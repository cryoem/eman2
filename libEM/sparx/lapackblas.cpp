#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "lapackblas.h"

integer ieeeck_(integer *ispec, real *zero, real *one)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1998   


    Purpose   
    =======   

    IEEECK is called from the ILAENV to verify that Infinity and   
    possibly NaN arithmetic is safe (i.e. will not trap).   

    Arguments   
    =========   

    ISPEC   (input) INTEGER   
            Specifies whether to test just for inifinity arithmetic   
            or whether to test for infinity and NaN arithmetic.   
            = 0: Verify infinity arithmetic only.   
            = 1: Verify infinity and NaN arithmetic.   

    ZERO    (input) REAL   
            Must contain the value 0.0   
            This is passed to prevent the compiler from optimizing   
            away this code.   

    ONE     (input) REAL   
            Must contain the value 1.0   
            This is passed to prevent the compiler from optimizing   
            away this code.   

    RETURN VALUE:  INTEGER   
            = 0:  Arithmetic failed to produce the correct answers   
            = 1:  Arithmetic produced the correct answers */
    /* System generated locals */
    integer ret_val;
    /* Local variables */
    static real neginf, posinf, negzro, newzro, nan1, nan2, nan3, nan4, nan5, 
	    nan6;


    ret_val = 1;

    posinf = *one / *zero;
    if (posinf <= *one) {
	ret_val = 0;
	return ret_val;
    }

    neginf = -(*one) / *zero;
    if (neginf >= *zero) {
	ret_val = 0;
	return ret_val;
    }

    negzro = *one / (neginf + *one);
    if (negzro != *zero) {
	ret_val = 0;
	return ret_val;
    }

    neginf = *one / negzro;
    if (neginf >= *zero) {
	ret_val = 0;
	return ret_val;
    }

    newzro = negzro + *zero;
    if (newzro != *zero) {
	ret_val = 0;
	return ret_val;
    }

    posinf = *one / newzro;
    if (posinf <= *one) {
	ret_val = 0;
	return ret_val;
    }

    neginf *= posinf;
    if (neginf >= *zero) {
	ret_val = 0;
	return ret_val;
    }

    posinf *= posinf;
    if (posinf <= *one) {
	ret_val = 0;
	return ret_val;
    }




/*     Return if we were only asked to check infinity arithmetic */

    if (*ispec == 0) {
	return ret_val;
    }

    nan1 = posinf + neginf;

    nan2 = posinf / neginf;

    nan3 = posinf / posinf;

    nan4 = posinf * *zero;

    nan5 = neginf * negzro;

    nan6 = nan5 * 0.f;

    if (nan1 == nan1) {
	ret_val = 0;
	return ret_val;
    }

    if (nan2 == nan2) {
	ret_val = 0;
	return ret_val;
    }

    if (nan3 == nan3) {
	ret_val = 0;
	return ret_val;
    }

    if (nan4 == nan4) {
	ret_val = 0;
	return ret_val;
    }

    if (nan5 == nan5) {
	ret_val = 0;
	return ret_val;
    }

    if (nan6 == nan6) {
	ret_val = 0;
	return ret_val;
    }

    return ret_val;
} /* ieeeck_ */




integer ilaenv_(integer *ispec, char *name__, char *opts, integer *n1, 
	integer *n2, integer *n3, integer *n4, ftnlen name_len, ftnlen 
	opts_len)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    ILAENV is called from the LAPACK routines to choose problem-dependent   
    parameters for the local environment.  See ISPEC for a description of   
    the parameters.   

    This version provides a set of parameters which should give good,   
    but not optimal, performance on many of the currently available   
    computers.  Users are encouraged to modify this subroutine to set   
    the tuning parameters for their particular machine using the option   
    and problem size information in the arguments.   

    This routine will not function correctly if it is converted to all   
    lower case.  Converting it to all upper case is allowed.   

    Arguments   
    =========   

    ISPEC   (input) INTEGER   
            Specifies the parameter to be returned as the value of   
            ILAENV.   
            = 1: the optimal blocksize; if this value is 1, an unblocked   
                 algorithm will give the best performance.   
            = 2: the minimum block size for which the block routine   
                 should be used; if the usable block size is less than   
                 this value, an unblocked routine should be used.   
            = 3: the crossover point (in a block routine, for N less   
                 than this value, an unblocked routine should be used)   
            = 4: the number of shifts, used in the nonsymmetric   
                 eigenvalue routines   
            = 5: the minimum column dimension for blocking to be used;   
                 rectangular blocks must have dimension at least k by m,   
                 where k is given by ILAENV(2,...) and m by ILAENV(5,...)   
            = 6: the crossover point for the SVD (when reducing an m by n   
                 matrix to bidiagonal form, if f2cmax(m,n)/min(m,n) exceeds   
                 this value, a QR factorization is used first to reduce   
                 the matrix to a triangular form.)   
            = 7: the number of processors   
            = 8: the crossover point for the multishift QR and QZ methods   
                 for nonsymmetric eigenvalue problems.   
            = 9: maximum size of the subproblems at the bottom of the   
                 computation tree in the divide-and-conquer algorithm   
                 (used by xGELSD and xGESDD)   
            =10: ieee NaN arithmetic can be trusted not to trap   
            =11: infinity arithmetic can be trusted not to trap   

    NAME    (input) CHARACTER*(*)   
            The name of the calling subroutine, in either upper case or   
            lower case.   

    OPTS    (input) CHARACTER*(*)   
            The character options to the subroutine NAME, concatenated   
            into a single character string.  For example, UPLO = 'U',   
            TRANS = 'T', and DIAG = 'N' for a triangular routine would   
            be specified as OPTS = 'UTN'.   

    N1      (input) INTEGER   
    N2      (input) INTEGER   
    N3      (input) INTEGER   
    N4      (input) INTEGER   
            Problem dimensions for the subroutine NAME; these may not all   
            be required.   

   (ILAENV) (output) INTEGER   
            >= 0: the value of the parameter specified by ISPEC   
            < 0:  if ILAENV = -k, the k-th argument had an illegal value.   

    Further Details   
    ===============   

    The following conventions have been used when calling ILAENV from the   
    LAPACK routines:   
    1)  OPTS is a concatenation of all of the character options to   
        subroutine NAME, in the same order that they appear in the   
        argument list for NAME, even if they are not used in determining   
        the value of the parameter specified by ISPEC.   
    2)  The problem dimensions N1, N2, N3, N4 are specified in the order   
        that they appear in the argument list for NAME.  N1 is used   
        first, N2 second, and so on, and unused problem dimensions are   
        passed a value of -1.   
    3)  The parameter value returned by ILAENV is checked for validity in   
        the calling subroutine.  For example, ILAENV is used to retrieve   
        the optimal blocksize for STRTRI as follows:   

        NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )   
        IF( NB.LE.1 ) NB = MAX( 1, N )   

    ===================================================================== */
    /* Table of constant values */
    static integer c__0 = 0;
    static real c_b162 = 0.f;
    static real c_b163 = 1.f;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer ret_val;
    /* Builtin functions   
       Subroutine */ void s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Local variables */
    static integer i__;
    static logical cname, sname;
    static integer nbmin;
    static char c1[1], c2[2], c3[3], c4[2];
    static integer ic, nb;
    extern integer ieeeck_(integer *, real *, real *);
    static integer iz, nx;
    static char subnam[6];




    switch (*ispec) {
	case 1:  goto L100;
	case 2:  goto L100;
	case 3:  goto L100;
	case 4:  goto L400;
	case 5:  goto L500;
	case 6:  goto L600;
	case 7:  goto L700;
	case 8:  goto L800;
	case 9:  goto L900;
	case 10:  goto L1000;
	case 11:  goto L1100;
    }

/*     Invalid value for ISPEC */

    ret_val = -1;
    return ret_val;

L100:

/*     Convert NAME to upper case if the first character is lower case. */

    ret_val = 1;
    s_copy(subnam, name__, (ftnlen)6, name_len);
    ic = *(unsigned char *)subnam;
    iz = 'Z';
    if (iz == 90 || iz == 122) {

/*        ASCII character set */

	if (ic >= 97 && ic <= 122) {
	    *(unsigned char *)subnam = (char) (ic - 32);
	    for (i__ = 2; i__ <= 6; ++i__) {
		ic = *(unsigned char *)&subnam[i__ - 1];
		if (ic >= 97 && ic <= 122) {
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
		}
/* L10: */
	    }
	}

    } else if (iz == 233 || iz == 169) {

/*        EBCDIC character set */

	if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 && 
		ic <= 169) {
	    *(unsigned char *)subnam = (char) (ic + 64);
	    for (i__ = 2; i__ <= 6; ++i__) {
		ic = *(unsigned char *)&subnam[i__ - 1];
		if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 
			162 && ic <= 169) {
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
		}
/* L20: */
	    }
	}

    } else if (iz == 218 || iz == 250) {

/*        Prime machines:  ASCII+128 */

	if (ic >= 225 && ic <= 250) {
	    *(unsigned char *)subnam = (char) (ic - 32);
	    for (i__ = 2; i__ <= 6; ++i__) {
		ic = *(unsigned char *)&subnam[i__ - 1];
		if (ic >= 225 && ic <= 250) {
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
		}
/* L30: */
	    }
	}
    }

    *(unsigned char *)c1 = *(unsigned char *)subnam;
    sname = *(unsigned char *)c1 == 'S' || *(unsigned char *)c1 == 'D';
    cname = *(unsigned char *)c1 == 'C' || *(unsigned char *)c1 == 'Z';
    if (! (cname || sname)) {
	return ret_val;
    }
    s_copy(c2, subnam + 1, (ftnlen)2, (ftnlen)2);
    s_copy(c3, subnam + 3, (ftnlen)3, (ftnlen)3);
    s_copy(c4, c3 + 1, (ftnlen)2, (ftnlen)2);

    switch (*ispec) {
	case 1:  goto L110;
	case 2:  goto L200;
	case 3:  goto L300;
    }

L110:

/*     ISPEC = 1:  block size   

       In these examples, separate code is provided for setting NB for   
       real and complex.  We assume that NB will take the same value in   
       single or double precision. */

    nb = 1;

    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	} else if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, 
		"RQF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)
		3, (ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) 
		== 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "PO", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	} else if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
	    nb = 32;
	} else if (sname && s_cmp(c3, "GST", (ftnlen)3, (ftnlen)3) == 0) {
	    nb = 64;
	}
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    nb = 64;
	} else if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
	    nb = 32;
	} else if (s_cmp(c3, "GST", (ftnlen)3, (ftnlen)3) == 0) {
	    nb = 64;
	}
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nb = 32;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nb = 32;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nb = 32;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nb = 32;
	    }
	}
    } else if (s_cmp(c2, "GB", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		if (*n4 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    } else {
		if (*n4 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    }
	}
    } else if (s_cmp(c2, "PB", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		if (*n2 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    } else {
		if (*n2 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    }
	}
    } else if (s_cmp(c2, "TR", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "LA", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "UUM", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (sname && s_cmp(c2, "ST", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "EBZ", (ftnlen)3, (ftnlen)3) == 0) {
	    nb = 1;
	}
    }
    ret_val = nb;
    return ret_val;

L200:

/*     ISPEC = 2:  minimum block size */

    nbmin = 2;
    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "RQF", (
		ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)3, (
		ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
		 {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	}
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nbmin = 8;
	    } else {
		nbmin = 8;
	    }
	} else if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
	    nbmin = 2;
	}
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
	    nbmin = 2;
	}
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nbmin = 2;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nbmin = 2;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nbmin = 2;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nbmin = 2;
	    }
	}
    }
    ret_val = nbmin;
    return ret_val;

L300:

/*     ISPEC = 3:  crossover point */

    nx = 0;
    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "RQF", (
		ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)3, (
		ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
		 {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	}
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
	if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
	    nx = 32;
	}
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
	    nx = 32;
	}
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nx = 128;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nx = 128;
	    }
	}
    }
    ret_val = nx;
    return ret_val;

L400:

/*     ISPEC = 4:  number of shifts (used by xHSEQR) */

    ret_val = 6;
    return ret_val;

L500:

/*     ISPEC = 5:  minimum column dimension (not used) */

    ret_val = 2;
    return ret_val;

L600:

/*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD) */

    ret_val = (integer) ((real) f2cmin(*n1,*n2) * 1.6f);
    return ret_val;

L700:

/*     ISPEC = 7:  number of processors (not used) */

    ret_val = 1;
    return ret_val;

L800:

/*     ISPEC = 8:  crossover point for multishift (used by xHSEQR) */

    ret_val = 50;
    return ret_val;

L900:

/*     ISPEC = 9:  maximum size of the subproblems at the bottom of the   
                   computation tree in the divide-and-conquer algorithm   
                   (used by xGELSD and xGESDD) */

    ret_val = 25;
    return ret_val;

L1000:

/*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap   

       ILAENV = 0 */
    ret_val = 1;
    if (ret_val == 1) {
	ret_val = ieeeck_(&c__0, &c_b162, &c_b163);
    }
    return ret_val;

L1100:

/*     ISPEC = 11: infinity arithmetic can be trusted not to trap   

       ILAENV = 0 */
    ret_val = 1;
    if (ret_val == 1) {
	ret_val = ieeeck_(&c__1, &c_b162, &c_b163);
    }
    return ret_val;

/*     End of ILAENV */

} /* ilaenv_ */



logical lsame_(char *ca, char *cb)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    LSAME returns .TRUE. if CA is the same letter as CB regardless of   
    case.   

    Arguments   
    =========   

    CA      (input) CHARACTER*1   
    CB      (input) CHARACTER*1   
            CA and CB specify the single characters to be compared.   

   ===================================================================== 
  


       Test if the characters are equal */
    /* System generated locals */
    logical ret_val;
    /* Local variables */
    static integer inta, intb, zcode;


    ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
    if (ret_val) {
	return ret_val;
    }

/*     Now test for equivalence if both characters are alphabetic. */

    zcode = 'Z';

/*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime   
       machines, on which ICHAR returns a value with bit 8 set.   
       ICHAR('A') on Prime machines returns 193 which is the same as   
       ICHAR('A') on an EBCDIC machine. */

    inta = *(unsigned char *)ca;
    intb = *(unsigned char *)cb;

    if (zcode == 90 || zcode == 122) {

/*        ASCII is assumed - ZCODE is the ASCII code of either lower o
r   
          upper case 'Z'. */

	if (inta >= 97 && inta <= 122) {
	    inta += -32;
	}
	if (intb >= 97 && intb <= 122) {
	    intb += -32;
	}

    } else if (zcode == 233 || zcode == 169) {

/*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower
 or   
          upper case 'Z'. */

	if (inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta 
		>= 162 && inta <= 169) {
	    inta += 64;
	}
	if (intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb 
		>= 162 && intb <= 169) {
	    intb += 64;
	}

    } else if (zcode == 218 || zcode == 250) {

/*        ASCII is assumed, on Prime machines - ZCODE is the ASCII cod
e   
          plus 128 of either lower or upper case 'Z'. */

	if (inta >= 225 && inta <= 250) {
	    inta += -32;
	}
	if (intb >= 225 && intb <= 250) {
	    intb += -32;
	}
    }
    ret_val = inta == intb;

/*     RETURN   

       End of LSAME */

    return ret_val;
} /* lsame_ */



#ifdef KR_headers
double pow_ri(ap, bp) real *ap; integer *bp;
#else
double pow_ri(real *ap, integer *bp)
#endif
{
double pow, x;
integer n;
unsigned long u;

pow = 1;
x = *ap;
n = *bp;

if(n != 0)
	{
	if(n < 0)
		{
		n = -n;
		x = 1/x;
		}
	for(u = n; ; )
		{
		if(u & 01)
			pow *= x;
		if(u >>= 1)
			x *= x;
		else
			break;
		}
	}
return(pow);
}


#ifdef KR_headers
double r_sign(a,b) real *a, *b;
#else
double r_sign(real *a, real *b)
#endif
{
double x;
x = (*a >= 0 ? *a : - *a);
return( *b >= 0 ? x : -x);
}



/* Subroutine */ int saxpy_(integer *n, real *sa, real *sx, integer *incx, 
	real *sy, integer *incy)
{
    /* System generated locals */
    integer i__1;
    /* Local variables */
    static integer i__, m, ix, iy, mp1;
/*     constant times a vector plus a vector.   
       uses unrolled loop for increments equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   
       Parameter adjustments */
    --sy;
    --sx;
    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*sa == 0.f) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }
/*        code for unequal increments or equal increments   
            not equal to 1 */
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sy[iy] += *sa * sx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;
/*        code for both increments equal to 1   
          clean-up loop */
L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sy[i__] += *sa * sx[i__];
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 4) {
	sy[i__] += *sa * sx[i__];
	sy[i__ + 1] += *sa * sx[i__ + 1];
	sy[i__ + 2] += *sa * sx[i__ + 2];
	sy[i__ + 3] += *sa * sx[i__ + 3];
/* L50: */
    }
    return 0;
} /* saxpy_ */



/* compare two strings */

#ifdef KR_headers
integer s_cmp(a0, b0, la, lb) char *a0, *b0; ftnlen la, lb;
#else
integer s_cmp(char *a0, char *b0, ftnlen la, ftnlen lb)
#endif
{
register unsigned char *a, *aend, *b, *bend;
a = (unsigned char *)a0;
b = (unsigned char *)b0;
aend = a + la;
bend = b + lb;

if(la <= lb)
	{
	while(a < aend)
		if(*a != *b)
			return( *a - *b );
		else
			{ ++a; ++b; }

	while(b < bend)
		if(*b != ' ')
			return( ' ' - *b );
		else	++b;
	}

else
	{
	while(b < bend)
		if(*a == *b)
			{ ++a; ++b; }
		else
			return( *a - *b );
	while(a < aend)
		if(*a != ' ')
			return(*a - ' ');
		else	++a;
	}
return(0);
}
/* Unless compiled with -DNO_OVERWRITE, this variant of s_copy allows the
 * target of an assignment to appear on its right-hand side (contrary
 * to the Fortran 77 Standard, but in accordance with Fortran 90),
 * as in  a(2:5) = a(4:7) .
 */



/* assign strings:  a = b */

#ifdef KR_headers
VOID s_copy(a, b, la, lb) register char *a, *b; ftnlen la, lb;
#else
void s_copy(char *a, char *b, ftnlen la, ftnlen lb)
#endif
{
	register char *aend, *bend;

	aend = a + la;

	if(la <= lb)
#ifndef NO_OVERWRITE
		if (a <= b || a >= b + la)
#endif
			while(a < aend)
				*a++ = *b++;
#ifndef NO_OVERWRITE
		else
			for(b += la; a < aend; )
				*--aend = *--b;
#endif

	else {
		bend = b + lb;
#ifndef NO_OVERWRITE
		if (a <= b || a >= bend)
#endif
			while(b < bend)
				*a++ = *b++;
#ifndef NO_OVERWRITE
		else {
			a += lb;
			while(b < bend)
				*--a = *--bend;
			a += lb;
			}
#endif
		while(a < aend)
			*a++ = ' ';
		}
	}



/* Subroutine */ int scopy_(integer *n, real *sx, integer *incx, real *sy, 
	integer *incy)
{
    /* System generated locals */
    integer i__1;
    /* Local variables */
    static integer i__, m, ix, iy, mp1;
/*     copies a vector, x, to a vector, y.   
       uses unrolled loops for increments equal to 1.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   
       Parameter adjustments */
    --sy;
    --sx;
    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }
/*        code for unequal increments or equal increments   
            not equal to 1 */
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sy[iy] = sx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;
/*        code for both increments equal to 1   
          clean-up loop */
L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sy[i__] = sx[i__];
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 7) {
	sy[i__] = sx[i__];
	sy[i__ + 1] = sx[i__ + 1];
	sy[i__ + 2] = sx[i__ + 2];
	sy[i__ + 3] = sx[i__ + 3];
	sy[i__ + 4] = sx[i__ + 4];
	sy[i__ + 5] = sx[i__ + 5];
	sy[i__ + 6] = sx[i__ + 6];
/* L50: */
    }
    return 0;
} /* scopy_ */




doublereal sdot_(integer *n, real *sx, integer *incx, real *sy, integer *incy)
{
    /* System generated locals */
    integer i__1;
    real ret_val;
    /* Local variables */
    static integer i__, m;
    static real stemp;
    static integer ix, iy, mp1;
/*     forms the dot product of two vectors.   
       uses unrolled loops for increments equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   
       Parameter adjustments */
    --sy;
    --sx;
    /* Function Body */
    stemp = 0.f;
    ret_val = 0.f;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }
/*        code for unequal increments or equal increments   
            not equal to 1 */
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stemp += sx[ix] * sy[iy];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    ret_val = stemp;
    return ret_val;
/*        code for both increments equal to 1   
          clean-up loop */
L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stemp += sx[i__] * sy[i__];
/* L30: */
    }
    if (*n < 5) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
	stemp = stemp + sx[i__] * sy[i__] + sx[i__ + 1] * sy[i__ + 1] + sx[
		i__ + 2] * sy[i__ + 2] + sx[i__ + 3] * sy[i__ + 3] + sx[i__ + 
		4] * sy[i__ + 4];
/* L50: */
    }
L60:
    ret_val = stemp;
    return ret_val;
} /* sdot_ */




/* Subroutine */ int sgemm_(char *transa, char *transb, integer *m, integer *
	n, integer *k, real *alpha, real *a, integer *lda, real *b, integer *
	ldb, real *beta, real *c__, integer *ldc)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;
    /* Local variables */
    static integer info;
    static logical nota, notb;
    static real temp;
    static integer i__, j, l, ncola;
    extern logical lsame_(char *, char *);
    static integer nrowa, nrowb;
    extern /* Subroutine */ int xerbla_(char *, integer *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
#define c___ref(a_1,a_2) c__[(a_2)*c_dim1 + a_1]
/*  Purpose   
    =======   
    SGEMM  performs one of the matrix-matrix operations   
       C := alpha*op( A )*op( B ) + beta*C,   
    where  op( X ) is one of   
       op( X ) = X   or   op( X ) = X',   
    alpha and beta are scalars, and A, B and C are matrices, with op( A )   
    an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.   
    Parameters   
    ==========   
    TRANSA - CHARACTER*1.   
             On entry, TRANSA specifies the form of op( A ) to be used in   
             the matrix multiplication as follows:   
                TRANSA = 'N' or 'n',  op( A ) = A.   
                TRANSA = 'T' or 't',  op( A ) = A'.   
                TRANSA = 'C' or 'c',  op( A ) = A'.   
             Unchanged on exit.   
    TRANSB - CHARACTER*1.   
             On entry, TRANSB specifies the form of op( B ) to be used in   
             the matrix multiplication as follows:   
                TRANSB = 'N' or 'n',  op( B ) = B.   
                TRANSB = 'T' or 't',  op( B ) = B'.   
                TRANSB = 'C' or 'c',  op( B ) = B'.   
             Unchanged on exit.   
    M      - INTEGER.   
             On entry,  M  specifies  the number  of rows  of the  matrix   
             op( A )  and of the  matrix  C.  M  must  be at least  zero.   
             Unchanged on exit.   
    N      - INTEGER.   
             On entry,  N  specifies the number  of columns of the matrix   
             op( B ) and the number of columns of the matrix C. N must be   
             at least zero.   
             Unchanged on exit.   
    K      - INTEGER.   
             On entry,  K  specifies  the number of columns of the matrix   
             op( A ) and the number of rows of the matrix op( B ). K must   
             be at least  zero.   
             Unchanged on exit.   
    ALPHA  - REAL            .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   
    A      - REAL             array of DIMENSION ( LDA, ka ), where ka is   
             k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.   
             Before entry with  TRANSA = 'N' or 'n',  the leading  m by k   
             part of the array  A  must contain the matrix  A,  otherwise   
             the leading  k by m  part of the array  A  must contain  the   
             matrix A.   
             Unchanged on exit.   
    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared   
             in the calling (sub) program. When  TRANSA = 'N' or 'n' then   
             LDA must be at least  f2cmax( 1, m ), otherwise  LDA must be at   
             least  f2cmax( 1, k ).   
             Unchanged on exit.   
    B      - REAL             array of DIMENSION ( LDB, kb ), where kb is   
             n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.   
             Before entry with  TRANSB = 'N' or 'n',  the leading  k by n   
             part of the array  B  must contain the matrix  B,  otherwise   
             the leading  n by k  part of the array  B  must contain  the   
             matrix B.   
             Unchanged on exit.   
    LDB    - INTEGER.   
             On entry, LDB specifies the first dimension of B as declared   
             in the calling (sub) program. When  TRANSB = 'N' or 'n' then   
             LDB must be at least  f2cmax( 1, k ), otherwise  LDB must be at   
             least  f2cmax( 1, n ).   
             Unchanged on exit.   
    BETA   - REAL            .   
             On entry,  BETA  specifies the scalar  beta.  When  BETA  is   
             supplied as zero then C need not be set on input.   
             Unchanged on exit.   
    C      - REAL             array of DIMENSION ( LDC, n ).   
             Before entry, the leading  m by n  part of the array  C must   
             contain the matrix  C,  except when  beta  is zero, in which   
             case C need not be set on entry.   
             On exit, the array  C  is overwritten by the  m by n  matrix   
             ( alpha*op( A )*op( B ) + beta*C ).   
    LDC    - INTEGER.   
             On entry, LDC specifies the first dimension of C as declared   
             in  the  calling  (sub)  program.   LDC  must  be  at  least   
             f2cmax( 1, m ).   
             Unchanged on exit.   
    Level 3 Blas routine.   
    -- Written on 8-February-1989.   
       Jack Dongarra, Argonne National Laboratory.   
       Iain Duff, AERE Harwell.   
       Jeremy Du Croz, Numerical Algorithms Group Ltd.   
       Sven Hammarling, Numerical Algorithms Group Ltd.   
       Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not   
       transposed and set  NROWA, NCOLA and  NROWB  as the number of rows   
       and  columns of  A  and the  number of  rows  of  B  respectively.   
       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;
    /* Function Body */
    nota = lsame_(transa, "N");
    notb = lsame_(transb, "N");
    if (nota) {
	nrowa = *m;
	ncola = *k;
    } else {
	nrowa = *k;
	ncola = *m;
    }
    if (notb) {
	nrowb = *k;
    } else {
	nrowb = *n;
    }
/*     Test the input parameters. */
    info = 0;
    if (! nota && ! lsame_(transa, "C") && ! lsame_(
	    transa, "T")) {
	info = 1;
    } else if (! notb && ! lsame_(transb, "C") && ! 
	    lsame_(transb, "T")) {
	info = 2;
    } else if (*m < 0) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*k < 0) {
	info = 5;
    } else if (*lda < f2cmax(1,nrowa)) {
	info = 8;
    } else if (*ldb < f2cmax(1,nrowb)) {
	info = 10;
    } else if (*ldc < f2cmax(1,*m)) {
	info = 13;
    }
    if (info != 0) {
	xerbla_("SGEMM ", &info);
	return 0;
    }
/*     Quick return if possible. */
    if (*m == 0 || *n == 0 || (*alpha == 0.f || *k == 0) && *beta == 1.f) {
	return 0;
    }
/*     And if  alpha.eq.zero. */
    if (*alpha == 0.f) {
	if (*beta == 0.f) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c___ref(i__, j) = 0.f;
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c___ref(i__, j) = *beta * c___ref(i__, j);
/* L30: */
		}
/* L40: */
	    }
	}
	return 0;
    }
/*     Start the operations. */
    if (notb) {
	if (nota) {
/*           Form  C := alpha*A*B + beta*C. */
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.f) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = 0.f;
/* L50: */
		    }
		} else if (*beta != 1.f) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = *beta * c___ref(i__, j);
/* L60: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (b_ref(l, j) != 0.f) {
			temp = *alpha * b_ref(l, j);
			i__3 = *m;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    c___ref(i__, j) = c___ref(i__, j) + temp * a_ref(
				    i__, l);
/* L70: */
			}
		    }
/* L80: */
		}
/* L90: */
	    }
	} else {
/*           Form  C := alpha*A'*B + beta*C */
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp = 0.f;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp += a_ref(l, i__) * b_ref(l, j);
/* L100: */
		    }
		    if (*beta == 0.f) {
			c___ref(i__, j) = *alpha * temp;
		    } else {
			c___ref(i__, j) = *alpha * temp + *beta * c___ref(i__,
				 j);
		    }
/* L110: */
		}
/* L120: */
	    }
	}
    } else {
	if (nota) {
/*           Form  C := alpha*A*B' + beta*C */
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.f) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = 0.f;
/* L130: */
		    }
		} else if (*beta != 1.f) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = *beta * c___ref(i__, j);
/* L140: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (b_ref(j, l) != 0.f) {
			temp = *alpha * b_ref(j, l);
			i__3 = *m;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    c___ref(i__, j) = c___ref(i__, j) + temp * a_ref(
				    i__, l);
/* L150: */
			}
		    }
/* L160: */
		}
/* L170: */
	    }
	} else {
/*           Form  C := alpha*A'*B' + beta*C */
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp = 0.f;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp += a_ref(l, i__) * b_ref(j, l);
/* L180: */
		    }
		    if (*beta == 0.f) {
			c___ref(i__, j) = *alpha * temp;
		    } else {
			c___ref(i__, j) = *alpha * temp + *beta * c___ref(i__,
				 j);
		    }
/* L190: */
		}
/* L200: */
	    }
	}
    }
    return 0;
/*     End of SGEMM . */
} /* sgemm_ */
#undef c___ref
#undef b_ref
#undef a_ref




/* Subroutine */ int sgemv_(char *trans, integer *m, integer *n, real *alpha, 
	real *a, integer *lda, real *x, integer *incx, real *beta, real *y, 
	integer *incy)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    static integer info;
    static real temp;
    static integer lenx, leny, i__, j;
    extern logical lsame_(char *, char *);
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_(char *, integer *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
/*  Purpose   
    =======   
    SGEMV  performs one of the matrix-vector operations   
       y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   
    where alpha and beta are scalars, x and y are vectors and A is an   
    m by n matrix.   
    Parameters   
    ==========   
    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the operation to be performed as   
             follows:   
                TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.   
                TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.   
                TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.   
             Unchanged on exit.   
    M      - INTEGER.   
             On entry, M specifies the number of rows of the matrix A.   
             M must be at least zero.   
             Unchanged on exit.   
    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   
    ALPHA  - REAL            .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   
    A      - REAL             array of DIMENSION ( LDA, n ).   
             Before entry, the leading m by n part of the array A must   
             contain the matrix of coefficients.   
             Unchanged on exit.   
    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared   
             in the calling (sub) program. LDA must be at least   
             f2cmax( 1, m ).   
             Unchanged on exit.   
    X      - REAL             array of DIMENSION at least   
             ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.   
             Before entry, the incremented array X must contain the   
             vector x.   
             Unchanged on exit.   
    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   
    BETA   - REAL            .   
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   
             Unchanged on exit.   
    Y      - REAL             array of DIMENSION at least   
             ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.   
             Before entry with BETA non-zero, the incremented array Y   
             must contain the vector y. On exit, Y is overwritten by the   
             updated vector y.   
    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   
    Level 2 Blas routine.   
    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   
       Test the input parameters.   
       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --x;
    --y;
    /* Function Body */
    info = 0;
    if (! lsame_(trans, "N") && ! lsame_(trans, "T") && ! lsame_(trans, "C")
	    ) {
	info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*lda < f2cmax(1,*m)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("SGEMV ", &info);
	return 0;
    }
/*     Quick return if possible. */
    if (*m == 0 || *n == 0 || *alpha == 0.f && *beta == 1.f) {
	return 0;
    }
/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set   
       up the start points in  X  and  Y. */
    if (lsame_(trans, "N")) {
	lenx = *n;
	leny = *m;
    } else {
	lenx = *m;
	leny = *n;
    }
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * *incy;
    }
/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A.   
       First form  y := beta*y. */
    if (*beta != 1.f) {
	if (*incy == 1) {
	    if (*beta == 0.f) {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = 0.f;
/* L10: */
		}
	    } else {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = *beta * y[i__];
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.f) {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = 0.f;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = *beta * y[iy];
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == 0.f) {
	return 0;
    }
    if (lsame_(trans, "N")) {
/*        Form  y := alpha*A*x + y. */
	jx = kx;
	if (*incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0.f) {
		    temp = *alpha * x[jx];
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			y[i__] += temp * a_ref(i__, j);
/* L50: */
		    }
		}
		jx += *incx;
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0.f) {
		    temp = *alpha * x[jx];
		    iy = ky;
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			y[iy] += temp * a_ref(i__, j);
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
/* L80: */
	    }
	}
    } else {
/*        Form  y := alpha*A'*x + y. */
	jy = ky;
	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp = 0.f;
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp += a_ref(i__, j) * x[i__];
/* L90: */
		}
		y[jy] += *alpha * temp;
		jy += *incy;
/* L100: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp = 0.f;
		ix = kx;
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp += a_ref(i__, j) * x[ix];
		    ix += *incx;
/* L110: */
		}
		y[jy] += *alpha * temp;
		jy += *incy;
/* L120: */
	    }
	}
    }
    return 0;
/*     End of SGEMV . */
} /* sgemv_ */
#undef a_ref




/* Subroutine */ int sger_(integer *m, integer *n, real *alpha, real *x, 
	integer *incx, real *y, integer *incy, real *a, integer *lda)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    static integer info;
    static real temp;
    static integer i__, j, ix, jy, kx;
    extern /* Subroutine */ int xerbla_(char *, integer *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
/*  Purpose   
    =======   
    SGER   performs the rank 1 operation   
       A := alpha*x*y' + A,   
    where alpha is a scalar, x is an m element vector, y is an n element   
    vector and A is an m by n matrix.   
    Parameters   
    ==========   
    M      - INTEGER.   
             On entry, M specifies the number of rows of the matrix A.   
             M must be at least zero.   
             Unchanged on exit.   
    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   
    ALPHA  - REAL            .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   
    X      - REAL             array of dimension at least   
             ( 1 + ( m - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the m   
             element vector x.   
             Unchanged on exit.   
    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   
    Y      - REAL             array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y.   
             Unchanged on exit.   
    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   
    A      - REAL             array of DIMENSION ( LDA, n ).   
             Before entry, the leading m by n part of the array A must   
             contain the matrix of coefficients. On exit, A is   
             overwritten by the updated matrix.   
    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared   
             in the calling (sub) program. LDA must be at least   
             f2cmax( 1, m ).   
             Unchanged on exit.   
    Level 2 Blas routine.   
    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   
       Test the input parameters.   
       Parameter adjustments */
    --x;
    --y;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    /* Function Body */
    info = 0;
    if (*m < 0) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < f2cmax(1,*m)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("SGER  ", &info);
	return 0;
    }
/*     Quick return if possible. */
    if (*m == 0 || *n == 0 || *alpha == 0.f) {
	return 0;
    }
/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A. */
    if (*incy > 0) {
	jy = 1;
    } else {
	jy = 1 - (*n - 1) * *incy;
    }
    if (*incx == 1) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (y[jy] != 0.f) {
		temp = *alpha * y[jy];
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    a_ref(i__, j) = a_ref(i__, j) + x[i__] * temp;
/* L10: */
		}
	    }
	    jy += *incy;
/* L20: */
	}
    } else {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*m - 1) * *incx;
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (y[jy] != 0.f) {
		temp = *alpha * y[jy];
		ix = kx;
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    a_ref(i__, j) = a_ref(i__, j) + x[ix] * temp;
		    ix += *incx;
/* L30: */
		}
	    }
	    jy += *incy;
/* L40: */
	}
    }
    return 0;
/*     End of SGER  . */
} /* sger_ */
#undef a_ref




/* Subroutine */ int slae2_(real *a, real *b, real *c__, real *rt1, real *rt2)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix   
       [  A   B  ]   
       [  B   C  ].   
    On return, RT1 is the eigenvalue of larger absolute value, and RT2   
    is the eigenvalue of smaller absolute value.   

    Arguments   
    =========   

    A       (input) REAL   
            The (1,1) element of the 2-by-2 matrix.   

    B       (input) REAL   
            The (1,2) and (2,1) elements of the 2-by-2 matrix.   

    C       (input) REAL   
            The (2,2) element of the 2-by-2 matrix.   

    RT1     (output) REAL   
            The eigenvalue of larger absolute value.   

    RT2     (output) REAL   
            The eigenvalue of smaller absolute value.   

    Further Details   
    ===============   

    RT1 is accurate to a few ulps barring over/underflow.   

    RT2 may be inaccurate if there is massive cancellation in the   
    determinant A*C-B*B; higher precision or correctly rounded or   
    correctly truncated arithmetic would be needed to compute RT2   
    accurately in all cases.   

    Overflow is possible only if RT1 is within a factor of 5 of overflow.   
    Underflow is harmless if the input data is 0 or exceeds   
       underflow_threshold / macheps.   

   =====================================================================   


       Compute the eigenvalues */
    /* System generated locals */
    real r__1;
    /* Builtin functions */
//    double sqrt(doublereal);
    /* Local variables */
    static real acmn, acmx, ab, df, tb, sm, rt, adf;


    sm = *a + *c__;
    df = *a - *c__;
    adf = dabs(df);
    tb = *b + *b;
    ab = dabs(tb);
    if (dabs(*a) > dabs(*c__)) {
	acmx = *a;
	acmn = *c__;
    } else {
	acmx = *c__;
	acmn = *a;
    }
    if (adf > ab) {
/* Computing 2nd power */
	r__1 = ab / adf;
	rt = adf * sqrt(r__1 * r__1 + 1.f);
    } else if (adf < ab) {
/* Computing 2nd power */
	r__1 = adf / ab;
	rt = ab * sqrt(r__1 * r__1 + 1.f);
    } else {

/*        Includes case AB=ADF=0 */

	rt = ab * sqrt(2.f);
    }
    if (sm < 0.f) {
	*rt1 = (sm - rt) * .5f;

/*        Order of execution important.   
          To get fully accurate smaller eigenvalue,   
          next line needs to be executed in higher precision. */

	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else if (sm > 0.f) {
	*rt1 = (sm + rt) * .5f;

/*        Order of execution important.   
          To get fully accurate smaller eigenvalue,   
          next line needs to be executed in higher precision. */

	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else {

/*        Includes case RT1 = RT2 = 0 */

	*rt1 = rt * .5f;
	*rt2 = rt * -.5f;
    }
    return 0;

/*     End of SLAE2 */

} /* slae2_ */




/* Subroutine */ int slaev2_(real *a, real *b, real *c__, real *rt1, real *
	rt2, real *cs1, real *sn1)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix   
       [  A   B  ]   
       [  B   C  ].   
    On return, RT1 is the eigenvalue of larger absolute value, RT2 is the   
    eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right   
    eigenvector for RT1, giving the decomposition   

       [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]   
       [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].   

    Arguments   
    =========   

    A       (input) REAL   
            The (1,1) element of the 2-by-2 matrix.   

    B       (input) REAL   
            The (1,2) element and the conjugate of the (2,1) element of   
            the 2-by-2 matrix.   

    C       (input) REAL   
            The (2,2) element of the 2-by-2 matrix.   

    RT1     (output) REAL   
            The eigenvalue of larger absolute value.   

    RT2     (output) REAL   
            The eigenvalue of smaller absolute value.   

    CS1     (output) REAL   
    SN1     (output) REAL   
            The vector (CS1, SN1) is a unit right eigenvector for RT1.   

    Further Details   
    ===============   

    RT1 is accurate to a few ulps barring over/underflow.   

    RT2 may be inaccurate if there is massive cancellation in the   
    determinant A*C-B*B; higher precision or correctly rounded or   
    correctly truncated arithmetic would be needed to compute RT2   
    accurately in all cases.   

    CS1 and SN1 are accurate to a few ulps barring over/underflow.   

    Overflow is possible only if RT1 is within a factor of 5 of overflow.   
    Underflow is harmless if the input data is 0 or exceeds   
       underflow_threshold / macheps.   

   =====================================================================   


       Compute the eigenvalues */
    /* System generated locals */
    real r__1;
    /* Builtin functions */
//    double sqrt(doublereal);
    /* Local variables */
    static real acmn, acmx, ab, df, cs, ct, tb, sm, tn, rt, adf, acs;
    static integer sgn1, sgn2;


    sm = *a + *c__;
    df = *a - *c__;
    adf = dabs(df);
    tb = *b + *b;
    ab = dabs(tb);
    if (dabs(*a) > dabs(*c__)) {
	acmx = *a;
	acmn = *c__;
    } else {
	acmx = *c__;
	acmn = *a;
    }
    if (adf > ab) {
/* Computing 2nd power */
	r__1 = ab / adf;
	rt = adf * sqrt(r__1 * r__1 + 1.f);
    } else if (adf < ab) {
/* Computing 2nd power */
	r__1 = adf / ab;
	rt = ab * sqrt(r__1 * r__1 + 1.f);
    } else {

/*        Includes case AB=ADF=0 */

	rt = ab * sqrt(2.f);
    }
    if (sm < 0.f) {
	*rt1 = (sm - rt) * .5f;
	sgn1 = -1;

/*        Order of execution important.   
          To get fully accurate smaller eigenvalue,   
          next line needs to be executed in higher precision. */

	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else if (sm > 0.f) {
	*rt1 = (sm + rt) * .5f;
	sgn1 = 1;

/*        Order of execution important.   
          To get fully accurate smaller eigenvalue,   
          next line needs to be executed in higher precision. */

	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else {

/*        Includes case RT1 = RT2 = 0 */

	*rt1 = rt * .5f;
	*rt2 = rt * -.5f;
	sgn1 = 1;
    }

/*     Compute the eigenvector */

    if (df >= 0.f) {
	cs = df + rt;
	sgn2 = 1;
    } else {
	cs = df - rt;
	sgn2 = -1;
    }
    acs = dabs(cs);
    if (acs > ab) {
	ct = -tb / cs;
	*sn1 = 1.f / sqrt(ct * ct + 1.f);
	*cs1 = ct * *sn1;
    } else {
	if (ab == 0.f) {
	    *cs1 = 1.f;
	    *sn1 = 0.f;
	} else {
	    tn = -cs / tb;
	    *cs1 = 1.f / sqrt(tn * tn + 1.f);
	    *sn1 = tn * *cs1;
	}
    }
    if (sgn1 == sgn2) {
	tn = *cs1;
	*cs1 = -(*sn1);
	*sn1 = tn;
    }
    return 0;

/*     End of SLAEV2 */

} /* slaev2_ */



doublereal slamch_(char *cmach)
{
/*  -- LAPACK auxiliary routine (version 3.0) --
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLAMCH determines single precision machine parameters.   

    Arguments   
    =========   

    CMACH   (input) CHARACTER*1   
            Specifies the value to be returned by SLAMCH:   
            = 'E' or 'e',   SLAMCH := eps   
            = 'S' or 's ,   SLAMCH := sfmin   
            = 'B' or 'b',   SLAMCH := base   
            = 'P' or 'p',   SLAMCH := eps*base   
            = 'N' or 'n',   SLAMCH := t   
            = 'R' or 'r',   SLAMCH := rnd   
            = 'M' or 'm',   SLAMCH := emin   
            = 'U' or 'u',   SLAMCH := rmin   
            = 'L' or 'l',   SLAMCH := emax   
            = 'O' or 'o',   SLAMCH := rmax   

            where   

            eps   = relative machine precision   
            sfmin = safe minimum, such that 1/sfmin does not overflow   
            base  = base of the machine   
            prec  = eps*base   
            t     = number of (base) digits in the mantissa   
            rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise   
            emin  = minimum exponent before (gradual) underflow   
            rmin  = underflow threshold - base**(emin-1)   
            emax  = largest exponent before overflow   
            rmax  = overflow threshold  - (base**emax)*(1-eps)   

   ===================================================================== 
*/
/* >>Start of File<<   
       Initialized data */
    static logical first = TRUE_;
    /* System generated locals */
    integer i__1;
    real ret_val;
    /* Builtin functions */
    double pow_ri(real *, integer *);
    /* Local variables */
    static real base;
    static integer beta;
    static real emin, prec, emax;
    static integer imin, imax;
    static logical lrnd;
    static real rmin, rmax, t, rmach;
    extern logical lsame_(char *, char *);
    static real small, sfmin;
    extern /* Subroutine */ int slamc2_(integer *, integer *, logical *, real 
	    *, integer *, real *, integer *, real *);
    static integer it;
    static real rnd, eps;



    if (first) {
	first = FALSE_;
	slamc2_(&beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax);
	base = (real) beta;
	t = (real) it;
	if (lrnd) {
	    rnd = 1.f;
	    i__1 = 1 - it;
	    eps = pow_ri(&base, &i__1) / 2;
	} else {
	    rnd = 0.f;
	    i__1 = 1 - it;
	    eps = pow_ri(&base, &i__1);
	}
	prec = eps * base;
	emin = (real) imin;
	emax = (real) imax;
	sfmin = rmin;
	small = 1.f / rmax;
	if (small >= sfmin) {

/*           Use SMALL plus a bit, to avoid the possibility of rou
nding   
             causing overflow when computing  1/sfmin. */

	    sfmin = small * (eps + 1.f);
	}
    }

    if (lsame_(cmach, "E")) {
	rmach = eps;
    } else if (lsame_(cmach, "S")) {
	rmach = sfmin;
    } else if (lsame_(cmach, "B")) {
	rmach = base;
    } else if (lsame_(cmach, "P")) {
	rmach = prec;
    } else if (lsame_(cmach, "N")) {
	rmach = t;
    } else if (lsame_(cmach, "R")) {
	rmach = rnd;
    } else if (lsame_(cmach, "M")) {
	rmach = emin;
    } else if (lsame_(cmach, "U")) {
	rmach = rmin;
    } else if (lsame_(cmach, "L")) {
	rmach = emax;
    } else if (lsame_(cmach, "O")) {
	rmach = rmax;
    }

    ret_val = rmach;
    return ret_val;

/*     End of SLAMCH */

} /* slamch_ */



/* Subroutine */ int slamc1_(integer *beta, integer *t, logical *rnd, logical 
	*ieee1)
{
/*  -- LAPACK auxiliary routine (version 3.0) --
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLAMC1 determines the machine parameters given by BETA, T, RND, and   
    IEEE1.   

    Arguments   
    =========   

    BETA    (output) INTEGER   
            The base of the machine.   

    T       (output) INTEGER   
            The number of ( BETA ) digits in the mantissa.   

    RND     (output) LOGICAL   
            Specifies whether proper rounding  ( RND = .TRUE. )  or   
            chopping  ( RND = .FALSE. )  occurs in addition. This may not 
  
            be a reliable guide to the way in which the machine performs 
  
            its arithmetic.   

    IEEE1   (output) LOGICAL   
            Specifies whether rounding appears to be done in the IEEE   
            'round to nearest' style.   

    Further Details   
    ===============   

    The routine is based on the routine  ENVRON  by Malcolm and   
    incorporates suggestions by Gentleman and Marovich. See   

       Malcolm M. A. (1972) Algorithms to reveal properties of   
          floating-point arithmetic. Comms. of the ACM, 15, 949-951.   

       Gentleman W. M. and Marovich S. B. (1974) More on algorithms   
          that reveal properties of floating point arithmetic units.   
          Comms. of the ACM, 17, 276-277.   

   ===================================================================== 
*/
    /* Initialized data */
    static logical first = TRUE_;
    /* System generated locals */
    real r__1, r__2;
    /* Local variables */
    static logical lrnd;
    static real a, b, c, f;
    static integer lbeta;
    static real savec;
    static logical lieee1;
    static real t1, t2;
    extern doublereal slamc3_(real *, real *);
    static integer lt;
    static real one, qtr;



    if (first) {
	first = FALSE_;
	one = 1.f;

/*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BE
TA,   
          IEEE1, T and RND.   

          Throughout this routine  we use the function  SLAMC3  to ens
ure   
          that relevant values are  stored and not held in registers, 
 or   
          are not affected by optimizers.   

          Compute  a = 2.0**m  with the  smallest positive integer m s
uch   
          that   

             fl( a + 1.0 ) = a. */

	a = 1.f;
	c = 1.f;

/* +       WHILE( C.EQ.ONE )LOOP */
L10:
	if (c == one) {
	    a *= 2;
	    c = slamc3_(&a, &one);
	    r__1 = -(doublereal)a;
	    c = slamc3_(&c, &r__1);
	    goto L10;
	}
/* +       END WHILE   

          Now compute  b = 2.0**m  with the smallest positive integer 
m   
          such that   

             fl( a + b ) .gt. a. */

	b = 1.f;
	c = slamc3_(&a, &b);

/* +       WHILE( C.EQ.A )LOOP */
L20:
	if (c == a) {
	    b *= 2;
	    c = slamc3_(&a, &b);
	    goto L20;
	}
/* +       END WHILE   

          Now compute the base.  a and c  are neighbouring floating po
int   
          numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and
 so   
          their difference is beta. Adding 0.25 to c is to ensure that
 it   
          is truncated to beta and not ( beta - 1 ). */

	qtr = one / 4;
	savec = c;
	r__1 = -(doublereal)a;
	c = slamc3_(&c, &r__1);
	lbeta = c + qtr;

/*        Now determine whether rounding or chopping occurs,  by addin
g a   
          bit  less  than  beta/2  and a  bit  more  than  beta/2  to 
 a. */

	b = (real) lbeta;
	r__1 = b / 2;
	r__2 = -(doublereal)b / 100;
	f = slamc3_(&r__1, &r__2);
	c = slamc3_(&f, &a);
	if (c == a) {
	    lrnd = TRUE_;
	} else {
	    lrnd = FALSE_;
	}
	r__1 = b / 2;
	r__2 = b / 100;
	f = slamc3_(&r__1, &r__2);
	c = slamc3_(&f, &a);
	if (lrnd && c == a) {
	    lrnd = FALSE_;
	}

/*        Try and decide whether rounding is done in the  IEEE  'round
 to   
          nearest' style. B/2 is half a unit in the last place of the 
two   
          numbers A and SAVEC. Furthermore, A is even, i.e. has last  
bit   
          zero, and SAVEC is odd. Thus adding B/2 to A should not  cha
nge   
          A, but adding B/2 to SAVEC should change SAVEC. */

	r__1 = b / 2;
	t1 = slamc3_(&r__1, &a);
	r__1 = b / 2;
	t2 = slamc3_(&r__1, &savec);
	lieee1 = t1 == a && t2 > savec && lrnd;

/*        Now find  the  mantissa, t.  It should  be the  integer part
 of   
          log to the base beta of a,  however it is safer to determine
  t   
          by powering.  So we find t as the smallest positive integer 
for   
          which   

             fl( beta**t + 1.0 ) = 1.0. */

	lt = 0;
	a = 1.f;
	c = 1.f;

/* +       WHILE( C.EQ.ONE )LOOP */
L30:
	if (c == one) {
	    ++lt;
	    a *= lbeta;
	    c = slamc3_(&a, &one);
	    r__1 = -(doublereal)a;
	    c = slamc3_(&c, &r__1);
	    goto L30;
	}
/* +       END WHILE */

    }

    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *ieee1 = lieee1;
    return 0;

/*     End of SLAMC1 */

} /* slamc1_ */



/* Subroutine */ int slamc2_(integer *beta, integer *t, logical *rnd, real *
	eps, integer *emin, real *rmin, integer *emax, real *rmax)
{
/*  -- LAPACK auxiliary routine (version 3.0) --
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLAMC2 determines the machine parameters specified in its argument   
    list.   

    Arguments   
    =========   

    BETA    (output) INTEGER   
            The base of the machine.   

    T       (output) INTEGER   
            The number of ( BETA ) digits in the mantissa.   

    RND     (output) LOGICAL   
            Specifies whether proper rounding  ( RND = .TRUE. )  or   
            chopping  ( RND = .FALSE. )  occurs in addition. This may not 
  
            be a reliable guide to the way in which the machine performs 
  
            its arithmetic.   

    EPS     (output) REAL   
            The smallest positive number such that   

               fl( 1.0 - EPS ) .LT. 1.0,   

            where fl denotes the computed value.   

    EMIN    (output) INTEGER   
            The minimum exponent before (gradual) underflow occurs.   

    RMIN    (output) REAL   
            The smallest normalized number for the machine, given by   
            BASE**( EMIN - 1 ), where  BASE  is the floating point value 
  
            of BETA.   

    EMAX    (output) INTEGER   
            The maximum exponent before overflow occurs.   

    RMAX    (output) REAL   
            The largest positive number for the machine, given by   
            BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point 
  
            value of BETA.   

    Further Details   
    ===============   

    The computation of  EPS  is based on a routine PARANOIA by   
    W. Kahan of the University of California at Berkeley.   

   ===================================================================== 
*/
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* Initialized data */
    static logical first = TRUE_;
    static logical iwarn = FALSE_;
    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3, r__4, r__5;
    /* Builtin functions */
    double pow_ri(real *, integer *);
    /* Local variables */
    static logical ieee;
    static real half;
    static logical lrnd;
    static real leps, zero, a, b, c;
    static integer i, lbeta;
    static real rbase;
    static integer lemin, lemax, gnmin;
    static real small;
    static integer gpmin;
    static real third, lrmin, lrmax, sixth;
    static logical lieee1;
    extern /* Subroutine */ int slamc1_(integer *, integer *, logical *, 
	    logical *);
    extern doublereal slamc3_(real *, real *);
    extern /* Subroutine */ int slamc4_(integer *, real *, integer *), 
	    slamc5_(integer *, integer *, integer *, logical *, integer *, 
	    real *);
    static integer lt, ngnmin, ngpmin;
    static real one, two;



    if (first) {
	first = FALSE_;
	zero = 0.f;
	one = 1.f;
	two = 2.f;

/*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values
 of   
          BETA, T, RND, EPS, EMIN and RMIN.   

          Throughout this routine  we use the function  SLAMC3  to ens
ure   
          that relevant values are stored  and not held in registers, 
 or   
          are not affected by optimizers.   

          SLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1. 
*/

	slamc1_(&lbeta, &lt, &lrnd, &lieee1);

/*        Start to find EPS. */

	b = (real) lbeta;
	i__1 = -lt;
	a = pow_ri(&b, &i__1);
	leps = a;

/*        Try some tricks to see whether or not this is the correct  E
PS. */

	b = two / 3;
	half = one / 2;
	r__1 = -(doublereal)half;
	sixth = slamc3_(&b, &r__1);
	third = slamc3_(&sixth, &sixth);
	r__1 = -(doublereal)half;
	b = slamc3_(&third, &r__1);
	b = slamc3_(&b, &sixth);
	b = dabs(b);
	if (b < leps) {
	    b = leps;
	}

	leps = 1.f;

/* +       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP */
L10:
	if (leps > b && b > zero) {
	    leps = b;
	    r__1 = half * leps;
/* Computing 5th power */
	    r__3 = two, r__4 = r__3, r__3 *= r__3;
/* Computing 2nd power */
	    r__5 = leps;
	    r__2 = r__4 * (r__3 * r__3) * (r__5 * r__5);
	    c = slamc3_(&r__1, &r__2);
	    r__1 = -(doublereal)c;
	    c = slamc3_(&half, &r__1);
	    b = slamc3_(&half, &c);
	    r__1 = -(doublereal)b;
	    c = slamc3_(&half, &r__1);
	    b = slamc3_(&half, &c);
	    goto L10;
	}
/* +       END WHILE */

	if (a < leps) {
	    leps = a;
	}

/*        Computation of EPS complete.   

          Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3
)).   
          Keep dividing  A by BETA until (gradual) underflow occurs. T
his   
          is detected when we cannot recover the previous A. */

	rbase = one / lbeta;
	small = one;
	for (i = 1; i <= 3; ++i) {
	    r__1 = small * rbase;
	    small = slamc3_(&r__1, &zero);
/* L20: */
	}
	a = slamc3_(&one, &small);
	slamc4_(&ngpmin, &one, &lbeta);
	r__1 = -(doublereal)one;
	slamc4_(&ngnmin, &r__1, &lbeta);
	slamc4_(&gpmin, &a, &lbeta);
	r__1 = -(doublereal)a;
	slamc4_(&gnmin, &r__1, &lbeta);
	ieee = FALSE_;

	if (ngpmin == ngnmin && gpmin == gnmin) {
	    if (ngpmin == gpmin) {
		lemin = ngpmin;
/*            ( Non twos-complement machines, no gradual under
flow;   
                e.g.,  VAX ) */
	    } else if (gpmin - ngpmin == 3) {
		lemin = ngpmin - 1 + lt;
		ieee = TRUE_;
/*            ( Non twos-complement machines, with gradual und
erflow;   
                e.g., IEEE standard followers ) */
	    } else {
		lemin = f2cmin(ngpmin,gpmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else if (ngpmin == gpmin && ngnmin == gnmin) {
	    if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1) {
		lemin = f2cmax(ngpmin,ngnmin);
/*            ( Twos-complement machines, no gradual underflow
;   
                e.g., CYBER 205 ) */
	    } else {
		lemin = f2cmin(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1 && gpmin == gnmin)
		 {
	    if (gpmin - f2cmin(ngpmin,ngnmin) == 3) {
		lemin = f2cmax(ngpmin,ngnmin) - 1 + lt;
/*            ( Twos-complement machines with gradual underflo
w;   
                no known machine ) */
	    } else {
		lemin = f2cmin(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else {
/* Computing MIN */
	    i__1 = f2cmin(ngpmin,ngnmin), i__1 = f2cmin(i__1,gpmin);
	    lemin = f2cmin(i__1,gnmin);
/*         ( A guess; no known machine ) */
	    iwarn = TRUE_;
	}
/* **   
   Comment out this if block if EMIN is ok */
	if (iwarn) {
	    first = TRUE_;
	    printf("\n\n WARNING. The value EMIN may be incorrect:- ");
	    printf("EMIN = %8i\n",lemin);
	    printf("If, after inspection, the value EMIN looks acceptable");
            printf("please comment out \n the IF block as marked within the"); 
            printf("code of routine SLAMC2, \n otherwise supply EMIN"); 
            printf("explicitly.\n");
	}
/* **   

          Assume IEEE arithmetic if we found denormalised  numbers abo
ve,   
          or if arithmetic seems to round in the  IEEE style,  determi
ned   
          in routine SLAMC1. A true IEEE machine should have both  thi
ngs   
          true; however, faulty machines may have one or the other. */

	ieee = ieee || lieee1;

/*        Compute  RMIN by successive division by  BETA. We could comp
ute   
          RMIN as BASE**( EMIN - 1 ),  but some machines underflow dur
ing   
          this computation. */

	lrmin = 1.f;
	i__1 = 1 - lemin;
	for (i = 1; i <= 1-lemin; ++i) {
	    r__1 = lrmin * rbase;
	    lrmin = slamc3_(&r__1, &zero);
/* L30: */
	}

/*        Finally, call SLAMC5 to compute EMAX and RMAX. */

	slamc5_(&lbeta, &lt, &lemin, &ieee, &lemax, &lrmax);
    }

    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *eps = leps;
    *emin = lemin;
    *rmin = lrmin;
    *emax = lemax;
    *rmax = lrmax;

    return 0;


/*     End of SLAMC2 */

} /* slamc2_ */



doublereal slamc3_(real *a, real *b)
{
/*  -- LAPACK auxiliary routine (version 3.0) --
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLAMC3  is intended to force  A  and  B  to be stored prior to doing 
  
    the addition of  A  and  B ,  for use in situations where optimizers 
  
    might hold one of these in a register.   

    Arguments   
    =========   

    A, B    (input) REAL   
            The values A and B.   

   ===================================================================== 
*/
/* >>Start of File<<   
       System generated locals */
    real ret_val;



    ret_val = *a + *b;

    return ret_val;

/*     End of SLAMC3 */

} /* slamc3_ */



/* Subroutine */ int slamc4_(integer *emin, real *start, integer *base)
{
/*  -- LAPACK auxiliary routine (version 3.0) --
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLAMC4 is a service routine for SLAMC2.   

    Arguments   
    =========   

    EMIN    (output) EMIN   
            The minimum exponent before (gradual) underflow, computed by 
  
            setting A = START and dividing by BASE until the previous A   
            can not be recovered.   

    START   (input) REAL   
            The starting point for determining EMIN.   

    BASE    (input) INTEGER   
            The base of the machine.   

   ===================================================================== 
*/
    /* System generated locals */
    integer i__1;
    real r__1;
    /* Local variables */
    static real zero, a;
    static integer i;
    static real rbase, b1, b2, c1, c2, d1, d2;
    extern doublereal slamc3_(real *, real *);
    static real one;



    a = *start;
    one = 1.f;
    rbase = one / *base;
    zero = 0.f;
    *emin = 1;
    r__1 = a * rbase;
    b1 = slamc3_(&r__1, &zero);
    c1 = a;
    c2 = a;
    d1 = a;
    d2 = a;
/* +    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.   
      $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP */
L10:
    if (c1 == a && c2 == a && d1 == a && d2 == a) {
	--(*emin);
	a = b1;
	r__1 = a / *base;
	b1 = slamc3_(&r__1, &zero);
	r__1 = b1 * *base;
	c1 = slamc3_(&r__1, &zero);
	d1 = zero;
	i__1 = *base;
	for (i = 1; i <= *base; ++i) {
	    d1 += b1;
/* L20: */
	}
	r__1 = a * rbase;
	b2 = slamc3_(&r__1, &zero);
	r__1 = b2 / rbase;
	c2 = slamc3_(&r__1, &zero);
	d2 = zero;
	i__1 = *base;
	for (i = 1; i <= *base; ++i) {
	    d2 += b2;
/* L30: */
	}
	goto L10;
    }
/* +    END WHILE */

    return 0;

/*     End of SLAMC4 */

} /* slamc4_ */



/* Subroutine */ int slamc5_(integer *beta, integer *p, integer *emin, 
	logical *ieee, integer *emax, real *rmax)
{
/*  -- LAPACK auxiliary routine (version 3.0) --
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLAMC5 attempts to compute RMAX, the largest machine floating-point   
    number, without overflow.  It assumes that EMAX + abs(EMIN) sum   
    approximately to a power of 2.  It will fail on machines where this   
    assumption does not hold, for example, the Cyber 205 (EMIN = -28625, 
  
    EMAX = 28718).  It will also fail if the value supplied for EMIN is   
    too large (i.e. too close to zero), probably with overflow.   

    Arguments   
    =========   

    BETA    (input) INTEGER   
            The base of floating-point arithmetic.   

    P       (input) INTEGER   
            The number of base BETA digits in the mantissa of a   
            floating-point value.   

    EMIN    (input) INTEGER   
            The minimum exponent before (gradual) underflow.   

    IEEE    (input) LOGICAL   
            A logical flag specifying whether or not the arithmetic   
            system is thought to comply with the IEEE standard.   

    EMAX    (output) INTEGER   
            The largest exponent before overflow   

    RMAX    (output) REAL   
            The largest machine floating-point number.   

   ===================================================================== 
  


       First compute LEXP and UEXP, two powers of 2 that bound   
       abs(EMIN). We then assume that EMAX + abs(EMIN) will sum   
       approximately to the bound that is closest to abs(EMIN).   
       (EMAX is the exponent of the required number RMAX). */
    /* Table of constant values */
    static real c_b5 = 0.f;
    
    /* System generated locals */
    integer i__1;
    real r__1;
    /* Local variables */
    static integer lexp;
    static real oldy;
    static integer uexp, i;
    static real y, z;
    static integer nbits;
    extern doublereal slamc3_(real *, real *);
    static real recbas;
    static integer exbits, expsum, try__;



    lexp = 1;
    exbits = 1;
L10:
    try__ = lexp << 1;
    if (try__ <= -(*emin)) {
	lexp = try__;
	++exbits;
	goto L10;
    }
    if (lexp == -(*emin)) {
	uexp = lexp;
    } else {
	uexp = try__;
	++exbits;
    }

/*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater   
       than or equal to EMIN. EXBITS is the number of bits needed to   
       store the exponent. */

    if (uexp + *emin > -lexp - *emin) {
	expsum = lexp << 1;
    } else {
	expsum = uexp << 1;
    }

/*     EXPSUM is the exponent range, approximately equal to   
       EMAX - EMIN + 1 . */

    *emax = expsum + *emin - 1;
    nbits = exbits + 1 + *p;

/*     NBITS is the total number of bits needed to store a   
       floating-point number. */

    if (nbits % 2 == 1 && *beta == 2) {

/*        Either there are an odd number of bits used to store a   
          floating-point number, which is unlikely, or some bits are 
  
          not used in the representation of numbers, which is possible
,   
          (e.g. Cray machines) or the mantissa has an implicit bit,   
          (e.g. IEEE machines, Dec Vax machines), which is perhaps the
   
          most likely. We have to assume the last alternative.   
          If this is true, then we need to reduce EMAX by one because 
  
          there must be some way of representing zero in an implicit-b
it   
          system. On machines like Cray, we are reducing EMAX by one 
  
          unnecessarily. */

	--(*emax);
    }

    if (*ieee) {

/*        Assume we are on an IEEE machine which reserves one exponent
   
          for infinity and NaN. */

	--(*emax);
    }

/*     Now create RMAX, the largest machine number, which should   
       be equal to (1.0 - BETA**(-P)) * BETA**EMAX .   

       First compute 1.0 - BETA**(-P), being careful that the   
       result is less than 1.0 . */

    recbas = 1.f / *beta;
    z = *beta - 1.f;
    y = 0.f;
    i__1 = *p;
    for (i = 1; i <= *p; ++i) {
	z *= recbas;
	if (y < 1.f) {
	    oldy = y;
	}
	y = slamc3_(&y, &z);
/* L20: */
    }
    if (y >= 1.f) {
	y = oldy;
    }

/*     Now multiply by BETA**EMAX to get RMAX. */

    i__1 = *emax;
    for (i = 1; i <= *emax; ++i) {
	r__1 = y * *beta;
	y = slamc3_(&r__1, &c_b5);
/* L30: */
    }

    *rmax = y;
    return 0;

/*     End of SLAMC5 */

} /* slamc5_ */




doublereal slanst_(char *norm, integer *n, real *d__, real *e)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    SLANST  returns the value of the one norm,  or the Frobenius norm, or   
    the  infinity norm,  or the  element of  largest absolute value  of a   
    real symmetric tridiagonal matrix A.   

    Description   
    ===========   

    SLANST returns the value   

       SLANST = ( f2cmax(abs(A(i,j))), NORM = 'M' or 'm'   
                (   
                ( norm1(A),         NORM = '1', 'O' or 'o'   
                (   
                ( normI(A),         NORM = 'I' or 'i'   
                (   
                ( normF(A),         NORM = 'F', 'f', 'E' or 'e'   

    where  norm1  denotes the  one norm of a matrix (maximum column sum),   
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and   
    normF  denotes the  Frobenius norm of a matrix (square root of sum of   
    squares).  Note that  f2cmax(abs(A(i,j)))  is not a  matrix norm.   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies the value to be returned in SLANST as described   
            above.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.  When N = 0, SLANST is   
            set to zero.   

    D       (input) REAL array, dimension (N)   
            The diagonal elements of A.   

    E       (input) REAL array, dimension (N-1)   
            The (n-1) sub-diagonal or super-diagonal elements of A.   

    =====================================================================   


       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer i__1;
    real ret_val, r__1, r__2, r__3, r__4, r__5;
    /* Builtin functions */
//    double sqrt(doublereal);
    /* Local variables */
    static integer i__;
    static real scale;
    extern logical lsame_(char *, char *);
    static real anorm;
    extern /* Subroutine */ int slassq_(integer *, real *, integer *, real *, 
	    real *);
    static real sum;


    --e;
    --d__;

    /* Function Body */
    if (*n <= 0) {
	anorm = 0.f;
    } else if (lsame_(norm, "M")) {

/*        Find f2cmax(abs(A(i,j))). */

	anorm = (r__1 = d__[*n], dabs(r__1));
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    r__2 = anorm, r__3 = (r__1 = d__[i__], dabs(r__1));
	    anorm = df2cmax(r__2,r__3);
/* Computing MAX */
	    r__2 = anorm, r__3 = (r__1 = e[i__], dabs(r__1));
	    anorm = df2cmax(r__2,r__3);
/* L10: */
	}
    } else if (lsame_(norm, "O") || *(unsigned char *)
	    norm == '1' || lsame_(norm, "I")) {

/*        Find norm1(A). */

	if (*n == 1) {
	    anorm = dabs(d__[1]);
	} else {
/* Computing MAX */
	    r__3 = dabs(d__[1]) + dabs(e[1]), r__4 = (r__1 = e[*n - 1], dabs(
		    r__1)) + (r__2 = d__[*n], dabs(r__2));
	    anorm = df2cmax(r__3,r__4);
	    i__1 = *n - 1;
	    for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing MAX */
		r__4 = anorm, r__5 = (r__1 = d__[i__], dabs(r__1)) + (r__2 = 
			e[i__], dabs(r__2)) + (r__3 = e[i__ - 1], dabs(r__3));
		anorm = df2cmax(r__4,r__5);
/* L20: */
	    }
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */

	scale = 0.f;
	sum = 1.f;
	if (*n > 1) {
	    i__1 = *n - 1;
	    slassq_(&i__1, &e[1], &c__1, &scale, &sum);
	    sum *= 2;
	}
	slassq_(n, &d__[1], &c__1, &scale, &sum);
	anorm = scale * sqrt(sum);
    }

    ret_val = anorm;
    return ret_val;

/*     End of SLANST */

} /* slanst_ */




doublereal slansy_(char *norm, char *uplo, integer *n, real *a, integer *lda, 
	real *work)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLANSY  returns the value of the one norm,  or the Frobenius norm, or   
    the  infinity norm,  or the  element of  largest absolute value  of a   
    real symmetric matrix A.   

    Description   
    ===========   

    SLANSY returns the value   

       SLANSY = ( f2cmax(abs(A(i,j))), NORM = 'M' or 'm'   
                (   
                ( norm1(A),         NORM = '1', 'O' or 'o'   
                (   
                ( normI(A),         NORM = 'I' or 'i'   
                (   
                ( normF(A),         NORM = 'F', 'f', 'E' or 'e'   

    where  norm1  denotes the  one norm of a matrix (maximum column sum),   
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and   
    normF  denotes the  Frobenius norm of a matrix (square root of sum of   
    squares).  Note that  f2cmax(abs(A(i,j)))  is not a  matrix norm.   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies the value to be returned in SLANSY as described   
            above.   

    UPLO    (input) CHARACTER*1   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is to be referenced.   
            = 'U':  Upper triangular part of A is referenced   
            = 'L':  Lower triangular part of A is referenced   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.  When N = 0, SLANSY is   
            set to zero.   

    A       (input) REAL array, dimension (LDA,N)   
            The symmetric matrix A.  If UPLO = 'U', the leading n by n   
            upper triangular part of A contains the upper triangular part   
            of the matrix A, and the strictly lower triangular part of A   
            is not referenced.  If UPLO = 'L', the leading n by n lower   
            triangular part of A contains the lower triangular part of   
            the matrix A, and the strictly upper triangular part of A is   
            not referenced.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= f2cmax(N,1).   

    WORK    (workspace) REAL array, dimension (LWORK),   
            where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,   
            WORK is not referenced.   

   =====================================================================   


       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    real ret_val, r__1, r__2, r__3;
    /* Builtin functions */
//    double sqrt(doublereal);
    /* Local variables */
    static real absa;
    static integer i__, j;
    static real scale;
    extern logical lsame_(char *, char *);
    static real value;
    extern /* Subroutine */ int slassq_(integer *, real *, integer *, real *, 
	    real *);
    static real sum;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --work;

    /* Function Body */
    if (*n == 0) {
	value = 0.f;
    } else if (lsame_(norm, "M")) {

/*        Find f2cmax(abs(A(i,j))). */

	value = 0.f;
	if (lsame_(uplo, "U")) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
		    r__2 = value, r__3 = (r__1 = a_ref(i__, j), dabs(r__1));
		    value = df2cmax(r__2,r__3);
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = j; i__ <= i__2; ++i__) {
/* Computing MAX */
		    r__2 = value, r__3 = (r__1 = a_ref(i__, j), dabs(r__1));
		    value = df2cmax(r__2,r__3);
/* L30: */
		}
/* L40: */
	    }
	}
    } else if (lsame_(norm, "I") || lsame_(norm, "O") || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is symmetric). */

	value = 0.f;
	if (lsame_(uplo, "U")) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		sum = 0.f;
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    absa = (r__1 = a_ref(i__, j), dabs(r__1));
		    sum += absa;
		    work[i__] += absa;
/* L50: */
		}
		work[j] = sum + (r__1 = a_ref(j, j), dabs(r__1));
/* L60: */
	    }
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
		r__1 = value, r__2 = work[i__];
		value = df2cmax(r__1,r__2);
/* L70: */
	    }
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		work[i__] = 0.f;
/* L80: */
	    }
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		sum = work[j] + (r__1 = a_ref(j, j), dabs(r__1));
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    absa = (r__1 = a_ref(i__, j), dabs(r__1));
		    sum += absa;
		    work[i__] += absa;
/* L90: */
		}
		value = df2cmax(value,sum);
/* L100: */
	    }
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */

	scale = 0.f;
	sum = 1.f;
	if (lsame_(uplo, "U")) {
	    i__1 = *n;
	    for (j = 2; j <= i__1; ++j) {
		i__2 = j - 1;
		slassq_(&i__2, &a_ref(1, j), &c__1, &scale, &sum);
/* L110: */
	    }
	} else {
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n - j;
		slassq_(&i__2, &a_ref(j + 1, j), &c__1, &scale, &sum);
/* L120: */
	    }
	}
	sum *= 2;
	i__1 = *lda + 1;
	slassq_(n, &a[a_offset], &i__1, &scale, &sum);
	value = scale * sqrt(sum);
    }

    ret_val = value;
    return ret_val;

/*     End of SLANSY */

} /* slansy_ */

#undef a_ref





doublereal slapy2_(real *x, real *y)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary   
    overflow.   

    Arguments   
    =========   

    X       (input) REAL   
    Y       (input) REAL   
            X and Y specify the values x and y.   

    ===================================================================== */
    /* System generated locals */
    real ret_val, r__1;
    /* Builtin functions */
//    double sqrt(doublereal);
    /* Local variables */
    static real xabs, yabs, w, z__;



    xabs = dabs(*x);
    yabs = dabs(*y);
    w = df2cmax(xabs,yabs);
    z__ = df2cmin(xabs,yabs);
    if (z__ == 0.f) {
	ret_val = w;
    } else {
/* Computing 2nd power */
	r__1 = z__ / w;
	ret_val = w * sqrt(r__1 * r__1 + 1.f);
    }
    return ret_val;

/*     End of SLAPY2 */

} /* slapy2_ */




/* Subroutine */ int slarfb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, real *v, integer *ldv, 
	real *t, integer *ldt, real *c__, integer *ldc, real *work, integer *
	ldwork)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    SLARFB applies a real block reflector H or its transpose H' to a   
    real m by n matrix C, from either the left or the right.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': apply H or H' from the Left   
            = 'R': apply H or H' from the Right   

    TRANS   (input) CHARACTER*1   
            = 'N': apply H (No transpose)   
            = 'T': apply H' (Transpose)   

    DIRECT  (input) CHARACTER*1   
            Indicates how H is formed from a product of elementary   
            reflectors   
            = 'F': H = H(1) H(2) . . . H(k) (Forward)   
            = 'B': H = H(k) . . . H(2) H(1) (Backward)   

    STOREV  (input) CHARACTER*1   
            Indicates how the vectors which define the elementary   
            reflectors are stored:   
            = 'C': Columnwise   
            = 'R': Rowwise   

    M       (input) INTEGER   
            The number of rows of the matrix C.   

    N       (input) INTEGER   
            The number of columns of the matrix C.   

    K       (input) INTEGER   
            The order of the matrix T (= the number of elementary   
            reflectors whose product defines the block reflector).   

    V       (input) REAL array, dimension   
                                  (LDV,K) if STOREV = 'C'   
                                  (LDV,M) if STOREV = 'R' and SIDE = 'L'   
                                  (LDV,N) if STOREV = 'R' and SIDE = 'R'   
            The matrix V. See further details.   

    LDV     (input) INTEGER   
            The leading dimension of the array V.   
            If STOREV = 'C' and SIDE = 'L', LDV >= f2cmax(1,M);   
            if STOREV = 'C' and SIDE = 'R', LDV >= f2cmax(1,N);   
            if STOREV = 'R', LDV >= K.   

    T       (input) REAL array, dimension (LDT,K)   
            The triangular k by k matrix T in the representation of the   
            block reflector.   

    LDT     (input) INTEGER   
            The leading dimension of the array T. LDT >= K.   

    C       (input/output) REAL array, dimension (LDC,N)   
            On entry, the m by n matrix C.   
            On exit, C is overwritten by H*C or H'*C or C*H or C*H'.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDA >= f2cmax(1,M).   

    WORK    (workspace) REAL array, dimension (LDWORK,K)   

    LDWORK  (input) INTEGER   
            The leading dimension of the array WORK.   
            If SIDE = 'L', LDWORK >= f2cmax(1,N);   
            if SIDE = 'R', LDWORK >= f2cmax(1,M).   

    =====================================================================   


       Quick return if possible   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    static real c_b14 = 1.f;
    static real c_b25 = -1.f;
    
    /* System generated locals */
    integer c_dim1, c_offset, t_dim1, t_offset, v_dim1, v_offset, work_dim1, 
	    work_offset, i__1, i__2;
    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, real *, real *, integer *, real *, integer *, real *, 
	    real *, integer *), scopy_(integer *, real *, 
	    integer *, real *, integer *), strmm_(char *, char *, char *, 
	    char *, integer *, integer *, real *, real *, integer *, real *, 
	    integer *);
    static char transt[1];
#define work_ref(a_1,a_2) work[(a_2)*work_dim1 + a_1]
#define c___ref(a_1,a_2) c__[(a_2)*c_dim1 + a_1]
#define v_ref(a_1,a_2) v[(a_2)*v_dim1 + a_1]


    v_dim1 = *ldv;
    v_offset = 1 + v_dim1 * 1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1 * 1;
    t -= t_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;
    work_dim1 = *ldwork;
    work_offset = 1 + work_dim1 * 1;
    work -= work_offset;

    /* Function Body */
    if (*m <= 0 || *n <= 0) {
	return 0;
    }

    if (lsame_(trans, "N")) {
	*(unsigned char *)transt = 'T';
    } else {
	*(unsigned char *)transt = 'N';
    }

    if (lsame_(storev, "C")) {

	if (lsame_(direct, "F")) {

/*           Let  V =  ( V1 )    (first K rows)   
                       ( V2 )   
             where  V1  is unit lower triangular. */

	    if (lsame_(side, "L")) {

/*              Form  H * C  or  H' * C  where  C = ( C1 )   
                                                    ( C2 )   

                W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)   

                W := C1' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    scopy_(n, &c___ref(j, 1), ldc, &work_ref(1, j), &c__1);
/* L10: */
		}

/*              W := W * V1 */

		strmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b14,
			 &v[v_offset], ldv, &work[work_offset], ldwork);
		if (*m > *k) {

/*                 W := W + C2'*V2 */

		    i__1 = *m - *k;
		    sgemm_("Transpose", "No transpose", n, k, &i__1, &c_b14, &
			    c___ref(*k + 1, 1), ldc, &v_ref(*k + 1, 1), ldv, &
			    c_b14, &work[work_offset], ldwork);
		}

/*              W := W * T'  or  W * T */

		strmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - V * W' */

		if (*m > *k) {

/*                 C2 := C2 - V2 * W' */

		    i__1 = *m - *k;
		    sgemm_("No transpose", "Transpose", &i__1, n, k, &c_b25, &
			    v_ref(*k + 1, 1), ldv, &work[work_offset], ldwork,
			     &c_b14, &c___ref(*k + 1, 1), ldc);
		}

/*              W := W * V1' */

		strmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b14, &
			v[v_offset], ldv, &work[work_offset], ldwork);

/*              C1 := C1 - W' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(j, i__) = c___ref(j, i__) - work_ref(i__, j);
/* L20: */
		    }
/* L30: */
		}

	    } else if (lsame_(side, "R")) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )   

                W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)   

                W := C1 */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    scopy_(m, &c___ref(1, j), &c__1, &work_ref(1, j), &c__1);
/* L40: */
		}

/*              W := W * V1 */

		strmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b14,
			 &v[v_offset], ldv, &work[work_offset], ldwork);
		if (*n > *k) {

/*                 W := W + C2 * V2 */

		    i__1 = *n - *k;
		    sgemm_("No transpose", "No transpose", m, k, &i__1, &
			    c_b14, &c___ref(1, *k + 1), ldc, &v_ref(*k + 1, 1)
			    , ldv, &c_b14, &work[work_offset], ldwork);
		}

/*              W := W * T  or  W * T' */

		strmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - W * V' */

		if (*n > *k) {

/*                 C2 := C2 - W * V2' */

		    i__1 = *n - *k;
		    sgemm_("No transpose", "Transpose", m, &i__1, k, &c_b25, &
			    work[work_offset], ldwork, &v_ref(*k + 1, 1), ldv,
			     &c_b14, &c___ref(1, *k + 1), ldc);
		}

/*              W := W * V1' */

		strmm_("Right", "Lower", "Transpose", "Unit", m, k, &c_b14, &
			v[v_offset], ldv, &work[work_offset], ldwork);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = c___ref(i__, j) - work_ref(i__, j);
/* L50: */
		    }
/* L60: */
		}
	    }

	} else {

/*           Let  V =  ( V1 )   
                       ( V2 )    (last K rows)   
             where  V2  is unit upper triangular. */

	    if (lsame_(side, "L")) {

/*              Form  H * C  or  H' * C  where  C = ( C1 )   
                                                    ( C2 )   

                W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)   

                W := C2' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    scopy_(n, &c___ref(*m - *k + j, 1), ldc, &work_ref(1, j), 
			    &c__1);
/* L70: */
		}

/*              W := W * V2 */

		strmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b14,
			 &v_ref(*m - *k + 1, 1), ldv, &work[work_offset], 
			ldwork);
		if (*m > *k) {

/*                 W := W + C1'*V1 */

		    i__1 = *m - *k;
		    sgemm_("Transpose", "No transpose", n, k, &i__1, &c_b14, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &c_b14, &
			    work[work_offset], ldwork);
		}

/*              W := W * T'  or  W * T */

		strmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - V * W' */

		if (*m > *k) {

/*                 C1 := C1 - V1 * W' */

		    i__1 = *m - *k;
		    sgemm_("No transpose", "Transpose", &i__1, n, k, &c_b25, &
			    v[v_offset], ldv, &work[work_offset], ldwork, &
			    c_b14, &c__[c_offset], ldc)
			    ;
		}

/*              W := W * V2' */

		strmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b14, &
			v_ref(*m - *k + 1, 1), ldv, &work[work_offset], 
			ldwork);

/*              C2 := C2 - W' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(*m - *k + j, i__) = c___ref(*m - *k + j, i__) 
				- work_ref(i__, j);
/* L80: */
		    }
/* L90: */
		}

	    } else if (lsame_(side, "R")) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )   

                W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)   

                W := C2 */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    scopy_(m, &c___ref(1, *n - *k + j), &c__1, &work_ref(1, j)
			    , &c__1);
/* L100: */
		}

/*              W := W * V2 */

		strmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b14,
			 &v_ref(*n - *k + 1, 1), ldv, &work[work_offset], 
			ldwork);
		if (*n > *k) {

/*                 W := W + C1 * V1 */

		    i__1 = *n - *k;
		    sgemm_("No transpose", "No transpose", m, k, &i__1, &
			    c_b14, &c__[c_offset], ldc, &v[v_offset], ldv, &
			    c_b14, &work[work_offset], ldwork);
		}

/*              W := W * T  or  W * T' */

		strmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - W * V' */

		if (*n > *k) {

/*                 C1 := C1 - W * V1' */

		    i__1 = *n - *k;
		    sgemm_("No transpose", "Transpose", m, &i__1, k, &c_b25, &
			    work[work_offset], ldwork, &v[v_offset], ldv, &
			    c_b14, &c__[c_offset], ldc)
			    ;
		}

/*              W := W * V2' */

		strmm_("Right", "Upper", "Transpose", "Unit", m, k, &c_b14, &
			v_ref(*n - *k + 1, 1), ldv, &work[work_offset], 
			ldwork);

/*              C2 := C2 - W */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, *n - *k + j) = c___ref(i__, *n - *k + j) 
				- work_ref(i__, j);
/* L110: */
		    }
/* L120: */
		}
	    }
	}

    } else if (lsame_(storev, "R")) {

	if (lsame_(direct, "F")) {

/*           Let  V =  ( V1  V2 )    (V1: first K columns)   
             where  V1  is unit upper triangular. */

	    if (lsame_(side, "L")) {

/*              Form  H * C  or  H' * C  where  C = ( C1 )   
                                                    ( C2 )   

                W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)   

                W := C1' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    scopy_(n, &c___ref(j, 1), ldc, &work_ref(1, j), &c__1);
/* L130: */
		}

/*              W := W * V1' */

		strmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b14, &
			v[v_offset], ldv, &work[work_offset], ldwork);
		if (*m > *k) {

/*                 W := W + C2'*V2' */

		    i__1 = *m - *k;
		    sgemm_("Transpose", "Transpose", n, k, &i__1, &c_b14, &
			    c___ref(*k + 1, 1), ldc, &v_ref(1, *k + 1), ldv, &
			    c_b14, &work[work_offset], ldwork);
		}

/*              W := W * T'  or  W * T */

		strmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - V' * W' */

		if (*m > *k) {

/*                 C2 := C2 - V2' * W' */

		    i__1 = *m - *k;
		    sgemm_("Transpose", "Transpose", &i__1, n, k, &c_b25, &
			    v_ref(1, *k + 1), ldv, &work[work_offset], ldwork,
			     &c_b14, &c___ref(*k + 1, 1), ldc);
		}

/*              W := W * V1 */

		strmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b14,
			 &v[v_offset], ldv, &work[work_offset], ldwork);

/*              C1 := C1 - W' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(j, i__) = c___ref(j, i__) - work_ref(i__, j);
/* L140: */
		    }
/* L150: */
		}

	    } else if (lsame_(side, "R")) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )   

                W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)   

                W := C1 */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    scopy_(m, &c___ref(1, j), &c__1, &work_ref(1, j), &c__1);
/* L160: */
		}

/*              W := W * V1' */

		strmm_("Right", "Upper", "Transpose", "Unit", m, k, &c_b14, &
			v[v_offset], ldv, &work[work_offset], ldwork);
		if (*n > *k) {

/*                 W := W + C2 * V2' */

		    i__1 = *n - *k;
		    sgemm_("No transpose", "Transpose", m, k, &i__1, &c_b14, &
			    c___ref(1, *k + 1), ldc, &v_ref(1, *k + 1), ldv, &
			    c_b14, &work[work_offset], ldwork);
		}

/*              W := W * T  or  W * T' */

		strmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - W * V */

		if (*n > *k) {

/*                 C2 := C2 - W * V2 */

		    i__1 = *n - *k;
		    sgemm_("No transpose", "No transpose", m, &i__1, k, &
			    c_b25, &work[work_offset], ldwork, &v_ref(1, *k + 
			    1), ldv, &c_b14, &c___ref(1, *k + 1), ldc);
		}

/*              W := W * V1 */

		strmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b14,
			 &v[v_offset], ldv, &work[work_offset], ldwork);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = c___ref(i__, j) - work_ref(i__, j);
/* L170: */
		    }
/* L180: */
		}

	    }

	} else {

/*           Let  V =  ( V1  V2 )    (V2: last K columns)   
             where  V2  is unit lower triangular. */

	    if (lsame_(side, "L")) {

/*              Form  H * C  or  H' * C  where  C = ( C1 )   
                                                    ( C2 )   

                W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)   

                W := C2' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    scopy_(n, &c___ref(*m - *k + j, 1), ldc, &work_ref(1, j), 
			    &c__1);
/* L190: */
		}

/*              W := W * V2' */

		strmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b14, &
			v_ref(1, *m - *k + 1), ldv, &work[work_offset], 
			ldwork);
		if (*m > *k) {

/*                 W := W + C1'*V1' */

		    i__1 = *m - *k;
		    sgemm_("Transpose", "Transpose", n, k, &i__1, &c_b14, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &c_b14, &
			    work[work_offset], ldwork);
		}

/*              W := W * T'  or  W * T */

		strmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - V' * W' */

		if (*m > *k) {

/*                 C1 := C1 - V1' * W' */

		    i__1 = *m - *k;
		    sgemm_("Transpose", "Transpose", &i__1, n, k, &c_b25, &v[
			    v_offset], ldv, &work[work_offset], ldwork, &
			    c_b14, &c__[c_offset], ldc);
		}

/*              W := W * V2 */

		strmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b14,
			 &v_ref(1, *m - *k + 1), ldv, &work[work_offset], 
			ldwork);

/*              C2 := C2 - W' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(*m - *k + j, i__) = c___ref(*m - *k + j, i__) 
				- work_ref(i__, j);
/* L200: */
		    }
/* L210: */
		}

	    } else if (lsame_(side, "R")) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )   

                W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)   

                W := C2 */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    scopy_(m, &c___ref(1, *n - *k + j), &c__1, &work_ref(1, j)
			    , &c__1);
/* L220: */
		}

/*              W := W * V2' */

		strmm_("Right", "Lower", "Transpose", "Unit", m, k, &c_b14, &
			v_ref(1, *n - *k + 1), ldv, &work[work_offset], 
			ldwork);
		if (*n > *k) {

/*                 W := W + C1 * V1' */

		    i__1 = *n - *k;
		    sgemm_("No transpose", "Transpose", m, k, &i__1, &c_b14, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &c_b14, &
			    work[work_offset], ldwork);
		}

/*              W := W * T  or  W * T' */

		strmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - W * V */

		if (*n > *k) {

/*                 C1 := C1 - W * V1 */

		    i__1 = *n - *k;
		    sgemm_("No transpose", "No transpose", m, &i__1, k, &
			    c_b25, &work[work_offset], ldwork, &v[v_offset], 
			    ldv, &c_b14, &c__[c_offset], ldc);
		}

/*              W := W * V2 */

		strmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b14,
			 &v_ref(1, *n - *k + 1), ldv, &work[work_offset], 
			ldwork);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, *n - *k + j) = c___ref(i__, *n - *k + j) 
				- work_ref(i__, j);
/* L230: */
		    }
/* L240: */
		}

	    }

	}
    }

    return 0;

/*     End of SLARFB */

} /* slarfb_ */

#undef v_ref
#undef c___ref
#undef work_ref





/* Subroutine */ int slarf_(char *side, integer *m, integer *n, real *v, 
	integer *incv, real *tau, real *c__, integer *ldc, real *work)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    SLARF applies a real elementary reflector H to a real m by n matrix   
    C, from either the left or the right. H is represented in the form   

          H = I - tau * v * v'   

    where tau is a real scalar and v is a real vector.   

    If tau = 0, then H is taken to be the unit matrix.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': form  H * C   
            = 'R': form  C * H   

    M       (input) INTEGER   
            The number of rows of the matrix C.   

    N       (input) INTEGER   
            The number of columns of the matrix C.   

    V       (input) REAL array, dimension   
                       (1 + (M-1)*abs(INCV)) if SIDE = 'L'   
                    or (1 + (N-1)*abs(INCV)) if SIDE = 'R'   
            The vector v in the representation of H. V is not used if   
            TAU = 0.   

    INCV    (input) INTEGER   
            The increment between elements of v. INCV <> 0.   

    TAU     (input) REAL   
            The value tau in the representation of H.   

    C       (input/output) REAL array, dimension (LDC,N)   
            On entry, the m by n matrix C.   
            On exit, C is overwritten by the matrix H * C if SIDE = 'L',   
            or C * H if SIDE = 'R'.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= f2cmax(1,M).   

    WORK    (workspace) REAL array, dimension   
                           (N) if SIDE = 'L'   
                        or (M) if SIDE = 'R'   

    =====================================================================   


       Parameter adjustments */
    /* Table of constant values */
    static real c_b4 = 1.f;
    static real c_b5 = 0.f;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer c_dim1, c_offset;
    real r__1;
    /* Local variables */
    extern /* Subroutine */ int sger_(integer *, integer *, real *, real *, 
	    integer *, real *, integer *, real *, integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int sgemv_(char *, integer *, integer *, real *, 
	    real *, integer *, real *, integer *, real *, real *, integer *);


    --v;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;
    --work;

    /* Function Body */
    if (lsame_(side, "L")) {

/*        Form  H * C */

	if (*tau != 0.f) {

/*           w := C' * v */

	    sgemv_("Transpose", m, n, &c_b4, &c__[c_offset], ldc, &v[1], incv,
		     &c_b5, &work[1], &c__1);

/*           C := C - v * w' */

	    r__1 = -(*tau);
	    sger_(m, n, &r__1, &v[1], incv, &work[1], &c__1, &c__[c_offset], 
		    ldc);
	}
    } else {

/*        Form  C * H */

	if (*tau != 0.f) {

/*           w := C * v */

	    sgemv_("No transpose", m, n, &c_b4, &c__[c_offset], ldc, &v[1], 
		    incv, &c_b5, &work[1], &c__1);

/*           C := C - w * v' */

	    r__1 = -(*tau);
	    sger_(m, n, &r__1, &work[1], &c__1, &v[1], incv, &c__[c_offset], 
		    ldc);
	}
    }
    return 0;

/*     End of SLARF */

} /* slarf_ */




/* Subroutine */ int slarfg_(integer *n, real *alpha, real *x, integer *incx, 
	real *tau)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SLARFG generates a real elementary reflector H of order n, such   
    that   

          H * ( alpha ) = ( beta ),   H' * H = I.   
              (   x   )   (   0  )   

    where alpha and beta are scalars, and x is an (n-1)-element real   
    vector. H is represented in the form   

          H = I - tau * ( 1 ) * ( 1 v' ) ,   
                        ( v )   

    where tau is a real scalar and v is a real (n-1)-element   
    vector.   

    If the elements of x are all zero, then tau = 0 and H is taken to be   
    the unit matrix.   

    Otherwise  1 <= tau <= 2.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the elementary reflector.   

    ALPHA   (input/output) REAL   
            On entry, the value alpha.   
            On exit, it is overwritten with the value beta.   

    X       (input/output) REAL array, dimension   
                           (1+(N-2)*abs(INCX))   
            On entry, the vector x.   
            On exit, it is overwritten with the vector v.   

    INCX    (input) INTEGER   
            The increment between elements of X. INCX > 0.   

    TAU     (output) REAL   
            The value tau.   

    =====================================================================   


       Parameter adjustments */
    /* System generated locals */
    integer i__1;
    real r__1;
    /* Builtin functions */
    double r_sign(real *, real *);
    /* Local variables */
    static real beta;
    extern doublereal snrm2_(integer *, real *, integer *);
    static integer j;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    static real xnorm;
    extern doublereal slapy2_(real *, real *), slamch_(char *);
    static real safmin, rsafmn;
    static integer knt;

    --x;

    /* Function Body */
    if (*n <= 1) {
	*tau = 0.f;
	return 0;
    }

    i__1 = *n - 1;
    xnorm = snrm2_(&i__1, &x[1], incx);

    if (xnorm == 0.f) {

/*        H  =  I */

	*tau = 0.f;
    } else {

/*        general case */

	r__1 = slapy2_(alpha, &xnorm);
	beta = -r_sign(&r__1, alpha);
	safmin = slamch_("S") / slamch_("E");
	if (dabs(beta) < safmin) {

/*           XNORM, BETA may be inaccurate; scale X and recompute them */

	    rsafmn = 1.f / safmin;
	    knt = 0;
L10:
	    ++knt;
	    i__1 = *n - 1;
	    sscal_(&i__1, &rsafmn, &x[1], incx);
	    beta *= rsafmn;
	    *alpha *= rsafmn;
	    if (dabs(beta) < safmin) {
		goto L10;
	    }

/*           New BETA is at most 1, at least SAFMIN */

	    i__1 = *n - 1;
	    xnorm = snrm2_(&i__1, &x[1], incx);
	    r__1 = slapy2_(alpha, &xnorm);
	    beta = -r_sign(&r__1, alpha);
	    *tau = (beta - *alpha) / beta;
	    i__1 = *n - 1;
	    r__1 = 1.f / (*alpha - beta);
	    sscal_(&i__1, &r__1, &x[1], incx);

/*           If ALPHA is subnormal, it may lose relative accuracy */

	    *alpha = beta;
	    i__1 = knt;
	    for (j = 1; j <= i__1; ++j) {
		*alpha *= safmin;
/* L20: */
	    }
	} else {
	    *tau = (beta - *alpha) / beta;
	    i__1 = *n - 1;
	    r__1 = 1.f / (*alpha - beta);
	    sscal_(&i__1, &r__1, &x[1], incx);
	    *alpha = beta;
	}
    }

    return 0;

/*     End of SLARFG */

} /* slarfg_ */




/* Subroutine */ int slarft_(char *direct, char *storev, integer *n, integer *
	k, real *v, integer *ldv, real *tau, real *t, integer *ldt)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    SLARFT forms the triangular factor T of a real block reflector H   
    of order n, which is defined as a product of k elementary reflectors.   

    If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;   

    If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.   

    If STOREV = 'C', the vector which defines the elementary reflector   
    H(i) is stored in the i-th column of the array V, and   

       H  =  I - V * T * V'   

    If STOREV = 'R', the vector which defines the elementary reflector   
    H(i) is stored in the i-th row of the array V, and   

       H  =  I - V' * T * V   

    Arguments   
    =========   

    DIRECT  (input) CHARACTER*1   
            Specifies the order in which the elementary reflectors are   
            multiplied to form the block reflector:   
            = 'F': H = H(1) H(2) . . . H(k) (Forward)   
            = 'B': H = H(k) . . . H(2) H(1) (Backward)   

    STOREV  (input) CHARACTER*1   
            Specifies how the vectors which define the elementary   
            reflectors are stored (see also Further Details):   
            = 'C': columnwise   
            = 'R': rowwise   

    N       (input) INTEGER   
            The order of the block reflector H. N >= 0.   

    K       (input) INTEGER   
            The order of the triangular factor T (= the number of   
            elementary reflectors). K >= 1.   

    V       (input/output) REAL array, dimension   
                                 (LDV,K) if STOREV = 'C'   
                                 (LDV,N) if STOREV = 'R'   
            The matrix V. See further details.   

    LDV     (input) INTEGER   
            The leading dimension of the array V.   
            If STOREV = 'C', LDV >= f2cmax(1,N); if STOREV = 'R', LDV >= K.   

    TAU     (input) REAL array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i).   

    T       (output) REAL array, dimension (LDT,K)   
            The k by k triangular factor T of the block reflector.   
            If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is   
            lower triangular. The rest of the array is not used.   

    LDT     (input) INTEGER   
            The leading dimension of the array T. LDT >= K.   

    Further Details   
    ===============   

    The shape of the matrix V and the storage of the vectors which define   
    the H(i) is best illustrated by the following example with n = 5 and   
    k = 3. The elements equal to 1 are not stored; the corresponding   
    array elements are modified but restored on exit. The rest of the   
    array is not used.   

    DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':   

                 V = (  1       )                 V = (  1 v1 v1 v1 v1 )   
                     ( v1  1    )                     (     1 v2 v2 v2 )   
                     ( v1 v2  1 )                     (        1 v3 v3 )   
                     ( v1 v2 v3 )   
                     ( v1 v2 v3 )   

    DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':   

                 V = ( v1 v2 v3 )                 V = ( v1 v1  1       )   
                     ( v1 v2 v3 )                     ( v2 v2 v2  1    )   
                     (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )   
                     (     1 v3 )   
                     (        1 )   

    =====================================================================   


       Quick return if possible   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    static real c_b8 = 0.f;
    
    /* System generated locals */
    integer t_dim1, t_offset, v_dim1, v_offset, i__1, i__2, i__3;
    real r__1;
    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int sgemv_(char *, integer *, integer *, real *, 
	    real *, integer *, real *, integer *, real *, real *, integer *), strmv_(char *, char *, char *, integer *, real *, 
	    integer *, real *, integer *);
    static real vii;
#define t_ref(a_1,a_2) t[(a_2)*t_dim1 + a_1]
#define v_ref(a_1,a_2) v[(a_2)*v_dim1 + a_1]


    v_dim1 = *ldv;
    v_offset = 1 + v_dim1 * 1;
    v -= v_offset;
    --tau;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1 * 1;
    t -= t_offset;

    /* Function Body */
    if (*n == 0) {
	return 0;
    }

    if (lsame_(direct, "F")) {
	i__1 = *k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (tau[i__] == 0.f) {

/*              H(i)  =  I */

		i__2 = i__;
		for (j = 1; j <= i__2; ++j) {
		    t_ref(j, i__) = 0.f;
/* L10: */
		}
	    } else {

/*              general case */

		vii = v_ref(i__, i__);
		v_ref(i__, i__) = 1.f;
		if (lsame_(storev, "C")) {

/*                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i) */

		    i__2 = *n - i__ + 1;
		    i__3 = i__ - 1;
		    r__1 = -tau[i__];
		    sgemv_("Transpose", &i__2, &i__3, &r__1, &v_ref(i__, 1), 
			    ldv, &v_ref(i__, i__), &c__1, &c_b8, &t_ref(1, 
			    i__), &c__1);
		} else {

/*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)' */

		    i__2 = i__ - 1;
		    i__3 = *n - i__ + 1;
		    r__1 = -tau[i__];
		    sgemv_("No transpose", &i__2, &i__3, &r__1, &v_ref(1, i__)
			    , ldv, &v_ref(i__, i__), ldv, &c_b8, &t_ref(1, 
			    i__), &c__1);
		}
		v_ref(i__, i__) = vii;

/*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i) */

		i__2 = i__ - 1;
		strmv_("Upper", "No transpose", "Non-unit", &i__2, &t[
			t_offset], ldt, &t_ref(1, i__), &c__1);
		t_ref(i__, i__) = tau[i__];
	    }
/* L20: */
	}
    } else {
	for (i__ = *k; i__ >= 1; --i__) {
	    if (tau[i__] == 0.f) {

/*              H(i)  =  I */

		i__1 = *k;
		for (j = i__; j <= i__1; ++j) {
		    t_ref(j, i__) = 0.f;
/* L30: */
		}
	    } else {

/*              general case */

		if (i__ < *k) {
		    if (lsame_(storev, "C")) {
			vii = v_ref(*n - *k + i__, i__);
			v_ref(*n - *k + i__, i__) = 1.f;

/*                    T(i+1:k,i) :=   
                              - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i) */

			i__1 = *n - *k + i__;
			i__2 = *k - i__;
			r__1 = -tau[i__];
			sgemv_("Transpose", &i__1, &i__2, &r__1, &v_ref(1, 
				i__ + 1), ldv, &v_ref(1, i__), &c__1, &c_b8, &
				t_ref(i__ + 1, i__), &c__1);
			v_ref(*n - *k + i__, i__) = vii;
		    } else {
			vii = v_ref(i__, *n - *k + i__);
			v_ref(i__, *n - *k + i__) = 1.f;

/*                    T(i+1:k,i) :=   
                              - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)' */

			i__1 = *k - i__;
			i__2 = *n - *k + i__;
			r__1 = -tau[i__];
			sgemv_("No transpose", &i__1, &i__2, &r__1, &v_ref(
				i__ + 1, 1), ldv, &v_ref(i__, 1), ldv, &c_b8, 
				&t_ref(i__ + 1, i__), &c__1);
			v_ref(i__, *n - *k + i__) = vii;
		    }

/*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i) */

		    i__1 = *k - i__;
		    strmv_("Lower", "No transpose", "Non-unit", &i__1, &t_ref(
			    i__ + 1, i__ + 1), ldt, &t_ref(i__ + 1, i__), &
			    c__1);
		}
		t_ref(i__, i__) = tau[i__];
	    }
/* L40: */
	}
    }
    return 0;

/*     End of SLARFT */

} /* slarft_ */

#undef v_ref
#undef t_ref





/* Subroutine */ int slartg_(real *f, real *g, real *cs, real *sn, real *r__)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SLARTG generate a plane rotation so that   

       [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.   
       [ -SN  CS  ]     [ G ]     [ 0 ]   

    This is a slower, more accurate version of the BLAS1 routine SROTG,   
    with the following other differences:   
       F and G are unchanged on return.   
       If G=0, then CS=1 and SN=0.   
       If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any   
          floating point operations (saves work in SBDSQR when   
          there are zeros on the diagonal).   

    If F exceeds G in magnitude, CS will be positive.   

    Arguments   
    =========   

    F       (input) REAL   
            The first component of vector to be rotated.   

    G       (input) REAL   
            The second component of vector to be rotated.   

    CS      (output) REAL   
            The cosine of the rotation.   

    SN      (output) REAL   
            The sine of the rotation.   

    R       (output) REAL   
            The nonzero component of the rotated vector.   

    ===================================================================== */
    /* Initialized data */
    static logical first = TRUE_;
    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    /* Builtin functions */
//    double log(doublereal), pow_ri(real *, integer *), sqrt(doublereal);
    double pow_ri(real *, integer *);
    /* Local variables */
    static integer i__;
    static real scale;
    static integer count;
    static real f1, g1, safmn2, safmx2;
    extern doublereal slamch_(char *);
    static real safmin, eps;



    if (first) {
	first = FALSE_;
	safmin = slamch_("S");
	eps = slamch_("E");
	r__1 = slamch_("B");
	i__1 = (integer) (log(safmin / eps) / log(slamch_("B")) / 
		2.f);
	safmn2 = pow_ri(&r__1, &i__1);
	safmx2 = 1.f / safmn2;
    }
    if (*g == 0.f) {
	*cs = 1.f;
	*sn = 0.f;
	*r__ = *f;
    } else if (*f == 0.f) {
	*cs = 0.f;
	*sn = 1.f;
	*r__ = *g;
    } else {
	f1 = *f;
	g1 = *g;
/* Computing MAX */
	r__1 = dabs(f1), r__2 = dabs(g1);
	scale = df2cmax(r__1,r__2);
	if (scale >= safmx2) {
	    count = 0;
L10:
	    ++count;
	    f1 *= safmn2;
	    g1 *= safmn2;
/* Computing MAX */
	    r__1 = dabs(f1), r__2 = dabs(g1);
	    scale = df2cmax(r__1,r__2);
	    if (scale >= safmx2) {
		goto L10;
	    }
/* Computing 2nd power */
	    r__1 = f1;
/* Computing 2nd power */
	    r__2 = g1;
	    *r__ = sqrt(r__1 * r__1 + r__2 * r__2);
	    *cs = f1 / *r__;
	    *sn = g1 / *r__;
	    i__1 = count;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		*r__ *= safmx2;
/* L20: */
	    }
	} else if (scale <= safmn2) {
	    count = 0;
L30:
	    ++count;
	    f1 *= safmx2;
	    g1 *= safmx2;
/* Computing MAX */
	    r__1 = dabs(f1), r__2 = dabs(g1);
	    scale = df2cmax(r__1,r__2);
	    if (scale <= safmn2) {
		goto L30;
	    }
/* Computing 2nd power */
	    r__1 = f1;
/* Computing 2nd power */
	    r__2 = g1;
	    *r__ = sqrt(r__1 * r__1 + r__2 * r__2);
	    *cs = f1 / *r__;
	    *sn = g1 / *r__;
	    i__1 = count;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		*r__ *= safmn2;
/* L40: */
	    }
	} else {
/* Computing 2nd power */
	    r__1 = f1;
/* Computing 2nd power */
	    r__2 = g1;
	    *r__ = sqrt(r__1 * r__1 + r__2 * r__2);
	    *cs = f1 / *r__;
	    *sn = g1 / *r__;
	}
	if (dabs(*f) > dabs(*g) && *cs < 0.f) {
	    *cs = -(*cs);
	    *sn = -(*sn);
	    *r__ = -(*r__);
	}
    }
    return 0;

/*     End of SLARTG */

} /* slartg_ */




/* Subroutine */ int slascl_(char *type__, integer *kl, integer *ku, real *
	cfrom, real *cto, integer *m, integer *n, real *a, integer *lda, 
	integer *info)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    SLASCL multiplies the M by N real matrix A by the real scalar   
    CTO/CFROM.  This is done without over/underflow as long as the final   
    result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that   
    A may be full, upper triangular, lower triangular, upper Hessenberg,   
    or banded.   

    Arguments   
    =========   

    TYPE    (input) CHARACTER*1   
            TYPE indices the storage type of the input matrix.   
            = 'G':  A is a full matrix.   
            = 'L':  A is a lower triangular matrix.   
            = 'U':  A is an upper triangular matrix.   
            = 'H':  A is an upper Hessenberg matrix.   
            = 'B':  A is a symmetric band matrix with lower bandwidth KL   
                    and upper bandwidth KU and with the only the lower   
                    half stored.   
            = 'Q':  A is a symmetric band matrix with lower bandwidth KL   
                    and upper bandwidth KU and with the only the upper   
                    half stored.   
            = 'Z':  A is a band matrix with lower bandwidth KL and upper   
                    bandwidth KU.   

    KL      (input) INTEGER   
            The lower bandwidth of A.  Referenced only if TYPE = 'B',   
            'Q' or 'Z'.   

    KU      (input) INTEGER   
            The upper bandwidth of A.  Referenced only if TYPE = 'B',   
            'Q' or 'Z'.   

    CFROM   (input) REAL   
    CTO     (input) REAL   
            The matrix A is multiplied by CTO/CFROM. A(I,J) is computed   
            without over/underflow if the final result CTO*A(I,J)/CFROM   
            can be represented without over/underflow.  CFROM must be   
            nonzero.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) REAL array, dimension (LDA,M)   
            The matrix to be multiplied by CTO/CFROM.  See TYPE for the   
            storage type.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= f2cmax(1,M).   

    INFO    (output) INTEGER   
            0  - successful exit   
            <0 - if INFO = -i, the i-th argument had an illegal value.   

    =====================================================================   


       Test the input arguments   

       Parameter adjustments */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    /* Local variables */
    static logical done;
    static real ctoc;
    static integer i__, j;
    extern logical lsame_(char *, char *);
    static integer itype, k1, k2, k3, k4;
    static real cfrom1;
    extern doublereal slamch_(char *);
    static real cfromc;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static real bignum, smlnum, mul, cto1;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]

    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

    /* Function Body */
    *info = 0;

    if (lsame_(type__, "G")) {
	itype = 0;
    } else if (lsame_(type__, "L")) {
	itype = 1;
    } else if (lsame_(type__, "U")) {
	itype = 2;
    } else if (lsame_(type__, "H")) {
	itype = 3;
    } else if (lsame_(type__, "B")) {
	itype = 4;
    } else if (lsame_(type__, "Q")) {
	itype = 5;
    } else if (lsame_(type__, "Z")) {
	itype = 6;
    } else {
	itype = -1;
    }

    if (itype == -1) {
	*info = -1;
    } else if (*cfrom == 0.f) {
	*info = -4;
    } else if (*m < 0) {
	*info = -6;
    } else if (*n < 0 || itype == 4 && *n != *m || itype == 5 && *n != *m) {
	*info = -7;
    } else if (itype <= 3 && *lda < f2cmax(1,*m)) {
	*info = -9;
    } else if (itype >= 4) {
/* Computing MAX */
	i__1 = *m - 1;
	if (*kl < 0 || *kl > f2cmax(i__1,0)) {
	    *info = -2;
	} else /* if(complicated condition) */ {
/* Computing MAX */
	    i__1 = *n - 1;
	    if (*ku < 0 || *ku > f2cmax(i__1,0) || (itype == 4 || itype == 5) && 
		    *kl != *ku) {
		*info = -3;
	    } else if (itype == 4 && *lda < *kl + 1 || itype == 5 && *lda < *
		    ku + 1 || itype == 6 && *lda < (*kl << 1) + *ku + 1) {
		*info = -9;
	    }
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SLASCL", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *m == 0) {
	return 0;
    }

/*     Get machine parameters */

    smlnum = slamch_("S");
    bignum = 1.f / smlnum;

    cfromc = *cfrom;
    ctoc = *cto;

L10:
    cfrom1 = cfromc * smlnum;
    cto1 = ctoc / bignum;
    if (dabs(cfrom1) > dabs(ctoc) && ctoc != 0.f) {
	mul = smlnum;
	done = FALSE_;
	cfromc = cfrom1;
    } else if (dabs(cto1) > dabs(cfromc)) {
	mul = bignum;
	done = FALSE_;
	ctoc = cto1;
    } else {
	mul = ctoc / cfromc;
	done = TRUE_;
    }

    if (itype == 0) {

/*        Full matrix */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a_ref(i__, j) = a_ref(i__, j) * mul;
/* L20: */
	    }
/* L30: */
	}

    } else if (itype == 1) {

/*        Lower triangular matrix */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = j; i__ <= i__2; ++i__) {
		a_ref(i__, j) = a_ref(i__, j) * mul;
/* L40: */
	    }
/* L50: */
	}

    } else if (itype == 2) {

/*        Upper triangular matrix */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = f2cmin(j,*m);
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a_ref(i__, j) = a_ref(i__, j) * mul;
/* L60: */
	    }
/* L70: */
	}

    } else if (itype == 3) {

/*        Upper Hessenberg matrix */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	    i__3 = j + 1;
	    i__2 = f2cmin(i__3,*m);
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a_ref(i__, j) = a_ref(i__, j) * mul;
/* L80: */
	    }
/* L90: */
	}

    } else if (itype == 4) {

/*        Lower half of a symmetric band matrix */

	k3 = *kl + 1;
	k4 = *n + 1;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	    i__3 = k3, i__4 = k4 - j;
	    i__2 = f2cmin(i__3,i__4);
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a_ref(i__, j) = a_ref(i__, j) * mul;
/* L100: */
	    }
/* L110: */
	}

    } else if (itype == 5) {

/*        Upper half of a symmetric band matrix */

	k1 = *ku + 2;
	k3 = *ku + 1;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	    i__2 = k1 - j;
	    i__3 = k3;
	    for (i__ = f2cmax(i__2,1); i__ <= i__3; ++i__) {
		a_ref(i__, j) = a_ref(i__, j) * mul;
/* L120: */
	    }
/* L130: */
	}

    } else if (itype == 6) {

/*        Band matrix */

	k1 = *kl + *ku + 2;
	k2 = *kl + 1;
	k3 = (*kl << 1) + *ku + 1;
	k4 = *kl + *ku + 1 + *m;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	    i__3 = k1 - j;
/* Computing MIN */
	    i__4 = k3, i__5 = k4 - j;
	    i__2 = f2cmin(i__4,i__5);
	    for (i__ = f2cmax(i__3,k2); i__ <= i__2; ++i__) {
		a_ref(i__, j) = a_ref(i__, j) * mul;
/* L140: */
	    }
/* L150: */
	}

    }

    if (! done) {
	goto L10;
    }

    return 0;

/*     End of SLASCL */

} /* slascl_ */

#undef a_ref





/* Subroutine */ int slaset_(char *uplo, integer *m, integer *n, real *alpha, 
	real *beta, real *a, integer *lda)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLASET initializes an m-by-n matrix A to BETA on the diagonal and   
    ALPHA on the offdiagonals.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies the part of the matrix A to be set.   
            = 'U':      Upper triangular part is set; the strictly lower   
                        triangular part of A is not changed.   
            = 'L':      Lower triangular part is set; the strictly upper   
                        triangular part of A is not changed.   
            Otherwise:  All of the matrix A is set.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    ALPHA   (input) REAL   
            The constant to which the offdiagonal elements are to be set.   

    BETA    (input) REAL   
            The constant to which the diagonal elements are to be set.   

    A       (input/output) REAL array, dimension (LDA,N)   
            On exit, the leading m-by-n submatrix of A is set as follows:   

            if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,   
            if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,   
            otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,   

            and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= f2cmax(1,M).   

   =====================================================================   


       Parameter adjustments */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]

    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

    /* Function Body */
    if (lsame_(uplo, "U")) {

/*        Set the strictly upper triangular or trapezoidal part of the   
          array to ALPHA. */

	i__1 = *n;
	for (j = 2; j <= i__1; ++j) {
/* Computing MIN */
	    i__3 = j - 1;
	    i__2 = f2cmin(i__3,*m);
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a_ref(i__, j) = *alpha;
/* L10: */
	    }
/* L20: */
	}

    } else if (lsame_(uplo, "L")) {

/*        Set the strictly lower triangular or trapezoidal part of the   
          array to ALPHA. */

	i__1 = f2cmin(*m,*n);
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
		a_ref(i__, j) = *alpha;
/* L30: */
	    }
/* L40: */
	}

    } else {

/*        Set the leading m-by-n submatrix to ALPHA. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a_ref(i__, j) = *alpha;
/* L50: */
	    }
/* L60: */
	}
    }

/*     Set the first f2cmin(M,N) diagonal elements to BETA. */

    i__1 = f2cmin(*m,*n);
    for (i__ = 1; i__ <= i__1; ++i__) {
	a_ref(i__, i__) = *beta;
/* L70: */
    }

    return 0;

/*     End of SLASET */

} /* slaset_ */

#undef a_ref





/* Subroutine */ int slasr_(char *side, char *pivot, char *direct, integer *m,
	 integer *n, real *c__, real *s, real *a, integer *lda)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLASR   performs the transformation   

       A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )   

       A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )   

    where A is an m by n real matrix and P is an orthogonal matrix,   
    consisting of a sequence of plane rotations determined by the   
    parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l'   
    and z = n when SIDE = 'R' or 'r' ):   

    When  DIRECT = 'F' or 'f'  ( Forward sequence ) then   

       P = P( z - 1 )*...*P( 2 )*P( 1 ),   

    and when DIRECT = 'B' or 'b'  ( Backward sequence ) then   

       P = P( 1 )*P( 2 )*...*P( z - 1 ),   

    where  P( k ) is a plane rotation matrix for the following planes:   

       when  PIVOT = 'V' or 'v'  ( Variable pivot ),   
          the plane ( k, k + 1 )   

       when  PIVOT = 'T' or 't'  ( Top pivot ),   
          the plane ( 1, k + 1 )   

       when  PIVOT = 'B' or 'b'  ( Bottom pivot ),   
          the plane ( k, z )   

    c( k ) and s( k )  must contain the  cosine and sine that define the   
    matrix  P( k ).  The two by two plane rotation part of the matrix   
    P( k ), R( k ), is assumed to be of the form   

       R( k ) = (  c( k )  s( k ) ).   
                ( -s( k )  c( k ) )   

    This version vectorises across rows of the array A when SIDE = 'L'.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            Specifies whether the plane rotation matrix P is applied to   
            A on the left or the right.   
            = 'L':  Left, compute A := P*A   
            = 'R':  Right, compute A:= A*P'   

    DIRECT  (input) CHARACTER*1   
            Specifies whether P is a forward or backward sequence of   
            plane rotations.   
            = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )   
            = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )   

    PIVOT   (input) CHARACTER*1   
            Specifies the plane for which P(k) is a plane rotation   
            matrix.   
            = 'V':  Variable pivot, the plane (k,k+1)   
            = 'T':  Top pivot, the plane (1,k+1)   
            = 'B':  Bottom pivot, the plane (k,z)   

    M       (input) INTEGER   
            The number of rows of the matrix A.  If m <= 1, an immediate   
            return is effected.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  If n <= 1, an   
            immediate return is effected.   

    C, S    (input) REAL arrays, dimension   
                    (M-1) if SIDE = 'L'   
                    (N-1) if SIDE = 'R'   
            c(k) and s(k) contain the cosine and sine that define the   
            matrix P(k).  The two by two plane rotation part of the   
            matrix P(k), R(k), is assumed to be of the form   
            R( k ) = (  c( k )  s( k ) ).   
                     ( -s( k )  c( k ) )   

    A       (input/output) REAL array, dimension (LDA,N)   
            The m by n matrix A.  On exit, A is overwritten by P*A if   
            SIDE = 'R' or by A*P' if SIDE = 'L'.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= f2cmax(1,M).   

    =====================================================================   


       Test the input parameters   

       Parameter adjustments */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    static integer info;
    static real temp;
    static integer i__, j;
    extern logical lsame_(char *, char *);
    static real ctemp, stemp;
    extern /* Subroutine */ int xerbla_(char *, integer *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]

    --c__;
    --s;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

    /* Function Body */
    info = 0;
    if (! (lsame_(side, "L") || lsame_(side, "R"))) {
	info = 1;
    } else if (! (lsame_(pivot, "V") || lsame_(pivot, 
	    "T") || lsame_(pivot, "B"))) {
	info = 2;
    } else if (! (lsame_(direct, "F") || lsame_(direct, 
	    "B"))) {
	info = 3;
    } else if (*m < 0) {
	info = 4;
    } else if (*n < 0) {
	info = 5;
    } else if (*lda < f2cmax(1,*m)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("SLASR ", &info);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }
    if (lsame_(side, "L")) {

/*        Form  P * A */

	if (lsame_(pivot, "V")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1.f || stemp != 0.f) {
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a_ref(j + 1, i__);
			    a_ref(j + 1, i__) = ctemp * temp - stemp * a_ref(
				    j, i__);
			    a_ref(j, i__) = stemp * temp + ctemp * a_ref(j, 
				    i__);
/* L10: */
			}
		    }
/* L20: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1.f || stemp != 0.f) {
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a_ref(j + 1, i__);
			    a_ref(j + 1, i__) = ctemp * temp - stemp * a_ref(
				    j, i__);
			    a_ref(j, i__) = stemp * temp + ctemp * a_ref(j, 
				    i__);
/* L30: */
			}
		    }
/* L40: */
		}
	    }
	} else if (lsame_(pivot, "T")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m;
		for (j = 2; j <= i__1; ++j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != 1.f || stemp != 0.f) {
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a_ref(j, i__);
			    a_ref(j, i__) = ctemp * temp - stemp * a_ref(1, 
				    i__);
			    a_ref(1, i__) = stemp * temp + ctemp * a_ref(1, 
				    i__);
/* L50: */
			}
		    }
/* L60: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m; j >= 2; --j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != 1.f || stemp != 0.f) {
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a_ref(j, i__);
			    a_ref(j, i__) = ctemp * temp - stemp * a_ref(1, 
				    i__);
			    a_ref(1, i__) = stemp * temp + ctemp * a_ref(1, 
				    i__);
/* L70: */
			}
		    }
/* L80: */
		}
	    }
	} else if (lsame_(pivot, "B")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1.f || stemp != 0.f) {
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a_ref(j, i__);
			    a_ref(j, i__) = stemp * a_ref(*m, i__) + ctemp * 
				    temp;
			    a_ref(*m, i__) = ctemp * a_ref(*m, i__) - stemp * 
				    temp;
/* L90: */
			}
		    }
/* L100: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1.f || stemp != 0.f) {
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a_ref(j, i__);
			    a_ref(j, i__) = stemp * a_ref(*m, i__) + ctemp * 
				    temp;
			    a_ref(*m, i__) = ctemp * a_ref(*m, i__) - stemp * 
				    temp;
/* L110: */
			}
		    }
/* L120: */
		}
	    }
	}
    } else if (lsame_(side, "R")) {

/*        Form A * P' */

	if (lsame_(pivot, "V")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1.f || stemp != 0.f) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a_ref(i__, j + 1);
			    a_ref(i__, j + 1) = ctemp * temp - stemp * a_ref(
				    i__, j);
			    a_ref(i__, j) = stemp * temp + ctemp * a_ref(i__, 
				    j);
/* L130: */
			}
		    }
/* L140: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *n - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1.f || stemp != 0.f) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a_ref(i__, j + 1);
			    a_ref(i__, j + 1) = ctemp * temp - stemp * a_ref(
				    i__, j);
			    a_ref(i__, j) = stemp * temp + ctemp * a_ref(i__, 
				    j);
/* L150: */
			}
		    }
/* L160: */
		}
	    }
	} else if (lsame_(pivot, "T")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n;
		for (j = 2; j <= i__1; ++j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != 1.f || stemp != 0.f) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a_ref(i__, j);
			    a_ref(i__, j) = ctemp * temp - stemp * a_ref(i__, 
				    1);
			    a_ref(i__, 1) = stemp * temp + ctemp * a_ref(i__, 
				    1);
/* L170: */
			}
		    }
/* L180: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *n; j >= 2; --j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != 1.f || stemp != 0.f) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a_ref(i__, j);
			    a_ref(i__, j) = ctemp * temp - stemp * a_ref(i__, 
				    1);
			    a_ref(i__, 1) = stemp * temp + ctemp * a_ref(i__, 
				    1);
/* L190: */
			}
		    }
/* L200: */
		}
	    }
	} else if (lsame_(pivot, "B")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1.f || stemp != 0.f) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a_ref(i__, j);
			    a_ref(i__, j) = stemp * a_ref(i__, *n) + ctemp * 
				    temp;
			    a_ref(i__, *n) = ctemp * a_ref(i__, *n) - stemp * 
				    temp;
/* L210: */
			}
		    }
/* L220: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *n - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1.f || stemp != 0.f) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a_ref(i__, j);
			    a_ref(i__, j) = stemp * a_ref(i__, *n) + ctemp * 
				    temp;
			    a_ref(i__, *n) = ctemp * a_ref(i__, *n) - stemp * 
				    temp;
/* L230: */
			}
		    }
/* L240: */
		}
	    }
	}
    }

    return 0;

/*     End of SLASR */

} /* slasr_ */

#undef a_ref





/* Subroutine */ int slasrt_(char *id, integer *n, real *d__, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    Sort the numbers in D in increasing order (if ID = 'I') or   
    in decreasing order (if ID = 'D' ).   

    Use Quick Sort, reverting to Insertion sort on arrays of   
    size <= 20. Dimension of STACK limits N to about 2**32.   

    Arguments   
    =========   

    ID      (input) CHARACTER*1   
            = 'I': sort D in increasing order;   
            = 'D': sort D in decreasing order.   

    N       (input) INTEGER   
            The length of the array D.   

    D       (input/output) REAL array, dimension (N)   
            On entry, the array to be sorted.   
            On exit, D has been sorted into increasing order   
            (D(1) <= ... <= D(N) ) or into decreasing order   
            (D(1) >= ... >= D(N) ), depending on ID.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    =====================================================================   


       Test the input paramters.   

       Parameter adjustments */
    /* System generated locals */
    integer i__1, i__2;
    /* Local variables */
    static integer endd, i__, j;
    extern logical lsame_(char *, char *);
    static integer stack[64]	/* was [2][32] */;
    static real dmnmx, d1, d2, d3;
    static integer start;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static integer stkpnt, dir;
    static real tmp;
#define stack_ref(a_1,a_2) stack[(a_2)*2 + a_1 - 3]

    --d__;

    /* Function Body */
    *info = 0;
    dir = -1;
    if (lsame_(id, "D")) {
	dir = 0;
    } else if (lsame_(id, "I")) {
	dir = 1;
    }
    if (dir == -1) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SLASRT", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 1) {
	return 0;
    }

    stkpnt = 1;
    stack_ref(1, 1) = 1;
    stack_ref(2, 1) = *n;
L10:
    start = stack_ref(1, stkpnt);
    endd = stack_ref(2, stkpnt);
    --stkpnt;
    if (endd - start <= 20 && endd - start > 0) {

/*        Do Insertion sort on D( START:ENDD ) */

	if (dir == 0) {

/*           Sort into decreasing order */

	    i__1 = endd;
	    for (i__ = start + 1; i__ <= i__1; ++i__) {
		i__2 = start + 1;
		for (j = i__; j >= i__2; --j) {
		    if (d__[j] > d__[j - 1]) {
			dmnmx = d__[j];
			d__[j] = d__[j - 1];
			d__[j - 1] = dmnmx;
		    } else {
			goto L30;
		    }
/* L20: */
		}
L30:
		;
	    }

	} else {

/*           Sort into increasing order */

	    i__1 = endd;
	    for (i__ = start + 1; i__ <= i__1; ++i__) {
		i__2 = start + 1;
		for (j = i__; j >= i__2; --j) {
		    if (d__[j] < d__[j - 1]) {
			dmnmx = d__[j];
			d__[j] = d__[j - 1];
			d__[j - 1] = dmnmx;
		    } else {
			goto L50;
		    }
/* L40: */
		}
L50:
		;
	    }

	}

    } else if (endd - start > 20) {

/*        Partition D( START:ENDD ) and stack parts, largest one first   

          Choose partition entry as median of 3 */

	d1 = d__[start];
	d2 = d__[endd];
	i__ = (start + endd) / 2;
	d3 = d__[i__];
	if (d1 < d2) {
	    if (d3 < d1) {
		dmnmx = d1;
	    } else if (d3 < d2) {
		dmnmx = d3;
	    } else {
		dmnmx = d2;
	    }
	} else {
	    if (d3 < d2) {
		dmnmx = d2;
	    } else if (d3 < d1) {
		dmnmx = d3;
	    } else {
		dmnmx = d1;
	    }
	}

	if (dir == 0) {

/*           Sort into decreasing order */

	    i__ = start - 1;
	    j = endd + 1;
L60:
L70:
	    --j;
	    if (d__[j] < dmnmx) {
		goto L70;
	    }
L80:
	    ++i__;
	    if (d__[i__] > dmnmx) {
		goto L80;
	    }
	    if (i__ < j) {
		tmp = d__[i__];
		d__[i__] = d__[j];
		d__[j] = tmp;
		goto L60;
	    }
	    if (j - start > endd - j - 1) {
		++stkpnt;
		stack_ref(1, stkpnt) = start;
		stack_ref(2, stkpnt) = j;
		++stkpnt;
		stack_ref(1, stkpnt) = j + 1;
		stack_ref(2, stkpnt) = endd;
	    } else {
		++stkpnt;
		stack_ref(1, stkpnt) = j + 1;
		stack_ref(2, stkpnt) = endd;
		++stkpnt;
		stack_ref(1, stkpnt) = start;
		stack_ref(2, stkpnt) = j;
	    }
	} else {

/*           Sort into increasing order */

	    i__ = start - 1;
	    j = endd + 1;
L90:
L100:
	    --j;
	    if (d__[j] > dmnmx) {
		goto L100;
	    }
L110:
	    ++i__;
	    if (d__[i__] < dmnmx) {
		goto L110;
	    }
	    if (i__ < j) {
		tmp = d__[i__];
		d__[i__] = d__[j];
		d__[j] = tmp;
		goto L90;
	    }
	    if (j - start > endd - j - 1) {
		++stkpnt;
		stack_ref(1, stkpnt) = start;
		stack_ref(2, stkpnt) = j;
		++stkpnt;
		stack_ref(1, stkpnt) = j + 1;
		stack_ref(2, stkpnt) = endd;
	    } else {
		++stkpnt;
		stack_ref(1, stkpnt) = j + 1;
		stack_ref(2, stkpnt) = endd;
		++stkpnt;
		stack_ref(1, stkpnt) = start;
		stack_ref(2, stkpnt) = j;
	    }
	}
    }
    if (stkpnt > 0) {
	goto L10;
    }
    return 0;

/*     End of SLASRT */

} /* slasrt_ */

#undef stack_ref





/* Subroutine */ int slassq_(integer *n, real *x, integer *incx, real *scale, 
	real *sumsq)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    SLASSQ  returns the values  scl  and  smsq  such that   

       ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,   

    where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is   
    assumed to be non-negative and  scl  returns the value   

       scl = f2cmax( scale, abs( x( i ) ) ).   

    scale and sumsq must be supplied in SCALE and SUMSQ and   
    scl and smsq are overwritten on SCALE and SUMSQ respectively.   

    The routine makes only one pass through the vector x.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of elements to be used from the vector X.   

    X       (input) REAL array, dimension (N)   
            The vector for which a scaled sum of squares is computed.   
               x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.   

    INCX    (input) INTEGER   
            The increment between successive values of the vector X.   
            INCX > 0.   

    SCALE   (input/output) REAL   
            On entry, the value  scale  in the equation above.   
            On exit, SCALE is overwritten with  scl , the scaling factor   
            for the sum of squares.   

    SUMSQ   (input/output) REAL   
            On entry, the value  sumsq  in the equation above.   
            On exit, SUMSQ is overwritten with  smsq , the basic sum of   
            squares from which  scl  has been factored out.   

   =====================================================================   


       Parameter adjustments */
    /* System generated locals */
    integer i__1, i__2;
    real r__1;
    /* Local variables */
    static real absxi;
    static integer ix;

    --x;

    /* Function Body */
    if (*n > 0) {
	i__1 = (*n - 1) * *incx + 1;
	i__2 = *incx;
	for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
	    if (x[ix] != 0.f) {
		absxi = (r__1 = x[ix], dabs(r__1));
		if (*scale < absxi) {
/* Computing 2nd power */
		    r__1 = *scale / absxi;
		    *sumsq = *sumsq * (r__1 * r__1) + 1;
		    *scale = absxi;
		} else {
/* Computing 2nd power */
		    r__1 = absxi / *scale;
		    *sumsq += r__1 * r__1;
		}
	    }
/* L10: */
	}
    }
    return 0;

/*     End of SLASSQ */

} /* slassq_ */




/* Subroutine */ int slatrd_(char *uplo, integer *n, integer *nb, real *a, 
	integer *lda, real *e, real *tau, real *w, integer *ldw)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLATRD reduces NB rows and columns of a real symmetric matrix A to   
    symmetric tridiagonal form by an orthogonal similarity   
    transformation Q' * A * Q, and returns the matrices V and W which are   
    needed to apply the transformation to the unreduced part of A.   

    If UPLO = 'U', SLATRD reduces the last NB rows and columns of a   
    matrix, of which the upper triangle is supplied;   
    if UPLO = 'L', SLATRD reduces the first NB rows and columns of a   
    matrix, of which the lower triangle is supplied.   

    This is an auxiliary routine called by SSYTRD.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is stored:   
            = 'U': Upper triangular   
            = 'L': Lower triangular   

    N       (input) INTEGER   
            The order of the matrix A.   

    NB      (input) INTEGER   
            The number of rows and columns to be reduced.   

    A       (input/output) REAL array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading   
            n-by-n upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading n-by-n lower triangular part of A contains the lower   
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   
            On exit:   
            if UPLO = 'U', the last NB columns have been reduced to   
              tridiagonal form, with the diagonal elements overwriting   
              the diagonal elements of A; the elements above the diagonal   
              with the array TAU, represent the orthogonal matrix Q as a   
              product of elementary reflectors;   
            if UPLO = 'L', the first NB columns have been reduced to   
              tridiagonal form, with the diagonal elements overwriting   
              the diagonal elements of A; the elements below the diagonal   
              with the array TAU, represent the  orthogonal matrix Q as a   
              product of elementary reflectors.   
            See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= (1,N).   

    E       (output) REAL array, dimension (N-1)   
            If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal   
            elements of the last NB columns of the reduced matrix;   
            if UPLO = 'L', E(1:nb) contains the subdiagonal elements of   
            the first NB columns of the reduced matrix.   

    TAU     (output) REAL array, dimension (N-1)   
            The scalar factors of the elementary reflectors, stored in   
            TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.   
            See Further Details.   

    W       (output) REAL array, dimension (LDW,NB)   
            The n-by-nb matrix W required to update the unreduced part   
            of A.   

    LDW     (input) INTEGER   
            The leading dimension of the array W. LDW >= f2cmax(1,N).   

    Further Details   
    ===============   

    If UPLO = 'U', the matrix Q is represented as a product of elementary   
    reflectors   

       Q = H(n) H(n-1) . . . H(n-nb+1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),   
    and tau in TAU(i-1).   

    If UPLO = 'L', the matrix Q is represented as a product of elementary   
    reflectors   

       Q = H(1) H(2) . . . H(nb).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),   
    and tau in TAU(i).   

    The elements of the vectors v together form the n-by-nb matrix V   
    which is needed, with W, to apply the transformation to the unreduced   
    part of the matrix, using a symmetric rank-2k update of the form:   
    A := A - V*W' - W*V'.   

    The contents of A on exit are illustrated by the following examples   
    with n = 5 and nb = 2:   

    if UPLO = 'U':                       if UPLO = 'L':   

      (  a   a   a   v4  v5 )              (  d                  )   
      (      a   a   v4  v5 )              (  1   d              )   
      (          a   1   v5 )              (  v1  1   a          )   
      (              d   1  )              (  v1  v2  a   a      )   
      (                  d  )              (  v1  v2  a   a   a  )   

    where d denotes a diagonal element of the reduced matrix, a denotes   
    an element of the original matrix that is unchanged, and vi denotes   
    an element of the vector defining H(i).   

    =====================================================================   


       Quick return if possible   

       Parameter adjustments */
    /* Table of constant values */
    static real c_b5 = -1.f;
    static real c_b6 = 1.f;
    static integer c__1 = 1;
    static real c_b16 = 0.f;
    
    /* System generated locals */
    integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3;
    /* Local variables */
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static integer i__;
    static real alpha;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *), 
	    sgemv_(char *, integer *, integer *, real *, real *, integer *, 
	    real *, integer *, real *, real *, integer *), saxpy_(
	    integer *, real *, real *, integer *, real *, integer *), ssymv_(
	    char *, integer *, real *, real *, integer *, real *, integer *, 
	    real *, real *, integer *);
    static integer iw;
    extern /* Subroutine */ int slarfg_(integer *, real *, real *, integer *, 
	    real *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define w_ref(a_1,a_2) w[(a_2)*w_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --e;
    --tau;
    w_dim1 = *ldw;
    w_offset = 1 + w_dim1 * 1;
    w -= w_offset;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }

    if (lsame_(uplo, "U")) {

/*        Reduce last NB columns of upper triangle */

	i__1 = *n - *nb + 1;
	for (i__ = *n; i__ >= i__1; --i__) {
	    iw = i__ - *n + *nb;
	    if (i__ < *n) {

/*              Update A(1:i,i) */

		i__2 = *n - i__;
		sgemv_("No transpose", &i__, &i__2, &c_b5, &a_ref(1, i__ + 1),
			 lda, &w_ref(i__, iw + 1), ldw, &c_b6, &a_ref(1, i__),
			 &c__1);
		i__2 = *n - i__;
		sgemv_("No transpose", &i__, &i__2, &c_b5, &w_ref(1, iw + 1), 
			ldw, &a_ref(i__, i__ + 1), lda, &c_b6, &a_ref(1, i__),
			 &c__1);
	    }
	    if (i__ > 1) {

/*              Generate elementary reflector H(i) to annihilate   
                A(1:i-2,i) */

		i__2 = i__ - 1;
		slarfg_(&i__2, &a_ref(i__ - 1, i__), &a_ref(1, i__), &c__1, &
			tau[i__ - 1]);
		e[i__ - 1] = a_ref(i__ - 1, i__);
		a_ref(i__ - 1, i__) = 1.f;

/*              Compute W(1:i-1,i) */

		i__2 = i__ - 1;
		ssymv_("Upper", &i__2, &c_b6, &a[a_offset], lda, &a_ref(1, 
			i__), &c__1, &c_b16, &w_ref(1, iw), &c__1);
		if (i__ < *n) {
		    i__2 = i__ - 1;
		    i__3 = *n - i__;
		    sgemv_("Transpose", &i__2, &i__3, &c_b6, &w_ref(1, iw + 1)
			    , ldw, &a_ref(1, i__), &c__1, &c_b16, &w_ref(i__ 
			    + 1, iw), &c__1);
		    i__2 = i__ - 1;
		    i__3 = *n - i__;
		    sgemv_("No transpose", &i__2, &i__3, &c_b5, &a_ref(1, i__ 
			    + 1), lda, &w_ref(i__ + 1, iw), &c__1, &c_b6, &
			    w_ref(1, iw), &c__1);
		    i__2 = i__ - 1;
		    i__3 = *n - i__;
		    sgemv_("Transpose", &i__2, &i__3, &c_b6, &a_ref(1, i__ + 
			    1), lda, &a_ref(1, i__), &c__1, &c_b16, &w_ref(
			    i__ + 1, iw), &c__1);
		    i__2 = i__ - 1;
		    i__3 = *n - i__;
		    sgemv_("No transpose", &i__2, &i__3, &c_b5, &w_ref(1, iw 
			    + 1), ldw, &w_ref(i__ + 1, iw), &c__1, &c_b6, &
			    w_ref(1, iw), &c__1);
		}
		i__2 = i__ - 1;
		sscal_(&i__2, &tau[i__ - 1], &w_ref(1, iw), &c__1);
		i__2 = i__ - 1;
		alpha = tau[i__ - 1] * -.5f * sdot_(&i__2, &w_ref(1, iw), &
			c__1, &a_ref(1, i__), &c__1);
		i__2 = i__ - 1;
		saxpy_(&i__2, &alpha, &a_ref(1, i__), &c__1, &w_ref(1, iw), &
			c__1);
	    }

/* L10: */
	}
    } else {

/*        Reduce first NB columns of lower triangle */

	i__1 = *nb;
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Update A(i:n,i) */

	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    sgemv_("No transpose", &i__2, &i__3, &c_b5, &a_ref(i__, 1), lda, &
		    w_ref(i__, 1), ldw, &c_b6, &a_ref(i__, i__), &c__1);
	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    sgemv_("No transpose", &i__2, &i__3, &c_b5, &w_ref(i__, 1), ldw, &
		    a_ref(i__, 1), lda, &c_b6, &a_ref(i__, i__), &c__1);
	    if (i__ < *n) {

/*              Generate elementary reflector H(i) to annihilate   
                A(i+2:n,i)   

   Computing MIN */
		i__2 = i__ + 2;
		i__3 = *n - i__;
		slarfg_(&i__3, &a_ref(i__ + 1, i__), &a_ref(f2cmin(i__2,*n), i__)
			, &c__1, &tau[i__]);
		e[i__] = a_ref(i__ + 1, i__);
		a_ref(i__ + 1, i__) = 1.f;

/*              Compute W(i+1:n,i) */

		i__2 = *n - i__;
		ssymv_("Lower", &i__2, &c_b6, &a_ref(i__ + 1, i__ + 1), lda, &
			a_ref(i__ + 1, i__), &c__1, &c_b16, &w_ref(i__ + 1, 
			i__), &c__1);
		i__2 = *n - i__;
		i__3 = i__ - 1;
		sgemv_("Transpose", &i__2, &i__3, &c_b6, &w_ref(i__ + 1, 1), 
			ldw, &a_ref(i__ + 1, i__), &c__1, &c_b16, &w_ref(1, 
			i__), &c__1);
		i__2 = *n - i__;
		i__3 = i__ - 1;
		sgemv_("No transpose", &i__2, &i__3, &c_b5, &a_ref(i__ + 1, 1)
			, lda, &w_ref(1, i__), &c__1, &c_b6, &w_ref(i__ + 1, 
			i__), &c__1);
		i__2 = *n - i__;
		i__3 = i__ - 1;
		sgemv_("Transpose", &i__2, &i__3, &c_b6, &a_ref(i__ + 1, 1), 
			lda, &a_ref(i__ + 1, i__), &c__1, &c_b16, &w_ref(1, 
			i__), &c__1);
		i__2 = *n - i__;
		i__3 = i__ - 1;
		sgemv_("No transpose", &i__2, &i__3, &c_b5, &w_ref(i__ + 1, 1)
			, ldw, &w_ref(1, i__), &c__1, &c_b6, &w_ref(i__ + 1, 
			i__), &c__1);
		i__2 = *n - i__;
		sscal_(&i__2, &tau[i__], &w_ref(i__ + 1, i__), &c__1);
		i__2 = *n - i__;
		alpha = tau[i__] * -.5f * sdot_(&i__2, &w_ref(i__ + 1, i__), &
			c__1, &a_ref(i__ + 1, i__), &c__1);
		i__2 = *n - i__;
		saxpy_(&i__2, &alpha, &a_ref(i__ + 1, i__), &c__1, &w_ref(i__ 
			+ 1, i__), &c__1);
	    }

/* L20: */
	}
    }

    return 0;

/*     End of SLATRD */

} /* slatrd_ */

#undef w_ref
#undef a_ref





doublereal snrm2_(integer *n, real *x, integer *incx)
{
/*        The following loop is equivalent to this call to the LAPACK   
          auxiliary routine:   
          CALL SLASSQ( N, X, INCX, SCALE, SSQ ) */
    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1;
    /* Builtin functions */
//    double sqrt(doublereal);
    /* Local variables */
    static real norm, scale, absxi;
    static integer ix;
    static real ssq;
/*  SNRM2 returns the euclidean norm of a vector via the function   
    name, so that   
       SNRM2 := sqrt( x'*x )   
    -- This version written on 25-October-1982.   
       Modified on 14-October-1993 to inline the call to SLASSQ.   
       Sven Hammarling, Nag Ltd.   
       Parameter adjustments */
    --x;
    /* Function Body */
    if (*n < 1 || *incx < 1) {
	norm = 0.f;
    } else if (*n == 1) {
	norm = dabs(x[1]);
    } else {
	scale = 0.f;
	ssq = 1.f;


	i__1 = (*n - 1) * *incx + 1;
	i__2 = *incx;
	for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
	    if (x[ix] != 0.f) {
		absxi = (r__1 = x[ix], dabs(r__1));
		if (scale < absxi) {
/* Computing 2nd power */
		    r__1 = scale / absxi;
		    ssq = ssq * (r__1 * r__1) + 1.f;
		    scale = absxi;
		} else {
/* Computing 2nd power */
		    r__1 = absxi / scale;
		    ssq += r__1 * r__1;
		}
	    }
/* L10: */
	}
	norm = scale * sqrt(ssq);
    }

    ret_val = norm;
    return ret_val;

/*     End of SNRM2. */

} /* snrm2_ */




/* Subroutine */ int sorg2l_(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    SORG2L generates an m by n real matrix Q with orthonormal columns,   
    which is defined as the last n columns of a product of k elementary   
    reflectors of order m   

          Q  =  H(k) . . . H(2) H(1)   

    as returned by SGEQLF.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix Q. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q. M >= N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines the   
            matrix Q. N >= K >= 0.   

    A       (input/output) REAL array, dimension (LDA,N)   
            On entry, the (n-k+i)-th column must contain the vector which   
            defines the elementary reflector H(i), for i = 1,2,...,k, as   
            returned by SGEQLF in the last k columns of its array   
            argument A.   
            On exit, the m by n matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= f2cmax(1,M).   

    TAU     (input) REAL array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by SGEQLF.   

    WORK    (workspace) REAL array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument has an illegal value   

    =====================================================================   


       Test the input arguments   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1;
    /* Local variables */
    static integer i__, j, l;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *), 
	    slarf_(char *, integer *, integer *, real *, integer *, real *, 
	    real *, integer *, real *);
    static integer ii;
    extern /* Subroutine */ int xerbla_(char *, integer *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0 || *n > *m) {
	*info = -2;
    } else if (*k < 0 || *k > *n) {
	*info = -3;
    } else if (*lda < f2cmax(1,*m)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SORG2L", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	return 0;
    }

/*     Initialise columns 1:n-k to columns of the unit matrix */

    i__1 = *n - *k;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (l = 1; l <= i__2; ++l) {
	    a_ref(l, j) = 0.f;
/* L10: */
	}
	a_ref(*m - *n + j, j) = 1.f;
/* L20: */
    }

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = *n - *k + i__;

/*        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left */

	a_ref(*m - *n + ii, ii) = 1.f;
	i__2 = *m - *n + ii;
	i__3 = ii - 1;
	slarf_("Left", &i__2, &i__3, &a_ref(1, ii), &c__1, &tau[i__], &a[
		a_offset], lda, &work[1]);
	i__2 = *m - *n + ii - 1;
	r__1 = -tau[i__];
	sscal_(&i__2, &r__1, &a_ref(1, ii), &c__1);
	a_ref(*m - *n + ii, ii) = 1.f - tau[i__];

/*        Set A(m-k+i+1:m,n-k+i) to zero */

	i__2 = *m;
	for (l = *m - *n + ii + 1; l <= i__2; ++l) {
	    a_ref(l, ii) = 0.f;
/* L30: */
	}
/* L40: */
    }
    return 0;

/*     End of SORG2L */

} /* sorg2l_ */

#undef a_ref





/* Subroutine */ int sorg2r_(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    SORG2R generates an m by n real matrix Q with orthonormal columns,   
    which is defined as the first n columns of a product of k elementary   
    reflectors of order m   

          Q  =  H(1) H(2) . . . H(k)   

    as returned by SGEQRF.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix Q. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q. M >= N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines the   
            matrix Q. N >= K >= 0.   

    A       (input/output) REAL array, dimension (LDA,N)   
            On entry, the i-th column must contain the vector which   
            defines the elementary reflector H(i), for i = 1,2,...,k, as   
            returned by SGEQRF in the first k columns of its array   
            argument A.   
            On exit, the m-by-n matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= f2cmax(1,M).   

    TAU     (input) REAL array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by SGEQRF.   

    WORK    (workspace) REAL array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument has an illegal value   

    =====================================================================   


       Test the input arguments   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    real r__1;
    /* Local variables */
    static integer i__, j, l;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *), 
	    slarf_(char *, integer *, integer *, real *, integer *, real *, 
	    real *, integer *, real *), xerbla_(char *, integer *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0 || *n > *m) {
	*info = -2;
    } else if (*k < 0 || *k > *n) {
	*info = -3;
    } else if (*lda < f2cmax(1,*m)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SORG2R", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	return 0;
    }

/*     Initialise columns k+1:n to columns of the unit matrix */

    i__1 = *n;
    for (j = *k + 1; j <= i__1; ++j) {
	i__2 = *m;
	for (l = 1; l <= i__2; ++l) {
	    a_ref(l, j) = 0.f;
/* L10: */
	}
	a_ref(j, j) = 1.f;
/* L20: */
    }

    for (i__ = *k; i__ >= 1; --i__) {

/*        Apply H(i) to A(i:m,i:n) from the left */

	if (i__ < *n) {
	    a_ref(i__, i__) = 1.f;
	    i__1 = *m - i__ + 1;
	    i__2 = *n - i__;
	    slarf_("Left", &i__1, &i__2, &a_ref(i__, i__), &c__1, &tau[i__], &
		    a_ref(i__, i__ + 1), lda, &work[1]);
	}
	if (i__ < *m) {
	    i__1 = *m - i__;
	    r__1 = -tau[i__];
	    sscal_(&i__1, &r__1, &a_ref(i__ + 1, i__), &c__1);
	}
	a_ref(i__, i__) = 1.f - tau[i__];

/*        Set A(1:i-1,i) to zero */

	i__1 = i__ - 1;
	for (l = 1; l <= i__1; ++l) {
	    a_ref(l, i__) = 0.f;
/* L30: */
	}
/* L40: */
    }
    return 0;

/*     End of SORG2R */

} /* sorg2r_ */

#undef a_ref





/* Subroutine */ int sorgql_(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    SORGQL generates an M-by-N real matrix Q with orthonormal columns,   
    which is defined as the last N columns of a product of K elementary   
    reflectors of order M   

          Q  =  H(k) . . . H(2) H(1)   

    as returned by SGEQLF.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix Q. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q. M >= N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines the   
            matrix Q. N >= K >= 0.   

    A       (input/output) REAL array, dimension (LDA,N)   
            On entry, the (n-k+i)-th column must contain the vector which   
            defines the elementary reflector H(i), for i = 1,2,...,k, as   
            returned by SGEQLF in the last k columns of its array   
            argument A.   
            On exit, the M-by-N matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= f2cmax(1,M).   

    TAU     (input) REAL array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by SGEQLF.   

    WORK    (workspace/output) REAL array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= f2cmax(1,N).   
            For optimum performance LWORK >= N*NB, where NB is the   
            optimal blocksize.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument has an illegal value   

    =====================================================================   


       Test the input arguments   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__3 = 3;
    static integer c__2 = 2;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    static integer i__, j, l, nbmin, iinfo;
    extern /* Subroutine */ int sorg2l_(integer *, integer *, integer *, real 
	    *, integer *, real *, real *, integer *);
    static integer ib, nb, kk, nx;
    extern /* Subroutine */ int slarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, real *, integer *, real *, 
	    integer *, real *, integer *, real *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int slarft_(char *, char *, integer *, integer *, 
	    real *, integer *, real *, real *, integer *);
    static integer ldwork, lwkopt;
    static logical lquery;
    static integer iws;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    nb = ilaenv_(&c__1, "SORGQL", " ", m, n, k, &c_n1, (ftnlen)6, (ftnlen)1);
    lwkopt = f2cmax(1,*n) * nb;
    work[1] = (real) lwkopt;
    lquery = *lwork == -1;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0 || *n > *m) {
	*info = -2;
    } else if (*k < 0 || *k > *n) {
	*info = -3;
    } else if (*lda < f2cmax(1,*m)) {
	*info = -5;
    } else if (*lwork < f2cmax(1,*n) && ! lquery) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SORGQL", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	work[1] = 1.f;
	return 0;
    }

    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code.   

   Computing MAX */
	i__1 = 0, i__2 = ilaenv_(&c__3, "SORGQL", " ", m, n, k, &c_n1, (
		ftnlen)6, (ftnlen)1);
	nx = f2cmax(i__1,i__2);
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked code. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and   
                determine the minimum value of NB. */

		nb = *lwork / ldwork;
/* Computing MAX */
		i__1 = 2, i__2 = ilaenv_(&c__2, "SORGQL", " ", m, n, k, &c_n1,
			 (ftnlen)6, (ftnlen)1);
		nbmin = f2cmax(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the first block.   
          The last kk columns are handled by the block method.   

   Computing MIN */
	i__1 = *k, i__2 = (*k - nx + nb - 1) / nb * nb;
	kk = f2cmin(i__1,i__2);

/*        Set A(m-kk+1:m,1:n-kk) to zero. */

	i__1 = *n - kk;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = *m - kk + 1; i__ <= i__2; ++i__) {
		a_ref(i__, j) = 0.f;
/* L10: */
	    }
/* L20: */
	}
    } else {
	kk = 0;
    }

/*     Use unblocked code for the first or only block. */

    i__1 = *m - kk;
    i__2 = *n - kk;
    i__3 = *k - kk;
    sorg2l_(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &work[1], &iinfo)
	    ;

    if (kk > 0) {

/*        Use blocked code */

	i__1 = *k;
	i__2 = nb;
	for (i__ = *k - kk + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
		i__2) {
/* Computing MIN */
	    i__3 = nb, i__4 = *k - i__ + 1;
	    ib = f2cmin(i__3,i__4);
	    if (*n - *k + i__ > 1) {

/*              Form the triangular factor of the block reflector   
                H = H(i+ib-1) . . . H(i+1) H(i) */

		i__3 = *m - *k + i__ + ib - 1;
		slarft_("Backward", "Columnwise", &i__3, &ib, &a_ref(1, *n - *
			k + i__), lda, &tau[i__], &work[1], &ldwork);

/*              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left */

		i__3 = *m - *k + i__ + ib - 1;
		i__4 = *n - *k + i__ - 1;
		slarfb_("Left", "No transpose", "Backward", "Columnwise", &
			i__3, &i__4, &ib, &a_ref(1, *n - *k + i__), lda, &
			work[1], &ldwork, &a[a_offset], lda, &work[ib + 1], &
			ldwork);
	    }

/*           Apply H to rows 1:m-k+i+ib-1 of current block */

	    i__3 = *m - *k + i__ + ib - 1;
	    sorg2l_(&i__3, &ib, &ib, &a_ref(1, *n - *k + i__), lda, &tau[i__],
		     &work[1], &iinfo);

/*           Set rows m-k+i+ib:m of current block to zero */

	    i__3 = *n - *k + i__ + ib - 1;
	    for (j = *n - *k + i__; j <= i__3; ++j) {
		i__4 = *m;
		for (l = *m - *k + i__ + ib; l <= i__4; ++l) {
		    a_ref(l, j) = 0.f;
/* L30: */
		}
/* L40: */
	    }
/* L50: */
	}
    }

    work[1] = (real) iws;
    return 0;

/*     End of SORGQL */

} /* sorgql_ */

#undef a_ref





/* Subroutine */ int sorgqr_(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    SORGQR generates an M-by-N real matrix Q with orthonormal columns,   
    which is defined as the first N columns of a product of K elementary   
    reflectors of order M   

          Q  =  H(1) H(2) . . . H(k)   

    as returned by SGEQRF.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix Q. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q. M >= N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines the   
            matrix Q. N >= K >= 0.   

    A       (input/output) REAL array, dimension (LDA,N)   
            On entry, the i-th column must contain the vector which   
            defines the elementary reflector H(i), for i = 1,2,...,k, as   
            returned by SGEQRF in the first k columns of its array   
            argument A.   
            On exit, the M-by-N matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= f2cmax(1,M).   

    TAU     (input) REAL array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by SGEQRF.   

    WORK    (workspace/output) REAL array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= f2cmax(1,N).   
            For optimum performance LWORK >= N*NB, where NB is the   
            optimal blocksize.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument has an illegal value   

    =====================================================================   


       Test the input arguments   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__3 = 3;
    static integer c__2 = 2;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    static integer i__, j, l, nbmin, iinfo, ib;
    extern /* Subroutine */ int sorg2r_(integer *, integer *, integer *, real 
	    *, integer *, real *, real *, integer *);
    static integer nb, ki, kk, nx;
    extern /* Subroutine */ int slarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, real *, integer *, real *, 
	    integer *, real *, integer *, real *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int slarft_(char *, char *, integer *, integer *, 
	    real *, integer *, real *, real *, integer *);
    static integer ldwork, lwkopt;
    static logical lquery;
    static integer iws;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    nb = ilaenv_(&c__1, "SORGQR", " ", m, n, k, &c_n1, (ftnlen)6, (ftnlen)1);
    lwkopt = f2cmax(1,*n) * nb;
    work[1] = (real) lwkopt;
    lquery = *lwork == -1;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0 || *n > *m) {
	*info = -2;
    } else if (*k < 0 || *k > *n) {
	*info = -3;
    } else if (*lda < f2cmax(1,*m)) {
	*info = -5;
    } else if (*lwork < f2cmax(1,*n) && ! lquery) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SORGQR", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	work[1] = 1.f;
	return 0;
    }

    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code.   

   Computing MAX */
	i__1 = 0, i__2 = ilaenv_(&c__3, "SORGQR", " ", m, n, k, &c_n1, (
		ftnlen)6, (ftnlen)1);
	nx = f2cmax(i__1,i__2);
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked code. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and   
                determine the minimum value of NB. */

		nb = *lwork / ldwork;
/* Computing MAX */
		i__1 = 2, i__2 = ilaenv_(&c__2, "SORGQR", " ", m, n, k, &c_n1,
			 (ftnlen)6, (ftnlen)1);
		nbmin = f2cmax(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the last block.   
          The first kk columns are handled by the block method. */

	ki = (*k - nx - 1) / nb * nb;
/* Computing MIN */
	i__1 = *k, i__2 = ki + nb;
	kk = f2cmin(i__1,i__2);

/*        Set A(1:kk,kk+1:n) to zero. */

	i__1 = *n;
	for (j = kk + 1; j <= i__1; ++j) {
	    i__2 = kk;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a_ref(i__, j) = 0.f;
/* L10: */
	    }
/* L20: */
	}
    } else {
	kk = 0;
    }

/*     Use unblocked code for the last or only block. */

    if (kk < *n) {
	i__1 = *m - kk;
	i__2 = *n - kk;
	i__3 = *k - kk;
	sorg2r_(&i__1, &i__2, &i__3, &a_ref(kk + 1, kk + 1), lda, &tau[kk + 1]
		, &work[1], &iinfo);
    }

    if (kk > 0) {

/*        Use blocked code */

	i__1 = -nb;
	for (i__ = ki + 1; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
	    i__2 = nb, i__3 = *k - i__ + 1;
	    ib = f2cmin(i__2,i__3);
	    if (i__ + ib <= *n) {

/*              Form the triangular factor of the block reflector   
                H = H(i) H(i+1) . . . H(i+ib-1) */

		i__2 = *m - i__ + 1;
		slarft_("Forward", "Columnwise", &i__2, &ib, &a_ref(i__, i__),
			 lda, &tau[i__], &work[1], &ldwork);

/*              Apply H to A(i:m,i+ib:n) from the left */

		i__2 = *m - i__ + 1;
		i__3 = *n - i__ - ib + 1;
		slarfb_("Left", "No transpose", "Forward", "Columnwise", &
			i__2, &i__3, &ib, &a_ref(i__, i__), lda, &work[1], &
			ldwork, &a_ref(i__, i__ + ib), lda, &work[ib + 1], &
			ldwork);
	    }

/*           Apply H to rows i:m of current block */

	    i__2 = *m - i__ + 1;
	    sorg2r_(&i__2, &ib, &ib, &a_ref(i__, i__), lda, &tau[i__], &work[
		    1], &iinfo);

/*           Set rows 1:i-1 of current block to zero */

	    i__2 = i__ + ib - 1;
	    for (j = i__; j <= i__2; ++j) {
		i__3 = i__ - 1;
		for (l = 1; l <= i__3; ++l) {
		    a_ref(l, j) = 0.f;
/* L30: */
		}
/* L40: */
	    }
/* L50: */
	}
    }

    work[1] = (real) iws;
    return 0;

/*     End of SORGQR */

} /* sorgqr_ */

#undef a_ref





/* Subroutine */ int sorgtr_(char *uplo, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    SORGTR generates a real orthogonal matrix Q which is defined as the   
    product of n-1 elementary reflectors of order N, as returned by   
    SSYTRD:   

    if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),   

    if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U': Upper triangle of A contains elementary reflectors   
                   from SSYTRD;   
            = 'L': Lower triangle of A contains elementary reflectors   
                   from SSYTRD.   

    N       (input) INTEGER   
            The order of the matrix Q. N >= 0.   

    A       (input/output) REAL array, dimension (LDA,N)   
            On entry, the vectors which define the elementary reflectors,   
            as returned by SSYTRD.   
            On exit, the N-by-N orthogonal matrix Q.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= f2cmax(1,N).   

    TAU     (input) REAL array, dimension (N-1)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by SSYTRD.   

    WORK    (workspace/output) REAL array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= f2cmax(1,N-1).   
            For optimum performance LWORK >= (N-1)*NB, where NB is   
            the optimal blocksize.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    =====================================================================   


       Test the input arguments   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *);
    static integer iinfo;
    static logical upper;
    static integer nb;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int sorgql_(integer *, integer *, integer *, real 
	    *, integer *, real *, real *, integer *, integer *), sorgqr_(
	    integer *, integer *, integer *, real *, integer *, real *, real *
	    , integer *, integer *);
    static logical lquery;
    static integer lwkopt;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < f2cmax(1,*n)) {
	*info = -4;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n - 1;
	if (*lwork < f2cmax(i__1,i__2) && ! lquery) {
	    *info = -7;
	}
    }

    if (*info == 0) {
	if (upper) {
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    i__3 = *n - 1;
	    nb = ilaenv_(&c__1, "SORGQL", " ", &i__1, &i__2, &i__3, &c_n1, (
		    ftnlen)6, (ftnlen)1);
	} else {
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    i__3 = *n - 1;
	    nb = ilaenv_(&c__1, "SORGQR", " ", &i__1, &i__2, &i__3, &c_n1, (
		    ftnlen)6, (ftnlen)1);
	}
/* Computing MAX */
	i__1 = 1, i__2 = *n - 1;
	lwkopt = f2cmax(i__1,i__2) * nb;
	work[1] = (real) lwkopt;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SORGTR", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	work[1] = 1.f;
	return 0;
    }

    if (upper) {

/*        Q was determined by a call to SSYTRD with UPLO = 'U'   

          Shift the vectors which define the elementary reflectors one   
          column to the left, and set the last row and column of Q to   
          those of the unit matrix */

	i__1 = *n - 1;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = j - 1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a_ref(i__, j) = a_ref(i__, j + 1);
/* L10: */
	    }
	    a_ref(*n, j) = 0.f;
/* L20: */
	}
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    a_ref(i__, *n) = 0.f;
/* L30: */
	}
	a_ref(*n, *n) = 1.f;

/*        Generate Q(1:n-1,1:n-1) */

	i__1 = *n - 1;
	i__2 = *n - 1;
	i__3 = *n - 1;
	sorgql_(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &work[1], 
		lwork, &iinfo);

    } else {

/*        Q was determined by a call to SSYTRD with UPLO = 'L'.   

          Shift the vectors which define the elementary reflectors one   
          column to the right, and set the first row and column of Q to   
          those of the unit matrix */

	for (j = *n; j >= 2; --j) {
	    a_ref(1, j) = 0.f;
	    i__1 = *n;
	    for (i__ = j + 1; i__ <= i__1; ++i__) {
		a_ref(i__, j) = a_ref(i__, j - 1);
/* L40: */
	    }
/* L50: */
	}
	a_ref(1, 1) = 1.f;
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    a_ref(i__, 1) = 0.f;
/* L60: */
	}
	if (*n > 1) {

/*           Generate Q(2:n,2:n) */

	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    i__3 = *n - 1;
	    sorgqr_(&i__1, &i__2, &i__3, &a_ref(2, 2), lda, &tau[1], &work[1],
		     lwork, &iinfo);
	}
    }
    work[1] = (real) lwkopt;
    return 0;

/*     End of SORGTR */

} /* sorgtr_ */

#undef a_ref





/* Subroutine */ int sscal_(integer *n, real *sa, real *sx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    /* Local variables */
    static integer i__, m, nincx, mp1;
/*     scales a vector by a constant.   
       uses unrolled loops for increment equal to 1.   
       jack dongarra, linpack, 3/11/78.   
       modified 3/93 to return if incx .le. 0.   
       modified 12/3/93, array(1) declarations changed to array(*)   
       Parameter adjustments */
    --sx;
    /* Function Body */
    if (*n <= 0 || *incx <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }
/*        code for increment not equal to 1 */
    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	sx[i__] = *sa * sx[i__];
/* L10: */
    }
    return 0;
/*        code for increment equal to 1   
          clean-up loop */
L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	sx[i__] = *sa * sx[i__];
/* L30: */
    }
    if (*n < 5) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i__ = mp1; i__ <= i__2; i__ += 5) {
	sx[i__] = *sa * sx[i__];
	sx[i__ + 1] = *sa * sx[i__ + 1];
	sx[i__ + 2] = *sa * sx[i__ + 2];
	sx[i__ + 3] = *sa * sx[i__ + 3];
	sx[i__ + 4] = *sa * sx[i__ + 4];
/* L50: */
    }
    return 0;
} /* sscal_ */




/* Subroutine */ int ssteqr_(char *compz, integer *n, real *d__, real *e, 
	real *z__, integer *ldz, real *work, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SSTEQR computes all eigenvalues and, optionally, eigenvectors of a   
    symmetric tridiagonal matrix using the implicit QL or QR method.   
    The eigenvectors of a full or band symmetric matrix can also be found   
    if SSYTRD or SSPTRD or SSBTRD has been used to reduce this matrix to   
    tridiagonal form.   

    Arguments   
    =========   

    COMPZ   (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only.   
            = 'V':  Compute eigenvalues and eigenvectors of the original   
                    symmetric matrix.  On entry, Z must contain the   
                    orthogonal matrix used to reduce the original matrix   
                    to tridiagonal form.   
            = 'I':  Compute eigenvalues and eigenvectors of the   
                    tridiagonal matrix.  Z is initialized to the identity   
                    matrix.   

    N       (input) INTEGER   
            The order of the matrix.  N >= 0.   

    D       (input/output) REAL array, dimension (N)   
            On entry, the diagonal elements of the tridiagonal matrix.   
            On exit, if INFO = 0, the eigenvalues in ascending order.   

    E       (input/output) REAL array, dimension (N-1)   
            On entry, the (n-1) subdiagonal elements of the tridiagonal   
            matrix.   
            On exit, E has been destroyed.   

    Z       (input/output) REAL array, dimension (LDZ, N)   
            On entry, if  COMPZ = 'V', then Z contains the orthogonal   
            matrix used in the reduction to tridiagonal form.   
            On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the   
            orthonormal eigenvectors of the original symmetric matrix,   
            and if COMPZ = 'I', Z contains the orthonormal eigenvectors   
            of the symmetric tridiagonal matrix.   
            If COMPZ = 'N', then Z is not referenced.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1, and if   
            eigenvectors are desired, then  LDZ >= f2cmax(1,N).   

    WORK    (workspace) REAL array, dimension (f2cmax(1,2*N-2))   
            If COMPZ = 'N', then WORK is not referenced.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  the algorithm has failed to find all the eigenvalues in   
                  a total of 30*N iterations; if INFO = i, then i   
                  elements of E have not converged to zero; on exit, D   
                  and E contain the elements of a symmetric tridiagonal   
                  matrix which is orthogonally similar to the original   
                  matrix.   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
    static real c_b9 = 0.f;
    static real c_b10 = 1.f;
    static integer c__0 = 0;
    static integer c__1 = 1;
    static integer c__2 = 2;
    
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;
    real r__1, r__2;
    /* Builtin functions */
//    double sqrt(doublereal), r_sign(real *, real *);
    double r_sign(real *, real *);
    /* Local variables */
    static integer lend, jtot;
    extern /* Subroutine */ int slae2_(real *, real *, real *, real *, real *)
	    ;
    static real b, c__, f, g;
    static integer i__, j, k, l, m;
    static real p, r__, s;
    extern logical lsame_(char *, char *);
    static real anorm;
    extern /* Subroutine */ int slasr_(char *, char *, char *, integer *, 
	    integer *, real *, real *, real *, integer *);
    static integer l1;
    extern /* Subroutine */ int sswap_(integer *, real *, integer *, real *, 
	    integer *);
    static integer lendm1, lendp1;
    extern /* Subroutine */ int slaev2_(real *, real *, real *, real *, real *
	    , real *, real *);
    extern doublereal slapy2_(real *, real *);
    static integer ii, mm, iscale;
    extern doublereal slamch_(char *);
    static real safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static real safmax;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, real *, 
	    real *, integer *, integer *, real *, integer *, integer *);
    static integer lendsv;
    extern /* Subroutine */ int slartg_(real *, real *, real *, real *, real *
	    ), slaset_(char *, integer *, integer *, real *, real *, real *, 
	    integer *);
    static real ssfmin;
    static integer nmaxit, icompz;
    static real ssfmax;
    extern doublereal slanst_(char *, integer *, real *, real *);
    extern /* Subroutine */ int slasrt_(char *, integer *, real *, integer *);
    static integer lm1, mm1, nm1;
    static real rt1, rt2, eps;
    static integer lsv;
    static real tst, eps2;
#define z___ref(a_1,a_2) z__[(a_2)*z_dim1 + a_1]


    --d__;
    --e;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1 * 1;
    z__ -= z_offset;
    --work;

    /* Function Body */
    *info = 0;

    if (lsame_(compz, "N")) {
	icompz = 0;
    } else if (lsame_(compz, "V")) {
	icompz = 1;
    } else if (lsame_(compz, "I")) {
	icompz = 2;
    } else {
	icompz = -1;
    }
    if (icompz < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ldz < 1 || icompz > 0 && *ldz < f2cmax(1,*n)) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SSTEQR", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (*n == 1) {
	if (icompz == 2) {
	    z___ref(1, 1) = 1.f;
	}
	return 0;
    }

/*     Determine the unit roundoff and over/underflow thresholds. */

    eps = slamch_("E");
/* Computing 2nd power */
    r__1 = eps;
    eps2 = r__1 * r__1;
    safmin = slamch_("S");
    safmax = 1.f / safmin;
    ssfmax = sqrt(safmax) / 3.f;
    ssfmin = sqrt(safmin) / eps2;

/*     Compute the eigenvalues and eigenvectors of the tridiagonal   
       matrix. */

    if (icompz == 2) {
	slaset_("Full", n, n, &c_b9, &c_b10, &z__[z_offset], ldz);
    }

    nmaxit = *n * 30;
    jtot = 0;

/*     Determine where the matrix splits and choose QL or QR iteration   
       for each block, according to whether top or bottom diagonal   
       element is smaller. */

    l1 = 1;
    nm1 = *n - 1;

L10:
    if (l1 > *n) {
	goto L160;
    }
    if (l1 > 1) {
	e[l1 - 1] = 0.f;
    }
    if (l1 <= nm1) {
	i__1 = nm1;
	for (m = l1; m <= i__1; ++m) {
	    tst = (r__1 = e[m], dabs(r__1));
	    if (tst == 0.f) {
		goto L30;
	    }
	    if (tst <= sqrt((r__1 = d__[m], dabs(r__1))) * sqrt((r__2 = d__[m 
		    + 1], dabs(r__2))) * eps) {
		e[m] = 0.f;
		goto L30;
	    }
/* L20: */
	}
    }
    m = *n;

L30:
    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m + 1;
    if (lend == l) {
	goto L10;
    }

/*     Scale submatrix in rows and columns L to LEND */

    i__1 = lend - l + 1;
    anorm = slanst_("I", &i__1, &d__[l], &e[l]);
    iscale = 0;
    if (anorm == 0.f) {
	goto L10;
    }
    if (anorm > ssfmax) {
	iscale = 1;
	i__1 = lend - l + 1;
	slascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n, 
		info);
	i__1 = lend - l;
	slascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, 
		info);
    } else if (anorm < ssfmin) {
	iscale = 2;
	i__1 = lend - l + 1;
	slascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n, 
		info);
	i__1 = lend - l;
	slascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, 
		info);
    }

/*     Choose between QL and QR iteration */

    if ((r__1 = d__[lend], dabs(r__1)) < (r__2 = d__[l], dabs(r__2))) {
	lend = lsv;
	l = lendsv;
    }

    if (lend > l) {

/*        QL Iteration   

          Look for small subdiagonal element. */

L40:
	if (l != lend) {
	    lendm1 = lend - 1;
	    i__1 = lendm1;
	    for (m = l; m <= i__1; ++m) {
/* Computing 2nd power */
		r__2 = (r__1 = e[m], dabs(r__1));
		tst = r__2 * r__2;
		if (tst <= eps2 * (r__1 = d__[m], dabs(r__1)) * (r__2 = d__[m 
			+ 1], dabs(r__2)) + safmin) {
		    goto L60;
		}
/* L50: */
	    }
	}

	m = lend;

L60:
	if (m < lend) {
	    e[m] = 0.f;
	}
	p = d__[l];
	if (m == l) {
	    goto L80;
	}

/*        If remaining matrix is 2-by-2, use SLAE2 or SLAEV2   
          to compute its eigensystem. */

	if (m == l + 1) {
	    if (icompz > 0) {
		slaev2_(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2, &c__, &s);
		work[l] = c__;
		work[*n - 1 + l] = s;
		slasr_("R", "V", "B", n, &c__2, &work[l], &work[*n - 1 + l], &
			z___ref(1, l), ldz);
	    } else {
		slae2_(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2);
	    }
	    d__[l] = rt1;
	    d__[l + 1] = rt2;
	    e[l] = 0.f;
	    l += 2;
	    if (l <= lend) {
		goto L40;
	    }
	    goto L140;
	}

	if (jtot == nmaxit) {
	    goto L140;
	}
	++jtot;

/*        Form shift. */

	g = (d__[l + 1] - p) / (e[l] * 2.f);
	r__ = slapy2_(&g, &c_b10);
	g = d__[m] - p + e[l] / (g + r_sign(&r__, &g));

	s = 1.f;
	c__ = 1.f;
	p = 0.f;

/*        Inner loop */

	mm1 = m - 1;
	i__1 = l;
	for (i__ = mm1; i__ >= i__1; --i__) {
	    f = s * e[i__];
	    b = c__ * e[i__];
	    slartg_(&g, &f, &c__, &s, &r__);
	    if (i__ != m - 1) {
		e[i__ + 1] = r__;
	    }
	    g = d__[i__ + 1] - p;
	    r__ = (d__[i__] - g) * s + c__ * 2.f * b;
	    p = s * r__;
	    d__[i__ + 1] = g + p;
	    g = c__ * r__ - b;

/*           If eigenvectors are desired, then save rotations. */

	    if (icompz > 0) {
		work[i__] = c__;
		work[*n - 1 + i__] = -s;
	    }

/* L70: */
	}

/*        If eigenvectors are desired, then apply saved rotations. */

	if (icompz > 0) {
	    mm = m - l + 1;
	    slasr_("R", "V", "B", n, &mm, &work[l], &work[*n - 1 + l], &
		    z___ref(1, l), ldz);
	}

	d__[l] -= p;
	e[l] = g;
	goto L40;

/*        Eigenvalue found. */

L80:
	d__[l] = p;

	++l;
	if (l <= lend) {
	    goto L40;
	}
	goto L140;

    } else {

/*        QR Iteration   

          Look for small superdiagonal element. */

L90:
	if (l != lend) {
	    lendp1 = lend + 1;
	    i__1 = lendp1;
	    for (m = l; m >= i__1; --m) {
/* Computing 2nd power */
		r__2 = (r__1 = e[m - 1], dabs(r__1));
		tst = r__2 * r__2;
		if (tst <= eps2 * (r__1 = d__[m], dabs(r__1)) * (r__2 = d__[m 
			- 1], dabs(r__2)) + safmin) {
		    goto L110;
		}
/* L100: */
	    }
	}

	m = lend;

L110:
	if (m > lend) {
	    e[m - 1] = 0.f;
	}
	p = d__[l];
	if (m == l) {
	    goto L130;
	}

/*        If remaining matrix is 2-by-2, use SLAE2 or SLAEV2   
          to compute its eigensystem. */

	if (m == l - 1) {
	    if (icompz > 0) {
		slaev2_(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2, &c__, &s)
			;
		work[m] = c__;
		work[*n - 1 + m] = s;
		slasr_("R", "V", "F", n, &c__2, &work[m], &work[*n - 1 + m], &
			z___ref(1, l - 1), ldz);
	    } else {
		slae2_(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2);
	    }
	    d__[l - 1] = rt1;
	    d__[l] = rt2;
	    e[l - 1] = 0.f;
	    l += -2;
	    if (l >= lend) {
		goto L90;
	    }
	    goto L140;
	}

	if (jtot == nmaxit) {
	    goto L140;
	}
	++jtot;

/*        Form shift. */

	g = (d__[l - 1] - p) / (e[l - 1] * 2.f);
	r__ = slapy2_(&g, &c_b10);
	g = d__[m] - p + e[l - 1] / (g + r_sign(&r__, &g));

	s = 1.f;
	c__ = 1.f;
	p = 0.f;

/*        Inner loop */

	lm1 = l - 1;
	i__1 = lm1;
	for (i__ = m; i__ <= i__1; ++i__) {
	    f = s * e[i__];
	    b = c__ * e[i__];
	    slartg_(&g, &f, &c__, &s, &r__);
	    if (i__ != m) {
		e[i__ - 1] = r__;
	    }
	    g = d__[i__] - p;
	    r__ = (d__[i__ + 1] - g) * s + c__ * 2.f * b;
	    p = s * r__;
	    d__[i__] = g + p;
	    g = c__ * r__ - b;

/*           If eigenvectors are desired, then save rotations. */

	    if (icompz > 0) {
		work[i__] = c__;
		work[*n - 1 + i__] = s;
	    }

/* L120: */
	}

/*        If eigenvectors are desired, then apply saved rotations. */

	if (icompz > 0) {
	    mm = l - m + 1;
	    slasr_("R", "V", "F", n, &mm, &work[m], &work[*n - 1 + m], &
		    z___ref(1, m), ldz);
	}

	d__[l] -= p;
	e[lm1] = g;
	goto L90;

/*        Eigenvalue found. */

L130:
	d__[l] = p;

	--l;
	if (l >= lend) {
	    goto L90;
	}
	goto L140;

    }

/*     Undo scaling if necessary */

L140:
    if (iscale == 1) {
	i__1 = lendsv - lsv + 1;
	slascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info);
	i__1 = lendsv - lsv;
	slascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &e[lsv], n, 
		info);
    } else if (iscale == 2) {
	i__1 = lendsv - lsv + 1;
	slascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info);
	i__1 = lendsv - lsv;
	slascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &e[lsv], n, 
		info);
    }

/*     Check for no convergence to an eigenvalue after a total   
       of N*MAXIT iterations. */

    if (jtot < nmaxit) {
	goto L10;
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (e[i__] != 0.f) {
	    ++(*info);
	}
/* L150: */
    }
    goto L190;

/*     Order eigenvalues and eigenvectors. */

L160:
    if (icompz == 0) {

/*        Use Quick Sort */

	slasrt_("I", n, &d__[1], info);

    } else {

/*        Use Selection Sort to minimize swaps of eigenvectors */

	i__1 = *n;
	for (ii = 2; ii <= i__1; ++ii) {
	    i__ = ii - 1;
	    k = i__;
	    p = d__[i__];
	    i__2 = *n;
	    for (j = ii; j <= i__2; ++j) {
		if (d__[j] < p) {
		    k = j;
		    p = d__[j];
		}
/* L170: */
	    }
	    if (k != i__) {
		d__[k] = d__[i__];
		d__[i__] = p;
		sswap_(n, &z___ref(1, i__), &c__1, &z___ref(1, k), &c__1);
	    }
/* L180: */
	}
    }

L190:
    return 0;

/*     End of SSTEQR */

} /* ssteqr_ */

#undef z___ref





/* Subroutine */ int ssterf_(integer *n, real *d__, real *e, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    SSTERF computes all eigenvalues of a symmetric tridiagonal matrix   
    using the Pal-Walker-Kahan variant of the QL or QR algorithm.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix.  N >= 0.   

    D       (input/output) REAL array, dimension (N)   
            On entry, the n diagonal elements of the tridiagonal matrix.   
            On exit, if INFO = 0, the eigenvalues in ascending order.   

    E       (input/output) REAL array, dimension (N-1)   
            On entry, the (n-1) subdiagonal elements of the tridiagonal   
            matrix.   
            On exit, E has been destroyed.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  the algorithm failed to find all of the eigenvalues in   
                  a total of 30*N iterations; if INFO = i, then i   
                  elements of E have not converged to zero.   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__0 = 0;
    static integer c__1 = 1;
    static real c_b32 = 1.f;
    
    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3;
    /* Builtin functions */
//    double sqrt(doublereal), r_sign(real *, real *);
    double r_sign(real *, real *);
    /* Local variables */
    static real oldc;
    static integer lend, jtot;
    extern /* Subroutine */ int slae2_(real *, real *, real *, real *, real *)
	    ;
    static real c__;
    static integer i__, l, m;
    static real p, gamma, r__, s, alpha, sigma, anorm;
    static integer l1;
    static real bb;
    extern doublereal slapy2_(real *, real *);
    static integer iscale;
    static real oldgam;
    extern doublereal slamch_(char *);
    static real safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static real safmax;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, real *, 
	    real *, integer *, integer *, real *, integer *, integer *);
    static integer lendsv;
    static real ssfmin;
    static integer nmaxit;
    static real ssfmax;
    extern doublereal slanst_(char *, integer *, real *, real *);
    extern /* Subroutine */ int slasrt_(char *, integer *, real *, integer *);
    static real rt1, rt2, eps, rte;
    static integer lsv;
    static real eps2;


    --e;
    --d__;

    /* Function Body */
    *info = 0;

/*     Quick return if possible */

    if (*n < 0) {
	*info = -1;
	i__1 = -(*info);
	xerbla_("SSTERF", &i__1);
	return 0;
    }
    if (*n <= 1) {
	return 0;
    }

/*     Determine the unit roundoff for this environment. */

    eps = slamch_("E");
/* Computing 2nd power */
    r__1 = eps;
    eps2 = r__1 * r__1;
    safmin = slamch_("S");
    safmax = 1.f / safmin;
    ssfmax = sqrt(safmax) / 3.f;
    ssfmin = sqrt(safmin) / eps2;

/*     Compute the eigenvalues of the tridiagonal matrix. */

    nmaxit = *n * 30;
    sigma = 0.f;
    jtot = 0;

/*     Determine where the matrix splits and choose QL or QR iteration   
       for each block, according to whether top or bottom diagonal   
       element is smaller. */

    l1 = 1;

L10:
    if (l1 > *n) {
	goto L170;
    }
    if (l1 > 1) {
	e[l1 - 1] = 0.f;
    }
    i__1 = *n - 1;
    for (m = l1; m <= i__1; ++m) {
	if ((r__3 = e[m], dabs(r__3)) <= sqrt((r__1 = d__[m], dabs(r__1))) * 
		sqrt((r__2 = d__[m + 1], dabs(r__2))) * eps) {
	    e[m] = 0.f;
	    goto L30;
	}
/* L20: */
    }
    m = *n;

L30:
    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m + 1;
    if (lend == l) {
	goto L10;
    }

/*     Scale submatrix in rows and columns L to LEND */

    i__1 = lend - l + 1;
    anorm = slanst_("I", &i__1, &d__[l], &e[l]);
    iscale = 0;
    if (anorm > ssfmax) {
	iscale = 1;
	i__1 = lend - l + 1;
	slascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n, 
		info);
	i__1 = lend - l;
	slascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, 
		info);
    } else if (anorm < ssfmin) {
	iscale = 2;
	i__1 = lend - l + 1;
	slascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n, 
		info);
	i__1 = lend - l;
	slascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, 
		info);
    }

    i__1 = lend - 1;
    for (i__ = l; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	r__1 = e[i__];
	e[i__] = r__1 * r__1;
/* L40: */
    }

/*     Choose between QL and QR iteration */

    if ((r__1 = d__[lend], dabs(r__1)) < (r__2 = d__[l], dabs(r__2))) {
	lend = lsv;
	l = lendsv;
    }

    if (lend >= l) {

/*        QL Iteration   

          Look for small subdiagonal element. */

L50:
	if (l != lend) {
	    i__1 = lend - 1;
	    for (m = l; m <= i__1; ++m) {
		if ((r__2 = e[m], dabs(r__2)) <= eps2 * (r__1 = d__[m] * d__[
			m + 1], dabs(r__1))) {
		    goto L70;
		}
/* L60: */
	    }
	}
	m = lend;

L70:
	if (m < lend) {
	    e[m] = 0.f;
	}
	p = d__[l];
	if (m == l) {
	    goto L90;
	}

/*        If remaining matrix is 2 by 2, use SLAE2 to compute its   
          eigenvalues. */

	if (m == l + 1) {
	    rte = sqrt(e[l]);
	    slae2_(&d__[l], &rte, &d__[l + 1], &rt1, &rt2);
	    d__[l] = rt1;
	    d__[l + 1] = rt2;
	    e[l] = 0.f;
	    l += 2;
	    if (l <= lend) {
		goto L50;
	    }
	    goto L150;
	}

	if (jtot == nmaxit) {
	    goto L150;
	}
	++jtot;

/*        Form shift. */

	rte = sqrt(e[l]);
	sigma = (d__[l + 1] - p) / (rte * 2.f);
	r__ = slapy2_(&sigma, &c_b32);
	sigma = p - rte / (sigma + r_sign(&r__, &sigma));

	c__ = 1.f;
	s = 0.f;
	gamma = d__[m] - sigma;
	p = gamma * gamma;

/*        Inner loop */

	i__1 = l;
	for (i__ = m - 1; i__ >= i__1; --i__) {
	    bb = e[i__];
	    r__ = p + bb;
	    if (i__ != m - 1) {
		e[i__ + 1] = s * r__;
	    }
	    oldc = c__;
	    c__ = p / r__;
	    s = bb / r__;
	    oldgam = gamma;
	    alpha = d__[i__];
	    gamma = c__ * (alpha - sigma) - s * oldgam;
	    d__[i__ + 1] = oldgam + (alpha - gamma);
	    if (c__ != 0.f) {
		p = gamma * gamma / c__;
	    } else {
		p = oldc * bb;
	    }
/* L80: */
	}

	e[l] = s * p;
	d__[l] = sigma + gamma;
	goto L50;

/*        Eigenvalue found. */

L90:
	d__[l] = p;

	++l;
	if (l <= lend) {
	    goto L50;
	}
	goto L150;

    } else {

/*        QR Iteration   

          Look for small superdiagonal element. */

L100:
	i__1 = lend + 1;
	for (m = l; m >= i__1; --m) {
	    if ((r__2 = e[m - 1], dabs(r__2)) <= eps2 * (r__1 = d__[m] * d__[
		    m - 1], dabs(r__1))) {
		goto L120;
	    }
/* L110: */
	}
	m = lend;

L120:
	if (m > lend) {
	    e[m - 1] = 0.f;
	}
	p = d__[l];
	if (m == l) {
	    goto L140;
	}

/*        If remaining matrix is 2 by 2, use SLAE2 to compute its   
          eigenvalues. */

	if (m == l - 1) {
	    rte = sqrt(e[l - 1]);
	    slae2_(&d__[l], &rte, &d__[l - 1], &rt1, &rt2);
	    d__[l] = rt1;
	    d__[l - 1] = rt2;
	    e[l - 1] = 0.f;
	    l += -2;
	    if (l >= lend) {
		goto L100;
	    }
	    goto L150;
	}

	if (jtot == nmaxit) {
	    goto L150;
	}
	++jtot;

/*        Form shift. */

	rte = sqrt(e[l - 1]);
	sigma = (d__[l - 1] - p) / (rte * 2.f);
	r__ = slapy2_(&sigma, &c_b32);
	sigma = p - rte / (sigma + r_sign(&r__, &sigma));

	c__ = 1.f;
	s = 0.f;
	gamma = d__[m] - sigma;
	p = gamma * gamma;

/*        Inner loop */

	i__1 = l - 1;
	for (i__ = m; i__ <= i__1; ++i__) {
	    bb = e[i__];
	    r__ = p + bb;
	    if (i__ != m) {
		e[i__ - 1] = s * r__;
	    }
	    oldc = c__;
	    c__ = p / r__;
	    s = bb / r__;
	    oldgam = gamma;
	    alpha = d__[i__ + 1];
	    gamma = c__ * (alpha - sigma) - s * oldgam;
	    d__[i__] = oldgam + (alpha - gamma);
	    if (c__ != 0.f) {
		p = gamma * gamma / c__;
	    } else {
		p = oldc * bb;
	    }
/* L130: */
	}

	e[l - 1] = s * p;
	d__[l] = sigma + gamma;
	goto L100;

/*        Eigenvalue found. */

L140:
	d__[l] = p;

	--l;
	if (l >= lend) {
	    goto L100;
	}
	goto L150;

    }

/*     Undo scaling if necessary */

L150:
    if (iscale == 1) {
	i__1 = lendsv - lsv + 1;
	slascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info);
    }
    if (iscale == 2) {
	i__1 = lendsv - lsv + 1;
	slascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info);
    }

/*     Check for no convergence to an eigenvalue after a total   
       of N*MAXIT iterations. */

    if (jtot < nmaxit) {
	goto L10;
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (e[i__] != 0.f) {
	    ++(*info);
	}
/* L160: */
    }
    goto L180;

/*     Sort eigenvalues in increasing order. */

L170:
    slasrt_("I", n, &d__[1], info);

L180:
    return 0;

/*     End of SSTERF */

} /* ssterf_ */




/* Subroutine */ int sswap_(integer *n, real *sx, integer *incx, real *sy, 
	integer *incy)
{
    /* System generated locals */
    integer i__1;
    /* Local variables */
    static integer i__, m;
    static real stemp;
    static integer ix, iy, mp1;
/*     interchanges two vectors.   
       uses unrolled loops for increments equal to 1.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   
       Parameter adjustments */
    --sy;
    --sx;
    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }
/*       code for unequal increments or equal increments not equal   
           to 1 */
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stemp = sx[ix];
	sx[ix] = sy[iy];
	sy[iy] = stemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;
/*       code for both increments equal to 1   
         clean-up loop */
L20:
    m = *n % 3;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stemp = sx[i__];
	sx[i__] = sy[i__];
	sy[i__] = stemp;
/* L30: */
    }
    if (*n < 3) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 3) {
	stemp = sx[i__];
	sx[i__] = sy[i__];
	sy[i__] = stemp;
	stemp = sx[i__ + 1];
	sx[i__ + 1] = sy[i__ + 1];
	sy[i__ + 1] = stemp;
	stemp = sx[i__ + 2];
	sx[i__ + 2] = sy[i__ + 2];
	sy[i__ + 2] = stemp;
/* L50: */
    }
    return 0;
} /* sswap_ */




/* Subroutine */ int ssyev_(char *jobz, char *uplo, integer *n, real *a, 
	integer *lda, real *w, real *work, integer *lwork, integer *info)
{
/*  -- LAPACK driver routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    SSYEV computes all eigenvalues and, optionally, eigenvectors of a   
    real symmetric matrix A.   

    Arguments   
    =========   

    JOBZ    (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only;   
            = 'V':  Compute eigenvalues and eigenvectors.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) REAL array, dimension (LDA, N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the   
            leading N-by-N upper triangular part of A contains the   
            upper triangular part of the matrix A.  If UPLO = 'L',   
            the leading N-by-N lower triangular part of A contains   
            the lower triangular part of the matrix A.   
            On exit, if JOBZ = 'V', then if INFO = 0, A contains the   
            orthonormal eigenvectors of the matrix A.   
            If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')   
            or the upper triangle (if UPLO='U') of A, including the   
            diagonal, is destroyed.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= f2cmax(1,N).   

    W       (output) REAL array, dimension (N)   
            If INFO = 0, the eigenvalues in ascending order.   

    WORK    (workspace/output) REAL array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The length of the array WORK.  LWORK >= f2cmax(1,3*N-1).   
            For optimal efficiency, LWORK >= (NB+2)*N,   
            where NB is the blocksize for SSYTRD returned by ILAENV.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the algorithm failed to converge; i   
                  off-diagonal elements of an intermediate tridiagonal   
                  form did not converge to zero.   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__0 = 0;
    static real c_b17 = 1.f;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    real r__1;
    /* Builtin functions */
//    double sqrt(doublereal);
    /* Local variables */
    static integer inde;
    static real anrm;
    static integer imax;
    static real rmin, rmax;
    static integer lopt;
    static real sigma;
    extern logical lsame_(char *, char *);
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    static logical lower, wantz;
    static integer nb, iscale;
    extern doublereal slamch_(char *);
    static real safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static real bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, real *, 
	    real *, integer *, integer *, real *, integer *, integer *);
    static integer indtau, indwrk;
    extern /* Subroutine */ int ssterf_(integer *, real *, real *, integer *);
    extern doublereal slansy_(char *, char *, integer *, real *, integer *, 
	    real *);
    static integer llwork;
    static real smlnum;
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int sorgtr_(char *, integer *, real *, integer *, 
	    real *, real *, integer *, integer *), ssteqr_(char *, 
	    integer *, real *, real *, real *, integer *, real *, integer *), ssytrd_(char *, integer *, real *, integer *, real *, 
	    real *, real *, real *, integer *, integer *);
    static real eps;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --w;
    --work;

    /* Function Body */
    wantz = lsame_(jobz, "V");
    lower = lsame_(uplo, "L");
    lquery = *lwork == -1;

    *info = 0;
    if (! (wantz || lsame_(jobz, "N"))) {
	*info = -1;
    } else if (! (lower || lsame_(uplo, "U"))) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < f2cmax(1,*n)) {
	*info = -5;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n * 3 - 1;
	if (*lwork < f2cmax(i__1,i__2) && ! lquery) {
	    *info = -8;
	}
    }

    if (*info == 0) {
	nb = ilaenv_(&c__1, "SSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
/* Computing MAX */
	i__1 = 1, i__2 = (nb + 2) * *n;
	lwkopt = f2cmax(i__1,i__2);
	work[1] = (real) lwkopt;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SSYEV ", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	work[1] = 1.f;
	return 0;
    }

    if (*n == 1) {
	w[1] = a_ref(1, 1);
	work[1] = 3.f;
	if (wantz) {
	    a_ref(1, 1) = 1.f;
	}
	return 0;
    }

/*     Get machine constants. */

    safmin = slamch_("Safe minimum");
    eps = slamch_("Precision");
    smlnum = safmin / eps;
    bignum = 1.f / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

    anrm = slansy_("M", uplo, n, &a[a_offset], lda, &work[1]);
    iscale = 0;
    if (anrm > 0.f && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm;
    }
    if (iscale == 1) {
	slascl_(uplo, &c__0, &c__0, &c_b17, &sigma, n, n, &a[a_offset], lda, 
		info);
    }

/*     Call SSYTRD to reduce symmetric matrix to tridiagonal form. */

    inde = 1;
    indtau = inde + *n;
    indwrk = indtau + *n;
    llwork = *lwork - indwrk + 1;
    ssytrd_(uplo, n, &a[a_offset], lda, &w[1], &work[inde], &work[indtau], &
	    work[indwrk], &llwork, &iinfo);
    lopt = (*n << 1) + work[indwrk];

/*     For eigenvalues only, call SSTERF.  For eigenvectors, first call   
       SORGTR to generate the orthogonal matrix, then call SSTEQR. */

    if (! wantz) {
	ssterf_(n, &w[1], &work[inde], info);
    } else {
	sorgtr_(uplo, n, &a[a_offset], lda, &work[indtau], &work[indwrk], &
		llwork, &iinfo);
	ssteqr_(jobz, n, &w[1], &work[inde], &a[a_offset], lda, &work[indtau],
		 info);
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

    if (iscale == 1) {
	if (*info == 0) {
	    imax = *n;
	} else {
	    imax = *info - 1;
	}
	r__1 = 1.f / sigma;
	sscal_(&imax, &r__1, &w[1], &c__1);
    }

/*     Set WORK(1) to optimal workspace size. */

    work[1] = (real) lwkopt;

    return 0;

/*     End of SSYEV */

} /* ssyev_ */

#undef a_ref





/* Subroutine */ int ssymv_(char *uplo, integer *n, real *alpha, real *a, 
	integer *lda, real *x, integer *incx, real *beta, real *y, integer *
	incy)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    static integer info;
    static real temp1, temp2;
    static integer i__, j;
    extern logical lsame_(char *, char *);
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_(char *, integer *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
/*  Purpose   
    =======   
    SSYMV  performs the matrix-vector  operation   
       y := alpha*A*x + beta*y,   
    where alpha and beta are scalars, x and y are n element vectors and   
    A is an n by n symmetric matrix.   
    Parameters   
    ==========   
    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the array A is to be referenced as   
             follows:   
                UPLO = 'U' or 'u'   Only the upper triangular part of A   
                                    is to be referenced.   
                UPLO = 'L' or 'l'   Only the lower triangular part of A   
                                    is to be referenced.   
             Unchanged on exit.   
    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   
    ALPHA  - REAL            .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   
    A      - REAL             array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper   
             triangular part of the symmetric matrix and the strictly   
             lower triangular part of A is not referenced.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower   
             triangular part of the symmetric matrix and the strictly   
             upper triangular part of A is not referenced.   
             Unchanged on exit.   
    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared   
             in the calling (sub) program. LDA must be at least   
             f2cmax( 1, n ).   
             Unchanged on exit.   
    X      - REAL             array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element vector x.   
             Unchanged on exit.   
    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   
    BETA   - REAL            .   
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   
             Unchanged on exit.   
    Y      - REAL             array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y. On exit, Y is overwritten by the updated   
             vector y.   
    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   
    Level 2 Blas routine.   
    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   
       Test the input parameters.   
       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --x;
    --y;
    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*lda < f2cmax(1,*n)) {
	info = 5;
    } else if (*incx == 0) {
	info = 7;
    } else if (*incy == 0) {
	info = 10;
    }
    if (info != 0) {
	xerbla_("SSYMV ", &info);
	return 0;
    }
/*     Quick return if possible. */
    if (*n == 0 || *alpha == 0.f && *beta == 1.f) {
	return 0;
    }
/*     Set up the start points in  X  and  Y. */
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (*n - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (*n - 1) * *incy;
    }
/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through the triangular part   
       of A.   
       First form  y := beta*y. */
    if (*beta != 1.f) {
	if (*incy == 1) {
	    if (*beta == 0.f) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = 0.f;
/* L10: */
		}
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = *beta * y[i__];
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.f) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = 0.f;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = *beta * y[iy];
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == 0.f) {
	return 0;
    }
    if (lsame_(uplo, "U")) {
/*        Form  y  when A is stored in upper triangle. */
	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[j];
		temp2 = 0.f;
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    y[i__] += temp1 * a_ref(i__, j);
		    temp2 += a_ref(i__, j) * x[i__];
/* L50: */
		}
		y[j] = y[j] + temp1 * a_ref(j, j) + *alpha * temp2;
/* L60: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[jx];
		temp2 = 0.f;
		ix = kx;
		iy = ky;
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    y[iy] += temp1 * a_ref(i__, j);
		    temp2 += a_ref(i__, j) * x[ix];
		    ix += *incx;
		    iy += *incy;
/* L70: */
		}
		y[jy] = y[jy] + temp1 * a_ref(j, j) + *alpha * temp2;
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    } else {
/*        Form  y  when A is stored in lower triangle. */
	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[j];
		temp2 = 0.f;
		y[j] += temp1 * a_ref(j, j);
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    y[i__] += temp1 * a_ref(i__, j);
		    temp2 += a_ref(i__, j) * x[i__];
/* L90: */
		}
		y[j] += *alpha * temp2;
/* L100: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[jx];
		temp2 = 0.f;
		y[jy] += temp1 * a_ref(j, j);
		ix = jx;
		iy = jy;
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    ix += *incx;
		    iy += *incy;
		    y[iy] += temp1 * a_ref(i__, j);
		    temp2 += a_ref(i__, j) * x[ix];
/* L110: */
		}
		y[jy] += *alpha * temp2;
		jx += *incx;
		jy += *incy;
/* L120: */
	    }
	}
    }
    return 0;
/*     End of SSYMV . */
} /* ssymv_ */
#undef a_ref




/* Subroutine */ int ssyr2_(char *uplo, integer *n, real *alpha, real *x, 
	integer *incx, real *y, integer *incy, real *a, integer *lda)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    static integer info;
    static real temp1, temp2;
    static integer i__, j;
    extern logical lsame_(char *, char *);
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_(char *, integer *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
/*  Purpose   
    =======   
    SSYR2  performs the symmetric rank 2 operation   
       A := alpha*x*y' + alpha*y*x' + A,   
    where alpha is a scalar, x and y are n element vectors and A is an n   
    by n symmetric matrix.   
    Parameters   
    ==========   
    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the array A is to be referenced as   
             follows:   
                UPLO = 'U' or 'u'   Only the upper triangular part of A   
                                    is to be referenced.   
                UPLO = 'L' or 'l'   Only the lower triangular part of A   
                                    is to be referenced.   
             Unchanged on exit.   
    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   
    ALPHA  - REAL            .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   
    X      - REAL             array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element vector x.   
             Unchanged on exit.   
    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   
    Y      - REAL             array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y.   
             Unchanged on exit.   
    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   
    A      - REAL             array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper   
             triangular part of the symmetric matrix and the strictly   
             lower triangular part of A is not referenced. On exit, the   
             upper triangular part of the array A is overwritten by the   
             upper triangular part of the updated matrix.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower   
             triangular part of the symmetric matrix and the strictly   
             upper triangular part of A is not referenced. On exit, the   
             lower triangular part of the array A is overwritten by the   
             lower triangular part of the updated matrix.   
    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared   
             in the calling (sub) program. LDA must be at least   
             f2cmax( 1, n ).   
             Unchanged on exit.   
    Level 2 Blas routine.   
    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   
       Test the input parameters.   
       Parameter adjustments */
    --x;
    --y;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < f2cmax(1,*n)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("SSYR2 ", &info);
	return 0;
    }
/*     Quick return if possible. */
    if (*n == 0 || *alpha == 0.f) {
	return 0;
    }
/*     Set up the start points in X and Y if the increments are not both   
       unity. */
    if (*incx != 1 || *incy != 1) {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*n - 1) * *incx;
	}
	if (*incy > 0) {
	    ky = 1;
	} else {
	    ky = 1 - (*n - 1) * *incy;
	}
	jx = kx;
	jy = ky;
    }
/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through the triangular part   
       of A. */
    if (lsame_(uplo, "U")) {
/*        Form  A  when A is stored in the upper triangle. */
	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != 0.f || y[j] != 0.f) {
		    temp1 = *alpha * y[j];
		    temp2 = *alpha * x[j];
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			a_ref(i__, j) = a_ref(i__, j) + x[i__] * temp1 + y[
				i__] * temp2;
/* L10: */
		    }
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0.f || y[jy] != 0.f) {
		    temp1 = *alpha * y[jy];
		    temp2 = *alpha * x[jx];
		    ix = kx;
		    iy = ky;
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			a_ref(i__, j) = a_ref(i__, j) + x[ix] * temp1 + y[iy] 
				* temp2;
			ix += *incx;
			iy += *incy;
/* L30: */
		    }
		}
		jx += *incx;
		jy += *incy;
/* L40: */
	    }
	}
    } else {
/*        Form  A  when A is stored in the lower triangle. */
	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != 0.f || y[j] != 0.f) {
		    temp1 = *alpha * y[j];
		    temp2 = *alpha * x[j];
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			a_ref(i__, j) = a_ref(i__, j) + x[i__] * temp1 + y[
				i__] * temp2;
/* L50: */
		    }
		}
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0.f || y[jy] != 0.f) {
		    temp1 = *alpha * y[jy];
		    temp2 = *alpha * x[jx];
		    ix = jx;
		    iy = jy;
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			a_ref(i__, j) = a_ref(i__, j) + x[ix] * temp1 + y[iy] 
				* temp2;
			ix += *incx;
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    }
    return 0;
/*     End of SSYR2 . */
} /* ssyr2_ */
#undef a_ref




/* Subroutine */ int ssyr2k_(char *uplo, char *trans, integer *n, integer *k, 
	real *alpha, real *a, integer *lda, real *b, integer *ldb, real *beta,
	 real *c__, integer *ldc)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;
    /* Local variables */
    static integer info;
    static real temp1, temp2;
    static integer i__, j, l;
    extern logical lsame_(char *, char *);
    static integer nrowa;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
#define c___ref(a_1,a_2) c__[(a_2)*c_dim1 + a_1]
/*  Purpose   
    =======   
    SSYR2K  performs one of the symmetric rank 2k operations   
       C := alpha*A*B' + alpha*B*A' + beta*C,   
    or   
       C := alpha*A'*B + alpha*B'*A + beta*C,   
    where  alpha and beta  are scalars, C is an  n by n  symmetric matrix   
    and  A and B  are  n by k  matrices  in the  first  case  and  k by n   
    matrices in the second case.   
    Parameters   
    ==========   
    UPLO   - CHARACTER*1.   
             On  entry,   UPLO  specifies  whether  the  upper  or  lower   
             triangular  part  of the  array  C  is to be  referenced  as   
             follows:   
                UPLO = 'U' or 'u'   Only the  upper triangular part of  C   
                                    is to be referenced.   
                UPLO = 'L' or 'l'   Only the  lower triangular part of  C   
                                    is to be referenced.   
             Unchanged on exit.   
    TRANS  - CHARACTER*1.   
             On entry,  TRANS  specifies the operation to be performed as   
             follows:   
                TRANS = 'N' or 'n'   C := alpha*A*B' + alpha*B*A' +   
                                          beta*C.   
                TRANS = 'T' or 't'   C := alpha*A'*B + alpha*B'*A +   
                                          beta*C.   
                TRANS = 'C' or 'c'   C := alpha*A'*B + alpha*B'*A +   
                                          beta*C.   
             Unchanged on exit.   
    N      - INTEGER.   
             On entry,  N specifies the order of the matrix C.  N must be   
             at least zero.   
             Unchanged on exit.   
    K      - INTEGER.   
             On entry with  TRANS = 'N' or 'n',  K  specifies  the number   
             of  columns  of the  matrices  A and B,  and on  entry  with   
             TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number   
             of rows of the matrices  A and B.  K must be at least  zero.   
             Unchanged on exit.   
    ALPHA  - REAL            .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   
    A      - REAL             array of DIMENSION ( LDA, ka ), where ka is   
             k  when  TRANS = 'N' or 'n',  and is  n  otherwise.   
             Before entry with  TRANS = 'N' or 'n',  the  leading  n by k   
             part of the array  A  must contain the matrix  A,  otherwise   
             the leading  k by n  part of the array  A  must contain  the   
             matrix A.   
             Unchanged on exit.   
    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared   
             in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'   
             then  LDA must be at least  f2cmax( 1, n ), otherwise  LDA must   
             be at least  f2cmax( 1, k ).   
             Unchanged on exit.   
    B      - REAL             array of DIMENSION ( LDB, kb ), where kb is   
             k  when  TRANS = 'N' or 'n',  and is  n  otherwise.   
             Before entry with  TRANS = 'N' or 'n',  the  leading  n by k   
             part of the array  B  must contain the matrix  B,  otherwise   
             the leading  k by n  part of the array  B  must contain  the   
             matrix B.   
             Unchanged on exit.   
    LDB    - INTEGER.   
             On entry, LDB specifies the first dimension of B as declared   
             in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'   
             then  LDB must be at least  f2cmax( 1, n ), otherwise  LDB must   
             be at least  f2cmax( 1, k ).   
             Unchanged on exit.   
    BETA   - REAL            .   
             On entry, BETA specifies the scalar beta.   
             Unchanged on exit.   
    C      - REAL             array of DIMENSION ( LDC, n ).   
             Before entry  with  UPLO = 'U' or 'u',  the leading  n by n   
             upper triangular part of the array C must contain the upper   
             triangular part  of the  symmetric matrix  and the strictly   
             lower triangular part of C is not referenced.  On exit, the   
             upper triangular part of the array  C is overwritten by the   
             upper triangular part of the updated matrix.   
             Before entry  with  UPLO = 'L' or 'l',  the leading  n by n   
             lower triangular part of the array C must contain the lower   
             triangular part  of the  symmetric matrix  and the strictly   
             upper triangular part of C is not referenced.  On exit, the   
             lower triangular part of the array  C is overwritten by the   
             lower triangular part of the updated matrix.   
    LDC    - INTEGER.   
             On entry, LDC specifies the first dimension of C as declared   
             in  the  calling  (sub)  program.   LDC  must  be  at  least   
             f2cmax( 1, n ).   
             Unchanged on exit.   
    Level 3 Blas routine.   
    -- Written on 8-February-1989.   
       Jack Dongarra, Argonne National Laboratory.   
       Iain Duff, AERE Harwell.   
       Jeremy Du Croz, Numerical Algorithms Group Ltd.   
       Sven Hammarling, Numerical Algorithms Group Ltd.   
       Test the input parameters.   
       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;
    /* Function Body */
    if (lsame_(trans, "N")) {
	nrowa = *n;
    } else {
	nrowa = *k;
    }
    upper = lsame_(uplo, "U");
    info = 0;
    if (! upper && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (! lsame_(trans, "N") && ! lsame_(trans, 
	    "T") && ! lsame_(trans, "C")) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*k < 0) {
	info = 4;
    } else if (*lda < f2cmax(1,nrowa)) {
	info = 7;
    } else if (*ldb < f2cmax(1,nrowa)) {
	info = 9;
    } else if (*ldc < f2cmax(1,*n)) {
	info = 12;
    }
    if (info != 0) {
	xerbla_("SSYR2K", &info);
	return 0;
    }
/*     Quick return if possible. */
    if (*n == 0 || (*alpha == 0.f || *k == 0) && *beta == 1.f) {
	return 0;
    }
/*     And when  alpha.eq.zero. */
    if (*alpha == 0.f) {
	if (upper) {
	    if (*beta == 0.f) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = 0.f;
/* L10: */
		    }
/* L20: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = *beta * c___ref(i__, j);
/* L30: */
		    }
/* L40: */
		}
	    }
	} else {
	    if (*beta == 0.f) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c___ref(i__, j) = 0.f;
/* L50: */
		    }
/* L60: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c___ref(i__, j) = *beta * c___ref(i__, j);
/* L70: */
		    }
/* L80: */
		}
	    }
	}
	return 0;
    }
/*     Start the operations. */
    if (lsame_(trans, "N")) {
/*        Form  C := alpha*A*B' + alpha*B*A' + C. */
	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.f) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = 0.f;
/* L90: */
		    }
		} else if (*beta != 1.f) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = *beta * c___ref(i__, j);
/* L100: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (a_ref(j, l) != 0.f || b_ref(j, l) != 0.f) {
			temp1 = *alpha * b_ref(j, l);
			temp2 = *alpha * a_ref(j, l);
			i__3 = j;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    c___ref(i__, j) = c___ref(i__, j) + a_ref(i__, l) 
				    * temp1 + b_ref(i__, l) * temp2;
/* L110: */
			}
		    }
/* L120: */
		}
/* L130: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.f) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c___ref(i__, j) = 0.f;
/* L140: */
		    }
		} else if (*beta != 1.f) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c___ref(i__, j) = *beta * c___ref(i__, j);
/* L150: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (a_ref(j, l) != 0.f || b_ref(j, l) != 0.f) {
			temp1 = *alpha * b_ref(j, l);
			temp2 = *alpha * a_ref(j, l);
			i__3 = *n;
			for (i__ = j; i__ <= i__3; ++i__) {
			    c___ref(i__, j) = c___ref(i__, j) + a_ref(i__, l) 
				    * temp1 + b_ref(i__, l) * temp2;
/* L160: */
			}
		    }
/* L170: */
		}
/* L180: */
	    }
	}
    } else {
/*        Form  C := alpha*A'*B + alpha*B'*A + C. */
	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp1 = 0.f;
		    temp2 = 0.f;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp1 += a_ref(l, i__) * b_ref(l, j);
			temp2 += b_ref(l, i__) * a_ref(l, j);
/* L190: */
		    }
		    if (*beta == 0.f) {
			c___ref(i__, j) = *alpha * temp1 + *alpha * temp2;
		    } else {
			c___ref(i__, j) = *beta * c___ref(i__, j) + *alpha * 
				temp1 + *alpha * temp2;
		    }
/* L200: */
		}
/* L210: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = j; i__ <= i__2; ++i__) {
		    temp1 = 0.f;
		    temp2 = 0.f;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp1 += a_ref(l, i__) * b_ref(l, j);
			temp2 += b_ref(l, i__) * a_ref(l, j);
/* L220: */
		    }
		    if (*beta == 0.f) {
			c___ref(i__, j) = *alpha * temp1 + *alpha * temp2;
		    } else {
			c___ref(i__, j) = *beta * c___ref(i__, j) + *alpha * 
				temp1 + *alpha * temp2;
		    }
/* L230: */
		}
/* L240: */
	    }
	}
    }
    return 0;
/*     End of SSYR2K. */
} /* ssyr2k_ */
#undef c___ref
#undef b_ref
#undef a_ref




/* Subroutine */ int ssytd2_(char *uplo, integer *n, real *a, integer *lda, 
	real *d__, real *e, real *tau, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SSYTD2 reduces a real symmetric matrix A to symmetric tridiagonal   
    form T by an orthogonal similarity transformation: Q' * A * Q = T.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is stored:   
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) REAL array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading   
            n-by-n upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading n-by-n lower triangular part of A contains the lower   
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   
            On exit, if UPLO = 'U', the diagonal and first superdiagonal   
            of A are overwritten by the corresponding elements of the   
            tridiagonal matrix T, and the elements above the first   
            superdiagonal, with the array TAU, represent the orthogonal   
            matrix Q as a product of elementary reflectors; if UPLO   
            = 'L', the diagonal and first subdiagonal of A are over-   
            written by the corresponding elements of the tridiagonal   
            matrix T, and the elements below the first subdiagonal, with   
            the array TAU, represent the orthogonal matrix Q as a product   
            of elementary reflectors. See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= f2cmax(1,N).   

    D       (output) REAL array, dimension (N)   
            The diagonal elements of the tridiagonal matrix T:   
            D(i) = A(i,i).   

    E       (output) REAL array, dimension (N-1)   
            The off-diagonal elements of the tridiagonal matrix T:   
            E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.   

    TAU     (output) REAL array, dimension (N-1)   
            The scalar factors of the elementary reflectors (see Further   
            Details).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    If UPLO = 'U', the matrix Q is represented as a product of elementary   
    reflectors   

       Q = H(n-1) . . . H(2) H(1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in   
    A(1:i-1,i+1), and tau in TAU(i).   

    If UPLO = 'L', the matrix Q is represented as a product of elementary   
    reflectors   

       Q = H(1) H(2) . . . H(n-1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),   
    and tau in TAU(i).   

    The contents of A on exit are illustrated by the following examples   
    with n = 5:   

    if UPLO = 'U':                       if UPLO = 'L':   

      (  d   e   v2  v3  v4 )              (  d                  )   
      (      d   e   v3  v4 )              (  e   d              )   
      (          d   e   v4 )              (  v1  e   d          )   
      (              d   e  )              (  v1  v2  e   d      )   
      (                  d  )              (  v1  v2  v3  e   d  )   

    where d and e denote diagonal and off-diagonal elements of T, and vi   
    denotes an element of the vector defining H(i).   

    =====================================================================   


       Test the input parameters   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    static real c_b8 = 0.f;
    static real c_b14 = -1.f;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    static real taui;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static integer i__;
    extern /* Subroutine */ int ssyr2_(char *, integer *, real *, real *, 
	    integer *, real *, integer *, real *, integer *);
    static real alpha;
    extern logical lsame_(char *, char *);
    static logical upper;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *), ssymv_(char *, integer *, real *, real *, 
	    integer *, real *, integer *, real *, real *, integer *), 
	    xerbla_(char *, integer *), slarfg_(integer *, real *, 
	    real *, integer *, real *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --d__;
    --e;
    --tau;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < f2cmax(1,*n)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SSYTD2", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	return 0;
    }

    if (upper) {

/*        Reduce the upper triangle of A */

	for (i__ = *n - 1; i__ >= 1; --i__) {

/*           Generate elementary reflector H(i) = I - tau * v * v'   
             to annihilate A(1:i-1,i+1) */

	    slarfg_(&i__, &a_ref(i__, i__ + 1), &a_ref(1, i__ + 1), &c__1, &
		    taui);
	    e[i__] = a_ref(i__, i__ + 1);

	    if (taui != 0.f) {

/*              Apply H(i) from both sides to A(1:i,1:i) */

		a_ref(i__, i__ + 1) = 1.f;

/*              Compute  x := tau * A * v  storing x in TAU(1:i) */

		ssymv_(uplo, &i__, &taui, &a[a_offset], lda, &a_ref(1, i__ + 
			1), &c__1, &c_b8, &tau[1], &c__1);

/*              Compute  w := x - 1/2 * tau * (x'*v) * v */

		alpha = taui * -.5f * sdot_(&i__, &tau[1], &c__1, &a_ref(1, 
			i__ + 1), &c__1);
		saxpy_(&i__, &alpha, &a_ref(1, i__ + 1), &c__1, &tau[1], &
			c__1);

/*              Apply the transformation as a rank-2 update:   
                   A := A - v * w' - w * v' */

		ssyr2_(uplo, &i__, &c_b14, &a_ref(1, i__ + 1), &c__1, &tau[1],
			 &c__1, &a[a_offset], lda);

		a_ref(i__, i__ + 1) = e[i__];
	    }
	    d__[i__ + 1] = a_ref(i__ + 1, i__ + 1);
	    tau[i__] = taui;
/* L10: */
	}
	d__[1] = a_ref(1, 1);
    } else {

/*        Reduce the lower triangle of A */

	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Generate elementary reflector H(i) = I - tau * v * v'   
             to annihilate A(i+2:n,i)   

   Computing MIN */
	    i__2 = i__ + 2;
	    i__3 = *n - i__;
	    slarfg_(&i__3, &a_ref(i__ + 1, i__), &a_ref(f2cmin(i__2,*n), i__), &
		    c__1, &taui);
	    e[i__] = a_ref(i__ + 1, i__);

	    if (taui != 0.f) {

/*              Apply H(i) from both sides to A(i+1:n,i+1:n) */

		a_ref(i__ + 1, i__) = 1.f;

/*              Compute  x := tau * A * v  storing y in TAU(i:n-1) */

		i__2 = *n - i__;
		ssymv_(uplo, &i__2, &taui, &a_ref(i__ + 1, i__ + 1), lda, &
			a_ref(i__ + 1, i__), &c__1, &c_b8, &tau[i__], &c__1);

/*              Compute  w := x - 1/2 * tau * (x'*v) * v */

		i__2 = *n - i__;
		alpha = taui * -.5f * sdot_(&i__2, &tau[i__], &c__1, &a_ref(
			i__ + 1, i__), &c__1);
		i__2 = *n - i__;
		saxpy_(&i__2, &alpha, &a_ref(i__ + 1, i__), &c__1, &tau[i__], 
			&c__1);

/*              Apply the transformation as a rank-2 update:   
                   A := A - v * w' - w * v' */

		i__2 = *n - i__;
		ssyr2_(uplo, &i__2, &c_b14, &a_ref(i__ + 1, i__), &c__1, &tau[
			i__], &c__1, &a_ref(i__ + 1, i__ + 1), lda)
			;

		a_ref(i__ + 1, i__) = e[i__];
	    }
	    d__[i__] = a_ref(i__, i__);
	    tau[i__] = taui;
/* L20: */
	}
	d__[*n] = a_ref(*n, *n);
    }

    return 0;

/*     End of SSYTD2 */

} /* ssytd2_ */

#undef a_ref





/* Subroutine */ int ssytrd_(char *uplo, integer *n, real *a, integer *lda, 
	real *d__, real *e, real *tau, real *work, integer *lwork, integer *
	info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    SSYTRD reduces a real symmetric matrix A to real symmetric   
    tridiagonal form T by an orthogonal similarity transformation:   
    Q**T * A * Q = T.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) REAL array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading   
            N-by-N upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of A contains the lower   
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   
            On exit, if UPLO = 'U', the diagonal and first superdiagonal   
            of A are overwritten by the corresponding elements of the   
            tridiagonal matrix T, and the elements above the first   
            superdiagonal, with the array TAU, represent the orthogonal   
            matrix Q as a product of elementary reflectors; if UPLO   
            = 'L', the diagonal and first subdiagonal of A are over-   
            written by the corresponding elements of the tridiagonal   
            matrix T, and the elements below the first subdiagonal, with   
            the array TAU, represent the orthogonal matrix Q as a product   
            of elementary reflectors. See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= f2cmax(1,N).   

    D       (output) REAL array, dimension (N)   
            The diagonal elements of the tridiagonal matrix T:   
            D(i) = A(i,i).   

    E       (output) REAL array, dimension (N-1)   
            The off-diagonal elements of the tridiagonal matrix T:   
            E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.   

    TAU     (output) REAL array, dimension (N-1)   
            The scalar factors of the elementary reflectors (see Further   
            Details).   

    WORK    (workspace/output) REAL array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= 1.   
            For optimum performance LWORK >= N*NB, where NB is the   
            optimal blocksize.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    If UPLO = 'U', the matrix Q is represented as a product of elementary   
    reflectors   

       Q = H(n-1) . . . H(2) H(1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in   
    A(1:i-1,i+1), and tau in TAU(i).   

    If UPLO = 'L', the matrix Q is represented as a product of elementary   
    reflectors   

       Q = H(1) H(2) . . . H(n-1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),   
    and tau in TAU(i).   

    The contents of A on exit are illustrated by the following examples   
    with n = 5:   

    if UPLO = 'U':                       if UPLO = 'L':   

      (  d   e   v2  v3  v4 )              (  d                  )   
      (      d   e   v3  v4 )              (  e   d              )   
      (          d   e   v4 )              (  v1  e   d          )   
      (              d   e  )              (  v1  v2  e   d      )   
      (                  d  )              (  v1  v2  v3  e   d  )   

    where d and e denote diagonal and off-diagonal elements of T, and vi   
    denotes an element of the vector defining H(i).   

    =====================================================================   


       Test the input parameters   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__3 = 3;
    static integer c__2 = 2;
    static real c_b22 = -1.f;
    static real c_b23 = 1.f;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *);
    static integer nbmin, iinfo;
    static logical upper;
    static integer nb, kk;
    extern /* Subroutine */ int ssytd2_(char *, integer *, real *, integer *, 
	    real *, real *, real *, integer *), ssyr2k_(char *, char *
	    , integer *, integer *, real *, real *, integer *, real *, 
	    integer *, real *, real *, integer *);
    static integer nx;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int slatrd_(char *, integer *, integer *, real *, 
	    integer *, real *, real *, real *, integer *);
    static integer ldwork, lwkopt;
    static logical lquery;
    static integer iws;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --d__;
    --e;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    lquery = *lwork == -1;
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < f2cmax(1,*n)) {
	*info = -4;
    } else if (*lwork < 1 && ! lquery) {
	*info = -9;
    }

    if (*info == 0) {

/*        Determine the block size. */

	nb = ilaenv_(&c__1, "SSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
	lwkopt = *n * nb;
	work[1] = (real) lwkopt;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SSYTRD", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	work[1] = 1.f;
	return 0;
    }

    nx = *n;
    iws = 1;
    if (nb > 1 && nb < *n) {

/*        Determine when to cross over from blocked to unblocked code   
          (last block is always handled by unblocked code).   

   Computing MAX */
	i__1 = nb, i__2 = ilaenv_(&c__3, "SSYTRD", uplo, n, &c_n1, &c_n1, &
		c_n1, (ftnlen)6, (ftnlen)1);
	nx = f2cmax(i__1,i__2);
	if (nx < *n) {

/*           Determine if workspace is large enough for blocked code. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  determine the   
                minimum value of NB, and reduce NB or force use of   
                unblocked code by setting NX = N.   

   Computing MAX */
		i__1 = *lwork / ldwork;
		nb = f2cmax(i__1,1);
		nbmin = ilaenv_(&c__2, "SSYTRD", uplo, n, &c_n1, &c_n1, &c_n1,
			 (ftnlen)6, (ftnlen)1);
		if (nb < nbmin) {
		    nx = *n;
		}
	    }
	} else {
	    nx = *n;
	}
    } else {
	nb = 1;
    }

    if (upper) {

/*        Reduce the upper triangle of A.   
          Columns 1:kk are handled by the unblocked method. */

	kk = *n - (*n - nx + nb - 1) / nb * nb;
	i__1 = kk + 1;
	i__2 = -nb;
	for (i__ = *n - nb + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
		i__2) {

/*           Reduce columns i:i+nb-1 to tridiagonal form and form the   
             matrix W which is needed to update the unreduced part of   
             the matrix */

	    i__3 = i__ + nb - 1;
	    slatrd_(uplo, &i__3, &nb, &a[a_offset], lda, &e[1], &tau[1], &
		    work[1], &ldwork);

/*           Update the unreduced submatrix A(1:i-1,1:i-1), using an   
             update of the form:  A := A - V*W' - W*V' */

	    i__3 = i__ - 1;
	    ssyr2k_(uplo, "No transpose", &i__3, &nb, &c_b22, &a_ref(1, i__), 
		    lda, &work[1], &ldwork, &c_b23, &a[a_offset], lda);

/*           Copy superdiagonal elements back into A, and diagonal   
             elements into D */

	    i__3 = i__ + nb - 1;
	    for (j = i__; j <= i__3; ++j) {
		a_ref(j - 1, j) = e[j - 1];
		d__[j] = a_ref(j, j);
/* L10: */
	    }
/* L20: */
	}

/*        Use unblocked code to reduce the last or only block */

	ssytd2_(uplo, &kk, &a[a_offset], lda, &d__[1], &e[1], &tau[1], &iinfo);
    } else {

/*        Reduce the lower triangle of A */

	i__2 = *n - nx;
	i__1 = nb;
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {

/*           Reduce columns i:i+nb-1 to tridiagonal form and form the   
             matrix W which is needed to update the unreduced part of   
             the matrix */

	    i__3 = *n - i__ + 1;
	    slatrd_(uplo, &i__3, &nb, &a_ref(i__, i__), lda, &e[i__], &tau[
		    i__], &work[1], &ldwork);

/*           Update the unreduced submatrix A(i+ib:n,i+ib:n), using   
             an update of the form:  A := A - V*W' - W*V' */

	    i__3 = *n - i__ - nb + 1;
	    ssyr2k_(uplo, "No transpose", &i__3, &nb, &c_b22, &a_ref(i__ + nb,
		     i__), lda, &work[nb + 1], &ldwork, &c_b23, &a_ref(i__ + 
		    nb, i__ + nb), lda);

/*           Copy subdiagonal elements back into A, and diagonal   
             elements into D */

	    i__3 = i__ + nb - 1;
	    for (j = i__; j <= i__3; ++j) {
		a_ref(j + 1, j) = e[j];
		d__[j] = a_ref(j, j);
/* L30: */
	    }
/* L40: */
	}

/*        Use unblocked code to reduce the last or only block */

	i__1 = *n - i__ + 1;
	ssytd2_(uplo, &i__1, &a_ref(i__, i__), lda, &d__[i__], &e[i__], &tau[
		i__], &iinfo);
    }

    work[1] = (real) lwkopt;
    return 0;

/*     End of SSYTRD */

} /* ssytrd_ */

#undef a_ref





/* Subroutine */ int strmm_(char *side, char *uplo, char *transa, char *diag, 
	integer *m, integer *n, real *alpha, real *a, integer *lda, real *b, 
	integer *ldb)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    /* Local variables */
    static integer info;
    static real temp;
    static integer i__, j, k;
    static logical lside;
    extern logical lsame_(char *, char *);
    static integer nrowa;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical nounit;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
/*  Purpose   
    =======   
    STRMM  performs one of the matrix-matrix operations   
       B := alpha*op( A )*B,   or   B := alpha*B*op( A ),   
    where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or   
    non-unit,  upper or lower triangular matrix  and  op( A )  is one  of   
       op( A ) = A   or   op( A ) = A'.   
    Parameters   
    ==========   
    SIDE   - CHARACTER*1.   
             On entry,  SIDE specifies whether  op( A ) multiplies B from   
             the left or right as follows:   
                SIDE = 'L' or 'l'   B := alpha*op( A )*B.   
                SIDE = 'R' or 'r'   B := alpha*B*op( A ).   
             Unchanged on exit.   
    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix A is an upper or   
             lower triangular matrix as follows:   
                UPLO = 'U' or 'u'   A is an upper triangular matrix.   
                UPLO = 'L' or 'l'   A is a lower triangular matrix.   
             Unchanged on exit.   
    TRANSA - CHARACTER*1.   
             On entry, TRANSA specifies the form of op( A ) to be used in   
             the matrix multiplication as follows:   
                TRANSA = 'N' or 'n'   op( A ) = A.   
                TRANSA = 'T' or 't'   op( A ) = A'.   
                TRANSA = 'C' or 'c'   op( A ) = A'.   
             Unchanged on exit.   
    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit triangular   
             as follows:   
                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   
                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   
             Unchanged on exit.   
    M      - INTEGER.   
             On entry, M specifies the number of rows of B. M must be at   
             least zero.   
             Unchanged on exit.   
    N      - INTEGER.   
             On entry, N specifies the number of columns of B.  N must be   
             at least zero.   
             Unchanged on exit.   
    ALPHA  - REAL            .   
             On entry,  ALPHA specifies the scalar  alpha. When  alpha is   
             zero then  A is not referenced and  B need not be set before   
             entry.   
             Unchanged on exit.   
    A      - REAL             array of DIMENSION ( LDA, k ), where k is m   
             when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.   
             Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k   
             upper triangular part of the array  A must contain the upper   
             triangular matrix  and the strictly lower triangular part of   
             A is not referenced.   
             Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k   
             lower triangular part of the array  A must contain the lower   
             triangular matrix  and the strictly upper triangular part of   
             A is not referenced.   
             Note that when  DIAG = 'U' or 'u',  the diagonal elements of   
             A  are not referenced either,  but are assumed to be  unity.   
             Unchanged on exit.   
    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared   
             in the calling (sub) program.  When  SIDE = 'L' or 'l'  then   
             LDA  must be at least  f2cmax( 1, m ),  when  SIDE = 'R' or 'r'   
             then LDA must be at least f2cmax( 1, n ).   
             Unchanged on exit.   
    B      - REAL             array of DIMENSION ( LDB, n ).   
             Before entry,  the leading  m by n part of the array  B must   
             contain the matrix  B,  and  on exit  is overwritten  by the   
             transformed matrix.   
    LDB    - INTEGER.   
             On entry, LDB specifies the first dimension of B as declared   
             in  the  calling  (sub)  program.   LDB  must  be  at  least   
             f2cmax( 1, m ).   
             Unchanged on exit.   
    Level 3 Blas routine.   
    -- Written on 8-February-1989.   
       Jack Dongarra, Argonne National Laboratory.   
       Iain Duff, AERE Harwell.   
       Jeremy Du Croz, Numerical Algorithms Group Ltd.   
       Sven Hammarling, Numerical Algorithms Group Ltd.   
       Test the input parameters.   
       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    /* Function Body */
    lside = lsame_(side, "L");
    if (lside) {
	nrowa = *m;
    } else {
	nrowa = *n;
    }
    nounit = lsame_(diag, "N");
    upper = lsame_(uplo, "U");
    info = 0;
    if (! lside && ! lsame_(side, "R")) {
	info = 1;
    } else if (! upper && ! lsame_(uplo, "L")) {
	info = 2;
    } else if (! lsame_(transa, "N") && ! lsame_(transa,
	     "T") && ! lsame_(transa, "C")) {
	info = 3;
    } else if (! lsame_(diag, "U") && ! lsame_(diag, 
	    "N")) {
	info = 4;
    } else if (*m < 0) {
	info = 5;
    } else if (*n < 0) {
	info = 6;
    } else if (*lda < f2cmax(1,nrowa)) {
	info = 9;
    } else if (*ldb < f2cmax(1,*m)) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("STRMM ", &info);
	return 0;
    }
/*     Quick return if possible. */
    if (*n == 0) {
	return 0;
    }
/*     And when  alpha.eq.zero. */
    if (*alpha == 0.f) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		b_ref(i__, j) = 0.f;
/* L10: */
	    }
/* L20: */
	}
	return 0;
    }
/*     Start the operations. */
    if (lside) {
	if (lsame_(transa, "N")) {
/*           Form  B := alpha*A*B. */
	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (k = 1; k <= i__2; ++k) {
			if (b_ref(k, j) != 0.f) {
			    temp = *alpha * b_ref(k, j);
			    i__3 = k - 1;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b_ref(i__, j) = b_ref(i__, j) + temp * a_ref(
					i__, k);
/* L30: */
			    }
			    if (nounit) {
				temp *= a_ref(k, k);
			    }
			    b_ref(k, j) = temp;
			}
/* L40: */
		    }
/* L50: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    for (k = *m; k >= 1; --k) {
			if (b_ref(k, j) != 0.f) {
			    temp = *alpha * b_ref(k, j);
			    b_ref(k, j) = temp;
			    if (nounit) {
				b_ref(k, j) = b_ref(k, j) * a_ref(k, k);
			    }
			    i__2 = *m;
			    for (i__ = k + 1; i__ <= i__2; ++i__) {
				b_ref(i__, j) = b_ref(i__, j) + temp * a_ref(
					i__, k);
/* L60: */
			    }
			}
/* L70: */
		    }
/* L80: */
		}
	    }
	} else {
/*           Form  B := alpha*A'*B. */
	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    for (i__ = *m; i__ >= 1; --i__) {
			temp = b_ref(i__, j);
			if (nounit) {
			    temp *= a_ref(i__, i__);
			}
			i__2 = i__ - 1;
			for (k = 1; k <= i__2; ++k) {
			    temp += a_ref(k, i__) * b_ref(k, j);
/* L90: */
			}
			b_ref(i__, j) = *alpha * temp;
/* L100: */
		    }
/* L110: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp = b_ref(i__, j);
			if (nounit) {
			    temp *= a_ref(i__, i__);
			}
			i__3 = *m;
			for (k = i__ + 1; k <= i__3; ++k) {
			    temp += a_ref(k, i__) * b_ref(k, j);
/* L120: */
			}
			b_ref(i__, j) = *alpha * temp;
/* L130: */
		    }
/* L140: */
		}
	    }
	}
    } else {
	if (lsame_(transa, "N")) {
/*           Form  B := alpha*B*A. */
	    if (upper) {
		for (j = *n; j >= 1; --j) {
		    temp = *alpha;
		    if (nounit) {
			temp *= a_ref(j, j);
		    }
		    i__1 = *m;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			b_ref(i__, j) = temp * b_ref(i__, j);
/* L150: */
		    }
		    i__1 = j - 1;
		    for (k = 1; k <= i__1; ++k) {
			if (a_ref(k, j) != 0.f) {
			    temp = *alpha * a_ref(k, j);
			    i__2 = *m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b_ref(i__, j) = b_ref(i__, j) + temp * b_ref(
					i__, k);
/* L160: */
			    }
			}
/* L170: */
		    }
/* L180: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = *alpha;
		    if (nounit) {
			temp *= a_ref(j, j);
		    }
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			b_ref(i__, j) = temp * b_ref(i__, j);
/* L190: */
		    }
		    i__2 = *n;
		    for (k = j + 1; k <= i__2; ++k) {
			if (a_ref(k, j) != 0.f) {
			    temp = *alpha * a_ref(k, j);
			    i__3 = *m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b_ref(i__, j) = b_ref(i__, j) + temp * b_ref(
					i__, k);
/* L200: */
			    }
			}
/* L210: */
		    }
/* L220: */
		}
	    }
	} else {
/*           Form  B := alpha*B*A'. */
	    if (upper) {
		i__1 = *n;
		for (k = 1; k <= i__1; ++k) {
		    i__2 = k - 1;
		    for (j = 1; j <= i__2; ++j) {
			if (a_ref(j, k) != 0.f) {
			    temp = *alpha * a_ref(j, k);
			    i__3 = *m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b_ref(i__, j) = b_ref(i__, j) + temp * b_ref(
					i__, k);
/* L230: */
			    }
			}
/* L240: */
		    }
		    temp = *alpha;
		    if (nounit) {
			temp *= a_ref(k, k);
		    }
		    if (temp != 1.f) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b_ref(i__, k) = temp * b_ref(i__, k);
/* L250: */
			}
		    }
/* L260: */
		}
	    } else {
		for (k = *n; k >= 1; --k) {
		    i__1 = *n;
		    for (j = k + 1; j <= i__1; ++j) {
			if (a_ref(j, k) != 0.f) {
			    temp = *alpha * a_ref(j, k);
			    i__2 = *m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b_ref(i__, j) = b_ref(i__, j) + temp * b_ref(
					i__, k);
/* L270: */
			    }
			}
/* L280: */
		    }
		    temp = *alpha;
		    if (nounit) {
			temp *= a_ref(k, k);
		    }
		    if (temp != 1.f) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b_ref(i__, k) = temp * b_ref(i__, k);
/* L290: */
			}
		    }
/* L300: */
		}
	    }
	}
    }
    return 0;
/*     End of STRMM . */
} /* strmm_ */
#undef b_ref
#undef a_ref




/* Subroutine */ int strmv_(char *uplo, char *trans, char *diag, integer *n, 
	real *a, integer *lda, real *x, integer *incx)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    static integer info;
    static real temp;
    static integer i__, j;
    extern logical lsame_(char *, char *);
    static integer ix, jx, kx;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical nounit;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
/*  Purpose   
    =======   
    STRMV  performs one of the matrix-vector operations   
       x := A*x,   or   x := A'*x,   
    where x is an n element vector and  A is an n by n unit, or non-unit,   
    upper or lower triangular matrix.   
    Parameters   
    ==========   
    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix is an upper or   
             lower triangular matrix as follows:   
                UPLO = 'U' or 'u'   A is an upper triangular matrix.   
                UPLO = 'L' or 'l'   A is a lower triangular matrix.   
             Unchanged on exit.   
    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the operation to be performed as   
             follows:   
                TRANS = 'N' or 'n'   x := A*x.   
                TRANS = 'T' or 't'   x := A'*x.   
                TRANS = 'C' or 'c'   x := A'*x.   
             Unchanged on exit.   
    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit   
             triangular as follows:   
                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   
                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   
             Unchanged on exit.   
    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   
    A      - REAL             array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper   
             triangular matrix and the strictly lower triangular part of   
             A is not referenced.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower   
             triangular matrix and the strictly upper triangular part of   
             A is not referenced.   
             Note that when  DIAG = 'U' or 'u', the diagonal elements of   
             A are not referenced either, but are assumed to be unity.   
             Unchanged on exit.   
    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared   
             in the calling (sub) program. LDA must be at least   
             f2cmax( 1, n ).   
             Unchanged on exit.   
    X      - REAL             array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element vector x. On exit, X is overwritten with the   
             tranformed vector x.   
    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   
    Level 2 Blas routine.   
    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   
       Test the input parameters.   
       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --x;
    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (! lsame_(trans, "N") && ! lsame_(trans, 
	    "T") && ! lsame_(trans, "C")) {
	info = 2;
    } else if (! lsame_(diag, "U") && ! lsame_(diag, 
	    "N")) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < f2cmax(1,*n)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    }
    if (info != 0) {
	xerbla_("STRMV ", &info);
	return 0;
    }
/*     Quick return if possible. */
    if (*n == 0) {
	return 0;
    }
    nounit = lsame_(diag, "N");
/*     Set up the start point in X if the increment is not unity. This   
       will be  ( N - 1 )*INCX  too small for descending loops. */
    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }
/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A. */
    if (lsame_(trans, "N")) {
/*        Form  x := A*x. */
	if (lsame_(uplo, "U")) {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[j] != 0.f) {
			temp = x[j];
			i__2 = j - 1;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    x[i__] += temp * a_ref(i__, j);
/* L10: */
			}
			if (nounit) {
			    x[j] *= a_ref(j, j);
			}
		    }
/* L20: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[jx] != 0.f) {
			temp = x[jx];
			ix = kx;
			i__2 = j - 1;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    x[ix] += temp * a_ref(i__, j);
			    ix += *incx;
/* L30: */
			}
			if (nounit) {
			    x[jx] *= a_ref(j, j);
			}
		    }
		    jx += *incx;
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (x[j] != 0.f) {
			temp = x[j];
			i__1 = j + 1;
			for (i__ = *n; i__ >= i__1; --i__) {
			    x[i__] += temp * a_ref(i__, j);
/* L50: */
			}
			if (nounit) {
			    x[j] *= a_ref(j, j);
			}
		    }
/* L60: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    if (x[jx] != 0.f) {
			temp = x[jx];
			ix = kx;
			i__1 = j + 1;
			for (i__ = *n; i__ >= i__1; --i__) {
			    x[ix] += temp * a_ref(i__, j);
			    ix -= *incx;
/* L70: */
			}
			if (nounit) {
			    x[jx] *= a_ref(j, j);
			}
		    }
		    jx -= *incx;
/* L80: */
		}
	    }
	}
    } else {
/*        Form  x := A'*x. */
	if (lsame_(uplo, "U")) {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = x[j];
		    if (nounit) {
			temp *= a_ref(j, j);
		    }
		    for (i__ = j - 1; i__ >= 1; --i__) {
			temp += a_ref(i__, j) * x[i__];
/* L90: */
		    }
		    x[j] = temp;
/* L100: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    temp = x[jx];
		    ix = jx;
		    if (nounit) {
			temp *= a_ref(j, j);
		    }
		    for (i__ = j - 1; i__ >= 1; --i__) {
			ix -= *incx;
			temp += a_ref(i__, j) * x[ix];
/* L110: */
		    }
		    x[jx] = temp;
		    jx -= *incx;
/* L120: */
		}
	    }
	} else {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[j];
		    if (nounit) {
			temp *= a_ref(j, j);
		    }
		    i__2 = *n;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			temp += a_ref(i__, j) * x[i__];
/* L130: */
		    }
		    x[j] = temp;
/* L140: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[jx];
		    ix = jx;
		    if (nounit) {
			temp *= a_ref(j, j);
		    }
		    i__2 = *n;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			ix += *incx;
			temp += a_ref(i__, j) * x[ix];
/* L150: */
		    }
		    x[jx] = temp;
		    jx += *incx;
/* L160: */
		}
	    }
	}
    }
    return 0;
/*     End of STRMV . */
} /* strmv_ */
#undef a_ref




/* Subroutine */ int xerbla_(char *srname, integer *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    XERBLA  is an error handler for the LAPACK routines.   
    It is called by an LAPACK routine if an input parameter has an   
    invalid value.  A message is printed and execution stops.   

    Installers may consider modifying the STOP statement in order to   
    call system-specific exception-handling facilities.   

    Arguments   
    =========   

    SRNAME  (input) CHARACTER*6   
            The name of the routine which called XERBLA.   

    INFO    (input) INTEGER   
            The position of the invalid parameter in the parameter list   

            of the calling routine.   

   ===================================================================== 
*/

    printf("** On entry to %6s, parameter number %2i had an illegal value\n",
		srname, *info);

/*     End of XERBLA */

    return 0;
} /* xerbla_ */

