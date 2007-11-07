// this file is taken from MacAnova by Wei Zhang 1/4/2007. MacAnova is
// also a GPL licensed software, so it is OK to use its source code here
// if things change in the future, we will need re-implement this part


/*
*(C)* The file is part of the source distribution of MacAnova
*(C)* version 4.12 or later
*(C)*
*(C)* Copyright (c) 2001 by Gary Oehlert and Christopher Bingham
*(C)* unless indicated otherwise
*(C)*
*(C)* You may give out copies of this software;  for conditions see the
*(C)* file COPYING included with this distribution
*(C)*
*(C)* This file is distributed WITHOUT ANY WARANTEE; without even
*(C)* the implied warantee of MERCHANTABILITY or FITNESS FOR
*(C)* A PARTICULAR PURPOSE
*(C)*
*/


#ifdef SEGMENTED
#if defined(MPW) || defined(MW_CW)
#pragma segment Rotfac
#endif /*MPW||MW_CW*/
#endif /*SEGMENTED*/

#include <math.h>
#include "varimax.h"

#ifndef ROTMAXITER
#define ROTMAXITER  100
#endif /*ROTMAXITER*/

#ifndef ROTEPSILON
#define ROTEPSILON 1e-5
#endif /*ROTEPSILON*/

/*
  010614 reorganized code and added capability for quartimax and
         orthomax rotation
*/
methodName     Methods[NMETHODS] =
{
	{"Varimax",    3},
	{"Quartimax",  5},
	{"Equimax",    4},
	{"Orthomax",   5},
	{"Oblimin",    7}
};

/*
  Code acquired from Doug Hawkins, University of Minnesota, 12/93

  Comments from Fortran original
	Routine to do a varimax rotation on the real array aload(nv,nnf).
	If nnf is positive, the routine feels free to reverse and reorder
	the factors;this is suppressed in nnf is entered as the negative
	of its actual value.  This suppression is desirable when doing
	a q-mode analysis.

  Translated to C by C. Bingham with the following modifications
     Argument nnf is now nf and is assumed positive and additional argument
     fnorm replaces the local variable fnorm.  If fnorm == (float *) 0,
	 reversing and reordering the factors is suppressed

	 Also argument itmax has been added to control the maximum number of
	 iterations and argument eps provides a convergence limit (originally
	 hard wired as 1e-4).
	 
	 Information on the iteration is printed if verbose != 0

	 varmx() returns a non-zero value if and only if fewer than itmax
	 iterations are required.

	 960919 Modified code to make it easier to add new rotation methods.
            except for slgihtly different error messages, it should have no
            effect on what it does.
     980303 Changed check before computation of rotaton angle to avoid
            atan2 domain error

     010614 added argument lambda to implement orthomax; lambda = 1 is
            varimax; lambda = 0 is quartimax
            Also computation of the criterion moved to separate function
            compcrit and names of certain variables were changed
*/

float compcrit(float *loadings, int nv, int nf, float lambda)
{
	float       crit = 0.0, *loadingj = loadings;
	float       fnv = (float) nv;
	int         i, j;
	
	for (j = 0; j < nf ;j++)
	{
		float       s2 = 0.0;

		for (i = 0;i < nv ;i++)
		{
			float        sq = loadingj[i]*loadingj[i];

			s2 += sq;
			crit += sq*sq;
		}
		crit -= lambda*s2*s2/fnv;
		loadingj += nv;
	} /*for (j = 0;j < nf ;j++)*/

	return (crit);
	
} /*compcrit()*/

/*
  010614 added arguments method and param and modified code so that it
         finds optimal orthomax rotation with parameter lambda = params[0]
         lambda == 1 <==> variamx
         lambda == 2 <==> quartimax
*/
int varmx(float *aload,int nv, int nf, int method, float *params,
				  float *fnorm,
				  int itmax, float eps, int )
/*float aload[nv][1];*/
{
	float         crit, startCrit, fnv = (float) nv;
	float        *aloadj, *aloadk;
	float         denominator, numerator, angl, trot;
	float         eps1 = eps, eps2 = eps;
	float         lambda = 0.0;	//avoid use lambda uninitialized
	int           inoim = 0, ict = 0, irot = 0;
	int           i, j, k, iflip, nf1 = nf - 1;

	if (method <= IORTHOMAX)
	{
		lambda = params[0];
	}

	startCrit = crit = compcrit(aload, nv, nf, lambda);

	do /*while (inoim < 2 && ict < itmax && iflip);*/
	{
		float      oldCrit = crit;

		iflip = 0;
		aloadj = aload;

		for (j = 0;j < nf1 ;j++)
		{
			aloadk = aloadj + nv;
			for (k = j + 1;k < nf ;k++)
			{
				float      a = 0.0, b = 0.0, c = 0.0, d = 0.0, s = 0.0;
				
				for (i = 0;i < nv ;i++)
				{
					float    c2 = aloadj[i]*aloadj[i] - aloadk[i]*aloadk[i];
					float    s2 = 2.0f*aloadj[i]*aloadk[i];

					a += c2;
					b += s2;
					c += c2*c2 - s2*s2;
					d += c2*s2;
				} /*for (i = 0;i < nv ;i++)*/

				denominator = fnv*c + lambda*(b*b - a*a);
				numerator = 2.0f*(fnv*d - lambda*a*b);

				if (fabs(numerator) > eps1*fabs(denominator))
				{
					iflip = 1;
					irot++;
					angl = 0.25f*atan2(numerator,denominator);
					
					c = cos(angl);
					s = sin(angl);
					for (i = 0;i < nv ;i++)
					{
						float   t = c*aloadj[i] + s*aloadk[i];

						aloadk[i] = -s*aloadj[i] + c*aloadk[i];
						aloadj[i] = t;
					} /*for (i = 0;i < nv ;i++)*/
				} /*if (fabs(numerator) >= eps1*fabs(denominator))*/
				aloadk += nv;
			} /*for (k = j + 1;k < nf ;k++)*/
			aloadj += nv;
		} /*for (j = 0;j < nf1 ;j++)*/
		ict++;

		crit = compcrit(aload, nv, nf, lambda);
		trot = (crit > 0.0f) ? (crit - oldCrit)/crit : 0.0f;
		inoim++;
		if (trot > eps2)
		{
			inoim = 0;
		}
	} while (inoim < 2 && ict < itmax && iflip);

	if (fnorm != (float *) 0)
	{
		aloadj = aload;
		for (j = 0;j < nf ;j++)
		{
			float     ssj = 0, sj = 0;
			
			for (i = 0;i < nv ;i++)
			{
				sj += aloadj[i];
				ssj += aloadj[i]*aloadj[i];
			} /*for (i = 0;i < nv ;i++)*/
			fnorm[j] = ssj;
			if (sj <= 0.0)
			{
				for (i = 0;i < nv ;i++)
				{
					aloadj[i] = -aloadj[i];
				}
			} /*if (sj <= 0.0)*/

			aloadk = aload;
			for (k = 0;k < j ;k++)
			{
				if (fnorm[k] < fnorm[j])
				{
					float       t = fnorm[k];

					fnorm[k] = fnorm[j];
					fnorm[j] = t;
					for (i = 0;i < nv ;i++)
					{
						t = aloadj[i];
						aloadj[i] = aloadk[i];
						aloadk[i] = t;
					} /*for (i = 0;i < nv ;i++)*/
				} /*if (fnorm[k] < fnorm[j])*/
				aloadk += nv;
			} /*for (k = 0;k < j ;k++)*/
			aloadj += nv;
		} /*for (j = 0;j < nf ;j++)*/
	} /*if (fnorm != (float *) 0)*/

	/*
	if (verbose)
	{
		char   *outstr = OUTSTR;
		
		outstr += formatChar(outstr, Methods[method].name, CHARASIS);
		outstr += formatChar(outstr, " starting criterion = ", CHARASIS);
		outstr += formatDouble(outstr, startCrit, DODEFAULT | TRIMLEFT);
		outstr += formatChar(outstr, ", final criterion = ", CHARASIS);
		outstr += formatDouble(outstr, crit, DODEFAULT | TRIMLEFT);
		putOUTSTR();
		sprintf(OUTSTR,	"%ld iterations and %ld rotations", ict, irot);
		putOUTSTR();
	} if (verbose)*/

	return (ict < itmax);
} /*varmx()*/

/*
  011125 added argument knorm. 
         knorm != (float *) 0 signals Kaiser normalization with knorm
         providing scratch for row norms
*/
int doRotation(int method, float *aload, int nv, int nf, 
					   float *params, float *fnorm, float *knorm,
					   int itmax, float eps, int verbose)
{
	int        reply = -1;

	if (method <= IORTHOMAX)
	{
		if (method < IORTHOMAX)
		{
			if (method == IVARIMAX)
			{
				params[0] = 1.0;
			}
			else if (method == IQUARTIMAX)
			{
				params[0] = 0.0;
			}
			else
			{
				/* equimax */
				params[0] = .5f*(float) nf;
			}
		}
		if (knorm != (float *) 0)
		{
			int         i, j, k;
			float       s;

			for (i = 0; i < nv; i++)
			{
				k = i;
				s = 0.0;
				for (j = 0; j < nf; j++, k += nv)
				{
					s += aload[k]*aload[k];
				}
				knorm[i] = s = (s > 0) ? sqrt(s) : 1.0f;

				k = i;
				for (j = 0; j < nf; j++, k += nv)
				{
					aload[k] /= s;
				}
			}
		} /*if (fcopy != (float *) 0)*/
		
		reply = varmx(aload, nv, nf, method, params, fnorm, itmax,
					  eps, verbose);
		if (knorm != (float *) 0)
		{
			int      i, j, k = 0;
			
			for (j = 0; j < nf; j++)
			{
				for (i = 0; i < nv; i++, k++)
				{
					aload[k] *= knorm[i];
				}
			}
		}
	} /*if (method <= IORTHOMAX)*/
	return (reply);
} /*doRotation()*/

