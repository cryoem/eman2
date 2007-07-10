#ifndef EMAN2_LIBEM_VARIMAX_H
#define EMAN2_LIBEM_VARIMAX_H

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



enum rotationScratch
{
	GFNORM = 0,
	GKNORM,
	GLABELS,
	NTRASH
};

enum rotationMethods
{
	IVARIMAX = 0,
	IQUARTIMAX,
	IEQUIMAX,	
	IORTHOMAX,
	IOBLIMIN,
	NMETHODS,
	NOTAVAILABLE = IOBLIMIN /* first non-implemented recognized method*/
};

typedef struct methodName
{
	const char      *name;
	short      length; /*minimum number of letters for recognition*/
} methodName;

int varmx(float *aload,int nv, int nf, int method, float *params,
				  float *fnorm,
				  int itmax, float eps, int verbose);

#endif

