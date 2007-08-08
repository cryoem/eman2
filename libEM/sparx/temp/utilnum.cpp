#include <math.h>
#include "utilnum.h"

int asta2(float *img, int nx, int ny, int ri, double *abaloc, int *klploc)
{
    // calculate the sum of the background pixels, and then set the intensity
    // of these pixels to zero. Also count the number of background pixels.
    // A background pixel is a pixel outside of the circle with radius ri
    // This is done for images assigned to a processor
    //
    int xcent = (nx / 2) + 1;
    int ycent = (ny / 2) + 1;
    int r_squared = ri*ri;

    int x_summand, y_summand;

    for ( int i = 0 ; i < nx ; ++i ) {
	x_summand = (i-xcent) * (i-xcent);
	for ( int j = 0 ; j < ny ; ++j ) {
	    y_summand = (j-ycent) * (j-ycent);
	    if ( x_summand + y_summand > r_squared ) {
		*abaloc += (double) img[j*nx + i];
		//chao set the background to zero
		img[j*nx+i]=0.0;
		++*klploc;
	    }
	}
    }

    return 0;
}

int ifix(float a)
{
    int ia = 0;
    if (a >= 0.0) {
       ia = (int)floor(a);
    }
    else {
       ia = (int)ceil(a);
    }
    return ia;
}


