/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "ctf.h"
#include "log.h"
#include "emdata.h"
#include "xydata.h"
#include "Assert.h"
#include "projector.h"

#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#include <sys/time.h>

#include <iostream>


using namespace EMAN;

using std::cout;
using std::cin;
using std::endl;
using std::string;

#include "spidutil.h"
#include "spidfft.h"

int  aprings(int  nimg, int nring, float  *imgstk, int nsam, 
             int  nrow, int *numr, float *refcstk, int lcirc, 
             char mode)
{
    int j, status;
    float *wr;

    status = 0;
    wr = (float*) calloc(nring, sizeof(float));
    if (!wr) {
	fprintf(stderr, "aprings: failed to allocate wr!\n");
        status = -1;
        goto EXIT;
    }


    // get wr weights
    ringwe(wr, numr, nring);
    if ( mode == 'H' || mode == 'h' )
	for (j=1;j<=nring;j++) wr(j) = wr(j)/2.0;
    for (j = 1; j<=nimg; j++) {
       apring1(&imgstk(1,1,j), nsam, nrow, &refcstk(1,j),
               lcirc, nring, mode, numr, wr);
    }

 EXIT:
    if (wr) free(wr);
    return status;
}

//-------------------------------------------------------------------

int  apring1(float *sqimg, int nsam , int  nrow, float *imgcirc, int lcirc, 
             int nring   , char mode, int *numr, float *wr)
{
   int    status = 0;
   int    nsb, nse, nrb, nre, maxrin, ltf, lt, lt2, lt3, ns2, nr2;
   int    inr, kcirc, jt, i, j;
   float  fnr2, fns2, yq, fi, dpi, x, y;
   double dfi;

   dpi    = 2.0*atan(1.0); 

   // calculate dimensions for normass
   nsb = -nsam/2;
   nse = -nsb-1+(nsam%2);
   nrb = -nrow/2;
   nre = -nrb-1+(nrow%2);

   //  normalize under the mask,  tried doing this on the
   //  polar rings but it gives different answers. al

   normass(sqimg,nsam,nsb,nse,nrb,nre,numr(1,1),numr(1,nring));

   maxrin = numr(3,nring);
   ns2    = nsam / 2 + 1;
   nr2    = nrow / 2 + 1;
   fnr2   = nr2;
   fns2   = ns2;

   // convert window from image into polar coordinates
        
   if (mode == 'f' || mode == 'F') {
      ltf = 4;
      for (i=1;i<=nring;i++) {
         //  radius of the ring
         inr             = numr(1,i);
         yq              = inr;
         lt              = numr(3,i) / ltf;
         lt2             = lt  + lt;
         lt3             = lt2 + lt;
         dfi             = dpi / lt;
         kcirc           = numr(2,i);
        
         imgcirc(kcirc)     = sqimg(ns2,     nr2+inr);
         imgcirc(lt+kcirc)  = sqimg(ns2+inr, nr2);
         imgcirc(lt2+kcirc) = sqimg(ns2,     nr2-inr);
         imgcirc(lt3+kcirc) = sqimg(ns2-inr, nr2);

         for (j=1;j<=lt - 1;j++) {
            fi = dfi     * j;
            x  = sin(fi) * yq;
            y  = cos(fi) * yq;
            jt = j + kcirc;

            imgcirc(jt)     = quadri(fns2+x,fnr2+y,nsam,nrow,sqimg);
            imgcirc(jt+lt)  = quadri(fns2+y,fnr2-x,nsam,nrow,sqimg);
            imgcirc(jt+lt2) = quadri(fns2-x,fnr2-y,nsam,nrow,sqimg);
            imgcirc(jt+lt3) = quadri(fns2-y,fnr2+x,nsam,nrow,sqimg);
         }
      }
   }
   else if (mode == 'h' || mode == 'H') {
      ltf = 2;
      for (i=1; i<=nring;i++) {
         // radius of the ring
         inr            = numr(1,i);
         yq             = inr;
         lt             = numr(3,i) / ltf;
         dfi            = dpi / lt;
         kcirc          = numr(2,i);

         imgcirc(kcirc)    = sqimg(ns2,     nr2+inr);
         imgcirc(lt+kcirc) = sqimg(ns2+inr, nr2);
 
         for (j=1;j<=lt - 1;j++) {
            fi          = dfi * j;
            x           = sin(fi) * yq;
            y           = cos(fi) * yq;
            jt          = j + kcirc;

            imgcirc(jt)    = quadri(fns2+x,fnr2+y,nsam,nrow,sqimg);
            imgcirc(jt+lt) = quadri(fns2+y,fnr2-x,nsam,nrow,sqimg);
         }
      }
   }

   //  fourier of circ 
   frngs(imgcirc,numr,nring);

   //    weight circ  using wr
   if (wr(1) > 0.0) {
      status = applyws(imgcirc,lcirc,numr,wr,nring);
   }

   return status;
}
    
// -----------------------------------------------------------------
float quadri(float xx, float yy, int nxdata, int nydata, float *fdata)
{
/*
c  purpose: quadratic interpolation 
c 
c  parameters:       xx,yy treated as circularly closed.
c                    fdata - image 1..nxdata, 1..nydata
c
c                    f3    fc       f0, f1, f2, f3 are the values
c                     +             at the grid points.  x is the
c                     + x           point at which the function
c              f2++++f0++++f1       is to be estimated. (it need
c                     +             not be in the first quadrant).
c                     +             fc - the outer corner point
c                    f4             nearest x.
c
c                                   f0 is the value of the fdata at
c                                   fdata(i,j), it is the interior mesh
c                                   point nearest  x.
c                                   the coordinates of f0 are (x0,y0),
c                                   the coordinates of f1 are (xb,y0),
c                                   the coordinates of f2 are (xa,y0),
c                                   the coordinates of f3 are (x0,yb),
c                                   the coordinates of f4 are (x0,ya),
c                                   the coordinates of fc are (xc,yc),
c
c                   o               hxa, hxb are the mesh spacings
c                   +               in the x-direction to the left
c                  hyb              and right of the center point.
c                   +
c            ++hxa++o++hxb++o       hyb, hya are the mesh spacings
c                   +               in the y-direction.
c                  hya
c                   +               hxc equals either  hxb  or  hxa
c                   o               depending on where the corner
c                                   point is located.
c
c                                   construct the interpolant
c                                   f = f0 + c1*(x-x0) +
c                                       c2*(x-x0)*(x-x1) +
c                                       c3*(y-y0) + c4*(y-y0)*(y-y1)
c                                       + c5*(x-x0)*(y-y0)
c
c
*/
    float x, y, dx0, dy0, f0, c1, c2, c3, c4, c5, dxb, dyb;
    float quadri;
    int   i, j, ip1, im1, jp1, jm1, ic, jc, hxc, hyc;

    x = xx;
    y = yy;

    // circular closure
    if (x < 1.0)               x = x+(1 - floor(x) / nxdata) * nxdata;
    if (x > (float)nxdata+0.5) x = fmod(x-1.0,(float)nxdata) + 1.0;
    if (y < 1.0)               y = y+(1 - floor(y) / nydata) * nydata;
    if (y > (float)nydata+0.5) y = fmod(y-1.0,(float)nydata) + 1.0;


    i   = (int) floor(x);
    j   = (int) floor(y);

    dx0 = x - i;
    dy0 = y - j;

    ip1 = i + 1;
    im1 = i - 1;
    jp1 = j + 1;
    jm1 = j - 1;

    if (ip1 > nxdata) ip1 = ip1 - nxdata;
    if (im1 < 1)      im1 = im1 + nxdata;
    if (jp1 > nydata) jp1 = jp1 - nydata;
    if (jm1 < 1)      jm1 = jm1 + nydata;

    f0  = fdata(i,j);
    c1  = fdata(ip1,j) - f0;
    c2  = (c1 - f0 + fdata(im1,j)) * 0.5;
    c3  = fdata(i,jp1) - f0;
    c4  = (c3 - f0 + fdata(i,jm1)) * 0.5;

    dxb = dx0 - 1;
    dyb = dy0 - 1;

    // hxc & hyc are either 1 or -1
    if (dx0 >= 0) {
       hxc = 1;
    }
    else {
       hxc = -1;
    }
    if (dy0 >= 0) {
       hyc = 1;
    }
    else {
       hyc = -1;
    }
 
    ic  = i + hxc;
    jc  = j + hyc;

    if (ic > nxdata) {
       ic = ic - nxdata;
    }
    else if (ic < 1) {
       ic = ic + nxdata;
    }

    if (jc > nydata) {
       jc = jc - nydata;
    }
    else if (jc < 1) {
       jc = jc + nydata;
    }

    c5  =  ( (fdata(ic,jc) - f0 - hxc * c1 - (hxc * (hxc - 1.0)) * c2 
            - hyc * c3 - (hyc * (hyc - 1.0)) * c4) * (hxc * hyc));

    quadri = f0 + dx0 * (c1 + dxb * c2 + dy0 * c5) + dy0 * (c3 + dyb * c4);

    return quadri; 
}

int alprbs(int *numr, int nring, int *lcirc, char mode)
{
/*
c  purpose: appears to circular rings, postitioned
c           in a linear array that holds rings concatenated together.
c           output is dependent on number of rings 
c                                                                      *
c  parameters:   numr(1,i) - ring number                      (sent)
c                numr(2,i) - beginning in circ                (ret.)
c                numr(3,i) - length in circ                   (ret.)
c                nring                                        (sent)
c                lcirc - total length of circ.                (ret.)
c
c image_processing_routine
c
*/
    int i, jp, ip;
    double dpi;
    int status = 0; 
    // hardwire for right now
    int MAXFFT = 32768;

    dpi = pi;
    if (mode == 'f' || mode == 'F') dpi=2*pi;

    *lcirc = 0;
    for (i=1;i<=nring;i++) {
       jp = (int)(dpi*numr(1,i));
       // original fortran code ip = 2**log2(jp), log2(jp) rounds up. 
       ip = (int)( pow(2,ceil(log2(jp))) );
       if (i < nring && jp > ip+ip/2)  ip=min0(MAXFFT,2*ip);

       //  last ring should be oversampled to allow higher accuracy
       //  of peak location (?).
       if (i == nring && jp > ip+ip/5) ip=min0(MAXFFT,2*ip);
       numr(3,i) = ip;
       if (i == 1) {
          numr(2,1) = 1;
       }
       else {
          numr(2,i) = numr(2,i-1)+numr(3,i-1);
       }
       *lcirc = *lcirc + ip;
    }
    return status;
}

// -----------------------------------------------------------------
void ringwe(float *wr, int *numr, int nring)
{
   int i, maxrin;
   float dpi;

   dpi = 8.0*atan(1.0);
   maxrin = numr(3,nring); 
   for  (i=1;i<=nring;i++) {
      wr(i)=(float)(numr(1,i))*dpi/(float)(numr(3,i))
           *(float)(maxrin)/(float)(numr(3,i));
	   cout << numr(1,i)<<"  "<< numr(2,i)<<"  "<< numr(3,i)<<"  " << wr(i)<<"  "<<endl;
   }
}

// -----------------------------------------------------------------
int setnring(int mr, int nr, int iskip)
{
    int nring = 0, i;

    i = mr;
    while (i<=nr) {
       nring = nring+1;
       i = i + iskip;
    } 
    return nring;
}

// -----------------------------------------------------------------
void numrinit(int mr, int nr, int iskip, int *numr)
{
    int nring = 0;
    int i;

    i = mr;
    while (i<=nr) {
      nring++;
      numr(1,nring) = i;
      i=i+iskip;
    }
}

// -----------------------------------------------------------------
void normass(float *sqimg, int nsam, int ns1, int ns2, int nr1, int nr2, 
             int ir1, int ir2)
{
/*
c serially normalizes x by variance
c
c  note   :    for parallel use normas instead al sept 01
c              difficult to add error flag due to use inside
c              parrallel region aug 05 al
*/
   // this is how sqimg is declared in spider
   // dimension  sqimg(ns1:ns2,nr1:nr2)

   double   av,vr,vrinv,dtemp;
   int      i1sq, i2sq, n, ir, jsq, i, j, irow, jcol;

   i1sq = ir1 * ir1;
   i2sq = ir2 * ir2;
   av   = 0.0;
   vr   = 0.0;
   n    = 0;

   for (j=nr1; j<=nr2; j++) {
      jsq = j * j;
      for (i=ns1;i<=ns2;i++) {
         ir = jsq + i * i;
         if (ir >= i1sq && ir <= i2sq) {
	    n  = n  + 1;
            irow = i-ns1+1;
            jcol = j-nr1+1; 
            av = av + sqimg(irow,jcol);
            vr = vr + sqimg(irow,jcol)*sqimg(irow,jcol);
         }
      }
   }

   av   = av / n;
   dtemp = (vr - n * av * av);
   if (dtemp > 0) {
      vr    = sqrt(dtemp / (n-1));
      vrinv = 1.0 / vr;

      //  array operation on x
      for ( j = nr1; j<=nr2;j++) 
         for (i = ns1; i <= ns2; i++) {
            irow = i - ns1 + 1; 
            jcol = j - nr1 + 1; 
            sqimg(irow,jcol) = (sqimg(irow,jcol) - av) * vrinv;
         }
   }
   else {
      // trap for blank image area
      // array operation on x
      for ( j = nr1; j<=nr2;j++) 
         for (i = ns1; i <= ns2; i++) {
            irow = i - ns1 + 1; 
            jcol = j - nr1 + 1; 
            sqimg(irow,jcol) = 0.0;
      }
   }
}

// -----------------------------------------------------------------

void frngs(float *circ, int *numr, int nring)
{
   int i, l; 
 
   for (i=1; i<=nring;i++) {
     l=(int)(log2(numr(3,i)));
     fftr_q(&circ(numr(2,i)),l);
   }
}

// -----------------------------------------------------------------
int  applyws(float *circ, int lcirc, int *numr, float *wr,
             int nring)
{
   int   maxrin, numr3i, numr2i;
   float w;
   int   i, j, jc;
   int   status = 0;

   maxrin = numr(3,nring);
 
   for (i=1;i<=nring;i++) {
      numr3i=numr(3,i);
      numr2i=numr(2,i);
      w=wr(i);
      circ(numr2i)=circ(numr2i)*w;
      if (numr3i == maxrin) {
         circ(numr2i+1)=circ(numr2i+1)*w;
      }
      else { 
         circ(numr2i+1)=circ(numr2i+1)*0.5*w;
      }
      for (j=3;j<=numr3i;j++) {
         jc=j+numr2i-1;
         circ(jc)=circ(jc)*w;
      }
   } 
   if (jc >= lcirc) status = -1;
   return status; 
}

//--------------------------------------------------
void normas(float *sqimg, int nsam, int ns1, int ns2, int nr1, int nr2, 
            int ir1, int ir2)
{
/*
c  purpose:    normalizes ring data.  covered area is: ir1....ir2      *
c                                                                      *
c  parameters:                                                         *
c
c  note   :    i think this is for parallel use only, because normass
c              is quicker for non_parallel use!! al sept 01
c                                                                      *
*/
    //  dimension  sqimg(ns1:ns2,nr1:nr2)

    double     av,vr;
    int        i1sq, i2sq, n, i, j, ir, j2, irow, jcol;

    i1sq = ir1 * ir1;
    i2sq = ir2 * ir2;

    av   = 0.0;
    vr   = 0.0;
    n    = 0;

    for (j=nr1;j<=nr2;j++) {
       j2 = j*j;
       for (i=ns1;i<=ns2;i++) {
          ir = j2 + i*i;
          jcol = j-nr1+1;
          if (ir >= i1sq && ir <= i2sq) {
             n++;
             irow = i-ns1+1;
             av = av + sqimg(irow,jcol);
             vr = vr + sqimg(irow,jcol)*sqimg(irow,jcol);
          }
       }
    } 

    av = av / n;

    //   multiplication is faster
    vr = 1.0 / (sqrt((vr-n*av*av) / (n-1)));

    for (j=nr1; j<=nr2; j++) {
       jcol = j-nr1+1;
       for (i=ns1;i<=ns2;i++) {
          irow = i-ns1+1;
          sqimg(irow,jcol) = (sqimg(irow,jcol) - av ) * vr;
       } 
    }
}

//-------------------------------------------------------
void alrq(float *xim,  int nsam , int nrow , int *numr,
          float *circ, int lcirc, int nring, char mode)
{
/* 
c                                                                     
c  purpose:                                                          
c                                                                   
c  parameters: convert to polar coordinates
c                                                                  
*/
   //  dimension         xim(nsam,nrow),circ(lcirc)
   //  integer           numr(3,nring)

   double dfi, dpi;
   int    ns2, nr2, i, inr, l, nsim, kcirc, lt, j;
   float  yq, xold, yold, fi, x, y;

   ns2 = nsam/2+1;
   nr2 = nrow/2+1;
   dpi = 2.0*atan(1.0);

//#pragma omp   parallel do private(i,j,inr,yq,l,lt,nsim,dfi,kcirc,
//#pragma omp&  xold,yold,fi,x,y)
   for (i=1;i<=nring;i++) {
     // radius of the ring
     inr = numr(1,i);
     yq  = inr;
     l   = numr(3,i);
     if (mode == 'h' || mode == 'H') {
        lt = l/2;
     }
     else if (mode == 'f' || mode == 'F' ) {
        lt = l/4;
     }

     nsim           = lt-1;
     dfi            = dpi/(nsim+1);
     kcirc          = numr(2,i);
     xold           = 0.0;
     yold           = inr;
     circ(kcirc)    = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
     xold           = inr;
     yold           = 0.0;
     circ(lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);

     if (mode == 'f' || mode == 'F') {
        xold              = 0.0;
        yold              = -inr;
        circ(lt+lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
        xold              = -inr;
        yold              = 0.0;
        circ(lt+lt+lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
     }

     for (j=1;j<=nsim;j++) {
        fi               = dfi*j;
        x                = sin(fi)*yq;
        y                = cos(fi)*yq;
        xold             = x;
        yold             = y;
        circ(j+kcirc)    = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
        xold             =  y;
        yold             = -x;
        circ(j+lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);

        if (mode == 'f' || mode == 'F')  {
           xold                = -x;
           yold                = -y;
           circ(j+lt+lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
           xold                = -y;
           yold                =  x;
           circ(j+lt+lt+lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
        };
     }
   }
//#pragma omp   end parallel do 
}

//---------------------------------------------------
int crosrng_ms(float *circ1, float *circ2, int  lcirc, int  nring,
               int   maxrin, int   *numr , double *qn, float *tot,
               double   *qm, double *tmt)
{
/*
c
c  checks both straight & mirrored positions
c
c  input - fourier transforms of rings!!
c  circ1 already multiplied by weights!
c
c  notes: aug 04 attempted speedup using 
c       premultiply  arrays ie( circ12 = circ1 * circ2) much slower
c       various  other attempts  failed to yield improvement
c       this is a very important compute demand in alignmen & refine.
c       optional limit on angular search should be added.
*/

   // dimension         circ1(lcirc),circ2(lcirc)

   // t(maxrin+2), q(maxrin+2), t7(-3:3)
   double *t, *q, t7[7];

   int   ip, jc, numr3i, numr2i, i, j, k, jtot;
   float t1, t2, t3, t4, c1, c2, d1, d2, pos;

   int   status = 0;

   *qn  = 0.0;
   *qm  = 0.0;
   *tot = 0.0;
   *tmt = 0.0; 

   ip = -(int)(log2(maxrin));

   //  c - straight  = circ1 * conjg(circ2)
   //  zero q array
  
   q = (double*)calloc(maxrin+2,sizeof(double));  
   if (!q) {
      status = -1;
      return status;
   }

   //   t - mirrored  = conjg(circ1) * conjg(circ2)
   //   zero t array
   t = (double*)calloc(maxrin+2,sizeof(double));
   if (!t) {
      status = -1;
      return status;
   } 

   //   premultiply  arrays ie( circ12 = circ1 * circ2) much slower

   for (i=1;i<=nring;i++) {

      numr3i = numr(3,i);
      numr2i = numr(2,i);

      t1   = circ1(numr2i) * circ2(numr2i);
      q(1) = q(1)+t1;
      t(1) = t(1)+t1;

      if (numr3i == maxrin)  {
         t1   = circ1(numr2i+1) * circ2(numr2i+1);
         q(2) = q(2)+t1;
         t(2) = t(2)+t1;
      }
      else {
	 t1          = circ1(numr2i+1) * circ2(numr2i+1);
	 q(numr3i+1) = q(numr3i+1)+t1;
      }

      for (j=3;j<=numr3i;j=j+2) {
	 jc     = j+numr2i-1;

 	 c1     = circ1(jc);
 	 c2     = circ1(jc+1);
         d1     = circ2(jc);
         d2     = circ2(jc+1);

  	 t1     = c1 * d1;
 	 t3     = c1 * d2;
 	 t2     = c2 * d2;
 	 t4     = c2 * d1;

	 q(j)   = q(j)   + t1 + t2;
	 q(j+1) = q(j+1) - t3 + t4;
	 t(j)   = t(j)   + t1 - t2;
	 t(j+1) = t(j+1) - t3 - t4;
      } 
  }

  fftr_d(q,ip);

  jtot = 0;
  *qn  = -1.0e20;
  for (j=1; j<=maxrin; j++) {
     if (q(j) >= *qn) {
        *qn  = q(j);
        jtot = j;
     }
  }

  if (jtot <= 0) {
     // some sort of error (probably compiler on mp on altix)
     fprintf(stderr, "no max in crosrng_ms, compiler error\n");
     status = -1;
     return status; 
  }
 
  for (k=-3;k<=3;k++) {
    j = ((jtot+k+maxrin-1)%maxrin)+1;
    t7(k+4) = q(j);
  }

  // this appears to interpolate? al
  prb1d(t7,7,&pos);
  *tot = (float)(jtot)+pos;

  // mirrored
  fftr_d(t,ip);

  // find angle
  *qm = -1.0e20;
  for (j=1; j<=maxrin;j++) {
     if ( t(j) >= *qm ) {
        *qm   = t(j);
        jtot = j;
     }
  }

  // find angle
  for (k=-3;k<=3;k++) {
    j       = ((jtot+k+maxrin-1)%maxrin) + 1;
    t7(k+4) = t(j);
  }

  // this appears to interpolate? al

  prb1d(t7,7,&pos);
  *tmt = float(jtot) + pos;

  free(t);
  free(q);

  return status;
}

//---------------------------------------------------
void prb1d(double *b, int npoint, float *pos)
{
   double  c2,c3;
   int     nhalf;

   nhalf = npoint/2 + 1;
   *pos  = 0.0;

   if (npoint == 7) {
      c2 = 49.*b(1) + 6.*b(2) - 21.*b(3) - 32.*b(4) - 27.*b(5)
         - 6.*b(6) + 31.*b(7);
      c3 = 5.*b(1) - 3.*b(3) - 4.*b(4) - 3.*b(5) + 5.*b(7);
   } 
   else if (npoint == 5) {
      c2 = (74.*b(1) - 23.*b(2) - 60.*b(3) - 37.*b(4)
         + 46.*b(5) ) / (-70.);
      c3 = (2.*b(1) - b(2) - 2.*b(3) - b(4) + 2.*b(5) ) / 14.0;
   }
   else if (npoint == 3) {
      c2 = (5.*b(1) - 8.*b(2) + 3.*b(3) ) / (-2.0);
      c3 = (b(1) - 2.*b(2) + b(3) ) / 2.0;
   }
   else if (npoint == 9) {
      c2 = (1708.*b(1) + 581.*b(2) - 246.*b(3) - 773.*b(4)
         - 1000.*b(5) - 927.*b(6) - 554.*b(7) + 119.*b(8)
         + 1092.*b(9) ) / (-4620.);
      c3 = (28.*b(1) + 7.*b(2) - 8.*b(3) - 17.*b(4) - 20.*b(5)
         - 17.*b(6) - 8.*b(7) + 7.*b(8) + 28.*b(9) ) / 924.0;
   }
   if (c3 != 0.0)  *pos = c2/(2.0*c3) - nhalf;
}
// -----------------------------------------------

void apmaster_1(char mode, float *divas, int nr, int *numth,
                int  lsam, int lrow, int *nsam, int *nrow)
{
/*
c parameters:
c       mode                degree mode                       (input)
c       divas               degrees                           (output)
c       numth               degrees                           (output)
c       lsam                orig size                         (input)
c       lrow                orig size                         (input)
c       nsam                new size                          (output)
c       nrow                new size                          (output)
*/
   int nra;

   if ( mode == 'h') {
      *divas = 180.0;
   }
   else {
      *divas = 360.0;
   }

   *numth = 1;
#ifdef sp_mp
//       find number of omp threads
//        call getthreads(numth)
#endif

   //  calculation of actual dimension of an image to be interpolated
   //  2*(no.of rings)+(0'th element)+2*(margin of 1)

   nra  = ((lsam-1)/2)*2+1;
   if ( ((lrow-1)/2)*2+1 < nra ) nra = ((lrow-1)/2)*2+1;
   if ( 2*nr+3 < nra ) nra = 2*nr+3;

   //  returns circular reduced nsam, nrow
   *nsam = nra;
   *nrow = nra;
}
//----------------------------------------------------------------
void win_resize(float *imgfrom, float *imgto, int lsam, int lrow, 
                int nsam, int nrow, int lr1, int lr2, int ls1, int ls2)
{
/*
c adpated from SPIDER ap_getdat.f
c purpose:       read read windowed image date into array x for 'ap' ops.
c
c parameters:
c       lsam,lrow           image dimensions                  (input)
c       nsam,nrow           output image dimensions           (input)
c       lr1,lr2,ls1,ls2     output image window               (input)
c       imgfrom             input image                       (input)
c       imgto               output output                     (output)
c
*/
   int window;

   int k3, k2, kt;

   //     real, dimension(lsam)                        :: bufin

    window = 0;
    if (lr1 != 1 || ls1 != 1 || lr2 != lrow || ls2 != lsam) window = 1;

    if (window) {
      // window from the whole image
      for (k2=lr1;k2<=lr2;k2++) {
         kt = k2-lr1+1;
         for (k3=ls1;k3<=ls2;k3++) 
            imgto(k3-ls1+1,kt) = imgfrom(k3,k2);
      }
    }
    else {
      // do a copy
      for (k2=1;k2<=lrow;k2++) 
         for (k3=1;k3<=lsam;k3++) 
            imgto(k3,k2) = imgfrom(k3,k2);
    }
}

//--------------------------------------------------
int apmd(EMData *refprj, Dict refparams, EMData *expimg, APMDopt options,
         float *fangles)
{
    // align      expimg using refprj as the reference
    // mr:        first ring 
    // nr:        last ring 
    // iskip:     number of rings skipped between the first and the last
    // mode:      ?
    // refparams: contains the angles used to generate refprj
    // fangles:   output angles assigned to expimg.

    int    mr, nr, iskip;
    char   mode; 
    int    nring, lcirc, nref, nimg, maxrin, status;
    int    nx, ny, nxr, nyr, nsam, nrow;
    int    nwsam, nwrow, ns1, ns2, nr1, nr2;
    int    lq, lr1, lr2, ls1, ls2; 
    float  *refstk, *refcstk, *imgstk, *imgcirc,*imgwindow;
    int    *numr;
    double *totmin, *totmir, *tmt;
    float  *tot;
    float  divas; 
    int    j, k, numth, idi, ldd, nang, mm;
    float  *dlist; 
    double eav;
    float  rang;
    vector <float> angrefs;

    status = 0;

    // retrieve APMD options
    mr = options.mr;
    nr = options.nr;
    iskip = options.iskip;
    mode = options.mode;

    nx    = expimg->get_xsize();
    ny    = expimg->get_ysize();
    nimg  = expimg->get_zsize();

    nxr   = refprj->get_xsize();
    nyr   = refprj->get_ysize();
    nref  = refprj->get_zsize();

    if (nx != nxr || ny != nyr) { 
	status = -2;   
        goto EXIT;
    }  

    nsam = nx;
    nrow = ny;

    // extract image data from the reference image EM object
    refstk = refprj->get_data();
    //  find number of reference-rings
    nring = setnring(mr, nr, iskip);

    numr = (int*) calloc(3*nring,sizeof(int));

    numrinit(mr, nr, iskip, numr);
    // calculates pointers for rings and stores them in numr
    status = alprbs(numr, nring, &lcirc, mode);
    if (status != 0) goto EXIT;

    // convert reference images to rings
    refcstk = (float*)calloc(lcirc*nref,sizeof(float));
    
    // apply ring weights to rings
    status = aprings(nref, nring, refstk, nsam, nrow, numr,
                     refcstk, lcirc, mode);
    if (status != 0) goto EXIT;

    // extract the image data from experimental EM image object
    imgstk = expimg->get_data();

    // set the window size of the exp image: 
    // nwsam,nwrow  is the new size of the image
    apmaster_1(mode, &divas, nr, &numth, nsam, nrow, &nwsam, &nwrow);

    // allocate work space
    imgwindow = (float*)calloc(nwsam*nwrow,sizeof(float));
    imgcirc   = (float*)calloc(lcirc, sizeof(float));

    totmin = (double*)calloc(nref,sizeof(double));
    totmir = (double*)calloc(nref,sizeof(double));
    tot    = (float*)calloc(nref,sizeof(float));
    tmt    = (double*)calloc(nref,sizeof(double));
    if (totmin == NULL || totmir == NULL || tot == NULL || tmt == NULL) {
       fprintf(stderr,"apmd: failed to allocate totmin, totmir, tot, tmt\n");
       status = -1;
       goto EXIT;
    }

    maxrin = numr(3,nring);
    ldd    = 3;
    dlist  = (float*)calloc(ldd*nimg, sizeof(float));
    if ( !dlist ) {
       status = -1;
       goto EXIT;
    }

    // calculate window paramters
    lq =nsam/2+1;
    lr1=(nwrow-1)/2;
    lr2=lq+lr1;
    lr1=lq-lr1;
    lq =nrow/2+1;
    ls1=(nwsam-1)/2;
    ls2=lq+ls1;
    ls1=lq-ls1;

    for (j = 1; j<=nimg; j++) {
       win_resize(&imgstk(1,1,j), imgwindow, nsam, nrow, nwsam, nwrow, 
                  lr1, lr2, ls1, ls2);

       // normalize the image 
       ns1 = -nwsam/2;
       ns2 =  nwsam/2;
       nr1 = -nwrow/2;
       nr2 =  nwrow/2;
       normas(imgwindow, nwsam, ns1, ns2, nr1, nr2, numr(1,1), numr(1,nring));

       // turn the image into rings
       alrq(imgwindow,nwsam,nwrow,numr,imgcirc,lcirc,nring,mode);

       // should we use frng instead??  frng is ||, but || level was moved up PAP
       frngs(imgcirc, numr, nring);

       for (k = 1; k<=nref; k++) {
          status = crosrng_ms(&refcstk(1,k), imgcirc, lcirc     , nring, 
                              maxrin       , numr   , &totmin(k), &tot(k), 
                              &totmir(k)   , &tmt(k));
       }

       eav = 1.0e-20; 
       for (k=1;k<=nref;k++) {
          if ( totmin(k) >= eav && totmin(k) >  totmir(k) ) {
             eav  = totmin(k);
             idi  = k;
             rang = tot(k);
          }
          else if ( totmir(k) >= eav ) {
             eav  = totmir(k);
             idi  = -k;
             rang = tmt(k);
          }
       }

       dlist(1,j) = idi;
       dlist(2,j) = eav;
       rang = (rang-1)/maxrin*divas;
       dlist(3,j) = rang;
       // printf("j = %d, %g %g %g\n", j, dlist(1,j), dlist(2,j), dlist(3,j));
    }

    // now turn dlist into output angles (spider VOMD)
    angrefs = refparams["anglelist"];
    nang    = angrefs.size()/3;
    if (nang != nref) {
       fprintf(stderr, "apmd: nang = %d, nref = %d\n", nang, nref);
       status = -3;
       goto EXIT;
    }

    for (j = 1; j<=nimg; j++) {
       mm = (int) dlist(1,j);
       if( mm != 0) {
          fangles(1,j)=-dlist(3,j)+360.0;
          if ( mm  < 0) {
             mm = -mm;
             fangles(1,j)=fangles(1,j)+180.0;
             fangles(2,j)=180.0-angrefs(2,mm);
             fangles(3,j)=angrefs(1,mm)+180.0;
             if ( fangles(1,j) >= 360.0) fangles(1,j)=fangles(1,j)-360.0;
             if ( fangles(3,j) >= 360.0) fangles(3,j)=fangles(3,j)-360.0;
          }
          else {
             fangles(2,j)=angrefs(2,mm);
             fangles(3,j)=angrefs(1,mm);
          }
       }
    } 
  
 EXIT:
    if (numr)      free(numr);
    if (refcstk)   free(refcstk);
    if (imgwindow) free(imgwindow);
    if (imgcirc)   free(imgcirc);
    if (totmin)    free(totmin);
    if (totmir)    free(totmir);
    if (tot)       free(tot);
    if (tmt)       free(tmt);
    if (dlist)     free(dlist);

    return status;
}

//----------------------------------------------------------------
int aprq2d(float   *sqimg, float     *bfc, int     *numr, int      nsam, 
           int       nrow, int   ishrange, int     istep, int       nsb, 
           int        nse, int        nrb, int       nre, int     lcirc, 
           int      nring, int     maxrin, int      nima, char     mode, 
           float *refdirs, float  *expdir, float   range, float  *diref, 
           float   *ccrot, float *rangnew, float *xshsum, float *yshsum,  
           int   *nimalcg, int   ckmirror, int limitrange)
/*
c 
c  parameters:
c                diref    number of  most similar ref. proj.  (output)
c                            (negative if mirrored)
c                ccrot    corr coeff.                         (output)
c                rangnew  inplane angle                       (output)
c                xshsum   shift                               (output)
c                yshsum   shift                               (output)
c                nimalcg                                      (output)
c
	dimension a(nsam,nrow),bfc(lcirc,nima),numr(3,nring) 
	double precision  fitp(-1:1,-1:1)
	double precision, dimension(*)    :: tt
	real, dimension(3,nima)           :: refdirs
	real, dimension(3)                :: expdir
        automatic arrays
	double precision  fit(-istep:istep,-istep:istep)
	dimension         rotmp(-istep:istep,-istep:istep)
        real, dimension(lcirc)             :: a_circ

*/
{
   float  *imgcirc;
   int    *lcg; 

   double ccrotd,peak,tota,tmta,tmt;
   int    mirrored;

   double quadpi=3.14159265358979323846;
   double dgr_to_rad = quadpi/180.0;

   int    imi, iend, mwant, jtma, itma;
   float  dt, dtabs, rangnewt;
   int    jt, it, irr, ir, ibe, isx, isy, status, idis;
   float  cnr2, cns2, co, so, afit, sx, sy, tot;

   double fitp[9], *fit;
   float  *rotmp;
   int    fitsize;

   status = 0;
 
   imgcirc = (float*)calloc(lcirc,sizeof(float));
   if (imgcirc == NULL) {
       fprintf(stderr,"aprq2d: failed to allocate imgcirc\n");
       status = -1;
       goto EXIT;
   }

   fitsize = (2*istep+1)*(2*istep+1);
   if ( istep >= 1) {
       fit   = (double*)calloc(fitsize,sizeof(double));
       rotmp = (float*) calloc(fitsize,sizeof(float));
   }
   else {
       status = -2;
       goto EXIT;
   }

   peak = 0.0;
   iend  = nima;

   if (limitrange) {
      // restricted range search
      lcg  = (int*) calloc(nima, sizeof(int));
      if (!lcg) {
         mwant = nima;
         status = -1;
         fprintf(stderr,"lcg: %d\n", mwant);
         fprintf(stderr, "range: %g, nima: %d\n", range, nima);
         goto EXIT;
      }

      *nimalcg = 0;
      for (imi=1; imi<=nima; imi++) {
         // dt near 1.0 = not-mirrored, dt near -1.0 = mirrored
         dt    = (expdir(1) * refdirs(1,imi) + 
                  expdir(2) * refdirs(2,imi) +
                  expdir(3) * refdirs(3,imi));
         dtabs = fabs(dt);

         if (dtabs >= range) {
            // mirored or non-mirrored is within range
            *nimalcg++;
            lcg(*nimalcg) = imi;
            if (dt < 0) lcg(*nimalcg) = -imi;
         }
      }

      if (*nimalcg <= 0) {
         // there is no reference within search range
         *xshsum  = 0;
         *yshsum  = 0;
         *diref   = 0;
         *rangnew = 0;
         *ccrot   = -1.0;
         goto EXIT; 
      }
      iend = *nimalcg;
      // end of restricted range search
   }

   ccrotd = -1.0e23;

   // go through centers for shift alignment
   for (jt=-ishrange;jt<=ishrange;jt=jt+istep) {
      cnr2 = nrow/2+1+jt;
      for (it=-ishrange;it<=ishrange;it=it+istep) {
         cns2 = nsam/2+1+it;

         // normalize under the mask
         //'normalize' image values over variance range
         normass(sqimg,nsam,nsb-it,nse-it,nrb-jt,nre-jt,numr(1,1),
                 numr(1,nring));

         // interpolation into polar coordinates
         // creates imgcirc (exp. image circles) for this position

         alrq_ms(sqimg,nsam,nrow,cns2,cnr2,numr,imgcirc,lcirc,nring,mode);

         // creates fourier of: a_circ
         frngs(imgcirc,numr,nring);

         // compare exp. image with all reference images
         for (irr=1;irr<=iend;irr++) {
            ir = irr;
            if (limitrange) ir = abs(lcg(irr));
              
            if (ckmirror) {
               if (limitrange) {
                  mirrored = 0;
		  if (lcg(irr) < 0) mirrored = 1;
                  // check either mirrored or non-mirrored position 
                  crosrng_e(&bfc(1,ir),imgcirc,lcirc,nring,
                            maxrin,numr,&tota,&tot,mirrored);
               }
               else {
                  // check both non-mirrored & mirrored positions 
	          status = crosrng_ms(&bfc(1,ir),imgcirc,lcirc,nring,
                                      maxrin,numr,&tota,&tot,&tmta,&tmt);
               }
            }
            else {
               // do not check mirrored position
		mirrored = 0;
	        crosrng_e(&bfc(1,ir),imgcirc,lcirc,nring,
                          maxrin,numr,&tota,&tot,mirrored);
            }

            if (tota >= ccrotd) {
               // good match with tota (mirrored or not)  position 
	       ccrotd  = tota;
	       ibe     = ir;
	       isx     = it;
	       isy     = jt;
	       *rangnew = ang_n(tot,mode,maxrin);
	       idis    = ir;
               if (limitrange && lcg(irr) < 0) idis = -ir;
	    }

            if (ckmirror && !limitrange) {
               // have to compare with mirrored position 
               if (tmta >= ccrotd) {
                  // good match with mirrored position 
	          ccrotd  = tmta;
	          ibe     = ir;
	          isx     = it;
	          isy     = jt;
	          *rangnew =  ang_n(tmt,mode,maxrin);
	          idis    = -ir;
	       }
            }
         } // endfor irr 
      } // endfor it
   } //  endfor jt

   // try to interpolate
   *ccrot = ccrotd;
   sx     = isx;
   sy     = isy;
   *diref = idis;

   // do not interpolate for point on the edge
   if ( abs(isx) != ishrange && abs(isy) != ishrange) {
      // have to find neighbouring values
      fit(0,0)   = ccrotd;
      rotmp(0,0) = *rangnew;

      for (jt=-istep;jt<=istep;jt++) {
         for (it=-istep;it<=istep;it++) {
            if (it !=0 || jt != 0) {
               cnr2 = nrow/2+1+jt+isy;
               cns2 = nsam/2+1+it+isx;

               normass(sqimg, nsam, nsb-(it+isx),nse-(it+isx),
                       nrb-(jt+isy),nre-(jt+isy), numr(1,1), numr(1,nring));

               alrq_ms(sqimg,nsam,nrow,cns2,cnr2,numr,imgcirc,lcirc,
                       nring,mode);

               frngs(imgcirc,numr,nring);

               //  if (idis .lt. 0)  check mirrored only
               mirrored = 0;
               if (idis < 0) mirrored = 1;
               crosrng_e(&bfc(1,ibe),imgcirc,lcirc,nring,
                         maxrin,numr,&fit(it,jt),&rotmp(it,jt),
                         mirrored);
               rotmp(it,jt) = ang_n(rotmp(it,jt),mode,maxrin);
            }
	 } // endfor it
      } //endfor jt 

      //  find the maximum within +/-istep
      //  maximum cannot be on the edge, i.e., it,jt/=istep
      afit     = fit(0,0);
      jtma     = 0;
      itma     = 0;
      rangnewt = rotmp(0,0);
      if ( istep > 1) {
         for (jt=-istep+1;jt<=istep-1;jt++) {
            for (it=-istep+1;it<=istep-1;it++) {
               if (fit(it,jt) > afit) {
                  afit     = fit(it,jt);
                  rangnewt = rotmp(it,jt);
                  itma     = it;
                  jtma     = jt;
               }
            }
         } 
      }
      //  temp variable overcomes compiler bug on opt 64 pgi 6.0
      *rangnew = rangnewt;

      //  copy values around the peak.
      for (jt=-1;jt<=1;jt++) 
         for (it=-1;it<=1;it++)
            fitp(it,jt) = fit(itma+it,jtma+jt);

      //  update location of the peak
      ccrotd = afit;
      isx    = isx+itma;
      isy    = isy+jtma;
      parabld(fitp,&sx,&sy,&peak);

      //  check whether interpolation is ok.
      if (fabs(sx) < 1.0 && fabs(sy) < 1.0) {
         //  not on edge of 3x3 area
         sx   = sx+isx;
         sy   = sy+isy;
         cnr2 = nrow/2+1+sy;
         cns2 = nsam/2+1+sx;

         normass(sqimg,nsam,nsb-isx,nse-isx,nrb-isy,nre-isy,numr(1,1),
                 numr(1,nring));

         alrq_ms(sqimg,nsam,nrow,cns2,cnr2,numr,imgcirc,lcirc,nring,mode);

         frngs(imgcirc,numr,nring);

         mirrored = 0;
         if (idis < 0) mirrored = 1;

         crosrng_e(&bfc(1,ibe),imgcirc,lcirc,nring,
                   maxrin,numr,&ccrotd,rangnew,mirrored);

         *ccrot   = ccrotd;
         *rangnew = ang_n(*rangnew,mode,maxrin);
      } 
      else {
         //  not on edge of 3x3 area
         sx = isx;
         sy = isy;
      }
   }

   sx = -sx;
   sy = -sy;

   // now have to change order of shift & rotation.
   // in this program image is shifted first, rotated second.
   // in 'rt sq' it is rotation first, shift second.
   // this part corresponds to 'sa p'.
   co      =  cos((*rangnew) * dgr_to_rad);
   so      = -sin((*rangnew) * dgr_to_rad);
   *xshsum = sx*co - sy*so;
   *yshsum = sx*so + sy*co;

   free(fit);
   free(rotmp);
   free(imgcirc);
   if (limitrange) free(lcg);

EXIT:
   return status;
}


//-----------------------------------------------
float ang_n(float rkk, char mode, int maxrin)
{
    float ang; 

    if (mode == 'H' || mode == 'h') {
	ang = fmod(((rkk-1.0) / maxrin+1.0)*180.0, 180.0);
    }
    else if ( mode == 'F' || mode == 'f') {
	ang = fmod(((rkk-1.0) / maxrin+1.0)*360.0, 360.0);
    }
    return ang;
}
//-----------------------------------------------
void alrq_ms(float *xim, int    nsam, int  nrow, float cns2, float cnr2,
             int  *numr, float *circ, int lcirc, int  nring, char  mode)
{
   double dpi, dfi;
   int    it, jt, inr, l, nsim, kcirc, lt;
   float  yq, xold, yold, fi, x, y;

   //     cns2 and cnr2 are predefined centers
   //     no need to set to zero, all elements are defined

   dpi = 2*atan(1.0);
   for (it=1;it<=nring;it++) {
      // radius of the ring
      inr = numr(1,it);
      yq  = inr;

      l = numr(3,it);
      if ( mode == 'h' || mode == 'H' ) { 
         lt = l / 2;
      }
      else if ( mode == 'f' || mode == 'F' ) {
         lt = l / 4;
      } 

      nsim  = lt - 1;
      dfi   = dpi / (nsim+1);
      kcirc = numr(2,it);
      xold  = 0.0;
      yold  = inr;

      circ(kcirc) = quadri(xold+cns2,yold+cnr2,nsam,nrow,xim);

      xold  = inr;
      yold  = 0.0;
      circ(lt+kcirc) = quadri(xold+cns2,yold+cnr2,nsam,nrow,xim);

      if ( mode == 'f' || mode == 'F' ) {
         xold = 0.0;
         yold = -inr;
         circ(lt+lt+kcirc) = quadri(xold+cns2,yold+cnr2,nsam,nrow,xim);

         xold = -inr;
         yold = 0.0;
         circ(lt+lt+lt+kcirc) = quadri(xold+cns2,yold+cnr2,nsam,nrow,xim);
      }
      
      for (jt=1;jt<=nsim;jt++) {
         fi   = dfi * jt;
         x    = sin(fi) * yq;
         y    = cos(fi) * yq;

         xold = x;
         yold = y;
         circ(jt+kcirc) = quadri(xold+cns2,yold+cnr2,nsam,nrow,xim);

         xold = y;
         yold = -x;
         circ(jt+lt+kcirc) = quadri(xold+cns2,yold+cnr2,nsam,nrow,xim);

         if ( mode == 'f' || mode == 'F' ) {
            xold = -x;
            yold = -y;
            circ(jt+lt+lt+kcirc) = quadri(xold+cns2,yold+cnr2,nsam,nrow,xim);

            xold = -y;
            yold = x;
            circ(jt+lt+lt+lt+kcirc) = quadri(xold+cns2,yold+cnr2,nsam,nrow,xim);
         }
      } // end for jt
   } //end for it
}
//-----------------------------------------------
void parabld(double *z33, float *xsh, float *ysh, double *peakv)
{
/*
c parabld  9/25/81 : parabolic fit to 3 by 3 peak neighborhood
c double precision version 
c
c the formula for paraboloid to be fiited into the nine points is:
c
c	f = c1 + c2*y + c3*y**2 + c4*x + c5*xy + c6*x**2
c
c the values of the coefficients c1 - c6 on the basis of the
c nine points around the peak, as evaluated by altran:
*/
   double c1,c2,c3,c4,c5,c6,denom;
   float  xmin, ymin;

   c1 = (26.*z33(1,1) - z33(1,2) + 2*z33(1,3) - z33(2,1) - 19.*z33(2,2)
         -7.*z33(2,3) + 2.*z33(3,1) - 7.*z33(3,2) + 14.*z33(3,3))/9.0;

   c2 = (8.* z33(1,1) - 8.*z33(1,2) + 5.*z33(2,1) - 8.*z33(2,2) + 3.*z33(2,3)
        +2.*z33(3,1) - 8.*z33(3,2) + 6.*z33(3,3))/(-6.);

   c3 = (z33(1,1) - 2.*z33(1,2) + z33(1,3) + z33(2,1) -2.*z33(2,2)
        + z33(2,3) + z33(3,1) - 2.*z33(3,2) + z33(3,3))/6.0;

   c4 = (8.*z33(1,1) + 5.*z33(1,2) + 2.*z33(1,3) -8.*z33(2,1) -8.*z33(2,2)
       - 8.*z33(2,3) + 3.*z33(3,2) + 6.*z33(3,3))/(-6.0);

   c5 = (z33(1,1) - z33(1,3) - z33(3,1) + z33(3,3))/4.0;

   c6 = (z33(1,1) + z33(1,2) + z33(1,3) - 2.*z33(2,1) - 2.*z33(2,2)
        -2.*z33(2,3) + z33(3,1) + z33(3,2) + z33(3,3))/6.0;

   // the peak coordinates of the paraboloid can now be evaluated as:

   *ysh=0.0;
   *xsh=0.0;
   denom=4.*c3*c6 - c5*c5;
   if (denom != 0.0) {
      *ysh=(c4*c5 - 2.*c2*c6) /denom-2.0;
      *xsh=(c2*c5 - 2.*c4*c3) /denom-2.0;
      *peakv= 4.*c1*c3*c6 - c1*c5*c5 -c2*c2*c6 + c2*c4*c5 - c4*c4*c3;
      *peakv= *peakv/denom;
      // limit interplation to +/- 1. range
      xmin = min(*xsh,1.0);
      ymin = min(*ysh,1.0);
      *xsh=max(xmin,-1.0);
      *ysh=max(ymin,-1.0);
   } 
}
//-----------------------------------------------
void crosrng_e(float *circ1, float *circ2, int lcirc,
               int    nring, int   maxrin, int *numr,
               double *qn, float *tot, int neg)
{
/*
c checks single position, neg is flag for checking mirrored position
c
c  input - fourier transforms of rings!
c  first set is conjugated (mirrored) if neg
c  circ1 already multiplied by weights!
c       automatic arrays
	dimension         t(maxrin+2)
	double precision  q(maxrin+2)
	double precision  t7(-3:3)
*/
   float *t;
   double t7[7], *q;
   int    i, j, k, ip, jc, numr3i, numr2i, jtot;
   float  pos; 

   ip = maxrin;
   q = (double*)calloc(maxrin+2, sizeof(double));
   t = (float*)calloc(maxrin+2, sizeof(float));
     
   for (i=1;i<=nring;i++) {
      numr3i = numr(3,i);
      numr2i = numr(2,i);

      t(1) = (circ1(numr2i)) * circ2(numr2i);

      if (numr3i != maxrin) {
         // test .ne. first for speed on some compilers
	 t(numr3i+1) = circ1(numr2i+1) * circ2(numr2i+1);
	 t(2)        = 0.0;

         if (neg) {
            // first set is conjugated (mirrored)
	    for (j=3;j<=numr3i;j=j+2) {
	       jc = j+numr2i-1;
	       t(j) =(circ1(jc))*circ2(jc)-(circ1(jc+1))*circ2(jc+1);
	       t(j+1) = -(circ1(jc))*circ2(jc+1)-(circ1(jc+1))*circ2(jc);
	    } 
         } 
         else {
	    for (j=3;j<=numr3i;j=j+2) {
	       jc = j+numr2i-1;
	       t(j) = (circ1(jc))*circ2(jc) + (circ1(jc+1))*circ2(jc+1);
	       t(j+1) = -(circ1(jc))*circ2(jc+1) + (circ1(jc+1))*circ2(jc);
	    }
         } 
         for (j=1;j<=numr3i+1;j++) q(j) = q(j) + t(j);
      }
      else {
	 t(2) = circ1(numr2i+1) * circ2(numr2i+1);
         if (neg) {
            // first set is conjugated (mirrored)
	    for (j=3;j<=maxrin;j=j+2) {
	       jc = j+numr2i-1;
	       t(j) = (circ1(jc))*circ2(jc) - (circ1(jc+1))*circ2(jc+1);
	       t(j+1) = -(circ1(jc))*circ2(jc+1) - (circ1(jc+1))*circ2(jc);
	    }
         }
         else {
	    for (j=3;j<=maxrin;j=j+2) {
	       jc = j+numr2i-1;
	       t(j) = (circ1(jc))*circ2(jc) + (circ1(jc+1))*circ2(jc+1);
	       t(j+1) = -(circ1(jc))*circ2(jc+1) + (circ1(jc+1))*circ2(jc);
	    } 
         }
         for (j = 1; j <= maxrin+2; j++) q(j) = q(j) + t(j);
      }
   }

   fftr_d(q,ip);

   *qn = -1.0e20;
   for (j=1;j<=maxrin;j++) {
      if (q(j) >= *qn) {
         *qn = q(j);
	 jtot = j;
      }
   } 

   for (k=-3;k<=3;k++) {
      j = (jtot+k+maxrin-1)%maxrin + 1;
      t7(k+4) = q(j);
   }

   prb1d(t7,7,&pos);

   *tot = (float)jtot + pos;

   if (q) free(q);
   if (t) free(t);
}
//-----------------------------------------------
int apmq(EMData *refprj, Dict refparams, EMData *expimg, APMQopt options,
         float  *angles, float *shifts)
{
    double quadpi = 3.14159265358979323846;
    double dgr_to_rad =  quadpi/180;

    int    mr, nr, iskip;
    float  range; 
    char   mode;

    int    nsam, nrow, numth, limitrange, maxrin, nimalcg, nring, status;
    int    nx, ny, nxr, nyr, nidi, nima, lcirc, nwsam, nwrow, nsb, nse,
           nrb, nre, it, ishrange, istep, ckmirror, nrad, nang;
    int    iref, ireft, mirrornew, imgref, ngotpar, mirrorold;
    float  ccrot, rangnew, xshnew, yshnew, peakv; 
    float  rangout, angdif, rangold, xshold, yshold, c, s;
 
    float  *refstk, *imgstk, *bfc, *imgwindow, *refdirs, *expdirs, 
           *angexps, *expdir;

    vector <float> angrefs;
    float  *dlist;
    int    ldd = 7;
    int    *numr;

    float  divas;

    status = 0;

    // retrieve APMQ options
    mr       = options.mr;
    nr       = options.nr;
    ishrange = options.shrange;
    istep    = options.istep;
    mode     = options.mode;
    range    = options.range; 

    angexps = options.angexps;

    nx    = expimg->get_xsize();
    ny    = expimg->get_ysize();
    nidi  = expimg->get_zsize();

    nxr   = refprj->get_xsize();
    nyr   = refprj->get_ysize();
    nima  = refprj->get_zsize();

    if (nx != nxr || ny != nyr) { 
	status = -2;   
        goto EXIT;
    }  

    nsam = nx;
    nrow = ny;

    nrad = min(nsam/2-1, nrow/2-1);
    if ( mr <=0 ) {
       fprintf(stderr,"first ring must be > 0\n");
       status = 10;
       goto EXIT;
    }
    if ( nr >= nrad) {
       fprintf(stderr,"last ring must be < %d\n", nrad);
       status = 10;
       goto EXIT;
    }
    if ( (ishrange+nr) > (nrad-1)) {
       fprintf(stderr,"last ring + translation must be < %d\n", nrad);
       status = 10;
       goto EXIT;
    }

    ngotpar   = 0;

    // reference angles are passed in from refparams
    angrefs = refparams["anglelist"];
    nang    = angrefs.size()/3;
    if (nang != nima) {
       fprintf(stderr, "apmd: nang = %d, nima = %d\n", nang, nima);
       status = -3;
       goto EXIT;
    }

    // extract image data from the reference image EM object
    refstk = refprj->get_data();
    //  find number of reference-rings
    iskip = 1; 
    nring = setnring(mr, nr, iskip);

    numr = (int*)calloc(3*nring,sizeof(int));
    if (!numr) {
	fprintf(stderr,"apmq: failed to allocate numr\n");
        status = -1;
        goto EXIT; 
    }

    numrinit(mr, nr, iskip, numr);
    
    //  Calculate pointers for rings and store them in numr
    status = alprbs(numr, nring, &lcirc, mode);
    if (status != 0) goto EXIT;

    maxrin = numr(3,nring);

    // convert reference images to rings
    bfc = (float*)calloc(lcirc*nima,sizeof(float));

    // read reference images into reference rings (bfc) array 
    status = aprings(nima, nring, refstk, nsam, nrow, numr, bfc, lcirc, mode);
    if (status !=0 ) goto EXIT;

    // extract the image data from experimental EM image object
    imgstk = expimg->get_data();

    // find divas, numth, nsam, & nrow
    apmaster_1(mode,&divas,nr,&numth,nsam,nrow,&nwsam,&nwrow);

    // *****STRANGELY, APMQ DOES NOT USE THE WINDOWED IMAGE*****
    //  Because it is not necessary - in the old APMD code that was done
    //     to save space  PAP  
    // allocate work space
    nwsam = nsam;
    nwrow = nrow; 
    imgwindow = (float*)calloc(nwsam*nwrow,sizeof(float));

    // calculate dimensions for normas
    nsb  = -nsam/2;
    nse  = -nsb-1+(nsam%2);
    nrb  = -nrow/2;
    nre  = -nrb-1+(nrow%2);

    limitrange = 0;
    ckmirror   = 1;
    if (range > 0.0 && range < 360.0) limitrange = 1;
    range      = cos(range*dgr_to_rad);
    nimalcg    = 1;

    if (limitrange) {
       // retrieve refangles for restricted angular search
       refdirs = (float*)calloc(3*nima,sizeof(float));
       if (!refdirs) {
          fprintf(stderr, "apmq: failed to allocate refdirs!\n");
          status = -1;
          goto EXIT;
       }
       // convert ref. angles to unitary directional vectors (refdirs).
       // call ap_getsata(angrefs,refdirs,3,nima,irtflg)


       // read previously determined exp. angles into angexps
       angexps = options.angexps;
       if (!angexps) {
          fprintf(stderr, "apmq: no angexps available !\n");
          status = -4;
          goto EXIT;
       }

       expdirs = (float*)calloc(3*nidi, sizeof(float));
       if (!expdirs) {
          fprintf(stderr, "apmq: failed to allocate expdirs!\n");
          status = -1;
          goto EXIT;
       }

       // convert exp. angles to unitary directional vectors(expdirs).
       // call ap_getsata(angexp,expdirs,7,nidi,irtflg)
    }

    dlist   = (float*)calloc(ldd*nidi,sizeof(float));
    if (dlist == NULL || expdirs == NULL) {
	fprintf(stderr, "apmq: failed to allocated dlist or expdirs...\n");
        status = -1;
        goto EXIT;
    }

    // loop over all sets of experimental (sample) images
    for (it=1;it<=nidi;it++) {
       // APMQ DOES NOT WINDOW, IT JUST COPIES imgstk to imgwindow
       win_resize(&imgstk(1,1,it), imgwindow, nsam, nrow, nwsam, nwrow,
                  1, nrow, 1, nsam); 
       printf("it = %d\n", it);

       if (limitrange) expdir = &expdirs(1,it); // otherwise expdir is a dummy

       status = aprq2d(imgwindow   , bfc         , numr ,  nsam   , 
                       nrow        , ishrange    , istep,  nsb    , 
                       nse         , nrb         , nre  ,  lcirc  , 
                       nring       , maxrin      , nima ,  mode   , 
                       refdirs     , expdir      , range,  &dlist(2,it),  
                       &dlist(3,it), &dlist(4,it), &dlist(5,it), 
                       &dlist(6,it), &nimalcg    , ckmirror    , 
                       limitrange);
//debug
//       printf("dlist(2,it) = %11.3e\n", dlist(2,it));
//       printf("dlist(3,it) = %11.3e\n", dlist(3,it));
//       printf("dlist(4,it) = %11.3e\n", dlist(4,it));
//       printf("dlist(5,it) = %11.3e\n", dlist(5,it));
//       printf("dlist(6,it) = %11.3e\n", dlist(6,it));
/*
c      output (in dlist position is increased by 1, no.1 is the key).
c      2 - number of the most similar reference projection.
c      3 - not-normalized correlation coefficient.
c      4 - psi angle. (in=plane rotation)
c      5 - sx shift
c      6 - sy shift
c      7 - input image number.
*/
       //  dlist(2,it) is list number of most similar ref. image 
       //  (<0 if mirrored, 0 if none )
       iref = (int)dlist(2,it);
       if (iref < 0) {
          //  mirrored reference image
          imgref = -iref;

          //      ireft is for refdirs index
          ireft     = -iref;
          mirrornew = 1;
       } 
       else if (iref == 0) {
          // no nearby reference image
          imgref = 0;
          // ireft is for refdirs index
          ireft  = 1;
          mirrornew = 0;
       }
       else {
          imgref = iref;
          //  ireft is for refdirs index
          ireft  = iref;
          mirrornew = 0;
       } 
 
       ccrot    = dlist(3,it);
       rangnew  = dlist(4,it);
       xshnew   = dlist(5,it);
       yshnew   = dlist(6,it);
       peakv    = 1.0;

       // imgref is number of most similar ref. image 
       if (imgref <= 0) {
          // no reference image selected
          ccrot = -1.0;
          peakv = 0.0;
       }
//  THE WHOLE NEXT SECTION IS SUSPICIOUS!  I WOULD REMOVE IT! PAP
       // set new projection angles
       angles(1,it) = 0.0; // default value
       angles(2,it) = 0.0; // 
       angles(2,it) = 0.0; //

       if (imgref > 0) {
          // use ref. angles as new projection angles
          angles(1,it) = angrefs(1,iref);
          angles(2,it) = angrefs(2,iref);
          angles(3,it) = angrefs(3,iref);

          if (mirrornew) {
             // ref. projection must be mirrored
             angles(1,it) = -angles(1,it);
             angles(2,it) = 180+angles(2,it);
          }
       }
       else if (ngotpar >= 3) {
          // keep old exp. proj. angles 
          angles(1,it) = angexps(1,it);
          angles(2,it) = angexps(2,it);
          angles(3,it) = angexps(3,it);
       }

       rangold   = 0.0;
       xshold    = 0.0;
       yshold    = 0.0;

       if (ngotpar >= 7 && ishrange > 0) {
          // use old inplane rot. & shift  
          rangold   = angexps(4,it);
          xshold    = angexps(5,it);
          yshold    = angexps(6,it);

          mirrorold = 0;
          if (angexps(7,it) > 0) mirrorold = 1;
          if (mirrorold) {
             status = 5;
             fprintf(stderr, "mirrorold = %d\n", mirrorold);
             goto EXIT;
          }
       }

      
       // combine rot. & shift with previous transformation
       c  =  cos(rangnew * dgr_to_rad);
       s  = -sin(rangnew * dgr_to_rad);

       shifts(1,it) = xshnew  + xshold*c - yshold*s;
       shifts(2,it) = yshnew  + xshold*s + yshold*c;
       rangout = rangold + rangnew;

       // list angles in range 0...360
       while ( rangout < 0.0) {
          rangout = rangout + 360.0;
       }
       while(rangout >= 360.0) {
          rangout = rangout - 360.0;
       } 

       // set flag for no angdif determined
       angdif = -1.0;

       if (imgref <= 0) {
          //  no relevant ref. image found
          angdif = 0.0;
       }
       else if (ngotpar >= 3) {
          //  can find angdif
          angdif = fabs(expdirs(1,it) * refdirs(1,iref) + 
                        expdirs(2,it) * refdirs(2,iref) + 
                        expdirs(3,it) * refdirs(3,iref));
          angdif = min(1.0,angdif);
          angdif = acos(angdif) / dgr_to_rad;
       }
    } // endfor it

EXIT:
    free(numr);
    if (limitrange) {
       free(refdirs);
       free(expdirs);
    }
    free(dlist);
    free(bfc);
    free(imgwindow);
    return status;
}
