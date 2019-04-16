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

using namespace EMAN;

using std::cout;
using std::cin;
using std::endl;
using std::string;

#include "spidfft.h"

void  fftr_q(float *xcmplx, int nv) 
{
   // dimension xcmplx(2,1); xcmplx(1,i) --- real, xcmplx(2,i) --- imaginary

   float tab1[15];
   int nu, inv, nu1, n, isub, n2, i1, i2, i;
   float ss, cc, c, s, tr, ti, tr1, tr2, ti1, ti2, t;

   tab1(1)=9.58737990959775e-5;
   tab1(2)=1.91747597310703e-4;
   tab1(3)=3.83495187571395e-4;
   tab1(4)=7.66990318742704e-4;
   tab1(5)=1.53398018628476e-3;
   tab1(6)=3.06795676296598e-3;
   tab1(7)=6.13588464915449e-3;
   tab1(8)=1.22715382857199e-2;
   tab1(9)=2.45412285229123e-2;
   tab1(10)=4.90676743274181e-2;
   tab1(11)=9.80171403295604e-2;
   tab1(12)=1.95090322016128e-1;
   tab1(13)=3.82683432365090e-1;
   tab1(14)=7.07106781186546e-1;
   tab1(15)=1.00000000000000;

   nu=abs(nv);
   inv=nv/nu;
   nu1=nu-1;
   n=(int)pow(2,nu1);
   isub=16-nu1;

   ss=-tab1(isub);
   cc=-2.0*pow(tab1(isub-1),2);
   c=1.0;
   s=0.0;
   n2=n/2;
   if ( inv > 0) {
      fftc_q(&xcmplx(1,1),&xcmplx(2,1),nu1,2);
      tr=xcmplx(1,1);
      ti=xcmplx(2,1);
      xcmplx(1,1)=tr+ti;
      xcmplx(2,1)=tr-ti;
      for (i=1;i<=n2;i++) {
         i1=i+1;
         i2=n-i+1;
         tr1=xcmplx(1,i1);
         tr2=xcmplx(1,i2);
         ti1=xcmplx(2,i1);
         ti2=xcmplx(2,i2);
         t=(cc*c-ss*s)+c;
         s=(cc*s+ss*c)+s;
         c=t;
         xcmplx(1,i1)=0.5*((tr1+tr2)+(ti1+ti2)*c-(tr1-tr2)*s);
         xcmplx(1,i2)=0.5*((tr1+tr2)-(ti1+ti2)*c+(tr1-tr2)*s);
         xcmplx(2,i1)=0.5*((ti1-ti2)-(ti1+ti2)*s-(tr1-tr2)*c);
         xcmplx(2,i2)=0.5*(-(ti1-ti2)-(ti1+ti2)*s-(tr1-tr2)*c);
     }
   }
   else {
     tr=xcmplx(1,1);
     ti=xcmplx(2,1);
     xcmplx(1,1)=0.5*(tr+ti);
     xcmplx(2,1)=0.5*(tr-ti);
     for (i=1; i<=n2; i++) {
        i1=i+1;
        i2=n-i+1;
        tr1=xcmplx(1,i1);
        tr2=xcmplx(1,i2);
        ti1=xcmplx(2,i1);
        ti2=xcmplx(2,i2);
        t=(cc*c-ss*s)+c;
        s=(cc*s+ss*c)+s;
        c=t;
        xcmplx(1,i1)=0.5*((tr1+tr2)-(tr1-tr2)*s-(ti1+ti2)*c);
        xcmplx(1,i2)=0.5*((tr1+tr2)+(tr1-tr2)*s+(ti1+ti2)*c);
        xcmplx(2,i1)=0.5*((ti1-ti2)+(tr1-tr2)*c-(ti1+ti2)*s);
        xcmplx(2,i2)=0.5*(-(ti1-ti2)+(tr1-tr2)*c-(ti1+ti2)*s);
     }
     fftc_q(&xcmplx(1,1),&xcmplx(2,1),nu1,-2);
   }
}

// -----------------------------------------------------------------
void fftc_q(float *br, float *bi, int ln, int ks)
{
   //  dimension  br(1),bi(1)

   int b3,b4,b5,b6,b7,b56;
   int n, k, l, j, i, ix0, ix1; 
   float rni, tr1, ti1, tr2, ti2, cc, c, ss, s, t, x2, x3, x4, x5, sgn;
   float tab1[15]; 
   int status=0;

   tab1(1)=9.58737990959775e-5;
   tab1(2)=1.91747597310703e-4;
   tab1(3)=3.83495187571395e-4;
   tab1(4)=7.66990318742704e-4;
   tab1(5)=1.53398018628476e-3;
   tab1(6)=3.06795676296598e-3;
   tab1(7)=6.13588464915449e-3;
   tab1(8)=1.22715382857199e-2;
   tab1(9)=2.45412285229123e-2;
   tab1(10)=4.90676743274181e-2;
   tab1(11)=9.80171403295604e-2;
   tab1(12)=1.95090322016128e-1;
   tab1(13)=3.82683432365090e-1;
   tab1(14)=7.07106781186546e-1;
   tab1(15)=1.00000000000000;

   n=(int)pow(2,ln);

   k=abs(ks);
   l=16-ln;
   b3=n*k;
   b6=b3;
   b7=k;
   if( ks > 0 ) {
      sgn=1.0;
   } 
   else {
      sgn=-1.0;
      rni=1.0/(float)n;
      j=1;
      for (i=1; i<=n;i++) {
         br(j)=br(j)*rni;
         bi(j)=bi(j)*rni;
         j=j+k;
      }
   }
L12:
   b6=b6/2;
   b5=b6;
   b4=2*b6;
   b56=b5-b6;
L14:
   tr1=br(b5+1);
   ti1=bi(b5+1);
   tr2=br(b56+1);
   ti2=bi(b56+1);

   br(b5+1)=tr2-tr1;
   bi(b5+1)=ti2-ti1;
   br(b56+1)=tr1+tr2;
   bi(b56+1)=ti1+ti2;

   b5=b5+b4;
   b56=b5-b6;
   if (b5 <= b3)  goto  L14;
   if (b6 == b7)  goto  L20;

   b4=b7;
   cc=2.0*pow(tab1(l),2);
   c=1.0-cc;
   l=l+1;
   ss=sgn*tab1(l);
   s=ss;
L16: 
   b5=b6+b4;
   b4=2*b6;
   b56=b5-b6;
L18:
   tr1=br(b5+1);
   ti1=bi(b5+1);
   tr2=br(b56+1);
   ti2=bi(b56+1);
   br(b5+1)=c*(tr2-tr1)-s*(ti2-ti1);
   bi(b5+1)=s*(tr2-tr1)+c*(ti2-ti1);
   br(b56+1)=tr1+tr2;
   bi(b56+1)=ti1+ti2;

   b5=b5+b4;
   b56=b5-b6;
   if(b5 <= b3)  goto L18;
   b4=b5-b6;
   b5=b4-b3;
   c=-c;
   b4=b6-b5;
   if(b5 < b4)  goto  L16;
   b4=b4+b7;
   if(b4 >= b5) goto  L12;

   t=c-cc*c-ss*s;
   s=s+ss*c-cc*s;
   c=t;
   goto  L16;
L20:
   ix0=b3/2;
   b3=b3-b7;
   b4=0;
   b5=0;
   b6=ix0;
   ix1=0;
   if ( b6 == b7) goto EXIT;
L22:
   b4=b3-b4;
   b5=b3-b5;
   x2=br(b4+1);
   x3=br(b5+1);
   x4=bi(b4+1);
   x5=bi(b5+1);
   br(b4+1)=x3;
   br(b5+1)=x2;
   bi(b4+1)=x5;
   bi(b5+1)=x4;
   if (b6 < b4) goto  L22;
L24:
   b4=b4+b7;
   b5=b6+b5;
   x2=br(b4+1);
   x3=br(b5+1);
   x4=bi(b4+1);
   x5=bi(b5+1);
   br(b4+1)=x3;
   br(b5+1)=x2;
   bi(b4+1)=x5;
   bi(b5+1)=x4;
   ix0=b6;
L26:
   ix0=ix0/2;
   ix1=ix1-ix0;
   if(ix1 >= 0)  goto  L26;

   ix0=2*ix0;
   b4=b4+b7;
   ix1=ix1+ix0;
   b5=ix1;
   if (b5 >= b4)  goto  L22;
   if (b4 < b6)   goto  L24;
EXIT:
   status = 0; 
}

// -------------------------------------------
void  fftr_d(double *xcmplx, int nv) 
{
   // double precision  x(2,1)
   int    i1, i2,  nu, inv, nu1, n, isub, n2, i;
   double tr1,tr2,ti1,ti2,tr,ti;
   double cc,c,ss,s,t;
   double tab1[15];

   tab1(1)=9.58737990959775e-5;
   tab1(2)=1.91747597310703e-4;
   tab1(3)=3.83495187571395e-4;
   tab1(4)=7.66990318742704e-4;
   tab1(5)=1.53398018628476e-3;
   tab1(6)=3.06795676296598e-3;
   tab1(7)=6.13588464915449e-3;
   tab1(8)=1.22715382857199e-2;
   tab1(9)=2.45412285229123e-2;
   tab1(10)=4.90676743274181e-2;
   tab1(11)=9.80171403295604e-2;
   tab1(12)=1.95090322016128e-1;
   tab1(13)=3.82683432365090e-1;
   tab1(14)=7.07106781186546e-1;
   tab1(15)=1.00000000000000;

   nu=abs(nv);
   inv=nv/nu;
   nu1=nu-1;
   n=(int)pow(2,nu1);
   isub=16-nu1;
   ss=-tab1(isub);
   cc=-2.0*pow(tab1(isub-1),2);
   c=1.0;
   s=0.0;
   n2=n/2;

   if ( inv > 0 ) {
      fftc_d(&xcmplx(1,1),&xcmplx(2,1),nu1,2);
      tr=xcmplx(1,1);
      ti=xcmplx(2,1);
      xcmplx(1,1)=tr+ti;
      xcmplx(2,1)=tr-ti;
      for (i=1;i<=n2;i++) {
         i1=i+1;
         i2=n-i+1;
         tr1=xcmplx(1,i1);
         tr2=xcmplx(1,i2);
         ti1=xcmplx(2,i1);
         ti2=xcmplx(2,i2);
         t=(cc*c-ss*s)+c;
         s=(cc*s+ss*c)+s;
         c=t;
         xcmplx(1,i1)=0.5*((tr1+tr2)+(ti1+ti2)*c-(tr1-tr2)*s);
         xcmplx(1,i2)=0.5*((tr1+tr2)-(ti1+ti2)*c+(tr1-tr2)*s);
         xcmplx(2,i1)=0.5*((ti1-ti2)-(ti1+ti2)*s-(tr1-tr2)*c);
         xcmplx(2,i2)=0.5*(-(ti1-ti2)-(ti1+ti2)*s-(tr1-tr2)*c);
      }
   }
   else {
      tr=xcmplx(1,1);
      ti=xcmplx(2,1);
      xcmplx(1,1)=0.5*(tr+ti);
      xcmplx(2,1)=0.5*(tr-ti);
      for (i=1;i<=n2;i++) {
         i1=i+1;
         i2=n-i+1;
         tr1=xcmplx(1,i1);
         tr2=xcmplx(1,i2);
         ti1=xcmplx(2,i1);
         ti2=xcmplx(2,i2);
         t=(cc*c-ss*s)+c;
         s=(cc*s+ss*c)+s;
         c=t;
         xcmplx(1,i1)=0.5*((tr1+tr2)-(tr1-tr2)*s-(ti1+ti2)*c);
         xcmplx(1,i2)=0.5*((tr1+tr2)+(tr1-tr2)*s+(ti1+ti2)*c);
         xcmplx(2,i1)=0.5*((ti1-ti2)+(tr1-tr2)*c-(ti1+ti2)*s);
         xcmplx(2,i2)=0.5*(-(ti1-ti2)+(tr1-tr2)*c-(ti1+ti2)*s);
      } 
      fftc_d(&xcmplx(1,1),&xcmplx(2,1),nu1,-2);
   } 
} 
//-----------------------------------------
void fftc_d(double *br, double *bi, int ln, int ks)
{
   double rni,sgn,tr1,tr2,ti1,ti2;
   double cc,c,ss,s,t,x2,x3,x4,x5;
   int    b3,b4,b5,b6,b7,b56;
   int    n, k, l, j, i, ix0, ix1, status=0; 
   double tab1[15];

   tab1(1)=9.58737990959775e-5;
   tab1(2)=1.91747597310703e-4;
   tab1(3)=3.83495187571395e-4;
   tab1(4)=7.66990318742704e-4;
   tab1(5)=1.53398018628476e-3;
   tab1(6)=3.06795676296598e-3;
   tab1(7)=6.13588464915449e-3;
   tab1(8)=1.22715382857199e-2;
   tab1(9)=2.45412285229123e-2;
   tab1(10)=4.90676743274181e-2;
   tab1(11)=9.80171403295604e-2;
   tab1(12)=1.95090322016128e-1;
   tab1(13)=3.82683432365090e-1;
   tab1(14)=7.07106781186546e-1;
   tab1(15)=1.00000000000000;

   n=(int)pow(2,ln);

   k=abs(ks);
   l=16-ln;
   b3=n*k;
   b6=b3;
   b7=k;
   if (ks > 0) {
      sgn=1.0;
   }
   else {
      sgn=-1.0;
      rni=1.0/(float)(n);
      j=1;
      for (i=1;i<=n;i++) {
         br(j)=br(j)*rni;
         bi(j)=bi(j)*rni;
         j=j+k;
      }
   }

L12:
   b6=b6/2;
   b5=b6;
   b4=2*b6;
   b56=b5-b6;

L14:
   tr1=br(b5+1);
   ti1=bi(b5+1);
   tr2=br(b56+1);
   ti2=bi(b56+1);

   br(b5+1)=tr2-tr1;
   bi(b5+1)=ti2-ti1;
   br(b56+1)=tr1+tr2;
   bi(b56+1)=ti1+ti2;

   b5=b5+b4;
   b56=b5-b6;
   if ( b5 <= b3 )  goto  L14;
   if ( b6 == b7 )  goto  L20;

   b4=b7;
   cc=2.0*pow(tab1(l),2);
   c=1.0-cc;
   l++;
   ss=sgn*tab1(l);
   s=ss;

L16:
   b5=b6+b4;
   b4=2*b6;
   b56=b5-b6;

L18:
   tr1=br(b5+1);
   ti1=bi(b5+1);
   tr2=br(b56+1);
   ti2=bi(b56+1);
   br(b5+1)=c*(tr2-tr1)-s*(ti2-ti1);
   bi(b5+1)=s*(tr2-tr1)+c*(ti2-ti1);
   br(b56+1)=tr1+tr2;
   bi(b56+1)=ti1+ti2;

   b5=b5+b4;
   b56=b5-b6;
   if ( b5 <= b3 )  goto  L18;
   b4=b5-b6;
   b5=b4-b3;
   c=-c;
   b4=b6-b5;
   if ( b5 < b4 )  goto  L16;
   b4=b4+b7;
   if ( b4 >= b5 ) goto  L12;

   t=c-cc*c-ss*s;
   s=s+ss*c-cc*s;
   c=t;
   goto  L16;

L20:
   ix0=b3/2;
   b3=b3-b7;
   b4=0;
   b5=0;
   b6=ix0;
   ix1=0;
   if (b6 == b7) goto EXIT;

L22:
   b4=b3-b4;
   b5=b3-b5;
   x2=br(b4+1);
   x3=br(b5+1);
   x4=bi(b4+1);
   x5=bi(b5+1);
   br(b4+1)=x3;
   br(b5+1)=x2;
   bi(b4+1)=x5;
   bi(b5+1)=x4;
   if(b6 < b4)  goto  L22;

L24:
   b4=b4+b7;
   b5=b6+b5;
   x2=br(b4+1);
   x3=br(b5+1);
   x4=bi(b4+1);
   x5=bi(b5+1);
   br(b4+1)=x3;
   br(b5+1)=x2;
   bi(b4+1)=x5;
   bi(b5+1)=x4;
   ix0=b6;

L26:
   ix0=ix0/2;
   ix1=ix1-ix0;
   if( ix1 >= 0)  goto L26;

   ix0=2*ix0;
   b4=b4+b7;
   ix1=ix1+ix0;
   b5=ix1;
   if ( b5 >= b4)  goto  L22;
   if ( b4 < b6)   goto  L24;

EXIT:
   status = 0;
} 

