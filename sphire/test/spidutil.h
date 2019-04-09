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

#define  angles(i,j)     angles [((j)-1)*3 + (i)-1]
#define  newangles(i,j)  newangles [((j)-1)*3 + (i)-1]
#define  angrefs(i,j)    angrefs [((j)-1)*3 + (i)-1]
#define  angexps(i,j)    angexps [((j)-1)*3 + (i)-1]
#define  shifts(i,j)     shifts[((j)-1)*2 + (i)-1]
#define  fangles(i,j)    fangles [((j)-1)*3 + (i)-1]
#define  numr(i,j)       numr   [((j)-1)*3 + (i)-1]
#define  imgcirc(i)      imgcirc[(i)-1]
#define  wr(i)           wr     [(i)-1]
#define  sqimg(i,j)      sqimg  [((j)-1)*nsam + (i)-1]
#define  xim(i,j)        xim    [((j)-1)*nsam + (i)-1]
#define  fdata(i,j)      fdata  [((j)-1)*nxdata + (i)-1] 
#define  min0(a, b)      ((a) >= (b) ? (b) : (a))
#define  min(a, b)       ((a) >= (b) ? (b) : (a))
#define  max(a, b)       ((a) <= (b) ? (b) : (a))
#define  tab1(i)         tab1   [(i)-1]
#define  xcmplx(i,j)     xcmplx [((j)-1)*2 + (i)-1]
#define  br(i)           br     [(i)-1]
#define  bi(i)           bi     [(i)-1]
#define  circ(i)         circ   [(i)-1]
#define  circ1(i)        circ1  [(i)-1]
#define  circ2(i)        circ2  [(i)-1]
#define  t(i)            t      [(i)-1]
#define  q(i)            q      [(i)-1]
#define  b(i)            b      [(i)-1]
#define  t7(i)           t7     [(i)-1]
#define  imgfrom(i,j)    imgfrom[((j)-1)*lsam + (i)-1]
#define  imgto(i,j)      imgto  [((j)-1)*nsam + (i)-1]
#define  imgstk(i,j,k)   imgstk[((k)-1)*nsam*nrow + ((j)-1)*nsam + (i)-1]
#define  refcstk(i,j)    refcstk[((j)-1)*lcirc + (i) - 1]
#define  imgwindow(i,j)  imgwindow [((j)-1)*nwsam + (i)-1]
#define  totmin(i)       totmin[(i)-1]
#define  totmir(i)       totmir[(i)-1]
#define  tot(i)          tot[(i)-1]
#define  tmt(i)          tmt[(i)-1]
#define  dlist(i,j)      dlist[((j)-1)*ldd + (i)-1]
#define  expdir(i)       expdir[(i)-1]
#define  expdirs(i,j)    expdirs[((j)-1)*3 + (i)-1]
#define  refdirs(i,j)    refdirs[((j)-1)*3 + (i)-1]
#define  refdir(i)       refdir[(i)-1]
#define  lcg(i)          lcg[(i)-1]
#define  bfc(i,j)        bfc[((j)-1)*lcirc + (i) - 1]
#define  fitp(i,j)       fitp[ ((j)+1)*3 + (i) + 1]
#define  fit(i,j)        fit[((j)+istep)*(2*istep+1) + (i) + istep]
#define  rotmp(i,j)      rotmp[((j)+istep)*(2*istep+1) + (i) + istep]
#define  z33(i,j)        z33[((j)-1)*3 + (i)-1]

struct APMQopt {
    int nr;      // first ring 
    int mr;      // last  ring
    int shrange; // shift search range
    int istep;   // shift search stepsize
    char mode;  
    float range; // angular search range
    float *angexps; // previously determined Euler angles
    float *angrefs; // reference angles for restricted angular search;
}; 

struct APMDopt {
    int nr;      // first ring
    int mr;      // last ring
    int iskip; 
    char mode;
    float range;
}; 

int  aprings(int nimg, int nring, float *imgstk, int nsam, 
             int nrow, int *numr, float *refcstk,int lcirc, 
             char mode);
int   apring1(float *sqimg, int  nsam, int  nrow, float *imgcirc, int lcirc, 
              int nring   , char mode, int *numr, float *wr);
float quadri(float xx, float yy, int nxdata, int nydata, float *fdata);
void  ringwe(float *wr, int *numr, int nring);
int   alprbs(int *numr, int nring, int *lcirc, char mode);
int   setnring(int mr, int nr, int iskip);
void  numrinit(int mr, int nr, int iskip, int *numr);
void  normass(float *sqimg, int nsam, int ns1, int ns2, int nr1, int nr2, 
              int ir1, int ir2);
void  normas(float *sqimg, int nsam, int ns1, int ns2, int nr1, int nr2,
             int ir1, int ir2);
void  frngs(float *circ, int *numr, int nring);
void  fftr_q(float *xcmplx, int nv);
void  fftc_q(float *br, float *bi, int ln, int ks);
int   applyws(float *circ, int lcirc, int *numr, float *wr,
              int nring);
void  alrq(float *xim,  int nsam , int nrow , int *numr,
           float *circ, int lcirc, int nring, char mode);
void  alrq_ms(float *xim, int    nsam, int  nrow, float cns2, float cnr2,
              int  *numr, float *circ, int lcirc, int  nring, char  mode);
int   crosrng_ms(float *circ1, float *circ2, int  lcirc, int   nring,
                 int   maxrin, int   *numr, double *qn, double *tot,
                 double   *qm, double *tmt);
void  prb1d(double *b, int npoint, float *pos);
void  apmaster_1(char mode, float *divas, int nr, int *numth,
                 int  lsam, int lrow, int *nsam, int *nrow);
void  win_resize(float *imgfrom, float *imgto, int lsam, int lrow, 
                 int nsam, int nrow, int lr1, int lr2, int ls1, int ls2);
int   apmd(EMData *refprj, Dict refparam, EMData *expimg, APMDopt options,
           float *fangle);
float ang_n(float rkk, char mode, int maxrin);
void parabld(double *z33, float *xsh, float *ysh, double *peakv);
void crosrng_e(float *circ1, float *circ2, int lcirc,
               int    nring, int   maxrin, int *numr,
               double *qn, float *tot, int neg);

int apmq(EMData *refprj, Dict refparams, EMData *expimg, APMQopt options,
         float *angles, float *shifts);

int aprq2d(float  *sqimg, float     *bfc, int     *numr, int      nsam, 
           int      nrow, int   ishrange, int     istep, int       nsb, 
           int       nse, int        nrb, int       nre, int     lcirc, 
           int     nring, int     maxrin, int      nima, char     mode, 
           float *refdir, float  *expdir, float   range, float  *diref, 
           float  *ccrot, float *rangnew, float *xshsum, float *yshsum, 
           int  *nimalcg, int   ckmirror, int limitrange);
