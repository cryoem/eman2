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
 
/** This file is a part of util.h, 
 * To use this file's functions, you should include "util.h"
 * NEVER directly include this file */
 
#ifndef util__sparx_h__
#define util__sparx_h__ 

public:

static int coveig(int n, float *covmat, float *eigval, float *eigvec);

static Dict ExpMinus4YSqr(float ymax,int nsamples);

static void WTM(EMData* PROJ, vector<float> SS,int DIAMETER,int NUMP);

static void WTF(EMData* PROJ,vector<float> SS,float SNR,int K,vector<float> exptable);

static Dict CANG(float PHI, float THETA, float PSI);
 
static void BPCQ(EMData* B, EMData *CUBE,vector<float> DM);

static vector<float> infomask(EMData* Vol, EMData* mask);

static void colreverse(float* beg, float* end, int nx);

static void slicereverse(float* beg, float* end, int nx,int ny);

static  void cyclicshift(EMData* image, Dict params);

static Dict im_diff(EMData* V1, EMData* V2, EMData* mask=0);

/** Creates a Two D Test Pattern
 * @param[in] Size must be odd
 * @param[in] p the x frequency
 * @param[in] q the y frequency
 * @param[in] a the x falloff
 * @param b
 * @param flag
 * @param[in] alphaDeg the projection angle
 * @return The 2D test pattern in real space, fourier space,
 *               or the projection in real or fourier space 
 *               or the FH of the pattern
*/
static EMData* TwoDTestFunc(int Size, float p, float q,  float a, float b, 
                   int flag=0, float alphaDeg=0); //PRB


/** Given a tabulated function y of x (n unordered points), and
 * Given the values of the m values xq to be interpolated
 * This routine returns the interpolated array yq, PRB
 * This function is called by splint 
 * @param x
 * @param[in] y of x is the tabulated function of length n
 * @param n
 * @param[in] xq is the x values to be splined: has m points.
 * @param[out]  yq are the splined values
 * @param m
*/
static void spline_mat(float *x, float *y, int n,  float *xq, float *yq, int m); //PRB

/** Given a tabulated function y of x (unordered), and
 * Given the values of the first derivatives at the end points
 * This routine returns an array y2, that contains the second derivatives
 * of the function at the tabulated points.   PRB
 * This function is called by splint 
 * @param x
 * @param[in] y of x is the tabulated function of length n
 * @param n
 * @param[in] yp1 : the derivatives of the first point.
 * @param[in] ypn : the derivatives of the last point.
 * @param[out]  y2 is the value of the second derivatives
*/
static void spline(float *x, float *y, int n, float yp1, float ypn, float *y2);

/** Given the arrays xa(ordered, ya of length n, which tabulate a function 
 *  and given the array y2a which is the output of spline and an unordered array xq,
 *  this routine returns a cubic-spline interpolated array yq.
 * @param xa
 * @param[in] ya of x is the tabulated function of length n
 * @param[in] y2a is returned from spline: second derivs
 * @param n
 * @param[in] xq is the x values to be splined: has m points.
 * @param[out]  yq are the splined values
 * @param m
*/  // PRB
static void splint( float *xa, float *ya, float *y2a, int n,  
                                     float *xq, float *yq, int m);


/** list the sorted lengths of the  integer lattice 
 * sites of a square sided image of  size Size. PRB
 * @param[out] PermMatTr  The matrix telling the ordering of the lattice 
 *                        sites wrt the array
 * @param [out] kValsSorted
 * @param[out] weightofkvalsSorted the number of sites at that distance
 * @param[in] Size The length of the image 
 * @param [out] SizeReturned
 * 
 */
 static void Radialize(int *PermMatTr,  float * kValsSorted,
            float *weightofkvalsSorted, int Size, int *SizeReturned);
		            

/** 1-D Kaiser-Bessel window function class.
 *  (It's a class so that the windowing parameters may be
 *   instantiated and held in the instance object.)
 *
 *  The I0 version can be tabulated and interpolated upon
 *  demand, but the max error needs to be checked.  The 
 *  "vtable" parameter corresponds to the maximum value of x
 *  for which the I0 window is non-zero.  Setting "vtable"
 *  different from "v" corresponds to a change in units of x.
 *  In practice, it is often handy to replace x in some sort
 *  of absolute units with x described in terms of grid
 *  intervals.

 *
 *  The get_kbsinh_win and get_kbi0_win functions return
 *  single-argument function objects, which is what a 
 *  generic routine is likely to want.
 *
 *  @see P. A. Penczek, R. Renka, and H. Schomberg,
 *      J. Opt. Soc. Am. _21_, 449 (2004)
 * 
 */
class KaiserBessel 
{
	protected:
		float alpha, v, r; /** Kaiser-Bessel parameters */
		int N; /** size in Ix-space */
		int K; /** I0 window size */
		float vtable; /** table I0 non-zero domain maximum */
		int ntable;
		vector<float> i0table;
		float dtable; /** table spacing */
		float alphar; /** alpha*r */
		float fac; /** 2*pi*alpha*r*v */
		float vadjust; 
		float facadj; /** 2*pi*alpha*r*vadjust */
		virtual void build_I0table(); /** Tabulate I0 window for speed */
		float fltb;
	public:
		KaiserBessel(float alpha_, int K, float r_,
				     float v_, int N_, float vtable_=0.f, 
					 int ntable_ = 5999);
		virtual ~KaiserBessel() {};
		/** Compute the maximum error in the table */
		float I0table_maxerror();
		vector<float> dump_table() {
			return i0table;
		}
		/** Kaiser-Bessel Sinh window function */
		virtual float sinhwin(float x) const;
		/** Kaiser-Bessel I0 window function */
		virtual float i0win(float x) const;
		/** Kaiser-Bessel I0 window function (uses table lookup) */
				inline float i0win_tab(float x) const {
					/*float absx = fabs(x);
					int loc = int(round(absx*fltb));
					return i0table[loc];*/
					float xt;
					if(x<0.f) xt = -x*fltb+0.5f; else xt = x*fltb+0.5f;
					return i0table[ (int) xt];
					/*return i0table[ (int) (fabs(x)*fltb+0.5f)];*/
#if 0 // old version
					if (absx > vtable) return 0.f;
					float loc = absx/dtable;
					return i0table[int(loc + 0.5f)];
#endif // 0
		}
		/** Return the size of the I0 window */
		int get_window_size() const { return K; }
		/** Sinh window function object */
		class kbsinh_win {
			KaiserBessel& kb;
			public:
			kbsinh_win(KaiserBessel& kb_) : kb(kb_) {}
			float operator()(float x) const {
				return kb.sinhwin(x);
			}
			int get_window_size() const {return kb.get_window_size();}
		};
		/** Sinh window function object factory */
		kbsinh_win get_kbsinh_win() { 
			return kbsinh_win(*this);
		}
		/** I0 window function object */
		class kbi0_win {
			KaiserBessel& kb;
			public:
			kbi0_win(KaiserBessel& kb_) : kb(kb_) {}
			float operator()(float x) const {
				return kb.i0win(x);
			}
			int get_window_size() const {return kb.get_window_size();}
		};
		/** I0 window function object factory */
		kbi0_win get_kbi0_win() { 
			return kbi0_win(*this);
		}
};

class FakeKaiserBessel : public KaiserBessel {
	public:
		FakeKaiserBessel(float alpha, int K, float r_,
				         float v_, int N_, float vtable_=0.f,
						 int ntable_ = 5999) 
        : KaiserBessel(alpha, K, r_, v_, N_, vtable_, ntable_) {
			build_I0table();
		}
		float sinhwin(float x) const;
		float i0win(float x) const;
		void build_I0table();
};

		/** Compute a vector containing quasi-evenly spaced Euler angles
		 *
		 * The order of angles in the vector is phi, theta, psi.
		 *
		 * @param[in] delta  Delta theta (spacing in theta).
		 * @param[in] t1  Starting (min) value of theta in degrees, default = 0.
		 * @param[in] t2  Ending (max) value of theta in degrees, default = 90.
		 * @param[in] p1  Starting (min) value of phi in degrees, default = 0.
		 * @param[in] p2  Ending (max) value of phi in degrees, default = 359.9.
		 * @return Vector of angles as a flat list of phi_0, theta_0, psi_0, ..., phi_N, theta_N, psi_N.
		 */
		static vector<float>
		even_angles(float delta, float t1=0, float t2=90, float p1=0, float p2=359.999);


		/** Quadratic interpolation (2D).
		 *
		 *  Note:  This routine starts counting from 1, not 0!
		 *
		 *	This routine uses six image points for interpolation:
		 *
		 *@see M. Abramowitz & I.E. Stegun, Handbook of Mathematical
		 *     Functions (Dover, New York, 1964), Sec. 25.2.67.
		 *     http://www.math.sfu.ca/~cbm/aands/page_882.htm.
		 *
		 *@see http://www.cl.cam.ac.uk/users/nad/pubs/quad.pdf
		 *
		 *@verbatim
                f3    fc
                |
                | x
         f2-----f0----f1
                |
                |
                f4
		 *@endverbatim
		 *
		 *	f0 - f4 are image values near the interpolated point X.
		 *	f0 is the interior mesh point nearest x.
		 *
		 *	Coords:  
		 *@li        f0 = (x0, y0)
		 *@li        f1 = (xb, y0)
		 *@li        f2 = (xa, y0)
		 *@li        f3 = (x0, yb)
		 *@li        f4 = (x0, ya)
		 *@li        fc = (xc, yc)
		 *
		 *	Mesh spacings: 
		 *@li              hxa -- x- mesh spacing to the left of f0
		 *@li              hxb -- x- mesh spacing to the right of f0
		 *@li              hyb -- y- mesh spacing above f0
		 *@li              hya -- y- mesh spacing below f0
		 *
		 *	Interpolant:
		 *	  f = f0 + c1*(x-x0) + c2*(x-x0)*(x-x1)
		 *			 + c3*(y-y0) + c4*(y-y0)*(y-y1)
		 *			 + c5*(x-x0)*(y-y0)
		 *
		 *	@param[in] x x-coord value
		 *	@param[in] y y-coord value
		 *  @param nx
		 *  @param ny
		 *	@param[in] image Image object (pointer)
		 *
		 *	@return Interpolated value
		 */
		static float quadri(float x, float y, int nx, int ny, float* image);
		/*static float quadris(float x, float y, int nx, int ny, float* image);*/
		static float bilinear(float xold, float yold, int nsam, int nrow, float* xim);


		/** Quadratic interpolation (3D)
		 * @param r
		 * @param s
		 * @param t
		 * @param fdata
		 *
		 * @return Interpolated value
		*/
		static float triquad(float r, float s, float t, float* fdata);

		/** Gaussian function class.
		 *
		 *  Usage:
		 *
		 *     Gaussian gauss(sigma);
		 *     float g = gauss(x);
		 */
		class Gaussian {
			float sigma;
			float rttwopisigma;
			float twosigma2;
			public:
			Gaussian(float sigma_ = 1.0) : sigma(sigma_) {
				rttwopisigma = sqrtf(static_cast<float>(twopi)*sigma);
				twosigma2 = 2*sigma*sigma;
			}
			inline float operator()(float x) const {
				return exp(-x*x/(twosigma2))/rttwopisigma;
			}
		};
		static void alrq(float *xim,  int nsam , int nrow , int *numr,
                                 float *circ, int lcirc, int nring, char mode);
        static EMData* Polar2D(EMData* image, vector<int> numr, string mode);
        static EMData* Polar2Dm(EMData* image, float cns2, float cnr2, vector<int> numr, string mode);
        static void alrq_ms(float *xim, int    nsam, int  nrow, float cns2, float cnr2,
                            int  *numr, float *circ, int lcirc, int  nring, char  mode);
        static void alrl_ms(float *xim, int    nsam, int  nrow, float cns2, float cnr2,
                            int  *numr, float *circ, int lcirc, int  nring, char  mode);
        static void alrq_msi(EMData* image,float cns2, float cnr2,
                           int  *numr, float *circ, int lcirc, int  nring, char  mode, Util::KaiserBessel& kb);
        static EMData* Polar2Dmi(EMData* image, float cns2, float cnr2, vector<int> numr, string mode, Util::KaiserBessel& kb);

        static void  fftr_q(float  *xcmplx, int nv);
        static void  fftr_d(double *xcmplx, int nv);
        static void  fftc_q(float  *br, float  *bi, int ln, int ks);
        static void  fftc_d(double *br, double *bi, int ln, int ks);
        static void  Frngs(EMData* circ, vector<int> numr);
        static void  frngs(float *circ, int *numr, int nring);
		static boost::tuple<double, float, int>
			Crosrng_e(EMData* circ1, EMData* circ2, vector<int> numr, int neg);
        static void  crosrng_e(float *circ1, float *circ2, int lcirc,
                               int    nring, int   maxrin, int *numr,
                               double *qn, float *tot, int neg);
			       
        static Dict Crosrng_ms(EMData* circ1, EMData* circ2, vector<int> numr);
        static void  crosrng_ms(float *circ1, float *circ2, int  lcirc, int  nring,
                                int   maxrin, int   *numr , double *qn, float *tot,
                                double   *qm, float *tmt);
        static Dict Crosrng_msr(EMData* circ1, EMData* circ2, vector<int> numr);
        static void  crosrng_msr(float *circ1, float *circ2, int  lcirc, int  nring,
                                int   maxrin, int   *numr , float *qn, float *tot,
                                float   *qm, float *tmt);
				
        static EMData* Crosrng_msg(EMData* circ1, EMData* circ2, vector<int> numr);
        static void  crosrng_msg(float *circ1, float *circ2, double *q, double *t, int  lcirc, int  nring,
                                int   maxrin, int   *numr );
	
        static void  prb1d(double *b, int npoint, float *pos);
	
	/* Decimates the image with respect to the image center.
	 * (i.e) the center of the original image is kept the same 
	 * and then the initial start pixel is calculated with respect to the 
	 * center of the image
	 * @params(image, x-pixel, y-pixel,z-pixel)
	 * works for all 3 dimensions
	**/
	static EMData* decimate(EMData* img, int x_step,int y_step=1,int z_step=1);
	
	static EMData* window(EMData* img,int new_nx ,int new_ny=1, int new_nz=1, int x_offset=0, int y_offset=0, int z_offset=0);

	static EMData* pad(EMData* img, int new_nx, int new_ny=1, int new_nz=1, int x_offset=0, int y_offset=0, int z_offset=0, char *params="average");

	/*static void histogram(EMData* image, EMData* mask);*/
	
	static Dict histc(EMData *ref,EMData *img,EMData *mask);
	
	static float hist_comp_freq(float PA,float PB,int size_img, int hist_len, EMData *img, vector<float> ref_freq_hist, EMData *mask, float ref_h_diff, float ref_h_min);
	/* The unit in the ctf function: dz: Angstrom, cs: CM  Ps: Angstrom, Voltage: Kv,dza: Angstrom, azz: degree wgh: None unit. b_factor: Angstrom^2 
	 The CTF function takes form of   *sin(-quadpi*(dz*lambda*ak^2-cs*lambda^3*ak^4/2.)-wgh)*exp(-b_factor*ak^2)*sign
          * sign can be set as +1 or -1 . The unit of frequency ak is 1/Angstrom
                  Attention: Envelope function in power spectrum has a form of exp(-b_factor*ak^2)
                                          */ 
	static EMData* ctf_img(int nx, int ny, int nz, float ps,float dz,float cs=2.0f,float voltage=100.0f,float dza=0.0f,float azz=0.0f,float wgh=.1,float b_factor=0.0f,float sign=-1.0f);
        static float tf(float dzz,float ak,float lambda,float cs,float wgh,float b_factor,float sign);
	static EMData *compress_image_mask(EMData* image, EMData* mask);
	static EMData *reconstitute_image_mask(EMData *image,EMData *mask);
	static vector<float> merge_peaks(vector<float> peak1, vector<float> peak2,float p_size);
	static vector<float> pw_extract(vector<float>pw, int n, int iswi,float ps);
	static vector<float> call_cl1(long int *k,long int *n, float *ps, long int *iswi, float *pw, float *q2, double *q, double *x, double *res, double *cu, double *s, long int *iu);
	static vector<float> lsfit(long int *ks,long int *n, long int *klm2d, long int *iswi, float *q1,double *q, double *x, double *res, double *cu, double *s,long int *iu);
	static void cl1(long int *k, long int *l, long int *m, long int *n, long int *klm2d,double *q, double *x, double *res, double *cu, long
	int *iu, double *s);
	static float eval(char * images,EMData * img, vector<int> S,int N, int K,int size);
	/*  VORONOI DIAGRAM */
	static vector<double> vrdg(EMData *th, EMData *ph);
	static void hsortd(double *theta,double *phi,int *key,int len,int option);
	static void voronoidiag(vector<double> theta,vector<double> phi,vector<double> weight,int n);
	static void voronoidiag(double *theta,double *phi,double* weight,int n);
	static void angstep(double* thetast,int len);
	static void voronoi(double *phi,double *theta,double *weight,int lenw,int low,int medium,int nt,int last);
	static void disorder2(double *x,double *y,double *z,int *key,int len);
	static void ang_to_xyz(double *x,double *y,double *z,int len);
	static double dot_product(double *x,double *y);
	static void flip23(double *x,double *y,double *z,int *key,int k,int len);
	static void swap(double x,double y);
	struct tmpstruct{
 		double theta1,phi1;
		int key1;
 		};
	static bool cmp1(tmpstruct tmp1,tmpstruct tmp2);
	static bool cmp2(tmpstruct tmp1,tmpstruct tmp2);
	/**********************************************************/
	/* ######### STRIDPACK USED COMMANDS FOR VORONOI #########################*/
	static int trmsh3_(int *n0, double *tol, double *x, double *y, double *z__, int *n, int *list, int *lptr, 
	       int *lend, int *lnew, int *indx, int *lcnt, int *near__, int *next, double *dist, int *ier);
	static double areav_(int *k, int *n, double *x, double *y, double *z__, int *list, int *lptr, int *lend, int *ier);
#endif	//util__sparx_h__
