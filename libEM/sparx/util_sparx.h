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

/** same function than Util::coveig but wrapped to use directly in python code */
/* Functions WTF, WTM and BPCQ accept as first parameter images in real and Fourier space. Output image is always in real space. */
static Dict coveig_for_py(int ncov, const vector<float>& covmatpy);

static void WTF(EMData* PROJ,vector<float> SS,float SNR,int K);

static void WTM(EMData* PROJ, vector<float> SS,int DIAMETER,int NUMP);

static Dict CANG(float PHI, float THETA, float PSI);

static void BPCQ(EMData* B, EMData *CUBE, const int radius);

static vector<float> infomask(EMData* Vol, EMData* mask, bool);

static vector<double> helixshiftali(vector<EMData*> ctx, vector<vector<float> > pcoords, int nsegms, float maxincline, int kang, int search_rng, int nxc);

static vector<double> snakeshiftali(vector<EMData*> sccf,  vector<vector<float> > pcoords, int nsegms, float maxincline, int kang, int search_rng, int nxc, int angnxc);

static vector<float> curhelixshiftali(vector<EMData*> ctx, vector<vector<float> > pcoords, int nsegms,  int search_rng, int nx, int ny);

static void colreverse(float* beg, float* end, int nx);

static void slicereverse(float* beg, float* end, int nx,int ny);

/**
 * Performs inplace integer cyclic shift as specified by the "dx","dy","dz" parameters on a 3d volume.
 * Implements the inplace swapping using reversals as descibed in  also:
 * http://www.csse.monash.edu.au/~lloyd/tildeAlgDS/Intro/Eg01/
 * @author  Phani Ivatury
 * @date 18-2006
 * @see http://www.csse.monash.edu.au/~lloyd/tildeAlgDS/Intro/Eg01/
 *
 *
 * A[0] A[1] A[2] A[3] A[4] A[5] A[6] A[7] A[8] A[9]
 *
 * 10   20   30   40   50   60   70   80   90   100
 * ------------
 *   m = 3 (shift left three places)
 *
 *   Reverse the items from 0..m-1 and m..N-1:
 *
 *   30   20   10   100  90   80   70   60   50   40
 *
 *   Now reverse the entire sequence:
 *
 *   40   50   60   70   80   90   100  10   20   30
 *
 *
 *   cycl_shift() in libpy/fundementals.py calls this function
 *
 *   Usage:
 *   EMData *im1 = new EMData();
 *   im1->set_size(70,80,85);
 *   im1->to_one();
 *   Dict params; params["dx"] = 10;params["dy"] = 10000;params["dz"] = -10;
 *   Utils::cyclicshift(im1,params);
 *   im1.peak_search(1,1)
 *   */
static void cyclicshift(EMData* image, Dict params);

static Dict im_diff(EMData* V1, EMData* V2, EMData* mask=0);

/** Creates a Two D Test Pattern
 * @param[in] Size must be odd
 * @param[in] p the x frequency
 * @param[in] q the y frequency
 * @param[in] a the x falloff
 * @param b the y falloff
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


static float bsplineBase(float u);
static float bsplineBasedu(float u);
static void convertTocubicbsplineCoeffs(vector<float> s, int DataLength, float EPSILON);

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


typedef struct
{
    // data structure used here NXCt3rIDWnWrG2bj
    float ims1;
    float ims2;
    float ims3;
} Ims;


class sincBlackman
{
	protected:
		int M; /** kernel size */
		float fc; /** cut-off frequency */
		int ntable;
		vector<float> sBtable;
		virtual void build_sBtable(); /** Tabulate kernel for speed */
		float fltb;
	public:
		sincBlackman(int M_, float fc_, int ntable_ = 1999);
		virtual ~sincBlackman() {};

		inline  float sBwin_tab(float x) const {
			float xt;
			if(x<0.0f) xt = -x*fltb+0.5f; else xt = x*fltb+0.5f;
			return sBtable[ (int) xt];
		}
		/** Return the size of the kernel */
		int get_sB_size() const { return M; }
};



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
			/*return i0table[ (int) (fabs(x)*fltb+0.5f)];
				if (absx > vtable) return 0.f;
				float loc = absx/dtable;
				return i0table[int(loc + 0.5f)]; */
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
		 *     http://www.math.sfu.ca/~cbm/aands/page_882.htm
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

        /** Quadratic interpolation (2D).
		 *
		 *  This is identical to quadri except the wrap around is not done circulantly.
		 *
		 *  Note:  This routine starts counting from 1, not 0!
		 *
		 *	This routine uses six image points for interpolation:
		 *
		 *@see M. Abramowitz & I.E. Stegun, Handbook of Mathematical
		 *     Functions (Dover, New York, 1964), Sec. 25.2.67.
		 *     http://www.math.sfu.ca/~cbm/aands/page_882.htm
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
		static float quadri_background(float x, float y, int nx, int ny, float* image, int xnew, int ynew);

		// Here counting is in C style, so coordinates of the pixel delx should be [0-nx-1]
		/* Commented by Zhengfan Yang on 04/20/07
		This function is written to replace get_pixel_conv(), which is too slow to use in practice.
		I made the following changes to get_pixel_conv():
		1. Use the same data passing scheme as quadri() and move the function from emdata_sparx.cpp to util_sparx.cpp
		2. Reduce usage of i0win_tab (from 98 calls to 14 calls in 2D case, from 1029 calls to 21 calls in 3D case!)
		3. Unfold the 'for' loop
		4. Reduce the usage of multiplications through some bracketing (from 98 times to 57 times in 2D case, from 1029 times to 400 times in 3D case)

		The shortcoming of this routine is that it only works for window size N=7. In case you want to use other window
		size, say N=5, you can easily modify it by referring my code.
		*/
		static float get_pixel_conv_new(int nx, int ny, int nz, float delx, float dely, float delz, float* data, Util::KaiserBessel& kb);

		// Here counting is in C style, so coordinates of the pixel delx should be [0-nx-1]
		/* Commented by Zhengfan Yang on 04/20/07
		This function is written to replace get_pixel_conv(), which is too slow to use in practice.
		I made the following changes to get_pixel_conv():
		1. Use the same data passing scheme as quadri() and move the function from emdata_sparx.cpp to util_sparx.cpp
		2. Reduce usage of i0win_tab (from 98 calls to 14 calls in 2D case, from 1029 calls to 21 calls in 3D case!)
		3. Unfold the 'for' loop
		4. Reduce the usage of multiplications through some bracketing (from 98 times to 57 times in 2D case, from 1029 times to 400 times in 3D case)

		The shortcoming of this routine is that it only works for window size N=7. In case you want to use other window
		size, say N=5, you can easily modify it by referring my code.
		*/
        static float get_pixel_conv_new_background(int nx, int ny, int nz, float delx, float dely, float delz, float* data, Util::KaiserBessel& kb, int xnew, int ynew);

		static std::complex<float> extractpoint2(int nx, int ny, float nuxnew, float nuynew, EMData *fimage, Util::KaiserBessel& kb);

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
	/*static void alrq(float *xim,  int nsam , int nrow , int *numr,
                             float *circ, int lcirc, int nring, char mode);*/
	static EMData* Polar2D(EMData* image, vector<int> numr, string mode);
	static EMData* Polar2Dm(EMData* image, float cns2, float cnr2, vector<int> numr, string cmode);
	/*static void alrq_ms(float *xim, int	 nsam, int  nrow, float cns2, float cnr2,
			    int  *numr, float *circ, int lcirc, int  nring, char  mode);*/
	static void alrl_ms(float *xim, int    nsam, int  nrow, float cns2, float cnr2,
			    int  *numr, float *circ, int lcirc, int  nring, char  mode);
	/*static void alrq_msi(EMData* image,float cns2, float cnr2,
			   int  *numr, float *circ, int lcirc, int  nring, char  mode, Util::KaiserBessel&
	    				       kb);*/
	static EMData* Polar2Dmi(EMData* image, float cns2, float cnr2, vector<int> numr, string cmode, Util::KaiserBessel& kb);

	static void  fftr_q(float  *xcmplx, int nv);
	static void  fftr_d(double *xcmplx, int nv);
	static void  fftc_q(float  *br, float  *bi, int ln, int ks);
	static void  fftc_d(double *br, double *bi, int ln, int ks);

	/** This function conducts the Single Precision Fourier Transform for a set of rings */
	static void  Frngs(EMData* circ, vector<int> numr);
	static void  Normalize_ring(EMData* ring, const vector<int>& numr);

	/** This function conducts the Single Precision Inverse Fourier Transform for a set of rings */
	static void  Frngs_inv(EMData* circ, vector<int> numr);

	/** This is a copy of Applyws routine from alignment.py */
	static void  Applyws(EMData* circ, vector<int> numr, vector<float> wr);

	/*
	  	A little note about different Crosrng:
		Basically, they all do cross-correlation function to two images in polar coordinates
		Crosrng_e is the original one
		Crosrng_ew is the one that you could apply weights to different rings
		Crosrng_ms assumes the user already apply weights to circ1, it also returns both
		           straight and mirrored positions simultaneously.
	        Crosrng_msg differs from the previous ones in that it returns the cross-correlation
			    function entirely instead of the peak value and position, thus makes it
			    possible to use the gridding method to determine the peak position
	        Crosrng_msg_s is same as Crosrng_msg except that it only checks straight position
	        Crosrng_msg_m is same as Crosrng_msg except that it only checks mirrored position
	  */
	static Dict Crosrng_e(EMData* circ1, EMData* circ2, vector<int> numr, int neg);
	static Dict Crosrng_rand_e(EMData* circ1, EMData* circ2, vector<int> numr, int neg, float previous_max, float an, int psi_pos);
	static Dict Crosrng_ew(EMData* circ1, EMData* circ2, vector<int> numr, vector<float> w, int neg);

	static Dict Crosrng_ms(EMData* circ1, EMData* circ2, vector<int> numr);
	static Dict Crosrng_ms_delta(EMData* circ1, EMData* circ2, vector<int> numr, float delta_start, float delta);

	/**
	 * checks either straight or mirrored position depending on flag
	 * input - fourier transforms of rings!!
	 * circ1 already multiplied by weights!
	*/
	static Dict Crosrng_sm_psi(EMData* circ1, EMData* circ2, vector<int> numr, float psi, int flag, float psimax);

    /**
	 * checks both straight & mirrored position
	 * input - fourier transforms of rings!!
	 * circ1 already multiplied by weights!
	*/
	static Dict Crosrng_psi(EMData* circ1, EMData* circ2, vector<int> numr, float psi, float psimax);
 
	/**
	 * checks both straight & mirrored positions
	 * Find only positions within +/-pis_max of 0 and 180 (for helcal use)
	 * input - fourier transforms of rings!
	 * circ1 already multiplied by weights!
	*/
	static Dict Crosrng_ns(EMData* circ1, EMData* circ2, vector<int> numr);

	/**
	 * checks both straight & mirrored positions
	 * input - fourier transforms of rings!!
	 * circ1 already multiplied by weights!
	 * returns EM object with 1D ccf
	*/
	static EMData* Crosrng_msg(EMData* circ1, EMData* circ2, vector<int> numr);

	/**
	 * checks both straight & mirrored positions
	 * input - fourier transforms of rings!!
	 * circ1 already multiplied by weights!
	 * returns EM object with 1D ccf
	*/
	static void Crosrng_msg_vec(EMData* circ1, EMData* circ2, vector<int> numr, float *q, float *t);
	static void Crosrng_msg_vec_snake(EMData* circ1, EMData* circ2, vector<int> numr, float *q);
	/**
	 * This program is half of the Crosrng_msg. It only checks straight position.
	 * input - fourier transforms of rings!!
	 * circ1 already multiplied by weights!
	 * returns EM object with 1D ccf
	*/
	static EMData* Crosrng_msg_s(EMData* circ1, EMData* circ2, vector<int> numr);

	/**
	 * This program is half of the Crosrng_msg. It only checks mirrored position.
	 * input - fourier transforms of rings!!
	 * circ1 already multiplied by weights!
	 * returns EM object with 1D ccf
	*/
	static EMData* Crosrng_msg_m(EMData* circ1, EMData* circ2, vector<int> numr);

	static vector<float> Crosrng_msg_vec_p(EMData* circ1, EMData* circ2, vector<int> numr );

	static void  prb3p(double *b, float *pos);
	static void  prb1d(double *b, int npoint, float *pos);

	static void update_fav(EMData* ave,EMData* dat, float tot, int mirror, vector<int> numr);
	static void sub_fav(EMData* ave,EMData* dat, float tot, int mirror, vector<int> numr);

	// helper functions for ali2d_ra
	static float ener(EMData* ave, vector<int> numr);

	static float ener_tot(const vector<EMData*>& data, vector<int> numr, vector<float> tot);

	/** k-means helper */
	static Dict min_dist_real(EMData* image, const vector<EMData*>& data);

	/** helper function for k-means */
	static Dict min_dist_four(EMData* image, const vector<EMData*>& data);

	/** helper to create the contengency table for partition matching (k-means)
	 * flag define is the list of stable obj must be store to stb, but the size st
	 * must be know before. The trick is first start wihtout the flag to get number
	 * of elements stable, then again with the flag to get the list. This avoid to
	 * have two differents functions for the same thing.
	 * */
	static int k_means_cont_table_(int* group1, int* group2, int* stb, long int s1, long int s2, int flag);

	// branch and bound matching algorithm

	
	/** initial_prune removes all classes C from Parts where there does not exist ANY feasible matching containing class C which has weight gt T.
	 * The first element of each class is its original index, and second is dummy variable
	*/
	static void initial_prune(vector <vector <int*> > & Parts, int* dimClasses, int nParts, int K, int T);

	/** Each class in Parts has its dummy variable set to 0 or 1. Only consider those with its dummy variable set to 1 (the 'active' ones)
	 * First element of each class is its original index, second is the dummy variable slot.
	*/
	static bool explore(vector <vector <int*> > & Parts, int* dimClasses, int nParts, int K, int T, int partref, int* curintx, int
			size_curintx, int* next, int size_next, int depth);

	
	

	/** make an intelligent "guess" at the largest weight of all possible feasible matches.
	 * we make "n_guesses" guesses and return the largest one.
	 * the largest weight of all feasible matches is guaranteed to be larger than or equal to the returned guess.
	*/
	static int generatesubmax(int* argParts, int* Indices, int* dimClasses, int nParts, int K, int T, int n_guesses, int LARGEST_CLASS);

	/** return the weight of the largest weighted feasible match, along with the match in the (preallocated) argument curmax.
	 * The returned weight has to be gt newT. If there is no such feasible matching, return 0 as *curmax
	*/
	static void search2(int* argParts, int* Indices, int* dimClasses, int nParts, int K, int T, int* matchlist, int* costlist, int J);
	
	static void explore2(int* argParts, int* Indices, int* dimClasses, int nParts, int K, int T, int* curintx, int size_curintx, int* next, int size_next, int depth, int J, int* matchlist, int*
costlist, int* curbranch);
	
	/** First element of output is total cost of the matches in the output
	 * Second element of output is the total number of matches in output
	 * So output has 2+(*(output+1))nParts elements.
	*/
	static bool sanitycheck(int* argParts, int* Indices, int* dimClasses, int nParts, int K, int T, int* output);

	/** K is the number of classes in each partition (should be the same for all partitions)
	 * the first element of each class is its original index in the partition, and second is dummy var
	 * MPI: if nTop <= 0, then initial prune is called, and the pruned partitions are returned in a 1D array. The first element is reserved for max_levels (the size of the smallest
	 * partition after pruning).
	 * if nTop > 0, then partitions are assumed to have been pruned, where only dummy variables of un-pruned partitions are set to 1, and findTopLargest is called
	 * to find the top weighted matches. The matches, where each match is preceded by its cost, is returned in a one dimensional vector.
	 *
	 * essentially the same as bb_enumerate but with the option to do mpi version.
	*/
	static vector<int> bb_enumerateMPI_(int* argParts, int* dimClasses, int nParts, int K, int T, int n_guesses, int LARGEST_CLASS, int J, int max_branching, float stmult,
	int branchfunc, int LIM);

	
	/** same as branch except the nFirst (=Levels[0]) possibilites for the first match are already chosen
	 * firstmatches stores the matches and corresponding cost, where each match is preceded by its cost....
	 * output is an int array, the first element is the cost of the output solution, the second element is the total number of matches in the solution
	 * and the rest is the list of matches. output is in one dimensional form.
	*/
	static int* branchMPI(int* argParts, int* Indices, int* dimClasses, int nParts, int K, int T, int curlevel,int n_guesses, int LARGEST_CLASS, int J, int max_branching,
	float stmult, int branchfunc, int LIM);
	
	static int branch_factor_0(int* costlist, int* matchlist, int J, int T, int nParts, int curlevel, int max_branching, int LIM);
	static int branch_factor_2(int* costlist, int* matchlist, int J, int T, int nParts, int curlevel, int max_branching, int LIM);
	static int branch_factor_3(int* costlist, int* matchlist, int J, int T, int nParts, int curlevel, int max_branching, int K, int LIM);
	static int branch_factor_4(int* costlist, int* matchlist, int J, int T, int nParts, int curlevel, int max_branching, float stmult);
	// new code common-lines

	//helper function for the weights calculation by Voronoi to Cml
	static vector<double> cml_weights(const vector<float>& cml);

	/** 2009-03-25 15:35:53 JB. This function calculates common-lines between sinogram */
	static vector<int> cml_line_insino(vector<float> Rot, int i_prj, int n_prj);

	/** 2009-03-30 15:35:07 JB. This function calculates all common-lines between sinogram */
	static vector<int> cml_line_insino_all(vector<float> Rot, vector<int> seq, int n_prj, int n_lines);

	/** 2009-03-25 15:35:05 JB. This function prepare rotation matrix for common-lines */
	static vector<double> cml_init_rot(vector<float> Ori);

	/** 2009-03-25 15:35:37 JB. this function update only one rotation amtrix according a new orientation*/
	static vector<float> cml_update_rot(vector<float> Rot, int iprj, float nph, float th, float nps);

	/** 2009-03-26 10:46:14 JB. This function calculate all common-lines in space for Voronoi */
	static vector<double> cml_line_in3d(vector<float> Ori, vector<int> seq, int nprj, int nlines);

	/** 2009-03-26 11:37:53 JB. This function spin all angle psi and evaluate the partial discrepancy belong common-lines */
	static vector<double> cml_spin_psi(const vector<EMData*>& data, vector<int> com, vector<float> weights, int iprj, vector<int> iw, int n_psi, int d_psi, int n_prj);
	static vector<double> cml_spin_psi_now(const vector<EMData*>& data, vector<int> com, int iprj, vector<int> iw, int n_psi, int d_psi, int n_prj);

	/** 2009-03-30 15:44:05 JB. Compute the discrepancy belong all common-lines */
	static double cml_disc(const vector<EMData*>& data, vector<int> com, vector<int> seq, vector<float> weights, int n_lines);

	/**  This function drop a line (line) to an 2D image (img).
	 *  The position of the line to the image is defined by (postline).
	 *  The part of the line paste is defined by (offset), the begin position
	 *  and (length) the size.
	 */
	static void set_line(EMData* img, int posline, EMData* line, int offset, int length);

	/** This function prepare the line from sinogram by cutting off some frequencies,
	 * and creating the mirror part (complexe conjugate of the first part). Then
	 * both lines (mirror and without) are drop to the sinogram.
	 * line is in Fourrier space, ilf low frequency, ihf high frequency, nblines
	 * number of lines of the half sinogram (the non miror part), sino the sinogram,
	 * pos_line the position of the line in the sino.
	 */
	static void cml_prepare_line(EMData* sino, EMData* line, int ilf, int ihf, int pos_line, int nblines);

	/* Decimates the image with respect to the image center.
	 * (i.e) the center of the original image is kept the same
	 * and then the initial start pixel is calculated with respect to the
	 * center of the image
	 * @params(image, x-pixel, y-pixel,z-pixel)
	 * works for all 3 dimensions
	**/
	static EMData* decimate(EMData* img, int x_step,int y_step=1,int z_step=1);

	static EMData* window(EMData* img,int new_nx ,int new_ny=1, int new_nz=1, int x_offset=0, int y_offset=0, int z_offset=0);

	static EMData* pad(EMData* img, int new_nx, int new_ny=1, int new_nz=1, int x_offset=0, int y_offset=0, int z_offset=0, const char *params="average");

	static vector<float> histogram(EMData* image, EMData* mask, int nbins = 128, float hmin =0.0f, float hmax = 0.0f );

	static Dict histc(EMData *ref,EMData *img,EMData *mask);

	static float hist_comp_freq(float PA,float PB,size_t size_img, int hist_len, EMData *img, vector<float> ref_freq_hist, EMData *mask, float ref_h_diff, float ref_h_min);


	/* The unit in the ctf function: dz: Angstrom, cs: CM  Ps: Angstrom, Voltage: Kv,dza: Angstrom, azz: degree wgh: None unit. b_factor: Angstrom^2
	 The CTF function takes form of   *sin(-quadpi*(dz*lambda*ak^2-cs*lambda^3*ak^4/2.)-wgh)*exp(-b_factor*ak^2)*sign
          * sign can be set as +1 or -1 . The unit of frequency ak is 1/Angstrom
                  Attention: Envelope function in power spectrum has a form of exp(-b_factor*ak^2)
                                          */
	static float   tf(float dzz, float ak, float voltage = 300.0f, float cs = 2.0f, float wgh = 0.1f, float b_factor = 0.0f, float sign = -1.0f);
	static EMData *compress_image_mask(EMData* image, EMData* mask);

	/** Recreates a n-d image using its compressed 1-D form and the mask */
	static EMData *reconstitute_image_mask(EMData *image,EMData *mask);

	static vector<float> merge_peaks(vector<float> peak1, vector<float> peak2,float p_size);
	static vector<float> pw_extract(vector<float>pw, int n, int iswi,float ps);
	static vector<float> call_cl1(long int *k,long int *n, float *ps, long int *iswi, float *pw, float *q2, double *q, double *x, double *res, double *cu, double *s, long int *iu);
	static vector<float> lsfit(long int *ks,long int *n, long int *klm2d, long int *iswi, float *q1,double *q, double *x, double *res, double *cu, double *s,long int *iu);
	static void cl1(long int *k, long int *l, long int *m, long int *n, long int *klm2d,double *q, double *x, double *res, double *cu, long
	int *iu, double *s);
	static float eval(char * images,EMData * img, vector<int> S,int N, int K,int size);

	/*  VORONOI DIAGRAM */
	static vector<double> vrdg(const vector<float>& ph, const vector<float>& th);
	static void hsortd(double *theta,double *phi,int *key,int len,int option);
	static void voronoidiag(double *theta,double *phi,double* weight,int n);
	/*static void angstep(double* thetast,int len);*/
	/*static void voronoi(double *phi,double *theta,double *weight,int lenw,int low,int medium,int nt,int last);*/
	static void voronoi(double *phi,double *theta,double *weight, int nt);
	static void disorder2(double *x,double *y,int *key,int len);
	static void ang_to_xyz(double *x,double *y,double *z,int len);
	static void flip23(double *x,double *y,double *z,int *key,int k,int len);
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
	/**********************************************************/
	/* ######### STRIDPACK USED COMMANDS FOR VORONOI #########################*/

    /*  Various operation on images */
	/* out = img + scalar * img1  */
	static EMData* madn_scalar(EMData* img, EMData* img1, float scalar);
	/* out = scalar * img  */
	static EMData* mult_scalar(EMData* img, float scalar);
	/* out = img + img1  */
	static EMData* addn_img(EMData* img, EMData* img1);
	/* out = img - img1  */
	static EMData* subn_img(EMData* img, EMData* img1);
	/* out = img * img1  */
	static EMData* muln_img(EMData* img, EMData* img1);
	/* out = img / img1  */
	static EMData* divn_img(EMData* img, EMData* img1);
	/* out = |img|^2  */
	static EMData* squaren_img(EMData* img);
	/* img /= Re(img1) with zero check  */
	static EMData* divn_filter(EMData* img, EMData* img1);

	/* img += scalar * img1 */
	static void mad_scalar(EMData* img, EMData* img1, float scalar);
	/* img *= scalar  */
	static void mul_scalar(EMData* img, float scalar);
	/* img += img1  */
	static void add_img(EMData* img, EMData* img1);
	/* img += abs(img1)  */
	static void add_img_abs(EMData* img, EMData* img1);
	/* img += img1**2  */
	static void add_img2(EMData* img, EMData* img1);
	/* img -= img1  */
	static void sub_img(EMData* img, EMData* img1);
	/* img *= img1  */
	static void mul_img(EMData* img, EMData* img1);
	/* img /= img1  */
	static void div_img(EMData* img, EMData* img1);
	/* img = |img|^2  */
	static void square_img(EMData* img);
	/* img /= Re(img1) with zero check  */
	static void div_filter(EMData* img, EMData* img1);

	//utility for sxlocres
	static void set_freq(EMData* freqvol, EMData* temp, EMData* mask, float cutoff, float freq);


	/* pack absolute values of complex image into  real image with addition of Friedel part  */
	static EMData* pack_complex_to_real(EMData* img);
private:
	static float ang_n(float peakp, string mode, int maxrin); //this function is used by apmq()
public:


	/** formerly known as apmq
	 * Determine shift and rotation between image and many reference
	 * images (crefim, weights have to be applied) quadratic
	 * interpolation
	 * */
	static vector<float> multiref_polar_ali_2d(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, string mode,
                vector< int >numr, float cnx, float cny);

	/* In this version, we return a list of peaks for all reference images */
	static vector<float> multiref_polar_ali_2d_peaklist(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, string mode,
                vector< int >numr, float cnx, float cny);

	/* In this version, we return a list of peaks for local subset of reference images */
	static vector<float> multiref_polar_ali_2d_peaklist_local(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, float ant, string mode,
                vector< int >numr, float cnx, float cny);

	/* This is used in ISAC program to assigning particles equally to grops */
	static vector<int> assign_groups(std::string matrix_address, int nref, int nima);

	static inline void getvec(float phi, float theta, float& x, float& y, float& z, int option=0) {
		float pi180 = M_PI/180.0f;
		
		if (theta > 180.0f) {
			theta -= 180.0f;
			phi += 180.0f;
		} else if (theta > 90.0f && option == 0) {
			theta = 180.0f - theta;
			phi += 180.0f;
		}

		phi   *= pi180;
		theta *= pi180;

		x = sin(theta)*cos(phi);
		y = sin(theta)*sin(phi);
		z = cos(theta);

		return;
	}
	
	static inline float ang_diff(float v11, float v12, float v13, float v21, float v22, float v23, int& mirror) {
		float v = v11*v21+v12*v22+v13*v23;
		if (v > 1) v = 1;
		if (v < -1) v = -1;
		if (v > 0) { mirror = 1; return acos(v)*180.0f/M_PI; }
		else { mirror = -1; return acos(-v)*180.0f/M_PI; }
	}
	
	/* Find nearest projection angles
		Notice: the input I use is different from python code, which I think is awkward.
	*/
	static int nearest_ang(const vector<float>& vecref, float x, float y, float z);

	/* Assign projection angles to nearest reference projections */
	static vector<int> assign_projangles(const vector<float>& projangles, const vector<float>& refangles); 

	/* Assign howmany projection angles to the nearest reference projection */
	static vector<int> nearestk_to_refdir(const vector<float>& projangles, const vector<float>& refangles, const int howmany); 

	/* Group projection angles by (phi, theta) */
	static vector<int> group_proj_by_phitheta(const vector<float>& projangles, const vector<float>& ref_ang, const int img_per_grp);

	/** formerly known as apmq
	 * Determine shift and rotation between image and many reference
	 * images (crefim, weights have to be applied) quadratic
	 * interpolation
	 * */
	static vector<float> multiref_polar_ali_2d_delta(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, string mode,
                vector< int >numr, float cnx, float cny, float delta_start, float delta);

	/** formerly known as apnq DO NOT CONSIDER MIRROR
	 * Determine shift and rotation between image and many reference
	 * images (crefim, weights have to be applied) quadratic
	 * interpolation
	 * */
	static vector<float> multiref_polar_ali_2d_nom(EMData* image, const vector< EMData* >& crefim,
                float xrng, float yrng, float step, string mode,
                vector< int >numr, float cnx, float cny);

	/** formerly known as apmq
	 * Determine shift and rotation between image and many reference
	 * images (crefim, weights have to be applied) quadratic
	 * interpolation
	 * */
	static vector<float> multiref_polar_ali_2d_local(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, float ant, string mode,
                vector< int >numr, float cnx, float cny, string sym);

	/* Returns first match with peak greater than previousmax or the best match in whole space (when there are no peaks > previousmax).
	 * The reference rings are checked in random order.
	 * */
	static vector<float> shc(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, float ant, string mode,
                vector< int >numr, float cnx, float cny, string sym);
	static vector<float> shc_multipeaks(EMData* image, const vector< EMData* >& crefim,
	            vector<float> xrng, vector<float> yrng, float step, float ant, string mode,
	            vector<int>numr, float cnx, float cny, int max_peaks_count);

	/** formerly known as apmq
	 * Determine shift and rotation between image and many reference
	 * images (crefim, weights have to be applied) quadratic
	 * interpolation
	 * Search for peaks only within +/-psi_max from 0 and 180 (helical)
	 * */
	static vector<float> multiref_polar_ali_helical(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, float psi_max, string mode,
                vector< int >numr, float cnx, float cny, int ynumber=-1);
	static vector<float> multiref_polar_ali_helical_local(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, float ant, float psi_max, string mode,
                vector< int >numr, float cnx, float cny, int ynumber=-1, float yrnglocal=-1.0);
	static vector<float> multiref_polar_ali_helical_90(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, float psi_max, string mode,
                vector< int >numr, float cnx, float cny, int ynumber=-1);
	static vector<float> multiref_polar_ali_helical_90_local(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, float ant, float psi_max, string mode,
                vector< int >numr, float cnx, float cny, int ynumber=-1, float yrnglocal=-1.0);
    /**  Next two for helicon **/
	static vector<float> multiref_polar_ali_helicon_local(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, float ant, float psi_max, string mode,
                vector< int >numr, float cnx, float cny, int ynumber=-1, float yrnglocal=-1.0);
	static vector<float> multiref_polar_ali_helicon_90_local(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, float ant, float psi_max, string mode,
                vector< int >numr, float cnx, float cny, int ynumber=-1, float yrnglocal=-1.0);

	/** formerly known as apmq
	 * Determine shift and rotation between image and many reference
	 * images (crefim, weights have to be applied) quadratic
	 * interpolation
	 * */
	static vector<float> multiref_polar_ali_2d_local_psi(EMData* image, const vector< EMData* >& crefim,
                vector<float> xrng, vector<float> yrng, float step, float ant, float psi_max, string mode,
                vector< int >numr, float cnx, float cny);
	/** Determine shift and rotation between image and one reference
	 * image (crefim, weights have to be applied) using quadratic
	 * interpolation, return a list of peaks  PAP  07/21/08
	 *
	 * ccf1d keeps 1d ccfs stored as (maxrin, -kx-1:kx+1, -ky-1:ky+1)
	 * margin is needed for peak search and both arrays are initialized with -1.0e20
	 * */
	static void multiref_peaks_ali2d(EMData* image, EMData* crefim,
                vector<float> xrng, vector<float> yrng, float step, string mode,
                vector< int >numr, float cnx, float cny, EMData* peaks, EMData* peakm);

	/** Determine shift and rotation between image and one reference
	 * image (crefim, weights have to be applied) using quadratic
	 * interpolation, return a list of peaks  PAP  07/21/08
	 *
	 * ccf1d keeps 1d ccfs stored as (maxrin, -kx-1:kx+1, -ky-1:ky+1)
	 * margin is needed for peak search and both arrays are initialized with -1.0e20
	 * */
	static void multiref_peaks_compress_ali2d(EMData* image, EMData* crefim, vector<float> xrng, vector<float> yrng,
	     float step, string mode, vector<int>numr, float cnx, float cny, EMData *peaks, EMData *peakm,
	     EMData *peaks_compress, EMData *peakm_compress);

	/** Determine shift and rotation between image and one reference
	 * image (crefim, weights have to be applied) using quadratic
	 * interpolation
	 * */
	static vector<float> ali2d_ccf_list(EMData* image, EMData* crefim, vector<float> xrng, vector<float> yrng,
	     float step, string mode, vector<int>numr, float cnx, float cny, double T);

	static vector<float> ali2d_ccf_list_snake(EMData* image, EMData* crefim, vector<float> wr, float xrng, float yrng,
	     float step, string mode, vector<int>numr, float cnx, float cny, double T);			     
     /*
	static void multiref_peaks_ali(EMData* image, const vector< EMData* >& crefim,
                float xrng, float yrng, float step, string mode,
                vector< int >numr, float cnx, float cny, EMData* peaks, EMData* peakm,
		    int nphi, int ntheta);
*/
	static vector<float> twoD_fine_ali(EMData* image, EMData *refim, EMData* mask, float ang, float sxs, float sys);

	static vector<float> twoD_fine_ali_G(EMData* image, EMData *refim, EMData* mask, Util::KaiserBessel& kb, float ang, float sxs, float sys);

	static vector<float> twoD_to_3D_ali(EMData* volft, Util::KaiserBessel& kb, EMData *refim, EMData* mask, float phi, float theta, float psi, float sxs, float sxy);

	static vector<float> twoD_fine_ali_SD(EMData* image, EMData *refim, EMData* mask, float ang, float sxs, float sys);

	static float ccc_images(EMData *, EMData *, EMData *, float , float , float );

	static vector<float> twoD_fine_ali_SD_G(EMData* image, EMData *refim, EMData* mask, Util::KaiserBessel& kb, float ang, float sxs, float sys);

	static float ccc_images_G(EMData* image, EMData* refim, EMData* mask, Util::KaiserBessel& kb, float ang, float sx, float sy);

	static float local_inner_product(EMData* image1, EMData* image2, int lx, int ly, int lz, int w);

	static EMData* move_points(EMData* img,  float qprob, int ri, int ro);

	static EMData* get_biggest_cluster( EMData* mg );

	//static EMData* ctf_img(int nx, int ny, int nz, float dz, float ps, float voltage=300.0f,float cs=2.0f,float wgh=0.1f,float b_factor=0.0f,float dza=0.0f,float azz=0.0f,float sign=-1.0f);
	static EMData* ctf_img(int nx, int ny, int nz, float dz, float ps, float voltage,float cs,float wgh,float b_factor,float dza,float azz,float sign);
	static EMData* ctf_rimg(int nx, int ny, int nz, float dz, float ps, float voltage,float cs,float wgh,float b_factor,float dza,float azz,float sign);
	static EMData* ctf2_rimg(int nx, int ny, int nz, float dz, float ps, float voltage,float cs,float wgh,float b_factor,float dza,float azz,float sign);

	static inline int mono(int k1, int k2) {
#ifdef _WIN32
		int  mk = _cpp_max(k1,k2);
		return  _cpp_min(k1,k2) + mk*(mk-1)/2;
#else
		int  mk = std::max(k1,k2);
		return  std::min(k1,k2) + mk*(mk-1)/2;
#endif	//_WIN32
	}

        static inline int nint180(float arg) {
	    int res = int(arg + 180.5) - 180;
	    return res;
        }
	
	static inline double mean(double *x, int n) {
		double s = 0.0;
		for (int i=0; i<n; i++) s+=x[i];
		return s/static_cast<double>(n);
	}

	static inline double var(double *x, int n) {
		double s = 0.0;
		double m = mean(x, n);
		for (int i=0; i<n; i++) s += (x[i]-m)*(x[i]-m);
		return s/static_cast<double>(n);
	}
	
	static vector<float> multi_align_error(vector<float> args, vector<float> all_ali_params, int d);
	static double multi_align_error_func(double* x, vector<float> all_ali_params, int nima, int num_ali, int d);
	static vector<double> multi_align_error_func2(double* x, vector<float> all_ali_params, int nima, int num_ali, int d);
	static void multi_align_error_dfunc(double* x, vector<float> all_ali_params, int nima, int num_ali, double* g, int d);
	
	static vector<float> cluster_pairwise(EMData* d, int K, float T, float F);
	//static vector<float> cluster_equalsize(EMData* d, int m);
	static vector<float> cluster_equalsize(EMData* d);
	static vector<float> vareas(EMData* d);

	/**
	 * This function returns a 2-D slice from a 3-D EMData object
	 * dim denotes the slice is perpendicular to which dimension
	 * 1 for x-dimension, 2 for y-dimension and 3 for z-dimension
	*/
	static EMData* get_slice(EMData *vol, int dim, int index);

	static void image_mutation(EMData *img, float mutation_rate);

	/** The purpose of this function is to convert a list to grey code and mutate them and convert them back */
	static void array_mutation(float* list, int len_list, float mutation_rate, float min_val, float max_val, int K, int is_mirror);

	static vector<float> list_mutation(vector<float> list, float mutation_rate, float min_val, float max_val, int K, int is_mirror);
	/*
			To restrict the value to [0, nx)
	*/
	static inline float restrict1(float x, int nx) {
		while ( x < 0.0f )        x += nx;
		while ( x >= (float)(nx) )  x -= nx;
		return x;
	}

	/** This function returns parameters from Transform object as a Dict object.
	 *  It allows to omit problem with passing Transform object through C++/Python layer (It causes memory leak). */
	static Dict get_transform_params(EMData* image, string xform, string convention);

	static void constrained_helix_exhaustive( vector<EMData*> data, vector<EMData*> fdata, vector<EMData*> refproj, vector<EMData*> rotproj
			, vector<float> dp_dphi_rise_delta, vector<int> nphi_phiwobble_range_ywobble_Dsym_nwx_nwy_nwxc_nwyc
			, bool FindPsi, float psi_max, vector<EMData*> crefim, vector<int> numr, int maxrin, string mode, int cnx, int cny);
	
	static std::vector<float> diff_between_matrix_of_3D_parameters_angles( std::vector<float> all_params, std::vector<float> rotations );

	/** Calculates max_clique in undirected graph, input: edges coded as pair of integers (integers correspond to vertices, must be >=0)
	 *  Returns: vertices (integers) belonging to max clique
	 * */
	static std::vector<int> max_clique(std::vector<int> edges);

#endif	//util__sparx_h__
