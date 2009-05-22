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

static Dict coveig_for_py(int ncov, const vector<float>& covmatpy);

static Dict ExpMinus4YSqr(float ymax,int nsamples);

static void WTM(EMData* PROJ, vector<float> SS,int DIAMETER,int NUMP);

static void WTF(EMData* PROJ,vector<float> SS,float SNR,int K,vector<float> exptable);

static Dict CANG(float PHI, float THETA, float PSI);

static void BPCQ(EMData* B, EMData *CUBE,vector<float> DM);

static vector<float> infomask(EMData* Vol, EMData* mask, bool);

static void colreverse(float* beg, float* end, int nx);

static void slicereverse(float* beg, float* end, int nx,int ny);

static void cyclicshift(EMData* image, Dict params);

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

		static float get_pixel_conv_new(int nx, int ny, int nz, float delx, float dely, float delz, float* data, Util::KaiserBessel& kb);

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
	static void  Frngs(EMData* circ, vector<int> numr);
	static void  Frngs_inv(EMData* circ, vector<int> numr);

	/*
	  	A little notes about different Crosrng:
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
	static Dict Crosrng_ew(EMData* circ1, EMData* circ2, vector<int> numr, vector<float> w, int neg);

	static Dict Crosrng_ms(EMData* circ1, EMData* circ2, vector<int> numr);
	static Dict Crosrng_ns(EMData* circ1, EMData* circ2, vector<int> numr);

	static EMData* Crosrng_msg(EMData* circ1, EMData* circ2, vector<int> numr);
	static void Crosrng_msg_vec(EMData* circ1, EMData* circ2, vector<int> numr, float *q, float *t);
	static EMData* Crosrng_msg_s(EMData* circ1, EMData* circ2, vector<int> numr);
	static EMData* Crosrng_msg_m(EMData* circ1, EMData* circ2, vector<int> numr);

	static vector<float> Crosrng_msg_vec_p(EMData* circ1, EMData* circ2, vector<int> numr );
	static void  prb1d(double *b, int npoint, float *pos);

	static void update_fav(EMData* ave,EMData* dat, float tot, int mirror, vector<int> numr);
	static void sub_fav(EMData* ave,EMData* dat, float tot, int mirror, vector<int> numr);

	static float ener(EMData* ave, vector<int> numr);

	static float ener_tot(const vector<EMData*>& data, vector<int> numr, vector<float> tot);

        // k-means helper
        static Dict min_dist_real(EMData* image, const vector<EMData*>& data);
        static Dict min_dist_four(EMData* image, const vector<EMData*>& data);

        // temp comment before removing 2009-03-31 09:35:20 JB
        //static float SqEuc_dist(EMData* image, EMData* width);

        static vector<float> cml_line_in3d_full(const vector<float>& Ori);
        static vector<double> cml_line_in3d_iagl(const vector<float>& Ori, float phi, float theta, int iprj);
        static vector<double> cml_weights(const vector<float>& cml);
        static vector<float> cml_spin(int n_psi, int i_prj, int n_prj, vector<float> weights, vector<int> com, const vector<EMData*>& data, int flag);
        static vector<int> cml_line_pos(float phi1, float theta1, float psi1, float phi2, float theta2, float psi2, int nangle);
        static vector<int> cml_list_line_pos(vector<float> Ori, float newphi, float newtheta, int i_prj, int n_prj, int nangle, int nlines);

        // new code common-lines
        static vector<int> cml_line_insino(vector<float> Rot, int i_prj, int n_prj);
        static vector<int> cml_line_insino_all(vector<float> Rot, vector<int> seq, int n_prj, int n_lines);
        static vector<double> cml_init_rot(vector<float> Ori);
        static vector<float> cml_update_rot(vector<float> Rot, int iprj, float nph, float th, float nps);
        static vector<double> cml_line_in3d(vector<float> Ori, vector<int> seq, int nprj, int nlines);
        static vector<double> cml_spin_psi(const vector<EMData*>& data, vector<int> com, vector<float> weights, int iprj, vector<int> iw, int n_psi, int d_psi, int n_prj); 
        static double cml_disc(const vector<EMData*>& data, vector<int> com, vector<int> seq, vector<float> weights, int n_lines);
static void set_line(EMData* img, int posline, EMData* line);

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

	static vector<float> histogram(EMData* image, EMData* mask, int nbins = 128, float hmin =0.0f, float hmax = 0.0f );

	static Dict histc(EMData *ref,EMData *img,EMData *mask);

	static float hist_comp_freq(float PA,float PB,int size_img, int hist_len, EMData *img, vector<float> ref_freq_hist, EMData *mask, float ref_h_diff, float ref_h_min);


	/* The unit in the ctf function: dz: Angstrom, cs: CM  Ps: Angstrom, Voltage: Kv,dza: Angstrom, azz: degree wgh: None unit. b_factor: Angstrom^2
	 The CTF function takes form of   *sin(-quadpi*(dz*lambda*ak^2-cs*lambda^3*ak^4/2.)-wgh)*exp(-b_factor*ak^2)*sign
          * sign can be set as +1 or -1 . The unit of frequency ak is 1/Angstrom
                  Attention: Envelope function in power spectrum has a form of exp(-b_factor*ak^2)
                                          */
	static float   tf(float dzz, float ak, float voltage = 300.0f, float cs = 2.0f, float wgh = 0.1f, float b_factor = 0.0f, float sign = -1.0f);
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
	/* img /= Re(img1) with zero check  */
	static EMData* divn_filter(EMData* img, EMData* img1);

	/* img += scalar * img1 */
	static void mad_scalar(EMData* img, EMData* img1, float scalar);
	/* img *= scalar  */
	static void mul_scalar(EMData* img, float scalar);
	/* img += img1  */
	static void add_img(EMData* img, EMData* img1);
	/* img += img1**2  */
	static void add_img2(EMData* img, EMData* img1);
	/* img -= img1  */
	static void sub_img(EMData* img, EMData* img1);
	/* img *= img1  */
	static void mul_img(EMData* img, EMData* img1);
	/* img /= img1  */
	static void div_img(EMData* img, EMData* img1);
	/* img /= Re(img1) with zero check  */
	static void div_filter(EMData* img, EMData* img1);
	/* pack absolute values of complex image into  real image with addition of Friedel part  */
	static EMData* pack_complex_to_real(EMData* img);
private:
	static float ang_n(float peakp, string mode, int maxrin); //this function is used by apmq()
public:
	static vector<float> multiref_polar_ali_2d(EMData* image, const vector< EMData* >& crefim,
                float xrng, float yrng, float step, string mode,
                vector< int >numr, float cnx, float cny);
	static vector<float> multiref_polar_ali_2d_nom(EMData* image, const vector< EMData* >& crefim,
                float xrng, float yrng, float step, string mode,
                vector< int >numr, float cnx, float cny);
	static vector<float> multiref_polar_ali_2d_local(EMData* image, const vector< EMData* >& crefim,
                float xrng, float yrng, float step, float ant, string mode,
                vector< int >numr, float cnx, float cny);

	static void multiref_peaks_ali2d(EMData* image, EMData* crefim,
                float xrng, float yrng, float step, string mode,
                vector< int >numr, float cnx, float cny, EMData* peaks, EMData* peakm);
	
	static void multiref_peaks_compress_ali2d(EMData* image, EMData* crefim, float xrng, float yrng, 
	     float step, string mode, vector<int>numr, float cnx, float cny, EMData *peaks, EMData *peakm, 
	     EMData *peaks_compress, EMData *peakm_compress);

	static vector<float> ali2d_ccf_list(EMData* image, EMData* crefim, float xrng, float yrng, 
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

	static EMData* move_points(EMData* img,  float qprob, int ri, int ro);

	static EMData* get_biggest_cluster( EMData* mg );

	//static EMData* ctf_img(int nx, int ny, int nz, float dz, float ps, float voltage=300.0f,float cs=2.0f,float wgh=0.1f,float b_factor=0.0f,float dza=0.0f,float azz=0.0f,float sign=-1.0f);
	static EMData* ctf_img(int nx, int ny, int nz, float dz, float ps, float voltage,float cs,float wgh,float b_factor,float dza,float azz,float sign);

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
	    float res = int(arg + 180.5) - 180.0;
	    return res;
        }

	static vector<float> cluster_pairwise(EMData* d, int K, float T, float F);
	//static vector<float> cluster_equalsize(EMData* d, int m);
	static vector<float> cluster_equalsize(EMData* d);
	static vector<float> vareas(EMData* d);
	static EMData* get_slice(EMData *vol, int dim, int index);
	/*
			To restrict the value to [0, nx)
	*/
	static inline float restrict1(float x, int nx) {
		while ( x < 0.0f )        x += nx;
		while ( x >= (float)(nx) )  x -= nx;
		return x;
	}

#endif	//util__sparx_h__
