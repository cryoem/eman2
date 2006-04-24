/**
 * $Id$
 */
 
/** This file is a part of util.h, 
 * To use this file's functions, you should #include "util.h"
 * NEVER directly include this file */
 
#ifndef util__sparx_h__
#define util__sparx_h__ 

public:

/** Creates a Two D Test Pattern
 * @param[in] Size, must be odd
 * @param[in] p the x frequency
 * @param[in] q the y frequency
 * @param[in] a the x falloff
 * @param[in] alpha the projection angle
 * @param[out] The 2D test pattern in real space, fourier space,
 *               or the projection in real or fourier space 
 *               or the FH of the pattern
*/



static vector<float> infomask(EMData* Vol, EMData* mask);

static void colreverse(float* beg, float* end, int nx);

static void slicereverse(float* beg, float* end, int nx,int ny);

static  void cyclicshift(EMData* image, Dict params);

static Dict im_diff(EMData* V1, EMData* V2, EMData* mask);

static EMData* TwoDTestFunc(int Size, float p, float q,  float a, float b, 
                   int flag=0, float alphaDeg=0); //PRB


/** Given a tabulated function y of x (n unordered points), and
 * Given the values of the m values xq to be interpolated
 * This routine returns the interpolated array yq, PRB
 * This function is called by splint 
 * @param[in] y of x is the tabulated function of length n
 * @param[in] xq is the x values to be splined: has m points.
 * @param[out]  yq are the splined values
*/
static void spline_mat(float *x, float *y, int n,  float *xq, float *yq, int m); //PRB

/** Given a tabulated function y of x (unordered), and
 * Given the values of the first derivatives at the end points
 * This routine returns an array y2, that contains the second derivatives
 * of the function at the tabulated points.   PRB
 * This function is called by splint 
 * @param[in] y of x is the tabulated function of length n
 * @param[in] yp1, ypn : the derivatives of the first and last point.
 * @param[out]  y2 is the value of the second derivatives
*/
static void spline(float *x, float *y, int n, float yp1, float ypn, float *y2);

/** Given the arrays xa(ordered, ya of length n, which tabulate a function 
 *  and given the array y2a which is the output of spline and an unordered array xq,
 *  this routine returns a cubic-spline interpolated array yq.
 * @param[in] y of x is the tabulated function of length n
 * @param[in] y2a is returned from spline: second derivs
 * @param[in] xq is the x values to be splined: has m points.
 * @param[out]  yq are the splined values
*/  // PRB
static void splint( float *xa, float *ya, float *y2a, int n,  
                                     float *xq, float *yq, int m);


/** list the sorted lengths of the  integer lattice 
 * sites of a square sided image of  size Size. PRB
 * @param[in] Size The length of the image 
 * @param[out] PermMatTr  The matrix telling the ordering of the lattice 
 *                        sites wrt the array
 * @param[out] weightofkvalsSorted the number of sites at that distance
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
				float i0win_tab(float x) const {
					float absx = fabs(x);
					int loc = int(round(absx*fltb));
					return i0table[loc];
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

		/** Compute a vector containing quais-evenly spaced Euler angles
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
		voea(float delta, float t1=0, float t2=90, 
			 float p1=0, float p2=359.9);

		/** Tri-Quadratic interpolation.
		 *
		 *	@param[in] r x-coord value
		 *	@param[in] s y-coord value
		 *	@param[in] t z-coord value
		 *	@param[in] f 3x3x3 grid of measured values
		 *
		 *	@return Interpolated value
		 */
		static float triquad(double r, double s, double t, float f[]);

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
		 *	@param[in] image Image object (pointer)
		 *
		 *	@return Interpolated value
		 */
		static float quadri(float x, float y, int nx, int ny, float* image);

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
       static EMData* alrq_msi(EMData* image,float cns2, float cnr2,
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
	
	
	/*static EMData* window(EMData* img,int new_nx,int new_ny=1, int new_nz=1, int shift_x=1, int shift_y=1, int
	shift_z=1);*/
	static EMData* window(EMData* img,int new_nx ,int new_ny=1, int new_nz=1, int x_offset=0, int y_offset=0, int z_offset=0);
	
	
	static EMData* pad(EMData* img, int new_nx, int new_ny=1, int new_nz=1, int x_offset=0, int y_offset=0, int z_offset=0,int background=0);
	
	
#endif	//util__sparx_h__
