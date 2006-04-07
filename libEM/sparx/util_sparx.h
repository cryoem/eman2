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

static EMData* im_diff(EMData* V1, EMData* V2, EMData* mask);
static vector<float> infomask(EMData* Vol, EMData* mask);


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

#endif	//util__sparx_h__
