/**
 * $Id$
 */
#ifndef eman__interp_h__
#define eman__interp_h__

#include <math.h>

namespace EMAN
{

	/** Interp defines the interpolation function used to generate
     * a e^-x^4 function in real space. Used for Fourier interpolation
     */
	class Interp
	{
	  public:
		static float get_hyperg(int i)
		{
			return HYPERGEOM[i];
		}

		static float hyperg(float v)
		{
			if (v < 0 || v > 4.998) {
				return 0;
			}
			float r = v / 0.001;
			int a = (int) floor(r);
			r -= a;
			return -(HYPERGEOM[a] * (1.0 - r) + HYPERGEOM[a + 1] * r);
		}

		static float *get_gimx()
		{
			if (!gimx) {
				init_gimx();
			}
			return gimx;
		}

	  private:
		static void init_gimx();

	  private:
		static float HYPERGEOM[];
		static float *gimx;
	};
}

#endif
