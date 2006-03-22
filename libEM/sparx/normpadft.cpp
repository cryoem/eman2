#include "emdata.h"

using namespace EMAN;
using std::vector;

namespace EMAN {
	/* #G2#
	Purpose: Create a new [normalized] [zero-padded] fft image. 
	Method: Normalize, pad with zeros, extend for fft, create new fft image,
	return new fft image.  
	Input: f real n-dimensional image
	flag specify normalize, pad, and/or extend 
	Output: fft image of normalized, zero-padded input image
	 */
	EMData* norm_pad_ft(EMData* fimage, bool donorm, bool dopad, int npad) {
		int nx = fimage->get_xsize();
		int ny = fimage->get_ysize();
		int nz = fimage->get_zsize();
		float mean = 0., stddev = 1.;
		if(donorm) { // Normalization requested
			mean = fimage->get_attr("mean");
			stddev = fimage->get_attr("sigma");
		}
		//  Padding requested
		if (dopad) {
			// sanity check
			if (npad < 2) npad = 2;
		} else {
			npad = 1;
		}
		EMData& fpimage = *fimage->pad_fft(npad);
		//  Perform the actual normalization (only on the
		//  non-zero section of the image)
		if(donorm) { // Normalization requested
			// indexing starts at 1
			vector<int> saved_offsets = fpimage.get_array_offsets();
			fpimage.set_array_offsets(1,1,1);
			for (int iz = 1; iz <= nz; iz++) {
				for (int iy = 1; iy <= ny; iy++) {
					for (int ix = 1; ix <= nx; ix++) {
						fpimage(ix,iy,iz) = (fpimage(ix,iy,iz)-mean)/stddev;
					}
				}
			}
			fpimage.set_array_offsets(saved_offsets); // reset
		}
		mean = fpimage.get_attr("mean");
		stddev = fpimage.get_attr("sigma");
		fpimage.do_fft_inplace();
		return &fpimage;
	}

}

/* vim: set ts=4 noet: */
