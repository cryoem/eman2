/**
 * $Id$
 */
#ifndef normpadft_h__
#define normpadft_h__ 

public:

/** Multiply a real image by (-1)**(ix+iy+iz) to center
 *  the fft version.
 *
 */
void center_origin();


/** Multiply a Fourier image by (-1)**(ix+iy+iz) to center it.
 *
 */
void center_origin_fft();


/** return an image object that has been padded with zeros to
 *  an integral (npad) factor (e.g., npad=2 means the new 
 *  image will be 2x larger in each direction).  The default
 *  is to pad 4x.
 * The current image is not changed.
 *
 * @return An image object that has been padded npad-times.
 */
EMData *zeropad_ntimes(int npad=4);


/** return an image object that has been fft-padded/unpadded.
 * The current image is not changed.
 *
 * @return An image object that has been fft-padded/unpadded.
 */
EMData *pad_fft(int npad = 1);


/** Remove padding, leaving a single corner of the image.
 *  The current image is changed in place.
 *
 *  The assumption is that after an in-place inverse fft
 *  the real-space image contains too much information 
 *  because it may have been zero-padded some integer factor
 *  of times and it has also been extended slightly along x
 *  for the fft.  Here we keep only the data corresponding 
 *  to ix=0,...,nxold-1, iy=0,...,nyold-1, iz=0,...,nzold-1,
 *  where nxold, nyold, nzold are the sizes of the original
 *  image.
 */
void postift_depad_corner_inplace();


/** Interpolate up image by padding with zeroes its Fourier transfrom.
 *  
 *  @param[in] nxni new x size (has to be larger/equal than the original x size)
 *  @param[in] nyni new y size (has to be larger/equal than the original y size)
 *  @param[in] nzni new z size (has to be larger/equal than the original z size)
 *  
 *  @return New interpolated up image.
 */
EMData* FourInterpol(int nxni, int nyni=0, int nzni=0);




#endif	//normpadft_h__
