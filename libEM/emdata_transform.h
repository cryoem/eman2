/**
 * $Id$
 */
#ifndef emdata__transform_h__
#define emdata__transform_h__ 

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


/** return the fast fourier transform (FFT) image of the current
 * image. the current image is not changed. The result is in
 * real/imaginary format.
 * @return The FFT of the current image in real/imaginary format.
 */
EMData *do_fft();


/** Do FFT inplace. And return the FFT image.
 * @return The FFT of the current image in real/imaginary format.
 */
EMData* do_fft_inplace();


/** return the inverse fourier transform (IFT) image of the current
 * image. the current image is not changed.
 *
 * @exception ImageFormatException If the image is not a complex image.
 * @return The current image's inverse fourier transform image.
 */
EMData *do_ift();


/* Do IFT inplace. And return the IFT image.
 * @return The IFT image.
 */
EMData* do_ift_inplace();


/** Render the image into an 8-bit image. 2D image only.
 *
 * @param x	origin of the area to render
 * @param y
 * @param xsize	size of the area to render in output pixels
 * @param ysize
 * @param bpl	bytes per line, if asrgb remember *3
 * @param scale	scale factor for rendering
 * @param min_gray	minimum gray value to render (0-255)
 * @param max_gray	maximum gray value to render (0-255)
 * @param min_render	float image density corresponding to min_gray
 * @param max_render	float image density corresponding to max_gray
 * @param asrgb	duplicate each output pixel 3x for RGB rendering
 * @exception ImageDimensionException If the image is not 2D.
 */
std::string render_amp8(int x, int y, int xsize, int ysize,
				 int bpl, float scale, int min_gray, int max_gray,
				 float min_render, float max_render,int asrgb);

		
/** Render the image into a 24-bit image. 2D image only.
 * @param x
 * @param y
 * @param xsize
 * @param ysize
 * @param bpl
 * @param scale
 * @param min_gray
 * @param max_gray
 * @param min_render
 * @param max_render
 * @param ref
 * @param cmap
 * @exception ImageDimensionException If the image is not 2D.
 */		 
void render_amp24(int x, int y, int xsize, int ysize,
				  int bpl, float scale, int min_gray, int max_gray,
				  float min_render, float max_render,
				  void *ref, void cmap(void *, int coord, unsigned char *tri));


/** convert the complex image from real/imaginary to amplitude/phase */
void ri2ap();


/** convert the complex image from amplitude/phase to real/imaginary */
void ap2ri();


/** Insert a clip into this image.
 * @param block An image block.
 * @param origin The origin location to insert the clip.
 * @exception ImageFormatException If clip is outside of the
 * destination image (i.e., this image).
 */
void insert_clip(EMData * block, const IntPoint & origin);


/** Add a scaled image into another image at a specified location.
 *  This is used, for example, to accumulate gaussians in
 *  programs like pdb2mrc.py. The center of 'block' will be positioned at
 *  'center' with scale factor 'scale'. Densities will be interpolated in
 *  'block' and multiplied by 'mult'.
 *
 * @param block The image to inserted.
 * @param center The center of the inserted block in 'this'.
 * @param scale  Scale factor.
 * @param mult_factor Number used to multiply the block's densities.
 * @exception ImageDimensionException If 'this' image is not 2D/3D.
 */
void insert_scaled_sum(EMData *block, const FloatPoint & center,
					   float scale=1.0, float mult_factor=1.0);

#endif	//emdata__transform_h__
