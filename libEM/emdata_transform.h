/**
 * $Id$
 */

/** This file is a part of "emdata.h", to use functions in this file,
 * you should "#include "emdata.h",
 * NEVER directly include this file. */ 

#ifndef emdata__transform_h__
#define emdata__transform_h__ 

public:

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
 * @param flags	1-duplicate each output pixel 3x for RGB rendering,2-add a 256 character greyscale histogram to the end of the image array
 * @exception ImageDimensionException If the image is not 2D.
 */
std::string render_amp8(int x, int y, int xsize, int ysize,
				 int bpl, float scale, int min_gray, int max_gray,
				 float min_render, float max_render,int flags);

		
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
