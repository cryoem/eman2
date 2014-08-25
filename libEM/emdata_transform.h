/**
 * $Id$
 */

/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
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
 * */

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
EMData *do_fft() const;


#ifdef EMAN2_USING_CUDA
/** return the fast fourier transform (FFT) image of the current
 * image. the current image is not changed. The result is in
 * real/imaginary format and exists only on the GPU
 * @return The FFT of the current image in real/imaginary format, existing on the GPU.
 * @ingroup CUDA_ENABLED
 */
EMData *do_fft_cuda(); // I am not using const b/c this function does modify cuda specific EMData members...

/** return the inplace fast fourier transform (FFT) image of the current
 * image. the current image is not changed. The result is in
 * real/imaginary format and exists only on the GPU
 * @return The FFT of the current image in real/imaginary format, existing on the GPU.
 * @ingroup CUDA_ENABLED
 */
void do_fft_inplace_cuda();

/**  return the inverse fourier transform (IFT) image of the current
 * image. The result exists only on the GPU
 * @param preserve_input whether or not this EMData object should be preserved. If this is unecessary than we can avoid a copy and run faster
 * @return The FFT of the current image in real/imaginary format, existing on the GPU.
 * @ingroup CUDA_ENABLED
 */
EMData *do_ift_cuda();

/**  return the inverse fourier transform (IFT) image of the current
 * image inplace. The result exists only on the GPU
 * @param preserve_input whether or not this EMData object should be preserved. If this is unecessary than we can avoid a copy and run faster
 * @return The FFT of the current image in real/imaginary format, existing on the GPU.
 * @ingroup CUDA_ENABLED
 */
void do_ift_inplace_cuda();
#endif //EMAN2_USING_CUDA



/** Do FFT inplace. And return the FFT image.
 * @return The FFT of the current image in real/imaginary format.
 */
void do_fft_inplace();


/** return the inverse fourier transform (IFT) image of the current
 * image. the current image may be changed if it is in amplitude/phase
 * format as opposed to real/imaginary format - if this change is
 * performed it is not undone.
 *
 * @exception ImageFormatException If the image is not a complex image.
 * @return The current image's inverse fourier transform image.
 */
EMData *do_ift();


/* Do IFT inplace. And return the IFT image.
 * @return The IFT image.
 */
void do_ift_inplace();


/** Render the image into an 8-bit image. 2D images only.
 * flags provide a way to do unusual things with this function, such
 * as calculating a histogram of the rendered area.
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
 * @param gamma
 * @param flags	1-duplicate each output pixel 3x for RGB rendering,2-add a 256 int greyscale histogram to the end of the image array,4-invert y axis,8-render 32 bit 0xffRRGGBB
 * @exception ImageDimensionException If the image is not 2D.
 */
std::string render_amp8(int x, int y, int xsize, int ysize,
				 int bpl, float scale, int min_gray, int max_gray,
				 float min_render, float max_render,float gamma,int flags);

/** Render the image into an 8-bit image. 2D images only.
 * flags provide a way to do unusual things with this function, such
 * as calculating a histogram of the rendered area.
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
 * @param gamma
 * @param flags	1-duplicate each output pixel 3x for RGB rendering,2-add a 256 int greyscale histogram to the end of the image array,4-invert y axis,8-render 32 bit 0xffRRGGBB
 * @exception ImageDimensionException If the image is not 2D.
 */
std::string render_ap24(int x, int y, int xsize, int ysize,
				 int bpl, float scale, int min_gray, int max_gray,
				 float min_render, float max_render,float gamma,int flags);

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


/** convert the complex image from real/imaginary to amplitude/phase
*/
void ri2ap();

/** convert the complex image from amplitude/phase to real/imaginary
 * */
void ap2ri();

/** convert the complex image from real/imaginary to Intensity/0.
This conversion cannot be reversed, and the image remains marked as R/I
*/
void ri2inten();

/**   This computes the rotational and translational bispectral
        invariants of an image. The invariants are labelled by the Fourier
    Harmonic label given by N.
     fVec is the real input image
      NK is the number of Fourier components one wishes to use in calculating this bispectrum
   the output is a single 2D image whose x,y labels are lengths, corresponding to the two lengths of sides of
         a triangle  */
EMData*   bispecRotTransInvN(int N, int NK);



/**   This computes the rotational and translational bispectral
        invariants of an image.
   the output is a single 3d Volume whose x,y labels are lengths, corresponding to the two lengths of sides of
         a triangle
         the z label is for the angle  */
EMData*  bispecRotTransInvDirect(int type=0);


/** Insert a clip into this image.
 * Very robust clip insertion code works in all way you might think possible
 * @param block An image block.
 * @param origin The origin location to insert the clip.
 * @ingroup CUDA_ENABLED (but slow)
 */
void insert_clip(const EMData * const block, const IntPoint & origin);


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
