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

#ifndef emdata__sparx_h__
#define emdata__sparx_h__ 

public:
/** returns the fourier harmonic transform (FH) image of the current
 * image (in real space). The current image is not changed. The result is in
 * real/imaginary format. The FH switch is set on.
 * @param OverSamplekB is a parameter controlling the fineness of the Fourier sampling
 * @return the Fourier Harmonic image
 */
EMData *real2FH(float OverSamplekB);

/** returns the fourier version of the image 
 * from the FH version. The current image is not changed. The result is in
 * real/imaginary format. The FH switch is set off.
 * @param Size is the size of the image to be returned
 * @param OverSamplekB is a parameter controlling the fineness of the Fourier sampling
 * @param IntensityFlag =0 is the usual; =1 means that the input was an intensity
 * @return the shuffled version of the FFT
 */
EMData *FH2F(int Size, float OverSamplekB, int IntensityFlag =0);

/** returns the real version of the image 
 * from the FH version. The current image is not changed. The result is in
 * real format.
 * @param Size is the size of the image to be returned
 * @param OverSamplekB is a parameter controlling the fineness of the Fourier sampling
 * @param IntensityFlag =0 is the usual; =1 means that the input was an intensity
 * @return the real version of the data
 */
EMData *FH2Real(int Size, float OverSamplekB, int IntensityFlag =0);
		
		
/** Create a (1-D) rotationally averaged image.
 * @exception ImageDimensionException If 'this' image is 1D.
 * @return 1-D rotationally-averaged image
 */					
EMData* rotavg();

/** Create a 2-D or 3-D rotationally averaged image.
 * @exception ImageDimensionException If 'this' image is 1D.
 * @return 2-D or 3-D rotationally-averaged image
 */					
EMData* rotavg_i();

/** Calculates the Center of Gravity 
 *  and the Radius of Gyration of the image.
 *  @returns the mass and the radius as vectors.
 *  
 */
 vector<float> cog();


		/** Calculate CCF in Fourier space as a function of spatial frequency
                 * between a pair of 2-3D images (corners not included).
		 * The input image 'with' must have the same size to 'this' image.
		 * Input images can be either real or Fourier in arbitrary combination.
		 *
		 * @param[in] with The image used to caculate the fourier shell
		 * @param[in] w Ring/shell width in Fourier space.
		 * @exception ImageFormatException If the 2 images are not
		 * same size.
		 * @exception NullPointerException if the input image is null
		 * @exception Cannot calculate FSC for 1D images
		 * @return  Vector of 3*k FSC results (frequencies, FSC values, error)
		 * k - length of FSC curve, depends on dimensions of the image and ring width
		 * 1 column - normalized frequency [0,0.5]
		 * 2 column - FSC,
		 * 3 column - error of the FSC = 1/sqrt(n), where n is the number of Fourier
		 *            coefficients within given shell.
		 */
		vector <float> calc_fourier_shell_correlation(EMData * with, float w = 1.0f);

		/** Subtract average outside of a circle
		 *  
		 *  
		 *  @return image with sbtracted average outside of a circle..
		 */
		EMData* average_circ_sub();


		/** Helper function for method nn.
		 *
		 * @param j y fourier index (frequency)
		 * @param n number of real elements.
		 * @param n2 Number of complex elements.
		 * @param wptr Normalization matrix [0:n2][1:n][1:n]
		 * @param bi Fourier transform matrix [0:n2][1:n]
		 * @param tf Transform3D reference
		 */
		void onelinenn(int j, int n, int n2, EMData* wptr, EMData* bi, const Transform3D& tf);

		void onelinenn_mult(int j, int n, int n2, EMData* wptr, EMData* bi, const Transform3D& tf, int mult);

		/** Nearest Neighbor interpolation.
		 *  Modifies the current object.
		 *
		 * @param wptr Normalization data.
		 * @param myfft FFT data.
		 * @param tf Transform3D reference
		 * @param mult
		 */
		void nn(EMData* wptr, EMData* myfft, const Transform3D& tf, int mult=1);

		/** Nearest Neighbor interpolation, meanwhile return necessary data such as
		 *  Kn, sum_k(F_k^n) ans sum_k(|F_k^n|^2)
		 *  Modifies the current object.
		 *
		 * @param wptr Normalization data.
		 * @param myfft FFT data.
		 * @param tf Transform3D reference
		 * @param mult
		 */
		void nn_SSNR(EMData* wptr, EMData* wptr2, EMData* myfft, const Transform3D& tf, int mult=1);		

		/** Nearest Neighbor interpolation, meanwhile return necessary data such as
		 *  Kn, sum_k(F_k^n) ans sum_k(|F_k^n|^2)
		 *  Modifies the current object.
		 *
		 * @param wptr Normalization data.
		 * @param myfft FFT data.
		 * @param tf Transform3D reference
		 * @param mult
		 */
		void nn_SSNR_ctf(EMData* wptr, EMData* wptr2, EMData* wptr3, EMData* myfft, const Transform3D& tf, int mult=1);
		
		/** Symmetrize plane 0
		 *  Modifies the current object.
		 *
		 * @param norm Normalization data.
		 */
		void symplane0(EMData* norm);

		/** Symmetrize plane 0
		 *  Modifies the current object.
		 *
		 * @param norm Normalization data.
		 */
		void symplane1(EMData* norm, EMData* norm2);
		
		/** Symmetrize plane 0
		 *  Modifies the current object.
		 *
		 * @param norm Normalization data.
		 */
		void symplane2(EMData* norm, EMData* norm2, EMData* norm3);
		
		/** Helper function for method nn4_ctf.
		 *
		 * @param j y fourier index (frequency)
		 * @param n number of real elements.
		 * @param n2 Number of complex elements.
		 * @param w Normalization matrix [0:n2][1:n][1:n]
		 * @param bi Fourier transform matrix [0:n2][1:n]
		 * @param tf Transform3D reference
		 * @param defocus
		 * @param mult
		 */
		void onelinenn_ctf(int j, int n, int n2, EMData* w, EMData* bi, const Transform3D& tf, float defocus, int mult);

		/** Nearest Neighbor interpolation.
		 *  Modifies the current object.
		 *
		 * @param w Normalization data.
		 * @param myfft FFT data.
		 * @param tf Transform3D reference
		 * @param mult
		 */
		void nn_ctf(EMData* w, EMData* myfft, const Transform3D& tf, int mult);

		/** Helper function for method nn4_ctf.
		 *  here it is assumed the projection data was already multiplied by the ctf...
		 *
		 * @param j y fourier index (frequency)
		 * @param n number of real elements.
		 * @param n2 Number of complex elements.
		 * @param w Normalization matrix [0:n2][1:n][1:n]
		 * @param bi Fourier transform matrix [0:n2][1:n]
		 * @param tf Transform3D reference
		 * @param defocus
		 * @param mult
		 */
		void onelinenn_ctf_applied(int j, int n, int n2, EMData* w, EMData* bi, const Transform3D& tf, float defocus, int mult);

		/** Nearest Neighbor interpolation.
		 *  Modifies the current object.
		 *  here it is assumed the projection data was already multiplied by the ctf...
		 *
		 * @param w Normalization data.
		 * @param myfft FFT data.
		 * @param tf Transform3D reference
		 * @param mult
		 */
		void nn_ctf_applied(EMData* w, EMData* myfft, const Transform3D& tf, int mult );

		/** Symmetrize plane 0
		 *  Modifies the current object.
		 *
		 * @param w Normalization data.
		 */
		void symplane0_ctf(EMData* w);

				
		/** Symmetrize volume in real space.
		 *  
		 *  @param[in] symmetry Point group of the target volume.
		 *  
		 *  @return New symmetrized volume object.
		 */
		EMData* symvol(string symmetry);

		
		/** Rotate-Shift-Scale-Circulantly image 
		 *
		 *  If the image is a volume, then all slices are
		 *  rotated/translated/scaled.
		 *  
		 *  @param[in] ang Rotation angle in degrees.
		 *  @param[in] delx Amount to shift rotation origin along x
		 *  @param[in] dely Amount to shift rotation origin along y
		 *  @param[in] scale Scaling factor (default=1.0)
		 *  @exception ImageDimensionException can not rotate 1 D image
		 *  @exception ImageDimensionException can not rotate 3 D image
		 *  @return New rotated/shifted/scaled image
		 */
		EMData* rot_scale_trans2D(float ang, float delx = 0.f, float dely = 0.f, float scale = 1.f);
		

		/** Rotate-Shift-Scale-Circulantly image 
		 *
		 *  If the image is a volume, then all slices are
		 *  rotated/translated/scaled.
		 *  
		 *  @param[in] RA Transform3D object
		 *  @exception ImageDimensionException can not rotate 1 D image
		 *  @return New rotated/shifted/scaled image
		 */
		EMData* rot_scale_trans(const Transform3D &RA);
		

                /** euclidean distance between two line
		 *
		 */

		float cm_euc(EMData* sinoj, int n1, int n2, float alpha1, float alpha2);


		
		/** Rotate-Shift-Scale-Circulantly image using convolution 
		 *
		 *  If the image is a volume, then all slices are
		 *  rotated/translated/scaled.
		 *  
		 *  @param[in] ang Rotation angle in degrees.
		 *  @param[in] delx Amount to shift rotation origin along x
		 *  @param[in] dely Amount to shift rotation origin along y
		 *  @param[in] kb convolution kernel
		 *  @param[in] scale
		 *  @exception ImageDimensionException can not rotate 1 D image
		 *  @exception ImageDimensionException can not rotate 3 D image
		 *  @return New rotated/shifted/scaled image
		 */
		EMData* rot_scale_conv(float ang, float delx, float dely, Util::KaiserBessel& kb, float scale = 1.0);

		
		/** Get pixel value image using convolution 
		 *
		 *  If the image is a volume, then all slices are
		 *  rotated/translated/scaled.
		 *  
		 *  @param[in] delx Amount to shift rotation origin along x
		 *  @param[in] dely Amount to shift rotation origin along y
		 *  @param[in] delz Amount to shift rotation origin along z
		 *  @param[in] kb convolution kernel
		 *  @exception ImageDimensionException can not rotate 1 D image
		 *  @return New rotated/shifted/scaled image
		 */
		float get_pixel_conv(float delx, float dely, float delz, Util::KaiserBessel& kb);


	
		/** Value of 2-D analytic masking (or 2-D convolution) at off-grid point.
		 *  
		 *  The only requirement for the window function object is that
		 *  it overload operator()(const float) and return a float.
		 *
		 *  @param[in] x x-value of the desired (potentially off-grid) point
		 *  @param[in] y y-value of the desired (potentially off-grid) point
		 *  @param[in] win Window (mask/kernel) function object.
		 *  @param[in] size Size of real-space kernel/mask.
		 *
		 *  @return Value of masked/convolved image at (x,y)
		 */
		/*template<class Win>
		float getconvpt2d(float x, float y, Win win, int size = 7);*/
		float getconvpt2d_kbi0(float x, float y, 
				Util::KaiserBessel::kbi0_win win, int size = 7);
		
		
		/** 2-D rotation using gridding convolution.
		 *  
		 *  The only requirement for the window function object is that
		 *  it overload operator()(const float) and return a float.
		 *
		 *  This routine does _not_ deconvolve out the window function
		 *  after rotation.
		 *
		 *  @param[in] x x-value of the desired (potentially off-grid) point
		 *  @param[in] y y-value of the desired (potentially off-grid) point
		 *  @param[in] win Window (mask/kernel) function object.
		 *  @param[in] size Size of real-space kernel/mask.
		 *  @exception ImageDimensionException only support 2D image
		 *  @return Rotated/convolved EMData image.
		 */
#if 0 // broken
		//template<class Win>
		//EMData* rotconvtrunc2d(float ang, Win win, int size = 7);
		//EMData* rotconvtrunc2d_kbi0(float ang, float alpha, int size) {
		//	Util::KaiserBessel kb(alpha, size-1);
		//	return rotconvtrunc2d(ang, kb.get_kbi0_win(), size);
		//}
		EMData* rotconvtrunc2d_kbi0(float ang, float alpha, int size);
		
		
		/** 2-D rotation using gridding convolution and deconvolution.
		 *  
		 *  This routine does deconvolve out the window function
		 *  after rotation.
		 *
		 *  @param[in] x x-value of the desired (potentially off-grid) point
		 *  @param[in] y y-value of the desired (potentially off-grid) point
		 *  @param[in] win Window (mask/kernel) function object.
		 *  @param[in] size Size of real-space kernel/mask.
		 *
		 *  @return Rotated/convolved EMData image.
		 */
		//template<class Win>
		//EMData* rotconvtrunc2d(float ang, Win win, int size = 7);
		//EMData* rotconvtrunc2d_kbi0(float ang, float alpha, int size) {
		//	Util::KaiserBessel kb(alpha, size-1);
		//	return rotconvtrunc2d(ang, kb.get_kbi0_win(), size);
		//}
		EMData* gridrot2d_kbi0(float ang, float alpha = 1., int size = 7) {
			EMData* rot = rotconvtrunc2d_kbi0(ang, alpha, size);
			Dict params;
			params["dopad"] = 1;
			params["alpha"] = alpha;
			rot->process_inplace("filter.kaisersinhinverse", params);
			return rot;
		}
#endif // 0
		
		
		/** fft_shuffle -- Shuffle a Fourier image to put the origin at (0,ny/2)
		 *  
		 *  Our usual FFT convention puts the origin at (0,0), but then
		 *  grid points corresponding to iy > ny/2 correspond to 
		 *  (unnormalized) frequencies iy-ny.  This routine rearranges
		 *  the columns of the Fourier image so that iy varies from
		 *  -ny/2 to ny/2 (or ny/2 - 1 for ny even).  This method acts
		 *  as a toggle, so to unshuffle a Fourier image just call
		 *  this method a second time.
		 */
		void fft_shuffle();
		
		
		/** center_padded -- Center a padded image
		 *  
		 *  Image padding leaves the image in the corner.  This method
		 *  moves that original image so that it is centered.
		 */
		void center_padded();
		
		
		/** extractpoint -- Gridding convolution
		 *
		 *  Note: Expected to be used in combination with fouriergridrot2d.
		 *
		 *  @param[in] xin x-position
		 *  @param[in] yin y-position
		 *  @param[in] kb  Kaiser-Bessel window 
		 *
		 *  @return Complex gridding result
		 *
		 *  @see P.A. Penczek, R. Renka, and H. Schomberg,
		 *       J. Opt. Soc. Am. A _21_, (2004)
		 *
		 */
		std::complex<float> extractpoint(float xin, float yin, Util::KaiserBessel& kb);
		
		
		/** extractplane -- Gridding convolution in 3D along a plane
		 *
		 *  Note: Expected to be used in combination with fourier gridding
		 *        projections.
		 *
		 *  @param[in] tf  transform matrix defining the intended plane.
		 *  @param[in] kb  Kaiser-Bessel window 
		 *
		 *  @return Complex gridding plane
		 *
		 *  @see P.A. Penczek, R. Renka, and H. Schomberg,
		 *       J. Opt. Soc. Am. A _21_, 499-509 (2004)
		 *
		 */
		EMData*  extractplane(const Transform3D& tf, Util::KaiserBessel& kb);
		
		EMData* fouriergridrot2d(float ang, Util::KaiserBessel& kb);
		
		
		/** divkbsinh -- Divide image by a Kaiser-Bessel sinh window.
		 *
		 *  @param[in] kb Kaiser-Bessel window object
		 *
		 *  Note: Ideally this method really should be a "processor"
		 *        instead, but at the moment a KaiserBessel object
		 *        cannot be passed as part of a Dict, making the usual
		 *        EMData::project() interface rather awkward here.
		 */
		void divkbsinh(const Util::KaiserBessel& kb);
		
		
		/** masked_stats -- Compute image statistics under a mask
		 *
		 *  Specifically, compute the average and standard deviation
		 *  under the mask.  Return the average, the standard deviation,
		 *  and the number of pixels under the mask.
		 *
		 *  @param[in] mask Mask image
		 *
		 *  @return dictionary containing "avg", "sigma", and "nmask" keys
		 */

		//  OBSOLETED  Dict masked_stats(const EMData* mask);


/** Search specified number peaks in 1D, 2D, or 3D real images.
* and output the peaks in descendent order: 
  The numbers coming out are: image dimension, then
  1D: pixel value, x coord, relative peak value, x coord( NX/2 center),
  ...
  2D: pixel value, x coord, y coord, realative peak value, x coord(NX/2 center) y coord(NY/2 center) 
  ...
  3D  pixel value, x coord, y coord, z coord, realative peak value, x coord(NX/2 center) y coord(NY/2 center) z coord(NZ/2 center)
  The function is supposed to return 0 dimension and first pixel value (0,0,0) when the image is contant.
  ...                   */
vector<float> peak_search(int ml, float invert);
/** Calculate the Phase approximation to center of gravity
 *  This operations works for 1-2-3-d images
 *  @returns both the center of gravity and the phase approximated center of gravity values.
 */
vector<float> phase_cog();

float find_3d_threshold(float mass,float pixel_size);

  
 /* Peak (with a radius of hf_p) search for particle picking:
                                                           */ 
vector<float> peak_ccf(float hf_p);
/* pixel power operation function */
EMData* get_pow(float n_pow);

private:
//utility function for peak_search()
static bool peakcmp(const Pixel& p1, const Pixel& p2);  
public:
EMData* extractline(Util::KaiserBessel& kb,float nuxnew,float nuynew);
static EMData* ctf_img(int nx, int ny, int nz, float dz, float ps, float voltage=300.0f,float cs=2.0f,float wgh=0.1f,float b_factor=0.0f,float dza=0.0f,float azz=0.0f,float sign=-1.0f);
#endif	//emdata__sparx_h__
