#
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#

from global_def import *

def absi(e):
	if e.is_complex():
		e.set_attr_dict({"is_complex_ri":1})
	return e.absi()

# Autocorrelation functions
def acf(e, center=True):
	"""
		Name
			acf - calculate the circulant autocorrelation function of an image
		Input
			e: input image, can be either real or Fourier
			center: if set to True (default), the origin of the result is at the center; if set to False, the origin is at (0,0).
		Output
			circulant autocorrelation function of the input image. Real.
	"""
	from EMAN2 import autocorrelation, fp_flag
	return autocorrelation(e, fp_flag.CIRCULANT, center)

def acfn(e, center=True):
	"""
		Name
			acfn - calculate the normalized circulant autocorrelation function of an image
		Input
			e: input image (real)
			center: if set to True (default), the origin of the result is at the center; if set to False, the origin is at (0,0).
		Output
			normalized circulant autocorrelation function of an input image. Real.
	"""
	from EMAN2 import autocorrelation, fp_flag
	return autocorrelation(e, fp_flag.CIRCULANT_NORMALIZED, center)

def acfp(e, center=True):
	from EMAN2 import autocorrelation, fp_flag
	return autocorrelation(e, fp_flag.PADDED, center)

def acfnp(e, center=True):
	"""
		Name
			acfnp - calculate the normalized autocorrelation function of an image
		Input
			e: input image (real)
			center: if set to True (default), the origin of the result is at the center; if set to False, the origin is at (0,0).
		Output
			normalized autocorrelation function of the input image. Real. 
	"""
	from EMAN2 import autocorrelation, fp_flag
	return autocorrelation(e, fp_flag.PADDED_NORMALIZED, center)

def acfpl(e, center=True):
	"""
		Name
			acfpl - calculate the normalized autocorrelation function of an image
		Input
			e: input image (real)
			center: if set to True (default), the origin of the result is at the center; if set to False, the origin is at (0,0).
		Output
			normalized autocorrelation function of the input image. Real. 
			
	"""
	from EMAN2 import autocorrelation, fp_flag
	return autocorrelation(e, fp_flag.PADDED_LAG, center)

def acfnpl(e, center=True):
	"""
		Name
			acfnpl - calculate the normalized autocorrelation function of an image
		Input
			e: input image (real)
			center: if set to True (default), the origin of the result is at the center; if set to False, the origin is at (0,0).
		Output
			autocorrelation function of the input image. Real. 
	"""
	from EMAN2 import autocorrelation, fp_flag
	return autocorrelation(e, fp_flag.PADDED_NORMALIZED_LAG, center)

def __buildweights(m, kb):
	weights = EMData()
	weights.set_size(m,m,1)
	for iy in xrange(m):
		wy = kb.sinhwin(iy-m//2)
		for ix in xrange(m):
			wx = kb.sinhwin(ix-m//2)
			weights.set_value_at(ix,iy,wx*wy)
	return weights

# shortcuts to Fourier product functions
# Correlation functions
def ccf(e, f, center=True):
	"""
	Return the circulant cross-correlation function of images e and f.
	Input images may be real or complex.  Output image is real.
	1-D, 2-D, or 3-D images supported.
	"""
	from EMAN2 import correlation, fp_flag
	return correlation(e,f,fp_flag.CIRCULANT, center)

def ccfn(e, f, center=True):
	"""
		Name
			ccfn - calculate the normalized circulant cross-correlation function between two images.
		Input
			e: input image (real)
			ref: second input image (real) (in the alignment problems, it should be the reference image).
			center: if set to True (default), the origin of the result is at the center
		Output
			normalized circulant cross-correlation function between image and ref. Real.
	"""
	from EMAN2 import correlation, fp_flag
	return correlation(e,f,fp_flag.CIRCULANT_NORMALIZED, center)

def ccfp(e, f, center=True):
	"""
		Name
			ccfp - calculate the cross-correlation function between two images
		Input
			e: input image (real)
			ref: second input image (real)
			center: if set to True (default), the origin of the result is at the center
		Output
			cross-correlation function between image and ref. Real.
	"""
	from EMAN2 import correlation, fp_flag
	return correlation(e,f,fp_flag.PADDED, center)

def ccfnp(e, f, center=True):
	"""
		Name
			ccfnp - calculate the normalized cross-correlation function between two images.
		Input
			e: input image (real)
			ref: second input image (real).
			center: if set to True (default), the origin of the result is at the center
		Output
			normalized cross-correlation function between image and ref. Real.
	"""
	from EMAN2 import correlation, fp_flag
	return correlation(e,f,fp_flag.PADDED_NORMALIZED, center)

def ccfpl(e, f, center=True):
	"""
		Name
			ccfpl - calculate the cross-correlation function between two images	
		Input
			e: input image (real)
			ref: second input image (real) 
			center: if set to True (default), the origin of the result is at the center
		Output
			cross-correlation function between image and ref. Real. 
	"""
	from EMAN2 import correlation, fp_flag
	return correlation(e,f, fp_flag.PADDED_LAG, center)

def ccfnpl(e, f, center=True):
	"""
		Name
			ccfnpl - calculate the normalized cross-correlation function between two images
		Input
			e: input image (real)
			ref: second input image (real) 
			center: if set to True (default), the origin of the result is at the center
	"""
	from EMAN2 import correlation, fp_flag
	return correlation(e,f,fp_flag.PADDED_NORMALIZED_LAG, center)
    
# Convolution functions
def cnv(e, f, center=True):
	"""
		Name
			cnv - calculate the circulant convolution function between two images
		Input
			e: input image, can be either real or Fourier
			ref: second input image, can be either real or Fourier.
			center: if set to True (default), the origin of the result is at the center
		Output
			circulant convolution function between image and ref. Real.
	"""
	from EMAN2 import convolution, fp_flag
	return convolution(e,f,fp_flag.CIRCULANT, center)

def cnvn(e, f, center=True):
	"""
		Name
			cnvn - calculate the normalized circulant convolution function between two images
		Input
			e: input image (real)
			ref: second input image (real).
			center: if set to True (default), the origin of the result is at the center	
		Output
			normalized circulant convolution function between image and ref. Real. 
	"""
	from EMAN2 import convolution, fp_flag
	return convolution(e,f,fp_flag.CIRCULANT_NORMALIZED, center)

def cnvp(e, f, center=True):
	"""
		Name
			cnvp - calculate the convolution function between two images 
		Input
			e: input image (real)
			ref: second input image (real).
			center: if set to True (default), the origin of the result is at the center			
		Output
			convolution function between image and ref. Real.
	"""
	from EMAN2 import convolution, fp_flag
	return convolution(e,f,fp_flag.PADDED, center)

def cnvnp(e, f, center=True):
	"""
		Name
			cnvnp - calculate the normalized convolution function between two images
		Input
			e: input image (real)
			ref: second input image (real) 
			center: if set to True (default), the origin of the result is at the center
		Output
			normalized convolution function between image and ref. Real. 
	"""
	from EMAN2 import convolution, fp_flag
	return convolution(e,f,fp_flag.PADDED_NORMALIZED, center)

def cnvpl(e, f, center=True):
	"""
		Name
			cnvpl - calculate the convolution function between two images
		Input
			e: input image (real)
			ref: second input image (real) 
			center: if set to True (default), the origin of the result is at the center
		Output
			convolution function between image and ref. Real. 
	"""
	from EMAN2 import convolution, fp_flag
	return convolution(e,f,fp_flag.PADDED_LAG, center)

def cnvnpl(e, f, center=True):
	"""
		Name
			cnvnpl - calculate the normalized convolution function between two images
		Input
			e: input image (real)
			ref:second input image (real) 
			center: if set to True (default), the origin of the result is at the center
		Output
			convolution function between image and ref. Real.
	"""
	from EMAN2 import convolution, fp_flag
	return convolution(e,f,fp_flag.PADDED_NORMALIZED_LAG, center)
    
    
# Selfcorrelation functions
def scf(e, center=True):
	"""
		Name
			scf - calculate the circulant self-correlation function of an image
		Input
			e: input image, can be either real or Fourier
			center: if set to True (default), the origin of the result is at the center
		Output
			circulant self-correlation function of the input image. Real.
	"""
	from EMAN2 import self_correlation, fp_flag
	return self_correlation(e, fp_flag.CIRCULANT, center)

def scfn(e, center=True):
	"""
		Name
			scfn - calculate the normalized circulant self-correlation function
		Input
			e: input image, can be either real or Fourier
			center: if set to True (default), the origin of the result is at the center
		Output
			normalized circulant self-correlation function of an input image. Real.
	"""
	from EMAN2 import self_correlation, fp_flag
	return self_correlation(e, fp_flag.CIRCULANT_NORMALIZED, center)

def scfp(e, center=True):
	"""
		Name
			scfp - calculate the self-correlation function of an image
		Input
			e: input image (real)
			center: if set to True (default), the origin of the result is at the center
		Output
			self-correlation function of the input image. Real. 
	"""
	from EMAN2 import self_correlation, fp_flag
	return self_correlation(e, fp_flag.PADDED, center)

def scfnp(e, center=True):
	"""
		Name
			scfnp - calculate the normalized self-correlation function of an image
		Input
			e: input image (real)
			center: if set to True (default), the origin of the result is at the center
		Output
			normalized self-correlation function of the input image. Real.
	"""
	from EMAN2 import self_correlation, fp_flag
	return self_correlation(e, fp_flag.PADDED_NORMALIZED, center)

def scfpl(e, center=True):
	"""
		Name
			scfpl - calculate the self-correlation function of an image
		Input
			e: input image (real)
			center: if set to True (default), the origin of the result is at the center
		Output
			self-correlation function of the input image. Real.
	"""
	from EMAN2 import self_correlation, fp_flag
	return self_correlation(e, fp_flag.PADDED_LAG, center)

def scfnpl(e, center=True):
	"""
		Name
			scfnpl - calculate the normalized self-correlation function of an image
		Input
			e: input image (real)
			center: if set to True (default), the origin of the result is at the center
		Output
			self-correlation function of the input image. Real.
	"""
	from EMAN2 import self_correlation, fp_flag
	return self_correlation(e, fp_flag.PADDED_NORMALIZED_LAG, center)
 
def cyclic_shift(img, dx=0, dy=0, dz=0):
	"""
		Name
			cyclic_shift - cyclically shift an image by an integer number of pixels
		Input
			image: input image
			ix: (optional) integer number of pixels by which the image has to be shifted in the x direction, can be negative
			iy: (optional) integer number of pixels by which the image has to be shifted in the y direction, can be negative
			iz: (optional) integer number of pixels by which the image has to be shifted in the z direction, can be negative
		Output
			output image
	"""
	e = img.copy()
	Util.cyclicshift(e,{"dx":dx,"dy":dy,"dz":dz})
	return e

def mirror(img, axis = 'x'):
	"""Mirror of an image about by changing the sign of the specified axis
	"""
	e = img.copy()
	e.process_inplace("xform.mirror", {"axis":axis})
	return e 

# FFT functions
def fft_(e, npad=1):
	"""Out-of-place fft / ift
		No padding performed, and fft-extension along x removed after ift. Zhong added in July,5,06
	"""
	if (e.is_complex()):
		return e.do_ift()
	else:
		if npad > 1:
			f = e.norm_pad(False, npad)
			f.do_fft_inplace()
			return f
		else:
			f = e.norm_pad(False, 1)
			f.do_fft_inplace()
			return f

def fft(e, npad=1):
	"""Out-of-place fft / ift
	   No padding performed, and fft-extension along x removed after ift.
	"""
	if (e.is_complex()):
		return  e.do_ift()
	else:
		# forward fft
		if npad > 1:
			f = e.norm_pad(False, npad)
			f.do_fft_inplace()
			return f
		else:
			return e.do_fft()

def fftip(e):
	"""In-place fft / ift
	   No padding performed, and fft-extension along x NOT removed after ift.
	"""
	if (e.is_complex()):
		# inverse fft
		e.do_ift_inplace()
		#e.postift_depad_corner_inplace()
	else:
		# forward fft
		e.do_fft_inplace()

def fpol(image, nnx, nny=1, nnz=1, RetReal = True):
	"""
		Name
			fpol -Interpolate image up by padding its Fourier transform with zeroes
		Input
			image: image to be interpolated.
			nnx: new nx dimension
			nny: new ny dimension (default = 0, meaning equal to the original ny)
			nnz: new nz dimension (default = 0, meaning equal to the original nz)
			RetReal: Logical flag, if True, the returned image is real, if False, it is Fourier
		Output
			the output interpolated up image
	"""
	from fundamentals import fft
	
	nx = image.get_xsize()
	ny = image.get_ysize()
	nz = image.get_zsize()
	
	if image.is_complex():
		nx -= (2-nx%2)
	
	if nx == nnx and ny == nny and nz == nnz:
		if image.is_complex() and RetReal: return fft(image)
		else: return image
	
	return  image.FourInterpol(nnx, nny, nnz, RetReal)

def fdecimate(image, nnx, nny=1, nnz=1, RetReal = True):
	"""
		Decimate image by truncating its Fourier transform
	"""
	return  image.FourTruncate(nnx, nny, nnz, RetReal)	
    
def fshift(e, delx=0, dely=0, delz=0):
	"""
		Name
			fshift - shift an image using phase multiplication in Fourier space
		Input
			image: input image
			sx: shift in pixels in the x direction, can be negative
			sy: shift in pixels in the y direction, can be negative
			sz: shift in pixels in the z direction, can be negative
		Output
			output image
	"""
	from EMAN2 import Processor

	params = {"filter_type" : Processor.fourier_filter_types.SHIFT,	"x_shift" : float(delx), "y_shift" : float(dely), "z_shift" : float(delz) }
	return Processor.EMFourierFilter(e, params)

def image_decimate(img, decimation=2, fit_to_fft = True, frequency_low=0, frequency_high=0):
	"""
		Window 2D image to FFT-friendly size, apply Butterworth low pass filter,
		and decimate image by integer factor
	"""
	from filter       import filt_btwl
	from fundamentals import smallprime, window2d
	from utilities    import get_image
	if type(img)     == str:	img=get_image(img)
	nz       = img.get_zsize()
	if( nz > 1):                    ERROR("This command works only for 2-D images", "image_decimate", 1)
	if decimation    <= 1  : 	ERROR("Improper decimation ratio", "image_decimate", 1)
	if(decimation    == 1.0): 	return  img.copy()
	if frequency_low <= 0  :	
		frequency_low     = 0.5/decimation-0.02
		if frequency_low <= 0 : ERROR("Butterworth pass-band frequency is too low","image_decimation",1)			
		frequency_high    = min(0.5/decimation + 0.02, 0.499)
	if fit_to_fft:
		nx       = img.get_xsize()
		ny       = img.get_ysize()
		nx_fft_m = smallprime(nx)
		ny_fft_m = smallprime(ny)
		e        = Util.window(img, nx_fft_m, ny_fft_m, 1, 0,0,0)
		e        = filt_btwl(e, frequency_low, frequency_high)
	else:
		e        = filt_btwl(img, frequency_low, frequency_high)
	return Util.decimate(e, int(decimation), int(decimation), 1)

def subsample(image, subsample_rate=1.0):
	if subsample_rate == 1.0: return  image.copy()
	template_min = 15
	frequency_cutoff = 0.5*subsample_rate
	sb = Util.sincBlackman(template_min, frequency_cutoff, 1999) # 1999 taken directly from util_sparx.h
	return image.downsample(sb, subsample_rate)

def resample(img, sub_rate=0.5):
	"""
		resample image based on the value of sub_rate.
		the input image can be either 2D image or 3D volume.
		sub_rate < 1.0, subsampling the image.
		sub_rate > 1.0, upsampling the image using new gridding interpolation.
		fit_to_fft will change the ouput image size to an fft_friendly size
	"""

	from fundamentals import subsample
	from utilities    import get_pixel_size, set_pixel_size

	if type(img) == str:
		from utilities    import get_image
		img = get_image(img)
	nx = img.get_xsize()
	ny = img.get_ysize()
	nz = img.get_zsize()
	if( ny == 1):  ERROR("Only 2D or 3D images allowed","resample",1)
	if sub_rate == 1.0: return  img.copy()
	elif sub_rate < 1.0:
		e = subsample(img, sub_rate)
	else:  #  sub_rate>1
		new_nx = int(nx*sub_rate+0.5)
		new_ny = int(ny*sub_rate+0.5)
		if nz==1:
			new_nz = 1
		else:
			new_nz = int(ny*sub_rate+0.5)
		if ( nx!=ny and nz==1 ):
			nn = max(new_nx, new_ny)
			e = Util.pad(img, nn, nn,  1, 0, 0, 0, "circumference")
			e, kb = prepi(e)
			e = Util.window( e.rot_scale_conv_new(0.0, 0.0, 0.0, kb, sub_rate), new_nx, new_ny, 1, 0,0,0)
		 
		elif ((nx!=ny or nx!=nz or ny!=nz) and nz>1):
			nn = max(new_nx, new_ny,new_nz)
			e = Util.pad(img, nn, nn,  nn, 0, 0, 0, "circumference")
			e, kb = prepi3D(e)
			e = Util.window( e.rot_scale_conv_new_3D(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, kb, sub_rate), new_nx, new_ny, new_nz, 0,0,0)
		else:
			if nz==1:
				e, kb = prepi(Util.pad(img, new_nx, new_ny, 1, 0, 0, 0, "circumference"))
				e = e.rot_scale_conv_new(0.0, 0.0, 0.0, kb, sub_rate)
			else:
				e, kb = prepi3D(Util.pad(img, new_nx, new_ny, new_nz, 0, 0, 0, "circumference"))
				e = e.rot_scale_conv_new_3D(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, kb, sub_rate)

	# Automatically adjust pixel size for ctf parameters
	from utilities import get_pixel_size, set_pixel_size
	apix = get_pixel_size(e)
	apix /= sub_rate
	set_pixel_size(e, apix)
	cc = e.get_attr_default("xform.projection", None)
	if cc:
		cp = cc.get_params("spider")
		cp["tx"] *= sub_rate
		cp["ty"] *= sub_rate
		from utilities import set_params_proj
		set_params_proj(e, [cp["phi"], cp["theta"], cp["psi"], -cp["tx"], -cp["ty"]]) # have to invert as set inverts them again
	cc = e.get_attr_default("xform.align2d", None)
	if cc:
		cp = cc.get_params("2D")
		cp["tx"] *= sub_rate
		cp["ty"] *= sub_rate
		from utilities import set_params2D
		set_params2D(e, [cp["alpha"], cp["tx"], cp["ty"], cp["mirror"], cp["scale"]])

	return 	e
	


def fdownsample(img, sub_rate=0.5, RetReal = True):
	"""
		resample image based on the value of sub_rate.
		the input image can be either 2D image or 3D volume.
		sub_rate < 1.0, subsampling the image.
		sub_rate > 1.0, upsampling the image using new gridding interpolation.
		fit_to_fft will change the ouput image size to an fft_friendly size
	"""

	from fundamentals import fdecimate
	from utilities    import get_pixel_size, set_pixel_size

	if type(img) == str:
		from utilities    import get_image
		img = get_image(img)
	nx = img.get_xsize()
	if img.is_complex():
		nx -= (2-nx%2)
	ny = img.get_ysize()
	nz = img.get_zsize()
	if( ny == 1):  ERROR("Only 2D or 3D images allowed","resample",1)
	if sub_rate == 1.0: return  img.copy()
	elif sub_rate < 1.0:
		nnx = int(nx*sub_rate+0.5)
		nny = int(ny*sub_rate+0.5)
		nnz = int(nz*sub_rate+0.5)
		e = fdecimate(img, nnx, nny, nnz, RetReal = RetReal)
	else:  #  sub_rate>1
		ERROR("fdownsample","upscaling not implemented",1)
		"""
		new_nx = int(nx*sub_rate+0.5)
		new_ny = int(ny*sub_rate+0.5)
		if nz==1:
			new_nz = 1
		else:
			new_nz = int(ny*sub_rate+0.5)
		if ( nx!=ny and nz==1 ):
			nn = max(new_nx, new_ny)
			e = Util.pad(img, nn, nn,  1, 0, 0, 0, "circumference")
			e, kb = prepi(e)
			e = Util.window( e.rot_scale_conv_new(0.0, 0.0, 0.0, kb, sub_rate), new_nx, new_ny, 1, 0,0,0)
		 
		elif ((nx!=ny or nx!=nz or ny!=nz) and nz>1):
			nn = max(new_nx, new_ny,new_nz)
			e = Util.pad(img, nn, nn,  nn, 0, 0, 0, "circumference")
			e, kb = prepi3D(e)
			e = Util.window( e.rot_scale_conv_new_3D(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, kb, sub_rate), new_nx, new_ny, new_nz, 0,0,0)
		else:
			if nz==1:
				e, kb = prepi(Util.pad(img, new_nx, new_ny, 1, 0, 0, 0, "circumference"))
				e = e.rot_scale_conv_new(0.0, 0.0, 0.0, kb, sub_rate)
			else:
				e, kb = prepi3D(Util.pad(img, new_nx, new_ny, new_nz, 0, 0, 0, "circumference"))
				e = e.rot_scale_conv_new_3D(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, kb, sub_rate)
		"""

	# Automatically adjust pixel size for ctf parameters
	from utilities import get_pixel_size, set_pixel_size
	apix = get_pixel_size(e)
	apix /= sub_rate
	set_pixel_size(e, apix)
	cc = e.get_attr_default("xform.projection", None)
	if cc:
		cp = cc.get_params("spider")
		cp["tx"] *= sub_rate
		cp["ty"] *= sub_rate
		from utilities import set_params_proj
		set_params_proj(e, [cp["phi"], cp["theta"], cp["psi"], -cp["tx"], -cp["ty"]]) # have to invert as set inverts them again
	cc = e.get_attr_default("xform.align2d", None)
	if cc:
		cp = cc.get_params("2D")
		cp["tx"] *= sub_rate
		cp["ty"] *= sub_rate
		from utilities import set_params2D
		set_params2D(e, [cp["alpha"], cp["tx"], cp["ty"], cp["mirror"], cp["scale"]])

	return 	e


def prepi(image):
	"""
		Name
			prepi - prepare 2-D image for rotation/shift
		Input
			image: input image that is going to be rotated and shifted using rtshgkb
		Output
			imageft: real space image prepared for gridding rotation and shift by convolution
			kb: interpolants (tabulated Kaiser-Bessel function)
	"""
	from EMAN2 import Processor

	M = image.get_xsize()
# pad two times
	npad = 2
	N = M*npad
# support of the window
	K = 6
	alpha = 1.75
	r = M/2
	v = K/2.0/N
	kb = Util.KaiserBessel(alpha, K, r, v, N)
	#out = rotshift2dg(image, angle*pi/180., sx, sy, kb,alpha)
# first pad it with zeros in Fourier space
	q = image.FourInterpol(N, N, 1, 0)
	params = {"filter_type": Processor.fourier_filter_types.KAISER_SINH_INVERSE,
	          "alpha":alpha, "K":K, "r":r, "v":v, "N":N}
	q = Processor.EMFourierFilter(q, params)
	params = {"filter_type" : Processor.fourier_filter_types.TOP_HAT_LOW_PASS,
		"cutoff_abs" : 0.25, "dopad" : False}
	q = Processor.EMFourierFilter(q, params)
	return fft(q), kb

def prepi3D(image):
	"""
		Name
			prepi3D - prepare 3-D image for rotation/shift
		Input
			image: input image that is going to be rotated and shifted using rot_shif3D_grid
		Output
			imageft: image prepared for gridding rotation and shift
			kb: interpolants (tabulated Kaiser-Bessel function)
	"""
	from EMAN2 import Processor

	M = image.get_xsize()
	# padding size:
	npad = 2
	N = M*npad
	# support of the window:
	K = 6
	alpha = 1.75
	r = M/2
	v = K/2.0/N
	kb = Util.KaiserBessel(alpha, K, r, v, N)
	# pad with zeros in Fourier space:
	q = image.FourInterpol(N, N, N, 0)
	params = {"filter_type": Processor.fourier_filter_types.KAISER_SINH_INVERSE,
	          "alpha":alpha, "K":K, "r":r, "v":v, "N":N}
	q = Processor.EMFourierFilter(q, params)
	params = {"filter_type" : Processor.fourier_filter_types.TOP_HAT_LOW_PASS,
		"cutoff_abs" : 0.25, "dopad" : False}
	q = Processor.EMFourierFilter(q, params)
	return fft(q), kb

def prepg(image, kb):
	"""
		Name
			prepg - prepare 2-D image for rotation/shift using gridding method.
		Input
			image: input image that is going to be rotated and shifted using rtshgkb
			kb: interpolants (tabulated Kaiser-Bessel function)
		Output
			imageft: image prepared for gridding rotation and shift
	"""
	from EMAN2 import Processor

	M = image.get_xsize()
# padd two times
	npad = 2
	N = M*npad
# support of the window
	K = 6
	alpha = 1.75
	r = M/2
	v = K/2.0/N
# first pad it with zeros in Fourier space
	o = image.FourInterpol(2*M,2*M,1,0)
	params = {"filter_type" : Processor.fourier_filter_types.KAISER_SINH_INVERSE,
		  "alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
	q = Processor.EMFourierFilter(o,params)
	return  fft(q)
	
def ramp(inputimage):
	"""
		Name
			ramp - remove linear trend from a 2D image
		Input
			image: 2D input image
		Output
			2D output image from which a 2D plane was subtracted
	"""
	e = inputimage.copy()
	e.process_inplace("filter.ramp")
	return e

def rot_avg(e):
	"""Rotational average.
	   Returns a 1-D image containing a rotational average of image e.
	"""
	return e.rotavg()

def rot_avg_table(e):
	"""Rotational average.
	   Returns a table containing a rotational average of image e.
	"""
	qt = e.rotavg()
	tab = []
	n = qt.get_xsize()
	for i in xrange(n):
		tab.append(qt.get_value_at(i))
	return tab

def rot_avg_image(image_to_be_averaged):
	"""
	Rotational average
	Returns a 2-D or 3-D image containing a rotational average of image e
	"""
	import types
	from utilities import get_im
	if type(image_to_be_averaged) is types.StringType: image_to_be_averaged = get_im(image_to_be_averaged)
	return image_to_be_averaged.rotavg_i()

def ro_textfile(e, filename, helpful_string=""):
	"""Rotational average stored as a text file.
	   Saves a text file (suitable for gnuplot) of the rotational average of e.
	"""
	out = open(filename, "w")
	out.write("#Rotational average: %s\n" % (helpful_string));
	f = e.rotavg()
	nr = f.get_xsize()
	for ir in xrange(nr):
		out.write("%d\t%12.5g\n" % (ir, f.get_value_at(ir)))
	out.close()

def rops(e):
	"""Rotational average of the power spectrum.
	   Returns a 1-D image containing a rotational average
	   of the periodogram of image e.
	"""
	from EMAN2 import periodogram
	ps = periodogram(e)
	return ps.rotavg()

def rops_textfile(e, filename, helpful_string="", lng = False):
	"""Rotational average of the periodogram stored as a text file.
	   Saves a text file (suitable for gnuplot) of the rotational average 
	   of the periodogram of image e.
	"""
	from EMAN2 import periodogram
	out = open(filename, "w")
	if helpful_string != "": out.write("#Rotational average: %s\n" % (helpful_string))
	ps = periodogram(e)
	f = ps.rotavg()
	nr = f.get_xsize()
	table = [0.0]*nr
	for ir in xrange(nr): table[ir] = f.get_value_at(ir)
	if lng:
		from math import log
		for ir in xrange(1,nr): table[ir] = log(table[ir])
		table[0] = table[1]
	for ir in xrange(nr): out.write("%d\t%12.5g\n" % (ir, table[ir]))
	out.close()
	
def rops_table(img, lng = False):

	""" 
		Calculate 1D rotationally averaged 
		power spectrum and save it in list
	"""
	from EMAN2 import periodogram
	e = periodogram(img)
	ro = e.rotavg()
	nr = ro.get_xsize()
	table = [0.0]*nr
	for ir in xrange(nr): table[ir] = ro.get_value_at(ir)
	if lng:
		from math import log10
		for ir in xrange(1,nr): table[ir] = log10(table[ir])
		table[0] = table[1]
	return table

def rops_dir(indir, output_dir = "1dpw2_dir"):
	"""
		Calculate 1D rotationally averaged power spectra from
		image stack listed in a directory
	"""
	from EMAN2 import periodogram
	import os
	flist = os.listdir(indir)
	print flist
	if os.path.exists(output_dir) is False: os.mkdir(output_dir)
	for i, v in enumerate(flist):
		(filename, filextension) = os.path.splitext(v)
		nima = EMUtil.get_image_count(os.path.join(indir,v))
		print nima
		for im in xrange(nima):
			e = EMData()
			file_name = os.path.join(indir,v)
			e.read_image(file_name, im)
			tmp1 = periodogram(e)
			tmp  = tmp1.rotavg()
			if im == 0:
				sum_ima  = model_blank(tmp.get_xsize())
				sum_ima += tmp
			else :  sum_ima += tmp
		table = []
		nr = sum_ima.get_xsize()
		for ir in xrange(nr):  table.append([sum_ima.get_value_at(ir)])
		drop_spider_doc(os.path.join(output_dir, "1dpw2_"+filename+".txt"), table)


def rotshift2dg(image, ang, dx, dy, kb, scale = 1.0):
	"""Rotate and shift an image using gridding
	"""
	from math import radians
	from EMAN2 import Processor

	M = image.get_xsize()
	alpha = 1.75
	K = 6
	N = M*2  # npad*image size
	r = M/2
	v = K/2.0/N
	# first pad it with zeros in Fourier space
	o = image.FourInterpol(N,N,1,0)
	params = {"filter_type" : Processor.fourier_filter_types.KAISER_SINH_INVERSE,
	          "alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
	q = Processor.EMFourierFilter(o,params)
	o = fft(q)

	# gridding rotation
	return o.rot_scale_conv(radians(ang), dx, dy, kb, scale)

def gridrot_shift2D(image, ang = 0.0, sx = 0.0, sy = 0.0, scale = 1.0):
	"""
		Rotate and shift an image using gridding in Fourier space.
	"""
	from EMAN2 import Processor
	from fundamentals import fftip, fft

	nx = image.get_xsize()
	# split shift into integer and fractional parts
	isx = int(sx)
	fsx = sx - isx
	isy = int(sy)
	fsy = sy - isy
	# prepare 
	npad = 2
	N = nx*npad
	K = 6
	alpha = 1.75
	r = nx/2
	v = K/2.0/N
	kb = Util.KaiserBessel(alpha, K, r, v, N)

	image1 = image.copy()  # This step is needed, otherwise image will be changed outside the function
	# divide out gridding weights
	image1.divkbsinh(kb)
	# pad and center image, then FFT
	image1 = image1.norm_pad(False, npad)
	fftip(image1)
	# Put the origin of the (real-space) image at the center
	image1.center_origin_fft()
	# gridding rotation
	image1 = image1.fouriergridrot2d(ang, scale, kb)
	if(fsx != 0.0 or fsy != 0.0):
		params = {"filter_type" : Processor.fourier_filter_types.SHIFT,	"x_shift" : float(fsx), "y_shift" : float(fsy), "z_shift" : 0.0 }
		image1 = Processor.EMFourierFilter(image1, params)
	# put the origin back in the corner
	image1.center_origin_fft()
	# undo FFT and remove padding (window)
	image1 = fft(image1)
	image1 = image1.window_center(nx)
	Util.cyclicshift(image1,{"dx":isx,"dy":isy,"dz":0})
	return image1

def rot_shift2D(img, angle = 0.0, sx = 0.0, sy = 0.0, mirror = 0, scale = 1.0, interpolation_method = None, mode = "background"):
	"""
		rotate/shift image using:
		1. linear    interpolation
		2. quadratic interpolation
		3. gridding
		mode specifies what will be put in corners, should they stick out:
		background - leave original values
		cyclic - use pixels from the image using wrap-around transformation
		
	"""

	if scale == 0.0 :  ERROR("0 scale factor encountered","rot_shift2D", 1)
	if(interpolation_method):  use_method = interpolation_method
	else:  use_method = interpolation_method_2D

	if(use_method == "linear" and mode == "cyclic"):
		T  = Transform({'type': 'SPIDER', 'psi': angle, 'tx': sx, 'ty': sy, 'scale':scale})
		img = img.rot_scale_trans(T)
		if  mirror: img.process_inplace("xform.mirror", {"axis":'x'})
		return img
	elif(use_method == "linear" and mode == "background"):
		T  = Transform({'type': 'SPIDER', 'psi': angle, 'tx': sx, 'ty': sy, 'scale':scale})
		img = img.rot_scale_trans_background(T)
		if  mirror: img.process_inplace("xform.mirror", {"axis":'x'})
		return img
	elif(use_method == "quadratic" and mode == "cyclic"):
		img = img.rot_scale_trans2D(angle, sx, sy, scale)
		if  mirror: img.process_inplace("xform.mirror", {"axis":'x'})
		return img
	elif(use_method == "quadratic" and mode == "background"):
		img = img.rot_scale_trans2D_background(angle, sx, sy, scale)
		if  mirror: img.process_inplace("xform.mirror", {"axis":'x'})
		return img
	elif(use_method == "gridding" and mode == "cyclic"): # previous rtshg
		from math import radians
		from fundamentals import prepi
		o, kb = prepi(img)
		# gridding rotation
		o = o.rot_scale_conv_new(radians(angle), sx, sy, kb, scale)
		if  mirror: o.process_inplace("xform.mirror", {"axis":'x'})
		return o
	elif(use_method == "gridding" and mode == "background"): # previous rtshg
		from math import radians
		from fundamentals import prepi
		o, kb = prepi(img)
		# gridding rotation
		o = o.rot_scale_conv_new_background(radians(angle), sx, sy, kb, scale)
		if  mirror: o.process_inplace("xform.mirror", {"axis":'x'})
		return o	
	elif(use_method == "ftgridding"): # previous rtshg
		from fundamentals import gridrot_shift2D
		img = gridrot_shift2D(img, angle, sx, sy, scale)
		if  mirror: img.process_inplace("xform.mirror", {"axis":'x'})
		return img
	else:	ERROR("rot_shift_2D interpolation method is incorrectly set", "rot_shift_2D", 1)

def rot_shift3D(image, phi = 0, theta = 0, psi = 0, sx = 0, sy = 0, sz = 0, scale = 1.0, mode="background"):
	"""
		Name
			rot_shift3D - rotate, shift, and scale a 3D image
		Input
			image:3D image to be rotated
			phi, theta, psi: Euler angles ( in SPIDER convention ) describing rotation
			sx, sy, sz: x, y and z shifts, respectively
			scale: rescaling factor. scale > 1 will magnify the input volume.
			mode: backgroud
		Output
			The rotated, shifted, and scaled output 3D volume
	"""

	if scale == 0.0 :  ERROR("0 scale factor encountered","rot_shift3D", 1)
	T1 = Transform({'scale':scale})
	T2 = Transform({'type': 'SPIDER', 'phi': phi, 'theta': theta, 'psi': psi, 'tx': sx, 'ty': sy, 'tz': sz})
	T  = T1*T2
	if (mode == "cyclic"): return image.rot_scale_trans(T)
	else: return image.rot_scale_trans_background(T)


def rot_shift3D_grid(img, phi=0.0, theta=0.0, psi=0.0, sx=0.0, sy=0.0, sz=0.0, scale=1.0, kb=None, mode="background", wrap=False):
	"""
		rotate/shift/scale image using the gridding method.
		if kb = None, the image is prepared in this function before the rot/shift operation.
		If kb is NOT None, then the supplied (img,kb) must be the output from prepi3D.
		'mode' specifies what will be put in corners, should they stick out:
			background - leave original values
			cyclic - use pixels from the image using wrap-around transformation
		'wrap': option for using wraparound pixels during translations
	"""

	if scale == 0.0 :  ERROR("scale=0 not allowed", "rot_shift3D_grid", 1)

	if mode == "cyclic":
		from math import radians
		if kb == None:
			o, kb = prepi3D(img)
		else:
			o = img
		# gridding rotation/shift:
		#if  mirror: o.process_inplace("xform.mirror", {"axis":'x'})
		return o.rot_scale_conv_new_3D(radians(phi), radians(theta), radians(psi), sx, sy, sz, kb, scale, wrap)
	elif mode == "background":
		from math import radians
		if kb == None:
			o, kb = prepi3D(img)
		else:
			o = img
		# gridding rotation
		#if  mirror: o.process_inplace("xform.mirror", {"axis":'x'})
		return o.rot_scale_conv_new_background_3D(radians(phi), radians(theta), radians(psi), sx, sy, sz, kb, scale, wrap)	
	else: ERROR("rot_shift3D_grid mode not valid", "rot_shift3D_grid", 1)


def rtshg(image, angle = 0.0, sx=0.0, sy=0.0, scale = 1.0):
	"""
		Name
			rtshg - rotate and shift a 2D image using the gridding method
		Input
			image: a 2D input image
			alpha: rotation angle
			sx: x shift (default = 0.0)
			sy: y shift (default = 0.0)
			scale: magnification change (default = 1.0)
		Output
			the output rotated and shifted image
	"""
	from math import radians
	o,kb = prepi(image)
	# gridding rotation
	return o.rot_scale_conv_new(radians(angle), sx, sy, kb, scale)

def rtshgkb(image, angle, sx, sy, kb, scale = 1.0):
	"""
		Name
			rtshgkb - rotate and shift a 2D image using the gridding method.
		Input
			sx: x shift
			sy: y shift
			kb: interpolants
		Output
			the output rotated and shifted image
	"""
	from math import radians
	# gridding rotation
	return image.rot_scale_conv_new(radians(angle), sx, sy, kb, scale)
	
def smallprime(arbit_num, numprime=3):
	primelist = [2,3,5,7,11,13,17,19,23]
	for i in xrange(1, arbit_num+1):
		x62 = arbit_num-i+1
		for k in xrange(1,arbit_num+1): # fake loop try to divide the arbit_num
			for j in xrange(0,numprime):
				x71 = primelist[j]*int(x62/primelist[j])
				if(x71 == x62):
					x62 = x62/primelist[j]
					if(x62 == 1):
						nicenum = arbit_num-i+1
						return nicenum
					else:
						break
		 
	nicenum = arbit_num-i+1
	return nicenum

def welch_pw2(img, win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0):
	""" 
		Calculate the power spectrum using Welch periodograms (overlapped periodogram)
	"""
	from fundamentals import window2d, ramp
	from EMAN2 import periodogram
	nx = img.get_xsize()
	ny = img.get_ysize()
	nx_fft = smallprime(nx)
	ny_fft = smallprime(ny)
	x_gaussian_hi = 1./win_size
	from filter    import filt_gaussh
	e_fil = filt_gaussh(window2d(img,nx_fft,ny_fft,"l"), x_gaussian_hi)
	x38 = 100/(100-overlp_x) # normalization of % of the overlap in x 
	x39 = 100/(100-overlp_y) # normalization of % of the overlap in y
	x26 = int(x38*((nx-2*edge_x)/win_size-1)+1)  # number of pieces horizontal dim.(X)
	x29 = int(x39*((ny-2*edge_y)/win_size-1)+1)  # number of pieces vertical dim.(Y)
	iz = 0	
	pw2 = EMData()
	for iy in xrange(1, x29+1):	
		x21 = (win_size/x39)*(iy-1) + edge_y  #  y-direction it should start from 0 if edge_y=0	      
		for ix in  xrange(1, x26+1):			 
			x22 = (win_size/x38)*(ix-1) + edge_x  # x-direction it should start from 0 if edge_x =0
			wi  = window2d(e_fil, win_size, win_size, "l", x22, x21)
			iz  = iz+1
			if (iz == 1): pw2  = periodogram(ramp(wi))
			else:         pw2 += periodogram(ramp(wi))
	return  pw2/float(iz)

def welch_pw2_tilt_band(img,theta,num_bnd=-1,overlp_y=50,edge_x=0,edge_y=0,win_s=256):
	""" 
		1. Calculate the power spectra of tilt bands
		2. The tilt micrograph is rotated such that the tilt axis is vertical (along Y axis)
		3. edge_x and edge_y are removed from the micrograph
	""" 
	from EMAN2 import periodogram
	nx = img.get_xsize()
	ny = img.get_ysize()
	num1 = int(nx-2*edge_x)
	num2 = int(ny-2*edge_y)
	nx_fft = smallprime(num1)
	ny_fft = smallprime(num2)
	img1 = window2d(img,nx_fft,ny_fft,"l",edge_x,edge_y)
	if(num_bnd == -1):
		num_bnd = int(nx_fft/win_s)
		win_x   = int(win_s)
	else:
		win_x = int(nx_fft/num_bnd)
		win_x = int(smallprime(win_x))
	win_y = win_x
	x_gaussian_hi = 1./win_x
	del img
	from filter import filt_gaussh
	from utilities import drop_image, rot_image
	# The input img is rotated such that tilt axis is vertical
	img2  = rot_image(img1,theta, 0, 0, 1.0,1.0)	
	e_fil = filt_gaussh(img2, x_gaussian_hi)
	del img1
	del img2
	x39 = 100/(100-overlp_y) # normalization of % of the overlap in y
	x29 = int(x39*((ny)/win_y-1)+1)  # number of pieces vertical dim.(Y)
	pw2 = EMData()
	pw2_band = []
	for ix in  xrange(1, num_bnd+1):
		x22 = (win_x)*(ix-1)# x-direction it should start from 0 if edge_x =0
		iz=0
		for iy in xrange(1, x29+1):	
			x21 = (win_y/x39)*(iy-1) #  y-direction it should start from 0 if edge_y=0	      			 
			wi = window2d(e_fil,win_x, win_y,"l",x22, x21)
			iz = iz+1
			if (iz == 1): pw2  = periodogram(ramp(wi))
			else:         pw2 += periodogram(ramp(wi))
		pw2/=float(iz)
		# drop_image(pw2,"band%03d"%(ix))
		pw2_band.append(pw2)	
	return 	pw2_band


def tilemic(img, win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0):
	""" 
		Calculate the power spectrum using Welch periodograms (overlapped periodogram)
	"""
	from fundamentals import window2d, ramp
	from EMAN2 import periodogram
	nx = img.get_xsize()
	ny = img.get_ysize()
	nx_fft = smallprime(nx)
	ny_fft = smallprime(ny)
	x_gaussian_hi = 1./win_size
	from filter    import filt_gaussh
	e_fil = filt_gaussh(window2d(img,nx_fft,ny_fft,"l"), x_gaussian_hi)
	x38 = 100/(100-overlp_x) # normalization of % of the overlap in x 
	x39 = 100/(100-overlp_y) # normalization of % of the overlap in y
	x26 = int(x38*((nx-2*edge_x)/win_size-1)+1)  # number of pieces horizontal dim.(X)
	x29 = int(x39*((ny-2*edge_y)/win_size-1)+1)  # number of pieces vertical dim.(Y)
	pw2 = []
	for iy in xrange(1, x29+1):	
		x21 = (win_size/x39)*(iy-1) + edge_y  #  y-direction it should start from 0 if edge_y=0	      
		for ix in  xrange(1, x26+1):			 
			x22 = (win_size/x38)*(ix-1) + edge_x  # x-direction it should start from 0 if edge_x =0
			wi  = ramp( window2d(e_fil, win_size, win_size, "l", x22, x21) )
			st = Util.infomask(wi, None, True)
			wi = (wi - st[0])/st[1]*win_size
			pw2.append(periodogram(wi))
	return  pw2


def window2d(img, isize_x, isize_y, opt="c", ix=0, iy=0):
	"""
		Three ways of windowing a portion of image from a large image field:
		1. "c" Get the central part: "c" ( used for reduce image size )
		2. "l" Get clip starts from the top left corner
		3. "a" Get clip with arbitrary point (ix, iy) as the image center point ( nx//2,ny//2 corresponds to image center )
	"""
	lx = img.get_xsize()
	ly = img.get_ysize()
	if(lx == isize_x and ly == isize_y):  return img.copy()
	from EMAN2 import Region
	if(opt == "l"): reg = Region(ix, iy, isize_x, isize_y)
	elif(opt == "c"):
		mx = lx//2-isize_x//2
		my = ly//2-isize_y//2
		reg = Region(mx, my, isize_x, isize_y)
	elif(opt == "a"):
		mx = ix-isize_x//2
		my = iy-isize_y//2
		reg = Region(mx, my, isize_x, isize_y)
	else:  ERROR("Unknown window2d option","window2d",1)
	return img.get_clip(reg)

"""
# GOLDEN SEARCH CODE
def bracket(f,x1,h):
	c = 1.618033989 
	f1 = f(x1)
	x2 = x1 + h; f2 = f(x2)
	# Determine downhill direction and change sign of h if needed
	if f2 > f1:
		h = -h
		x2 = x1 + h; f2 = f(x2)
		# Check if minimum between x1 - h and x1 + h
		if f2 > f1: return x2,x1 - h 
	# Search loop
	for i in range (100):    
		h = c*h
		x3 = x2 + h; f3 = f(x3)
		if f3 > f2: return x1,x3
		x1 = x2; x2 = x3
		f1 = f2; f2 = f3
	print "Bracket did not find a mimimum"        
 
def goldsearch(f,a,b,tol=1.0e-9):
	from math import log, ceil
	nIter = int(ceil(-2.078087*log(tol/abs(b-a)))) # Eq. (10.4)
	R = 0.618033989
	C = 1.0 - R
	# First telescoping
	x1 = R*a + C*b; x2 = C*a + R*b
	f1 = f(x1); f2 = f(x2)
	# Main loop
	for i in range(nIter):
		if f1 > f2:
			a = x1
			x1 = x2; f1 = f2
			x2 = C*a + R*b; f2 = f(x2)
		else:
			b = x2
			x2 = x1; f2 = f1
			x1 = R*a + C*b; f1 = f(x1)
	if f1 < f2: return x1,f1
	else: return x2,f2
"""
