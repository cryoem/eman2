#
from __future__ import print_function
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

from builtins import range
from builtins import object
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
	for iy in range(m):
		wy = kb.sinhwin(iy-m//2)
		for ix in range(m):
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

def fpol(image, nnx, nny=1, nnz=1, RetReal = True, normalize = True):
	"""
		Name
			fpol -Interpolate image up by padding its Fourier transform with zeroes
		Input
			image: image to be interpolated.
			nnx: new nx dimension
			nny: new ny dimension (default = 0, meaning equal to the original ny)
			nnz: new nz dimension (default = 0, meaning equal to the original nz)
			RetReal:     Logical flag, if True, the returned image is real, if False, it is Fourier
			normalize:   Logical flag, if True, the returned image is normalized such that the power (sum of amplitudes squared) is preserved
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

	return  image.FourInterpol(nnx, nny, nnz, RetReal, normalize)

def fdecimate(image, nnx, nny=1, nnz=1, RetReal = True, normalize = True):
	"""
		Decimate image by truncating its Fourier transform
		Input
			image: image to be interpolated.
			nnx: new nx dimension
			nny: new ny dimension (default = 0, meaning equal to the original ny)
			nnz: new nz dimension (default = 0, meaning equal to the original nz)
			RetReal:     Logical flag, if True, the returned image is real, if False, it is Fourier
			normalize:   Logical flag, if True, the returned image is normalized such that the power (sum of amplitudes squared) is preserved
		Output
			the output decimated image
	"""
	return  image.FourTruncate(nnx, nny, nnz, RetReal, normalize)
    
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
	from fundamentals import smallprime
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


def prepf(image, npad = 2):
	"""
		Name
			prepf - prepare 2-D image for Fourier interpolation rotation/shift
		Input
			image: input image that is going to be rotated and shifted using fourier_rotate_shift2d
		Output
			imageft: Fourier space image prepared for Fourier interpolation rotation/shift
	"""
	
	cimage = Util.pad(image,2*(image.get_xsize()),2*(image.get_ysize()),1,0,0,0,"0.0")
	cimage.set_attr("npad",npad)
	cimage.div_sinc(1)
	cimage = cimage.norm_pad(False, 1)
	cimage.do_fft_inplace()
	cimage.center_origin_fft()
	cimage.fft_shuffle()
	cimage.set_attr_dict({"npad":npad,"prepf":1})
	return cimage

def prepi(image, RetReal = True):
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
	# first pad it with zeros in Fourier space
	q = image.FourInterpol(N, N, 1, 0)
	params = {"filter_type": Processor.fourier_filter_types.KAISER_SINH_INVERSE,
	          "alpha":alpha, "K":K, "r":r, "v":v, "N":N}
	q = Processor.EMFourierFilter(q, params)
	params = {"filter_type" : Processor.fourier_filter_types.TOP_HAT_LOW_PASS,
		"cutoff_abs" : 0.25, "dopad" : False}
	q = Processor.EMFourierFilter(q, params)
	if RetReal: return fft(q), kb
	else:  return q, kb


def prep_refim_gridding(refim, wr, numr, mode = "F"):
	from fundamentals import prepi
	nx = refim.get_xsize()
	ny = refim.get_ysize()
	cnx = nx//2+1
	cny = ny//2+1
	#precalculate rings
	temp,kb = prepi(refim)
	crefim = Util.Polar2Dmi(temp, cnx, cny, numr, mode, kb)
	Util.Normalize_ring(cimage, numr, 0 )
	Util.Frngs(crefim, numr)
	Util.Applyws(crefim, numr, wr)
	return  crefim,kb

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
	for i in range(n):
		tab.append(qt.get_value_at(i))
	return tab

def rot_avg_image(image_to_be_averaged):
	"""
	Rotational average
	Returns a 2-D or 3-D image containing a rotational average of image e
	"""
	import types
	from utilities import get_im
	if type(image_to_be_averaged) is bytes: image_to_be_averaged = get_im(image_to_be_averaged)
	return image_to_be_averaged.rotavg_i()

def ro_textfile(e, filename, helpful_string=""):
	"""Rotational average stored as a text file.
	   Saves a text file (suitable for gnuplot) of the rotational average of e.
	"""
	out = open(filename, "w")
	out.write("#Rotational average: %s\n" % (helpful_string));
	f = e.rotavg()
	nr = f.get_xsize()
	for ir in range(nr):
		out.write("%d\t%12.5g\n" % (ir, f.get_value_at(ir)))
	out.close()

def rops(e):
	"""Rotational average of the power spectrum.
	   Returns a 1-D image containing a rotational average
	   of the periodogram of image e.
		Input image can be real or Fourier, can be rectangular
		output length mapped onto x-dimension length
	"""
	from utilities import model_blank
	table = Util.rotavg_fourier(img)
	table = table[:len(table)//2]
	scale = (img.get_xsize() - 2*img.is_complex())*img.get_ysize()*img.get_zsize()
	scale = 4.0/scale/scale
	for i in range(len(table)): table[i] *= scale
	if lng:
		from math import log10
		for ir in range(1,nr): table[ir] = log10(table[ir])
		table[0] = table[1]
	ps = model_blank(len(table))
	for i in range(len(table)): ps[i] = table[i]
	return ps

def rops_textfile(e, filename, lng = False):
	"""Rotational average of the periodogram stored as a text file.
	   Saves a text file (suitable for gnuplot) of the rotational average 
	   of the periodogram of image e.
		Input image can be real or Fourier, can be rectangular
		output length mapped onto x-dimension length
	"""
	from utilities import write_text_file
	table = Util.rotavg_fourier(img)
	table = table[:len(table)//2]
	scale = (img.get_xsize() - 2*img.is_complex())*img.get_ysize()*img.get_zsize()
	scale = 4.0/scale/scale
	for i in range(len(table)): table[i] *= scale
	if lng:
		from math import log10
		for ir in range(1,nr): table[ir] = log10(table[ir])
		table[0] = table[1]
	write_text_file([list(range(nr)),table], filename)
	
def rops_table(img, lng = False):

	""" 
		Calculate 1D rotationally averaged 
		power spectrum and save it in list
		Input image can be real or Fourier, can be rectangular
		output length mapped onto x-dimension length
	"""
	table = Util.rotavg_fourier(img)
	table = table[:len(table)//2]
	scale = (img.get_xsize() - 2*img.is_complex())*img.get_ysize()*img.get_zsize()
	scale = 4.0/scale/scale
	for i in range(len(table)): table[i] *= scale
	if lng:
		from math import log10
		for ir in range(1,nr): table[ir] = log10(table[ir])
		table[0] = table[1]
	return table

'''
It is not used anywhere, so I commented it out  02/03/2019 PAP
def rops_dir(indir, output_dir = "1dpw2_dir"):
	"""
		Calculate 1D rotationally averaged power spectra from
		image stack listed in a directory
	"""
	from utilities import get_im, write_text_file
	import os
	flist = os.listdir(indir)
	if os.path.exists(output_dir) is False: os.mkdir(output_dir)
	for i, v in enumerate(flist):
		(filename, filextension) = os.path.splitext(v)
		nima = EMUtil.get_image_count(os.path.join(indir,v))
		for im in range(nima):
			e = get_im(os.path.join(indir,v), im)
			temp = Util.rotavg_fourier(img)
			temp = table[:len(temp)//2]
			if im == 0:
				table= temp[:]
				scale = (img.get_xsize() - 2*img.is_complex())*img.get_ysize()*img.get_zsize()
				scale = 4.0/scale/scale
			else :  for i in range(len(table)): table[i] += temp[i]
		for i in range(len(table)): table[i] *= scale/nima
		write_text_file(table, os.path.join(output_dir, "1dpw2_"+filename+".txt"))
'''

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

'''
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
'''

def gridrot_shift2D(image, ang = 0.0, sx = 0.0, sy = 0.0, scale = 1.0):
	"""
		Rotate and shift an image using gridding in Fourier space.
	"""
	from EMAN2 import Processor
	from fundamentals import fftip, fft
	from utilities import compose_transform2

	nx = image.get_xsize()
	if( (ang != 0.0) or (scale != 1.0) ):
		image1 = image.copy()  # This step is needed, otherwise image will be changed outside the function
		# prepare 
		npad = 2
		N = nx*npad
		K = 6
		alpha = 1.75
		r = nx/2
		v = K/2.0/N
		kb = Util.KaiserBessel(alpha, K, r, v, N)

		# invert shifts
		_, tsx, tsy, _ = compose_transform2(0.,sx,sy,1,-ang,0,0,1)
		# split shift into integer and fractional parts
		isx = int(tsx)
		fsx = tsx - isx
		isy = int(tsy)
		fsy = tsy - isy
		# shift image to the center
		Util.cyclicshift(image1,{"dx":isx,"dy":isy,"dz":0})
		# bring back fractional shifts
		_, tsx, tsy, _ = compose_transform2(0.,fsx,fsy,1,ang,0,0,1)
		# divide out gridding weights
		image1.divkbsinh(kb)
		# pad and center image, then FFT
		image1 = image1.norm_pad(False, npad)
		fftip(image1)
		# Put the origin of the (real-space) image at the center
		image1.center_origin_fft()
		# gridding rotation
		image1 = image1.fouriergridrot2d(ang, scale, kb)
		if( (fsx != 0.0) or (fsy != 0.0) ):
			params = {"filter_type" : Processor.fourier_filter_types.SHIFT,	"x_shift" : float(tsx), "y_shift" : float(tsy), "z_shift" : 0.0 }
			image1 = Processor.EMFourierFilter(image1, params)
		# put the origin back in the corner
		image1.center_origin_fft()
		# undo FFT and remove padding (window)
		image1 = fft(image1)
		image1 = image1.window_center(nx)
	else:
		if( (sx != 0.0) or (sy != 0.0) ):
			params = {"filter_type" : Processor.fourier_filter_types.SHIFT,	"x_shift" : float(sx), "y_shift" : float(sy), "z_shift" : 0.0 }
			image1 = Processor.EMFourierFilter(image, params)
		else: return image
		
	return image1


def ft2polargrid(image, ring_length, nb, ne):
	"""
		resample to polar coordinates using gridding in Fourier space.
	"""
	from EMAN2 import Processor
	from fundamentals import fftip, fft

	nx = image.get_xsize()
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
	return image1.ft2polargrid(ring_length, nb, ne, kb)


def rot_shift2D(img, angle = 0.0, sx = 0.0, sy = 0.0, mirror = 0, scale = 1.0, interpolation_method = None, mode = "background"):
	"""
		rotate/shift image using (interpolation_method):
		1. linear    interpolation
		2. quadratic interpolation
		3. gridding
		4. ftgridding
		5. fourier
		mode specifies what will be put in corners, should they stick out:
		background - leave original values
		cyclic - use pixels from the image using wrap-around transformation
		
	"""

	if scale == 0.0 :  ERROR("0 scale factor encountered","rot_shift2D", 1)
	if(interpolation_method):  use_method = interpolation_method
	else:  use_method = interpolation_method_2D

	if(use_method == "linear" and mode == "cyclic"):
		T  = Transform({'type': 'SPIDER', 'psi': angle, 'tx': sx, 'ty': sy, 'scale':scale})
		img = img.rot_scale_trans(T, None)
		if  mirror: img.process_inplace("xform.mirror", {"axis":'x'})
		return img
	elif(use_method == "linear" and mode == "background"):
		T  = Transform({'type': 'SPIDER', 'psi': angle, 'tx': sx, 'ty': sy, 'scale':scale})
		img = img.rot_scale_trans_background(T, None)
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
	elif(use_method == "fourier"):
		if img.is_complex() :
			if img.get_attr_default("prepf",0) != 1:  ERROR("Incorrect input image Fourier format","rot_shift2D fourier method",1)
		img = img.fourier_rotate_shift2d(angle, sx, sy, 2)
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
	if (mode == "cyclic"): return image.rot_scale_trans(T, None)
	else: return image.rot_scale_trans_background(T, None)


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
	primelist = [2,3,5,7,11,13,17,19,23,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97]
	lip = min(numprime,len(primelist))
	for i in range(1, arbit_num+1):
		x62 = arbit_num-i+1
		for k in range(1,arbit_num+1): # fake loop try to divide the arbit_num
			for j in range(0,lip):
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

def sinc2inv(nx):
	from math import sqrt
	s = sincinv(nx)
	return [i*i for i in s]

def sincinv(nx):
	from math import pi,sin
	cdf =pi/nx
	npad = 1
	nxb = nx/2/npad
	nxe = nxb + (nx/npad)%2
	s = [1.0]*nx
	for i in range( -nxb, nxe):
		if( i != 0 ):
			rrr=abs(i)
			s[i+nxb] = (rrr*cdf)/sin(rrr*cdf)
	return s

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
	for iy in range(1, x29+1):	
		x21 = (win_size/x39)*(iy-1) + edge_y  #  y-direction it should start from 0 if edge_y=0	      
		for ix in  range(1, x26+1):			 
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
	for ix in  range(1, num_bnd+1):
		x22 = (win_x)*(ix-1)# x-direction it should start from 0 if edge_x =0
		iz=0
		for iy in range(1, x29+1):	
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
		Calculate set of periodograms for tiles.  Returns a list.
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
	for iy in range(1, x29+1):	
		x21 = (win_size/x39)*(iy-1) + edge_y  #  y-direction it should start from 0 if edge_y=0	      
		for ix in  range(1, x26+1):			 
			x22 = (win_size/x38)*(ix-1) + edge_x  # x-direction it should start from 0 if edge_x =0
			wi  = ramp( window2d(e_fil, win_size, win_size, "l", x22, x21) )
			st = Util.infomask(wi, None, True)
			wi = (wi - st[0])/st[1]*win_size
			pw2.append(periodogram(wi))
	return  pw2


def window2d(img, isize_x, isize_y, opt="c", ix=0, iy=0):
	"""
		Three ways of windowing out a portion of an image:
		1. "c" Get the central part: "c" ( default setting )
		2. "l" Get clip starts from the top left corner
		3. "a" Get clip with arbitrary point (ix, iy) used as the input image center point 
				in case 1, 0,0 corresponds to image center
				in case 2, nx//2,ny//2 corresponds to image center
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
	print("Bracket did not find a mimimum")        
 
def goldsearch(f,a,b,tol=1.0e-9):
	from math import log, ceil
	nIter = int(ceil(-2.078087*log(tol/abs(b-a))))
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


#  01/06/2016 - This is my recoding of old FORTRAN code with the hope that python's double precission
#                 will fix the problem of rotation of a 0,0,0 direction.  It does not as one neeeds psi
#                 in this case as well.  So, the only choice is to use small theta instead of exact 0,0,0 direction
def rotate_params(params, transf):
	from fundamentals import rotmatrix, recmat, rotmatrix
	matinv = rotmatrix( -transf[2], -transf[1], -transf[0] )
	n = len(params)
	cpar = [None]*n
	for i in range(n):
		d = rotmatrix( params[i][0], params[i][1], params[i][2] )
		phi, theta, psi = recmat(mulmat(d,matinv))
		cpar[i] = [phi, theta, psi]
	return cpar

"""
def rotmatrix(phi,theta,psi):
	from math import sin,cos,radians
	rphi   = radians(phi)
	rtheta = radians(theta)
	rpsi   = radians(psi)
	mat = [[0.0]*3,[0.0]*3,[0.0]*3]

	mat[0][0] =  cos(rpsi)*cos(rtheta)*cos(rphi) - sin(rpsi)*sin(rphi)
	mat[1][0] = -sin(rpsi)*cos(rtheta)*cos(rphi) - cos(rpsi)*sin(rphi)
	mat[2][0] =            sin(rtheta)*cos(rphi)


	mat[0][1] =  cos(rpsi)*cos(rtheta)*sin(rphi) + sin(rpsi)*cos(rphi)
	mat[1][1] = -sin(rpsi)*cos(rtheta)*sin(rphi) + cos(rpsi)*cos(rphi)
	mat[2][1] =            sin(rtheta)*sin(rphi)


	mat[0][2] = -cos(rpsi)*sin(rtheta)
	mat[1][2] =  sin(rpsi)*sin(rtheta)
	mat[2][2] =            cos(rtheta)
	return mat
"""

def rotmatrix(phi,theta,psi):
	from math import sin,cos,radians
	rphi   = radians(phi)
	rtheta = radians(theta)
	rpsi   = radians(psi)
	cosphi = cos(rphi)
	sinphi = sin(rphi)
	costheta = cos(rtheta)
	sintheta = sin(rtheta)
	cospsi = cos(rpsi)
	sinpsi = sin(rpsi)
	mat = [[0.0]*3,[0.0]*3,[0.0]*3]

	mat[0][0] =  cospsi*costheta*cosphi - sinpsi*sinphi
	mat[1][0] = -sinpsi*costheta*cosphi - cospsi*sinphi
	mat[2][0] =            sintheta*cosphi


	mat[0][1] =  cospsi*costheta*sinphi + sinpsi*cosphi
	mat[1][1] = -sinpsi*costheta*sinphi + cospsi*cosphi
	mat[2][1] =            sintheta*sinphi


	mat[0][2] = -cospsi*sintheta
	mat[1][2] =  sinpsi*sintheta
	mat[2][2] =            costheta
	return mat

def mulmat(m1,m2):
	mat = [[0.0]*3,[0.0]*3,[0.0]*3]
	"""
	for i in range(3):
		for j in range(3):
			for k in range(3):
				mat[i][j] += m1[i][k]*m2[k][j]
			#mat[i][j] = round(mat[i][j],8)
	"""
	mat[0][0] = m1[0][0]*m2[0][0] + m1[0][1]*m2[1][0] + m1[0][2]*m2[2][0]
	mat[0][1] = m1[0][0]*m2[0][1] + m1[0][1]*m2[1][1] + m1[0][2]*m2[2][1]
	mat[0][2] = m1[0][0]*m2[0][2] + m1[0][1]*m2[1][2] + m1[0][2]*m2[2][2]
	mat[1][0] = m1[1][0]*m2[0][0] + m1[1][1]*m2[1][0] + m1[1][2]*m2[2][0]
	mat[1][1] = m1[1][0]*m2[0][1] + m1[1][1]*m2[1][1] + m1[1][2]*m2[2][1]
	mat[1][2] = m1[1][0]*m2[0][2] + m1[1][1]*m2[1][2] + m1[1][2]*m2[2][2]
	mat[2][0] = m1[2][0]*m2[0][0] + m1[2][1]*m2[1][0] + m1[2][2]*m2[2][0]
	mat[2][1] = m1[2][0]*m2[0][1] + m1[2][1]*m2[1][1] + m1[2][2]*m2[2][1]
	mat[2][2] = m1[2][0]*m2[0][2] + m1[2][1]*m2[1][2] + m1[2][2]*m2[2][2]

	return mat

def recmat(mat):
	from math import acos,asin,atan2,degrees,pi
	def sign(x):
		if( x >= 0.0 ): return 1
		else:  return -1
	"""
	mat = [[0.0]*3,[0.0]*3,[0.0]*3]
	# limit precision
	for i in range(3):
		for j in range(3):
			mat[i][j] = inmat[i][j]
			#if(abs(inmat[i][j])<1.0e-8):  mat[i][j] = 0.0
			#else: mat[i][j] = inmat[i][j]
	for i in range(3):
		for j in range(3):  print  "     %14.8f"%mat[i][j],
		print ""
	"""
	if(mat[2][2] == 1.0):
		theta = 0.0
		psi = 0.0
		if( mat[0][0] == 0.0 ):
			phi = asin(mat[0][1])
		else:
			phi = atan2(mat[0][1],mat[0][0])
	elif(mat[2][2] == -1.0):
		theta = pi
		psi = 0.0
		if(mat[0][0] == 0.0):
			phi = asin(-mat[0][1])
		else:
			phi = atan2(-mat[0][1],-mat[0][0])
	else:
		theta = acos(mat[2][2])
		st = sign(theta)
		#print theta,st,mat[2][0]
		if(mat[2][0] == 0.0):
			if( st != sign(mat[2][1]) ):
				phi = 1.5*pi
			else:
				phi = 0.5*pi
		else:
			phi = atan2(st*mat[2][1], st*mat[2][0])

		#print theta,st,mat[0][2],mat[1][2]
		if(mat[0][2] == 0.0):
			if( st != sign(mat[1][2]) ):
				psi = 1.5*pi
			else:
				psi = 0.5*pi
		else:
			psi = atan2(st*mat[1][2], -st*mat[0][2])
	#pi2 = 2*pi
	#return  degrees(round(phi%pi2,8)),degrees(round(theta%pi2,8)),degrees(round(psi%pi2,8))
	#return  degrees(round(phi,10)%pi2)%360.0,degrees(round(theta,10)%pi2)%360.0,degrees(round(psi,10)%pi2)%360.0
	return  degrees(phi)%360.0,degrees(theta)%360.0,degrees(psi)%360.0


"""

def mulmat_np(m1,m2):
	import numpy as np
	mat1 = np.matrix(m1,dtype="f8")
	mat2 = np.matrix(m2,dtype="f8")
	mat1 = np.array(mat1*mat2)
	return [list(q) for q in mat1]

def rotmatrix_np(phi,theta,psi):
	import numpy as np
	mat = np.matrix(((0.,0.,0.),(0.,0.,0.),(0.,0.,0.)), dtype = "f8")
	rphi   = np.radians(np.float64(phi))
	rtheta = np.radians(np.float64(theta))
	rpsi   = np.radians(np.float64(psi))
	cosphi = np.cos(rphi)
	sinphi = np.sin(rphi)
	costheta = np.cos(rtheta)
	sintheta = np.sin(rtheta)
	cospsi = np.cos(rpsi)
	sinpsi = np.sin(rpsi)

	mat[0,0] =  cospsi*costheta*cosphi - sinpsi*sinphi
	mat[1,0] = -sinpsi*costheta*cosphi - cospsi*sinphi
	mat[2,0] =            sintheta*cosphi


	mat[0,1] =  cospsi*costheta*sinphi + sinpsi*cosphi
	mat[1,1] = -sinpsi*costheta*sinphi + cospsi*cosphi
	mat[2,1] =            sintheta*sinphi


	mat[0,2] = -cospsi*sintheta
	mat[1,2] =  sinpsi*sintheta
	mat[2,2] =            costheta
	return mat


def recmat_np(mat):
	#from math import np.arccos,np.np.arcsin,np.arctan2,degrees,pi
	import numpy as np
	'''
	def sign(x):
		if( x >= 0.0 ): return 1
		else:  return -1
	mat = [[0.0]*3,[0.0]*3,[0.0]*3]
	# limit precision
	for i in range(3):
		for j in range(3):
			mat[i,j] = inmat[i,j]
			#if(abs(inmat[i,j])<1.0e-8):  mat[i,j] = 0.0
			#else: mat[i,j] = inmat[i,j]
	for i in range(3):
		for j in range(3):  print  "     %14.8f"%mat[i,j],
		print ""
	'''
	if(mat[2,2] == 1.0):
		theta = 0.0
		psi = 0.0
		if( mat[0,0] == 0.0 ):
			phi = np.np.arcsin(mat[0,1])
		else:
			phi = np.arctan2(mat[0,1],mat[0,0])
	elif(mat[2,2] == -1.0):
		theta = pi
		psi = 0.0
		if(mat[0,0] == 0.0):
			phi = np.np.arcsin(-mat[0,1])
		else:
			phi = np.arctan2(-mat[0,1],-mat[0,0])
	else:
		theta = np.arccos(mat[2,2])
		st = np.sign(theta)
		#print theta,st,mat[2,0]
		if(mat[2,0] == 0.0):
			if( st != np.sign(mat[2,1]) ):
				phi = 1.5*pi
			else:
				phi = 0.5*pi
		else:
			phi = np.arctan2(st*mat[2,1], st*mat[2,0])

		#print theta,st,mat[0,2],mat[1,2]
		if(mat[0,2] == 0.0):
			if( st != np.sign(mat[1,2]) ):
				psi = 1.5*pi
			else:
				psi = 0.5*pi
		else:
			psi = np.arctan2(st*mat[1,2], -st*mat[0,2])
	pi2 = 2*np.pi
	#return  degrees(round(phi,10)%pi2),degrees(round(theta,10)%pi2),degrees(round(psi,10)%pi2)
	#return  degrees(round(phi,10)%pi2)%360.0,degrees(round(theta,10)%pi2)%360.0,degrees(round(psi,10)%pi2)%360.0
	return  np.degrees(np.mod(phi,pi2)),np.degrees(np.mod(theta,pi2)),np.degrees(np.mod(psi,pi2))
"""

class symclass(object):
	import numpy as np
	def __init__(self, sym):
		"""
		  sym: cn, dn, oct, tet, icos
		"""
		from math import degrees, radians, sin, cos, tan, atan, acos, sqrt, pi
		#from utilities import get_sym, get_symt
		from string import lower
		self.sym = sym.lower()
		if(self.sym[0] == "c"):
			self.nsym = int(self.sym[1:])
			if(self.nsym<1):  ERROR("For Cn symmetry, we need n>0","symclass",1)
			self.brackets = [[360./self.nsym,90.0,360./self.nsym,90.0],[360./self.nsym,180.0,360./self.nsym,180.0]]
			self.symangles = []
			for i in range(self.nsym):
				self.symangles.append([0.0, 0.0, i*360./self.nsym])

		elif(self.sym[0] == "d"):
			self.nsym = 2*int(self.sym[1:])
			if(self.nsym<1):  ERROR("For Dn symmetry, we need n>0","symclass",1)
			self.brackets = [[360./self.nsym,90.0,360./self.nsym,90.0],[360./self.nsym*2,90.0,360./self.nsym*2,90.0]]
			self.symangles = []
			for i in range(self.nsym/2):
				self.symangles.append([0.0, 0.0, i*360./self.nsym*2])
			for i in range(self.nsym/2):
				self.symangles.append([0.0, 180.0, (i*360./self.nsym*2+180.0*(int(self.sym[1:])%2))%360.0])

		elif(self.sym[:3] == "oct"):
			self.nsym = 24
			ncap = 4
			cap_sig = 360.0/ncap  # also called platonic_params["az_max"]
			alpha = degrees(acos(1.0/(sqrt(3.0)*tan(2*pi/ncap/2.0)))) # also platonic_params["alt_max"]
			theta = degrees(0.5*acos( cos(radians(cap_sig))/(1.0-cos(radians(cap_sig))) ))  #  also platonic_params["theta_c_on_two"]
			self.brackets = [[180./ncap,theta,cap_sig,alpha],[360./ncap,theta,cap_sig,alpha]]
			self.symangles = [[0.0,0.0,float(i)] for i in range(0,271,90)]
			for i in range(0,271,90):
				for j in range(0,271,90):
					self.symangles.append([float(j),90.0,float(i)])
			for i in range(0,271,90):  self.symangles.append([0.0,180.0,float(i)])

		elif(self.sym[:3] == "tet"):
			self.nsym = 12
			ncap = 3
			cap_sig = 360.0/ncap  # also called platonic_params["az_max"]
			alpha = degrees(acos(1.0/(sqrt(3.0)*tan(2*pi/ncap/2.0)))) # also platonic_params["alt_max"]
			theta = degrees(0.5*acos( cos(radians(cap_sig))/(1.0-cos(radians(cap_sig))) ))  #  also platonic_params["theta_c_on_two"]
			self.brackets = [[360.0/ncap,theta,cap_sig,alpha],[360.0/ncap,theta,cap_sig,alpha]]
			lvl1 = degrees(acos(-1.0/3.0)) # There  are 3 faces at this angle
			self.symangles = [ [0.,0.,0.], [0., 0., 120.], [0., 0., 240.]]
			for l1 in range(0,241,120):
				for l2 in range(60,301,120):
					self.symangles.append([float(l1),lvl1,float(l2)])
			
			"""
			#  These angles were translated from eman to spider, but the do not agree with definitions of subunit above
			for l1 in range(30,271,120):
				for l2 in range(30,271,120):
					self.symangles.append([float(l1),lvl1,float(l2)])
			"""

		elif(self.sym[:4] == "icos"):
			self.nsym = 60
			ncap = 5
			cap_sig = 360.0/ncap  # also called platonic_params["az_max"]
			alpha = degrees(acos(1.0/(sqrt(3.0)*tan(2*pi/ncap/2.0)))) # also platonic_params["alt_max"]
			theta = degrees(0.5*acos( cos(radians(cap_sig))/(1.0-cos(radians(cap_sig))) ))  #  also platonic_params["theta_c_on_two"]
			self.brackets = [[36.,theta,cap_sig,alpha],[72.,theta,cap_sig,alpha]]
			lvl1= degrees(atan(2.0))  #there are 5 pentagons with centers at this height (angle)
			lvl2 = 180.0 - lvl1      #there are 5 pentagons with centers at this height (angle)
			self.symangles = [[0.0,0.0,float(i)] for i in range(0,288+1,72)]
			for l1 in range(0,288+1,72):
				for l2 in range(36,324+1,72):
					self.symangles.append([float(l1),lvl1,float(l2)])
			for l1 in range(36,324+1,72):
				for l2 in range(0,288+1,72):
					self.symangles.append([float(l1),lvl2,float(l2)])
			for i in range(0,288+1,72):  self.symangles.append([0.0,180.0,float(i)])
		
		else:  ERROR("Unknown symmetry","symclass",1)

		#
		self.transform = []
		self.symatrix = []
		for args in self.symangles:
			self.transform.append(Transform({"type":"spider", "phi":args[0], "theta":args[1], "psi":args[2]}))
			self.symatrix.append(rotmatrix(args[0],args[1],args[2]))

	def is_in_subunit(self, phi, theta, inc_mirror=1):
		"""
		Input: a projection direction specified by (phi, theta).
				inc_mirror = 1 consider mirror directions as unique
				inc_mirror = 0 consider mirror directions as outside of unique range.
		Output: True if input projection direction is in the first asymmetric subunit,
		        False otherwise.
		"""
		from math import degrees, radians, sin, cos, tan, atan, acos, sqrt
		if( (self.sym[0] == "c")  or  (self.sym[0] == "d" and (self.nsym//2)%2 == 0) ):
			if((phi>= 0.0 and phi<self.brackets[inc_mirror][0]) and (theta<=self.brackets[inc_mirror][1])):  return True
			else:  return False
		elif( self.sym[0] == "d" and (self.nsym//2)%2 == 1 ):
			if(theta<=self.brackets[inc_mirror][1]):
				phib = 360.0/self.nsym
				if( phi>=0.0 and phi<self.brackets[1][0] ):
					if(inc_mirror==1):  return True
					elif( (phi>= 0.0 and phi<phib/2) or (phi>= phib and phi<(phib+phib/2)) ): return True
			return False
			
		elif( (self.sym[:3] == "oct")  or  (self.sym[:4] == "icos") ):
			if( phi>= 0.0 and phi<self.brackets[inc_mirror][0] and theta<=self.brackets[inc_mirror][3] ):
				tmphi = min(phi, self.brackets[inc_mirror][2]-phi)
				baldwin_lower_alt_bound = \
				(sin(radians(self.brackets[inc_mirror][2]/2.0-tmphi))/tan(radians(self.brackets[inc_mirror][1])) + \
					sin(radians(tmphi))/tan(radians(self.brackets[inc_mirror][3])))/sin(radians(self.brackets[inc_mirror][2]/2.0))
				baldwin_lower_alt_bound = degrees(atan(1.0/baldwin_lower_alt_bound))
				#print(  "  baldwin_lower_alt_bound ",self.brackets,baldwin_lower_alt_bound,theta)
				if(baldwin_lower_alt_bound>theta): return True
				else: return False
			else:
				#print "phi",self.brackets
				return False
		elif( self.sym[:3] == "tet" ):
			if( phi>= 0.0 and phi<self.brackets[inc_mirror][0] and theta<=self.brackets[inc_mirror][3] ):
				tmphi = min(phi, self.brackets[inc_mirror][2]-phi)
				baldwin_lower_alt_bound = \
				(sin(radians(self.brackets[inc_mirror][2]/2.0-tmphi))/tan(radians(self.brackets[inc_mirror][1])) + \
					sin(radians(tmphi))/tan(radians(self.brackets[inc_mirror][3])))/sin(radians(self.brackets[inc_mirror][2]/2.0))
				baldwin_lower_alt_bound = degrees(atan(1.0/baldwin_lower_alt_bound))
				#print(  "  baldwin_lower_alt_bound ",phi,theta,baldwin_lower_alt_bound,self.brackets[inc_mirror])
				if(baldwin_lower_alt_bound>theta):
					if( inc_mirror == 1 ):
						return True
					else:
						baldwin_upper_alt_bound = \
						(sin(radians(self.brackets[inc_mirror][2]/2.0-tmphi))/tan(radians(self.brackets[inc_mirror][1])) + \
							sin(radians(tmphi))/tan(radians(self.brackets[inc_mirror][3]/2.0)))/sin(radians(self.brackets[inc_mirror][2]/2.0))
						baldwin_upper_alt_bound = degrees(atan(1.0/baldwin_upper_alt_bound))
						#print(  "  baldwin_upper_alt_bound ",phi,theta,baldwin_upper_alt_bound,self.brackets[inc_mirror])
						if(baldwin_upper_alt_bound<theta): return False
						else:  return True
				else: return False
			else:
				#print "phi",self.brackets
				return False
		else:  ERROR("unknown symmetry","symclass: is_in_subunit",1)

	def symmetry_related(self, angles):
		"""
		Generate all symmetry related instances of input angles
		Input:  is a list of list of n triplets  
			[[phi0,theta0,psi0],[phi1,theta1,psi1],...]
		Output: is a list of list triplets whose length is k times larger than that of input,
		        where k is symmetry multiplicity.  First k are symmetry related version of the first triplet, 
		        second k are symmetry related version of the second triplet and so on...
		   [[phi0,theta0,0],[phi0,theta0,0]_SYM0,...,[phi1,theta1,],[phi1,theta1,]_SYM1,...]
		"""
		redang = [angles[:]]
		if(self.sym[0] == "c"):
			qt = 360.0/self.nsym
			for l in range(1,self.nsym):
				redang.append([(angles[0]+l*qt)%360.0, angles[1], angles[2]])
		elif(self.sym[0] == "d"):
			nsm = self.nsym/2
			qt = 360.0/nsm
			for l in range(1,nsm):
				redang.append([(angles[0]+l*qt)%360.0, angles[1], angles[2]])
			for l in range(nsm,self.nsym):
				redang.append([(360.0-redang[l-nsm][0])%360.0, 180.0-angles[1], (angles[2]+180.0*(nsm%2))%360.0])
		else:
			from fundamentals import rotmatrix, recmat, mulmat
			mat = rotmatrix(angles[0],angles[1],angles[2])
			for l in range(1,self.nsym):
				p1,p2,p3 = recmat( mulmat( mat , self.symatrix[l]) )
				redang.append([p1,p2,p3])

		return redang


	def symmetry_neighbors(self, angles):
		"""
		Generate symmetry related instances of input angles only for asymmetric regions adjacent to the zero's one
		  input is a list of lists  [[phi0,theta0,0],[phi1,theta1,0],...]
		Output: is a list of list triplets whose length is k times larger than that of input,
		        where k is the number of adjacent regions for a given point group symmetry.  
		        First k are symmetry related version of the first triplet, 
		        second k are symmetry related version of the second triplet and so on...
		  output is [[phi0,theta0,0],[phi0,theta0,0]_SYM1,...,[phi1,theta1,],[phi1,theta1,]_SYM1,...]
		"""
		if( self.sym[0] == "c" or self.sym[0] == "d" ):
			temp = Util.symmetry_neighbors(angles, self.sym)
			nt = len(temp)/3
			return [[temp[l*3],temp[l*3+1],0.0] for l in range(nt) ]
		#  Note symmetry neighbors below refer to the particular order 
		#   in which this class generates symmetry matrices
		neighbors = {}
		neighbors["oct"]  = [1,2,3,8,9,12,13]
		neighbors["tet"]  = [1,2,3,4,6,7]
		neighbors["icos"] = [1,2,3,4,6,7,11,12]
		sang = [[] for l in range(len(angles)*(len(neighbors[self.sym])+1))]
		for i,q in enumerate(angles):  sang[i*(len(neighbors[self.sym])+1)] = angles[i][:]
		from fundamentals import rotmatrix, recmat, mulmat
		for i,q in enumerate(angles):
			mat = rotmatrix(q[0],q[1],q[2])
			for j,l in enumerate(neighbors[self.sym]):
				p1,p2,p3 = recmat( mulmat( mat , self.symatrix[l]) )
				sang[i*(len(neighbors[self.sym])+1) + (j+1)] = [p1,p2,0.0]
		return sang

	def reduce_anglesets(self, angles, inc_mirror=1):
		"""
		  Input is either list ot lists [[phi,thet,psi],[],[]] or a triplet [phi,thet,psi]
				inc_mirror = 1 consider mirror directions as unique
				inc_mirror = 0 consider mirror directions as outside of unique range.
		  It will map all triplets to the first asymmetric subunit.
		"""
		from math import degrees, radians, sin, cos, tan, atan, acos, sqrt
		import types
		is_platonic_sym = self.sym[0] == "o" or self.sym[0] == "i"
		if(self.sym[0] == "c"): qs = 360.0/self.nsym
		elif(self.sym[0] == "d"): qs = 720.0/self.nsym
		if( type(angles[0]) is list):
			toprocess = angles
			lis = True
		else:
			toprocess = [angles]
			lis = False
		redang = []
		for q in toprocess:
			phi = q[0]; theta = q[1]; psi = q[2]
			if is_platonic_sym:
				if(not self.is_in_subunit(phi, theta, 1)):
					mat = rotmatrix(phi,theta,psi)
					for l in range(self.nsym):
						p1,p2,p3 = recmat( mulmat( mat , self.symatrix[l]) )
						#print(p1,p2,p3)
						if(self.is_in_subunit(p1, p2, 1)):
							phi=p1; theta=p2; psi=p3
							#print("  FOUND ")
							break
					if( inc_mirror == 0 ):
						if(phi>=self.brackets[0][0]):
							phi = self.brackets[1][0]-phi
							if(l>0): psi = (360.0-psi)%360.0
				elif(inc_mirror == 0):
					if(phi>=self.brackets[0][0]):
						phi = self.brackets[1][0]-phi
						psi = (360.0-psi)%360.0
					elif(inc_mirror == 0):
						if(phi>=self.brackets[0][0]):  phi = self.brackets[1][0]-phi
			elif( self.sym[0] == "t" ):
				if(not self.is_in_subunit(phi, theta, inc_mirror)):
					mat = rotmatrix(phi,theta,psi)
					fifi = False
					for l in range(self.nsym):
						p1,p2,p3 = recmat( mulmat( mat , self.symatrix[l]) )
						if(self.is_in_subunit(p1, p2, inc_mirror)):
							if(inc_mirror):
								phi=p1; theta=p2; psi=p3
								fifi = True
								break
							else:
								if(self.is_in_subunit(p1, p2, 0)):
									phi=p1; theta=p2; psi=p3
									fifi = True
									break
					if(inc_mirror == 1 and not fifi): print("  FAILED no mirror ")
					if( not fifi ):
						phi = (180.0+phi)%360.0; theta = 180.0 - theta; psi = (180.0 - psi)%360.0
						mat = rotmatrix(phi,theta,psi)
						for l in range(self.nsym):
							p1,p2,p3 = recmat( mulmat( mat , self.symatrix[l]) )
							if(self.is_in_subunit(p1, p2, 0)):
								phi=p1; theta=p2; psi=p3
								fifi = True
								break

					if( not fifi ):  print("  FAILED mirror ")
			else:
				if( theta>90.0  and inc_mirror == 0 ):
					phi = (180.0+phi)%360.0; theta = 180.0 - theta; psi = (180.0 - psi)%360.0
				phi = phi%qs
				if(self.sym[0] == "d"):
					if( inc_mirror == 0 ):
						if((self.nsym//2)%2 == 0):
							if(phi>=qs/2):
								phi = qs-phi
								psi = (360.0-psi)%360.0
						else:
							if(phi>=360.0/self.nsym/2 and phi<360.0/self.nsym):
								phi = 360.0/self.nsym-phi
								psi = 360.0 - psi
							elif(phi>=360.0/self.nsym+360.0/self.nsym/2 and phi<720.0/self.nsym):
								phi = 720.0/self.nsym-phi+360.0/self.nsym
								psi = (360.0-psi)%360.0
					else:
						if(theta>90.0):
							mat = rotmatrix(phi,theta,psi)
							for l in range(self.nsym):
								p1,p2,p3 = recmat( mulmat( mat , self.symatrix[l]) )
								#print(p1,p2,p3)
								if(self.is_in_subunit(p1, p2, 1)):
									phi=p1; theta=p2; psi=p3
									#print("  FOUND ")
									break

			redang.append([phi, theta, psi])

		if lis: return redang
		else: return redang[0]

	def reduce_angles(self, phiin, thetain, psiin, inc_mirror=1):
		"""
		  It will map three input angles to the first asymmetric subunit.
				inc_mirror = 1 consider mirror directions as unique
				inc_mirror = 0 consider mirror directions as outside of unique range.
		"""
		from math import degrees, radians, sin, cos, tan, atan, acos, sqrt
		is_platonic_sym = self.sym[0] == "o" or self.sym[0] == "i"
		if(self.sym[0] == "c"): qs = 360.0/self.nsym
		elif(self.sym[0] == "d"): qs = 720.0/self.nsym
		if is_platonic_sym:
			if(not self.is_in_subunit(phiin, thetain, 1)):
				mat = rotmatrix(phiin,thetain,psiin)
				phi=phiin; theta=thetain; psi=psiin
				for l in range(self.nsym):
					p1,p2,p3 = recmat( mulmat( mat , self.symatrix[l]) )
					#print(p1,p2,p3)
					if(self.is_in_subunit(p1, p2, 1)):
						phi=p1; theta=p2; psi=p3
						#print("  FOUND ")
						break
				if( inc_mirror == 0 ):
					if(phi>=self.brackets[0][0]):
						phi = self.brackets[1][0]-phi
						if(l>0): psi = (360.0-psi)%360.0
			elif(inc_mirror == 0):
				phi = phiin; theta = thetain; psi = psiin
				if(phi>=self.brackets[0][0]):
					phi = self.brackets[1][0]-phi
					psi = (360.0-psi)%360.0
			else:
				phi = phiin; theta = thetain; psi = psiin
		elif( self.sym[0] == "t" ):
			phi=phiin; theta=thetain; psi=psiin
			if(not self.is_in_subunit(phi, theta, inc_mirror)):
				mat = rotmatrix(phi,theta,psi)
				for l in range(self.nsym):
					p1,p2,p3 = recmat( mulmat( mat , self.symatrix[l]) )
					#print(" ASYM  ",l,p1,p2,p3 )
					if(self.is_in_subunit(p1, p2, inc_mirror)):
						if(inc_mirror):  return p1,p2,p3
						else:
							if(self.is_in_subunit(p1, p2, 0)):  return p1,p2,p3
				if(inc_mirror == 1): print("  FAILED no mirror ")
				phi = (180.0+phi)%360.0; theta = 180.0 - theta; psi = (180.0 - psi)%360.0
				mat = rotmatrix(phi,theta,psi)
				for l in range(self.nsym):
					p1,p2,p3 = recmat( mulmat( mat , self.symatrix[l]) )
					#print(" MIR  ",l,p1,p2,p3 )
					if(self.is_in_subunit(p1, p2, 0)):  return p1,p2,p3

				print("  FAILED mirror ")
		else:
			if( thetain>90.0 and inc_mirror == 0 ):
				phi = (180.0+phiin)%360.0; theta = 180.0 - thetain; psi = (180.0 - psiin)%360.0
			else:
				phi = phiin; theta = thetain; psi = psiin
			phi = phi%qs
			if(self.sym[0] == "d"):
				if( inc_mirror == 0 ):
					if((self.nsym//2)%2 == 0):
						if(phi>=qs/2):
							phi = qs-phi
							psi = 360.0 - psi
					else:
						if(phi>=360.0/self.nsym/2 and phi<360.0/self.nsym):
							phi = 360.0/self.nsym-phi
							psi = 360.0 - psi
						elif(phi>=360.0/self.nsym+360.0/self.nsym/2 and phi<720.0/self.nsym):
							phi = 720.0/self.nsym-phi+360.0/self.nsym
							psi = 360.0 - psi
				else:
					if(theta>90.0):
						mat = rotmatrix(phi,theta,psi)
						for l in range(self.nsym):
							p1,p2,p3 = recmat( mulmat( mat , self.symatrix[l]) )
							#print(p1,p2,p3)
							if(self.is_in_subunit(p1, p2, 1)):
								phi=p1; theta=p2; psi=p3
								#print("  FOUND ")
								break
						

		return phi, theta, psi

	"""
	def reduce_normal(self, phi, theta, psi, inc_mirror):
		from math import degrees, radians, sin, cos, tan, atan, acos, sqrt
		return False
	"""

	def even_angles(self, delta = 15.0, theta1=-1.0, theta2=-1.0, phi1=-1.0, phi2=-1.0, \
					method = 'S', phiEqpsi = "Zero", inc_mirror = 1):
		"""Create a list of Euler angles suitable for projections.
		   method is either 'S' - for Saff algorithm
					   or   'P' - for Penczek '94 algorithm
						 'S' assumes phi1<phi2 and phi2-phi1>> delta ;
		   symmetry  - if this is set to point-group symmetry (cn or dn) or helical symmetry with point-group symmetry (scn or sdn), \
						 it will yield angles from the asymmetric unit, not the specified range;
		"""

		from math      import pi, sqrt, cos, acos, tan, sin, radians, degrees
		from utilities import even_angles_cd
		angles = []
		phi2_org = phi2
		if(phi2_org < 0.0):  phi2_org = self.brackets[1][0] - 1.0e-7 # exclude right border of unit
		theta2_org = theta2
		if(theta2_org < 0.0): theta2_org = self.brackets[1][3]
		if(phi2<phi1 or theta2<theta1 or delta <= 0.0):  ERROR("even_angles","incorrect parameters (phi1,phi2,theta1,theta2,delta): %f   %f   %f   %f   %f"%(phi1,phi2,theta1,theta2,delta),1)
		if(phi1 < 0.0):  phi1 = 0.0
		if(phi2 < 0.0):  phi2 = self.brackets[inc_mirror][0] - 1.0e-7 # exclude right border of unit
		if(theta1 < 0.0): theta1 = 0.0
		if(theta2 < 0.0): theta2 = self.brackets[inc_mirror][3]
		#print " parameters (phi1,phi2,theta,theta2,delta): %f   %f   %f   %f   %f"%(phi1,phi2,theta1,theta2,delta)
		#
		if(self.sym[0] != "s"):
			"""Create a list of Euler angles suitable for projections.
			   method is either 'S' - for Saff algorithm
							  or   'P' - for Penczek '94 algorithm
					  'S' assumes phi1<phi2 and phi2-phi1>> delta ;
			   phiEqpsi  - set this to 'Minus', if you want psi=-phi;
			"""
			angles = []
			is_platonic_sym = self.sym[0] == "o" or self.sym[0] == "t" or self.sym[0] == "i"
			if (method == 'P'):
				theta = theta1
				while(theta <= theta2):
					phi = phi1
					if(theta==0.0 or theta==180.0): detphi = 2*phi2
					else:  detphi = delta/sin(radians(theta))
					while(phi<phi2):
						if(self.is_in_subunit(phi, theta, inc_mirror)): 	angles.append([phi, theta, 0.0])
						else:  	angles.append([phi, theta, 0.0])
						phi += detphi
					theta += delta
			else:
				# I have to use original phi2 and theta2 to compute Deltaz and wedgeFactor as otherwise
				# points for include mirror differ from do not include mirror.
				Deltaz  = cos(radians(theta2_org))-cos(radians(theta1))
				s       = delta*pi/180.0
				NFactor = 3.6/s
				wedgeFactor = abs(Deltaz*(phi2_org-phi1)/720.0)
				NumPoints   = int(NFactor*NFactor*wedgeFactor)
				angles.append([phi1, theta1, 0.0])
				# initialize loop
				phistep = phi2_org-phi1
				z1 = cos(radians(theta1))
				phi = phi1
				for k in range(1, NumPoints-1):
					z = z1 + Deltaz*k/(NumPoints-1)
					r = sqrt(1.0-z*z)
					phi = phi1+(phi + delta/r - phi1)%phistep
					theta = degrees(acos(z))
					if(theta>180.0):  break
					if(not self.is_in_subunit(phi, theta, inc_mirror)): continue
					angles.append([phi, theta, 0.0])
				#angles.append([p2,t2,0])  # This is incorrect, as the last angle is really the border, not the element we need. PAP 01/15/07
			if (phiEqpsi == 'Minus'):
				for k in range(len(angles)): angles[k][2] = (720.0 - angles[k][0])%360.0
			if( (self.sym[0] == "c" or self.sym[0] == "d") and ((theta2 == 180.) or (theta2 >= 180. and delta == 180.0))):  angles.append( [0.0, 180.0, 0.0] )

		"""
		#  Helical symmetry should not be here
		elif(self.sym[0]  == "s"):

			#if symetry is "s", deltphi=delta, theata intial=theta1, theta end=90, delttheta=theta2
			# for helical, theta1 cannot be 0.0
			if theta1 > 90.0:
				ERROR('theta1 must be less than 90.0 for helical symmetry', 'even_angles', 1)
			if theta1 == 0.0: theta1 =90.0
			theta_number = int((90.0 - theta1)/theta2)
			#for helical, symmetry = s or scn
			cn = int(self.sym[2:])
			for j in range(theta_number,-1, -1):

				if( j == 0):
					if (self.sym[1] =="c"):
						if cn%2 == 0:
							k=int(359.99/cn/delta)
						else:
							k=int(359.99/2/cn/delta)
					elif (self.sym[1] =="d"):
						if cn%2 == 0:
							k=int(359.99/2/cn/delta)
						else:
							k=int(359.99/4/cn/delta)
					else:
						ERROR("For helical strucutre, we only support scn and sdn symmetry","even_angles",1)

				else:
					if (self.sym[1] =="c"):
						k=int(359.99/cn/delta)
					elif (self.sym[1] =="d"):
						k=int(359.99/2/cn/delta)

				for i in range(k+1):
						angles.append([i*delta,90.0-j*theta2,90.0])
			"""
		return angles
