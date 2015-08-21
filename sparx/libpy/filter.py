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

def filt_median(f, nx, ny, nz = 1, kernelshape = "BLOCK"):
	"""
	Name
		filt_median - Calculate the median filtered image.
	Input
		input image
		nxk, nyk, nzk
			Size of the kernel on each dimension.
			- They have to be odd so that the center of kernal can be well defined.
			- Note: For CIRCULAR kernel, we currently only allow circle or sphere kernel, i.e., nxk and nyk must be same for 2 dimensional image and nxk, nyk and nzk must be same for 3 dimensional image.
		kernelshape
			Shape of the kernel
			- BLOCK is for block kernel (DEFAULT);
			- CIRCULAR is for circular kernel;
			- CROSS is for cross kernal;
			- Note: For 1-dimensional image, these three kernels degenerate into the same shape.
	Output
		median filtered image
	"""
	from EMAN2 import kernel_shape, filt_median_

	if kernelshape=="BLOCK":         return filt_median_(f,nx,ny,nz,kernel_shape.BLOCK)
	elif kernelshape=="CIRCULAR":    return filt_median_(f,nx,ny,nz,kernel_shape.CIRCULAR)
	elif kernelshape=="CROSS":       return filt_median_(f,nx,ny,nz,kernel_shape.CROSS)
	else: print "Unknown kernel shape."

# Fourier filters
def filt_tophatl(e, freq, pad = False):
	"""
		Name
			filt_tophatl - top-hat low-pass Fourier filter (truncation of a Fourier series)
		Input
			e: input image (can be either real or Fourier)
			freq: stop-band frequency
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			All frequencies are in absolute frequency units fa and their valid range is [0:0.5].
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
	from EMAN2 import Processor
	params = {"filter_type" : Processor.fourier_filter_types.TOP_HAT_LOW_PASS,
		"cutoff_abs" : freq, "dopad" : pad}
	return Processor.EMFourierFilter(e, params)
    
def filt_tophath(e, freq, pad = False):
	"""
		Name
			filt_tophath - top-hat high-pass Fourier filter (truncation of a Fourier series)
		Input
			e: input image (can be either real or Fourier)
			freq: pass-band frequency
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			All frequencies are in absolute frequency units fa and their valid range is [0:0.5].
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
	from EMAN2 import Processor
	params = {"filter_type" : Processor.fourier_filter_types.TOP_HAT_HIGH_PASS,
		  "cutoff_abs" : freq, "dopad" : pad}
	return Processor.EMFourierFilter(e, params)
    
def filt_tophatb(e, freql, freqh, pad = False):
	"""
		Name
			filt_tophatb - top-hat band-pass Fourier filter (truncation of a Fourier series)
		Input
			e: input image (can be either real or Fourier)
			freql: low-end frequency of the filter pass-band
			freqh: high-end frequency of the filter pass-band
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			All frequencies are in absolute frequency units fa and their valid range is [0:0.5].
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
	from EMAN2 import Processor
	params = {"filter_type" : Processor.fourier_filter_types.TOP_HAT_BAND_PASS,
		  "low_cutoff_frequency" : freql, "high_cutoff_frequency" : freqh, "dopad" : pad}
	return Processor.EMFourierFilter(e, params)
    
def filt_tophato(e, freql, freqh, value, pad = False):
	"""
		Name
			filt_tophato - top-hat homomorphic Fourier filter (truncation of a Fourier series)
		Input
			e: input image (can be either real or Fourier)
			freql: low-end frequency of the filter pass-band
			freqh: high-end frequency of the filter pass-band
			value: value of the filter in spatial frequencies lower than freql
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			All frequencies are in absolute frequency units fa and their valid range is [0:0.5].
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier

	"""
	from EMAN2 import Processor
	params = {"filter_type" : Processor.fourier_filter_types.TOP_HOMOMORPHIC,
		  "low_cutoff_frequency" : freql, "high_cutoff_frequency" : freqh, "value_at_zero_frequency" : value, "dopad" : pad}
	return Processor.EMFourierFilter(e, params)

    

def filt_gaussl(e, sigma, pad = False):
	"""
		Name
			filt_gaussl - Gaussian low-pass Fourier filter
		Input
			e: input image (can be either real or Fourier)
			sigma: standard deviation of the Gaussian function in absolute frequency units fa.
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			Note: for sigma = 0.5 (the Nyquist frequency) the value of the filter at the maximum frequency is G(fN)=1e=0.61.
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
	from EMAN2 import Processor
	params = {"filter_type" : Processor.fourier_filter_types.GAUSS_LOW_PASS,
		  "cutoff_abs" : sigma, "dopad" : pad}
	return Processor.EMFourierFilter(e, params)
    
def filt_gaussinv(e, sigma, pad = False):
	"""
		Name
			filt_gaussinv - inverse Gaussian (high-pass) Fourier filter (division by a Gaussian function in Fourier space)
		Input
			e: input image (can be either real or Fourier)
			sigma: standard deviation of the Gaussian function in absolute frequency units fa.
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			Note: for sigma = 0.5 (the Nyquist frequency) the value of the filter at the maximum frequency is G(fN)=e=1.65.
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
	from EMAN2 import Processor
	params = {"filter_type" : Processor.fourier_filter_types.GAUSS_INVERSE,
		  "cutoff_abs" : sigma, "dopad" : pad}
	return Processor.EMFourierFilter(e, params)
 
def filt_gaussh(e, sigma, pad = False):
	"""
		Name
			filt_gaussh - Gaussian high-pass Fourier filter
		Input
			e: input image (can be either real or Fourier)
			sigma: standard deviation of the Gaussian function in absolute frequency units fa.
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			Note: for sigma = 0.5 (the Nyquist frequency) the value of the filter at the maximum frequency is 1-G(0)=1-1e=0.39
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
	from EMAN2 import Processor
	params = {"filter_type" : Processor.fourier_filter_types.GAUSS_HIGH_PASS,
		  "cutoff_abs" : sigma, "dopad" : pad}
	return Processor.EMFourierFilter(e, params)
    
def filt_gaussb(e, sigma, center, pad = False):
	"""
		Name
			filt_gaussb - Gaussian band-pass Fourier filter
		Input
			e: input image (can be either real or Fourier)
			sigma: standard deviation of the Gaussian function in absolute frequency units fa.
			center: frequency in absolute frequency units fa at which the filter reaches maximum (=1).
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
	from EMAN2 import Processor
	params = {"filter_type" : Processor.fourier_filter_types.GAUSS_BAND_PASS,
		  "cutoff_abs" : sigma, "center" : center, "dopad" : pad}
	return Processor.EMFourierFilter(e, params)
    
def filt_gausso(e, sigma, value, pad = False):
	"""
	Name
		filt_gausso - Gaussian homomorphic Fourier filter
	Input
		e: input image (can be either real or Fourier)
		sigma: standard deviation of the Gaussian function in absolute frequency units fa.
		value: center - value of the filter at zero freqeuncy
		pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
	Output
		filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
	from EMAN2 import Processor
	params = {"filter_type" : Processor.fourier_filter_types.GAUSS_HOMOMORPHIC,
		  "cutoff_abs" : sigma, "value_at_zero_frequency" : value, "dopad" : pad}
	return Processor.EMFourierFilter(e, params)



def filt_btwl(e, freql, freqh, pad = False):
	"""
		Name
			filt_btwl - Butterworth low-pass Fourier filter
		Input
			e: input image (can be either real or Fourier)
			freql: low - pass-band frequency
			freqh: high - stop-band frequency
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			All frequencies are in absolute frequency units fa and their valid range is [0:0.5].
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
	from EMAN2 import Processor
	params = {"filter_type" : Processor.fourier_filter_types.BUTTERWORTH_LOW_PASS,
		  "low_cutoff_frequency" : freql, "high_cutoff_frequency" : freqh, "dopad" : pad}
	return Processor.EMFourierFilter(e, params)
    
def filt_btwh(e, freql, freqh, pad = False):
	"""
		Name
			filt_btwh - Butterworth high-pass Fourier filter

		Input
			e: input image (can be either real or Fourier)
			freql: low - stop-band frequency
			freqh: high - pass-band frequency
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			All frequencies are in absolute frequency units fa and their valid range is [0:0.5].
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
	from EMAN2 import Processor
	params = {"filter_type" : Processor.fourier_filter_types.BUTTERWORTH_HIGH_PASS,
		  "low_cutoff_frequency" : freql, "high_cutoff_frequency" : freqh, "dopad" : pad}
	return Processor.EMFourierFilter(e, params)
    
def filt_btwo(e, freql, freqh, value, pad = False):
	"""
		Name
			filt_btwo - Butterworth homomorphic Fourier filter
		Input
			e: input image (can be either real or Fourier)
			freql: low - stop-band frequency
			freqh: high - pass-band frequency
			value: value of the filter at zero frequency
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			All frequencies are in absolute frequency units fa and their valid range is [0:0.5].
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier

	"""
	from EMAN2 import Processor
	params = {"filter_type" : Processor.fourier_filter_types.BUTTERWORTH_HOMOMORPHIC,
		  "low_cutoff_frequency" : freql, "high_cutoff_frequency" : freqh,
	    	  "value_at_zero_frequency" : value, "dopad" : pad}
	return Processor.EMFourierFilter(e, params)
   


def filt_tanl(e, freq, fall_off, pad = False):
	"""
		Name
			filt_tanl - hyperbolic tangent low-pass Fourier filter
		Input
			e: input image (can be either real or Fourier)
			freq: stop-band frequency
			fall_off: fall off of the filter
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			All frequencies are in absolute frequency units fa and their valid range is [0:0.5].
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
	from EMAN2 import Processor
	params = {"filter_type" : Processor.fourier_filter_types.TANH_LOW_PASS,
		  "cutoff_abs" : freq, "fall_off": fall_off, "dopad" : pad}
	return Processor.EMFourierFilter(e, params)
    
def filt_tanh(e, freq, fall_off, pad = False):
	"""
		Name
			filt_tanh - hyperbolic tangent high-pass Fourier filter
		Input
			e: input image (can be either real or Fourier)
			freq: pass-band frequency
			fall_off: fall off of the filter
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			All frequencies are in absolute frequency units fa and their valid range is [0:0.5].
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
	from EMAN2 import Processor
	params = {"filter_type" : Processor.fourier_filter_types.TANH_HIGH_PASS,
		  "cutoff_abs" : freq, "fall_off": fall_off, "dopad" : pad}
	return Processor.EMFourierFilter(e, params)
    
def filt_tanb(e, freql, low_fall_off, freqh, high_fall_off, pad = False):
	"""
		Name
			filt_tanb - hyperbolic tangent band-pass Fourier filter
		Input
			e: input image (can be either real or Fourier)
			freql: low-end frequency of the filter pass-band
			low_fall_off: fall off of the filter at the low-end frequency
			freqh: high-end frequency of the filter pass-band
			high_fall_off: fall off of the filter at the high-and frequency
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			All frequencies are in absolute frequency units fa and their valid range is [0:0.5].
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
	from EMAN2 import Processor
	params = {"filter_type" : Processor.fourier_filter_types.TANH_BAND_PASS,
		  "low_cutoff_frequency" : freql, "Low_fall_off": low_fall_off,
	    	  "high_cutoff_frequency" : freqh, "high_fall_off": high_fall_off, "dopad" : pad}
	return Processor.EMFourierFilter(e, params)
    
def filt_tano(e, freq, fall_off, value, pad = False):
	"""
		Name
			filt_tano - hyperbolic tangent homomorphic Fourier filter
		Input
			e: input image (can be either real or Fourier)
			freq: pass-band frequency
			fall_off: fall off of the filter
			value: value of the filter at zero frequency
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			All frequencies are in absolute frequency units fa and their valid range is [0:0.5].
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier.
	"""
	from EMAN2 import Processor
	params = {"filter_type" : Processor.fourier_filter_types.TANH_HOMOMORPHIC,
		  "cutoff_abs" : freq, "fall_off": fall_off,
	    	  "value_at_zero_frequency" : value, "dopad" : pad}
	return Processor.EMFourierFilter(e, params)


  
def filt_kaisersinh(e, alpha):
	from EMAN2 import Processor
	M = e.get_xsize()
	K = 6
	N = M*2  # npad*image size
	r=M/2
	v=K/2.0/N
	from EMAN2 import Processor
	params = {"filter_type":Processor.fourier_filter_types.KAISER_SINH,
		  "alpha":alpha, "K":K,"r":r,"v":v,"N":N}
	return Processor.EMFourierFilter(e, params)
    
def filt_kaisersinhp(e, alpha):
	from EMAN2 import Processor
	M = e.get_xsize()
	K = 6
	N = M*2  # npad*image size
	r=M/2
	v=K/2.0/N
	from EMAN2 import Processor
	params = {"filter_type":Processor.fourier_filter_types.KAISER_SINH,
		  "dopad" : 1, "alpha":alpha, "K":K,"r":r,"v":v,"N":N}
	return Processor.EMFourierFilter(e, params)
    
def filt_kaisersinhinv(e, alpha):
	from EMAN2 import Processor
	M = e.get_xsize()
	K = 6
	N = M*2  # npad*image size
	r=M/2
	v=K/2.0/N
	from EMAN2 import Processor
	params = {"filter_type":Processor.fourier_filter_types.KAISER_SINH_INVERSE,
		  "alpha":alpha, "K":K,"r":r,"v":v,"N":N}
	return Processor.EMFourierFilter(e, params)
    
def filt_kaisersinhinvp(e, alpha):
	from EMAN2 import Processor
	M = e.get_xsize()
	K = 6
	N = M*2  # npad*image size
	r=M/2
	v=K/2.0/N
	from EMAN2 import Processor
	params = {"filter_type":Processor.fourier_filter_types.KAISER_SINH_INVERSE,
		  "dopad" : 1, "alpha":alpha, "K":K,"r":r,"v":v,"N":N}
	return Processor.EMFourierFilter(e, params)

def filt_table(e, table):
	"""
		Name
			filt_table - filter image using a user-provided (as a list) of Fourier filter values
		Input
			e: input image (real or Fourier)
			table: user-provided 1-D table of filter values.
		Output
			image Fourier-filtered using coefficients in table (real or Fourier, according the input image)
		Options
			fast: use the fast method; may combust certain computers.
			huge: gobble memory; there is plenty of it, anyway.
	"""
	from EMAN2 import Processor
	params = {"filter_type" : Processor.fourier_filter_types.RADIAL_TABLE,
			"table" : table}
	return Processor.EMFourierFilter(e, params)

def filt_ctf(img, ctf, dopad=True, sign=1, binary = 0):
	"""
		Name
			filt_ctf - apply Contrast Transfer Function (CTF) to an image in Fourier space
		Input
			image: input image, it can be real or complex
			ctf: an CTF object, please see CTF_info for description.
			pad: apply padding with zeroes in real space to twice the size before CTF application (Default is True, if set to False, no padding, if input image is Fourier, the flag has no effect).
			sign: sign of the CTF. If cryo data had inverted contrast, i.e., images were multiplied by -1 and particle projections appear bright on dark background, it has to be set to -1). (Default is 1).
			binary: phase flipping if set to 1 (default is 0).
		Output
			image multiplied in Fourier space by the CTF, the output image has the same format as the input image.
	"""
	from EMAN2 import Processor
	assert img.get_ysize() > 1
	dict        = ctf.to_dict()
	dz          = dict["defocus"]
	cs          = dict["cs"]
	voltage     = dict["voltage"]
	pixel_size  = dict["apix"]
	b_factor    = dict["bfactor"]
	ampcont     = dict["ampcont"]
	dza         = dict["dfdiff"]
	azz         = dict["dfang"]

	if dopad and not img.is_complex():  ip = 1
	else:                               ip = 0

	params = {"filter_type": Processor.fourier_filter_types.CTF_,
		"defocus" : dz,
		"Cs": cs,
		"voltage": voltage,
		"Pixel_size": pixel_size,
		"B_factor": b_factor,
		"amp_contrast": ampcont,
		"dopad": ip,
		"binary": binary,
		"sign": sign,
		"dza": dza,
		"azz":azz}
	tmp = Processor.EMFourierFilter(img, params)
	tmp.set_attr_dict({"ctf":ctf})
	return tmp

def filt_unctf(e, dz, cs, voltage, pixel, wgh=0.1, b_factor=0.0, sign=-1.0, dza=0.0, azz=0.0):
	"""
		defocus: Angstrom
		cs: mm
		voltage: Kv
		pixel size: Angstrom
		b_factor: Angstrom^2
		wgh: Unitless
	"""	
	from EMAN2 import Processor
	params = {"filter_type": Processor.fourier_filter_types.CTF_,
		"defocus" : dz,
		"Cs" : cs,
		"voltage" : voltage,
		"Pixel_size" : pixel,
		"B_factor" : b_factor,
		"amp_contrast": wgh,
		"sign": sign,
		"dza": dza,
		"azz":azz,
		"undo": 1}
	return Processor.EMFourierFilter(e, params)

def filt_params(dres, high = 0.95, low = 0.1):
	"""
		Name
			filt_params - using the FSC curve establish two frequencies: stop-band, corresponding to the point where the curves drops below predefined high value, and pass-band, corresponding to the point where the curves drops below predefined low value
		Input
			dres - list produced by the fsc funcion
			dres[0] - absolute frequencies
			dres[1] - fsc, because it was calculated from half the dataset, convert it to full using rn = 2r/(1+r)
			dres[2] - number of point use to calculate fsc coeff
			low: cutoff of the fsc curve
			high: cutoff of the fsc curve. Note: high and low values are optional.
			return parameters of the Butterworth filter (low and high frequencies)
		Output
			freql: stop-band frequency
			freqh: pass-band frequency
	"""
	n = len(dres[1])
	if ( (2*dres[1][n-1]/(1.0+dres[1][n-1])) > low):
		# the curve never falls below preset level, most likely something's wrong; however, return reasonable values
		#  First, find out whether it actually does fall
		nend = 0
		for i in xrange(n-1,1,-1):
			if ( (2*dres[1][i]/(1.0+dres[1][i]) ) < low):
				nend = i
				break
		if( nend == 0 ):
			#it does not fall anywhere
			freql = 0.4
			freqh = 0.4999
			return freql,freqh
		else:
			n = nend +1
	#  normal case
	freql = 0.001 # this is in case the next loop does not work and the curve never drops below high
	for i in xrange(1,n-1):
		if ( (2*dres[1][i]/(1.0+dres[1][i]) ) < high):
			freql = dres[0][i]
			break
	freqh = 0.1 # this is in case the next loop does not work and the curve never rises above low
	for i in xrange(n-2,0,-1):
		if ( (2*dres[1][i]/(1.0+dres[1][i]) ) > low):
			freqh =min( max(dres[0][i], freql+0.1), 0.45)
			break
	return freql,freqh

def filt_from_fsc(dres, low = 0.1):
	"""
		Name
			filt_from_fsc - generate a filter using the FSC curve
		Input
			dres: - list produced by the fsc funcion
			dres[0] - absolute frequencies
			dres[1] - fsc, because it was calculated from the dataset split into halves, convert it to full using rn = 2r/(1+r)
			dres[2] - number of point use to calculate fsc coeff
			low (low is option): cutoff of the fsc curve
			return filter curve rn = 2r/(1+r)
		Output
			filtc: a list of filter values
			Note: it is assumed that the FSC curve was computed from the dataset split into halves, so the filter values are calculated using rf = 2r1+r

	"""
	n = len(dres[1])
	filtc = [0.0]*n
	filtc[0] = 1.0
	last = n
	for i in xrange(1, n-1):
		if dres[1][i] < 0.0:
			qt = 0.0
		else:
			qt = 2*(dres[1][i-1]/(1.0+dres[1][i-1]) + dres[1][i]/(1.0+dres[1][i]) + dres[1][i+1]/(1.0+dres[1][i+1]))/3.0
		if qt < low:
			filtc[i] = low
			last = i
			break
		else:
			filtc[i] = qt
	if last < n-1:
		for i in xrange(last, min(last+4, n)):
			qt /= 1.2
			filtc[i] = qt
	return filtc

def filt_from_fsc2(dres, low = 0.1):
	"""
		dres - list produced by the fsc funcion
		dres[0] - absolute frequencies
		dres[1] - fsc, because it was calculated from the dataset split into halves, convert it to full using rn = 2r/(1+r)
		dres[2] - number of point use to calculate fsc coeff
		low cutoff of the fsc curve
		return filter curve rn = 2r/(1+r)
	"""
	n = len(dres[1])
	filtc = [0.0]*n
	filtc[0] = 1.0
	last = n
	for i in xrange(1,n-1):
		if(dres[1][i]<0.0):
			qt = 0.0
		else:
			qt = (dres[1][i-1]/(1.0+dres[1][i-1]) + dres[1][i]/(1.0+dres[1][i]) + dres[1][i+1]/(1.0+dres[1][i+1]))/3.0
		if ( qt < low ):
			filtc[i] = low
			last = i
			break
		else:
			filtc[i] = qt
	if(last<n-1):
		for i in xrange(last, min(last+4,n)):
			qt /= 1.2
			filtc[i] = qt
	return filtc

def filt_from_fsc_bwt(dres, low = 0.1):
	"""
		Name
			filt_from_fsc_bwt - Create filter using Fourier Ring correlation coefficients with Butterworth filter (bwf) cut off.
		Input
			dres - list produced by the fsc funcion
			dres[i][0] - absolute frequencies
			dres[i][1] - fsc, because it was calculated from the dataset split into halves, convert it to full using rn = 2r/(1+r)
			dres[i][2] - number of point use to calculate fsc coeff
			return filter curve rn = 2r/(1+r) when frc>0.9
		                    rn =1./sqrt(1.+(dres[i][0]/rad)**order) when frc<=0.9
			low: pass band
			ssnr calculated from frc with Butterworth fitler-like cutoff
		Output
			filtc: filter created from FRC
	"""
	from math import log,sqrt
	n = len(dres[0])
	filtc = [0.0]*n
	eps=0.882 # copied from spider fq np command 
	a=10.624 # copied from spider fq np command 
	lowf=.45
	for i in xrange(n):
		if(dres[1][i]<low):
			lowf=dres[0][i]
			break
	highf=lowf+.05
	order=2.*log(eps/sqrt(a**2-1))/log(lowf/highf)
	rad=lowf/eps**(2./order)
	for i in xrange(n):		
		if(dres[1][i]<low): 
			qt = 1./sqrt(1.+(dres[0][i]/rad)**order)
		else:
			qt = 2*dres[1][i]/(1.0+dres[1][i])
		filtc[i] = qt
	return  filtc


def fit_tanh(dres, low = 0.1):
	"""
		dres - list produced by the fsc funcion
		dres[0] - absolute frequencies
		dres[1] - fsc, because it was calculated from the dataset split into halves, convert it to full using rn = 2r/(1+r)
		dres[2] - number of points use to calculate fsc coeff
		low cutoff of the fsc curve
		return parameters of the tanh filter: freq - cutoff frequency at which filter value is 0.5, and fall_off, the 'width' of the filter
	"""
	def fit_tanh_func(args, data):
		from math import pi, tanh
		v = 0.0

		if(data[1][0] < 0.0 ):
			data[1][0] *= -1.0

		for i in xrange(len(data[0])):
			fsc =  2*data[1][i]/(1.0+data[1][i])
			if args[0]==0 or args[1]==0: qt=0
			else: qt  = fsc - 0.5*( tanh(pi*(data[0][i]+args[0])/2.0/args[1]/args[0]) - tanh(pi*(data[0][i]-args[0])/2.0/args[1]/args[0]) )
			v  -= qt*qt
		#print args,v
		return v

	setzero = False
	for i in xrange(1,len(dres[0])):
		if not setzero:
			if(2*dres[1][i]/(1.0+dres[1][i]) <low):  setzero = True
		if setzero:  dres[1][i] = 0.0

	freq = -1.0
	for i in xrange(1,len(dres[0])-1):
		if ( (2*dres[1][i]/(1.0+dres[1][i]) ) < 0.5):
			freq = dres[0][i-1]
			break
	if freq < 0.0:
		# the curve never falls below 0.5, most likely something's wrong; however, return reasonable values
		freq = 0.4
		fall_off = 0.2
		return freq, fall_off

	from utilities import amoeba
	args   = [freq, 0.1]
	scale  = [0.05, 0.05]
	result = amoeba(args, scale, fit_tanh_func, data = dres)

	'''
	args[0] = result[0][0]
	args[1] = result[0][1]
	#print  args
	from math import pi, tanh
	for i in xrange(len(dres[0])):
		fsc = 2*dres[1][i]/(1.0+dres[1][i])
		print i, dres[0][i],fsc , 0.5*( tanh(pi*(dres[0][i]+args[0])/2.0/args[1]/args[0]) - tanh(pi*(dres[0][i]-args[0])/2.0/args[1]/args[0]) )
		#qt  = fsc - 0.5*( tanh(pi*(dres[0][i]+args[0])/2.0/args[1]/args[0]) - tanh(pi*(dres[0][i]-args[0])/2.0/args[1]/args[0]) )
	'''
	return result[0][0], result[0][1]
	

def fit_tanh1(dres, low = 0.1):
	"""
		dres - list produced by the fsc funcion
		dres[0] - absolute frequencies
		dres[1] - fsc, to be conservative, do not use factor of 2.
		dres[2] - number of points use to calculate fsc coeff
		low cutoff of the fsc curve
		return parameters of the tanh filter: freq - cutoff frequency at which filter value is 0.5, and fall_off, the 'width' of the filter
	"""
	def fit_tanh_func(args, data):
		from math import pi, tanh
		v = 0.0
		for i in xrange(len(data[0])):
			fsc =  data[1][i]
			if args[0]==0 or args[1]==0: qt=0
			else: qt  = fsc - 0.5*( tanh(pi*(data[0][i]+args[0])/2.0/args[1]/args[0]) - tanh(pi*(data[0][i]-args[0])/2.0/args[1]/args[0]) )
			v  -= qt*qt
		#print args,v
		return v

	setzero = False
	for i in xrange(1,len(dres[0])):
		if not setzero:
			if(dres[1][i] <low):  setzero = True
		if setzero:  dres[1][i] = 0.0

	freq = -1.0
	for i in xrange(1,len(dres[0])-1):
		if ( dres[1][i] < 0.5):
			freq = dres[0][i-1]
			break
	if(freq < 0.0):
		# the curve never falls below 0.5, most likely something's wrong; however, return reasonable values
		freq = 0.2
		fall_off = 0.2
		return freq, fall_off
	from utilities import amoeba
	args   = [freq, 0.1]
	scale  = [0.05, 0.05]
	result = amoeba(args, scale, fit_tanh_func, data = dres)

	'''
	args[0] = result[0][0]
	args[1] = result[0][1]
	#print  args
	from math import pi, tanh
	for i in xrange(len(dres[0])):
		fsc = 2*dres[1][i]/(1.0+dres[1][i])
		print i, dres[0][i],fsc , 0.5*( tanh(pi*(dres[0][i]+args[0])/2.0/args[1]/args[0]) - tanh(pi*(dres[0][i]-args[0])/2.0/args[1]/args[0]) )
		#qt  = fsc - 0.5*( tanh(pi*(dres[0][i]+args[0])/2.0/args[1]/args[0]) - tanh(pi*(dres[0][i]-args[0])/2.0/args[1]/args[0]) )
	'''
	return result[0][0], result[0][1]

def tanhfilter(nx, fl, aa):
	#  generate discretized tanh filter
	from math import pi, tanh
	n = nx//2 + nx%2
	f = [0.0]*n
	for i in xrange(n):
		x = float(i)/nx
		f[i] = 0.5*( tanh(pi*(x+fl)/2.0/aa/fl) - tanh(pi*(x-fl)/2.0/aa/fl) )
	return  [[float(i)/nx for i in xrange(n)],f]

def filt_matched(ima, SNR, Pref):
	""" 
		Calculate Matched filter. The reference image is 
		assumed to be a Wiener average
		See paper Alignment under noise by PRB et al
	"""

	from filter import filt_from_fsc_bwt, filt_table
	from math import sqrt
	from EMAN2 import EMData
	from fundamentals import rops_table
	
	ctf_2 = ima.get_attr('ctf_2')
	PU = ima.get_attr('PU')
	Pn1 = ima.get_attr('Pn1')
	Pn2 = ima.get_attr('Pn2')
	TE= ima.get_attr('TE')
	HM=[]
	TMP1=[]
	TMP2=[]
	for j in xrange(len(Pref)-1):
		if(SNR[j]>.05): 
			thm=SNR[j]*(SNR[j]+1.)*PU[j]/Pref[j]
			print thm
			hm=sqrt(thm)
			deno=(SNR[j]+1)*(ctf_2[j]*Pn2[j]*TE[j]**2+Pn1[j])+ctf_2[j]*PU[j]*TE[j]**2
			xval=hm/deno 
		else: 
			hm=0.0
			xval=0.0
		HM.append(xval)
		TMP1.append(Pref[j])
		TMP2.append(hm)
		
	img=filt_table(ima,HM)
	res=[]
	res.append(img)
	res.append(TMP1)
	res.append(TMP2)
	del HM
	return res	

def filt_vols( vols, fscs, mask3D ):
	from math          import sqrt
	from filter        import fit_tanh, filt_tanl, filt_table
	from fundamentals  import rops_table
	from morphology    import threshold

	flmin = 1.0
	flmax = -1.0
	nvol = len(vols)
	for i in xrange(nvol):
		fl, aa = fit_tanh( fscs[i] )
		if (fl < flmin):
			flmin = fl
			aamin = aa
		if (fl > flmax):
			flmax = fl
			idmax = i
	print " Filter tanl, parameters: ",flmin-0.05, "  ",  aamin
	volmax = vols[idmax]
	volmax = filt_tanl( volmax, flmin-0.05, aamin )
	pmax = rops_table( volmax )

	for i in xrange(nvol):
		ptab = rops_table( vols[i] )
		for j in xrange( len(ptab) ):
			ptab[j] = sqrt( pmax[j]/ptab[j] )

		vols[i] = filt_table( vols[i], ptab )
		#stat = Util.infomask( vols[i], mask3D, False )
		#volf -= stat[0]
		Util.mul_img( vols[i], mask3D )
		#volf = threshold( volf )

	return vols

def filterlocal(ui, vi, m, falloff, myid, main_node, number_of_proc):
	from mpi 	  	  import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	  	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_gatherv, mpi_send, mpi_recv
	from mpi 	  	  import MPI_SUM, MPI_FLOAT, MPI_INT, MPI_TAG_UB
	from utilities import bcast_number_to_all, bcast_list_to_all, model_blank, bcast_EMData_to_all, reduce_EMData_to_root
	from morphology import threshold_outside
	from filter import filt_tanl
	from fundamentals import fft, fftip

	if(myid == main_node):

		nx = vi.get_xsize()
		ny = vi.get_ysize()
		nz = vi.get_zsize()
		#  Round all resolution numbers to two digits
		for x in xrange(nx):
			for y in xrange(ny):
				for z in xrange(nz):
					ui.set_value_at_fast( x,y,z, round(ui.get_value_at(x,y,z), 2) )
		dis = [nx,ny,nz]
	else:
		falloff = 0.0
		radius  = 0
		dis = [0,0,0]
	falloff = bcast_number_to_all(falloff, main_node)
	dis = bcast_list_to_all(dis, myid, source_node = main_node)

	if(myid != main_node):
		nx = int(dis[0])
		ny = int(dis[1])
		nz = int(dis[2])

		vi = model_blank(nx,ny,nz)
		ui = model_blank(nx,ny,nz)

	bcast_EMData_to_all(vi, myid, main_node)
	bcast_EMData_to_all(ui, myid, main_node)

	fftip(vi)  #  volume to be filtered

	st = Util.infomask(ui, m, True)


	filteredvol = model_blank(nx,ny,nz)
	cutoff = max(st[2] - 0.01,0.0)
	while(cutoff < st[3] ):
		cutoff = round(cutoff + 0.01, 2)
		#if(myid == main_node):  print  cutoff,st
		pt = Util.infomask( threshold_outside(ui, cutoff - 0.00501, cutoff + 0.005), m, True)  # Ideally, one would want to check only slices in question...
		if(pt[0] != 0.0):
			#print cutoff,pt[0]
			vovo = fft( filt_tanl(vi, cutoff, falloff) )
			for z in xrange(myid, nz, number_of_proc):
				for x in xrange(nx):
					for y in xrange(ny):
						if(m.get_value_at(x,y,z) > 0.5):
							if(round(ui.get_value_at(x,y,z),2) == cutoff):
								filteredvol.set_value_at_fast(x,y,z,vovo.get_value_at(x,y,z))

	mpi_barrier(MPI_COMM_WORLD)
	reduce_EMData_to_root(filteredvol, myid, main_node, MPI_COMM_WORLD)
	return filteredvol
