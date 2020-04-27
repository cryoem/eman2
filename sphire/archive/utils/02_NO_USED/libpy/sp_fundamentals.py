







































from __future__ import print_function
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
	pass#IMPORTIMPORTIMPORT from EMAN2 import autocorrelation, fp_flag
	return EMAN2_cppwrap.autocorrelation(e, EMAN2_cppwrap.fp_flag.CIRCULANT, center)

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
	pass#IMPORTIMPORTIMPORT from EMAN2 import autocorrelation, fp_flag
	return EMAN2_cppwrap.autocorrelation(e, EMAN2_cppwrap.fp_flag.CIRCULANT_NORMALIZED, center)

def acfp(e, center=True):
	pass#IMPORTIMPORTIMPORT from EMAN2 import autocorrelation, fp_flag
	return EMAN2_cppwrap.autocorrelation(e, EMAN2_cppwrap.fp_flag.PADDED, center)

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
	pass#IMPORTIMPORTIMPORT from EMAN2 import autocorrelation, fp_flag
	return EMAN2_cppwrap.autocorrelation(e, EMAN2_cppwrap.fp_flag.PADDED_NORMALIZED, center)

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
	pass#IMPORTIMPORTIMPORT from EMAN2 import autocorrelation, fp_flag
	return EMAN2_cppwrap.autocorrelation(e, EMAN2_cppwrap.fp_flag.PADDED_LAG, center)

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
	pass#IMPORTIMPORTIMPORT from EMAN2 import autocorrelation, fp_flag
	return EMAN2_cppwrap.autocorrelation(e, EMAN2_cppwrap.fp_flag.PADDED_NORMALIZED_LAG, center)

def __buildweights(m, kb):
	weights = EMAN2_cppwrap.EMData()
	weights.set_size(m,m,1)
	for iy in range(m):
		wy = kb.sinhwin(iy-m//2)
		for ix in range(m):
			wx = kb.sinhwin(ix-m//2)
			weights.set_value_at(ix,iy,wx*wy)
	return weights

# shortcuts to Fourier product functions
# Correlation functions









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
	pass#IMPORTIMPORTIMPORT from EMAN2 import correlation, fp_flag
	return EMAN2_cppwrap.correlation(e,f,EMAN2_cppwrap.fp_flag.CIRCULANT_NORMALIZED, center)

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
	pass#IMPORTIMPORTIMPORT from EMAN2 import correlation, fp_flag
	return EMAN2_cppwrap.correlation(e,f,EMAN2_cppwrap.fp_flag.PADDED, center)

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
	pass#IMPORTIMPORTIMPORT from EMAN2 import correlation, fp_flag
	return EMAN2_cppwrap.correlation(e,f,EMAN2_cppwrap.fp_flag.PADDED_NORMALIZED, center)

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
	pass#IMPORTIMPORTIMPORT from EMAN2 import correlation, fp_flag
	return EMAN2_cppwrap.correlation(e,f, EMAN2_cppwrap.fp_flag.PADDED_LAG, center)

def ccfnpl(e, f, center=True):
	"""
		Name
			ccfnpl - calculate the normalized cross-correlation function between two images
		Input
			e: input image (real)
			ref: second input image (real) 
			center: if set to True (default), the origin of the result is at the center
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import correlation, fp_flag
	return EMAN2_cppwrap.correlation(e,f,EMAN2_cppwrap.fp_flag.PADDED_NORMALIZED_LAG, center)
    
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
	pass#IMPORTIMPORTIMPORT from EMAN2 import convolution, fp_flag
	return EMAN2_cppwrap.convolution(e,f,EMAN2_cppwrap.fp_flag.CIRCULANT, center)

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
	pass#IMPORTIMPORTIMPORT from EMAN2 import convolution, fp_flag
	return EMAN2_cppwrap.convolution(e,f,EMAN2_cppwrap.fp_flag.CIRCULANT_NORMALIZED, center)

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
	pass#IMPORTIMPORTIMPORT from EMAN2 import convolution, fp_flag
	return EMAN2_cppwrap.convolution(e,f,EMAN2_cppwrap.fp_flag.PADDED, center)

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
	pass#IMPORTIMPORTIMPORT from EMAN2 import convolution, fp_flag
	return EMAN2_cppwrap.convolution(e,f,EMAN2_cppwrap.fp_flag.PADDED_NORMALIZED, center)

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
	pass#IMPORTIMPORTIMPORT from EMAN2 import convolution, fp_flag
	return EMAN2_cppwrap.convolution(e,f,EMAN2_cppwrap.fp_flag.PADDED_LAG, center)

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
	pass#IMPORTIMPORTIMPORT from EMAN2 import convolution, fp_flag
	return EMAN2_cppwrap.convolution(e,f,EMAN2_cppwrap.fp_flag.PADDED_NORMALIZED_LAG, center)
    
    
# Selfcorrelation functions













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
	pass#IMPORTIMPORTIMPORT from EMAN2 import self_correlation, fp_flag
	return EMAN2_cppwrap.self_correlation(e, EMAN2_cppwrap.fp_flag.CIRCULANT_NORMALIZED, center)

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
	pass#IMPORTIMPORTIMPORT from EMAN2 import self_correlation, fp_flag
	return EMAN2_cppwrap.self_correlation(e, EMAN2_cppwrap.fp_flag.PADDED, center)

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
	pass#IMPORTIMPORTIMPORT from EMAN2 import self_correlation, fp_flag
	return EMAN2_cppwrap.self_correlation(e, EMAN2_cppwrap.fp_flag.PADDED_NORMALIZED, center)

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
	pass#IMPORTIMPORTIMPORT from EMAN2 import self_correlation, fp_flag
	return EMAN2_cppwrap.self_correlation(e, EMAN2_cppwrap.fp_flag.PADDED_LAG, center)

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
	pass#IMPORTIMPORTIMPORT from EMAN2 import self_correlation, fp_flag
	return EMAN2_cppwrap.self_correlation(e, EMAN2_cppwrap.fp_flag.PADDED_NORMALIZED_LAG, center)
 
























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

























































































def image_decimate(img, decimation=2, fit_to_fft = True, frequency_low=0, frequency_high=0):
	"""
		Window 2D image to FFT-friendly size, apply Butterworth low pass filter,
		and decimate image by integer factor
	"""
	pass#IMPORTIMPORTIMPORT from sp_filter       import filt_btwl
	pass#IMPORTIMPORTIMPORT from sp_fundamentals import smallprime
	pass#IMPORTIMPORTIMPORT from sp_utilities    import get_image
	if type(img)     == str:	img=sp_utilities.get_image(img)
	nz       = img.get_zsize()
	if( nz > 1):                    sp_global_def.ERROR("This command works only for 2-D images", "image_decimate", 1)
	if decimation    <= 1  : 	sp_global_def.ERROR("Improper decimation ratio", "image_decimate", 1)
	if(decimation    == 1.0): 	return  img.copy()
	if frequency_low <= 0  :	
		frequency_low     = 0.5/decimation-0.02
		if frequency_low <= 0 : sp_global_def.ERROR("Butterworth pass-band frequency is too low","image_decimation",1)			
		frequency_high    = min(0.5/decimation + 0.02, 0.499)
	if fit_to_fft:
		nx       = img.get_xsize()
		ny       = img.get_ysize()
		nx_fft_m = smallprime(nx)
		ny_fft_m = smallprime(ny)
		e        = EMAN2_cppwrap.Util.window(img, nx_fft_m, ny_fft_m, 1, 0,0,0)
		e        = sp_filter.filt_btwl(e, frequency_low, frequency_high)
	else:
		e        = sp_filter.filt_btwl(img, frequency_low, frequency_high)
	return EMAN2_cppwrap.Util.decimate(e, int(decimation), int(decimation), 1)
















































































def fdownsample(img, sub_rate=0.5, RetReal = True):
	"""
		resample image based on the value of sub_rate.
		the input image can be either 2D image or 3D volume.
		sub_rate < 1.0, subsampling the image.
		sub_rate > 1.0, upsampling the image using new gridding interpolation.
		fit_to_fft will change the ouput image size to an fft_friendly size
	"""

	pass#IMPORTIMPORTIMPORT from sp_fundamentals import fdecimate
	pass#IMPORTIMPORTIMPORT from sp_utilities    import get_pixel_size, set_pixel_size

	if type(img) == str:
		pass#IMPORTIMPORTIMPORT from sp_utilities    import get_image
		img = sp_utilities.get_image(img)
	nx = img.get_xsize()
	if img.is_complex():
		nx -= (2-nx%2)
	ny = img.get_ysize()
	nz = img.get_zsize()
	if( ny == 1):  sp_global_def.ERROR("Only 2D or 3D images allowed","resample",1)
	if sub_rate == 1.0: return  img.copy()
	elif sub_rate < 1.0:
		nnx = int(nx*sub_rate+0.5)
		nny = int(ny*sub_rate+0.5)
		nnz = int(nz*sub_rate+0.5)
		e = fdecimate(img, nnx, nny, nnz, RetReal = RetReal)
	else:  #  sub_rate>1
		sp_global_def.ERROR("fdownsample","upscaling not implemented",1)
		"""Multiline Comment0"""
		#MULTILINEMULTILINEMULTILINE 0
		#MULTILINEMULTILINEMULTILINE 0
		#MULTILINEMULTILINEMULTILINE 0
			#MULTILINEMULTILINEMULTILINE 0
		#MULTILINEMULTILINEMULTILINE 0
			#MULTILINEMULTILINEMULTILINE 0
		#MULTILINEMULTILINEMULTILINE 0
			#MULTILINEMULTILINEMULTILINE 0
			#MULTILINEMULTILINEMULTILINE 0
			#MULTILINEMULTILINEMULTILINE 0
			#MULTILINEMULTILINEMULTILINE 0

		#MULTILINEMULTILINEMULTILINE 0
			#MULTILINEMULTILINEMULTILINE 0
			#MULTILINEMULTILINEMULTILINE 0
			#MULTILINEMULTILINEMULTILINE 0
			#MULTILINEMULTILINEMULTILINE 0
		#MULTILINEMULTILINEMULTILINE 0
			#MULTILINEMULTILINEMULTILINE 0
				#MULTILINEMULTILINEMULTILINE 0
				#MULTILINEMULTILINEMULTILINE 0
			#MULTILINEMULTILINEMULTILINE 0
				#MULTILINEMULTILINEMULTILINE 0
				#MULTILINEMULTILINEMULTILINE 0
		#MULTILINEMULTILINEMULTILINE 0

	# Automatically adjust pixel size for ctf parameters
	pass#IMPORTIMPORTIMPORT from sp_utilities import get_pixel_size, set_pixel_size
	apix = sp_utilities.get_pixel_size(e)
	apix /= sub_rate
	sp_utilities.set_pixel_size(e, apix)
	cc = e.get_attr_default("xform.projection", None)
	if cc:
		cp = cc.get_params("spider")
		cp["tx"] *= sub_rate
		cp["ty"] *= sub_rate
		pass#IMPORTIMPORTIMPORT from sp_utilities import set_params_proj
		sp_utilities.set_params_proj(e, [cp["phi"], cp["theta"], cp["psi"], -cp["tx"], -cp["ty"]]) # have to invert as set inverts them again
	cc = e.get_attr_default("xform.align2d", None)
	if cc:
		cp = cc.get_params("2D")
		cp["tx"] *= sub_rate
		cp["ty"] *= sub_rate
		pass#IMPORTIMPORTIMPORT from sp_utilities import set_params2D
		sp_utilities.set_params2D(e, [cp["alpha"], cp["tx"], cp["ty"], cp["mirror"], cp["scale"]])

	return 	e
























































def prep_refim_gridding(refim, wr, numr, mode = "F"):
	pass#IMPORTIMPORTIMPORT from sp_fundamentals import prepi
	nx = refim.get_xsize()
	ny = refim.get_ysize()
	cnx = nx//2+1
	cny = ny//2+1
	#precalculate rings
	temp,kb = prepi(refim)
	crefim = EMAN2_cppwrap.Util.Polar2Dmi(temp, cnx, cny, numr, mode, kb)
	EMAN2_cppwrap.Util.Normalize_ring(cimage, numr, 0 )
	EMAN2_cppwrap.Util.Frngs(crefim, numr)
	EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)
	return  crefim,kb










































































def rot_avg(e):
	"""Rotational average.
	   Returns a 1-D image containing a rotational average of image e.
	"""
	return e.rotavg()












def rot_avg_image(image_to_be_averaged):
	"""
	Rotational average
	Returns a 2-D or 3-D image containing a rotational average of image e
	"""
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT from sp_utilities import get_im
	if type(image_to_be_averaged) is bytes: image_to_be_averaged = sp_utilities.get_im(image_to_be_averaged)
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
	pass#IMPORTIMPORTIMPORT from sp_utilities import model_blank
	table = EMAN2_cppwrap.Util.rotavg_fourier(img)
	table = table[:len(table)//2]
	scale = (img.get_xsize() - 2*img.is_complex())*img.get_ysize()*img.get_zsize()
	scale = 4.0/scale/scale
	for i in range(len(table)): table[i] *= scale
	if lng:
		pass#IMPORTIMPORTIMPORT from math import log10
		for ir in range(1,len(table)): table[ir] = numpy.log10(table[ir])
		table[0] = table[1]
	ps = sp_utilities.model_blank(len(table))
	for i in range(len(table)): ps[i] = table[i]
	return ps

def rops_textfile(img, filename, lng = False):
	"""Rotational average of the periodogram stored as a text file.
	   Saves a text file (suitable for gnuplot) of the rotational average 
	   of the periodogram of img.
		Input image can be real or Fourier, can be rectangular
		output length mapped onto x-dimension length
	"""
	pass#IMPORTIMPORTIMPORT from sp_utilities import write_text_file
	table = EMAN2_cppwrap.Util.rotavg_fourier(img)
	table = table[:len(table)//2]
	scale = (img.get_xsize() - 2*img.is_complex())*img.get_ysize()*img.get_zsize()
	scale = 4.0/scale/scale
	for i in range(len(table)): table[i] *= scale
	if lng:
		pass#IMPORTIMPORTIMPORT from math import log10
		for ir in range(1,len(table)): table[ir] = numpy.log10(table[ir])
		table[0] = table[1]
	sp_utilities.write_text_file([list(range(len(table) ) ),table], filename)
	














































def rotshift2dg(image, ang, dx, dy, kb, scale = 1.0):
	"""Rotate and shift an image using gridding
	"""
	pass#IMPORTIMPORTIMPORT from math import radians
	pass#IMPORTIMPORTIMPORT from EMAN2 import Processor

	M = image.get_xsize()
	alpha = 1.75
	K = 6
	N = M*2  # npad*image size
	r = M/2
	v = K/2.0/N
	# first pad it with zeros in Fourier space
	o = image.FourInterpol(N,N,1,0)
	params = {"filter_type" : EMAN2_cppwrap.Processor.fourier_filter_types.KAISER_SINH_INVERSE,
	          "alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
	q = EMAN2_cppwrap.Processor.EMFourierFilter(o,params)
	o = fft(q)

	# gridding rotation
	return o.rot_scale_conv(numpy.radians(ang), dx, dy, kb, scale)







































































































def ft2polargrid(image, ring_length, nb, ne):
	"""
		resample to polar coordinates using gridding in Fourier space.
	"""
	pass#IMPORTIMPORTIMPORT from EMAN2 import Processor
	pass#IMPORTIMPORTIMPORT from sp_fundamentals import fftip, fft

	nx = image.get_xsize()
	# prepare 
	npad = 2
	N = nx*npad
	K = 6
	alpha = 1.75
	r = nx/2
	v = K/2.0/N
	kb = EMAN2_cppwrap.Util.KaiserBessel(alpha, K, r, v, N)

	image1 = image.copy()  # This step is needed, otherwise image will be changed outside the function
	# divide out gridding weights
	image1.divkbsinh(kb)
	# pad and center image, then FFT
	image1 = image1.norm_pad(False, npad)
	fftip(image1)
	# Put the origin of the (real-space) image at the center
	image1.center_origin_fft()
	return image1.ft2polargrid(ring_length, nb, ne, kb)

























































































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

	if scale == 0.0 :  sp_global_def.ERROR("scale=0 not allowed", "rot_shift3D_grid", 1)

	if mode == "cyclic":
		pass#IMPORTIMPORTIMPORT from math import radians
		if kb == None:
			o, kb = prepi3D(img)
		else:
			o = img
		# gridding rotation/shift:
		#if  mirror: o.process_inplace("xform.mirror", {"axis":'x'})
		return o.rot_scale_conv_new_3D(numpy.radians(phi), numpy.radians(theta), numpy.radians(psi), sx, sy, sz, kb, scale, wrap)
	elif mode == "background":
		pass#IMPORTIMPORTIMPORT from math import radians
		if kb == None:
			o, kb = prepi3D(img)
		else:
			o = img
		# gridding rotation
		#if  mirror: o.process_inplace("xform.mirror", {"axis":'x'})
		return o.rot_scale_conv_new_background_3D(numpy.radians(phi), numpy.radians(theta), numpy.radians(psi), sx, sy, sz, kb, scale, wrap)	
	else: sp_global_def.ERROR("rot_shift3D_grid mode not valid", "rot_shift3D_grid", 1)






















































def sinc2inv(nx):
	pass#IMPORTIMPORTIMPORT from math import sqrt
	s = sincinv(nx)
	return [i*i for i in s]

def sincinv(nx):
	pass#IMPORTIMPORTIMPORT from math import pi,sin
	cdf =numpy.pi/nx
	npad = 1
	nxb = nx/2/npad
	nxe = nxb + (nx/npad)%2
	s = [1.0]*nx
	for i in range( -nxb, nxe):
		if( i != 0 ):
			rrr=abs(i)
			s[i+nxb] = (rrr*cdf)/numpy.sin(rrr*cdf)
	return s

def welch_pw2(img, win_size=512, overlp_x=50, overlp_y=50, edge_x=0, edge_y=0):
	""" 
		Calculate the power spectrum using Welch periodograms (overlapped periodogram)
	"""
	pass#IMPORTIMPORTIMPORT from sp_fundamentals import window2d, ramp
	pass#IMPORTIMPORTIMPORT from EMAN2 import periodogram
	nx = img.get_xsize()
	ny = img.get_ysize()
	nx_fft = smallprime(nx)
	ny_fft = smallprime(ny)
	x_gaussian_hi = 1./win_size
	pass#IMPORTIMPORTIMPORT from sp_filter    import filt_gaussh
	e_fil = sp_filter.filt_gaussh(window2d(img,nx_fft,ny_fft,"l"), x_gaussian_hi)
	x38 = 100/(100-overlp_x) # normalization of % of the overlap in x 
	x39 = 100/(100-overlp_y) # normalization of % of the overlap in y
	x26 = int(x38*((nx-2*edge_x)/win_size-1)+1)  # number of pieces horizontal dim.(X)
	x29 = int(x39*((ny-2*edge_y)/win_size-1)+1)  # number of pieces vertical dim.(Y)
	iz = 0	
	pw2 = EMAN2_cppwrap.EMData()
	for iy in range(1, x29+1):	
		x21 = (win_size/x39)*(iy-1) + edge_y  #  y-direction it should start from 0 if edge_y=0	      
		for ix in  range(1, x26+1):			 
			x22 = (win_size/x38)*(ix-1) + edge_x  # x-direction it should start from 0 if edge_x =0
			wi  = window2d(e_fil, win_size, win_size, "l", x22, x21)
			iz  = iz+1
			if (iz == 1): pw2  = EMAN2_cppwrap.periodogram(ramp(wi))
			else:         pw2 += EMAN2_cppwrap.periodogram(ramp(wi))
	return  pw2/float(iz)

def welch_pw2_tilt_band(img,theta,num_bnd=-1,overlp_y=50,edge_x=0,edge_y=0,win_s=256):
	""" 
		1. Calculate the power spectra of tilt bands
		2. The tilt micrograph is rotated such that the tilt axis is vertical (along Y axis)
		3. edge_x and edge_y are removed from the micrograph
	""" 
	pass#IMPORTIMPORTIMPORT from EMAN2 import periodogram
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
	pass#IMPORTIMPORTIMPORT from sp_filter import filt_gaussh
	pass#IMPORTIMPORTIMPORT from sp_utilities import drop_image, rot_image
	# The input img is rotated such that tilt axis is vertical
	img2  = rot_image(img1,theta, 0, 0, 1.0,1.0)	
	e_fil = sp_filter.filt_gaussh(img2, x_gaussian_hi)
	del img1
	del img2
	x39 = 100/(100-overlp_y) # normalization of % of the overlap in y
	x29 = int(x39*((ny)/win_y-1)+1)  # number of pieces vertical dim.(Y)
	pw2 = EMAN2_cppwrap.EMData()
	pw2_band = []
	for ix in  range(1, num_bnd+1):
		x22 = (win_x)*(ix-1)# x-direction it should start from 0 if edge_x =0
		iz=0
		for iy in range(1, x29+1):	
			x21 = (win_y/x39)*(iy-1) #  y-direction it should start from 0 if edge_y=0	      			 
			wi = window2d(e_fil,win_x, win_y,"l",x22, x21)
			iz = iz+1
			if (iz == 1): pw2  = EMAN2_cppwrap.periodogram(ramp(wi))
			else:         pw2 += EMAN2_cppwrap.periodogram(ramp(wi))
		pw2/=float(iz)
		# drop_image(pw2,"band%03d"%(ix))
		pw2_band.append(pw2)	
	return 	pw2_band

























































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
	sp_global_def.sxprint("Bracket did not find a mimimum")        
 
















































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































