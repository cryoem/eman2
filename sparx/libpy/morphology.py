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
import numpy

from builtins import range
from global_def import *
import global_def
import scipy.ndimage.morphology as snm
import EMAN2_cppwrap


def fill_soft_edge_kernel_mask(kernel_mask, length, mode):
	"""
	Get the soft edge kernel value at the specified position.
	If the position is greater than the length, the value is zero.

	Arguments:
	position - Kernel position
	length - Maximum edge width
	mode - Soft edge mode: 'c' - cosine else gaussian

	Returns:
	Edge kernel value
	"""
	if mode.lower() == 'c':
		numpy.add(0.5, numpy.multiply(0.5, numpy.cos(numpy.pi * kernel_mask / float(length), out=kernel_mask), out=kernel_mask), out=kernel_mask)
	else:
		Q = -4.605170185988091
		numpy.exp(Q * (kernel_mask / float(length))**2, out=kernel_mask)


def soft_edge(img, length, mode='c', do_approx=False):
	"""
	Add a soft edge mask to a 2D/3D image binary image.

	Arguments:
	img - Input 2D or 3D image
	length - Length of the edge in pixel.
	mode - Mode of the mask 'c' - cosine; 'g' - gaussian

	Returns:
	Expanded mask
	"""
	if isinstance(img, EMAN2_cppwrap.EMData):
		img_data = EMAN2_cppwrap.EMNumPy.em2numpy(img)
		return_object = EMAN2_cppwrap.EMData(*img_data.shape)
		return_data = EMAN2_cppwrap.EMNumPy.em2numpy(return_object)
		out_eman = True
	else:
		img_data = img
		out_eman = False

	if length <= 0:
		return img.copy()

	# Get the mask shape for the soft edge kernel
	kernel_mask_dim = 2 * length + 1
	dimension = len(img_data.shape)
	mask_shape = tuple([kernel_mask_dim]*dimension)
	if dimension not in (2, 3):
		global_def.ERROR('morphology/soft_edge', 'Only 2D and 3D images are supported!', 1)

	# Create the outline for the array by erosing it once.
	# Pad the outline with the edge mask to avoid edge effects later.
	outline = img_data - snm.binary_erosion(img_data)
	outline = numpy.pad(outline, length + 1, mode='constant', constant_values=0)
	outline_index = numpy.where(outline == 1)

	# Fill the kernel with the soft edge values
	edge_norm = length**2
	cosine_falloff = 100
	if dimension == 2:
		x, y = numpy.ogrid[0:kernel_mask_dim, 0:kernel_mask_dim]
		kernel_mask = numpy.sqrt(((x - length)**2 + (y - length)**2) / float(edge_norm))*cosine_falloff
	elif dimension == 3:
		x, y, z = numpy.ogrid[0:kernel_mask_dim, 0:kernel_mask_dim, 0:kernel_mask_dim]
		kernel_mask = numpy.sqrt(((x - length)**2 + (y - length)**2 + (z - length)**2) / float(edge_norm))*cosine_falloff
	else:
		assert False

	if do_approx:
		numpy.add(kernel_mask, numpy.copysign(0.5, kernel_mask), kernel_mask)
		numpy.trunc(kernel_mask, kernel_mask)
	kernel_mask[kernel_mask >= cosine_falloff] = cosine_falloff
	fill_soft_edge_kernel_mask(kernel_mask, cosine_falloff, mode)

	# Replace the region around every outline pixel with the gaussian kernel.
	if dimension == 2:
		for x, y in zip(*outline_index):
			x_start = x - length
			x_stop = x + length + 1
			y_start = y - length
			y_stop = y + length + 1
			mask_slice = outline[
				x_start:x_stop,
				y_start:y_stop,
				]
			numpy.maximum(kernel_mask, mask_slice, mask_slice)
		outline = outline[
			length+1:outline.shape[0]-length-1,
			length+1:outline.shape[1]-length-1,
			]
	elif dimension == 3:
		for x, y, z in zip(*outline_index):
			x_start = x - length
			x_stop = x + length + 1
			y_start = y - length
			y_stop = y + length + 1
			z_start = z - length
			z_stop = z + length + 1
			mask_slice = outline[
				x_start:x_stop,
				y_start:y_stop,
				z_start:z_stop,
				]
			numpy.maximum(kernel_mask, mask_slice, mask_slice)
		outline = outline[
			length+1:outline.shape[0]-length-1,
			length+1:outline.shape[1]-length-1,
			length+1:outline.shape[2]-length-1,
			]
	else:
		assert False

	# Return a EMData object if an EMData object was the input
	combined_mask = numpy.maximum(img_data, outline)
	if out_eman:
		return_data[...] = combined_mask
	else:
		return_object = combined_mask
	return return_object


def binarize(img, minval = 0.0):
	"""
		Name
			binarize - create a binary image from an input image
		Input
			img: input image
			minval: value below which image pixels will be set to zero.
		Output
			binarized image.
	"""
	return img.process( "threshold.binary", {"value": minval} )

def collapse(img, minval = -1.0, maxval = 1.0):
	"""
		Name
			collapse - binarize image by setting to zero values within predefined range, and to one outside of this range
		Input
			input: input image
			minval: minimum bracket value (default is -1.0).
			maxval: maximum bracket value (default is 1.0).
		Output
			collapsed image.
	"""
	# for values between minval and maxval set to one, to zero outside
	return img.process( "threshold.binaryrange", {"low": minval, "high":maxval} )
			
def dilation(f, mask = None, morphtype="BINARY"):
	"""
		Name
			dilation - Calculate the dilated image.
		Input
			The first input image
			mask: The second input image used as the mask.
			The size of the mask has to be odd so that the center of mask can be well defined.
			The size of the mask should be smaller than the size of the first input image.
			morph_type: Type of the dilation
			BINARY is for binary dilation;
			GRAYLEVEL is for graylevel dilation.
		Output
			dilated image
	"""
	from EMAN2 import morph_type, filt_dilation_

	if not mask:
		from utilities import model_circle
		nx = f.get_xsize()
		ny = f.get_ysize()
		nz = f.get_zsize()
		if(nz == 1):	mask = model_circle(2,5,5)
		elif(nz >1):  mask = model_circle(2,5,5,5)
		else:  ERROR("Command does not work for 1D images","dilation",1)

	if morphtype=="BINARY":
		return filt_dilation_(f, mask, morph_type.BINARY)
	elif morphtype=="GRAYLEVEL":
		return filt_dilation_(f, mask, morph_type.GRAYLEVEL)
	else: ERROR("Unknown dilation type","dilation",1)

def erosion(f, mask = None, morphtype="BINARY"):
	"""
		Name
			erosion - Calculate the eroded image.
		Input
			The first input image
			mask: The second input image used as the mask.
			The size of the mask has to be odd so that the center of mask can be well defined.
			The size of the mask should be smaller than the size of the first input image.
			morph_type: Type of the erosion
			BINARY is for binary erosion (DEFAULT);
			GRAYLEVEL is for graylevel erosion.
		Output
			eroded image
	"""
	from EMAN2 import morph_type, filt_erosion_

	if not mask:
		from utilities import model_blank
		nx = f.get_xsize()
		ny = f.get_ysize()
		nz = f.get_zsize()
		from utilities import model_circle
		nx = f.get_xsize()
		ny = f.get_ysize()
		nz = f.get_zsize()
		if(nz == 1):	mask = model_circle(2,5,5)
		elif(nz >1):  mask = model_circle(2,5,5,5)
		else:  ERROR("Command does not work for 1D images","erosion",1)

	if morphtype=="BINARY":
		return filt_erosion_(f, mask, morph_type.BINARY)
	elif morphtype=="GRAYLEVEL":
		return filt_erosion_(f, mask, morph_type.GRAYLEVEL)
	else: ERROR("Unknown erosion type","erosion",1)

def invert(im):
	"""
	 Invert contrast of an image (while keeping the average)
	"""
	p = Util.infomask(im, None, True)[0]
	return ((-1.0*im) + 2*p)

#def compress(img, value = 0.0, frange=1.0):
#	return img.process( "threshold.compress", {"value": value, "range": frange } )

def expn(img, a = 1.0, b=0.0):
	"""
		Name
			expn - generate image whose pixels are generated of raising to a given power pixels of the input image
		Input
			image: input real image
		Output
			the output image whose pixels are given by o=ir
			r: exponent
	"""
	return img.process( "math.exp", {"low": 1.0/a, "high":b})

def power(img, x = 3.0):
	"""
		Name
			power - generate image whose pixels are generated of raising to a given power pixels of the input image
		Input
			image: input real image
		Output
			the output image whose pixels are given by o=ir
			x: power
	"""
	return img.process( "math.pow", {"pow": x})

def alog10(img):
	return img.process( "math.log")

def square_root(img):
	"""
		Name
			square_root - create image whose pixels will be square roots of pixels of the input image
		Input
			input image
		Output
			output image.
	"""
	[a,b,c,d] = Util.infomask(img, None, False)
	if(c<0.0):  ERROR("Cannot calculate square root of negative pixels","square_root",1)
	return img.process( "math.sqrt" )

def square(img):
	"""
		Name
			square - create image whose pixels will be squared values of pixels of the input image
		Input
			input image
		Output
			output image.
	"""
	return img.process( "math.squared")

def threshold(img, minval = 0.0):
	"""
		Name
			threshold - replace values below given threshold by zero
		Input
			img: input image
			minval: value below which image pixels will be set to zero. 
		Output
			thresholded image.
	"""
	return img.process( "threshold.belowtozero", {"minval": minval} )

def threshold_to_zero(img, minval = 0.0):
	"""
		Name
			threshold_to_zero - replace values below given threshold by zero and values above by (value-threshold)
		Input
			img: input image
			minval: value below which image pixels will be set to zero
		Output
			thresholded image.
	"""
	return img.process( "threshold.belowtozero_cut", {"minval": minval } )

def threshold_to_minval(img, minval = 0.0):
	"""
		Name
			threshold_to_minval - replace values below given threshold by the threshold value
		Input
			img: input image
			minval: value below which image pixels will be set to this value
		Output
			thresholded image.
	"""
	return img.process( "threshold.belowtominval", {"minval": minval } )

def threshold_outside(img, minval, maxval):
	"""
		Name
			threshold_outside - replace values outside given thresholds by respective threshold values
		Input
			img: input image
			minval: value below which image pixels will be set to this value.
			maxval: value above which image pixels will be set to this value.
	"""
	return img.process( "threshold.clampminmax", {"minval": minval, "maxval": maxval } )

def threshold_inside(img, minval, maxval):
	"""
		Name
			threshold_inside - replace values inside given thresholds by zeroes
		Input
			img: input image
			minval, maxval: image pixels that have values between these thresholds will be set to zero.
		WARNING: This function is written in python and thus very slow for large images/volumes
	"""
	im = img.copy()
	nx = im.get_xsize()
	ny = im.get_ysize()
	nz = im.get_zsize()
	for z in range(nz):
		for y in range(ny):
			for x in range(nx):
				q = im.get_value_at(x,y,z)
				if( q>minval and q < maxval): im.set_value_at(x,y,z,0.0)
	return im

def threshold_maxval(img, maxval = 0.0):
	"""
		Name
			threshold_maxval - replace values above given threshold by the threshold value
		Input
			img: input image
			maxval: value above which image pixels will be set to this value
		Output
			thresholded image.
	"""
	st = Util.infomask(img, None, True)
	return img.process( "threshold.clampminmax", {"minval": st[2], "maxval": maxval } )

def notzero(img):
	"""
		Name
			notzero - replace values that are not zero by 1.0
		Input
			img: input image
		Output
			binary image.
	"""
	return img.process( "threshold.notzero" )

def linchange(a, fct):
	"""
	reinterpolate a line given as a list by a factor of fct.
	Useful for adjusting 1D power spectra, uses linear interplation
	"""
	fctf = float(fct)
	n = len(a)
	m = int(n*fctf+0.5)
	o = [0.0]*m
	for i in range(m):
		x = i/fctf
		j = min(int(x), n-2)
		dx = x-j
		o[i] = (1.0-dx)*a[j] + dx*a[j+1]
	return o

## Fitting of "ideal" 3D FSC, as elaborated in Resolution, Meth Enz, 2010, Eq.3.25
"""
def ideal_fsc(n,a,N,B,g):
	#  Note number of image N and scaling factor for SSNR are equivalent
	o=[]
	for i in xrange(n/2):
		s=float(i)/n/a
		#  For 2D float(i) has to be replaced by 1.0.
		o.append(N*exp(-B*s*s/4.)/(N*exp(-B*s*s/4.)+float(i)/g))
	return o

def difsc(args,data):
	z=ideal_fsc(360,1.34,50000,args[0],args[1])
	vv = 0.0
	for l in xrange(len(z)):
		vv -= (data[0][l] - z[l])**2
	print(args,vv)
	return vv

##amoeba(arg, [1.0,0.1], difsc, data=data)
##z=ideal_fsc(360,1.34,50000,202.,.05);write_text_file(z,"zz.txt")
"""


## CTF related functions
def rotavg_ctf(img, defocus, Cs, voltage, Pixel_size, amp = 0.0, ang = 0.0):
	"""1D rotational average of a 2D power spectrum (img)
	   based on estimated CTF parameters, including astigmatism amplitude and angle
	"""
	from math import sqrt,atan2,pi,sin,radians
	defc = defocus*10000
	astigmamp = amp*10000
	lam = 12.398/sqrt(voltage*(1022.0+voltage))
	angrad = radians(ang)
	nx = img.get_xsize()
	lr = [0.0]*2*(nx//2+1)
	cnt = [0.0]*2*(nx//2+1)
	nc = nx//2
	nr2 = nc*nc + 1
	if( Cs == 0.0 ):
		for ix in range(nx):
			x = ix - nc
			for iy in range(nx):
				y = iy - nc
				r2 = x*x + y*y
				if( r2 < nr2 ):
					dfa = defc - astigmamp/2*sin(2*(atan2(y,x) + angrad))
					try:
						u = sqrt( dfa/defc ) * sqrt(r2)
						iu = int(u)
						du = u - iu
						lr[iu]    += (1.0-du)*img.get_value_at(ix,iy)
						lr[iu+1]  +=       du*img.get_value_at(ix,iy)
						cnt[iu]   += 1.0-du
						cnt[iu+1] +=     du
					except: pass
	else:
		Cst = Cs*1.e7
		for ix in range(nx):
			x = ix - nc
			for iy in range(nx):
				y = iy - nc
				r2 = x*x + y*y
				if( r2 < nr2 ):
					s = sqrt(r2)/(nc*2*Pixel_size)
					dfa = defc - astigmamp/2*sin(2*(atan2(y,x) + angrad))
					#u = sqrt(r2)*sqrt(1.0 -  astigmamp/2./defc*sin(2*(atan2(y,x) - angrad)))
					#print ix,iy,sqrt(r2),defc,dfa,lam,s,u
					#print  ix,iy,sqrt(r2),defc**2 + Cst**2*lam**4*s**4 - 2*Cst*lam**2*s**2*dfa
					#print  defc
					#print  defc - sqrt( defc**2 + Cst**2*lam**4*s**4 - 2*Cst*lam**2*s**2*dfa)
					try:
						u = sqrt( Cst*(defc - sqrt( defc**2 + Cst**2*lam**4*s**4 - 2*Cst*lam**2*s**2*dfa)))/(Cst*lam) * nc*2*Pixel_size
						iu = int(u)
						du = u - iu
						lr[iu]    += (1.0-du)*img.get_value_at(ix,iy)
						lr[iu+1]  +=       du*img.get_value_at(ix,iy)
						cnt[iu]   += 1.0-du
						cnt[iu+1] +=     du
					except:  pass
	for ix in range(nc):  lr[ix] /= max(cnt[ix], 1.0)
	return lr[:nc]


def ctf_1d(nx, ctf, sign = 1, doabs = False):
	"""
		Generate a list of 1D CTF values 
		Input
			nx: image size to which CTF will be applied.
			ctf: CTF object created using generate_ctf
			sign: sign of the CTF.
		Output
			a list of CTF values.
	"""
	dict = ctf.to_dict()
	dz = dict["defocus"]
	cs = dict["cs"]
	voltage = dict["voltage"]
	pixel_size = dict["apix"]
	bfactor = dict["bfactor"]
	ampcont = dict["ampcont"]


	ctf_1 = []
	scl    = 1./pixel_size/nx
	length = int(1.41*float(nx/2)) + 1
	ctf_1 = [0.0]*length
	if doabs:
		for i in range(length): ctf_1[i] = abs(Util.tf(dz, i*scl, voltage, cs, ampcont, bfactor, sign))
	else:
		for i in range(length): ctf_1[i] = Util.tf(dz, i*scl, voltage, cs, ampcont, bfactor, sign)
	return ctf_1

def ctf_2(nx, ctf):
	"""
		Generate a list of 1D CTF^2 values 
		Input
			nx: image size to which CTF will be applied.
			ctf: ctf object created using generate_ctf
		Output
			a list of CTF2 values.
	"""
	dict       = ctf.to_dict()
	dz         = dict["defocus"]
	cs         = dict["cs"]
	voltage    = dict["voltage"]
	pixel_size = dict["apix"]
	b_factor   = dict["bfactor"]
	ampcont    = dict["ampcont"]

	ctf_2  = []
	scl    = 1.0/pixel_size/nx
	length = int(1.7321*float(nx/2)) + 2
	ctf_2 = [0.0]*length
	for i in range(length):
		ctf_val = Util.tf(dz, i*scl, voltage, cs, ampcont, b_factor)
		ctf_2[i] = ctf_val*ctf_val
	return ctf_2


def ctf_img(nx, ctf, sign = 1, ny = 0, nz = 1):
	"""
		Generate a 1-2-3-D complex image containing the CTF.
	 	Default is 2D output.
	  	Input
			nx: x image size.
			ctf: ctf object, see CTF_info for description.
			sign: sign of the CTF.
			ny: y image size
			nz: z image size 
		Output
			ctfimg: complex image containing CTF.
	"""
	dict = ctf.to_dict()
	dz = dict["defocus"]
	cs = dict["cs"]
	voltage = dict["voltage"]
	pixel_size = dict["apix"]
	b_factor = dict["bfactor"]
	ampcont = dict["ampcont"]
	dza = dict["dfdiff"]
	azz = dict["dfang"]

	if(ny < 1):  ny = nx
	return  Util.ctf_img(nx, ny, nz, dz, pixel_size, voltage, cs, ampcont, b_factor, dza, azz, sign)

def ctf_img_real(nx, ctf, sign = 1, ny = 0, nz = 1):
	"""
		Generate a 1-2-3-D real image containing the CTF.
	 	Default is 2D output.
	  	Input
			nx: x image size.
			ctf: ctf object, see CTF_info for description.
			sign: sign of the CTF.
			ny: y image size
			nz: z image size 
		Output
			ctfimg: real image containing CTF, x-size half of the complex
	"""
	dict = ctf.to_dict()
	dz = dict["defocus"]
	cs = dict["cs"]
	voltage = dict["voltage"]
	pixel_size = dict["apix"]
	b_factor = dict["bfactor"]
	ampcont = dict["ampcont"]
	dza = dict["dfdiff"]
	azz = dict["dfang"]

	if(ny < 1):  ny = nx
	return  Util.ctf_img_real(nx, ny, nz, dz, pixel_size, voltage, cs, ampcont, b_factor, dza, azz, sign)


def ctf_rimg(nx, ctf, sign = 1, ny = 0, nz = 1):
	"""
		Generate a 1-2-3-D real image containing the CTF.
	 	Default is 2D output.
	  	Input
			nx: x image size.
			ctf: ctf object, see CTF_info for description.
			sign: sign of the CTF, if sign=0 compute |CTF|
			ny: y image size
			nz: z image size 
		Output
			ctfimg: image containing CTF^2.
	"""
	dict       = ctf.to_dict()
	dz         = dict["defocus"]
	cs         = dict["cs"]
	voltage    = dict["voltage"]
	pixel_size = dict["apix"]
	b_factor   = dict["bfactor"]
	ampcont    = dict["ampcont"]
	dza        = dict["dfdiff"]
	azz        = dict["dfang"]

	if(ny < 1):  ny = nx
	return  Util.ctf_rimg(nx, ny, nz, dz, pixel_size, voltage, cs, ampcont, b_factor, dza, azz, sign)

def ctf2_rimg(nx, ctf, sign = 1, ny = 0, nz = 1):
	"""
		Generate a 1-2-3-D real image containing the CTF^2.
	 	Default is 2D output.
	  	Input
			nx: x image size.
			ctf: ctf object, see CTF_info for description.
			sign: sign of the CTF.
			ny: y image size
			nz: z image size 
		Output
			ctfimg: image containing CTF^2.
	"""
	dict       = ctf.to_dict()
	dz         = dict["defocus"]
	cs         = dict["cs"]
	voltage    = dict["voltage"]
	pixel_size = dict["apix"]
	b_factor   = dict["bfactor"]
	ampcont    = dict["ampcont"]
	dza        = dict["dfdiff"]
	azz        = dict["dfang"]

	if(ny < 1):  ny = nx
	return  Util.ctf2_rimg(nx, ny, nz, dz, pixel_size, voltage, cs, ampcont, b_factor, dza, azz, sign)


def ctflimit(nx, defocus, cs, voltage, pix):
	"""
	 Find aliasing limit given 
	   nx - window size in pixels
	   defocus
	   cs
	   voltage
	   pix - pixel size in A
	 Ouput:
	  Fourier pixel at which aliasing will occur and corresponding Fourier frequency
	  Note for window size nx maximum Fourier pixel number is nx/2.
	"""
	import numpy as np
	def ctfperiod(defocus, Cs, lam, freq):
		# find local "period" T by solving fourth order polynomial resulting from equation:
		#  sin( 2pi (gamma(freq) + 1) ) = sin( 2pi (gamma(freq+T) )
		cis = Cs*1.e7
		A = 0.5*defocus*10000.0*lam
		B = 0.25*cis*lam*lam*lam
		f2 = freq*freq
		"""
		for i in xrange(800):
			ff = freq+(i-400)*0.00002
			print  ff,Util.tf(defocus, ff, voltage, Cs, 10., 0.0, 1.0)
		"""
		rot = np.roots([B, 4*B*freq, 6*B*f2-A, 4*B*f2*freq-2*A*freq, -1.0])
		#print np.roots([A,2*A*freq,1.0]),-freq-np.sqrt(f2/2-1./A),-freq+np.sqrt(f2-1./A)
		#print  Util.tf(defocus, freq, voltage, Cs, 10., 0.0, 1.0),Util.tf(defocus, freq+min(rot), voltage, Cs, 10., 0.0, 1.0)
		return min(abs(rot))

	n = nx//2+1
	#  Width of Fourier pixel
	fwpix = 1./(2*pix)/n
	
	# Fourier cycle
	fcycle = 1./(2*fwpix)
	#Fourier period
	fper = 1.0/fcycle
	#print "Image size %6d,   pixel size  %7.4f  Width of Fourier pixel %7.5f   Fourier period  %8.5f "%(nx,pix,fwpix,fper)

	
	#CTF
	lam = 12.398/np.sqrt(voltage*(1022.0+voltage))  #  All units in A
	z1 = defocus*10000.0
	ctfper = ctfperiod(defocus, cs, lam, 1./(2*pix))
	#print " At Nyquist, the CTF period is ",ctfper
	for ii in range(n-1,1,-1):
		#print ii
		xr = ii/float(n-1)/(2*pix)
		ctfper = ctfperiod(defocus, cs, lam, xr)
		#print ii,xr,ctfper
		if(ctfper >  fper):
			#print  " Limiting frequency is:",xr,"  limiting resolution is:",1.0/xr
			return  int(xr/fwpix+0.5),xr
	return nx//2,1.0/(2*pix)

def compare_ctfs(nx, ctf1, ctf2):
	sign = 1.0
	dict       = ctf1.to_dict()
	dz         = dict["defocus"]
	cs         = dict["cs"]
	voltage    = dict["voltage"]
	pixel_size = dict["apix"]
	b_factor   = dict["bfactor"]
	ampcont    = dict["ampcont"]
	dza        = dict["dfdiff"]
	azz        = dict["dfang"]
	cim1 = Util.ctf_rimg(nx, 1, 1, dz, pixel_size, voltage, cs, ampcont, b_factor, dza, azz, sign)
	dict       = ctf2.to_dict()
	dz         = dict["defocus"]
	cs         = dict["cs"]
	voltage    = dict["voltage"]
	pixel_size2= dict["apix"]
	b_factor   = dict["bfactor"]
	ampcont    = dict["ampcont"]
	dza        = dict["dfdiff"]
	azz        = dict["dfang"]
	if(pixel_size != pixel_size2):
		ERROR("CTFs have different pixel sizes, pixel size from the first one used", "compare_ctfs", 0)
	cim2 = Util.ctf_rimg(nx, 1, 1, dz, pixel_size, voltage, cs, ampcont, b_factor, dza, azz, sign)
	for i in range(nx//2,nx):
		if(cim1.get_value_at(i)*cim2.get_value_at(i) < 0.0):
			limi = i-nx//2
			break
	return limi, pixel_size*nx/limi
	

###----D-----------------------		
def defocus_env_baseline_fit(roo, i_start, i_stop, nrank, iswi):
	"""
		    iswi = 2 using polynomial n rank to fit envelope function
			iswi = 3 using polynomial n rank to fit background
	"""
	TMP = imf_params_cl1(roo[i_start:i_stop], nrank, iswi)
	curve   = [0]*len(roo)
	curve[i_start:i_stop] = TMP[1][:i_stop-i_start]
	return curve

def defocus_get(fnam_roo, volt=300, Pixel_size=1, Cs=2, wgh=.1, f_start=0, f_stop=-1, docf="a", skip="#", round_off=1, nr1=3, nr2=6):
	"""
	
		1.Estimating envelope function and baseline noise using constrained simplex method
		  so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get 
		   defocus which matches the extracted CTF imprints 
	"""

	from math 	import sqrt, atan
	from utilities 	import read_text_row
	roo     = []
	res     = []
	if(docf == "a"):
		TMP_roo = read_text_row(fnam_roo, "a", skip)
		for i in range(len(TMP_roo)):	roo.append(TMP_roo[i][1])
	else:
		TMP_roo=read_text_row(fnam_roo,"s",";")
		for i in range(len(TMP_roo)):	roo.append(TMP_roo[i][2])
	Res_roo = []
	Res_TE  = []	
	if f_start == 0 : 	i_start=0
	else: 			i_start=int(Pixel_size*2.*len(roo)/f_start)
	if f_stop <= i_start : 	i_stop=len(roo)
	else: 			i_stop=int(Pixel_size*2.*len(roo)/f_stop)

	TE  = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr1), 4)
	Pn1 = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr2), 3)
	Res_roo = []
	Res_TE  = []	
	for i in range(len(roo)):
		Res_roo.append(roo[i] - Pn1[i])
		Res_TE.append( TE[i]  - Pn1[i])
	#
	defocus=defocus_guess(Res_roo, Res_TE, volt, Cs, Pixel_size, wgh, i_start, i_stop, 2, round_off)
	del    roo
	del    TE
	del    Pn1
	del    Res_TE
	del    Res_roo	
	return defocus

def defocus_gett(roo, voltage=300.0, Pixel_size=1.0, Cs=2.0, wgh=0.1, f_start=0.0, f_stop=-1.0, round_off=1.0, nr1=3, nr2=6, parent=None):
	"""
	
		1. Estimate envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get 
		   defocus which matches the extracted CTF imprints 
	"""
	from utilities  import generate_ctf
	#print "CTF params:", voltage, Pixel_size, Cs, wgh, f_start, f_stop, round_off, nr1, nr2, parent
	if f_start == 0 : 	i_start = 0
	else: 			    i_start = int(Pixel_size*2.*len(roo)*f_start)
	if f_stop <= f_start : 	i_stop  = len(roo)
	else:                   i_stop  = min(len(roo), int(Pixel_size*2.*len(roo)*f_stop))

	#print "f_start, i_start, f_stop, i_stop:", f_start, i_start, f_stop, i_stop
	TE  = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr1), 4)
	Pn1 = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr2), 3)
	Res_roo = []
	Res_TE  = []
	for i in range(len(roo)):
		Res_roo.append(roo[i] - Pn1[i])
		Res_TE.append( TE[i]  - Pn1[i])

	defocus = defocus_guess(Res_roo, Res_TE, voltage, Cs, Pixel_size, wgh, i_start, i_stop, 2, round_off)

	nx  = int(len(Res_roo)*2)
	ctf = ctf_2(nx, generate_ctf([defocus, Cs, voltage, Pixel_size, 0.0, wgh]))
	if (parent is not None):
		parent.ctf_data=[roo, Res_roo, Res_TE]
		parent.i_start = i_start
		parent.i_stop = i_stop
		from utilities import write_text_file
		write_text_file([roo, Res_roo, Res_TE, ctf], "procpw.txt")
	else:
		from utilities import write_text_file
		write_text_file([roo, Res_roo, Res_TE, ctf], "procpw.txt")
	return defocus

def defocus_get_Eudis(fnam_roo, volt=300, Pixel_size=1, Cs=2, wgh=.1, f_start=0, f_stop=-1, docf="a" ,skip="#", round_off=1, nr1=3, nr2=6):
	"""
		1. Estimating envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get 
		   defocus which matches the extracted CTF imprints
		3. It returns Euclidean distance for defocus selection 
	"""
	from math 	import sqrt, atan
	from utilities 	import read_text_row, generate_ctf
	roo     = []
	res     = []
	if docf == "a":
		TMP_roo = read_text_row(fnam_roo, "a", skip)
		for i in range(len(TMP_roo)): # remove first record
			roo.append(TMP_roo[i][1])
	else:
		skip = ";"
		TMP_roo = read_text_row(fnam_roo, "s", skip)
		for i in range(len(TMP_roo)): # remove first record
			roo.append(TMP_roo[i][2])
	Res_roo = []
	Res_TE  = []	
	if f_start == 0 : 	i_start=0
	else: 			i_start=int(Pixel_size*2.*len(roo)/f_start)
	if f_stop <= i_start :	i_stop=len(roo)
	else: 			i_stop=int(Pixel_size*2.*len(roo)/f_stop)	
	TE  = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr1), 4)
	Pn1 = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr2), 3)
	Res_roo = []
	Res_TE  = []
	for i in range(len(roo)):
		Res_roo.append( roo[i] - Pn1[i] )
		Res_TE.append(  TE[i]  - Pn1[i] )
	#
	defocus=defocus_guess(Res_roo, Res_TE, volt, Cs, Pixel_size, wgh, i_start, i_stop, 2, round_off)
	nx  = int(len(roo)*2)
	ctf = ctf_2(nx, generate_ctf([defocus,Cs,voltage,Pixel_size, 0.0, wgh]))
	for i in range(len(Res_TE)):
		ctf[i]=ctf[i]*Res_TE[i]
	dis = defocus_L2_euc(ctf, Res_roo, i_start, i_stop)
	return [defocus, dis]

def defocus_L2_euc(v1,v2, ist,istp):
	from math import sqrt
	dis    = 0.0
	pw_sum = 0.0
	if ist == istp :	ERROR("No pw2 curve is included  ", "defocus_L2_euc", 0)
	else:			tfeq = istp-ist
	for i in range(ist,istp,1):
		dis+=    (v1[i]-v2[2])**2
		pw_sum+= (v1[i])**2
	if pw_sum <= 0:		ERROR("negative or zero power ", "defocus_L2_euc", 1)
	if dis    <= 0:		ERROR("bad fitting, change options settings and try again  ", "defocus_L2_euc", 0)
	else:
		res = sqrt(dis)/sqrt(pw_sum)/tfeq	
		return res

def defocus_guess(Res_roo, Res_TE, volt, Cs, Pixel_size, ampcont=10.0, istart=0, istop=-1, defocus_estimation_method=2, round_off=1, dz_low=1000., dz_high=200000., nloop=100):
	"""
		Use specified frequencies area (istart-istop)to estimate defocus
		1.  The searching range is limited to dz_low (.1um) ~ dz_high (20 um).
		    The user can modify this limitation accordingly
		2.  changing nloop can speed up the estimation
		3.  defocus_estimation_method = 1  use squared error
		    defocus_estimation_method = 2  use normalized inner product
		Input:
		  Res_roo - background-subtracted Power Spectrum
		  Res_TE  - background-subtracted Envelope
	"""
	
	from math import sqrt
	from utilities import generate_ctf

	if istop <= istart : 			istop=len(Res_roo)
	step = (dz_high-dz_low)/nloop
	if step > 10000.   : 			step     =  10000.     # Angstrom

	xval_e = 0.0
	for ifreq in range(len(Res_TE)):
		xval_e += Res_TE[ifreq]**2
	if (xval_e == 0.0): return xvale_e  #  This is strange, returns defocus=0.

	if round_off >= 1: 			cut_off  =  1.
	else: 					    cut_off  =  round_off # do extreme fitting

	length = len(Res_roo)
	nx     = int(length*2)
	if defocus_estimation_method == 1 : 	diff_min =  1.e38
	else: 					                diff_min = -1.e38
	while (step >= cut_off):
		for i_dz in range(nloop):
			dz     = dz_low + step*i_dz
			ctf    = ctf_2(nx, generate_ctf([dz, Cs, volt, Pixel_size, 0.0, ampcont]))
			diff   = 0.0
			if defocus_estimation_method == 1:
				for ifreq in range(istart, istop, 1):
					diff += (ctf[ifreq]*Res_TE[ifreq] - Res_roo[ifreq])**2
					if diff < diff_min :
						defocus  = dz
						diff_min = diff
			else:
				diff  = 0.0
				sum_a = 0.0
				sum_b = 0.0
				for ifreq in range(istart, istop, 1):
					xval   =  ctf[ifreq]*Res_TE[ifreq]
					diff  +=  Res_roo[ifreq]*xval
					sum_a +=  Res_roo[ifreq]**2
					sum_b +=  xval**2
				diff /= (sqrt(sum_a*sum_b)*( istop - istart + 1 ))
				if diff > diff_min :
					defocus  = dz
					diff_min = diff

		dz_low = defocus-step*2
		if( dz_low < 0 ): 	dz_low=0.0
		dz_high = defocus + step*2
		step /= 10.
	defocus = int( defocus/round_off )*round_off
	return defocus


def defocus_guess1(Res_roo, Res_TE, volt, Cs, Pixel_size, ampcont=10.0, istart=0, istop=-1, defocus_estimation_method=2, round_off=1, dz_low=1000., dz_high=200000., nloop=100):
	"""
		Use specified frequencies area (istart-istop) to estimate defocus from crossresolution curve
		1.  The searching range is limited to dz_low (.1um) ~ dz_high (20 um).
		    The user can modify this limitation accordingly
		2.  changing nloop can speed up the estimation
		3.  defocus_estimation_method = 1  use squared error
		    defocus_estimation_method = 2  use normalized inner product
		Input:
		  Res_roo - Cross-resolution
		  Res_TE  - Envelope
	"""
	
	from math import sqrt
	from utilities import generate_ctf
	from morphology import ctf_1d

	if istop <= istart : 			istop=len(Res_roo)
	step = (dz_high-dz_low)/nloop
	if step > 10000.   : 			step     =  10000.     # Angstrom

	xval_e = 0.0
	for ifreq in range(len(Res_TE)):
		xval_e += Res_TE[ifreq]**2
	if (xval_e == 0.0): return xvale_e  #  This is strange, returns defocus=0.

	if round_off >= 1: 			cut_off  =  1.
	else: 					    cut_off  =  round_off # do extreme fitting

	length = len(Res_roo)
	nx     = int(length*2)
	if defocus_estimation_method == 1 : 	diff_min =  1.e38
	else: 					                diff_min = -1.e38
	while (step >= cut_off):
		for i_dz in range(nloop):
			dz     = dz_low + step*i_dz
			ctf    = ctf_1d(nx, generate_ctf([dz, Cs, volt, Pixel_size, 0.0, ampcont]))
			diff   = 0.0
			if defocus_estimation_method == 1:
				for ifreq in range(istart, istop, 1):
					diff += (ctf[ifreq]*Res_TE[ifreq] - Res_roo[ifreq])**2
					if diff < diff_min :
						defocus  = dz
						diff_min = diff
			else:
				diff  = 0.0
				sum_a = 0.0
				sum_b = 0.0
				for ifreq in range(istart, istop, 1):
					xval   =  ctf[ifreq]*Res_TE[ifreq]
					diff  +=  Res_roo[ifreq]*xval
					sum_a +=  Res_roo[ifreq]**2
					sum_b +=  xval**2
				diff /= (sqrt(sum_a*sum_b)*( istop - istart + 1 ))
				if diff > diff_min :
					defocus  = dz
					diff_min = diff

		dz_low = defocus-step*2
		if( dz_low < 0 ): 	dz_low=0.0
		dz_high = defocus + step*2
		step /= 10.
	defocus = int( defocus/round_off )*round_off
	return defocus

def defocus_get_fast(indir, writetodoc="w", Pixel_size=1, volt=120, Cs=2, wgh=.1, round_off=100, dz_max0=50000, f_l0=30, f_h0=5, nr_1=5, nr_2=5, prefix="roo", docf="a",skip="#", micdir="no", print_screen="p"):
	"""
		Estimate defocus using user supplied 1D power spectrum area
		writetodoc="a" return the estimated defoci in a list, and write them down also in a text file
		writetodoc="l" output estimated defocus in a list
		writetodoc="w" output estimated defocus in a text file
	"""
	import os
	import types
	from utilities import set_arb_params, get_image
	if writetodoc[0]   != "a" and writetodoc[0]   != "l" and writetodoc[0] != "a": 	writetodoc= "a"
	if print_screen[0] != "p" and print_screen[0] != "n"			     : 	print_screen = "n"
	if os.path.exists(indir) == False: 	ERROR("roodir doesn't exist", "defocus_get_fast",1)
	ctf_dicts = ["defocus", "Pixel_size", "voltage", "Cs", "amp_contrast", "B_factor", "sign"] 
	flist = os.listdir(indir)
	res   = []
	f_l   = f_l0
	f_h   = f_h0
	if(f_l <= 1 and f_l> 0)	:
		 f_l = 1./f_l
		 f_h = 1./f_h
	if(f_h > f_l or f_l <= 0 or f_h <= 0): 
		f_h  = 8
		f_l  = 30
	if nr_1       <=  1 	:	nr_1      =  5.
	if nr_2       <=  1 	:	nr_2      =  5.
	if round_off  <=  0	: 	round_off =  100.
	if dz_max0    <=  1	: 	dz_max0   =  100000.
	dz_max=dz_max0
	if writetodoc[0] == "w" or writetodoc == "a":
		fdefo_nam = "defocus.txt"
		out       =  open(fdefo_nam, "w")
		out.write("#defocus: %s\n")
	defocus = 0
	ncount  = 0	
	for i, v in enumerate(flist):
		(fnam, fext) = os.path.splitext(v)
		if(fnam[0:len(prefix)] == prefix):
			ncount   += 1
			fnam_root = fnam[len(prefix):]
			nr1       = int(nr_1)
			nr2       = int(nr_2)
			istart    = int(f_l)
			istop     = int(f_h)
			fnam_roo  = os.path.join(indir, v)
			defocus = defocus_get(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
			if(defocus > dz_max):
				while(nr1 <= 7 or nr2 <= 7):
					nr1 += 1
					nr2 += 1
					defocus = defocus_get(fnam_roo, volt, Pixel_size, Cs, wgh,istart, istop, docf,skip, round_off, nr1, nr2)
					#if(print_screen[0] == "p" or print_screen[0] == "P" ): print "defocus",defocus,"Euclidean distance",dis,"starting feq",istart,"stop freq",istop,"P R E", nr1,"P R B", nr2
					if(defocus<dz_max): break
			if(defocus > dz_max):
				while(nr1 >= 2 and nr2 >= 2):
					nr1 -= 1
					nr2 -= 1
					defocus = defocus_get(fnam_roo, volt,Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
					#if(print_sreen[0] == "p" or print_screen=="P"): print "defocus",defocus,"Euclidean distance",dis,"starting feq",istart,"stop freq",istop,"P R E", nr1,"P R B", nr2
					if(defocus < dz_max): break
			if(defocus > dz_max):
				while(istart > istop):
					nr1    =  5
					nr2    =  5
					istart -=.5
					defocus = defocus_get(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf,skip, round_off, nr1, nr2)
					#if(print_screen[0] == "p" or print_screen == "P"): print "defocus",defocus,"Euclidean distance",dis,"starting feq",istart,"stop freq",istop,"P R E", nr1,"P R B", nr2
					if(defocus < dz_max): break
			if(defocus > dz_max):
				while(istart > istop):
					nr1     = 5										    	
					nr2     = 5
					istop  += 0.5
					defocus = defocus_get(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
					if(print_screen == "p" or print_screen == "P"): print("defocus",defocus,"Euclidean distance", dis, "starting feq", istart, "stop freq", istop,"P R E", nr1,"P R B", nr2)
					if(defocus < dz_max): 				break
			if(defocus >= dz_max): 					ERROR("defocus_get_fast fails at estimating defocus", fnam, action = 0)
			print("", flist[i], '%5d'%(defocus)) 	# screen output, give the user a general impression about estimated defoci
			if(writetodoc[0] == "w" or writetodoc[0] != "l"):	out.write("%d\t%5d\t%s\n" % (ncount,defocus,flist[i]))
			if(writetodoc[0] == "l"):				res.append(defocus)
			if type(micdir) is bytes : 
				ctf_param = [defocus, Pixel_size, volt, Cs, wgh, 0, 1]
				mic_name  = os.path.join(micdir,""+ fnam_root+ ".hdf")
				if os.path.exists(mic_name) :
					e = get_image (mic_name)
					U______set_arb_params(e, ctf_param, ctf_dicts)  # THIS IS INCORRECT< PLEASE CHANGE
					e.write_image(mic_name,0, EMUtil.ImageType.IMAGE_HDF, True)
					print("ctf parameters is written back into headers of ", mic_name)
				#else :  print  mic_name, " Not found"
	if(len(res) == 0 and  writetodoc == "l" ):				ERROR("No input file is found, check the input directory of file prefix", indir, 1)
	else:
		if(writetodoc[0] == "a"):
			out.close()
			return res
		if(writetodoc[0] == "l"): 	return res
		if(writetodoc[0] == "w"): 	out.close()

def defocus_get_fast_MPI(indir, writetodoc="w", Pixel_size=1, volt=120, Cs=2, wgh=.1, round_off=100, dz_max0=50000, f_l0=30, f_h0=5, nr_1=5, nr_2=5, prefix_of_="roo", docf="a",skip="#",print_screen="no"):
	"""
		Estimate defocus using user defined 1D power spectrum area
		writetodoc="a" return the estimated defoci in a list, and write them down also in a text file
		writetodoc="l" output estimated defocus in a list
		writetodoc="w" output estimated defocus in a text file
	"""
	import os
	import sys
	if os.path.exists(indir) == False: 	ERROR("roodir doesn't exist", "defocus_get_fast",1)
	flist = os.listdir(indir)
	for i, v in enumerate(flist):
		micname                  = os.path.join(indir,v)
		(filename, filextension) = os.path.splitext(v)
		if(filename[0:len(prefix_of_)] == prefix_of_):
			mic_name_list.append(micname)
			nima += 1
	if nima < 1: 	ERROR("No  is found, check either directory or prefix of s is correctly given","pw2sp",1)
	
	sys.argv       = mpi_init(len(sys.argv),sys.argv)
	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	#  chose a random node as a main one...
	main_node      = 0
	if(myid == 0): main_node = randint(0,number_of_proc-1)
	main_node      = mpi_bcast(main_node, 1, MPI_INT, 0, MPI_COMM_WORLD)

	if(myid == main_node):
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)
	if(number_of_proc <= nima ):	nimage_per_node = nima/number_of_proc
	else: 				nimage_per_node = 1 
	image_start    = myid * nimage_per_node
	if(myid == number_of_proc-1):  image_end = nima
	else:                          image_end = image_start + nimage_per_node
	
	if writetodoc[0]   != "a" and writetodoc[0]   != "l" and writetodoc[0] != "a": 	writetodoc= "a"
	if print_screen[0] != "p" and print_screen[0] != "n"			     : 	print_screen = "n"
	res   = []
	f_l   = f_l0
	f_h   = f_h0
	if(f_l <= 1 and f_l> 0)	:
		 f_l = 1./f_l
		 f_h = 1./f_h
	if(f_h > f_l or f_l <= 0 or f_h <= 0): 
		f_h  = 8
		f_l  = 30
	if nr_1       <=  1 	:	nr_1      =  5.
	if nr_2       <=  1 	:	nr_2      =  5.
	if round_off  <=  0	: 	round_off =  100.
	if dz_max0    <=  1	: 	dz_max0   =  100000.
	dz_max = dz_max0
	if writetodoc[0] == "w" or writetodoc == "a":
		fdefo_nam = "defocus.txt"
		out       =  open(fdefo_nam, "w")
		out.write("#defocus: %s\n")
	defocus = 0
	ncount  = 0
	nr1	= int(nr_1)
	nr2	= int(nr_2)
	istart	= int(f_l )
	istop	= int(f_h )
	for i in range(image_start,image_end):
		filename=mic_name_list[i] 
		print('%-15s%-30s'%("s # ",filename))
		(f_nam, filextension) = os.path.splitext(filename)
		fnam_roo     = "particle_"+f_nam[len(prefix_of_)+len(indir)+2:]+filextension	
#	for i, v in enumerate(flist):
		ncount   += 1
		defocus = defocus_get(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
		if(defocus > dz_max):
			while(nr1 <= 7 or nr2 <= 7):
				nr1 += 1
				nr2 += 1
				defocus = defocus_get(fnam_roo, volt, Pixel_size, Cs, wgh,istart, istop, docf,skip, round_off, nr1, nr2)
				if(print_screen[0] == "p" or print_screen[0] == "P" ): print("defocus",defocus,"Euclidean distance",dis,"starting feq",istart,"stop freq",istop,"P R E", nr1,"P R B", nr2)
				if(defocus<dz_max): break
		if(defocus > dz_max):
			while(nr1 >= 2 and nr2 >= 2):
				nr1 -= 1
				nr2 -= 1
				defocus = defocus_get(fnam_roo, volt,Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
				if(print_sreen[0] == "p" or print_screen=="P"): print("defocus",defocus,"Euclidean distance",dis,"starting feq",istart,"stop freq",istop,"P R E", nr1,"P R B", nr2)
				if(defocus < dz_max): break
		if(defocus > dz_max):
			while(istart > istop):
				nr1    =  5
				nr2    =  5
				istart -=.5
				defocus = defocus_get(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf,skip, round_off, nr1, nr2)
				if(print_screen[0] == "p" or print_screen == "P"): print("defocus",defocus,"Euclidean distance",dis,"starting feq",istart,"stop freq",istop,"P R E", nr1,"P R B", nr2)
				if(defocus < dz_max): break
		if(defocus > dz_max):
			while(istart > istop):
				nr1     = 5										    	
				nr2     = 5
				istop  += 0.5
				defocus = defocus_get(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
				if(print_screen == "p" or print_screen == "P"): print("defocus",defocus,"Euclidean distance", dis, "starting feq", istart, "stop freq", istop,"P R E", nr1,"P R B", nr2)
				if(defocus < dz_max): 				break
		if(defocus >= dz_max): 					ERROR("defocus_get_fast fails at estimating defocus", fnam, action = 0)
		print("", flist[i], '%10.3g'(defocus)) 	# screen output, give the user a general impression about estimated defoci
		if(writetodoc[0] == "w" or writetodoc[0] != "l"):	out.write("%d\t%f\t%s\n" % (ncount,defocus,flist[i]))
		if(writetodoc[0] == "l"):				res.append(defocus)
	if(len(res) == 0 and  writetodoc == "l" ):				ERROR("No input file is found, check the input directory of file prefix", indir, 1)
	else:
		if writetodoc[0] == "a":
			out.close()
			return res
	if(writetodoc[0] == "l"): 	return res
	if(writetodoc[0] == "w"): 	out.close()

def defocus_get_slow(indir, writetodoc="w", Pixel_size=1, volt=120, Cs=2, wgh=.1, round_off=100, dz_max0=50000, f_l0=30, f_h0=5, prefix="roo", docf="s", skip=";",micdir="", print_screen="p"):
	"""
		Estimate defocus using user provided 1D power spectrum
		mode=1 return the estimated defoci in a list, and writes them down also in a text file
		mode=2 output estimated defocus in a list
		mode=3 output estimated defocus in a text file
		This is a slow version, more accurate than no s version
	"""
	from morphology import defocus_get_Eudis
	import os
	if writetodoc[0]   != "a" and writetodoc[0]   != "l" and writetodoc[0] != "a" : writetodoc   = "a"
	if print_screen[0] != "p" and print_screen[0] != "n": 				print_screen = "n" 
	if os.path.exists(indir) == False: 	ERROR("roodir doesn't exist", "defocus_get_slow",1)
	flist=os.listdir(indir)
	res  = []
	f_l  = f_l0
	f_h  = f_h0
	if f_l <= 1 and f_l > 0:
		 f_l = 1./f_l
		 f_h = 1./f_h
	if f_h > f_l or f_l <= 0 or f_h <= 0: 
		f_h=8.  # angstrom
		f_l=30. # angstrom 
	if round_off <= 0: 	round_off = 100.
	if dz_max0   <= 1: 	dz_max    = 100000.
	dz_max = dz_max0
	if( writetodoc[0] == "w" or writetodoc == "a" ):
		fdefo_nam = "defocus.txt"
		out = open(fdefo_nam, "w")
		out.write("#Coordinates: %s\n")
	ncount = 0	
	for i, v in enumerate(flist):
		(fnam, fext) = os.path.splitext(v)
		if(fnam[0:len(prefix)] == prefix):
			istart   = int(f_l)
			istop    = int(f_h)
			fnam_roo = os.path.join(indir,v)
			Mdis     = 1.e22
			defo     = 0.0
			for nr1 in range(2,7,1):
				for nr2 in range(2,7,1):
					[defocus, dis]     = defocus_get_Eudis(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
					if(print_screen[0]=="p"): print("defocus",defocus,"Euclidean distance",dis,"starting feq",istart,"stop freq",istop,"P R E", nr1,"P R B", nr2)
					if(Mdis > dis):
						defo = defocus
						Mdis = dis
			if(defo > dz_max):
				istart-= 1.
				for nr1 in range(3,5,1):
					for nr2 in range(2,4,1):
						[defocus, dis] = defocus_get_Eudis(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
						if(Mdis>dis):
							defo = defocus
							Mdis = dis
			if(defo >= dz_max): 	ERROR("defo_get_s fails at estimating defocus from ", fnam, 0)
			else:				print("", flist[i], defo) # screen output, give the user a general impression about estimated defoci		
			if writetodoc    == "w" or writetodoc[0] == "a":out.write("%d\t%f\t%s\n" % (ncount, defo, fdefo_nam))
			if writetodoc[0] == "l" : 	res.append(defo)
	if  len(res) == 0 and writetodoc == "l" :  ERROR("No input file, check the input directory", indir, 1)
	else:
		if writetodoc[0] == "a":
			out.close()
			return res
		if writetodoc[0] == "l": 	return res
		if writetodoc[0] == "w":	out.close()

def flcc(t, e):
	"""
		Fast local cross correlation function 
		See Alan Roseman's paper in Ultramicroscopy
	"""
	from utilities import model_blank
	from fundamentals import ccf
	tmp        = EMData()
	mic_avg_sq = EMData()
	mic_sq     = EMData()
	mask       = model_blank(t.get_xsize(), t.get_ysize(), 1)	
	mask       +=1. 
	[mean_t, sigma_t, imin_t, imax_t] = Util.infomask(t,None,False)
	nx         = e.get_xsize()
	ny         = e.get_ysize()		
	n_pixelt   = t.get_xsize()*t.get_ysize()  # get total pixels in template   
	n_pixele   = nx*ny  # get total pixels in mic
	t          = (t-mean_t)/sigma_t # normalize the template such that the average of template is zero.
	t_pad      = Util.pad(t,    nx, ny, 1, {"background":0}, 0, 0, 0)
	m_pad      = Util.pad(mask, nx, ny, 1, {"background":0}, 0, 0, 0) # create a mask (blank, value=1 )file and pad to size of mic   	 	
	tmp        = ccf(e, m_pad)/n_pixele # calculate the local average
	mic_avg_sq = tmp*tmp    # calculate average square
	tmp        = e*e
	mic_sq     = ccf(tmp,m_pad)/n_pixelt 	  # calculate the average of squared mic	       
	tmp        = mic_sq-mic_avg_sq*n_pixelt   #  
	mic_var    = tmp.get_pow(.5)              # Calculate the local variance of the image 
	cc_map     = ccf(e,t_pad)
	cc_map    /= (mic_var*n_pixelt) # Normalize the cross correlation map 
	return cc_map

##-----------------------------img formation parameters related functions---------------------------------
def imf_params_cl1(pw, n=2, iswi=3, Pixel_size=1):
	"""
		Extract image formation parameters using constrained simplex method
		The output is a list of list, which contains the following four elements:
		1. frequencies in 1/Angstrom
		2. fitted curve, either background noise or envelope function
		3. original power spectrum to be fitted
		4. The parameters
		Attention:
		    iswi= 2 using polynomial n rank to fit no-Gaussian envelope function
			iswi =3 using polynomial n rank to fit background
			n = the polynomial rank +1
			The optimization tend to fail when the polynomial rank is higher than 6 
	"""
	feq  = []
	cur  = []
	parm = []
	t    = Util.pw_extract(pw, n, iswi, Pixel_size)
	for i in range(len(pw)):
		cur.append(t[i*2])
		feq.append(t[i*2+1])
	npam = len(t)-2*len(pw)
	for i in range(npam):
		k    = 2*len(pw)+i
		parm.append(t[k])
	return [feq, cur, pw, parm]

def imf_B_factor_get(res_N, x, ctf_params):
	from scipy.optimize import fmin
	nx    = len(res_N)*2
	ctf   = ctf_1d(nx, ctf_params)
	p     = [1,1]
	xopt  = fmin(residuals_B1, p, (res_N,x))
	p     = xopt
	xopt1 = fmin(residuals_B2, p, (res_N,ctf[1][0:nx-1], x))
	print(xopt)
	return xopt

def imf_residuals_B1(p,y,x):
	"""
		Give the initial guess of B-factor
	"""
	from numpy import exp
	C,B = p
	err = 0.0
	for i in range(len(y)):
		err+= abs(y[i] - C*exp(-B*x[i]*x[i]))  # should be 4*B
	return err

def imf_residuals_B2(p,y,ctf,x):
	"""
		fit B-factor in case of considering CTF effect
	""" 
	from numpy import exp
	C,B = p
	err = 0.0
	for i in range(len(y)):
		err+= abs(y[i] - ctf[i]*C*exp(-B*x[i]*x[i]))  # should be 4*B
	return err

def imf_params_get(fstrN, fstrP, ctf_params, pu, nrank, q, lowf=0.01):
	"""
		Extract image formation parameters using optimization method
		Output params: 1. freq; 2.Pn1; 3.B factor.4. C; 5. C*Pu; 6. Pn2
	"""
	params = []
	w      = []
	pw_N   = get_1dpw_list(fstrN)
	pw_P   = get_1dpw_list(fstrP)
	t_N    = imf_params_cl1(pw_N,nrank,3,ctf_params[0])
	t_P    = imf_params_cl1(pw_P,nrank,3,ctf_params[0])
	res_N  = []
	res_P  = []
	for i in range(len(t_N[0])):
		res_N.append(t_N[2][i] - t_N[1][i])
		res_P.append(t_P[2][i] - t_N[1][i])
	params.append(t_N[0]) # freq
	params.append(t_N[1]) # baseline
#	params.append(t_N[1])
	parm1  = imf_B_factor_get(res_N,t_N[0],ctf_params)
	params.append(parm1[1])
	n_lowf = lowf*ctf_params[0]*len(res_P)*2
	n_lowf = int(n_lowf)
	for i in range(len(res_P)):
		if(i <= n_lowf): w.append(0.)
		else:            w.append(1.)
	parm2 = imf_fit_pu(res_P,t_N[0],ctf_params,pu,parm1[0],parm1[1],q,w)
	params.append(parm2[1])
	params.append(parm2[0])
	for i in range(len(res_N)):
		res_N[i] *= q
	params.append(res_N)
	return params

def imf_fit_pu(res_P, x, ctf_params, pu, C, B, q, w):
	from scipy.optimize import fmin
	res   = []
	nx    = len(res_P)*2
	ctf   = ctf_1d(nx, ctf_params)
	for i in range(len(pu)):
		res_P[i] = res_P[i]-q*C*ctf[1][i]*w[i]
		pu[i]   *= ctf[1][i]
	p     = [1]
	xopt  = fmin(residuals_pu,p,(res_P,pu,x))
	res.append(pu)
	res.append(xopt[0])
	return res

def imf_residuals_pu(p,y,pu,x):
	"""
		fit B-factor in case of considering CTF effect
	""" 
	from numpy import exp
	C   = p
	err = 0.0
	for i in range(len(y)):
		err+= abs(y[i] - C*pu[i])
	return err

def residuals_simplex(args, data):
	err      = 0.0
	for i in range(len(data[0])):  err -= (data[0][i] - (args[0] + (args[1]/(data[1][i]/args[2]+1.0)**2)))**2
	return err

def residuals_lsq(p,y,x):
	c1,c2,c3 = p
	err	 = []
	for i in range(len(y)):
		err.append(abs(y[i] - c1-c2/(x[i]+c3)**2))
	return err

def residuals_lsq_peak(p,y,x,c):
	from numpy import exp
	d1,d2,d3 = p
	c1,c2,c3 = c
	err	 = []
	for i in range(len(y)):
		tmp1 = exp(-(x[i] - d2)**2/d3)
		tmp2 = exp(c1)*exp(c2/(x[i] + c3)**2)
		err.append(abs(y[i] - tmp2 - d1*tmp1))
	return err

def residual_1dpw2(list_1dpw2, polynomial_rankB = 2, Pixel_size = 1, cut_off = 0):
	"""
		calculate signal residual from 1D rotationally averaged power spectra 
	"""
	background = []
	freq       = []
	out = Util.pw_extract(list_1dpw2[0:cut_off + 1], polynomial_rankB, 3, Pixel_size )
	for i in range(len(list_1dpw2)):
		j = i*2
		k = i*2+1
		if i <= cut_off:
			res.append(list_1dpw2[i]-background[i])
			freq.append(i/(2*Pixel_size*len(list_1dpw2)))
		else : 
			res.append(0.0)
			freq.append(i/(2*Pixel_size*len(list_1dpw2)))
	return res, freq

def adaptive_mask(vol, nsigma = 1.0, threshold = -9999.0, ndilation = 3, edge_width = 5, mode = "C"):
	"""
		Name
			adaptive_mask - create a mask from a given image.
		Input
			img: input image
			nsigma: value for initial thresholding of the image.
		Output
			mask: The mask will have values one, zero, with cosine smooth transition between two regions.
	"""
	from utilities  import model_circle
	from morphology import binarize, dilation
	nx = vol.get_xsize()
	ny = vol.get_ysize()
	nz = vol.get_zsize()
	mc = model_circle(nx//2, nx, ny, nz) - model_circle(nx//3, nx, ny, nz)
	s1 = Util.infomask(vol, mc, True) # flip true: find statistics under the mask (mask >0.5)
	if threshold <= -9999.0:
		# Use automatic mode
		s1 = [s1[0] + s1[1] * nsigma, s1[0], s1[1], nsigma]
		# new s1[0] is calculated threshold for binarize
	else: 
		# use the user-provided threshold
		if s1[1] != 0.0:
			s1 = [threshold, s1[0], s1[1], (threshold - s1[0])/s1[1]] 
		else:
			s1 = [threshold, s1[0], s1[1], 0.0]
		# new s1[3] is calculated nsigma corresponding to user-provided threshold
	mask = Util.get_biggest_cluster(binarize(vol, s1[0]))
	for i in range(ndilation):   mask = dilation(mask)
	mask = Util.soft_edge(mask, edge_width, mode)
	return mask

def adaptive_mask_scipy(vol, nsigma = 1.0, threshold = -9999.0, ndilation = 3, edge_width = 5, mode = "C", allow_disconnected=False, nerosion = 0, do_approx=False):
	"""
		Name
			adaptive_mask - create a mask from a given image.
		Input
			img: input image
			nsigma: value for initial thresholding of the image.
		Output
			mask: The mask will have values one, zero, with cosine smooth transition between two regions.
	"""
	from utilities  import model_circle
	from morphology import binarize, dilation
	nx = vol.get_xsize()
	ny = vol.get_ysize()
	nz = vol.get_zsize()
	mc = model_circle(nx//2, nx, ny, nz) - model_circle(nx//3, nx, ny, nz)
	s1 = Util.infomask(vol, mc, True) # flip true: find statistics under the mask (mask >0.5)
	if threshold <= -9999.0:
		# Use automatic mode
		bin_threshold = s1[0] + s1[1] * nsigma
		#s1 = [s1[0] + s1[1] * nsigma, s1[0], s1[1], nsigma]
		# new s1[0] is calculated threshold for binarize
	else: 
		# use the user-provided threshold
		bin_threshold = threshold
		#if s1[1] != 0.0:
		#	s1 = [threshold, s1[0], s1[1], (threshold - s1[0])/s1[1]] 
		#else:
		#	s1 = [threshold, s1[0], s1[1], 0.0]
		# new s1[3] is calculated nsigma corresponding to user-provided threshold

	mask = binarize(vol, bin_threshold)
	if not allow_disconnected:
		mask = Util.get_biggest_cluster(mask)
	for i in range(ndilation):
		mask = dilation(mask)
	for i in range(nerosion):
		mask = erosion(mask)
	if edge_width > 0:
		mask = soft_edge(mask, edge_width, mode, do_approx)
	return mask

'''
def adaptive_mask2D(img, nsigma = 1.0, ndilation = 3, kernel_size = 11, gauss_standard_dev =9):
	"""
		Name
			adaptive_mask - create a mask from a given image.
		Input
			img: input image
			nsigma: value for initial thresholding of the image.
		Output
			mask: The mask will have values one, zero, with Gaussian smooth transition between two regions.
	"""
	from utilities  import gauss_edge, model_circle
	from morphology import binarize, dilation
	nx = img.get_xsize()
	ny = img.get_ysize()
	mc = model_circle(nx//2, nx, ny) - model_circle(nx//3, nx, ny)
	s1 = Util.infomask(img, mc, True)
	mask = Util.get_biggest_cluster(binarize(img, s1[0]+s1[1]*nsigma))
	for i in range(ndilation):   mask = dilation(mask)
	#mask = gauss_edge(mask, kernel_size, gauss_standard_dev)
	return mask
'''

def cosinemask(im, radius = -1, cosine_width = 5, bckg = None, s=999999.0):
	"""
		Apply mask with a cosine fall-off setting values outside of radius_cosine_width to the average computed outside.
		The fall-off begins from pixel at a distance radius from the center,
		i.e., mask(radius) = 1 and mask(radius+cosine_width)=0.
		if s=999999.0 using average else program takes in user-provided s  
	"""
	return  Util.cosinemask(im, radius, cosine_width, bckg, s)
'''
	from utilities import model_blank
	from math import cos, sqrt, pi
	nx = im.get_xsize()
	ny = im.get_ysize()
	nz = im.get_zsize()
	if(radius < 0):
		if(ny == 1):    radius = nx//2 - cosine_width
		elif(nz == 1):  radius = min(nx,ny)//2 - cosine_width
		else:           radius = min(nx,ny,nz)//2 - cosine_width
	radius_p = radius + cosine_width
	om = im.copy()
	cz = nz//2
	cy = ny//2
	cx = nx//2
	if bckg:
		for z in xrange(nz):
			tz = (z-cz)**2
			for y in xrange(ny):
				ty = tz + (y-cy)**2
				for x in xrange(nx):
					r = sqrt(ty + (x-cx)**2)
					if(r > radius_p):
						om.set_value_at_fast(x,y,z, bckg.get_value_at(x,y,z))
					elif(r>=radius):
						temp = (0.5 + 0.5 * cos(pi*(radius_p - r)/cosine_width ))
						om.set_value_at_fast(x,y,z, om.get_value_at(x,y,z) + temp*(bckg.get_value_at(x,y,z)-om.get_value_at(x,y,z)))
	else:
		u = 0.0
		s = 0.0
		for z in xrange(nz):
			tz = (z-cz)**2
			for y in xrange(ny):
				ty = tz + (y-cy)**2
				for x in xrange(nx):
					r = sqrt(ty + (x-cx)**2)
					if(r > radius_p):
						u += 1.0
						s += om.get_value_at(x,y,z)
					elif(r>=radius):
						temp = (0.5 + 0.5 * cos(pi*(radius_p - r)/cosine_width ))
						u += temp
						s += om.get_value_at(x,y,z)*temp
		s /= u
		for z in xrange(nz):
			tz = (z-cz)**2
			for y in xrange(ny):
				ty = tz + (y-cy)**2
				for x in xrange(nx):
					r = sqrt(ty + (x-cx)**2)
					if(r > radius_p):
						om.set_value_at_fast(x,y,z, s)
					elif(r>=radius):
						temp = (0.5 + 0.5 * cos(pi*(radius_p - r)/cosine_width ))
						om.set_value_at_fast(x,y,z, om.get_value_at(x,y,z) + temp*(s-om.get_value_at(x,y,z)))
						#om.set_value_at_fast(x,y,z, om.get_value_at(x,y,z)*(0.5 + 0.5 * cos(pi*(radius_p - r)/cosine_width )))
	return om
	"""
'''


def get_shrink_3dmask(nxinit, mask_file_name):
	from utilities import get_im
	from fundamentals import resample
	from morphology   import binarize
	mask3d = get_im(mask_file_name)
	nx2 = nxinit
	nx1 = mask3d.get_xsize()
	if nx1 == nx2:
		return mask3d
	else:
		shrinkage = float(nx2)/nx1
		mask3d    = binarize(resample(mask3d,shrinkage),0.5) #added 0.5 here to fix binarization problem
		return mask3d


def get_biggest_cluster(mg):
	"""
	  Input: binary image
	  Output: image that contains the largest connected subset in the input image
	  This code was written by Wei in C and put in util_sparx.cpp
	"""
	from utilities import model_blank
	from morphology import collapse

	nx = mg.get_xsize()
	ny = mg.get_ysize()
	nz = mg.get_zsize()

	lg = mg.copy()
	s = Util.infomask(lg, None, True)
	nnc = int(s[0]*nx*ny*nz)

	cls = model_blank(nx,ny,nz)

	l = []
	grp = 0

	# endless loop
	while(grp < 1):
		grp -= 1

		fif = True
		for k in range(nz):
			if(fif):
				for j in range(ny):
					if(fif):
						for i in range(nx):
							if(fif and lg[i,j,k]):
								lg[i,j,k]=0
								l.append([i,j,k])
								fif = False

		if(fif):
			#  All points checked, extract the largest group
			grp = -1-grp
			# check if any groups found
			if(grp == 0):  return cls
			cls *= -1
			# if one group, simply return it
			if(grp == 1):  return cls
			st = 0.0
			for ig in range(1,grp):
				s = Util.infomask(collapse(cls, ig-0.5, ig+0.5), None, True)
				if(s[0] > st):
					st = s[0]
					tig = ig
			return  collapse(cls, tig-0.5, tig+0.5)

		while(len(l) > 0):
			cr = l[0]
			del l[0]
			cls[cr[0], cr[1], cr[2]] = grp
			lg[cr[0], cr[1], cr[2]] = 0
			for k in range(-1,2,1):
				kq = cr[2]+k
				if(kq>-1 and kq<nz):
					for j in range(-1,2,1):
						jq = cr[1]+j
						if(jq>-1 and jq<ny):
							for i in range(-1,2,1):
								iq = cr[0]+i
								if(iq>-1 and iq<nx):
									if(lg[iq,jq,kq]):
										lg[iq,jq,kq]=0
										l.append([iq,jq,kq])

def adaptive_mask_mass(vol, mass=2000, Pixel_size=3.6):
	from utilities  import gauss_edge, model_blank
	from morphology import binarize, threshold, dilation
	from filter     import filt_gaussl
	nx = vol.get_xsize()
	a = filt_gaussl(vol, 0.15, True)
	TH = a.find_3d_threshold(mass, Pixel_size)
	a = binarize(a,TH)
	d = a.delete_disconnected_regions(0,0,0)

	d = dilation(d, model_blank(3,3,3,1.0), "BINARY")
	#d = filt_dilation(d, model_blank(3,3,3,1.0), "BINARY")
	d = gauss_edge(d)
	return d
	#Util.mul_img(vol, d)
	#return threshold(vol, 0.0)

"""
def refine_with_mask(vol):
	from filter     import filt_dilation
	from utilities  import model_circle, model_gauss, drop_image
	from morphology import collapse
	# does not seem to be working all that well
	nx = vol.get_xsize()
	outer_radius = nx/2-2
	inner_radius = nx/2-4
	outer_sphere = model_circle(outer_radius, nx, nx, nx)
	inner_sphere = model_circle(inner_radius, nx, nx, nx)
	shell = outer_sphere - inner_sphere
	avg,sigma,amin,amax = Util.infomask(vol, shell, False)
	print  avg,sigma,amin,amax
	vol -= avg
	mask = collapse(vol, -1.5*sigma, 1.5*sigma)
	#from utilities import drop_image
	#drop_image(mask,"m1.spi","s")
	mask -= 1.0
	mask *= -1.0

	att = model_circle(1,3,3,3)
	mask = filt_dilation(filt_dilation(filt_dilation(mask, att, "BINARY"), att, "BINARY"), att, "BINARY")

	gauss = model_gauss(2.0, 9, 9, 9)
	avg,sigma,amin,amax = Util.infomask(gauss, None, False)
	gauss /= (avg*9*9*9)
	mask = rsconvolution(mask, gauss)
	#drop_image(mask,"m2.spi","s")
	vol *= mask
	return vol
"""


def compute_bfactor(pws, freq_min, freq_max, pixel_size = 1.0):
	"""
		Estimate B-factor from power spectrum
		pws          : 1D rotational average of power spectrum, length should be half of the image size
		idx_freq_min : the index of the minimum frequency of fitting range
		idx_freq_max : the index of the maximum frequency of fitting range
		pixel_size   :  in A
	"""
	from math import log, sqrt
	from statistics import linreg
	nr = len(pws)
	"""
	if (idx_freq_min < 0):
		ERROR("compute_bfactor", "Invalid value of idx_freq_min. Setting to 0", 0)
		idx_freq_min = 0
	if (idx_freq_max >= nr):
		pERROR("compute_bfactor", "Invalid value of idx_freq_max. Setting to %d" % (nr - 1), 0)
		idx_freq_max = (nr - 1)
	"""

	pws_log = [0.0]*nr
	x = [0.0]*nr
	q = min(pws)
	for i in range(1,nr):
		pws_log[i] = log(pws[i])  #/q)
		x[i] = (float(i)/(2*nr)/pixel_size)**2
	idx_freq_min = 1
	for i in range(1,nr):
		if(x[i] > freq_min**2):
			idx_freq_min = i
			break

	idx_freq_max = 1
	for i in range(1,nr):
		idx_freq_max = i
		if(x[i] > freq_max**2):
			break

	# Linear regression will crash if min & max frequencies are only one apart
	if idx_freq_max-idx_freq_min <= 1: 
		ERROR("B_start is too high a resolution! Decrease it (under Advanced) and re-run the program! ", "compute_bfactor")
		
	B, s = linreg(x[idx_freq_min:idx_freq_max], pws_log[idx_freq_min:idx_freq_max])
	#print  B,s

	ff = [0.0]*nr
	from math import exp
	for i in range(nr):  ff[i] = B*x[i] + s

	return -B, [x,ff,pws_log], idx_freq_min, idx_freq_max


################
#
#  CTER code (new version since 2016/03/16)
#
################
# 
# NOTE: 2016/03/16 Toshio Moriya
# In this version, the IO-related interface is simplified for sxcter.py and sxgui.py
# Since cter() was used in not only sxcter.py but also e2boxer.py and sxhelixboxer.py,
# This new version is added to avoid breaking e2boxer.py and sxhelixboxer.py
#  
# NOTE: 2016/03/16 Toshio Moriya
# To get a single  file name from a GUI application, 
# there must be a better way than using guimic...
# 
# NOTE: 2016/11/16 Toshio Moriya
# Now, this function assume the MPI setup and clean up is done by caller, such as mpi_init, and mpi_finalize
# 
# NOTE: 07/11/2017  PAP
#       This is "exact" copy of mrk version with a switch to amplitudes (square root of PW)
# 
def cter_mrk(input_image_path, output_directory, selection_list = None, wn = 512, pixel_size = -1.0, \
			Cs = 2.0, voltage = 300.0, wgh = 10.0, f_start = -1.0, f_stop = -1.0, \
			kboot = 16, overlap_x = 50, overlap_y = 50, edge_x = 0, edge_y = 0, \
			check_consistency = False, stack_mode = False, debug_mode = False, \
			program_name = "cter_mrk() in morphology.py", \
			RUNNING_UNDER_MPI = False, main_mpi_proc = 0, my_mpi_proc_id = 0, n_mpi_procs = 1):
	"""
	Arguments
		input_image_path  :  file name pattern for Micrographs Modes (e.g. 'Micrographs/mic*.mrc') or particle stack file path for Stack Mode (e.g. 'bdb:stack'; must be stack_mode = True).
		output_directory  : output directory
	"""
	from   EMAN2 import periodogram
	from   EMAN2db import db_check_dict, db_parse_path
	from   applications import MPI_start_end
	from   utilities import read_text_file, write_text_file, get_im, model_blank, model_circle, amoeba, generate_ctf
	from   utilities import if_error_then_all_processes_exit_program
	from   utilities import wrap_mpi_bcast
	from   sys import exit
	import numpy as np
	import os
	import glob
	from   fundamentals import tilemic, rot_avg_table, resample
	from   morphology   import threshold, bracket_def, bracket, goldsearch_astigmatism
	from   morphology   import defocus_baseline_fit, simpw1d, movingaverage, localvariance, defocusgett
	from   morphology   import defocus_guessn, defocusget_from_crf, make_real
	from   morphology   import fastigmatism, fastigmatism1, fastigmatism2, fastigmatism3, simctf, simctf2, simctf2out, fupw,ctf2_rimg
	from   alignment    import Numrinit, ringwe
	from   statistics   import table_stat
	from   pixel_error  import angle_ave
	from   inspect      import currentframe, getframeinfo
	from   global_def   import ERROR
	import global_def
	from   time import time
	from   mpi import MPI_COMM_WORLD, mpi_barrier

	# ====================================================================================
	# Prepare processing
	# ====================================================================================

	
	# ------------------------------------------------------------------------------------
	# Find the CTER Running Mode before checking error conditions
	# ------------------------------------------------------------------------------------
	i_enum = -1; idx_cter_mode_invalid       = i_enum; 
	i_enum += 1; idx_cter_mode_all_mics      = i_enum  # All Micrographs Mode - Process all s in a directory
	i_enum += 1; idx_cter_mode_selected_mics = i_enum  # Selected Micrographs Mode - Process all s in a selection list file
	i_enum += 1; idx_cter_mode_single_mic    = i_enum  # Single Micrograph Mode - Process a single 
	i_enum += 1; idx_cter_mode_stack         = i_enum  # Stack Mode - Process a stack (Advanced Option)
	i_enum += 1; idx_cter_mode_counts        = i_enum
	
	cter_mode_idx = idx_cter_mode_invalid
	cter_mode_name = None
	if stack_mode == False:
		# One of three Micrograph Modes
		# For any of Micrograph Modes, input image file name must be a file path pattern containing wild card "*" 
		if selection_list == None:
			# User did not use selection list option 
			# -> All Micrographs Mode
			cter_mode_idx = idx_cter_mode_all_mics
			cter_mode_name = "All Micrographs Mode"
		else:
			if os.path.splitext(selection_list)[1] == ".txt":
				# User specified a selection list text file path containing".txt" extension through selection list option
				# -> Selected Micrographs Mode
				cter_mode_idx = idx_cter_mode_selected_mics
				cter_mode_name = "Selected Micrographs Mode"
			else: 
				# User specified an image file path (a non-text file path) through selection list option
				# -> Single Micrograph Mode
				cter_mode_idx = idx_cter_mode_single_mic
				cter_mode_name = "Single Micrograph Mode"
	else: 
		# (Particle) Stack Mode
		cter_mode_idx = idx_cter_mode_stack
		cter_mode_name = "Stack Mode"
	
	if my_mpi_proc_id == main_mpi_proc:
		print(" ")
		print(("----- Running with %s -----" % (cter_mode_name)))
	
	# ------------------------------------------------------------------------------------
	# Check mode-dependent error conditions of input arguments and options if abort is necessary. All nodes do this checking
	# ------------------------------------------------------------------------------------
	error_message_list = [] # List of error messages. If no error is found, the length should be zero
	if not stack_mode:
		
		# Check error conditions applicable to any of Micrograph Mode 
		if input_image_path.find("*") == -1:
			error_message_list.append("Input image file path (%s) for %s must be a  path pattern containing wild card (*). Please check input_image_path argument." % (input_image_path, cter_mode_name))
		
		if input_image_path[:len("bdb:")].lower() == "bdb:":
			error_message_list.append("BDB file can not be selected as input image file path (%s) for %s. Please check input_image_path argument and convert the image format." % (input_image_path, cter_mode_name))
		
		# Check error conditions applicable to Selected Micrographs Mode 
		if cter_mode_idx == idx_cter_mode_selected_mics:
			if not os.path.exists(selection_list): 
				error_message_list.append("Selection list text file specified by selection_list option (%s) for %s does not exists. Please check selection_list option." % (selection_list, cter_mode_name))
		
		if cter_mode_idx == idx_cter_mode_single_mic:
			if not os.path.exists(os.path.join(os.path.dirname(input_image_path), os.path.basename(selection_list))): 
				error_message_list.append("Micrograph specified by selection_list option (%s) for %s does not exist. Please check selection_list option." % (selection_list, cter_mode_name))
			# 
			if RUNNING_UNDER_MPI and n_mpi_procs != 1:
				error_message_list.append("%s supports only a single processor version. Please change MPI settings." % (cter_mode_name))
		
	else: 
		# Check error conditions
		if input_image_path.find("*") != -1:
			error_message_list.append("Stack file path specified by input_image_path (%s) for %s should not contain wild card (*). Please check input_image_path argument." % (input_image_path, cter_mode_name))
		
		is_not_found_input_image_file = False
		if input_image_path[:len("bdb:")].lower() == "bdb:":
			if not db_check_dict(input_image_path): 
				is_not_found_input_image_file = True
		else:
			if not os.path.exists(input_image_path): 
				is_not_found_input_image_file = True
		if is_not_found_input_image_file:
			error_message_list.append("Stack file specified by input_image_path (%s) for %s does not exist. Please check input_image_path argument." % (input_image_path, cter_mode_name))
		
		if RUNNING_UNDER_MPI and n_mpi_procs != 1:
			error_message_list.append("%s supports only a single processor version. Please change MPI settings." % (cter_mode_name))
	
	# --------------------------------------------------------------------------------
	# check output-related error conditions (mode-independent). All nodes do this checking
	# --------------------------------------------------------------------------------
	if os.path.exists(output_directory):
		error_message_list.append("Output directory (%s) exists already. Please check output_directory argument." % (output_directory))
	
	# --------------------------------------------------------------------------------
	# Check error conditions of options (mode-independent). All nodes do this checking
	# --------------------------------------------------------------------------------
	if pixel_size <= 0.0:
		error_message_list.append("Pixel size (%f) must not be negative. Please set a pasitive value larger than 0.0 to pixel_size option." % (pixel_size))

	if wn <= 0.0:
		error_message_list.append("CTF window size (%d) must not be negative. Please set a valid value larger than 0 to wn option." % (wn))
	
	# --------------------------------------------------------------------------------
	# Print all error messages and abort the process if necessary.
	# --------------------------------------------------------------------------------
	error_status = None
	if len(error_message_list) > 0:
		# Detected error! Print all error messages
		if my_mpi_proc_id == main_mpi_proc:
			print(" ")
			for error_message in error_message_list:  
				print(("ERROR!!! %s" % (error_message)))
		error_status = ("Detected %d error(s) related to arguments and options. Run %s -h for help. Exiting..." % (len(error_message_list), program_name), getframeinfo(currentframe()))
	if_error_then_all_processes_exit_program(error_status)
	if RUNNING_UNDER_MPI:
		# Wait for all mpi processes to check error conditions, especially existence of output directory
		# Without this barrier, main mpi process can create output directory before some child mpi process check this error.
		mpi_barrier(MPI_COMM_WORLD)

	del error_message_list # Don't need this anymore
	
	# ------------------------------------------------------------------------------------
	# Check warning conditions of options
	# ------------------------------------------------------------------------------------
	if my_mpi_proc_id == main_mpi_proc:
		if stack_mode:
			if selection_list != None:
				print(" ")
				print(("WARNING!!! --selection_list option will be ignored in %s." % (cter_mode_name)))
			if wn != 512:
				print(" ")
				print(("WARNING!!! --wn option will be ignored in %s." % (cter_mode_name)))
			if overlap_x != 50:
				print(" ")
				print(("WARNING!!! --overlap_x option will be ignored in %s." % (cter_mode_name)))
			if overlap_y != 50:
				print(" ")
				print(("WARNING!!! --overlap_y option will be ignored in %s." % (cter_mode_name)))
			if edge_x != 0:
				print(" ")
				print(("WARNING!!! --edge_x option will be ignored in %s." % (cter_mode_name)))
			if edge_y != 0:
				print(" ")
				print(("WARNING!!! --edge_y option will be ignored in %s." % (cter_mode_name)))
			if check_consistency:
				print(" ")
				print(("WARNING!!! --check_consistency option will be ignored in %s." % (cter_mode_name)))

	# ====================================================================================
	# Create the input file path list and also check input-related error conditions if abort is necessary.
	# ====================================================================================
	input_file_path_list = []
	if not stack_mode:
		# --------------------------------------------------------------------------------
		# Prepare the variables for all sections in  mode case 
		# --------------------------------------------------------------------------------
		# Micrograph basename pattern (directory path is removed from  path pattern)
		mic_pattern = input_image_path
		mic_basename_pattern = os.path.basename(mic_pattern)
	
		# Global entry dictionary (all possible entries from all lists) for all mic id substring
		global_entry_dict = {} # mic id substring is the key
		subkey_input_mic_path = "Input Micrograph Path"
		subkey_selected_mic_basename = "Selected Micrograph Basename"
	
		# List keeps only id substrings of s whose all necessary information are available
		valid_mic_id_substr_list = [] 
		
		# --------------------------------------------------------------------------------
		# Obtain the list of  id sustrings using a single CPU (i.e. main mpi process)
		# --------------------------------------------------------------------------------
		# NOTE: Toshio Moriya 2016/11/15
		# The below is not a real while.  
		# It gives if-statements an opportunity to use break when errors need to be reported
		# However, more elegant way is to use 'raise' statement of exception mechanism...
		# 
		error_status = None
		while my_mpi_proc_id == main_mpi_proc:
			
			# --------------------------------------------------------------------------------
			# Prepare variables for this section
			# --------------------------------------------------------------------------------
			# Prefix and suffix of  basename pattern 
			# to find the head/tail indices of  id substring
			mic_basename_tokens = mic_basename_pattern.split('*')
			# Find head index of  id substring
			mic_id_substr_head_idx = len(mic_basename_tokens[0])
		
			# --------------------------------------------------------------------------------
			# Register  id substrings found in the input directory (specified by  path pattern)
			# to the global entry dictionary
			# --------------------------------------------------------------------------------
			# Generate the list of  paths in the input directory
			print(" ")
			print("Checking the input directory...")
			input_mic_path_list = glob.glob(mic_pattern)
			# Check error condition of input  file path list
			print(("Found %d micrographs in %s." % (len(input_mic_path_list), os.path.dirname(mic_pattern))))
			if len(input_mic_path_list) == 0:
				# The result shouldn't be empty if the specified  file name pattern is invalid
				error_status = ("There are no micrographs whose paths match with the specified file path pattern (%s) for %s. Please check input_image_path. Run %s -h for help." % (mic_pattern, cter_mode_name, program_name), getframeinfo(currentframe()))
				break
		
			# Register  id substrings to the global entry dictionary
			for input_mic_path in input_mic_path_list:
				# Find tail index of  id substring and extract the substring from the  name
				input_mic_basename = os.path.basename(input_mic_path)
				mic_id_substr_tail_idx = input_mic_basename.index(mic_basename_tokens[1])
				mic_id_substr = input_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
				if not mic_id_substr in global_entry_dict:
					# print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from input_mic_path_list " % (mic_id_substr))
					global_entry_dict[mic_id_substr] = {}
				global_entry_dict[mic_id_substr][subkey_input_mic_path] = input_mic_path
		
			# --------------------------------------------------------------------------------
			# Register  id substrings found in the selection list
			# to the global entry dictionary
			# --------------------------------------------------------------------------------
			# Generate the list of selected  paths in the selection file
			selected_mic_path_list = []
			# Generate  lists according to the execution mode
			if cter_mode_idx == idx_cter_mode_all_mics:
				# Treat all s in the input directory as selected ones
				selected_mic_path_list = input_mic_path_list
			else:
				if os.path.splitext(selection_list)[1] == ".txt":
					print(" ")
					print("Checking the selection list...")
					selected_mic_path_list = read_text_file(selection_list)
				
					# Check error condition of  entry lists
					print(("Found %d microgarph entries in %s." % (len(selected_mic_path_list), selection_list)))
					if len(selected_mic_path_list) == 0:
						error_status = ("The provided  list file (%s) for %s mode contains no entries. Please check selection_list option and make sure the file contains a  list. Run %s -h for help." % (selection_list, cter_mode_name, program_name), getframeinfo(currentframe()))
						break
				else:
					print(" ")
					print(("Processing a single micorgprah: %s..." % (selection_list)))
					selected_mic_path_list = [selection_list]
			
				selected_mic_directory = os.path.dirname(selected_mic_path_list[0])
				if selected_mic_directory != "":
					print(("    NOTE: Program disregards the directory paths in the selection list (%s)." % (selected_mic_directory)))
		
			# Register  id substrings to the global entry dictionary
			for selected_mic_path in selected_mic_path_list:
				# Find tail index of  id substring and extract the substring from the  name
				selected_mic_basename = os.path.basename(selected_mic_path)
				mic_id_substr_tail_idx = selected_mic_basename.index(mic_basename_tokens[1])
				mic_id_substr = selected_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
				if not mic_id_substr in global_entry_dict:
					# print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from selected_mic_path_list " % (mic_id_substr))
					global_entry_dict[mic_id_substr] = {}
				global_entry_dict[mic_id_substr][subkey_selected_mic_basename] = selected_mic_basename
			
			# --------------------------------------------------------------------------------
			# Clean up variables related to registration to the global entry dictionary
			# --------------------------------------------------------------------------------
			del mic_basename_tokens
			del mic_id_substr_head_idx
		
			# --------------------------------------------------------------------------------
			# Create the list containing only valid  id substrings
			# --------------------------------------------------------------------------------
			# Prepare lists to keep track of invalid (rejected) s 
			no_input_mic_id_substr_list = []
		
			print(" ")
			print("Checking the input datasets consistency...")
		
			# Loop over substring id list
			for mic_id_substr in global_entry_dict:
				mic_id_entry = global_entry_dict[mic_id_substr]
			
				warinnig_messages = []
				# selected  basename must have been registed always .
				if subkey_selected_mic_basename in mic_id_entry: 
					# Check if associated input  exists
					if not subkey_input_mic_path in mic_id_entry:
						input_mic_path = mic_pattern.replace("*", mic_id_substr)
						warinnig_messages.append("    associated input  %s." % (input_mic_path))
						no_input_mic_id_substr_list.append(mic_id_substr)
				
					if len(warinnig_messages) > 0:
						print(("WARNING!!! Micrograph ID %s does not have:" % (mic_id_substr)))
						for warinnig_message in warinnig_messages:
							print(warinnig_message)
						print("    Ignores this as an invalid entry.")
					else:
						# print("MRK_DEBUG: adding mic_id_substr := ", mic_id_substr)
						valid_mic_id_substr_list.append(mic_id_substr)
				# 	# This entry is not in the selection list. Do nothing
			
			# Check the input dataset consistency and save the result to a text file, if necessary.
			if check_consistency:
				# Create output directory
				os.mkdir(output_directory)
			
				# Open the consistency check file
				inconsist_mic_list_path = os.path.join(output_directory,"inconsist_mic_id_file.txt")
				print(" ")
				print(("Generating the input datasets consistency report in %s..." % (inconsist_mic_list_path)))
				inconsist_mic_list_file = open(inconsist_mic_list_path, "w")
				inconsist_mic_list_file.write("# The information about inconsistent  IDs\n")
				# Loop over substring id list
				for mic_id_substr in global_entry_dict:
					mic_id_entry = global_entry_dict[mic_id_substr]
				
					consistency_messages = []
					# Check if associated input  path exists
					if not subkey_input_mic_path in mic_id_entry:
						input_mic_path = mic_pattern.replace("*", mic_id_substr)
						consistency_messages.append("    associated input  %s." % (input_mic_path))
				
					# Check if associated selected  basename exists
					if not subkey_selected_mic_basename in mic_id_entry:
						input_mic_path = mic_pattern.replace("*", mic_id_substr)
						consistency_messages.append("    associated selected  %s." % (input_mic_path))
				
					if len(consistency_messages) > 0:
						inconsist_mic_list_file.write("Micrograph ID %s does not have:\n" % (mic_id_substr))
						for consistency_message in consistency_messages:
							inconsist_mic_list_file.write(consistency_message)
							inconsist_mic_list_file.write("\n")
			
				# Close the consistency check file, if necessary
				inconsist_mic_list_file.flush()
				inconsist_mic_list_file.close()
			
			# Since mic_id_substr is once stored as the key of global_entry_dict and extracted with the key order
			# we need sort the valid_mic_id_substr_list here
			if debug_mode: print(("BEFORE SORT: valid_mic_id_substr_list := ", valid_mic_id_substr_list))
			valid_mic_id_substr_list.sort(key=str.lower) # Sort list of  IDs using case insensitive string comparison
			if debug_mode: print(("AFTER SORT: valid_mic_id_substr_list := ", valid_mic_id_substr_list))
			
			# --------------------------------------------------------------------------------
			# Print out the summary of input consistency
			# --------------------------------------------------------------------------------
			print(" ")
			print("Summary of dataset consistency check...")
			print(("  Detected  IDs               : %6d" % (len(global_entry_dict))))
			print(("  Entries in input directory  : %6d" % (len(input_mic_path_list))))
			print(("  Entries in selection list   : %6d" % (len(selected_mic_path_list))))
			print(("  Rejected by no input        : %6d" % (len(no_input_mic_id_substr_list))))
			print(("  Valid Entries               : %6d" % (len(valid_mic_id_substr_list))))
			
			# --------------------------------------------------------------------------------
			# Check MPI error condition
			# --------------------------------------------------------------------------------
			if len(valid_mic_id_substr_list) < n_mpi_procs:
				error_status = ("Number of MPI processes (%d) supplied by --np in mpirun cannot be greater than %d (number of valid s that satisfy all criteria to be processed). Run %s -h for help." % (n_mpi_procs, len(valid_mic_id_substr_list, program_name)), getframeinfo(currentframe()))
				break
			
			# --------------------------------------------------------------------------------
			# Create input file path list
			# --------------------------------------------------------------------------------
			for mic_id_substr in valid_mic_id_substr_list:
				mic_path = global_entry_dict[mic_id_substr][subkey_input_mic_path]
				input_file_path_list.append(mic_path)
			
			# --------------------------------------------------------------------------------
			# Clean up variables related to tracking of invalid (rejected) s 
			# --------------------------------------------------------------------------------
			del input_mic_path_list
			del selected_mic_path_list
			del no_input_mic_id_substr_list
			
			break
		
		
		# --------------------------------------------------------------------------------
		# Clean up the variables for all sections in  mode case 
		# --------------------------------------------------------------------------------
		del mic_pattern
		del mic_basename_pattern
		del global_entry_dict
		del subkey_input_mic_path
		del subkey_selected_mic_basename
		del valid_mic_id_substr_list
		
		# --------------------------------------------------------------------------------
		# Print all error messages and abort the process if necessary.
		# --------------------------------------------------------------------------------
		# NOTE: Toshio Moriya 2016/11/15
		# The following function takes care of the case when an if-statement uses break for occurence of an error.
		# However, more elegant way is to use 'exception' statement of exception mechanism...
		# 
		if_error_then_all_processes_exit_program(error_status)
		
	else:
		input_file_path_list.append(input_image_path)
	
	if RUNNING_UNDER_MPI:
		# Wait for main mpi process to create the input file path list
		mpi_barrier(MPI_COMM_WORLD)
		
		# All mpi processes should know input file path list
		input_file_path_list = wrap_mpi_bcast(input_file_path_list, main_mpi_proc)
	
	# ====================================================================================
	# Prepare input file path(s)
	# ====================================================================================
	# 
	# NOTE: 2016/03/17 Toshio Moriya
	# From here on, stack (and namics) will be used to distinguish stack mode and  mode.
	# However, a single input_file_path_list should be sufficient since we already know the mode.
	# Let's consider this refactoring in the future.
	# 
	stack = None # (particle) stack file name: if it is not None, cter runs with stack mode. Otherwise, runs with  mode
	namics = []  #  file name list
	if not stack_mode:
		namics = input_file_path_list
		if debug_mode: print(("BEFORE SORT: namics := ", namics))
		namics.sort(key=str.lower) # Sort list of s using case insensitive string comparison
		if debug_mode: print(("AFTER SORT: namics := ", namics))
	else:
		stack = input_file_path_list[0]
	
	del input_file_path_list # Don't need this anymore
	
	# Make output directory
	outpwrot = "%s/pwrot" % (output_directory)
	if stack == None: 
		outmicthumb = "%s/micthumb" % (output_directory)
	if debug_mode:  
		outravg = "%s/ravg" % (output_directory)
	if my_mpi_proc_id == main_mpi_proc:
		# Make output directory
		if not os.path.exists(output_directory):
			os.mkdir(output_directory)
		os.mkdir(outpwrot)
		if stack == None: 
			os.mkdir(outmicthumb)
		if debug_mode:
			os.mkdir(outravg)
	
	if RUNNING_UNDER_MPI:
		# Make all mpi processes wait for main mpi process to create output directory
		mpi_barrier(MPI_COMM_WORLD)
	
	# Set up loop variables depending on the cter mode
	if stack == None:
		if RUNNING_UNDER_MPI:
			set_start, set_end = MPI_start_end(len(namics), n_mpi_procs, my_mpi_proc_id)
		else:
			set_start = 0
			set_end = len(namics)
	else: 
		pw2 = []
		data = EMData.read_images(stack)
		nima = len(data)
		for i in range(nima):
			pw2.append(periodogram(data[i]))
		wn = pw2[0].get_xsize()
		set_start = 0
		set_end = 1
	
	# Set up progress message
	if my_mpi_proc_id == main_mpi_proc:
		print(" ")
		print("Estimating CTF parameters...")
		if stack == None:
			print("  Micrographs processed by main process (including percent of progress):")
			progress_percent_step = (set_end - set_start)/100.0 # the number of micrograms for main mpi processer divided by 100
	
	totresi = []
	missing_img_names = []
	rejected_img_names = []
	for ifi in range(set_start, set_end):
		img_type = ""
		img_name = ""
		img_basename_root = ""
		
		if stack == None:
			img_type = "Micrograph"
			img_name = namics[ifi]
			
			if my_mpi_proc_id == main_mpi_proc:
				print(("    Processing %s ---> %6.2f%%" % (img_name, (ifi - set_start) / progress_percent_step)))
			
			if not os.path.exists(img_name):
				missing_img_names.append(img_name)
				print("    %s %s: Can not find this file. Skipping the estimation and CTF parameters are not stored..." % (img_type, img_name))
				continue

			mic = get_im(img_name)
			try:
				pw2 = tilemic(mic, win_size = wn, overlp_x = overlap_x, overlp_y = overlap_y, edge_x = edge_x, edge_y = edge_y)
			except:
				print("MRK_DEBUG: tilemic() in cter_mrk() raised an exception. The micrographs {} might have a problem. Please check it and remove it if necessary.".format(img_name))
				raise
			del mic
		else:
			img_type = "Stack"
			img_name = stack
			
			numFM = EMUtil.get_image_count(img_name)
			pw2 = [None]*numFM
			for i in range(numFM):
				pw2.append(periodogram(get_im(img_name,i)))

		for i in range(len(pw2)):
			pw2[i] = square_root(pw2[i])
		if debug_mode: print("    %s %s: Process %04d started the processing. Detected %d image(s) in this %s file." % (img_type, img_name, ifi, numFM, img_type.lower()))

		if db_check_dict(img_name) == False:
			img_basename_root = os.path.splitext(os.path.basename(img_name))[0]
		else:
			path, dictname, keys = db_parse_path(img_name)
			img_basename_root = dictname

		
		nimi = len(pw2)
		adefocus = [0.0] * kboot
		aamplitu = [0.0] * kboot
		aangle   = [0.0] * kboot
		
		allroo = []
		for imi in range(nimi):
			allroo.append(rot_avg_table(pw2[imi]))
		lenroo = len(allroo[0])
		#print time(),nimi
		
		for nboot in range(kboot):
			#at = time()
			if(nboot == 0): boot = list(range(nimi))
			else:
				from random import randint
				for imi in range(nimi): boot[imi] = randint(0, nimi - 1)
			qa = model_blank(wn, wn)
			roo  = np.zeros(lenroo, np.float32)
			sroo = np.zeros(lenroo, np.float32)
			aroo = np.zeros(lenroo, np.float32)
			
			for imi in range(nimi):
				Util.add_img(qa, pw2[boot[imi]])
				temp1 = np.array(allroo[boot[imi]])
				roo += temp1
				temp2 = movingaverage(temp1, 10)
				aroo += temp2
				sroo += temp2**2
			sroo[0] = sroo[1]
			aroo[0] = aroo[1]
			sroo = (sroo-aroo**2 / nimi) / nimi
			aroo /= nimi
			roo  /= nimi
			qa   /= nimi
			
			if f_start < 0:
				#  Find a break point
				bp = 1.e23
				for i in range(5, lenroo - 5):
					#t1 = linreg(sroo[:i])
					#t2 = linreg(sroo[i:])
					#tt = t1[1][0] + t2[1][0]
					xtt = np.array(list(range(i)), np.float32)
					zet = np.poly1d(np.polyfit(xtt,sroo[:i], 2))
					t1 = sum((sroo[:i] - zet(xtt))**2)
					xtt = np.array(list(range(i, lenroo)), np.float32)
					zet = np.poly1d(np.polyfit(xtt, sroo[i:], 2) )
					tt = t1 + sum((sroo[i:] - zet(xtt))**2)
					if tt < bp:
						bp = tt
						istart = i
				#istart = 25
				#print istart
				f_start = istart / (pixel_size * wn)
			"""
			hi = hist_list(sroo,2)
			# hi[0][1] is the threshold
			for i in xrange(1,len(sroo)):
				if(sroo[i] < hi[0][1]):
					istart = i
					break
			"""
			#write_text_file([roo.tolist(),aroo.tolist(),sroo.tolist()], "sroo%03d.txt"%ifi)
			rooc = roo.tolist()

			#print namics[ifi],istart,f_start

			defc, subpw, ctf2, baseline, envelope, istart, istop = defocusgett_pap(rooc, wn, voltage = voltage, Pixel_size = pixel_size, Cs = Cs, ampcont = wgh, f_start = f_start, f_stop = f_stop, round_off = 1.0, nr1 = 3, nr2 = 6, parent = None, DEBug = debug_mode)
			#defc, subpw, ctf2, baseline, envelope, istart, istop = defocusgett(rooc, wn, voltage = voltage, Pixel_size = pixel_size, Cs = Cs, ampcont = wgh, f_start = f_start, f_stop = f_stop, round_off = 1.0, nr1 = 3, nr2 = 6, parent = None, DEBug = debug_mode)
			if debug_mode:
				print("  RESULT %s" % (img_name), defc, istart, istop)
				
				freq = list(range(len(subpw)))
				for i in range(len(freq)):  freq[i] = float(i) / wn / pixel_size
#				write_text_file([freq, subpw.tolist(), ctf2, envelope.tolist(), baseline.tolist()], "%s/ravg%05d.txt" % (output_directory, ifi))
				fou = os.path.join(outravg, "%s_ravg_%02d.txt" % (img_basename_root, nboot))
				write_text_file([freq, subpw.tolist(), ctf2, envelope.tolist(), baseline.tolist()], fou)
			#mpi_barrier(MPI_COMM_WORLD)
			
			#exit()
			bg = baseline.tolist()
			en = envelope.tolist()
			
			bckg = model_blank(wn, wn, 1, 1)
			envl = model_blank(wn, wn, 1, 1)
			
			from math import sqrt
			nc = wn // 2
			bg.append(bg[-1])
			en.append(en[-1])
			for i in range(wn):
				for j in range(wn):
					r = sqrt((i - nc)**2 + (j - nc)**2)
					ir = int(r)
					if(ir < nc):
						dr = r - ir
						bckg.set_value_at(i, j, (1. - dr) * bg[ir] + dr * bg[ir + 1] )
						envl.set_value_at(i, j, (1. - dr) * en[ir] + dr * en[ir + 1] )
			
			#qa.write_image("rs1.hdf")
			
			mask = model_circle(istop - 1, wn, wn) * (model_blank(wn, wn, 1, 1.0) - model_circle(istart, wn, wn))
			qse = threshold((qa - bckg))#*envl
			#(qse*mask).write_image("rs2.hdf")
			#qse.write_image("rs3.hdf")
			##  SIMULATION
			#bang = 0.7
			#qse = ctf2_rimg(wn, generate_ctf([defc,Cs,voltage,pixel_size,0.0,wgh, bang, 37.0]) )
			#qse.write_image("rs3.hdf")
			
			cnx = wn // 2 + 1
			cny = cnx
			mode = "H"
			istop = min(wn // 2 - 2, istop)    #2-26-2015@ming
			numr = Numrinit(istart, istop, 1, mode)
			wr = ringwe(numr, mode)
			
			crefim = Util.Polar2Dm(qse*mask, cnx, cny, numr, mode)
			Util.Frngs(crefim, numr)
			Util.Applyws(crefim, numr, wr)

			#pc = ctf2_rimg(wn,generate_ctf([defc,Cs,voltage,pixel_size,0.0,wgh]))
			#print ccc(pc*envl, subpw, mask)

			bang = 0.0
			bamp = 0.0
			bdef = defc
			bold = 1.e23
			while(True):
				#  in simctf2 data[3] is astigmatism amplitude
				"""
				data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
				#astdata = [crefim, numr, wn, bdef, Cs, voltage, pixel_size, wgh, bang]
				for qqq in xrange(200):
					qbdef = 1.0 + qqq*0.001
					print " VALUE AT THE BEGGINING OF while LOOP  ",qbdef,simctf2(qbdef, data)#,fastigmatism3(bamp,astdata)
				"""
				"""
				bamp = 0.7
				bang = 37.0
				
				data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
				astdata = [crefim, numr, wn, bdef, Cs, voltage, Pixel_size, wgh, bang]
				print " VALUE AT THE BEGGINING OF while LOOP  ",bdef,bamp,bang,simctf2(bdef, data),fastigmatism3(bamp,astdata,mask)
				#print  simctf2out(1.568,data)
				#exit()
				
				for kdef in xrange(14000,17000,10):
					dz = kdef/10000.0
					ard = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
					#print ard
					aqd = [crefim, numr, wn, dz, Cs, voltage, pixel_size, wgh, bang]
					#print aqd
					print  dz,simctf2(dz,ard),fastigmatism3(bamp,aqd,mask)
					#print aqd[-1]
				exit()
				"""
				data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
				h = 0.05 * bdef
				amp1, amp2 = bracket_def(simctf2_pap, data, bdef * 0.9, h)
				#amp1, amp2 = bracket_def(simctf2, data, bdef * 0.9, h)
				#print "bracketing of the defocus  ",amp1, amp2
				#print " ttt ",time()-srtt
				#print "bracketing of the defocus  ",amp1,amp2,simctf2(amp1, data),simctf2(amp2, data),h
				amp1, val2 = goldsearch_astigmatism(simctf2_pap, data, amp1, amp2, tol = 1.0e-3)
				#amp1, val2 = goldsearch_astigmatism(simctf2, data, amp1, amp2, tol = 1.0e-3)
				#print "golden defocus ",amp1, val2,simctf2(amp1, data)
				#bdef, bcc = goldsearch_astigmatism(simctf2, data, amp1, amp2, tol=1.0e-3)
				#print "correction of the defocus  ",bdef,bcc
				#print " ttt ",time()-srtt
				"""
				crot2 = rotavg_ctf(ctf2_rimg(wn,generate_ctf([bdef, Cs, voltage, pixel_size, 0.0, wgh, bamp, bang])), bdef, Cs, voltage, pixel_size, bamp, bang)
				pwrot = rotavg_ctf(qa-bckg, bdef, Cs, voltage, pixel_size, bamp, bang)
				write_text_file([range(len(subroo)),asubroo, ssubroo, sen, pwrot, crot2],"rotinf%04d.txt"%ifi)
				qse.write_image("qse.hdf")
				mask.write_image("mask.hdf")
				exit()
				"""

				astdata = [crefim, numr, wn, bdef, Cs, voltage, pixel_size, wgh, bang, mask]
				h = 0.01
				amp1, amp2 = bracket(fastigmatism3_pap, astdata, h)
				#amp1, amp2 = bracket(fastigmatism3, astdata, h)
				#print "  astigmatism bracket  ",amp1,amp2,astdata[-1]
				#print " ttt ",time()-srtt
				bamp, bcc = goldsearch_astigmatism(fastigmatism3_pap, astdata, amp1, amp2, tol = 1.0e-3)
				junk = fastigmatism3_pap(bamp,astdata)
				#bamp, bcc = goldsearch_astigmatism(fastigmatism3, astdata, amp1, amp2, tol = 1.0e-3)
				#junk = fastigmatism3(bamp,astdata)
				bang = astdata[8]
				
				#print astdata[8]
				#print  fastigmatism3(0.0,astdata)
				#print astdata[8]
				#temp = 0.0
				#print bdef, Cs, voltage, pixel_size, temp, wgh, bamp, bang, -bcc
				#data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
				#astdata = [crefim, numr, wn, bdef, Cs, voltage, pixel_size, wgh, bang]
				#print " VALUE WITHIN the while LOOP  ",bdef,bamp,bang,simctf2(bdef, data),fastigmatism3(bamp,astdata)
				#print "  golden search ",bamp,data[-1], fastigmatism3(bamp,data), fastigmatism3(0.0,data)
				#print " ttt ",time()-srtt
				#bamp = 0.5
				#bang = 277
				
				dama = amoeba([bdef, bamp], [0.2, 0.2], fupw_pap, 1.e-4, 1.e-4, 500, astdata)
				#dama = amoeba([bdef, bamp], [0.2, 0.2], fupw, 1.e-4, 1.e-4, 500, astdata)
				if debug_mode:  print("AMOEBA    ", dama)
				bdef = dama[0][0]
				bamp = dama[0][1]
				astdata = [crefim, numr, wn, bdef, Cs, voltage, pixel_size, wgh, bang, mask]
				junk = fastigmatism3_pap(bamp, astdata)
				#junk = fastigmatism3(bamp, astdata)
				bang = astdata[8]
				if debug_mode:  print(" after amoeba ", bdef, bamp, bang)
				#  The looping here is blocked as one shot at amoeba is good enough.  To unlock it, remove - from bold.
				if(bcc < -bold): bold = bcc
				else:           break

			#data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
			#print " VALUE AFTER the while LOOP  ",bdef,bamp,bang,simctf2(bdef, data),fastigmatism3(bamp,astdata)
			#temp = 0.0
			#print ifi,bdef, Cs, voltage, pixel_size, temp, wgh, bamp, bang, -bcc
			#freq = range(len(subpw))
			#for i in xrange(len(freq)):  freq[i] = float(i)/wn/pixel_size
			#ctf2 = ctf_2(wn, generate_ctf([bdef,Cs,voltage,pixel_size,0.0,wgh]))[:len(freq)]
			#write_text_file([freq, subpw.tolist(), ctf2, envelope.tolist(), baseline.tolist()],"ravg/ravg%05d.txt"%ifi)
			#print " >>>> ",wn, bdef, bamp, Cs, voltage, pixel_size, wgh, bang
			#data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
			#print  simctf2out(bdef, data)
			#exit()
			adefocus[nboot] = bdef
			aamplitu[nboot] = bamp
			aangle[nboot]   = bang
			#from sys import exit
			#exit()

		#print " ttt ",time()-srtt
		#from sys import exit
		#exit()
		ad1, ad2, ad3, ad4 = table_stat(adefocus) # return values: average, variance, minimum, maximum
		if ad2 <= 0.0:
			print(("    %s %s: Detected the variance less than zero (defocus statistics: avg = %f, var = %f, min = %f, max = %f)." % (img_type, img_name, ad1, ad2, ad3, ad4)))
			print(("           The program ignores this estimate..."))
			continue
		
		reject = []
		thr = 3 * sqrt(ad2)
		for i in range(len(adefocus)):
			if(abs(adefocus[i] - ad1) > thr):
				print(("    %s %s: Rejected an outlier defocus estimate (defocus = %f, average defocus = %f, threshold = %f)." % (img_type, img_name, adefocus[i], ad1, thr)))
				reject.append(i)
		
		if(len(reject) > 0):
			print(("    %s %s: Total number of rejects %s" % (img_type, img_name, len(reject))))
			for i in range(len(reject) - 1, -1, -1):
				del adefocus[i]
				del aamplitu[i]
				del aangle[i]
		
		if(len(adefocus) < 2):
			print(("    %s %s: After rejection of outliers, there is too few estimated defocus values. Skipping the estimation and CTF parameters are not stored..." % (img_type, img_name)))
		else:
			#print "adefocus",adefocus
			#print  "aamplitu",aamplitu
			#print "aangle",aangle
			ad1, ad2, ad3, ad4 = table_stat(adefocus)
			bd1, bd2, bd3, bd4 = table_stat(aamplitu)
			cd1,cd2 = angle_ave([2*q for q in aangle])  # Have to use this trick as the function works for range [0,360]
			cd1/=2
			cd2/=2
			temp = 0.0
			stdavad1 = np.sqrt(kboot * max(0.0, ad2))
			stdavbd1 = np.sqrt(kboot * max(0.0, bd2))
			cd2 *= np.sqrt(kboot)
			
			# Adjust value ranges of astig. amp. and angle.
			if( bd1 < 0.0 ):
				bd1 = -bd1
				cd1 += 90.0
			cd1 = cd1%180
			
			#  SANITY CHECK, do not produce anything if defocus abd astigmatism amplitude are out of whack
			reject_img_messages = []
			try:
				pwrot2 = rotavg_ctf( model_blank(wn, wn), ad1, Cs, voltage, pixel_size, bd1, cd1)
			except:
				reject_img_messages.append("    - Astigmatism amplitude (%f) is larger than defocus (%f) or defocus (%f) is negative." % (bd1, ad1, ad1))
			
			if len(reject_img_messages) > 0:
				rejected_img_names.append(img_name)
				print("    %s %s: Rejected the CTF estimate - " % (img_type, img_name), ad1, Cs, voltage, pixel_size, wgh, bd1, cd1, "(def, Cs, vol, apix, amp_contrast, astig_amp, astig_angle)")
				print("    %s %s: because... " % (img_type, img_name))
				for reject_img_message in reject_img_messages:
					print(reject_img_message)
				print("    %s %s: Skipping the estimation and CTF parameters are not stored..." % (img_type, img_name))
			else:
				#  Estimate the point at which (sum_errordz ctf_1(dz+errordz))^2 falls to 0.5
				import random as rqt
				
				supe = model_blank(wn, wn)
				niter = 1000
				for it in range(niter):
					Util.add_img(supe, Util.ctf_rimg(wn, wn, 1, ad1 + rqt.gauss(0.0,stdavad1), pixel_size, voltage, Cs, 0.0, wgh, bd1 + rqt.gauss(0.0,stdavbd1), cd1 + rqt.gauss(0.0,cd2), 1))
				ni = wn // 2
				supe /= niter
				pwrot2 = rotavg_ctf(supe, ad1, Cs, voltage, pixel_size, bd1, cd1)
				for i in range(ni):  pwrot2[i] = pwrot2[i]**2
				
				ibec = 0
				for it in range(ni - 1, 0, -1):
					if pwrot2[it] > 0.5 :
						ibec = it
						break
				from morphology import ctf_1d
				ct = generate_ctf([ad1, Cs, voltage, pixel_size, temp, wgh, 0.0, 0.0])
				cq = ctf_1d(wn, ct)
				
				supe = [0.0] * ni
				niter = 1000
				for i in range(niter):
					cq = generate_ctf([ad1 + rqt.gauss(0.0,stdavad1), Cs, voltage, pixel_size, 0.0, wgh, 0.0, 0.0])
					ci = ctf_1d(wn, cq)[:ni]
					for l in range(ni):  supe[l] +=ci[l]
				
				for l in range(ni):  supe[l] = (supe[l] / niter)**2
				
				ib1 = 0
				for it in range(ni - 1, 0, -1):
					if supe[it] > 0.5:
						ib1 = it
						break
				ibec = ibec / (pixel_size * wn)  #  with astigmatism
				ib1  = ib1 / (pixel_size * wn)   #  no astigmatism
				#from utilities import write_text_file
				#write_text_file([range(ni), supe[:ni],pwrot2[:ni]],"fifi.txt")
				
				# Compute defocus CV and astig. amp. CV (CV: coefficient of variation; ratio of error (SD) relative to average (mean))
				cvavad1 = stdavad1 / ad1 * 100 # use percentage
				
				cvavbd1 = stdavbd1 / bd1 * 100 # use percentage
				
				# Compute CTF limit (theoretical resolution limit based on the oscillations of CTF) 
				# For output, use ctflim (relative frequency limit [1/A]), not ctflim_abs (absolute frequency limit)
				# 
				# NOTE: 2016/03/23 Toshio Moriya
				# xr is limiting frequency [1/A]. Max is Nyquist frequency = 1.0/(2*apix[A/pixel]). <UNIT: [1/(A/pixel)/[pixel])] => [(pixel)/(A*pixel] => [1/A]>
				# 1.0/xr is limiting period (Angstrom resolution) [A]. Min is Nyquist period = (2*apix[A/pixel]). <UNIT: [1/(1/A)] = [A]>
				# fwpix is width of Fourier pixel [pixel/A] := 1.0[pixel]/(2*apix[A/pixel])/box_half[pixel] = 1[pixel]/fullsize[A]). <UNIT: [pixel/(A/pixel)/(pixel)] = [pixel*(pixel/A)*(1/pixel) = [pixel/A]>
				# int(xr/fwpix+0.5) is limiting_absolute_frequency [1/pixel]. <Unit:[(1/A)/(pixel/A)] = [(1/A)*(A/pixel)] = [1/pixel]>
				# return  int(xr/fwpix+0.5),xr, which is limiting_abs_frequency [1/pixel], and Limiting_frequency[1/A]
				#
				ctflim_abs, ctflim = ctflimit(wn, ad1, Cs, voltage, pixel_size)
				
				"""
				for i in xrange(len(ssubroo)):
					asubroo[i] /= kboot
					ssubroo[i]  = sqrt(max(0.0, ssubroo[i]-kboot*asubroo[i]**2)/kboot)
					sen[i]     /= kboot
				"""
				lnsb = len(subpw)
				try:		crot1 = rotavg_ctf(ctf_rimg(wn, generate_ctf([ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1]), sign=0), ad1, Cs, voltage, pixel_size, 0.0, 0.0)[:lnsb]
				except:		crot1 = [0.0] * lnsb
				try:		pwrot1 = rotavg_ctf(threshold(qa - bckg), ad1, Cs, voltage, pixel_size, 0.0, 0.0)[:lnsb]
				except:		pwrot1 = [0.0] * lnsb
				try:		crot2 = rotavg_ctf(ctf_rimg(wn, generate_ctf([ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1]), sign=0), ad1, Cs, voltage, pixel_size, bd1, cd1)[:lnsb]
				except:		crot2 = [0.0] * lnsb
				try:		pwrot2 = rotavg_ctf(threshold(qa - bckg), ad1, Cs, voltage, pixel_size, bd1, cd1)[:lnsb]
				except:		pwrot2 = [0.0] * lnsb
				#  #1 - rotational averages without astigmatism, #2 - with astigmatism
				lnsb = min(len(crot2),len(pwrot1),len(crot2),len(pwrot2))
				write_text_file([list(range(lnsb)), [float(i)/wn/pixel_size for i in range(lnsb)], pwrot1, crot1, pwrot2, crot2], os.path.join(outpwrot, "%s_rotinf.txt"%(img_basename_root)))
				
				#
				# NOTE: 2016/03/23 Toshio Moriya
				# Compute mean of extrema differences (differences at peak & trough) between 
				# (1) experimental rotational average with astigmatism (pwrot2)
				# (2) experimental rotational average without astigmatism (pwrot1), and
				# as a indication of goodness of astigmatism estimation by cter.
				# The peak & trough detection uses fitted rotational average with astigmatism (crot2) 
				# Start from 1st trough while ignoring 1st peak.
				# End at astigmatism frequency limit.
				# 
				is_peak_target = True
				pre_crot2_val = crot2[0]
				extremum_counts = 0
				extremum_diff_sum = 0
				for i in range(1, len(crot2)):
					cur_crot2_val = crot2[i]
					if is_peak_target == True and pre_crot2_val > cur_crot2_val:
						# peak search state
						extremum_i = i - 1
						extremum_counts += 1
						extremum_diff_sum += pwrot2[extremum_i] - pwrot1[extremum_i] # This should be positive if astigmatism estimation is good
						# print "MRK_DEBUG: Peak Search  : extremum_i = %03d, freq[extremum_i] = %12.5g, extremum_counts = %03d, (pwrot2[extremum_i] - pwrot1[extremum_i]) = %12.5g, extremum_diff_sum = %12.5g " % (extremum_i, freq[extremum_i] , extremum_counts, (pwrot2[extremum_i] - pwrot1[extremum_i]), extremum_diff_sum)
						is_peak_target = False
					elif is_peak_target == False and pre_crot2_val < cur_crot2_val:
						# trough search state
						extremum_i = i - 1
						extremum_counts += 1
						extremum_diff_sum += pwrot1[extremum_i] - pwrot2[extremum_i] # This should be positive if astigmatism estimation is good
						# print "MRK_DEBUG: Trough Search: extremum_i = %03d, freq[extremum_i] = %12.5g, extremum_counts = %03d, (pwrot1[extremum_i] - pwrot2[extremum_i]) = %12.5g, extremum_diff_sum = %12.5g " % (extremum_i, freq[extremum_i] , extremum_counts, (pwrot1[extremum_i] - pwrot2[extremum_i]), extremum_diff_sum)
						is_peak_target = True
					pre_crot2_val = cur_crot2_val
#				#if extremum_counts == 0: ERROR("Logical Error: Encountered unexpected zero extremum counts. Consult with the developer." % (bd1), "%s in %s" % (__name__, os.path.basename(__file__))) # MRK_ASSERT
#				extremum_diff_avg = 1.1#extremum_diff_sum / extremum_counts
				
				# print "MRK_DEBUG: extremum_avg = %12.5g, extremum_diff_sum = %12.5g, extremum_counts = %03d," % (extremum_avg, extremum_diff_sum, extremum_counts)
				
#				if stack == None:     cmd = "echo " + "    " + namics[ifi] + "  >>  " + fou
#				else:                 cmd = "echo " + "    " + "  >>  " + fou
#				os.system(cmd)
				
				ed1      = wgh            # total amplitude contrast. Here, it is also constant amplitude contrast since Volta phase shift is not estimated with this function. 
				stdaved1 = 0.0            # dummy value for error of total amplitude contrast estimation
				max_freq = 0.5/pixel_size # dummy value for maximum frequency. set to Nyquist frequency for now. let's add the implementation in near future (Toshio 2017/12/06)
				reserved = 0.0            # dummy value for reserved spot, which might be used for parameter of external programs (e.g. CTFFIND4, GCTF, and etc.)
				# wgh                     # constant amplitude contrast provided by user (default 10%). Here, it is also total amplitude contrast since Volta phase shift is not estimated with this function.
				phase_shift = ampcont2angle(ed1) - ampcont2angle(wgh) # Volta phase shift [deg] = total amplitude contrast phase shift [deg] (ed1) -  constant amplitude contrast phase shift [deg]; ed1 is boot strap average of total amplitude contrast [%]
				
				if debug_mode: print(("    %s %s: Process %04d finished the processing. Estimated CTF parmaters are stored in %s." % (img_type, img_name, ifi, os.path.join(output_directory, "partres.txt"))))
#				if debug_mode: print((ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, cvavad1, cvavbd1, extremum_diff_avg, ib1, ibec, ctflim))
#				# totresi.append([img_name, ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, cvavad1, cvavbd1, extremum_diff_avg, ib1, ibec, ctflim])
#				totresi.append([img_name, ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdaved1, stdavbd1, cd2, cvavad1, cvavbd1, extremum_diff_avg, ib1, ibec, ctflim])
				if debug_mode: print((ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdaved1, stdavbd1, cd2, cvavad1, cvavbd1, ib1, ibec, ctflim, max_freq, reserved, wgh, phase_shift))
				totresi.append([img_name, ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1, stdavad1, stdaved1, stdavbd1, cd2, cvavad1, cvavbd1, ib1, ibec, ctflim, max_freq, reserved, wgh, phase_shift])
				
#				if stack == None:
#					print  namics[ifi], ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec
#				else:
#					print               ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec
#				if stack == None:
#					totresi.append( [ namics[ifi], ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec])
#				else:
#					totresi.append( [ 0, ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec])
#				#if ifi == 4 : break

		if stack == None:
			img_mic = get_im(namics[ifi])
			# create  thumbnail
			nx = img_mic.get_xsize()
			if nx > 512:
				img_micthumb = resample(img_mic, 512.0/nx)
			else:
				img_micthumb = img_mic
			fou = os.path.join(outmicthumb, "%s_thumb.hdf" % (img_basename_root))
			img_micthumb.write_image(fou)

	if RUNNING_UNDER_MPI:
		from utilities import wrap_mpi_gatherv
		totresi = wrap_mpi_gatherv(totresi, 0, MPI_COMM_WORLD)
		missing_img_names = wrap_mpi_gatherv(missing_img_names, 0, MPI_COMM_WORLD)
		rejected_img_names = wrap_mpi_gatherv(rejected_img_names, 0, MPI_COMM_WORLD)
	
	if my_mpi_proc_id == main_mpi_proc:
		outf = open(os.path.join(output_directory, "partres.txt"), "w")
		for i in range(len(totresi)):
			for k in range(1, len(totresi[i])):
				outf.write("  %12.5g" % totresi[i][k])
			outf.write("  %s\n" % totresi[i][0])
		outf.close()
		
		print(" ")
		print(("Summary of %s processing..." % (img_type.lower())))
		missing_counts = len(missing_img_names)
		print(("  Missing  : %d" % (missing_counts)))
		if missing_counts > 0:
			outfile_path = os.path.join(output_directory, "missing_%s_list.txt" % (img_type.lower()))
			print(("    Saving list of missing in %s..." % (outfile_path)))
			outf = open(outfile_path, "w")
			for missing_img_name in missing_img_names:
				outf.write("%s\n" % missing_img_name)
			outf.close()
		
		rejected_counts = len(rejected_img_names)
		print(("  Rejected : %d" % (rejected_counts)))
		if rejected_counts > 0:
			outfile_path = os.path.join(output_directory, "rejected_%s_list.txt" % (img_type.lower()))
			print(("    Saving list of rejected in %s..." % (outfile_path)))
			outf = open(outfile_path, "w")
			for rejected_img_name in rejected_img_names:
				outf.write("%s\n" % rejected_img_name)
			outf.close()
	
	if cter_mode_idx == idx_cter_mode_stack:
		return totresi[0][1], totresi[0][7], totresi[0][8], totresi[0][9], totresi[0][10], totresi[0][11]



################
#
#  CTER code (new version since 2016/03/16)
#
################
# 
# NOTE: 2016/03/16 Toshio Moriya
# In this version, the IO-related interface is simplified for sxcter.py and sxgui.py
# Since cter() was used in not only sxcter.py but also e2boxer.py and sxhelixboxer.py,
# This new version is added to avoid breaking e2boxer.py and sxhelixboxer.py
#  
# NOTE: 2016/03/16 Toshio Moriya
# To get a single  file name from a GUI application, 
# there must be a better way than using guimic...
# 
# NOTE: 2016/11/16 Toshio Moriya
# Now, this function assume the MPI setup and clean up is done by caller, such as mpi_init, and mpi_finalize
# 07/11/2017  This is power spectrum version of cter_mrk
def cter_pap(input_image_path, output_directory, selection_list = None, wn = 512, \
			pixel_size = -1.0, Cs = 2.0, voltage = 300.0, wgh = 10.0, f_start = -1.0, f_stop = -1.0, \
			kboot = 16, overlap_x = 50, overlap_y = 50, edge_x = 0, edge_y = 0, \
			check_consistency = False, stack_mode = False, debug_mode = False, \
			program_name = "cter_pap() in morphology.py", \
			RUNNING_UNDER_MPI = False, main_mpi_proc = 0, my_mpi_proc_id = 0, n_mpi_procs = 1):
	"""
	Arguments
		input_image_path  :  file name pattern for Micrographs Modes (e.g. 'Micrographs/mic*.mrc') or particle stack file path for Stack Mode (e.g. 'bdb:stack'; must be stack_mode = True).
		output_directory  : output directory
	"""
	from   EMAN2 import periodogram
	from   EMAN2db import db_check_dict, db_parse_path
	from   applications import MPI_start_end
	from   utilities import read_text_file, write_text_file, get_im, model_blank, model_circle, amoeba, generate_ctf
	from   utilities import if_error_then_all_processes_exit_program
	from   utilities import wrap_mpi_bcast
	from   sys import exit
	import numpy as np
	import os
	import glob
	from   fundamentals import tilemic, rot_avg_table, resample
	from   morphology   import threshold, bracket_def, bracket, goldsearch_astigmatism
	from   morphology   import defocus_baseline_fit, simpw1d, movingaverage, localvariance, defocusgett
	from   morphology   import defocus_guessn, defocusget_from_crf, make_real
	from   morphology   import fastigmatism, fastigmatism1, fastigmatism2, fastigmatism3, simctf, simctf2, simctf2out, fupw,ctf2_rimg
	from   alignment    import Numrinit, ringwe
	from   statistics   import table_stat
	from   pixel_error  import angle_ave
	from   inspect      import currentframe, getframeinfo
	from   global_def   import ERROR
	import global_def
	from   time import time
	from   mpi import MPI_COMM_WORLD, mpi_barrier

	# ====================================================================================
	# Prepare processing
	# ====================================================================================

	
	# ------------------------------------------------------------------------------------
	# Find the CTER Running Mode before checking error conditions
	# ------------------------------------------------------------------------------------
	i_enum = -1; idx_cter_mode_invalid       = i_enum; 
	i_enum += 1; idx_cter_mode_all_mics      = i_enum  # All Micrographs Mode - Process all s in a directory
	i_enum += 1; idx_cter_mode_selected_mics = i_enum  # Selected Micrographs Mode - Process all s in a selection list file
	i_enum += 1; idx_cter_mode_single_mic    = i_enum  # Single Micrograph Mode - Process a single 
	i_enum += 1; idx_cter_mode_stack         = i_enum  # Stack Mode - Process a stack (Advanced Option)
	i_enum += 1; idx_cter_mode_counts        = i_enum

	cter_mode_idx = idx_cter_mode_invalid
	cter_mode_name = None
	if stack_mode == False:
		# One of three Micrograph Modes
		# For any of Micrograph Modes, input image file name must be a file path pattern containing wild card "*" 
		if selection_list == None:
			# User did not use selection list option 
			# -> All Micrographs Mode
			cter_mode_idx = idx_cter_mode_all_mics
			cter_mode_name = "All Micrographs Mode"
		else:
			if os.path.splitext(selection_list)[1] == ".txt":
				# User specified a selection list text file path containing".txt" extension through selection list option
				# -> Selected Micrographs Mode
				cter_mode_idx = idx_cter_mode_selected_mics
				cter_mode_name = "Selected Micrographs Mode"
			else: 
				# User specified an image file path (a non-text file path) through selection list option
				# -> Single Micrograph Mode
				cter_mode_idx = idx_cter_mode_single_mic
				cter_mode_name = "Single Micrograph Mode"
	else: 
		# (Particle) Stack Mode
		cter_mode_idx = idx_cter_mode_stack
		cter_mode_name = "Stack Mode"

	if my_mpi_proc_id == main_mpi_proc:
		print(" ")
		print(("----- Running with %s -----" % (cter_mode_name)))

	# ------------------------------------------------------------------------------------
	# Check mode-dependent error conditions of input arguments and options if abort is necessary. All nodes do this checking
	# ------------------------------------------------------------------------------------
	error_message_list = [] # List of error messages. If no error is found, the length should be zero
	if not stack_mode:

		# Check error conditions applicable to any of Micrograph Mode 
		if input_image_path.find("*") == -1:
			error_message_list.append("Input image file path (%s) for %s must be a  path pattern containing wild card (*). Please check input_image_path argument." % (input_image_path, cter_mode_name))

		if input_image_path[:len("bdb:")].lower() == "bdb:":
			error_message_list.append("BDB file can not be selected as input image file path (%s) for %s. Please check input_image_path argument and convert the image format." % (input_image_path, cter_mode_name))

		# Check error conditions applicable to Selected Micrographs Mode 
		if cter_mode_idx == idx_cter_mode_selected_mics:
			if not os.path.exists(selection_list): 
				error_message_list.append("Selection list text file specified by selection_list option (%s) for %s does not exists. Please check selection_list option." % (selection_list, cter_mode_name))

		if cter_mode_idx == idx_cter_mode_single_mic:
			if not os.path.exists(os.path.join(os.path.dirname(input_image_path), os.path.basename(selection_list))): 
				error_message_list.append("Micrograph specified by selection_list option (%s) for %s does not exist. Please check selection_list option." % (selection_list, cter_mode_name))
			# 
			if RUNNING_UNDER_MPI and n_mpi_procs != 1:
				error_message_list.append("%s supports only a single processor version. Please change MPI settings." % (cter_mode_name))

	else: 
		# Check error conditions
		if input_image_path.find("*") != -1:
			error_message_list.append("Stack file path specified by input_image_path (%s) for %s should not contain wild card (*). Please check input_image_path argument." % (input_image_path, cter_mode_name))

		is_not_found_input_image_file = False
		if input_image_path[:len("bdb:")].lower() == "bdb:":
			if not db_check_dict(input_image_path): 
				is_not_found_input_image_file = True
		else:
			if not os.path.exists(input_image_path): 
				is_not_found_input_image_file = True
		if is_not_found_input_image_file:
			error_message_list.append("Stack file specified by input_image_path (%s) for %s does not exist. Please check input_image_path argument." % (input_image_path, cter_mode_name))

		if RUNNING_UNDER_MPI and n_mpi_procs != 1:
			error_message_list.append("%s supports only a single processor version. Please change MPI settings." % (cter_mode_name))

	# --------------------------------------------------------------------------------
	# check output-related error conditions (mode-independent). All nodes do this checking
	# --------------------------------------------------------------------------------
	if os.path.exists(output_directory):
		error_message_list.append("Output directory (%s) exists already. Please check output_directory argument." % (output_directory))
	
	# --------------------------------------------------------------------------------
	# Check error conditions of options (mode-independent). All nodes do this checking
	# --------------------------------------------------------------------------------
	if pixel_size <= 0.0:
		error_message_list.append("Pixel size (%f) must not be negative. Please set a pasitive value larger than 0.0 to pixel_size option." % (pixel_size))

	if wn <= 0.0:
		error_message_list.append("CTF window size (%d) must not be negative. Please set a valid value larger than 0 to wn option." % (wn))
	
	# --------------------------------------------------------------------------------
	# Print all error messages and abort the process if necessary.
	# --------------------------------------------------------------------------------
	error_status = None
	if len(error_message_list) > 0:
		# Detected error! Print all error messages
		if my_mpi_proc_id == main_mpi_proc:
			print(" ")
			for error_message in error_message_list:  
				print(("ERROR!!! %s" % (error_message)))
		error_status = ("Detected %d error(s) related to arguments and options. Run %s -h for help. Exiting..." % (len(error_message_list), program_name), getframeinfo(currentframe()))
	if_error_then_all_processes_exit_program(error_status)
	if RUNNING_UNDER_MPI:
		# Wait for all mpi processes to check error conditions, especially existence of output directory
		# Without this barrier, main mpi process can create output directory before some child mpi process check this error.
		mpi_barrier(MPI_COMM_WORLD)

	del error_message_list # Don't need this anymore
	
	# ------------------------------------------------------------------------------------
	# Check warning conditions of options
	# ------------------------------------------------------------------------------------
	if my_mpi_proc_id == main_mpi_proc:
		if stack_mode:
			if selection_list != None:
				print(" ")
				print(("WARNING!!! --selection_list option will be ignored in %s." % (cter_mode_name)))
			if wn != 512:
				print(" ")
				print(("WARNING!!! --wn option will be ignored in %s." % (cter_mode_name)))
			if overlap_x != 50:
				print(" ")
				print(("WARNING!!! --overlap_x option will be ignored in %s." % (cter_mode_name)))
			if overlap_y != 50:
				print(" ")
				print(("WARNING!!! --overlap_y option will be ignored in %s." % (cter_mode_name)))
			if edge_x != 0:
				print(" ")
				print(("WARNING!!! --edge_x option will be ignored in %s." % (cter_mode_name)))
			if edge_y != 0:
				print(" ")
				print(("WARNING!!! --edge_y option will be ignored in %s." % (cter_mode_name)))
			if check_consistency:
				print(" ")
				print(("WARNING!!! --check_consistency option will be ignored in %s." % (cter_mode_name)))

	# ====================================================================================
	# Create the input file path list and also check input-related error conditions if abort is necessary.
	# ====================================================================================
	input_file_path_list = []
	if not stack_mode:
		# --------------------------------------------------------------------------------
		# Prepare the variables for all sections in  mode case 
		# --------------------------------------------------------------------------------
		# Micrograph basename pattern (directory path is removed from  path pattern)
		mic_pattern = input_image_path
		mic_basename_pattern = os.path.basename(mic_pattern)
	
		# Global entry dictionary (all possible entries from all lists) for all mic id substring
		global_entry_dict = {} # mic id substring is the key
		subkey_input_mic_path = "Input Micrograph Path"
		subkey_selected_mic_basename = "Selected Micrograph Basename"
	
		# List keeps only id substrings of s whose all necessary information are available
		valid_mic_id_substr_list = [] 
		
		# --------------------------------------------------------------------------------
		# Obtain the list of  id sustrings using a single CPU (i.e. main mpi process)
		# --------------------------------------------------------------------------------
		# NOTE: Toshio Moriya 2016/11/15
		# The below is not a real while.  
		# It gives if-statements an opportunity to use break when errors need to be reported
		# However, more elegant way is to use 'raise' statement of exception mechanism...
		# 
		error_status = None
		while my_mpi_proc_id == main_mpi_proc:
			
			# --------------------------------------------------------------------------------
			# Prepare variables for this section
			# --------------------------------------------------------------------------------
			# Prefix and suffix of  basename pattern 
			# to find the head/tail indices of  id substring
			mic_basename_tokens = mic_basename_pattern.split('*')
			# Find head index of  id substring
			mic_id_substr_head_idx = len(mic_basename_tokens[0])
		
			# --------------------------------------------------------------------------------
			# Register  id substrings found in the input directory (specified by  path pattern)
			# to the global entry dictionary
			# --------------------------------------------------------------------------------
			# Generate the list of  paths in the input directory
			print(" ")
			print("Checking the input directory...")
			input_mic_path_list = glob.glob(mic_pattern)
			# Check error condition of input  file path list
			print(("Found %d micrographs in %s." % (len(input_mic_path_list), os.path.dirname(mic_pattern))))
			if len(input_mic_path_list) == 0:
				# The result shouldn't be empty if the specified  file name pattern is invalid
				error_status = ("There are no micrographs whose paths match with the specified file path pattern (%s) for %s. Please check input_image_path. Run %s -h for help." % (mic_pattern, cter_mode_name, program_name), getframeinfo(currentframe()))
				break
		
			# Register  id substrings to the global entry dictionary
			for input_mic_path in input_mic_path_list:
				# Find tail index of  id substring and extract the substring from the  name
				input_mic_basename = os.path.basename(input_mic_path)
				mic_id_substr_tail_idx = input_mic_basename.index(mic_basename_tokens[1])
				mic_id_substr = input_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
				if not mic_id_substr in global_entry_dict:
					# print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from input_mic_path_list " % (mic_id_substr))
					global_entry_dict[mic_id_substr] = {}
				global_entry_dict[mic_id_substr][subkey_input_mic_path] = input_mic_path
		
			# --------------------------------------------------------------------------------
			# Register  id substrings found in the selection list
			# to the global entry dictionary
			# --------------------------------------------------------------------------------
			# Generate the list of selected  paths in the selection file
			selected_mic_path_list = []
			# Generate  lists according to the execution mode
			if cter_mode_idx == idx_cter_mode_all_mics:
				# Treat all s in the input directory as selected ones
				selected_mic_path_list = input_mic_path_list
			else:
				if os.path.splitext(selection_list)[1] == ".txt":
					print(" ")
					print("Checking the selection list...")
					selected_mic_path_list = read_text_file(selection_list)
				
					# Check error condition of  entry lists
					print(("Found %d microgarph entries in %s." % (len(selected_mic_path_list), selection_list)))
					if len(selected_mic_path_list) == 0:
						error_status = ("The provided  list file (%s) for %s mode contains no entries. Please check selection_list option and make sure the file contains a  list. Run %s -h for help." % (selection_list, cter_mode_name, program_name), getframeinfo(currentframe()))
						break
				else:
					print(" ")
					print(("Processing a single micorgprah: %s..." % (selection_list)))
					selected_mic_path_list = [selection_list]
			
				selected_mic_directory = os.path.dirname(selected_mic_path_list[0])
				if selected_mic_directory != "":
					print(("    NOTE: Program disregards the directory paths in the selection list (%s)." % (selected_mic_directory)))
		
			# Register  id substrings to the global entry dictionary
			for selected_mic_path in selected_mic_path_list:
				# Find tail index of  id substring and extract the substring from the  name
				selected_mic_basename = os.path.basename(selected_mic_path)
				mic_id_substr_tail_idx = selected_mic_basename.index(mic_basename_tokens[1])
				mic_id_substr = selected_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
				if not mic_id_substr in global_entry_dict:
					# print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from selected_mic_path_list " % (mic_id_substr))
					global_entry_dict[mic_id_substr] = {}
				global_entry_dict[mic_id_substr][subkey_selected_mic_basename] = selected_mic_basename
			
			# --------------------------------------------------------------------------------
			# Clean up variables related to registration to the global entry dictionary
			# --------------------------------------------------------------------------------
			del mic_basename_tokens
			del mic_id_substr_head_idx
		
			# --------------------------------------------------------------------------------
			# Create the list containing only valid  id substrings
			# --------------------------------------------------------------------------------
			# Prepare lists to keep track of invalid (rejected) s 
			no_input_mic_id_substr_list = []
		
			print(" ")
			print("Checking the input datasets consistency...")
		
			# Loop over substring id list
			for mic_id_substr in global_entry_dict:
				mic_id_entry = global_entry_dict[mic_id_substr]
			
				warinnig_messages = []
				# selected  basename must have been registed always .
				if subkey_selected_mic_basename in mic_id_entry: 
					# Check if associated input  exists
					if not subkey_input_mic_path in mic_id_entry:
						input_mic_path = mic_pattern.replace("*", mic_id_substr)
						warinnig_messages.append("    associated input  %s." % (input_mic_path))
						no_input_mic_id_substr_list.append(mic_id_substr)
				
					if len(warinnig_messages) > 0:
						print(("WARNING!!! Micrograph ID %s does not have:" % (mic_id_substr)))
						for warinnig_message in warinnig_messages:
							print(warinnig_message)
						print("    Ignores this as an invalid entry.")
					else:
						# print("MRK_DEBUG: adding mic_id_substr := ", mic_id_substr)
						valid_mic_id_substr_list.append(mic_id_substr)
				# 	# This entry is not in the selection list. Do nothing
			
			# Check the input dataset consistency and save the result to a text file, if necessary.
			if check_consistency:
				# Create output directory
				os.mkdir(output_directory)
			
				# Open the consistency check file
				inconsist_mic_list_path = os.path.join(output_directory,"inconsist_mic_id_file.txt")
				print(" ")
				print(("Generating the input datasets consistency report in %s..." % (inconsist_mic_list_path)))
				inconsist_mic_list_file = open(inconsist_mic_list_path, "w")
				inconsist_mic_list_file.write("# The information about inconsistent  IDs\n")
				# Loop over substring id list
				for mic_id_substr in global_entry_dict:
					mic_id_entry = global_entry_dict[mic_id_substr]
				
					consistency_messages = []
					# Check if associated input  path exists
					if not subkey_input_mic_path in mic_id_entry:
						input_mic_path = mic_pattern.replace("*", mic_id_substr)
						consistency_messages.append("    associated input  %s." % (input_mic_path))
				
					# Check if associated selected  basename exists
					if not subkey_selected_mic_basename in mic_id_entry:
						input_mic_path = mic_pattern.replace("*", mic_id_substr)
						consistency_messages.append("    associated selected  %s." % (input_mic_path))
				
					if len(consistency_messages) > 0:
						inconsist_mic_list_file.write("Micrograph ID %s does not have:\n" % (mic_id_substr))
						for consistency_message in consistency_messages:
							inconsist_mic_list_file.write(consistency_message)
							inconsist_mic_list_file.write("\n")
			
				# Close the consistency check file, if necessary
				inconsist_mic_list_file.flush()
				inconsist_mic_list_file.close()
			
			# Since mic_id_substr is once stored as the key of global_entry_dict and extracted with the key order
			# we need sort the valid_mic_id_substr_list here
			if debug_mode: print(("BEFORE SORT: valid_mic_id_substr_list := ", valid_mic_id_substr_list))
			valid_mic_id_substr_list.sort(key=str.lower) # Sort list of  IDs using case insensitive string comparison
			if debug_mode: print(("AFTER SORT: valid_mic_id_substr_list := ", valid_mic_id_substr_list))
			
			# --------------------------------------------------------------------------------
			# Print out the summary of input consistency
			# --------------------------------------------------------------------------------
			print(" ")
			print("Summary of dataset consistency check...")
			print(("  Detected  IDs               : %6d" % (len(global_entry_dict))))
			print(("  Entries in input directory  : %6d" % (len(input_mic_path_list))))
			print(("  Entries in selection list   : %6d" % (len(selected_mic_path_list))))
			print(("  Rejected by no input        : %6d" % (len(no_input_mic_id_substr_list))))
			print(("  Valid Entries               : %6d" % (len(valid_mic_id_substr_list))))
			
			# --------------------------------------------------------------------------------
			# Check MPI error condition
			# --------------------------------------------------------------------------------
			if len(valid_mic_id_substr_list) < n_mpi_procs:
				error_status = ("Number of MPI processes (%d) supplied by --np in mpirun cannot be greater than %d (number of valid s that satisfy all criteria to be processed). Run %s -h for help." % (n_mpi_procs, len(valid_mic_id_substr_list, program_name)), getframeinfo(currentframe()))
				break
			
			# --------------------------------------------------------------------------------
			# Create input file path list
			# --------------------------------------------------------------------------------
			for mic_id_substr in valid_mic_id_substr_list:
				mic_path = global_entry_dict[mic_id_substr][subkey_input_mic_path]
				input_file_path_list.append(mic_path)
			
			# --------------------------------------------------------------------------------
			# Clean up variables related to tracking of invalid (rejected) s 
			# --------------------------------------------------------------------------------
			del input_mic_path_list
			del selected_mic_path_list
			del no_input_mic_id_substr_list
			
			break
		
		
		# --------------------------------------------------------------------------------
		# Clean up the variables for all sections in  mode case 
		# --------------------------------------------------------------------------------
		del mic_pattern
		del mic_basename_pattern
		del global_entry_dict
		del subkey_input_mic_path
		del subkey_selected_mic_basename
		del valid_mic_id_substr_list
		
		# --------------------------------------------------------------------------------
		# Print all error messages and abort the process if necessary.
		# --------------------------------------------------------------------------------
		# NOTE: Toshio Moriya 2016/11/15
		# The following function takes care of the case when an if-statement uses break for occurence of an error.
		# However, more elegant way is to use 'exception' statement of exception mechanism...
		# 
		if_error_then_all_processes_exit_program(error_status)
		
	else:
		input_file_path_list.append(input_image_path)
	
	if RUNNING_UNDER_MPI:
		# Wait for main mpi process to create the input file path list
		mpi_barrier(MPI_COMM_WORLD)
		
		# All mpi processes should know input file path list
		input_file_path_list = wrap_mpi_bcast(input_file_path_list, main_mpi_proc)
	
	# ====================================================================================
	# Prepare input file path(s)
	# ====================================================================================
	# 
	# NOTE: 2016/03/17 Toshio Moriya
	# From here on, stack (and namics) will be used to distinguish stack mode and  mode.
	# However, a single input_file_path_list should be sufficient since we already know the mode.
	# Let's consider this refactoring in the future.
	# 
	stack = None # (particle) stack file name: if it is not None, cter runs with stack mode. Otherwise, runs with  mode
	namics = []  #  file name list
	if not stack_mode:
		namics = input_file_path_list
		if debug_mode: print(("BEFORE SORT: namics := ", namics))
		namics.sort(key=str.lower) # Sort list of s using case insensitive string comparison
		if debug_mode: print(("AFTER SORT: namics := ", namics))
	else:
		stack = input_file_path_list[0]
	
	del input_file_path_list # Don't need this anymore
	
	# Make output directory
	outpwrot = "%s/pwrot" % (output_directory)
	if stack == None: 
		outmicthumb = "%s/micthumb" % (output_directory)
	if debug_mode:  
		outravg = "%s/ravg" % (output_directory)
	if my_mpi_proc_id == main_mpi_proc:
		# Make output directory
		if not os.path.exists(output_directory):
			os.mkdir(output_directory)
		os.mkdir(outpwrot)
		if stack == None: 
			os.mkdir(outmicthumb)
		if debug_mode:
			os.mkdir(outravg)
	
	if RUNNING_UNDER_MPI:
		# Make all mpi processes wait for main mpi process to create output directory
		mpi_barrier(MPI_COMM_WORLD)
	
	# Set up loop variables depending on the cter mode
	if stack == None:
		if RUNNING_UNDER_MPI:
			set_start, set_end = MPI_start_end(len(namics), n_mpi_procs, my_mpi_proc_id)
		else:
			set_start = 0
			set_end = len(namics)
	else:
		pw2 = []
		data = EMData.read_images(stack)
		nima = len(data)
		for i in range(nima):
			pw2.append(periodogram(data[i]))
		wn = pw2[0].get_xsize()
		set_start = 0
		set_end = 1

	# Set up progress message
	if my_mpi_proc_id == main_mpi_proc:
		print(" ")
		print("Estimating CTF parameters...")
		if stack == None:
			print("  Micrographs processed by main process (including percent of progress):")
			progress_percent_step = (set_end - set_start)/100.0 # the number of micrograms for main mpi processer divided by 100

	totresi = []
	missing_img_names = []
	rejected_img_names = []
	for ifi in range(set_start, set_end):
		img_type = ""
		img_name = ""
		img_basename_root = ""
		
		if stack == None:
			img_type = "Micrograph"
			img_name = namics[ifi]
			
			if my_mpi_proc_id == main_mpi_proc:
				print(("    Processing %s ---> %6.2f%%" % (img_name, (ifi - set_start) / progress_percent_step)))
			
			if not os.path.exists(img_name):
				missing_img_names.append(img_name)
				print("    %s %s: Can not find this file. Skipping the estimation and CTF parameters are not stored..." % (img_type, img_name))
				continue

			mic = get_im(img_name)
			try:
				pw2 = tilemic(mic, win_size = wn, overlp_x = overlap_x, overlp_y = overlap_y, edge_x = edge_x, edge_y = edge_y)
			except:
				print("MRK_DEBUG: tilemic() in cter_pap() raised an exception. The micrographs {} might have a problem. Please check it and remove it if necessary.".format(img_name))
				raise
			del mic
		else:
			img_type = "Stack"
			img_name = stack
			
			numFM = EMUtil.get_image_count(img_name)
			pw2 = [None]*numFM
			for i in range(numFM):
				pw2.append(periodogram(get_im(img_name,i)))

		if debug_mode: print("    %s %s: Process %04d started the processing. Detected %d image(s) in this %s file." % (img_type, img_name, ifi, numFM, img_type.lower()))
		
		if db_check_dict(img_name) == False:
			img_basename_root = os.path.splitext(os.path.basename(img_name))[0]
		else:
			path, dictname, keys = db_parse_path(img_name)
			img_basename_root = dictname

		nimi = len(pw2)
		adefocus = [0.0] * kboot
		aamplitu = [0.0] * kboot
		aangle   = [0.0] * kboot
		
		allroo = []
		for imi in range(nimi):
			allroo.append(rot_avg_table(pw2[imi]))
		lenroo = len(allroo[0])
		#print time(),nimi
		
		for nboot in range(kboot):
			#at = time()
			if(nboot == 0): boot = list(range(nimi))
			else:
				from random import randint
				for imi in range(nimi): boot[imi] = randint(0, nimi - 1)
			qa = model_blank(wn, wn)
			roo  = np.zeros(lenroo, np.float32)
			sroo = np.zeros(lenroo, np.float32)
			aroo = np.zeros(lenroo, np.float32)

			for imi in range(nimi):
				Util.add_img(qa, pw2[boot[imi]])
				temp1 = np.array(allroo[boot[imi]])
				roo += temp1
				temp2 = movingaverage(temp1, 10)
				aroo += temp2
				sroo += temp2**2
			sroo[0] = sroo[1]
			aroo[0] = aroo[1]
			sroo = (sroo-aroo**2 / nimi) / nimi
			aroo /= nimi
			roo  /= nimi
			qa   /= nimi

			if f_start < 0:
				#  Find a break point
				bp = 1.e23
				for i in range(5, lenroo - 5):
					#t1 = linreg(sroo[:i])
					#t2 = linreg(sroo[i:])
					#tt = t1[1][0] + t2[1][0]
					xtt = np.array(list(range(i)), np.float32)
					zet = np.poly1d(np.polyfit(xtt,sroo[:i], 2))
					t1 = sum((sroo[:i] - zet(xtt))**2)
					xtt = np.array(list(range(i, lenroo)), np.float32)
					zet = np.poly1d(np.polyfit(xtt, sroo[i:], 2) )
					tt = t1 + sum((sroo[i:] - zet(xtt))**2)
					if tt < bp:
						bp = tt
						istart = i
				#istart = 25
				#print istart
				f_start = istart / (pixel_size * wn)
			"""
			hi = hist_list(sroo,2)
			# hi[0][1] is the threshold
			for i in xrange(1,len(sroo)):
				if(sroo[i] < hi[0][1]):
					istart = i
					break
			"""
			#write_text_file([roo.tolist(),aroo.tolist(),sroo.tolist()], "sroo%03d.txt"%ifi)
			rooc = roo.tolist()

			#print namics[ifi],istart,f_start
			
			defc, subpw, ctf2, baseline, envelope, istart, istop = defocusgett(rooc, wn, voltage = voltage, Pixel_size = pixel_size, Cs = Cs, ampcont = wgh, f_start = f_start, f_stop = f_stop, round_off = 1.0, nr1 = 3, nr2 = 6, parent = None, DEBug = debug_mode)
			if debug_mode:
				print("  RESULT %s" % (img_name), defc, istart, istop)
				
				freq = list(range(len(subpw)))
				for i in range(len(freq)):  freq[i] = float(i) / wn / pixel_size
#				write_text_file([freq, subpw.tolist(), ctf2, envelope.tolist(), baseline.tolist()], "%s/ravg%05d.txt" % (output_directory, ifi))
				fou = os.path.join(outravg, "%s_ravg_%02d.txt" % (img_basename_root, nboot))
				write_text_file([freq, subpw.tolist(), ctf2, envelope.tolist(), baseline.tolist()], fou)
			#mpi_barrier(MPI_COMM_WORLD)
			
			#exit()
			bg = baseline.tolist()
			en = envelope.tolist()
			
			bckg = model_blank(wn, wn, 1, 1)
			envl = model_blank(wn, wn, 1, 1)
			
			from math import sqrt
			nc = wn // 2
			bg.append(bg[-1])
			en.append(en[-1])
			for i in range(wn):
				for j in range(wn):
					r = sqrt((i - nc)**2 + (j - nc)**2)
					ir = int(r)
					if(ir < nc):
						dr = r - ir
						bckg.set_value_at(i, j, (1. - dr) * bg[ir] + dr * bg[ir + 1] )
						envl.set_value_at(i, j, (1. - dr) * en[ir] + dr * en[ir + 1] )
			
			#qa.write_image("rs1.hdf")
			
			mask = model_circle(istop - 1, wn, wn) * (model_blank(wn, wn, 1, 1.0) - model_circle(istart, wn, wn))
			qse = threshold((qa - bckg))#*envl
			#(qse*mask).write_image("rs2.hdf")
			#qse.write_image("rs3.hdf")
			##  SIMULATION
			#bang = 0.7
			#qse = ctf2_rimg(wn, generate_ctf([defc,Cs,voltage,pixel_size,0.0,wgh, bang, 37.0]) )
			#qse.write_image("rs3.hdf")
			
			cnx = wn // 2 + 1
			cny = cnx
			mode = "H"
			istop = min(wn // 2 - 2, istop)    #2-26-2015@ming
			numr = Numrinit(istart, istop, 1, mode)
			wr = ringwe(numr, mode)
			
			crefim = Util.Polar2Dm(qse*mask, cnx, cny, numr, mode)
			Util.Frngs(crefim, numr)
			Util.Applyws(crefim, numr, wr)
			
			#pc = ctf2_rimg(wn,generate_ctf([defc,Cs,voltage,pixel_size,0.0,wgh]))
			#print ccc(pc*envl, subpw, mask)
			
			bang = 0.0
			bamp = 0.0
			bdef = defc
			bold = 1.e23
			while(True):
				#  in simctf2 data[3] is astigmatism amplitude
				"""
				data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
				#astdata = [crefim, numr, wn, bdef, Cs, voltage, pixel_size, wgh, bang]
				for qqq in xrange(200):
					qbdef = 1.0 + qqq*0.001
					print " VALUE AT THE BEGGINING OF while LOOP  ",qbdef,simctf2(qbdef, data)#,fastigmatism3(bamp,astdata)
				"""
				"""
				bamp = 0.7
				bang = 37.0
				
				data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
				astdata = [crefim, numr, wn, bdef, Cs, voltage, Pixel_size, wgh, bang]
				print " VALUE AT THE BEGGINING OF while LOOP  ",bdef,bamp,bang,simctf2(bdef, data),fastigmatism3(bamp,astdata,mask)
				#print  simctf2out(1.568,data)
				#exit()
				
				for kdef in xrange(14000,17000,10):
					dz = kdef/10000.0
					ard = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
					#print ard
					aqd = [crefim, numr, wn, dz, Cs, voltage, pixel_size, wgh, bang]
					#print aqd
					print  dz,simctf2(dz,ard),fastigmatism3(bamp,aqd,mask)
					#print aqd[-1]
				exit()
				"""
				data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
				h = 0.05 * bdef
				amp1, amp2 = bracket_def(simctf2, data, bdef * 0.9, h)
				#print "bracketing of the defocus  ",amp1, amp2
				#print " ttt ",time()-srtt
				#print "bracketing of the defocus  ",amp1,amp2,simctf2(amp1, data),simctf2(amp2, data),h
				amp1, val2 = goldsearch_astigmatism(simctf2, data, amp1, amp2, tol = 1.0e-3)
				#print "golden defocus ",amp1, val2,simctf2(amp1, data)
				#bdef, bcc = goldsearch_astigmatism(simctf2, data, amp1, amp2, tol=1.0e-3)
				#print "correction of the defocus  ",bdef,bcc
				#print " ttt ",time()-srtt
				"""
				crot2 = rotavg_ctf(ctf2_rimg(wn,generate_ctf([bdef, Cs, voltage, pixel_size, 0.0, wgh, bamp, bang])), bdef, Cs, voltage, pixel_size, bamp, bang)
				pwrot = rotavg_ctf(qa-bckg, bdef, Cs, voltage, pixel_size, bamp, bang)
				write_text_file([range(len(subroo)),asubroo, ssubroo, sen, pwrot, crot2],"rotinf%04d.txt"%ifi)
				qse.write_image("qse.hdf")
				mask.write_image("mask.hdf")
				exit()
				"""
				
				astdata = [crefim, numr, wn, bdef, Cs, voltage, pixel_size, wgh, bang, mask]
				h = 0.01
				amp1, amp2 = bracket(fastigmatism3, astdata, h)
				#print "  astigmatism bracket  ",amp1,amp2,astdata[-1]
				#print " ttt ",time()-srtt
				bamp, bcc = goldsearch_astigmatism(fastigmatism3, astdata, amp1, amp2, tol = 1.0e-3)
				junk = fastigmatism3(bamp,astdata)
				bang = astdata[8]
				
				#print astdata[8]
				#print  fastigmatism3(0.0,astdata)
				#print astdata[8]
				#temp = 0.0
				#print bdef, Cs, voltage, pixel_size, temp, wgh, bamp, bang, -bcc
				#data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
				#astdata = [crefim, numr, wn, bdef, Cs, voltage, pixel_size, wgh, bang]
				#print " VALUE WITHIN the while LOOP  ",bdef,bamp,bang,simctf2(bdef, data),fastigmatism3(bamp,astdata)
				#print "  golden search ",bamp,data[-1], fastigmatism3(bamp,data), fastigmatism3(0.0,data)
				#print " ttt ",time()-srtt
				#bamp = 0.5
				#bang = 277
				
				dama = amoeba([bdef, bamp], [0.2, 0.2], fupw, 1.e-4, 1.e-4, 500, astdata)
				if debug_mode:  print("AMOEBA    ", dama)
				bdef = dama[0][0]
				bamp = dama[0][1]
				astdata = [crefim, numr, wn, bdef, Cs, voltage, pixel_size, wgh, bang, mask]
				junk = fastigmatism3(bamp, astdata)
				bang = astdata[8]
				if debug_mode:  print(" after amoeba ", bdef, bamp, bang)
				#  The looping here is blocked as one shot at amoeba is good enough.  To unlock it, remove - from bold.
				if(bcc < -bold): bold = bcc
				else:           break
			
			#data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
			#print " VALUE AFTER the while LOOP  ",bdef,bamp,bang,simctf2(bdef, data),fastigmatism3(bamp,astdata)
			#temp = 0.0
			#print ifi,bdef, Cs, voltage, pixel_size, temp, wgh, bamp, bang, -bcc
			#freq = range(len(subpw))
			#for i in xrange(len(freq)):  freq[i] = float(i)/wn/pixel_size
			#ctf2 = ctf_2(wn, generate_ctf([bdef,Cs,voltage,pixel_size,0.0,wgh]))[:len(freq)]
			#write_text_file([freq, subpw.tolist(), ctf2, envelope.tolist(), baseline.tolist()],"ravg/ravg%05d.txt"%ifi)
			#print " >>>> ",wn, bdef, bamp, Cs, voltage, pixel_size, wgh, bang
			#data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
			#print  simctf2out(bdef, data)
			#exit()
			adefocus[nboot] = bdef
			aamplitu[nboot] = bamp
			aangle[nboot]   = bang
			#from sys import exit
			#exit()
		
		#print " ttt ",time()-srtt
		#from sys import exit
		#exit()
		ad1, ad2, ad3, ad4 = table_stat(adefocus) # return values: average, variance, minimum, maximum
		if ad2 <= 0.0:
			print(("    %s %s: Detected the variance less than zero (defocus statics: avg = %f, var = %f, min = %f, max = %f)." % (img_type, img_name, ad1, ad2, ad3, ad4)))
			print(("           The program ignores this estimate..."))
			continue
		
		reject = []
		thr = 3 * sqrt(ad2)
		for i in range(len(adefocus)):
			if(abs(adefocus[i] - ad1) > thr):
				print(("    %s %s: Rejected an outlier defocus estimate (defocus = %f, average defocus = %f, threshold = %f)." % (img_type, img_name, adefocus[i], ad1, thr)))
				reject.append(i)
		
		if(len(reject) > 0):
			print(("    %s %s: Total number of rejects %s" % (img_type, img_name, len(reject))))
			for i in range(len(reject) - 1, -1, -1):
				del adefocus[i]
				del aamplitu[i]
				del aangle[i]
		
		if(len(adefocus) < 2):
			print(("    %s %s: After rejection of outliers, there is too few estimated defocus values. Skipping the estimation and CTF parameters are not stored..." % (img_type, img_name)))
		else:
			#print "adefocus",adefocus
			#print  "aamplitu",aamplitu
			#print "aangle",aangle
			ad1, ad2, ad3, ad4 = table_stat(adefocus)
			bd1, bd2, bd3, bd4 = table_stat(aamplitu)
			cd1,cd2 = angle_ave([2*q for q in aangle])  # Have to use this trick as the function works for range [0,360]
			cd1/=2
			cd2/=2
			temp = 0.0
			stdavad1 = np.sqrt(kboot * max(0.0, ad2))
			stdavbd1 = np.sqrt(kboot * max(0.0, bd2))
			cd2 *= np.sqrt(kboot)
			
			# Adjust value ranges of astig. amp. and angle.
			if( bd1 < 0.0 ):
				bd1 = -bd1
				cd1 = 90.0 + cd1
			cd1 = cd1 % 180
			
			#  SANITY CHECK, do not produce anything if defocus abd astigmatism amplitude are out of whack
			reject_img_messages = []
			try:
				pwrot2 = rotavg_ctf( model_blank(wn, wn), ad1, Cs, voltage, pixel_size, bd1, cd1)
			except:
				reject_img_messages.append("    - Astigmatism amplitude (%f) is larger than defocus (%f) or defocus (%f) is negative." % (bd1, ad1, ad1))
			
			if len(reject_img_messages) > 0:
				rejected_img_names.append(img_name)
				print("    %s %s: Rejected the CTF estimate - " % (img_type, img_name), ad1, Cs, voltage, pixel_size, wgh, bd1, cd1, "(def, Cs, vol, apix, amp_contrast, astig_amp, astig_angle)")
				print("    %s %s: because... " % (img_type, img_name))
				for reject_img_message in reject_img_messages:
					print(reject_img_message)
				print("    %s %s: Skipping the estimation and CTF parameters are not stored..." % (img_type, img_name))
			else: # assert(len(img_reject_messages) == 0)
				#  Estimate the point at which (sum_errordz ctf_1(dz+errordz))^2 falls to 0.5
				import random as rqt
				
				supe = model_blank(wn, wn)
				niter = 1000
				for it in range(niter):
					Util.add_img(supe, Util.ctf_rimg(wn, wn, 1, ad1 + rqt.gauss(0.0,stdavad1), pixel_size, voltage, Cs, 0.0, wgh, bd1 + rqt.gauss(0.0,stdavbd1), cd1 + rqt.gauss(0.0,cd2), 1))
				ni = wn // 2
				supe /= niter
				pwrot2 = rotavg_ctf(supe, ad1, Cs, voltage, pixel_size, bd1, cd1)
				for i in range(ni):  pwrot2[i] = pwrot2[i]**2
				
				ibec = 0
				for it in range(ni - 1, 0, -1):
					if pwrot2[it] > 0.5 :
						ibec = it
						break
				from morphology import ctf_1d
				ct = generate_ctf([ad1, Cs, voltage, pixel_size, temp, wgh, 0.0, 0.0])
				cq = ctf_1d(wn, ct)
				
				supe = [0.0] * ni
				niter = 1000
				for i in range(niter):
					cq = generate_ctf([ad1 + rqt.gauss(0.0,stdavad1), Cs, voltage, pixel_size, 0.0, wgh, 0.0, 0.0])
					ci = ctf_1d(wn, cq)[:ni]
					for l in range(ni):  supe[l] +=ci[l]
				
				for l in range(ni):  supe[l] = (supe[l] / niter)**2
				
				ib1 = 0
				for it in range(ni - 1, 0, -1):
					if supe[it] > 0.5:
						ib1 = it
						break
				ibec = ibec / (pixel_size * wn)  #  with astigmatism
				ib1  = ib1 / (pixel_size * wn)   #  no astigmatism
				#from utilities import write_text_file
				#write_text_file([range(ni), supe[:ni],pwrot2[:ni]],"fifi.txt")
				
				# Compute defocus CV and astig. amp. CV (CV: coefficient of variation; ratio of error (SD) relative to average (mean))
				#if ad1 < max(0.0, valid_min_defocus): ERROR("Logical Error: Encountered unexpected defocus value (%f). Consult with the developer." % (ad1), "%s in %s" % (__name__, os.path.basename(__file__))) # MRK_ASSERT
				#  TOSHIO - no need to ceck, it was computed as sqrt above, so it cannot be <0
				#if stdavad1 < 0.0: ERROR("Logical Error: Encountered unexpected defocus SD value (%f). Consult with the developer." % (stdavad1), "%s in %s" % (__name__, os.path.basename(__file__))) # MRK_ASSERT
				cvavad1 = stdavad1 / ad1 * 100 # use percentage
				
				if bd1 < 0.0: ERROR("Logical Error: Encountered unexpected astig. amp. value (%f). Consult with the developer." % (bd1), "%s in %s" % (__name__, os.path.basename(__file__))) # MRK_ASSERT

				bd1 = max(bd1, 1.0e-15)
				
				cvavbd1 = stdavbd1 / bd1 * 100 # use percentage
				
				# Compute CTF limit (theoretical resolution limit based on the oscillations of CTF) 
				# For output, use ctflim (relative frequency limit [1/A]), not ctflim_abs (absolute frequency limit)
				# 
				# NOTE: 2016/03/23 Toshio Moriya
				# xr is limiting frequency [1/A]. Max is Nyquist frequency = 1.0/(2*apix[A/pixel]). <UNIT: [1/(A/pixel)/[pixel])] => [(pixel)/(A*pixel] => [1/A]>
				# 1.0/xr is limiting period (Angstrom resolution) [A]. Min is Nyquist period = (2*apix[A/pixel]). <UNIT: [1/(1/A)] = [A]>
				# fwpix is width of Fourier pixel [pixel/A] := 1.0[pixel]/(2*apix[A/pixel])/box_half[pixel] = 1[pixel]/fullsize[A]). <UNIT: [pixel/(A/pixel)/(pixel)] = [pixel*(pixel/A)*(1/pixel) = [pixel/A]>
				# int(xr/fwpix+0.5) is limiting_absolute_frequency [1/pixel]. <Unit:[(1/A)/(pixel/A)] = [(1/A)*(A/pixel)] = [1/pixel]>
				# return  int(xr/fwpix+0.5),xr, which is limiting_abs_frequency [1/pixel], and Limiting_frequency[1/A]
				#
				ctflim_abs, ctflim = ctflimit(wn, ad1, Cs, voltage, pixel_size)
				
				"""
				for i in xrange(len(ssubroo)):
					asubroo[i] /= kboot
					ssubroo[i]  = sqrt(max(0.0, ssubroo[i]-kboot*asubroo[i]**2)/kboot)
					sen[i]     /= kboot
				"""
				lnsb = len(subpw)
				try:		crot1 = rotavg_ctf(ctf_rimg(wn, generate_ctf([ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1]), sign=0), ad1, Cs, voltage, pixel_size, 0.0, 0.0)[:lnsb]
				except:		crot1 = [0.0] * lnsb
				try:		pwrot1 = rotavg_ctf(threshold(qa - bckg), ad1, Cs, voltage, pixel_size, 0.0, 0.0)[:lnsb]
				except:		pwrot1 = [0.0] * lnsb
				try:		crot2 = rotavg_ctf(ctf_rimg(wn, generate_ctf([ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1]), sign=0), ad1, Cs, voltage, pixel_size, bd1, cd1)[:lnsb]
				except:		crot2 = [0.0] * lnsb
				try:		pwrot2 = rotavg_ctf(threshold(qa - bckg), ad1, Cs, voltage, pixel_size, bd1, cd1)[:lnsb]
				except:		pwrot2 = [0.0] * lnsb
				#  #1 - rotational averages without astigmatism, #2 - with astigmatism
				lnsb = min(lnsb,len(crot2),len(pwrot1),len(crot2),len(pwrot2))
				write_text_file([list(range(lnsb)), [float(i)/wn/pixel_size for i in range(lnsb)], pwrot1[:lnsb], crot1[:lnsb], pwrot2[:lnsb], crot2[:lnsb]], os.path.join(outpwrot, "%s_rotinf.txt"%(img_basename_root)))
				#
				# NOTE: 2016/03/23 Toshio Moriya
				# Compute mean of extrema differences (differences at peak & trough) between 
				# (1) experimental rotational average with astigmatism (pwrot2)
				# (2) experimental rotational average without astigmatism (pwrot1), and
				# as a indication of goodness of astigmatism estimation by cter.
				# The peak & trough detection uses fitted rotational average with astigmatism (crot2) 
				# Start from 1st trough while ignoring 1st peak.
				# End at astigmatism frequency limit.
				# 
				is_peak_target = True
				pre_crot2_val = crot2[0]
				extremum_counts = 0
				extremum_diff_sum = 0
				for i in range(1, len(crot2)):
					cur_crot2_val = crot2[i]
					if is_peak_target == True and pre_crot2_val > cur_crot2_val:
						# peak search state
						extremum_i = i - 1
						extremum_counts += 1
						extremum_diff_sum += pwrot2[extremum_i] - pwrot1[extremum_i] # This should be positive if astigmatism estimation is good
						# print "MRK_DEBUG: Peak Search  : extremum_i = %03d, freq[extremum_i] = %12.5g, extremum_counts = %03d, (pwrot2[extremum_i] - pwrot1[extremum_i]) = %12.5g, extremum_diff_sum = %12.5g " % (extremum_i, freq[extremum_i] , extremum_counts, (pwrot2[extremum_i] - pwrot1[extremum_i]), extremum_diff_sum)
						is_peak_target = False
					elif is_peak_target == False and pre_crot2_val < cur_crot2_val:
						# trough search state
						extremum_i = i - 1
						extremum_counts += 1
						extremum_diff_sum += pwrot1[extremum_i] - pwrot2[extremum_i] # This should be positive if astigmatism estimation is good
						# print "MRK_DEBUG: Trough Search: extremum_i = %03d, freq[extremum_i] = %12.5g, extremum_counts = %03d, (pwrot1[extremum_i] - pwrot2[extremum_i]) = %12.5g, extremum_diff_sum = %12.5g " % (extremum_i, freq[extremum_i] , extremum_counts, (pwrot1[extremum_i] - pwrot2[extremum_i]), extremum_diff_sum)
						is_peak_target = True
					pre_crot2_val = cur_crot2_val
#				#if extremum_counts == 0: ERROR("Logical Error: Encountered unexpected zero extremum counts. Consult with the developer." % (bd1), "%s in %s" % (__name__, os.path.basename(__file__))) # MRK_ASSERT
#				extremum_diff_avg = 1.1#extremum_diff_sum / extremum_counts
				
				# print "MRK_DEBUG: extremum_avg = %12.5g, extremum_diff_sum = %12.5g, extremum_counts = %03d," % (extremum_avg, extremum_diff_sum, extremum_counts)
				
#				if stack == None:     cmd = "echo " + "    " + namics[ifi] + "  >>  " + fou
#				else:                 cmd = "echo " + "    " + "  >>  " + fou
#				os.system(cmd)
				
				ed1      = wgh            # total amplitude contrast. Here, it is also constant amplitude contrast since Volta phase shift is not estimated with this function. 
				stdaved1 = 0.0            # dummy value for error of total amplitude contrast estimation
				max_freq = 0.5/pixel_size # dummy value for maximum frequency. set to Nyquist frequency for now. let's add the implementation in near future (Toshio 2017/12/06)
				reserved = 0.0            # dummy value for reserved spot, which might be used for parameter of external programs (e.g. CTFFIND4, GCTF, and etc.)
				# wgh                     # constant amplitude contrast provided by user (default 10%). Here, it is also total amplitude contrast since Volta phase shift is not estimated with this function.
				phase_shift = ampcont2angle(ed1) - ampcont2angle(wgh) # Volta phase shift [deg] = total amplitude contrast phase shift [deg] (ed1) -  constant amplitude contrast phase shift [deg]; ed1 is boot strap average of total amplitude contrast [%]
				
				if debug_mode: print(("    %s %s: Process %04d finished the processing. Estimated CTF parmaters are stored in %s." % (img_type, img_name, ifi, os.path.join(output_directory, "partres.txt"))))
#				if debug_mode: print((ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, cvavad1, cvavbd1, extremum_diff_avg, ib1, ibec, ctflim))
#				# totresi.append([img_name, ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, cvavad1, cvavbd1, extremum_diff_avg, ib1, ibec, ctflim])
#				totresi.append([img_name, ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdaved1, stdavbd1, cd2, cvavad1, cvavbd1, extremum_diff_avg, ib1, ibec, ctflim])
				if debug_mode: print((ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1, stdavad1, stdavbd1, stdavbd1, cd2, cvavad1, cvavbd1, ib1, ibec, ctflim, max_freq, reserved, wgh, phase_shift))
				totresi.append([img_name, ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1, stdavad1, stdaved1, stdavbd1, cd2, cvavad1, cvavbd1, ib1, ibec, ctflim, max_freq, reserved, wgh, phase_shift])
				
#				if stack == None:
#					print  namics[ifi], ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec
#				else:
#					print               ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec
#				if stack == None:
#					totresi.append( [ namics[ifi], ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec])
#				else:
#					totresi.append( [ 0, ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec])
#				#if ifi == 4 : break

		if stack == None:
			img_mic = get_im(namics[ifi])
			# create  thumbnail
			nx = img_mic.get_xsize()
			if nx > 512:
				img_micthumb = resample(img_mic, 512.0/nx)
			else:
				img_micthumb = img_mic
			fou = os.path.join(outmicthumb, "%s_thumb.hdf" % (img_basename_root))
			img_micthumb.write_image(fou)

	if RUNNING_UNDER_MPI:
		from utilities import wrap_mpi_gatherv
		totresi = wrap_mpi_gatherv(totresi, 0, MPI_COMM_WORLD)
		missing_img_names = wrap_mpi_gatherv(missing_img_names, 0, MPI_COMM_WORLD)
		rejected_img_names = wrap_mpi_gatherv(rejected_img_names, 0, MPI_COMM_WORLD)
	
	if my_mpi_proc_id == main_mpi_proc:
		outf = open(os.path.join(output_directory, "partres.txt"), "w")
		for i in range(len(totresi)):
			for k in range(1, len(totresi[i])):
				outf.write("  %12.5g" % totresi[i][k])
			outf.write("  %s\n" % totresi[i][0])
		outf.close()
		
		print(" ")
		print(("Summary of %s processing..." % (img_type.lower())))
		missing_counts = len(missing_img_names)
		print(("  Missing  : %d" % (missing_counts)))
		if missing_counts > 0:
			outfile_path = os.path.join(output_directory, "missing_%s_list.txt" % (img_type.lower()))
			print(("    Saving list of missing in %s..." % (outfile_path)))
			outf = open(outfile_path, "w")
			for missing_img_name in missing_img_names:
				outf.write("%s\n" % missing_img_name)
			outf.close()
		
		rejected_counts = len(rejected_img_names)
		print(("  Rejected : %d" % (rejected_counts)))
		if rejected_counts > 0:
			outfile_path = os.path.join(output_directory, "rejected_%s_list.txt" % (img_type.lower()))
			print(("    Saving list of rejected in %s..." % (outfile_path)))
			outf = open(outfile_path, "w")
			for rejected_img_name in rejected_img_names:
				outf.write("%s\n" % rejected_img_name)
			outf.close()
	
	if cter_mode_idx == idx_cter_mode_stack:
		return totresi[0][1], totresi[0][7], totresi[0][8], totresi[0][9], totresi[0][10], totresi[0][11]



########################################
# functions used by cter
# Later on make sure these functions don't conflict with those used
# in the cross resolution program getastcrfNOE.py
########################################

def ampcont2angle(A):
	#  convert amplitude contrast to phase shift
	from math import sqrt, atan, degrees
	if(A == 100.0):  return 90.0
	elif(A == -100.0):  return 90.0
	elif(A<0.0):  return degrees(atan(A/sqrt(1.0e4-A**2)))+180.0
	else:  return degrees(atan(A/sqrt(1.0e4-A**2)))

def angle2ampcont(phi):
	#  convert phase shift to amplitude contrast
	from math import sqrt, tan, radians
	return tan(radians(phi))/sqrt(1.0+tan(radians(phi))**2)*100.0

def bracket_original(f, x1, h):
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


 
def bracket_def(f, dat, x1, h):
	c = 1.618033989 
	f1 = f(x1, dat)
	x2 = x1 + h
	f2 = f(x2, dat)
	#print x1,f1,x2,f2
	# Determine downhill direction and change sign of h if needed
	if f2 > f1:
		#print  h
		h = -h
		x2 = x1 + h
		f2 = f(x2, dat)
		# Check if minimum between x1 - h and x1 + h
		if f2 > f1: return x2,x1 - h 
	# Search loop
	for i in range (100):
		h = c*h
		x3 = x2 + h; f3 = f(x3, dat)
		#print i,x1,f1,x2,f2,x3,f3
		if f3 > f2: return x1,x3
		x1 = x2; x2 = x3
		f1 = f2; f2 = f3
	print("Bracket did not find a mimimum")
	return None, x3


def bracket(f, dat, h):
	c = 1.618033989
	x1 = 0.0
	f1 = f(x1, dat)
	x2 = x1 + h
	f2 = f(x2, dat)
	# Search loop
	for i in range (100):    
		h = c*h
		x3 = x2 + h
		f3 = f(x3, dat)
		#print i,x1,f1,x2,f2,x3,f3
		if f3 > f2: return x1,x3
		x1 = x2; x2 = x3
		f1 = f2; f2 = f3
	print("Bracket did not find a mimimum")
 
def goldsearch_astigmatism(f, dat, a, b, tol=1.0e-3):
	from math import log, ceil
	nIter = int(ceil(-2.078087*log(tol/abs(b-a)))) # Eq. (10.4)
	R = 0.618033989
	C = 1.0 - R
	# First telescoping
	x1 = R*a + C*b; x2 = C*a + R*b
	f1 = f(x1, dat); f2 = f(x2, dat)
	# Main loop
	for i in range(nIter):
		if f1 > f2:
			a = x1
			x1 = x2; f1 = f2
			x2 = C*a + R*b; f2 = f(x2, dat)
		else:
			b = x2
			x2 = x1; f2 = f1
			x1 = R*a + C*b; f1 = f(x1, dat)
	if f1 < f2: return x1,f1
	else: return x2,f2

def defocus_baseline_fit(roo, i_start, i_stop, nrank, iswi):
	"""
		    iswi = 2 using polynomial n rank to fit envelope function
			iswi = 3 using polynomial n rank to fit background
			The background fit is done between i_start, i_stop, but the entire baseline curve is evaluated and subtracted
	"""
	import numpy as np
	from morphology import imf_params_cl1
	
	TMP = imf_params_cl1(roo[i_start:i_stop], nrank, iswi)
	nroo = len(roo)
	baseline   = np.zeros(nroo, np.float32)
	ord = len(TMP[-1])
	freqstep = TMP[0][1]
	if(iswi == 3):
		for i in range(len(roo)):
			freq = freqstep*(i-i_start)
			tmp = TMP[-1][-1]
			for j in range(1,ord): tmp += TMP[-1][j-1]*freq**j
			baseline[i] = tmp
	else:
		for i in range(i_start, i_stop):
			freq = freqstep*(i-i_start)
			tmp = TMP[-1][-1]
			for j in range(1,ord): tmp += TMP[-1][j-1]*freq**j
			baseline[i] = tmp
	return np.exp(baseline)

def simpw1d(defocus, data):
	import numpy as np
	from morphology import ctf_2
	from utilities import generate_ctf
	from math import sqrt
	
	#[defocus, cs, voltage, apix, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle]
	#  data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
	# data[1] - envelope
	#ct = data[1]*np.array( ctf_1d(data[2], generate_ctf([defocus, data[4], data[5], data[6], 0.0, data[7], 0.0, 0.0]), doabs= True)[data[8]:data[9]], np.float32)
	ct = data[1]*np.array( ctf_2(data[2], generate_ctf([defocus, data[4], data[5], data[6], 0.0, data[7], 0.0, 0.0]))[data[8]:data[9]], np.float32)
	#print  " 1d  ",sum(data[0]*ct),np.linalg.norm(ct,2)
	return  -sum(data[0]*ct)/np.linalg.norm(ct,2)


def simpw1d_pap(defocus, data):
	import numpy as np
	from morphology import ctf_1d
	from utilities import generate_ctf
	from math import sqrt
	
	#[defocus, cs, voltage, apix, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle]
	#  data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
	# data[1] - envelope
	#ct = data[1]*np.array( ctf_1d(data[2], generate_ctf([defocus, data[4], data[5], data[6], 0.0, data[7], 0.0, 0.0]), doabs= True)[data[8]:data[9]], np.float32)
	ct = np.array( ctf_1d(data[2], generate_ctf([defocus, data[4], data[5], data[6], 0.0, data[7], 0.0, 0.0]), doabs= True)[data[8]:data[9]], np.float32)
	#print  " 1d  ",sum(data[0]*ct),np.linalg.norm(ct,2)
	return  -sum(data[0]*ct/data[1])/np.linalg.norm(ct,2)

def simpw1d_print(defocus, data):
	import numpy as np
	from morphology import ctf_1d
	from utilities import generate_ctf
	from math import sqrt
	
	#[defocus, cs, voltage, apix, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle]
	#  data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
	# data[1] - envelope
	#ct = data[1]*np.array( ctf_1d(data[2], generate_ctf([defocus, data[4], data[5], data[6], 0.0, data[7], 0.0, 0.0]), doabs= True)[data[8]:data[9]], np.float32)
	ct = np.array( ctf_1d(data[2], generate_ctf([defocus, data[4], data[5], data[6], 0.0, data[7], 0.0, 0.0]), doabs= True)[data[8]:data[9]], np.float32)
	#print  " 1d  ",sum(data[0]*ct),np.linalg.norm(ct,2)
	for i in range(len(data[0])):  print(i,i+data[8],data[0][i],ct[i],data[1][i],data[0][i]/data[1][i])
	return  -sum(data[0]*ct/data[1])/np.linalg.norm(ct,2)

def simpw2d(defocus, data2d):
	from utilities import generate_ctf
	from morphology import ctf_rimg
	from math import sqrt
	
	#             0        1     2      3     4         5             6                      7           
	#           [defocus, cs, voltage, apix, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle]
	#
	#             0        1             2      3    4         5           6        7            8                     9            10
	#  data2d = [nx, experimental_pw, defocus, Cs, voltage, Pixel_size, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle, mask]
	
	defocust = max(min(defocus, 6.0), 0.01)
	data2d[7] = max(min(data2d[7],99.0), 1.0)
	ct = ctf_rimg(data2d[0], generate_ctf([defocust, data2d[3], data2d[4], data2d[5], data2d[6], data2d[7], data2d[8], data2d[9]]), sign=0, ny=data2d[0])
	q2 = ct.cmp("dot", ct, dict(negative = 0, mask = data2d[10], normalize = 0))#Util.infomask(ct*ct, data2d[10], True)[0]
	q1 = ct.cmp("dot", data2d[1], dict(negative = 0, mask = data2d[10], normalize = 0))
	'''
	from utilities import info
	print  info(data2d[1], data2d[10])
	print  info(ct, data2d[10])
	print q1,q2
	'''
	return  -q1/q2


def simpw1dc(defocus, data):
	import numpy as np
	from morphology import ctf_2
	from utilities import generate_ctf
	
	#[defocus, cs, voltage, apix, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle]
	#  data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
	# data[1] - envelope
	ct = data[1]*np.array( ctf_2(data[2], generate_ctf([defocus, data[4], data[5], data[6], 0.0, data[7], 0.0, 0.0]))[data[8]:data[9]], np.float32)
	print(" 1d  ",sum(data[0]*ct),np.linalg.norm(ct,2))
	return  2.0-sum(data[0]*ct)/np.linalg.norm(ct,2),ctf_2(data[2], generate_ctf([defocus, data[4], data[5], data[6], 0.0, data[7], 0.0, 0.0]))

def simpw2dc(defocus, data2d):
	from utilities import generate_ctf
	from morphology import ctf2_rimg
	from math import sqrt
	
	#             0        1     2      3     4         5             6                      7           
	#           [defocus, cs, voltage, apix, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle]
	#
	#             0        1             2      3    4         5           6        7            8                     9            10
	#  data2d = [nx, experimental_pw, defocus, Cs, voltage, Pixel_size, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle, mask]
	

	ct = ctf2_rimg(data2d[0], generate_ctf([defocus, data2d[3], data2d[4], data2d[5], data2d[6], data2d[7], data2d[8], data2d[9]]), ny=data2d[0])
	from utilities import info
	q1 = ct.cmp("dot", data2d[1], dict(negative = 0, mask = data2d[10], normalize = 0))
	q2 = sqrt(ct.cmp("dot", ct, dict(negative = 0, mask = data2d[10], normalize = 0)))
	'''
	print  info(data2d[1], data2d[10])
	print  info(ct, data2d[10])
	'''
	print(" 2d  ",q1,q2)
	return  2.0-q1/q2,ct

def movingaverage(data, window_size, skip = 3):
	import numpy as np
	ld = len(data)
	qt = sum(data[skip:skip+4])/3.0
	tt = type(data[0])
	qt = np.concatenate( ( np.array([qt]*(window_size+skip), tt), data[skip:], np.tile(data[-1],(window_size)) ))
	out = np.empty(ld, np.float32)
	nc1 = window_size - window_size//2
	nc2 = window_size + window_size//2 +1
	for i in range(ld):   out[i] = sum(qt[i+nc1:i+nc2])
	return out*np.float32(1.0/window_size)

def localvariance(data, window_size, skip = 3):
	import numpy as np
	ld = len(data)
	qt = sum(data[skip:skip+4])/3.0
	tt = type(data[0])
	qt = np.concatenate( ( np.array([qt]*(window_size+skip), tt), data[skip:], np.tile(data[-1],(window_size)) ))
	out = np.empty(ld, np.float32)
	nc1 = window_size - window_size//2
	nc2 = window_size + window_size//2 +1
	qnorm = np.float32(1.0/window_size)
	for i in range(ld):
		sav = sum(qt[i+nc1:i+nc2])*qnorm
		sdv = sum(qt[i+nc1:i+nc2]**2)
		out[i] = (qt[i] - sav)/np.sqrt(sdv*qnorm - sav*sav)
	out += min(out)
	return out

def defocusgett(roo, nx, voltage=300.0, Pixel_size=1.0, Cs=2.0, ampcont=0.1, f_start=-1.0, f_stop=-1.0, round_off=1.0, nr1=3, nr2=6, parent=None, DEBug=False):
	"""
	
		1. Estimate envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get 
		   defocus which matches the extracted CTF imprints 
	"""
	from utilities  import generate_ctf
	import numpy as np
	from morphology import ctf_2, bracket_def, defocus_baseline_fit, ctflimit, simpw1d, goldsearch_astigmatism

	#print "CTF params:", voltage, Pixel_size, Cs, wgh, f_start, f_stop, round_off, nr1, nr2, parent

	if f_start == 0 : 	    i_start = 0
	else: 			        i_start = int(Pixel_size*nx*f_start+0.5)
	if f_stop <= f_start :
		i_stop  = len(roo)
		adjust_fstop = True
	else:
		i_stop  = min(len(roo), int(Pixel_size*nx*f_stop+0.5))
		adjust_fstop = False

	nroo = len(roo)

	if DEBug:  print("f_start, i_start, f_stop, i_stop:", f_start, i_start, f_stop, i_stop-1)
	#TE  = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr1), 4)
	#baseline = defocus_baseline_fit(roo, i_start, i_stop, int(nr2), 3)

	baseline = defocus_baseline_fit(roo, i_start,nroo, int(nr2), 3)
	subpw = np.array(roo, np.float32) - baseline
	subpw[0] = subpw[1]
	#write_text_file([roo,baseline,subpw],"dbg.txt")
	#print "IN defocusgett  ",np.min(subpw),np.max(subpw)
	for i in range(len(subpw)):  subpw[i] = max(subpw[i],0.0)
	#print "IN defocusgett  ",np.min(subpw),np.max(subpw)
	#envelope = movingaverage(  subpw   , nroo//4, 3)
	envelope = np.array([1.0]*len(subpw), np.float32)
	#write_text_file([roo,baseline,subpw],"dbgt.txt")

	#print "IN defocusgett  ",np.min(subpw),np.max(subpw),np.min(envelope)
	#envelope = np.ones(nroo, np.float32)
	defocus = 0.0
	data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
	#for i in xrange(nroo):
	#	print  i,"   ",roo[i],"   ",baseline[i],"   ",subpw[i],"   ",envelope[i]
	h = 0.1
	#def1, def2 = bracket(simpw1d, data, h)
	#if DEBug:  print "first bracket ",def1, def2,simpw1d(def1, data),simpw1d(def2, data)
	#def1=0.1
	ndefs = 18
	defound = []
	for  idef in range(ndefs):
		def1 = (idef+1)*0.5
		def1, def2 = bracket_def(simpw1d, data, def1, h)
		if DEBug:  print("second bracket ",idef,def1, def2,simpw1d(def1, data),simpw1d(def2, data),h)
		def1, val2 = goldsearch_astigmatism(simpw1d, data, def1, def2, tol=1.0e-3)
		if DEBug:  print("golden ",idef,def1, val2,simpw1d(def1, data))
		if def1>0.0:  defound.append([val2,def1])
	defound.sort()
	del defound[3:]
	def1 = defound[0][1]
	if adjust_fstop:
		from morphology import ctflimit
		newstop, fnewstop = ctflimit(nx, def1, Cs, voltage, Pixel_size)
		if DEBug:  print("newstop  ", int(newstop*0.7), fnewstop*0.7, i_stop, newstop)
		if( newstop != i_stop and (newstop-i_start)>min(10,(i_stop-i_start))):
			i_stop = newstop
			data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
			"""
			def1, def2 = bracket_def(simpw1d, data, def1, h)
			if(def1 > def2):
				temp = def1
				def1 = def2
				def2 = temp
			print "adjusted bracket ",def1, def2,simpw1d(def1, data)
			"""
			h = 0.05
			for idef in range(3):
				def1, def2 = bracket_def(simpw1d, data, defound[idef][1], h)
				if DEBug:  print(" adjusted def ",def1,def2)
				def1, val2 = goldsearch_astigmatism(simpw1d, data, def1, def2, tol=1.0e-3)
				if DEBug:  print("adjusted golden ",def1, val2,simpw1d(def1, data))
				if def1>0.0:  defound[idef] = [val2,def1]
			defound.sort()
			def1 = defound[0][1]
	if DEBug: print(" ultimate defocus",def1,defound)

	#defocus = defocus_guessn(Res_roo, voltage, Cs, Pixel_size, ampcont, i_start, i_stop, 2, round_off)
	#print simpw1d(def1, data),simpw1d(4.372, data)
	"""
	def1 = 0.02
	def2 = 10.
	def1, def2 = goldsearch_astigmatism(simpw1d, data, def1, def2, tol=1.0e-3)
	print "golden ",def1, def2,simpw1d(def1, data)
	"""
	if DEBug:
		qm = 1.e23
		toto = []
		for i in range(1000,100000,5):
			dc = float(i)/10000.0
			qt = simpw1d(dc, data)
			toto.append([dc,qt])
			if(qt<qm):
				qm=qt
				defi = dc
		from utilities import write_text_row
		write_text_row(toto,"toto1.txt")
		print(" >>>>>>>>>  ",defi,simpw1d(defi, data))#,generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont])
		#def1 = defi
	#exit()
	ctf2 = ctf_2(nx, generate_ctf([def1, Cs, voltage, Pixel_size, 0.0, ampcont]))

	return def1, subpw, ctf2, baseline, envelope, i_start, i_stop


def defocusgett_pap(roo, nx, voltage=300.0, Pixel_size=1.0, Cs=2.0, ampcont=0.1, f_start=-1.0, f_stop=-1.0, round_off=1.0, nr1=3, nr2=6, parent=None, DEBug=False):
	"""
	
		1. Estimate envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get 
		   defocus which matches the extracted CTF imprints 
	"""
	from utilities  import generate_ctf
	import numpy as np
	from morphology import ctf_1d, bracket_def, defocus_baseline_fit, ctflimit, simpw1d_pap, goldsearch_astigmatism

	#print "CTF params:", voltage, Pixel_size, Cs, wgh, f_start, f_stop, round_off, nr1, nr2, parent

	if f_start == 0 : 	    i_start = 0
	else: 			        i_start = int(Pixel_size*nx*f_start+0.5)
	if f_stop <= f_start :
		i_stop  = len(roo)
		adjust_fstop = True
	else:
		i_stop  = min(len(roo), int(Pixel_size*nx*f_stop+0.5))
		adjust_fstop = False

	nroo = len(roo)

	if DEBug:  print("f_start, i_start, f_stop, i_stop:", f_start, i_start, f_stop, i_stop-1)
	#TE  = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr1), 4)
	#baseline = defocus_baseline_fit(roo, i_start, i_stop, int(nr2), 3)

	baseline = defocus_baseline_fit(roo, i_start,nroo, int(nr2), 3)
	subpw = np.array(roo, np.float32) - baseline
	subpw[0] = subpw[1]
	#write_text_file([roo,baseline,subpw],"dbg.txt")
	#print "IN defocusgett  ",np.min(subpw),np.max(subpw)
	for i in range(len(subpw)):  subpw[i] = max(subpw[i],0.0)
	#print "IN defocusgett  ",np.min(subpw),np.max(subpw)
	#envelope = movingaverage(  subpw   , nroo//4, 3)
	envelope = np.array([1.0]*len(subpw), np.float32)
	#write_text_file([roo,baseline,subpw],"dbgt.txt")

	#print "IN defocusgett  ",np.min(subpw),np.max(subpw),np.min(envelope)
	#envelope = np.ones(nroo, np.float32)
	defocus = 0.0
	data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
	#for i in xrange(nroo):
	#	print  i,"   ",roo[i],"   ",baseline[i],"   ",subpw[i],"   ",envelope[i]
	h = 0.1
	#def1, def2 = bracket(simpw1d, data, h)
	#if DEBug:  print "first bracket ",def1, def2,simpw1d(def1, data),simpw1d(def2, data)
	#def1=0.1
	ndefs = 18
	defound = []
	for  idef in range(ndefs):
		def1 = (idef+1)*0.5
		def1, def2 = bracket_def(simpw1d_pap, data, def1, h)
		if DEBug:  print("second bracket ",idef,def1, def2,simpw1d(def1, data),simpw1d_pap(def2, data),h)
		def1, val2 = goldsearch_astigmatism(simpw1d_pap, data, def1, def2, tol=1.0e-3)
		if DEBug:  print("golden ",idef,def1, val2,simpw1d_pap(def1, data))
		if def1>0.0:  defound.append([val2,def1])
	defound.sort()
	del defound[3:]
	def1 = defound[0][1]
	if adjust_fstop:
		from morphology import ctflimit
		newstop, fnewstop = ctflimit(nx, def1, Cs, voltage, Pixel_size)
		if DEBug:  print("newstop  ", int(newstop*0.7), fnewstop*0.7, i_stop, newstop)
		if( newstop != i_stop and (newstop-i_start)>min(10,(i_stop-i_start))):
			i_stop = newstop
			data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
			"""
			def1, def2 = bracket_def(simpw1d, data, def1, h)
			if(def1 > def2):
				temp = def1
				def1 = def2
				def2 = temp
			print "adjusted bracket ",def1, def2,simpw1d(def1, data)
			"""
			h = 0.05
			for idef in range(3):
				def1, def2 = bracket_def(simpw1d_pap, data, defound[idef][1], h)
				if DEBug:  print(" adjusted def ",def1,def2)
				def1, val2 = goldsearch_astigmatism(simpw1d_pap, data, def1, def2, tol=1.0e-3)
				if DEBug:  print("adjusted golden ",def1, val2,simpw1d_pap(def1, data))
				if( def1>0.0 ):  defound[idef] = [val2,def1]
			defound.sort()
			def1 = defound[0][1]
	if DEBug: print(" ultimate defocus",def1,defound)

	#defocus = defocus_guessn(Res_roo, voltage, Cs, Pixel_size, ampcont, i_start, i_stop, 2, round_off)
	#print simpw1d(def1, data),simpw1d(4.372, data)
	"""
	def1 = 0.02
	def2 = 10.
	def1, def2 = goldsearch_astigmatism(simpw1d, data, def1, def2, tol=1.0e-3)
	print "golden ",def1, def2,simpw1d(def1, data)
	"""
	if DEBug:
		qm = 1.e23
		toto = []
		for i in range(1000,100000,5):
			dc = float(i)/10000.0
			qt = simpw1d_pap(dc, data)
			toto.append([dc,qt])
			if(qt<qm):
				qm=qt
				defi = dc
		from utilities import write_text_row
		write_text_row(toto,"toto1.txt")
		print(" >>>>>>>>>  ",defi,simpw1d(defi, data))#,generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont])
		#def1 = defi
	#exit()
	ctf2 = ctf_1d(nx, generate_ctf([def1, Cs, voltage, Pixel_size, 0.0, ampcont]), doabs = True)

	return def1, subpw, ctf2, baseline, envelope, i_start, i_stop


def defocus_guessn(roo, volt, Cs, Pixel_size, ampcont, istart, i_stop):
	"""
		Use specified frequencies area (istart-istop)to estimate defocus
		1.  The searching range is limited to dz_low (.1um) ~ dz_high (20 um).
		    The user can modify this limitation accordingly
		2.  changing nloop can speed up the estimation
		3.  defocus_estimation_method = 1  use squared error
		    defocus_estimation_method = 2  use normalized inner product
		Input:
		  Res_roo - background-subtracted Power Spectrum
		  Res_TE  - background-subtracted Envelope
	"""
	
	from math import sqrt
	from utilities import generate_ctf
	from morphology import ctf_2
	
	import numpy as np

	data = np.array(roo,np.float32)

	envelope = movingaverage(data, 60)
	nx  = int(len(roo)*2)
	nn = len(data)
	goal = -1.e23
	for d in range(20000,56000,10):
		dz = d/10000.
		ct = np.array( ctf_2(nx, generate_ctf([dz, Cs, volt, Pixel_size, 0.0, ampcont]))[:nn], np.float32)
		ct = (ct - sum(ct)/nn)*envelope
		g = sum(ct[istart:]*sub[istart:])/sum(ct[istart:])
		#print d,dz,g
		if(g>goal):
			defocus = d
			goal = g
			#print " ****************************************** ", defocus, goal,istart
	#from utilities import write_text_file
	#write_text_file([sub,envelope,ct,temp],"oto.txt")
	ct = np.array( ctf_2(nx, generate_ctf([defocus, Cs, volt, Pixel_size, 0.0, ampcont]))[:nn], np.float32)
	temp = ct
	ct = (ct - sum(ct)/nn)*envelope
	for i in range(nn):  print(sub[i],envelope[i],ct[i],temp[i])
	from sys import exit
	exit()
	#defocus = int( defocus/round_off )*round_off
	return defocus


'''
def defocusgett_(roo, voltage=300.0, Pixel_size=1.0, Cs=2.0, wgh=0.1, f_start=0.0, f_stop=-1.0, round_off=1.0, nr1=3, nr2=6, parent=None):
	"""
		1. Estimate envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get 
		   defocus which matches the extracted CTF imprints 
	"""
	from utilities  import generate_ctf
	from morphology import defocus_env_baseline_fit, defocus_guess, ctf_2, defocus_guessn
	
	#print "CTF params:", voltage, Pixel_size, Cs, wgh, f_start, f_stop, round_off, nr1, nr2, parent

	defocus = defocus_guessn(roo, voltage, Cs, Pixel_size, wgh, f_start)

	nx  = int(len(roo)*2)
	#ctf = ctf_1d(nx, generate_ctf([defocus, Cs, voltage, Pixel_size, 0.0, wgh]))
	ctf2 = ctf_2(nx, generate_ctf([defocus, Cs, voltage, Pixel_size, 0.0, wgh]))
	"""
	if (parent is not None):
		parent.ctf_data=[roo, Res_roo, Res_TE]
		parent.i_start = i_start
		parent.i_stop  = i_stop
		from utilities import write_text_file
		write_text_file([range(len(roo)), roo, Res_roo, Res_TE, ctf], "procpw.txt")
	else:
		from utilities import write_text_file
		write_text_file([range(len(roo)), roo, Res_roo, Res_TE, ctf, TE, Pn1], "procpw.txt")
	"""
	return defocus, [], ctf2, [],[],0,0
'''


def defocusget_from_crf(roo, voltage=300.0, Pixel_size=1.0, Cs=2.0, ampcont=10., f_start=0.0, f_stop=-1.0, round_off=1.0, nr1=3, nr2=6):
	"""
	
		1. Estimate envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get 
		   defocus which matches the extracted CTF imprints 
	"""
	from utilities  import generate_ctf
	from morphology import defocus_env_baseline_fit, defocus_guess, defocus_guess1, ctf_1d
	
	#print "CTF params:", voltage, Pixel_size, Cs, wgh, f_start, f_stop, round_off, nr1, nr2, parent
	if f_start == 0 : 	i_start = 0
	else: 			    i_start = int(Pixel_size*2.*len(roo)*f_start)
	if f_stop <= f_start :
		i_stop  = len(roo)
	else: 
		i_stop  = int(Pixel_size*2.*len(roo)*f_stop)
		if i_stop > len(roo): i_stop  = len(roo)

	#print "f_start, i_start, f_stop, i_stop:", f_start, i_start, f_stop, i_stop-1
	rot2 = [0.0]*len(roo)
	for i in range(len(roo)):  rot2[i] = abs(roo[i])
	TE  = defocus_env_baseline_fit(rot2, i_start, i_stop, int(nr1), 4)

	defocus = defocus_guess1(roo, TE, voltage, Cs, Pixel_size, wgh, i_start, i_stop, 2, round_off)

	nx  = int(len(roo)*2)
	ctf = ctf_1d(nx, generate_ctf([defocus, Cs, voltage, Pixel_size, 0.0, ampcont]))

	#from utilities import write_text_file
	#write_text_file([range(len(roo)), roo, ctf, TE], "procrf.txt")

	return defocus, TE, ctf, i_start, i_stop



def make_real(t):
	from utilities import model_blank

	nx = t.get_ysize()
	ny2 = nx//2
	q = model_blank(nx,nx)
	for iy in range(0,nx):
		jy = ny2-iy
		if(jy<0): jy += nx
		jm = nx-jy
		if( jy == 0 ): jm = 0
		for ix in range(0,nx,2):
			tr = t.get_value_at(ix,iy)
			jx = ix//2
			q.set_value_at(jx+ny2, jm, tr)
			q.set_value_at(ny2-jx, jy, tr)
	return q


def fastigmatism(amp, data):
	from morphology import ctf2_rimg
	from utilities import generate_ctf
	
	nx = data[0].get_xsize()
	qt = 0.5*nx**2
	bcc = -1.0
	for j in range(90):
		ang = j
		pc = ctf2_rimg(nx, generate_ctf([data[3], data[4], data[5], data[6], 0.0, data[7], amp, ang]) )
		cuc = (pc*data[1]).cmp("dot", data[0], {"mask":data[2], "negative":0, "normalize":1})
		if( cuc > bcc ):
			bcc = cuc
			bamp = amp
			bang = ang
			#print bdef,bamp,bang,bcc
		#else:
		#	print bdef,amp,ang,cuc
	data[-1] = bang
	return  -bcc


def fastigmatism1(amp, data):
	
	from morphology import ctf_rimg
	from utilities import generate_ctf
	
	nx = data[0].get_xsize()
	qt = 0.5*nx**2
	bcc = -1.0
	for j in range(90):
		ang = j
		pc = ctf_rimg(nx, generate_ctf([data[3], data[4], data[5], data[6], 0.0, data[7], amp, ang]) )
		cuc = (pc*data[1]).cmp("dot", data[0], {"mask":data[2], "negative":0, "normalize":1})
		if( cuc > bcc ):
			bcc = cuc
			bamp = amp
			bang = ang
			#print bdef,bamp,bang,bcc
		#else:
		#	print bdef,amp,ang,cuc
	#pc.write_image("pc.hdf")
	#data[0].write_image("true.hdf")
	data[-1] = bang
	return  -bcc


"""
def fastigmatism2(amp, data):
	
	from morphology import ctf_rimg
	from utilities import generate_ctf
	from alignment import ornq
	
	cnx = data[2]//2+1
	#qt = 0.5*nx**2
	pc = ctf_rimg(data[2], generate_ctf([data[3], data[4], data[5], data[6], 0.0, data[7], amp, 0.0]) )
	ang, sxs, sys, mirror, peak = ornq(pc, crefim, 0.0, 0.0, 1, "H", numr, cnx, cnx)# shold be data[0]
	data[-1] = ang
	return  -peak
"""

def fastigmatism3(amp, data):
	from morphology import ctf2_rimg
	from utilities  import generate_ctf
	from alignment  import ornq
	from math       import sqrt
	#  data[0] - crefim
	#  data[1] - numr
	#  data[2] - nx (image is square)
	#  data[8] - astigmatism amplitude
	#  data[9] - mask defining the region of interest

	cnx = data[2]//2+1
	#qt = 0.5*nx**2
	#B = 0.0
	pc = ctf2_rimg(data[2], generate_ctf([data[3], data[4], data[5], data[6], 0.0, data[7], amp, 0.0]) )
	#st = Util.infomask(pc, data[9], True)
	#Util.mul_scalar(pc, 1.0/st[0])
	ang, sxs, sys, mirror, peak = ornq(pc, data[0], [0.0,0.0], [0.0,0.0], 1, "H", data[1], cnx, cnx)
	#print  ang, sxs, sys, mirror, peak
	#exit()
	data[8] = ang
	return  -peak

fastigmatism2 = fastigmatism3

def fastigmatism3_pap(amp, data):
	from morphology import ctf_rimg
	from utilities  import generate_ctf
	from math       import sqrt
	#  data[0] - crefim
	#  data[1] - numr
	#  data[2] - nx (image is square)
	#  data[8] - astigmatism amplitude
	#  data[9] - mask defining the region of interest

	cnx = data[2]//2+1
	#qt = 0.5*nx**2
	#B = 0.0
	pc = ctf_rimg(data[2], generate_ctf([data[3], data[4], data[5], data[6], 0.0, data[7], amp, 0.0]), sign = 0)
	#st = Util.infomask(pc, data[9], True)
	#Util.mul_scalar(pc, 1.0/st[0])
	ang, sxs, sys, mirror, peak = ornq_vpp(pc, data[0], [0.0,0.0], [0.0,0.0], 1, "H", data[1], cnx, cnx)
	#print  ang, sxs, sys, mirror, peak
	#exit()
	data[8] = ang
	return  -peak

def simctf(amp, data):
	from morphology import ctf2_rimg
	from utilities import generate_ctf
	
	nx = data[2]
	qt = 0.5*nx**2
	pc = ctf_rimg(nx, generate_ctf([amp, data[4], data[5], data[6], 0.0, data[7], data[3], data[8]]) )
	bcc = pc.cmp("dot", data[0], {"mask":data[1], "negative":0, "normalize":1})
	return  -bcc

def simctf2(dz, data):
	from morphology import ctf2_rimg
	from utilities import generate_ctf
	
	#nx = data[2]
	#qt = 0.5*nx**2
	#print  data
	#print dz, data[4], data[5], data[6], 0.0, data[7], data[3], data[8]
	pc = ctf_rimg(data[2], generate_ctf([dz, data[4], data[5], data[6], 0.0, data[7], data[3], data[8]]) )
	bcc = pc.cmp("dot", data[0], {"mask":data[1], "negative":0, "normalize":1})
	#print " simctf2   ",amp,-bcc
	return  -bcc

def simctf2_pap(dz, data):
	from morphology import ctf_rimg
	from utilities import generate_ctf
	
	#nx = data[2]
	#qt = 0.5*nx**2
	#print  data
	#print dz, data[4], data[5], data[6], 0.0, data[7], data[3], data[8]
	pc = ctf_rimg(data[2], generate_ctf([dz, data[4], data[5], data[6], 0.0, data[7], data[3], data[8]]), sign=0)
	bcc = pc.cmp("dot", data[0], {"mask":data[1], "negative":0, "normalize":1})
	#print " simctf2   ",amp,-bcc
	return  -bcc

def simctf2out(dz, data):
	from morphology import ctf2_rimg, localvariance
	from utilities import generate_ctf, model_blank, pad
	
	nx = data[2]
	qt = 0.5*nx**2
	pc = ctf2_rimg(nx, generate_ctf([dz, data[4], data[5], data[6], 0.0, data[7], data[3], data[8]]) )
	pc.write_image("ocou2.hdf")
	normpw = localvariance(  data[0]   , nx//8, 0)#has to be changed

	(normpw*data[1]).write_image("ocou1.hdf")
	mm = pad(model_blank(nx//2,nx,1,1.0),nx,nx,1,0.0,-nx//4)
	s = Util.infomask(pc, None, True)
	pc -= s[0]
	pc /= s[1]
	dout = model_blank(nx,nx)
	dout += pc*mm
	s = Util.infomask(normpw, data[1], True)
	dout += ((normpw-s[0])/s[1])*(model_blank(nx,nx,1,1)-mm)*data[1]
	dout.write_image("ocou3.hdf")
	bcc = pc.cmp("dot", data[0], {"mask":data[1], "negative":0, "normalize":1})
	#print " simctf2   ",amp,-bcc
	return  -bcc


def fupw(args, data):
	from morphology import fastigmatism3
	return -fastigmatism3(args[1],[data[0], data[1], data[2], args[0], data[4], data[5], data[6], data[7], data[8], data[9]])


def fupw_pap(args, data):
	from morphology import fastigmatism3_pap
	return -fastigmatism3_pap(args[1],[data[0], data[1], data[2], args[0], data[4], data[5], data[6], data[7], data[8], data[9]])

########################################
# end of functions used by ctfer
########################################

########################################################################
# start of code used for estimation of cross resolution
########################################################################

def simpw1d_crf(defocus, data):
	from morphology import ctf_1d
	import numpy as np
	from utilities import generate_ctf
	
	#[defocus, Cs, volt, Pixel_size, 0.0, ampcont]
	# data[1] - envelope
	ct = data[1]*np.array( ctf_1d(data[2], generate_ctf([defocus, data[4], data[5], data[6], 0.0, data[7], 0.0, 0.0]))[data[8]:data[9]], np.float32)
	return  2.0-sum(data[0]*ct)/np.linalg.norm(ct,2)

def linregnp(y):
	import numpy as np
	ny = len(y)
	ff = type(y[0])
	x = np.array(list(range(ny)), ff)
	return  np.linalg.lstsq(np.vstack([x, np.ones(ny,ff)]).T, y)

def defocusgett_crf(roo, nx, voltage=300.0, Pixel_size=1.0, Cs=2.0, ampcont=0.1, f_start=-1.0, f_stop=-1.0, round_off=1.0, nr1=3, nr2=6, parent=None, DEBug=False):
	#(roo, nx, voltage=300.0, Pixel_size=1.0, Cs=2.0, ampcont=0.1, f_start=-1.0, f_stop=-1.0, round_off=1.0, nr1=3, nr2=6, parent=None):
	"""
	
		1. Estimate envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get 
		   defocus which matches the extracted CTF imprints 
	"""
	from utilities  import generate_ctf
	from morphology import bracket_def, goldsearch_astigmatism, ctflimit, simpw1d_crf, ctf_1d
	import numpy as np


	#print "CTF params:", voltage, Pixel_size, Cs, wgh, f_start, f_stop, round_off, nr1, nr2, parent

	if f_start == 0 : 	    i_start = 0
	else: 			        i_start = int(Pixel_size*nx*f_start+0.5)
	if f_stop <= f_start :
		i_stop  = len(roo)
		adjust_fstop = True
	else:
		i_stop  = min(len(roo), int(Pixel_size*nx*f_stop+0.5))
		adjust_fstop = False

	nroo = len(roo)
	"""
	#  THIS COULD ALSO WORK< TRY IT
	rot2 = [0.0]*len(roo)
	for i in xrange(len(roo)):  rot2[i] = abs(roo[i])
	TE  = defocus_env_baseline_fit(rot2, i_start, i_stop, int(nr1), 4)  # This does envelope
	"""
	#print "f_start, i_start, f_stop, i_stop:", f_start, i_start, f_stop, i_stop-1
	#  THERE IS NO NEED FOR BASELINE!!!
	#write_text_file([roo,baseline,subpw],"dbg.txt")
	#print "IN defocusgett  ",np.min(subpw),np.max(subpw)
	envelope = np.array([1.0]*i_stop*2, np.float32) # movingaverage(  abs( np.array(roo, np.float32) )   , nroo//4, 3)
	#write_text_file([roo,baseline,subpw],"dbgt.txt")

	#print "IN defocusgett  ",np.min(subpw),np.max(subpw),np.min(envelope)
	#envelope = np.ones(nroo, np.float32)
	defocus = 0.0
	data = [roo[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
	#for i in xrange(nroo):
	#	print  i,"   ",roo[i],"   ",baseline[i],"   ",subpw[i],"   ",envelope[i]
	h = 0.1
	#def1, def2 = bracket(simpw1d, data, h)
	#if DEBug:  print "first bracket ",def1, def2,simpw1d(def1, data),simpw1d(def2, data)
	#def1=0.1
	ndefs = 18
	defound = []
	for  idef in range(ndefs):
		def1 = (idef+1)*0.5
		def1, def2 = bracket_def(simpw1d_crf, data, def1, h)
		#if DEBug:  print "second bracket ",idef,def1, def2,simpw1d_crf(def1, data),simpw1d(def2, data),h
		def1, val2 = goldsearch_astigmatism(simpw1d_crf, data, def1, def2, tol=1.0e-3)
		#if DEBug:  print "golden ",idef,def1, val2,simpw1d_crf(def1, data)
		if def1>0.0:  defound.append([val2,def1])
	defound.sort()
	del defound[3:]
	if DEBug:  print(" BEST DEF CANDIDATES",defound)
	if adjust_fstop:
		from morphology import ctflimit
		newstop,fnewstop = ctflimit(nx, defound[0][1], Cs, voltage, Pixel_size)
		if DEBug:  
			print("newstop  ",int(newstop),fnewstop,i_stop,newstop,nx, defound[0][1])
		if( newstop != i_stop):
			i_stop = newstop
			data = [roo[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
			h = 0.05
			for idef in range(3):
				def1, def2 = bracket_def(simpw1d_crf, data, defound[idef][1], h)
				if DEBug:  print(" adjusted def ",def1,def2)
				def1, val2 = goldsearch_astigmatism(simpw1d_crf, data, def1, def2, tol=1.0e-3)
				if DEBug:  print("adjusted golden ",def1, val2,simpw1d_crf(def1, data))
				if def1>0.0:  defound[idef] = [val2,def1]
			defound.sort()
	def1 = defound[0][1]
	if DEBug: print(" ultimate defocus",def1,defound)

	#defocus = defocus_guessn(Res_roo, voltage, Cs, Pixel_size, ampcont, i_start, i_stop, 2, round_off)
	#print simpw1d(def1, data),simpw1d(4.372, data)
	"""
	def1 = 0.02
	def2 = 10.
	def1, def2 = goldsearch_astigmatism(simpw1d, data, def1, def2, tol=1.0e-3)
	print "golden ",def1, def2,simpw1d(def1, data)
	"""
	if DEBug and False:
		qm = 1.e23
		toto = []
		for i in range(1000,100000,5):
			dc = float(i)/10000.0
			qt = simpw1d_crf(dc, data)
			toto.append([dc,qt])
			if(qt<qm):
				qm=qt
				defi = dc
		write_text_row(toto,"toto1.txt")
		print(" >>>>>>>>>  ",defi,simpw1d_crf(defi, data))#,generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont])
		#def1 = defi
	#exit()
	ctf1d = ctf_1d(nx, generate_ctf([def1, Cs, voltage, Pixel_size, 0.0, ampcont]))

	return def1, ctf1d, None, envelope, i_start, i_stop

def envelopegett_crf(defold, roo, nx, voltage=300.0, Pixel_size=1.0, Cs=2.0, ampcont=0.1, f_start=-1.0, f_stop=-1.0, round_off=1.0, nr1=3, nr2=6, parent=None, DEBug=False):
	"""
	
		
	"""
	from utilities  import generate_ctf
	from morphology import ctf_1d
	import numpy as np


	#print "CTF params:", voltage, Pixel_size, Cs, wgh, f_start, f_stop, round_off, nr1, nr2, parent

	if f_start == 0 : 	    i_start = 0
	else: 			        i_start = int(Pixel_size*nx*f_start+0.5)
	if f_stop <= f_start :
		i_stop  = len(roo)
		adjust_fstop = True
	else:
		i_stop  = min(len(roo), int(Pixel_size*nx*f_stop+0.5))
		adjust_fstop = False

	nroo = len(roo)
	"""
	#  THIS COULD ALSO WORK< TRY IT
	rot2 = [0.0]*len(roo)
	for i in xrange(len(roo)):  rot2[i] = abs(roo[i])
	TE  = defocus_env_baseline_fit(rot2, i_start, i_stop, int(nr1), 4)  # This does envelope
	"""
	#print "f_start, i_start, f_stop, i_stop:", f_start, i_start, f_stop, i_stop-1
	#  THERE IS NO NEED FOR BASELINE!!!
	#write_text_file([roo,baseline,subpw],"dbg.txt")
	#print "IN defocusgett  ",np.min(subpw),np.max(subpw)
	envelope = []#movingaverage(  abs( np.array(roo, np.float32) )   , nroo//4, 3)

	if adjust_fstop:
		from morphology import ctflimit
		newstop,fnewstop = ctflimit(nx, defold, Cs, voltage, Pixel_size)
		if DEBug:  print("newstop  ",int(newstop*0.7),fnewstop*0.7,i_stop)
	
	ctf1d = ctf_1d(nx, generate_ctf([defold, Cs, voltage, Pixel_size, 0.0, ampcont]))

	return envelope, i_start, i_stop

def fufu(args,data):
	from morphology import fastigmatism2
	return -fastigmatism2(args[1],[data[0], data[1], data[2], args[0], data[4], data[5], data[6], data[7], data[8]])

#  
# NOTE: 2016/03/21 Toshio Moriya
# getastcrfNOE() function does not work with the new output format of cter_mrk()
# 
def getastcrfNOE(refvol, datfilesroot, voltage=300.0, Pixel_size= 1.264, Cs = 2.0, wgh = 7.0, kboot=16, DEBug = False):
	"""

	#####################   Estimation from crossresolution   ############################


	#  mpirun -np 1 python /Volumes/pawel//ED/getastcrf.py <reference volume>  <name string of the data file>


	"""

	from applications import MPI_start_end
	from utilities import read_text_file, write_text_file, get_im, model_blank, model_circle, amoeba, generate_ctf
	from sys import exit
	import numpy as np
	import os
	from mpi  import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD, mpi_barrier
	from fundamentals import tilemic, rot_avg_table
	from morphology import threshold, bracket_def, bracket, goldsearch_astigmatism, defocus_baseline_fit, simpw1d, movingaverage, localvariance, defocusgett, defocus_guessn, defocusget_from_crf, make_real, fastigmatism, fastigmatism1, fastigmatism2, fastigmatism3, simctf, simctf2, simctf2out, fupw,ctf2_rimg
	from alignment import Numrinit, ringwe
	from statistics import table_stat
	from pixel_error import angle_ave

	myid = mpi_comm_rank(MPI_COMM_WORLD)
	ncpu = mpi_comm_size(MPI_COMM_WORLD)
	main_node = 0
	
	UseOldDef = True

	#f_start = 0.022
	#f_stop  = 0.24
	
	#lenroot = len(nameroot)


	#ll = read_text_file("lookup",-1)
	#ll = [[367], [12031], [25]]
	ll = read_text_file("lookup")

	#set_start, set_end = MPI_start_end(len(ll[0]), ncpu, myid)
	#for k in xrange(3):
	#	ll[k] = map(int,ll[k])[set_start:set_end]

	volft,kb = prep_vol(get_im(refvol))


	totresi = []
	for ifi in range(len(ll)):
		#  There should be something here that excludes sets with too few images
		#namics = datfilesroot+"%05d_%06d"%(ll[0][ifi],ll[1][ifi])
		namics = datfilesroot+"%05d"%ll[ifi]
		d = EMData.read_images( namics )
		nimi = len(d)
		nx = d[0].get_xsize()
		ny2 = nx//2
		if UseOldDef:
			defold, csi, volti, apixi, bfcti, ampconti, astampi, astangi = get_ctf(d[0])
			if DEBug:  print(" USING OLD CTF  ",defold, csi, volti, apixi, bfcti, ampconti, astampi, astangi)

		fa  = [None]*nimi
		fb  = [None]*nimi
		fbc = [None]*nimi

		for imi in range(nimi):
			phi,theta,psi,tx,ty = get_params_proj(d[imi])
			# next is test
			#psi = 0.0
			#d[imi] = prgs(volft, kb, [phi,theta,psi,-tx,-ty])
			####
			fa[imi] = fft( d[imi] )
			#prgs(volft, kb, [phi,theta,psi,-tx,-ty]).write_image("bdb:projs",imi)
			#fb[imi] = filt_ctf(fa[imi] , generate_ctf([1.9, Cs, voltage, Pixel_size, 0.0, wgh, 0.9, 177.]),False) + fft(model_gauss_noise(2., nx,nx)) #!!!!!!!!!!fft(get_im("bdb:projs",imi))  #
	
			#  next modified for test
			fb[imi] = fft( prgs(volft, kb, [phi,theta,psi,-tx,-ty]) )   #fa[imi].copy()#
			#fbc[imi] = fb[imi].conjg()  #  THIS IS WRONG PAP 06/19/2018
			# next is test
			#fa[imi] = filt_ctf(fa[imi] , generate_ctf([defold, Cs, voltage, Pixel_size, 0.0, wgh, 0.9, 77.]),False) + fft(model_gauss_noise(2., nx,nx)) 

		del d  # I do not think projections are needed anymore
		adefocus = [0.0]*kboot
		aamplitu = [0.0]*kboot
		aangle   = [0.0]*kboot
	
		if True:  #try:
			for nboot in range(kboot):
				if(nboot == 0): boot = list(range(nimi))
				else:
					from random import randint
					for imi in range(nimi): boot[imi] = randint(0,nimi-1)
	
				qs = model_blank(nx,nx)
				qa = model_blank(nx,nx)
				qb = model_blank(nx,nx)
				crf1d = []
				for imboot in range(nimi):
					imi = boot[imboot]
		
					temp = fsc(fa[imi],fb[imi])[1]
					if( len(crf1d) == 0 ): crf1d = [0.0]*len(temp)
					for k in range(len(temp)):  crf1d[k] += temp[k]
					t  = make_real( Util.muln_img(fa[imi], fbc[imi]) )
					Util.mul_scalar(t, 1.0/(float(nx)**4))
					Util.add_img(qs , t)

					Util.add_img(qa, periodogram(fa[imi]))
					Util.add_img(qb, periodogram(fb[imi]))

				for k in range(len(temp)):  crf1d[k] /= nimi
				"""
				sroo[0] = sroo[1]
				aroo[0] = aroo[1]
				sroo = (sroo-aroo**2/nimi)/nimi
				aroo /= nimi
				roo  /= nimi
				qa /= nimi
				"""
	
				from math import sqrt
				nc = nx//2
				tqa = [0.0]*(nc+1)
				for i in range(nx):
					for j in range(nx):
						r = sqrt((i-nc)**2 + (j-nc)**2)
						ir = int(r)
						if(ir<nc):
							dr = r - ir
							qqa = qa.get_value_at(i,j)
							tqa[ir]   += (1.0-dr)*qqa
							tqa[ir+1] +=       dr*qqa
				for i in range(nc+1): tqa[i] = sqrt(max(tqa[i],0.0))

				divs = model_blank(nx, nx, 1, 1.0)
				for i in range(nx):
					for j in range(nx):
						r = sqrt((i-nc)**2 + (j-nc)**2)
						ir = int(r)
						if(ir<nc):
							dr = r - ir
							divs.set_value_at(i,j,  (1.-dr)*tqa[ir] + dr*tqa[ir+1] )
				#if(nboot == 0): qs.write_image("rs1.hdf")
				Util.div_img(qs, divs)
				qs.set_value_at(ny2,ny2,1.0)
				#if(nboot == 0): write_text_file(crf1d,"crf1d.txt")
				#if(nboot == 0): qs.write_image("rs2.hdf")


				sroo = rot_avg_table(qs)
				lenroo = len(sroo)
				#  Find a break point
				bp = 1.e23
				for i in range(1,max(3,lenroo//4)):
					if( sroo[i] <0.0 ):
						istart = max(3,i//2)
						break

				#istart = 25
				#print istart
	
				f_start = istart/(Pixel_size*nx)
				#print namics[ifi],istart,f_start


				if UseOldDef:
					envelope, istart, istop = envelopegett_crf(defold, crf1d, nx, voltage=voltage, Pixel_size=Pixel_size, Cs=Cs, ampcont=wgh, f_start=f_start, f_stop=-1.0, round_off=1.0, nr1=3, nr2=6, parent=None, DEBug=DEBug)
					defc = defold
				else:
					defc, ctf1d, baseline, envelope, istart, istop = defocusgett_crf(crf1d, nx, voltage=voltage, Pixel_size=Pixel_size, Cs=Cs, ampcont=wgh, f_start=f_start, f_stop=-1.0, round_off=1.0, nr1=3, nr2=6, parent=None, DEBug=DEBug)
					if DEBug:  print("  RESULT ",namics,defc, istart, istop)
					if DEBug:
						freq = list(range(len(crf1d)))
						for i in range(len(crf1d)):  freq[i] = float(i)/nx/Pixel_size
						write_text_file([freq, crf1d, ctf1d, envelope.tolist()],"ravg%05d.txt"%ifi)
				#mpi_barrier(MPI_COMM_WORLD)
				#   NOT USING ENVELOPE!
				"""
				en = envelope.tolist()
				en.append(en[-1])
				envl = model_blank(nx, nx, 1, 1.0)
				for i in xrange(nx):
					for j in xrange(nx):
						r = sqrt((i-nc)**2 + (j-nc)**2)
						ir = int(r)
						if(ir<nc):
							dr = r - ir
							envl.set_value_at(i,j,  (1.-dr)*en[ir] + dr*en[ir+1] )
				"""

				#exit()
				#istop = nx//4
				#mask = model_circle(istop-1,nx,nx)*(model_blank(nx,nx,1,1.0)-model_circle(istart,nx,nx))
				#qse = qa*envl
				#(qse*mask).write_image("rs2.hdf")
				#qse.write_image("rs3.hdf")
				##  SIMULATION
				#bang = 0.7
				#qse = ctf2_rimg(nx, generate_ctf([defc,Cs,voltage,Pixel_size,0.0,wgh, bang, 37.0]) )
				#qse.write_image("rs3.hdf")

				mask = model_circle(istop-1,nx,nx)*(model_blank(nx,nx,1,1.0)-model_circle(istart,nx,nx))
				qse = qs #* envl
				#if(nboot == 0): (qs*mask).write_image("rs5.hdf")


				cnx = nx//2+1
				cny = cnx
				mode = "H"
				numr = Numrinit(istart, istop, 1, mode)
				wr   = ringwe(numr, mode)

				crefim = Util.Polar2Dm(qse, cnx, cny, numr, mode)
				Util.Frngs(crefim, numr)
				Util.Applyws(crefim, numr, wr)


				"""
				astdata = [crefim, numr, nx, 1.1, Cs, voltage, Pixel_size, wgh, 0.]
				junk = fastigmatism2(0.5,astdata)
				bang = astdata[-1]
				print "   IHIHIHIHIHI   ",bang,junk
				#exit()
				"""

				#pc = ctf2_rimg(nx,generate_ctf([defc,Cs,voltage,Pixel_size,0.0,wgh]))
				#print ccc(pc*envl, subpw, mask)

				bang = 0.0
				bamp = 0.0
				bdef = defc
				bold = 1.e23
				dstep = 0.1
				while( True):
					data = [qse, mask, nx, bamp, Cs, voltage, Pixel_size, wgh, bang]
					h = 0.05*bdef
					#print "  bdef  at the beginning of while loop   ",nboot,bdef
					amp1, amp2 = bracket_def(simctf, data, bdef*0.9, h)
					#print "bracketing of the defocus  ",nboot,amp1, amp2
					amp1, val2 = goldsearch_astigmatism(simctf, data, amp1, amp2, tol=1.0e-3)
					#print "golden defocus ",amp1, val2,simctf(amp1, data)
					#,simctf(1.1, [qse, mask, nx, 0.5, Cs, voltage, Pixel_size, wgh, 77.])
					#print " ttt ",time()-srtt
					#bdef = 1.1
					#bamp = 0.2
					#write_text_file(ctf_1d(data[2], generate_ctf([amp1, data[4], data[5], data[6], 0.0, data[7], 0.0, 0.0])), "gctf1d.txt")
					#exit()
	
					astdata = [crefim, numr, nx, bdef, Cs, voltage, Pixel_size, wgh, bang]
	
					h = 0.01
					amp1,amp2 = bracket(fastigmatism2, astdata, h)
					#print "  astigmatism bracket  ",nboot,amp1,amp2,astdata[-1]
					#print " ttt ",time()-srtt
	
					bamp, bcc = goldsearch_astigmatism(fastigmatism2, astdata, amp1, amp2, tol=1.0e-3)
	
					#junk = fastigmatism2(bamp,astdata)
					#bang = astdata[-1]
					#print " ang within the loop   ",bdef, bamp,astdata[-1],junk, fastigmatism2(0.5,[crefim, numr, nx, 1.1, Cs, voltage, Pixel_size, wgh, bang])
					#print  fastigmatism2(0.0,astdata)
					#print astdata[-1]
					#print "  golden search ",bamp,data[-1], fastigmatism2(bamp,data), fastigmatism2(0.0,data)
					#print " ttt ",time()-srtt
					#bamp = 0.5
					#bang = 277
					dama = amoeba([bdef, bamp],[0.2,0.2], fufu, 1.e-4,1.e-4,500, astdata)

					if DEBug:  print("AMOEBA    ",nboot,dama)
					bdef = dama[0][0]
					bamp = dama[0][1]
					astdata = [crefim, numr, nx, bdef, Cs, voltage, Pixel_size, wgh, bang]
					junk = fastigmatism2(bamp, astdata)
					bang = astdata[-1]
					if DEBug:  print(" after amoeba ", nboot,bdef, bamp, bang)
					#  The looping here is blocked as one shot at amoeba is good enough.  To unlock it, remove - from bold.
					if(bcc < -bold): bold = bcc
					else:           break

				adefocus[nboot] = bdef
				aamplitu[nboot] = bamp
				aangle[nboot]   = bang
				#from sys import exit
				if DEBug:  print("this is what I found  ",nboot,bdef,bamp,bang)
				#exit()

			#print " ttt ",time()-srtt
			#from sys import exit
			#exit()
			ad1,ad2,ad3,ad4 = table_stat(adefocus)
			reject = []
			thr = 3*sqrt(ad3)
			for i in range(len(adefocus)):
				if(abs(adefocus[i]-ad1)>thr):
					if DEBug:  print(adefocus[i],ad1,thr)
					reject.append(i)
			if(len(reject)>0):
				if DEBug:  print("  Number of rejects  ",namics,len(reject))
				for i in range(len(reject)-1,-1,-1):
					del adefocus[i]
					del aamplitu[i]
					del aangle[i]
			if(len(adefocus)<2):
				print("  After rejection of outliers too few estimated defocus values for :",namics)
			else:
				#print "adefocus",adefocus
				#print  "aamplitu",aamplitu
				#print "aangle",aangle
				ad1,ad2,ad3,ad4 = table_stat(adefocus)
				bd1,bd2,bd3,bd4 = table_stat(aamplitu)
				cd1,cd2 = angle_ave([2*q for q in aangle])  # Have to use this trick as the function works for range [0,360]
				cd1/=2
				cd2/=2
				temp = 0.0
				print(namics,ad1, Cs, voltage, Pixel_size, temp, wgh, bd1, cd1, sqrt(max(0.0,ad2)),sqrt(max(0.0,bd2)),cd2) 
				totresi.append( [ namics, ad1, Cs, voltage, Pixel_size, temp, wgh, bd1, cd1, sqrt(max(0.0,ad2)),sqrt(max(0.0,bd2)),cd2 ])
				#if ifi == 4 : break
				"""
				for i in xrange(len(ssubroo)):
					asubroo[i] /= kboot
					ssubroo[i]  = sqrt(max(0.0, ssubroo[i]-kboot*asubroo[i]**2)/kboot)
					sen[i]     /= kboot
				"""
				#ctf_rimg(nx,generate_ctf([ad1, Cs, voltage, Pixel_size, temp, wgh, bd1, cd1])).write_image("ctf1.hdf")
				lnsb = len(crf1d)
	
				try:		crot2 = rotavg_ctf(ctf_rimg(nx,generate_ctf([ad1, Cs, voltage, Pixel_size, temp, wgh, bd1, cd1])), ad1, Cs, voltage, Pixel_size, bd1, cd1)[:lnsb]
				except:     crot2 = [0.0]*lnsb
				try:		pwrot2 = rotavg_ctf(qs, ad1, Cs, voltage, Pixel_size, bd1, cd1)[:lnsb]
				except:     pwrot2 = [0.0]*lnsb
				try:		crot1 = rotavg_ctf(ctf_rimg(nx,generate_ctf([ad1, Cs, voltage, Pixel_size, temp, wgh, bd1, cd1])), ad1, Cs, voltage, Pixel_size, 0.0, 0.0)[:lnsb]
				except:     crot1 = [0.0]*lnsb
				try:		pwrot1 = rotavg_ctf(qs, ad1, Cs, voltage, Pixel_size, 0.0, 0.0)[:lnsb]
				except:     pwrot1 = [0.0]*lnsb
				freq = list(range(lnsb))
				for i in range(len(freq)):  freq[i] = float(i)/nx/Pixel_size
				#fou = "crfrot/rotinf%05d_%06d.txt"%(ll[0][ifi],ll[1][ifi])
				fou = "crfrot/rotinf%05d.txt"%(ll[ifi])
				#  #1 - rotational averages without astigmatism, #2 - with astigmatism
				write_text_file([list(range(len(crot1))), freq, pwrot1, crot1, pwrot2, crot2],fou)
				cmd = "echo "+"    "+namics+"  >>  "+fou
				os.system(cmd)
		else:  #except:
			print(namics,"     FAILED")
	#from utilities import write_text_row
	outf = open( "partcrf/partcrf_%05d"%myid, "w")
	for i in range(len(totresi)):
		for k in range(1,len(totresi[i])): outf.write("  %12.5g"%totresi[i][k])
		outf.write("  %s\n"%totresi[i][0])
	outf.close()		

########################################################################
# end of code used for estimation of cross resolution
########################################################################


################
#
#  CTER code (07/10/2017)
#
################
# 
# 
def cter_vpp(input_image_path, output_directory, selection_list = None, wn = 512, 
			pixel_size = -1.0, Cs = 2.0, voltage = 300.0, wgh = 10.0, f_start = -1.0, f_stop = -1.0, \
			kboot = 16, overlap_x = 50, overlap_y = 50, edge_x = 0, edge_y = 0, \
			check_consistency = False, stack_mode = False, debug_mode = False, \
			program_name = "cter_vpp() in morphology.py", vpp_options = [], \
			RUNNING_UNDER_MPI = False, main_mpi_proc = 0, my_mpi_proc_id = 0, n_mpi_procs = 1):
	"""
	Arguments
		input_image_path  :  file name pattern for Micrographs Modes (e.g. 'Micrographs/mic*.mrc') or particle stack file path for Stack Mode (e.g. 'bdb:stack'; must be stack_mode = True).
		output_directory  : output directory
	"""
	from   EMAN2 import periodogram
	from   EMAN2db import db_check_dict, db_parse_path
	from   applications import MPI_start_end
	from   utilities import read_text_file, write_text_file, get_im, model_blank, model_circle, amoeba, generate_ctf
	from   utilities import if_error_then_all_processes_exit_program
	from   utilities import wrap_mpi_bcast
	from   fundamentals import tilemic, rot_avg_table, resample
	from   morphology   import threshold, bracket_def, bracket, goldsearch_astigmatism
	from   morphology   import defocus_baseline_fit, simpw1d, movingaverage, localvariance, defocusgett
	from   morphology   import defocus_guessn, defocusget_from_crf, make_real
	from   morphology   import fastigmatism, fastigmatism1, fastigmatism2, fastigmatism3, simctf, simctf2, simctf2out, fupw,ctf2_rimg
	from   alignment    import Numrinit, ringwe
	from   statistics   import table_stat
	from   pixel_error  import angle_ave
	from   global_def   import ERROR
	import global_def

	from   sys import exit
	import numpy as np
	import os
	import glob
	from   time import time
	from   inspect      import currentframe, getframeinfo
	from   mpi import MPI_COMM_WORLD, mpi_barrier

	# ====================================================================================
	# Prepare processing
	# ====================================================================================
	#  vpp_options = [defocus_min,  defocus_max,  defocus_step,  phase_min,  phase_max,  phase_step]
	# ------------------------------------------------------------------------------------

	# ------------------------------------------------------------------------------------
	# Find the CTER Running Mode before checking error conditions
	# ------------------------------------------------------------------------------------
	i_enum = -1; idx_cter_mode_invalid       = i_enum; 
	i_enum += 1; idx_cter_mode_all_mics      = i_enum  # All Micrographs Mode - Process all s in a directory
	i_enum += 1; idx_cter_mode_selected_mics = i_enum  # Selected Micrographs Mode - Process all s in a selection list file
	i_enum += 1; idx_cter_mode_single_mic    = i_enum  # Single Micrograph Mode - Process a single 
	i_enum += 1; idx_cter_mode_stack         = i_enum  # Stack Mode - Process a stack (Advanced Option)
	i_enum += 1; idx_cter_mode_counts        = i_enum

	cter_mode_idx = idx_cter_mode_invalid
	cter_mode_name = None
	if stack_mode == False:
		# One of three Micrograph Modes
		# For any of Micrograph Modes, input image file name must be a file path pattern containing wild card "*" 
		if selection_list == None:
			# User did not use selection list option 
			# -> All Micrographs Mode
			cter_mode_idx = idx_cter_mode_all_mics
			cter_mode_name = "All Micrographs Mode"
		else:
			if os.path.splitext(selection_list)[1] == ".txt":
				# User specified a selection list text file path containing".txt" extension through selection list option
				# -> Selected Micrographs Mode
				cter_mode_idx = idx_cter_mode_selected_mics
				cter_mode_name = "Selected Micrographs Mode"
			else: 
				# User specified an image file path (a non-text file path) through selection list option
				# -> Single Micrograph Mode
				cter_mode_idx = idx_cter_mode_single_mic
				cter_mode_name = "Single Micrograph Mode"
	else: 
		# (Particle) Stack Mode
		cter_mode_idx = idx_cter_mode_stack
		cter_mode_name = "Stack Mode"
	
	if my_mpi_proc_id == main_mpi_proc:
		print(" ")
		print(("----- Running with %s -----" % (cter_mode_name)))
	
	# ------------------------------------------------------------------------------------
	# Check mode-dependent error conditions of input arguments and options if abort is necessary. All nodes do this checking
	# ------------------------------------------------------------------------------------
	error_message_list = [] # List of error messages. If no error is found, the length should be zero
	if not stack_mode:
		
		# Check error conditions applicable to any of Micrograph Mode 
		if input_image_path.find("*") == -1:
			error_message_list.append("Input image file path (%s) for %s must be a  path pattern containing wild card (*). Please check input_image_path argument." % (input_image_path, cter_mode_name))
		
		if input_image_path[:len("bdb:")].lower() == "bdb:":
			error_message_list.append("BDB file can not be selected as input image file path (%s) for %s. Please check input_image_path argument and convert the image format." % (input_image_path, cter_mode_name))
		
		# Check error conditions applicable to Selected Micrographs Mode 
		if cter_mode_idx == idx_cter_mode_selected_mics:
			if not os.path.exists(selection_list): 
				error_message_list.append("Selection list text file specified by selection_list option (%s) for %s does not exists. Please check selection_list option." % (selection_list, cter_mode_name))
		
		if cter_mode_idx == idx_cter_mode_single_mic:
			if not os.path.exists(os.path.join(os.path.dirname(input_image_path), os.path.basename(selection_list))): 
				error_message_list.append("Micrograph specified by selection_list option (%s) for %s does not exist. Please check selection_list option." % (selection_list, cter_mode_name))
			# 
			if RUNNING_UNDER_MPI and n_mpi_procs != 1:
				error_message_list.append("%s supports only a single processor version. Please change MPI settings." % (cter_mode_name))
		
	else: 
		# Check error conditions
		if input_image_path.find("*") != -1:
			error_message_list.append("Stack file path specified by input_image_path (%s) for %s should not contain wild card (*). Please check input_image_path argument." % (input_image_path, cter_mode_name))
		
		is_not_found_input_image_file = False
		if input_image_path[:len("bdb:")].lower() == "bdb:":
			if not db_check_dict(input_image_path): 
				is_not_found_input_image_file = True
		else:
			if not os.path.exists(input_image_path): 
				is_not_found_input_image_file = True
		if is_not_found_input_image_file:
			error_message_list.append("Stack file specified by input_image_path (%s) for %s does not exist. Please check input_image_path argument." % (input_image_path, cter_mode_name))
		
		if RUNNING_UNDER_MPI and n_mpi_procs != 1:
			error_message_list.append("%s supports only a single processor version. Please change MPI settings." % (cter_mode_name))
	
	# --------------------------------------------------------------------------------
	# check output-related error conditions (mode-independent). All nodes do this checking
	# --------------------------------------------------------------------------------
	if os.path.exists(output_directory):
		error_message_list.append("Output directory (%s) exists already. Please check output_directory argument." % (output_directory))
	
	# --------------------------------------------------------------------------------
	# Check error conditions of options (mode-independent). All nodes do this checking
	# --------------------------------------------------------------------------------
	if pixel_size <= 0.0:
		error_message_list.append("Pixel size (%f) must not be negative. Please set a pasitive value larger than 0.0 to pixel_size option." % (pixel_size))

	if wn <= 0.0:
		error_message_list.append("CTF window size (%d) must not be negative. Please set a valid value larger than 0 to wn option." % (wn))
	
	# --------------------------------------------------------------------------------
	# Print all error messages and abort the process if necessary.
	# --------------------------------------------------------------------------------
	error_status = None
	if len(error_message_list) > 0:
		# Detected error! Print all error messages
		if my_mpi_proc_id == main_mpi_proc:
			print(" ")
			for error_message in error_message_list:  
				print(("ERROR!!! %s" % (error_message)))
		error_status = ("Detected %d error(s) related to arguments and options. Run %s -h for help. Exiting..." % (len(error_message_list), program_name), getframeinfo(currentframe()))
	if_error_then_all_processes_exit_program(error_status)
	if RUNNING_UNDER_MPI:
		# Wait for all mpi processes to check error conditions, especially existence of output directory
		# Without this barrier, main mpi process can create output directory before some child mpi process check this error.
		mpi_barrier(MPI_COMM_WORLD)
	del error_message_list # Don't need this anymore
	
	# ------------------------------------------------------------------------------------
	# Check warning conditions of options
	# ------------------------------------------------------------------------------------
	if my_mpi_proc_id == main_mpi_proc:
		if stack_mode:
			if selection_list != None:
				print(" ")
				print(("WARNING!!! --selection_list option will be ignored in %s." % (cter_mode_name)))
			if wn != 512:
				print(" ")
				print(("WARNING!!! --wn option will be ignored in %s." % (cter_mode_name)))
			if overlap_x != 50:
				print(" ")
				print(("WARNING!!! --overlap_x option will be ignored in %s." % (cter_mode_name)))
			if overlap_y != 50:
				print(" ")
				print(("WARNING!!! --overlap_y option will be ignored in %s." % (cter_mode_name)))
			if edge_x != 0:
				print(" ")
				print(("WARNING!!! --edge_x option will be ignored in %s." % (cter_mode_name)))
			if edge_y != 0:
				print(" ")
				print(("WARNING!!! --edge_y option will be ignored in %s." % (cter_mode_name)))
			if check_consistency:
				print(" ")
				print(("WARNING!!! --check_consistency option will be ignored in %s." % (cter_mode_name)))

	# ====================================================================================
	# Create the input file path list and also check input-related error conditions if abort is necessary.
	# ====================================================================================
	input_file_path_list = []
	if not stack_mode:
		# --------------------------------------------------------------------------------
		# Prepare the variables for all sections in  mode case 
		# --------------------------------------------------------------------------------
		# Micrograph basename pattern (directory path is removed from  path pattern)
		mic_pattern = input_image_path
		mic_basename_pattern = os.path.basename(mic_pattern)

		# Global entry dictionary (all possible entries from all lists) for all mic id substring
		global_entry_dict = {} # mic id substring is the key
		subkey_input_mic_path = "Input Micrograph Path"
		subkey_selected_mic_basename = "Selected Micrograph Basename"
	
		# List keeps only id substrings of s whose all necessary information are available
		valid_mic_id_substr_list = [] 
		
		# --------------------------------------------------------------------------------
		# Obtain the list of  id sustrings using a single CPU (i.e. main mpi process)
		# --------------------------------------------------------------------------------
		# NOTE: Toshio Moriya 2016/11/15
		# The below is not a real while.  
		# It gives if-statements an opportunity to use break when errors need to be reported
		# However, more elegant way is to use 'raise' statement of exception mechanism...
		# 
		error_status = None
		while my_mpi_proc_id == main_mpi_proc:			
			# --------------------------------------------------------------------------------
			# Prepare variables for this section
			# --------------------------------------------------------------------------------
			# Prefix and suffix of  basename pattern 
			# to find the head/tail indices of  id substring
			mic_basename_tokens = mic_basename_pattern.split('*')
			# Find head index of  id substring
			mic_id_substr_head_idx = len(mic_basename_tokens[0])
		
			# --------------------------------------------------------------------------------
			# Register  id substrings found in the input directory (specified by  path pattern)
			# to the global entry dictionary
			# --------------------------------------------------------------------------------
			# Generate the list of  paths in the input directory
			print(" ")
			print("Checking the input directory...")
			input_mic_path_list = glob.glob(mic_pattern)
			# Check error condition of input  file path list
			print(("Found %d micrographs in %s." % (len(input_mic_path_list), os.path.dirname(mic_pattern))))
			if len(input_mic_path_list) == 0:
				# The result shouldn't be empty if the specified  file name pattern is invalid
				error_status = ("There are no micrographs whose paths match with the specified file path pattern (%s) for %s. Please check input_image_path. Run %s -h for help." % (mic_pattern, cter_mode_name, program_name), getframeinfo(currentframe()))
				break
		
			# Register  id substrings to the global entry dictionary
			for input_mic_path in input_mic_path_list:
				# Find tail index of  id substring and extract the substring from the  name
				input_mic_basename = os.path.basename(input_mic_path)
				mic_id_substr_tail_idx = input_mic_basename.index(mic_basename_tokens[1])
				mic_id_substr = input_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
				if not mic_id_substr in global_entry_dict:
					# print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from input_mic_path_list " % (mic_id_substr))
					global_entry_dict[mic_id_substr] = {}
				global_entry_dict[mic_id_substr][subkey_input_mic_path] = input_mic_path
		
			# --------------------------------------------------------------------------------
			# Register  id substrings found in the selection list
			# to the global entry dictionary
			# --------------------------------------------------------------------------------
			# Generate the list of selected  paths in the selection file
			selected_mic_path_list = []
			# Generate  lists according to the execution mode
			if cter_mode_idx == idx_cter_mode_all_mics:
				# Treat all s in the input directory as selected ones
				selected_mic_path_list = input_mic_path_list
			else:
				if os.path.splitext(selection_list)[1] == ".txt":
					print(" ")
					print("Checking the selection list...")
					selected_mic_path_list = read_text_file(selection_list)
				
					# Check error condition of  entry lists
					print(("Found %d microgarph entries in %s." % (len(selected_mic_path_list), selection_list)))
					if len(selected_mic_path_list) == 0:
						error_status = ("The provided  list file (%s) for %s mode contains no entries. Please check selection_list option and make sure the file contains a  list. Run %s -h for help." % (selection_list, cter_mode_name, program_name), getframeinfo(currentframe()))
						break
				else:
					print(" ")
					print(("Processing a single micorgprah: %s..." % (selection_list)))
					selected_mic_path_list = [selection_list]
			
				selected_mic_directory = os.path.dirname(selected_mic_path_list[0])
				if selected_mic_directory != "":
					print(("    NOTE: Program disregards the directory paths in the selection list (%s)." % (selected_mic_directory)))
		
			# Register  id substrings to the global entry dictionary
			for selected_mic_path in selected_mic_path_list:
				# Find tail index of  id substring and extract the substring from the  name
				selected_mic_basename = os.path.basename(selected_mic_path)
				mic_id_substr_tail_idx = selected_mic_basename.index(mic_basename_tokens[1])
				mic_id_substr = selected_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
				if not mic_id_substr in global_entry_dict:
					# print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from selected_mic_path_list " % (mic_id_substr))
					global_entry_dict[mic_id_substr] = {}
				global_entry_dict[mic_id_substr][subkey_selected_mic_basename] = selected_mic_basename
			
			# --------------------------------------------------------------------------------
			# Clean up variables related to registration to the global entry dictionary
			# --------------------------------------------------------------------------------
			del mic_basename_tokens
			del mic_id_substr_head_idx
		
			# --------------------------------------------------------------------------------
			# Create the list containing only valid  id substrings
			# --------------------------------------------------------------------------------
			# Prepare lists to keep track of invalid (rejected) s 
			no_input_mic_id_substr_list = []
		
			print(" ")
			print("Checking the input datasets consistency...")
		
			# Loop over substring id list
			for mic_id_substr in global_entry_dict:
				mic_id_entry = global_entry_dict[mic_id_substr]
			
				warinnig_messages = []
				# selected  basename must have been registed always .
				if subkey_selected_mic_basename in mic_id_entry: 
					# Check if associated input  exists
					if not subkey_input_mic_path in mic_id_entry:
						input_mic_path = mic_pattern.replace("*", mic_id_substr)
						warinnig_messages.append("    associated input  %s." % (input_mic_path))
						no_input_mic_id_substr_list.append(mic_id_substr)
				
					if len(warinnig_messages) > 0:
						print(("WARNING!!! Micrograph ID %s does not have:" % (mic_id_substr)))
						for warinnig_message in warinnig_messages:
							print(warinnig_message)
						print("    Ignores this as an invalid entry.")
					else:
						# print("MRK_DEBUG: adding mic_id_substr := ", mic_id_substr)
						valid_mic_id_substr_list.append(mic_id_substr)
				# 	# This entry is not in the selection list. Do nothing
			
			# Check the input dataset consistency and save the result to a text file, if necessary.
			if check_consistency:
				# Create output directory
				os.mkdir(output_directory)
			
				# Open the consistency check file
				inconsist_mic_list_path = os.path.join(output_directory,"inconsist_mic_id_file.txt")
				print(" ")
				print(("Generating the input datasets consistency report in %s..." % (inconsist_mic_list_path)))
				inconsist_mic_list_file = open(inconsist_mic_list_path, "w")
				inconsist_mic_list_file.write("# The information about inconsistent  IDs\n")
				# Loop over substring id list
				for mic_id_substr in global_entry_dict:
					mic_id_entry = global_entry_dict[mic_id_substr]
				
					consistency_messages = []
					# Check if associated input  path exists
					if not subkey_input_mic_path in mic_id_entry:
						input_mic_path = mic_pattern.replace("*", mic_id_substr)
						consistency_messages.append("    associated input  %s." % (input_mic_path))
				
					# Check if associated selected  basename exists
					if not subkey_selected_mic_basename in mic_id_entry:
						input_mic_path = mic_pattern.replace("*", mic_id_substr)
						consistency_messages.append("    associated selected  %s." % (input_mic_path))
				
					if len(consistency_messages) > 0:
						inconsist_mic_list_file.write("Micrograph ID %s does not have:\n" % (mic_id_substr))
						for consistency_message in consistency_messages:
							inconsist_mic_list_file.write(consistency_message)
							inconsist_mic_list_file.write("\n")
			
				# Close the consistency check file, if necessary
				inconsist_mic_list_file.flush()
				inconsist_mic_list_file.close()
			
			# Since mic_id_substr is once stored as the key of global_entry_dict and extracted with the key order
			# we need sort the valid_mic_id_substr_list here
			if debug_mode: print(("BEFORE SORT: valid_mic_id_substr_list := ", valid_mic_id_substr_list))
			valid_mic_id_substr_list.sort(key=str.lower) # Sort list of  IDs using case insensitive string comparison
			if debug_mode: print(("AFTER SORT: valid_mic_id_substr_list := ", valid_mic_id_substr_list))
			
			# --------------------------------------------------------------------------------
			# Print out the summary of input consistency
			# --------------------------------------------------------------------------------
			print(" ")
			print("Summary of dataset consistency check...")
			print(("  Detected  IDs               : %6d" % (len(global_entry_dict))))
			print(("  Entries in input directory  : %6d" % (len(input_mic_path_list))))
			print(("  Entries in selection list   : %6d" % (len(selected_mic_path_list))))
			print(("  Rejected by no input        : %6d" % (len(no_input_mic_id_substr_list))))
			print(("  Valid Entries               : %6d" % (len(valid_mic_id_substr_list))))
			
			# --------------------------------------------------------------------------------
			# Check MPI error condition
			# --------------------------------------------------------------------------------
			if len(valid_mic_id_substr_list) < n_mpi_procs:
				error_status = ("Number of MPI processes (%d) supplied by --np in mpirun cannot be greater than %d (number of valid s that satisfy all criteria to be processed). Run %s -h for help." % (n_mpi_procs, len(valid_mic_id_substr_list, program_name)), getframeinfo(currentframe()))
				break
			
			# --------------------------------------------------------------------------------
			# Create input file path list
			# --------------------------------------------------------------------------------
			for mic_id_substr in valid_mic_id_substr_list:
				mic_path = global_entry_dict[mic_id_substr][subkey_input_mic_path]
				input_file_path_list.append(mic_path)
			
			# --------------------------------------------------------------------------------
			# Clean up variables related to tracking of invalid (rejected) s 
			# --------------------------------------------------------------------------------
			del input_mic_path_list
			del selected_mic_path_list
			del no_input_mic_id_substr_list
			
			break
		
		
		# --------------------------------------------------------------------------------
		# Clean up the variables for all sections in  mode case 
		# --------------------------------------------------------------------------------
		del mic_pattern
		del mic_basename_pattern
		del global_entry_dict
		del subkey_input_mic_path
		del subkey_selected_mic_basename
		del valid_mic_id_substr_list
		
		# --------------------------------------------------------------------------------
		# Print all error messages and abort the process if necessary.
		# --------------------------------------------------------------------------------
		# NOTE: Toshio Moriya 2016/11/15
		# The following function takes care of the case when an if-statement uses break for occurence of an error.
		# However, more elegant way is to use 'exception' statement of exception mechanism...
		# 
		if_error_then_all_processes_exit_program(error_status)
		
	else:
		input_file_path_list.append(input_image_path)
	
	if RUNNING_UNDER_MPI:
		# Wait for main mpi process to create the input file path list
		mpi_barrier(MPI_COMM_WORLD)
		
		# All mpi processes should know input file path list
		input_file_path_list = wrap_mpi_bcast(input_file_path_list, main_mpi_proc)

	
	# ====================================================================================
	# Prepare input file path(s)
	# ====================================================================================
	# 
	# NOTE: 2016/03/17 Toshio Moriya
	# From here on, stack (and namics) will be used to distinguish stack mode and  mode.
	# However, a single input_file_path_list should be sufficient since we already know the mode.
	# Let's consider this refactoring in the future.
	# 
	stack = None # (particle) stack file name: if it is not None, cter runs with stack mode. Otherwise, runs with  mode
	namics = []  #  file name list
	if not stack_mode:
		namics = input_file_path_list
		if debug_mode: print(("BEFORE SORT: namics := ", namics))
		namics.sort(key=str.lower) # Sort list of s using case insensitive string comparison
		if debug_mode: print(("AFTER SORT: namics := ", namics))
	else:
		stack = input_file_path_list[0]
	
	del input_file_path_list # Don't need this anymore
	
	# Make output directory
	outpwrot = "%s/pwrot" % (output_directory)
	if stack == None: 
		outmicthumb = "%s/micthumb" % (output_directory)
	if debug_mode:  
		outravg = "%s/ravg" % (output_directory)
	if my_mpi_proc_id == main_mpi_proc:
		# Make output directory
		if not os.path.exists(output_directory):
			os.mkdir(output_directory)
		os.mkdir(outpwrot)
		if stack == None: 
			os.mkdir(outmicthumb)
		if debug_mode:
			os.mkdir(outravg)
	
	if RUNNING_UNDER_MPI:
		# Make all mpi processes wait for main mpi process to create output directory
		mpi_barrier(MPI_COMM_WORLD)
	
	# Set up loop variables depending on the cter mode
	if stack == None:
		if RUNNING_UNDER_MPI:
			set_start, set_end = MPI_start_end(len(namics), n_mpi_procs, my_mpi_proc_id)
		else:
			set_start = 0
			set_end = len(namics)
	else: 
		pw2 = []
		data = EMData.read_images(stack)
		nima = len(data)
		for i in range(nima):
			pw2.append(periodogram(data[i]))
		wn = pw2[0].get_xsize()
		set_start = 0
		set_end = 1
	
	# Set up progress message
	if my_mpi_proc_id == main_mpi_proc:
		print(" ")
		print("Estimating CTF parameters...")
		if stack == None:
			print("  Micrographs processed by main process (including percent of progress):")
			progress_percent_step = (set_end - set_start)/100.0 # the number of micrograms for main mpi processer divided by 100
	
	totresi = []
	missing_img_names = []
	rejected_img_names = []
	for ifi in range(set_start, set_end):
		if stack == None:
			img_type = "Micrograph"
			img_name = namics[ifi]
			img_basename_root = os.path.splitext(os.path.basename(img_name))[0]

			if my_mpi_proc_id == main_mpi_proc:
				print(("    Processing %s ---> %6.2f%%" % (img_name, (ifi - set_start) / progress_percent_step)))

			if not os.path.exists(img_name):
				missing_img_names.append(img_name)
				print("    %s %s: Can not find this file. Skipping the estimation and CTF parameters are not stored..." % (img_type, img_name))
				continue
			mic = get_im(img_name)
			try:
				pw2 = tilemic(mic, win_size = wn, overlp_x = overlap_x, overlp_y = overlap_y, edge_x = edge_x, edge_y = edge_y)
			except:
				print("MRK_DEBUG: tilemic() in cter_vpp() raised an exception. The micrographs {} might have a problem. Please check it and remove it if necessary.".format(img_name))
				raise
			del mic

		else:
			img_type = "Stack"
			img_name = stack

			numFM = EMUtil.get_image_count(img_name)
			pw2 = [None]*numFM
			for i in range(numFM):
				pw2[i] = periodogram(get_im(img_name,i))

		for i in range(len(pw2)):
			pw2[i] = square_root(pw2[i])

		if debug_mode: print("    %s %s: Process %04d started the processing. Detected %d image(s) in this %s file." % (img_type, img_name, ifi, img_type.lower()))

		if db_check_dict(img_name) == False:
			img_basename_root = os.path.splitext(os.path.basename(img_name))[0]
		else:
			path, dictname, keys = db_parse_path(img_name)
			img_basename_root = dictname

		#  VPP code starts here  03/08/2017
		nimi = len(pw2)
		adefocus = [0.0] * kboot
		aampcont = [0.0] * kboot
		aamplitu = [0.0] * kboot
		aangle   = [0.0] * kboot
		
		allroo = []
		for imi in range(nimi):
			allroo.append(rot_avg_table(pw2[imi]))
		lenroo = len(allroo[0])

		for nboot in range(kboot):
			#at = time()
			if(nboot == 0): boot = list(range(nimi))
			else:
				from random import randint
				for imi in range(nimi): boot[imi] = randint(0, nimi - 1)
			qa = model_blank(wn, wn)
			roo  = np.zeros(lenroo, np.float32)
			sroo = np.zeros(lenroo, np.float32)
			aroo = np.zeros(lenroo, np.float32)

			for imi in range(nimi):
				Util.add_img(qa, pw2[boot[imi]])
				temp1 = np.array(allroo[boot[imi]])
				roo += temp1
				temp2 = movingaverage(temp1, 10)
				aroo += temp2
				sroo += temp2**2
			sroo[0] = sroo[1]
			aroo[0] = aroo[1]
			sroo = (sroo-aroo**2 / nimi) / nimi
			aroo /= nimi
			roo  /= nimi
			qa   /= nimi

			if f_start < 0:
				#  Find a break point
				bp = 1.e23
				for i in range(5, lenroo - 5):
					#t1 = linreg(sroo[:i])
					#t2 = linreg(sroo[i:])
					#tt = t1[1][0] + t2[1][0]
					xtt = np.array(list(range(i)), np.float32)
					zet = np.poly1d(np.polyfit(xtt,sroo[:i], 2))
					t1 = sum((sroo[:i] - zet(xtt))**2)
					xtt = np.array(list(range(i, lenroo)), np.float32)
					zet = np.poly1d(np.polyfit(xtt, sroo[i:], 2) )
					tt = t1 + sum((sroo[i:] - zet(xtt))**2)
					if tt < bp:
						bp = tt
						istart = i
				f_start = istart / (pixel_size * wn)
			"""
			hi = hist_list(sroo,2)
			# hi[0][1] is the threshold
			for i in xrange(1,len(sroo)):
				if(sroo[i] < hi[0][1]):
					istart = i
					break
			"""
			#write_text_file([roo.tolist(),aroo.tolist(),sroo.tolist()], "sroo%03d.txt"%ifi)
			rooc = roo.tolist()

			defc, ampcont, subpw, baseline, envelope, istart, istop = defocusgett_vpp(rooc, wn, voltage = voltage, Pixel_size = pixel_size, Cs = Cs, \
										f_start = f_start, f_stop = f_stop, vpp_options = vpp_options, nr1 = 3, nr2 = 6, parent = None, DEBug = debug_mode)
			if debug_mode:
				print("  RESULT 1 %s" % (img_name), nboot, defc, ampcont, istart, istop)#, (time()-at)/60.
				'''
				freq = range(len(subpw))
				for i in xrange(len(freq)):  freq[i] = float(i) / wn / pixel_size
				#write_text_file([freq, subpw.tolist(), ctf2, envelope.tolist(), baseline.tolist()], "%s/ravg%05d.txt" % (output_directory, ifi))
				#fou = os.path.join(outravg, "%s_ravg_%02d.txt" % (img_basename_root, nboot))
				fou = os.path.join(".", "%s_ravg_%02d.txt" % (img_basename_root, nboot))
				write_text_file([freq, subpw.tolist(), ctf2, envelope.tolist(), baseline.tolist()], fou)
				'''
			#mpi_barrier(MPI_COMM_WORLD)
			
			#exit()
			bg = baseline.tolist()
			en = envelope.tolist()

			bckg = model_blank(wn, wn, 1, 1)
			envl = model_blank(wn, wn, 1, 1)

			from math import sqrt
			nc = wn // 2
			ne = istop
			ns = istart - 1
			bg.append(bg[-1])
			en.append(en[-1])
			for i in range(wn):
				for j in range(wn):
					r = sqrt((i - nc)**2 + (j - nc)**2)
					ir = int(r+0.5)
					if( ir < nc ):
						# This is awkward but it is needed for storing of results
						dr = r - ir
						bckg.set_value_at(i, j, (1. - dr) * bg[ir] + dr * bg[ir + 1] )
						if(ir > ns and ir < ne):
							dr = r - ir
							envl.set_value_at(i, j, (1. - dr) * en[ir] + dr * en[ir + 1] )

			qse = threshold((qa - bckg))/envl
			#print  "  fit1  ", nboot,(time()-at)/60.0
			#at = time()
			#(qse*mask).write_image("rs2.hdf")
			#bckg.write_image("bckg.hdf")
			#envl.write_image("envl.hdf")
			#qse.write_image("qse.hdf")
			del envl, baseline, envelope
			#exit()
			##  SIMULATION
			#bang = 0.7
			#qse = ctf2_rimg(wn, generate_ctf([defc,Cs,voltage,pixel_size,0.0,wgh, bang, 37.0]) )
			#qse.write_image("rs3.hdf")
			#at = time()
			defc, ampcont, astamp, astang, score =  defocusgett_vpp2(qse, wn, defc, ampcont, voltage = voltage, Pixel_size = pixel_size, Cs = Cs, i_start=istart, i_stop=istop, parent = None, DEBug = debug_mode)
			if debug_mode:
				print("  RESULT 2 %s" % (img_name), nboot, defc, ampcont, astamp, astang, score)#, (time()-at)/60.
				'''
				freq = range(len(subpw))
				for i in xrange(len(freq)):  freq[i] = float(i) / wn / pixel_size
				#write_text_file([freq, subpw.tolist(), ctf2, envelope.tolist(), baseline.tolist()], "%s/ravg%05d.txt" % (output_directory, ifi))
				#fou = os.path.join(outravg, "%s_ravg_%02d.txt" % (img_basename_root, nboot))
				fou = os.path.join(".", "%s_ravg22_%02d.txt" % (img_basename_root, nboot))
				write_text_file([freq, subpw.tolist(), ctf2, envelope.tolist(), baseline.tolist()], fou)
				'''
			adefocus[nboot] = defc
			aampcont[nboot] = ampcont
			aamplitu[nboot] = astamp
			aangle[nboot]   = astang
			#from sys import exit
			#exit()
			#from morphology import ctf_rimg, ctf_1d
			#cq = ctf_1d(wn, generate_ctf([defc, Cs, voltage, pixel_size, 0.0, ampcont,astamp,astang]), doabs = True)[20:150]
			#write_text_file([subpw[20:150],cq],"pwds%02d.txt"%nboot)
			#print  "  fit2  ", nboot,(time()-at)/60.0
		
		#print  "  xxx2  ", nboot,(time()-at)/60.0

		ad1, ad2, ad3, ad4 = table_stat(adefocus) # return values: average, variance, minimum, maximum
		if ad2 <= 0.0:
			print(("    %s %s: Detected the variance less than zero (defocus statics: avg = %f, var = %f, min = %f, max = %f)." % (img_type, img_name, ad1, ad2, ad3, ad4)))
			print(("           The program ignores this estimate..."))
			continue
		
		reject = []
		thr = 3 * sqrt(ad2)
		for i in range(len(adefocus)):
			if( abs(adefocus[i] - ad1) > thr ):
				print(("    %s %s: Rejected an outlier defocus estimate (defocus = %f, average defocus = %f, threshold = %f)." % (img_type, img_name, adefocus[i], ad1, thr)))
				reject.append(i)

		if(len(reject) > 0):
			print(("    %s %s: Total number of rejects %s" % (img_type, img_name, len(reject))))
			for i in range(len(reject) - 1, -1, -1):
				del adefocus[i]
				del aampcont[i]
				del aamplitu[i]
				del aangle[i]

		#print  "  xxx2  ", nboot,(time()-at)/60.0
		if(len(adefocus) < 2):
			print(("    %s %s: After rejection of outliers, there is too few estimated defocus values. Skipping the estimation and CTF parameters are not stored..." % (img_type, img_name)))
		else:
			#print "adefocus",adefocus
			#print  "aamplitu",aamplitu
			#print "aangle",aangle
			ad1, ad2, ad3, ad4 = table_stat(adefocus)
			#  compute statistics of ampcont using phase shifts instead
			ed1, ed2 = angle_ave([2*ampcont2angle(q) for q in aampcont]) # Have to use this trick as the function works for range [0,360]
			ed1/=2
			ed2/=2
			bd1, bd2, bd3, bd4 = table_stat(aamplitu)
			cd1,cd2 = angle_ave([2*q for q in aangle])  # Have to use this trick as the function works for range [0,360]
			cd1/=2
			cd2/=2
			temp = 0.0
			stdavad1 = np.sqrt(kboot * max(0.0, ad2))
			ed2 *= np.sqrt(kboot)
			stdavbd1 = np.sqrt(kboot * max(0.0, bd2))
			cd2 *= np.sqrt(kboot)
			
			# Adjust value ranges of astig. amp. and angle.
			if( bd1 < 0.0 ):
				bd1 = -bd1
				cd1 += 90.0
			cd1 = cd1%180
			
			#  SANITY CHECK, do not produce anything if defocus abd astigmatism amplitude are out of whack
			reject_img_messages = []
			try:
				pwrot2 = rotavg_ctf( model_blank(wn, wn), ad1, Cs, voltage, pixel_size, bd1, cd1)
			except:
				reject_img_messages.append("    - Astigmatism amplitude (%f) is larger than defocus (%f) or defocus (%f) is negative." % (bd1, ad1, ad1))
			"""
			valid_min_defocus = 0.05
			if ad1 < valid_min_defocus:
				reject_img_messages.append("    - Defocus (%f) is smaller than valid minimum value (%f)." % (ad1, valid_min_defocus))
			"""

			if len(reject_img_messages) > 0:
				rejected_img_names.append(img_name)
				print("    %s %s: Rejected the CTF estimate - " % (img_type, img_name), ad1, Cs, voltage, pixel_size, angle2ampcont(ed1), bd1, cd1, "(def, Cs, vol, apix, amp_contrast, astig_amp, astig_angle)")
				print("    %s %s: because... " % (img_type, img_name))
				for reject_img_message in reject_img_messages:
					print(reject_img_message)
				print("    %s %s: Skipping the estimation and CTF parameters are not stored..." % (img_type, img_name))
			else:
				#  Estimate the point at which (sum_errordz ctf_1(dz+errordz))^2 falls to 0.5
				import random as rqt
				#  NOW WE SWITCH PHASE SHIFT TO AMPLITUDE CONTRAST
				ed1 = angle2ampcont(ed1)
				#at = time()
				supe = model_blank(wn, wn)
				niter = 1000
				for it in range(niter):
					Util.add_img(supe, Util.ctf_rimg(wn, wn, 1, ad1 + rqt.gauss(0.0,stdavad1), pixel_size, voltage, Cs, 0.0, ed1 + angle2ampcont(rqt.gauss(0.0,ed2)), bd1 + rqt.gauss(0.0,stdavbd1), cd1 + rqt.gauss(0.0,cd2), 1))
				ni = wn // 2
				supe /= niter
				pwrot2 = rotavg_ctf(supe, ad1, Cs, voltage, pixel_size, bd1, cd1)
				for i in range(ni):  pwrot2[i] = pwrot2[i]**2
				ibec = 0
				for it in range(ni - 1, 0, -1):
					if pwrot2[it] > 0.5 :
						ibec = it
						break
				#print  "  uuu1 ",(time()-at)/60.
				#  NOW WE SWITCH PHASE SHIFT TO AMPLITUDE CONTRAST
				ed2 = angle2ampcont(ed2)
				if(ed2<0.0): ed2 = 180.0 - ed2

				from morphology import ctf_1d
				ct = generate_ctf([ad1, Cs, voltage, pixel_size, temp, ed1, 0.0, 0.0])
				cq = ctf_1d(wn, ct)
				#at = time()
				supe = [0.0] * ni
				niter = 1000
				for i in range(niter):
					cq = generate_ctf([ad1 + rqt.gauss(0.0,stdavad1), Cs, voltage, pixel_size, 0.0, ed1, 0.0, 0.0])
					ci = ctf_1d(wn, cq)[:ni]
					for l in range(ni):  supe[l] +=ci[l]
				
				for l in range(ni):  supe[l] = (supe[l] / niter)**2
				
				ib1 = 0
				for it in range(ni - 1, 0, -1):
					if supe[it] > 0.5:
						ib1 = it
						break
				ibec = ibec / (pixel_size * wn)  #  with astigmatism
				ib1  = ib1 / (pixel_size * wn)   #  no astigmatism
				#print  " error est  ",(time()-at)/60.0
				#from utilities import write_text_file
				#write_text_file([range(ni), supe[:ni],pwrot2[:ni]],"fifi.txt")
				
				# Compute defocus CV and astig. amp. CV (CV: coefficient of variation; ratio of error (SD) relative to average (mean))
				#  I blocked them, make no sense here.
				#if ad1 < max(0.0, valid_min_defocus): ERROR("Logical Error: Encountered unexpected defocus value (%f). Consult with the developer." % (ad1), "%s in %s" % (__name__, os.path.basename(__file__))) # MRK_ASSERT
				#if stdavad1 < 0.0: ERROR("Logical Error: Encountered unexpected defocus SD value (%f). Consult with the developer." % (stdavad1), "%s in %s" % (__name__, os.path.basename(__file__))) # MRK_ASSERT
				cvavad1 = stdavad1 / ad1 * 100 # use percentage
				
				#if bd1 < 0.0: ERROR("Logical Error: Encountered unexpected astig. amp. value (%f). Consult with the developer." % (bd1), "%s in %s" % (__name__, os.path.basename(__file__))) # MRK_ASSERT
				if stdavbd1 < 0.0: ERROR("Logical Error: Encountered unexpected astig. amp. SD value (%f). Consult with the developer." % (stdavbd1), "%s in %s" % (__name__, os.path.basename(__file__))) # MRK_ASSERT

				bd1 = max(bd1, 1.0e-15)
				cvavbd1 = stdavbd1 / bd1 * 100 # use percentage
				
				# Compute CTF limit (theoretical resolution limit based on the oscillations of CTF) 
				# For output, use ctflim (relative frequency limit [1/A]), not ctflim_abs (absolute frequency limit)
				# 
				# NOTE: 2016/03/23 Toshio Moriya
				# xr is limiting frequency [1/A]. Max is Nyquist frequency = 1.0/(2*apix[A/pixel]). <UNIT: [1/(A/pixel)/[pixel])] => [(pixel)/(A*pixel] => [1/A]>
				# 1.0/xr is limiting period (Angstrom resolution) [A]. Min is Nyquist period = (2*apix[A/pixel]). <UNIT: [1/(1/A)] = [A]>
				# fwpix is width of Fourier pixel [pixel/A] := 1.0[pixel]/(2*apix[A/pixel])/box_half[pixel] = 1[pixel]/fullsize[A]). <UNIT: [pixel/(A/pixel)/(pixel)] = [pixel*(pixel/A)*(1/pixel) = [pixel/A]>
				# int(xr/fwpix+0.5) is limiting_absolute_frequency [1/pixel]. <Unit:[(1/A)/(pixel/A)] = [(1/A)*(A/pixel)] = [1/pixel]>
				# return  int(xr/fwpix+0.5),xr, which is limiting_abs_frequency [1/pixel], and Limiting_frequency[1/A]
				#
				ctflim_abs, ctflim = ctflimit(wn, ad1, Cs, voltage, pixel_size)

				"""
				for i in xrange(len(ssubroo)):
					asubroo[i] /= kboot
					ssubroo[i]  = sqrt(max(0.0, ssubroo[i]-kboot*asubroo[i]**2)/kboot)
					sen[i]     /= kboot
				"""
				#print  " error est2  ",(time()-at)/60.0
				#print " ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1 ",ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1
				lnsb = len(subpw)
				try:		crot1 = rotavg_ctf(ctf_rimg(wn, generate_ctf([ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1]), sign=0), ad1, Cs, voltage, pixel_size, 0.0, 0.0)[:lnsb]
				except:		crot1 = [0.0] * lnsb
				try:		pwrot1 = rotavg_ctf(threshold(qa - bckg), ad1, Cs, voltage, pixel_size, 0.0, 0.0)[:lnsb]
				except:		pwrot1 = [0.0] * lnsb
				try:		crot2 = rotavg_ctf(ctf_rimg(wn, generate_ctf([ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1]), sign=0), ad1, Cs, voltage, pixel_size, bd1, cd1)[:lnsb]
				except:		crot2 = [0.0] * lnsb
				try:		pwrot2 = rotavg_ctf(threshold(qa - bckg), ad1, Cs, voltage, pixel_size, bd1, cd1)[:lnsb]
				except:		pwrot2 = [0.0] * lnsb
				#  #1 - rotational averages without astigmatism, #2 - with astigmatism
				lnsb = min(lnsb,len(crot2),len(pwrot1),len(crot2),len(pwrot2))
				write_text_file([list(range(lnsb)), [float(i)/wn/pixel_size for i in range(lnsb)], pwrot1[:lnsb], crot1[:lnsb], pwrot2[:lnsb], crot2[:lnsb]], os.path.join(outpwrot, "%s_rotinf.txt"%(img_basename_root)))

				#
				#print  " error est3  ",(time()-at)/60.0
				# NOTE: 2016/03/23 Toshio Moriya
				# Compute mean of extrema differences (differences at peak & trough) between 
				# (1) experimental rotational average with astigmatism (pwrot2)
				# (2) experimental rotational average without astigmatism (pwrot1), and
				# as a indication of goodness of astigmatism estimation by cter.
				# The peak & trough detection uses fitted rotational average with astigmatism (crot2) 
				# Start from 1st trough while ignoring 1st peak.
				# End at astigmatism frequency limit.
				# 
				'''
				is_peak_target = True
				pre_crot2_val = crot2[0]
				extremum_counts = 0
				extremum_diff_sum = 0
				for i in xrange(1, len(crot2)):
					cur_crot2_val = crot2[i]
					if is_peak_target == True and pre_crot2_val > cur_crot2_val:
						# peak search state
						extremum_i = i - 1
						extremum_counts += 1
						extremum_diff_sum += pwrot2[extremum_i] - pwrot1[extremum_i] # This should be positive if astigmatism estimation is good
						# print "MRK_DEBUG: Peak Search  : extremum_i = %03d, freq[extremum_i] = %12.5g, extremum_counts = %03d, (pwrot2[extremum_i] - pwrot1[extremum_i]) = %12.5g, extremum_diff_sum = %12.5g " % (extremum_i, freq[extremum_i] , extremum_counts, (pwrot2[extremum_i] - pwrot1[extremum_i]), extremum_diff_sum)
						is_peak_target = False
					elif is_peak_target == False and pre_crot2_val < cur_crot2_val:
						# trough search state
						extremum_i = i - 1
						extremum_counts += 1
						extremum_diff_sum += pwrot1[extremum_i] - pwrot2[extremum_i] # This should be positive if astigmatism estimation is good
						# print "MRK_DEBUG: Trough Search: extremum_i = %03d, freq[extremum_i] = %12.5g, extremum_counts = %03d, (pwrot1[extremum_i] - pwrot2[extremum_i]) = %12.5g, extremum_diff_sum = %12.5g " % (extremum_i, freq[extremum_i] , extremum_counts, (pwrot1[extremum_i] - pwrot2[extremum_i]), extremum_diff_sum)
						is_peak_target = True
					pre_crot2_val = cur_crot2_val
				'''
#				#if extremum_counts == 0: ERROR("Logical Error: Encountered unexpected zero extremum counts. Consult with the developer." % (bd1), "%s in %s" % (__name__, os.path.basename(__file__))) # MRK_ASSERT
#				extremum_diff_avg = 1.1#extremum_diff_sum / extremum_counts
				
#				#print "MRK_DEBUG: extremum_avg = %12.5g, extremum_diff_sum = %12.5g, extremum_counts = %03d," % (extremum_avg, extremum_diff_sum, extremum_counts)
#				#print "MRK_DEBUG: extremum_diff_avg = %12.5g, extremum_diff_sum = %12.5g, extremum_counts = %03d," % (extremum_diff_avg, extremum_diff_sum, extremum_counts)
				
#				if stack == None:     cmd = "echo " + "    " + namics[ifi] + "  >>  " + fou
#				else:                 cmd = "echo " + "    " + "  >>  " + fou
#				os.system(cmd)
				cvavbd1 = stdavbd1 / bd1 * 100 # use percentage
				
				max_freq = 0.5/pixel_size # dummy value for maximum frequency. set to Nyquist frequency for now. let's add the implementation in near future (Toshio 2017/12/06)
				reserved = 0.0            # dummy value for reserved spot, which might be used for parameter of external programs (e.g. CTFFIND4, GCTF, and etc.)
				# wgh                     # constant amplitude contrast provided by user (default 10%)
				phase_shift = ampcont2angle(ed1) - ampcont2angle(wgh) # Volta phase shift [deg] = total amplitude contrast phase shift [deg] (ed1) -  constant amplitude contrast phase shift [deg]; ed1 is boot strap average of total amplitude contrast [%]
				
				if debug_mode: print(("    %s %s: Process %04d finished the processing. Estimated CTF parmaters are stored in %s." % (img_type, img_name, ifi, os.path.join(output_directory, "partres.txt"))))
#				if debug_mode: print((ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1, stdavad1, ed2, stdavbd1, cd2, cvavad1, cvavbd1, extremum_diff_avg, ib1, ibec, ctflim))
#				totresi.append([img_name, ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1, stdavad1, ed2, stdavbd1, cd2, cvavad1, cvavbd1, extremum_diff_avg, ib1, ibec, ctflim])
				if debug_mode: print((ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1, stdavad1, ed2, stdavbd1, cd2, cvavad1, cvavbd1, ib1, ibec, ctflim, max_freq, reserved, wgh, phase_shift))
				totresi.append([img_name, ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1, stdavad1, ed2, stdavbd1, cd2, cvavad1, cvavbd1, ib1, ibec, ctflim, max_freq, reserved, wgh, phase_shift])
				
#				if stack == None:
#					print  namics[ifi], ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec
#				else:
#					print               ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec
#				if stack == None:
#					totresi.append( [ namics[ifi], ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec])
#				else:
#					totresi.append( [ 0, ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec])
#				#if ifi == 4 : break
				#print  " error est4  ",(time()-at)/60.0

		if stack == None:
			img_mic = get_im(namics[ifi])
			# create  thumbnail
			nx = img_mic.get_xsize()
			if nx > 512:
				img_micthumb = resample(img_mic, 512.0/nx)
			else:
				img_micthumb = img_mic
			img_basename_root = os.path.splitext(os.path.basename(img_name))[0]
			fou = os.path.join(outmicthumb, "%s_thumb.hdf" % (img_basename_root))
			img_micthumb.write_image(fou)

	if RUNNING_UNDER_MPI:
		from utilities import wrap_mpi_gatherv
		totresi = wrap_mpi_gatherv(totresi, 0, MPI_COMM_WORLD)
		missing_img_names = wrap_mpi_gatherv(missing_img_names, 0, MPI_COMM_WORLD)
		rejected_img_names = wrap_mpi_gatherv(rejected_img_names, 0, MPI_COMM_WORLD)
	
	if my_mpi_proc_id == main_mpi_proc:
		outf = open(os.path.join(output_directory, "partres.txt"), "w")
		for i in range(len(totresi)):
			for k in range(1, len(totresi[i])):
				outf.write("  %12.5g" % totresi[i][k])
			outf.write("  %s\n" % totresi[i][0])
		outf.close()
		
		print(" ")
		print(("Summary of %s processing..." % (img_type.lower())))
		missing_counts = len(missing_img_names)
		print(("  Missing  : %d" % (missing_counts)))
		if missing_counts > 0:
			outfile_path = os.path.join(output_directory, "missing_%s_list.txt" % (img_type.lower()))
			print(("    Saving list of missing in %s..." % (outfile_path)))
			outf = open(outfile_path, "w")
			for missing_img_name in missing_img_names:
				outf.write("%s\n" % missing_img_name)
			outf.close()
		
		rejected_counts = len(rejected_img_names)
		print(("  Rejected : %d" % (rejected_counts)))
		if rejected_counts > 0:
			outfile_path = os.path.join(output_directory, "rejected_%s_list.txt" % (img_type.lower()))
			print(("    Saving list of rejected in %s..." % (outfile_path)))
			outf = open(outfile_path, "w")
			for rejected_img_name in rejected_img_names:
				outf.write("%s\n" % rejected_img_name)
			outf.close()
	
	if cter_mode_idx == idx_cter_mode_stack:
		return totresi[0][1], totresi[0][7], totresi[0][8], totresi[0][9], totresi[0][10], totresi[0][11]
	
########################################
# functions used by cter_vpp
########################################

def defocusgett_vpp(roo, nx, voltage=300.0, Pixel_size=1.0, Cs=2.0, f_start=-1.0, f_stop=-1.0, vpp_options = [], nr1 = 3, nr2 = 6, parent=None, DEBug=False):
	"""
		1. Estimate envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get 
		   defocus which matches the extracted CTF imprints 
	"""
	from utilities  import generate_ctf, write_text_file
	import numpy as np
	from morphology import defocus_baseline_fit, simpw1d, angle2ampcont

	#print "CTF params:", voltage, Pixel_size, Cs, wgh, f_start, f_stop, round_off, nr1, nr2, parent

	if f_start == 0 : 	    i_start = 0
	else: 			        i_start = int(Pixel_size*nx*f_start+0.5)
	if f_stop <= f_start :
		i_stop  = len(roo)
		adjust_fstop = True
	else:
		i_stop  = min(len(roo), int(Pixel_size*nx*f_stop+0.5))
		adjust_fstop = False

	i_stop = min(nx//2-2, i_stop)  #  this is for resampling into polar

	nroo = len(roo)

	if DEBug:  print("f_start, i_start, f_stop, i_stop:", f_start, i_start, f_stop, i_stop-1)
	#TE  = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr1), 4)
	#baseline = defocus_baseline_fit(roo, i_start, i_stop, int(nr2), 3)
	baseline = defocus_baseline_fit(roo, i_start, nroo, int(nr2), 3)
	subpw = np.array(roo, np.float32) - baseline
	subpw[0] = subpw[1]
	#write_text_file([roo,baseline,subpw],"dbg.txt")
	#print "IN defocusgett  ",np.min(subpw),np.max(subpw)
	for i in range(len(subpw)):  subpw[i] = max(subpw[i],0.0)
	#print "IN defocusgett  ",np.min(subpw),np.max(subpw)
	#envelope = movingaverage(  subpw, nroo//8, 3)
	
	#envelope = defocus_baseline_fit(roo, i_start, nroo, int(nr2), 2) - baseline
	envelope = defocus_baseline_fit(roo, i_start, min(int(i_stop*1.45),nx//2-2), int(nr2), 2) - baseline
	#  Process envelope
	qm = np.max(envelope[5:])
	dc = max(np.min(envelope[5:]), qm/1000.)
	for i in range(len(envelope)):
		if(envelope[i]<dc): envelope[i] = qm
	
	#envelope = np.array([1.0]*len(subpw), np.float32)
	#write_text_file([roo,baseline,subpw,envelope],"dbgt.txt")

	#print "IN defocusgett  ",np.min(subpw),np.max(subpw),np.min(envelope)
	#envelope = np.ones(nroo, np.float32)
	defocus = 0.0
	ampcont = 0.0
	data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
	qm = 1.e23
	#toto = []
	#  vpp_options = [defocus_min,  defocus_max,  defocus_step,  phase_min,  phase_max,  phase_step]
	#  This is in degrees
	if( vpp_options[4] < vpp_options[3] ): vpp_options[4] += 180.0
	a = vpp_options[3]
	while( a<=vpp_options[4] ):
		data[7] = angle2ampcont(a%180.0)
		dc = vpp_options[0]
		while( dc <= vpp_options[1] ):
			qt = simpw1d_pap(dc, data)
			#toto.append([a,data[7],dc,qt])
			if(qt<qm):
				qm = qt
				defi = dc
				ampcont = data[7]
				#print  a,dc,qt
			dc += vpp_options[2]
		a += vpp_options[5]
	#'''
	if DEBug:
		from utilities import write_text_row
		#write_text_row(toto,"toto1.txt")
		data[7] = ampcont
		print(" >>>>>>>>>  ",defi,data[7],ampcont2angle(data[7]),simpw1d_print(defi, data))#,generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont])
		#data[7]=10.
		#defi = 4.5
		#print " >>>>>>>>>  ",defi,data[7],simpw1d_print(defi, data)#,generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont])
		#def1 = defi
	#exit()
	#'''
	#ctf2 = ctf_1d(nx, generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont]), doabs= True)
	'''
	from utilities import write_text_file
	foki = subpw.tolist()
	write_text_file([foki,ctf2[:len(foki)]],"toto1.txt")
	'''
	return defi, ampcont, subpw.tolist(), baseline, envelope, i_start, i_stop  #, ctf2



def defocusgett_vpp2(qse, wn, xdefc, xampcont, voltage=300.0, Pixel_size=1.0, Cs=2.0, i_start=0, i_stop=0, parent=None, DEBug=False):
	"""
		1. Estimate envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get
		   defocus which matches the extracted CTF imprints
		Switching to phase shift altogether 07/07/2017
	"""
	#from utilities  import generate_ctf
	#import numpy as np
	from utilities import amoeba
	from alignment import Numrinit, ringwe

	cnx = wn // 2 + 1
	cny = cnx
	mode = "H"
	numr = Numrinit(i_start, i_stop-1, 1, mode)
	wr = ringwe(numr, mode)

	crefim = Util.Polar2Dm(qse, cnx, cny, numr, mode)
	Util.Frngs(crefim, numr)
	Util.Applyws(crefim, numr, wr)

	bdef = 0.
	bphs = 0.0  	#  phase_shift
	bamp = 0.0      #  initial astigmatism amplitude
	bang = 0.0      #  initial astigmatism angle
	astdata = [crefim, numr, wn, bdef, Cs, voltage, Pixel_size, bphs, bamp, bang]
	initial_ast_ang = 0.0
	xphase_shift = ampcont2angle(xampcont)
	dama = amoeba([xdefc,xphase_shift,initial_ast_ang], [0.1, 2.0, 0.05], fupw_vpp, 1.e-4, 1.e-4, 500, astdata)
	qma = -dama[-2]
	if DEBug: print(" amoeba  %7.2f  %7.2f  %12.6g  %12.6g"%(dama[0][0],dama[0][1]%180.0,dama[0][2],qma))
	dpefi = dama[0][0]
	dphshift = dama[0][1]%180.0
	dastamp = dama[0][2]
	astdata = [crefim, numr, wn, dpefi, Cs, voltage, Pixel_size, dphshift, dastamp, bang]
	junk = fastigmatism3_vpp(dama[0][2], astdata)
	#  Corrected here PAP 06/16/2018
	dastang = astdata[9]


	'''

	dp = 1.0e23
	toto = []

	for aa in xrange(0,20,5):
		a = xampcont + aa - 10.
		print "  fdasfdsfa  ",a
		for i in xrange(0,2000,500):
			dc = xdefc + float(i-1000)/10000.0
			ju1 = dc # defocus
			ju2 = float(a) # amp contrast
			ju3 = 0.0  # astigma amp
			dama = amoeba([ju1,ju2,ju3], [0.1, 2.0, 0.05], fupw_vpp, 1.e-4, 1.e-4, 500, astdata)
			qma = -dama[-2]
			print  " amoeba  %7.2f  %7.2f  %12.6g  %12.6g"%(dama[0][0],dama[0][1],dama[0][2],qma)
			toto.append([dama[0][0],dama[0][1],dama[0][2],qma])
			if(qma<dp):
				dp = qma
				dpefi = dama[0][0]
				dpmpcont = dama[0][1]
				dastamp = dama[0][2]
				astdata = [crefim, numr, wn, dpefi, Cs, voltage, Pixel_size, dpmpcont, dastamp, bang]
				junk = fastigmatism3_vpp(dama[0][2], astdata)
				dastang = astdata[8]
				print " FOUND ANGLE",junk, qma, dpefi,dpmpcont,dastamp,dastang
		#from sys import exit
		#exit()
	'''
	if DEBug:
		#from utilities import write_text_row
		#write_text_row(toto,"toto1.txt")
		print(" repi3  ", dpefi, dphshift, dastamp, dastang, junk)

	return dpefi, angle2ampcont(dphshift), dastamp, dastang, qma#dp


def fupw_vpp(args, data):
	from morphology import fastigmatism3_vpp
	#  args = [defocus, phaseshift, astigma-amp]
	#                                   0       1     2   3     4    5         6          7     8     9 
	#            (astdata) =          [crefim, numr, wn, bdef, Cs, voltage, Pixel_size, bphs, bamp, bang]
	#
	#   [data[0], data[1], data[2], args[0], data[4], data[5], data[6], args[1], data[8], data[9]]
	#   [crefim,   numr,   wn, (args)defocus, Cs,   voltage, Pixel_size,(a)phshift, (a)astamp, ang, mask]
	#
	#print  " fuw_vpp           ",args[0],args[1],args[2]
	args[0] = max(min(args[0], 6.0), 0.01)
	args[1] = args[1]#  %180.0  #  Phase shift within valid range
	args[2] = max(min(args[2], 3.0), 0.0)
	#                        (a)astamp
	return fastigmatism3_vpp(args[2],[data[0], data[1], data[2], args[0], data[4], data[5], data[6], args[1], data[8], data[9]])


def fastigmatism3_vpp(amp, data):
	from morphology import ctf2_rimg
	from utilities  import generate_ctf
	from alignment  import ornq
	from math       import sqrt
	#  data[0] - crefim
	#  data[1] - numr
	#  data[2] - nx (image is square)
	#  data[8] - astigmatism amplitude
	#  data[9] - mask defining the region of interest
	#
	#      0        1          2       3        4       5         6         7      8        9 
	#   [data[0], data[1], data[2], args[0], data[4], data[5], data[6], args[1], data[8], data[9]]
	#   [crefim,   numr,   wn, (args)defocus, Cs,   voltage, Pixel_size,(a)phshift, (a)astamp, ang]
	#
	#  generate_ctf
	#      0      1    2       3       4        5        6                      7
	#  [defocus, cs, voltage, apix, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle]
	#  [ microns, mm, kV, Angstroms, A^2, microns, radians]

	cnx = data[2]//2+1
	#qt = 0.5*nx**2
	#B = 0.0
	pc = ctf_rimg(data[2], generate_ctf([data[3], data[4], data[5], data[6], 0.0, angle2ampcont(data[7]%180.0), amp, 0.0]), sign = 0. )
	#st = pc.cmp("dot", pc, dict(negative = 0, mask = data[10], normalize = 0))
	#Util.mul_scalar(pc, 1.0/st)
	ang, sxs, sys, mirror, peak = ornq_vpp(pc, data[0], [0.0,0.0], [0.0,0.0], 1, "H", data[1], cnx, cnx)
	#print  ang, sxs, sys, mirror, peak
	#print  " fastigmatism3_vpp ",round(data[3],3), data[7], amp,round(ang,2),round(peak,3)
	# PAP 06/16/2018 - there seems to be a mistake here as amplitude angle should be #9 on the list
	# So I corrected it
	data[9] = ang
	return  peak


def ornq_vpp(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi = 0.0):
	"""Determine shift and rotation between image and reference image (refim)
	   no mirror
		quadratic interpolation
		cnx, cny in FORTRAN convention
	"""
	from math import pi, cos, sin, radians
	from alignment import ang_n
	#from utilities import info
	#print "ORNQ"
	peak = -1.0E23

	lkx = int(xrng[0]/step)
	rkx = int(xrng[-1]/step)

	lky = int(yrng[0]/step)
	rky = int(yrng[-1]/step)

	for i in range(-lky, rky+1):
		iy = i*step
		for j in range(-lkx, rkx+1):
			ix = j*step
			cimage = Util.Polar2Dm(image, cnx+ix, cny+iy, numr, mode)
			Util.Frngs(cimage, numr)
			Util.Normalize_ring(cimage, numr, 0)
			retvals = Util.Crosrng_e(crefim, cimage, numr, 0, deltapsi)
			qn = retvals["qn"]
			if qn >= peak:
				sx = -ix
				sy = -iy
				ang = ang_n(retvals["tot"], mode, numr[-1])
				peak = qn
	# mirror is returned as zero for consistency
	mirror = 0
	co =  cos(radians(ang))
	so = -sin(radians(ang))
	sxs = sx*co - sy*so
	sys = sx*so + sy*co
	return  ang, sxs, sys, mirror, peak

def Xdefocusgett_vpp2(qse, roo, nx, xdefc, xampcont, voltage=300.0, Pixel_size=1.0, Cs=2.0, f_start=-1.0, f_stop=-1.0, round_off=1.0, nr1=3, nr2=6, parent=None, DEBug=False):
	"""
		1. Estimate envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get
		   defocus which matches the extracted CTF imprints
	"""
	from utilities  import generate_ctf
	import numpy as np
	from morphology import ctf_2, bracket_def, defocus_baseline_fit, ctflimit, simpw1d, goldsearch_astigmatism

	#print "CTF params:", voltage, Pixel_size, Cs, wgh, f_start, f_stop, round_off, nr1, nr2, parent

	if f_start == 0 : 	    i_start = 0
	else: 			        i_start = int(Pixel_size*nx*f_start+0.5)
	if f_stop <= f_start :
		i_stop  = len(roo)
		adjust_fstop = True
	else:
		i_stop  = min(len(roo), int(Pixel_size*nx*f_stop+0.5))
		adjust_fstop = False

	nroo = len(roo)

	if DEBug:  print("f_start, i_start, f_stop, i_stop:", f_start, i_start, f_stop, i_stop-1)
	#TE  = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr1), 4)
	#baseline = defocus_baseline_fit(roo, i_start, i_stop, int(nr2), 3)

	baseline = defocus_baseline_fit(roo, i_start,nroo, int(nr2), 3)
	subpw = np.array(roo, np.float32) - baseline
	subpw[0] = subpw[1]
	#write_text_file([roo,baseline,subpw],"dbg.txt")
	#print "IN defocusgett  ",np.min(subpw),np.max(subpw)
	for i in range(len(subpw)):  subpw[i] = max(subpw[i],0.0)
	#print "IN defocusgett  ",np.min(subpw),np.max(subpw)
	#envelope = movingaverage(  subpw   , nroo//4, 3)
	envelope = np.array([1.0]*len(subpw), np.float32)
	#write_text_file([roo,baseline,subpw],"dbgt.txt")

	#print "IN defocusgett  ",np.min(subpw),np.max(subpw),np.min(envelope)
	#envelope = np.ones(nroo, np.float32)
	defocus = 0.0
	ampcont = 0.0
	data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
	wn = 512
	from utilities import model_circle, model_blank, amoeba
	from alignment import Numrinit, ringwe
	mask = model_circle(i_stop - 1, wn, wn) * (model_blank(wn, wn, 1, 1.0) - model_circle(i_start, wn, wn))
	from fundamentals import rot_avg_table
	zizi = rot_avg_table(qse)[i_start:i_stop]
	from utilities import write_text_file
	dudi = subpw[i_start:i_stop]
	#print dudi.tolist()
	#print zizi

	cnx = wn // 2 + 1
	cny = cnx
	mode = "H"
	numr = Numrinit(i_start, i_stop, 1, mode)
	wr = ringwe(numr, mode)
	
	crefim = Util.Polar2Dm(qse*mask, cnx, cny, numr, mode)
	print("  CREFIM    ",Util.infomask(qse*mask,None,True),Util.infomask(crefim,None,True))
	Util.Frngs(crefim, numr)
	Util.Applyws(crefim, numr, wr)
	bdef = 0.
	baco = 0.0  #  amplitude contrast
	bamp = 0.0      #  initial astigmatism amplitude
	bang = 0.0      #  initial astigmatism angle
	astdata = [crefim, numr, wn, bdef, Cs, voltage, Pixel_size, baco, bamp, bang, mask]
	data2d = [nx, qse, bdef, Cs, voltage, Pixel_size, 0.0, baco, bamp, bang, mask]

	print(" i_start:i_stop",i_start,i_stop)

	qm = 1.e23
	dm = 1.e23
	dp = 1.0e23
	toto = []
	'''
	#for a in xrange(0,101,10):
	for a in xrange(20,21,10):
		data[7] = float(a)
		print "  fdasfdsfa  ",a
		#for i in xrange(1000,100000,50):
		for i in xrange(5000,5001,5):
			dc = float(i)/10000.0
			qt,ct1 = simpw1dc(dc, data)
			write_text_file(ct1,"testi1.txt")
			write_text_file(data[0],"testd1.txt")
			ju1 = dc # defocus
			ju2 = float(a) # amp contrast
			ju3 = 0.0  # astigma amp
			#dama = amoeba([ju1,ju2,ju3], [0.02, 1.0, 0.02], fupw_vpp, 1.e-4, 1.e-4, 1, astdata)
			data2d[7] = float(a)
			zigi,ct2 = simpw2dc(dc, data2d)
			ct2.write_image("testi2.hdf")
			data2d[1].write_image("testd2.hdf")
			#print  dc,data[7],qt,dama
			toto.append([dc,data[7],qt,zigi])#,dama[-2]])
	'''
	for aa in range(0,20,4):
		a = xampcont + aa
		data[7] = float(a)
		print("  fdasfdsfa  ",a)
		for i in range(0,2000,200):
			dc = xdefc + float(i-1000)/10000.0
			qt = simpw1d(dc, data)
			ju1 = dc # defocus
			ju2 = float(a) # amp contrast
			ju3 = 0.0  # astigma amp
			dama = amoeba([ju1,ju2,ju3], [0.005, 2.0, 0.002], fupw_vpp, 1.e-4, 1.e-4, 200, astdata)
			data2d[7] = float(a)
			zigi = simpw2d(dc, data2d)
			qma = -dama[-2]
			print(" amoeba  %7.2f  %7.2f  %12.6g  %12.6g  %12.6g  %7.2f  %7.2f  %7.2f "%(dc,data[7],qma,zigi,qt,dama[0][0],dama[0][1],dama[0][2]), dama)
			toto.append([dc,data[7],qt,zigi,qma])
			if(qma<dp):
				dp = qma
				dpefi = dama[0][0]
				dpmpcont = dama[0][1]
			if(zigi<dm):
				dm = zigi
				ddefi = dc
				dampcont = data[7]
			if(qt<qm):
				qm = qt
				defi = dc
				ampcont = data[7]
	if DEBug:
		from utilities import write_text_row
		write_text_row(toto,"toto1.txt")
		print(" repi3  ",dp,dpefi,dpmpcont)
		print(" resi2  ",qm,defi,ampcont)
		print(" resi1  ",dm,ddefi,dampcont)
		
		#print " >>>>>>>>>  ",defi,simpw1d(defi, data)#,generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont])
		#def1 = defi
	#exit()
	from morphology import ctf2_rimg, ctf_rimg, square_root
	ctf2 = ctf_rimg(nx, generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont]), sign=0)
	cq = ctf_1d(nx, generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont]), doabs = True)[20:150]
	qse.write_image("qse.hdf")
	ctf2.write_image("c1.hdf")
	ctf22 = ctf_rimg(nx, generate_ctf([ddefi, Cs, voltage, Pixel_size, 0.0, dampcont]), sign=0)
	ci = ctf_1d(nx, generate_ctf([ddefi, Cs, voltage, Pixel_size, 0.0, dampcont]), doabs = True)[20:150]
	dq = ctf_1d(nx, generate_ctf([dpefi, Cs, voltage, Pixel_size, 0.0, dpmpcont]), doabs = True)[20:150]
	write_text_file([dudi.tolist(),zizi,cq,ci,dq],"pwds.txt")
	ctf22.write_image("c2.hdf")
	'''
	from utilities import write_text_file
	foki = subpw.tolist()
	write_text_file([foki,ctf2[:len(foki)]],"toto1.txt")
	'''
	return defi, ampcont, subpw, ctf2, baseline, envelope, i_start, i_stop


def Xdefocusgett_vpp22(qse, roo, nx, voltage=300.0, Pixel_size=1.0, Cs=2.0, f_start=-1.0, f_stop=-1.0, round_off=1.0, nr1=3, nr2=6, parent=None, DEBug=False):
	"""
		1. Estimate envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get
		   defocus which matches the extracted CTF imprints
	"""
	from utilities  import generate_ctf
	import numpy as np
	from morphology import ctf_2, bracket_def, defocus_baseline_fit, ctflimit, simpw1d, goldsearch_astigmatism

	#print "CTF params:", voltage, Pixel_size, Cs, wgh, f_start, f_stop, round_off, nr1, nr2, parent

	if f_start == 0 : 	    i_start = 0
	else: 			        i_start = int(Pixel_size*nx*f_start+0.5)
	if f_stop <= f_start :
		i_stop  = len(roo)
		adjust_fstop = True
	else:
		i_stop  = min(len(roo), int(Pixel_size*nx*f_stop+0.5))
		adjust_fstop = False

	nroo = len(roo)

	if DEBug:  print("f_start, i_start, f_stop, i_stop:", f_start, i_start, f_stop, i_stop-1)
	#TE  = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr1), 4)
	#baseline = defocus_baseline_fit(roo, i_start, i_stop, int(nr2), 3)

	baseline = defocus_baseline_fit(roo, i_start,nroo, int(nr2), 3)
	subpw = np.array(roo, np.float32) - baseline
	subpw[0] = subpw[1]
	#write_text_file([roo,baseline,subpw],"dbg.txt")
	#print "IN defocusgett  ",np.min(subpw),np.max(subpw)
	for i in range(len(subpw)):  subpw[i] = max(subpw[i],0.0)
	#print "IN defocusgett  ",np.min(subpw),np.max(subpw)
	#envelope = movingaverage(  subpw   , nroo//4, 3)
	envelope = np.array([1.0]*len(subpw), np.float32)
	#write_text_file([roo,baseline,subpw],"dbgt.txt")

	#print "IN defocusgett  ",np.min(subpw),np.max(subpw),np.min(envelope)
	#envelope = np.ones(nroo, np.float32)
	defocus = 0.0
	ampcont = 0.0
	data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
	wn = 512
	from utilities import model_circle, model_blank, amoeba
	from alignment import Numrinit, ringwe
	mask = model_circle(i_stop - 1, wn, wn) * (model_blank(wn, wn, 1, 1.0) - model_circle(i_start, wn, wn))
	from fundamentals import rot_avg_table
	zizi = rot_avg_table(qse)[i_start:i_stop]
	from utilities import write_text_file
	dudi = subpw[i_start:i_stop]
	#print dudi.tolist()
	#print zizi

	cnx = wn // 2 + 1
	cny = cnx
	mode = "H"
	numr = Numrinit(i_start, i_stop, 1, mode)
	wr = ringwe(numr, mode)
	
	crefim = Util.Polar2Dm(qse*mask, cnx, cny, numr, mode)
	Util.Frngs(crefim, numr)
	Util.Applyws(crefim, numr, wr)
	bdef = 0.
	baco = 0.0  #  amplitude contrast
	bamp = 0.0      #  initial astigmatism amplitude
	bang = 0.0      #  initial astigmatism angle
	astdata = [crefim, numr, wn, bdef, Cs, voltage, Pixel_size, baco, bamp, bang, mask]
	data2d = [nx, qse, bdef, Cs, voltage, Pixel_size, 0.0, baco, bamp, bang, mask]

	print(" i_start:i_stop",i_start,i_stop)

	qm = 1.e23
	dm = 1.e23
	dp = 1.0e23
	toto = []
	'''
	#for a in xrange(0,101,10):
	for a in xrange(20,21,10):
		data[7] = float(a)
		print "  fdasfdsfa  ",a
		#for i in xrange(1000,100000,50):
		for i in xrange(5000,5001,5):
			dc = float(i)/10000.0
			qt,ct1 = simpw1dc(dc, data)
			write_text_file(ct1,"testi1.txt")
			write_text_file(data[0],"testd1.txt")
			ju1 = dc # defocus
			ju2 = float(a) # amp contrast
			ju3 = 0.0  # astigma amp
			#dama = amoeba([ju1,ju2,ju3], [0.02, 1.0, 0.02], fupw_vpp, 1.e-4, 1.e-4, 1, astdata)
			data2d[7] = float(a)
			zigi,ct2 = simpw2dc(dc, data2d)
			ct2.write_image("testi2.hdf")
			data2d[1].write_image("testd2.hdf")
			#print  dc,data[7],qt,dama
			toto.append([dc,data[7],qt,zigi])#,dama[-2]])
	'''
	for a in range(5,96,10):
		data[7] = float(a)
		print("  fdasfdsfa  ",a)
		for i in range(1000,100000,5000):
			dc = float(i)/10000.0
			qt = simpw1d(dc, data)
			ju1 = dc # defocus
			ju2 = float(a) # amp contrast
			ju3 = 0.0  # astigma amp
			dama = amoeba([ju1,ju2,ju3], [0.002, 0.001, 0.002], fupw_vpp, 1.e-4, 1.e-4, 1, astdata)
			data2d[7] = float(a)
			zigi = simpw2d(dc, data2d)
			qma = dama[-2]/42.
			print(" amoeba  %7.2f  %7.2f  %12.6g  %12.6g  %12.6g  %7.2f  %7.2f  %7.2f "%(dc,data[7],qma,zigi,qt,dama[0][0],dama[0][1],dama[0][2]), dama)
			toto.append([dc,data[7],qt,zigi,qma])
			if(qma<dp):
				dp = qma
				dpefi = dc
				dpmpcont = data[7]
			if(zigi<dm):
				dm = zigi
				ddefi = dc
				dampcont = data[7]
			if(qt<qm):
				qm = qt
				defi = dc
				ampcont = data[7]
	if DEBug:
		from utilities import write_text_row
		write_text_row(toto,"toto1.txt")
		print(" repi3  ",dp,dpefi,dpmpcont)
		print(" resi2  ",qm,defi,ampcont)
		print(" resi1  ",dm,ddefi,dampcont)
		
		#print " >>>>>>>>>  ",defi,simpw1d(defi, data)#,generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont])
		#def1 = defi
	#exit()
	from morphology import ctf2_rimg, ctf_rimg, square_root
	ctf2 = ctf_rimg(nx, generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont]), sign=0)
	cq = ctf_1d(nx, generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont]), doabs = True)[20:150]
	qse.write_image("qse.hdf")
	ctf2.write_image("c1.hdf")
	ctf22 = ctf_rimg(nx, generate_ctf([ddefi, Cs, voltage, Pixel_size, 0.0, dampcont]), sign=0)
	ci = ctf_1d(nx, generate_ctf([ddefi, Cs, voltage, Pixel_size, 0.0, dampcont]), doabs = True)[20:150]
	dq = ctf_1d(nx, generate_ctf([dpefi, Cs, voltage, Pixel_size, 0.0, dpmpcont]), doabs = True)[20:150]
	write_text_file([dudi.tolist(),zizi,cq,ci,dq],"pwds.txt")
	ctf22.write_image("c2.hdf")
	'''
	from utilities import write_text_file
	foki = subpw.tolist()
	write_text_file([foki,ctf2[:len(foki)]],"toto1.txt")
	'''
	return defi, ampcont, subpw, ctf2, baseline, envelope, i_start, i_stop

