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
		from utilities import model_blank
		nx = f.get_xsize()
		ny = f.get_ysize()
		nz = f.get_zsize()
		if(nz == 1):	mask = model_blank(3,3,bckg = 1.0)
		elif(nz >1):  mask = model_blank(3,3,3,bckg = 1.0)
		else:  ERROR("Command does not work for 1D images","dilation",1)

	if morphtype=="BINARY":
		return filt_dilation_(f, mask, morph_type.BINARY)
	elif morphtype=="GRAYLEVEL":
		return filt_dilation_(f, mask, morph_type.GRAYLEVEL)
	else: print "Unknown dilation type."

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
		if(nz == 1):	mask = model_blank(3,3,bckg = 1.0)
		elif(nz >1):  mask = model_blank(3,3,3,bckg = 1.0)
		else:  ERROR("Command does not work for 1D images","dilation",1)

	if morphtype=="BINARY":
		return filt_erosion_(f, mask, morph_type.BINARY)
	elif morphtype=="GRAYLEVEL":
		return filt_erosion_(f, mask, morph_type.GRAYLEVEL)
	else: print "Unknown erosion type."

def invert(im):
	"""
	 Invert contrast of an image (while keeping the average)
	"""
	p = Util.infomask(im, None, True)
	return ((-1.0*im) + 2*p[0])

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

## CTF related functions
def ctf_1d(nx, ctf, sign = 1):
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
	for i in xrange(length):
		ctf_1[i] = Util.tf(dz, i*scl, voltage, cs, ampcont, bfactor, sign)
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
	dict = ctf.to_dict()
	dz = dict["defocus"]
	cs = dict["cs"]
	voltage = dict["voltage"]
	pixel_size = dict["apix"]
	b_factor = dict["bfactor"]
	ampcont = dict["ampcont"]

	ctf_2  = []
	scl    = 1.0/pixel_size/nx
	length = int(1.41*float(nx/2)) + 1
	ctf_2 = [0.0]*length
	for i in xrange(length):
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

###----D-----------------------		
def defocus_env_baseline_fit(roo, i_start, i_stop, nrank, iswi):
	TMP_roo = []
	curve   = [0]*len(roo)
	for i in xrange(i_start,i_stop,1):	TMP_roo.append(roo[i])
	nc     = -1								 
	TMP = imf_params_cl1(TMP_roo, nrank, iswi)
	for i in xrange(i_start, i_stop, 1):
		nc      += 1
		curve[i] = TMP[1][nc]
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
		for i in xrange(len(TMP_roo)):	roo.append(TMP_roo[i][1])
	else:
		TMP_roo=read_text_row(fnam_roo,"s",";")
 		for i in xrange(len(TMP_roo)):	roo.append(TMP_roo[i][2])
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
	for i in xrange(len(roo)):
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
	res     = []
	Res_roo = []
	Res_TE  = []	
	if f_start == 0 : 	i_start = 0
	else: 			i_start = int(Pixel_size*2.*len(roo)*f_start)
	if f_stop <= f_start : 	i_stop  = len(roo)
		
	else: 			
		i_stop  = int(Pixel_size*2.*len(roo)*f_stop)
		if i_stop > len(roo): i_stop  = len(roo)

	#print "f_start, i_start, f_stop, i_stop:", f_start, i_start, f_stop, i_stop
	TE  = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr1), 4)
	Pn1 = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr2), 3)
	Res_roo = []
	Res_TE  = []
	for i in xrange(len(roo)):
		Res_roo.append(roo[i] - Pn1[i])
		Res_TE.append( TE[i]  - Pn1[i])

	defocus = defocus_guess(Res_roo, Res_TE, voltage, Cs, Pixel_size, wgh, i_start, i_stop, 2, round_off)

	nx  = int(len(Res_roo)*2)
	#ctf = ctf_1d(nx, generate_ctf([defocus, Cs, voltage, Pixel_size, 0.0, wgh]))
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
		for i in xrange(len(TMP_roo)): # remove first record
			roo.append(TMP_roo[i][1])
	else:
		skip = ";"
		TMP_roo = read_text_row(fnam_roo, "s", skip)
 		for i in xrange(len(TMP_roo)): # remove first record
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
	for i in xrange(len(roo)):
		Res_roo.append( roo[i] - Pn1[i] )
		Res_TE.append(  TE[i]  - Pn1[i] )
	#
	defocus=defocus_guess(Res_roo, Res_TE, volt, Cs, Pixel_size, wgh, i_start, i_stop, 2, round_off)
	nx  = int(len(roo)*2)
	ctf = ctf_2(nx, generate_ctf([defocus,Cs,voltage,Pixel_size, 0.0, wgh]))
	for i in xrange(len(Res_TE)):
		ctf[i]=ctf[i]*Res_TE[i]
	dis = defocus_L2_euc(ctf, Res_roo, i_start, i_stop)
	return [defocus, dis]

def defocus_L2_euc(v1,v2, ist,istp):
	from math import sqrt
	dis    = 0.0
	pw_sum = 0.0
	if ist == istp :	ERROR("No pw2 curve is included  ", "defocus_L2_euc", 0)
	else:			tfeq=istp-ist
	for i in xrange(ist,istp,1):
		dis+=    (v1[i]-v2[2])**2
		pw_sum+= (v1[i])**2
	if pw_sum <= 0:		ERROR("negative or zero power ", "defocus_L2_euc", 1)
	if dis    <= 0:		ERROR("bad fitting, change options settings and try again  ", "defocus_L2_euc", 0)
	else:	
		res = sqrt(dis)/sqrt(pw_sum)/tfeq	
		return res

def defocus_guess(Res_roo, Res_TE, volt, Cs, Pixel_size, wgh, istart=0, istop=-1, defocus_estimation_method=2, round_off=1, dz_low=1000., dz_high=200000., nloop=100, ampcont=0.0):
	"""
		Use specified frequencies area (istart-istop)to estimate defocus
		1.  The searching range is limited to dz_low (.1um) ~ dz_high (20 um).
		    The user can modify this limitation accordingly
		2.  changing nloop can speed up the estimation
		3.  mode=1 compare squared error 
		    mode=2; compare cross correlation coefficients 
	"""
	
	from math import sqrt
	from utilities import generate_ctf

	if istop <= istart : 			istop=len(Res_roo)
	step = (dz_high-dz_low)/nloop
	if step > 10000.   : 			step     =  10000.     # angstrom 
	if defocus_estimation_method == 1 : 	diff_min =  1.e38
	else: 					diff_min = -1.e38
	xval_e = 0
	for ifreq in xrange(len(Res_TE)):
		xval_e += Res_TE[ifreq]**2
	if (xval_e == 0):
		defocus = 0 
		return defocus
	if round_off >= 1: 			cut_off  =  1.
	else: 					cut_off  =  round_off # do extreme fitting	
	while (step >= cut_off):
		for i_dz in xrange(nloop):
			dz     = dz_low + step*i_dz
			diff   = 0
			length = len(Res_roo)
			nx     = int(length*2)
			ctf    = ctf_2(nx, generate_ctf([dz,Cs,volt,Pixel_size,ampcont,wgh]))
			if defocus_estimation_method == 1:
	        		for ifreq in xrange(istart, istop, 1):
	        			diff += (ctf[ifreq]*Res_TE[ifreq] - Res_roo[ifreq])**2
	        	       	if diff < diff_min or dz == dz_low:
	        		       defocus  = dz
	        		       diff_min = diff
			else:
				diff  = 0.0
				sum_a = 0.0
				sum_b = 0.0
	        		for ifreq in xrange(istart, istop, 1):
	        		       xval   =  ctf[ifreq]**2*Res_TE[ifreq]
	        		       diff  +=  Res_roo[ifreq]**2*xval
	        		       sum_a +=  Res_roo[ifreq]**4
	        		       sum_b +=  xval**2
	        	       	diff/=(sqrt(sum_a*sum_b)*( istop - istart + 1 ))
	        	       	if diff > diff_min or dz == dz_low:
	        		       defocus  = dz
	        		       diff_min = diff
	        		       xscore   = diff_min

		dz_low = defocus-step*2
		if( dz_low < 0 ): 	dz_low=0.0
		dz_high = defocus + step*2
		step/=10.
	defocus=int( defocus/round_off )*round_off
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
	if os.path.exists(indir) == False					     : 	ERROR("roodir doesn't exist", "defocus_get_fast",1)
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
					if(print_screen == "p" or print_screen == "P"): print "defocus",defocus,"Euclidean distance", dis, "starting feq", istart, "stop freq", istop,"P R E", nr1,"P R B", nr2
					if(defocus < dz_max): 				break
			if(defocus >= dz_max): 					ERROR("defocus_get_fast fails at estimating defocus", fnam, action = 0)
			print "micrograph", flist[i], '%5d'%(defocus) 	# screen output, give the user a general impression about estimated defoci
			if(writetodoc[0] == "w" or writetodoc[0] != "l"):	out.write("%d\t%5d\t%s\n" % (ncount,defocus,flist[i]))
			if(writetodoc[0] == "l"):				res.append(defocus)
			if type(micdir) is types.StringType : 
				ctf_param = [defocus, Pixel_size, volt, Cs, wgh, 0, 1]
				mic_name  = os.path.join(micdir,"micrograph"+ fnam_root+ ".hdf")
				if os.path.exists(mic_name) :
					e = get_image (mic_name)
					U______set_arb_params(e, ctf_param, ctf_dicts)  # THIS IS INCORRECT< PLEASE CHANGE
					e.write_image(mic_name,0, EMUtil.ImageType.IMAGE_HDF, True)
					print "ctf parameters is written back into headers of ", mic_name
				#else :  print  mic_name, " Not found"
	if(len(res) == 0 and  writetodoc == "l" ):				ERROR("No input file is found, check the input directory of file prefix", indir, 1)
	else:
		if(writetodoc[0] == "a"):
			out.close()
			return res
		if(writetodoc[0] == "l"): 	return res
		if(writetodoc[0] == "w"): 	out.close()

def defocus_get_fast_MPI(indir, writetodoc="w", Pixel_size=1, volt=120, Cs=2, wgh=.1, round_off=100, dz_max0=50000, f_l0=30, f_h0=5, nr_1=5, nr_2=5, prefix_of_micrograph="roo", docf="a",skip="#",print_screen="no"):
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
		if(filename[0:len(prefix_of_micrograph)] == prefix_of_micrograph):
			mic_name_list.append(micname)
			nima += 1
	if nima < 1: 	ERROR("No micrograph is found, check either directory or prefix of micrographs is correctly given","pw2sp",1)
	
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
	for i in xrange(image_start,image_end):
		filename=mic_name_list[i] 
		print '%-15s%-30s'%("micrographs # ",filename)
		(f_nam, filextension) = os.path.splitext(filename)
		fnam_roo     = "particle_"+f_nam[len(prefix_of_micrograph)+len(indir)+2:]+filextension	
#	for i, v in enumerate(flist):
		ncount   += 1
		defocus = defocus_get(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
		if(defocus > dz_max):
			while(nr1 <= 7 or nr2 <= 7):
				nr1 += 1
				nr2 += 1
				defocus = defocus_get(fnam_roo, volt, Pixel_size, Cs, wgh,istart, istop, docf,skip, round_off, nr1, nr2)
				if(print_screen[0] == "p" or print_screen[0] == "P" ): print "defocus",defocus,"Euclidean distance",dis,"starting feq",istart,"stop freq",istop,"P R E", nr1,"P R B", nr2
				if(defocus<dz_max): break
		if(defocus > dz_max):
			while(nr1 >= 2 and nr2 >= 2):
				nr1 -= 1
				nr2 -= 1
				defocus = defocus_get(fnam_roo, volt,Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
				if(print_sreen[0] == "p" or print_screen=="P"): print "defocus",defocus,"Euclidean distance",dis,"starting feq",istart,"stop freq",istop,"P R E", nr1,"P R B", nr2
				if(defocus < dz_max): break
		if(defocus > dz_max):
			while(istart > istop):
		 		nr1    =  5
				nr2    =  5
				istart -=.5
				defocus = defocus_get(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf,skip, round_off, nr1, nr2)
				if(print_screen[0] == "p" or print_screen == "P"): print "defocus",defocus,"Euclidean distance",dis,"starting feq",istart,"stop freq",istop,"P R E", nr1,"P R B", nr2
				if(defocus < dz_max): break
		if(defocus > dz_max):
			while(istart > istop):
				nr1     = 5										    	
				nr2     = 5
				istop  += 0.5
				defocus = defocus_get(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
				if(print_screen == "p" or print_screen == "P"): print "defocus",defocus,"Euclidean distance", dis, "starting feq", istart, "stop freq", istop,"P R E", nr1,"P R B", nr2
				if(defocus < dz_max): 				break
		if(defocus >= dz_max): 					ERROR("defocus_get_fast fails at estimating defocus", fnam, action = 0)
		print "micrograph", flist[i], '%10.3g'(defocus) 	# screen output, give the user a general impression about estimated defoci
		if(writetodoc[0] == "w" or writetodoc[0] != "l"):	out.write("%d\t%f\t%s\n" % (ncount,defocus,flist[i]))
		if(writetodoc[0] == "l"):				res.append(defocus)
	if(len(res) == 0 and  writetodoc == "l" ):				ERROR("No input file is found, check the input directory of file prefix", indir, 1)
	else:
		if writetodoc[0] == "a":
			out.close()
			return res
	if(writetodoc[0] == "l"): 	return res
	if(writetodoc[0] == "w"): 	out.close()

def defocus_get_slow(indir, writetodoc="w", Pixel_size=1, volt=120, Cs=2, wgh=.1, round_off=100, dz_max0=50000, f_l0=30, f_h0=5, prefix="roo", docf="s", skip=";",micdir="micrograph", print_screen="p"):
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
			for nr1 in xrange(2,7,1):
				for nr2 in xrange(2,7,1):
					[defocus, dis]     = defocus_get_Eudis(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
					if(print_screen[0]=="p"): print "defocus",defocus,"Euclidean distance",dis,"starting feq",istart,"stop freq",istop,"P R E", nr1,"P R B", nr2
					if(Mdis > dis):
						defo = defocus
						Mdis = dis
			if(defo > dz_max):
				istart-= 1.
				for nr1 in xrange(3,5,1):
					for nr2 in xrange(2,4,1):
						[defocus, dis] = defocus_get_Eudis(fnam_roo, volt, Pixel_size, Cs, wgh, istart, istop, docf, skip, round_off, nr1, nr2)
						if(Mdis>dis):
							defo = defocus
							Mdis = dis
			if(defo >= dz_max): 	ERROR("defo_get_s fails at estimating defocus from ", fnam, 0)
			else:				print "micrograph", flist[i], defo # screen output, give the user a general impression about estimated defoci		
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
	mic_sq     = ccf(tmp,m_pad)/n_pixelt 	 # calculate the average of squared mic	       
	tmp        = mic_sq-mic_avg_sq*n_pixelt   #  
	mic_var    = tmp.get_pow(.5)          # Calculate the local variance of the image 
	cc_map     = ccf(e,t_pad)
	cc_map    /= (mic_var*n_pixelt) # Normalize the cross correlation map 
	return cc_map

##-----------------------------img formation parameters related functions---------------------------------
def imf_params_cl1(pw,n=2,iswi=3,Pixel_size=1):
	"""
		Extract image formation parameters using contrained simplex method
		The output is a list of list, which contains the following four elements:
		1. frequencies in 1/Angstrom
		2. fitted curve, either background noise or envelope function
		3. original power spectrum to be fitted
		4. The parameters
		Attension:
		        iswi= 2 using polynomial n rank to fit no-Gaussian envelope function
			iswi =3 using polynomail n rank to fit background
			n = the polynomial rank +1
			The optimization tend to fail when the polynomial rank is higher than 6 
	"""
	res  = []
	feq  = []
	cur  = []
	parm = []
	t    = Util.pw_extract(pw, n, iswi, Pixel_size)
	for i in xrange(len(pw)):
		j = i*2
		k = i*2+1
		cur.append(t[j])
		feq.append(t[k])
	npam = len(t)-2*len(pw)
	for i in xrange(npam):
		k    = 2*len(pw)+i
		parm.append(t[k])
	res.append(feq )
	res.append(cur )
	res.append(pw  )
	res.append(parm)	
	return res

def imf_get_1dpw_list(fstr):
	pw   = []
	data = read_spider_doc(fstr)
	for i in xrange(len(data)):
		pw.append(data[i][0])
	return pw
	
def imf_B_factor_get(res_N,x,ctf_params):
	from scipy.optimize import fmin
	nx    = len(res_N)*2
	ctf   = ctf_1d(nx, ctf_params)
	p     = [1,1]	
	xopt  = fmin(residuals_B1, p, (res_N,x))
	p     = xopt
	xopt1 = fmin(residuals_B2, p, (res_N,ctf[1][0:nx-1], x))
	print  xopt
	return xopt

def imf_residuals_B1(p,y,x):
	"""
		Give the initial guess of B-factor
	"""
	from numpy import exp
	C,B = p
	err = 0.0
	for i in xrange(len(y)):
		err+= abs(y[i] - C*exp(-B*x[i]*x[i]))  # should be 4*B
	return err

def imf_residuals_B2(p,y,ctf,x):
	"""
		fit B-factor in case of considering CTF effect
	""" 
	from numpy import exp
	C,B = p
	err = 0.0
	for i in xrange(len(y)):
		err+= abs(y[i] - ctf[i]*C*exp(-B*x[i]*x[i]))  # should be 4*B
	return err

def imf_params_get(fstrN, fstrP, ctf_params, pu, nrank, q, lowf=0.01):
	"""
		Extract image formation paramters using optimization method
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
	for i in xrange(len(t_N[0])):
		res_N.append(t_N[2][i] - t_N[1][i])
		res_P.append(t_P[2][i] - t_N[1][i])
	params.append(t_N[0]) # freq
	params.append(t_N[1]) # baseline
#	params.append(t_N[1])
	parm1  = imf_B_factor_get(res_N,t_N[0],ctf_params)
	params.append(parm1[1])
	n_lowf = lowf*ctf_params[0]*len(res_P)*2
	n_lowf = int(n_lowf)
	for i in xrange(len(res_P)):
		if(i <= n_lowf):
			w.append(0.)
		else:
			w.append(1.)
	parm2 = imf_fit_pu(res_P,t_N[0],ctf_params,pu,parm1[0],parm1[1],q,w)
	params.append(parm2[1])
	params.append(parm2[0])
	for i in xrange(len(res_N)):
		res_N[i]*= q
	params.append(res_N)
	return params

def imf_fit_pu(res_P, x, ctf_params, pu, C, B, q, w):
	from scipy.optimize import fmin
	res   = []
	nx    = len(res_P)*2
	ctf   = ctf_1d(nx, ctf_params)
	for i in xrange(len(pu)):
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
	for i in xrange(len(y)):
		err+= abs(y[i] - C*pu[i])
	return err

def residuals_simplex(args, data):
	err      = 0.0
	for i in xrange(len(data[0])):  err -= (data[0][i] - (args[0] + (args[1]/(data[1][i]/args[2]+1.0)**2)))**2
	return err

def residuals_lsq(p,y,x):
	c1,c2,c3 = p
	err	 = []
	for i in xrange(len(y)):
		err.append(abs(y[i] - c1-c2/(x[i]+c3)**2))
	return err

def residuals_lsq_peak(p,y,x,c):
	from numpy import exp
	d1,d2,d3 = p
	c1,c2,c3 = c
	err	 = []
	for i in xrange(len(y)):
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
	for i in xrange(len(list_1dpw2)):
		j = i*2
		k = i*2+1
		if i <= cut_off:
			res.append(list_1dpw2[i]-background[i])
			freq.append(i/(2*Pixel_size*len(list_1dpw2)))
		else : 
			res.append(0.0)
			freq.append(i/(2*Pixel_size*len(list_1dpw2)))
	return res, freq

def adaptive_mask(vol, nsigma = 1.0, ndilation = 3, kernel_size = 11, gauss_standard_dev =9):
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
	nx = vol.get_xsize()
	ny = vol.get_ysize()
	nz = vol.get_zsize()
	mc = model_circle(nx//2, nx, ny, nz) - model_circle(nx//3, nx, ny, nz)
	s1 = Util.infomask(vol, mc, True)
	mask = Util.get_biggest_cluster(binarize(vol, s1[0]+s1[1]*nsigma))
	for i in xrange(ndilation):   mask = dilation(mask)
	mask = gauss_edge(mask, kernel_size, gauss_standard_dev)
	return mask

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
	for i in xrange(ndilation):   mask = dilation(mask)
	#mask = gauss_edge(mask, kernel_size, gauss_standard_dev)
	return mask


'''

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
		for k in xrange(nz):
			if(fif):
				for j in xrange(ny):
					if(fif):
						for i in xrange(nx):
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
			for ig in xrange(1,grp):
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
			for k in xrange(-1,2,1):
				kq = cr[2]+k
				if(kq>-1 and kq<nz):
					for j in xrange(-1,2,1):
						jq = cr[1]+j
						if(jq>-1 and jq<ny):
							for i in xrange(-1,2,1):
								iq = cr[0]+i
								if(iq>-1 and iq<nx):
									if(lg[iq,jq,kq]):
										lg[iq,jq,kq]=0
										l.append([iq,jq,kq])

def adaptive_mask(vol, mass=2000, Pixel_size=3.6):
	from utilities  import gauss_edge, model_blank
	from morphology import binarize, threshold
	from filter     import filt_gaussl, filt_dilation
	nx = vol.get_xsize()
	a = filt_gaussl(vol, 0.15, True)
	TH = a.find_3d_threshold(mass, Pixel_size)
	a = binarize(a,TH)
	d = a.delete_disconnected_regions(0,0,0)

	d = filt_dilation(d, model_blank(3,3,3,1.0), "BINARY")
	#d = filt_dilation(d, model_blank(3,3,3,1.0), "BINARY")
	d = gauss_edge(d)
	return d
	#Util.mul_img(vol, d)
	#return threshold(vol, 0.0)

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
'''


'''
Start code for helical consistency
'''

# global variables
helical_ref_weight = 200.0
helical_nonref_weight = 0.5
helical_dphi = -166.5
helical_dz = 27.6 # Angstroms
helical_pixel_size = 1.84
helical_sgnfil=-1000
helical_ref = -1
helical_ysgn = -1000
helical_filtheta=90
helical_ptclcoords=[]
helical_thetapsi1 = []
helical_THR_CONS_PHI=1.5
helical_w0=[]
helical_THR_CONS_Y=0.5
helical_w0y=[]

def S6(w):
        tot = 0
        np = len(helical_thetapsi1)
        for ii in xrange(0,np):
                jj = helical_ref 
                wij = get_iphi(helical_thetapsi1[ii], helical_thetapsi1[jj], w[jj], helical_dz, helical_dphi, helical_pixel_size,helical_ptclcoords,helical_filtheta,helical_sgnfil)
                conscost = (abs( (w[ii]%360.0) - (wij%360.0)))%360.
                conscost = min(conscost, 360.0-conscost)
                        
                dfcost = 0
                if ii == helical_ref:
                        dfi = (abs((w[ii]%360.0) - (helical_w0[ii]%360.0)))%360.0
                        wid = min(dfi, 360.-dfi)
                        dfcost = wid*helical_ref_weight 
                        
                tot += (dfcost + conscost)
        
        return tot        

def S6y(w):
        tot = 0
        np = len(helical_thetapsi1)
        for ii in xrange(0,np):
                jj = helical_ref
                ycons = get_iy(helical_thetapsi1[ii],helical_thetapsi1[jj],w[jj], helical_dz, helical_pixel_size,helical_ptclcoords,helical_filtheta,helical_ysgn)
               
                minycost = -1
                yii_min = -1
                for yii in ycons:
                        ycost = abs(yii - helical_w0y[ii])
                        if (minycost < 0) or (ycost < minycost):
                                minycost = ycost
                                yii_min = yii
                
                ycost = abs(yii_min - w[ii])
             
                dfty = 0
                if ii == helical_ref:
                        dfyi = abs(w[ii] - helical_w0y[ii])
                        dfty = dfyi*helical_ref_weight
                        
                tot += (ycost + dfty)
        
        return tot   

def get_delta_p_y(aD_p, aip, aiy, apref):
        ymin = aiy[0]
        for y in aiy:
                if abs(ymin - aip) > abs(y - aip):
                        ymin = y
        bestD = abs(ymin - apref)
        delta_p = bestD - aD_p                
        return delta_p, bestD
       
def get_delta_p(aD_p, aiphi, apref):
        D = abs(aiphi - apref) # ideal distnace between ref and seg calculated wrt D_p
        D_1 = (aiphi-apref)%360.
        D_2 = min(D_1, 360.-D_1)
        D_3 = max(D_1, 360.-D_1)
        bestD = D
        delta_p = bestD - aD_p
        if abs(delta_p) > abs(D_1 - aD_p):
                bestD = D_1
                delta_p = D_1 - aD_p
        if abs(delta_p) > abs(D_2 - aD_p):
                bestD = D_2
                delta_p = D_2 - aD_p
        if abs(delta_p) > abs(D_3 - aD_p):
                bestD = D_3
                delta_p = D_3 - aD_p
        return delta_p, bestD          
              
def get_neighborhoods(refseg, ps, newparams,dz, dphi, pixel_size, ptclcoords,filtheta,sgnfil,delta_phi,nbrphi):
        pref = newparams[refseg][0]
        
        aref ={}
        aref_g0 ={}
        aref_e0 ={}
        
        D={}
        delta_p={}
        
        for iseg in ps:
                if iseg == refseg:
                        continue
                
                ip = newparams[iseg][0]
                
                D_p = abs(ip - pref) # absolute value between ref and iseg
                iphi = get_iphi(iseg,refseg,pref, dz, dphi, pixel_size,ptclcoords,filtheta,sgnfil)
                idelta_p, iD = get_delta_p(D_p, iphi, pref)
                D[iseg] = iD
                delta_p[iseg]= idelta_p
                
                if abs(delta_p[iseg]) >= delta_phi:
                        print "1 enforced level of consistency is too not strict enough cmpared to desired level of consistency"
                        
                        sys.exit()
                
                dsgn = 1.0      
                if (ip - pref) == D_p:
                        dsgn=-1.0
                a1 = dsgn*delta_p[iseg] - delta_phi
                a2 = dsgn*delta_p[iseg] + delta_phi
                
                if a1 >=0 or a2 <= 0:
                        print "something wrong with a1 and a2, phi"
                        sys.exit()
                        
                if abs(2*delta_phi) <= D_p:
                        aref[iseg] = [a1/2.0,a2/2.0]
                        continue
                
                if D_p > 0:
                        aref_g0[iseg] = [a1/2.0, a2/2.0]
                
                        if abs(aref_g0[iseg][0]) > D_p:
                                aref_g0[iseg][0] = -D_p/2.0
                        if abs(aref_g0[iseg][1]) > D_p:
                                aref_g0[iseg][1] = D_p/2.0
                        
                delta_pp = D[iseg] + D_p
                
                if abs(delta_pp) <= delta_phi:
                       
                        dsgn2 = 1.0      
                        if (pref- ip) == D_p:
                                dsgn2=-1.0
                        a1e0 = dsgn2*delta_pp - delta_phi
                        a2e0 = dsgn2*delta_pp + delta_phi
                        if a1e0 > a1:
                                a1 = a1e0
                        if a2e0 < a2:
                                a2 = a2e0
                        
                        aref_e0[iseg] = [a1/2.0,a2/2.0]         
               
        aseg ={}
        aseg_g0 ={}
        aseg_e0 ={}
        allkeys = aref.keys() + aref_e0.keys() + aref_g0.keys()
        for iseg in ps:
                if iseg == refseg:
                        continue
                if not(iseg in allkeys):
                        print "iseg is not covered in ref keys"
                        print D_p
                        sys.exit()
        for iseg in ps:
                if iseg == refseg:
                        continue
                ip = newparams[iseg][0]
               
                D_p = abs(ip - pref) # absolute value between ref and iseg
                
                if abs(delta_p[iseg]) >= delta_phi:
                        print "abs(delta_p) >= delta_phi"
                        sys.exit()
                        
                dsgn = 1.0      
                if (pref - ip) == D_p:
                        dsgn=-1.0
                        
                if iseg in aref.keys():
                
                        a1ref = aref[iseg][0]
                        a2ref = aref[iseg][1]
                        a1 = dsgn*delta_p[iseg] - delta_phi + a2ref
                        a2 = dsgn*delta_p[iseg] + delta_phi + a1ref
                
                        if a1 > 0 or a2 < 0:
                                print "1 something wrong with a1, a2"
                                sys.exit()
                        
                        aseg[iseg]=[a1,a2]
                        continue
                        
                if iseg in aref_e0.keys():
                
                        a1ref = aref_e0[iseg][0]
                        a2ref = aref_e0[iseg][1]
                        
                        a1 = dsgn*delta_p[iseg] - delta_phi + a2ref
                        a2 = dsgn*delta_p[iseg] + delta_phi + a1ref
                        
                        if a1 > 0 or a2 < 0:
                                print "1 something wrong with a1, a2"
                                sys.exit()
                        
                        aseg_e0[iseg] = [a1, a2]    
                        delta_pp = D[iseg] + D_p
                        
                        dsgn2 = 1.0      
                        if (ip-pref) == D_p:
                                dsgn2=-1.0 
                        
                        a1e0 = dsgn2*delta_pp - delta_phi + a2ref
                        a2e0 = dsgn2*delta_pp + delta_phi + a1ref
                
                        if a1e0 > 0 or a2e0 < 0:
                                print "2 something wrong with a1, a2"
                                sys.exit()
                        
                        if a1e0 > aseg_e0[iseg][0] :
                                aseg_e0[iseg][0] = a1e0
                        if a2e0 < aseg_e0[iseg][1]:
                                aseg_e0[iseg][1] = a2e0
                                
                if iseg in aref_g0.keys():
                        a1ref = aref_g0[iseg][0]
                        a2ref = aref_g0[iseg][1]
                        a1 = dsgn*delta_p[iseg] - delta_phi + a2ref
                        a2 = dsgn*delta_p[iseg] + delta_phi + a1ref
                
                        if a1 > 0 or a2 < 0:
                                print "1 something wrong with a1, a2"
                                sys.exit()
                                
                        # adjust a1, a2 so that |max(|a1|,|a2|)| + |max(|a1ref|,|a2ref|)| < D_p
                        small_enough = False
                        aa1ref =abs(a1ref)
                        aa2ref =abs(a2ref)
                        aa1 = abs(a1)
                        aa2 = abs(a2)
                        
                        while not(small_enough):
                                
                                if max(aa1, aa2) + max(aa1ref, aa2ref) < D_p:
                                        small_enough = True
                                else:
                                        if aa1 > aa2:
                                                aa1 = 0.95*aa1
                                        else:
                                                aa2 = 0.95*aa2
                        aseg_g0[iseg]=[-aa1, aa2]
                # set aseg[iseg] to the larger of the two intervals in aseg_e0 and aseg_g0      
                e0len = -1
                g0len = -1
                e0reflen = -1
                g0reflen = -1
                if iseg in aref_e0.keys():
                        e0len = aseg_e0[iseg][1] - aseg_e0[iseg][0]
                        e0reflen = aref_e0[iseg][1] - aref_e0[iseg][0]
                        
                if iseg in aref_g0.keys():
                        g0len = aseg_g0[iseg][1] - aseg_g0[iseg][0]
                        g0reflen = aref_g0[iseg][1] - aref_g0[iseg][0]
                        
                if e0len < 0 and g0len < 0 and e0reflen < 0 and g0reflen < 0:
                        print "something wrong"
                        sys.exit()
                
                if ((e0len + e0reflen) > (g0len + g0reflen)):
                        aseg[iseg] = [aseg_e0[iseg][0],aseg_e0[iseg][1]]
                        aref[iseg]=[aref_e0[iseg][0],aref_e0[iseg][1]]
                else:
                        aseg[iseg] = [aseg_g0[iseg][0],aseg_g0[iseg][1]]
                        aref[iseg] = [aref_g0[iseg][0],aref_g0[iseg][1]]
              
        # set ref's neighborhood to intersection of aref[iseg] for all iseg
        a1ref = 1000
        a2ref = -1000
        for iseg in ps:
                if iseg != refseg:
                        
                        nbrphi[iseg]=[aseg[iseg][0], aseg[iseg][1]]
               
                        if aref[iseg][0] > a1ref or a1ref > 0:
                                a1ref = aref[iseg][0]
                
                        if aref[iseg][1] < a2ref or a2ref < 0:
                                a2ref = aref[iseg][1]  
        nbrphi[refseg] = [a1ref, a2ref]
        
                
def get_neighborhoods_y(refseg, ps, newparams, dz, pixel_size,ptclcoords,filtheta,ysgn,delta_y,nbry):
        
        dpp_half = dz/pixel_size/2.0
        pref = newparams[refseg][4]
        
        aref ={}
        aref_g0 ={}
        aref_e0 ={}
        
        D={}
        delta_p={}
        
        for iseg in ps:
                if iseg == refseg:
                        continue
                
                ip = newparams[iseg][4]
                
                D_p = abs(ip - pref) # absolute value between ref and iseg
                iy = get_iy(iseg,refseg,pref, dz, pixel_size,ptclcoords,filtheta,ysgn)
                idelta_p, iD = get_delta_p_y(D_p, ip, iy, pref)
                D[iseg] = iD
                delta_p[iseg]=idelta_p 
                
                if abs(delta_p[iseg]) >= delta_y:
                        print "1 enforced level of consistency is too not strict enough cmpared to desired level of consistency"
                        sys.exit()
                          
                dsgn = 1.0      
                if (ip - pref) == D_p:
                        dsgn=-1.0
                a1 = dsgn*delta_p[iseg] - delta_y
                a2 = dsgn*delta_p[iseg] + delta_y
                
                if a1 >=0 or a2 <= 0:
                        print "something wrong with a1 and a2, phi"
                        sys.exit()
                        
                if abs(2*delta_y) <= D_p:
                        aref[iseg] = [a1/2.0,a2/2.0]
                        if pref+aref[iseg][0] < -1*dpp_half:
                                aref[iseg][0] = (-1*dpp_half) - pref
                        if pref + aref[iseg][1] > dpp_half:
                                aref[iseg][1] = dpp_half - pref        
                        continue
                
                if D_p > 0:
                        aref_g0[iseg] = [a1/2.0, a2/2.0]
                
                        if abs(aref_g0[iseg][0]) > D_p:
                                aref_g0[iseg][0] = -D_p/2.0
                        if abs(aref_g0[iseg][1]) > D_p:
                                aref_g0[iseg][1] = D_p/2.0
                
                        if pref+aref_g0[iseg][0] < -1*dpp_half:
                                aref_g0[iseg][0] = (-1*dpp_half) - pref
                        if pref + aref_g0[iseg][1] > dpp_half:
                                aref_g0[iseg][1] = dpp_half - pref     
                                
                delta_pp = D[iseg] + D_p
                
                if abs(delta_pp) <= delta_y:
                       
                        dsgn2 = 1.0      
                        if (pref- ip) == D_p:
                                dsgn2=-1.0
                        a1e0 = dsgn2*delta_pp - delta_y
                        a2e0 = dsgn2*delta_pp + delta_y
                        if a1e0 > a1:
                                a1 = a1e0
                        if a2e0 < a2:
                                a2 = a2e0
                        
                        aref_e0[iseg] = [a1/2.0,a2/2.0]         
                        
                        if pref+aref_e0[iseg][0] < -1*dpp_half:
                                aref_e0[iseg][0] = (-1*dpp_half) - pref
                        if pref + aref_e0[iseg][1] > dpp_half:
                                aref_e0[iseg][1] = dpp_half - pref    
                                
        aseg ={}
        aseg_g0 ={}
        aseg_e0 ={}
        allkeys = aref.keys() + aref_e0.keys() + aref_g0.keys()
        for iseg in ps:
                if iseg == refseg:
                        continue
                if not(iseg in allkeys):
                        print "iseg is not covered in ref keys"
                        print D_p
                        sys.exit()
                
        for iseg in ps:
                if iseg == refseg:
                        continue
                ip = newparams[iseg][4]
               
                D_p = abs(ip - pref) # absolute value between ref and iseg
               
                if abs(delta_p[iseg]) >= delta_y:
                        print "abs(delta_p) >= delta_y"
                        sys.exit()
                        
                dsgn = 1.0      
                if (pref - ip) == D_p:
                        dsgn=-1.0
                
                if iseg in aref.keys():
                
                        a1ref = aref[iseg][0]
                        a2ref = aref[iseg][1]
                        a1 = dsgn*delta_p[iseg] - delta_y + a2ref
                        a2 = dsgn*delta_p[iseg] + delta_y + a1ref
                
                        if a1 > 0 or a2 < 0:
                                print "1 something wrong with a1, a2"
                                sys.exit()
                        
                        aseg[iseg]=[a1,a2]
                        
                        if ip+aseg[iseg][0] < -1*dpp_half:
                                aseg[iseg][0] = (-1*dpp_half) - ip
                        if ip + aseg[iseg][1] > dpp_half:
                                aseg[iseg][1] = dpp_half - ip
                        
                        continue
            
                if iseg in aref_e0.keys():
                
                        a1ref = aref_e0[iseg][0]
                        a2ref = aref_e0[iseg][1]
                        
                        a1 = dsgn*delta_p[iseg] - delta_y + a2ref
                        a2 = dsgn*delta_p[iseg] + delta_y + a1ref
                        
                        if a1 > 0 or a2 < 0:
                                print "1 something wrong with a1, a2"
                                sys.exit()
                        
                        aseg_e0[iseg] = [a1, a2]    
                        delta_pp = D[iseg] + D_p
                        
                        dsgn2 = 1.0      
                        if (ip-pref) == D_p:
                                dsgn2=-1.0 
                        
                        a1e0 = dsgn2*delta_pp - delta_y + a2ref
                        a2e0 = dsgn2*delta_pp + delta_y + a1ref
                
                        if a1e0 > 0 or a2e0 < 0:
                                print "2 something wrong with a1, a2"
                                sys.exit()
                        
                        if a1e0 > aseg_e0[iseg][0] :
                                aseg_e0[iseg][0] = a1e0
                        if a2e0 < aseg_e0[iseg][1]:
                                aseg_e0[iseg][1] = a2e0
        
                if iseg in aref_g0.keys():
                        a1ref = aref_g0[iseg][0]
                        a2ref = aref_g0[iseg][1]
                        a1 = dsgn*delta_p[iseg] - delta_y + a2ref
                        a2 = dsgn*delta_p[iseg] + delta_y + a1ref
                
                        if a1 > 0 or a2 < 0:
                                print "1 something wrong with a1, a2"
                                sys.exit()
                                
                        # adjust a1, a2 so that |max(|a1|,|a2|)| + |max(|a1ref|,|a2ref|)| < D_p
                        small_enough = False
                        aa1ref =abs(a1ref)
                        aa2ref =abs(a2ref)
                        aa1 = abs(a1)
                        aa2 = abs(a2)
                        while not(small_enough):
                                
                                if max(aa1, aa2) + max(aa1ref, aa2ref) < D_p:
                                        small_enough = True
                                else:
                                        if aa1 > aa2:
                                                aa1 = 0.95*aa1
                                        else:
                                                aa2 = 0.95*aa2
                                                
                        aseg_g0[iseg]=[-aa1, aa2]
                        
                # set aseg[iseg] to the larger of the two intervals in aseg_e0 and aseg_g0      
                e0len = -1
                g0len = -1
                e0reflen = -1
                g0reflen = -1
                if iseg in aref_e0.keys():
                        e0len = aseg_e0[iseg][1] - aseg_e0[iseg][0]
                        e0reflen = aref_e0[iseg][1] - aref_e0[iseg][0]
                        
                if iseg in aref_g0.keys():
                        g0len = aseg_g0[iseg][1] - aseg_g0[iseg][0]
                        g0reflen = aref_g0[iseg][1] - aref_g0[iseg][0]
                        
                if e0len < 0 and g0len < 0 and e0reflen < 0 and g0reflen < 0:
                        print "something wrong"
                        sys.exit()
                
                if ((e0len + e0reflen) > (g0len + g0reflen)):
                        aseg[iseg] = [aseg_e0[iseg][0],aseg_e0[iseg][1]]
                        aref[iseg]=[aref_e0[iseg][0],aref_e0[iseg][1]]
                else:
                        aseg[iseg] = [aseg_g0[iseg][0],aseg_g0[iseg][1]]
                        aref[iseg] = [aref_g0[iseg][0],aref_g0[iseg][1]]
        
                if ip+aseg[iseg][0] < -1*dpp_half:
                        aseg[iseg][0] = (-1*dpp_half) - ip
                if ip + aseg[iseg][1] > dpp_half:
                        aseg[iseg][1] = dpp_half - ip
         
         # set ref's neighborhood to intersection of aref[iseg] for all iseg
        a1ref = 1000
        a2ref = -1000
        for iseg in ps:
                if iseg != refseg:
                        
                        nbry[iseg]=[aseg[iseg][0], aseg[iseg][1]]
               
                        if aref[iseg][0] > a1ref or a1ref > 0:
                                a1ref = aref[iseg][0]
                
                        if aref[iseg][1] < a2ref or a2ref < 0:
                                a2ref = aref[iseg][1]  
        nbry[refseg] = [a1ref, a2ref]
                
def get_dist(ix, iy, jx, jy, theta):
        from math import pi, sqrt, sin
        qv = pi/180
        d = sqrt((ix - jx)**2 + (iy - jy)**2)
	d = (d/abs(sin(theta*qv)))
	return d

def get_iphi(iseg, refseg, phiref, dz, dphi, pixel_size,ptclcoords,filtheta,sgnfil):
        ix = ptclcoords[iseg][0]
        iy = ptclcoords[iseg][1]
        jx = ptclcoords[refseg][0]
        jy = ptclcoords[refseg][1]
        d = pixel_size * get_dist(ix, iy, jx, jy, filtheta)
        if iseg > refseg:
                iphi = (((phiref%360.0) + sgnfil*round(d/dz)*dphi)%360.0)
        else:
                iphi = (((phiref%360.0) - sgnfil*round(d/dz)*dphi)%360.0)   
        return iphi
        
def get_iy(iseg, refseg, yref, dz, pixel_size,ptclcoords,filtheta,ysgn):
        
        ix = ptclcoords[iseg][0]
        iy = ptclcoords[iseg][1]
        refx = ptclcoords[refseg][0]
        refy = ptclcoords[refseg][1]
        d = pixel_size * get_dist(ix, iy, refx, refy, filtheta)
      
        dbar = d%dz
        dbar1 = dz-dbar
      
        ycons=[] # all possible helical consistent y-shifts
             
        if iseg < refseg:
                yii_1 = ysgn*dbar
                yii_2 = -1*ysgn*dbar1
        else:
                yii_1 = -1*ysgn*dbar
                yii_2 = ysgn*dbar1
        
        yii_1 += (yref*pixel_size)
        yii_2 += (yref*pixel_size)
        
        if yii_1 >= dz:
                yii_1 = yii_1%dz
        if yii_1 <= -dz:
                yii_1 = -1.0*((abs(yii_1))%dz)
        
        if yii_2 >= dz:
                yii_2 = yii_2%dz
        if yii_2 <= -dz:
                yii_2 = -1.0*((abs(yii_2))%dz)
       
        if abs(yii_1) <= 0.5*dz:
                ycons.append(yii_1/pixel_size)
        if abs(yii_2) <= 0.5*dz:
                ycons.append(yii_2/pixel_size)
        
        return ycons
        
def find_params_phi(w0, ref, ps, dz, dphi, pixel_size,ptclcoords,filtheta,sgnfil):
        tot = 0
        np = len(ps)
        for ii in xrange(0,np):
                jj = ref 
                wij = get_iphi(ps[ii], ps[jj], w0[jj], dz, dphi, pixel_size,ptclcoords,filtheta,sgnfil)
                conscost = (abs( (w0[ii]%360.0) - (wij%360.0)))%360.
                conscost = min(conscost, 360.0-conscost)
                tot += (conscost)
        
        return tot        

def find_params_y(w0y, ref, ps, dz,pixel_size,ptclcoords,filtheta,ysgn):
        tot = 0
        np = len(ps)
        rng = dz/pixel_size
        for ii in xrange(0,np):
                jj = ref
                ycons = get_iy(ps[ii],ps[jj],w0y[jj], dz, pixel_size,ptclcoords,filtheta,ysgn)
               
                minycost = -1
                yii_min = -1
                for yii in ycons:
                        ycost = abs(yii - w0y[ii])
                        if (minycost < 0) or (ycost < minycost):
                                minycost = ycost
                                yii_min = yii
                
                ycost = (abs(yii_min - w0y[ii]))*(360./rng)
                tot += (ycost)
        
        return tot   
        

def num_cons_segs(w,ps,ref, dz, dphi, pixel_size ,ptclcoords, sgnfil, THR=3.5, STRICT=False,filtheta=90):
        ncons=0
        nt = 0
        np = len(ps)
        for ii in xrange(0,np):
                jj = ref
                
                wij = get_iphi(ps[ii], ps[jj], w[jj], dz, dphi, pixel_size,ptclcoords,filtheta,sgnfil)
                conscost = (abs( (w[ii]%360.0) - (wij%360.0)))%360.
                conscost = min(conscost, 360.0-conscost)
                    
                nt += 1
                
                if not(STRICT):
                        if conscost <= THR:
                                ncons += 1
                else:
                        if conscost < THR:
                                ncons += 1
                #else:
                #        print "not consistent: ", ps[ii]
                #        print wij, w[ii]
        return nt, ncons

def num_cons_segs_y(w, ps, ref, dz, pixel_size, ptclcoords, ysgn, THR=1.5, STRICT=False,filtheta=90):
        ncons=0
        nt = 0
        np = len(ps)
        for ii in xrange(0,np):
                jj = ref
                
                ycons = get_iy(ps[ii],ps[jj],w[jj], dz, pixel_size,ptclcoords,filtheta,ysgn)
                minycost = -1
                yii_min = -1
                for yii in ycons:
                        ycost = abs(yii - w[ii])
                        if (minycost < 0) or (ycost < minycost):
                                minycost = ycost
                                yii_min = yii
                
                ycost = abs(yii_min - w[ii])
                
                nt += 1
                
                if not(STRICT):
                        if abs(ycost) <= THR:
                                ncons += 1
                else:
                        if abs(ycost) < THR:
                                ncons += 1
                #else: 
                #        print "y not consistent: ", ps[ii]
                #        print yii_min, w[ii],w0y[ii]
                #        sys.exit()
        return nt, ncons   


# ps2 are the IDs of segments whose phi angle need to be predicted
def predict_phi(ps2, refseg, refphi, filtheta,dz, dphi, pixel_size,ptclcoords,sgnfil):
        np = len(ps2)
        prphi=[0 for iiiii in xrange(np)]
        for ii in xrange(np):
                iphi = get_iphi(ps2[ii], refseg, refphi, dz, dphi, pixel_size,ptclcoords,filtheta,sgnfil)
                prphi[ii] = iphi
        return prphi
        
def predict_y(ps2, refseg, refy, filtheta,dz, pixel_size,ptclcoords,ysgn):
        np = len(ps2)
        pry=[0 for iiiii in xrange(np)]
        for ii in xrange(np):
                iyall = get_iy(ps2[ii], refseg, refy, dz, pixel_size,ptclcoords,filtheta,ysgn)
                ymin = iyall[0]
                for yshift in iyall:
                        if abs(yshift) < abs(ymin):
                                ymin = yshift
                pry[ii] = ymin
                if abs(ymin - refy) > 0.16:
                        print "ymin not equal to refy: ", ps2[ii], refseg, iyall, refy
                        sys.exit()
        return pry
'''
 Level of consistency enforeced on parameters
        THR_CONS_PHI=1.5
        THR_CONS_Y=1.0

 Desired level of consistency AFTER refinement
        delta_phi = 3.5
        delta_y = 1.5
        
 dz is in Angstroms
'''
def helical_consistency(parmfile, miclistfile, ptclcoordsfile, THR_CONS_PHI=1.5,THR_CONS_Y=1.0,delta_phi = 3.5, delta_y = 1.5, dphi = -166.5,dz = 27.6, pixel_size = 1.84, fileNewIDs='newIDs_morph.txt', fileNewParams='newparams2_morph.txt', fileNbrPhi='nbrphi_morph.txt', fileNbrY='nbry_morph.txt'):        
        
        from utilities import read_text_row, write_text_file, write_text_row
        from scipy.optimize import minimize
        
        global helical_dphi
        global helical_dz
        global helical_pixel_size
        global helical_sgnfil
        global helical_ref
        global helical_filtheta
        global helical_ptclcoords
        global helical_thetapsi1
        global helical_ysgn
        global helical_THR_CONS_PHI
        global helical_w0
        global helical_THR_CONS_Y
        global helical_w0y
        
        helical_dphi = dphi
        helical_dz = dz
        helical_pixel_size =pixel_size
        
        miclist = read_text_row(miclistfile)
        ptclcoords=read_text_row(ptclcoordsfile)
        helical_ptclcoords = ptclcoords    
        params=read_text_row(parmfile)  
        
        helical_THR_CONS_PHI=THR_CONS_PHI
        helical_THR_CONS_Y = THR_CONS_Y
        MAXIT = 20

        dpp_half = dz/pixel_size/2.0

        # For C1, so psi angles for a particular filament should be ~same if consistent
        # theta = 90
        nima = len(params)
        newparams=[[] for iiiii in xrange(nima)]
        nbrphi=[[] for iiiii in xrange(nima)]
        nbry=[[] for iiiii in xrange(nima)]
        N = len(miclist)
        newIDs = []

        for i in xrange(N):
                print "filament ", i
                mic = miclist[i][6:]
                mic = map(int, mic)
                a90=[]
                a270=[]
                for iseg in mic:
                        ipsi = params[iseg][2]
                        if ipsi - 90 > 90:
                                a90.append(iseg)
                        else:
                                a270.append(iseg)
                if len(a90) > len(a270):
                        thetapsi1 = a90
                else:
                        thetapsi1 = a270
       
                if (len(a90) == len(a270)) or len(mic)==1:
                        # these cannot be predicted. Don't include them in final stack
                        continue
                newIDs.extend(mic)     
                helical_thetapsi1= thetapsi1
                filtheta=params[thetapsi1[0]][1]
                helical_filtheta = filtheta
                w0 = []
                for iseg in thetapsi1:
                        w0.append(params[iseg][0])
                helical_w0=w0
                w0y=[]
                for iseg in thetapsi1:
                        w0y.append(params[iseg][4])  
                helical_w0y = w0y
                best_cost = -1
                best_sgnfil = 100
                best_ref=-1
                best_ysgn = -100
        
                for ref in xrange(len(thetapsi1)):
                        if len(thetapsi1) > 1:
                                if abs(params[thetapsi1[ref]][4] * pixel_size) > 0.5*dz:
                                        continue
                        for sgnfil in [-1,1]:
                                for ysgn in [-1,1]:
                                        startphicost = find_params_phi(w0, ref, thetapsi1, dz, dphi,pixel_size,ptclcoords,filtheta,sgnfil)
                                        startycost = find_params_y(w0y, ref, thetapsi1,dz,pixel_size,ptclcoords,filtheta,ysgn)
                                        #startphicost=0.0*startphicost
                                        #startycost=0.0*startycost
                                        cost = startphicost + startycost
                               
                                        if (best_cost < 0) or (cost < best_cost):
                                                best_cost = cost
                                                best_sgnfil = sgnfil
                                                best_ref = ref
                                                best_ysgn = ysgn
                sgnfil = best_sgnfil 
                helical_sgnfil = sgnfil  
                if best_ref < 0:
                        print "no references round, all segments had shifts > dpp_half"
                        sys.exit()
                ref = best_ref
                helical_ref =ref
                ysgn=best_ysgn
                helical_ysgn = ysgn
        
                res = w0
        
                # iterate until we get to the desired level of consistency
                #print "optimize filament ", i
                h1, h2 = num_cons_segs(res,thetapsi1, ref, dz, dphi, pixel_size,ptclcoords, sgnfil,THR=THR_CONS_PHI, STRICT=True,filtheta=filtheta)
                counter = 0
                while h2 < h1:
                        res = minimize(S6, w0, method='Powell')   
                        h1, h2 = num_cons_segs(res,thetapsi1, ref, dz, dphi, pixel_size,ptclcoords,sgnfil, THR=THR_CONS_PHI, STRICT=True,filtheta=filtheta)
                        counter = counter + 1
                        if counter > MAXIT:
                                print "cannot reach desired level of consistency for phi"
                                sys.exit()
        
                resy = w0y  
                h1y, h2y = num_cons_segs_y(resy,thetapsi1, ref,dz,pixel_size,ptclcoords, ysgn,THR=THR_CONS_Y, STRICT=True,filtheta=filtheta)
        
                counter = 0
                while h2y < h1y:
                        resy = minimize(S6y, resy, method='Powell')  
                        h1y, h2y = num_cons_segs_y(resy,thetapsi1, ref, dz, pixel_size, ptclcoords,ysgn,THR=THR_CONS_Y, STRICT=True,filtheta=filtheta)
                        counter = counter + 1
                        if counter > MAXIT:
                                print "cannot reach desired level of consistency for y"
                                sys.exit()
                    
                for ii in xrange(len(thetapsi1)):
                        newparams[thetapsi1[ii]] = [(res[ii])%360.,params[thetapsi1[ii]][1],params[thetapsi1[ii]][2],params[thetapsi1[ii]][3], resy[ii]]
        
                # predict what consistent params should be for everything not in thetapsi1
        
                thetapsi2 = []
                for iseg in mic:
                        if not(iseg in thetapsi1):
                                thetapsi2.append(iseg)
        
                # theta nad psi of segments in thetapsi2 will be set to that of the reference's phi and theta
        
                reftheta = params[thetapsi1[ref]][1]
                refpsi = params[thetapsi1[ref]][2]
        
                # phi angles and y-shifts of segments in thetapsi2 will be predicted using the reference's adjusted params (newparams)
        
                predphi = predict_phi(thetapsi2, thetapsi1[ref], newparams[thetapsi1[ref]][0], filtheta, dz, dphi, pixel_size,ptclcoords,sgnfil)
                predy = predict_y(thetapsi2, thetapsi1[ref], newparams[thetapsi1[ref]][4], filtheta,dz, pixel_size,ptclcoords,ysgn)
        
                for ii in xrange(len(thetapsi2)):
                        newparams[thetapsi2[ii]] = [(predphi[ii])%360., reftheta, refpsi,params[thetapsi2[ii]][3], predy[ii]]
        
                # First determine allowed neighborhood for the reference
                refseg = thetapsi1[ref]
                thetapsi3 = thetapsi1 + thetapsi2
                get_neighborhoods(refseg,thetapsi3, newparams, dz, dphi, pixel_size,ptclcoords,filtheta,sgnfil,delta_phi,nbrphi)
                get_neighborhoods_y(refseg,thetapsi3, newparams, dz, pixel_size,ptclcoords,filtheta,ysgn,delta_y,nbry)   
        
      
        newIDs.sort()
        nnima=len(newIDs)
        newparams2=[[] for jj in xrange(nnima)]
        testparams=[[] for jj in xrange(nnima)]
        nbrphi2 = [[] for jj in xrange(nnima)]
        nbry2 = [[] for jj in xrange(nnima)]
        for ii in xrange(nnima):
                iseg = newIDs[ii]
                newparams2[ii] = [newparams[iseg][0],newparams[iseg][1],newparams[iseg][2],newparams[iseg][3],newparams[iseg][4]]
                testparams[ii] = [newparams[iseg][0],newparams[iseg][1],newparams[iseg][2],newparams[iseg][3],newparams[iseg][4]]

                nbrphi2[ii] = [nbrphi[iseg][0], nbrphi[iseg][1]]
                nbry2[ii] = [nbry[iseg][0], nbry[iseg][1]]

        write_text_file(newIDs,fileNewIDs)  
        write_text_row(newparams2,fileNewParams)  # consistency enforced on small_stack_parameters (from elmar)
        write_text_row(nbrphi2,fileNbrPhi)  
        write_text_row(nbry2,fileNbrY)  

def helical_consistency_predict(parmfile, miclistfile, ptclcoordsfile, dphi = -166.5,dz = 27.6, pixel_size = 1.84, fileNewIDs='newIDs.txt', fileNewParams='pred_params.txt'):        
        
        from utilities import read_text_row, write_text_file, write_text_row
        
        miclist = read_text_row(miclistfile)
        ptclcoords=read_text_row(ptclcoordsfile)
        helical_ptclcoords = ptclcoords    
        params=read_text_row(parmfile)  

        dpp_half = dz/pixel_size/2.0

        # For C1, so psi angles for a particular filament should be ~same if consistent
        # theta = 90
        nima = len(params)
        newparams=[[] for iiiii in xrange(nima)]
        N = len(miclist)
        newIDs = []
        
        for i in xrange(N):
                mic = miclist[i][6:]
                mic = map(int, mic)
                a90=[]
                a270=[]
                for iseg in mic:
                        ipsi = params[iseg][2]
                        if ipsi - 90 > 90:
                                a90.append(iseg)
                        else:
                                a270.append(iseg)
                if len(a90) > len(a270):
                        thetapsi1 = a90
                else:
                        thetapsi1 = a270
       
                if (len(a90) == len(a270)) or len(mic)==1:
                        # these cannot be predicted. Don't include them in final stack
                        continue
                newIDs.extend(mic)     
                filtheta=params[thetapsi1[0]][1]
                w0 = []
                for iseg in thetapsi1:
                        w0.append(params[iseg][0])
                
                w0y=[]
                for iseg in thetapsi1:
                        w0y.append(params[iseg][4])  
                
                best_cost = -1
                best_sgnfil = 100
                best_ref=-1
                best_ysgn = -100
        
                for ref in xrange(len(thetapsi1)):
                        if len(thetapsi1) > 1:
                                if abs(params[thetapsi1[ref]][4] * pixel_size) > 0.5*dz:
                                        continue
                        for sgnfil in [-1,1]:
                                for ysgn in [-1,1]:
                                        startphicost = find_params_phi(w0, ref, thetapsi1, dz, dphi,pixel_size,ptclcoords,filtheta,sgnfil)
                                        startycost = find_params_y(w0y, ref, thetapsi1,dz,pixel_size,ptclcoords,filtheta,ysgn)
                                        #startphicost=0.0*startphicost
                                        #startycost=0.0*startycost
                                        cost = startphicost + startycost
                               
                                        if (best_cost < 0) or (cost < best_cost):
                                                best_cost = cost
                                                best_sgnfil = sgnfil
                                                best_ref = ref
                                                best_ysgn = ysgn
                sgnfil = best_sgnfil  
                if best_ref < 0:
                        print "no references round, all segments had shifts > dpp_half"
                        sys.exit()
                ref = best_ref
                ysgn=best_ysgn
        
                # predict what consistent params should be for everything using ref
                reftheta = params[thetapsi1[ref]][1]
                refpsi = params[thetapsi1[ref]][2]
        
                # phi angles and y-shifts of segments in thetapsi2 will be predicted using the reference's adjusted params (newparams)
        
                predphi = predict_phi(mic, thetapsi1[ref], params[thetapsi1[ref]][0], filtheta, dz, dphi, pixel_size,ptclcoords,sgnfil)
                predy = predict_y(mic, thetapsi1[ref], params[thetapsi1[ref]][4], filtheta,dz, pixel_size,ptclcoords,ysgn)
        
                for ii in xrange(len(mic)):
                        newparams[mic[ii]] = [(predphi[ii])%360., reftheta, refpsi,params[mic[ii]][3], predy[ii]]
        
        newIDs.sort()
        nnima=len(newIDs)
        newparams2=[[] for jj in xrange(nnima)]
        for ii in xrange(nnima):
                iseg = newIDs[ii]
                newparams2[ii] = [newparams[iseg][0],newparams[iseg][1],newparams[iseg][2],newparams[iseg][3],newparams[iseg][4]]

        write_text_file(newIDs,fileNewIDs)  
        write_text_row(newparams2,fileNewParams)  # consistency enforced on small_stack_parameters (from elmar)

def helical_consistency_test(parmfile, miclistfile, ptclcoordsfile,THR_CONS_PHI=1.5,THR_CONS_Y=1.0, dphi = -166.5,dz = 27.6, pixel_size = 1.84):        
        
        from utilities import read_text_row, write_text_file, write_text_row
        
        miclist = read_text_row(miclistfile)
        ptclcoords=read_text_row(ptclcoordsfile)
        params=read_text_row(parmfile)  

        dpp_half = dz/pixel_size/2.0

        # For C1, so psi angles for a particular filament should be ~same if consistent
        # theta = 90
        nima = len(params)
        newparams=[[] for iiiii in xrange(nima)]
        N = len(miclist)
        newIDs = []
        
        totphitested=0.0
        totphicons=0.0
        totytested = 0.0
        totycons=0.0
        notest=0.0
        
        for i in xrange(N):
                mic = miclist[i][6:]
                mic = map(int, mic)
                a90=[]
                a270=[]
                for iseg in mic:
                        ipsi = params[iseg][2]
                        if ipsi - 90 > 90:
                                a90.append(iseg)
                        else:
                                a270.append(iseg)
                if len(a90) > len(a270):
                        thetapsi1 = a90
                else:
                        thetapsi1 = a270
       
                if (len(a90) == len(a270)) or len(mic)==1:
                        # these cannot be predicted. Don't include them in final stack
                        notest += len(mic)
                        continue
                        
                newIDs.extend(mic)     
                notest += (len(mic)-len(thetapsi1))
                filtheta=params[thetapsi1[0]][1]
                w0 = []
                for iseg in thetapsi1:
                        w0.append(params[iseg][0])
                w0y=[]
                for iseg in thetapsi1:
                        w0y.append(params[iseg][4])  
                
                best_cost = -1
                best_sgnfil = 100
                best_ref=-1
                best_ysgn = -100
        
                for ref in xrange(len(thetapsi1)):
                        if len(thetapsi1) > 1:
                                if abs(params[thetapsi1[ref]][4] * pixel_size) > 0.5*dz:
                                        continue
                        for sgnfil in [-1,1]:
                                for ysgn in [-1,1]:
                                        startphicost = find_params_phi(w0, ref, thetapsi1, dz, dphi,pixel_size,ptclcoords,filtheta,sgnfil)
                                        startycost = find_params_y(w0y, ref, thetapsi1,dz,pixel_size,ptclcoords,filtheta,ysgn)
                                        #startphicost=0.0*startphicost
                                        #startycost=0.0*startycost
                                        cost = startphicost + startycost
                               
                                        if (best_cost < 0) or (cost < best_cost):
                                                best_cost = cost
                                                best_sgnfil = sgnfil
                                                best_ref = ref
                                                best_ysgn = ysgn
                sgnfil = best_sgnfil  
                if best_ref < 0:
                        print "no references round, all segments had shifts > dpp_half"
                        sys.exit()
                ref = best_ref
                ysgn=best_ysgn
        
               
                h1, h2 = num_cons_segs(w0,thetapsi1,ref, dz, dphi, pixel_size,ptclcoords,sgnfil,THR=THR_CONS_PHI,filtheta=filtheta)
                h1y, h2y = num_cons_segs_y(w0y, thetapsi1,ref,dz,pixel_size,ptclcoords,ysgn,THR=THR_CONS_Y,filtheta=filtheta)
                totphicons = totphicons + h2
                totphitested = totphitested + h1
                totycons = totycons + h2y
                totytested = totytested + h1y
        
        print "Consistency for segments which had psi agreement: "
        print "         phi: ", totphitested, totphicons,totphicons/totphitested
        print "         y: ", totytested, totycons, totycons/totytested
        print "Inconsistent because lack of psi agreement: ", notest
'''
End code for helical consistency
'''