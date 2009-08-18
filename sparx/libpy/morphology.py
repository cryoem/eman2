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
#
"""
Morphology functions
"""
from EMAN2_cppwrap import *
from global_def import *

def binarize(img, minval = 0.0):
	return img.process( "threshold.binary", {"value": minval} )

def collapse(img, minval = -1.0, maxval = 1.0):
	# for values between minval and maxval set to one, to zero outside
	return img.process( "threshold.binaryrange", {"low": minval, "high":maxval} )
			
def dilation(f, mask = None, morphtype="BINARY"):
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
	return img.process( "math.exp", {"low": 1.0/a, "high":b})

def power(img, x = 3.0):
	return img.process( "math.pow", {"pow": x})

def alog10(img):
	return img.process( "math.log")

def square_root(img):
	[a,b,c,d] = Util.infomask(img, None, False)
	if(c<0.0):  ERROR("Cannot calculate square root of negative pixels","square_root",1)
	return img.process( "math.sqrt" )

def square(img):
	return img.process( "math.squared")

def threshold(img, minval = 0.0):
	return img.process( "threshold.belowtozero", {"minval": minval} )

def threshold_to_zero(img, minval = 0.0):
	return img.process( "threshold.belowtozero_cut", {"minval": minval } )	

def threshold_to_minval(img, minval = 0.0):
	return img.process( "threshold.belowtominval", {"minval": minval } )	

def threshold_outside(img, minval, maxval):
	return img.process( "threshold.clampminmax", {"minval": minval, "maxval": maxval } )	

def threshold_maxval(img, maxval = 0.0):
	st = Util.infomask(img, None, True)	
	return img.process( "threshold.clampminmax", {"minval": st[2], "maxval": maxval } )	

## CTF related functions
def ctf_1d(nx, ctf, sign = 1):
	"""
		Generate a list of 1D CTF values 
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
	  
	"""
	dict = ctf.to_dict()
	dz = dict["defocus"]
	cs = dict["cs"]
	voltage = dict["voltage"]
	pixel_size = dict["apix"]
	b_factor = dict["bfactor"]
	ampcont = dict["ampcont"]


	if(ny < 1):  ny = nx
	return  Util.ctf_img(nx, ny, nz, dz, pixel_size, voltage, cs, ampcont, b_factor, 0.0, 0.0, sign)

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
	from utilities 	import read_txt_col
	roo     = []
	res     = []
	if(docf == "a"):
		TMP_roo = read_txt_col(fnam_roo, "a", skip)
		for i in xrange(len(TMP_roo)):	roo.append(TMP_roo[i][1])
	else:
		TMP_roo=read_txt_col(fnam_roo,"s",";")
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
			
def defocus_gett(roo, voltage=300.0, Pixel_size=1.0, Cs=2.0, wgh=0.1, f_start=0.0, f_stop=-1.0, docf="a", skip="#", round_off=1.0, nr1=3, nr2=6,parent=None):
	"""
	
		1. Estimating envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get 
		   defocus which matches the extracted CTF imprints 
	"""
	from utilities import generate_ctf
	from math 	import sqrt, atan
	res     = []
	Res_roo = []
	Res_TE  = []	
	if f_start == 0 : 	i_start = 0
	else: 			i_start = int(Pixel_size*2.*len(roo)/f_start)
	if f_stop <= i_start : 	i_stop  = len(roo)
	else: 			i_stop  = int(Pixel_size*2.*len(roo)/f_stop)
	#print f_start, i_start, f_stop, i_stop
	TE  = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr1), 4)
	Pn1 = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr2), 3)
	Res_roo = []
	Res_TE  = []	
	for i in xrange(len(roo)):
		Res_roo.append(roo[i] - Pn1[i])
		Res_TE.append( TE[i]  - Pn1[i])
	#
	defocus = defocus_guess(Res_roo, Res_TE, voltage, Cs, Pixel_size, wgh, i_start, i_stop, 2, round_off)
	nx  = int(len(Res_roo)*2)
	ctf = ctf_1d(nx, generate_ctf([defocus,Cs,voltage,Pixel_size, 0.0, wgh]))
	if (parent is not None):
		parent.ctf_data=[roo, Res_roo, Res_TE]
		parent.i_start = i_start
		parent.i_stop = i_stop
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
	from utilities 	import read_txt_col, generate_ctf
	roo     = []
	res     = []
	if docf == "a":
		TMP_roo = read_txt_col(fnam_roo, "a", skip)
		for i in xrange(len(TMP_roo)): # remove first record
			roo.append(TMP_roo[i][1])
	else:
		skip = ";"
		TMP_roo = read_txt_col(fnam_roo, "s", skip)
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
					set_arb_params(e, ctf_param, ctf_dicts)
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
	if os.path.exists(indir) == False					     : 	ERROR("roodir doesn't exist", "defocus_get_fast",1)
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
	t_pad      = Util.pad(t,    {"background":0}, nx, ny, 1, 0, 0, 0)
	m_pad      = Util.pad(mask, {"background":0}, nx, ny, 1, 0, 0, 0) # create a mask (blank, value=1 )file and pad to size of mic   	 	
	tmp        = ccf(e, m_pad)/n_pixele # calculate the local average
	mic_avg_sq = tmp*tmp    # calculate average square
	tmp        = e*e
	mic_sq     = ccf(tmp,m_pad)/n_pixelt 	 # calculate the average of squared mic	       
	tmp        = mic_sq-mic_avg_sq*n_pixelt   #  
	mic_var    = tmp.get_pow(.5)          # Calculate the local variance of the image 
	cc_map     = ccf(e,t_pad)
	cc_map    /= (mic_var*n_pixelt) # Normalize the cross correlation map 
	del tmp
	del mic_var
	del mic_sq
	del mic_avg_sq
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
		err+= abs(y[i] - C*exp(-B*x[i]*x[i]))
	return err
	
def imf_residuals_B2(p,y,ctf,x):
	"""
		fit B-factor in case of considering CTF effect
	""" 
	from numpy import exp
	C,B = p
	err = 0.0
	for i in xrange(len(y)):
		err+= abs(y[i] - ctf[i]*C*exp(-B*x[i]*x[i]))
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
