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
from EMAN2_cppwrap import *
from global_def import *
	
def ali2d_single_iter(data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, mode, CTF=False, random_method="", T=1.0, ali_params="xform.align2d"):
	"""
		single iteration of 2D alignment using ormq
		if CTF = True, apply CTF to data (not to reference!)
	"""
	from utilities import model_circle, combine_params2, inverse_transform2, get_params2D, set_params2D
	from alignment import Applyws, ormq

	if CTF:
		from filter  import filt_ctf

	# 2D alignment using rotational ccf in polar coords and quadratic interpolation
	cimage = Util.Polar2Dm(tavg, cnx, cny, numr, mode)
	Util.Frngs(cimage, numr)
	Applyws(cimage, numr, wr)
	
	maxrin = numr[-1]
	sx_sum = 0.0
	sy_sum = 0.0
	for im in xrange(len(data)):
		if CTF:
			#Apply CTF to image
			ctf_params = data[im].get_attr("ctf")
			ima = filt_ctf(data[im], ctf_params, True)
		else:
			ima = data[im]
		alpha, sx, sy, mirror, dummy = get_params2D(ima, ali_params)
		alpha, sx, sy, mirror        = combine_params2(alpha, sx, sy, mirror, 0.0, -cs[0], -cs[1], 0)
		alphai, sxi, syi, scalei     = inverse_transform2(alpha, sx, sy, 1.0)

		# align current image to the reference
		if random_method == "SA":
			peaks = ormq_peaks(ima, cimage, xrng, yrng, step, mode, numr, cnx+sxi, cny+syi)
			[angt, sxst, syst, mirrort, peakt, select] = sim_anneal(peaks, T, step, mode, maxrin)
			[alphan, sxn, syn, mn] = combine_params2(0.0, -sxi, -syi, 0, angt, sxst, syst, mirrort)
			data[im].set_attr_dict({"select":select, "peak":peakt})
			set_params2D(data[im], [alphan, sxn, syn, mn, 1.0], ali_params)
		else:
			[angt, sxst, syst, mirrort, peakt] = ormq(ima, cimage, xrng, yrng, step, mode, numr, cnx+sxi, cny+syi)
			# combine parameters and set them to the header, ignore previous angle and mirror
			[alphan, sxn, syn, mn] = combine_params2(0.0, -sxi, -syi, 0, angt, sxst, syst, mirrort)
			set_params2D(data[im], [alphan, sxn, syn, mn, 1.0], ali_params)

		if mn == 0: sx_sum += sxn
		else:       sx_sum -= sxn
		sy_sum += syn

	return sx_sum, sy_sum


def ali2d_random_ccf(data, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, mode, CTF=False, random_method="", T=1.0, ali_params="xform.align2d"):
	"""
		single iteration of 2D alignment using ormq
		if CTF = True, apply CTF to data (not to reference!)
	"""
	from utilities import model_circle, combine_params2, inverse_transform2, get_params2D, set_params2D
	from alignment import Applyws, ormq

	if CTF:
		from filter  import filt_ctf

	# 2D alignment using rotational ccf in polar coords and quadratic interpolation
	cimage = Util.Polar2Dm(tavg, cnx, cny, numr, mode)
	Util.Frngs(cimage, numr)
	Applyws(cimage, numr, wr)
	
	maxrin = numr[-1]
	sx_sum = 0.0
	sy_sum = 0.0
	for im in xrange(len(data)):
		if CTF:
			#Apply CTF to image
			ctf_params = data[im].get_attr("ctf")
			ima = filt_ctf(data[im], ctf_params, True)
		else:
			ima = data[im]
		alpha, sx, sy, mirror, dummy = get_params2D(ima, ali_params)
		alpha, sx, sy, mirror        = combine_params2(alpha, sx, sy, mirror, 0.0, -cs[0], -cs[1], 0)
		alphai, sxi, syi, scalei     = inverse_transform2(alpha, sx, sy, 1.0)

		# align current image to the reference
		if random_method == "SA":
			peak = Util.ali2d_ccf_list(ima, cimage, xrng, yrng, step, mode, numr, cnx+sxi, cny+syi, T)
			[angt, sxst, syst, mirrort, peakt, select] = sim_ccf(peak, T, step, mode, maxrin)
			[alphan, sxn, syn, mn] = combine_params2(0.0, -sxi, -syi, 0, angt, sxst, syst, mirrort)
			data[im].set_attr_dict({"select":select, "peak":peakt})
			set_params2D(data[im], [alphan, sxn, syn, mn, 1.0], ali_params)
		else:
			[angt, sxst, syst, mirrort, peakt] = ormq(ima, cimage, xrng, yrng, step, mode, numr, cnx+sxi, cny+syi)
			# combine parameters and set them to the header, ignore previous angle and mirror
			[alphan, sxn, syn, mn] = combine_params2(0.0, -sxi, -syi, 0, angt, sxst, syst, mirrort)
			set_params2D(data[im], [alphan, sxn, syn, mn, 1.0], ali_params)

		if mn == 0: sx_sum += sxn
		else:       sx_sum -= sxn
		sy_sum += syn

	return sx_sum, sy_sum

def ang_n(tot, mode, maxrin):
	"""
	  Calculate angle based on the position of the peak
	"""
	from math import fmod
	if (mode == 'f' or mode == 'F'): return fmod(((tot-1.0) / maxrin+1.0)*360.0, 360.0)
	else:                            return fmod(((tot-1.0) / maxrin+1.0)*180.0, 180.0)

def Applyws(circ, numr, wr):
	"""
	  Apply weights to FTs of rings
	"""
	nring = len(numr)/3
	maxrin = numr[len(numr)-1]
	for i in xrange(nring):
		numr3i = numr[2+i*3]
		numr2i = numr[1+i*3]-1
		w = wr[i]
		circ[numr2i] *= w
		if (numr3i == maxrin): circ[numr2i+1] *= w
		else: circ[numr2i+1] *= 0.5*w
		for j in xrange(2,numr3i):
			jc = j+numr2i
			circ[jc] *= w

def crit2d(args, data):
	#print  " AMOEBA ",args
	#  data: 0 - kb,  1 - mask, 2 - nima,  3 - current ave, 4 - current image in the gridding format
	#from utilities import info
	from fundamentals import rtshgkb
	mn = data[4].get_attr('mirror')
	temp = rtshgkb(data[4], args[0], args[1], args[2], data[0])
	if  mn: temp.process_inplace("mirror", {"axis":'x'})
	#temp2 = data[3] + temp/data[2]
	temp2 = Util.madn_scalar(data[3], temp, 1.0/data[2]) 
	v = temp2.cmp("dot", temp2, {"negative":0, "mask":data[1]})
	#print  " AMOEBA ",args,mn,v
	return v

def eqproj(args, data):
	from projection import prgs
	#from fundamentals import cyclic_shift
	#from utilities import info
	#print  " AMOEBA ",args
	#  data: 0 - volkb,  1 - kb, 2 - image,  3 - mask,
	#  args: 0 - phi, 1 - theta, 2 - psi, 3 - sx, 4 - sy
	prj = prgs(data[0], data[1], args)

	# the idea is for the mask to "follow" the projection
	#isx = int(args[3]+100000.5)-100000 # this is a strange trick to take care of negative sx
	#isy = int(args[4]+100000.5)-100000
	#shifted_mask = cyclic_shift(data[3], isx, isy)
	#info(proj)
	#info(data[2])
	#info(proj,None,"proj")
	#info(data[2],None,"data[2")
	#info(data[3],None,"data[3")

	#info(shifted_mask,None,"shifted mask")
	#v = -proj.cmp("SqEuclidean", data[2], {"mask":shifted_mask})
	#        CURRENTLY THE DISTANCE IS cross-correlation coefficient
	#v = -prj.cmp("SqEuclidean", data[2], {"mask":data[3]})
	v = prj.cmp("ccc", data[2], {"mask":data[3], "negative":0})
	#v = proj.cmp("ccc", data[2], {"mask":shifted_mask, "negative":0})
	#print  " AMOEBA o", args, v
	return v

def eqprojDot(args, data):
	from projection import project
	from filter import filt_ctf
	phi = args[0]
	tht = args[1]
	psi = args[2]
	vol = data[0]
	img = data[1]
	s2x = data[2]
	s2y = data[3]
	msk = data[4]
	CTF = data[5]
        ou  = data[6]

	tmp = img.process( "normalize.mask", {"mask":msk, "no_sigma":0} )
	ref = project( vol, [phi,tht,psi,-s2x,-s2y], ou )
	if CTF:
		ctf = img.get_attr( "ctf" )
		ref = filt_ctf( ref, ctf )
	return ref.cmp( "dot", tmp, {"mask":msk, "negative":0} )

def eqprojEuler(args, data):
	from projection import prgs
	prj = prgs(data[0], data[1], [args[0], args[1], args[2], data[3], data[4]])
	v = prj.cmp("ccc", data[2], {"mask":data[5], "negative":0})
	return v

def find_symm(vol, mask, sym_gp, phi, theta, psi, scale, ftolerance, xtolerance):
	
	from utilities import amoeba, model_circle
	from alignment import symm_func
	args   = [phi, theta, psi]
	data   = [vol, mask, sym_gp]
	result = amoeba(args, scale, symm_func, ftolerance, xtolerance, 500, data)

	return result

#   Implemented in c in utli_sparx
#  Helper functions for ali2d_r
def kbt(nx,npad=2):
	# padd two times
	N=nx*npad
	# support of the window
	K=6
	alpha=1.75
	r=nx/2
	v=K/2.0/N
	return Util.KaiserBessel(alpha, K, r, K/(2.*N), N)
     

#  AP stuff  01/18/06
    
    
def log2(n):
	""" Returns the smallest power by which 2 has to be raised to obtain
	    an integer less equal n
	"""
	m = 1
	k =-1
	while (m <= n):
		i = m
		k +=1
		m = 2*i
	return k
    
def Numrinit(first_ring, last_ring, skip=1, mode="F"):
	"""This function calculates the necessary information for the 2D 
	   polar interpolation. For each ring, three elements are recorded:
	   numr[i*3]:  Radius of this ring
	   numr[i*3+1]: Total number of samples of all inner rings+1
	   		(Or, the beginning point of this ring)
	   numr[i*3+2]: Number of samples of this ring. This number is an 
	   		FFT-friendly power of the 2.
			
	   "F" means a full circle interpolation
	   "H" means a half circle interpolation
	"""
	MAXFFT=32768
	from math import pi

	if (mode == 'f' or mode == 'F'): dpi = 2*pi
	else:                            dpi = pi

	numr = []
	lcirc = 1
	for k in xrange(first_ring, last_ring+1, skip):
		numr.append(k)
		jp = int(dpi * k+0.5)
		ip = 2**(log2(jp)+1)  # two times oversample each ring
		if (k+skip <= last_ring and jp > ip+ip//2): ip=min(MAXFFT,2*ip)
		if (k+skip  > last_ring and jp > ip+ip//5):  ip=min(MAXFFT,2*ip)

		numr.append(lcirc)
		numr.append(ip)
		lcirc += ip

	lcirc -= 1
	return  numr
	
def ringwe(numr, mode="F"):
	"""Calculate ring weights for rotational alignment
	   Roughly speaking, this weight is inversely proportional to the radius of the ring.
	"""
	from math import pi
	if (mode == 'f' or mode == 'F'):
		dpi = 2*pi
	else:
		dpi = pi
	nring = len(numr)/3
	wr=[0.0]*nring
	maxrin = float(numr[len(numr)-1])
	for i in xrange(0,nring): wr[i] = numr[i*3]*dpi/float(numr[2+i*3])*maxrin/float(numr[2+i*3])
	return wr

def normalize_ring(circ, numr):
	"""
	  Apply weights to FTs of rings
	"""
	from math import sqrt,pi
	av=0.0
	sq=0.0
	nn=0
	nring = len(numr)/3
	maxrin = numr[len(numr)-1]
	for i in xrange(nring):
		numr3i = numr[2+i*3]
		numr2i = numr[1+i*3]-1
		w = numr[3*i]*2*pi/float(numr[3*i+2])
		for j in xrange(numr3i):
			jc = j+numr2i
			av += circ[jc] * w
			sq += circ[jc] * circ[jc]* w
			nn += w
        avg = av/nn
        sgm = sqrt( (sq-av*av/nn)/nn )

        circ -= avg
        circ /= sgm



def ringw_real(numr, mode="F"):
	"""Calculate ring weights for rotational alignment
	   Roughly speaking, this weight is inversely proportional to the radius of the ring.
	"""
	from math import pi
	if (mode == 'f' or mode == 'F'):
		dpi = 2*pi
	else:
		dpi = pi
	nring = len(numr)/3
	wr=[0.0]*nring
	for i in xrange(nring): wr[i] = numr[i*3]*dpi/float(numr[2+i*3])
	return wr



def ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny):
	"""Determine shift and rotation between image and reference image (crefim)
		crefim should be as FT of polar coords with applied weights
	        consider mirror
		quadratic interpolation
		cnx, cny in FORTRAN convention
	"""
	from math import pi, cos, sin
	#print "ORMQ"
	peak = -1.0E23
	ky = int(2*yrng/step+0.5)//2
	kx = int(2*xrng/step+0.5)//2
	for i in xrange(-ky, ky+1):
		iy = i*step
	 	for  j in xrange(-kx, kx+1):
			ix = j*step
			cimage = Util.Polar2Dm(image, cnx+ix, cny+iy, numr, mode)
			Util.Frngs(cimage, numr)
			retvals = Util.Crosrng_ms(crefim, cimage, numr)
			qn = retvals["qn"]
			qm = retvals["qm"]
			if (qn >= peak or qm >= peak):
		 		sx = -ix
		 		sy = -iy
				if (qn >= qm):
			 		ang = ang_n(retvals["tot"], mode, numr[len(numr)-1])
			 		peak = qn
			 		mirror = 0
		   		else:
			 		ang = ang_n(retvals["tmt"], mode, numr[len(numr)-1])
			 		peak = qm
			 		mirror = 1
	co =  cos(ang*pi/180.0)
	so = -sin(ang*pi/180.0)
	sxs = sx*co - sy*so
	sys = sx*so + sy*co
	return  ang, sxs, sys, mirror, peak
			

def ormq_peaks(image, crefim, xrng, yrng, step, mode, numr, cnx, cny):
	"""
	Determine shift and rotation between image and reference image (crefim)
	crefim should be as FT of polar coords with applied weights
	consider mirror
	quadratic interpolation
	cnx, cny in FORTRAN convention
	"""
	from utilities import peak_search

	ccfs = EMData()
	ccfm = EMData()
	Util.multiref_peaks_ali2d(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, ccfs, ccfm)

	peaks = peak_search(ccfs, 1000)
	for i in xrange(len(peaks)):	peaks[i].append(0)

	peakm = peak_search(ccfm, 1000)
	for i in xrange(len(peakm)):	peakm[i].append(1)
	peaks += peakm

	return peaks


'''
def process_peak_1d_pad(peaks, step, mode, numr, nx):
	from math import pi, cos, sin
	
	peak_num = len(peaks)
	maxrin = numr[-1]
	for i in xrange(peak_num):
		peaks[i][1]= ang_n(peaks[i][1]+1-nx/2, mode, maxrin)
	peaks.sort(reverse=True)
	
	return peaks

def find_position(list_a, t):
	"""
	The function determines how many elements of list_a is larger or equal than t.
	Here we assume list_a is sorted reversely.
	"""
	if list_a[0] < t:
		return 0
	elif list_a[-1] >= t:
		return len(list_a)
	else:
		K = len(list_a)
		k_min = 1
		k_max = K-1
		while k_min != k_max:
			k = (k_min+k_max)/2
			if list_a[k] < t:
				if list_a[k-1] >= t:
					k_min = k
					k_max = k
				else:
					k_max = k-1
			else:
				k_min = k+1
		return k_min


def select_major_peaks(g, max_major_peaks, min_height, dim):

	from filter import filt_gaussl
	from fundamentals import fft
	from utilities import peak_search
	
	G = fft(g)
	
	found = False
	min_fl = 0.001
	max_fl = 0.5
	
	while found == False:
		fl = (min_fl+max_fl)/2
		peakg = peak_search(fft(filt_gaussl(G, fl)), 1000)
		K = len(peakg)
		list_a = [0.0]*K
		for i in xrange(K):  list_a[i] = peakg[i][dim+1]
		k = find_position(list_a, min_height)
		if k > max_major_peaks: 
			max_fl = fl
		elif k < max_major_peaks:
			min_fl = fl
		else:
			found = True
		if max_fl - min_fl < 0.001: found = True
		 
	return peakg[0:k] 


def select_major_peaks_Gaussian_fitting(peak):

	# Generate the histogram of the angle
	ang_bin = [0]*30
	for i in xrange(len(angle)):
		#angle.append(peak[i][1])
		bin_num = int(angle[i]/12)
		ang_bin[bin_num] += 1
	ang_bin_index = []
	for i in xrange(30):
		ang_bin_index.append([ang_bin[i], i])
	ang_bin_index.sort(reverse=True)
	print ang_bin
	print ang_bin_index
	
	K = 5
	A = [0.0]*K
	mu = [0.0]*K
	sigma = [0.0]*K
	
	for k in xrange(K):
		A[k] = ang_bin_index[k][0]
		mu[k] = ang_bin_index[k][1]*12+6
		sigma[k] = 5.0
	
	print A, mu, sigma 
	
	
	return []


def ormq_peaks_major(image, crefim, xrng, yrng, step, mode, numr, cnx, cny):
	"""
	Determine shift and rotation between image and reference image (crefim)
	crefim should be as FT of polar coords with applied weights
	consider mirror
	quadratic interpolation
	cnx, cny in FORTRAN convention
	"""
	from utilities import peak_search, pad
	
	ccfs = EMData()
	ccfm = EMData()
	ccfs_compress = EMData()
	ccfm_compress = EMData()

	Util.multiref_peaks_compress_ali2d(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, ccfs, ccfm, ccfs_compress, ccfm_compress)
	
	nx = ccfs.get_xsize()
	ny = ccfs.get_ysize()
	nz = ccfs.get_zsize()

	peaks = peak_search(ccfs, 1000)
	peakm = peak_search(ccfm, 1000)

	peaks = process_peak(peaks, step, mode, numr)
	peakm = process_peak(peakm, step, mode, numr)

	max_major_peaks = 5
	min_height = 0.7
	
	peaks_major = select_major_peaks(pad(ccfs_compress, nx*2), max_major_peaks, min_height, 1)
	peakm_major = select_major_peaks(pad(ccfm_compress, nx*2), max_major_peaks, min_height, 1)

	peaks_major = process_peak_1d_pad(peaks_major, step, mode, numr, nx)
	peakm_major = process_peak_1d_pad(peakm_major, step, mode, numr, nx)	

	"""
	ccfs_compress = EMData(nx, 1, 1, True)
	ccfm_compress = EMData(nx, 1, 1, True)

	for x in xrange(nx):
		slices = Util.window(ccfs, 1, ny-2, nz-2, x-nx/2, 0, 0)
		slicem = Util.window(ccfm, 1, ny-2, nz-2, x-nx/2, 0, 0)
		
		[means, dummy, dummy, dummy] = Util.infomask(slices, None, True)
		[meanm, dummy, dummy, dummy] = Util.infomask(slicem, None, True)
		
		ccfs_compress.set_value_at(x, 0, 0, means)
		ccfm_compress.set_value_at(x, 0, 0, meanm)

	peaks_major = select_major_peaks_Gaussian_fitting(peaks)
	peakm_major = select_major_peaks_Gaussian_fitting(peakm)
	
	fs = Util.window(ccfs, nx, ny-2, nz-2, 0, 0, 0)
	fm = Util.window(ccfm, nx, ny-2, nz-2, 0, 0, 0)

	dummy1, dummy2, mins, dummy3 = Util.infomask(fs, None, True)
	dummy1, dummy2, minm, dummy3 = Util.infomask(fm, None, True)

	gs = threshold_to_minval(ccfs, mins)
	gm = threshold_to_minval(ccfm, minm)

	peaks_major = select_major_peaks(gs, max_major_peaks, min_height, 3)
	peakm_major = select_major_peaks(gm, max_major_peaks, min_height, 3)
	
	peaks_major = select_major_peaks(pad(ccfs_compress, nx*2), max_major_peaks, min_height, 1)
	peakm_major = select_major_peaks(pad(ccfm_compress, nx*2), max_major_peaks, min_height, 1)
	time4 = time()
	peaks = peak_search(ccfs, 1000)
	peakm = peak_search(ccfm, 1000)
	time5 = time()
	peaks = process_peak(peaks, step, mode, numr)
	peakm = process_peak(peakm, step, mode, numr)
	peaks_major = process_peak_1d_pad(peaks_major, step, mode, numr)
	peakm_major = process_peak_1d_pad(peakm_major, step, mode, numr)	
	time6 = time()
	"""
	
	return peaks, peakm, peaks_major, peakm_major
'''


def select_k(dJe, T):
	"""
	This routine is used in simulated annealing to select a random path
	based on the weight of the each path and the temperature.
	"""
	from random import random

	K = len(dJe)

	p  = [0.0] * K
	ut = 1.0/T
	for k in xrange(K): p[k] = dJe[k]**ut

	sumq = float(sum(p))
	for k in xrange(K): p[k] /= sumq
	#print  p

	for k in xrange(1, K-1): p[k] += p[k-1]
	# the next line looks strange, but it assures that at least the lst element is selected
	p[K-1] = 2.0

	pb = random()
	select = 0

	while(p[select] < pb):  select += 1
	#select = 0
	return select

def sim_anneal(peaks, T, step, mode, maxrin):
	from random import random
	from math import pi, cos, sin

	peaks.sort(reverse=True)

	if T < 0.0:
		select = int(-T)
		ang = ang_n(peaks[select][1]+1, mode, maxrin)
		sx  = -peaks[select][6]*step
		sy  = -peaks[select][7]*step

		co =  cos(ang*pi/180.0)
		so = -sin(ang*pi/180.0)
		sxs = sx*co - sy*so
		sys = sx*so + sy*co

		mirror = peaks[select][8]
		peak   = peaks[select][0]/peaks[0][0]
	elif T == 0.0:
		select = 0
	
		ang = ang_n(peaks[select][1]+1, mode, maxrin)
		sx  = -peaks[select][6]*step
		sy  = -peaks[select][7]*step

		co =  cos(ang*pi/180.0)
		so = -sin(ang*pi/180.0)
		sxs = sx*co - sy*so
		sys = sx*so + sy*co

		mirror = peaks[select][8]
		peak   = peaks[select][0]/peaks[0][0]
	else:
		K = len(peaks)
		qt = peaks[0][0]
		p  = [0.0] * K
		ut = 1.0/T
		for k in xrange(K): p[k] = (peaks[k][0]/qt)**ut

		sumq = float(sum(p))
		cp  = [0.0] * K
		for k in xrange(K):
			p[k] /= sumq
			cp[k] = p[k]
		#print  p

		for k in xrange(1, K-1): cp[k] += cp[k-1]
		# the next line looks strange, but it assures that at least the lst element is selected
		cp[K-1] = 2.0

		pb = random()
		select = 0
		while(cp[select] < pb):  select += 1

		ang = ang_n(peaks[select][1]+1, mode, maxrin)
		sx  = -peaks[select][6]*step
		sy  = -peaks[select][7]*step

		co =  cos(ang*pi/180.0)
		so = -sin(ang*pi/180.0)
		sxs = sx*co - sy*so
		sys = sx*so + sy*co

		mirror = peaks[select][8]
		peak   = p[select]
		
	return  ang, sxs, sys, mirror, peak, select

def sim_ccf(peaks, T, step, mode, maxrin):
	from random import random
	from math import pi, cos, sin

	if T < 0.0:
		select = int(-T)
		ang = ang_n(peaks[select][1]+1, mode, maxrin)
		sx  = -peaks[select][2]*step
		sy  = -peaks[select][3]*step

		co =  cos(ang*pi/180.0)
		so = -sin(ang*pi/180.0)
		sxs = sx*co - sy*so
		sys = sx*so + sy*co

		mirror = peaks[select][4]
		peak   = peaks[select][0]/peaks[0][0]
	elif T == 0.0:
		select = 0
	
		ang = ang_n(peaks[select][1]+1, mode, maxrin)
		sx  = -peaks[select][2]*step
		sy  = -peaks[select][3]*step

		co =  cos(ang*pi/180.0)
		so = -sin(ang*pi/180.0)
		sxs = sx*co - sy*so
		sys = sx*so + sy*co

		mirror = peaks[select][4]
		peak   = peaks[select][0]/peaks[0][0]
	else:
		select = int(peaks[5])
		ang = ang_n(peaks[1]+1, mode, maxrin)
		sx  = -peaks[2]*step
		sy  = -peaks[3]*step

		co =  cos(ang*pi/180.0)
		so = -sin(ang*pi/180.0)
		sxs = sx*co - sy*so
		sys = sx*so + sy*co

		mirror = int(peaks[4])
		peak   = peaks[0]

	return  ang, sxs, sys, mirror, peak, select


def sim_anneal2(peaks, Iter, T0, F, SA_stop):
	from math import exp, pow
	from random import random

	# Determine the current temperature
	T = T0*pow(F, Iter)	

	K = len(peaks)
	p = [0.0] * K

	if T > 0.0001 and Iter < SA_stop:
	
		dJe = [0.0]*K
		for k in xrange(K):
			dJe[k] = peaks[k][0]/peaks[0][0]

		# q[k]
		q      = [0.0] * K
		arg    = [0.0] * K
		maxarg = 0
		for k in xrange(K):
			arg[k] = dJe[k] / T
			if arg[k] > maxarg: maxarg = arg[k]
		limarg = 200
		if maxarg > limarg:
			sumarg = float(sum(arg))
			for k in xrange(K): q[k] = exp(arg[k] * limarg / sumarg)
		else:
			for k in xrange(K): q[k] = exp(arg[k])

		sumq = float(sum(q))
		for k in xrange(K):
			p[k] = q[k] / sumq
	else:
		p[0] = 1.0
	
	return p


def sim_anneal3(peaks, peakm, peaks_major, peakm_major, Iter, T0, F, SA_stop):
	from math import pow, sin, sqrt, pi
	from random import random

	# Determine the current temperature
	T = T0*pow(F, Iter)
	max_peak = 5
	DEG_to_RAD = pi/180.0

	dim = 1

	if T > 0.001 and Iter < SA_stop:

		K = len(peaks_major)
		dJe = [0.0]*K
		for k in xrange(K):	dJe[k] = peaks_major[k][dim+1]
		
		select_major = select_k(dJe, T)
		
		ang_m = peaks_major[select_major][1]
		#sx_m = peaks_major[select_major][6]
		#sy_m = peaks_major[select_major][7]
		
		neighbor = []
		for i in xrange(len(peaks)):
			ang = peaks[i][1]
			#sx = peaks[i][6]
			#sy = peaks[i][7]		
			dist = 64*abs(sin((ang-ang_m)/2*DEG_to_RAD))#+sqrt((sx-sx_m)**2+(sy-sy_m)**2)
			neighbor.append([dist, i])
		neighbor.sort()

		dJe = [0.0]*max_peak
		for k in xrange(max_peak):   dJe[k] = peaks[neighbor[k][1]][4]
		select_s = neighbor[select_k(dJe, T)][1]
			
		#############################################################################################################

		K = len(peakm_major)
		dJe = [0.0]*K
		for k in xrange(K): 	dJe[k] = peakm_major[k][dim+1]

		select_major = select_k(dJe, T)
				
		ang_m = peakm_major[select_major][1]
		#sx_m = peakm_major[select_major][6]
		#sy_m = peakm_major[select_major][7]
		
		neighbor = []
		for i in xrange(len(peakm)):
			ang = peakm[i][1]
			#sx = peakm[i][6]
			#sy = peakm[i][7]		
			dist = 64*abs(sin((ang-ang_m)/2*DEG_to_RAD))#+sqrt((sx-sx_m)**2+(sy-sy_m)**2)
			neighbor.append([dist, i])
		neighbor.sort()

		dJe = [0.0]*max_peak
		for k in xrange(max_peak):   dJe[k] = peakm[neighbor[k][1]][4]
		select_m = neighbor[select_k(dJe, T)][1]

		ps = peaks[select_s][0]
		pm = peakm[select_m][0]
		pk = select_k([1.0, min(ps/pm, pm/ps)], T)
		
		if ps > pm and pk == 0 or ps < pm and pk == 1: use_mirror = 0
		else: use_mirror = 1
	else:
		select_s = 0
		select_m = 0
		ps = peaks[select_s][0]
		pm = peakm[select_m][0]
		if ps > pm:
			use_mirror = 0
		else:
			use_mirror = 1
	
	if use_mirror == 0:
		select = select_s	
		ang = peaks[select][1]
		sx  = peaks[select][6]
		sy  = peaks[select][7]
		mirror = 0
		peak = peaks[select][0]
	else:
		select = select_m
		ang = peakm[select][1]
		sx  = peakm[select][6]
		sy  = peakm[select][7]
		mirror = 1
		peak = peakm[select][0]
		
	return  ang, sx, sy, mirror, peak, select


def prep_vol_kb(vol, kb, npad=2):
	# prepare the volume
	volft = vol.copy()
	volft.divkbsinh(kb)
	volft = volft.norm_pad(False, npad)
	volft.do_fft_inplace()
	volft.center_origin_fft()
	volft.fft_shuffle()
	return  volft

def prepare_refrings( volft, kb, delta, ref_a, sym, numr, MPI=False):
        from projection import prep_vol, prgs
        from math       import sin, cos, pi
	from applications import MPI_start_end
	from utilities  import even_angles
	# generate list of Eulerian angles for reference projections
	#  phi, theta, psi
	mode = "F"
	ref_angles = even_angles(delta, symmetry=sym, method = ref_a, phiEqpsi = "Minus")
	wr_four  = ringwe(numr, "F")
	ny = volft.get_ysize()/2
        nx = ny
	cnx = nx//2 + 1
	cny = ny//2 + 1
	qv = pi/180.
        num_ref = len(ref_angles)

	if MPI:
		from mpi import mpi_comm_rank, mpi_comm_size, MPI_COMM_WORLD
		myid = mpi_comm_rank( MPI_COMM_WORLD )
		ncpu = mpi_comm_size( MPI_COMM_WORLD )
	else:
		ncpu = 1
		myid = 0
	from applications import MPI_start_end
	ref_start,ref_end = MPI_start_end( num_ref, ncpu, myid )

	refrings = []     # list of (image objects) reference projections in Fourier representation

	sizex = numr[ len(numr)-2 ] + numr[ len(numr)-1 ] - 1

        for i in xrange(num_ref):
		prjref = EMData()
		prjref.set_size(sizex, 1, 1)
		refrings.append(prjref)

        for i in xrange(ref_start, ref_end):
		prjref = prgs(volft, kb, [ref_angles[i][0], ref_angles[i][1], ref_angles[i][2], 0.0, 0.0])
		cimage = Util.Polar2Dm(prjref, cnx, cny, numr, mode)  # currently set to quadratic....
		normalize_ring(cimage, numr)

		Util.Frngs(cimage, numr)
		Applyws(cimage, numr, wr_four)
		refrings[i] = cimage

	if MPI:
		from utilities import bcast_EMData_to_all
		for i in xrange(num_ref):
			for j in xrange(ncpu):
				ref_start,ref_end = MPI_start_end(num_ref,ncpu,j)
				if i >= ref_start and i < ref_end: rootid = j
			
			bcast_EMData_to_all( refrings[i], myid, rootid )

	
	for i in xrange(len(ref_angles)):
		n1 = sin(ref_angles[i][1]*qv)*cos(ref_angles[i][0]*qv)
		n2 = sin(ref_angles[i][1]*qv)*sin(ref_angles[i][0]*qv)
		n3 = cos(ref_angles[i][1]*qv)
		refrings[i].set_attr_dict( {"n1":n1, "n2":n2, "n3":n3} )
		refrings[i].set_attr("phi", ref_angles[i][0])
		refrings[i].set_attr("theta", ref_angles[i][1])
		refrings[i].set_attr("psi", ref_angles[i][2])




	return refrings

def refprojs( volft, kb, ref_angles, last_ring, mask2D, cnx, cny, numr, mode, wr ):
        from projection import prgs

	ref_proj_rings = []     # list of (image objects) reference projections in Fourier representation
        for i in xrange(len(ref_angles)):
		#prjref = project(volref, [ref_angles[i][0], ref_angles[i][1], ref_angles[i][2], 0.0, 0.0], last_ring)
		prjref = prgs(volft, kb, [ref_angles[i][0], ref_angles[i][1], ref_angles[i][2], 0.0, 0.0])
		prjref.process_inplace("normalize.mask", {"mask":mask2D, "no_sigma":1})  # (i-a)/s
		cimage = Util.Polar2Dm(prjref, cnx, cny, numr, mode)  # currently set to quadratic....
		Util.Frngs(cimage, numr)
		Applyws(cimage, numr, wr)
		ref_proj_rings.append(cimage)

	return ref_proj_rings

def proj_ali_incore(data, refrings, numr, mask2D, xrng, yrng, step):
	from utilities    import even_angles, model_circle, compose_transform2, print_msg
	from alignment    import prepare_refprojs
	from utilities    import get_params_proj, set_params_proj

	mode = "F"
	#  center is in SPIDER convention
	nx   = data.get_xsize()
	ny   = data.get_ysize()
	cnx  = nx//2 + 1
	cny  = ny//2 + 1
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)

	# prepare 2-D mask for normalization
	mask2D = model_circle(last_ring, nx, ny)
	if(first_ring > 0): mask2D -= model_circle(first_ring, nx, ny)


	phi, theta, psi, sxo, syo = get_params_proj( projdata[imn] )
	[ang, sxs, sys, mirror, iref, peak] = Util.multiref_polar_ali_2d(data, refrings, xrng, yrng, step, mode, numr, cnx-sxo, cny-syo)
	iref=int(iref)
	#[ang,sxs,sys,mirror,peak,numref] = apmq(projdata[imn], ref_proj_rings, xrng, yrng, step, mode, numr, cnx-sxo, cny-syo)
	#ang = (ang+360.0)%360.0
	# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
	#  What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
	angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1, -ang, 0.0, 0.0, 1)
	if mirror:
		phi   = (refrings[iref].get_attr("phi")+540.0)%360.0
		theta = 180.0-ref_angles[iref].get_attr("theta")
		psi   = (540.0-ref_angles[iref].get_attr("psi")+angb)%360.0
		s2x   = sxb + sxo
		s2y   = syb + syo
	else:
		phi   = refrings[iref].get_attr("phi")
		theta = refrings[iref].get_attr("theta")
		psi   = (refrings[iref].get_attr("psi")+angb+360.0)%360.0
		s2x   = sxb + sxo
		s2y   = syb + syo
	set_params_proj(data, [phi, theta, psi, s2x, s2y])
	return peak

def proj_ali_incore_local(data, refrings, numr, mask2D, xrng, yrng, step, an):
	from utilities    import even_angles, model_circle, compose_transform2, bcast_EMData_to_all
	from utilities import set_params_proj, get_params_proj
	from math         import cos, sin, pi

	mode = "F"
	#  begin from applying the mask, i.e., subtract the average outside the mask and multiply by the mask
	nx   = data.get_xsize()
	ny   = data.get_ysize()
	#  center is in SPIDER convention
	cnx  = nx//2 + 1
	cny  = ny//2 + 1


	qv = pi/180.
	ant = abs(cos(an*qv))
	phi, theta, psi, sxo, syo = get_params_proj(data)

	[ang, sxs, sys, mirror, iref, peak] = Util.multiref_polar_ali_2d_local(data,refrings,xrng,yrng,step,ant,mode,numr,cnx-sxo,cny-syo)
	iref=int(iref)
	#[ang,sxs,sys,mirror,peak,numref] = apmq_local(projdata[imn], ref_proj_rings, xrng, yrng, step, ant, mode, numr, cnx-sxo, cny-syo)
	#ang = (ang+360.0)%360.0
	if iref > -1:
		# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
		# What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
		angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1,  -ang, 0.,0.,1)
		if  mirror:
			phi   = (refrings[iref].get_attr("phi")+540.0)%360.0
			theta = 180.0-refrings[iref].get_attr("theta")
			psi   = (540.0-refrings[iref].get_attr("psi")+angb)%360.0
			s2x   = sxb + sxo
			s2y   = syb + syo
		else:
			phi   = refrings[iref].get_attr("phi")
			theta = refrings[iref].get_attr("theta")
			psi   = (refrings[iref].get_attr("psi")+angb+360.0)%360.0
			s2x   = sxb+sxo
			s2y   = syb+syo

		from utilities import set_params_proj, get_params_proj
		set_params_proj( data, [phi, theta, psi, s2x, s2y] )
		return peak
	else:
		return -1.0e23


def proj_ali_incore_peaks(volref, mask3D, projdata, first_ring, last_ring, rstep, xrng, yrng, step, delta, ref_a, symmetry, CTF = False, finfo=None, MPI=False):
	from utilities    import even_angles, model_circle, compose_transform2, print_msg
	from alignment    import prepare_refprojs
	from utilities    import get_params_proj, set_params_proj

	mode    = "F"
	# generate list of Eulerian angles for reference projections
	#  phi, theta, psi
	ref_angles = []#even_angles(delta, symmetry = symmetry, method = ref_a, phiEqpsi = "Minus")
	nphi = 0
	for i in xrange(0,360.0,delta): nphi += 1
	ntheta = 0
	for j in xrange(0,91.0+delta, delta):
		ntheta += 1
		for i in xrange(0,360.0,delta):
			ref_angles.append([float(i), float(j), 0.0])
	for i in xrange(len(ref_angles)):
		print ref_angles[i]
	print len(ref_angles)

	#  begin from applying the mask, i.e., subtract the average outside the mask and multiply by the mask
	if mask3D:
		[mean, sigma, xmin, xmax ] =  Util.infomask(volref, mask3D, False)
		volref -= mean
		Util.mul_img(volref, mask3D)
	if CTF:
		from filter    import filt_ctf

	nx   = volref.get_xsize()
	ny   = volref.get_ysize()
	#  center is in SPIDER convention
	cnx  = nx//2 + 1
	cny  = ny//2 + 1
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr   = ringwe(numr ,mode)

	# prepare 2-D mask for normalization
	mask2D = model_circle(last_ring, nx, ny)
	if(first_ring > 0): mask2D -= model_circle(first_ring, nx, ny)

	# generate reference projections in polar coords
	ref_proj_rings = prepare_refprojs( volref, ref_angles, last_ring, mask2D, cnx, cny, numr, mode, wr, MPI )
	print  "  projections prepared"
	for imn in xrange(len(projdata)):
		#if(imn%10 == 0):  print_msg("%d  "%(imn))
		phi, theta, psi, sxo, syo = get_params_proj( projdata[imn] )
		if not(finfo is None):
			finfo.write( "proj %4d old params: %8.3f %8.3f %8.3f %8.3f %8.3f\n" %(imn, phi, theta, psi, sxo, syo) )
			finfo.flush()
		from utilities import ttime
		print ttime()
		ccfs = EMData()
		ccfm = EMData()
		if CTF:
			ctf_params = projdata[imn].get_attr("ctf")
			ima = filt_ctf(projdata[imn], ctf_params, True)
			Util.multiref_peaks_ali(ima.process("normalize.mask", {"mask":mask2D, "no_sigma":1}), ref_proj_rings, xrng, yrng, step, mode, numr, cnx-sxo, cny-syo, ccfs, ccfm, nphi, ntheta)
		else:
			Util.multiref_peaks_ali(projdata[imn].process("normalize.mask", {"mask":mask2D, "no_sigma":1}), ref_proj_rings, xrng, yrng, step, mode, numr, cnx-sxo, cny-syo, ccfs, ccfm, nphi, ntheta)
		print ttime()
		from utilities import info
		info(ccfs)
		info(ccfm)
		from utilities import peak_search
		pp = peak_search(ccfs,1000)
		pm = peak_search(ccfm,1000)
		print ttime()
		print len(pp),pp[0],len(pm),pm[0]
		"""
		numref=int(nref)
		#[ang,sxs,sys,mirror,peak,numref] = apmq(projdata[imn], ref_proj_rings, xrng, yrng, step, mode, numr, cnx-sxo, cny-syo)
		#ang = (ang+360.0)%360.0
		# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
		#  What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
		angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1,  -ang, 0.,0.,1)
		if  mirror:
                        phi   = (ref_angles[numref][0]+540.0)%360.0
                        theta = 180.0-ref_angles[numref][1]
                        psi   = (540.0-ref_angles[numref][2]+angb)%360.0
                        s2x   = sxb + sxo
                        s2y   = syb + syo
                else:
                        phi   = ref_angles[numref][0]
                        theta = ref_angles[numref][1]
                        psi   = (ref_angles[numref][2]+angb+360.0)%360.0
                        s2x   = sxb + sxo
                        s2y   = syb + syo
		set_params_proj( projdata[imn], [phi, theta, psi, s2x, s2y] )
		projdata[imn].set_attr('peak', peak)
		if not(finfo is None):
			finfo.write( "proj %4d new params: %8.3f %8.3f %8.3f %8.3f %8.3f\n" %(imn, phi, theta, psi, s2x, s2y) )
			finfo.flush()
		"""

def ali_vol_func(params, data):
	from utilities    import compose_transform3
	from fundamentals import rot_shift3D
	#print  params
	#print  data[3]
	#cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale= compose_transform3(data[3][0], data[3][1], data[3][2], data[3][3], data[3][4], data[3][5], data[3][6], params[0], params[1], params[2],params[3], params[4], params[5],1.0)
	#rint  cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale
	x = rot_shift3D(data[0], params[0], params[1], params[2], params[3], params[4], params[5], 1.0)
	res = -x.cmp("ccc", data[1], {"mask":data[2]})
	print  " %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f  %10.5f" %(params[0], params[1], params[2],params[3], params[4], params[5], -res)
	return res

def ali_vol_func_rotate(params, data):
	from utilities    import compose_transform3
	from fundamentals import rot_shift3D
	cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale= compose_transform3(data[3][0], data[3][1], data[3][2], data[3][3], data[3][4], data[3][5], data[3][6], params[0], params[1], params[2],0.0,0.0,0.0,1.0)
	x = rot_shift3D(data[0], cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale)
	res = -x.cmp(data[4], data[1], {"mask":data[2]})
	print  " %9.3f %9.3f %9.3f  %12.3e" %(params[0], params[1], params[2], -res)
	return res

def ali_vol_func_shift(params, data):
	from utilities    import compose_transform3
	from fundamentals import rot_shift3D
	cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale= compose_transform3(data[3][0], data[3][1], data[3][2], data[3][3], data[3][4], data[3][5], data[3][6], 0.0,0.0,0.0, params[0], params[1], params[2],1.0)
	x = rot_shift3D(data[0], cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale)
	res = -x.cmp(data[4], data[1], {"mask":data[2]})
	print  " %9.3f %9.3f %9.3f  %12.3e" %(params[0], params[1], params[2], -res)
	return res

def ali_vol_func_scale(params, data):
	from utilities    import compose_transform3
	from fundamentals import rot_shift3D
	cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale= compose_transform3(data[3][0], data[3][1], data[3][2], data[3][3], data[3][4], data[3][5], data[3][6], params[0], params[1], params[2], params[3], params[4], params[5], params[6])
	x = rot_shift3D(data[0], cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale)
	res = -x.cmp(data[4], data[1], {"mask":data[2]})
	print  " %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f  %12.3e" %(params[0], params[1], params[2],params[3], params[4], params[5], params[6], -res)
	return res

def ali_vol_func_only_scale(params, data):
	from utilities    import compose_transform3
	from fundamentals import rot_shift3D
	cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale= compose_transform3(data[3][0], data[3][1], data[3][2], data[3][3], data[3][4], data[3][5], data[3][6], 0.0,0.0,0.0,0.0,0.0,0.0, params[0])
	x = rot_shift3D(data[0], cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale)
	res = -x.cmp(data[4], data[1], {"mask":data[2]})
	print  " %9.3f  %12.3e" %(params[0], -res)
	return res

def helios_func(params, data):
	sm = data[0].helicise(data[2], params[0], params[1], data[3], data[4])
	q = sm.cmp("dot", sm, {"negative":0})
	#print  params,q
	return  q

def helios(vol, pixel_size, dp, dphi, section_use = 1.0, radius = 0.0):
	from alignment    import helios_func
	from utilities    import amoeba
	nx = vol.get_xsize()
	ny = vol.get_ysize()
	nz = vol.get_zsize()
	if(radius <= 0.0):    radius = nx//2-1
	params = [dp, dphi]
	#print  " input params ",params
	data=[vol, params, pixel_size, section_use, radius]
	new_params = [dp, dphi]
	new_params = amoeba(new_params, [0.05*dp, 0.05*abs(dphi)], helios_func, 1.0e-3, 1.0e-3, 500, data)
	#print  " new params ", new_params[0], new_params[1]
	return  vol.helicise(pixel_size, new_params[0][0], new_params[0][1], section_use, radius), new_params[0][0], new_params[0][1]

def sub_favj(ave, data, jtot, mirror, numr):
	'''
		Subtract FT of rings from the average
	'''
	from math import pi,sin,cos
	#from utilities  import print_col
	# trig functions in radians
	pi2 = pi*2
	nring = len(numr)/3
	maxrin = numr[len(numr)-1]
	#print  "update",psi
	#print_col(ave)
	#print_col(data)
	#print numr
	if(mirror):
		# for mirrored data has to be conjugated
		for i in xrange(nring):
			numr3i = numr[2+i*3]
			np = numr[1+i*3]-1
			ave[np]   -= data[np]
			ave[np+1] -= data[np+1]*cos(pi2*(jtot-1)/2.0*numr3i/maxrin)
			for j in xrange(2, numr3i, 2):
				arg = pi2*(jtot-1)*int(j/2)/maxrin
				cs = complex(data[np + j],data[np + j +1])*complex(cos(arg),sin(arg))
				ave[np + j]    -= cs.real
				ave[np + j +1] += cs.imag
	else:
		for i in xrange(nring):
			numr3i = numr[2+i*3]
			np = numr[1+i*3]-1
			ave[np]   -= data[np]
			ave[np+1] -= data[np+1]*cos(pi2*(jtot-1)/2.0*numr3i/maxrin)
			for j in xrange(2, numr3i, 2):
				arg = pi2*(jtot-1)*int(j/2)/maxrin
				cs = complex(data[np + j],data[np + j +1])*complex(cos(arg),sin(arg))
				ave[np + j]    -= cs.real
				ave[np + j +1] -= cs.imag
	#print_col(ave)

def symm_func(args, data):
	from utilities import sym_vol
	from fundamentals  import  rot_shift3D
	sym = sym_vol(rot_shift3D(data[0], args[0], args[1], args[2]), data[2])
	avg = sym.cmp("dot",sym,{"mask":data[1], "negative":0})
	print avg, args
	return avg

def update_favj(ave, data, jtot, mirror, numr):
	'''
		Add FT of rings to the average
	'''
	from math import pi,sin,cos
	#from utilities  import print_col
	# trig functions in radians
	pi2 = pi*2
	nring = len(numr)/3
	maxrin = numr[len(numr)-1]
	#print  "update",psi
	#print_col(ave)
	#print_col(data)
	#print numr
	if(mirror):
		# for mirrored data has to be conjugated
		for i in xrange(nring):
			numr3i = numr[2+i*3]
			np = numr[1+i*3]-1
			ave[np]   += data[np]
			ave[np+1] += data[np+1]*cos(pi2*(jtot-1)/2.0*numr3i/maxrin)
			for j in xrange(2, numr3i, 2):
				arg = pi2*(jtot-1)*int(j/2)/maxrin
				cs = complex(data[np + j],data[np + j +1])*complex(cos(arg),sin(arg))
				ave[np + j]    += cs.real
				ave[np + j +1] -= cs.imag
	else:
		for i in xrange(nring):
			numr3i = numr[2+i*3]
			np = numr[1+i*3]-1
			ave[np]   += data[np]
			ave[np+1] += data[np+1]*cos(pi2*(jtot-1)/2.0*numr3i/maxrin)
			for j in xrange(2, numr3i, 2):
				arg = pi2*(jtot-1)*int(j/2)/maxrin
				cs = complex(data[np + j],data[np + j +1])*complex(cos(arg),sin(arg))
				ave[np + j]    += cs.real
				ave[np + j +1] += cs.imag
	#print_col(ave)

def fine_2D_refinement(data, br, mask, tavg, group = -1):
	from utilities import amoeba
	from fundamentals 	import rtshgkb, prepg

	# IMAGES ARE SQUARES!
	nx = data[0].get_xsize()
	#  center is in SPIDER convention
	cnx = int(nx/2)+1
	cny = cnx

	if(group > -1):
		nima = 0
		for im in xrange(len(data)):
			if(data[im].get_attr('ref_num') == group):  nima += 1
	else:  nima = len(data)

	# prepare KB interpolants
	kb = kbt(nx)
	# load stuff for amoeba
	stuff = []
	stuff.insert(0, kb)
	stuff.insert(1, mask)
	stuff.insert(2, nima)
	#stuff.insert(3,tave)  # current average
	#stuff.insert(4,data)  # current image in the gridding format
	weights = [br]*3 # weights define initial bracketing, so one would have to figure how to set them correctly

	for im in xrange(len(data)):
		if(group > -1):
			if(data[im].get_attr('ref_num') != group):  continue
		# subtract current image from the average
		alpha  = data[im].get_attr('alpha')
		sx     = data[im].get_attr('sx')
		sy     = data[im].get_attr('sy')
		mirror = data[im].get_attr('mirror')
		ddata  = prepg(data[im], kb)
		ddata.set_attr_dict({'alpha': alpha, 'sx':sx, 'sy':sy, 'mirror': mirror})
		temp   = rtshgkb(ddata, alpha, sx, sy, kb)
		if  mirror: temp.process_inplace("mirror", {"axis":'x'})
		#  Subtract current image from the average
		refim = Util.madn_scalar(tavg, temp, -1.0/float(nima)) 
		stuff.append(refim)  # curent ave-1
		stuff.append(ddata)  # curent image
		# perform amoeba alignment
		params = [alpha, sx, sy]
		outparams =  amoeba(params, weights, crit2d, 1.e-4, 1.e-4, 500, stuff)
		del stuff[3]
		del stuff[3]
		# set parameters to the header
		data[im].set_attr_dict({'alpha':outparams[0][0], 'sx':outparams[0][1], 'sy':outparams[0][2],'mirror': mirror})
		# update the average
		temp = rtshgkb(ddata, outparams[0][0], outparams[0][1], outparams[0][2], kb)
		if  mirror: temp.process_inplace("mirror",{"axis":'x'})
		#check whether the criterion actually increased
		# add current image to the average
		tavg = Util.madn_scalar(refim, temp, 1.0/float(nima))
		#print  im,tave.cmp("dot", tave, {"negative":0,"mask":mask}),params,outparams[0],outparams[2]
		#tave,tvar = ave_var_series_g(data,kb)
		#print  " Criterium on the fly ", tave.cmp("dot", tave, {"negative":0,"mask":mask})

def ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny):
	"""Determine shift and rotation between image and reference image (refim)
	   no mirror
		quadratic interpolation
		cnx, cny in FORTRAN convention
	"""
	from math import pi, cos, sin
	from alignment import ang_n
	#from utilities import info
	#print "ORNQ"
	peak = -1.0E23
	ky = int(2*yrng/step+0.5)//2
	kx = int(2*xrng/step+0.5)//2
	
	
	for i in xrange(-ky,ky+1):
		iy=i*step
		for j in xrange(-kx,kx+1):
			ix=j*step
			cimage=Util.Polar2Dm(image, cnx+ix, cny+iy, numr, mode)
			Util.Frngs(cimage, numr)
			retvals = Util.Crosrng_e(crefim, cimage, numr, 0)
			qn = retvals["qn"]
			if(qn >= peak):
				sx = -ix
				sy = -iy
				ang = ang_n(retvals["tot"], mode, numr[len(numr)-1])
				peak=qn
	# mirror is returned as zero for consistency
	mirror=0
	co =  cos(ang*pi/180.0)
	so = -sin(ang*pi/180.0)
	sxs = sx*co - sy*so
	sys = sx*so + sy*co
	return  ang, sxs, sys, mirror, peak
       
def align2d(image, refim, xrng=0, yrng=0, step=1, first_ring=1, last_ring=0, rstep=1, mode = "F"):
	"""  Determine shift and rotation between image and reference image
	     quadratic interpolation
	"""
	#from utilities import print_col
	from alignment import Numrinit, ringwe, Applyws
	step = float(step)
	nx = refim.get_xsize()
	ny = refim.get_ysize()
	if(last_ring == 0):  last_ring = nx/2-2-int(max(xrng,yrng))
	# center in SPIDER convention
	cnx = nx//2+1
	cny = ny//2+1
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr   = ringwe(numr, mode)
	#cimage=Util.Polar2Dmi(refim, cnx, cny, numr, mode, kb)
	crefim = Util.Polar2Dm(refim, cnx, cny, numr, mode)
	#crefim = Util.Polar2D(refim, numr, mode)
	#print_col(crefim)
	Util.Frngs(crefim, numr)
	Applyws(crefim, numr, wr)
	return ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny)

def align2d_peaks(image, refim, xrng=0, yrng=0, step=1, first_ring=1, last_ring=0, rstep=1, mode = "F"):
	"""  Determine shift and rotation between image and reference image
	     quadratic interpolation
	"""
	#from utilities import print_col
	from alignment import Numrinit, ringwe, Applyws
	step = float(step)
	nx = refim.get_xsize()
	ny = refim.get_ysize()
	if(last_ring == 0):  last_ring = nx/2-2-int(max(xrng,yrng))
	# center in SPIDER convention
	cnx = nx//2+1
	cny = ny//2+1
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr   = ringwe(numr, mode)
	#cimage=Util.Polar2Dmi(refim, cnx, cny, numr, mode, kb)
	crefim = Util.Polar2Dm(refim, cnx, cny, numr, mode)
	#crefim = Util.Polar2D(refim, numr, mode)
	#print_col(crefim)
	Util.Frngs(crefim, numr)
	Applyws(crefim, numr, wr)
	return ormq_peaks(image, crefim, xrng, yrng, step, mode, numr, cnx, cny)

def align2d_g(image, refim, xrng=0, yrng=0, step=1, first_ring=1, last_ring=0, rstep=1, mode = "F"):
	"""  Determine shift and rotation between image and reference image
	     quadratic interpolation
	"""
	from development import ormy2
	from alignment import Numrinit, ringwe, Applyws
	from fundamentals import fft
	
	step = float(step)
	nx = refim.get_xsize()
	ny = refim.get_ysize()
	if(last_ring == 0):  last_ring = nx/2-2-int(max(xrng,yrng))
	# center in SPIDER convention
	cnx = int(nx/2)+1
	cny = int(ny/2)+1
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr = ringwe(numr, mode)

	N = nx*2
	K = 6
	alpha = 1.75
	r = nx/2
	v = K/2.0/N
	kb = Util.KaiserBessel(alpha, K, r, v, N)
	refi = refim.FourInterpol(N,N,1,0)  
	params = {"filter_type" : Processor.fourier_filter_types.KAISER_SINH_INVERSE,"alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
	q = Processor.EMFourierFilter(refi,params)
	refi = fft(q)
	crefim = Util.Polar2Dmi(refi,cnx,cny,numr,mode,kb)

	Util.Frngs(crefim, numr)
	Applyws(crefim, numr, wr)
	numr = Numrinit(first_ring, last_ring, rstep, mode)

	return ormy2(image,refim,crefim,xrng,yrng,step,mode,numr,cnx,cny,"gridding")

def max_pixel_error(alpha1, sx1, sy1, alpha2, sx2, sy2, d):
	from math import sin, pi, sqrt
	
	return abs(sin((alpha1-alpha2)/180.0*pi/2))*d+sqrt((sx1-sx2)**2+(sy1-sy2)**2)
