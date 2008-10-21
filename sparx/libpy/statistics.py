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

def add_oe_series(data):
	"""
		Calculate odd and even sum of an image series using current alignment parameters
	"""
	from utilities    import model_blank, get_params2D
	from fundamentals import rot_shift2D
	n = len(data)
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	ave1 = model_blank(nx,ny)
	ave2 = model_blank(nx,ny)
	for i in xrange(n):
		alpha, sx, sy, mirror, scale = get_params2D(data[i])
		temp = rot_shift2D(data[i], alpha, sx, sy, mirror, scale, "quadratic")
		if i%2 == 0: Util.add_img(ave1, temp)
		else:          Util.add_img(ave2, temp)
	return ave1, ave2

def add_ave_varf(data, mask = None, mode = "a", CTF = False, ctf_2_sum = None):
	"""
		Calculate average of an image series and variance, sum of squares in Fourier space
		mode - "a": use current alignment parameters
		CTF  - if True, use CTF for calculations of both average and variance.
	"""
	from utilities    import    model_blank, get_params2D
	from fundamentals import    rot_shift2D, fft

	n = len(data)
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	ave   = model_blank(nx, ny)
	sumsq = EMData(nx, ny, 1, False)
	var   = EMData(nx, ny, 1, False)

	if CTF:
		from morphology   import ctf_img
		from filter       import filt_ctf, filt_table
		from utilities    import get_arb_params
		parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
		if data[0].get_attr_default('ctf_applied', 1) == 1:
			ERROR("data cannot be ctf-applied","add_ave_varf",1)
		if ctf_2_sum:  get_ctf2 = False
		else:          get_ctf2 = True
		if get_ctf2: ctf_2_sum = EMData(nx, ny, 1, False)
	 	for i in xrange(n):
	 		ima = data[i].copy()
	 		if mode == "a":
				alpha, sx, sy, mirror, scale = get_params2D(ima)
				ima = rot_shift2D(ima, alpha, sx, sy, mirror, scale, "quadratic")
				#  Here we have a possible problem: varf works only if CTF is applied after rot/shift
				#    while calculation of average (and in general principle) CTF should be applied before rot/shift
				#    here we use the first possibility
	 		ctf_params = get_arb_params(ima, parnames)
	 		oc = filt_ctf(ima, ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5], pad=False)
			Util.add_img(ave, oc)
 			Util.add_img2(var, fft(ima))
	 		if get_ctf2: Util.add_img2(ctf_2_sum, ctf_img(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5]))
		sumsq = fft(ave)
		ave = fft(Util.divn_img(sumsq, ctf_2_sum))
		Util.mul_img(sumsq, sumsq.conjg())
		Util.div_img(sumsq, ctf_2_sum)
	 	Util.sub_img(var, sumsq)
	else:
		for i in xrange(n):
			if mode == "a":
				alpha, sx, sy, mirror, scale = get_params2D(data[i])
				ima = rot_shift2D(data[i], alpha, sx, sy, mirror, scale, "quadratic")
			else:
				ima = data[i].copy()
			Util.add_img(ave, ima)
			fim = fft(ima)
			Util.add_img2(var, fim)
		sumsq = fft(ave)
		Util.mul_scalar(ave, 1.0/float(n))
		Util.mul_img(sumsq, sumsq.conjg())
		Util.mul_scalar(sumsq, 1.0/float(n))
		Util.sub_img(var, sumsq)
		
	Util.mul_scalar(var, 1.0/float(n-1))	
	st = Util.infomask(var, None, True)
	if st[2]<0.0:  ERROR("Negative variance!", "add_oe_ave_varf", 1)
	return ave, var, sumsq
	

def add_ave_varf_MPI(data, mask = None, mode = "a", CTF = False):
	"""
		Calculate sum of an image series and sum of squares in Fourier space
		Since this is the MPI version, we need to reduce sum and sum of squares 
		on the main node and calculate variance there.
		mode - "a": use current alignment parameters
		CTF  - if True, use CTF for calculations of the sum.
	"""
	from utilities    import    model_blank, get_params2D
	from fundamentals import    rot_shift2D, fft

	n = len(data)
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	ave = model_blank(nx, ny)
	var = EMData(nx, ny, 1, False)

	if CTF:
		from morphology   import ctf_img
		from filter       import filt_ctf, filt_table
		from utilities    import get_arb_params
		parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
		if data[0].get_attr_default('ctf_applied', 1) == 1:
			ERROR("data cannot be ctf-applied","add_ave_varf_MPI",1)
	 	for i in xrange(n):
	 		if mode == "a":
				alpha, sx, sy, mirror, scale = get_params2D(data[i])
				ima = rot_shift2D(data[i], alpha, sx, sy, mirror, scale, "quadratic")
				#  Here we have a possible problem: varf works only if CTF is applied after rot/shift
				#    while calculation of average (and in general principle) CTF should be applied before rot/shift
				#    here we use the first possibility
			else:
				ima = data[i].copy()
	 		ctf_params = get_arb_params(ima, parnames)
	 		oc = filt_ctf(ima, ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5], pad=False)
			Util.add_img(ave, oc)
 			Util.add_img2(var, fft(ima))
	else:
		for i in xrange(n):
			if mode == "a":
				alpha, sx, sy, mirror, scale = get_params2D(data[i])
				ima = rot_shift2D(data[i], alpha, sx, sy, mirror, scale, "quadratic")
			else:
				ima = data[i].copy()
			Util.add_img(ave, ima)
			Util.add_img2(var, fft(ima))
	
	return ave, var
	

def add_ave_varf_ML_MPI(data, mask = None, mode = "a", CTF = False):
	"""
		Calculate sum of an image series and sum of squares in Fourier space
		Since this is the MPI version, we need to reduce sum and sum of squares 
		on the main node and calculate variance there.
		mode - "a" use current alignment parameters
		CTF  - if True, use CTF for calculations of the sum.
	"""
	from utilities    import    model_blank, get_params2D
	from fundamentals import    rot_shift2D, fft

	n = len(data)
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	ave = model_blank(nx, ny)
	var = EMData(nx, ny, 1, False)

	if CTF:
		from morphology   import ctf_img
		from filter       import filt_ctf, filt_table
		from utilities    import get_arb_params
		parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
		if data[0].get_attr_default('ctf_applied', 1) == 1: 
			ERROR("data cannot be ctf-applied","add_oe_ave_varf",1)
	 	for i in xrange(n):
	 		if mode == "a":
				attr_list = data[i].get_attr_dict()	
				if attr_list.has_key('Prob'):
					p = data[i].get_attr('Prob')
					K = len(p)
		 			alpha  = data[i].get_attr('alpha_l')
		 			sx     = data[i].get_attr('sx_l')
	 				sy     = data[i].get_attr('sy_l')
	 				mirror = data[i].get_attr('mirror_l')
		 			ima_tmp = rot_shift2D(data[i], alpha[0], sx[0], sy[0], mirror[0])
					Util.mul_scalar(ima_tmp, p[0])
					for k in xrange(1, K):
						ima_tmp2 = rot_shift2D(data[i], alpha[k], sx[k], sy[k], mirror[k])
						Util.mul_scalar(ima_tmp2, p[k])
						Util.add_img(ima_tmp, ima_tmp2)
					ima = ima_tmp.copy()
				else:
	 				alpha, sx, sy, mirror, scale = get_params2D(data[i])
		 			ima = rot_shift2D(data[i], alpha, sx, sy, mirror)
				#  Here we have a possible problem: varf works only if CTF is applied after rot/shift
				#    while calculation of average (and in general principle) CTF should be applied before rot/shift
				#    here we use the first possibility
			else:
		 		ima = data[i].copy()
	 		ctf_params = get_arb_params(ima, parnames)
	 		oc = filt_ctf(ima, ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5], pad=False)
			Util.add_img(ave, oc)
 			Util.add_img2(var, fft(ima))
	else:
		for i in xrange(n):
			if mode == "a":
	 			alpha, sx, sy, mirror, scale = get_params2D(data[i])
		 		ima = rot_shift2D(data[i], alpha, sx, sy, mirror)
			else:
				ima = data[i].copy()
			Util.add_img(ave, ima)
			Util.add_img2(var, fft(ima))
	
	return ave, var

def ave_var_s(data):
	"""
		Calculate average and variance of an image series
	"""
	from utilities import model_blank
	n = len(data)
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	nz = data[0].get_zsize()
	ave = model_blank(nx,ny,nz)
	var = model_blank(nx,ny,nz)
	for i in xrange(n):
		Util.add_img(ave, data[i])
		Util.add_img2(var, data[i])
	ave /= n
	return ave, (var - ave*ave*n)/(n-1)

def add_oe(data):
	"""
		Calculate odd and even sum of an image series
	"""
	from utilities import model_blank
	n = len(data)
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	nz = data[0].get_zsize()
	ave1 = model_blank(nx,ny,nz)
	ave2 = model_blank(nx,ny,nz)
	for i in xrange(n):
		if i%2 == 0: Util.add_img(ave1, data[i])
		else:        Util.add_img(ave2, data[i])
	return ave1, ave2

def ave_series(data):
	"""
		Calculate average of a image series using current alignment parameters
		data - real space image series
	"""
	from utilities    import model_blank, get_params2D
	from fundamentals import rot_shift2D
	n = len(data)
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	ave = model_blank(nx,ny)
	for i in xrange(n):
	 	alpha, sx, sy, mirror, scale = get_params2D(data[i])
		temp = rot_shift2D(data[i], alpha, sx, sy, mirror)
		Util.add_img(ave, temp)

	return ave/n

def ave_series_ctf(data, ctf2):
	"""
		Calculate average of an image series using current alignment parameters and ctf
		data - real space image series premultiplied by the CTF
	"""
	from utilities    import model_blank, get_params2D
	from filter       import filt_table
	from fundamentals import rot_shift2D
	n = len(data)
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	ave = model_blank(nx,ny)
	for i in xrange(n):
	 	alpha, sx, sy, mirror, scale = get_params2D(data[i])
		temp = rot_shift2D(data[i], alpha, sx, sy, mirror)
		Util.add_img(ave, temp)

	return filt_table(ave, ctf2)

def ave_var_series(data, kb):
	"""
		Calculate average and variance of an image series using current alignment parameters
	"""
	from fundamentals import rotshift2dg
	from utilities    import model_blank
	n = len(data)
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	ave = model_blank(nx,ny)
	var = model_blank(nx,ny)
	for i in xrange(n):
		alpha = data[i].get_attr('alpha')
		sx    =  data[i].get_attr('sx')
		sy    =  data[i].get_attr('sy')
		mirror =  data[i].get_attr('mirror')
		temp = rotshift2dg(data[i], alpha, sx, sy, kb)
		if  mirror: temp.process_inplace("mirror",{"axis":'x'})
		Util.add_img(ave, temp)
		Util.add_img2(var, temp)

	ave /= n
	return ave, (var - ave*ave*n)/(n-1)

def ave_var_series_g(data, kb):
	"""
		Calculate average and variance of a image series using current alignment parameters,
		data contains images prepared for gridding
	"""
	from fundamentals import rtshgkb
	from utilities    import model_blank
	n = len(data)
	ny = data[0].get_ysize()/2
	nx = ny
	ave = model_blank(nx,ny)
	var = model_blank(nx,ny)
	for i in xrange(n):
		alpha = data[i].get_attr('alpha')
		sx =  data[i].get_attr('sx')
		sy =  data[i].get_attr('sy')
		mirror =  data[i].get_attr('mirror')
		temp = rtshgkb(data[i], alpha, sx, sy, kb)
		if  mirror: temp.process_inplace("mirror",{"axis":'x'})
		Util.add_img(ave, temp)
		Util.add_img2(var, temp)

	ave /= n
	return ave, (var - ave*ave*n)/(n-1)
	
def ave_oe_series_g(data, kb):
	"""
		Calculate odd and even averages of a image series using current alignment parameters,
		      data contains images prepared for gridding
	"""
	from fundamentals import rtshgkb
	from utilities    import model_blank
	n  = len(data)
	ny = data[0].get_ysize()/2
	nx = ny
	ave1 = model_blank(nx, ny)
	ave2 = model_blank(nx, ny)
	for i in xrange(n):
		alpha  =  data[i].get_attr('alpha')
		sx     =  data[i].get_attr('sx')
		sy     =  data[i].get_attr('sy')
		mirror =  data[i].get_attr('mirror')
		temp = rtshgkb(data[i], alpha, sx, sy, kb)
		if  mirror: temp.process_inplace("mirror", {"axis":'x'})
		if i%2 == 0: Util.add_img(ave1, temp)
		else:        Util.add_img(ave2, temp)
	ave1 /= (n/2+(n%2))
	ave2 /= (n/2)
	return ave1, ave2


def ave_oe_series_d(data):
	"""
		Calculate odd and even averages of an image series		      
	"""
	n  = len(data)
	ave1 = data[0].copy()
	ave2 = data[1].copy()
	for i in xrange(2,n):
		if i%2 == 0: Util.add_img(ave1, data[i])
 		else:        Util.add_img(ave2, data[i])
	return ave1/(n//2+(n%2)), ave2/(n//2)


def ave_oe_series(stack):
	"""
		Calculate odd and even averages of an image stack using current alignment parameters
	"""
	from utilities import model_blank, get_params2D
	from fundamentals import rot_shift2D
	n = EMUtil.get_image_count(stack)
	ima = EMData()
	ima.read_image(stack, 0, True)
	nx = ima.get_xsize()
	ny = ima.get_ysize()
	ave1 = model_blank(nx,ny)
	ave2 = model_blank(nx,ny)
	for i in xrange(n):
		ima = EMData()
		ima.read_image(stack,i)
		alpha, sx, sy, mirror, scale = get_params2D(ima)
		temp = rot_shift2D(ima, alpha, sx, sy, mirror)
		if i%2 == 0: Util.add_img(ave1, temp)
		else:        Util.add_img(ave2, temp)
	return ave1/(n//2+(n%2)), ave2/(n//2)


def ave_oe_series_indexed(stack, idx_ref):
	"""
		Calculate odd and even averages of an image series using current alignment parameters,
	"""
	from utilities import model_blank
	from fundamentals import rot_shift2D
	
	ntot = 0
	n = EMUtil.get_image_count(stack)
	ima = EMData()
	ima.read_image(stack,0)
	nx = ima.get_xsize()
	ny = ima.get_ysize()
	ave1 = model_blank(nx,ny)
	ave2 = model_blank(nx,ny)
	for i in xrange(n):
		if i == 0: ima = EMData()
		ima.read_image(stack, i)
		if idx_ref == ima.get_attr('ref_num'):
			ntot+=1
			alpha, sx, sy, mirror, scale = get_params2D(ima)
		 	temp = rot_shift2D(ima, alpha, sx, sy, mirror)
			if i%2 == 0: Util.add_img(ave1, temp)
			else:        Util.add_img(ave2, temp)
	if ntot >= 0:	return ave1/(ntot/2+(ntot%2)), ave2/(ntot/2), ntot
	else:		return ave1, ave2, ntot
	
def ave_var_series_one(data, skip, kb):
	"""
		Calculate average and variance of an image series using current alignment parameters
	"""
	from fundamentals import rotshift2dg
	n = len(data)
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	ave = model_blank(nx,ny)
	var = model_blank(nx,ny)
	for i in xrange(n):
		if( i!=skip):
			alpha  = data[i].get_attr('alpha')
			sx     = data[i].get_attr('sx')
			sy     = data[i].get_attr('sy')
			mirror = data[i].get_attr('mirror')
			temp = rotshift2dg(data[i], alpha, sx, sy, kb)
			if  mirror: temp.process_inplace("mirror", {"axis":'x'})
			Util.add_img(ave, temp)
			Util.add_img2(var, temp)

	ave /= n-1
	return ave, (var - ave*ave*(n-1))/(n-2)

def add_series(stack, i1=0 ,i2=0):
	""" Calculate average and variance files for an image series

	Usage:  average,variance = add_series(stack,i1,i2)
	  i1 - first file in image series
	  i2 - last file in image series
	  average and variance are output objects
	  
	"""
	from utilities import model_blank, get_im

	if(i2==0): i2 = EMUtil.get_image_count(stack)-1
	ave = get_im(stack, i1)
	var = ave*ave  #pow(ave,2.0)
	nx = ave.get_xsize()
	ny = ave.get_ysize()
	nz = ave.get_zsize()

	# process the remaining files
	for index in xrange(i1+1,i2+1):
		e = get_im(stack, index)
		Util.add_img(ave, e)        #ave += e
		Util.add_img2(var, e)       #var += e*e  #pow(e,2.0)

	ii=i2-i1+1
	ave = Util.mult_scalar(ave, 1.0/float(ii))  
	#ave /= ii
	#var = (var - pow(ave,2)/ii)/(ii-1)
	#var = (var-ave*ave*ii)/(ii-1)
	e = model_blank(nx, ny, nz)
	Util.add_img2(e, ave)
	var = Util.madn_scalar(var, e, -float(ii))
	Util.mul_scalar(var, 1.0/float(ii-1))
	#var -= ave*ave*ii
	#var /= (ii-1)
	
	return ave, var

def add_series_class(stack, i1 = 0, i2 = 0):
	""" Calculate average and variance files for each group in an image series

	Usage:  average,variance = add_series(stack,i1,i2)
	  i1 - first file in image series
	  i2 - last file in image series
	  average and variance are output objects
	  
	"""
	from utilities import model_blank, get_im
	if(i2==0): i2 = EMUtil.get_image_count(stack)-1
	e = get_im(stack, i1)
	kc = e.get_attr('nclass')
	nx = e.get_xsize()
	ny = e.get_ysize()
	nz = e.get_zsize()
	ave = []
	var = []
	e = model_blank(nx,ny,nz)
	for k in xrange(kc):
		ave.append(e.copy())
		var.append(e.copy())
	
	nclass = [0]*kc
	# process files
	for index in xrange(i1,i2+1):
		e = get_im(stack, index)
		g = e.get_attr('ref_num')
		nclass[g] += 1
		Util.add_img(ave[g], e)
		#ave[g] += e
		#ave[g] = ave[g] + e
		Util.add_img2(var[g], e)
		#var[g] = var[g] + e*e  #pow(e,2.0)

	for k in xrange(kc):
		ii = nclass[k]
		if(ii > 0):
			ave[k] = Util.mult_scalar(ave[k], 1.0/float(ii))         #ave[k] = ave[k]/ii
			if(ii > 1):
				#var[k] = (var[k] - ave[k]*ave[k]*ii) / (ii-1)
				temp = model_blank(nx, ny, nz)
				Util.add_img2(temp, ave[k])
				var[k] = Util.madn_scalar(var[k], temp, -float(ii))
				Util.mul_scalar(var[k], 1.0/float(ii-1))
			else:
				var[k] = model_blank(nx,ny,nz)

	return ave, var, nclass

def add_series_class_mem(data, assign, kc):
	""" Calculate average and variance files for each group in an image series

	Usage:  average,variance = add_series(data, assign, kc)
		data   - list of images
		assign - list of group assignments
		kc     - number of groups

	  average and variance are output objects
	  
	"""
	from utilities import model_blank

	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	nz = data[0].get_zsize()
	ave = []
	var = []
	e = model_blank(nx,ny,nz)
	for k in xrange(kc):
		ave.append(e.copy())
		var.append(e.copy())
	
	nclass = [0]*kc
	# process files
	for index in xrange(len(data)):
		g = assign[index]
		nclass[g] += 1
		Util.add_img(ave[g], data[index])
		#ave[g] += e
		#ave[g] = ave[g] + e
		Util.add_img2(var[g], data[index])        
		#var[g] = var[g] + e*e  #pow(e,2.0)

	for k in xrange(kc):
		ii = nclass[k]
		ave[k] = Util.mult_scalar(ave[k], 1.0/float(ii))         #ave[k] = ave[k]/ii
		if(ii > 1):
			#var[k] = (var[k] - ave[k]*ave[k]*ii) / (ii-1)
			temp = model_blank(nx, ny, nz)
			Util.add_img2(temp, ave[k])
			var[k] = Util.madn_scalar(var[k], temp, -float(ii))
			Util.mul_scalar(var[k], 1.0/float(ii-1))
		else:
			var[k].to_zero()

	return ave, var, nclass

def add_series_class_ctf(images, ctf1, ctf2, snr, assign, kc):
	""" Calculate average and variance files for each group in an image series
		taking into account CTF information

	Usage:  average,variance = add_series(images, assign, kc)
		images   - list of images
		assign - list of group assignments
		kc     - number of groups

	  average and variance are output objects
	  
	"""
	from fundamentals import  fftip, fft
	from filter       import  filt_table
	from utilities    import  model_blank #, info
	nx = images[0].get_xsize()
	ny = images[0].get_ysize()
	nz = images[0].get_zsize()
	ave = []
	var = []
	tota = model_blank(nx - 2 + images[0].get_attr('is_fftodd'), ny, nz)
	for k in xrange(kc):
		var.append(tota.copy())   # these are real!
	fftip(tota)
	for k in xrange(kc):
		ave.append(tota.copy())   # Fourier
	
	nclass = [0]*kc
	lctf = len(ctf2[0])
	ctf_2 = [[0.0]*lctf for k in xrange(kc)]

	# First the averages
	for im in xrange(len(images)):
		ctf_x = filt_table( images[im], ctf1[im] )
		k = assign[im]
		Util.add_img(ave[k], ctf_x)
		#Cls[k].C += ctf_x
		nclass[k] += 1
		for i in xrange(lctf): ctf_2[k][i] += ctf2[im][i]

	# and the total average
	tcft2 = [0.0]*lctf
	for k in xrange(kc):
		Util.add_img(tota, ave[k])
		for i in xrange(lctf): tcft2[i] += ctf_2[k][i]
	for i in xrange(lctf):  tcft2[i] = 1.0/(tcft2[i] + 1.0/snr)
	tota = filt_table( tota, tcft2 )

	for k in xrange(kc):
		for i in xrange(lctf):  ctf_2[k][i] = 1.0/(ctf_2[k][i] + 1.0/snr)
		ave[k] = filt_table( ave[k], ctf_2[k] )

	# Variance files have to be in real space - nobody wants to look at them in Fourier space!
	
	totv = model_blank(nx - 2 + images[0].get_attr('is_fftodd'), ny, nz) 
	for im in xrange(len(images)):
		# so we have to compute inverse FTs of all images
		fim = fft(images[im])
		#  first total variance
		temp = fft( filt_table( tota, ctf1[im] ) )
		temp = Util.subn_img( fim, temp)  # in real space now!
		Util.add_img2(totv, temp)
	
		# now variance images for each group
		k = assign[im]
		temp = fft( filt_table( ave[k], ctf1[im] ) )
		temp = Util.subn_img( fim, temp)  # in real space now!
		Util.add_img2(var[k], temp)

	Util.mul_scalar(totv, 1.0/float(len(images)-1))

	for k in xrange(kc):
		if(nclass[k] > 1):
			Util.mul_scalar(var[k], 1.0/float(nclass[k]-1))
		else:
			var[k].to_zero()

	# finally compute inverse FFT of averages
	for k in xrange(kc):
		ave[k] = fft(ave[k])

	return fft(tota), totv, ave, var, nclass
	
def aves(stack, mode="a", i1 = 0, i2 = 0):
	"""
		Calculate the average and variance for
		1. mode="a" for alignment
		2. mode=else for normal summation
	"""
	from utilities    import model_blank, get_params2D
	from fundamentals import rot_shift2D
	if (i2==0): i2 = EMUtil.get_image_count(stack)-1
	nima = i2 - i1 +1
	ima = EMData()
	ima.read_image(stack, i1)
	nx  = ima.get_xsize()
	ny  = ima.get_xsize()
	ave = model_blank(nx,ny)
	var = model_blank(nx,ny)
	for i in xrange(i1, i2 + 1):
		if (i > i1):
			ima = EMData()
			ima.read_image(stack, i)
		if (mode=="a"):
			alpha, sx, sy, mirror, scale = get_params2D(ima)
 			out = rot_shift2D(ima, alpha, sx, sy, mirror)
			Util.add_img(ave, out)
			Util.add_img2(var, out)
		else: 
			Util.add_img(ave, ima)
			Util.add_img2(var, ima)
	#var[k] = (var[k] - ave[k]*ave[k]*ii) / (ii-1)

	ave = Util.mult_scalar(ave, 1.0/float(nima))
	temp = model_blank(nx, ny)
	Util.add_img2(temp, ave)
	var = Util.madn_scalar(var, temp, -float(nima))
	Util.mul_scalar(var, 1.0/float(nima-1))

	return ave, var

def ave_ali(name_stack, name_out = None, start = 0, end = -1, ali = False):
	from statistics import aves
	from sys        import exit
	
	"""
	   Calculate the average and variance for
	   1. mode="a" for alignment
	   2. mode=else for normal summation

	   Export the result into two separate stacks

	   This function is called by sxave_ali.py
	"""
	
	# check var
	N = EMUtil.get_image_count(name_stack)
	if end == -1:
		end = N - 1
	elif end >= N:
		print 'Error: the last image is %d' % (N-1)
		exit()

	# choose the mode with or without alignment parameters
	if ali:
		mode = 'a'
	else:
		mode = 'nothing'
		
	
	# compute the ave and var
	ave, var = aves(name_stack, mode, start, end)

	if name_out is None:
		ave.write_image('ave_' + name_stack)
		var.write_image('var_' + name_stack)
	else:
		ave.write_image('ave_' + name_out)
		var.write_image('var_' + name_out)


def aves_w(stack, mode="a"):
	"""
		Apply alignment parameters, and calculate Wiener 
		Summation using CTF and SNR saved in header
		mode="a" will apply alignment parameters to the input image.
	"""
	
	from filter       import filt_table
	from utilities    import model_blank, get_params2D
	from fundamentals import rot_shift2D
	
	
	e    = EMData()
	wie  = EMData()
	nima = EMUtil.get_image_count(stack)
	e.read_image(stack,0, True)
	nx = e.get_xsize()
	ny = e.get_ysize()
	twie = model_blank(nx, ny)
	ctf_1     = e.get_attr('ctf_1')
	ctf_2_sum = [0]*len(ctf_1)
	for i in xrange(nima):
		e.read_image(stack,i)
		active = e.get_attr('active')
		if(active):
			if(mode == "a"): # for alignment 
				alpha, sx, sy, mirror, scale = get_params2D(e)
	 			out = rot_shift2D(e, alpha, sx, sy, mirror)
			else:
				out=e.copy()
			ctf_mul = e.get_attr('ctf_applied')
			ctf_1   = e.get_attr('ctf_1')
			ctf_2   = e.get_attr('ctf_2')
			TE      = e.get_attr('TE')
			snr = e.get_attr('SNR') # SNR can be either float or list
			import types
			if(type(snr) is types.ListType):
				for k in xrange(len(snr)):
					if(ctf_mul): ctf_1[k]  =snr[k]*ctf_2[k]*TE[k]*TE[k]*TE[k]
					else:        ctf_1[k] *= snr[k]*TE[k]
				Util.add_img(twie, filt_table(out, ctf_1))
				for j in xrange(len(snr)):
					if(ctf_mul==1): ctf_2_sum[j]+=ctf_2[j]*snr[j]*ctf_2[j]*TE[j]**4
					else:	       	ctf_2_sum[j]+=ctf_2[j]*snr[j]*TE[j]**2
			else:
				for k in xrange(len(TE)):
					if(ctf_mul==1): ctf_1[k]=snr*ctf_2[k]*TE[k]*TE[k]*TE[k]
					else: ctf_1[k]*=snr*TE[k]
				wie = filt_table(out, ctf_1)*snr
				Util.add_img(twie, wie)
				for j in xrange(len(TE)):
					if(ctf_mul==1): ctf_2_sum[j]+=ctf_2[j]*snr*ctf_2[j]*TE[j]**4
					else: ctf_2_sum[j]+=ctf_2[j]*snr*TE[j]**2

	for i in xrange(len(TE)):  ctf_2_sum[i] = 1./(ctf_2_sum[i]+1.)
	return filt_table(twie, ctf_2_sum)
	
def aves_wiener(input_stack, mode="a", SNR=1.0):
	"""
		Apply alignment parameters, and calculate Wiener 
		Summation using CTF info and SNR saved in header
		mode="a" will apply alignment parameters to the input image.
	"""
	from utilities import get_arb_params
	
	n = EMUtil.get_image_count(input_stack)
	ima = EMData()
	ima.read_image(input_stack, 0)
	nx = ima.get_xsize()
	ny = ima.get_xsize()
	
	if(ima.get_attr_default('ctf_applied', 2) > 0):
		ERROR("data cannot be ctf-applied","prepare_2d_forPCA",1)
	from fundamentals import fft, rot_shift2D
	from morphology   import ctf_img
	from filter 	  import filt_ctf
	from utilities    import pad, get_params2D
	parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]

	nx2 = 2*nx
	ny2 = 2*ny
	ave       = EMData(nx2, ny2, 1, False)
	ctf_2_sum = EMData(nx2, ny2, 1, False)
	from math import sqrt
	snrsqrt = sqrt(SNR)

	for i in xrange(n):
		ima = EMData()
		ima.read_image(input_stack, i)
		ctf_params = get_arb_params(ima, parnames)
		if mode == "a":
	 		alpha, sx, sy, mirror, scale = get_params2D(ima)
		 	ima = rot_shift2D(ima, alpha, sx, sy, mirror)
		oc = filt_ctf(fft(pad(ima, nx2, ny2, background = 0.0)), ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
		Util.mul_scalar(oc, SNR)
		Util.add_img(ave, oc)
		Util.add_img2(ctf_2_sum, snrsqrt*ctf_img(nx2, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5], ny = ny2, nz = 1))
	ctf_2_sum += 1.0
	Util.div_filter(ave, ctf_2_sum)
	# variance
	var = EMData(nx,ny)
	var.to_zero()
	for i in xrange(n):
		ima = EMData()
		ima.read_image(input_stack, i)
		ctf_params = get_arb_params(ima, parnames)
		if mode == "a":
			alpha, sx, sy, mirror, scale = get_params2D(ima)
 			ima = rot_shift2D(ima, alpha, sx, sy, mirror)
		oc = filt_ctf(ave, ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5], pad= True)
		Util.sub_img(ima, Util.window(fft(oc),nx,ny,1,0,0,0))
		Util.add_img2(var, ima)
	ave = Util.window(fft(ave),nx,ny,1,0,0,0)
	Util.mul_scalar(var, 1.0/(n-1))
	return ave, var

def ssnr2d(data, mask = None, mode=""):
	'''
	Calculate ssnr and variance in Fourier space for 2D or 3D images
	If mode = "a" apply alignment parameters
	'''
	from morphology   import threshold
	from utilities    import get_params2D
	from fundamentals import fft, rot_shift2D
	import  types
	if (type(data) is types.StringType):
		n = EMUtil.get_image_count(data)
		ima = EMData()
		ima.read_image(data, 0, True)
		nx = ima.get_xsize()
		ny = ima.get_ysize()
		nz = ima.get_zsize()
	else:
		n = len(data)
		nx = data[0].get_xsize()
		ny = data[0].get_ysize()
		nz = data[0].get_zsize()

	sumsq = EMData(nx, ny, nz, False)
	var   = EMData(nx, ny, nz, False)

	for i in xrange(n):
		if (type(data) is types.StringType):
			ima = EMData()
			ima.read_image(data, i)
			if(mode == "a"):
 				alpha, sx, sy, mirror, scale = get_params2D(ima)
	 			ima = rot_shift2D(ima, alpha, sx, sy, mirror)
			if(mask):  Util.mul_img(ima, mask)
			fim = fft(ima)
		else:
			if(mode == "a"):
 				alpha, sx, sy, mirror, scale = get_params2D(data[i])
	 			ima = rot_shift2D(data[i], alpha, sx, sy, mirror)
				if(mask):  fim = fft(Util.muln_img(ima, mask))
				else    :  fim = fft(ima)
			else:
				if(mask):  fim = fft(Util.muln_img(data[i], mask))
				else:      fim = fft(data[i])
		Util.add_img(sumsq, fim)
		Util.add_img2(var, fim)
	Util.mul_img(sumsq, sumsq.conjg())
	var   = Util.pack_complex_to_real(var)
	sumsq = Util.pack_complex_to_real(sumsq)
	var = (var - sumsq/n)/(n-1)
	# convert to real images
	ssnr   = sumsq/var/n - 1.0
	from fundamentals import rot_avg_table
	rvar = rot_avg_table(var)
	rsumsq = rot_avg_table(sumsq)
	rssnr = []
	for i in xrange(len(rvar)):
		if(rvar[i] > 0.0): qt = max(0.0, rsumsq[i]/rvar[i]/n - 1.0)
		else:              ERROR("ssnr2d","rvar negative",1)
		rssnr.append(qt)

	return rssnr, rsumsq, rvar, ssnr, sumsq, var

def ssnr2d_ctf(data, mask = None, mode=""):
	'''
	Calculate ssnr and variance in Fourier space for 2D images including CTF information
	If mode = "a" apply alignment parameters
	'''
	from fundamentals import fft, rot_shift2D, rot_avg_table
	from morphology   import ctf_img, threshold
	from filter       import filt_ctf
	from utilities    import get_arb_params, get_params2D
	import  types
	
	parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
	if (type(data) is types.StringType):
		n = EMUtil.get_image_count(data)
		ima = EMData()
		ima.read_image(data, 0, True)
		if(ima.get_attr_default('ctf_applied', 1) == 1):
			ERROR("data cannot be ctf-applied","varfctf",1)
		nx = ima.get_xsize()
		ny = ima.get_ysize()
		nz = ima.get_zsize()
	else:
		if(data[0].get_attr_default('ctf_applied', 1) == 1):
			ERROR("data cannot be ctf-applied","varfctf",1)
		n = len(data)
		nx = data[0].get_xsize()
		ny = data[0].get_ysize()
		nz = data[0].get_zsize()

	ctf_2_sum = EMData(nx, ny, nz, False)
	sumsq     = EMData(nx, ny, nz, False)
	var       = EMData(nx, ny, nz, False)

	for i in xrange(n):
		if (type(data) is types.StringType):
			ima = EMData()
			ima.read_image(data, i)
		else:
			ima = data[i].copy()
		ctf_params = get_arb_params(ima, parnames)
		if(mode == "a"):
			alpha, sx, sy, mirror, scale = get_params2D(ima)
 			ima = rot_shift2D(ima, alpha, sx, sy, mirror)
		if(mask):  Util.mul_img(ima, mask)
		ctfimg = ctf_img(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5], ny = ny, nz = nz)
		Util.add_img2(ctf_2_sum, ctfimg)
		ima = fft(ima)
		Util.add_img2(var, ima)
		Util.mul_img(ima, ctfimg)
		Util.add_img(sumsq, ima)
	Util.div_filter(sumsq, ctf_2_sum)
	Util.mul_img(sumsq, sumsq.conjg())
	Util.mul_img(sumsq, ctf_2_sum)
	Util.sub_img(var, sumsq)
	Util.mul_scalar(var, 1.0/float(n-1))
	
	var   = Util.pack_complex_to_real(var)
	sumsq = Util.pack_complex_to_real(sumsq)
	ssnr   = sumsq/var - 1.0
	rvar = rot_avg_table(var)
	rsumsq = rot_avg_table(sumsq)
	rssnr = []
	for i in xrange(len(rvar)):
		if(rvar[i] > 0.0): qt = max(0.0, rsumsq[i]/rvar[i] - 1.0)
		else:              ERROR("ssnr2d","rvar negative",1)
		rssnr.append(qt)
	return rssnr, rsumsq, rvar, ssnr, sumsq, var

def varf(data, mask = None, mode=""):
	'''
	Calculate variance in Fourier space for 2D or 3D images
	If mode = "a" apply alignment parameters
	'''
	from fundamentals import fft, rot_shift2D
	from utilities    import get_params2D
	import  types
	if (type(data) is types.StringType):
		n = EMUtil.get_image_count(data)
		ima = EMData()
		ima.read_image(data, 0, True)
		nx = ima.get_xsize()
		ny = ima.get_ysize()
		nz = ima.get_zsize()
	else:
		n = len(data)
		nx = data[0].get_xsize()
		ny = data[0].get_ysize()
		nz = data[0].get_zsize()

	sumsq = EMData(nx, ny, nz, False)
	var   = EMData(nx, ny, nz, False)

	for i in xrange(n):
		if (type(data) is types.StringType):
			ima = EMData()
			ima.read_image(data, i)
		else:
			ima = data[i].copy()
		if(mode == "a"):
		 	alpha, sx, sy, mirror, scale = get_params2D(ima)
	 		ima = rot_shift2D(ima, alpha, sx, sy, mirror)
		if(mask):  Util.mul_img(ima, mask)
		fim = fft(ima)
		Util.add_img(sumsq, fim)
		Util.add_img2(var, fim)

	Util.mul_img(sumsq, sumsq.conjg())
	Util.mad_scalar(var, sumsq, -1.0/float(n))
	Util.mul_scalar(var, 1.0/float(n-1))
	st = Util.infomask(var, None, True)
	if(st[2]<0.0):  ERROR("Negative variance!","varf",1)
	from fundamentals import rot_avg_table

	return var, rot_avg_table(Util.pack_complex_to_real(var))

def varfctf(data, mask = None, mode=""):
	'''
	Calculate variance in Fourier space for 2D or 3D images including ctf information
	If mode = "a" apply alignment parameters
	'''
	from fundamentals import fft, rot_shift2D
	from morphology   import ctf_img
	from filter       import filt_ctf
	from utilities    import get_arb_params, get_params2D
	import  types

	parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
	if (type(data) is types.StringType):
		n = EMUtil.get_image_count(data)
		ima = EMData()
		ima.read_image(data, 0, True)
		if(ima.get_attr_default('ctf_applied', 1) == 1):
			ERROR("data cannot be ctf-applied","varfctf",1)
		nx = ima.get_xsize()
		ny = ima.get_ysize()
		nz = ima.get_zsize()
	else:
		if(data[0].get_attr_default('ctf_applied', 1) == 1):
			ERROR("data cannot be ctf-applied","varfctf",1)
		n = len(data)
		nx = data[0].get_xsize()
		ny = data[0].get_ysize()
		nz = data[0].get_zsize()

	ctf_2_sum = EMData(nx, ny, nz, False)
	sumsq     = EMData(nx, ny, nz, False)
	var       = EMData(nx, ny, nz, False)

	for i in xrange(n):
		if (type(data) is types.StringType):
			ima = EMData()
			ima.read_image(data, i)
		else:
			ima = data[i].copy()
		ctf_params = get_arb_params(ima, parnames)
		if(mode == "a"):
			alpha, sx, sy, mirror, scale = get_params2D(ima)
 			ima = rot_shift2D(ima, alpha, sx, sy, mirror)
		if(mask): Util.mul_img(ima, mask)
		oc = filt_ctf(ima, ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5], pad=True)
		Util.add_img(sumsq, fft(oc))
		Util.add_img2(var, fft(ima))
		Util.add_img2(ctf_2_sum, ctf_img(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5], ny = ny, nz = nz))
	Util.mul_img(sumsq, sumsq.conjg())
	Util.div_img(sumsq, ctf_2_sum)
	Util.sub_img(var, sumsq)
	Util.mul_scalar(var, 1.0/float(n-1))
	st = Util.infomask(var, None, True)
	if(st[2]<0.0):  ERROR("Negative variance!","varfctf",1)

	from fundamentals import rot_avg_table

	return var, rot_avg_table(Util.pack_complex_to_real(var))


def ccc(img1, img2, mask=None):
	"""Cross-correlation coefficient.
	   
	   Usage: result = ccc(image1, image2 [, mask])
	"""
	if mask: return img1.cmp("ccc", img2, {"mask":mask,"negative":0})
	else:    return img1.cmp("ccc", img2, {"negative":0})

def fsc(img1, img2, w = 1.0, filename=None):
	"""Fourier Shell (or Ring) Correlation.

	   Usage: [frsc =] fsc(image1, image2 [, w, filename])

	   Computes the Fourier Shell (3d) or Fourier Ring (2d) correlation
	   function of two images.  If a filename is provided, then the 
	   result is saved using that filename.
	"""
	result = img1.calc_fourier_shell_correlation(img2, w)
	# repack results as a list of (freq, fsc, n) triplets
	size = len(result)/3
	frsc = []
	for i in xrange(3):
		frsc.append(result[i*size:(i+1)*size])
	if filename:
		outf = file(filename, "w")
		for i in xrange(size):
			datstrings = []
			datstrings.append("  %12f" % (frsc[0][i]))
			datstrings.append("  %12f" % (frsc[1][i]))
			datstrings.append("  %12f" % (frsc[2][i]))
			datstrings.append("\n")
			outf.write("".join(datstrings))
		outf.close()
	return frsc

def fsc_mask(img1, img2, mask = None, w = 1.0, filename=None):
	""" 
	        Compute fsc using mask and normalization.

		Usage: [frsc =] fsc(image1, image2 [, mask, w, filename])

		If no mask is provided, using circular mask with R=nx//2

	"""
	from statistics import fsc
	from utilities import model_circle
	from morphology import binarize, erosion
	nx = img1.get_xsize()
	ny = img1.get_ysize()
	nz = img1.get_zsize()
	if( mask == None):  mask = model_circle(nx//2, nx, ny, nz)
	m = binarize(mask, 0.5)
	Util.sub_img(m, erosion(m))
	s1 = Util.infomask(img1, m, True)
	s2 = Util.infomask(img2, m, True)
	return fsc((img1-s1[0])*mask, (img2-s2[0])*mask, w, filename)

def get_refstack(imgstack,params,nref,refstack,cs,mask,center,Iter):

	"""
		Calculate multiple references from imgstack using aligment parameter table
	"""
	from filter import fshift
			 	
 	refimg=EMData()
	refimgo=EMData()
	refimge=EMData()
 	nima = EMUtil.get_image_count(imgstack)
 	ima = EMData()
	tc=[]
	print len(params)
	for ir in xrange(nref):
		ncnt=0
		for im in xrange(nima):
			if(ir+1==int(params[im][4])):
				ncnt+=1
				ima.read_image(imgstack,im)
				out = rot_shift2D(ima, params[im][0], params[im][1], params[im][2], params[im][3])
				if(ncnt<=2):
					if(ncnt%2==0): refimge=out
					if(ncnt%2==1): refimgo=out
				else:
					if(ncnt%2==0): refimge+=out
					if(ncnt%2==1): refimgo+=out
		if(center):
			if(ncnt>=2):
				tavg= (refimgo*(int(ncnt/2)+(ncnt%2)) + refimge*int(ncnt/2))/ncnt
				dropImage(tavg,"tavg.spi")
				cs[ir] = tavg.phase_cog()
				refimg = fshift(tavg, -cs[ir][0], -cs[ir][1])
			else:
				cs[ir] = refimgo.phase_cog()
				refimg = fshift(refimgo, -cs[ir][0], -cs[ir][1])				
		else:
			if(ncnt>=2):
				refimg= (refimgo*(int(ncnt/2)+(ncnt%2)) + refimge*int(ncnt/2))/ncnt
			else:
				refimg=refimgo
		a0 = refimg.cmp("dot", refimg, {"negative":0, "mask":mask}) # tave' * tave
		print " ITERATION #",'%5d'%(Iter),"  criterion = ",'%11.4g'%(a0), "reference number", '%5d'%(Iter)
		refimg.write_image(refstack,ir)
	return cs
	
def get_1dpw_table_stack(stack):
	"""
		calculate the 1D rotationally averaged power spectrum of one stack file
	"""
	table=[]
	nima  = EMUtil.get_image_count(stack)
	img   = EMData()
	rosum = EMData()
	for i in xrange(nima):
		img.read_image(stack,i)
		e  = periodogram(img)
		ro = e.rotavg()
		if(i==0): rosum = ro.copy()
		else: rosum += ro
	rosum=rosum/nima
	nr = ro.get_xsize()
	for ir in xrange(nr):  table.append([ir/float(nr-1)/2, rosum.get_value_at(ir)])
	return table

def histogram(image, mask = None, nbins = 0, hmin = 0.0, hmax = 0.0):
	return  Util.histogram(image, mask, nbins, hmin, hmax)

def im_diff(im1, im2, mask = None):
	import types
	from utilities import model_circle, get_im
	if type(im1) == types.StringType : im1 = get_im(im1)
	if type(im2) == types.StringType : im2 = get_im(im2) 
	if type(mask) == types.ClassType : l = Util.im_diff(im1, im2, mask)
	else:   
		nx = im1.get_xsize()
		ny = im1.get_ysize()
		nz = im1.get_zsize()
		if   type(mask) == types.FloatType or type(mask) == types.IntType: m = model_circle(mask, nx, ny, nz)
		elif type(mask) == types.StringType:   m = get_im(mask)
		else:
			if   im1.get_ndim() == 3: radius = min(nx,ny,nz)//2 - 1
			elif im1.get_ndim() == 2: radius = min(nx,ny)//2    - 1
			else:                     radius = int(nx)//2       - 1
			m = model_circle(radius, nx, ny, nz)
		l = Util.im_diff(im1, im2, m)
	return  l["imdiff"], l["A"], l["B"]


##############################################################################################
### K-MEANS ##################################################################################

# k-means open and prepare images
def k_means_open_im(stack, maskname, N_start, N_stop, N, CTF):
	from utilities     import get_params2D, getImage
	from fundamentals  import rot_shift2D, rot_shift3D
	
	if CTF:
		from morphology		import ctf_2, ctf_1d
		from filter		import filt_ctf, filt_table
		from fundamentals 	import fftip

	im_M = [0] * N
	im = EMData()
	im.read_image(stack, N_start, True)
	nx = im.get_xsize()
	ny = im.get_ysize()
	nz = im.get_zsize()
	
	if CTF:
		parnames    = ('Pixel_size', 'defocus', 'voltage', 'Cs', 'amp_contrast', 'B_factor',  'ctf_applied')
		ctf	    = []
		ctf2        = []
		ctf_params  = get_arb_params(im, parnames)
		if ctf_params[6]: ERROR('K-means cannot be performed on CTF-applied images', 'k_means', 1)

	if maskname != None:
		if isinstance(maskname, basestring):
			mask = getImage(maskname)
	else:
		mask = None

	DATA = im.read_images(stack, range(N_start, N_stop))
	ct   = 0
	for i in xrange(N_start, N_stop):
		image = DATA[ct].copy()
		# 3D object
		if nz > 1:
			try:	phi, theta, psi, s3x, s3y, s3z, mirror = get_params3D(image)
			except:	phi, theta, psi, s3x, s3y, s3z, mirror = 0, 0, 0, 0, 0, 0, 0
			image  = rot_shift3D(image, phi, theta, psi, s3x, s3y, s3z, scale)
			if mirror: image.process_inplace('mirror', {'axis':'x'})
		# 2D object
		elif ny > 1:
			try:	alpha, sx, sy, mirror = get_params2D(image)
			except: alpha, sx, sy, mirror  = 0, 0, 0, 0
			image = rot_shift2D(image, alpha, sx, sy, mirror)
		# obtain ctf
		if CTF:
			ctf_params = get_arb_params(image, parnames)
			ctf.append(ctf_1d(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5]))
			ctf2.append(ctf_2(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5]))

		# apply mask
		if mask != None:
			if CTF: Util.mul_img(image, mask)
			else: image = Util.compress_image_mask(image, mask)

		# fft
		if CTF: fftip(image)

		# mem the original size
		if im == N_start:
			image.set_attr('or_nx', nx)
			image.set_attr('or_ny', ny)
			image.set_attr('or_nz', nz)

		# store image
		im_M[i] = image.copy()
		ct += 1

	if CTF: return im_M, mask, ctf, ctf2
	else:   return im_M, mask, None, None

# k-means init for MPI version
def k_means_init_MPI(stack):
	from mpi 	  import mpi_init, mpi_comm_size, mpi_comm_rank, mpi_barrier, MPI_COMM_WORLD
	from mpi          import mpi_bcast, MPI_INT
	import sys
	
	# init
	sys.argv       = mpi_init(len(sys.argv),sys.argv)
	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)

	# chose a main one
	main_node = 0

	# define the node affected to the images
	N = 0
	if myid == main_node: N = EMUtil.get_image_count(stack)
	mpi_barrier(MPI_COMM_WORLD)
	N = mpi_bcast(N, 1, MPI_INT, main_node, MPI_COMM_WORLD)
	N = N.tolist()[0]

	im_per_node = max(N / number_of_proc, 1)
	N_start     = myid * im_per_node
	if myid == (number_of_proc -1):
		N_stop = N
	else:
		N_stop = min(N_start + im_per_node, N)

	return main_node, myid, number_of_proc, N_start, N_stop, N

# k-means write the head of the logfile
def k_means_headlog(stackname, outname, method, N, K, crit, maskname, trials, maxit, CTF, T0, F, SA2, rnd, ncpu):
	from utilities import print_msg

	if F != 0: SA = True
	else:      SA = False

	if isinstance(K, list):
		if len(K) == 2: Ks = 'from ' + str(K[0]) + ' to ' + str(K[1])
	else:
		Ks = str(K)

	if method == 'cla': method = 'Classical'

	if ncpu > 1: methodhead = method + ' MPI'
	else:        methodhead = method

	print_msg('\n************* k-means %s *************\n' % methodhead)
	print_msg('Input stack                 : %s\n'     % stackname)
	print_msg('Number of images            : %i\n'     % N)
	print_msg('Maskfile                    : %s\n'     % maskname)
	print_msg('Number of clusters          : %s\n'     % Ks)
	print_msg('Number of trials            : %i\n'     % trials)
	print_msg('Maximum iteration           : %i\n'     % maxit)
	print_msg('Data with CTF               : %s\n'     % CTF)
	print_msg('Criterion                   : %s\n'     % crit)
	print_msg('Optimization method         : %s\n'     % method)
	if SA:
		print_msg('Simulate annealing          : ON\n')
		if SA2: print_msg('   select neighbour         : closer according T\n')
		else:	print_msg('   select neighbour         : randomly\n')
		print_msg('   T0                       : %f\n' % T0)
		print_msg('   F                        : %f\n' % F)
	else:
		print_msg('Simulate annealing          : OFF\n')
	print_msg('Random seed                 : %i\n'     % rnd)
	print_msg('Number of cpus              : %i\n'     % ncpu)
	print_msg('Output seed names           : %s\n\n'   % outname)

# K-means write results output directory
def k_means_export(Cls, crit, assign, out_seedname):
	from utilities import print_msg
	
	if out_seedname.split(':')[0] == 'bdb':
		BDB = True
	else:
		BDB = False
		import os
		if os.path.exists(out_seedname):  os.system('rm -rf ' + out_seedname)
		os.mkdir(out_seedname)

	# write the report on the logfile
	Je = 0
	for k in xrange(Cls['k']): Je += Cls['Ji'][k]

	print_msg('\n\n_Details____________________________________________________\n')
	print_msg('\n\t%s\t%11.6e\n\n' % ('The total Sum of Squares Error (Je) = ', Je))

	for name in crit['name']:
		if name == 'C':   print_msg('\t%s\t%11.4e\n' % ('Criteria Coleman', crit['C']))
		elif name == 'H': print_msg('\t%s\t%11.4e\n' % ('Criteria Harabasz', crit['H']))
		elif name == 'D': print_msg('\t%s\t%11.4e\n' % ('Criteria Davies-Bouldin', crit['D']))
		else:             ERROR('Kind of criterion k-means unknown', 'k_means_out_res', 0)	
	print_msg('\n')

	for k in xrange(Cls['k']):
		print_msg('\t%s\t%d\t%s\t%d' % ('Cluster no:', k, 'No of Objects = ', Cls['n'][k]))
		if(Cls['n'][k] > 1): print_msg('\t%s\t%11.6e\t%s\t%11.6e\n' % ('Sum of Squares Error Ji', Cls['Ji'][k], ' Variance', Cls['Ji'][k] / float(Cls['n'][k]-1)))
		else:                print_msg('\t%s\t%11.6e\n' % ('Sum of Squares Error Ji', Cls['Ji'][k]))

		# limitation of hdf file in the numbers of attributes
		if Cls['n'][k] > 16000 and not BDB:
			print 'WARNING: limitation of number attributes in hdf file, the results will be export in separate files \n'
			outfile = open('%d_kmeans' % k, 'w')
			list_images = []
			for i in xrange(len(assign)):
				if assign[i] == k:
					list_images.append(i)
					outfile.write(str(i) +'\n')
			outfile.close()
			Cls['ave'][k].set_attr_dict({'Class_average':1.0, 'nobjects':float(Cls['n'][k])})
		else:
			lassign = []
			for i in xrange(len(assign)):
				if(assign[i] == k):  lassign.append(float(i))
			Cls['ave'][k].set_attr('Class_average', 1.0)
			Cls['ave'][k].set_attr('nobjects', float(Cls['n'][k]))
			Cls['ave'][k].set_attr('members', lassign)
			Cls['ave'][k].set_attr('Ji', Cls['Ji'][k])
			Cls['ave'][k].set_attr('Je', Je)

			Cls['var'][k].set_attr('Class_variance', 1.0)
			Cls['var'][k].set_attr('nobjects', float(Cls['n'][k]))
			Cls['var'][k].set_attr('members', lassign)
			Cls['var'][k].set_attr('Ji', Cls['Ji'][k])
			Cls['var'][k].set_attr('Je', Je)

		if BDB:
			Cls['ave'][k].write_image(out_seedname + '_ave', k)
			Cls['var'][k].write_image(out_seedname + '_var', k)
		else:
			Cls['ave'][k].write_image(out_seedname + "/average.hdf", k)
			Cls['var'][k].write_image(out_seedname + "/variance.hdf", k)


# K-means compute criterion in order to validate the number of groups
def k_means_criterion(Cls, crit_name=''):
	from utilities		import model_blank
	from fundamentals	import fftip
	
	if crit_name == 'all':	crit_name = 'CHD'
	
	# Informations about images
	nx   = Cls['ave'][0].get_xsize()
	ny   = Cls['ave'][0].get_ysize()
	nz   = Cls['ave'][0].get_zsize()
	N    = Cls['N']
	
	# if complex need buf complex
	if Cls['ave'][0].is_complex():
		buf  = model_blank(nx, ny, nz)
		buf.set_complex(1)
	else:	
		buf  = model_blank(nx, ny, nz)
				
	# Criterion
	Crit         = {}
	Crit['name'] = crit_name
	Crit['C']    = 0
	Crit['H']    = 0
	Crit['D']    = 0
			
	# Compute the criterion required
	for name in Crit['name']:
		
		# Coleman criterion C = Tr(Cls) * S Tr(Cls(im)) -> C = S (ave-m_ave)**2  *  Je
		if name == 'C':
			Tr, Je = 0, 0
			buf.to_zero()
			for k in xrange(Cls['k']):
				Util.add_img(buf, Cls['ave'][k])
				Je += Cls['Ji'][k]
			Util.mul_scalar(buf, 1.0/float(Cls['k']))
			for k in xrange(Cls['k']):	Tr += buf.cmp("SqEuclidean", Cls['ave'][k])
			Crit['C'] = Tr * Je
		
		# Harabasz criterion H = [Tr(Cls)/(K-1)] / [S Tr(Cls(im))/(N-K)]
		elif name == 'H':
			Tr, Je = 0, 0
			
			if Crit['C'] == 0:			
				buf.to_zero()
				for k in xrange(Cls['k']):
					Util.add_img(buf, Cls['ave'][k])
					Je += Cls['Ji'][k]
				Util.mul_scalar(buf, 1.0/float(Cls['k']))
				for k in xrange(Cls['k']):	Tr += buf.cmp("SqEuclidean", Cls['ave'][k])
				if Je > 0:
					Crit['H'] = (Tr / (Cls['k'] - 1)) / (Je / (N - Cls['k']))
				else:
					Crit['H'] = 1e20
			else:
				# Tr already compute in Crit['C']
				for k in xrange(Cls['k']):	Je += Cls['Ji'][k]							
				if Je > 0:
					Crit['H'] = ((Crit['C'] / Je) / (Cls['k'] - 1)) / (Je / (N - Cls['k']))
				else:
					Crit['H'] = 1e20				
		
		# Davies-Bouldin criterion DB = 1/K S_i[max all(i>j) ((Ji/n) + (Jj/n))/(ave_j-ave_i)**2]
		elif name == 'D':
			db, DB, err, ji, jj = 0, 0, 0, 0, 0
			for i in xrange(Cls['k']):
				val_max = 0
				for j in xrange(Cls['k']):
					if i != j:		
						err = Cls['ave'][j].cmp("SqEuclidean",Cls['ave'][i])
						if err > 0:
							if Cls['n'][i] > 0:	ji = Cls['Ji'][i] / Cls['n'][i]
							else:			ji = 0
							if Cls['n'][j] > 0:	jj = Cls['Ji'][j] / Cls['n'][j]
							else:			jj = 0
							db = (ji + jj) / err
						else:
							db = 1e20													
					if db > val_max:	val_max = db					
				DB += val_max				
			Crit['D'] = DB / Cls['k']			
			
		else:
			ERROR("Kind of criterion k-means unknown","k_means_criterion",1)
	
	# return the results
	return Crit

# K-means with classical method
def k_means_classical(im_M, mask, K, rand_seed, maxit, trials, CTF, F=0, T0=0, SA2=False, DEBUG=False):
	from utilities 		import model_blank, get_im, running_time
	from random    		import seed, randint
	from utilities 		import print_msg
	from copy		import deepcopy
	from sys		import exit
	import time
	if CTF[0]:
		from filter	        import filt_ctf, filt_table
		from fundamentals 	import fftip

		ctf  = deepcopy(CTF[1])
		ctf2 = deepcopy(CTF[2])
		CTF  = True
	else:
		CTF  = False 

	# Simulate annealing use or not
	if F != 0: SA = True
	else:      SA = False

	if SA:
		# for simulate annealing
		from math   import exp
		from random import random

	if mask != None:
		if isinstance(mask, basestring):
			ERROR('Mask must be an image, not a file name!', 'k-means', 1)

	N = len(im_M)

	t_start = time.time()
		
	# Informations about images
	if CTF:
		nx  = im_M[0].get_attr('or_nx')
		ny  = im_M[0].get_attr('or_ny')
		nz  = im_M[0].get_attr('or_nz')
		buf = model_blank(nx, ny, nz)
		fftip(buf)		
		nx   = im_M[0].get_xsize()
		ny   = im_M[0].get_ysize()
		nz   = im_M[0].get_zsize()
		norm = nx * ny * nz
	else:
		nx   = im_M[0].get_xsize()
		ny   = im_M[0].get_ysize()
		nz   = im_M[0].get_zsize()
		norm = nx * ny * nz
		buf  = model_blank(nx, ny, nz)

	# Variables			
	if rand_seed > 0:  seed(rand_seed)
	else:              seed()
	Cls        = {}
	Cls['n']   = [0]*K   # number of objects in a given cluster
	Cls['ave'] = [0]*K   # value of cluster average
	Cls['var'] = [0]*K   # value of cluster variance
	Cls['Ji']  = [0]*K   # value of Ji
	Cls['k']   =  K	     # value of number of clusters
	Cls['N']   =  N
	assign     = [0]*N 
	
	if CTF:
		Cls_ctf2    = {}
		len_ctm	    = len(ctf2[0])
			
	# TRIALS
	if trials > 1:
		MemCls, MemJe, MemAssign = {}, {}, {}
	else:
		trials = 1
	ntrials = 0
	while ntrials < trials:
		ntrials  += 1

		# for simulate annealing
		if SA: T = T0
		
		# Init the cluster by an image empty
		buf.to_zero()
		for k in xrange(K):
			Cls['ave'][k] = buf.copy()
			Cls['var'][k] = buf.copy()
			Cls['n'][k]   = 0
			Cls['Ji'][k]  = 0
		
		## Random method
		retrial = 20
		while retrial > 0:
			retrial -= 1
			i = 0
			for im in xrange(N):
				assign[im] = randint(0, K-1)
				Cls['n'][assign[im]] += 1
			
			flag, k = 1, K
			while k>0 and flag:
				k -= 1
				if Cls['n'][k] == 0:
					flag = 0
					if retrial == 0:
						ERROR('Empty class in the initialization', 'k_means_classical', 1)
					for k in xrange(K):
						Cls['n'][k] = 0
			
			if flag == 1:	retrial = 0
		
		## Calculate averages, if CTF: ave = S CTF.F / S CTF**2
		if CTF:
			# first init ctf2
			for k in xrange(K):	Cls_ctf2[k] = [0] * len_ctm
						
			for im in xrange(N):
				# compute ctf2				
				for i in xrange(len_ctm):	Cls_ctf2[assign[im]][i] += ctf2[im][i]
				
				# compute average first step
				CTFxF = filt_table(im_M[im], ctf[im])
				Util.add_img(Cls['ave'][assign[im]], CTFxF)
						
			for k in xrange(K):
				for i in xrange(len_ctm):	Cls_ctf2[k][i] = 1.0 / float(Cls_ctf2[k][i])
				Cls['ave'][k] = filt_table(Cls['ave'][k], Cls_ctf2[k])

			# compute Ji and Je
			for n in xrange(N):
				CTFxAve               = filt_table(Cls['ave'][assign[n]], ctf[n])
				Cls['Ji'][assign[n]] += CTFxAve.cmp("SqEuclidean", im_M[n]) / norm
			Je = 0
			for k in xrange(K):        Je = Cls['Ji'][k]
																			
		else:
			# compute average
			for im in xrange(N):	Util.add_img(Cls['ave'][assign[im]], im_M[im])
			for k in xrange(K):	Cls['ave'][k] = Util.mult_scalar(Cls['ave'][k], 1.0 / float(Cls['n'][k]))

			# compute Ji and Je
			Je = 0
			for n in xrange(N):	Cls['Ji'][assign[n]] += im_M[n].cmp("SqEuclidean",Cls['ave'][assign[n]])/norm
			for k in xrange(K):	Je += Cls['Ji'][k]	
		
		## Clustering		
		ite       = 0
		watch_dog = 0
		old_Je    = 0
		change    = True

		if DEBUG: print 'init Je', Je
		
		print_msg('\n__ Trials: %2d _________________________________%s\n'%(ntrials, time.strftime('%a_%d_%b_%Y_%H_%M_%S', time.localtime())))
		while change and watch_dog < maxit:
			ite       += 1
			watch_dog += 1
			change     = False
			Je	   = 0
			if SA: ct_pert = 0
		
			for im in xrange(N):
		
				if CTF:
					CTFxAVE = []
					for k in xrange(K): CTFxAVE.append(filt_table(Cls['ave'][k], ctf[im]))
					res = Util.min_dist(im_M[im], CTFxAVE)
				else:
					res = Util.min_dist(im_M[im], Cls['ave'])
				
				# Simulate annealing
				if SA:
					if SA2:
						dJe = [0.0] * K
						ni  = float(Cls['n'][assign[im]])
						di  = res['dist'][assign[im]]
											
						for k in xrange(K):
							if k != assign[im]:
								nj  = float(Cls['n'][k])
								dj  = res['dist'][k]
								
								dJe[k] = -( (nj/(nj+1))*(dj/norm) - (ni/(ni-1))*(di/norm) )
															
							else:
								dJe[k] = 0

						# norm <0 [-1;0], >=0 [0;+1], if just 0 norm to 1
						nbneg  =  0
						nbpos  =  0
						minneg =  0
						maxpos =  0
						for k in xrange(K):
							if dJe[k] < 0.0:
								nbneg += 1
								if dJe[k] < minneg: minneg = dJe[k]
							else:
								nbpos += 1
								if dJe[k] > maxpos: maxpos = dJe[k]
						if nbneg != 0:                   dneg = -1.0 / minneg
						if nbpos != 0 and maxpos != 0:   dpos =  1.0 / maxpos
						for k in xrange(K):
							if dJe[k] < 0.0: dJe[k] = dJe[k] * dneg
							else:
								if maxpos != 0: dJe[k] = dJe[k] * dpos
								else:           dJe[k] = 1.0

						# q[k]
						q      = [0.0] * K
						arg    = [0.0] * K
						maxarg = 0
						for k in xrange(K):
							arg[k] = dJe[k] / T
							if arg[k] > maxarg: maxarg = arg[k]
						limarg = 17
						if maxarg > limarg:
							sumarg = float(sum(arg))
							for k in xrange(K): q[k] = exp(arg[k] * limarg / sumarg)
						else:
							for k in xrange(K): q[k] = exp(arg[k])
										
						# p[k]
						p = [[0.0, 0] for i in xrange(K)]
						sumq = float(sum(q))
						for k in xrange(K):
							p[k][0] = q[k] / sumq
							p[k][1] = k
											
						p.sort()
						c = [0.0] * K
						c[0] = p[0][0]
						for k in xrange(1, K): c[k] = c[k-1] + p[k][0]

						pb = random()
						select = -1
						for k in xrange(K):
							if c[k] > pb:
								select = p[k][1]
								break


						if select != res['pos']:
							ct_pert    += 1
							res['pos']  = select


					else:
						if exp( -(1) / float(T) ) > random():
							res['pos']  = randint(0, K - 1)
							ct_pert    += 1

				# update assign
				if res['pos'] != assign[im]:
					Cls['n'][assign[im]] -= 1
					if Cls['n'][assign[im]] < 1: Cls['n'][assign[im]] = 0
					assign[im]            = res['pos']
					Cls['n'][assign[im]] += 1
					change                = True

							
			# manage empty cluster
			for k in xrange(K):
				if Cls['n'][k] < 1: ERROR('Empty groups in kmeans_classical', 'k_means_classical', 1)
													
			# Update clusters
			for k in xrange(K):
				Cls['ave'][k] = buf.copy()
				Cls['Ji'][k]  = 0
				Je = 0
			
			if CTF:
				# first init ctf2
				for k in xrange(K):	Cls_ctf2[k] = [0] * len_ctm

				for im in xrange(N):
					# compute ctf2				
					for i in xrange(len_ctm):	Cls_ctf2[assign[im]][i] += ctf2[im][i]

					# compute average first step
					CTFxF = filt_table(im_M[im], ctf[im])
					Util.add_img(Cls['ave'][assign[im]], CTFxF)

				for k in xrange(K):
					for i in xrange(len_ctm):	Cls_ctf2[k][i] = 1.0 / float(Cls_ctf2[k][i])
					Cls['ave'][k] = filt_table(Cls['ave'][k], Cls_ctf2[k])

				# compute Ji and Je
				for n in xrange(N):
					CTFxAve               = filt_table(Cls['ave'][assign[n]], ctf[n])
					Cls['Ji'][assign[n]] += CTFxAve.cmp("SqEuclidean", im_M[n]) / norm
				Je = 0
				for k in xrange(K):       Je += Cls['Ji'][k]
			
			else:
				for im in xrange(N): Util.add_img(Cls['ave'][assign[im]], im_M[im])
				for k in xrange(K):  Cls['ave'][k] = Util.mult_scalar(Cls['ave'][k], 1.0/float(Cls['n'][k]))

				# compute Ji and Je
				Je = 0
				for n in xrange(N):	Cls['Ji'][assign[n]] += im_M[n].cmp("SqEuclidean",Cls['ave'][assign[n]]) / norm
				for k in xrange(K):	Je += Cls['Ji'][k]	
									
				
			# threshold convergence control
			if Je != 0:
				thd = abs(Je - old_Je) / Je
			else:
				thd = 0
			if SA:
				if thd < 1.0e-12 and ct_pert == 0: watch_dog =maxit
			else:
				if thd < 1.0e-8:                   watch_dog = maxit
			old_Je = Je

			if SA:
				# Simulate annealing, update temperature
				T *= F
				if SA2:
					if T < 0.09: SA = False

				print_msg('> iteration: %5d    criterion: %11.6e    T: %13.8f  ct disturb: %5d\n' % (ite, Je, T, ct_pert))
				if DEBUG: print '> iteration: %5d    criterion: %11.6e    T: %13.8f  ct disturb: %5d' % (ite, Je, T, ct_pert)
			else:
				print_msg('> iteration: %5d    criterion: %11.6e\n' % (ite, Je))
				if DEBUG: print '> iteration: %5d    criterion: %11.6e' % (ite, Je)
						
		# memorize the result for this trial	
		if trials > 1:
			MemCls[ntrials-1]    = deepcopy(Cls)
			MemJe[ntrials-1]     = deepcopy(Je)
			MemAssign[ntrials-1] = deepcopy(assign)
						
	# if severals trials choose the best
	if trials > 1:
		val_min = 1.0e20
		best    = -1
		for n in xrange(trials):
			if MemJe[n] < val_min:
				val_min = MemJe[n]
				best    = n
		
		# affect the best
		Cls    = MemCls[best]
		Je     = MemJe[best]
		assign = MemAssign[best]		
	
	if CTF:
		# compute Ji and the variance S (F - CTF * Ave)**2
		for n in xrange(N):
			CTFxAve   	      = filt_table(Cls['ave'][assign[n]], ctf[n])	
			Cls['Ji'][assign[n]] += CTFxAve.cmp("SqEuclidean", im_M[n]) / norm
			
			buf.to_zero()
			buf = Util.subn_img(CTFxAve, im_M[n])
			Util.add_img(Cls['var'][assign[n]], buf) # **2
			
	else:
		# compute Ji
		for n in xrange(N): 	Cls['Ji'][assign[n]] += im_M[n].cmp("SqEuclidean",Cls['ave'][assign[n]]) / norm
		
		# compute the variance 1/n S(im-ave)**2 -> 1/n (Sim**2 - n ave**2)
		for im in xrange(N):	Util.add_img2(Cls['var'][assign[im]],im_M[im])
		for k in xrange(K):
			buf.to_zero()
			Util.add_img2(buf, Cls['ave'][k])
			Util.mad_scalar(Cls['var'][k], buf, -float(Cls['n'][k]))
			Util.mul_scalar(Cls['var'][k], 1.0/float(Cls['n'][k]))
						
			# Uncompress ave and var images if the mask is used
			if mask != None:
				Cls['ave'][k] = Util.reconstitute_image_mask(Cls['ave'][k], mask)
				Cls['var'][k] = Util.reconstitute_image_mask(Cls['var'][k], mask)
	
	# prepare the results
	if CTF:
		# ifft
		for k in xrange(K):
			fftip(Cls['ave'][k])
			fftip(Cls['var'][k])
			Cls['ave'][k].postift_depad_corner_inplace()
			Cls['var'][k].postift_depad_corner_inplace()
	
	# information display
	running_time(t_start)
	print_msg('Criterion = %11.6e \n' % Je)
	for k in xrange(K):	print_msg('Cls[%i]: %i\n'%(k, Cls['n'][k]))

	# to debug
	if DEBUG: print Cls['n']
	
	# return Cls, assign
	return Cls, assign


# K-means with SSE method
def k_means_SSE(im_M, mask, K, rand_seed, maxit, trials, CTF, F=0, T0=0, SA2=False, DEBUG=False):
	from utilities    import model_blank, get_im, running_time
	from utilities    import print_begin_msg, print_end_msg, print_msg
	from random       import seed, randint, shuffle
	from copy         import deepcopy
	from sys          import exit
	import time
	if CTF[0]:
		from filter		import filt_ctf, filt_table
		from fundamentals	import fftip

		ctf  = deepcopy(CTF[1])
		ctf2 = deepcopy(CTF[2])
		CTF  = True
	else:
		CTF  = False

	# Simulate annealing use or not
	if T0 != 0: SA = True
	else:       SA = False

	if SA:
		# for simulate annealing
		from math   import exp
		from random import random

	if mask != None:
		if isinstance(mask, basestring):
			ERROR('Mask must be an image, not a file name!', 'k-means', 1)

	N = len(im_M)

	t_start = time.time()	
	
	# Information about images
	if CTF:
		nx  = im_M[0].get_attr('or_nx')
		ny  = im_M[0].get_attr('or_ny')
		nz  = im_M[0].get_attr('or_nz')
		buf = model_blank(nx, ny, nz)
		fftip(buf)		
		nx   = im_M[0].get_xsize()
		ny   = im_M[0].get_ysize()
		nz   = im_M[0].get_zsize()
		norm = nx * ny * nz
	else:
		nx   = im_M[0].get_xsize()
		ny   = im_M[0].get_ysize()
		nz   = im_M[0].get_zsize()
		norm = nx * ny * nz
		buf  = model_blank(nx, ny, nz)

	# Variables
	if(rand_seed > 0):  seed(rand_seed)
	else:               seed()
	Cls = {}
	Cls['n']   = [0]*K     # number of objects in a given cluster
	Cls['ave'] = [0]*K     # value of cluster average
	Cls['var'] = [0]*K     # value of cluster variance
	Cls['Ji']  = [0]*K     # value of Ji
	Cls['k']   =  K	       # value of number of clusters
	Cls['N']   =  N
	assign     = [0]*N
	
        if CTF:
		Cls_ctf2   = {}
		len_ctm	   = len(ctf2[0])
	
	## TRIALS
	if trials > 1:
		MemCls, MemJe, MemAssign = {}, {}, {}
	else:
		trials = 1
	ntrials = 0
	while ntrials < trials:
		ntrials += 1

		# for simulate annealing
		if SA: T = T0
	
		# Init the cluster by an image empty
		buf.to_zero()
		for k in xrange(K):
			Cls['ave'][k] = buf.copy()
			Cls['var'][k] = buf.copy()
			Cls['n'][k]   = 0
			Cls['Ji'][k]  = 0

		## Random method
		retrial = 20
		while retrial > 0:
			retrial -= 1
			i = 0
			for im in xrange(N):
				assign[im] = randint(0, K-1)
				Cls['n'][assign[im]] += 1
			flag,k = 1,K
			while k>0 and flag:
				k -= 1 
				if Cls['n'][k] == 0:
					flag = 0
					if retrial == 0:
						ERROR('Empty class in the initialization', 'k_means_SSE', 1)
					for k in xrange(K):
						Cls['n'][k] = 0
			if flag == 1:
				retrial = 0
		
		if CTF:
			## Calculate averages ave = S CTF.F / S CTF**2, first init ctf2
			for k in xrange(K):	Cls_ctf2[k] = [0] * len_ctm
			
			for im in xrange(N):
				# compute Sum ctf2
				for i in xrange(len_ctm):	Cls_ctf2[assign[im]][i] += ctf2[im][i]
				
				# compute average first step
				CTFxF = filt_table(im_M[im], ctf[im])
				Util.add_img(Cls['ave'][assign[im]], CTFxF)
			
			for k in xrange(K):
				valCTF = [0] * len_ctm
				for i in xrange(len_ctm):	valCTF[i] = 1.0 / float(Cls_ctf2[k][i])
				Cls['ave'][k] = filt_table(Cls['ave'][k], valCTF)
			
			## Compute Ji = S(im - CTFxAve)**2 and Je = S Ji
			for n in xrange(N):
				CTFxAve		      = filt_table(Cls['ave'][assign[n]], ctf[n])
				Cls['Ji'][assign[n]] += CTFxAve.cmp("SqEuclidean", im_M[n]) / norm
			Je = 0
			for k in xrange(K):	  Je += Cls['Ji'][k]
		else:
			## Calculate averages
			for im in xrange(N):	Util.add_img(Cls['ave'][assign[im]], im_M[im])
			for k in xrange(K):	Cls['ave'][k] = Util.mult_scalar(Cls['ave'][k], 1.0/float(Cls['n'][k]))
				
			# Compute Ji = S(im - ave)**2 and Je = S Ji
			Je = 0
			for n in xrange(N):	Cls['Ji'][assign[n]] += im_M[n].cmp("SqEuclidean",Cls['ave'][assign[n]])/norm
			for k in xrange(K):	Je += Cls['Ji'][k]	
				
		## Clustering		
		ite       = 0
		watch_dog = 0
		old_Je    = 0
		change    = True
		order     = range(N)

		if DEBUG: print 'init Je', Je
		
		print_msg('\n__ Trials: %2d _________________________________%s\n'%(ntrials, time.strftime('%a_%d_%b_%Y_%H_%M_%S', time.localtime())))

		while change and watch_dog < maxit:
			ite       += 1
			watch_dog += 1
			change     = False
			shuffle(order)
			if SA: ct_pert = 0

			for imn in xrange(N):
				# to select random image
				im = order[imn]
				assign_to = -1
				
				# compute SqEuclidean (objects and centroids)
				if CTF:
					# compute the minimum distance with centroids
					# CTF: (F - CTFxAve)**2
					CTFxAve = []
					for k in xrange(K):
						tmp = filt_table(Cls['ave'][k], ctf[im])
						CTFxAve.append(tmp.copy())
					res = Util.min_dist(im_M[im], CTFxAve)
				else:
					# compute the minimum distance with centroids
					res = Util.min_dist(im_M[im], Cls['ave'])

					
				# Simulate Annealing
				if SA:
					if SA2:
						dJe = [0.0] * K
						ni  = float(Cls['n'][assign[im]])
						di  = res['dist'][assign[im]]
						for k in xrange(K):
							if k != assign[im]:
								nj  = float(Cls['n'][k])
								dj  = res['dist'][k]
								
								dJe[k] = -( (nj/(nj+1))*(dj/norm) - (ni/(ni-1))*(di/norm) )
															
							else:
								dJe[k] = 0

						# norm <0 [-1;0], >=0 [0;+1], if just 0 norm to 1
						nbneg  =  0
						nbpos  =  0
						minneg =  0
						maxpos =  0
						for k in xrange(K):
							if dJe[k] < 0.0:
								nbneg += 1
								if dJe[k] < minneg: minneg = dJe[k]
							else:
								nbpos += 1
								if dJe[k] > maxpos: maxpos = dJe[k]
						if nbneg != 0:                   dneg = -1.0 / minneg
						if nbpos != 0 and maxpos != 0:   dpos =  1.0 / maxpos
						for k in xrange(K):
							if dJe[k] < 0.0: dJe[k] = dJe[k] * dneg
							else:
								if maxpos != 0: dJe[k] = dJe[k] * dpos
								else:           dJe[k] = 1.0

						# q[k]
						q      = [0.0] * K
						arg    = [0.0] * K
						maxarg = 0
						for k in xrange(K):
							arg[k] = dJe[k] / T
							if arg[k] > maxarg: maxarg = arg[k]
						limarg = 17
						if maxarg > limarg:
							sumarg = float(sum(arg))
							for k in xrange(K): q[k] = exp(arg[k] * limarg / sumarg)
						else:
							for k in xrange(K): q[k] = exp(arg[k])
										
						# p[k]
						p = [[0.0, 0] for i in xrange(K)]
						sumq = float(sum(q))
						for k in xrange(K):
							p[k][0] = q[k] / sumq
							p[k][1] = k
											
						p.sort()
						c = [0.0] * K
						c[0] = p[0][0]
						for k in xrange(1, K): c[k] = c[k-1] + p[k][0]

						pb = random()
						select = -1
						for k in xrange(K):
							if c[k] > pb:
								select = p[k][1]
								break
					

						if select != res['pos']:
							ct_pert    += 1
							res['pos']  = select

						
					else:
						if exp( -(1) / float(T) ) > random():
							res['pos']  = randint(0, K - 1)
							ct_pert    += 1
							
			
				# moving object and update iteratively
				if res['pos'] != assign[im]:
					assign_from = assign[im]
					assign_to   = res['pos']
								
					if CTF:
						# Update average
						
						# compute valCTF = CTFi / (S ctf2 - ctf2i)
						valCTF = [0] * len_ctm
						for i in xrange(len_ctm):
							valCTF[i] = Cls_ctf2[assign_from][i] - ctf2[im][i]
							valCTF[i] = ctf[im][i] / valCTF[i]
						# compute CTFxAve
						CTFxAve = filt_table(Cls['ave'][assign_from], ctf[im])
						# compute F - CTFxAve
						buf.to_zero()
						buf = Util.subn_img(im_M[im], CTFxAve) 
						# compute valCTF * (F - CTFxAve)
						buf = filt_table(buf, valCTF)
						# sub the value at the average
						Util.sub_img(Cls['ave'][assign_from], buf)

						# compute valCTF = CTFi / (S ctf2 + ctf2i)
						valCTF = [0] * len_ctm
						for i in xrange(len_ctm):
							valCTF[i] = ctf[im][i] / (Cls_ctf2[assign_to][i] + ctf2[im][i])
						# compute CTFxAve
						CTFxAve = filt_table(Cls['ave'][assign_to], ctf[im])
						# compute F - CTFxAve
						buf.to_zero()
						buf = Util.subn_img(im_M[im], CTFxAve) 
						# compute valCTF * (F - CTFxAve)
						buf = filt_table(buf, valCTF)
						# add the value at the average
						Util.add_img(Cls['ave'][assign_to], buf)	
					else:
						# Update average
						buf.to_zero()
						buf = Util.mult_scalar(Cls['ave'][assign_from], float(Cls['n'][assign_from]))
						Util.sub_img(buf,im_M[im])
						Cls['ave'][assign_from] = Util.mult_scalar(buf, 1.0/float(Cls['n'][assign_from]-1))

						buf.to_zero()
						buf = Util.mult_scalar(Cls['ave'][assign_to], float(Cls['n'][assign_to]))
						Util.add_img(buf, im_M[im])
						Cls['ave'][assign_to] = Util.mult_scalar(buf, 1.0/float(Cls['n'][assign_to]+1))
				
					
					# new number of objects in clusters
					Cls['n'][assign_from] -= 1
					assign[im]             = assign_to
					Cls['n'][assign_to]   += 1
					if CTF:
						# update Sum ctf2
						for i in xrange(len_ctm):
							Cls_ctf2[assign_from][i] -= ctf2[im][i]
							Cls_ctf2[assign_to][i]   += ctf2[im][i]
										
					# empty cluster control
					if Cls['n'][assign_from] == 1:
						print_msg('\nEmpty clusters!!\n')
						for k in xrange(K):	print_msg('Cls[%i]: %i\n'%(k, Cls['n'][k]))
						ERROR("Empty groups in kmeans_SSE","k_means_SSE",1)						
						exit()							
								
					change = True
			
			if CTF:
				## Compute Ji = S(im - CTFxAve)**2 and Je = S Ji
				for k in xrange(K): Cls['Ji'][k] = 0
				for n in xrange(N):
					CTFxAve		      = filt_table(Cls['ave'][assign[n]], ctf[n])
					Cls['Ji'][assign[n]] += CTFxAve.cmp("SqEuclidean", im_M[n]) / norm
				Je = 0
				for k in xrange(K):	  Je += Cls['Ji'][k]
			else:
				# Compute Je
				Je = 0
				for k in xrange(K):     Cls['Ji'][k] = 0
				for n in xrange(N):	Cls['Ji'][assign[n]] += im_M[n].cmp("SqEuclidean",Cls['ave'][assign[n]]) / norm
				for k in xrange(K):	Je += Cls['Ji'][k]
	
			# threshold convergence control
			if Je != 0:
				thd = abs(Je - old_Je) / Je
			else:
				thd = 0

			if SA:
				if thd < 1e-12 and ct_pert == 0: watch_dog = maxit
			else:
				if thd < 1e-8:	                 watch_dog = maxit
			old_Je = Je

			if SA:
				# Simulate annealing, update temperature
				T *= F
				if SA2:
					if T < 0.09: SA = False

				print_msg('> iteration: %5d    criterion: %11.6e    T: %13.8f  ct disturb: %5d\n' % (ite, Je, T, ct_pert))
				if DEBUG: print '> iteration: %5d    criterion: %11.6e    T: %13.8f  ct disturb: %5d' % (ite, Je, T, ct_pert)
			else:
				print_msg('> iteration: %5d    criterion: %11.6e\n'%(ite, Je))
				if DEBUG: print '> iteration: %5d    criterion: %11.6e'%(ite, Je)

		if CTF:
			## Calculate averages ave = S CTF.F / S CTF**2, first init ctf2
			for k in xrange(K):	Cls_ctf2[k] = [0] * len_ctm
			for im in xrange(N):
				# compute Sum ctf2
				for i in xrange(len_ctm):	Cls_ctf2[assign[im]][i] += ctf2[im][i]
				# compute average first step
				CTFxF = filt_table(im_M[im], ctf[im])
				Util.add_img(Cls['ave'][assign[im]], CTFxF)
			for k in xrange(K):
				valCTF = [0] * len_ctm
				for i in xrange(len_ctm):	valCTF[i] = 1.0 / float(Cls_ctf2[k][i])
				Cls['ave'][k] = filt_table(Cls['ave'][k], valCTF)
			## Compute Ji = S(im - CTFxAve)**2 and Je = S Ji
			for k in xrange(K): Cls['Ji'][k] = 0
			for n in xrange(N):
				CTFxAve		      = filt_table(Cls['ave'][assign[n]], ctf[n])
				Cls['Ji'][assign[n]] += CTFxAve.cmp("SqEuclidean", im_M[n]) / norm
			Je = 0
			for k in xrange(K):	  Je += Cls['Ji'][k]
		else:
			# Calculate the real averages, because the iterations method cause approximation
			buf.to_zero()
			for k in xrange(K):     Cls['ave'][k] = buf.copy()
			for im in xrange(N):	Util.add_img(Cls['ave'][assign[im]], im_M[im])
			for k in xrange(K):	Cls['ave'][k] = Util.mult_scalar(Cls['ave'][k], 1.0/float(Cls['n'][k]))

			# Compute the accurate Je, because during the iterations Je is aproximated from average
			Je = 0
			for k in xrange(K):     Cls['Ji'][k] = 0
			for n in xrange(N):	Cls['Ji'][assign[n]] += im_M[n].cmp("SqEuclidean",Cls['ave'][assign[n]]) / norm
			for k in xrange(K):	Je += Cls['Ji'][k]	

		# memorize the result for this trial	
		if trials > 1:
			MemCls[ntrials-1]    = deepcopy(Cls)
			MemJe[ntrials-1]     = deepcopy(Je)
			MemAssign[ntrials-1] = deepcopy(assign)

			
	# if severals trials choose the best
	if trials > 1:
		val_min = 1.0e20
		best    = -1
		for n in xrange(trials):
			if MemJe[n] < val_min:
				val_min = MemJe[n]
				best    = n
		# affect the best
		Cls    = MemCls[best]
		Je     = MemJe[best]
		assign = MemAssign[best]
		
	if CTF:
		# compute the variance S (F - CTF * Ave)**2
		buf.to_zero()
		for k in xrange(K): Cls['var'][k] = buf.copy()
		
		for n in xrange(N):
			CTFxAve = filt_table(Cls['ave'][assign[n]], ctf[n])
			
			buf.to_zero()
			buf     = Util.subn_img(im_M[n], CTFxAve)
			Util.add_img(Cls['var'][assign[n]], buf) ## **2
		
	else:
		# compute the variance 1/n S(im-ave)**2 -> 1/n (Sim**2 - n ave**2)
		for im in xrange(N):	Util.add_img2(Cls['var'][assign[im]], im_M[im])
		for k in xrange(K):
			buf.to_zero()
			Util.add_img2(buf, Cls['ave'][k])
			Cls['var'][k] = Util.madn_scalar(Cls['var'][k], buf, -float(Cls['n'][k]))
			Util.mul_scalar(Cls['var'][k], 1.0/float(Cls['n'][k]))
			
			# Uncompress ave and var images if the mask is used
			if mask != None:
				Cls['ave'][k] = Util.reconstitute_image_mask(Cls['ave'][k], mask)
				Cls['var'][k] = Util.reconstitute_image_mask(Cls['var'][k], mask)

	# write the results if out_dire is defined
	if CTF:
		# ifft
		for k in xrange(K):
			fftip(Cls['ave'][k])
			fftip(Cls['var'][k])
			Cls['ave'][k].postift_depad_corner_inplace()
			Cls['var'][k].postift_depad_corner_inplace()

	# information display
	running_time(t_start)
	print_msg('Criterion = %11.6e \n' % Je)
	for k in xrange(K):	print_msg('Cls[%i]: %i\n'%(k, Cls['n'][k]))
	
	# to debug
	if DEBUG: print Cls['n']
		
	# return Cls, assign and Je
	return Cls, assign


# K-means MPI with classical method
def k_means_cla_MPI(im_M, mask, K, rand_seed, maxit, trials, CTF, myid, main_node, N_start, N_stop, F=0, T0=0, SA2=False):
	from utilities    import model_blank, get_im
	from utilities    import bcast_EMData_to_all, reduce_EMData_to_root
	from utilities    import print_msg, running_time
	from random       import seed, randint
	from copy	  import deepcopy
	from mpi 	  import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_recv, mpi_send
	from mpi 	  import MPI_SUM, MPI_FLOAT, MPI_INT, MPI_LOR
	import time
	import sys
	if CTF[0]:
		from filter		import filt_ctf, filt_table
		from fundamentals	import fftip

		ctf  = deepcopy(CTF[1])
		ctf2 = deepcopy(CTF[2])
		CTF  = True
	else:
		CTF  = False

	# Simulate annealing
	if F != 0: SA = True
	else:      SA = False

	if SA:
		from math   import exp
		from random import random

	if mask != None:
		if isinstance(mask, basestring):
			ERROR('Mask must be an image, not a file name!', 'k-means', 1)

	N = len(im_M)

	# [id]   part of code different for each node
	# [sync] synchronise each node
	# [main] part of code just for the main node
	# [all]  code write for all node

	t_start = time.time()
	
	# [all] Informations on images or mask for the norm
	if CTF:
		nx  = im_M[N_start].get_attr('or_nx')
		ny  = im_M[N_start].get_attr('or_ny')
		nz  = im_M[N_start].get_attr('or_nz')
		buf = model_blank(nx, ny, nz)
		fftip(buf)		
		nx   = im_M[N_start].get_xsize()
		ny   = im_M[N_start].get_ysize()
		nz   = im_M[N_start].get_zsize()
		norm = nx * ny * nz
	else:
		nx   = im_M[N_start].get_xsize()
		ny   = im_M[N_start].get_ysize()
		nz   = im_M[N_start].get_zsize()
		norm = nx * ny * nz
		buf  = model_blank(nx, ny, nz)
	
	# [all] define parameters
	if rand_seed > 0: seed(rand_seed)
	else:             seed()
	Cls={}
	Cls['n']   = [0]*K   # number of objects in a given cluster
	Cls['ave'] = [0]*K   # value of cluster average
	Cls['var'] = [0]*K   # value of cluster variance
	Cls['Ji']  = [0]*K   # value of ji
	Cls['k']   =  K	     # value of number of clusters
	Cls['N']   =  N
	assign     = [0]*N
	
	if CTF:
		Cls_ctf2    = {}
		len_ctm	    = len(ctf2[0])
	
	# TRIALS
	if trials > 1:
		MemCls, MemJe, MemAssign = {}, {}, {}
	else:
		trials = 1
	ntrials   = 0
	while ntrials < trials:
		ntrials  += 1

		# Simulate annealing
		if SA: T = T0
		
		# [all] Init the cluster by an image empty
		buf.to_zero()
		for k in xrange(K):
			Cls['ave'][k] = buf.copy()
			Cls['var'][k] = buf.copy()
			Cls['Ji'][k]  = 0
			Cls['n'][k]   = 0
			OldClsn       = [0] * K

		## [main] Random method
		if myid == main_node:
			retrial = 20
			while retrial > 0:
				retrial -= 1
				i = 0
				for im in xrange(N):
					assign[im] = randint(0, K-1)
					Cls['n'][int(assign[im])] += 1
				flag,k = 1,0
				while k < K and flag:
					if Cls['n'][k] == 0:
						flag = 0
						if retrial == 0:
							ERROR('Empty class in the initialization', 'k_means_cla_MPI', 1)
						for k in xrange(K):
							Cls['n'][k] = 0
					k += 1
				if flag == 1: retrial = 0

		# [sync] waiting the assign is finished
		mpi_barrier(MPI_COMM_WORLD)

		# [all] send assign to the others proc and the number of objects in each clusters
		assign = mpi_bcast(assign, N, MPI_INT, main_node, MPI_COMM_WORLD)
		assign = assign.tolist()     # convert array gave by MPI to list
		Cls['n'] = mpi_bcast(Cls['n'], K, MPI_FLOAT, main_node, MPI_COMM_WORLD)
		Cls['n'] = Cls['n'].tolist() # convert array gave by MPI to list
		
		## 
		if CTF:
			# [all] first init ctf2
			for k in xrange(K):	Cls_ctf2[k] = [0] * len_ctm
			
			# [id] compute local S ctf2 and local S ave	
			for im in xrange(N_start, N_stop):
				# ctf2
				for i in xrange(len_ctm):
					Cls_ctf2[int(assign[im])][i] += ctf2[im][i]
				# ave
				CTFxF = filt_table(im_M[im], ctf[im])
				Util.add_img(Cls['ave'][int(assign[im])], CTFxF)
			
			# [sync] waiting the result
			mpi_barrier(MPI_COMM_WORLD)
			
			# [all] compute global sum, broadcast the results and obtain the average ave = S CTF.F / S CTF**2
			for k in xrange(K):
				Cls_ctf2[k] = mpi_reduce(Cls_ctf2[k], len_ctm, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
				Cls_ctf2[k] = mpi_bcast(Cls_ctf2[k],  len_ctm, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				Cls_ctf2[k] = Cls_ctf2[k].tolist()    # convert array gave by MPI to list
				
				reduce_EMData_to_root(Cls['ave'][k], myid, main_node)
				bcast_EMData_to_all(Cls['ave'][k], myid, main_node)
				
				for i in xrange(len_ctm):	Cls_ctf2[k][i] = 1.0 / Cls_ctf2[k][i]
				Cls['ave'][k] = filt_table(Cls['ave'][k], Cls_ctf2[k])
								
		else:
			# [id] Calculates averages, first calculate local sum
			for im in xrange(N_start, N_stop):	Util.add_img(Cls['ave'][int(assign[im])], im_M[im])

			# [sync] waiting the result
			mpi_barrier(MPI_COMM_WORLD)

			# [all] compute global sum, broadcast the results and obtain the average
			for k in xrange(K):
				reduce_EMData_to_root(Cls['ave'][k], myid, main_node) 
				bcast_EMData_to_all(Cls['ave'][k], myid, main_node)
				Cls['ave'][k] = Util.mult_scalar(Cls['ave'][k], 1.0/float(Cls['n'][k]))
		
		## Clustering		
		ite       = 0
		watch_dog = 0
		old_Je    = 0
		change    = 1
		if myid == main_node: print_msg('\n__ Trials: %2d _________________________________%s\n'%(ntrials, time.strftime('%a_%d_%b_%Y_%H_%M_%S', time.localtime())))
		
		while change and watch_dog < maxit:
			ite       += 1
			watch_dog += 1
			change     = 0
			Je         = 0
			if SA:
			   ct_pert = 0
			
			# [all] init the number of objects
			for k in xrange(K):	Cls['n'][k]=0

			# [id] assign each images with err_min between all clusters averages
			for im in xrange(N_start, N_stop):

				# [all] compute min dist between object and centroids
				if CTF:
					CTFxAve = []
					for k in xrange(K):
						tmp = filt_table(Cls['ave'][k], ctf[im])
						CTFxAve.append(tmp.copy())
					res = Util.min_dist(im_M[im], CTFxAve)
				else:
					res = Util.min_dist(im_M[im], Cls['ave'])

				# [all] Simulate annealing
				if SA:
					if SA2:
						dJe = [0.0] * K
						ni  = float(Cls['n'][assign[im]])
						di  = res['dist'][assign[im]]
						for k in xrange(K):
							if k != assign[im]:
								nj  = float(Cls['n'][k])
								dj  = res['dist'][k]
								
								dJe[k] = -( (nj/(nj+1))*(dj/norm) - (ni/(ni-1))*(di/norm) )
															
							else:
								dJe[k] = 0

						# norm <0 [-1;0], >=0 [0;+1], if just 0 norm to 1
						nbneg  =  0
						nbpos  =  0
						minneg =  0
						maxpos =  0
						for k in xrange(K):
							if dJe[k] < 0.0:
								nbneg += 1
								if dJe[k] < minneg: minneg = dJe[k]
							else:
								nbpos += 1
								if dJe[k] > maxpos: maxpos = dJe[k]
						if nbneg != 0:                   dneg = -1.0 / minneg
						if nbpos != 0 and maxpos != 0:   dpos =  1.0 / maxpos
						for k in xrange(K):
							if dJe[k] < 0.0: dJe[k] = dJe[k] * dneg
							else:
								if maxpos != 0: dJe[k] = dJe[k] * dpos
								else:           dJe[k] = 1.0

						# q[k]
						q      = [0.0] * K
						arg    = [0.0] * K
						maxarg = 0
						for k in xrange(K):
							arg[k] = dJe[k] / T
							if arg[k] > maxarg: maxarg = arg[k]
						limarg = 17
						if maxarg > limarg:
							sumarg = float(sum(arg))
							for k in xrange(K): q[k] = exp(arg[k] * limarg / sumarg)
						else:
							for k in xrange(K): q[k] = exp(arg[k])
										
						# p[k]
						p = [[0.0, 0] for i in xrange(K)]
						sumq = float(sum(q))
						for k in xrange(K):
							p[k][0] = q[k] / sumq
							p[k][1] = k
											
						p.sort()
						c = [0.0] * K
						c[0] = p[0][0]
						for k in xrange(1, K): c[k] = c[k-1] + p[k][0]

						pb = random()
						select = -1
						for k in xrange(K):
							if c[k] > pb:
								select = p[k][1]
								break
					

						if select != res['pos']:
							ct_pert    += 1
							res['pos']  = select

						
					else:
						if exp( -(1) / float(T) ) > random():
							res['pos']  = randint(0, K - 1)
							ct_pert    += 1
				
				# [all] move object
				if res['pos'] != assign[im]:
					assign[im]  = res['pos']
					change      = 1
					
				# [id] calculate the new local number of objects
				Cls['n'][int(assign[im])] += 1
				
			# [sync] waiting the result
			mpi_barrier(MPI_COMM_WORLD)

			# [all] sum the number of objects in each node and broadcast
			Cls['n'] = mpi_reduce(Cls['n'], K, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
			Cls['n'] = mpi_bcast(Cls['n'], K, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			Cls['n'] = Cls['n'].tolist() # convert array gave by MPI to list
			
			# [all] init average and ctf2
			for k in xrange(K):
				if Cls['n'][k] < 1:
					ERROR('Empty groups in kmeans_classical', 'k_means_cla MPI', 1)
					exit()			
				
				Cls['ave'][k].to_zero()
				Cls['Ji'][k] = 0
				if CTF:	Cls_ctf2[k] = [0] * len_ctm			
			
			if CTF:
				# [id] compute local S ctf2 and local S ave	
				for im in xrange(N_start, N_stop):
					# ctf2
					for i in xrange(len_ctm):
						Cls_ctf2[int(assign[im])][i] += ctf2[im][i]
					# ave
					CTFxF = filt_table(im_M[im], ctf[im])
					Util.add_img(Cls['ave'][int(assign[im])], CTFxF)
				
				# [sync] waiting the result
				mpi_barrier(MPI_COMM_WORLD)
				
				# [all] compute global sum, broadcast the results and obtain the average ave = S CTF.F / S CTF**2
				for k in xrange(K):
					Cls_ctf2[k] = mpi_reduce(Cls_ctf2[k], len_ctm, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
					Cls_ctf2[k] = mpi_bcast(Cls_ctf2[k], len_ctm, MPI_FLOAT, main_node, MPI_COMM_WORLD)
					Cls_ctf2[k] = Cls_ctf2[k].tolist() # convert array gave by MPI to list

					reduce_EMData_to_root(Cls['ave'][k], myid, main_node)
					bcast_EMData_to_all(Cls['ave'][k], myid, main_node)
					
					for i in xrange(len_ctm):	Cls_ctf2[k][i] = 1.0 / float(Cls_ctf2[k][i])
					Cls['ave'][k] = filt_table(Cls['ave'][k], Cls_ctf2[k])

				# [id] compute Ji
				for im in xrange(N_start, N_stop):
					CTFxAve = filt_table(Cls['ave'][int(assign[im])], ctf[im])
					Cls['Ji'][int(assign[im])] += CTFxAve.cmp("SqEuclidean", im_M[im]) / norm

				# [all] waiting the result
				mpi_barrier(MPI_COMM_WORLD)

				# [all] global sum Ji
				Cls['Ji'] = mpi_reduce(Cls['Ji'], K, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
				Cls['Ji'] = mpi_bcast(Cls['Ji'],  K, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				Cls['Ji'] = Cls['Ji'].tolist()
			
			else:			
				# [id] Update clusters averages
				for im in xrange(N_start, N_stop):	Util.add_img(Cls['ave'][int(assign[im])], im_M[im])

				# [sync] waiting the result
				mpi_barrier(MPI_COMM_WORLD)

				# [all] compute global sum, broadcast the results and obtain the average
				for k in xrange(K):
					reduce_EMData_to_root(Cls['ave'][k], myid, main_node) 
					bcast_EMData_to_all(Cls['ave'][k], myid, main_node)
					Cls['ave'][k] = Util.mult_scalar(Cls['ave'][k], 1.0/float(Cls['n'][k]))

				# [id] compute Ji
				for im in xrange(N_start, N_stop): Cls['Ji'][int(assign[im])] += im_M[im].cmp("SqEuclidean", Cls['ave'][int(assign[im])])/norm


			# [all] compute Je
			Je = 0
			for k in xrange(K): Je += Cls['Ji'][k]

			# [all] waiting the result
			mpi_barrier(MPI_COMM_WORLD)

			# [all] calculate Je global sum and broadcast
			Je = mpi_reduce(Je, 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
			Je = mpi_bcast(Je, 1, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			Je = Je.tolist()
			Je = Je[0]
			
			if SA:
				# Simulate annealing, update temperature
				T *= F
				if SA2:
					if T < 0.09: SA = False

				#[id] informations display
				if myid == main_node: print_msg('> iteration: %5d    criterion: %11.6e   T: %13.8f  disturb:  %5d\n' % (ite, Je, T, ct_pert))
			else:
				# [id] informations display
				if myid == main_node: print_msg('> iteration: %5d    criterion: %11.6e\n' % (ite, Je))
						
			# threshold convergence control
			if Je != 0: thd = abs(Je - old_Je) / Je
			else:       thd = 0
			if SA:
				if thd < 1e-12 and ct_pert == 0: change = 0
			else:
				if thd < 1e-8:	change = 0
			old_Je = Je
			
			# [all] Need to broadcast this value because all node must run together
			change = mpi_reduce(change, 1, MPI_INT, MPI_LOR, main_node, MPI_COMM_WORLD)
			change = mpi_bcast(change, 1, MPI_INT, main_node, MPI_COMM_WORLD)
			change = change.tolist()
			change = change[0]
		
		# [all] waiting the result
		mpi_barrier(MPI_COMM_WORLD)

		# [id] memorize the result for this trial	
		if trials > 1:
			MemCls[ntrials-1]    = deepcopy(Cls)
			MemJe[ntrials-1]     = deepcopy(Je)
			MemAssign[ntrials-1] = deepcopy(assign)
			
	# if severals trials choose the best
	if trials > 1:
		val_min = 1.0e20
		best    = -1
		for n in xrange(trials):
			if MemJe[n] < val_min:
				val_min = MemJe[n]
				best    = n
		# affect the best
		Cls    = MemCls[best]
		Je     = MemJe[best]
		assign = MemAssign[best]
	
	if CTF:
		# [id] compute Ji and the variance S (F - CTFxAve)**2
		for im in xrange(N_start, N_stop):
			CTFxAve = filt_table(Cls['ave'][int(assign[im])], ctf[im])
			Cls['Ji'][int(assign[im])] += CTFxAve.cmp("SqEuclidean", im_M[im]) / norm
			
			buf.to_zero()
			buf = Util.subn_img(CTFxAve, im_M[im])
			Util.add_img(Cls['var'][int(assign[im])], buf) # **2
		
		# [all] waiting the result
		mpi_barrier(MPI_COMM_WORLD)
		
		# [all] global sum Ji and var
		Cls['Ji'] = mpi_reduce(Cls['Ji'], K, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
		Cls['Ji'] = mpi_bcast(Cls['Ji'],  K, MPI_FLOAT, main_node, MPI_COMM_WORLD)
		Cls['Ji'] = Cls['Ji'].tolist()
		for k in xrange(K):
			reduce_EMData_to_root(Cls['var'][k], myid, main_node)
			
	else:
		# [id] compute Ji and the variance 1/n S(im-ave)**2 -> 1/n (Sim**2 - n ave**2)	
		for im in xrange(N_start, N_stop):
			Cls['Ji'][int(assign[im])] += im_M[im].cmp("SqEuclidean", Cls['ave'][int(assign[im])])/norm		
			Util.add_img2(Cls['var'][int(assign[im])], im_M[im])
		
		# [all] waiting the result
		mpi_barrier(MPI_COMM_WORLD)

		# [all] global sum ji and im**2
		Cls['Ji'] = mpi_reduce(Cls['Ji'], K, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
		Cls['Ji'] = mpi_bcast(Cls['Ji'],  K, MPI_FLOAT, main_node, MPI_COMM_WORLD)
		Cls['Ji'] = Cls['Ji'].tolist()
		
		for k in xrange(K): reduce_EMData_to_root(Cls['var'][k], myid, main_node)	
		
		
		# [main] caclculate the variance for each cluster
		if myid == main_node:
			for k in xrange(K):
				buf.to_zero()
				Util.add_img2(buf, Cls['ave'][k])
				Cls['var'][k] = Util.madn_scalar(Cls['var'][k], buf, -float(Cls['n'][k]))
				Util.mul_scalar(Cls['var'][k], 1.0/float(Cls['n'][k]))
				
				# Uncompress ave and var images if the mask is used
				if mask != None:
					Cls['ave'][k] = Util.reconstitute_image_mask(Cls['ave'][k], mask)
					Cls['var'][k] = Util.reconstitute_image_mask(Cls['var'][k], mask)

	# [id] prepare assign to update
	v = range(N_start, N_stop)
	for n in xrange(N_start):
		assign[n] = 0
	for n in xrange(N_stop, N):
		assign[n] = 0
	# [all] gather in main_node
	assign = mpi_reduce(assign, N, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
	assign = assign.tolist() # convert array gave by MPI to list

	# [main_node] write the result
	if myid == main_node and CTF:
		# ifft
		for k in xrange(K):
			fftip(Cls['ave'][k])
			fftip(Cls['var'][k])
			Cls['ave'][k].postift_depad_corner_inplace()
			Cls['var'][k].postift_depad_corner_inplace()

			
	# [main_node] information display
	if myid == main_node:
		running_time(t_start)
		print_msg('Criterion = %11.4e \n' % Je)
		for k in xrange(K):	print_msg('Cls[%i]: %i\n'%(k, Cls['n'][k]))

	# [all] waiting all nodes
	mpi_barrier(MPI_COMM_WORLD)
		
	if myid == main_node: return Cls, assign
	else:                 return None, None

# K-means MPI with SSE method
def k_means_SSE_MPI(im_M, mask, K, rand_seed, maxit, trials, CTF, myid, main_node, ncpu, N_start, N_stop, F=0, T0=0, SA2=False):
	from utilities    import model_blank, get_im
	from utilities    import bcast_EMData_to_all, reduce_EMData_to_root
	from utilities    import print_msg, running_time
	from random       import seed, randint, shuffle
	from copy	  import deepcopy
	from mpi 	  import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_recv, mpi_send
	from mpi 	  import MPI_SUM, MPI_FLOAT, MPI_INT, MPI_LOR
	import time
	import sys
	if CTF[0]:
		from filter		import filt_ctf, filt_table
		from fundamentals	import fftip

		ctf  = deepcopy(CTF[1])
		ctf2 = deepcopy(CTF[2])
		CTF  = True
	else:
		CTF  = False

	# Simulate annealing
	if F != 0: SA = True
	else:      SA = False

	if SA:
		from math   import exp
		from random import random

	if mask != None:
		if isinstance(mask, basestring):
			ERROR('Mask must be an image, not a file name!', 'k-means', 1)

	N = len(im_M)
	
	# [id]   part of code different for each node
	# [sync] synchronise each node
	# [main] part of code just for the main node
	# [all]  code write for all node
	
	t_start = time.time()	

	# [all] Informations on images or mask for the norm
	if CTF:
		nx  = im_M[N_start].get_attr('or_nx')
		ny  = im_M[N_start].get_attr('or_ny')
		nz  = im_M[N_start].get_attr('or_nz')
		buf = model_blank(nx, ny, nz)
		fftip(buf)		
		nx   = im_M[N_start].get_xsize()
		ny   = im_M[N_start].get_ysize()
		nz   = im_M[N_start].get_zsize()
		norm = nx * ny * nz
	else:
		nx   = im_M[N_start].get_xsize()
		ny   = im_M[N_start].get_ysize()
		nz   = im_M[N_start].get_zsize()
		norm = nx * ny * nz
		buf  = model_blank(nx, ny, nz)
	
	# [all] define parameters	
	if rand_seed > 0: seed(rand_seed)
	else:             seed()
	Cls={}
	Cls['n']   = [0]*K   # number of objects in a given cluster
	Cls['ave'] = [0]*K   # value of cluster average
	Cls['var'] = [0]*K   # value of cluster variance
	Cls['Ji']  = [0]*K   # value of ji
	Cls['k']   =  K	     # value of number of clusters
	Cls['N']   =  N
	assign     = [0]*N
	
	if CTF:
		Cls_ctf2    = {}
		len_ctm	    = len(ctf2[0])
	
	# TRIALS
	if trials > 1:
		MemCls, MemJe, MemAssign = {}, {}, {}
	else:
		trials = 1
	ntrials   = 0
	while ntrials < trials:
		ntrials += 1

		# Simulate annealing
		if SA: T = T0
		
		# [all] Init the cluster by an image empty
		buf.to_zero()
		for k in xrange(K):
			Cls['ave'][k] = buf.copy()
			Cls['var'][k] = buf.copy()
			Cls['Ji'][k]  = 0
			Cls['n'][k]   = 0
			OldClsn       = [0] * K
		
		## [main] Random method
		if myid == main_node:
			retrial = 20
			while retrial > 0:
				retrial -= 1
				i = 0
				for im in xrange(N):
					assign[im] = randint(0, K-1)
					Cls['n'][int(assign[im])] += 1
				flag,k = 1,0
				while k < K and flag:
					if Cls['n'][k] == 0:
						flag = 0
						if retrial == 0:
							ERROR('Empty class in the initialization', 'k_means_SSE_MPI', 1)
						for k in xrange(K):
							Cls['n'][k] = 0
					k += 1
				if flag == 1:	retrial = 0
		
		# [sync] waiting the assign is finished
		mpi_barrier(MPI_COMM_WORLD)
		
		# [all] send assign to the others proc and the number of objects in each cluster
		assign   = mpi_bcast(assign, N, MPI_INT, main_node, MPI_COMM_WORLD)
		assign   = assign.tolist()     # convert array gave by MPI to list
		Cls['n'] = mpi_bcast(Cls['n'], K, MPI_FLOAT, main_node, MPI_COMM_WORLD)
		Cls['n'] = Cls['n'].tolist() # convert array gave by MPI to list
			
		## 
		if CTF:
			# [all] first init ctf2
			for k in xrange(K):	Cls_ctf2[k] = [0] * len_ctm
			
			# [id] compute local S ctf2 and local S ave	
			for im in xrange(N_start, N_stop):
				# ctf2
				for i in xrange(len_ctm):
					Cls_ctf2[int(assign[im])][i] += ctf2[im][i]
				# ave
				CTFxF = filt_table(im_M[im], ctf[im])
				Util.add_img(Cls['ave'][int(assign[im])], CTFxF)
			
			# [sync] waiting the result
			mpi_barrier(MPI_COMM_WORLD)
			
			# [all] compute global sum, broadcast the results and obtain the average ave = S CTF.F / S CTF**2
			for k in xrange(K):
				Cls_ctf2[k] = mpi_reduce(Cls_ctf2[k], len_ctm, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
				Cls_ctf2[k] = mpi_bcast(Cls_ctf2[k],  len_ctm, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				Cls_ctf2[k] = Cls_ctf2[k].tolist()    # convert array gave by MPI to list
				
				reduce_EMData_to_root(Cls['ave'][k], myid, main_node)
				bcast_EMData_to_all(Cls['ave'][k], myid, main_node)
				
				valCTF = [0] * len_ctm
				for i in xrange(len_ctm):	valCTF[i] = 1.0 / Cls_ctf2[k][i]
				Cls['ave'][k] = filt_table(Cls['ave'][k], valCTF)
			
			# [all] waiting the result
			mpi_barrier(MPI_COMM_WORLD)
										
		else:
			# [id] Calculates averages, first calculate local sum
			for im in xrange(N_start, N_stop):	Util.add_img(Cls['ave'][int(assign[im])], im_M[im])

			# [sync] waiting the result
			mpi_barrier(MPI_COMM_WORLD)

			# [all] compute global sum, broadcast the results and obtain the average
			for k in xrange(K):
				reduce_EMData_to_root(Cls['ave'][k], myid, main_node) 
				bcast_EMData_to_all(Cls['ave'][k], myid, main_node)
				Cls['ave'][k] = Util.mult_scalar(Cls['ave'][k], 1.0/float(Cls['n'][k]))

			# [sync] waiting the result
			mpi_barrier(MPI_COMM_WORLD)

	
		## Clustering
		ite       = 0
		watch_dog = 0
		old_Je    = 0
		sway_Je   = [0] * 2
		sway_ct   = 0
		change    = 1
		loc_order = range(N_start, N_stop)
		if myid == main_node: print_msg('\n__ Trials: %2d _________________________________%s\n'%(ntrials, time.strftime('%a_%d_%b_%Y_%H_%M_%S', time.localtime())))
		
		th_update = 200 #N / ncpu / 5
		while change and watch_dog < maxit:
			ite       += 1
			watch_dog += 1
			change     = 0
			err	   = {}
			order	   = [0]*N
			index	   = 0
			Je	   = 0
			if SA:
			   ct_pert = 0
			shuffle(loc_order)	
				
			# [id] random number to image
			for n in xrange(N_start, N_stop):
				order[n] = loc_order[index]
				index   += 1
			imn = N_start
			all_proc_end = 0

			ct_update = 0
			while imn < N_stop:
				# to select random image
				im        = order[imn]
				imn      +=  1
				
				# [all] compute min dist between object and centroids
				if CTF:
					CTFxAve = []
					for k in xrange(K):
						tmp = filt_table(Cls['ave'][k], ctf[im])
						CTFxAve.append(tmp.copy())
					res = Util.min_dist(im_M[im], CTFxAve)
					
				else:
					res = Util.min_dist(im_M[im], Cls['ave'])

					
				# [all] Simulate annealing
				if SA:
					if SA2:
						dJe = [0.0] * K
						ni  = float(Cls['n'][assign[im]])
						di  = res['dist'][assign[im]]											
						for k in xrange(K):
							if k != assign[im]:
								nj  = float(Cls['n'][k])
								dj  = res['dist'][k]
								
								dJe[k] = -( (nj/(nj+1))*(dj/norm) - (ni/(ni-1))*(di/norm) )
															
							else:
								dJe[k] = 0

						# norm <0 [-1;0], >=0 [0;+1], if just 0 norm to 1
						nbneg  =  0
						nbpos  =  0
						minneg =  0
						maxpos =  0
						for k in xrange(K):
							if dJe[k] < 0.0:
								nbneg += 1
								if dJe[k] < minneg: minneg = dJe[k]
							else:
								nbpos += 1
								if dJe[k] > maxpos: maxpos = dJe[k]
						if nbneg != 0:                   dneg = -1.0 / minneg
						if nbpos != 0 and maxpos != 0:   dpos =  1.0 / maxpos
						for k in xrange(K):
							if dJe[k] < 0.0: dJe[k] = dJe[k] * dneg
							else:
								if maxpos != 0: dJe[k] = dJe[k] * dpos
								else:           dJe[k] = 1.0

						# q[k]
						q      = [0.0] * K
						arg    = [0.0] * K
						maxarg = 0
						for k in xrange(K):
							arg[k] = dJe[k] / T
							if arg[k] > maxarg: maxarg = arg[k]
						limarg = 17
						if maxarg > limarg:
							sumarg = float(sum(arg))
							for k in xrange(K): q[k] = exp(arg[k] * limarg / sumarg)
						else:
							for k in xrange(K): q[k] = exp(arg[k])

						# p[k]
						p = [[0.0, 0] for i in xrange(K)]
						sumq = float(sum(q))
						for k in xrange(K):
							p[k][0] = q[k] / sumq
							p[k][1] = k
											
						p.sort()
						c = [0.0] * K
						c[0] = p[0][0]
						for k in xrange(1, K): c[k] = c[k-1] + p[k][0]

						pb = random()
						select = -1
						for k in xrange(K):
							if c[k] > pb:
								select = p[k][1]
								break
					

						if select != res['pos']:
							ct_pert    += 1
							res['pos']  = select

						
					else:
						if exp( -(1) / float(T) ) > random():
							res['pos']  = randint(0, K - 1)
							ct_pert    += 1
					
				# [all] moving object
				if res['pos'] != assign[im]:
					assign[im] = res['pos']
					change     = 1
				'''
				# count update
				ct_update += 1

				# if reach th update all (SSE)
				if ct_update > th_update:
					ct_update = 0
					
					# [id] compute the number of objects
					for k in xrange(K): 		  Cls['n'][k] = 0
					for n in xrange(N_start, N_stop): Cls['n'][int(assign[n])] += 1			

					# [all] sum the number of objects in each node and broadcast
					Cls['n'] = mpi_reduce(Cls['n'], K, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
					Cls['n'] = mpi_bcast(Cls['n'], K, MPI_FLOAT, main_node, MPI_COMM_WORLD)
					Cls['n'] = Cls['n'].tolist() # convert array gave by MPI to list

					# [all] init average, Ji and ctf2
					for k in xrange(K):
						if Cls['n'][k] < 2:
							ERROR('Empty groups in kmeans_classical', 'k_means_SSE_MPI', 1)
							exit()	

						Cls['ave'][k].to_zero()
						if CTF:	Cls_ctf2[k] = [0] * len_ctm	

					if CTF:
						# [id] compute local S ctf2 and local S ave	
						for im in xrange(N_start, N_stop):
							# ctf2
							for i in xrange(len_ctm):
								Cls_ctf2[int(assign[im])][i] += ctf2[im][i]
							# ave
							CTFxF = filt_table(im_M[im], ctf[im])
							Util.add_img(Cls['ave'][int(assign[im])], CTFxF)

						# [sync] waiting the result
						mpi_barrier(MPI_COMM_WORLD)

						# [all] compute global sum, broadcast the results and obtain the average ave = S CTF.F / S CTF**2
						for k in xrange(K):
							Cls_ctf2[k] = mpi_reduce(Cls_ctf2[k], len_ctm, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
							Cls_ctf2[k] = mpi_bcast(Cls_ctf2[k], len_ctm, MPI_FLOAT, main_node, MPI_COMM_WORLD)
							Cls_ctf2[k] = Cls_ctf2[k].tolist() # convert array gave by MPI to list

							reduce_EMData_to_root(Cls['ave'][k], myid, main_node)
							bcast_EMData_to_all(Cls['ave'][k], myid, main_node)

							valCTF = [0] * len_ctm
							for i in xrange(len_ctm):	valCTF[i] = 1.0 / float(Cls_ctf2[k][i])
							Cls['ave'][k] = filt_table(Cls['ave'][k], valCTF)
						
					else:			
						# [id] Update clusters averages
						for im in xrange(N_start, N_stop): Util.add_img(Cls['ave'][int(assign[im])], im_M[im])

						# [sync] waiting the result
						mpi_barrier(MPI_COMM_WORLD)

						# [all] compute global sum, broadcast the results and obtain the average
						for k in xrange(K):
							reduce_EMData_to_root(Cls['ave'][k], myid, main_node) 
							bcast_EMData_to_all(Cls['ave'][k], myid, main_node)
							Cls['ave'][k] = Util.mult_scalar(Cls['ave'][k], 1.0/float(Cls['n'][k]))

					# [all] waiting the result
					mpi_barrier(MPI_COMM_WORLD)
				'''		
														
			# [sync]
			mpi_barrier(MPI_COMM_WORLD)
			
			# [id] compute the number of objects
			for k in xrange(K): 		  Cls['n'][k] = 0
			for n in xrange(N_start, N_stop): Cls['n'][int(assign[n])] += 1			

			# [all] sum the number of objects in each node and broadcast
			Cls['n'] = mpi_reduce(Cls['n'], K, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
			Cls['n'] = mpi_bcast(Cls['n'], K, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			Cls['n'] = Cls['n'].tolist() # convert array gave by MPI to list
			
			# [all] init average, Ji and ctf2
			for k in xrange(K):
				if Cls['n'][k] < 2:
					ERROR('Empty groups in kmeans_classical', 'k_means_SSE_MPI', 1)
					exit()	
				
				Cls['ave'][k].to_zero()
				Cls['Ji'][k] = 0
				if CTF:	Cls_ctf2[k] = [0] * len_ctm	
			
			if CTF:
				# [id] compute local S ctf2 and local S ave	
				for im in xrange(N_start, N_stop):
					# ctf2
					for i in xrange(len_ctm):
						Cls_ctf2[int(assign[im])][i] += ctf2[im][i]
					# ave
					CTFxF = filt_table(im_M[im], ctf[im])
					Util.add_img(Cls['ave'][int(assign[im])], CTFxF)
				
				# [sync] waiting the result
				mpi_barrier(MPI_COMM_WORLD)
				
				# [all] compute global sum, broadcast the results and obtain the average ave = S CTF.F / S CTF**2
				for k in xrange(K):
					Cls_ctf2[k] = mpi_reduce(Cls_ctf2[k], len_ctm, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
					Cls_ctf2[k] = mpi_bcast(Cls_ctf2[k], len_ctm, MPI_FLOAT, main_node, MPI_COMM_WORLD)
					Cls_ctf2[k] = Cls_ctf2[k].tolist() # convert array gave by MPI to list

					reduce_EMData_to_root(Cls['ave'][k], myid, main_node)
					bcast_EMData_to_all(Cls['ave'][k], myid, main_node)
					
					valCTF = [0] * len_ctm
					for i in xrange(len_ctm):	valCTF[i] = 1.0 / float(Cls_ctf2[k][i])
					Cls['ave'][k] = filt_table(Cls['ave'][k], valCTF)

				# [id] compute Ji
				for im in xrange(N_start, N_stop):
					CTFxAve = filt_table(Cls['ave'][int(assign[im])], ctf[im])
					Cls['Ji'][int(assign[im])] += CTFxAve.cmp("SqEuclidean", im_M[im]) / norm

				# [all] waiting the result
				mpi_barrier(MPI_COMM_WORLD)

				# [all] global sum Ji
				Cls['Ji'] = mpi_reduce(Cls['Ji'], K, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
				Cls['Ji'] = mpi_bcast(Cls['Ji'],  K, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				Cls['Ji'] = Cls['Ji'].tolist()

				
			else:			
				# [id] Update clusters averages
				for im in xrange(N_start, N_stop): Util.add_img(Cls['ave'][int(assign[im])], im_M[im])

				# [sync] waiting the result
				mpi_barrier(MPI_COMM_WORLD)

				# [all] compute global sum, broadcast the results and obtain the average
				for k in xrange(K):
					reduce_EMData_to_root(Cls['ave'][k], myid, main_node) 
					bcast_EMData_to_all(Cls['ave'][k], myid, main_node)
					Cls['ave'][k] = Util.mult_scalar(Cls['ave'][k], 1.0/float(Cls['n'][k]))

				# [id] compute Ji
				for im in xrange(N_start, N_stop): Cls['Ji'][int(assign[im])] += im_M[im].cmp("SqEuclidean", Cls['ave'][int(assign[im])])/norm

			# [all] compute Je
			Je = 0
			for k in xrange(K): Je += Cls['Ji'][k]

			# [all] waiting the result
			mpi_barrier(MPI_COMM_WORLD)

			# [all] calculate Je global sum and broadcast
			Je = mpi_reduce(Je, 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
			Je = mpi_bcast(Je, 1, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			Je = Je.tolist()
			Je = Je[0]
			
			if SA:
				# Simulate annealing, update temperature
				T *= F
				if SA2:
					if T < 0.09: SA = False
				
				#[id] informations display
				if myid == main_node: print_msg('> iteration: %5d    criterion: %11.6e   T: %13.8f  disturb:  %5d\n' % (ite, Je, T, ct_pert))
			else:
				#[id] informations display
				if myid == main_node: print_msg('> iteration: %5d    criterion: %11.6e\n' % (ite, Je))
			
			# Convergence control: threshold on the criterion value
			if Je != 0: thd = abs(Je - old_Je) / Je
			else:       thd = 0
			if SA:
				if thd < 1.0e-12 and ct_pert == 0: change = 0
			else:
				if thd < 1.0e-8: change = 0

			# Convergence control: if Je sway, means clusters unstable, due to the parallel version of k-means
			# store Je_n, Je_(n-1), Je_(n-2)
			sway_Je[1] = sway_Je[0]
			sway_Je[0] = old_Je
			old_Je = Je

			# count the consecutive number of Je periods egal at 2 (if Je_n == Je_(n-2))
			if Je == sway_Je[1]: sway_ct += 1
			else:                sway_ct == 0
			# if more of 4 periods consecutive, stop the program, detection of clusters unstable
			if sway_ct == 4 and not SA:
				change = 0
				if myid == main_node: print_msg('> Stop iteration, unstable cluster detected: Criterion sway.\n')
			
			# [all] Need to broadcast this value because all node must run together
			change = mpi_reduce(change, 1, MPI_INT, MPI_LOR, main_node, MPI_COMM_WORLD)
			change = mpi_bcast(change, 1, MPI_INT, main_node, MPI_COMM_WORLD)
			change = change.tolist()
			change = change[0]
			
		# [all] waiting the result
		mpi_barrier(MPI_COMM_WORLD)
				
		# [id] memorize the result for this trial	
		if trials > 1:
			MemCls[ntrials-1]    = deepcopy(Cls)
			MemJe[ntrials-1]     = deepcopy(Je)
			MemAssign[ntrials-1] = deepcopy(assign)
			
	# if severals trials choose the best
	if trials > 1:
		val_min = 1.0e20
		best    = -1
		for n in xrange(trials):
			if MemJe[n] < val_min:
				val_min = MemJe[n]
				best    = n
		# affect the best
		Cls    = MemCls[best]
		Je     = MemJe[best]
		assign = MemAssign[best]
	
	
	if CTF:
		# [id] the variance S (F - CTFxAve)**2
		for im in xrange(N_start, N_stop):
			CTFxAve = filt_table(Cls['ave'][int(assign[im])], ctf[im])
			buf.to_zero()
			buf = Util.subn_img(CTFxAve, im_M[im])
			Util.add_img(Cls['var'][int(assign[im])], buf) # **2
		
		# [all] waiting the result
		mpi_barrier(MPI_COMM_WORLD)
		
		# [all] global sum var
		for k in xrange(K): reduce_EMData_to_root(Cls['var'][k], myid, main_node)
		
		# compute criterion
		crit = k_means_criterion(Cls, assign, crit_name)
		
	else:
		# [id] compute the variance 1/n S(im-ave)**2 -> 1/n (Sim**2 - n ave**2)	
		for im in xrange(N_start, N_stop): Util.add_img2(Cls['var'][int(assign[im])], im_M[im])
		
		# [all] waiting the result
		mpi_barrier(MPI_COMM_WORLD)

		# [all] global sum var
		for k in xrange(K): reduce_EMData_to_root(Cls['var'][k], myid, main_node)	
		
		# [main] caclculate the variance for each cluster
		if myid == main_node:
			for k in xrange(K):
				buf.to_zero()
				Util.add_img2(buf, Cls['ave'][k])
				Cls['var'][k] = Util.madn_scalar(Cls['var'][k], buf, -float(Cls['n'][k]))
				Util.mul_scalar(Cls['var'][k], 1.0/float(Cls['n'][k]))
				
				# Uncompress ave and var images if the mask is used
				if mask != None:
					Cls['ave'][k] = Util.reconstitute_image_mask(Cls['ave'][k], mask)
					Cls['var'][k] = Util.reconstitute_image_mask(Cls['var'][k], mask)
	
	# [id] prepare assign to update
	for n in xrange(N_start):
		assign[n] = 0
	for n in xrange(N_stop, N):
		assign[n] = 0
		
	# [all] gather in main_node
	assign = mpi_reduce(assign, N, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
	assign = assign.tolist() # convert array gave by MPI to list

	"""
	# [main] Watch dog
	if myid == main_node:
		import pickle
		f = open('Assign', 'w')
		pickle.dump(assign, f)
		f.close()
	"""
	# compute Ji global
	for k in xrange(K): Cls['Ji'][k] = mpi_reduce(Cls['Ji'][k], 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
	if myid == main_node:
		for k in xrange(K): Cls['Ji'][k] = Cls['Ji'][k].tolist()[0]

	# [main_node] write the result
	if myid == main_node and CTF:
		# ifft
		for k in xrange(K):
			fftip(Cls['ave'][k])
			fftip(Cls['var'][k])
			Cls['ave'][k].postift_depad_corner_inplace()
			Cls['var'][k].postift_depad_corner_inplace()

	# [main_node] information display
	if myid == main_node:
		running_time(t_start)
		print_msg('Criterion = %11.6e \n' % Je)
		for k in xrange(K): print_msg('Cls[%i]: %i\n' % (k, Cls['n'][k]))

	# [all] waiting all nodes
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node: return Cls, assign
	else:                 return None, None

# to figure out the number of clusters
def k_means_groups_serial(stack, out_file, maskname, opt_method, K1, K2, rand_seed, maxit, trials, crit_name, CTF, F, T0, SA2, DEBUG=False):
	from utilities   import print_begin_msg, print_end_msg, print_msg, running_time
	from statistics  import k_means_open_im, k_means_criterion, k_means_headlog
	from statistics  import k_means_classical, k_means_SSE
	import time
	import string
	import os
	import sys

	if os.path.exists(out_file):  os.system('rm -rf ' + out_file)
	os.mkdir(out_file)

	if stack.split(':')[0] == 'bdb': BDB = True
	else:                            BDB = False
	
	t_start = time.time()
	N       = EMUtil.get_image_count(stack)								

	print_begin_msg('k-means groups')
	k_means_headlog(stack, '', opt_method, N, [K1, K2], crit_name, maskname, trials, maxit, CTF, T0, F, SA2, rand_seed, 1)
	
	[im_M, mask, ctf, ctf2] = k_means_open_im(stack, maskname, 0, N, N, CTF, BDB)
	
	# watch file
	out = open(out_file + '/WATCH_GRP_KMEANS', 'w')
	out.write('Watch grp k-means ' + time.ctime() +'\n')
	out.close()
	
	Crit = {}
	KK   = range(K1, K2 + 1)	# Range of works
	if crit_name == 'all': crit_name='CHD'
	
	# init the criterion
	for name in crit_name:
		if name == 'C':		
			Crit['C']=[0] * len(KK)
			txt_C = ''
		elif name == 'H':
			Crit['H']=[0] * len(KK)
			txt_H = ''
		elif name == 'D':
			Crit['D']=[0] * len(KK)
			txt_DB = ''
		else:
			ERROR('Kind of criterion k-means unknown', 'k_means_groups', 1)
			sys.exit()
	
	# init the file result
	file_crit = open(out_file + '/' + out_file, 'w')
	file_crit.write('# Criterion of k-means group\n')
	txt = '# N  '
	for name in crit_name:
		if name == 'C':
			txt += '  Coleman'.ljust(13, ' ')
		elif name == 'H':
			txt += '  Harabasz'.ljust(13, ' ')
		elif name == 'D':
			txt += '  Davies-Bouldin'.ljust(13, ' ')

	txt += '\n'
	file_crit.write(txt)
		
	# Compute the criterion and format
	index = -1
	stop  = False 
	for K in KK:
		
		index += 1
		print_msg('\n')
		print_msg('| K=%d |====================================================================\n' % K)

		flag_run = True
		ct_rnd   = 0
		while flag_run:
			try:
				if opt_method   == 'cla': [Cls, assign] = k_means_classical(im_M, mask, K, rand_seed, maxit, trials, [CTF, ctf, ctf2], F, T0, SA2, DEBUG)
				elif opt_method == 'SSE': [Cls, assign] = k_means_SSE(im_M, mask, K, rand_seed, maxit, trials, [CTF, ctf, ctf2], F, T0, SA2, DEBUG)
				else:			  ERROR('Kind of k-means unknown', 'k_means_groups', 1)
				flag_run = False
			except SystemExit:
				rand_seed += 100
				ct_rnd    += 1
				print_msg('Empty cluster, restart with random seed: %d\n' % rand_seed)
				if ct_rnd >= 5:
					stop = True
					break
				else: flag_run = True

		if stop:
			print_msg('\n=== STOP number of groups to high. ===\n\n')
			break
		
		crit = k_means_criterion(Cls, crit_name)

		# watch file
		out = open(out_file + '/WATCH_GRP_KMEANS', 'a')
		out.write('%3d  C: %11.4e  H: %11.4e  DB: %11.4e | %s\n' % (K, crit['C'], crit['H'], crit['D'], time.ctime())) 
		out.close()
		
		# result file
		txt = '%3d ' % K
		ipos, pos = 2, {}
		for name in crit_name:
			if name == 'C':
				txt += '  %11.4e' % crit['C']
				pos['C'] = ipos
				Crit['C'][index]  = crit['C']
			elif name == 'H':
				txt += '  %11.4e' % crit['H']
				pos['H'] = ipos
				Crit['H'][index]  = crit['H']
			elif name == 'D':
				txt += '  %11.4e' % crit['D']
				pos['D'] = ipos
				Crit['D'][index]  = crit['D']
			ipos += 1
		txt += '\n'
		file_crit.write(txt)
		
	file_crit.close()	
	
	# make script file to gnuplot
	out = open(out_file + '/' + out_file + '.p', 'w')
	out.write('# Gnuplot script file for plotting result in kmeans groups\n')
	out.write('# cmd: gnuplot> load \x22%s.p\x22\n' % out_file)
	out.write('reset\n')
	out.write('set autoscale\n')
	txt = 'plot'
	
	for name in crit_name:
		if name == 'C':
			# norm plot [0;1]
			d = max(Crit['C']) - min(Crit['C'])
			a = 1 / d
			b = 0.5 - ((max(Crit['C']) * a + min(Crit['C']) * a) / 2)
			txt += ' \x22%s\x22 using 1:($%d*(%11.4e)+(%11.4e)) title \x22Coleman\x22 w l,' % (out_file, pos['C'], a, b)
		elif name == 'H':
			# norm plot [0;1]
			d = max(Crit['H']) - min(Crit['H'])
			a = 1 / d
			b = 0.5 - ((max(Crit['H']) * a + min(Crit['H']) * a) / 2)
			txt += ' \x22%s\x22 using 1:($%d*(%11.4e)+(%11.4e)) title \x22Harabasz\x22 w l,' % (out_file, pos['H'], a, b)
		elif name == 'D':
			# norm plot [0;1]
			d = max(Crit['D']) - min(Crit['D'])
			a = 1 / d
			b = 0.5 - ((max(Crit['D']) * a + min(Crit['D']) * a) / 2)
			txt += ' \x22%s\x22 using 1:($%d*(%11.4e)+(%11.4e)) title \x22Davies-Bouldin\x22 w l,' % (out_file, pos['D'], a, b)
	txt = txt.rstrip(',') + '\n'
	out.write(txt)
	out.close()
		
	# runtime
	running_time(t_start)
	print_end_msg('k-means groups')
		
	return Crit, KK

# to figure out the number of clusters MPI version
def k_means_groups_MPI(stack, out_file, maskname, opt_method, K1, K2, rand_seed, maxit, trials, crit_name, CTF, F, T0, SA2):
	from utilities   import print_begin_msg, print_end_msg, print_msg, running_time
	from mpi 	 import MPI_COMM_WORLD, mpi_barrier
	from statistics  import k_means_init_MPI, k_means_open_im, k_means_criterion, k_means_headlog
	from statistics  import k_means_cla_MPI, k_means_SSE_MPI
	import sys
	import time
	import os
	
	# [id]   part of different code for each node
	# [sync] synchronise each node
	# [main] part of code just for the main node
	# [all]  code write for all node
	N = EMUtil.get_image_count(stack)
	main_node, myid, ncpu, N_start, N_stop = k_means_init_MPI(N)

	if myid == main_node:
		if os.path.exists(out_file): os.system('rm -rf ' + out_file)
		os.mkdir(out_file)

	if stack.split(':')[0] == 'bdb': BDB = True
	else:                            BDB = False

	if myid == main_node:
		t_start = time.time()

		print_begin_msg('k-means groups')
		k_means_headlog(stack, '', opt_method, N, [K1, K2], crit_name, maskname, trials, maxit, CTF, T0, F, SA2, rand_seed, ncpu)
		
		# watch file
		out = open(out_file + '/WATCH_GRP_KMEANS', 'w')
		out.write('Watch grp k-means ' + time.ctime() +'\n')
		out.close()

	[im_M, mask, ctf, ctf2] = k_means_open_im(stack, maskname, N_start, N_stop, N, CTF, BDB)

	# define the node affected to the numbers tests
	Crit = {}
	KK = range(K1, K2 + 1)	# Range of works
	if crit_name == 'all': crit_name='CHD'
	
	#[all] init the criterion
	# init the criterion
	for name in crit_name:
		if name == 'C':		
			Crit['C']=[0] * len(KK)
			txt_C = ''
		elif name == 'H':
			Crit['H']=[0] * len(KK)
			txt_H = ''
		elif name == 'D':
			Crit['D']=[0] * len(KK)
			txt_DB = ''
		else:	ERROR('Kind of criterion k-means unknown', 'k_means_groups', 1)
	
	# init the file result
	if myid == main_node:
		file_crit = open(out_file + '/' + out_file, 'w')
		file_crit.write('# Criterion of k-means group\n')
		txt = '# N  '
		for name in crit_name:
			if name == 'C':
				txt += '  Coleman'.ljust(13, ' ')
			elif name == 'H':
				txt += '  Harabasz'.ljust(13, ' ')
			elif name == 'D':
				txt += '  Davies-Bouldin'.ljust(13, ' ')
			else:
				ERROR('Kind of criterion k-means unknown', 'k_means_groups', 1)
		txt += '\n'
		file_crit.write(txt)
	
	index = -1

	stop  = False
	# [id] Measure the criterions
	for K in KK:

		index += 1
		if myid == main_node:
			print_msg('\n')
			print_msg('| K=%d |====================================================================\n' % K)

		flag_run = True
		ct_rnd   = 0
		while flag_run:
			try:
				if   opt_method == 'cla':
					[Cls, assign] = k_means_cla_MPI(im_M, mask, K, rand_seed, maxit, trials, [CTF, ctf, ctf2], myid, main_node, N_start, N_stop, F, T0, SA2)
				elif opt_method == 'SSE':
					[Cls, assign] = k_means_SSE_MPI(im_M, mask, K, rand_seed, maxit, trials, [CTF, ctf, ctf2], myid, main_node, ncpu, N_start, N_stop, F, T0, SA2)
				else:   ERROR('kind of k-means unknown', 'k_means_groups', 1)
				flag_run = False
			except SystemExit:
				rand_seed += 100
				ct_rnd    += 1
				if myid == main_node: print_msg('Empty cluster, restart with random seed: %d\n' % rand_seed)
				if ct_rnd >= 5:
					stop = True
					break
				else: flag_run = True

		if stop:
			if myid == main_node: print_msg('\n=== STOP number of groups to high. ===\n\n')
			break
					

		if myid == main_node: crit = k_means_criterion(Cls, crit_name)

		mpi_barrier(MPI_COMM_WORLD)

		# [main] Format and export the result
		if myid == main_node:
			# watch file
			out = open(out_file + '/WATCH_GRP_KMEANS', 'a')
			out.write('%d  C: %11.4e  H: %11.4e  DB: %11.4e | %s\n' % (K, crit['C'], crit['H'], crit['D'], time.ctime())) 
			out.close()

			# result file
			txt = '%3d ' % K
			ipos, pos = 2, {}
			for name in crit_name:
				if name == 'C':
					txt += '  %11.4e' % crit['C']
					pos['C'] = ipos
					Crit['C'][index]  = crit['C']	# Coleman
				elif name == 'H':
					txt += '  %11.4e' % crit['H']
					pos['H'] = ipos
					Crit['H'][index]  = crit['H']	# Harabasz
				elif name == 'D':
					txt += '  %11.4e' % crit['D']
					pos['D'] = ipos
					Crit['D'][index]  = crit['D']	# Davies & Bouldin
				ipos += 1
			txt += '\n'
			file_crit.write(txt)

	if myid == main_node:
		# make script file to gnuplot
		out = open(out_file + '/' + out_file + '.p', 'w')
		out.write('# Gnuplot script file for plotting result in kmeans groups\n')
		out.write('# cmd: gnuplot> load \x22%s.p\x22\n' % out_file)
		out.write('reset\n')
		out.write('set autoscale\n')
		txt = 'plot'

		for name in crit_name:
			if name == 'C':
				# norm plot [0;1]
				d = max(Crit['C']) - min(Crit['C'])
				a = 1 / d
				b = 0.5 - ((max(Crit['C']) * a + min(Crit['C']) * a) / 2)
				txt += ' \x22%s\x22 using 1:($%d*(%11.4e)+(%11.4e)) title \x22Coleman\x22 w l,' % (out_file, pos['C'], a, b)
			elif name == 'H':
				# norm plot [0;1]
				d = max(Crit['H']) - min(Crit['H'])
				a = 1 / d
				b = 0.5 - ((max(Crit['H']) * a + min(Crit['H']) * a) / 2)
				txt += ' \x22%s\x22 using 1:($%d*(%11.4e)+(%11.4e)) title \x22Harabasz\x22 w l,' % (out_file, pos['H'], a, b)
			elif name == 'D':
				# norm plot [0;1]
				d = max(Crit['D']) - min(Crit['D'])
				a = 1 / d
				b = 0.5 - ((max(Crit['D']) * a + min(Crit['D']) * a) / 2)
				txt += ' \x22%s\x22 using 1:($%d*(%11.4e)+(%11.4e)) title \x22Davies-Bouldin\x22 w l,' % (out_file, pos['D'], a, b)

		txt = txt.rstrip(',') + '\n'
		out.write(txt)
		out.close()
		
		# logfile
		file_crit.close()
		running_time(t_start)
		print_end_msg('k-means groups')

	if myid == main_node: return Crit, KK
	else:                 return None, None

### END K-MEANS ##############################################################################
##############################################################################################

# helper functions for ali2d_ra and ali2d_rac
def kmn(data, numr, wr, cm = 0, max_iter = 10, this_seed = 1000):
	from utilities    import model_blank, print_msg
	from random       import seed, randint, shuffle
	seed(this_seed)
	if (not cm):  mirror = 0
	nima = len(data)
	# Randomization/shuffling
	rnd = range(nima)
	shuffle(rnd)
	# calculate random approximation of the total average
	tave = data[rnd[0]].copy()
	for imi in xrange(1, nima):
		# align current image to the reference minus this image
		im = rnd[imi]
		if (cm):
			retvals = Util.Crosrng_ew(tave, data[im], numr, wr, 0)
			qn  = retvals["qn"]
			tot = retvals["tot"]
			retvals = Util.Crosrng_ew(tave, data[im], numr, wr, 1)
			qm  = retvals["qn"]
			tmt = retvals["tot"]
			if (qn >= qm):
				alpha  = tot
				mirror = 0
			else:
				alpha  = tmt
				mirror = 1
		else:
			retvals = Util.Crosrng_ew(tave, data[im], numr, wr, 0)
			alpha = retvals["tot"]
		data[im].set_attr_dict({'alpha':alpha, 'mirror': mirror})
		Util.update_fav(tave, data[im], alpha, mirror, numr)
	
	nx = tave.get_xsize()
	a0 = Util.ener(tave, numr)
	msg = "Initial criterion = %20.7e\n"%(a0)
	print_msg(msg)

	# do the alignment
	for Iter in xrange(max_iter):
		again = False
		for im in xrange(nima):
			it = randint(0, nima-1)
			tmp = rnd[im]; rnd[im] = rnd[it]; rnd[it] = tmp;
		for imi in xrange(nima):
			# subtract current image from the average
			im = rnd[imi]
			alpha  = data[im].get_attr('alpha')
			mirror =  data[im].get_attr('mirror')
			# Subtract current image from the average
			refim = tave.copy()
			Util.sub_fav(refim, data[im], alpha, mirror, numr)
			# align current image to the reference minus this image
			if(cm):
				retvals = Util.Crosrng_ew(refim, data[im], numr, wr, 0)
				qn  = retvals["qn"]
				tot = retvals["tot"]
				retvals = Util.Crosrng_ew(refim, data[im], numr, wr, 1)
				qm  = retvals["qn"]
				tmt = retvals["tot"]
		   		if (qn >= qm):
					alpha  = tot
					mirror = 0
		   		else:
					alpha  = tmt
					mirror = 1
			else:
				retvals = Util.Crosrng_ew(refim, data[im], numr, wr, 0)
				alpha = retvals["tot"]
			Util.update_fav(refim, data[im], alpha, mirror, numr)
			# calculate the criterion
			a1 = Util.ener(refim, numr)
			if(a1 > a0):
				# replace the average by the improved average and set the new parameters to the image, otherwise, do nothing
				tave = refim.copy()
				data[im].set_attr_dict({'alpha':alpha, 'mirror': mirror})
				a0 = a1
				again = True
		if (again):
			# calculate total average using current alignment parameters
			tave = model_blank(nx)
			for im in xrange(nima):  
				alpha  = data[im].get_attr('alpha')
				mirror = data[im].get_attr('mirror')
				Util.update_fav(tave, data[im], alpha, mirror, numr)
			a1 = Util.ener(tave, numr)
			msg = "ITERATION #%3d        criterion = %20.7e\n"%(Iter+1,a1)
			print_msg(msg)
			if (a1 <= a0):  break
			else:          a0 = a1
		else:  break


def kmn_a(data, numr, wr, cm = 0, max_iter = 10, this_seed = 1000):
	from utilities    import model_blank, print_msg, amoeba
	from random       import seed, randint, shuffle
	seed(this_seed)
	if (not cm):  mirror = 0
	nima = len(data)
	# Randomization/shuffling
	rnd = range(nima)
	shuffle(rnd)
	# calculate random approximation of the total average
	tave = data[rnd[0]].copy()
	for imi in xrange(1, nima):
		# align current image to the reference minus this image
		im = rnd[imi]
		if (cm):
			retvals = Util.Crosrng_ew(tave, data[im], numr, wr, 0)
			qn  = retvals["qn"]
			tot = retvals["tot"]
			retvals = Util.Crosrng_ew(tave, data[im], numr, wr, 1)
			qm  = retvals["qn"]
			tmt = retvals["tot"]
			if (qn >= qm):
				alpha  = tot
				mirror = 0
			else:
				alpha  = tmt
				mirror = 1
		else:
			retvals = Util.Crosrng_ew(tave, data[im], numr, wr, 0)
			alpha = retvals["tot"]
		data[im].set_attr_dict({'alpha':alpha, 'mirror': mirror})
		Util.update_fav(tave, data[im], alpha, mirror, numr)
	
	nx = tave.get_xsize()
	a0 = Util.ener(tave, numr)
	msg = "Initial criterion = %-20.7e\n"%(a0)
	print_msg(msg)

	# do the alignment
	for Iter in xrange(max_iter):
		again = False
		for im in xrange(nima):
			it = randint(0, nima-1)
			tmp = rnd[im]; rnd[im] = rnd[it]; rnd[it] = tmp;
		for imi in xrange(nima):
			# subtract current image from the average
			im = rnd[imi]
			alpha  = data[im].get_attr('alpha')
			mirror =  data[im].get_attr('mirror')
			# Subtract current image from the average
			refim = tave.copy()
			Util.sub_fav(refim, data[im], alpha, mirror, numr)
			# align current image to the reference minus this image
			if(cm):
				retvals = Util.Crosrng_ew(refim, data[im], numr, wr, 0)
				qn  = retvals["qn"]
				tot = retvals["tot"]
				retvals = Util.Crosrng_ew(refim, data[im], numr, wr, 1)
				qm  = retvals["qn"]
				tmt = retvals["tot"]
		   		if (qn >= qm):
					alpha  = tot
					mirror = 0
		   		else:
					alpha  = tmt
					mirror = 1
			else:
				retvals = Util.Crosrng_ew(refim, data[im], numr, wr, 0)
				alpha = retvals["tot"]
			Util.update_fav(refim, data[im], alpha, mirror, numr)
			# calculate the criterion
			a1 = Util.ener(refim, numr)
			if(a1 > a0):
				# replace the average by the improved average and set the new parameters to the image, otherwise, do nothing
				tave = refim.copy()
				data[im].set_attr_dict({'alpha':alpha, 'mirror': mirror})
				a0 = a1
				again = True
		if (again):
			# calculate total average using current alignment parameters
			tave = model_blank(nx)
			for im in xrange(nima):  
				alpha  = data[im].get_attr('alpha')
				mirror = data[im].get_attr('mirror')
				Util.update_fav(tave, data[im], alpha, mirror, numr)
			a1 = Util.ener(tave, numr)
			msg = "ITERATION #%3d        criterion = %20.7e\n"%(Iter+1,a1)
			print_msg(msg)
			if (a1 <= a0):  break
			else:          a0 = a1
		else:  break
	amoeba_data = []
	amoeba_data.append(nx)
	amoeba_data.append(data)
	amoeba_data.append(numr)
	amoeba_data.append(nima)
	all_alpha = []
	for im in xrange(nima):
		alpha = data[im].get_attr('alpha')
		all_alpha.append(alpha)
	alpha_range = [2.0]*nima
	ps = amoeba(all_alpha, alpha_range, multi_search_func, 1.e-4, 1.e-4, 1000, amoeba_data)
	for im in xrange(nima):
		data[im].set_attr_dict({'alpha': ps[0][im]})
	msg = "Final criterion = %20.7e\n"%(ps[1])
	print_msg(msg)


def kmn_g(data, numr, wr, stack, check_mirror = False, max_iter = 10, this_seed = 1000):

	from utilities    import   model_blank, print_msg, amoeba, combine_params2
	from random       import   seed, randint
	from alignment    import   Applyws, ang_n
	from development  import   oned_search_func
	from random       import   gauss
	
	mode = "F"
	seed(this_seed)
	if (not check_mirror):  mirror = 0
	nima = len(data)
	# Randomization/shuffling
	rnd = range(nima)
	for im in xrange(nima):
		it = randint(0,nima-1)
		tmp = rnd[im]; rnd[im] = rnd[it]; rnd[it] = tmp;
	# calculate random approximation of the total average
	tave = data[rnd[0]].copy()
	tave_w = tave.copy()
	Applyws(tave_w, numr, wr)
	maxrin = numr[len(numr)-1]
		
	line = EMData()
	line.set_size(maxrin,1,1)
	M=maxrin
	# do not pad
	npad=1
	N=M*npad
	# support of the window
	K=6
	alpha=1.75
	r=M/2
	v=K/2.0/N
	kbline = Util.KaiserBessel(alpha, K, r, K/(2.*N), N)
	parline = {"filter_type" : Processor.fourier_filter_types.KAISER_SINH_INVERSE,"alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
	amoeba_data = []
	amoeba_data.append(kbline)
		
	for imi in xrange(1, nima):
		# align current image to the reference minus this image
		im = rnd[imi]
		if (check_mirror):
			qt = Util.Crosrng_msg(tave_w, data[im], numr)
				
			# straight
			for i in xrange(0,maxrin): line[i]=qt[i,0]					
			#  find the current maximum and its location
			ps=line.peak_search(1,1)				
			qn=ps[1]
			jtot=ps[2]/2
			q=Processor.EMFourierFilter(line,parline)
			amoeba_data.insert(0,q)
			ps = amoeba([jtot], [2.0], oned_search_func, 1.e-4, 1.e-4, 500, amoeba_data)
			del amoeba_data[0]
			jtot=ps[0][0]*2
			qn=ps[1]

			# mirror
			for i in xrange(0,maxrin): line[i]=qt[i,1]
			#  find the current maximum and its location
			ps=line.peak_search(1,1)				
			qm=ps[1]
			mtot=ps[2]/2
			q=Processor.EMFourierFilter(line,parline)
			amoeba_data.insert(0,q)
			ps = amoeba([mtot], [2.0], oned_search_func, 1.e-4, 1.e-4, 500, amoeba_data)
			del amoeba_data[0]
			mtot=ps[0][0]*2
			qm=ps[1]

			if (qn >= qm):
				#alpha = ang_n(jtot+1, mode, maxrin)
				alpha = jtot+1
				mirror = 0
			else:
				#alpha = ang_n(mtot+1, mode, maxrin)
				alpha = mtot+1
				mirror = 1
		else:		
			line_s = Util.Crosrng_msg_s(tave_w, data[im], numr)
			
			# straight
			#  find the current maximum and its location
			ps=line_s.peak_search(1,1)
			qn=ps[1]
			jtot=ps[2]/2
			q=Processor.EMFourierFilter(line_s,parline)
			amoeba_data.insert(0,q)
			ps = amoeba([jtot], [2.0], oned_search_func, 1.e-4, 1.e-4, 500, amoeba_data)
			del amoeba_data[0]
			jtot=ps[0][0]*2
			qn=ps[1]

			#alpha = ang_n(jtot+1, mode, maxrin)
			alpha = jtot+1
			mirror = 0
			
		data[im].set_attr_dict({'alpha':alpha, 'mirror': mirror})		
		Util.update_fav(tave, data[im], alpha, mirror, numr)
		tave_w = tave.copy()
		Applyws(tave_w, numr, wr)
	
	nx = tave.get_xsize()
	a0 = Util.ener(tave, numr)
	msg = "Initial criterion = %-20.7e\n"%(a0)
	print_msg(msg)
	
	# do the alignment
	for Iter in xrange(max_iter):
		again = False
		for im in xrange(nima):
			it = randint(0, nima-1)
			tmp = rnd[im]; rnd[im] = rnd[it]; rnd[it] = tmp;
		for imi in xrange(nima):
			# subtract current image from the average
			im = rnd[imi]
			alpha  = data[im].get_attr('alpha')
			mirror = data[im].get_attr('mirror')
			# Subtract current image from the average
			refim = tave.copy()
			Util.sub_fav(refim, data[im], alpha, mirror, numr)
			refim_w = refim.copy()
			Applyws(refim_w, numr, wr)
			# align current image to the reference minus this image
			if (check_mirror):
				qt = Util.Crosrng_msg(refim_w, data[im], numr)
					
				# straight
				for i in xrange(0,maxrin): line[i]=qt[i,0]					
				#  find the current maximum and its location
				ps=line.peak_search(1,1)				
				qn=ps[1]
				jtot=ps[2]/2
				q=Processor.EMFourierFilter(line,parline)
				amoeba_data.insert(0,q)
				ps = amoeba([jtot], [2.0], oned_search_func, 1.e-4, 1.e-4, 500, amoeba_data)
				del amoeba_data[0]
				jtot=ps[0][0]*2
				qn=ps[1]

				# mirror
				for i in xrange(0,maxrin): line[i]=qt[i,1]
				#  find the current maximum and its location
				ps=line.peak_search(1,1)				
				qm=ps[1]
				mtot=ps[2]/2
				q=Processor.EMFourierFilter(line,parline)
				amoeba_data.insert(0,q)
				ps = amoeba([mtot], [2.0], oned_search_func, 1.e-4, 1.e-4, 500, amoeba_data)
				del amoeba_data[0]
				mtot=ps[0][0]*2
				qm=ps[1]

				if (qn >= qm):
					#alpha = ang_n(jtot+1, mode, maxrin)
					alpha = jtot+1
					mirror = 0
				else:
					#alpha = ang_n(mtot+1, mode, maxrin)
					alpha = mtot+1
					mirror = 1
			else:
				line_s = Util.Crosrng_msg_s(refim_w, data[im], numr)
				
				# straight
				#  find the current maximum and its location
				ps=line_s.peak_search(1,1)				
				qn=ps[1]
				jtot=ps[2]/2
				q=Processor.EMFourierFilter(line_s,parline)
				amoeba_data.insert(0,q)
				ps = amoeba([jtot], [2.0], oned_search_func, 1.e-4, 1.e-4, 500, amoeba_data)
				del amoeba_data[0]
				jtot=ps[0][0]*2
				qn=ps[1]
	
				#alpha = ang_n(jtot+1, mode, maxrin)
				alpha =  jtot+1
				mirror = 0
			Util.update_fav(refim, data[im], alpha, mirror, numr)
			# calculate the criterion
			a1 = Util.ener(refim, numr)
			if (a1 > a0):
				# replace the average by the improved average and set the new parameters to the image, otherwise, do nothing
				tave = refim.copy()
				data[im].set_attr_dict({'alpha':alpha, 'mirror': mirror})
				a0 = a1
				again = True
		if (again):
			# calculate total average using current alignment parameters
			tave = model_blank(nx)
			for im in xrange(nima):  
				alpha  = data[im].get_attr('alpha')
				mirror = data[im].get_attr('mirror')
				Util.update_fav(tave, data[im], alpha, mirror, numr)
			a1 = Util.ener(tave, numr)
			msg = "ITERATION #%3d        criterion = %20.7e\n"%(Iter+1,a1)
			print_msg(msg)
			if (a1 <= a0):  break
			else:          a0 = a1
		else:  break
	
	
	temp = EMData()
	for im in xrange(nima):
		alpha_original   = data[im].get_attr('alpha_original')
		alpha = data[im].get_attr('alpha')
		sx    =  data[im].get_attr('sx')
		sy    =  data[im].get_attr('sy')
		mirror =  data[im].get_attr('mirror')
		alpha = ang_n(alpha+1, mode, maxrin)
		#  here the original angle is irrelevant, used only to determine proper shifts
		alpha_original_n, sxn, syn, mir = combine_params2(0, -sx, -sy, 0, -alpha_original, 0,0,0)
		alphan, sxn, syn, mir           = combine_params2(0, -sxn, -syn, 0, alpha, 0,0, mirror)
		temp.read_image(stack, im, True)
		#if(CTF and data_had_ctf == 0):   temp.set_attr('ctf_applied', 0)
		temp.set_attr_dict({'alpha':alphan, 'sx':sxn, 'sy':syn, 'mirror': mir})
		temp.write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)
	
	ave1, ave2 =  ave_oe_series(stack)	
	fsc(ave1, ave2, 1.0, "fsc_before_amoeba")
	
	
	amoeba_data = []
	amoeba_data.append(data)
	# temp test
 	#amoeba_data.append(stack)
	
	amoeba_data.append(numr)
	
	# temp measure
	#amoeba_data.append(mode)
	#amoeba_data.append(maxrin)
	
	alpha_range = [1.0]*nima
	new_alpha = []

	for Iter in xrange(10):
		tave = model_blank(nx)
		all_alpha = []
		for im in xrange(nima):
			alpha = data[im].get_attr('alpha')
			if Iter==0: shake = 0.0
			else: shake = gauss(0, 0.5)
			all_alpha.append(alpha+shake)
			Util.update_fav(tave, data[im], alpha+shake, 0, numr)
		amoeba_data.append(tave)
		amoeba_data.append(all_alpha)		
		ps = amoeba(all_alpha, alpha_range, multi_search_func, 1.e-4, 1.e-4, 10000, amoeba_data)
		dummy = amoeba_data.pop()
		dummy = amoeba_data.pop()
		msg = "Trial %2d of amoeba:   criterion = %20.7e    Iteration = %4d\n"%(Iter+1, ps[1], ps[2])
		print_msg(msg)				
		if ps[1]>a0:
			a0 = ps[1]			
			new_alpha = ps[0]
		else:
			if new_alpha == []:	new_alpha = all_alpha		
				
		temp = EMData()
		for im in xrange(nima):
			alpha_original   = data[im].get_attr('alpha_original')
			alpha = ps[0][im]
			sx    =  data[im].get_attr('sx')
			sy    =  data[im].get_attr('sy')
			mirror =  data[im].get_attr('mirror')
			alpha = ang_n(alpha+1, mode, maxrin)
			#  here the original angle is irrelevant, used only to determine proper shifts
			alpha_original_n, sxn, syn, mir = combine_params2(0, -sx, -sy, 0, -alpha_original, 0,0,0)
			alphan, sxn, syn, mir           = combine_params2(0, -sxn, -syn, 0, alpha, 0,0, mirror)
			temp.read_image(stack, im, True)
			#if(CTF and data_had_ctf == 0):   temp.set_attr('ctf_applied', 0)
			temp.set_attr_dict({'alpha':alphan, 'sx':sxn, 'sy':syn, 'mirror': mir})
			temp.write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)

		ave1, ave2 =  ave_oe_series(stack)	
		fsc(ave1, ave2, 1.0, "fsc_trial_%02d"%(Iter+1))
		
	for im in xrange(nima):
		data[im].set_attr_dict({'alpha': new_alpha[im]})
		
	msg = "Final criterion = %-20.7e\n"%(a0)
	print_msg(msg)
	

def multi_search_func(args, data):
	
	from utilities import model_blank
	
	fdata = data[0]
	numr = data[1]
	tave = data[2].copy()
	ori_alpha = data[3]
	
	nima = len(fdata)

	total_change = 0
       	update_list = []
	for im in xrange(nima):
		if args[im]!=ori_alpha[im]: 
			total_change += 1 
			update_list.append(im)
	
	if total_change < nima*0.3:
		for im in update_list:
			Util.sub_fav(tave, fdata[im], ori_alpha[im], 0, numr) 
			Util.update_fav(tave, fdata[im], args[im], 0, numr)
	else:
		nx = tave.get_xsize()
		tave = model_blank(nx)
		for im in xrange(nima):
			Util.update_fav(tave, fdata[im], args[im], 0, numr)
			
	energy = Util.ener(tave, numr)	
	#energy = Util.ener_tot(fdata, numr, args)
	
	return energy

"""
def multi_search_func2(args, data):
	
	from utilities import model_blank, model_circle
	from alignment import ang_n
	
	stack = data[0]
	numr = data[1]
	mode = data[2]
	maxrin = data[3]
	
	nima = EMUtil.get_image_count(stack)	
	for im in xrange(nima):
		alpha = args[im]
		alpha = ang_n(alpha+1, mode, maxrin)
		img = EMData()
		img.read_image(stack, im)
		img.set_attr_dict({"alpha":alpha})
		img.write_image(stack, im)
	wie, var = aves_wiener(stack, "a")
	mask = model_circle(72, 150, 150)
	h = Util.infomask(var, mask, True)
	print h
	
	return -h[0]
"""

def kmn_ctf(data, ref_data, numr, wr, cm = 0, max_iter = 10, this_seed = 1000):
	from utilities    import model_blank, print_msg
	from random       import seed, randint, shuffle
	#  This version is not quite correct, as ctf2 is not modified in the avergae, but it cannot be done any other simple way :-(
	seed(this_seed)
	if(not cm):  mirror = 0
	nima = len(data)
	# Randomization/shuffling
	rnd = range(nima)
	shuffle(rnd)
	# calculate random approximation of the total average
	tave = ref_data[rnd[0]].copy()
	for imi in xrange(1, nima):
		# align current image to the reference minus this image
		im = rnd[imi]
		if(cm):
			retvals = Util.Crosrng_ew(tave, data[im], numr, wr, 0)
			qn  = retvals["qn"]
			tot = retvals["tot"]
			retvals = Util.Crosrng_ew(tave, data[im], numr, wr, 1)
			qm  = retvals["qn"]
			tmt = retvals["tot"]
			if (qn >= qm):
				alpha  = tot
				mirror = 0
			else:
				alpha  = tmt
				mirror = 1
		else:
			retvals = Util.Crosrng_ew(tave, data[im], numr, wr, 0)
			alpha = retvals["tot"]
		data[im].set_attr_dict({'alpha':alpha, 'mirror': mirror})
		Util.update_fav(tave, ref_data[im], alpha, mirror, numr)
	
	nx = tave.get_xsize()
	a0 = Util.ener(tave, numr)
	msg = "Initial criterion = %-20.7e\n"%(a0)
	print_msg(msg)

	# do the alignment
	for Iter in xrange(max_iter):
		again = False
		for im in xrange(nima):
			it = randint(0,nima-1)
			tmp = rnd[im]; rnd[im] = rnd[it]; rnd[it] = tmp;
		for imi in xrange(nima):
			# subtract current image from the average
			im = rnd[imi]
			alpha  = data[im].get_attr('alpha')
			mirror =  data[im].get_attr('mirror')
			# Subtract current image from the average
			refim = tave.copy()
			Util.sub_fav(refim, ref_data[im], alpha, mirror, numr)
			# align current image to the reference minus this image
			if(cm):
				retvals = Util.Crosrng_ew(refim, data[im], numr, wr, 0)
				qn  = retvals["qn"]
				tot = retvals["tot"]
				retvals = Util.Crosrng_ew(refim, data[im], numr, wr, 1)
				qm  = retvals["qn"]
				tmt = retvals["tot"]
		   		if (qn >= qm):
					alpha  = tot
					mirror = 0
		   		else:
					alpha  = tmt
					mirror = 1
			else:
				retvals = Util.Crosrng_ew(refim, data[im], numr, wr, 0)
				alpha = retvals["tot"]
			Util.update_fav(refim, ref_data[im], alpha, mirror, numr)
			# calculate the criterion
			a1 = Util.ener(refim, numr)
			if(a1 > a0):
				# replace the average by the improved average and set the new parameters to the image, otherwise, do nothing
				tave = refim.copy()
				data[im].set_attr_dict({'alpha':alpha, 'mirror': mirror})
				a0 = a1
				again = True
		if(again):
			# calculate total average using current alignment parameters
			tave = model_blank(nx)
			for im in xrange(nima):  
				alpha  = data[im].get_attr('alpha')
				mirror =  data[im].get_attr('mirror')
				Util.update_fav(tave, ref_data[im], alpha, mirror, numr)
			a1 = Util.ener(tave, numr)
			msg = "ITERATION #%3d        criterion = %20.7e\n"%(Iter+1,a1)
			print_msg(msg)
			if(a1 <= a0):  break
			else:          a0 = a1
		else:  break

def kmnr(data, assign, nima, k, numr, wr, cm = 0, max_iter = 10, this_seed = 1000, MPI = False):
	from random       import seed, randint, random, shuffle
	from utilities    import model_blank
	seed(this_seed)
	if(not cm):  mirror = 0
	ntot = len(data)

	if not MPI:
		# create list of objects in this group
		lob = [0]*nima
		im = -1
		for i in xrange(ntot):
			if(assign[i] == k):
				im +=1
				lob[im] = i
	else:
		# use the list of index images given (used in ali2d_rac_MPI)
		im  = nima - 1

	# random / shuffle
	rnd = range(nima)
	shuffle(rnd)
				
	# calculate random approximation of the total average
	if not MPI: tave = data[lob[rnd[0]]].copy()
	else:       tave = data[rnd[0]].copy()
		
	for imi in xrange(1, nima):
		# align current image to the reference minus this image
		if not MPI: im = lob[rnd[imi]]
		else:       im = rnd[imi]
		
		if(cm):
			retvals = Util.Crosrng_ew(tave, data[im], numr, wr, 0)
			qn  = retvals["qn"]
			tot = retvals["tot"]
			retvals = Util.Crosrng_ew(tave, data[im], numr, wr, 1)
			qm  = retvals["qn"]
			tmt = retvals["tot"]
			if (qn >= qm):
			   alpha  = tot
			   mirror = 0
			else:
			   alpha  = tmt
			   mirror = 1
		else:
			retvals = Util.Crosrng_ew(tave, data[im], numr, wr, 0)
			alpha = retvals["tot"]
		
		data[im].set_attr_dict({'alpha':alpha, 'mirror': mirror})
		Util.update_fav(tave, data[im], alpha, mirror, numr)
	
	nr = tave.get_xsize()
	at = Util.ener(tave, numr)
	a0 = at
	#print  "Initial criterion = ",a0
	temp = EMData()
	# do the alignment
	for iter in xrange(max_iter):
		again = False
		rnd = range(nima)
		shuffle(rnd)
		for imi in xrange(nima):
			# subtract current image from the average
			if not MPI: im = lob[rnd[imi]]
			else:       im = rnd[imi]
			
			alpha  = data[im].get_attr('alpha')
			mirror =  data[im].get_attr('mirror')
			#  Subtract current image from the average
			refim = tave.copy()
			Util.sub_fav(refim, data[im], alpha, mirror, numr)
			# align current image to the reference minus this image
			if(cm):
				retvals = Util.Crosrng_ew(refim, data[im], numr, wr, 0)
				qn  = retvals["qn"]
				tot = retvals["tot"]
				retvals = Util.Crosrng_ew(refim, data[im], numr, wr, 1)
				qm  = retvals["qn"]
				tmt = retvals["tot"]
		   		if (qn >= qm):
				   alpha  = tot
				   mirror = 0
		   		else:
				   alpha  = tmt
				   mirror = 1
			else:
				retvals = Util.Crosrng_ew(refim, data[im], numr, wr, 0)
				alpha = retvals["tot"]
			Util.update_fav(refim, data[im], alpha, mirror, numr)
			# calculate the criterion
			a1 = Util.ener(refim, numr)
			#print  im, a1, psi, mirror
			if(a1 > a0):
				# replace the average by the improved average and set the new parameters to the image, otherwise, do nothing
				tave = refim.copy()
				data[im].set_attr_dict({'alpha':alpha, 'mirror': mirror})
				a0 = a1
				again = True
		if(again):
			# calculate total average using current alignment parameters
			tave = model_blank(nr)
			for imi in xrange(nima):
				if not MPI: im = lob[imi]
				else:       im = rnd[imi]
				
				alpha  = data[im].get_attr('alpha')
				mirror =  data[im].get_attr('mirror')
				Util.update_fav(tave, data[im], alpha, mirror, numr)
			a0 = Util.ener(tave, numr)
			#print " ITERATION #",iter+1,"  criterion = ",a0
			if(a0 <= at):  break
			else:          at = a0
				
		else:  break

	return  tave
	
def Wiener_CSQ(data, K, assign, Cls, ctf1, ctf2, snr = 1.0):
	from filter       import filt_table

	N = len(data)
	lctf = len(ctf2[0])
	ctf_2 = [[0.0]*lctf for k in xrange(K)]

	for k in xrange(K):
		Cls[k].C.to_zero()
		Cls[k].n = 0
		Cls[k].SSE = 0.0

	#Calculate averages
	for im in xrange(N):
		k = assign[im]
		avec = filt_table( data[im], ctf1[im] )  # multiply data by the CTF
		Util.add_img(Cls[k].C, avec)
		for i in xrange(lctf): ctf_2[k][i] += ctf2[im][i]
		Cls[k].n +=1

	# check whether there are empty classes
	for k in xrange(K):
		if(Cls[k].n == 0):
			ERROR("empty class", "Wiener_CSQ", 0)
			return  [],[]

	ctf_temp = [0.0]*lctf
	for k in xrange(K):
		for i in xrange(lctf): ctf_temp[i] = 1.0/(ctf_2[k][i] + 1.0/snr)
		avek = filt_table( Cls[k].C, ctf_temp )
		for im in xrange(len(data)):
			if(k == assign[im]):
				Cls[k].SSE += data[im].cmp("SqEuclidean", filt_table( avek, ctf1[im] ))

	return Cls, ctf_2
	
def Wiener_sse(data, K, assign, Cls, ctf1, ctf2, snr = 1.0):
	from filter       import filt_table
	from fundamentals import fft
	N = len(data)
	lctf = len(ctf2[0])
	ctf_2 = [[0.0]*lctf for k in xrange(K)]

	for k in xrange(K):
		Cls[k].C.to_zero()
		Cls[k].n = 0
		Cls[k].SSE = 0.0

	#Calculate averages
	for im in xrange(N):
		ctf_x = filt_table( data[im], ctf1[im] )
		k = assign[im]
		Util.add_img(Cls[k].C, ctf_x)
		#Cls[k].C += ctf_x
		for i in xrange(lctf): ctf_2[k][i] += ctf2[im][i]
		Cls[k].n +=1

	# check whether there are empty classes
	for k in xrange(K):
		if(Cls[k].n == 0):
			return  []
			ERROR("empty class", "Wiener_sse", 1)

	#  Calculate partial SSEs and total SSE
	for k in xrange(K):
		for i in xrange(lctf):  ctf_2[k][i] = 1.0/(ctf_2[k][i] + 1.0/snr)
		Cls[k].C = filt_table( Cls[k].C, ctf_2[k] )

	for im in xrange(len(data)):
		k = assign[im]
		ctf_x = filt_table( Cls[k].C, ctf1[k] )
		Cls[k].SSE += data[im].cmp("SqEuclidean", ctf_x)

	return Cls

def k_means_aves(images, N, K, rand_seed, outdir):
	from utilities import model_blank
	from random       import seed, randint
	from sys import exit
	
	#  This version is not quite correct, as ctf2 is not modified in the avergae, but it cannot be done any other simple way :-(
	seed(rand_seed)
	assign = [0]*N
	nx = images[0].get_xsize()
	ny = images[0].get_ysize()
	nz = images[0].get_zsize()
	norm = nx*ny*nz
	e = model_blank(nx, ny, nz)

	Cls  = [0]*K #  number of elements in each class
	suma = []
	for k in xrange(K): suma.append(e.copy())
	#initialize randomly
	init_method = "Random"
	for i in xrange(3,6): assign[i]=1
	if(init_method == "Random"):
		retrial = 20
		while (retrial > 0):
			retrial -= 1
			#for im in xrange(N):  assign[im] = randint(0, K-1)
			for im in xrange(N):  Cls[assign[im]] +=1
			ier = 1
			k = K
			while( k > 0 and ier):
				k -= 1
				if(Cls[k] == 0):
					ier =0
					if(retrial == 0): ERROR(" Empty class in the initialization", "init_kmeans_aves", 1)
					for k in xrange(K): Cls[k] = 0
			if(ier == 1):  retrial = 0
	#  Calculate sums of images
	for im in xrange(N): Util.add_img(suma[assign[im]], images[im])
	# calculate criterion
	CRIT = [0.0]*K
	CT = 0.0
	for k in xrange(K):
		suma[k] = Util.mult_scalar(suma[k], 1.0/float(Cls[k]))
		CRIT[k] += norm*suma[k].cmp("dot",suma[k],{"negative":0})
		CT += CRIT[k]
		print  k,"  ",CRIT[k]
	print  " TOTAL = ",CT
	print  " Cls ",Cls
	# clustering
	from random import shuffle
	order = range(N)
	cnt = 0
	change = True
	maxit = 100
	while (change and cnt < maxit):
		#shuffle(order)
		cnt += 1
		change = False
		for imn in xrange(N):
			im = order[imn]
			#current assignment
			ca = assign[im]
			e = Util.subn_img( Util.mult_scalar(suma[ca], float(Cls[ca])), images[im] )
			JMINUS = norm*e.cmp("dot",e,{"negative":0})
			print  "JMINUS  ",im,ca, JMINUS
			JMINUS /= float(Cls[ca]-1)**2
			print  "JMINUS  ",im,ca, JMINUS
			Gain = [0.0]*K
			Gain_max = 0.0
			assign_to = -1
			for k in xrange(K):
				#  consider moving object assign[im] to remining K-1 groups
				if(k != ca):
					e = Util.subn_img( Util.mult_scalar(suma[k], float(Cls[k])), images[im] )
					JPLUS = norm*e.cmp("dot",e,{"negative":0})
					JPLUS /= float(Cls[k]+1)**2
					Gain[k] = JMINUS + JPLUS - CRIT[ca] - CRIT[k]
					print  "gain  ",k,"  ",Gain[k], JPLUS
					if(Gain[k] > Gain_max):
						Gain_max = Gain[k]
						assign_to = k
						JPLUSk = JPLUS
						print  "ASSIGNED to ",k,Gain_max
			# if improvement, move object
			if(assign_to > -1):
				print  "   >>>>>>>>>>>>>>>>>    MOVING OBJECT",im,Cls,CRIT
				CT += Gain[assign_to]
				assign[im] = assign_to
				Cls[ca] -= 1
				CRIT[ca] = JMINUS
				Cls[assign_to] += 1	
				CRIT[assign_to] = JPLUSk
				print  "   >>>>>>>>>>>>>>>>>    MOVING OBJECT",im,assign_to,CT, Cls,CRIT
				change = True
		print " ITERATION ",cnt
		print  CT
		print  CRIT
		print  assign
	exit()

def var_bydef(vol_stack, vol_list, info):
	"""  var_bydef calculates variance from definition
	"""
	if type(vol_stack)==type(''):
	    average = EMData()
	    average.read_image( vol_stack, vol_list[0] )
	else:
	    average = vol_stack[ vol_list[0] ].copy()
	average.to_zero()

	nimg = 0
	if not(info is None):
	    info.write( "Calculating average:" ) 
	    info.flush()
	for ivol in vol_list: 
		curt = EMData()
		if type(vol_stack)==type(''): curt.read_image(vol_stack, ivol )
		else:                       curt = vol_stack[ivol]
		average += curt
		if(nimg % 50==0 and  not(info is None) ) :
			info.write( "\n" )
			info.write( " %4d " % nimg )
			info.flush()
		nimg=nimg+1
		if not(info is None):
			info.write( "." )
			info.flush( )
	if not(info is None):
	    info.write( "\n" )
	average /= len(vol_list)
	vars = average.copy()
	vars.to_zero()
	nimg=0
	if not(info is None):
		info.write( "Calculating variance:\n" )
		info.flush()
	for ivol in vol_list:
		if type(vol_stack)==type(''):
			curt = EMData()
			curt.read_image(vol_stack,ivol)
		else:
			curt = vol_stack[ivol].copy()
		curt -= average
		curt *= curt
		vars += curt
		if(nimg % 50==0 and not(info is None) ) :
			info.write( "\n" )
			info.write( " %4d " % nimg )
			info.flush()
		nimg=nimg+1
		if not(info is None):
		        info.write( "." )
		        info.flush( )
	if not(info is None):
		info.write( "\n" )
	return vars/(len(vol_list)-1)


def histogram2d( datai, dataj, nbini, nbinj ) :
	fmaxi = max( datai )
	fmini = min( datai )
	fmaxj = max( dataj )
	fminj = min( dataj )

	binsize_i = (fmaxi - fmini)/(nbini-1)
	binsize_j = (fmaxj - fminj)/(nbinj-1)
	start_i = fmini-binsize_i/2.0
	start_j = fminj-binsize_j/2.0

	region = [None]*nbini
	hist = [None]*nbinj
	for i in xrange(nbini):
		region[i] = [None]*nbinj
		hist[i] = [None]*nbinj
		for j in xrange(nbinj) :
			region[i][j] = (start_i + i*binsize_i, start_j + j*binsize_j)
			hist[i][j] = 0

	assert len(datai)==len(dataj)
	for k in xrange( len(datai) ):
		id = int( (datai[k]-start_i)/binsize_i )
		jd = int( (dataj[k]-start_j)/binsize_j )
		hist[id][jd]+=1

	return region,hist

def linreg(X, Y):
	"""
	  Linear regression y=ax+b
	"""
	Sx = Sy = Sxx = Syy = Sxy = 0.0
	N = len(X)
	for x, y in map(None, X, Y):
		Sx  += x
		Sy  += y
		Sxx += x*x
		Syy += y*y
		Sxy += x*y
	det = Sxx * N - Sx * Sx
	a, b = (Sxy * N - Sy * Sx)/det, (Sxx * Sy - Sx * Sxy)/det
	"""
	from math import sqrt
	meanerror = residual = 0.0
	for x, y in map(None, X, Y):
		meanerror += (y - Sy/N)**2
		residual  += (y - a * x - b)**2
	RR = 1.0 - residual/meanerror
	ss = residual / (N-2)
	Var_a, Var_b = ss * N / det, ss * Sxx / det
	print "y=ax+b"
	print "N= %d" % N
	print "a= %g \\pm t_{%d;\\alpha/2} %g" % (a, N-2, sqrt(Var_a)) 
	print "b= %g \\pm t_{%d;\\alpha/2} %g" % (b, N-2, sqrt(Var_b))
	print "R^2= %g" % RR
	print "s^2= %g" % ss
	"""
	return a,b

def pearson(X, Y):
	"""
	  Pearson correlation coefficient
	"""
	from math import sqrt
	Sx = Sy = Sxx = Syy = Sxy = 0.0
	N = len(X)
	for x, y in map(None, X, Y):
		Sx  += x
		Sy  += y
		Sxx += x*x
		Syy += y*y
		Sxy += x*y
	return (Sxy - Sx * Sy / N) / sqrt((Sxx - Sx*Sx/N)*(Syy - Sy*Sy/N))

def table_stat(X):
	"""
	  Basic statistics of numbers stored in a list
	"""
	av = X[0]
	va = X[0]
	mi = X[0]
	ma = X[0]
	for i in xrange(1,len(X)):
		av += X[i]
		va += X[i]*X[i]
		mi = min(mi, X[i])
		ma = max(ma, X[i])
	return  av/len(X),(va - av*av/len(X))/float(len(X)-1) , mi, ma

def get_power_spec(stack_file, start_particle, end_particle):
	# computes the rotationally averaged power spectra of a series of images, e.g. a defocus group
	# and averages these spectra into one spectrum for the whole set of images
	# returns a list
	from fundamentals import rops_table
	from utilities import get_im

	ima = get_im(stack_file, start_particle)
	PSrot_avg = rops_table(ima)
	nnn = len(PSrot_avg)
	q= end_particle-start_particle+1
	for i in xrange(start_particle+1, end_particle+1):
		ima = get_im(stack_file, i)
		PSrot = rops_table(ima)
		for y in xrange(nnn): PSrot_avg[y] += PSrot[y]
	for y in xrange(nnn): PSrot_avg[y] /= q
	return PSrot_avg
	
def noise_corrected_PW(pw, lo_limit, hi_limit, abs_limit):
	from math import log, sqrt, exp
	# pw 	    : a list containing the values of a power spectrum		
	# lo_limit  : lower frequency threshold for minima search -> lt in freq. units
	# hi_limit  : upper frequency threshold for minima search -> ut in freq. units
	# abs_limit : upper limit of the power spectrum

	# returns a noise corrected power spectrum, a list of corresponding frequency values and the factors a and b
	# which were used to substract the function f(x)=exp( a*x*x+b ) from the original power spectrum pw 


	def fit_min(lt,ut,pw):
		# finds minima of given rot. averaged power spectrum of lengh nnn
		# searches between frequency domain boundaries
		xx = []
		yy = []
		nnn = len(pw)
		for k in xrange(lt,ut+1):
			if(pw[k-1] > pw[k] and pw[k+1] > pw[k]):
				xx.append(pow(float(k)*0.5/nnn,2))
				yy.append(log(pw[k]))	
		if (len(xx)==0):
			for k in xrange(lt,ul+1):
				xx.append(pow(float(k)*0.5/nnn,2))
				yy.append(log(pw[k]))
		else:
			for k in xrange(ut,ul+1):
				xx.append(pow(float(k)*0.5/nnn,2))
				yy.append(log(pw[k]))
		return 	linreg(xx, yy)

	def get_ns_pw(pw_in,a,b):
		# computes a noise substracted power spectrum 
		# based on model f(x)=a*x*x+b
		pw_sub = []
		nnn = len(pw_in)
		freq =[]
		for k in xrange(nnn):
			f = float(k)*0.5/nnn
			pw_sub.append(pw_in[k]-exp(a*f*f+b))
			freq.append(f)
		return pw_sub,freq
	
	nnn=len(pw)
	lt = int(lo_limit*nnn/0.5)
	ut = int(hi_limit*nnn/0.5)
	ul = int(abs_limit*nnn/0.5)

	aa,bb = fit_min(lt,ut,pw)
	pw_n,freq = get_ns_pw(pw,aa,bb)

	x2 = []
	y2 = []
	xc = []

	for k in xrange(lt,ut+1):
		if(pw_n[k-1] > pw_n[k] and pw_n[k+1] > pw_n[k]):
			x2.append(float(k))
	
	#if(len(x2) < 3):
		#print "only ",len(xc)," minima found, can only estimate noise-substracted power spectrum"
		#for k in xrange(lt,ut+1): 
			#xc.append(pow(float(k)*0.5/nnn,2))
			#y2.append(log(pw[k]))
	#else:
	for k in xrange(len(x2)): 
		x_value = int(x2[k])		
		xc.append(pow(float(x_value)*0.5/nnn,2))
		y2.append(log(pw[x_value]))

	a2,b2 = linreg(xc,y2)
	pw_n2,freq2 = get_ns_pw(pw,a2,b2) # pw_n2 is the noise corrected power spectrum after 2nd run
	return freq2, pw_n2, a2, b2


class def_variancer:
	def __init__(self, nx, ny, nz):
		from utilities import model_blank
		self.nimg = 0
		self.sum1 = model_blank(nx, ny, nz)
		self.imgs = []

	def insert(self, img):
		self.nimg += 1
		Util.add_img(self.sum1, img)

		self.imgs.append(img)

        def mpi_getvar(self, myid, rootid):
		from utilities import reduce_EMData_to_root, bcast_EMData_to_all
		from mpi import mpi_reduce, MPI_INT, MPI_SUM, MPI_COMM_WORLD
		avg = self.sum1.copy()

		reduce_EMData_to_root( avg, myid, rootid )
		nimg = mpi_reduce( self.nimg, 1, MPI_INT, MPI_SUM, rootid, MPI_COMM_WORLD)
		
		if myid==rootid:
   		    nimg = int(nimg[0])
		    avg /= nimg

		bcast_EMData_to_all( avg, myid, rootid )

		var = avg.copy()
		var.to_zero()
		for img in self.imgs:
			s = img - avg
			Util.add_img2( var, s )
	
		reduce_EMData_to_root( var, myid, rootid )
		if myid==rootid:
			var /= (nimg-1)
			var.set_attr( "nimg", nimg )
			return var, avg

		return None, None


        def mpi_getavg(self, myid, rootid ):
		from mpi import mpi_reduce, MPI_INT, MPI_SUM, MPI_COMM_WORLD
		from utilities import reduce_EMData_to_root

		cpy1 = self.sum1.copy()

		reduce_EMData_to_root( cpy1, myid, rootid )
		
		nimg = mpi_reduce( self.nimg, 1, MPI_INT, MPI_SUM, rootid, MPI_COMM_WORLD)
		
		if myid==rootid:
			nimg = int( nimg[0] )
			return cpy1/nimg

		return None

	def getvar(self):
		avg1 = self.sum1/self.nimg

		tmp = avg1.copy()
		Util.mul_img( tmp, avg1 )
		avg2 -= tmp

		avg2 *= (float(self.nimg)/float(self.nimg-1))
		 
		return avg2


	def getavg(self):
		return self.sum1/self.nimg


class inc_variancer:
	def __init__(self):
		self.nimg = 0
		self.sum1 = None
		self.sum2 = None

	def insert(self, img):
		self.nimg += 1

		if self.sum1 is None:
			assert self.sum2 is None
			self.sum1 = img.copy()
			self.sum2 = img.copy()
			Util.mul_img(self.sum2, img)
		else:
			Util.add_img(self.sum1, img)
			Util.add_img2(self.sum2, img)


        def mpi_getvar(self, myid, rootid):
		from utilities import reduce_EMData_to_root
		from mpi import mpi_reduce, MPI_INT, MPI_SUM, MPI_COMM_WORLD
		cpy1 = self.sum1.copy()
		cpy2 = self.sum2.copy()

		reduce_EMData_to_root( cpy1, myid, rootid )
		reduce_EMData_to_root( cpy2, myid, rootid )
		sum_nimg = mpi_reduce( self.nimg, 1, MPI_INT, MPI_SUM, rootid, MPI_COMM_WORLD)
		if myid==rootid:

   		    sum_nimg = int(sum_nimg[0])


		    avg1 = cpy1/sum_nimg
		    avg2 = cpy2/sum_nimg


		    tmp = avg1.copy()
		    Util.mul_img( tmp, avg1 )
		    avg2 -= tmp
		    avg2 *= (float(sum_nimg)/float(sum_nimg-1))
		    return avg2
		return None


        def mpi_getavg(self, myid, rootid ):
		from mpi import mpi_reduce, MPI_INT, MPI_SUM, MPI_COMM_WORLD
		from utilities import reduce_EMData_to_root

		cpy1 = self.sum1.copy()

		reduce_EMData_to_root( cpy1, myid, rootid )
		
		sum_nimg = mpi_reduce( self.nimg, 1, MPI_INT, MPI_SUM, rootid, MPI_COMM_WORLD)
		
		if myid==rootid:
			sum_nimg = int( sum_nimg[0] )
			return cpy1/sum_nimg

		return None

	def getvar(self):
		avg1 = self.sum1/self.nimg
		avg2 = self.sum2/self.nimg

		tmp = avg1.copy()
		Util.mul_img( tmp, avg1 )
		avg2 -= tmp

		avg2 *= (float(self.nimg)/float(self.nimg-1))
		 
		return avg2


	def getavg(self):
		return self.sum1/self.nimg

def mono(k1,k2):
	"""
	get index of a square nxn matrix stored in a triangular form
	for i in xrange(1,n):
	    for j in xrange(i):
		print  i,j,mono(i,j)

	"""
	mk = max(k1,k2)
	return  min(k1,k2) + mk*(mk-1)/2
	
def cluster_pairwise(d, K):
	"""
	  d  - lower half of the square matrix of pairwsie distances
	  K  - number of classes
	"""
	from statistics import mono
	from random import randint, shuffle
	from math import sqrt
	N = 1 + int((sqrt(1.0 + 8.0*len(d))-1.0)/2.0)
	if(N*(N-1)/2 != len(d)):
		print  "  incorrect dimension"
		return
	cent = [0]*K
	assign = range(N)
	shuffle(assign)
	for k in xrange(K):  cent[k] = assign[k]
	# assign
	assign = [0]*N
	change = True
	it = -1
	while(change):
		change = False
		it += 1
		#print  "Iteration  ", it
		# dispersion is a sum of distances from objects to group centers
		disp = 0.0
		for i in xrange(N):
			qm = 1.0e23
			for k in xrange(K):
				if(i == cent[k]):
					qm = 0.0
					na = k
				else:
					dt = d[mono(i,cent[k])]
					if(dt < qm):
						qm = dt
						na = k
			disp += qm
			if(na != assign[i]):
				assign[i] = na
				change = True
		#print disp
		#print  assign
		# find centers
		for k in xrange(K):
			qm = 1.0e23
			for i in xrange(N):
				if(assign[i] == k):
					q = 0.0
					for j in xrange(N):
						if(assign[j] == k):
							if(i != j):
								#it cannot be the same object
								q += d[mono(i,j)]
					if(q < qm):
						qm = q
						cent[k] = i
	return  assign,cent,disp,it
	"""
	# write out information
	from utilities import write_text_file
	write_text_file(cent, "cent")
	for k in xrange(K):
		cent = []
		for i in xrange(N):
			if(assign[i] == k):  cent.append(i)
		write_text_file(assign, "assign%03d"%(k))
	"""	
	
	
def cluster_equalsize(d, m):
	"""
	  d  - lower half of the square matrix of pairwsie distances
	  m  - number of objects per class
	"""
	from statistics import mono
	from random import randint, shuffle
	from math import sqrt
	nd = d.get_xsize()
	N = 1 + int((sqrt(1.0 + 8.0*nd)-1.0)/2.0)
	if(N*(N-1)/2 != nd):
		print  "  incorrect dimension"
		return
	K = N/m
	active = [True]*N
	groupping = [None]*K
	for k in xrange(K):
		# find two most similiar objects among active
		dm = 1.0e23
		print 'K:', k
		f = open('WATCH', 'a')
		f.write('K: %d\n' % k)
		f.close()
		for i in xrange(1,N):
			if(active[i]):
				for j in xrange(i):
					if(active[j]):
						qd = d.get_value_at(mono(i,j))
						if(qd < dm):
							dm = qd
							p = [i,j]
		groupping[k] = p
		active[p[0]] = False
		active[p[1]] = False
		#print  dm
		#print  groupping[k]
		#print active
		# find progressively objects most similar to those in the current list
		for l in xrange(2,m):
			dm = 1.0e23
			for i in xrange(N):
				if(active[i]):
					qt = 0.0
					for j in groupping[k]:
						qt += d.get_value_at(mono(i,j))
					if(qt < dm):
						dm = qt
						p = i
			groupping[k].append(p)
			active[p] = False
			#print  dm
			#print  groupping[k]
			#print active
	# there might be remaining objects when N is not divisible by m, simply put them in one group
	if(N%m != 0):
		K += 1
		groupping.append([])
		for i in xrange(N):
			if(active[i]):
				groupping[K-1].append(i)
	# find centers
	cent = [0]*K
	for k in xrange(K):
		dm = 1.0e23
		for i in groupping[k]:
				q = 0.0
				for j in groupping[k]:
					if(i != j):
						#it cannot be the same object
						q += d.get_value_at(mono(i,j))
				if(q < dm):
					dm = q
					cent[k] = i
	# dispersion is a sum of distances from objects to group centers
	disp = 0.0
	for k in xrange(K):
		for i in groupping[k]:
			if(i != cent[k]):
				disp += d.get_value_at(mono(i,cent[k]))
	print  disp
	# try swapping elements
	for k1 in xrange(1,K):
		for k2 in xrange(k1):
			print  " trying to swap",k1,k2
			for i in groupping[k1]:
				if(i != cent[k1]): d1 = d.get_value_at(mono(i,cent[k1]))
				else:              d1 = 0.0
				d2 = d.get_value_at(mono(i,cent[k2]))
				for j in groupping[k2]:
					if(j != cent[k2]): e1 = d.get_value_at(mono(j,cent[k2]))
					else:              e1 = 0.0
					e2 = d.get_value_at(mono(j,cent[k1]))
					if(d1+e1 > d2+e2):
						#  swap!
						print  " SWAP"
						l1 = groupping[k1].index(i)
						l2 = groupping[k2].index(j)
						temp = groupping[k1][l1]
						groupping[k1][l1] = groupping[k2][l2]
						groupping[k1][l2] = temp
	# dispersion is a sum of distance from objects to group centers
	disp = 0.0
	for k in xrange(K):
		for i in groupping[k]:
			if(i != cent[k]):
				disp += d.get_value_at(mono(i,cent[k]))
	print  disp
	for k in xrange(K):
		groupping[k].sort()
	
	return  groupping,cent,disp
