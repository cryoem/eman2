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

def avgvar(data, mode='a', interp='quadratic', i1=0, i2=0, use_odd=True, use_even=True):
	'''
	
	INPUT
	
	data: image stack, can be 2D or 3D, must be in real space
	mode: whether to apply alignment parameters. Default mode='a' means apply parameters
	rot_method: specifies the function by which images are rotated/shifted if alignment parameters are to be applied. This is only relevant for the case where images are 2D, in which case rot_method can be either rot_shift2D or rotshift2dg, with the default being rot_shift2D. If images are 3D, rot_shift3D will be used to rotate/shift the images.
	interp: interpolation method to use for rot_method when applying alignment parameters.
	i1: index of first image to be used.
	i2: index of last image to be used. If i2 = 0, then i2 defaults to one less than number of images in the data
	use_odd: images with indices between i1 and i2 which are odd are used if and only if use_odd is set to True. Default is True.
	use_even: images with indices between i1 and i2 which are even are used if and only if use_even is set to True. Default is True.
	
	OUTPUT
		
	ave: the average of the image series in real space
	var: the variance of the image series in real space

	'''
	from utilities    import model_blank
	from alignment    import kbt

	inmem = True
	if type(data) == type(""):
		inmem = False
		from utilities    import get_im	

	img2D = True
	if inmem:
		img = data[0]
	else:
		img = get_im(data,0)
	nx = img.get_xsize()
	ny = img.get_ysize()
	nz = img.get_zsize()
	if nz > 1:
		img2D = False

	if mode == 'a':
		if img2D:
			from utilities import get_params2D
			from fundamentals import rot_shift2D
		else:
			from utilities import get_params3D
			from fundamentals import rot_shift3D

	if inmem:
		data_nima = len(data)
	else:
		data_nima = EMUtil.get_image_count(data)
	if i2 == 0: i2 = data_nima-1

	ave = model_blank(nx,ny,nz)
	var = model_blank(nx,ny,nz)
	nima = 0
	for i in xrange(i1, i2+1):
		if not(use_odd) and i%2 == 1:
			continue
		if not(use_even) and i%2 == 0:
			continue
		nima += 1
		if inmem:
			img = data[i]
		else:
			img = get_im(data, i)
		if (mode == 'a'):
			if img2D:
				angle, sx, sy, mirror, scale = get_params2D(img)
				img = rot_shift2D(img, angle, sx, sy, mirror, scale, interp)
			else:
				phi, theta, psi, s3x, s3y, s3z, mirror, scale = get_params3D(img)
				img = rot_shift3D(img, phi, theta, psi, s3x, s3y, s3z, scale)
		Util.add_img(ave, img)
		Util.add_img2(var, img)

	Util.mul_scalar(ave, 1.0 /float(nima) )
	return ave, (var - ave*ave*nima)/(nima-1)

def avgvar_ctf(data, mode='a', interp='quadratic', i1=0, i2=0, use_odd=True, use_even=True, snr=1.0e23, dopa = True):
	'''
	
	INPUT
	
	data: image stack, must be 2D, must be in real space
	mode: whether to apply alignment parameters. Default mode='a' means apply parameters
	rot_method: specifies the function by which images are rotated/shifted if alignment parameters are to be applied. This is only relevant for the case where images are 2D, in which case rot_method can be either rot_shift2D or rotshift2dg, with the default being rot_shift2D. If images are 3D, rot_shift3D will be used to rotate/shift the images.
	interp: interpolation method to use for rot_method when applying alignment parameters.
	i1: index of first image to be used.
	i2: index of last image to be used. If i2 = 0, then i2 defaults to one less than number of images in the data
	use_odd: images with indices between i1 and i2 which are odd are used if and only if use_odd is set to True. Default is True.
	use_even: images with indices between i1 and i2 which are even are used if and only if use_even is set to True. Default is True.
	snr: signal to noise ratio, default 1.0
	
	OUTPUT
	
	tavg: The best estimate (Wiener filter) given the image series and estimated CTF parms, in real space.
			
		var: Variance (in real space) calculated as follows, [1/(n-1)]*[sum_j { F[(H_j*(O_j - F^{-1}(H_j*tavg))^2]/SUM_CTF^2} where O_j is the j-th image in real space, F^{-1} denotes inverse fourier transform operator, and H_j is the CTF of the j-th image
	
	'''
	
	from utilities    import model_blank, pad
	from alignment    import kbt
	from fundamentals import fft, fftip, window2d
	from filter       import filt_ctf
	from morphology   import ctf_img

	inmem = True
	if type(data) == type(""):
		inmem = False
		from utilities    import get_im	

	if inmem:
		img = data[0]
	else:
		img = get_im(data,0)
	nx = img.get_xsize()
	ny = img.get_ysize()
	nz = img.get_zsize()
	if nz > 1:
		print "images must be 2D for CTF correction.....exiting"
		sys.exit()

	if img.get_attr_default('ctf_applied', 0) == 1:
		print "data cannot be ctf-applied....exiting"
		sys.exit()

	if mode == 'a':
		from utilities import get_params2D
		from fundamentals import rot_shift2D

	if inmem:
		data_nima = len(data)
	else:
		data_nima = EMUtil.get_image_count(data)

	if i2 == 0: i2 = data_nima-1
	if dopa:
		nx2 = nx*2
		ny2 = ny*2
	else:
		nx2 = nx
		ny2 = ny
	ave = EMData(nx2, ny2, 1, False)
	ctf_2_sum = EMData(nx2, ny2, 1, False)
	nima = 0
	for i in xrange(i1, i2+1):
		if not(use_odd) and i%2 == 1:
			continue
		if not(use_even) and i%2 == 0:
			continue
		nima += 1
		if inmem:
			img = data[i].copy()
		else:
			img = get_im(data, i)

		if (mode == 'a'):
			angle, sx, sy, mirror, scale = get_params2D(img)
			img = rot_shift2D(img, angle, sx, sy, mirror, scale, interp)

		img = pad(img, nx2, ny2, 1, background = "circumference")
		ctf_params = img.get_attr("ctf")
		fftip(img)
		Util.add_img(ave, filt_ctf(img, ctf_params))
		Util.add_img2(ctf_2_sum, ctf_img(nx2, ctf_params))

	snr_img = EMData(nx2, ny2, 1, False)
	fnx = snr_img.get_xsize()
	fny = snr_img.get_ysize()
	for i in xrange(fnx):
		for j in xrange(fny):
			snr_img.set_value_at(i,j, 1.0/snr)
	Util.add_img(ctf_2_sum, snr_img)
	Util.div_filter(ave, ctf_2_sum)

	# calculate variance in real space
	#totv = model_blank(nx2, ny2, nz)
	tvar = model_blank(nx, ny, nz)
	for i in xrange(i1, i2+1):
		if not(use_odd) and i%2 == 1:
			continue
		if not(use_even) and i%2 == 0:
			continue
		if inmem:
			img = data[i].copy()
		else:
			img = get_im(data, i)

		if (mode == 'a'):
			angle, sx, sy, mirror, scale = get_params2D(img)
			img = rot_shift2D(img, angle, sx, sy, mirror, scale, interp)

		img = pad(img, nx2, ny2, 1, background = "circumference")
		fftip(img)
		ctf_params = img.get_attr("ctf")
		img = filt_ctf(img-filt_ctf(ave, ctf_params, dopa), ctf_params, dopa)
		Util.div_filter(img, ctf_2_sum)
		img = window2d(fft(img),nx,ny)
		#Util.add_img(totv, img)
		Util.add_img2(tvar, img)
	Util.mul_scalar(tvar, float(nima)*nima/(nima-1)) # the strange factor is due to the fact that division by ctf^2 is equivalent to division by nima
	return  window2d(fft(ave),nx,ny) , tvar#,(tvar - totv*totv/nima), tvar, totv,tavg

def add_oe_series(data, ali_params="xform.align2d"):
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
		alpha, sx, sy, mirror, scale = get_params2D(data[i], ali_params)
		temp = rot_shift2D(data[i], alpha, sx, sy, mirror, scale, "quadratic")
		if i%2 == 0: Util.add_img(ave1, temp)
		else:         Util.add_img(ave2, temp)
	return ave1, ave2

def add_ave_varf(data, mask = None, mode = "a", CTF = False, ctf_2_sum = None, ali_params = "xform.align2d"):
	"""
		Calculate average of an image series and variance, sum of squares in Fourier space
		mode - "a": use current alignment parameters
		CTF  - if True, use CTF for calculations of both average and variance.
	"""
	from utilities    import    model_blank, get_params2D, info
	from fundamentals import    rot_shift2D, fft, fftip

	n = len(data)
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	ave1 = EMData(nx, ny, 1, False)
	ave2 = EMData(nx, ny, 1, False)
	var  = EMData(nx, ny, 1, False)
	
	if CTF:
		from morphology   import ctf_img
		from filter       import filt_ctf, filt_table
		if data[0].get_attr_default('ctf_applied', 1) == 1:
			ERROR("data cannot be ctf-applied", "add_ave_varf", 1)
		if ctf_2_sum:  get_ctf2 = False
		else:          get_ctf2 = True
		if get_ctf2: ctf_2_sum = EMData(nx, ny, 1, False)
	 	for i in xrange(n):
	 		if mode == "a":
				alpha, sx, sy, mirror, scale = get_params2D(data[i], ali_params)
				ima = rot_shift2D(data[i], alpha, sx, sy, mirror, scale, "quadratic")
				if mask:  Util.mul_img(ima, mask)
				fftip(ima)
				#  Here we have a possible problem: varf works only if CTF is applied after rot/shift
				#    while calculation of average (and in general principle) CTF should be applied before rot/shift
				#    here we use the first possibility
			else:
				if  mask:   ima = fft(Util.muln_img(data[i], mask))
				else:       ima = fft(data[i])
	 		ctf_params = data[i].get_attr("ctf")
	 		ima_filt = filt_ctf(ima, ctf_params, dopad=False)
			if(i%2 == 0):  Util.add_img(ave1, ima_filt)
			else:          Util.add_img(ave2, ima_filt)
 			Util.add_img2(var, ima)
	 		if get_ctf2: Util.add_img2(ctf_2_sum, ctf_img(nx, ctf_params))
		sumsq = Util.addn_img(ave1, ave2)
		tavg = Util.divn_img(sumsq, ctf_2_sum)
		Util.mul_img(sumsq, sumsq.conjg())
		Util.div_img(sumsq, ctf_2_sum)
	 	Util.sub_img(var, sumsq)
	else:
		for i in xrange(n):
			if mode == "a":
				alpha, sx, sy, mirror, scale = get_params2D(data[i], ali_params)
				ima = rot_shift2D(data[i], alpha, sx, sy, mirror, scale, "quadratic")
				if mask:  Util.mul_img(ima, mask)
				fftip(ima)
			else:
				if  mask:   ima = fft(Util.muln_img(data[i], mask))
				else:       ima = fft(data[i])
			if(i%2 == 0):   Util.add_img(ave1, ima)
			else:           Util.add_img(ave2, ima)
			Util.add_img2(var, ima)
		sumsq = Util.addn_img(ave1, ave2)
		tavg = Util.mult_scalar(sumsq, 1.0/float(n))
		Util.mul_img(sumsq, sumsq.conjg())
		Util.mul_scalar(sumsq, 1.0/float(n))
		Util.sub_img(var, sumsq)

	Util.mul_scalar(var, 1.0/float(n-1))
	var.set_value_at(0, 0, 1.0)
	st = Util.infomask(var, None, True)
	if st[2] < 0.0:  ERROR("Negative variance!", "add_ave_varf", 1)
	return tavg, ave1, ave2, var, sumsq

def add_ave_varf_MPI(myid, data, mask = None, mode = "a", CTF = False, ctf_2_sum = None, ali_params = "xform.align2d", main_node = 0, comm = -1):
	"""
		Calculate sum of an image series and sum of squares in Fourier space
		Since this is the MPI version, we need to reduce sum and sum of squares 
		on the main node and calculate variance there.
		mode - "a": use current alignment parameters
		CTF  - if True, use CTF for calculations of the sum.
	"""
	from utilities    import    model_blank, get_params2D
	from fundamentals import    rot_shift2D, fft, fftip
	from utilities    import    reduce_EMData_to_root
	from mpi          import    mpi_reduce, MPI_INT, MPI_SUM
	
	if comm == -1:
		from mpi import MPI_COMM_WORLD
		comm = MPI_COMM_WORLD

	n = len(data)
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	ave1 = EMData(nx, ny, 1, False)
	ave2 = EMData(nx, ny, 1, False)
	var  = EMData(nx, ny, 1, False)
	
	if CTF:
		from filter       import filt_ctf
		from morphology   import ctf_img
		if data[0].get_attr_default('ctf_applied', 1) == 1:
	 		ERROR("data cannot be ctf-applied", "add_ave_varf_MPI", 1)
		if ctf_2_sum:  get_ctf2 = False
		else:          get_ctf2 = True
		if get_ctf2: ctf_2_sum = EMData(nx, ny, 1, False)
	 	for i in xrange(n):
	 		if mode == "a":
				alpha, sx, sy, mirror, scale = get_params2D(data[i], ali_params)
				ima = rot_shift2D(data[i], alpha, sx, sy, mirror, scale, "quadratic")
				if mask:  Util.mul_img(ima, mask)
				fftip(ima)
			else:
				if  mask:   ima = fft(Util.muln_img(data[i], mask))
				else:       ima = fft(data[i])
	 		ctf_params = data[i].get_attr("ctf")
	 		ima_filt = filt_ctf(ima, ctf_params, dopad=False)
			if(i%2 == 0):   Util.add_img(ave1, ima_filt)
			else:           Util.add_img(ave2, ima_filt)
 			Util.add_img2(var, ima)
	 		if get_ctf2: Util.add_img2(ctf_2_sum, ctf_img(nx, ctf_params))
	else:
		get_ctf2 = False
		for i in xrange(n):
			if mode == "a":
				alpha, sx, sy, mirror, scale = get_params2D(data[i], ali_params)
				ima = rot_shift2D(data[i], alpha, sx, sy, mirror, scale, "quadratic")
				if mask:  Util.mul_img(ima, mask)
				fftip(ima)
			else:
				if  mask:   ima = fft(Util.muln_img(data[i], mask))
				else:       ima = fft(data[i])
			if(i%2 == 0):   Util.add_img(ave1, ima)
			else:           Util.add_img(ave2, ima)
			Util.add_img2(var, ima)
	reduce_EMData_to_root(ave1, myid, main_node, comm)
	reduce_EMData_to_root(ave2, myid, main_node, comm)
	reduce_EMData_to_root(var, myid, main_node, comm)
	if get_ctf2: reduce_EMData_to_root(ctf_2_sum, myid, main_node, comm)
	nima = n
	nima = mpi_reduce(nima, 1, MPI_INT, MPI_SUM, main_node, comm)
	if myid == main_node:
		nima = int(nima)
		sumsq = Util.addn_img(ave1, ave2)
		if CTF:
			tavg = Util.divn_img(sumsq, ctf_2_sum)
			Util.mul_img(sumsq, sumsq.conjg())
			Util.div_img(sumsq, ctf_2_sum)
		else:
			tavg = Util.mult_scalar(sumsq, 1.0/float(nima))
			Util.mul_img(sumsq, sumsq.conjg())
			Util.mul_scalar(sumsq, 1.0/float(nima))
		Util.sub_img(var, sumsq)
		Util.mul_scalar(var, 1.0/float(nima-1))
		var.set_value_at(0, 0, 1.0)
		st = Util.infomask(var, None, True)
		if st[2] < 0.0:  ERROR("Negative variance!", "add_ave_varf_MPI", 1)
	else:
		tavg  = EMData()
		sumsq = EMData()
	return tavg, ave1, ave2, var, sumsq

def sum_oe(data, mode = "a", CTF = False, ctf_2_sum = None):
	"""
		Calculate average of an image series
		mode - "a": use current alignment parameters
		CTF  - if True, use CTF for calculations of the average.
		In addition, calculate odd and even sums, these are not divided by the ctf^2
		If ctf_2_sum not provided, sum of ctf^2 will be calculated and returned
	"""
	from utilities    import    model_blank, get_params2D
	from fundamentals import    rot_shift2D, fft

	n      = len(data)
	nx     = data[0].get_xsize()
	ny     = data[0].get_ysize()
	ave1   = model_blank(nx, ny)
	ave2   = model_blank(nx, ny)

	if CTF:
		from morphology   import ctf_img
		from filter       import filt_ctf
		if data[0].get_attr_default('ctf_applied', 1) == 1:  ERROR("data cannot be ctf-applied", "sum_oe", 1)
		if ctf_2_sum:  get_ctf2 = False
		else:          get_ctf2 = True
		if get_ctf2: ctf_2_sum = EMData(nx, ny, 1, False)
	 	for i in xrange(n):
	 		ctf_params = data[i].get_attr("ctf")
	 		if mode == "a":
				alpha, sx, sy, mirror, scale = get_params2D(data[i])
				ima = rot_shift2D(data[i], alpha, sx, sy, mirror, scale, "quadratic")
			else:
				ima = data[i]
	 		ima_filt = filt_ctf(ima, ctf_params, dopad=True)
			if i%2 == 0:	Util.add_img(ave1, ima_filt)
			else:	        Util.add_img(ave2, ima_filt)
	 		if get_ctf2: Util.add_img2(ctf_2_sum, ctf_img(nx, ctf_params))
	else:
		for i in xrange(n):
			if mode == "a":
				alpha, sx, sy, mirror, scale = get_params2D(data[i])
				ima = rot_shift2D(data[i], alpha, sx, sy, mirror, scale, "quadratic")
			else:
				ima = data[i]
			if i%2 == 0:	Util.add_img(ave1, ima)
			else:	        Util.add_img(ave2, ima)

	if  CTF:
		if get_ctf2: return ave1, ave2, ctf_2_sum
		else:        return  ave1, ave2
	else:        return  ave1, ave2

def ave_var(data, mode = "a", listID=None):
	"""
		Calculate average and variance of a 2D or 3D image series
		with optional application of orientation parameters
		data can be either in-core stack or a disk file
	"""
	from utilities import model_blank, get_im
	if  type(data) == type(""): n = EMUtil.get_image_count(data)
	else:                       n = len(data)
	if listID == None:
		listID = range(n)
	img = get_im(data, 0)
	nx = img.get_xsize()
	ny = img.get_ysize()
	nz = img.get_zsize()
	if(mode == "a"):
		if(nz > 1):
			ali_params = "xform.align3d"
			from fundamentals import rot_shift3D
			from utilities import get_params3D
		else:
			ali_params = "xform.align2d"
			from fundamentals import rot_shift2D
			from utilities import get_params2D

	ave = model_blank(nx,ny,nz)
	var = model_blank(nx,ny,nz)
	nlistID = len(listID)
	for i in xrange(nlistID):
		img = get_im(data,listID[i])
		if(mode == "a"):
			if(nz > 1):
				phi, theta, psi, s3x, s3y, s3z, mirror, scale = get_params3D(img)
				img = rot_shift3D(img, phi, theta, psi, s3x, s3y, s3z, scale)
			else:
				angle, sx, sy, mirror, scale = get_params2D(img)
				img = rot_shift2D(img, angle, sx, sy, mirror, scale)
		Util.add_img(ave, img)
		Util.add_img2(var, img)
	Util.mul_scalar(ave, 1.0 /float(nlistID) )

	return ave, (var - ave*ave*nlistID)/(nlistID-1)

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

def ave_series(data, pave = True):
	"""
		Calculate average of a image series using current alignment parameters
		data - real space image series
	"""
	from utilities    import model_blank, get_params2D
	from fundamentals import rot_shift2D
	n = len(data)
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	ave = model_blank(nx, ny)
	for i in xrange(n):
	 	alpha, sx, sy, mirror, scale = get_params2D(data[i])
		temp = rot_shift2D(data[i], alpha, sx, sy, mirror)
		Util.add_img(ave, temp)
	if pave:  Util.mul_scalar(ave, 1.0/float(n))
	return ave

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

'''
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
		if  mirror: temp.process_inplace("xform.mirror",{"axis":'x'})
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
		if  mirror: temp.process_inplace("xform.mirror",{"axis":'x'})
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
		if  mirror: temp.process_inplace("xform.mirror", {"axis":'x'})
		if i%2 == 0: Util.add_img(ave1, temp)
		else:        Util.add_img(ave2, temp)
	ave1 /= (n/2+(n%2))
	ave2 /= (n/2)
	return ave1, ave2
'''

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


def ave_oe_series_textfile(stack, textfile):
	"""
		Calculate odd and even averages of an image stack using alignment parameters in a text file
	"""
	from utilities import model_blank, read_text_file
	from fundamentals import rot_shift2D
	
	n = EMUtil.get_image_count(stack)
	ima = EMData()
	ima.read_image(stack, 0, True)
	nx = ima.get_xsize()
	ny = ima.get_ysize()
	ave1 = model_blank(nx, ny)
	ave2 = model_blank(nx, ny)
	params = read_text_file(textfile, -1)
	for i in xrange(n):
		ima = EMData()
		ima.read_image(stack, i)
		alpha = params[0][i]
		sx = params[1][i]
		sy = params[2][i]
		mirror = params[3][i]
		temp = rot_shift2D(ima, alpha, sx, sy, mirror)
		if i%2 == 0: Util.add_img(ave1, temp)
		else:        Util.add_img(ave2, temp)
	return ave1/(n/2+n%2), ave2/(n/2)


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
	from utilities import model_blank
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
			if  mirror: temp.process_inplace("xform.mirror", {"axis":'x'})
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

	if(i2==0):
		if  type(stack) == type(""): i2 = EMUtil.get_image_count(stack)-1
		else:                       i2 = len(stack)-1
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
	e = model_blank(nx, ny, nz)
	Util.add_img2(e, ave)
	var = Util.madn_scalar(var, e, -float(ii))
	Util.mul_scalar(var, 1.0/float(ii-1))
	
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
	from utilities    import get_im, model_blank, get_params2D
	from fundamentals import rot_shift2D

	if i2 == 0:
		if type(stack) == type(""):  i2 = EMUtil.get_image_count(stack)-1
		else:  i2 = len(stack)-1
	nima = i2-i1+1

	ima = get_im(stack, i1)
	nx  = ima.get_xsize()
	ny  = ima.get_ysize()
	ave = model_blank(nx,ny)
	var = model_blank(nx,ny)
	for i in xrange(i1, i2 + 1):
		if i > i1:
			ima = get_im(stack, i)
		if mode=="a":
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
	
def aveq(stack, mode="a", i1 = 0, i2 = 0):
	"""
		Calculate the average and variance for
		1. mode="a" for alignment
		2. mode=else for normal summation
	"""
	from utilities    import get_im, model_blank, get_params2D
	from fundamentals import rot_shift2D

	if i2 == 0:
		if type(stack) == type(""):  i2 = EMUtil.get_image_count(stack)-1
		else:  i2 = len(stack)-1
	nima = i2-i1+1

	ima = get_im(stack, i1)
	nx  = ima.get_xsize()
	ny  = ima.get_ysize()
	ave = model_blank(nx,ny)
	for i in xrange(i1, i2 + 1):
		if i > i1:
			ima = get_im(stack, i)
		if mode=="a":
			alpha, sx, sy, mirror, scale = get_params2D(ima)
			Util.add_img(ave, rot_shift2D(ima, alpha, sx, sy, mirror))
		else: 
			Util.add_img(ave, ima)

	ave = Util.mult_scalar(ave, 1.0/float(nima))
	return ave

def aves_w(stack, mode="a"):
	"""
		Apply alignment parameters, and calculate Wiener average
		using CTF and SNR saved in header
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
		# # horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = e.get_attr('active')
		# if(active):
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
		Apply alignment parameters, and calculate Wiener average using CTF info
		mode="a" will apply alignment parameters to the input image.
	"""
	
	from  fundamentals import fft, rot_shift2D
	from  morphology   import ctf_img
	from  filter 	   import filt_ctf
	from  utilities    import pad, get_params2D, get_im
	from  math 	   import sqrt
	
	if type(input_stack) == type(""):	n = EMUtil.get_image_count(input_stack)
	else:  n = len(input_stack)
	ima = get_im(input_stack, 0)
	nx = ima.get_xsize()
	ny = ima.get_xsize()

	if ima.get_attr_default('ctf_applied', 2) > 0:	ERROR("data cannot be ctf-applied", "aves_wiener", 1)

	nx2 = nx*2
	ny2 = ny*2
	ave       = EMData(nx2, ny2, 1, False)
	ctf_2_sum = EMData(nx2, ny2, 1, False)
	snrsqrt = sqrt(SNR)

	for i in xrange(n):
		ima = get_im(input_stack, i)
		ctf_params = ima.get_attr("ctf")
		if mode == "a":
	 		alpha, sx, sy, mirror, scale = get_params2D(ima)
		 	ima = rot_shift2D(ima, alpha, sx, sy, mirror)
		oc = filt_ctf(fft(pad(ima, nx2, ny2, background = 0.0)), ctf_params, dopad=False)
		Util.mul_scalar(oc, SNR)
		Util.add_img(ave, oc)
		Util.add_img2(ctf_2_sum, snrsqrt*ctf_img(nx2, ctf_params, ny = ny2, nz = 1))
	ctf_2_sum += 1.0
	Util.div_filter(ave, ctf_2_sum)
	# variance
	var = EMData(nx,ny)
	var.to_zero()
	for i in xrange(n):
		ima = get_im(input_stack, i)
		ctf_params = ima.get_attr("ctf")
		if mode == "a":
			alpha, sx, sy, mirror, scale = get_params2D(ima)
 			ima = rot_shift2D(ima, alpha, sx, sy, mirror)
		oc = filt_ctf(ave, ctf_params, dopad=False)
		Util.sub_img(ima, Util.window(fft(oc),nx,ny,1,0,0,0))
		Util.add_img2(var, ima)
	ave = Util.window(fft(ave),nx,ny,1,0,0,0)
	Util.mul_scalar(var, 1.0/(n-1))
	#  The variance is incorrect, so I replaced it by a blank image, will fix later
	var.to_zero()
	return ave, var


def aves_adw(input_stack, mode="a", SNR=1.0, Ng = -1):
	"""
		Apply alignment parameters, and calculate Wiener average using CTF info
		mode="a" will apply alignment parameters to the input image.
	"""
	
	from  fundamentals import fft, rot_shift2D
	from  morphology   import ctf_img, ctf_1d, ctf_2
	from  filter 	   import filt_ctf, filt_table
	from  utilities    import pad, get_params2D, get_im
	from  math 	   import sqrt
	
	if type(input_stack) == type(""):	n = EMUtil.get_image_count(input_stack)
	else:  n = len(input_stack)
	ima = get_im(input_stack, 0)
	nx  = ima.get_xsize()

	if ima.get_attr_default('ctf_applied', 2) > 0:	ERROR("data cannot be ctf-applied", "aves_wiener", 1)

	ctf_abs_sum = EMData(nx, nx, 1, False)
	ctf_2_sum = EMData(nx, nx, 1, False)

	Ave = EMData(nx, nx, 1, False)

	if Ng == -1: Ng = n
	for i in xrange(n):
		ima = get_im(input_stack, i)
		ctf_params = ima.get_attr("ctf")
		ctfimg = ctf_img(nx, ctf_params)
		Util.add_img2(ctf_2_sum, ctfimg)
		Util.add_img_abs(ctf_abs_sum, ctfimg)
		if mode == "a":
	 		alpha, sx, sy, mirror, scale = get_params2D(ima)
		 	ima = rot_shift2D(ima, alpha, sx, sy, mirror)
		oc = filt_ctf(fft(ima), ctf_params, dopad=False)
		Util.add_img(Ave, oc)

	adw_img = Util.mult_scalar(ctf_2_sum, SNR)
	#adw_img += 1.0
	Util.div_filter(adw_img, ctf_abs_sum)
	#Util.mul_scalar(adw_img, float(Ng-1)/(n-1)/SNR)
	Util.mul_scalar(adw_img, float(Ng-1)/(n-1))
	adw_img += float(n-Ng)/(n-1)
	#Util.mul_scalar(adw_img, SNR)
	#Util.mul_scalar(ctf_2_sum, SNR)
	#ctf_2_sum += 1.0

	ave = fft(Util.divn_filter(Util.muln_img(Ave, adw_img), ctf_2_sum))

	# variance
	var = EMData(nx, nx)
	var.to_zero()
	for i in xrange(n):
		ima = get_im(input_stack, i)
		ctf_params = ima.get_attr("ctf")
		if mode == "a":
			alpha, sx, sy, mirror, scale = get_params2D(ima)
 			ima = rot_shift2D(ima, alpha, sx, sy, mirror)
		oc = filt_ctf(ave, ctf_params, dopad=False)
		Util.sub_img(ima, oc)
		Util.add_img2(var, ima)
	Util.mul_scalar(var, 1.0/(n-1))
	#  The variance is incorrect, so I replaced it by a blank image, will fix later
	var.to_zero()
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
	# convert to real images
	var   = Util.pack_complex_to_real(var)
	sumsq = Util.pack_complex_to_real(sumsq)
	var = (var - sumsq/n)/(n-1)
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

def ssnr2d_ctf(data, mask = None, mode="", dopa=True):
	'''
	Calculate ssnr and variance in Fourier space for 2D images including CTF information
	If mode = "a" apply alignment parameters
	'''
	from fundamentals import fft, rot_shift2D, rot_avg_table, fftip
	from morphology   import ctf_img, threshold
	from filter       import filt_ctf
	from utilities    import get_params2D, pad
	import  types
	
	if type(data) is types.StringType:
		n = EMUtil.get_image_count(data)
		ima = EMData()
		ima.read_image(data, 0, True)
		if ima.get_attr_default('ctf_applied', 1) == 1:
			ERROR("data cannot be ctf-applied","ssnr2d",1)
		nx = ima.get_xsize()
		ny = ima.get_ysize()
	else:
		if data[0].get_attr_default('ctf_applied', 1) == 1:
			ERROR("data cannot be ctf-applied","ssnr2d",1)
		n = len(data)
		nx = data[0].get_xsize()
		ny = data[0].get_ysize()
	if  dopa:
		nx2 = nx*2
		ny2 = ny*2
	else:
		nx2 = nx
		ny2 = ny
	ctf_2_sum = EMData(nx2, ny2, 1, False)
	sumsq     = EMData(nx2, ny2, 1, False)

	for i in xrange(n):
		if type(data) is types.StringType:
			ima = EMData()
			ima.read_image(data, i)
		else:
			ima = data[i].copy()
		ctf_params = ima.get_attr('ctf')
		if mode == "a":
			alpha, sx, sy, mirror, scale = get_params2D(ima)
 			ima = rot_shift2D(ima, alpha, sx, sy, mirror)
		if mask:  Util.mul_img(ima, mask)
		if  dopa:  ima = pad(ima, nx2, ny2, 1, background = "circumference")
		fftip(ima)
		Util.add_img(sumsq, filt_ctf(ima, ctf_params, dopa))
		Util.add_img2(ctf_2_sum, ctf_img(nx2, ctf_params))
	print  "   NEW "

	ave = Util.divn_filter(sumsq, ctf_2_sum)

	var       = EMData(nx2, ny2, 1, False)
	for i in xrange(n):
		if type(data) is types.StringType:
			ima = EMData()
			ima.read_image(data, i)
		else:
			ima = data[i].copy()
		ctf_params = ima.get_attr('ctf')
		if mode == "a":
			alpha, sx, sy, mirror, scale = get_params2D(ima)
 			ima = rot_shift2D(ima, alpha, sx, sy, mirror)
		if mask:  Util.mul_img(ima, mask)
		if dopa:  ima = pad(ima, nx2, ny2, 1, background = "circumference")
		fftip(ima)
		"""
		ima = ima-filt_ctf(ave, ctf_params, dopa)
		Util.add_img2(var, ima)
		"""
		ima = filt_ctf(ima-filt_ctf(ave, ctf_params, dopa), ctf_params, dopa)
		#Util.div_filter(ima, ctf_2_sum)
		Util.add_img2(var, ima)

	"""
	Util.mul_scalar(var, 1.0/float(n-1))
	Util.mul_img(ave, ave.conjg())
	ave *= n
	"""
	Util.mul_img(sumsq, sumsq.conjg())
	#sumsq  = fft(window2d(fft(sumsq),nx,ny))

	#Util.div_filter(sumsq, ctf_2_sum)
	#Util.sub_img(var, sumsq)
	#Util.mul_scalar(var, 1.0/float(n-1))
	#Util.div_filter(sumsq, ctf_2_sum)

	#Util.mul_img(ave, ave.conjg())

	#Util.mul_img(ave, ctf_2_sum)
	#Util.mul_img(ave, ctf_2_sum)
	#ave *= n
	#var /= (n-1)

	nn = nx2//2
	nm = ny2//2

	from utilities import info
	from fundamentals import resample
	sumsq = Util.pack_complex_to_real(sumsq)
	sumsq[nn,nm] = sumsq[nn+1,nm]
	#if dopa:  sumsq = resample(sumsq,0.5)
	#info(var,None, "   tvar in ssnr2d")
	tvar = Util.divn_filter(var, ctf_2_sum)
	#info(tvar,None, "   tvar in ssnr2d first div")
	Util.div_filter(tvar, ctf_2_sum)
	#info(tvar,None, "   tvar in ssnr2d second div")
	tvar =  Util.pack_complex_to_real(tvar)
	#info(tvar,None, "   tvar in ssnr2d pack")
	tvar[nn,nm] = tvar[nn+1,nm]
	#if dopa:  tvar  = resample(tvar,0.5)
	#info(tvar,None, "   tvar in ssnr2d after resample")


	var = Util.pack_complex_to_real(var)
	var[nn,nm] = var[nn+1,nm]
	#info(var,None, "   var in ssnr2d")
	#if dopa:  var  = resample(var,0.5)
	#info(var,None, "   var in ssnr2d after resample")
	ssnr = sumsq/var - 1.0
	rave = rot_avg_table(sumsq)
	rvar = rot_avg_table(var)

	rssnr = []
	for i in xrange(len(rvar)):
		if rvar[i] > 0.0: qt = max(0.0, rave[i]/rvar[i] - 1.0)
		else:              ERROR("ssnr2d","rvar negative",1)
		rssnr.append(qt)
	return rssnr, rave, rvar, ssnr, sumsq, ave, tvar

def ssnr2d_ctf_OLD(data, mask = None, mode=""):
	'''
	Calculate ssnr and variance in Fourier space for 2D images including CTF information
	If mode = "a" apply alignment parameters
	'''
	from fundamentals import fft, fftip, rot_shift2D, rot_avg_table
	from morphology   import ctf_img, threshold
	from filter       import filt_ctf
	from utilities    import get_params2D
	import  types
	
	if type(data) is types.StringType:
		n = EMUtil.get_image_count(data)
		ima = EMData()
		ima.read_image(data, 0, True)
		if ima.get_attr_default('ctf_applied', 1) == 1:
			ERROR("data cannot be ctf-applied","ssnr2d",1)
		nx = ima.get_xsize()
		ny = ima.get_ysize()
		nz = ima.get_zsize()
	else:
		if data[0].get_attr_default('ctf_applied', 1) == 1:
			ERROR("data cannot be ctf-applied","ssnr2d",1)
		n = len(data)
		nx = data[0].get_xsize()
		ny = data[0].get_ysize()
		nz = data[0].get_zsize()

	ctf_2_sum = EMData(nx, ny, nz, False)
	sumsq     = EMData(nx, ny, nz, False)
	var       = EMData(nx, ny, nz, False)

	for i in xrange(n):
		if type(data) is types.StringType:
			ima = EMData()
			ima.read_image(data, i)
		else:
			ima = data[i].copy()
		ctf_params = ima.get_attr('ctf')
		if mode == "a":
			alpha, sx, sy, mirror, scale = get_params2D(ima)
 			ima = rot_shift2D(ima, alpha, sx, sy, mirror)
		if mask:  Util.mul_img(ima, mask)
		fftip(ima)
		oc = filt_ctf(ima, ctf_params)
		Util.add_img(sumsq, oc)
		Util.add_img2(var, ima)
		Util.add_img2(ctf_2_sum, ctf_img(nx, ctf_params, ny=ny))
	Util.mul_img(sumsq, sumsq.conjg())
	Util.div_filter(sumsq, ctf_2_sum)
	Util.sub_img(var, sumsq)
	Util.mul_scalar(var, 1.0/float(n-1))
	Util.div_filter(sumsq, ctf_2_sum)

	var   = Util.pack_complex_to_real(var)
	sumsq = Util.pack_complex_to_real(sumsq)
	sumsq *= n
	ssnr   = sumsq/var - 1.0
	rvar = rot_avg_table(var)
	rsumsq = rot_avg_table(sumsq)
	rssnr = []
	for i in xrange(len(rvar)):
		if rvar[i] > 0.0: qt = max(0.0, rsumsq[i]/rvar[i] - 1.0)
		else:              ERROR("ssnr2d","rvar negative",1)
		rssnr.append(qt)
	return rssnr, rsumsq, rvar, ssnr, sumsq, var

def varf(data, mask = None, mode="a"):
	'''
	Calculate variance in Fourier space for 2D or 3D images, (no CTF correction)
	If mode = "a" apply alignment parameters
	'''
	from fundamentals import fftip, rot_shift2D
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
		fftip(ima)
		Util.add_img(sumsq, ima)
		Util.add_img2(var, ima)

	Util.mul_img(sumsq, sumsq.conjg())
	Util.mad_scalar(var, sumsq, -1.0/float(n))
	Util.mul_scalar(var, 1.0/float(n-1))
	st = Util.infomask(var, None, True)
	if(st[2]<0.0):  ERROR("Negative variance!","varf",1)
	from fundamentals import rot_avg_table

	return var, rot_avg_table(Util.pack_complex_to_real(var))

def varfctf(data, mask = None, mode="a", dopad = True):
	'''
	Calculate variance in Fourier space for 2D or 3D images including ctf correction
	If mode = "a" apply alignment parameters
	This command is for ML average, i.e., A = sum_k (CTF_k F_k) / sum_k ( CTF_k^2 )
	'''
	from fundamentals import fftip, fft, rot_shift2D, window2d, cyclic_shift
	from morphology   import ctf_img
	from filter       import filt_ctf
	from utilities    import get_arb_params, get_params2D
	import  types

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
	if dopad:
		nx2 = 2*nx
		ny2 = 2*ny
		if( nz>1 ): nz2 = 2*nz
		else:       nz2 = nz
		from utilities import pad
	else:
		nx2 = nx
		ny2 = ny
		nz2 = nz
	ctf_2_sum = EMData(nx2, ny2, nz2, False)
	sumsq     = EMData(nx2, ny2, nz2, False)
	var       = EMData(nx2, ny2, nz, False)

	for i in xrange(n):
		if (type(data) is types.StringType):
			ima = EMData()
			ima.read_image(data, i)
		else:
			ima = data[i].copy()
		ctf_params = ima.get_attr("ctf")
		if(mode == "a"):
			alpha, sx, sy, mirror, scale = get_params2D(ima)
 			ima = rot_shift2D(ima, alpha, sx, sy, mirror)
		if(mask): Util.mul_img(ima, mask)
		if dopad:  ima = pad(ima, nx2, ny2, nz2, background = "circumference")
		fftip(ima)
		oc = filt_ctf(ima, ctf_params)
		Util.add_img(sumsq, oc)
		Util.add_img2(var, ima)
		Util.add_img2(ctf_2_sum, ctf_img(nx2, ctf_params, ny = ny2, nz = nz2))
	Util.mul_img(sumsq, sumsq.conjg())
	Util.div_filter(sumsq, ctf_2_sum)
	Util.sub_img(var, sumsq)
	Util.mul_scalar(var, 1.0/float(n-1))
	st = Util.infomask(var, None, True)
	if(st[2]<0.0):  ERROR("Negative variance!","varfctf",1)
	if dopad:  #  CHECK THIS< CAN IT BE DONE BETTER??
		var = fft( cyclic_shift(window2d(cyclic_shift(fft(var), nx, ny, nz), nx, ny), -nx//2, -ny//2, -nz//2) )

	from fundamentals import rot_avg_table

	return var, rot_avg_table(Util.pack_complex_to_real(var))

def varf2d(data, ave, mask = None, mode="a"):
	'''
	Calculate variance in Fourier space for 2D images including ctf correction
	ave is the Wiener average of data
	If mode = "a" apply alignment parameters
	This command is for Wiener average, i.e., A = sum_k (CTF_k SNR_k F_k) / [sum_k ( CTF_k^2 SNR_k) + 1]
	'''
	from fundamentals import fft, rot_shift2D
	from morphology   import ctf_img
	from filter       import filt_ctf
	from utilities    import get_arb_params, get_params2D, get_im
	from fundamentals import fft
	import  types

	if (type(data) is types.StringType):
		n = EMUtil.get_image_count(data)
	else:
		n = len(data)
	ima = get_im(data)
	nx = ima.get_xsize()
	ny = ima.get_ysize()
	nz = ima.get_zsize()
	if(ima.get_attr_default('ctf_applied', 1) == 1):
		ERROR("data cannot be ctf-applied","varf2d",1)
	if(nz > 1): ERROR("data cannot be 3D","varf2d",1)

	var = EMData(nx, ny, nz, False)

	defocus = -1.0
	for i in xrange(n):
		ima = get_im(data, i)
		ctf = ima.get_attr("ctf")
		if(mode == "a"):
			alpha, sx, sy, mirror, scale = get_params2D(ima)
 			ima = rot_shift2D(ima, alpha, sx, sy, mirror)
		if(mask): Util.mul_img(ima, mask)
		if(defocus != ctf.defocus):
			oc = filt_ctf(ave, ctf, dopad=True)
			defocus = ctf.defocus
			#print i, "  changed  ctf  ",defocus,Util.infomask(oc, None, True)
		Util.add_img2(var, fft(Util.subn_img(ima, oc)))

	Util.mul_scalar(var, 1.0/float(n-1))
	st = Util.infomask(var, None, True)
	if(st[2]<0.0):  ERROR("Negative variance!","varf2d",1)

	from fundamentals import rot_avg_table

	return var, rot_avg_table(Util.pack_complex_to_real(var))

def varf2d_MPI(myid, data, ave, mask = None, mode = "a", CTF = False, main_node = 0, comm = -1):
	"""
	Calculate variance in Fourier space for 2D images optionally including ctf correction
	ave is the Wiener average of data
	If mode = "a" apply alignment parameters
	This command is for Wiener average, i.e., A = sum_k (CTF_k SNR_k F_k) / [sum_k ( CTF_k^2 SNR_k) + 1]
	"""
	from utilities    import    model_blank, get_params2D
	from fundamentals import    rot_shift2D, fft, fftip
	from utilities    import    reduce_EMData_to_root
	from mpi          import    mpi_reduce, MPI_INT, MPI_SUM
	
	if comm == -1:
		from mpi import MPI_COMM_WORLD
		comm = MPI_COMM_WORLD

	n = len(data)
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	var  = EMData(nx, ny, 1, False)
	
	if CTF:
		from filter       import filt_ctf
		from morphology   import ctf_img
		if data[0].get_attr_default('ctf_applied', 1) == 1:
	 		ERROR("data cannot be ctf-applied", "add_ave_varf_MPI", 1)
		defocus = -1.0
	 	for i in xrange(n):
			ctf = data[i].get_attr("ctf")
			if(mode == "a"):
				alpha, sx, sy, mirror, scale = get_params2D(data[i])
 				ima = rot_shift2D(data[i], alpha, sx, sy, mirror)
			else:
				ima = data[i].copy()
			if(mask): Util.mul_img(ima, mask)
			if(defocus != ctf.defocus):
				oc = filt_ctf(ave, ctf, dopad=True)
				defocus = ctf.defocus
			Util.add_img2(var, fft(Util.subn_img(ima, oc)))
	else:
		for i in xrange(n):
			if mode == "a":
				alpha, sx, sy, mirror, scale = get_params2D(data[i])
				ima = rot_shift2D(data[i], alpha, sx, sy, mirror, scale, "quadratic")
				if mask:  Util.mul_img(ima, mask)
				fftip(ima)
			else:
				if  mask:   ima = fft(Util.muln_img(data[i], mask))
				else:       ima = fft(data[i])
			Util.add_img2(var, ima)
	reduce_EMData_to_root(var, myid, main_node, comm)
	nima = n
	nima = mpi_reduce(nima, 1, MPI_INT, MPI_SUM, main_node, comm)
	if myid == main_node:
		Util.mul_scalar(var, 1.0/float(nima-1))
		if var.get_value_at(0, 0) < 0.0:	var.set_value_at(0, 0, 0.0)		
		st = Util.infomask(var, None, True)
		if st[2] < 0.0:  ERROR("Negative variance!", "varf2_MPI", 1)
		from fundamentals import rot_avg_table
		return var, rot_avg_table(Util.pack_complex_to_real(var))
	else:
		return  EMData(), [0]  # return minimum what has to be returned, but no meaning.

def varf3d(prjlist,ssnr_text_file = None, mask2D = None, reference_structure = None, ou = -1, rw = 1.0, npad = 1, CTF = False, sign = 1, sym ="c1"):
	"""
	  Calculate variance in Fourier space of an object reconstructed from projections
	  
	  Known problems: properly speaking, the SSNR has to be calculated using snr=inf and this is what recons3d_nn_SSNR_MPI does.
	  So, when one computes reference structure, snr should be 1.0e20.  However, when the reference structure is passed
	  from the reconstruction program, it was computed using different snr.  I tested it and in practice there is no difference,
	  as this only changes the background variance due to reconstruction algorithm, which is much lower anyway.  PAP.
	"""
	from reconstruction import   recons3d_nn_SSNR, recons3d_4nn, recons3d_4nn_ctf
	from utilities      import   model_blank
	from projection     import   prep_vol, prgs

	[ssnr1, vol_ssnr1] = recons3d_nn_SSNR(prjlist, mask2D, rw, npad, sign, sym, CTF)

	nx  = prjlist[0].get_xsize()
	if ou == -1: radius = int(nx/2) - 1
	else:        radius = int(ou)
	if(reference_structure == None):
		if CTF :
			snr = 1.0#e20
			reference_structure = recons3d_4nn_ctf(prjlist, range(prjlist), snr, sign, sym, 0, npad)
		else  :
			reference_structure = recons3d_4nn(prjlist, range(prjlist), sym, npad)

	volft,kb = prep_vol(reference_structure)
	del reference_structure
	from utilities import get_params_proj
	if CTF: from filter import filt_ctf
	re_prjlist = []
	for prj in prjlist:
		phi,theta,psi,tx,ty = get_params_proj(prj)
		proj = prgs(volft, kb, [phi,theta,psi,-tx,-ty])
		if CTF:
			ctf_params = prj.get_attr("ctf")			
			proj = filt_ctf(proj, ctf_params)
			proj.set_attr('sign', 1)
		re_prjlist.append(proj)
	del volft
	[ssnr2, vol_ssnr2] = recons3d_nn_SSNR(re_prjlist, mask2D, rw, npad, sign, sym, CTF)

	outf = file(ssnr_text_file, "w")
	for i in xrange(len(ssnr2[0])):
		datstrings = []
		datstrings.append("  %15f" % ssnr1[0][i])    #  have to subtract 0.5 as in C code there is round.
		datstrings.append("  %15e" % ssnr1[1][i])    # SSNR
		datstrings.append("  %15e" % ssnr1[2][i])    # variance
		datstrings.append("  %15f" % ssnr1[3][i])    # number of points in the shell
		datstrings.append("  %15f" % ssnr1[4][i])    # number of added Fourier points
		datstrings.append("  %15e" % ssnr1[5][i])    # square of signal
		datstrings.append("  %15e" % ssnr2[1][i])    # SSNR
		datstrings.append("  %15e" % ssnr2[2][i])    # variance
		datstrings.append("  %15e" % ssnr2[5][i])    # square of signal
		datstrings.append("\n")
		outf.write("".join(datstrings))
	outf.close()
	vol_ssnr1 = Util.subn_img(Util.pack_complex_to_real(vol_ssnr1), Util.pack_complex_to_real(vol_ssnr2))
	del  vol_ssnr2
	nc = nx//2
	r2 = radius**2
	for i in xrange(nx):
		for j in xrange(nx):
			for k in xrange(nx):
				if( ( (i-nc)**2+(j-nc)**2+(k-nc)**2) <r2):
					if( vol_ssnr1.get_value_at(i,j,k) <= 0.0):
						bm = -1.0
						for i1 in xrange(-1,2):
							for i2 in xrange(-1,2):
								for i3 in xrange(-1,2):
									tm = vol_ssnr1.get_value_at(i+i1,j+i2,k+i3)
									if(tm > bm):
										bm = tm
						vol_ssnr1.set_value_at(i,j,k,bm)
									
				else:  vol_ssnr1.set_value_at(i,j,k, 1.0)
	return vol_ssnr1
	#from morphology import threshold_to_minval
	#return  threshold_to_minval(Util.subn_img(Util.pack_complex_to_real(vol_ssnr1), Util.pack_complex_to_real(vol_ssnr2)), 1.0)

def varf3d_MPI(prjlist, ssnr_text_file = None, mask2D = None, reference_structure = None, ou = -1, rw = 1.0, npad = 1, CTF = False, sign = 1, sym ="c1", myid = 0, mpi_comm = None):
	"""
	  Calculate variance in Fourier space of an object reconstructed from projections

	  Known problems: properly speaking, the SSNR has to be calculated using snr=inf and this is what recons3d_nn_SSNR_MPI does.
	  So, when one computes reference structure, snr should be 1.0e20.  However, when the reference structure is passed
	  from the reconstruction program, it was computed using different snr.  I tested it and in practice there is no difference,
	  as this only changes the background variance due to reconstruction algorithm, which is much lower anyway.  PAP.
	"""
	from reconstruction import   recons3d_nn_SSNR_MPI, recons3d_4nn_MPI, recons3d_4nn_ctf_MPI
	from utilities      import   bcast_EMData_to_all, model_blank
	from projection     import   prep_vol, prgs
	from mpi import MPI_COMM_WORLD
	
	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD
	
	if myid == 0: [ssnr1, vol_ssnr1] = recons3d_nn_SSNR_MPI(myid, prjlist, mask2D, rw, npad, sign, sym, CTF, mpi_comm=mpi_comm)
	else:                              recons3d_nn_SSNR_MPI(myid, prjlist, mask2D, rw, npad, sign, sym, CTF, mpi_comm=mpi_comm)

	nx  = prjlist[0].get_xsize()
	if ou == -1: radius = int(nx/2) - 2
	else:        radius = int(ou)
	if(reference_structure == None):
		if CTF :
			snr = 1.0#e20
			if myid == 0 :
				reference_structure = recons3d_4nn_ctf_MPI(myid, prjlist, snr, sign, sym, mpi_comm=mpi_comm)
			else :
				recons3d_4nn_ctf_MPI(myid, prjlist, snr, sign, sym, mpi_comm=mpi_comm)
				reference_structure = model_blank(nx, nx, nx)
		else  :
			if myid == 0 :
				reference_structure = recons3d_4nn_MPI(myid, prjlist, sym, mpi_comm=mpi_comm)
			else :
				recons3d_4nn_MPI(myid, prjlist, sym, mpi_comm=mpi_comm)
				reference_structure = model_blank(nx, nx, nx)
		bcast_EMData_to_all(reference_structure, myid, 0, mpi_comm)
	#if myid == 0:  reference_structure.write_image("refer.hdf",0)
	#vol *= model_circle(radius, nx, nx, nx)
	volft,kb = prep_vol(reference_structure)
	del reference_structure
	from utilities import get_params_proj
	if CTF: from filter import filt_ctf
	re_prjlist = []
	for prj in prjlist:
		phi,theta,psi,tx,ty = get_params_proj(prj)
		proj = prgs(volft, kb, [phi,theta,psi,-tx,-ty])
		if CTF:
			ctf_params = prj.get_attr("ctf")
			proj = filt_ctf(proj, ctf_params)
			proj.set_attr('sign', 1)
		re_prjlist.append(proj)
	del volft
	if myid == 0: [ssnr2, vol_ssnr2] = recons3d_nn_SSNR_MPI(myid, re_prjlist, mask2D, rw, npad, sign, sym, CTF, mpi_comm=mpi_comm)
	else:                              recons3d_nn_SSNR_MPI(myid, re_prjlist, mask2D, rw, npad, sign, sym, CTF, mpi_comm=mpi_comm)
	del re_prjlist

	if myid == 0 and ssnr_text_file != None:
		outf = file(ssnr_text_file, "w")
		for i in xrange(len(ssnr2[0])):
			datstrings = []
			datstrings.append("  %15f" % ssnr1[0][i])    #  have to subtract 0.5 as in C code there is round.
			datstrings.append("  %15e" % ssnr1[1][i])    # SSNR
			datstrings.append("  %15e" % ssnr1[2][i])    # variance
			datstrings.append("  %15f" % ssnr1[3][i])    # number of points in the shell
			datstrings.append("  %15f" % ssnr1[4][i])    # number of added Fourier points
			datstrings.append("  %15e" % ssnr1[5][i])    # square of signal
			datstrings.append("  %15e" % ssnr2[1][i])    # SSNR
			datstrings.append("  %15e" % ssnr2[2][i])    # variance
			datstrings.append("  %15e" % ssnr2[5][i])    # square of signal
			datstrings.append("\n")
			outf.write("".join(datstrings))
		outf.close()
	if myid == 0:
		vol_ssnr1 = Util.subn_img(Util.pack_complex_to_real(vol_ssnr1), Util.pack_complex_to_real(vol_ssnr2))
		del  vol_ssnr2
		# what follows is a risky business.  There should be a better way to deal with negative values. but for the time being...
		nc = nx//2
		r2 = radius**2
		for i in xrange(nx):
			for j in xrange(nx):
				for k in xrange(nx):
					if( ( (i-nc)**2+(j-nc)**2+(k-nc)**2) <r2):
						if( vol_ssnr1.get_value_at(i,j,k) <= 0.0):
							bm = -1.0
							for i1 in xrange(-1,2):
								for i2 in xrange(-1,2):
									for i3 in xrange(-1,2):
										tm = vol_ssnr1.get_value_at(i+i1,j+i2,k+i3)
										if(tm > bm):
											bm = tm
							vol_ssnr1.set_value_at(i,j,k,bm)
										
					else:  vol_ssnr1.set_value_at(i,j,k, 1.0)
		return  vol_ssnr1
		#from morphology import threshold_to_minval
		#return  threshold_to_minval( Util.subn_img(Util.pack_complex_to_real(vol_ssnr1), Util.pack_complex_to_real(vol_ssnr2)), 1.0)
	else:  return  model_blank(2,2,2)


def ccc(img1, img2, mask=None):
	"""Cross-correlation coefficient.	   
	   Usage: result = ccc(image1, image2 [, mask])
	"""
	return img1.cmp("ccc", img2, {"mask":mask,"negative":0})

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


def locres(vi, ui, m, nk, cutoff, step, myid, main_node, number_of_proc):
	from mpi 	  	  import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	  	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_gatherv, mpi_send, mpi_recv
	from mpi 	  	  import MPI_SUM, MPI_FLOAT, MPI_INT, MPI_TAG_UB
	from fundamentals import fft
	from utilities import model_blank, bcast_EMData_to_all, recv_EMData, send_EMData, bcast_number_to_all, info
	from filter import filt_tophatb
	from EMAN2 import rsconvolution
	from morphology import square_root
	

	nx = m.get_xsize()
	ny = m.get_ysize()
	nz = m.get_zsize()

	mc = model_blank(nx,ny,nz,1.0)-m

	if(myid == main_node):
		st = Util.infomask(vi,m,True)
		vi -= st[0]

		st = Util.infomask(ui,m,True)
		ui -= st[1]

	bcast_EMData_to_all(vi, myid, main_node)
	bcast_EMData_to_all(ui, myid, main_node)

	vf = fft(vi)
	uf = fft(ui)

	if(myid == 0):
		freqvol = model_blank(nx,ny,nz)
		resolut = []
	lp = int(max(nx,ny,nz)/2/step+0.5)
	step = 0.5/lp
	lt = lp//number_of_proc
	lp = (lt+1)*number_of_proc
	bailout = 0
	for i in xrange(myid,lp,number_of_proc):
		fl = step*i
		fh = fl+step
		freq=(fl+fh)/2.0
		#print number_of_proc,myid,lp,i,step,fl,fh,freq

		if i>0 :
			v = fft(filt_tophatb( vf, fl, fh))
			u = fft(filt_tophatb( uf, fl, fh))
			tmp1 = Util.muln_img(v,v)
			tmp2 = Util.muln_img(u,u)
			tmp3 = Util.muln_img(u,v)
			do = Util.infomask(square_root(Util.muln_img(tmp1,tmp2)),m,True)[0]
			dp = Util.infomask(tmp3,m,True)[0]
			#print "dpdo   ",myid,dp,do
			if do == 0.0: dis = [freq, 0.0]
			else:  dis = [freq, dp/do]
		else:
			tmp1 = model_blank(nx,ny,nz,1.0)
			tmp2 = model_blank(nx,ny,nz,1.0)
			tmp3 = model_blank(nx,ny,nz,1.0)
			dis = [freq, 1.0]


		tmp1 = Util.box_convolution(tmp1, nk)
		tmp2 = Util.box_convolution(tmp2, nk)
		tmp3 = Util.box_convolution(tmp3, nk)

		Util.mul_img(tmp1,tmp2)

		tmp1 = square_root(tmp1)

		Util.mul_img(tmp1,m)
		Util.add_img(tmp1,mc)


		Util.mul_img(tmp3,m)
		Util.add_img(tmp3,mc)

		Util.div_img(tmp3,tmp1)

		Util.mul_img(tmp3,m)

		mpi_barrier(MPI_COMM_WORLD)

		if(myid == main_node):
			for k in xrange(number_of_proc):
				if(k != main_node):
					#print " start receiving",myid,i
					tag_node = k+1001
					dis = mpi_recv(2, MPI_FLOAT, k, MPI_TAG_UB, MPI_COMM_WORLD)
					#print  "received ",myid,dis
					tmp3 = recv_EMData(k, tag_node)
					#print  "received ",myid
				if(dis[0] <=0.5):  resolut.append(dis)
				fl = step*(i+k)
				fh = fl+step
				freq=(fl+fh)/2.0
				#print k,dis,Util.infomask(tmp3,m,True)
				#if(k == number_of_proc-1):  bailout = 1
				bailout = 0
				#print  "setting freqvol  ",k
				Util.set_freq(freqvol,tmp3,m,cutoff,freq)
				"""
				for x in xrange(nx):
					for y in xrange(ny):
						for z in xrange(nz):
							if(m.get_value_at(x,y,z) > 0.5):
								if(freqvol.get_value_at(x,y,z) == 0.0):
									if(tmp3.get_value_at(x,y,z) < cutoff):
										freqvol.set_value_at(x,y,z,freq)
										bailout = 0
									else:
										if(k == number_of_proc-1):
											bailout = 0
				"""
				#print k,freq,Util.infomask(freqvol,m,True)

		else:
			tag_node = myid+1001
			#print   "sent from", myid,dis
			mpi_send(dis, 2, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)
			#print   "sending EMD from", myid
			send_EMData(tmp3, main_node, tag_node)
			#print   "sent EMD from",myid

		bailout = bcast_number_to_all(bailout, main_node)
		if(bailout == 1):  break

	mpi_barrier(MPI_COMM_WORLD)
	if( myid == main_node ):  return freqvol, resolut
	else:  return None, None


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
				drop_image(tavg,"tavg.spi")
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
		calculate 1D rotationally averaged power spectrum of images in stack file
		Input
			stack
		Output
			a list containing 1D rotationally averaged power spectrum
	"""
	from utilities import get_im
	from EMAN2 import periodogram
	if  type(stack) == type(""): nima = EMUtil.get_image_count(stack)
	else:                       nima = len(stack)
	for i in xrange(nima):
		img = get_im(stack,i)
		e  = periodogram(img)
		ro = e.rotavg()
		if(i==0): rosum = ro.copy()
		else: Util.add_img(rosum, ro)
	rosum /= nima
	nr = rosum.get_xsize()
	table = [0.0]*nr
	for ir in xrange(nr):  table[ir] = rosum.get_value_at(ir)
	return table

def histogram(image, mask = None, nbins = 0, hmin = 0.0, hmax = 0.0):
	"""
		Name
			histogram - calculate a histogram of the image pixel values
		Input
			input: input image
			mask: optional binary mask
			nbins:number of bins. Optional. 
			hmin, hmax: Optional. 
		Output
			h:a list containining 2*nbins elements.
	"""
	return  Util.histogram(image, mask, nbins, hmin, hmax)

def im_diff(im1, im2, mask = None):
	import types
	from utilities import model_circle, get_im
	if type(im1) == types.StringType : im1 = get_im(im1)
	if type(im2) == types.StringType : im2 = get_im(im2)
	nx = im1.get_xsize()
	ny = im1.get_ysize()
	nz = im1.get_zsize()
	if mask != None :
		if   type(mask) == types.FloatType or type(mask) == types.IntType: m = model_circle(mask, nx, ny, nz)
		elif type(mask) == types.StringType:   m = get_im(mask)
		else: m = mask
	else:
		if   im1.get_ndim() == 3: radius = min(nx,ny,nz)//2 - 1
		elif im1.get_ndim() == 2: radius = min(nx,ny)//2    - 1
		else:                     radius = int(nx)//2       - 1
		m = model_circle(radius, nx, ny, nz)
	l = Util.im_diff(im1, im2, m)
	return  l["imdiff"], l["A"], l["B"]


##############################################################################################
### K-MEANS ##################################################################################

# init cluster assignment randomly
def k_means_init_asg_rnd(N, K):
	from random import randint
	assign  = [0] * N
	nc      = [0] * K
	retrial = 20
	while retrial > 0:
		retrial -= 1
		i = 0
		for im in xrange(N):
			assign[im] = randint(0, K-1)
			nc[assign[im]] += 1
		flag,k = 1,K
		while k>0 and flag:
			k -= 1 
			if nc[k] <= 1:
				flag = 0
				if retrial == 0: ERROR('Empty class in the initialization', 'k_means_SSE', 1)
				for k in xrange(K): nc[k] = 0
		if flag == 1:
			retrial = 0

	return assign, nc

# init cluster assignment by D2 weighting
def k_means_init_asg_d2w(im, N, K):
	from random    import randrange, gauss
	from numpy     import ones
	#from utilities import print_msg
	from sys       import exit
	
	C = [im[randrange(N)]]
	#print_msg('\n')
	d = ones((N)) * 1e4
	for k in xrange(K-1):
		#print_msg('\rInitialization with D2 weighting method... %i / %i' % (k+2, K))
		for n in xrange(N):
			val = im[n].cmp("SqEuclidean", C[k])
			if val < d[n]: d[n] = val

		SD2  = d.sum()
		p    = d / float(SD2)
		p    = map(float, p)
		ps   = zip(p, range(N))
		ps.sort(reverse = True)
		ind  = gauss(0, n // 6) # 6 is a cst define empirically
		ind  = int(abs(ind))
		C.append(im[ps[ind][1]])

	#print_msg('\n')
	assign = [0] * N
	nc     = [0] * K
	for n in xrange(N):
		res            = Util.min_dist_real(im[n], C)
		assign[n]      = res['pos']
		nc[assign[n]] += 1

	#print nc

	return assign, nc

# Convert local assignment to absolute assignment
def k_means_locasg2glbasg(ASG, LUT, N):
	Nloc = len(ASG)
	GASG = [-1] * N
	for n in xrange(Nloc): GASG[LUT[n]] = ASG[n]

	return GASG

# # # horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
# # k-means open and prepare images, only unstable objects (active = 1)
# def k_means_list_active(stack):
# 	from utilities     import file_type
# 	from EMAN2db import db_open_dict
# 	
# 	N    = EMUtil.get_image_count(stack)
# 
# 	image = EMData()
# 	image.read_image(stack, 0, True)
# 	flagactive = True
# 	try:	active = image.get_attr('active')
# 	except: flagactive = False
# 	del image
# 
# 	if flagactive:
# 		LUT  = []
# 		ext  = file_type(stack)
# 		if ext == 'bdb':
# 			DB = db_open_dict(stack)
# 			for n in xrange(N):
# 				if DB.get_attr(n, 'active'): LUT.append(n)
# 			DB.close()
# 		else:
# 			im = EMData()
# 			for n in xrange(N):
# 				im.read_image(stack, n, True)
# 				if im.get_attr('active'): LUT.append(n)
# 
# 		N = len(LUT)
# 	else:
# 		LUT = range(N)
# 
# 	return LUT, N

# k-means, prepare to open images later
def k_means_init_open_im(stack, maskname):
	from utilities import get_image, get_im, model_blank, file_type

	ext = file_type(stack)
	if ext == 'txt': TXT = True
	else:            TXT = False

	# open mask if defined
	if maskname != None:
		mask = get_image(maskname)
		im   = Util.compress_image_mask(mask, mask)
		m    = im.get_xsize()
		del im
	else:
		mask = None
		if TXT:
			line = open(stack, 'r').readline()
			m    = len(line.split())
		else:
			im = get_im(stack, 0)
			m  = im.get_xsize() * im.get_ysize() * im.get_zsize()
			del im

	# get some params
	if TXT:
		Ntot = len(open(stack, 'r').readlines())
	else:   Ntot = EMUtil.get_image_count(stack)

	# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1	
	# # check if the flag active is used, in the case where k-means will run for the stability
	# if TXT:
	# 	flagactive = False
	# else:
	# 	image = EMData()
	# 	image.read_image(stack, 0, True)
	# 	flagactive = True
	# 	try:	active = image.get_attr('active')
	# 	except: flagactive = False
	# 	del image
	# 
	# # if flag active used, prepare the list of images
	# if flagactive:
	# 	active = EMUtil.get_all_attributes(stack, "active")
	# 	lim  = []
	# 	for n in xrange(len(active)):
	# 		if active[n]: lim.append(n)
	# 	"""
	# 	from utilities import write_text_file
	# 	from sys import exit
	# 	write_text_file(active,'active')
	# 	write_text_file(lim,'lim')
	# 	exit()
	# 	"""
	# 	del active
	# 	N = len(lim)
	# else:
	# 	lim = range(Ntot)
	# 	N = Ntot

	lim = range(Ntot)
	N = Ntot

	return lim, mask, N, m, Ntot

# k-means open and prepare images
def k_means_open_im(stack, mask, CTF, lim, flagnorm = False):
	from utilities     import get_params2D, get_image, get_params3D, file_type, model_blank, print_msg
	from fundamentals  import rot_shift2D, rot_shift3D
	from sys           import exit
	if CTF:
		from morphology		import ctf_2, ctf_1d
		from filter		import filt_ctf, filt_table
		from fundamentals 	import fftip
		from utilities          import get_arb_params

	ext = file_type(stack)
	if ext == 'txt': TXT = True
	else:            TXT = False
	N   = len(lim)

	# to manage coord fact in text file format
	if TXT:
		IM   = [None] * N
		c    = 0
		data = open(stack, 'r').readlines()
		nx   = len(data[0].split())
		for idi in lim:
			line = data[idi]
			im   = model_blank(nx)
			line = line.split()
			for i in xrange(nx):
				val = float(line[i])
				im.set_value_at_fast(i, 0, val)

			if mask != None: im = Util.compress_image_mask(im, mask)

			IM[c] = im.copy()
			c += 1

	else:
	        im = EMData()
	        im.read_image(stack, 0)
	        nx = im.get_xsize()
	        ny = im.get_ysize()
	        nz = im.get_zsize()
	        if CTF:
	        	ctf	    = [[] for i in xrange(N)]
	        	ctf2	    = [[] for i in xrange(N)]
	        	ctf_params  = im.get_attr( "ctf" )
	        	if im.get_attr("ctf_applied")>0.0: ERROR('K-means cannot be performed on CTF-applied images', 'k_means', 1)

	        IM = im.read_images(stack, lim)
	        for i in xrange(N):
	        	# 3D object
	        	if nz > 1:
	        		try:
	        			phi, theta, psi, s3x, s3y, s3z, mirror, scale = get_params3D(IM[i])
	        			IM[i]  = rot_shift3D(IM[i], phi, theta, psi, s3x, s3y, s3z, scale)
	        			if mirror: IM[i].process_inplace('xform.mirror', {'axis':'x'})
	        		except:
	        			#ERROR('K-MEANS no 3D alignment parameters found', "k_means_open_im", 1)
	        			#sys.exit()
					pass
	        	# 2D object
	        	elif ny > 1:
	        		try:
	        			alpha, sx, sy, mirror, scale = get_params2D(IM[i])
	        			IM[i] = rot_shift2D(IM[i], alpha, sx, sy, mirror, scale)
	        		except: 
	        			#ERROR('K-MEANS no 2D alignment parameters found', "k_means_open_im", 1)
	        			#sys.exit()
					pass

	        	# obtain ctf
	        	if CTF:
	        		ctf_params = IM[i].get_attr( "ctf" )
	        		ctf[i]     = ctf_1d(nx, ctf_params)
	        		ctf2[i]    = ctf_2(nx, ctf_params)

	        	if flagnorm:
	        		# normalize
	        		ave, std, mi, mx = Util.infomask(IM[i], mask, True)
	        		IM[i] -= ave
	        		IM[i] /= std

	        	# apply mask
	        	if mask != None:
	        		if CTF: Util.mul_img(IM[i], mask)
	        		else:	IM[i] = Util.compress_image_mask(IM[i], mask)

	        	# fft
	        	if CTF: fftip(IM[i])

	        	# mem the original size
	        	if i == 0:
	        		IM[i].set_attr('or_nx', nx)
	        		IM[i].set_attr('or_ny', ny)
	        		IM[i].set_attr('or_nz', nz)

	if CTF: return IM, ctf, ctf2
	else:   return IM, None, None

# k-means write the head of the logfile
def k_means_headlog(stackname, outname, method, N, K, crit, maskname, trials, maxit, CTF, T0, F, rnd, ncpu, m, init_method='random'):
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
	print_msg('Number of pixels under mask : %i\n'     % m)
	print_msg('Number of clusters          : %s\n'     % Ks)
	print_msg('Number of trials            : %i\n'     % trials)
	print_msg('Maximum iteration           : %i\n'     % maxit)
	print_msg('Data with CTF               : %s\n'     % CTF)
	print_msg('Criterion                   : %s\n'     % crit)
	print_msg('Optimization method         : %s\n'     % method)
	if SA:
		print_msg('Simulated annealing          : ON\n')
		#if SA2: print_msg('   select neighbour         : closer according T\n')
		#else: 	 print_msg('   select neighbour         : randomly\n')
		print_msg('   T0                       : %f\n' % T0)
		print_msg('   F                        : %f\n' % F)
	else:
		print_msg('Simulated annealing          : OFF\n')
	print_msg('Random seed                 : %i\n'     % rnd)
	print_msg('Initialization method       : %s\n'     % init_method)
	print_msg('Number of CPUs              : %i\n'     % ncpu)
	print_msg('Output seed names           : %s\n\n'   % outname)

# K-means write results output directory
def k_means_export(Cls, crit, assign, out_seedname, part = -1, TXT = False):
	from utilities import print_msg
	import os

	if not os.path.exists(out_seedname): os.mkdir(out_seedname)

	# write the report on the logfile
	Je = 0
	flagHDF = False
	for k in xrange(Cls['k']):
		Je += Cls['Ji'][k]
		if Cls['n'][k] > 16000:
			flagHDF = True
			print_msg('\nWARNING: limitation of number attributes in hdf format, the results will be export in separate text files\n')

	print_msg('\n\n_Details____________________________________________________\n')
	print_msg('\n\t%s\t%11.6e\n\n' % ('The total Sum of Squares Error (Je) = ', Je))

	for name in crit['name']:
		if name   == 'C': print_msg('\t%s\t%11.4e\n' % ('Criteria Coleman', crit['C']))
		elif name == 'H': print_msg('\t%s\t%11.4e\n' % ('Criteria Harabasz', crit['H']))
		elif name == 'D': print_msg('\t%s\t%11.4e\n' % ('Criteria Davies-Bouldin', crit['D']))
		else:             ERROR('Kind of criterion k-means unknown', 'k_means_out_res', 0)	
	print_msg('\n')

	for k in xrange(Cls['k']):
		print_msg('\t%s\t%d\t%s\t%d' % ('Cluster no:', k, 'No of Objects = ', Cls['n'][k]))
		if(Cls['n'][k] > 1): print_msg('\t%s\t%11.6e\t%s\t%11.6e\n' % ('Sum of Squares Error Ji', Cls['Ji'][k], ' Variance', Cls['Ji'][k] / float(Cls['n'][k]-1)))
		else:               print_msg('\t%s\t%11.6e\n' % ('Sum of Squares Error Ji', Cls['Ji'][k]))

		lassign = []
		for i in xrange(len(assign)):
			if(assign[i] == k):  lassign.append(float(i))

		# limitation of hdf file in the numbers of attributes
		if flagHDF or TXT:
			
			if part != -1:
				outfile = open(os.path.join(out_seedname, 'kmeans_part_%02i_grp_%03i.txt' % (part + 1, k + 1)), 'w')
			else:
				outfile = open(os.path.join(out_seedname, 'kmeans_grp_%03i.txt' % (k + 1)), 'w')
			list_images = []
			for i in xrange(len(assign)):
				if assign[i] == k:
					list_images.append(i)
					outfile.write(str(i) +'\n')
			outfile.close()
			Cls['ave'][k].set_attr_dict({'Class_average':1.0, 'nobjects':float(Cls['n'][k])})
		else:
			
			Cls['ave'][k].set_attr('Class_average', 1.0)
			Cls['ave'][k].set_attr('nobjects', float(Cls['n'][k]))
			Cls['ave'][k].set_attr('members', lassign)
			Cls['ave'][k].set_attr('Ji', float(Cls['Ji'][k]) )
			Cls['ave'][k].set_attr('Je', Je)

			Cls['var'][k].set_attr('Class_variance', 1.0)
			Cls['var'][k].set_attr('nobjects', float(Cls['n'][k]))
			Cls['var'][k].set_attr('members', lassign)
			Cls['var'][k].set_attr('Ji', float( Cls['Ji'][k] ) )
			Cls['var'][k].set_attr('Je', Je)

		if part == -1:
			Cls['ave'][k].write_image(os.path.join(out_seedname, "averages.hdf"), k)
			Cls['var'][k].write_image(os.path.join(out_seedname, "variances.hdf"), k)
		else:
			Cls['ave'][k].write_image(os.path.join(out_seedname, "averages_%02i.hdf" % part), k)
			Cls['var'][k].write_image(os.path.join(out_seedname, "variances_%02i.hdf" % part), k)

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
	Crit['C']    = 0.0
	Crit['H']    = 0.0
	Crit['D']    = 0.0
			
	ch = True
	try:
		name = Crit['name'].index("C")
	except:
		try:
			name = Crit['name'].index("H")
		except:
			ch = False
	if ch:
		Tr, Je = 0.0, 0.0
		buf.to_zero()
		nob = 0
		for k in xrange(Cls['k']):
			Util.add_img(buf, Cls['ave'][k]*Cls['n'][k])
			nob += Cls['n'][k]
			Je += Cls['Ji'][k]
		Util.mul_scalar(buf, 1.0/float(nob))
		for k in xrange(Cls['k']):	Tr += buf.cmp("SqEuclidean", Cls['ave'][k])
	
	# Compute the criterion required
	for name in Crit['name']:
	
		# Coleman criterion C = Tr(Cls) * S Tr(Cls(im)) -> C = S (ave-m_ave)**2  *  Je
		if name == 'C':
			Crit['C'] = Tr * Je

		# Harabasz criterion H = [Tr(Cls)/(K-1)] / [S Tr(Cls(im))/(N-K)]
		elif name == 'H':

			if Je > 0.0:
				Crit['H'] = (Tr / (Cls['k'] - 1)) / (Je / (N - Cls['k']))
			else:
				Crit['H'] = 1.0e20

		# Davies-Bouldin criterion DB = 1/K S_i[max all(i>j) ((Ji/n) + (Jj/n))/(ave_j-ave_i)**2]
		elif name == 'D':
			db, DB, err, ji, jj = 0.0, 0.0, 0.0, 0, 0
			for i in xrange(Cls['k']):
				val_max = 0.0
				for j in xrange(Cls['k']):
					if i != j:		
						err = Cls['ave'][j].cmp("SqEuclidean",Cls['ave'][i])
						if err > 0.0:
							if Cls['n'][i] > 0:	ji = Cls['Ji'][i] / Cls['n'][i]
							else:			ji = 0
							if Cls['n'][j] > 0:	jj = Cls['Ji'][j] / Cls['n'][j]
							else:			jj = 0
							db = (ji + jj) / err
						else:
							db = 1.0e20													
					if db > val_max:	val_max = db					
				DB += val_max
			Crit['D'] = DB / Cls['k']			

		else:
			ERROR("Criterion type for K-means unknown","k_means_criterion",1)
	# return the results
	return Crit

# K-means SA selection, modification of the function alignment.py/select_K which didn't work
# very well for k-means, error of selection, insufficiently smooth, ... JB 2009-01-16 12:30:43
def select_kmeans(dJe, T):
	from random import random
	from math   import exp

	K    = len(dJe)
	p    = [[0.0, k] for k in xrange(K)]
	pw   = 1.0 / T
	sump = 0.0
	for k in xrange(K):
	    arg    = -dJe[k] * pw
	    arg    = min( max(arg, -30.0), 30.0)
	    p[k][0] = exp(arg)
	    sump  += p[k][0]

	for k in xrange(K): p[k][0] /= sump
	p.sort()
	for k in xrange(1, K): p[k][0] += p[k - 1][0]
	rnd = random()
	s   = 0
	while p[s][0] < rnd: s += 1
	s   = p[s][1]

	return s

# K-means with classical method
def k_means_cla(im_M, mask, K, rand_seed, maxit, trials, CTF, F=0, T0=0, DEBUG=False, rnd_method = 'rnd'):
	from utilities 		import model_blank, get_im, running_time
	from random    		import seed, randint
	from utilities 		import print_msg
	from copy		import deepcopy
	import sys
	import time
	if CTF[0]:
		from filter	        import filt_ctf, filt_table
		from fundamentals 	import fftip

		ctf  = deepcopy(CTF[1])
		ctf2 = deepcopy(CTF[2])
		CTF  = True
	else:
		CTF  = False 

	# Simulated annealing use or not
	if F != 0: SA = True
	else:      SA = False

	if SA:
		# for simulated annealing
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
		
	flag_empty = False
	ntrials    = 0
	wd_trials  = 0
	SA_run     = SA
	
	ALL_EMPTY = True
	while ntrials < trials:
		ntrials  += 1

		# for simulated annealing
		SA = SA_run
		if SA: T = T0

		# Init the cluster by an image empty
		buf.to_zero()
		for k in xrange(K):
			Cls['ave'][k] = buf.copy()
			Cls['var'][k] = buf.copy()
			Cls['n'][k]   = 0
			Cls['Ji'][k]  = 0

		if rnd_method == 'd2w': assign, Cls['n'] = k_means_init_asg_d2w(im_M, N, K)
		else:                   assign, Cls['n'] = k_means_init_asg_rnd(N, K)

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
		print_msg('Criterion: %11.6e \n' % Je)

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
					res = Util.min_dist_four(im_M[im], CTFxAVE)
				else:
					res = Util.min_dist_real(im_M[im], Cls['ave'])

				# Simulated annealing
				if SA:
					dJe = [0.0] * K
					ni  = float(Cls['n'][assign[im]])
					di  = res['dist'][assign[im]]
					for k in xrange(K):
						if k != assign[im]:
							nj     = float(Cls['n'][k])
							dj     = res['dist'][k]
							dJe[k] = (ni/(ni-1))*(di/norm) - (nj/(nj+1))*(dj/norm)
						else:
							dJe[k] = 0

					# normalize and select
					mindJe = min(dJe)
					scale  = max(dJe) - mindJe
					for k in xrange(K): dJe[k] = 1 - (dJe[k] - mindJe) / scale
					select = select_kmeans(dJe, T)

					if select != res['pos']:
						ct_pert    += 1
						res['pos']  = select
					
				# update assign
				if res['pos'] != assign[im]:
					Cls['n'][assign[im]] -= 1
					if Cls['n'][assign[im]] < 1: Cls['n'][assign[im]] = 0
					assign[im]            = res['pos']
					Cls['n'][assign[im]] += 1
					change                = True
				
			# manage empty cluster
			for k in xrange(K):
				if Cls['n'][k] <= 1:
					print_msg('>>> WARNING: Empty cluster, restart with new partition.\n\n')
					flag_empty = True
					break
			if flag_empty: break
													
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
			if Je != 0: thd = abs(Je - old_Je) / Je
			else:	    thd = 0

			# Simulated annealing, update temperature
			if SA:
				if thd < 1.0e-12 and ct_pert == 0: watch_dog = maxit
				T *= F
				if T < 0.009 and ct_pert < 5: SA = False
				print_msg('> iteration: %5d    criterion: %11.6e    T: %13.8f  ct disturb: %5d\n' % (ite, Je, T, ct_pert))
				if DEBUG: print '> iteration: %5d    criterion: %11.6e    T: %13.8f  ct disturb: %5d' % (ite, Je, T, ct_pert)
			else:
				if thd < 1.0e-8: watch_dog = maxit
				print_msg('> iteration: %5d    criterion: %11.6e\n' % (ite, Je))
				if DEBUG: print '> iteration: %5d    criterion: %11.6e' % (ite, Je)

			old_Je = Je

		if not flag_empty:
			# memorize the result for this trial	
			if trials > 1:
				MemCls[ntrials-1]    = deepcopy(Cls)
				MemJe[ntrials-1]     = deepcopy(Je)
				MemAssign[ntrials-1] = deepcopy(assign)
				print_msg('# Criterion: %11.6e \n' % Je)
				ALL_EMPTY = False
			# set to zero watch dog trials
			wd_trials = 0
		else:
			flag_empty  = False
			wd_trials  += 1
			if wd_trials > 10:
				if trials > 1:
					MemJe[ntrials-1] = 1e10
					if ntrials == trials:
						#print_msg('>>> WARNING: After ran 10 times with different partitions, one cluster is still empty, STOP k-means.\n\n')
						#sys.exit()
						print_msg('>>> WARNING: After ran 10 times with different partitions, one cluster is still empty. \n\n')
						
					else:	print_msg('>>> WARNING: After ran 10 times with different partitions, one cluster is still empty, start the next trial.\n\n')
				else:
					print_msg('>>> WARNING: After ran 10 times with different partitions, one cluster is still empty, STOP k-means.\n\n')
					sys.exit()
				wd_trials = 0
			else:
				ntrials -= 1
	
	# if all trials resulted in empty cluster, exit!
	
	if trials > 1:
		if ALL_EMPTY:
			print_msg('>>> WARNING: All trials resulted in empty clusters, STOP k-means.\n\n')
			sys.exit()
						
	# if severals trials choose the best
	if trials > 1:
		val_min = 1.0e20
		best    = -1
		for n in xrange(trials):
			if MemJe[n] < val_min:
				val_min = MemJe[n]
				best    = n
		
		# affect the besta
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
			Cls['ave'][k].do_ift_inplace()
			Cls['var'][k].do_ift_inplace()
			Cls['ave'][k].depad()
			Cls['var'][k].depad()

	# information display
	running_time(t_start)
	print_msg('Criterion = %11.6e \n' % Je)
	for k in xrange(K):	print_msg('Cls[%i]: %i\n'%(k, Cls['n'][k]))

	# to debug
	if DEBUG: print Cls['n']
	
	# return Cls, assign
	return Cls, assign


# K-means with SSE method
def k_means_SSE(im_M, mask, K, rand_seed, maxit, trials, CTF, F=0, T0=0, DEBUG=False, rnd_method = 'rnd'):
	from utilities    import model_blank, get_im, running_time
	from utilities    import print_begin_msg, print_end_msg, print_msg
	from random       import seed, randint, shuffle
	from copy         import deepcopy
	import sys
	import time
	
	if CTF[0]:
		from filter		import filt_ctf, filt_table
		from fundamentals	import fftip

		ctf  = deepcopy(CTF[1])
		ctf2 = deepcopy(CTF[2])
		CTF  = True
	else:
		CTF  = False

	# Simulated annealing use or not
	if T0 != 0: SA = True
	else:       SA = False

	if SA:
		# for simulated annealing
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
	flag_empty = False	
	ntrials    = 0
	wd_trials  = 0
	SA_run     = SA
	
	ALL_EMPTY = True	
	while ntrials < trials:
		ntrials += 1

		# for simulated annealing
		SA = SA_run
		if SA: T = T0
	
		# Init the cluster by an image empty
		buf.to_zero()
		for k in xrange(K):
			Cls['ave'][k] = buf.copy()
			Cls['var'][k] = buf.copy()
			Cls['n'][k]   = 0
			Cls['Ji'][k]  = 0

		if rnd_method == 'd2w': assign, Cls['n'] = k_means_init_asg_d2w(im_M, N, K)
		else:	                assign, Cls['n'] = k_means_init_asg_rnd(N, K)
		
		
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
		print_msg('Criterion: %11.6e \n' % Je)

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
					res = Util.min_dist_four(im_M[im], CTFxAve)
				else:
					# compute the minimum distance with centroids
					res = Util.min_dist_real(im_M[im], Cls['ave'])

				dJe = [0.0] * K
				ni  = float(Cls['n'][assign[im]])
				di  = res['dist'][assign[im]]
				for k in xrange(K):
					if k != assign[im]:
						nj  = float(Cls['n'][k])
						dj  = res['dist'][k]
						dJe[k] =  (ni/(ni-1))*(di/norm) - (nj/(nj+1))*(dj/norm)
					else:
						dJe[k] = 0	
				# Simulate Annealing
				if SA:
					
					
					# normalize and select
					mindJe = min(dJe)
					scale  = max(dJe) - mindJe
					for k in xrange(K): dJe[k] = 1 - (dJe[k] - mindJe) / scale
					select = select_kmeans(dJe, T)
					
					if select != res['pos']:
						ct_pert    += 1
						res['pos']  = select
				else:
					max_value = -1.e30
					for i in xrange( len(dJe) ):
						if( dJe[i] >= max_value) :
							max_value = dJe[i]
							res['pos'] = i

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
					
					if Cls['n'][assign_from] <= 1:
						print_msg('>>> WARNING: Empty cluster, restart with new partition %d.\n\n' % wd_trials)
						flag_empty = True
												
					change = True

				# empty cluster
				if flag_empty: break

			# empty cluster
			if flag_empty: break
			
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
			if Je != 0: thd = abs(Je - old_Je) / Je
			else:       thd = 0

			# Simulated annealing, update temperature
			if SA:
				if thd < 1e-12 and ct_pert == 0: watch_dog = maxit
				T *= F
				if T < 0.009: SA = False
				print_msg('> iteration: %5d    criterion: %11.6e    T: %13.8f  ct disturb: %5d\n' % (ite, Je, T, ct_pert))
				if DEBUG: print '> iteration: %5d    criterion: %11.6e    T: %13.8f  ct disturb: %5d' % (ite, Je, T, ct_pert)
			else:
				if thd < 1e-8: watch_dog = maxit
				print_msg('> iteration: %5d    criterion: %11.6e\n'%(ite, Je))
				if DEBUG: print '> iteration: %5d    criterion: %11.6e'%(ite, Je)

			old_Je = Je

		# if no empty cluster
		if not flag_empty:

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
				print_msg('# Criterion: %11.6e \n' % Je)
				ALL_EMPTY = False
			# set to zero watch dog trials
			wd_trials = 0
		else:
			flag_empty  = False
			wd_trials  += 1
			if wd_trials > 10:
				
				if trials > 1:
					MemJe[ntrials-1] = 1e10
					if ntrials == trials:
						#print_msg('>>> WARNING: After ran 10 times with different partitions, one cluster is still empty, STOP k-means.\n\n')
						#sys.exit()
						print_msg('>>> WARNING: After ran 10 times with different partitions, one cluster is still empty. \n\n')
					else:	print_msg('>>> WARNING: After ran 10 times with different partitions, one cluster is still empty, start the next trial.\n\n')
				else:
					print_msg('>>> WARNING: After ran 10 times with different partitions, one cluster is still empty, STOP k-means.\n\n')
					sys.exit()
				wd_trials = 0
			else:
				ntrials -= 1

	if trials > 1:
		if ALL_EMPTY:
			print_msg('>>> WARNING: All trials resulted in empty clusters, STOP k-means.\n\n')
			sys.exit()

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
			Cls['ave'][k].do_ift_inplace()
			Cls['var'][k].do_ift_inplace()
			Cls['ave'][k].depad()
			Cls['var'][k].depad()

	# information display
	running_time(t_start)
	print_msg('Criterion = %11.6e \n' % Je)
	for k in xrange(K):	print_msg('Cls[%i]: %i\n'%(k, Cls['n'][k]))
	
	# to debug
	if DEBUG: print Cls['n']
		
	# return Cls, assign and Je
	return Cls, assign

def k_means_SSE_combine(Cls, assign, Je, N, K, ncpu, myid, main_node):
	# Common
	from utilities   import print_begin_msg, print_end_msg, print_msg, file_type, running_time
	from statistics  import k_means_locasg2glbasg
	from time        import time
	import sys, os
	
	from mpi        import mpi_init, mpi_comm_size, mpi_comm_rank, mpi_barrier
	from mpi        import MPI_COMM_WORLD, MPI_INT, mpi_bcast
	from mpi	import MPI_FLOAT, MPI_INT, mpi_recv, mpi_send, MPI_TAG_UB
	from utilities  import bcast_number_to_all, recv_EMData, send_EMData

	#print "my id ==", myid, " assign [10:20] ", assign[10:20], " Je===", Je, "Cls==",Cls[ 'n' ], " Ji==", Cls['Ji']

	if myid == main_node:
		je_return = [0.0]*(ncpu)
		for n1 in xrange(ncpu):
			if n1 != main_node: je_return[n1]	=	mpi_recv(1,MPI_FLOAT, n1, MPI_TAG_UB, MPI_COMM_WORLD)
 			else:               je_return[main_node]  = Je
	else:
		mpi_send(Je, 1, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)
	n_best = -1
	if(myid == main_node):
		je_return = map(float, je_return)
		#print "recived je", je_return
		min_je= 1.0e30

		#for i in xrange(ncpu):
			#print "i ==", i, "Je==", je_return[i]

		for i in xrange(ncpu):
			if( je_return[i] < min_je ):
				if Je > 0:
					min_je = je_return[i]
					n_best = i
		#print "main_node n_best===", n_best

	n_best = bcast_number_to_all(n_best,   source_node = main_node)

	if( n_best >=0):
		
		if myid == main_node:
			assign_return = [0]*(N)
			if n_best == main_node: assign_return[0:N-1] = assign[0:N-1] 
			else: assign_return = mpi_recv(N,MPI_INT, n_best, MPI_TAG_UB, MPI_COMM_WORLD)

		else:
			if n_best == myid:
				mpi_send(assign, N, MPI_INT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)

		if myid == main_node:	
			r_Cls={}
			r_Cls['n']   = [0]*K   # number of objects in a given cluster
			r_Cls['ave'] = [0]*K   # value of cluster average
			r_Cls['var'] = [0]*K   # value of cluster variance
			r_Cls['Ji']  = [0]*K   # value of ji
			r_Cls['k']   =  K	     # value of number of clusters
			r_Cls['N']   =  N
			#get 'n'	
		if myid == main_node:

			if n_best == main_node: r_Cls['n'] = Cls['n'] 
			else:
				r_Cls['n'] = mpi_recv(K,MPI_INT, n_best, MPI_TAG_UB, MPI_COMM_WORLD)
				r_Cls['n'] = map(int, r_Cls['n'])

		else:
			if n_best == myid:
				mpi_send(Cls['n'], K, MPI_INT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)
		# get 'Ji'
		if myid == main_node:

			if n_best == main_node: r_Cls['Ji'] = Cls['Ji'] 
			else: r_Cls['Ji'] = mpi_recv(K,MPI_FLOAT, n_best, MPI_TAG_UB, MPI_COMM_WORLD)

		else:
			if n_best == myid:
				mpi_send(Cls['Ji'], K, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)


		for k in xrange( K):
			tag_cls_ave = 1000
			tag_cls_var = 1010	
			if myid == main_node:

				if n_best == main_node: 
					r_Cls['ave'][k] = Cls['ave'][k] 
					r_Cls['var'][k] = Cls['var'][k]
				else: 
					r_Cls['ave'][k] = recv_EMData( n_best, tag_cls_ave )
					r_Cls['var'] [k]= recv_EMData( n_best, tag_cls_var )

			else:
				if n_best == myid:
					send_EMData( Cls['ave'][k], main_node, tag_cls_ave )	
					send_EMData( Cls['var'][k], main_node, tag_cls_var )	


			mpi_barrier(MPI_COMM_WORLD)
		mpi_barrier(MPI_COMM_WORLD)


		if myid == main_node: return assign_return, r_Cls, je_return, n_best
		else: return 0.0, 0.0, 0.0, 0.0
	else:
		if myid == main_node: return 0.0, 0.0, je_return, n_best
		else: return 0.0, 0.0, 0.0, 0.0



def k_means_SSE_collect(Cls, assign, Je, N, K, ncpu, myid, main_node):
	# Common
	from utilities   import print_begin_msg, print_end_msg, print_msg, file_type, running_time
	from statistics  import k_means_locasg2glbasg
	from time        import time
	import sys, os
	
	from mpi        import mpi_init, mpi_comm_size, mpi_comm_rank, mpi_barrier
	from mpi        import MPI_COMM_WORLD, MPI_INT, mpi_bcast
	from mpi	import MPI_FLOAT, MPI_INT, mpi_recv, mpi_send, MPI_TAG_UB
	from utilities  import bcast_number_to_all, recv_EMData, send_EMData
	
	
	
	#print "my id ==", myid, " assign [10:20] ", assign[10:20], " Je===", Je, "Cls==",Cls[ 'n' ], " Ji==", Cls['Ji']
	
			
	'''if myid == main_node:
		je_return = [0.0]*(ncpu)
		for n1 in xrange(ncpu):
			if n1 != main_node: je_return[n1]	=	mpi_recv(1,MPI_FLOAT, n1, MPI_TAG_UB, MPI_COMM_WORLD)
 			else:               je_return[main_node]  = Je
	else:
		mpi_send(Je, 1, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)'''
		
	if  myid == main_node:
		r_assign = []
		for i in xrange(ncpu):
			r_assign.append( [0]*N )

	if myid == main_node:
		for n in xrange( ncpu ):
			if n == main_node:
				r_assign[ n ] = assign
			else:
				r_assign[ n ] = mpi_recv(N,MPI_INT, n, MPI_TAG_UB, MPI_COMM_WORLD)
	else:
		mpi_send(assign, N, MPI_INT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)
		
	
	if  myid == main_node:
		r_cls = []
		for i in xrange(ncpu):
			
			t_Cls={}
			t_Cls['n']   = [0]*K   # number of objects in a given cluster
			t_Cls['ave'] = [0]*K   # value of cluster average
			t_Cls['var'] = [0]*K   # value of cluster variance
			t_Cls['Ji']  = [0]*K   # value of ji
			t_Cls['k']   =  K	     # value of number of clusters
			t_Cls['N']   =  N
			r_cls.append( t_Cls )
			del t_Cls
	#print " myiid ==", myid, "Cls['n'] before ==", Cls['n']
	mpi_barrier(MPI_COMM_WORLD)
	

	for n in xrange( ncpu ):
		if myid == main_node:
			if n == main_node:
				(r_cls[n])['n'] = Cls['n']
			else:
				(r_cls[n])['n'] = mpi_recv(K,MPI_INT, n, MPI_TAG_UB, MPI_COMM_WORLD)
		
		else:
			if myid == n:
				mpi_send(Cls['n'], K, MPI_INT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)	
		mpi_barrier(MPI_COMM_WORLD)
	'''if myid == main_node:
		for n in xrange( ncpu ):
			print " n==", n, "Cls['n'] after ==", (r_cls[ n ])['n']	'''
	if myid == main_node:
		for n in xrange( ncpu ):
			if n == main_node:
				(r_cls[n])['Ji'] = Cls['Ji']
			else:
				(r_cls[n])['Ji'] = mpi_recv(K,MPI_FLOAT, n, MPI_TAG_UB, MPI_COMM_WORLD)
	else:
		mpi_send(Cls['Ji'], K, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)
		

	for k in xrange( K):	
		tag_cls_ave = 1000
		tag_cls_var = 1010	
		if myid == main_node:
			for n in xrange( ncpu):
				if n == main_node:
					(r_cls[n])['ave'][k] = Cls['ave'][k]
					(r_cls[n])['var'][k] = Cls['var'][k]
				else:
					(r_cls[n])['ave'][k] = recv_EMData( n, tag_cls_ave )
					(r_cls[n])['var'][k] = recv_EMData( n, tag_cls_var )

		else:
			send_EMData( Cls['ave'][k], main_node, tag_cls_ave )	
			send_EMData( Cls['var'][k], main_node, tag_cls_var )	


		mpi_barrier(MPI_COMM_WORLD)
	mpi_barrier(MPI_COMM_WORLD)

	

	if myid == main_node: return r_assign, r_cls
	else: return 0.0, 0.0


	

def k_means_SSE_MPI(im_M, mask, K, rand_seed, maxit, trials, CTF, F=0, T0=0, DEBUG=False, rnd_method = 'rnd', myid = 0, main_node =0, jumping = 1):
	from utilities    import model_blank, get_im, running_time
	#from utilities    import print_begin_msg, print_end_msg, print_msg
	from random       import seed, randint, shuffle
	from copy         import deepcopy
	import sys
	import time
	#using jumping to change method for initialization
	if CTF[0]:
		from filter		import filt_ctf, filt_table
		from fundamentals	import fftip

		ctf  = deepcopy(CTF[1])
		ctf2 = deepcopy(CTF[2])
		CTF  = True
	else:
		CTF  = False

	# Simulated annealing use or not
	if T0 != 0: SA = True
	else:       SA = False

	if SA:
		# for simulated annealing
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
	if jumping ==1:
		from random import jumpahead
		if(myid != main_node):  jumpahead(17*myid+123)
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
	flag_empty = False	
	ntrials    = 0
	wd_trials  = 0
	SA_run     = SA
	
	ALL_EMPTY = True	
	while ntrials < trials:
		ntrials += 1

		# for simulated annealing
		SA = SA_run
		if SA: T = T0
	
		# Init the cluster by an image empty
		buf.to_zero()
		for k in xrange(K):
			Cls['ave'][k] = buf.copy()
			Cls['var'][k] = buf.copy()
			Cls['n'][k]   = 0
			Cls['Ji'][k]  = 0

		if rnd_method == 'd2w': assign, Cls['n'] = k_means_init_asg_d2w(im_M, N, K)
		else:	                assign, Cls['n'] = k_means_init_asg_rnd(N, K)


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

		#if DEBUG: print 'init Je', Je

		#print_msg('\n__ Trials: %2d _________________________________%s\n'%(ntrials, time.strftime('%a_%d_%b_%Y_%H_%M_%S', time.localtime())))
		#print_msg('Criterion: %11.6e \n' % Je)

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
					res = Util.min_dist_four(im_M[im], CTFxAve)
				else:
					# compute the minimum distance with centroids
					res = Util.min_dist_real(im_M[im], Cls['ave'])

				dJe = [0.0] * K
				ni  = float(Cls['n'][assign[im]])
				di  = res['dist'][assign[im]]
				for k in xrange(K):
					if k != assign[im]:
						nj  = float(Cls['n'][k])
						dj  = res['dist'][k]
						dJe[k] =  (ni/(ni-1))*(di/norm) - (nj/(nj+1))*(dj/norm)
					else:
						dJe[k] = 0	
				# Simulate Annealing
				if SA:
					
					
					# normalize and select
					mindJe = min(dJe)
					scale  = max(dJe) - mindJe
					for k in xrange(K): dJe[k] = 1 - (dJe[k] - mindJe) / scale
					select = select_kmeans(dJe, T)
					
					if select != res['pos']:
						ct_pert    += 1
						res['pos']  = select
				else:
					max_value = -1.e30
					for i in xrange( len(dJe) ):
						if( dJe[i] >= max_value) :
							max_value = dJe[i]
							res['pos'] = i

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
					if Cls['n'][assign_from] <= 1:
						#print_msg('>>> WARNING: Empty cluster, restart with new partition %d.\n\n' % wd_trials)
						flag_empty = True
												
					change = True

				# empty cluster
				if flag_empty: break

			# empty cluster
			if flag_empty: break
			
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
			if Je != 0: thd = abs(Je - old_Je) / Je
			else:       thd = 0

			# Simulated annealing, update temperature
			if SA:
				if thd < 1e-12 and ct_pert == 0: watch_dog = maxit
				T *= F
				if T < 0.009: SA = False
				#print_msg('> iteration: %5d    criterion: %11.6e    T: %13.8f  ct disturb: %5d\n' % (ite, Je, T, ct_pert))
				#if DEBUG: print '> iteration: %5d    criterion: %11.6e    T: %13.8f  ct disturb: %5d' % (ite, Je, T, ct_pert)
			else:
				if thd < 1e-8: watch_dog = maxit
				#print_msg('> iteration: %5d    criterion: %11.6e\n'%(ite, Je))
				#if DEBUG: print '> iteration: %5d    criterion: %11.6e'%(ite, Je)

			old_Je = Je

		# if no empty cluster
		if not flag_empty:

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
				#print_msg('# Criterion: %11.6e \n' % Je)
				ALL_EMPTY = False
			# set to zero watch dog trials
			wd_trials = 0
		else:
			flag_empty  = False
			wd_trials  += 1
			if wd_trials > 10:
				
				if trials > 1:
					MemJe[ntrials-1] = 1e10
					#if ntrials == trials:
						#print_msg('>>> WARNING: After ran 10 times with different partitions, one cluster is still empty. \n\n')
					#else:	print_msg('>>> WARNING: After ran 10 times with different partitions, one cluster is still empty, start the next trial.\n\n')
				else:
					#print_msg('>>> WARNING: After ran 10 times with different partitions, one cluster is still empty, STOP k-means.\n\n')
					#sys.exit()
					return 0.0, 0.0, -1e30
				wd_trials = 0
			else:
				ntrials -= 1

	'''if trials > 1:
		if ALL_EMPTY:
			#print_msg('>>> WARNING: All trials resulted in empty clusters, STOP k-means.\n\n')
			sys.exit()

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
		assign = MemAssign[best]'''
		
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
			Cls['ave'][k].do_ift_inplace()
			Cls['var'][k].do_ift_inplace()
			Cls['ave'][k].depad()
			Cls['var'][k].depad()

	# information display
	#running_time(t_start)
	#print_msg('Criterion = %11.6e \n' % Je)
	#for k in xrange(K):	print_msg('Cls[%i]: %i\n'%(k, Cls['n'][k]))
	
	# to debug
	if DEBUG: print Cls['n']
		
	# return Cls, assign and Je
	return Cls, assign, Je
	
	



# K-means MPI with classical method
def k_means_cla_MPI(IM, mask, K, rand_seed, maxit, trials, CTF, F, T0, myid, main_node, N_start, N_stop, N):
	from utilities    import model_blank, get_im
	from utilities    import bcast_EMData_to_all, reduce_EMData_to_root
	from utilities    import print_msg, running_time
	from random       import seed, randint, jumpahead
	from copy	  import deepcopy
	from mpi 	  import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_recv, mpi_send
	from mpi 	  import MPI_SUM, MPI_FLOAT, MPI_INT, MPI_LOR
	import time
	import sys
	if CTF[0]:
		from filter		import filt_ctf, filt_table
		from fundamentals	import fftip

		tmpctf  = deepcopy(CTF[1])
		tmpctf2 = deepcopy(CTF[2])
		CTF     = True
		ctf     = [None] * N
		ctf[N_start:N_stop]  = tmpctf
		ctf2    = [None] * N
		ctf2[N_start:N_stop] = tmpctf2
	else:
		CTF  = False

	## TODO change data struct to works directly with IM
	im_M = [None] * N
	im_M[N_start:N_stop] = IM

	# Simulated annealing
	if F != 0: SA = True
	else:      SA = False

	if SA:
		from math   import exp
		from random import random

	if mask != None:
		if isinstance(mask, basestring):
			ERROR('Mask must be an image, not a file name!', 'k-means', 1)

	
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
	if(myid != main_node):  jumpahead(17*myid+123)
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
		len_ctm	    = len(ctf2[N_start])

	# TRIALS
	if trials > 1: MemCls, MemJe, MemAssign = {}, {}, {}
	else: trials = 1

	ntrials    = 0
	wd_trials  = 0
	SA_run     = SA
	ALL_EMPTY = True
	while ntrials < trials:
		ntrials  += 1

		# Simulated annealing
		SA = SA_run
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
		FLAG_EXIT = 0
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
					if Cls['n'][k] <= 1:
						flag = 0
						if retrial == 0:
							print_msg('Empty class in the initialization k_means_cla_MPI\n')
							FLAG_EXIT = 1
							flag      = 1
						for k in xrange(K):
							Cls['n'][k] = 0
					k += 1
				if flag == 1: retrial = 0

		# [sync] waiting the assign is finished
		mpi_barrier(MPI_COMM_WORLD)

		# [all] check if need to exit due to initialization
		FLAG_EXIT = mpi_reduce(FLAG_EXIT, 1, MPI_INT, MPI_LOR, main_node, MPI_COMM_WORLD)
		FLAG_EXIT = mpi_bcast(FLAG_EXIT, 1, MPI_INT, main_node, MPI_COMM_WORLD)
		FLAG_EXIT = map(int, FLAG_EXIT)[0]
		if FLAG_EXIT: sys.exit()

		# [all] send assign to the others proc and the number of objects in each clusters
		assign = mpi_bcast(assign, N, MPI_INT, main_node, MPI_COMM_WORLD)
		assign = map(int, assign)     # convert array gave by MPI to list
		Cls['n'] = mpi_bcast(Cls['n'], K, MPI_FLOAT, main_node, MPI_COMM_WORLD)
		Cls['n'] = map(int, Cls['n']) # convert array gave by MPI to list
		
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
				Cls_ctf2[k] = map(float, Cls_ctf2[k])    # convert array gave by MPI to list

				reduce_EMData_to_root(Cls['ave'][k], myid, main_node)
				bcast_EMData_to_all(Cls['ave'][k], myid, main_node)

				for i in xrange(len_ctm):	Cls_ctf2[k][i] = 1.0 / Cls_ctf2[k][i]
				Cls['ave'][k] = filt_table(Cls['ave'][k], Cls_ctf2[k])

			# [id] compute Ji
			for im in xrange(N_start, N_stop):
				CTFxAve = filt_table(Cls['ave'][int(assign[im])], ctf[im])
				Cls['Ji'][int(assign[im])] += CTFxAve.cmp("SqEuclidean", im_M[im]) / norm

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
		Je = map(float, Je)[0]
		
		## Clustering		
		ite       = 0
		watch_dog = 0
		old_Je    = 0
		change    = 1
		if myid == main_node:
			print_msg('\n__ Trials: %2d _________________________________%s\n'%(ntrials, time.strftime('%a_%d_%b_%Y_%H_%M_%S', time.localtime())))
			print_msg('Criterion: %11.6e \n' % Je)
		
		while change and watch_dog < maxit:
			ite       += 1
			watch_dog += 1
			change     = 0
			Je         = 0
			if SA:
			   ct_pert = 0

			# [id] assign each images with err_min between all clusters averages
			for im in xrange(N_start, N_stop):

				# [all] compute min dist between object and centroids
				if CTF:
					CTFxAve = []
					for k in xrange(K):
						tmp = filt_table(Cls['ave'][k], ctf[im])
						CTFxAve.append(tmp.copy())
					res = Util.min_dist_four(im_M[im], CTFxAve)
				else:
					res = Util.min_dist_real(im_M[im], Cls['ave'])

				# [all] Simulated annealing
				if SA:
					dJe = [0.0] * K
					ni  = float(Cls['n'][assign[im]])
					di  = res['dist'][assign[im]]
					for k in xrange(K):
						if k != assign[im]:
							nj  = float(Cls['n'][k])
							dj  = res['dist'][k]
							dJe[k] = (ni/(ni-1))*(di/norm) - (nj/(nj+1))*(dj/norm)
						else:
							dJe[k] = 0

					# normalize and select
					mindJe = min(dJe)
					scale  = max(dJe) - mindJe
					for k in xrange(K): dJe[k] = 1 - (dJe[k] - mindJe) / scale
					select = select_kmeans(dJe, T)

					if select != res['pos']:
						ct_pert    += 1
						res['pos']  = select
				
				# [all] move object
				if res['pos'] != assign[im]:
					assign[im] = res['pos']
					change     = 1

			# [id] compute the number of objects
			for k in xrange(K): 		  Cls['n'][k] = 0
			for n in xrange(N_start, N_stop): Cls['n'][int(assign[n])] += 1			
				
			# [sync] waiting the result
			mpi_barrier(MPI_COMM_WORLD)

			# [all] sum the number of objects in each node and broadcast
			Cls['n'] = mpi_reduce(Cls['n'], K, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
			Cls['n'] = mpi_bcast(Cls['n'], K, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			Cls['n'] = map(int, Cls['n']) # convert array gave by MPI to list
			
			# [all] init average and ctf2
			FLAG_EMPTY = 0
			for k in xrange(K):
				if Cls['n'][k] <= 1:
					if myid == main_node: print_msg('>>> WARNING: Empty cluster, restart with new partition.\n\n')
					FLAG_EMPTY = 1
					break
				
				Cls['ave'][k].to_zero()
				Cls['Ji'][k] = 0
				if CTF:	Cls_ctf2[k] = [0] * len_ctm

			# [all] broadcast empty cluster information
			mpi_barrier(MPI_COMM_WORLD)
			FLAG_EMPTY = mpi_reduce(FLAG_EMPTY, 1, MPI_INT, MPI_LOR, main_node, MPI_COMM_WORLD)
			FLAG_EMPTY = mpi_bcast(FLAG_EMPTY, 1, MPI_INT, main_node, MPI_COMM_WORLD)
			FLAG_EMPTY = map(int, FLAG_EMPTY)[0]
			if FLAG_EMPTY: break
			
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
					Cls_ctf2[k] = map(float, Cls_ctf2[k]) # convert array gave by MPI to list

					reduce_EMData_to_root(Cls['ave'][k], myid, main_node)
					bcast_EMData_to_all(Cls['ave'][k], myid, main_node)
					
					for i in xrange(len_ctm):	Cls_ctf2[k][i] = 1.0 / float(Cls_ctf2[k][i])
					Cls['ave'][k] = filt_table(Cls['ave'][k], Cls_ctf2[k])

				# [id] compute Ji
				for im in xrange(N_start, N_stop):
					CTFxAve = filt_table(Cls['ave'][int(assign[im])], ctf[im])
					Cls['Ji'][int(assign[im])] += CTFxAve.cmp("SqEuclidean", im_M[im]) / norm
			
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
			Je = map(float, Je)[0]

			# threshold convergence control
			if Je != 0: thd = abs(Je - old_Je) / Je
			else:       thd = 0

			# Simulated annealing, update temperature
			if SA:
				if thd < 1e-12 and ct_pert == 0: change = 0
				T *= F
				if T < 0.009 and ct_pert < 5: SA = False
				#[id] informations display
				if myid == main_node: print_msg('> iteration: %5d    criterion: %11.6e   T: %13.8f  disturb:  %5d\n' % (ite, Je, T, ct_pert))
			else:
				if thd < 1e-8:	change = 0
				# [id] informations display
				if myid == main_node: print_msg('> iteration: %5d    criterion: %11.6e\n' % (ite, Je))
				
			old_Je = Je
			
			# [all] Need to broadcast this value because all node must run together
			change = mpi_reduce(change, 1, MPI_INT, MPI_LOR, main_node, MPI_COMM_WORLD)
			change = mpi_bcast(change, 1, MPI_INT, main_node, MPI_COMM_WORLD)
			change = map(int, change)[0]
		
		# [all] waiting the result
		mpi_barrier(MPI_COMM_WORLD)

		if not FLAG_EMPTY:
			# [id] memorize the result for this trial	
			if trials > 1:
				MemCls[ntrials-1]    = deepcopy(Cls)
				MemJe[ntrials-1]     = deepcopy(Je)
				MemAssign[ntrials-1] = deepcopy(assign)
				if myid == main_node: print_msg('# Criterion: %11.6e \n' % Je)
				ALL_EMPTY = False
			# set to zero watch dog trials
			wd_trials = 0
		else:
			FLAG_EMPTY  = 0
			wd_trials  += 1
			if wd_trials > 10:
				if trials > 1:
					MemJe[ntrials-1] = 1e10
					if ntrials == trials:
						#if myid == main_node: print_msg('>>> WARNING: After ran 10 times with different partitions, one cluster is still empty, STOP k-means.\n\n')
						#sys.exit()
						if myid == main_node: print_msg('>>> WARNING: After ran 10 times with different partitions, one cluster is still empty.\n\n')
					else:
						if myid == main_node: print_msg('>>> WARNING: After ran 10 times with different partitions, one cluster is still empty, start the next trial.\n\n')
				else:
					
					if myid == main_node: print_msg('>>> WARNING: After ran 10 times with different partitions, one cluster is still empty, STOP k-means.\n\n')
					sys.exit()
				wd_trials = 0
			else:
				ntrials -= 1
	
	if trials > 1:
		if ALL_EMPTY:
			if myid == main_node: print_msg('>>> WARNING: After ran 10 times with different partitions, one cluster is still empty, STOP k-means.\n\n')
			sys.exit()
			
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
		Cls['Ji'] = map(float, Cls['Ji'])
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
		Cls['Ji'] = map(float, Cls['Ji'])
		
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
	assign = map(int, assign) # convert array given by MPI to list

	# [main_node] write the result
	if myid == main_node and CTF:
		# ifft
		for k in xrange(K):
			Cls['ave'][k].do_ift_inplace()
			Cls['var'][k].do_ift_inplace()
			Cls['ave'][k].depad()
			Cls['var'][k].depad()
			
	# [main_node] information display
	if myid == main_node:
		running_time(t_start)
		print_msg('Criterion = %11.6e \n' % Je)
		for k in xrange(K):	print_msg('Cls[%i]: %i\n'%(k, Cls['n'][k]))

	# [all] waiting all nodes
	mpi_barrier(MPI_COMM_WORLD)
		
	if myid == main_node: return Cls, assign
	else:                 return None, None

# K-means CUDA
def k_means_CUDA(stack, mask, LUT, m, N, Ntot, K, maxit, F, T0, rand_seed, outdir, TXT, nbpart, logging = -1, flagnorm = False):
	from statistics import k_means_cuda_error, k_means_cuda_open_im
	from statistics import k_means_locasg2glbasg, k_means_cuda_export
	from utilities  import print_msg, running_time
	from time       import time
	
	# Init memory
	Kmeans = MPICUDA_kmeans()
	status = Kmeans.setup(m, N, N, K, 0)
	if status:
		k_means_cuda_error(status)
		sys.exit()
	k_means_cuda_open_im(Kmeans, stack, LUT, mask, flagnorm)
	Kmeans.compute_im2()
	status = Kmeans.init_mem(-1) # select free device
	if status:
		k_means_cuda_error(status)
		sys.exit()

	if isinstance(rand_seed, list): rnd = rand_seed
	else:                           rnd = [rand_seed]
	
	for ipart in xrange(nbpart):
		if logging != -1: logging.info('...... Start partition: %d' % (ipart + 1))
		
		Kmeans.random_ASG(rnd[ipart])
		Kmeans.compute_AVE()

		# K-means iteration
		t_start = time()
		if F != 0:
			switch_SA = True
			Kmeans.set_T(T0)
		else:   switch_SA = False
		ferror = 0
		ite    = 0
		memct  = 0
		stopct = 5
		while ite < maxit:
			if switch_SA:
				status = Kmeans.one_iter_SA()
				T      = Kmeans.get_T()
				ct     = Kmeans.get_ct_im_mv()

				print_msg('> iteration: %5d    T: %13.8f    ct disturb: %5d\n' % (ite, T, ct))
				if ct == 0: memct += 1
				else:       memct  = 0
				T *= F
				if T < 0.00001 or memct >= stopct : switch_SA = False
				Kmeans.set_T(T)
			else:   
				status = Kmeans.one_iter()
				ct     = Kmeans.get_ct_im_mv()
				print_msg('> iteration: %5d                        ct disturb: %5d\n' % (ite, ct))
				if status == 255: break

			ite += 1
			if status != 0 and status != 255:
				ferror = 1
				break

			# update
			Kmeans.compute_AVE()

		if ferror:
			k_means_cuda_error(status)
			exit()

		running_time(t_start)
		print_msg('Number of iterations        : %i\n' % ite)	
		Ji   = Kmeans.compute_ji()
		crit = Kmeans.compute_criterion(Ji)
		AVE  = Kmeans.get_AVE()
		ASG  = Kmeans.get_ASG()
		GASG = k_means_locasg2glbasg(ASG, LUT, Ntot)
		if nbpart > 1: k_means_cuda_export(GASG, AVE, outdir, mask, crit, ipart, TXT)
		else:          k_means_cuda_export(GASG, AVE, outdir, mask, crit,    -1, TXT)

	Kmeans.shutdown()
	del Kmeans

	return crit
	
# K-mean CUDA
def k_means_SSE_CUDA(stack, mask, LUT, m, N, Ntot, K, maxit, F, T0, rand_seed, outdir, TXT, nbpart, logging = -1, flagnorm = False):
	from statistics import k_means_cuda_error, k_means_cuda_open_im
	from statistics import k_means_locasg2glbasg, k_means_cuda_export
	from utilities  import print_msg, running_time
	from time       import time
	
	# Init memory
	Kmeans = MPICUDA_kmeans()
	status = Kmeans.setup(m, N, N, K, 0)
	if status:
		k_means_cuda_error(status)
		sys.exit()
	k_means_cuda_open_im(Kmeans, stack, LUT, mask, flagnorm)
	Kmeans.compute_im2() #h_im2
	status = Kmeans.init_mem(-1) # select free device
	if status:
		k_means_cuda_error(status)
		sys.exit()

	if isinstance(rand_seed, list): rnd = rand_seed
	else:                           rnd = [rand_seed]
	
	for ipart in xrange(nbpart):
		print "nbpart ==", nbpart
		if logging != -1: logging.info('...... Start partition: %d' % (ipart + 1))
		
		Kmeans.random_ASG(rnd[ipart])
		Kmeans.compute_AVE()  #get h_AVE, h_AVE2

		# K-means iteration
		t_start = time()
		if F != 0:
			switch_SA = True
			Kmeans.set_T(T0)
		else:   switch_SA = False
		ferror = 0
		ite    = 0
		memct  = 0
		stopct = 5
		while ite < maxit:
			
			if switch_SA:
				status = Kmeans.one_iter_SA()
				T      = Kmeans.get_T()
				ct     = Kmeans.get_ct_im_mv()

				print_msg('> iteration: %5d    T: %13.8f    ct disturb: %5d\n' % (ite, T, ct))
				if ct == 0: memct += 1
				else:       memct  = 0
				T *= F
				if T < 0.00001 or memct >= stopct : switch_SA = False
				Kmeans.set_T(T)
			else:   
				status = Kmeans.one_iter_SSE()
				ct     = Kmeans.get_ct_im_mv()
				print_msg('> iteration: %5d                        ct disturb: %5d\n' % (ite, ct))
				
				if status == 255: break

			ite += 1
			if status != 0 and status != 255:
				ferror = 1
				break

			# update
			Kmeans.compute_AVE()

		if ferror:
			k_means_cuda_error(status)
			exit()

		#Kmeans.AVE_to_host()
		running_time(t_start)
		print_msg('Number of iterations        : %i\n' % ite)	
		Ji   = Kmeans.compute_ji()
		crit = Kmeans.compute_criterion(Ji)
		AVE  = Kmeans.get_AVE()
		ASG  = Kmeans.get_ASG()
		GASG = k_means_locasg2glbasg(ASG, LUT, Ntot)
		if nbpart > 1: k_means_cuda_export(GASG, AVE, outdir, mask, crit, ipart, TXT)
		else:          k_means_cuda_export(GASG, AVE, outdir, mask, crit,    -1, TXT)

	Kmeans.shutdown()
	del Kmeans

	return crit


## tmp
def dump_AVE(AVE, mask, myid, ite = 0):
    #mask = get_im(maskname, 0)
    K = 256
    for k in xrange(K):
        NEWAVE = Util.reconstitute_image_mask(AVE[k], mask)
        NEWAVE.write_image('ave_from_%02i_ite_%02i.hdf' % (myid, ite), k)

# K-mean CUDA
def k_means_CUDA_MPI(stack, mask, LUT, m, N, Ntot, K, maxit, F, T0, rand_seed, myid, main_node, ncpu, outdir, TXT, nbpart, logging = -1, flagnorm = False):
	from applications import MPI_start_end
	from statistics   import k_means_cuda_error, k_means_cuda_open_im
	from statistics   import k_means_locasg2glbasg, k_means_cuda_export
	from mpi          import mpi_bcast, mpi_reduce, mpi_barrier, mpi_gatherv
	from mpi          import MPI_COMM_WORLD, MPI_INT, MPI_SUM, MPI_LOR, MPI_FLOAT
	from utilities    import print_msg, running_time
	from time         import time, sleep
	import sys

	# CST
	NGPU_PER_NODES = 4

	#if myid == main_node: t1 = time()
	# Init memory
	Kmeans                = MPICUDA_kmeans()
	N_start, N_stop       = MPI_start_end(N, ncpu, myid)
	lut                   = LUT[N_start:N_stop]
	n                     = len(lut)

	#  this is needed for gathering of ASG
	disps     = []
	recvcount = []
	for im in xrange(ncpu):
		if im == main_node:  disps.append(0)
		else:                disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(N, ncpu, im)
		recvcount.append(ie - ib)

	status = Kmeans.setup(m, N, n, K, N_start) 
	if status:
		k_means_cuda_error(status)
		sys.exit()
	k_means_cuda_open_im(Kmeans, stack, LUT, mask, flagnorm)
	Kmeans.compute_im2()
	status = Kmeans.init_mem(myid % NGPU_PER_NODES)
	if status:
		k_means_cuda_error(status)
		sys.exit()

	if myid == main_node:
		if isinstance(rand_seed, list): rnd = rand_seed
		else:                           rnd = [rand_seed]

	for ipart in xrange(nbpart):
		if logging != -1 and myid == main_node: logging.info('...... Start partition: %d' % (ipart + 1))

		# Init averages
		if myid == main_node:
			Kmeans.random_ASG(rnd[ipart])
			ASG = Kmeans.get_ASG()
		else:   ASG = None
		mpi_barrier(MPI_COMM_WORLD)
		ASG = mpi_bcast(ASG, N, MPI_INT, main_node, MPI_COMM_WORLD)
		ASG = map(int, ASG)
		Kmeans.set_ASG(ASG)
		Kmeans.compute_NC()
		Kmeans.compute_AVE()

		#if myid == main_node: print 'Init: ', time() - t1, 's'
		# K-means iterations
		if myid == main_node: tstart = time()
		if F  != 0:
			switch_SA = True
			Kmeans.set_T(T0)
		else:   switch_SA = False

		ite    = 0
		fsync  = 0
		ferror = 0
		ctconv = 0
		while ite < maxit:
			
			stop = 0
			
			#if myid == main_node: ts1 = time()
			if switch_SA:
				status = Kmeans.one_iter_SA()
				T      = Kmeans.get_T()
				ct     = Kmeans.get_ct_im_mv()
				if ct == 0: ctconv += 1
				else:       ctconv  = 0
				if myid == main_node:
					print_msg('> iteration: %5d    T: %13.8f    ct disturb: %5d %5d\n' % (ite, T, ct, ctconv))
				T *= F
				Kmeans.set_T(T)
				
				if T < 0.00001: switch_SA = False
				if ctconv >= 10: stop = 1
			else:
				status = Kmeans.one_iter()
				
				ct     = Kmeans.get_ct_im_mv()
				if myid == main_node:
					print_msg('> iteration: %5d                        ct disturb: %5d\n' % (ite, ct))
				if status != 0: stop = 1
			stop = mpi_reduce(stop, 1, MPI_INT, MPI_LOR, main_node, MPI_COMM_WORLD)
			stop = mpi_bcast(stop, 1, MPI_INT, main_node, MPI_COMM_WORLD)
			stop = int(stop[0])
			
			if stop: break

			#if myid == main_node: print 'ite cuda', time() - ts1, 's'
			ite += 1
			#if myid == main_node: ts2 = time()		
			# update
			asg = Kmeans.get_asg()
			ASG = mpi_gatherv(asg, n, MPI_INT, recvcount, disps, MPI_INT, main_node, MPI_COMM_WORLD)
			ASG = mpi_bcast(ASG, N, MPI_INT, main_node, MPI_COMM_WORLD)
			ASG = map(int, ASG)
			Kmeans.set_ASG(ASG)
			#if myid == main_node: print 'com asg', time() - ts2, 's'
			
			#if myid == main_node: ts3 = time()
			Kmeans.compute_NC()
			Kmeans.compute_AVE()
			#if myid == main_node: print 'new ave', time() - ts3, 's'
			
		#if myid == main_node:
		#	print 'Iteration time:', time() - tstart, 's'

		if status != 255 and status != 0: error = 1
		else:             error = 0

		not_empty_class_error=0
		if status != 5:
			not_empty_class_error=1
		
		not_empty_class_error = mpi_reduce(not_empty_class_error, 1, MPI_INT, MPI_LOR, main_node, MPI_COMM_WORLD)
		if myid == main_node:
			not_empty_class_error = int(not_empty_class_error[0])
			if logging != -1 and not_empty_class_error == 0:
				logging.info("EMPTY_CLASS_ERROR_K=%d"%K)	

		error = mpi_reduce(error, 1, MPI_INT, MPI_LOR, main_node, MPI_COMM_WORLD)
		error = mpi_bcast(error, 1, MPI_INT, main_node, MPI_COMM_WORLD)
		error = int(error[0])
		if error:
			k_means_cuda_error(status)
			exit()

		if myid == main_node:
			running_time(tstart)
			print_msg('Number of iterations        : %i\n' % ite)
		ji   = Kmeans.compute_ji()
		Ji   = mpi_reduce(ji, K, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
		Ji   = mpi_bcast(Ji, K, MPI_FLOAT, main_node, MPI_COMM_WORLD)
		Ji   = map(float, Ji)
		crit = Kmeans.compute_criterion(Ji)
		AVE  = Kmeans.get_AVE()
		ASG  = Kmeans.get_ASG()
		if myid == main_node:
			GASG = k_means_locasg2glbasg(ASG, LUT, Ntot)
			if nbpart > 1: k_means_cuda_export(GASG, AVE, outdir, mask, crit, ipart, TXT)
			else:          k_means_cuda_export(GASG, AVE, outdir, mask, crit,    -1, TXT)

	Kmeans.shutdown()
	del Kmeans
	
	return crit


def k_means_CUDA_MPI_YANG(stack, mask, LUT, m, N, Ntot, K, maxit, F, T0, rand_seed, myid, main_node, ncpu, outdir, TXT, ipart, logging = -1, flagnorm = False, comm = -1, gpuid = 0):
	from applications import MPI_start_end
	from statistics   import k_means_cuda_error, k_means_cuda_open_im
	from statistics   import k_means_locasg2glbasg, k_means_cuda_export
	from mpi          import mpi_bcast, mpi_reduce, mpi_barrier, mpi_gatherv
	from mpi          import MPI_COMM_WORLD, MPI_INT, MPI_SUM, MPI_LOR, MPI_FLOAT
	from utilities    import print_msg, running_time
	from time         import time, sleep
	import sys

	if comm == -1:  comm = MPI_COMM_WORLD	

	# Init memory
	Kmeans                = MPICUDA_kmeans()
	N_start, N_stop       = MPI_start_end(N, ncpu, myid)
	lut                   = LUT[N_start:N_stop]
	n                     = len(lut)

	#  this is needed for gathering of ASG
	disps     = []
	recvcount = []
	for im in xrange(ncpu):
		if im == main_node:  disps.append(0)
		else:                disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(N, ncpu, im)
		recvcount.append(ie - ib)

	status = Kmeans.setup(m, N, n, K, N_start) 
	if status:
		k_means_cuda_error(status)
		sys.exit()
	k_means_cuda_open_im(Kmeans, stack, LUT, mask, flagnorm)
	Kmeans.compute_im2()
	status = Kmeans.init_mem(gpuid)
	if status:
		k_means_cuda_error(status)
		sys.exit()

	if myid == main_node:
		if isinstance(rand_seed, list): rnd = rand_seed
		else:                           rnd = [rand_seed]

	if logging != -1 and myid == main_node: logging.info('...... Start partition: %d' % (ipart + 1))

	# Init averages
	if myid == main_node:
		Kmeans.random_ASG(rnd[ipart])
		ASG = Kmeans.get_ASG()
	else:   ASG = None
	mpi_barrier(comm)
	ASG = mpi_bcast(ASG, N, MPI_INT, main_node, comm)
	ASG = map(int, ASG)
	Kmeans.set_ASG(ASG)
	Kmeans.compute_NC()
	Kmeans.compute_AVE()

	#if myid == main_node: print 'Init: ', time() - t1, 's'
	# K-means iterations
	if myid == main_node: tstart = time()
	if F  != 0:
		switch_SA = True
		Kmeans.set_T(T0)
	else:   switch_SA = False

	ite    = 0
	fsync  = 0
	ferror = 0
	ctconv = 0
	while ite < maxit:
		stop = 0
		
		#if myid == main_node: ts1 = time()
		if switch_SA:
			status = Kmeans.one_iter_SA()
			T      = Kmeans.get_T()
			ct     = Kmeans.get_ct_im_mv()
			if ct == 0: ctconv += 1
			else:       ctconv  = 0
			if myid == main_node:
				print_msg('> iteration: %5d    T: %13.8f    ct disturb: %5d %5d\n' % (ite, T, ct, ctconv))
			T *= F
			Kmeans.set_T(T)
			
			if T < 0.00001: switch_SA = False
			if ctconv >= 10: stop = 1
		else:
			status = Kmeans.one_iter()
			ct     = Kmeans.get_ct_im_mv()
			if myid == main_node:
				print_msg('> iteration: %5d                        ct disturb: %5d\n' % (ite, ct))
			if status != 0: stop = 1
		stop = mpi_reduce(stop, 1, MPI_INT, MPI_LOR, main_node, comm)
		stop = mpi_bcast(stop, 1, MPI_INT, main_node, comm)
		stop = int(stop[0])
		if stop: break

		#if myid == main_node: print 'ite cuda', time() - ts1, 's'
		ite += 1
		#if myid == main_node: ts2 = time()		
		# update
		asg = Kmeans.get_asg()
		ASG = mpi_gatherv(asg, n, MPI_INT, recvcount, disps, MPI_INT, main_node, comm)
		ASG = mpi_bcast(ASG, N, MPI_INT, main_node, comm)
		ASG = map(int, ASG)
		Kmeans.set_ASG(ASG)
		#if myid == main_node: print 'com asg', time() - ts2, 's'
		
		#if myid == main_node: ts3 = time()
		Kmeans.compute_NC()
		Kmeans.compute_AVE()
		#if myid == main_node: print 'new ave', time() - ts3, 's'
		
	#if myid == main_node:
	#	print 'Iteration time:', time() - tstart, 's'

	if status != 255 and status != 0: error = 1
	else:             error = 0

	not_empty_class_error=0
	if status != 5:
		not_empty_class_error=1
	
	not_empty_class_error = mpi_reduce(not_empty_class_error, 1, MPI_INT, MPI_LOR, main_node, comm)
	if myid == main_node:
		not_empty_class_error = int(not_empty_class_error[0])
		if logging != -1 and not_empty_class_error == 0:
			logging.info("EMPTY_CLASS_ERROR_K=%d"%K)	

	error = mpi_reduce(error, 1, MPI_INT, MPI_LOR, main_node, MPI_COMM_WORLD)
	error = mpi_bcast(error, 1, MPI_INT, main_node, MPI_COMM_WORLD)
	error = int(error[0])
	if error:
		k_means_cuda_error(status)
		exit(5)

	if myid == main_node:
		running_time(tstart)
		print_msg('Number of iterations        : %i\n' % ite)
	ji   = Kmeans.compute_ji()
	Ji   = mpi_reduce(ji, K, MPI_FLOAT, MPI_SUM, main_node, comm)
	Ji   = mpi_bcast(Ji, K, MPI_FLOAT, main_node, comm)
	Ji   = map(float, Ji)
	crit = Kmeans.compute_criterion(Ji)
	AVE  = Kmeans.get_AVE()
	ASG  = Kmeans.get_ASG()
	if myid == main_node:
		GASG = k_means_locasg2glbasg(ASG, LUT, Ntot)
		k_means_cuda_export(GASG, AVE, outdir, mask, crit, ipart, TXT)

	Kmeans.shutdown()
	del Kmeans
	
	return crit

## K-MEANS GROUPS ######################################################################

# make script file to gnuplot
def k_means_groups_gnuplot(file, src, C, DB, H):
	out = open(file, 'w')
	out.write('# Gnuplot script file for plotting result in kmeans groups\n')
	out.write('# $ gnuplot %s.p\n' % src)
	out.write('reset\n')
	out.write('set autoscale\n')
	txt = 'plot'

	WORLD = [C, DB, H]
	name  = ['Coleman', 'Davies-Bouldin', 'Harabasz']
	pos   = [3, 5, 7]

	# norm plot [0;1]
	for i in xrange(3):
		minv = min(WORLD[i])
		maxv = max(WORLD[i])
		d = maxv - minv
		a = 1 / d
		b = 0.5 - ((maxv + minv) * a / 2)
		txt += ' \x22%s\x22 u 1:($%d*(%11.4e)+(%11.4e)) ti \x22%s\x22 w l,' % (src, pos[i], a, b, name[i])

	out.write(txt.rstrip(',') + '\n')
	out.close()

# to figure out the number of clusters
def k_means_groups_serial(stack, outdir, maskname, opt_method, K1, K2, rand_seed, maxit, trials, CTF, F, T0, DEBUG = False, flagnorm = False):
	from utilities   import print_begin_msg, print_end_msg, print_msg, running_time, file_type
	from statistics  import k_means_open_im, k_means_criterion, k_means_headlog
	from statistics  import k_means_cla, k_means_SSE, k_means_groups_gnuplot
	from statistics  import k_means_init_open_im
	import os, sys, time

	if os.path.exists(outdir): ERROR('Output directory exists, please change the name and restart the program', "k_means_groups_serial", 1)
	os.mkdir(outdir)

	t_start = time.time()
	print_begin_msg('k-means groups')	

	ext = file_type(stack)
	if ext == 'txt': TXT = True
	else:            TXT = False
	LUT, mask, N, m, Ntot = k_means_init_open_im(stack, maskname)
	IM, ctf, ctf2         = k_means_open_im(stack, mask, CTF, LUT, flagnorm)
	k_means_headlog(stack, outdir, opt_method, N, [K1, K2], 'CHD', maskname, trials, maxit,\
				CTF, T0, F, rand_seed, 1, m)
	
	# init
	KK       = range(K1, K2 + 1)	# Range of works
	C, DB, H = [], [], []
	sp       = 15                   # cst space to format file

	# init the file result
	file_crit = open(outdir + '/' + outdir, 'w')
	file_crit.write('# Criterion of k-means group\n')
	file_crit.write('# %s %s %s  %s\n' % ('N ', 'Coleman'.ljust(sp), 'Davies-Bouldin'.ljust(sp), 'Harabasz'.ljust(sp)))
	file_crit.close()
		
	# Compute the criterion and format
	for K in KK:
		
		print_msg('\n')
		print_msg('| K=%d |====================================================================\n' % K)

		try:
			if opt_method   == 'cla': [Cls, assign] = k_means_cla(IM, mask, K, rand_seed, maxit, trials, [CTF, ctf, ctf2], F, T0, DEBUG)
			elif opt_method == 'SSE': [Cls, assign] = k_means_SSE(IM, mask, K, rand_seed, maxit, trials, [CTF, ctf, ctf2], F, T0, DEBUG)
			else:			  ERROR('Kind of k-means unknown', 'k_means_groups', 1)

		except SystemExit:
			ERROR('Empty cluster, number of groups too high %d'%K, 'k_means_groups', 1)
		
		crit = k_means_criterion(Cls, 'CHD')

		# res file
		file_crit = open(outdir + '/' + outdir, 'a')
		file_crit.write('%3d  C: %11.4e  DB: %11.4e  H: %11.4e | %s\n' % (K, crit['C'], crit['D'], crit['H'], time.ctime())) 
		file_crit.close()
		
		# mem
		C.append(crit['C'])
		H.append(crit['H'])
		DB.append(crit['D'])
		
	# gnuplot script
	k_means_groups_gnuplot(outdir + '/' + outdir + '.p', outdir, C, DB, H)
	running_time(t_start)
	print_end_msg('k-means groups')

# to figure out the number of clusters CUDA version
def k_means_groups_CUDA(stack, outdir, maskname, K1, K2, rand_seed, maxit, F, T0):
	from utilities   import print_begin_msg, print_end_msg, print_msg, running_time, file_type
	from statistics  import k_means_cuda_init_open_im, k_means_cuda_headlog
	from statistics  import k_means_groups_gnuplot, k_means_CUDA
	import time, os, sys

	if os.path.exists(outdir): ERROR('Output directory exists, please change the name and restart the program', "k_means_groups_CUDA", 1)
	os.mkdir(outdir)
	t_start = time.time()
	
	ext = file_type(stack)
	if ext == 'txt': TXT = True
	else:            TXT = False

	# init to open images
	LUT, mask, N, m, Ntot = k_means_cuda_init_open_im(stack, maskname)
	# write logfile
	print_begin_msg('k-means groups')
	k_means_cuda_headlog(stack, outdir, 'cuda', N, [K1, K2], maskname, maxit, T0, F, rand_seed, 1, m)

	# init
	KK       = range(K1, K2 + 1)	# Range of works
	C, DB, H = [], [], []
	sp       = 15                   # cst space to format file result
		
	# init the file result
	file_crit = open(outdir + '/' + outdir, 'w')
	file_crit.write('# Criterion of k-means group\n')
	file_crit.write('# %s %s %s  %s\n' % ('N ', 'Coleman'.ljust(sp), 'Davies-Bouldin'.ljust(sp), 'Harabasz'.ljust(sp)))
	file_crit.close()
		
	# Compute the criterion and format
	for K in KK:
		print_msg('\n')
		print_msg('| K=%d |====================================================================\n' % K)

		#try:
		crit = k_means_CUDA(stack, mask, LUT, m, N, Ntot, K, maxit, F, T0, rand_seed, outdir, TXT, 1)
		#except SystemExit:
		#	ERROR('Empty cluster or device error', 'k_means_groups_CUDA', 1)

		Je, Ci, Hi, DBi = crit
		# result file
		file_crit = open(outdir + '/' + outdir, 'a')
		file_crit.write('%3d  C: %11.4e  DB: %11.4e  H: %11.4e | %s\n' % (K, Ci, DBi, Hi, time.ctime())) 
		file_crit.close()

		# mem for latter
		C.append(Ci)
		DB.append(DBi)
		H.append(Hi)

	k_means_groups_gnuplot(outdir + '/' + outdir + '.p', outdir, C, DB, H)
		
	# runtime
	running_time(t_start)
	print_end_msg('k-means groups')

# to figure out the number of clusters MPI version
def k_means_groups_MPI(stack, outdir, maskname, opt_method, K1, K2, rand_seed, maxit, trials, CTF, F, T0, flagnorm):
	from utilities    import print_begin_msg, print_end_msg, print_msg, running_time, file_type
	from statistics   import k_means_open_im, k_means_criterion, k_means_headlog
	from statistics   import k_means_cla_MPI, k_means_SSE_MPI
	from applications import MPI_start_end
	from mpi 	  import mpi_init, mpi_comm_size, mpi_comm_rank, mpi_barrier, MPI_COMM_WORLD, mpi_bcast, MPI_INT, mpi_send, mpi_recv, MPI_TAG_UB
	import sys, os, time
	from utilities import bcast_number_to_all

	sys.argv  = mpi_init(len(sys.argv), sys.argv)
	ncpu      = mpi_comm_size(MPI_COMM_WORLD)
	myid      = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	ext = file_type(stack)
	if ext == 'txt': TXT = True
	else:            TXT = False
	
	if os.path.exists(outdir): ERROR('Output directory exists, please change the name and restart the program', "k_means_groups_MPI", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)
	
	if myid == main_node:
		#print_begin_msg('k-means groups_MPI')
		t_start = time.time()
		os.mkdir(outdir)
	
	LUT, mask, N, m, Ntot = k_means_init_open_im(stack, maskname)
	
	IM, ctf, ctf2         = k_means_open_im(stack, mask, CTF, LUT, flagnorm)

	if myid == main_node: k_means_headlog(stack, outdir, opt_method, N, [K1, K2], 'CHD', maskname, ncpu, maxit, CTF, T0, F, rand_seed, ncpu, m)

	KK = range(K1, K2 + 1)	# Range of works
	if myid == main_node:
		# init
		C, DB, H = [], [], []
		sp       = 15                   # cst space to format file

		# init the file result
		file_crit = open(outdir + '/' + outdir, 'w')
		file_crit.write('# Criterion of k-means group\n')
		file_crit.write('# %s %s %s  %s\n' % ('N ', 'Coleman'.ljust(sp), 'Davies-Bouldin'.ljust(sp), 'Harabasz'.ljust(sp)))
		file_crit.close()

		#k_means_headlog(stack, outdir, opt_method, N, [K1, K2], 'CHD', maskname, trials, maxit, CTF, T0, F, rand_seed, ncpu, m)

	# get some criterion
	for K in KK:
		if myid == main_node:
			print_msg('\n')
			print_msg('| K=%d |====================================================================\n' % K)
			t_start1 = time.time()

		

		[Cls, assign, Je] = k_means_SSE_MPI(IM, mask, K, rand_seed, maxit, 
					1, [CTF, ctf, ctf2], F, T0, False, "rnd", myid = myid, main_node = main_node, jumping = 1)


		from statistics import k_means_SSE_combine
		[ assign_return, r_Cls, je_return, n_best] = k_means_SSE_combine(Cls, assign, Je, N, K, ncpu, myid, main_node)
		if myid == main_node: running_time(t_start1)
		n_best_get = 0
		mpi_barrier(MPI_COMM_WORLD)
		if myid == main_node:
			for n1 in xrange(ncpu):
				if n1 != main_node: mpi_send(n_best, 1, MPI_INT, n1, MPI_TAG_UB, MPI_COMM_WORLD) 
 				else:               n_best_get  = n_best
		else: n_best_get	=	mpi_recv(1,MPI_INT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)
		n_best_get = int(n_best_get)
		mpi_barrier(MPI_COMM_WORLD)
		#print "myid==",myid," n_best==", n_best_get
		mpi_barrier(MPI_COMM_WORLD)
		if myid == main_node:

			if n_best == -1:
				print_msg('>>> WARNING: All trials resulted in empty clusters, STOP k-means.\n\n')
				#print "assign_return===", assign_return[10:20], "cls_n return==", r_Cls['n'], "Ji==", r_Cls['Ji'], "ave size ==", r_Cls['ave'][0].get_xsize()
			else:
				for i in xrange( ncpu ):
					if( je_return[i] <0 ):
						print_msg('> Trials: %5d    resulted in empty clusters  \n' % (i) )
					else:
						print_msg('> Trials: %5d    criterion: %11.6e  \n' % (i, je_return[i]) )
				crit = k_means_criterion(r_Cls, 'CHD')
				#glb_assign = k_means_locasg2glbasg(assign_return, LUT, Ntot)
				#k_means_export(r_Cls, crit, glb_assign, outdir, -1, TXT)
				#print_end_msg('k-means MPI end')
			
		if n_best_get== -1:
			sys.exit()
		
		
		if myid == main_node:
			crit = k_means_criterion(r_Cls, 'CHD')

			# res file
			file_crit = open(outdir + '/' + outdir, 'a')
			file_crit.write('%3d  C: %11.4e  DB: %11.4e  H: %11.4e | %s\n' % (K, crit['C'], crit['D'], crit['H'], time.ctime())) 
			file_crit.close()

			# mem
			C.append(crit['C'])
			H.append(crit['H'])
			DB.append(crit['D'])

		mpi_barrier(MPI_COMM_WORLD)
		
	# gnuplot script
	if myid == main_node:
		k_means_groups_gnuplot(outdir + '/' + outdir + '.p', outdir, C, DB, H)
		running_time(t_start)
		print_end_msg('k-means groups')

## K-MEANS CUDA ###########################################################################
# 2009-02-20 15:39:43

# k-means print out error given by the cuda code
def k_means_cuda_error(status):
	from utilities import print_msg
	# status info:
	#   0 - all is ok
	#   1 - error to init host memory
	#   2 - error to init device memory
	#   3 - error system on device
	#   4 - init assignment with empty class
	#   5 - classification return empty class
	#   6 - error to select the device
	# 255 - k-means done
	if status == 0 or status == 255: return
	print_msg('============================================\n')
	if   status == 1: print_msg('* ERROR: allocation host memory            *\n')
	elif status == 2: print_msg('* ERROR: allocation device memory          *\n')
	elif status == 3: print_msg('* ERROR: system device                     *\n')
	elif status == 4: print_msg('* ERROR: random assignment (empty class)   *\n')
	elif status == 5: print_msg('* ERROR: classification return empty class *\n')
	elif status == 6: print_msg('* ERROR: fail to select the device         *\n')
	print_msg('============================================\n')

# k-means write the head of the logfile for CUDA
def k_means_cuda_headlog(stackname, outname, method, N, K, maskname, maxit, T0, F, rnd, ncpu, m):
	from utilities import print_msg
	from math import log

	if F != 0: SA = True
	else:      SA = False

	if method == 'cla': method = 'No optimisation'

	if ncpu > 1: methodhead = 'CUDA MPI'
	else:        methodhead = 'CUDA'

	if isinstance(K, list): 
		txtK = '%i to %i' % (K[0], K[1])
		k    = K
	else:	
		txtK = str(K)
		k    = [K]
	if isinstance(rnd, list):
		txtrnd = '%i ' * len(rnd)
		txtrnd = txtrnd % tuple(rnd)
	else:
		txtrnd = str(rnd)

	# memory estimation
	#        IM                         AVE              DIST
	device = N * m * 4 / float(ncpu) + max(k) * m * 4 + N * max(k) * 4 / float(ncpu)
	#        IM          AVE              DIST
	host   = N * m * 4 + max(k) * m * 4 + N * max(k) * 4 / float(ncpu)      
	#       ASG     NC           IM2                   AVE2
	host  += N * 2 + max(k) * 4 + N * 4 / float(ncpu) + max(k) * 4
	ie_device  = int(log(device) // log(1e3))
	ie_host    = int(log(host)   // log(1e3))
	device    /= (1e3 ** ie_device)
	host      /= (1e3 ** ie_host)
	txt        = ['', 'k', 'M', 'G', 'T']
	device     = '%5.2f %sB' % (device, txt[ie_device])
	host       = '%5.2f %sB' % (host,   txt[ie_host])

	print_msg('\n************* k-means %s *************\n' % methodhead)
	print_msg('Input stack                 : %s\n'     % stackname)
	print_msg('Number of images            : %i\n'     % N)
	print_msg('Maskfile                    : %s\n'     % maskname)
	print_msg('Number of pixels under mask : %i\n'     % m)
	print_msg('Number of clusters          : %s\n'     % txtK)
	print_msg('Maximum iteration           : %i\n'     % maxit)
	print_msg('Criterion                   : CHD\n'    )
	print_msg('Optimization method         : %s\n'     % method)
	if SA:
		print_msg('Simulated annealing          : ON\n')
		print_msg('   F                        : %f\n' % F)
		if T0 != -1: print_msg('   T0                       : %f\n' % T0)
		else:        print_msg('   T0                       : AUTO\n')

	else:
		print_msg('Simulated annealing          : OFF\n')
	print_msg('Random seed                 : %s\n'     % txtrnd)
	print_msg('Number of Cs              : %i\n'     % ncpu)
	print_msg('Output seed names           : %s\n'     % outname)
	print_msg('Memory on device            : %s\n'     % device)
	print_msg('Memory on host              : %s\n\n'   % host)

# k-means, prepare to open images later
def k_means_cuda_init_open_im(stack, maskname):
	from utilities import get_image, get_im, model_blank, file_type
	from EMAN2db import db_open_dict

	ext = file_type(stack)
	if ext == 'txt': TXT = True
	else:            TXT = False

	# open mask if defined
	if maskname != None: mask = get_image(maskname)
	else:
		# anyway image must be a flat image
		if TXT:
			line = open(stack, 'r').readline()
			nx   = len(line.split())
			mask = model_blank(nx)
			mask.to_one()
		else:
			im = get_im(stack, 0)
			mask = model_blank(im.get_xsize(), im.get_ysize(), im.get_zsize())
			mask.to_one()
			del im

	# get some params
	if TXT: Ntot = len(open(stack, 'r').readlines())
	else:   Ntot = EMUtil.get_image_count(stack)
	im = Util.compress_image_mask(mask, mask)
	m  = im.get_xsize()
	del im

	# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
	# # check if the flag active is used, in the case where k-means will run for the stability
	# if TXT:
	# 	flagactive = False
	# else:
	# 	image = EMData()
	# 	image.read_image(stack, 0, True)
	# 	flagactive = True
	# 	try:	active = image.get_attr('active')
	# 	except: flagactive = False
	# 	del image
	# 
	# 
	# # if flag active used, prepare the list of images
	# if flagactive:
	# 	lim  = []
	# 	ext  = file_type(stack)
	# 	if ext == 'bdb':
	# 		DB = db_open_dict(stack)
	# 		for n in xrange(Ntot):
	# 			if DB.get_attr(n, 'active'): lim.append(n)
	# 		DB.close()
	# 	else:
	# 		im = EMData.read_images(stack, range(Ntot), True)
	# 		for n in xrange(Ntot):
	# 			if im[n].get_attr('active'): lim.append(n)
	# 		del im
	# 	N = len(lim)
	# else:
	# 	lim = range(Ntot)
	# 	N = Ntot

	lim = range(Ntot)
	N = Ntot

	return lim, mask, N, m, Ntot

# k-means open, prepare, and load images for CUDA k-means
def k_means_cuda_open_im(KmeansCUDA, stack, lim, mask, flagnorm = False):
	from utilities     import get_params2D, get_params3D, get_im, file_type, model_blank
	from fundamentals  import rot_shift2D, rot_shift3D
	
	ext = file_type(stack)
	if ext == 'txt': TXT = True
	else:            TXT = False

	# to manage coord fact in text file format
	if TXT:
		c    = 0
		data = open(stack, 'r').readlines()
		nx   = len(data[0].split())
		for line in data:
			im   = model_blank(nx)
			line = line.split()
			for i in xrange(nx):
				val = float(line[i])
				im.set_value_at_fast(i, 0, val)
			im = Util.compress_image_mask(im, mask)
			KmeansCUDA.append_flat_image(im, c)
			c += 1
			
		return

	# some parameters
	image = get_im(stack, 0)
	nx = image.get_xsize()
	ny = image.get_ysize()
	nz = image.get_zsize()
	del image

	# open one by one to avoid twice allocation of memory (python/C)
	# even if it takes more time
	c = 0
	for i in lim:
		image = get_im(stack, i)
		
		# 3D object
		if nz > 1:
			try:
				phi, theta, psi, s3x, s3y, s3z, mirror, scale = get_params3D(image)
				image = rot_shift3D(image, phi, theta, psi, s3x, s3y, s3z, scale)
				if mirror: image.process_inplace('xfrom.mirror', {'axis':'x'})
			except:
				ERROR('K-MEANS no 3D alignment parameters found', "k_means_cuda_open_im", 1)
				sys.exit()
		# 2D object
		elif ny > 1:
			try:
				alpha, sx, sy, mirror, scale = get_params2D(image)
				image = rot_shift2D(image, alpha, sx, sy, mirror, scale)
			except: 
				ERROR('K_MEANS no 2D alignment parameters found', "k_means_cuda_open_im", 1)
				sys.exit()

		if flagnorm:
			# normalize
			ave, std, mi, mx = Util.infomask(image, mask, True)
			image -= ave
			image /= std

		# apply mask 
		image = Util.compress_image_mask(image, mask)

		# load to C function through the kmeansCUDA object
		KmeansCUDA.append_flat_image(image, c)
		c += 1

# K-means write only the major info to the header, call by the stability process
def k_means_cuda_info(INFO):
	from utilities import print_msg
	
	# write the report on the logfile
	time_run = int(INFO['time'])
	time_h   = time_run / 3600
	time_m   = (time_run % 3600) / 60
	time_s   = (time_run % 3600) % 60
	
	print_msg('Running time is             : %s h %s min %s s\n' % (str(time_h).rjust(2, '0'), str(time_m).rjust(2, '0'), str(time_s).rjust(2, '0')))
	print_msg('Number of iterations        : %i\n' % INFO['noi'])
	print_msg('Partition criterion is      : %11.6e (total sum of squares error)\n' % INFO['Je'])
	print_msg('Criteria Coleman is         : %11.6e\n' % INFO['C'])
	print_msg('Criteria Harabasz is        : %11.6e\n' % INFO['H'])
	print_msg('Criteria Davies-Bouldin is  : %11.6e\n' % INFO['DB'])

# K-means write results output directory
def k_means_cuda_export(PART, FLATAVE, out_seedname, mask, crit, part = -1, TXT = False):
	from utilities import print_msg
	import os
	if not os.path.exists(out_seedname): os.mkdir(out_seedname)

	Je, C, H, DB = crit
	print_msg('Partition criterion is      : %11.6e (total sum of squares error)\n' % Je)
	print_msg('Criteria Coleman is         : %11.6e\n' % C)
	print_msg('Criteria Harabasz is        : %11.6e\n' % H)
	print_msg('Criteria Davies-Bouldin is  : %11.6e\n' % DB)

	# prepare list of images id for each group
	K   = max(PART) + 1
	N   = len(PART)
	GRP = [[] for i in xrange(K)]
	for n in xrange(N):
		# if image are assigned somewhere (active)
		if int(PART[n]) != -1: GRP[int(PART[n])].append(n)

	flagHDF = False
	for k in xrange(K):
		if len(GRP[k]) > 16000: flagHDF = True
	if flagHDF: print_msg('\nWARNING: limitation of number attributes in hdf format, the results will be export in separate text files\n')

	# write the details of the clustering
	print_msg('\n-- Details ----------------------------\n')
	for k in xrange(K):
		print_msg('\t%s\t%d\t%s\t%d\n' % ('Cluster no:', k, 'No of Objects = ', len(GRP[k])))

		# reconstitute averages
		AVE = Util.reconstitute_image_mask(FLATAVE[k], mask)
		
		# limitation of hdf format
		if flagHDF or TXT:
			if part != -1:
				outfile = open(os.path.join(out_seedname, 'k_means_part_%02i_grp_%03i.txt' % (part, k)), 'w')
			else:
				outfile = open(os.path.join(out_seedname, 'k_means_grp_%03i.txt' % (k)), 'w')
			for id in GRP[k]: outfile.write('%i\n' % int(id))
			outfile.close()
			AVE.set_attr_dict({'Class_average':1.0, 'nobjects': len(GRP[k])})
		else:
			AVE.set_attr('Class_average', 1.0)
			AVE.set_attr('nobjects', len(GRP[k]))
			AVE.set_attr('members', GRP[k])

		if part == -1: AVE.write_image(os.path.join(out_seedname, 'averages.hdf'), k)
		else:          AVE.write_image(os.path.join(out_seedname, 'averages_%02i.hdf' % part), k)
	print_msg('\n')

## K-MEANS STABILITY ######################################################################
# 2008-12-18 11:35:11 

# K-means SA define the first temperature T0 with a couple of testing values
def k_means_SA_T0(im_M, mask, K, rand_seed, CTF, F):
	from utilities 		import model_blank, print_msg
	from alignment          import select_k
	from random    		import seed, randint
	import sys
	import time
	if CTF[0]:
		from filter	        import filt_ctf, filt_table
		from fundamentals 	import fftip

		ctf  = deepcopy(CTF[1])
		ctf2 = deepcopy(CTF[2])
		CTF  = True
	else:
		CTF  = False

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
			if Cls['n'][k] <= 1:
				flag = 0
				if retrial == 0: sys.exit()
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
	th = int(float(N)*0.8)
	T0 = -1
	lT = []
	Tm = 40
	for i in xrange(1, 10): lT.append(i/10.)
	lT.extend(range(1, 5))
	lT.extend(range(5, Tm, 2))
	for T in lT:
		ct_pert = 0
		for rep in xrange(2):
			for im in xrange(N):
				if CTF:
					CTFxAVE = []
					for k in xrange(K): CTFxAVE.append(filt_table(Cls['ave'][k], ctf[im]))
					res = Util.min_dist_four(im_M[im], CTFxAVE)
				else:
					res = Util.min_dist_real(im_M[im], Cls['ave'])
		
				# Simulated annealing
				dJe = [0.0] * K
				ni  = float(Cls['n'][assign[im]])
				di  = res['dist'][assign[im]]											
				for k in xrange(K):
					if k != assign[im]:
						nj  = float(Cls['n'][k])
						dj  = res['dist'][k]
						dJe[k] = (ni/(ni-1))*(di/norm) - (nj/(nj+1))*(dj/norm)
					else:
						dJe[k] = 0

				# normalize and select
				mindJe = min(dJe)
				scale  = max(dJe) - mindJe
				for k in xrange(K): dJe[k] = (dJe[k] - mindJe) / scale
				select = select_k(dJe, T)

				if select != res['pos']:
					ct_pert    += 1
					res['pos']  = select

		ct_pert /= 2.0

		# select the first temperature if > th
		if ct_pert > th:
			T0 = T
			break

	# if not found, set to the max value
	if T0 == -1: T0 = Tm
	
	# return Cls, assign
	return T0, ct_pert

# K-means SA define the first temperature T0 (MPI version) with a couple of testing values
def k_means_SA_T0_MPI(im_M, mask, K, rand_seed, CTF, F, myid, main_node, N_start, N_stop):
	from utilities 		import model_blank, print_msg, bcast_EMData_to_all, reduce_EMData_to_root
	from random    		import seed, randint
	from alignment          import select_k
	from mpi                import mpi_reduce, mpi_bcast, mpi_barrier, mpi_recv, mpi_send
	from mpi                import MPI_SUM, MPI_FLOAT, MPI_INT, MPI_LOR, MPI_COMM_WORLD
	from copy               import deepcopy
	import sys
	import time
	if CTF[0]:
		from filter	        import filt_ctf, filt_table
		from fundamentals 	import fftip

		ctf  = deepcopy(CTF[1])
		ctf2 = deepcopy(CTF[2])
		CTF  = True
	else:
		CTF  = False

	from math   import exp
	from random import random

	if mask != None:
		if isinstance(mask, basestring):
			ERROR('Mask must be an image, not a file name!', 'k-means', 1)

	N = len(im_M)

	t_start = time.time()
		
	# Informations about images
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
		len_ctm	    = len(ctf2[N_start])

	# Init the cluster by an image empty
	buf.to_zero()
	for k in xrange(K):
		Cls['ave'][k] = buf.copy()
		Cls['var'][k] = buf.copy()
		Cls['n'][k]   = 0
		Cls['Ji'][k]  = 0

	## [main] Random method
	FLAG_EXIT = 0
	if myid == main_node:
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
				if Cls['n'][k] <= 1:
					flag = 0
					if retrial == 0: FLAG_EXIT = 1
					for k in xrange(K):
						Cls['n'][k] = 0

			if flag == 1:	retrial = 0

	# if need all node quit
	mpi_barrier(MPI_COMM_WORLD)
	FLAG_EXIT = mpi_reduce(FLAG_EXIT, 1, MPI_INT, MPI_LOR, main_node, MPI_COMM_WORLD)
	FLAG_EXIT = mpi_bcast(FLAG_EXIT, 1, MPI_INT, main_node, MPI_COMM_WORLD)
	FLAG_EXIT = map(int, FLAG_EXIT)[0]
	mpi_barrier(MPI_COMM_WORLD)
	if FLAG_EXIT: sys.exit()

	# [sync] waiting assignment
	assign = mpi_bcast(assign, N, MPI_INT, main_node, MPI_COMM_WORLD)
	assign = map(int, assign)     # convert array gave by MPI to list
	Cls['n'] = mpi_bcast(Cls['n'], K, MPI_FLOAT, main_node, MPI_COMM_WORLD)
	Cls['n'] = map(float, Cls['n']) # convert array gave by MPI to list

	## Calculate averages, if CTF: ave = S CTF.F / S CTF**2
	if CTF:
		# first init ctf2
		for k in xrange(K):	Cls_ctf2[k] = [0] * len_ctm

		for im in xrange(N_start, N_stop):
			# compute ctf2				
			for i in xrange(len_ctm):	Cls_ctf2[assign[im]][i] += ctf2[im][i]

			# compute average first step
			CTFxF = filt_table(im_M[im], ctf[im])
			Util.add_img(Cls['ave'][assign[im]], CTFxF)

		# sync
		mpi_barrier(MPI_COMM_WORLD)
		for k in xrange(K):
			Cls_ctf2[k] = mpi_reduce(Cls_ctf2[k], len_ctm, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
			Cls_ctf2[k] = mpi_bcast(Cls_ctf2[k],  len_ctm, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			Cls_ctf2[k] = map(float, Cls_ctf2[k])    # convert array gave by MPI to list
			reduce_EMData_to_root(Cls['ave'][k], myid, main_node)
			bcast_EMData_to_all(Cls['ave'][k], myid, main_node)

			for i in xrange(len_ctm):	Cls_ctf2[k][i] = 1.0 / float(Cls_ctf2[k][i])
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
	th = int(float(N)*0.8)
	T0 = -1
	lT = []
	Tm = 40
	for i in xrange(1, 10): lT.append(i/10.)
	lT.extend(range(1, 5))
	lT.extend(range(5, Tm, 2))
	for T in lT:
		ct_pert = 0
		for rep in xrange(2):
			for im in xrange(N_start, N_stop):
				if CTF:
					CTFxAVE = []
					for k in xrange(K): CTFxAVE.append(filt_table(Cls['ave'][k], ctf[im]))
					res = Util.min_dist_four(im_M[im], CTFxAVE)
				else:
					res = Util.min_dist_real(im_M[im], Cls['ave'])

				# Simulated annealing
				dJe = [0.0] * K
				ni  = float(Cls['n'][assign[im]])
				di  = res['dist'][assign[im]]											
				for k in xrange(K):
					if k != assign[im]:
						nj  = float(Cls['n'][k])
						dj  = res['dist'][k]
						dJe[k] = (ni/(ni-1))*(di/norm) - (nj/(nj+1))*(dj/norm)
					else:
						dJe[k] = 0

				# normalize and select
				mindJe = min(dJe)
				scale  = max(dJe) - mindJe
				for k in xrange(K): dJe[k] = (dJe[k] - mindJe) / scale
				select = select_k(dJe, T)

				if select != res['pos']:
					ct_pert    += 1
					res['pos']  = select

		# sync
		mpi_barrier(MPI_COMM_WORLD)
		ct_pert = mpi_reduce(ct_pert, 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
		ct_pert = mpi_bcast(ct_pert, 1, MPI_INT, main_node, MPI_COMM_WORLD)
		ct_pert = map(int, ct_pert)[0]
		ct_pert /= 2.0

		# select the first temperature if > th
		if ct_pert > th:
			T0 = T
			break

	# sync
	mpi_barrier(MPI_COMM_WORLD)

	# if not found, set to the max value
	if T0 == -1: T0 = Tm
	
	return T0, ct_pert






'''
-- Munkres algorithm (or Hungarian algorithm) ----------------------------------

Copyright and License
=====================

Copyright (c) 2008 Brian M. Clapper

This is free software, released under the following BSD-like license:

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.

  2. The end-user documentation included with the redistribution, if any,
     must include the following acknowlegement:

     This product includes software developed by Brian M. Clapper
     (bmc@clapper.org, http://www.clapper.org/bmc/). That software is
     copyright (c) 2008 Brian M. Clapper.

     Alternately, this acknowlegement may appear in the software itself, if
     and wherever such third-party acknowlegements normally appear.

THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL BRIAN M. CLAPPER BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 


'''
import sys
class Munkres:
    """
Calculate the Munkres solution to the classical assignment problem.
See the module documentation for usage.
    """   
    def __init__(self):
        """Create a new instance"""
        self.C = None
        self.row_covered = []
        self.col_covered = []
        self.n = 0
        self.Z0_r = 0
        self.Z0_c = 0
        self.marked = None
        self.path = None

    def make_cost_matrix(profit_matrix, inversion_function):
        """
        Create a cost matrix from a profit matrix by calling
        'inversion_function' to invert each value. The inversion
        function must take one numeric argument (of any type) and return
        another numeric argument which is presumed to be the cost inverse
        of the original profit.

        This is a static method. Call it like this:

        cost_matrix = Munkres.make_cost_matrix(matrix, inversion_func)

        For example:

        cost_matrix = Munkres.make_cost_matrix(matrix, lambda x : sys.maxint - x)
        """
        cost_matrix = []
        for row in profit_matrix:
            cost_row = []
            for value in row:
                cost_row += [inversion_function(value)]
            cost_matrix += [cost_row]
        return cost_matrix

    make_cost_matrix = staticmethod(make_cost_matrix)

    def compute(self, cost_matrix):
        """
        Compute the indexes for the lowest-cost pairings between rows and
        columns in the database. Returns a list of (row, column) tuples
        that can be used to traverse the matrix.

        The matrix must be square.
        """
        self.C = self.__copy_matrix(cost_matrix)
        self.n = len(cost_matrix)
        self.row_covered = [False for i in range(self.n)]
        self.col_covered = [False for i in range(self.n)]
        self.Z0_r = 0
        self.Z0_c = 0
        self.path = self.__make_matrix(self.n * 2, 0)
        self.marked = self.__make_matrix(self.n, 0)

        done = False
        step = 1

        steps = { 1 : self.__step1,
                  2 : self.__step2,
                  3 : self.__step3,
                  4 : self.__step4,
                  5 : self.__step5,
                  6 : self.__step6 }

        while not done:
            try:
                func = steps[step]
                #print 'calling ' + str(func)
                step = func()
            except KeyError:
                done = True

        # Look for the starred columns
        results = []
        for i in range(self.n):
            for j in range(self.n):
                if self.marked[i][j] == 1:
                    results += [(i, j)]
        assert(len(results) == self.n)

        return results

    def __copy_matrix(self, matrix):
        """Return an exact copy of the supplied matrix"""
        copy = []
        for row in matrix:
            new_row = []
            for item in row:
                new_row += [item]
            copy += [new_row]
        return copy

    def __make_matrix(self, n, val):
        """Create an NxN matrix, populating it with the specific value."""
        matrix = []
        for i in range(n):
            matrix += [[val for j in range(n)]]
        return matrix

    def __step1(self):
        """
        For each row of the matrix, find the smallest element and
        subtract it from every element in its row. Go to Step 2.
        """
        C = self.C
        n = self.n
        for i in range(n):
            minval = self.C[i][0]
            # Find the minimum value for this row
            for j in range(n):
                if minval > self.C[i][j]:
                    minval = self.C[i][j]

            # Subtract that minimum from every element in the row.
            for j in range(n):
                self.C[i][j] -= minval

        return 2

    def __step2(self):
        """
        Find a zero (Z) in the resulting matrix. If there is no starred
        zero in its row or column, star Z. Repeat for each element in the
        matrix. Go to Step 3.
        """
        n = self.n
        for i in range(n):
            for j in range(n):
                if (self.C[i][j] == 0) and \
                   (not self.col_covered[j]) and \
                   (not self.row_covered[i]):
                    self.marked[i][j] = 1
                    self.col_covered[j] = True
                    self.row_covered[i] = True

        self.__clear_covers()
        return 3

    def __step3(self):
        """
        Cover each column containing a starred zero. If K columns are
        covered, the starred zeros describe a complete set of unique
        assignments. In this case, Go to DONE, otherwise, Go to Step 4.
        """
        n = self.n
        count = 0
        for i in range(n):
            for j in range(n):
                if self.marked[i][j] == 1:
                    self.col_covered[j] = True
                    count += 1

        if count >= n:
            step = 7 # done
        else:
            step = 4

        return step

    def __step4(self):
        """
        Find a noncovered zero and prime it. If there is no starred zero
        in the row containing this primed zero, Go to Step 5. Otherwise,
        cover this row and uncover the column containing the starred
        zero. Continue in this manner until there are no uncovered zeros
        left. Save the smallest uncovered value and Go to Step 6.
        """
        step = 0
        done = False
        row = -1
        col = -1
        star_col = -1
        while not done:
            (row, col) = self.__find_a_zero()
            if row < 0:
                done = True
                step = 6
            else:
                self.marked[row][col] = 2
                star_col = self.__find_star_in_row(row)
                if star_col >= 0:
                    col = star_col
                    self.row_covered[row] = True
                    self.col_covered[col] = False
                else:
                    done = True
                    self.Z0_r = row
                    self.Z0_c = col
                    step = 5

        return step

    def __step5(self):
        """
        Construct a series of alternating primed and starred zeros as
        follows. Let Z0 represent the uncovered primed zero found in Step 4.
        Let Z1 denote the starred zero in the column of Z0 (if any).
        Let Z2 denote the primed zero in the row of Z1 (there will always
        be one). Continue until the series terminates at a primed zero
        that has no starred zero in its column. Unstar each starred zero
        of the series, star each primed zero of the series, erase all
        primes and uncover every line in the matrix. Return to Step 3
        """
        count = 0
        path = self.path
        path[count][0] = self.Z0_r
        path[count][1] = self.Z0_c
        done = False
        while not done:
            row = self.__find_star_in_col(path[count][1])
            if row >= 0:
                count += 1
                path[count][0] = row
                path[count][1] = path[count-1][1]
            else:
                done = True

            if not done:
                col = self.__find_prime_in_row(path[count][0])
                count += 1
                path[count][0] = path[count-1][0]
                path[count][1] = col

        self.__convert_path(path, count)
        self.__clear_covers()
        self.__erase_primes()
        return 3

    def __step6(self):
        """
        Add the value found in Step 4 to every element of each covered
        row, and subtract it from every element of each uncovered column.
        Return to Step 4 without altering any stars, primes, or covered
        lines.
        """
        minval = self.__find_smallest()
        for i in range(self.n):
            for j in range(self.n):
                if self.row_covered[i]:
                    self.C[i][j] += minval
                if not self.col_covered[j]:
                    self.C[i][j] -= minval
        return 4

    def __find_smallest(self):
        """Find the smallest uncovered value in the matrix."""
        minval = sys.maxint
        for i in range(self.n):
            for j in range(self.n):
                if (not self.row_covered[i]) and (not self.col_covered[j]):
                    if minval > self.C[i][j]:
                        minval = self.C[i][j]
        return minval

    def __find_a_zero(self):
        """Find the first uncovered element with value 0"""
        row = -1
        col = -1
        i = 0
        n = self.n
        done = False

        while not done:
            j = 0
            while True:
                if (self.C[i][j] == 0) and \
                   (not self.row_covered[i]) and \
                   (not self.col_covered[j]):
                    row = i
                    col = j
                    done = True
                j += 1
                if j >= n:
                    break
            i += 1
            if i >= n:
                done = True

        return (row, col)

    def __find_star_in_row(self, row):
        """
        Find the first starred element in the specified row. Returns
        the column index, or -1 if no starred element was found.
        """
        col = -1
        for j in range(self.n):
            if self.marked[row][j] == 1:
                col = j
                break

        return col

    def __find_star_in_col(self, col):
        """
        Find the first starred element in the specified row. Returns
        the row index, or -1 if no starred element was found.
        """
        row = -1
        for i in range(self.n):
            if self.marked[i][col] == 1:
                row = i
                break

        return row

    def __find_prime_in_row(self, row):
        """
        Find the first prime element in the specified row. Returns
        the column index, or -1 if no starred element was found.
        """
        col = -1
        for j in range(self.n):
            if self.marked[row][j] == 2:
                col = j
                break

        return col

    def __convert_path(self, path, count):
        for i in range(count+1):
            if self.marked[path[i][0]][path[i][1]] == 1:
                self.marked[path[i][0]][path[i][1]] = 0
            else:
                self.marked[path[i][0]][path[i][1]] = 1

    def __clear_covers(self):
        """Clear all covered matrix cells"""
        for i in range(self.n):
            self.row_covered[i] = False
            self.col_covered[i] = False

    def __erase_primes(self):
        """Erase all prime markings"""
        for i in range(self.n):
            for j in range(self.n):
                if self.marked[i][j] == 2:
                    self.marked[i][j] = 0

# NEED TO BE FIX
'''
# Hungarian algorithm between two partitions
def Hungarian(part1, part2):
	from statistics import Munkres
	from numpy      import zeros, array
	import sys

	K = len(part1)
	# prepare matrix
	MAT = [[0] * K for i in xrange(K)]
	for k1 in xrange(K):
		for k2 in xrange(K):
			for index in part1[k1]:
				if index in part2[k2]:
					MAT[k1][k2] += 1
	cost_MAT = []
	for row in MAT:
		cost_row = []
		for col in row:
			cost_row += [sys.maxint - col]
		cost_MAT += [cost_row]

	m = Munkres()
	indexes = m.compute(cost_MAT)

	return indexes, MAT
'''

# K-means main stability stream command line
def k_means_stab_stream(stack, outdir, maskname, K, npart = 5, F = 0, T0 = 0, th_nobj = 0, rand_seed = 0, opt_method = 'cla', CTF = False, maxit = 1e9, flagnorm = False):
	from utilities 	 import print_begin_msg, print_end_msg, print_msg
	from utilities   import model_blank, get_image, get_im, file_type
	from statistics  import k_means_stab_update_tag, k_means_headlog, k_means_export, k_means_init_open_im
	from statistics  import k_means_cla, k_means_SSE, k_means_criterion, k_means_locasg2glbasg, k_means_open_im
	from statistics  import k_means_stab_asg2part, k_means_stab_pwa, k_means_stab_export, k_means_stab_export_txt, k_means_stab_H, k_means_stab_bbenum
	import sys, logging, os, pickle
	
	ext = file_type(stack)
	if ext == 'txt': TXT = True
	else:            TXT = False

	# create a directory
	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', " ", 1)
	os.mkdir(outdir)

	# create main log
	f = open(os.path.join(outdir, 'main_log.txt'), 'w')
	f.close()
		
	logging.basicConfig(filename = os.path.join(outdir, 'main_log.txt'), format = '%(asctime)s     %(message)s', level = logging.INFO)
	logging.info('Clustering')

	# manage random seed
	rnd = []
	if(rand_seed > 0):
		for n in xrange(1, npart + 1): rnd.append(n * (2**n) + rand_seed)
	else:
		from random import randint
		for n in xrange(1, npart + 1): rnd.append(randint(1,91234567))
	logging.info('... Init list random seed: %s' % rnd)

	trials       = 1
	critname     = ''
	logging.info('... K = %03d' % K)

	# open unstable images
	logging.info('... Open images')
	LUT, mask, N, m, Ntot = k_means_init_open_im(stack, maskname)
	IM, ctf, ctf2         = k_means_open_im(stack, mask, CTF, LUT, flagnorm)

	logging.info('... %d unstable images found' % N)
	if N < 2:
		logging.info('[STOP] Not enough images')
		sys.exit()

	# loop over partition
	print_begin_msg('k-means')
	for n in xrange(npart):
		# info
		logging.info('...... Start partition: %d' % (n + 1))
		k_means_headlog(stack, 'partition %d' % (n + 1), opt_method, N, K, critname, maskname, trials, maxit, CTF, T0, F, rnd[n], 1, m)
		
		# classification
		flag_cluster = False
		if   opt_method == 'cla':
			try:			[Cls, assign] = k_means_cla(IM, mask, K, rnd[n], maxit, trials, [CTF, ctf, ctf2], F, T0, False)
			except SystemExit:	flag_cluster  = True
		elif opt_method == 'SSE':
			try:			[Cls, assign] = k_means_SSE(IM, mask, K, rnd[n], maxit, trials, [CTF, ctf, ctf2], F, T0, False)
			except SystemExit:      flag_cluster  = True
		if flag_cluster:
			logging.info('[ERROR] Empty cluster')
			sys.exit()

		# export partition
		crit       = k_means_criterion(Cls, critname)
		glb_assign = k_means_locasg2glbasg(assign, LUT, Ntot)
		k_means_export(Cls, crit, glb_assign, outdir, n, TXT)
				
	# end of classification
	print_end_msg('k-means')

	# convert all assignment to partition
	logging.info('... Matching')
	ALL_PART = k_means_stab_asg2part(outdir, npart)

	# calculate the stability
	if npart == 2:
		stb, nb_stb, STB_PART = k_means_stab_H(ALL_PART)
		logging.info('... Stability: %5.2f %% (%d objects)' % (stb, nb_stb))
	else:
		MATCH, STB_PART, CT_s, CT_t, ST, st = k_means_stab_bbenum(ALL_PART,T=th_nobj)
		logging.info('... Stability: %5.2f %% (%d objects)' % (sum(ST) / float(len(ST)), sum(CT_s)))
	
	if TXT:
		count_k, id_rejected = k_means_stab_export_txt(STB_PART, outdir, th_nobj)
		logging.info('... Export %i stable class averages: averages_grp_i.txt (rejected %i images)' % (count_k, len(id_rejected)))
	else:
		# export the stable class averages
		count_k, id_rejected = k_means_stab_export(STB_PART, stack, outdir, th_nobj, CTF)
		logging.info('... Export %i stable class averages: averages.hdf (rejected %i images)' % (count_k, len(id_rejected)))

		# tag informations to the header
		logging.info('... Update info to the header')
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# k_means_stab_update_tag(stack, STB_PART, id_rejected)
	
	logging.info('... Done')

# K-means main stability stream command line
# added argument num_first_matches (jia)
def k_means_stab_MPI_stream(stack, outdir, maskname, K, npart = 5, F = 0, T0 = 0, th_nobj = 0, rand_seed = 0, opt_method = 'cla', CTF = False, maxit = 1e9, flagnorm = False, num_first_matches=1):
	from mpi         import mpi_init, mpi_comm_size, mpi_comm_rank, mpi_barrier, MPI_COMM_WORLD
	from mpi         import mpi_bcast, MPI_FLOAT, MPI_INT, mpi_send, mpi_recv, MPI_TAG_UB
	from utilities 	 import print_begin_msg, print_end_msg, print_msg, running_time
	from utilities   import model_blank, get_image, get_im, file_type
	from statistics  import k_means_stab_update_tag, k_means_headlog, k_means_init_open_im, k_means_open_im
	from statistics  import k_means_cla_MPI, k_means_SSE_MPI, k_means_criterion, k_means_locasg2glbasg
	from statistics  import k_means_stab_asg2part, k_means_stab_pwa, k_means_stab_export, k_means_stab_H, k_means_export, k_means_stab_export_txt, k_means_stab_bbenum, k_means_stab_getinfo
	from applications import MPI_start_end
	import sys, logging, os, pickle, time
	
	ncpu      = mpi_comm_size(MPI_COMM_WORLD)
	myid      = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	mpi_barrier(MPI_COMM_WORLD)
	npart = ncpu
	ext = file_type(stack)
	if ext == 'txt': TXT = True
	else:            TXT = False

	nx = 0
        if myid == main_node:
		if os.path.exists(outdir):
			nx = 1
			ERROR('Output directory exists, please change the name and restart the program', " k_means_mpi", 0)
		else:
			os.system( "mkdir " + outdir )
	nx = mpi_bcast(nx, 1, MPI_INT, 0, MPI_COMM_WORLD)
	nx = int(nx[0])
	if(nx != 0):
		import sys
		exit()

	mpi_barrier( MPI_COMM_WORLD )

	if myid == main_node:
		# create main log
		f = open(os.path.join(outdir, 'main_log.txt'), 'w')
		f.close()
		
	logging.basicConfig(filename = os.path.join(outdir, 'main_log.txt'), format = '%(asctime)s     %(message)s', level = logging.INFO)
	if myid == main_node: logging.info('Clustering')
	rnd = [ ]
	if(rand_seed > 0):
		for n in xrange( 1,ncpu+1):
			rnd.append( n * (2**n) + rand_seed)
	else:
		from random import randint
		for n in xrange( 1,ncpu+1):
			rnd.append( randint(1,91234567) )
	if myid == main_node: logging.info('... Init list random seed: %s' % rnd)

	trials       = 1
	if myid == main_node: logging.info('... K = %03d' % K)

	# open unstable images
	if myid == main_node: logging.info('... Open images')

	LUT, mask, N, m, Ntot = k_means_init_open_im(stack, maskname)

	IM, ctf, ctf2         = k_means_open_im(stack, mask, CTF, LUT, flagnorm)

	if myid == main_node:
		logging.info('... %d images found' % N)
	if N < K:
		logging.info('[STOP] Not enough images')
		sys.exit()

	# loop over partition
	if myid == main_node: 
		print_begin_msg('k-means')
		t_start = time.time()
	
	[Cls, assign, Je] = k_means_SSE_MPI(IM, mask, K, rnd[myid], maxit,1, [CTF, ctf, ctf2], F, T0, False, "rnd", myid = myid, main_node = main_node, jumping = 0) # no jumping
	
	
	if myid == main_node:
		je_return = [0.0]*(ncpu)
		for n1 in xrange(ncpu):
			if n1 != main_node: je_return[n1]	=	mpi_recv(1,MPI_FLOAT, n1, MPI_TAG_UB, MPI_COMM_WORLD)
 			else:               je_return[main_node]  = Je
	else:
		mpi_send(Je, 1, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)

	
	n_best = 0
	#if any of je , 0, n_bes =-1 and broadcast to all
	if myid == main_node:
		for n in xrange( ncpu):
			je_return[n] = float(je_return[n])
			if( je_return[n] < 0): n_best = -1
	mpi_barrier(MPI_COMM_WORLD)
	if myid == main_node:
		for n1 in xrange(ncpu):
			if n1 != main_node: mpi_send(n_best, 1, MPI_INT, n1, MPI_TAG_UB, MPI_COMM_WORLD) 
 			else:               n_best  = n_best
	else: n_best	=	mpi_recv(1,MPI_INT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)
	n_best = int(n_best)		
	#print "myid ==", myid, "n_best==", n_best
	mpi_barrier( MPI_COMM_WORLD )
	if n_best == -1: 
		if myid == main_node:
			print_msg('>> K is too big and resulted empty cluters, stop k-means  \n' )
		sys.exit()
	
	[r_assign, r_cls] = k_means_SSE_collect(Cls, assign, Je, N, K, ncpu, myid, main_node)
	if myid == main_node:
		for n in xrange( ncpu ):
			r_cls[n]['n'] = map(int, r_cls[n]['n'] )
			r_assign[n] = map(int, r_assign[n] )
			k_means_headlog(stack, outdir, opt_method, N, K, 
						      '', maskname, 1, maxit, CTF, T0, 
						      F, rnd[n], ncpu, m)
			print_msg('\n>>>>>>>>partion:       %5d'%(n+1))
			crit       = k_means_criterion(r_cls[n], 'CHD')
			glb_assign = k_means_locasg2glbasg(r_assign[n], LUT, Ntot)
			k_means_export(r_cls[n], crit, glb_assign, outdir, n, TXT)
			del crit, glb_assign
			
			
	if myid == main_node:
		# end of classification
		running_time(t_start)
		print_end_msg('k-means')

		# convert all assignments to partition
		logging.info('... Matching')
		ALL_PART = k_means_stab_asg2part(outdir, ncpu )

		# calculate the stability
		if npart == 2:
			stb, nb_stb, STB_PART = k_means_stab_H(ALL_PART)
			logging.info('... Stability: %5.2f %% (%d objects)' % (stb, nb_stb))
		# To do the non-mpi version of bbenum......
		else:
			MATCH, STB_PART, CT_s, CT_t, ST, st = k_means_stab_bbenum(ALL_PART,T=th_nobj)
			logging.info('... Stability: %5.2f %% (%d objects)' % (sum(ST) / float(len(ST)), sum(CT_s)))

	if myid == main_node:	
		# export the stable class averages
		if TXT:	count_k, id_rejected = k_means_stab_export_txt(STB_PART, outdir, th_nobj)
		else:   count_k, id_rejected = k_means_stab_export(STB_PART, stack, outdir, th_nobj, CTF)

		logging.info('... Export %i stable class averages: averages.hdf (rejected %i images)' % (count_k, len(id_rejected)))

		if not TXT:
		        # tag informations to the header
			logging.info('... Update info to the header')
			# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1				
			# k_means_stab_update_tag(stack, STB_PART, id_rejected)

		logging.info('... Done')




# Match two partitions with hungarian algorithm
def k_means_match_clusters_asg(asg1, asg2):
	# asg1 and asg2 are numpy array
	from numpy      import zeros, array
	from statistics import Munkres
	import sys

	K        = len(asg1)
	MAT      = [[0] * K for i in xrange(K)] 
	cost_MAT = [[0] * K for i in xrange(K)]
	dummy    = array([0], 'int32')
	for k1 in xrange(K):
		for k2 in xrange(K):
			MAT[k1][k2] = Util.k_means_cont_table(asg1[k1], asg2[k2], dummy, asg1[k1].size, asg2[k2].size, 0)

	for i in xrange(K):
		for j in xrange(K):
			cost_MAT[i][j] = sys.maxint - MAT[i][j]
	m = Munkres()
	indexes = m.compute(cost_MAT)
	list_stable = []
	nb_tot_objs = 0
	for r, c in indexes:
		cont = MAT[r][c]
		if cont == 0:
			list_stable.append(array([], 'int32'))
			continue
		nb_tot_objs += cont
		objs = zeros(cont, 'int32')
		dummy = Util.k_means_cont_table(asg1[r], asg2[c], objs, asg1[r].size, asg2[c].size, 1)
		list_stable.append(objs)

	return list_stable, nb_tot_objs

# Match two partitions with hungarian algorithm and also return the matches whose corresponding stable sets have size larger than T. 
#  Is otherwise identical
# to k_means_match_clusters_asg
# 10/11/11: If the number of elements in common between two classes, e.g., asg1[i] and asg2[j], is not greater than threshold T, 
#     then their "cost" is set to 0 in MAT, i.e., the cost table which is input to Munkres. 
def k_means_match_clusters_asg_new(asg1, asg2, T=0):
	# asg1 and asg2 are numpy array
	from numpy      import zeros, array
	from statistics import Munkres
	import sys

	K        = len(asg1)
	MAT      = [[0] * K for i in xrange(K)] 
	cost_MAT = [[0] * K for i in xrange(K)]
	dummy    = array([0], 'int32')
	for k1 in xrange(K):
		for k2 in xrange(K):
			MAT[k1][k2] = Util.k_means_cont_table(asg1[k1], asg2[k2], dummy, asg1[k1].size, asg2[k2].size, 0)
			if MAT[k1][k2] <= T:
				MAT[k1][k2] = 0
	for i in xrange(K):
		for j in xrange(K):
			cost_MAT[i][j] = sys.maxint - MAT[i][j]
	m = Munkres()
	indexes = m.compute(cost_MAT)
	newindexes = []
	list_stable = []
	nb_tot_objs = 0
	for r, c in indexes:
		cont = MAT[r][c]
		if cont <= T:
			list_stable.append(array([], 'int32'))
			continue
		nb_tot_objs += cont
		objs = zeros(cont, 'int32')
		dummy = Util.k_means_cont_table(asg1[r], asg2[c], objs, asg1[r].size, asg2[c].size, 1)
		list_stable.append(objs)
		newindexes.append([r,c])

	return newindexes, list_stable, nb_tot_objs

# Hierarchical stability between partitions given by k-means
def k_means_stab_H(ALL_PART):
	from copy       import deepcopy
	from statistics import k_means_match_clusters_asg

	nb_part = len(ALL_PART)
	K       = len(ALL_PART[0])
	tot_gbl = 0
	for i in xrange(K): tot_gbl += len(ALL_PART[0][i])
	
	for h in xrange(0, nb_part - 1):
		newPART = []
		for n in xrange(1, nb_part - h):
			LIST_stb, tot_n = k_means_match_clusters_asg(ALL_PART[0], ALL_PART[n])
			newPART.append(LIST_stb)

			nb_stb = 0
			for i in xrange(K): nb_stb += len(LIST_stb[i])


		ALL_PART = []
		ALL_PART = deepcopy(newPART)

	nb_stb = 0
	for i in xrange(K): nb_stb += len(ALL_PART[0][i])
	stability = (float(nb_stb) / float(tot_gbl)) * 100

	STB_PART = []
	for part in ALL_PART[0]: STB_PART.append(map(int, part))

	return stability, nb_stb, STB_PART

# Pairwise recurence agreement matching between paritions given by k-means
def k_means_match_pwa(PART, lim = -1):
	from numpy import zeros, array

	# get table contengenci between two partitions
	def get_mat(part1, part2):
	    K = len(part1)

	    MAT  = zeros((K, K), 'int32')
	    dummy = array([0], 'int32')
	    for k1 in xrange(K):
		    for k2 in xrange(K):
			    MAT[k1][k2] = Util.k_means_cont_table(part1[k1], part2[k2], dummy, part1[k1].size, part2[k2].size, 0)

	    return MAT

	# find alone maximum (maximum along col and same maximum along row)
	def max_g(mat):
            K   = len(mat)
            m1  = mat.argmax(axis=0)
	    m2  = mat.argmax(axis=1)
	    val = []
	    for k in xrange(K):
		    # is a max along the col and the row
		    if m2[m1[k]] == k:
			    val.append([mat[m1[k]][k], m1[k], k])
	    val.sort(reverse = True)
	    # change to flat format [l0, c0, l1, c1, ..., li, ci]
	    res = zeros((2 * len(val)), 'int32')
	    ct  = 0
	    for obj in val:
		    res[ct] = obj[1]
		    res[ct+1] = obj[2]
		    ct += 2

	    return res

	# recursive agreement
	def agree(level, MAX, lmax, np, pos, Nmax, res):
		res[level] = lmax[0]
		if level == (np - 1): return res
		for k in xrange(level + 1, np):
			k2 = lmax[k - level]
			ind1 = mono(level, k)
			ind2 = pos[ind1]
			flag = False
			for i in xrange(Nmax[ind1]):
				i *= 2
				if MAX[ind2 + i] == res[level] and MAX[ind2 + i + 1] == k2:
					flag = True
					lmax[k - level - 1] = k2
					break
			if not flag:
				res[0] = -1
				return res

		res = agree(level + 1, MAX, lmax, np, pos, Nmax, res)

		return res

	# mono 
	def mono(i, j):
		a = max(i, j)
		return min(i, j) + a * (a - 1) // 2

	#====== main ==================
	# prepare table
	np   = len(PART)
	if lim == -1: lim = len(PART[0])               # number of groups
	MAX  = []
	Nmax = zeros((np - 1), 'int32')                # number of maximum per pairwise table
	pos  = zeros((np * (np - 1) / 2 + 1), 'int32') # position list of maximum in MAX
	for i in xrange(1, np):
		for j in xrange(i):
			mat  = get_mat(PART[j], PART[i])
			lmax = max_g(mat)
			nb   = min((len(lmax) // 2), lim)
			lmax = lmax[:2 * nb]
			pos[mono(i, j) + 1] = 2 * nb
			MAX.extend(lmax)

	MAX = array(MAX, 'int32')

	# matching
	Nmax  = pos[1:] 
	Nmax  = Nmax / 2
	pos   = pos.cumsum()
	res   = zeros((np), 'int32')
	lmax  = zeros((np - 1), 'int32')
	MATCH = []
	if np > 2:
		# for each maximum in p0p1
		for k in xrange(0, 2*Nmax[0], 2):
			k1      = MAX[k]
			lmax[0] = MAX[k + 1]
			# search for the same maximum in p0pi
			for pwpart in xrange(1, np):
				flag = False
				for ki in xrange(0, 2*Nmax[pwpart], 2):
					ind1 = pos[mono(pwpart, 0)] + ki
					ki   = MAX[ind1]
					if ki == k1:
						lmax[pwpart - 1] = MAX[ind1 + 1]
						flag = True
						break
				if not flag: break
			# recursive agreement with all remain pairwise tables for the hypotetical list of maximum lmax
			res[0] = k1
			res = agree(1, MAX, lmax, np, pos, Nmax, res)
			if res[0] != -1: MATCH.append(res.copy())
	else:
		# if only two partitions return the list of maximum
		for i in xrange(0, 2*Nmax[0], 2):
			MATCH.append(array([MAX[i], MAX[i+1]], 'int32'))

	return MATCH

# Stability with pairwise agreement matching
def k_means_stab_pwa(PART, lim = -1):
	from statistics import k_means_match_pwa
	from copy       import deepcopy

	MATCH    = k_means_match_pwa(PART, lim)
	np       = len(PART)
	K        = len(PART[0])
	STB_PART = [[] for i in xrange(K)]
	nm       = len(MATCH)
	CT_t     = [0] * K
	CT_s     = [0] * K
	ST       = [0] * K
	ct_t     = 0
	ct_s     = 0
	st       = 0

	for k in xrange(nm):
		kk   = int(MATCH[k][0]) # due to numpy obj
		vmax = [0] * np
		vmin = [0] * np
		for i in xrange(np):
		    vmax[i] = max(PART[i][int(MATCH[k][i])])
		    vmin[i] = min(PART[i][int(MATCH[k][i])])

		vmax = int(max(vmax))
		vmin = int(min(vmin))
		vd   = vmax - vmin + 1

		asg = [0] * vd
		for i in xrange(np):
		    for item in PART[i][int(MATCH[k][i])]: asg[int(item) - vmin] += 1

		stb  = []
		for i in xrange(vd):
		    if asg[i] != 0:
			CT_t[kk] += 1
			if asg[i] == np:
			    CT_s[kk] += 1
			    stb.append(i + vmin)

		STB_PART[kk] = deepcopy(stb)

	for k in xrange(K):
		if CT_t[k] == 0: continue
		ST[k] = 100.0 * CT_s[k] / float(CT_t[k])

	if sum(CT_t) == 0:
		st = 0
	else:   st = 100.0 * sum(CT_s) / float(sum(CT_t))

	return MATCH, STB_PART, CT_s, CT_t, ST, st

# Export stable averages to text file
def k_means_stab_export_txt(PART, outdir, th_nobj):
	K  = len(PART)
	RK = 0
	for k in xrange(K):
		if len(PART[k]) >= th_nobj:
			f = open(outdir + '/stable_grp_%03i.txt' % k, 'w')
			for id in PART[k]: f.write('%i\n' % id)
			f.close()
			RK += 1
	return RK, []

# Build and export the stable class averages 
def k_means_stab_export(PART, stack, outdir, th_nobj, CTF = False):
	from utilities    import model_blank, get_params2D, get_im
	from fundamentals import rot_shift2D, fftip
	
	K    = len(PART)
	im   = EMData()
	im.read_image(stack, 0)
	nx   = im.get_xsize()
	ny   = im.get_ysize()
	imbk = model_blank(nx, ny)
	AVE  = []
	lrej = []
	Kr   = 0
	ck   = 0
	for k in xrange(K):
		if len(PART[k]) >= th_nobj:
			Kr += 1
	if Kr == 0:
		imbk.write_image(outdir + '/averages.hdf', 0)
		return 0, []
	for k in xrange(Kr): AVE.append(imbk.copy())
	for k in xrange(K):
		nobjs = len(PART[k])
		if nobjs >= th_nobj:
			if CTF:
				data = []
				for ID in PART[k]: data.append(get_im(stack, int(ID)))
				AVE[ck], dum, dum, dum, dum = add_ave_varf(data, None, 'a', True)
				fftip(AVE[ck])
				AVE[ck].depad()
			else:
				for ID in PART[k]:
					im.read_image(stack, int(ID))
					if im.get_ysize() > 1:
						alpha, sx, sy, mirror, scale = get_params2D(im)
						im = rot_shift2D(im, alpha, sx, sy, mirror)

					Util.add_img(AVE[ck], im)
				Util.mul_scalar(AVE[ck], 1.0 / float(nobjs))

			AVE[ck].set_attr('Class_average', 1.0)
			AVE[ck].set_attr('nobjects', nobjs)
			AVE[ck].set_attr('members', PART[k])
			AVE[ck].set_attr('k_ref', k)
			AVE[ck].write_image(outdir + '/averages.hdf', ck)
			ck += 1
		else:
			lrej.extend(PART[k])

	return ck, lrej

# Init the header for the stack file
# TODO this function need to be removed (not used)
def k_means_stab_init_tag(stack):
	from utilities import file_type, write_header
	from EMAN2db import db_open_dict
	
	N   = EMUtil.get_image_count(stack)
	ext = file_type(stack)
	if ext == 'bdb':
		DB = db_open_dict(stack)
		for n in xrange(N):
			DB.set_attr(n, 'stab_active', 1)
			DB.set_attr(n, 'stab_part', -2)
		DB.close()
	else:
		im = EMData()
		for n in xrange(N):
			im.read_image(stack, n, True)
			im.set_attr('stab_active', 1)
			im.set_attr('stab_part', -2)
			write_header(stack, im, n) 

# Convert local all assignment to absolute all partition
def k_means_stab_asg2part(outdir, npart):
	from numpy     import array
	from utilities import get_im
	from os        import path, listdir

	# first check special case, when membership is export as txt file
	# due to the limitation of hdf header or the data is TXT file
	ALL_PART = []
	if path.isfile(path.join(outdir, 'kmeans_part_01_grp_001.txt')):
		from utilities import read_text_file
		lname = listdir(outdir)
		K     = 0
		for name in lname:
			if name.find('kmeans_part_01_grp_') == 0: K += 1
		for n in xrange(npart):
			part = []
			for k in xrange(K):
				lid = read_text_file(path.join(outdir, 'kmeans_part_%02i_grp_%03i.txt' % (n+1, k+1)), 0)
				lid = array(lid, 'int32')
				lid.sort()
				part.append(lid.copy())
			ALL_PART.append(part)
	else:
		for n in xrange(npart):
			name = path.join(outdir, 'averages_%02i.hdf' % n)			
			K    = EMUtil.get_image_count(name)
			part = []
			for k in xrange(K):
				im  = get_im(name, k)
				lid = im.get_attr('members')
				lid = array(lid, 'int32')
				lid.sort()
				part.append(lid.copy())
			ALL_PART.append(part)

	return ALL_PART

# Convert local assignment to absolute partion
def k_means_asg_locasg2glbpart(ASG, LUT):
	K = max(ASG) + 1
	N = len(ASG)
	PART = [[] for i in xrange(K)]
	for n in xrange(N): PART[ASG[n]].append(LUT[n])

	return PART

# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
# # Update information to the header of the stack file
# def k_means_stab_update_tag(stack, STB_PART, lrej):
# 	from utilities import file_type, write_header
# 	from EMAN2db import db_open_dict
# 
# 	N  = EMUtil.get_image_count(stack)
# 	# prepare for active images
# 	list_stb = []
# 	for part in STB_PART: list_stb.extend(part)
# 
# 	ext = file_type(stack)
# 	if ext == 'bdb':
# 		DB = db_open_dict(stack)
# 		# these are stable
# 		for ID in list_stb: DB.set_attr(ID, 'active', 0)
# 		# these are still unstable
# 		for ID in lrej: DB.set_attr(ID, 'active', 1)
# 		DB.close()
# 	else:
# 		im = EMData()
# 		# these are stable
# 		for ID in list_stb:
# 			im.read_image(stack, ID, True)
# 			im.set_attr('active', 0)
# 			write_header(stack, im, ID)
# 		# these are still unstable
# 		for ID in lrej:
# 			im.read_image(stack, ID, True)
# 			im.set_attr('active', 1)
# 			write_header(stack, im, ID)

# Gather all stable class averages in the same stack
def k_means_stab_gather(nb_run, maskname, outdir):
	from utilities import get_image
	from os        import path

	ct   = 0
	im   = EMData()
	for nr in xrange(1, nb_run + 1):
		name = outdir + '/average_stb_run%02d.hdf' % nr
		if path.exists(name):
			N = EMUtil.get_image_count(name)
			if nr == 1:
				if maskname != None: mask = get_image(maskname)
				else: mask   = None
			for n in xrange(N):
				im.read_image(name, n)
				try:
					nobjs = im.get_attr('nobjects') # check if ave not empty
					ret = Util.infomask(im, mask, False) # 
					im  = (im - ret[0]) / ret[1]        # normalize
					im.write_image(outdir + '/averages.hdf', ct)
					ct += 1
				except: pass

	return ct

# extract group to a stack of images for each classe, and apply alignment
def k_means_extract_class_ali(stack_name, ave_name, dir):
	from   utilities  import get_im, get_params2D, set_params2D
	from fundamentals import rot_shift2D
	import os

	K = EMUtil.get_image_count(ave_name)
	# all images are not open, I assume we pick up only few images (from stable_averages)
	for k in xrange(K):
		im  = get_im(ave_name, k)
		lid = im.get_attr('members')
		# if there are only one image in member 
		if not isinstance(lid, list): lid = [lid]
		trg = os.path.join(dir, 'class_%03i.hdf' % k)
		for i, item in enumerate(lid):
			im = get_im(stack_name, item)
			alpha, sx, sy, mir, scale = get_params2D(im, 'xform.align2d')
			im = rot_shift2D(im, alpha, sx, sy, mir, scale)
			set_params2D(im, [0.0, 0.0, 0.0, 0, 1.0], 'xform.align2d')
			# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
			# im.set_attr('active', 1)
			im.write_image(trg, i)
	
# compute pixel error for a class given	
def k_means_class_pixerror(class_name, dir, ou, xr, ts, maxit, fun, CTF=False, snr=1.0, Fourvar=False):
	from applications import header, ali2d
	from utilities    import estimate_stability
	from statistics   import aves
	import os
	name = class_name.split('.')[0]
	file = os.path.join(dir, class_name)
	header(file, 'xform.align2d', randomize=True)
	ali2d(file, os.path.join(dir, '%s_01' % name), ou=ou, xr=xr, ts=ts, maxit=maxit,
		CTF=CTF, snr=snr, Fourvar=Fourvar, user_func_name=fun, MPI=False)
	header(file, 'xform.align2d', backup=True, suffix='_round1')
	header(file, 'xform.align2d', randomize=True)
	ali2d(file, os.path.join(dir, '%s_02' % name), ou=ou, xr=xr, ts=ts, maxit=maxit,
		CTF=CTF, snr=snr, Fourvar=Fourvar, user_func_name=fun, MPI=False)

	data1 = EMData.read_images(file)
	header(file, 'xform.align2d_round1', restore=True)
	data2 = EMData.read_images(file)

	stab_mirror, list_pix_err, ccc = estimate_stability(data1, data2, CTF, snr, ou)
	ave_pix_err = sum(list_pix_err) / float(len(list_pix_err))

	ave, var = aves(file, 'a')
	ave.set_attr('err_mir', stab_mirror)
	ave.set_attr('err_pix', ave_pix_err)
	ave.set_attr('err_ccc', ccc)

	return ave

# ISC procedure, update configuration file with ite
def isc_update_ite_conf(conf_file, ite):
	import ConfigParser
	config = ConfigParser.ConfigParser()
	config.read(conf_file)
	config.set('main', 'ite', ite)
	config.write(open(conf_file, 'w'))

# ISC procedure, read configuration file
def isc_read_conf(conf_file):
	import ConfigParser

	# read config file
	config  = ConfigParser.ConfigParser()
	config.read(conf_file)
	cfgmain = dict(config.items('main'))
	cfgali  = dict(config.items('alignment'))
	cfgclu  = dict(config.items('clustering'))
	cfgrali = dict(config.items('realignment'))
	# convert data if need

	#   alignment
	#cfgali['n_ite']      = int(cfgali['n_ite'])
	try:	cfgali['fourvar']    = eval(cfgali['fourvar'])
	except:	cfgali['fourvar']    = False
	cfgali['maxit']      = int(cfgali['maxit'])
	cfgali['maxit']      = int(cfgali['maxit'])
	try:	cfgali['ctf']        = eval(cfgali['ctf'])
	except:	cfgali['ctf']        = False
	try:	cfgali['snr']        = float(cfgali['snr'])
	except:	cfgali['snr']        = 1.0
	cfgali['ou']         = int(cfgali['ou'])
	cfgali['nb_cpu']     = int(cfgali['nb_cpu'])
	try:       cfgali['cuda']       = eval(cfgali['cuda'])
	except:     cfgali['cuda']       = False
	try:       cfgali['ng'] 	     = int(cfgali['ng'])
	except:   cfgali['ng'] 	     = -1
	try:	cfgali['dst']	     = eval(cfgali['dst'])	
	except:	cfgali['dst']	     = 0.0	
	try:	cfgali['center']     = int(cfgali['center'])
	except:	cfgali['center']     = -1

	#   clustering
	try:      cfgclu['f']          = float(cfgclu['f'])
	except:    cfgclu['f']          = 0.9
	try:	cfgclu['th_nobj']    = int(cfgclu['th_nobj'])
	except:	cfgclu['th_nobj']    = 1
	try:      cfgclu['t0']         = float(cfgclu['t0'])
	except:    cfgclu['t0']         = 2.0
	try:	cfgclu['ctf']        = eval(cfgclu['ctf'])
	except:	cfgclu['ctf']        = False
	cfgclu['maxit']      = int(cfgclu['maxit'])
	cfgclu['rand_seed']  = int(cfgclu['rand_seed'])
	cfgclu['im_per_grp'] = int(cfgclu['im_per_grp'])
	cfgclu['nb_part']    = int(cfgclu['nb_part'])
	cfgclu['nb_cpu']     = int(cfgclu['nb_cpu'])
	try:	cfgclu['flag_tiny']	= eval(cfgclu['flag_tiny'])
	except:	cfgclu['flag_tiny']	= False
	try:	cfgclu['cuda']	= eval(cfgclu['cuda'])
	except:	cfgclu['cuda']	= False
	try:	cfgclu['filter_cutoff']        = float(cfgclu['filter_cutoff'])
	except:	cfgclu['filter_cutoff']        = 0.0
	try:	cfgclu['filter_falloff']        = float(cfgclu['filter_falloff'])
	except:	cfgclu['filter_falloff']        = 0.0
	try:	cfgclu['new_nx']        = int(cfgclu['new_nx'])
	except:	cfgclu['new_nx']        = -1
	try:	cfgclu['new_ny']        = int(cfgclu['new_ny'])
	except:	cfgclu['new_ny']        = -1

	#   realignment
	cfgrali['fourvar']   = eval(cfgrali['fourvar'])
	cfgrali['maxit']     = int(cfgrali['maxit'])
	cfgrali['ctf']       = eval(cfgrali['ctf'])
	try:    cfgrali['snr']       = float(cfgrali['snr'])
	except:  cfgrali['snr']       = 1.0
	cfgrali['ou']        = int(cfgrali['ou'])
	cfgrali['nb_cpu']    = int(cfgrali['nb_cpu'])
	try:	cfgrali['cuda']	     = eval(cfgrali['cuda'])
	except:	cfgrali['cuda']	     = False
	try:	cfgrali['ng']	     = int(cfgrali['ng'])
	except:	cfgrali['ng']	     = -1
	try:	cfgrali['num_ali']   = int(cfgrali['num_ali'])
	except:	cfgrali['num_ali']   = 4
	try:	cfgrali['th_mir']    = float(cfgrali['th_mir'])
	except:	cfgrali['th_mir']    = 0.5
	try:	cfgrali['th_err']    = float(cfgrali['th_err'])
	except:	cfgrali['th_err']    = 1.0
	try:	cfgrali['center']    = int(cfgrali['center'])
	except:	cfgrali['center']    = -1
	try:	cfgrali['dst']       = float(cfgrali['dst'])
	except:	cfgrali['dst']       = 0.0

	#   main
	cfgmain['ite']       = int(cfgmain['ite'])
	cfgmain['maxit']     = int(cfgmain['maxit'])

	return cfgmain, cfgali, cfgclu, cfgrali

def isc_ave_huge(ave_tiny, org_data, ave_huge):
	from utilities    import get_im, get_params2D, model_blank
	from fundamentals import rot_shift2D
	from sys import exit
	im = get_im(org_data, 0)
	nx = im.get_xsize()
	ny = im.get_ysize()
	D  = EMData.read_images(org_data)
	A  = EMData.read_images(ave_tiny)
	K  = len(A)
	for k in xrange(K):
		asg = A[k].get_attr('members')
		asg = map(int, asg)
		ave = model_blank(nx, ny)
		for id in asg:
			a, sx, sy, mir, sc = get_params2D(D[id])
			im = rot_shift2D(D[id], a, sx, sy, mir, sc)
			Util.add_img(ave, im)
		Util.mul_scalar(ave, 1.0 / float(len(asg)))

		ave.set_attr('members', asg)
		ave.set_attr('nobjects', len(asg))
		ave.set_attr('Class_average', 1)
		ave.write_image(ave_huge, k)





	

### END K-MEANS ##############################################################################
##############################################################################################

##############################################################################################
### PY CLUSTER ###############################################################################
# 2008-12-08 12:39:54 JB
#
# This is part of "python-cluster". A library to group similar items together.
# Copyright (C) 2006   Michel Albert
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#

def py_cluster_median(numbers):
   """
   Return the median of the list of numbers.
   found at: http://mail.python.org/pipermail/python-list/2004-December/253517.html
   """
   # Sort the list and take the middle element.
   n = len(numbers)
   copy = numbers[:] # So that "numbers" keeps its original order
   copy.sort()
   if n & 1:         # There is an odd number of elements
      return copy[n // 2]
   else:
      return (copy[n // 2 - 1] + copy[n // 2]) / 2.0

def py_cluster_mean(numbers):
   """
   Returns the arithmetic mean of a numeric list.
   found at: http://mail.python.org/pipermail/python-list/2004-December/253517.html
   """
   return float(sum(numbers)) / float(len(numbers))

def py_cluster_genmatrix(list, combinfunc, symmetric=False, diagonal=None):
   """
   Takes a list and generates a 2D-matrix using the supplied combination
   function to calculate the values.

   PARAMETERS
      list        - the list of items
      combinfunc  - the function that is used to calculate teh value in a cell.
                    It has to cope with two arguments.
      symmetric   - Whether it will be a symmetric matrix along the diagonal.
                    For example, it the list contains integers, and the
                    combination function is abs(x-y), then the matrix will be
                    symmetric.
                    Default: False
      diagonal    - The value to be put into the diagonal. For some functions,
                    the diagonal will stay constant. An example could be the
                    function "x-y". Then each diagonal cell will be "0".
                    If this value is set to None, then the diagonal will be
                    calculated.
                    Default: None
   """
   matrix = []
   row_index = 0
   for item in list:
      row = []
      col_index = 0
      for item2 in list:
         if diagonal is not None and col_index == row_index:
            # if this is a cell on the diagonal
            row.append(diagonal)
         elif symmetric and col_index < row_index:
            # if the matrix is symmetric and we are "in the lower left triangle"
            row.append( matrix[col_index][row_index] )
         else:
            # if this cell is not on the diagonal
            row.append(combinfunc(item, item2))
         col_index += 1
      matrix.append(row)
      row_index += 1
   return matrix

class py_Cluster:
   """
   A collection of items. This is internally used to detect clustered items in
   the data so we could distinguish other collection types (lists, dicts, ...)
   from the actual clusters. This means that you could also create clusters of
   lists with this class.
   """

   def __str__(self):
      return "<Cluster@%s(%s)>" % (self.__level, self.__items)

   def __repr__(self):
      return self.__str__()

   def __init__(self, level, *args):
      """
      Constructor

      PARAMETERS
         level - The level of this cluster. This is used in hierarchical
                 clustering to retrieve a specific set of clusters. The higher
                 the level, the smaller the count of clusters returned. The
                 level depends on the difference function used.
         *args - every additional argument passed following the level value
                 will get added as item to the cluster. You could also pass a
                 list as second parameter to initialise the cluster with that
                 list as content
      """
      self.__level = level
      if len(args) == 0: self.__items = []
      else:              self.__items = list(args)

   def append(self, item):
      """
      Appends a new item to the cluster

      PARAMETERS
         item  -  The item that is to be appended
      """
      self.__items.append(item)

   def items(self, newItems = None):
      """
      Sets or gets the items of the cluster

      PARAMETERS
         newItems (optional) - if set, the items of the cluster will be
                               replaced with that argument.
      """
      if newItems is None: return self.__items
      else:                self.__items = newItems

   def fullyflatten(self, *args):
      """
      Completely flattens out this cluster and returns a one-dimensional list
      containing the cluster's items. This is useful in cases where some items
      of the cluster are clusters in their own right and you only want the
      items.

      PARAMETERS
         *args - only used for recursion.
      """
      flattened_items = []
      if len(args) == 0: collection = self.__items
      else:              collection = args[0].items()

      for item in collection:
         if isinstance(item, py_Cluster):
            flattened_items = flattened_items + self.fullyflatten(item)
         else:
            flattened_items.append(item)

      return flattened_items

   def level(self):
      """
      Returns the level associated with this cluster
      """
      return self.__level

   def display(self, depth=0):
      """
      Pretty-prints this cluster. Useful for debuging
      """
      print depth*"   " + "[level %s]" % self.__level
      for item in self.__items:
         if isinstance(item, py_Cluster):
            item.display(depth+1)
         else:
            print depth*"   "+"%s" % item

   def topology(self):
      """
      Returns the structure (topology) of the cluster as tuples.

      Output from cl.data:
          [<Cluster@0.833333333333(['CVS', <Cluster@0.818181818182(['34.xls',
          <Cluster@0.789473684211([<Cluster@0.555555555556(['0.txt',
          <Cluster@0.181818181818(['ChangeLog', 'ChangeLog.txt'])>])>,
          <Cluster@0.684210526316(['20060730.py',
          <Cluster@0.684210526316(['.cvsignore',
          <Cluster@0.647058823529(['About.py',
          <Cluster@0.625(['.idlerc', '.pylint.d'])>])>])>])>])>])>])>]

      Corresponding output from cl.topo():
          ('CVS', ('34.xls', (('0.txt', ('ChangeLog', 'ChangeLog.txt')),
          ('20060730.py', ('.cvsignore', ('About.py',
          ('.idlerc', '.pylint.d')))))))
      """

      left  = self.__items[0]
      right = self.__items[1]
      if isinstance(left, py_Cluster):
          first = left.topology()
      else:
          first = left
      if isinstance(right, py_Cluster):
          second = right.topology()
      else:
          second = right
      return first, second

   def getlevel(self, threshold):
      """
      Retrieve all clusters up to a specific level threshold. This
      level-threshold represents the maximum distance between two clusters. So
      the lower you set this threshold, the more clusters you will receive and
      the higher you set it, you will receive less but bigger clusters.

      PARAMETERS
         threshold - The level threshold

      NOTE
         It is debatable whether the value passed into this method should
         really be as strongly linked to the real cluster-levels as it is right
         now. The end-user will not know the range of this value unless s/he
         first inspects the top-level cluster. So instead you might argue that
         a value ranging from 0 to 1 might be a more useful approach.
      """

      left  = self.__items[0]
      right = self.__items[1]

      # if this object itself is below the threshold value we only need to
      # return it's contents as a list
      if self.level() <= threshold:
         return [self.fullyflatten()]

      # if this cluster's level is higher than the threshold we will investgate
      # it's left and right part. Their level could be below the threshold
      if isinstance(left, py_Cluster) and left.level() <= threshold:
         if isinstance(right, py_Cluster):
            return [left.fullyflatten()] + right.getlevel(threshold)
         else:
            return [left.fullyflatten()] + [[right]]
      elif isinstance(right, py_Cluster) and right.level() <= threshold:
         if isinstance(left, py_Cluster):
            return left.getlevel(threshold) + [right.fullyflatten()]
         else:
            return [[left]] + [right.fullyflatten()]

      # Alright. We covered the cases where one of the clusters was below the
      # threshold value. Now we'll deal with the clusters that are above by
      # recursively applying the previous cases.
      if isinstance(left, py_Cluster) and isinstance(right, py_Cluster):
         return left.getlevel(threshold) + right.getlevel(threshold)
      elif isinstance(left, py_Cluster):
         return left.getlevel(threshold) + [[right]]
      elif isinstance(right, py_Cluster):
         return [[left]] + right.getlevel(threshold)
      else:
         return [[left], [right]]

class py_cluster_BaseClusterMethod:
   """
   The base class of all clustering methods.
   """

   def __init__(self, input, distance_function):
      """
      Constructs the object and starts clustering

      PARAMETERS
         input             - a list of objects
         distance_function - a function returning the distance - or opposite of
                             similarity ( distance = -similarity ) - of two
                             items from the input. In other words, the closer
                             the two items are related, the smaller this value
                             needs to be. With 0 meaning they are exactly the
                             same.

      NOTES
         The distance function should always return the absolute distance
         between two given items of the list. Say,

         distance(input[1], input[4]) = distance(input[4], input[1])

         This is very important for the clustering algorithm to work!
         Naturally, the data returned by the distance function MUST be a
         comparable datatype, so you can perform arithmetic comparisons on
         them (< or >)! The simplest examples would be floats or ints. But as
         long as they are comparable, it's ok.
      """
      self.distance = distance_function
      self._input = input    # the original input
      self._data  = input[:] # clone the input so we can work with it

   def topo(self):
      """
      Returns the structure (topology) of the cluster.

      See Cluster.topology() for information.
      """
      return self.data[0].topology()

   def __get_data(self):
      """
      Returns the data that is currently in process.
      """
      return self._data
   data = property(__get_data)

   def __get_raw_data(self):
      """
      Returns the raw data (data without being clustered).
      """
      return self._input
   raw_data = property(__get_raw_data)

class py_cluster_HierarchicalClustering(py_cluster_BaseClusterMethod):
   """
   Implementation of the hierarchical clustering method as explained in
   http://www.elet.polimi.it/upload/matteucc/Clustering/tutorial_html/hierarchical.html

   USAGE
      >>> from cluster import HierarchicalClustering
      >>> # or: from cluster import *
      >>> cl = HierarchicalClustering([123,334,345,242,234,1,3], lambda x,y: float(abs(x-y)))
      >>> cl.getlevel(90)
      [[345, 334], [234, 242], [123], [3, 1]]

      Note that all of the returned clusters are more that 90 apart

   """

   def __init__(self, data, distance_function, linkage='single'):
      """
      Constructor

      See BaseClusterMethod.__init__ for more details.
      """
      py_cluster_BaseClusterMethod.__init__(self, data, distance_function)

      # set the linkage type to single
      self.setLinkageMethod(linkage)
      self.__clusterCreated = False

   def setLinkageMethod(self, method):
      """
      Sets the method to determine the distance between two clusters.

      PARAMETERS:
         method - The name of the method to use. It must be one of 'single',
                  'complete', 'average' or 'uclus'
      """
      if method == 'single':
         self.linkage = self.singleLinkageDistance
      elif method == 'complete':
         self.linkage = self.completeLinkageDistance
      elif method == 'average':
         self.linkage = self.averageLinkageDistance
      elif method == 'uclus':
         self.linkage = self.uclusDistance
      else:
         raise ValueError, 'distance method must be one of single, complete, average of uclus'

   def uclusDistance(self, x, y):
      """
      The method to determine the distance between one cluster an another
      item/cluster. The distance equals to the *average* (median) distance from
      any member of one cluster to any member of the other cluster.

      PARAMETERS
         x  -  first cluster/item
         y  -  second cluster/item
      """
      # create a flat list of all the items in <x>
      if not isinstance(x, py_Cluster): x = [x]
      else: x = x.fullyflatten()

      # create a flat list of all the items in <y>
      if not isinstance(y, py_Cluster): y = [y]
      else: y = y.fullyflatten()

      distances = []
      for k in x:
         for l in y:
            distances.append(self.distance(k,l))
      return py_cluster_median(distances)

   def averageLinkageDistance(self, x, y):
      """
      The method to determine the distance between one cluster an another
      item/cluster. The distance equals to the *average* (mean) distance from
      any member of one cluster to any member of the other cluster.

      PARAMETERS
         x  -  first cluster/item
         y  -  second cluster/item
      """
      # create a flat list of all the items in <x>
      if not isinstance(x, py_Cluster): x = [x]
      else: x = x.fullyflatten()

      # create a flat list of all the items in <y>
      if not isinstance(y, py_Cluster): y = [y]
      else: y = y.fullyflatten()

      distances = []
      for k in x:
         for l in y:
            distances.append(self.distance(k,l))
      return py_cluster_mean(distances)

   def completeLinkageDistance(self, x, y):
      """
      The method to determine the distance between one cluster an another
      item/cluster. The distance equals to the *longest* distance from any
      member of one cluster to any member of the other cluster.

      PARAMETERS
         x  -  first cluster/item
         y  -  second cluster/item
      """

      # create a flat list of all the items in <x>
      if not isinstance(x, py_Cluster): x = [x]
      else: x = x.fullyflatten()

      # create a flat list of all the items in <y>
      if not isinstance(y, py_Cluster): y = [y]
      else: y = y.fullyflatten()

      # retrieve the minimum distance (single-linkage)
      maxdist = self.distance(x[0], y[0])
      for k in x:
         for l in y:
            maxdist = max(maxdist, self.distance(k,l))

      return maxdist

   def singleLinkageDistance(self, x, y):
      """
      The method to determine the distance between one cluster an another
      item/cluster. The distance equals to the *shortest* distance from any
      member of one cluster to any member of the other cluster.

      PARAMETERS
         x  -  first cluster/item
         y  -  second cluster/item
      """

      # create a flat list of all the items in <x>
      if not isinstance(x, py_Cluster): x = [x]
      else: x = x.fullyflatten()

      # create a flat list of all the items in <y>
      if not isinstance(y, py_Cluster): y = [y]
      else: y = y.fullyflatten()

      # retrieve the minimum distance (single-linkage)
      mindist = self.distance(x[0], y[0])
      for k in x:
         for l in y:
            mindist = min(mindist, self.distance(k,l))

      return mindist

   def cluster(self, matrix=None, level=None, sequence=None):
      """
      Perform hierarchical clustering. This method is automatically called by
      the constructor so you should not need to call it explicitly.

      PARAMETERS
         matrix   -  The 2D list that is currently under processing. The matrix
                     contains the distances of each item with each other
         level    -  The current level of clustering
         sequence -  The sequence number of the clustering
      """

      if matrix is None:
         # create level 0, first iteration (sequence)
         level    = 0
         sequence = 0
         matrix   = []

      # if the matrix only has two rows left, we are done
      while len(matrix) > 2 or matrix == []:

         matrix = py_cluster_genmatrix(self._data, self.linkage, True, 0)

         smallestpair = None
         mindistance  = None
         rowindex = 0   # keep track of where we are in the matrix
         # find the minimum distance
         for row in matrix:
            cellindex = 0 # keep track of where we are in the matrix
            for cell in row:
               # if we are not on the diagonal (which is always 0)
               # and if this cell represents a new minimum...
               if (rowindex != cellindex) and ( cell < mindistance or smallestpair is None ):
                  smallestpair = ( rowindex, cellindex )
                  mindistance  = cell
               cellindex += 1
            rowindex += 1

         sequence += 1
         level     = matrix[smallestpair[1]][smallestpair[0]]
         cluster   = py_Cluster(level, self._data[smallestpair[0]], self._data[smallestpair[1]])

         # maintain the data, by combining the the two most similar items in the list
         # we use the min and max functions to ensure the integrity of the data.
         # imagine: if we first remove the item with the smaller index, all the
         # rest of the items shift down by one. So the next index will be
         # wrong. We could simply adjust the value of the second "remove" call,
         # but we don't know the order in which they come. The max and min
         # approach clarifies that
         self._data.remove(self._data[max(smallestpair[0], smallestpair[1])]) # remove item 1
         self._data.remove(self._data[min(smallestpair[0], smallestpair[1])]) # remove item 2
         self._data.append(cluster)               # append item 1 and 2 combined

      # all the data is in one single cluster. We return that and stop
      self.__clusterCreated = True
      return

   def getlevel(self, threshold):
      """
      Returns all clusters with a maximum distance of <threshold> in between
      each other

      PARAMETERS
         threshold - the maximum distance between clusters

      SEE-ALSO
         Cluster.getlevel(threshold)
      """

      # if it's not worth clustering, just return the data
      if len(self._input) <= 1: return self._input

      # initialize the cluster if not yet done
      if not self.__clusterCreated: self.cluster()

      return self._data[0].getlevel(threshold)

### END PY_CLUSTER ###########################################################################
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
	from alignment    import   ang_n
	from development  import   oned_search_func
	from random       import   gauss
	from EMAN2 import Processor
	
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
	Util.Applyws(tave_w, numr, wr)
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
		Util.Applyws(tave_w, numr, wr)
	
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
			Util.Applyws(refim_w, numr, wr)
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

'''
# This code looks obsolete, I do not write it. JB.. So I commented it out, PAP 12/31/09
def k_means_aves(images, N, K, rand_seed, outdir):
	from utilities import model_blank
	from random       import seed, randint
	from sys import exit
	
	#  This version is not quite correct, as ctf2 is not modified in the average, but it cannot be done any other simple way :-(
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
				print  "   #####################    MOVING OBJECT",im,Cls,CRIT
				CT += Gain[assign_to]
				assign[im] = assign_to
				Cls[ca] -= 1
				CRIT[ca] = JMINUS
				Cls[assign_to] += 1	
				CRIT[assign_to] = JPLUSk
				print  "   #####################    MOVING OBJECT",im,assign_to,CT, Cls,CRIT
				change = True
		print " ITERATION ",cnt
		print  CT
		print  CRIT
		print  assign
	exit()
'''

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

def histogram2d(datai, dataj, nbini, nbinj):
	fmaxi = max(datai)
	fmini = min(datai)
	fmaxj = max(dataj)
	fminj = min(dataj)

	binsize_i = (fmaxi-fmini)/nbini
	binsize_j = (fmaxj-fminj)/nbinj
	start_i = fmini
	start_j = fminj

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
		idd = min(int((datai[k]-start_i)/binsize_i), nbini-1) 
		jdd = min(int((dataj[k]-start_j)/binsize_j), nbinj-1)
		hist[idd][jdd]+=1

	return region,hist

def hist_list(data, nbin = -1, fminiu = None, fmaxiu = None):
	"""
	  Calculate histogram of the list elements
	  nbin - number of bins, if not provided it will be set such that in average there is 10 elements per bin
	  fminiu - user provided minimum value for the histogram, it has to be smaller than the smallest element in data
	  fmaxiu - user provided maximum value for the histogram, it has to be larger than the largest element in data
	"""
	if nbin < 0:  nbin = len(data)/10
	fmaxi = max(data)
	fmini = min(data)

	if fmaxi == fmini:
		hist = [0]*nbin
		hist[0] = len(data)
		return [fmaxi]*nbin, hist
	if fminiu != None:
		if fminiu < fmini : fmini = fminiu
	if fmaxiu != None:
		if fmaxiu > fmaxi : fmaxi = fmaxiu

	binsize_i = (fmaxi-fmini)/float(nbin)
	start_i = fmini

	region = [None]*nbin
	hist = [None]*nbin
	for i in xrange(nbin):
		region[i] = start_i + i*binsize_i
		hist[i] = 0

	for d in data:
		i = min(int((d-start_i)/binsize_i), nbin-1)
		hist[i] += 1

	return region, hist

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
	  Pearson correlation coefficient between two lists
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
	  Basic statistics of numbers stored in a list: average, variance, minimum, maximum
	"""
	av = X[0]
	va = X[0]*X[0]
	mi = X[0]
	ma = X[0]
	N = len(X)
	for i in xrange(1,N):
		av += X[i]
		va += X[i]*X[i]
		mi = min(mi, X[i])
		ma = max(ma, X[i])
	return  av/N,(va - av*av/N)/float(N-1) , mi, ma

def get_power_spec(stack_file, start_particle, end_particle):
	"""
		Name
			get_power_spec - computes the rotationally averaged power spectra of a series of images
		Input
			stack_file: stack images in hdf format
			start_particle: first particle to use
			end_particle: last particle to use
		Output
			PSrot_avg: a rotationally averaged 76-bin power spectrum for a set of images
	"""
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
	"""
		Name
			noise_corrected_PW - returns a noise corrected power spectrum and the factors a and b which were used to subtract the function f(x)=exp( a*x*x+b ) from the original power spectrum
		Input
			ps: a list containing the values of a power spectrum
			lo_limit: lower frequency threshold for minima search
			hi_limit: upper frequency threshold for minima search
			abs_limit: upper limit of the power spectrum, no usable information is contained above this threshold
		Output
			freq: list of corresponding frequency values
			pw_ns: noise corrected power spectrum
			b_factor: exponential factor
			norm: normalization factor
	"""
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
			Util.add_img2( var, Util.subn_img(img, avg) )

		reduce_EMData_to_root( var, myid, rootid )
		if myid == rootid:
			var /= (nimg-1)
			var.set_attr( "nimg", nimg )
			return var, avg
		else:
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
		else:
			return None

	def getvar(self):
		avg1 = self.sum1/self.nimg

		tmp = avg1.copy()
		Util.mul_img( tmp, avg1 )
		Util.sub_img(avg2 , tmp)

		avg2 *= (float(self.nimg)/float(self.nimg-1))
		 
		return avg2

	def getavg(self):
		return self.sum1/self.nimg


class inc_variancer:
	def __init__(self, nx, ny, nz):
		import numpy
		self.nx = nx
		self.ny = ny
		self.nz = nz
		self.ntot = nx*ny*nz
		self.sum1 = numpy.array( [0.0]*self.ntot )
		self.sum2 = numpy.array( [0.0]*self.ntot )
		self.nimg = 0

	def insert(self, img):
		from numpy import reshape
		from utilities import get_image_data
		data = reshape( get_image_data(img), (self.ntot,) )
		self.sum1 += data
		self.sum2 += (data*data)
		self.nimg += 1
		data = None


        def mpi_getvar(self, myid, rootid):
		from utilities import memory_usage, get_image_data, model_blank
		from mpi import mpi_reduce, MPI_DOUBLE, MPI_INT, MPI_SUM, MPI_COMM_WORLD
		from numpy import reshape

		cpy1 = self.sum1.copy()
		cpy2 = self.sum2.copy()
		sum1 = mpi_reduce( cpy1, self.ntot, MPI_DOUBLE, MPI_SUM, rootid, MPI_COMM_WORLD )
		sum2 = mpi_reduce( cpy2, self.ntot, MPI_DOUBLE, MPI_SUM, rootid, MPI_COMM_WORLD )
		sum_nimg = mpi_reduce( self.nimg, 1, MPI_INT, MPI_SUM, rootid, MPI_COMM_WORLD)
		if myid==rootid:

			nimg = int(sum_nimg[0])

			avg = model_blank( self.nx, self.ny, self.nz )
			var = model_blank( self.nx, self.ny, self.nz )
			vdata = reshape( get_image_data(var), (self.ntot,) )
			adata = reshape( get_image_data(avg), (self.ntot,) )
	
			
			for i in xrange(self.ntot):
				t1 = sum1[i]/nimg
				t2 = sum2[i]/nimg
				t2 = (t2 - t1*t1)*float(nimg)/float(nimg-1)
				adata[i] = t1
				vdata[i] = t2
			
			avg.set_attr( "nimg", nimg )
			var.set_attr( "nimg", nimg )

			del sum_nimg
			del sum1
			del sum2
			del cpy1
			del cpy2
			return var,avg

		del sum_nimg
		del sum1
		del sum2
		del cpy1
		del cpy2
		return model_blank(self.nx,self.ny,self.nz), model_blank(self.nx,self.ny,self.nz)


        def mpi_getavg(self, myid, rootid ):
		from mpi import mpi_reduce, MPI_INT, MPI_SUM, MPI_COMM_WORLD
		import numpy

		cpy1 = self.sum1.copy()
		cpy1 = mpi_reduce(cpy1, ntot, MPI_DOUBLE, MPI_SUM, rootid, MPI_COMM_WORLD)
		sum_nimg = mpi_reduce( self.nimg, 1, MPI_INT, MPI_SUM, rootid, MPI_COMM_WORLD)
		
		if myid==rootid:
			sum_nimg = int( sum_nimg[0] )
			cpy1 /= sum_nimg

			avg = model_blank( self.nx, self.ny, self.nz )
			adata = get_image_data(var).reshape( [self.ntot] ).astype( numpy.float32 )

			adata[0:self.ntot] = cpy1[0:self.ntot]
			avg.set_attr( "nimg", sum_nimg )
			return avg

		return None

	def getvar(self):
		avg1 = self.sum1/self.nimg
		avg2 = self.sum2/self.nimg

		tmp = avg1.copy()
		Util.mul_img( tmp, avg1 )
		Util.sub_img(avg2 , tmp)

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


class pcanalyzer:
	def __init__(self, mask, nvec=3, incore=False, MPI=False, scratch=None):
		import os
		self.mask = mask.copy()
		if MPI:
			from mpi import mpi_comm_rank, MPI_COMM_WORLD
			self.myid = mpi_comm_rank( MPI_COMM_WORLD )
			if( scratch == None):
				self.file = os.path.join("." , "maskedimg%04d.bin" % self.myid )
			else:
				self.file = os.path.join(scratch , "maskedimg%04d.bin" % self.myid )
			self.MPI  = True
		else:
			if( scratch == None):
				self.file = os.path.join("." , "maskedimg.bin" )
			else:
				self.file = os.path.join(scratch , "maskedimg.bin" )
			self.MPI  = False
			self.myid = 0
		#self.sdir   = "."
		self.nimg   = 0
		self.nvec   = nvec
		self.fw     = None
		self.fr     = None
		self.avgdat = None
		self.myBuff = []
		self.myBuffPos = 0
		self.incore = incore

	def writedat( self, data ):
		import array
		if not self.incore:
			if self.fw is None:
				self.fw = open( self.file, "wb" )
			data.tofile( self.fw )
		else:
			if len(self.myBuff) <= self.myBuffPos:
				self.myBuff.append(data.copy())
				self.myBuffPos = len(self.myBuff)
			else:
				self.myBuff[self.myBuffPos] = data.copy()
				self.myBuffPos += 1

	def read_dat( self, data ):
		from numpy import fromfile, float32
		if not self.incore:
			if not(self.fw is None) and not( self.fw.closed ):
				self.fw.close()
			if self.fr is None:
				self.fr = open( self.file, "rb" )
			assert not(self.fr is None) and not self.fr.closed
			Util.readarray( self.fr, data, self.ncov )
		else:
			data[:] = self.myBuff[self.myBuffPos]
			self.myBuffPos += 1
		if not(self.avgdat is None):
			data -= self.avgdat

	def close_dat( self ):
		if not self.incore:
			if not(self.fw is None) and not( self.fw.closed ):
				self.fw.close()
			self.fw = None
			if not(self.fr is None) and not( self.fr.closed ):
				self.fr.close()
			self.fr = None	
		else:
			self.myBuffPos = 0

	def usebuf( self ):
		nx = self.mask.get_xsize()
		ny = self.mask.get_ysize()
		nz = self.mask.get_zsize()

		self.ncov = 0
		for ix in xrange(nx):
			for iy in xrange(ny):
				for iz in xrange(nz):
					if( self.mask.get_value_at(ix,iy,iz) >= 0.5 ):
						self.ncov += 1
		import os
		size = os.stat( self.file )[6]
		self.nimg = size/(self.ncov*4)
		assert self.nimg * self.ncov*4 == size
		self.bufused = True

	def shuffle( self ):
		assert self.bufused
		from random import shuffle, seed
		from numpy  import zeros, float32, array
		from string import replace
		seed( 10000 + 10*self.myid )

		shfflfile = replace( self.file, "masked", "shuffled" )

		#print self.myid, "shuffling"
		sumdata = zeros( (self.ncov), float32 )
		imgdata = zeros( (self.ncov), float32 )
		if not self.incore: 
			self.fr = open( self.file, "rb" )
		self.avgdata = None

		if not self.incore: 
			fw = open( shfflfile, "wb" )
		for i in xrange(self.nimg):
			self.read_dat( imgdata )
			shuffle( imgdata )
			sumdata += imgdata
			if not self.incore:
				imgdata.tofile( fw )
			else:
				self.myBuff[self.myBuffPos-1] = imgdata.copy()

		if not self.incore: 
			self.fr.close()
			fw.close()
		else:
			self.close_dat()
		
		if self.MPI:
			from mpi import mpi_reduce, mpi_bcast, MPI_FLOAT, MPI_INT, MPI_SUM, MPI_COMM_WORLD
			sumdata = mpi_reduce( sumdata, self.ncov, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD )
			sumdata = mpi_bcast(  sumdata, self.ncov, MPI_FLOAT, 0, MPI_COMM_WORLD )
			sumdata = array(sumdata, float32)
 
			sumnimg = mpi_reduce( self.nimg, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD )
			sumnimg = mpi_bcast(  sumnimg,   1, MPI_INT, 0, MPI_COMM_WORLD )
		else:
			sumnimg = self.nimg

		self.file = shfflfile
		self.avgdat = sumdata[:]/float(sumnimg)
		#print "done shuffling,nimg:", float(sumnimg)



	def setavg( self, avg ):
		from numpy import zeros, float32
		from utilities import get_image_data
		tmpimg = Util.compress_image_mask( avg, self.mask )
		avgdat = get_image_data(tmpimg)
		self.avgdat = zeros( (len(avgdat)), float32 )
		self.avgdat[:] = avgdat[:]

	def insert( self, img ):
		assert self.mask.get_xsize()==img.get_xsize()
		assert self.mask.get_ysize()==img.get_ysize()
		assert self.mask.get_zsize()==img.get_zsize()

		from utilities import get_image_data
		tmpimg = Util.compress_image_mask( img, self.mask )
		tmpdat = get_image_data(tmpimg)
		if self.incore:
			self.myBuffPos = len(self.myBuff)
		self.writedat( tmpdat )
		if self.incore:
			self.close_dat()                                   #   WRITEDAT
		self.nimg +=1
		self.ncov = tmpimg.get_xsize()

	def analyze( self ):
		#if self.myid==0:
		#	print "analyze: ", self.ncov, " nvec: ", self.nvec
		from time import time
		from numpy import zeros, float32, int32, int64
		ncov = self.ncov
		kstep = self.nvec + 20 # the choice of kstep is purely heuristic

		diag    = zeros( (kstep), float32 )
		subdiag = zeros( (kstep), float32 )
		vmat    = zeros( (kstep, ncov), float32 )

		lanczos_start = time()
		kstep = self.lanczos( kstep, diag, subdiag, vmat )
		#print 'time for lanczos: ', time() - lanczos_start
		if not self.MPI or self.myid==0:
			qmat = zeros( (kstep,kstep), float32 )
			lfwrk = 100 + 4*kstep + kstep*kstep
			liwrk =   3 + 5*kstep

			fwork = zeros( (lfwrk), float32 )
			iwork = zeros( (liwrk), int32 )
			info = Util.sstevd( "V", kstep, diag, subdiag, qmat, kstep, fwork, lfwrk, iwork, liwrk)

			eigval = zeros( (self.nvec), float32 )
			for j in xrange(self.nvec):
				eigval[j] = diag[kstep-j-1]

			from utilities import model_blank, get_image_data
			eigimgs = []
			for j in xrange(self.nvec):
				tmpimg = model_blank(ncov, 1, 1)
				eigvec = get_image_data( tmpimg )
				trans = 'N'
				Util.sgemv( trans, ncov, kstep, 1.0, vmat, ncov, qmat[kstep-j-1], 1, 0.0, eigvec, 1 );

				eigimg = Util.reconstitute_image_mask(tmpimg, self.mask)
				eigimg.set_attr( "eigval", float(eigval[j])/(self.nimg - 1) )
				eigimgs.append( eigimg )

			return eigimgs


	def lanczos( self, kstep, diag, subdiag, V ):
		from numpy import zeros, float32, array
		from time import time

		all_start = time()

		ncov = self.ncov
		v0 = zeros( (ncov), float32)
		Av = zeros( (ncov), float32)

		hvec = zeros( (kstep), float32 )
		htmp = zeros( (kstep), float32 )
		imgdata = zeros( (ncov), float32 )

		for i in xrange(ncov):
			v0[i] = 1.0

		beta = Util.snrm2(ncov, v0, 1)
		for i in xrange(ncov):
			V[0][i] = v0[i]/beta

		for i in xrange(self.nimg):
			self.read_dat(imgdata)                                     #  READ_DAT			
			alpha = Util.sdot( ncov, imgdata, 1, V[0], 1 )
			Util.saxpy( ncov, alpha, imgdata, 1, Av, 1 )
		self.close_dat()

		if self.MPI:
			from mpi import mpi_reduce, mpi_bcast, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD
			Av = mpi_reduce( Av, ncov, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD )
			Av = mpi_bcast(  Av, ncov, MPI_FLOAT, 0, MPI_COMM_WORLD )
			Av = array(Av, float32)
		#print 'iter 0: ', time() - all_start


		diag[0] = Util.sdot( ncov, V[0], 1, Av, 1 )
		alpha = -diag[0]
		Util.saxpy( ncov, float(alpha), V[0], 1, Av, 1 )

		TOL = 1.0e-7
		for iter in xrange(1, kstep):
			iter_start = time()
			beta = Util.snrm2(ncov, Av, 1)
			if( beta < TOL ):
				kstep = iter+1
				break

			subdiag[iter-1] = beta
			for i in xrange(ncov):
				V[iter][i] = Av[i]/beta

			Av[:] = 0.0

			for i in xrange(self.nimg):
				self.read_dat( imgdata )                                #READ_DAT
				alpha = Util.sdot( ncov, imgdata, 1, V[iter], 1 )
				Util.saxpy( ncov, float(alpha), imgdata, 1, Av, 1 )
			self.close_dat()


			if self.MPI:
				Av = mpi_reduce( Av, ncov, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD )
				Av = mpi_bcast(  Av, ncov, MPI_FLOAT, 0, MPI_COMM_WORLD )
				Av = array(Av, float32)

			trans = 'T'
			Util.sgemv( trans, ncov, iter+1,  1.0, V, ncov, Av, 1,
			              0.0, hvec, 1 )

			trans = 'N'
			Util.sgemv( trans, ncov, iter+1, -1.0, V, ncov, hvec, 1,
			              1.0,     Av, 1 )

			trans = 'T'
			Util.sgemv( trans, ncov, iter+1,  1.0, V, ncov, Av, 1,
			              0.0,   htmp, 1 )

			Util.saxpy(iter+1, 1.0, htmp, 1, hvec, 1)

			trans = 'N'
			Util.sgemv( trans, ncov, iter+1, -1.0, V, ncov, htmp, 1,
			              1.0,     Av, 1 )

			diag[iter] = hvec[iter]

			#print 'iter, time, overall_time: ', iter, time()-iter_start, time()-all_start
		return kstep



class pcanalyzebck:
	def __init__(self, mask, nvec, dataw, list_of_particles, dm, variance, fl, aa, MPI=False ):
		import os
		self.mask = mask.copy()
		if MPI:
			from mpi import mpi_comm_rank, MPI_COMM_WORLD
			self.myid = mpi_comm_rank( MPI_COMM_WORLD )
			"""
			if( scratch == None):
				self.file = os.path.join(sdir , "maskedimg%04d.bin" % self.myid )
			else:
				self.file = os.path.join(scratch , "maskedimg%04d.bin" % self.myid )
			"""
			self.MPI  = True
		else:
			"""
			if( scratch == None):
				self.file = os.path.join(sdir , "maskedimg.bin" )
			else:
				self.file = os.path.join(scratch , "maskedimg.bin" )
			"""
			self.MPI  = False
			self.myid = 0
		#self.sdir   = sdir
		self.dataw  = dataw
		self.list_of_particles = list_of_particles
		self.dm     = dm
		self.variance = variance
		self.nimg   = len(dataw)
		self.nvec   = nvec
		self.fw     = None
		self.fr     = None
		self.fl     = fl
		self.aa     = aa
		self.avgdat = None
		self.getncov()

	"""
	def writedat( self, data ):
		import array

		if self.fw is None:
			self.fw = open( self.file, "wb" )

		data.tofile( self.fw )

	def read_dat( self, data ):
		from numpy import fromfile, float32
		if not(self.fw is None) and not( self.fw.closed ):
			self.fw.close()

		assert not(self.fr is None) and not self.fr.closed
		Util.readarray( self.fr, data, self.ncov )
		if not(self.avgdat is None):
			data -= self.avgdat
	"""
	def get_dat( self, k ):
		from reconstruction import backproject_swbp
		from filter  import filt_tanl
		vb = Util.divn_img(backproject_swbp(self.dataw[k], self.list_of_particles[k], self.dm), self.variance)
		if(self.fl > 0.0):
			vb = filt_tanl(vb, self.fl, self.aa)
		#vb -= pc[0]
		#vb *= (refstat[1]/pc[1])
		from numpy import zeros, float32
		from utilities import get_image_data

		tmpimg = Util.compress_image_mask( vb, self.mask )
		data = get_image_data(tmpimg)
		if not(self.avgdat is None):
			data -= self.avgdat
		return data

	def getncov( self ):  # used to be called usebuf
		nx = self.mask.get_xsize()
		ny = self.mask.get_ysize()
		nz = self.mask.get_zsize()

		self.ncov = 0
		for ix in xrange(nx):
			for iy in xrange(ny):
				for iz in xrange(nz):
					if( self.mask.get_value_at(ix,iy,iz) >= 0.5 ):
						self.ncov += 1
		#import os
		#size = os.stat( self.file )[6]
		#self.nimg = size/(self.ncov*4)
		#assert self.nimg * self.ncov*4 == size
		#self.bufused = True

	"""
	def shuffle( self ):
		assert self.bufused
		from random import shuffle, seed
		from numpy  import zeros, float32, array
		from string import replace
		seed( 10000 + 10*self.myid )

		shfflfile = replace( self.file, "masked", "shuffled" )

		#print self.myid, "shuffling"
		sumdata = zeros( (self.ncov), float32 )
		imgdata = zeros( (self.ncov), float32 )
		self.fr = open( self.file, "rb" )
		self.avgdata = None

		fw = open( shfflfile, "wb" )
		for i in xrange(self.nimg):
			self.read_dat( imgdata )
			shuffle( imgdata )
			sumdata += imgdata
			imgdata.tofile( fw )

		self.fr.close()
		fw.close()

		if self.MPI:
			from mpi import mpi_reduce, mpi_bcast, MPI_FLOAT, MPI_INT, MPI_SUM, MPI_COMM_WORLD
			sumdata = mpi_reduce( sumdata, self.ncov, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD )
			sumdata = mpi_bcast(  sumdata, self.ncov, MPI_FLOAT, 0, MPI_COMM_WORLD )
			sumdata = array(sumdata, float32)
 
			sumnimg = mpi_reduce( self.nimg, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD )
			sumnimg = mpi_bcast(  sumnimg,   1, MPI_INT, 0, MPI_COMM_WORLD )

		self.file = shfflfile
		self.avgdat = sumdata[:]/float(sumnimg)
		#print "done shuffling,nimg:", float(sumnimg)
	"""

	def setavg( self, avg ):
		from numpy import zeros, float32
		from utilities import get_image_data
		tmpimg = Util.compress_image_mask( avg, self.mask )
		avgdat = get_image_data(tmpimg)
		self.avgdat = zeros( (len(avgdat)), float32 )
		self.avgdat[:] = avgdat[:]

	"""
	def insert( self, img ):
		assert self.mask.get_xsize()==img.get_xsize()
		assert self.mask.get_ysize()==img.get_ysize()
		assert self.mask.get_zsize()==img.get_zsize()

		from utilities import get_image_data
		tmpimg = Util.compress_image_mask( img, self.mask )
		tmpdat = get_image_data(tmpimg)
		self.writedat( tmpdat )                                   #   WRITEDAT
		self.nimg +=1
		self.ncov = tmpimg.get_xsize()
	"""
	def analyze( self ):
		if self.myid==0:
			print "analyze: ", self.ncov, " nvec: ", self.nvec
		from time import time
		from numpy import zeros, float32, int32, int64
		ncov = self.ncov
		kstep = self.nvec + 20 # the choice of kstep is purely heuristic

		diag    = zeros( (kstep), float32 )
		subdiag = zeros( (kstep), float32 )
		vmat    = zeros( (kstep, ncov), float32 )

		lanczos_start = time()
		kstep = self.lanczos( kstep, diag, subdiag, vmat )
		print 'time for lanczos: ', time() - lanczos_start

		if not self.MPI or self.myid==0:
			qmat = zeros( (kstep,kstep), float32 )
			lfwrk = 100 + 4*kstep + kstep*kstep
			liwrk =   3 + 5*kstep

			fwork = zeros( (lfwrk), float32 )
			iwork = zeros( (liwrk), int32 )
			info = Util.sstevd( "V", kstep, diag, subdiag, qmat, kstep, fwork, lfwrk, iwork, liwrk)

			eigval = zeros( (self.nvec), float32 )
			for j in xrange(self.nvec):
				eigval[j] = diag[kstep-j-1]

			from utilities import model_blank, get_image_data
			eigimgs = []
			for j in xrange(self.nvec):
				tmpimg = model_blank(ncov, 1, 1)
				eigvec = get_image_data( tmpimg )
				trans = 'N'
				Util.sgemv( trans, ncov, kstep, 1.0, vmat, ncov, qmat[kstep-j-1], 1, 0.0, eigvec, 1 );

				eigimg = Util.reconstitute_image_mask(tmpimg, self.mask)
				eigimg.set_attr( "eigval", float(eigval[j]) )
				eigimgs.append( eigimg )

			return eigimgs


	def lanczos( self, kstep, diag, subdiag, V ):
		from numpy import zeros, float32, array
		from time import time
		all_start = time()
	
		ncov = self.ncov
		v0 = zeros( (ncov), float32)
		Av = zeros( (ncov), float32)

		hvec = zeros( (kstep), float32 )
		htmp = zeros( (kstep), float32 )
		imgdata = zeros( (ncov), float32 )
	
		for i in xrange(ncov):
			v0[i] = 1.0

		beta = Util.snrm2(ncov, v0, 1)
		for i in xrange(ncov):
			V[0][i] = v0[i]/beta

		#self.fr = open( self.file, "rb" )
		for i in xrange(self.nimg):
			#self.read_dat(imgdata)                                     #  READ_DAT
			imgdata = self.get_dat(i)
			alpha = Util.sdot( ncov, imgdata, 1, V[0], 1 )
			Util.saxpy( ncov, alpha, imgdata, 1, Av, 1 )
		#self.fr.close()

		if self.MPI:
			from mpi import mpi_reduce, mpi_bcast, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD
			Av = mpi_reduce( Av, ncov, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD )
			Av = mpi_bcast(  Av, ncov, MPI_FLOAT, 0, MPI_COMM_WORLD )
			Av = array(Av, float32)
		#print 'iter 0: ', time() - all_start


		diag[0] = Util.sdot( ncov, V[0], 1, Av, 1 )
		alpha = -diag[0]
		Util.saxpy( ncov, float(alpha), V[0], 1, Av, 1 )

		TOL = 1.0e-7
		for iter in xrange(1, kstep):
			iter_start = time()
			beta = Util.snrm2(ncov, Av, 1)
			if( beta < TOL ):
				kstep = iter+1
				break

			subdiag[iter-1] = beta
			for i in xrange(ncov):
				V[iter][i] = Av[i]/beta

			Av[:] = 0.0

			#self.fr = open( self.file, "rb" )
			for i in xrange(self.nimg):
				#self.read_dat( imgdata )                                #READ_DAT
				imgdata = self.get_dat( i)
				alpha = Util.sdot( ncov, imgdata, 1, V[iter], 1 )
				Util.saxpy( ncov, float(alpha), imgdata, 1, Av, 1 )
			#self.fr.close()


			if self.MPI:
				Av = mpi_reduce( Av, ncov, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD )
				Av = mpi_bcast(  Av, ncov, MPI_FLOAT, 0, MPI_COMM_WORLD )
				Av = array(Av, float32)

			trans = 'T'
			Util.sgemv( trans, ncov, iter+1,  1.0, V, ncov, Av, 1,
			              0.0, hvec, 1 )

			trans = 'N'
			Util.sgemv( trans, ncov, iter+1, -1.0, V, ncov, hvec, 1,
			              1.0,     Av, 1 )

			trans = 'T'
			Util.sgemv( trans, ncov, iter+1,  1.0, V, ncov, Av, 1,
			              0.0,   htmp, 1 )

			Util.saxpy(iter+1, 1.0, htmp, 1, hvec, 1)

			trans = 'N'
			Util.sgemv( trans, ncov, iter+1, -1.0, V, ncov, htmp, 1,
			              1.0,     Av, 1 )

			diag[iter] = hvec[iter]

			print 'iter, time, overall_time: ', iter, time()-iter_start, time()-all_start
		return kstep



def k_means_stab_bbenum(PART, T=10, nguesses=5, J=50, max_branching=40, stmult=0.25, branchfunc=2, LIM=-1, doMPI_init=False, njobs=-1,do_mpi=False, K=-1,cdim=[],nstart=-1, nstop=-1, top_Matches=[]):
	from statistics import k_means_match_clusters_asg_new
	"""
		
		Input:
			PART:   list of list of arrays. So PART[0] corresponds to the first partition, PART[1] the second partition and so on. 
				Each partition is a list of arrays of distinct integers sorted in ascending order. The arrays corresponds
				to classes, so PART[i][j] is the j-th class of the i-th partition. 
				
				Example of how to construct a partition:
					K    = EMUtil.get_image_count(stack_name)
					part0 = []
					for k in xrange(K):
						im  = get_im(stack_name, k)
						lid = im.get_attr('members')
						lid = array(lid, 'int32')
						lid.sort()
						part0.append(lid.copy())
						
			K:	Number of classes in the first partition.
	
			T:      An integer. It is the threshold such that the stable sets corresponding to the matches in the output have size larger than T. 
				Specifically, if there are say four partitions, and the i-th match in the output, i.e., MATCH[i], is [2,12,1,5], then
				the 2nd class of the first partition, 12th class of the second, first class of the third, and fifth class of the
				fourth have at least T elements in common.
			
			nguesses: Not used anymore. I'm going to remove it.
			
			J:	An integer. In the branching matching algorithm we use, each step corresponds to choosing a match to add to the 
				collection of matches we will eventually output. 
				Since at each step there are many different possibilities and 
				it is generally not feasible to branch on them all, we consider only J matches with the largest weights in choosing which matches to branch on.
				See branchfunc below for more details.
				Intuitively, a larger J should give better results. 
				If J=1, the algorithm defaults to the greedy algorithm, i.e, at each step, we choose the match which has the largest cost and add it to the collection of matches to output.
			        
			
			branchfunc: An integer. Determines which branching function should be used. 
			
				    Roughly speaking, the algorithm we use builds up the collection of matches iteratively. During each step, we determine
				    the next match that should be added to the collection. Since there are many possibilities, the role of the branching function
				    is to determine both how many and which possibilities to explore.
				    
				    The branching functions chooses the possibilities to explore from the 
				    J matches with the largest weights.
				    
				    There are two possible choices for branchfunc: 2 or 4
				    
				    branchfunc = 2:  The J matches are considered in order according to weight, beginning with the match with the largest weight. 
				    	             Always branches on the current match with the largest weight. 
						     Branch on the match with the second largest weight only if it is infeasible with the match with the largest weight.
				    	             Similarly, for each subsequent match, only branch on it if it is infeasible with at least LIM of the matches which have already been chosen to branch on.
						    
						    
				    branchfunc = 4:  The J matches are considered in order according to weight, beginning with the match with the largest weight. 
				    		     We branch based on distribution of the cost of the J largest matches. 
				    		     Let stdev be the standard deviation of the weights of the J matches. 
						     We always branch on the match with the largest weight.
						     For the remaining J-1 matches, we branch on those which are within stmult*stdev of the weight of match with the largest weight.
						     	 	
				    	
			LIM: 	An integer smaller than J. See explanation for branchfunc above.
				  	
			max_branching: This is an upper bound on the maximum number of branches to explore. See explanation for branchfunc above. This is to ensure we get a result in reasonable time. Intuitively, the larger max_branching is, the likelihood of getting a better result (i.e, a collection of matches with a large total weight) is increased.  
	
			stmult: An integer. See explanation for branchfunc above.
	
		Output: MATCH, STB_PART, CT_s, CT_t, ST, st
			
			If there are exactly two partitions, i.e., len(PART) == 2, then the matching will be computed using Hungarian algorithm, and those matches for which the corresponding stable set has size greater than T will be output. If T==0, then the output matches will be optimal.
			
			(EDIT 10/11/11: If there are exactly two partitions, i.e., len(PART) == 2, then the matching will be computed using Hungarian algorithm such that ONLY matches with weights greater than T are considered during the course of the algorithm. This means the matching computed using the Hungarian algorithm will contain only matches with weight greater than T. See k_means_match_clusters_asg_new.) 
			
			MATCH: A list of arrays. Each array has len(PART) elements, and the i-th element corresponds to a class in the i-th partition. 
			       So MATCH[0][0] is an index into PART[0], MATCH[0][1] is an index into PART[1], etc.
			    
			STB_PART: A list of lists. The stable set corresponding to the i-th MATCH, i.e., MATCH[i], is stored in STB_PART[MATCH[i][0]]. 
			
			CT_s:   A list of integers. The size of the stable set corresponding to  the i-th MATCH, i.e., MATCH[i], is stored in CT_s[MATCH[i][0]].
				The quality of the output collection of matches, i.e., MATCH, can be determined by summing CT_s. The larger the better.
			
			CT_t:	A list of integers. The number of elements in the union of the classes in the i-th match is stored in CT_t[MATCH[i][0]].
			
			st:     st = 100.0 * sum(CT_s) / float(sum(CT_t)) 
			
			ST:	ST[k] = 100.0 * CT_s[k] / float(CT_t[k])	
	
	"""
	
	from copy import deepcopy
	from numpy import array, append

	MATCH=[]
	np = len(PART)
	# do MPI_init: compute pruned partitions and Njobs top matches and return
	#if doMPI_init:
	#	newParts,topMatches,class_dim, num_found=k_means_match_bbenum(PART,T=T, J=J, nguesses=nguesses, DoMPI_init=True,Njobs=njobs)
	#	return newParts, topMatches, class_dim,num_found
	
	#if do_mpi:
	#	MATCH,cost= k_means_match_bbenum(PART, T=T, J=J, nguesses=nguesses, DoMPI=True, K=K, np=np,c_dim=cdim, N_start=nstart,N_stop=nstop,topMatches=top_Matches)
	#	return MATCH,cost
			
	if np > 2:
		MATCH= k_means_match_bbenum(PART, T=T, J=J, max_branching=max_branching, stmult=stmult, branchfunc=branchfunc, LIM=LIM, nguesses=nguesses, DoMPI=False)

	
	list_stb=[]
	tot_n=0
	if np == 2:
		# First make sure the two partitions have the same number of classes. If not, pad the one with less with junk.
		K1 = len(PART[0])
		K2 = len(PART[1])
		maxK = max(K1, K2)
		if K1 < maxK or K2 < maxK:
			# ffind garbage value to pad empty classes with
			garbage_value = 923456
			garbage_incr=1
			for i in xrange(np):
				K = len(PART[i])
				for j in xrange(K):
					pSize = PART[i][j].size
					if pSize >= 1:
						for p in xrange(pSize):
							if PART[i][j][p] >= garbage_value:
								garbage_value = PART[i][j][p] + garbage_incr
			for i in xrange(np):
				if len(PART[i]) < maxK:
					# pad with empty arrays
					df = maxK - len(PART[i])
					for pd in xrange(df):
						PART[i].append(array([garbage_value, garbage_value+1],'int32'))
						garbage_value = garbage_value + 2
						
		# now call 
		MATCH, list_stb, tot_n = k_means_match_clusters_asg_new(PART[0], PART[1], T)
			
							
	K=len(PART[0])		
	
	STB_PART = [[] for i in xrange(K)]
	nm       = len(MATCH)
	CT_t     = [0] * K
	CT_s     = [0] * K
	ST       = [0] * K
	ct_t     = 0
	ct_s     = 0
	st       = 0


	for k in xrange(nm):
		kk   = int(MATCH[k][0]) # due to numpy obj
		vmax = [0] * np
		vmin = [0] * np
		for i in xrange(np):
		    vmax[i] = max(PART[i][int(MATCH[k][i])])
		    vmin[i] = min(PART[i][int(MATCH[k][i])])

		vmax = int(max(vmax))
		vmin = int(min(vmin))
		vd   = vmax - vmin + 1

		asg = [0] * vd
		for i in xrange(np):
		    for item in PART[i][int(MATCH[k][i])]: asg[int(item) - vmin] += 1

		stb  = []
		for i in xrange(vd):
		    if asg[i] != 0:
			CT_t[kk] += 1
			if asg[i] == np:
			    CT_s[kk] += 1
			    stb.append(i + vmin)

		STB_PART[kk] = deepcopy(stb)

	for k in xrange(K):
		if CT_t[k] == 0: continue
		ST[k] = 100.0 * CT_s[k] / float(CT_t[k])

	if sum(CT_t) == 0:
		st = 0
	else:   st = 100.0 * sum(CT_s) / float(sum(CT_t))

	if np == 2:
		if sum(CT_s) != tot_n:
			print sum(CT_s), tot_n
			print "something wrong!!"
			sys.exit()

	return MATCH, STB_PART, CT_s, CT_t, ST, st

# DO NOT copy memory - could lead to crashes	
# This is the wrapper function for bb_enumerateMPI. It packages the arguments and formats the output....
def k_means_match_bbenum(PART, T=10, J=1, max_branching=40, stmult=0.25, nguesses=5, branchfunc=2, LIM=-1, DoMPI_init=False, Njobs=-1, DoMPI=False, K=-1, np=-1, c_dim=[],N_start=-1, N_stop=-1, topMatches=[]):	
	from numpy import array, append, insert,sort
	
	MATCH=[]
	output=[]
	
	
	if DoMPI==True and DoMPI_init==False:
		print "Not supporting MPI currently"
		#if len(levels)<1:	
		#	for i in xrange(K):
		#		levels.append(1)

		#ar_levels = array(levels, 'int32')
		#ar_class_dim	= array(c_dim,'int32')	
		#ar_newParts = array(PART, 'int32')
		#firstmatches = topMatches[N_start*(np+1):(N_stop+1)*(np+1)]
		#output=Util.branchMPIpy(ar_newParts,ar_class_dim,np,K,T,ar_levels,K,nguesses,N_stop-N_start+1, array(firstmatches, 'int32'))
	else:
		LARGEST_CLASS = 0
		np=len(PART)
		
		#ar_argParts = array([],'int32')
		onedParts = []
		class_dim=[]
		
		# figure out the garbage value to pad empty classes with
		
		garbage_value = 923456
		garbage_incr = 10
		for i in xrange(np):
			K = len(PART[i])
			for j in xrange(K):
				pSize = PART[i][j].size
				if pSize > LARGEST_CLASS:
					LARGEST_CLASS = pSize
				if pSize >= 1:
					for p in xrange(pSize):
						if PART[i][j][p] >= garbage_value:
							garbage_value = PART[i][j][p] + garbage_incr
		
		# deal with the case where not all the partitions have the same number of classes
		max_K = len(PART[0])
		for i in xrange(np):
			if len(PART[i]) > max_K:
				max_K = len(PART[i])
		for i in xrange(np):
			if not(len(PART[i]) == max_K):
				# pad with empty arrays
				df = max_K - len(PART[i])
				for j in xrange(df):
					PART[i].append(array([garbage_value, garbage_value+1],'int32'))
					garbage_value = garbage_value + 2	
					
		
		K = len(PART[0])		
		# serialize all the arguments in preparation for passing into c++ function
		
		
		for i in xrange(np):
			for j in xrange(K):
				pSize = PART[i][j].size
					
				onedParts.append(j)
				onedParts.append(0)
				zero_pad = 0
				if pSize > 0:
					for p in xrange(pSize):
						onedParts.append(PART[i][j][p])	
				else:
					# pad with garbage value
					zero_pad=2
					onedParts.append(garbage_value)
					onedParts.append(garbage_value+1)
					garbage_value = garbage_value + zero_pad
						
				class_dim.append(pSize + zero_pad + 2)
				
				#ar_argParts = append(ar_argParts,[j,0])
				#ar_argParts=append(ar_argParts,PART[i][j])
		ar_argParts = array(onedParts,'int32')
		ar_class_dim = array(class_dim,'int32')
		

		# Single processor version
		output = Util.bb_enumerateMPI(ar_argParts, ar_class_dim,np,K,T, nguesses,LARGEST_CLASS, J, max_branching, stmult, branchfunc, LIM)
		
	# first element of output is the total  cost of the solution, second element is the number of matches
	# in the output solution, and then follows the list of matches.
	
	# convert each match into an array and insert into MATCH
	num_matches = output[1]
	
	for j in xrange(num_matches):
		# get the j-th match
		ar_match = array(output[j*np + 2: j*np + 2+np],'int32')
		MATCH.append(ar_match)
		
	# order Matches in Match by group in first partition
	outMATCH=[]
	for j in xrange(num_matches):
		amatch=[]
		for js in xrange(len(MATCH[j])):
			amatch.append(MATCH[j][js])
		outMATCH.append(amatch)
	outMATCH.sort()
			
	if DoMPI:
		print "Not supporting MPI currently"
		return outMATCH, output[0]
	else:
		return outMATCH
# match is a list, where every five tuple corresponds to a match
def k_means_stab_getinfo(PART, match):
	
	from copy import deepcopy
	from numpy import array
	K=len(PART[0])
	np = len(PART)
	
	MATCH=[]
	
	# convert argument match to a list of arrays, where each array is a match
	len_match = len(match)
	if (len_match % np) != 0:
		print "something wrong in k_means_stab_getinfo"
		sys.exit()
	num_matches = len_match/np
	for i in xrange(num_matches):
		MATCH.append(array(match[i*np:(i+1)*np]))
	
	STB_PART = [[] for i in xrange(K)]
	nm       = len(MATCH)
	CT_t     = [0] * K
	CT_s     = [0] * K
	ST       = [0] * K
	ct_t     = 0
	ct_s     = 0
	st       = 0


	for k in xrange(nm):
		kk   = int(MATCH[k][0]) # due to numpy obj
		vmax = [0] * np
		vmin = [0] * np
		for i in xrange(np):
		    vmax[i] = max(PART[i][int(MATCH[k][i])])
		    vmin[i] = min(PART[i][int(MATCH[k][i])])

		vmax = int(max(vmax))
		vmin = int(min(vmin))
		vd   = vmax - vmin + 1

		asg = [0] * vd
		for i in xrange(np):
		    for item in PART[i][int(MATCH[k][i])]: asg[int(item) - vmin] += 1

		stb  = []
		for i in xrange(vd):
		    if asg[i] != 0:
			CT_t[kk] += 1
			if asg[i] == np:
			    CT_s[kk] += 1
			    stb.append(i + vmin)

		STB_PART[kk] = deepcopy(stb)

	for k in xrange(K):
		if CT_t[k] == 0: continue
		ST[k] = 100.0 * CT_s[k] / float(CT_t[k])

	if sum(CT_t) == 0:
		st = 0
	else:   st = 100.0 * sum(CT_s) / float(sum(CT_t))

	return MATCH, STB_PART, CT_s, CT_t, ST, st

def match_lists(l1, l2):
	count = 0
	for l in l1:
		try:
			i = l2.index(l)
			count += 1
		except:
			pass
	return count

def center_of_gravity(a):
	return a.cog()


def center_of_gravity_phase(a):
	return a.phase_cog()

def fit_ctf(crossresolution, ctf_params, rangedef = -1.0, i1 = 0, i2 = 0, chisquare=False):
	"""
		ctf_params = [defocus, cs, voltage, apix, bfactor, ampcont]
	"""
	from math import copysign
	from morphology import ctf_1d
	from utilities import generate_ctf
	n = len(crossresolution[1])
	if(i2 <= i1):  i2 = n
	nx = 2*n
	if(rangedef == -1.0): rangedef = ctf_params[0]*0.1
	sgncrs = [0.0]*n
	for i in xrange(n):  sgncrs[i] = copysign(1.0, crossresolution[1][i])
	
	qt = 1.0e23
	nstep = 21
	for j in xrange(21):
		defi = ctf_params[0]-rangedef + rangedef*0.1*j
		ctf = ctf_1d(nx, generate_ctf([defi]+ctf_params[1:]))
		disc = 0.0
		if chisquare:
			for k in xrange(i1,i2):
				disc += ((sgncrs[k] - copysign(1.0, ctf[k]))/crossresolution[2][k])**2
		else:
			for k in xrange(i1,i2):
				disc += (sgncrs[k] - copysign(1.0, ctf[k]))**2
		if( disc < qt):
			best_def = defi
			qt = disc
	return best_def

def randprojdir(ang, sigma):
	""" 
		Randomize projection angles
		ang - projection directions in rows
		output - same as ang, but with first two positions (phi, theta) dispersed using gaussian noise with sigma
	"""
	import random as rdq
	l = len(ang[0])
	aout = []
	for n in xrange(len(ang)):
		t = Transform({"type":"spider","phi":ang[n][0],"theta":ang[n][1],"psi":0.0})
		r = Transform({"type":"spider","phi":rdq.random()*360.0,"theta":abs(rdq.gauss(0.0,sigma)),"psi":0.0})
		r = r*t
		d = r.get_params("spider")
		aout.append([d["phi"],d["theta"]] + [ang[n][k] for k in xrange(2,l)])

	return aout
