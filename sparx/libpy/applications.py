#
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holfds
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

def ali2d(stack, outdir, maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", dst=0.0, center=-1, maxit=0, \
		CTF=False, snr=1.0, Fourvar=False, Ng=-1, user_func_name="ref_ali2d", CUDA=False, GPUID="", MPI=False):
	"""
		Name
			ali2d - Perform 2-D reference-free alignment of an image series
		Input
			stack: set of 2-D images in a stack file (format hdf or bdb), images have to be square (nx=ny)
			maskfile: optional mask file to be used internally during alignment
			inner_radius: inner radius for rotational correlation > 0 (set to 1)
			outer_radius: outer radius for rotational correlation < nx/2-1 (set to nx/2-2, should be set to the radius of the particle)
			ring_step: step between rings in rotational correlation > 0 (set to 1)
			x_range: range for translation search in x direction, search is +/xr 
			y_range: range for translation search in y direction, search is +/yr
			translation_step: step of translation search in both directions, 
			center_type
			max_iteration: maximum number of iterations the program will perform 
			CTF: if this flag is set, the program will use CTF information provided in file headers (for details see I_O). (default no CTF)
			snr: signal-to-noise ratio of the data (default SNR=1.0)
			Fourvar: use Fourier variance to weight the reference (recommended, default False)
			function: name of the user-supplied-function that prepares reference image for each iteration (default ref_ali2d)
			rand_alpha: if this flag is set, the program start with random alpha angle (default no rand_alpha)
			MPI: if this flag is set, use MPI version
		Output
			output_directory: directory name into which the output files will be written.
			header: the alignment parameters are stored in the headers of input files as 'xform.align2d'.
	"""
	if MPI:
		ali2d_MPI(stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, dst, center, maxit, CTF, snr, Fourvar, Ng, user_func_name, CUDA, GPUID)
		return

	from utilities    import print_begin_msg, print_end_msg, print_msg
	from utilities    import file_type
	import os

	print_begin_msg("ali2d")
	
	# Comment by Zhengfan Yang on 09/03/10
	# I have decided that outdir should be able to take None as parameter, if the user does not need an output directory.
	# This is particularly useful in the ISAC program, nobody will care what is in the alignment and realignment directory.
	# It takes a lot of disk space and time to write them, too. It may not be much for one iteration, but will be enormous 
	# for several hundred iterations. This change will not affect the normal use.
	if outdir:
		if os.path.exists(outdir):   ERROR('Output directory exists, please change the name and restart the program', "ali2d", 1)
		os.mkdir(outdir)

	if file_type(stack) == "bdb":
		from EMAN2db import db_open_dict
		dummy = db_open_dict(stack, True)
	active = EMUtil.get_all_attributes(stack, 'active')
	list_of_particles = []
	for im in xrange(len(active)):
		if active[im]:  list_of_particles.append(im)
	del active
	nima = len(list_of_particles)
	
	data = EMData.read_images(stack, list_of_particles)
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
	
	print_msg("Input stack                 : %s\n"%(stack))
	
	ali2d_data(data, outdir, maskfile, ir, ou, rs, xr, yr, ts, dst, center, maxit, CTF, snr, Fourvar, Ng, user_func_name, CUDA, GPUID, True)
	
	# write out headers
	from utilities import write_headers
	write_headers(stack, data, list_of_particles)
	print_end_msg("ali2d")


def ali2d_data(data, outdir, maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", dst=0.0, center=-1, maxit=0, \
		CTF=False, snr=1.0, Fourvar=False, Ng=-1, user_func_name="ref_ali2d", CUDA=False, GPUID="", from_ali2d=False):

	# Comment by Zhengfan Yang 02/25/11
	# This is where ali2d() actually runs, the reason I divided it into two parts is that
	# this alignment program supports list of EMData not a stack.

	from utilities    import drop_image, get_image, get_input_from_string, get_params2D, set_params2D
	from statistics   import fsc_mask, sum_oe, hist_list
	from alignment    import Numrinit, ringwe, ali2d_single_iter
	from pixel_error  import max_2D_pixel_error
	from fundamentals import fshift, fft, rot_avg_table
	from utilities    import print_begin_msg, print_end_msg, print_msg
	from utilities    import model_blank, model_circle
	import os

	if from_ali2d == False: print_begin_msg("ali2d_data")
	
	if outdir and from_ali2d == False:
		if os.path.exists(outdir):   ERROR('Output directory exists, please change the name and restart the program', "ali2d", 1)
		os.mkdir(outdir)

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)

	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit);

	if max_iter == 0:
		max_iter  = 10
		auto_stop = True
	else:
		auto_stop = False

	nima = len(data)

	if Ng == -1:
		Ng = nima
	elif Ng == -2:
		Ng = int(0.98*nima)

	nx = data[0].get_xsize()

	# default value for the last ring
	if last_ring == -1:  last_ring = nx/2-2

	import user_functions
	user_func = user_functions.factory[user_func_name]
	print_msg("Number of active images     : %s\n"%(nima))
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Inner radius                : %i\n"%(first_ring))
	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Ring step                   : %i\n"%(rstep))
	print_msg("X search range              : %s\n"%(xrng))
	print_msg("Y search range              : %s\n"%(yrng))
	print_msg("Translational step          : %s\n"%(step))
	print_msg("Discrete angle used         : %d\n"%(dst))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("Use Fourier variance        : %s\n"%(Fourvar))
	print_msg("Number of groups            : %d\n"%(Ng))
	print_msg("CTF correction              : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	if auto_stop:
		print_msg("Stop iteration with         : criterion\n")
	else:
		print_msg("Stop iteration with         : maxit\n")
	print_msg("User function               : %s\n"%(user_func_name))
	print_msg("Using CUDA                  : %s\n"%(CUDA))
	if CUDA:
		GPUID = get_input_from_string(GPUID)
		GPUID = int(GPUID[0])
		R = CUDA_Aligner(GPUID)
		print_msg("GPU ID                      : %d\n"%(GPUID))

	if maskfile:
		import	types
		if type(maskfile) is types.StringType:
			print_msg("Maskfile                    : %s\n\n"%(maskfile))
			mask = get_image(maskfile)
		else:
			print_msg("Maskfile                    : user provided in-core mask\n\n")
			mask = maskfile
	else :
		print_msg("Maskfile                    : default, a circle with radius %i\n\n"%(last_ring))
		mask = model_circle(last_ring, nx, nx)

	cnx = nx/2+1
 	cny = cnx
 	mode = "F"
	if CTF:
		if data[0].get_attr_default('ctf_applied', 2) > 0:	ERROR("data cannot be ctf-applied", "ali2d", 1)
		from filter import filt_ctf
		from morphology import ctf_img
		ctf_abs_sum = EMData(nx, nx, 1, False)
		ctf_2_sum = EMData(nx, nx, 1, False)
	else:
		ctf_2_sum = None

	if Fourvar:
		if CTF:
			from statistics   import varf2d
		else:
			from statistics   import varf

	if CUDA:
		all_ali_params = []
		all_ctf_params = []

	for im in xrange(nima):
		#  Subtract averages outside of mask from all input data
		st = Util.infomask(data[im], mask, False)
		data[im] -= st[0]
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			ctfimg = ctf_img(nx, ctf_params)
			if CUDA:
				all_ctf_params.extend([ctf_params.defocus, ctf_params.cs, ctf_params.voltage, ctf_params.apix, ctf_params.bfactor, ctf_params.ampcont])
			Util.add_img2(ctf_2_sum, ctfimg)
			Util.add_img_abs(ctf_abs_sum, ctfimg)
		if CUDA:
			alpha, sx, sy, mirror, scale = get_params2D(data[im])
			all_ali_params.extend([alpha, sx, sy, mirror])
	if CTF: 
		adw_img = Util.mult_scalar(ctf_2_sum, snr)
		Util.div_filter(adw_img, ctf_abs_sum)
		Util.mul_scalar(adw_img, float(Ng-1)/(nima-1))
		adw_img += float(nima-Ng)/(nima-1)

	# startup
	numr = Numrinit(first_ring, last_ring, rstep, mode) 	#precalculate rings
 	wr = ringwe(numr, mode)

	ref_data = [mask, center, None, None]
	sx_sum = 0.0
	sy_sum = 0.0
	a0 = -1.0e22

	cs = [0.0]*2
	total_iter = 0

	if CUDA:
		from math import log, pi
		RING_LENGTH = 2**(int(log(2*pi*last_ring)/log(2))+1)
		NRING = 2**(int(log(last_ring)/log(2))+1)

	for N_step in xrange(len(xrng)):

		if CUDA:
			R.setup(len(data), nx, nx, RING_LENGTH, NRING, last_ring, step[N_step], int(xrng[N_step]/step[N_step]+0.5), int(yrng[N_step]/step[N_step]+0.5), CTF)
			for im in xrange(len(data)):	R.insert_image(data[im], im)
			if CTF:  R.filter_stack(all_ctf_params)

		msg = "\nX range = %5.2f   Y range = %5.2f   Step = %5.2f\n"%(xrng[N_step], yrng[N_step], step[N_step])
		print_msg(msg)
		for Iter in xrange(max_iter):
			total_iter += 1
			print_msg("Iteration #%4d\n"%(total_iter))
			if CUDA:
				ave1 = model_blank(nx, nx)
				ave2 = model_blank(nx, nx)
				R.sum_oe(all_ctf_params, all_ali_params, ave1, ave2)
				# Comment by Zhengfan Yang on 02/01/10
				# The reason for this step is that in CUDA 2-D FFT, the image is multipled by NX*NY times after
				# FFT and IFFT, so we want to decrease it such that the criterion is in line with non-CUDA version
				# However, this step is not mandatory.
				if CTF:
					ave1 /= (nx*2)**2
					ave2 /= (nx*2)**2
			else:
				ave1, ave2 = sum_oe(data, "a", CTF, EMData())  # pass empty object to prevent calculation of ctf^2
			if CTF:
				tavg_Ng = fft(Util.divn_filter(Util.muln_img(fft(Util.addn_img(ave1, ave2)), adw_img), ctf_2_sum))
				tavg = fft(Util.divn_filter(fft(Util.addn_img(ave1, ave2)), ctf_2_sum))
			else: tavg = (ave1+ave2)/nima

			if outdir:
				tavg.write_image(os.path.join(outdir, "aqc.hdf"), total_iter-1)
				if CTF:
					tavg_Ng.write_image(os.path.join(outdir, "aqc_view.hdf"), total_iter-1)
				frsc = fsc_mask(ave1, ave2, mask, 1.0, os.path.join(outdir, "resolution%03d"%(total_iter)))
			else:
				frsc = fsc_mask(ave1, ave2, mask, 1.0)

			if Fourvar:
				if CTF: vav, rvar = varf2d(data, tavg, mask, "a")
				else: vav, rvar = varf(data, tavg, mask, "a")
				tavg = fft(Util.divn_img(fft(tavg), vav))
				vav_r	= Util.pack_complex_to_real(vav)
				if outdir:
					vav_r.write_image(os.path.join(outdir, "varf.hdf"), total_iter-1)

			ref_data[2] = tavg
			ref_data[3] = frsc
			#  call user-supplied function to prepare reference image, i.e., center and filter it
			if center == -1:
				# When center = -1, which is by default, we use the average center method
				ref_data[1] = 0
				tavg, cs = user_func(ref_data)
				cs[0] = sx_sum/float(nima)
				cs[1] = sy_sum/float(nima)
				tavg = fshift(tavg, -cs[0], -cs[1])
				msg = "Average center x =      %10.3f        Center y       = %10.3f\n"%(cs[0], cs[1])
				print_msg(msg)
			else:
				tavg, cs = user_func(ref_data)

			# write the current filtered average
			if outdir:
				tavg.write_image(os.path.join(outdir, "aqf.hdf"), total_iter-1)

			# a0 should increase; stop algorithm when it decreases.  However, it will depend on filtration, so it is not quite right.
			a1 = tavg.cmp("dot", tavg, dict(negative = 0, mask = ref_data[0]))
			msg = "Criterion %d = %15.8e\n"%(total_iter, a1)
			print_msg(msg)
			if total_iter == len(xrng)*max_iter: break
			if a1 < a0:
				if auto_stop == True: break
			else:	a0 = a1

			if CUDA:
				old_ali_params = all_ali_params[:]
			else:
				old_ali_params = []
			        for im in xrange(nima):
			        	alphan, sxn, syn, mirror, scale = get_params2D(data[im])
		        		old_ali_params.extend([alphan, sxn, syn, mirror])

			if Iter%4 != 0 or total_iter > max_iter*len(xrng)-10: delta = 0.0 
			else: delta = dst
			if CUDA:
				all_ali_params = R.ali2d_single_iter(tavg, all_ali_params, cs[0], cs[1], 1, delta)
				sx_sum = all_ali_params[-2]
				sy_sum = all_ali_params[-1]
				for im in xrange(len(data)):  all_ali_params[im*4+3] = int(all_ali_params[im*4+3])
			else:
				sx_sum, sy_sum = ali2d_single_iter(data, numr, wr, cs, tavg, cnx, cny, xrng[N_step], yrng[N_step], step[N_step], mode, CTF=CTF, delta=delta)

		        pixel_error = 0.0
		        mirror_consistent = 0
			pixel_error_list = []
		        for im in xrange(nima):
		        	if CUDA:
					alpha = all_ali_params[im*4]
					sx = all_ali_params[im*4+1]
					sy = all_ali_params[im*4+2]
					mirror = all_ali_params[im*4+3]
				else:
			        	alpha, sx, sy, mirror, scale = get_params2D(data[im]) 
		        	if old_ali_params[im*4+3] == mirror:
		        		this_error = max_2D_pixel_error(old_ali_params[im*4:im*4+3], [alpha, sx, sy], last_ring)
		        		pixel_error += this_error
					pixel_error_list.append(this_error)
					mirror_consistent += 1
			print_msg("Mirror consistent rate = %6.4f%%\n"%(float(mirror_consistent)/nima*100))
			if mirror_consistent != 0:
				print_msg("Among the mirror consistent images, average pixel error is %0.4f, their distribution is:\n"%(float(pixel_error)/float(mirror_consistent)))
 				region, hist = hist_list(pixel_error_list, 20)	
				for p in xrange(20):
					print_msg("      %8.4f: %5d\n"%(region[p], hist[p]))
			print_msg("\n\n\n")
		if CUDA: R.finish()

	if CUDA:
		for im in xrange(nima):
			set_params2D(data[im], [all_ali_params[im*4], all_ali_params[im*4+1], all_ali_params[im*4+2], all_ali_params[im*4+3], 1.0])

	if outdir:
		drop_image(tavg, os.path.join(outdir, "aqfinal.hdf"))

	if from_ali2d == False: print_end_msg("ali2d_data")

def ali2d_MPI(stack, outdir, maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", dst=0.0, center=-1, maxit=0, CTF=False, snr=1.0, \
			Fourvar=False, Ng=-1, user_func_name="ref_ali2d", CUDA=False, GPUID=""):

	from utilities    import model_circle, model_blank, drop_image, get_image, get_input_from_string
	from utilities    import reduce_EMData_to_root, bcast_EMData_to_all, send_attr_dict, file_type, bcast_number_to_all, bcast_list_to_all
	from statistics   import fsc_mask, sum_oe, hist_list, varf2d_MPI
	from alignment    import Numrinit, ringwe, ali2d_single_iter
	from pixel_error  import max_2D_pixel_error
	from numpy        import reshape, shape
	from fundamentals import fshift, fft, rot_avg_table
	from utilities    import get_params2D, set_params2D
	from utilities    import print_msg, print_begin_msg, print_end_msg
	import os
	import sys
	from mpi 	  import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_gatherv
	from mpi 	  import MPI_SUM, MPI_FLOAT, MPI_INT

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	ftp = file_type(stack)
	
	if outdir:
		if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ali2d_MPI", 1, myid)
		mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		print_begin_msg("ali2d_MPI")
		if outdir:	os.mkdir(outdir)

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	
	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit);
	
	if max_iter == 0:
		max_iter = 10
		auto_stop = True
	else:
		auto_stop = False

	if myid == main_node:
       		if ftp == "bdb":
			from EMAN2db import db_open_dict
			dummy = db_open_dict(stack, True)
		active = EMUtil.get_all_attributes(stack, 'active')
		list_of_particles = []
		for im in xrange(len(active)):
			if active[im]:  list_of_particles.append(im)
		del active
		nima = len(list_of_particles)
	else:
		nima = 0
	nima = bcast_number_to_all(nima, source_node = main_node)
	
	if myid != main_node:
		list_of_particles = [-1]*nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
	list_of_particles = list_of_particles[image_start: image_end]
	
	if Ng == -1:
		Ng = nima
	elif Ng == -2:
		Ng = int(0.98*nima)	

	# read nx and ctf_app (if CTF) and broadcast to all nodes
	if myid == main_node:
		ima = EMData()
		ima.read_image(stack, list_of_particles[0], True)
		nx = ima.get_xsize()
		if CTF:	ctf_app = ima.get_attr_default('ctf_applied', 2)
		del ima
	else:
		nx = 0
		if CTF:	ctf_app = 0
	nx = bcast_number_to_all(nx, source_node = main_node)
	if CTF:
		ctf_app = bcast_number_to_all(ctf_app, source_node = main_node)
		if ctf_app > 0:	ERROR("data cannot be ctf-applied", "ali2d_MPI", 1, myid)

	# default value for the last ring
	if last_ring == -1: last_ring = nx/2-2

	if CUDA:
		GPUID = get_input_from_string(GPUID)
		GPUID = map(int, GPUID)
	
	if myid == main_node:
		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Number of active images     : %d\n"%(nima))
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Inner radius                : %i\n"%(first_ring))
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Ring step                   : %i\n"%(rstep))
		print_msg("X search range              : %s\n"%(xrng))
		print_msg("Y search range              : %s\n"%(yrng))
		print_msg("Translational step          : %s\n"%(step))
		print_msg("Discrete angle used         : %d\n"%(dst))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Use Fourier variance        : %s\n"%(Fourvar))
		print_msg("Number of groups            : %d\n"%(Ng))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		if auto_stop:
			print_msg("Stop iteration with         : criterion\n")
		else:
			print_msg("Stop iteration with         : maxit\n")

		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("User function               : %s\n"%(user_func_name))
		print_msg("Number of processors used   : %d\n"%(number_of_proc))
		print_msg("Using CUDA                  : %s\n"%(CUDA))
		if CUDA:
			print_msg("GPU IDs                     : %s\n"%(GPUID))

	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  
			if myid == main_node:		print_msg("Maskfile                    : %s\n\n"%(maskfile))
			mask = get_image(maskfile)
		else:
			if myid == main_node: 		print_msg("Maskfile                    : user provided in-core mask\n\n")
			mask = maskfile
	else:
		if myid == main_node: 	print_msg("Maskfile                    : default, a circle with radius %i\n\n"%(last_ring))
		mask = model_circle(last_ring, nx, nx)

	cnx  = nx/2+1
 	cny  = cnx
 	mode = "F"
	data = []
	if CTF:
		from filter import filt_ctf
		from morphology   import ctf_img
		ctf_abs_sum = EMData(nx, nx, 1, False)
		ctf_2_sum = EMData(nx, nx, 1, False)
	else:
		ctf_2_sum = None

	from global_def import CACHE_DISABLE
	if CACHE_DISABLE:
		data = EMData.read_images(stack, list_of_particles)
	else:
		for i in xrange(number_of_proc):
			if myid == i:
				data = EMData.read_images(stack, list_of_particles)
			if ftp == "bdb": mpi_barrier(MPI_COMM_WORLD)
		
	if CUDA:
		nGPU = len(GPUID)
		GPUID = GPUID[myid%nGPU]
		R = CUDA_Aligner(GPUID)
		all_ali_params = []
		all_ctf_params = []

	for im in xrange(len(data)):
		data[im].set_attr('ID', list_of_particles[im])
		st = Util.infomask(data[im], mask, False)
		data[im] -= st[0]
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			ctfimg = ctf_img(nx, ctf_params)
			if CUDA:
				all_ctf_params.extend([ctf_params.defocus, ctf_params.cs, ctf_params.voltage, ctf_params.apix, ctf_params.bfactor, ctf_params.ampcont])
			Util.add_img2(ctf_2_sum, ctfimg)
			Util.add_img_abs(ctf_abs_sum, ctfimg)
		if CUDA:
			alpha, sx, sy, mirror, scale = get_params2D(data[im])
			all_ali_params.extend([alpha, sx, sy, mirror])

	if CTF:
		reduce_EMData_to_root(ctf_2_sum, myid, main_node)
		reduce_EMData_to_root(ctf_abs_sum, myid, main_node)
		if myid == main_node:
			adw_img = Util.mult_scalar(ctf_2_sum, snr)
			Util.div_filter(adw_img, ctf_abs_sum)
			Util.mul_scalar(adw_img, float(Ng-1)/(nima-1))
			adw_img += float(nima-Ng)/(nima-1)
	else:  ctf_2_sum = None
	# startup
	numr = Numrinit(first_ring, last_ring, rstep, mode) 	#precalculate rings
 	wr = ringwe(numr, mode)
	
	if myid == main_node:
		# initialize data for the reference preparation function
		ref_data = [mask, center, None, None]
		sx_sum = 0.0
		sy_sum = 0.0
		a0 = -1.0e22
		
	recvcount = []
	disp = []
	for i in xrange(number_of_proc):
		ib, ie = MPI_start_end(nima, number_of_proc, i)
		recvcount.append(ie-ib)
		if i == 0:
			disp.append(0)
		else:
			disp.append(disp[i-1]+recvcount[i-1])

	again = 1
	total_iter = 0
	cs = [0.0]*2

	if CUDA:
		from math import log, pi
		RING_LENGTH = 2**(int(log(2*pi*last_ring)/log(2))+1)
		NRING = 2**(int(log(last_ring)/log(2))+1)

	for N_step in xrange(len(xrng)):

		if CUDA:
			R.setup(len(data), nx, nx, RING_LENGTH, NRING, last_ring, step[N_step], int(xrng[N_step]/step[N_step]+0.5), int(yrng[N_step]/step[N_step]+0.5), CTF)
			for im in xrange(len(data)):	R.insert_image(data[im], im)
			if CTF:  R.filter_stack(all_ctf_params)

		msg = "\nX range = %5.2f   Y range = %5.2f   Step = %5.2f\n"%(xrng[N_step], yrng[N_step], step[N_step])
		if myid == main_node: print_msg(msg)
		for Iter in xrange(max_iter):
			total_iter += 1
			if CUDA:
				ave1 = model_blank(nx, nx)
				ave2 = model_blank(nx, nx)
				R.sum_oe(all_ctf_params, all_ali_params, ave1, ave2)
				# Comment by Zhengfan Yang on 02/01/10
				# The reason for this step is that in CUDA 2-D FFT, the image is multipled by NX*NY times after
				# FFT and IFFT, so we want to decrease it such that the criterion is in line with non-CUDA version
				# However, this step is not mandatory.
				if CTF:
					ave1 /= (nx*2)**2
					ave2 /= (nx*2)**2
			else:
				ave1, ave2 = sum_oe(data, "a", CTF, EMData())  # pass empty object to prevent calculation of ctf^2
			reduce_EMData_to_root(ave1, myid, main_node)
			reduce_EMData_to_root(ave2, myid, main_node)
			if myid == main_node:
				print_msg("Iteration #%4d\n"%(total_iter))
				if CTF: 
					tavg_Ng = fft(Util.divn_filter(Util.muln_img(fft(Util.addn_img(ave1, ave2)), adw_img), ctf_2_sum))
					tavg = fft(Util.divn_filter(fft(Util.addn_img(ave1, ave2)), ctf_2_sum))
				else:	 tavg = (ave1+ave2)/nima
				if outdir:
					tavg.write_image(os.path.join(outdir, "aqc.hdf"), total_iter-1)
					if CTF:
						tavg_Ng.write_image(os.path.join(outdir, "aqc_view.hdf"), total_iter-1)
					frsc = fsc_mask(ave1, ave2, mask, 1.0, os.path.join(outdir, "resolution%03d"%(total_iter)))
				else:
					frsc = fsc_mask(ave1, ave2, mask, 1.0)
			else:
				tavg =  model_blank(nx, nx)
			del ave1, ave2
				
			if Fourvar:  
				bcast_EMData_to_all(tavg, myid, main_node)
				vav, rvar = varf2d_MPI(myid, data, tavg, mask, "a", CTF)

			if myid == main_node:
				if Fourvar:
					tavg    = fft(Util.divn_img(fft(tavg), vav))
					vav_r	= Util.pack_complex_to_real(vav)
					if outdir:
						vav_r.write_image(os.path.join(outdir, "varf.hdf"), total_iter-1)

				ref_data[2] = tavg
				ref_data[3] = frsc

				#  call user-supplied function to prepare reference image, i.e., center and filter it
				if center == -1:
					# When center = -1, which is by default, we use the average center method
					ref_data[1] = 0
					tavg, cs = user_func(ref_data)
					cs[0] = float(sx_sum)/nima
					cs[1] = float(sy_sum)/nima
					tavg = fshift(tavg, -cs[0], -cs[1])
					msg = "Average center x =      %10.3f        Center y       = %10.3f\n"%(cs[0], cs[1])
					print_msg(msg)
				else:
					tavg, cs = user_func(ref_data)

				# write the current filtered average
				if outdir:
					tavg.write_image(os.path.join(outdir, "aqf.hdf"), total_iter-1)

				# a0 should increase; stop algorithm when it decreases.    However, the result will depend on filtration, so it is not quite right.
				a1 = tavg.cmp("dot", tavg, dict(negative = 0, mask = ref_data[0]))
				msg = "Criterion %d = %15.8e\n"%(total_iter, a1)
				print_msg(msg)
				
				if a1 < a0:
					if auto_stop: 	again = 0
				else:	a0 = a1
			else:
				tavg = model_blank(nx, nx)
				cs = [0.0]*2

			if auto_stop:
				again = mpi_bcast(again, 1, MPI_INT, main_node, MPI_COMM_WORLD)
				if int(again[0]) == 0: break

			if Fourvar:  del vav
			bcast_EMData_to_all(tavg, myid, main_node)
			cs = mpi_bcast(cs, 2, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			cs = map(float, cs)
			if total_iter != max_iter*len(xrng):
				if CUDA:
					old_ali_params = all_ali_params[:]
				else:
					old_ali_params = []
				        for im in xrange(len(data)):  
						alpha, sx, sy, mirror, scale = get_params2D(data[im])
						old_ali_params.extend([alpha, sx, sy, mirror])

				if Iter%4 != 0 or total_iter > max_iter*len(xrng)-10: delta = 0.0 
				else: delta = dst
				if CUDA:
					all_ali_params = R.ali2d_single_iter(tavg, all_ali_params, cs[0], cs[1], 1, delta)
					sx_sum = all_ali_params[-2]
					sy_sum = all_ali_params[-1]
					for im in xrange(len(data)):  all_ali_params[im*4+3] = int(all_ali_params[im*4+3])
				else:	
					sx_sum, sy_sum = ali2d_single_iter(data, numr, wr, cs, tavg, cnx, cny, xrng[N_step], yrng[N_step], step[N_step], mode, CTF=CTF, delta=delta)
					
				sx_sum = mpi_reduce(sx_sum, 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
				sy_sum = mpi_reduce(sy_sum, 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)

			        pixel_error       = 0.0
			        mirror_consistent = 0
				pixel_error_list  = []
			        for im in xrange(len(data)):
			        	if CUDA:
						alpha = all_ali_params[im*4]
						sx = all_ali_params[im*4+1]
						sy = all_ali_params[im*4+2]
						mirror = all_ali_params[im*4+3]
					else:
						alpha, sx, sy, mirror, scale = get_params2D(data[im])
			        	if old_ali_params[im*4+3] == mirror:
		        			this_error = max_2D_pixel_error(old_ali_params[im*4:im*4+3], [alpha, sx, sy], last_ring)
		        			pixel_error += this_error
						pixel_error_list.append(this_error)
						mirror_consistent += 1
					else:
						pixel_error_list.append(-1)
				mirror_consistent = mpi_reduce(mirror_consistent, 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
				pixel_error       = mpi_reduce(pixel_error, 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
				pixel_error_list  = mpi_gatherv(pixel_error_list, len(data), MPI_FLOAT, recvcount, disp, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				if myid == main_node:
					print_msg("Mirror consistency rate = %8.4f%%\n"%(float(mirror_consistent)/nima*100))
					if mirror_consistent!=0:
						print_msg("Among the mirror-consistent images, average of pixel errors is %0.4f, and their distribution is:\n"%(float(pixel_error)/float(mirror_consistent)))
						pixel_error_list = map(float, pixel_error_list)
						for i in xrange(nima-1, -1, -1):
							if pixel_error_list[i] < 0:  del pixel_error_list[i]
						region, hist = hist_list(pixel_error_list, 20)	
						for p in xrange(20):
							print_msg("      %10.6f: %5d\n"%(region[p], hist[p]))
					print_msg("\n\n\n")
		if CUDA: R.finish()

	if CUDA:
		for im in xrange(len(data)):
			set_params2D(data[im], [all_ali_params[im*4], all_ali_params[im*4+1], all_ali_params[im*4+2], all_ali_params[im*4+3], 1.0])

	if myid == main_node and outdir:  drop_image(tavg, os.path.join(outdir, "aqfinal.hdf"))
	# write out headers and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
	par_str = ["xform.align2d", "ID"]
	if myid == main_node:
		from utilities import file_type
		if(file_type(stack) == "bdb"):
			from utilities import recv_attr_dict_bdb
			recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
		else:
			from utilities import recv_attr_dict
			recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	else:           send_attr_dict(main_node, data, par_str, image_start, image_end)
	if myid == main_node: print_end_msg("ali2d_MPI")

'''
def ORGali2d_c(stack, outdir, maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", center=-1, maxit=0, \
		CTF=False, snr=1.0, Fourvar = False, user_func_name="ref_ali2d", CUDA=False, GPU=0, MPI=False):
	if MPI:
		ali2d_c_MPI(stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, center, maxit, CTF, snr, Fourvar, user_func_name, CUDA, GPU)
		return

	from utilities    import model_circle, drop_image, get_image, get_input_from_string, get_params2D
	from statistics   import fsc_mask, sum_oe, hist_list
	from alignment    import Numrinit, ringwe, ali2d_single_iter
	from pixel_error  import max_2D_pixel_error
	from filter       import filt_ctf, filt_table, filt_tophatb
	from fundamentals import fshift
	from utilities    import print_begin_msg, print_end_msg, print_msg
	from fundamentals import fft, rot_avg_table
	from utilities    import write_text_file, file_type
	import os
		
	print_begin_msg("ali2d_c")

	if os.path.exists(outdir):   ERROR('Output directory exists, please change the name and restart the program', " ORGali2d_c", 1)
	os.mkdir(outdir)

	import user_functions
	user_func = user_functions.factory[user_func_name]

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)

	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit);

	if max_iter == 0:
		max_iter  = 10
		auto_stop = True
	else:
		auto_stop = False

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Inner radius                : %i\n"%(first_ring))

	if(file_type(stack) == "bdb"):
		from EMAN2db import db_open_dict
		dummy = db_open_dict(stack, True)
	active = EMUtil.get_all_attributes(stack, 'active')
	list_of_particles = []
	for im in xrange(len(active)):
		if active[im]:  list_of_particles.append(im)
	del active
	nima = len(list_of_particles)
	ima  = EMData()
	ima.read_image(stack, list_of_particles[0], True)
	nx = ima.get_xsize()

	# default value for the last ring
	if last_ring == -1:  last_ring = nx/2-2

	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Ring step                   : %i\n"%(rstep))
	print_msg("X search range              : %s\n"%(xrng))
	print_msg("Y search range              : %s\n"%(yrng))
	print_msg("Translational step          : %s\n"%(step))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("Use Fourier variance        : %s\n"%(Fourvar))
	print_msg("CTF correction              : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	if auto_stop:  print_msg("Stop iteration with         : criterion\n")
	else:           print_msg("Stop iteration with         : maxit\n")
	print_msg("User function               : %s\n"%(user_func_name))

	if maskfile:
		import	types
		if type(maskfile) is types.StringType:
			print_msg("Maskfile                    : %s\n\n"%(maskfile))
			mask=get_image(maskfile)
		else:
			print_msg("Maskfile                    : user provided in-core mask\n\n")
			mask = maskfile
	else :
		print_msg("Maskfile                    : default, a circle with radius %i\n\n"%(last_ring))
		mask = model_circle(last_ring, nx, nx)

	cnx = nx/2+1
 	cny = cnx
 	mode = "F"
	data = []
	if CTF:
		ctf_params = ima.get_attr("ctf")
		if ima.get_attr_default('ctf_applied', 2) > 0:	ERROR("data cannot be ctf-applied", "ORGali2d_c", 1)
		from morphology   import ctf_img
		ctf_2_sum = EMData(nx, nx, 1, False)
	else:
		ctf_2_sum = None
	if  Fourvar:
		from statistics   import add_ave_varf

	del ima
	data = EMData.read_images(stack, list_of_particles)
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		st = Util.infomask(data[im], mask, False)
		data[im] -= st[0]
		if CTF:
			ctf_params = data[im].get_attr("ctf")
	 		Util.add_img2(ctf_2_sum, ctf_img(nx, ctf_params))
	if(CTF and (not Fourvar)):  ctf_2_sum += 1.0/snr  # note this is complex addition (1.0/snr,0.0)
	# startup
	numr = Numrinit(first_ring, last_ring, rstep, mode) 	#precalculate rings
 	wr = ringwe(numr, mode)
	
	# initialize data for the reference preparation function
	#  mask can be modified in user_function
	ref_data = [mask, center, None, None]

	cs = [0.0]*2
	# iterate
	total_iter = 0
	a0 = -1e22
	sx_sum = 0.0
	sy_sum = 0.0
	for N_step in xrange(len(xrng)):
		msg = "\nX range = %5.2f   Y range = %5.2f   Step = %5.2f\n"%(xrng[N_step], yrng[N_step], step[N_step])
		print_msg(msg)
		for Iter in xrange(max_iter):
			total_iter += 1
			print_msg("Iteration #%4d\n"%(total_iter))
			if  Fourvar:  
				tavg, ave1, ave2, vav, sumsq = add_ave_varf(data, mask, "a", CTF, ctf_2_sum)
				# write the current average
				fft(tavg).write_image(os.path.join(outdir, "aqc.hdf"), total_iter-1)
				tavg    = fft(Util.divn_img(tavg, vav))

				vav_r	= Util.pack_complex_to_real(vav)
				vav_r.write_image(os.path.join(outdir, "varf.hdf"), total_iter-1)
				sumsq_r = Util.pack_complex_to_real(sumsq)
				rvar    = rot_avg_table(vav_r)
				rsumsq  = rot_avg_table(sumsq_r)
				frsc = []
				freq = []
				for i in xrange(len(rvar)):
					qt = max(0.0, rsumsq[i]/rvar[i] - 1.0)
					frsc.append(qt/(qt+1.0))
					freq.append(float(i)/(len(rvar)-1)*0.5)
				frsc = [freq, frsc]
				del freq
				write_text_file(frsc, os.path.join(outdir, "resolution%03d"%(total_iter)))
			else:
				ave1, ave2 = sum_oe(data, "a", CTF, ctf_2_sum)
				if CTF:  tavg = fft(Util.divn_img(fft(Util.addn_img(ave1, ave2)), ctf_2_sum))
				else:	 tavg = (ave1+ave2)/nima
				tavg.write_image(os.path.join(outdir, "aqc.hdf"), total_iter-1)
				frsc = fsc_mask(ave1, ave2, mask, 1.0, os.path.join(outdir, "resolution%03d"%(total_iter)))

			ref_data[2] = tavg
			ref_data[3] = frsc
			#  call user-supplied function to prepare reference image, i.e., center and filter it
			if center == -1:
				# When center = -1, which is by default, we use the average center method
				ref_data[1] = 0
				tavg, cs = user_func(ref_data)
				cs[0] = sx_sum/float(nima)
				cs[1] = sy_sum/float(nima)
				tavg = fshift(tavg, -cs[0], -cs[1])
				msg = "Average center x =      %10.3f        Center y       = %10.3f\n"%(cs[0], cs[1])
				print_msg(msg)
			else:
				tavg, cs = user_func(ref_data)

			# write the current filtered average
			tavg.write_image(os.path.join(outdir, "aqf.hdf"), total_iter-1)

			# a0 should increase; stop algorithm when it decreases.
			if Fourvar:  
				Util.div_filter(sumsq, vav)
				sumsq = filt_tophatb(sumsq, 0.01, 0.49)
				a1 = Util.infomask(sumsq, None, True)
				a1 = a1[0]
			else:
				a1 = tavg.cmp("dot", tavg, dict(negative = 0, mask = ref_data[0]))
			msg = "Criterion = %15.8e\n"%(a1)
			print_msg(msg)
			if total_iter == len(xrng)*max_iter: break
			if a1 < a0:
				if auto_stop == True: break
			else:	a0 = a1

			old_ali_params = []
		        for im in xrange(nima):
		        	alphan, sxn, syn, mirror, scale = get_params2D(data[im])
		        	old_ali_params.append([alphan, sxn, syn, mirror, scale])

			sx_sum, sy_sum = ali2d_single_iter(data, numr, wr, cs, tavg, cnx, cny, xrng[N_step], yrng[N_step], step[N_step], mode, CTF)

		        pixel_error = 0.0
		        mirror_changed = 0
			pixel_error_list = []
		        for im in xrange(nima):
		        	alphan, sxn, syn, mirror, scale = get_params2D(data[im]) 
		        	if old_ali_params[im][3] == mirror:
		        		this_error = max_2D_pixel_error(old_ali_params[im][0], old_ali_params[im][1], old_ali_params[im][2], alphan, sxn, syn, last_ring)
		        		pixel_error += this_error
					pixel_error_list.append(this_error)
		        	else:
		        		mirror_changed += 1
			print_msg("Mirror changed = %6.4f%%\n"%(float(mirror_changed)/nima*100))
			print_msg("Among the mirror consistent images, average pixel error is %0.4f, their distribution is:\n"%(pixel_error/float(nima-mirror_changed)))
 			region, hist = hist_list(pixel_error_list, 20)	
			for p in xrange(20):
				print_msg("      %8.4f: %5d\n"%(region[p], hist[p]))
			print_msg("\n\n\n")
			
	drop_image(tavg, os.path.join(outdir, "aqfinal.hdf"))
	# write out headers
	from utilities import write_headers
	write_headers(stack, data, list_of_particles)
	print_end_msg("ali2d_c")

def ORGali2d_c_MPI(stack, outdir, maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", center=-1, maxit=0, CTF=False, snr=1.0, \
			Fourvar = False, user_func_name="ref_ali2d", CUDA=False, GPU=0):

	from utilities    import model_circle, model_blank, drop_image, get_image, get_input_from_string
	from utilities    import reduce_EMData_to_root, bcast_EMData_to_all, send_attr_dict, file_type, bcast_number_to_all, bcast_list_to_all
	from statistics   import fsc_mask, sum_oe, add_ave_varf_MPI, hist_list
	from alignment    import Numrinit, ringwe, ali2d_single_iter
	from pixel_error  import max_2D_pixel_error
	from filter       import filt_table, filt_ctf, filt_tophatb
	from numpy        import reshape, shape
	from fundamentals import fshift, fft, rot_avg_table
	from utilities    import write_text_file, get_params2D, set_params2D
	from utilities    import print_msg, print_begin_msg, print_end_msg
	import os
	import sys
	from mpi 	  import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_gatherv
	from mpi 	  import MPI_SUM, MPI_FLOAT, MPI_INT

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	ftp = file_type(stack)

	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ORGali2d_c_MPI ", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		print_begin_msg("ali2d_c_MPI")
		os.mkdir(outdir)

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	
	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit);
	
	if max_iter == 0:
		max_iter = 10
		auto_stop = True
	else:
		auto_stop = False

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]
		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Inner radius                : %i\n"%(first_ring))
	
	if myid == main_node:
       		if(file_type(stack) == "bdb"):
			from EMAN2db import db_open_dict
			dummy = db_open_dict(stack, True)
		active = EMUtil.get_all_attributes(stack, 'active')
		list_of_particles = []
		for im in xrange(len(active)):
			if active[im]:  list_of_particles.append(im)
		del active
		nima = len(list_of_particles)
	else:
		nima = 0
	nima = bcast_number_to_all(nima, source_node = main_node)
	
	if myid != main_node:
		list_of_particles = [-1]*nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)
	
	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
	list_of_particles = list_of_particles[image_start: image_end]

	ima = EMData()
	for i in xrange(number_of_proc):
		if myid == i:
			ima.read_image(stack, list_of_particles[0], True)
		if ftp == "bdb": mpi_barrier(MPI_COMM_WORLD)
	nx = ima.get_xsize()
	# default value for the last ring
	if last_ring == -1: last_ring = nx/2-2

	if myid == main_node:
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Ring step                   : %i\n"%(rstep))
		print_msg("X search range              : %s\n"%(xrng))
		print_msg("Y search range              : %s\n"%(yrng))
		print_msg("Translational step          : %s\n"%(step))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Use Fourier variance        : %s\n"%(Fourvar))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		if auto_stop:
			print_msg("Stop iteration with         : criterion\n")
		else:
			print_msg("Stop iteration with         : maxit\n")
		print_msg("User function               : %s\n"%(user_func_name))
		print_msg("Number of processors used   : %d\n"%(number_of_proc))
		print_msg("Using CUDA                  : %s\n"%(CUDA))
		print_msg("Number of GPUs              : %d\n"%(GPU))

		
	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  
			if myid == main_node:		print_msg("Maskfile                    : %s\n\n"%(maskfile))
			mask = get_image(maskfile)
		else:
			if myid == main_node: 		print_msg("Maskfile                    : user provided in-core mask\n\n")
			mask = maskfile
	else : 
		if myid == main_node: 	print_msg("Maskfile                    : default, a circle with radius %i\n\n"%(last_ring))
		mask = model_circle(last_ring, nx, nx)

	cnx  = nx/2+1
 	cny  = cnx
 	mode = "F"
	data = []
	if CTF:
		ctf_params = ima.get_attr("ctf")
		if ima.get_attr_default('ctf_applied', 2) > 0:	ERROR("data cannot be ctf-applied", "ORGali2d_c_MPI", 1,myid)
		from filter import filt_ctf
		from morphology   import ctf_img
		ctf_2_sum = EMData(nx, nx, 1, False)
	else:
		ctf_2_sum = None
	if  Fourvar:
		from statistics   import add_ave_varf

	del ima

	for i in xrange(number_of_proc):
		if myid == i: 
			data = EMData.read_images(stack, list_of_particles)
		if ftp == "bdb": mpi_barrier(MPI_COMM_WORLD)
	
	if CUDA:
		GPUID = myid%GPU
		all_ali_params = []
		all_ctf_params = []
	for im in xrange(len(data)):
		data[im].set_attr('ID', list_of_particles[im])
		st = Util.infomask(data[im], mask, False)
		data[im] -= st[0]
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			Util.add_img2(ctf_2_sum, ctf_img(nx, ctf_params))
			if CUDA:
				all_ctf_params.append(ctf_params.defocus)
				all_ctf_params.append(ctf_params.cs)
				all_ctf_params.append(ctf_params.voltage)
				all_ctf_params.append(ctf_params.apix)
				all_ctf_params.append(ctf_params.bfactor)
				all_ctf_params.append(ctf_params.ampcont)
		if CUDA:
			alpha, sx, sy, mirror, scale = get_params2D(data[im])
			all_ali_params.append(alpha)
			all_ali_params.append(sx)
			all_ali_params.append(sy)
			all_ali_params.append(mirror)

	if CTF:
		reduce_EMData_to_root(ctf_2_sum, myid, main_node)
		if(not Fourvar):  ctf_2_sum += 1.0/snr  # note this is complex addition (1.0/snr,0.0)
	else:  ctf_2_sum = None
	# startup
	numr = Numrinit(first_ring, last_ring, rstep, mode) 	#precalculate rings
 	wr = ringwe(numr, mode)
	
	if myid == main_node:
		# initialize data for the reference preparation function
		ref_data = [mask, center, None, None]
		sx_sum = 0.0
		sy_sum = 0.0
		a0 = -1.0e22
		
	recvcount = []
	disp = []
	for i in xrange(number_of_proc):
		ib, ie = MPI_start_end(nima, number_of_proc, i)
		recvcount.append(ie-ib)
		if i == 0:
			disp.append(0)
		else:
			disp.append(disp[i-1]+recvcount[i-1])

	again = True
	total_iter = 0
	cs = [0.0]*2

	for N_step in xrange(len(xrng)):

		if CUDA:
			R = CUDA_Aligner()
			R.setup(len(data), nx, nx, 256, 32, last_ring, step[N_step], int(xrng[N_step]/step[N_step]+0.5), int(yrng[N_step]/step[N_step]+0.5), CTF)
			for im in xrange(len(data)):	R.insert_image(data[im], im)
			if CTF:  R.filter_stack(all_ctf_params, GPUID)
					
		msg = "\nX range = %5.2f   Y range = %5.2f   Step = %5.2f\n"%(xrng[N_step], yrng[N_step], step[N_step])
		if myid == main_node: print_msg(msg)
		for Iter in xrange(max_iter):
			total_iter += 1
			if  Fourvar:  
				tavg, ave1, ave2, vav, sumsq = add_ave_varf_MPI(myid, data, mask, "a", CTF, ctf_2_sum)
			else:
				if CUDA:
					ave1 = model_blank(nx, nx)
					ave2 = model_blank(nx, nx)
					R.sum_oe(all_ctf_params, all_ali_params, ave1, ave2, GPUID)
				else:
					ave1, ave2 = sum_oe(data, "a", CTF, EMData())  # pass empty object to prevent calculation of ctf^2
				reduce_EMData_to_root(ave1, myid, main_node)
				reduce_EMData_to_root(ave2, myid, main_node)

			if myid == main_node:
				print_msg("Iteration #%4d\n"%(total_iter))
				if Fourvar:
					fft(tavg).write_image(os.path.join(outdir, "aqc.hdf"), total_iter-1)
					tavg    = fft(Util.divn_img(tavg, vav))
					vav_r	= Util.pack_complex_to_real(vav)
					vav_r.write_image(os.path.join(outdir, "varf.hdf"), total_iter-1)
					rvar	= rot_avg_table(vav_r)
					rsumsq  = rot_avg_table(Util.pack_complex_to_real(sumsq))
					frsc = []
					freq = []
					for i in xrange(len(rvar)):
						qt = max(0.0, rsumsq[i]/rvar[i] - 1.0)
						frsc.append(qt/(qt+1.0))
						freq.append(float(i)/(len(rvar)-1)*0.5)
					frsc = [freq, frsc]
					del freq
					write_text_file(frsc, os.path.join(outdir, "resolution%03d"%(total_iter)) )
				else:
					if CTF:  tavg = fft(Util.divn_img(fft(Util.addn_img(ave1, ave2)), ctf_2_sum))
					else:	 tavg = (ave1+ave2)/nima
					tavg.write_image(os.path.join(outdir, "aqc.hdf"), total_iter-1)
					frsc = fsc_mask(ave1, ave2, mask, 1.0, os.path.join(outdir, "resolution%03d"%(total_iter)))

				ref_data[2] = tavg
				ref_data[3] = frsc
				
				#  call user-supplied function to prepare reference image, i.e., center and filter it
				if center == -1:
					# When center = -1, which is by default, we use the average center method
					ref_data[1] = 0
					tavg, cs = user_func(ref_data)
					cs[0] = float(sx_sum)/nima
					cs[1] = float(sy_sum)/nima
					tavg = fshift(tavg, -cs[0], -cs[1])
					msg = "Average center x =      %10.3f        Center y       = %10.3f\n"%(cs[0], cs[1])
					print_msg(msg)
				else:
					tavg, cs = user_func(ref_data)

				# write the current filtered average
				tavg.write_image(os.path.join(outdir, "aqf.hdf"), total_iter-1)

				# a0 should increase; stop algorithm when it decreases.    
				if Fourvar:  
					Util.div_filter(sumsq, vav)
					sumsq = filt_tophatb(sumsq, 0.01, 0.49)
					a1 = Util.infomask(sumsq, None, True)
					a1 = a1[0]
				else:
					a1 = tavg.cmp("dot", tavg, dict(negative = 0, mask = ref_data[0]))
				msg = "Criterion %d = %15.8e\n"%(total_iter, a1)
				print_msg(msg)
				
				if a1 < a0:
					if auto_stop: 
						again = False
						break
				else:	a0 = a1
			else:
				tavg = EMData(nx, nx, 1, True)
				cs = [0.0]*2

			bcast_EMData_to_all(tavg, myid, main_node)
			cs = mpi_bcast(cs, 2, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			cs = map(float, cs)
			if auto_stop:
				again = mpi_bcast(again, 1, MPI_INT, main_node, MPI_COMM_WORLD)
				if not again: break
			if total_iter != max_iter*len(xrng):
				if CUDA:
					old_ali_params = all_ali_params[:]
				else:
					old_ali_params = []
				        for im in xrange(len(data)):  
						alpha, sx, sy, mirror, scale = get_params2D(data[im])
						old_ali_params.append(alpha)
						old_ali_params.append(sx)
						old_ali_params.append(sy)
						old_ali_params.append(mirror)

				if CUDA:
					all_ali_params = R.ali2d_single_iter(tavg, all_ali_params, cs[0], cs[1], GPUID, 1)
					sx_sum = all_ali_params[-2]
					sy_sum = all_ali_params[-1]
					for im in xrange(len(data)):  all_ali_params[im*4+3] = int(all_ali_params[im*4+3])
				else:	
					sx_sum, sy_sum = ali2d_single_iter(data, numr, wr, cs, tavg, cnx, cny, xrng[N_step], yrng[N_step], step[N_step], mode, CTF)
					
				sx_sum = mpi_reduce(sx_sum, 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
				sy_sum = mpi_reduce(sy_sum, 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)

			        pixel_error = 0.0
			        mirror_consistent = 0
				pixel_error_list = []
			        for im in xrange(len(data)):
			        	if CUDA:
						alpha = all_ali_params[im*4]
						sx = all_ali_params[im*4+1]
						sy = all_ali_params[im*4+2]
						mirror = all_ali_params[im*4+3]
					else:
						alpha, sx, sy, mirror, scale = get_params2D(data[im]) 
			        	if old_ali_params[im*4+3] == mirror:
		        			this_error = max_2D_pixel_error(old_ali_params[im*4:im*4+3], [alpha, sx, sy], last_ring)
		        			pixel_error += this_error
						pixel_error_list.append(this_error)
						mirror_consistent += 1
					else:
						pixel_error_list.append(-1)
				mirror_consistent = mpi_reduce(mirror_consistent, 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
				pixel_error = mpi_reduce(pixel_error, 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
				pixel_error_list = mpi_gatherv(pixel_error_list, len(data), MPI_FLOAT, recvcount, disp, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				if myid == main_node:
					print_msg("Mirror consistent rate = %8.4f%%\n"%(float(mirror_consistent)/nima*100))
					print_msg("Among the mirror consistent images, average pixel error is %0.4f, their distribution is:\n"%(float(pixel_error)/mirror_consistent))
					pixel_error_list = map(float, pixel_error_list)
					for i in xrange(nima-1, -1, -1):
						if pixel_error_list[i] < 0:  del pixel_error_list[i]
					region, hist = hist_list(pixel_error_list, 20)	
					for p in xrange(20):
						print_msg("      %10.6f: %5d\n"%(region[p], hist[p]))
					print_msg("\n\n\n")
		if CUDA: R.finish()

	if CUDA:
		for im in xrange(len(data)):
			set_params2D(data[im], [all_ali_params[im*4], all_ali_params[im*4+1], all_ali_params[im*4+2], all_ali_params[im*4+3], 1.0])

	if myid == main_node:  drop_image(tavg, os.path.join(outdir, "aqfinal.hdf"))
	# write out headers and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
	par_str = ["xform.align2d", "ID"]
	if myid == main_node:
		from utilities import file_type
		if(file_type(stack) == "bdb"):
			from utilities import recv_attr_dict_bdb
			recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
		else:
			from utilities import recv_attr_dict
			recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	else:           send_attr_dict(main_node, data, par_str, image_start, image_end)
	if myid == main_node: print_end_msg("ali2d_c_MPI")
'''

def local_ali2d(stack, outdir, maskfile = None, ou = -1, br = 1.75, center = 1, eps = 0.001, maxit = 10, CTF = False, snr = 1.0, user_func_name="ref_ali2d"):
	"""
		Name
			local_ali2d - Perform local refinement of 2-D alignment of an image series
		Input
			stack: set of 2-D images in a stack file (format hdf), images have to be squares 
			maskfile: optional mask file to use internally during alignment
			outer_radius: outer radius for rotational correlation < int(nx/2)-1 
			br: brackets for the search of orientation parameters
			center_type: center the average
			epsilon: stopping criterion
			max_iter: maximum number of iterations the program will perform 
			CTF: default no CTF
			snr: Signal-to-Noise Ratio of the data
			function: name of the user-supplied function
		Output
			output_directory: directory name into which the output files will be written.
			header: the alignment parameters are stored in the headers of input files as 'alpha', 'sx', 'sy', 'mirror'
	"""
# 2D alignment using amoeba and gridding interpolation
	from alignment    	import kbt
	from utilities    	import model_circle, amoeba, compose_transform2, drop_image, get_arb_params, get_image, get_params2D, set_params2D
	from alignment    	import fine_2D_refinement, crit2d
	from statistics   	import add_oe_series, fsc_mask
	from filter 		import filt_from_fsc_bwt,filt_table
	from morphology         import ctf_2, ctf_1d
	import os
	import sys
	import types
	output = sys.stdout
	
	from utilities import print_begin_msg, print_end_msg, print_msg
	
	# create the output directory, if it does not existm
	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "local_ali2d", 1)
	os.mkdir(outdir)
	
	max_iter  = int(maxit)
	last_ring = int(ou)

	print_begin_msg("local_ali2d")	
	
	import user_functions
	user_func = user_functions.factory[user_func_name]
	
	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Output directory            : %s\n"%(outdir))

	ima = EMData()
	ima.read_image(stack, 0, True)
	nx = ima.get_xsize()
	# default value for the last ring
	if (last_ring == -1): last_ring=nx//2-2

	if maskfile:
		if(type(maskfile) is types.StringType):
			print_msg("Maskfile                    : %s\n"%(maskfile))
			mask = get_image(maskfile)
		else:
			print_msg("Maskfile                    : user provided in-core mask\n")
			mask = maskfile
	else :
		print_msg("Maskfile                    : default, a circle with radius %i\n"%(last_ring))
		mask = model_circle(last_ring, nx, nx)
	# initialize data for the reference preparation function
	ref_data = [mask, center]

	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Search range                : %-5.2f\n"%(br))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("Error tolerance             : %f\n"%(eps))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("CTF correction              : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n\n"%(snr))
	
	low = 0.5

	# read images
	del ima
	data = EMData.read_images(stack)
	nima = len(data)
	if(CTF):
		ctf_params = data[0].get_attr("ctf")
		data_had_ctf = data[0].get_attr("ctf_applied")
		ctm = ctf_1d(nx, ctf_params)
		lctf = len(ctm)
		ctf2 = []
		ctf2.append([0.0]*lctf)
		ctf2.append([0.0]*lctf)
		ctfb2 = [0.0]*lctf
		for im in xrange(nima):
			ctf_params = data[im].get_attr( "ctf" )
			ctm = ctf_2(nx, ctf_params)
			k = im%2
			for i in xrange(lctf):  ctf2[k][i] += ctm[i]
			if(data[im].get_attr("ctf_applied") == 0):
				st = Util.infomask(data[im], mask, False)
				data[im] -= st[0]
				from filter import filt_ctf
				data[im] = filt_ctf(data[im], ctf_params)
				data[im].set_attr('ctf_applied', 1)
		for i in xrange(lctf):
			ctfb2[i] = 1.0/(ctf2[0][i] + ctf2[1][i] + 1.0/snr)
			for k in xrange(2):
				ctf2[k][i] = 1.0/(ctf2[k][i] + 1.0/snr)
	#calculate averages
	av1, av2 = add_oe_series(data)
	tavg = Util.addn_img(av1, av2)
	if(CTF):
		from filter import filt_table
		tavg = filt_table(tavg, ctfb2)
	else:
		tavg /= nima

	Util.mul_img( tavg, mask)
	a0 = tavg.cmp("dot", tavg, {"negative":0, "mask":mask})
	msg = "Initial criterion = %-20.7e\n"%(a0)	
	print_msg(msg)
	# do the alignment

	for Iter in xrange(max_iter):
		again = False
		fine_2D_refinement(data, br, ref_data[0], tavg)

		# calculate total average using current alignment parameters
		av1, av2 = add_oe_series(data)
		if(CTF):
			tavg = filt_table(Util.addn_img(av1, av2), ctfb2)
			av1  = filt_table(av1, ctf2[0])
			av2  = filt_table(av2, ctf2[1])
		else:
			tavg = (av1 + av2)/nima
		drop_image(tavg, os.path.join(outdir, "aqe_%03d.hdf"%(Iter)))

		frsc = fsc_mask(av1, av2, ref_data[0], 1.0, os.path.join(outdir, "dre%03d"%Iter))

		ref_data.append(tavg)
		ref_data.append(frsc)
		#  call user-supplied function to prepare reference image, i.e., center and filter it
		tavg, cs = user_func( ref_data )
		del ref_data[3]
		del ref_data[2]

		if center:
			#  apply centering parameters to shifts
			for im in xrange(nima):
				alpha, sx, sy, mirror, scale    = get_params2D(data[im])
				alphan, sxn, syn, scale = compose_transform2(alpha, sx, sy, 1.0, 0.0, -cs[0], -cs[1], 1.0)
				set_params2D( data[im], [alphan, sxn, syn, mirror, scale])

		# write current average
		a1 = tavg.cmp("dot", tavg, {"negative":0, "mask":ref_data[0]})
		msg = "ITERATION #%3d        criterion = %20.7e\n"%(Iter+1,a1)
		print_msg(msg)
		# write the current average
		drop_image(tavg, os.path.join(outdir, "aqf_%03d.hdf"%(Iter)))

	if(CTF and data_had_ctf == 0):
		for im in xrange(nima): data[im].set_attr('ctf_applied', 0)
	from utilities import write_headers
	write_headers(stack, data, range(nima))
	print_end_msg("local_ali2d")


def mref_ali2d(stack, refim, outdir, maskfile=None, ir=1, ou=-1, rs=1, xrng=0, yrng=0, step=1, center=1, maxit=0, CTF=False, snr=1.0, user_func_name="ref_ali2d", rand_seed=1000, MPI=False):
	"""
		Name
			mref_ali2d - Perform 2-D multi-reference alignment of an image series
		Input
			stack: set of 2-D images in a stack file, images have to be squares
			refim: set of initial reference 2-D images in a stack file 
			maskfile: optional maskfile to be used in the alignment
			inner_radius: inner radius for rotational correlation > 0
			outer_radius: outer radius for rotational correlation < nx/2-1
			ring_step: step between rings in rotational correlation >0
			x_range: range for translation search in x direction, search is +/xr 
			y_range: range for translation search in y direction, search is +/yr 
			translation_step: step of translation search in both directions
			center: center the average
			max_iter: maximum number of iterations the program will perform
			CTF: if this flag is set, the program will use CTF information provided in file headers
			snr: signal-to-noise ratio of the data
			rand_seed: the seed used for generating random numbers
			MPI: whether to use MPI version
		Output
			output_directory: directory name into which the output files will be written.
			header: the alignment parameters are stored in the headers of input files as 'xform.align2d'.
	"""
# 2D multi-reference alignment using rotational ccf in polar coordinates and quadratic interpolation
	if MPI:
		mref_ali2d_MPI(stack, refim, outdir, maskfile, ir, ou, rs, xrng, yrng, step, center, maxit, CTF, snr, user_func_name, rand_seed)
		return

	from utilities      import   model_circle, combine_params2, inverse_transform2, drop_image, get_image
	from utilities	    import   center_2D, get_im, get_params2D, set_params2D
	from statistics     import   fsc
	from alignment      import   Numrinit, ringwe, Applyws, fine_2D_refinement
	from fundamentals   import   rot_shift2D, fshift
	from morphology     import   ctf_2
	from filter         import   filt_btwl, filt_params
	from random         import   seed, randint
	import os
	import sys

	from utilities      import   print_begin_msg, print_end_msg, print_msg
	
	# create the output directory, if it does not exist
	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "mref_ali2d", 1)
	os.mkdir(outdir)
	
	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit)
	if max_iter == 0:
		max_iter  = 10
		auto_stop = True
	else:
		auto_stop = False

	print_begin_msg("mref_ali2d")

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Reference stack             : %s\n"%(refim))	
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Maskfile                    : %s\n"%(maskfile))
	print_msg("Inner radius                : %i\n"%(first_ring))

	ima = EMData()
	ima.read_image(stack, 0)
	nx = ima.get_xsize()
	# default value for the last ring
	if last_ring == -1: last_ring = nx/2-2

	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Ring step                   : %i\n"%(rstep))
	print_msg("X search range              : %i\n"%(xrng))
	print_msg("Y search range              : %i\n"%(yrng))
	print_msg("Translational step          : %i\n"%(step))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("CTF correction              : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	print_msg("Random seed                 : %i\n\n"%(rand_seed))

	output = sys.stdout

	import user_functions
	user_func = user_functions.factory[user_func_name]

	if maskfile:
		import types
		if type(maskfile) is types.StringType:  mask = get_image(maskfile)
		else: mask = maskfile
	else: mask = model_circle(last_ring, nx, nx)
	#  references
	refi = []
	numref = EMUtil.get_image_count(refim)
	#  CTF stuff
	if CTF:
		ctf_params = ima.get_attr("ctf")
		data_had_ctf = ima.get_attr("ctf_applied")
		ctm = ctf_2(nx, ctf_params)
		lctf = len(ctm)
		ctf2 = [[[0.0]*lctf for k in xrange(2)] for j in xrange(numref)]

	# IMAGES ARE SQUARES! center is in SPIDER convention
	cnx = nx/2+1
	cny = cnx

	mode = "F"
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr = ringwe(numr, mode)
	# reference images
	params = []
	#read all data
	data = EMData.read_images(stack)
	nima = len(data)
	# prepare the reference
	ima.to_zero()
	for j in xrange(numref):
		temp = EMData()
		temp.read_image(refim, j)
		#  eve, odd, numer of even, number of images.  After frc, totav
		refi.append([temp, ima.copy(), 0])
	seed(rand_seed)
	a0 = -1.
	again = True
	Iter = 0

	ref_data = [mask, center, None, None]

	while Iter < max_iter and again:
		#again = False
		ringref = []
		#print "numref",numref
		for j in xrange(numref):
			refi[j][0].process_inplace("normalize.mask", {"mask":mask, "no_sigma":1})
			cimage = Util.Polar2Dm(refi[j][0], cnx, cny, numr, mode)
			Util.Frngs(cimage, numr)
			Applyws(cimage, numr, wr)
			ringref.append(cimage)
			# zero refi
			refi[j][0].to_zero()
			refi[j][1].to_zero()
			refi[j][2] = 0
			if CTF:
				for i in xrange(lctf): 
					ctf2[j][0][i] = 0.0
					ctf2[j][1][i] = 0.0
		assign = [[] for i in xrange(numref)]
		sx_sum = [0.0]*numref
		sy_sum = [0.0]*numref
		for im in xrange(nima):
			if CTF:
				ctf_params = data[im].get_attr("ctf")
				if data[im].get_attr("ctf_applied") == 0:
					st = Util.infomask(data[im], mask, False)
					data[im] -= st[0]
					from filter import filt_ctf
					data[im] = filt_ctf(data[im], ctf_params)
					data[im].set_attr('ctf_applied', 1)
			alpha, sx, sy, mirror, scale = get_params2D(data[im])
			alphai, sxi, syi, scalei = inverse_transform2(alpha, sx, sy)
			# normalize
			data[im].process_inplace("normalize.mask", {"mask":mask, "no_sigma":0})
			# align current image to the reference
			[angt, sxst, syst, mirrort, xiref, peakt] = Util.multiref_polar_ali_2d(data[im], 
				ringref, xrng, yrng, step, mode, numr, cnx+sxi, cny+syi)
			iref = int(xiref)
			# combine parameters and set them to the header, ignore previous angle and mirror
			[alphan, sxn, syn, mn] = combine_params2(0.0, -sxi, -syi, 0, angt, sxst, syst, int(mirrort))
			set_params2D(data[im], [alphan, sxn, syn, int(mn), scale])
			if mn == 0: sx_sum[iref] += sxn
			else: sx_sum[iref] -= sxn
			sy_sum[iref] += syn
			data[im].set_attr('assign', iref)
			# apply current parameters and add to the average
			temp = rot_shift2D(data[im], alphan, sxn, syn, mn)
			it = im%2
			Util.add_img(refi[iref][it], temp)
			if CTF:
				ctm = ctf_2(nx, ctf_params)
				for i in xrange(lctf):  ctf2[iref][it][i] += ctm[i]
			assign[iref].append(im)
			refi[iref][2] += 1
		del ringref
		if again:
			a1 = 0.0
			for j in xrange(numref):
				msg = "   group #%3d   number of particles = %7d\n"%(j, refi[j][2])
				print_msg(msg)
				if refi[j][2] < 4:
					#ERROR("One of the references vanished","mref_ali2d",1)
					#  if vanished, put a random image there
					assign[j] = []
					assign[j].append(randint(0, nima-1))
					refi[j][0] = data[assign[j][0]].copy()
				else:
					max_inter = 0  # switch off fine refi.
					br = 1.75
					#  the loop has to 
					for INter in xrange(max_inter+1):
						# Calculate averages at least ones, meaning even if no within group refinement was requested
						if CTF:
							for i in xrange(lctf):  ctm[i] = 1.0 / (ctf2[j][0][i] + 1.0/snr)
							from filter import filt_table
							av1 = filt_table(refi[j][0], ctm)
							for i in xrange(lctf):  ctm[i] = 1.0 / (ctf2[j][1][i] + 1.0/snr)
							av2 = filt_table(refi[j][1], ctm)
							frsc = fsc(av1, av2, 1.0, os.path.join(outdir,"drm_%03d_%04d.txt"%(Iter, j)))
							#Now the total average
							for i in xrange(lctf):  ctm[i] = 1.0 / (ctf2[j][0][i] + ctf2[j][1][i] + 1.0/snr)
							refi[j][0] = filt_table(Util.addn_img(refi[j][0], refi[j][1]), ctm)
						else:
							frsc = fsc(refi[j][0], refi[j][1], 1.0, os.path.join(outdir,"drm_%03d_%04d.txt"%(Iter, j)))
							Util.add_img(refi[j][0], refi[j][1])
							Util.mul_scalar(refi[j][0], 1.0/float(refi[j][2]))
							
						ref_data[2] = refi[j][0]
						ref_data[3] = frsc						
						refi[j][0], cs = user_func(ref_data)
						if center == -1:
							cs[0] = sx_sum[j]/len(assign[j])
							cs[1] = sy_sum[j]/len(assign[j])
							refi[j][0] = fshift(refi[j][0], -cs[0], -cs[1])
						for i in xrange(len(assign[j])):
							im = assign[j][i]
							alpha, sx, sy, mirror, scale =  get_params2D(data[im])
							alphan, sxn, syn, mirrorn = combine_params2(alpha, sx, sy, mirror, 0.0, -cs[0], -cs[1], 0)
							set_params2D(data[im], [alphan, sxn, syn, int(mirrorn), scale])
						# refine images within the group
						#  Do the refinement only if max_inter>0, but skip it for the last iteration.
						if INter < max_inter:
							fine_2D_refinement(data, br, mask, refi[j][0], j)
							#  Calculate updated average
							refi[j][0].to_zero()
							refi[j][1].to_zero()
							for i in xrange(len(assign[j])):
								im = assign[j][i]
								alpha, sx, sy, mirror, scale = get_params2D(data[im])
								# apply current parameters and add to the average
								temp = rot_shift2D(data[im], alpha, sx, sy, mn)
								it = im%2
								Util.add_img(refi[j][it], temp)

				# write the current average
				TMP = []
				for i_tmp in xrange(len(assign[j])):  TMP.append(float(assign[j][i_tmp]))
				TMP.sort()
				refi[j][0].set_attr_dict({'ave_n': refi[j][2], 'members': TMP })
				del TMP
				# replace the name of the stack with reference with the current one
				newrefim = os.path.join(outdir,"aqm%03d.hdf"%Iter)
				refi[j][0].write_image(newrefim, j)
				a1 += refi[j][0].cmp("dot", refi[j][0], {"negative":0, "mask":mask})
			Iter += 1
			msg = "ITERATION #%3d        criterion = %20.7e\n"%(Iter,a1)
			print_msg(msg)
			if a1 < a0:
				if auto_stop == True:	break
			else:	a0 = a1

	newrefim = os.path.join(outdir,"multi_ref.hdf")
	for j in xrange(numref):  refi[j][0].write_image(newrefim, j)
	if CTF:
		if data_had_ctf == 0:
			for im in xrange(nima): data[im].set_attr('ctf_applied', 0)
	from utilities import write_headers
	write_headers(stack, data, range(nima))
	print_end_msg("mref_ali2d")


def mref_ali2d_MPI(stack, refim, outdir, maskfile = None, ir=1, ou=-1, rs=1, xrng=0, yrng=0, step=1, center=1, maxit=10, CTF=False, snr=1.0, user_func_name="ref_ali2d", rand_seed=1000):
# 2D multi-reference alignment using rotational ccf in polar coordinates and quadratic interpolation

	from utilities      import   model_circle, combine_params2, inverse_transform2, drop_image, get_image, get_im
	from utilities      import   reduce_EMData_to_root, bcast_EMData_to_all, bcast_number_to_all
	from utilities      import   send_attr_dict
	from utilities	    import   center_2D
	from statistics     import   fsc_mask
	from alignment      import   Numrinit, ringwe, Applyws
	from fundamentals   import   rot_shift2D, fshift
	from utilities      import   get_params2D, set_params2D
	from random         import   seed, randint
	from morphology     import   ctf_2
	from filter         import   filt_btwl, filt_params
	from numpy          import   reshape, shape
	from utilities      import   print_msg, print_begin_msg, print_end_msg
	import os
	import sys
	from mpi 	  import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_recv, mpi_send
	from mpi 	  import MPI_SUM, MPI_FLOAT, MPI_INT, MPI_TAG_UB

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	# create the output directory, if it does not exist
	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "mref_ali2d_MPI ", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		print_begin_msg("mref_ali2d_MPI")
		os.mkdir(outdir)

	nima = EMUtil.get_image_count(stack)
	
	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)

	nima = EMUtil.get_image_count(stack)
	ima = EMData()
	ima.read_image(stack, image_start)

	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit)

	if max_iter == 0:
		max_iter  = 10
		auto_stop = True
	else:
		auto_stop = False

	if myid == main_node:
		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Reference stack             : %s\n"%(refim))	
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Maskfile                    : %s\n"%(maskfile))
		print_msg("Inner radius                : %i\n"%(first_ring))

	nx = ima.get_xsize()
	# default value for the last ring
	if last_ring == -1: last_ring=nx/2-2
	
	if myid == main_node:
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Ring step                   : %i\n"%(rstep))
		print_msg("X search range              : %f\n"%(xrng))
		print_msg("Y search range              : %f\n"%(yrng))
		print_msg("Translational step          : %f\n"%(step))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Random seed                 : %i\n\n"%(rand_seed))	

	import user_functions
	user_func = user_functions.factory[user_func_name]

	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  mask = get_image(maskfile)
		else: mask = maskfile
	else : mask = model_circle(last_ring, nx, nx)
	#  references, do them on all processors...
	refi = []
	numref = EMUtil.get_image_count(refim)
	#  CTF stuff
	if CTF:
		ctf_params = ima.get_attr("ctf")
		data_had_ctf = ima.get_attr("ctf_applied")
		ctm = ctf_2(nx, ctf_params)
		lctf = len(ctm)

	# IMAGES ARE SQUARES! center is in SPIDER convention
	cnx = nx/2+1
	cny = cnx

	mode = "F"
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr = ringwe(numr, mode)
	# reference images
	again = True
	params = []
	# prepare reference images on all nodes
	ima.to_zero()
	for j in xrange(numref):
		#  even, odd, numer of even, number of images.  After frc, totav
		refi.append([get_im(refim,j), ima.copy(), 0])
	#  for each node read its share of data
	data = EMData.read_images(stack, range(image_start, image_end))
	for im in xrange(image_start, image_end):
		data[im-image_start].set_attr('ID', im)
		if CTF:
			ctf_params = data[im-image_start].get_attr( "ctf" )
			if data[im-image_start].get_attr("ctf_applied") == 0:
				st = Util.infomask(data[im-image_start], mask, False)
				data[im-image_start] -= st[0]
				from filter import filt_ctf
				data[im-image_start] = filt_ctf(data[im-image_start], ctf_params)
				data[im-image_start].set_attr('ctf_applied', 1)
	if myid == main_node:  seed(rand_seed)

	a0 = -1.0
	again = True
	Iter = 0

	ref_data = [mask, center, None, None]

	while Iter < max_iter and again:
		ringref = []
		for j in xrange(numref):
			refi[j][0].process_inplace("normalize.mask", {"mask":mask, "no_sigma":1}) # normalize reference images to N(0,1)
			cimage = Util.Polar2Dm(refi[j][0] , cnx, cny, numr, mode)
			Util.Frngs(cimage, numr)
			Applyws(cimage, numr, wr)
			ringref.append(cimage)
			# zero refi
			refi[j][0].to_zero()
			refi[j][1].to_zero()
			refi[j][2] = 0
		if CTF: ctf2 = [[[0.0]*lctf for k in xrange(2)] for j in xrange(numref)]
		assign = [[] for i in xrange(numref)]
		# begin MPI section
		for im in xrange(image_start, image_end):
			alpha, sx, sy, mirror, scale = get_params2D(data[im-image_start])
			alphai, sxi, syi, scalei = inverse_transform2(alpha, sx, sy)
			# normalize
			data[im-image_start].process_inplace("normalize.mask", {"mask":mask, "no_sigma":0}) # subtract average under the mask
			# align current image to the reference
			[angt, sxst, syst, mirrort, xiref, peakt] = Util.multiref_polar_ali_2d(data[im-image_start], 
				ringref, xrng, yrng, step, mode, numr, cnx+sxi, cny+syi)
			iref = int(xiref)
			# combine parameters and set them to the header, ignore previous angle and mirror
			[alphan, sxn, syn, mn] = combine_params2(0.0, -sxi, -syi, 0, angt, sxst, syst, (int)(mirrort))
			set_params2D(data[im-image_start], [alphan, sxn, syn, int(mn), scale])
			data[im-image_start].set_attr('assign',iref)
			# apply current parameters and add to the average
			temp = rot_shift2D(data[im-image_start], alphan, sxn, syn, mn)
			it = im%2
			Util.add_img( refi[iref][it], temp)
			assign[iref].append(im)
			if CTF:
				#  I wonder whether params are still there....
				ctf_params = data[im-image_start].get_attr("ctf")
				ctm = ctf_2(nx, ctf_params)
				for i in xrange(lctf):  ctf2[iref][it][i] += ctm[i]
			#assign[im] = iref
			refi[iref][2] += 1.0
		del ringref
		# end MPI section, bring partial things together, calculate new reference images, broadcast them back
		if CTF:
		# bring ctf2 together on main node
			s = shape(ctf2)
			ctf2  = mpi_reduce(ctf2, 2*lctf*numref, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
			if myid == main_node: ctf2 = reshape(ctf2, s)
		for j in xrange(numref):
			reduce_EMData_to_root(refi[j][0], myid, main_node)
			reduce_EMData_to_root(refi[j][1], myid, main_node)
			refi[j][2] = mpi_reduce(refi[j][2], 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
			if(myid == main_node): refi[j][2] = int(refi[j][2][0])
		# gather assignements
		for j in xrange(numref):
			if myid == main_node:
				for n in xrange(number_of_proc):
					if n != main_node:
						ln =  mpi_recv(1, MPI_INT, n, MPI_TAG_UB, MPI_COMM_WORLD)
						lis = mpi_recv(ln[0], MPI_INT, n, MPI_TAG_UB, MPI_COMM_WORLD)
						for l in xrange(ln[0]): assign[j].append(int(lis[l]))
			else:
				mpi_send(len(assign[j]), 1,MPI_INT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)
				mpi_send(assign[j], len(assign[j]), MPI_INT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)

		if myid == main_node:
			# replace the name of the stack with reference with the current one
			refim = os.path.join(outdir,"aqm%03d.hdf"%Iter)
			a1 = 0.0
			ave_fsc = []
			for j in xrange(numref):
				if refi[j][2] < 4:
					#ERROR("One of the references vanished","mref_ali2d_MPI",1)
					#  if vanished, put a random image (only from main node!) there
					assign[j] = []
					assign[j].append( randint(image_start, image_end-1) - image_start )
					refi[j][0] = data[assign[j][0]].copy()
					#print 'ERROR', j
				else:
					if CTF:
						for i in xrange(lctf):  ctm[i] = 1.0 / (ctf2[j][0][i] + 1.0/snr)
						from filter import filt_table
						av1 = filt_table( refi[j][0], ctm)
						for i in xrange(lctf):  ctm[i] = 1.0 / (ctf2[j][1][i] + 1.0/snr)
						av2 = filt_table( refi[j][1], ctm)
						from statistics import fsc
						#frsc = fsc_mask(av1, av2, mask, 1.0, os.path.join(outdir,"drm%03d%04d"%(Iter, j)))
						frsc = fsc(av1, av2, 1.0, os.path.join(outdir,"drm%03d%04d.txt"%(Iter, j)))
						#Now the total average
						for i in xrange(lctf):  ctm[i] = 1.0 / (ctf2[j][0][i] + ctf2[j][1][i] + 1.0/snr)
						refi[j][0] = filt_table( Util.addn_img( refi[j][0], refi[j][1] ), ctm)
					else:
						#frsc = fsc_mask(refi[j][0], refi[j][1], mask, 1.0, os.path.join(outdir,"drm%03d%04d"%(Iter, j)))
						from statistics import fsc
						frsc = fsc(refi[j][0], refi[j][1], 1.0, os.path.join(outdir,"drm%03d%04d.txt"%(Iter,j)))
						Util.add_img( refi[j][0], refi[j][1] )
						Util.mul_scalar( refi[j][0], 1.0/float(refi[j][2]) )
				        	
					if ave_fsc == []:
						for i in xrange(len(frsc[1])): ave_fsc.append(frsc[1][i])
						c_fsc = 1
					else:
						for i in xrange(len(frsc[1])): ave_fsc[i] += frsc[1][i]
						c_fsc += 1
					#print 'OK', j, len(frsc[1]), frsc[1][0:5], ave_fsc[0:5]			


			#print 'sum', sum(ave_fsc)
			if sum(ave_fsc) != 0:		
				for i in xrange(len(ave_fsc)):
					ave_fsc[i] /= float(c_fsc)
					frsc[1][i]  = ave_fsc[i]
			
			for j in xrange(numref):
				ref_data[2]    = refi[j][0]
				ref_data[3]    = frsc
				refi[j][0], cs = user_func(ref_data)	

				# write the current average
				TMP = []
				for i_tmp in xrange(len(assign[j])): TMP.append(float(assign[j][i_tmp]))
				TMP.sort()
				refi[j][0].set_attr_dict({'ave_n': refi[j][2],  'members': TMP })
				del TMP
				refi[j][0].process_inplace("normalize.mask", {"mask":mask, "no_sigma":1})
				refi[j][0].write_image(refim, j)
				a1 += refi[j][0].cmp("dot", refi[j][0], {"negative":0, "mask":mask})
				# send refi[j][0]  back!

			Iter += 1
			msg = "ITERATION #%3d        criterion = %20.7e\n\n"%(Iter,a1)
			print_msg(msg)
			for j in xrange(numref):
				msg = "   group #%3d   number of particles = %7d\n"%(j, refi[j][2])
				print_msg(msg)
			
			if a1 < a0:
				if (auto_stop == True):	again = False
			else:	a0 = a1
		#again = mpi_bcast(again, 1, MPI_INT, main_node, MPI_COMM_WORLD)
		Iter  = bcast_number_to_all(Iter, main_node)
		if CTF:  del  ctf2
		if again:
			for j in xrange(numref):
				bcast_EMData_to_all(refi[j][0], myid, main_node)

	#  clean up
	del assign
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
	if CTF and data_had_ctf == 0:
		for im in xrange(len(data)): data[im].set_attr('ctf_applied', 0)
	par_str = ['xform.align2d', 'assign', 'ID']
	if myid == main_node:
		from utilities import file_type
		if(file_type(stack) == "bdb"):
			from utilities import recv_attr_dict_bdb
			recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
		else:
			from utilities import recv_attr_dict
			recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	else:           send_attr_dict(main_node, data, par_str, image_start, image_end)

	if myid == main_node:
		newrefim = os.path.join(outdir, "multi_ref.hdf")
		for j in xrange(numref):
			refi[j][0].write_image(newrefim, j)
		print_end_msg("mref_ali2d_MPI")


def ali2d_ra(stack, maskfile = None, ir = 1, ou = -1, rs = 1, maxit = 10, check_mirror = False, CTF = False, rand_seed = 1000):
# 2D rotational alignment using ccf in polar coordinates

	from utilities    import model_circle, compose_transform2, combine_params2, drop_image, get_im, get_arb_params, get_params2D, set_params2D
	from alignment    import Numrinit, ringwe, ang_n
	from statistics   import kmn, kmn_ctf
	from morphology   import ctf_2
	from statistics   import add_series
	from applications import transform2d
	from random       import random
	from utilities    import print_begin_msg, print_end_msg, print_msg

	print_begin_msg("ali2d_ra")

	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit); 

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Maskfile                    : %s\n"%(maskfile))
	print_msg("Inner radius                : %i\n"%(first_ring))

	temp = EMData()
	temp.read_image(stack, 0)
	nx = temp.get_xsize()
	ny = nx
	# default value for the last ring
	if (last_ring == -1): last_ring=nx//2-2

	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Ring step                   : %i\n"%(rstep))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("Consider Mirror             : %s\n"%(check_mirror))	
	print_msg("CTF correction              : %s\n"%(CTF))
	print_msg("Random seed                 : %i\n\n"%(rand_seed))
	
	# consider mirror
	# somebody dedicated could write a version with option "H" for half rings that would work for ACF functions.
	mode = "F"

	nima = EMUtil.get_image_count(stack)

	# precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr   = ringwe(numr, mode)
	lnumr = numr[len(numr)-1]
	# prepare 2-D mask for normalization
	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask2D = get_im(maskfile)
		else: mask2D = maskfile
	else : mask2D = model_circle(last_ring, nx, nx)
	if (first_ring > 0):
		tave = model_circle(first_ring-1, nx, ny)
		mask2D -= tave

	# read images and resample them into polar coordinates
	data = []
	#  center is in SPIDER convention
	cnx = int(nx/2) + 1
	cny = int(ny/2) + 1

	if(CTF):
		# for rotational alignment with CTF correction, we have more complicated strategies.
		#   We need two series of images in polar coordinates: multiplied by the ctf, and 'ctf-corrected', i.e., also divided by ctf^2+1/snr
		#   Thus, is alignment, the reference is always computed from the second series, while the first is used for alignment 
		parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
		#                        0                1              2            3              4             5                6
		for im in xrange(nima):
			if(im>0):
				temp = EMData()
				temp.read_image(stack, im)
			ctf_params = temp.get_attr( "ctf" )
			if(im == 0):  data_had_ctf = temp.get_attr('ctf_applied')
			ctf = ctf_2(nx, ctf_params)
			if(im == 0):
				lctf = len(ctf)
				ctf2 = [0.0]*lctf
			for i in xrange(lctf):  ctf2[i] += ctf[i]
		del ctf
		ref_data = []
		for im in xrange(nima):
			temp = EMData()
			temp.read_image(stack, im)
			st = Util.infomask(temp, mask2D, False)
			ctf_params = temp.get_attr( "ctf" )
			temp -= st[0]
			if(temp.get_attr("ctf_applied") == 0):
				from filter import filt_ctf
				temp = filt_ctf(temp, ctf_params)
				temp.set_attr('ctf_applied', 1)
			from filter       import filt_table
			refc = filt_table(temp, ctf2)
			
			alpha_original, sx, sy, miri, scale = get_params2D(temp)
			#tempg = prepg(temp, kb)
			#cimage = Util.Polar2Dmi(tempg, cnx+sx, cny+sy, numr, mode, kb)
			alphan, sxn, syn, mir = combine_params2(0, -sx, -sy, 0, -alpha_original, 0,0,0)
			cimage = Util.Polar2Dm(temp, cnx+sxn, cny+syn, numr, mode)
			Util.Frngs(cimage, numr)
			data.append(cimage)

			#  Here alpha is position of the peak (i.e., starts from 1), just put as a place holder, will be determined in kmn
			data[im].set_attr_dict({'alpha':1.0, 'alpha_original':alpha_original, 'sx':sx, 'sy':sy, 'mirror': 0})
			cimage = Util.Polar2Dm(refc, cnx+sxn, cny+syn, numr, mode)
			Util.Frngs(cimage, numr)
			#  We do not need any attributes for ref_data, as they are going to be taken from data
			ref_data.append(cimage)

		del temp
		del refc
		del mask2D
		del ctf2
		kmn_ctf(data, ref_data, numr, wr, check_mirror, max_iter, rand_seed)
	else:
		for im in xrange(nima):
			if (im>0):
				temp = EMData()
				temp.read_image(stack, im)
			alpha_original, sx, sy, miri, scale = get_params2D(temp)
			temp.process_inplace("normalize.mask", {"mask":mask2D, "no_sigma":0})
			#tempg = prepg(temp, kb)
			#cimage = Util.Polar2Dmi(tempg, cnx+sx, cny+sy, numr, mode, kb)
			alphan, sxn, syn, mir = combine_params2(0.0, -sx, -sy, 0, -alpha_original, 0.0, 0.0, 0)
			cimage = Util.Polar2Dm(temp, cnx+sxn, cny+syn, numr, mode)
			Util.Frngs(cimage, numr)
			data.append(cimage)
			#  Here alpha is position of the peak (i.e., starts from 1), just put as a place holder, will be determined in kmn
			data[im].set_attr_dict({'alpha':1.0, 'alpha_original':alpha_original, 'sx':sx, 'sy':sy, 'mirror':0})
		del temp
		#del tempg
		del mask2D
		kmn(data, numr, wr, check_mirror, max_iter, rand_seed)
	#  write out the alignment parameters to headers
	from utilities import write_header, file_type
	temp = EMData()
	for im in xrange(nima):
		alpha_original   = data[im].get_attr('alpha_original')
		alpha = data[im].get_attr('alpha')
		sx    =  data[im].get_attr('sx')
		sy    =  data[im].get_attr('sy')
		mirror =  data[im].get_attr('mirror')
		alpha = ang_n(alpha, mode, lnumr)
		#  here the original angle is irrelevant, used only to determine proper shifts
		alpha_original_n, sxn, syn, mir = combine_params2(0, -sx, -sy, 0, -alpha_original, 0,0,0)
		alphan, sxn, syn, mir           = combine_params2(0, -sxn, -syn, 0, alpha, 0,0, mirror)
		temp.read_image(stack, im, True)
		if(CTF and data_had_ctf == 0):   temp.set_attr('ctf_applied', 0)
		set_params2D(temp, [alphan, sxn, syn, mir, 1.0])
		write_header(stack, temp, im)
		#temp.write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)
	print_end_msg("ali2d_ra")

def ali2d_rag(stack, maskfile = None, ir = 1, ou = -1, rs = 1, maxit = 10, check_mirror = False, CTF = False, rand_seed = 1000):
# 2D rotational alignment using ccf in polar coordinates and gridding-based interpolation

	from utilities    import model_circle, compose_transform2, combine_params2, drop_image, get_im, get_arb_params
	from alignment    import Numrinit, ringwe, ang_n
	from statistics   import kmn_g, kmn_ctf
	from morphology   import ctf_2
	from statistics   import add_series
	from applications import transform2d
	from random       import random
	from utilities    import print_begin_msg, print_end_msg, print_msg
	from fundamentals import prepi

	print_begin_msg("ali2d_rag")

	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit); 

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Maskfile                    : %s\n"%(maskfile))
	print_msg("Inner radius                : %i\n"%(first_ring))

	temp = EMData()
	temp.read_image(stack, 0)
	nx = temp.get_xsize()
	ny = nx
	# default value for the last ring
	if (last_ring == -1): last_ring=nx//2-2

	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Ring step                   : %i\n"%(rstep))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("Consider Mirror             : %s\n"%(check_mirror))	
	print_msg("CTF correction              : %s\n"%(CTF))
	print_msg("Random seed                 : %i\n\n"%(rand_seed))
	
	# consider mirror
	# somebody dedicated could write a version with option "H" for half rings that would work for ACF functions.
	mode = "F"

	nima = EMUtil.get_image_count(stack)

	# precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr   = ringwe(numr, mode)
	maxrin = numr[len(numr)-1]
	# prepare 2-D mask for normalization
	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask2D = get_im(maskfile)
		else: mask2D = maskfile
	else : mask2D = model_circle(last_ring, nx, nx)
	if (first_ring > 0):
		tave = model_circle(first_ring-1, nx, ny)
		mask2D -= tave

	# read images and resample them into polar coordinates
	data = []
	#  center is in SPIDER convention
	cnx = int(nx/2) + 1
	cny = int(ny/2) + 1

	if(CTF):
		# for rotational alignment with CTF correction, we have more complicated strategies.
		#   We need two series of images in polar coordinates: multiplied by the ctf, and 'ctf-corrected', i.e., also divided by ctf^2+1/snr
		#   Thus, is alignment, the reference is always computed from the second series, while the first is used for alignment 
		parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
		#                        0                1              2            3              4             5                6
		for im in xrange(nima):
			if(im>0):
				temp = EMData()
				temp.read_image(stack, im)
			ctf_params = temp.get_attr( "ctf" )
			if(im == 0):  data_had_ctf = temp.get_attr('ctf_applied')
			ctf = ctf_2(nx, ctf_params)
			if(im == 0):
				lctf = len(ctf)
				ctf2 = [0.0]*lctf
			for i in xrange(lctf):  ctf2[i] += ctf[i]
		del ctf
		ref_data = []
		for im in xrange(nima):
			temp = EMData()
			temp.read_image(stack, im)
			st = Util.infomask(temp, mask2D, False)
					
			temp -= st[0]
			if(temp.get_attr("ctf_applied") == 0):
				from filter import filt_ctf
				ctf_params = temp.get_attr( "ctf" )
				temp = filt_ctf(temp, ctf_params)
				temp.set_attr('ctf_applied', 1)
			from filter       import filt_table
			refc = filt_table(temp, ctf2)
			
			alpha_original = temp.get_attr('alpha')
			sx =  temp.get_attr('sx')
			sy =  temp.get_attr('sy')
			miri = temp.get_attr('mirror')
			alphan, sxn, syn, mir = combine_params2(0, -sx, -sy, 0, -alpha_original, 0,0,0)
			tempg, kb = prepi(temp)
			cimage = Util.Polar2Dmi(tempg, cnx+sxn, cny+syn, numr, mode, kb)
			Util.Frngs(cimage, numr)
			data.append(cimage)

			#  Here alpha is position of the peak (i.e., starts from 1), just put as a place holder, will be determined in kmn
			data[im].set_attr_dict({'alpha':1.0, 'alpha_original':alpha_original, 'sx':sx, 'sy':sy, 'mirror': 0})
			cimage = Util.Polar2Dmi(refc, cnx+sxn, cny+syn, numr, mode, kb)
			Util.Frngs(cimage, numr)
			#  We do not need any attributes for ref_data, as they are going to be taken from data
			ref_data.append(cimage)

		del temp
		del refc
		del mask2D
		del ctf2
		kmn_ctf(data, ref_data, numr, wr, check_mirror, max_iter, rand_seed)
	else:
		for im in xrange(nima):
			if (im>0):
				temp = EMData()
				temp.read_image(stack, im)
			alpha_original = temp.get_attr('alpha')
			sx =  temp.get_attr('sx')
			sy =  temp.get_attr('sy')
			miri = temp.get_attr('mirror')
			temp.process_inplace("normalize.mask", {"mask":mask2D, "no_sigma":0})
			alphan, sxn, syn, mir = combine_params2(0.0, -sx, -sy, 0, -alpha_original, 0.0, 0.0, 0)
			tempg, kb = prepi(temp)
			cimage = Util.Polar2Dmi(tempg, cnx+sxn, cny+syn, numr, mode, kb)
			Util.Frngs(cimage, numr)
			data.append(cimage)
			#  Here alpha is position of the peak (i.e., starts from 1), just put as a place holder, will be determined in kmn
			data[im].set_attr_dict({'alpha':1.0, 'alpha_original':alpha_original, 'sx':sx, 'sy':sy, 'mirror':0})
		del temp
		del tempg
		del mask2D
		kmn_g(data, numr, wr, stack, check_mirror, max_iter, rand_seed)
	#  write out the alignment parameters to headers
	from utilities import write_header, file_type
	ext = file_type(stack)
	if(ext == "bdb"):
		from EMAN2db import EMAN2DB
		DB = EMAN2DB()
		DB = EMAN2DB.open_db(ipath)
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
		if(CTF and data_had_ctf == 0):   temp.set_attr('ctf_applied', 0)
		temp.set_attr_dict({'alpha':alphan, 'sx':sxn, 'sy':syn, 'mirror': mir})
		write_header(stack, temp, im)
		#temp.write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)
	if(ext == "bdb"):
		DB.close_dict(ipath)
	print_end_msg("ali2d_rag")
	
def ali2d_rac(stack, maskfile = None, ir = 1, ou = -1, rs = 1, nclass = 2, maxit = 10, maxin = 10, check_mirror = False, rand_seed = 1000, MPI=False):
# 2D rotational classification and alignment using ccf in polar coords, no CTF
	test = True
	if MPI:
		if test:
			from development import ali2d_rac_MPI
			ali2d_rac_MPI(stack, maskfile, ir, ou, rs, nclass, maxit, maxin, check_mirror, rand_seed)
			return
		else:
			print 'ali2d_rac: no mpi version'
			return

	from utilities    import model_circle, combine_params2, drop_image
	from alignment    import Numrinit, ringwe, ang_n
	from statistics   import kmnr, kmn, add_series_class
	from random       import seed, randint

	from utilities    import info, ttime, print_list_format
	import time

	seed(rand_seed)
	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); kc=int(nclass); max_iter=int(maxit); max_internal=int(maxin);

	# somebody dedicated could write a version with option "H" for half rings that would work for ACF functions.
	mode = 'F'
	nima = EMUtil.get_image_count(stack)

	from utilities import print_begin_msg, print_end_msg, print_msg
	print_begin_msg('ali2d_rac')
	print_msg("Input stack                 : %s\n"%(stack))
	print_msg('Number of images            : %i\n' % nima)
	print_msg("Maskfile                    : %s\n"%(maskfile))
	print_msg("Inner radius                : %i\n"%(first_ring))
	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Ring step                   : %i\n"%(rstep))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg('Maximum intern iteration    : %i\n' % max_internal)
	print_msg("Consider Mirror             : %s\n"%(check_mirror))
	print_msg('Number of classes           : %i\n' % kc)
	print_msg("Random seed                 : %i\n"%(rand_seed))
	print_msg('Ouput stack                 : %s\n\n' % stack)
	t_start = time.time()

	# create the output directory, if it does not exist
	#if os.path.exists(outdir) is False: os.mkdir(outdir)
	temp = EMData()
	temp.read_image(stack, 0)
	nx = temp.get_xsize()
	ny = nx
	
	# default value for the last ring
	if(last_ring==-1): last_ring=nx//2-2
	
	# precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	# norm to criterion
	norm_rsd = 0
	for n in xrange(1, len(numr), 3): norm_rsd += numr[n]
		
	wr    = ringwe(numr ,mode)
	lnumr = numr[len(numr)-1]
	# prepare 2-D ask for normalization
	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask2D = get_image(maskfile)
		else: mask2D = maskfile
	else : mask2D = model_circle(last_ring,nx,nx)
	if(first_ring > 0):
		tave = model_circle(first_ring-1, nx, ny)
		mask2D -= tave


	# read images and resample them into polar coordinates
	data = []
	org  = []
	#  center 
	cnx = int(nx/2)
	cny = int(ny/2)
	
	for im in xrange(nima):
		if(im>0):
			temp = EMData()
			temp.read_image(stack, im)
		sx =  temp.get_attr('sx')
		sy =  temp.get_attr('sy')
		alpha_original = temp.get_attr('alpha')
		miri = temp.get_attr('mirror')
		[mean, sigma, qn, qm] = Util.infomask(temp, mask2D, True)
		temp = (temp - mean)/sigma
		alpha_original_n,sxn,syn,mir = combine_params2(0, -sx, -sy, 0, -alpha_original,0,0,0)
		cimage = Util.Polar2Dm(temp, cnx+sxn, cny+syn, numr, mode)
		Util.Frngs(cimage, numr)
		data.append(cimage)
		#  Here alpha is postion of the pick (i.e., starts from 1)
		data[im].set_attr_dict({'alpha':1.0, 'alpha_original':alpha_original, 'sx':sx, 'sy':sy, 'mirror': 0})

	del temp
	del mask2D
	
	# initial classification
	assign = [0]*nima
	nclass = [0]*kc
	init = "Random"
	#init = "Even"
	#init = "Batch"
	if(init == "Random"):
		again = True
		while(again):
			for i in xrange (nima):
				assign[i] = randint(0,kc-1)
				nclass[assign[i]] += 1
			# make sure than no class is zero
			again = False
			for k in xrange(kc):
				if(nclass[k] == 0): again = True
	elif(init == "Even"):
		for i in xrange (nima):
			assign[i] = i%kc
			nclass[assign[i]] += 1
	elif(init == "Batch"):
		kt = (nima+kc)//kc
		for i in xrange (nima):
			assign[i] = i//kt
			nclass[assign[i]] += 1

	
	print  '%-20s'%("Initialization")
	print_list_format(nclass)

	again = True
	it = -1

	t= time.time()
	while(again and it < (max_iter - 1)):
		it += 1
		# averages
		tave = []
		for k in xrange(kc):
			temp = kmnr(data, assign, nclass[k], k, numr, wr, check_mirror, max_internal, rand_seed+(it * 10))
			temp /= nclass[k]
			tave.append(temp)
					
		#  assign objects to averages, angles are irrelevant
		for k in xrange(kc):  nclass[k] = 0
		again = False

		# compute norm to distance
		norm = []
		for k in xrange(kc):
			retval = Util.Crosrng_ew(tave[k], tave[k], numr, wr, 0)
			sq1    = Util.ener(tave[k], numr)
			norm.append(2*sq1 / float(retval['qn']))
			
		Je_rsd = 0
		for im in xrange (nima):
			dmin = 1.0e20
			g    = -1
			for k in xrange(kc):
				
				if (check_mirror):
					retvals = Util.Crosrng_ew(tave[k], data[im], numr, wr, 0)
					qn  = retvals["qn"]
					retvals = Util.Crosrng_ew(tave[k], data[im], numr, wr, 1)
					qm  = retvals["qn"]
					qn = max(qn,qm)
				else:
					retvals = Util.Crosrng_ew(tave[k], data[im], numr, wr, 0)
					qn      = retvals["qn"]

				q1 = Util.ener(tave[k], numr)
				q2 = Util.ener(data[im], numr)
				qn = q1 + q2 - (qn * norm[k])
									
				if(qn < dmin):
					dmin = qn
					g = k

			if( g <0 ):
				print  "  Error in assignment of objects to averages "
				break
			if(assign[im] != g):
				again = True
				assign[im] = g

			nclass[g] += 1
			Je_rsd    += dmin / float(norm_rsd)

		#print 'before remove class'
		#print_list_format(nclass)
			
		# eliminate empty classes
		rep = True
		while(kc > 1 and rep):
			rep = False
			for k in xrange(kc):
				if(nclass[k] == 0):
					del nclass[k]
					del tave[k]
					#  fix the assign list
					for i in xrange(nima):
						if(assign[i] > k):   assign[i] -= 1
						elif assign[i] == k: assign[i] = randint(0, kc - 2)
					kc -= 1
					print  "  Empty class was eliminated, new number of classes =",kc
					if(kc == 1): again = False
					rep = True
					break

		print '%-20s %5d'%("ITERATION #",it+1)
		print_list_format(nclass)

		print_msg('> iteration %d      criterion %5.3e\n' % (it+1, Je_rsd))
		print '> iteration %d      criterion %5.3e\n' % (it+1, Je_rsd)

	print_msg('\n')
	for k in xrange(kc): print_msg('Cls[%3d]: %d\n' % (k, nclass[k]))
	print_msg('\n')

	print 'time ite:', time.time() - t
	print '%30s %5d %10s'%("Numbers of objects in ",kc," classes:")
	print_list_format(nclass)
	#align class averages and transfer parameters to individual images
	for k in xrange(kc):
		tave[k].set_attr_dict({'alpha':1.0, 'mirror':0})
	
	kmn(tave, numr, wr, check_mirror, max_iter)

	#print ttime()

	talpha = [0]*kc
	tmir   = [0]*kc
	for k in xrange(kc):
		talpha[k] = ang_n(tave[k].get_attr('alpha'), mode, lnumr)
		tmir[k]   = tave[k].get_attr('mirror')

	del tave
	#  write out the alignment parameters to headers
	del temp
	from utilities import write_header, file_type
	ext = file_type(stack)
	if(ext == "bdb"):
		from EMAN2db import EMAN2DB
		DB = EMAN2DB()
		DB = EMAN2DB.open_db(ipath)
	temp = EMData()
	for im in xrange(nima):
				
		#  First combine with angle of the average
		alpha = ang_n(data[im].get_attr('alpha'), mode, lnumr)
		mirror =  data[im].get_attr('mirror')
		alpha, k,it, mirror = combine_params2(alpha, 0,0, mirror, talpha[assign[im]], 0,0, tmir[assign[im]])

		#  Second combine with given alignment
		alpha_original   =  data[im].get_attr('alpha_original')
		sx    =  data[im].get_attr('sx')
		sy    =  data[im].get_attr('sy')
		alpha_original_n, sxn, syn, mir = combine_params2(0, -sx, -sy, 0, -alpha_original,0,0,0)
		alphan, sxn, syn, mir           = combine_params2(0, -sxn, -syn, 0, alpha, 0,0,mirror)
		temp.read_image(stack, im, True)
		temp.set_attr_dict({'alpha':alphan, 'sx':sxn, 'sy':syn, 'mirror': mir, 'nclass':kc, 'assign':assign[im]})
		
		#if(data_had_ctf == 0):   temp.set_attr('ctf_applied', 0)
		write_header(stack, temp, im)
		#temp.write_image(stack ,im, EMUtil.ImageType.IMAGE_HDF, True)
	if(ext == "bdb"):
		DB.close_dict(ipath)
	del temp
	del data
	#for im in xrange(nima): data[im].write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)
	#transform2d(stack, "dala.hdf")
	#ave, var, nclass = add_series_class("dala.hdf")
	#for k in xrange(kc):
	#		ave[k].set_attr_dict({'Class_average':1, 'nobjects':nclass[k]})
	#		var[k].set_attr_dict({'Class_variance':1, 'nobjects':nclass[k]})
	#		ave[k].write_image("class_ave.hdf", k)
	#		var[k].write_image("class_var.hdf", k)

	print_msg("\nTime: %f s\n" % (time.time() - t_start))
	print_end_msg("ali2d_rac")

def ali2d_cross_res(stack, outdir, maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", center=1, maxit=0, CTF=False, snr=1.0, user_func_name="ref_ali2d"):
 	"""
		Split data into odd and even data sets and align them seperately
		Cross resolution alignment
	"""
	import os
	from utilities 		import model_circle, combine_params2, drop_image
	from utilities          import get_input_from_string, get_image, get_arb_params, set_arb_params
	from fundamentals 	import rot_shift2D
	from statistics 	      import add_oe_series, ave_series_ctf, ave_series, fsc_mask
	from alignment 		import Numrinit, ringwe, ali2d_single_iter, align2d
	from filter 		import filt_table, filt_ctf
	from morphology     import ctf_2

	from utilities import print_begin_msg, print_end_msg, print_msg
	import	types
	
	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ali2d_cross_res", 1)
	os.mkdir(outdir)
		
	print_begin_msg("ali2d_cross_res")

	import user_functions
	user_func = user_functions.factory[user_func_name]

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)

	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit);

	if max_iter == 0:
		max_iter  = 10
		auto_stop = True
	else:
		auto_stop = False

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Inner radius                : %i\n"%(first_ring))

	nima = EMUtil.get_image_count(stack)
	ima = EMData()
	ima.read_image(stack, 0)
	nx = ima.get_xsize()
	# default value for the last ring
	if (last_ring == -1):  last_ring = nx//2-2

	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Ring step                   : %i\n"%(rstep))
	print_msg("X search range              : %s\n"%(xrng))
	print_msg("Y search range              : %s\n"%(yrng))
	print_msg("Translational step          : %s\n"%(step))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("CTF correction              : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	if auto_stop:  print_msg("Stop iteration with         : criterion\n")
	else:           print_msg("Stop iteration with         : maxit\n")


	if maskfile:
		if(type(maskfile) is types.StringType):
			print_msg("Maskfile                    : %s\n\n"%(maskfile))
			mask=get_image(maskfile)
		else:
			print_msg("Maskfile                    : user provided in-core mask\n\n")
			mask = maskfile
	else :
		print_msg("Maskfile                    : default, a circle with radius %i\n\n"%(last_ring))
		mask = model_circle(last_ring, nx, nx)

	#  ODD-EVEN is controlled by setting NG to 2
	NG = 2
	cnx = int(nx/2)+1
 	cny = cnx
 	mode = "F"
	if(CTF):
		ctf_params = ima.get_attr( "ctf" )
		data_had_ctf = ima.get_attr( "ctf_applied" )
		ctm = ctf_2(nx, ctf_params)
		lctf = len(ctm)
		ctf2 = [[[0.0]*lctf for j in xrange(2)] for i in xrange(NG)]
		ctfb2 = [[0.0]*lctf for i in xrange(NG)]
	del ima
	all_data = EMData.read_images(stack)
	for im in xrange(nima):
		all_data[im].set_attr('ID', im)
		k = im%NG
		if(CTF):
			ctf_params = all_data[im].get_attr( "ctf" )
			ctm = ctf_2(nx, ctf_params)

			kl = (im//2)%NG  # not sure it will work for NG>2
			for i in xrange(lctf):
				ctf2[k][kl][i] += ctm[i]
	
  			if(all_data[im].get_attr("ctf_applied") == 0):
				st = Util.infomask(all_data[im], mask, False)
				all_data[im] -= st[0]
				all_data[im] = filt_ctf(all_data[im], ctf_params)
				all_data[im].set_attr('ctf_applied', 1)

	#  create to lists of images in groups.
	data = [[] for i in xrange(NG)]
	for im in xrange(nima):
		k = im%NG
		data[k].append(all_data[im])

	if(CTF):
		ctf_tot = [0.0]*lctf
		for i in xrange(lctf):
			for k in xrange(NG):
				ctf_tot[i] += ctf2[k][0][i] + ctf2[k][1][i]
			ctf_tot[i] = 1.0/(ctf_tot[i] + 1.0/snr)
		for k in xrange(NG):
			for i in xrange(lctf):
				ctfb2[k][i] = 1.0/(ctf2[k][0][i] + ctf2[k][1][i] + 1.0/snr)
				for kl in xrange(2):
					ctf2[k][kl][i] = 1.0/(ctf2[k][kl][i] + 1.0/snr)
 	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
 	wr = ringwe(numr, mode)
	tavg = [None]*NG
	for k in xrange(NG):
		av1, av2 = add_oe_series(data[k])
		Util.add_img(av1, av2)
		if(CTF):  tavg[k] = filt_table(av1, ctfb2[k])
		else:	  tavg[k] = av1/len(data[k])
		drop_image(tavg[k],os.path.join(outdir, "aqc_%03d_%03d.hdf"%(k, 0)))
	fscross = fsc_mask(tavg[0], tavg[1], mask, 1.0, os.path.join(outdir, "drcross_%03d"%(0)))

	if(CTF): total_ave = ave_series_ctf(all_data, ctf_tot)
	else:    total_ave = ave_series(all_data)
	drop_image(total_ave, os.path.join(outdir, "total_ave_%03d.hdf"%(0)))
	a0 = total_ave.cmp("dot", total_ave, dict(negative = 0, mask = mask))
	msg = "Initial criterion = %-20.7e\n"%(a0)
	print_msg(msg)
	params = ["alpha", "sx", "sy", "mirror"]
	# initialize data for the reference preparation function
	#  mask can be modified in user_function
	ref_data = [mask, center, None, None]
	cs=[[0.0,0.0]]*NG
	total_iter = 0
	for N_step in xrange(len(xrng)):
		msg = "\nX range = %5.2f   Y range = %5.2f   Step = %5.2f\n"%(xrng[N_step], yrng[N_step], step[N_step])
		print_msg(msg)
		for Iter in xrange(max_iter):
			total_iter += 1
			frsc = []
			ktavg = [None]*NG
			for k in xrange(NG):
				ali2d_single_iter(data[k], numr, wr, cs[k], tavg[k], cnx, cny, xrng[N_step], yrng[N_step], step[N_step], mode)
				av1, av2 = add_oe_series(data[k])
				if(CTF):
					tavg[k] = filt_table(Util.addn_img(av1, av2), ctfb2[k])
					av1    = filt_table(av1, ctf2[k][0])
					av2    = filt_table(av2, ctf2[k][1])
				else:
					tavg[k] = (av1+av2)/len(data[k])
				drop_image(tavg[k], os.path.join(outdir, "aqc_%03d_%03d.hdf"%(k, total_iter)))

				frsc.append(fsc_mask(av1, av2, ref_data[0], 1.0, os.path.join(outdir, "resolution_%03d_%03d"%(k, total_iter))))
				#  prepare averages for alignment
				kref_data = [mask, 0, tavg[k], frsc[k]]
				#  call the user-supplied function to prepare reference image, i.e., filter it, but do not center!
				ktavg[k], cs[k] = user_func( kref_data )
				del kref_data
			#  This should be done only for estimation of resolution, nothing else!
			alpha, sx, sy, mirror, peak = align2d(ktavg[0], ktavg[1], xrng[0], yrng[0], step=0.25, first_ring = first_ring, last_ring = last_ring, rstep=1, mode = mode)
			#  apply parameters to the original average
			favg2 = rot_shift2D(tavg[0], alpha, sx, sy, mn, interpolation_method="gridding")
			fscross = fsc_mask(favg2, tavg[1], ref_data[0], 1.0, os.path.join(outdir, "drcross_%03d"%(total_iter)))
			del favg2
			# Here one may want to apply rot-shift of the first average to all images in its group
			for im in xrange(len(data[0])):
				ps = get_arb_params(data[0][im], params)
				an,sxn,syn,mrn = combine_params2(ps[0], ps[1], ps[2], ps[3], alpha, sx, sy, mirror)
				set_arb_params(data[0][im], [an,sxn,syn,mrn], params)
			k = 0
			av1, av2 = add_oe_series(data[k])
			if(CTF):
				tavg[k] = filt_table(Util.addn_img(av1, av2), ctfb2[k])
			else:
				tavg[k] = (av1+av2)/len(data[k])
			#  Here we have to change fsc values.  The reason is that we have crossresolution, so snr can be calculated directly,
			#        while in user function the fit to fsc is done under assumption that is was calculated by splitting the dataset, so it has a factor of 2
			for i in xrange(len(fscross[1])):   fscross[1][i] = fscross[1][i]/(2.0-fscross[1][i])
			for k in xrange(NG):
				#  Apply the same filtration to all averages
				ref_data[2] = tavg[k]
				ref_data[3] = fscross
				#  call the user-supplied function to prepare reference image, i.e., filter and center!
				tavg[k], cs[k] = user_func( ref_data )
				drop_image(tavg[k], os.path.join(outdir, "aqf_%03d_%03d.hdf"%(k, total_iter)))

			if(CTF): total_ave = ave_series_ctf(all_data, ctf_tot)
			else:    total_ave = ave_series(all_data)
			drop_image(total_ave, os.path.join(outdir, "total_ave_%03d.hdf"%(total_iter)))
			# a0 should increase; stop algorithm when it decreases.
			a1 = total_ave.cmp("dot", total_ave, dict(negative = 0, mask = ref_data[0]))
			msg = "ITERATION #%3d	     criterion = %20.7e\n"%(total_iter,a1)
			print_msg(msg)
			if(a1 < a0):
				if (auto_stop == True): break
			else:	a0 = a1
	# write out headers
	if(CTF and data_had_ctf == 0):
		for k in xrange(NG):
			for im in xrange(len(data[k])):
				data[k][im].set_attr('ctf_applied', 0)
	from utilities import write_header, file_type
	ext = file_type(stack)
	if(ext == "bdb"):
		from EMAN2db import EMAN2DB
		DB = EMAN2DB()
		DB = EMAN2DB.open_db(ipath)
	for im in xrange(nima):
		k=im%NG
		imm = im//NG
		write_header(stack, data[k][imm], im)
		#data[k][imm].write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)
	if(ext == "bdb"):
		DB.close_dict(ipath)
	print_end_msg("ali2d_cross_res")

def ali3d_a(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta="10 6 4 4", an="-1", 
	    center = 1.0, maxit = 5, CTF = False, ref_a = "S", sym="c1", user_func_name="ref_ali3d"):

	from utilities      import model_circle, drop_image
	from utilities      import get_image, get_input_from_string
	from utilities      import get_arb_params, set_arb_params
	from filter	    import filt_params, filt_btwl, filt_from_fsc, filt_table, fit_tanh, filt_tanl
	from alignment	    import proj_ali_incore, proj_ali_incore_local
	from statistics     import fsc_mask
	from fundamentals   import fft
	from reconstruction import recons3d_nn_SSNR
	import os
	import types
	from utilities      import print_begin_msg, print_end_msg, print_msg
	print_begin_msg("ali3d_a")

	import user_functions
	user_func = user_functions.factory[user_func_name]

	if CTF : from reconstruction import recons3d_4nn_ctf
	else   : from reconstruction import recons3d_4nn

	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', " ali3d_a", 1)
	os.mkdir(outdir)

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min( len(xrng), len(yrng), len(step), len(delta) )
	if (an == "-1"):
		an = []
		for i in xrange(lstp):   an.append(-1)
	else:
		an = get_input_from_string(an)
	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Reference volume            : %s\n"%(ref_vol))	
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Maskfile                    : %s\n"%(maskfile))
	print_msg("Inner radius                : %i\n"%(first_ring))

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	if (last_ring == -1):	last_ring = nx//2 - 2

	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Ring step                   : %i\n"%(rstep))
	print_msg("X search range              : %s\n"%(xrng))
	print_msg("Y search range              : %s\n"%(yrng))
	print_msg("Translational step          : %s\n"%(step))
	print_msg("Angular step                : %s\n"%(delta))
	print_msg("Angular search range        : %s\n"%(an))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("CTF correction              : %s\n"%(CTF))
	print_msg("Reference projection method : %s\n"%(ref_a))
	print_msg("Symmetry group              : %s\n\n"%(sym))

	if (maskfile) :
		if (type(maskfile) is types.StringType): mask3D = get_image(maskfile)
		else                                  : mask3D = maskfile
	else          :   mask3D = model_circle(last_ring, nx, nx, nx)
	mask = model_circle(last_ring, nx, nx)
	#drop_image(vol, os.path.join(outdir,"ref_vol00.hdf"))

	del ima
	data = EMData.read_images(stack)
	nima = len(data)
	for im in xrange(nima):
		data[im].set_attr('ID', im)
		if(CTF):
			if(data[im].get_attr_default('ctf_applied', 2) > 0):
				ERROR("data cannot be ctf-applied","ali2d_a",1)
			st = Util.infomask(data[im], mask, False)
			data[im] -= st[0]

	# initialize data for the reference preparation function
	ref_data = [mask3D, ceter, None, None]

	# do the projection matching
	for N_step in xrange(lstp):
 		for Iter in xrange(max_iter):
			print_msg("ITERATION #%3d\n"%(N_step*max_iter + Iter+1))

			if(an[N_step] == -1):	peak, pixel_error = proj_ali_incore(vol, mask3D, data, first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], ref_a, sym, CTF, MPI=False)
			else:	                peak, pixel_error = proj_ali_incore_local(vol, mask3D, data, first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], an[N_step], ref_a, sym, CTF, MPI=False)
			#  3D stuff
			if(CTF): vol1 = recons3d_4nn_ctf(data, range(0,nima,2), 1.0e6, 1, sym)
			else:	 vol1 = recons3d_4nn(data, range(0,nima,2), sym)
			if(CTF): vol2 = recons3d_4nn_ctf(data, range(1,nima,2), 1.0e6, 1, sym)
			else:	 vol2 = recons3d_4nn(data, range(1,nima,2), sym)

			fscc = fsc_mask(vol1, vol2, mask3D, 1.0, os.path.join(outdir, "resolution%04d"%(N_step*max_iter+Iter+1)))
			del vol1
			del vol2

			# calculate new and improved 3D
			if(CTF): vol = recons3d_4nn_ctf(data, range(nima), 1.0e6, 1, sym)
			else:	 vol = recons3d_4nn(data, range(nima), sym)
			# store the reference volume
			drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(N_step*max_iter+Iter+1)))
			did, vav = recons3d_nn_SSNR(data,  mask2D = None, ring_width=1, npad =1, sign=1, symmetry = sym, CTF =CTF)
			vol = fft(Util.divn_filter(fft(vol), vav))
			ref_data[2] = vol
			ref_data[3] = fscc
			#  call user-supplied function to prepare reference image, i.e., center and filter it
			vol, cs = user_func( ref_data )
			if center == 1:
				from utilities import rotate_3D_shift
				rotate_3D_shift(data, [-cs[0], -cs[1], -cs[2]])
			drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(N_step*max_iter+Iter+1)))
	#  here we  write header info
	from utilities import write_headers
	write_headers( stack, data, range(nima))
	print_end_msg("ali3d_a")


def ali3d(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", deltapsi = "-1", startpsi = "-1",
	    center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",
	    user_func_name = "ref_ali3d", fourvar = True, npad = 4, debug = False, MPI = False, termprec = 0.0):
	"""
		Name
			ali3d - Perform 3-D projection matching given initial reference volume and image series
		Input
			stack: set of 2-D images in a stack file, images have to be squares
			ref_vol: initial reference volume
			outdir: directory name into which the results will be written
			maskfile: filename of the file containing 3D mask.
			ir: inner radius for rotational correlation > 0 
			ou: outer radius for rotational correlation <int(nx/2)-1 
			rs: steps between rings in rotational correlation >0
			xr: range for translation search in x direction in each iteration, search is +/xr
			yr: range for translation search in y direction in each iteration, search is +/yr
			ts: step size of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional.
			delta: angular step for the reference projections in respective iterations
			an: angular neighborhood for local searches
			center: average center method
			max_iter: maximum iterations at each angle step
			CTF: if the flag is present, program will use the CTF information stored in file headers
			snr: signal noise ratio used in the 3D reconstruction
			ref_a: method for creating quasi-uniform distribution of the projection directions of reference projections: "S" - spiral
			sym: symmetry of the refined structure
			function: name of the user-supplied-function
			MPI: if presetm use MPI version
		Output
			output_directory: directory name into which the output files will be written.
			header: the alignment parameters are stored in the headers of input files as Transform Object xform.proj
	"""
	if MPI:
		ali3d_MPI(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, yr, ts,
	        	delta, an, deltapsi, startpsi, center, maxit, CTF, snr, ref_a, sym, user_func_name,
			fourvar, npad, debug, termprec)
		return

	from alignment      import proj_ali_incore, proj_ali_incore_local
	from utilities      import model_circle, drop_image, get_image, get_input_from_string
	from utilities      import get_params_proj
	from utilities      import estimate_3D_center, rotate_3D_shift
	from filter         import filt_params, fit_tanh, filt_tanl, filt_ctf
	from statistics     import fsc_mask
	import os
	import types
	from utilities      import print_begin_msg, print_end_msg, print_msg
	print_begin_msg("ali3d")

	from alignment      import Numrinit, prepare_refrings
	from projection     import prep_vol

	import user_functions
	user_func = user_functions.factory[user_func_name]

	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ali3d", 1)
	os.mkdir(outdir)

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(step), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)
	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Reference volume            : %s\n"%(ref_vol))	
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Maskfile                    : %s\n"%(maskfile))
	print_msg("Inner radius                : %i\n"%(first_ring))

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	if last_ring == -1:	last_ring = nx/2 - 2

	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Ring step                   : %i\n"%(rstep))
	print_msg("X search range              : %s\n"%(xrng))
	print_msg("Y search range              : %s\n"%(yrng))
	print_msg("Translational step          : %s\n"%(step))
	print_msg("Angular step                : %s\n"%(delta))
	print_msg("Angular search range        : %s\n"%(an))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("CTF correction              : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	print_msg("Reference projection method : %s\n"%(ref_a))
	print_msg("Symmetry group              : %s\n\n"%(sym))

	if maskfile :
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else                                  : mask3D = maskfile
	else          :   mask3D = model_circle(last_ring, nx, nx, nx)
	mask2D = model_circle(last_ring, nx, nx) - model_circle(first_ring, nx, nx)
	numr   = Numrinit(first_ring, last_ring, rstep, "F")

	if CTF:
		from reconstruction import recons3d_4nn_ctf
		from filter         import filt_ctf
	else: from reconstruction import recons3d_4nn

	if debug:  outf = file(os.path.join(outdir, "progress"), "w")
	else:      outf = None

	active = EMUtil.get_all_attributes(stack, 'active')
	list_of_particles = []
	for im in xrange(len(active)):
		if(active[im]):  list_of_particles.append(im)
	del active
	data = EMData.read_images(stack, list_of_particles)
        for im in xrange(len(data)):
                data[im].set_attr('ID', list_of_particles[im])
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)

	nima = len(data)
	# initialize data for the reference preparation function
	ref_data = [ mask3D, max(center,0), None, None ]#  for center -1 switch of centering by user function

	cs = [0.0]*3
	# do the projection matching
	for N_step in xrange(lstp):
 		for Iter in xrange(max_iter):
			print_msg("\nITERATION #%3d\n"%(N_step*max_iter+Iter+1))

			volft,kb = prep_vol( vol )
			refrings = prepare_refrings( volft, kb, nx, delta[N_step], ref_a, sym, numr, MPI=False)
			del volft,kb

			for im in xrange( nima ):

				if an[N_step] == -1:	
					peak, pixel_error = proj_ali_incore(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step])
				else:
					peak, pixel_error = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step])
			if center == -1:
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center(data)
				msg = "Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
				print_msg(msg)
				rotate_3D_shift(data, [-cs[0], -cs[1], -cs[2]])

			if CTF:   vol1 = recons3d_4nn_ctf(data, range(0, nima, 2), snr, 1, sym)
			else:	   vol1 = recons3d_4nn(data, range(0, nima, 2), sym)
			if CTF:   vol2 = recons3d_4nn_ctf(data, range(1, nima, 2), snr, 1, sym)
			else:	   vol2 = recons3d_4nn(data, range(1, nima, 2), sym)

			fscc = fsc_mask(vol1, vol2, mask3D, 1.0, os.path.join(outdir, "resolution%04d"%(N_step*max_iter+Iter+1)))
			del vol1
			del vol2

			# calculate new and improved 3D
			if CTF:  vol = recons3d_4nn_ctf(data, range(nima), snr, 1, sym)
			else:	 vol = recons3d_4nn(data, range(nima), sym)
			# store the reference volume
			drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(N_step*max_iter+Iter+1)))
			ref_data[2] = vol
			ref_data[3] = fscc

			#  call user-supplied function to prepare reference image, i.e., center and filter it
			vol, dummy = user_func(ref_data)

			drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(N_step*max_iter+Iter+1)))
			#  here we write header info
			from utilities import write_headers
			#from utilities import write_select_headers
			if CTF:
				for dat in data:  dat.set_attr('ctf_applied',0)
			write_headers(stack, data, list_of_particles)
			#list_params= ['ID','xform.projection']
			#write_select_headers(stack, data, list_params)
			if CTF:
				for dat in data:  dat.set_attr('ctf_applied', 1)
	print_end_msg("ali3d")

def ali3d_MPI(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", deltapsi = "-1", startpsi = "-1",
	    center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",  user_func_name = "ref_ali3d",
	    fourvar = True, npad = 4, debug = False, termprec = 0.0):

	from alignment       import Numrinit, prepare_refrings, proj_ali_incore, proj_ali_incore_local
	from utilities       import model_circle, get_image, drop_image, get_input_from_string
	from utilities       import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all
	from utilities       import send_attr_dict
	from utilities       import get_params_proj, file_type
	from fundamentals    import rot_avg_image
	import os
	import types
	from utilities       import print_begin_msg, print_end_msg, print_msg
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi             import mpi_reduce, MPI_INT, MPI_SUM
	from filter          import filt_ctf
	from projection      import prep_vol, prgs
	from statistics      import hist_list, varf3d_MPI


	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ali3d_MPI", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)
	
	if myid == main_node:
		print_begin_msg("ali3d_MPI")
		os.mkdir(outdir)
	mpi_barrier(MPI_COMM_WORLD)

	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(step), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)

	if deltapsi == "-1":
		deltapsi = [-1] * lstp
	else:
		deltapsi = get_input_from_string(deltapsi)

	if startpsi == "-1":
		startpsi = [-1] * lstp
	else:
		startpsi = get_input_from_string(startpsi)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	if last_ring < 0:	last_ring = int(nx/2) - 2

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Reference volume            : %s\n"%(ref_vol))	
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Maskfile                    : %s\n"%(maskfile))
		print_msg("Inner radius                : %i\n"%(first_ring))
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Ring step                   : %i\n"%(rstep))
		print_msg("X search range              : %s\n"%(xrng))
		print_msg("Y search range              : %s\n"%(yrng))
		print_msg("Translational step          : %s\n"%(step))
		print_msg("Angular step                : %s\n"%(delta))
		print_msg("Angular search range        : %s\n"%(an))
		print_msg("Delta psi                   : %s\n"%(deltapsi))
		print_msg("Start psi                   : %s\n"%(startpsi))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Percentage of change for termination: %f\n"%(termprec))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Reference projection method : %s\n"%(ref_a))
		print_msg("Symmetry group              : %s\n\n"%(sym))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else: mask3D = model_circle(last_ring, nx, nx, nx)

	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	mask2D  = model_circle(last_ring,nx,nx) - model_circle(first_ring,nx,nx)

	fscmask = mask3D  #model_circle(last_ring,nx,nx,nx)  For a fancy mask circle would work better  PAP 7/21/11
	if CTF:
		from reconstruction import rec3D_MPI
		from filter         import filt_ctf
	else:	 from reconstruction import rec3D_MPI_noCTF

	if myid == main_node:
       		if(file_type(stack) == "bdb"):
			from EMAN2db import db_open_dict
			dummy = db_open_dict(stack, True)
		active = EMUtil.get_all_attributes(stack, 'active')
		list_of_particles = []
		for im in xrange(len(active)):
			if active[im]:  list_of_particles.append(im)
		del active
		nima = len(list_of_particles)
	else:
		nima = 0
	total_nima = bcast_number_to_all(nima, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()

	data = EMData.read_images(stack, list_of_particles)
	if fourvar:  original_data = []
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		if fourvar: original_data.append(data[im].copy())
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)

	if debug:
		finfo.write( '%d loaded  \n' % nima )
		finfo.flush()
	if myid == main_node:
		# initialize data for the reference preparation function
		ref_data = [ mask3D, max(center,0), None, None, None, None ]
		# for method -1, switch off centering in user function

	from time import time	

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )

	pixer = [0.0]*nima
	cs = [0.0]*3
	total_iter = 0
	# do the projection matching
	for N_step in xrange(lstp):
		terminate = 0
		Iter = -1
 		while(Iter < max_iter-1 and terminate == 0):
			Iter += 1
			total_iter += 1
			if myid == main_node:
				start_time = time()
				print_msg("\nITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f, delta psi = %5.2f, start psi = %5.2f\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step],deltapsi[N_step],startpsi[N_step]))

			volft,kb = prep_vol( vol )
			refrings = prepare_refrings( volft, kb, nx, delta[N_step], ref_a, sym, numr, True)
			del volft,kb
			if myid== main_node:
				print_msg( "Time to prepare rings: %d\n" % (time()-start_time) )
				start_time = time()

			for im in xrange( nima ):

				if deltapsi[N_step] > 0.0:
					from alignment import proj_ali_incore_delta
					peak, pixer[im] = proj_ali_incore_delta(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],startpsi[N_step],deltapsi[N_step],finfo)						
				elif an[N_step] == -1:
					peak, pixer[im] = proj_ali_incore(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],finfo)
				else:
					peak, pixer[im] = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step],finfo)

			if myid == main_node:
				print_msg("Time of alignment = %d\n"%(time()-start_time))
				start_time = time()

			#output pixel errors
			from mpi import mpi_gatherv
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			terminate = 0
			if(myid == main_node):
				recvbuf = map(float, recvbuf)
				from statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				if(region[0] < 0.0):  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				# Terminate if 95% within 1 pixel error
				im = 0
				for lhx in xrange(lhist):
					if(region[lhx] > 1.0): break
					im += histo[lhx]
				precn = 100*float(total_nima-im)/float(total_nima)
				msg = " Number of particles that changed orientations %7d, percentage of total: %5.1f\n"%(total_nima-im, precn)
				print_msg(msg)
				if(precn <= termprec):  terminate = 1
				del region, histo
			del recvbuf
			terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
			terminate = int(terminate[0])

			if(center == -1):
				from utilities      import estimate_3D_center_MPI, rotate_3D_shift
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					print_msg(msg)
				cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(data, cs)

			if CTF: vol, fscc = rec3D_MPI(data, snr, sym, fscmask, os.path.join(outdir, "resolution%04d"%(total_iter)), myid, main_node, npad = npad)
			else:    vol, fscc = rec3D_MPI_noCTF(data, sym, fscmask, os.path.join(outdir, "resolution%04d"%(total_iter)), myid, main_node, npad = npad)

			if myid == main_node:
				print_msg("3D reconstruction time = %d\n"%(time()-start_time))
				start_time = time()

			if fourvar:
			#  Compute Fourier variance
				for im in xrange(nima):
					original_data[im].set_attr( 'xform.projection', data[im].get_attr('xform.projection') )
				varf = varf3d_MPI(original_data, ssnr_text_file = os.path.join(outdir, "ssnr%04d"%(total_iter)), mask2D = None, reference_structure = vol, ou = last_ring, rw = 1.0, npad = 1, CTF = CTF, sign = 1, sym =sym, myid = myid)
				if myid == main_node:
					print_msg("Time to calculate 3D Fourier variance= %d\n"%(time()-start_time))
					start_time = time()
					varf = 1.0/varf
			else:  varf = None

			if myid == main_node:
				drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(total_iter)))
				ref_data[2] = vol
				ref_data[3] = fscc
				ref_data[4] = varf
				#  call user-supplied function to prepare reference image, i.e., center and filter it
				vol, cs = user_func(ref_data)
				drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(total_iter)))

			del varf
			bcast_EMData_to_all(vol, myid, main_node)
			# write out headers, under MPI writing has to be done sequentially
			mpi_barrier(MPI_COMM_WORLD)
			par_str = ['xform.projection', 'ID']
			if myid == main_node:
	   			if(file_type(stack) == "bdb"):
	        			from utilities import recv_attr_dict_bdb
	        			recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        		else:
	        			from utilities import recv_attr_dict
	        			recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
				print_msg("Time to write header information= %d\n"%(time()-start_time))
				start_time = time()
	        	else:	       send_attr_dict(main_node, data, par_str, image_start, image_end)
	if myid == main_node: print_end_msg("ali3d_MPI")


def mref_ali3d(stack, ref_vol, outdir, maskfile=None, focus = None, maxit=1, ir=1, ou=-1, rs=1, 
           xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta="10 6 4 4", an="-1", 
	     center = 1.0, nassign = 3, nrefine = 1, CTF = False, snr = 1.0,  ref_a = "S", sym="c1",
	     user_func_name="ref_ali3d", MPI=False, npad = 4, debug = False, fourvar=False, termprec = 0.0):
	if MPI:
		mref_ali3d_MPI(stack, ref_vol, outdir, maskfile, focus, maxit, ir, ou, rs, xr, yr, ts,
		 delta, an, center, nassign, nrefine, CTF, snr, ref_a, sym, user_func_name, npad, debug, fourvar, termprec)
		return
	from utilities      import model_circle, drop_image, get_image, get_input_from_string
	from utilities      import get_arb_params, set_arb_params, get_im, write_headers
	from projection     import prep_vol, prgs
	from utilities      import get_params_proj, estimate_3D_center
	from alignment      import proj_ali_incore, proj_ali_incore_local, Numrinit, prepare_refrings
	from filter	    import filt_params, filt_tanl
	from fundamentals   import fshift
	from statistics     import fsc_mask
	from utilities      import print_begin_msg, print_end_msg, print_msg

	import os
	import types
	# 2D alignment using rotational ccf in polar coords and linear
	# interpolation	
	print_begin_msg("mref_ali3d")

	import user_functions
	user_func = user_functions.factory[user_func_name]

	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "mref_ali3d", 1)
	os.mkdir(outdir)

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min( len(xrng), len(yrng), len(step), len(delta) )
	if (an == "-1"):
		an = []
		for i in xrange(len(xrng)):   an.append(-1)
		from alignment	  import proj_ali_incore
	else:
		an = get_input_from_string(an)
		from alignment	  import proj_ali_incore_local

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)

	numref = EMUtil.get_image_count(ref_vol)
	for  iref in xrange(numref):
		volref     = EMData()
		volref.read_image(ref_vol, iref)
		volref.write_image(os.path.join(outdir, "volf0000.hdf"), iref)

	nx      = volref.get_xsize()
	if last_ring < 0:	last_ring = nx//2 - 2

	fscmask = model_circle(last_ring, nx, nx, nx)

	import user_functions
	user_func = user_functions.factory[user_func_name]

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Reference volume            : %s\n"%(ref_vol))	
	print_msg("Number of reference volumes : %i\n"%(numref))
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Maskfile                    : %s\n"%(maskfile))
	print_msg("Inner radius                : %i\n"%(first_ring))
	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Ring step                   : %i\n"%(rstep))
	print_msg("X search range              : %s\n"%(xrng))
	print_msg("Y search range              : %s\n"%(yrng))
	print_msg("Translational step          : %s\n"%(step))
	print_msg("Angular step                : %s\n"%(delta))
	print_msg("Angular search range        : %s\n"%(an))
	print_msg("Maximum number of reassignment iterations   : %i\n"%(nassign))
	print_msg("Maximum number of alignment iterations      : %i\n"%(nrefine))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("CTF correction              : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	print_msg("Reference projection method : %s\n"%(ref_a))
	print_msg("Symmetry group              : %s\n\n"%(sym))

	if(maskfile):
		if(type(maskfile) is types.StringType):	 mask3D = get_image(maskfile)
		else: 	                                 mask3D = maskfile
	else        :   mask3D = model_circle(last_ring, nx, nx, nx)
	
	numr = Numrinit(first_ring, last_ring, rstep, "F")
	mask2D = model_circle(last_ring,nx,nx) - model_circle(first_ring,nx,nx)


	if debug:  finfo = file(os.path.join(outdir, "progress"), "w")
	else:      finfo = None

	active = EMUtil.get_all_attributes(stack, 'active')
	list_of_particles = []
	for im in xrange(len(active)):
		if(active[im]):  list_of_particles.append(im)
	del active
	data = EMData.read_images(stack, list_of_particles)
	nima = len(data)

	if CTF :
		#  ERROR if ctf applied
		if data[0].get_attr("ctf_applied") > 0:  ERROR("mref_ali3d does not work for CTF-applied data", "mref_ali3d", 1)
		from reconstruction import recons3d_4nn_ctf
		from filter import filt_ctf
	else   : from reconstruction import recons3d_4nn

	# initialize data for the reference preparation function
	ref_data = [mask3D, center, None, None]

	# do the projection matching
	total_iter = 0
	tr_dummy = Transform({"type":"spider"})

	Niter = int(lstp*maxit*(nassign+nrefine))
	for Iter in xrange(Niter):
		N_step = (Iter%(lstp*(nassign+nrefine)))/(nassign+nrefine)
		if Iter%(nassign+nrefine) < nassign:
			runtype = "ASSIGNMENT"
		else:
			runtype = "REFINEMENT"


		total_iter += 1
		print_msg("%s ITERATION #%3d\n"%(runtype,total_iter))
		peaks = [-1.0e23]*nima
		trans = [tr_dummy]*nima

		cs = [0.0]*3
		for iref in xrange(numref):
			vol = get_im(os.path.join(outdir, "volf%04d.hdf"%( total_iter-1)), iref)
			if(CTF):
				previous_defocus = -1.0
			else:
				volft, kb = prep_vol(vol)
				if runtype=="REFINEMENT":
					refrings = prepare_refrings(volft,kb,nx,delta[N_step],ref_a,sym,numr)
			for im in xrange(nima):
				if(CTF):
					ctf_params = data[im].get_attr("ctf")
					if(ctf_params.defocus != previous_defocus):
						previous_defocus = ctf_params.defocus
						volft,kb = prep_vol(filt_ctf(vol, ctf_params))
					if runtype=="REFINEMENT":
						refrings = prepare_refrings(volft,kb,nx,delta[N_step],ref_a,sym,numr)

				if runtype=="ASSIGNMENT":
					phi,tht,psi,s2x,s2y = get_params_proj(data[im])
					ref = prgs( volft,kb,[phi,tht,psi,-s2x,-s2y] )
					peak = ref.cmp("ccc", data[im], {"mask":mask2D, "negative":0})
					if not(finfo is None):
						finfo.write( "ID,iref,peak: %6d %d %8.5f" % (list_of_particles[im],iref,peak) )
						finfo.flush()
				else:		
					if(an[N_step] == -1):	
						peak, pixel_error = proj_ali_incore(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step])
					else:	           
						peak, pixel_error = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step])
					if not(finfo is None):
						phi,tht,psi,s2x,s2y = get_params_proj(data[im])
						finfo.write( "ID,iref,peak,trans: %6d %d %f %f %f %f %f %f"%(list_of_particles[im],iref,peak,phi,tht,psi,s2x,s2y) )
						finfo.flush()


				if(peak > peaks[im]):
					peaks[im] = peak
					data[im].set_attr('group', iref)
					if runtype=="REFINEMENT":
						trans[im] = data[im].get_attr( "xform.projection" )
					if not(finfo is None):
						finfo.write( " current best\n" )
						finfo.flush()
				else:
					if not(finfo is None):
						finfo.write( "\n" )
						finfo.flush()

		if runtype=="REFINEMENT":
			for im in xrange(nima):
				data[im].set_attr('xform.projection', trans[im])

			if center == -1:
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center(data, total_nima, myid, number_of_proc, main_node)
				msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
				print_msg(msg)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(data, cs)


		fscc = [None]*numref
		for iref in xrange(numref):
			list_p = []
			for im in xrange(nima):
				if(iref == data[im].get_attr('group')):
					list_p.append(im)
			print_msg("Group number : %i"%(iref) + ",  number of objects: %i\n"%(len(list_p)))
			#  3D stuff
			if(CTF): vodd = recons3d_4nn_ctf(data, [list_p[im] for im in xrange(0,len(list_p), 2)], snr, 1, sym)
			else:    vodd = recons3d_4nn(data, [list_p[im] for im in xrange(1,len(list_p), 2)], sym)
			if(CTF): veve = recons3d_4nn_ctf(data, [list_p[im] for im in xrange(0,len(list_p), 2)], snr, 1, sym)
			else:    veve = recons3d_4nn(data,[list_p[im] for im in xrange(1,len(list_p), 2)], sym)

			fscc[iref] = fsc_mask(vodd, veve, mask3D, 1.0, os.path.join(outdir, "resolution_%02d_%04d"%(iref, total_iter)))
			
			# calculate new and improved 3D
			if(CTF): volref = recons3d_4nn_ctf(data, list_p, snr, 1, sym)
			else:	 volref = recons3d_4nn(data, list_p, sym)
			volref.write_image(os.path.join(outdir, "vol%04d.hdf"%( total_iter)), iref)
			

		refdata = [None]*7
		refdata[0] = numref
		refdata[1] = outdir
		refdata[2] = fscc
		refdata[3] = total_iter
		refdata[4] = None #varf
		refdata[5] = fscmask
		refdata[6] = False
		user_func( refdata )

		#  here we  write header info
		write_headers( stack, data, list_of_particles)

	print_end_msg("mref_ali3d")


"""
# This is version with the same number of images per group.
def mref_ali3d_MPI(stack, ref_vol, outdir, maskfile=None, focus = None, maxit=1, ir=1, ou=-1, rs=1, 
            xr ="4 2  2  1", yr="-1", ts="1 1 0.5 0.25",   delta="10  6  4  4", an="-1",
	      center = -1, nassign = 3, nrefine= 1, CTF = False, snr = 1.0,  ref_a="S", sym="c1",
	      user_func_name="ref_ali3d", npad = 4, debug = False, fourvar=False, termprec = 0.0):
	from utilities      import model_circle, reduce_EMData_to_root, bcast_EMData_to_all, bcast_number_to_all, drop_image
	from utilities      import bcast_string_to_all, bcast_list_to_all, get_image, get_input_from_string, get_im
	from utilities      import get_arb_params, set_arb_params, drop_spider_doc, send_attr_dict
	from utilities      import get_params_proj, set_params_proj, model_blank
	from filter         import filt_params, filt_btwl, filt_ctf, filt_table, fit_tanh, filt_tanl
	from utilities      import rotate_3D_shift,estimate_3D_center_MPI
	from alignment      import Numrinit, prepare_refrings, proj_ali_incore
	from random         import randint, random
	from filter         import filt_ctf
	from utilities      import print_begin_msg, print_end_msg, print_msg
	from projection     import prep_vol, prgs, project, prgq, gen_rings_ctf
	from morphology     import binarize

	import os
	import types
	from mpi            import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi            import mpi_reduce, mpi_gatherv, mpi_scatterv, MPI_INT, MPI_SUM

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	if myid == 0:
		if os.path.exists(outdir):  nx = 1
		else:  nx = 0
	else:  nx = 0
	ny = bcast_number_to_all(nx, source_node = main_node)
	
	if ny == 1:  ERROR('Output directory exists, please change the name and restart the program', "mref_ali3d_MPI", 1,myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:	
		os.mkdir(outdir)
	mpi_barrier(MPI_COMM_WORLD)

	from time import time	

	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		finfo = open(os.path.join(outdir, "progress%04d"%myid), 'w')
		frec  = open( os.path.join(outdir, "recons%04d"%myid), "w" )
	else:
		finfo = None
		frec  = None

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min( len(xrng), len(yrng), len(step), len(delta) )
	if (an == "-1"):
		an = []
		for i in xrange(len(xrng)):   an.append(-1)
	else:
		from  alignment	    import proj_ali_incore_local
		an      = get_input_from_string(an)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	center      = int(center)

	numref = EMUtil.get_image_count(ref_vol)
	volref     = EMData()
	volref.read_image(stack, 0)
	nx      = volref.get_xsize()
	if last_ring < 0:	last_ring = nx//2 - 2

	if (myid == main_node):
		import user_functions
		user_func = user_functions.factory[user_func_name]
		print_begin_msg("mref_ali3d_MPI")
		print_msg("Input stack                               : %s\n"%(stack))
		print_msg("Reference volumes                         : %s\n"%(ref_vol))	
		print_msg("Number of reference volumes               : %i\n"%(numref))
		print_msg("Output directory                          : %s\n"%(outdir))
		print_msg("User function                             : %s\n"%(user_func_name))
		if(focus != None):  \
		print_msg("Maskfile 3D for focused clustering        : %s\n"%(focus))
		print_msg("Maskfile 2D                               : %s\n"%(maskfile))
		print_msg("Inner radius                              : %i\n"%(first_ring))
		print_msg("Outer radius                              : %i\n"%(last_ring))
		print_msg("Ring step                                 : %i\n"%(rstep))
		print_msg("X search range                            : %s\n"%(xrng))
		print_msg("Y search range                            : %s\n"%(yrng))
		print_msg("Translational step                        : %s\n"%(step))
		print_msg("Angular step                              : %s\n"%(delta))
		print_msg("Angular search range                      : %s\n"%(an))
		print_msg("Number of assignments in each iteration   : %i\n"%(nassign))
		print_msg("Number of alignments in each iteration    : %i\n"%(nrefine))
		print_msg("Number of iterations                      : %i\n"%(lstp*maxit) )
		print_msg("Center type                               : %i\n"%(center))
		print_msg("CTF correction                            : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio                     : %f\n"%(snr))
		print_msg("Reference projection method               : %s\n"%(ref_a))
		print_msg("Symmetry group                            : %s\n\n"%(sym))
		print_msg("Percentage of change for termination      : %f\n"%(termprec))

	if(maskfile):
		if(type(maskfile) is types.StringType): mask3D = get_image(maskfile)
		else: 	                                mask3D = maskfile
	else        :  mask3D = model_circle(last_ring, nx, nx, nx)

	numr     = Numrinit(first_ring, last_ring, rstep, "F")
	mask2D   = model_circle(last_ring, nx, nx) - model_circle(first_ring, nx, nx)

	if(myid == main_node):
		active = EMUtil.get_all_attributes(stack, 'active')
		list_of_particles = []
		for im in xrange(len(active)):
			if(active[im]):  list_of_particles.append(im)
		del active
		total_nima = len(list_of_particles)
	else:
		total_nima =0

	total_nima = bcast_number_to_all(total_nima, source_node = main_node)

	if(myid != main_node):
		list_of_particles = [-1]*total_nima

	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)

	if debug:
		finfo.write( "Image_start, image_end: %d %d\n" %(image_start, image_end) )
		finfo.flush()

	start_time = time()
	data = EMData.read_images(stack, list_of_particles)
	if(myid == 0):
		print_msg( "Time to read data: %d\n" % (time()-start_time) );start_time = time()
	#  Initialize Particle ID and set group number to non-existant -1
 	for im in xrange(nima):
 		data[im].set_attr_dict({'ID':list_of_particles[im], 'group':-1})

	if fourvar:
		from reconstruction import rec3D_MPI
		from statistics     import varf3d_MPI
		#  Compute Fourier variance
		vol, fscc = rec3D_MPI(data, snr, sym, model_circle(last_ring, nx, nx, nx), os.path.join(outdir, "resolution0000"), myid, main_node, finfo=frec, npad=npad)
		varf = varf3d_MPI(data, os.path.join(outdir, "ssnr0000"), None, vol, last_ring, 1.0, 1, CTF, 1, sym, myid)
		if myid == main_node:   
			varf = 1.0/varf
			varf.write_image( os.path.join(outdir,"varf0000.hdf") )
	else:
		varf = None

	if myid == main_node:
		for  iref in xrange(numref):
			get_im(ref_vol, iref).write_image(os.path.join(outdir, "volf0000.hdf"), iref)
	mpi_barrier( MPI_COMM_WORLD )

	if CTF:
		if(data[0].get_attr("ctf_applied") > 0.0):  ERROR("mref_ali3d_MPI does not work for CTF-applied data", "mref_ali3d_MPI", 1, myid)
		from reconstruction import rec3D_MPI
	else:
		from reconstruction import rec3D_MPI_noCTF

	if debug:
		finfo.write( '%d loaded  \n' % len(data) )
		finfo.flush()

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )

	total_iter = 0
	tr_dummy = Transform({"type":"spider"})

	if(focus != None):
		if(myid == main_node):
			vol = get_im(focus)
		else:
			vol =  model_blank(nx, nx, nx)
		bcast_EMData_to_all(vol, myid, main_node)
		focus, kb = prep_vol(vol)

	Niter = int(lstp*maxit*(nassign + nrefine) )
	for Iter in xrange(Niter):
		N_step = (Iter%(lstp*(nassign+nrefine)))/(nassign+nrefine)
		if Iter%(nassign+nrefine) < nassign:
			runtype = "ASSIGNMENT"
		else:
			runtype = "REFINEMENT"

		total_iter += 1
		if(myid == main_node):
			print_msg("\n%s ITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f\n"%(runtype, total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step]))
			start_ime = time()
	
		peaks =  [ [ -1.0e23 for im in xrange(nima) ] for iref in xrange(numref) ]
		if runtype=="REFINEMENT":
 			trans = [ [ tr_dummy for im in xrange(nima) ] for iref in xrange(numref) ]
			pixer = [ [  0.0     for im in xrange(nima) ] for iref in xrange(numref) ]

		cs = [0.0]*3
		for iref in xrange(numref):
			if(myid == main_node):
				vol = get_im(os.path.join(outdir, "volf%04d.hdf"%(total_iter-1)), iref)
			else:
				vol =  model_blank(nx, nx, nx)
			bcast_EMData_to_all(vol, myid, main_node)

			volft, kb = prep_vol(vol)
			del vol
			if CTF:
				previous_defocus = -1.0
				if runtype=="REFINEMENT":
					start_time = time()
					prjref = prgq( volft, kb, nx, delta[N_step], ref_a, sym, MPI=True)
					if(myid == 0):
						print_msg( "Calculation of projections: %d\n" % (time()-start_time) );start_time = time()
					del volft, kb

			else:
				if runtype=="REFINEMENT":
					start_time = time()
					refrings = prepare_refrings( volft, kb, nx, delta[N_step], ref_a, sym, numr)
					if(myid == 0):
						print_msg( "Initial time to prepare rings: %d\n" % (time()-start_time) );start_time = time()
					del volft, kb


			start_time = time()
			for im in xrange(nima):
				if(CTF):
					ctf = data[im].get_attr( "ctf" )
					if runtype=="REFINEMENT":
						if(ctf.defocus != previous_defocus):
							previous_defocus = ctf.defocus
							rstart_time = time()
							refrings = gen_rings_ctf( prjref, nx, ctf, numr)
							if(myid == 0):
								print_msg( "Repeated time to prepare rings: %d\n" % (time()-rstart_time) );rstart_time = time()

				if runtype=="ASSIGNMENT":
					phi,tht,psi,s2x,s2y = get_params_proj(data[im])
					ref = prgs( volft, kb, [phi,tht,psi,-s2x,-s2y])
					if CTF:  ref = filt_ctf( ref, ctf )
					if(focus != None):  mask2D = binarize( prgs( focus, kb, [phi,tht,psi,-s2x,-s2y]) )
					peak = ref.cmp("ccc",data[im],{"mask":mask2D, "negative":0})
					if not(finfo is None):
						finfo.write( "ID,iref,peak: %6d %d %8.5f" % (list_of_particles[im],iref,peak) )
				else:
					if(an[N_step] == -1):
						peak, pixel_error = proj_ali_incore(data[im], refrings, numr, xrng[N_step], yrng[N_step], step[N_step])
					else:
						peak, pixel_error = proj_ali_incore_local(data[im], refrings, numr, xrng[N_step], yrng[N_step], step[N_step], an[N_step])
					if not(finfo is None):
						phi,tht,psi,s2x,s2y = get_params_proj(data[im])
						finfo.write( "ID,iref,peak,trans: %6d %d %f %f %f %f %f %f"%(list_of_particles[im],iref,peak,phi,tht,psi,s2x,s2y) )
						finfo.flush()

				peaks[iref][im] = peak
				if runtype=="REFINEMENT":
					pixer[iref][im] = pixel_error
					trans[iref][im] = data[im].get_attr( "xform.projection" )
			if(myid == 0):
				print_msg( "Time to process particles for reference %3d: %d\n" % (iref, time()-start_time) );start_time = time()


		if runtype=="ASSIGNMENT":  del volft, kb, ref
		else:
			if CTF: del prjref
			del refrings

		#  send peak values to the main node, do the assignments, and bring them back
		from numpy import float32, empty, inner, abs
		if( myid == 0 ):
			dtot = empty( (numref, total_nima), dtype = float32)
		for  iref in xrange(numref):
			recvbuf = mpi_gatherv(peaks[iref], nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			if( myid == 0 ): dtot[iref] = recvbuf
		#  prepare reference directions
		from utilities import even_angles, getvec
		refa = even_angles(15.0)
		numrefang = len(refa)
		refanorm = empty( (numrefang, 3), dtype = float32)
		for i in xrange(numrefang):
			tmp = getvec(refa[i][0], refa[i][1])
			for j in xrange(3):
				refanorm[i][j] = tmp[j]
		del  refa, tmp

		transv = empty( (nima, 3), dtype = float32)
		if runtype=="ASSIGNMENT":
			for im in xrange(nima):
				trns = data[im].get_attr( "xform.projection" )
				for j in xrange(3):
					transv[im][j] = trns.at(2,j)
		else:
			# For REFINEMENT we have a problem, as the exact angle is not known only after the next step of assigning projections.
			# So, we will assume it is the one with max peak
			for im in xrange(nima):
				qt = -1.0e23
				it = -1
				for iref in xrange(numref):
					pt = peaks[iref][im]
					if(pt > qt):
						qt = pt
						it = iref
				for j in xrange(3):
					transv[im][j] = trans[it][im].at(2,j)
		#  We have all vectors, now create a list of assignments of images to references
		refassign = [-1]*nima
		for im in xrange(nima):
			refassign[im] = abs(inner(refanorm,transv[im])).argmax()
		assigntorefa = mpi_gatherv(refassign, nima, MPI_INT, recvcount, disps, MPI_INT, main_node, MPI_COMM_WORLD)
		assigntorefa = map(int, assigntorefa)

		del refassign, refanorm, transv
		del recvbuf


		if myid == main_node:
			SA = False
			asi = [[] for iref in xrange(numref)]
			report_error = 0
			for imrefa in xrange(numrefang):
				from utilities import findall
				N = findall(assigntorefa, imrefa)
				current_nima = len(N)
				if( current_nima >= numref and report_error == 0):
					tasi = [[] for iref in xrange(numref)]
					maxasi = current_nima//numref
					nt = current_nima
					kt = numref
					K = range(numref)

					d = empty( (numref, current_nima), dtype = float32)
					for ima in xrange(current_nima):
						for iref in xrange(numref):  d[iref][ima] = dtot[iref][N[ima]]

					while nt > 0 and kt > 0:
						l = d.argmax()
						group = l//current_nima
						ima   = l-current_nima*group
						if SA:
							J = [0.0]*numref
							sJ = 0
							Jc = [0.0]*numref
							for iref in xrange(numref):
								J[iref] = exp(d[iref][ima]/T)
								sJ += J[iref]
							for iref in xrange(numref):
								J[iref] /= sJ
							Jc[0] = J[0]
							for iref in xrange(1, numref):
								Jc[iref] = Jc[iref-1]+J[iref]
							sss = random()
							for group in xrange(numref):
								if( sss <= Jc[group]): break
						tasi[group].append(N[ima])
						N[ima] = -1
						for iref in xrange(numref):  d[iref][ima] = -1.e10
						nt -= 1
						masi = len(tasi[group])
						if masi == maxasi:
							for im in xrange(current_nima):  d[group][im] = -1.e10
							kt -= 1
					else:
						for ima in xrange(current_nima):
							if N[ima] > -1:
								qm = -1.e10
								for iref in xrange(numref):
									qt = dtot[iref][N[ima]]
									if( qt > qm ):
										qm = qt
										group = iref
								tasi[group].append(N[ima])

					del d, N, K
					if  SA:  del J, Jc
					for iref in xrange(numref):
						asi[iref] += tasi[iref]
					del tasi
				else:
					report_error = 1
			del dtot

		else:
			assignment = []
			report_error = 0

		report_error = bcast_number_to_all(report_error, source_node = main_node)
		if report_error == 1:  ERROR('Number of images within a group too small', "mref_ali3d_MPI", 1, myid)
		if myid == main_node:
			assignment = [0]*total_nima
			for iref in xrange(numref):
				for im in xrange(len(asi[iref])):
					assignment[asi[iref][im]] = iref
			del asi

		'''
		if myid == main_node:
			SA = False
			maxasi = total_nima//numref
			asi = [[] for iref in xrange(numref)]
			nt = total_nima
			kt = numref
			K = range(numref)
			N = range(total_nima)

			while nt > 0 and kt > 1:
				l = d.argmax()
				group = l//total_nima
				ima   = l-total_nima*group
				if SA:
					J = [0.0]*numref
					sJ = 0
					Jc = [0.0]*numref
					for iref in xrange(numref):
						J[iref] = exp(d[iref][ima]/T)
						sJ += J[iref]
					for iref in xrange(numref):
						J[iref] /= sJ
					Jc[0] = J[0]
					for iref in xrange(1, numref):
						Jc[iref] = Jc[iref-1]+J[iref]
					sss = random()
					for group in xrange(numref):
						if( sss <= Jc[group]): break
				asi[group].append(N[ima])
				for iref in xrange(numref):  d[iref][ima] = -1.e10
				nt -= 1
				masi = len(asi[group])
				if masi == maxasi:
					for im in xrange(total_nima):  d[group][im] = -1.e10
					kt -= 1
			else:
				mas = [len(asi[iref]) for iref in xrange(numref)]
				group = mas.index(min(mas))
				del mas
				for im in xrange(total_nima):
					kt = 0
					go = True
					while(go and kt < numref):
						if d[kt][im] > -1.e10:
							asi[group].append(im)
							go = False
						kt += 1

			assignment = [0]*total_nima
			for iref in xrange(numref):
				for im in xrange(len(asi[iref])):
					assignment[asi[iref][im]] = iref

			del asi, d, N, K
			if  SA:  del J, Jc


		else:
			assignment = []
		'''

		assignment = mpi_scatterv(assignment, recvcount, disps, MPI_INT, recvcount[myid], MPI_INT, main_node, MPI_COMM_WORLD)
		assignment = map(int, assignment)


		#  compute number of particles that changed assignment and how many are in which group
		nchng = 0
		npergroup = [0]*numref
		for im in xrange(nima):
			iref = data[im].get_attr('group')
			npergroup[assignment[im]] += 1
			if( iref != assignment[im]): nchng += 1
			data[im].set_attr('group', assignment[im])
		nchng = mpi_reduce(nchng, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
		npergroup = mpi_reduce(npergroup, numref, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
		npergroup = map(int, npergroup)
		terminate = 0
		if( myid == 0 ):
			nchng = int(nchng[0])
			precn = 100*float(nchng)/float(total_nima)
			msg = " Number of particles that changed assignments %7d, percentage of total: %5.1f\n"%(nchng, precn)
			print_msg(msg)
			msg = " Group       number of particles\n"
			print_msg(msg)
			for iref in xrange(numref):
				msg = " %5d       %7d\n"%(iref+1, npergroup[iref])
				print_msg(msg)
			if(precn <= termprec):  terminate = 1
		terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
		terminate = int(terminate[0])

		if runtype=="REFINEMENT":
			for im in xrange(nima):
				data[im].set_attr('xform.projection', trans[assignment[im]][im])
				pixer[0][im] = pixer[assignment[im]][im]
			pixer = pixer[0]

			if(center == -1):
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)				
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					print_msg(msg)
				cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(data, cs)
			#output pixel errors
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			if(myid == main_node):
				recvbuf = map(float, recvbuf)
				from statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				if(region[0] < 0.0):  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.3f      %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				del region, histo
			del recvbuf

		fscc = [None]*numref

		if fourvar and runtype=="REFINEMENT":
			sumvol = model_blank(nx, nx, nx)

		start_time = time()
		for iref in xrange(numref):
			#  3D stuff
			from time import localtime, strftime
			if(CTF): volref, fscc[iref] = rec3D_MPI(data, snr, sym, model_circle(last_ring, nx, nx, nx), os.path.join(outdir, "resolution_%02d_%04d"%(iref, total_iter)), myid, main_node, index = iref, npad = npad, finfo=frec)
			else:    volref, fscc[iref] = rec3D_MPI_noCTF(data, sym, model_circle(last_ring, nx, nx, nx), os.path.join(outdir, "resolution_%02d_%04d"%(iref, total_iter)), myid, main_node, index = iref, npad = npad, finfo=frec)
			if(myid == 0):
				print_msg( "Time to compute 3D: %d\n" % (time()-start_time) );start_time = time()

			if(myid == main_node):
				volref.write_image(os.path.join(outdir, "vol%04d.hdf"%( total_iter)), iref)
				if fourvar and runtype=="REFINEMENT":
					sumvol += volref
			del volref

		if runtype=="REFINEMENT":
			if fourvar:
				varf = varf3d_MPI(data, os.path.join(outdir, "ssnr%04d"%total_iter), None,sumvol,last_ring, 1.0, 1, CTF, 1, sym, myid)
				if myid == main_node:   
					varf = 1.0/varf
					varf.write_image( os.path.join(outdir,"varf%04d.hdf"%total_iter) )


		if(myid == main_node):
			refdata = [None]*7
			refdata[0] = numref
			refdata[1] = outdir
			refdata[2] = None
			refdata[3] = total_iter
			refdata[4] = varf
			refdata[5] = mask3D
			refdata[6] = (runtype=="REFINEMENT") # whether to align on 50S, this only happens at refinement step
			user_func( refdata )

		#  here we  write header info
		mpi_barrier(MPI_COMM_WORLD)
		start_time = time()
		if runtype=="REFINEMENT":
			par_str = ['xform.projection', 'ID', 'group']
		else:
			par_str = ['group', 'ID' ]
	        if myid == main_node:
			from utilities import file_type
	        	if(file_type(stack) == "bdb"):
	        		from utilities import recv_attr_dict_bdb
	        		recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        	else:
	        		from utilities import recv_attr_dict
	        		recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        else:		send_attr_dict(main_node, data, par_str, image_start, image_end)
		if(terminate == 1):
			if myid==main_node:
				print_end_msg("mref_ali3d_MPI terminated due to small number of objects changing assignments")
			from sys import exit
			exit()
		if(myid == 0):
			print_msg( "Time to write headers: %d\n" % (time()-start_time) );start_time = time()
	if myid==main_node:
		print_end_msg("mref_ali3d_MPI")
"""

def mref_ali3d_MPI(stack, ref_vol, outdir, maskfile=None, focus = None, maxit=1, ir=1, ou=-1, rs=1, 
            xr ="4 2  2  1", yr="-1", ts="1 1 0.5 0.25",   delta="10  6  4  4", an="-1",
	      center = -1, nassign = 3, nrefine= 1, CTF = False, snr = 1.0,  ref_a="S", sym="c1",
	      user_func_name="ref_ali3d", npad = 4, debug = False, fourvar=False, termprec = 0.0):
	from utilities      import model_circle, reduce_EMData_to_root, bcast_EMData_to_all, bcast_number_to_all, drop_image
	from utilities      import bcast_string_to_all, bcast_list_to_all, get_image, get_input_from_string, get_im
	from utilities      import get_arb_params, set_arb_params, drop_spider_doc, send_attr_dict
	from utilities      import get_params_proj, set_params_proj, model_blank
	from filter         import filt_params, filt_btwl, filt_ctf, filt_table, fit_tanh, filt_tanl
	from utilities      import rotate_3D_shift,estimate_3D_center_MPI
	from alignment      import Numrinit, prepare_refrings, proj_ali_incore
	from random         import randint
	from filter         import filt_ctf
	from utilities      import print_begin_msg, print_end_msg, print_msg
	from projection     import prep_vol, prgs, project, prgq, gen_rings_ctf
	import os
	import types
	from mpi            import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi            import mpi_reduce, MPI_INT, MPI_SUM

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	if os.path.exists(outdir): ERROR('Output directory exists, please change the name and restart the program', "mref_ali3d_MPI ", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:	
		print_begin_msg("mref_ali3d_MPI")
		os.mkdir(outdir)
	mpi_barrier(MPI_COMM_WORLD)

	from time import time	

	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		finfo = open(os.path.join(outdir, "progress%04d"%myid), 'w')
		frec  = open( os.path.join(outdir, "recons%04d"%myid), "w" )
	else:
		finfo = None
		frec  = None

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min( len(xrng), len(yrng), len(step), len(delta) )
	if (an == "-1"):
		an = []
		for i in xrange(len(xrng)):   an.append(-1)
	else:
		from  alignment	    import proj_ali_incore_local
		an      = get_input_from_string(an)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	center      = int(center)

	numref = EMUtil.get_image_count(ref_vol)
	volref     = EMData()
	volref.read_image(stack, 0)
	nx      = volref.get_xsize()
	if last_ring < 0:	last_ring = nx//2 - 2

	fscmask = model_circle(last_ring, nx, nx, nx)

	if (myid == main_node):
		import user_functions
		user_func = user_functions.factory[user_func_name]
		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Reference volumes           : %s\n"%(ref_vol))	
		print_msg("Number of reference volumes : %i\n"%(numref))
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("User function               : %s\n"%(user_func_name))
		print_msg("Maskfile                    : %s\n"%(maskfile))
		print_msg("Inner radius                : %i\n"%(first_ring))
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Ring step                   : %i\n"%(rstep))
		print_msg("X search range              : %s\n"%(xrng))
		print_msg("Y search range              : %s\n"%(yrng))
		print_msg("Translational step          : %s\n"%(step))
		print_msg("Angular step                : %s\n"%(delta))
		print_msg("Angular search range        : %s\n"%(an))
		print_msg("Number of assignments in each iteration   : %i\n"%(nassign))
		print_msg("Number of alignments in each iteration    : %i\n"%(nrefine))
		print_msg("Number of iterations                      : %i\n"%(lstp*maxit) )
		print_msg("Center type                 : %i\n"%(center))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Reference projection method : %s\n"%(ref_a))
		print_msg("Symmetry group              : %s\n\n"%(sym))
		print_msg("Percentage of change for termination: %f\n"%(termprec))

	if(maskfile):
		if(type(maskfile) is types.StringType): mask3D = get_image(maskfile)
		else: 	                                mask3D = maskfile
	else        :  mask3D = model_circle(last_ring, nx, nx, nx)

	numr     = Numrinit(first_ring, last_ring, rstep, "F")
	mask2D   = model_circle(last_ring, nx, nx) - model_circle(first_ring, nx, nx)

	if(myid == main_node):
		active = EMUtil.get_all_attributes(stack, 'active')
		list_of_particles = []
		for im in xrange(len(active)):
			if(active[im]):  list_of_particles.append(im)
		del active
		nima = len(list_of_particles)
	else:
		nima =0

	nima = bcast_number_to_all(nima, source_node = main_node)

	if(myid != main_node):
		list_of_particles = [-1]*nima

	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
	# create a list of images for each node
	total_nima = nima
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)

	if debug:
		finfo.write( "Image_start, image_end: %d %d\n" %(image_start, image_end) )
		finfo.flush()

	start_time = time()
	data = EMData.read_images(stack, list_of_particles)
	if(myid == 0):
		print_msg( "Time to read data: %d\n" % (time()-start_time) );start_time = time()
	#  Initialize Particle ID and set group number to non-existant -1
	assignment = [-1]*len(data)
 	for im in xrange(len(data)):
 		data[im].set_attr_dict({'ID':list_of_particles[im], 'group':-1})

	if fourvar:
		from reconstruction import rec3D_MPI
		from statistics     import varf3d_MPI
		#  Compute Fourier variance
		vol, fscc = rec3D_MPI(data, snr, sym, fscmask, os.path.join(outdir, "resolution0000"), myid, main_node, finfo=frec, npad=npad)
		varf = varf3d_MPI(data, os.path.join(outdir, "ssnr0000"), None, vol, last_ring, 1.0, 1, CTF, 1, sym, myid)
		if myid == main_node:   
			varf = 1.0/varf
			varf.write_image( os.path.join(outdir,"varf0000.hdf") )
	else:
		varf = None


	if myid == main_node:
		for  iref in xrange(numref):
			get_im(ref_vol, iref).write_image(os.path.join(outdir, "volf0000.hdf"), iref)
	mpi_barrier( MPI_COMM_WORLD )

	if CTF:
		if(data[0].get_attr("ctf_applied") > 0.0):  ERROR("mref_ali3d_MPI does not work for CTF-applied data", "mref_ali3d_MPI", 1, myid)
		from reconstruction import rec3D_MPI
	else:
		from reconstruction import rec3D_MPI_noCTF

	if debug:
		finfo.write( '%d loaded  \n' % len(data) )
		finfo.flush()

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )

	total_iter = 0
	tr_dummy = Transform({"type":"spider"})

	Niter = int(lstp*maxit*(nassign + nrefine) )
	for Iter in xrange(Niter):
		N_step = (Iter%(lstp*(nassign+nrefine)))/(nassign+nrefine)
		if Iter%(nassign+nrefine) < nassign:
			runtype = "ASSIGNMENT"
		else:
			runtype = "REFINEMENT"

		total_iter += 1
		if(myid == main_node):
			print_msg("\n%s ITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f\n"%(runtype, total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step]))
			start_ime = time()
	
		peaks = [ -1.0e23]*nima
 		trans = [tr_dummy]*nima
		if runtype=="REFINEMENT":  pixer = [0.0]*nima

		cs = [0.0]*3
		for iref in xrange(numref):
			vol = get_im(os.path.join(outdir, "volf%04d.hdf"%(total_iter-1)), iref)
			volft, kb = prep_vol(vol)
			if CTF:
				previous_defocus = -1.0
				if runtype=="REFINEMENT":
					start_time = time()
					prjref = prgq( volft, kb, nx, delta[N_step], ref_a, sym, MPI=True)
					if(myid == 0):
						print_msg( "Calculation of projections: %d\n" % (time()-start_time) );start_time = time()
					del volft, kb

			else:
				if runtype=="REFINEMENT":
					start_time = time()
					refrings = prepare_refrings( volft, kb, nx, delta[N_step], ref_a, sym, numr)
					if(myid == 0):
						print_msg( "Initial time to prepare rings: %d\n" % (time()-start_time) );start_time = time()
					del volft, kb


			start_time = time()
			for im in xrange(nima):
				if(CTF):
					ctf = data[im].get_attr( "ctf" )
					if runtype=="REFINEMENT":
						if(ctf.defocus != previous_defocus):
							previous_defocus = ctf.defocus
							rstart_time = time()
							refrings = gen_rings_ctf( prjref, nx, ctf, numr)
							if(myid == 0):
								print_msg( "Repeated time to prepare rings: %d\n" % (time()-rstart_time) );rstart_time = time()

				if runtype=="ASSIGNMENT":
					phi,tht,psi,s2x,s2y = get_params_proj(data[im])
					ref = prgs( volft, kb, [phi,tht,psi,-s2x,-s2y])
					if CTF:  ref = filt_ctf( ref, ctf )
					peak = ref.cmp("ccc",data[im],{"mask":mask2D, "negative":0})
					if not(finfo is None):
						finfo.write( "ID,iref,peak: %6d %d %8.5f" % (list_of_particles[im],iref,peak) )
				else:
					if(an[N_step] == -1):
						peak, pixel_error = proj_ali_incore(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step])
					else:
						peak, pixel_error = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step])
					if not(finfo is None):
						phi,tht,psi,s2x,s2y = get_params_proj(data[im])
						finfo.write( "ID,iref,peak,trans: %6d %d %f %f %f %f %f %f"%(list_of_particles[im],iref,peak,phi,tht,psi,s2x,s2y) )
						finfo.flush()

				if(peak > peaks[im]):
					peaks[im] = peak
					data[im].set_attr('group', iref)
					if runtype=="REFINEMENT":
						pixer[im] = pixel_error
						trans[im] = data[im].get_attr( "xform.projection" )
					if not(finfo is None):
						finfo.write( " current best\n" )
						finfo.flush()
				else:
					if not(finfo is None):
						finfo.write( "\n" )
						finfo.flush()
			if(myid == 0):
				print_msg( "Time to process particles for reference %3d: %d\n" % (iref, time()-start_time) );start_time = time()


		del peaks
		if runtype=="ASSIGNMENT":  del volft, kb, ref
		else:
			if CTF: del prjref
			del refrings
		#  compute number of particles that changed assignment and how man are in which group
		nchng = 0
		npergroup = [0]*numref
		for im in xrange(nima):
			iref = data[im].get_attr('group')
			npergroup[iref] += 1
			if( iref != assignment[im]):
				assignment[im] = iref
				nchng += 1
		nchng = mpi_reduce(nchng, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
		npergroup = mpi_reduce(npergroup, numref, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
		npergroup = map(int, npergroup)
		terminate = 0
		if( myid == 0 ):
			nchng = int(nchng[0])
			precn = 100*float(nchng)/float(total_nima)
			msg = " Number of particles that changed assignments %7d, percentage of total: %5.1f\n"%(nchng, precn)
			print_msg(msg)
			msg = " Group       number of particles\n"
			print_msg(msg)
			for iref in xrange(numref):
				msg = " %5d       %7d\n"%(iref+1, npergroup[iref])
				print_msg(msg)
			if(precn <= termprec):  terminate = 1
		terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
		terminate = int(terminate[0])

		if runtype=="REFINEMENT":
			for im in xrange(nima):
				data[im].set_attr('xform.projection', trans[im])

			if(center == -1):
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)				
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					print_msg(msg)
				cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(data, cs)
			#output pixel errors
			from mpi import mpi_gatherv
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			if(myid == main_node):
				recvbuf = map(float, recvbuf)
				from statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				if(region[0] < 0.0):  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.3f      %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				del region, histo
			del recvbuf

		if CTF: del vol
		fscc = [None]*numref

		if fourvar and runtype=="REFINEMENT":
			sumvol = model_blank(nx, nx, nx)

		sart_time = time()
		for iref in xrange(numref):
			#  3D stuff
			from time import localtime, strftime
			if(CTF): volref, fscc[iref] = rec3D_MPI(data, snr, sym, fscmask, os.path.join(outdir, "resolution_%02d_%04d"%(iref, total_iter)), myid, main_node, index = iref, npad = npad, finfo=frec)
			else:    volref, fscc[iref] = rec3D_MPI_noCTF(data, sym, fscmask, os.path.join(outdir, "resolution_%02d_%04d"%(iref, total_iter)), myid, main_node, index = iref, npad = npad, finfo=frec)
			if(myid == 0):
				print_msg( "Time to compute 3D: %d\n" % (time()-start_time) );start_time = time()

			if(myid == main_node):
				volref.write_image(os.path.join(outdir, "vol%04d.hdf"%( total_iter)), iref)
				if fourvar and runtype=="REFINEMENT":
					sumvol += volref
			del volref

		if runtype=="REFINEMENT":
			if fourvar:
				varf = varf3d_MPI(data, os.path.join(outdir, "ssnr%04d"%total_iter), None,sumvol,last_ring, 1.0, 1, CTF, 1, sym, myid)
				if myid == main_node:   
					varf = 1.0/varf
					varf.write_image( os.path.join(outdir,"varf%04d.hdf"%total_iter) )


		if(myid == main_node):
			refdata = [None]*7
			refdata[0] = numref
			refdata[1] = outdir
			refdata[2] = None
			refdata[3] = total_iter
			refdata[4] = varf
			refdata[5] = mask3D
			refdata[6] = (runtype=="REFINEMENT") # whether align on 50S, this only happens at refinement step
			user_func( refdata )

		#  here we  write header info
		mpi_barrier(MPI_COMM_WORLD)
		start_time = time()
		if runtype=="REFINEMENT":
			par_str = ['xform.projection', 'ID', 'group']
		else:
			par_str = ['group', 'ID' ]
	        if myid == main_node:
			from utilities import file_type
	        	if(file_type(stack) == "bdb"):
	        		from utilities import recv_attr_dict_bdb
	        		recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        	else:
	        		from utilities import recv_attr_dict
	        		recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        else:		send_attr_dict(main_node, data, par_str, image_start, image_end)
		if(terminate == 1):
			if myid==main_node:
				print_end_msg("mref_ali3d_MPI terminated due to small number of objects changing assignments")
			from sys import exit
			exit()
		if(myid == 0):
			print_msg( "Time to write headers: %d\n" % (time()-start_time) );start_time = time()
	if myid==main_node:
		print_end_msg("mref_ali3d_MPI")

def get_refiparams(nx):
	from EMAN2 import Processor
	M = nx
	npad = 2
	N = M*npad
	K = 6
	alpha = 1.75
	r = M/2
	v = K/2.0/N
	return {"filter_type": Processor.fourier_filter_types.KAISER_SINH_INVERSE, "alpha":alpha, "K":K, "r":r, "v":v, "N":N}

def local_ali3dm_MPI_(stack, refvol, outdir, maskfile, ou=-1,  delta=2, ts=0.25, maxit=10, nassign=4, nrefine=1, CTF = None,
                snr=1.0, sym="c1", user_func_name="ref_ali3d", fourvar=False, debug=False, termprec = 0.0 ):
	"""
	  Focus on intersubunit region	
	"""
	from alignment	    import eqproj_cascaded_ccc
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl, filt_tanl, filt_vols
	from fundamentals   import fshift, rot_avg_image
	from projection     import prep_vol, prgs, project
	from utilities      import amoeba_multi_level, bcast_string_to_all, model_circle, get_arb_params, set_arb_params, drop_spider_doc
	from utilities      import bcast_number_to_all, bcast_list_to_all,get_image, drop_image, bcast_EMData_to_all, send_attr_dict
	from utilities      import get_params_proj, set_params_proj, get_im
	from utilities      import model_blank, print_begin_msg, print_msg, print_end_msg, file_type
	from reconstruction import rec3D_MPI
	from statistics     import ccc
	from pixel_error    import max_3D_pixel_error
	from math           import pi, sqrt
	from string         import replace
	from mpi            import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi            import mpi_reduce, MPI_INT, MPI_SUM
	from utilities      import estimate_3D_center_MPI, rotate_3D_shift
	from EMAN2 import Processor
	import os
	import sys
	from EMAN2 import Vec2f


	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "local_ali3dm_MPI_ ", 1,myid)
	mpi_barrier(MPI_COMM_WORLD)
	
	if myid == main_node:
		print_begin_msg("local_ali3dm_MPI")
		os.mkdir(outdir)
	
	mpi_barrier(MPI_COMM_WORLD)



	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		finfo = open(os.path.join(outdir, "progress%04d"%myid), 'w')
	else:
		finfo = None

	from time import time	


        # refine step define on which step refinement will be carried
        # if set to -1, no (??) refinement only assignment 

	nx  = get_image( refvol ).get_xsize()
	ou = int(ou)
	if(ou <= 0):  ou = nx//2-2
	if maskfile:
		import  types
		if(type(maskfile) is types.StringType): 
			mask3D = get_image(maskfile)
		else:   
			mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)

	mask2D  = model_circle(ou, nx, nx)
	fscmask = model_circle(ou, nx, nx, nx)

	numref = EMUtil.get_image_count(refvol)
	if myid==main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]
		for krf in xrange(numref):
			vol = get_im(refvol, krf)
			vol.write_image( os.path.join(outdir, "volf0000.hdf"), krf )
			vol = None
		print_begin_msg("local_ali3dm_MPI")
		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Reference volume            : %s\n"%(refvol))	
		print_msg("Number of reference volumes : %i\n"%(numref))
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Maskfile                    : %s\n"%(maskfile))
		print_msg("Angular step                : %f\n"%(delta))
		print_msg("Shift search range          : %f\n"%(ts))

	if(myid == main_node):
       		if(file_type(stack) == "bdb"):
			from EMAN2db import db_open_dict
			dummy = db_open_dict(stack, True)
		active = EMUtil.get_all_attributes(stack, "active")
		list_of_particles = []
		for im in xrange( len(active) ):
			if( active[im] ) : list_of_particles.append(im)
		del active
		nima = len( list_of_particles )
		start_time = time()
	else:
		nima = 0

	nima = bcast_number_to_all( nima, source_node = main_node )

	if(myid != main_node):
		list_of_particles = [-1]*nima

	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
	# create a list of images for each node
	total_nima = nima
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)

	if debug:
		finfo.write( "image_start, image_end: %d %d\n" % (image_start, image_end) )
		finfo.flush()

	data = EMData.read_images(stack, list_of_particles)	
	if(debug) :
		finfo.write( '%d loaded  \n' % len(data) )	
		finfo.flush()

	#  Initialize Particle ID and set group number to non-existant -1
	assignment = [-1]*len(data)
	for im in xrange(len(data)):
		data[im].set_attr('ID', list_of_particles[im])

	if fourvar:
		#  I am not sure why it is here!  PAP 09/26/09
		from reconstruction import rec3D_MPI
		from statistics     import varf3d_MPI
		#  Compute Fourier variance
		vol, fscc = rec3D_MPI(data, snr, sym, fscmask, os.path.join(outdir, "resolution0000"), myid, main_node, finfo=finfo)
		varf = varf3d_MPI(data, os.path.join(outdir, "ssnr0000"), None, vol, int(ou), 1.0, 1, CTF, 1, sym, myid)
		if myid == main_node:
			varf = 1.0/varf
			varf.write_image( os.path.join(outdir,"varf0000.hdf") )
			print_msg("Time to calculate 3D Fourier variance = %d\n"%(time()-start_time))
			start_time = time()
	else:
		varf = None


	if(CTF):
		if(data[0].get_attr("ctf_applied") > 0):
			ERROR( "local_ali3dm does not work on ctf_applied data", "local_ali3dm_MPI_", 1,myid)
		from reconstruction import rec3D_MPI
	else:
		from reconstruction import rec3D_MPI_noCTF
	

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )



	if(myid == main_node):
		momm = get_im("varmask.hdf")
	else:
		momm =  model_blank(nx, nx, nx)
	bcast_EMData_to_all(momm, myid, main_node)
	from projection import project
	


        refiparams = get_refiparams(nx)

	maxit = maxit*(nassign+nrefine)
	Iter = -1
	iteration = 0
	while(Iter < maxit - 1):
		#for Iter in xrange(maxit):
		Iter += 1
		if Iter%(nassign+nrefine) < nassign :
			runtype = "ASSIGNMENT"
		else:
			runtype = "REFINEMENT"

		iteration += 1
		if(myid == main_node) :
			start_time = time()
			print_msg( runtype + (" ITERATION #%3d\n"%iteration) )

		peaks = [-1.0e23] * nima
		if runtype=="REFINEMENT":  pixer = [0.0]*nima

		for krf in xrange(numref):
			vol = get_im(os.path.join(outdir, "volf%04d.hdf"%(iteration-1)), krf)
			if CTF:
				previous_defocus = -1
			else:
				volft,kb = prep_vol(vol)

			for im in xrange(nima):
				img = data[im]
				if CTF:
					ctf = img.get_attr( "ctf" )
					if ctf.defocus != previous_defocus:
						ctfvol = filt_ctf( vol, ctf )
						volft, kb = prep_vol( ctfvol )
						previous_defocus = ctf.defocus

				phi,tht,psi,s2x,s2y = get_params_proj(img)
				t1 = img.get_attr("xform.projection")
				dp = t1.get_params("spider")
				phi =  dp["phi"]
				tht =  dp["theta"]
				psi =  dp["psi"]
				s2x = -dp["tx"]
				s2y = -dp["ty"]
				if runtype=="ASSIGNMENT":
					pmomm = Util.addn_img(project(momm, [phi,tht,psi,-s2x,-s2y], ou), mask2D)
					ref   = Util.muln_img(prgs( volft, kb, [phi,tht,psi,-s2x,-s2y] ), pmomm)
					peak = ref.cmp("ccc",img,{"mask":mask2D, "negative":0})
					if not(finfo is None):
						finfo.write( "ID,iref,peak: %6d %d %f"%(list_of_particles[im],krf,peak) )
						finfo.flush()
				else:
					refi = img.process( "normalize.mask", {"mask":mask2D, "no_sigma":0} )
					refi = refi.FourInterpol(nx*2,nx*2,0,True)
					refi = Processor.EMFourierFilter(refi, refiparams)
					refdata = [None]*7
					refdata[0] = volft
					refdata[1] = kb
					refdata[2] = img
					refdata[3] = mask2D
					refdata[4] = refi
					refdata[5] = [-s2x,-s2y]
					refdata[6] = ts
					weight_phi = max(delta, delta*abs((tht-90.0)/180.0*pi))
					[ang,peak,qiter,sft] = amoeba_multi_level([phi,tht,psi],[weight_phi,delta,weight_phi],eqproj_cascaded_ccc, 1.0,1.e-2, 500, refdata)
					if not(finfo is None):
						finfo.write( "ID,iref,peak,trans: %6d %d %f %f %f %f %f %f"%(list_of_particles[im],krf,peak,ang[0],ang[1],ang[2],-sft[0],-sft[1]) )
						finfo.flush()

				if(peak > peaks[im]):
					peaks[im] = peak
					data[im].set_attr( "group",  krf )
					if runtype=="REFINEMENT":
						t2 = Transform({"type":"spider","phi":ang[0],"theta":ang[1],"psi":ang[2]})
						t2.set_trans(Vec2f(sft[0], sft[1]))
						data[im].set_attr("xform.projection", t2)
						pixer[im] = max_3D_pixel_error(t1, t2, ou)

					if not(finfo is None):
						finfo.write( " current best" )

				if not(finfo is None):
					finfo.write( "\n" )
				if( myid== main_node and (im>0) and ( ((im)%(nima//2) == 0) or (im == nima-1) ) ):
					print_msg( "Time to process %6d particles : %d\n" % (nima//2, time()-start_time) )
					start_time = time()

		del peaks
		del vol, volft
		#  compute number of particles that changed assignments and how many are in which group
		nchng = 0
		npergroup = [0]*numref
		for im in xrange(nima):
			iref = data[im].get_attr('group')
			npergroup[iref] += 1
			if( iref != assignment[im]):
				assignment[im] = iref
				nchng += 1
		nchng = mpi_reduce(nchng, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
		npergroup = mpi_reduce(npergroup, numref, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
		npergroup = map(int, npergroup)
		terminate = 0
		if( myid == 0 ):
			nchng = int(nchng[0])
			precn = 100*float(nchng)/float(total_nima)
			msg = " Number of particles that changed assignment %7d, percentage of total: %5.1f\n"%(nchng, precn)
			print_msg(msg)
			msg = " Group       number of particles\n"
			print_msg(msg)
			for iref in xrange(numref):
				msg = " %5d       %7d\n"%(iref+1, npergroup[iref])
				print_msg(msg)
			if(precn <= termprec):  terminate = 1
		terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
		terminate = int(terminate[0])


		if runtype=="REFINEMENT":
			if(True):
				cs = [0.0]*3
				cs[0],cs[1],cs[2],dummy,dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)				
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					print_msg(msg)
				cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(data, cs)
			#output pixel errors
			from mpi import mpi_gatherv
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			if(myid == main_node):
				recvbuf = map(float, recvbuf)
				from statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				if(region[0] < 0.0):  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				del region, histo
			del recvbuf


		fscc = [None]*numref
		if fourvar and runtype=="REFINEMENT":
			sumvol = model_blank(nx, nx, nx)
		for krf in xrange(numref):
			if CTF:
				vol, fscc[krf] = rec3D_MPI(data, snr, sym, fscmask, os.path.join(outdir, "resolution%02d_%04d"%(krf, iteration)), myid, main_node, index = krf)
			else:
				vol, fscc[krf] = rec3D_MPI_noCTF(data, snr, sym, fscmask, os.path.join(outdir, "resolution%02d_%04d"%(krf, iteration)), myid, main_node, index = krf)
				

			if(myid==main_node):
				vol.write_image(os.path.join(outdir,"vol%04d.hdf"%iteration),krf)
				print_msg("3D reconstruction time = %d\n"%(time()-start_time))
				start_time = time()
				if fourvar and runtype=="REFINEMENT":
					sumvol += vol

		if runtype=="REFINEMENT":
			if fourvar:
				varf = varf3d_MPI(data, os.path.join(outdir, "ssnr%04d"%iteration), None, sumvol, int(ou), 1.0, 1, CTF, 1, sym, myid)
				if myid == main_node:   
					varf = 1.0/varf
					varf.write_image( os.path.join(outdir,"varf%04d.hdf"%iteration) )
					print_msg("Time to calculate 3D Fourier variance= %d\n"%(time()-start_time))
					start_time = time()


		if(myid == main_node):
			refdata = [None]*7
			refdata[0] = numref
			refdata[1] = outdir
			refdata[2] = None # fscc
			refdata[3] = iteration
			refdata[4] = varf
			refdata[5] = mask3D
			refdata[6] = (runtype=="REFINEMENT") # whether align on 50S, this only happens at refinement step
			user_func( refdata )


		mpi_barrier(MPI_COMM_WORLD)
		# here we should write header info, just in case the program crashes...
		# write out headers  and STOP, under MPI writing has to be done sequentially

		if runtype=="REFINEMENT":
			par_str = ["xform.projection", "group", "ID"]
		else:
			par_str = ["group", "ID"]

	        if myid == main_node:
			from utilities import file_type
	        	if(file_type(stack) == "bdb"):
	        		from utilities import recv_attr_dict_bdb
	        		recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        	else:
	        		from utilities import recv_attr_dict
	        		recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        else:		send_attr_dict(main_node, data, par_str, image_start, image_end)
		if myid == main_node:
			print_msg("Time to write header information= %d\n"%(time()-start_time))
			start_time = time()
		if(terminate == 1  and runtype=="ASSIGNMENT"):
			if myid==main_node:
				#print_end_msg("local_ali3dm_MPI terminated due to small number of objects changing assignments")
				print_msg("local_ali3dm_MPI abandoned assignments due to small number of objects changing assignments\n")
			from sys import exit
			exit()
			while(runtype == "ASSIGNMENT"):
				Iter += 1
				if Iter%(nassign+nrefine) < nassign :
					runtype = "ASSIGNMENT"
				else:
					runtype = "REFINEMENT"
			Iter += -1

	if myid==main_node:
		print_end_msg("local_ali3dm_MPI")


def local_ali3dm_MPI(stack, refvol, outdir, maskfile, ou=-1,  delta=2, ts=0.25, maxit=10, nassign=4, nrefine=1, CTF = None,
                snr=1.0, sym="c1", user_func_name="ref_ali3d", fourvar=False, npad = 4, debug=False, termprec = 0.0 ):
	"""
	  The original, fully operational version	
	"""
	from alignment	    import eqproj_cascaded_ccc
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl, filt_tanl, filt_vols
	from fundamentals   import fshift, rot_avg_image
	from projection     import prep_vol, prgs, project
	from utilities      import amoeba_multi_level, bcast_string_to_all, model_circle, get_arb_params, set_arb_params, drop_spider_doc
	from utilities      import bcast_number_to_all, bcast_list_to_all,get_image, drop_image, bcast_EMData_to_all, send_attr_dict
	from utilities      import get_params_proj, set_params_proj, get_im
	from utilities      import model_blank, print_begin_msg, print_msg, print_end_msg, file_type
	from reconstruction import rec3D_MPI
	from statistics     import ccc
	from pixel_error    import max_3D_pixel_error
	from math           import pi, sqrt
	from string         import replace
	from mpi            import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi            import mpi_reduce, MPI_INT, MPI_SUM
	from utilities      import estimate_3D_center_MPI, rotate_3D_shift
	from EMAN2 import Processor
	import os
	import sys
	from EMAN2 import Vec2f


	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "local_ali3dm_MPI ", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if(myid == main_node):
		print_begin_msg("local_ali3dm_MPI")
		os.mkdir(outdir)
	
	mpi_barrier(MPI_COMM_WORLD)



	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		finfo = open(os.path.join(outdir, "progress%04d"%myid), 'w')
	else:
		finfo = None

	from time import time	


        # refine step define on which step refinement will be carried
        # if set to -1, no (??) refinement only assignment 

	nx  = get_image( refvol ).get_xsize()
	ou = int(ou)
	if(ou <= 0):  ou = nx//2-2
	if maskfile:
		import  types
		if(type(maskfile) is types.StringType): 
			mask3D = get_image(maskfile)
		else:   
			mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)

	mask2D  = model_circle(ou, nx, nx)
	fscmask = model_circle(ou, nx, nx, nx)

	numref = EMUtil.get_image_count(refvol)
	if myid==main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]
		for krf in xrange(numref):
			vol = get_im(refvol, krf)
			vol.write_image( os.path.join(outdir, "volf0000.hdf"), krf )
			vol = None
		print_begin_msg("local_ali3dm_MPI")
		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Reference volume            : %s\n"%(refvol))	
		print_msg("Number of reference volumes : %i\n"%(numref))
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Maskfile                    : %s\n"%(maskfile))
		print_msg("Angular step                : %f\n"%(delta))
		print_msg("Shift search range          : %f\n"%(ts))

	if(myid == main_node):
       		if(file_type(stack) == "bdb"):
			from EMAN2db import db_open_dict
			dummy = db_open_dict(stack, True)
		active = EMUtil.get_all_attributes(stack, "active")
		list_of_particles = []
		for im in xrange( len(active) ):
			if( active[im] ) : list_of_particles.append(im)
		del active
		nima = len( list_of_particles )
		start_time = time()
	else:
		nima = 0

	nima = bcast_number_to_all( nima, source_node = main_node )

	if(myid != main_node):
		list_of_particles = [-1]*nima

	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
	# create a list of images for each node
	total_nima = nima
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)

	if debug:
		finfo.write( "image_start, image_end: %d %d\n" % (image_start, image_end) )
		finfo.flush()

	data = EMData.read_images(stack, list_of_particles)	
	if(debug) :
		finfo.write( '%d loaded  \n' % len(data) )	
		finfo.flush()

	#  Initialize Particle ID and set group number to non-existant -1
	assignment = [-1]*len(data)
	for im in xrange(len(data)):
		data[im].set_attr('ID', list_of_particles[im])

	if fourvar:
		#  I am not sure why it is here!  PAP 09/26/09
		from reconstruction import rec3D_MPI
		from statistics     import varf3d_MPI
		#  Compute Fourier variance
		vol, fscc = rec3D_MPI(data, snr, sym, fscmask, os.path.join(outdir, "resolution0000"), myid, main_node, finfo=finfo, npad = npad)
		varf = varf3d_MPI(data, os.path.join(outdir, "ssnr0000"), None, vol, int(ou), 1.0, 1, CTF, 1, sym, myid)
		if myid == main_node:
			varf = 1.0/varf
			varf.write_image( os.path.join(outdir,"varf0000.hdf") )
			print_msg("Time to calculate 3D Fourier variance = %d\n"%(time()-start_time))
			start_time = time()
	else:
		varf = None


	if(CTF):
		if(data[0].get_attr("ctf_applied") > 0):
			ERROR( "local_ali3dm does not work on ctf_applied data", "local_ali3dm_MPI", 1,myid)
		from reconstruction import rec3D_MPI
	else:
		from reconstruction import rec3D_MPI_noCTF
	

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )
	
        refiparams = get_refiparams(nx)

	maxit = maxit*(nassign+nrefine)
	Iter = -1
	iteration = 0
	while(Iter < maxit - 1):
		#for Iter in xrange(maxit):
		Iter += 1
		if Iter%(nassign+nrefine) < nassign :
			runtype = "ASSIGNMENT"
		else:
			runtype = "REFINEMENT"

		iteration += 1
		if(myid == main_node) :
			start_time = time()
			print_msg( runtype + (" ITERATION #%3d\n"%iteration) )

		peaks = [-1.0e23] * nima
		if runtype=="REFINEMENT":  pixer = [0.0]*nima

		for krf in xrange(numref):
			vol = get_im(os.path.join(outdir, "volf%04d.hdf"%(iteration-1)), krf)
			if CTF:
				previous_defocus = -1
			else:
				volft,kb = prep_vol(vol)

			for im in xrange(nima):
				img = data[im]
				if CTF:
					ctf = img.get_attr( "ctf" )
					if ctf.defocus != previous_defocus:
						ctfvol = filt_ctf( vol, ctf )
						volft, kb = prep_vol( ctfvol )
						previous_defocus = ctf.defocus

				phi,tht,psi,s2x,s2y = get_params_proj(img)
				t1 = img.get_attr("xform.projection")
				dp = t1.get_params("spider")
				phi =  dp["phi"]
				tht =  dp["theta"]
				psi =  dp["psi"]
				s2x = -dp["tx"]
				s2y = -dp["ty"]
				if runtype=="ASSIGNMENT":
					ref  = prgs( volft, kb, [phi,tht,psi,-s2x,-s2y] )
					peak = ref.cmp("ccc",img,{"mask":mask2D, "negative":0})
					if not(finfo is None):
						finfo.write( "ID,iref,peak: %6d %d %f"%(list_of_particles[im],krf,peak) )
						finfo.flush()
				else:
					refi = img.process( "normalize.mask", {"mask":mask2D, "no_sigma":0} )
					refi = refi.FourInterpol(nx*2,nx*2,0,True)
					refi = Processor.EMFourierFilter(refi, refiparams)
					refdata = [None]*7
					refdata[0] = volft
					refdata[1] = kb
					refdata[2] = img
					refdata[3] = mask2D
					refdata[4] = refi
					refdata[5] = [-s2x,-s2y]
					refdata[6] = ts
					weight_phi = max(delta, delta*abs((tht-90.0)/180.0*pi))
					[ang,peak,qiter,sft] = amoeba_multi_level([phi,tht,psi],[weight_phi,delta,weight_phi],eqproj_cascaded_ccc, 1.0,1.e-2, 500, refdata)
					if not(finfo is None):
						finfo.write( "ID,iref,peak,trans: %6d %d %f %f %f %f %f %f"%(list_of_particles[im],krf,peak,ang[0],ang[1],ang[2],-sft[0],-sft[1]) )
						finfo.flush()

				if(peak > peaks[im]):
					peaks[im] = peak
					data[im].set_attr( "group",  krf )
					if runtype=="REFINEMENT":
						t2 = Transform({"type":"spider","phi":ang[0],"theta":ang[1],"psi":ang[2]})
						t2.set_trans(Vec2f(sft[0], sft[1]))
						data[im].set_attr("xform.projection", t2)
						pixer[im] = max_3D_pixel_error(t1, t2, ou)

					if not(finfo is None):
						finfo.write( " current best" )

				if not(finfo is None):
					finfo.write( "\n" )
				if( myid== main_node and (im>0) and ( ((im)%(nima//2) == 0) or (im == nima-1) ) ):
					print_msg( "Time to process %6d particles : %d\n" % (nima//2, time()-start_time) )
					start_time = time()

		del peaks
		del vol, volft
		#  compute number of particles that changed assignments and how many are in which group
		nchng = 0
		npergroup = [0]*numref
		for im in xrange(nima):
			iref = data[im].get_attr('group')
			npergroup[iref] += 1
			if( iref != assignment[im]):
				assignment[im] = iref
				nchng += 1
		nchng = mpi_reduce(nchng, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
		npergroup = mpi_reduce(npergroup, numref, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
		npergroup = map(int, npergroup)
		terminate = 0
		if( myid == 0 ):
			nchng = int(nchng[0])
			precn = 100*float(nchng)/float(total_nima)
			msg = " Number of particles that changed assignment %7d, percentage of total: %5.1f\n"%(nchng, precn)
			print_msg(msg)
			msg = " Group       number of particles\n"
			print_msg(msg)
			for iref in xrange(numref):
				msg = " %5d       %7d\n"%(iref+1, npergroup[iref])
				print_msg(msg)
			if(precn <= termprec):  terminate = 1
		terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
		terminate = int(terminate[0])


		if runtype=="REFINEMENT":
			if(True):
				cs = [0.0]*3
				cs[0],cs[1],cs[2],dummy,dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)				
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					print_msg(msg)
				cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(data, cs)
			#output pixel errors
			from mpi import mpi_gatherv
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			if(myid == main_node):
				recvbuf = map(float, recvbuf)
				from statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				if(region[0] < 0.0):  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				del region, histo
			del recvbuf


		fscc = [None]*numref
		if fourvar and runtype=="REFINEMENT":
			sumvol = model_blank(nx, nx, nx)
		for krf in xrange(numref):
			if CTF:
				vol, fscc[krf] = rec3D_MPI(data, snr, sym, fscmask, os.path.join(outdir, "resolution%02d_%04d"%(krf, iteration)), myid, main_node, index = krf, npad = npad)
			else:
				vol, fscc[krf] = rec3D_MPI_noCTF(data, snr, sym, fscmask, os.path.join(outdir, "resolution%02d_%04d"%(krf, iteration)), myid, main_node, index = krf, npad = npad)
				

			if(myid==main_node):
				vol.write_image(os.path.join(outdir,"vol%04d.hdf"%iteration),krf)
				print_msg("3D reconstruction time = %d\n"%(time()-start_time))
				start_time = time()
				if fourvar and runtype=="REFINEMENT":
					sumvol += vol

		if runtype=="REFINEMENT":
			if fourvar:
				varf = varf3d_MPI(data, os.path.join(outdir, "ssnr%04d"%iteration), None, sumvol, int(ou), 1.0, 1, CTF, 1, sym, myid)
				if myid == main_node:   
					varf = 1.0/varf
					varf.write_image( os.path.join(outdir,"varf%04d.hdf"%iteration) )
					print_msg("Time to calculate 3D Fourier variance= %d\n"%(time()-start_time))
					start_time = time()


		if(myid == main_node):
			refdata = [None]*7
			refdata[0] = numref
			refdata[1] = outdir
			refdata[2] = None # fscc
			refdata[3] = iteration
			refdata[4] = varf
			refdata[5] = mask3D
			refdata[6] = (runtype=="REFINEMENT") # whether align on 50S, this only happens at refinement step
			user_func( refdata )


		mpi_barrier(MPI_COMM_WORLD)
		# here we should write header info, just in case the program crashes...
		# write out headers  and STOP, under MPI writing has to be done sequentially

		if runtype=="REFINEMENT":
			par_str = ["xform.projection", "group", "ID"]
		else:
			par_str = ["group", "ID"]

	        if myid == main_node:
			from utilities import file_type
	        	if(file_type(stack) == "bdb"):
	        		from utilities import recv_attr_dict_bdb
	        		recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        	else:
	        		from utilities import recv_attr_dict
	        		recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        else:		send_attr_dict(main_node, data, par_str, image_start, image_end)
		if myid == main_node:
			print_msg("Time to write header information= %d\n"%(time()-start_time))
			start_time = time()
		if(terminate == 1  and runtype=="ASSIGNMENT"):
			if myid==main_node:
				#print_end_msg("local_ali3dm_MPI terminated due to small number of objects changing assignments")
				print_msg("local_ali3dm_MPI abandoned assignments due to small number of objects changing assignments\n")

			while(runtype == "ASSIGNMENT"):
				Iter += 1
				if Iter%(nassign+nrefine) < nassign :
					runtype = "ASSIGNMENT"
				else:
					runtype = "REFINEMENT"
			Iter += -1

	if myid==main_node:
		print_end_msg("local_ali3dm_MPI")


def local_ali3d(stack, outdir, maskfile = None, ou = -1,  delta = 2, ts=0.25, center = -1, maxit = 10, 
           CTF = False, snr = 1.0, sym = "c1", chunk = -1.0, user_func_name = "ref_ali3d",
	     fourvar = True, npad = 4, debug = False, MPI = False):
	"""
		
	"""

	if MPI:
		local_ali3d_MPI(stack, outdir, maskfile, ou, delta, ts, center, maxit,
				CTF, snr, sym, chunk, user_func_name, 
				fourvar, npad, debug)
		return

	from alignment      import eqproj_cascaded_ccc
	from projection     import prep_vol
	from utilities      import model_circle, get_params_proj, set_params_proj
	from utilities      import get_image, drop_image
	from utilities      import amoeba_multi_level, rotate_3D_shift, estimate_3D_center
	from math           import pi
	from statistics     import fsc_mask
	from utilities      import print_begin_msg, print_end_msg, print_msg
	from EMAN2 import Processor
	import os 
	import sys

	if os.path.exists(outdir): ERROR('Output directory exists, please change the name and restart the program', "local_ali3d", 1)
	os.mkdir(outdir)
	
	print_begin_msg('local_ali3d')

	import user_functions
	user_func = user_functions.factory[user_func_name]

	if CTF:
		ima = EMData()
		ima.read_image(stack, 0)
		ctf_applied = ima.get_attr("ctf_applied")
		del ima
		if ctf_applied == 1:  ERROR("local_ali3d does not work for CTF-applied data", "local_ali3d", 1)
		from reconstruction import recons3d_4nn_ctf
		from filter         import filt_ctf
	else   : from reconstruction import recons3d_4nn

	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	ima     = EMData()
	ima.read_image(stack, 0)
	nx      = ima.get_xsize()
	del ima
	if last_ring == -1:	last_ring = nx//2 - 2

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Maskfile                    : %s\n"%(maskfile))
	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Angular search range        : %s\n"%(delta))
	print_msg("Translation search range    : %f\n"%(ts))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("CTF correction              : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	print_msg("Symmetry group              : %s\n"%(sym))
	if chunk <= 0.0:  chunk = 1.0
	print_msg("Chunk size                  : %f\n\n"%(chunk))
	
	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else:
		mask3D = model_circle(last_ring, nx, nx, nx)
	mask2D = model_circle(last_ring, nx, nx)


	if debug:  outf = file(os.path.join(outdir, "progress"), "w")
	else:      outf = None

	active = EMUtil.get_all_attributes(stack, 'active')
	list_of_particles = []
	for im in xrange(len(active)):
		if(active[im]):  list_of_particles.append(im)
	del active
	dataim = EMData.read_images(stack, list_of_particles)
	nima = len(dataim)

	if debug:
		outf.write("  data read")
		outf.write("\n")
		outf.flush()

	n_in_chunk  = max(int(chunk * nima), 1)
	n_of_chunks = nima//n_in_chunk + min(nima%n_in_chunk, 1)
	
	if debug:
		outf = file(os.path.join(outdir, "progress"), "w")
		outf.write("  chunk = "+str(chunk)+"   ")
		outf.write("\n")
		outf.flush()
		outf.write("  chunk = "+str(n_in_chunk)+"   ")
		outf.write("  chunk = "+str(n_of_chunks)+"   ")
		outf.write("\n")
		outf.flush()

	# initialize data for the reference preparation function
	ref_data = [mask3D, center, None, None]

	M = nx
	npad = 2
	N = M*npad
	K = 6
	alpha = 1.75
	r = M/2
	v = K/2.0/N
	params = {"filter_type": Processor.fourier_filter_types.KAISER_SINH_INVERSE, "alpha":alpha, "K":K, "r":r, "v":v, "N":N}

	data = [None]*7
	data[3] = mask2D
	cs = [0.0]*3

	for iteration in xrange(maxit+1):
		print_msg("ITERATION #%3d\n"%(iteration+1))
		for ic in xrange(n_of_chunks):
			if(center == -1):
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center(dataim)				
				rotate_3D_shift(dataim, [-cs[0], -cs[1], -cs[2]])
				msg = "Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
				print_msg(msg)				
			# compute updated 3D at the beginning of each chunk
			#  3D stuff
			if CTF: vol1 = recons3d_4nn_ctf(dataim, range(0, nima, 2), snr, 1, sym)
			else:	   vol1 = recons3d_4nn(dataim, range(0, nima, 2), sym)

			if CTF: vol2 = recons3d_4nn_ctf(dataim, range(1, nima, 2), snr, 1, sym)
			else:	   vol2 = recons3d_4nn(dataim, range(1, nima, 2), sym)

	    		# resolution
			fscc = fsc_mask(vol1, vol2, mask3D, 1.0, os.path.join(outdir, "resolution%04d"%(iteration*n_of_chunks+ic+1)))
			del vol1
			del vol2

			# calculate new and improved 3D
			if CTF: vol = recons3d_4nn_ctf(dataim, range(nima), snr, 1, sym)
			else:	   vol = recons3d_4nn(dataim, range(nima), sym)

			# store the reference volume
			drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(iteration*n_of_chunks+ic+1)))
			ref_data[2] = vol
			ref_data[3] = fscc

			#  call user-supplied function to prepare reference image, i.e., center and filter it
			vol, dummy = user_func(ref_data)

			drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(iteration*n_of_chunks+ic+1)))
			if(iteration == maxit):
				#  in last iteration quit here
				print_end_msg("local_ali3d")
				return

			if not CTF:
				data[0], data[1] = prep_vol(vol)

			image_start_in_chunk = ic*n_in_chunk
			image_end_in_chunk   = min(image_start_in_chunk + n_in_chunk, nima)
			if debug:
				outf.write("image_start_in_chunk "+str(image_start_in_chunk)+"\n")
				outf.write("\n")
				outf.write("image_end_in_chunk "+str(image_end_in_chunk)+"\n")
				outf.write("\n")
				outf.flush()
			if CTF:  previous_defocus = -1.0
			for imn in xrange(image_start_in_chunk, image_end_in_chunk):
				if CTF:
					ctf_params = dataim[imn].get_attr( "ctf" )
					if ctf_params.defocus != previous_defocus:
						previous_defocus = ctf_params.defocus
						data[0], data[1] = prep_vol(filt_ctf(vol, ctf_params))

				data[2] = dataim[imn]

				if ts > 0.0:
					refi = dataim[imn].FourInterpol(nx*2, nx*2, 0, False)
					data[4] = Processor.EMFourierFilter(refi, params)
				
				phi, theta, psi, tx, ty = get_params_proj(dataim[imn])
				atparams = [phi, theta, psi]
				data[5] = [tx, ty]
				data[6] = ts
				data[5][0] *= -1
				data[5][1] *= -1

				if debug:
					initial, dummy  = eqproj_cascaded_ccc(atparams, data)  # this is if we need initial discrepancy
					outf.write("Image "+str(imn)+"\n")
					outf.write('Old  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %11.4f'%(phi,theta,psi,tx,ty,initial))
					outf.write("\n")
					
				# change signs of shifts for projections
			
				weight_phi = max(delta, delta*abs((atparams[1]-90.0)/180.0*pi))

				optm_params = amoeba_multi_level(atparams, [weight_phi, delta, weight_phi], eqproj_cascaded_ccc, 1.e-4, 1.e-4, 500, data)
				optm_params[0].append(optm_params[3][0])
				optm_params[0].append(optm_params[3][1])
				optm_params[0][3] *= -1
				optm_params[0][4] *= -1
				
				if debug:
					outf.write('New  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %11.4f  %4d'%(optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4], optm_params[1], optm_params[2]))
					outf.write("\n")
					outf.flush()

				set_params_proj(dataim[imn], optm_params[0])

			#  here we write header infomation
			from utilities import write_headers
			#write_headers(stack, dataim, list_of_particles)


def local_ali3d_MPI(stack, outdir, maskfile, ou = -1,  delta = 2, ts=0.25, center = -1, maxit = 10, 
                CTF = False, snr = 1.0, sym = "c1", chunk = -1.0, user_func_name = "ref_ali3d",
		    fourvar = True, npad = 4, debug = False):
	"""
		
	"""
	from alignment        import eqproj_cascaded_ccc
	from filter           import filt_ctf
	from projection       import prep_vol
	from utilities        import bcast_string_to_all, bcast_number_to_all, model_circle, get_params_proj, set_params_proj
	from utilities        import bcast_EMData_to_all, bcast_list_to_all, send_attr_dict
	from utilities        import get_image, drop_image, file_type
	from utilities        import amoeba_multi_level, rotate_3D_shift, estimate_3D_center_MPI
	from utilities        import print_begin_msg, print_end_msg, print_msg
	from reconstruction   import rec3D_MPI, rec3D_MPI_noCTF
	from statistics       import varf3d_MPI
	from math             import pi
	from mpi              import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi              import mpi_reduce, MPI_INT, MPI_SUM
	from EMAN2 import Processor
	import os
	import sys
	from EMAN2 import Vec2f


	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)

	if CTF:
		from filter import filt_ctf

	main_node = 0
	
	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "local_ali3d_MPI ", 1,myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		print_begin_msg("local_ali3d_MPI")
		os.mkdir(outdir)
		import user_functions
		user_func = user_functions.factory[user_func_name]
		if CTF:
			ima = EMData()
			ima.read_image(stack, 0)
			ctf_applied = ima.get_attr("ctf_applied")
			del ima
			if ctf_applied == 1:  ERROR("local_ali3d does not work for CTF-applied data", "local_ali3d_MPI", 1,myid)
	mpi_barrier(MPI_COMM_WORLD)

	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		finfo = open(os.path.join(outdir, "progress%04d"%myid), 'w')
	else:
		finfo = None

	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	if myid == main_node:
       		if(file_type(stack) == "bdb"):
			from EMAN2db import db_open_dict
			dummy = db_open_dict(stack, True)
		active = EMUtil.get_all_attributes(stack, 'active')
		list_of_particles = []
		for im in xrange(len(active)):
			if(active[im]):  list_of_particles.append(im)
		del active
		nima = len(list_of_particles)
		ima     = EMData()
		ima.read_image(stack, 0)
		nx      = ima.get_xsize()
		del ima
	else:
		nima = 0
		nx = 0
	total_nima = bcast_number_to_all(nima, source_node = main_node)
	nx = bcast_number_to_all(nx, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)

	if last_ring < 0:	last_ring = int(nx/2) - 2

	if chunk <= 0.0:  chunk = 1.0

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_begin_msg("local_ali3d_MPI")
		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Maskfile                    : %s\n"%(maskfile))
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Angular search range        : %s\n"%(delta))
		print_msg("Shift search range          : %f\n"%(ts))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Symmetry group              : %s\n"%(sym))
		print_msg("Chunk size                  : %f\n\n"%(chunk))

	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else:
		mask3D = model_circle(last_ring, nx, nx, nx)
	mask2D = model_circle(last_ring, nx, nx)

	if debug:
		finfo.write( "image_start, image_end: %d %d\n" %(image_start, image_end) )
		finfo.flush()

	n_in_chunk  = max(int(chunk*(image_end-image_start)), 1)
	n_of_chunks = (image_end-image_start)//n_in_chunk + min((image_end-image_start)%n_in_chunk, 1)

	if debug:
		finfo.write("  chunk = "+str(chunk)+"   ")
		finfo.write("\n")
		finfo.flush()
		finfo.write("  Number of images in a chunk = "+str(n_in_chunk)+"   ")
		finfo.write("  Number of chunks = "+str(n_of_chunks)+"   ")
		finfo.write("\n")
		finfo.flush()

	dataim = EMData.read_images(stack, list_of_particles)
	for im in xrange(len(dataim)):
		dataim[im].set_attr('ID', list_of_particles[im])
	del list_of_particles

	if debug:
		finfo.write("  First image on this processor: "+str(image_start)+"   ")
		finfo.write("  Last  image on this processor: "+str(image_end)+"   ")
		finfo.write("\n")
		finfo.flush()

	if myid == main_node:
		# initialize data for the reference preparation function
		ref_data = [ mask3D, max(center,0), None, None, None ]
		# for method -1, switch off centering in user function
		ref_data.append( None )

	from time import time	
		
	M = nx
	npad = 2
	N = M*npad
	K = 6
	alpha = 1.75
	r = M/2
	v = K/2.0/N
	params = {"filter_type": Processor.fourier_filter_types.KAISER_SINH_INVERSE, "alpha":alpha, "K":K, "r":r, "v":v, "N":N}


	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )

	#  this is needed for gathering of pixel errors
	pixer = [0.0]*nima
	data = [None]*7
	data[3] = mask2D
	cs = [0.0]*3

	for iteration in xrange(maxit+1):
		if myid == main_node:
			start_time = time()
			print_msg("ITERATION #%3d\n"%(iteration+1))
		if debug:
			finfo.write("  iteration = "+str(iteration)+"   ")
			finfo.write("\n")
			finfo.flush()
		for ic in xrange(n_of_chunks):
			if(center == -1):
				if debug:
					finfo.write("  begin centering \n")
					finfo.flush()
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(dataim, total_nima, myid, number_of_proc, main_node)
				cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(dataim, cs)
				if myid == main_node:
					msg = "Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					print_msg(msg)		
					print_msg("Time to center = %d\n"%(time()-start_time))
					start_time = time()
			# compute updated 3D before each chunk
 	    		# resolution
			if debug:
				finfo.write("  begin reconstruction = "+str(image_start))
				finfo.write("\n")
				finfo.flush()

			if CTF: vol, fscc = rec3D_MPI(dataim, snr, sym, mask3D, os.path.join(outdir, "resolution%03d_%03d"%(iteration, ic)), myid, main_node, npad = npad)
			else:   vol, fscc = rec3D_MPI_noCTF(dataim, sym, mask3D, os.path.join(outdir, "resolution%03d_%03d"%(iteration, ic)), myid, main_node, npad = npad)

			if myid == main_node:
				drop_image(vol, os.path.join(outdir, "vol%03d_%03d.hdf"%(iteration, ic) ))
				print_msg("3D reconstruction time = %d\n"%(time()-start_time))
				start_time = time()
			if debug:
				finfo.write("  done reconstruction = "+str(image_start))
				finfo.write("\n")
				finfo.flush()

			if fourvar:
			#  Compute Fourier variance
				varf = varf3d_MPI(dataim, ssnr_text_file = os.path.join(outdir, "ssnr%03d_%03d"%(iteration, ic)), mask2D = None, reference_structure = vol, ou = ou, rw = 1.0, npad = 1, CTF = CTF, sign = 1, sym =sym, myid = myid)
				if myid == main_node:
					varf = 1.0/varf
					print_msg("Time to calculate 3D Fourier variance= %d\n"%(time()-start_time))
					start_time = time()

			else:  varf = None
			if myid == main_node:
				ref_data[2] = vol
				ref_data[3] = fscc
				ref_data[4] = varf
				#  call user-supplied function to prepare reference image, i.e., filter it
				# When center = -1, which is by default, we use the average center method
				ref_data[1] = 0
				vol, dummy = user_func(ref_data)
				drop_image(vol, os.path.join(outdir, "volf%03d_%03d.hdf"%(iteration, ic)))
			del varf
			
			# in last iteration return here
			if(iteration == maxit):
				if myid == main_node: print_end_msg("local_ali3d_MPI")
				return
			bcast_EMData_to_all(vol, myid, main_node)

			if not CTF:
				data[0], data[1] = prep_vol(vol)

			image_start_in_chunk = image_start + ic*n_in_chunk
			image_end_in_chunk   = min(image_start_in_chunk + n_in_chunk, image_end)
			if debug:
				finfo.write("Chunk "+str(ic)+"   Number of images in this chunk: "+str(n_in_chunk)+"\n")
				finfo.write("First image in this chunk: "+str(image_start_in_chunk)+"   Last image in this chunk: "+str(image_end_in_chunk-1)+"\n")
				finfo.flush()
			if CTF:  previous_defocus = -1.0
			for imn in xrange(image_start_in_chunk, image_end_in_chunk):
				if CTF:
					ctf_params = dataim[imn-image_start].get_attr( "ctf" )
					if ctf_params.defocus != previous_defocus:
						previous_defocus = ctf_params.defocus
						data[0], data[1] = prep_vol(filt_ctf(vol, ctf_params))

				data[2] = dataim[imn-image_start]
				if ts > 0.0:
					refi = dataim[imn-image_start].FourInterpol(nx*2, nx*2, 0, True)
					data[4] = Processor.EMFourierFilter(refi, params)

				#phi, theta, psi, tx, ty = get_params_proj(dataim[imn-image_start])
				t1 = dataim[imn-image_start].get_attr("xform.projection")
				dp = t1.get_params("spider")
				atparams = [dp["phi"], dp["theta"], dp["psi"]]
				data[5] = [dp["tx"], dp["ty"]]
				if debug:
					# we have to distiguish between no shift situation, which is done through ccc, and shift, which is done using gridding in 2D
					if(ts == 0.0):  data[6] = 0.0
					else:  data[6] = -1.0#ts#-1.0
					initial, dummy = eqproj_cascaded_ccc(atparams, data)  # this is if we need initial discrepancy
					finfo.write("Image "+str(imn)+"\n")
					finfo.write('Old  %6.1f  %6.1f  %6.1f   %5.2f  %5.2f  %11.4e\n'%(atparams[0],atparams[1],atparams[2], -dummy[0], -dummy[1], initial))
				# change signs of shifts for projections
				data[6] = ts
				#from random import random
				#data[5] = [(random()-0.5)*2,(random()-0.5)*2]  #  HERE !!!!!!!!!!!

				weight_phi = max(delta, delta*abs((atparams[1]-90.0)/180.0*pi))

				optm_params = amoeba_multi_level(atparams, [weight_phi, delta, weight_phi], eqproj_cascaded_ccc, 1.0, 1.e-2, 500, data)
				optm_params[0].append(optm_params[3][0])
				optm_params[0].append(optm_params[3][1])
				optm_params[0][3] *= -1
				optm_params[0][4] *= -1

				if debug:
					finfo.write('New  %6.1f  %6.1f  %6.1f   %5.2f  %5.2f  %11.4e  %4d\n'%(optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4], optm_params[1], optm_params[2]))
					finfo.flush()

				#from sys import exit
				#exit()
				t2 = Transform({"type":"spider","phi":optm_params[0][0],"theta":optm_params[0][1],"psi":optm_params[0][2]})
				t2.set_trans(Vec2f(-optm_params[0][3], -optm_params[0][4]))
				dataim[imn-image_start].set_attr("xform.projection", t2)
				from pixel_error import max_3D_pixel_error
				pixer[imn-image_start] = max_3D_pixel_error(t1, t2, last_ring)
				#set_params_proj(dataim[imn-image_start], optm_params[0])
				if( myid== main_node and (imn>image_start_in_chunk) and ( ((imn-image_start_in_chunk)%(n_in_chunk//2) == 0) or (imn == image_end_in_chunk-1) ) ):
					print_msg( "Time to process %6d particles : %d\n" % (n_in_chunk//2, time()-start_time) )
					start_time = time()

			# release memory of volft
			data[0] = None

			# write out headers, under MPI writing has to be done sequentially
			mpi_barrier(MPI_COMM_WORLD)
			par_str = ['xform.projection', 'ID']
			if myid == main_node:
				from utilities import file_type
				if(file_type(stack) == "bdb"):
					from utilities import recv_attr_dict_bdb
					recv_attr_dict_bdb(main_node, stack, dataim, par_str, image_start, image_end, number_of_proc)
				else:
					from utilities import recv_attr_dict
					recv_attr_dict(main_node, stack, dataim, par_str, image_start, image_end, number_of_proc)
			else:	        send_attr_dict(main_node, dataim, par_str, image_start, image_end)
			if myid == main_node:
				print_msg("Time to write header information= %d\n"%(time()-start_time))
				start_time = time()

		#output pixel errors after all headers were processed
		from mpi import mpi_gatherv
		recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
		mpi_barrier(MPI_COMM_WORLD)
		terminate = 0
		if(myid == main_node):
			recvbuf = map(float, recvbuf)
			from statistics import hist_list
			lhist = 20
			region, histo = hist_list(recvbuf, lhist)
			if(region[0] < 0.0):  region[0] = 0.0
			msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
			print_msg(msg)
			for lhx in xrange(lhist):
				msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
				print_msg(msg)
			# Terminate if 95% within 1 pixel error
			im = 0
			for lhx in xrange(lhist):
				if(region[lhx] > 1.0): break
				im += histo[lhx]
			if(im/float(total_nima) > 0.95):  terminate = 1
			del region, histo
		del recvbuf
		terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
		terminate = int(terminate[0])

def autowin(indir,outdir, noisedoc, noisemic, templatefile, deci, CC_method, p_size, sigma, hf_p, n_peak_max, contrast_invert=1, CTF = False, prm = "micrograph", MPI=False):
	""" 
		automated particle picking program
	   	CC_method=1  Applying high pass filter to detect protein particles
		sigma can be adjusted to perform particle picking. 
		CC_method=2  Using fast local normalization method to detect protein particles
		The output files are detected coordinates and windowed particles
	"""
	if MPI:
		autowin_MPI(indir, outdir, noisedoc, noisemic, templatefile, deci, CC_method, p_size, sigma, hf_p, n_peak_max, contract_invert, prm)
		return
	from utilities 		import get_image, model_circle, ce_fit, drop_image,info
	from utilities          import get_arb_params, set_arb_params
	from fundamentals 	import smallprime, window2d, ccf, ramp, fft
	from filter 		import filt_gaussh, filt_tanl
	from string 		import split
	from morphology 	import flcc
	import os
	if os.path.exists(indir)  is False: ERROR("micrograph directory does not exsit", "autowin",1)
	else                              : flist=os.listdir(indir)
	if os.path.exists(outdir)         : ERROR('Output directory exists, please change the name and restart the program', "autowin", 1) 
	os.mkdir(outdir)
	t       = EMData()
	e_n     = EMData()
	e       = get_image(noisemic)
	i_tem   = EMUtil.get_image_count(templatefile)
	e_n.read_image(templatefile, 0)
	nx      = e_n.get_xsize()
	ny      = e_n.get_ysize()
	rad     = int(nx/2)-1
	mask    = model_circle(rad, nx, ny)
	f       = open(noisedoc, "r")
	rstring = f.readlines() 
	f.close
	xs      = rstring[0]
	tmp     = split(xs)
	x       = tmp[0]
	y       = tmp[1]
	x_noi   = int(x)-p_size/2
	y_noi   = int(y)-p_size/2
	reg     = Region(x_noi,y_noi, p_size, p_size)# Get the reference noise from the input mic and noise coordinates
	e_n     = e.get_clip(reg)
	if CTF : ctf_dicts = ["defocus", "Pixel_size", "voltage", "Cs", "amp_contrast", "B_factor", "sign"]
	for i, v in enumerate(flist):
		(filename, filextension) = os.path.splitext(v)
		out1_file     = "coord_" + filename + ".txt"
		out2_file     = "particle_" + filename + ".hdf"
		print "micrograph being processed:", filename
		f_coord       = os.path.join(outdir, out1_file)
		file_particle = os.path.join(outdir, out2_file) #  parse output file names 			
		peaks         = []
		file_input    = os.path.join(indir, v)
		img1 = EMData()
		img1.read_image(file_input)
		if CTF : ctf_params    = img1.get_attr( "ctf" )
		img1         *= contrast_invert
		nx            = img1.get_xsize()
		ny            = img1.get_ysize()
		N_ptl         = int(nx*ny/p_size/p_size) # number of possible particles
		if(N_ptl >= n_peak_max):	N_ptl = n_peak_max
		sigma_win     = float(float(sigma)/float(p_size)) # filter radius
		nx_d          = int(nx/deci)
		ny_d          = int(ny/deci)
		nx_fft_p      = smallprime(nx_d)
		ny_fft_p      = smallprime(ny_d)
		nx_fft_m      = nx_fft_p*int(deci)
 		ny_fft_m      = ny_fft_p*int(deci)
		img1          = window2d(img1, nx_fft_m, ny_fft_m, "l")
		if(CC_method == 1):
			if(int(deci) == 1):
				img1         = filt_gaussh(fft(img1), sigma_win)
			else:
				feq_deci = 0.5/deci
				img1       = Util.decimate(fft(filt_tanl(filt_gaussh(fft(img1), sigma_win), feq_deci, 0.04)), int(deci), int(deci),1)
				img1       = fft(img1)
		else:
			feq_deci           = 0.5/deci
			img2               = Util.decimate(filt_tanl(img1, feq_deci, 0.04), int(deci), int(deci), 1)
			img2               = fft(img2)
			img1               = Util.decimate(fft(filt_tanl(filt_gaussh(fft(img1,sigma_win)), feq_deci, 0.04)), int(deci),int(deci), 1)
		for j in xrange(i_tem):
			t.read_image(templatefile, j)
			if(int(CC_method) == 2):	cc_map = flcc(t, img2)
			else:
				t_pad  = Util.pad(t, nx_fft_p, ny_fft_p, 1,0,0,0, "circumference")
				cc_map = ccf(img1, t_pad)
				del t_pad
			peaks.insert(0, cc_map.peak_ccf(p_size/2-1.0))
		if(int(CC_method) == 2): del img2
		peak = peaks[0]
		for j in xrange(1,i_tem):#output results
			peak1 = Util.merge_peaks(peak, peaks[j], hf_p)
			del peak
			peak  = peak1
		n_peak = int(len(peak)/3)
		if n_peak <= N_ptl :	N_wi=int(n_peak)
		else:			N_wi=int(N_ptl )
		out = open(f_coord, "w")
		out.write("#Coordinates: %s\n")
		if N_wi == 0 :	ERROR("Number of particles is zero", "autowin", 0)
		if(CC_method == 1):  img1 = fft(img1)
		for k in xrange(N_wi):
			x       = peak[k*3+1] -p_size/2
			y       = peak[k*3+2] -p_size/2
			# print "x==",x, "y===",y, " ccf==",peak[k*3]
			out.write("%d\t%f\t%f\n" % (k+1,x,y))
			reg     = Region(x,y, p_size, p_size)
			wi      = img1.get_clip(reg)
			ra      = ramp(wi)
			outlist = ce_fit(ra, e_n, mask)
			outlist[2].set_attr_dict({'xp':peak[k*3+1], 'yp': peak[k*3+2], 'mic':filename})
			if CTF: outlist[2].set( "ctf", ctf_params)
			outlist[2].write_image(file_particle, k)
		out.close()

def autowin_MPI(indir,outdir, noisedoc, noisemic, templatefile, deci, CC_method, p_size, sigma, hf_p,n_peak_max, contrast_invert=1, prefix_of_micrograph="micrograph"):
	""" 
		automated particle picking program
	   	CC_method=1  Applying high pass filter to detect protein particles
		sigma can be adjusted to perform particle picking. 
		CC_method=2  Using fast local normalization method to detect protein particles
		The output files are detected coordinates and windowed particles
	"""
	from utilities 		import get_image, model_circle,ce_fit,drop_image,info
	from fundamentals 	import smallprime, window2d, ccf, ramp, fft
	from filter 		import filt_gaussh, filt_tanl
	from string 		import split
	from morphology 	import flcc
	from random     	import randint
	import sys
	import os
	from mpi 	    import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	    import mpi_barrier, mpi_bcast
	from mpi 	    import MPI_INT

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	#  chose a random node as a main one...
	main_node      = 0
	if(myid == 0): main_node = randint(0,number_of_proc-1)
	main_node      = mpi_bcast(main_node, 1, MPI_INT, 0, MPI_COMM_WORLD)
	if os.path.exists(indir)  is False: ERROR("micrograph directory does not exsit", "autowin_MPI",1,myid)	
	flist = os.listdir(indir)
	nima          = 0
	mic_name_list = []
	for i, v in enumerate(flist):
		micname                  = os.path.join(indir,v)
		(filename, filextension) = os.path.splitext(v  )
		if(filename[0:len(prefix_of_micrograph)] == prefix_of_micrograph):
			mic_name_list.append(micname)
			nima += 1
	
	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "autowin_MPI ", 1,myid)
	mpi_barrier(MPI_COMM_WORLD)
	
	if myid == int(main_node): # directory cleaning only performed by main node
		print_begin_msg("autowin_MPI")
		os.mkdir(outdir)
	mpi_barrier(MPI_COMM_WORLD)
	
	#if(number_of_proc <= nima ):	nimage_per_node = nima/number_of_proc
	#else: 				nimage_per_node = 1 
	#image_start     = myid * nimage_per_node
	#if(myid == number_of_proc-1):  image_end = nima
	#else:                          image_end = image_start + nimage_per_node
	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
	if image_start == image_end: 
		if image_end == nima:	image_start -= 1
		else: image_end += 1
	
	t       = EMData()
	e_n     = EMData()
	e       = get_image(noisemic)
	i_tem   = EMUtil.get_image_count(templatefile)
	e_n.read_image(templatefile,0)
	nx      = e_n.get_xsize()
	ny      = e_n.get_ysize()
	rad     = int(nx/2)-1
	mask    = model_circle(rad, nx, ny)
	f       = open(noisedoc, "r")
	rstring = f.readlines() 
	f.close
	xs      = rstring[0]
	tmp     = split(xs)
	x       = tmp[0]
	y       = tmp[1]
	x_noi   = int(x)-p_size/2
	y_noi   = int(y)-p_size/2
	reg     = Region(x_noi,y_noi, p_size, p_size)# Get the reference noise from the input mic and noise coordinates
	e_n     = e.get_clip(reg)
	
	for i in xrange(image_start,image_end):
		filename=mic_name_list[i] 
		print '%-15s%-30s'%("micrographs # ",filename)
		(f_nam, filextension) = os.path.splitext(filename)
		out1_file     = "coord_"+f_nam[len(prefix_of_micrograph)+len(indir)+2:]+".txt"
		out2_file     = "particle_"+f_nam[len(prefix_of_micrograph)+len(indir)+2:]+filextension
		f_coord       = os.path.join(outdir, out1_file)
		file_particle = os.path.join(outdir, out2_file) #  parse output file names 			
		peaks         = []
		img1 = EMData()
		img1.read_image(filename)
		img1          *= contrast_invert
		nx            = img1.get_xsize()
		ny            = img1.get_ysize()
		N_ptl         = int(nx*ny/p_size/p_size) # number of possible particles
		if(N_ptl >= n_peak_max):	N_ptl = n_peak_max
		sigma_win     = float(float(sigma)/float(p_size)) # filter radius
		nx_d          = int(nx/deci)
		ny_d          = int(ny/deci)
		nx_fft_p      = smallprime(nx_d)
		ny_fft_p      = smallprime(ny_d)
		nx_fft_m      = nx_fft_p*int(deci)
 		ny_fft_m      = ny_fft_p*int(deci)
		img1          = window2d(img1, nx_fft_m, ny_fft_m, "l")
		if(CC_method == 1):
			if(int(deci) == 1):
				img1         = filt_gaussh(fft(img1), sigma_win)
			else:
				feq_deci = 0.5/deci
				img1       = Util.decimate(fft(filt_tanl(filt_gaussh(fft(img1), sigma_win), feq_deci, 0.04)), int(deci), int(deci),1)
				img1       = fft(img1)
		else:
			feq_deci           = 0.5/deci
			img2               = Util.decimate(filt_tanl(img1, feq_deci, 0.04), int(deci), int(deci), 1)
			img2               = fft(img2)
			img1               = Util.decimate(fft(filt_tanl(filt_gaussh(fft(img1,sigma_win)), feq_deci, 0.04)), int(deci),int(deci), 1)
		for j in xrange(i_tem):
			t.read_image(templatefile, j)
			if(int(CC_method) == 2):	cc_map = flcc(t, img2)
			else:
				t_pad  = Util.pad(t, nx_fft_p, ny_fft_p, 1,0,0,0, "circumference")
				cc_map = ccf(img1, t_pad)
				del t_pad
			peaks.insert(0,cc_map.peak_ccf(p_size/2-1.0))
		if(int(CC_method) == 2): del img2
		peak = peaks[0]
		for j in xrange(1,i_tem):#output results
			peak1 = Util.merge_peaks(peak, peaks[j], hf_p)
			del peak
			peak  = peak1
		n_peak = int(len(peak)/3)
		if n_peak <= N_ptl :	N_wi=int(n_peak)
		else:			N_wi=int(N_ptl )			
		out = open(f_coord, "w")
		out.write("#Coordinates: %s\n")
		if N_wi == 0 :	ERROR("Number of particles is zero","autowin.py",0,myid)
		if(CC_method == 1):img1 = fft(img1)			
		for k in xrange(N_wi):
			x       = peak[k*3+1] -p_size/2
			y       = peak[k*3+2] -p_size/2
			# print "x==",x, "y===",y, " ccf==",peak[k*3]
			out.write("%d\t%f\t%f\n" % (k+1,x,y))
			reg     = Region(x, y, p_size, p_size)
			wi      = img1.get_clip(reg)
			ra      = ramp(wi)
			outlist = ce_fit(ra,e_n,mask)
			outlist[2].set_attr_dict({'xp':peak[k*3+1], 'yp': peak[k*3+2],'mic':filename})
			outlist[2].write_image(file_particle, k)
		out.close()


def ihrsr(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, ynumber, 
          txs, delta, initial_theta, delta_theta, an, maxit, CTF, snr, dp, ndp, dp_step, dphi, ndhpi, dphi_step, psi_max,
	  rmin, rmax, fract, nise, npad, sym, user_func_name, datasym,
	  fourvar, debug = False, MPI = False, chunk = -1.0, WRAP = 1):
	if MPI:
		if (chunk <= 0.0):
			if WRAP == 1:
				ihrsr_MPI(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, ynumber, 
			txs, delta, initial_theta, delta_theta, an, maxit, CTF, snr, dp, ndp, dp_step, dphi, ndhpi, dphi_step, psi_max,
			rmin, rmax, fract, nise, npad, sym, user_func_name, datasym,
			fourvar, debug)
			
			else:
				ihrsr_MPI_no_wrap(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, ynumber, 
			txs, delta, initial_theta, delta_theta, an, maxit, CTF, snr, dp, ndp, dp_step, dphi, ndhpi, dphi_step, psi_max,
			rmin, rmax, fract, nise, npad, sym, user_func_name, datasym,
			fourvar, debug)
		else:
			ihrsr_chunk_MPI(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, ynumber, 
			txs, delta, initial_theta, delta_theta, an, maxit, CTF, snr, dp, ndp, dp_step, dphi, ndhpi, dphi_step, psi_max,
			rmin, rmax, fract, nise, npad, sym, user_func_name, datasym,
			fourvar, debug, chunk)
		
		return
	'''
	from utilities      import model_circle, drop_image
	from utilities      import get_image, get_input_from_string
	from utilities      import get_params_proj, set_params_proj
	from alignment	    import proj_ali_helical, helios, Numrinit, prepare_refrings
	from projection     import prep_vol
	from statistics     import ccc
	from fundamentals   import cyclic_shift, rot_shift3D
	import os
	import types
	from utilities      import print_begin_msg, print_end_msg, print_msg
	print_begin_msg("ihrsr")

	import user_functions
	user_func = user_functions.factory[user_func_name]

	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ihrsr", 1)
	os.mkdir(outdir)
        
	xrng        = get_input_from_string(xr)
	#Guozhi Tao changed----since I cannot test this code, just modify the code to make sure no crash--begin
	if  ny == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(xr)
	#since I cannot test this code, just modify the code to make sure no crash--end
	step        = get_input_from_string(txs)
	delta       = get_input_from_string(delta)
	lstp = min( len(xrng), len(yrng), len(step), len(delta) )
	if (an == "-1"):
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)
	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)

	print_msg("Input stack                               : %s\n"%(stack))
	print_msg("Reference volume                          : %s\n"%(ref_vol))	
	print_msg("Output directory                          : %s\n"%(outdir))
	print_msg("Maskfile                                  : %s\n"%(maskfile))
	print_msg("Inner radius                              : %i\n"%(first_ring))

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	if (last_ring == -1):	last_ring = nx//2 - 2

	print_msg("Outer radius                              : %i\n"%(last_ring))
	print_msg("Ring step                                 : %i\n"%(rstep))
	print_msg("X search range                            : %s\n"%(xrng))
	print_msg("Y search range                            : %s\n"%(yrng))
	print_msg("Translational step                        : %s\n"%(step))
	print_msg("Angular step                              : %s\n"%(delta))
	print_msg("Angular search range                      : %s\n"%(an))
	print_msg("max radius for helical search (in Ang)    : %f\n"%(rmax))
	print_msg("fraction of volume used for helical search: %f\n"%(fract))
	print_msg("initial symmetry - angle                  : %f\n"%(dphi))
	print_msg("initial symmetry - axial rise             : %f\n"%(dp))
	print_msg("Maximum number of iterations              : %i\n"%(max_iter))
	print_msg("CTF correction                            : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio                     : %f\n"%(snr))
	print_msg("Symmetry group                            : %s\n"%(sym))
	print_msg("symmetry doc file                         : %s\n"%(datasym))
	print_msg("npad                                      : %i\n"%(npad))
	print_msg("User function                             : %s\n"%(user_func_name))

	if (maskfile) :
		if (type(maskfile) is types.StringType): mask3D = get_image(maskfile)
		else                                  : mask3D = maskfile
	else          :   mask3D = model_circle(last_ring, nx, nx, nx)
	mask2D = model_circle(last_ring, nx, nx) - model_circle(first_ring, nx, nx)
	numr   = Numrinit(first_ring, last_ring, rstep, "F")

	if CTF:
		from reconstruction import recons3d_4nn_ctf
		from filter import filt_ctf

	else: from reconstruction import recons3d_4nn


	#drop_image(vol, os.path.join(outdir,"ref_vol00.hdf"))
	#sym = "c1"
	symref = "s"+sym
	ref_a= "P"


	active = EMUtil.get_all_attributes(stack, 'active')
	list_of_particles = []
	for im in xrange(len(active)):
		if(active[im]):  list_of_particles.append(im)
	del active
	data = EMData.read_images(stack, list_of_particles)
        for im in xrange(len(data)):
                data[im].set_attr('ID', list_of_particles[im])
	nima = len(data)
	if(data[0].get_attr_default('ctf_applied', 2) > 0):  ctf_applied = True
	else:   ctf_applied = False
	if CTF:
		pixel_size = data[0].get_attr('ctf').apix
	else:
		pixel_size = data[0].get_attr('pixel_size')
	print_msg("Pixel size in Angstroms                   : %f\n\n"%(pixel_size))

	finfo = None#open("desperado", 'w')
	# do the projection matching
	drop_image(vol, os.path.join(outdir, "aligned0000.hdf"))
	for N_step in xrange(lstp):
		for Iter in xrange(max_iter):
			print_msg("ITERATION #%3d\n"%(N_step*max_iter + Iter+1))

			if CTF:
				if ctf_applied == False:
					previous_defocus = -1.0
				else:
					volft,kb = prep_vol( vol )
					refrings = prepare_refrings( volft, kb, nx, delta[N_step], ref_a, symref, numr, MPI=False, phiEqpsi = "Zero")
					del volft, kb
			else:
				volft,kb = prep_vol( vol )
				refrings = prepare_refrings( volft, kb, nx, delta[N_step], ref_a, symref, numr, MPI=False, phiEqpsi = "Zero")
				del volft, kb
			sx = 0.0
			sy = 0.0
			phihi = [0.0]*nima
			for im in xrange( nima ):
				if CTF and ctf_applied == False:
					ctf = data[im].get_attr( "ctf" )
					if ctf.defocus != previous_defocus:
						previous_defocus = ctf.defocus
						ctfvol = filt_ctf(vol, ctf)
						volft,kb = prep_vol( ctfvol )
						refrings = prepare_refrings( volft, kb, nx, delta[N_step], ref_a, symref, numr, MPI=False, phiEqpsi = "Zero")
						del volft, kb

				#if an[N_step] == -1:
				#	peak, pixel_error = proj_ali_helical(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step])
				#else:
				peak, phihi[im], sxi, syi, pixel_error = proj_ali_helical(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step], psi_max)

				sx += sxi
				sy += syi
				#print_msg("Image %i,  psi %9.2f,    s2x  %9.2f, s2y  %9.2f,  peak  %10.3e \n"%(im, paramali[2], paramali[3], paramali[4], peak))

			# histogram of phi's
			from statistics import hist_list
			lhist = 30
			region, histo = hist_list(phihi, lhist)
			msg = "      Histogram of phi angles\n      phi         number of particles\n"
			print_msg(msg)
			for lhx in xrange(lhist):
				msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
				print_msg(msg)
			del phihi, region, histo, msg
			#center projections
			sx /= nima
			sy /= nima
			for im in xrange( nima ):
				paramali = get_params_proj(data[im])
				set_params_proj(data[im], [paramali[0], paramali[1], paramali[2], paramali[3] - sx, paramali[4] - sy])

			#  3D stuff
			#  calculate new and improved 3D
			if(CTF): vol = recons3d_4nn_ctf(data, range(nima), snr, npad = npad)
			else:	 vol = recons3d_4nn(data, range(nima), npad = npad)

			ref_data = [vol]
			vol = user_func(ref_data)

			# store the reference volume
			drop_image(vol, os.path.join(outdir, "unsymmetrized%04d.hdf"%(N_step*max_iter+Iter+1)))
			if(N_step*max_iter+Iter+1 > 2):
				vol, dp, dphi = helios(vol, pixel_size, dp, dphi, fract, rmax, rmin)
			else:
				#  in the first two steps the symmetry is imposed
				vol = vol.helicise(pixel_size,dp, dphi, fract, rmax, rmin)
			print_msg("New delta z and delta phi      : %s,    %s\n\n"%(dp,dphi))
			fofo = open(os.path.join(outdir,datasym),'a')
			fofo.write('  %12.4f   %12.4f\n'%(dp,dphi))
			fofo.close()
			drop_image(vol, os.path.join(outdir, "aligned%04d.hdf"%(N_step*max_iter+Iter+1)) )
			#  here we  write header info
			from utilities import write_headers
			write_headers( stack, data, list_of_particles)
	print_end_msg("ihrsr")
	'''
def ihrsr_MPI(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, ynumber, 
	txs, delta, initial_theta, delta_theta, an, maxit, CTF, snr, dp, ndp, dp_step, dphi, ndphi, dphi_step, psi_max,
	rmin, rmax, fract, nise, npad, sym, user_func_name, datasym,
	fourvar, debug):

	from alignment      import Numrinit, prepare_refrings, proj_ali_helical, proj_ali_helical_90, proj_ali_helical_local, proj_ali_helical_90_local, helios,helios7
	from utilities      import model_circle, get_image, drop_image, get_input_from_string, pad, model_blank
	from utilities      import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all
	from utilities      import send_attr_dict
	from utilities      import get_params_proj, set_params_proj, file_type
	from fundamentals   import rot_avg_image
	from pixel_error    import max_3D_pixel_error
	import os
	import types
	from utilities      import print_begin_msg, print_end_msg, print_msg
	from mpi            import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi            import mpi_recv,  mpi_send, MPI_TAG_UB
	from mpi            import mpi_reduce, MPI_INT, MPI_SUM
	from filter         import filt_ctf
	from projection     import prep_vol, prgs
	from statistics     import hist_list, varf3d_MPI
	from applications   import MPI_start_end
	from EMAN2 import Vec2f
	from string    import lower,split
	from math import cos, pi

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	if myid == 0:
		if os.path.exists(outdir):  nx = 1
		else:  nx = 0
	else:  nx = 0
	ny = bcast_number_to_all(nx, source_node = main_node)
	
	if ny == 1:  ERROR('Output directory exists, please change the name and restart the program', "ihrsr_MPI", 1,myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		os.mkdir(outdir)
	mpi_barrier(MPI_COMM_WORLD)
	
	'''ndp    = 12
	sndp   = 0.1
	ndphi  = 12
	sndphi = 0.1'''
	nlprms = (2*ndp+1)*(2*ndphi+1)
	if nlprms< number_of_proc:
		ERROR('number of cpus is larger than the number of helical search, please reduce it or at this moment modify ndp,dphi in the program', "ihrsr_MPI", 1,myid)
	


	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None


	#sym = "c1"
	symref = "s"+sym

	ref_a= "P"
	symmetryLower = sym.lower()
	symmetry_string = split(symmetryLower)[0]

	xrng        = get_input_from_string(xr)
	ynumber	    = get_input_from_string(ynumber)
	for i in xrange(len(ynumber)):
		if(ynumber[i]%2==1):
			ynumber[i]=ynumber[i]+1
	yrng =[]

	for i in xrange(len(xrng)):
		yrng.append(dp/2)

	stepx        = get_input_from_string(txs)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(stepx), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	ny      = vol.get_ysize()
	nz      = vol.get_zsize()
	nmax = max(nx, ny, nz)
	if ( nx == ny ):
		if( nx == nz):
			xysize = -1
			zsize = -1
		elif( nx < nz):
			xysize = nx
			zsize = -1
		else:
			zsize = nz
			xysize = -1
	
	else:
		ERROR('the x and y size have to be same, please change the reference volume and restart the program', "ihrsr_MPI", 1,myid)
	

	if last_ring < 0:	last_ring = int(nx/2) - 2

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                               : %s\n"%(stack))
		print_msg("Reference volume                          : %s\n"%(ref_vol))	
		print_msg("Output directory                          : %s\n"%(outdir))
		print_msg("Maskfile                                  : %s\n"%(maskfile))
		print_msg("Inner radius                              : %i\n"%(first_ring))
		print_msg("Outer radius                              : %i\n"%(last_ring))
		print_msg("Ring step                                 : %i\n"%(rstep))
		print_msg("X search range                            : %s\n"%(xrng))
		print_msg("Y number                                  : %s\n"%(ynumber))
		print_msg("Translational stepx                       : %s\n"%(stepx))
		print_msg("Angular step                              : %s\n"%(delta))
		print_msg("Angular search range                      : %s\n"%(an))
		print_msg("Initial Theta                             : %s\n"%(initial_theta))
		print_msg("min radius for helical search (in pix)    : %5.4f\n"%(rmin))
		print_msg("max radius for helical search (in pix)    : %5.4f\n"%(rmax))
		print_msg("fraction of volume used for helical search: %5.4f\n"%(fract))
		print_msg("initial symmetry - angle                  : %5.4f\n"%(dphi))
		print_msg("initial symmetry - axial rise             : %5.4f\n"%(dp))
		print_msg("Maximum iteration                         : %i\n"%(max_iter))
		print_msg("Data with CTF                             : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio                     : %5.4f\n"%(snr))
		print_msg("symmetry output doc file                  : %s\n"%(datasym))
		print_msg("number of times initial symmetry is imposed: %i\n"%(nise))
		print_msg("npad                                      : %i\n"%(npad))
		print_msg("User function                             : %s\n"%(user_func_name))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	#else: mask3D = model_circle(last_ring, nx, nx, nx)

	numr	= Numrinit(first_ring, last_ring, rstep, "F")

	if CTF:
		from reconstruction import recons3d_4nn_ctf_MPI
		from filter         import filt_ctf
	else:	 from reconstruction import recons3d_4nn_MPI

	if myid == main_node:
       		if(file_type(stack) == "bdb"):
			from EMAN2db import db_open_dict
			dummy = db_open_dict(stack, True)
		active = EMUtil.get_all_attributes(stack, 'active')
		list_of_particles = []
		for im in xrange(len(active)):
			if active[im]:  list_of_particles.append(im)
		del active
		nima = len(list_of_particles)
	else:
		nima = 0
	total_nima = bcast_number_to_all(nima, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()
	
	mask2D = pad( model_blank( int(nmax-20),nmax,1,bckg=1.0), nmax, nmax, 1,0.0)

	data = EMData.read_images(stack, list_of_particles)
	if fourvar:  original_data = []
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		sttt = Util.infomask(data[im], mask2D, False)
		data[im] = data[im] - sttt[0]
		if fourvar: original_data.append(data[im].copy())
		if CTF:
			st = data[im].get_attr_default("ctf_applied", 0)
			if(st == 0):
				ctf_params = data[im].get_attr("ctf")
				data[im] = filt_ctf(data[im], ctf_params)
				data[im].set_attr('ctf_applied', 1)
		#Util.mul_img(data[im], mask2D)  #?????
	del mask2D

	if debug:
		finfo.write( '%d loaded  \n' % nima )
		finfo.flush()
	if CTF:
		pixel_size = data[0].get_attr('ctf').apix
	else:
		pixel_size = data[0].get_attr('pixel_size')
	for i in xrange(len(xrng)):
		yrng[i]=dp/(2*pixel_size)
	from math import sin, pi
	ou_max = ( nmax/2.0)*sin( initial_theta*pi/180) - dp/2.0/pixel_size -1.0
	if ( ou > ou_max):
		ERROR('ou should be less than----( nmax/2.0)*sin( initial_theta*pi/180) - dp/2.0/pixel_size -2.0 ', "ihrsr_MPI", 1,myid)
	del ou_max
	if myid == main_node:
		print_msg("Pixel size in Angstroms                   : %5.4f\n\n"%(pixel_size))
		print_msg("Y search range (pix) initialized as       : %s\n\n"%(yrng))


	from time import time	

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )

	#jeanmod
	total_iter = 0
	# do the projection matching
	
	for N_step in xrange(lstp):
		terminate = 0
		Iter = 0
 		while(Iter < max_iter and terminate == 0):
			yrng[N_step]=float(dp)/(2*pixel_size) #will change it later according to dp
			if(ynumber[N_step]==0):
				stepy = 0.0
			else:
				stepy = (2*yrng[N_step]/ynumber[N_step])

			pixer  = [0.0]*nima
			modphi = [0.0]*nima
			Iter += 1
			total_iter += 1
			if myid == main_node:
				start_time = time()
				print_msg("\nITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.4f, xrange = %5.4f,stepx = %5.4f, yrange = %5.4f,  stepy = %5.4f, ynumber = %3d\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],stepx[N_step],yrng[N_step],stepy,ynumber[N_step]))
			if( xysize == -1 and zsize==-1 ):
				volft,kb = prep_vol( vol )
				refrings = prepare_refrings( volft, kb, nmax, delta[N_step], ref_a, symref, numr, MPI = True, phiEqpsi = "Zero", initial_theta =initial_theta, delta_theta = delta_theta)
				del volft,kb
			else:
				volft, kbx, kby, kbz = prep_vol( vol )
				refrings = prepare_refrings( volft, kbz, nmax, delta[N_step], ref_a, symref, numr, MPI = True, phiEqpsi = "Zero", kbx = kbx, kby = kby, initial_theta =initial_theta, delta_theta = delta_theta)
				del volft, kbx, kby, kbz
			if myid== main_node:
				print_msg( "Time to prepare rings: %d\n" % (time()-start_time) )
				start_time = time()
			#split refrings to two list: refrings1 (theta =90), and refrings2( theat not 90)
			refrings1= []
			refrings2= []
			sn = int(symmetry_string[1:])
			for i in xrange( len(refrings) ):
				if( sn%2 ==0 and abs( refrings[i].get_attr('n3') ) <1.0e-6 and (symmetry_string[0] == "c" or symmetry_string[0] =="d" ) ):
					refrings1.append( refrings[i])
		
				else:
					refrings2.append( refrings[i])
					'''if myid == main_node:
						print_msg("\nphi = %5.2f, theta = %5.2f, psi=%5.2f\n"%( refrings[i].get_attr('phi'), refrings[i].get_attr('theta'), refrings[i].get_attr('psi') ) )
			if myid == main_node:
				print_msg("\nlen(ref1) = %4d, len(ref2) = %4d\n"%(len(refrings1), len(refrings2)) )'''		
			del refrings
			from numpy import float32
			dpp = float32(float(dp)/pixel_size)
			dpp = float( dpp )
			dpp_half = dpp/2.0
			
			for im in xrange( nima ):
				peak1 = None
				peak2 = None
				if ( len(refrings1) > 0):
					if  an[N_step] == -1:
						peak1, phihi1, theta1, psi1, sxi1, syi1, t11 = proj_ali_helical_90(data[im],refrings1,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step],psi_max,finfo,)
					else:
						peak1, phihi1, theta1, psi1, sxi1, syi1, t11 = proj_ali_helical_90_local(data[im],refrings1,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step], an[N_step], psi_max,  finfo,)
				if( len(refrings2) > 0):
					if  an[N_step] == -1:
						peak2, phihi2, theta2, psi2, sxi2, syi2, t12 = proj_ali_helical(data[im],refrings2,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step],psi_max,finfo,)
					else:
						peak2, phihi2, theta2, psi2, sxi2, syi2, t12 = proj_ali_helical_local(data[im],refrings2,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step], an[N_step], psi_max,finfo,)
				if peak1 is None: 
					peak = peak2
					phihi = phihi2
					theta = theta2
					psi = psi2
					sxi = sxi2
					syi = syi2
					t1 = t12
				elif peak2 is None:
					peak = peak1
					phihi = phihi1
					theta = theta1
					psi = psi1
					sxi = sxi1
					syi = syi1
					t1 = t11
				else:
					if(peak1 >= peak2):
						peak = peak1
						phihi = phihi1
						theta = theta1
						psi = psi1
						sxi = sxi1
						syi = syi1
						t1 = t11
					else:
						peak = peak2
						phihi = phihi2
						theta = theta2
						psi = psi2
						sxi = sxi2
						syi = syi2
						t1 = t12
				#peak, phihi, theta, psi, sxi, syi, t1 = proj_ali_helical(data[im],refrings,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step],psi_max,finfo,)
				if(peak > -1.0e22):
					#Guozhi Tao: wrap y-shifts back into box within rise of one helical unit by changing phi
					
					tp = Transform({"type":"spider","phi":phihi,"theta":theta,"psi":psi})
					tp.set_trans( Vec2f( -sxi, -syi ) )
					dtp = tp.get_params("spider")
					dtp_ty = float( dtp["ty"] )
					del dtp
					if( abs(dtp_ty) >dpp_half):
						dtp_ty_temp = dtp_ty
						if( abs(psi-90) < 90  ):
							sign_psi = 1
						else:
							sign_psi = -1
						if( dtp_ty > 0):
							period_step = -1*sign_psi
						else:
							period_step = 1*sign_psi
						nperiod = 0
						while( abs( dtp_ty_temp ) > dpp_half ):
							nperiod += period_step
							th = Transform({"type":"spider","phi": -nperiod*dphi, "tz":nperiod*dpp})
							tfinal = tp*th
							df = tfinal.get_params("spider")
							dtp_ty_temp = float( df["ty"] )
							
						phihi = float(df["phi"])
						sxi   = float(-df["tx"])
						syi   = float(-df["ty"])
		
					# unique range identified by [k0,k1], [k2,k3]
					tp = Transform({"type":"spider","phi":phihi,"theta":theta,"psi":psi})
					tp.set_trans( Vec2f( -sxi, -syi ) )
					k0 = 0.0
					k2 = k0+180

					if( abs( tp.at(2,2) )<1.0e-6 ):
						if (symmetry_string[0] =="c"):
							if sn%2 == 0:
								k1=360.0/sn
							else:
								k1=360.0/2/sn
						elif (symmetry_string[0] =="d"):
							if sn%2 == 0:
								k1=360.0/2/sn
							else:
								k1=360.0/4/sn
					else:
						k1=360.0/sn
					k3 = k1 +180
					from utilities import get_sym
					T = get_sym(symmetry_string[0:])
					
					for i in xrange( len(T) ):
						ttt = tp*Transform({"type":"spider","phi":T[i][0],"theta":T[i][1],"psi":T[i][2]})
						d1 = ttt.get_params("spider")
						
						if ( abs( tp.at(2,2) )<1.0e-6 ):
							if( sn%2==1 ): # theta=90 and n odd, only one of the two region match

								if( ( d1['phi'] <= float(k1) and d1['phi'] >= float(k0) ) or ( d1['phi'] < float(k3) and d1['phi'] >= float(k2) )):
									
									sxnew = - d1["tx"]
									synew = - d1["ty"]
									phinew = d1['phi']
									thetanew = d1["theta"]
									psinew = d1["psi"]
							else: #for theta=90 and n even, there is no mirror version during aligment, so only consider region [k0,k1]

								if( d1['phi'] <= float(k1) and d1['phi'] >= float(k0) ) :
									
									sxnew = - d1["tx"]
									synew = - d1["ty"]
									phinew = d1['phi']
									thetanew = d1["theta"]
									psinew = d1["psi"]
									
						else: #theta !=90, # if theta >90, put the projection into [k2,k3]. Otherwise put it into the region [k0,k1]
							if( sn==1):
								sxnew = sxi
								synew = syi
								phinew = phihi
								thetanew = theta
								psinew = psi
							else:

								if (tp.at(2,2) >0.0): #theta <90
									
									if( d1['phi'] <= float(k1) and d1['phi'] >= float(k0) ):
										if( cos( pi*float( d1['theta'] )/180.0 )>0.0 ):
											
											sxnew = - d1["tx"]
											synew = - d1["ty"]
											phinew = d1['phi']
											thetanew = d1["theta"]
											psinew = d1["psi"]
											
								else:
									if(  d1['phi'] <= float(k3) and d1['phi'] >= float(k2) ):
										if( cos( pi*float( d1['theta'] )/180.0 )<0.0 ):
											
											sxnew = - d1["tx"]
											synew = - d1["ty"]
											phinew = d1['phi']
											thetanew = d1["theta"]
											psinew = d1["psi"]
											
						del ttt,d1


					t2 = Transform({"type":"spider","phi":phinew,"theta":thetanew,"psi":psinew})
					t2.set_trans(Vec2f(-sxnew, -synew))
					data[im].set_attr("xform.projection", t2)
					pixer[im]  = max_3D_pixel_error(t1, t2, numr[-3])
					modphi[im] = phinew
				else:
					# peak not found, parameters not modified
					pixer[im]  = 0.0
					phihi, theta, psi, sxi, syi = get_params_proj(data[im])
					modphi[im] = phihi

			del refrings1, refrings2
			if myid == main_node:
				print_msg("Time of alignment = %d\n"%(time()-start_time))
				start_time = time()

			#output pixel errors
			from mpi import mpi_gatherv
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			del pixer		
			mpi_barrier(MPI_COMM_WORLD)
			terminate = 0
			if(myid == main_node):
				recvbuf = map(float, recvbuf)
				from utilities import write_text_file
				write_text_file([range(len(recvbuf)), recvbuf], os.path.join(outdir, "pixer_%04d_%04d.txt"%(N_step+1,Iter)) )
				from statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				if(region[0] < 0.0):  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.2f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				# Terminate if 95% within 1 pixel error
				im = 0
				for lhx in xrange(lhist):
					if(region[lhx] > 1.0): break
					im += histo[lhx]
				#if(im/float(total_nima) > 0.95):  terminate = 1
				del region, histo
			#output distribution of phi
			#jeanmod
			recvbuf = mpi_gatherv(modphi, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			#end jeanmod
			mpi_barrier(MPI_COMM_WORLD)
			del modphi
			if(myid == main_node):
				recvbuf = map(float, recvbuf)
				phi_value_0 = []
				phi_value_180 = []
				for i in xrange ( len ( recvbuf ) ):
					if ( recvbuf[i] < 180.0):
						phi_value_0.append( recvbuf[i] )
					else:
						phi_value_180.append( recvbuf[i] ) 
 				lhist = int( round(max(phi_value_0)/delta[N_step]) )
								# if delta is big, number of bins (lhist) will be small, leave it as it is
				# if delta is small, number of bins (lhist) will be big, adjust lhist = lhist/n such as the total 
				# number of bins close to 30, thus most likely we can see each bin contains particles.
				from math import ceil
				if ( len( phi_value_180) > 0):
					if lhist > 15:
						lhist = int(   lhist/ceil((lhist/15.0))  ) 
				else:
					if lhist > 30:
						lhist = int(   lhist/ceil((lhist/30.0))  )  
				region, histo = hist_list(phi_value_0, lhist)
				msg = "\n      Distribution of phi\n      phi         number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.2f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				del region, histo, phi_value_0
				if ( len( phi_value_180) > 0):
					region, histo = hist_list(phi_value_180, lhist)
					for lhx in xrange(lhist):
						msg = " %10.2f     %7d\n"%(region[lhx], histo[lhx])
						print_msg(msg)
					del region, histo, phi_value_180			
			del recvbuf
			terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
			terminate = int(terminate[0])
			if myid == main_node:
				print_msg("Time to compute pixer = %d\n"%(time()-start_time))
				start_time = time()
			# write out headers, under MPI writing has to be done sequentially
			mpi_barrier(MPI_COMM_WORLD)
			m = 5
			from mpi import mpi_recv, mpi_send, MPI_TAG_UB, MPI_COMM_WORLD, MPI_FLOAT
			if myid == main_node:
				
				fexp = open(os.path.join(outdir, "parameters_%04d_%04d.txt"%(N_step+1,Iter)),"w")
				for n in xrange(number_of_proc):
					if n!=main_node:
						t = mpi_recv(recvcount[n]*m,MPI_FLOAT, n, MPI_TAG_UB, MPI_COMM_WORLD)
						for i in xrange(recvcount[n]):
							for j in xrange(m):
								fexp.write(" %15.5f  "%t[j+i*m])
							fexp.write("\n")
					else:
						t = [0.0]*m
						for i in xrange(recvcount[myid]):
							t = get_params_proj(data[i])
							for j in xrange(m):
								fexp.write(" %15.5f  "%t[j])
							fexp.write("\n")
				fexp.close()
				del t
	        	else:
				nvalue = [0.0]*m*recvcount[myid]
				t = [0.0]*m
				for i in xrange(recvcount[myid]):
					t = get_params_proj(data[i])
					for j in xrange(m):
						nvalue[j + i*m] = t[j]
				mpi_send(nvalue, recvcount[myid]*m, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)
				del nvalue
			if myid == main_node:
				print_msg("Time to write parameters = %d\n"%(time()-start_time))
				start_time = time()

			if CTF: vol = recons3d_4nn_ctf_MPI(myid, data, symmetry=sym, snr = snr, npad = npad, xysize = xysize, zsize = zsize)
			else:    vol = recons3d_4nn_MPI(myid, data, symmetry=sym, npad = npad, xysize = xysize, zsize = zsize)

			if myid == main_node:
				print_msg("\n3D reconstruction time = %d\n"%(time()-start_time))
				start_time = time()

			if fourvar:
			#  Compute Fourier variance
				for im in xrange(nima):
					original_data[im].set_attr( 'xform.projection', data[im].get_attr('xform.projection') )
				varf = varf3d_MPI(original_data, ssnr_text_file = os.path.join(outdir, "ssnr%04d"%(total_iter)), mask2D = None, reference_structure = vol, ou = last_ring, rw = 1.0, npad = 1, CTF = CTF, sign = 1, sym =sym, myid = myid)
				if myid == main_node:
					print_msg("Time to calculate 3D Fourier variance= %d\n"%(time()-start_time))
					start_time = time()
					varf = 1.0/varf
					drop_image(varf, os.path.join(outdir, "varf%04d.hdf"%(total_iter)))
			else:  varf = None

			#search for helical symmetry
			if myid == main_node:
				drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(total_iter)))
							
			if(total_iter > nise):
				bcast_EMData_to_all(vol, myid, main_node)
				#from filter import filt_gaussl
				#vol = filt_gaussl(vol, 0.25)

				if myid == main_node:
					lprms = []
					for i in xrange(-ndp,ndp+1,1):
						for j in xrange(-ndphi,ndphi+1,1):
							lprms.append( dp   + i*dp_step)
							lprms.append( dphi + j*dphi_step)
					#print "lprms===",lprms
					recvpara = []
					for im in xrange(number_of_proc):
						helic_ib, helic_ie = MPI_start_end(nlprms, number_of_proc, im)
						recvpara.append(helic_ib )
						recvpara.append(helic_ie )

				para_start, para_end = MPI_start_end(nlprms, number_of_proc, myid)

				list_dps     = [0.0]*((para_end-para_start)*2)
				list_fvalues = [-1.0]*((para_end-para_start)*1)

				if myid == main_node:
					for n in xrange(number_of_proc):
						if n!=main_node: mpi_send(lprms[2*recvpara[2*n]:2*recvpara[2*n+1]], 2*(recvpara[2*n+1]-recvpara[2*n]), MPI_FLOAT, n, MPI_TAG_UB, MPI_COMM_WORLD)
						else:    list_dps = lprms[2*recvpara[2*0]:2*recvpara[2*0+1]]
				else:
					list_dps = mpi_recv((para_end-para_start)*2, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)

				list_dps = map(float, list_dps)

				local_pos = [0.0, 0.0, -1.0e20]
				for i in xrange(para_end-para_start):
					fvalue = helios7(vol, pixel_size, list_dps[i*2], list_dps[i*2+1], fract, rmax, rmin)
					if(fvalue >= local_pos[2]):
						local_pos = [list_dps[i*2], list_dps[i*2+1], fvalue ]
				if myid == main_node:
					list_return = [0.0]*(3*number_of_proc)
					for n in xrange(number_of_proc):
						if n != main_node: list_return[3*n:3*n+3]                 = mpi_recv(3,MPI_FLOAT, n, MPI_TAG_UB, MPI_COMM_WORLD)
 						else:              list_return[3*main_node:3*main_node+3]  = local_pos[:]
				else:
					mpi_send(local_pos, 3, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)

				if myid == main_node:	
					maxvalue = list_return[2]
					for i in xrange(number_of_proc):
						if( list_return[i*3+2] >= maxvalue ):
							maxvalue = list_return[i*3+2]
							dp       = list_return[i*3+0]
							dphi     = list_return[i*3+1]
					dp   = float(dp)
					dphi = float(dphi)
					#print  "  GOT dp dphi",dp,dphi

					vol  = vol.helicise(pixel_size,dp, dphi, fract, rmax, rmin)
					print_msg("New delta z and delta phi      : %s,    %s\n\n"%(dp,dphi))		
			
			else:
				if myid==main_node:
					#  in the first nise steps the symmetry is imposed
					vol = vol.helicise(pixel_size,dp, dphi, fract, rmax, rmin)
					print_msg("Imposed delta z and delta phi      : %s,    %s\n\n"%(dp,dphi))
			if(myid==main_node):
				fofo = open(os.path.join(outdir,datasym),'a')
				fofo.write('  %12.4f   %12.4f\n'%(dp,dphi))
				fofo.close()
				ref_data = [vol]
				if  fourvar:  ref_data.append(varf)
				vol = user_func(ref_data)
				vol = vol.helicise(pixel_size,dp, dphi, fract, rmax, rmin)

				drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(total_iter)))
				print_msg("\nSymmetry search and user function time = %d\n"%(time()-start_time))
				start_time = time()

			bcast_EMData_to_all(vol, myid, main_node)
			dp   = bcast_number_to_all(dp,   source_node = main_node)
			dphi = bcast_number_to_all(dphi, source_node = main_node)
			#
			del varf
	par_str = ["xform.projection"]
	if myid == main_node:
	   	if(file_type(stack) == "bdb"):
	        	from utilities import recv_attr_dict_bdb
	        	recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        else:
	        	from utilities import recv_attr_dict
	        	recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
		print_msg("Time to write header information= %d\n"%(time()-start_time))
		start_time = time()
	else:	       send_attr_dict(main_node, data, par_str, image_start, image_end)
	if myid == main_node: print_end_msg("ihrsr_MPI")
def ihrsr_MPI_no_wrap(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, ynumber, 
	txs, delta, initial_theta, delta_theta, an, maxit, CTF, snr, dp, ndp, dp_step, dphi, ndphi, dphi_step, psi_max,
	rmin, rmax, fract, nise, npad, sym, user_func_name, datasym,
	fourvar, debug):

	from alignment      import Numrinit, prepare_refrings, proj_ali_helical, proj_ali_helical_90, proj_ali_helical_local, proj_ali_helical_90_local, helios,helios7
	from utilities      import model_circle, get_image, drop_image, get_input_from_string, pad, model_blank
	from utilities      import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all
	from utilities      import send_attr_dict
	from utilities      import get_params_proj, set_params_proj, file_type
	from fundamentals   import rot_avg_image
	from pixel_error    import max_3D_pixel_error
	import os
	import types
	from utilities      import print_begin_msg, print_end_msg, print_msg
	from mpi            import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi            import mpi_recv,  mpi_send, MPI_TAG_UB
	from mpi            import mpi_reduce, MPI_INT, MPI_SUM
	from filter         import filt_ctf
	from projection     import prep_vol, prgs
	from statistics     import hist_list, varf3d_MPI
	from applications   import MPI_start_end
	from EMAN2 import Vec2f
	from string    import lower,split
	from math import cos, pi

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	if myid == 0:
		if os.path.exists(outdir):  nx = 1
		else:  nx = 0
	else:  nx = 0
	ny = bcast_number_to_all(nx, source_node = main_node)
	
	if ny == 1:  ERROR('Output directory exists, please change the name and restart the program', "ihrsr_MPI", 1,myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		os.mkdir(outdir)
	mpi_barrier(MPI_COMM_WORLD)
	
	'''ndp    = 12
	sndp   = 0.1
	ndphi  = 12
	sndphi = 0.1'''
	nlprms = (2*ndp+1)*(2*ndphi+1)
	if nlprms< number_of_proc:
		ERROR('number of cpus is larger than the number of helical search, please reduce it or at this moment modify ndp,dphi in the program', "ihrsr_MPI", 1,myid)
	


	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None


	#sym = "c1"
	symref = "s"+sym

	ref_a= "P"
	symmetryLower = sym.lower()
	symmetry_string = split(symmetryLower)[0]

	xrng        = get_input_from_string(xr)
	ynumber	    = get_input_from_string(ynumber)
	for i in xrange(len(ynumber)):
		if(ynumber[i]%2==1):
			ynumber[i]=ynumber[i]+1
	yrng =[]

	for i in xrange(len(xrng)):
		yrng.append(dp/2)

	stepx        = get_input_from_string(txs)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(stepx), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	ny      = vol.get_ysize()
	nz      = vol.get_zsize()
	nmax = max(nx, ny, nz)
	if ( nx == ny ):
		if( nx == nz):
			xysize = -1
			zsize = -1
		elif( nx < nz):
			xysize = nx
			zsize = -1
		else:
			zsize = nz
			xysize = -1
	
	else:
		ERROR('the x and y size have to be same, please change the reference volume and restart the program', "ihrsr_MPI", 1,myid)
	

	if last_ring < 0:	last_ring = int(nx/2) - 2

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]
		print_msg("no helical wrapping \n")
		print_msg("Input stack                               : %s\n"%(stack))
		print_msg("Reference volume                          : %s\n"%(ref_vol))	
		print_msg("Output directory                          : %s\n"%(outdir))
		print_msg("Maskfile                                  : %s\n"%(maskfile))
		print_msg("Inner radius                              : %i\n"%(first_ring))
		print_msg("Outer radius                              : %i\n"%(last_ring))
		print_msg("Ring step                                 : %i\n"%(rstep))
		print_msg("X search range                            : %s\n"%(xrng))
		print_msg("Y number                                  : %s\n"%(ynumber))
		print_msg("Translational stepx                       : %s\n"%(stepx))
		print_msg("Angular step                              : %s\n"%(delta))
		print_msg("Angular search range                      : %s\n"%(an))
		print_msg("min radius for helical search (in pix)    : %5.4f\n"%(rmin))
		print_msg("max radius for helical search (in pix)    : %5.4f\n"%(rmax))
		print_msg("fraction of volume used for helical search: %5.4f\n"%(fract))
		print_msg("initial symmetry - angle                  : %5.4f\n"%(dphi))
		print_msg("initial symmetry - axial rise             : %5.4f\n"%(dp))
		print_msg("Maximum iteration                         : %i\n"%(max_iter))
		print_msg("Data with CTF                             : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio                     : %5.4f\n"%(snr))
		print_msg("symmetry output doc file                  : %s\n"%(datasym))
		print_msg("number of times initial symmetry is imposed: %i\n"%(nise))
		print_msg("npad                                      : %i\n"%(npad))
		print_msg("User function                             : %s\n"%(user_func_name))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	#else: mask3D = model_circle(last_ring, nx, nx, nx)

	numr	= Numrinit(first_ring, last_ring, rstep, "F")

	if CTF:
		from reconstruction import recons3d_4nn_ctf_MPI
		from filter         import filt_ctf
	else:	 from reconstruction import recons3d_4nn_MPI

	if myid == main_node:
       		if(file_type(stack) == "bdb"):
			from EMAN2db import db_open_dict
			dummy = db_open_dict(stack, True)
		active = EMUtil.get_all_attributes(stack, 'active')
		list_of_particles = []
		for im in xrange(len(active)):
			if active[im]:  list_of_particles.append(im)
		del active
		nima = len(list_of_particles)
	else:
		nima = 0
	total_nima = bcast_number_to_all(nima, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()
	
	mask2D = pad( model_blank( int(nmax-20),nmax,1,bckg=1.0), nmax, nmax, 1,0.0)

	data = EMData.read_images(stack, list_of_particles)
	if fourvar:  original_data = []
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		sttt = Util.infomask(data[im], mask2D, False)
		data[im] = data[im] - sttt[0]
		if fourvar: original_data.append(data[im].copy())
		if CTF:
			st = data[im].get_attr_default("ctf_applied", 0)
			if(st == 0):
				ctf_params = data[im].get_attr("ctf")
				data[im] = filt_ctf(data[im], ctf_params)
				data[im].set_attr('ctf_applied', 1)
		#Util.mul_img(data[im], mask2D)  #?????
	del mask2D

	if debug:
		finfo.write( '%d loaded  \n' % nima )
		finfo.flush()
	if CTF:
		pixel_size = data[0].get_attr('ctf').apix
	else:
		pixel_size = data[0].get_attr('pixel_size')
	for i in xrange(len(xrng)):
		yrng[i]=dp/(2*pixel_size)

	if myid == main_node:
		print_msg("Pixel size in Angstroms                   : %5.4f\n\n"%(pixel_size))
		print_msg("Y search range (pix) initialized as       : %s\n\n"%(yrng))


	from time import time	

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )

	#jeanmod
	total_iter = 0
	# do the projection matching
	
	for N_step in xrange(lstp):
		terminate = 0
		Iter = 0
 		while(Iter < max_iter and terminate == 0):
			yrng[N_step]=float(dp)/(2*pixel_size) #will change it later according to dp
			if(ynumber[N_step]==0):
				stepy = 0.0
			else:
				stepy = (2*yrng[N_step]/ynumber[N_step])

			pixer  = [0.0]*nima
			modphi = [0.0]*nima
			Iter += 1
			total_iter += 1
			if myid == main_node:
				start_time = time()
				print_msg("\nITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.4f, xrange = %5.4f,stepx = %5.4f, yrange = %5.4f,  stepy = %5.4f, ynumber = %3d\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],stepx[N_step],yrng[N_step],stepy,ynumber[N_step]))
			if( xysize == -1 and zsize==-1 ):
				volft,kb = prep_vol( vol )
				refrings = prepare_refrings( volft, kb, nmax, delta[N_step], ref_a, symref, numr, MPI = True, phiEqpsi = "Zero", initial_theta =initial_theta, delta_theta = delta_theta)
				del volft,kb
			else:
				volft, kbx, kby, kbz = prep_vol( vol )
				refrings = prepare_refrings( volft, kbz, nmax, delta[N_step], ref_a, symref, numr, MPI = True, phiEqpsi = "Zero", kbx = kbx, kby = kby, initial_theta =initial_theta, delta_theta = delta_theta)
				del volft, kbx, kby, kbz
			if myid== main_node:
				print_msg( "Time to prepare rings: %d\n" % (time()-start_time) )
				start_time = time()
			#split refrings to two list: refrings1 (theta =90), and refrings2( theat not 90)
			refrings1= []
			refrings2= []
			sn = int(symmetry_string[1:])
			for i in xrange( len(refrings) ):
				if( sn%2 ==0 and abs( refrings[i].get_attr('n3') ) <1.0e-6 and (symmetry_string[0] == "c" or symmetry_string[0] =="d" ) ):
					refrings1.append( refrings[i])
		
				else:
					refrings2.append( refrings[i])
					'''if myid == main_node:
						print_msg("\nphi = %5.2f, theta = %5.2f, psi=%5.2f\n"%( refrings[i].get_attr('phi'), refrings[i].get_attr('theta'), refrings[i].get_attr('psi') ) )
			if myid == main_node:
				print_msg("\nlen(ref1) = %4d, len(ref2) = %4d\n"%(len(refrings1), len(refrings2)) )'''		
			del refrings
			for im in xrange( nima ):
				peak1 = None
				peak2 = None
				if ( len(refrings1) > 0):
					if  an[N_step] == -1:
						peak1, phihi1, theta1, psi1, sxi1, syi1, t11 = proj_ali_helical_90(data[im],refrings1,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step],psi_max,finfo,)
					else:
						peak1, phihi1, theta1, psi1, sxi1, syi1, t11 = proj_ali_helical_90_local(data[im],refrings1,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step], an[N_step], psi_max,  finfo,)
				if( len(refrings2) > 0):
					if  an[N_step] == -1:
						peak2, phihi2, theta2, psi2, sxi2, syi2, t12 = proj_ali_helical(data[im],refrings2,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step],psi_max,finfo,)
					else:
						peak2, phihi2, theta2, psi2, sxi2, syi2, t12 = proj_ali_helical_local(data[im],refrings2,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step], an[N_step], psi_max,finfo,)
				if peak1 is None: 
					peak = peak2
					phihi = phihi2
					theta = theta2
					psi = psi2
					sxi = sxi2
					syi = syi2
					t1 = t12
				elif peak2 is None:
					peak = peak1
					phihi = phihi1
					theta = theta1
					psi = psi1
					sxi = sxi1
					syi = syi1
					t1 = t11
				else:
					if(peak1 >= peak2):
						peak = peak1
						phihi = phihi1
						theta = theta1
						psi = psi1
						sxi = sxi1
						syi = syi1
						t1 = t11
					else:
						peak = peak2
						phihi = phihi2
						theta = theta2
						psi = psi2
						sxi = sxi2
						syi = syi2
						t1 = t12
				#peak, phihi, theta, psi, sxi, syi, t1 = proj_ali_helical(data[im],refrings,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step],psi_max,finfo,)
				if(peak > -1.0e22):
						
					# unique range identified by [k0,k1], [k2,k3]
					tp = Transform({"type":"spider","phi":phihi,"theta":theta,"psi":psi})
					tp.set_trans( Vec2f( -sxi, -syi ) )
					k0 = 0.0
					k2 = k0+180

					if( abs( tp.at(2,2) )<1.0e-6 ):
						if (symmetry_string[0] =="c"):
							if sn%2 == 0:
								k1=360.0/sn
							else:
								k1=360.0/2/sn
						elif (symmetry_string[0] =="d"):
							if sn%2 == 0:
								k1=360.0/2/sn
							else:
								k1=360.0/4/sn
					else:
						k1=360.0/sn
					k3 = k1 +180
					from utilities import get_sym
					T = get_sym(symmetry_string[0:])
					
					for i in xrange( len(T) ):
						ttt = tp*Transform({"type":"spider","phi":T[i][0],"theta":T[i][1],"psi":T[i][2]})
						d1 = ttt.get_params("spider")
						
						if ( abs( tp.at(2,2) )<1.0e-6 ):
							if( sn%2==1 ): # theta=90 and n odd, only one of the two region match

								if( ( d1['phi'] <= float(k1) and d1['phi'] >= float(k0) ) or ( d1['phi'] < float(k3) and d1['phi'] >= float(k2) )):
									
									sxnew = - d1["tx"]
									synew = - d1["ty"]
									phinew = d1['phi']
									thetanew = d1["theta"]
									psinew = d1["psi"]
							else: #for theta=90 and n even, there is no mirror version during aligment, so only consider region [k0,k1]

								if( d1['phi'] <= float(k1) and d1['phi'] >= float(k0) ) :
									
									sxnew = - d1["tx"]
									synew = - d1["ty"]
									phinew = d1['phi']
									thetanew = d1["theta"]
									psinew = d1["psi"]
									
						else: #theta !=90, # if theta >90, put the projection into [k2,k3]. Otherwise put it into the region [k0,k1]
							if( sn==1):
								sxnew = sxi
								synew = syi
								phinew = phihi
								thetanew = theta
								psinew = psi
							else:

								if (tp.at(2,2) >0.0): #theta <90
									
									if( d1['phi'] <= float(k1) and d1['phi'] >= float(k0) ):
										if( cos( pi*float( d1['theta'] )/180.0 )>0.0 ):
											
											sxnew = - d1["tx"]
											synew = - d1["ty"]
											phinew = d1['phi']
											thetanew = d1["theta"]
											psinew = d1["psi"]
											
								else:
									if(  d1['phi'] <= float(k3) and d1['phi'] >= float(k2) ):
										if( cos( pi*float( d1['theta'] )/180.0 )<0.0 ):
											
											sxnew = - d1["tx"]
											synew = - d1["ty"]
											phinew = d1['phi']
											thetanew = d1["theta"]
											psinew = d1["psi"]
											
						del ttt,d1


					t2 = Transform({"type":"spider","phi":phinew,"theta":thetanew,"psi":psinew})
					t2.set_trans(Vec2f(-sxnew, -synew))
					data[im].set_attr("xform.projection", t2)
					pixer[im]  = max_3D_pixel_error(t1, t2, numr[-3])
					modphi[im] = phinew
				else:
					# peak not found, parameters not modified
					pixer[im]  = 0.0
					phihi, theta, psi, sxi, syi = get_params_proj(data[im])
					modphi[im] = phihi

			del refrings1, refrings2
			if myid == main_node:
				print_msg("Time of alignment = %d\n"%(time()-start_time))
				start_time = time()

			#output pixel errors
			from mpi import mpi_gatherv
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			del pixer		
			mpi_barrier(MPI_COMM_WORLD)
			terminate = 0
			if(myid == main_node):
				recvbuf = map(float, recvbuf)
				from utilities import write_text_file
				write_text_file([range(len(recvbuf)), recvbuf], os.path.join(outdir, "pixer_%04d_%04d.txt"%(N_step+1,Iter)) )
				from statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				if(region[0] < 0.0):  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.2f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				# Terminate if 95% within 1 pixel error
				im = 0
				for lhx in xrange(lhist):
					if(region[lhx] > 1.0): break
					im += histo[lhx]
				#if(im/float(total_nima) > 0.95):  terminate = 1
				del region, histo
			#output distribution of phi
			#jeanmod
			recvbuf = mpi_gatherv(modphi, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			#end jeanmod
			mpi_barrier(MPI_COMM_WORLD)
			del modphi
			if(myid == main_node):
				recvbuf = map(float, recvbuf)
				phi_value_0 = []
				phi_value_180 = []
				for i in xrange ( len ( recvbuf ) ):
					if ( recvbuf[i] < 180.0):
						phi_value_0.append( recvbuf[i] )
					else:
						phi_value_180.append( recvbuf[i] ) 
 				lhist = int( round(max(phi_value_0)/delta[N_step]) )
								# if delta is big, number of bins (lhist) will be small, leave it as it is
				# if delta is small, number of bins (lhist) will be big, adjust lhist = lhist/n such as the total 
				# number of bins close to 30, thus most likely we can see each bin contains particles.
				from math import ceil
				if ( len( phi_value_180) > 0):
					if lhist > 15:
						lhist = int(   lhist/ceil((lhist/15.0))  ) 
				else:
					if lhist > 30:
						lhist = int(   lhist/ceil((lhist/30.0))  )  
				region, histo = hist_list(phi_value_0, lhist)
				msg = "\n      Distribution of phi\n      phi         number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.2f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				del region, histo, phi_value_0
				if ( len( phi_value_180) > 0):
					region, histo = hist_list(phi_value_180, lhist)
					for lhx in xrange(lhist):
						msg = " %10.2f     %7d\n"%(region[lhx], histo[lhx])
						print_msg(msg)
					del region, histo, phi_value_180			
			del recvbuf
			terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
			terminate = int(terminate[0])
			if myid == main_node:
				print_msg("Time to compute pixer = %d\n"%(time()-start_time))
				start_time = time()
			# write out headers, under MPI writing has to be done sequentially
			mpi_barrier(MPI_COMM_WORLD)
			m = 5
			from mpi import mpi_recv, mpi_send, MPI_TAG_UB, MPI_COMM_WORLD, MPI_FLOAT
			if myid == main_node:
				
				fexp = open(os.path.join(outdir, "parameters_%04d_%04d.txt"%(N_step+1,Iter)),"w")
				for n in xrange(number_of_proc):
					if n!=main_node:
						t = mpi_recv(recvcount[n]*m,MPI_FLOAT, n, MPI_TAG_UB, MPI_COMM_WORLD)
						for i in xrange(recvcount[n]):
							for j in xrange(m):
								fexp.write(" %15.5f  "%t[j+i*m])
							fexp.write("\n")
					else:
						t = [0.0]*m
						for i in xrange(recvcount[myid]):
							t = get_params_proj(data[i])
							for j in xrange(m):
								fexp.write(" %15.5f  "%t[j])
							fexp.write("\n")
				fexp.close()
				del t
	        	else:
				nvalue = [0.0]*m*recvcount[myid]
				t = [0.0]*m
				for i in xrange(recvcount[myid]):
					t = get_params_proj(data[i])
					for j in xrange(m):
						nvalue[j + i*m] = t[j]
				mpi_send(nvalue, recvcount[myid]*m, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)
				del nvalue
			if myid == main_node:
				print_msg("Time to write parameters = %d\n"%(time()-start_time))
				start_time = time()

			if CTF: vol = recons3d_4nn_ctf_MPI(myid, data, symmetry=sym, snr = snr, npad = npad, xysize = xysize, zsize = zsize)
			else:    vol = recons3d_4nn_MPI(myid, data, symmetry=sym, npad = npad, xysize = xysize, zsize = zsize)

			if myid == main_node:
				print_msg("\n3D reconstruction time = %d\n"%(time()-start_time))
				start_time = time()

			if fourvar:
			#  Compute Fourier variance
				for im in xrange(nima):
					original_data[im].set_attr( 'xform.projection', data[im].get_attr('xform.projection') )
				varf = varf3d_MPI(original_data, ssnr_text_file = os.path.join(outdir, "ssnr%04d"%(total_iter)), mask2D = None, reference_structure = vol, ou = last_ring, rw = 1.0, npad = 1, CTF = CTF, sign = 1, sym =sym, myid = myid)
				if myid == main_node:
					print_msg("Time to calculate 3D Fourier variance= %d\n"%(time()-start_time))
					start_time = time()
					varf = 1.0/varf
					drop_image(varf, os.path.join(outdir, "varf%04d.hdf"%(total_iter)))
			else:  varf = None

			#search for helical symmetry
			if myid == main_node:
				drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(total_iter)))
							
			if(total_iter > nise):
				bcast_EMData_to_all(vol, myid, main_node)
				#from filter import filt_gaussl
				#vol = filt_gaussl(vol, 0.25)

				if myid == main_node:
					lprms = []
					for i in xrange(-ndp,ndp+1,1):
						for j in xrange(-ndphi,ndphi+1,1):
							lprms.append( dp   + i*dp_step)
							lprms.append( dphi + j*dphi_step)
					#print "lprms===",lprms
					recvpara = []
					for im in xrange(number_of_proc):
						helic_ib, helic_ie = MPI_start_end(nlprms, number_of_proc, im)
						recvpara.append(helic_ib )
						recvpara.append(helic_ie )

				para_start, para_end = MPI_start_end(nlprms, number_of_proc, myid)

				list_dps     = [0.0]*((para_end-para_start)*2)
				list_fvalues = [-1.0]*((para_end-para_start)*1)

				if myid == main_node:
					for n in xrange(number_of_proc):
						if n!=main_node: mpi_send(lprms[2*recvpara[2*n]:2*recvpara[2*n+1]], 2*(recvpara[2*n+1]-recvpara[2*n]), MPI_FLOAT, n, MPI_TAG_UB, MPI_COMM_WORLD)
						else:    list_dps = lprms[2*recvpara[2*0]:2*recvpara[2*0+1]]
				else:
					list_dps = mpi_recv((para_end-para_start)*2, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)

				list_dps = map(float, list_dps)

				local_pos = [0.0, 0.0, -1.0e20]
				for i in xrange(para_end-para_start):
					fvalue = helios7(vol, pixel_size, list_dps[i*2], list_dps[i*2+1], fract, rmax, rmin)
					if(fvalue >= local_pos[2]):
						local_pos = [list_dps[i*2], list_dps[i*2+1], fvalue ]
				if myid == main_node:
					list_return = [0.0]*(3*number_of_proc)
					for n in xrange(number_of_proc):
						if n != main_node: list_return[3*n:3*n+3]                 = mpi_recv(3,MPI_FLOAT, n, MPI_TAG_UB, MPI_COMM_WORLD)
 						else:              list_return[3*main_node:3*main_node+3]  = local_pos[:]
				else:
					mpi_send(local_pos, 3, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)

				if myid == main_node:	
					maxvalue = list_return[2]
					for i in xrange(number_of_proc):
						if( list_return[i*3+2] >= maxvalue ):
							maxvalue = list_return[i*3+2]
							dp       = list_return[i*3+0]
							dphi     = list_return[i*3+1]
					dp   = float(dp)
					dphi = float(dphi)
					#print  "  GOT dp dphi",dp,dphi

					vol  = vol.helicise(pixel_size,dp, dphi, fract, rmax, rmin)
					print_msg("New delta z and delta phi      : %s,    %s\n\n"%(dp,dphi))		
			
			else:
				if myid==main_node:
					#  in the first nise steps the symmetry is imposed
					vol = vol.helicise(pixel_size,dp, dphi, fract, rmax, rmin)
					print_msg("Imposed delta z and delta phi      : %s,    %s\n\n"%(dp,dphi))
			if(myid==main_node):
				fofo = open(os.path.join(outdir,datasym),'a')
				fofo.write('  %12.4f   %12.4f\n'%(dp,dphi))
				fofo.close()
				ref_data = [vol]
				if  fourvar:  ref_data.append(varf)
				vol = user_func(ref_data)
				vol = vol.helicise(pixel_size,dp, dphi, fract, rmax, rmin)

				drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(total_iter)))
				print_msg("\nSymmetry search and user function time = %d\n"%(time()-start_time))
				start_time = time()

			bcast_EMData_to_all(vol, myid, main_node)
			dp   = bcast_number_to_all(dp,   source_node = main_node)
			dphi = bcast_number_to_all(dphi, source_node = main_node)
			#
			del varf
	par_str = ["xform.projection"]
	if myid == main_node:
	   	if(file_type(stack) == "bdb"):
	        	from utilities import recv_attr_dict_bdb
	        	recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        else:
	        	from utilities import recv_attr_dict
	        	recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
		print_msg("Time to write header information= %d\n"%(time()-start_time))
		start_time = time()
	else:	       send_attr_dict(main_node, data, par_str, image_start, image_end)
	if myid == main_node: print_end_msg("ihrsr_MPI")


def ihrsr_chunk_MPI(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, ynumber, 
	txs, delta, initial_theta, delta_theta, an, maxit, CTF, snr, dp, ndp, dp_step, dphi, ndphi, dphi_step, psi_max,
	rmin, rmax, fract, nise, npad, sym, user_func_name, datasym,
	fourvar, debug, chunk):

	from alignment      import Numrinit, prepare_refrings, proj_ali_helical, proj_ali_helical_90, proj_ali_helical_local, proj_ali_helical_90_local, helios,helios7
	from utilities      import model_circle, get_image, drop_image, get_input_from_string, pad, model_blank
	from utilities      import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all
	from utilities      import send_attr_dict
	from utilities      import get_params_proj, set_params_proj, file_type
	from fundamentals   import rot_avg_image
	from pixel_error    import max_3D_pixel_error
	import os
	import types
	from utilities      import print_begin_msg, print_end_msg, print_msg
	from mpi            import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi            import mpi_recv,  mpi_send, MPI_TAG_UB
	from mpi            import mpi_reduce, MPI_INT, MPI_SUM
	from filter         import filt_ctf
	from projection     import prep_vol, prgs
	from statistics     import hist_list, varf3d_MPI
	from applications   import MPI_start_end
	from EMAN2 import Vec2f
	from string    import lower,split
	from math import cos, pi
	from random import shuffle

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	if myid == 0:
		if os.path.exists(outdir):  nx = 1
		else:  nx = 0
	else:  nx = 0
	ny = bcast_number_to_all(nx, source_node = main_node)
	
	if ny == 1:  ERROR('Output directory exists, please change the name and restart the program', "ihrsr_MPI", 1,myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		os.mkdir(outdir)
	mpi_barrier(MPI_COMM_WORLD)
	
	'''ndp    = 12
	sndp   = 0.1
	ndphi  = 12
	sndphi = 0.1'''
	nlprms = (2*ndp+1)*(2*ndphi+1)
	if nlprms< number_of_proc:
		ERROR('number of cpus is larger than the number of helical search, please reduce it or at this moment modify ndp,dphi in the program', "ihrsr_MPI", 1,myid)
	


	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None


	#sym = "c1"
	symref = "s"+sym

	ref_a= "P"
	symmetryLower = sym.lower()
	symmetry_string = split(symmetryLower)[0]

	xrng        = get_input_from_string(xr)
	ynumber	    = get_input_from_string(ynumber)
	for i in xrange(len(ynumber)):
		if(ynumber[i]%2==1):
			ynumber[i]=ynumber[i]+1
	yrng =[]
	
	for i in xrange(len(xrng)):
		yrng.append(dp/2)
	
	stepx        = get_input_from_string(txs)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(stepx), len(delta))
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	ny      = vol.get_ysize()
	nz      = vol.get_zsize()
	nmax = max(nx, ny, nz)
	if ( nx == ny ):
		if( nx == nz):
			xysize = -1
			zsize = -1
		elif( nx < nz):
			xysize = nx
			zsize = -1
		else:
			zsize = nz
			xysize = -1
	
	else:
		ERROR('the x and y size have to be same, please change the reference volume and restart the program', "ihrsr_MPI", 1,myid)
	

	if last_ring < 0:	last_ring = int(nx/2) - 2

	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                               : %s\n"%(stack))
		print_msg("Reference volume                          : %s\n"%(ref_vol))	
		print_msg("Output directory                          : %s\n"%(outdir))
		print_msg("Maskfile                                  : %s\n"%(maskfile))
		print_msg("Inner radius                              : %i\n"%(first_ring))
		print_msg("Outer radius                              : %i\n"%(last_ring))
		print_msg("Ring step                                 : %i\n"%(rstep))
		print_msg("X search range                            : %s\n"%(xrng))
		print_msg("Y number                                  : %s\n"%(ynumber))
		print_msg("Translational stepx                        : %s\n"%(stepx))
		print_msg("Angular step                              : %s\n"%(delta))
		print_msg("Angular search range                      : %s\n"%(an))
		print_msg("min radius for helical search (in pix)    : %5.4f\n"%(rmin))
		print_msg("max radius for helical search (in pix)    : %5.4f\n"%(rmax))
		print_msg("fraction of volume used for helical search: %5.4f\n"%(fract))
		print_msg("initial symmetry - angle                  : %5.4f\n"%(dphi))
		print_msg("initial symmetry - axial rise             : %5.4f\n"%(dp))
		print_msg("Maximum iteration                         : %i\n"%(max_iter))
		print_msg("Data with CTF                             : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio                     : %5.4f\n"%(snr))
		print_msg("symmetry output doc file                  : %s\n"%(datasym))
		print_msg("number of times initial symmetry is imposed: %i\n"%(nise))
		print_msg("npad                                      : %i\n"%(npad))
		print_msg("User function                             : %s\n"%(user_func_name))
		print_msg("chunk                                     : %f\n"%(chunk))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	#else: mask3D = model_circle(last_ring, nx, nx, nx)

	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	if CTF:
		from reconstruction import recons3d_4nn_ctf_MPI
		from filter         import filt_ctf
	else:	 from reconstruction import recons3d_4nn_MPI

	if myid == main_node:
       		if(file_type(stack) == "bdb"):
			from EMAN2db import db_open_dict
			dummy = db_open_dict(stack, True)
		active = EMUtil.get_all_attributes(stack, 'active')
		list_of_particles = []
		for im in xrange(len(active)):
			if active[im]:  list_of_particles.append(im)
		del active
		nima = len(list_of_particles)
	else:
		nima = 0
	total_nima = bcast_number_to_all(nima, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()
	mask2D = pad( model_blank( int(2*rmax),nmax,1,bckg=1.0), nmax, nmax, 1,0.0)
	data = EMData.read_images(stack, list_of_particles)
	if fourvar:  original_data = []
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		sttt = Util.infomask(data[im], mask2D, False)
		data[im] = data[im] - sttt[0]
		if fourvar: original_data.append(data[im].copy())
		if CTF:
			st = data[im].get_attr_default("ctf_applied", 0)
			if(st == 0):
				ctf_params = data[im].get_attr("ctf")
				data[im] = filt_ctf(data[im], ctf_params)
				data[im].set_attr('ctf_applied', 1)
		Util.mul_img(data[im], mask2D)  #?????
	del mask2D

	if debug:
		finfo.write( '%d loaded  \n' % nima )
		finfo.flush()
	if CTF:
		pixel_size = data[0].get_attr('ctf').apix
	else:
		pixel_size = data[0].get_attr('pixel_size')
	for i in xrange(len(xrng)):
		yrng[i]=dp/(2*pixel_size)
		
	if myid == main_node:
		print_msg("Pixel size in Angstroms                   : %5.4f\n\n"%(pixel_size))
		print_msg("Y search range (pix) initialized as       : %s\n\n"%(yrng))
	#  this is needed for determining chunk, each chunk data used for aligment--anima, data used for reconstruction--rnima
	nchunk = int(1.0/chunk+0.5) 
	chunk_list=[0]*2*nchunk
	anima = [0]*nchunk
	rnima =[0]*nchunk
	for i in xrange(nchunk):
		chunk_list[2*i+0] = int(round(float(nima)/nchunk*i))
		chunk_list[2*i+1]   = int(round(float(nima)/nchunk*(i+1)))
		anima[i] = chunk_list[2*i+1]-chunk_list[2*i+0]
		rnima[i] = nima - anima[i]
	

	from time import time	

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )

	#jeanmod
	total_iter = 0
	# do the projection matching
	
	for N_step in xrange(lstp):
		terminate = 0
		Iter = 0
 		while(Iter < max_iter and terminate == 0):
			yrng[N_step]=float(dp)/(2*pixel_size) #will change it later according to dp
			if(ynumber[N_step]==0):
				stepy=0.0
			else:
				stepy=(2*yrng[N_step]/ynumber[N_step])
				
			pixer  = [0.0]*nima
			modphi = [0.0]*nima
			Iter += 1
			total_iter += 1
			if myid == main_node:
				start_time = time()
				print_msg("\nITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.4f, xrange = %5.4f,stepx = %5.4f, yrange = %5.4f,  stepy = %5.4f, ynumber = %3d\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],stepx[N_step],yrng[N_step],stepy,ynumber[N_step]))
			c1 = range(nima)
			shuffle(c1)
			c1_set = set(c1)
			for nch in xrange(nchunk):
				# data is all the input data of each processor, len(data) = nima
				# rdata is data used for reconstruction, len(rdata) = rnima
				# adata is the data used for aligment, len(data) = anima
				list1 = []
				adata = [None]* anima[nch]
				rdata = [None]*rnima[nch]
	
				for k in xrange( len(adata) ):
					adata[k] = data[c1[chunk_list[2*nch]+k]]
					list1.append( c1[chunk_list[2*nch]+k] )
	
				list1_set = set( list1 )
				list2 = list(  c1_set.difference(list1_set) )
				for k in xrange( len(rdata) ):   
					rdata[k]= data[ list2[k] ]
				if myid == main_node:
					start_time = time()
				if CTF: vol = recons3d_4nn_ctf_MPI(myid, rdata, symmetry=sym, snr = snr, npad = npad, xysize = xysize, zsize = zsize)
				else:    vol = recons3d_4nn_MPI(myid, rdata, symmetry=sym, npad = npad, xysize = xysize, zsize = zsize)
				del rdata, list2
				if myid==main_node:
					#  for each chunk, we directly impose helicla parameters
					vol = vol.helicise(pixel_size,dp, dphi, fract, rmax, rmin)
					ref_data = [vol]
					vol = user_func(ref_data)
					vol = vol.helicise(pixel_size,dp, dphi, fract, rmax, rmin)


				bcast_EMData_to_all(vol, myid, main_node)

				if( xysize == -1 and zsize==-1 ):
					volft,kb = prep_vol( vol )
					refrings = prepare_refrings( volft, kb, nmax, delta[N_step], ref_a, symref, numr, MPI = True, phiEqpsi = "Zero", initial_theta =initial_theta, delta_theta = delta_theta)
					del volft,kb
				else:
					volft, kbx, kby, kbz = prep_vol( vol )
					refrings = prepare_refrings( volft, kbz, nmax, delta[N_step], ref_a, symref, numr, MPI = True, phiEqpsi = "Zero", kbx = kbx, kby = kby, initial_theta =initial_theta, delta_theta = delta_theta)
					del volft, kbx, kby, kbz
				
				#split refrings to two list: refrings1 (theta =90), and refrings2( theat not 90)
				refrings1= []
				refrings2= []
				sn = int(symmetry_string[1:])
				for i in xrange( len(refrings) ):
					if( sn%2 ==0 and abs( refrings[i].get_attr('n3') ) <1.0e-6 and (symmetry_string[0] == "c" or symmetry_string[0] =="d" ) ):
						refrings1.append( refrings[i])

					else:
						refrings2.append( refrings[i])
								
				del refrings
				for im in xrange( anima[nch] ):
					peak1 = None
					peak2 = None
					if ( len(refrings1) > 0):
						if  an[N_step] == -1:
							peak1, phihi1, theta1, psi1, sxi1, syi1, t11 = proj_ali_helical_90(adata[im],refrings1,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step],psi_max,finfo,)
						else:
							peak1, phihi1, theta1, psi1, sxi1, syi1, t11 = proj_ali_helical_90_local(adata[im],refrings1,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step], an[N_step], psi_max,  finfo,)
					if( len(refrings2) > 0):
						if  an[N_step] == -1:
							peak2, phihi2, theta2, psi2, sxi2, syi2, t12 = proj_ali_helical(adata[im],refrings2,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step],psi_max,finfo,)
						else:
							peak2, phihi2, theta2, psi2, sxi2, syi2, t12 = proj_ali_helical_local(adata[im],refrings2,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step], an[N_step], psi_max,finfo,)
					if peak1 is None: 
						peak = peak2
						phihi = phihi2
						theta = theta2
						psi = psi2
						sxi = sxi2
						syi = syi2
						t1 = t12
					elif peak2 is None:
						peak = peak1
						phihi = phihi1
						theta = theta1
						psi = psi1
						sxi = sxi1
						syi = syi1
						t1 = t11
					else:
						if(peak1 >= peak2):
							peak = peak1
							phihi = phihi1
							theta = theta1
							psi = psi1
							sxi = sxi1
							syi = syi1
							t1 = t11
						else:
							peak = peak2
							phihi = phihi2
							theta = theta2
							psi = psi2
							sxi = sxi2
							syi = syi2
							t1 = t12
					#peak, phihi, theta, psi, sxi, syi, t1 = proj_ali_helical(data[im],refrings,numr,xrng[N_step],yrng[N_step],stepx[N_step],ynumber[N_step],psi_max,finfo,)
					if(peak > -1.0e22):
						#Guozhi Tao: wrap y-shifts back into box within rise of one helical unit by changing phi
						jdelta=0.0
						dpp = (float(dp)/pixel_size)

						tp = Transform({"type":"spider","phi":phihi,"theta":theta,"psi":psi})
						tp.set_trans( Vec2f( -sxi, -syi ) )
						trans1 = tp.get_pre_trans()
						dyi = (float(trans1[2])/dpp)-int(trans1[2]/dpp)

						if dyi < -0.5-jdelta:  eyi = dyi+1.0
						elif dyi > 0.5+jdelta:  eyi = dyi-1.0
						else:                   eyi = dyi

						nperiod = float(eyi*dpp-trans1[2])/dpp
						th = Transform({"type":"spider","phi": -nperiod*dphi, "tz":nperiod*dpp})
						tfinal = tp*th
						ddd = tfinal.get_params("spider")
						sxnew = - ddd["tx"]
						synew = - ddd["ty"]
						phinew  = ddd["phi"]

						phihi = phinew
						sxi   = sxnew
						syi   = synew

						# unique range identified by [k0,k1], [k2,k3]
						tp = Transform({"type":"spider","phi":phihi,"theta":theta,"psi":psi})
						tp.set_trans( Vec2f( -sxi, -syi ) )
						k0 = 0.0
						k2 = k0+180

						if( abs( tp.at(2,2) )<1.0e-6 ):
							if (symmetry_string[0] =="c"):
								if sn%2 == 0:
									k1=360.0/sn
								else:
									k1=360.0/2/sn
							elif (symmetry_string[0] =="d"):
								if sn%2 == 0:
									k1=360.0/2/sn
								else:
									k1=360.0/4/sn
						else:
							k1=360.0/sn
						k3 = k1 +180
						from utilities import get_sym
						T = get_sym(symmetry_string[0:])

						for i in xrange( len(T) ):
							ttt = tp*Transform({"type":"spider","phi":T[i][0],"theta":T[i][1],"psi":T[i][2]})
							d1 = ttt.get_params("spider")

							if ( abs( tp.at(2,2) )<1.0e-6 ):
								if( sn%2==1 ): # theta=90 and n odd, only one of the two region match

									if( ( d1['phi'] <= float(k1) and d1['phi'] >= float(k0) ) or ( d1['phi'] < float(k3) and d1['phi'] >= float(k2) )):

										sxnew = - d1["tx"]
										synew = - d1["ty"]
										phinew = d1['phi']
										thetanew = d1["theta"]
										psinew = d1["psi"]
								else: #for theta=90 and n even, there is no mirror version during aligment, so only consider region [k0,k1]

									if( d1['phi'] <= float(k1) and d1['phi'] >= float(k0) ) :

										sxnew = - d1["tx"]
										synew = - d1["ty"]
										phinew = d1['phi']
										thetanew = d1["theta"]
										psinew = d1["psi"]

							else: #theta !=90, # if theta >90, put the projection into [k2,k3]. Otherwise put it into the region [k0,k1]
								if( sn==1):
									sxnew = sxi
									synew = syi
									phinew = phihi
									thetanew = theta
									psinew = psi
								else:

									if (tp.at(2,2) >0.0): #theta <90

										if( d1['phi'] <= float(k1) and d1['phi'] >= float(k0) ):
											if( cos( pi*float( d1['theta'] )/180.0 )>0.0 ):

												sxnew = - d1["tx"]
												synew = - d1["ty"]
												phinew = d1['phi']
												thetanew = d1["theta"]
												psinew = d1["psi"]

									else:
										if(  d1['phi'] <= float(k3) and d1['phi'] >= float(k2) ):
											if( cos( pi*float( d1['theta'] )/180.0 )<0.0 ):

												sxnew = - d1["tx"]
												synew = - d1["ty"]
												phinew = d1['phi']
												thetanew = d1["theta"]
												psinew = d1["psi"]

							del ttt,d1


						t2 = Transform({"type":"spider","phi":phinew,"theta":thetanew,"psi":psinew})
						t2.set_trans(Vec2f(-sxnew, -synew))
						adata[im].set_attr("xform.projection", t2)
						pixer[list1[im]]  = max_3D_pixel_error(t1, t2, numr[-3])
						modphi[list1[im]] = phinew
					else:
						# peak not found, parameters not modified
						pixer[list1[im]]  = 0.0
						phihi, theta, psi, sxi, syi = get_params_proj(adata[im])
						modphi[list1[im]] = phihi
				del adata, list1, list1_set
				del refrings1, refrings2
				
			#output pixel errors after all chunks are used
			from mpi import mpi_gatherv
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			del pixer		
			mpi_barrier(MPI_COMM_WORLD)
			terminate = 0
			if(myid == main_node):
				recvbuf = map(float, recvbuf)
				from utilities import write_text_file
				write_text_file([range(len(recvbuf)), recvbuf], os.path.join(outdir, "pixer_%04d_%04d.txt"%(N_step+1,Iter)) )
				from statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				if(region[0] < 0.0):  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.2f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				# Terminate if 95% within 1 pixel error
				im = 0
				for lhx in xrange(lhist):
					if(region[lhx] > 1.0): break
					im += histo[lhx]
				#if(im/float(total_nima) > 0.95):  terminate = 1
				del region, histo
			#output distribution of phi
			#jeanmod
			recvbuf = mpi_gatherv(modphi, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			#end jeanmod
			mpi_barrier(MPI_COMM_WORLD)
			del modphi
			if(myid == main_node):
				recvbuf = map(float, recvbuf)
				lhist = 30
				region, histo = hist_list(recvbuf, lhist)
				msg = "\n      Distribution of phi\n      phi         number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.2f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				del region, histo
			del recvbuf
			terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
			terminate = int(terminate[0])
			if myid == main_node:
				print_msg("Time to compute pixer = %d\n"%(time()-start_time))
				start_time = time()
			# write out headers, under MPI writing has to be done sequentially
			mpi_barrier(MPI_COMM_WORLD)
			m = 5
			from mpi import mpi_recv, mpi_send, MPI_TAG_UB, MPI_COMM_WORLD, MPI_FLOAT
			if myid == main_node:
				
				fexp = open(os.path.join(outdir, "parameters_%04d_%04d.txt"%(N_step+1,Iter)),"w")
				for n in xrange(number_of_proc):
					if n!=main_node:
						t = mpi_recv(recvcount[n]*m,MPI_FLOAT, n, MPI_TAG_UB, MPI_COMM_WORLD)
						for i in xrange(recvcount[n]):
							for j in xrange(m):
								fexp.write(" %15.5f  "%t[j+i*m])
							fexp.write("\n")
					else:
						t = [0.0]*m
						for i in xrange(recvcount[myid]):
							t = get_params_proj(data[i])
							for j in xrange(m):
								fexp.write(" %15.5f  "%t[j])
							fexp.write("\n")
				fexp.close()
				del t
	        	else:
				nvalue = [0.0]*m*recvcount[myid]
				t = [0.0]*m
				for i in xrange(recvcount[myid]):
					t = get_params_proj(data[i])
					for j in xrange(m):
						nvalue[j + i*m] = t[j]
				mpi_send(nvalue, recvcount[myid]*m, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)
				del nvalue
			if myid == main_node:
				print_msg("Time to write parameters = %d\n"%(time()-start_time))
				start_time = time()

			if CTF: vol = recons3d_4nn_ctf_MPI(myid, data, symmetry=sym, snr = snr, npad = npad, xysize = xysize, zsize = zsize )
			else:    vol = recons3d_4nn_MPI(myid, data, symmetry=sym, npad = npad, xysize = xysize, zsize = zsize )

			if myid == main_node:
				print_msg("\n3D reconstruction time = %d\n"%(time()-start_time))
				start_time = time()

			if fourvar:
			#  Compute Fourier variance
				for im in xrange(nima):
					original_data[im].set_attr( 'xform.projection', data[im].get_attr('xform.projection') )
				varf = varf3d_MPI(original_data, ssnr_text_file = os.path.join(outdir, "ssnr%04d"%(total_iter)), mask2D = None, reference_structure = vol, ou = last_ring, rw = 1.0, npad = 1, CTF = CTF, sign = 1, sym =sym, myid = myid)
				if myid == main_node:
					print_msg("Time to calculate 3D Fourier variance= %d\n"%(time()-start_time))
					start_time = time()
					varf = 1.0/varf
					drop_image(varf, os.path.join(outdir, "varf%04d.hdf"%(total_iter)))
			else:  varf = None

			#search for helical symmetry
			if myid == main_node:
				drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(total_iter)))
							
			if(total_iter > nise):
				bcast_EMData_to_all(vol, myid, main_node)
				#from filter import filt_gaussl
				#vol = filt_gaussl(vol, 0.25)

				if myid == main_node:
					lprms = []
					for i in xrange(-ndp,ndp+1,1):
						for j in xrange(-ndphi,ndphi+1,1):
							lprms.append( dp   + i*dp_step)
							lprms.append( dphi + j*dphi_step)
					#print "lprms===",lprms
					recvpara = []
					for im in xrange(number_of_proc):
						helic_ib, helic_ie = MPI_start_end(nlprms, number_of_proc, im)
						recvpara.append(helic_ib )
						recvpara.append(helic_ie )

				para_start, para_end = MPI_start_end(nlprms, number_of_proc, myid)

				list_dps     = [0.0]*((para_end-para_start)*2)
				list_fvalues = [-1.0]*((para_end-para_start)*1)

				if myid == main_node:
					for n in xrange(number_of_proc):
						if n!=main_node: mpi_send(lprms[2*recvpara[2*n]:2*recvpara[2*n+1]], 2*(recvpara[2*n+1]-recvpara[2*n]), MPI_FLOAT, n, MPI_TAG_UB, MPI_COMM_WORLD)
						else:    list_dps = lprms[2*recvpara[2*0]:2*recvpara[2*0+1]]
				else:
					list_dps = mpi_recv((para_end-para_start)*2, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)

				list_dps = map(float, list_dps)

				local_pos = [0.0, 0.0, -1.0e20]
				for i in xrange(para_end-para_start):
					fvalue = helios7(vol, pixel_size, list_dps[i*2], list_dps[i*2+1], fract, rmax, rmin)
					if(fvalue >= local_pos[2]):
						local_pos = [list_dps[i*2], list_dps[i*2+1], fvalue ]
				if myid == main_node:
					list_return = [0.0]*(3*number_of_proc)
					for n in xrange(number_of_proc):
						if n != main_node: list_return[3*n:3*n+3]                 = mpi_recv(3,MPI_FLOAT, n, MPI_TAG_UB, MPI_COMM_WORLD)
 						else:              list_return[3*main_node:3*main_node+3]  = local_pos[:]
				else:
					mpi_send(local_pos, 3, MPI_FLOAT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)

				if myid == main_node:	
					maxvalue = list_return[2]
					for i in xrange(number_of_proc):
						if( list_return[i*3+2] >= maxvalue ):
							maxvalue = list_return[i*3+2]
							dp       = list_return[i*3+0]
							dphi     = list_return[i*3+1]
					dp   = float(dp)
					dphi = float(dphi)
					#print  "  GOT dp dphi",dp,dphi

					vol  = vol.helicise(pixel_size,dp, dphi, fract, rmax, rmin)
					print_msg("New delta z and delta phi      : %s,    %s\n\n"%(dp,dphi))		
			
			else:
				if myid==main_node:
					#  in the first nise steps the symmetry is imposed
					vol = vol.helicise(pixel_size,dp, dphi, fract, rmax, rmin)
					print_msg("Imposed delta z and delta phi      : %s,    %s\n\n"%(dp,dphi))
			if(myid==main_node):
				fofo = open(os.path.join(outdir,datasym),'a')
				fofo.write('  %12.4f   %12.4f\n'%(dp,dphi))
				fofo.close()
				ref_data = [vol]
				if  fourvar:  ref_data.append(varf)
				vol = user_func(ref_data)
				vol = vol.helicise(pixel_size,dp, dphi, fract, rmax, rmin)

				drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(total_iter)))
				print_msg("\nSymmetry search and user function time = %d\n"%(time()-start_time))
				start_time = time()

			bcast_EMData_to_all(vol, myid, main_node)
			dp   = bcast_number_to_all(dp,   source_node = main_node)
			dphi = bcast_number_to_all(dphi, source_node = main_node)
			#
			del varf
	par_str = ["xform.projection"]
	if myid == main_node:
	   	if(file_type(stack) == "bdb"):
	        	from utilities import recv_attr_dict_bdb
	        	recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        else:
	        	from utilities import recv_attr_dict
	        	recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
		print_msg("Time to write header information= %d\n"%(time()-start_time))
		start_time = time()
	else:	       send_attr_dict(main_node, data, par_str, image_start, image_end)
	if myid == main_node: print_end_msg("ihrsr_chunk_MPI")

def copyfromtif(indir, outdir=None, input_extension="tif", film_or_CCD="f", output_extension="hdf", contrast_invert=1, Pixel_size=1, scanner_param_a=1,scanner_param_b=1, scan_step=63.5, magnification=40, MPI=False):
	"""
		purpose: convert raw image into SPIDER format images
		1. film_or_CCD =f (scanned film, density will be converted into OD )or c (CCD)
		2. contrast_invert =1 or -1  invert contrast or not 
		3. The default vale of input_extension is tiff.
		4. Pixel_size denotes the final converted image pixel size. When it set as 
		   negative, the program will switch back to image_decimate reduce image size
		   integer times. 
	"""
	if MPI:
		copyfromtif_MPI(indir, outdir, input_extension, film_or_CCD, output_extension, contrast_invert, Pixel_size, scanner_param_a, scanner_param_b, scan_step, magnification)
		return

	from utilities 		import get_image, drop_image
	from fundamentals 	import smallprime, window2d, resample, image_decimate
	from filter 		import filt_btwl
	import types
	import os
	if os.path.exists(indir) is False: ERROR("Input directory doesn't exist","copyfromtif",1)
	else                             : flist          = os.listdir(indir)
	if(type(outdir)          is types.StringType):
		if os.path.exists(outdir) :   ERROR("Output directory exists, please change the name and restart the program","copyfromtif",1)
		os.mkdir(outdir)	
	else: 	
		outdir = ("micrographs")# default output directory
		os.mkdir(outdir)
		
	gridding       = False
	Pixel_size_raw = scan_step/magnification
	if Pixel_size == 0 : Pixel_size = Pixel_size_raw
	if Pixel_size <  0 : 
		scaling_ratio  = - int(Pixel_size)
		gridding       = False
		Pixel_size     = Pixel_size_raw*scaling_ratio
	else               : scaling_ratio   = Pixel_size/Pixel_size_raw	
	e              = EMData()
	e1             = EMData()
	X19            = 2**16-1.0
	for i, v in enumerate(flist):
		tifname                 = os.path.join(indir, v)
		(rawname, rawextension) = os.path.splitext(v)
		if(rawextension == "."+input_extension):
			print "The raw file under processing is", tifname
			e = get_image(tifname)
			if(film_or_CCD == "f"):
				e  -= 1.
				e  /= X19
				e1  = e.log10()
				if(scanner_param_a == 1):
					e1  *= -1.883	# scanned by the scanner in Pawel's lab
					e1  += 0.06698
				elif(scanner_param_a == 2):
					e1  *= -1.6090  # scanned by the scanner in Dowhan's lab
					e1  += 0.012
				else:
					e1  *= scanner_param_a # Check the local scanner
					e1  += scanner_param_b
			else:   e1   =e*contrast_invert
			#if  gridding : e1    = resample(e, scaling_ratio, 1) # resample will pad image to four times 
			#else         : e1    = image_decimate(e, scaling_ratio, 1)
			f_micrograph         = "micrograph_"+rawname+"."+ output_extension
			f_micname            = os.path.join(outdir, f_micrograph)
			if output_extension == "spi": drop_image(e1,f_micname,"s")
			else:   e1.write_image(f_micname)

def copyfromtif_MPI(indir, outdir=None, input_extension="tif", film_or_CCD="f", output_extension="hdf", contrast_invert=1, Pixel_size=1, scanner_param_a=1,scanner_param_b=1, scan_step=63.5, magnification=40):
	"""
		purpose: convert raw image into SPIDER format images
		1. film_or_CCD =f (scanned film, density will be converted into OD )or c (CCD)
		2. contrast_invert =1 or -1  invert contrast or not 
		3. The default vale of input_extension is tiff.
		4. Pixel_size denotes the final converted image pixel size. When it set as 
		   negative, the program will switch back to image_decimate reduce image size
		   integer times. 
	"""
	from utilities 		import get_image, drop_image
	from fundamentals 	import smallprime, window2d, resample, image_decimate
	from filter 		import filt_btwl
	from random     	import randint
	import types
	import os
	import sys
	from mpi 	    import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	    import mpi_barrier, mpi_bcast
	from mpi 	    import MPI_INT

	if os.path.exists(indir) is False: ERROR("Input directory doesn't exist","copyfromtif_MPI",1,myid)
	else                             : flist = os.listdir(indir)
	if(type(outdir)          is types.StringType):
		if os.path.exists(outdir): ERROR("Output directory exists, please change the name and restart the program","copyfromtif_MPI",1,myid)
		os.mkdir(outdir)	
	else: 	
		os.system(" rm -rf micrograph ")
		outdir ="micrograph"# default output directory
		os.mkdir(outdir)	
	gridding       = True
	Pixel_size_raw = scan_step/magnification
	if Pixel_size == 0 : Pixel_size = Pixel_size_raw
	if Pixel_size <  0 : 
		scaling_ratio  = - int(Pixel_size)
		gridding       = False
		Pixel_size     = Pixel_size_raw*scaling_ratio
	else               : scaling_ratio   = Pixel_size/Pixel_size_raw
	#	
	nima           = 0
	mic_name_list  = []
	for i, v in enumerate(flist):
		micname  		= os.path.join(indir,v)
		(filename, filextension)= os.path.splitext(v  )
		if(filextension == "."+input_extension):
		       mic_name_list.append(micname)
		       nima += 1
	if nima < 1:    ERROR("No micrograph is found, check either directory or prefix of micrographs is correctly given","copyfromtif_MPI",1,myid)
       
	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid	      = mpi_comm_rank(MPI_COMM_WORLD)
       #  chose a random node as a main one...
	main_node      = 0
	if(myid == 0): main_node = randint(0,number_of_proc-1)
	main_node      = mpi_bcast(main_node, 1, MPI_INT, 0, MPI_COMM_WORLD)

	#nimage_per_node = nima/number_of_proc
	#image_start     = myid * nimage_per_node
	#if(myid == number_of_proc-1):  image_end = nima
	#else:			      image_end = image_start + nimage_per_node
	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)

	e              = EMData()
	e1             = EMData()
	X19            = 2**16-1.0
	#
	for i in xrange(image_start,image_end):
		filename               = mic_name_list[i]
		(tifname, filextension)= os.path.splitext(filename)
		print '%-30s%-30s'%("micrographs under proccessing : ",filename[len(indir)+1:])
		e = get_image(filename)
		if(film_or_CCD == "f"):
			e  -= 1.
			e  /= X19
			e1  = e.log10()
			if(scanner_param_a == 1):
				e1  *= -1.883	# scanned by the scanner in Pawel's lab
				e1  += 0.06698
			elif(scanner_param_a == 2):
				e1  *= -1.6090  # scanned by the scanner in Dowhan's lab
				e1  += 0.012
			else:
				e1  *= scanner_param_a # Check the local scanner
				e1  += scanner_param_b
		e                    = e1*contrast_invert
		if  gridding : e1    = resample(e, scaling_ratio) # resample will pad image two times 
		else         : e1    = image_decimate(e, scaling_ratio, 1)
		f_micrograph         = "micrograph_"+tifname[len(indir)+1:]+"."+ output_extension
		f_micname            = os.path.join(outdir, f_micrograph)
		e1.write_image(f_micname)

# reworked cpy is able to process lists instead of single files. ins_list can
#    be either single file names like "test.hdf" or "bdb:image1", but also lists
#    of these names. note that mixed lists (containing both raw file names and
#    db objects) also work.

def cpy(ins_list, ous):

	# reworked to include lists, since we want to be able to copy lists of images
	#    into one file, concatenating.
	if isinstance(ins_list,list):
		# got a list of input files
		image_list = ins_list
	else:
		# got a single input file
		image_list = [ins_list,]

	gl_index = 0

	from utilities import file_type

	oextension = file_type(ous)	

	if oextension == "bdb":
		from EMAN2db import db_open_dict
		DB = db_open_dict(ous)

	# iterate over all images in the list, even if it's only one...
	for ins in image_list:

		#print ins
		nima = EMUtil.get_image_count(ins)
		data = EMData()
		iextension = file_type(ins)

		if iextension == "bdb":
			from EMAN2db import db_open_dict

		if nima == 1 and oextension == "spi":
			data.read_image(ins)
			data.write_image(ous, 0, EMUtil.ImageType.IMAGE_SINGLE_SPIDER)
			
		elif iextension == "bdb" and oextension == "bdb":
			
			OB = db_open_dict(ins)
			for i in range(nima):
				DB[gl_index] = OB[i]
				gl_index += 1
			OB.close()

		elif iextension == "bdb":
			
			DB = db_open_dict(ins)
			for i in range(nima):
				a = DB[i]
				a.write_image(ous, gl_index)
				gl_index += 1
			DB.close()
			
		elif oextension == "bdb":
			
			for i in range(nima):
				a = EMData()
				a.read_image(ins, i)
				DB[gl_index] = a
				gl_index += 1
			
		else:
			for im in xrange(nima):
				data.read_image(ins, im)
				data.write_image(ous, gl_index)
				gl_index += 1

	if oextension == "bdb":
		DB.close()

def dele_flist(flist):
	""" 
		delete a series files listed in a document file
	"""
	delist=[]
	inf = file(flist, "r")
	strg=inf.readline()
	while (len(strg)>0):
		sh_com="rm -f"+strg
		print sh_com
		delist.append(sh_com)
		strg=inf.readline()
	for i in xrange(len(delist)):  os.system(delist[i])
		
def defocus_calc(roodir, method, writetodoc="w", Pixel_size=1, voltage=120, Cs=1, amp_contrast=.1, round_off=100, dz_max=50000., frequency_low=30, frequency_high=5, polynomial_rank_baseline=5, polynomial_rank_envelope=5, prefix="roo", format="spider", skip_comment="#", micdir = "no", print_screen="no"):	
	from morphology import defocus_get_slow, defocus_get_fast
	if( method == "s"): 	defocus_get_slow(roodir, writetodoc, Pixel_size, voltage, Cs, amp_contrast, round_off, dz_max, frequency_low, frequency_high, prefix, format, skip_comment, micdir, print_screen)
	else: 			defocus_get_fast(roodir, writetodoc, Pixel_size, voltage, Cs, amp_contrast, round_off, dz_max, frequency_low, frequency_high, polynomial_rank_baseline,polynomial_rank_envelope, prefix, format, skip_comment,micdir, print_screen)

'''
def iso_kmeans(images, out_dir, parameter, K=None, mask=None, init_method="Random"):
	from statistics import init_Kmeans,Kmeans_step,kmeans_ave_var,iso_kmeans_rm_cluster,iso_kmeans_split,iso_kmeans_merge
	import os
	
	e=EMData()
	if(type(images) == str):
		flag = 1
		N = EMUtil.get_image_count(images)
		e.read_image(images,0)
	else:
		if(type(images) == list):
			flag = 0
			N = len(images)
			e = images[0].copy()
		else:
			return 'Invalid Input file format --- The input images can only be a stack file or a list of images'

	nx = e.get_xsize()
	ny = e.get_ysize()
	nz = e.get_zsize()
	size = nx*ny*nz
	if(mask == None):
		mask = model_blank(nx, ny, nz, 1.0)

	class Cluster:
		def __init__(self,a,n,C,SSE,V):
			self.a   = a
			self.n   = n
			self.C   = C
			self.SSE = SSE
			self.V   = V

	m = "%s" % (parameter)
	parameter *= size

	out  = open(out_dir+"/Isodata_kmeans_chart_"+m+".txt",'w')

	[assign,TE,Cls] = init_Kmeans(images, N, K, e, flag, mask, init_method)

	[assign,TE,Cls] = Kmeans_step(images,N,len(Cls),assign,TE,Cls,size,mask,flag)

	out.write("%s\t%d\n" % ("After K means the Total Error is",TE))

	for k in xrange(len(Cls)):
		out.write("\n%s\t%d\t%s\t%g\t%s\t%d\n" % ("Cluster",k,"variance = ",Cls[k].V/size,"No of Object = ",Cls[k].n))
		out.write("\n%s\t%d\t\t%s\n" % ("Cluster",k,"assignment"))
		Cls[k].a.sort()
		for kk in xrange(len(Cls[k].a)):
			out.write("\t%d" % (Cls[k].a[kk]))
			if(float(kk+1)/10==int(kk+1)/10):
				out.write("\n")
		out.write("\n")

	cnt = 0
	change = True
	while(change):

		print " ITERATION #",cnt,"  Total Error = ",TE
		change=False

		[assign,Cls] = iso_kmeans_rm_cluster(images,N,assign,Cls)

		[assign,TE,Cls] = iso_kmeans_split(images,N,assign,TE,Cls,Cluster,size,parameter,mask,init_method,flag)

		length_split = len(Cls)

		[assign,Cls] = iso_kmeans_rm_cluster(images,N,assign,Cls)

		[assign,TE,Cls] = iso_kmeans_merge(images,N,assign,TE,Cls,Cluster,size,parameter,mask,init_method,flag)

		length_merge =  len(Cls)
		
		if(length_merge != length_split):
			change = True
		cnt +=1

	out.write("%s\t%d\n" % ("After Isodata Clustering the Total Error is",TE))
		
	for k in xrange(len(Cls)):
		out.write("\n%s\t%d\t%s\t%g\t%s\t%d\n" % ("Cluster",k,"variance = ",Cls[k].V/size,"No of Object = ",Cls[k].n))
		out.write("\n%s\t%d\t\t%s\n" % ("Cluster",k,"assignment"))
		Cls[k].a.sort()
		for kk in xrange(len(Cls[k].a)):
			out.write("\t%d" % (Cls[k].a[kk]))
			if(float(kk+1)/10==int(kk+1)/10):
				out.write("\n")
		out.write("\n")
	
	os.system('rm -f '+out_dir+'/Isodata_kmeans_average_'+m+'.spi')
	os.system('rm -f '+out_dir+'/Isodata_kmeans_variance_'+m+'.spi')
	for k in xrange(len(Cls)):
		res = kmeans_ave_var(Cls[k].a,images,flag)
		res[0].write_image(out_dir+"/Isodata_kmeans_average_"+m+".spi",k)
		res[1].write_image(out_dir+"/Isodata_kmeans_variance_"+m+".spi",k)
'''

def project3d(volume, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None):
	from projection    import   prgs, prep_vol
	from utilities     import   even_angles, read_text_row, set_params_proj, model_gauss_noise, info
	from string        import   split
	from filter        import   filt_ctf,filt_gaussl
	import os
	import types

	if listagls is None:
		angles = even_angles(delta, symmetry = symmetry, method = method, phiEqpsi = phiEqpsi)
	elif(type(listagls) is types.StringType):
		angles = read_text_row(listagls, "", "")
	else:
		angles = listagls

	# try to parse the CTFs list. this is either not set (None), a filename or a list of values
	if listctfs is None:
		# not set, so simply ignore it 
		pass
	elif (type(listctfs) is types.StringType):
		# a string, so assume this is a filename and try to open the file
		try:
			ctfs = read_text_row(listctfs, "", "")
		except:
			ctfs = [None for ii in xrange(len(angles))]
	else:
		# assume this a list of len(angles)
		ctfs = listctfs

	if not noise is None:
		# try to convert noise string to float. ignore noise if this fails
		try:
			noise_level = float(noise)
		except:
			noise_level = None
	# ignore noise, since it was not requested
	else:
		noise_level = None

	if(type(volume) is types.StringType):
		vol = EMData()
		vol.read_image(volume)
		if(mask):
			if(type(mask) is types.StringType):
				maski = EMData()
				maski.read_image(volume)
				Util.mul_img(vol, maski)
				del maski
			else:
				Util.mul_img(vol, mask)
		nx = vol.get_xsize()
		ny = vol.get_ysize()
		nz = vol.get_zsize()
		
		if(nx==nz&ny==nz):
			volft, kb = prep_vol(vol)
		else:
			volft, kbx,kby,kbz = prep_vol(vol)
	else:
		if(mask):
			if(type(mask) is types.StringType):
				maski = EMData()
				maski.read_image(volume)
				Util.mul_img(vol, maski)
				del maski
			else:
				Util.mul_img(vol, mask)
		nx = vol.get_xsize()
		ny = vol.get_ysize()
		nz = vol.get_zsize()
		
		if(nx==nz&ny==nz):
			volft, kb = prep_vol(vol)
		else:
			volft, kbx,kby,kbz = prep_vol(vol)

	if(type(stack) is types.StringType):
		Disk = True
		os.system("rm -f  "+stack)	
	else:
		out = []
		Disk = False
	
	s2x=0
	s2y=0
	
	for i in xrange(len(angles)):
		if(len(angles[i]) == 3):
			if(nx==nz&ny==nz):
				proj = prgs(volft, kb, [angles[i][0], angles[i][1], angles[i][2], 0.0, 0.0])
			else:
				proj = prgs(volft, kbz, [angles[i][0], angles[i][1], angles[i][2], 0.0, 0.0],kbx,kby)
			set_params_proj(proj, [angles[i][0], angles[i][1], angles[i][2], 0.0, 0.0])
		else:
			if(nx==nz&ny==nz):
				proj = prgs(volft, kb, [angles[i][0], angles[i][1], angles[i][2], -angles[i][3], -angles[i][4]])
			else:
				proj = prgs(volft, kbz, [angles[i][0], angles[i][1], angles[i][2], -angles[i][3], -angles[i][4]],kbx,kby)
			set_params_proj(proj, angles[i])
		proj.set_attr_dict({'active':1})

		# add noise, if noise is set. this is two-fold: application of noise before
		#    ctf filtering and after it.
		if noise is not None:
			try:
				# no mask, so call w/ false
				noise_ima = model_gauss_noise(noise_level,proj.get_xsize(),
							      proj.get_ysize())
			except:
				pass
			else:
				proj += noise_ima

		# apply ctf, if ctf option is set and if we can create a valid CTF object
		try:
			ctf = EMAN2Ctf()
			# order of values is the one applied in sxheader for import / export!
			ctf.from_dict({ "defocus":ctfs[i][0], "cs":ctfs[i][1], "voltage":ctfs[i][2], 
					"apix":ctfs[i][3], "bfactor":ctfs[i][4], "ampcont":ctfs[i][5] })
		except:
			# there are no ctf values, so ignore this and set no values
			proj.set_attr( "error",1)
		else:
			# setting of values worked, so apply ctf and set the header info correctly
			proj = filt_ctf(proj,ctf)
			proj.set_attr( "ctf",ctf)
			proj.set_attr( "ctf_applied",0)

		# add second noise level that is not affected by CTF
		if noise is not None:
			try:
				noise_ima = model_gauss_noise(noise_level,proj.get_xsize(),
							      proj.get_ysize())
			except:
				pass
			else:
				proj += noise_ima

		if(Disk):
			proj.write_image(stack, i)
		else: 
			out.append(proj)
	if(not Disk):  return out

def pw2sp(indir, outdir = None, w =256, xo =50, yo = 50, xd = 0, yd = 0, r = 0, prefix_of_micrograph="micrograph", MPI=False):
	""" 
		Purpose: 
		Calculate power spectra of a list of micrographs in a given directory using Welch's periodogram
		The input options enable one selects area in micrographs to calculate overlapped periodogram.
	"""
	if MPI:
		pw2sp_MPI(indir, outdir, w, xo, yo, xd, yd, r, prefix_of_micrograph)
		return

	from utilities    import get_image,drop_image, model_circle, info
	from fundamentals import welch_pw2, ro_textfile
	import sys
	import os
	import types
	if os.path.exists(indir) is False: ERROR("Input directory doesn't exist","pw2sp",1)
	else				 : flist    = os.listdir(indir)
	if(type(outdir)          is types.StringType):
		if os.path.exists(outdir) is True: ERROR("Output directory exists, please change the name and restart the program","pw2sp",1)
		os.mkdir(outdir)	
	else: 	os.system("mkdir power")
	mask     = model_circle(int(r), int(w), int(w), nz=1)
	mask    -= 1
	mask    *= -1
	pw2_mask = EMData()
	ncount   = 0
	for i, v in enumerate(flist):
		micname=os.path.join(indir,v)
		(filename, filextension) = os.path.splitext(v)
		if(filename[0:len(prefix_of_micrograph)] == prefix_of_micrograph):
			print "micrographs :",filename
			ncount   += 1
			roofile   = "roo"+filename[len(prefix_of_micrograph):]+".txt"
			pw2file   = "pw2"+filename[len(prefix_of_micrograph):]+filextension
			e         = get_image(micname)
 			pw2       = welch_pw2(e,int(w),int(xo),int(yo),int(xd),int(yd))
			pw2_mask  = pw2*mask
			pw2name   = os.path.join(outdir,pw2file)
			drop_image(pw2_mask, pw2name)
			rotxtname = os.path.join(outdir,roofile)
			ro_textfile(pw2, rotxtname)
	if ncount < 1: 	ERROR("No micrograph is found, check either directory or prefix of micrographs is correctly given","pw2sp",1)
 
def pw2sp_MPI(indir, outdir, w =256, xo =50, yo = 50, xd = 0, yd = 0, r = 0, prefix_of_micrograph="micrograph"):
	""" 
		Purpose: 
		Calculate power spectra of a list of micrographs in a given directory using Welch's periodogram
		The input options enable one selects area in micrographs to calculate overlapped periodogram.
	"""
	from utilities    	import get_image,drop_image, model_circle, info
	from fundamentals 	import welch_pw2, ro_textfile
	from random     	import randint
	import sys
	import os
	import types
	from mpi 	    import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	    import mpi_barrier, mpi_bcast
	from mpi 	    import MPI_INT

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	#  chose a random node as a main one...
	main_node      = 0
	if(myid == 0): main_node = randint(0,number_of_proc-1)
	main_node      = mpi_bcast(main_node, 1, MPI_INT, 0, MPI_COMM_WORLD)	
	
	if os.path.exists(outdir): ERROR("Output directory exists, please change the name and restart the program","pw2sp_MPI ",1,myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == int(main_node):  # only main node do cleaning & creating jobs
		print_begin_msg("pw2sp_MPI")
		os.mkdir(outdir)
	mpi_barrier(MPI_COMM_WORLD)
	
	# get the total micrograph number 
	if os.path.exists(indir) is False: ERROR("Input directory doesn't exist","pw2sp_MPI",1,myid)
	else:	                           flist = os.listdir(indir)
	nima          = 0
	mic_name_list = []
	for i, v in enumerate(flist):
		micname                  = os.path.join(indir,v)
		(filename, filextension) = os.path.splitext(v)
		if(filename[0:len(prefix_of_micrograph)] == prefix_of_micrograph):
			mic_name_list.append(micname)
			nima += 1
	if nima < 1: 	ERROR("No micrograph is found, check spelling of either directory or micrographs","pw2sp_MPI",1,myid)	

	#if nima < number_of_proc: nimage_per_node = 1
	#else                    : nimage_per_node = nima/number_of_proc
	#image_start     = myid *  nimage_per_node
	#if(myid == number_of_proc-1):  image_end = nima
	#else:                          image_end = image_start + nimage_per_node				
	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
	if image_start == image_end: 
		if image_end == nima:	image_start -= 1
		else: image_end += 1

	mask     =  model_circle(int(r), int(w), int(w), nz=1)
	mask    -=  1
	mask    *= -1
	pw2_mask =  EMData() 	
	for i in xrange(image_start,image_end):
		filename=mic_name_list[i] 
		print '%-15s%-30s'%("micrographs : ",filename)
		(f_nam, filextension) = os.path.splitext(filename)
		roofile   = "roo_"+f_nam[len(prefix_of_micrograph)+len(indir)+2:]+".txt"
		pw2file   = "pw2_"+f_nam[len(prefix_of_micrograph)+len(indir)+2:]+filextension
		e         = get_image(filename)
 		pw2       = welch_pw2(e,int(w),int(xo),int(yo),int(xd),int(yd))
		pw2_mask  = pw2*mask
		pw2name   = os.path.join(outdir,pw2file)
		drop_image(pw2_mask, pw2name)
		rotxtname = os.path.join(outdir,roofile)
		ro_textfile(pw2, rotxtname)
	
def ra_cef(indir, noise, outdir, prf, num):
	"""
		Apply ramp and cefit to detected particles
		1. Fits a least-squares plane to the picture, 
		   and subtracts the plane from the picture.(ramp) 
		2. Fit the histogram of the input image under 
		   mask with the reference image (ce fit)
	"""
	flist  = os.listdir(indir)
	e      = EMData()
	e_n    = get_image(noise)
	nx     = e_n.get_xsize()
	ny     = e_n.get_ysize()
	radius = nx//2-1
	mask   = model_circle(radius, nx, ny)
	for i, v in enumerate(flist):
		(filename,filextension) = os.path.splitext(v)
		if(filename[0:4] == "ptl_"):
			infile  = os.path.join(indir, v)
			outfile = os.path.join(outdir,prf+v)
			for j in xrange(num):
				e.read_image(infile, j)
				ra=ramp(e)
				e=ce_fit(ra, e_n, mask)
				e.write_image(outfile, j)

def ali_vol_2(vol, refv, ang_scale, shift_scale, radius=None, discrepancy = "ccc"):
	#rotation and shift
	from alignment    import ali_vol_func
	from utilities    import get_im, model_circle
	from utilities    import amoeba
	from fundamentals import rot_shift3D

	nx = refv.get_xsize()
	ny = refv.get_ysize()
	nz = refv.get_zsize()
	if(radius != None):   mask = model_circle(radius, nx, ny, nz)
	else:                 mask = model_circle(float(min(nx, ny, nz)//2-2), nx, ny, nz)

	scale = [ang_scale, ang_scale, ang_scale, shift_scale, shift_scale, shift_scale]
	data  = [vol, refv, mask]
	new_params = [0.0]*6
	new_params = amoeba(new_params, scale, ali_vol_func, 1.e-4, 1.e-4, 500, data)
	vol = rot_shift3D(vol, new_params[0][0], new_params[0][1], new_params[0][2], new_params[0][3], new_params[0][4], new_params[0][5],1.0)
	return vol

def ali_vol_3(vol, refv, ang_scale, shift_scale, radius=None, discrepancy = "ccc", mask=None):
	#rotation and shift
	from alignment    import ali_vol_func
	from utilities    import model_circle, amoeba

	nx = refv.get_xsize()
	ny = refv.get_ysize()
	nz = refv.get_zsize()
	if mask is None:
		if(radius != None): mask = model_circle(radius, nx, ny, nz)
		else:               mask = model_circle(float(min(nx, ny, nz)//2-2), nx, ny, nz)

	scale = [ang_scale, ang_scale, ang_scale, shift_scale, shift_scale, shift_scale]
	data=[vol, refv, mask, discrepancy]
	new_params = [0.0]*6
	opt_params,funval,niter = amoeba(new_params, scale, ali_vol_func, 1.e-4, 1.e-4, 500, data)
	return opt_params

def ali_vol(vol, refv, ang_scale, shift_scale, radius=None, discrepancy = "ccc"):
	"""
		Name
			sxali_vol - align a 3D structure with respect of a 3D reference structure
		Input
			aligned_volume.hdf: 3D structure to be aligned.
		reference_volume.hdf
			3D reference structure.
		ang_scale
			correct angles are expected to be within +/-ang_scale of the values preset in the header of the structure to be aligned
		shift_scale
			correct shifts are expected to be within +/-shift_scale of the values preset in the header of the structure to be aligned
		mag_scale
			correct magnification is expected to be within +/-mag_scale of the value preset in the header of the structure to be aligned
		r
			radius of a spherical mask centered at nx/2, ny/2, nz/2
		Note - there are no defaults for three scale parameters. At least one has to appear.
	"""


	#rotation and shift
	from alignment    import ali_vol_func
	from utilities    import get_image, model_circle, get_params3D, set_params3D
	from utilities    import amoeba, compose_transform3
	from fundamentals import rot_shift3D
	
	ref = get_image(refv)
	nx = ref.get_xsize()
	ny = ref.get_ysize()
	nz = ref.get_zsize()
	if(radius != None):    mask = model_circle(radius, nx, ny, nz)
	else:                  mask = model_circle(float(min(nx, ny, nz)//2-2), nx, ny, nz)

	#names_params = ["phi", "theta", "psi", "s3x", "s3y", "s3z", "scale"]
	phi, theta, psi, s3x, s3y, s3z, mirror, scale = get_params3D(ref)
	paramsr = [phi, theta, psi, s3x, s3y, s3z, mirror, scale]
	# print  " params of the reference volume", paramsr
	ref = rot_shift3D(ref, paramsr[0], paramsr[1], paramsr[2], paramsr[3], paramsr[4], paramsr[5], paramsr[7])

	e = get_image(vol)
	phi, theta, psi, s3x, s3y, s3z, mirror, scale =  get_params3D(e)
	paramsv = [phi, theta, psi, s3x, s3y, s3z, mirror, scale]
	e = rot_shift3D(e, phi, theta, psi, s3x, s3y, s3z, scale)
	# print  " input params ", paramsv
	params = [phi, theta, psi, s3x, s3y, s3z]
	data = [e, ref, mask, params, discrepancy]
	new_params = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	
	new_params = amoeba(new_params, [ang_scale, ang_scale, ang_scale, shift_scale, shift_scale, shift_scale], ali_vol_func, 1.e-4, 1.e-4, 500, data)
	print "amoeba: func_value =",new_params[1], "iter =",new_params[2]

	cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale= compose_transform3(paramsv[0], paramsv[1], paramsv[2], paramsv[3], paramsv[4], paramsv[5], paramsv[7], new_params[0][0], new_params[0][1], new_params[0][2], new_params[0][3], new_params[0][4], new_params[0][5],1.0)
	# print  " new params ", cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale, new_params[1]
	set_params3D(e, [cphi, ctheta, cpsi, cs2x, cs2y, cs2z, 0, cscale])
	if type(vol)==type(""):
		from utilities import write_headers
		write_headers( vol, [e], [0])
	else:
		return e

def ali_vol_n(vol, refv, ang_scale, shift_scale, radius=None, discrepancy="ccc", rsdec=1):
	"""
		Name
			sxali_vol_n - align a 3D structure with respect of a 3D reference structure.
				Like sxali_vol, but bypassing the composition of transformations.
		Input
			aligned_volume.hdf: 3D structure to be aligned.
		reference_volume.hdf
			3D reference structure.
		ang_scale
			correct angles are expected to be within +/-ang_scale of the values preset in the header of the structure to be aligned
		shift_scale
			correct shifts are expected to be within +/-shift_scale of the values preset in the header of the structure to be aligned
		mag_scale
			correct magnification is expected to be within +/-mag_scale of the value preset in the header of the structure to be aligned
		radius
			radius of a spherical mask centered at nx/2, ny/2, nz/2
		rsdec
			(int) if given and >1, the map, after being transformed, is convolved with a gaussian kernel, and then decimated rsdec times.
		Note - there are no defaults for three scale parameters. At least one has to appear.
	"""


	#rotation and shift
	from alignment    import ali_vol_func
	from utilities    import get_image, model_circle
	from utilities    import amoeba, get_params3D, set_params3D
	from utilities    import get_arb_params, set_arb_params
	
	ref = get_image(refv)
	nx = ref.get_xsize()
	ny = ref.get_ysize()
	nz = ref.get_zsize()
	if(radius != None):    mask = model_circle(radius, nx, ny, nz)
	else:                  mask = model_circle(float(min(nx, ny, nz)//2-2), nx, ny, nz)

	names_params = ["phi", "theta", "psi", "s3x", "s3y", "s3z"]

	e = get_image(vol)
	params =  get_arb_params(e, names_params)
	data = [e, ref, mask, None, discrepancy, rsdec]
	
	new_params = amoeba(params, [ang_scale, ang_scale, ang_scale, shift_scale, shift_scale, shift_scale], ali_vol_func, 1.e-5, 1.e-4, 500, data)

	set_arb_params(e, [new_params[0][0], new_params[0][1], new_params[0][2], new_params[0][3], new_params[0][4], new_params[0][5]], names_params)
	if type(vol)==type(""):
		from utilities import write_headers
		write_headers( vol, [e], [0])
	else:
		return e

def ali_vol_grid(vol, params, refv, ang_scale, shift_scale, radius=None, discrepancy="dot", kb=None, wrap=False):
	"""
		Name
			sxali_vol_grid - align a 3D structure with respect of a 3D reference structure.
				Like ali_vol_n, but using gridding transformations.

		Arguments:
		vol
			3D structure to be aligned
		params
			starting parameters to be applied to vol
		refv
			3D reference structure
		ang_scale
			correct angles are expected to be within +/-ang_scale of the values preset in the header of the structure to be aligned
		shift_scale
			correct shifts are expected to be within +/-shift_scale of the values preset in the header of the structure to be aligned
		radius
			radius of a spherical mask centered at nx/2, ny/2, nz/2
		discrepancy
			discrepancy (or similarity) measure
		kb
			if given, then (vol,kb) must be the output from prepi3D. If not given, then prepi3D is called here.
		wrap
			if True, use wraparound pixels during translations
	"""


	#rotation and shift
	from alignment    import ali_vol_func_grid
	from utilities    import get_image, model_circle
	from utilities    import amoeba, get_params3D, set_params3D
	from utilities    import get_arb_params, set_arb_params
	
	ref = get_image(refv)
	nx = ref.get_xsize()
	ny = ref.get_ysize()
	nz = ref.get_zsize()
	if(radius != None):    mask = model_circle(radius, nx, ny, nz)
	else:                  mask = model_circle(float(min(nx, ny, nz)//2-2), nx, ny, nz)

	e = get_image(vol)
	
	if kb==None:
		ee, kb = prepi3D(e)
		data = [ee, ref, mask, None, discrepancy, kb, wrap]
	else:
		data = [e, ref, mask, None, discrepancy, kb, wrap]
	
	new_params = amoeba(params, [ang_scale, ang_scale, ang_scale, shift_scale, shift_scale, shift_scale], ali_vol_func_grid, 1.e-5, 1.e-4, 500, data)

	return [new_params[0][0], new_params[0][1], new_params[0][2], new_params[0][3], new_params[0][4], new_params[0][5]]


def ali_vol_M(vol, refv, ang_scale, shift_scale, mask=None, discrepancy = "ccc"):
	"""
		Name
			sxali_vol_m - align a 3D structure with respect of a 3D reference structure.
				Like sxali_vol_n, but taking a specified mask instead of a radius.
		Input
			aligned_volume.hdf: 3D structure to be aligned.
		reference_volume.hdf
			3D reference structure.
		ang_scale
			correct angles are expected to be within +/-ang_scale of the values preset in the header of the structure to be aligned
		shift_scale
			correct shifts are expected to be within +/-shift_scale of the values preset in the header of the structure to be aligned
		mag_scale
			correct magnification is expected to be within +/-mag_scale of the value preset in the header of the structure to be aligned
		radius
			radius of a spherical mask centered at nx/2, ny/2, nz/2
		Note - there are no defaults for three scale parameters. At least one has to appear.
	"""


	#rotation and shift
	from alignment    import ali_vol_func
	from utilities    import get_image, model_circle
	from utilities    import amoeba, get_params3D, set_params3D
	from utilities    import get_arb_params, set_arb_params
	
	ref = get_image(refv)
	nx = ref.get_xsize()
	ny = ref.get_ysize()
	nz = ref.get_zsize()

	names_params = ["phi", "theta", "psi", "s3x", "s3y", "s3z"]

	e = get_image(vol)

	if(mask == None):
		mask = model_circle(float(min(nx, ny, nz)//2-2), nx, ny, nz)
		minval = None
	elif(mask == "tight"):
		minval = 0.1*Util.infomask(e,None,True)[0]   # threshold for binarizing
	elif(isinstance(mask, (int))):
		minval = -mask    # negative denotes mask is a moving sphere
		mask = model_circle(-minval, nx, ny, nz)
	else:
		minval = None

	params = get_arb_params(e, names_params)
	data = [e, ref, mask, minval, discrepancy]
	
	maxiter = 500
	new_params = amoeba(params, [ang_scale, ang_scale, ang_scale, shift_scale, shift_scale, shift_scale], ali_vol_func, 1.e-5, 1.e-4, maxiter, data)

	if new_params[2]>=maxiter:
		print "Warning: amoeba reached the max number of iterations allowed."

	set_arb_params(e, [new_params[0][0], new_params[0][1], new_params[0][2], new_params[0][3], new_params[0][4], new_params[0][5]], names_params)
	if type(vol)==type(""):
		from utilities import write_headers
		write_headers( vol, [e], [0])
	else:
		return e


def ali_vol_nopsi(vol, refv, ang_scale, shift_scale, radius=None, discrepancy = "ccc"):
	"""
		Name
			sxali_vol_nopsi - align a 3D structure with respect of a 3D reference structure,
				like sxali_vol, but keeping psi fixed at 0.
		Input
			aligned_volume.hdf: 3D structure to be aligned.
		reference_volume.hdf
			3D reference structure.
		ang_scale
			correct angles are expected to be within +/-ang_scale of the values preset in the header of the structure to be aligned
		shift_scale
			correct shifts are expected to be within +/-shift_scale of the values preset in the header of the structure to be aligned
		mag_scale
			correct magnification is expected to be within +/-mag_scale of the value preset in the header of the structure to be aligned
		radius
			radius of a spherical mask centered at nx/2, ny/2, nz/2
		Note - there are no defaults for three scale parameters. At least one has to appear.
	"""


	#rotation and shift
	from alignment    import ali_vol_func_nopsi
	from utilities    import get_image, model_circle
	from utilities    import amoeba, get_params3D, set_params3D
	from utilities    import get_arb_params, set_arb_params
	
	ref = get_image(refv)
	nx = ref.get_xsize()
	ny = ref.get_ysize()
	nz = ref.get_zsize()
	if(radius != None):    mask = model_circle(radius, nx, ny, nz)
	else:                  mask = model_circle(float(min(nx, ny, nz)//2-2), nx, ny, nz)

	names_params = ["phi", "theta", "s3x", "s3y", "s3z"]

	e = get_image(vol)
	params = get_arb_params(e, names_params)
	data = [e, ref, mask, None, discrepancy]
	
	new_params = amoeba(params, [ang_scale, ang_scale, shift_scale, shift_scale, shift_scale], ali_vol_func_nopsi, 1.e-5, 1.e-4, 500, data)

	set_arb_params(e, [new_params[0][0], new_params[0][1], new_params[0][2], new_params[0][3], new_params[0][4]], names_params)
	if type(vol)==type(""):
		from utilities import write_headers
		write_headers( vol, [e], [0])
	else:
		return e


def ali_vol_rotate(vol, refv, ang_scale, radius=None, discrepancy = "ccc"):
	#rotation 
	from alignment    import ali_vol_func_rotate
	from utilities    import get_image, model_circle, get_params3D, set_params3D
	from utilities    import amoeba, compose_transform3
	from fundamentals import rot_shift3D

	ref = get_image(refv)
	nx = ref.get_xsize()
	ny = ref.get_ysize()
	nz = ref.get_zsize()
	if(radius != None):    mask = model_circle(radius, nx, ny, nz)
	else:                  mask = model_circle(float(min(nx, ny, nz)//2-2), nx, ny, nz)

	#names_params = ["phi", "theta", "psi", "s3x", "s3y", "s3z", "scale"]
	params = get_params3D(ref)
	#print  " params of the reference volume",params
	ref = rot_shift3D(ref, params[0], params[1], params[2], params[3], params[4], params[5], params[7])

	e = get_image(vol)
	params = get_params3D(e)
	#e = rot_shift3D(e, params[0], params[1], params[2], params[3], params[4], params[5], params[7])
	#print  " input params ", params
	data = [e, ref, mask, params, discrepancy]
	new_params = [0.0, 0.0, 0.0]
	new_params = amoeba(new_params, [ang_scale, ang_scale, ang_scale], ali_vol_func_rotate, 1.e-4, 1.e-4, 500, data)
	cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale= compose_transform3(params[0], params[1], params[2], params[3], params[4], params[5], params[7], new_params[0][0], new_params[0][1], new_params[0][2],0.0,0.0,0.0,1.0)
	#print  " new params ", cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale, new_params[1]
	set_params3D(e, [cphi, ctheta, cpsi, cs2x, cs2y, cs2z, 0, cscale])
	if type(vol)==type(""):
		from utilities import write_headers
		write_headers( vol, [e], [0])
	else:
		return e

def ali_vol_shift(vol, refv, shift_scale, radius=None, discrepancy = "ccc"):
	# shift
	from alignment    import ali_vol_func_shift
	from utilities    import get_image, model_circle, get_params3D, set_params3D
	from utilities    import amoeba, compose_transform3
	from fundamentals import rot_shift3D

	ref = get_image(refv)
	nx = ref.get_xsize()
	ny = ref.get_ysize()
	nz = ref.get_zsize()
	if(radius != None):    mask = model_circle(radius, nx, ny, nz)
	else:                  mask = model_circle(float(min(nx, ny, nz)//2-2), nx, ny, nz)

	#names_params = ["phi", "theta", "psi", "s3x", "s3y", "s3z", "scale"]
	params = get_params3D(ref)
	#print  " params of the reference volume",params
	ref = rot_shift3D(ref, params[0], params[1], params[2], params[3], params[4], params[5], params[7])

	e = get_image(vol)
	params = get_params3D(e)
	#e = rot_shift3D(e, params[0], params[1], params[2], params[3], params[4], params[5], params[7])
	#print  " input params ",params
	data = [e, ref, mask, params, discrepancy]
	new_params = [0.0, 0.0, 0.0]
	new_params = amoeba(new_params, [shift_scale, shift_scale, shift_scale], ali_vol_func_shift, 1.e-4, 1.e-4, 500, data)
	cphi, ctheta, cpsi, cs3x, cs3y, cs3z, cscale= compose_transform3(params[0], params[1], params[2], params[3], params[4], params[5], params[7], 0.0,0.0,0.0, new_params[0][0], new_params[0][1], new_params[0][2],1.0)
	#print  " new params ", cphi, ctheta, cpsi, cs3x, cs3y, cs3z, cscale, new_params[1]
	set_params3D(e, [cphi, ctheta, cpsi, cs3x, cs3y, cs3z, 0, cscale])
	if type(vol)==type(""):
		from utilities import write_headers
		write_headers( vol, [e], [0])
	else:
		return e

def ali_vol_scale(vol, refv, ang_scale, shift_scale, mag_scale, radius=None, discrepancy = "ccc"):
	# rotation shift and scale
	from alignment    import ali_vol_func_scale
	from utilities    import get_image, model_circle, get_params3D, set_params3D
	from utilities    import amoeba, compose_transform3
	from fundamentals import rot_shift3D
	ref = get_image(refv)
	nx = ref.get_xsize()
	ny = ref.get_ysize()
	nz = ref.get_zsize()
	if(radius != None):    mask = model_circle(radius, nx, ny, nz)
	else:                  mask = model_circle(float(min(nx, ny, nz)//2-2), nx, ny, nz)

	#names_params = ["phi", "theta", "psi", "s3x", "s3y", "s3z", "scale"]
	params = get_params3D(ref)
	print  " params of the reference volume",params
	ref = rot_shift3D(ref, params[0], params[1], params[2], params[3], params[4], params[5], params[7])

	e = get_image(vol)
	params = get_params3D(e)
	print  " input params ",params
	data = [e, ref, mask, params, discrepancy]
	new_params = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
	new_params = amoeba(new_params, [ang_scale, ang_scale, ang_scale, shift_scale, shift_scale, shift_scale, mag_scale], ali_vol_func_scale, 1.e-4, 1.e-4, 500, data)
	cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale= compose_transform3(params[0], params[1], params[2], params[3], params[4], params[5], params[7], new_params[0][0], new_params[0][1], new_params[0][2], new_params[0][3], new_params[0][4], new_params[0][5], new_params[0][6])
	print  " new params ", cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale, new_params[1]
	set_params3D(e, [cphi, ctheta, cpsi, cs2x, cs2y, cs2z, 0, cscale])
	
	if type(vol)==type(""):
		from utilities import write_headers
		write_headers( vol, [e], [0])
	else:
		return e

def ali_vol_only_scale(vol, refv, mag_scale, radius=None, discrepancy = "ccc"):
	# scale
	from alignment    import ali_vol_func_only_scale
	from utilities    import get_image, model_circle, get_params3D, set_params3D
	from utilities    import amoeba, compose_transform3
	from fundamentals import rot_shift3D
	ref = get_image(refv)
	nx = ref.get_xsize()
	ny = ref.get_ysize()
	nz = ref.get_zsize()
	if(radius != None):    mask = model_circle(radius, nx, ny, nz)
	else:                  mask = model_circle(float(min(nx, ny, nz)//2-2), nx, ny, nz)


	#names_params = ["phi", "theta", "psi", "s3x", "s3y", "s3z", "scale"]
	params = get_params3D(ref)
	print  " params of the reference volume",params
	ref = rot_shift3D(ref, params[0], params[1], params[2], params[3], params[4], params[5], params[6])

	e = get_image(vol)
	params = get_params3D(e)
	print  " input params ",params
	data = [e, ref, mask, params, discrepancy]
	new_params = [1.0]
	new_params = amoeba(new_params, [mag_scale], ali_vol_func_only_scale, 1.e-4, 1.e-4, 500, data)
	cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale= compose_transform3(params[0], params[1], params[2], params[3], params[4], params[5], params[7], 0.0,0.0,0.0,0.0,0.0,0.0, new_params[0][0])
	print  " new params ", cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale, new_params[1]
	set_params3D(e, [cphi, ctheta, cpsi, cs2x, cs2y, cs2z, 0, cscale])
	
	if type(vol) == type(""):
		from utilities import write_headers
		write_headers( vol, [e], [0])
	else:
		return e

def rot_sym(infile, outfile, sym_gp="d4", \
			radius=None, phi=0, theta=0, psi=0, phirange=20, thetarange=20, psirange=20, ftolerance=1.e-4, xtolerance=1.e-4):
	
	from alignment     import find_symm
	from utilities     import drop_image, model_circle, sym_vol
	from fundamentals  import  rot_shift3D
	
	e=EMData()
	e.read_image(infile)
	mask = EMData()
	if radius == None: radius = e.get_xsize()/2.0
	mask  = model_circle(radius, e.get_xsize(), e.get_ysize(), e.get_zsize())
	scale = [phirange, thetarange, psirange]
	res = find_symm(e, mask, sym_gp, phi, theta, psi, scale, ftolerance, xtolerance)
	print  "  RESULTS: ",res

	sym = sym_vol(rot_shift3D(e, res[0][0], res[0][1], res[0][2] ), sym_gp)

	drop_image(sym, outfile)
	
def transform2d(stack_data, stack_data_ali):
# apply 2D alignment parameters stored in the header of the input stack file using gridding interpolation and create an output stack file
	from fundamentals   import rot_shift2D
	from utilities 	    import set_arb_params, set_params2D, get_params2D
	from utilities      import print_begin_msg, print_end_msg, print_msg
	import os
	
	print_begin_msg("transform2d")	
	print_msg("Input stack                 : %s\n"%(stack_data))
	print_msg("Output stack                : %s\n\n"%(stack_data_ali))

	if os.path.exists(stack_data_ali): os.system("rm -f "+stack_data_ali)

	attributes = ['nclass', 'assign']
	t = Transform({"type":"2D"})
	data = EMData()
	nima = EMUtil.get_image_count(stack_data)
	data.read_image(stack_data, 0)
	for im in xrange(nima):
		if im>0:
			data = EMData()
			data.read_image(stack_data, im)
		l = data.get_attr_dict()
		params = []
		for ia in xrange(len(attributes)):
			if(attributes[ia] in l):
				params.append(data.get_attr(attributes[ia]))
			else:
				params.append(0)
		al2d = get_params2D(data)
		# apply params to the image
		temp = rot_shift2D(data, al2d[0], al2d[1], al2d[2], al2d[3])
		temp.set_attr("xform.align2d", t)
		set_arb_params(temp, params, attributes)
		temp.write_image(stack_data_ali, im)
	print_end_msg("transform2d")

def recons3d_n(prj_stack, pid_list, vol_stack, CTF=False, snr=1.0, sign=1, npad=4, sym="c1", listfile = "", group = -1, verbose=0, MPI=False,xysize=-1, zsize = -1):
	if MPI:
		recons3d_n_MPI(prj_stack, pid_list, vol_stack, CTF, snr, 1, npad, sym, listfile, group, verbose,xysize, zsize)
		return

	from reconstruction import recons3d_4nn_ctf, recons3d_4nn
	from utilities import drop_image
	from utilities import print_begin_msg, print_end_msg, print_msg

	print_begin_msg("recons3d_n")
	print_msg("Input stack                 : %s\n"%(prj_stack))
	print_msg("Output volume               : %s\n"%(vol_stack))
	print_msg("Padding factor              : %i\n"%(npad))
	print_msg("CTF correction              : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	print_msg("CTF sign                    : %i\n"%(sign))
	print_msg("Symmetry group              : %s\n\n"%(sym))
	if(listfile):
		from utilities import read_text_file
		pid_list = read_text_file(listfile, 0)
		pid_list = map(int, pid_list)
		print_msg("Reconstruction for images listed in file : %s\n\n"%(listfile))
	elif(group > -1):
		print_msg("Reconstruction for group             : %i\n\n"%(group))
		tmp_list = EMUtil.get_all_attributes(prj_stack, 'group')
		pid_list = []
		for i in xrange(len(tmp_list)):
			if(tmp_list[i] == group):  pid_list.append(i)
		del tmp_list

	if CTF: vol = recons3d_4nn_ctf(prj_stack, pid_list, snr, 1, sym, verbose, npad, xysize=xysize, zsize=zsize)
	else:   vol = recons3d_4nn(prj_stack,  pid_list, sym, npad, xysize=xysize, zsize = zsize)
	if(vol_stack[-3:] == "spi"):
		drop_image(vol, vol_stack, "s")
	else:
		drop_image(vol, vol_stack)
	print_end_msg("recons3d_n")

def recons3d_n_MPI(prj_stack, pid_list, vol_stack, CTF, snr, sign, npad, sym, listfile, group, verbose,xysize, zsize):
	from reconstruction import recons3d_4nn_ctf_MPI, recons3d_4nn_MPI
	from utilities      import get_im, drop_image, bcast_number_to_all
	from utilities      import print_begin_msg, print_end_msg, print_msg
	from string         import replace
	from time           import time
	from mpi 	    import mpi_comm_size, mpi_comm_rank, mpi_bcast, MPI_INT, MPI_COMM_WORLD

	myid  = mpi_comm_rank(MPI_COMM_WORLD)
	nproc = mpi_comm_size(MPI_COMM_WORLD)
	time_start = time()

	if(myid == 0):
		print_begin_msg("recons3d_n_MPI")
		print_msg("Input stack  	       : %s\n"%(prj_stack))
		print_msg("Output volume	       : %s\n"%(vol_stack))
		print_msg("Padding factor	       : %i\n"%(npad))
		print_msg("CTF correction	       : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("CTF sign		       : %i\n"%(sign))
		print_msg("Symmetry group	       : %s\n\n"%(sym))
		if(listfile):
			from utilities import read_text_file
			pid_list = read_text_file(listfile, 0)
			pid_list = map(int, pid_list)
			print_msg("Reconstruction for images listed in file : %s\n\n"%(listfile))
		elif(group > -1):
			print_msg("Reconstruction for group		: %i\n\n"%(group))
			tmp_list = EMUtil.get_all_attributes(prj_stack, 'group')
			pid_list = []
			for i in xrange(len(tmp_list)):
				if(tmp_list[i] == group):  pid_list.append(i)
			del tmp_list
		nima = len(pid_list)
	else:
		nima = 0

	nima = bcast_number_to_all(nima, source_node = 0)

	if(listfile or group > -1):
		if myid != 0:
			pid_list = [-1]*nima
		pid_list = mpi_bcast(pid_list, nima, MPI_INT, 0, MPI_COMM_WORLD)
		pid_list = map(int, pid_list)
	else:
		if(not pid_list):  pid_list = range(nima)

	if verbose==0:
		finfo = None
	else:
		infofile = "progress%04d.txt"%(myid+1)
		finfo = open( infofile, 'w' )

	image_start, image_end = MPI_start_end(nima, nproc, myid)

	prjlist = EMData.read_images(prj_stack, pid_list[image_start:image_end])
	del pid_list

	if CTF: vol = recons3d_4nn_ctf_MPI(myid, prjlist, snr, sign, sym, finfo, npad,xysize, zsize)
	else:	vol = recons3d_4nn_MPI(myid, prjlist, sym, finfo, npad,xysize, zsize)
	if myid == 0 :
		if(vol_stack[-3:] == "spi"):
			drop_image(vol, vol_stack, "s")
		else:
			drop_image(vol, vol_stack)
		if not(finfo is None):
			finfo.write( "result written to " + vol_stack + "\n")
			finfo.write( "Total time: %10.3f\n" % (time()-time_start) )
			finfo.flush()

def recons3d_f(prj_stack, vol_stack, fsc_file, mask=None, CTF=True, snr=1.0, sym="c1", listfile = "", group = -1, npad = 4, verbose=1, MPI=False):
	if MPI:
		recons3d_f_MPI(prj_stack, vol_stack, fsc_file, mask, CTF, snr, sym, listfile, group, npad, verbose)
		return

	nima = EMUtil.get_image_count( prj_stack )

	from reconstruction import recons3d_4nn_ctf, recons3d_4nn
	from statistics     import fsc_mask
	from utilities      import drop_image
	if(listfile):
		from utilities import read_text_file
		pid_list = read_text_file(listfile, 0)
		pid_list = map(int, pid_list)
	elif(group > -1):
			tmp_list = EMUtil.get_all_attributes(prj_stack, 'group')
			pid_list = []
			for i in xrange(len(tmp_list)):
				if(tmp_list[i] == group):  pid_list.append(i)
			del tmp_list
	else:
		pid_list = range(nima)
	if CTF:
		volodd = recons3d_4nn_ctf(prj_stack, [ pid_list[i] for i in xrange(0, len(pid_list), 2) ], snr, 1, sym, verbose, npad)
		voleve = recons3d_4nn_ctf(prj_stack, [ pid_list[i] for i in xrange(1, len(pid_list), 2) ], snr, 1, sym, verbose, npad)
		t = fsc_mask( volodd, voleve, mask, filename=fsc_file)
		del volodd, voleve
		volall = recons3d_4nn_ctf(prj_stack, pid_list,                                          snr, 1, sym, verbose, npad)
	else:
		volodd = recons3d_4nn(prj_stack, [ pid_list[i] for i in xrange(0, len(pid_list), 2) ], sym, npad)
		voleve = recons3d_4nn(prj_stack, [ pid_list[i] for i in xrange(1, len(pid_list), 2) ], sym, npad)
		t = fsc_mask( volodd, voleve, mask, filename=fsc_file)
		del volodd, voleve
		volall = recons3d_4nn(prj_stack, pid_list,                                          sym, npad)
	if(vol_stack[-3:] == "spi"):
		drop_image(volall, vol_stack, "s")
	else:
		drop_image(volall, vol_stack)

def recons3d_f_MPI(prj_stack, vol_stack, fsc_file, mask, CTF=True, snr=1.0, sym="c1", listfile="", group=-1, npad = 4, verbose=1):

	from mpi       import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD, mpi_bcast, MPI_INT
	from utilities import drop_image, bcast_number_to_all
	nproc = mpi_comm_size( MPI_COMM_WORLD )
	myid  = mpi_comm_rank( MPI_COMM_WORLD )
	
	if(myid == 0):
		if(listfile):
			from utilities import read_text_file
			pid_list = read_text_file(listfile, 0)
			pid_list = map(int, pid_list)
			nima = len(pid_list)
		elif(group > -1):
			tmp_list = EMUtil.get_all_attributes(prj_stack, 'group')
			pid_list = []
			for i in xrange(len(tmp_list)):
				if(tmp_list[i] == group):  pid_list.append(i)
			del tmp_list
			nima = len(pid_list)
		else:
			nima = EMUtil.get_image_count(prj_stack)
			pid_list = range(nima)
	else:
		nima = 0

	nima = bcast_number_to_all(nima, source_node = 0)

	if myid != 0:
		pid_list = [-1]*nima
	pid_list = mpi_bcast(pid_list, nima, MPI_INT, 0, MPI_COMM_WORLD)
	pid_list = map(int, pid_list)

	image_start, image_end = MPI_start_end(nima, nproc, myid)

	imgdata = EMData.read_images(prj_stack, pid_list[image_start:image_end])
	del pid_list

	if verbose==0:
		finfo = None
	else:
		infofile = "progress%04d.txt" % (myid)
		finfo = open( infofile, 'w' )

	odd_start = image_start%2
	eve_start = (odd_start+1)%2
	if CTF:
		from reconstruction import rec3D_MPI

		vol,fsc = rec3D_MPI(imgdata, snr, sym, mask, fsc_file, myid, 0, 1.0, odd_start, eve_start, finfo, npad = npad)
	else :
		from reconstruction import rec3D_MPI_noCTF
		vol,fsc = rec3D_MPI_noCTF(imgdata, sym, mask, fsc_file, myid, 0, 1.0, odd_start, eve_start, finfo, npad = npad)
	if myid == 0:
		if(vol_stack[-3:] == "spi"):
			drop_image(vol, vol_stack, "s")
		else:
			drop_image(vol, vol_stack)

def ssnr3d(stack, output_volume = None, ssnr_text_file = None, mask = None, reference_structure = None, ou = -1, rw = 1.0,  npad = 1, CTF = False, sign = 1, sym ="c1", MPI = False, random_angles = 0):
	'''
	        Perform 3D reconstruction using selected particle images, 
	        and calculate spectrum signal noise ratio (SSNR).
	        1. The selection file is supposed to be in SPIDER text format.
	        2. The 3D alignment parameters have been written in headers of the particle images.
	''' 
	if MPI:
		ssnr3d_MPI(stack, output_volume, ssnr_text_file, mask, reference_structure, ou, rw, npad, CTF, sign, sym, random_angles)
		return

	from utilities               import model_circle, get_im
	from filter                  import filt_ctf
	from reconstruction          import recons3d_nn_SSNR, recons3d_4nn, recons3d_4nn_ctf
	from projection              import prep_vol, prgs
	from utilities               import print_begin_msg, print_end_msg, print_msg
	
	print_begin_msg("ssnr3d")

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Output volume               : %s\n"%(output_volume))	
	print_msg("SSNR text file              : %s\n"%(ssnr_text_file))
	print_msg("Outer radius                : %i\n"%(ou))
	print_msg("Ring width                  : %i\n"%(rw))
	print_msg("Padding factor              : %i\n"%(npad))
	print_msg("CTF correction              : %s\n"%(CTF))
	print_msg("CTF sign                    : %i\n"%(sign))
	print_msg("Symmetry group              : %s\n"%(sym))
	
	fring_width = float(rw)
	if mask:
		import  types
		if type(mask) is types.StringType:
			print_msg("Maskfile                    : %s\n"%(mask))
			mask2D=get_im(mask)
		else:
			print_msg("Maskfile                    : user provided in-core mask\n")
			mask2D = mask
	else:
		print_msg("Maskfile                    : None\n")
		mask2D = None

	[ssnr1, vol_ssnr1] = recons3d_nn_SSNR(stack, mask2D, rw, npad, sign, sym, CTF, random_angles)
	vol_ssnr1.write_image(output_volume, 0)
	del vol_ssnr1
	# perform 3D reconstruction
	if reference_structure == None:
		nima = EMUtil.get_image_count(stack)
		if CTF:
			snr = 1.0e20
			vol = recons3d_4nn_ctf(stack, range(nima), snr, sign, sym)
		else:   vol = recons3d_4nn(stack, range(nima), sym)
	else:
		vol = get_im(reference_structure)

	# re-project the reconstructed volume
	nx = vol.get_xsize()
	if int(ou) == -1: radius = nx//2 - 1
	else :            radius = int(ou)
	#
	vol *= model_circle(radius, nx, nx, nx)
	volft, kb = prep_vol(vol)
	del vol
	prjlist = []
	from utilities import get_params_proj
	for i in xrange(nima):
		e = EMData()
		e.read_image(stack, i, True)
		e.set_attr('sign', 1)
		phi, theta, psi, tx, ty = get_params_proj(e)
		proj = prgs(volft, kb, [phi, theta, psi, -tx, -ty])
		if CTF :
			ctf_params = e.get_attr("ctf")			
			proj = filt_ctf(proj, ctf_params)
		prjlist.append(proj)
	del volft
	[ssnr2, vol_ssnr2] = recons3d_nn_SSNR(prjlist, mask2D, rw, npad, sign, sym, CTF, random_angles)
	vol_ssnr2.write_image(output_volume, 1)
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

	print_end_msg('ssnr3d')
	
	'''
	qt = 0.0
	for i in xrange(len(ssnr1)):
		tqt = ssnr1[i][1] - ssnr2[i][1]
		if( tqt<qt ): qt = tqt
	for i in xrange(len(ssnr1)): ssnr1[i][1] -= (ssnr2[i][1] + qt)
	from utilities import dropSpiderDoc19289
	dropSpiderDoc(ssnr_text_file+".doc", ssnr1)
	dropImage(vol_ssnr2, output_volume+"2.spi", "s")
	'''

def ssnr3d_MPI(stack, output_volume = None, ssnr_text_file = None, mask = None, reference_structure = None, ou = -1, rw = 1.0, npad = 1, CTF = False, sign = 1, sym ="c1", random_angles = 0):
	from reconstruction import recons3d_nn_SSNR_MPI, recons3d_4nn_MPI, recons3d_4nn_ctf_MPI
	from utilities      import bcast_EMData_to_all, model_blank, model_circle, get_im
	from projection     import prep_vol, prgs
	from mpi            import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD

	nima = EMUtil.get_image_count(stack)
	nproc = mpi_comm_size(MPI_COMM_WORLD)
	myid  = mpi_comm_rank(MPI_COMM_WORLD)

	image_start, image_end = MPI_start_end(nima, nproc, myid)

	if mask:
		import  types
		if type(mask) is types.StringType: mask2D = get_im(mask)
		else: mask2D = mask
	else:
		mask2D = None

	prjlist = EMData.read_images(stack, range(image_start, image_end))
	if random_angles > 0:
		for prj in prjlist:
			active = prj.get_attr_default('active', 1)
			if active == 1:
				if random_angles == 2:
					from  random import  random
					phi	 = 360.0*random()
					theta	 = 180.0*random()
					psi	 = 360.0*random()
					xform_proj = Transform( {"type":"spider", "phi":phi, "theta":theta, "psi":psi} )
					prj.set_attr("xform.projection", xform_proj)
				elif random_angles == 3:
					from  random import  random
					phi    = 360.0*random()
					theta  = 180.0*random()
					psi    = 360.0*random()
					tx     = 6.0*(random() - 0.5)
					ty     = 6.0*(random() - 0.5)
					xform_proj = Transform( {"type":"spider", "phi":phi, "theta":theta, "psi":psi, "tx":tx, "ty":ty} )
					prj.set_attr("xform.projection", xform_proj)
				elif random_angles  == 1:
					from  random import  random
					old_xform_proj = prj.get_attr( "xform.projection" )
					dict = old_xform_proj.get_rotation( "spider" )
					dict["psi"] = 360.0*random()
					xform_proj = Transform( dict )
					prj.set_attr("xform.projection", xform_proj)
		random_angles = 0
	if myid == 0: [ssnr1, vol_ssnr1] = recons3d_nn_SSNR_MPI(myid, prjlist, mask2D, rw, npad, sign, sym, CTF, random_angles)  
	else:	                           recons3d_nn_SSNR_MPI(myid, prjlist, mask2D, rw, npad, sign, sym, CTF, random_angles)
	if myid == 0:
		vol_ssnr1.write_image( output_volume, 0)
		del vol_ssnr1

	nx  = prjlist[0].get_xsize()
	if ou == -1: radius = int(nx/2) - 1
	else:        radius = int(ou)
	if(reference_structure == None):
		vol = model_blank(nx, nx, nx)
		if CTF:
			snr = 1.0e20
			if myid == 0 : vol = recons3d_4nn_ctf_MPI(myid, prjlist, snr, sign, sym)
			else :  	     recons3d_4nn_ctf_MPI(myid, prjlist, snr, sign, sym)
		else:
			if myid == 0 : vol = recons3d_4nn_MPI(myid, prjlist, sym)
			else:		     recons3d_4nn_MPI(myid, prjlist, sym)
	else:
		if myid == 0: vol = get_im(reference_structure)

	bcast_EMData_to_all(vol, myid, 0)
	re_prjlist = []
	#vol *= model_circle(radius, nx, nx, nx)
	volft, kb = prep_vol(vol)
	del vol
	from utilities import get_params_proj
	if CTF: from filter import filt_ctf
	for prj in prjlist:
		phi, theta, psi, tx, ty = get_params_proj(prj)
		proj = prgs(volft, kb, [phi, theta, psi, -tx, -ty])
		if CTF:
			ctf_params = prj.get_attr("ctf")			
			proj = filt_ctf(proj, ctf_params)
			proj.set_attr('sign', 1)
		re_prjlist.append(proj)
	del volft, prjlist
	if myid == 0: [ssnr2, vol_ssnr2] = recons3d_nn_SSNR_MPI(myid, re_prjlist, mask2D, rw, npad, sign, sym, CTF, random_angles)
	else:                              recons3d_nn_SSNR_MPI(myid, re_prjlist, mask2D, rw, npad, sign, sym, CTF, random_angles)
	if myid == 0:
		vol_ssnr2.write_image( output_volume, 1)
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
		"""
		qt = 0.0
		for i in xrange(len(ssnr2)):
			tqt = ssnr1[i][1] - ssnr2[i][1]
			if( tqt<qt ): qt = tqt
		for i in xrange(len(ssnr1)): ssnr1[i][1] -= (ssnr2[i][1] + qt)
		
		dropSpiderDoc(ssnr_text_file+".doc", ssnr1)
		vol_ssnr2, output_volume+"2.spi", "s")
		"""

def pca(input_stacks, output_stack, subavg, mask_radius, sdir, nvec, shuffle, genbuf, maskfile="", MPI=False, verbose=False):
	from utilities import get_image, get_im, model_circle, model_blank
	from statistics import pcanalyzer

	if len(input_stacks)==0:
		print "Error: no input file."
		return

	if mask_radius > 0 and maskfile !="":
		print "Error: mask radius and mask file cannot be used at the same time"
		return

	if mask_radius >0:

		if(verbose):
			print "Using spherical mask, rad=", mask_radius

		if maskfile!="":
			print "Error: mask radius and mask file cannot be given at the same time"
			return

		data = get_im( input_stacks[0] )
		mask = model_circle(mask_radius, data.get_xsize(), data.get_ysize(), data.get_zsize())

	elif(maskfile!="") :
		if(verbose):
			print "Using mask: ", maskfile
		mask = get_image( maskfile )
	else:
		data = EMData()
		data.read_image( input_stacks[0], 0, True)
		mask = model_blank(data.get_xsize(), data.get_ysize(), data.get_zsize(), bckg=1.0)

	pca = pcanalyzer(mask, sdir, nvec, MPI)


	if subavg != "":
		if(verbose):
			print "Subtracting ", subavg, " from each image"
		avg = get_image( subavg )
		pca.setavg( avg )

	files = file_set( input_stacks )
	if MPI:
		from mpi import mpi_comm_rank, mpi_comm_size, MPI_COMM_WORLD
		myid = mpi_comm_rank( MPI_COMM_WORLD )
		ncpu = mpi_comm_size( MPI_COMM_WORLD )
	else:
		myid = 0
		ncpu = 1

	if genbuf:
		if shuffle:
			print "Error: shuffle works only with usebuf"
			return

		bgn,end = MPI_start_end( files.nimg(), ncpu, myid )
		for i in xrange(bgn,end):
			fname, imgid = files.get( i )
			
			data = get_im( fname, imgid)
			pca.insert( data )
			if(verbose):
				 print "Inserting image %s, %4d" % (fname, imgid)
	else:
		pca.usebuf( )
		print myid, "using existing buff, nimg: ", pca.nimg
		if shuffle:
			pca.shuffle()

	vecs = pca.analyze()
	if myid==0:
		for i in xrange( len(vecs) ):
			vecs[i].write_image( output_stack, i)


def prepare_2d_forPCA(input_stack, output_stack, average, avg = False, CTF = False):
	"""
		Prepare 2D images for PCA
		Average of all images is calculated using header alignment information, 
		  subtracted from each image and the difference is written to the output_stack
		If CTF, average is calculated as
		   Av = sum(CTF_k*Im_k)/sum(CTF_k^2
		and the difference as
		   Im_k - CTF_k*Av
		If avg = True, average outside of a circle r = nx//2-1 is subtracted from each image
	"""
	from utilities    import model_blank, model_circle, get_arb_params, set_arb_params
	from fundamentals import rot_shift2D
	pali = ["alpha", "sx", "sy", "mirror"]
	n = EMUtil.get_image_count(input_stack)
	ima = EMData()
	ima.read_image(input_stack, 0)
	nx = ima.get_xsize()
	ny = ima.get_xsize()
	
	if avg:  mask = model_circle( nx//2-1, nx, ny)
	if  CTF:
		if(ima.get_attr_default('ctf_applied', 2) > 0):
			ERROR("data cannot be ctf-applied","prepare_2d_forPCA",1)
		from fundamentals import fft
		from morphology   import ctf_img
		from filter 	  import filt_ctf
		from utilities    import pad
		parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]

		nx2 = 2*nx
		ny2 = 2*ny
		ave       = EMData(nx2, ny2, 1, False)
		ctf_2_sum = EMData(nx2, ny2, 1, False)

		for i in xrange(n):
			ima = EMData()
			ima.read_image(input_stack, i)
			ctf_params = ima.get_arb_params("ctf")
			ali = get_arb_params(ima, pali)
			ima    = rot_shift2D(ima, ali[0], ali[1], ali[2], ali[3])
			if avg:
				st = Util.infomask(ima, mask, False)
				ima -= st[0]
			oc = filt_ctf(fft(pad(ima, nx2, ny2, background = 0.0)), ctf_params)
			Util.add_img(ave, oc)
			Util.add_img2(ctf_2_sum, ctf_img(nx2, ctf_params, ny = ny2, nz = 1))
		Util.div_filter(ave, ctf_2_sum)
		for i in xrange(n):
			ima = EMData()
			ima.read_image(input_stack, i)
			ctf_params = ima.get_attr( "ctf" )
			ali = get_arb_params(ima, pali)
			ima    = rot_shift2D(ima, ali[0], ali[1], ali[2], ali[3])
			if avg:
				st = Util.infomask(ima, mask, False)
				ima -= st[0]
			oc = filt_ctf(ave, ctf_params, pad= True)
			Util.sub_img(ima, Util.window(fft(oc),nx,ny,1,0,0,0))
			set_arb_params(ima, [0.0,0.0,0.0,0], pali)
			ima.write_image(output_stack, i)
		Util.window(fft(ave),nx,ny,1,0,0,0).write_image(average)
	else:
		ave  = model_blank( nx, ny)
		for i in xrange(n):
			ima = EMData()
			ima.read_image(input_stack, i)
			ali = get_arb_params(ima, pali)
			ima    = rot_shift2D(ima, ali[0], ali[1], ali[2], ali[3])
			if avg:
				st = Util.infomask(ima, mask, False)
				ima -= st[0]
			Util.add_img(ave, ima)
		ave /= n
		ave.write_image(average)
		for i in xrange(n):
			ima = EMData()
			ima.read_image(input_stack, i)
			ali = get_arb_params(ima, pali)
			ima    = rot_shift2D(ima, ali[0], ali[1], ali[2], ali[3])
			if avg:
				st = Util.infomask(ima, mask, False)
				ima -= st[0]
			Util.sub_img(ima, ave)
			set_arb_params(ima, [0.0,0.0,0.0,0], pali)
			ima.write_image(output_stack, i)

def varimax(input_stack, imglist, output_stack, maskfile, mask_radius, verbose ) :
	from utilities import get_im, model_circle
	from EMAN2     import Analyzers

	data = get_im( input_stack )

	if maskfile:
		import types
		if type(maskfile) is types.StringType: mask = get_im(maskfile)
		else:                                  mask = maskfile
	else:
		if(mask_radius < 1):  mask_radius = data.get_xsize()//2-2
		mask = model_circle( mask_radius, data.get_xsize(), data.get_ysize(), data.get_zsize() )

	ana = Analyzers.get( "varimax", {"mask":mask} )
	sumeig =0.0
	#from utilities import info
	#from math import sqrt
	for i in imglist:
		data = get_im( input_stack, i)
		eigval = data.get_attr_default('eigval', 1.0)
		sumeig += eigval
		#Util.mul_scalar(data, sqrt(eigval))
		#info(data)
		ana.insert_image( data )
		#print "Inserting image %4d" % i
	del data
	vecs = ana.analyze()

	sumeig /= len(vecs)
	#print  sumeig
	for iout in xrange(len(vecs)):
		#info(vecs[iout])
		vecs[iout].set_attr('eigval', sumeig)
		vecs[iout].write_image( output_stack, iout)

def bootstrap_genbuf(prj_stack, buf_prefix, npad, verbose, CTF=False):
	from EMAN2 import newfile_store
	import os
	size = 1
	myid = 0
	print  
	if os.path.exists( buf_prefix + ".bin" ):
		ERROR('Output file exists, please change the name and restart the program', "bootstrap_genbuf", 1)

	if(verbose == 1):  finfo=open( os.path.join(outdir, "progress%04d.txt" % myid), "w" )
	else:              finfo = None

	store = newfile_store(buf_prefix, npad, CTF)

	nimage = EMUtil.get_image_count( prj_stack )
	for i in xrange(nimage):
		proj = EMData()
		proj.read_image( prj_stack, i )
		store.add_image( proj, proj.get_attr("xform.projection") )

		if( verbose == 1 and ((i+1) % 100 == 0  or i==nimage-1)) :
			finfo.write( "projection %6d buffered\n" % (i+1) )
			finfo.flush()
 
def bootstrap_run(prj_stack, media, outdir, nvol, CTF, snr, sym, verbose, MPI=False):

	import string
	from mpi import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD, mpi_init, mpi_barrier
	import os

	if MPI:
		size = mpi_comm_size(MPI_COMM_WORLD)
		myid = mpi_comm_rank(MPI_COMM_WORLD)
	else:
		size = 1
		myid = 0

	if myid==0:
		if os.path.exists(outdir): ERROR('Output directory exists, please change the name and restart the program', "bootstrap_run", 1,myid)
		os.system( "mkdir " + outdir )
	if MPI:
		mpi_barrier( MPI_COMM_WORLD )	


	myvolume_file = "%s/bsvol%04d.hdf" % (outdir, myid)
	if verbose != 0 :
		mystatus_file = "%s/status%04d.inf" % (outdir, myid)
		mystatus = open( mystatus_file, 'w' )
	else:
		mystatus = None

	mynvol = nvol/size

	if myid==(size-1) : mynvol = mynvol + (nvol%size)


	nproj = EMUtil.get_image_count(prj_stack)
	if verbose != 0 : mystatus.write( "# of projs: %d\n" % nproj )

	npad = 4
	sign = 1
	list_proj = range(nproj)
	from reconstruction import bootstrap_nn
	bootstrap_nn( prj_stack, myvolume_file, list_proj, mynvol, media, npad, sym, mystatus, CTF, snr, sign)
	
def params_2D_to_3D(stack):
	from utilities import params_2D_3D, print_begin_msg, print_end_msg, print_msg, get_params2D, set_params_proj, write_header
	
	#print_begin_msg("params_2D_to_3D")
	#print_msg("Input stack                 : %s\n\n"%(stack))

	nima = EMUtil.get_image_count(stack)
	ima = EMData()
	for im in xrange(nima):
		ima.read_image(stack, im, True)
		p = get_params2D(ima)
		p = params_2D_3D(p[0], p[1], p[2], int(p[3]))
		set_params_proj(ima, p)
		write_header(stack, ima, im)
	#print_end_msg("params_2D_to_3D")
	
def params_3D_to_2D(stack):
	from utilities import params_3D_2D, print_begin_msg, print_end_msg, print_msg, set_params2D, write_header
	
	#print_begin_msg("params_3D_to_2D")
	#print_msg("Input stack                 : %s\n\n"%(stack))

	nima = EMUtil.get_image_count(stack)
	ima = EMData()
	for im in xrange(nima):
		ima.read_image(stack, im, True)
		from utilities import set_params_proj, get_params_proj
		phi,theta,psi,s2x,s2y = get_params_proj( ima )
		alpha, sx, sy, mirror = params_3D_2D(phi, theta, psi, s2x, s2y)
		set_params2D(ima, [alpha, sx, sy, mirror, 1.0])
		write_header(stack, ima, im)
	#print_end_msg("params_3D_to_2D")


# application find structure
def cml_find_structure_main(stack, out_dir, ir, ou, delta, dpsi, lf, hf, rand_seed, maxit, given = False, first_zero = False, flag_weights = False, debug = False, trials = 1):
	from projection import cml_open_proj, cml_init_global_var, cml_head_log, cml_disc, cml_export_txtagls
	from projection import cml_find_structure, cml_export_struc, cml_end_log
	from utilities  import print_begin_msg, print_msg, print_end_msg, start_time, running_time
	from copy       import deepcopy
	from random     import seed, random
	import time, sys, os

	# logfile
	t_start = start_time()
	print_begin_msg('find_struct')
	print_msg('\n')

	out_dir = out_dir.rstrip('/')
	if os.path.exists(out_dir): ERROR('Output directory exists, please change the name and restart the program', "cml_find_structure_main", 1)
	os.mkdir(out_dir)

	if rand_seed > 0: seed(rand_seed)
	else:             seed()

	# Open and transform projections
	Prj, Ori = cml_open_proj(stack, ir, ou, lf, hf, dpsi)
	# Init the global vars
	cml_init_global_var(dpsi, delta, len(Prj), debug)
	# Update logfile
	cml_head_log(stack, out_dir, delta, ir, ou, lf, hf, rand_seed, maxit, given, flag_weights, trials, 1)

	ibest    = -1
	bestdisc = 1.0e20
	MEM      = []
	for itrial in xrange(trials):
		# if not angles given select randomly orientation for each projection
		if not given:
			j = 0
			for n in xrange(len(Prj)):
				if first_zero and n == 0:
					Ori[j]   = 0.0
					Ori[j+1] = 0.0
					Ori[j+2] = 0.0
				else:
					Ori[j]   = random() * 360  # phi
					Ori[j+1] = random() * 180  # theta
					Ori[j+2] = random() * 360  # psi
				j += 4

		# prepare rotation matrix
		Rot = Util.cml_init_rot(Ori)
		# Compute the first disc
		disc_init = cml_disc(Prj, Ori, Rot)
		# Update progress file
		cml_export_txtagls(out_dir, 'angles_%03i' % itrial, Ori, disc_init, 'Init')
		# Find structure
		Ori, disc, ite = cml_find_structure(Prj, Ori, Rot, out_dir, 'angles_%03i' % itrial, maxit, first_zero, flag_weights)
		print_msg('Trial %03i\tdiscrepancy init: %10.7f\tnb ite: %i\tdiscrepancy end: %10.7f\n' % (itrial, disc_init, ite + 1, disc))
		if disc < bestdisc:
			bestdisc = disc
			ibest    = itrial
			MEM      = deepcopy(Ori)

		# Export structure
		cml_export_struc(stack, out_dir, itrial, Ori)

	print_msg('\n Selected trial #%03i with disc %10.7f\n' % (ibest, bestdisc))
	os.system('cp %s/structure_%03i.hdf %s/structure.hdf' % (out_dir, ibest, out_dir))
	cml_end_log(MEM)
	running_time(t_start)
	print_end_msg('find_struct')

# application find structure
def cml_find_structure_MPI2(stack, out_dir, ir, ou, delta, dpsi, lf, hf, rand_seed, maxit, given = False, first_zero = False, flag_weights = False, debug = False, trials = 1):
	from projection import cml_open_proj, cml_init_global_var, cml_head_log, cml_disc, cml_export_txtagls
	from projection import cml_find_structure2, cml_export_struc, cml_end_log
	from utilities  import print_begin_msg, print_msg, print_end_msg, start_time, running_time
	from copy       import deepcopy
	from random     import seed, random
	import time, sys, os
	from mpi 	  import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_gatherv
	from mpi 	  import MPI_SUM, MPI_FLOAT, MPI_INT

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	# logfile
	if myid == main_node:
		t_start = start_time()
		print_begin_msg('find_struct')
		print_msg('\n')

		out_dir = out_dir.rstrip('/')
		if os.path.exists(out_dir): ERROR('Output directory exists, please change the name and restart the program', "cml_find_structure_main", 1)
		os.mkdir(out_dir)

	if rand_seed > 0: seed(rand_seed)
	else:             seed()

	# Open and transform projections
	Prj, Ori = cml_open_proj(stack, ir, ou, lf, hf, dpsi)
	# Init the global vars
	cml_init_global_var(dpsi, delta, len(Prj), debug)
	# Update logfile
	if myid == main_node:
		cml_head_log(stack, out_dir, delta, ir, ou, lf, hf, rand_seed, maxit, given, flag_weights, trials, number_of_proc)

	ibest    = -1
	bestdisc = 1.0e20
	MEM      = []
	for itrial in xrange(trials):
		# if not angles given select randomly orientation for each projection
		if not given:
			j = 0
			for n in xrange(len(Prj)):
				if first_zero and n == 0:
					Ori[j]   = 0.0
					Ori[j+1] = 0.0
					Ori[j+2] = 0.0
				else:
					Ori[j]   = random() * 360  # phi
					Ori[j+1] = random() * 180  # theta
					Ori[j+2] = random() * 360  # psi
				j += 4

		# prepare rotation matrix
		Rot = Util.cml_init_rot(Ori)
		# Compute the first disc
		disc_init = cml_disc(Prj, Ori, Rot)
		# Update progress file
		if myid == main_node:
			cml_export_txtagls(out_dir, 'angles_%03i' % itrial, Ori, disc_init, 'Init')
		# Find structure
		Ori, disc, ite = cml_find_structure2(Prj, Ori, Rot, out_dir, 'angles_%03i' % itrial, maxit, first_zero, flag_weights, myid, main_node, number_of_proc)
		if myid == main_node:
			print_msg('Trial %03i\tdiscrepancy init: %10.7f\tnb ite: %i\tdiscrepancy end: %10.7f\n' % (itrial, disc_init, ite + 1, disc))
		if disc < bestdisc:
			bestdisc = disc
			ibest    = itrial
			MEM      = deepcopy(Ori)

		# Export structure
		if myid == main_node:
			cml_export_struc(stack, out_dir, itrial, Ori)

	if myid == main_node:
		print_msg('\n Selected trial #%03i with disc %10.7f\n' % (ibest, bestdisc))
		os.system('cp %s/structure_%03i.hdf %s/structure.hdf' % (out_dir, ibest, out_dir))
		cml_end_log(MEM)
		running_time(t_start)
		print_end_msg('find_struct')

# application find structure
def cml_find_structure_MPI(stack, out_dir, ir, ou, delta, dpsi, lf, hf, rand_seed, maxit, given = False, first_zero = False, flag_weights = False, debug = False, trials = 10):
	from projection import cml_open_proj, cml_init_global_var, cml_head_log, cml_disc, cml_export_txtagls
	from projection import cml_find_structure, cml_export_struc, cml_end_log, cml_init_rnd, cml2_ori_collinearity
	from utilities  import print_begin_msg, print_msg, print_end_msg, start_time, running_time
	from random     import seed, random
	from mpi        import mpi_init, mpi_comm_size, mpi_comm_rank, mpi_bcast
	from mpi        import mpi_barrier, MPI_COMM_WORLD, mpi_reduce, MPI_FLOAT, MPI_INT, MPI_SUM
	import time, sys, os

	# init
	sys.argv  = mpi_init(len(sys.argv),sys.argv)
	ncpu      = mpi_comm_size(MPI_COMM_WORLD)
	myid      = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	mpi_barrier(MPI_COMM_WORLD)

	if os.path.exists(out_dir): ERROR('Output directory exists, please change the name and restart the program', "cml_find_structure_MPI ", 1,myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		t_start = start_time()
		print_begin_msg('find_struct')
		os.mkdir(out_dir)


	flag = 0
	if myid == main_node:
		if ncpu > trials:
			print '** WARNING **'
			print 'Find structure MPI: number of trials must be superior than the number of processors.'
			flag = 1
	flag = mpi_bcast(flag, 1, MPI_INT, 0, MPI_COMM_WORLD)
	flag = int(flag[0])
	if flag != 0: sys.exit()

	N_start, N_stop = MPI_start_end(trials, ncpu, myid)
	lrnd    = cml_init_rnd(trials, rand_seed)
	out_dir = out_dir.rstrip('/')

	# Open and transform projections
	Prj, Ori = cml_open_proj(stack, ir, ou, lf, hf, dpsi)

	# Init the global vars
	cml_init_global_var(dpsi, delta, len(Prj), debug)

	# Update logfile
	if myid == main_node: cml_head_log(stack, out_dir, delta, ir, ou, lf, hf, rand_seed, maxit, given, flag_weights, trials, ncpu)

	disc_init = [0.0] * trials
	disc_end  = [0.0] * trials
	ite       = [0]   * trials
	coll      = [0.0] * trials
	for itrial in xrange(N_start, N_stop):

		# if not angles given select randomly orientation for each projection
		if not given:
			seed(lrnd[itrial])
			j = 0
			for n in xrange(len(Prj)):
				if first_zero and n == 0:
					Ori[j]   = 0.0
					Ori[j+1] = 0.0
					Ori[j+2] = 0.0
				else:
					Ori[j]   = random() * 360  # phi
					Ori[j+1] = random() * 180  # theta
					Ori[j+2] = random() * 360  # psi
				j += 4

		# prepare rotation matrix
		Rot = Util.cml_init_rot(Ori)
		# Compute the first disc
		disc_init[itrial] = cml_disc(Prj, Ori, Rot)
		# Update progress file
		cml_export_txtagls(out_dir, 'angles_%03i' % itrial, Ori, disc_init[itrial], 'Init')
		# Find structure
		Ori, disc_end[itrial], ite[itrial] = cml_find_structure(Prj, Ori, Rot, out_dir, 'angles_%03i' % itrial, maxit, first_zero, flag_weights)
		# Export structure
		cml_export_struc(stack, out_dir, itrial, Ori)

	#from development import cml2_ori_collinearity
	coll[itrial] = cml2_ori_collinearity(Ori)

	mpi_barrier(MPI_COMM_WORLD)
	disc_init = mpi_reduce(disc_init, trials, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
	disc_init = disc_init.tolist()
	disc_end  = mpi_reduce(disc_end, trials, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
	disc_end  = disc_end.tolist()
	ite       = mpi_reduce(ite, trials, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
	ite       = ite.tolist()
	coll      = mpi_reduce(coll, trials, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
	coll      = coll.tolist()

	if myid == main_node:
		score = [0.0] * trials
		print_msg('\n')
		for i in xrange(trials):
			score[i] = disc_end[i] * (1 - coll[i])
			print_msg('Trial  %03i\trnd %10i discrepnacy init: %10.7f\tnb ite: %i\tdiscrepancy end: %10.7f\tcollinearity: %f\tscore: %f\n' % (i, lrnd[i], disc_init[i], ite[i] + 1, disc_end[i], coll[i], score[i]))
			
		ibest = disc_end.index(min(disc_end))
		#ibest = score.index(min(score))
		print_msg('\n Selected trial #%03i with discrepancy %10.7f\n' % (ibest, disc_end[ibest]))
		os.system('cp %s/structure_%03i.hdf %s/structure.hdf' % (out_dir, ibest, out_dir))

		running_time(t_start)
		print_end_msg('find_struct')

def extract_value( s ):
	from string import atoi, atof

	try:
		i = atoi( s )
		return i
	except:
		pass

	try:
		f = atof( s )
		return f
	except:
		pass
	
	return s 

def header(stack, params, zero=False, one=False, randomize=False, rand_alpha=False, fimport=None, 
	   fexport=None, fprint=False, backup=False, suffix='_backup', restore=False, delete=False):
	from string    import split
	from utilities import write_header, file_type,generate_ctf
	from random    import random, randint
	from utilities import set_params2D, get_params2D, set_params3D, get_params3D, set_params_proj, get_params_proj, set_ctf, get_ctf
	from EMAN2 import Vec2f


	op = zero+one+randomize+rand_alpha+(fimport!=None)+(fexport!=None)+fprint+backup+restore+delete
	if op == 0:
		print "Error: no operation selected!"
		return
	elif op > 1:
		print "Error: more than one operation at the same time!"
		return


	params = split(params)

	if fimport != None: fimp = open(fimport, 'r')
	if fexport != None: fexp = open(fexport, 'w')

	nimage = EMUtil.get_image_count(stack)
	ext = file_type(stack)
	if ext == "bdb":
		from EMAN2db import db_open_dict
		DB = db_open_dict(stack)
	for i in xrange(nimage):
		if fimport != None:
			line = fimp.readline()
			if len(line)==0 :
				print "Error: file " + fimport + " has only " + str(i) + " lines, while there are " + str(nimage) + " images in the file."
				return
			parmvalues = split(line)
			il=0
			for p in params:
				if p[:13] == "xform.align2d":
					if len(parmvalues) < il+3:
						print "Not enough parameters!"
						return
					alpha = extract_value(parmvalues[il])
					sx = extract_value(parmvalues[il+1])
					sy = extract_value(parmvalues[il+2])
					if len(parmvalues) > (il+3):
						mirror = int(extract_value(parmvalues[il+3]))
					else:
						mirror = 0
					if len(parmvalues) > (il+4):
						scale = extract_value(parmvalues[il+4])
					else:
						scale = 1.0
					#set_params2D(img, [alpha, sx, sy, mirror, scale], params[0])
					t = Transform({"type":"2D","alpha":alpha,"tx":sx,"ty":sy,"mirror":mirror,"scale":scale})
					if ext == "bdb":
						DB.set_attr(i, "xform.align2d", t)
					elif ext == "hdf":					
						EMUtil.write_hdf_attribute(stack, "xform.align2d", t, i)
					il+=5
			
				elif p[:16] == "xform.projection":
					if len(parmvalues) < il+5:
						print "Not enough parameters!"
						return
					phi = extract_value(parmvalues[il])
					theta = extract_value(parmvalues[il+1])
					psi = extract_value(parmvalues[il+2])
					s2x = extract_value(parmvalues[il+3])
					s2y = extract_value(parmvalues[il+4])
					#set_params_proj(img, [phi, theta, psi, s2x, s2y], params[0])
					t = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
					t.set_trans(Vec2f(-s2x, -s2y))
					if ext == "bdb":
						DB.set_attr(i, "xform.projection", t)
					elif ext == "hdf":
						EMUtil.write_hdf_attribute(stack, "xform.projection", t, i)
					il+=5				
				elif p[:13] == "xform.align3d":
					if len(parmvalues) < il+8:
						print "Not enough parameters!"
						return
					phi = extract_value(parmvalues[il])
					theta = extract_value(parmvalues[il+1])
					psi = extract_value(parmvalues[il+2])
					s3x = extract_value(parmvalues[il+3])
					s3y = extract_value(parmvalues[il+4])
					s3z = extract_value(parmvalues[il+5])
					mirror = int(extract_value(parmvalues[il+6]))
					scale = extract_value(parmvalues[il+7])
					#set_params3D(img, [phi, theta, psi, s3x, s3y, s3z, mirror, scale], params[0])
					t = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi,"tx":s3x,"ty":s3y,"tz":s3z,"mirror":mirror,"scale":scale})
					if ext == "bdb":
						DB.set_attr(i, "xform.align3d", t)
					elif ext == "hdf":
						EMUtil.write_hdf_attribute(stack, "xform.align3d", t, i)
					il+=8	
				elif p == "ctf":
					if len(parmvalues) < il+6:
						print "Not enough parameters!"
						return			
					defocus = extract_value(parmvalues[il])
					cs = extract_value(parmvalues[il+1])
					voltage = extract_value(parmvalues[il+2])
					apix = extract_value(parmvalues[il+3])
					bfactor = extract_value(parmvalues[il+4])
					ampcont = extract_value(parmvalues[il+5])
					dfdiff = extract_value(parmvalues[il+6])
					dfang = extract_value(parmvalues[il+7])
					#set_ctf(img, [defocus, cs, voltage, apix, bfactor, ampcont])
					ctf = generate_ctf([defocus, cs, voltage, apix, bfactor, ampcont, dfdiff, dfang]) 
					if ext == "bdb":
						DB.set_attr(i, "ctf", ctf)
					elif ext == "hdf":
						EMUtil.write_hdf_attribute(stack, "ctf", ctf, i)
					il+=6	
				else:
					#if len(params)!=len(parmvalues):
						#print "Error: %d params need to be set, while %d values are provided in line %d of file." % ( len(params), len(parmvalues), i )
						#return
					if ext == "bdb":
						DB.set_attr(i, p, extract_value(parmvalues[il]))
					elif ext == "hdf":
						EMUtil.write_hdf_attribute(stack, p, extract_value(parmvalues[il]), i)
					il+=1
			
		else:
			for p in params:
			
				if zero:
					if p[:13] == "xform.align2d":
						#set_params2D(img, [0.0, 0.0, 0.0, 0, 1.0], p)
						t = Transform({"type":"2D","alpha":0.0,"tx":0.0,"ty":0.0,"mirror":0,"scale":1.0})
						if ext == "bdb":
							DB.set_attr(i, "xform.align2d", t)
						elif ext == "hdf":
							EMUtil.write_hdf_attribute(stack, "xform.align2d", t, i)	
					elif p[:16] == "xform.projection":
						#set_params_proj(img, [0.0, 0.0, 0.0, 0.0, 0.0], p)
						t = Transform({"type":"spider","phi":0.0,"theta":0.0,"psi":0.0})
						t.set_trans(Vec2f(0.0, 0.0))
						if ext == "bdb":
							DB.set_attr(i, "xform.projection", t)
						elif ext == "hdf":
							EMUtil.write_hdf_attribute(stack, "xform.projection", t, i)	
					elif p[:13] == "xform.align3d":
						#set_params3D(img, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 1.0], p)
						t = Transform({"type":"spider","phi":0.0,"theta":0.0,"psi":0.0,"tx":0.0,"ty":0.0,"tz":0.0,"mirror":0,"scale":1.0})
						if ext == "bdb":
							DB.set_attr(i, "xform.align3d", t)
						elif ext == "hdf":
							EMUtil.write_hdf_attribute(stack, "xform.align3d", t, i)
					elif p == "ctf":
						print "Invalid operation!"
						return
					else:
						#img.set_attr(p, 0)
						if ext == "bdb":
							DB.set_attr(i, p, 0)
						elif ext == "hdf":
							EMUtil.write_hdf_attribute(stack,p,0,i)
				elif one:
					if p[:6] == "xform." or p == "ctf":
						print "Invalid operation!"
						return
					else:
						#img.set_attr(p, 1)
						if ext == "bdb":
							DB.set_attr(i, p, 1)
						elif ext == "hdf":
							EMUtil.write_hdf_attribute(stack, p, 1, i)
				elif randomize:
					if p[:13] == "xform.align2d":
						alpha = random()*360.0
						sx = random()*2.0-1.0
						sy = random()*2.0-1.0
						mirror = randint(0, 1)
						scale = 1.0
						#set_params2D(img, [alpha, sx, sy, mirror, scale], p)
						t = Transform({"type":"2D","alpha":alpha,"tx":sx,"ty":sy,"mirror":mirror,"scale":scale})
						if ext == "bdb":
							DB.set_attr(i, "xform.align2d", t)
						elif ext == "hdf":
							EMUtil.write_hdf_attribute(stack, "xform.align2d", t, i)
					elif p[:16] == "xform.projection":
						phi = random()*360.0
						theta = random()*180.0
						psi = random()*360.0
						s2x = random()*4.0-2.0
						s2y = random()*4.0-2.0
						#set_params_proj(img, [phi, theta, psi, s2x, s2y], p)
						t = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
						t.set_trans(Vec2f(-s2x, -s2y))
						if ext == "bdb":
							DB.set_attr(i, "xform.projection", t)
						elif ext == "hdf":
							EMUtil.write_hdf_attribute(stack, "xform.projection", t, i)
					elif p[:13] == "xform.align3d":
						phi = random()*360.0
						theta = random()*180.0
						psi = random()*360.0
						s3x = random()*4.0-2.0
						s3y = random()*4.0-2.0
						s3z = random()*4.0-2.0
						mirror = randint(0, 1)
						scale = 1.0
						#set_params3D(img, [phi, theta, psi, s3x, s3y, s3z, mirror, scale], p)	
						t = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi,"tx":s3x,"ty":s3y,"tz":s3z,"mirror":mirror,"scale":scale})
						if ext == "bdb":
							DB.set_attr(i, "xform.align3d", t)
						elif ext == "hdf":
							EMUtil.write_hdf_attribute(stack, "xform.align3d", t, i)					
					else:
						print "Invalid operation!"
						return						
				elif rand_alpha:
					if p[:13] == "xform.align2d":
						alpha = random()*360.0
						sx = 0.0
						sy = 0.0
						mirror = randint(0, 1)
						scale = 1.0
						#set_params2D(img, [alpha, sx, sy, mirror, scale], p)
						t = Transform({"type":"2D","alpha":alpha,"tx":sx,"ty":sy,"mirror":mirror,"scale":scale})
						if ext == "bdb":
							DB.set_attr(i, "xform.align2d",t)
						elif ext == "hdf":
							EMUtil.write_hdf_attribute(stack, "xform.align2d", t, i)	
					elif p[:16] == "xform.projection":
						phi = random()*360.0
						theta = random()*180.0
						psi = random()*360.0
						s2x = 0.0
						s2y = 0.0
						#set_params_proj(img, [phi, theta, psi, s2x, s2y], p)
					
						t = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
						t.set_trans(Vec2f(-s2x, -s2y))
						if ext == "bdb":
							DB.set_attr(i, "xform.projection",t)
						elif ext == "hdf":
							EMUtil.write_hdf_attribute(stack, "xform.projection", t, i)
					elif p[:13] == "xform.align3d":
						phi = random()*360.0
						theta = random()*180.0
						psi = random()*360.0
						s3x = 0.0
						s3y = 0.0
						s3z = 0.0
						mirror = randint(0, 1)
						scale = 1.0
						#set_params3D(img, [phi, theta, psi, s3x, s3y, s3z, mirror, scale], p)	
						t = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi,"tx":s3x,"ty":s3y,"tz":s3z,"mirror":mirror,"scale":scale})
						if ext == "bdb":
							DB.set_attr(i, "xform.align3d",t)
						elif ext == "hdf":
							EMUtil.write_hdf_attribute(stack, "xform.align3d", t, i)					
					else:
						print "Invalid operation!"
						return	
											
				elif fexport != None:
					if p[:13] == "xform.align2d":
						#alpha, sx, sy, mirror, scale = get_params2D(img, p)
						if ext == "bdb":
							t = DB.get_attr(i, "xform.align2d")
							d = t.get_params("2D")
							fexp.write("%15.5f %15.5f %15.5f %10d %10.3f "%(d["alpha"],d["tx"],d["ty"],d["mirror"],d["scale"]))
														
						elif ext == "hdf":
							t = EMUtil.read_hdf_attribute(stack, "xform.align2d",i)
							d = t.get_params("2D")
							fexp.write("%15.5f %15.5f %15.5f %10d %10.3f "%(d["alpha"],d["tx"],d["ty"],d["mirror"],d["scale"]))
		
					elif p[:16] == "xform.projection":
						#phi, theta, psi, s2x, s2y = get_params_proj(img, p)
						if ext == "bdb":
							t = DB.get_attr(i,"xform.projection")
							d = t.get_params("spider")
							fexp.write("%15.5f %15.5f %15.5f %15.5f %15.5f "%(d["phi"],d["theta"],d["psi"],-d["tx"],-d["ty"]))
							
						elif ext == "hdf":
							t = EMUtil.read_hdf_attribute(stack, "xform.projection",i)
							d = t.get_params("spider")
							fexp.write("%15.5f %15.5f %15.5f %15.5f %15.5f "%(d["phi"],d["theta"],d["psi"],-d["tx"],-d["ty"]))
							
					elif p[:13] == "xform.align3d":
						#phi, theta, psi, s3x, s3y, s3z, mirror, scale = get_params3D(img, p)
						if ext == "bdb":
							t = DB.get_attr(i,"xform.align3d")
							d = t.get_params("spider")
							fexp.write("%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %10d %10.3f "%(d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"],d["mirror"],d["scale"]))
							
						elif ext == "hdf":
							t =EMUtil.read_hdf_attribute(stack, "xform.align3d",i)
							d = t.get_params("spider")	
							fexp.write("%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %10d %10.3f "%(d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"],d["mirror"],d["scale"]))
							
					elif p == "ctf":
						#defocus, cs, voltage, apix, bfactor, ampcont = get_ctf(img)
						if ext == "bdb":
							t = DB.get_attr(i,"ctf")
							fexp.write("%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f"%(t.defocus, t.cs, t.voltage, t.apix, t.bfactor, t.ampcont, t.dfdiff, t.dfang))
							
						elif ext == "hdf":
							t = EMUtil.read_hdf_attribute(stack, "ctf",i)
							fexp.write("%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f"%(t.defocus, t.cs, t.voltage, t.apix, t.bfactor, t.ampcont, t.dfdiff, t.dfang))
							
					else:
						if ext == "bdb":
							fexp.write("%15s "%str(DB.get_attr(i, p)))

						elif ext == "hdf":
							fexp.write("%15s "%str(EMUtil.read_hdf_attribute(stack, p, i)))
				elif fprint:
					if p[:13] == "xform.align2d":
						#alpha, sx, sy, mirror, scale = get_params2D(img, p)
						if ext == "bdb":
							t = DB.get_attr(i,"xform.align2d")
							d = t.get_params("2D")
							print "%15.5f %15.5f %15.5f %10d %10.3f"%(d["alpha"],d["tx"],d["ty"],d["mirror"],d["scale"]),
							
						elif ext == "hdf":
							t = EMUtil.read_hdf_attribute(stack, "xform.align2d",i)
							d = t.get_params("2D")
							print "%15.5f %15.5f %15.5f %10d %10.3f"%(d["alpha"],d["tx"],d["ty"],d["mirror"],d["scale"]),
					elif p[:16] == "xform.projection":
						#phi, theta, psi, s2x, s2y = get_params_proj(img, p)
						if ext == "bdb":
							t = DB.get_attr(i,"xform.projection")
							d = t.get_params("spider")
							print "%15.5f %15.5f %15.5f %15.5f %15.5f"%(d["phi"],d["theta"],d["psi"],-d["tx"],-d["ty"]),
						elif ext == "hdf":
							t = EMUtil.read_hdf_attribute(stack, "xform.projection", i)
							d = t.get_params("spider")
							print "%15.5f %15.5f %15.5f %15.5f %15.5f"%(d["phi"],d["theta"],d["psi"],-d["tx"],-d["ty"]),
					elif p[:13] == "xform.align3d":
						#phi, theta, psi, s3x, s3y, s3z, mirror, scale = get_params3D(img, p)
						if ext == "bdb":
							t = DB.get_attr(i, "xform.align3d")
							d = t.get_params("spider")
							print "%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %10d %10.3f"%(d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"],d["mirror"],d["scale"]),
						elif ext == "hdf":
							t = EMUtil.read_hdf_attribute(stack, "xform.align3d", i)
							d = t.get_params("spider")
							print "%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %10d %10.3f"%(d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"],d["mirror"],d["scale"]),
					elif p == "ctf":
						#defocus, cs, voltage, apix, bfactor, ampcont = get_ctf(img)
						if ext == "bdb":
							t = DB.get_attr(i, "ctf")
							print "%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f"%(t.defocus, t.cs, t.voltage, t.apix, t.bfactor, t.ampcont),
						elif ext == "hdf":
							t = EMUtil.read_hdf_attribute(stack,"ctf", i)
							print "%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f"%(t.defocus, t.cs, t.voltage, t.apix, t.bfactor, t.ampcont),
					else:
						if ext == "bdb":
							print "%15s"%str(DB.get_attr(i, p)),
						elif ext == "hdf":
							print "%15s"%str(EMUtil.read_hdf_attribute(stack, p, i)),
				elif backup:
					#t = img.get_attr(p)
					#img.set_attr(p+suffix, t)
					if ext == "bdb":
						t= DB.get_attr(i, p)
						DB.set_attr(i, p+suffix, t)
					elif ext == "hdf":
						t= EMUtil.read_hdf_attribute(stack, p, i)
						EMUtil.write_hdf_attribute(stack,p+suffix, t, i)
				
				elif restore:
					if p == "xform.align2d" or p == "xform.align3d" or p == "xform.projection":
						print  "ERROR, no suffix in xform!"
						return
					#t = img.get_attr(p)
					#if ext == "bdb":
						#for i in xrange(nimage):
							#t= DB.get_attr(i,p)
					#elif ext == "hdf":
						#for i in xrange(nimage):
							#t= EMUtil.read_hdf_attribute(stack,p,i)
					if p[:13] == "xform.align2d":
						#img.set_attr(p[:13], t)
						if ext == "bdb":
							t = DB.get_attr(i, p)
							DB.set_attr(i, "xform.align2d", t)
						elif ext == "hdf":
							t = EMUtil.read_hdf_attribute(stack, p, i)
							EMUtil.write_hdf_attribute(stack, "xform.align2d", t, i)
					elif p[:16] == "xform.projection":
						#img.set_attr(p[:10], t)
						if ext == "bdb":
							t = DB.get_attr(i, p)
							DB.set_attr(i, "xform.projection", t)
						elif ext == "hdf":
							t = EMUtil.read_hdf_attribute(stack, p, i)
							EMUtil.write_hdf_attribute(stack, "xform.projection", t, i)
					elif p[:13] == "xform.align3d":
						#img.set_attr(p[:13], t)
						if ext == "bdb":
							t = DB.get_attr(i, p)
							DB.set_attr(i, "xform.align3d", t)
					elif ext == "hdf":
						for i in xrange(nimage):
							t = EMUtil.read_hdf_attribute(stack, p, i)
							EMUtil.write_hdf_attribute(stack, "xform.align3d", t, i)
					else:
						#img.set_attr(p[:-len(suffix)], t)
						if ext == "bdb":
							t = DB.get_attr(i, p)
							DB.set_attr(i, p[:-len(suffix)],t)
						elif ext == "hdf":
							t = EMUtil.read_hdf_attribute(stack, p, i)
							EMUtil.write_hdf_attribute(stack, p[:-len(suffix)], t, i)
				elif delete:
					img = EMData()
					img.read_image(stack, i, True)
					img.del_attr(p)
					write_header(stack, img, i)
		#if zero or one or randomize or rand_alpha or backup or restore or delete:
			#write_header(stack, img, i)
			if fexport != None:
				fexp.write( "\n" )
			if fprint:
				print " "
	if ext == "bdb": DB.close()

def imgstat_ccc( stacks, rad ):
	from EMAN2 import EMUtil
	from utilities import get_im, model_circle
	from statistics import ccc
	from projection import prep_vol,prgs
	from utilities	import get_params_proj

	if len(stacks)>3: ERROR("Error: ccc should be run on two stacks","imgstat_ccc",1)

	nimg1 = EMUtil.get_image_count( stacks[0] )
	nimg2 = EMUtil.get_image_count( stacks[1] )


	if nimg2==1 and get_im(stacks[0]).get_zsize()==1 and get_im(stacks[1]).get_zsize() > 1:
		print "ccc between prj and volume"
		volccc = True
		volft,kb = prep_vol( get_im(stacks[1]) )
	else:
		volccc = False


	
	nimg = max(nimg1,nimg2)
	if(nimg2<nimg1): nimg2 = 1

        imgtmp = get_im( stacks[0] )

	if rad==-1:
		if len(stacks) == 3:  mask = get_im(stacks[2])
		else:                 mask = None
	else:
		if len(stacks) == 3:    ERROR("Error: Mask radius and mask file canot be given simultaneously","imgstat_ccc",1)
		else:
			nx = imgtmp.get_xsize()
			ny = imgtmp.get_ysize()
			nz = imgtmp.get_zsize()
			mask = model_circle( rad, nx, ny, nz )

	for i in xrange(nimg):
		img1 = get_im( stacks[0], i )

		if nimg2==1:
			if volccc:
				phi,tht,psi,s2x,s2y = get_params_proj( img1 )
				img2 = prgs( volft,kb, [phi,tht,psi,-s2x,-s2y] )
			else:
				img2 = get_im( stacks[1] )
		else:
			img2 = get_im( stacks[1], i )

		val = ccc(img1, img2, mask)

		print "%6d: %10.5f" % (i, val)

def imgstat_fsc( stacks, fscfile, rad ):
	from utilities import get_im, model_circle
	from statistics import fsc_mask

	if len(stacks)>3: ERROR("Error: fsc should be run on two images","imgstat_fsc",1)

	img1 = get_im( stacks[0] )
	img2 = get_im( stacks[1] )

	nx = img1.get_xsize()
	ny = img1.get_ysize()
	nz = img1.get_zsize()

	if  EMUtil.get_image_count(stacks[0])>1: ERROR("Error: %s is an stack, fsc should be run on images","imgstat_fsc",1)

	if img2.get_xsize() != nx or img2.get_ysize() != ny or img2.get_zsize() != nz: ERROR("Error: input images has different sizes","imgstat_fsc",1)

	if rad==-1:
		if len(stacks) == 3: mask = get_im(stacks[2])
		else:                mask = None
	else:
		if len(stacks) == 3:  ERROR("Error: Mask radius and mask file canot be given simultaneously","imgstat_fsc",1)
		else:    mask = model_circle( rad, nx, ny, nz )

	fsc_mask( img1, img2, mask, filename=fscfile )
	
def imgstat_inf( stacks, rad ):
	from EMAN2 import EMUtil
	from utilities import get_im, model_circle
	if len(stacks)>2: ERROR("Error: inf should be run on one file","imgstat_inf",1)

	nimg = EMUtil.get_image_count( stacks[0] )
	img1 = get_im( stacks[0] )
	
	nx = img1.get_xsize()
	ny = img1.get_ysize()
	nz = img1.get_zsize()

	if rad==-1:
		if len(stacks) == 2:  mask = get_im(stacks[1])
		else:                 mask = None
	else:
		if len(stacks) == 2:    ERROR("Error: Mask radius and mask file canot be given simultaneously","imgstat_inf",1)
		else:			mask = model_circle( rad, nx, ny, nz )


	for i in xrange(nimg):

		img = get_im( stacks[0], i )

		[avg,sigma,fmin,fmax] = Util.infomask( img, mask, True )

		print "nx,ny,nz,avg,sigma,min,max: %6d %6d %6d %11.4e %10.5f %10.5f %10.5f" % (nx, ny, nz, avg, sigma, fmin, fmax )

def imgstat( stacks, ifccc, fscfile, pinf, rad ):
	if ifccc:
		imgstat_ccc( stacks, rad )
		return

	if len(fscfile)>0:
		imgstat_fsc( stacks, fscfile, rad )
		return

	if pinf:
		imgstat_inf( stacks, rad )
		return

def MPI_start_end(nima, nproc, myid):
	image_start = int(round(float(nima)/nproc*myid))
	image_end   = int(round(float(nima)/nproc*(myid+1)))
	return image_start, image_end

def normal_prj( prj_stack, outdir, refvol, weights, r, niter, snr, sym, verbose = 0, CTF = False, MPI=False ):
	def peak_range( nx, ctf ):
		"""
		  Find first maximum of the CTF, use CTF^2, so the sign will be ignored
		"""
		from morphology import ctf_2
		ctf = ctf_2( nx, ctf )

		for i in xrange( 1, len(ctf)-1 ):
			prev = ctf[i-1]
			curt = ctf[i]
			next = ctf[i+1]

			if curt > prev and curt > next:
				freq = float(i)/nx
				return [freq-0.03, freq+0.02]

		assert false

	from utilities     import get_image, get_im, model_circle, drop_spider_doc
	from utilities     import drop_image, get_params_proj
	from projection    import prep_vol, prgs
	from filter        import filt_ctf, filt_btwo, filt_tophatb
	from fundamentals  import fft
	from statistics    import ccc
	import os

	if(MPI and not (weights is None)):
		ERROR('Application of weights does not have MPI version', "normal_prj", 1,myid)

	if MPI:
		from mpi import mpi_comm_size, mpi_comm_rank, mpi_barrier, mpi_init, mpi_reduce, mpi_bcast, MPI_COMM_WORLD, MPI_FLOAT, MPI_SUM
		from utilities     import bcast_EMData_to_all
		nproc = mpi_comm_size( MPI_COMM_WORLD )
		myid  = mpi_comm_rank( MPI_COMM_WORLD )
	else:
		nproc = 1
		myid  = 0
		if( not (weights is None) ):
			#  This section is application of weights
			from utilities import read_text_file
			s = read_text_file(weights)
			img_number     = EMUtil.get_image_count( prj_stack )
			if(len(s) != img_number):  ERROR('Number of images does not agree with number of weights', "normal_prj", 1,myid)
			for i in xrange(img_number):
				img = get_im(prj_stack, i)
				Util.mul_scalar(img, s[i])
				img.write_image(outdir, i)
			return
	
	if myid== 0:
		if os.path.exists(outdir): ERROR('Output directory exists, please change the name and restart the program', "normal_prj", 1,myid)
    		os.mkdir(outdir)

	if MPI:
		mpi_barrier( MPI_COMM_WORLD )


	img = get_image( prj_stack )
	nx  = img.get_xsize()
	ny  = img.get_ysize()
	del img

	if r <= 0:  r = nx//2-1

	img_number     = EMUtil.get_image_count( prj_stack )
	img_node_start, img_node_end = MPI_start_end(img_number, nproc, myid )


	if(verbose == 1):  info=open( os.path.join(outdir, "progress%04d.txt" % myid), "w" )
	else:              info = None

	imgdata = EMData.read_images(prj_stack, range(img_node_start, img_node_end))

	if(verbose == 1):
		info.write( ' all images loaded\n' )
		info.flush( )

	pred = [1.0]* len(imgdata)

	from reconstruction import rec3D_MPI,rec3D_MPI_noCTF,rec3D_MPI_noCTF, recons3d_4nn_ctf, recons3d_4nn
        if refvol is None:
		fsc_file = os.path.join(outdir, "fsc_init.dat")
		vol_file = os.path.join(outdir, "vol_init.hdf")
		if  MPI:
			if  CTF:  refvol, fscc = rec3D_MPI( imgdata, snr, sym, None, fsc_file, myid )
			else:     refvol, fscc = rec3D_MPI_noCTF( imgdata, sym, None, fsc_file, myid )
			bcast_EMData_to_all( refvol, myid )
		else:
			if CTF:   refvol = recons3d_4nn_ctf( imgdata, range(len(imgdata)), snr, 1, sym)
			else:	   refvol = recons3d_4nn( imgdata, range(len(imgdata)), sym)
		if myid==0:
			refvol.write_image( vol_file )
		if(verbose == 1):
			info.write( "inital reconstructed volume written to " + vol_file + "\n" )
			info.flush()


    	mask = model_circle( r, nx, ny )
	for iter in xrange(niter) :
		refvol, kb = prep_vol( refvol )

		scales = []
		for i in xrange( len(imgdata) ) :
			exp_prj = imgdata[i].copy()

			phi,theta,psi,s2x,s2y = get_params_proj( exp_prj )

			ref_prj = filt_btwo( fft( prgs( refvol, kb, [phi, theta, psi, -s2x, -s2y] ) ), 0.01, 0.1, 0.2)

			if  CTF:
				ctf = exp_prj.get_attr( "ctf" )
				ref_prj = filt_ctf( filt_ctf( ref_prj, ctf ), ctf )
				frange = peak_range( nx, ctf)

				if exp_prj.get_attr('ctf_applied')==0.0:
					exp_prj = filt_ctf( fft(exp_prj), ctf )
				else:
					exp_prj = fft(exp_prj)
				ref_prj = filt_tophatb( ref_prj, frange[0], frange[1], False )
				exp_prj = filt_tophatb( exp_prj, frange[0], frange[1], False )
			else:
				exp_prj = fft(exp_prj)

			ref_prj = fft(ref_prj)
			exp_prj = fft(exp_prj)
			Util.mul_img(ref_prj, mask)
			Util.mul_img(exp_prj, mask)
			curtccc = ccc( ref_prj, exp_prj, mask )

			try:
				a = exp_prj.dot( ref_prj ) / exp_prj.dot(exp_prj)
			except:
				print 'exception at myid, i:', myid, i
				a = 1.0

		        scales.append( a )
			if(verbose == 1):
				info.write( "i, a, ccc:  %4d %10.5f  %10.5f\n" %(i, a, curtccc) )
				info.flush()

 	   	sum_scale = sum( scales )

		if  MPI:
			total_sum_scale = mpi_reduce( sum_scale, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD )
			total_sum_scale = mpi_bcast( total_sum_scale, 1, MPI_FLOAT, 0, MPI_COMM_WORLD)
			sum_scale = float(total_sum_scale[0])

		avg_scale = sum_scale/img_number

		assert( len(imgdata)==len(scales) )

		for i in xrange( len(imgdata) ):
			s = scales[i] / avg_scale
			imgdata[i] *= s
			pred[i] *= s

    		scale_file = os.path.join(outdir, "newscale%04d_%04d.txt" % (myid, iter))
    		drop_spider_doc( scale_file, pred )

		fsc_file = os.path.join(outdir, ( "fsc_%04d.dat" % iter ))
		vol_file = os.path.join(outdir, ( "vol_%04d.hdf" % iter ))
		if(verbose == 1):
			info.write( 'running reconstruction\n' )
			info.flush()
		if  MPI:
			if  CTF:  refvol, fscc = rec3D_MPI( imgdata, snr, sym, None, fsc_file, myid )
			else:     refvol, fscc = rec3D_MPI_noCTF( imgdata, sym, None, fsc_file, myid )
			bcast_EMData_to_all( refvol, myid )
		else:
			if CTF:   refvol = recons3d_4nn_ctf( imgdata, range(len(imgdata)), snr, 1, sym)
			else:	   refvol = recons3d_4nn( imgdata, range(len(imgdata)), sym)
		if(verbose == 1):
			info.write( 'reconstruction finished\n' )
			info.flush()

		if myid==0:
			drop_image( refvol, vol_file )
			if(verbose == 1):
				info.write( "reconstructed volume written to " + vol_file  + "\n")
				info.flush()

		if  MPI:  bcast_EMData_to_all( refvol, myid )
		[mean,sigma,fmin,fmax] = Util.infomask( refvol, None, True )
		if(verbose == 1):
			info.write( 'vol all after reconstruction, myid: %d %10.3e %10.3e %10.3e %10.3e\n' % ( myid, mean, sigma, fmin, fmax ) )
			info.flush()

	if(myid == 0):
		foutput = open( os.path.join(outdir,"weights.txt"), "w" )
	if  MPI:
		from mpi import MPI_INT, MPI_TAG_UB, mpi_recv, mpi_send
		if myid == 0:
			ltot = 0
			base = 0
			for iq in xrange( nproc ):
				if(iq == 0):
					ltot = spill_out(ltot, base, pred, 1, foutput)
				else:
					lend = mpi_recv(1, MPI_INT, iq, MPI_TAG_UB, MPI_COMM_WORLD)
					lend = int(lend[0])
					pred = mpi_recv(lend, MPI_FLOAT, iq, MPI_TAG_UB, MPI_COMM_WORLD)
					ltot = spill_out(ltot, base, pred, 1, foutput)
				base += len(pred)
		else:
			mpi_send([len(pred)], 1, MPI_INT, 0, MPI_TAG_UB, MPI_COMM_WORLD)
			mpi_send(pred, len(pred), MPI_FLOAT, 0, MPI_TAG_UB, MPI_COMM_WORLD)
	else:
		ltot = 0
		base = 0
		ltot = spill_out(ltot, base, pred, 1, foutput)

"""
def incvar(prefix, nfile, nprj, output, fl, fh, radccc, writelp, writestack):
	from statistics import variancer, ccc
	from string     import atoi, replace, split, atof
	from utilities  import get_im, circumference, model_circle, drop_image
	from filter     import filt_btwl
	from math       import sqrt
	import os

	all_varer = variancer()
	odd_varer = variancer()
	eve_varer = variancer()

	filname = prefix + "0000.hdf"
	img = get_im(filname, 0)
	n = img.get_xsize()
	
	if os.path.exists(output):		os.system("rm -f "+output)
	if os.path.exists('stack_'+output):	os.system("rm -f stack_"+output)
	if os.path.exists('odd_stack_'+output):	os.system("rm -f odd_stack_"+output)
	if os.path.exists('eve_stack_'+output):	os.system("rm -f eve_stack_"+output)
	if os.path.exists('avg_'+output):	os.system("rm -f avg_"+output)

	cccmask = model_circle(radccc, n, n, n)
	scale = sqrt( nprj )
	totnimg = 0
	iwrite = 0
	radcir = n/2-3
	for i in xrange(nfile):
		filename = prefix + '%04d.hdf'%i 

		print 'loading file ', filename
		nimg = EMUtil.get_image_count( filename )
		for j in xrange(nimg):
			
			img = get_im( filename, j )
			img *= scale
			img = circumference( img, radcir, radcir+2 )
			img = filt_btwl(img, fl, fh)
						
			if writelp:
				img.write_image("btwl_cir_"+filename, j)

			if totnimg%2==0: odd_varer.insert(img)
			else: eve_varer.insert(img)

			all_varer.insert(img)

			totnimg += 1

			if totnimg%100==0:
				odd_var = odd_varer.getvar()
				eve_var = eve_varer.getvar()
				all_var = all_varer.getvar()
			
				if writestack:
					odd_var.write_image( 'odd_stack_' + output, iwrite )
					eve_var.write_image( 'eve_stack_' + output, iwrite )
					all_var.write_image( 'stack_' + output, iwrite )
					iwrite += 1  
				print 'ntot, ccc: %6d %10.3f' % (totnimg, ccc(odd_var, eve_var, cccmask))  

	all_var = all_varer.getvar()
	odd_var = odd_varer.getvar()
	eve_var = eve_varer.getvar()
	print 'ntot, ccc: %6d %10.3f' % (totnimg, ccc(odd_var, eve_var, cccmask))  

	avg = all_varer.getavg()
	avg.write_image( 'avg_' + output, 0 )

	if writestack:
		all_var.write_image( 'stack_' + output, iwrite )
		odd_var.write_image( 'odd_stack_' + output, iwrite )
		eve_var.write_image( 'eve_stack_' + output, iwrite )

	all_var = circumference( all_var, radcir, radcir+2 )
	all_var.write_image( output, 0 )
"""

class file_set :

	def __init__( self, files ):
		nfile = len(files)
		self.files = files
		self.fends = [None] * nfile

		totimg = 0
		for i in xrange(nfile):
			totimg += EMUtil.get_image_count( self.files[i] )
			self.fends[i] = totimg		

	def nimg(self):
		return self.fends[-1]

	def get(self, imgid):
		assert imgid < self.fends[-1]

		ifile = 0
		while imgid >= self.fends[ifile]:
			ifile += 1

		if ifile==0:
			return self.files[0], imgid

		return self.files[ifile], imgid - self.fends[ifile-1]

def defvar(files, outdir, fl, aa, radccc, frepa = "default", pca=False, pcamask=None, pcanvec=None):
	from utilities  import get_im, get_image, circumference, model_blank
	from filter     import filt_tanl
	from math       import sqrt
	import os
	
	if os.path.exists(outdir): ERROR('Output directory exists, please change the name and restart the program', " defvar", 1)
	os.mkdir(outdir)

	#finf = open( outdir + "/var_progress.txt", "w" )
	print "  START "

	img = get_image( files[0] )
	nx = img.get_xsize()
	ny = img.get_ysize()
	nz = img.get_zsize()
	if(frepa == "None"):  repair = False
	else:
		repair = True
		if(frepa == "default"):
			#from math import sqrt
			#from utilities import model_gauss
			#rota = model_gauss(sqrt(2.0)*nx, nx,ny,nz)
			#Util.mul_scalar( rota, 1.0/(rota.get_value_at(nx//2, ny//2, nz//2)) )
			from utilities import model_blank
			rota = model_blank(nx, ny, nz, 1.0)
		else:   rota = get_im(frepa)

	radcir = min(nx,ny,nz)//2 - 2

	if pca :
		from statistics import pcanalyzer
		pcamask = get_im( pcamask)
		pcaer = pcanalyzer(pcamask, outdir, pcanvec, False)

	avgfile  = os.path.join(outdir, "avg.hdf")
	varfile  = os.path.join(outdir, "var.hdf")
	varfileE = os.path.join(outdir, "varE.hdf")
	avgfileE = os.path.join(outdir, "avgE.hdf")
	varfileO = os.path.join(outdir, "varO.hdf")
	avgfileO = os.path.join(outdir, "avgO.hdf")

	if(radccc < 1):  radcir = min(nx,ny,nz)//2-2
	else:            radcir = radccc

        nfiles = len( files )

	avg1 = model_blank(nx,ny,nz)
	avg2 = model_blank(nx,ny,nz)

	total_img = 0
	mf = 0
	for f in files:
		nimg = EMUtil.get_image_count( f )
		#print f," A  ",nimg
		mf += 1
		for i in xrange(nimg):
			img = get_im( f, i )
			if(fl > 0.0):
				img = filt_tanl( img, fl, aa )
			if(repair):
				Util.div_img(img, rota) #img = circumference(Util.divn_img(img, rota), radcir)
				if pca:
					pc   = Util.infomask(img, pcamask, True)
					img -= pc[0]
					img *= (refstat[1]/pc[1])
			if(total_img%2 == 0):	Util.add_img(avg1, img)
			else:			Util.add_img(avg2, img)
			total_img += 1

	avg = Util.addn_img(avg1, avg2)
	Util.mul_scalar(avg, 1.0/float(total_img))
	"""
	Util.mul_scalar(avg1, 1.0/float(total_img//2+total_img%2 - 1 ))
	avg1.write_image(avgfileE)
	Util.mul_scalar(avg2, 1.0/float(total_img//2 - 1) )
	avg2.write_image(avgfileO)
	"""
	avg.write_image( avgfile)

	del avg1, avg2

	from utilities import model_circle
	#cccmask = model_circle(radccc, nx, ny, nz)
	var1 = model_blank(nx,ny,nz)
	var2 = model_blank(nx,ny,nz)
	for f in files:
		nimg = EMUtil.get_image_count( f )
		for i in xrange(nimg):
			img = get_im( f, i )
			#img = circumference( img, radcir )
			if(fl > 0.0): img = filt_tanl( img, fl, aa)
			if(repair):  Util.div_img(img, rota) #img = circumference(Util.divn_img(img, rota), radcir)
			if pca:
				pc = Util.infomask(img, pcamask, True)
				img -= pc[0]
				img *= (refstat[1]/pc[1])
				#img += refstat[1]
			if pca : pcaer.insert(img)
			Util.sub_img(img, avg)
			if(total_img%2 == 0): Util.add_img2(var1, img)
			else:                 Util.add_img2(var2 , img)

	var = Util.addn_img(var1, var2)
	Util.mul_scalar(var, 1.0/float(total_img-1) )
	"""
	Util.mul_scalar(var1, 1.0/float(total_img//2+total_img%2 - 1 ))
	var1.write_image(varfileE)
	Util.mul_scalar(var2, 1.0/float(total_img//2 - 1) )
	var2.write_image(varfileO)
	"""
	circumference(var, radcir-1).write_image( varfile )
	del var, var1, var2#, cccmask

	if pca:
		assert not(avg is None)

		pcaer.setavg( avg )

		eigs = pcaer.analyze()
		eigfile = os.path.join(outdir, "eigvol.hdf")
		for i in xrange( len(eigs) ):
			eigs[i].write_image( eigfile, i )

def var_mpi(files, outdir, fl, aa, radccc, frepa = "default", pca=False, pcamask=None, pcanvec=None):
	from string     import atoi, replace, split, atof
	from utilities  import get_im, circumference, model_circle, model_blank
	from utilities  import bcast_EMData_to_all, reduce_EMData_to_root
	from filter     import filt_tanl
	from mpi        import mpi_comm_rank, mpi_comm_size, mpi_barrier, mpi_bcast, mpi_reduce
	from mpi        import MPI_COMM_WORLD, MPI_INT, MPI_FLOAT, MPI_SUM
	import os
	"""
	  writelp means overwrite original stacks with repaired ones
	  writestack means write new stacks with low-pass filtered data
	"""

        myid = mpi_comm_rank( MPI_COMM_WORLD )
        ncpu = mpi_comm_size( MPI_COMM_WORLD )
	main_node=0
	
	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "mref_ali2d_MPI ", 1)
	mpi_barrier(MPI_COMM_WORLD)
        if myid== main_node:   os.mkdir(outdir)

	mpi_barrier( MPI_COMM_WORLD )

	if( myid == main_node ):
		print "  START "
		img = get_im(files[0])
		nx = img.get_xsize()
		ny = img.get_ysize()
		nz = img.get_zsize()
		#del img
	else:
		nx = 0
		ny = 0
		nz = 0
	nx = mpi_bcast(nx, 1, MPI_INT, 0, MPI_COMM_WORLD)
	nx = int(nx[0])
	ny = mpi_bcast(ny, 1, MPI_INT, 0, MPI_COMM_WORLD)
	ny = int(ny[0])
	nz = mpi_bcast(nz, 1, MPI_INT, 0, MPI_COMM_WORLD)
	nz = int(nz[0])
	if(frepa == "None"):  repair = False
	else:
		repair = True
		if(frepa == "default"):
			#from math import sqrt
			#from utilities import model_gauss
			#rota = model_gauss(sqrt(2.0)*nx, nx,ny,nz)
			#Util.mul_scalar( rota, 1.0/(rota.get_value_at(nx//2, ny//2, nz//2)) )
			from utilities import model_blank
			rota = model_blank(nx, ny, nz, 1.0)
		else:   rota = get_im(frepa)

	if pca:
		from statistics import pcanalyzer
		if(myid == 0):  pcamask = get_im( pcamask)
		else:           pcamask = model_blank(nx,ny,nz)
		bcast_EMData_to_all(pcamask, myid)
		pcaer = pcanalyzer(pcamask, outdir, pcanvec, True)
		if( myid == 0 ):  refstat = Util.infomask(img, pcamask, True)
		else:             refstat = [0.0,0.0,0.0,0.0]
		refstat = mpi_bcast(refstat, 4, MPI_FLOAT, 0, MPI_COMM_WORLD)
		refstat = map(float, refstat)

	avgfile  = os.path.join(outdir, "avg.hdf")
	varfile  = os.path.join(outdir, "var.hdf")
	varfileE = os.path.join(outdir, "varE.hdf")
	avgfileE = os.path.join(outdir, "avgE.hdf")
	varfileO = os.path.join(outdir, "varO.hdf")
	avgfileO = os.path.join(outdir, "avgO.hdf")

	if(radccc < 1):  radcir = min(nx,ny,nz)//2-2
	else:            radcir = radccc

        nfiles = len( files )
	if(nfiles < ncpu):
		ERROR('Number of files less than number of processors specified, reduce number of processors', " var_mpi", 1, myid)
		
	file_start, file_end = MPI_start_end(nfiles, ncpu, myid)

	#ndump = 100 # write out ccc after each 100 steps

	iwrite = 0
	istack = 0
	iprint = 0
	iadded = 0

	avg1 = model_blank(nx,ny,nz)
	avg2 = model_blank(nx,ny,nz)

	total_img = 0
	for ifile in xrange(file_start, file_end):
		nimg = EMUtil.get_image_count( files[ifile] )
		#print myid," A  ",files[ifile],"   ",nimg
		for i in xrange(nimg):
			img = get_im( files[ifile], i )
			#img = circumference( img, radcir )
			if(fl > 0.0):
				img = filt_tanl( img, fl, aa )
			if(repair):
				Util.div_img(img, rota) #img = circumference(Util.divn_img(img, rota), radcir)
				if pca:
					pc   = Util.infomask(img, pcamask, True)
					img -= pc[0]
					img *= (refstat[1]/pc[1])
			if(total_img%2 == 0):	Util.add_img(avg1, img)
			else:			Util.add_img(avg2, img)
			total_img += 1
	reduce_EMData_to_root(avg1, myid)
	reduce_EMData_to_root(avg2, myid)
	total_img = mpi_reduce(total_img, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
	if( myid == 0) :
		total_img = int(total_img[0])
		avg = Util.addn_img(avg1, avg2)
		Util.mul_scalar(avg, 1.0/float(total_img))
	else:    avg = model_blank(nx,ny,nz)
	bcast_EMData_to_all( avg, myid )
	if( myid == 0 ):
		#Util.mul_scalar(avg1, 1.0/float(total_img//2+total_img%2 - 1 ))
		#avg1.write_image(avgfileE)
		#Util.mul_scalar(avg2, 1.0/float(total_img//2 - 1) )
		#avg2.write_image(avgfileO)
		avg.write_image( avgfile)

	del avg1, avg2

	var1 = model_blank(nx,ny,nz)
	var2 = model_blank(nx,ny,nz)
	for ifile in xrange(file_start, file_end):
		nimg = EMUtil.get_image_count( files[ifile] )
		#print myid," V  ",files[ifile],"   ",nimg
		for i in xrange(nimg):
			img = get_im( files[ifile], i )
			#img = circumference( img, radcir )
			if(fl > 0.0): img = filt_tanl( img, fl, aa)
			if(repair):  Util.div_img(img, rota) #img = circumference(Util.divn_img(img, rota), radcir)
			if pca:
				pc = Util.infomask(img, pcamask, True)
				img -= pc[0]
				img *= (refstat[1]/pc[1])
				#img += refstat[1]
			if pca : pcaer.insert(img)
			Util.sub_img(img, avg)
			if(total_img%2 == 0): Util.add_img2(var1, img)
			else:                 Util.add_img2(var2 , img)

	reduce_EMData_to_root(var1, myid)
	reduce_EMData_to_root(var2, myid)
	if( myid == 0):
		var = Util.addn_img(var1, var2)
		Util.mul_scalar(var, 1.0/float(total_img-1) )
	else:    var = model_blank(nx,ny,nz)
	bcast_EMData_to_all( var, myid )
	if(  (myid == 0)):
		#Util.mul_scalar(var1, 1.0/float(total_img//2+total_img%2 - 1 ))
		#circumference(var1, radcir-1).write_image(varfileE)
		#Util.mul_scalar(var2, 1.0/float(total_img//2 - 1) )
		#circumference(var2, radcir-1).write_image(varfileO)

		circumference(var, radcir-1).write_image( varfile )
		del var1, var2
	del var

	if pca:
		assert not(avg is None)

		pcaer.setavg( avg )

		eigs = pcaer.analyze()

		if myid==0:
			eigfile = os.path.join(outdir, "eigvol.hdf")
			for i in xrange( len(eigs) ):
				eigs[i].write_image( eigfile, i )

def factcoords_vol( vol_stacks, avgvol_stack, eigvol_stack, prefix, rad = -1, neigvol = -1, fl=0.0, aa=0.0, MPI=False):
	from utilities import get_im, model_circle, model_blank

	if MPI:
		from mpi import mpi_comm_rank, mpi_comm_size, MPI_COMM_WORLD
		ncpu = mpi_comm_size(MPI_COMM_WORLD)
		myid = mpi_comm_rank(MPI_COMM_WORLD)
	else:
		ncpu = 1
		myid = 0


        if( myid == 0 ):
		foutput = open( prefix+".txt", "w" )

	
	if(neigvol < 0):
		eigvols = EMData.read_images(eigvol_stack)
		neigvol = len(eigvols)
	else:
		eigvols = EMData.read_images(eigvol_stack,range(neigvol))

	from math import sqrt
	eigvals = [0.0]*neigvol
	for j in xrange(neigvol):
		eigvals[j] = sqrt( eigvols[j].get_attr_default('eigval',1.0) )
		Util.mul_scalar(eigvols[j] , eigvals[j])
	if( avgvol_stack != None):
		avgvol = get_im( avgvol_stack )

	nx = eigvols[0].get_xsize()
	ny = eigvols[0].get_ysize()
	nz = eigvols[0].get_zsize()

	m = model_circle( rad, nx, ny, nz )
	files = file_set( vol_stacks )
	vol_bgn,vol_end = MPI_start_end( files.nimg(), ncpu, myid )

	d = []
	for i in xrange( vol_bgn, vol_end ):
		fname,imgid = files.get( i )
		exp_vol = get_im( fname, imgid )
		if(avgvol_stack != None):  
			Util.sub_img(exp_vol, avgvol)

		for j in xrange( neigvol ):
			d.append( exp_vol.cmp( "ccc", eigvols[j], {"negative":0, "mask":m} )*eigvals[j] )

	if  MPI:
		from mpi import MPI_INT, MPI_FLOAT, MPI_TAG_UB, MPI_COMM_WORLD, mpi_recv, mpi_send
		if myid == 0:
			ltot = 0
			base = 0
			for iq in xrange( ncpu ):
				if(iq == 0):
					ltot = spill_out(ltot, base, d, neigvol, foutput)
				else:
					lend = mpi_recv(1, MPI_INT, iq, MPI_TAG_UB, MPI_COMM_WORLD)
					lend = int(lend[0])
					d = mpi_recv(lend, MPI_FLOAT, iq, MPI_TAG_UB, MPI_COMM_WORLD)
					ltot = spill_out(ltot, base, d, neigvol, foutput)
				base += len(d)/neigvol
		else:
			mpi_send([len(d)], 1, MPI_INT, 0, MPI_TAG_UB, MPI_COMM_WORLD)
			mpi_send(d, len(d), MPI_FLOAT, 0, MPI_TAG_UB, MPI_COMM_WORLD)
	else:
		ltot = 0
		base = 0
		ltot = spill_out(ltot, base, d, neigvol, foutput)

def factcoords_prj( prj_stacks, avgvol_stack, eigvol_stack, prefix, rad, neigvol, fl=0.0, aa=0.0, CTF = False, MPI=False):
	from utilities    import get_im, get_image, model_circle, model_blank, get_params_proj
	from projection   import prgs, prep_vol
	from filter       import filt_ctf, filt_tanl
	from statistics   import im_diff

	if MPI:
		from mpi import mpi_comm_rank, mpi_comm_size, MPI_COMM_WORLD
		ncpu = mpi_comm_size( MPI_COMM_WORLD )
		myid = mpi_comm_rank( MPI_COMM_WORLD )
	else:
		ncpu = 1
		myid = 0

        if(myid == 0):
		foutput = open( prefix+".txt", "w" )

	nx = get_im( prj_stacks[0] ).get_xsize()
	ny = nx

	avgvol = get_im( avgvol_stack )
	a = Util.infomask(avgvol, model_circle(int(rad), nx, nx, nx), False)
	avgvol -= a[0]
	avgvol,kb = prep_vol( avgvol )

	if neigvol==-1:  neigvol = EMUtil.get_image_count( eigvol_stack )
	# average volumes and eigen volumes.
	eigvols = [None]*neigvol
	eigvals = [0.0]*neigvol
	from math import sqrt
	for j in xrange(neigvol):
		eigvols[j] = get_im(eigvol_stack, j)
		eigvals[j] = sqrt( eigvols[j].get_attr('eigval') )
		eigvols[j], kb = prep_vol( eigvols[j] )

	m = model_circle( int(rad), nx, ny )

	files = file_set( prj_stacks )
	nprj  = files.nimg()
	img_bgn, img_end = MPI_start_end( nprj, ncpu, myid )
	#ltot = -1
	d = []
	for i in xrange( img_bgn, img_end ):
		fname,imgid = files.get(i)
		#if(i%1000 == 0):  print  "  ",myid,"   ",i
		exp_prj = get_im( fname, imgid )

		phi,theta,psi,s2x,s2y = get_params_proj(exp_prj)
		if CTF:  ctf = exp_prj.get_attr("ctf")

		ref_prj = prgs( avgvol, kb, [phi, theta, psi, -s2x, -s2y] )
		if  CTF:  ref_prj = filt_ctf( ref_prj, ctf )
		if(fl > 0.0):  exp_prj = filt_tanl(exp_prj, fl, aa)
		#ltot += 1
		#ref_prj.write_image("projection.hdf",ltot)
		diff,a,b = im_diff( ref_prj, exp_prj, m)
		#nrmd = diff.cmp( "dot", diff, {"negative":0, "mask":m} )  #CHANGED HERE
		#diff.write_image("difference.hdf",ltot)

		for j in xrange( neigvol ) :

			ref_eigprj = prgs( eigvols[j], kb, [phi, theta, psi, -s2x, -s2y] )
			if  CTF:  ref_eigprj = filt_ctf( ref_eigprj, ctf )
			#ref_eigprj.write_image("eigen.hdf",ltot)

			#d.append( diff.cmp( "dot", ref_eigprj, {"negative":0, "mask":m} )*eigvals[j]/nrmd )   CHANGED HERE
			d.append( diff.cmp( "ccc", ref_eigprj, {"negative":0, "mask":m} ))#*eigvals[j] )
        		#print  i,j,d[-1]
	if  MPI:
		from mpi import MPI_INT, MPI_FLOAT, MPI_TAG_UB, MPI_COMM_WORLD, mpi_recv, mpi_send
		if myid == 0:
			ltot = 0
			base = 0
			for iq in xrange( ncpu ):
				if(iq == 0):
					ltot = spill_out(ltot, base, d, neigvol, foutput)
				else:
					lend = mpi_recv(1, MPI_INT, iq, MPI_TAG_UB, MPI_COMM_WORLD)
					lend = int(lend[0])
					d = mpi_recv(lend, MPI_FLOAT, iq, MPI_TAG_UB, MPI_COMM_WORLD)
					ltot = spill_out(ltot, base, d, neigvol, foutput)
				base += len(d)/neigvol
		else:
			mpi_send([len(d)], 1, MPI_INT, 0, MPI_TAG_UB, MPI_COMM_WORLD)
			mpi_send(d, len(d), MPI_FLOAT, 0, MPI_TAG_UB, MPI_COMM_WORLD)
	else:
		ltot = 0
		base = 0
		ltot = spill_out(ltot, base, d, neigvol, foutput)

def spill_out(ltot, base, d, neigvol, foutput):
	loc = 0
	for i in xrange(len(d)//neigvol):
		for j in xrange( neigvol):
			foutput.write( "    %e" % d[loc] )
			ltot += 1
			loc += 1
		foutput.write( "    %d" % (i+base) )
		foutput.write( "\n" )
	foutput.flush()
	return  ltot

def refvol( vollist, fsclist, output, mask ):
	from utilities     import get_image, read_fsc
	from fundamentals  import rops_table
	from math          import sqrt
	from filter        import filt_tanl, fit_tanh, filt_table, filt_vols
	from morphology    import threshold

	nvol = len(vollist)
	assert len(fsclist)==nvol

	fscs = [None]*nvol
	vols = [None]*nvol
	for i in xrange(nvol):
		fscs[i] = read_fsc( fsclist[i] )
		vols[i] = get_image( vollist[i] )
		print 'rawvol, resolution: ', vollist[i], fsclist[i]

	m    = get_image( mask )
	volfs = filt_vols( vols, fscs, m )

	for i in xrange(nvol):
		volfs[i].write_image( output, i )

# -- K-means main ---------------------------------------------------------------------------
# K-means main driver
def k_means_main(stack, out_dir, maskname, opt_method, K, rand_seed, maxit, trials, critname,
		 CTF = False, F = 0, T0 = 0, MPI = False, CUDA = False, DEBUG = False, flagnorm = False,
		 init_method = 'rnd'):
	# Common
	from utilities   import print_begin_msg, print_end_msg, print_msg, file_type, running_time
	from statistics  import k_means_locasg2glbasg
	from time        import time
	import sys, os
	#import time
	if MPI:
		from mpi        import mpi_init, mpi_comm_size, mpi_comm_rank, mpi_barrier
		from mpi        import MPI_COMM_WORLD, MPI_INT, mpi_bcast
		from mpi	import MPI_FLOAT, MPI_INT, mpi_recv, mpi_send, MPI_TAG_UB
		from utilities  import bcast_number_to_all, recv_EMData, send_EMData
		
	if CUDA:
		from statistics import k_means_cuda_init_open_im, k_means_cuda_headlog
		from statistics import k_means_cuda_export
		if MPI: from statistics import k_means_CUDA_MPI
		else:   from statistics import k_means_CUDA, k_means_SSE_CUDA
	else:
		from statistics import k_means_init_open_im, k_means_open_im, k_means_headlog
		from statistics import k_means_criterion, k_means_export
		if MPI: from statistics import k_means_cla, k_means_SSE_MPI
		else:   from statistics import k_means_cla, k_means_SSE

	ext = file_type(stack)
	if ext == 'txt': TXT = True
	else:            TXT = False

	if (T0 == 0 and F != 0) or (T0 != 0 and F == 0):
		ERROR('Ambigues parameters F=%f T0=%f' % (F, T0), 'k_means_main', 1)
		sys.exit()

	if MPI:
		sys.argv  = mpi_init(len(sys.argv), sys.argv)
		ncpu      = mpi_comm_size(MPI_COMM_WORLD)
		myid      = mpi_comm_rank(MPI_COMM_WORLD)
		main_node = 0
		mpi_barrier(MPI_COMM_WORLD)

		if os.path.exists(out_dir): ERROR('Output directory exists, please change the name and restart the program', "k_means_main ", 1)
		mpi_barrier(MPI_COMM_WORLD)

	else:
		if os.path.exists(out_dir): ERROR('Output directory exists, please change the name and restart the program', "k_means_main ", 1)

	if MPI and not CUDA:
		
		if myid == main_node:	print_begin_msg('k-means')
		mpi_barrier(MPI_COMM_WORLD)
		
		LUT, mask, N, m, Ntot = k_means_init_open_im(stack, maskname)
		N_min = N
		
		
		
		IM, ctf, ctf2         = k_means_open_im(stack, mask, CTF, LUT, flagnorm)
		
		
		if myid == main_node: 
			k_means_headlog(stack, out_dir, opt_method, N, K, 
						      critname, maskname, ncpu, maxit, CTF, T0, 
						      F, rand_seed, ncpu, m)
			t_start = time()
		
		[Cls, assign, Je] = k_means_SSE_MPI(IM, mask, K, rand_seed, maxit, 
					1, [CTF, ctf, ctf2], F, T0, DEBUG, init_method, myid = myid, main_node = main_node, jumping = 1)
					
				
		
		from statistics import k_means_SSE_combine
		[ assign_return, r_Cls, je_return, n_best] = k_means_SSE_combine(Cls, assign, Je, N, K, ncpu, myid, main_node)
		mpi_barrier(MPI_COMM_WORLD)
		if myid == main_node:
		
			if n_best == -1:
				print_msg('>>> WARNING: All trials resulted in empty clusters, STOP k-means.\n\n')
				print_end_msg('k-means MPI end')
				running_time(t_start)	
			#print "assign_return===", assign_return[10:20], "cls_n return==", r_Cls['n'], "Ji==", r_Cls['Ji'], "ave size ==", r_Cls['ave'][0].get_xsize()
			else:
				for i in xrange( ncpu ):
					if( je_return[i] <0 ):
						print_msg('> Trials: %5d    resulted in empty clusters  \n' % (i) )
					else:
						print_msg('> Trials: %5d    criterion: %11.6e  \n' % (i, je_return[i]) )
				running_time(t_start)
				crit = k_means_criterion(r_Cls, critname)
				glb_assign = k_means_locasg2glbasg(assign_return, LUT, Ntot)
				k_means_export(r_Cls, crit, glb_assign, out_dir, -1, TXT)
				print_end_msg('k-means MPI end')
	
	
	
	#don't touch below code

	elif CUDA and not MPI: # added 2009-02-20 16:27:26 # modify 2009-09-23 13:52:29
		print_begin_msg('k-means')
		LUT, mask, N, m, Ntot = k_means_cuda_init_open_im(stack, maskname)
		k_means_cuda_headlog(stack, out_dir, 'cla', N, K, maskname, maxit, T0, F, rand_seed, 1, m)
		if   opt_method == 'cla':
			k_means_CUDA(stack, mask, LUT, m, N, Ntot, K, maxit, F, T0, rand_seed, out_dir, TXT, 1, flagnorm=flagnorm)
		else:
			k_means_SSE_CUDA(stack, mask, LUT, m, N, Ntot, K, maxit, F, T0, rand_seed, out_dir, TXT, 1, flagnorm=flagnorm)
		print_end_msg('k-means')
	#don't touch below code
	elif MPI and CUDA: # added 2009-09-22 14:34:45
		print "tao mpi and cuda"
		LUT, mask, N, m, Ntot = k_means_cuda_init_open_im(stack, maskname)
		if myid == main_node:
			print_begin_msg('k-means')
			k_means_cuda_headlog(stack, out_dir, 'cuda', N, K, maskname, maxit, T0, F, rand_seed, ncpu, m)
		k_means_CUDA_MPI(stack, mask, LUT, m, N, Ntot, K, maxit, F, T0, rand_seed, myid, main_node, ncpu, out_dir, TXT, 1, flagnorm=flagnorm)
		if myid == main_node:
			print_end_msg('k-means')

		mpi_barrier(MPI_COMM_WORLD)

	else:
		print_begin_msg('k-means')
		LUT, mask, N, m, Ntot = k_means_init_open_im(stack, maskname)
		IM, ctf, ctf2         = k_means_open_im(stack, mask, CTF, LUT, flagnorm)
		k_means_headlog(stack, out_dir, opt_method, N, K, critname, maskname, trials, maxit, 
					CTF, T0, F, rand_seed, 1, m, init_method)
		
		if   opt_method == 'cla':
			[Cls, assign] = k_means_cla(IM, mask, K, rand_seed, maxit, 
					trials, [CTF, ctf, ctf2], F, T0, DEBUG, init_method)
		elif opt_method == 'SSE':
			[Cls, assign] = k_means_SSE(IM, mask, K, rand_seed, maxit, 
					trials, [CTF, ctf, ctf2], F, T0, DEBUG, init_method)
		else:
			ERROR('opt_method %s unknown!' % opt_method, 'k_means_main', 1,myid)
			sys.exit()
		crit = k_means_criterion(Cls, critname)
		glb_assign = k_means_locasg2glbasg(assign, LUT, Ntot)
		k_means_export(Cls, crit, glb_assign, out_dir, -1, TXT)

		print_end_msg('k-means')


# -- K-means groups ---------------------------------------------------------------------------
			
# K-means groups driver
def k_means_groups(stack, out_file, maskname, opt_method, K1, K2, rand_seed, maxit, trials, CTF=False, F=0.0, T0=0.0, MPI=False, CUDA=False, DEBUG=False, flagnorm=False):

	#import os
	#if os.path.exists(out_file):
		#ERROR('Output directory exists, please change the name and restart the program', "k_means_groups", 1)
	
	if MPI:
	 	#print "MPI version of kmeans group is under development"
		#sys.exit()
		from statistics import k_means_groups_MPI
		k_means_groups_MPI(stack, out_file, maskname, opt_method, K1, K2, rand_seed, maxit, trials, CTF, F, T0, flagnorm)
	elif CUDA:
		import os
		if os.path.exists(out_file):
			ERROR('Output directory exists, please change the name and restart the program', "k_means_groups_CUDA", 1)
		from statistics import k_means_groups_CUDA
		k_means_groups_CUDA(stack, out_file, maskname, K1, K2, rand_seed, maxit, F, T0)
	else:
		import os
		if os.path.exists(out_file):
			ERROR('Output directory exists, please change the name and restart the program', "k_means_groups_serial", 1)
		from statistics import k_means_groups_serial
		k_means_groups_serial(stack, out_file, maskname, opt_method, K1, K2, rand_seed, maxit, trials, CTF, F, T0, DEBUG, flagnorm)



# 2008-12-08 12:46:11 JB
# Plot angles distribution on a hemisphere from a list of given projection
def plot_projs_distrib(stack, outplot):
	from projection import plot_angles
	from utilities  import get_params_proj, file_type
	import sys

	N    = EMUtil.get_image_count(stack)
	ext  = file_type(stack)
	if ext == 'bdb':
		from EMAN2db import db_open_dict
		DB = db_open_dict(stack)
	agls = []
	for n in xrange(N):
		im = EMData()
		im.read_image(stack, n, True)
		try:
			p0, p1, p2, p3, p4 = get_params_proj(im)
		except RuntimeError:
			print 'Projection #%d from %s has no angles set!' % (n, stack)
			sys.exit()
		agls.append([p0, p1])

	if ext == 'bdb': DB.close()

	im = plot_angles(agls)
	im.write_image(outplot, 0)

# 2008-12-08 12:46:46 JB
# Wrap for the HAC part of py_cluster in the statistics.py file
def HAC_clustering(stack, dendoname, maskname, kind_link, kind_dist, flag_diss):
	from statistics   import ccc, py_cluster_HierarchicalClustering
	from copy         import deepcopy
	from utilities    import get_im, get_params2D, get_params3D
	from fundamentals import rot_shift2D, rot_shift3D

	N    = EMUtil.get_image_count(stack)
	if maskname != None: mask = get_im(maskname)
	else:                mask = None

	IM = EMData.read_images(stack)
	ny = IM[0].get_ysize()
	nz = IM[0].get_zsize()
	for n in xrange(N):
		# 3D object
		if nz > 1:
			phi, theta, psi, s3x, s3y, s3z, mirror, scale = get_params3D(IM[n])
			IM[n]  = rot_shift3D(IM[n], phi, theta, psi, s3x, s3y, s3z, scale)
			if mirror: IM[n].process_inplace('xform.mirror', {'axis':'x'})
		# 2D object
		elif ny > 1:
			alpha, sx, sy, mirror, scale = get_params2D(IM[n])
			IM[n] = rot_shift2D(IM[n], alpha, sx, sy, mirror, scale)

		if mask != None: Util.mul_img(IM[n], mask)
		IM[n].set_attr('ID_hclus', n)

	if kind_dist   == 'SqEuc':
		if flag_diss: cl = py_cluster_HierarchicalClustering(IM, lambda x,y: -x.cmp("SqEuclidean", y), linkage = kind_link)
		else:        cl = py_cluster_HierarchicalClustering(IM, lambda x,y: x.cmp("SqEuclidean", y), linkage = kind_link)
	elif kind_dist == 'CCC':
		if flag_diss: cl = py_cluster_HierarchicalClustering(IM, lambda x,y: -ccc(x, y, mask), linkage = kind_link)
		else:        cl = py_cluster_HierarchicalClustering(IM, lambda x,y: ccc(x, y, mask), linkage = kind_link)

	k     = N
	Dendo = {}
	doc   = open(dendoname + '.txt', 'w')
	for val in xrange(0, 10000):
	    if flag_diss: th  = -(val / 1000.0)
	    else:         th  = val / 1000.0

	    res = cl.getlevel(th)

	    GP = []
	    for gp in res:
		OBJ = []
		for obj in gp:
		    OBJ.append(obj.get_attr('ID_hclus'))
		GP.append(OBJ)

	    newk = len(GP)
	    if newk != k:
		k = newk
		doc.write('%7.3f %d\n' % (th, k))
		Dendo[k] = deepcopy(GP)
	doc.close()

	import pickle
	f = open(dendoname + '.dendo', 'w')
	pickle.dump(Dendo, f)
	f.close()

# 2008-12-08 15:20:24 JB
# Compute the averages from the dendogram given by the function HAC_clustering
def HAC_averages(stack, dendoname, avename, K):
	from utilities    import get_im, get_params2D, get_params3D
	from fundamentals import rot_shift2D, rot_shift3D
	from utilities    import model_blank
	import sys
	
	N    = EMUtil.get_image_count(stack)

	try:
		import pickle
		f = open(dendoname, 'r')
		Dendo = pickle.load(f)
		f.close()
	except:
		print 'Impossible to read dendogram structure.'
		sys.exit()

	list_k = Dendo.keys()
	if K not in list_k:
		print 'Dendogram does not contain the draw for K=%d' % K
		sys.exit()

	IM = EMData.read_images(stack)
	nx = IM[0].get_xsize()
	ny = IM[0].get_ysize()
	nz = IM[0].get_zsize()
	for n in xrange(N):
		# 3D object
		if nz > 1:
			phi, theta, psi, s3x, s3y, s3z, mirror, scale = get_params3D(IM[n])
			IM[n]  = rot_shift3D(IM[n], phi, theta, psi, s3x, s3y, s3z, scale)
			if mirror: IM[n].process_inplace('xform.mirror', {'axis':'x'})
		# 2D object
		elif ny > 1:
			alpha, sx, sy, mirror, scale = get_params2D(IM[n])
			IM[n] = rot_shift2D(IM[n], alpha, sx, sy, mirror, scale)

	part = Dendo[K]
	AVE  = []
	ct   = 0
	for k in xrange(K):
		AVE.append(model_blank(nx, ny, nz))
		nobj = len(part[k])
		if nobj > 1:
			for id in part[k]: Util.add_img(AVE[k], IM[id])
			Util.mul_scalar(AVE[k], 1 / float(len(part[k])))
			AVE[k].set_attr('nobjects', len(part[k]))
			AVE[k].set_attr('members',  part[k])
			AVE[k].write_image(avename, k)
		else:
			im = IM[part[k][0]].copy()
			im.set_attr('nobjects', 1)
			im.set_attr('members', part[k][0])
			im.write_image(avename, k)


def tomo(box):
	print "box: ", box
	cmd = "your_executable -box " + box

	print "cmd: ", cmd
	from os import system
	status = system( cmd )
	print 'status: ', status

# Calculate averages to a given stack (wrap for ave_var in statistics)
def ave_ali(name_stack, name_out = None, ali = False, active = False):
	from statistics import ave_var, add_ave_varf, k_means_list_active
	from utilities  import file_type
	"""
	   This function is called by sxave_ali.py
	"""
	N = EMUtil.get_image_count(name_stack)
	if ali:	mode = 'a'
	else:	mode = ''
	if active: listID, N = k_means_list_active(name_stack)
	else:      listID    = range(N)
	data = EMData().read_images(name_stack, listID)
	
	ave, var = ave_var(data, mode)

	ext = file_type(name_stack)
	if name_out is None:
		if ext == 'bdb': name = name_stack.split(':')[1] + '.hdf'
		else:            name = name_stack
		ave.write_image('ave_' + name, 0)
	else:
		ave.write_image(name_out, 0)
#####################
#      ISAC CODE
#####################

def iter_isac(stack, ir, ou, rs, xr, yr, ts, maxit, CTF, snr, dst, FL, FH, FF, init_iter, main_iter, iter_reali, match_first, max_round, match_second, stab_ali, thld_err, indep_run, thld_grp, img_per_grp, generation):

	from filter 	  import filt_tanl
	from pixel_error  import multi_align_stability
	from statistics   import ave_series
	from utilities    import model_blank, write_text_file, set_params2D, get_params2D
	from utilities    import gather_EMData, bcast_EMData_to_all, bcast_list_to_all
	from utilities    import send_EMData, recv_EMData

	from mpi	  import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_send, mpi_recv, mpi_gatherv, mpi_scatterv
	from mpi	  import MPI_FLOAT, MPI_SUM, MPI_INT, MPI_TAG_UB
	from mpi	  import mpi_comm_split

	from random       import randint
	from sys          import exit
	from time         import localtime
	from subprocess   import call
	import os

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	if main_iter%iter_reali != 0:
		ERROR("main_iter should be a multiple of iter_reali, please reset them and restart and program", "iter_isac", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if generation == 0:
		ERROR("Generation should begin from 1, please reset it and restart and program", "iter_isac", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	ali_params_dir = "ali_params_generation_%d"%generation
	if os.path.exists(ali_params_dir):  
		ERROR('Output directory for alignment parameters exists, please either change its name or delete it and restart the program', "iter_isac", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		print "****************************************************************************************************"
		print "*                                                                                                  *"
		print "*                 Beginning of the ISAC program               "+strftime("%a, %d %b %Y %H:%M:%S", localtime())+"            *"
		print "*                                                                                                  *"
		print "* Iterative Stable Alignment and Clustering                                                        *"
		print "* By Zhengfan Yang, Jia Fang, Francisco Asturias and Pawel A. Penczek                   July 2011  *"
		print "*                                                                                                  *"
		print '* REFERENCE: Z. Yang, J. Fang, J. Chittuluru, F. J. Asturias and P. A. Penczek, "Iterative Stable  *'
		print '*            Alignment and Clustering of 2D Transmission Electron Microscope Images",              *' 
		print '*            Structure 20, 237-247, February 8, 2012.                                              *'
		print "*                                                                                                  *"
		print "* Last updated: 02/24/12                                                                           *"
		print "****************************************************************************************************"            
		print "*                                       Generation %3d                                             *"%(generation)
		print "****************************************************************************************************"
		print ""
	                                  

	color = myid%indep_run
	key = myid/indep_run
	group_comm = mpi_comm_split(MPI_COMM_WORLD, color, key)
	group_number_of_proc = mpi_comm_size(group_comm)
	group_main_node = 0

	# Read data on each processor, there are two ways, one is read on main_node and send them to all other nodes
	# The other way is all nodes reading it one by one, we have to test to determine which way is better.
	# The test shows that way 1 (18s) is way faster then way 2 (197s) on the test on 16 nodes.
	# The drawback of way 1 is it cannot have all attibutes, but I assume this is not important.
	# Method 1:
	if myid == main_node:
		alldata = EMData.read_images(stack)
		ndata = len(alldata)
		# alldata_n stores the original index of the particle (i.e., the index before running Generation 1)  
		alldata_n = [0]*ndata
		if generation > 1:
			for i in xrange(ndata): alldata_n[i] = alldata[i].get_attr('data_n')
		else:
			for i in xrange(ndata): alldata_n[i] = i
		nx = alldata[0].get_xsize()
	else:
		ndata = 0
		nx = 0
	ndata = mpi_bcast(ndata, 1, MPI_INT, main_node, MPI_COMM_WORLD)
	ndata = int(ndata[0])
	nx = mpi_bcast(nx, 1, MPI_INT, main_node, MPI_COMM_WORLD)
	nx = int(nx[0])

	if myid != main_node:
		alldata = [model_blank(nx, nx) for i in xrange(ndata)]
	mpi_barrier(MPI_COMM_WORLD)
	data = [None]*ndata
	tdummy = Transform({"type":"2D"})
	for im in xrange(ndata):
		bcast_EMData_to_all(alldata[im], myid, main_node)
		mpi_barrier(MPI_COMM_WORLD)  # has to be here, otherwise it chokes on our cluster.  PAP
		# This is the absolute ID, the only time we use it is
		# when setting the members of 4-way output. All other times, the id in 'members' is 
		# the relative ID.
		alldata[im].set_attr_dict({"xform.align2d": tdummy, "ID": im})
		data[im] = alldata[im]
	mpi_barrier(MPI_COMM_WORLD)
	'''
	# Method 2:
	for i in xrange(number_of_proc):
		if myid == i:
			print "Reading on node", i
			data = EMData.read_images(stack)
		mpi_barrier(MPI_COMM_WORLD)	
	'''

	ali_params_filename = "ali_params_%d"%color

	if myid == main_node:
		print "******************************************************************************************"
		print "*            Beginning of the first phase           "+strftime("%a, %d %b %Y %H:%M:%S", localtime())+"            *"
		print "*                                                                                        *"
		print "* The first phase is an exploratory phase. In this phase, we set the criteria to be very *"
		print "* loose and try to find as much candidate class averages as possible. This phase         *"
		print "* typically should have 10 to 20 rounds (default = 20). The candidate class averages are *"
		print "* stored in class_averages_candidate_generation_n.hdf.                                   *"
		print "******************************************************************************************"
		 
	avg_num = 0
	Iter = 1
	match_initialization = False
	avg_first_stage = "class_averages_candidate_generation_%d.hdf"%generation

	# I am adding here Artificial Intelligence for stopping 
	#  The program should stop if
	#    (a)  three times in a row it could not find new stable groups
	couldnt_find_stable = 0
	#    (b)  if number of groups to process is less than three
	K = ndata/img_per_grp
	while Iter <= max_round and couldnt_find_stable < 3 and K > 3:
		if myid == main_node: 
			print "################################################################################"
			print "#           Beginning of Round %2d           "%Iter+strftime("%a, %d %b %Y %H:%M:%S", localtime())+"          #"
			print "################################################################################"
			print "**********************************************************************"
			print "               Initialization of averages using EQ-mref               "
			print "**********************************************************************"
			print "     Number of current images = %d            Number of groups = %d"%(ndata, K)

		# Generate random averages for each group
		if key == group_main_node:
			refim = generate_random_averages(data, K, 9023)
			#for j in xrange(len(refim)):  refim[j].write_image("refim_%d.hdf"%color, j)
		else:
			refim = [model_blank(nx, nx) for i in xrange(K)]

		for i in xrange(K):
			bcast_EMData_to_all(refim[i], key, group_main_node, group_comm)

		# Generate inital averages
		#if myid == main_node: print "     Generating initial averages ",localtime()[:5]
		refi = isac_MPI(data, refim, maskfile=None, outname=None, ir=ir, ou=ou, rs=rs, xrng=xr, yrng=yr, step=ts, 
				maxit=maxit, isac_iter=init_iter, CTF=CTF, snr=snr, rand_seed=-1, color=color, comm=group_comm, 
				stability=False, FL=FL, FH=FH, FF=FF, dst=dst)
		del refim

		# gather the data on main node
		if match_initialization:
			if key == group_main_node:
				#print "Begin gathering ...", myid, len(refi)
				refi = gather_EMData(refi, indep_run, myid, main_node)
			if myid == main_node:
				# Match all averages in the initialization and select good ones
				#print "before matching, len = ", len(refi)
				current_refim = get_unique_averages(refi, indep_run)
				# If data_good is too few, add some random ones, otherwise, cut to K
				print " found data good = ", len(current_refim)
				if len(current_refim) > K:
					current_refim = current_refim[:K]
				elif len(current_refim) < K:
					defi = K - len(current_refim)
					for i in xrange(defi):
						current_refim.append(refi[randint(0, indep_run*K-1)].copy())
			else:
				current_refim = [model_blank(nx, nx) for i in xrange(K)]

			mpi_barrier(MPI_COMM_WORLD)	
		else:
			current_refim = refi

		# broadcast current_refim to all nodes
		for i in xrange(K):
			bcast_EMData_to_all(current_refim[i], myid, main_node)
		mpi_barrier(MPI_COMM_WORLD)

		#if key == group_main_node:
		#	for i in xrange(K):
		#		current_refim[i].write_image("init_group%d_round%d.hdf"%(color, Iter), i)

		# Run ISAC
		if myid == main_node:
			print "**********************************************************************"
			print "                 Run the main part of ISAC program                    "
			print "**********************************************************************"
	
		for mloop in xrange(1, match_first+1):
			if myid == main_node:
				print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
				print "                  Loop %3d for 2-way matching               "%mloop
				print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			refi = isac_MPI(data, current_refim, maskfile=None, outname=None, ir=ir, ou=ou, rs=rs, xrng=xr, yrng=yr, step=ts, 
					maxit=maxit, isac_iter=main_iter, CTF=CTF, snr=snr, rand_seed=-1, color=color, comm=group_comm, 
					stability=True, stab_ali=stab_ali, iter_reali=iter_reali, thld_err=thld_err, FL=FL, FH=FH, FF=FF, dst=dst)

			all_ali_params = [[] for i in xrange(5)]
			for im in data:
				alpha, sx, sy, mirror, scale = get_params2D(im)
				all_ali_params[0].append(alpha)
				all_ali_params[1].append(sx)			
				all_ali_params[2].append(sy)
				all_ali_params[3].append(mirror)
				all_ali_params[4].append(scale)
			if key == group_main_node:
				call(['rm', '-f', ali_params_filename])
				write_text_file(all_ali_params, ali_params_filename)
			del all_ali_params

			# gather the data from the group main node to the main node
			if key == group_main_node:
				refi = gather_EMData(refi, indep_run, myid, main_node)

			if mloop != match_first:
				if myid == main_node:
					current_refim = match_2_way(data, refi, indep_run, thld_grp, FH, FF)
				else:
					current_refim = [model_blank(nx, nx) for i in xrange(K)]
				for k in xrange(K):
					bcast_EMData_to_all(current_refim[k], myid, main_node)
			mpi_barrier(MPI_COMM_WORLD)
		del current_refim

		# Run Matching
		if myid == main_node:
			#print " Before matching ...  ", localtime()[:5] #len(data), len(refi), indep_run
			matched_data = match_2_way(data, refi, indep_run, thld_grp, FH, FF)
			members = []
			for im in matched_data:
				im.write_image(avg_first_stage, avg_num)
				avg_num += 1
				members.extend(im.get_attr('members'))

			# Because it's 2-way matching, it is possible some members are accounted for twice, we must delete the duplicate ones.   Yang 03/28/11
			members.sort()
			for i in xrange(len(members)-1, 0, -1):
				if members[i] == members[i-1]: del members[i]
			for i in xrange(len(members)-1): assert members[i]!=members[i+1]
			mem_len = len(members)
			print "In Round #%d, we found %d stable and reproducible averages, accounted for %d particles.  "%(Iter, len(matched_data), mem_len)
		else:
			mem_len = 0
			members = []
		mem_len = mpi_bcast(mem_len, 1, MPI_INT, main_node, MPI_COMM_WORLD)
		mem_len = int(mem_len[0])

		if mem_len > 0:
			# In members we have absolute ID
			members = mpi_bcast(members, mem_len, MPI_INT, main_node, MPI_COMM_WORLD)
			members = map(int, members)

			# Take out the good ones and use the remaining ones for initialization again
			nndata = ndata-len(members)
			newdata = [-1]*nndata
			ll = 0
			for i in xrange(ndata):
				abs_id = data[i].get_attr("ID")
				if abs_id not in members:
					newdata[ll] = abs_id
					ll += 1
			data = [alldata[im] for im in newdata]
			del newdata
			for im in data:
				im.set_attr("xform.align2d", tdummy)
			ndata = nndata

			couldnt_find_stable = 0
			K = ndata/img_per_grp
		else:
			couldnt_find_stable += 1
		Iter += 1
		mpi_barrier(MPI_COMM_WORLD)
	
	del data
	refim_stack = avg_first_stage

	if myid == main_node:
		print "******************************************************************************************"
		print "*              End of the first phase             "+strftime("%a, %d %b %Y %H:%M:%S", localtime())+"              *"
		print "******************************************************************************************"
		print ""
		print "******************************************************************************************"
		print "*           Beginning of the second phase           "+strftime("%a, %d %b %Y %H:%M:%S", localtime())+"            *"
		print "*                                                                                        *"
		print "* The second phase is where the actual class averages are generated, it typically have   *"
		print "* 3~9 iterations (default = 5) of matching. The first half of iterations are 2-way       *"
		print "* matching, the second half of iterations are 3-way matching, and the last iteration is  *"
		print "* 4-way matching. In the second phase, three files will be generated:                    *"
		print "* class_averages_generation_n.hdf : class averages generated in this generation          *"
		print "* generation_n_accounted.txt      : IDs of accounted particles in this generation        *"
		print "* generation_n_unaccounted.txt    : IDs of unaccounted particles in this generation      *"
		print "******************************************************************************************"
		try:
			refim = EMData.read_images(refim_stack)
		except:
			refim = []
		nrefim = len(refim)
	else:
		nrefim = 0
	nrefim = mpi_bcast(nrefim, 1, MPI_INT, main_node, MPI_COMM_WORLD)
	nrefim = int(nrefim[0])

	if myid != main_node:
		refim = [model_blank(nx, nx) for i in xrange(nrefim)]
	mpi_barrier(MPI_COMM_WORLD)

	nn = [0]*nrefim
	for i in xrange(nrefim):
		bcast_EMData_to_all(refim[i], myid, main_node)
		if myid == main_node: n_objects = refim[i].get_attr('n_objects')
		else: n_objects = 0
		n_objects = mpi_bcast(n_objects, 1, MPI_INT, main_node, MPI_COMM_WORLD)
		nn[i] = int(n_objects[0])
		if myid != main_node:  refim[i].set_attr('n_objects', int(n_objects[0]))
		mpi_barrier(MPI_COMM_WORLD)
	nn.sort()

	if len(nn) > 0: img_per_grp = nn[-1]
	refim_all = refim

	two_way_loop = match_second/2
	ndata = len(alldata)
	K = ndata/img_per_grp
	for mloop in xrange(1, match_second+1):
		if mloop <= two_way_loop:
			wayness = 2
			wayness_english = "two"
		elif mloop != match_second:
			wayness = 3
			wayness_english = "three"
		else:
			wayness = 4
			wayness_english = "four"
		if myid == main_node:		
			print "################################################################################"
			print "#       Iteration %2d for %d-way matching       "%(mloop, wayness)+strftime("%a, %d %b %Y %H:%M:%S", localtime())+"        #"
			print "################################################################################"

			members = []
			for im in refim_all:
				if im.get_attr('n_objects') > 1:
					members.extend(im.get_attr('members'))
			members.sort()
			for i in xrange(len(members)-1, 0, -1):
				if members[i] == members[i-1]: del members[i]
			n_members = len(members)
		else:
			n_members = 0
		n_members = mpi_bcast(n_members, 1, MPI_INT, main_node, MPI_COMM_WORLD)
		n_members = int(n_members[0])
		if myid != main_node:
			members = [0]*n_members
		members = mpi_bcast(members, n_members, MPI_INT, main_node, MPI_COMM_WORLD)
		members = map(int, members)

		ndata = len(alldata)
		nleft = ndata-n_members
		data_left = [None]*nleft
		c = 0
		for i in xrange(ndata):
			if i not in members:
				data_left[c] = alldata[i]
				c += 1

		K_left = nleft/img_per_grp
		if K_left > 0:
			if myid == main_node: 
				print "**********************************************************************"
		 		print "        Generating initial averages for unaccounted for images        "
				print "**********************************************************************"
				print "   Number of images unaccounted for = %d     Number of groups = %d"%(nleft, K_left)

			# Generate random averages for each group
			if key == group_main_node:
				refim_left = generate_random_averages(data_left, K_left)
				#for j in xrange(K_left):  refim_left[j].write_image("refim_left_%d.hdf"%color, j)
			else:
				refim_left = [model_blank(nx, nx) for i in xrange(K_left)]

			for i in xrange(K_left):
				bcast_EMData_to_all(refim_left[i], key, group_main_node, group_comm)

			# Generate initial averages for the unaccouted images
			refim_left = isac_MPI(data_left, refim_left, maskfile=None, outname=None, ir=ir, ou=ou, rs=rs, xrng=xr, yrng=yr, step=ts, 
					maxit=maxit, isac_iter=init_iter, CTF=CTF, snr=snr, rand_seed=-1, color=color, comm=group_comm, stability=False, 
					FL=FL, FH=FH, FF=FF, dst=dst)

			if len(refim) < K:
				# This will only happen in the first iteration, if applicable
				for k in xrange(K_left):
					refim.append(refim_left[k])
				refim = refim[:K]
			else:
				refim = refim[:K]
				ileft = 0
				for k in xrange(K):
					if refim[k].get_attr('n_objects') == 1:
						refim[k] = refim_left[ileft]
						ileft += 1
						if ileft >= K_left:  break
			mpi_barrier(MPI_COMM_WORLD)

			#if key == group_main_node:
			#	for i in xrange(K):
			#		refim[i].write_image("init_group%d_2nd_phase_round%d.hdf"%(color, mloop), i)

		# Run ISAC
		if myid == main_node:
			print "**********************************************************************"
			print "                 Run the main part of ISAC program                    "
			print "**********************************************************************"
			print "    Number of images = %d               Number of groups = %d"%(ndata, K)

		refim = isac_MPI(alldata, refim, maskfile=None, outname=None, ir=ir, ou=ou, rs=rs, xrng=xr, yrng=yr, step=ts, 
				maxit=maxit, isac_iter=main_iter, CTF=CTF, snr=snr, rand_seed=-1, color=color, comm=group_comm, 
				stability=True, stab_ali=stab_ali, iter_reali=iter_reali, thld_err=thld_err, FL=FL, FH=FH, FF=FF, dst=dst)

		all_ali_params = [[] for i in xrange(5)]
		for im in alldata:
			alpha, sx, sy, mirror, scale = get_params2D(im)
			all_ali_params[0].append(alpha)
			all_ali_params[1].append(sx)			
			all_ali_params[2].append(sy)
			all_ali_params[3].append(mirror)
			all_ali_params[4].append(scale)
		if key == group_main_node:
			call(['rm', '-f', ali_params_filename])
			write_text_file(all_ali_params, ali_params_filename)

		# gather refim to the main node
		if key == group_main_node:
			refim = gather_EMData(refim, indep_run, myid, main_node)

		if mloop != match_second:
			if myid == main_node:
 				print "**********************************************************************"
 				print "                   Run the %d-way matching algorithm                   "%wayness
 				print "**********************************************************************"
				# In this last two-way loop, we find all unique 2-way matches and use it as the starting
				# point of three-way match
				if mloop == two_way_loop:
					refim_all = match_2_way(alldata, refim, indep_run, thld_grp, FH, FF)
					# If they are enough, good; otherwise, add some random images into it.
					if len(refim_all) > K:
						print "Since it is more than the number of groups (%d), we only use the first %d of them."%(K, K)
						refim_all = refim_all[:K]
					elif len(refim_all) < K:
						defi = K - len(refim_all)
						print "Since it is less than the number of groups (%d), we have to append %d random images."%(K, defi)
						for i in xrange(defi):
							# put some dummy avgs here
							temp_id = randint(0, ndata-1)
							ave = alldata[temp_id].copy()
							ave.set_attr_dict({"members": [temp_id], "n_objects": 1})
							refim_all.append(ave)
					for i in xrange(K*(indep_run-1)):
						refim_all.append(refim[i%K])
				else:
					refim_all = match_2_way(alldata, refim, indep_run, thld_grp, FH, FF, find_unique = False, wayness = wayness)
			else:
				refim_all = [model_blank(nx, nx) for i in xrange(K*indep_run)]
			for k in xrange(K*indep_run):
				bcast_EMData_to_all(refim_all[k], myid, main_node)
				if myid == main_node: n_objects = refim_all[k].get_attr('n_objects')
				else: n_objects = 0
				n_objects = mpi_bcast(n_objects, 1, MPI_INT, main_node, MPI_COMM_WORLD)
				if myid != main_node:  refim_all[k].set_attr('n_objects', int(n_objects[0]))
			refim = refim_all[color*K:(color+1)*K]
			#if key == group_main_node:
			#	for k in xrange(K): refim[k].write_image("%s_way_match_%02d_%02d.hdf"%(wayness_english, mloop, color), k)
		mpi_barrier(MPI_COMM_WORLD)

	if key == group_main_node:
		call(['rm', '-f', ali_params_filename])

	if myid == main_node:
		print "**********************************************************************"
		print "                   Run the 4-way matching algorithm                   "
		print "**********************************************************************"
	
	# Run 4-way Matching
	# Comment by Zhengfan Yang on 6/20/11
	# The original design was way too slow, we have to send all stable sets to each node and let each node to do the realignment
	# and send them back, even though the code will be much more complicated.
	# I have decided that main node should not do realignment, otherwise it could clog the whole operation if it happened to have
	# a very large group.  The main node is used to send and collect information.
	
	if myid == main_node:
		STB_PART = match_independent_runs(alldata, refim, indep_run, thld_grp)
		l_STB = len(STB_PART)
		call(['mkdir', ali_params_dir])
	else:
		l_STB = 0
	mpi_barrier(MPI_COMM_WORLD)
	l_STB = mpi_bcast(l_STB, 1, MPI_INT, main_node, MPI_COMM_WORLD)
	l_STB = int(l_STB[0])

	if myid == main_node:
		for i in xrange(l_STB):
			node_to_run = i%(number_of_proc-1)+1
			mpi_send(len(STB_PART[i]), 1, MPI_INT, node_to_run, i+10000, MPI_COMM_WORLD)
			mpi_send(STB_PART[i], len(STB_PART[i]), MPI_INT, node_to_run, i+20000, MPI_COMM_WORLD)

		members_acc = []
		ave_num = 0	
		for i in xrange(l_STB):
			node_to_run = i%(number_of_proc-1)+1
			l_stable_members = mpi_recv(1, MPI_INT, node_to_run, i+30000, MPI_COMM_WORLD)
			l_stable_members = int(l_stable_members[0])
			stable_members = mpi_recv(l_stable_members, MPI_INT, node_to_run, i+40000, MPI_COMM_WORLD)
			stable_members = map(int, stable_members)
			mirror_consistent_rate = mpi_recv(1, MPI_FLOAT, node_to_run, i+50000, MPI_COMM_WORLD)
			mirror_consistend_rate = float(mirror_consistent_rate[0])
			pix_err = mpi_recv(1, MPI_FLOAT, node_to_run, i+60000, MPI_COMM_WORLD)
			pix_err = float(pix_err[0])
	
			print  "Group %d ...... Mirror consistent rate = %f"%(i, mirror_consistent_rate)
			print  "Group %d ...... Average pixel error = %f"%(i, pix_err)
			print  "Group %d ...... Size of stable subset = %d"%(i, l_stable_members)
			print  "Group %d ......"%i,
		
			if l_stable_members <= thld_grp:
				print "Size of stable subset smaller than the threshold, discarded\n"
				continue
			print "Size of stable subset larger than the threshold, kept\n"

			ave = recv_EMData(node_to_run, i+70000)
			stable_members_ori = [0]*l_stable_members
			for j in xrange(l_stable_members): stable_members_ori[j] = alldata_n[stable_members[j]]
			ave.set_attr_dict({"members": stable_members_ori, "n_objects": l_stable_members})
			ave.write_image("class_averages_generation_%d.hdf"%generation, ave_num)
			mpi_send(ave_num, 1, MPI_INT, node_to_run, i+80000, MPI_COMM_WORLD)
			ave_num += 1
			members_acc.extend(stable_members_ori)
			
		members_acc.sort()
		for i in xrange(len(members_acc)-1): assert members_acc[i] != members_acc[i+1]

		members_unacc = [0]*(ndata-len(members_acc))
		c = 0
		for i in xrange(ndata):
			if alldata_n[i] in members_acc: continue
			members_unacc[c] = alldata_n[i]
			c += 1
		print "In the second phase, we found %d stable and reproducible averages that account for %d particles.  "%(ave_num, len(members_acc))
		if members_acc != []:
			write_text_file(members_acc, "generation_%d_accounted.txt"%generation)
		if members_unacc != []:
			write_text_file(members_unacc, "generation_%d_unaccounted.txt"%generation)
		print "******************************************************************************************"
		print "*             End of the second phase             "+strftime("%a, %d %b %Y %H:%M:%S", localtime())+"              *"
		print "******************************************************************************************"
	else:
		for i in xrange(l_STB):
			node_to_run = i%(number_of_proc-1)+1
			if myid != node_to_run: continue
			l_STB_PART = mpi_recv(1, MPI_INT, main_node, i+10000, MPI_COMM_WORLD)
			l_STB_PART = int(l_STB_PART[0])
			STB_PART = mpi_recv(l_STB_PART, MPI_INT, main_node, i+20000, MPI_COMM_WORLD)
			STB_PART = map(int, STB_PART)
			STB_PART.sort()

			class_data = [None]*l_STB_PART
			members_id = [0]*l_STB_PART
 		        for im in xrange(l_STB_PART):
 	                        class_data[im] = alldata[STB_PART[im]]
 				members_id[im] = alldata[STB_PART[im]].get_attr('ID')
 			for im in xrange(l_STB_PART-1):
 				assert members_id[im] != members_id[im+1]

			ali_params = [[] for j in xrange(stab_ali)]
			for ii in xrange(stab_ali):
				ave = within_group_refinement(class_data, None, True, ir, ou, rs, [xr], [yr], [ts], dst, maxit, FH, FF)
				for im in xrange(l_STB_PART):
					alpha, sx, sy, mirror, scale = get_params2D(class_data[im])
					ali_params[ii].extend([alpha, sx, sy, mirror])
			if ou == -1:  ou = nx/2-2
			stable_set, mirror_consistent_rate, pix_err = multi_align_stability(ali_params, 0.0, 10000.0, thld_err, False, ou*2)
			
			l_stable_set = len(stable_set)
			stable_set_id = [0]*l_stable_set
			all_alpha = [0]*l_stable_set
			all_sx = [0]*l_stable_set
			all_sy = [0]*l_stable_set
			all_mirror = [0]*l_stable_set
			all_scale = [1.0]*l_stable_set
			for j in xrange(l_stable_set): 
				stable_set_id[j] = members_id[stable_set[j][1]]
				all_alpha[j] = ali_params[stab_ali-1][stable_set[j][1]*4]
				all_sx[j] = ali_params[stab_ali-1][stable_set[j][1]*4+1]
				all_sy[j] = ali_params[stab_ali-1][stable_set[j][1]*4+2]
				all_mirror[j] = ali_params[stab_ali-1][stable_set[j][1]*4+3]

			mpi_send(l_stable_set, 1, MPI_INT, main_node, i+30000, MPI_COMM_WORLD)
			mpi_send(stable_set_id, l_stable_set, MPI_INT, main_node, i+40000, MPI_COMM_WORLD)
			mpi_send(mirror_consistent_rate, 1, MPI_FLOAT, main_node, i+50000, MPI_COMM_WORLD)
			mpi_send(pix_err, 1, MPI_FLOAT, main_node, i+60000, MPI_COMM_WORLD)
			
			if l_stable_set > thld_grp:
				send_EMData(ave, main_node, i+70000)		
				ave_num = mpi_recv(1, MPI_INT, main_node, i+80000, MPI_COMM_WORLD)
				ave_num = int(ave_num[0])
				write_text_file([all_alpha, all_sx, all_sy, all_mirror, all_scale], "%s/ali_params_%03d"%(ali_params_dir, ave_num))
	
	mpi_barrier(MPI_COMM_WORLD)
	
	if myid == main_node:
		print "****************************************************************************************************"
		print "*                                                                                                  *"
		print "*                   End of the ISAC program                 "+strftime("%a, %d %b %Y %H:%M:%S", localtime())+"              *"
		print "*                                                                                                  *"
		print "****************************************************************************************************"            


def isac_MPI(stack, refim, maskfile = None, outname = "avim", ir=1, ou=-1, rs=1, xrng=0, yrng=0, step=1, 
            maxit=30, isac_iter=10, CTF=False, snr=1.0, rand_seed=-1, color=0, comm=-1, 
	    stability=False, stab_ali=5, iter_reali=1, thld_err=1.732, FL=0.1, FH=0.3, FF=0.2, dst=90.0):

	from alignment      import   Numrinit, ringwe, Applyws
	from applications   import   MPI_start_end
	from filter         import   filt_tanl, filt_params
	from fundamentals   import   rot_shift2D, fshift
	from morphology     import   ctf_2
	from pixel_error    import   multi_align_stability
	from statistics     import   ave_series
	from utilities      import   model_circle, combine_params2, inverse_transform2, get_image, get_im
	from utilities      import   reduce_EMData_to_root, bcast_EMData_to_all
	from utilities      import   get_params2D, set_params2D

	import os
	import sys
	from random         import   seed, randint, jumpahead
	from subprocess     import   call
	from mpi 	    import   mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD, mpi_barrier
	from mpi 	    import   mpi_reduce, mpi_bcast, mpi_barrier, mpi_recv, mpi_send
	from mpi 	    import   MPI_SUM, MPI_FLOAT, MPI_INT, MPI_TAG_UB

	if comm == -1: comm = MPI_COMM_WORLD		

	number_of_proc = mpi_comm_size(comm)
	myid = mpi_comm_rank(comm)
	main_node = 0

	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(isac_iter)

	if type(stack) == type(""):
		#read all data
		alldata = EMData.read_images(stack)
	else:
		alldata = stack
	nx = alldata[0].get_xsize()
	# default value for the last ring
	if last_ring == -1: last_ring = nx/2-2

	nima = len(alldata)
	
	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)

	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  mask = get_image(maskfile)
		else: mask = maskfile
	else : mask = model_circle(last_ring, nx, nx)
	if type(refim) == type(""):
		refi = EMData.read_images(refim)
	else:
		# It's safer to make a hard copy here. Although I am not sure, I believe a shallow copy
		# has messed up the program.
		# refi = refim
		refi = [None for i in xrange(len(refim))]
		for i in xrange(len(refim)):  refi[i] = refim[i].copy()
	numref = len(refi)

	#  CTF stuff
	if CTF:
		ctf_params = ima.get_attr("ctf")
		data_had_ctf = ima.get_attr("ctf_applied")
		ctm = ctf_2(nx, ctf_params)
		lctf = len(ctm)

	# IMAGES ARE SQUARES! center is in SPIDER convention
	cnx = nx/2+1
	cny = cnx

	mode = "F"
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr = ringwe(numr, mode)
	# reference images
	params = []
	#  for each node read its share of data
	#data = EMData.read_images(stack, range(image_start, image_end))
	for im in xrange(image_start, image_end):
		#data[im-image_start].set_attr('ID', im)
		if CTF:
			ctf_params = alldata[im].get_attr( "ctf" )
			if alldata[im].get_attr("ctf_applied") == 0:
				st = Util.infomask(alldata[im], mask, False)
				alldata[im] -= st[0]
				from filter import filt_ctf
				alldata[im] = filt_ctf(alldata[im], ctf_params)
				alldata[im].set_attr('ctf_applied', 1)

	if rand_seed > -1: seed(rand_seed)
	else:  seed()
	if myid != main_node: jumpahead(17*myid + 12345)

	fl = FL
	Iter = -1
	main_iter = 0

	while main_iter < max_iter:
		Iter += 1
		#if myid == main_node: print "Iter = ", Iter, "    main_iter = ", main_iter, "    len data = ", image_end-image_start, localtime()[0:5]
		ringref = []
		for j in xrange(numref):
			refi[j].process_inplace("normalize.mask", {"mask":mask, "no_sigma":1}) # normalize reference images to N(0,1)
			cimage = Util.Polar2Dm(refi[j] , cnx, cny, numr, mode)
			Util.Frngs(cimage, numr)
			Applyws(cimage, numr, wr)
			ringref.append(cimage)
		if CTF: ctf2 = [[[0.0]*lctf for k in xrange(2)] for j in xrange(numref)]
		peak_list = [[] for i in xrange(numref)]
		#  nima is the total number of images, not the one on this node, tha latter is (image_end-image_start)
		d = [0.0]*(numref*nima)
		# begin MPI section
		for im in xrange(image_start, image_end):
			alpha, sx, sy, mirror, scale = get_params2D(alldata[im])
			alphai, sxi, syi, scalei = inverse_transform2(alpha, sx, sy)
			# normalize
			alldata[im].process_inplace("normalize.mask", {"mask":mask, "no_sigma":0}) # subtract average under the mask

			# align current image to the reference
			temp = Util.multiref_polar_ali_2d_peaklist(alldata[im],
				ringref, xrng, yrng, step, mode, numr, cnx+sxi, cny+syi)
			for iref in xrange(numref):
				[alphan, sxn, syn, mn] = \
				   combine_params2(0.0, -sxi, -syi, 0, temp[iref*5+1], temp[iref*5+2], temp[iref*5+3], int(temp[iref*5+4]))
				peak_list[iref].extend([alphan, sxn, syn, mn])
				d[iref*nima+im] = temp[iref*5]
		del ringref
		del temp

		d = mpi_reduce(d, numref*nima, MPI_FLOAT, MPI_SUM, main_node, comm)
		if myid == main_node:
			d = map(float, d)
			id_list_long = Util.assign_groups(d, numref, nima)
			id_list = [[] for i in xrange(numref)]
			maxasi = nima/numref
			for i in xrange(maxasi*numref):
				id_list[i/maxasi].append(id_list_long[i])
			for i in xrange(nima%maxasi):
				id_list[id_list_long[-1]].append(id_list_long[maxasi*numref+i])
			for iref in xrange(numref):
				id_list[iref].sort()

			belongsto = [0]*nima
			for iref in xrange(numref):
				for im in id_list[iref]: belongsto[im] = iref
		else:
			belongsto = [0]*nima
		del d
		mpi_barrier(comm)
		belongsto = mpi_bcast(belongsto, nima, MPI_INT, main_node, comm)
		belongsto = map(int, belongsto)

		#  Compute partial averages
		members = [0]*numref
		sx_sum = [0.0]*numref
		sy_sum = [0.0]*numref
		for j in xrange(numref):  refi[j].to_zero()
		for im in xrange(image_start, image_end):
			matchref = belongsto[im]
			alphan = peak_list[matchref][(im-image_start)*4+0]
			sxn = peak_list[matchref][(im-image_start)*4+1]
			syn = peak_list[matchref][(im-image_start)*4+2]
			mn = peak_list[matchref][(im-image_start)*4+3]
			if mn == 0: sx_sum[matchref] += sxn
			else:       sx_sum[matchref] -= sxn
			sy_sum[matchref] += syn
			# apply current parameters and add to the average
			Util.add_img(refi[matchref], rot_shift2D(alldata[im], alphan, sxn, syn, mn))
			if CTF:
				ctm = ctf_2(nx, ctf_params)
				for i in xrange(lctf):  ctf2[matchref][it][i] += ctm[i]
			members[matchref] += 1
		sx_sum = mpi_reduce(sx_sum, numref, MPI_FLOAT, MPI_SUM, main_node, comm)
		sy_sum = mpi_reduce(sy_sum, numref, MPI_FLOAT, MPI_SUM, main_node, comm)
		members = mpi_reduce(members, numref, MPI_INT, MPI_SUM, main_node, comm)
		if myid != main_node:
			sx_sum = [0.0]*numref
			sy_sum = [0.0]*numref
			members = [0.0]*numref
		sx_sum = mpi_bcast(sx_sum, numref, MPI_FLOAT, main_node, comm)
		sy_sum = mpi_bcast(sy_sum, numref, MPI_FLOAT, main_node, comm)
		members = mpi_bcast(members, numref, MPI_INT, main_node, comm)
		sx_sum = map(float, sx_sum)
		sy_sum = map(float, sy_sum)
		members = map(int, members)

		for j in xrange(numref):
			sx_sum[j] /= float(members[j])
			sy_sum[j] /= float(members[j])

		for im in xrange(image_start, image_end):
			matchref = belongsto[im]
			alphan = peak_list[matchref][(im-image_start)*4+0]
			sxn = peak_list[matchref][(im-image_start)*4+1]
			syn = peak_list[matchref][(im-image_start)*4+2]
			mn = int(peak_list[matchref][(im-image_start)*4+3])
			if mn == 0:
				set_params2D(alldata[im], [alphan, sxn-sx_sum[matchref], syn-sy_sum[matchref], mn, scale])
			else:
				set_params2D(alldata[im], [alphan, sxn+sx_sum[matchref], syn-sy_sum[matchref], mn, scale])

		for j in xrange(numref):
			reduce_EMData_to_root(refi[j], myid, main_node, comm)
			if myid == main_node:
				# Golden rule when to do within group refinement
				Util.mul_scalar(refi[j], 1.0/float(members[j]))
				refi[j] = filt_tanl(refi[j], fl, FF)
				refi[j] = fshift(refi[j], -sx_sum[j], -sy_sum[j])
				set_params2D(refi[j], [0.0, 0.0, 0.0, 0, 1.0])	
		
		if myid == main_node:
			dummy = within_group_refinement(refi, mask, True, first_ring, last_ring, rstep, [xrng], [yrng], [step], dst, maxit, FH, FF)
			ref_ali_params = []
			for j in xrange(numref):
				alpha, sx, sy, mirror, scale = get_params2D(refi[j])
				refi[j] = rot_shift2D(refi[j], alpha, sx, sy, mirror)
				ref_ali_params.extend([alpha, sx, sy, mirror])
		else:
			ref_ali_params = [0.0]*(numref*4)
		ref_ali_params = mpi_bcast(ref_ali_params, numref*4, MPI_FLOAT, main_node, comm)
		ref_ali_params = map(float, ref_ali_params)

		for j in xrange(numref):
			bcast_EMData_to_all(refi[j], myid, main_node, comm)

		# Compensate the centering to averages
		for im in xrange(image_start, image_end):
			matchref = belongsto[im]
			alpha, sx, sy, mirror, scale = get_params2D(alldata[im])
			alphan, sxn, syn, mirrorn = combine_params2(alpha, sx, sy, mirror, ref_ali_params[matchref*4], ref_ali_params[matchref*4+1], \
				ref_ali_params[matchref*4+2], int(ref_ali_params[matchref*4+3]))
			set_params2D(alldata[im], [alphan, sxn, syn, int(mirrorn), 1.0])

		do_within_group = 0
		fl += 0.05
		if fl >= FH:
			fl = FL
			do_within_group = 1
			
		if do_within_group == 1:
			#if myid == main_node: print "Doing within group alignment .......", localtime()[0:5]

			# Broadcast the alignment parameters to all nodes
			for i in xrange(number_of_proc):
				im_start, im_end = MPI_start_end(nima, number_of_proc, i)
				if myid == i:
					ali_params = []
					for im in xrange(image_start, image_end):
						alpha, sx, sy, mirror, scale = get_params2D(alldata[im])
						ali_params.extend([alpha, sx, sy, mirror])
				else:
					ali_params = [0.0]*((im_end-im_start)*4)
				ali_params = mpi_bcast(ali_params, len(ali_params), MPI_FLOAT, i, comm)
				ali_params = map(float, ali_params)
				for im in xrange(im_start, im_end):
					alpha = ali_params[(im-im_start)*4]
					sx = ali_params[(im-im_start)*4+1]
					sy = ali_params[(im-im_start)*4+2]
					mirror = int(ali_params[(im-im_start)*4+3])
					set_params2D(alldata[im], [alpha, sx, sy, mirror, 1.0])
	
			# Here the assumption is that numref > number_of_proc, should be checked.
			main_iter += 1
			for j in xrange(myid, numref, number_of_proc):
				assign = []
				for im in xrange(nima):
					if j == belongsto[im]:  assign.append(im)

				randomize = True  # I think there is no reason not to be True
				class_data = [alldata[im] for im in assign]
				refi[j] = within_group_refinement(class_data, mask, randomize, first_ring, last_ring, rstep, [xrng], [yrng], [step], dst, maxit, FH, FF)

				# Here stability does not need to be checked for each main iteration, it only need to
				# be done for every 'iter_reali' iterations. If one really want it to be checked each time
				# simple set iter_reali to 1, which is the default value right now.
				if stability and main_iter%iter_reali==0:
					ali_params = [[] for qq in xrange(stab_ali)]
					for ii in xrange(stab_ali):
						if ii > 0:  # The first one does not have to be repeated
							dummy = within_group_refinement(class_data, mask, randomize, first_ring, last_ring, rstep, [xrng], [yrng], [step], dst, maxit, FH, FF)
						for im in xrange(len(class_data)):
							alpha, sx, sy, mirror, scale = get_params2D(class_data[im])
							ali_params[ii].extend([alpha, sx, sy, mirror])

					stable_set, mirror_consistent_rate, err = multi_align_stability(ali_params, 0.0, 10000.0, thld_err, False, last_ring*2)

					#print  "Class %d ...... Size of stable subset = %d"%(j, len(stable_set))
					#print  "Class %d ...... Mirror consistent rate = %f"%(j, mirror_consistent_rate)
					#print  "Class %d ...... Average error = %f"%(j, err)

					# If the size of stable subset is too small (say 1, 2), it will cause many problems, so we manually increase it to 5
					while len(stable_set) < 5:
						duplicate = True
						while duplicate:
							duplicate = False
							p = randint(0, len(class_data)-1)
							for ss in stable_set:
								if p == ss[1]: duplicate = True
						stable_set.append([100.0, p])

					stable_data = []
					stable_members = []
					for err in stable_set:
						im = err[1]
						stable_members.append(assign[im])
						stable_data.append(class_data[im])
					stable_members.sort()

					refi[j] = filt_tanl(ave_series(stable_data), FH, FF)
					refi[j].set_attr('members', stable_members)
					refi[j].set_attr('n_objects', len(stable_members))
					del stable_members
				# end of stability
				del assign
			mpi_barrier(comm)

			for j in xrange(numref):
				bcast_EMData_to_all(refi[j], myid, j%number_of_proc, comm)

			for im in xrange(nima):
				done_on_node = belongsto[im]%number_of_proc
				if myid == done_on_node:
					alpha, sx, sy, mirror, scale = get_params2D(alldata[im])
					ali_params = [alpha, sx, sy, mirror]
				else:
					ali_params = [0.0]*4
				ali_params = mpi_bcast(ali_params, 4, MPI_FLOAT, done_on_node, comm)
				ali_params = map(float, ali_params)
				set_params2D(alldata[im], [ali_params[0], ali_params[1], ali_params[2], int(ali_params[3]), 1.0])

			if stability and main_iter%iter_reali == 0:
				# In this case, we need to set the 'members' attr using stable members from the stability test
				for j in xrange(numref):
					done_on_node = j%number_of_proc
					if done_on_node != main_node:
						if myid == main_node:
							mem_len = mpi_recv(1, MPI_INT, done_on_node, MPI_TAG_UB, comm)
							mem_len = int(mem_len[0])
							members = mpi_recv(mem_len, MPI_INT, done_on_node, MPI_TAG_UB, comm)
							members = map(int, members)
							refi[j].set_attr_dict({'members': members,'n_objects': mem_len})
						elif myid == done_on_node:
							members = refi[j].get_attr('members')
							mpi_send(len(members), 1, MPI_INT, main_node, MPI_TAG_UB, comm)
							mpi_send(members, len(members), MPI_INT, main_node, MPI_TAG_UB, comm)
			#if myid == main_node:  print "within group alignment done. ", localtime()[0:5]
		# end of do_within_group
		mpi_barrier(comm)

		if myid == main_node:
			#  I added a switch here.  I do not think we will need those in the future.  PAP 03/26
			if outname != None:
				call(['rm', '-f', outname+'%02d_%03d.hdf'%(color, Iter)])
			if stability and main_iter%iter_reali == 0:
				# In this case, the attr 'members' is defined as the stable members, its setting is done
				# in the code before
				if outname != None:
					for j in xrange(numref):
						refi[j].write_image(outname+"%02d_%03d.hdf"%(color, Iter), j)
			else:
				for j in xrange(numref):
					refi[j].set_attr_dict({'members': id_list[j], 'n_objects': len(id_list[j])})
					if outname != None:
						refi[j].write_image(outname+"%02d_%03d.hdf"%(color, Iter), j)
		mpi_barrier(comm)

	return refi


def within_group_refinement(data, maskfile, randomize, ir, ou, rs, xrng, yrng, step, dst, maxit, FH, FF):

	# Comment by Zhengfan Yang 03/11/11
	# This is a simple version of ali2d_data (down to the bone), no output dir, no logfile, no CTF, no MPI or CUDA, no Fourvar, no auto stop, no user function
	
	from alignment    import Numrinit, ringwe, ali2d_single_iter
	from filter	  import filt_tanl
	from fundamentals import fshift
	from random	  import randint, random
	from statistics   import ave_series
	from utilities    import get_input_from_string, model_circle, set_params2D, get_params2D, combine_params2, inverse_transform2

	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit);
	nima = len(data)
	nx = data[0].get_xsize()
	if last_ring == -1:  last_ring = nx/2-2
	if maskfile: mask = maskfile
	else: mask = model_circle(last_ring, nx, nx)

	if randomize:
		for im in data:
			alpha, sx, sy, mirror, scale = get_params2D(im)
			alphai, sxi, syi, mirrori = inverse_transform2(alpha, sx, sy)
			alphan, sxn, syn, mirrorn = combine_params2(0.0, -sxi, -syi, 0, random()*360.0, 0.0, 0.0, randint(0, 1))
			#set_params2D(im, [random()*360.0, 0.0, 0.0, randint(0, 1), 1.0])
			set_params2D(im, [alphan, sxn, syn, mirrorn, 1.0])

	cnx = nx/2+1
 	cny = cnx
 	mode = "F"
	numr = Numrinit(first_ring, last_ring, rstep, mode)
 	wr = ringwe(numr, mode)

	sx_sum = 0.0
	sy_sum = 0.0
	cs = [0.0]*2
	total_iter = 0
	for N_step in xrange(len(xrng)):
		for Iter in xrange(max_iter):
			total_iter += 1
			tavg = ave_series(data)
		        fl = 0.1+(FH-0.1)*Iter/float(max_iter-1)
			tavg = filt_tanl(tavg, fl, FF)
			if total_iter == len(xrng)*max_iter:  return tavg
			cs[0] = sx_sum/float(nima)
			cs[1] = sy_sum/float(nima)
			tavg = fshift(tavg, -cs[0], -cs[1])
			if Iter%4 != 0 or total_iter > max_iter*len(xrng)-10: delta = 0.0
			else: delta = dst
			sx_sum, sy_sum = ali2d_single_iter(data, numr, wr, cs, tavg, cnx, cny, xrng[N_step], yrng[N_step], step[N_step], mode, CTF=False, delta=delta)


def match_independent_runs(data, refi, n_group, T):

	from numpy       import array
	from statistics  import k_means_stab_bbenum
	from time        import localtime

	K = len(refi)/n_group
	Parts = []
	for i in xrange(n_group):
                part = [] 
                for k in xrange(K): 
                        lid = refi[i*K+k].get_attr('members') 
                        lid = array(lid, 'int32') 
                        lid.sort() 
                        part.append(lid)
		Parts.append(part)

	#print Parts
	#print "Before matching = ", localtime()[:5]
        MATCH, STB_PART, CT_s, CT_t, ST, st = k_means_stab_bbenum(Parts, T=T, J=50, max_branching=40, stmult=0.1, branchfunc=2)
	#print "After matching = ", localtime()[:5]

        # I commented out next three, not much use printing them,  PAP.
	#print MATCH
        #print STB_PART
        #print CT_s
        #print CT_t
        #print ST
 	#print st

        cost_by_match_thresh = []
        for i in xrange(len(CT_s)):
                if CT_s[i] > T:
                        cost_by_match_thresh.append(CT_s[i])

        print "%d-way match: total cost of matches over threshold: "%(len(Parts)), sum(cost_by_match_thresh)
        print "%d-way match: total number of matches over threshold: "%(len(Parts)), len(cost_by_match_thresh)
        print "%d-way match: cost by match over threshold: "%(len(Parts)), cost_by_match_thresh
        print " "

	STB_PART_cleaned = []
	for i in xrange(len(STB_PART)):
		if len(STB_PART[i]) > T:
			STB_PART_cleaned.append(STB_PART[i])
	
	return STB_PART_cleaned


def match_2_way(data, refi, indep_run, thld_grp, FH, FF, find_unique=True, wayness=2):

	from utilities import read_text_row, set_params2D

	from statistics import ave_series, k_means_stab_bbenum
	from random     import randint, shuffle
	from filter 	import filt_tanl
	from numpy      import array

	K = len(refi)/indep_run
	run = range(indep_run)
	shuffle(run)

	reproducible_avgs = []
        
	for irun in xrange(indep_run):
		filename = "ali_params_%d"%run[irun]
		all_ali_params = read_text_row(filename)
	
		Parts = []
	        part = [] 
        	for k in xrange(K): 
                	lid = refi[run[irun]*K+k].get_attr('members') 
                        lid = array(lid, 'int32') 
                        lid.sort() 
                        part.append(lid)
		Parts.append(part)

	        part = [] 
        	for k in xrange(K): 
                	lid = refi[run[(irun+1)%indep_run]*K+k].get_attr('members') 
                        lid = array(lid, 'int32') 
                        lid.sort() 
                        part.append(lid)
		Parts.append(part)

		if wayness == 3:
		        part = [] 
	        	for k in xrange(K): 
        	        	lid = refi[run[(irun+2)%indep_run]*K+k].get_attr('members') 
                	        lid = array(lid, 'int32') 
                        	lid.sort() 
	                        part.append(lid)
			Parts.append(part)

	        MATCH, STB_PART, CT_s, CT_t, ST, st = k_means_stab_bbenum(Parts, T=thld_grp, J=50, max_branching=40, stmult=0.1, branchfunc=2)

       		cost_by_match_thresh = []
	        for i in xrange(len(CT_s)):
        	        if CT_s[i] > thld_grp:
                	        cost_by_match_thresh.append(CT_s[i])

	        print "%d-way match: total cost of matches over threshold: "%(len(Parts)), sum(cost_by_match_thresh)
	        print "%d-way match: total number of matches over threshold: "%(len(Parts)), len(cost_by_match_thresh)
        	print "%d-way match: cost by match over threshold: "%(len(Parts)), cost_by_match_thresh
	        print " "

		for i in xrange(len(STB_PART)):
                	if len(STB_PART[i]) > 0:
                        	class_data = []
				members_id = []
                        	for im in STB_PART[i]:
					set_params2D(data[im], [all_ali_params[im][0], all_ali_params[im][1], all_ali_params[im][2], int(all_ali_params[im][3]), 1.0])
                                	class_data.append(data[im])
					members_id.append(data[im].get_attr('ID'))
        			ave = ave_series(class_data)
                		ave.set_attr_dict({"members": members_id, "n_objects": len(members_id)})
      	        		ave = filt_tanl(ave, FH, FF)
        	        	reproducible_avgs.append(ave)
			else:
				# put some dummy avgs here
				temp_id = randint(0, len(data)-1)
				ave = data[temp_id].copy()
	                	ave.set_attr_dict({"members": [temp_id], "n_objects": 1})
				reproducible_avgs.append(ave)

	if find_unique:
		# Here the idea is like this: in all reproducible averages, find the ones that are unique
		reproducible_avgs_unique = get_unique_averages(reproducible_avgs, indep_run)
		print "Found %d unique class averages through %d-way matching"%(len(reproducible_avgs_unique), wayness)
		return reproducible_avgs_unique
	else:
		return reproducible_avgs


def generate_random_averages(data, K, rand_seed = -1):

	from statistics import aveq
	from random import randint, shuffle, seed
	#  I prefer to take random images....  PAP
	if rand_seed == -1: seed()
	else: seed(rand_seed)
	ndata = len(data)
	ll = range(ndata)
	shuffle(ll)
	avgs = [data[ll[i]].copy() for i in xrange(K)]
	return avgs

	'''
	from alignment import align2d
	from fundamentals import rot_shift2D
	for im in xrange(K,ndata):
		wnmr = im%K
		alpha,sx,sy,mirror,peak = align2d(data[ll[im]], avgs[wnmr], 1,1,0.5,1,30)
		Util.add_img(avgs[wnmr], rot_shift2D(data[ll[im]], alpha, sx, sy, mirror))

	l = [[] for i in xrange(K)]
	for im in xrange(ndata):
		l[randint(0, K-1)].append(im)

	avgs = []
	for k in xrange(K):
		if len(l[k]) > 1:
			temp_stack = [None]*len(l[k])
			for i in xrange(len(l[k])):
				temp_stack[i] = data[l[k][i]]
			ave = aveq(temp_stack, mode = "")
		elif len(l[k]) == 1: ave = data[l[k][0]]
		else: ave = data[randint(0, ndata)]
		avgs.append(ave)
	'''
	return avgs


def get_unique_averages(data, indep_run, m_th=0.45):
	
	size_all = len(data)
	size = size_all/indep_run
	assert size_all%indep_run == 0	

	# Meaning of flag
	# 0 - not matched yet
	# 1 - matched, kept
	# 2 - matched, discard, as we only want to keep a unique copy
	# 3 - not matched, kept
	flag = [0]*size_all
	
	for i in xrange(size_all-size):
		m1 = map(int, data[i].get_attr('members'))
		if len(m1) < 5: continue
		for j in xrange((i/size+1)*size, size_all):
			m2 = map(int, data[j].get_attr('members'))
			if len(m2) < 5: continue
			m = set(m1).intersection(set(m2))
			if float(len(m)) > m_th*(len(m1)+len(m2))/2:
				if flag[i] == 0 and flag[j] == 0:
					flag[i] = 1
					flag[j] = 2
					break
				elif flag[i] == 2 or flag[j] == 2:
					flag[i] = 2
					flag[j] = 2
				else:
					print "Impossible: Something is wrong!"
		if flag[i] == 0: flag[i] = 3

	data_good = []
	# Give priority to the averages that appeared more than once
	for im in xrange(size_all):
		if flag[im] == 1: 	data_good.append(data[im])
	for im in xrange(size_all):
		if flag[im] == 3: 	data_good.append(data[im])
	
	return data_good

