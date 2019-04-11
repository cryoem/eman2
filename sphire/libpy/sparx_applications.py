#
from __future__ import print_function
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

import EMAN2_cppwrap
import EMAN2db
import sparx_alignment
import copy
import sparx_filter
import sparx_fundamentals
import sparx_global_def
import sparx_logger
import math
import sparx_morphology
import mpi
import sparx_multi_shc
import numpy
import numpy.random
import os
import pickle
import sparx_pixel_error
import sparx_projection
import random
import sparx_reconstruction
import scipy.optimize
import sparx_statistics
import string
import string as sting
import sys
import time
import sparx_user_functions
import sparx_utilities
pass#IMPORTIMPORTIMPORT import EMAN2
pass#IMPORTIMPORTIMPORT import EMAN2_cppwrap
pass#IMPORTIMPORTIMPORT import EMAN2db
pass#IMPORTIMPORTIMPORT import alignment
pass#IMPORTIMPORTIMPORT import applications
pass#IMPORTIMPORTIMPORT import copy
pass#IMPORTIMPORTIMPORT import datetime
pass#IMPORTIMPORTIMPORT import development
pass#IMPORTIMPORTIMPORT import filter
pass#IMPORTIMPORTIMPORT import fundamentals
pass#IMPORTIMPORTIMPORT import global_def
pass#IMPORTIMPORTIMPORT import isac
pass#IMPORTIMPORTIMPORT import logger
pass#IMPORTIMPORTIMPORT import math
pass#IMPORTIMPORTIMPORT import morphology
pass#IMPORTIMPORTIMPORT import mpi
pass#IMPORTIMPORTIMPORT import multi_shc
pass#IMPORTIMPORTIMPORT import numpy
pass#IMPORTIMPORTIMPORT import numpy.random
pass#IMPORTIMPORTIMPORT import os
pass#IMPORTIMPORTIMPORT import pickle
pass#IMPORTIMPORTIMPORT import pixel_error
pass#IMPORTIMPORTIMPORT import projection
pass#IMPORTIMPORTIMPORT import random
pass#IMPORTIMPORTIMPORT import reconstruction
pass#IMPORTIMPORTIMPORT import scipy.optimize
pass#IMPORTIMPORTIMPORT import statistics
pass#IMPORTIMPORTIMPORT import string
pass#IMPORTIMPORTIMPORT import string as sting
pass#IMPORTIMPORTIMPORT import sys
pass#IMPORTIMPORTIMPORT import time
pass#IMPORTIMPORTIMPORT import types
pass#IMPORTIMPORTIMPORT import user_functions
pass#IMPORTIMPORTIMPORT import utilities
from builtins import range
from builtins import object
pass#IMPORTIMPORTIMPORT from global_def import *
pass#IMPORTIMPORTIMPORT import numpy.random

def ali2d_MPI(stack, outdir, maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", \
			nomirror = False, dst=0.0, center=-1, maxit=0, CTF=False, snr=1.0, \
			Fourvar=False, Ng=-1, user_func_name="ref_ali2d", CUDA=False, GPUID="", random_method = ""):

	pass#IMPORTIMPORTIMPORT from utilities    import model_circle, model_blank, drop_image, get_image, get_input_from_string
	pass#IMPORTIMPORTIMPORT from utilities    import reduce_EMData_to_root, bcast_EMData_to_all, send_attr_dict, file_type
	pass#IMPORTIMPORTIMPORT from utilities    import bcast_number_to_all, bcast_list_to_all
	pass#IMPORTIMPORTIMPORT from statistics   import fsc_mask, sum_oe, hist_list, varf2d_MPI
	pass#IMPORTIMPORTIMPORT from alignment    import Numrinit, ringwe, ali2d_single_iter
	pass#IMPORTIMPORTIMPORT from pixel_error  import pixel_error_2D
	pass#IMPORTIMPORTIMPORT from numpy        import reshape, shape
	pass#IMPORTIMPORTIMPORT from fundamentals import fshift, fft, rot_avg_table
	pass#IMPORTIMPORTIMPORT from utilities    import get_params2D, set_params2D
	pass#IMPORTIMPORTIMPORT from utilities    import print_msg, print_begin_msg, print_end_msg
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT import sys
	pass#IMPORTIMPORTIMPORT from mpi 	  import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from mpi 	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_gatherv
	pass#IMPORTIMPORTIMPORT from mpi 	  import MPI_SUM, MPI_FLOAT, MPI_INT

	number_of_proc = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
	myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
	main_node = 0
	
	ftp = sparx_utilities.file_type(stack)
	
	if outdir:
		if os.path.exists(outdir):  sparx_global_def.ERROR('Output directory exists, please change the name and restart the program', "ali2d_MPI", 1, myid)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

	if myid == main_node:
		if outdir:
			os.mkdir(outdir)
		pass#IMPORTIMPORTIMPORT import global_def
		sparx_global_def.LOGFILE =  os.path.join(outdir, sparx_global_def.LOGFILE)
		sparx_utilities.print_begin_msg("ali2d_MPI")

	xrng        = sparx_utilities.get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = sparx_utilities.get_input_from_string(yr)
	step        = sparx_utilities.get_input_from_string(ts)
	
	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit);
	
	if max_iter == 0:
		max_iter = 10
		auto_stop = True
	else:
		auto_stop = False

	if myid == main_node:
		if ftp == "bdb":
			pass#IMPORTIMPORTIMPORT from EMAN2db import db_open_dict
			dummy = EMAN2db.db_open_dict(stack, True)
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = EMUtil.get_all_attributes(stack, 'active')
		# list_of_particles = []
		# for im in xrange(len(active)):
		# 	if active[im]:  list_of_particles.append(im)
		# del active
		# nima = len(list_of_particles)
		
		nima = EMAN2_cppwrap.EMUtil.get_image_count(stack)
		list_of_particles = list(range(nima))
	
	else:
		nima = 0
	nima = sparx_utilities.bcast_number_to_all(nima, source_node = main_node)
	
	if myid != main_node:
		list_of_particles = [-1]*nima
	list_of_particles = sparx_utilities.bcast_list_to_all(list_of_particles, myid,  source_node = main_node)

	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
	list_of_particles = list_of_particles[image_start: image_end]
	
	if Ng == -1:
		Ng = nima
	elif Ng == -2:
		Ng = int(0.98*nima)	

	# read nx and ctf_app (if CTF) and broadcast to all nodes
	if myid == main_node:
		ima = EMAN2_cppwrap.EMData()
		ima.read_image(stack, list_of_particles[0], True)
		nx = ima.get_xsize()
		if CTF:	ctf_app = ima.get_attr_default('ctf_applied', 0)
		del ima
	else:
		nx = 0
		if CTF:	ctf_app = 0
	nx = sparx_utilities.bcast_number_to_all(nx, source_node = main_node)
	if CTF:
		ctf_app = sparx_utilities.bcast_number_to_all(ctf_app, source_node = main_node)
		if ctf_app > 0:	sparx_global_def.ERROR("data cannot be ctf-applied", "ali2d_MPI", 1, myid)
		phase_flip = True
		pass#IMPORTIMPORTIMPORT from filter import filt_ctf
	else:
		phase_flip = False
	CTF = False

	# default value for the last ring
	if last_ring == -1: last_ring = nx/2-2

	if last_ring + max([max(xrng), max(yrng)]) > (nx-1) // 2:
		sparx_global_def.ERROR('Shift or radius is too large - particle crosses image boundary', "ali2d_MPI", 1)

	if CUDA:
		GPUID = sparx_utilities.get_input_from_string(GPUID)
		GPUID = list(map(int, GPUID))
	
	if myid == main_node:
		sparx_utilities.print_msg("Input stack                 : %s\n"%(stack))
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# print_msg("Number of active images     : %d\n"%(nima))
		sparx_utilities.print_msg("Number of images            : %d\n"%(nima))
		sparx_utilities.print_msg("Output directory            : %s\n"%(outdir))
		sparx_utilities.print_msg("Inner radius                : %i\n"%(first_ring))
		sparx_utilities.print_msg("Outer radius                : %i\n"%(last_ring))
		sparx_utilities.print_msg("Ring step                   : %i\n"%(rstep))
		sparx_utilities.print_msg("X search range              : %s\n"%(xrng))
		sparx_utilities.print_msg("Y search range              : %s\n"%(yrng))
		sparx_utilities.print_msg("Translational step          : %s\n"%(step))
		sparx_utilities.print_msg("Disable checking mirror     : %s\n"%(nomirror))
		sparx_utilities.print_msg("Discrete angle used         : %d\n"%(dst))
		sparx_utilities.print_msg("Center type                 : %i\n"%(center))
		sparx_utilities.print_msg("Maximum iteration           : %i\n"%(max_iter))
		sparx_utilities.print_msg("Use Fourier variance        : %s\n"%(Fourvar))
		#print_msg("Number of groups            : %d\n"%(Ng))
		sparx_utilities.print_msg("CTF correction              : %s\n"%(CTF))
		sparx_utilities.print_msg("Phase flip                  : %s\n"%(phase_flip))
		sparx_utilities.print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		if auto_stop:
			sparx_utilities.print_msg("Stop iteration with         : criterion\n")
		else:
			sparx_utilities.print_msg("Stop iteration with         : maxit\n")

		pass#IMPORTIMPORTIMPORT import user_functions
		user_func = sparx_user_functions.factory[user_func_name]

		sparx_utilities.print_msg("User function               : %s\n"%(user_func_name))
		sparx_utilities.print_msg("Number of processors used   : %d\n"%(number_of_proc))
		sparx_utilities.print_msg("Using CUDA                  : %s\n"%(CUDA))
		if CUDA:
			sparx_utilities.print_msg("GPU IDs                     : %s\n"%(GPUID))

	if maskfile:
		pass#IMPORTIMPORTIMPORT import  types
		if type(maskfile) is bytes:  
			if myid == main_node:		sparx_utilities.print_msg("Maskfile                    : %s\n\n"%(maskfile))
			mask = sparx_utilities.get_image(maskfile)
		else:
			if myid == main_node: 		sparx_utilities.print_msg("Maskfile                    : user provided in-core mask\n\n")
			mask = maskfile
	else:
		if myid == main_node: 	sparx_utilities.print_msg("Maskfile                    : default, a circle with radius %i\n\n"%(last_ring))
		mask = sparx_utilities.model_circle(last_ring, nx, nx)

	cnx  = nx/2+1
	cny  = cnx
	if  random_method == "SCF":		mode = "H"
	else: 							mode = "F"
	data = []
	if CTF:
		pass#IMPORTIMPORTIMPORT from filter import filt_ctf
		pass#IMPORTIMPORTIMPORT from morphology   import ctf_img
		ctf_abs_sum = EMAN2_cppwrap.EMData(nx, nx, 1, False)
		ctf_2_sum = EMAN2_cppwrap.EMData(nx, nx, 1, False)
	else:
		ctf_2_sum = None

	pass#IMPORTIMPORTIMPORT from global_def import CACHE_DISABLE
	if sparx_global_def.CACHE_DISABLE:
		data = EMAN2_cppwrap.EMData.read_images(stack, list_of_particles)
	else:
		for i in range(number_of_proc):
			if myid == i:
				data = EMAN2_cppwrap.EMData.read_images(stack, list_of_particles)
			if ftp == "bdb": mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		
	if CUDA:
		nGPU = len(GPUID)
		GPUID = GPUID[myid%nGPU]
		R = CUDA_Aligner(GPUID)
		all_ali_params = []
		all_ctf_params = []

	for im in range(len(data)):
		data[im].set_attr('ID', list_of_particles[im])
		st = EMAN2_cppwrap.Util.infomask(data[im], mask, False)
		data[im] -= st[0]
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			ctfimg = sparx_morphology.ctf_img(nx, ctf_params)
			if CUDA:
				all_ctf_params.extend([ctf_params.defocus, ctf_params.cs, ctf_params.voltage, ctf_params.apix, ctf_params.bfactor, ctf_params.ampcont])
			EMAN2_cppwrap.Util.add_img2(ctf_2_sum, ctfimg)
			EMAN2_cppwrap.Util.add_img_abs(ctf_abs_sum, ctfimg)
		if CUDA:
			alpha, sx, sy, mirror, scale = sparx_utilities.get_params2D(data[im])
			all_ali_params.extend([alpha, sx, sy, mirror])
		if( random_method == "SHC" ):  data[im].set_attr('previousmax',1.0e-23)
		if phase_flip:  data[im] = sparx_filter.filt_ctf(data[im], data[im].get_attr("ctf"), binary = True)

	if CTF:
		sparx_utilities.reduce_EMData_to_root(ctf_2_sum, myid, main_node)
		sparx_utilities.reduce_EMData_to_root(ctf_abs_sum, myid, main_node)
		if myid == main_node:
			adw_img = EMAN2_cppwrap.Util.mult_scalar(ctf_2_sum, snr)
			EMAN2_cppwrap.Util.div_filter(adw_img, ctf_abs_sum)
			EMAN2_cppwrap.Util.mul_scalar(adw_img, float(Ng-1)/(nima-1))
			adw_img += float(nima-Ng)/(nima-1)
	else:  ctf_2_sum = None
	# startup
	numr = sparx_alignment.Numrinit(first_ring, last_ring, rstep, mode)  #precalculate rings
	wr = sparx_alignment.ringwe(numr, mode)
	
	if myid == main_node:
		# initialize data for the reference preparation function
		ref_data = [mask, center, None, None]
		sx_sum = 0.0
		sy_sum = 0.0
		a0 = -1.0e22
		
	recvcount = []
	disp = []
	for i in range(number_of_proc):
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
		pass#IMPORTIMPORTIMPORT from math import log, pi
		RING_LENGTH = 2**(int(numpy.log(2*numpy.pi*last_ring)/numpy.log(2))+1)
		NRING       = 2**(int(numpy.log(last_ring)/numpy.log(2))+1)

	for N_step in range(len(xrng)):

		if CUDA:
			R.setup(len(data), nx, nx, RING_LENGTH, NRING, last_ring, step[N_step], int(xrng[N_step]/step[N_step]+0.5), int(yrng[N_step]/step[N_step]+0.5), CTF)
			for im in range(len(data)):	R.insert_image(data[im], im)
			if CTF:  R.filter_stack(all_ctf_params)

		msg = "\nX range = %5.2f   Y range = %5.2f   Step = %5.2f\n"%(xrng[N_step], yrng[N_step], step[N_step])
		if myid == main_node: sparx_utilities.print_msg(msg)
		for Iter in range(max_iter):
			total_iter += 1
			if CUDA:
				ave1 = sparx_utilities.model_blank(nx, nx)
				ave2 = sparx_utilities.model_blank(nx, nx)
				R.sum_oe(all_ctf_params, all_ali_params, ave1, ave2)
				# Comment by Zhengfan Yang on 02/01/10
				# The reason for this step is that in CUDA 2-D FFT, the image is multipled by NX*NY times after
				# FFT and IFFT, so we want to decrease it such that the criterion is in line with non-CUDA version
				# However, this step is not mandatory.
				if CTF:
					ave1 /= (nx*2)**2
					ave2 /= (nx*2)**2
			else:
				ave1, ave2 = sparx_statistics.sum_oe(data, "a", CTF, EMAN2_cppwrap.EMData())  # pass empty object to prevent calculation of ctf^2
			sparx_utilities.reduce_EMData_to_root(ave1, myid, main_node)
			sparx_utilities.reduce_EMData_to_root(ave2, myid, main_node)
			if myid == main_node:
				sparx_utilities.print_msg("Iteration #%4d\n"%(total_iter))
				if CTF: 
					tavg_Ng = sparx_fundamentals.fft(EMAN2_cppwrap.Util.divn_filter(EMAN2_cppwrap.Util.muln_img(sparx_fundamentals.fft(EMAN2_cppwrap.Util.addn_img(ave1, ave2)), adw_img), ctf_2_sum))
					tavg    = sparx_fundamentals.fft(EMAN2_cppwrap.Util.divn_filter(sparx_fundamentals.fft(EMAN2_cppwrap.Util.addn_img(ave1, ave2)), ctf_2_sum))
				else:	 tavg = (ave1+ave2)/nima
				if outdir:
					tavg.write_image(os.path.join(outdir, "aqc.hdf"), total_iter-1)
					if CTF:
						tavg_Ng.write_image(os.path.join(outdir, "aqc_view.hdf"), total_iter-1)
					frsc = sparx_statistics.fsc_mask(ave1, ave2, mask, 1.0, os.path.join(outdir, "resolution%03d"%(total_iter)))
				else:
					frsc = sparx_statistics.fsc_mask(ave1, ave2, mask, 1.0)
			else:
				tavg =  sparx_utilities.model_blank(nx, nx)
			del ave1, ave2
				
			if Fourvar:  
				sparx_utilities.bcast_EMData_to_all(tavg, myid, main_node)
				vav, rvar = sparx_statistics.varf2d_MPI(myid, data, tavg, mask, "a", CTF)

			if myid == main_node:
				if Fourvar:
					tavg    = sparx_fundamentals.fft(EMAN2_cppwrap.Util.divn_img(sparx_fundamentals.fft(tavg), vav))
					vav_r	= EMAN2_cppwrap.Util.pack_complex_to_real(vav)
					if outdir:
						vav_r.write_image(os.path.join(outdir, "varf.hdf"), total_iter-1)


				# a0 should increase; stop algorithm when it decreases.    
				#     However, the result will depend on filtration, so it is not quite right.
				#  moved it here, so it is for unfiltered average and thus hopefully makes more sense
				a1 = tavg.cmp("dot", tavg, dict(negative = 0, mask = ref_data[0]))
				msg = "Criterion %d = %15.8e"%(total_iter, a1)
				# numpy.log.add(msg)


				ref_data[2] = tavg
				ref_data[3] = frsc

				#  call user-supplied function to prepare reference image, i.e., center and filter it
				if center == -1:
					# When center = -1, which is by default, we use the average center method
					ref_data[1] = 0
					tavg, cs = user_func(ref_data)
					cs[0] = float(sx_sum)/nima
					cs[1] = float(sy_sum)/nima
					tavg = sparx_fundamentals.fshift(tavg, -cs[0], -cs[1])
					msg = "Average center x =      %10.3f        Center y       = %10.3f\n"%(cs[0], cs[1])
					sparx_utilities.print_msg(msg)
				else:
					tavg, cs = user_func(ref_data)

				# write the current filtered average
				if outdir:
					tavg.write_image(os.path.join(outdir, "aqf.hdf"), total_iter-1)

				if a1 < a0:
					if auto_stop: 	again = 0
				else:	a0 = a1
			else:
				tavg = sparx_utilities.model_blank(nx, nx)
				cs = [0.0]*2

			if auto_stop:
				again = mpi.mpi_bcast(again, 1, mpi.MPI_INT, main_node, mpi.MPI_COMM_WORLD)
				if int(again[0]) == 0: break

			if Fourvar:  del vav
			sparx_utilities.bcast_EMData_to_all(tavg, myid, main_node)
			cs = mpi.mpi_bcast(cs, 2, mpi.MPI_FLOAT, main_node, mpi.MPI_COMM_WORLD)
			cs = list(map(float, cs))
			if total_iter != max_iter*len(xrng):
				if CUDA:
					old_ali_params = all_ali_params[:]
				else:
					old_ali_params = []
					for im in range(len(data)):  
						alpha, sx, sy, mirror, scale = sparx_utilities.get_params2D(data[im])
						old_ali_params.extend([alpha, sx, sy, mirror])

				if Iter%4 != 0 or total_iter > max_iter*len(xrng)-10: delta = 0.0 
				else: delta = dst
				if CUDA:
					all_ali_params = R.ali2d_single_iter(tavg, all_ali_params, cs[0], cs[1], 1, delta)
					sx_sum = all_ali_params[-2]
					sy_sum = all_ali_params[-1]
					for im in range(len(data)):  all_ali_params[im*4+3] = int(all_ali_params[im*4+3])
				else:
					sx_sum, sy_sum, nope = sparx_alignment.ali2d_single_iter(data, numr, wr, cs, tavg, cnx, cny, \
												xrng[N_step], yrng[N_step], step[N_step], \
												nomirror=nomirror, mode=mode, CTF=CTF, delta=delta, \
												random_method = random_method)

				sx_sum = mpi.mpi_reduce(sx_sum, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, main_node, mpi.MPI_COMM_WORLD)
				sy_sum = mpi.mpi_reduce(sy_sum, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, main_node, mpi.MPI_COMM_WORLD)
				#  for SHC
				if  random_method == "SHC":
					nope   = mpi.mpi_reduce(nope, 1, mpi.MPI_INT, mpi.MPI_SUM, main_node, mpi.MPI_COMM_WORLD)
					nope   = mpi.mpi_bcast(nope, 1, mpi.MPI_INT, main_node, mpi.MPI_COMM_WORLD)
					if int(nope[0]) == nima: break

				pixel_error       = 0.0
				mirror_consistent = 0
				pixel_error_list  = []
				for im in range(len(data)):
					if CUDA:
						alpha = all_ali_params[im*4]
						sx = all_ali_params[im*4+1]
						sy = all_ali_params[im*4+2]
						mirror = all_ali_params[im*4+3]
					else:
						alpha, sx, sy, mirror, scale = sparx_utilities.get_params2D(data[im])
						if old_ali_params[im*4+3] == mirror:
							this_error = sparx_pixel_error.pixel_error_2D(old_ali_params[im*4:im*4+3], [alpha, sx, sy], last_ring)
							pixel_error += this_error
							pixel_error_list.append(this_error)
							mirror_consistent += 1
						else:
							pixel_error_list.append(-1)
				mirror_consistent = mpi.mpi_reduce(mirror_consistent, 1, mpi.MPI_INT, mpi.MPI_SUM, main_node, mpi.MPI_COMM_WORLD)
				pixel_error       = mpi.mpi_reduce(pixel_error, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, main_node, mpi.MPI_COMM_WORLD)
				pixel_error_list  = mpi.mpi_gatherv(pixel_error_list, len(data), mpi.MPI_FLOAT, recvcount, disp, mpi.MPI_FLOAT, main_node, mpi.MPI_COMM_WORLD)
				if myid == main_node:
					sparx_utilities.print_msg("Mirror consistency rate = %8.4f%%\n"%(float(mirror_consistent)/nima*100))
					if mirror_consistent!=0:
						sparx_utilities.print_msg("Among the mirror-consistent images, average of pixel errors is %0.4f, and their distribution is:\n"%(float(pixel_error)/float(mirror_consistent)))
						pixel_error_list = list(map(float, pixel_error_list))
						for i in range(nima-1, -1, -1):
							if pixel_error_list[i] < 0:  del pixel_error_list[i]
						region, hist = sparx_statistics.hist_list(pixel_error_list, 20)	
						for p in range(20):
							sparx_utilities.print_msg("      %10.6f: %5d\n"%(region[p], hist[p]))
					sparx_utilities.print_msg("\n\n\n")
		if CUDA: R.finish()

	if CUDA:
		for im in range(len(data)):
			sparx_utilities.set_params2D(data[im], [all_ali_params[im*4], all_ali_params[im*4+1], all_ali_params[im*4+2], all_ali_params[im*4+3], 1.0])

	if myid == main_node and outdir:  sparx_utilities.drop_image(tavg, os.path.join(outdir, "aqfinal.hdf"))
	# write out headers and STOP, under MPI writing has to be done sequentially
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	par_str = ["xform.align2d", "ID"]
	if myid == main_node:
		pass#IMPORTIMPORTIMPORT from utilities import file_type
		if(sparx_utilities.file_type(stack) == "bdb"):
			pass#IMPORTIMPORTIMPORT from utilities import recv_attr_dict_bdb
			sparx_utilities.recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
		else:
			pass#IMPORTIMPORTIMPORT from utilities import recv_attr_dict
			sparx_utilities.recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	else:           sparx_utilities.send_attr_dict(main_node, data, par_str, image_start, image_end)
	if myid == main_node: sparx_utilities.print_end_msg("ali2d_MPI")



def ali2d_base(stack, outdir, maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", \
			nomirror = False, dst=0.0, center=-1, maxit=0, CTF=False, snr=1.0, \
			Fourvar=False, user_func_name="ref_ali2d", random_method = "", log = None, \
			number_of_proc = 1, myid = 0, main_node = 0, mpi_comm = None, write_headers = False):

	pass#IMPORTIMPORTIMPORT from utilities    import model_circle, model_blank, drop_image, get_image, get_input_from_string
	pass#IMPORTIMPORTIMPORT from utilities    import reduce_EMData_to_root, bcast_EMData_to_all, send_attr_dict, file_type
	pass#IMPORTIMPORTIMPORT from utilities    import bcast_number_to_all, bcast_list_to_all
	pass#IMPORTIMPORTIMPORT from statistics   import fsc_mask, sum_oe, hist_list, varf2d_MPI
	pass#IMPORTIMPORTIMPORT from alignment    import Numrinit, ringwe, ali2d_single_iter
	pass#IMPORTIMPORTIMPORT from pixel_error  import pixel_error_2D
	pass#IMPORTIMPORTIMPORT from numpy        import reshape, shape
	pass#IMPORTIMPORTIMPORT from fundamentals import fshift, fft, rot_avg_table
	pass#IMPORTIMPORTIMPORT from utilities    import get_params2D, set_params2D
	pass#IMPORTIMPORTIMPORT from utilities    import wrap_mpi_gatherv
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT import sys
	pass#IMPORTIMPORTIMPORT from mpi 	  import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	pass#IMPORTIMPORTIMPORT from mpi 	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_gatherv
	pass#IMPORTIMPORTIMPORT from mpi 	  import MPI_SUM, MPI_FLOAT, MPI_INT

	if log == None:
		pass#IMPORTIMPORTIMPORT from logger import Logger
		log = sparx_logger.Logger()

	if mpi_comm == None:
		mpi_comm = mpi.MPI_COMM_WORLD

	# ftp = file_type(stack)

	if myid == main_node:
		pass#IMPORTIMPORTIMPORT import global_def
		sparx_global_def.LOGFILE =  os.path.join(outdir, sparx_global_def.LOGFILE)
		log.add("Start  ali2d_MPI")

	xrng        = sparx_utilities.get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = sparx_utilities.get_input_from_string(yr)
	step        = sparx_utilities.get_input_from_string(ts)
	
	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit);
	
	if max_iter == 0:
		max_iter = 10
		auto_stop = True
	else:
		auto_stop = False

	pass#IMPORTIMPORTIMPORT import types
	if( type(stack) is bytes ):
		if myid == main_node:
			total_nima = EMAN2_cppwrap.EMUtil.get_image_count(stack)
		else:
			total_nima = 0
		total_nima = sparx_utilities.bcast_number_to_all(total_nima)
		list_of_particles = list(range(total_nima))

		image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
		list_of_particles = list_of_particles[image_start:image_end]
		nima = len(list_of_particles)
		data = EMAN2_cppwrap.EMData.read_images(stack, list_of_particles)

	else:
		data = stack
		total_nima = len(data)
		total_nima = mpi.mpi_reduce(total_nima, 1, mpi.MPI_INT, mpi.MPI_SUM, 0, mpi.MPI_COMM_WORLD)
		total_nima = mpi.mpi_bcast(total_nima, 1, mpi.MPI_INT, main_node, mpi.MPI_COMM_WORLD)[0]
		list_of_particles = list(range(total_nima))
		image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
		list_of_particles = list_of_particles[image_start:image_end]
		nima = len(list_of_particles)

		print('nima is' , nima)
		print('data length is ', len(data))
		assert( nima == len(data))

	# read nx and ctf_app (if CTF) and broadcast to all nodes
	if myid == main_node:
		nx = data[0].get_xsize()
		if CTF:	ctf_app = data[0].get_attr_default('ctf_applied', 0)
		# del ima
	else:
		nx = 0
		if CTF:	ctf_app = 0
	nx = sparx_utilities.bcast_number_to_all(nx, source_node = main_node)
	if CTF:
		ctf_app = sparx_utilities.bcast_number_to_all(ctf_app, source_node = main_node)
		if ctf_app > 0:	sparx_global_def.ERROR("data cannot be ctf-applied", "ali2d_MPI", 1, myid)
		phase_flip = True
		pass#IMPORTIMPORTIMPORT from filter import filt_ctf
	else:
		phase_flip = False
	CTF = False

	# default value for the last ring
	if last_ring == -1: last_ring = nx/2-2

	if last_ring + max([max(xrng), max(yrng)]) > (nx-1) // 2:
		sparx_global_def.ERROR('Shift or radius is too large - particle crosses image boundary', "ali2d_MPI", 1)
	
	if myid == main_node:
		# log.add("Input stack                 : %s"%(stack))
		log.add("Number of images            : %d"%(total_nima))
		log.add("Output directory            : %s"%(outdir))
		log.add("Inner radius                : %i"%(first_ring))
		log.add("Outer radius                : %i"%(last_ring))
		log.add("Ring step                   : %i"%(rstep))
		log.add("X search range              : %s"%(xrng))
		log.add("Y search range              : %s"%(yrng))
		log.add("Translational step          : %s"%(step))
		log.add("Disable checking mirror     : %s"%(nomirror))
		log.add("Discrete angle used         : %d"%(dst))
		log.add("Center type                 : %i"%(center))
		log.add("Maximum iteration           : %i"%(max_iter))
		#log.add("Use Fourier variance        : %s\n"%(Fourvar))
		log.add("CTF correction              : %s"%(CTF))
		log.add("Phase flip                  : %s"%(phase_flip))
		#log.add("Signal-to-Noise Ratio       : %f\n"%(snr))
		if auto_stop:
			log.add("Stop iteration with         : criterion")
		else:
			log.add("Stop iteration with         : maxit")

		pass#IMPORTIMPORTIMPORT import user_functions
		user_func = sparx_user_functions.factory[user_func_name]

		log.add("User function               : %s"%(user_func_name))
		log.add("Number of processors used   : %d"%(number_of_proc))

	if maskfile:
		pass#IMPORTIMPORTIMPORT import  types
		if type(maskfile) is bytes:  
			if myid == main_node:		log.add("Maskfile                    : %s"%(maskfile))
			mask = sparx_utilities.get_image(maskfile)
		else:
			if myid == main_node: 		log.add("Maskfile                    : user provided in-core mask")
			mask = maskfile
	else:
		if myid == main_node: 	log.add("Maskfile                    : default, a circle with radius %i"%(last_ring))
		mask = sparx_utilities.model_circle(last_ring, nx, nx)

	cnx  = nx/2+1
	cny  = cnx
	if  random_method == "SCF":		mode = "H"
	else: 							mode = "F"

	if CTF:
		pass#IMPORTIMPORTIMPORT from filter import filt_ctf
		pass#IMPORTIMPORTIMPORT from morphology   import ctf_img
		ctf_abs_sum = EMAN2_cppwrap.EMData(nx, nx, 1, False)
		ctf_2_sum = EMAN2_cppwrap.EMData(nx, nx, 1, False)
	else:
		ctf_2_sum = None

	for im in range(nima):
		data[im].set_attr('ID', list_of_particles[im])
		sparx_utilities.set_params2D(data[im], [0.0, 0.0, 0.0, 0, 1.0], 'xform.align2d')
		st = EMAN2_cppwrap.Util.infomask(data[im], mask, False)
		data[im] -= st[0]
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			ctfimg = sparx_morphology.ctf_img(nx, ctf_params)
			EMAN2_cppwrap.Util.add_img2(ctf_2_sum, ctfimg)
			EMAN2_cppwrap.Util.add_img_abs(ctf_abs_sum, ctfimg)
		if( random_method == "SHC" ):  data[im].set_attr('previousmax',1.0e-23)
		if phase_flip:  data[im] = sparx_filter.filt_ctf(data[im], data[im].get_attr("ctf"), binary = True)

	if CTF:
		sparx_utilities.reduce_EMData_to_root(ctf_2_sum, myid, main_node)
		sparx_utilities.reduce_EMData_to_root(ctf_abs_sum, myid, main_node)
		if myid == main_node:
			adw_img = EMAN2_cppwrap.Util.mult_scalar(ctf_2_sum, snr)
			EMAN2_cppwrap.Util.div_filter(adw_img, ctf_abs_sum)
			EMAN2_cppwrap.Util.mul_scalar(adw_img, float(Ng-1)/(nima-1))
			adw_img += float(nima-Ng)/(nima-1)
	else:  ctf_2_sum = None

	# startup
	numr = sparx_alignment.Numrinit(first_ring, last_ring, rstep, mode) 	#precalculate rings
	wr = sparx_alignment.ringwe(numr, mode)

	if myid == main_node:
		# initialize data for the reference preparation function
		ref_data = [mask, center, None, None]
		sx_sum = 0.0
		sy_sum = 0.0
		a0 = -1.0e22
		
	recvcount = []
	disp = []
	for i in range(number_of_proc):
		ib, ie = MPI_start_end(total_nima, number_of_proc, i)
		recvcount.append(ie-ib)
		if i == 0:
			disp.append(0)
		else:
			disp.append(disp[i-1]+recvcount[i-1])

	again = 1
	total_iter = 0
	cs = [0.0]*2
	delta = 0.0
	for N_step in range(len(xrng)):

		for Iter in range(max_iter):
			total_iter += 1
			ave1, ave2 = sparx_statistics.sum_oe(data, "a", CTF, EMAN2_cppwrap.EMData())  # pass empty object to prevent calculation of ctf^2
			sparx_utilities.reduce_EMData_to_root(ave1, myid, main_node)
			sparx_utilities.reduce_EMData_to_root(ave2, myid, main_node)
			sys.stdout.flush()
			if myid == main_node:
				log.add("Iteration #%4d"%(total_iter))
				msg = "X range = %5.2f   Y range = %5.2f   Step = %5.2f"%(xrng[N_step], yrng[N_step], step[N_step])
				log.add(msg)
				if CTF: 
					tavg_Ng = sparx_fundamentals.fft(EMAN2_cppwrap.Util.divn_filter(EMAN2_cppwrap.Util.muln_img(sparx_fundamentals.fft(EMAN2_cppwrap.Util.addn_img(ave1, ave2)), adw_img), ctf_2_sum))
					tavg    = sparx_fundamentals.fft(EMAN2_cppwrap.Util.divn_filter(sparx_fundamentals.fft(EMAN2_cppwrap.Util.addn_img(ave1, ave2)), ctf_2_sum))
				else:
					tavg = (ave1+ave2)/total_nima
				if outdir:
					print('without CTF address', outdir)
					tavg.write_image(os.path.join(outdir, "aqc.hdf"), total_iter-1)
					
					if CTF:
						print('with CTF address', outdir)
						tavg_Ng.write_image(os.path.join(outdir, "aqc_view.hdf"), total_iter-1)
					frsc = sparx_statistics.fsc_mask(ave1, ave2, mask, 1.0, os.path.join(outdir, "resolution%03d"%(total_iter)))
				else:
					frsc = sparx_statistics.fsc_mask(ave1, ave2, mask, 1.0)
			else:
				tavg =  sparx_utilities.model_blank(nx, nx)
			del ave1, ave2
				
			if Fourvar:  
				sparx_utilities.bcast_EMData_to_all(tavg, myid, main_node)
				vav, rvar = sparx_statistics.varf2d_MPI(myid, data, tavg, mask, "a", CTF)

			if myid == main_node:
				if Fourvar:
					tavg    = sparx_fundamentals.fft(EMAN2_cppwrap.Util.divn_img(sparx_fundamentals.fft(tavg), vav))
					vav_r	= EMAN2_cppwrap.Util.pack_complex_to_real(vav)
					if outdir:
						vav_r.write_image(os.path.join(outdir, "varf.hdf"), total_iter-1)


				# a0 should increase; stop algorithm when it decreases.    
				#     However, the result will depend on filtration, so it is not quite right.
				#  moved it here, so it is for unfiltered average and thus hopefully makes more sense
				a1 = tavg.cmp("dot", tavg, dict(negative = 0, mask = ref_data[0]))
				msg = "Criterion %d = %15.8e"%(total_iter, a1)
				log.add(msg)


				ref_data[2] = tavg
				ref_data[3] = frsc

				#  call user-supplied function to prepare reference image, i.e., center and filter it
				if center == -1:
					# When center = -1, which is by default, we use the average center method
					ref_data[1] = 0
					tavg, cs = user_func(ref_data)
					cs[0] = float(sx_sum)/total_nima
					cs[1] = float(sy_sum)/total_nima
					tavg = sparx_fundamentals.fshift(tavg, -cs[0], -cs[1])
					msg = "Average center x =      %10.3f        Center y       = %10.3f"%(cs[0], cs[1])
					log.add(msg)
				else:
					if delta != 0.0:
						cnt = ref_data[1]
						ref_data[1] = 0
					tavg, cs = user_func(ref_data)
					if delta != 0.0:
						ref_data[1] = cnt
				# write the current filtered average
				if outdir:
					tavg.write_image(os.path.join(outdir, "aqf.hdf"), total_iter-1)

				if a1 < a0:
					if auto_stop: 	again = 0
				else:	a0 = a1
			else:
				tavg = sparx_utilities.model_blank(nx, nx)
				cs = [0.0]*2

			if auto_stop:
				again = mpi.mpi_bcast(again, 1, mpi.MPI_INT, main_node, mpi_comm)
				if int(again[0]) == 0: break

			if Fourvar:  del vav
			sparx_utilities.bcast_EMData_to_all(tavg, myid, main_node)
			cs = mpi.mpi_bcast(cs, 2, mpi.MPI_FLOAT, main_node, mpi_comm)
			cs = list(map(float, cs))
			if total_iter != max_iter*len(xrng):
				old_ali_params = []
				for im in range(nima):  
					alpha, sx, sy, mirror, scale = sparx_utilities.get_params2D(data[im])
					old_ali_params.extend([alpha, sx, sy, mirror])

				if Iter%4 != 0 or total_iter > max_iter*len(xrng)-10: delta = 0.0
				else: delta = dst
				sx_sum, sy_sum, nope = sparx_alignment.ali2d_single_iter(data, numr, wr, cs, tavg, cnx, cny, \
												xrng[N_step], yrng[N_step], step[N_step], \
												nomirror=nomirror, mode=mode, CTF=CTF, delta=delta, \
												random_method = random_method)

				sx_sum = mpi.mpi_reduce(sx_sum, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, main_node, mpi_comm)
				sy_sum = mpi.mpi_reduce(sy_sum, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, main_node, mpi_comm)
				#  for SHC
				if  random_method == "SHC":
					nope   = mpi.mpi_reduce(nope, 1, mpi.MPI_INT, mpi.MPI_SUM, main_node, mpi_comm)
					nope   = mpi.mpi_bcast(nope, 1, mpi.MPI_INT, main_node, mpi_comm)
					if int(nope[0]) == total_nima: break

				pixel_error       = 0.0
				mirror_consistent = 0
				pixel_error_list  = [-1.0]*nima
				for im in range(nima):
					alpha, sx, sy, mirror, scale = sparx_utilities.get_params2D(data[im])
					if old_ali_params[im*4+3] == mirror:
						this_error = sparx_pixel_error.pixel_error_2D(old_ali_params[im*4:im*4+3], [alpha, sx, sy], last_ring)
						pixel_error += this_error
						pixel_error_list[im] = this_error
						mirror_consistent += 1
				del old_ali_params
				mirror_consistent = mpi.mpi_reduce(mirror_consistent, 1, mpi.MPI_INT, mpi.MPI_SUM, main_node, mpi_comm)
				pixel_error       = mpi.mpi_reduce(pixel_error, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, main_node, mpi_comm)
				pixel_error_list  = mpi.mpi_gatherv(pixel_error_list, nima, mpi.MPI_FLOAT, recvcount, disp, mpi.MPI_FLOAT, main_node, mpi_comm)
				if myid == main_node:
					log.add("Mirror consistency rate = %8.4f%%"%(float(mirror_consistent)/total_nima*100))
					if mirror_consistent!=0:
						log.add("Among the mirror-consistent images, average of pixel errors is %0.4f, and their distribution is:"%(float(pixel_error)/float(mirror_consistent)))
						pixel_error_list = list(map(float, pixel_error_list))
						for i in range(total_nima-1, -1, -1):
							if pixel_error_list[i] < 0:  del pixel_error_list[i]
						region, hist = sparx_statistics.hist_list(pixel_error_list, 20)
						for p in range(20):
							log.add("      %14.2f: %6d"%(region[p], hist[p]))
					log.add("\n\n")

	if myid == main_node and outdir:  sparx_utilities.drop_image(tavg, os.path.join(outdir, "aqfinal.hdf"))
	# write out headers and STOP, under MPI writing has to be done sequentially
	mpi.mpi_barrier(mpi_comm)
	if write_headers:
		par_str = ["xform.align2d", "ID"]
		if myid == main_node:
			pass#IMPORTIMPORTIMPORT from utilities import file_type
			if(sparx_utilities.file_type(stack) == "bdb"):
				pass#IMPORTIMPORTIMPORT from utilities import recv_attr_dict_bdb
				sparx_utilities.recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
			else:
				pass#IMPORTIMPORTIMPORT from utilities import recv_attr_dict
				sparx_utilities.recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
		else:           sparx_utilities.send_attr_dict(main_node, data, par_str, image_start, image_end)
	params = []
	for im in range(nima):  
		alpha, sx, sy, mirror, scale = sparx_utilities.get_params2D(data[im])
		params.append([alpha, sx, sy, mirror])
	#params = wrap_mpi_gatherv(params, main_node, mpi_comm)

	if myid == main_node: log.add("Finished ali2d_base")

	return params #, data

"""Multiline Comment1"""

def mref_ali3d_MPI(stack, ref_vol, outdir, maskfile=None, focus = None, maxit=1, ir=1, ou=-1, rs=1, \
            xr ="4 2  2  1", yr="-1", ts="1 1 0.5 0.25",   delta="10  6  4  4", an="-1", center = -1, \
            nassign = 3, nrefine= 1, CTF = False, snr = 1.0,  ref_a="S", sym="c1",
			user_func_name="ref_ali3d", npad = 2, debug = False, fourvar=False, termprec = 0.0,\
			mpi_comm = None, log = None):
	pass#IMPORTIMPORTIMPORT from utilities      import model_circle, reduce_EMData_to_root, bcast_EMData_to_all, bcast_number_to_all, drop_image
	pass#IMPORTIMPORTIMPORT from utilities      import bcast_list_to_all, get_image, get_input_from_string, get_im
	pass#IMPORTIMPORTIMPORT from utilities      import get_arb_params, set_arb_params, drop_spider_doc, send_attr_dict
	pass#IMPORTIMPORTIMPORT from utilities      import get_params_proj, set_params_proj, model_blank, wrap_mpi_bcast
	pass#IMPORTIMPORTIMPORT from filter         import filt_params, filt_btwl, filt_ctf, filt_table, fit_tanh, filt_tanl
	pass#IMPORTIMPORTIMPORT from utilities      import rotate_3D_shift,estimate_3D_center_MPI
	pass#IMPORTIMPORTIMPORT from alignment      import Numrinit, prepare_refrings, proj_ali_incore
	pass#IMPORTIMPORTIMPORT from random         import randint, random
	pass#IMPORTIMPORTIMPORT from filter         import filt_ctf
	pass#IMPORTIMPORTIMPORT from utilities      import print_begin_msg, print_end_msg, print_msg
	pass#IMPORTIMPORTIMPORT from projection     import prep_vol, prgs, project, prgq, gen_rings_ctf
	pass#IMPORTIMPORTIMPORT from morphology     import binarize

	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT from mpi            import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	pass#IMPORTIMPORTIMPORT from mpi            import mpi_reduce, mpi_gatherv, mpi_scatterv, MPI_INT, MPI_SUM


	if mpi_comm == None: mpi_comm = mpi.MPI_COMM_WORLD

	if log == None:
		pass#IMPORTIMPORTIMPORT from logger import Logger
		log = sparx_logger.Logger()

	number_of_proc = mpi.mpi_comm_size(mpi_comm)
	myid           = mpi.mpi_comm_rank(mpi_comm)
	main_node = 0

	if myid == 0:
		if os.path.exists(outdir):  nx = 1
		else:  nx = 0
	else:  nx = 0
	ny = sparx_utilities.bcast_number_to_all(nx, source_node = main_node)
	
	if ny == 1:  sparx_global_def.ERROR('Output directory exists, please change the name and restart the program', "mref_ali3d_MPI", 1,myid)
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

	if myid == main_node:	
		os.mkdir(outdir)
		pass#IMPORTIMPORTIMPORT import global_def
		sparx_global_def.LOGFILE =  os.path.join(outdir, sparx_global_def.LOGFILE)
		log.add("Equal Kmeans-modified K-means  ")
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

	pass#IMPORTIMPORTIMPORT from time import time

	if debug:
		pass#IMPORTIMPORTIMPORT from time import sleep
		while not os.path.exists(outdir):
			print("Node ",myid,"  waiting...")
			time.sleep(5)

		finfo = open(os.path.join(outdir, "progress%04d"%myid), 'w')
		frec  = open( os.path.join(outdir, "recons%04d"%myid), "w" )
	else:
		finfo = None
		frec  = None

	xrng        = sparx_utilities.get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = sparx_utilities.get_input_from_string(yr)
	step        = sparx_utilities.get_input_from_string(ts)
	delta       = sparx_utilities.get_input_from_string(delta)
	lstp = min( len(xrng), len(yrng), len(step), len(delta) )
	if (an == "-1"):
		an = []
		for i in range(len(xrng)):   an.append(-1)
	else:
		pass#IMPORTIMPORTIMPORT from  alignment	    import proj_ali_incore_local
		an      = sparx_utilities.get_input_from_string(an)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	center      = int(center)

	numref = EMAN2_cppwrap.EMUtil.get_image_count(ref_vol)
	volref     = EMAN2_cppwrap.EMData()
	volref.read_image(stack, 0)
	nx      = volref.get_xsize()
	if last_ring < 0:	last_ring = nx//2 - 2

	if (myid == main_node):
		pass#IMPORTIMPORTIMPORT import user_functions
		user_func = sparx_user_functions.factory[user_func_name]
		log.add("mref_ali3d_MPI")
		log.add("Input stack                               : %s"%(stack))
		log.add("Reference volumes                         : %s"%(ref_vol))	
		log.add("Number of reference volumes               : %i"%(numref))
		log.add("Output directory                          : %s"%(outdir))
		log.add("User function                             : %s"%(user_func_name))
		if(focus != None):  \
		log.add("Maskfile 3D for focused clustering        : %s"%(focus))
		log.add("Overall 3D mask applied in user function  : %s"%(maskfile))
		log.add("Inner radius                              : %i"%(first_ring))
		log.add("Outer radius                              : %i"%(last_ring))
		log.add("Ring step                                 : %i"%(rstep))
		log.add("X search range                            : %s"%(xrng))
		log.add("Y search range                            : %s"%(yrng))
		log.add("Translational step                        : %s"%(step))
		log.add("Angular step                              : %s"%(delta))
		log.add("Angular search range                      : %s"%(an))
		log.add("Number of assignments in each iteration   : %i"%(nassign))
		log.add("Number of alignments in each iteration    : %i"%(nrefine))
		log.add("Number of iterations                      : %i"%(lstp*maxit) )
		log.add("Center type                               : %i"%(center))
		log.add("CTF correction                            : %s"%(CTF))
		log.add("Signal-to-Noise Ratio                     : %f"%(snr))
		log.add("Reference projection method               : %s"%(ref_a))
		log.add("Symmetry group                            : %s"%(sym))
		log.add("Percentage of change for termination      : %f"%(termprec))
		log.add("User function                             : %s"%(user_func_name))

	if(maskfile):
		if(type(maskfile) is bytes): mask3D = sparx_utilities.get_image(maskfile)
		else: 	                                mask3D = maskfile
	else        :  mask3D = sparx_utilities.model_circle(last_ring, nx, nx, nx)

	numr     = sparx_alignment.Numrinit(first_ring, last_ring, rstep, "F")
	mask2D   = sparx_utilities.model_circle(last_ring, nx, nx)
	if(first_ring > 1):  mask2D -= sparx_utilities.model_circle(first_ring, nx, nx)


	if( type(stack) is bytes ):
		if myid == main_node:
			total_nima = EMAN2_cppwrap.EMUtil.get_image_count( stack )
		else:
			total_nima = 0
		total_nima = sparx_utilities.wrap_mpi_bcast(total_nima, main_node, mpi_comm)
		list_of_particles = list(range(total_nima))
		image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
		# create a list of images for each node
		list_of_particles = list_of_particles[image_start: image_end]
		nima = len(list_of_particles)

	else:
		list_of_particles = list(range(len(stack)))
		nima = len(list_of_particles)
		total_nima = len(list_of_particles)
		total_nima = mpi.mpi_reduce(total_nima, 1, mpi.MPI_INT, mpi.MPI_SUM, 0, mpi.MPI_COMM_WORLD)
		total_nima = mpi.mpi_bcast(total_nima, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD)
		total_nima = int(total_nima[0])
	"""Multiline Comment18"""
	if debug:
		finfo.write( "Image_start, image_end: %d %d\n" %(image_start, image_end) )
		finfo.flush()

	start_time = time.time()
	data = [None]*nima
	#  Here the assumption is that input are always volumes.  It should be most likely be changed so optionally these are group assignments.
	#  Initialize Particle ID and set group number to non-existant -1
	for im in range(nima):
		if( type(stack) is bytes ):
			data[im] = sparx_utilities.get_im(stack, list_of_particles[im])
			data[im].set_attr_dict({'ID':list_of_particles[im], 'group':-1})
		else:
			data[im] = stack[list_of_particles[im]]
			#  NOTE: in case data comes in, it would have to have ID set as there is no way to tell here what was the original ordering.
			data[im].set_attr_dict({ 'group':-1})
	if(myid == 0):
		log.add( "Time to read data: %d" % (time.time()-start_time) );start_time = time.time()

	if fourvar:
		pass#IMPORTIMPORTIMPORT from reconstruction import rec3D_MPI
		pass#IMPORTIMPORTIMPORT from statistics     import varf3d_MPI
		#  Compute Fourier variance
		vol, fscc = sparx_reconstruction.rec3D_MPI(data, snr, sym, sparx_utilities.model_circle(last_ring, nx, nx, nx), os.path.join(outdir, "resolution0000"), myid, main_node, finfo=frec, npad=npad)
		varf = sparx_statistics.varf3d_MPI(data, os.path.join(outdir, "ssnr0000"), None, vol, last_ring, 1.0, 1, CTF, 1, sym, myid)
		if myid == main_node:   
			varf = 1.0/varf
			varf.write_image( os.path.join(outdir,"varf0000.hdf") )
	else:
		varf = None

	if myid == main_node:
		refdata = [None]*7
		for  iref in range(numref):
			vol = sparx_utilities.get_im(ref_vol, iref).write_image(os.path.join(outdir, "vol0000.hdf"), iref)
		refdata[0] = numref
		refdata[1] = outdir
		refdata[2] = None
		refdata[3] = 0
		#refdata[4] = varf
		refdata[5] = mask3D
		refdata[6] = False # whether to align on 50S, this only happens at refinement step
		user_func( refdata )
		#vol.write_image(os.path.join(outdir, "volf0000.hdf"), iref)
	mpi.mpi_barrier( mpi.MPI_COMM_WORLD )

	if CTF:
		if(data[0].get_attr_default("ctf_applied",0) > 0):  sparx_global_def.ERROR("mref_ali3d_MPI does not work for CTF-applied data", "mref_ali3d_MPI", 1, myid)
		pass#IMPORTIMPORTIMPORT from reconstruction import rec3D_MPI
	else:
		pass#IMPORTIMPORTIMPORT from reconstruction import rec3D_MPI_noCTF

	if debug:
		finfo.write( '%d loaded  \n' % len(data) )
		finfo.flush()

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in range(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                   disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )

	total_iter = 0
	tr_dummy = EMAN2_cppwrap.Transform({"type":"spider"})

	if(focus != None):
		if(myid == main_node):
			vol = sparx_utilities.get_im(focus)
		else:
			vol =  sparx_utilities.model_blank(nx, nx, nx)
		sparx_utilities.bcast_EMData_to_all(vol, myid, main_node)
		focus, kb = sparx_projection.prep_vol(vol)

	Niter = int(lstp*maxit*(nassign + nrefine) )
	for Iter in range(Niter):
		N_step = (Iter%(lstp*(nassign+nrefine)))/(nassign+nrefine)
		if Iter%(nassign+nrefine) < nassign:
			runtype = "ASSIGNMENT"
		else:
			runtype = "REFINEMENT"

		total_iter += 1
		if(myid == main_node):
			log.add("\n%s ITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f"%(runtype, total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step]))
			start_ime = time.time()
	
		peaks =  [ [ -1.0e23 for im in range(nima) ] for iref in range(numref) ]
		if runtype=="REFINEMENT":
			trans = [ [ tr_dummy for im in range(nima) ] for iref in range(numref) ]
			pixer = [ [  0.0     for im in range(nima) ] for iref in range(numref) ]
			if(an[N_step] > 0):
				pass#IMPORTIMPORTIMPORT from utilities    import even_angles
				ref_angles = sparx_utilities.even_angles(delta[N_step], symmetry=sym, method = ref_a, phiEqpsi = "Zero")
				# generate list of angles
				pass#IMPORTIMPORTIMPORT from alignment import generate_list_of_reference_angles_for_search
				list_of_reference_angles = \
				sparx_alignment.generate_list_of_reference_angles_for_search(ref_angles, sym=sym)
				del ref_angles
			else:  list_of_reference_angles = [[1.0,1.0]]

		cs = [0.0]*3
		for iref in range(numref):
			if(myid == main_node):
				volft = sparx_utilities.get_im(os.path.join(outdir, "volf%04d.hdf"%(total_iter-1)), iref)
			else:
				volft =  sparx_utilities.model_blank(nx, nx, nx)
			sparx_utilities.bcast_EMData_to_all(volft, myid, main_node)

			volft, kb = sparx_projection.prep_vol(volft)
			if CTF:
				previous_defocus = -1.0
				if runtype=="REFINEMENT":
					start_time = time.time()
					prjref = sparx_projection.prgq( volft, kb, nx, delta[N_step], ref_a, sym, MPI=True)
					if(myid == 0):
						log.add( "Calculation of projections: %d" % (time.time()-start_time) );start_time = time.time()
					del volft, kb

			else:
				if runtype=="REFINEMENT":
					start_time = time.time()
					refrings = sparx_alignment.prepare_refrings( volft, kb, nx, delta[N_step], ref_a, sym, numr)
					if(myid == 0):
						log.add( "Initial time to prepare rings: %d" % (time.time()-start_time) );start_time = time.time()
					del volft, kb


			start_time = time.time()
			for im in range(nima):
				if(CTF):
					ctf = data[im].get_attr( "ctf" )
					if runtype=="REFINEMENT":
						if(ctf.defocus != previous_defocus):
							previous_defocus = ctf.defocus
							rstart_time = time.time()
							refrings = sparx_projection.gen_rings_ctf( prjref, nx, ctf, numr)
							if(myid == 0):
								log.add( "Repeated time to prepare rings: %d" % (time.time()-rstart_time) );rstart_time = time.time()

				if runtype=="ASSIGNMENT":
					phi,tht,psi,s2x,s2y = sparx_utilities.get_params_proj(data[im])
					ref = sparx_projection.prgs( volft, kb, [phi,tht,psi,-s2x,-s2y])
					if CTF:  ref = sparx_filter.filt_ctf( ref, ctf )
					if(focus != None):  mask2D = sparx_morphology.binarize( sparx_projection.prgs( focus, kb, [phi,tht,psi,-s2x,-s2y]) )  #  Should be precalculated!!
					peak = ref.cmp("ccc",data[im],{"mask":mask2D, "negative":0})
					if not(finfo is None):
						finfo.write( "ID, iref, peak: %6d %d %8.5f\n" % (list_of_particles[im],iref,peak) )
				else:
					if(an[N_step] == -1):
						peak, pixel_error = sparx_alignment.proj_ali_incore(data[im], refrings, numr, xrng[N_step], yrng[N_step], step[N_step])
					else:
						peak, pixel_error = sparx_alignment.proj_ali_incore_local(data[im], refrings, list_of_reference_angles, numr, \
																	xrng[N_step], yrng[N_step], step[N_step], an[N_step],sym=sym)
					if not(finfo is None):
						phi,tht,psi,s2x,s2y = sparx_utilities.get_params_proj(data[im])
						finfo.write( "ID, iref, peak,t rans: %6d %d %f %f %f %f %f %f\n"%(list_of_particles[im],iref,peak,phi,tht,psi,s2x,s2y) )
						finfo.flush()

				peaks[iref][im] = peak
				if runtype=="REFINEMENT":
					pixer[iref][im] = pixel_error
					trans[iref][im] = data[im].get_attr( "xform.projection" )

			if(myid == 0):
				log.add( "Time to process particles for reference %3d: %d" % (iref, time.time()-start_time) );start_time = time.time()


		if runtype=="ASSIGNMENT":  del volft, kb, ref
		else:
			if CTF: del prjref
			del refrings
			if(an[N_step] > 0): del list_of_reference_angles


		#  send peak values to the main node, do the assignments, and bring them back
		pass#IMPORTIMPORTIMPORT from numpy import float32, empty, inner, abs
		if( myid == 0 ):
			dtot = numpy.empty( (numref, total_nima), dtype = numpy.float32)
		for  iref in range(numref):
			recvbuf = mpi.mpi_gatherv(peaks[iref], nima, mpi.MPI_FLOAT, recvcount, disps, mpi.MPI_FLOAT, main_node, mpi.MPI_COMM_WORLD)
			if( myid == 0 ): dtot[iref] = recvbuf
		del recvbuf


		#  The while loop over even angles delta should start here.
		#  prepare reference directions
		pass#IMPORTIMPORTIMPORT from utilities import even_angles, getvec
		refa = sparx_utilities.even_angles(60.0)
		numrefang = len(refa)
		refanorm = numpy.empty( (numrefang, 3), dtype = numpy.float32)
		for i in range(numrefang):
			tmp = sparx_utilities.getvec(refa[i][0], refa[i][1])
			for j in range(3):
				refanorm[i][j] = tmp[j]
		del  refa, tmp

		transv = numpy.empty( (nima, 3), dtype = numpy.float32)
		if runtype=="ASSIGNMENT":
			for im in range(nima):
				trns = data[im].get_attr( "xform.projection" )
				for j in range(3):
					transv[im][j] = trns.at(2,j)
		else:
			# For REFINEMENT we have a problem, as the exact angle is known only after the next step of assigning projections.
			# So, we will assume it is the one with max peak
			for im in range(nima):
				qt = -1.0e23
				it = -1
				for iref in range(numref):
					pt = peaks[iref][im]
					if(pt > qt):
						qt = pt
						it = iref
				for j in range(3):
					transv[im][j] = trans[it][im].at(2,j)
		#  We have all vectors, now create a list of assignments of images to references
		refassign = [-1]*nima
		for im in range(nima):
			refassign[im] = abs(numpy.inner(refanorm,transv[im])).argmax()
		assigntorefa = mpi.mpi_gatherv(refassign, nima, mpi.MPI_INT, recvcount, disps, mpi.MPI_INT, main_node, mpi.MPI_COMM_WORLD)
		assigntorefa = list(map(int, assigntorefa))

		del refassign, refanorm, transv


		"""Multiline Comment19"""


		if myid == main_node:
			SA = False
			asi = [[] for iref in range(numref)]
			report_error = 0
			for imrefa in range(numrefang):
				pass#IMPORTIMPORTIMPORT from utilities import findall
				N = sparx_utilities.findall(imrefa, assigntorefa)
				current_nima = len(N)
				if( current_nima >= numref and report_error == 0):
					tasi = [[] for iref in range(numref)]
					maxasi = current_nima//numref
					nt = current_nima
					kt = numref
					K = list(range(numref))

					d = numpy.empty( (numref, current_nima), dtype = numpy.float32)
					for ima in range(current_nima):
						for iref in range(numref):  d[iref][ima] = dtot[iref][N[ima]]

					while nt > 0 and kt > 0:
						l = d.argmax()
						group = l//current_nima
						ima   = l-current_nima*group
						if SA:
							J = [0.0]*numref
							sJ = 0
							Jc = [0.0]*numref
							for iref in range(numref):
								J[iref] = numpy.exp(d[iref][ima]/T)
								sJ += J[iref]
							for iref in range(numref):
								J[iref] /= sJ
							Jc[0] = J[0]
							for iref in range(1, numref):
								Jc[iref] = Jc[iref-1]+J[iref]
							sss = numpy.random.random()
							for group in range(numref):
								if( sss <= Jc[group]): break
						tasi[group].append(N[ima])
						N[ima] = -1
						for iref in range(numref):  d[iref][ima] = -1.e10
						nt -= 1
						masi = len(tasi[group])
						if masi == maxasi:
							for im in range(current_nima):  d[group][im] = -1.e10
							kt -= 1
					else:
						for ima in range(current_nima):
							if N[ima] > -1:
								qm = -1.e10
								for iref in range(numref):
									qt = dtot[iref][N[ima]]
									if( qt > qm ):
										qm = qt
										group = iref
								tasi[group].append(N[ima])

					del d, N, K
					if  SA:  del J, Jc
					for iref in range(numref):
						asi[iref] += tasi[iref]
					del tasi
				else:
					report_error = 1
			#  This should be deleted only once we know that the number of images is sufficiently large, see below.
			del dtot

		else:
			assignment = []
			report_error = 0

		report_error = sparx_utilities.bcast_number_to_all(report_error, source_node = main_node)
		if report_error == 1:  sparx_global_def.ERROR('Number of images within a group too small', "mref_ali3d_MPI", 1, myid)
		if myid == main_node:
			assignment = [0]*total_nima
			for iref in range(numref):
				for im in range(len(asi[iref])):
					assignment[asi[iref][im]] = iref
			del asi
		
		"""Multiline Comment20"""

		"""Multiline Comment21"""

		assignment = mpi.mpi_scatterv(assignment, recvcount, disps, mpi.MPI_INT, recvcount[myid], mpi.MPI_INT, main_node, mpi.MPI_COMM_WORLD)
		assignment = list(map(int, assignment))


		#  compute number of particles that changed assignment and how many are in which group
		nchng = 0
		npergroup = [0]*numref
		for im in range(nima):
			iref = data[im].get_attr('group')
			npergroup[assignment[im]] += 1
			if( iref != assignment[im]): nchng += 1
			data[im].set_attr('group', assignment[im])
		nchng = mpi.mpi_reduce(nchng, 1, mpi.MPI_INT, mpi.MPI_SUM, 0, mpi.MPI_COMM_WORLD)
		npergroup = mpi.mpi_reduce(npergroup, numref, mpi.MPI_INT, mpi.MPI_SUM, 0, mpi.MPI_COMM_WORLD)
		npergroup = list(map(int, npergroup))
		terminate = 0
		if( myid == 0 ):
			nchng = int(nchng[0])
			precn = 100*float(nchng)/float(total_nima)
			msg = " Number of particles that changed assignments %7d, percentage of total: %5.1f"%(nchng, precn)
			log.add(msg)
			msg = " Group       number of particles"
			log.add(msg)
			for iref in range(numref):
				msg = " %5d       %7d"%(iref+1, npergroup[iref])
				log.add(msg)
			if(precn <= termprec):  terminate = 1
		terminate = mpi.mpi_bcast(terminate, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD)
		terminate = int(terminate[0])

		if runtype=="REFINEMENT":
			for im in range(nima):
				data[im].set_attr('xform.projection', trans[assignment[im]][im])
				pixer[0][im] = pixer[assignment[im]][im]
			pixer = pixer[0]

			if(center == -1):
				cs[0], cs[1], cs[2], dummy, dummy = sparx_utilities.estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)				
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f"%(cs[0], cs[1], cs[2])
					log.add(msg)
				cs = mpi.mpi_bcast(cs, 3, mpi.MPI_FLOAT, main_node, mpi.MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				sparx_utilities.rotate_3D_shift(data, cs)
			#output pixel errors
			recvbuf = mpi.mpi_gatherv(pixer, nima, mpi.MPI_FLOAT, recvcount, disps, mpi.MPI_FLOAT, main_node, mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			if(myid == main_node):
				recvbuf = list(map(float, recvbuf))
				pass#IMPORTIMPORTIMPORT from statistics import hist_list
				lhist = 20
				region, histo = sparx_statistics.hist_list(recvbuf, lhist)
				if(region[0] < 0.0):  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles"
				log.add(msg)
				for lhx in range(lhist):
					msg = " %10.3f      %7d"%(region[lhx], histo[lhx])
					log.add(msg)
				del region, histo
			del recvbuf

		fscc = [None]*numref

		if fourvar and runtype=="REFINEMENT":
			sumvol = sparx_utilities.model_blank(nx, nx, nx)

		start_time = time.time()
		for iref in range(numref):
			#  3D stuff
			pass#IMPORTIMPORTIMPORT from time import localtime, strftime
			if(CTF): volref, fscc[iref] = sparx_reconstruction.rec3D_MPI(data, snr, sym, sparx_utilities.model_circle(last_ring, nx, nx, nx), os.path.join(outdir, "resolution_%02d_%04d"%(iref, total_iter)), myid, main_node, index = iref, npad = npad, finfo=frec)
			else:    volref, fscc[iref] = sparx_reconstruction.rec3D_MPI_noCTF(data, sym, sparx_utilities.model_circle(last_ring, nx, nx, nx), os.path.join(outdir, "resolution_%02d_%04d"%(iref, total_iter)), myid, main_node, index = iref, npad = npad, finfo=frec)
			if(myid == 0):
				log.add( "Time to compute 3D: %d" % (time.time()-start_time) );start_time = time.time()

			if(myid == main_node):
				volref.write_image(os.path.join(outdir, "vol%04d.hdf"%( total_iter)), iref)
				if fourvar and runtype=="REFINEMENT":
					sumvol += volref
			del volref

		if runtype=="REFINEMENT":
			if fourvar:
				varf = sparx_statistics.varf3d_MPI(data, os.path.join(outdir, "ssnr%04d"%total_iter), None,sumvol,last_ring, 1.0, 1, CTF, 1, sym, myid)
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

		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		if terminate ==1: # headers are only updated when the program is going to terminate
			start_time = time.time()
			if runtype=="REFINEMENT":
				par_str = ['xform.projection', 'ID', 'group']
			else:
				par_str = ['group', 'ID' ]
			if myid == main_node:
				pass#IMPORTIMPORTIMPORT from utilities import file_type
				if(sparx_utilities.file_type(stack) == "bdb"):
					pass#IMPORTIMPORTIMPORT from utilities import recv_attr_dict_bdb
					sparx_utilities.recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
				else:
					pass#IMPORTIMPORTIMPORT from utilities import recv_attr_dict
					sparx_utilities.recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
			else:		sparx_utilities.send_attr_dict(main_node, data, par_str, image_start, image_end)
			if(myid == 0):
				log.add( "Time to write headers: %d\n" % (time.time()-start_time) );start_time = time.time()
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			if myid==main_node:
				log.add("mref_ali3d_MPI terminated due to small number of objects changing assignments")
			break
	if myid==main_node:
		log.add("mref_ali3d_MPI finishes")



def Kmref_ali3d_MPI(stack, ref_vol, outdir, maskfile=None, focus = None, maxit=1, ir=1, ou=-1, rs=1, 
            xr ="4 2  2  1", yr="-1", ts="1 1 0.5 0.25",   delta="10  6  4  4", an="-1",
	      center = -1, nassign = 3, nrefine= 1, CTF = False, snr = 1.0,  ref_a="S", sym="c1",
	      user_func_name="ref_ali3d", npad = 4, debug = False, fourvar=False, termprec = 0.0, mpi_comm = None, log = None): 
	pass#IMPORTIMPORTIMPORT from utilities      import model_circle, reduce_EMData_to_root, bcast_EMData_to_all, bcast_number_to_all, drop_image
	pass#IMPORTIMPORTIMPORT from utilities      import bcast_list_to_all, get_image, get_input_from_string, get_im
	pass#IMPORTIMPORTIMPORT from utilities      import get_arb_params, set_arb_params, drop_spider_doc, send_attr_dict
	pass#IMPORTIMPORTIMPORT from utilities      import get_params_proj, set_params_proj, model_blank, write_text_row, write_text_file
	pass#IMPORTIMPORTIMPORT from filter         import filt_params, filt_btwl, filt_ctf, filt_table, fit_tanh, filt_tanl
	pass#IMPORTIMPORTIMPORT from utilities      import rotate_3D_shift,estimate_3D_center_MPI
	pass#IMPORTIMPORTIMPORT from alignment      import Numrinit, prepare_refrings, proj_ali_incore
	pass#IMPORTIMPORTIMPORT from random         import randint
	pass#IMPORTIMPORTIMPORT from filter         import filt_ctf
	pass#IMPORTIMPORTIMPORT from utilities      import print_begin_msg, print_end_msg, print_msg
	pass#IMPORTIMPORTIMPORT from projection     import prep_vol, prgs, project, prgq, gen_rings_ctf
	pass#IMPORTIMPORTIMPORT from utilities      import wrap_mpi_recv, wrap_mpi_send
	pass#IMPORTIMPORTIMPORT from copy           import deepcopy
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT from mpi            import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	pass#IMPORTIMPORTIMPORT from mpi            import mpi_reduce, MPI_INT, MPI_SUM
	
	if mpi_comm == None: mpi_comm = mpi.MPI_COMM_WORLD
	number_of_proc = mpi.mpi_comm_size(mpi_comm)
	myid           = mpi.mpi_comm_rank(mpi_comm)
	main_node = 0
	if log == None:
		pass#IMPORTIMPORTIMPORT from logger import Logger
		log =sparx_logger.Logger()

	if os.path.exists(outdir): sparx_global_def.ERROR('Output directory exists, please change the name and restart the program', "Kmref_ali3d_MPI ", 1, myid)
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

	if myid == main_node:	
		os.mkdir(outdir)
		pass#IMPORTIMPORTIMPORT import global_def
		sparx_global_def.LOGFILE =  os.path.join(outdir, sparx_global_def.LOGFILE)
		log.add("Kmref_ali3d_MPI - Traditional Kmeans clustering  !")
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

	pass#IMPORTIMPORTIMPORT from time import time

	if debug:
		pass#IMPORTIMPORTIMPORT from time import sleep
		while not os.path.exists(outdir):
			print("Node ",myid,"  waiting...")
			time.sleep(5)

		finfo = open(os.path.join(outdir, "progress%04d"%myid), 'w')
		frec  = open( os.path.join(outdir, "recons%04d"%myid), "w" )
	else:
		finfo = None
		frec  = None

	xrng        = sparx_utilities.get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = sparx_utilities.get_input_from_string(yr)
	step        = sparx_utilities.get_input_from_string(ts)
	delta       = sparx_utilities.get_input_from_string(delta)
	lstp = min( len(xrng), len(yrng), len(step), len(delta) )
	if an == "-1":
		an = []
		for i in range(len(xrng)):   an.append(-1)
	else:
		pass#IMPORTIMPORTIMPORT from  alignment	    import proj_ali_incore_local
		an      = sparx_utilities.get_input_from_string(an)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	center      = int(center)

	numref = EMAN2_cppwrap.EMUtil.get_image_count(ref_vol)
	volref     = EMAN2_cppwrap.EMData()
	volref.read_image(stack, 0)
	nx      = volref.get_xsize()
	if last_ring < 0:	last_ring = nx//2 - 2

	fscmask = sparx_utilities.model_circle(last_ring, nx, nx, nx)

	if myid == main_node:
		pass#IMPORTIMPORTIMPORT import user_functions
		user_func = sparx_user_functions.factory[user_func_name]
		log.add("Input stack                 : %s"%(stack))
		log.add("Reference volumes           : %s"%(ref_vol))	
		log.add("Number of reference volumes : %i"%(numref))
		log.add("Output directory            : %s"%(outdir))
		log.add("User function               : %s"%(user_func_name))
		log.add("Maskfile                    : %s"%(maskfile))
		log.add("Inner radius                : %i"%(first_ring))
		log.add("Outer radius                : %i"%(last_ring))
		log.add("Ring step                   : %i"%(rstep))
		log.add("X search range              : %s"%(xrng))
		log.add("Y search range              : %s"%(yrng))
		log.add("Translational step          : %s"%(step))
		log.add("Angular step                : %s"%(delta))
		log.add("Angular search range        : %s"%(an))
		log.add("Number of assignments in each iteration   : %i"%(nassign))
		log.add("Number of alignments in each iteration    : %i"%(nrefine))
		log.add("Number of iterations                      : %i"%(lstp*maxit) )
		log.add("Center type                 : %i"%(center))
		log.add("CTF correction              : %s"%(CTF))
		log.add("Signal-to-Noise Ratio       : %f"%(snr))
		log.add("Reference projection method : %s"%(ref_a))
		log.add("Symmetry group              : %s"%(sym))
		log.add("Percentage of change for termination: %f"%(termprec))
		log.add("User function               : %s"%(user_func_name))

	if maskfile:
		if type(maskfile) is bytes: mask3D = sparx_utilities.get_image(maskfile)
		else: 	                                mask3D = maskfile
	else        :  mask3D = sparx_utilities.model_circle(last_ring, nx, nx, nx)

	numr     = sparx_alignment.Numrinit(first_ring, last_ring, rstep, "F")
	mask2D   = sparx_utilities.model_circle(last_ring, nx, nx) - sparx_utilities.model_circle(first_ring, nx, nx)

	if myid == main_node:
		nima =EMAN2_cppwrap.EMUtil.get_image_count( stack )
		list_of_particles=list(range(nima))
	else:
		nima = 0

	nima = sparx_utilities.bcast_number_to_all(nima, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*nima

	list_of_particles = sparx_utilities.bcast_list_to_all(list_of_particles, myid,  source_node = main_node)

	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
	# create a list of images for each node
	total_nima = nima
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)

	if debug:
		finfo.write( "Image_start, image_end: %d %d\n" %(image_start, image_end) )
		finfo.flush()

	start_time = time.time()
	data = EMAN2_cppwrap.EMData.read_images(stack, list_of_particles)
	if myid == main_node:
		log.add( "Time to read data: %d\n" % (time.time()-start_time) );start_time = time.time()
	#  Initialize Particle ID and set group number to non-existant -1
	assignment = [-1]*len(data)
	for im in range(len(data)):
		data[im].set_attr_dict({'ID':list_of_particles[im], 'group':-1})

	if fourvar:
		pass#IMPORTIMPORTIMPORT from reconstruction import rec3D_MPI
		pass#IMPORTIMPORTIMPORT from statistics     import varf3d_MPI
		#  Compute Fourier variance
		vol, fscc = sparx_reconstruction.rec3D_MPI(data, snr, sym, fscmask, os.path.join(outdir, "resolution0000"), myid, main_node, finfo=frec, npad=npad)
		varf = sparx_statistics.varf3d_MPI(data, os.path.join(outdir, "ssnr0000"), None, vol, last_ring, 1.0, 1, CTF, 1, sym, myid)
		if myid == main_node:   
			varf = 1.0/varf
			varf.write_image( os.path.join(outdir,"varf0000.hdf") )
	else:
		varf = None

	if myid == main_node:
		for  iref in range(numref):
			sparx_utilities.get_im(ref_vol, iref).write_image(os.path.join(outdir, "volf0000.hdf"), iref)
	mpi.mpi_barrier( mpi.MPI_COMM_WORLD )

	if CTF:
		if(data[0].get_attr("ctf_applied") > 0.0):  sparx_global_def.ERROR("Kmref_ali3d_MPI does not work for CTF-applied data", "Kmref_ali3d_MPI", 1, myid)
		pass#IMPORTIMPORTIMPORT from reconstruction import rec3D_MPI
	else:
		pass#IMPORTIMPORTIMPORT from reconstruction import rec3D_MPI_noCTF

	if debug:
		finfo.write( '%d loaded  \n' % len(data) )
		finfo.flush()

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in range(number_of_proc):
		if im == main_node:  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )

	total_iter = 0
	tr_dummy = EMAN2_cppwrap.Transform({"type":"spider"})

	Niter = int(lstp*maxit*(nassign + nrefine) )
	for Iter in range(Niter):
		N_step = (Iter%(lstp*(nassign+nrefine)))/(nassign+nrefine)
		if Iter%(nassign+nrefine) < nassign:
			runtype = "ASSIGNMENT"
		else:
			runtype = "REFINEMENT"

		total_iter += 1
		if myid == main_node:
			log.add("\n%s ITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f"%(runtype, total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step]))
			start_ime = time.time()
	
		peaks = [ -1.0e23]*nima
		if runtype=="REFINEMENT":
			trans = [tr_dummy]*nima
			pixer = [0.0]*nima
			if(an[N_step] > 0):
				pass#IMPORTIMPORTIMPORT from utilities    import even_angles
				ref_angles = sparx_utilities.even_angles(delta[N_step], symmetry=sym, method = ref_a, phiEqpsi = "Zero")
				# generate list of angles
				pass#IMPORTIMPORTIMPORT from alignment import generate_list_of_reference_angles_for_search
				list_of_reference_angles = \
				sparx_alignment.generate_list_of_reference_angles_for_search(ref_angles, sym=sym)
				del ref_angles
			else:  list_of_reference_angles = [[1.0,1.0]]
 
		cs = [0.0]*3
		for iref in range(numref):
			if myid==main_node:
				volft = sparx_utilities.get_im(os.path.join(outdir, "volf%04d.hdf"%(total_iter-1)), iref)
			else:
				volft=sparx_utilities.model_blank(nx,nx,nx)
			sparx_utilities.bcast_EMData_to_all(volft, myid, main_node)
			volft, kb = sparx_projection.prep_vol(volft)

			if CTF:
				previous_defocus = -1.0
				if runtype=="REFINEMENT":
					start_time = time.time()
					prjref = sparx_projection.prgq( volft, kb, nx, delta[N_step], ref_a, sym, MPI=True)
					if myid == main_node:
						log.add( "Calculation of projections: %d" % (time.time()-start_time) );start_time = time.time()
					del volft, kb
			else:
				if runtype=="REFINEMENT":
					start_time = time.time()
					refrings = sparx_alignment.prepare_refrings( volft, kb, nx, delta[N_step], ref_a, sym, numr)
					if myid == main_node:
						log.add( "Initial time to prepare rings: %d" % (time.time()-start_time) );start_time = time.time()
					del volft, kb

			start_time = time.time()
			for im in range(nima):
				if CTF:
					ctf = data[im].get_attr( "ctf" )
					if runtype=="REFINEMENT":
						if ctf.defocus != previous_defocus:
							previous_defocus = ctf.defocus
							rstart_time = time.time()
							refrings = sparx_projection.gen_rings_ctf( prjref, nx, ctf, numr)
							if myid == main_node:
								log.add( "Repeated time to prepare rings: %d" % (time.time()-rstart_time) );rstart_time = time.time()

				if runtype=="ASSIGNMENT":
					phi,tht,psi,s2x,s2y = sparx_utilities.get_params_proj(data[im])
					ref = sparx_projection.prgs( volft, kb, [phi,tht,psi,-s2x,-s2y])
					if CTF:  ref = sparx_filter.filt_ctf( ref, ctf )
					peak = ref.cmp("ccc",data[im],{"mask":mask2D, "negative":0})
					if not(finfo is None):
						finfo.write( "ID,iref,peak: %6d %d %8.5f\n" % (list_of_particles[im],iref,peak) )
				else:
					if an[N_step] == -1:
						peak, pixel_error = sparx_alignment.proj_ali_incore(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step])
					else:
						peak, pixel_error = sparx_alignment.proj_ali_incore_local(data[im], refrings, list_of_reference_angles, numr,\
																	xrng[N_step], yrng[N_step], step[N_step], an[N_step])
					if not(finfo is None):
						phi,tht,psi,s2x,s2y = sparx_utilities.get_params_proj(data[im])
						finfo.write( "ID,iref,peak,trans: %6d %d %f %f %f %f %f %f\n"%(list_of_particles[im],iref,peak,phi,tht,psi,s2x,s2y) )
						finfo.flush()

				if peak > peaks[im]:
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
			if myid == main_node:
				log.add( "Time to process particles for reference %3d: %d" % (iref, time.time()-start_time) );start_time = time.time()


		del peaks
		if runtype=="ASSIGNMENT":  del volft, kb, ref
		else:
			if CTF: del prjref
			del refrings
			if an[N_step] > 0: del list_of_reference_angles


		#  compute number of particles that changed assignment and how man are in which group
		nchng = 0
		npergroup = [0]*numref
		for im in range(nima):
			iref = data[im].get_attr('group')
			npergroup[iref] += 1
			if iref != assignment[im]:
				assignment[im] = iref
				nchng += 1
		nchng = mpi.mpi_reduce(nchng, 1, mpi.MPI_INT, mpi.MPI_SUM, 0, mpi.MPI_COMM_WORLD)
		npergroup = mpi.mpi_reduce(npergroup, numref, mpi.MPI_INT, mpi.MPI_SUM, 0, mpi.MPI_COMM_WORLD)
		npergroup = list(map(int, npergroup))
		terminate  = 0
		empty_group =0
		if myid == main_node:
			nchng = int(nchng[0])
			precn = 100*float(nchng)/float(total_nima)
			msg = " Number of particles that changed assignments %7d, percentage of total: %5.1f"%(nchng, precn)
			log.add(msg)
			msg = " Group       number of particles"
			log.add(msg)
			for iref in range(numref):
				msg = " %5d       %7d"%(iref+1, npergroup[iref])
				log.add(msg)
				if npergroup[iref]==0:
					empty_group =1
			if precn <= termprec:  
				terminate = 1
			if empty_group ==1:
				terminate = 1
		terminate = mpi.mpi_bcast(terminate, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD)
		terminate = int(terminate[0])
		empty_group = mpi.mpi_bcast(empty_group, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD)
		empty_group = int(empty_group[0])
		if empty_group ==1: break # program stops whenever empty_group appears!
		if runtype=="REFINEMENT":
			for im in range(nima):
				data[im].set_attr('xform.projection', trans[im])


			if center == -1:
				cs[0], cs[1], cs[2], dummy, dummy = sparx_utilities.estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)				
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f"%(cs[0], cs[1], cs[2])
					log.add(msg)
				cs = mpi.mpi_bcast(cs, 3, mpi.MPI_FLOAT, main_node, mpi.MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				sparx_utilities.rotate_3D_shift(data, cs)
			#output pixel errors
			pass#IMPORTIMPORTIMPORT from mpi import mpi_gatherv
			recvbuf = mpi.mpi_gatherv(pixer, nima, mpi.MPI_FLOAT, recvcount, disps, mpi.MPI_FLOAT, main_node, mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			if myid == main_node:
				recvbuf = list(map(float, recvbuf))
				pass#IMPORTIMPORTIMPORT from statistics import hist_list
				lhist = 20
				region, histo = sparx_statistics.hist_list(recvbuf, lhist)
				if region[0] < 0.0:  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles"
				log.add(msg)
				for lhx in range(lhist):
					msg = " %10.3f      %7d"%(region[lhx], histo[lhx])
					log.add(msg)
				del region, histo
			del recvbuf


			for im in range(nima):
				phi,theta,psi,tx,ty = sparx_utilities.get_params_proj(data[im])
				trans[im] = [phi,theta,psi,tx,ty]
			if myid == main_node:
				all_trans = []
				for klm in range(number_of_proc):
					if(klm == main_node):  all_trans.extend(copy.deepcopy(trans))
					else:  all_trans.extend(sparx_utilities.wrap_mpi_recv(klm, mpi.MPI_COMM_WORLD))
			else:  sparx_utilities.wrap_mpi_send(trans, main_node, mpi.MPI_COMM_WORLD)
			if myid == main_node:
				sparx_utilities.write_text_row(all_trans, os.path.join(outdir, "params_%04d.txt"%(total_iter)) )
				del all_trans


		if myid == main_node:
			all_trans = []
			for klm in range(number_of_proc):
				if(klm == main_node):  all_trans.extend(copy.deepcopy(assignment))
				else:  all_trans.extend(sparx_utilities.wrap_mpi_recv(klm, mpi.MPI_COMM_WORLD))
		else:  sparx_utilities.wrap_mpi_send(assignment, main_node, mpi.MPI_COMM_WORLD)
		if myid == main_node:
			sparx_utilities.write_text_file(all_trans, os.path.join(outdir, "assignment_%04d.txt"%(total_iter)) )
			del all_trans

		#if CTF: del vol
		fscc = [None]*numref

		if fourvar and runtype=="REFINEMENT":
			sumvol = sparx_utilities.model_blank(nx, nx, nx)

		sart_time = time.time()
		for iref in range(numref):
			#  3D stuff
			pass#IMPORTIMPORTIMPORT from time import localtime, strftime
			if CTF: volref, fscc[iref] = sparx_reconstruction.rec3D_MPI(data, snr, sym, fscmask, os.path.join(outdir, "resolution_%02d_%04d.txt"%(iref, total_iter)), myid, main_node, index = iref, npad = npad, finfo=frec)
			else:    volref, fscc[iref] = sparx_reconstruction.rec3D_MPI_noCTF(data, sym, fscmask, os.path.join(outdir, "resolution_%02d_%04d.txt"%(iref, total_iter)), myid, main_node, index = iref, npad = npad, finfo=frec)
			if myid == main_node:
				log.add( "Time to compute 3D: %d" % (time.time()-start_time) );start_time = time.time()

			if myid == main_node:
				volref.write_image(os.path.join(outdir, "vol%04d.hdf"%( total_iter)), iref)
				if fourvar and runtype=="REFINEMENT":
					sumvol += volref
			del volref

		if runtype=="REFINEMENT":
			if fourvar:
				varf = sparx_statistics.varf3d_MPI(data, os.path.join(outdir, "ssnr%04d"%total_iter), None,sumvol,last_ring, 1.0, 1, CTF, 1, sym, myid)
				if myid == main_node:   
					varf = 1.0/varf
					varf.write_image( os.path.join(outdir,"varf%04d.hdf"%total_iter) )

		if myid == main_node:
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
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		#start_time = time()
		#if runtype=="REFINEMENT":
		#	par_str = ['xform.projection', 'ID', 'group']
		#else:
		#	par_str = ['group', 'ID' ]
	        #if myid == main_node:
		#	from utilities import file_type
	        #	if file_type(stack) == "bdb":
	        #		from utilities import recv_attr_dict_bdb
	        #		recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        #	else:
	        #		from utilities import recv_attr_dict
	        #		recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        #else:		send_attr_dict(main_node, data, par_str, image_start, image_end)
		if terminate == 1:
			if myid == main_node:
				log.add("Kmref_ali3d_MPI terminated due to small number of objects changing assignments")
			break
		#if myid == main_node:
		#	log.add( "Time to write headers: %d\n" % (time()-start_time) );start_time = time()
	######writing paritition only in the end of the program
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	"""Multiline Comment22"""
	if myid == main_node:
		log.add("Kmref_ali3d_MPI is done!")
	return empty_group


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

	pass#IMPORTIMPORTIMPORT from utilities import file_type

	oextension = sparx_utilities.file_type(ous)	

	if oextension == "bdb":
		pass#IMPORTIMPORTIMPORT from EMAN2db import db_open_dict
		DB = EMAN2db.db_open_dict(ous)

	# iterate over all images in the list, even if it's only one...
	for ins in image_list:
		nima = EMAN2_cppwrap.EMUtil.get_image_count(ins)
		data = EMAN2_cppwrap.EMData()
		iextension = sparx_utilities.file_type(ins)

		if iextension == "bdb":
			pass#IMPORTIMPORTIMPORT from EMAN2db import db_open_dict

		if nima == 1 and oextension == "spi":
			data.read_image(ins)
			data.write_image(ous, 0, EMAN2_cppwrap.EMUtil.ImageType.IMAGE_SINGLE_SPIDER)
			
		elif iextension == "bdb" and oextension == "bdb":
			
			OB = EMAN2db.db_open_dict(ins)
			for i in range(nima):
				DB[gl_index] = OB[i]
				gl_index += 1
			OB.close()

		elif iextension == "bdb":
			
			DB = EMAN2db.db_open_dict(ins)
			for i in range(nima):
				a = DB[i]
				a.write_image(ous, gl_index)
				gl_index += 1
			DB.close()
			
		elif oextension == "bdb":
			
			for i in range(nima):
				a = EMAN2_cppwrap.EMData()
				a.read_image(ins, i)
				DB[gl_index] = a
				gl_index += 1
			
		else:
			for im in range(nima):
				data.read_image(ins, im)
				data.write_image(ous, gl_index)
				gl_index += 1

	if oextension == "bdb":
		DB.close()

def project3d(volume, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False):
	pass#IMPORTIMPORTIMPORT from projection    import   prgs, prep_vol, project
	pass#IMPORTIMPORTIMPORT from utilities     import   even_angles, read_text_row, set_params_proj, model_gauss_noise, info
	pass#IMPORTIMPORTIMPORT from string        import   split
	pass#IMPORTIMPORTIMPORT from filter        import   filt_ctf,filt_gaussl
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT import types
	if trillinear:
		pass#IMPORTIMPORTIMPORT from fundamentals import fft
		pass#IMPORTIMPORTIMPORT from projection   import prgl
		pass#IMPORTIMPORTIMPORT from morphology   import ctf_img_real
	if trillinear and realsp:
		sparx_global_def.ERROR("Both relion mode and realsp mode are specified","project3d", 1)
	
	if listagls is None:
		angles = sparx_utilities.even_angles(delta, symmetry = symmetry, method = method, phiEqpsi = phiEqpsi)
	elif(type(listagls) is bytes):
		angles = sparx_utilities.read_text_row(listagls, "", "")
	else:
		angles = listagls

	# try to parse the CTFs list. this is either not set (None), a filename or a list of values
	if listctfs is None:
		# not set, so simply ignore it 
		ctfs = None
	elif (type(listctfs) is bytes):
		# a string, so assume this is a filename and try to open the file
		try:
			ctfs = sparx_utilities.read_text_row(listctfs, "", "")
		except:
			ctfs = [None for ii in range(len(angles))]
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

	if(type(volume) is bytes):
		vol = EMAN2_cppwrap.EMData()
		vol.read_image(volume)
		if(mask):
			if(type(mask) is bytes):
				maski = EMAN2_cppwrap.EMData()
				maski.read_image(volume)
				EMAN2_cppwrap.Util.mul_img(vol, maski)
				del maski
			else:
				EMAN2_cppwrap.Util.mul_img(vol, mask)
		nx = vol.get_xsize()
		ny = vol.get_ysize()
		nz = vol.get_zsize()
		
		if realsp:
			volft = vol
		elif trillinear:
			volft = sparx_projection.prep_vol(vol, npad = 2, interpolation_method = 1)
		else:
			if(nx==nz&ny==nz):
				volft, kb = sparx_projection.prep_vol(vol)
			else:
				volft, kbx,kby,kbz = sparx_projection.prep_vol(vol)
	else:
		vol = volume
		if(mask):
			vol = vol.copy()
			if(type(mask) is bytes):
				maski = EMAN2_cppwrap.EMData()
				maski.read_image(volume)
				EMAN2_cppwrap.Util.mul_img(vol, maski)
				del maski
			else:
				EMAN2_cppwrap.Util.mul_img(vol, mask)
		nx = vol.get_xsize()
		ny = vol.get_ysize()
		nz = vol.get_zsize()
		
		if realsp:
			volft = vol
		elif trillinear:
			volft = sparx_projection.prep_vol(vol, npad = 2, interpolation_method = 1)
		else:
			if(nx==nz & ny==nz):
				volft, kb = sparx_projection.prep_vol(vol)
			else:
				volft, kbx,kby,kbz = sparx_projection.prep_vol(vol)

	if(type(stack) is bytes):
		Disk = True
		os.system("rm -f  "+stack)	
	else:
		out = []
		Disk = False
	
	s2x=0
	s2y=0
	
	for i in range(len(angles)):
		if(len(angles[i]) == 3):
			if realsp:
				proj = sparx_projection.project(volft, [angles[i][0], angles[i][1], angles[i][2], 0.0, 0.0], 10*nx)
			elif trillinear:
				if ctfs is not None: proj = sparx_projection.prgl(volft, [angles[i][0], angles[i][1], angles[i][2], 0.0, 0.0], 1, False)
				else:                proj = sparx_projection.prgl(volft, [angles[i][0], angles[i][1], angles[i][2], 0.0, 0.0], 1, True)
			else:
				if(nx==nz & ny==nz):
					proj = sparx_projection.prgs(volft, kb, [angles[i][0], angles[i][1], angles[i][2], 0.0, 0.0])
				else:
					proj = sparx_projection.prgs(volft, kbz, [angles[i][0], angles[i][1], angles[i][2], 0.0, 0.0], kbx, kby)
			sparx_utilities.set_params_proj(proj, [angles[i][0], angles[i][1], angles[i][2], 0.0, 0.0])
		else:
			if realsp:
				proj = sparx_projection.project(volft, [angles[i][0], angles[i][1], angles[i][2], -angles[i][3], -angles[i][4]], 10*nx)
			elif trillinear:
				if ctfs is not None: proj = sparx_projection.prgl(volft, [angles[i][0], angles[i][1], angles[i][2], -angles[i][3], -angles[i][4]], 1, False)
				else:                proj = sparx_projection.prgl(volft, [angles[i][0], angles[i][1], angles[i][2], -angles[i][3], -angles[i][4]], 1, True)
			else:
				if(nx==nz&ny==nz):
					proj = sparx_projection.prgs(volft, kb, [angles[i][0], angles[i][1], angles[i][2], -angles[i][3], -angles[i][4]])
				else:
					proj = sparx_projection.prgs(volft, kbz, [angles[i][0], angles[i][1], angles[i][2], -angles[i][3], -angles[i][4]],kbx,kby)
			if trillinear:sparx_utilities.set_params_proj(proj, angles[i][0:5])
			else:         sparx_utilities.set_params_proj(proj, angles[i])
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# proj.set_attr_dict({'active':1})
		

		# add noise, if noise is set. this is two-fold: application of noise before
		#    ctf filtering and after it.
		if noise is not None:
			try:
				# no mask, so call w/ false
				noise_ima = sparx_utilities.model_gauss_noise(noise_level, proj.get_xsize(), proj.get_ysize())
			except:
				pass
			else:
				proj += noise_ima

		# apply ctf, if ctf option is set and if we can create a valid CTF object
		if ctfs is not None:
			try:
				pass#IMPORTIMPORTIMPORT from utilities import generate_ctf
				if(len(ctfs[i]) == 6):    ctf = sparx_utilities.generate_ctf([ctfs[i][0], ctfs[i][1], ctfs[i][2], ctfs[i][3], ctfs[i][4], ctfs[i][5]])
				elif(len(ctfs[i]) == 8):  ctf = sparx_utilities.generate_ctf([ctfs[i][0], ctfs[i][1], ctfs[i][2], ctfs[i][3], ctfs[i][4], ctfs[i][5], ctfs[i][6], ctfs[i][7]])
				else:  1.0/0.0
			except:
				# there are no ctf values, so ignore this and set no values
				sparx_global_def.ERROR("Incorrect ctf values","project3d",1)
			# setting of values worked, so apply ctf and set the header info correctly
			if trillinear:
				if ctfs is not None:
					proj.set_attr_dict({"is_complex": 0})
					EMAN2_cppwrap.Util.mulclreal(proj, sparx_morphology.ctf_img_real(proj.get_ysize(), ctf))
				proj.set_attr_dict({"padffted":1, "is_complex":1})
				proj = sparx_fundamentals.fft(proj)
			else: proj = sparx_filter.filt_ctf(proj,ctf)
			proj.set_attr( "ctf",ctf)
			proj.set_attr( "ctf_applied",0)

		# add second noise level that is not affected by CTF
		if noise is not None:
			try:
				noise_ima = sparx_utilities.model_gauss_noise(noise_level, proj.get_xsize(), proj.get_ysize())
			except:
				pass
			else:
				proj += noise_ima

		if(Disk):
			proj.write_image(stack, i)
		else: 
			out.append(proj)
	if(not Disk):  return out

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
	pass#IMPORTIMPORTIMPORT from alignment    import ali_vol_func
	pass#IMPORTIMPORTIMPORT from utilities    import get_image, model_circle, get_params3D, set_params3D
	pass#IMPORTIMPORTIMPORT from utilities    import amoeba, compose_transform3
	pass#IMPORTIMPORTIMPORT from fundamentals import rot_shift3D
	
	ref = sparx_utilities.get_image(refv)
	nx = ref.get_xsize()
	ny = ref.get_ysize()
	nz = ref.get_zsize()
	if(radius != None):    mask = sparx_utilities.model_circle(radius, nx, ny, nz)
	else:                  mask = sparx_utilities.model_circle(float(min(nx, ny, nz)//2-2), nx, ny, nz)

	#names_params = ["phi", "theta", "psi", "s3x", "s3y", "s3z", "scale"]
	phi, theta, psi, s3x, s3y, s3z, mirror, scale = sparx_utilities.get_params3D(ref)
	paramsr = [phi, theta, psi, s3x, s3y, s3z, mirror, scale]
	# print  " params of the reference volume", paramsr
	ref = sparx_fundamentals.rot_shift3D(ref, paramsr[0], paramsr[1], paramsr[2], paramsr[3], paramsr[4], paramsr[5], paramsr[7])

	e = sparx_utilities.get_image(vol)
	phi, theta, psi, s3x, s3y, s3z, mirror, scale =  sparx_utilities.get_params3D(e)
	paramsv = [phi, theta, psi, s3x, s3y, s3z, mirror, scale]
	e = sparx_fundamentals.rot_shift3D(e, phi, theta, psi, s3x, s3y, s3z, scale)
	# print  " input params ", paramsv
	params = [phi, theta, psi, s3x, s3y, s3z]
	data = [e, ref, mask, params, discrepancy]
	new_params = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

	new_params = sparx_utilities.amoeba(new_params, [ang_scale, ang_scale, ang_scale, shift_scale, shift_scale, shift_scale], sparx_alignment.ali_vol_func, 1.e-4, 1.e-4, 500, data)
	print("amoeba: func_value =",new_params[1], "iter =",new_params[2])

	cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale= sparx_utilities.compose_transform3(paramsv[0], paramsv[1], paramsv[2], paramsv[3], paramsv[4], paramsv[5], paramsv[7], new_params[0][0], new_params[0][1], new_params[0][2], new_params[0][3], new_params[0][4], new_params[0][5],1.0)
	# print  " new params ", cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale, new_params[1]
	sparx_utilities.set_params3D(e, [cphi, ctheta, cpsi, cs2x, cs2y, cs2z, 0, cscale])
	if type(vol)==type(""):
		pass#IMPORTIMPORTIMPORT from utilities import write_headers
		sparx_utilities.write_headers( vol, [e], [0])
	else:
		return e

def recons3d_n_trl_MPI_one_node(prjlist, CTF, snr, sign, npad, sym, group, niter, verbose, upweighted, compensate, chunk_id):
	pass#IMPORTIMPORTIMPORT from reconstruction import recons3d_4nn_ctf_MPI, recons3d_4nn_MPI, recons3d_4nnf_MPI
	pass#IMPORTIMPORTIMPORT from utilities      import get_im, drop_image, bcast_number_to_all, write_text_file, read_text_file, info
	pass#IMPORTIMPORTIMPORT from string         import replace
	pass#IMPORTIMPORTIMPORT from time           import time
	pass#IMPORTIMPORTIMPORT from mpi            import mpi_comm_size, mpi_comm_rank, mpi_bcast, MPI_INT, MPI_COMM_WORLD, mpi_barrier
	pass#IMPORTIMPORTIMPORT from EMAN2      import Reconstructors
	pass#IMPORTIMPORTIMPORT from fundamentals import fftip, fft
	
	myid       = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
	nproc      = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
	mpi_comm   = mpi.MPI_COMM_WORLD
	time_start = time.time()
	nnxo       = 0
	if(myid == 0): nnxo = prjlist[0].get_ysize()
	else:          nnxo = 0
	nnxo = sparx_utilities.bcast_number_to_all(nnxo, source_node = 0)
	
	if verbose==0:
		finfo = None
	else:
		infofile = "progress%04d.txt"%(myid+1)
		finfo = open( infofile, 'w' )

	nnnx = ((prjlist[0].get_ysize())*2+3)

	pass#IMPORTIMPORTIMPORT from utilities      import read_text_file, read_text_row, write_text_file, info, model_blank, get_im
	pass#IMPORTIMPORTIMPORT from fundamentals   import fft,fshift
	pass#IMPORTIMPORTIMPORT from reconstruction import insert_slices, insert_slices_pdf
	pass#IMPORTIMPORTIMPORT from utilities      import reduce_EMData_to_root, model_blank
	pass#IMPORTIMPORTIMPORT from filter         import filt_table
	pass#IMPORTIMPORTIMPORT from filter		    import filt_ctf
	# reconstruction step 
	refvol = sparx_utilities.model_blank(nnnx)
	refvol.set_attr("fudge", 1.0)
	if CTF: do_ctf = 1
	else:   do_ctf = 0
	if not (finfo is None): nimg = 0
	
	fftvol = EMAN2_cppwrap.EMData()
	weight = EMAN2_cppwrap.EMData()
	
	params = {"size":nnnx, "npad":2, "snr":1.0, "sign":1, "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = EMAN2_cppwrap.Reconstructors.get( "nn4_ctfw", params )
	r.setup()
	m = [1.0]*nnnx
	is_complex = prjlist[0].get_attr("is_complex")
	if chunk_id== -1:
		for image in prjlist:
			if not is_complex: image = sparx_fundamentals.fft(image)
			image = sparx_filter.filt_ctf(image, image.get_attr("ctf"))
			image.set_attr("padffted",1)
			image.set_attr("npad",1)
			image.set_attr("bckgnoise",m)
			if (image.get_attr("group") == group): 
				if not upweighted:  sparx_reconstruction.insert_slices_pdf(r, sparx_filter.filt_table(image, image.get_attr("bckgnoise")) )
				else: sparx_reconstruction.insert_slices_pdf(r, image)
	else:
		for image in prjlist:
			if not is_complex: image = sparx_fundamentals.fft(image)
			image =sparx_filter.filt_ctf(image, image.get_attr("ctf"))
			image.set_attr("padffted",1)
			image.set_attr("npad",1)
			image.set_attr("bckgnoise",m)
			if (image.get_attr("group") == group) and (image.get_attr("chunk_id") == chunk_id): 
				if not upweighted:  sparx_reconstruction.insert_slices_pdf(r, sparx_filter.filt_table(image, image.get_attr("bckgnoise")) )
				else: sparx_reconstruction.insert_slices_pdf(r, image)		

	if not (finfo is None): 
		finfo.write( "begin reduce\n" )
		finfo.flush()

	sparx_utilities.reduce_EMData_to_root(fftvol, myid, 0, comm=mpi_comm)
	sparx_utilities.reduce_EMData_to_root(weight, myid, 0, comm=mpi_comm)

	if not (finfo is None): 
		finfo.write( "after reduce\n" )
		finfo.flush()

	if myid == 0: dummy = r.finish(compensate)
	mpi.mpi_barrier(mpi_comm)
	dummy  = EMAN2_cppwrap.EMData()
	if myid == 0 : # post-insertion operations, done only in main_node
		fftvol = EMAN2_cppwrap.Util.shrinkfvol(fftvol,2)
		weight = EMAN2_cppwrap.Util.shrinkfvol(weight,2) # reduce size 
		if( sym != "c1" ):
			fftvol    = fftvol.symfvol(sym, -1)
			weight    = weight.symfvol(sym, -1)  # symmetrize if not asymmetric
		#fftvol = Util.divn_cbyr(fftvol, weight)
		nz     = weight.get_zsize()
		ny     = weight.get_ysize()
		nx     = weight.get_xsize()
		pass#IMPORTIMPORTIMPORT from utilities  import tabessel
		pass#IMPORTIMPORTIMPORT from morphology import notzero
		beltab = sparx_utilities.tabessel(ny, nnxo) # iterative process
		nwe    = sparx_morphology.notzero(weight)
		#Util.save_slices_on_disk(weight,"slices.hdf")
		for i in range(niter):
			cvv = EMAN2_cppwrap.Util.mulreal(nwe, weight)
			#cvv = Util.read_slice_and_multiply(nwe,weight)
			cvv = sparx_fundamentals.fft(cvv)
			EMAN2_cppwrap.Util.mul_img_tabularized(cvv, nnxo, beltab)
			cvv = sparx_fundamentals.fft(cvv)
			EMAN2_cppwrap.Util.divabs(nwe, cvv)
		pass#IMPORTIMPORTIMPORT import os
		#os.system(" rm slices.hdf")
		del  beltab
		pass#IMPORTIMPORTIMPORT from morphology   import cosinemask, threshold_outside
		pass#IMPORTIMPORTIMPORT from fundamentals import fshift, fpol
		
		nwe    = sparx_morphology.threshold_outside(nwe, 0.0, 1.0e20)
		nx     = fftvol.get_ysize()
		fftvol = sparx_fundamentals.fshift(fftvol,nx//2,nx//2,nx//2)
		EMAN2_cppwrap.Util.mulclreal(fftvol, nwe)
		fftvol = sparx_fundamentals.fft(fftvol) 
		fftvol = EMAN2_cppwrap.Util.window(fftvol, nnxo, nnxo, nnxo)
		fftvol = sparx_fundamentals.fpol(fftvol, nnxo, nnxo, nnxo, True, False)
		fftvol = sparx_morphology.cosinemask(fftvol, nnxo//2-1,5,None)
		fftvol.div_sinc(1)
		#fftvol.write_image(vol_stack)
		return fftvol
	else:
		return dummy



def pca(input_stacks, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="", MPI=False, verbose=False):
	"""
		PCA of a set of images (can be 1-2-3-D).
		input_stacks - 
		subavg       - file name containing the average of the input stack.  If None, average will not be subtracted
		mask_radius  - radius of a spherical mask, cannot be specified if maskfile provided
		nvec         - number of egeinimages to be computed
		incore       - do in-core calculations, preferable for small datasets (default False)
		shuffle      - Shuffle test (default False)
		genbuf       - generate disk buffer (default True), to use the disk buffer with data set to False
		maskfile     - name of the mask file 
	"""
	pass#IMPORTIMPORTIMPORT from utilities import get_image, get_im, model_circle, model_blank
	pass#IMPORTIMPORTIMPORT from statistics import pcanalyzer
	pass#IMPORTIMPORTIMPORT import types

	if type(input_stacks[0]) is bytes: data_on_disk = True	 # input_stacks is a file name
	else:
		data_on_disk = False # input_stacks is a list of images not a file name
		if MPI:
			sparx_global_def.ERROR('MPI version for data in memory version is not implemented', "pca", 1)

	if mask_radius > 0 and maskfile !="":
		sparx_global_def.ERROR('Error: mask radius and mask file cannot be used at the same time', "pca", 1)

	if mask_radius >0:

		if(verbose): print("Using spherical mask, rad=", mask_radius)

		if maskfile!="":   sparx_global_def.ERROR('mask radius and mask file cannot be used at the same time', "pca", 1)
		if data_on_disk:
			data = sparx_utilities.get_im( input_stacks[0] )
		else:
			data = input_stacks[0]
		mask = sparx_utilities.model_circle(mask_radius, data.get_xsize(), data.get_ysize(), data.get_zsize())

	elif(maskfile!="") :
		if(verbose): print("Using mask: ", maskfile)
		mask = sparx_utilities.get_image( maskfile )
	else:
		data = EMAN2_cppwrap.EMData()
		if data_on_disk:
			data.read_image( input_stacks[0], 0, True)
		else:
			data = input_stacks[0]
		mask = sparx_utilities.model_blank(data.get_xsize(), data.get_ysize(), data.get_zsize(), bckg=1.0)

	pca = sparx_statistics.pcanalyzer(mask, nvec, incore, MPI)

	if subavg != "":
		if(verbose): print("Subtracting ", subavg, " from each image")
		avg = sparx_utilities.get_image( subavg )
		pca.setavg( avg )

	if data_on_disk:
		files = file_set( input_stacks )
	if MPI:
		pass#IMPORTIMPORTIMPORT from mpi import mpi_comm_rank, mpi_comm_size, MPI_COMM_WORLD
		myid = mpi.mpi_comm_rank( mpi.MPI_COMM_WORLD )
		ncpu = mpi.mpi_comm_size( mpi.MPI_COMM_WORLD )
	else:
		myid = 0
		ncpu = 1

	if genbuf:
		if shuffle: sparx_global_def.ERROR('Shuffle works only with usebuf', "pca", 1)

		if data_on_disk:
			bgn,end = MPI_start_end( files.nimg(), ncpu, myid )
		else:
			bgn,end = MPI_start_end( len(input_stacks), ncpu, myid )
		for i in range(bgn,end):
			if data_on_disk:
				fname, imgid = files.get( i )
				data = sparx_utilities.get_im( fname, imgid)
				if(verbose):  print("Inserting image %s, %4d" % (fname, imgid))
			else:
				data = input_stacks[i]
			pca.insert( data )

	else:
		pca.usebuf( )
		if(verbose):  print(myid, "using existing buff, nimg: ", pca.nimg)
		if shuffle:
			pca.shuffle()

	vecs = pca.analyze()
	return vecs
	
def prepare_2d_forPCA(data, mode = "a", output_stack = None, CTF = False):
	"""
		Prepare 2D images for PCA
		Average of all images is calculated using header alignment information, 
		  subtracted from each image and the difference is written to the output_stack
		If CTF, the average is calculated as
		   Av = sum(CTF_k*Im_k)/sum(CTF_k^2)
		and the difference as
		   CTF_k(Im_k - CTF_k*Av)/sum(CTF_k^2)
		average outside of a circle r = nx//2-1 is subtracted from each image
	"""
	pass#IMPORTIMPORTIMPORT from utilities    import model_blank, model_circle, set_params2D, get_params2D
	pass#IMPORTIMPORTIMPORT from fundamentals import rot_shift2D
	dopa = True
	if type(data) == type(""):
		inmem = False
		pass#IMPORTIMPORTIMPORT from utilities    import get_im
	else:
		inmem = True

	if inmem:
		n = len(data)
	else:
		n = EMAN2_cppwrap.EMUtil.get_image_count(data)

	if inmem:
		img = data[0]
	else:
		img = sparx_utilities.get_im(data,0)

	nx = img.get_xsize()
	ny = img.get_ysize()
	nz = img.get_zsize()
	
	if( output_stack == None):  outstack = [None]*n

	mask = sparx_utilities.model_circle( nx//2-2, nx, ny)
	if  CTF:
		if(img.get_attr_default('ctf_applied', 0) > 0):
			sparx_global_def.ERROR("data cannot be ctf-applied","prepare_2d_forPCA",1)
		pass#IMPORTIMPORTIMPORT from fundamentals import fft, fftip, window2d
		pass#IMPORTIMPORTIMPORT from morphology   import ctf_img
		pass#IMPORTIMPORTIMPORT from filter 	  import filt_ctf
		pass#IMPORTIMPORTIMPORT from utilities    import pad

		nx2 = 2*nx
		ny2 = 2*ny
		ave       = EMAN2_cppwrap.EMData(nx2, ny2, 1, False)
		ctf_2_sum = EMAN2_cppwrap.EMData(nx2, ny2, 1, False)

		for i in range(n):
			if inmem:
				img = data[i].copy()
			else:
				img = sparx_utilities.get_im(data, i)
			ctf_params = img.get_attr("ctf")
			if (mode == 'a'):
				angle, sx, sy, mirror, scale = sparx_utilities.get_params2D(img)
				img = sparx_fundamentals.rot_shift2D(img, angle, sx, sy, mirror, scale)
			st = EMAN2_cppwrap.Util.infomask(img, mask, False)
			img -= st[0]
			img = sparx_utilities.pad(img, nx2,ny2, 1, background = "circumference")
			sparx_fundamentals.fftip(img)
			EMAN2_cppwrap.Util.add_img(ave, sparx_filter.filt_ctf(img, ctf_params))
			EMAN2_cppwrap.Util.add_img2(ctf_2_sum, sparx_morphology.ctf_img(nx2, ctf_params))
		EMAN2_cppwrap.Util.div_filter(ave, ctf_2_sum)
		for i in range(n):
			if inmem:
				img = data[i].copy()
			else:
				img = sparx_utilities.get_im(data, i)
			ctf_params = img.get_attr("ctf")
			if (mode == 'a'):
				angle, sx, sy, mirror, scale = sparx_utilities.get_params2D(img)
				img = sparx_fundamentals.rot_shift2D(img, angle, sx, sy, mirror, scale)
			st = EMAN2_cppwrap.Util.infomask(img, mask, False)
			img -= st[0]
			img = sparx_utilities.pad(img, nx2,ny2, 1, background = "circumference")
			sparx_fundamentals.fftip(img)
			img = sparx_filter.filt_ctf(img-sparx_filter.filt_ctf(ave, ctf_params, dopa), ctf_params, dopa)
			EMAN2_cppwrap.Util.div_filter(img, ctf_2_sum)
			img = sparx_fundamentals.window2d(sparx_fundamentals.fft(img),nx,ny)
			sparx_utilities.set_params2D(img, [0.0,0.0,0.0,0,1.0])
			if( output_stack == None):  outstack[i] = img
			else:                       img.write_image(output_stack, i)
	else:
		ave  = sparx_utilities.model_blank( nx, ny)
		for i in range(n):
			if inmem:
				img = data[i].copy()
			else:
				img = sparx_utilities.get_im(data, i)
			angle, sx, sy, mirror, scale = sparx_utilities.get_params2D(img)
			img = sparx_fundamentals.rot_shift2D(img, angle, sx, sy, mirror, scale)
			st = EMAN2_cppwrap.Util.infomask(img, mask, False)
			img -= st[0]
			EMAN2_cppwrap.Util.add_img(ave, img)
		ave /= n
		for i in range(n):
			if inmem:
				img = data[i].copy()
			else:
				img = sparx_utilities.get_im(data, i)
			angle, sx, sy, mirror, scale = sparx_utilities.get_params2D(img)
			img = sparx_fundamentals.rot_shift2D(img, angle, sx, sy, mirror, scale)
			st = EMAN2_cppwrap.Util.infomask(img, mask, False)
			img -= st[0]
			EMAN2_cppwrap.Util.sub_img(img, ave)
			sparx_utilities.set_params2D(img, [0.0,0.0,0.0,0,1.0])
			if( output_stack == None):  outstack[i] = img
			else:                       img.write_image(output_stack, i)
	if( output_stack == None):  return ave, outstack
	else:                       return None

def extract_value( s ):
	pass#IMPORTIMPORTIMPORT from string import atoi, atof

	try:
		i = string.atoi( s )
		return i
	except:
		pass

	try:
		f = string.atof( s )
		return f
	except:
		pass

	return s 

def header(stack, params, zero=False, one=False, set = 0.0, randomize=False, rand_alpha=False, fimport=None, 
	   fexport=None, fprint=False, backup=False, suffix='_backup', restore=False, delete=False, consecutive=False):
	pass#IMPORTIMPORTIMPORT from string    import split
	pass#IMPORTIMPORTIMPORT from utilities import write_header, file_type, generate_ctf
	pass#IMPORTIMPORTIMPORT from random    import random, randint
	pass#IMPORTIMPORTIMPORT from utilities import set_params2D, get_params2D, set_params3D, get_params3D, set_params_proj, get_params_proj, set_ctf, get_ctf
	pass#IMPORTIMPORTIMPORT from EMAN2 import Vec2f

	if set == 0.0: doset = False
	else:          doset = True

	op = zero+one++consecutive+randomize+rand_alpha+(fimport!=None)+(fexport!=None)+fprint+backup+restore+delete+doset
	if op == 0:
		print("Error: no operation selected!")
		return
	elif op > 1:
		print("Error: more than one operation at the same time!")
		return


	params = string.split(params)

	if fimport != None: fimp = open(fimport, 'r')
	if fexport != None: fexp = open(fexport, 'w')

	nimage = EMAN2_cppwrap.EMUtil.get_image_count(stack)
	ext = sparx_utilities.file_type(stack)
	if ext == "bdb":
		pass#IMPORTIMPORTIMPORT from EMAN2db import db_open_dict
		DB = EMAN2db.db_open_dict(stack)
	for i in range(nimage):
		if fimport != None:
			line = fimp.readline()
			if len(line)==0 :
				print("Error: file " + fimport + " has only " + str(i) + " lines, while there are " + str(nimage) + " images in the file.")
				return
			parmvalues = string.split(line)
			il=0
			for p in params:
				if p[:13] == "xform.align2d":
					if len(parmvalues) < il+3:
						print("Not enough parameters!")
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
					t = EMAN2_cppwrap.Transform({"type":"2D","alpha":alpha,"tx":sx,"ty":sy,"mirror":mirror,"scale":scale})
					if ext == "bdb":
						DB.set_attr(i, "xform.align2d", t)
					elif ext == "hdf":
						EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, "xform.align2d", t, i)
					il+=5

				elif p[:16] == "xform.projection":
					if len(parmvalues) < il+3:
						print("Not enough parameters!")
						return
					phi = extract_value(parmvalues[il])
					theta = extract_value(parmvalues[il+1])
					psi = extract_value(parmvalues[il+2])
					if len(parmvalues) > il+3:
						s2x = extract_value(parmvalues[il+3])
					else:
						s2x = 0.0
					if len(parmvalues) > il+4:
						s2y = extract_value(parmvalues[il+4])
					else:
						s2y = 0.0
					#set_params_proj(img, [phi, theta, psi, s2x, s2y], params[0])
					t = EMAN2_cppwrap.Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
					t.set_trans(EMAN2_cppwrap.Vec2f(-s2x, -s2y))
					if ext == "bdb":
						DB.set_attr(i, "xform.projection", t)
					elif ext == "hdf":
						EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, "xform.projection", t, i)
					il = min(il+5, len(parmvalues))
				elif p[:13] == "xform.align3d":
					if len(parmvalues) < il+8:
						print("Not enough parameters!")
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
					t = EMAN2_cppwrap.Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi,"tx":s3x,"ty":s3y,"tz":s3z,"mirror":mirror,"scale":scale})
					if ext == "bdb":
						DB.set_attr(i, "xform.align3d", t)
					elif ext == "hdf":
						EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, "xform.align3d", t, i)
					il+=8	
				elif p == "ctf":
					if len(parmvalues) < il+6:
						print("Not enough parameters!")
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
					ctf = sparx_utilities.generate_ctf([defocus, cs, voltage, apix, bfactor, ampcont, dfdiff, dfang]) 
					if ext == "bdb":
						DB.set_attr(i, "ctf", ctf)
					elif ext == "hdf":
						EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, "ctf", ctf, i)
					il+=6	
				else:
					#if len(params)!=len(parmvalues):
						#print "Error: %d params need to be set, while %d values are provided in line %d of file." % ( len(params), len(parmvalues), i )
						#return
					if ext == "bdb":
						DB.set_attr(i, p, extract_value(parmvalues[il]))
					elif ext == "hdf":
						EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, p, extract_value(parmvalues[il]), i)
					il+=1

		else:
			for p in params:

				if zero:
					if p[:13] == "xform.align2d":
						#set_params2D(img, [0.0, 0.0, 0.0, 0, 1.0], p)
						t = EMAN2_cppwrap.Transform({"type":"2D","alpha":0.0,"tx":0.0,"ty":0.0,"mirror":0,"scale":1.0})
						if ext == "bdb":
							DB.set_attr(i, "xform.align2d", t)
						elif ext == "hdf":
							EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, "xform.align2d", t, i)	
					elif p[:16] == "xform.projection":
						#set_params_proj(img, [0.0, 0.0, 0.0, 0.0, 0.0], p)
						t = EMAN2_cppwrap.Transform({"type":"spider"})
						if ext == "bdb":
							DB.set_attr(i, "xform.projection", t)
						elif ext == "hdf":
							EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, "xform.projection", t, i)	
					elif p[:13] == "xform.align3d":
						#set_params3D(img, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 1.0], p)
						t = EMAN2_cppwrap.Transform({"type":"spider"})
						if ext == "bdb":
							DB.set_attr(i, "xform.align3d", t)
						elif ext == "hdf":
							EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, "xform.align3d", t, i)
					elif p == "ctf":
						print("Invalid operation!")
						return
					else:
						#img.set_attr(p, 0)
						if ext == "bdb":
							DB.set_attr(i, p, 0)
						elif ext == "hdf":
							EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack,p,0,i)
				elif one:
					if p[:6] == "xform." or p == "ctf":
						print("Invalid operation!")
						return
					else:
						#img.set_attr(p, 1)
						if ext == "bdb":
							DB.set_attr(i, p, 1)
						elif ext == "hdf":
							EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, p, 1, i)
				elif doset:
					if p[:6] == "xform." or p == "ctf":
						print("Invalid operation!")
						return
					else:
						#img.set_attr(p, 1)
						if ext == "bdb":
							DB.set_attr(i, p, set)
						elif ext == "hdf":
							EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, p, set, i)
				elif consecutive:
					if ext == "bdb":
						DB.set_attr(i, p, i)
					elif ext == "hdf":
						EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, p, i, i)					
				elif randomize:
					if p[:13] == "xform.align2d":
						alpha = numpy.random.random()*360.0
						sx = numpy.random.random()*2.0-1.0
						sy = numpy.random.random()*2.0-1.0
						mirror = random.randint(0, 1)
						scale = 1.0
						#set_params2D(img, [alpha, sx, sy, mirror, scale], p)
						t = EMAN2_cppwrap.Transform({"type":"2D","alpha":alpha,"tx":sx,"ty":sy,"mirror":mirror,"scale":scale})
						if ext == "bdb":
							DB.set_attr(i, "xform.align2d", t)
						elif ext == "hdf":
							EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, "xform.align2d", t, i)
					elif p[:16] == "xform.projection":
						phi = numpy.random.random()*360.0
						theta = numpy.random.random()*180.0
						psi = numpy.random.random()*360.0
						s2x = numpy.random.random()*4.0-2.0
						s2y = numpy.random.random()*4.0-2.0
						#set_params_proj(img, [phi, theta, psi, s2x, s2y], p)
						t = EMAN2_cppwrap.Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
						t.set_trans(EMAN2_cppwrap.Vec2f(-s2x, -s2y))
						if ext == "bdb":
							DB.set_attr(i, "xform.projection", t)
						elif ext == "hdf":
							EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, "xform.projection", t, i)
					elif p[:13] == "xform.align3d":
						phi = numpy.random.random()*360.0
						theta = numpy.random.random()*180.0
						psi = numpy.random.random()*360.0
						s3x = numpy.random.random()*4.0-2.0
						s3y = numpy.random.random()*4.0-2.0
						s3z = numpy.random.random()*4.0-2.0
						mirror = random.randint(0, 1)
						scale = 1.0
						#set_params3D(img, [phi, theta, psi, s3x, s3y, s3z, mirror, scale], p)	
						t = EMAN2_cppwrap.Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi,"tx":s3x,"ty":s3y,"tz":s3z,"mirror":mirror,"scale":scale})
						if ext == "bdb":
							DB.set_attr(i, "xform.align3d", t)
						elif ext == "hdf":
							EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, "xform.align3d", t, i)					
					else:
						print("Invalid operation!")
						return						
				elif rand_alpha:
					if p[:13] == "xform.align2d":
						alpha = numpy.random.random()*360.0
						sx = 0.0
						sy = 0.0
						mirror = random.randint(0, 1)
						scale = 1.0
						#set_params2D(img, [alpha, sx, sy, mirror, scale], p)
						t = EMAN2_cppwrap.Transform({"type":"2D","alpha":alpha,"tx":sx,"ty":sy,"mirror":mirror,"scale":scale})
						if ext == "bdb":
							DB.set_attr(i, "xform.align2d",t)
						elif ext == "hdf":
							EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, "xform.align2d", t, i)	
					elif p[:16] == "xform.projection":
						phi = numpy.random.random()*360.0
						theta = numpy.random.random()*180.0
						psi = numpy.random.random()*360.0
						s2x = 0.0
						s2y = 0.0
						#set_params_proj(img, [phi, theta, psi, s2x, s2y], p)
					
						t = EMAN2_cppwrap.Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
						t.set_trans(EMAN2_cppwrap.Vec2f(-s2x, -s2y))
						if ext == "bdb":
							DB.set_attr(i, "xform.projection",t)
						elif ext == "hdf":
							EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, "xform.projection", t, i)
					elif p[:13] == "xform.align3d":
						phi = numpy.random.random()*360.0
						theta = numpy.random.random()*180.0
						psi = numpy.random.random()*360.0
						s3x = 0.0
						s3y = 0.0
						s3z = 0.0
						mirror = random.randint(0, 1)
						scale = 1.0
						#set_params3D(img, [phi, theta, psi, s3x, s3y, s3z, mirror, scale], p)	
						t = EMAN2_cppwrap.Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi,"tx":s3x,"ty":s3y,"tz":s3z,"mirror":mirror,"scale":scale})
						if ext == "bdb":
							DB.set_attr(i, "xform.align3d",t)
						elif ext == "hdf":
							EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, "xform.align3d", t, i)					
					else:
						print("Invalid operation!")
						return	
											
				elif fexport != None:
					if p[:13] == "xform.align2d":
						#alpha, sx, sy, mirror, scale = get_params2D(img, p)
						if ext == "bdb":
							t = DB.get_attr(i, "xform.align2d")
							d = t.get_params("2D")
							fexp.write("%15.5f %15.5f %15.5f %10d %10.3f "%(d["alpha"],d["tx"],d["ty"],d["mirror"],d["scale"]))
														
						elif ext == "hdf":
							t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, "xform.align2d",i)
							d = t.get_params("2D")
							fexp.write("%15.5f %15.5f %15.5f %10d %10.3f "%(d["alpha"],d["tx"],d["ty"],d["mirror"],d["scale"]))
		
					elif p[:16] == "xform.projection":
						#phi, theta, psi, s2x, s2y = get_params_proj(img, p)
						if ext == "bdb":
							t = DB.get_attr(i,"xform.projection")
							d = t.get_params("spider")
							fexp.write("%15.5f %15.5f %15.5f %15.5f %15.5f "%(d["phi"],d["theta"],d["psi"],-d["tx"],-d["ty"]))
							
						elif ext == "hdf":
							t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, "xform.projection",i)
							d = t.get_params("spider")
							fexp.write("%15.5f %15.5f %15.5f %15.5f %15.5f "%(d["phi"],d["theta"],d["psi"],-d["tx"],-d["ty"]))
							
					elif p[:13] == "xform.align3d":
						#phi, theta, psi, s3x, s3y, s3z, mirror, scale = get_params3D(img, p)
						if ext == "bdb":
							t = DB.get_attr(i,"xform.align3d")
							d = t.get_params("spider")
							fexp.write("%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %10d %10.3f "%(d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"],d["mirror"],d["scale"]))
							
						elif ext == "hdf":
							t =EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, "xform.align3d",i)
							d = t.get_params("spider")	
							fexp.write("%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %10d %10.3f "%(d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"],d["mirror"],d["scale"]))
							
					elif p == "ctf":
						#defocus, cs, voltage, apix, bfactor, ampcont = get_ctf(img)
						if ext == "bdb":
							t = DB.get_attr(i,"ctf")							
						elif ext == "hdf":
							t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, "ctf",i)
						fexp.write("%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f"%(t.defocus, t.cs, t.voltage, t.apix, t.bfactor, t.ampcont, t.dfdiff, t.dfang))

					else:
						if ext == "bdb":
							fexp.write("%15s "%str(DB.get_attr(i, p)))

						elif ext == "hdf":
							fexp.write("%15s "%str(EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, p, i)))
				elif fprint:
					if p[:13] == "xform.align2d":
						#alpha, sx, sy, mirror, scale = get_params2D(img, p)
						if ext == "bdb":
							t = DB.get_attr(i,"xform.align2d")
							d = t.get_params("2D")
							print("%15.5f %15.5f %15.5f %10d %10.3f"%(d["alpha"],d["tx"],d["ty"],d["mirror"],d["scale"]), end=' ')
							
						elif ext == "hdf":
							t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, "xform.align2d",i)
							d = t.get_params("2D")
							print("%15.5f %15.5f %15.5f %10d %10.3f"%(d["alpha"],d["tx"],d["ty"],d["mirror"],d["scale"]), end=' ')
					elif p[:16] == "xform.projection":
						#phi, theta, psi, s2x, s2y = get_params_proj(img, p)
						if ext == "bdb":
							t = DB.get_attr(i,"xform.projection")
							d = t.get_params("spider")
							print("%15.5f %15.5f %15.5f %15.5f %15.5f"%(d["phi"],d["theta"],d["psi"],-d["tx"],-d["ty"]), end=' ')
						elif ext == "hdf":
							t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, "xform.projection", i)
							d = t.get_params("spider")
							print("%15.5f %15.5f %15.5f %15.5f %15.5f"%(d["phi"],d["theta"],d["psi"],-d["tx"],-d["ty"]), end=' ')
					elif p[:13] == "xform.align3d":
						#phi, theta, psi, s3x, s3y, s3z, mirror, scale = get_params3D(img, p)
						if ext == "bdb":
							t = DB.get_attr(i, "xform.align3d")
							d = t.get_params("spider")
							print("%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %10d %10.3f"%(d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"],d["mirror"],d["scale"]), end=' ')
						elif ext == "hdf":
							t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, "xform.align3d", i)
							d = t.get_params("spider")
							print("%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %10d %10.3f"%(d["phi"],d["theta"],d["psi"],d["tx"],d["ty"],d["tz"],d["mirror"],d["scale"]), end=' ')
					elif p == "ctf":
						#defocus, cs, voltage, apix, bfactor, ampcont = get_ctf(img)
						if ext == "bdb":
							t = DB.get_attr(i, "ctf")
						elif ext == "hdf":
							t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack,"ctf", i)
						print("%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f"%(t.defocus, t.cs, t.voltage, t.apix, t.bfactor, t.ampcont, t.dfdiff, t.dfang), end=' ')

					else:
						if ext == "bdb":
							print("%15s"%str(DB.get_attr(i, p)), end=' ')
						elif ext == "hdf":
							print("%15s"%str(EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, p, i)), end=' ')
				elif backup:
					#t = img.get_attr(p)
					#img.set_attr(p+suffix, t)
					if ext == "bdb":
						t= DB.get_attr(i, p)
						DB.set_attr(i, p+suffix, t)
					elif ext == "hdf":
						t= EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, p, i)
						EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack,p+suffix, t, i)
				
				elif restore:
					if p == "xform.align2d" or p == "xform.align3d" or p == "xform.projection":
						print("ERROR, no suffix in xform!")
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
							t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, p, i)
							EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, "xform.align2d", t, i)
					elif p[:16] == "xform.projection":
						#img.set_attr(p[:10], t)
						if ext == "bdb":
							t = DB.get_attr(i, p)
							DB.set_attr(i, "xform.projection", t)
						elif ext == "hdf":
							t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, p, i)
							EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, "xform.projection", t, i)
					elif p[:13] == "xform.align3d":
						#img.set_attr(p[:13], t)
						if ext == "bdb":
							t = DB.get_attr(i, p)
							DB.set_attr(i, "xform.align3d", t)
					elif ext == "hdf":
						for i in range(nimage):
							t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, p, i)
							EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, "xform.align3d", t, i)
					else:
						#img.set_attr(p[:-len(suffix)], t)
						if ext == "bdb":
							t = DB.get_attr(i, p)
							DB.set_attr(i, p[:-len(suffix)],t)
						elif ext == "hdf":
							t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, p, i)
							EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, p[:-len(suffix)], t, i)
				elif delete:
					img = EMAN2_cppwrap.EMData()
					img.read_image(stack, i, True)
					img.del_attr(p)
					sparx_utilities.write_header(stack, img, i)
		#if zero or one or randomize or rand_alpha or backup or restore or delete:
			#write_header(stack, img, i)
			if fexport != None:
				fexp.write( "\n" )
			if fprint:
				print(" ")
	if ext == "bdb": DB.close()

def MPI_start_end(nima, nproc, myid):
	image_start = int(round(float(nima)/nproc*myid))
	image_end   = int(round(float(nima)/nproc*(myid+1)))
	return image_start, image_end

class file_set(object) :

	def __init__( self, files ):
		nfile = len(files)
		self.files = files
		self.fends = [None] * nfile

		totimg = 0
		for i in range(nfile):
			totimg += EMAN2_cppwrap.EMUtil.get_image_count( self.files[i] )
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

def refvol( vollist, fsclist, output, mask ):
	pass#IMPORTIMPORTIMPORT from utilities     import get_image, read_fsc
	pass#IMPORTIMPORTIMPORT from fundamentals  import rops_table
	pass#IMPORTIMPORTIMPORT from math          import sqrt
	pass#IMPORTIMPORTIMPORT from filter        import filt_tanl, fit_tanh, filt_table, filt_vols
	pass#IMPORTIMPORTIMPORT from morphology    import threshold

	nvol = len(vollist)
	assert len(fsclist)==nvol

	fscs = [None]*nvol
	vols = [None]*nvol
	for i in range(nvol):
		fscs[i] = sparx_utilities.read_fsc( fsclist[i] )
		vols[i] = sparx_utilities.get_image( vollist[i] )
		print('rawvol, resolution: ', vollist[i], fsclist[i])

	m    = sparx_utilities.get_image( mask )
	volfs = sparx_filter.filt_vols( vols, fscs, m )

	for i in range(nvol):
		volfs[i].write_image( output, i )

# -- K-means main ---------------------------------------------------------------------------
# K-means main driver
def within_group_refinement(data, maskfile, randomize, ir, ou, rs, xrng, yrng, step, dst, maxit, FH, FF, method = "", CTF = False):

	# Comment by Zhengfan Yang 03/11/11
	# This is a simple version of ali2d_data (down to the bone), no output dir, no logfile, no MPI or CUDA, no Fourvar, no auto stop, no user function
	#  Added CTF fot simple version PAP 03/30/2017

	pass#IMPORTIMPORTIMPORT from alignment    import Numrinit, ringwe, ali2d_single_iter
	pass#IMPORTIMPORTIMPORT from filter	      import filt_tanl
	pass#IMPORTIMPORTIMPORT from fundamentals import fshift, fft
	pass#IMPORTIMPORTIMPORT from random	      import randint, random
	pass#IMPORTIMPORTIMPORT from statistics   import ave_series
	pass#IMPORTIMPORTIMPORT from utilities    import get_input_from_string, model_circle, set_params2D, get_params2D, combine_params2, inverse_transform2

	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit);
	nima = len(data)
	nx = data[0].get_xsize()
	if last_ring == -1:  last_ring = nx/2-2
	if maskfile: mask = maskfile
	else: mask = sparx_utilities.model_circle(last_ring, nx, nx)

	if randomize :
		for im in data:
			alpha, sx, sy, mirror, scale = sparx_utilities.get_params2D(im)
			alphai, sxi, syi, mirrori    = sparx_utilities.inverse_transform2(alpha, sx, sy)
			alphan, sxn, syn, mirrorn    = sparx_utilities.combine_params2(0.0, -sxi, -syi, 0, numpy.random.random()*360.0, 0.0, 0.0, random.randint(0, 1))
			#alphan, sxn, syn, mirrorn    = combine_params2(0.0, -sxi+randint(-xrng[0],xrng[0]), -syi+randint(-xrng[0],xrng[0]), 0, random()*360.0, 0, 0, randint(0, 1))
			sparx_utilities.set_params2D(im, [alphan, sxn, syn, mirrorn, 1.0])


	cnx = nx/2+1
	cny = cnx
	mode = "F"
	numr = sparx_alignment.Numrinit(first_ring, last_ring, rstep, mode)
	wr = sparx_alignment.ringwe(numr, mode)

	sx_sum = 0.0
	sy_sum = 0.0
	cs = [0.0]*2
	total_iter = 0
	if(method == "SHC"):
		#  This is my failed attempt to use SHC for 2D alignment.  
		#    Inexplicably, it did not do all that well.  While initially it converges fast
		#     and generally yields a very good solution, it converges to cluster of images scattered
		#     around the 'best' solution, i.e., the shifts are within a fraction of a pixel of what they
		#     should be and, as a result, some are in wrong positions and overall pixel error is large.
		#     Overall, Yang's method works much better, so I am leaving it at that.  PAP 01/22/2015
		for im in data:  im.set_attr('previousmax', -1.0e23)
		tavg = sparx_statistics.ave_series(data)
		for N_step in range(len(xrng)):
			nope = 0
			Iter = 0
			while(nope < len(data)//1 and Iter < max_iter ):
				total_iter += 1
				Iter += 1
				if( FH > 0.0):
					tavg = sparx_filter.filt_tanl(sparx_fundamentals.fft(tavg), FH, FF)
					if( xrng[0] > 0.0 ): cs[0] = sx_sum/float(nima)
					if( yrng[0] > 0.0 ): cs[1] = sy_sum/float(nima)
					tavg = sparx_fundamentals.fft(sparx_fundamentals.fshift(tavg, -cs[0], -cs[1]))
				else:
					tavg = sparx_filter.filt_tanl(tavg, FH, FF)
				sx_sum, sy_sum, nope = sparx_alignment.ali2d_single_iter(data, numr, wr, cs, tavg, cnx, cny, \
															xrng[N_step], yrng[N_step], step[N_step], \
															mode=mode, CTF=False, random_method = method)
				#print  "  iteration  shc   %03d   %03d   %7.2f    %7.2f  "%(total_iter,nope,cs[0],cs[1])
				#print total_iter,nope
				#for i in data:  print "  ",i.get_attr('previousmax'),
				#print "  "
				#tavg.write_image('tata.hdf',total_iter-1)
				tavg = sparx_statistics.ave_series(data)
		"""Multiline Comment49"""


	elif( method == "PCP"):
		pass#IMPORTIMPORTIMPORT from alignent import prepref
		pass#IMPORTIMPORTIMPORT from utilities import model_circle
		stp = step[-1]
		rings = sparx_alignment.prepref(data, sparx_utilities.model_circle(nx//2-1,nx,nx), cnx, cnx, numr, mode, xrng[0], xrng[0], stp)
		print(" rings  ",len(rings))
		for im in range(len(data)):
			rings[im][0][0].set_attr("inx",nx)
		tavg = sparx_statistics.ave_series(data)
		for N_step in range(len(xrng)):
			print(" xrng ",xrng[N_step])
			for Iter in range(max_iter):
				total_iter += 1
				if( FH > 0.0):
					fl = 0.1+(FH-0.1)*Iter/float(max_iter-1)
					tavg = sparx_filter.filt_tanl(tavg, fl, FF)
					"""Multiline Comment50"""
				else:
					"""Multiline Comment51"""
				cs = [0,0]
				#print  "  iteration  std   %03d   %7.2f    %7.2f  "%(total_iter,cs[0],cs[1])
				if Iter%4 != 0 or total_iter > max_iter*len(xrng)-10: delta = 0.0
				else:                                                 delta = dst
				sx_sum, sy_sum, nope = sparx_alignment.ali2d_single_iter(rings, numr, wr, cs, tavg, cnx, cny, \
															xrng[N_step], yrng[N_step], step[N_step], \
															mode=mode, CTF=False, delta=delta, random_method = method)
				for im in range(len(data)):
					alpha, tx, ty, mirror, scale = sparx_utilities.get_params2D(rings[im][0][0])
					sparx_utilities.set_params2D(data[im],[alpha, tx, ty, mirror, scale])
				tavg = sparx_statistics.ave_series(data)
				#print  "tata ",total_iter-1,Util.infomask(tavg,None,True)
				#tavg.write_image('tata.hdf',total_iter-1)
	else:
		if CTF:
			pass#IMPORTIMPORTIMPORT from morphology   import ctf_img
			pass#IMPORTIMPORTIMPORT from fundamentals import fft
			ctf2 = EMAN2_cppwrap.EMData(nx, nx, 1, False)
			cdata = []
			for im in data:
				ctt = sparx_morphology.ctf_img(nx, im.get_attr("ctf"))
				EMAN2_cppwrap.Util.add_img2(ctf2, ctt)
				cdata.append(sparx_fundamentals.fft(EMAN2_cppwrap.Util.muln_img(sparx_fundamentals.fft(im), ctt)))
		else:
			cdata = [None]*len(data)
			for i in range(len(data)):
				cdata[i] = data[i]

		for N_step in range(len(xrng)):
			for Iter in range(max_iter):
				total_iter += 1
				tavg = sparx_statistics.ave_series(cdata)
				if CTF: tavg = sparx_fundamentals.fft(EMAN2_cppwrap.Util.divn_img(sparx_fundamentals.fft(tavg), ctf2))
				if( FH > 0.0):
					fl = 0.1+(FH-0.1)*Iter/float(max_iter-1)
					tavg = sparx_filter.filt_tanl(tavg, fl, FF)
				"""Multiline Comment52"""
				if total_iter == len(xrng)*max_iter:
					if CTF:
						for i in range(len(data)): data[i].set_attr("xform.align2d",cdata[i].get_attr("xform.align2d"))
					return tavg
				cs = [0,0]
				#print  "  iteration  std   %03d   %7.2f    %7.2f  "%(total_iter,cs[0],cs[1])
				if Iter%4 != 0 or total_iter > max_iter*len(xrng)-10: delta = 0.0
				else:                                                 delta = dst
				sx_sum, sy_sum, nope = sparx_alignment.ali2d_single_iter(cdata, numr, wr, cs, tavg, cnx, cny, \
															xrng[N_step], yrng[N_step], step[N_step], \
															mode=mode, CTF=False, delta=delta)
				#for im in data:  print get_params2D(im)
				#print  "tata ",total_iter-1,Util.infomask(tavg,None,True)
				#tavg.write_image('tata.hdf',total_iter-1)

	return tavg

"""Multiline Comment53"""
def ali3d_mref_Kmeans_MPI(ref_list, outdir, this_data_list_file, Tracker): 
	pass#IMPORTIMPORTIMPORT from utilities      import model_circle, reduce_EMData_to_root, bcast_EMData_to_all, bcast_number_to_all, drop_image
	pass#IMPORTIMPORTIMPORT from utilities      import bcast_list_to_all, get_image, get_input_from_string, get_im, read_text_file
	pass#IMPORTIMPORTIMPORT from utilities      import get_arb_params, set_arb_params, drop_spider_doc, send_attr_dict
	pass#IMPORTIMPORTIMPORT from utilities      import get_params_proj, set_params_proj, model_blank, write_text_file, get_shrink_data_huang
	pass#IMPORTIMPORTIMPORT from filter         import filt_params, filt_btwl, filt_ctf, filt_table, fit_tanh, filt_tanl
	pass#IMPORTIMPORTIMPORT from utilities      import rotate_3D_shift,estimate_3D_center_MPI,get_im
	###-------
	pass#IMPORTIMPORTIMPORT from utilities      import get_attr_stack, get_sorting_attr_stack, get_sorting_params, get_sorting_params_refine
	pass#IMPORTIMPORTIMPORT from utilities      import parsing_sorting_params, fill_in_mpi_list, wrap_mpi_bcast, get_groups_from_partition
	pass#IMPORTIMPORTIMPORT from utilities      import remove_small_groups, set_filter_parameters_from_adjusted_fsc

	###------- 
	pass#IMPORTIMPORTIMPORT from alignment      import Numrinit, prepare_refrings, proj_ali_incore
	pass#IMPORTIMPORTIMPORT from random         import randint
	pass#IMPORTIMPORTIMPORT from filter         import filt_ctf
	pass#IMPORTIMPORTIMPORT from utilities      import print_begin_msg, print_end_msg, print_msg
	pass#IMPORTIMPORTIMPORT from projection     import prep_vol, prgs, prgl, project, prgq, gen_rings_ctf
	pass#IMPORTIMPORTIMPORT from applications   import MPI_start_end
	pass#IMPORTIMPORTIMPORT from reconstruction import rec3D_MPI_noCTF,rec3D_two_chunks_MPI
	pass#IMPORTIMPORTIMPORT from morphology     import binarize, get_shrink_3dmask
	pass#IMPORTIMPORTIMPORT from fundamentals   import fftip, rops_table, fft
	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT from mpi            import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	pass#IMPORTIMPORTIMPORT from mpi            import mpi_reduce, MPI_INT, MPI_SUM
	pass#IMPORTIMPORTIMPORT from time import sleep
	### - reconstruction parameters - No need to change
	fourvar   = False
	debug     = False
	snr       = 1.0
	ref_a     = "S"
	npad      = 2
	############################################################
	pass#IMPORTIMPORTIMPORT from logger import Logger,BaseLogger_Files
	log       = sparx_logger.Logger()
	log   =sparx_logger.Logger(sparx_logger.BaseLogger_Files())
	log.prefix=outdir+"/"
	myid      = Tracker["constants"]["myid"]
	main_node = Tracker["constants"]["main_node"]
	number_of_proc     = Tracker["constants"]["nproc"]
	shrinkage = Tracker["shrinkage"]
	### input parameters
	maxit    = Tracker["constants"]["maxit"]
	ou       = Tracker["constants"]["radius"]
	ir       = Tracker["constants"]["ir"]
	rs       = Tracker["constants"]["rs"]
	xr       = Tracker["constants"]["xr"]
	yr       = Tracker["constants"]["yr"]
	ts       = Tracker["constants"]["ts"]
	delta    = Tracker["constants"]["delta"]
	an       = Tracker["constants"]["an"]
	center   = Tracker["constants"]["center"]
	nassign            = Tracker["constants"]["nassign"]
	nrefine            = Tracker["constants"]["nrefine"]
	CTF                = Tracker["constants"]["CTF"]
	sym                = Tracker["constants"]["sym"]
	termprec           = Tracker["constants"]["stoprnct"]
	if Tracker["constants"]["mask3D"]:
		maskfile = Tracker["mask3D"]
	else:
		maskfile = None
	user_func_name     = Tracker["constants"]["user_func"]
	Tracker["lowpass"] = Tracker["low_pass_filter"]
	Tracker["falloff"] = .1
	if Tracker["constants"]["PWadjustment"]:
		Tracker["PWadjustment"]=Tracker["PW_dict"][Tracker["constants"]["nxinit"]]
	else:
		Tracker["PWadjustment"]=Tracker["constants"]["PWadjustment"]	
	mpi_comm = mpi.MPI_COMM_WORLD
	###--------------------------
	if os.path.exists(outdir): sparx_global_def.ERROR('Output directory exists, please change the name and restart the program', "Kmref_ali3d_MPI ", 1, myid)
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	###
	if myid == main_node:	
		os.mkdir(outdir)
		pass#IMPORTIMPORTIMPORT import global_def
		sparx_global_def.LOGFILE =  os.path.join(outdir, sparx_global_def.LOGFILE)
		log.add("Kmref_ali3d_MPI - Traditional Kmeans clustering  !")
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	######
	#Tracker["applyctf"] = False
	while not os.path.exists(this_data_list_file):
		#print  " my_id",myid
		time.sleep(2)
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)	
	data, old_shifts    = sparx_utilities.get_shrink_data_huang(Tracker,Tracker["nxinit"],this_data_list_file,Tracker["constants"]["partstack"],myid, main_node, number_of_proc, preshift = True)
	if myid ==main_node:
		list_of_particles     = sparx_utilities.read_text_file(this_data_list_file)
		total_nima            = len(list_of_particles)
	else:   
		total_nima            = 0
		list_of_particles     = 0
	total_nima = sparx_utilities.bcast_number_to_all(total_nima,main_node)
	list_of_particles     = sparx_utilities.wrap_mpi_bcast(list_of_particles, main_node)
	
	
	pass#IMPORTIMPORTIMPORT from time import time

	if debug:
		pass#IMPORTIMPORTIMPORT from time import sleep
		while not os.path.exists(outdir):
			print("Node ",myid,"  waiting...")
			time.sleep(5)
		finfo = open(os.path.join(outdir, "progress%04d"%myid), 'w')
		frec  = open( os.path.join(outdir, "recons%04d"%myid), "w" )
	else:
		finfo = None
		frec  = None

	xrng        = sparx_utilities.get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = sparx_utilities.get_input_from_string(yr)
	step        = sparx_utilities.get_input_from_string(ts)
	delta       = sparx_utilities.get_input_from_string(delta)
	lstp = min( len(xrng), len(yrng), len(step), len(delta) )
	if an == "-1":
		an = []
		for i in range(len(xrng)):   an.append(-1)
	else:
		pass#IMPORTIMPORTIMPORT from  alignment	    import proj_ali_incore_local
		an      = sparx_utilities.get_input_from_string(an)
	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(int(ou)*shrinkage+.5)
	center      = int(center)
	image_start, image_end = MPI_start_end(len(Tracker["this_data_list"]), number_of_proc, myid)
	numref      = len(ref_list)
	nx          = ref_list[0].get_xsize()
	if last_ring < 0:       last_ring = nx//2 - 2

	if Tracker["constants"]["focus3Dmask"]:
		focus               = Tracker["focus3D"]
	else:
		focus =      None


	fscmask     = sparx_utilities.model_circle(last_ring, nx, nx, nx)
	stack       = Tracker["constants"]["stack"]
	pass#IMPORTIMPORTIMPORT import user_functions
	user_func = sparx_user_functions.factory[user_func_name]
	if myid == main_node:
		#import user_functions
		#user_func = sparx_user_functions.factory[user_func_name]
		log.add("Input stack                 : %s"%(stack))
		#log.add("Reference volumes           : %s"%(ref_vol))	
		log.add("Number of reference volumes : %i"%(numref))
		log.add("Output directory            : %s"%(outdir))
		log.add("User function               : %s"%(user_func_name))
		log.add("Maskfile                    : %s"%(maskfile))
		if(focus != None):  \
		log.add("Maskfile 3D for focused clustering        : %s"%(focus))
		log.add("Inner radius                : %i"%(first_ring))
		log.add("Outer radius                : %i"%(last_ring))
		log.add("Ring step                   : %i"%(rstep))
		log.add("X search range              : %s"%(xrng))
		log.add("Y search range              : %s"%(yrng))
		log.add("Translational step          : %s"%(step))
		log.add("Angular step                : %s"%(delta))
		log.add("Angular search range        : %s"%(an))
		log.add("Number of assignments in each iteration   : %i"%(nassign))
		log.add("Number of alignments in each iteration    : %i"%(nrefine))
		log.add("Number of iterations                      : %i"%(lstp*maxit) )
		log.add("Center type                 : %i"%(center))
		log.add("CTF correction              : %s"%(CTF))
		log.add("Reference projection method : %s"%(ref_a))
		log.add("Symmetry group              : %s"%(sym))
		log.add("Percentage of change for termination: %f"%(termprec))
		log.add("User function               : %s"%(user_func_name))
		log.add("total number of particles                 : %d"%len(Tracker["this_data_list"]))
		log.add("shrinkage is                              : %f"%shrinkage)
		log.add("the text file for get_shrink_data is %s"%this_data_list_file)
	if maskfile:
		if type(maskfile) is bytes:  mask3D = sparx_morphology.get_shrink_3dmask(Tracker["nxinit"],maskfile)
		else: 	                                mask3D = maskfile
	else:  mask3D = sparx_utilities.model_circle(last_ring, nx, nx, nx)
	numr       = sparx_alignment.Numrinit(first_ring, last_ring, rstep, "F")
	mask2D     = sparx_utilities.model_circle(last_ring, nx, nx) - sparx_utilities.model_circle(first_ring, nx, nx)
	total_nima = len(Tracker["this_data_list"])
	nima       = len(data)
	#####
	if debug:
		finfo.write( "Image_start, image_end: %d %d\n" %(image_start, image_end) )
		finfo.flush()
	start_time = time.time()
	if myid == main_node:
		log.add( "Time to read data: %d\n" % (time.time()-start_time) );start_time = time.time()
	#  Initialize Particle ID and set group number to non-existant -1
	assignment = [-1]*len(data)
	for im in range(len(data)):
		data[im].set_attr_dict({'group':-1})
	"""Multiline Comment71"""
	# volume already filtered before, here is for using userfunc
	#lowpass = 0.5 
	#for  iref in xrange(numref):
	#	set_filter_parameters_from_adjusted_fsc(Tracker["constants"]["total_stack"],Tracker["number_of_ref_class"][iref],Tracker)
	#	lowpass = min(lowpass, Tracker["lowpass"])
	#Tracker["lowpass"] = lowpass
	
	res = 0.5
	for i in range(len(Tracker["global_fsc"][0])-1,0,-1):
		if Tracker["global_fsc"][1][i]>0.5:
			res = fsc_in[0][i]
			break

	Tracker["lowpass"] = min(0.45, res-0.05)
	Tracker["falloff"] = 0.1
	highres = []
	for  iref in range(numref): highres.append(int(res*Tracker["nxinit"]+.5))	
	if myid ==main_node:
	
		for  iref in range(numref):
			ref_list[iref].write_image(os.path.join(outdir, "vol0000.hdf"), iref)
			#set_filter_parameters_from_adjusted_fsc(Tracker["constants"]["total_stack"],Tracker["number_of_ref_class"][iref],Tracker)
			log.add("%d reference low pass filter is %f  %f  %d"%(iref, Tracker["lowpass"],Tracker["falloff"],Tracker["number_of_ref_class"][iref]))
			log.add("%d highres                   %f"%(iref, highres[iref]))
			
			if Tracker["mask3D"]:
				mask3D          = sparx_utilities.get_im(Tracker["mask3D"])
				stat            = EMAN2_cppwrap.Util.infomask(ref_list[iref], mask3D, False)
				ref_list[iref] -= stat[0]
				if stat[1]!=0.0: EMAN2_cppwrap.Util.mul_scalar(ref_list[iref], 1.0/stat[1])
				else:
					pass#IMPORTIMPORTIMPORT from morphology import erosion
					bv = sparx_utilities.model_blank(3, 3, 3)
					bv +=1.
					while stat[1]==0:
						ermask = sparx_morphology.erosion(mask3D, bv)
						stat   = EMAN2_cppwrap.Util.infomask(ref_list[iref], ermask, False)
					EMAN2_cppwrap.Util.mul_scalar(ref_list[iref], 1.0/stat[1])
				
			if(Tracker["constants"]["PWadjustment"]):
				rt = sparx_utilities.read_text_file(Tracker["PW_dict"][Tracker["constants"]["nxinit"]])
				ro = sparx_fundamentals.rops_table(ref_list[iref])
				for i in range(1,len(ro)):  ro[i] = (rt[i]/ro[i])**Tracker["constants"]["upscale"]
				ref_list[iref] = sparx_filter.filt_table(ref_list[iref],ro)
				
			if (Tracker["constants"]["low_pass_filter"]==-1.):  ref_list[iref] = sparx_filter.filt_tanl(ref_list[iref], Tracker["lowpass"], Tracker["falloff"])                                       # low pass from resolution 
			else:                                               ref_list[iref] = sparx_filter.filt_tanl(ref_list[iref], min(Tracker["constants"]["low_pass_filter"]/Tracker["shrinkage"],0.45), Tracker["falloff"]) # user define filter
				
			if Tracker["mask3D"]: EMAN2_cppwrap.Util.mul_img(ref_list[iref], mask3D)
			ref_list[iref].write_image(os.path.join(outdir, "volf0000.hdf"), iref)
			
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

	if CTF:
		#if(data[0].get_attr("ctf_applied") > 0.0):  ERROR("mref_ali3d_MPI does not work for CTF-applied data", "mref_ali3d_MPI", 1, myid)
		pass#IMPORTIMPORTIMPORT from reconstruction import rec3D_MPI
	else:
		pass#IMPORTIMPORTIMPORT from reconstruction import rec3D_MPI_noCTF

	if debug:
		finfo.write( '%d loaded  \n' % len(data) )
		finfo.flush()
	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in range(number_of_proc):
		if im == main_node:  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )
	total_iter = 0
	tr_dummy = EMAN2_cppwrap.Transform({"type":"spider"})

	if(focus != None):
		if(myid == main_node):
			focus = sparx_morphology.get_shrink_3dmask(Tracker["nxinit"],focus)
		else:
			focus =  sparx_utilities.model_blank(nx, nx, nx)
		sparx_utilities.bcast_EMData_to_all(focus, myid, main_node)
		st = EMAN2_cppwrap.Util.infomask(focus, None, True)
		if( st[0] == 0.0 ):  sparx_global_def.ERROR("sxrsort3d","incorrect focused mask, after binarize all values zero",1, myid)
		focus = sparx_projection.prep_vol(focus, 1, 1)


	Niter = int(lstp*maxit*(nassign + nrefine) )
	for Iter in range(Niter):
		N_step = (Iter%(lstp*(nassign+nrefine)))/(nassign+nrefine)
		if Iter%(nassign+nrefine) < nassign:
			runtype = "ASSIGNMENT"
		else:
			runtype = "REFINEMENT"

		total_iter += 1
		if myid == main_node:
			log.add("\n%s ITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f" \
                        %(runtype, total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step]))
			start_ime = time.time()	
		peaks = [ -1.0e23]*nima
		if runtype=="REFINEMENT":
			trans = [tr_dummy]*nima
			pixer = [0.0]*nima
			if(an[N_step] > 0):
				pass#IMPORTIMPORTIMPORT from utilities    import even_angles
				ref_angles = sparx_utilities.even_angles(delta[N_step], symmetry=sym, method = ref_a, phiEqpsi = "Zero")
				# generate list of angles
				pass#IMPORTIMPORTIMPORT from alignment import generate_list_of_reference_angles_for_search
				list_of_reference_angles = \
				sparx_alignment.generate_list_of_reference_angles_for_search(ref_angles, sym=sym)
				del ref_angles
			else:  list_of_reference_angles = [[1.0,1.0]]
		cs = [0.0]*3
		pass#IMPORTIMPORTIMPORT from fundamentals import fft
		if( not focus ):
			for im in range(nima):  data[im] = sparx_fundamentals.fft(data[im])

		for iref in range(numref):
			if myid==main_node: volft = sparx_utilities.get_im(os.path.join(outdir, "volf%04d.hdf"%(total_iter-1)), iref)
			else:				volft=sparx_utilities.model_blank(nx,nx,nx)
			sparx_utilities.bcast_EMData_to_all(volft, myid, main_node)
			volft= sparx_projection.prep_vol(volft,1,1)
			if CTF:
				previous_defocus = -1.0
				if runtype=="REFINEMENT":
					start_time = time.time()
					prjref = sparx_projection.prgq( volft, kb, nx, delta[N_step], ref_a, sym, MPI=True)
					if myid == main_node:
						log.add( "Calculation of projections: %d" % (time.time()-start_time) );start_time = time.time()
					del volft, kb
			else:
				if runtype=="REFINEMENT":
					start_time = time.time()
					refrings = sparx_alignment.prepare_refrings( volft, kb, nx, delta[N_step], ref_a, sym, numr)
					if myid == main_node:
						log.add( "Initial time to prepare rings: %d" % (time.time()-start_time) );start_time = time.time()
					del volft, kb

			start_time = time.time()
			for im in range(nima):
				if CTF:
					ctf = data[im].get_attr( "ctf" )
					if runtype=="REFINEMENT":
						if ctf.defocus != previous_defocus:
							previous_defocus = ctf.defocus
							rstart_time = time.time()
							refrings = sparx_projection.gen_rings_ctf( prjref, nx, ctf, numr)
							if myid == main_node: log.add( "Repeated time to prepare rings: %d" % (time.time()-rstart_time) );rstart_time = time.time()
				if runtype=="ASSIGNMENT":
					phi,tht,psi,s2x,s2y = sparx_utilities.get_params_proj(data[im])
					"""Multiline Comment72"""
					#  Standard distance
					#  Ref is in reciprocal space
					if CTF:
						ref = sparx_filter.filt_ctf( sparx_projection.prgl( volft, [phi,tht,psi,-s2x,-s2y], 1, False), ctf )
					else:
						ref = sparx_projection.prgl( volft, [phi,tht,psi,-s2x,-s2y], 1, False)
					pass#IMPORTIMPORTIMPORT from filter import filt_tophatl
					pass#IMPORTIMPORTIMPORT from math import sqrt
					ref = sparx_filter.filt_tophatl(ref, float(highres[iref])/(ref.get_ysize()))
					ref.set_attr("is_complex",0)
					ref.set_value_at(0,0,0.0)
					nrmref = numpy.sqrt(EMAN2_cppwrap.Util.innerproduct(ref, ref, None))
					if(focus):
						mask2D = sparx_morphology.binarize( sparx_projection.prgl( focus, [phi,tht,psi,-s2x,-s2y]), 1)
						tempx = sparx_fundamentals.fft(data[im]*mask2D)
						tempx.set_attr("is_complex",0)
						peak = EMAN2_cppwrap.Util.innerproduct(ref, tempx, None)
					else:
						data[im].set_attr("is_complex",0)
						peak = EMAN2_cppwrap.Util.innerproduct(ref, data[im], None)
						data[im].set_attr("is_complex",1)
					peak /= nrmref


					"""Multiline Comment73"""

				if peak > peaks[im]:
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
			if myid == main_node:log.add( "Time to process particles for reference %3d: %d" % (iref, time.time()-start_time) );start_time = time.time()
		del peaks
		if runtype=="ASSIGNMENT":  del volft, ref #kb, ref
		else:
			if CTF: del prjref
			del refrings
			if an[N_step] > 0: del list_of_reference_angles


		#  compute number of particles that changed assignment and how many are in which group
		nchng = 0
		npergroup = [0]*numref
		for im in range(nima):
			iref = data[im].get_attr('group')
			npergroup[iref] += 1
			if iref != assignment[im]:
				assignment[im] = iref
				nchng += 1
		nchng            = mpi.mpi_reduce(nchng, 1, mpi.MPI_INT, mpi.MPI_SUM, 0, mpi.MPI_COMM_WORLD)
		npergroup        = mpi.mpi_reduce(npergroup, numref, mpi.MPI_INT, mpi.MPI_SUM, 0, mpi.MPI_COMM_WORLD)
		npergroup        = list(map(int, npergroup))
		terminate        = 0
		empty_group      = 0
		empty_group_list = []
		if myid == main_node:
			ngroup = []
			nchng = int(nchng[0])
			precn = 100*float(nchng)/float(total_nima)
			msg = " Number of particles that changed assignments %7d, percentage of total: %5.1f"%(nchng, precn)
			log.add(msg)
			msg = " Group       number of particles"
			log.add(msg)
			for iref in range(numref):
				msg = " %5d       %7d"%(iref+1, npergroup[iref])
				log.add(msg)
				if npergroup[iref]==0:
					empty_group =1
					empty_group_list.append(iref)
				ngroup.append(int(npergroup[iref]))
			if precn <= termprec:
				terminate = 1
			if empty_group ==1:
				terminate = 1
		else:
			ngroup  = 0
		terminate   = mpi.mpi_bcast(terminate, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD)
		terminate   = int(terminate[0])
		ngroup      = sparx_utilities.wrap_mpi_bcast(ngroup,main_node)
		empty_group = mpi.mpi_bcast(empty_group, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD)
		empty_group = int(empty_group[0])
		if empty_group ==1: break # program stops whenever empty_group appears!
		if runtype=="REFINEMENT":
			for im in range(nima):  data[im].set_attr('xform.projection', trans[im])
			if center == -1:
				cs[0], cs[1], cs[2], dummy, dummy = sparx_utilities.estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)				
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f"%(cs[0], cs[1], cs[2])
					log.add(msg)
				cs = mpi.mpi_bcast(cs, 3, mpi.MPI_FLOAT, main_node, mpi.MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				sparx_utilities.rotate_3D_shift(data, cs)
			#output pixel errors
			pass#IMPORTIMPORTIMPORT from mpi import mpi_gatherv
			recvbuf = mpi.mpi_gatherv(pixer, nima, mpi.MPI_FLOAT, recvcount, disps, mpi.MPI_FLOAT, main_node, mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			if myid == main_node:
				recvbuf = list(map(float, recvbuf))
				pass#IMPORTIMPORTIMPORT from statistics import hist_list
				lhist = 20
				region, histo = sparx_statistics.hist_list(recvbuf, lhist)
				if region[0] < 0.0:  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles"
				log.add(msg)
				for lhx in range(lhist):
					msg = " %10.3f      %7d"%(region[lhx], histo[lhx])
					log.add(msg)
				del region, histo
			del recvbuf
		#if CTF: del vol
		fscc = [None]*numref
		if fourvar and runtype=="REFINEMENT":
			sumvol = sparx_utilities.model_blank(nx, nx, nx)
		sart_time = time.time()
	
		#lowpass = 0.5
		#for  iref in xrange(numref):
		#	set_filter_parameters_from_adjusted_fsc(Tracker["constants"]["total_stack"],Tracker["number_of_ref_class"][iref],Tracker)
		#lowpass = min(lowpass, Tracker["lowpass"])
		#Tracker["lowpass"] = lowpass

		if( not focus ):
			for im in range(nima):  data[im] = sparx_fundamentals.fft(data[im])

		highres = []
		lowpass_tmp =[]
		tmpref =[]
		pass#IMPORTIMPORTIMPORT from statistics import fsc
		for iref in range(numref):
			#  3D stuff
			pass#IMPORTIMPORTIMPORT from time import localtime, strftime
			if Tracker["constants"]["3d-interpolation"]==" ":
				chunk_id   = 0
				niter      = 10
				upweighted = False
				compensate = False
				verbose    =  0
				sign       = 1
				volref0    = recons3d_n_trl_MPI_one_node(data, CTF, snr, sign, npad, sym, iref, niter, verbose, upweighted, compensate, chunk_id)
				chunk_id   = 1
				volref1    = recons3d_n_trl_MPI_one_node(data, CTF, snr, sign, npad, sym, iref, niter, verbose, upweighted, compensate, chunk_id)
				if myid == main_node:
					fscc[iref] = sparx_statistics.fsc(volref0, volref1)
					volref  = volref0 + volref1
					del volref1
					del volref0
				#chunk_id   = -1
				#volref = recons3d_n_trl_MPI_one_node(data, CTF, snr, sign, npad, sym, iref, niter, verbose, upweighted, compensate, chunk_id)
			else:
				if CTF: volref, fscc[iref] = sparx_reconstruction.rec3D_two_chunks_MPI(data, snr, sym, fscmask, os.path.join(outdir, "resolution_%02d_%04d"%(iref, total_iter)), myid, main_node, index = iref, npad = npad, finfo=frec)
				else:   volref, fscc[iref] = sparx_reconstruction.rec3D_MPI_noCTF(data, sym, fscmask, os.path.join(outdir, "resolution_%02d_%04d"%(iref, total_iter)), myid, main_node, index = iref, npad = npad, finfo=frec)
			if myid == main_node:
				log.add( "Time to compute 3D: %d" % (time.time()-start_time) );start_time = time.time()
				volref.write_image(os.path.join(outdir, "vol%04d.hdf"%( total_iter)), iref)
				if fourvar and runtype=="REFINEMENT": sumvol += volref
				res = 0.5
				for ifreq in range(len(fscc[iref][0])-1,0,-1):
					if fscc[iref][1][ifreq] > 0.5 : # always use .5 as cutoff
						res = fscc[iref][0][ifreq]
						break
				Tracker["lowpass"] = min(0.45, res)
				#Tracker["lowpass"] = max(0.11, res)
				Tracker["falloff"] = 0.1
				log.add(" low pass filter is %f    %f   %d"%(Tracker["lowpass"], Tracker["falloff"], ngroup[iref]))
				tmpref.append(volref)
			else:
				Tracker["lowpass"] = 0.0
				Tracker["falloff"] = 0.0
				res   = 0.0
			Tracker["lowpass"] = sparx_utilities.wrap_mpi_bcast(Tracker["lowpass"], main_node, mpi_comm)
			Tracker["falloff"] = sparx_utilities.wrap_mpi_bcast(Tracker["falloff"], main_node, mpi_comm)
			res = sparx_utilities.wrap_mpi_bcast(res, main_node, mpi_comm)
			highres.append(int(res*Tracker["nxinit"]+ 0.5))
			lowpass_tmp.append(Tracker["lowpass"])
		Tracker["lowpass"]=min(lowpass_tmp)
		if myid ==main_node:
				log.add(" the adopted low pass filter is %f    %f   "%(Tracker["lowpass"], Tracker["falloff"]))
		for iref in range(numref):
			if myid == main_node:
				log.add("%d highres                   %f"%(iref, highres[iref]))
				if Tracker["mask3D"]:
					mask3D = sparx_utilities.get_im(Tracker["mask3D"])
					stat = EMAN2_cppwrap.Util.infomask(tmpref[iref], mask3D, False)
					tmpref[iref] -= stat[0]
					EMAN2_cppwrap.Util.mul_scalar(tmpref[iref], 1.0/stat[1])
					
				if(Tracker["constants"]["PWadjustment"]):
				
					rt = sparx_utilities.read_text_file(Tracker["PW_dict"][Tracker["constants"]["nxinit"]])
					ro = sparx_fundamentals.rops_table(tmpref[iref])
					for i in range(1,len(ro)):  ro[i] = (rt[i]/ro[i])**Tracker["constants"]["upscale"]
					tmpref[iref] =sparx_filter.filt_table(tmpref[iref],ro)

				if (Tracker["constants"]["low_pass_filter"]==-1.):  tmpref[iref] = sparx_filter.filt_tanl(tmpref[iref], Tracker["lowpass"], Tracker["falloff"])                                       # low pass from resolution 
				else:                                               tmpref[iref] = sparx_filter.filt_tanl(tmpref[iref], min(Tracker["constants"]["low_pass_filter"]/Tracker["shrinkage"],0.45), Tracker["falloff"]) # user define filter			
					
				if Tracker["mask3D"]: EMAN2_cppwrap.Util.mul_img(tmpref[iref], mask3D)
				tmpref[iref].write_image(os.path.join(outdir, "volf%04d.hdf"%( total_iter)), iref)
		del tmpref

		if runtype=="REFINEMENT":
			if fourvar:
				varf = sparx_statistics.varf3d_MPI(data, os.path.join(outdir, "ssnr%04d"%total_iter), None,sumvol,last_ring, 1.0, 1, CTF, 1, sym, myid)
				if myid == main_node:   
					varf = 1.0/varf
					varf.write_image( os.path.join(outdir,"varf%04d.hdf"%total_iter) )                            		
		"""Multiline Comment74"""
		#  here we  write header info
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		#start_time = time()
		if runtype=="REFINEMENT": par_str = ['xform.projection', 'ID', 'group']
		else: par_str = ['group', 'ID' ]
		#if myid == main_node:
		#	from utilities import file_type
	        #	if file_type(stack) == "bdb":
	        #		from utilities import recv_attr_dict_bdb
	        #		recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        #	else:
	        # 		from utilities import recv_attr_dict
	        #		recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        #else:		send_attr_dict(main_node, data, par_str, image_start, image_end)
		if terminate == 1:
			if myid == main_node: log.add("Kmref_ali3d_MPI terminated due to small number of objects changing assignments")
			#final_list = get_sorting_params(Tracker,data)
			#res_groups = get_groups_from_partition(final_list, Tracker["this_data_list"], numref)
			#if myid ==main_node:
			#	nc = 0
			#	final_list_saved_file =os.path.join(outdir,"list2.txt")
			#	write_text_file(final_list,final_list_saved_file)
			#	for igrp in xrange(len(res_groups)):
			#			if len(res_groups[igrp])>0:
			#				saved_file = os.path.join(outdir,"Class%d.txt"%nc)
			#				write_text_file(res_groups[igrp],saved_file)
			#				nc +=1
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			#Tracker["this_partition"]=final_list
			break
		if myid == main_node:log.add( "Time to write headers: %d\n" % (time.time()-start_time) )
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	######writing partition only in the end of the program
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	if nrefine!=0: par_str = ['xform.projection', 'ID', 'group']
	else: par_str = ['group', 'ID' ]
	"""Multiline Comment75"""
	if myid == main_node:log.add("Kmref_ali3d_MPI is done!")
	final_list = sparx_utilities.get_sorting_params_refine(Tracker,data, total_nima)
	group_list, ali3d_params_list = sparx_utilities.parsing_sorting_params(final_list)
	res_groups = sparx_utilities.get_groups_from_partition(group_list, list_of_particles, numref)
	final_group_list, res_groups = sparx_utilities.remove_small_groups(res_groups, Tracker["constants"]["smallest_group"])
	if myid ==main_node:
		nc = 0
		sparx_utilities.write_text_file(group_list, os.path.join(outdir,"list2.txt"))
		for igrp in range(len(res_groups)):
			if len(res_groups[igrp])>0:
				sparx_utilities.write_text_file(res_groups[igrp], os.path.join(outdir,"Class%d.txt"%nc))
				nc +=1
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	Tracker["this_partition"]=final_list
	return empty_group,res_groups,final_group_list
##########################


def mref_ali3d_EQ_Kmeans(ref_list, outdir, particle_list_file, Tracker):
	pass#IMPORTIMPORTIMPORT from utilities      import model_circle, reduce_EMData_to_root, bcast_EMData_to_all, bcast_number_to_all, drop_image
	pass#IMPORTIMPORTIMPORT from utilities      import  bcast_list_to_all, get_image, get_input_from_string, get_im
	pass#IMPORTIMPORTIMPORT from utilities      import get_arb_params, set_arb_params, drop_spider_doc, send_attr_dict
	pass#IMPORTIMPORTIMPORT from utilities      import get_params_proj, set_params_proj, model_blank, wrap_mpi_bcast, write_text_file
	pass#IMPORTIMPORTIMPORT from filter         import filt_params, filt_btwl, filt_ctf, filt_table, fit_tanh, filt_tanl
	pass#IMPORTIMPORTIMPORT from utilities      import rotate_3D_shift,estimate_3D_center_MPI, get_shrink_data_huang, get_im
	####-------
	pass#IMPORTIMPORTIMPORT from utilities      import get_attr_stack, get_sorting_attr_stack, get_sorting_params, get_sorting_params_refine
	pass#IMPORTIMPORTIMPORT from utilities      import parsing_sorting_params, fill_in_mpi_list, wrap_mpi_bcast, get_groups_from_partition
	pass#IMPORTIMPORTIMPORT from utilities      import remove_small_groups, set_filter_parameters_from_adjusted_fsc
	###-----------
	pass#IMPORTIMPORTIMPORT from alignment      import Numrinit, prepare_refrings, proj_ali_incore
	pass#IMPORTIMPORTIMPORT from random         import randint, random
	pass#IMPORTIMPORTIMPORT from filter         import filt_ctf
	pass#IMPORTIMPORTIMPORT from utilities      import print_begin_msg, print_end_msg, print_msg, read_text_file
	pass#IMPORTIMPORTIMPORT from projection     import prep_vol, prgs, prgl, project, prgq, gen_rings_ctf
	pass#IMPORTIMPORTIMPORT from morphology     import binarize, get_shrink_3dmask
	pass#IMPORTIMPORTIMPORT from statistics		import fsc

	pass#IMPORTIMPORTIMPORT import os
	pass#IMPORTIMPORTIMPORT import types
	pass#IMPORTIMPORTIMPORT from mpi            import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	pass#IMPORTIMPORTIMPORT from mpi            import mpi_reduce, mpi_gatherv, mpi_scatterv, MPI_INT, MPI_SUM
	pass#IMPORTIMPORTIMPORT from applications   import MPI_start_end
	pass#IMPORTIMPORTIMPORT from reconstruction import rec3D_two_chunks_MPI, rec3D_MPI_noCTF
	pass#IMPORTIMPORTIMPORT from fundamentals   import fftip, rops_table, fft
	mpi_comm = mpi.MPI_COMM_WORLD
	#####  reconstruction parameters, no need to change.
	fourvar   = False
	snr       = 1.0
	debug     = False
	ref_a     = "S"
	npad      = 2
	#####################################################
	pass#IMPORTIMPORTIMPORT from logger import Logger,BaseLogger_Files
	log       = sparx_logger.Logger()
	log   =sparx_logger.Logger(sparx_logger.BaseLogger_Files())
	log.prefix=outdir+"/"
	myid           = Tracker["constants"]["myid"]
	main_node      = Tracker["constants"]["main_node"]
	number_of_proc = Tracker["constants"]["nproc"]
	shrinkage      = Tracker["shrinkage"]
	### input parameters
	maxit          = Tracker["constants"]["maxit"]
	ou             = Tracker["constants"]["radius"]
	ir             = Tracker["constants"]["ir"]
	rs             = Tracker["constants"]["rs"]
	xr             = Tracker["constants"]["xr"]
	yr             = Tracker["constants"]["yr"]
	ts             = Tracker["constants"]["ts"]
	delta          = Tracker["constants"]["delta"]
	an             = Tracker["constants"]["an"]
	center         = Tracker["constants"]["center"]
	nassign        = Tracker["constants"]["nassign"]
	nrefine        = Tracker["constants"]["nrefine"]
	CTF            = Tracker["constants"]["CTF"]
	sym            = Tracker["constants"]["sym"]
	termprec       = Tracker["constants"]["stoprnct"]
	if Tracker["constants"]["mask3D"]:
		maskfile            = Tracker["mask3D"]
	else:
		maskfile            = None
	if Tracker["constants"]["focus3Dmask"]: focus = Tracker["constants"]["focus3Dmask"]
	else:                                   focus = None
	partstack           = Tracker["constants"]["partstack"]
	user_func_name      = Tracker["constants"]["user_func"]
	Tracker["lowpass"]  = Tracker["low_pass_filter"]
	Tracker["falloff"]  = .1
	if Tracker["constants"]["PWadjustment"] !='':
		Tracker["PWadjustment"] = Tracker["PW_dict"][Tracker["constants"]["nxinit"]]
	else:
		Tracker["PWadjustment"] = None	
	####################################################
	pass#IMPORTIMPORTIMPORT from time import sleep
	#Tracker["applyctf"] = True # 
	data, old_shifts = sparx_utilities.get_shrink_data_huang(Tracker, Tracker["nxinit"], particle_list_file, partstack, myid,main_node,number_of_proc,preshift=True)
	if myid ==main_node:
		list_of_particles = sparx_utilities.read_text_file(particle_list_file)
		total_nima        =len(list_of_particles)
	else:
		total_nima = 0
		list_of_particles = 0
	total_nima            = sparx_utilities.bcast_number_to_all(total_nima,main_node)
	list_of_particles     = sparx_utilities.wrap_mpi_bcast(list_of_particles, main_node)
	
	if myid == main_node:
		if os.path.exists(outdir):  nx = 1
		else:  nx = 0
	else:  nx = 0
	ny = sparx_utilities.bcast_number_to_all(nx, source_node = main_node)
	#image_start, image_end = MPI_start_end(total_data,number_of_proc, myid)
	
	if ny == 1:  sparx_global_def.ERROR('Output directory exists, please change the name and restart the program', "mref_ali3d_iter", 1,myid)
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	if myid == main_node:	
		os.mkdir(outdir)
		pass#IMPORTIMPORTIMPORT import global_def
		sparx_global_def.LOGFILE =  os.path.join(outdir, sparx_global_def.LOGFILE)
		log.add("Equal K-means  ")
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	pass#IMPORTIMPORTIMPORT from time import time
	if debug:
		pass#IMPORTIMPORTIMPORT from time import sleep
		while not os.path.exists(outdir):
			print("Node ",myid,"  waiting...")
			time.sleep(5)
		finfo = open(os.path.join(outdir, "progress%04d"%myid), 'w')
		frec  = open( os.path.join(outdir, "recons%04d"%myid), "w" )
	else:
		finfo = None
		frec  = None
	xrng        = sparx_utilities.get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = sparx_utilities.get_input_from_string(yr)
	step        = sparx_utilities.get_input_from_string(ts)
	delta       = sparx_utilities.get_input_from_string(delta)
	lstp = min( len(xrng), len(yrng), len(step), len(delta) )
	if (an == "-1"):
		an = []
		for i in range(len(xrng)):   an.append(-1)
	else:
		pass#IMPORTIMPORTIMPORT from  alignment	    import proj_ali_incore_local
		an      = sparx_utilities.get_input_from_string(an)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(int(ou)*shrinkage+.5)
	center      = int(center)
	numref  = len(ref_list)
	nx      = ref_list[0].get_xsize()
	if last_ring < 0:last_ring = nx//2 - 2
	pass#IMPORTIMPORTIMPORT import user_functions
	user_func = sparx_user_functions.factory[user_func_name]
	if (myid == main_node):
		#import user_functions
		#user_func = sparx_user_functions.factory[user_func_name]
		log.add("mref_ali3d_MPI")
		log.add("Input stack                               : %s"%(Tracker["constants"]["stack"]))
		#log.add("Reference volumes                         : %s"%(ref_vol))	
		log.add("Number of reference volumes               : %i"%(numref))
		log.add("Output directory                          : %s"%(outdir))
		log.add("User function                             : %s"%(user_func_name))
		if(focus != None):  \
		log.add("Maskfile 3D for focused clustering        : %s"%(focus))
		log.add("Overall 3D mask applied in user function  : %s"%(maskfile))
		log.add("Inner radius                              : %i"%(first_ring))
		log.add("Outer radius                              : %i"%(last_ring))
		log.add("Ring step                                 : %i"%(rstep))
		log.add("X search range                            : %s"%(xrng))
		log.add("Y search range                            : %s"%(yrng))
		log.add("Translational step                        : %s"%(step))
		log.add("Angular step                              : %s"%(delta))
		log.add("Angular search range                      : %s"%(an))
		log.add("Number of assignments in each iteration   : %i"%(nassign))
		log.add("Number of alignments in each iteration    : %i"%(nrefine))
		log.add("Number of iterations                      : %i"%(lstp*maxit) )
		log.add("Center type                               : %i"%(center))
		log.add("CTF correction                            : %s"%(CTF))
		log.add("Symmetry group                            : %s"%(sym))
		log.add("Percentage of change for termination      : %f"%(termprec))
		log.add("User function                             : %s"%(user_func_name))
		log.add("total number of particles                 : %d"%len(Tracker["this_data_list"]))
		log.add("shrinkage is                              : %f"%shrinkage)
		log.add("the particle id files for get_shrink_dat is %s"%particle_list_file) 
	if(maskfile):
		if(type(maskfile) is bytes): mask3D = sparx_morphology.get_shrink_3dmask(Tracker["nxinit"],maskfile) 
		else: 	                                mask3D = maskfile
	else        :  mask3D = sparx_utilities.model_circle(last_ring, nx, nx, nx)
	numr     = sparx_alignment.Numrinit(first_ring, last_ring, rstep, "F")
	mask2D   = sparx_utilities.model_circle(last_ring, nx, nx)
	if(first_ring > 1):  mask2D -= sparx_utilities.model_circle(first_ring,nx,nx)
	nima       = len(data) # local number of particles
	
	"""Multiline Comment76"""
	if debug:
		finfo.write( "Image_start, image_end: %d %d\n" %(image_start,image_end))
		finfo.flush()
	start_time = time.time()
	#  Here the assumption is that input are always volumes.  It should be most likely be changed so optionally these are group assignments.
	#  Initialize Particle ID and set group number to non-existant -1
	for im in range(nima):
		data[im].set_attr_dict({'group':-1})
	if(myid == 0):
		log.add( "Time to read data: %d" % (time.time()-start_time) );start_time = time.time()

	if myid == main_node:
		refdata = [None]*7
		for  iref in range(numref):
			ref_list[iref].write_image(os.path.join(outdir, "vol0000.hdf"), iref)
		#refdata[0] = 
		#user_func(refdata)
	refdata        =[None]*4
	res = 0.5
	for i in range(len(Tracker["global_fsc"][0])-1,0,-1):
		if Tracker["global_fsc"][1][i] > 0.5:
				res=fsc_in[0][i]
				break
	highres = []
	for  iref in range(numref): highres.append(int(res*Tracker["nxinit"] + 0.5))
	Tracker["lowpass"] = min(0.45, res-0.05)
	Tracker["lowpass"] = max(0.11,res)
	Tracker["falloff"] = 0.1
	##-----------------------------------------------
	if myid == main_node:  ### 3-D mask, low pass filter, and power spectrum adjustment
		for iref in range(numref):
		
			#set_filter_parameters_from_adjusted_fsc(Tracker["constants"]["total_stack"],Tracker["number_of_ref_class"][iref],Tracker)
			log.add("%d  low pass filter   %f %f  %d"%(iref,Tracker["lowpass"],Tracker["falloff"],Tracker["number_of_ref_class"][iref]))
			log.add("%d highres                   %f"%(iref, highres[iref]))
			
			if Tracker["mask3D"]:
				mask3D = sparx_utilities.get_im(Tracker["mask3D"])
				stat   = EMAN2_cppwrap.Util.infomask(ref_list[iref], mask3D, False)
				ref_list[iref] -= stat[0]
				if stat[1]!=0.0: EMAN2_cppwrap.Util.mul_scalar(ref_list[iref], 1.0/stat[1])
				else:
					pass#IMPORTIMPORTIMPORT from morphology import erosion
					bv = sparx_utilities.model_blank(3, 3, 3)
					bv +=1.
					while stat[1]==0:
						ermask = sparx_morphology.erosion(mask3D, bv)
						stat   = EMAN2_cppwrap.Util.infomask(ref_list[iref], ermask, False)
					EMAN2_cppwrap.Util.mul_scalar(ref_list[iref], 1.0/stat[1])
					
			if(Tracker["constants"]["PWadjustment"] != ''):
				rt = sparx_utilities.read_text_file(Tracker["PW_dict"][Tracker["constants"]["nxinit"]])
				ro = sparx_fundamentals.rops_table(ref_list[iref])
				for i in range(1,len(ro)):  ro[i] = (rt[i]/ro[i])**Tracker["constants"]["upscale"]
				ref_list[iref] =sparx_filter.filt_table(ref_list[iref],ro)
				
			if (Tracker["constants"]["low_pass_filter"]==-1.):  ref_list[iref] = sparx_filter.filt_tanl(ref_list[iref], Tracker["lowpass"], Tracker["falloff"])                                       # low pass from resolution 
			else:                                               ref_list[iref] = sparx_filter.filt_tanl(ref_list[iref], min(Tracker["constants"]["low_pass_filter"]/Tracker["shrinkage"],0.45), Tracker["falloff"]) # user define filter
							

			if Tracker["mask3D"]: EMAN2_cppwrap.Util.mul_img(ref_list[iref], mask3D)
			ref_list[iref].write_image(os.path.join(outdir,"volf0000.hdf"),iref)
			
	mpi.mpi_barrier( mpi.MPI_COMM_WORLD )
	if CTF:
		#if(data[0].get_attr_default("ctf_applied",0) > 0):  ERROR("mref_ali3d_MPI does not work for CTF-applied data", "mref_ali3d_MPI", 1, myid)
		pass#IMPORTIMPORTIMPORT from reconstruction import rec3D_MPI
	else:
		pass#IMPORTIMPORTIMPORT from reconstruction import rec3D_MPI_noCTF
	if debug:
		finfo.write( '%d loaded  \n' % len(data) )
		finfo.flush()
	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in range(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                   disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )
	total_iter = 0
	tr_dummy = EMAN2_cppwrap.Transform({"type":"spider"})

	if(focus != None):
		if(myid == main_node):
			focus = sparx_morphology.get_shrink_3dmask(Tracker["nxinit"], focus)
		else:
			focus =  sparx_utilities.model_blank(Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"])
		sparx_utilities.bcast_EMData_to_all(focus, myid, main_node)
		st = EMAN2_cppwrap.Util.infomask(focus, None, True)
		if( st[0] == 0.0 ):  sparx_global_def.ERROR("sxrsort3d","incorrect focused mask, after binarize all values zero",1, myid)
		focus = sparx_projection.prep_vol(focus, 1, 1)

	Niter = int(lstp*maxit*(nassign + nrefine) )
	for Iter in range(Niter):
		N_step = (Iter%(lstp*(nassign+nrefine)))/(nassign+nrefine)
		if Iter%(nassign+nrefine) < nassign:
			runtype = "ASSIGNMENT"
		else:
			runtype = "REFINEMENT"

		total_iter += 1
		if(myid == main_node):
			log.add("\n%s ITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f  \
			"%(runtype, total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step]))
			start_ime = time.time()
		peaks =  [ [ -1.0e23 for im in range(nima) ] for iref in range(numref) ]
		if runtype=="REFINEMENT":
			trans = [ [ tr_dummy for im in range(nima) ] for iref in range(numref) ]
			pixer = [ [  0.0     for im in range(nima) ] for iref in range(numref) ]
			if(an[N_step] > 0):
				pass#IMPORTIMPORTIMPORT from utilities    import even_angles
				ref_angles = sparx_utilities.even_angles(delta[N_step], symmetry=sym, method = ref_a, phiEqpsi = "Zero")
				# generate list of angles
				pass#IMPORTIMPORTIMPORT from alignment import generate_list_of_reference_angles_for_search
				list_of_reference_angles = \
				sparx_alignment.generate_list_of_reference_angles_for_search(ref_angles, sym=sym)
				del ref_angles
			else:  list_of_reference_angles = [[1.0,1.0]]
		cs = [0.0]*3
		pass#IMPORTIMPORTIMPORT from fundamentals import fft
		if( not focus ):
			for im in range(nima):  data[im] = sparx_fundamentals.fft(data[im])

		for iref in range(numref):
			if(myid == main_node):	volft = sparx_utilities.get_im(os.path.join(outdir, "volf%04d.hdf"%(total_iter-1)), iref)
			else: 					volft =  sparx_utilities.model_blank(nx, nx, nx)
			sparx_utilities.bcast_EMData_to_all(volft, myid, main_node)
			volft = sparx_projection.prep_vol(volft, 1, 1)
			if CTF:
				previous_defocus=-1.0
				if runtype=="REFINEMENT":
					start_time = time.time()
					prjref = sparx_projection.prgq( volft, kb, nx, delta[N_step], ref_a, sym, MPI=True)
					if(myid == 0): log.add( "Calculation of projections: %d" % (time.time()-start_time) );start_time = time.time()
					del volft, kb
			else:
				if runtype=="REFINEMENT":
					start_time = time.time()
					refrings = sparx_alignment.prepare_refrings( volft, kb, nx, delta[N_step], ref_a, sym, numr)
					if(myid == 0): log.add( "Initial time to prepare rings: %d" % (time.time()-start_time) );start_time = time.time()
					del volft, kb
			start_time = time.time()
			for im in range(nima):
				if(CTF):
					ctf = data[im].get_attr("ctf")
					if runtype=="REFINEMENT":
						if(ctf.defocus != previous_defocus):
							previous_defocus = ctf.defocus
							rstart_time = time.time()
							refrings = sparx_projection.gen_rings_ctf( prjref, nx, ctf, numr)
							if(myid == 0): log.add( "Repeated time to prepare rings: %d" % (time.time()-rstart_time) );rstart_time = time.time()

				if runtype=="ASSIGNMENT":
					phi,tht,psi,s2x,s2y = sparx_utilities.get_params_proj(data[im])
					"""Multiline Comment77"""
					#  Standard distance
					#  Ref is in reciprocal space
					if CTF:
						ref = sparx_filter.filt_ctf( sparx_projection.prgl( volft, [phi,tht,psi,-s2x,-s2y], 1, False), ctf )
					else:
						ref = sparx_projection.prgl( volft, [phi,tht,psi,-s2x,-s2y], 1, False)
					pass#IMPORTIMPORTIMPORT from filter import filt_tophatl
					pass#IMPORTIMPORTIMPORT from math import sqrt
					ref = sparx_filter.filt_tophatl(ref, float(highres[iref])/(ref.get_ysize()))
					ref.set_attr("is_complex",0)
					ref.set_value_at(0,0,0.0)
					nrmref = numpy.sqrt(EMAN2_cppwrap.Util.innerproduct(ref, ref, None))
					if(focus):
						mask2D = sparx_morphology.binarize( sparx_projection.prgl( focus, [phi,tht,psi,-s2x,-s2y]), 1)
						tempx = sparx_fundamentals.fft(data[im]*mask2D)
						tempx.set_attr("is_complex",0)
						peak = EMAN2_cppwrap.Util.innerproduct(ref, tempx, None)
					else:
						data[im].set_attr("is_complex",0)
						peak = EMAN2_cppwrap.Util.innerproduct(ref, data[im], None)
						data[im].set_attr("is_complex",1)
					peak /= nrmref

					"""Multiline Comment78"""



					if not(finfo is None):
						finfo.write( "ID, iref, peak: %6d %d %8.5f\n" % (list_of_particles[im],iref,peak) )
				else:
					if(an[N_step] == -1):
						peak, pixel_error = sparx_alignment.proj_ali_incore(data[im], refrings, numr, xrng[N_step], yrng[N_step], step[N_step])
					else:
						peak, pixel_error = sparx_alignment.proj_ali_incore_local(data[im], refrings, list_of_reference_angles, numr, \
																	xrng[N_step], yrng[N_step], step[N_step], an[N_step],sym=sym)
					if not(finfo is None):
						phi,tht,psi,s2x,s2y = sparx_utilities.get_params_proj(data[im])
						finfo.write( "ID, iref, peak,t rans: %6d %d %f %f %f %f %f %f\n"%(list_of_particles[im],iref,peak,phi,tht,psi,s2x,s2y) )
						finfo.flush()

				peaks[iref][im] = peak
				if runtype=="REFINEMENT":
					pixer[iref][im] = pixel_error
					trans[iref][im] = data[im].get_attr( "xform.projection" )
			if(myid == 0):log.add( "Time to process particles for reference %3d: %d" % (iref, time.time()-start_time) );start_time = time.time()
		if runtype=="ASSIGNMENT":  del volft, ref #kb, ref
		else:
			if CTF: del prjref
			del refrings
			if(an[N_step] > 0): del list_of_reference_angles
		#  send peak values to the main node, do the assignments, and bring them back
		pass#IMPORTIMPORTIMPORT from numpy import float32, empty, inner, abs
		if( myid == 0 ):
			dtot = numpy.empty( (numref, total_nima), dtype = numpy.float32)
		for  iref in range(numref):
			recvbuf = mpi.mpi_gatherv(peaks[iref], nima, mpi.MPI_FLOAT, recvcount, disps, mpi.MPI_FLOAT, main_node, mpi.MPI_COMM_WORLD)
			if( myid == 0 ): dtot[iref] = recvbuf
		del recvbuf
		#  The while loop over even angles delta should start here.
		#  prepare reference directions
		pass#IMPORTIMPORTIMPORT from utilities import even_angles, getvec
		if   Tracker["constants"]["protein_shape"]=="g"  :refa = sparx_utilities.even_angles(60.0)     # globular proteins
		elif Tracker["constants"]["protein_shape"]=="f"  :refa = sparx_utilities.even_angles(40.0, theta1=65, theta2=115) # filament proteins
		numrefang = len(refa)
		refanorm = numpy.empty( (numrefang, 3), dtype = numpy.float32)
		for i in range(numrefang):
			tmp = sparx_utilities.getvec(refa[i][0], refa[i][1])
			for j in range(3):
				refanorm[i][j] = tmp[j]
		del  refa, tmp
		transv = numpy.empty( (nima, 3), dtype = numpy.float32)
		if runtype=="ASSIGNMENT":
			for im in range(nima):
				trns = data[im].get_attr("xform.projection")
				for j in range(3):
					transv[im][j] = trns.at(2,j)
		else:
			# For REFINEMENT we have a problem, as the exact angle is known only after the next step of assigning projections.
			# So, we will assume it is the one with max peak
			for im in range(nima):
				qt = -1.0e23
				it = -1
				for iref in range(numref):
					pt = peaks[iref][im]
					if(pt > qt):
						qt = pt
						it = iref
				for j in range(3):
					transv[im][j] = trans[it][im].at(2,j)
		#  We have all vectors, now create a list of assignments of images to references
		refassign = [-1]*nima
		for im in range(nima):
			refassign[im] = abs(numpy.inner(refanorm,transv[im])).argmax()
		assigntorefa = mpi.mpi_gatherv(refassign, nima, mpi.MPI_INT, recvcount, disps, mpi.MPI_INT, main_node, mpi.MPI_COMM_WORLD)
		assigntorefa = list(map(int, assigntorefa))
		del refassign, refanorm, transv
		"""Multiline Comment79"""
		if myid == main_node:
			SA = False
			asi = [[] for iref in range(numref)]
			report_error = 0
			for imrefa in range(numrefang):
				pass#IMPORTIMPORTIMPORT from utilities import findall
				N = sparx_utilities.findall(imrefa, assigntorefa)
				current_nima = len(N)
				if( current_nima >= numref and report_error == 0):
					tasi = [[] for iref in range(numref)]
					maxasi = current_nima//numref
					nt = current_nima
					kt = numref
					K = list(range(numref))

					d = numpy.empty( (numref, current_nima), dtype = numpy.float32)
					for ima in range(current_nima):
						for iref in range(numref):  d[iref][ima] = dtot[iref][N[ima]]

					while nt > 0 and kt > 0:
						l = d.argmax()
						group = l//current_nima
						ima   = l-current_nima*group
						if SA:
							J = [0.0]*numref
							sJ = 0
							Jc = [0.0]*numref
							for iref in range(numref):
								J[iref] = numpy.exp(d[iref][ima]/T)
								sJ += J[iref]
							for iref in range(numref):
								J[iref] /= sJ
							Jc[0] = J[0]
							for iref in range(1, numref):
								Jc[iref] = Jc[iref-1]+J[iref]
							sss = numpy.random.random()
							for group in range(numref):
								if( sss <= Jc[group]): break
						tasi[group].append(N[ima])
						N[ima] = -1
						for iref in range(numref):  d[iref][ima] = -1.e10
						nt -= 1
						masi = len(tasi[group])
						if masi == maxasi:
							for im in range(current_nima):  d[group][im] = -1.e10
							kt -= 1
					else:
						for ima in range(current_nima):
							if N[ima] > -1:
								qm = -1.e10
								for iref in range(numref):
									qt = dtot[iref][N[ima]]
									if( qt > qm ):
										qm = qt
										group = iref
								tasi[group].append(N[ima])
					del d, N, K
					if  SA:  del J, Jc
					for iref in range(numref):
						asi[iref] += tasi[iref]
					del tasi
				else:
					report_error = 1
			#  This should be deleted only once we know that the number of images is sufficiently large, see below.
			del dtot
		else:
			assignment = []
			report_error = 0
		report_error = sparx_utilities.bcast_number_to_all(report_error, source_node = main_node)
		if report_error == 1:  sparx_global_def.ERROR('Number of images within a group too small', "mref_ali3d_MPI", 1, myid)
		if myid == main_node:
			assignment = [0]*total_nima
			for iref in range(numref):
				for im in range(len(asi[iref])):
					assignment[asi[iref][im]] = iref
			del asi
		
		"""Multiline Comment80"""

		"""Multiline Comment81"""
		assignment = mpi.mpi_scatterv(assignment, recvcount, disps, mpi.MPI_INT, recvcount[myid], mpi.MPI_INT, main_node, mpi.MPI_COMM_WORLD)
		assignment = list(map(int, assignment))
		#  compute number of particles that changed assignment and how many are in which group
		nchng = 0
		npergroup = [0]*numref
		for im in range(nima):
			iref = data[im].get_attr('group')
			npergroup[assignment[im]] += 1
			if( iref != assignment[im]): nchng += 1
			data[im].set_attr('group', assignment[im])
		nchng = mpi.mpi_reduce(nchng, 1, mpi.MPI_INT, mpi.MPI_SUM, 0, mpi.MPI_COMM_WORLD)
		npergroup = mpi.mpi_reduce(npergroup, numref, mpi.MPI_INT, mpi.MPI_SUM, 0, mpi.MPI_COMM_WORLD)
		npergroup = list(map(int, npergroup))
		terminate = 0
		if( myid == 0 ):
			ngroup=[]
			nchng = int(nchng[0])
			precn = 100*float(nchng)/float(total_nima)
			msg = " Number of particles that changed assignments %7d, percentage of total: %5.1f"%(nchng, precn)
			log.add(msg)
			msg = " Group       number of particles"
			log.add(msg)
			for iref in range(numref):
				msg = " %5d       %7d"%(iref+1, npergroup[iref])
				log.add(msg)
				ngroup.append(int(npergroup[iref]))
			if(precn <= termprec):  terminate = 1
		else:
			ngroup = 0
		terminate = mpi.mpi_bcast(terminate, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD)
		terminate = int(terminate[0])
		ngroup = sparx_utilities.wrap_mpi_bcast(ngroup,main_node)

		if runtype=="REFINEMENT":
			for im in range(nima):
				data[im].set_attr('xform.projection', trans[assignment[im]][im])
				pixer[0][im] = pixer[assignment[im]][im]
			pixer = pixer[0]
			if(center == -1):
				cs[0], cs[1], cs[2], dummy, dummy = sparx_utilities.estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)				
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f"%(cs[0], cs[1], cs[2])
					log.add(msg)
				cs = mpi.mpi_bcast(cs, 3, mpi.MPI_FLOAT, main_node, mpi.MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				sparx_utilities.rotate_3D_shift(data, cs)
			#output pixel errors
			recvbuf = mpi.mpi_gatherv(pixer, nima, mpi.MPI_FLOAT, recvcount, disps, mpi.MPI_FLOAT, main_node, mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			if(myid == main_node):
				recvbuf = list(map(float, recvbuf))
				pass#IMPORTIMPORTIMPORT from statistics import hist_list
				lhist = 20
				region, histo = sparx_statistics.hist_list(recvbuf, lhist)
				if(region[0] < 0.0):  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles"
				log.add(msg)
				for lhx in range(lhist):
					msg = " %10.3f      %7d"%(region[lhx], histo[lhx])
					log.add(msg)
				del region, histo
			del recvbuf
		fscc = [None]*numref
		if fourvar and runtype=="REFINEMENT":
			sumvol = sparx_utilities.model_blank(nx, nx, nx)
		start_time = time.time()

		if( not focus ):
			for im in range(nima):  data[im] = sparx_fundamentals.fft(data[im])
		highres = []
		lowpass_tmp = []
		tmpref =[]
		for iref in range(numref):
			#  3D stuff
			pass#IMPORTIMPORTIMPORT from time import localtime, strftime
			if Tracker["constants"]["3d-interpolation"]=="trl":
				chunk_id   = 0
				niter      = 10
				upweighted = False
				compensate = False
				verbose    =  0
				sign       = 1 
				volref0    = recons3d_n_trl_MPI_one_node(data, CTF, snr, sign, npad, sym, iref, niter, verbose, upweighted, compensate, chunk_id)
				chunk_id   = 1
				volref1    = recons3d_n_trl_MPI_one_node(data, CTF, snr, sign, npad, sym, iref, niter, verbose, upweighted, compensate, chunk_id)
				if myid == main_node:
					fscc[iref] = sparx_statistics.fsc(volref0, volref1)
					volref = volref1+volref0
					del volref1
					del volref0
				#chunk_id   = -1
				#volref = recons3d_n_trl_MPI_one_node(data, CTF, snr, sign, npad, sym, iref, niter, verbose, upweighted, compensate, chunk_id)
			else:
				if(CTF): volref, fscc[iref] = sparx_reconstruction.rec3D_two_chunks_MPI(data, snr, sym, mask3D,\
					os.path.join(outdir, "resolution_%02d_%04d"%(iref, total_iter)), myid, main_node, index=iref, npad=npad, finfo=frec)
				else:    volref, fscc[iref] = sparx_reconstruction.rec3D_MPI_noCTF(data, sym,mask3D,\
					os.path.join(outdir, "resolution_%02d_%04d"%(iref, total_iter)), myid, main_node, index=iref, npad=npad, finfo=frec)
			if(myid == 0):
				log.add( "Time to compute 3D: %d" % (time.time()-start_time) );start_time = time.time()
				volref.write_image(os.path.join(outdir, "vol%04d.hdf"%( total_iter)), iref)
				tmpref.append(volref)			
				if fourvar and runtype=="REFINEMENT": sumvol += volref
			#set_filter_parameters_from_adjusted_fsc(Tracker["constants"]["total_stack"],ngroup[iref],Tracker)
			if myid ==main_node:
				res = 0.5
				for ifreq in range(len(fscc[iref][0])-1,0,-1):
					if fscc[iref][1][ifreq] > 0.5 : # always use .5 as cutoff
						res = fscc[iref][0][ifreq]
						break

				Tracker["lowpass"] = min(0.45, res)
				#Tracker["lowpass"] = max(0.11, res)
				Tracker["falloff"] = 0.1
				log.add(" low pass filter is %f    %f   %d"%(Tracker["lowpass"], Tracker["falloff"], ngroup[iref]))
			else:
				Tracker["lowpass"] = 0.0
				Tracker["lowpass"] = 0.0
				res = 0.0
			res = sparx_utilities.bcast_number_to_all(res, main_node)
			highres.append(int(res*Tracker["nxinit"]+ 0.5))
			if myid ==main_node:
				log.add("%d highres                   %f"%(iref, highres[iref]))
			Tracker["lowpass"] = sparx_utilities.bcast_number_to_all(Tracker["lowpass"], main_node)
			Tracker["falloff"] = sparx_utilities.bcast_number_to_all(Tracker["falloff"], main_node)
			lowpass_tmp.append(Tracker["lowpass"])
		Tracker["lowpass"] = min(lowpass_tmp)
		if myid ==main_node:
			log.add(" the adopted  low pass filter is %f    %f   "%(Tracker["lowpass"], Tracker["falloff"]))
		for iref in range(numref):
			
			if myid ==main_node:
				volref =tmpref[iref]
				if Tracker["mask3D"]: 
				
					mask3D = sparx_utilities.get_im(Tracker["mask3D"])
					stat = EMAN2_cppwrap.Util.infomask(volref, mask3D, False)
					volref -= stat[0]
					EMAN2_cppwrap.Util.mul_scalar(volref, 1.0/stat[1])
					
				if(Tracker["constants"]["PWadjustment"]):
				
					rt = sparx_utilities.read_text_file(Tracker["PW_dict"][Tracker["constants"]["nxinit"]])
					ro = sparx_fundamentals.rops_table(volref)
					for i in range(1,len(ro)):  ro[i] = (rt[i]/ro[i])**Tracker["constants"]["upscale"]
					volref =sparx_filter.filt_table(volref,ro)
							
				if (Tracker["constants"]["low_pass_filter"]==-1.):  volref = sparx_filter.filt_tanl(volref, Tracker["lowpass"], Tracker["falloff"])                                       # low pass from resolution 
				else:                                               volref = sparx_filter.filt_tanl(volref, min(Tracker["constants"]["low_pass_filter"]/Tracker["shrinkage"],0.45), Tracker["falloff"]) # user define filter
			
	
				if Tracker["mask3D"]: EMAN2_cppwrap.Util.mul_img(volref, mask3D)
				volref.write_image( os.path.join(outdir,"volf%04d.hdf"%(total_iter)), iref)
				del volref
				
		"""Multiline Comment82"""
		#if(myid == main_node):
		"""Multiline Comment83"""
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		if terminate == 0: # headers are only updated when the program is going to terminate
			start_time = time.time()
			if nrefine!=0: par_str = ['xform.projection', 'ID', 'group']
			else: par_str = ['group', 'ID' ]
			"""Multiline Comment84"""
			if(myid == 0):log.add( "Time to write headers: %d\n" % (time.time()-start_time) );start_time = time.time()
		else:
			final_list = sparx_utilities.get_sorting_params_refine(Tracker, data, total_nima)
			group_list, ali3d_params_list = sparx_utilities.parsing_sorting_params(final_list)
			if myid ==main_node:
				group_list_saved_file =os.path.join(outdir, "list2.txt")
				sparx_utilities.write_text_file(group_list,group_list_saved_file)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			Tracker["this_partition"]=group_list
			break
	if terminate != 1:
		final_list = sparx_utilities.get_sorting_params_refine(Tracker,data, total_nima)
		group_list, ali3d_params_list = sparx_utilities.parsing_sorting_params(final_list)  
		if myid ==main_node:
			group_list_saved_file =os.path.join(outdir, "list2.txt")
			sparx_utilities.write_text_file(group_list,group_list_saved_file)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		Tracker["this_partition"]=group_list
	### program finishes
	if myid ==main_node:  log.add(" mref_ali3d_EQ_Kmeans finishes !")
		
######
 
