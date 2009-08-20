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

#  This file contains code under development or not currently used.

def ali2d_g(stack, outdir, maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", center=-1, maxit=0, \
		CTF=False, snr=1.0, Fourvar = False, user_func_name="ref_ali2d", rand_alpha = False, MPI=False):
	if MPI:
		ali2d_g_MPI(stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, center, maxit, CTF, snr, Fourvar, user_func_name, rand_alpha)
		return

	from utilities    import model_circle, drop_image, get_image, get_input_from_string
	from statistics   import fsc_mask, sum_oe
	from alignment    import Numrinit, ringwe
	from development  import ali2d_single_iter_g
	from filter       import filt_ctf, filt_table, filt_tophatb
	from utilities    import print_begin_msg, print_end_msg, print_msg
	from fundamentals import fft, rot_avg_table, fshift, prepi
	from utilities    import write_text_file, file_type
	from applications import MPI_start_end
	import os
		
	print_begin_msg("ali2d_g")

	if os.path.exists(outdir):   ERROR('Output directory exists, please change the name and restart the program', " ", 1)
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
	print_msg("Data with CTF               : %s\n"%(CTF))
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
		if ima.get_attr_default('ctf_applied', 2) > 0:	ERROR("data cannot be ctf-applied", "ali2d_c_MPI", 1)
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
	if CTF:  ctf_2_sum += 1.0/snr  # note this is complex addition (1.0/snr,0.0)
	# startup
	numr = Numrinit(first_ring, last_ring, rstep, mode) 	#precalculate rings
 	wr = ringwe(numr, mode)
	
	data_prep = []
	for im in xrange(nima):
		img_prep, kb = prepi(data[im])
		data_prep.append(img_prep)

	# initialize data for the reference preparation function
	#  mask can be modified in user_function
	ref_data = []
	ref_data.append( mask )
	ref_data.append( center )
	ref_data.append( None )
	ref_data.append( None )
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
			if  Fourvar:  
				tavg, ave1, ave2, vav, sumsq = add_ave_varf(data, mask, "a", CTF, ctf_2_sum)
				# write the current average
				drop_image(fft(tavg), os.path.join(outdir, "aqc_%03d.hdf"%(total_iter)))
				tavg    = fft(Util.divn_img(tavg, vav))

				vav_r	= Util.pack_complex_to_real(vav)
				drop_image(vav_r, os.path.join(outdir, "varf_%03d.hdf"%(total_iter)))
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
				write_text_file(frsc, os.path.join(outdir, "resolution%03d"%(total_iter)) )
			else:
				ave1, ave2 = sum_oe(data, "a", CTF, ctf_2_sum)
				if CTF:  tavg = fft(Util.divn_img(fft(Util.addn_img(ave1, ave2)), ctf_2_sum))
				else:	 tavg = (ave1+ave2)/nima
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
			drop_image(tavg, os.path.join(outdir, "aqf_%03d.hdf"%(total_iter)))

			# a0 should increase; stop algorithm when it decreases.
			if Fourvar:  
				Util.div_filter(sumsq, vav)
				sumsq = filt_tophatb(sumsq, 0.01, 0.49)
				a1 = Util.infomask(sumsq, None, True)
				a1 = a1[0]
			else:
				a1 = tavg.cmp("dot", tavg, dict(negative = 0, mask = ref_data[0]))
			msg = "Iteration   #%5d	     criterion = %20.7e\n"%(total_iter,a1)
			print_msg(msg)
			if total_iter == len(xrng)*max_iter: break
			if a1 < a0:
				if auto_stop == True: break
			else:	a0 = a1
			sx_sum, sy_sum = ali2d_single_iter_g(data, data_prep, kb, numr, wr, cs, tavg, cnx, cny, xrng[N_step], yrng[N_step], step[N_step], mode, CTF)
			
	drop_image(tavg, os.path.join(outdir, "aqfinal.hdf"))
	# write out headers
	from utilities import write_headers
	write_headers(stack, data, list_of_particles)
	print_end_msg("ali2d_g")


def ali2d_g_MPI(stack, outdir, maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", center=-1, maxit=0, CTF=False, snr=1.0, \
			Fourvar = False, user_func_name="ref_ali2d", rand_alpha=False):

	from utilities    import model_circle, model_blank, drop_image, get_image, get_input_from_string
	from utilities    import reduce_EMData_to_root, bcast_EMData_to_all, send_attr_dict, file_type, bcast_number_to_all, bcast_list_to_all
	from statistics   import fsc_mask, sum_oe, add_ave_varf_MPI
	from alignment    import Numrinit, ringwe, ali2d_single_iter
	from development  import ali2d_single_iter_g
	from filter       import filt_table, filt_ctf, filt_tophatb
	from numpy        import reshape, shape
	from fundamentals import fshift, fft, rot_avg_table, prepi
	from utilities    import write_text_file
	from utilities    import print_msg, print_begin_msg, print_end_msg
	from applications import MPI_start_end
	import os
	import sys
	from mpi 	  import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	  import mpi_reduce, mpi_bcast, mpi_barrier
	from mpi 	  import MPI_SUM, MPI_FLOAT, MPI_INT

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	ftp = file_type(stack)

	if myid == main_node:
		print_begin_msg("ali2d_g_MPI")
		if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', " ", 1)
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
		nima =0
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
		print_msg("Data with CTF               : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		if auto_stop:
			print_msg("Stop iteration with         : criterion\n")
		else:
			print_msg("Stop iteration with         : maxit\n")
		print_msg("User function               : %s\n"%(user_func_name))
		print_msg("Number of processors used   : %d\n"%(number_of_proc))

		
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
		if ima.get_attr_default('ctf_applied', 2) > 0:	ERROR("data cannot be ctf-applied", "ali2d_c_MPI", 1)
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
	
	for im in xrange(len(data)):
		data[im].set_attr('ID', list_of_particles[im])
		st = Util.infomask(data[im], mask, False)
		data[im] -= st[0]
		if CTF:
			ctf_params = data[im].get_attr("ctf")
			Util.add_img2(ctf_2_sum, ctf_img(nx, ctf_params))
	if CTF:
		reduce_EMData_to_root(ctf_2_sum, myid, main_node)

		ctf_2_sum += 1.0/snr # this is complex addition (1.0/snr,0)
	else:  ctf_2_sum = None
	# startup
	numr = Numrinit(first_ring, last_ring, rstep, mode) 	#precalculate rings
 	wr = ringwe(numr, mode)
	
	data_prep = []
	for im in xrange(image_start, image_end):
		img_prep, kb = prepi(data[im-image_start])
		data_prep.append(img_prep)

	if myid == main_node:
		# initialize data for the reference preparation function
		ref_data = []
		ref_data.append( mask )
		ref_data.append( center )
		ref_data.append( None )
		ref_data.append( None )
		sx_sum = 0.0
		sy_sum = 0.0
		a0 = -1.0e22

	again = True
	total_iter = 0
	cs = [0.0]*2

	for N_step in xrange(len(xrng)):
		msg = "\nX range = %5.2f   Y range = %5.2f   Step = %5.2f\n"%(xrng[N_step], yrng[N_step], step[N_step])
		if myid == main_node: print_msg(msg)
		for Iter in xrange(max_iter):
			total_iter += 1
			if Fourvar:  
				tavg, ave1, ave2, vav, sumsq = add_ave_varf_MPI(myid, data, mask, "a", CTF, ctf_2_sum)
			else:
				ave1, ave2 = sum_oe(data, "a", CTF, EMData())  # pass empty object to prevent calculation of ctf^2
				reduce_EMData_to_root(ave1, myid, main_node)
				reduce_EMData_to_root(ave2, myid, main_node)

			if myid == main_node:

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
					write_text_file(frsc, os.path.join(outdir, "resolution%03d"%(total_iter)))
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
				msg = "Iteration   #%5d	     criterion = %20.7e\n"%(total_iter, a1)
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
			cs = [float(cs[0]), float(cs[1])]
			if auto_stop:
				again = mpi_bcast(again, 1, MPI_INT, main_node, MPI_COMM_WORLD)
				if not again: break
			if total_iter != max_iter*len(xrng):
				if total_iter >= max_iter*len(xrng)-5: 
					sx_sum, sy_sum = ali2d_single_iter_g(data, data_prep, kb, numr, wr, cs, tavg, cnx, cny, xrng[N_step], yrng[N_step], step[N_step], mode, CTF)
				else:
					sx_sum, sy_sum = ali2d_single_iter(data, numr, wr, cs, tavg, cnx, cny, xrng[N_step], yrng[N_step], step[N_step], mode, CTF)
				sx_sum = mpi_reduce(sx_sum, 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
				sy_sum = mpi_reduce(sy_sum, 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)

	if myid == main_node:  drop_image(tavg, os.path.join(outdir, "aqfinal.hdf"))
	# write out headers  and STOP, under MPI writing has to be done sequentially
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
	if myid == main_node: print_end_msg("ali2d_g_MPI")


def ali2d_single_iter_g(data, data_prep, kb, numr, wr, cs, tavg, cnx, cny, xrng, yrng, step, mode, CTF=False, ali_params="xform.align2d"):
	"""
		single iteration of 2D alignment using ormy3g
		if CTF = True, apply CTF to data (not to reference!)
	"""
	from utilities import combine_params2, inverse_transform2, get_params2D, set_params2D
	from alignment import Applyws
	from development import ormy3g
	from fundamentals import prepi

	if CTF:
		from filter  import filt_ctf

	# 2D alignment using rotational ccf in polar coords and NUFFT interpolation
	tavgi, kb = prepi(tavg)
	cimage = Util.Polar2Dmi(tavgi, cnx, cny, numr, mode, kb)
	Util.Frngs(cimage, numr)
	Applyws(cimage, numr, wr)

	maxrin = numr[-1]
	sx_sum = 0.0
	sy_sum = 0.0
	for im in xrange(len(data)):
		if CTF:
			#Apply CTF to image
			ctf_params = data[im].get_attr("ctf")
			ima = filt_ctf(data_prep[im], ctf_params, True)
		else:
			ima = data_prep[im].copy()
		alpha, sx, sy, mirror, dummy = get_params2D(data[im], ali_params)
		alpha, sx, sy, mirror        = combine_params2(alpha, sx, sy, mirror, 0.0, -cs[0], -cs[1], 0)
		alphai, sxi, syi, scalei     = inverse_transform2(alpha, sx, sy, 1.0)

		# align current image to the reference
		[angt, sxst, syst, mirrort, peakt] = ormy3g(ima, kb, cimage, xrng, yrng, step, mode, numr, cnx+sxi, cny+syi)
		# combine parameters and set them to the header, ignore previous angle and mirror
		[alphan, sxn, syn, mn] = combine_params2(0.0, -sxi, -syi, 0, angt, sxst, syst, mirrort)
		set_params2D(data[im], [alphan, sxn, syn, mn, 1.0], ali_params)

		if mn == 0: sx_sum += sxn
		else:       sx_sum -= sxn
		sy_sum += syn

	return sx_sum, sy_sum


def ali_param_3d_to_2d(stack):
	from utilities import get_arb_params, set_arb_params, ali_params_3D_2D
	import types 
	parameters_2D_name=["alpha", "sx",    "sy",  "mirror" ]
	parameters_3D_name=[ "phi",  "theta", "psi", "s2x", "s2y" ]
	if(type(stack) == types.StringType): 
		nima = EMUtil.get_image_count(stack)
		for im in xrange(nima):
			if( type(stack) == types.StringType ): ima.read_image(stack,im)
			parameters_3D = get_arb_params(ima, parameters_3D_name)
			alpha, sx, sy, mirror = ali_params_3D_2D(parameters_3D[0],parameters_3D[1],parameters_3D[2],parameters_3D[3],parameters_3D[4],parameters_3D[5])
			set_arb_params( ima, [ alpha, sx, sy, mirror ], parameters_2D_name)
			ima.write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)
	else:
		for im in xrange(len(stack)):
			parameters_3D = get_arb_params(stack[im], parameters_3D_name)
			alpha, sx, sy, mirror = ali_params_3D_2D(parameters_3D[0],parameters_3D[1],parameters_3D[2],parameters_3D[3],parameters_3D[4],parameters_3D[5])
			set_arb_params( stack[im], [ alpha, sx, sy, mirror ], parameters_2D_name)

def ali_param_2d_to_3d(stack):
	import os
	from utilities import compose_transform2
	from utilities import get_arb_params, set_arb_params
	
	""" 
		transform 2D alignment parameters into 3D alignment parameters.	
		has to be used if no 3D aligment parameters have been determined yet.

	"""
	dummy = EMData()
	header_3d_ID=["alpha","sx","sy","mirror", "psi","theta","phi","s2x","s2y"]
	header_2d_ID=["alpha","sx","sy","mirror"]
	nima = EMUtil.get_image_count(stack)
	for i in xrange(nima):
		dummy.read_image(stack,i, True)
		param_2d=get_arb_params(dummy,header_2d_ID)
		alphan,s2x, s2y,scalen = compose_transform2(0,param_2d[1], param_2d[2], 1, -param_2d[0], 0, 0, 1)
		if( param_2d[3] > 0.5 ) :
			phi   = 180.0
			theta = 180.0
			psi   = (540.0 + alphan)%360.0
		else:
			psi   = (alphan)%360.0
			phi   = 0.0
			theta = 0.0
		params=[ 0., 0., 0., 0, psi, theta, phi, s2x, s2y]
		set_arb_params(dummy, params, header_3d_ID)
		dummy.write_image(stack,i, EMUtil.ImageType.IMAGE_HDF, True)


'''	NOT USED
def ali2d_r(stack, outdir, ir=1, ou=-1, rs=1, center=1, maxit=10):
# 2D alignment using rotational ccf in polar coords and quadratic interpolation
	from utilities import model_circle, compose_transform2, combine_params2, drop_image, start_time, finish_time
	from fundamentals import rot_shift2D
	from statistics import aves, fsc, ave_oe_series
	from alignment import Numrinit, ringwe, Applyws, ang_n
	from filter import fshift
	import os
	from EMAN2 import EMUtil,EMData
	tt =start_time()
	first_ring=int(ir)
	last_ring=int(ou)
	rstep=int(rs)
	max_iter = int(maxit)
	# create the output directory, if it does not exist
	if os.path.exists(outdir) is False: os.mkdir(outdir)
	# reading image and gethering the information
	nima = EMUtil.get_image_count(stack)
	ima = EMData()
	ima.read_image(stack,0)
	nx=ima.get_xsize()
	# default value for the last ring
	if(last_ring==-1): last_ring=nx-2
	mask = model_circle(last_ring,nx,nx)
	# calculate total average using current alignment parameters - 'aves' function removed
	tave = aves(stack)
	a0 = tave.cmp("dot", tave, {"negative":0,"mask":mask})
	print  "Initial criterion=",a0
	#  center is in SPIDER convention
	cnx = int(nx/2)+18
	cny = cnx
	mode = "F"
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr   = ringwe(numr, mode)
	again = True
	cs = [0.0]*2  # additional shifts from the centering of the average (if requested) from the previous iteration
	#iteration
	for iter in xrange(max_iter):
		#prepare the reference
		cimage=Util.Polar2Dm(tave, cnx, cny, numr, mode)
		Util.Frngs(cimage, numr)
		Applyws(cimage, numr, wr)
		# prepare the images and align them
		for im in xrange(nima):
			# get image & params.
			ima = EMData()
			ima.read_image(stack,im)
			alpah = ima.get_attr('alpha')
			sx =  ima.get_attr('sx')
			sy =  ima.get_attr('sy')
			mirror =  ima.get_attr('mirror')
			#  add centering from the previous iteration, if image was mirrored, filp the sign of x center
			if mirror: alpha, sx, sy, scale = compose_transform2(alpha, sx, sy, 1.0, 0.0,  cs[0], -cs[1], 1.0)
			else:      alpha, sx, sy, scale = compose_transform2(alpha, sx, sy, 1.0, 0.0, -cs[0], -cs[1], 1.0)
			# align current image to the reference
			peak = -1.0E23
			cimage_o=Util.Polar2Dm(ima, cnx-sx, cny-sy, numr, mode)
			Util.Frngs(cimage_o, numr)
			retvals = Util.Crosrng_ms(cimage, cimage_o, numr)
			qn = retvals["qn"]
			qm = retvals["qm"]
			if(qn >=peak or qm >=peak):
				if (qn >= qm):
					angt = ang_n(retvals["tot"], mode, numr[len(numr)-1])
					peak=qn
					mirror=0
				else:
					angt = ang_n(retvals["tmt"], mode, numr[len(numr)-1])
					peak=qm
					mirror=1
			angt = ang_n(retvals["tot"], mode, numr[len(numr)-1])		
			# combine parameters and set them to the header, ignore previous angle and mirror
			[alphan, sxn, syn, mn] = combine_params2(0.0, sx, sy, 0, angt, 0, 0, mirror)
			#print  "combine with previous ",iter,im,psin, sxn, syn, mn
			ima.set_attr_dict({'alpha':alphan, 'sx':sxn, 'sy':syn,'mirror':mn})
			ima.write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)
			
		av1,av2 = ave_oe_series(stack) # odd, even averages.
		frsc = fsc(av1, av2, 1.0, os.path.join(outdir, "drc%03d"%iter))
		if center:
			# Center.
			av1 = (av1*(int(nima/2)+(nima%2)) + av2*int(nima/2))/nima
			# Potentially replace this:
			cs = av1.phase_cog()
			#print  "SHIFT ",cs
			tave = fshift(av1, -cs[0], -cs[1])
		else:
			tave = (av1*(int(nima/2)+(nima%2)) + av2*int(nima/2))/nima
		# a0 should increase; stop algorithm when it decreases.
		a0 = tave.cmp("dot", tave, {"negative":0,"mask":mask}) # tave' * tave
		print'%-15 %5d %15 = %10.g'%(" ITERATION #",iter+1,"criterion",a0)
		# write the current average
		drop_image(tave, os.path.join(outdir, "aqc%03d.spi"%iter))
	print finish_time(tt)
'''
'''	Not used anywhere
def ali3d_b(stack, ref_vol, outdir, maskfile = None,ir=1, ou=-1, rs=1, xr=0, yr=0,
	    ts=1, dtheta=15, maxit=10, ref_a="S", sym="c1", CTF = None):  
    # 2D alignment using rotational ccf in polar coords and gridding
    # interpolation 
    
	from utilities      import model_circle, drop_image
	from alignment      import proj_ali
	from reconstruction import recons3d_4nn_ctf, recons3d_4nn
	from statistics     import fsc
	from filter         import filt_params, filt_btwl
	import os
	snr=1.
	sign=1

	first_ring = int(ir);    last_ring = int(ou); rstep = int(rs); 
	xrng       = int(xr);    yrng      = int(yr); step  = int(ts);
	max_iter   = int(maxit); 
	# create the output directory, if it does not exist
	if os.path.exists(outdir) is False: os.mkdir(outdir)
	nima = EMUtil.get_image_count(stack)
	#  compare the image size
	vol = EMData()
	vol.read_image(stack,0)
	nx = vol.get_xsize()
	#vol.read_image(ref_vol)
	# default value for the last ring
	if(last_ring==-1): last_ring=nx-2
	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask3D=get_image(maskfile)
		else: mask3D = maskfile
	else: mask3D = model_circle(last_ring, nx, nx, nx)
	#mask2D = model_circle(last_ring, nx, nx, 1)
    
	# do the projection matching
	# begin a refinement loop, slowly decrease dtheta inside the loop
	for iter in xrange(max_iter):
		print " ali3d_b: ITERATION #",iter
		proj_ali(ref_vol, mask3D, stack, first_ring, last_ring, rstep, xrng, yrng, step, dtheta, ref_a, sym, CTF)

		# resolution
		list_p = range(0,nima,2)
		if(CTF): vol1 = recons3d_4nn_ctf(stack, list_p, snr, sign, sym)
		else:    vol1 = recons3d_4nn(stack, list_p, sym)

		list_p = range(1,nima,2)
		if(CTF): vol2 = recons3d_4nn_ctf(stack, list_p, snr, sign, sym)
		else:    vol2 = recons3d_4nn(stack, list_p, sym)

		[mean, sigma, xmin, xmax ] =  Util.infomask(vol1, mask3D, false)
		vol1 = (vol1-mean)*mask3D                                    #use processor instead
		[mean, sigma, xmin, xmax ] =  Util.infomask(vol2, mask3D)
		vol2 = (vol2-mean)*mask3D

		res = fsc(vol1, vol2, 1.0, os.path.join(outdir,"dres%03d"%iter))
		del vol1
		del vol2
		
		# calculate new and improved 3D
		list_p = range(nima)
		if(CTF): vol = recons3d_4nn_ctf(stack, list_p, snr, sign, sym)
		else:    vol = recons3d_4nn(stack, list_p, sym)
		# store the reference volume
		drop_image(vol,os.path.join(outdir,"vol%03d.spi"%iter))

		# here figure the filtration parameters and filter vol for the
		# next iteration
		fl, fh = filt_params(res)
		vol    = filt_btwl(vol, fl, fh)
		# store the filtred reference volume
		drop_image(vol, os.path.join(outdir, "volf%03d.spi"%iter))
		#  here the loop should end
	return vol
'''

''' 	        
def ali3d_e_iter(stack, mask3d, iteration, defocus, ptl_defgrp, radius, lowf, highf, snr, dtheta, symmetry):

	"""
		3D alignment inner loop 
	"""
	
	import os
	from statistics import fsc
	from filter     import filt_btwl, filt_ctf
	from projection import prep_vol
	from utilities  import amoeba, info
	from math       import pi,sin
	
	img=EMData()	
	firstp=0
	lastp = 0
	sign=1.
	par_str=["s2x", "s2y", "psi", "theta", "phi"]
	par_str1=["psi", "theta", "phi", "s2x", "s2y"]
	fdres_pat="dres{***}"
	nl = EMUtil.get_image_count(stack)
	img.read_image(stack,0)
	ctf_params = img.get_attr('ctf')
	Pixel_size = ctf_params.apix
	cs = ctf_params.cs
	voltage = ctf_params.voltage
	for mic in range(len(ptl_defgrp)):
		iteration += 1
		lastp += ptl_defgrp[mic]
		list_of_particles=range(firstp,lastp)
		# Calculate reconstruction using particles NOT in the current group					  	  
		list_p = range(nl)
		for e in list_of_particles:  list_p.remove(e)								  	  
		vol = recons3d_4nn_ctf(stack,list_p, snr, sign, symmetry)
		ctf_params.defocus = defocus[mic]
		vol_ctf = filt_ctf(vol, ctf_params)
		vol = filt_btwl(vol_ctf, lowf, highf)
		vol*=mask3d
		data=[]
		volft,kb=prep_vol(vol)
		data.insert(0,volft)
		data.insert(1,kb)
		data.insert(3,mask3d)
		for imn in range(ptl_defgrp[mic]):
			img.read_image(stack, list_of_particles[imn])
			data.insert(2,img)
			atparams = get_arb_params(img,par_str1)
			initial = eqproj(atparams,data)
			stheta  = abs(sin(atparams[1]*pi/180.0))
			if(stheta >0.0): weight_phi=min(180.0,dtheta/stheta/2.0)
			else: weight_phi=180.0
			from utilities import start_time, finish_time
			t3=start_time()
			optm_params =  amoeba(atparams, [weight_phi,dtheta,weight_phi,1.0,1.0],eqproj, 1.e-4,1.e-4,500,data)
			del  data[2]
			tmp_params = [optm_params[0][3],optm_params[0][4],optm_params[0][0],optm_params[0][1],optm_params[0][2]]
			update_attrt(img,par_str,tmp_params,"r")
			t4 = finish_time(t3)
			img.write_image(stack, list_of_particles[imn], EMUtil.ImageType.IMAGE_HDF, True)
		firstp=lastp
		# after all groups have been processed, calculate the current resolution
		# for the time being, do it after each iteration
		list_p = range(0,nl,2)
		vol1=recons3d_4nn_ctf(stack, list_p, snr, sign, symmetry)
		#drop_image(vol1,"Cnvolcmm1.spi")
		list_p = range(1,nl,2)
		vol2=recons3d_4nn_ctf(stack, list_o, snr, sign, symmetry)
		vol1 *= mask3d
		vol2 *= mask3d
		fdres = parse_spider_fname(fdres_pat, iteration)							  	  
		fsc(vol1, vol2, 1.0, fdres)		 										  	  											  	  
	return  iteration

	
def update_attrt(img, par_str, params, opt="s"):	
	"""
		Update attributes stored in the headers
		1. subtract new
		2. replaced by new
		3. add new
	"""
	old_params = get_arb_params(img,par_str)
	for i in xrange(len(par_str)):
		if(opt == "s"): old_params[i]-=params[i]
		if(opt == "r"): old_params[i]=params[i]
		if(opt == "a"): old_params[i]+=params[i]
	set_arb_params(img, old_params, par_str)
	del old_params


def ali3d_e(stack, maskfile, radius=25, snr=1.0, lowf=.25, highf=.3, max_it=10, dtheta=2, symmetry="c1", CTF = None):
	"""
		
	"""
	nl = EMUtil.get_image_count(stack) # check the data, and get defocus
	img = EMData()
	defocus = []
	ptl_defgrp = []  	 
	par_str = ["defocus","defgrp"]
	ptl_num_def = 0
	def_num = -1
	for i in xrange(nl):
		img.read_image(stack,i)		 
		tmp_params = get_arb_params(img, par_str)
		if(i==0):
			nx = img.get_xsize()		 
			defocus.append(tmp_params[0])
			def_num += 1
		if(tmp_params[1]!=def_num): 
			defocus.append(tmp_params[0])
			ptl_defgrp.append(ptl_num_def) 
			def_num+=1
			ptl_num_def=0
		ptl_num_def+=1
	ptl_defgrp.append(ptl_num_def)
	iteration=-1
	if os.path.exists(maskfile): mask3D=get_image(maskfile)
	else : mask3d=model_circle(radius, nx, nx, nx)
	for dummy in xrange(1000): # fake loop
		if(iteration <= max_it): iteration=ali3d_e_iter(stack,mask3d,iteration,defocus,ptl_defgrp,radius,lowf,highf,snr,dtheta,symmetry)		
'''

"""
This is not used, mainly because it is not parallelizable.
def ali2d_b(stack, outdir, maskfile=None, ir=1, ou=-1, rs=1, xr=0, yr=0, ts=1, center=1, maxit=10):
# 2D alignment using rotational ccf in polar coords and gridding interpolation

	import os
	from utilities import model_circle, compose_transform2, combine_params2, drop_image, info
	from statistics import aves, fsc, ave_var_series_g, ave_oe_series_g
	from alignment import Numrinit, ringwe, Applyws, ormqip, kbt, ormy3
	from fundamentals import prepg, rtshgkb, fshift
	from filter import filt_from_fsc_bwt, filt_table
	#from string import *

	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); xrng=int(xr); yrng=int(yr); step=int(ts); max_iter=int(maxit);
	
	if os.path.exists(outdir) : os.system('rm -rf '+outdir)
	os.mkdir(outdir)
	
	nima = EMUtil.get_image_count(stack)
	temp = EMData()
	temp.read_image(stack,0)

	# prepare KB interpolants
	nx=temp.get_xsize()
	kb = kbt(nx)

	# default value for the last ring
	#if(last_ring==-1): last_ring=nx-2
	if last_ring == -1 : last_ring=int(nx/2)-2

	# read images and prepare them for gridding
	data = []
	for im in xrange(nima):
		if(im>0):
			temp = EMData()
			temp.read_image(stack,im)
		alpha = temp.get_attr('alpha')
		sx =  temp.get_attr('sx')
		sy =  temp.get_attr('sy')
		mn =  temp.get_attr('mirror')
		data.append(prepg(temp,kb))
		data[im].set_attr_dict({'alpha':alpha, 'sx':sx, 'sy':sy, 'mirror': mn})
	if maskfile==None : mask = model_circle(last_ring,nx,nx)
	else              : get_image(maskfile)
	
	# calculate total average using current alignment parameters
	tave,tvar = ave_var_series_g(data,kb)
	#drop_image(tave,"a1.spi")
	#drop_image(tvar,"a2.spi")
	a0 = tave.cmp("dot", tave, {"negative":0, "mask":mask})
	print  "Initial criterion=",a0
	
	# do the alignment
	# IMAGES ARE SQUARES!
	#  center is in SPIDER convention
	cnx = int(nx/2)+1
	cny = cnx

	mode = "F"
	#precalculate rings
	numr=Numrinit(first_ring,last_ring,rstep,mode)
	wr=ringwe(numr,mode)
	
	for iter in xrange(max_iter):
		again = False
		for im in xrange(nima):
			# subtract current image from the average
			alpha = data[im].get_attr('alpha')
			sx =  data[im].get_attr('sx')
			sy =  data[im].get_attr('sy')
			mirror =  data[im].get_attr('mirror')
			temp = rtshgkb(data[im], alpha, sx, sy, kb)
			if  mirror: temp.process_inplace("mirror",{"axis":'x'})

			#  Subtract current image from the average
			refim = tave - temp/nima
			refim_p = prepg(refim,kb)
			cimage=Util.Polar2Dmi(refim_p,cnx,cny,numr,mode,kb)
			Util.Frngs(cimage, numr)
			Applyws(cimage, numr, wr)

			# align current image to the reference minus this image
			[angt,sxst,syst,mirrort,peakt]=ormy3(temp,cimage,xrng,yrng,step,mode,numr,cnx,cny,"gridding")			
			# print angt, sxst, syst, mirrort, peakt
			#[angt,sxst,syst,mirrort,peakt]=ormqip(prepg(temp,kb),cimage,xrng,yrng,step,mode,numr,kb,cnx,cny,nx)			
			# print angt, sxst, syst, mirrort, peakt

			# combine parameters and set them to the header
			[angn,sxn,syn,mn]=combine_params2(alpha, sx, sy, mirror, angt, sxst, syst, mirrort)

			# apply params to the image
			temp = rtshgkb(data[im], angn, sxn, syn, kb)
			if  mn: temp.process_inplace("mirror",{"axis":'x'})

			# check whether the criterion actually increased
			# add current image to the average
			temp = refim + temp/nima

			# calculate the criterion
			a1 = temp.cmp("dot", temp, {"negative":0, "mask":mask})
			print a1, angn, sxn, syn, mn
			if a1 > a0:
				#print  im,"  ",a1,"  ",mirror,"  ",mirrort,"  ",psi,"  ",angt,"  ",sxst,"  ",syst
				# replace the average by the improved average and set the new parameters to the image, otherwise, do nothing
				tave = temp.copy()
				data[im].set_attr_dict({'alpha':angn, 'sx':sxn, 'sy':syn,'mirror': mn})
				a0 = a1
				again = True
		if again :
			# calculate total average using current alignment parameters
			print " ITERATION #",iter,"  criterion = ",a0
			av1,av2 = ave_oe_series_g(data,kb)
			frsc = fsc(av1,av2,1.0,os.path.join(outdir,"drm%03d"%iter))
			if center:
				#center
				av1 = (av1*(int(nima/2)+(nima%2)) + av2*int(nima/2))/nima
				s = av1.phase_cog()
				tave = fshift(av1, -s[0], -s[1])
				for im in xrange(nima):
					alpha = data[im].get_attr('alpha')
					sx =  data[im].get_attr('sx')
					sy =  data[im].get_attr('sy')
					alphan,sxn,syn,scale = compose_transform2(alpha, sx, sy, 1.0, 0.0, -s[0], -s[1], 1.0)
					data[im].set_attr_dict({'alpha':alphan, 'sx':sxn, 'sy':syn})
			else:
				tave = (av1*(int(nima/2)+(nima%2)) + av2*int(nima/2))/nima
			# write the current average
			# fil_table=filt_from_fsc_bwt(frsc, low=0.5)
			# tave=filt_table(tave,fil_table)
			drop_image(tave,os.path.join(outdir,"aqm%03d.spi"%iter),"s")
		else:
			break

	for im in xrange(nima):
		temp.read_image(stack,im)
		alpha = data[im].get_attr('alpha')
		sx =  data[im].get_attr('sx')
		sy =  data[im].get_attr('sy')
		mn =  data[im].get_attr('mirror')
		temp.set_attr_dict({'alpha':alpha, 'sx':sx, 'sy':sy,'mirror': mn})
		temp.write_image(stack,im)
	#for im in xrange(nima): data[im].write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)
	return	

def ali2d_b2(stack, outdir, maskfile=None, ir=1, ou=-1, rs=1, xr=0, yr=0, ts=1, center=1, maxit=10):
# 2D alignment using rotational ccf in polar coords and gridding interpolation

	import os
	from utilities import model_circle, compose_transform2, combine_params2, drop_image, info
	from statistics import aves, fsc, ave_var_series_g, ave_oe_series_g
	from alignment import Numrinit, ringwe, Applyws, ormqip, kbt, ormy3
	from fundamentals import prepg, rtshgkb, fshift
	from filter import filt_from_fsc_bwt, filt_table
	#from string import *

	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); xrng=int(xr); yrng=int(yr); step=int(ts); max_iter=int(maxit);
	
	if os.path.exists(outdir) : os.system('rm -rf '+outdir)
	os.mkdir(outdir)
	
	nima = EMUtil.get_image_count(stack)
	temp = EMData()
	temp.read_image(stack,0)

	# prepare KB interpolants
	nx=temp.get_xsize()
	kb = kbt(nx)

	# default value for the last ring
	#if(last_ring==-1): last_ring=nx-2
	if last_ring == -1 : last_ring=int(nx/2)-2

	# read images and prepare them for gridding
	data = []
	for im in xrange(nima):
		if(im>0):
			temp = EMData()
			temp.read_image(stack,im)
		alpha = temp.get_attr('alpha')
		sx =  temp.get_attr('sx')
		sy =  temp.get_attr('sy')
		mn =  temp.get_attr('mirror')
		data.append(prepg(temp,kb))
		data[im].set_attr_dict({'alpha':alpha, 'sx':sx, 'sy':sy, 'mirror': mn})
	if maskfile==None : mask = model_circle(last_ring,nx,nx)
	else              : get_image(maskfile)
	
	# calculate total average using current alignment parameters
	tave,tvar = ave_var_series_g(data,kb)
	#drop_image(tave,"a1.spi")
	#drop_image(tvar,"a2.spi")
	a0 = tave.cmp("dot", tave, {"negative":0, "mask":mask})
	print  "Initial criterion=",a0
	
	# do the alignment
	# IMAGES ARE SQUARES!
	#  center is in SPIDER convention
	cnx = int(nx/2)+1
	cny = cnx

	mode = "F"
	#precalculate rings
	numr=Numrinit(first_ring,last_ring,rstep,mode)
	wr=ringwe(numr,mode)
	
	for iter in xrange(max_iter):
		again = False
		refim = tave
		refim_p = prepg(refim,kb)
		cimage=Util.Polar2Dmi(refim_p,cnx,cny,numr,mode,kb)
		Util.Frngs(cimage, numr)
		Applyws(cimage, numr, wr)
		a2 = a0
		for im in xrange(nima):
			# subtract current image from the average
			alpha = data[im].get_attr('alpha')
			sx =  data[im].get_attr('sx')
			sy =  data[im].get_attr('sy')
			mirror =  data[im].get_attr('mirror')
			temp = rtshgkb(data[im], alpha, sx, sy, kb)
			if  mirror: temp.process_inplace("mirror",{"axis":'x'})

			# align current image to the reference minus this image
			[angt,sxst,syst,mirrort,peakt]=ormy3(temp,cimage,xrng,yrng,step,mode,numr,cnx,cny,"gridding")			

			# combine parameters and set them to the header
			[angn,sxn,syn,mn]=combine_params2(alpha, sx, sy, mirror, angt, sxst, syst, mirrort)

			# apply params to the image
			temp = rtshgkb(data[im], angn, sxn, syn, kb)
			if  mn: temp.process_inplace("mirror",{"axis":'x'})

			# calculate the criterion
			a1 = temp.cmp("dot", tave, {"negative":0, "mask":mask})
			print a1, angn, sxn, syn, mn
			if a1 > a0:
				# replace the average by the improved average and set the new parameters to the image, otherwise, do nothing
				tave = temp.copy()
				data[im].set_attr_dict({'alpha':angn, 'sx':sxn, 'sy':syn,'mirror': mn})
				if a1>a2 : a2 = a1
				again = True
		a0 = a2
		if again :
			# calculate total average using current alignment parameters
			print " ITERATION #",iter,"  criterion = ",a0
			av1,av2 = ave_oe_series_g(data,kb)
			frsc = fsc(av1,av2,1.0,os.path.join(outdir,"drm%03d"%iter))
			if center:
				#center
				av1 = (av1*(int(nima/2)+(nima%2)) + av2*int(nima/2))/nima
				s = av1.phase_cog()
				tave = fshift(av1, -s[0], -s[1])
				for im in xrange(nima):
					alpha = data[im].get_attr('alpha')
					sx =  data[im].get_attr('sx')
					sy =  data[im].get_attr('sy')
					alphan,sxn,syn,scale = compose_transform2(alpha, sx, sy, 1.0, 0.0, -s[0], -s[1], 1.0)
					data[im].set_attr_dict({'alpha':alphan, 'sx':sxn, 'sy':syn})
			else:
				tave = (av1*(int(nima/2)+(nima%2)) + av2*int(nima/2))/nima
			# write the current average
			# fil_table=filt_from_fsc_bwt(frsc, low=0.5)
			# tave=filt_table(tave,fil_table)
			drop_image(tave,os.path.join(outdir,"aqm%03d.spi"%iter),"s")
		else:
			break

	for im in xrange(nima):
		temp.read_image(stack,im)
		alpha = data[im].get_attr('alpha')
		sx =  data[im].get_attr('sx')
		sy =  data[im].get_attr('sy')
		mn =  data[im].get_attr('mirror')
		temp.set_attr_dict({'alpha':alpha, 'sx':sx, 'sy':sy,'mirror': mn})
		temp.write_image(stack,im)
	#for im in xrange(nima): data[im].write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)
	return	
"""

def eqprojGSq(args, data):
	from projection import prgs
	from utilities import info

	proj = prgs(data[0], data[1], args)
	v = -proj.cmp("SqEuclidean", data[2], {"mask":data[3]})
	return v
	
def eqprojGccc(args, data):
	from projection import prgs
	from utilities import info

	proj = prgs(data[0], data[1], args)
	v = proj.cmp("ccc", data[2], {"mask":data[3], "negative":0})
	return v


def eqprojL(args, data):
	from projection import project

	proj = project(data[0], args, data[1])
	v = -proj.cmp("SqEuclidean", data[2], {"mask":data[3]})
	return v

def ali3d_e_G(stack, ref_vol, maskfile=None, radius=-1, dtheta=2, crit="SqEuclidean", opti_algorithm="amoeba"):
	"""
	An experimental version of ali3d_e, comparing SqEuclidean (or ccc) in the real space		
	"""
	from projection import prep_vol
	from utilities  import amoeba, model_circle, get_params_proj, set_params_proj
	from utilities  import get_image
	from math       import pi
	from sys        import exit
	from time       import time

	nima = EMUtil.get_image_count(stack)

	vol = EMData()
	vol.read_image(ref_vol)
	nx = vol.get_xsize()
	if radius <= 0:  radius = nx/2-1
	
	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  mask3D = get_image(maskfile)
		else: mask3D = maskfile
	else:
		mask3D = model_circle(radius, nx, nx, nx)
	mask2D = model_circle(radius, nx, nx)

	dataim = []
	for im in xrange(nima):
		ima = EMData()
		ima.read_image(stack, im)
		dataim.append(ima)
	
	data = [None]*4
	data[3] = mask2D

	total_time = 0.0
	total_time2 = 0.0

	Util.mul_img(vol, mask3D)
	volft, kb = prep_vol(vol)
	data[0] = volft
	data[1] = kb

	for imn in xrange(nima):
		begin_time = time()
	
		data[2] = dataim[imn]
		phi, theta, psi, s2x, s2y = get_params_proj(dataim[imn])
		atparams = [phi, theta, psi, s2x, s2y]
		#  change signs of shifts for projections
		atparams[3] *= -1
		atparams[4] *= -1
		weight_phi = max(dtheta, dtheta*abs((atparams[1]-90.0)/180.0*pi))
		
		if opti_algorithm == "amoeba":
			# For downhill simplex method			
			if crit == "SqEuclidean":
				optm_params =  amoeba(atparams, [weight_phi, dtheta, weight_phi, 1.0, 1.0], eqprojGSq, 1.e-5, 1.e-5, 500, data)
			elif crit == "ccc":
				optm_params =  amoeba(atparams, [weight_phi, dtheta, weight_phi, 1.0, 1.0], eqprojGccc, 1.e-5, 1.e-5, 500, data)
			else: 
				print "Unknown criterion!"
			optm_params[0][3] *= -1
			optm_params[0][4] *= -1
			set_params_proj(dataim[imn], optm_params[0])
		elif opti_algorithm == "LBFGSB":
			# For LBFGSB method
			optm_params = Util.twoD_to_3D_ali(volft, kb, dataim[imn], mask2D, atparams[0], atparams[1], atparams[2], atparams[3], atparams[4])
			optm_params[3] *= -1
			optm_params[4] *= -1
			set_params_proj(dataim[imn], [optm_params[0], optm_params[1], optm_params[2], optm_params[3], optm_params[4]])
		else:
			print "Unknown optimization method!"
			exit()
		used_time = time()-begin_time
		total_time += used_time
		total_time2 += used_time**2
	for im in xrange(nima):
		dataim[im].write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)
	
	return total_time, total_time2

def eqprojG2(args, data):
	from fundamentals import fft
	from sys import exit
	
	R = Transform({"type":"spider", "phi":args[0], "theta":args[1], "psi":args[2], "tx":0.0, "ty":0.0, "tz":0.0, "mirror":0, "scale":1.0})
	temp = data[0].extract_plane(R, data[1])
	temp.fft_shuffle()
	temp.center_origin_fft()

	if args[3]!=0.0 or args[4]!=0.0:
		params = {"filter_type":Processor.fourier_filter_types.SHIFT, "x_shift":args[3], "y_shift":args[4], "z_shift":0.0}
		temp = Processor.EMFourierFilter(temp, params)
	
	v = -temp.cmp("SqEuclidean", data[4])
	return v


def ali3d_e_G2(stack, ref_vol, maskfile=None, radius=-1, dtheta=2, opti_algorithm="amoeba"):
	"""
	An experimental version of ali3d_e, comparing SqEuclidean in the Fourier space		
	"""
	from projection import prep_vol
	from utilities  import amoeba, model_circle, get_params_proj, set_params_proj
	from utilities  import get_image
	from math       import pi
	from sys 	import exit
	from time       import time

	nima = EMUtil.get_image_count(stack)

	vol = EMData()
	vol.read_image(ref_vol)
	nx = vol.get_xsize()
	if radius <= 0:  radius = nx/2-1
	
	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  mask3D = get_image(maskfile)
		else: mask3D = maskfile
	else:
		mask3D = model_circle(radius, nx, nx, nx)
	mask2D = model_circle(radius, nx, nx)

	dataim = []
	for im in xrange(nima):
		ima = EMData()
		ima.read_image(stack, im)
		dataim.append(ima)

	data = [None]*5
	data[3] = mask2D
	
	total_time = 0.0
	total_time2 = 0.0

	Util.mul_img(vol, mask3D)
	volft, kb  = prep_vol(vol)
	data[0] = volft
	data[1] = kb
	for imn in xrange(nima):
		begin_time = time()
	
		data[2] = dataim[imn].copy()
		refi = dataim[imn].norm_pad(False, 2)
		refi.do_fft_inplace()
		data[4] = refi.copy()
		
		phi, theta, psi, sx, sy = get_params_proj(dataim[imn])
		atparams = [phi, theta, psi, sx, sy]
		#  change signs of shifts for projections
		atparams[3] *= -1
		atparams[4] *= -1
		weight_phi = max(dtheta, dtheta*abs((atparams[1]-90.0)/180.0*pi))			
		
		if opti_algorithm == "amoeba":
			# For downhill simplex method 
			optm_params = amoeba(atparams, [weight_phi, dtheta, weight_phi, 1.0, 1.0], eqprojG2, 1.e-5, 1.e-5, 500, data)
			optm_params[0][3] *= -1
			optm_params[0][4] *= -1
			set_params_proj(dataim[imn], optm_params[0])
		elif opti_algorithm == "LBFGSB":
			# For LBFGSB method
			optm_params = Util.twoD_to_3D_ali(volft, kb, dataim[imn], mask2D, atparams[0], atparams[1], atparams[2], atparams[3], atparams[4])
			optm_params[3] *= -1
			optm_params[4] *= -1
			set_params_proj(dataim[imn], [optm_params[0], optm_params[1], optm_params[2], optm_params[3], optm_params[4]])
		else:
			print "Unknown optimization algorithm!"
			exit()
		used_time = time()-begin_time
		total_time += used_time
		total_time2 += used_time**2

	for im in xrange(nima):
		dataim[im].write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)
	
	return total_time, total_time2
			
def eqprojG3(args, data):
	from utilities import amoeba
	from fundamentals import ccf
	from applications import twoD_fine_search
	
	R = Transform({"type":"spider", "phi":args[0], "theta":args[1], "psi":args[2], "tx":0.0, "ty":0.0, "tz":0.0, "mirror":0, "scale":1.0})
	temp = data[0].extract_plane(R, data[1])
	M    = temp.get_ysize()	
	"""
	temp.fft_shuffle()
	temp = temp.Four_ds(M/2, M/2, 1, False)
	temp.center_origin_fft()
	temp = temp.FourInterpol_i(M, M, 1, False)
	"""
	temp = temp.Four_shuf_ds_cen_us(M, M, 1, False)

	nx = M/2
	sx = (nx-data[7][0]*2)/2.0
	sy = (nx-data[7][1]*2)/2.0

	data2 = []
	data2.append(ccf(temp, data[5]))
	data2.append(data[6])
	#  search for shift
	ps = amoeba([sx, sy], [0.05, 0.05], twoD_fine_search, 1.e-5, 1.e-5, 500, data2)
	
	s2x = (nx-ps[0][0]*2)/2
	s2y = (nx-ps[0][1]*2)/2
	
	params2 = {"filter_type":Processor.fourier_filter_types.SHIFT, "x_shift":s2x*2, "y_shift":s2y*2, "z_shift":0.0}
	temp2 = Processor.EMFourierFilter(temp, params2)
	v = -temp2.cmp("SqEuclidean", data[4])
	
	return v, [s2x, s2y]

def prepij(image):
	M = image.get_ysize()
	npad = 2
	N = M*npad
	K = 6
	alpha = 1.75
	r = M/2
	v = K/2.0/N
	kb = Util.KaiserBessel(alpha, K, r, v, N)
	params = {"filter_type": Processor.fourier_filter_types.KAISER_SINH_INVERSE, "alpha":alpha, "K":K, "r":r, "v":v, "N":N}
	o = image.FourInterpol(2*M, 2*M, 1, 0)
	q = Processor.EMFourierFilter(o, params)
	return  o, q, kb


def ali3d_e_G3(stack, ref_vol, maskfile=None, radius=-1, dtheta=2):
	"""
	An experimental version of ali3d_e, using ccf in the Fourier space		
	"""
	from projection import prep_vol
	from utilities  import amoeba_multi_level, model_circle, get_params_proj, set_params_proj
	from utilities  import get_image
	from math       import pi
	from time       import time

	nima = EMUtil.get_image_count(stack)

	vol = EMData()
	vol.read_image(ref_vol)
	nx = vol.get_xsize()
	if radius <= 0:  radius = nx/2-1
	
	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  mask3D = get_image(maskfile)
		else: mask3D = maskfile
	else:
		mask3D = model_circle(radius, nx, nx, nx)
	mask2D = model_circle(radius, nx, nx)

	dataim = []
	for im in xrange(nima):
		ima = EMData()
		ima.read_image(stack, im)
		dataim.append(ima)
	
	data = [None]*8
	data[3] = mask2D

	total_time = 0.0
	total_time2 = 0.0

	Util.mul_img(vol, mask3D)
	volft, kb  = prep_vol(vol)
	data[0] = volft
	data[1] = kb

	for imn in xrange(nima):
		begin_time = time()
	
		data[2] = dataim[imn]
		refi = dataim[imn].copy()
		oo, qq, kb2 = prepij(refi)
		data[4] = oo
		data[5] = qq
		data[6] = kb2 
		
		phi, theta, psi, sx, sy = get_params_proj(dataim[imn])
		atparams = [phi, theta, psi]
		data[7] = [-sx, -sy]
		weight_phi = max(dtheta, dtheta*abs((atparams[1]-90.0)/180.0*pi))			
		
		# For downhill simplex method 
		optm_params = amoeba_multi_level(atparams, [weight_phi, dtheta, weight_phi], eqprojG3, 1.e-5, 1.e-5, 500, data)
		set_params_proj(dataim[imn], [optm_params[0][0], optm_params[0][1], optm_params[0][2], -optm_params[3][0], -optm_params[3][1]])
		
		used_time = time()-begin_time
		total_time += used_time
		total_time2 += used_time**2
	for im in xrange(nima):
		dataim[im].write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)
			
	return total_time, total_time2


def ali_G3(data, atparams, delta):
	from utilities import amoeba_multi_level
	from development import prepij
	oo, qq, kb2 = prepij(data[2])
	data.insert(4, oo)
	data.insert(5, qq)
	data.insert(6, kb2)
	data.insert(7, [-atparams[3], -atparams[4]])

	weight_phi = max(delta, delta*abs((atparams[1]-90.0)/180.0*pi))	
	optm_params =  amoeba_multi_level(atparams, [weight_phi, delta, weight_phi], eqprojG3, 1.e-5,1.e-5,500,data)

	del data[4]
	del data[4]
	del data[4]
	del data[4]

	return optm_params


def ali3d_e_G4(stack, ref_vol, maskfile=None, radius=-1, dtheta=2):
	"""
	An experimental version of ali3d_e, using ccf in the real space		
	"""
	from projection   import prep_vol
	from utilities    import amoeba_multi_level, model_circle, get_params_proj, set_params_proj
	from utilities    import get_image
	from math         import pi
	from time         import time

	nima = EMUtil.get_image_count(stack)

	vol = EMData()
	vol.read_image(ref_vol)
	nx = vol.get_xsize()
	if radius <= 0:  radius = nx/2-1
	
	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  mask3D = get_image(maskfile)
		else: mask3D = maskfile
	else:
		mask3D = model_circle(radius, nx, nx, nx)
	mask2D = model_circle(radius, nx, nx)

	dataim = []
	for im in xrange(nima):
		ima = EMData()
		ima.read_image(stack, im)
		ima *= mask2D
		ima.process_inplace("normalize.mask", {"mask":mask2D, "no_sigma":1})
		dataim.append(ima)
	
	npad = 2
	N = nx*npad
	K = 6
	alpha = 1.75
	r = nx/2
	v = K/2.0/N
	params = {"filter_type":Processor.fourier_filter_types.KAISER_SINH_INVERSE, "alpha":alpha, "K":K, "r":r, "v":v, "N":N}
	del alpha, r, v, K, N
	
	data = [None]*6
	data[3] = mask2D

	total_time = 0.0
	total_time2 = 0.0

	Util.mul_img(vol, mask3D)
	volft, kb  = prep_vol(vol)
	data[0] = volft
	data[1] = kb

	for imn in xrange(nima):
		begin_time = time()

		data[2] = dataim[imn]
		refi = dataim[imn].copy()
		refi = refi.FourInterpol(nx*2, nx*2, 0, True)
		data[4] = Processor.EMFourierFilter(refi, params)

		phi, theta, psi, sx, sy = get_params_proj(dataim[imn])
		atparams = [phi, theta, psi]
		data[5] = [-sx, -sy]
		weight_phi = max(dtheta, dtheta*abs((atparams[1]-90.0)/180.0*pi))

		# For downhill simplex method
		optm_params = amoeba_multi_level(atparams, [weight_phi, dtheta, weight_phi], eqprojG4, 1.e-5, 1.e-5, 500, data)
		set_params_proj(dataim[imn], [optm_params[0][0], optm_params[0][1], optm_params[0][2], -optm_params[3][0], -optm_params[3][1]])

		used_time = time()-begin_time
		total_time += used_time
		total_time2 += used_time**2
	for im in xrange(nima):
		dataim[im].write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)
		
	return total_time, total_time2


def eqprojG4(args, data):
	from utilities import peak_search, amoeba
	from fundamentals import fft, ccf, fpol
	from applications import twoD_fine_search

	volft 	= data[0]
	kb	= data[1]
	prj	= data[2]
	mask2D	= data[3]
	refi	= data[4]
	shift	= data[5]

	R = Transform({"type":"spider", "phi":args[0], "theta":args[1], "psi":args[2], "tx":0.0, "ty":0.0, "tz":0.0, "mirror":0, "scale":1.0})
	temp = volft.extract_plane(R, kb)
	temp.fft_shuffle()
	temp.center_origin_fft()
	temp.do_ift_inplace()
	M = temp.get_ysize()/2
	refprj = temp.window_center(M)

	refprj.process_inplace("normalize.mask", {"mask":mask2D, "no_sigma":1})
	refprj *= mask2D
	
	nx = refprj.get_ysize()
	sx = (nx-shift[0]*2)/2
	sy = (nx-shift[1]*2)/2
	
	proj2x = fpol(refprj, 2*M, 2*M, 0, False)
	product = ccf(proj2x, refi)

	data2 = [0]*2
	data2[0] = product
	data2[1] = kb
	ps = amoeba([sx, sy], [0.05, 0.05], twoD_fine_search, 1.e-4, 1.e-4, 500, data2)

	v = ps[1]
	s2x = nx/2-ps[0][0]
	s2y = nx/2-ps[0][1]

	return v, [s2x, s2y]

def ali3d_e_G_DUD(stack, ref_vol, maskfile=None, radius=-1, dtheta=2):
	"""
	An experimental version of ali3d_e, using DUD algorithm for optimization
	"""
	from projection   import prep_vol, prgs
	from utilities    import get_params_proj, set_params_proj
	from utilities    import model_circle, get_image
	from math         import pi
	from time         import time
	from numpy        import matrix, diag, concatenate
	from utilities    import get_image_data

	nima = EMUtil.get_image_count(stack)

	vol = EMData()
	vol.read_image(ref_vol)
	nx = vol.get_xsize()
	if radius <= 0:  radius = nx/2-1
	
	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  mask3D = get_image(maskfile)
		else: mask3D = maskfile
	else:
		mask3D = model_circle(radius, nx, nx, nx)
	mask2D = model_circle(radius, nx, nx)

	dataim = []
	for im in xrange(nima):
		ima = EMData()
		ima.read_image(stack, im)
		dataim.append(ima)
	
	volft, kb  = prep_vol(vol)

	total_time = 0.0
	total_time2 = 0.0

	for im in xrange(nima):
		begin_time = time()

		refim = dataim[im].copy()

		phi, theta, psi, s2x, s2y = get_params_proj(refim)

		s2x = -s2x
		s2y = -s2y

		P = 5
		max_P = 50
		converge = False
		beta = []
		delta = 0.1
		beta.append([phi, theta, psi, s2x, s2y])
		beta.append([phi+delta, theta, psi, s2x, s2y])
		beta.append([phi, theta+delta, psi, s2x, s2y])
		beta.append([phi, theta, psi+delta, s2x, s2y])
		beta.append([phi, theta, psi, s2x+delta, s2y])
		beta.append([phi, theta, psi, s2x, s2y+delta])

		Y = matrix(get_image_data(refim)).flatten(1)
		omega = []
		proj_0 = prgs(volft, kb, beta[0])
		for p in xrange(1, P+1):
			proj_p = prgs(volft, kb, beta[p])
			diffim = Util.subn_img(proj_p, proj_0)
			P2 = matrix(get_image_data(diffim)).flatten(1)
			if p == 1: F = P2
			else: F = concatenate((F, P2))
			omega.append(proj_p.cmp("SqEuclidean", proj_0))
			
					
		proj_p = prgs(volft, kb, beta[P])
		Fs = matrix(get_image_data(proj_p)).flatten(1)

		A = [[delta, 0, 0, 0, 0], [0, delta, 0, 0, 0], [0, 0, delta, 0, 0], [0, 0, 0, delta, 0], [0, 0, 0, 0, delta]]
		while converge == False and P < max_P:
			A_matrix = matrix(A)
			omega_matrix_inv = matrix(diag(omega))**(-1)
			G = F.T*omega_matrix_inv*A_matrix*(A_matrix.T*omega_matrix_inv*A_matrix)**(-1)
					
			H = (G.T*G)**(-1)*G.T*(Y.T-Fs.T)

			dphi = float(H[0][0])
			dtheta = float(H[1][0])
			dpsi = float(H[2][0])
			ds2x = float(H[3][0])
			ds2y = float(H[4][0])
			d = dphi**2+dtheta**2+dpsi**2+ds2x**2+ds2y**2 
			if d < 1e-6: converge = True
			if converge == False:
				A.append([dphi, dtheta, dpsi, ds2x, ds2y])
		 		beta.append([beta[P][0]+dphi, beta[P][1]+dtheta, beta[P][2]+dpsi, beta[P][3]+ds2x, beta[P][4]+ds2y])
				P += 1
				proj_p = prgs(volft, kb, beta[P])
				diffim = Util.subn_img(proj_p, proj_0)
				P2 = matrix(get_image_data(diffim)).flatten(1)
				F = concatenate((F, P2))
				Fs = matrix(get_image_data(proj_p)).flatten(1)
				omega.append(proj_p.cmp("SqEuclidean", proj_0))

		set_params_proj(dataim[im], [beta[P][0], beta[P][1], beta[P][2], -beta[P][3], -beta[P][4]])

		used_time = time()-begin_time
		total_time += used_time
		total_time2 += used_time**2
		
	for im in xrange(nima):
		dataim[im].write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)
		
	return total_time, total_time2


			
def ali3d_e_L(stack, ref_vol, outdir, maskfile, radius=-1, snr=1.0, dtheta=2, max_it=10, symmetry="c1", CTF = False):
	"""
	A temporary replacement of ali3d_e		
	"""
	from alignment	import Numrinit, ringwe, Applyws
	from filter     import filt_ctf, filt_params, filt_table, filt_from_fsc
	from fundamentals   import fshift, fft
	from projection import prep_vol
	from utilities  import amoeba, model_circle, get_arb_params, set_arb_params, drop_spider_doc
	from utilities  import drop_image, amoeba_multi_level
	from math       import pi,sin
	from string     import replace
	from reconstruction import recons3d_4nn, recons3d_4nn_ctf
	from statistics import fsc
	import os 
	import sys
	from time import time

	nima = EMUtil.get_image_count(stack)

	vol = EMData()
	vol.read_image(ref_vol)
	nx  = vol.get_xsize()
	if radius <= 0:  radius = nx//2-1
	if os.path.exists(outdir):  os.system('rm -rf '+outdir)
	os.mkdir(outdir)
	vol.write_image(os.path.join(outdir,"ref_volf00.hdf"))
	
	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask3D=get_image(maskfile)
		else: mask3D = maskfile
	else:
		mask3D = model_circle(radius, nx, nx, nx)
	mask2D = model_circle(radius, nx, nx)
	dataim = []
	parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
	for im in xrange(nima):
		ima = EMData()
		ima.read_image(stack, im)
		if(CTF):
			ctf_params = get_arb_params(ima, parnames)
			if(ctf_params[6] == 0):
				from filter import filt_ctf
				if(im==0): print  " APPLYING CTF"
				ima = filt_ctf(ima, ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
				set_arb_params(ima, ctf_params, parnames)  #I am not sure whether this is needed
		dataim.append(ima)
	

	par_str=["phi", "theta", "psi", "s2x", "s2y"]
	total_time = 0.0
	total_time2 = 0.0
	for iteration in xrange(max_it):
		Util.mul_img(vol, mask3D)
		data = []
		data.insert(0, vol)
		data.insert(1, radius)
		data.insert(3, mask2D)
		#new_params = []
		#step = 0
		for imn in xrange(nima):
			begin_time = time()
			#print iteration, imn
			data.insert(2, dataim[imn])
			
			atparams = get_arb_params(dataim[imn], par_str)
			#  change signs of shifts for projections
			#atparams[3] *= -1
			#atparams[4] *= -1
			#initial  = eqproj(atparams, data)  # this is if we need initial discrepancy
			
			#stheta   = abs(sin(atparams[1]*pi/180.0))
			#if stheta > 0.0: weight_phi=min(180.0,dtheta/stheta/2.0)
			#else           : weight_phi=180.0
			weight_phi = max(dtheta, dtheta*abs((atparams[1]-90.0)/180.0*pi))			
			optm_params =  amoeba(atparams, [weight_phi,dtheta,weight_phi,1.0,1.0], eqprojL, 1.e-5,1.e-5,500,data)
			
			optm_params[0].append(imn)
			#optm_params[0][3] *= -1
			#optm_params[0][4] *= -1
			#new_params.append(optm_params[0])
			del  data[2]
			set_arb_params(dataim[imn], [optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4]], par_str)
			end_time = time()
			total_time += (end_time-begin_time) 
			total_time2 += (end_time-begin_time)**2
			#step += optm_params[2]
			#print "Projection", imn, "  ", optm_params[2], "steps"
		#print "Average step = ", step/float(nima)
		"""
		#  3D stuff
		list_p = range(0,nima,2)
		if(CTF): vol1 = recons3d_4nn_ctf(stack, list_p, snr, 1, symmetry)
		else:	 vol1 = recons3d_4nn(stack, list_p, symmetry)

		list_p = range(1,nima,2)
		if(CTF): vol2 = recons3d_4nn_ctf(stack, list_p, snr, 1, symmetry)
		else:	 vol2 = recons3d_4nn(stack, list_p, symmetry)

		[mean, sigma, xmin, xmax ] =  Util.infomask(vol1, mask3D, False)
		vol1 -= mean
		Util.mul_img(vol1, mask3D)
		[mean, sigma, xmin, xmax ] =  Util.infomask(vol2, mask3D, False)
		vol2 -= mean
		Util.mul_img(vol2, mask3D)
		
		filename = replace("fsc%3d"%(iteration+1)," ","0")
		res = fsc(vol1, vol2, 1.0, os.path.join(outdir, filename))
		del vol1
		del vol2

		# calculate new and improved 3D
		list_p = range(nima)
		if(CTF): vol = recons3d_4nn_ctf(stack, list_p, snr, 1, symmetry)
		else:	 vol = recons3d_4nn(stack, list_p, symmetry)
		# store the reference volume
		filename =  replace("vol%3d.spi"%(iteration+1)," ","0")
		drop_image(vol,os.path.join(outdir, filename), "s")
		# lk = 0
		# while(fscc[1][lk] >0.9 and fscc[0][lk]<0.25):	lk+=1
		# fl = fscc[0][lk]
		# fh = min(fl+0.1,0.49)
		# vol = filt_btwl(vol, fl, fh)
		cs   = vol.phase_cog()
		vol  = fshift(vol, -cs[0], -cs[1] -cs[2])
		filename = replace("volf%3d.spi"%(iteration+1)," ","0")
		drop_image(vol,os.path.join(outdir, filename), "s")
		"""
		
		for im in xrange(nima):
			dataim[im].write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)
	return total_time, total_time2


def max_3D_pixel_error(t1, t2, r):
	"""
	  Compute maximum pixel error between two projection directions
	  assuming object has radius r, t1 is the projection transformation
	  of the first projection and t2 of the second one, respectively:
		t = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
		t.set_trans(Vec2f(-tx, -ty))
	  Note the function is symmetric in t1, t2.
	"""
	from math import sin, cos, pi, sqrt
	t3 = t2*t1.inverse()
	tm = t3.get_matrix()

	ddmax = 0.0
	for i in xrange(1, r+1):
		for ang in xrange(int(2*pi*i+0.5)):
			sx = i*cos(ang)
			sy = i*sin(ang)
			dx = (tm[0]-1)*sx+tm[1]*sy+tm[3]
			dy = tm[4]*sx+(tm[5]-1)*sy+tm[7]
			dz = tm[8]*sx+tm[9]*sy+tm[8]
			dd = dx*dx+dy*dy+dz*dz
			if dd > ddmax: ddmax=dd
	return sqrt(ddmax)


def max_3D_pixel_error1(phi1, theta1, psi1, sx1, sy1, phi2, theta2, psi2, sx2, sy2, r):
	from numpy import mat, linalg
	from math import cos, sin, sqrt, pi

	phi1 = phi1/180.0*pi
	theta1 = theta1/180.0*pi
	psi1 = psi1/180.0*pi
	phi2 = phi2/180.0*pi
	theta2 = theta2/180.0*pi
	psi2 = psi2/180.0*pi	
	R1a = mat([[cos(phi1), sin(phi1), 0],[-sin(phi1), cos(phi1), 0],[0, 0, 1]])
	R1b = mat([[cos(theta1), 0, -sin(theta1)],[0, 1, 0],[sin(theta1), 0, cos(theta1)]])
	R1c = mat([[cos(psi1), sin(psi1), 0],[-sin(psi1), cos(psi1), 0],[0, 0, 1]])
	R1  = R1a*R1b*R1c	
	R1inv = linalg.inv(R1) 
	R2a = mat([[cos(phi2), sin(phi2), 0],[-sin(phi2), cos(phi2), 0],[0, 0, 1]])
	R2b = mat([[cos(theta2), 0, -sin(theta2)],[0, 1, 0],[sin(theta2), 0, cos(theta2)]])
	R2c = mat([[cos(psi2), sin(psi2), 0],[-sin(psi2), cos(psi2), 0],[0, 0, 1]])
	R2  = R2a*R2b*R2c
	s1  = mat([[sx1],[sy1],[0]])
	s2  = mat([[sx2],[sy2],[0]])
	dmax = 0
	for i in range(-r,r+1):
		for j in range(-r,r+1):
			for k in range(-r,r+1):
				dis = sqrt(i*i+j*j+k*k)
				if (r-dis >= 0) and (r-dis <= 1):
					X = mat([[i],[j],[k]])
					d = R1inv*(R2*X+s2-s1)-X
					dd = sqrt(d[0,0]*d[0,0]+d[1,0]*d[1,0]+d[2,0]*d[2,0])
					if dd > dmax: 
						dmax=dd
	return dmax

def max_3D_pixel_error2(phi1, theta1, psi1, sx1, sy1, phi2, theta2, psi2, sx2, sy2, r):
	from numpy import mat, linalg
	from math import cos, sin, sqrt, pi
	from utilities import even_angles

	phi1 = phi1/180.0*pi
	theta1 = theta1/180.0*pi
	psi1 = psi1/180.0*pi
	phi2 = phi2/180.0*pi
	theta2 = theta2/180.0*pi
	psi2 = psi2/180.0*pi	
	R1a = mat([[cos(phi1), sin(phi1), 0],[-sin(phi1), cos(phi1), 0],[0, 0, 1]])
	R1b = mat([[cos(theta1), 0, -sin(theta1)],[0, 1, 0],[sin(theta1), 0, cos(theta1)]])
	R1c = mat([[cos(psi1), sin(psi1), 0],[-sin(psi1), cos(psi1), 0],[0, 0, 1]])
	R1  = R1a*R1b*R1c	
	R1inv = linalg.inv(R1) 
	R2a = mat([[cos(phi2), sin(phi2), 0],[-sin(phi2), cos(phi2), 0],[0, 0, 1]])
	R2b = mat([[cos(theta2), 0, -sin(theta2)],[0, 1, 0],[sin(theta2), 0, cos(theta2)]])
	R2c = mat([[cos(psi2), sin(psi2), 0],[-sin(psi2), cos(psi2), 0],[0, 0, 1]])
	R2  = R2a*R2b*R2c
	s1  = mat([[sx1],[sy1],[0]])
	s2  = mat([[sx2],[sy2],[0]])
	dmax = 0
	angles = even_angles(0.5, 0.0, 180.0, 0.0, 359.9, "S", symmetry='c1', phiEqpsi='Minus')
	nangles = len(angles)
	for i in range(nangles):
		X = mat([[r*sin(angles[i][1]*pi/180.0)*cos(angles[i][0]*pi/180.0)],[r*sin(angles[i][1]*pi/180.0)*sin(angles[i][0]*pi/180.0)],[r*cos(angles[i][1]*pi/180.0)]])
		d = R1inv*(R2*X+s2-s1)-X
		dd = sqrt(d[0,0]*d[0,0]+d[1,0]*d[1,0]+d[2,0]*d[2,0])
		if dd > dmax: 	
			dmax=dd
	return dmax

"""	
def max_3D_pixel_error3(phi1, theta1, psi1, sx1, sy1, phi2, theta2, psi2, sx2, sy2, r):
	from numpy import mat, linalg, concatenate, eye
	from math import cos, sin, sqrt, pi
	from utilities import even_angles

	phi1 = phi1/180.0*pi
	theta1 = theta1/180.0*pi
	psi1 = psi1/180.0*pi
	phi2 = phi2/180.0*pi
	theta2 = theta2/180.0*pi
	psi2 = psi2/180.0*pi	
	R1a = mat([[cos(phi1), sin(phi1), 0],[-sin(phi1), cos(phi1), 0],[0, 0, 1]])
	R1b = mat([[cos(theta1), 0, -sin(theta1)],[0, 1, 0],[sin(theta1), 0, cos(theta1)]])
	R1c = mat([[cos(psi1), sin(psi1), 0],[-sin(psi1), cos(psi1), 0],[0, 0, 1]])
	R1  = R1a*R1b*R1c
	R1e = concatenate((concatenate((R1,[[sx1],[sy1],[0]]),1), [[0,0,0,1]]))
	R1inv = linalg.inv(R1e) 
	R2a = mat([[cos(phi2), sin(phi2), 0],[-sin(phi2), cos(phi2), 0],[0, 0, 1]])
	R2b = mat([[cos(theta2), 0, -sin(theta2)],[0, 1, 0],[sin(theta2), 0, cos(theta2)]])
	R2c = mat([[cos(psi2), sin(psi2), 0],[-sin(psi2), cos(psi2), 0],[0, 0, 1]])
	R2  = R2a*R2b*R2c
	R2e = concatenate((concatenate((R2,[[sx2],[sy2],[0]]),1), [[0,0,0,1]]))
	dmax = 0
	angles = even_angles(0.5, 0.0, 180.0, 0.0, 359.9, "S", symmetry='c1', phiEqpsi='Minus')
	nangles = len(angles)
	for i in range(nangles):
		X = mat([[r*sin(angles[i][1]*pi/180.0)*cos(angles[i][0]*pi/180.0)],[r*sin(angles[i][1]*pi/180.0)*sin(angles[i][0]*pi/180.0)],[r*cos(angles[i][1]*pi/180.0)], [1.0]])
		d = (R1inv*R2e-eye(4))*X
		dd = sqrt(d[0,0]*d[0,0]+d[1,0]*d[1,0]+d[2,0]*d[2,0])
		if dd > dmax: 	
			dmax=dd
	return dmax
"""

def max_3D_pixel_error3(phi1, theta1, psi1, sx1, sy1, phi2, theta2, psi2, sx2, sy2, r):
	from numpy import mat, linalg, concatenate, eye
	from math import cos, sin, sqrt, pi
	from utilities import even_angles

	PI180 = pi/180.0
	phi1 = phi1*PI180
	theta1 = theta1*PI180
	psi1 = psi1*PI180
	phi2 = phi2*PI180
	theta2 = theta2*PI180
	psi2 = psi2*PI180
	R1a = mat([[cos(phi1), sin(phi1), 0], [-sin(phi1), cos(phi1), 0], [0, 0, 1]])
	R1b = mat([[cos(theta1), 0, -sin(theta1)], [0, 1, 0], [sin(theta1), 0, cos(theta1)]])
	R1c = mat([[cos(psi1), sin(psi1), 0], [-sin(psi1), cos(psi1), 0], [0, 0, 1]])
	R1  = R1a*R1b*R1c
	R1e = concatenate((concatenate((R1, [[sx1], [sy1], [0]]), 1), [[0, 0, 0, 1]]))
	R1inv = linalg.inv(R1e) 
	R2a = mat([[cos(phi2), sin(phi2), 0], [-sin(phi2), cos(phi2), 0], [0, 0, 1]])
	R2b = mat([[cos(theta2), 0, -sin(theta2)], [0, 1, 0], [sin(theta2), 0, cos(theta2)]])
	R2c = mat([[cos(psi2), sin(psi2), 0], [-sin(psi2), cos(psi2), 0], [0, 0, 1]])
	R2  = R2a*R2b*R2c
	R2e = concatenate((concatenate((R2, [[sx2], [sy2], [0]]), 1), [[0, 0, 0, 1]]))
	L = R1inv*R2e-eye(4)
	dmax = 0
	angles = even_angles(1.0, 0.0, 180.0, 0.0, 359.9, "S", symmetry='c1', phiEqpsi='Minus')
	nangles = len(angles)
	for i in range(nangles):
		rad0 = angles[i][0]*PI180
		rad1 = angles[i][1]*PI180
		X = mat([[r*sin(rad1)*cos(rad0)], [r*sin(rad1)*sin(rad0)], [r*cos(rad1)], [1.0]])
		d = L*X
		dd = d[0,0]*d[0,0]+d[1,0]*d[1,0]+d[2,0]*d[2,0]
		if dd > dmax: 	
			dmax=dd
	return sqrt(dmax)
	

def ali_with_mask(vol, mdl, mask):
	vol_mask = vol*mask
	mdl_mask = mdl*mask

	vol_cs = vol_mask.phase_cog()
	mdl_cs = mdl_mask.phase_cog()

	vol_mask = fshift(vol, -vol_cs[0], -vol_cs[1], -vol_cs[2])
	mdl_mask = fshift(mdl, -mdl_cs[0], -mdl_cs[1], -mdl_cs[2])
	trans = ali_vol_3(vol_mask, mdl_mask, 2.0, 1.0)
	newvol = fshift(vol, -vol_cs[0], -vol_cs[1], -vol_cs[2])
	newvol = rot_shift3D(newvol, trans[0], trans[1], trans[2], trans[3], trans[4], trans[5], 1.0)
	newvol = fshift(vol,  mdl_cs[0],  mdl_cs[1],  mdl_cs[2])
	return newvol

#def ali_vol_rotation(vol, refv, radius=None):
def rot_3D_add_n(vol, refv, phi=0, theta=0, psi=0, radius=None):
	# rotation only
	from alignment    import ali_vol_func
	from utilities    import get_im, model_circle, model_blank, get_arb_params, set_arb_params, drop_image
	from utilities    import amoeba
	from fundamentals import rot_shift3D
	ang_scale=5.0
	ref = get_im(refv)
	nx = ref.get_xsize()
	ny = ref.get_ysize()
	nz = ref.get_zsize()
	if(radius != None):    mask = model_circle(radius, nx, ny, nz)
	else:                  mask = model_circle(float(min(nx, ny, nz)//2-2), nx, ny, nz)

	names_params = ["phi", "theta", "psi"]
	params = get_arb_params(ref, names_params)
	ref = rot_shift3D(ref, params[0], params[1], params[2])

	e = get_im(vol)
	params = get_arb_params(ref, names_params)
	print  params
	scale = [ang_scale, ang_scale, ang_scale]
	data=[e, ref, mask]
	params = amoeba(params, scale, ali_vol_func, 1.e-4, 1.e-4,500, data)
	print  params
	set_arb_params(e, params, names_params)
	e.write_image(stack, n,  EMUtil.ImageType.IMAGE_HDF, True)

def rot_3D_add_n(stack, outdir, phi=0, theta=0, psi=0, radius=None):
	
	from alignment import find_rotation
	from utilities import get_im, model_circle, model_blank, get_arb_params, set_arb_params, drop_image
	from fundamentals import rot_shift3D
	import sys
	import os

	N = EMUtil.get_image_count(stack)

	ref = get_im(stack, 0)
	nx = ref.get_xsize()
	ny = ref.get_ysize()
	nz = ref.get_zsize()
	if(radius != None):
		mask = model_circle(radius, nx, ny, nz)
	else:
		mask = model_circle(float(min(nx, ny, nz)//2-2), nx, ny, nz)

	names_params = ["phi", "theta", "psi", "s3x", "s3y", "s3z"]
	params = get_arb_params(ref, names_params)
	ref = rot_shift3D(ref, params[0], params[1], params[2], params[3], params[4], params[5])

	for n in xrange(1,N):
		e = get_im(stack, n)
		params = get_arb_params(ref, names_params)
		print  n,params
		params = find_rotation(e, ref, mask, params)
		e = rot_shift3D(e, params[0], params[1], params[2], params[3], params[4], params[5])
		Util.add_img(ref, e)
		print  n,params
		set_arb_params(e, params, names_params)
		e.write_image(stack, n,  EMUtil.ImageType.IMAGE_HDF, True)
	Criter = ref.cmp("dot", ref, {"negative":0, "mask":mask})
	

	angle_change = True
	cnt=0
	while(angle_change):
		print "Iteration: ",cnt,Criter
		nref = model_blank(nx, ny, nz)
		for n in xrange(N):
			e = get_im(stack, n)
			params = get_arb_params(ref, names_params)
			temp = rot_shift3D(e, params[0], params[1], params[2], params[3], params[4], params[5])
			Util.sub_img(ref, temp)
			Util.add_img(nref, temp)
			params = find_rotation(e, ref, mask, params)
			e = rot_shift3D(e, params[0], params[1], params[2], params[3], params[4], params[5])
			Util.add_img(ref, e)
			set_arb_params(e, params, names_params)
			e.write_image(stack, n,  EMUtil.ImageType.IMAGE_HDF, True)
		cnt += 1
		newt = nref.cmp("dot", nref, {"negative":0, "mask":mask})
		if(newt <= Criter):  angle_change = False
		Criter = newt
		ref = nref.copy()

	drop_image(ref,outdir+"/vols.spi")



def pufit(infile, Pixel_size):
	""" 
		Taken from morphology
		Fit power spectrum model from 
		Power spectrum of the structure 
		created by X-ray crystallography
	"""
	from utilities  	import get_image, amoeba, read_spider_doc, get_im
	from math               import log
	from fundamentals 	import rops
	#from numpy 		import exp
	#from scipy.optimize 	import fmin,leastsq, fmin
	#from optparse 		import OptionParser
	#from scipy 		import arange
	#from numpy 		import log
	#from sparx 		import TFList
	fpw2         = rops(get_im(infile))
	nr           = fpw2.get_xsize()
	fpw = []
	frq = []
	for i in xrange(1,nr):
		frq.append(i/float(2*nr)/Pixel_size)
		fpw.append(fpw2.get_value_at(i))
	scale = fpw[0]
	norm = fpw[1]
	nr = len(fpw)
	for i in xrange(nr):	fpw[i] /= norm
	flo = []
	for i in xrange(nr):	flo.append(log(fpw[i]))
	data = []
	data.insert(0,flo)
	data.insert(1,frq)
	params = [-10.0, 10.0, 0.05]
	ooo = amoeba(params, [0.1, 1.0, 0.001], residuals_simplex, 1.e-4, 1.e-4, 500, data)
	print ooo
	for i in xrange(1,nr):  print  frq[i], flo[i], ooo[0][0] + ooo[0][1]/(data[1][i]/ooo[0][2]+1.0)**2
	return  ooo[0]
	'''
	# fit Gaussian like peaks 
	#  Fitting the peak around 1/10 Angstrom
	if 2.*Pixel_size+1. < 10.:
		xpeak_1   = 1./(5.) # 7
		xpeak_2   = 1./(11.)  # 16
		ix_1      = int(xpeak_1/step)
		ix_2      = int(xpeak_2/step)
		freq_carb = .1
		sig_carb  = sgm1**2
		p1        = [y_true[ix_1], freq_carb, sig_carb]
		print '%10.3g %10.3g'%(ix_1,ix_2)
		y1        = y_true[ix_2:ix_1]
		x1        = x[ix_2:ix_1]
		print plsq[0]
		plsq     = leastsq(residuals_lsq_peak,p1,(y1,x1,c0))
		d1,d2,d3 = plsq[0]
		c1,c2,c3 = c0
		for i in xrange(len(pw2)):	pw2[i]=exp(c1)*exp(c2/(x[i]+c3)**2)+d1*exp(-(x[i]-d2)**2/d3)
       # fitting the peak around 1/5 Angstrom
		if (2.*Pixel_size + 1. < 5.):
			xpeak_1   = 1./(2.*Pixel_size+1.)
			xpeak_2   = 1./(6.)
			ix_1      = int(xpeak_1/step)
			ix_2      = int(xpeak_2/step)
			freq_carb = .2
			sig_carb  = sgm2**2
			p1        = [y_true[ix_1],freq_carb,sig_carb]
			y1        = y_true[ix_2:ix_1]
			x1        = x[ix_2:ix_1]
			plsq      = leastsq(residuals_lsq_peak,p1,(y1,x1,c0))
			p1        = plsq[0]
			plsq      = leastsq(residuals_lsq_peak,p1,(y1,x1,c0))
			h1,h2,h3  = plsq[0]
			c1,c2,c3  = c0
			for i in xrange(len(pw2)):	pw2[i] = exp(c1)*exp(c2/(x[i] + c3)**2)+d1*exp(-(x[i]-d2)**2/d3) + h1*exp(-(x[i]-h2)**2/h3)
	else: # Fitting those power spectrum whose nyquist frequency less than 1/10 Angstrom 
		pw2 = exp(c1)*exp(c2/(x+c3)**2)
	out = TFList(nrows,3)
	for i in xrange(1,nrows+1):
		out.SetVal(i,2,y_true[i-1]/y_scale)
		out.SetVal(i,1,x[i-1])
		out.SetVal(i,2,pw2[i-1]/y_scale)	
	out.write(outfile)
	'''


def extract_signal_from_1dpw2_cl1(indir, outdir, Pixel_size = 1, polynomial_rankB = 2, cut_off = 0):
	"""
		Extract signal from 1D rotationally averaged 
		power spectra saved in a directory 
	"""
	import os
	from morphorlogy import residual_1dpw2
	flist = os.listdir(outdir)
	if os.path.exists(indir)  is False: os.mkdir(outdir)
	for i, single_file in xrange(len(flist)):
		(filename, filextension) = os.path.splitext(single_file)
		input_file    = os.path.join(indir, single_file)
		pw = read_spider_doc(input_file)
		out_file     = os.path.join(outdir, "residual_of_" + filename + ".txt")
		res = residual_1dpw2(list_1dpw2, polynomial_rankB, Pixel_size, cut_off)
		drop_spider_doc(out_file, res)



"##################   WHAT FOLLOWS ARE VARIOUS ALIGNMENT PROCEDURES     #####################"


def oned_search_func(args, data):
	#print  " AMOEBA i",args
	v = data[0].get_pixel_conv7(args[0]*2, 0, 0, data[1])
	#print  " AMOEBA o",args,v
	return v
	
def apmq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny):
	"""Determine shift and rotation between image and many reference images (crefim, weights have to be applied)
		quadratic interpolation
	"""
	from math import pi, cos, sin
	#from utilities import ttime#,info
	#print "APMQ"
	#print  ttime()
	#cimage=Util.Polar2Dm(image, cnx, cny, numr, mode)
	peak = -1.0E23
	ky = int(2*yrng/step+0.5)//2
	kx = int(2*xrng/step+0.5)//2
	for i in xrange(-ky, ky+1):
		iy = i*step
		for  j in xrange(-kx, kx+1):
			ix = j*step
			cimage=Util.Polar2Dm(image, cnx+ix, cny+iy, numr, mode)
			#for l in xrange(20000):   cimage=Util.Polar2Dm(image, cnx+ix, cny+iy, numr, mode)
			#print  "Polar2Dm ",ttime()
			Util.Frngs(cimage, numr)
			#for l in xrange(20000):   retvals = Util.Crosrng_ms(cimage, cimage, numr)
			#print  "Crosrng_ms ",ttime()
			#from sys import exit
			#exit()
			#imref = retvals["qn"]
			#for l in xrange(6000):   Util.Frngs(cimage, numr)
			#print  "Frngs    ",ttime()
			#  compare with all reference images
			for iref in xrange(len(crefim)):
				retvals = Util.Crosrng_ms(crefim[iref], cimage, numr)
				#print  retvals
				#from sys import exit
				#exit()
				#for l in xrange(6000):  retvals = Util.Crosrng_ms(crefim[iref], cimage, numr)
				#print  "Crosrng  ",ttime()
				qn = retvals["qn"]
				qm = retvals["qm"]
				if(qn >=peak or qm >=peak):
					sx = -ix
					sy = -iy
					nref = iref
					if (qn >= qm):
						ang = ang_n(retvals["tot"], mode, numr[len(numr)-1])
						peak=qn
						mirror=0
					else:
						ang = ang_n(retvals["tmt"], mode, numr[len(numr)-1])
						peak=qm
						mirror=1
	co =  cos(ang*pi/180.0)
	so = -sin(ang*pi/180.0)
	sxs = sx*co - sy*so
	sys = sx*so + sy*co
	return  ang, sxs, sys, mirror, peak, nref

def diff_projdir(n1, phi2, theta2):
	from math import cos, sin, pi
	qv = pi/180.
	#n1 = [ sin(theta1*qv)*cos(phi1) ,  sin(theta1*qv)*sin(phi1*qv) , cos(theta1*qv)]
	n2 = [ sin(theta2*qv)*cos(phi2*qv) ,  sin(theta2*qv)*sin(phi2*qv) , cos(theta2*qv)]
	#return acos(abs(n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2]))*180.0/pi
	return abs(n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2])
	
def apmq_local(image, crefim, xrng, yrng, step, ant, mode, numr, cnx, cny):
	"""Determine shift and rotation between image and many reference images (crefim, weights have to be applied)
		quadratic interpolation
	"""
	from math import pi, cos, sin
	from alignment import diff_projdir
	#from utilities import ttime#,info
	#print "APMQ"
	#print  ttime()
	#cimage=Util.Polar2Dm(image, cnx, cny, numr, mode)
	phi,theta,psi,s2x,s2y = get_params_proj( image )
	peak = -1.0E23
	ky = int(2*yrng/step+0.5)//2
	kx = int(2*xrng/step+0.5)//2
	for i in xrange(-ky, ky+1):
		iy = i*step
		for  j in xrange(-kx, kx+1):
			ix = j*step
			cimage = Util.Polar2Dm(image, cnx+ix, cny+iy, numr, mode)
			#for l in xrange(20000):   cimage=Util.Polar2Dm(image, cnx+ix, cny+iy, numr, mode)
			#print  "Polar2Dm ",ttime()
			Util.Frngs(cimage, numr)
			#for l in xrange(20000):   retvals = Util.Crosrng_ms(cimage, cimage, numr)
			#print  "Crosrng_ms ",ttime()
			#from sys import exit
			#exit()
			#imref = retvals["qn"]
			#for l in xrange(6000):   Util.Frngs(cimage, numr)
			#print  "Frngs    ",ttime()
			#  compare with all reference images
			for iref in xrange(len(crefim)):
				n1 = []
				n1.append(crefim[iref].get_attr('n1'))
				n1.append(crefim[iref].get_attr('n2'))
				n1.append(crefim[iref].get_attr('n3'))
				if(diff_projdir(n1, phi, theta) >= ant):
					retvals = Util.Crosrng_ms(crefim[iref], cimage, numr)
					#print  retvals
					#from sys import exit
					#exit()
					#for l in xrange(6000):  retvals = Util.Crosrng_ms(crefim[iref], cimage, numr)
					#print  "Crosrng  ",ttime()
					qn = retvals["qn"]
					qm = retvals["qm"]
					if(qn >=peak or qm >=peak):
						sx = -ix
						sy = -iy
						nref = iref
						if (qn >= qm):
							ang = ang_n(retvals["tot"], mode, numr[len(numr)-1])
							peak=qn
							mirror=0
						else:
							ang = ang_n(retvals["tmt"], mode, numr[len(numr)-1])
							peak=qm
							mirror=1
	if(peak == -1.0E23):
		ang = 0; sxs = 0; sys = 0; mirror = 0;
		nref = -1
	else:
		co =  cos(ang*pi/180.0)
		so = -sin(ang*pi/180.0)
		sxs = sx*co - sy*so
		sys = sx*so + sy*co
	return  ang, sxs, sys, mirror, peak, nref

def apmqs(image, refim, first_ring, last_ring, rstep, xrng, yrng, step, mode = "F"):
	"""Determine shift and rotation between images and many reference images (refim)
		quadratic interpolation
	"""
	from utilities import model_circle
	#from utilities import info, ttime
	#print "ORMQ"
	#print  ttime()
	mode = "F"
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr = ringwe(numr, mode)
	crefim = []
	nx=refim[0].get_xsize()
	ny=refim[0].get_ysize()
	mask = model_circle(last_ring, nx, ny)
	cnx = int(nx/2)+1
	cny = int(ny/2)+1
	for j in xrange(len(refim)):
		cimage = Util.Polar2Dm(refim[j].process_inplace("normalize.mask", {"mask":mask, "no_sigma":1}), cnx, cny, numr, mode)
		#cimage = Util.Polar2Dm(refim[j], cnx, cny, numr, mode)
		Util.Frngs(cimage, numr)
		Applyws(cimage, numr, wr)
		crefim.append(cimage)

	outparams = []
	for im in xrange(len(image)):	
		outparams.append(apmq(image[im], crefim, xrng, yrng, step, mode, numr, cnx, cny))
	return  outparams
	
def apnq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny):
	"""Determine shift and rotation between image and many reference images (crefim, weights have to be applied)
		quadratic interpolation
		no MIRROR
	"""
	from math import pi, cos, sin
	#print "APNQ"
	peak = -1.0E23
	ky = int(2*yrng/step+0.5)//2
	kx = int(2*xrng/step+0.5)//2
	for i in xrange(-ky, ky+1):
		iy = i*step
		for  j in xrange(-kx, kx+1):
			ix = j*step
			cimage=Util.Polar2Dm(image, cnx+ix, cny+iy, numr, mode)
			Util.Frngs(cimage, numr)
			#  compare with all reference images
			for iref in xrange(len(crefim)):
				retvals = Util.Crosrng_e(crefim[iref], cimage, numr, 0)
				qn = retvals["qn"]
				if(qn >= peak):
					sx = -ix
					sy = -iy
					nref = iref
					ang = ang_n(retvals["tot"], mode, numr[len(numr)-1])
					peak=qn
	mirror=0
	co =  cos(ang*pi/180.0)
	so = -sin(ang*pi/180.0)
	sxs = sx*co - sy*so
	sys = sx*so + sy*co
	return  ang, sxs, sys, mirror, peak, nref

def apply_trans2d(data, kb):
	"""
		Apply alignment parameters
	"""
	from fundamentals import rotshift2dg
	n = len(data)
	nx = data[0].get_xsize()
	ny = data[0].get_ysize()
	dala = []
	for i in xrange(n):
		alpha  = data[i].get_attr('alpha')
		sx     = data[i].get_attr('sx')
		sy     = data[i].get_attr('sy')
		mirror = data[i].get_attr('mirror')
		temp = rotshift2dg(data[i], alpha, sx, sy, kb)
		if  mirror: temp.process_inplace("mirror",{"axis":'x'})
		dala.append(temp)
	return dala


def aptrans(image, params):
	""""
	  apply spider style transformation parameters to an image
	"""
	from fundamentals import rtshg
	output = rtshg(image, params[0], params[1], params[2])
	if params[3]: output.process_inplace("mirror", {"axis":'x'})
	return  output


	
def combine_params_new(spider_trans,ref_ali_params):
	"""
	  combine spider style transformation doc file (alpha,sx,sy,mirrory)
	  with the results of ref_ali (phi,theta,psi,nsx,nsy)
	  Output:
		angles: in spider style (psi,theta,phi)
	"""
	angles=[]
	new_trans=[]
	for i in xrange(0,len(ref_ali_params)/5):
		angles.append([ref_ali_params[i*5],ref_ali_params[i*5+1],ref_ali_params[i*5+2]])
		if(spider_trans[i*4+3]==0):
			#no mirror, simply add shifts
			new_trans.append([spider_trans[i*4],spider_trans[i*4+1]+ref_ali_params[i*5+2],spider_trans[i*4+2]+ref_ali_params[i*5+3],spider_trans[i*4+3]])
		else:
			#mirror, subtract shifts
			new_trans.append([spider_trans[i*4],spider_trans[i*4+1]-ref_ali_params[i*5+2],spider_trans[i*4+2]-ref_ali_params[i*5+3],spider_trans[i*4+3]])
	return  angles,new_trans

def dfunc(args, data):
	from math import pi
	res = data[0].rot_scale_conv7(args[0]*pi/180., args[1], args[2], data[2], 1.0)
	#v = res.cmp("ccc", data[1], {"mask":data[3], "negative":0})
	v = -res.cmp("SqEuclidean", data[1], {"mask":data[3]})
	return v

def dfunc2(args, data):
	res = data[0].rot_scale_trans2D(args[0], args[1], args[2], 1.0)
	#v = res.cmp("ccc", data[1], {"mask":data[2], "negative":0})
	v = -res.cmp("SqEuclidean", data[1], {"mask":data[2]})
	return v

def dfunc_i(args, data):
	from math import pi
	res = data[0].rot_scale_conv7(args[0]*pi/180., args[1], args[2], data[2], 1.0)
	v = res.cmp("SqEuclidean", data[1], {"mask":data[3]})
	return v

def dfunc2_i(args, data):
	res = data[0].rot_scale_trans2D(args[0], args[1], args[2], 1.0)
	v = res.cmp("SqEuclidean", data[1], {"mask":data[2]})
	return v
		
def eqproj_spin(args, data):
	from projection import prgs
	from math import pi
	#print  " AMOEBA ",args
	#  data: 0 - volkb,  1 - kb, 2 - image,  3 - mask, 
	EULER_SPIN=Transform3D.EulerType.SPIN
	EULER_SPIDER = Transform3D.EulerType.SPIDER
	dict_spin = dict()
	dict_spin["Omega"]=args[0] 
	dict_spin["n1"]=sin(args[1]*pi/180.)*cos(args[2]*pi/180.)
	dict_spin["n2"]=sin(args[1]*pi/180.)*sin(args[2]*pi/180.)
	dict_spin["n3"]=cos(args[1]*pi/180.)
	Rspin=Transform3D(EULER_SPIN,dict_spin)
	vv=Rspin.get_rotation(EULER_SPIDER)
	params=5*[0.0]
	params[0]=vv["phi"]
	params[1]=vv["theta"]
	params[2]=vv["psi"]
	params[3]=args[3]
	params[4]=args[4]

	proj=prgs(data[0], data[1], params)
	v = -proj.cmp("SqEuclidean", data[2], {"mask":data[3]})
	#print  " AMOEBA o",params,v
	return v

def eqproj_(args, data):
	from projection import prgs
	#print  " AMOEBA ",args
	#  data: 0 - volkb,  1 - kb, 2 - image,  3 - mask, 
	proj = prgs(data[0], data[1], args)
	v = -proj.cmp("SqEuclidean", data[2], {"mask":data[3]})
	#print  " AMOEBA o",args,v
	return v

def eqproj_xyx(args,data):
	from projection import prgs
	#print  " AMOEBA ",args
	#  data: 0 - volkb,  1 - kb, 2 - image,  3 - mask, 
	proj=prgs(data[0], data[1], args)
	v = -proj.cmp("SqEuclidean", data[2], {"mask":data[3]})
	#print  " AMOEBA o",args,v
	return v
    

def eqproj2(args, data):
	from projection import prgs
	from fundamentals import fft
	from utilities import amoeba_multi_level
	from math import cos, sin, pi
	
	proj = prgs(data[0], data[1], args+[0.0, 0.0, 0.0])
	
	M = data[2].get_xsize()
	sx = data[4][1]
	sy = data[4][2]
	crefim = data[5]
	mode = data[6][0]	
	numr = data[6][1]
	cnx = data[6][2]
	cny = data[6][3]
	kb = data[6][4]
	params = data[6][5]
	maxrin = data[6][6]
	kbline = data[6][7]
	parline = data[6][8]
	
	imali=proj.FourInterpol(2*M,2*M,1,0)
	q=Processor.EMFourierFilter(imali,params)
	imali=fft(q)
	
	is_mirror = 0
		
	data0 = []
	data0.insert(0,imali)
	data0.insert(1,cnx)
	data0.insert(2,cny)
	data0.insert(3,numr)
	data0.insert(4,mode)
	data0.insert(5,kb)
	data0.insert(6,crefim)
	data0.insert(7,parline)
	data0.insert(8,maxrin)
	data0.insert(9,kbline)
	data0.insert(10,is_mirror)
	
	ps2 = amoeba_multi_level([sx,sy],[1,1],func_loop2,1.e-4,1.e-4,500,data0)	
		
	sx = -ps2[0][0]
	sy = -ps2[0][1]
	ang = ps2[3]
			
	co =  cos(ang*pi/180.0)
	so = -sin(ang*pi/180.0)
	sxs = sx*co - sy*so
	sys = sx*so + sy*co

	proj = prgs(data[0], data[1], args+[ang, sxs, sys])
	v = proj.cmp("ccc", data[2], {"mask":data[3], "negative":0})
	return v, [ang, sxs, sys]

def eqproj_xyz(args,data):
    from projection import prgs
    #print  " AMOEBA ",args
    #  data: 0 - volkb,  1 - kb, 2 - image,  3 - mask, 
    EULER_XYZ=Transform3D.EulerType.XYZ
    EULER_SPIDER = Transform3D.EulerType.SPIDER
    dict_xyz = dict()
    dict_xyz["xtilt"]=args[0]
    dict_xyz["ytilt"]=args[1]
    dict_xyz["ztilt"]=args[2]
    Rxyz=Transform3D(EULER_XYZ,dict_xyz)
    vv=Rxyz.get_rotation(EULER_SPIDER)
    params=5*[0.0]
    params[0]=vv["phi"]
    params[1]=vv["theta"]
    params[2]=vv["psi"]
    params[3]=args[3]
    params[4]=args[4]

    proj=prgs(data[0], data[1], params)
    v = -proj.cmp("SqEuclidean", data[2], {"mask":data[3]})
    #print  " AMOEBA o",args,v
    return v

def orm(image, refim, xrng=0, yrng=0, step=1, first_ring=1, last_ring=0, rstep=1, mode = "F"):
	"""  Determine shift and rotation between image and reference image
	     quadratic interpolation
	"""
	nx=refim.get_xsize()
	ny=refim.get_ysize()
	if(last_ring == 0):  last_ring = nx//2-2-max(xrng,yrng)
	# center in SPIDER convention
	cnx = int(nx/2)+1
	cny = int(ny/2)+1
	#precalculate rings
	numr=Numrinit(first_ring, last_ring, rstep, mode)
	wr=ringwe(numr ,mode)
	#for iring in xrange(len(wr)):
	#	print wr[iring], iring
	#cimage=Util.Polar2Dmi(refim, cnx, cny, numr, mode, kb)
	cimage=Util.Polar2Dm(refim, cnx, cny, numr, mode)
	Util.Frngs(cimage, numr)
	Applyws(cimage, numr, wr)
	return ormq(image, cimage, xrng, yrng, step, mode, numr, cnx, cny)
       
def ormqi(image, refim, first_ring, last_ring, rstep=1, xrng=0, yrng=0, step=1, mode="F"):
	"""Determine shift and rotation between image and reference image (refim)
		using gridding resampling
	"""
	from math import pi, cos, sin
	from fundamentals import fft
	#from utilities import info
	#print "ORMQI"
	numr=Numrinit(first_ring, last_ring, rstep, mode)
	wr = ringwe(numr,mode)
	nx = image.get_xsize()
	ny = image.get_ysize()
	cnx = nx/2+1
	cny = ny/2+1

	M=image.get_xsize()
	# padd two times
	npad=2
	N=M*npad
	# support of the window
	K=6
	alpha=1.75
	r=M/2
	v=K/2.0/N
	kb = Util.KaiserBessel(alpha, K, r, K/(2.*N), N)
	oo=refim.FourInterpol(2*M,2*M,1,0)
	params = {"filter_type" : Processor.fourier_filter_types.KAISER_SINH_INVERSE,"alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
	q=Processor.EMFourierFilter(oo,params)
	oo=fft(q)
	crefim=Util.Polar2Dmi(oo, cnx, cny, numr, mode ,kb)
	#info(crefim)
	Util.Frngs(crefim, numr)
	Applyws(crefim, numr, wr)

	oo=image.FourInterpol(2*M,2*M,1,0)
	params = {"filter_type" : Processor.fourier_filter_types.KAISER_SINH_INVERSE,"alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
	q=Processor.EMFourierFilter(oo,params)
	oo=fft(q)
	 
	peak = -1.0E23
	for iy in xrange(-yrng, yrng+1, step):
		for  ix in xrange(-xrng, xrng+1, step):
			cimage=Util.Polar2Dmi(oo, cnx+ix, cny+iy, numr, mode, kb)
			Util.Frngs(cimage, numr)

			#retvals = Util.Crosrng_e(crefim, cimage, numr, 0)
			#print  "  e 0 ",retvals["tot"],retvals["qn"]
			#retvals = Util.Crosrng_e(crefim, cimage, numr, 1)
			#print  "  e 1 ",retvals["tot"],retvals["qn"]
			retvals = Util.Crosrng_ms(crefim, cimage, numr)
			#print  "  ms  ",retvals["tot"],retvals["qn"],retvals["tmt"],retvals["qm"]
			qn = retvals["qn"]
			qm = retvals["qm"]
			if(qn >=peak or qm >=peak):
				sx = -ix
				sy = -iy
				if (qn >= qm):
					ang = ang_n(retvals["tot"], mode, numr[len(numr)-1])
					peak=qn
					mirror=0
				else:
					ang = ang_n(retvals["tmt"], mode, numr[len(numr)-1])
					peak=qm
					mirror=1
	co =  cos(ang*pi/180.0)
	so = -sin(ang*pi/180.0)
	sxs = sx*co - sy*so
	sys = sx*so + sy*co
	return  ang,sxs,sys,mirror,peak
       
def ormqip(image,crefim,xrng,yrng,step,mode,numr,kb,cnx,cny,nx):
	"""Determine shift and rotation between image (kb prepared) and reference image (crefim)
	"""
	from math import pi,cos,sin
	from fundamentals import fft
	#from utilities import info
	#print "ORMQIP"
	peak = -1.0E23
	for iy in xrange(-yrng,yrng+1,step):
		for  ix in xrange(-xrng,xrng+1,step):
			cimage = Util.Polar2Dmi(image, cnx+ix, cny+iy, numr, mode, kb)
			#cimage=Util.Polar2Dm(image,cnx+ix,cny+iy,numr,mode)
			#info(cimage)
			Util.Frngs(cimage, numr)
			#info(cimage)
			retvals = Util.Crosrng_ms(crefim,cimage,numr)
			#print  retvals["qn"],retvals["tot"],retvals["qm"],retvals["tmt"]
			qn = retvals["qn"]
			qm = retvals["qm"]
			if(qn >=peak or qm >=peak):
				sx = -ix
				sy = -iy
				if (qn >= qm):
					ang = ang_n(retvals["tot"], mode, numr[len(numr)-1])
					peak=qn
					mirror=0
				else:
					ang = ang_n(retvals["tmt"], mode, numr[len(numr)-1])
					peak=qm
					mirror=1
	co =  cos(ang*pi/180.0)
	so = -sin(ang*pi/180.0)
	sxs = sx*co - sy*so
	sys = sx*so + sy*co
	return  ang,sxs,sys,mirror,peak
       
def ormqg(image,refim,first_ring,last_ring,rstep,xrng,yrng,step,mode):
     """Determine shift and rotation between image and reference image (refim)
     """
     from math import pi
     from fundamentals import fft
     #from utilities import info
     #print "ORMQG"
     numr=Numrinit(first_ring,last_ring,rstep,mode)
     wr = ringwe(numr,mode)
     nx = image.get_xsize()
     ny = image.get_ysize()
     cnx = nx/2+1
     cny = ny/2+1
     
     M=image.get_xsize()
     # padd two times
     npad=2
     N=M*npad
     # support of the window
     K=6
     alpha=1.75
     r=M/2
     v=K/2.0/N
     kb = Util.KaiserBessel(alpha, K, r, K/(2.*N), N)
     oo=refim.FourInterpol(N,N,1,0)
     params = {"filter_type" : Processor.fourier_filter_types.KAISER_SINH_INVERSE,"alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
     q=Processor.EMFourierFilter(oo,params)
     oo=fft(q)
     crefim=Util.Polar2Dmi(oo,cnx,cny,numr,mode,kb)
     #info(crefim)
     Util.Frngs(crefim, numr)
     Applyws(crefim, numr,wr)
     
     oo=image.FourInterpol(2*M,2*M,1,0)
     params = {"filter_type" : Processor.fourier_filter_types.KAISER_SINH_INVERSE,"alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
     q=Processor.EMFourierFilter(oo,params)
     oo=fft(q)
     maxrin=numr[len(numr)-1]
     #print  maxrin
     #  peak_s is a 3D ccf array straight
     peak_s = EMData()
     peak_s.set_size(maxrin,2*xrng+1,2*yrng+1)
     #  peak_m is a 3D ccf array mirrored
     peak_m = EMData()
     peak_m.set_size(maxrin,2*xrng+1,2*yrng+1)
     for iy in xrange(-yrng,yrng+1,step):
       for  ix in xrange(-xrng,xrng+1,step):
	 cimage=Util.Polar2Dmi(oo,cnx+ix,cny+iy,numr,mode,kb)
	 #info(cimage)
	 Util.Frngs(cimage, numr)
	 qt = Util.Crosrng_msg(crefim,cimage,numr)
	 for k in xrange(0,maxrin):
	   peak_s[(k+maxrin/2)%maxrin,ix+xrng,iy+yrng]=qt[k,0]
	   peak_m[(k+maxrin/2)%maxrin,ix+xrng,iy+yrng]=qt[k,1]
     return  peak_s,peak_m

"""
def ormql(image,refim,first_ring,last_ring,rstep,xrng,yrng,step,mode):
	
	#Determine shift and rotation between image and reference image (refim)
	
	from math import pi, cos, sin
	from fundamentals import fft
	from utilities import amoeba, model_circle
	#print "ORMQL"
	numr=Numrinit(first_ring,last_ring,rstep,mode)
	wr = ringwe(numr,mode)
	nx = image.get_xsize()
	ny = image.get_ysize()
	cnx = nx/2+1
	cny = ny/2+1

	M=image.get_xsize()
	# padd two times
	npad=2
	N=M*npad
	# support of the window
	K=6
	alpha=1.75
	r=M/2
	v=K/2.0/N
	kb = Util.KaiserBessel(alpha, K, r, K/(2.*N), N)
	refi=refim.FourInterpol(N,N,1,0)
	params = {"filter_type" : Processor.fourier_filter_types.KAISER_SINH_INVERSE,"alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
	q=Processor.EMFourierFilter(refi,params)
	refi=fft(q)
	crefim=Util.Polar2Dmi(refi,cnx,cny,numr,mode,kb)
	#info(crefim)
	Util.Frngs(crefim, numr)
	Applyws(crefim, numr,wr)

	imali=image.FourInterpol(2*M,2*M,1,0)
	params = {"filter_type" : Processor.fourier_filter_types.KAISER_SINH_INVERSE,"alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
	q=Processor.EMFourierFilter(imali,params)
	imali=fft(q)
	maxrin=numr[len(numr)-1]

	#print  maxrin
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
	data=[]
	data.insert(1,kbline)
	#dama=[]
	#dama.insert(1,kbline)

	peak = -1.0E23
	#qm=peak  # remove later
	for iy in xrange(-yrng,yrng+1,step):
		for  ix in xrange(-xrng,xrng+1,step):
			cimage=Util.Polar2Dmi(imali,cnx+ix,cny+iy,numr,mode,kb)
			Util.Frngs(cimage, numr)
			qt = Util.Crosrng_msg(crefim,cimage,numr)
			# straight
			for i in xrange(0,maxrin): line[i]=qt[i,0]

			#  find the current maximum and its location
			ps=line.peak_search(1,1)
			#print  ps
			qn=ps[1]
			jtot=ps[2]  #/2
			q=Processor.EMFourierFilter(line,parline)
			data.insert(0,q)
			#info(q)
			#print  q.peak_search(1,1)
			#print  ix,iy
			#print  " before amoeba",ps[2]/2,q.get_pixel_conv(ps[2]/2, 0, 0, data[1])
			ps = amoeba([jtot], [2.0],afunc,1.e-4,1.e-4,500,data)
			del data[0]
			#print  ps
			jtot=ps[0][0]*2
			qn=ps[1]
			#print  qn,2*jtot*360/maxrin
			# mirror
			for i in xrange(0,maxrin): line[i]=qt[i,1]
			#print " MIRROR"
			#  find the current maximum and its location
			ps=line.peak_search(1,1)
			qm=ps[1]
			mtot=ps[3]/2
			q=Processor.EMFourierFilter(line,parline)
			data.insert(0,q)
			#print qm,mtot
			ps = amoeba([mtot], [2.0],afunc,1.e-4,1.e-4,500,data)
			del data[0]
			#print  qm,ps
			mtot=ps[0][0]*2
			qm=ps[1]
			#print  qm,mtot
	 
			if(qn >=peak or qm >=peak):
				if (qn >= qm):
					ang = ang_n(jtot+1, mode, numr[len(numr)-1])
					sx = -ix
					sy = -iy
					peak=qn
					mirror=0
				else:
					ang = ang_n(mtot+1, mode, numr[len(numr)-1])
					sx = -ix
					sy = -iy
					peak=qm
					mirror=1

      
	co =  cos(ang*pi/180.0)
	so = -sin(ang*pi/180.0)
	sxs = sx*co - sy*so
	sys = sx*so + sy*co
	return  ang,sxs,sys,mirror,peak    
	
	#postprocessing by amoeba search	
	del data[0]
	mask=model_circle(last_ring,nx,ny)
	data.insert(0,imali)
	data.insert(1,refim)
	data.insert(2,kb)
	data.insert(3,mask)
	ps = amoeba([ang,sxs,sys], [360.0/(2.0*pi*last_ring),1.0,1.0],dfunc,1.e-4,1.e-4,500,data)
	#print  ps
	#return  ps[0][0],ps[0][1],ps[0][2],mirror,-ps[1]
"""

def ormql(image,crefim,xrng,yrng,step,mode,numr,cnx,cny):
	"""
	Determine shift and rotation between image and reference image (refim)
	"""
	from math import pi, cos, sin
	from fundamentals import fft
	from utilities import amoeba, model_circle
	
	M=image.get_xsize()
	npad=2
	N=M*npad
	#support of the window
	K=6
	alpha=1.75
	r=M/2
	v=K/2.0/N
	kb = Util.KaiserBessel(alpha, K, r, K/(2.*N), N)

	imali=image.FourInterpol(2*M,2*M,1,0)
	params = {"filter_type" : Processor.fourier_filter_types.KAISER_SINH_INVERSE,"alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
	q=Processor.EMFourierFilter(imali,params)
	imali=fft(q)
	
	maxrin=numr[len(numr)-1]

	#print  maxrin
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
	data=[]
	data.insert(1,kbline)

	peak = -1.0E23
	#qm=peak  # remove later
	for iy in xrange(-yrng,yrng+1,step):
		for  ix in xrange(-xrng,xrng+1,step):
			cimage=Util.Polar2Dmi(imali,cnx+ix,cny+iy,numr,mode,kb)
			Util.Frngs(cimage, numr)
			qt = Util.Crosrng_msg(crefim,cimage,numr)
			
			# straight
			for i in xrange(0,maxrin): line[i]=qt[i,0]								
			#  find the current maximum and its location
			ps=line.peak_search(1,1)
			qn=ps[1]
			jtot=ps[2]/2
			q=Processor.EMFourierFilter(line,parline)
			data.insert(0,q)
			#info(q)
			#print  q.peak_search(1,1)
			#print  ix,iy
			#print  " before amoeba",ps[2]/2,q.get_pixel_conv(ps[2]/2, 0, 0, data[1])
			#print  " before amoeba",ps[2] ,q.get_pixel_conv(ps[2] , 0, 0, data[1])
			ps = amoeba([jtot], [2.0], oned_search_func,1.e-4,1.e-4,500,data)
			del data[0]
			#print  ps
			jtot=ps[0][0]*2
			qn=ps[1]
			#print  qn,2*jtot*360/maxrin
			
			# mirror
			for i in xrange(0,maxrin): line[i]=qt[i,1]
			#print " MIRROR"
			#  find the current maximum and its location
			ps=line.peak_search(1,1)
			# print ps
			qm=ps[1]
			mtot=ps[2]/2
			q=Processor.EMFourierFilter(line,parline)
			data.insert(0,q)
			#print qm,mtot
			ps = amoeba([mtot], [2.0],oned_search_func,1.e-4,1.e-4,500,data)
			del data[0]
			#print  qm,ps
			mtot=ps[0][0]*2
			qm=ps[1]
			#print  qm,mtot	 

			if(qn >=peak or qm >=peak):
				if (qn >= qm):
					ang = ang_n(jtot+1, mode, numr[len(numr)-1])
					sx = -ix
					sy = -iy
					peak=qn
					mirror=0
				else:
					ang = ang_n(mtot+1, mode, numr[len(numr)-1])
					sx = -ix
					sy = -iy
					peak=qm
					mirror=1
     
	co =  cos(ang*pi/180.0)
	so = -sin(ang*pi/180.0)
	sxs = sx*co - sy*so
	sys = sx*so + sy*co
	return  ang,sxs,sys,mirror,peak    	


def ormq_all(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, interpolation_method=None, opti_algorithm=None, opti_strategy=None, refim=None):
	"""
	Determine shift and rotation between image and reference image (refim)
	One can specify the interpolation method, optimization algorithm and optimization
	strategy for the routine.
	It should be noted that crefim should be preprocessed before using this	function.
	Also, if using "concurrent" as optimization strategy, one should provide the original
	reference image
	"""
	from math import pi, cos, sin
	from fundamentals import fft, mirror
	from utilities import amoeba, amoeba_multi_level, model_circle
	from alignment import ang_n
	from sys import exit
	
	if interpolation_method == None:
		interpolation_method = interpolation_method_2D #use global setting
	
	nx = image.get_xsize()
	maxrin = numr[-1]

	peak = -1.0E23
	ky = int(2*yrng/step+0.5)/2
	kx = int(2*xrng/step+0.5)/2

	if interpolation_method == "gridding":
		npad = 2
		N = nx*npad
		#support of the window
		K = 6
		alpha = 1.75
		r = nx/2
		v = K/2.0/N
		kb = Util.KaiserBessel(alpha, K, r, v, N)
		imali = image.FourInterpol(2*nx, 2*nx, 1, 0)
		params = {"filter_type":Processor.fourier_filter_types.KAISER_SINH_INVERSE, "alpha":alpha, "K":K, "r":r, "v":v, "N":N}
		del alpha, K, r, v, N
		q = Processor.EMFourierFilter(imali,params)
		imali = fft(q)
	
		# do not pad
		npad = 1
		N = maxrin*npad
		# support of the window
		K = 6
		alpha = 1.75
		r = maxrin/2
		v = K/2.0/N
		kbline = Util.KaiserBessel(alpha, K, r, v, N)
		parline = {"filter_type":Processor.fourier_filter_types.KAISER_SINH_INVERSE, "alpha":alpha, "K":K, "r":r, "v":v, "N":N}
		del alpha, K, r, v, N
		data = [kbline]

		for i in xrange(-ky, ky+1):
			iy = i*step
			for j in xrange(-kx, kx+1):
				ix = j*step
				cimage = Util.Polar2Dmi(imali, cnx+ix, cny+iy, numr, mode, kb)
				Util.Frngs(cimage, numr)
				
				# straight
				line_s = Util.Crosrng_msg_s(crefim, cimage, numr)
				#  find the current maximum and its location
				ps = line_s.max_search()
				jtot = ps[0]/2
				q = Processor.EMFourierFilter(line_s, parline)
				data.insert(0, q)
				ps = amoeba([jtot], [2.0], oned_search_func, 1.e-4, 1.e-4, 500, data)
				del data[0]
				jtot = ps[0][0]*2
				qn = ps[1]
				
				# mirror
				line_m = Util.Crosrng_msg_m(crefim, cimage, numr)
				# find the current maximum and its location
				ps = line_m.max_search()
				mtot = ps[0]/2
				q = Processor.EMFourierFilter(line_m, parline)
				data.insert(0, q)
				ps = amoeba([mtot], [2.0], oned_search_func, 1.e-4, 1.e-4, 500, data)
				del data[0]
				mtot = ps[0][0]*2
				qm = ps[1]

				if qn >= peak or qm >= peak:
					if qn >= qm:
						ang = ang_n(jtot+1, mode, maxrin)
						"""
						NOTICE: This is important!!!!
						
						The reason that there is a "+1" here while there is no "+1" in quadratic is that
						in program such as Crosrng_ms (used in quadratic), the returned position begins
						from 1 insteat of 0. Although this is kind of inconvenient, we decide to keep this 
						convention. Therefore, in function ang_n, there is a "-1" for offsetting. However, 
						for gridding, the program such as Crosrng_msg returns ccf. Also, we choose to have 
						the found position begin from 0, otherwise, it will become very complicated
						since it involves "*2" and "/2". Therefore, when using ang_n, we should have another
						"+1" to offset back.
						"""
						sx = -ix
						sy = -iy
						peak = qn
						is_mirror = 0
					else:
						ang = ang_n(mtot+1, mode, maxrin)
						sx = -ix
						sy = -iy
						peak = qm
						is_mirror = 1
		if opti_algorithm == None:
			co =  cos(ang*pi/180.0)
			so = -sin(ang*pi/180.0)
			sxs = sx*co - sy*so
			sys = sx*so + sy*co
			return  ang, sxs, sys, is_mirror, peak
		elif opti_strategy == "concurrent":
			co =  cos(ang*pi/180.0)
			so = -sin(ang*pi/180.0)
			sxs = sx*co - sy*so
			sys = sx*so + sy*co
			last_ring = numr[-3]
			mask = model_circle(last_ring, nx, nx)
		
			if refim == None:
				print "For concurrent strategy, you need to provide referenece image as well as reference image in polar coordinates!"
				exit()
			if opti_algorithm == "amoeba":
				if is_mirror:	
					image2 = mirror(image)
					imali = image2.FourInterpol(2*nx, 2*nx, 1, 0)
					q = Processor.EMFourierFilter(imali, params)
					imali = fft(q)
					ang = -ang
					sxs = -sxs
				data = []
				data.insert(0, imali)
				data.insert(1, refim)
				data.insert(2, kb)
				data.insert(3, mask)
				ps = amoeba([ang, sxs, sys], [360.0/(2.0*pi*last_ring), step, step], dfunc, 1.e-4, 1.e-4, 500, data)

				ps[0][0] = ps[0][0]%360.0
				if is_mirror:
					ps[0][0] = -ps[0][0]
					ps[0][1] = -ps[0][1]
				return  ps[0][0], ps[0][1], ps[0][2], is_mirror, -ps[1]
			elif opti_algorithm == "LBFGSB":
				from scipy.optimize.lbfgsb import fmin_l_bfgs_b
				if is_mirror:
					image2 = mirror(image)
					imali = image2.FourInterpol(2*nx, 2*nx, 1, 0)
					q = Processor.EMFourierFilter(imali, params)
					imali = fft(q)
					x0 = [-ang, -sxs, sys]
				else:
					x0 = [ang, sxs, sys]
				bounds = [(x0[0]-2.0, x0[0]+2.0), (x0[1]-1.5, x0[1]+1.5), (x0[2]-1.5, x0[2]+1.5)]
				data = []
				data.append(imali)
				data.append(refim)
				data.append(kb)
				data.append(mask)
				ps, val, d = fmin_l_bfgs_b(dfunc_i, x0, args=[data], approx_grad=1, bounds=bounds, m=10, factr=1e1, pgtol=1e-5, epsilon=1e-3, iprint=-1, maxfun=15000)
				ps[0] = ps[0]%360.0
				if is_mirror:
					ps[0] = -ps[0]
					ps[1] = -ps[1]
	
				return  ps[0], ps[1], ps[2], is_mirror, val
			elif opti_algorithm == "SD":
				if is_mirror:	
					image2 = mirror(image)
					imali = image2.FourInterpol(2*nx, 2*nx, 1, 0)
					q = Processor.EMFourierFilter(imali, params)
					imali = fft(q)
					ang = -ang
					sxs = -sxs
				ps = Util.twoD_fine_ali_SD_G(imali, refim, mask, kb, ang, sxs, sys)
				ps[0] = ps[0]%360.0
				if is_mirror:
					ps[0] = -ps[0]
					ps[1] = -ps[1]
				return ps[0], ps[1], ps[2], is_mirror, ps[3]
			elif opti_algorithm == "DUD":
				if is_mirror:
					image2 = mirror(image)
					imali = image2.FourInterpol(2*nx, 2*nx, 1, 0)
					q = Processor.EMFourierFilter(imali, params)
					imali = fft(q)
					ang = -ang
					sxs = -sxs
				
				P = 3
				max_P = 20
				converge = False
				beta = []
				delta = 0.1
				beta.append([ang, sxs, sys])
				beta.append([ang+delta, sxs, sys])
				beta.append([ang, sxs+delta, sys])
				beta.append([ang, sxs, sys+delta])

				from numpy import matrix, diag, concatenate
				from utilities import get_image_data

				Y = matrix(get_image_data(refim)).flatten(1)
				omega = []
				imali0 = imali.rot_scale_conv7(beta[0][0]*pi/180, beta[0][1], beta[0][2], kb, 1.0)
				for p in xrange(1, P+1):
					imalip = imali.rot_scale_conv7(beta[p][0]*pi/180, beta[p][1], beta[p][2], kb, 1.0)
					diffim = Util.subn_img(imalip, imali0)
					P2 = matrix(get_image_data(diffim)).flatten(1)
					if p == 1: F = P2
					else: F = concatenate((F, P2))
					omega.append(imalip.cmp("SqEuclidean", imali0))
					
				imalip = imali.rot_scale_conv7(beta[P][0]*pi/180, beta[P][1], beta[P][2], kb, 1.0)
				Fs = matrix(get_image_data(imalip)).flatten(1)

				A = [[delta, 0, 0], [0, delta, 0], [0, 0, delta]]
				while converge == False and P < max_P:
					A_matrix = matrix(A)
					omega_matrix_inv = matrix(diag(omega))**(-1)
					G = F.T*omega_matrix_inv*A_matrix*(A_matrix.T*omega_matrix_inv*A_matrix)**(-1)
					
					H = (G.T*G)**(-1)*G.T*(Y.T-Fs.T)

					dalpha = float(H[0][0])
					dsx = float(H[1][0])
					dsy = float(H[2][0])
					if dalpha**2+dsx**2+dsy**2 < 1e-6: converge = True
					beta.append([float(beta[P][0])+dalpha, float(beta[P][1])+dsx, float(beta[P][2])+dsy])
					P += 1
					if converge == False:
						imali_new = imali.rot_scale_conv7(float(beta[P][0])*pi/180, float(beta[P][1]), float(beta[P][2]), kb, 1.0)
						diffim = Util.subn_img(imali_new, imali0)
						P2 = matrix(get_image_data(diffim)).flatten(1)
						F = concatenate((F, P2))
						Fs = matrix(get_image_data(imali_new)).flatten(1)
						omega.append(imali_new.cmp("SqEuclidean", imali0))
						A.append([dalpha, dsx, dsy])
				#print "P =", P, "converge =", converge

				ang = float(beta[P][0])
				sxs = float(beta[P][1])
				sys = float(beta[P][2])
				ang = ang%360.0
				if is_mirror:
					ang = -ang
					sxs = -sxs
				return ang, sxs, sys, is_mirror, 0
			else:
				print "Unknown optimization algorithm!"
		elif opti_strategy == "cascaded":
			data0 = []
			data0.insert(0, imali)
			data0.insert(1, cnx)
			data0.insert(2, cny)
			data0.insert(3, numr)
			data0.insert(4, mode)
			data0.insert(5, kb)
			data0.insert(6, crefim)
			data0.insert(7, parline)
			data0.insert(8, maxrin)
			data0.insert(9, kbline)
			data0.insert(10, is_mirror)

			ps2 = amoeba_multi_level([-sx, -sy], [1, 1], func_loop2, 1.e-4, 1.e-4, 50, data0)
		
			if ps2[1] > peak:
				sx = -ps2[0][0]
				sy = -ps2[0][1]
				ang = ps2[3]
				peak = ps2[1]
			co =  cos(ang*pi/180.0)
			so = -sin(ang*pi/180.0)
			sxs = sx*co - sy*so
			sys = sx*so + sy*co
			return  ang, sxs, sys, is_mirror, peak	
		else: print "Unknown optimization strategy!"
		
	elif interpolation_method == "quadratic":
		for i in xrange(-ky, ky+1):
			iy = i*step
			for  j in xrange(-kx, kx+1):
				ix = j*step
				cimage = Util.Polar2Dm(image, cnx+ix, cny+iy, numr, mode)
				Util.Frngs(cimage, numr)
				retvals = Util.Crosrng_ms(crefim, cimage, numr)
				qn = retvals["qn"]
				qm = retvals["qm"]
				if qn >= peak or qm >= peak:
					sx = -ix
					sy = -iy
					if qn >= qm:
						ang = ang_n(retvals["tot"], mode, maxrin)
						peak = qn
						is_mirror = 0
					else:
						ang = ang_n(retvals["tmt"], mode, maxrin)
						peak = qm
						is_mirror = 1		
		if opti_algorithm == None:
			co =  cos(ang*pi/180.0)
			so = -sin(ang*pi/180.0)
			sxs = sx*co - sy*so
			sys = sx*so + sy*co
			return  ang, sxs, sys, is_mirror, peak
		elif opti_strategy == "concurrent":
			co =  cos(ang*pi/180.0)
			so = -sin(ang*pi/180.0)
			sxs = sx*co - sy*so
			sys = sx*so + sy*co
			last_ring = numr[-3]
			mask = model_circle(last_ring, nx, nx)

			if refim == None:
				print "For concurrent strategy, you need to provide referenece image as well as reference image in polar coordinates!"
				exit()
			if opti_algorithm == "amoeba":
				data = []
				if is_mirror:
					data.insert(0, mirror(image))
					ang = -ang
					sxs = -sxs
				else:
					data.insert(0, image)
				data.insert(1, refim)
				data.insert(2, mask)
				ps = amoeba([ang, sxs, sys], [360.0/(2.0*pi*last_ring), step, step], dfunc2, 1.e-4, 1.e-4, 500, data)
			
				ps[0][0] = ps[0][0]%360.0
				if is_mirror:
					ps[0][0] = -ps[0][0]
					ps[0][1] = -ps[0][1]
	
				return  ps[0][0], ps[0][1], ps[0][2], is_mirror, -ps[1]
			elif opti_algorithm == "LBFGSB":
				from scipy.optimize.lbfgsb import fmin_l_bfgs_b
				if is_mirror:
					image = mirror(image)
					x0 = [-ang, -sxs, sys] 
				else:
					x0 = [ang, sxs, sys] 
				bounds = [(x0[0]-2.0, x0[0]+2.0), (x0[1]-1.5, x0[1]+1.5), (x0[2]-1.5, x0[2]+1.5)]
				data = []
				data.append(image)
				data.append(refim)
				data.append(mask)			
				ps, val, d = fmin_l_bfgs_b(dfunc2_i, x0, args=[data], approx_grad=1, bounds=bounds, m=10, factr=1e1, pgtol=1e-5, epsilon=1e-3, iprint=-1, maxfun=15000)
				ps[0] = ps[0]%360.0
				if is_mirror:
					ps[0] = -ps[0]
					ps[1] = -ps[1]
				return  ps[0], ps[1], ps[2], is_mirror, val
			elif opti_algorithm == "SD":
				if is_mirror:
					image = mirror(image)
					ang = -ang
					sxs = -sxs
				ps = Util.twoD_fine_ali_SD(image, refim, mask, ang, sxs, sys)
				ps[0] = ps[0]%360.0
				if is_mirror:
					ps[0] = -ps[0]
					ps[1] = -ps[1]
				return ps[0], ps[1], ps[2], is_mirror, ps[3]
			elif opti_algorithm == "DUD":
				if is_mirror:
					image = mirror(image)
					ang = -ang
					sxs = -sxs
				P = 3
				max_P = 20
				converge = False
				beta = []
				delta = 0.1
				beta.append([ang, sxs, sys])
				beta.append([ang+delta, sxs, sys])
				beta.append([ang, sxs+delta, sys])
				beta.append([ang, sxs, sys+delta])

				from numpy import matrix, diag, concatenate
				from utilities import get_image_data

				Y = matrix(get_image_data(refim)).flatten(1)
				omega = []
				imali0 = image.rot_scale_trans2D(beta[0][0], beta[0][1], beta[0][2], 1.0)
				for p in xrange(1, P+1):
					imalip = image.rot_scale_trans2D(beta[p][0], beta[p][1], beta[p][2], 1.0)
					diffim = Util.subn_img(imalip, imali0)
					P2 = matrix(get_image_data(diffim)).flatten(1)
					if p == 1: F = P2
					else: F = concatenate((F, P2))
					omega.append(imalip.cmp("SqEuclidean", imali0))
					
				imalip = image.rot_scale_trans2D(beta[P][0], beta[P][1], beta[P][2], 1.0)
				Fs = matrix(get_image_data(imalip)).flatten(1)

				A = [[delta, 0, 0], [0, delta, 0], [0, 0, delta]]
				while converge == False and P < max_P:
					A_matrix = matrix(A)
					omega_matrix_inv = matrix(diag(omega))**(-1)
					G = F.T*omega_matrix_inv*A_matrix*(A_matrix.T*omega_matrix_inv*A_matrix)**(-1)
					
					H = (G.T*G)**(-1)*G.T*(Y.T-Fs.T)

					dalpha = float(H[0][0])
					dsx = float(H[1][0])
					dsy = float(H[2][0])
					if dalpha**2+dsx**2+dsy**2 < 1e-6: converge = True
					beta.append([float(beta[P][0])+dalpha, float(beta[P][1])+dsx, float(beta[P][2])+dsy])
					P += 1
					if converge == False:
						imali_new = image.rot_scale_trans2D(float(beta[P][0]), float(beta[P][1]), float(beta[P][2]), 1.0)
						diffim = Util.subn_img(imali_new, imali0)
						P2 = matrix(get_image_data(diffim)).flatten(1)
						F = concatenate((F, P2))
						Fs = matrix(get_image_data(imali_new)).flatten(1)
						omega.append(imali_new.cmp("SqEuclidean", imali0))
						A.append([dalpha, dsx, dsy])

				ang = float(beta[P][0])
				sxs = float(beta[P][1])
				sys = float(beta[P][2])
				ang = ang%360.0
				if is_mirror:
					ang = -ang
					sxs = -sxs
				return ang, sxs, sys, is_mirror, 0
			
			else:
				print "Unknown optimization algorithm!"
		elif opti_strategy == "cascaded":
			data = []
			data.insert(0, image)
			data.insert(1, cnx)
			data.insert(2, cny)
			data.insert(3, numr)
			data.insert(4, mode)
			data.insert(5, crefim)
			ps2 = amoeba_multi_level([-sx, -sy], [1, 1], func_loop, 1.e-4, 1.e-4, 50, data)
			if ps2[1] > peak:
				sx = -ps2[0][0]
				sy = -ps2[0][1]
				ang = ps2[3]
				peak = ps2[1]
	
			co =  cos(ang*pi/180.0)
			so = -sin(ang*pi/180.0)
			sxs = sx*co - sy*so
			sys = sx*so + sy*co
			return  ang, sxs, sys, is_mirror, peak
		else: print "Unknown optimization strategy!"
		
	elif interpolation_method == "linear":
		print "This interpolation method is not supported!"

	else: 
		print "Error: Unknown interpolation method!"


def ormy2(image, refim, crefim, xrng, yrng, step, mode, numr, cnx, cny, interpolation_method=None):
	"""
	Determine shift and rotation between image and reference image (refim)
	It should be noted that crefim should be preprocessed before using this	function.
	This function is mostly same as the ormy, the only difference is that it 
	does fine refinement in the end
	"""
	from math import pi, cos, sin, sqrt
	from fundamentals import fft, mirror
	from utilities import amoeba, model_circle
	from sys import exit
	from alignment import ang_n
	from random import random
	
	if interpolation_method==None:
		interpolation_method = interpolation_method_2D #use global setting
		
	maxrin = numr[-1]
	nx = image.get_xsize()

	if interpolation_method=="gridding":
		npad = 2
		N = nx*npad
		#support of the window
		K = 6
		alpha = 1.75
		r = nx/2
		v = K/2.0/N
		kb = Util.KaiserBessel(alpha, K, r, v, N)
		imali = image.FourInterpol(2*nx, 2*nx, 1, 0)
		params = {"filter_type":Processor.fourier_filter_types.KAISER_SINH_INVERSE, "alpha":alpha, "K":K, "r":r, "v":v, "N":N}
		del K, alpha, r, v
		q = Processor.EMFourierFilter(imali, params)
		imali = fft(q)
	
		line = EMData()
		line.set_size(maxrin, 1, 1)
		# do not pad
		npad = 1
		N = maxrin*npad
		# support of the window
		K = 6
		alpha = 1.75
		r = maxrin/2
		v = K/2.0/N
		kbline = Util.KaiserBessel(alpha, K, r, v, N)
		parline = {"filter_type":Processor.fourier_filter_types.KAISER_SINH_INVERSE, "alpha":alpha, "K":K, "r":r, "v":v, "N":N}
		del K, alpha, r, v
		data = []
		data.insert(1, kbline)
	elif interpolation_method=="quadratic":
		pass
	elif interpolation_method=="linear":
		pass
	else: 
		print "Error: Unknown interpolation method!"
		exit()

	peak = -1.0E23
	ky = int(2*yrng/step+0.5)/2
	kx = int(2*xrng/step+0.5)/2
	
	if interpolation_method=="gridding":
		for i in xrange(-ky, ky+1):
			iy = i*step
			for j in xrange(-kx, kx+1):
				ix = j*step
				cimage = Util.Polar2Dmi(imali, cnx+ix, cny+iy, numr, mode, kb)
				Util.Frngs(cimage, numr)
				qt = Util.Crosrng_msg(crefim, cimage, numr)
				
				# straight
				for i in xrange(0, maxrin): line[i] = qt[i, 0]					
				#  find the current maximum and its location
				ps = line.max_search()
				qn  = ps[1]
				jtot = ps[0]/2
				q = Processor.EMFourierFilter(line, parline)
				data.insert(0, q)
				ps = amoeba([jtot], [2.0], oned_search_func, 1.e-4, 1.e-4, 500, data)
				del data[0]
				jtot = ps[0][0]*2
				qn = ps[1]

				# mirror
				for i in xrange(0, maxrin): line[i] = qt[i, 1]
				#  find the current maximum and its location
				ps = line.max_search()
				qm = ps[1]
				mtot = ps[0]/2
				q = Processor.EMFourierFilter(line, parline)
				data.insert(0, q)
				ps = amoeba([mtot], [2.0], oned_search_func, 1.e-4, 1.e-4, 500, data)
				del data[0]
				mtot = ps[0][0]*2
				qm = ps[1]

				if qn >= peak or qm >= peak:
					if qn >= qm:
						ang = ang_n(jtot+1, mode, maxrin)
						sx = -ix
						sy = -iy
						peak = qn
						is_mirror = 0
					else:
						ang = ang_n(mtot+1, mode, maxrin)
						sx = -ix
						sy = -iy
						peak = qm
						is_mirror = 1
	elif interpolation_method=="quadratic":
		for i in xrange(-ky, ky+1):
			iy = i*step
			for j in xrange(-kx, kx+1):
				ix = j*step
				cimage = Util.Polar2Dm(image, cnx+ix, cny+iy, numr, mode)
				Util.Frngs(cimage, numr)
				retvals = Util.Crosrng_ms(crefim, cimage, numr)
				qn = retvals["qn"]
				qm = retvals["qm"]
				if qn >= peak or qm >= peak:
					sx = -ix
					sy = -iy
					if qn >= qm:
						ang = ang_n(retvals["tot"], mode, maxrin)
						peak = qn
						is_mirror = 0
					else:
						ang = ang_n(retvals["tmt"], mode, maxrin)
						peak = qm
						is_mirror = 1		
	
	co =  cos(ang*pi/180.0)
	so = -sin(ang*pi/180.0)
	sxs = sx*co - sy*so
	sys = sx*so + sy*co
	
	last_ring = numr[-3]
	mask = model_circle(last_ring, nx, nx)

	if interpolation_method=="gridding":
		if is_mirror:	
			image2 = mirror(image)
			imali = image2.FourInterpol(2*nx, 2*nx, 1, 0)
			q = Processor.EMFourierFilter(imali, params)
			imali = fft(q)
			ang = -ang
			sxs = -sxs
		data.insert(0, imali)
		data.insert(1, refim)
		data.insert(2, kb)
		data.insert(3, mask)
		ps = amoeba([ang, sxs, sys], [360.0/(2.0*pi*last_ring), step, step], dfunc, 1.e-4, 1.e-4, 500, data)
	elif interpolation_method=="quadratic":
		data = []
		if is_mirror:
			data.insert(0, mirror(image))
			ang = -ang
			sxs = -sxs
		else:
			data.insert(0, image)
		data.insert(1, refim)
		data.insert(2, mask)
		ps = amoeba([ang, sxs, sys], [360.0/(2.0*pi*last_ring), step, step], dfunc2, 1.e-4, 1.e-4, 500, data)
	else: pass
	
	ps[0][0] = ps[0][0]%360.0
	if is_mirror:
		ps[0][0] = -ps[0][0]
		ps[0][1] = -ps[0][1]
	
	return  ps[0][0], ps[0][1], ps[0][2], is_mirror, -ps[1]


def ormy2lbfgsb(image,refim,crefim,xrng,yrng,step,mode,numr,cnx,cny,interpolation_method=None):
	"""
	Determine shift and rotation between image and reference image (refim)
	It should be noted that crefim should be preprocessed before using this
	function.
	This function is mostly same as the ormy2, the only difference is that 
	it use LBFGSB for parameter optimization
	"""
	from math import pi, cos, sin
	from fundamentals import fft, mirror
	from utilities import amoeba, model_circle
	from sys import exit
	from alignment import ang_n
	from scipy.optimize.lbfgsb import fmin_l_bfgs_b
	
	if interpolation_method==None:
		interpolation_method = interpolation_method_2D #use global setting
		
	if interpolation_method=="gridding":
		M=image.get_xsize()
		npad=2
		N=M*npad
		#support of the window
		K=6
		alpha=1.75
		r=M/2
		v=K/2.0/N
		kb = Util.KaiserBessel(alpha, K, r, K/(2.*N), N)

		imali=image.FourInterpol(2*M,2*M,1,0)
		params = {"filter_type" : Processor.fourier_filter_types.KAISER_SINH_INVERSE,"alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
		q=Processor.EMFourierFilter(imali,params)
		imali=fft(q)
	
		maxrin=numr[len(numr)-1]
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
		data=[]
		data.insert(1,kbline)
	elif interpolation_method=="quadratic":
		pass
	elif interpolation_method=="linear":
		pass
	else: 
		print "Error: Unknown interpolation method!"
		exit()

	peak = -1.0E23
	ky = int(2*yrng/step+0.5)//2
	kx = int(2*xrng/step+0.5)//2
	
	if interpolation_method=="gridding":
		for i in xrange(-ky,ky+1):
			iy=i*step
			for  j in xrange(-kx,kx+1):
				ix=j*step
				cimage=Util.Polar2Dmi(imali,cnx+ix,cny+iy,numr,mode,kb)
				Util.Frngs(cimage, numr)
				qt = Util.Crosrng_msg(crefim,cimage,numr)
				
				# straight
				for i in xrange(0,maxrin): line[i]=qt[i,0]					
				#  find the current maximum and its location
				ps=line.peak_search(1,1)				
				qn=ps[1]
				jtot=ps[2]/2
				q=Processor.EMFourierFilter(line,parline)
				data.insert(0,q)
				ps = amoeba([jtot], [2.0], oned_search_func, 1.e-4,1.e-4,500,data)
				del data[0]
				jtot=ps[0][0]*2
				qn=ps[1]

				# mirror
				for i in xrange(0,maxrin): line[i]=qt[i,1]
				#  find the current maximum and its location
				ps=line.peak_search(1,1)				
				qm=ps[1]
				mtot=ps[2]/2
				q=Processor.EMFourierFilter(line,parline)
				data.insert(0,q)
				ps = amoeba([mtot], [2.0], oned_search_func,1.e-4,1.e-4,500,data)
				del data[0]
				mtot=ps[0][0]*2
				qm=ps[1]

				if(qn >=peak or qm >=peak):
					if (qn >= qm):
						ang = ang_n(jtot+1, mode, numr[len(numr)-1])
						sx = -ix
						sy = -iy
						peak=qn
						is_mirror=0
					else:
						ang = ang_n(mtot+1, mode, numr[len(numr)-1])
						sx = -ix
						sy = -iy
						peak=qm
						is_mirror=1
	elif interpolation_method=="quadratic":
		for i in xrange(-ky, ky+1):
			iy = i*step
			for  j in xrange(-kx, kx+1):
				ix = j*step
				cimage=Util.Polar2Dm(image, cnx+ix, cny+iy, numr, mode)
				Util.Frngs(cimage, numr)
				retvals = Util.Crosrng_ms(crefim, cimage, numr)
				qn = retvals["qn"]
				qm = retvals["qm"]
				if(qn >=peak or qm >=peak):
					sx = -ix
					sy = -iy
					if (qn >= qm):
						ang = ang_n(retvals["tot"], mode, numr[len(numr)-1])
						peak=qn
						is_mirror=0
					else:
						ang = ang_n(retvals["tmt"], mode, numr[len(numr)-1])
						peak=qm
						is_mirror=1		
	
	co =  cos(ang*pi/180.0)
	so = -sin(ang*pi/180.0)
	sxs = sx*co - sy*so
	sys = sx*so + sy*co	
	
	last_ring = numr[len(numr)-3]
	nx = image.get_xsize()
	ny = image.get_ysize()
	mask = model_circle(last_ring,nx,ny)

	if interpolation_method=="gridding":
		if is_mirror:	
			M=image.get_xsize()
			image2 = mirror(image)
			imali=image2.FourInterpol(2*M,2*M,1,0)
			params = {"filter_type" : Processor.fourier_filter_types.KAISER_SINH_INVERSE,"alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
			q=Processor.EMFourierFilter(imali,params)
			imali=fft(q)
			#ps = Util.twoD_fine_ali_G(imali, refim, mask, kb, -ang, -sxs, sys);
			x0 = [-ang, -sxs, sys] 
			bounds = [(x0[0]-2.0, x0[0]+2.0), (x0[1]-1.5, x0[1]+1.5), (x0[2]-1.5, x0[2]+1.5)]
			data = []
			data.append(imali)
			data.append(refim)
			data.append(kb)
			data.append(mask)
			ps, val, d = fmin_l_bfgs_b(dfunc_i, x0, args=[data], approx_grad=1, bounds=bounds, m=10, factr=1e1, pgtol=1e-5, epsilon=1e-3, iprint=-1, maxfun=15000)
		else:
			#ps = Util.twoD_fine_ali_G(imali, refim, mask, kb, ang, sxs, sys);			
			x0 = [ang, sxs, sys] 
			bounds = [(x0[0]-2.0, x0[0]+2.0), (x0[1]-1.5, x0[1]+1.5), (x0[2]-1.5, x0[2]+1.5)]
			data = []
			data.append(imali)
			data.append(refim)
			data.append(kb)
			data.append(mask)
			ps, val, d = fmin_l_bfgs_b(dfunc_i, x0, args=[data], approx_grad=1, bounds=bounds, m=10, factr=1e1, pgtol=1e-5, epsilon=1e-3, iprint=-1, maxfun=15000)
	elif interpolation_method=="quadratic":
		if is_mirror:
			image=mirror(image)
			#ps = Util.twoD_fine_ali(image, refim, mask, -ang, -sxs, sys)
			x0 = [-ang, -sxs, sys] 
			bounds = [(x0[0]-2.0, x0[0]+2.0), (x0[1]-1.5, x0[1]+1.5), (x0[2]-1.5, x0[2]+1.5)]
			data = []
			data.append(image)
			data.append(refim)
			data.append(mask)
			ps, val, d = fmin_l_bfgs_b(dfunc2_i, x0, args=[data], approx_grad=1, bounds=bounds, m=10, factr=1e1, pgtol=1e-5, epsilon=1e-3, iprint=-1, maxfun=15000)
		else:
			#ps = Util.twoD_fine_ali(image, refim, mask, ang, sxs, sys)
			x0 = [ang, sxs, sys] 
			bounds = [(x0[0]-2.0, x0[0]+2.0), (x0[1]-1.5, x0[1]+1.5), (x0[2]-1.5, x0[2]+1.5)]
			data = []
			data.append(image)
			data.append(refim)
			data.append(mask)			
			ps, val, d = fmin_l_bfgs_b(dfunc2_i, x0, args=[data], approx_grad=1, bounds=bounds, m=10, factr=1e1, pgtol=1e-5, epsilon=1e-3, iprint=-1, maxfun=15000)
	else: pass

	ps[0] = ps[0]%360.0
	if is_mirror:
		ps[0] = -ps[0]
		ps[1] = -ps[1]
	return  ps[0], ps[1], ps[2], is_mirror


def ormy2sd(image,refim,crefim,xrng,yrng,step,mode,numr,cnx,cny,interpolation_method=None):
	"""
	Determine shift and rotation between image and reference image (refim)
	It should be noted that crefim should be preprocessed before using this
	function.
	This function is mostly same as the ormy2, the only difference is that 
	it use steepest descent algorithm for parameter optimization
	"""
	from math import pi, cos, sin
	from fundamentals import fft, mirror
	from utilities import amoeba, model_circle
	from sys import exit
	from alignment import ang_n
	
	if interpolation_method==None:
		interpolation_method = interpolation_method_2D #use global setting
		
	if interpolation_method=="gridding":
		M=image.get_xsize()
		npad=2
		N=M*npad
		#support of the window
		K=6
		alpha=1.75
		r=M/2
		v=K/2.0/N
		kb = Util.KaiserBessel(alpha, K, r, K/(2.*N), N)

		imali=image.FourInterpol(2*M,2*M,1,0)
		params = {"filter_type" : Processor.fourier_filter_types.KAISER_SINH_INVERSE,"alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
		q=Processor.EMFourierFilter(imali,params)
		imali=fft(q)
	
		maxrin=numr[len(numr)-1]
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
		data=[]
		data.insert(1,kbline)
	elif interpolation_method=="quadratic":
		pass
	elif interpolation_method=="linear":
		pass
	else: 
		print "Error: Unknown interpolation method!"
		exit()

	peak = -1.0E23
	ky = int(2*yrng/step+0.5)//2
	kx = int(2*xrng/step+0.5)//2
	
	if interpolation_method=="gridding":
		for i in xrange(-ky,ky+1):
			iy=i*step
			for  j in xrange(-kx,kx+1):
				ix=j*step
				cimage=Util.Polar2Dmi(imali,cnx+ix,cny+iy,numr,mode,kb)
				Util.Frngs(cimage, numr)
				qt = Util.Crosrng_msg(crefim,cimage,numr)
				
				# straight
				for i in xrange(0,maxrin): line[i]=qt[i,0]					
				#  find the current maximum and its location
				ps=line.peak_search(1,1)				
				qn=ps[1]
				jtot=ps[2]/2
				q=Processor.EMFourierFilter(line,parline)
				data.insert(0,q)
				ps = amoeba([jtot], [2.0], oned_search_func, 1.e-4,1.e-4,500,data)
				del data[0]
				jtot=ps[0][0]*2
				qn=ps[1]

				# mirror
				for i in xrange(0,maxrin): line[i]=qt[i,1]
				#  find the current maximum and its location
				ps=line.peak_search(1,1)				
				qm=ps[1]
				mtot=ps[2]/2
				q=Processor.EMFourierFilter(line,parline)
				data.insert(0,q)
				ps = amoeba([mtot], [2.0], oned_search_func,1.e-4,1.e-4,500,data)
				del data[0]
				mtot=ps[0][0]*2
				qm=ps[1]

				if(qn >=peak or qm >=peak):
					if (qn >= qm):
						ang = ang_n(jtot+1, mode, numr[len(numr)-1])
						sx = -ix
						sy = -iy
						peak=qn
						is_mirror=0
					else:
						ang = ang_n(mtot+1, mode, numr[len(numr)-1])
						sx = -ix
						sy = -iy
						peak=qm
						is_mirror=1
	elif interpolation_method=="quadratic":
		for i in xrange(-ky, ky+1):
			iy = i*step
			for  j in xrange(-kx, kx+1):
				ix = j*step
				cimage=Util.Polar2Dm(image, cnx+ix, cny+iy, numr, mode)
				Util.Frngs(cimage, numr)
				retvals = Util.Crosrng_ms(crefim, cimage, numr)
				qn = retvals["qn"]
				qm = retvals["qm"]
				if(qn >=peak or qm >=peak):
					sx = -ix
					sy = -iy
					if (qn >= qm):
						ang = ang_n(retvals["tot"], mode, numr[len(numr)-1])
						peak=qn
						is_mirror=0
					else:
						ang = ang_n(retvals["tmt"], mode, numr[len(numr)-1])
						peak=qm
						is_mirror=1		
	
	co =  cos(ang*pi/180.0)
	so = -sin(ang*pi/180.0)
	sxs = sx*co - sy*so
	sys = sx*so + sy*co	
	
	last_ring = numr[len(numr)-3]
	nx = image.get_xsize()
	ny = image.get_ysize()
	mask = model_circle(last_ring,nx,ny)

	if interpolation_method=="gridding":
		if is_mirror:	
			M=image.get_xsize()
			image2 = mirror(image)
			imali=image2.FourInterpol(2*M,2*M,1,0)
			params = {"filter_type" : Processor.fourier_filter_types.KAISER_SINH_INVERSE,"alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
			q=Processor.EMFourierFilter(imali,params)
			imali=fft(q)
			ps = Util.twoD_fine_ali_SD_G(imali, refim, mask, kb, -ang, -sxs, sys);
		else:
			ps = Util.twoD_fine_ali_SD_G(imali, refim, mask, kb, ang, sxs, sys);			
	elif interpolation_method=="quadratic":
		if is_mirror:
			image=mirror(image)
			ps = Util.twoD_fine_ali_SD(image, refim, mask, -ang, -sxs, sys)
		else:
			ps = Util.twoD_fine_ali_SD(image, refim, mask, ang, sxs, sys)
	else: pass

	if ps[0] < 180: ps[0] = ps[0]+360
	if ps[0] > 180: ps[0] = ps[0]-360
	return  ps[0], ps[1], ps[2], is_mirror, ps[3]


def ormy3(image,crefim,xrng,yrng,step,mode,numr,cnx,cny,interpolation_method=None):
	"""
	Determine shift and rotation between image and reference image (refim)
	
	It should be noted that crefim should be preprocessed before using this	function.
	
	This function is mostly same as the ormy2, the difference is that is use a different 
	strategy for refinement. To be specific, after the coarse search, it refines [sx, sy] 
	instead of refining [ang, sx, sy]. 
	
	In the future, we can delete ormy and ormy2 if this version is proved to be stable.
	"""
	from math import pi, cos, sin
	from fundamentals import fft, mirror
	from utilities import amoeba_multi_level, model_circle, amoeba
	from sys import exit
	from alignment import ang_n
	
	if interpolation_method==None:
		interpolation_method = interpolation_method_2D #use global setting
		
	if interpolation_method=="gridding":
		M=image.get_xsize()
		npad=2
		N=M*npad
		#support of the window
		K=6
		alpha=1.75
		r=M/2
		v=K/2.0/N
		kb = Util.KaiserBessel(alpha, K, r, K/(2.*N), N)

		imali=image.FourInterpol(2*M,2*M,1,0)
		params = {"filter_type" : Processor.fourier_filter_types.KAISER_SINH_INVERSE,"alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
		q=Processor.EMFourierFilter(imali,params)
		imali=fft(q)
	
		maxrin=numr[len(numr)-1]
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
		data=[]
		data.insert(1,kbline)
	elif interpolation_method=="quadratic":
		pass
	elif interpolation_method=="linear":
		pass
	else: 
		print "Error: Unknown interpolation method!"
		exit()

	peak = -1.0E23
	ky = int(2*yrng/step+0.5)//2
	kx = int(2*xrng/step+0.5)//2

	if interpolation_method=="gridding":
		for i in xrange(-ky,ky+1):
			iy=i*step
			for  j in xrange(-kx,kx+1):
				ix=j*step
				cimage=Util.Polar2Dmi(imali,cnx+ix,cny+iy,numr,mode,kb)
				Util.Frngs(cimage, numr)
				
				# straight, find the current maximum and its location
				line_s = Util.Crosrng_msg_s(crefim,cimage,numr)
				ps=line_s.max_search()				
				qn=ps[1]
				jtot=ps[0]/2
				q=Processor.EMFourierFilter(line_s,parline)
				data.insert(0,q)
				ps = amoeba([jtot], [2.0], oned_search_func, 1.e-4,1.e-4,50,data)
				del data[0]
				jtot=ps[0][0]*2
				qn=ps[1]

				# mirror, find the current maximum and its location
				line_m = Util.Crosrng_msg_m(crefim,cimage,numr)
				ps=line_m.max_search()				
				qm=ps[1]
				mtot=ps[0]/2
				q=Processor.EMFourierFilter(line_m,parline)
				data.insert(0,q)
				ps = amoeba([mtot], [2.0], oned_search_func,1.e-4,1.e-4,50,data)
				del data[0]
				mtot=ps[0][0]*2
				qm=ps[1]

				if(qn >=peak or qm >=peak):
					if (qn >= qm):
						ang = ang_n(jtot+1, mode, numr[len(numr)-1])
						sx = -ix
						sy = -iy
						peak = qn
						is_mirror=0
					else:
						ang = ang_n(mtot+1, mode, numr[len(numr)-1])
						sx = -ix
						sy = -iy
						peak = qm
						is_mirror=1
		data0 = []
		data0.insert(0,imali)
		data0.insert(1,cnx)
		data0.insert(2,cny)
		data0.insert(3,numr)
		data0.insert(4,mode)
		data0.insert(5,kb)
		data0.insert(6,crefim)
		data0.insert(7,parline)
		data0.insert(8,maxrin)
		data0.insert(9,kbline)
		data0.insert(10,is_mirror)
		
		ps2 = amoeba_multi_level([-sx,-sy],[1,1],func_loop2,1.e-4,1.e-4,50,data0)
		
		if ps2[1] > peak:
			sx = -ps2[0][0]
			sy = -ps2[0][1]
			ang = ps2[3]
	elif interpolation_method=="quadratic":
		for i in xrange(-ky, ky+1):
			iy = i*step
			for  j in xrange(-kx, kx+1):
				ix = j*step
				cimage=Util.Polar2Dm(image, cnx+ix, cny+iy, numr, mode)
				Util.Frngs(cimage, numr)
				retvals = Util.Crosrng_ms(crefim, cimage, numr)
				qn = retvals["qn"]
				qm = retvals["qm"]
				if(qn >=peak or qm >=peak):
					sx = -ix
					sy = -iy
					if (qn >= qm):
						ang = ang_n(retvals["tot"], mode, numr[len(numr)-1])
						peak=qn
						is_mirror=0
					else:
						ang = ang_n(retvals["tmt"], mode, numr[len(numr)-1])
						peak=qm
						is_mirror=1	 
		
		data = []
		data.insert(0,image)
		data.insert(1,cnx)
		data.insert(2,cny)
		data.insert(3,numr)
		data.insert(4,mode)
		data.insert(5,crefim)
		ps2 = amoeba_multi_level([-sx,-sy],[1,1],func_loop,1.e-4,1.e-4,50,data)
		if ps2[1]>peak:
			sx = -ps2[0][0]
			sy = -ps2[0][1]
			ang = ps2[3]
		
	co =  cos(ang*pi/180.0)
	so = -sin(ang*pi/180.0)
	sxs = sx*co - sy*so
	sys = sx*so + sy*co
	return  ang,sxs,sys,is_mirror,peak	


def ormy3g(imali, kb, crefim, xrng, yrng, step, mode, numr, cnx, cny):
	"""
	Determine shift and rotation between image and reference image (refim)
	
	It should be noted that crefim should be preprocessed before using this	function.
	
	This function is mostly same as the ormy2, the difference is that is use a different 
	strategy for refinement. To be specific, after the coarse search, it refines [sx, sy] 
	instead of refining [ang, sx, sy]. 
	
	In the future, we can delete ormy and ormy2 if this version is proved to be stable.
	
	This version is only used for gridding interpolation and assume imali is alread prepared, which
	makes it suitable for alignment of a set of images.
	Also, according to my experiment, ormy2 is not good for alignment of a set of images.	
	"""
	from math import pi, cos, sin
	from fundamentals import fft, mirror
	from utilities import amoeba, model_circle, amoeba_multi_level
	from alignment import ang_n
	
	maxrin = numr[-1]
	# do not pad
	npad = 1
	N = maxrin*npad
	# support of the window
	K = 6
	alpha = 1.75
	r = maxrin/2
	v = K/2.0/N
	kbline = Util.KaiserBessel(alpha, K, r, v, N)
	parline = {"filter_type":Processor.fourier_filter_types.KAISER_SINH_INVERSE, "alpha":alpha, "K":K, "r":r, "v":v, "N":N}
	del K, r, v, N
	data = []
	data.insert(1, kbline)

	peak = -1.0E23
	ky = int(2*yrng/step+0.5)//2
	kx = int(2*xrng/step+0.5)//2

	for i in xrange(-ky, ky+1):
		iy = i*step
		for j in xrange(-kx, kx+1):
			ix = j*step
			cimage = Util.Polar2Dmi(imali, cnx+ix, cny+iy, numr, mode, kb)
			Util.Frngs(cimage, numr)
			
			# straight, find the current maximum and its location
			line_s = Util.Crosrng_msg_s(crefim, cimage, numr)
			ps = line_s.max_search()	
			qn = ps[1]
			jtot = ps[0]/2
			q = Processor.EMFourierFilter(line_s, parline)
			data.insert(0, q)
			ps = amoeba([jtot], [2.0], oned_search_func, 1.e-4, 1.e-4, 50, data)
			del data[0]
			jtot = ps[0][0]*2
			qn = ps[1]

			# mirror, find the current maximum and its location
			line_m = Util.Crosrng_msg_m(crefim, cimage, numr)
			ps = line_m.max_search()				
			qm = ps[1]
			mtot = ps[0]/2
			q = Processor.EMFourierFilter(line_m, parline)
			data.insert(0, q)
			ps = amoeba([mtot], [2.0], oned_search_func, 1.e-4, 1.e-4, 50, data)
			del data[0]
			mtot = ps[0][0]*2
			qm = ps[1]

			if qn >= peak or qm >= peak:
				if qn >= qm:
					ang = ang_n(jtot+1, mode, maxrin)
					sx = -ix
					sy = -iy
					peak = qn
					is_mirror = 0
				else:
					ang = ang_n(mtot+1, mode, maxrin)
					sx = -ix
					sy = -iy
					peak = qm
					is_mirror = 1

	data0 = []
	data0.insert(0, imali)
	data0.insert(1, cnx)
	data0.insert(2, cny)
	data0.insert(3, numr)
	data0.insert(4, mode)
	data0.insert(5, kb)
	data0.insert(6, crefim)
	data0.insert(7, parline)
	data0.insert(8, maxrin)
	data0.insert(9, kbline)
	data0.insert(10, is_mirror)
	
	ps2 = amoeba_multi_level([-sx, -sy], [1.0, 1.0], func_loop2, 1.e-4, 1.e-4, 50, data0)
	if ps2[1] > peak and abs(ps2[0][0])<xrng and abs(ps2[0][1])<yrng:
		sx = -ps2[0][0]
		sy = -ps2[0][1]
		ang = ps2[3]
				
	co =  cos(ang*pi/180.0)
	so = -sin(ang*pi/180.0)
	sxs = sx*co - sy*so
	sys = sx*so + sy*co
	return ang, sxs, sys, is_mirror, peak	


def func_loop(args, data):
	from alignment import ang_n

	image = data[0]
	cnx = data[1]
	cny = data[2]
	numr = data[3]
	mode = data[4]
	crefim = data[5]
	ix = args[0]
	iy = args[1]
	cimage = Util.Polar2Dm(image, cnx+ix, cny+iy, numr, mode)
	Util.Frngs(cimage, numr)
	retvals = Util.Crosrng_ms(crefim, cimage, numr)
	qn = retvals["qn"]
	qm = retvals["qm"]
	if qn >= qm:
		ang = ang_n(retvals["tot"], mode, numr[-1])
		return qn, ang
	else :
		ang = ang_n(retvals["tmt"], mode, numr[-1])
		return qm, ang

		
def func_loop2(args, data):
	from utilities import amoeba
	from alignment import ang_n
	
	imali = data[0]
	cnx = data[1]
	cny = data[2]
	numr = data[3]
	mode = data[4]
	kb = data[5]
	crefim = data[6]
	parline = data[7]
	maxrin = data[8]
	kbline = data[9]
	is_mirror = data[10]
	
	ix = args[0]
	iy = args[1]
	
	cimage = Util.Polar2Dmi(imali, cnx+ix, cny+iy, numr, mode, kb)
	Util.Frngs(cimage, numr)
	
	if is_mirror == 0:
		line = Util.Crosrng_msg_s(crefim, cimage, numr)
		#  Straight, find the current maximum and its location
		ps = line.max_search()
		jtot = ps[0]/2
		q = Processor.EMFourierFilter(line, parline)
		data1 = [q, kbline]
		ps = amoeba([jtot], [2.0], oned_search_func, 1.e-4, 1.e-4, 50, data1)
		jtot = ps[0][0]*2
		ang = ang_n(jtot+1, mode, maxrin)
		return ps[1], ang
	else:
		line = Util.Crosrng_msg_m(crefim, cimage, numr)
		#  Mirror, find the current maximum and its location
		ps = line.max_search()
		mtot = ps[0]/2
		q = Processor.EMFourierFilter(line, parline)
		data1 = [q, kbline]
		ps = amoeba([mtot], [2.0], oned_search_func, 1.e-4, 1.e-4, 50, data1)
		mtot = ps[0][0]*2
		ang = ang_n(mtot+1, mode, maxrin)
		return ps[1], ang

'''
def proj_ali(volref, mask, projdata, first_ring, last_ring, rstep, xrng, yrng, step, dtheta, ref_a, sym, CTF = False):
	from utilities    import even_angles, model_circle, drop_spider_doc, info
	from projection   import prgs
	from fundamentals import prepg
	from statistics   import ttime
	from math import cos
	from string import split

	print  ttime()
	mode    = "F"
	# generate list of Eulerian angles for reference projections
	#  phi, theta, psi
	tmp_sym = split(sym)
	if(tmp_sym[0][0]=="d" or tmp_sym[0][0]=="D"): nfold = int(tmp_sym[0][1])*2
	else:                                       nfold = int(tmp_sym[0][1])
	ref_angles = even_angles(dtheta, 0.0, 90.0, 0.0, 359.99/nfold, ref_a, "Zero")
	#drop_spider_doc("angles.txt",ref_angles)
	num_ref = len(ref_angles)
	# generate reference projections in polar coords
	#  begin from applying the mask, i.e., subtract the average under the mask and multiply by the mask
	#print "volref"
	#info(volref)
	if(mask):
		[mean, sigma, xmin, xmax ] =  Util.infomask(volref, mask)

		volref -= mean
		Util.mul_img(volref, mask)

	nx=volref.get_xsize()
	ny=volref.get_ysize()
	#  center is in SPIDER convention
	cnx = int(nx/2)+1
	cny = int(ny/2)+1
	#precalculate rings
	numr=Numrinit(first_ring, last_ring, rstep, mode)
	wr = ringwe(numr ,mode)
	# K-B
	kb=kbt(nx)
	# prepare the reference volume
	volrefft = prep_vol_kb(volref,kb)
	print  ttime()

	# prepare 2-D ask for normalization
	mask2D = model_circle(last_ring, nx, ny)
	if(first_ring > 0):
		prjref = model_circle(first_ring, nx, ny)
		mask2D -= prjref

	ref_proj_rings = []     # list of (image objects) reference projections in Fourier representation
	# loop over reference angles
	for i in xrange(num_ref):
		params = [ref_angles[i][0],ref_angles[i][1],ref_angles[i][2],0.0,0.0]
		#print  params
		prjref = prgs(volrefft, kb, params)

		prjref.process_inplace("normalize.mask", {"mask":mask2D, "no_sigma":1})  # (i-a)/s
		#prjref.write_image("prj.spi",i)
		#info(prjref)
		#temp   = prepg(prjref, kb)
		#cimage = Util.Polar2Dmi(temp, cnx, cny, numr, mode, kb)
		cimage=Util.Polar2Dm(prjref, cnx, cny, numr, mode)
		Util.Frngs(cimage, numr)  # currently set to gridding....
		Applyws(cimage, numr, wr)
		ref_proj_rings.insert(i, cimage)
		#info(ref_proj_rings[i])

	exptimage = EMData()
	#  find out how many images there are in the stack
	nprojdata = EMUtil.get_image_count(projdata)
	for imn in xrange(nprojdata):
		print  imn, ttime()
		exptimage.read_image(projdata, imn)
		if(CTF):
			ctf_params = e.get_attr('ctf')
			if e.get_attr('ctf_applied') == 0:
				from filter import filt_ctf
				exptimage = filt_ctf(exptimage, ctf_params)
				e.set_attr('ctf_applied', 1)
		#pg = prepg(exptimage, kb)
		#phi1 = exptimage.get_attr('phi')
		#theta1 = exptimage.get_attr('theta')
		#psi1 = exptimage.get_attr('psi')
		sxo = exptimage.get_attr('s2x')
		syo = exptimage.get_attr('s2y')
		peak = -1.0E23
		for i in xrange(num_ref):
			#print  "REFERENCE IMAGE:  ",i,"  ",ttime()
			#info(ref_proj_rings[i])
			[angt,sxst,syst,mirrort,peakt] = ormq(exptimage, ref_proj_rings[i], xrng, yrng, step, mode, numr, cnx-sxo, cny-syo)
			#print angt,sxst,syst,mirrort,peakt

			#print  "RESULT:  ",i, angt, sxst, syst, mirrort, peakt
			if(peakt>peak):
				peak   = peakt
				numref = i
				ang    = -angt
				sxs    = sxst
				sys    = syst
				mirror = mirrort
		ang = (ang+360.0)%360.0
		if  mirror:
			# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
			#  What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
			#print  imn," M ",ang,numref,ref_angles[numref],-sxs+sxo,-sys+syo
			exptimage.set_attr_dict({'phi':(ref_angles[numref][0]+180.0)%360.0, 'theta':180-ref_angles[numref][1], 'psi':(180+ref_angles[numref][2]-ang+360.0)%360.0, 'sx':-sxs+sxo, 'sy':-sys+syo})
		else:
			#print  imn,"   ",ang,numref,ref_angles[numref],sxs+sxo,sys+syo
			exptimage.set_attr_dict({'phi':ref_angles[numref][0], 'theta':ref_angles[numref][1], 'psi':(ref_angles[numref][2]+ang+360.0)%360.0, 's2x':sxs+sxo, 's2y':sys+syo})
		exptimage.write_image(projdata, imn, EMUtil.ImageType.IMAGE_HDF, True)
'''



'''
#  This version has gridding projections
def proj_ali_incore(volref, mask3D, projdata, first_ring, last_ring, rstep, xrng, yrng, step, delta, ref_a, symmetry, CTF = False):
	from alignment    import apmq
	from utilities    import even_angles, model_circle, compose_transform2
	from projection   import prgs
	from fundamentals import prepg

	#from sys import exit
	from utilities    import drop_image,info
	from statistics import memory
	initmem = memory()
	print  initmem/1.0e6
	mode    = "F"
	# generate list of Eulerian angles for reference projections
	#  phi, theta, psi
	ref_angles = even_angles(delta, symmetry = symmetry, method = ref_a, phiEqpsi = "Zero")
	#from utilities import drop_spider_doc
	#drop_spider_doc("angles.txt",ref_angles)
	num_ref = len(ref_angles)
	# generate reference projections in polar coords
	#  begin from applying the mask, i.e., subtract the average under the mask and multiply by the mask
	if(mask3D):
		[mean, sigma, xmin, xmax ] =  Util.infomask(volref, mask3D, True)
		volref -= mean
		Util.mul_img(volref, mask3D)
	#drop_image(volref, "volref.spi", "s")
	#exit()
	nx   = volref.get_xsize()
	ny   = volref.get_ysize()
	#  center is in SPIDER convention
	cnx  = int(nx/2)+1
	cny  = int(ny/2)+1
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr   = ringwe(numr ,mode)
	# K-B
	kb   = kbt(nx)
	# prepare the reference volume
	volrefft = prep_vol_kb(volref,kb)
	#print  ttime()

	# prepare 2-D ask for normalization
	mask2D = model_circle(last_ring, nx, ny)
	if(first_ring > 0):
		prjref = model_circle(first_ring, nx, ny)
		mask2D -= prjref

	print  " pro_ali before refproj ",memory(initmem)/1.0e6,num_ref
	ref_proj_rings = []     # list of (image objects) reference projections in Fourier representation
	# loop over reference angles
	for i in xrange(num_ref):
		params = [ref_angles[i][0], ref_angles[i][1], ref_angles[i][2], 0.0,0.0]
		#print  params
		prjref = prgs(volrefft, kb, params)
		#from utilities import get_im
		#prjref = get_im('projspi.spi',i)
		prjref.process_inplace("normalize.mask", {"mask":mask2D, "no_sigma":1})  # (i-a)/s
		#prjref.write_image("prjsparx_n.spi",i)
		#info(prjref)
		#temp   = prepg(prjref, kb)
		#cimage = Util.Polar2Dmi(temp, cnx, cny, numr, mode, kb)
		cimage=Util.Polar2Dm(prjref, cnx, cny, numr, mode)
		Util.Frngs(cimage, numr)  # currently set to quadratic/gridding....
		Applyws(cimage, numr, wr)
		ref_proj_rings.append(cimage)
		if(i == 0):  info(cimage)
		#info(ref_proj_rings[i])
	print  " pro_ali after refproj ",memory(initmem)/1.0e6
	#exit()
	#from utilities  import ttime
	#print  ttime()
	#soto = []
	for imn in xrange(len(projdata)):
		#pg = prepg(exptimage, kb)
		sxo = projdata[imn].get_attr('s2x')
		syo = projdata[imn].get_attr('s2y')
		#peak = -1.0E23
		# new
		[ang, sxs, sys, mirror, nref, peak] = Util.multiref_polar_ali_2d(projdata[imn], ref_proj_rings, xrng, yrng, step, mode, numr, cnx-sxo, cny-syo)
		numref=int(nref)
		#[ang,sxs,sys,mirror,peak,numref] = apmq(projdata[imn], ref_proj_rings, xrng, yrng, step, mode, numr, cnx-sxo, cny-syo)
		#ang = (ang+360.0)%360.0
		# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
		#  What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
		angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1,  -ang, 0.,0.,1)
		if  mirror:
			#print  imn," M ",ang,numref,ref_angles[numref],-sxs+sxo,-sys+syo
			projdata[imn].set_attr_dict({'phi':(ref_angles[numref][0]+540.0)%360.0, 'theta':180.0-ref_angles[numref][1], 'psi':(540.0-ref_angles[numref][2]+angb)%360.0, 's2x':sxb+sxo, 's2y':syb+syo})
		else:
			#print  imn,"   ",ang,numref,ref_angles[numref],sxs+sxo,sys+syo
			projdata[imn].set_attr_dict({'phi':ref_angles[numref][0], 'theta':ref_angles[numref][1], 'psi':(ref_angles[numref][2]+angb+360.0)%360.0, 's2x':sxb+sxo, 's2y':syb+syo})
	print  " pro_ali after apmq ",memory(initmem)/1.0e6
	#print  ttime()
	#from sys import exit
	#exit()
'''

def ref_3d_ali(image, mask, volref, currentparams):
	from projection import prep_vol
	from utilities  import amoeba
	# currentparams: 0 - phi, 1 - theta, 2 - psi, 3 - sx, 4 - sy.
	data=[]
	volft,kb=prep_vol(volref)
	data.insert(0,volft)
	data.insert(1,kb)
	data.insert(2,image)
	data.insert(3,mask)
	#print  eqproj(currentparams,data)
	#return amoeba(currentparams, [1.0,1.0,1.0,0.5,0.5],eqproj,1.e-4,1.e-4,500,data)
	return amoeba(currentparams, [2.0,2.0,2.0,1.0,1.0],eqproj,1.e-4,1.e-4,500,data)

def ref_ali(volref, mask3D, mask, stack_name, list_of_particles, trans_params, angles, delta):
	from projection import prep_vol
	from utilities  import amoeba,info
	from math       import pi
	from time       import time
	# currentparams: 0 - phi, 1 - theta, 2 - psi, 3 - sx, 4 - sy, mirror, peak
	# 
	#[mean, sigma, xmin, xmax, nx, ny, nz ] =  info(volref,mask3D)
	#volref=(volref-mean)*mask3D

	data=[]
	volft,kb=prep_vol(volref)
	data.insert(0,volft)
	data.insert(1,kb)
	data.insert(3,mask)
	image=EMData()
	newp = []
	for i in xrange (0,len(list_of_particles)):
		imn=list_of_particles[i]
		#print  "  Stack: ",stack_name,"  image number =",imn," time: ",time()
		image.read_image(stack_name,imn)
		#info(image)
		#temp=aptrans(image,trans_params[4*imn:4*imn+4])
		#temp=Util.window(image,81,81)
		data.insert(2,image)
		#info(temp)
		#data.insert(2,aptrans(image,trans_params[4*imn:4*imn+4]))
		
		# atparams: 0 - phi, 1 - theta, 2 - psi, 3 - sx, 4 - sy, mirror, peak
		atparams = [angles[imn][0],angles[imn][1],angles[imn][2],0.0,0.0]
		initial = eqproj(atparams,data) # DDDDDDDDDDDDDDDDDDD
		stheta = abs(sin(atparams[1]*pi/180.0))
		if(stheta >0.0):
			weight_phi=min(180.0,delta/stheta/2.0)
		else:
			weight_phi=180.0
		#print  delta,weight_phi,time()
		#print  atparams,eqproj(atparams,data)
		params =  amoeba(atparams, [weight_phi,delta,weight_phi,1.0,1.0],eqproj,1.e-4,1.e-4,500,data)
		del  data[2]
		print "ref_ali -- imn:", imn, "initial:", initial, "params[1]:", params[1]  # DDDDDDDDDDDDDDDDDDD
		newp.append([params[0][0],params[0][1],params[0][2],params[0][3],params[0][4]])
	return  newp

def ref_3d_ali(image,mask,volref,currentparams):
	from projection import prep_vol
	from utilities import amoeba
	data=[]
	volft,kb=prep_vol(volref)
	data.insert(0,volft)
	data.insert(1,kb)
	data.insert(2,image)
	data.insert(3,mask)
	print  eqproj(currentparams,data)
	#return amoeba(currentparams, [1.0,1.0,1.0,0.5,0.5],eqproj,1.e-4,1.e-4,500,data)
	return amoeba(currentparams, [2.0,2.0,2.0,1.0,1.0],eqproj,1.e-4,1.e-4,500,data)


# this is a new version of ref_ali where params are in the header, no CTF for the time being
def ref_alis(volref, mask3D, mask, stack_name, list_of_particles, delta):
	from projection import prep_vol
	from utilities import amoeba,info
	from math import pi
	from time import time
	# 
	#[mean, sigma, xmin, xmax, nx, ny, nz ] =  info(volref,mask3D)
	#volref=(volref-mean)*mask3D

	data=[]
	volft,kb=prep_vol(volref)
	data.insert(0,volft)
	data.insert(1,kb)
	data.insert(3,mask)
	image=EMData()
	for i in xrange (0,len(list_of_particles)):
		imn=list_of_particles[i]
		#print  "  Stack: ",stack_name,"  image number =",imn," time: ",time()
		image.read_image(stack_name,imn)
		phi   = image.get_attr('phi')
		theta = image.get_attr('theta')
		psi   = image.get_attr('psi')
		#  have to invert the sign, as now the shift is applied to the projection, not to the image
		sx = -image.get_attr('sx')
		sy = -image.get_attr('sy')
		data.insert(2,image)
		# atparams: 0 - phi, 1 - theta, 2 - psi, 3 - sx, 4 - sy, mirror, peak
		atparams = [phi,theta,psi,sx,sy]
		print atparams
		initial = eqproj(atparams,data) # DDDDDDDDDDDDDDDDDDD
		stheta = abs(sin(atparams[1]*pi/180.0))
		if(stheta >0.0):
			weight_phi=min(180.0,delta/stheta/2.0)
		else:
			weight_phi=180.0
		#print  delta,weight_phi,time()
		#print  atparams,eqproj(atparams,data)
		params =  amoeba(atparams, [weight_phi,delta,weight_phi,1.0,1.0],eqproj,1.e-4,1.e-4,500,data)
		print params[0][0],params[0][1],params[0][2],params[0][3],params[0][4]
		del  data[2]
		print "ref_ali -- imn:", imn, "initial:", initial, "params[1]:", params[1],params[2]  # DDDDDDDDDDDDDDDDDDD
		image.set_attr_dict({'phi':params[0][0], 'theta':params[0][1], 'psi':params[0][2], 'sx':-params[0][3], 'sy':-params[0][4]})
		image.write_image(stack_name, i, EMUtil.ImageType.IMAGE_HDF, True)
	return


def ref_3d_ali_xyz(image,mask,volref,currentparams):
	from projection import prep_vol
	from utilities import amoeba
	EULER_XYZ=Transform3D.EulerType.XYZ
	EULER_SPIDER = Transform3D.EulerType.SPIDER
	ss=dict()
	ss["phi"]=currentparams[0]
	ss["theta"]=currentparams[1]
	ss["psi"]=currentparams[2]
	Rss=Transform3D(EULER_SPIDER,ss)
	vv=Rss.get_rotation(EULER_XYZ)
	params=5*[0.0]
	params[0]=vv["xtilt"]
	params[1]=vv["ytilt"]
	params[2]=vv["ztilt"]
	params[3]=currentparams[3]
	params[4]=currentparams[4]
	#print  params
	data=[]
	volft,kb=prep_vol(volref)
	data.insert(0,volft)
	data.insert(1,kb)
	data.insert(2,image)
	data.insert(3,mask)
	#print  eqproj_xyz(params,data)
	#return amoeba(currentparams, [1.0,1.0,1.0,0.5,0.5],eqproj,1.e-4,1.e-4,500,data)
	args=amoeba(params, [4.0,4.0,4.0,1.0,1.0],eqproj_xyz,1.e-4,1.e-4,500,data)
	dict_xyz = dict()
	dict_xyz["xtilt"]=args[0][0]
	dict_xyz["ytilt"]=args[0][1]
	dict_xyz["ztilt"]=args[0][2]
	Rxyz=Transform3D(EULER_XYZ,dict_xyz)
	vv=Rxyz.get_rotation(EULER_SPIDER)
	#print  vv
	return  args

def ref_3d_ali_spin(image,mask,volref,currentparams):
    from projection import prep_vol
    from utilities import amoeba
    #print  currentparams
    EULER_SPIN=Transform3D.EulerType.SPIN
    EULER_SPIDER = Transform3D.EulerType.SPIDER
    ss=dict()
    ss["phi"]   = currentparams[0]
    ss["theta"] = currentparams[1]
    ss["psi"]   = currentparams[2]
    Rss=Transform3D(EULER_SPIDER,ss)
    vv=Rss.get_rotation(EULER_SPIN)
    #print  vv
    params=5*[0.0]
    from math import acos,pi
    thetan=acos(vv["n3"])*180.0/pi
    sinn=sin(thetan*pi/180.0)
    phin=acos(vv["n1"]/sinn)*180.0/pi
    params[0]=vv["Omega"]
    params[1]=thetan  #theta
    params[2]=phin  #phi
    params[3]=currentparams[3]
    params[4]=currentparams[4]
    #print  params


    dict_spin = dict()
    dict_spin["Omega"]=params[0]
    dict_spin["n1"]=sin(params[1]*pi/180.)*cos(params[2]*pi/180.)
    dict_spin["n2"]=sin(params[1]*pi/180.)*sin(params[2]*pi/180.)
    dict_spin["n3"]=cos(params[1]*pi/180.)
    Rspin=Transform3D(EULER_SPIN,dict_spin)
    ch=Rspin.get_rotation(EULER_SPIDER)
    #print  ch

    data=[]
    volft,kb=prep_vol(volref)
    data.insert(0,volft)
    data.insert(1,kb)
    data.insert(2,image)
    data.insert(3,mask)
    #print  eqproj_spin(params,data)
    #return amoeba(currentparams, [1.0,1.0,1.0,0.5,0.5],eqproj,1.e-4,1.e-4,500,data)
    args=amoeba(params, [4.0,4.0,4.0,1.0,1.0],eqproj_spin,1.e-4,1.e-4,500,data)
    #print  args
    dict_spin = dict()
    dict_spin["Omega"]=args[0][0]
    dict_spin["n1"]=sin(args[0][1]*pi/180.)*cos(args[0][2]*pi/180.)
    dict_spin["n2"]=sin(args[0][1]*pi/180.)*sin(args[0][2]*pi/180.)
    dict_spin["n3"]=cos(args[0][1]*pi/180.)
    Rspin=Transform3D(EULER_SPIN,dict_spin)
    vv=Rspin.get_rotation(EULER_SPIDER)
    params[0]=vv["phi"]
    params[1]=vv["theta"]
    params[2]=vv["psi"]
    params[3]=args[0][3]
    params[4]=args[0][4]
    #print  vv
    return  args  #,params


def rot_func(args,data):
	
	from utilities import drop_image,info

	'''The shift(float value) is divided into translational movement and an integer shift, dx & sx_ respectively'''
	sx_ = int(args[3]);sy_ = int(args[4]);sz_ = int(args[4]); #integer shift
	dx = args[3]-sx_;dy = args[4]-sy_;dz = args[5]-sz_;       #decimal shift/translational motion
	'''----------------------------------------------'''
	x=data[0].copy()
	x.rotate_translate(args[0],args[1],args[2],dx,dy,dz)
	#data[0]=x.copy()
	#data[0].rotate(args[0],args[1],args[2])
	print data[0].get_xsize(),data[0].get_ysize(),data[0].get_zsize(),sx_,sy_,sz_,args[0],args[1],args[2]
	temp_ref_image = Util.window(data[1],data[0].get_xsize(),data[0].get_ysize(),data[0].get_zsize(),sx_,sy_,sz_)
	rot_shift3D(image_in, phi, theta=0, psi=0, sx=0, sy=0,sz=0, scale=1, method='linear')
	res = -x.cmp("SqEuclidean",temp_ref_image,{"mask":data[2]})
	print res
	info(temp_ref_image)
	info(data[0])
	drop_image(x,"out.tfc")
	del temp_ref_image
	return res

def rot_align_func(inp_image,ref_image,mask_image,az,alt,phi,sx,sy,sz):

	from utilities import amoeba

	args=[az,alt,phi,sx,sy,sz]
	scale = [1.349, -0.349, 1.349, 1.1, -1.1, 1.1]
	data=[inp_image,ref_image,mask_image]
	result=amoeba(args, scale, rot_func, 1.e-4, 1.e-4, 500 ,data)
	print result
	#info(result[1])
	del data[:],args[:],scale[:]
	return result

"""
#  written in C
def enerp(ave, numr):
	'''
		Squared norm of FTs of rings
	'''
	from math import pi
	#from utilities  import print_col
	#print  "  ENER "
	#print_col(ave)
	pi2 = pi*2
	nring = len(numr)/3
	maxrin = numr[len(numr)-1]
	ener = 0.0
	#print  numr
	for i in xrange(nring):
		numr3i = numr[2+i*3]
		np = numr[1+i*3]-1
		tq = 2*pi*numr[i*3]/numr3i
		en = tq*(ave[np]**2+ave[np+1]**2)*0.5
		for j in xrange(np+2, np+numr3i-1):  en += tq*ave[j]**2
		ener += en/numr3i
	#print  ener, ave.cmp("dot", ave, {"negative":0})
	return  ener
"""


'''
def proj_ali_old(volref, mask3D, projdata, first_ring, last_ring, rstep, xrng, yrng, step, dtheta, ref_a):
	from utilities import even_angles, info
	from projection import prgs
	from fundamentals import prepg
	from statistics import ttime
	print  ttime()
	mode="F"
	# generate list of Eulerian angles for reference projections
	#  phi, theta, psi
	ref_angles = even_angles(dtheta,0,90,0,359.99,ref_a)  # add symmetries - PHIL
	num_ref=len(ref_angles)
	# generate reference projections in polar coords
	#  begin from applying the mask, i.e., subtract the average under the mask and multiply by the mask
	if(mask3D):
		[mean, sigma, xmin, xmax, nx, ny, nz ] =  info(volref,mask3D)
		volref=(volref-mean)*mask3D

	nx=volref.get_xsize()
	ny=volref.get_ysize()
	#  center is in SPIDER convention
	cnx = int(nx/2)+1
	cny = int(ny/2)+1
	#precalculate rings
	numr=Numrinit(first_ring, last_ring, rstep, mode)
	wr=ringwe(numr ,mode)
	# K-B
	kb=kbt(nx)
	# prepare the reference volume
	volrefft = prep_vol_kb(volref,kb)
	print  ttime()

	ref_proj_rings = []     # list of (image objects) reference projections in Fourier representation
	# loop over reference angles
	for i in xrange(num_ref):
		params = [ref_angles[i][0],ref_angles[i][1],ref_angles[i][2],0.0,0.0]
		#print  params
		prjref = prgs(volrefft, kb, params)
		#info(prjref)
		temp = prepg(prjref, kb)
		cimage=Util.Polar2Dmi(temp, cnx, cny, numr, mode, kb)
		#cimage=Util.Polar2Dm(prjref,cnx,cny,numr,mode)
		Util.Frngs(cimage, numr)
		Applyws(cimage, numr, wr)
		ref_proj_rings.insert(i,cimage)
		#info(ref_proj_rings[i])
	 
	exptimage = EMData()
	#  find out how many images there are in the stack
	nprojdata = EMUtil.get_image_count(projdata)
	for imn in xrange(nprojdata):
		print  imn, ttime()
		exptimage.read_image(projdata, imn)
		pg = prepg(exptimage, kb)
		#phi1 = exptimage.get_attr('phi')
		#theta1 = exptimage.get_attr('theta')
		#psi1 = exptimage.get_attr('psi')
		sxo = exptimage.get_attr('sx')
		syo = exptimage.get_attr('sy')
		peak = -1.0E23
		for i in xrange(num_ref):
			#print  "REFERENCE IMAGE:  ",i,"  ",ttime()
			#info(ref_proj_rings[i])
			[angt,sxst,syst,mirrort,peakt]=ormqip(pg,ref_proj_rings[i],xrng,yrng,step,mode,numr,kb,cnx-sxo,cny-syo,nx)
			#print  "RESULT:  ",i,angt,sxst,syst,mirrort,peakt
			if(peakt>peak):
				peak=peakt
				numref=i
				ang=-angt
				sxs=sxst
				sys=syst
				mirror=mirrort
		#print  imn,newangles,sxs,sys,mirror,peak
		ang = (ang+360.0)%360.0
		if  mirror:
			# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
			#  What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
			print  imn," M ",ang,numref,ref_angles[numref],-sxs+sxo,-sys+syo
			exptimage.set_attr_dict({'phi':(ref_angles[numref][0]+180.0)%360.0, 'theta':180-ref_angles[numref][1], 'psi':(180+ref_angles[numref][2]-ang+360.0)%360.0, 'sx':-sxs+sxo, 'sy':-sys+syo})
		else:
			print  imn,"   ",ang,numref,ref_angles[numref],sxs+sxo,sys+syo
			exptimage.set_attr_dict({'phi':ref_angles[numref][0], 'theta':ref_angles[numref][1], 'psi':(ref_angles[numref][2]+ang+360.0)%360.0, 'sx':sxs+sxo, 'sy':sys+syo})
		exptimage.write_image(projdata, imn, EMUtil.ImageType.IMAGE_HDF, True)

	return
'''


'''
def rec3D_MPI_noCTF_index(data, index, symmetry, myid, main_node = 0):
'''
  #This function is to be called within an MPI program to do a reconstruction on a dataset kept in the memory 
  #Computes reconstruction from projections whose attribute 'group' = index...
'''
	from utilities  import model_blank, reduce_EMData_to_root, get_image,send_EMData, recv_EMData
	from random import randint

	Ttype = Transform3D.EulerType.SPIDER
	fftvol = EMData()
	weight = EMData()
	nx = data[0].get_xsize()
	rec_params = {"size":nx, "npad":4, "symmetry":symmetry, "fftvol":fftvol, "weight":weight}
	rec = Reconstructors.get( "nn4", rec_params )
	rec.setup()

	for i in xrange(len(data)):
		group = data[i].get_attr('group')
		if(group == index):
			phi   = data[i].get_attr('phi')
			theta = data[i].get_attr('theta')
			psi   = data[i].get_attr('psi')
			rec.insert_slice(data[i], Transform3D(Ttype, phi, theta, psi) )

	reduce_EMData_to_root(fftvol, myid, main_node)
	reduce_EMData_to_root(weight, myid, main_node)

	if myid == main_node:   return rec.finish()
	else:                   return model_blank(nx,nx,nx)
'''
	
def combine_ref_ali_to_trans(ref_ali, spider_trans):

	from utilities import combine_params2
	"""
		Serve for ali3d_e
	"""
	res=[]
	new_angle=[ref_ali[2],ref_ali[1],ref_ali[0]]
	res.append(new_angles)
	new_spider_trans=combine_params2(spider_trans,ref_ali[3:6])
	res.append(new_spider_trans)	
	return res
	
'''  Not used anywhere
def Damerau_Levenshtein_Distance(str1, str2):
	"""
	December 11th 2006
	    computes the edit distance
	     which includes deletion, insertion, substitution
		as well as transposition
	"""
	lenStr1 = len(str1);
	lenStr2 = len(str2);
	
	d=[[0 for j in range(lenStr2)] for i in range(lenStr1) ];
	
	for i in range(lenStr1):
		for j in range(lenStr2):
			cost=1;   
			if (str1[i]==str2[j]) :
				cost=0;       # deletion, insertion, substitution
			d[i][j] = min( d[i-1][j]+1 , d[i][j-1]+1 , d[i-1][ j-1] + cost);
	
			if (i>0)&(j>0)&(str1[i]==str2[j-1])& (str1[i-1]==str2[j]):
				d[i][j] = min( d[i][j] , d[i-2][ j-2] + cost); # transposition
	return d[lenStr1-1][lenStr2-1]
'''
'''
def randperm(n):
	"""
	December 11th 2006
	    computes a random permutation of 1 through n
		This was very naive, rewritten by PAP
		WARNING - the list returned in (1 ... n)!
		doesn't seem to be used anywhere, so I commented it out PAP.
	"""
	from random import seed, randint
	seed()
	
	assigned = range(1, n+1) ; # the assignments
	for i in xrange(n):
		k = randint(0, n-1)
		t = assigned[i]
		assigned[i] = assigned[k]
		assigned[k] = t

	return assigned
'''
'''
###########-------------------------------------------------------------------------------------------------###########
#### SPIDER 3D reconstruction related parameters conversion functions   ###
def copy_trans_angles_s(transtr,ptlstr,angstr,defo_fname,stack,Pixel_size,voltage,cs=2,ctf_mul=1,symmetry="c1",deci=1,verify=1,m_radius=0,amp_contrast=.1):
	"""
	1. This command reads output parameters from a SPIDER generated 3D reconstruction, copy all 2D images 
	  (classified in defocus groups) into one hdf stack file, and store the 3D alignment paramters in the headers. The 
	   related headers are all itemized in par_str.   
	2. The program assumes that in the SPIDER 3D reconstruction directory there is a defocus doc file which includes all defoci value of all data
	   sets.
	3. The verification will perform 3D reconstrcutions of a. all particles; b. even number particls; c. odd number particles, and calculate
	   Fourier ring correlation function between volumes created by odd and even number particles.
	4. The user is surpposed to have sufficent experience about performing SPIDER 3D reconstruction and its data structure.	
	"""
	from filt import filt_params,filt_ctf,filt_btwl
	from reconstruction import recons3d_4nn_ctf
	from utilities import read_txt_col,parse_spider_fname,set_arb_params,model_circle,info
	par_str=["psi","theta","phi","sx","sy", "mirror","active","ctf_applied","defocus", "amp_contrast", "voltage", "Cs", "Pixel_size","defgrp"] 
	
	data=read_txt_col(defo_fname,"s")
	img=EMData()
	ptl_defgrp=[]
	snr=1. 
	sign=1.
	n_ptl=-1
	total=0
	angles=[]
	for igrp in xrange(len(data)):
		ftran=parse_spider_fname(transtr,igrp+1)
		fang=parse_spider_fname(angstr,igrp+1)
		ptl_stack=parse_spider_fname(ptlstr,igrp+1)
		trans=read_txt_col(ftran,"s")
		angles=read_txt_col(fang,"s")
		nima = EMUtil.get_image_count(ptl_stack)
		ptl_defgrp.append(nima)
		total+=nima
		for im in xrange(nima):	
			params=[trans[im][2],angles[im][3],angles[im][4],-trans[im][3],-trans[im][4],trans[im][5],1,ctf_mul,data[igrp][2],amp_contrast,voltage,cs,Pixel_size,igrp]		
			n_ptl+=1
			img.read_image(ptl_stack, im)
			if(ctf_mul==1):
				ctf_params = EMAN2Ctf()
				ctf_params.from_dict("defocus": data[igrp][2], "cs": cs, "voltage": voltage, "apix": Pixel_size, "ampcont": amp_contrast)
				img=filt_ctf(img, ctf_params)
			if(deci>1):
				lowf=.5/deci-.05
				highf=.5/deci+.05
				img=filt_btwl(img,lowf,highf)
				img=Util.decimate(img,deci,deci)
			set_arb_params(img,params,par_str)
			ctf_params = EMAN2Ctf()
			ctf_params.from_dict("defocus": data[igrp][2], "cs": cs, "voltage": voltage, "apix": Pixel_size, "ampcont": amp_contrast)
			img.set_attr('ctf', ctf_params)
			img.set_attr('ctf_applied', 1)
			img.write_image(stack,n_ptl)
	m_radius=0
	if(verify==1): get_3Dfsc_stack(stack,m_radius,sign=1,snr=1,symmetry="c1")	

def copy_trans_angles_a(par_list,sp_stack,hdf_stack,new_dim,voltage,Pixel_size,ctf_mul,amp_contrast=.1,cs=2,deci=1,verify=1,symmetry="c1",m_radius=0,fmat="a",skip="C"):
	"""
		1. Read 3D reconstruction parameters from one arbitary format text file
		and store them as headers in one stack hdf file. 
		2. The parameter text file records the parameters with such sequential order:
			psi,theta,phi,sx,sy,defocus1, defocus2 
	"""
	from utilities import read_txt_col,set_arb_params,shift,window2d
	from filt import fshift,filt_ctf
	
	data=read_txt_col(par_list,fmat,skip)
	defocus=[]
	angles=[]
	shift=[]
	par_str=["psi","theta","phi","sx","sy","active","ctf_applied","defocus", "amp_contrast", "voltage", "Cs", "Pixel_size"]
	for i in xrange(len(data)):
		angles.append(data[i][1:4])
		shift.append(data[i][4:6])
		xval=(data[i][8]+data[i][9])/2.
		defocus.append(xval)
	img=EMData()
	for i in xrange(len(data)):
		img.read_image(sp_stack,i)
		if(ctf_mul>0):
			ctf_params = EMAN2Ctf()
			ctf_params.from_dict("defocus": defocus[i], "cs": cs, "voltage": voltage, "apix": Pixel_size, "ampcont": amp_contrast)
			img=filt_ctf(img, ctf_params)
		img=fshift(img,-shift[i][0],-shift[i][1])
		if(deci>1): 
			lowf=.5/deci-.05
			highf=.5/deci+.05
			img=filt_btwl(img,lowf,highf)
			img=Util.decimate(img,deci,deci)
		if(new_dim>0): img=window2d(img,int(new_dim),int(new_dim))
		params=[angles[i][0],angles[i][1],angles[i][2],0,0,1,ctf_mul,defocus[i],amp_contrast,voltage,cs,Pixel_size] # normal order
		set_arb_params(img,params,par_str)
		ctf_params = EMAN2Ctf()
		ctf_params.from_dict("defocus": defocus[i], "cs": cs, "voltage": voltage, "apix": Pixel_size, "ampcont": amp_contrast)
		img.set_attr('ctf', ctf_params)
		img.set_attr('ctf_applied', 1)
		img.write_image(hdf_stack,i)
	sign=1
	snr=1
	if(verify==1): get_3Dfsc_stack(hdf_stack,m_radius,sign,snr,symmetry)

def get_3Dfsc_stack(stack,m_radius=0,sign=1,snr=1,symmetry="c1",dres="fsc001",fval="val.spi",fvol="vol.spi",fvole="vole.spi",fvolo="volo.spi"):
	"""
		Calculate 3D vlomes from a stack file, and calculate the 
		Fourier ring correlation curve
	"""
	from filt import filt_params,filt_ctf,filt_btwl
	from reconstruction import recons3d_4nn_ctf
	from utilities import model_circle,info
	
	total=EMUtil.get_image_count(stack)
	list_p=[]
	for i in xrange(total):						 
		list_p.append(i)
	print "start reconstruction"						 
	val=recons3d_4nn_ctf(stack,list_p, snr, sign, symmetry)
	drop_image(val,"val.spi")
	print "finished"
	plist=[]
	for i in xrange(0,total,2):
		plist.append(i) 
	vole=recons3d_4nn_ctf(stack, plist, snr, sign, symmetry)		 
	drop_image(vole,"vole.spi")
	print "even finished"
	del plist
	plist=[]
	for i in xrange(1,total,2):
		plist.append(i) 
	volo=recons3d_4nn_ctf(stack, plist, snr, sign, symmetry)		 
	drop_image(volo,"volo.spi")
	print "odd finished"
	nx=vole.get_xsize()
	if(m_radius<=1): m_radius=int(nx/2)-10
	mask=model_circle(m_radius, nx, nx, nx)
	s=info(vole,mask)
	vole-=s[0]
	s=info(volo)
	volo-=s[0]
	vole*=mask
	volo*=mask
	res=fsc(volo,vole,1,"fsc001")
	fl, fh = filt_params(res)
	vol=filt_btwl(val,fl,fh)
	drop_image(vol,"vol.spi")
##############

'''


" From statistics "
def memory(since=0.0):
	'''Return memory usage in bytes.
	'''
	return _VmB('VmSize:') - since

def resident(since=0.0):
	'''Return resident memory usage in bytes.
	'''
	return _VmB('VmRSS:') - since

def stacksize(since=0.0):
	'''Return stack size in bytes.
	'''
	return _VmB('VmStk:') - since

def var_bydef(vol_stack, vol_list, info=None):
	"""  var_bydef calculates variances from the definition
	"""
	import sys
	if(info == None): info = sys.stdout
	average = EMData()
	average.read_image( vol_stack, vol_list[0] )
	average.to_zero()
	nimg = 0
	info.write( "Calculating average:" )
	info.flush()
	for ivol in vol_list:
		curt = EMData()
		curt.read_image(vol_stack, ivol )
		Util.add_img(average, curt)
		#average += curt
		if(nimg % 50==0 ) :
			info.write( "\n" )
			info.write( " %4d " % nimg )
			info.flush()
		nimg += 1
		info.write( "." )
		info.flush( )
	info.write( "\n" )
	    
	average /= len(vol_list)

	varis = average.copy()
	varis.to_zero()

	nimg=0
	info.write( "Calculating variance:" )
	info.flush()
	for ivol in vol_list:
		curt = EMData()
		curt.read_image(vol_stack,ivol)
		Util.sub_img(curt, average)
		#curt -= average
		#curt *= curt
		#varis += curt
		Util.add_img2(varis, curt)
		if(nimg % 50==0 ) :
			info.write( "\n" )
			info.write( " %4d " % nimg )
			info.flush()
		nimg=nimg+1
		info.write( "." )
	info.flush( )
	info.write( "\n" )
	
	return varis / (len(vol_list)-1)
	

def _VmB(VmKey):
	'''Private.
	'''
	
	import os
	_proc_status = '/proc/%d/status' % os.getpid()

	_scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
		  'KB': 1024.0, 'MB': 1024.0*1024.0}

	# get pseudo file  /proc/<pid>/status
	#print _scale
	try:
		t = open(_proc_status)
		v = t.read()
		t.close()
	except:
		return 0.0  # non-Linux?
	# get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
	i = v.index(VmKey)
	v = v[i:].split(None, 3)  # whitespace
	if len(v) < 3:
		return 0.0  # invalid format?
	# convert Vm value to bytes
	#print _scale[v[2]],v[2],v[1]
	return float(v[1]) * _scale[v[2]]

'''####################################################################################################
################################  K-MEANS WITH PROBABILITY OF WEIGHTS CALCULATION #####################
####################################################################################################'''
'''
def K_means_prob(images, K, mask=None, init_method="Random"):
	
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

	e.to_zero()
	
	[Cls, p] = init_Kmeans_prob(images,N,K,e,size,mask,init_method,flag)
	
	[Cls, p] = Kmeans_step_prob(images,N,K,p,Cls,size,mask,flag)
	
	for k in xrange(K):
		out2 = open("class"+str(k),"w")
		nk=0
		lst = []
		for n in xrange(N):
			dm = p[n][k]
			fl = True
			for  k1 in xrange(K):
				if(k!=k1):
					if(dm < p[n][k1]):  fl = False
			if(fl):
				nk +=1
				lst.append(n)
		out2.write("\n%s\t%d\t%s\t%d\n" % ("Cluster no:",k,"No of Object = ",nk))
		for kk in xrange(nk):
			out2.write("\t%d" % (lst[kk]))
			if(float(kk+1)/10==int(kk+1)/10):
				out2.write("\n")
		out2.write("\n")
		out2.close()
		del  lst
	

def init_Kmeans_prob(images,N,K,e,size,mask,init_method,flag):
	
	from random import seed,randint
	class Cluster:
		def __init__(self,C):
			self.C   = C
			#self.prevC = prevC
	
	Cls = []	
	e.to_zero()
	
	no = [0]*K
	p = [[0.0]*K for n in xrange(N)]
	
	seed()
	if(init_method == "Random"):

		for k in xrange(K):
			Cls.append(Cluster(None))
			Cls[k].C = e.copy()
			
		for n in xrange(N):
			
			k = randint(0,(K-1))
			
			if(flag):
				e.read_image(images,n)
			else:
				e = images[n].copy()
			
			Cls[k].C += e
			no[k] +=1
	
	for k in xrange(K):
		Cls[k].C /= no[k]
		Cls[k].C.write_image("init_K_means_average.spi",k)

	del Cluster
	
	return [Cls, p]
	
def Kmeans_step_prob(images,N,K,p,Cls,size,mask,flag):
	from math import exp
	alpha = 0.0
	e=EMData()
	e = Cls[0].C.copy()
	e.to_zero()
	c = [e.copy() for k in xrange(K)]
	cnt=0
	d = [0]*K
	for iterat in xrange(100):
		alpha += 0.1
		change =  False
		while (change==False):
			print "Iteration Number :",cnt, alpha
			cnt+=1
			#change = False
			out = open("kmns_"+str(cnt)+"_"+str(alpha)+".txt","w")
			#out1 = open("kmns_nornorm_"+str(alpha)+"_"+str(cnt)+"_"+str(MAX)+".txt","w")
			out.write("%s\t%d\n" % ("ITERATION NUMBER IS : ",cnt))
			x1 = 0
			for n in xrange(N):
				if(flag):
					e.read_image(images,n)
				else:
					e = images[n].copy()

				sumd = 0.0
				for k in xrange(K):
					if(mask==None):
						dist2 = e.cmp("SqEuclidean",Cls[k].C)
					else:
						dist2 = e.cmp("SqEuclidean",Cls[k].C,{"mask":mask})
					d[k] = exp(-alpha*dist2)				
					sumd += d[k]

				for k in xrange(K): p[n][k] = d[k]/sumd
				#out1.write("%f\t" % (p[n][k]))
				#out1.write("\n")

				#for  likmn  in xrange(MAX):
				#	sump = 0.0
				#	for k in xrange(K):
				#		p[n][k] = pow(p[n][k],alpha)
				#		sump += p[n][k]
				#	for k in xrange(K): p[n][k] /= sump

				for k in xrange(K): out.write("%f\t" % (p[n][k]))
				out.write("\n")
			for k in xrange(K): c[k] = Cls[k].C.copy()

			Cls = centroid_prob(images,N,K,Cls,p,mask,flag)
			change = True
			for k in xrange(K):
				dot = Cls[k].C.cmp("dot", c[k], {"normalize":1,"negative":0})
				if(dot < 0.999): change = False

			print "ITERATION " ,alpha,cnt,change

			for k in xrange(K): Cls[k].C.write_image("K_means_average_"+str(cnt)+"_"+str(alpha)+".spi",k)
			#out1.close()
			out.close()
			#alpha += 0.05
			
	return [Cls,p]


def centroid_prob(images,N,K,Cls,p,mask,flag):
	
	from utilities import ttime
	print ttime()
	e=EMData()
	for k in xrange(K): Cls[k].C.to_zero()

	sumw = [0.0]*K
	for n in xrange(N):
		if(flag):
			e.read_image(images,n)
		else:
			e = images[n].copy()
		
		for k in xrange(K):
			Cls[k].C = Util.madn_scalar(e, Cls[k].C, p[n][k])
			sumw[k] += p[n][k]
	for k in xrange(K): Cls[k].C /= sumw[k]
	print ttime()
		
	return Cls
#######################################################################################################################
'''


def create_write_gridprojections(volume, filepattern, anglelist,
				     kb_K=6, kb_alpha=1.75,
				     start=1, npad=2, dofakealign=False,
				     noise_sigma=0.0):
	"""Given a volume, a set of projection angles, and Kaiser-Bessel
	   window parameters, use gridding to generate projections along
	   the various angles and write them out to disk using the passed-in
	   filepattern.

	   Usage:
	     anglist = voea(60.)
	     vol = get_image("model.tcp")
	     kb_K = 6
	     kb_alpha = 1.75
	     create_write_gridprojections(vol,"prj{***}.tcp",anglist,kb_K,kb_alpha)
	"""
	from utilities import parse_spider_fname, model_gauss_noise
	from fundmentals import grotshift2d
	myparams = {"angletype":"SPIDER",
		    "anglelist":anglelist,
		    "kb_alpha":kb_alpha,
		    "kb_K":kb_K,
		    "npad":npad}
	projvol = volume.project("fourier_gridding",myparams)
	nangles = len(anglelist) / 3
	nx = volume.get_xsize()
	ny = volume.get_ysize()
	m = volume.get_xsize()
	n = npad*m
	kb = Util.KaiserBessel(kb_alpha, kb_K, m/2, kb_K/(2.*n),n)
	for i in xrange(nangles):
		if projvol.get_zsize() > 1:
			# extract the i'th z-plane
			proj = projvol.get_clip(Region(0,0,i,nx,ny,1))
		else:
			proj = projvol
		if noise_sigma > 0.0:
			proj += model_gauss_noise(noise_sigma, m, m)
		if dofakealign:
			# Randomly rotate and translate forward and back
			# to simulate an alignment
			from random import uniform
			from math import pi
			angle = uniform(0,2*pi)
			dx = uniform(0,10.)
			dy = uniform(0,10.)
			p2 = grotshift2d(proj, angle, kb, npad, dx, dy)
			p2parback = grotshift2d(p2, 0., kb, npad, -dx, -dy)
			p2back = grotshift2d(p2parback, -angle, kb, npad, 0., 0.)
			proj = p2back
		projname = parse_spider_fname(filepattern, i+start)
		proj.write_image(projname, 0, EMUtil.ImageType.IMAGE_SINGLE_SPIDER)

'''  From project '''
"""
def cre(volume, filepattern, anglelist,
				 kb_K=6, kb_alpha=1.75,
				 start=1, npad=2, dofakealign=False,
				 noise_sigma=0.0):
#    Given a volume, a set of projection angles, and Kaiser-Bessel
#       window parameters, use gridding to generate projections along
#       the various angles and write them out to disk using the passed-in
#       filepattern.
#
#       Usage:
#         anglist = voea(60.)
#         vol = get_image("model.tcp")
#         kb_K = 6
#         kb_alpha = 1.75
#         create_write_gridprojections(vol,"prj{***}.tcp",anglist,kb_K,kb_alpha)
#
    myparams = {"angletype":"SPIDER",
		"anglelist":anglelist,
		"kb_alpha":kb_alpha,
		"kb_K":kb_K,
		"npad":npad}
    #proj = volume.project("fourier_gridding",myparams)
    #projname = parse_spider_fname(filepattern, start)
    #proj.write_image(projname, 0, EMUtil.ImageType.IMAGE_SINGLE_SPIDER)
    " This does not work in C code, DO NOT USE!"
    return volume.project("fourier_gridding",myparams)
"""

'''

def pw2_adjustment(vol_in, vol_model, fscdoc, q, shell_of_noise = 20, cut_off_fsc = .2, save_curves = False):
	"""
		Enforce power spectrum of a 3D volume to that of corresponding 
		X-ray structure. The adjustment is equivalent to applying a 
		band_pass filter to the volume.
	"""
	import os
	import types
	from utilities import get_image, read_text_file
	from fundamentals import rops_table
	from filter import filt_table
	if type(vol_in) is types.StringType : vol_to_be_corrected = get_im( vol_in )
	else :                                vol_to_be_corrected = vol_in
	nx    = vol_to_be_corrected.get_xsize()
	r1    = nx//2 - 1 - shell_of_noise
	r2    = nx//2 - 1
	m1    = model_circle(r1, nx, nx, nx)
	m1   -= 1.
	m1   *=-1.
	m2    = model_circle(r2, nx, nx, nx)
	m_t   = m1*m2
	tmp_l = info(m_t)
	from utilities import gauss_edge
	mask  = gauss_edge(m_t)
	vol_n = vol_to_be_corrected*mask
	t = info(vol_n)
	vol_n -= t[0]
	t = info(vol_to_be_corrected)
	vol_to_be_corrected -= t[0]
	p3o   = rops_table(vol_to_be_corrected)
	p3b   = rops_table(vol_n)
	for i in xrange(len(p3b)): p3b[i] /= tmp_l[0]
	if type(vol_model) is types.StringType : vol_model = get_im( vol_model )
	pu    = rops_table(vol_model)
	fsc   = read_text_file(fscdoc, 1)
	p_n   = []
	out = Util.pw_extract(p3b[1:], 3, 3)
	p_n.append(0.)
	for i in xrange(1,len(p3b)):
		j = i*2
		p_n.append(out[j])
	tmp_snr = []
	eps   = 0.882
	a     = 10.624
	for j in xrange(len(fsc)):
		if fsc[j] < cut_off_fsc:
			lowf  = float(j)/len(fsc)/2.
			break
	highf = lowf+.05
	order=2.*log(eps/sqrt(a**2-1))/log(lowf/highf)
	rad   = lowf/eps**(2./order)
	H     = []	       
	for j in xrange(1,len(fsc)):
		if fsc[j] < cut_off_fsc  and fsc[j] > 0. :
			freq = float(j)/len(fsc)/2.
			snr = 1./sqrt(1.+(freq/rad)**order)
		elif fsc[j] <= 0.:               snr = 0.0 
		else:                            snr = fsc[j]/(fsc[j]+1.)
		tmp_snr.append(snr)
		tmpx = (p3o[j]-(q**2*p3b[j]+p_n[j]*(1.-q**2)))
		if    tmpx < 0. and float(j)/len(fsc) < highf : ERROR("q is too large", "pw2_adjustment")
		elif  tmpx < 0. : H.append(0.0)
		else:
			x0   = pu[j]**.5*snr
			x1   = tmpx**.5
			H.append(x0/x1)
	out_vol = filt_table(vol_to_be_corrected, H)
	if save_curves :
		from fundamentals import rops_textfile
		rops_textfile(out_vol,  "roo_adjusted.txt")
		rops_textfile(vol_model,"roo_model.txt")
		rops_textfile(vol_to_be_corrected,"roo_original.txt")
		drop_spider_doc("snr.txt", tmp_snr)
		drop_spider_doc("H.txt", H)
		drop_spider_doc("bg.txt",p_n )
		drop_spider_doc("roo_noise.txt",p3b)
		drop_spider_doc("roo_original.txt",p3o)
		drop_image(out_vol,"vol_adjusted.spi", "s")
		drop_image(mask,"vol_mask.spi", "s")
	return  out_vol
'''


'''  From fundamentals  '''
	
def fpol_stack(stack, stack_copy, nnx, nny=0, nnz=1, RetReal = True):
	nima = EMUtil.get_image_count(stack)
	ima = EMData()
	imao = EMData()
	for im in xrange(nima):
		ima.read_image(stack,im)
		imao=ima.FourInterpol(nnx, nny, nnz, RetReal)
		imao.write_image(stack_copy, im)

def select_stack(stack, select_stack):
	nima = EMUtil.get_image_count(stack)
	ima = EMData()
	ima.read_image(stack,0)
	nclass = ima.get_attr('nclass')
	class2 = []
	for im in xrange(nima):
		ima.read_image(stack,im)
		n = ima.get_attr('assign')
		if n not in class2:
			class2.append(n)		
			ima.write_image(select_stack, n)
				
def rot_shift(image, RA, method='quadrilinear'):
	return image.rot_scale_trans(RA)

'''
def rtsh(image, angle, sx=0., sy=0.):
	m=image.get_xsize()
	# padd two times
	npad=2
	n=m*npad
	# support of the window
	k=6
	alpha=1.75
	kb = Util.KaiserBessel(alpha, k, m/2, k/(2.*n), n)
	return grotshift2d(image, angle, kb, npad, sx, sy)
'''

###############################################################################################
# draft ali2d_rac_MPI
###############################################################################################

def ali2d_rac_MPI(stack, maskfile = None, kmeans = 'None', ir = 1, ou = -1, rs = 1, nclass = 2, maxit = 10, maxin = 10, check_mirror = False, rand_seed = 10):
	from global_def import MPI
	from utilities  import bcast_EMData_to_all, reduce_EMData_to_root, bcast_number_to_all
	from utilities  import send_EMData, recv_EMData
	from utilities  import model_circle, combine_params2
	from utilities  import get_arb_params, set_arb_params
	from statistics import MPIlogfile_init, MPIlogfile_print, MPIlogfile_end
	from statistics import kmnr, kmn
	from random     import randint, seed, shuffle
	from alignment  import Numrinit, ringwe, ang_n
	from copy       import deepcopy
	import time
	import sys
	from mpi 	  import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_recv, mpi_send
	from mpi 	  import MPI_SUM, MPI_FLOAT, MPI_INT, MPI_LOR

	# [id]   part of code different for each node
	# [sync] synchronise each node
	# [main] part of code just for the main node
	# [all]  code write for all node

	# debug mode
	DEBUG = False
	K_th  = 1

	# To work
	if kmeans == 'None':
		start_k = False
	else:
		start_k = True
	
	# ------------------- MPI init ----------------------------- #
	# init
	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)

	# chose a random node as a main one
	main_node = 0
	if myid  == 0:	main_node = randint(0,number_of_proc-1)
	main_node = bcast_number_to_all(main_node,0)
	mpi_barrier(MPI_COMM_WORLD)

	if number_of_proc > nclass:
		if myid == main_node:
			print 'need number of cpus > K'
		return	
		
	t_start = time.time()
	
	# ---------------------------------------------------------- #
	
	seed(rand_seed)
	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); K=int(nclass); max_iter=int(maxit); max_internal=int(maxin);
	nima = EMUtil.get_image_count(stack)
	temp = EMData()
	temp.read_image(stack, 0)
	nx = temp.get_xsize()
	ny = nx

	# default value for the last ring
	if(last_ring==-1): last_ring=nx//2-2

	# somebody dedicated could write a version with option "H" for half rings that would work for ACF functions.
	mode = "F"

	# [all] data structure
	grp_asg    = [] # assignment of groups to node
	GRP        = [] # list of local groups, each group contain index images
	IM         = [] # list of local groups, each group contain images
	AVE        = [] # list of averages (for all groups)

	# [main] 
	if myid == main_node:
		# define assignment of groups to node
		grp_asg = []
		for k in xrange(K): grp_asg.append(k % 	number_of_proc)
		grp_asg.sort()

		if start_k == False:
			# randomly list of index
			tmp = range(0, nima)
			shuffle(tmp)

			list_i = []
			for n in xrange(nima): list_i.append(n % K)
			list_i.sort()

			list_index = [[] for n in xrange(K)]
			n = 0
			for i in list_i:
			    list_index[i].append(tmp[n])
			    n += 1
		else:
			# pick-up index from ave of k-means
			list_index = [[] for n in xrange(K)]
			im = EMData()
			for k in xrange(K):
				im.read_image(kmeans, k, True)
				listim = im.get_attr('members')
				for index in listim: list_index[k].append(int(index))
				if k == 0: print list_index

	## to test
	#return

	# [main] broadcast list of assignment of groups to node
	grp_asg = mpi_bcast(grp_asg, K, MPI_INT, main_node, MPI_COMM_WORLD)
	grp_asg = grp_asg.tolist()

	# [all] display infos
	MPIlogfile_init(myid)
	MPIlogfile_print(myid, '************* ali2d_rac MPI *************\n')
	MPIlogfile_print(myid, 'Input stack                          : %s\n' % stack)
	MPIlogfile_print(myid, 'Mask file                            : %s\n' % maskfile)
	MPIlogfile_print(myid, 'Inner radius                         : %i\n' % first_ring)
	MPIlogfile_print(myid, 'Outer radius                         : %i\n' % last_ring)
	MPIlogfile_print(myid, 'Ring step                            : %i\n' % rstep)
	MPIlogfile_print(myid, 'Maximum iteration                    : %i\n' % max_iter)
	MPIlogfile_print(myid, 'Maximum internal iteration           : %i\n' % max_internal)
	MPIlogfile_print(myid, 'Consider mirror                      : %s\n' % check_mirror)
	MPIlogfile_print(myid, 'Random seed                          : %i\n' % rand_seed)
	MPIlogfile_print(myid, 'Number of classes                    : %i\n' % K)
	MPIlogfile_print(myid, 'Number of images                     : %i\n' % nima)
	MPIlogfile_print(myid, 'Number of cpus                       : %i\n' % number_of_proc)
	MPIlogfile_print(myid, 'ID cpu                               : %i\n' % myid)

	nb_grp = 0
	for k in xrange(K):
		if grp_asg[k] == myid: nb_grp += 1

	MPIlogfile_print(myid, 'Number of classes in this cpu        : %i\n' % nb_grp)
	MPIlogfile_print(myid, 'Output file                          : %s\n\n' % stack)

	# [all] generate LUT group, ex global indice grp 4 5 6 -> local indice in the node 0 1 2
	# [4, 5, 6] -> [-1, -1, -1, -1, 0, 1, 2]
	lut_grp = [-1] * K
	i       = 0
	for k in xrange(K):
		if myid == grp_asg[k]:
			lut_grp[k] = i
			i += 1
	
	# [main] send list_index to each individual group
	for k in xrange(K):
		size = []
		dst  = -1
		tmp  = []

		if myid == main_node:
			if grp_asg[k] != main_node:
				size = len(list_index[k])
				mpi_send(size, 1, MPI_INT, grp_asg[k], 0, MPI_COMM_WORLD)
				mpi_send(list_index[k], size, MPI_INT, grp_asg[k], 0, MPI_COMM_WORLD)
			else:
				GRP.append(list_index[k])
		else:
			if myid == grp_asg[k]:
				size = mpi_recv(1, MPI_INT, main_node, 0, MPI_COMM_WORLD)
				tmp  = mpi_recv(size, MPI_INT, main_node, 0, MPI_COMM_WORLD)
				tmp  = tmp.tolist()
				GRP.append(tmp)

		mpi_barrier(MPI_COMM_WORLD)

	# [all] precalculate rings
	numr     = Numrinit(first_ring, last_ring, rstep, mode)
	wr       = ringwe(numr ,mode)
	lnumr    = numr[len(numr)-1]
	norm_rsd = 0
	for n in xrange(1, len(numr), 3): norm_rsd += numr[n]
	
	# [all] prepare 2-D mask for normalization
	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):
			mask2D = get_image(maskfile)
		else:
			mask2D = maskfile
	else:
		mask2D = model_circle(last_ring,nx,nx)
	if(first_ring > 0):
		tave = model_circle(first_ring-1, nx, ny)
		mask2D -= tave

	#  center is in SPIDER convention
	cnx = int(nx/2) + 1
	cny = int(ny/2) + 1
	
	# [id] read images in each local groups and resample them into polar coordinates
	if DEBUG: MPIlogfile_print(myid, 'prepare images\n')
	for grp in xrange(nb_grp):
		data = []
		im   = 0
		for index in GRP[grp]:
			temp = EMData()
			temp.read_image(stack, index)
			sx = temp.get_attr('sx')
			sy = temp.get_attr('sy')
			alpha_original = temp.get_attr('alpha')
			miri = temp.get_attr('mirror')
			[mean, sigma, qn, qm] = Util.infomask(temp, mask2D, True)
			temp = (temp - mean)/sigma
			alpha_original_n,sxn,syn,mir = combine_params2(0, -sx, -sy, 0, -alpha_original,0,0,0)
			cimage = Util.Polar2Dm(temp, cnx+sxn, cny+syn, numr, mode)
			Util.Frngs(cimage, numr)
			data.append(cimage)
			data[im].set_attr_dict({'alpha':1.0, 'alpha_original':alpha_original, 'sx':sx, 'sy':sy, 'mirror': 0})

			im += 1

		IM.append(deepcopy(data))

	del temp
	del data
	del mask2D

	# use image to 0 for define a blank image (use to reset the average)
	blank = IM[0][0]
	blank.to_zero()

	# prepare the syncronize sequence (for the com between the nodes)
	sync_org = []
	sync_dst = []
	for i in xrange(number_of_proc):
		tmp = range(number_of_proc)
		sync_dst.extend(tmp)
		for j in xrange(number_of_proc):
			sync_org.append(i)

	# debug
	if DEBUG:
		flow_ctrl = [0] * nima
		flow_cnt  = 0
	
	# [id] classification
	again =  1
	it    = -1
	while again and it < (max_iter - 1):
		it += 1
	
		mpi_barrier(MPI_COMM_WORLD)

		if DEBUG:
			flow_ite = 0

			if myid == main_node:
				tg = time.time()
				print '== ITE %d ==' % it

		# [id] compute average send to main_node
		AVE = []
		for k in xrange(K): AVE.append(blank)
		
		if DEBUG:
			MPIlogfile_print(myid, 'after init AVE\n')
			if myid == main_node: t1 = time.time()
		
		for k in xrange(K):
			if myid == grp_asg[k]:
				temp    = kmnr(IM[lut_grp[k]], -1, len(GRP[lut_grp[k]]), -1, numr, wr, check_mirror, max_internal, rand_seed, myid)
				temp   /= len(GRP[lut_grp[k]])
				AVE[k]  = temp.copy()

		if DEBUG and myid == main_node: print 'time average: %d s' % (time.time() - t1)

		if DEBUG: MPIlogfile_print(myid, 'compute ave\n')

		mpi_barrier(MPI_COMM_WORLD)

		for k in xrange(K):
			if myid != main_node:
				if myid == grp_asg[k]:
					send_EMData(AVE[k], main_node, 0)
			else:
				if grp_asg[k] != main_node and grp_asg[k] != -1:
					AVE[k] = recv_EMData(grp_asg[k], 0)

			mpi_barrier(MPI_COMM_WORLD)

		if DEBUG: MPIlogfile_print(myid, 'recv average\n')

		# [main] broadcast each average to other node
		for k in xrange(K):
			tmp = AVE[k].copy()
			bcast_EMData_to_all(tmp, myid, main_node)
			AVE[k] = tmp.copy() # need to use .copy() to store the new im
			mpi_barrier(MPI_COMM_WORLD)

		if DEBUG: MPIlogfile_print(myid, 'after send AVE\n')

		# [all] compute norm for each average
		norm = []
		for k in xrange(K):
			q   = Util.ener(AVE[k], numr)
			res = Util.Crosrng_ew(AVE[k], AVE[k], numr, wr, 0)
			norm.append((2 * q) / float(res['qn']))

		# [id] info display
		MPIlogfile_print(myid, '\n___________ Iteration %d _____________%s\n' % (it + 1, time.strftime('%a_%d_%b_%Y_%H_%M_%S', time.localtime())))
		for k in xrange(K):
			if lut_grp[k] != -1:
				MPIlogfile_print(myid, 'objects in class %3d : %d \n' % (k, len(GRP[lut_grp[k]])))

		MPIlogfile_print(myid, '\n')

		if DEBUG and myid == main_node: t1 = time.time()
		
		# [id] classification
		list_exg = []
		list_dst = []
		again    = 0
		Je_rsd   = 0
		for ngrp in xrange(K):
			if myid == grp_asg[ngrp]:
				exg = []
				dst = []
				for index in xrange(len(GRP[lut_grp[ngrp]])):
					dmin    = 1.0e20
					old_grp = ngrp
					new_grp = -1
					for k in xrange(K):
						if grp_asg[k] != -1: # if the group k is not empty
							if (check_mirror):
								retvals = Util.Crosrng_ew(AVE[k], IM[lut_grp[ngrp]][index], numr, wr, 0)
								qn      = retvals["qn"]
								retvals = Util.Crosrng_ew(AVE[k], IM[lut_grp[ngrp]][index], numr, wr, 1)
								qm      = retvals["qn"]
								q1      = Util.ener(AVE[k], numr)
								q2      = Util.ener(IM[lut_grp[ngrp]][index], numr)
								qn      = max(qn,qm)
								qn      = q1 + q2 - (qn * norm[k])
							else:
								retvals = Util.Crosrng_ew(AVE[k], IM[lut_grp[ngrp]][index], numr, wr, 0)
								qn      = retvals["qn"]
								q1      = Util.ener(AVE[k], numr)
								q2      = Util.ener(IM[lut_grp[ngrp]][index], numr)
								qn      = q1 + q2 - (qn * norm[k])
							if(qn < dmin):
								dmin    = qn
								new_grp = k

					Je_rsd += (dmin / float(norm_rsd))
					
					if new_grp < 0:
						print  'Error in assignment of objects to averages'
						break
					
					## TO TEST
					#if ngrp    == 0: new_grp = 1
					#if new_grp == 0: new_grp = 1
				
					if old_grp != new_grp:
						again = 1
						exg.append(GRP[lut_grp[ngrp]][index])
						dst.append(new_grp)
					
				# store list exg and dst for each loc groups
				list_exg.append(exg)
				list_dst.append(dst)

		if DEBUG and myid == main_node: print 'time classification: %d s' % (time.time() - t1)

		mpi_barrier(MPI_COMM_WORLD)

		if DEBUG: MPIlogfile_print(myid, 'after classification\n')

		# [sync] gather value of criterion Je_rsd
		Je_rsd = mpi_reduce(Je_rsd, 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
		Je_rsd = mpi_bcast(Je_rsd, 1, MPI_FLOAT, main_node, MPI_COMM_WORLD)
		Je_rsd = Je_rsd.tolist()
		Je_rsd = Je_rsd[0]

		MPIlogfile_print(myid, 'Criterion %11.4e\n' % Je_rsd)

		# [sync] with the other node
		again = mpi_reduce(again, 1, MPI_INT, MPI_LOR, main_node, MPI_COMM_WORLD)
		again = mpi_bcast(again, 1, MPI_INT, main_node, MPI_COMM_WORLD)
		again = again.tolist()
		again = again[0]
		
		if DEBUG:
			for grp in xrange(len(list_exg)):
				for index in list_exg[grp]:
					flow_ctrl[index] += 1
			if myid == main_node: t1 = time.time()

		# COMMUNICATION
		mpi_barrier(MPI_COMM_WORLD)
		for node in xrange(number_of_proc * number_of_proc):
			# define the origin and destination node
			org_node = sync_org[node]
			dst_node = sync_dst[node]

			# define list of group inside org node and dst node
			org_grp = []
			dst_grp = []
			for k in xrange(K):
				if grp_asg[k] == org_node:
					org_grp.append(k)
				if grp_asg[k] == dst_node:
					dst_grp.append(k)

			bag           = []
			pack          = []
			pack_size     = []
			bag_im        = []
			pack_im       = []
			lst_remove    = [[] for n in xrange(len(org_grp))]
			if myid == org_node:
				# prepare data and pack
				for outgrp in dst_grp:
					bag        = []
					bag_im     = []
					for ingrp in org_grp:
						remove   = []
						size_exg = len(list_exg[lut_grp[ingrp]])
						for n in xrange(size_exg):
							if list_dst[lut_grp[ingrp]][n] == outgrp:
								index = list_exg[lut_grp[ingrp]][n]
								bag.append(index)
								pos = GRP[lut_grp[ingrp]].index(index)
								bag_im.append(IM[lut_grp[ingrp]][pos])
								remove.append(index)

						lst_remove[lut_grp[ingrp]].extend(remove)
			
					pack.extend(bag)
					pack_im.extend(bag_im)
					pack_size.append(len(bag))
						
			# first send the pack_size to all
			pack_size = mpi_bcast(pack_size, len(dst_grp), MPI_INT, org_node, MPI_COMM_WORLD)
			pack_size = pack_size.tolist()

			mpi_barrier(MPI_COMM_WORLD)

			if DEBUG:
				data_send = 0
				data_recv = 0
				data_in   = 0
			
			# send-recv
			if myid != org_node:
				pack    = []
				pack_im = []
			for n in xrange(sum(pack_size)):
				if org_node != dst_node:
					if myid == org_node:
						# send index
						mpi_send(pack[n], 1, MPI_INT, dst_node, 0, MPI_COMM_WORLD)
						# send header [alpha, sx, sy, mirror]
						head = get_arb_params(pack_im[n], ['alpha', 'alpha_original', 'sx', 'sy', 'mirror'])
						mpi_send(head, 5, MPI_FLOAT, dst_node, 0, MPI_COMM_WORLD)
						# send image
						send_EMData(pack_im[n], dst_node, 1)
						
						if DEBUG: data_send += 1
					
					if myid == dst_node:
						# recv index
						index = mpi_recv(1, MPI_INT, org_node, 0, MPI_COMM_WORLD)
						index = index.tolist()
						pack.append(index[0])
						# recv header
						head    = mpi_recv(5, MPI_FLOAT, org_node, 0, MPI_COMM_WORLD)
						head    = head.tolist()
						head[4] = int(head[4]) # mirror is integer
						# recv image
						img   = recv_EMData(org_node, 1)
						set_arb_params(img, head, ['alpha', 'alpha_original', 'sx', 'sy', 'mirror'])
						pack_im.append(img)
						
						if DEBUG: data_recv += 1
				else:
					if DEBUG: data_in += 1

				mpi_barrier(MPI_COMM_WORLD)
	
			if DEBUG:
				flow_cnt += data_send
				flow_ite += data_send
				if dst_node != org_node:
					if myid == org_node:
						MPIlogfile_print(myid, 'send data %s %d ---> %d %2s: exp [%d] int [%d]\n' % (org_grp, org_node, dst_node, str(dst_grp).ljust(50, ' '), data_send, data_in))
					if myid == dst_node:
						MPIlogfile_print(myid, 'recv data %s %d <--- %d %2s: exp [%d] int [%d]\n' % (dst_grp, dst_node, org_node, str(org_grp).ljust(50, ' '), data_recv, data_in))
				else:
					if myid == org_node:
						MPIlogfile_print(myid, 'loca data %s %d <--> %d %2s: exp [%d] int [%d]\n' % (dst_grp, dst_node, org_node, str(org_grp).ljust(50, ' '), data_recv, data_in))

			if myid == org_node:
				# remove index and images in each group
				for ingrp in org_grp:
					for index in lst_remove[lut_grp[ingrp]]:
						pos = GRP[lut_grp[ingrp]].index(index)
						trash = GRP[lut_grp[ingrp]].pop(pos)
						trash = IM[lut_grp[ingrp]].pop(pos)
						del trash
				
			if myid == dst_node:
				i = 0
				j = 0
				for ingrp in dst_grp:
					bag    = []
					bag_im = []
					for n in xrange(j, j + pack_size[i]):
						bag.append(pack[n])
						bag_im.append(pack_im[n])
						j += 1
						
					i   += 1

					# add index and images in each group
					GRP[lut_grp[ingrp]].extend(bag)
					IM[lut_grp[ingrp]].extend(bag_im)

			# [sync]
		        mpi_barrier(MPI_COMM_WORLD)

		if DEBUG and myid == main_node:
			print 'time comm: %d s' % (time.time() - t1)
			print '\n'
			print 'time iteration: %d s\n\n' % (time.time() - tg)
		
		# [sync]
		mpi_barrier(MPI_COMM_WORLD)

		if DEBUG:
			flow_ite = mpi_reduce(flow_ite, 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
			flow_ite = flow_ite.tolist()

		flag = True
		if flag:
			# search the groups to remove
			grp_rmv = [0] * K
			nbr_obj = [0] * K
			for k in xrange(K):
				if myid == grp_asg[k]:
					if len(GRP[lut_grp[k]]) < K_th:
						grp_rmv[k] = 1
						nbr_obj[k] = len(GRP[lut_grp[k]])

			for k in xrange(K):
				# if the group k need to removed
				if grp_rmv[k] > 0:
					if nbr_obj[k] != 0:
						MPIlogfile_print(myid, 'Warning: class %d not enough images, remove this group\n' % k)
						org_mv = k
						dst_mv = -1

						## TO TEST
						#MPIlogfile_print(myid, 'list lut_grp: %s\n' % lut_grp)
						#MPIlogfile_print(myid, 'list grp rmv: %s\n' % grp_rmv)
						
						# search if in the node there is another group to move objects
						for q in xrange(K):
							if lut_grp[q] != -1 and grp_rmv[q] == 0:
								dst_mv = q

						## TO TEST
						#MPIlogfile_print(myid, 'grp destination %d\n' % dst_mv)
								
						# move object if there is a group destination
						if dst_mv != -1:
							#MPIlogfile_print(myid, 'from %d -> %d\n' % (org_mv, dst_mv))
							GRP[lut_grp[dst_mv]].extend(GRP[lut_grp[org_mv]])
							IM[lut_grp[dst_mv]].extend(IM[lut_grp[dst_mv]])
							MPIlogfile_print(myid, 'move objects grp %d -> grp %d\n' % (org_mv, dst_mv))
					else:
						MPIlogfile_print(myid, 'Warning: class %d is empty, remove this group.\n' % k)
	
			# clean groups
			for k in xrange(K - 1, -1, -1):
				if grp_rmv[k] > 0:
					#MPIlogfile_print(myid, 'lut_grp[k] and k: %d %d\n' % (lut_grp[k], k))
					del GRP[lut_grp[k]]
					del IM[lut_grp[k]]
							
			# broadcast the grp remove
			grp_rmv = mpi_reduce(grp_rmv, K, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
			grp_rmv = mpi_bcast(grp_rmv, K, MPI_INT, main_node, MPI_COMM_WORLD)
			grp_rmv = grp_rmv.tolist()

			# new group assign node
			for k in xrange(K):
				if grp_rmv[k] > 0: grp_asg[k] = -1			

			# new list loc group assign
			lut_grp = [-1] * K
			i       = 0
			for k in xrange(K):
				if myid == grp_asg[k]:
					lut_grp[k] = i
					i += 1

			# check empty node
			kill = 0
			if len(GRP) == 0:
				print 'Error: cpu %d is empty (no group), kill all process.\n' % myid
				MPIlogfile_print(myid, 'Error: cpu %d is empty (no group), kill all process.\n' % myid)
				kill = 1
			
			mpi_barrier(MPI_COMM_WORLD)
			kill = mpi_reduce(kill, 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
			kill = mpi_bcast(kill, 1, MPI_INT, main_node, MPI_COMM_WORLD)
			kill = kill.tolist()
			kill = kill[0]

			if kill > 0: return

		else:
			# check empty class
			kill = 0
			for k in xrange(K):
				if myid == grp_asg[k]:
					if len(GRP[lut_grp[k]]) < 1:
						print 'Error: class %d is empty' % k
						MPIlogfile_print(myid, 'Error: class %d is empty, kill all process.\n' % k)
						kill = 1

				mpi_barrier(MPI_COMM_WORLD)
			kill = mpi_reduce(kill, 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
			kill = mpi_bcast(kill, 1, MPI_INT, main_node, MPI_COMM_WORLD)
			kill = kill.tolist()
			kill = kill[0]

			if kill > 0: return

	if DEBUG:
		flow_ctrl = mpi_reduce(flow_ctrl, nima, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
		flow_ctrl = flow_ctrl.tolist()

		print myid, ': flow count [%d]' % flow_cnt

		if myid == main_node:
			print '\n ############ flow control ####################\n'
			#print TXT_flow
			#print flow_ctrl
			print '\n'
			print 'sum: %d' % sum(flow_ctrl)
	
	mpi_barrier(MPI_COMM_WORLD)

	# display the number of objects in each group
	ngrp = [0] * K
	for k in xrange(K):
		if myid == grp_asg[k]: ngrp[k] = len(GRP[lut_grp[k]])
	ngrp = mpi_reduce(ngrp, K, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
	ngrp = mpi_bcast(ngrp, K, MPI_INT, main_node, MPI_COMM_WORLD)
	ngrp = ngrp.tolist()
	MPIlogfile_print(myid, '\n== Object in each Class ==\n')
	count = 0
	for k in xrange(K):
		if ngrp[k] == 0:
			MPIlogfile_print(myid, 'class %3d : %d [removed]\n' % (k, ngrp[k]))
		else:
			MPIlogfile_print(myid, 'class %3d : %d\n' % (k, ngrp[k]))
			count += 1

	MPIlogfile_print(myid, 'number of classes: %d\n' % count)
	
	MPIlogfile_print(myid, '\nPREPARE RESULTS\n')

	'''
	# update ave according the empty group
	newAVE = []
	L      = 0
	for k in xrange(K):
		if ngrp[k] != 0:
			tmp = AVE[k].copy()
			newAVE.append(tmp)
			L += 1
	K = L
	'''

	# [main] align class averages and transfer parameters to individual images
	talpha = [0] * K
	tmir   = [0] * K
	if myid == main_node:
		for k in xrange(K): AVE[k].set_attr_dict({'alpha':1.0, 'mirror':0})
		kmn(AVE, numr, wr, check_mirror, max_iter, rand_seed)
		for k in xrange(K):
			talpha[k] = ang_n(AVE[k].get_attr('alpha'), mode, lnumr)
			tmir[k]   = AVE[k].get_attr('mirror')

		del AVE
	
	# [all] broadcast talpha and tmir
	talpha = mpi_bcast(talpha, K, MPI_FLOAT, main_node, MPI_COMM_WORLD)
	talpha = talpha.tolist()
	tmir   = mpi_bcast(tmir, K, MPI_INT, main_node, MPI_COMM_WORLD)
	tmir   = tmir.tolist()
	
	mpi_barrier(MPI_COMM_WORLD)

	#  write out the alignment parameters to headers
	del temp
	temp = EMData()

	if DEBUG and myid == main_node: t1 = time.time()

	for k in xrange(K):
		if myid == grp_asg[k]:
			for n in xrange(len(GRP[lut_grp[k]])):
				#  First combine with angle of the average
				alpha  = ang_n(IM[lut_grp[k]][n].get_attr('alpha'), mode, lnumr)
				mirror =  IM[lut_grp[k]][n].get_attr('mirror')
				alpha, tmp, it, mirror = combine_params2(alpha, 0, 0, mirror, talpha[k], 0, 0, tmir[k])
			
				#  Second combine with given alignment
				alpha_original = IM[lut_grp[k]][n].get_attr('alpha_original')
				sx    =  IM[lut_grp[k]][n].get_attr('sx')
				sy    =  IM[lut_grp[k]][n].get_attr('sy')
				alpha_original_n, sxn, syn, mir = combine_params2(0, -sx, -sy, 0, -alpha_original, 0, 0, 0)
				alphan, sxn, syn, mir           = combine_params2(0, -sxn, -syn, 0, alpha, 0, 0, mirror)
				temp.read_image(stack, GRP[lut_grp[k]][n], True)
				temp.set_attr_dict({'alpha':alphan, 'sx':sxn, 'sy':syn, 'mirror': mir, 'nclass':K, 'ref_num':k})
				temp.write_image(stack, GRP[lut_grp[k]][n], EMUtil.ImageType.IMAGE_HDF, True)

		mpi_barrier(MPI_COMM_WORLD)

	if DEBUG and myid == main_node: print 'time write header: %d s' % (time.time() - t1)

	del temp
	del IM

	MPIlogfile_print(myid, '\n Time: %d s\n' % (time.time() - t_start))
	MPIlogfile_end(myid)



'''

def ali2d_cross_res(stack, outdir, maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", center=1, maxit=0, CTF=False, snr=1.0, user_func_name="ref_ali2d"):
 	"""
		Split data into odd and even data sets and align them seperately
		Cross resolution alignment
	"""
	import os
	from utilities 		import model_circle, model_blank, combine_params2, drop_image
	from utilities      import get_input_from_string, get_image, get_arb_params, set_arb_params
	from fundamentals 	import rot_shift2D
	from statistics 	import add_oe_series, ave_series_ctf, ave_series, fsc_mask
	from alignment 		import Numrinit, ringwe, ali2d_single_iter, align2d
	from filter 		import filt_table, filt_ctf
	from morphology     import ctf_2

	from utilities import print_begin_msg, print_end_msg, print_msg
	import	types
		
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
	print_msg("Data with CTF               : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	if auto_stop:  print_msg("Stop iteration with         : criterion\n")
	else:           print_msg("Stop iteration with         : maxit\n")

	if os.path.exists(outdir):
		os.system('rm -rf '+outdir)
	os.mkdir(outdir)

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
	data = [[] for i in xrange(NG)]
	all_data = []
	if(CTF):
		ctf_params = ima.get_attr( "ctf" )
		data_had_ctf = ima.get_attr( "ctf_applied" )
		ctm = ctf_2(nx, ctf_params)
		lctf = len(ctm)
		ctf2 = [[[0.0]*lctf for j in xrange(2)] for i in xrange(NG)]
		ctfb2 = [[0.0]*lctf for i in xrange(NG)]
	for im in xrange(nima):
		k = im%NG
		ima = EMData()
		ima.read_image(stack, im)
		if(CTF):
			ctf_params = ima.get_attr( "ctf" )
			ctm = ctf_2(nx, ctf_params)
			kl = (im//2)%NG  # not sure it will work for NG>2
			for i in xrange(lctf):
				ctf2[k][kl][i] += ctm[i]
			if(ima.get_attr("ctf_applied") == 0):
				st = Util.infomask(ima, mask, False)
				ima -= st[0]
				#      filt_ctf(image,	     dz,		  cs,		   voltage,	     pixel_size,     amp_contrast=0.1,	  b_factor=0.0):
				ima = filt_ctf(ima, ctf_params)
				ima.set_attr('ctf_applied', 1)
		data[k].append(ima)

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

	total_ave = model_blank(nx,nx)
	if(CTF):
		for k in xrange(NG):
			Util.add_img(total_ave, ave_series_ctf(data[k], ctf_tot))
	else:
		for k in xrange(NG):
			Util.add_img(total_ave, ave_series(data[k])*len(data[k]))
		total_ave /= nima
	drop_image(total_ave, os.path.join(outdir, "total_ave_%03d.hdf"%(0)))
	a0 = total_ave.cmp("dot", total_ave, dict(negative = 0, mask = mask))
	msg = "Initial criterion = %-20.7e\n"%(a0)
	print_msg(msg)
	params = ["alpha", "sx", "sy", "mirror"]
	# initialize data for the reference preparation function
	#  mask can be modified in user_function
	ref_data = []
	ref_data.append( mask )
	ref_data.append( center )
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
				ali2d_single_iter(data[k], numr, wr, cs[k], tavg[k], cnx, cny, xrng[N_step], yrng[N_step], step[N_step], mode, range(len(data[k])))
				av1, av2 = add_oe_series(data[k])
				if(CTF):
					tavg[k] = filt_table(Util.addn_img(av1, av2), ctfb2[k])
					av1    = filt_table(av1, ctf2[k][0])
					av2    = filt_table(av2, ctf2[k][1])
				else:
					tavg[k] = (av1+av2)/len(data[k])
				drop_image(tavg[k], os.path.join(outdir, "aqc_%03d_%03d.hdf"%(k, total_iter)))

				frsc.append(fsc_mask(av1, av2, ref_data[0], 1.0, os.path.join(outdir, "drc_%03d_%03d"%(k, total_iter))))
				#  prepare averages for alignment
				kref_data = []
				kref_data.append( mask )
				kref_data.append( 0 )
				kref_data.append( tavg[k] )
				kref_data.append( frsc[k] )
				#  call the user-supplied function to prepare reference image, i.e., filter it, but do not center!
				ktavg[k], cs[k] = user_func( kref_data )
				del kref_data
			#  This should be done only for estimation of resolution, nothing else!
			alpha, sx, sy, mirror, peak = align2d(ktavg[0], ktavg[1], xrng[0], yrng[0], step=0.25, first_ring = first_ring, last_ring = last_ring, rstep=1, mode = mode)
			#  apply parameters to the original average
			favg2 = rot_shift2D(tavg[0], alpha, sx, sy, mirrorm interpolation_method="gridding")
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
				ref_data.append( tavg[k] )
				ref_data.append( fscross  )
				#  call the user-supplied function to prepare reference image, i.e., filter and center!
				tavg[k], cs[k] = user_func( ref_data )
				del ref_data[2]
				del ref_data[2]
				drop_image(tavg[k], os.path.join(outdir, "aqf_%03d_%03d.hdf"%(k, total_iter)))
			total_ave = model_blank(nx,nx)
			if(CTF):
				for k in xrange(NG):
					Util.add_img(total_ave, ave_series_ctf(data[k], ctf_tot))
			else:
				for k in xrange(NG):
					Util.add_img(total_ave, ave_series(data[k])*len(data[k]))
				total_ave /= nima
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
	for im in xrange(nima):
		k=im%NG
		imm = im//NG
		data[k][imm].write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)
	print_end_msg("ali2d_cross_res")
'''



# I'll put everything used in alignment based on SSNR here.
def err_func(args, data):
	from utilities import compose_transform2
	from math import sin, pi, sqrt
	
	v1 = data[0]
	v2 = data[1]
	nima = data[2]
	mirror_same = data[3]
	mirror_flip = data[4]

	err = 0.0
	alpha2 = args[0]
	sx2 = args[1]
	sy2 = args[2]
	for im in xrange(nima):
		mirror = v1[im*4+3]
		mirror_i = v2[im*4+3]
		if mirror == mirror_i:
			alpha = v1[im*4]
			sx = v1[im*4+1]
			sy = v1[im*4+2]
			alpha_i2 = v2[im*4]
			sx_i2 = v2[im*4+1]
			sy_i2 = v2[im*4+2]
			sign = 1-2*mirror
			if mirror_flip: sign = -sign
			alpha3, sx3, sy3, scale3 = compose_transform2(alpha*sign, sx*sign, sy, 1.0, alpha2, sx2, sy2, 1.0)
			this_err = 50*abs(sin((alpha3-alpha_i2)/180.0*pi/2))+sqrt((sx3-sx_i2)**2+(sy3-sy_i2)**2) 
			err += this_err

	return -err/mirror_same

def ave_ali_err(stack):
	from utilities import inverse_transform2, amoeba, get_params2D
	from math import sqrt
	from random import random

	nima = EMUtil.get_image_count(stack)
	v1 = []
	v2 = []
	mirror_same = 0 
	for im in xrange(nima):
		img = EMData()
		img.read_image(stack, im)
		alpha, sx, sy, mirror, dummy = get_params2D(img)
		v1.append(alpha)
		v1.append(sx)
		v1.append(sy)
		v1.append(mirror)
		alpha_i, sx_i, sy_i, mirror_i, dummy = get_params2D(img, "xform.align2d_ideal")
		if mirror==mirror_i:	mirror_same += 1
		alpha_i2, sx_i2, sy_i2, scale2 = inverse_transform2(alpha_i, sx_i, sy_i)
		v2.append(alpha_i2)
		v2.append(sx_i2)
		v2.append(sy_i2)
		v2.append(mirror_i)

	same_rate = float(mirror_same)/nima
	mirror_flip = False
	if same_rate < 0.5:
		mirror_flip = True
		for im in xrange(nima):
			v2[im*4+3]=1-v2[im*4+3]
		mirror_same = nima-mirror_same
		same_rate = 1-same_rate
	
	data = []
	data.append(v1)
	data.append(v2)
	data.append(nima)
	data.append(mirror_same)
	data.append(mirror_flip)

	min_err = 1e10
	trial = 20
	for i in xrange(trial):
		alpha = random()*360.0
		sx = random()*6.0-3.0
		sy = random()*6.0-3.0
		if i==0: args = [0.0, 0.0, 0.0]
		else: args = [alpha, sx, sy]
		ps = amoeba(args, [1.0, 1.0, 1.0], err_func, 1.e-5, 1.e-5, 1000, data)
		err = -ps[1]
		if err < min_err:
			min_err = err
			min_ali = ps[0]
	return min_err
	

def ave_ali_err_MPI(stack):
	from utilities import inverse_transform2, amoeba, get_params2D
	from math import sqrt
	from utilities import reduce_EMData_to_root, print_msg
	from mpi import mpi_comm_size, mpi_comm_rank, mpi_bcast, mpi_reduce, MPI_COMM_WORLD, MPI_FLOAT, MPI_MIN
	from random import random

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	nima = EMUtil.get_image_count(stack)
	v1 = []
	v2 = []
	mirror_same = 0 
	img = EMData()
	for im in xrange(nima):
		img.read_image(stack, im)
		alpha, sx, sy, mirror, dummy = get_params2D(img)
		v1.append(alpha)
		v1.append(sx)
		v1.append(sy)
		v1.append(mirror)
		alpha_i, sx_i, sy_i, mirror_i, dummy = get_params2D(img, "xform.align2d_ideal")
		if mirror==mirror_i:	mirror_same += 1
		alpha_i2, sx_i2, sy_i2, scale2 = inverse_transform2(alpha_i, sx_i, sy_i)
		v2.append(alpha_i2)
		v2.append(sx_i2)
		v2.append(sy_i2)
		v2.append(mirror_i)
	
	same_rate = float(mirror_same)/nima
	mirror_flip = False
	if same_rate < 0.5:
		mirror_flip = True
		for im in xrange(nima):
			v2[im*4+3]=1-v2[im*4+3]
		mirror_same = nima-mirror_same
		same_rate = 1-same_rate
	
	data = []
	data.append(v1)
	data.append(v2)
	data.append(nima)
	data.append(mirror_same)
	data.append(mirror_flip)

	ic = 0
	min_err = 1e10
	trial = max(20/number_of_proc, 1)
	for i in xrange(trial):
		alpha = random()*360.0
		sx = random()*6.0-3.0
		sy = random()*6.0-3.0
		if i==0 and myid==0: args = [0.0, 0.0, 0.0]
		else: args = [alpha, sx, sy]
		ps = amoeba(args, [1.0, 1.0, 1.0], err_func, 1.e-5, 1.e-5, 1000, data)
		err = -ps[1]
		if err < min_err:	min_err = err
	min_err = mpi_reduce(min_err, 1, MPI_FLOAT, MPI_MIN, main_node, MPI_COMM_WORLD)
	if myid==main_node:	min_err = min_err[0]
	min_err = mpi_bcast(min_err, 1, MPI_FLOAT, main_node, MPI_COMM_WORLD)
	min_err = min_err[0]
	return min_err
	
	
"""
def gridrot_shift2D(image, kb, ang = 0.0, sx = 0.0, sy = 0.0):
	from math import pi, sin, cos
	
	N = image.get_ysize()
	image1 = image.copy()
	image1 = image1.fouriergridrot2d(ang, 1.0, kb)
	image1.center_origin_fft()
	for x in xrange((N+2)/2):
		for y in xrange(N):
			if y > N/2-1: yy = y-N
			else: yy = y
			Re = image1.get_value_at(x*2,y)
			Im = image1.get_value_at(x*2+1,y)
			ang = -2*pi*(sx*x+sy*yy)/N
			image1.set_value_at(x*2, y, Re*cos(ang)-Im*sin(ang)) 
			image1.set_value_at(x*2+1, y, Re*sin(ang)+Im*cos(ang))
	return image1
"""
	
def SSNR_func(args, data):
	from utilities import print_msg
	
	img_data = data[0]
	nima = data[1]
	maskI = data[2]
	kb = data[3]
	CTF = data[6]
	SSNR_fit = data[7]
	if CTF:
		index_list = data[8]
		ctfimg_list = data[9]
		ctfimg2 = data[10]
		ctf_data = 3
	else:
		ctf_data = 0
	if SSNR_fit:
		SSNR_img = data[8+ctf_data]
		SSNR_sig = data[9+ctf_data]
		
	N = img_data[0].get_ysize()
	avgimg  = EMData(N, N, 1, False)
	var     = EMData(N, N, 1, False)
	
	for im in xrange(nima):
		imgft = img_data[im].copy()
		alpha = args[im*3]
		sx = args[im*3+1]
		sy = args[im*3+2]
		#imgft = gridrot_shift2D(imgft, kb, alpha, sx, sy)
		imgft = imgft.fouriergridrot_shift2d(alpha, sx, sy, kb)
		imgft.center_origin_fft()
		Util.add_img2(var, imgft)
		if CTF:
			ctfimg = ctfimg_list[index_list[im]]
			Util.mul_img(imgft, ctfimg)
		Util.add_img(avgimg, imgft)

		#print "%2d %20.12e %20.12e"%(im, imgft.get_value_at(2, 0), imgft.get_value_at(3, 0))

	'''
	yyy1 = avgimg.get_value_at(2, 0)
	yyy2 = avgimg.get_value_at(3, 0)
	yyy3 = var.get_value_at(2, 0)
	print yyy1, yyy2, yyy3
	yyy1 /= nima
	yyy2 /= nima
	yyy = yyy1**2+yyy2**2
	print yyy*nima, (yyy3-nima*yyy)/(nima-1), yyy*nima/((yyy3-nima*yyy)/(nima-1))
	'''
	
	if CTF:
		Util.div_filter(avgimg, ctfimg2)
	else:
		Util.mul_scalar(avgimg, 1.0/float(nima))		
	avgimg_2 = avgimg.copy()
	Util.mul_img(avgimg_2, avgimg_2.conjg())
	if CTF:
		Util.mul_img(avgimg_2, ctfimg2)
	else:
		Util.mul_scalar(avgimg_2, float(nima))
	Util.sub_img(var, avgimg_2)
	Util.mul_scalar(var, 1.0/float(nima-1))
	
	SSNR = avgimg_2.copy()
	Util.div_filter(SSNR, var)
	
	#print avgimg_2.get_value_at(2, 0), var.get_value_at(2, 0), SSNR.get_value_at(2, 0)
	#print " "
	
	a0 = Util.infomask(SSNR, maskI, True)
	sum_SSNR = a0[0]
	print_msg("SSNR = %20.7f\n"%(sum_SSNR))
	
	if SSNR_fit:
		beta = 10.0
		Util.sub_img(SSNR, maskI)
		SSNR = SSNR.process("threshold.belowtozero", {"minval": 0.0})
		Util.sub_img(SSNR, SSNR_img)
		Util.div_filter(SSNR, SSNR_sig)
		Util.mul_img(SSNR, SSNR)
		a0 = Util.infomask(SSNR, maskI, True)
		sum_SSNR -= beta*a0[0]
	
	return -sum_SSNR


def SSNR_grad(args, data):
	from numpy import zeros, array, float64

	img_data = data[0]
	nima = data[1]
	maskI = data[2]
	kb = data[3]
	_jX = data[4]
	_jY = data[5]
	CTF = data[6]
	SSNR_fit = data[7]
	if CTF:
		index_list = data[8]
		ctfimg_list = data[9]
		ctfimg2 = data[10]
		ctf_data = 3
	else:
		ctf_data = 0
	if SSNR_fit:
		SSNR_img = data[8+ctf_data]
		SSNR_sig = data[9+ctf_data]
	
	N = img_data[0].get_ysize()
	avgimg  = EMData(N, N, 1, False)
	var     = EMData(N, N, 1, False)
	img_data_new = []
	d_img = []

	for im in xrange(nima):
		imgft = img_data[im].copy()
		alpha = args[im*3]
		sx = args[im*3+1]
		sy = args[im*3+2]
		
		dalpha = 0.05
		#imgft0 = gridrot_shift2D(imgft, kb, alpha, sx, sy)
		#imgft1 = gridrot_shift2D(imgft, kb, alpha-dalpha, sx, sy)
		#imgft2 = gridrot_shift2D(imgft, kb, alpha+dalpha, sx, sy)
		imgft0 = imgft.fouriergridrot_shift2d(alpha, sx, sy, kb)
		imgft0.center_origin_fft()
		imgft1 = imgft.fouriergridrot_shift2d(alpha-dalpha, sx, sy, kb)
		imgft1.center_origin_fft()
		imgft2 = imgft.fouriergridrot_shift2d(alpha+dalpha, sx, sy, kb)
		imgft2.center_origin_fft()
		Util.sub_img(imgft2, imgft1)
		Util.mul_scalar(imgft2, 1/(2*dalpha))
				
		img_data_new.append(imgft0)
		d_img.append(imgft2)
		Util.add_img2(var, imgft0)
		if CTF:
			ctfimg = ctfimg_list[index_list[im]]
			Util.mul_img(imgft0, ctfimg)
		Util.add_img(avgimg, imgft0)
		
	if CTF:
		Util.div_filter(avgimg, ctfimg2)
	else:
		Util.mul_scalar(avgimg, 1.0/float(nima))
	avgimg_2 = avgimg.copy()
	Util.mul_img(avgimg_2, avgimg_2.conjg())
	if CTF:
		Util.mul_img(avgimg_2, ctfimg2)
	else:
		Util.mul_scalar(avgimg_2, float(nima))
	sumimg2 = var.copy()
	Util.sub_img(var, avgimg_2)
	Util.mul_scalar(var, 1.0/float(nima-1))
	
	avgimg_conj = avgimg.conjg()
	dSSNR = avgimg_conj.copy()
	Util.div_filter(dSSNR, var)

	if SSNR_fit:
		beta = 10.0
		SSNR = avgimg_2.copy()
		Util.div_filter(SSNR, var)
		Util.sub_img(SSNR, maskI)
		SSNR = SSNR.process("threshold.belowtozero", {"minval": 0.0})
		Util.sub_img(SSNR, SSNR_img)
		Util.div_filter(SSNR, SSNR_sig)
		Util.div_filter(SSNR, SSNR_sig)
		Util.mul_scalar(SSNR, 2*beta)
		C = maskI.copy()
		Util.sub_img(C, SSNR)
		Util.mul_img(dSSNR, C)
	
	g = zeros(args.shape, float64)
	accurate = True

	for im in xrange(nima):
		img_new = img_data_new[im].copy()
		dSSNR_copy = dSSNR.copy()

		if accurate: 
			
			img_new_copy = img_new.copy()
			Util.sub_img(img_new_copy, avgimg)
			Util.mul_img(img_new_copy, dSSNR)
			img_new_copy = img_new_copy.conjg()
			Util.mul_scalar(img_new_copy, nima/float(nima-1))
			C = maskI.copy()
			Util.sub_img(C, img_new_copy)
			Util.mul_img(dSSNR_copy, C)
			'''
			C = sumimg2.copy()
			img_new_congj = img_new.conjg()
			Util.mul_img(img_new_congj, avgimg)
			Util.mul_scalar(img_new_congj, float(nima))
			Util.sub_img(C, img_new_congj)
			Util.div_filter(C, var)
			Util.mul_scalar(C, 1.0/float(nima-1))
			Util.mul_img(dSSNR_copy, C)
			'''
			
		Util.mul_img(dSSNR_copy, d_img[im])
		if CTF:
			Util.mul_img(dSSNR_copy, ctfimg_list[index_list[im]])
		Util.add_img(dSSNR_copy, dSSNR_copy.conjg())
		a0 = Util.infomask(dSSNR_copy, maskI, True)
		g[im*3] = -a0[0]
		
		dSSNR_copy = dSSNR.copy()
		if accurate:
			Util.mul_img(dSSNR_copy, C)
		Util.mul_img(dSSNR_copy, img_new)
		if CTF:
			Util.mul_img(dSSNR_copy, ctfimg_list[index_list[im]])
		dSSNR_fft = dSSNR_copy.copy()		
		Util.mul_img(dSSNR_fft, _jX)
		Util.add_img(dSSNR_fft, dSSNR_fft.conjg())
		
		#if im == 0: 	print dSSNR_fft.get_value_at(2, 0)
		
		a0 = Util.infomask(dSSNR_fft, maskI, True)
		g[im*3+1] = -a0[0]
		
		dSSNR_fft = dSSNR_copy.copy()		
		Util.mul_img(dSSNR_fft, _jY)
		Util.add_img(dSSNR_fft, dSSNR_fft.conjg())
		a0 = Util.infomask(dSSNR_fft, maskI, True)
		g[im*3+2] = -a0[0]
	return g


def SSNR_func_MPI(args, data):
	from applications import MPI_start_end
	from utilities import reduce_EMData_to_root, print_msg
	from mpi import mpi_comm_size, mpi_comm_rank, mpi_bcast, MPI_COMM_WORLD, MPI_FLOAT

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	img_data = data[0]
	nima = data[1]
	maskI = data[2]
	kb = data[3]
	CTF = data[6]
	SSNR_fit = data[7]
	if CTF:
		index_list = data[8]
		ctfimg_list = data[9]
		ctfimg2 = data[10]
		ctf_data = 3
	else:
		ctf_data = 0
	if SSNR_fit:
		"""
		SSNR_img = data[8+ctf_data]
		SSNR_sig = data[9+ctf_data]
		"""
		SSNR_r = data[8+ctf_data]
	
	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
	
	N = img_data[0].get_ysize()
	avgimg  = EMData(N, N, 1, False)
	var     = EMData(N, N, 1, False)

	for im in xrange(image_start, image_end):
		imgft = img_data[im-image_start].copy()
		alpha = args[im*3]
		sx = args[im*3+1]
		sy = args[im*3+2]
		#imgft = gridrot_shift2D(imgft, kb, alpha, sx, sy)			
		imgft = imgft.fouriergridrot_shift2d(alpha, sx, sy, kb)
		imgft.center_origin_fft()
		Util.add_img2(var, imgft)
		if CTF:
			ctfimg = ctfimg_list[index_list[im]]
			Util.mul_img(imgft, ctfimg)
		Util.add_img(avgimg, imgft)
		
	reduce_EMData_to_root(avgimg, myid, main_node)
	reduce_EMData_to_root(var, myid, main_node)
	if myid == main_node:
		if CTF:
			Util.div_filter(avgimg, ctfimg2)
		else:
			Util.mul_scalar(avgimg, 1.0/float(nima))
		avgimg_2 = avgimg.copy()
		Util.mul_img(avgimg_2, avgimg_2.conjg())
		if CTF:
			Util.mul_img(avgimg_2, ctfimg2)
		else:
			Util.mul_scalar(avgimg_2, float(nima))
		Util.sub_img(var, avgimg_2)
		Util.mul_scalar(var, 1.0/float(nima-1))
	
		SSNR = avgimg_2.copy()
		Util.div_filter(SSNR, var)
		
		a0 = Util.infomask(SSNR, maskI, True)
		sum_SSNR = a0[0]
		print_msg("SSNR = %20.7f\n"%(sum_SSNR))

		if SSNR_fit:
			from fundamentals import rot_avg_table
			beta = 10.0
			avgimg_2 = Util.pack_complex_to_real(avgimg_2)
			var   = Util.pack_complex_to_real(var)
			ravgimg_2 = rot_avg_table(avgimg_2)
			rvar = rot_avg_table(var)
			SSNR_diff = 0.0
			for i in xrange(N/2+1):
				qt = max(0.0, ravgimg_2[i]/rvar[i] - 1.0)
				diff = (qt-SSNR_r[i])/max(1.0, SSNR_r[i])
				print "In bin %3d  Target SSNR = %10.4f  Actual SSNR= %10.4f difference in pct is %7.3f"%(i,SSNR_r[i], qt, diff)
				print "Numerator: ", ravgimg_2[i], "Denominator: ", rvar[i]
				SSNR_diff += diff**2
			sum_SSNR -= beta*SSNR_diff				
			"""
			Util.sub_img(SSNR, maskI)
			SSNR = SSNR.process("threshold.belowtozero", {"minval": 0.0})
			Util.sub_img(SSNR, SSNR_img)
			Util.div_filter(SSNR, SSNR_sig)
			Util.mul_img(SSNR, SSNR)
			a0 = Util.infomask(SSNR, maskI, True)
			sum_SSNR -= beta*a0[0]
			"""
	else: 	
		sum_SSNR = 0.0
	
	sum_SSNR = mpi_bcast(sum_SSNR, 1, MPI_FLOAT, main_node, MPI_COMM_WORLD)
	return -sum_SSNR[0]


def SSNR_grad_MPI(args, data):
	from numpy import zeros, array, float32, float64
	from applications import MPI_start_end
	from utilities import reduce_EMData_to_root, bcast_EMData_to_all 
	from utilities import reduce_array_to_root, bcast_array_to_all
	from mpi import mpi_comm_size, mpi_comm_rank, mpi_bcast, MPI_COMM_WORLD, MPI_FLOAT

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	img_data = data[0]
	nima = data[1]
	maskI = data[2]
	kb = data[3]
	_jX = data[4]
	_jY = data[5]
	CTF = data[6]
	SSNR_fit = data[7]
	if CTF:
		index_list = data[8]
		ctfimg_list = data[9]
		ctfimg2 = data[10]
		ctf_data = 3
	else:
		ctf_data = 0
	if SSNR_fit:
		"""
		SSNR_img = data[8+ctf_data]
		SSNR_sig = data[9+ctf_data]
		"""
		SSNR_r = data[8+ctf_data]
	
	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)

	N = img_data[0].get_ysize()
	avgimg  = EMData(N, N, 1, False)
	var     = EMData(N, N, 1, False)
	img_data_new = []
	d_img = []

	for im in xrange(image_start, image_end):
		imgft = img_data[im-image_start].copy()
		alpha = args[im*3]
		sx = args[im*3+1]
		sy = args[im*3+2]
		
		dalpha = 0.05
		#imgft0 = gridrot_shift2D(imgft, kb, alpha, sx, sy)
		#imgft1 = gridrot_shift2D(imgft, kb, alpha-dalpha, sx, sy)
		#imgft2 = gridrot_shift2D(imgft, kb, alpha+dalpha, sx, sy)
		imgft0 = imgft.fouriergridrot_shift2d(alpha, sx, sy, kb)
		imgft0.center_origin_fft()		
		imgft1 = imgft.fouriergridrot_shift2d(alpha-dalpha, sx, sy, kb)
		imgft1.center_origin_fft()
		imgft2 = imgft.fouriergridrot_shift2d(alpha+dalpha, sx, sy, kb)
		imgft2.center_origin_fft()
		Util.sub_img(imgft2, imgft1)
		Util.mul_scalar(imgft2, 1/(2*dalpha))
				
		img_data_new.append(imgft0)
		d_img.append(imgft2)
		Util.add_img2(var, imgft0)
		if CTF:
			ctfimg = ctfimg_list[index_list[im]]
			Util.mul_img(imgft0, ctfimg)
		Util.add_img(avgimg, imgft0)
		
	reduce_EMData_to_root(avgimg, myid, main_node)
	reduce_EMData_to_root(var, myid, main_node)
	bcast_EMData_to_all(avgimg, myid, main_node)
	bcast_EMData_to_all(var, myid, main_node)

	if CTF:
		Util.div_filter(avgimg, ctfimg2)
	else:
		Util.mul_scalar(avgimg, 1.0/float(nima))
	avgimg_2 = avgimg.copy()
	Util.mul_img(avgimg_2, avgimg_2.conjg())
	if CTF:
		Util.mul_img(avgimg_2, ctfimg2)
	else:
		Util.mul_scalar(avgimg_2, float(nima))
	Util.sub_img(var, avgimg_2)
	Util.mul_scalar(var, 1.0/float(nima-1))
	
	avgimg_conj = avgimg.conjg()
	dSSNR = avgimg_conj.copy()
	Util.div_filter(dSSNR, var)

	if SSNR_fit:
		from fundamentals import rot_avg_table
		from math import sqrt, pi
		beta = 10.0
		avgimg_2 = Util.pack_complex_to_real(avgimg_2)
		var   = Util.pack_complex_to_real(var)
		ravgimg_2 = rot_avg_table(avgimg_2)
		rvar = rot_avg_table(var)
		C = EMData(N, N, 1, False)
		S = pi*(0.49*N)**2
		for x in xrange((N+2)/2):
			for y in xrange(N):
 				if y > N/2-1: yy = y-N
				else: yy = y
				r = sqrt(x**2+yy**2)
				if r < 0.49*N:
					i = int(r+0.5)
					qt = max(0.0, ravgimg_2[i]/rvar[i] - 1.0)
					temp = 1-2*beta*(qt-SSNR_r[i])/max(1.0, SSNR_r[i])**2*S/max(1, 2*pi*i)
					C.set_value_at(x*2, y, temp)
		Util.mul_img(dSSNR, C)
		
		"""
		SSNR = avgimg_2.copy()
		Util.div_filter(SSNR, var)
		Util.sub_img(SSNR, maskI)
		SSNR = SSNR.process("threshold.belowtozero", {"minval": 0.0})
		Util.sub_img(SSNR, SSNR_img)
		Util.div_filter(SSNR, SSNR_sig)
		Util.div_filter(SSNR, SSNR_sig)
		Util.mul_scalar(SSNR, 2*beta)
		C = maskI.copy()
		Util.sub_img(C, SSNR)
		Util.mul_img(dSSNR, C)
		"""
	
	h = zeros(args.shape, float32)
	accurate = True
	
	for im in xrange(image_start, image_end):
		img_new = img_data_new[im-image_start].copy()
		dSSNR_copy = dSSNR.copy()

		if accurate: 
			img_new_copy = img_new.copy()
			Util.sub_img(img_new_copy, avgimg)
			Util.mul_img(img_new_copy, dSSNR)
			img_new_copy = img_new_copy.conjg()
			Util.mul_scalar(img_new_copy, nima/float(nima-1))
			C = maskI.copy()
			Util.sub_img(C, img_new_copy)
			Util.mul_img(dSSNR_copy, C)

		Util.mul_img(dSSNR_copy, d_img[im-image_start])
		if CTF:
			Util.mul_img(dSSNR_copy, ctfimg_list[index_list[im]])
		Util.add_img(dSSNR_copy, dSSNR_copy.conjg())
		a0 = Util.infomask(dSSNR_copy, maskI, True)
		h[im*3] = -a0[0]
		
		dSSNR_copy = dSSNR.copy()
		if accurate:
			Util.mul_img(dSSNR_copy, C)
		Util.mul_img(dSSNR_copy, img_new)
		if CTF:
			Util.mul_img(dSSNR_copy, ctfimg_list[index_list[im]])
		
		dSSNR_fft = dSSNR_copy.copy()		
		Util.mul_img(dSSNR_fft, _jX)
		Util.add_img(dSSNR_fft, dSSNR_fft.conjg())
		a0 = Util.infomask(dSSNR_fft, maskI, True)
		h[im*3+1] = -a0[0]
		
		dSSNR_fft = dSSNR_copy.copy()		
		Util.mul_img(dSSNR_fft, _jY)
		Util.add_img(dSSNR_fft, dSSNR_fft.conjg())
		a0 = Util.infomask(dSSNR_fft, maskI, True)
		h[im*3+2] = -a0[0]
		
	reduce_array_to_root(h, myid, main_node)
	bcast_array_to_all(h, myid, main_node)
	
	g = zeros(args.shape, float64)
	g = float64(h)
	return g


def ali_SSNR(stack, maskfile=None, ou=-1, maxit=10, CTF=False, opti_method="CG", SSNR_fit=False, SSNR=[], MPI=False):

	if MPI:
		ali_SSNR_MPI(stack, maskfile, ou, maxit, CTF, opti_method, SSNR_fit, SSNR)
		return
		
	from math import pi, sqrt
	from fundamentals import fftip, mirror
	from numpy import Inf
	from scipy.optimize.lbfgsb import fmin_l_bfgs_b
	from scipy.optimize.optimize import fmin_cg
	from utilities import get_image, get_params2D, set_params2D, print_begin_msg, print_end_msg, print_msg
	
	if CTF:
		from utilities import get_arb_params
		from morphology import ctf_img

	print_begin_msg("ali_SSNR")

	if opti_method!="CG" and opti_method!="LBFGSB":
		 print "Unknown optimization method!"
		 return
		 
	nima = EMUtil.get_image_count(stack)
	ima = EMData()
	ima.read_image(stack, 0)
	nx = ima.get_xsize()
	
	last_ring = int(ou);	max_iter = int(maxit)
	if last_ring == -1:	last_ring = nx//2-2
	
	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("Data with CTF               : %s\n"%(CTF))
	print_msg("Optimization method         : %s\n"%(opti_method))

	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  
			print_msg("Maskfile                    : %s\n\n"%(maskfile))
			mask = get_image(maskfile)
		else:
			print_msg("Maskfile                    : user provided in-core mask\n\n")
			mask = maskfile
	else :
		from utilities import model_circle 
		print_msg("Maskfile                    : default, a circle with radius %i\n\n"%(last_ring))
		mask = model_circle(last_ring, nx, nx)
	
	
	# prepare kb for gridding interpolation
	npad = 2
	N = nx*npad
	K = 6
	alpha = 1.75
	r = nx/2
	v = K/2.0/N
	kb = Util.KaiserBessel(alpha, K, r, v, N)
	del alpha, K, r, v

	# generate the mask in Fourier space
	maskI = EMData(N, N, 1, False)
	for x in xrange((N+2)/2):
		for y in xrange(N):
 			if y > N/2-1: yy = y-N
			else: yy = y
			if x**2+yy**2 < (N*0.49)**2:
				maskI.set_value_at(x*2, y, 1) 
	maskI.set_value_at(0, 0, 0)
	maskI.set_value_at(1, 0, 0)
	
	if CTF:
		defocus_list = []
		ctfimg_list = []
		index_list = []
		ctfimg2 = EMData(N, N, 1, False)

	# There two matrixes in Fourier space is necessary for calculating derivatives
	_jX = EMData(N, N, 1, False)
	_jY = EMData(N, N, 1, False)
	for x in xrange((N+2)/2):
		for y in xrange(N):
 			if y > N/2-1: yy = y-N
			else: yy = y
		 	_jX.set_value_at(x*2+1, y, -x*2*pi/N) 
 			_jY.set_value_at(x*2+1, y, -yy*2*pi/N) 

	img_data= []
	x0 = [0.0]*(nima*3)
	mir = [0]*nima	
	if opti_method == "LBFGSB":	bounds = []

	# pre-processing, get initial parameters and boundaries (for LBFGSB)
	for im in xrange(nima):
		if im>0:
			ima = EMData()
			ima.read_image(stack, im)
		if CTF:
			ctf_params = ima.get_attr('ctf')

		x0[im*3], x0[im*3+1], x0[im*3+2], mir[im], dummy = get_params2D(ima)
		
		if mir[im] == 1:
			x0[im*3] = 360.0-x0[im*3]
			x0[im*3+1] = -x0[im*3+1]
			ima = mirror(ima)

		st = Util.infomask(ima, mask, False)
		ima -= st[0]	
		ima.divkbsinh(kb)
		ima = ima.norm_pad(False, npad)
		fftip(ima)
		ima.center_origin_fft()
		img_data.append(ima)
		if CTF:
			if not (ctf_params.defocus in defocus_list):
				defocus_list.append(ctf_params.defocus)
				ctfimg = ctf_img(N, ctf_params, ny = N, nz = 1)
				ctfimg_list.append(ctfimg)
				index_list.append(len(defocus_list)-1)
			else:
				index = defocus_list.index(ctf_params.defocus)
				ctfimg = ctfimg_list[index]
				index_list.append(index)
			Util.add_img2(ctfimg2, ctfimg)			

		if opti_method == "LBFGSB":
			bounds.append((x0[im*3]-2.0, x0[im*3]+2.0))
			bounds.append((x0[im*3+1]-1.0, x0[im*3+1]+1.0))
			bounds.append((x0[im*3+2]-1.0, x0[im*3+2]+1.0))

	# Use a gradient method here
	data = []
	data.append(img_data)
	data.append(nima)
	data.append(maskI)
	data.append(kb)
	data.append(_jX)
	data.append(_jY)
	data.append(CTF)
	data.append(SSNR_fit)

	if CTF:
		data.append(index_list)
		data.append(ctfimg_list)
		data.append(ctfimg2)
	if SSNR_fit:
		"""
		SSNR_img = EMData(N, N, 1, False)
		SSNR_sig = EMData(N, N, 1, False)
		for x in xrange((N+2)/2):
			for y in xrange(N):
 				if y > N/2-1: yy = y-N
				else: yy = y
				r = sqrt(x**2+yy**2)
				if r < N*0.49:
					i = int(r/2)
					j = r/2-i
					temp = (SSNR[i]*(1-j)+SSNR[i+1]*j)*nima
					SSNR_img.set_value_at(x*2, y, temp) 
					if temp < 1.0: temp = 1.0
					SSNR_sig.set_value_at(x*2, y, temp)					
		data.append(SSNR_img)
		data.append(SSNR_sig)
		"""
		SSNR_r = []
		for i in xrange(N/2+1):
			if i%2==0: SSNR_r.append(SSNR[i/2]*nima)
			else:	SSNR_r.append((SSNR[i/2]+SSNR[i/2+1])*0.5*nima)
		data.append(SSNR_r)
	
	from numpy import array
	x1 = array(x0)
	aaa = SSNR_func(x1, data)
	bbb = SSNR_grad(x1, data)

	ccc = [0.0]*(nima*3)
	for i in xrange(nima*3):
		x1[i] += 0.1
		ccc[i] = SSNR_func(x1, data)
		x1[i] -= 0.1
	t1 = []
	t2 = []
	t3 = []
	for i in xrange(nima*3):
		f1 = ccc[i]-aaa
		f2 = bbb[i]*0.1
		if f1!=0.0: pct = abs(f1-f2)/abs(f1)
		else: pct = 0.0
		if i%3==0: t1.append(pct)
		elif i%3==1: t2.append(pct)
		else: t3.append(pct)
		print i, f1, f2, pct
	t1.sort(), t2.sort(), t3.sort()
	print "Median error = ", t1[nima/2], t2[nima/2], t3[nima/2]
	exit()
	
	if opti_method == "CG":
		ps = fmin_cg(SSNR_func, x0, fprime=SSNR_grad, args=([data]), gtol=1e-3, norm=Inf, epsilon=1e-5,
        	      maxiter=max_iter, full_output=0, disp=0, retall=0, callback=None)		
	else:
		ps, val, d = fmin_l_bfgs_b(SSNR_func, x0, args=[data], fprime=SSNR_grad, bounds=bounds, m=10, 
			factr=1e3, pgtol=1e-4, epsilon=1e-2, iprint=-1, maxfun=max_iter)
	
	# write the result into the header
	for im in xrange(nima):
		ima = EMData()
		ima.read_image(stack, im)
		if mir[im] == 0:
			set_params2D(ima, [ps[im*3], ps[im*3+1], ps[im*3+2], 0, 1.0])
		else:
			set_params2D(ima, [360.0-ps[im*3], -ps[im*3+1], ps[im*3+2], 1, 1.0])
		ima.write_image(stack, im)	

	print_end_msg("ali_SSNR")


def ali_SSNR_MPI(stack, maskfile=None, ou=-1, maxit=10, CTF=False, opti_method="CG", SSNR_fit=False, SSNR=[]):

	from applications import MPI_start_end
	from math import pi, sqrt
	from fundamentals import fftip, mirror
	from numpy import Inf
	from scipy.optimize.lbfgsb import fmin_l_bfgs_b
	from scipy.optimize.optimize import fmin_cg
	from mpi import mpi_comm_size, mpi_comm_rank, mpi_bcast, MPI_COMM_WORLD, mpi_barrier
	from utilities import get_image, set_params2D, get_params2D, print_begin_msg, print_end_msg, print_msg
	from utilities import model_circle
	
	if CTF:
		from utilities import get_arb_params
		from morphology import ctf_img

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	if myid == main_node:	print_begin_msg("ali_SSNR_MPI")

	if opti_method!="CG" and opti_method!="LBFGSB":
		 print "Unknown optimization method!"
		 return
		 
	nima = EMUtil.get_image_count(stack)
	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
	ima = EMData()
	ima.read_image(stack, 0)
	nx = ima.get_xsize()
	
	last_ring = int(ou);	max_iter = int(maxit)
	if last_ring == -1:	last_ring = nx//2-2
	
	if myid == main_node:
		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Data with CTF               : %s\n"%(CTF))
		print_msg("Optimization method         : %s\n"%(opti_method))
		print_msg("Number of processors used   : %d\n"%(number_of_proc))

	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  
			if myid == main_node:		print_msg("Maskfile                    : %s\n\n"%(maskfile))
			mask = get_image(maskfile)
		else:
			if myid == main_node: 		print_msg("Maskfile                    : user provided in-core mask\n\n")
			mask = maskfile
	else : 
		if myid==main_node: 	print_msg("Maskfile                    : default, a circle with radius %i\n\n"%(last_ring))
		mask = model_circle(last_ring, nx, nx)
	
	
	# prepare kb for gridding interpolation
	npad = 2
	N = nx*npad
	K = 6
	alpha = 1.75
	r = nx/2
	v = K/2.0/N
	kb = Util.KaiserBessel(alpha, K, r, v, N)
	del alpha, K, r, v

	# generate the mask in Fourier space
	maskI = EMData(N, N, 1, False)
	for x in xrange((N+2)/2):
		for y in xrange(N):
 			if y > N/2-1: yy = y-N
			else: yy = y
			if x**2+yy**2 < (N*0.49)**2:
				maskI.set_value_at(x*2, y, 1) 
	maskI.set_value_at(0, 0, 0)
	maskI.set_value_at(1, 0, 0)
	
	if CTF:
		defocus_list = []
		ctfimg_list = []
		index_list = []
		ctfimg2 = EMData(N, N, 1, False)

	# There two matrixes in Fourier space is necessary for calculating derivatives
	_jX = EMData(N, N, 1, False)
	_jY = EMData(N, N, 1, False)
	for x in xrange((N+2)/2):
		for y in xrange(N):
 			if y > N/2-1: yy = y-N
			else: yy = y
		 	_jX.set_value_at(x*2+1, y, -x*2*pi/N) 
 			_jY.set_value_at(x*2+1, y, -yy*2*pi/N) 

	img_data= []
	x0 = [0.0]*(nima*3)
	mir = [0]*nima
	if opti_method == "LBFGSB":	bounds = []

	# pre-processing, get initial parameters and boundaries (for LBFGSB)
	for im in xrange(nima):
		if im>0:
			ima = EMData()
			ima.read_image(stack, im)
		if CTF:
			ctf_params = ima.get_attr('ctf')

		x0[im*3], x0[im*3+1], x0[im*3+2], mir[im], dummy = get_params2D(ima)
		
		if mir[im] == 1:
			x0[im*3] = 360-x0[im*3]
			x0[im*3+1] = -x0[im*3+1]
			ima = mirror(ima)		

		if (im >= image_start) and (im < image_end):
			st = Util.infomask(ima, mask, False)
			ima -= st[0]	
			ima.divkbsinh(kb)
			ima = ima.norm_pad(False, npad)
			fftip(ima)
			ima.center_origin_fft()
			img_data.append(ima)
		if CTF:
			if not (ctf_params.defocus in defocus_list):
				defocus_list.append(ctf_params.defocus)
				ctfimg = ctf_img(N, ctf_params, ny = N, nz = 1)
				ctfimg_list.append(ctfimg)
				index_list.append(len(defocus_list)-1)
			else:
				index = defocus_list.index(ctf_params.defocus)
				ctfimg = ctfimg_list[index]
				index_list.append(index)
			Util.add_img2(ctfimg2, ctfimg)			
		
		if opti_method == "LBFGSB":
			bounds.append((x0[im*3]-2.0, x0[im*3]+2.0))
			bounds.append((x0[im*3+1]-1.0, x0[im*3+1]+1.0))
			bounds.append((x0[im*3+2]-1.0, x0[im*3+2]+1.0))

	# Use a gradient method here
	data = []
	data.append(img_data)
	data.append(nima)
	data.append(maskI)
	data.append(kb)
	data.append(_jX)
	data.append(_jY)
	data.append(CTF)
	data.append(SSNR_fit)

	if CTF:
		data.append(index_list)
		data.append(ctfimg_list)
		data.append(ctfimg2)
	if SSNR_fit:
		"""
		SSNR_img = EMData(N, N, 1, False)
		SSNR_sig = EMData(N, N, 1, False)
		for x in xrange((N+2)/2):
			for y in xrange(N):
 				if y > N/2-1: yy = y-N
				else: yy = y
				r = sqrt(x**2+yy**2)
				if r < N*0.49:
					i = int(r/2)
					j = r/2-i
					temp = (SSNR[i]*(1-j)+SSNR[i+1]*j)*nima
					SSNR_img.set_value_at(x*2, y, temp) 
					if temp < 1.0: temp = 1.0
					SSNR_sig.set_value_at(x*2, y, temp)					
		data.append(SSNR_img)
		data.append(SSNR_sig)
		"""
		SSNR_r = []
		for i in xrange(N/2+1):
			if i%2==0: SSNR_r.append(SSNR[i/2]*nima)
			else:	SSNR_r.append((SSNR[i/2]+SSNR[i/2+1])*0.5*nima)
		data.append(SSNR_r)
	
	'''
	from numpy import array
	x0 = array(x0)
	aaa = SSNR_func_MPI(x0, data)
	bbb = SSNR_grad_MPI(x0, data)
	ccc = [0.0]*(len(x0))
	for i in xrange(len(x0)):
		x0[i] += 0.1
		ccc[i] = SSNR_func_MPI(x0, data)
		x0[i] -= 0.1
		mpi_barrier(MPI_COMM_WORLD)
	
	if myid == 0:
		t = 0
		for i in xrange(len(x0)):
			x1 = ccc[i]-aaa
			x2 = bbb[i]*0.1
			if x1!=0.0: pct = abs(x1-x2)/abs(x1)
			else: pct = 0.0
			t += pct	
			print i, x1, x2, pct
		print "Average error = ", t/len(x0)
	exit()
	'''

	if opti_method == "CG":
		ps = fmin_cg(SSNR_func_MPI, x0, fprime=SSNR_grad_MPI, args=([data]), gtol=1e-3, norm=Inf, epsilon=1e-5,
        	      maxiter=max_iter, full_output=0, disp=0, retall=0, callback=None)
	else:		
		ps, val, d = fmin_l_bfgs_b(SSNR_func_MPI, x0, args=[data], fprime=SSNR_grad_MPI, bounds=bounds, m=10, 
			factr=1e3, pgtol=1e-4, epsilon=1e-2, iprint=-1, maxfun=max_iter)
	
	# write the result into the header
	if myid == main_node: 
		for im in xrange(nima):
			ima = EMData()
			ima.read_image(stack, im)
			if mir[im] == 0:
				set_params2D(ima, [ps[im*3], ps[im*3+1], ps[im*3+2], 0, 1.0])
			else:
				set_params2D(ima, [360.0-ps[im*3], -ps[im*3+1], ps[im*3+2], 1, 1.0])
			ima.write_image(stack, im)

	mpi_barrier(MPI_COMM_WORLD)	
	if myid == main_node:	print_end_msg("ali_SSNR_MPI")


###############################################################################################
## GA CLUSTERING ##############################################################################

'''
#-------------------------------------#
| Draft Clustering Genetic Algorithm  |
#-------------------------------------#
'''

def GA_watch(gen, FIT, tstart, msg = None):
	from time      import time, ctime
	from utilities import print_msg

	init = False
	end  = False

	if msg ==  'init':  init = True
	elif msg == 'end':  end = True
	
	if init: f = open('WATCH', 'w')
	else:    f = open('WATCH', 'a')

	max_fit = max(FIT)
	min_fit = min(FIT)
	mean = 0
	for n in range(len(FIT)): mean += FIT[n]
	mean /= len(FIT)
	std  = 0
	for n in xrange(len(FIT)): std += (FIT[n] - mean)**2
	std /= len(FIT)

	if not end:
		#f.write('|gen %5d|  max fit: %11.4e    min fit: %11.4e    mean fit: %11.4e    std fit: %11.4e      %s\n' % (gen, max_fit, min_fit, mean, std, ctime()))
		f.write('|gen %5d|   off-line fit: %11.6e (min)   on-line fit: %11.6e (mean)              %s   %s\n' % (gen, min_fit, mean, ctime(), str(int(time() - tstart)).rjust(8, '0')))
	if gen % 100 == 0 or end:
		#print_msg('|gen %5d|  max fit: %11.4e    min fit: %11.4e    mean fit: %11.4e    std fit: %11.4e      %s\n' % (gen, max_fit, min_fit, mean, std, ctime()))
		print_msg('|gen %5d|   off-line fit: %11.6e (min)   on-line fit: %11.6e (mean)              %s   %s\n' % (gen, min_fit, mean, ctime(), str(int(time() - tstart)).rjust(8, '0')))

	if not end: f.close()

def list_send_im(list, dst):
	from utilities import send_EMData, model_blank
	N    = len(list)
	data = model_blank(N)
	for n in xrange(N): data.set_value_at(n, float(list[n]))
	send_EMData(data, dst, 0)

def list_recv_im(org):
	from utilities import recv_EMData
	data = recv_EMData(org, 0)
	N    = data.get_xsize()
	list = [0] * N
	for n in xrange(N): list[n] = int(data.get_value_at(n))
	return list

# list must be allocated
def list_bcast_im(list, myid, main_node):
	from utilities import bcast_EMData_to_all, model_blank
	N   = len(list)
	data = model_blank(N)
	for n in xrange(N): data.set_value_at(n, float(list[n]))
	bcast_EMData_to_all(data, myid, main_node)
	for n in xrange(N): list[n] = int(data.get_value_at(n))
	return list

def GA(stack, rad, K, sizepop, maxgen, pcross, pmut, pselect, elit = False, POP = None, debug = False, control = False, rand_seed = 10, PART = None, SSE = False, MPI = False):
	from random    import randint, randrange, random, gauss
	from utilities import get_im, model_blank, model_circle
	from utilities import print_begin_msg, print_end_msg, print_msg
	from copy      import deepcopy
	from sys       import exit
	from time      import ctime
	if MPI:
		from statistics_jb import GA_MPI
		GA_MPI(stack, rad, K, sizepop, maxgen, pcross, pmut, pselect, elit, POP, debug, control, rand_seed, PART, SSE)
		return

	seed(rand_seed)

	if POP != None:
		restart = True
		sizepop = len(POP)
	else:
		restart = False
		
	# Informations about images
	image = get_im(stack, 0)
	nx    = image.get_xsize()
	ny    = image.get_ysize()
	N     = EMUtil.get_image_count(stack)
	mask  = model_circle(rad, nx, ny)
	del image

	print_begin_msg('Genetic Algorithm')
	print_msg('Stack name                      : %s\n' % stack)
	print_msg('Number of images                : %i\n' % N)
	print_msg('Number of clusters              : %i\n' % K)
	print_msg('Radius maskfile                 : %d\n' % rad)
	print_msg('Size of population              : %d\n' % sizepop)
	print_msg('Start from new population       : %s\n' % (not restart))
	print_msg('Maximum number of generations   : %d\n' % maxgen)
	print_msg('Crossover rate                  : %f\n' % pcross)
	print_msg('Mutation rate                   : %f\n\n' % pmut)

	# data structure
	if not restart:	POP  = [[] for i in xrange(sizepop)]
	FIT  = [0]    * sizepop
	IM   = [None] * N
	

	# open images
	for n in xrange(N):
		im = get_im(stack, n)
		im = Util.compress_image_mask(im, mask)
		IM[n] = im.copy()

		if n == 0:
			# new params according the mask
			nx    = im.get_xsize()
			ny    = im.get_ysize()
			norm  = nx * ny
			buf   = model_blank(nx,ny)
			maskv = buf.copy()
			maskv.to_one()
	del im

	if not restart:
		# first population (randomly)
		for npop in xrange(sizepop):
			for im in xrange(N):
				POP[npop].append(randint(0, K-1))
	
	# cst for fitness
	cst = buf.copy()
	for n in xrange(N): Util.add_img2(cst, IM[n])
	
	# compute fitness
	for npop in xrange(sizepop):
		tmp        = [0] * K
		ct         = [0] * K
		for k in xrange(K): tmp[k] = buf.copy()

		for im in xrange(N):
			Util.add_img(tmp[POP[npop][im]], IM[im])
			ct[POP[npop][im]] += 1

		val = buf.copy()
		for k in xrange(K):
			if ct[k] != 0: Util.add_img(val, Util.mult_scalar(Util.muln_img(tmp[k], tmp[k]), 1 / float(ct[k])))

		res = Util.subn_img(cst, val)
                ret = Util.infomask(res, maskv, True)
		FIT[npop] = ret[0] / float(norm)

	# loop of generation
	tstart = time()
	for ngen in xrange(maxgen):
		if ngen == 0: GA_watch(ngen, FIT, tstart, 'init')
		else:         GA_watch(ngen, FIT, tstart)

		print 'gen: %6d            %s' % (ngen, ctime())

		# Elitism method
		if elit:
			listfit = []
			for n in xrange(sizepop): listfit.append([FIT[n], n])
			listfit.sort()
			ELIT1     = POP[listfit[0][1]]
			ELIT1_fit = FIT[listfit[0][1]]
			ELIT2     = POP[listfit[1][1]]
			ELIT2_fit = FIT[listfit[1][1]]

		NEWPOP  = [[] for i in xrange(sizepop)]
		NEWFIT  = [0] * sizepop
		
		# select chromosome
		#for nselect in xrange(int(ratio_select * sizepop)):
		for nselect in xrange(int(sizepop/2)):
			
			# Selection Gauss Method
			chr1 = int(min(abs(gauss(0, pselect)), 1.0) * sizepop)
			chr2 = int(min(abs(gauss(0, pselect)), 1.0) * sizepop)
			if chr2 == chr1:    chr2 += 1
			if chr2 >= sizepop: chr2 = 0
			listfit = []
			for n in xrange(sizepop): listfit.append([FIT[n], n])
			listfit.sort()
			chr1 = listfit[chr1][1]
			chr2 = listfit[chr2][1]

					
			# cross-over
			if random() <= pcross:
				s1 = int(random() * N)
								
				child1 = [0] * N
				child2 = [0] * N
				
				child1[0:s1]  = POP[chr1][0:s1]
				child1[s1:N]  = POP[chr2][s1:N]

				child2[0:s1]  = POP[chr2][0:s1]
				child2[s1:N]  = POP[chr1][s1:N]
			else:
				child1 = POP[chr1]
				child2 = POP[chr2]
			

			'''
			# multi-cross-over
			if random() <= pcross:
				nb       = int(random() * N)
				template = []
				child1   = [0] * N
				child2   = [0] * N
				for n in xrange(nb): template.append(int(random() * N))
				for n in xrange(N):
					if n in template:
						child1[n] = POP[chr2][n]
						child2[n] = POP[chr1][n]
					else:
						child1[n] = POP[chr1][n]
						child2[n] = POP[chr2][n]
			else:
				child1 = POP[chr1]
				child2 = POP[chr2]
			'''
								
			
			# mutation child1
			for n in xrange(N):
				if random() <= pmut:
					child1[n] = randrange(0, K)

			# mutation child2
			for n in xrange(N):
				if random() <= pmut:
					child2[n] = randrange(0, K)
								
			# fitness 1
			tmp        = [0] * K
			ct         = [0] * K
			for k in xrange(K): tmp[k] = buf.copy()

			for im in xrange(N):
				Util.add_img(tmp[child1[im]], IM[im])
				ct[child1[im]] += 1

			val = buf.copy()
			for k in xrange(K):
				if ct[k] != 0: Util.add_img(val, Util.mult_scalar(Util.muln_img(tmp[k], tmp[k]), 1 / float(ct[k])))

			res = Util.subn_img(cst, val)
			ret = Util.infomask(res, maskv, True)
			fitchild1 = ret[0] / float(norm)
	 		
			# fitness 2
			tmp        = [0] * K
			ct         = [0] * K
			for k in xrange(K): tmp[k] = buf.copy()

			for im in xrange(N):
				Util.add_img(tmp[child2[im]], IM[im])
				ct[child2[im]] += 1

			val = buf.copy()
			for k in xrange(K):
				if ct[k] != 0: Util.add_img(val, Util.mult_scalar(Util.muln_img(tmp[k], tmp[k]), 1 / float(ct[k])))

			res = Util.subn_img(cst, val)
			ret = Util.infomask(res, maskv, True)
			fitchild2 = ret[0] / float(norm)
		
			# Make new POP
			NEWPOP[2*nselect]     = child1
			NEWPOP[2*nselect + 1] = child2
			NEWFIT[2*nselect]     = fitchild1
			NEWFIT[2*nselect + 1] = fitchild2

		# replace pop by new pop
		POP = deepcopy(NEWPOP)
		del NEWPOP
		FIT = deepcopy(NEWFIT)
		del NEWFIT

		# Elitism method
		if elit:
			listfit = []
			for n in xrange(sizepop): listfit.append([FIT[n], n])
			listfit.sort()
			POP[listfit[-1][1]] = ELIT1
			FIT[listfit[-1][1]] = ELIT1_fit
			POP[listfit[-2][1]] = ELIT2
			FIT[listfit[-2][1]] = ELIT2_fit


	# export result
	best = []
	for n in xrange(sizepop): best.append([FIT[n], n])
	best.sort()
	AVE = [0] * K
	CT  = [0] * K
	for k in xrange(K):
		AVE[k] = buf.copy()
		CT[k]  = 0

	listim = [[] for i in xrange(K)]
	for im in xrange(N):
		Util.add_img(AVE[POP[best[0][1]][im]], IM[im])
		listim[POP[best[0][1]][im]].append(im)
		CT[POP[best[0][1]][im]] += 1
	for k in xrange(K):
		if CT[k] > 0: Util.mul_scalar(AVE[k], 1.0 / float(CT[k]))

	for k in xrange(K):
		AVE[k] = Util.reconstitute_image_mask(AVE[k], mask)
		AVE[k].set_attr('members', listim[k])
		AVE[k].set_attr('nobjects', CT[k])
		AVE[k].write_image('average.hdf', k)

	# Save POP
	import pickle
	f = open('Pop_gen%d' % (ngen+1), 'w')
	pickle.dump(POP, f)
	f.close()

	GA_watch(ngen, FIT, tstart, 'end')

	print_msg('\n')
	print_end_msg('Genetic Algorithm')


def GA_MPI(stack, rad, K, sizepop, maxgen, pcross, pmut, pselect, elit = False, gPOP = None, debug = False, control = False, rand_seed = 10, gPART = None, SSE = False):
	from random     import randint, randrange, random, gauss, seed
	from utilities  import get_im, model_blank, model_circle
	from utilities  import print_begin_msg, print_end_msg, print_msg
	from utilities  import bcast_number_to_all, get_im
	from time       import ctime, time
	from mpi 	import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	import mpi_reduce, mpi_bcast, mpi_barrier, mpi_recv, mpi_send
	from mpi 	import MPI_SUM, MPI_FLOAT, MPI_INT, MPI_LOR
	import sys

	# ------------------- MPI init ----------------------------- #
	# init
	sys.argv       = mpi_init(len(sys.argv),sys.argv)
	nb_nodes       = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	seed(rand_seed)

	# choose a random node as a main one
	main_node = 0
	if myid  == 0:	main_node = randrange(0, nb_nodes)
	main_node = bcast_number_to_all(main_node,0)
	mpi_barrier(MPI_COMM_WORLD)

	# Informations about images
	image = get_im(stack, 0)
	nx    = image.get_xsize()
	ny    = image.get_ysize()
	N     = EMUtil.get_image_count(stack)
	mask  = model_circle(rad, nx, ny)
	del image

	# Start from new POP?
	if gPOP is None:
		restart   = False
	elif isinstance(gPOP, basestring):
		import pickle
		try:
			f   = open(gPOP, 'r')
			POP = pickle.load(f)
			f.close()
			restart = True
			sizepop = len(POP)
			if myid != main_node: del POP
		except IOError:
			if myid == main_node: print '*** ERROR: file %s not exist' % gPOP
			sys.exit()
	else:
		if myid == main_node: print '*** ERROR: error with POP file name'
		sys.exit()

	# Start from partition?
	if gPART is None:
		flagpart = False
	elif isinstance(gPART, basestring):
		try:
			flagpart = True
			gK  = EMUtil.get_image_count(gPART)
			tmp = [[] for i in xrange(gK)]
			im  = EMData()
			gN  = 0
			for k in xrange(gK):
				im.read_image(gPART, k, True)
				tmp[k] = im.get_attr('members')
				gN += len(tmp[k])
			part = [0] * gN
			for k in xrange(gK):
				for i in tmp[k]: part[int(i)] = k
			if myid != main_node: del part
		except:
			if myid == main_node: print '*** ERROR: error to load partition %s' % gPART
			raise
			sys.exit()
	else:
		if myid == main_node: print '*** ERROR: error with partition name'
		sys.exit()
			
	# Check value
	if sizepop % 2 != 0:
		if myid == main_node: print '*** ERROR: Size of population must be even: %d' % sizepop
		sys.exit()
	elif sizepop < 2 * nb_nodes:
		if myid == main_node: print '*** ERROR: The minimum size of population is %d (2 * nb nodes), actually is %d.' % (nb_nodes * 2, sizepop)
		sys.exit()
	if gPOP is not None and gPART is not None:
		if myid == main_node: print '*** ERROR: Can\'t start from population and partition given in same time, choose one.'
		sys.exit()
	if flagpart and gN != N:
		if myid == main_node: print '*** ERROR: partition and stack file don\'t have the same number of images: %d and %d' % (gN, N)
		sys.exit()
	if flagpart and gK != K:
		if myid == main_node: print '*** ERROR: the value of K=%d and the number of groups in the partition is different: %d' % (K, gK)
		sys.exit()
	
	# number of chromosomes pair per Node
	nsubpop = [0] * nb_nodes
	ct      = 0
	for n in xrange(int(sizepop/2)):
		nsubpop[ct] += 1
		ct          += 1
		if ct >= nb_nodes: ct = 0
		
	i_subpop = [[] for i in xrange(nb_nodes)]
	ct       = 0
	for n in xrange(nb_nodes):
		for i in xrange(nsubpop[n]):
			i_subpop[n].append(ct)
			ct += 1

	if myid == main_node:
		print nsubpop
		print i_subpop

	if myid == main_node:
		start_time = time()
		print_begin_msg('Genetic Algorithm MPI')
		print_msg('Stack name                      : %s\n'   % stack)
		print_msg('Number of images                : %i\n'   % N)
		print_msg('Number of clusters              : %i\n'   % K)
		print_msg('Radius maskfile                 : %d\n'   % rad)
		print_msg('Size of population              : %d\n'   % sizepop)
		print_msg('Number of cpus                  : %d\n'   % nb_nodes)
		print_msg('Start from given population     : %s\n'   % gPOP)
		print_msg('Start from given partition      : %s\n'   % gPART)
		print_msg('Elitism method                  : %s\n'   % elit)
		print_msg('Heuristics (k-means SSE)        : %s\n'   % SSE)
		print_msg('Maximum number of generations   : %d\n'   % maxgen)
		print_msg('Crossover rate                  : %f\n'   % pcross)
		print_msg('Mutation rate                   : %f\n\n' % pmut)

	# data structure
	IM = [None] * N
	if myid == main_node:
		if not restart: POP  = [[] for i in xrange(sizepop)]
		FIT  = [0] * sizepop
	
	# [ALL] open images
	for n in xrange(N):
		im = get_im(stack, n)
		im = Util.compress_image_mask(im, mask)
		IM[n] = im.copy()

		if n == 0:
			# new params according the mask
			nx    = im.get_xsize()
			ny    = im.get_ysize()
			norm  = nx * ny
			buf   = model_blank(nx,ny)
			maskv = buf.copy()
			maskv.to_one()
	del im

	# [ALL] cst for compute the fitness
	cst = buf.copy()
	for n in xrange(N): Util.add_img2(cst, IM[n])

	# [MAIN] create first pop
	if myid == main_node:
		# first population (randomly) or start from given population\
		if not restart: # new POP
			if not flagpart: # no partition given
				for npop in xrange(sizepop):
					for im in xrange(N): POP[npop].append(randint(0, K-1))
			else:
				# duplicate the same chromosome
				for npop in xrange(sizepop): POP[npop] = part
				del part
				
		# [MAIN] compute fitness
		for npop in xrange(sizepop):
			tmp = [0] * K
			ct  = [0] * K
			for k in xrange(K): tmp[k] = buf.copy()

			for im in xrange(N):
				Util.add_img(tmp[POP[npop][im]], IM[im])
				ct[POP[npop][im]] += 1

			val = buf.copy()
			for k in xrange(K):
				if ct[k] != 0: Util.add_img(val, Util.mult_scalar(Util.muln_img(tmp[k], tmp[k]), 1 / float(ct[k])))

			res = Util.subn_img(cst, val)
			ret = Util.infomask(res, maskv, True)
			FIT[npop] = ret[0] / float(norm)

	if debug: MPIlogfile_init(myid)

	## TO WORKS
	#if myid == main_node: print FIT[0]
	#sys.exit()
	if control and myid == main_node: mem_ctrl = ''
	

	# [ALL] loop of generation
	for ngen in xrange(maxgen):
		if debug: MPIlogfile_print(myid, '|start gen %6d|_____________________________%s\n' % (ngen, ctime()))

		# [ALL]
		if debug: MPIlogfile_print(myid, ' ' * 20 + 'Selection:\n')
		if control and myid == main_node:
			f = open('params', 'w')
			f.write('pmut:    %f\n' % pmut)
			f.write('pcross:  %f\n' % pcross)
			f.write('pselect: %f\n' % pselect)
			f.write('_____history_____\n')
			f.write(mem_ctrl)
			f.close()
		
		# [MAIN] logfile
		if myid == main_node:
			
			listchr1 = []
			listchr2 = []
			for npop in xrange(int(sizepop/2)): 
				# [MAIN]
				# Selection Gauss Method
				chr1 = int(min(abs(gauss(0, pselect)), 1.0) * sizepop)
				chr2 = int(min(abs(gauss(0, pselect)), 1.0) * sizepop)
				if chr2 == chr1:    chr2 += 1
				if chr2 >= sizepop: chr2 = 0
				listfit = []
				for n in xrange(sizepop): listfit.append([FIT[n], n])
				listfit.sort()
				if debug: MPIlogfile_print(myid, ' ' * 25 + 'for %d: %d %d' % (npop, chr1, chr2))
				chr1 = listfit[chr1][1]
				chr2 = listfit[chr2][1]
	
				listchr1.append(chr1)
				listchr2.append(chr2)
				
				if debug: MPIlogfile_print(myid, ' #%d #%d\n' % (chr1, chr2))

			# Elitism method
			if elit:
				ELIT1     = POP[listfit[0][1]]
				ELIT1_fit = FIT[listfit[0][1]]
				ELIT2     = POP[listfit[1][1]]
				ELIT2_fit = FIT[listfit[1][1]]
				

		mpi_barrier(MPI_COMM_WORLD)

		# [MAIN] send parents to the others node
		for node in xrange(nb_nodes):
			if myid == main_node and node != main_node:
				for i in i_subpop[node]:
					list_send_im(POP[listchr1[i]], node)
					list_send_im(POP[listchr2[i]], node)
					if debug: MPIlogfile_print(myid, ' ' * 25 + 'parents send to %d\n' % node)

			if myid == node:
				SUBPOP1 = [[] for i in xrange(nsubpop[node])]
				SUBPOP2 = [[] for i in xrange(nsubpop[node])]
				if myid != main_node:
					for i in xrange(nsubpop[node]):
						SUBPOP1[i] = list_recv_im(main_node)
						SUBPOP2[i] = list_recv_im(main_node)
						
						if debug: MPIlogfile_print(myid, ' ' * 25 + 'parents recv from %d\n' % main_node)
				else:
					for i in xrange(nsubpop[node]):
						SUBPOP1[i] = POP[listchr1[node]]
						SUBPOP2[i] = POP[listchr2[node]]
						if debug: MPIlogfile_print(myid, ' ' * 25 + 'parents recv from me\n')

			mpi_barrier(MPI_COMM_WORLD)

		# [MAIN] init new POP
		if myid == main_node:
			POP  = [[0] * N for i in xrange(sizepop)]
			FIT  = [0]  * sizepop

		# loop of subpop for each node
		SUBFIT1 = [0.0] * nsubpop[myid]
		SUBFIT2 = [0.0] * nsubpop[myid]
		for nsub in xrange(nsubpop[myid]):
			parent1 = SUBPOP1[nsub]
			parent2 = SUBPOP2[nsub]

			# [ID] cross-over
			if random() <= pcross:
				s1 = int(random() * N)

				child1 = [0] * N
				child2 = [0] * N

				child1[0:s1] = parent1[0:s1]
				child1[s1:N] = parent2[s1:N]

				child2[0:s1] = parent2[0:s1]
				child2[s1:N] = parent1[s1:N]

				if debug: MPIlogfile_print(myid, ' ' * 20 + 'Cross-over: cross %d\n' % s1)
			else:
				child1 = parent1
				child2 = parent2

				if debug: MPIlogfile_print(myid, ' ' * 20 + 'Cross-over: no cross\n')

			del parent1, parent2

			# [ID] mutation child1
			ct_mut = 0
			for n in xrange(N):
				if random() <= pmut:
					ct_mut += 1
					child1[n] = randrange(0, K)

			if debug: MPIlogfile_print(myid, ' ' * 20 + 'Mutation: child1 %d\n' % ct_mut)

			# [ID] mutation child2
			ct_mut = 0
			for n in xrange(N):
				if random() <= pmut:
					ct_mut += 1
					child2[n] = randrange(0, K)

			if debug: MPIlogfile_print(myid, ' ' * 20 + 'Mutation: child2 %d\n' % ct_mut)

			# [ID] Heuristics child1, k-means SSE
			#if SSE and ngen % 10 == 0:
			if SSE:
				AVE = [0] * K
				CT  = [0] * K
				for k in xrange(K):     AVE[k] = buf.copy()
				for im in xrange(N):
					Util.add_img(AVE[child1[im]], IM[im])
					CT[child1[im]] += 1
				for k in xrange(K):	AVE[k] = Util.mult_scalar(AVE[k], 1.0/float(CT[k]))
				for im in xrange(N):	child1[im] = Util.min_dist(IM[im], AVE)['pos']

			# [ID] Heuristics child2, k-means SSE
			#if SSE and ngen % 10 == 0:
			if SSE:
				AVE = [0] * K
				CT  = [0] * K
				for k in xrange(K):     AVE[k] = buf.copy()
				for im in xrange(N):
					Util.add_img(AVE[child2[im]], IM[im])
					CT[child2[im]] += 1
				for k in xrange(K):	AVE[k] = Util.mult_scalar(AVE[k], 1.0/float(CT[k]))
				for im in xrange(N):	child2[im] = Util.min_dist(IM[im], AVE)['pos']

			# [ID] compute fitness for the child1
			tmp        = [0] * K
			ct         = [0] * K
			for k in xrange(K): tmp[k] = buf.copy()

			for im in xrange(N):
				Util.add_img(tmp[child1[im]], IM[im])
				ct[child1[im]] += 1

			val = buf.copy()
			for k in xrange(K):
				if ct[k] != 0: Util.add_img(val, Util.mult_scalar(Util.muln_img(tmp[k], tmp[k]), 1 / float(ct[k])))

			res = Util.subn_img(cst, val)
			ret = Util.infomask(res, maskv, True)
			fitchild1 = ret[0] / float(norm)

			if debug: MPIlogfile_print(myid, ' ' * 20 + 'Fitness child1: %11.4e\n' % fitchild1)

			# [ID] compute fitness for the child2
			tmp        = [0] * K
			ct         = [0] * K
			for k in xrange(K): tmp[k] = buf.copy()

			for im in xrange(N):
				Util.add_img(tmp[child2[im]], IM[im])
				ct[child2[im]] += 1

			val = buf.copy()
			for k in xrange(K):
				if ct[k] != 0: Util.add_img(val, Util.mult_scalar(Util.muln_img(tmp[k], tmp[k]), 1 / float(ct[k])))

			res = Util.subn_img(cst, val)
			ret = Util.infomask(res, maskv, True)
			fitchild2 = ret[0] / float(norm)

			if debug: MPIlogfile_print(myid, ' ' * 20 + 'Fitness child2: %11.4e\n' % fitchild2)

			# Replacement
			if debug: MPIlogfile_print(myid, ' ' * 20 + 'Replacement:\n')
			SUBPOP1[nsub] = child1
			SUBPOP2[nsub] = child2
			SUBFIT1[nsub] = fitchild1
			SUBFIT2[nsub] = fitchild2
			del child1, child2

		mpi_barrier(MPI_COMM_WORLD)
		
		# [ID] send SUBPOP to main_node
		for node in xrange(nb_nodes):
			if myid == node and myid != main_node:
				for n in xrange(nsubpop[node]):
					list_send_im(SUBPOP1[n], main_node)
					list_send_im(SUBPOP2[n], main_node)
					mpi_send(SUBFIT1[n], 1, MPI_FLOAT, main_node, 1, MPI_COMM_WORLD)
					mpi_send(SUBFIT2[n], 1, MPI_FLOAT, main_node, 1, MPI_COMM_WORLD)
					if debug: MPIlogfile_print(myid, ' ' * 25 + 'childs send to %d\n' % main_node)

			if myid == main_node:
				if node != main_node:
					for i in i_subpop[node]:
						POP[2*i]     = list_recv_im(node)
						POP[2*i + 1] = list_recv_im(node)
						FIT[2*i]     = mpi_recv(1, MPI_FLOAT, node, 1, MPI_COMM_WORLD).tolist()[0]
						FIT[2*i + 1] = mpi_recv(1, MPI_FLOAT, node, 1, MPI_COMM_WORLD).tolist()[0]
						
						if debug: MPIlogfile_print(myid, ' ' * 25 + 'childs recv from %d\n' % node)
				else:
					ct = 0
					for i in i_subpop[node]:
						POP[2*i]     = SUBPOP1[ct]
						POP[2*i + 1] = SUBPOP2[ct]
						FIT[2*i]     = SUBFIT1[ct]
						FIT[2*i + 1] = SUBFIT2[ct]
						ct += 1

						if debug: MPIlogfile_print(myid, ' ' * 25 + 'childs recv from me\n')

			mpi_barrier(MPI_COMM_WORLD)

		del SUBPOP1, SUBPOP2

		# Elitism method, swap the worse by the elit
		if myid == main_node and elit:
			listfit = []
			for n in xrange(sizepop): listfit.append([FIT[n], n])
			listfit.sort()
			POP[listfit[-1][1]] = ELIT1
			FIT[listfit[-1][1]] = ELIT1_fit
			POP[listfit[-2][1]] = ELIT2
			FIT[listfit[-2][1]] = ELIT2_fit

		# logfile
		if myid == main_node:
			if ngen == 0: GA_watch(ngen, FIT, 'init')
			else:         GA_watch(ngen, FIT)

		# TO WORKS
		if control:
			try:
				cmd = open('control', 'r').readline().strip(' \n').split(':')
				if cmd[0] == 'exit':
					#f = open('control', 'w')
					#f.write('')
					#f.close()
					break
				elif cmd[0] == 'pselect':
					pselect = float(cmd[1])
					savecmd = True
				elif cmd[0] == 'pmut':
					pmut    = float(cmd[1])
					savecmd = True
				elif cmd[0] == 'pcross':
					pcross  = float(cmd[1])
					savecmd = True
				elif cmd[0] == 'sse':
					sse     = int(cmd[1])
					savecmd = True
				else:
					savecmd = False

				#f = open('control', 'w')
				#f.write('')
				#f.close()

				if myid == main_node and savecmd: mem_ctrl += 'gen %d: %s %s\n' % (ngen, cmd[0], cmd[1])
				
			except IOError:
				f   = open('control', 'w')
				f.close()
		

	if myid == main_node:
		GA_watch(ngen, FIT, 'end')
		print_msg('\n')
		time_run = int(time() - start_time)
		time_h   = time_run / 3600
		time_m   = (time_run % 3600) / 60
		time_s   = (time_run % 3600) % 60
		print_msg('Time: %s h %s min %s s\n' % (str(time_h).rjust(2, '0'), str(time_m).rjust(2, '0'), str(time_s).rjust(2, '0')))
	if debug: MPIlogfile_end(myid)

	# [MAIN] export result
	if myid == main_node:
		best = []
		for n in xrange(sizepop): best.append([FIT[n], n])
		best.sort()
		AVE = [0] * K
		CT  = [0] * K
		for k in xrange(K):
			AVE[k] = buf.copy()
			CT[k]  = 0

		listim = [[] for i in xrange(K)]
		for im in xrange(N):
			Util.add_img(AVE[POP[best[0][1]][im]], IM[im])
			listim[POP[best[0][1]][im]].append(im)
			CT[POP[best[0][1]][im]] += 1
		for k in xrange(K):
			if CT[k] > 0: Util.mul_scalar(AVE[k], 1.0 / float(CT[k]))

		for k in xrange(K):
			AVE[k] = Util.reconstitute_image_mask(AVE[k], mask)
			AVE[k].set_attr('members', listim[k])
			AVE[k].set_attr('nobjects', CT[k])
			AVE[k].write_image('average.hdf', k)

		# Save POP
		import pickle
		f = open('Pop_gen%d' % (ngen+1), 'w')
		pickle.dump(POP, f)
		f.close()
		
		print_end_msg('Genetic Algorithm MPI')

	mpi_barrier(MPI_COMM_WORLD)
	return

## END GA CLUSTERING ##########################################################################
###############################################################################################

###############################################################################################
## COMMON LINES OLD CODE ############ 2009-05-11 11:55:18 ########################## JB #######

# transform an image to sinogram (mirror include)
def cml_sinogram(image2D, diameter, d_psi = 1):
	from math         import cos, sin
	from fundamentals import fft
	
	M_PI  = 3.141592653589793238462643383279502884197
	
	# prepare 
	M = image2D.get_xsize()
	# padd two times
	npad  = 2
	N     = M * npad
	# support of the window
	K     = 6
	alpha = 1.75
	r     = M / 2
	v     = K / 2.0 / N

	kb     = Util.KaiserBessel(alpha, K, r, K / (2. * N), N)
	volft  = image2D.average_circ_sub()  	# ASTA - in spider
	volft.divkbsinh(kb)		  	# DIVKB2 - in spider
	volft  = volft.norm_pad(False, npad)
	volft.do_fft_inplace()
	volft.center_origin_fft()
	volft.fft_shuffle()

	# get line projection
	nangle = int(360 / d_psi)     
	dangle = 2 * M_PI / float(nangle)
	data   = []
	for j in xrange(nangle):
		nuxnew =  cos(dangle * j)
		nuynew = -sin(dangle * j)
		line   = volft.extractline(kb, nuxnew, nuynew)
		rlines = fft(line)
		data.append(rlines.copy())

	# copy each line in the same im
	e = EMData()
	e.set_size(data[0].get_xsize() ,len(data), 1)
	for n in xrange(len(data)):
		nx = data[n].get_xsize()
		for i in xrange(nx): e.set_value_at(i, n, data[n].get_value_at(i))

	Util.cyclicshift(e, {"dx":M, "dy":0, "dz":0} )

	return Util.window(e, diameter, len(data), 1, 0, 0, 0)

'''
# NOT USED WILL BE REMOVED jb 2009-03-23 11:41:13
def common_line_in3D(ph1, th1, ph2, th2):
	from math import pi, sqrt, cos, sin, asin, acos

	deg_rad = pi / 180.0
	ph1 *= deg_rad 
	th1 *= deg_rad 
	ph2 *= deg_rad 
	th2 *= deg_rad

	# cross-product between normal vector of projections
	nx = sin(th1)*sin(ph1)*cos(th2) - cos(th1)*sin(th2)*sin(ph2)
	ny = cos(th1)*sin(th2)*cos(ph2) - cos(th2)*sin(th1)*cos(ph1)
	nz = sin(th1)*cos(ph1)*sin(th2)*sin(ph2) - sin(th1)*sin(ph1)*sin(th2)*cos(ph2)

	# normalize
	norm    = nx**2 + ny**2 + nz**2
	rt_norm = sqrt(norm)
	nx /= rt_norm
	ny /= rt_norm
	nz /= rt_norm

	# if theta > 90, apply mirror 
	if nz < 0: nx = -nx; ny = -ny; nz = -nz
	
	# calculate phi and theta (deg)
	thetaCom  = acos(nz)

	if    thetaCom == 0: phiCom = 0
	else:
		val = ny / sin(thetaCom)
		if val > 1.0:  val = 1.0
		if val < -1.0: val = -1.0
		phiCom = asin(val)
	
		phiCom    = (phiCom * 180 / pi + 360)%360
		thetaCom *= (180 / pi)

	return phiCom , thetaCom
'''

'''
# OBSOLETE WILL BE REMOVED jb 2009-03-23 11:40:30
def cml_weights_full(Ori):
	from projection   import common_line_in3D

	# gbl vars
	global g_n_prj, g_n_lines, g_anglst

	
	# gbl vars
	l_phs  = [0.0] * g_n_lines  # angle phi of the common lines
	l_ths  = [0.0] * g_n_lines  # angle theta of the common lines
	n      = 0
	for i in xrange(g_n_prj - 1):
		for j in xrange(i + 1, g_n_prj):
			l_phs[n], l_ths[n] = common_line_in3D(Ori[4*i], Ori[4*i+1], Ori[4*j], Ori[4*j+1])
			n+= 1

	tol = 3

	# search the closer cml lines
	ocp_same   = [-1] * g_n_lines
	num_agl    = 0
	for i in xrange(g_n_lines):
	    if ocp_same[i] == -1:
		ocp_same[i] = num_agl
		for j in xrange(i + 1, g_n_lines):
		    if ocp_same[j] == -1:
			dist = (l_phs[i] - l_phs[j])**2 + (l_ths[i] - l_ths[j])**2
			if dist < tol: ocp_same[j] = num_agl

		num_agl += 1

	if num_agl > 2:

		# create the new vector n_phi n_theta without closer
		n_phi   = [0.0] * num_agl
		n_theta = [0.0] * num_agl
		nb_same = [0]   * num_agl
		num_agl = 0
		for n in xrange(g_n_lines):
		    nb_same[ocp_same[n]] += 1
		    if ocp_same[n] == num_agl:
			n_phi[num_agl]   = l_phs[n]
			n_theta[num_agl] = l_ths[n]
			num_agl += 1

		# Voronoi
		n_weights = Util.vrdg(n_phi, n_theta)

		weights = [0.0] * g_n_lines
		for i in xrange(g_n_lines):
			if nb_same[ocp_same[i]] > 1:
				weights[i] = n_weights[ocp_same[i]] / float(nb_same[ocp_same[i]])
			else:
				weights[i] = n_weights[ocp_same[i]]

	else:
		weights = [6.28/float(g_n_lines)] * g_n_lines

	return weights
'''
	
'''
# OBSOLETE WILL BE REMOVED jb 2009-03-23 11:39:54
# compute the weight of the common lines
def cml_weights_iagl(Ori, iagl, iprj):
	from projection   import common_line_in3D

	# gbl vars
	global g_n_prj, g_n_lines, g_anglst

	
	# gbl vars
	l_phs  = [0.0] * g_n_lines  # angle phi of the common lines
	l_ths  = [0.0] * g_n_lines  # angle theta of the common lines
	n      = 0
	for i in xrange(g_n_prj - 1):
		for j in xrange(i + 1, g_n_prj):
			if i == iprj:   l_phs[n], l_ths[n] = common_line_in3D(g_anglst[iagl][0], g_anglst[iagl][1], Ori[4*j], Ori[4*j+1])
			elif j == iprj:	l_phs[n], l_ths[n] = common_line_in3D(Ori[4*i], Ori[4*i+1], g_anglst[iagl][0], g_anglst[iagl][1])
			else:		l_phs[n], l_ths[n] = common_line_in3D(Ori[4*i], Ori[4*i+1], Ori[4*j], Ori[4*j+1])
			n+= 1

	tol = 3

	# search the closer cml lines
	ocp_same   = [-1] * g_n_lines
	num_agl    = 0
	for i in xrange(g_n_lines):
	    if ocp_same[i] == -1:
		ocp_same[i] = num_agl
		for j in xrange(i + 1, g_n_lines):
		    if ocp_same[j] == -1:
			dist = (l_phs[i] - l_phs[j])**2 + (l_ths[i] - l_ths[j])**2
			#print i, j, dist
			if dist < tol: ocp_same[j] = num_agl

		num_agl += 1

	if num_agl > 2:

		# create the new vector n_phi n_theta without closer
		n_phi   = [0.0] * num_agl
		n_theta = [0.0] * num_agl
		nb_same = [0]   * num_agl
		num_agl = 0
		for n in xrange(g_n_lines):
		    nb_same[ocp_same[n]] += 1
		    if ocp_same[n] == num_agl:
			n_phi[num_agl]   = l_phs[n]
			n_theta[num_agl] = l_ths[n]
			num_agl += 1

		# Voronoi
		n_weights = Util.vrdg(n_phi, n_theta)

		weights = [0.0] * g_n_lines
		for i in xrange(g_n_lines):
			if nb_same[ocp_same[i]] > 1:
				weights[i] = n_weights[ocp_same[i]] / float(nb_same[ocp_same[i]])
			else:
				weights[i] = n_weights[ocp_same[i]]

	else:
		weights = [6.28/float(g_n_lines)] * g_n_lines

	return weights
'''

# open and transform projections
def cml_open_proj(stack, ir, ou, d_psi, lf, hf):
	from projection   import cml_sinogram
	from utilities    import model_circle, get_params_proj, model_blank, get_im
	from fundamentals import fft

	nprj = EMUtil.get_image_count(stack)                # number of projections
	Prj = []                                            # list of projections
	Ori = [-1] * 4 * nprj                              # orientation intial (phi, theta, psi, index) for each projection

	for i in xrange(nprj):
		image = get_im(stack, i)

		# read initial angles if given
		try:	Ori[4*i], Ori[4*i+1], Ori[4*i+2], s2x, s2y = get_params_proj(image)
		except:	pass
		
		if(i == 0):
			nx = image.get_xsize()
			if(ou < 1):
				ou = nx // 2 - 1
			#else:
			#	ou = int(ou) // 2
			#	ou = 2 * ou +1
			diameter = 2 * ou + 1
			diameter = int(diameter)
			mask2D   = model_circle(ou, nx, nx)
			circ     = mask2D.copy()
			if ou > 1:  circ   -= model_circle(ou - 1, nx, nx)
			if ir > 0:  mask2D -= model_circle(ir, nx, nx)

		# normalize under the mask
		[mean_a, sigma, imin, imax] = Util.infomask(image, circ, True)
		image = (image - mean_a) / sigma
		Util.mul_img(image, mask2D)

		# sinogram
		sino = cml_sinogram(image, diameter, d_psi)

		# prepare the cut positions in order to filter (lf: low freq; hf: high freq)
		ihf = min(int(2 * hf * diameter), diameter + (diameter + 1) % 2)
		ihf = ihf + (ihf + 1) % 2    # index ihf must be odd to take the img part
		ilf = max(int(2 * lf * diameter), 0)
		ilf = ilf + ilf % 2          # index ilf must be even to fall in the real part
		bdf = ihf - ilf + 1

		# process lines
		nxe = sino.get_xsize()
		nye = sino.get_ysize()
		prj = model_blank(bdf, nye)
		prj.set_complex(True)
		for li in xrange(nye):

			# get the line li
			line = model_blank(nxe)
			for ci in xrange(nxe): line.set_value_at(ci, 0, sino.get_value_at(ci, li))

			# normalize this line
			[mean_l, sigma_l, imin, imax] = Util.infomask(line, None, True)
			line = (line - mean_l) / sigma_l

			# fft
			line = fft(line)
	
			# filter (cut part of coef)
			ct = 0
			for ci in xrange(ilf, ihf + 1):
				prj.set_value_at(ct, li, line.get_value_at(ci, 0))
				ct += 1
	
		# store the projection
		Prj.append(prj)

	return Prj, Ori


# export result obtain by the function find_struct
def cml_export_struc(stack, outseed, Ori, BDB):
	from projection import plot_angles
	from utilities  import set_params_proj, get_im

	global g_n_prj
	
	pagls = []
	for i in xrange(g_n_prj):
		data = get_im(stack, i)
		p = [Ori[4*i], Ori[4*i+1], Ori[4*i+2], 0.0, 0.0]
		set_params_proj(data, p)
		data.set_attr('active', 1)
		if BDB:	data.write_image(stack, i)
		else:	data.write_image(stack, i)

		# prepare angles to plot
		pagls.append([Ori[4*i], Ori[4*i+1], Ori[4*i+2]])

	# plot angles
	im = plot_angles(pagls)
	if BDB: im.write_image('bdb:%s_plot_agls' % outseed, 0)
	else:   im.write_image(outseed + 'plot_agls.hdf')

# init the global average used for lot of function to cml
def cml_init_global_var(dpsi, delta, nprj, debug):
	from utilities import even_angles
	global g_anglst, g_d_psi, g_n_psi, g_i_prj, g_n_lines, g_n_prj, g_n_anglst, g_debug
	
	g_anglst   = even_angles(delta, 0.0, 179.9, 0.0, 359.9, 'P')
	g_n_anglst = len(g_anglst)
	g_d_psi    = dpsi
	g_n_psi    = int(360 / dpsi)
	g_i_prj    = -1
	g_n_lines  = ((nprj - 1) * nprj) / 2
	g_n_prj    = nprj
	g_debug    = debug
	

# write the head of the logfile
def cml_head_log(stack, outdir, delta, ir, ou, lf, hf, rand_seed, maxit, given):
	from utilities import print_msg

	# call global var
	global g_anglst, g_n_prj, g_d_psi, g_n_anglst
	
	print_msg('Input stack                  : %s\n'     % stack)
	print_msg('Number of projections        : %d\n'     % g_n_prj)
	print_msg('Output directory             : %s\n'     % outdir)
	print_msg('Angular step                 : %5.2f\n'  % delta)
	print_msg('Sinogram angle accuracy      : %5.2f\n'  % g_d_psi)
	print_msg('Inner particle radius        : %5.2f\n'  % ir)	
	print_msg('Outer particle radius        : %5.2f\n'  % ou)
	print_msg('Filter, minimum frequency    : %5.3f\n'  % lf)
	print_msg('Filter, maximum frequency    : %5.3f\n'  % hf)
	print_msg('Random seed                  : %i\n'     % rand_seed)
	print_msg('Number of maximum iterations : %d\n'     % maxit)
	print_msg('Start from given orientations: %s\n'     % given)
	
	
	#print_msg('Number of trials            : %d\n'     % trials)
	#print_msg('Number of cpus              : %i\n'     % ncpu)
	#if refine:
	#	print_msg('Refinement                  : True\n')
	#else:
	#	print_msg('Refinement                  : False\n')
	print_msg('Number of angles             : %i\n\n'   % g_n_anglst)

# write the end of the logfile
def cml_end_log(Ori, disc, disc_nw, ite):
	from utilities import print_msg
	global g_n_prj
	print_msg('\n\n')
	for i in xrange(g_n_prj): print_msg('Projection #%s: phi %10.5f    theta %10.5f    psi %10.5f\n' % (str(i).rjust(3, '0'), Ori[4*i], Ori[4*i+1], Ori[4*i+2]))
	print_msg('\nNumber of iterations: %d\n' % ite)
	print_msg('Discrepancy: %10.3f\n' % abs(disc))
	print_msg('Discrepancy without weigths: %10.3f\n' % abs(disc_nw))

# display the list of angles for each iterations
def cml_export_txtagls(outdir, Ori, disc, title):
	import time
	global g_n_prj, g_i_prj

	angfile = open(outdir + 'angles', 'a')

	angfile.write('|%s|-----------------------------------------------%s---------\n' % (title, time.ctime()))
	for i in xrange(g_n_prj): angfile.write('%10.3f\t%10.3f\t%10.3f\n' % (Ori[4*i], Ori[4*i+1], Ori[4*i+2]))
			
	angfile.write('\nDiscrepancy: %10.3f\n\n' % disc)
	angfile.close()

# export the progress of the find_struc function
def cml_export_progress(outdir, ite, iprj, iagl, psi, disc, cmd):
	infofile = open(outdir + 'progress', 'a')
	global g_anglst

	if cmd == 'progress':
		txt_ite = str(ite).rjust(3, '0')
		txt_i   = str(iprj).rjust(3, '0')
		txt_a   = str(iagl).rjust(3, '0')
		txt     = 'Ite: %s Prj: %s Agls: %s >> Agls (phi, theta, psi): %10.3f %10.3f %10.3f   Disc: %10.7f' % (txt_ite, txt_i, txt_a, g_anglst[iagl][0], g_anglst[iagl][1], psi, disc)
		

	elif cmd == 'choose':
		txt   = 'Ite: %s  Select Agls: %s >> Agls (phi, theta, psi): %10.3f %10.3f %10.3f   Disc: %10.7f\n' % (str(ite).rjust(3, '0'), str(iagl).rjust(3, '0'), g_anglst[iagl][0], g_anglst[iagl][1], psi, disc)

	infofile.write(txt + '\n')
	infofile.close()


'''
# WILL BE REMOVED jb 2009-03-23 11:43:30
# compute the common lines in sino
def get_common_line_angles(phi1, theta1, psi1, phi2, theta2, psi2, nangle, STOP=False):
	from math import fmod
	R1    = Transform({"type":"spider", "phi":phi1, "theta":theta1, "psi":psi1})
	R2    = Transform({"type":"spider", "phi":phi2, "theta":theta2, "psi":psi2})
	R2T   = R2.inverse()
	R2to1 = R1*R2T

	eulerR2to1 = R2to1.get_rotation("spider")
	phiR2to1   = eulerR2to1["phi"]
	thetaR2to1 = eulerR2to1["theta"]
	psiR2to1   = eulerR2to1["psi"]

	alphain1 = fmod(psiR2to1  + 270.0, 360.0)
	alphain2 = fmod(-phiR2to1 + 270.0, 360.0)

	n1 = int(nangle * fmod(alphain1 + 360, 360) / 360.0)
	n2 = int(nangle * fmod(alphain2 + 360, 360) / 360.0)
	
	return n1, n2
'''


# compute discrepancy according the projections and orientations
def cml_disc(Prj, Ori, flag_weights):
	from math        import pi, fmod

	# gbl vars
	global g_n_prj, g_n_psi, g_n_lines

	if flag_weights:
		cml = Util.cml_line_in3d_full(Ori)    # c-code
		weights = Util.cml_weights(cml)       # c-code
		weights = cml_norm_weights(weights)
	else:   weights = [1.0] * g_n_lines

	com = [0] * 2 * g_n_lines

	# compute the common lines
	count = 0
	for i in xrange(g_n_prj - 1):
		for j in xrange(i + 1, g_n_prj):
			[com[count], com[count + 1]] = Util.cml_line_pos(Ori[4*i], Ori[4*i+1], Ori[4*i+2], Ori[4*j], Ori[4*j+1], Ori[4*j+2], g_n_psi)        # c  code
			count += 2

	n = 0
	L_tot = 0.0

	# compute the discrepancy for all sinograms
	for i in xrange(g_n_prj - 1):
		for j in xrange(i + 1, g_n_prj):
			L      = Prj[i].cm_euc(Prj[j], com[n], com[n + 1])
			L_tot += (L * weights[int(n/2)])
			n     += 2

	return L_tot

# cml spin function for one orientation
def cml_spin(Prj, iprj, Ori, iagl, weights, flag):
	import time
	# gbl vars
	global g_n_prj, g_n_psi, g_n_lines, g_anglst
	com = [0] * 2 * g_n_lines

	#t1 = time.time()
	# compute the common line (only for iprj)
	com = Util.cml_list_line_pos(Ori, g_anglst[iagl][0], g_anglst[iagl][1], iprj, g_n_prj, g_n_psi, g_n_lines)
	#print 'line pos', time.time() - t1, 's'
	
	#t2 = time.time()
	res = Util.cml_spin(g_n_psi, iprj, g_n_prj, weights, com, Prj, flag)
	#print 'spin', time.time() - t2, 's'

	return res[0], int(res[1])

def cml_norm_weights(w):
	wm = max(w)
	nw = []
	for i in xrange(len(w)): nw.append(wm - w[i])
	sw = sum(nw)
	for i in xrange(len(w)): nw[i] /= sw
	return nw

# find structure
def cml_find_structure(Prj, Ori, outdir, maxit, first_zero, flag_weights):
	from projection import cml_spin, cml_export_progress
	import time
	import sys
	
	# global vars
	global g_i_prj, g_n_prj, g_n_anglst, g_anglst, g_d_psi, g_debug, g_n_lines

	# list of free orientation
	ocp = [-1] * g_n_anglst

	if first_zero:
		listprj = range(1, g_n_prj)
		ocp[0]  = 0 
	else:   listprj = range(g_n_prj)

	# iteration loop
	for ite in xrange(maxit):
		t_start = time.time()

		# loop over i prj
		change = False
		for iprj in listprj:

			# Store current index of angles assign
			cur_agl      = Ori[4*iprj+3]
			ocp[cur_agl] = -1

			# loop over all angles
			best_disc = 1e20
			best_psi  = -1
			best_iagl = -1
			#t1 = time.time()
			#c  = 0
			
	
			for iagl in xrange(g_n_anglst):
				# if agls free
				if ocp[iagl] == -1:
					#t1 = time.time()

					# weights
					if flag_weights:
						cml = Util.cml_line_in3d_iagl(Ori, g_anglst[iagl][0], g_anglst[iagl][1], iprj)   # c-code
						#print cml
						weights = Util.cml_weights(cml)                                                  # c-code
						weights = cml_norm_weights(weights)
					else:   weights = [1.0] * g_n_lines

					#print weights
					#print 'get weights', time.time() - t1

					#t2 = time.time()
					disc, ind_psi = cml_spin(Prj, iprj, Ori, iagl, weights, 0)
					#print 'spin', time.time() - t2
				
					# select the best
					if disc < best_disc:
						best_disc = disc
						best_psi  = ind_psi
						best_iagl = iagl

					#c += 1
					if g_debug: cml_export_progress(outdir, ite, iprj, iagl, ind_psi * g_d_psi, disc, 'progress')
				else:
					if g_debug: cml_export_progress(outdir, ite, iprj, iagl, -1, -1, 'progress')

				#if iprj == 1: print disc

			# if change, assign
			if best_iagl != cur_agl:
				ocp[best_iagl] = iprj
				Ori[4*iprj]    = g_anglst[best_iagl][0] # phi
				Ori[4*iprj+1]  = g_anglst[best_iagl][1] # theta
				Ori[4*iprj+2]  = best_psi * g_d_psi     # psi
				Ori[4*iprj+3]  = best_iagl              # index

				change = True
			else:
				ocp[cur_agl]   = iprj

			#if iprj == 1: sys.exit()

			if g_debug: cml_export_progress(outdir, ite, iprj, best_iagl, best_psi * g_d_psi, best_disc, 'choose')

		print 'time', time.time() - t_start, 's'

		# if one change, compute new full disc
		disc = cml_disc(Prj, Ori, flag_weights)

		if g_debug:
			disc2 = cml_disc(Prj, Ori, False)
			print 'Ite: ', disc, '           %6.2f s     ' % (time.time() - t_start), disc2

		# display in the progress file
		cml_export_txtagls(outdir, Ori, disc, 'Ite: %s' % str(ite + 1).rjust(3, '0'))

		if not change: break

	return Ori, disc, ite

# application find structure
def find_struct(stack, out_seedname, ir, ou, delta, dpsi, lf, hf, rand_seed, maxit, given = False, first_zero = False, flag_weights = False, debug = False):
	from utilities import print_begin_msg, print_msg, print_end_msg, start_time, running_time
	import time
	import sys

	# logfile
	t_start = start_time()
	print_begin_msg('find_struct')

	if out_seedname.split(':')[0] == 'bdb':
		BDB = True
		outdir = ''
		out_seedname = out_seedname[4:]
	else:
		BDB = False
		import os
		if os.path.exists(out_seedname):  os.system('rm -rf ' + out_seedname)
		os.mkdir(out_seedname)
		out_seedname += '/'
		outdir = out_seedname

	# import
	from projection  import cml_open_proj, cml_init_global_var, cml_head_log
	from projection  import cml_disc, cml_export_txtagls, cml_export_struc
	from projection  import cml_end_log, cml_find_structure
	from random      import seed, random

	# Open and transform projections
	Prj, Ori = cml_open_proj(stack, ir, ou, dpsi, lf, hf)
	
	# if not angles given select randomly orientation for each projection
	if not given:
		if rand_seed > 0: seed(rand_seed)
		else:             seed()
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

	# Init the global vars
	cml_init_global_var(dpsi, delta, len(Prj), debug)
	
	# Update logfile
	cml_head_log(stack, outdir, delta, ir, ou, lf, hf, rand_seed, maxit, given)

	# Compute the first disc
	disc = cml_disc(Prj, Ori, flag_weights)

	# Update progress file
	cml_export_txtagls(outdir, Ori, disc, 'Init')

	# Find structure
	Ori, disc, ite = cml_find_structure(Prj, Ori, outdir, maxit, first_zero, flag_weights)

	# Export structure
	cml_export_struc(stack, out_seedname, Ori, BDB)

	# Compute disc without weights
	disc_now = cml_disc(Prj, Ori, False)

	# Update logfile
	cml_end_log(Ori, disc, disc_now, ite)
	running_time(t_start)
	print_end_msg('find_struct')

	return disc_now



## END CML OLD CODE ###########################################################################
###############################################################################################

###############################################################################################
## COMMON LINES NEW VERSION ############ 2009-03-26 14:06:22 ####################### JB #######


# export result obtain by the function find_struct MPI GA
def cml2_export_struc_GA(stack, outdir, Ori, igen):
	from projection import plot_angles
	from utilities  import set_params_proj, get_im

	#global g_n_prj
	g_n_prj = len(Ori) // 4
	
	pagls = []
	for i in xrange(g_n_prj):
		data = get_im(stack, i)
		p = [Ori[4*i], Ori[4*i+1], Ori[4*i+2], 0.0, 0.0]
		set_params_proj(data, p)
		data.set_attr('active', 1)
		data.write_image(outdir + '/structure_gen_%03i.hdf' % igen, i)

		# prepare angles to plot
		pagls.append([Ori[4*i], Ori[4*i+1], Ori[4*i+2]])

	# plot angles
	im = plot_angles(pagls)
	im.write_image(outdir + '/plot_agls_gen_%03i.hdf' % igen)

# given two set of orientations (structure), determine the error between them after 3D LSQ registration
def cml2_error_ori(Ori1, Ori2):
	from math  import sin, cos, pi, sqrt
	from numpy import array, linalg, matrix

	rad2deg = 180 / pi
	deg2rad = 1 / rad2deg
	nori    = len(Ori1) // 4

	U1 = []
	for n in xrange(nori):
		ind = 4*n
		phi, theta, psi = Ori1[ind], Ori1[ind+1], Ori1[ind+2]
		map = False
		if theta > 90:
			theta = theta - 2 * (theta - 90)
			map   = True
		phi   *= deg2rad
		theta *= deg2rad
		x      = sin(theta) * sin(phi)
		y      = sin(theta) * cos(phi)
		val = 1 - x*x - y*y
		if val < 0: val = 0
		z = sqrt(val)
		if map: z = -z
		
		U1.append([x, y, z])

	U2 = []
	for n in xrange(nori):
		ind = 4*n
		phi, theta, psi = Ori2[ind], Ori2[ind+1], Ori2[ind+2]
		map = False
		if theta > 90:
			theta = theta - 2 * (theta - 90)
			map   = True
		phi   *= deg2rad
		theta *= deg2rad
		x      = sin(theta) * sin(phi)
		y      = sin(theta) * cos(phi)
		val = 1 - x*x - y*y
		if val < 0: val = 0
		z = sqrt(val)
		if map: z = -z
		
		U2.append([x, y, z])

	# compute all Suv with uv = {xx, xy, xz, yx, ..., zz}
	Suv   = [0] * 9
	c     = 0
	for i in xrange(3):
		for j in xrange(3):
			for s in xrange(nori):
				Suv[c] += (U2[s][i] * U1[s][j])
			c += 1

        # create matrix N
	N = array([[Suv[0]+Suv[4]+Suv[8], Suv[5]-Suv[7],         Suv[6]-Suv[2],                 Suv[1]-Suv[3]], 
		   [Suv[5]-Suv[7],        Suv[0]-Suv[4]-Suv[8],  Suv[1]+Suv[3],                 Suv[6]+Suv[2]], 
		   [Suv[6]-Suv[2],        Suv[1]+Suv[3],        -Suv[0]+Suv[4]-Suv[8],          Suv[5]+Suv[7]],
		   [Suv[1]-Suv[3],        Suv[6]+Suv[2],         Suv[5]+Suv[7],         -Suv[0]-Suv[4]+Suv[8]]])

        # eigenvector corresponding to the most positive eigenvalue
	val, vec = linalg.eig(N)
	q0, qx, qy, qz = vec[:, val.argmax()]

        # create quaternion Rot matrix
	Rot = matrix([[q0*q0+qx*qx-qy*qy-qz*qz, 2*(qx*qy-q0*qz),                 2*(qx*qz+q0*qy)],
		      [2*(qy*qx+q0*qz),         q0*q0-qx*qx+qy*qy-qz*qz,         2*(qy*qz-q0*qx)],
		      [2*(qz*qx-q0*qy),         2*(qz*qy+q0*qx),         q0*q0-qx*qx-qy*qy+qz*qz]])

        # calculate error between the two sets of marker (err = sum p1t * Rot * p2)
	sum_err = 0
	for n in xrange(nori):
		p1  = matrix([U1[n]]) # row
		p2  = matrix([U2[n]])
		p2t = p2.transpose()  # column
		err = p1 * Rot * p2t
		sum_err += err.getA()[0][0]
	
	# reverse the error, like that if is the same orientation err = 0 and not err = nori
	sum_err = nori - sum_err

	return sum_err

# this function return the degree of colinearity of the orientations found (if colinear the value is close to zero)
def cml2_ori_collinearity(Ori):
	from math  import sin, cos, pi
	from numpy import array, linalg, matrix, zeros, power

	# ori 3d sphere map to 2d plan
	rad2deg = 180 / pi
	deg2rad = 1 / rad2deg
	nori    = len(Ori) // 4
	lx, ly  = [], []
	for n in xrange(nori):
		ind     = n * 4
		ph, th  = Ori[ind:ind+2]
		ph     *= deg2rad
		th     *= deg2rad
		lx.append(sin(th) * sin(ph))
		ly.append(sin(th) * cos(ph))

	# direct least square fitting of ellipse (IEEE ICIP Pilu96)
	D = zeros((nori, 6))
	for c in xrange(nori):
		D[c, :] = [lx[c]*lx[c], lx[c]*ly[c], ly[c]*ly[c], lx[c], ly[c], 1.0]
	D = matrix(D)
	S = D.transpose() * D
	C = zeros((6, 6))
	C[5, 5] = 0
	C[0, 2] = -2
	C[1, 2] = 1
	C[2, 0] = -2
	C = matrix(C)
	val, vec = linalg.eig(S.getI() * C)
	ell = vec[:, val.argmin()]
	verr = D * ell
	verr = power(verr, 2)
	serr = sum(verr)

	# sum squared error
	return serr.getA()[0][0]

# this function is used in find_structure MPI version, to merge and mutatted solution, GA string digit
def cml2_GA_digit(Ori1, Ori2, npop, pcross, pmut):
	from random import randint, random
	import sys

	nori = len(Ori1) // 4
	ngen = nori * 3

	# encode each gen to 3 digit 000
	parent1 = ''
	parent2 = ''
	for n in xrange(nori):
		ind1 = n * 4
		ind2 = ind1 + 1
		ind3 = ind1 + 2
		parent1 += '%03i%03i%03i' % (Ori1[ind1], Ori1[ind2], Ori1[ind3])
		parent2 += '%03i%03i%03i' % (Ori2[ind1], Ori2[ind2], Ori2[ind3])

	#print 'p1 ', parent1
	#print 'p2 ', parent2

	# create new pop
	POP = []
	for ipop in xrange(npop // 2):
		# uniform cross-over
		if random() < pcross:
			child1 = ''
			child2 = ''
			#pat    = ''
			for iori in xrange(nori):
				ind1 = 9 * iori
				ind2 = ind1 + 9
				if random() < 0.6:
					child1 += parent1[ind1:ind2]
					child2 += parent2[ind1:ind2]
					#pat += '111111111'
				else:
					child1 += parent2[ind1:ind2]
					child2 += parent1[ind1:ind2]
					#pat += '000000000'
		else:
			child1 = parent1
			child2 = parent2
		
		#print '\npa ', pat
		#print 'c1 ', child1
		#print 'c2 ', child2

		# mutation
		mut1 = ''
		mut2 = ''
		#pat  = ''
		for igen in xrange(ngen):
			ind1 = 3 * igen
			ind2 = ind1 + 3
			if random() < pmut:
				v = randint(0, 359)
				mut1 += '%03i' % v
				#pat += 'xxx'
			else:	
				mut1 += child1[ind1:ind2]
				#pat += '___'
			if random() < pmut:
				v = randint(0, 359)
				mut2 += '%03i' % v
			else:	mut2 += child2[ind1:ind2]
		
		#print '\npa ', pat
		#print 'm1 ', mut1
		#print 'm2 ', mut2
			
		# decode
		pop1 = []
		pop2 = []

		for iori in xrange(nori):
			ind1 = 9 * iori
			ind2 = ind1 + 3
			for i in xrange(3):
				pop1.append(float(mut1[ind1:ind2]))
				pop2.append(float(mut2[ind1:ind2]))
				ind1 += 3
				ind2 += 3
			pop1.append(-1)
			pop2.append(-1)

		#print 'p1 ', pop1
		#print 'p2 ', pop2

		POP.append(pop1)
		POP.append(pop2)
	
		#sys.exit()
	return POP

# this function is used in find_structure MPI version, to merge and mutatted solution, GA string binary
def cml2_GA(Ori1, Ori2, npop, pcross, pmut):
	from random import randint, random
	import sys

	# reduce to 9 bits (value max of phi and theta is 360)
	def int2bin(n):
		n = int(n)
		res = ''
		while n > 0:
			res = str(n % 2) + res
			n   = n >> 1
		return res.rjust(9, '0')

	nori = len(Ori1) // 4
	ngen = nori * 2 * 9

	# first map Ori to binary string
	chromo1 = ''
	chromo2 = ''
	for n in xrange(nori):
		ind      = n*4
		chromo1 += (int2bin(Ori1[ind]) + int2bin(Ori1[ind+1]))
		chromo2 += (int2bin(Ori2[ind]) + int2bin(Ori2[ind+1]))

	POP = []
	for ipop in xrange(npop):
		# Cycle crossover
		#co = ''
		#co = co.join([str(randint(0, 1)) for i in xrange(ngen)])
		#child = ''
		#for i in xrange(ngen):
		#	if co[i] == '0': child += chromo1[i]
		#	else:            child += chromo2[i]

		if random() <= pcross:
			# One-point crossover
			pt = randint(0, ngen)
			if   pt == 0:    child = chromo1
			elif pt == ngen: child = chromo2
			else:		 child = chromo1[:pt] + chromo2[pt:] 
		else:
			child = chromo1

		# mutation
		for n in xrange(ngen):
			if random() <= pmut:
				pt = randint(0, ngen - 1)
				if child[pt] == '0': child = child[:pt] + '1' + child[pt+1:]
				else:                child = child[:pt] + '0' + child[pt+1:]

		Ori = []
		for n in xrange(nori):
			ind = n*18
			Ori.extend([int(child[ind:ind+9], 2), int(child[ind+9:ind+18], 2), 0, -1]) 

		POP.append(Ori)

	#seq  = range(int(vmin), int(vmax), 10)
	#p1   = choice(seq)
	#seq.remove(p1)
	#p2   = choice(seq)

	return POP



# plot for each individu in the population all orientation on the top hemisphere
def cml2_plot_POP(POP, out_dir, igen):
	from projection import plot_angles
	npop  = len(POP)
	nagls = len(POP[0]) // 4
	for n in xrange(npop):
		agls = []
		for i in xrange(nagls):	agls.append([POP[n][4*i], POP[n][4*i+1], POP[n][4*i+2]])
		im = plot_angles(agls)
		im.write_image(out_dir + '/plot_angles_pop_%03i.hdf' % igen, n)

# application find structure MPI version
def cml2_main_mpi(stack, out_dir, ir, ou, delta, dpsi, lf, hf, rand_seed, maxit, given = False, first_zero = False, flag_weights = False, debug = False, maxgen = 10, pcross = 0.95, pmut = 0.05):
	from mpi        import mpi_init, mpi_comm_size, mpi_comm_rank, mpi_barrier, MPI_COMM_WORLD
	from mpi        import mpi_reduce, mpi_bcast, MPI_INT, MPI_LOR, MPI_FLOAT, MPI_SUM, mpi_send, mpi_recv
	from utilities  import print_begin_msg, print_msg, print_end_msg, start_time, running_time
	from random     import seed, random
	from copy       import deepcopy
	from projection import cml_open_proj, cml_init_global_var, cml_disc, cml_export_txtagls, cml_find_structure
	import time, sys, os, cPickle, logging
	
	# init
	sys.argv = mpi_init(len(sys.argv), sys.argv)
	ncpu     = mpi_comm_size(MPI_COMM_WORLD)
	myid     = mpi_comm_rank(MPI_COMM_WORLD)

	# chose a main one
	main_node = 0

	if myid == main_node:
		# logfile
		#print_begin_msg('find_struct')

		out_dir = out_dir.rstrip('/')
		if os.path.exists(out_dir): os.system('rm -rf ' + out_dir)
		os.mkdir(out_dir)

	mpi_barrier(MPI_COMM_WORLD)
	logging.basicConfig(filename = out_dir + '/main_log.txt', format = '%(asctime)s     %(message)s', level = logging.INFO)

	# manage random seed
	rnd = []
	for n in xrange(1, ncpu + 1): rnd.append(n * (10**n))
	if rand_seed > 0: seed(rnd[myid])
	else:             seed()
	
	# Open and transform projections
	Prj, Ori = cml_open_proj(stack, ir, ou, lf, hf)
	
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

	# Init the global vars
	cml_init_global_var(dpsi, delta, len(Prj), debug)
	
	# Update logfile
	#cml_head_log(stack, outdir, delta, ir, ou, lf, hf, rand_seed, maxit, given)

	POP       = [[] for i in xrange(ncpu)]
	POP[myid] = deepcopy(Ori)

	#pmutmin = 0.01
	#pmutmax = 0.3
	#F       = 0.5
	g_n_prj = len(Ori) // 4
	flag_stop = 0
	for igen in xrange(maxgen):
		t_start = start_time()

		if myid == main_node: logging.info('----- GEN %03i -----' % igen)

		Ori = deepcopy(POP[myid])
		Rot = Util.cml_init_rot(Ori)
		disc = cml_disc(Prj, Ori, Rot, flag_weights)
		cml_export_txtagls(out_dir, 'angles_node_%02i' % myid, Ori, disc, 'Init')

		# gather first disc
		first_disc = [0.0] * ncpu
		first_disc[myid] = disc
		mpi_barrier(MPI_COMM_WORLD)
		first_disc = mpi_reduce(first_disc, ncpu, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
		first_disc = first_disc.tolist()
		if myid == main_node:
			meandisc = sum(first_disc) / float(ncpu)
			logging.info('>>> FIRST DISC: min %7.2f max %7.2f mean %7.2f' % (min(first_disc), max(first_disc), meandisc))
			first_disc = tuple(first_disc)
			txt = '%6.1f ' * len(first_disc)
			logging.info(txt % first_disc)

		#logging.info('[gen: %03i node: %02i]  first disc: %f' % (igen, myid, disc))

		# Find structure
		Ori, disc, ite = cml_find_structure(Prj, Ori, Rot, out_dir, 'angles_node_%02i' % myid, maxit, first_zero, flag_weights)
		#logging.info('[gen: %03i node: %02i]  disc: %f  nb ite: %i' % (igen, myid, disc, ite))
		#f = open(out_dir + '/gen_%03i_node_%02i' % (igen, myid), 'w')
		#cPickle.dump(Ori, f)
		#f.close()

		# gather disc
		all_disc = [0.0] * ncpu
		all_disc[myid] = disc
		mpi_barrier(MPI_COMM_WORLD)
		all_disc = mpi_reduce(all_disc, ncpu, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
		all_disc = all_disc.tolist()

		# gather Ori
		POP[myid] = deepcopy(Ori)
		for icpu in xrange(ncpu):
			if icpu != main_node:
				if myid == icpu:
					mpi_send(Ori, len(Ori), MPI_FLOAT, main_node, 0, MPI_COMM_WORLD)
				if myid == main_node:
					loc_ORI   = mpi_recv(len(Ori), MPI_FLOAT, icpu, 0, MPI_COMM_WORLD)
					loc_ORI   = loc_ORI.tolist()
					POP[icpu] = deepcopy(loc_ORI)
			mpi_barrier(MPI_COMM_WORLD)

		if myid == main_node:
			## TO TEST
			'''
			from projection import plot_angles
			for iori in xrange(ncpu):
				cml = Util.cml_line_in3d(POP[iori], g_seq, g_n_prj, g_n_lines)
				tmp_agls = []
				for i in xrange(len(cml) // 2):
					ind = i * 2
					tmp_agls.append([cml[ind], cml[ind+1], 1])
				im = plot_angles(tmp_agls)
				im.write_image(out_dir + '/cml_pop_%03i.hdf' % igen, iori)
			'''
			# degree of collinearity
			colli = []
			for i in xrange(ncpu): colli.append(cml2_ori_collinearity(POP[i]))
			colli = tuple(colli)
			txt = '%5.3f ' * ncpu
			logging.info(txt % colli)

			meandisc = sum(all_disc) / float(ncpu)
			minerr = 1e6
			mem    = [-1, -1]
			serr   = 0
			for i in xrange(ncpu):
				for j in xrange(i+1, ncpu):
					err   = cml2_error_ori(POP[i], POP[j])
					err  /= (colli[i] * colli[j]) 
					serr += err
					if err < minerr:
						minerr = err
						mem    = [i, j]

			serr /= float((ncpu-1)*ncpu/2.0)
			logging.info('>>> END DISC:   min %7.2f max %7.2f mean %7.2f' % (min(all_disc), max(all_disc), meandisc))
			all_disc = tuple(all_disc)
			txt = '%6.1f ' * ncpu
			logging.info(txt % all_disc)
			
			logging.info('>>> best pair is %s (disc %6.1f %6.1f) with err = %f   err mean: %f' % (mem, all_disc[mem[0]], all_disc[mem[1]], minerr, serr))
			if all_disc[mem[0]] < all_disc[mem[1]]:
				best_mem = mem[0]
			else:
				best_mem = mem[1]
			cml2_export_struc_GA(stack, out_dir, POP[best_mem], igen)
			logging.info('>>> export best structure given by node %i (disc %6.1f)' % (best_mem, all_disc[best_mem]))

			cml2_plot_POP(POP, out_dir, igen)

			POP = cml2_GA_digit(POP[mem[0]], POP[mem[1]], ncpu, pcross, pmut)
			#logging.info('>>> new populations pcross: %f pmut: %f' % (pcross, pmut))
		
			#if minerr == 0: flag_stop = 1

		# Broadcast new POP
		mpi_barrier(MPI_COMM_WORLD)
		if myid != main_node: POP = [[] for i in xrange(ncpu)]
		for icpu in xrange(ncpu):
			Ori       = deepcopy(POP[icpu])
			Ori       = mpi_bcast(Ori, g_n_prj * 4, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			Ori       = Ori.tolist()
			POP[icpu] = deepcopy(Ori)
			
		mpi_barrier(MPI_COMM_WORLD)
		flag_stop = mpi_reduce(flag_stop, 1, MPI_INT, MPI_LOR, main_node, MPI_COMM_WORLD)
		flag_stop = mpi_bcast(flag_stop, 1, MPI_INT, main_node, MPI_COMM_WORLD)
		flag_stop = flag_stop.tolist()[0]
		if flag_stop: break

	
	#if myid == main_node:
	        # Export structure
		#cml2_export_struc(stack, out_dir, Ori)
	if myid == main_node:
		logging.info('=== END ===')

	# Update logfile
	#cml_end_log(Ori, disc, disc_now, ite)

	#running_time(t_start)
	#print_end_msg('find_struct')

	#return disc_now








## END COMMON LINES NEW VERSION ###############################################################
###############################################################################################

"""

def ali3d_dB(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an="-1",
	    center = 1, maxit = 5, CTF = False, snr = 1.0,  ref_a="S", sym="c1"):
	#  THIS IS FOR PROCESSING BERLIN DATASET
	from utilities      import model_circle, reduce_EMData_to_root, bcast_EMData_to_all, drop_image
	from utilities      import bcast_string_to_all, get_image, get_input_from_string
	from utilities      import get_arb_params, set_arb_params, drop_spider_doc,recv_attr_dict, send_attr_dict
	from utilities      import drop_spider_doc,get_im
	from filter	    import filt_params, filt_btwl, filt_from_fsc, filt_table, filt_gaussl
	from alignment	    import proj_ali_incore
	from random	    import randint
	from fundamentals   import fshift,rot_avg_image
	import os
	import types
	from string         import replace
	from reconstruction import rec3D_MPI, rec3D_MPI_noCTF
	from mpi 	    import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	    import mpi_barrier
	
	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	if(myid == main_node):
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)
        finfo = None
	#info_file = replace("progress%4d"%myid, ' ', '0')
        #finfo = open(info_file, 'w')

	max_iter    = int(maxit)
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
		from  alignment	    import proj_ali_incore_localB
		an      = get_input_from_string(an)
	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	if last_ring < 0:	last_ring = int(nx/2) - 2
	if(maskfile):
		if(type(maskfile) is types.StringType): mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else:         mask3D = model_circle(last_ring, nx, nx, nx)
	nima            = EMUtil.get_image_count(stack)
	mask = model_circle(last_ring, nx, nx)
	
	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)

	if(myid == main_node)       :  drop_image(vol, os.path.join(outdir,"ref_volf00.spi"), "s")
	if(CTF):
		ctf_dicts = ["defocus", "Cs", "voltage", "Pixel_size", "amp_contrast", "B_factor", "ctf_applied" ]
		#                 0      1         2            3               4                  5               6
		from reconstruction import rec3D_MPI
	else:
		from reconstruction import rec3D_MPI_noCTF
	from utilities import read_spider_doc, set_arb_params
	prm = read_spider_doc("params_new.doc")
	prm_dict = ["phi", "theta", "psi", "s2x", "s2y"]
	for im in xrange(image_start, image_end):
		data[im-image_start].set_attr('ID', im)
		set_arb_params(ima, prm[im], prm_dict)
		if(CTF):
			ctf_params = get_arb_params(data[im], ctf_dicts)
			if(im == image_start): data_had_ctf = ctf_params[6]
			if(ctf_params[6] == 0):
				st = Util.infomask(data[im-image_start], mask, False)
				data[im-image_start] -= st[0]
				from filter import filt_ctf
				data[im-image_start] = filt_ctf(data[im], ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
				data[im-image_start].set_attr('ctf_applied', 1)
	del prm
	# do the projection matching
	from  string        import replace
	from utilities      import bcast_number_to_all
	from morphology     import adaptive_mask, threshold, threshold_to_minval
	from statistics     import histogram
	from reconstruction import recons3d_nn_SSNR_MPI
	for N_step in xrange(lstp):
 		for Iter in xrange(max_iter):
			#print " ali3d_d_MPI: ITERATION #",N_step*max_iter + Iter+1
			#  We have to find a mask on a filtered volume, apply it to non filtered, filter, 
			#          derive adjustment of rotational average on thresholded filtered volume
			#          finally apply mask and adjustment to not-filtered and fitler it
			#  vol   - filtered
			#  volo  - not-filtered
			#if(N_step == 0 and Iter == 0):
			#	adam = adaptive_mask(vol, 3000, 2.22)
			#	vol = threshold( adam*vol)
			#	h = histogram( vol )
			#	vol *= get_im("rotav_tteftu_with_tRNA.spi") / threshold_to_minval( rot_avg_image(vol), h[len(h)//2])
			#	vol = filt_btwl(vol, 0.3, 0.4)
			if myid == 0:
				[ssnr1, vol_ssnr] = recons3d_nn_SSNR_MPI(myid, data, model_circle(last_ring, nx, nx), ctf = CTF)
				del ssnr1
				#  change volume to fsc
				vol_ssnr /= vol_ssnr + 1.0
				#  filter it somewhat
				vol_ssnr = threshold(filt_gaussl(vol_ssnr, 0.1))
			else:  recons3d_nn_SSNR_MPI(myid, data, model_circle(last_ring, nx, nx), ctf = CTF)

			if(CTF): volo, fscc = rec3D_MPI(data, snr, sym, mask3D, os.path.join(outdir, replace("resolution%4d"%(N_step*max_iter+Iter+1),' ','0')), myid, main_node)
			else:    volo, fscc = rec3D_MPI_noCTF(data, sym, mask3D, os.path.join(outdir, replace("resolution%4d"%(N_step*max_iter+Iter+1),' ','0')), myid, main_node)

			mpi_barrier(MPI_COMM_WORLD)
			if(myid == main_node):
				drop_image(volo,  os.path.join(outdir, replace("vol%4d.spi"%(N_step*max_iter+Iter+1),' ','0')), "s")
				from fundamentals import rot_avg
				filt = filt_from_fsc(fscc, 0.1)
				# This is a cheat - it moves the resolution curve to the right making the filter more lenient
				#filt = [ 1.0, 1.0, 1.0] + filt
				#filt = [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0] + filt
				rfsc = rot_avg(vol_ssnr)
				rfsc.set_value_at(0, 1.0)
				nrfsc = rfsc.get_xsize()
				for irf in xrange(nrfsc):
					qtq = rfsc.get_value_at(irf)
					if(qtq > 0.0):   rfsc.set_value_at(irf, filt[irf]/qtq)
					else:            rfsc.set_value_at(irf, 0.0)
				vol_ssnr = vol_ssnr.mult_radial(rfsc)  # adjust radially fsc filter to fsc curve
				drop_image(vol_ssnr,  os.path.join(outdir, replace("filter%4d.spi"%(N_step*max_iter+Iter+1),' ','0')), "s")
				#  Filter by volume
				vol = volo.filter_by_image(vol_ssnr)
				#vol  = filt_table(volo, filt)
				if(center == 1):
					cs   = vol.phase_cog()
					vol  = fshift(vol, -cs[0], -cs[1] -cs[2])
					volo = fshift(volo, -cs[0], -cs[1] -cs[2])
				drop_image(vol, os.path.join(outdir, replace("volf%4d.spi"%(N_step*max_iter+Iter+1),' ','0')), "s")
				adam = adaptive_mask(vol, 3000, 2.22)
				vol = threshold( adam*vol )
				h = histogram( vol )
				vol = threshold( adam*volo ) * get_im("rotav_tteftu_with_tRNA.spi") / threshold_to_minval( rot_avg_image(vol), h[len(h)//2])
				del volo
				#  Filter by volume
				vol = vol.filter_by_image(vol_ssnr)
				#vol  = filt_table(vol, filt)
				#vol = filt_btwl(vol, fl, fh)
				drop_image(vol, os.path.join(outdir, replace("vhlf%4d.spi"%(N_step*max_iter+Iter+1),' ','0')), "s")
				del adam, vol_ssnr, rfsc, filt
				del h
			bcast_EMData_to_all(vol, myid, main_node)
			if(an[N_step] == -1): proj_ali_incore(vol, mask3D, data, first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], ref_a, sym, finfo, MPI=True)
			else:                 proj_ali_incore_localB(vol, mask3D, data, first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], an[N_step], ref_a, sym, finfo, MPI=True)
			prm = []
			for im in xrange(image_start, image_end):
				prml = get_arb_params(data[im-image_start], prm_dict)
				prml.append(im)
				prm.append(prml)
			drop_spider_doc(os.path.join(outdir, replace("new_params%8d"%((N_step*max_iter+Iter+1)*1000+myid),' ','0')), prm," phi, theta, psi, s2x, s2y, image number")
			del prm, prml

	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
	if(CTF and data_had_ctf == 0):
		for im in xrange(len(data)): data[im].set_attr('ctf_applied', 0)
	#par_str = ['phi', 'theta', 'psi', 's2x', 's2y']
	if(myid == main_node): recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	else: send_attr_dict(main_node, data, par_str, image_start, image_end)






def ssnr3d(stack, output_volume = None, ssnr_text_file = None, mask = None, ou = -1, rw = 1.0,  npad = 1, CTF = False, sign = 1, sym ="c1", random_angles = 0):

'''
	Perform 3D reconstruction using selected particle images, 
	and calculate spectrum signal noise ratio (SSNR).
	1. The selection file is supposed to be in SPIDER text format.
	2. The 3D alignment parameters have been written in headers of the particle images.
''' 


	from global_def import MPI
	if MPI:
		ssnr3d_MPI(stack, output_volume, ssnr_text_file, mask, ou, rw, npad, CTF, sign, sym, random_angles)
		return

	from utilities import read_spider_doc, drop_image, get_arb_params, set_arb_params, model_circle, get_im
	from filter import filt_ctf
	from reconstruction import recons3d_nn_SSNR, recons3d_4nn, recons3d_4nn_ctf
	from projection import prep_vol, prgs
	from utilities import print_begin_msg, print_end_msg, print_msg
	
	print_begin_msg("ssnr3d")

	if output_volume  is None: output_volume  = "SSNR"
	if ssnr_text_file is None: ssnr_text_file = "ssnr"

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Output volume               : %s\n"%(output_volume))	
	print_msg("SSNR text file              : %s\n"%(ssnr_text_file))
	print_msg("Outer radius                : %i\n"%(ou))
	print_msg("Ring width                  : %i\n"%(rw))
	print_msg("Padding factor              : %i\n"%(npad))
	print_msg("data with CTF               : %s\n"%(CTF))
	print_msg("CTF sign                    : %i\n"%(sign))
	print_msg("Symmetry group              : %s\n\n"%(sym))
	
	nima = EMUtil.get_image_count(stack)
	fring_width = float(rw)
	if mask:
		import  types
		if(type(mask) is types.StringType):
			print_msg("Maskfile                    : %s\n"%(mask))
			mask2D=get_im(mask)
		else:
			print_msg("Maskfile                    : user provided in-core mask")
			mask2D = mask
	else:
		print_msg("Maskfile                    : None")
		mask2D = None

	[ssnr1, vol_ssnr1] = recons3d_nn_SSNR(stack, mask2D, rw, npad, sign, sym, CTF, random_angles)
	drop_image(vol_ssnr1, output_volume+"1.spi", "s")
	outf = file(ssnr_text_file+"1.txt", "w")
	for i in xrange(len(ssnr1)):
		datstrings = []
		datstrings.append("  %15f" % ssnr1[0][i])    #  have to subtract 0.5 as in C code there is round.
		datstrings.append("  %15e" % ssnr1[1][i])   # SSNR
		datstrings.append("  %15e" % ssnr1[2][i])   # variance divided by two numbers
		datstrings.append("  %15f" % ssnr1[3][i])		 # number of points in the shell
		datstrings.append("  %15f" % ssnr1[4][i])		 # number of added Fourier points
		datstrings.append("  %15f" % ssnr1[5][i])		 # square of signal
		datstrings.append("\n")
		outf.write("".join(datstrings))
	outf.close()
	print_end_msg("ssnr3d")
	return ssnr1, vol_ssnr1
	'''
	# perform 3D reconstruction
	if CTF:
		snr = 1.0e20
		vol = recons3d_4nn_ctf(stack, range(nima), snr, sign, sym)
	else :   vol = recons3d_4nn(stack, range(nima), sym)
	# re-project the reconstructed volume
	if CTF :img_dicts = ["phi", "theta", "psi", "s2x", "s2y", "defocus", "Pixel_size",\
                  "voltage", "Cs", "amp_contrast", "sign", "B_factor", "active", "ctf_applied"]
	else   :img_dicts = ["phi", "theta", "psi", "s2x", "s2y", "active"]
	nx = vol.get_xsize()
	if int(ou) == -1: radius = nx//2 - 1
	else :            radius = int(ou)
	#
	prjlist = []
	vol *= model_circle(radius, nx, nx, nx)
	volft,kb = prep_vol(vol)
	del vol
	for i in xrange(nima):
		e = EMData()
		e.read_image(stack, i, True)
		e.set_attr('sign', 1)
		params = get_arb_params(e, img_dicts)
		#proj = project(vol,[params[0], params[1], params[2], params[3], params[4]] , radius)
		proj = prgs(volft, kb, [params[0], params[1], params[2], params[3], params[4]])
		if CTF :  proj = filt_ctf(proj, params[5], params[8], params[7], params[6], params[9], params[11])
		set_arb_params(proj, params, img_dicts)
		if(CTF):	 proj.set_attr('ctf_applied', 1)
		else:		 proj.set_attr('ctf_applied', 0)
		prjlist.append(proj)
	del volft
	[ssnr2, vol_ssnr2] = recons3d_nn_SSNR(prjlist, mask2D, CTF, sym, npad, sign, fring_width, filename=ssnr_text_file+"2.txt")
	qt = 0.0
	for i in xrange(len(ssnr1)):
		tqt = ssnr1[i][1] - ssnr2[i][1]
		if( tqt<qt ): qt = tqt
	for i in xrange(len(ssnr1)): ssnr1[i][1] -= (ssnr2[i][1] + qt)
	from utilities import drop_spider_doc19289
	drop_spider_doc(ssnr_text_file+".doc", ssnr1)
	drop_image(vol_ssnr2, output_volume+"2.spi", "s")
	

def ssnr3d_MPI(stack, output_volume = None, ssnr_text_file = None, mask = None, ou = -1, rw = 1.0, npad = 1, CTF = False, sign = 1, sym ="c1", random_angles = 0):
	from reconstruction import recons3d_nn_SSNR_MPI, recons3d_4nn_MPI, recons3d_4nn_ctf_MPI
	from string import replace
	from time import time
	from utilities import info, get_arb_params, set_arb_params, bcast_EMData_to_all, model_blank, model_circle, get_im, drop_image
	from filter import filt_ctf
	from projection import prep_vol, prgs

	if output_volume  is None: output_volume  = "SSNR.spi"
	if ssnr_text_file is None: ssnr_text_file = "ssnr"
	from utilities import read_spider_doc

	#   TWO NEXT LINEs FOR EXTRA PROJECT
	#params_ref = read_spider_doc("params_new.txt")
	#headers = ["phi", "theta", "psi", "s2x", "s2y"]

	nima = EMUtil.get_image_count(stack)
	nproc = mpi_comm_size(MPI_COMM_WORLD)
	myid  = mpi_comm_rank(MPI_COMM_WORLD)

	#nimage_per_node = nima/nproc
	#image_start     = myid * nimage_per_node
	#if (myid == nproc-1): image_end = nima
	#else:	              image_end = image_start + nimage_per_node
	image_start, image_end = MPI_start_end(nima, nproc, myid)

	if mask:
		import  types
		if(type(mask) is types.StringType):  mask2D=get_im(mask)
		else: mask2D = mask
	else:
		mask2D = None

	prjlist = []
	for i in range(image_start, image_end):
		prj = EMData()
		prj.read_image( stack, i)
		#  TWO NEXT LINEs FOR EXTRA PROJECT
		#tmp_par = [params_ref[i][0], params_ref[i][1], params_ref[i][2], params_ref[i][3], params_ref[i][4], 0]
		#set_arb_params(prj, params_ref[i], headers)
		prjlist.append( prj )
	if myid == 0: [ssnr1, vol_ssnr1] = recons3d_nn_SSNR_MPI(myid, prjlist, mask2D, rw, npad, sign, sym, CTF, random_angles)  
	else:	                           recons3d_nn_SSNR_MPI(myid, prjlist, mask2D, rw, npad, sign, sym, CTF, random_angles)
	if myid == 0:
		drop_image(vol_ssnr1, output_volume+"1.spi", "s")
		outf = file(ssnr_text_file+"1.txt", "w")
		for i in xrange(len(ssnr1)):
			datstrings = []
			datstrings.append("  %15f" % ssnr1[0][i])    #  have to subtract 0.5 as in C code there is round.
			datstrings.append("  %15e" % ssnr1[1][i])    # SSNR
			datstrings.append("  %15e" % ssnr1[2][i])    # variance divided by two numbers
			datstrings.append("  %15f" % ssnr1[3][i])    # number of points in the shell
			datstrings.append("  %15f" % ssnr1[4][i])    # number of added Fourier points
			datstrings.append("  %15f" % ssnr1[5][i])    # square of signal
			datstrings.append("\n")
			outf.write("".join(datstrings))
		outf.close()
		return ssnr1, vol_ssnr1

	'''
	nx  = prjlist[0].get_xsize()
	vol = model_blank(nx,nx,nx)
	if ou == -1: radius = int(nx/2) - 1
	else:        radius = int(ou)
	if CTF :
		snr = 1.0e20
		if myid == 0 : vol = recons3d_4nn_ctf_MPI(myid, prjlist, snr, sign, sym)
		else :  	     recons3d_4nn_ctf_MPI(myid, prjlist, snr, sign, sym)
	else  :
		if myid == 0 : vol = recons3d_4nn_MPI(myid, prjlist, sym)
		else:		     recons3d_4nn_MPI(myid, prjlist, sym)
	bcast_EMData_to_all(vol, myid, 0)
	if CTF: img_dicts = ["phi", "theta", "psi", "s2x", "s2y", "defocus", "Pixel_size",\
	                    "voltage", "Cs", "amp_contrast", "sign", "B_factor", "active", "ctf_applied"]
	else   : img_dicts = ["phi", "theta", "psi", "s2x", "s2y"]
	re_prjlist = []
	vol *= model_circle(radius, nx, nx, nx)
	volft,kb = prep_vol(vol)
	del vol
	for i in xrange(image_start, image_end):
		prjlist[i-image_start].set_attr('sign', 1)
		params = get_arb_params(prjlist[i-image_start], img_dicts)
		#proj   = project(vol, [params[0], params[1], params[2], params[3], params[4]] , radius)
		proj = prgs(volft, kb, [params[0], params[1], params[2], params[3], params[4]])
		if CTF: proj = filt_ctf(proj, params[5], params[8], params[7], params[6], params[9], params[11])
		set_arb_params(proj,params,img_dicts)
		if(CTF):	 proj.set_attr('ctf_applied', 1)
		else:		 proj.set_attr('ctf_applied', 0)
		re_prjlist.append(proj)
	del volft
	if myid == 0: [ssnr2, vol_ssnr2] = recons3d_nn_SSNR_MPI(myid, re_prjlist, ssnr_text_file+"2.txt", mask2D, rw, npad, sign, sym, CTF)
	else:                              recons3d_nn_SSNR_MPI(myid, re_prjlist, ssnr_text_file+"2.txt", mask2D, rw, npad, sign, sym, CTF)
	if myid == 0 :
		qt = 0.0
		for i in xrange(len(ssnr1)):
			tqt = ssnr1[i][1] - ssnr2[i][1]
			if( tqt<qt ): qt = tqt
		for i in xrange(len(ssnr1)): ssnr1[i][1] -= (ssnr2[i][1] + qt)
		from utilities import drop_spider_doc
		drop_spider_doc(ssnr_text_file+".doc", ssnr1)
		drop_image(vol_ssnr2, output_volume+"2.spi", "s")
"""


def ali3d_eB_MPI_LAST_USED(stack, ref_vol, outdir, maskfile, ou=-1,  delta=2, maxit=10, CTF = None, snr=1.0, sym="c1", chunk = -1.0, user_func_name="ref_aliB_cone"):
	'''
		Cone
	'''
	from alignment	    import proj_ali_incore_cone
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl, filt_tanl
	from fundamentals   import fshift, rot_avg_image
	from projection     import prep_vol, prgs
	from utilities      import bcast_string_to_all, model_circle, get_arb_params, set_arb_params, drop_spider_doc
	from utilities      import get_image, drop_image, bcast_EMData_to_all, send_attr_dict, recv_attr_dict
	from utilities      import read_spider_doc, get_im
	from reconstruction import rec3D_MPI
	from statistics     import ccc
	from math           import pi, sqrt
	from string         import replace
	import os
	import sys
	from mpi 	    import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	    import mpi_barrier

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	if(myid == main_node):
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)
		import user_functions
		user_func = user_functions.factory[user_func_name]
	mpi_barrier(MPI_COMM_WORLD)

	nima = EMUtil.get_image_count(stack)

	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)

	# figure the size of the chunk (3D is updated after each chunk).  Chunk should be given as 0.0< chunk <= 1.0.  1.0 means all projections
	if(chunk <= 0.0):  chunk = 1.0
	n_in_chunk  = max(int(chunk * (image_end-image_start+1)), 1)
	n_of_chunks = (image_end-image_start+1)//n_in_chunk + min((image_end-image_start+1)%n_in_chunk,1)
	
	outf = file(replace("progress%4d"%myid,' ','0'), "w")
	outf.write("  chunk = "+str(chunk)+"   ")
	outf.write("\n")
	outf.flush()
	outf.write("  n_in_chunk = "+str(n_in_chunk)+"   ")
	outf.write("  n_of_chunks = "+str(n_of_chunks)+"   ")
	outf.write("\n")
	outf.flush()
	#  Here we assume that reference volume exists
	vol = EMData()
	vol.read_image(ref_vol)
	nx  = vol.get_xsize()
	#ima = EMData()
	#ima.read_image(stack)
	#nx  = ima.get_xsize()
	if(ou <= 0):  ou = nx//2-2
	#if(myid == main_node):  vol.write_image(os.path.join(outdir,"ref_volf00.hdf"))
	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask3D = get_image(maskfile)
		else:                                   mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)
	mask2D = model_circle(ou, nx, nx)



	#from utilities import read_spider_doc, set_arb_params
	#prm = read_spider_doc("params_new.doc")
	#prm_dict = ["phi", "theta", "psi", "s2x", "s2y"]
	#prm_dict = ["phi", "theta", "psi"]
	#from random import seed,gauss
	#seed()
	dataim = EMData.read_images(image_start, image_end)
	for im in xrange(image_start, image_end):
		dataim[im].set_attr('ID', im)
		#angn = get_arb_params(ima, prm_dict)
		#for ian in xrange(3):  angn[ian] += gauss(0.0, 1.0)
		#set_arb_params(ima, angn, prm_dict)
		#set_arb_params(ima, prm[im], prm_dict)
		if(CTF):
			ctf_params = data[im].get_attr('ctf')
			if(im == image_start): data_had_ctf = data[im].get_attr('ctf_applied')
			if data[im].get_attr('ctf_applied') == 0:
				st = Util.infomask(data[im], mask2D, False)
				data[im] -= st[0]
				from filter import filt_ctf
				data[im] = filt_ctf(data[im], ctf_params)
				data[im].set_attr('ctf_applied', 1)

	outf.write("  data read = "+str(image_start)+"   ")
	outf.write("\n")
	outf.flush()


	if (myid == main_node):
		# initialize data for the reference preparation function
		from utilities import read_text_file
		ref_data = []
		ref_data.append( mask3D )
		ref_data.append( read_text_file("pwpdb.txt", 1) )

	from utilities      import bcast_number_to_all
	from morphology     import threshold, threshold_to_minval

	par_str=["phi", "theta", "psi", "s2x", "s2y"]

	#  this is needed for gathering and scattering of cccfs
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(nima, number_of_proc, im)
		recvcount.append( ie - ib )

	from utilities import even_angles
	xrng = 0.0
	yrng = 0.0
	step = 1.0
	template_angles = even_angles(0.2, 0.0, 0.4, phiEqpsi = 'Zero', method = 'P')
	from filter import  filt_tophatb, filt_gaussl
	fifi = True
	for iteration in xrange(maxit):
		outf.write("  iteration = "+str(iteration)+"   ")
		outf.write("\n")
		outf.flush()
		for  ic  in xrange(n_of_chunks):
			# compute updated 3D after each chunk
			if(fifi):
				# resolution
				outf.write("  begin reconstruction = "+str(image_start)+"   ")
				outf.write("\n")
				outf.flush()
				vol, fscc = rec3D_MPI(dataim, snr, sym, mask3D, os.path.join(outdir, replace("resolution%3d_%3d"%(iteration, ic),' ','0') ), myid, main_node)
				outf.write("  done reconstruction = "+str(image_start)+"   ")
				outf.write("\n")
				outf.flush()

				mpi_barrier(MPI_COMM_WORLD)
				if(myid == main_node):
					drop_image(vol, os.path.join(outdir, replace("vol%3d_%3d.hdf"%(iteration, ic),' ','0') ))
					
					ref_data.append( vol )
					ref_data.append( fscc )
					#  call user-supplied function to prepare reference image, i.e., filter it
					vol = user_func( ref_data )
					#  HERE CS SHOULD BE USED TO MODIFY PROJECTIONS' PARAMETERS  !!!
					del ref_data[2]
					del ref_data[2]
					drop_image(vol, os.path.join(outdir, replace("volf%3d_%3d.hdf"%(iteration, ic),' ','0') ))
				#from sys import exit
				#exit()

				bcast_EMData_to_all(vol, myid, main_node)
			fifi = True
			volft,kb  = prep_vol(vol)

			image_start_in_chunk = image_start + ic*n_in_chunk
			image_end_in_chunk   = min(image_start_in_chunk + n_in_chunk, image_end)
			outf.write("ic "+str(ic)+"   image_start "+str(image_start)+"   n_in_chunk "+str(n_in_chunk)+"   image_end "+str(image_end)+"\n")
			outf.write("image_start_in_chunk "+str(image_start_in_chunk)+"  image_end_in_chunk "+str(image_end_in_chunk)+"\n")
			outf.flush()
			for imn in xrange(image_start_in_chunk, image_end_in_chunk):
				#from utilities import start_time, finish_time
				#t3=start_time()
				proj_ali_incore_cone(volft, kb, template_angles, dataim[imn-image_start], 1, ou, 1, xrng, yrng, step, outf)

			soto = []
			for imn in xrange(image_start, image_end):
				from utilities import set_params_proj, get_params_proj
				phi,theta,psi,s2x,s2y = get_params_proj( dataim[imn-image_start] )
				soto.append([phi,theta,psi,s2x,s2y,imn])
			drop_spider_doc(os.path.join(outdir, replace("new_params%3d_%3d_%3d"%(iteration, ic, myid),' ','0')), soto," phi, theta, psi, s2x, s2y, image number")
			del soto

			#from sys import exit
			#exit()
		#  here we should write header info, just in case the program crashes...
	del vol
	del volft
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)

def ali3d_eB_CCC(stack, ref_vol, outdir, maskfile, ou=-1,  delta=2, maxit=10, CTF = None, snr=1.0, sym="c1", chunk = -1.0, user_func_name="ref_aliB_cone"):
	'''
		Cone, modified version to test CCC
		single processor version
	'''
	from utilities      import print_begin_msg, print_end_msg, print_msg

	from alignment	    import proj_ali_incore_cone
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl, filt_tanl
	from fundamentals   import fshift, rot_avg_image
	from projection     import prep_vol, prgs
	from utilities      import model_circle, get_arb_params, set_arb_params, drop_spider_doc
	from utilities      import get_image, drop_image, send_attr_dict, recv_attr_dict
	from utilities      import read_spider_doc, get_im
	from statistics     import ccc
	from statistics     import fsc_mask
	from math           import pi, sqrt
	from string         import replace
	import os
	import sys

	number_of_proc = 1
	myid = 0
	main_node = 0
	
	if os.path.exists(outdir):  os.system('rm -rf '+outdir)
	os.mkdir(outdir)
	import user_functions
	user_func = user_functions.factory[user_func_name]

	if CTF :from reconstruction import recons3d_4nn_ctf
	else   : from reconstruction import recons3d_4nn


	nima = EMUtil.get_image_count(stack)

	image_start = 0
	image_end   = nima

	# figure the size of the chunk (3D is updated after each chunk).  Chunk should be given as 0.0< chunk <= 1.0.  1.0 means all projections
	if(chunk <= 0.0):  chunk = 1.0
	n_in_chunk  = max(int(chunk * (image_end-image_start+1)), 1)
	n_of_chunks = (image_end-image_start+1)//n_in_chunk + min((image_end-image_start+1)%n_in_chunk,1)
	
	
	outf = file(os.path.join(outdir, "progress"), "w")
	outf.write("  chunk = "+str(chunk)+"   ")
	outf.write("\n")
	outf.flush()
	outf.write("  n_in_chunk = "+str(n_in_chunk)+"   ")
	outf.write("  n_of_chunks = "+str(n_of_chunks)+"   ")
	outf.write("\n")
	outf.flush()
	#  Here we assume that reference volume exists
	vol = EMData()
	vol.read_image(ref_vol)
	nx  = vol.get_xsize()
	#ima = EMData()
	#ima.read_image(stack)
	#nx  = ima.get_xsize()
	ou = int(ou)
	if(ou <= 0):  ou = nx//2-2
	#if(myid == main_node):  vol.write_image(os.path.join(outdir,"ref_volf00.hdf"))
	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask3D = get_image(maskfile)
		else:                                   mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)
	mask2D = model_circle(ou, nx, nx)

	dataim = []
	#from utilities import read_spider_doc, set_arb_params
	#prm = read_spider_doc("params_new.doc")
	#prm_dict = ["phi", "theta", "psi", "s2x", "s2y"]
	#prm_dict = ["phi", "theta", "psi"]
	#from random import seed,gauss
	#seed()
	for im in xrange(image_start, image_end):
		ima = EMData()
		ima.read_image(stack, im)
		#angn = get_arb_params(ima, prm_dict)
		#for ian in xrange(3):  angn[ian] += gauss(0.0, 1.0)
		#set_arb_params(ima, angn, prm_dict)
		#set_arb_params(ima, prm[im], prm_dict)
		if(CTF):
			ctf_params = ima.get_attr('ctf')
			if(im == image_start): data_had_ctf = ima.get_attr('ctf_applied')
			if ima.get_attr('ctf_applied') == 0:
				st = Util.infomask(ima, mask2D, False)
				ima -= st[0]
				from filter import filt_ctf
				ima = filt_ctf(ima, ctf_params)
				ima.set_attr('ctf_applied', 1)
		dataim.append(ima)
	outf.write("  data read = "+str(image_start)+"   ")
	outf.write("\n")
	outf.flush()


	# initialize data for the reference preparation function
	from utilities import read_text_file
	ref_data = []
	ref_data.append( mask3D )
	ref_data.append( read_text_file("pwpdb.txt", 1) )
	from utilities import read_text_file
	fscc = [read_text_file("resolution000_000",0), read_text_file("resolution000_000",1)]
 	jtep = 0

	par_str=["phi", "theta", "psi", "s2x", "s2y"]

	from utilities import even_angles
	xrng = 0.0
	yrng = 0.0
	step = 1.0
	template_angles = even_angles(0.2, 0.0, 3.5, phiEqpsi = 'Zero', method = 'P')
	print  len(template_angles)
	for iteration in xrange(maxit):
		msg = "ITERATION #%3d\n"%(iteration+1)
		print_msg(msg)
		for  ic  in xrange(n_of_chunks):
			jtep += 1
			drop_image(vol, os.path.join(outdir, replace("vol%3d_%3d.hdf"%(iteration, ic),' ','0') ))
			ref_data.append( vol )
			ref_data.append( fscc )
			#  call user-supplied function to prepare reference image, i.e., filter it
			vol = user_func( ref_data )
			#  HERE CS SHOULD BE USED TO MODIFY PROJECTIONS' PARAMETERS  !!!
			del ref_data[2]
			del ref_data[2]
			drop_image(vol, os.path.join(outdir, replace("volf%3d_%3d.hdf"%(iteration, ic),' ','0') ))

			volft,kb  = prep_vol(vol)

			image_start_in_chunk = image_start + ic*n_in_chunk
			image_end_in_chunk   = min(image_start_in_chunk + n_in_chunk, image_end)
			outf.write("ic "+str(ic)+"   image_start "+str(image_start)+"   n_in_chunk "+str(n_in_chunk)+"   image_end "+str(image_end)+"\n")
			outf.write("image_start_in_chunk "+str(image_start_in_chunk)+"  image_end_in_chunk "+str(image_end_in_chunk)+"\n")
			outf.flush()
			for imn in xrange(image_start_in_chunk, image_end_in_chunk):
			#for imn in xrange(2,3):
				proj_ali_incore_cone(volft, kb, template_angles, dataim[imn-image_start], 1, ou, 1, xrng, yrng, step, outf)

			soto = []
			for imn in xrange(image_start, image_end):
				from utilities import set_params_proj, get_params_proj
				phi,theta,psi,s2x,s2y = get_params_proj( dataim[imn-image_start] )
				soto.append([phi,theta,psi,s2x,s2y,imn])
			drop_spider_doc(os.path.join(outdir, replace("new_params%3d_%3d"%(iteration, ic),' ','0')), soto," phi, theta, psi, s2x, s2y, image number")
			del soto
			from sys import exit
			exit()

			# compute updated 3D after each chunk
 	    		# resolution
			#print  " start reconstruction",image_start,image_end
			#  3D stuff
			list_p = range(0,nima,2)
 			if(CTF): vol1 = recons3d_4nn_ctf(stack, list_p, snr, 1, sym)
			else:	 vol1 = recons3d_4nn(stack, list_p, sym)

			list_p = range(1,nima,2)
			if(CTF): vol2 = recons3d_4nn_ctf(stack, list_p, snr, 1, sym)
			else:	 vol2 = recons3d_4nn(stack, list_p, sym)

			fscc = fsc_mask(vol1, vol2, mask3D, 1.0, os.path.join(outdir, replace("resolution%4d"%(iteration*n_of_chunks+ic+1),' ','0')))
			del vol1
			del vol2

			# calculate new and improved 3D
			list_p = range(nima)
			if(CTF): vol = recons3d_4nn_ctf(stack, list_p, snr, 1, sym)
			else:	 vol = recons3d_4nn(stack, list_p, sym)
			# store the reference volume
			#drop_image(vol,os.path.join(outdir, replace("vol%4d.spi"%(N_step*max_iter+Iter+1),' ','0')), "s")
			drop_image(vol,os.path.join(outdir, replace("vol%4d.hdf"%(iteration*n_of_chunks+ic+1),' ','0')), "s")

def ali3d_eB_MPI_conewithselect(stack, ref_vol, outdir, maskfile, ou=-1,  delta=2, maxit=10, CTF = None, snr=1.0, sym="c1", chunk = -1.0):
	'''
		Cone
	'''
	from alignment	    import proj_ali_incore_cone
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl, filt_tanl
	from fundamentals   import fshift, rot_avg_image
	from projection     import prep_vol, prgs
	from utilities      import bcast_string_to_all, model_circle, get_arb_params, set_arb_params, drop_spider_doc
	from utilities      import get_image, drop_image, bcast_EMData_to_all, send_attr_dict, recv_attr_dict
	from utilities      import read_spider_doc, get_im
	from reconstruction import rec3D_MPI
	from statistics     import ccc
	from math           import pi, sqrt
	from string         import replace
	import os
	import sys
	#from development    import ali_G3
	from mpi 	    import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	    import mpi_barrier, mpi_gatherv, mpi_scatterv
	from mpi 	    import MPI_FLOAT, MPI_INT, MPI_SUM

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	if(myid == main_node):
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)
		from utilities import read_text_file
		pwpdb = read_text_file("pwpdb.txt", 1)
	mpi_barrier(MPI_COMM_WORLD)

	nima = EMUtil.get_image_count(stack)

	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)

	# figure the size of the chunk (3D is updated after each chunk).  Chunk should be given as 0.0< chunk <= 1.0.  1.0 means all projections
	if(chunk <= 0.0):  chunk = 1.0
	n_in_chunk  = max(int(chunk * (image_end-image_start+1)), 1)
	n_of_chunks = (image_end-image_start+1)//n_in_chunk + min((image_end-image_start+1)%n_in_chunk,1)
	
	outf = file(replace("progress%4d"%myid,' ','0'), "w")
	outf.write("  chunk = "+str(chunk)+"   ")
	outf.write("\n")
	outf.flush()
	outf.write("  n_in_chunk = "+str(n_in_chunk)+"   ")
	outf.write("  n_of_chunks = "+str(n_of_chunks)+"   ")
	outf.write("\n")
	outf.flush()
	#  Here we assume that reference volume exists
	vol = EMData()
	vol.read_image(ref_vol)
	nx  = vol.get_xsize()
	#ima = EMData()
	#ima.read_image(stack)
	#nx  = ima.get_xsize()
	if(ou <= 0):  ou = nx//2-2
	#if(myid == main_node):  vol.write_image(os.path.join(outdir,"ref_volf00.hdf"))
	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask3D = get_image(maskfile)
		else:                                   mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)
	mask2D = model_circle(ou, nx, nx)

	dataim = []
	#from utilities import read_spider_doc, set_arb_params
	#prm = read_spider_doc("params_new.doc")
	#prm_dict = ["phi", "theta", "psi", "s2x", "s2y"]
	#prm_dict = ["phi", "theta", "psi"]
	#from random import seed,gauss
	#seed()
	for im in xrange(image_start, image_end):
		ima = EMData()
		ima.read_image(stack, im)
		#angn = get_arb_params(ima, prm_dict)
		#for ian in xrange(3):  angn[ian] += gauss(0.0, 1.0)
		#set_arb_params(ima, angn, prm_dict)
		#set_arb_params(ima, prm[im], prm_dict)
		if(CTF):
			ctf_params = ima.get_attr('ctf')
			if(im == image_start): data_had_ctf = ima.get_attr('ctf_applied')
			if ima.get_attr('ctf_applied') == 0:
				st = Util.infomask(ima, mask2D, False)
				ima -= st[0]
				from filter import filt_ctf
				ima = filt_ctf(ima, ctf_params)
				ima.set_attr('ctf_applied', 1)
		dataim.append(ima)
	outf.write("  data read = "+str(image_start)+"   ")
	outf.write("\n")
	outf.flush()

	from utilities      import bcast_number_to_all
	from morphology     import threshold, threshold_to_minval

	par_str=["phi", "theta", "psi", "s2x", "s2y"]

	#  this is needed for gathering and scattering of cccfs
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(nima, number_of_proc, im)
		recvcount.append( ie - ib )

	from utilities import even_angles
	xrng = 0.0
	yrng = 0.0
	step = 1.0
	template_angles = even_angles(0.2, 0.0, 0.4, phiEqpsi = 'Zero', method = 'P')
	from filter import  filt_tophatb, filt_gaussl
	volftb = vol.copy()
	for iteration in xrange(maxit):
		outf.write("  iteration = "+str(iteration)+"   ")
		outf.write("\n")
		outf.flush()
		for  ic  in xrange(n_of_chunks):
			# compute updated 3D after each chunk
			
			outf.write("  generate projections = "+str(image_start)+"   ")
			outf.write("\n")
			outf.flush()
			# SORT PROJECTIONS
			volftb,kb  = prep_vol( filt_tophatb(volftb, 0.28, 0.44, False) )

			qcc = []
			for imn in xrange(image_start, image_end):
				atparams = get_arb_params(dataim[imn-image_start], par_str)
				projt = prgs(volftb, kb, [atparams[0], atparams[1], atparams[2], -atparams[3], -atparams[4]])
				qcc.append(ccc(projt, dataim[imn-image_start], mask2D))
				#qqcc=ccc(projt, dataim[imn-image_start], mask2D)
				#dataim[imn-image_start] /= (1.0-qqcc*qqcc)
			del projt
			del volftb

			recvbuf = mpi_gatherv(qcc, len(dataim), MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			del qcc
			mpi_barrier(MPI_COMM_WORLD)
			if(myid == main_node):
				templ = []
				for im in xrange(len(recvbuf)):    templ.append([float(recvbuf[im]), im])
				del recvbuf
				'''
				for im in xrange(len(templ)):
					outf.write('ccc, image %12.5f  %07d'%( templ[im][0], templ[im][1]  ))
					outf.write("\n")
				outf.flush()
				'''
				
				templ.sort()
				ilow = int(0.25*len(templ))  # reject 25% worst images.
				for im in xrange(ilow):  templ[im] = [templ[im][1], 0]
				for im in xrange(ilow, len(templ)):  templ[im] = [templ[im][1], 1]
				templ.sort()
				sendbuf = []
				for im in xrange(len(templ)):	sendbuf.append(templ[im][1])
				del templ
				'''
				qb = -1.0
				qs = 10.0
				for im in xrange(len(recvbuf)):
					qt = float(recvbuf[im])
					qb = max(qb,qt)
					qs = min(qs,qt)
				qs -= 1.0e-3
				qb -= qs
				sendbuf = []
				for im in xrange(len(recvbuf)):
					sendbuf.append((float(recvbuf[im])-qs)/qb)
				del recvbuf
				'''
			else:
				sendbuf = []
			mpi_barrier(MPI_COMM_WORLD)
			#recvbuf = mpi_scatterv(sendbuf, recvcount, disps, MPI_FLOAT, recvcount[myid], MPI_FLOAT, main_node, MPI_COMM_WORLD)
			recvbuf = mpi_scatterv(sendbuf, recvcount, disps, MPI_INT, recvcount[myid], MPI_INT, main_node, MPI_COMM_WORLD)
			del sendbuf

			for imn in xrange(image_start, image_end):
				dataim[imn-image_start].set_attr_dict({'active': int(recvbuf[imn-image_start])})
				#dataim[imn-image_start] /= float(recvbuf[imn-image_start])
			'''
			nact = 0
			for imn in xrange(image_start, image_end):
				nact += dataim[imn-image_start].get_attr('active')
			nact = float(nact)
			tn = mpi_reduce(nact, 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
			if(myid == main_node):
				outf.write('total number of used images %12.2f  '%(float(tn)))
				outf.write("\n")
				outf.flush()
			'''
			# resolution
			outf.write("  begin reconstruction = "+str(image_start)+"   ")
			outf.write("\n")
			outf.flush()
			vol, fscc = rec3D_MPI(dataim, snr, sym, mask3D, os.path.join(outdir, replace("resolution%3d_%3d"%(iteration, ic),' ','0') ), myid, main_node)
			outf.write("  done reconstruction = "+str(image_start)+"   ")
			outf.write("\n")
			outf.flush()
			volftb = vol.copy()

			#  restore original normalization
			#for imn in xrange(image_start, image_end):
			#	#dataim[imn-image_start].set_attr_dict({'active': int(recvbuf[imn-image_start])})
			#	dataim[imn-image_start] *= float(recvbuf[imn-image_start])
			del recvbuf
			mpi_barrier(MPI_COMM_WORLD)
			if(myid == main_node):
				drop_image(vol, os.path.join(outdir, replace("vol%3d_%3d.hdf"%(iteration, ic),' ','0') ))
				stat = Util.infomask(vol, mask3D, False)
				vol -= stat[0]
				vol /= stat[1]
				vol = threshold(vol)
				from  fundamentals  import  rops_table
				pwem = rops_table(vol)
				ftb = []
				for idum in xrange(len(pwem)):
					ftb.append(sqrt(pwpdb[idum]/pwem[idum]))
				from filter import filt_table, fit_tanh
				vol = filt_table(vol, ftb)
				del ftb, stat
				Util.mul_img(vol, mask3D)
				fl, aa = fit_tanh(fscc)
				vol = filt_tanl(vol, fl, aa)
				#vol = filt_gaussl(filt_tanl(vol, fl, aa),  0.2)
				outf.write('tanh params %8.4f  %8.4f '%(fl, aa))
				outf.write("\n")
				outf.flush()
				drop_image(vol, os.path.join(outdir, replace("volf%3d_%3d.hdf"%(iteration, ic),' ','0') ))
			#from sys import exit
			#exit()

			bcast_EMData_to_all(volftb, myid, main_node)
			bcast_EMData_to_all(vol, myid, main_node)
			volft,kb  = prep_vol(vol)

			image_start_in_chunk = image_start + ic*n_in_chunk
			image_end_in_chunk   = min(image_start_in_chunk + n_in_chunk, image_end)
			outf.write("ic "+str(ic)+"   image_start "+str(image_start)+"   n_in_chunk "+str(n_in_chunk)+"   image_end "+str(image_end)+"\n")
			outf.write("image_start_in_chunk "+str(image_start_in_chunk)+"  image_end_in_chunk "+str(image_end_in_chunk)+"\n")
			outf.flush()
			for imn in xrange(image_start_in_chunk, image_end_in_chunk):
				#from utilities import start_time, finish_time
				#t3=start_time()
				proj_ali_incore_cone(volft, kb, template_angles, dataim[imn-image_start], 1, ou, 1, xrng, yrng, step, outf)

			soto = []
			for imn in xrange(image_start, image_end):
				from utilities import set_params_proj, get_params_proj
				phi,theta,psi,s2x,s2y = get_params_proj( dataim[imn-image_start] )
				soto.append([phi,theta,psi,s2x,s2y,imn])
			drop_spider_doc(os.path.join(outdir, replace("new_params%3d_%3d_%3d"%(iteration, ic, myid),' ','0')), soto," phi, theta, psi, s2x, s2y, image number")
			del soto

			#from sys import exit
			#exit()
		#  here we should write header info, just in case the program crashes...
	del vol
	del volft
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)


"""
"""
def proj_ali_incore_index(volref, iref, mask3D, projdata, first_ring, last_ring, rstep, xrng, yrng, step, delta, ref_a, symmetry, MPI):
	from utilities    import even_angles, model_circle, compose_transform2, get_params_proj, set_params_proj
	from alignment    import prepare_refprojs
	#  DO NOT USE THIS ONE< WILL BE OBSOLETED SOON  PAP 01/25/08
	mode    = "F"
	# generate list of Eulerian angles for reference projections
	#  phi, theta, psi
	ref_angles = even_angles(delta, symmetry = symmetry, method = ref_a, phiEqpsi = "Minus")
	#from utilities import drop_spider_doc
	#drop_spider_doc("angles.txt",ref_angles)
	#  begin from applying the mask, i.e., subtract the average outside the mask and multiply by the mask
	if(mask3D):
		[mean, sigma, xmin, xmax ] =  Util.infomask(volref, mask3D, False)
		volref -= mean
		Util.mul_img(volref, mask3D)
	#drop_image(volref, "volref.spi", "s")
	#exit()
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
	#soto = []
	for imn in xrange(len(projdata)):
		peako = projdata[imn].get_attr('peak')
		from utilities import set_params_proj, get_params_proj
		phi,theta,psi,sxo,syo = get_params_proj( projdata[imn] )
		[ang, sxs, sys, mirror, nref, peak] = Util.multiref_polar_ali_2d(projdata[imn].process("normalize.mask", {"mask":mask2D, "no_sigma":1}), ref_proj_rings, xrng, yrng, step, mode, numr, cnx-sxo, cny-syo)
		if(peak > peako):
			numref=int(nref)
			projdata[imn].set_attr_dict({'peak':peak})
			projdata[imn].set_attr_dict({'group':iref})
			#[ang,sxs,sys,mirror,peak,numref] = apmq(projdata[imn], ref_proj_rings, xrng, yrng, step, mode, numr, cnx-sxo, cny-syo)
			#ang = (ang+360.0)%360.0
			# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
			#  What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
			angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1,  -ang, 0.,0.,1)
			if  mirror:
                		phi = (ref_angles[numref][0]+540.0)%360.0
                		theta = 180.0-ref_angles[numref][1]
                		psi = (540.0-ref_angles[numref][2]+angb)%360.0
                		s2x = sxb + sxo
                		s2y = syb + syo
                	else:
                		phi = ref_angles[numref][0]
                		theta = ref_angles[numref][1]
                		psi = (ref_angles[numref][2]+angb+360.0)%360.0
                		s2x = sxb + sxo
                		s2y = syb + syo
			from utilities import set_params_proj, get_params_proj
                	set_params_proj( projdata[imn], [phi, theta, psi, s2x, s2y] )
			#soto.append( [phi, theta,psi,s2x, s2y, mirror, numref, peak] )
	#from utilities import drop_spider_doc
	#drop_spider_doc("ali_s_params.txt",soto)

def proj_ali_incore_localB(volref, mask3D, projdata, first_ring, last_ring, rstep, xrng, yrng, step, delta, an, ref_a, symmetry, info=None, MPI=False):
	#This is for Berlin only
	from utilities    import even_angles, model_circle, compose_transform2, bcast_EMData_to_all
	from alignment    import prepare_refprojs
	from math         import cos, sin, pi
	qv = pi/180.
        
	mode    = "F"
	# generate list of Eulerian angles for reference projections
	#  phi, theta, psi
	ref_angles = even_angles(delta, symmetry = symmetry, method = ref_a, phiEqpsi = "Minus")
	#from utilities import drop_spider_doc
	#drop_spider_doc("angles.txt",ref_angles)
	#  begin from applying the mask, i.e., subtract the average outside the mask and multiply by the mask
	if(mask3D):
		[mean, sigma, xmin, xmax ] =  Util.infomask(volref, mask3D, False)
		volref -= mean
		Util.mul_img(volref, mask3D)
	#drop_image(volref, "volref.spi", "s")
	#exit()
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

	for i in xrange(len(ref_angles)):
		n1 = sin(ref_angles[i][1]*qv)*cos(ref_angles[i][0]*qv)
		n2 = sin(ref_angles[i][1]*qv)*sin(ref_angles[i][0]*qv)
		n3 = cos(ref_angles[i][1]*qv)
		ref_proj_rings[i].set_attr_dict( {"n1":n1, "n2":n2, "n3":n3} )

	ant = abs(cos(an*qv))
	for imn in xrange(len(projdata)):
		from utilities import set_params_proj, get_params_proj
		phi,theta,psi,sxo,syo = get_params_proj( projdata[imn] )
		if not(info is None):
			info.write( "prj %4d old params: %8.3f %8.3f %8.3f %8.3f %8.3f\n" %(imn, phi, theta, psi, sxo, syo) )
			info.flush()
		# This is for Berlin only
		from utilities import get_arb_params
		ctf_params = projdata[imn].get_attr('ctf')
		from morphology import ctf_2
		ctf2 = ctf_2(nx, ctf_params)
		nct = len(ctf2)
		from math import exp
		envt = []
		for i in xrange(nct):
			# e(x)*h(x)/(bckg(x)+e(x)**2*ctf(x)**2/20)
			xs = float(i)/2.22/nx
			et = exp(-70.0*xs**2)
			bckgt = exp(-0.8-120.*xs**2)+0.01
			ht = 1.0-0.6*exp(-xs**2/2.0/0.012**2)
			fmt = et/(bckgt + ctf2[i]*et**2/10.0)*ht
			envt.append(fmt)
		from filter import filt_table
		ima = filt_table(projdata[imn], envt)
		#from utilities import drop_spider_doc
		#if(myid == main_node):
		#	if(im == image_start):  drop_spider_doc("matc.doc", envt)
		[ang, sxs, sys, mirror, nref, peak] = Util.multiref_polar_ali_2d_local(ima.process("normalize.mask", {"mask":mask2D, "no_sigma":1}), ref_proj_rings, xrng, yrng, step, ant, mode, numr, cnx-sxo, cny-syo)
		#[ang, sxs, sys, mirror, nref, peak] = Util.multiref_polar_ali_2d_local(projdata[imn].process("normalize.mask", {"mask":mask2D, "no_sigma":1}), ref_proj_rings, xrng, yrng, step, ant, mode, numr, cnx-sxo, cny-syo)
		numref=int(nref)
		#[ang,sxs,sys,mirror,peak,numref] = apmq_local(projdata[imn], ref_proj_rings, xrng, yrng, step, ant, mode, numr, cnx-sxo, cny-syo)
		#ang = (ang+360.0)%360.0
		if(numref > -1):
			# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
			#  What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
			angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1,  -ang, 0.,0.,1)
			if  mirror:
				phi = (ref_angles[numref][0]+540.0)%360.0
				theta = 180.0-ref_angles[numref][1]
				psi = (540.0-ref_angles[numref][2]+angb)%360.0
				s2x = sxb + sxo
				s2y = syb + syo
			else:
				phi = ref_angles[numref][0]
				theta = ref_angles[numref][1]
				psi = (ref_angles[numref][2]+angb+360.0)%360.0
				s2x = sxb+sxo
				s2y = syb+syo

			from utilities import set_params_proj, get_params_proj
			set_params_proj( projdata[imn], [phi, theta, psi, s2x, s2y])

			# if -1 local search did not have any neighbors, simply skip it
			if not(info is None):
				info.write( "prj %4d new params: %8.3f %8.3f %8.3f %8.3f %8.3f\n" %(imn, phi, theta, psi, s2x, s2y) )
				info.flush()
	return

def proj_ali_incore_cone(volref, kb, template_angles, projdata, first_ring, last_ring, rstep, xrng, yrng, step, finfo=None):
	#alignment within a cone, no mirror considered

	from utilities    import even_angles, model_circle, compose_transform2, print_msg
	from alignment    import refprojs
	mode    = "F"
	#  Volume is prepared earlier
	nx   = projdata.get_xsize()
	ny   = projdata.get_ysize()
	#  center is in SPIDER convention
	cnx  = nx//2 + 1
	cny  = ny//2 + 1
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr   = ringwe(numr ,mode)

	# prepare 2-D mask for normalization
	mask2D = model_circle(last_ring, nx, ny)
	if(first_ring > 0): mask2D -= model_circle(first_ring, nx, ny)
	#  Generate ref_angles according to angles of the current projection

	from utilities import set_params_proj, get_params_proj
	phi, theta, psi, sxo, syo = get_params_proj( projdata )
	R2  = Transform({"type":"spider", "phi":phi,"theta":theta,"psi":0.0})
	ref_angles = []
	for i in xrange(len(template_angles)):
		R1  = Transform({"type":"spider", "phi":template_angles[i][0], "theta":template_angles[i][1], "psi":0.0})
		RR = R1*R2
		Euler = RR.get_rotation("spider")
		ref_angles.append( [ Euler['phi'], Euler['theta'], 0.0])

	# generate reference projections in polar coords
	ref_proj_rings = refprojs( volref, kb, ref_angles, last_ring, mask2D, cnx, cny, numr, mode, wr )

	#if(imn%10 == 0):  print_msg("%d  "%(imn))
	from utilities import set_params_proj, get_params_proj
	phi,theta,psi,sxo,syo = get_params_proj( projdata )
	if not(finfo is None):
		finfo.write( "old params: %8.3f %8.3f %8.3f %8.3f %8.3f\n" %(phi, theta, psi, sxo, syo) )
		finfo.flush()
	[ang, sxs, sys, nref, peak] = Util.multiref_polar_ali_2d_nom(projdata.process("normalize.mask", {"mask":mask2D, "no_sigma":1}), ref_proj_rings, xrng, yrng, step, mode, numr, cnx-sxo, cny-syo)
	numref=int(nref)
	#[ang,sxs,sys,mirror,peak,numref] = apmq(projdata, ref_proj_rings, xrng, yrng, step, mode, numr, cnx-sxo, cny-syo)
	#ang = (ang+360.0)%360.0
	# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
	#  What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
	angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1,  -ang, 0.,0.,1)
        phi   = ref_angles[numref][0]
        theta = ref_angles[numref][1]
        psi   = (ref_angles[numref][2]+angb+360.0)%360.0
        s2x   = sxb + sxo
        s2y   = syb + syo
	projdata.set_attr( "peak", peak )

	from utilities import set_params_proj, get_params_proj
	set_params_proj( projdata, [phi, theta, psi, s2x, s2y] )
        if not(finfo is None):
		finfo.write( "new params: %8.3f %8.3f %8.3f %8.3f %8.3f\n" %(phi, theta, psi, s2x, s2y) )
		finfo.flush()

###############################################################################################
## OLD CODE K-MEANS STABLE ####################################################################
# 2009-07-29 14:04:48 All should be remove after a while, until to be sure nothing
#                     is missing in the new code. JB.

# 2008-12-18 11:43:11

# K-means main stability
def k_means_stab(stack, outdir, maskname, opt_method, K, npart = 5, CTF = False, F = 0, FK = 0, maxrun = 50, th_nobj = 0, th_stab_max = 6.0, th_stab_min = 3.0, th_dec = 5, rand_seed = 0, MPI = False, CUDA = False):
	if MPI:
		k_means_stab_MPI(stack, outdir, maskname, opt_method, K, npart, CTF, F, FK, maxrun, th_nobj, th_stab_max, th_stab_min, th_dec, rand_seed)
		return
	elif CUDA:
		k_means_stab_CUDA(stack, outdir, maskname, K, npart, F, FK, maxrun, th_nobj, th_stab_max, th_stab_min, th_dec, rand_seed)
		return
	
	from utilities 	 import print_begin_msg, print_end_msg, print_msg, file_type
	from statistics  import k_means_stab_update_tag, k_means_headlog, k_means_open_unstable, k_means_stab_gather
	from statistics  import k_means_classical, k_means_SSE, k_means_SA_T0, k_means_stab_init_tag
	from statistics  import k_means_headlog, k_means_stab_asg2part, k_means_stab_pwa, k_means_stab_export
	import sys, logging, os

	# set params
	trials   = 1
	maxit    = 1e9
	critname = ''

	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', " ", 1)
	os.mkdir(outdir)

	# create main log
	f = open(outdir + '/main_log.txt', 'w')
	f.close()
	logging.basicConfig(filename = outdir + '/main_log.txt', format = '%(asctime)s     %(message)s', level = logging.INFO)
	logging.info('::: Start k-means stability :::')

	# manage random seed
	rnd = []
	for n in xrange(1, npart + 1): rnd.append(n * (2**n) + rand_seed)
	logging.info('Init list random seed: %s' % rnd)

	# loop over run
	stb          = 6.0
	flag_run     = True
	num_run      = 0
	while flag_run:
		flag_cluster = False
		num_run += 1
		logging.info('RUN %02d K = %03d %s' % (num_run, K, 40 * '-'))

		# open unstable images
		logging.info('... Open images')
		[im_M, mask, ctf, ctf2, LUT, N] = k_means_open_unstable(stack, maskname, CTF)
		logging.info('... %d unstable images found' % N)
		if N < 2:
			logging.info('[STOP] Not enough images')
			num_run -= 1
			break
		if F != 0:
			try:
				T0, ct_pert = k_means_SA_T0(im_M, mask, K, rand_seed, [CTF, ctf, ctf2], F)
				logging.info('... Select first temperature T0: %4.2f (dst %d)' % (T0, ct_pert))
			except SystemExit:
				logging.info('[STOP] Not enough images')
				num_run -= 1
				break
		else: T0 = 0

		# loop over partition
		ALL_ASG = []
		print_begin_msg('k-means')
		for n in xrange(npart):
			logging.info('...... Start partition: %d' % (n + 1))
			k_means_headlog(stack, 'partition %d' % (n+1), opt_method, N, K, critname, maskname, trials, maxit, CTF, T0, F, rnd[n], 1)
			if   opt_method == 'cla':
				try:			[Cls, assign] = k_means_classical(im_M, mask, K, rnd[n], maxit, trials, [CTF, ctf, ctf2], F, T0, False)
				except SystemExit:	flag_cluster  = True
			elif opt_method == 'SSE':
				try:			[Cls, assign] = k_means_SSE(im_M, mask, K, rnd[n], maxit, trials, [CTF, ctf, ctf2], F, T0, False)
				except SystemExit:      flag_cluster  = True
			if flag_cluster:
				logging.info('[WARNING] Empty cluster')
				break
			
			ALL_ASG.append(assign)

		print_end_msg('k-means')

		if flag_cluster:
			num_run -= 1
			if   FK == 0: K -= th_dec
			else:
				newK = int(float(K) * FK)
				if (K - newK) < th_dec: K -= th_dec
				else: K = newK

			if K > 1:
				logging.info('[WARNING] Restart the run with K = %d' % K)
				continue
			else:
				logging.info('[STOP] Not enough number of clusters ')
				break

		# convert local assignment to absolute partition
		logging.info('... Convert local assign to abs partition')
		ALL_PART = k_means_stab_asg2part(ALL_ASG, LUT)

		# calculate the stability
		MATCH, STB_PART, CT_s, CT_t, ST, st = k_means_stab_pwa(ALL_PART, 100)
		logging.info('... Stability: %5.2f %% (%d objects)' % (st, sum(CT_s)))

		# manage the stability
		if stb < th_stab_max:
			if   FK == 0: K -= th_dec
			elif FK != 1:
				newK = int(float(K) * FK)
				if (K - newK) < th_dec: K -= th_dec
				else: K = newK

		if K < 2:
			logging.info('[STOP] Not enough number of clusters')
			break
		if stb < th_stab_min:
			logging.info('[WARNING] Stability too low, restart the run with K = %d' % K)
			num_run -= 1
			continue

		# export the stable class averages
		count_k, id_rejected = k_means_stab_export(STB_PART, stack, num_run, outdir, th_nobj, CTF)
		logging.info('... Export %i stable class averages: average_stb_run%02d.hdf (rejected %i images)' % (count_k, num_run, len(id_rejected)))
		
		# tag informations to the header
		logging.info('... Update info to the header')
		k_means_stab_update_tag(stack, ALL_PART, STB_PART, num_run, id_rejected)

		# stop if max run is reach
		if num_run >= maxrun:
			flag_run = False
			logging.info('[STOP] Max number of runs is reached (%d)' % maxrun)

	# merge all stable averages
	ct = k_means_stab_gather(num_run, maskname, outdir)
	logging.info('Gather and normalize all stable class averages: averages.hdf (%d images)' % ct)
	
	logging.info('::: END k-means stability :::')

# K-means main stability
def k_means_stab_CUDA(stack, outdir, maskname, K, npart = 5, F = 0, FK = 0, maxrun = 50, th_nobj = 0, th_stab_max = 6.0, th_stab_min = 3.0, th_dec = 5, rand_seed = 0):
	
	from utilities 	 import print_begin_msg, print_end_msg, print_msg
	from utilities   import model_blank, get_image, get_im
	from statistics  import k_means_cuda_init_open_im, k_means_cuda_open_im, k_means_stab_init_tag
	from statistics  import k_means_cuda_headlog, k_means_cuda_error, k_means_cuda_info

	from statistics  import k_means_stab_update_tag, k_means_stab_gather, k_means_stab_init_tag
	from statistics  import k_means_stab_asg2part, k_means_stab_H, k_means_stab_export
	import sys, logging, os, pickle

	# create a directory
	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', " ", 1)
	os.mkdir(outdir)

	# create main log
	f = open(outdir + '/main_log.txt', 'w')
	f.close()
		
	logging.basicConfig(filename = outdir + '/main_log.txt', format = '%(asctime)s     %(message)s', level = logging.INFO)
	logging.info('::: Start k-means stability :::')

	# manage random seed
	rnd = []
	for n in xrange(1, npart + 1): rnd.append(n * (2**n) + rand_seed)
	logging.info('Init list random seed: %s' % rnd)

	# loop over run
	flag_run     = True
	num_run      = 0
	maxit        = int(1e9)
	T0           = float(-1) # auto
	while flag_run:
		flag_cluster = False
		num_run += 1
		logging.info('RUN %02d K = %03d %s' % (num_run, K, 40 * '-'))

		# create k-means obj
		KmeansCUDA = CUDA_kmeans()

		# open unstable images
		logging.info('... Open images')
		LUT, mask, N, m = k_means_cuda_init_open_im(stack, maskname)
		logging.info('... %d unstable images found' % N)
		if N < 2:
			logging.info('[STOP] Not enough images')
			num_run -= 1
			break
		KmeansCUDA.setup(m, N, K, F, T0, maxit, 0)
		k_means_cuda_open_im(KmeansCUDA, stack, LUT, mask)

		# loop over partition
		ALL_ASG = []
		print_begin_msg('k-means')
		for n in xrange(npart):
			# info
			logging.info('...... Start partition: %d' % (n + 1))
			partdir = 'partition %d' % (n + 1)
			ncpu   = 1
			method = 'cla'
			k_means_cuda_headlog(stack, partdir, method, N, K, maskname, maxit, T0, F, rnd[n], ncpu, m)

			# classification
			KmeansCUDA.set_rnd(rnd[n])
			status = KmeansCUDA.kmeans()
			if   status == 0:
				pass
			elif status == 5 or status == 4:
				logging.info('[WARNING] Empty cluster')
				k_means_cuda_error(status)
				flag_cluster = True
				break
			else:
				k_means_cuda_error(status)
				logging.info('[ERROR] CUDA status %i' % status)
				sys.exit()

			# get back the partition and its infos
			PART = KmeansCUDA.get_partition()
			ALL_ASG.append(PART)
			INFO = KmeansCUDA.get_info()
			k_means_cuda_info(INFO)

			import cPickle
			f = open('part%i.pck' % n, 'w')
			cPickle.dump(PART, f)
			f.close()

		# end of classification
		print_end_msg('k-means')

		# destroy k-means
		del KmeansCUDA

		# flow control if empty cluster
		if flag_cluster:
			num_run -= 1
			if   FK == 0: K -= th_dec
			else:
				newK = int(float(K) * FK)
				if (K - newK) < th_dec: K -= th_dec
				else: K = newK

			if K > 1:
				logging.info('[WARNING] Restart the run with K = %d' % K)
				continue
			else:
				logging.info('[STOP] Not enough number of clusters ')
				break

		# convert local assignment to absolute partition
		logging.info('... Convert local assign to abs partition')
		ALL_PART = k_means_stab_asg2part(ALL_ASG, LUT)

		# glooton control
		try:
			cmd = open('control', 'r').readline()
			cmd = cmd.strip(' \n')
			if cmd == 'stop':
				logging.info('[STOP] request by the user')
				break
		except: pass

		# calculate the stability
		stb, nb_stb, STB_PART = k_means_stab_H(ALL_PART)
		logging.info('... Stability: %5.2f %% (%d objects)' % (stb, nb_stb))

		# security to stop k-means stable cuda properly
		if os.path.exists('control'):
			try:
				cmd = open('control', 'r').readline()
				cmd = cmd.strip(' \n')
				if cmd == 'stop':
					logging.info('[STOP] request by the user')

					# export the stable class averages
					count_k, id_rejected = k_means_stab_export(STB_PART, stack, num_run, outdir, th_nobj, CTF)
					logging.info('... Export %i stable class averages: average_stb_run%02d.hdf (rejected %i images)' % (count_k, num_run, len(id_rejected)))

					# tag informations to the header
					logging.info('... Update info to the header')
					k_means_stab_update_tag(stack, ALL_PART, STB_PART, num_run, id_rejected)
					try:
						import os
						os.system('rm control')
					except: pass

					break
			except: pass

		# manage the stability
		if stb < th_stab_max:
			if   FK == 0: K -= th_dec
			elif FK != 1:
				newK = int(float(K) * FK)
				if (K - newK) < th_dec: K -= th_dec
				else: K = newK

		if K < 2:
			logging.info('[STOP] Not enough number of clusters')
			break
		if stb < th_stab_min:
			logging.info('[WARNING] Stability too low, restart the run with K = %d' % K)
			num_run -= 1
			continue
	
		# export the stable class averages
		count_k, id_rejected = k_means_stab_export(STB_PART, stack, num_run, outdir, th_nobj, CTF)
		logging.info('... Export %i stable class averages: average_stb_run%02d.hdf (rejected %i images)' % (count_k, num_run, len(id_rejected)))

		# tag informations to the header
		logging.info('... Update info to the header')
		k_means_stab_update_tag(stack, ALL_PART, STB_PART, num_run, id_rejected)

		# stop if max run is reach
		if num_run >= maxrun:
			flag_run = False
			logging.info('[STOP] Max number of runs is reached (%d)' % maxrun)

	# merge stable averages
	ct = k_means_stab_gather(num_run, maskname, outdir)
	logging.info('Gather and normalize all stable class averages: averages.hdf (%d images)' % ct)
	
	logging.info('::: END k-means stability :::')


# K-means main stability
def k_means_stab_MPI(stack, outdir, maskname, opt_method, K, npart = 5, CTF = False, F = 0, FK = 0, maxrun = 50, th_nobj = 0, th_stab_max = 6.0, th_stab_min = 3.0, th_dec = 5, rand_seed = 0):
	from utilities 	 import print_begin_msg, print_end_msg, print_msg, file_type
	from statistics  import k_means_stab_update_tag, k_means_headlog, k_means_open_unstable_MPI, k_means_stab_gather
	from statistics  import k_means_cla_MPI, k_means_SSE_MPI, k_means_SA_T0_MPI, k_means_stab_init_tag
	from statistics  import k_means_headlog, k_means_stab_asg2part, k_means_stab_pwa, k_means_stab_export
	from mpi         import mpi_init, mpi_comm_size, mpi_comm_rank, mpi_barrier, MPI_COMM_WORLD
	from mpi         import mpi_bcast, MPI_FLOAT
	import sys, logging, os

	sys.argv  = mpi_init(len(sys.argv), sys.argv)
	nb_cpu    = mpi_comm_size(MPI_COMM_WORLD)
	myid      = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	mpi_barrier(MPI_COMM_WORLD)

	# set params
	trials   = 1
	maxit    = 1e9
	critname = ''

	if myid == main_node:
		if os.path.exists(outdir): ERROR('Output directory exists, please change the name and restart the program', " ", 1)
		os.mkdir(outdir)

	mpi_barrier(MPI_COMM_WORLD)
		
	logging.basicConfig(filename = outdir + '/main_log.txt', format = '%(asctime)s     %(message)s', level = logging.INFO)
	if myid == main_node: logging.info('::: Start k-means stability :::')

	# manage random seed
	rnd = []
	for n in xrange(1, npart + 1): rnd.append(n * (2**n) + rand_seed)
	if myid == main_node: logging.info('Init list random seed: %s' % rnd)

	# loop over run
	stb          = 6.0
	flag_run     = True
	num_run      = 0
	while flag_run:
		flag_cluster = False
		num_run += 1
		if myid == main_node: logging.info('RUN %02d K = %03d %s' % (num_run, K, 40 * '-'))

		# open unstable images
		if myid == main_node: logging.info('... Open images')
		[im_M, mask, ctf, ctf2, LUT, N, N_start, N_stop] = k_means_open_unstable_MPI(stack, maskname, CTF, nb_cpu, main_node, myid)
		if myid == main_node: logging.info('... %d unstable images found' % N)
		if N < nb_cpu:
			logging.info('[STOP] Node %02d - Not enough images' % myid)
			num_run -= 1
			break
		if F != 0:
			try:
				T0, ct_pert = k_means_SA_T0_MPI(im_M, mask, K, rand_seed, [CTF, ctf, ctf2], F, myid, main_node, N_start, N_stop)
				if myid == main_node: logging.info('... Select first temperature T0: %4.2f (dst %d)' % (T0, ct_pert))
			except SystemExit:
				if myid == main_node: logging.info('[STOP] Not enough images')
				num_run -= 1
				mpi_barrier(MPI_COMM_WORLD)
				break
		else: T0 = 0

		# loop over partition
		ALL_ASG = []
		if myid == main_node: print_begin_msg('k-means')
		for n in xrange(npart):
			flag_cluster = False
			if myid == main_node:
				logging.info('...... Start partition: %d' % (n + 1))
				k_means_headlog(stack, 'partition %d' % (n+1), opt_method, N, K, critname, maskname, trials, maxit, CTF, T0, F, rnd[n], nb_cpu)
			if   opt_method == 'cla':
				try: [Cls, assign] = k_means_cla_MPI(im_M, mask, K, rnd[n], maxit, trials, [CTF, ctf, ctf2], myid, main_node, N_start, N_stop, F, T0)
				except SystemExit:
					if myid == main_node: logging.info('[WARNING] Empty cluster')
					mpi_barrier(MPI_COMM_WORLD)
					flag_cluster = True
					break
			elif opt_method == 'SSE':
				try: [Cls, assign] = k_means_SSE_MPI(im_M, mask, K, rnd[n], maxit, trials, [CTF, ctf, ctf2], myid, main_node, nb_cpu, N_start, N_stop, F, T0)
				except SystemExit:
					if myid == main_node: logging.info('[WARNING] Empty cluster')
					mpi_barrier(MPI_COMM_WORLD)
					flag_cluster = True
					break
			
			ALL_ASG.append(assign)
				
		if myid == main_node: print_end_msg('k-means')

		if flag_cluster:
			num_run -= 1
			if   FK == 0: K -= th_dec
			else:
				newK = int(float(K) * FK)
				if (K - newK) < th_dec: K -= th_dec
				else: K = newK
				
			if K > 1:
				if myid == main_node: logging.info('[WARNING] Restart the run with K = %d' % K)
				mpi_barrier(MPI_COMM_WORLD)
				continue
			else:
				if myid == main_node: logging.info('[STOP] Not enough number of clusters ')
				mpi_barrier(MPI_COMM_WORLD)
				break
	
		if myid == main_node:
			# convert local assignment to absolute partition
			logging.info('... Convert local assign to abs partition')
			ALL_PART = k_means_stab_asg2part(ALL_ASG, LUT)

			# calculate the stability
			MATCH, STB_PART, CT_s, CT_t, ST, st = k_means_stab_pwa(ALL_PART, 50)
			logging.info('... Stability: %5.2f %% (%d objects)' % (st, sum(CT_s)))
			

		mpi_barrier(MPI_COMM_WORLD)
		stb = mpi_bcast(stb, 1, MPI_FLOAT, main_node, MPI_COMM_WORLD)
		stb = stb.tolist()[0]

		# manage the stability
		if stb < th_stab_max:
			if   FK == 0: K -= th_dec
			elif FK != 1:
				newK = int(float(K) * FK)
				if (K - newK) < th_dec: K -= th_dec
				else: K = newK

		if K < 2:
			if myid == main_node: logging.info('[STOP] Not enough number of clusters')
			mpi_barrier(MPI_COMM_WORLD)
			break
		if stb < th_stab_min:
			if myid == main_node: logging.info('[WARNING] Stability too low, restart the run with K = %d' % K)
			mpi_barrier(MPI_COMM_WORLD)
			num_run -= 1
			continue

		if myid == main_node:
			# export the stable class averages
			count_k, id_rejected = k_means_stab_export(STB_PART, stack, num_run, outdir, th_nobj, CTF)
			logging.info('... Export %i stable class averages: average_stb_run%02d.hdf (rejected %i images)' % (count_k, num_run, len(id_rejected)))
			
			# tag informations to the header
			logging.info('... Update info to the header')
			k_means_stab_update_tag(stack, ALL_PART, STB_PART, num_run, id_rejected)

		mpi_barrier(MPI_COMM_WORLD)

		# stop if max run is reach
		if num_run >= maxrun:
			flag_run = False
			if myid == main_node: logging.info('[STOP] Max number of runs is reached (%d)' % maxrun)

	if myid == main_node:
		# merge and clean all stable averages
		ct = k_means_stab_gather(num_run, maskname, outdir)
		logging.info('Gather and normalize all stable class averages: averages.hdf (%d images)' % ct)

		logging.info('::: END k-means stability :::')

	mpi_barrier(MPI_COMM_WORLD)


'''
# sxk_means_stability.py

import os
import global_def
from   global_def import *
from   optparse import OptionParser
import sys
def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack outdir <maskfile> --K=2 --nb_part=5 --opt_method='SSE' --CTF --F=0.9 --max_run=10 --th_nobj=10 --th_stab_max=6.0 --th_stab_min=3.0 --min_dec_K=5 --FK=0 --rand_seed=0 --MPI --CUDA"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--K",              type="int",          default=2,         help="Number of classes for k-means (default 2)")
	parser.add_option("--nb_part",        type="int",          default=5,         help="Number of partitions used to calculate the stability (default 5)")
	parser.add_option("--opt_method",     type='string',       default="SSE",     help="K-means method: SSE (default), cla")
	parser.add_option("--CTF",            action="store_true", default=False,     help="Perform classification using CTF information")
	parser.add_option("--F",              type="float",        default=0.0,       help="Factor to decrease temperature in simulate annealing, ex.: 0.95")
	parser.add_option("--FK",             type="float",        default=0,         help="Cooling factor to decrease K according the number of runs")
	parser.add_option("--max_run",        type="int",          default=50,        help="Maximum number of runs (default 50)")
	parser.add_option("--th_nobj",        type="int",          default=10,        help="Cleanning threshold, classes with number of images < th_nobj are removed (default 10)")
	parser.add_option("--th_stab_max",    type="float",        default=6.0,       help="Stability threshold where K stop to decrease (default 6.0)")
	parser.add_option("--th_stab_min",    type="float",        default=0.0,       help="Stability threshold required to start the next run in percent (default 3.0)")
	parser.add_option("--min_dec_K",      type="int",          default=5,         help="Minimum decrement value to K when the stability is to low (default 5)")
	parser.add_option("--rand_seed",      type="int",          default=0,         help="Random seed")
	parser.add_option("--MPI",            action="store_true", default=False,     help="MPI version")
	parser.add_option("--CUDA",           action="store_true", default=False,     help="CUDA version")
	
	(options, args) = parser.parse_args()
    	if len(args) < 2 or len(args) > 3:
				print "usage: " + usage
        			print "Please run '" + progname + " -h' for detailed options"
	elif(options.opt_method != "cla"  and options.opt_method != "SSE"):
			sys.stderr.write("ERROR: unknown method\n\n")
			sys.exit()
	else:
		if len(args) == 2: mask = None
		else:              mask = args[2]

		if options.K < 2:
			sys.stderr.write('ERROR: K must be > 1 group\n\n')
			sys.exit()

		if options.nb_part < 2:
			sys.stderr.write('ERROR: nb_part must be > 1 partition\n\n')
			sys.exit()

		from  applications  import  k_means_stab
		global_def.BATCH = True
		k_means_stab(args[0], args[1], mask, options.opt_method, options.K, options.nb_part, options.CTF, options.F, options.FK, options.max_run, options.th_nobj, options.th_stab_max, options.th_stab_min, options.min_dec_K, options.rand_seed, options.MPI, options.CUDA)
		global_def.BATCH = False

if __name__ == "__main__":
	        main()



'''



## END OLD CODE K-MEANS STABLE ################################################################
###############################################################################################
