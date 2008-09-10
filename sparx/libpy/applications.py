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

def ali2d_reduce(stack, outdir, maskfile=None, ir=1, ou=-1, rs=1, center=1, maxit=0, CTF=False, user_func_name="ref_ali2d", randomize = False):
	from fundamentals import resample
	from utilities    import model_circle, get_arb_params, set_arb_params, getImage, get_params2D, set_params_2D
	from applications import ali2d_a
	import os
	
	from utilities import print_begin_msg, print_end_msg, print_msg
	import	types
		
	print_begin_msg("ali2d_reduce")

	import user_functions
	user_func = user_functions.factory[user_func_name]

	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit);

	nima = EMUtil.get_image_count(stack)
	ima = EMData()
	ima.read_image(stack, 0, True)
	if CTF:
		if(ima.get_attr_default('ctf_applied', 2) > 0):
			ERROR("data cannot be ctf-applied","ali2d_ar",1)
		pixel_size = ima.get_attr('Pixel_size')
	nx = ima.get_xsize()

	downsmpl = []
	qt = max(min(nx//2-ou, nx//32), 1)
	downsmpl.append(qt)
	if(qt > 2):
		while qt > 0:
			qt = qt//2
			downsmpl.append(qt)
		del downsmpl[len(downsmpl)-1]
	if(downsmpl[len(downsmpl)-1] > 1):  downsmpl.append(1)

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Inner radius                : %i\n"%(first_ring))
	# default value for the last ring
	if (last_ring == -1):  last_ring = nx//2-2

	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Ring step                   : %i\n"%(rstep))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("Data with CTF               : %s\n"%(CTF))

	if os.path.exists(outdir):
		os.system('rm -rf '+outdir)
	os.mkdir(outdir)

	if maskfile:
		if(type(maskfile) is types.StringType):
			print_msg("Maskfile                    : %s\n\n"%(maskfile))
			mask=getImage(maskfile)
		else:
			print_msg("Maskfile                    : user provided in-core mask\n\n")
			mask = maskfile
	else :
		print_msg("Maskfile                    : default, a circle with radius %i\n\n"%(last_ring))
		mask = model_circle(last_ring, nx, nx)

	if(len(downsmpl) == 1):
		outdir_temp = os.path.join(outdir,"001")
		ali2d_a(stack, outdir_temp, maskfile, ir, ou, rs, xr="1.0", yr="1.0", ts="0.2", center=center, maxit=maxit, CTF=CTF, user_func_name=user_func_name)
	else:
		attributes = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor"]

		for nd in downsmpl:
			outdir_temp = os.path.join(outdir,"%03d"%(downsmpl.index(nd)))
			os.mkdir(outdir_temp)
			if(nd > 1):
				sub_rate = 1.0/float(nd)
				ou_temp = max(int(ou*sub_rate),2)
				if CTF:  npx = pixel_size/sub_rate
				print_msg("New level of reduction	  x: %i\n"%(nd))
				print_msg("Pixel size used		   : %f\n"%(npx))
				maskfile_temp = resample(mask, sub_rate, fit_to_fft = True, num_prime = 5)
				tnx = maskfile_temp.get_xsize()
				print_msg("Image size used		   : %i\n"%(tnx))
				stack_temp = os.path.join(outdir,"temp%03d.hdf"%(downsmpl.index(nd)))
				for im in xrange(nima):
					ima = EMData()
					ima.read_image(stack, im)
					t = get_params2D(ima)
					if CTF:
						prm = get_arb_params(ima, attributes)
						prm[0] = npx
					ima = resample(ima, sub_rate, fit_to_fft=True, num_prime=5)
					if(nd == downsmpl[0]):
						t[1] *= sub_rate
						t[2] *= sub_rate
						set_params2D(ima, t)
					else:
						imold = EMData()
						imold.read_image(stack_previous ,im, True)
						prm_previous = get_params2D(imold)
						prm_previous[1] *= sub_rate/sub_rate_previous
						prm_previous[2] *= sub_rate/sub_rate_previous
						set_params2D(ima, prm_previous)
					if CTF: set_arb_params(ima, prm, attributes)
					ima.write_image(stack_temp, im)
				if(nd != downsmpl[0]):  os.system('rm -rf '+stack_previous)
			else:
				print_msg("Original image size     : %i\n"%(nx))
				print_msg("Pixel size              : %f\n"%(pixel_size))
				stack_temp = stack
				sub_rate = 1.0
				maskfile_temp = maskfile
				ou_temp = ou
				if(sub_rate_previous != 1.0):
					for im in xrange(nima):
						ima = EMData()
						ima.read_image(stack, im, True)
						imold = EMData()
						imold.read_image(stack_previous ,im, True)
						prm_previous = get_params2D(imold)
						prm_previous[1] *= sub_rate_previous
						prm_previous[2] *= sub_rate_previous
						set_params2D(ima, prm_previous)
						ima.write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)
					os.system('rm -rf '+stack_previous)

			ali2d_a(stack_temp, outdir_temp, maskfile_temp, ir, ou_temp, rs, xr="1.0", yr="1.0", ts="0.2", center=center, maxit=maxit, CTF=CTF, user_func_name=user_func_name)
			stack_previous = stack_temp
			sub_rate_previous = sub_rate

	print_end_msg("ali2d_reduce")

def ali2d_a(stack, outdir, maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", center=1, maxit=0, CTF=False, user_func_name="ref_ali2d", random_method="", F=0.996, MPI=False):

	if MPI:
		ali2d_a_MPI(stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, center, maxit, CTF, user_func_name, random_method, F)
		return
	from utilities    import model_circle, combine_params2, dropImage, getImage, get_arb_params, get_input_from_string
	from statistics   import add_oe_ave_varf
	from alignment    import Numrinit, ringwe, ali2d_s
	from filter       import filt_tophatb
	from morphology   import ctf_2
	from random       import random
	from utilities    import print_begin_msg, print_end_msg, print_msg
	from fundamentals import fft, rot_avg_table
	import os
		
	print_begin_msg("ali2d_a")

	if os.path.exists(outdir): os.system('rm -rf '+outdir)
	os.mkdir(outdir)

	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit);

	if max_iter == 0:
		max_iter  = 10
		auto_stop = True
	else:
		auto_stop = False

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)

	import user_functions
	user_func = user_functions.factory[user_func_name]

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Inner radius                : %i\n"%(first_ring))

	ima = EMData()
	ima.read_image(stack, 0, True)
	if CTF:
		if(ima.get_attr_default('ctf_applied', 2) > 0):
			ERROR("data cannot be ctf-applied","ali2d_a",1)
	nx = ima.get_xsize()
	# default value for the last ring
	if last_ring == -1:  last_ring = nx//2-2

	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Ring step                   : %i\n"%(rstep))
	print_msg("X search range              : %s\n"%(xrng))
	print_msg("Y search range              : %s\n"%(yrng))
	print_msg("Translational step          : %s\n"%(step))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("Data with CTF               : %s\n"%(CTF))
	print_msg("Simulated Annealing         : %s\n"%(randomize))
	if randomize: print_msg("Cooling Rate                : %f\n"%(F))
	if auto_stop:   print_msg("Stop iteration with         : criterion\n")
	else:           print_msg("Stop iteration with         : maxit\n")

	if maskfile:
		import	types
		if(type(maskfile) is types.StringType):
			print_msg("Maskfile                    : %s\n\n"%(maskfile))
			mask=getImage(maskfile)
		else:
			print_msg("Maskfile                    : user provided in-core mask\n\n")
			mask = maskfile
	else :
		print_msg("Maskfile                    : default, a circle with radius %i\n\n"%(last_ring))
		mask = model_circle(last_ring, nx, nx)

	cnx = int(nx/2)+1
 	cny = cnx
 	mode = "F"

	if CTF:
		from morphology   import ctf_img
		parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor"]
		ctf_2_sum = EMData(nx, nx, 1, False)
	else:
		ctf_2_sum = None

	del ima
	data = EMData.read_images(stack)
	for im in xrange(nima):
		data[im].set_attr('ID', im)
		if CTF:
			st = Util.infomask(data[im], mask, False)
			data[im] -= st[0]
			ctf_params = get_arb_params(data[im], parnames)
	 		Util.add_img2(ctf_2_sum, ctf_img(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5]))

	numr = Numrinit(first_ring, last_ring, rstep, mode) 	#precalculate rings
 	wr = ringwe(numr, mode)
	a0 = -1.0e22

	# initialize data for the reference preparation function
	#  mask can be modified in user_function
	ref_data = []
	ref_data.append( mask )
	ref_data.append( center )
	ref_data.append( None )
	ref_data.append( None )
	
	cs = [0.0]*2
	total_iter = 0
	for N_step in xrange(len(xrng)):
		msg = "\nX range = %5.2f   Y range = %5.2f   Step = %5.2f\n"%(xrng[N_step], yrng[N_step], step[N_step])
		print_msg(msg)
		for Iter in xrange(max_iter):
					
			tavg, vav, sumsq = add_oe_ave_varf(data, mask, mode="a", CTF=CTF, ctf_2_sum=ctf_2_sum)

			if total_iter==0 or randomize==False: select = 0
			else : 	
				select = 0
				for im in xrange(nima):
					this_select = data[im].get_attr("select")
					select += this_select

			total_iter += 1 
			
			if randomize and total_iter%10 == 0:
				dropImage(tavg, os.path.join(outdir, "aqc_%03d.hdf"%(total_iter)))
				dropImage(vav, os.path.join(outdir, "vav_%03d.hdf"%(total_iter)))

			tavg = fft(Util.divn_img(fft(tavg), vav))

			vav_r   = Util.pack_complex_to_real(vav)
			sumsq_r = Util.pack_complex_to_real(sumsq)
			rvar = rot_avg_table(vav_r)
			rsumsq = rot_avg_table(sumsq_r)
			frsc = []
			for i in xrange(len(rvar)):
				qt = max(0.0, rsumsq[i]/rvar[i] - 1.0)
				frsc.append([i/(len(rvar)-1)*0.5, qt/(qt+1)])

			ref_data[2] = tavg
			ref_data[3] = frsc
			#  call user-supplied function to prepare reference image, i.e., center and filter it
			tavg, cs = user_func( ref_data )
			
			Util.div_filter(sumsq, vav)
			sumsq = filt_tophatb(sumsq, 0.01, 0.49)
			a1 = Util.infomask(sumsq, None, True)
			a1 = a1[0]

			msg = "ITERATION   #%5d    criterion = %15.7e    average select = %5.3f\n\n"%(total_iter, a1, float(select)/nima)
			print_msg(msg)
			# write the current average
			if randomize and total_iter%10 == 0:
				dropImage(tavg, os.path.join(outdir, "aqf_%03d.hdf"%(total_iter)))
			# a0 should increase; stop algorithm when it decreases.    
			if(a1 < a0):
				if (auto_stop == True): break
			else:	a0 = a1
			if(N_step == (len(xrng)-1) and Iter == (max_iter-1)):  break
			ali2d_s(data, numr, wr, cs, tavg, cnx, cny, xrng[N_step], yrng[N_step], step[N_step], mode, range(nima), CTF=CTF, randomize=randomize, Iter=total_iter, F=F)
	# write out headers
	from utilities import write_headers
	write_headers(stack, data, range(nima))
	print_end_msg("ali2d_a")


def ali2d_a_MPI(stack, outdir, maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", center=1, maxit=0, CTF=False, user_func_name="ref_ali2d", random_method="", F=0.996):

	from utilities    import model_circle, combine_params2, dropImage, getImage, get_arb_params, get_input_from_string
	from utilities    import reduce_EMData_to_root, bcast_EMData_to_all, send_attr_dict, recv_attr_dict
	from statistics   import add_oe_ave_varf_MPI, add_oe_ave_varf_ML_MPI
	from alignment    import Numrinit, ringwe, ali2d_s
	from filter       import filt_tophatb
	from morphology   import ctf_2
	from numpy        import reshape, shape
	from utilities    import print_msg, print_begin_msg, print_end_msg
	from fundamentals import fft, rot_avg_table
	import os

	from mpi 	  import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	  import mpi_reduce, mpi_bcast, mpi_barrier
	from mpi 	  import MPI_SUM, MPI_FLOAT, MPI_INT

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	if myid == main_node:
		print_begin_msg("ali2d_a_MPI")
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)

	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit);
	
	if max_iter == 0:
		max_iter = 10
		auto_stop = True
	else:
		auto_stop = False

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	
	if myid == main_node:
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Inner radius                : %i\n"%(first_ring))

	nima = EMUtil.get_image_count(stack)	
	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)	
	ima = EMData()
	ima.read_image(stack, image_start, True)
	if CTF:
		if(ima.get_attr_default('ctf_applied', 2) > 0):
			ERROR("data cannot be ctf-applied","ali2d_a_MPI",1)
	nx = ima.get_xsize()
	# default value for the last ring
	if last_ring == -1: last_ring = nx//2-2

	if myid == main_node:
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Ring step                   : %i\n"%(rstep))
		print_msg("X search range              : %s\n"%(xrng))
		print_msg("Y search range              : %s\n"%(yrng))
		print_msg("Translational step          : %s\n"%(step))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Data with CTF               : %s\n"%(CTF))
		if random_method != "": 	
			print_msg("Random method               : %s\n"%(random_method))
		if random_method == "SA": 
			print_msg("Cooling Rate                : %f\n"%(F))
		if auto_stop: print_msg("Stop iteration with         : criterion\n")
		else:         print_msg("Stop iteration with         : maxit\n")
		print_msg("Number of processors used   : %d\n"%(number_of_proc))

	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  
			if myid == main_node:		print_msg("Maskfile                    : %s\n\n"%(maskfile))
			mask = getImage(maskfile)
		else:
			if myid == main_node: 		print_msg("Maskfile                    : user provided in-core mask\n\n")
			mask = maskfile
	else : 
		if myid==main_node: 	print_msg("Maskfile                    : default, a circle with radius %i\n\n"%(last_ring))
		mask = model_circle(last_ring,nx,nx)

	cnx  = nx//2+1
 	cny  = cnx
 	mode = "F"
	data = []
	if CTF:
		from morphology   import ctf_img
		parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
		ctf_2_sum = EMData(nx, nx, 1, False)

	for im in xrange(image_start, image_end):
		ima = EMData()
		ima.read_image(stack, im)
		if CTF:
			ctf_params = get_arb_params(ima, parnames)
			st = Util.infomask(ima, mask, False)
			ima -= st[0]
	 		Util.add_img2(ctf_2_sum, ctf_img(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5]))
		data.append(ima)
	
	numr = Numrinit(first_ring, last_ring, rstep, mode) 	#precalculate rings
 	wr = ringwe(numr, mode)
	a0 = -1.0e22
	
	if CTF: reduce_EMData_to_root(ctf_2_sum, myid, main_node)
	if myid == main_node:
		# initialize data for the reference preparation function
		ref_data = []
		ref_data.append( mask )
		ref_data.append( center )
		ref_data.append( None )
		ref_data.append( None )
		if CTF:	ctf_2sum = ctf_2_sum.copy()
	if CTF:	del ctf_2_sum	

	again = True
	cs = [0.0]*2
	total_iter = 0
	for N_step in xrange(len(xrng)):
		msg = "\nX range = %5.2f   Y range = %5.2f   Step = %5.2f\n"%(xrng[N_step], yrng[N_step], step[N_step])
		if myid == main_node: print_msg(msg)
		for Iter in xrange(max_iter):
		
			tavg, vav = add_oe_ave_varf_MPI(data, mask, mode="a", CTF=CTF)
			if random_method == "ML": 
				tavg_ML, vav_ML = add_oe_ave_varf_ML_MPI(data, mask, mode="a", CTF=CTF)
			#  bring all partial sums together
			reduce_EMData_to_root(tavg, myid, main_node)
			reduce_EMData_to_root(vav, myid, main_node)
			if random_method == "ML": 
				reduce_EMData_to_root(tavg_ML, myid, main_node)
				reduce_EMData_to_root(vav_ML, myid, main_node)			
			
			attr_list = data[0].get_attr_dict()	
			if attr_list.has_key("select") == False:
				select = 0
			else : 	
				select = 0
				for im in xrange(image_start, image_end):
					this_select = data[im-image_start].get_attr("select")
					select += this_select
			select  = mpi_reduce(select, 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
			
			total_iter += 1
			
			if myid == main_node:
				sumsq = fft(tavg)
				if random_method == "ML":
					sumsq_ML = fft(tavg_ML)				
				if CTF:	
					tavg = fft(Util.divn_img(sumsq, ctf_2sum))
				 	Util.mul_img(sumsq, sumsq.conjg())
				 	Util.div_img(sumsq, ctf_2sum)
			 		Util.sub_img(vav, sumsq)
					if random_method == "ML":
						tavg_ML = fft(Util.divn_img(sumsq_ML, ctf_2sum))
					 	Util.mul_img(sumsq_ML, sumsq_ML.conjg())
					 	Util.div_img(sumsq_ML, ctf_2sum)
			 			Util.sub_img(vav_ML, sumsq_ML)					
				else:
					Util.mul_scalar(tavg, 1.0/float(nima))
				 	Util.mul_img(sumsq, sumsq.conjg())
					Util.mul_scalar(sumsq, 1.0/float(nima))
					Util.sub_img(vav, sumsq)
					if random_method == "ML":
						Util.mul_scalar(tavg_ML, 1.0/float(nima))
					 	Util.mul_img(sumsq_ML, sumsq_ML.conjg())
						Util.mul_scalar(sumsq_ML, 1.0/float(nima))
						Util.sub_img(vav_ML, sumsq_ML)
												
				Util.mul_scalar(vav, 1.0/float(nima-1))
				if random_method == "ML":
					Util.mul_scalar(vav_ML, 1.0/float(nima-1))				

				if random_method=="" or total_iter%10 == 0:
					dropImage(tavg, os.path.join(outdir, "aqc_%03d.hdf"%(total_iter)))
					dropImage(vav, os.path.join(outdir, "vav_%03d.hdf"%(total_iter)))

				tavg = fft(Util.divn_img(fft(tavg), vav))
				if random_method == "ML":
					tavg_ML = fft(Util.divn_img(fft(tavg_ML), vav_ML))
					
 				if random_method != "ML":
					vav_r   = Util.pack_complex_to_real(vav)
		 		 	sumsq_r = Util.pack_complex_to_real(sumsq)
					rvar = rot_avg_table(vav_r)
					rsumsq = rot_avg_table(sumsq_r)
					frsc = []
					for i in xrange(len(rvar)):
						qt = max(0.0, rsumsq[i]/rvar[i] - 1.0)
						frsc.append([i/(len(rvar)-1)*0.5, qt/(qt+1)])
						
					ref_data[2] = tavg
					ref_data[3] = frsc
				else:
					vav_r   = Util.pack_complex_to_real(vav_ML)
		 		 	sumsq_r = Util.pack_complex_to_real(sumsq_ML)
					rvar = rot_avg_table(vav_r)
					rsumsq = rot_avg_table(sumsq_r)
					frsc = []
					for i in xrange(len(rvar)):
						qt = max(0.0, rsumsq[i]/rvar[i] - 1.0)
						frsc.append([i/(len(rvar)-1)*0.5, qt/(qt+1)])
						
					ref_data[2] = tavg_ML
					ref_data[3] = frsc
					
				#  call user-supplied function to prepare reference image, i.e., center and filter it
				tavg, cs = user_func( ref_data )

				Util.div_filter(sumsq, vav)
				sumsq = filt_tophatb(sumsq, 0.01, 0.49)
				a1 = Util.infomask(sumsq, None, True)
				a1 = a1[0]

				msg = "ITERATION   #%5d    criterion = %15.7e    average select = %5.3f\n\n"%(total_iter, a1, float(select)/nima)
				print_msg(msg)
				# write the current average
				if random_method=="" or total_iter%10 == 0:
					dropImage(tavg, os.path.join(outdir, "aqf_%03d.hdf"%(total_iter)))
				# a0 should increase; stop algorithm when it decreases.    
				if a1 < a0:
					if (auto_stop == True): break
				else:	a0 = a1
			else:
				tavg = EMData(nx, nx, 1, True)
				cs = [0.0]*2

			again = mpi_bcast(again, 1, MPI_INT, main_node, MPI_COMM_WORLD)
			bcast_EMData_to_all(tavg, myid, main_node)
			if not again: break
			if N_step == len(xrng)-1 and Iter == max_iter-1:  break
			cs = mpi_bcast(cs, 2, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			ali2d_s(data, numr, wr, cs, tavg, cnx, cny, xrng[N_step], yrng[N_step], step[N_step], mode, CTF=CTF, random_method=random_method, Iter=total_iter, F=F)
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
	#if(CTF and data_had_ctf == 0):
	#	for im in xrange(len(data)):  data[im].set_attr('ctf_applied', 0)
	par_str = ["xform.align2d"]
	if myid == main_node: recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	else: send_attr_dict(main_node, data, par_str, image_start, image_end)
	if myid == main_node:  print_end_msg("ali2d_a_MPI")


"""
#* ali2d_c
def ali2d_c(stack, outdir, maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", center=1, maxit=0, CTF=False, snr=1.0, user_func_name="ref_ali2d", rand_alpha = False):
	#** prepare
	from global_def import MPI
	if MPI:
		ali2d_c_MPI(stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, center, maxit, CTF, snr, user_func_name, rand_alpha)
		return

	from utilities    import model_circle, combine_params2, dropImage, getImage, get_arb_params, get_input_from_string
	from utilities    import set_arb_params
	from statistics   import aves, ave_oe_series, fsc_mask, add_oe_series
	from alignment    import Numrinit, ringwe, ali2d_s
	from filter       import filt_ctf, filt_table
	from morphology   import ctf_2
	from random       import random
	import os
	
	from utilities import print_begin_msg, print_end_msg, print_msg
	import	types
		
	print_begin_msg("ali2d_c")

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
	ima.read_image(stack, 0, True)
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
			mask=getImage(maskfile)
		else:
			print_msg("Maskfile                    : user provided in-core mask\n\n")
			mask = maskfile
	else :
		print_msg("Maskfile                    : default, a circle with radius %i\n\n"%(last_ring))
		mask = model_circle(last_ring, nx, nx)

	cnx = int(nx/2)+1
 	cny = cnx
 	mode = "F"
	data = []
	od = []
	if(CTF):
		parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
		#                        0                1              2     3              4                5                6
		ctf_params = get_arb_params(ima, parnames)
		data_had_ctf = ctf_params[6]
		ctm = ctf_2(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
		lctf = len(ctm)
		ctf2 = []
		ctf2.append([0.0]*lctf)
		ctf2.append([0.0]*lctf)
		ctfb2 = [0.0]*lctf
	for im in xrange(nima):
		ima = EMData()
		ima.read_image(stack, im)
		if(CTF):
			ctf_params = get_arb_params(ima, parnames)
			ctm = ctf_2(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
			k = im%2
			for i in xrange(lctf):
				ctf2[k][i] += ctm[i]
			if(ctf_params[6] == 0):
				st = Util.infomask(ima, mask, False)
				ima -= st[0]
				from utilities import info
				print im
				info(ima)
				od.append(ima)
				#      filt_ctf(image,	     dz,		  cs,		   voltage,	     pixel_size,     amp_contrast=0.1,	  b_factor=0.0):
				ima = filt_ctf(ima, ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5], pad = True)
				ima.set_attr('ctf_applied', 1)
				info(ima)

		data.append(ima)
	
	if(CTF):
		for i in xrange(lctf):
			ctfb2[i] = 1.0/(ctf2[0][i] + ctf2[1][i] + 1.0/snr)
			for k in xrange(2):
				ctf2[k][i] = 1.0/(ctf2[k][i] + 1.0/snr)

	#** startup
	numr = Numrinit(first_ring, last_ring, rstep, mode) 	#precalculate rings
 	wr = ringwe(numr, mode)
	a0 = 1.0e22

	# initialize data for the reference preparation function
	#  mask can be modified in user_function
	ref_data = []
	ref_data.append( mask )
	ref_data.append( center )
	cs = [0.0]*2
	#** iterate
	total_iter = 0
	for N_step in xrange(len(xrng)):
		msg = "\nX range = %5.2f   Y range = %5.2f   Step = %5.2f\n"%(xrng[N_step], yrng[N_step], step[N_step])
		print_msg(msg)
		for Iter in xrange(max_iter):
					
			av1, av2 = add_oe_series(data)
			if(CTF):
				tavg = filt_table(Util.addn_img(av1, av2), ctfb2)
				av1  = filt_table(av1, ctf2[0])
				av2  = filt_table(av2, ctf2[1])
			else:
				tavg = (av1+av2)/nima
			total_iter += 1 

			from fundamentals import fft
			from statistics import varfctf
			frsc = fsc_mask(av1, av2, ref_data[0], 1.0, os.path.join(outdir, "drc%03d"%(total_iter)))

			dropImage(tavg, os.path.join(outdir, "aqc_%03d.hdf"%(total_iter)))
			par_str = ["alpha", "sx", "sy", "mirror"]
			for i in xrange(len(data)):
				set_arb_params(od[i],get_arb_params(data[i],par_str),par_str)
			vav, rf = varfctf(od, mask, "a")
			dropImage(vav, os.path.join(outdir, "vav_%03d.hdf"%(total_iter)))
			tavg = fft(Util.divn_img(fft(tavg), vav))
			dropImage(tavg, os.path.join(outdir, "aqv_%03d.hdf"%(total_iter)))

			ref_data.append( tavg )
			ref_data.append( frsc )
			#  call user-supplied function to prepare reference image, i.e., center and filter it
			tavg, cs = user_func( ref_data )
			del ref_data[2]
			del ref_data[2]
			# a0 should increase; stop algorithm when it decreases.    
			a1 = tavg.cmp("dot", tavg, dict(negative = 0, mask = ref_data[0]))
			msg = "ITERATION #%3d	     criterion = %20.7e\n"%(total_iter,a1)
			print_msg(msg)
			# write the current average
			dropImage(tavg, os.path.join(outdir, "aqf_%03d.hdf"%(total_iter)))
			if(a1 < a0):
				if (auto_stop == True): break
			else:	a0 = a1
			ali2d_s(data, numr, wr, cs, tavg, cnx, cny, xrng[N_step], yrng[N_step], step[N_step], mode, range(nima))
	# write out headers
	if(CTF and data_had_ctf == 0):
		for im in xrange(nima):
			data[im].set_attr('ctf_applied', 0)
	for im in xrange(nima):
		data[im].write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)
	print_end_msg("ali2d_c")


"""

#* ali2d_c
def ali2d_c(stack, outdir, maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", center=1, maxit=0, CTF=False, snr=1.0, user_func_name="ref_ali2d", rand_alpha = False, MPI=False):
	if MPI:
		ali2d_c_MPI(stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, center, maxit, CTF, snr, user_func_name, rand_alpha)
		return

	from utilities    import model_circle, combine_params2, dropImage, getImage, get_arb_params, get_input_from_string
	from statistics   import fsc_mask, add_oe_series
	from alignment    import Numrinit, ringwe, ali2d_s
	from filter       import filt_ctf, filt_table
	from morphology   import ctf_2
	from utilities    import print_begin_msg, print_end_msg, print_msg
	import os
		
	print_begin_msg("ali2d_c")

	if os.path.exists(outdir):   os.system('rm -rf '+outdir)
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

	nima = EMUtil.get_image_count(stack)
	ima = EMData()
	ima.read_image(stack, 0, True)
	nx = ima.get_xsize()
	# default value for the last ring
	if last_ring == -1:  last_ring = nx//2-2

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

	if maskfile:
		import	types
		if type(maskfile) is types.StringType:
			print_msg("Maskfile                    : %s\n\n"%(maskfile))
			mask=getImage(maskfile)
		else:
			print_msg("Maskfile                    : user provided in-core mask\n\n")
			mask = maskfile
	else :
		print_msg("Maskfile                    : default, a circle with radius %i\n\n"%(last_ring))
		mask = model_circle(last_ring, nx, nx)

	cnx = nx//2+1
 	cny = cnx
 	mode = "F"
	data = []
	if CTF:
		parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
		#                        0                1              2     3              4                5                6
		ctf_params = get_arb_params(ima, parnames)
		data_had_ctf = ctf_params[6]
		ctm = ctf_2(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
		lctf = len(ctm)
		ctf2 = []
		ctf2.append([0.0]*lctf)
		ctf2.append([0.0]*lctf)
		ctfb2 = [0.0]*lctf
	del ima
	data = EMData.read_images(stack)
	for im in xrange(nima):
		data[im].set_attr('ID', im)
		if CTF:
			ctf_params = get_arb_params(data[im], parnames)
			ctm = ctf_2(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
			k = im%2
			for i in xrange(lctf):
				# ctf2[k][i] += ctm[i]
				ctf2[k][i] = ctf2[k][i] + ctm[i]
			if ctf_params[6] == 0:
				st = Util.infomask(data[im], mask, False)
				data[im] -= st[0]
				data[im] = filt_ctf(data[im], ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
				data[im].set_attr('ctf_applied', 1)

	
	if CTF:
		for i in xrange(lctf):
			ctfb2[i] = 1.0/(ctf2[0][i] + ctf2[1][i] + 1.0/snr)
			for k in xrange(2):
				ctf2[k][i] = 1.0/(ctf2[k][i] + 1.0/snr)

	#** startup
	numr = Numrinit(first_ring, last_ring, rstep, mode) 	#precalculate rings
 	wr = ringwe(numr, mode)
	
	av1, av2 = add_oe_series(data)
	Util.add_img(av1, av2)
	if CTF:  tavg = filt_table(av1, ctfb2)
	else:     tavg = av1/nima
	a0 = tavg.cmp("dot", tavg, dict(negative = 0, mask = mask))
	msg = "Initial criterion = %-20.7e\n"%(a0)
	print_msg(msg)

	# initialize data for the reference preparation function
	#  mask can be modified in user_function
	ref_data = []
	ref_data.append( mask )
	ref_data.append( center )
	ref_data.append( None )
	ref_data.append( None )
	cs = [0.0]*2
	#** iterate
	total_iter = 0
	for N_step in xrange(len(xrng)):
		msg = "\nX range = %5.2f   Y range = %5.2f   Step = %5.2f\n"%(xrng[N_step], yrng[N_step], step[N_step])
		print_msg(msg)
		for Iter in xrange(max_iter):
					
			ali2d_s(data, numr, wr, cs, tavg, cnx, cny, xrng[N_step], yrng[N_step], step[N_step], mode, range(nima))
			av1, av2 = add_oe_series(data)
			if CTF:
				tavg = filt_table(Util.addn_img(av1, av2), ctfb2)
				av1  = filt_table(av1, ctf2[0])
				av2  = filt_table(av2, ctf2[1])
			else:
				tavg = (av1+av2)/nima
			total_iter += 1 

			dropImage(tavg, os.path.join(outdir, "aqc_%03d.hdf"%(total_iter)))

			frsc = fsc_mask(av1, av2, ref_data[0], 1.0, os.path.join(outdir, "drc%03d"%(total_iter)))

			ref_data[2] = tavg
			ref_data[3] = frsc
			#  call user-supplied function to prepare reference image, i.e., center and filter it
			tavg, cs = user_func( ref_data )

			# a0 should increase; stop algorithm when it decreases.    
			a1 = tavg.cmp("dot", tavg, dict(negative = 0, mask = ref_data[0]))
			msg = "ITERATION   #%5d	     criterion = %20.7e\n\n"%(total_iter,a1)
			print_msg(msg)
			# write the current average
			dropImage(tavg, os.path.join(outdir, "aqf_%03d.hdf"%(total_iter)))
			if a1 < a0:
				if auto_stop == True: break
			else:	a0 = a1
	# write out headers
	if CTF and data_had_ctf == 0:
		for im in xrange(nima):	data[im].set_attr('ctf_applied', 0)
	from utilities import write_headers
	write_headers(stack, data, range(nima))
	print_end_msg("ali2d_c")


def ali2d_c_MPI(stack, outdir, maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", center=1, maxit=0, CTF=False, snr=1.0, user_func_name="ref_ali2d", rand_alpha=False):

	from utilities    import model_circle, model_blank, combine_params2, dropImage, getImage, get_arb_params, get_input_from_string
	from utilities    import reduce_EMData_to_root, bcast_EMData_to_all, send_attr_dict, recv_attr_dict
	from statistics   import fsc_mask, add_oe_series
	from alignment    import Numrinit, ringwe, ali2d_s
	from filter       import filt_table, filt_ctf
	from morphology   import ctf_2
	from numpy        import reshape, shape
	from utilities    import print_msg, print_begin_msg, print_end_msg
	import os
	import sys
	from mpi 	  import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	  import mpi_reduce, mpi_bcast, mpi_barrier
	from mpi 	  import MPI_SUM, MPI_FLOAT, MPI_INT

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	if myid == main_node:
		print_begin_msg("ali2d_c_MPI")
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
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

	nima = EMUtil.get_image_count(stack)
	
	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
	
	ima = EMData()
	ima.read_image(stack, image_start, True)
	nx = ima.get_xsize()
	# default value for the last ring
	if last_ring == -1: last_ring = nx//2-2

	if myid == main_node:
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Ring step                   : %i\n"%(rstep))
		print_msg("X search range              : %s\n"%(xrng))
		print_msg("Y search range              : %s\n"%(yrng))
		print_msg("Translational step          : %s\n"%(step))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("data with CTF               : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		if auto_stop:
			print_msg("Stop iteration with         : criterion\n")
		else:
			print_msg("Stop iteration with         : maxit\n")
		print_msg("Number of processors used   : %d\n"%(number_of_proc))

		
	if maskfile:
		import  types
		if type(maskfile) is types.StringType:  
			if myid == main_node:		print_msg("Maskfile                    : %s\n\n"%(maskfile))
			mask=getImage(maskfile)
		else:
			if myid == main_node: 		print_msg("Maskfile                    : user provided in-core mask\n\n")
			mask = maskfile
	else : 
		if myid==main_node: 	print_msg("Maskfile                    : default, a circle with radius %i\n\n"%(last_ring))
		mask = model_circle(last_ring, nx, nx)

	cnx  = nx//2+1
 	cny  = cnx
 	mode = "F"
	data = []
	if CTF:
		parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
		#                        0                1              2     3              4              5               6
		ctf_params = get_arb_params(ima, parnames)
		data_had_ctf = ctf_params[6]
		ctm = ctf_2(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
		lctf = len(ctm)
		ctf2 = []
		ctf2.append([0.0]*lctf)
		ctf2.append([0.0]*lctf)
		# the next only required for the total average
		if myid == main_node:  ctfb2 = [0.0]*lctf

	del ima
	data = EMData.read_images(stack, range(image_start, image_end))
	for im in xrange(image_start, image_end):
		data[im-image_start].set_attr('ID', im)
		if CTF:
			ctf_params = get_arb_params(data[im-image_start], parnames)
			ctm = ctf_2(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
			k = im%2
			for i in xrange(lctf):  ctf2[k][i] += ctm[i]
			if ctf_params[6] == 0:
				st = Util.infomask(data[im-image_start], mask, False)
				data[im-image_start] -= st[0]
				from filter import filt_ctf
				data[im-image_start] = filt_ctf(data[im-image_start], ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
				data[im-image_start].set_attr('ctf_applied', 1)
		
	
	if CTF:
		# bring ctf2 together, keep them only on main node,  strange trick required because mpi_reduce changes the nformat to numarray
		s = shape(ctf2)
		ctf2  = mpi_reduce(ctf2, 2*lctf, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
		if myid == main_node:
			ctf2 = reshape(ctf2, s)
			ctf2n = []
			ctf2n.append([0.0]*lctf)
			ctf2n.append([0.0]*lctf)
			for i in xrange(lctf):
				ctfb2[i] = 1.0/(ctf2[0][i] + ctf2[1][i] + 1.0/snr)
				for k in xrange(2):
					ctf2n[k][i] = 1.0/(ctf2[k][i] + 1.0/snr)
		del  ctf2

	numr = Numrinit(first_ring, last_ring, rstep, mode) 	#precalculate rings
 	wr = ringwe(numr, mode)
	
	av1, av2 = add_oe_series(data)
	Util.add_img(av1, av2)
	tavg = model_blank(nx, nx)
	#  get the total sum
	reduce_EMData_to_root(av1, myid, main_node)

	if myid == main_node:
		if CTF: tavg = filt_table(av1, ctfb2)
		else:    tavg = av1/nima
		a0 = tavg.cmp("dot", tavg, {"negative":0, "mask":mask})
		msg = "Initial criterion = %-20.7e\n"%(a0)
		print_msg(msg)
		# initialize data for the reference preparation function
		ref_data = []
		ref_data.append( mask )
		ref_data.append( center )
		ref_data.append( None )
		ref_data.append( None )

	bcast_EMData_to_all(tavg, myid, main_node)

	again = True
	cs = [0.0]*2
	total_iter = 0
	for N_step in xrange(len(xrng)):
		msg = "\nX range = %5.2f   Y range = %5.2f   Step = %5.2f\n"%(xrng[N_step], yrng[N_step], step[N_step])
		if myid == main_node: print_msg(msg)
		for Iter in xrange(max_iter):
			ali2d_s(data, numr, wr, cs, tavg, cnx, cny, xrng[N_step], yrng[N_step], step[N_step], mode)
			av1, av2 = add_oe_series(data)
			#  bring all partial sums together
			reduce_EMData_to_root(av1, myid, main_node)
			reduce_EMData_to_root(av2, myid, main_node)
			if myid == main_node:
				if CTF:
					tavg = filt_table(Util.addn_img(av1, av2), ctfb2)
					av1  = filt_table(av1, ctf2n[0])
					av2  = filt_table(av2, ctf2n[1])
				else:
					tavg = (av1 + av2)/nima
				total_iter += 1

				dropImage(tavg, os.path.join(outdir, "aqc_%03d.hdf"%(total_iter)))

				frsc = fsc_mask(av1, av2, ref_data[0], 1.0, os.path.join(outdir, "drc%03d"%(total_iter)))
				ref_data[2] = tavg
				ref_data[3] = frsc
				#  call user-supplied function to prepare reference image, i.e., center and filter it
				tavg, cs = user_func( ref_data )

				# a0 should increase; stop algorithm when it decreases.    
				a1 = tavg.cmp("dot", tavg, dict(negative = 0, mask = ref_data[0]))
				msg = "ITERATION   #%5d	     criterion = %20.7e\n\n"%(total_iter,a1)
				print_msg(msg)
				# write the current average
				dropImage(tavg, os.path.join(outdir, "aqf_%03d.hdf"%(total_iter)))
				if a1 < a0:
					if auto_stop == True: break
				else:	a0 = a1

			again = mpi_bcast(again, 1, MPI_INT, main_node, MPI_COMM_WORLD)
			bcast_EMData_to_all(tavg, myid, main_node)
			cs = mpi_bcast(cs, 2, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			if not again: break
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
	if CTF and data_had_ctf == 0:
		for im in xrange(len(data)):  data[im].set_attr('ctf_applied', 0)
	par_str = ["xform.align2d"]
	if myid == main_node: recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	else: send_attr_dict(main_node, data, par_str, image_start, image_end)
	if myid == main_node:	print_end_msg("ali2d_c_MPI")


def ali2d_e(stack, outdir, maskfile = None, ou = -1, br = 1.75, center = 1, eps = 0.001, maxit = 10, CTF = False, snr = 1.0, user_func_name="ref_ali2d"):
# 2D alignment using amoeba and gridding interpolation
	from alignment    	import kbt
	from utilities    	import model_circle, amoeba, compose_transform2, dropImage, get_arb_params, getImage, get_params2D, set_params2D
	from alignment    	import fine_2D_refinement, crit2d
	from statistics   	import add_oe_series, ave_var_series, fsc_mask
	from filter 		import filt_from_fsc_bwt,filt_table
	from morphology         import ctf_2, ctf_1d
	import os
	import sys
	import types
	output = sys.stdout
	
	from utilities import print_begin_msg, print_end_msg, print_msg
	max_iter  = int(maxit)
	last_ring = int(ou)

	print_begin_msg("ali2d_e")	
	
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
			mask = getImage(maskfile)
		else:
			print_msg("Maskfile                    : user provided in-core mask\n")
			mask = maskfile
	else :
		print_msg("Maskfile                    : default, a circle with radius %i\n"%(last_ring))
		mask = model_circle(last_ring, nx, nx)
	# initialize data for the reference preparation function
	ref_data = []
	ref_data.append( mask )
	ref_data.append( center )
	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Search range                : %-5.2f\n"%(br))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("Error tolerance             : %f\n"%(eps))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("data with CTF               : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n\n"%(snr))
	
	
	# create the output directory, if it does not existm
	if os.path.exists(outdir):  
	    os.system('rm -rf '+outdir)
	os.mkdir(outdir)
	low = 0.5

	# read images
	del ima
	data = EMData.read_images(stack)
	nima = len(data)
	if(CTF):
		parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
		#                        0                1              2     3              4             5                6
		ctf_params = get_arb_params(data[0], parnames)
		data_had_ctf = ctf_params[6]
		ctm = ctf_1d(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
		lctf = len(ctm)
		ctf2 = []
		ctf2.append([0.0]*lctf)
		ctf2.append([0.0]*lctf)
		ctfb2 = [0.0]*lctf
		for im in xrange(nima):
			ctf_params = get_arb_params(data[im], parnames)
			ctm = ctf_2(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
			k = im%2
			for i in xrange(lctf):  ctf2[k][i] += ctm[i]
			if(ctf_params[6] == 0):
				st = Util.infomask(data[im], mask, False)
				data[im] -= st[0]
				from filter import filt_ctf
				data[im] = filt_ctf(data[im], ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
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
		dropImage(tavg, os.path.join(outdir, "aqe_%03d.hdf"%(Iter)))

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
		dropImage(tavg, os.path.join(outdir, "aqf_%03d.hdf"%(Iter)))

	if(CTF and data_had_ctf == 0):
		for im in xrange(nima): data[im].set_attr('ctf_applied', 0)
	from utilities import write_headers
	write_headers(stack, data, range(nima))
	print_end_msg("ali2d_e")	

def ali2d_m(stack, refim, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, xr = 0, yr = 0, ts = 1, center = 1, maxit = 10, CTF = False, snr = 1.0, rand_seed = 1000, MPI=False):
# 2D multi-reference alignment using rotational ccf in polar coords and quadratic interpolation
	if MPI:
		ali2d_m_MPI(stack, refim, outdir, maskfile, ir, ou, rs, xr, yr, ts, center, maxit, CTF, snr, rand_seed )
		return

	from utilities      import   model_circle, compose_transform2, combine_params2, dropImage, getImage
	from utilities	    import   center_2D, get_arb_params, get_im, get_params2D, set_params2D
	from statistics     import   fsc
	from alignment      import   Numrinit, ringwe, Applyws, fine_2D_refinement
	from fundamentals   import   rot_shift2D, fshift
	from morphology     import   ctf_2
	from filter         import   filt_btwl, filt_params
	from random         import   seed, randint
	import os
	import sys
	
	from utilities      import   print_begin_msg, print_end_msg, print_msg
	print_begin_msg("ali2d_m")

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Reference stack             : %s\n"%(refim))	
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Maskfile                    : %s\n"%(maskfile))
	print_msg("Inner radius                : %i\n"%(first_ring))

	ima  = EMData()
	ima.read_image(stack, 0)
	nx = ima.get_xsize()
	# default value for the last ring
	if (last_ring == -1): last_ring = nx//2-2

	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Ring step                   : %i\n"%(rstep))
	print_msg("X search range              : %i\n"%(xrng))
	print_msg("Y search range              : %i\n"%(yrng))
	print_msg("Translational step          : %i\n"%(step))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("data with CTF               : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	print_msg("Random seed                 : %i\n\n"%(rand_seed))

	# create the output directory, if it does not exist
	if os.path.exists(outdir):  os.system('rm -rf '+outdir)
	os.mkdir(outdir)
	output = sys.stdout

	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask=getImage(maskfile)
		else: mask = maskfile
	else : mask = model_circle(last_ring, nx, nx)
	#  references
	refi = []
	numref = EMUtil.get_image_count(refim)
	#  CTF stuff
	if(CTF):
		parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
		ctf_params = get_arb_params(ima, parnames)
		data_had_ctf = ctf_params[6]
		ctm = ctf_2(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
		lctf = len(ctm)
		ctf2 = [[[0.0]*lctf for k in xrange(2)] for j in xrange(numref)]
	# do the alignment
	# IMAGES ARE SQUARES!
	#  center is in SPIDER convention
	cnx = int(nx/2)+1
	cny = cnx

	mode = "F"
	#precalculate rings
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr = ringwe(numr, mode)
	# reference images
	again = True
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
	while(Iter < max_iter and again):
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
			if(CTF):
				for i in xrange(lctf): 
					ctf2[j][0][i] = 0.0
					ctf2[j][1][i] = 0.0
		assign = [[] for i in xrange(numref)]
		for im in xrange(nima):
			#
			if(CTF):
				ctf_params = get_arb_params(data[im], parnames)
				if(ctf_params[6] == 0):
					st = Util.infomask(data[im], mask, False)
					data[im] -= st[0]
					from filter import filt_ctf
					data[im] = filt_ctf(data[im], ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
					data[im].set_attr('ctf_applied', 1)
			#normalize
			alpha, sx, sy, mirror, scale =  get_params2D(data[im])
			# align current image to the reference
			data[im].process_inplace("normalize.mask", {"mask":mask, "no_sigma":0})
			[angt, sxst, syst, mirrort, xiref, peakt] = \
				Util.multiref_polar_ali_2d(data[im], ringref, xrng, yrng, step, mode, numr, cnx-sx, cny-sy)
			iref = int(xiref)
			# combine parameters and set them to the header, ignore previous angle and mirror
			[alphan, sxn, syn, mn] = combine_params2(0.0, sx, sy, 0, angt, sxst, syst, mirrort)
			#print  "combine with previous ",Iter,im,psin, sxn, syn, mn, iter
			set_params2D(data[im], [alphan, sxn, syn, mn, scale])
			data[im].set_attr('assign',iref)
			# apply current parameters and add to the average
			temp = rot_shift2D(data[im], alphan, sxn, syn)
			if mn: temp.process_inplace("mirror", {"axis":'x'})
			it = im%2
			Util.add_img( refi[iref][it], temp)
			if(CTF):
				ctm = ctf_2(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
				for i in xrange(lctf):  ctf2[iref][it][i] += ctm[i]
			assign[iref].append(im)
			refi[iref][2] += 1

		del ringref
		if(again):
			a1 = 0.0
			for j in xrange(numref):
				msg = "   group #%3d   number of particles = %7d\n"%(j,refi[j][2])
				print_msg(msg)
				if(refi[j][2] < 4):
					#ERROR("One of the references vanished","ali2d_m",1)
					#  if vanished, put a random image there
					assign[j] = []
					assign[j].append( randint(0, nima-1) )
					refi[j][0] = data[assign[j][0]].copy()
				else:
					max_inter=0  # switch off fine refi.
					br = 1.75
					#  the loop has to 
					for INter in xrange(max_inter+1):
						# Calculate averages at least ones, meaning even if no within group refinement was requested
						if(CTF):
							for i in xrange(lctf):  ctm[i] = 1.0 / (ctf2[j][0][i] + 1.0/snr)
							from filter import filt_table
							av1 = filt_table( refi[j][0], ctm)
							for i in xrange(lctf):  ctm[i] = 1.0 / (ctf2[j][1][i] + 1.0/snr)
							av2 = filt_table( refi[j][1], ctm)
							frsc = fsc(av1, av2, 1.0, os.path.join(outdir,"drm_%03d_%04d"%(Iter, j)))
							#Now the total average
							for i in xrange(lctf):  ctm[i] = 1.0 / (ctf2[j][0][i] + ctf2[j][1][i] + 1.0/snr)
							refi[j][0] = filt_table( Util.addn_img( refi[j][0], refi[j][1] ), ctm)
						else:
							frsc = fsc(refi[j][0], refi[j][1], 1.0, os.path.join(outdir,"drm_%03d_%04d"%(Iter, j)))
							Util.add_img( refi[j][0], refi[j][1] )
							Util.mul_scalar( refi[j][0], 1.0/float(refi[j][2]) )
						#  low pass filter references
						lowfq, highfq = filt_params(frsc, low = 0.1)
						refi[j][0]  = filt_btwl(refi[j][0], lowfq, highfq)
						refi[j][0], csx, csy = center_2D(refi[j][0], center)
						msg = "   group #%3d   filter parameters = %6.4f, %6.4f,  center parameters (x,y) = %10.3e, %10.3e\n"%(j, lowfq, highfq, csx,csy)
						print_msg(msg)
						for i in xrange(len(assign[j])):
							im = assign[j][i]
							alpha, sx, sy, mirror, scale =  get_params2D(data[im])
							alphan, sxn, syn, scale = compose_transform2(alpha, sx, sy, 1.0, 0.0, -csx, -csy, 1.0)
							set_params2D(data[im], [alphan, sxn, syn, mirror, scale])
						# refine images within the group
						#  Do the refinement only if max_inter>0, but skip it for the last iteration.
						if(INter < max_inter):
							fine_2D_refinement(data, br, mask, refi[j][0], j)
							#  Calculate updated average
							refi[j][0].to_zero()
							refi[j][1].to_zero()
							for i in xrange(len(assign[j])):
								im = assign[j][i]
								alpha, sx, sy, mirror, scale =  get_params2D(data[im])
								# apply current parameters and add to the average
								temp = rot_shift2D(data[im], alpha, sx, sy)
								if mn: temp.process_inplace("mirror", {"axis":'x'})
								it = im%2
								Util.add_img( refi[j][it], temp)

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
			if (a1 > a0 and Iter < max_iter) :  a0 = a1
			else:
				os.system('cp '+newrefim+' multi_ref.hdf')
				break

	if(CTF):
		if(data_had_ctf == 0):
			for im in xrange(nima): data[im].set_attr('ctf_applied', 0)
	from utilities import write_headers
	write_headers(stack, data, range(nima))
	print_end_msg("ali2d_m")

def ali2d_m_MPI(stack, refim, outdir, maskfile = None, ir=1, ou=-1, rs=1, xr=0, yr=0, ts=1, center=1, maxit=10, CTF = False, snr = 1., rand_seed = 1000):
# 2D multi-reference alignment using rotational ccf in polar coords and quadratic interpolation

	from utilities      import   model_circle, compose_transform2, combine_params2, dropImage, getImage
	from utilities      import   reduce_EMData_to_root, bcast_EMData_to_all, bcast_number_to_all
	from utilities      import   send_attr_dict, recv_attr_dict
	from utilities	    import   center_2D
	from statistics     import   fsc_mask
	from alignment      import   Numrinit, ringwe, Applyws
	from fundamentals   import   rot_shift2D, fshift
	from utilities      import   get_arb_params, get_params2D, set_params2D
	from random         import   seed, randint
	from morphology     import   ctf_2
	from filter         import   filt_btwl, filt_params
	from Numeric        import   reshape, shape
	from utilities      import   print_msg, print_begin_msg, print_end_msg
	import os
	import sys
	from mpi 	  import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	  import mpi_reduce, mpi_bcast, mpi_barrier, mpi_recv, mpi_send
	from mpi 	  import MPI_SUM, MPI_FLOAT, MPI_INT, MPI_TAG_UB

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	#  chose a random node as a main one...
	main_node = 0
	if(myid == 0): main_node = randint(0,number_of_proc-1)
	main_node = mpi_bcast(main_node, 1, MPI_INT, 0, MPI_COMM_WORLD)
	main_node = main_node[0][0]
	# create the output directory, if it does not exist
	if(myid == main_node):
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)

	if (myid == main_node):
		print_begin_msg("ali2d_m_MPI")

	nima = EMUtil.get_image_count(stack)
	
	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)

	nima = EMUtil.get_image_count(stack)
	ima = EMData()
	ima.read_image(stack, image_start)
	# 

	if (myid == main_node):
		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Reference stack             : %s\n"%(refim))	
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Maskfile                    : %s\n"%(maskfile))
		print_msg("Inner radius                : %i\n"%(first_ring))

	nx=ima.get_xsize()
	# default value for the last ring
	if(last_ring == -1): last_ring=nx//2-2
	
	if (myid == main_node):
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Ring step                   : %i\n"%(rstep))
		print_msg("X search range              : %i\n"%(xrng))
		print_msg("Y search range              : %i\n"%(yrng))
		print_msg("Translational step          : %i\n"%(step))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("data with CTF               : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Random seed                 : %i\n\n"%(rand_seed))	

	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask=getImage(maskfile)
		else: mask = maskfile
	else : mask = model_circle(last_ring, nx, nx)
	#  references, do them on all processors...
	refi = []
	numref = EMUtil.get_image_count(refim)
	#  CTF stuff
	if(CTF):
		parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
		ctf_params = get_arb_params(ima, parnames)
		data_had_ctf = ctf_params[6]
		ctm = ctf_2(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
		lctf = len(ctm)
	# do the alignment
	# IMAGES ARE SQUARES!
	#  center is in SPIDER convention
	cnx = int(nx/2)+1
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
		temp = EMData()
		temp.read_image(refim, j)
		#  even, odd, numer of even, number of images.  After frc, totav
		refi.append([temp, ima.copy(), 0])
	#  for each node read its share of data
	data = EMData.read_images(stack, range(image_start, image_end))
	for im in xrange(image_start, image_end):
		data[im].set_attr('ID', im)
		if(CTF):
			ctf_params = get_arb_params(data[im], parnames)
			if(ctf_params[6] == 0):
				st = Util.infomask(data[im], mask, False)
				data[im] -= st[0]
				from filter import filt_ctf
				data[im] = filt_ctf(data[im], ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
				data[im].set_attr('ctf_applied', 1)
	if(myid == main_node):  seed(rand_seed)
	a0 = -1.
	again = True
	Iter = 0
	while(Iter < max_iter and again):
		ringref = []
		for j in xrange(numref):
			refi[j][0].process_inplace("normalize.mask", {"mask":mask, "no_sigma":1}) # normalize reference images to N(0,1)
			cimage =Util.Polar2Dm(refi[j][0] , cnx, cny, numr, mode)
			Util.Frngs(cimage, numr)
			Applyws(cimage, numr, wr)
			ringref.append(cimage)
			# zero refi
			refi[j][0].to_zero()
			refi[j][1].to_zero()
			refi[j][2] = 0
		if(CTF): ctf2 = [[[0.0]*lctf for k in xrange(2)] for j in xrange(numref)]
		#from utilities import ttime
		#print  ttime()
		assign = [[] for i in xrange(numref)]
		# begin MPI section
		for im in xrange(image_start, image_end):
			#
			alpha, sx, sy, mirror, scale =  get_params2D(data[im-image_start])
			#normalize
			# align current image to the reference
			data[im-image_start].process_inplace("normalize.mask", {"mask":mask, "no_sigma":0}) # subtract average under the mask
			[angt, sxst, syst, mirrort, xiref, peakt] = \
				 Util.multiref_polar_ali_2d(data[im-image_start], ringref, xrng, yrng, step, mode, numr, cnx-sx, cny-sy)
			iref = int(xiref)
			# combine parameters and set them to the header, ignore previous angle and mirror
			[alphan, sxn, syn, mn] = combine_params2(0.0, sx, sy, 0, angt, sxst, syst, mirrort)
			#if(iref  ==0):print  "combine with previous ",Iter,im,alphan, sxn, syn, mn, iref
			set_params2D(data[im-image_start], [alphan, sxn, syn, mn, scale])
			data[im-image_start].set_attr('assign',iref)
			# apply current parameters and add to the average
			temp = rot_shift2D(data[im-image_start], alphan, sxn, syn)
			if mn: temp.process_inplace("mirror", {"axis":'x'})
			it = im%2
			Util.add_img( refi[iref][it], temp)
			assign[iref].append(im)
			if(CTF):
				#  I wonder whether params are still there....
				ctf_params = get_arb_params(data[im-image_start], parnames)
				ctm = ctf_2(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
				for i in xrange(lctf):  ctf2[iref][it][i] += ctm[i]
			#assign[im] = iref
			refi[iref][2] += 1.0
			#if(im%50 == 0):
			#	output.write( "\n" )
			#	output.write( " %6d " % im )
			#	output.flush()
			#output.write( "." )
			#output.write( str (refi[iref][2]) )
			#output.flush( )
		#output.write( "\n" )
		#print  ttime()
		del ringref
		# end MPI section, bring partial things together, calculate new reference images, broadcast them back
		if(CTF):
		# bring ctf2 together on main node
			s = shape(ctf2)
			ctf2  = mpi_reduce(ctf2, 2*lctf*numref, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
			if(myid == main_node): ctf2 = reshape(ctf2, s)
		for j in xrange(numref):
			reduce_EMData_to_root(refi[j][0], myid, main_node)
			reduce_EMData_to_root(refi[j][1], myid, main_node)
			refi[j][2] = mpi_reduce(refi[j][2], 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
			if(myid == main_node): refi[j][2] = int(refi[j][2][0])
		# gather assignements
		for j in xrange(numref):
			if(myid == main_node):
				for n in xrange(number_of_proc):
					if(n != main_node):
						ln =  mpi_recv(1, MPI_INT, n, MPI_TAG_UB, MPI_COMM_WORLD)
						lis = mpi_recv(ln[0], MPI_INT, n, MPI_TAG_UB, MPI_COMM_WORLD)
						for l in xrange(ln[0]): assign[j].append(int(lis[l]))
			else:
				mpi_send(len(assign[j]), 1,MPI_INT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)
				mpi_send(assign[j], len(assign[j]), MPI_INT, main_node, MPI_TAG_UB, MPI_COMM_WORLD)

		if(myid == main_node):
			# replace the name of the stack with reference with the current one
			refim = os.path.join(outdir,"aqm%03d.hdf"%Iter)
			a1 = 0.0
			for j in xrange(numref):
				if(refi[j][2]<4):
					#ERROR("One of the references vanished","ali2d_m_MPI",1)
					#  if vanished, put a random image (only from main node!) there
					assign[j] = []
					assign[j].append( randint(image_start, image_end-1) - image_start )
					refi[j][0] = data[assign[j][0]].copy()
				else:
					if(CTF):
						for i in xrange(lctf):  ctm[i] = 1.0 / (ctf2[j][0][i] + 1.0/snr)
						from filter import filt_table
						av1 = filt_table( refi[j][0], ctm)
						for i in xrange(lctf):  ctm[i] = 1.0 / (ctf2[j][1][i] + 1.0/snr)
						av2 = filt_table( refi[j][1], ctm)
						frsc = fsc_mask(av1, av2, mask, 1.0, os.path.join(outdir,"drm%03d%04d"%(Iter, j)))
						#Now the total average
						for i in xrange(lctf):  ctm[i] = 1.0 / (ctf2[j][0][i] + ctf2[j][1][i] + 1.0/snr)
						refi[j][0] = filt_table( Util.addn_img( refi[j][0], refi[j][1] ), ctm)
					else:
						frsc = fsc_mask(refi[j][0], refi[j][1], mask, 1.0, os.path.join(outdir,"drm%03d%04d"%(Iter, j)))
						Util.add_img( refi[j][0], refi[j][1] )
						Util.mul_scalar( refi[j][0], 1.0/float(refi[j][2]) )
					#  low pass filter references
					lowfq, highfq = filt_params(frsc, low = 0.3)
					lowfq = max(0.05, lowfq)
					highfq = min(0.4,max(highfq,lowfq)+0.1)
					#from filter import filt_from_fsc, filt_table
					#filt = filt_from_fsc(frsc, 0.1)
					#refi[j][0] = filt_table(refi[j][0], filt)
					#del filt
					#from morphology import threshold
					#refi[j][0]  = filt_btwl(threshold(refi[j][0]), lowfq, highfq)
					refi[j][0]  = filt_btwl(refi[j][0], lowfq, highfq)
					#refi[j][0], csx, csy = center_2D(refi[j][0], center, searching_range = max(xrng, yrng))
					refi[j][0], csx, csy = center_2D(refi[j][0], center)
					msg = "   group #%3d   filter parameters = %6.4f, %6.4f,  center parameters (x,y) = %10.3e, %10.3e\n"%(j, lowfq, highfq, csx, csy)
					print_msg(msg)
				# write the current average
				TMP = []
				for i_tmp in xrange(len(assign[j])):
					TMP.append(float(assign[j][i_tmp]))
				TMP.sort()
				refi[j][0].set_attr_dict({'ave_n': refi[j][2],  'members': TMP })
				del TMP
				#  high pass filter references  (activate two lines below)
				#from filter import filt_btwo
				#refi[j][0] = filt_btwo(refi[j][0], 0.01, 0.1, 0.1)
				refi[j][0].process_inplace("normalize.mask", {"mask":mask, "no_sigma":1})
				refi[j][0].write_image(refim, j)
				a1 += refi[j][0].cmp("dot", refi[j][0], {"negative":0, "mask":mask})
				# send refi[j][0]  back!
			Iter += 1
			msg = "ITERATION #%3d        criterion = %20.7e\n\n"%(Iter,a1)
			print_msg(msg)
			#if(a1 > a0):  a0 = a1
			#else:         again = False
		#again = mpi_bcast(again, 1, MPI_INT, main_node, MPI_COMM_WORLD)
		Iter  = bcast_number_to_all( Iter, main_node )
		if(CTF):  del  ctf2
		if(again):
			for j in xrange(numref):
				bcast_EMData_to_all(refi[j][0], myid, main_node)
	#  clean up
	del refi
	del assign
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
	if(CTF and data_had_ctf == 0):
		for im in xrange(len(data)): data[im].set_attr('ctf_applied', 0)
	list_params = ['xform.align2d', 'assign']
	if(myid == main_node): recv_attr_dict(main_node, stack, data, list_params, image_start, image_end, number_of_proc)
	else: send_attr_dict(main_node, data, list_params, image_start, image_end)

	if (myid==main_node):	print_end_msg("ali2d_m_MPI")

def ali2d_ra(stack, maskfile = None, ir = 1, ou = -1, rs = 1, maxit = 10, check_mirror = False, CTF = False, rand_seed = 1000):
# 2D rotational alignment using ccf in polar coordinates

	from utilities    import model_circle, compose_transform2, combine_params2, dropImage, get_im, get_arb_params
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
	print_msg("data with CTF               : %s\n"%(CTF))
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
			ctf_params = get_arb_params(temp, parnames)
			if(im == 0):  data_had_ctf = temp.get_attr('ctf_applied')
			ctf = ctf_2(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
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
			if(ctf_params[6] == 0):
				from filter import filt_ctf
				temp = filt_ctf(temp, ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
				temp.set_attr('ctf_applied', 1)
			from filter       import filt_table
			refc = filt_table(temp, ctf2)
			
			alpha_original = temp.get_attr('alpha')
			sx =  temp.get_attr('sx')
			sy =  temp.get_attr('sy')
			miri = temp.get_attr('mirror')
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
			alpha_original = temp.get_attr('alpha')
			sx =  temp.get_attr('sx')
			sy =  temp.get_attr('sy')
			miri = temp.get_attr('mirror')
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
		temp.set_attr_dict({'alpha':alphan, 'sx':sxn, 'sy':syn, 'mirror': mir})
		write_header(stack, temp, im)
		#temp.write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)
	print_end_msg("ali2d_ra")

def ali2d_rag(stack, maskfile = None, ir = 1, ou = -1, rs = 1, maxit = 10, check_mirror = False, CTF = False, rand_seed = 1000):
# 2D rotational alignment using ccf in polar coordinates and gridding-based interpolation

	from utilities    import model_circle, compose_transform2, combine_params2, dropImage, get_im, get_arb_params
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
	print_msg("data with CTF               : %s\n"%(CTF))
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
			ctf_params = get_arb_params(temp, parnames)
			if(im == 0):  data_had_ctf = temp.get_attr('ctf_applied')
			ctf = ctf_2(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
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
			if(ctf_params[6] == 0):
				from filter import filt_ctf
				temp = filt_ctf(temp, ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
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
		DB=EMAN2DB()
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

	from utilities    import model_circle, combine_params2, dropImage
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
		if(type(maskfile) is types.StringType):  mask2D = getImage(maskfile)
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
		DB=EMAN2DB()
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
	from utilities 		import model_circle, combine_params2, dropImage
	from utilities      import get_input_from_string, getImage, get_arb_params, set_arb_params
	from fundamentals 	import rot_shift2D
	from statistics 	import add_oe_series, ave_series_ctf, ave_series, fsc_mask
	from alignment 		import Numrinit, ringwe, ali2d_s, align2d
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
			mask=getImage(maskfile)
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
		parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
		#                        0                1              2     3              4                5                6
		ctf_params = get_arb_params(ima, parnames)
		data_had_ctf = ctf_params[6]
		ctm = ctf_2(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
		lctf = len(ctm)
		ctf2 = [[[0.0]*lctf for j in xrange(2)] for i in xrange(NG)]
		ctfb2 = [[0.0]*lctf for i in xrange(NG)]
	del ima
	all_data = EMData.read_images(stack)
	for im in xrange(nima):
		all_data[im].set_attr('ID', im)
		k = im%NG
		if(CTF):
			ctf_params = get_arb_params(all_data[im], parnames)
			ctm = ctf_2(nx, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
			kl = (im//2)%NG  # not sure it will work for NG>2
			for i in xrange(lctf):
				ctf2[k][kl][i] += ctm[i]
			if(ctf_params[6] == 0):
				st = Util.infomask(all_data[im], mask, False)
				all_data[im] -= st[0]
				#      filt_ctf(all_data[im]ge,	     dz,		  cs,		   voltage,	     pixel_size,     amp_contrast=0.1,	  b_factor=0.0):
				all_data[im] = filt_ctf(all_data[im], ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
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
		dropImage(tavg[k],os.path.join(outdir, "aqc_%03d_%03d.hdf"%(k, 0)))
	fscross = fsc_mask(tavg[0], tavg[1], mask, 1.0, os.path.join(outdir, "drcross_%03d"%(0)))

	if(CTF): total_ave = ave_series_ctf(all_data, ctf_tot)
	else:    total_ave = ave_series(all_data)
	dropImage(total_ave, os.path.join(outdir, "total_ave_%03d.hdf"%(0)))
	a0 = total_ave.cmp("dot", total_ave, dict(negative = 0, mask = mask))
	msg = "Initial criterion = %-20.7e\n"%(a0)
	print_msg(msg)
	params = ["alpha", "sx", "sy", "mirror"]
	# initialize data for the reference preparation function
	#  mask can be modified in user_function
	ref_data = []
	ref_data.append( mask )
	ref_data.append( center )
	ref_data.append( None )
	ref_data.append( None )
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
				ali2d_s(data[k], numr, wr, cs[k], tavg[k], cnx, cny, xrng[N_step], yrng[N_step], step[N_step], mode, range(len(data[k])))
				av1, av2 = add_oe_series(data[k])
				if(CTF):
					tavg[k] = filt_table(Util.addn_img(av1, av2), ctfb2[k])
					av1    = filt_table(av1, ctf2[k][0])
					av2    = filt_table(av2, ctf2[k][1])
				else:
					tavg[k] = (av1+av2)/len(data[k])
				dropImage(tavg[k], os.path.join(outdir, "aqc_%03d_%03d.hdf"%(k, total_iter)))

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
			favg2 = rot_shift2D(tavg[0], alpha, sx, sy, interpolation_method="gridding")
			if  mirror: favg2.process_inplace("mirror",{"axis":'x'})
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
				dropImage(tavg[k], os.path.join(outdir, "aqf_%03d_%03d.hdf"%(k, total_iter)))

			if(CTF): total_ave = ave_series_ctf(all_data, ctf_tot)
			else:    total_ave = ave_series(all_data)
			dropImage(total_ave, os.path.join(outdir, "total_ave_%03d.hdf"%(total_iter)))
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
		DB=EMAN2DB()
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

	from utilities      import model_circle, dropImage
	from utilities      import getImage, get_input_from_string
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

	if CTF :
		from reconstruction import recons3d_4nn_ctf
		ctf_dicts = ["defocus", "Cs", "voltage", "Pixel_size", "amp_contrast", "B_factor", "ctf_applied" ]
		#                     0          1         2             3                 4                   5               6
	else   : from reconstruction import recons3d_4nn

	if os.path.exists(outdir):  os.system('rm -rf '+outdir)
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
	print_msg("data with CTF               : %s\n"%(CTF))
	print_msg("Reference projection method : %s\n"%(ref_a))
	print_msg("Symmetry group              : %s\n\n"%(sym))

	if (maskfile) :
		if (type(maskfile) is types.StringType): mask3D = getImage(maskfile)
		else                                  : mask3D = maskfile
	else          :   mask3D = model_circle(last_ring, nx, nx, nx)
	mask = model_circle(last_ring, nx, nx)
	#dropImage(vol, os.path.join(outdir,"ref_vol00.hdf"))

	del ima
	data = EMData.read_images(stack)
	nima = len(data)
	for im in xrange(nima):
		data[im].set_attr('ID', im)
		if(CTF):
			if(data[im].get_attr_default('ctf_applied', 2) > 0):
				ERROR("data cannot be ctf-applied","ali2d_ar",1)
			st = Util.infomask(data[im], mask, False)
			data[im] -= st[0]

	# initialize data for the reference preparation function
	ref_data = []
	ref_data.append( mask3D )
	ref_data.append( center )
	ref_data.append( None )
	ref_data.append( None )

	# do the projection matching
	for N_step in xrange(lstp):
 		for Iter in xrange(max_iter):
			print_msg("ITERATION #%3d\n"%(N_step*max_iter + Iter+1))

			if(an[N_step] == -1):	proj_ali_incore(vol, mask3D, data, first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], ref_a, sym, CTF, MPI=False)
			else:	                proj_ali_incore_local(vol, mask3D, data, first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], an[N_step], ref_a, sym, CTF, MPI=False)
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
			dropImage(vol, os.path.join(outdir, "vol%04d.hdf"%(N_step*max_iter+Iter+1)))
			did, vav = recons3d_nn_SSNR(data,  mask2D = None, ring_width=1, npad =1, sign=1, symmetry = sym, CTF =CTF)
			vol = fft(Util.divn_filter(fft(vol), vav))
			ref_data[2] = vol
			ref_data[3] = fscc
			#  call user-supplied function to prepare reference image, i.e., center and filter it
			vol, cs = user_func( ref_data )
			if center == 1:
				from utilities import rotate_3D_shift
				rotate_3D_shift(data, cs)
			dropImage(vol, os.path.join(outdir, "volf%04d.hdf"%(N_step*max_iter+Iter+1)))
	#  here we  write header info
	from utilities import write_headers
	write_headers( stack, data, range(nima))
	print_end_msg("ali3d_a")

def ali3d_d(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta="10 6 4 4", an="-1", 
	    center = 1.0, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym="c1", user_func_name="ref_ali3d", pinfo = False, MPI=False):
	if MPI:
		ali3d_d_MPI(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, yr, ts, delta, an, center, maxit, CTF, snr, ref_a, sym, user_func_name)
		return

	from utilities      import model_circle, dropImage
	from utilities      import getImage, get_input_from_string
	from utilities      import get_arb_params, set_arb_params
	from filter         import filt_params, filt_btwl, filt_from_fsc, filt_table, fit_tanh, filt_tanl
	from alignment	    import proj_ali_incore, proj_ali_incore_local
	from statistics     import fsc_mask
	import os
	import types
	from utilities      import print_begin_msg, print_end_msg, print_msg
	print_begin_msg("ali3d_d")

	import user_functions
	user_func = user_functions.factory[user_func_name]

	if CTF :
		from reconstruction import recons3d_4nn_ctf
		ctf_dicts = ["defocus", "Cs", "voltage", "Pixel_size", "amp_contrast", "B_factor", "ctf_applied" ]
		#                     0          1         2             3                 4                   5               6
	else   : from reconstruction import recons3d_4nn

	if os.path.exists(outdir):  os.system('rm -rf '+outdir)
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
	print_msg("data with CTF               : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	print_msg("Reference projection method : %s\n"%(ref_a))
	print_msg("Symmetry group              : %s\n\n"%(sym))

	if (maskfile) :
		if (type(maskfile) is types.StringType): mask3D = getImage(maskfile)
		else                                  : mask3D = maskfile
	else          :   mask3D = model_circle(last_ring, nx, nx, nx)
	mask = model_circle(last_ring, nx, nx)
	#dropImage(vol, os.path.join(outdir,"ref_vol00.hdf"))
	if pinfo:  outf = file(os.path.join(outdir, "progress"), "w")
	else:     outf = None

	data = EMData.read_images(stack)
	nima = len(data)
	for im in xrange(nima):
		data[im].set_attr('ID', im)
		if(CTF):
			ctf_params = get_arb_params(data[im], ctf_dicts)
			if(im == 0): data_had_ctf = ctf_params[6]
			if(ctf_params[6] == 0):
				st = Util.infomask(data[im], mask, False)
				data[im] -= st[0]
				from filter import filt_ctf
				data[im] = filt_ctf(data[im], ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
				data[im].set_attr('ctf_applied', 1)

	# initialize data for the reference preparation function
	ref_data = []
	ref_data.append( mask3D )
	ref_data.append( center )
	ref_data.append( None )
	ref_data.append( None )

	# do the projection matching
	for N_step in xrange(lstp):
 		for Iter in xrange(max_iter):
			print_msg("ITERATION #%3d\n"%(N_step*max_iter + Iter+1))
			#from filter import filt_gaussinv
			#vol = filt_gaussinv(vol, 0.175, True)
			#dropImage(vol, os.path.join(outdir, replace("vhl%4d.spi"%(N_step*max_iter+Iter+1),' ','0')), "s")
			if(an[N_step] == -1):	proj_ali_incore(vol, mask3D, data, first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], ref_a, sym, finfo = outf, MPI=False)
			else:	               proj_ali_incore_local(vol, mask3D, data, first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], an[N_step], ref_a, sym, finfo = outf, MPI=False)
			#  3D stuff
			if(CTF): vol1 = recons3d_4nn_ctf(data, range(0,nima,2), snr, 1, sym)
			else:	 vol1 = recons3d_4nn(data, range(0,nima,2), sym)
			if(CTF): vol2 = recons3d_4nn_ctf(data, range(1,nima,2), snr, 1, sym)
			else:	 vol2 = recons3d_4nn(data, range(1,nima,2), sym)

			fscc = fsc_mask(vol1, vol2, mask3D, 1.0, os.path.join(outdir, "resolution%04d"%(N_step*max_iter+Iter+1)))
			del vol1
			del vol2

			# calculate new and improved 3D
			if(CTF): vol = recons3d_4nn_ctf(data, range(nima), snr, 1, sym)
			else:	 vol = recons3d_4nn(data, range(nima), sym)
			# store the reference volume
			dropImage(vol, os.path.join(outdir, "vol%04d.hdf"%(N_step*max_iter+Iter+1)))
			ref_data[2] = vol
			ref_data[3] = fscc
			#  call user-supplied function to prepare reference image, i.e., center and filter it
			vol, cs = user_func( ref_data )

			if center == 1:
				from utilities import rotate_3D_shift
				rotate_3D_shift(data, cs)
			dropImage(vol, os.path.join(outdir, "volf%04d.hdf"%(N_step*max_iter+Iter+1)))
	#  here we  write header info
	if(CTF and data_had_ctf == 0):
		for im in xrange(nima): data[im].set_attr('ctf_applied', 0)
	from utilities import write_headers
	write_headers( stack, data, range(nima))
	print_end_msg("ali3d_d")

def ali3d_d_MPI(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an="-1",
	    center = 1.0, maxit = 5, CTF = False, snr = 1.0,  ref_a="S", sym="c1", user_func_name="ref_ali3d",debug=False):

	from utilities      import model_circle, reduce_EMData_to_root, bcast_EMData_to_all, dropImage
	from utilities      import bcast_list_to_all, bcast_string_to_all, getImage, get_input_from_string
	from utilities      import get_arb_params, set_arb_params, dropSpiderDoc,recv_attr_dict, send_attr_dict
	from utilities      import dropSpiderDoc, get_im
	from alignment	    import proj_ali_incore
	from random	    import randint
	from fundamentals   import rot_avg_image
	import os
	import types
	from reconstruction import rec3D_MPI, rec3D_MPI_noCTF
	from utilities      import print_begin_msg, print_end_msg, print_msg
	from mpi 	    import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	    import mpi_barrier

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	if(myid == main_node):
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)
	from string import replace
	if debug:
        	info_file = outdir+("/progress%04d"%myid)
        	finfo = open(info_file, 'w')
	else:
		finfo = None

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
		from  alignment	    import proj_ali_incore_local
		an      = get_input_from_string(an)

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	if (myid == main_node):
		import user_functions
		user_func = user_functions.factory[user_func_name]

		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Reference volume            : %s\n"%(ref_vol))	
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Maskfile                    : %s\n"%(maskfile))
		print_msg("Inner radius                : %i\n"%(first_ring))

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	if last_ring < 0:	last_ring = int(nx/2) - 2

	if (myid == main_node):
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Ring step                   : %i\n"%(rstep))
		print_msg("X search range              : %s\n"%(xrng))
		print_msg("Y search range              : %s\n"%(yrng))
		print_msg("Translational step          : %s\n"%(step))
		print_msg("Angular step                : %s\n"%(delta))
		print_msg("Angular search range        : %s\n"%(an))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("data with CTF               : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Reference projection method : %s\n"%(ref_a))
		print_msg("Symmetry group              : %s\n\n"%(sym))

	if(maskfile):
		if(type(maskfile) is types.StringType): mask3D = getImage(maskfile)
		else:                                  mask3D = maskfile
	else:         mask3D = model_circle(last_ring, nx, nx, nx)
	nima            = EMUtil.get_image_count(stack)
	mask = model_circle(last_ring, nx, nx)

	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
	if debug:
		finfo.write( "image_start, image_end: %d %d\n" %(image_start, image_end) )
		finfo.flush()

	if(CTF):
		ctf_dicts = ["defocus", "Cs", "voltage", "Pixel_size", "amp_contrast", "B_factor", "ctf_applied" ]
		#                 0      1         2            3               4                  5               6
		from reconstruction import rec3D_MPI
	else:
		from reconstruction import rec3D_MPI_noCTF

	data = EMData.read_images(stack, range(image_start, image_end))
	for im in xrange(image_start, image_end):
		data[im].set_attr('ID', im)
		if(CTF):
			ctf_params = get_arb_params(data[im], ctf_dicts)
			if(im == image_start): data_had_ctf = ctf_params[6]
			if(ctf_params[6] == 0):
				st = Util.infomask(data[im], mask, False)
				data[im] -= st[0]
				from filter import filt_ctf
				data[im] = filt_ctf(data[im], ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
				data[im].set_attr('ctf_applied', 1)
		if(im%100==0 and debug) :
			finfo.write( '%d loaded  ' % im )
			finfo.flush()

	if (myid == main_node):
		# initialize data for the reference preparation function
		ref_data = []
		ref_data.append( mask3D )
		ref_data.append( center )
		ref_data.append( None )
		ref_data.append( None )

	# do the projection matching
	for N_step in xrange(lstp):
 		for Iter in xrange(max_iter):
			if(an[N_step] == -1): proj_ali_incore(vol, mask3D, data, first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], ref_a, sym, finfo = finfo, MPI=True)
			else:                proj_ali_incore_local(vol, mask3D, data, first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], an[N_step], ref_a, sym, finfo = finfo, MPI=True)

			if(CTF): vol, fscc = rec3D_MPI(data, snr, sym, mask3D, os.path.join(outdir, "resolution%04d"%(N_step*max_iter+Iter+1)), myid, main_node)
			else:    vol, fscc = rec3D_MPI_noCTF(data, sym, mask3D, os.path.join(outdir, "resolution%04d"%(N_step*max_iter+Iter+1)), myid, main_node)
	
			mpi_barrier(MPI_COMM_WORLD)
			if(myid == main_node):
				dropImage(vol,  os.path.join(outdir, "vol%04d.hdf"%(N_step*max_iter+Iter+1)))
				ref_data[2] = vol
				ref_data[3] = fscc
				#  call user-supplied function to prepare reference image, i.e., center and filter it

				vol, cs = user_func( ref_data )
				if center == 1:
					from utilities import rotate_3D_shift
					rotate_3D_shift(data, cs)
				dropImage(vol,  os.path.join(outdir, "volf%04d.hdf"%(N_step*max_iter+Iter+1)))
			bcast_EMData_to_all(vol, myid, main_node)
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
	if(CTF and data_had_ctf == 0):
		for im in xrange(len(data)): data[im].set_attr('ctf_applied', 0)
	par_str = ['phi', 'theta', 'psi', 's2x', 's2y']
	if(myid == main_node): recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	else: send_attr_dict(main_node, data, par_str, image_start, image_end)
	if (myid == main_node): print_end_msg("ali3d_d_MPI")

def ali3d_dB(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an="-1",
	    center = 1, maxit = 5, CTF = False, snr = 1.0,  ref_a="S", sym="c1"):
	#  THIS IS FOR PROCESSING BERLIN DATASET
	from utilities      import model_circle, reduce_EMData_to_root, bcast_EMData_to_all, dropImage
	from utilities      import bcast_list_to_all, bcast_string_to_all, getImage, get_input_from_string
	from utilities      import get_arb_params, set_arb_params, dropSpiderDoc,recv_attr_dict, send_attr_dict
	from utilities      import dropSpiderDoc,get_im
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
		if(type(maskfile) is types.StringType): mask3D = getImage(maskfile)
		else:                                  mask3D = maskfile
	else:         mask3D = model_circle(last_ring, nx, nx, nx)
	nima            = EMUtil.get_image_count(stack)
	mask = model_circle(last_ring, nx, nx)
	
	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)

	if(myid == main_node)       :  dropImage(vol, os.path.join(outdir,"ref_volf00.spi"), "s")
	if(CTF):
		ctf_dicts = ["defocus", "Cs", "voltage", "Pixel_size", "amp_contrast", "B_factor", "ctf_applied" ]
		#                 0      1         2            3               4                  5               6
		from reconstruction import rec3D_MPI
	else:
		from reconstruction import rec3D_MPI_noCTF
	from utilities import readSpiderDoc, set_arb_params
	prm = readSpiderDoc("params_new.doc")
	prm_dict = ["phi", "theta", "psi", "s2x", "s2y"]
	for im in xrange(image_start, image_end):
		data[im].set_attr('ID', im)
		set_arb_params(ima, prm[im], prm_dict)
		if(CTF):
			ctf_params = get_arb_params(data[im], ctf_dicts)
			if(im == image_start): data_had_ctf = ctf_params[6]
			if(ctf_params[6] == 0):
				st = Util.infomask(data[im], mask, False)
				data[im] -= st[0]
				from filter import filt_ctf
				data[im] = filt_ctf(data[im], ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
				data[im].set_attr('ctf_applied', 1)
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
			else:    recons3d_nn_SSNR_MPI(myid, data, model_circle(last_ring, nx, nx), ctf = CTF)

			if(CTF): volo, fscc = rec3D_MPI(data, snr, sym, mask3D, os.path.join(outdir, replace("resolution%4d"%(N_step*max_iter+Iter+1),' ','0')), myid, main_node)
			else:    volo, fscc = rec3D_MPI_noCTF(data, sym, mask3D, os.path.join(outdir, replace("resolution%4d"%(N_step*max_iter+Iter+1),' ','0')), myid, main_node)

			mpi_barrier(MPI_COMM_WORLD)
			if(myid == main_node):
				dropImage(volo,  os.path.join(outdir, replace("vol%4d.spi"%(N_step*max_iter+Iter+1),' ','0')), "s")
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
				dropImage(vol_ssnr,  os.path.join(outdir, replace("filter%4d.spi"%(N_step*max_iter+Iter+1),' ','0')), "s")
				#  Filter by volume
				vol = volo.filter_by_image(vol_ssnr)
				#vol  = filt_table(volo, filt)
				if(center == 1):
					cs   = vol.phase_cog()
					vol  = fshift(vol, -cs[0], -cs[1] -cs[2])
					volo = fshift(volo, -cs[0], -cs[1] -cs[2])
				dropImage(vol, os.path.join(outdir, replace("volf%4d.spi"%(N_step*max_iter+Iter+1),' ','0')), "s")
				adam = adaptive_mask(vol, 3000, 2.22)
				vol = threshold( adam*vol )
				h = histogram( vol )
				vol = threshold( adam*volo ) * get_im("rotav_tteftu_with_tRNA.spi") / threshold_to_minval( rot_avg_image(vol), h[len(h)//2])
				del volo
				#  Filter by volume
				vol = vol.filter_by_image(vol_ssnr)
				#vol  = filt_table(vol, filt)
				#vol = filt_btwl(vol, fl, fh)
				dropImage(vol, os.path.join(outdir, replace("vhlf%4d.spi"%(N_step*max_iter+Iter+1),' ','0')), "s")
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
			dropSpiderDoc(os.path.join(outdir, replace("new_params%8d"%((N_step*max_iter+Iter+1)*1000+myid),' ','0')), prm," phi, theta, psi, s2x, s2y, image number")
			del prm, prml

	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
	if(CTF and data_had_ctf == 0):
		for im in xrange(len(data)): data[im].set_attr('ctf_applied', 0)
	#par_str = ['phi', 'theta', 'psi', 's2x', 's2y']
	if(myid == main_node): recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	else: send_attr_dict(main_node, data, par_str, image_start, image_end)

def ali3d_mN(stack, ref_vol, outdir, maskfile = None, ir=1, ou=-1, rs=1, 
            xr  ="4 2  2  1",      yr="-1",
            ts  ="1 1 0.5 0.25",   delta="10  6  4  4",
	    center = 0, maxit= 5, CTF = False, snr = 1.0,  ref_a="S", symmetry="c1", MPI=False):
	if MPI:
		ali3d_m_MPI(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, yr, ts, 
		            delta, center, maxit, CTF, snr, ref_a, symmetry)
		return
	from utilities      import model_circle, reduce_EMData_to_root, bcast_EMData_to_all, dropImage
	from utilities      import bcast_list_to_all, bcast_string_to_all, getImage, get_input_from_string
	from utilities      import get_arb_params, set_arb_params, dropSpiderDoc
	from filter	    import filt_params,filt_btwl, filt_from_fsc, filt_table
	from alignment	    import proj_ali_incore_index
	from fundamentals   import fshift
	from utilities      import print_begin_msg, print_end_msg, print_msg
	import os
	import types
	from string         import replace
	from random         import shuffle, seed, randint
	# 2D alignment using rotational ccf in polar coords and linear
	# interpolation
	seed()
	print_begin_msg("ali3d_d")
	if CTF :
		from reconstruction import recons3d_4nn_ctf
		ctf_dicts = ["defocus", "Cs", "voltage", "Pixel_size", "amp_contrast", "B_factor", "ctf_applied" ]
		#                     0          1         2             3                 4                   5               6
	else   : from reconstruction import recons3d_4nn

	if os.path.exists(outdir):  os.system('rm -rf '+outdir)
	os.mkdir(outdir)

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	#if (an == "-1"):
	#	an = []
	#	for i in xrange(len(xrng)):   an.append(-1)
	#else:
	#	an = get_input_from_string(an)
	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit);

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Reference volume            : %s\n"%(ref_vol))	
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Maskfile                    : %s\n"%(maskfile))
	print_msg("Inner radius                : %i\n"%(first_ring))

	numref            = 2

	ima = EMData()
	ima.read_image(stack, 0)
	nx      = ima.get_xsize()
	if last_ring < 0:	last_ring = int(nx/2) - 2
	if(maskfile):
		if(type(maskfile) is types.StringType):	 mask3D = getImage(maskfile)
		else: 	                                mask3D = maskfile
	else        :   mask3D = model_circle(last_ring, nx, nx, nx)
	mask = model_circle(last_ring, nx, nx)

	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Ring step                   : %i\n"%(rstep))
	print_msg("X search range              : %s\n"%(xrng))
	print_msg("Y search range              : %s\n"%(yrng))
	print_msg("Translational step          : %s\n"%(step))
	print_msg("Angular step                : %s\n"%(delta))
	#print_msg("Angular search range        : %s\n"%(an))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("data with CTF               : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	print_msg("Reference projection method : %s\n"%(ref_a))
	print_msg("Symmetry group              : %s\n\n"%(symmetry))

	for  iref in xrange(numref):  dropImage(volref[iref],os.path.join(outdir, replace("ref_vol%2d.spi"%iref,' ','0')), "s")

	if(CTF):
		ctf_dicts = ["defocus", "Cs", "voltage", "Pixel_size", "amp_contrast", "B_factor", "ctf_applied" ]
		#                        0      1              2            3               4                  5               6

	del ima
	data = EMData.read_images(stack)
	nima = len(data)
	asi = []
	for im in xrange(nima):
		data[im].set_attr('ID', im)
		data[im].set_attr_dict({'group':im%numref})  # this pretends to be general.  However, for the time being no CTF
		if(CTF):
			ctf_params = get_arb_params(data[im], ctf_dicts)
			if(im == 0): data_had_ctf = ctf_params[6]
			if(ctf_params[6] == 0):
				st = Util.infomask(data[im], mask, False)
				data[im] -= st[0]
				from filter import filt_ctf
				data[im] = filt_ctf(data[im], ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
				data[im].set_attr('ctf_applied', 1)
		asi. append(randint(0,numref))
		data[im].set_attr('group',asi[im])
	volref = []
	for  iref in xrange(numref):
		# calculate initial 3D
		list_p = []
		for im in xrange(nima):
			if(data[im].get_attr('group') == iref):  list_p.append(im)
		print  iref,len(list_p)
		if(CTF): volref.append( recons3d_4nn_ctf(data, list_p, snr, 1, symmetry) )
		else:	 volref.append( recons3d_4nn(data, list_p, symmetry) )
		del list_p

	# do the projection matching
	for N_step in xrange(len(xrng)):
 		for Iter in xrange(max_iter):
			for im in xrange(nima):
				data[im].set_attr_dict({'peak':-1.0e23})
				order = range(nima)
				shuffle(order)
				for iref in xrange(numref):
					#print " ali3d_m_MPI: ITERATION #",N_step*max_iter + Iter+1
					proj_ali_incore_index(volref[iref], iref, mask3D, data[order[im]], first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], ref_a, symmetry, MPI=False)
				list_p = []
				for im in xrange(nima):
					if(data[im].get_attr('group') == asi[order[im]]):  list_p.append(im)
				print  asi[order[im]],len(list_p)
				if(CTF): volref[asi[order[im]]] = recons3d_4nn_ctf(data, list_p, snr, 1, symmetry)
				else:	 volref[asi[order[im]]] = recons3d_4nn(data, list_p, symmetry)
				del list_p
				volref[asi[order[im]]] = filt_btwl(volref[asi[order[im]]], 0.3, 0.4)
				#  changed assignment?
				group = data[order[im]].get_attr('group')
				if( asi[order[im]] != group ):
					list_p = []
					for im in xrange(nima):
						if(data[im].get_attr('group') == group):  list_p.append(im)
					print  group,len(list_p)
					if(CTF): volref[group] = recons3d_4nn_ctf(data, list_p, snr, 1, symmetry)
					else:	 volref[group] = recons3d_4nn(data, list_p, symmetry)
					del list_p
					volref[group] = filt_btwl(volref[group], 0.3, 0.4)

			soto = []
			for im in xrange(nima):
				peak = data[im].get_attr('peak')
				s2x = data[im].get_attr('s2x')
				s2y = data[im].get_attr('s2y')
				phi = data[im].get_attr('phi')
				theta = data[im].get_attr('theta')
				psi = data[im].get_attr('psi')
				group = data[im].get_attr('group')
				soto.append([phi,theta,psi,s2x,s2y,peak,group])
			from utilities import dropSpiderDoc
			dropSpiderDoc(os.path.join(outdir, replace("params%4d"%(N_step*max_iter+Iter+1),' ','0')),soto)
			for iref in xrange(numref):
				ilm = 0
				for im in xrange(nima):
					if(asi[im] == iref):  ilm +=1
				print  "  REF ",iref,"  elements:",ilm
	del  volref
	#  here we  write header info
	if(CTF and data_had_ctf == 0):
		for im in xrange(len(data)): data[im].set_attr('ctf_applied', 0)
	from utilities import write_headers
	write_headers( stack, data, range(nima))
	print_end_msg("ali3d_m")

def ali3d_m(stack, ref_vol, outdir, maskfile = None, ir=1, ou=-1, rs=1, 
            xr  ="4 2  2  1",      yr="-1",
            ts  ="1 1 0.5 0.25",   delta="10  6  4  4",
	    center = 0, maxit= 5, CTF = False, snr = 1.0,  ref_a="S", symmetry="c1", MPI=False):
	if MPI:
		ali3d_m_MPI(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, yr, ts, 
		            delta, center, maxit, CTF, snr, ref_a, symmetry)
		return
	from utilities      import model_circle, reduce_EMData_to_root, bcast_EMData_to_all, dropImage
	from utilities      import bcast_list_to_all, bcast_string_to_all, getImage, get_input_from_string
	from utilities      import get_arb_params, set_arb_params, dropSpiderDoc
	from filter	    import filt_params,filt_btwl, filt_from_fsc, filt_table
	from alignment	    import proj_ali_incore_index
	from fundamentals   import fshift
	from utilities      import print_begin_msg, print_end_msg, print_msg
	import os
	import types
	from string         import replace
	# 2D alignment using rotational ccf in polar coords and linear
	# interpolation	
	print_begin_msg("ali3d_d")
	if CTF :
		from reconstruction import recons3d_4nn_ctf
		ctf_dicts = ["defocus", "Cs", "voltage", "Pixel_size", "amp_contrast", "B_factor", "ctf_applied" ]
		#                     0          1         2             3                 4                   5               6
	else   : from reconstruction import recons3d_4nn

	if os.path.exists(outdir):  os.system('rm -rf '+outdir)
	os.mkdir(outdir)

	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	#if (an == "-1"):
	#	an = []
	#	for i in xrange(len(xrng)):   an.append(-1)
	#else:
	#	an = get_input_from_string(an)
	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit);

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Reference volume            : %s\n"%(ref_vol))	
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Maskfile                    : %s\n"%(maskfile))
	print_msg("Inner radius                : %i\n"%(first_ring))

	numref            = EMUtil.get_image_count(ref_vol)
	volref = []
	for  iref in xrange(numref):
		vol     = EMData()
		vol.read_image(ref_vol, iref)
		volref.append(vol)
	del vol

	nx      = volref[0].get_xsize()
	if last_ring < 0:	last_ring = int(nx/2) - 2
	if(maskfile):
		if(type(maskfile) is types.StringType):	 mask3D = getImage(maskfile)
		else: 	                                mask3D = maskfile
	else        :   mask3D = model_circle(last_ring, nx, nx, nx)
	mask = model_circle(last_ring, nx, nx)

	print_msg("Outer radius                : %i\n"%(last_ring))
	print_msg("Ring step                   : %i\n"%(rstep))
	print_msg("X search range              : %s\n"%(xrng))
	print_msg("Y search range              : %s\n"%(yrng))
	print_msg("Translational step          : %s\n"%(step))
	print_msg("Angular step                : %s\n"%(delta))
	#print_msg("Angular search range        : %s\n"%(an))
	print_msg("Maximum iteration           : %i\n"%(max_iter))
	print_msg("Center type                 : %i\n"%(center))
	print_msg("data with CTF               : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	print_msg("Reference projection method : %s\n"%(ref_a))
	print_msg("Symmetry group              : %s\n\n"%(symmetry))

	if(CTF):
		ctf_dicts = ["defocus", "Cs", "voltage", "Pixel_size", "amp_contrast", "B_factor", "ctf_applied" ]
		#                        0      1              2            3               4                  5               6

	data = EMData.read_images(stack)
	nima = len(data)
	for im in xrange(nima):
		data[im].set_attr('ID', im)
		data[im].set_attr_dict({'group':im%numref})  # this pretends to be general.  However, for the time being no CTF
		if(CTF):
			ctf_params = get_arb_params(data[im], ctf_dicts)
			if(im == 0): data_had_ctf = ctf_params[6]
			if(ctf_params[6] == 0):
				st = Util.infomask(data[im], mask, False)
				data[im] -= st[0]
				from filter import filt_ctf
				data[im] = filt_ctf(data[im], ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
				data[im].set_attr('ctf_applied', 1)

	soto = []
	for im in xrange(nima):
		s2x = data[im].get_attr('s2x')
		s2y = data[im].get_attr('s2y')
		phi = data[im].get_attr('phi')
		theta = data[im].get_attr('theta')
		psi = data[im].get_attr('psi')
		soto.append([phi,theta,psi,s2x,s2y])
	from utilities import dropSpiderDoc
	dropSpiderDoc(os.path.join(outdir, replace("params%4d"%(0),' ','0')),soto)
	volref = []
	for  iref in xrange(numref):
		# calculate initial 3D
		list_p = []
		for im in xrange(nima):
			if(data[im].get_attr('group') == iref):  list_p.append(im)
		print  iref,len(list_p)
		if(CTF): volref.append( recons3d_4nn_ctf(data, list_p, snr, 1, symmetry) )
		else:	 volref.append( recons3d_4nn(data, list_p, symmetry) )
		del list_p

	for  iref in xrange(numref):  dropImage(volref[iref],os.path.join(outdir, replace("ref_vol%2d.spi"%iref,' ','0')), "s")

	from utilities import ttime
	# do the projection matching
	for N_step in xrange(len(xrng)):
 		for Iter in xrange(max_iter):
			for im in xrange(nima):  data[im].set_attr_dict({'peak':-1.0e23})
			for iref in xrange(numref):
				print  ttime()
				#print " ali3d_m_MPI: ITERATION #",N_step*max_iter + Iter+1
				proj_ali_incore_index(volref[iref], iref, mask3D, data, first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], ref_a, symmetry, MPI=False)
			soto = []
			for im in xrange(nima):
				peak = data[im].get_attr('peak')
				s2x = data[im].get_attr('s2x')
				s2y = data[im].get_attr('s2y')
				phi = data[im].get_attr('phi')
				theta = data[im].get_attr('theta')
				psi = data[im].get_attr('psi')
				group = data[im].get_attr('group')
				soto.append([phi,theta,psi,s2x,s2y,peak,group])
			from utilities import dropSpiderDoc
			dropSpiderDoc(os.path.join(outdir, replace("params%4d"%(N_step*max_iter+Iter+1),' ','0')),soto)
			for iref in xrange(numref):
				"""
				print  ttime()
				list_p = []
				get = 1
				for im in xrange(nima):
					if(data[im].get_attr('group') == iref):
						if(get):
							get = 0
							list_p.append(im)
						else:    get = 1
 				if(CTF): vol1 = recons3d_4nn_ctf(data, list_p, snr, 1, symmetry)
				else:    vol1 = recons3d_4nn(data, list_p, symmetry)
				list_p = []
				get = 0
				for im in xrange(nima):
					if(data[im].get_attr('group') == iref):
						if(get):
							get = 0
							list_p.append(im)
						else:    get = 1
				if(CTF): vol2 = recons3d_4nn_ctf(data, list_p, snr, 1, symmetry)
				else:	 vol2 = recons3d_4nn(data, list_p, symmetry)
				if(mask3D == None):
					fscc = fsc( vol1, vol2, 1.0, os.path.join(outdir, replace("resolution%2d%4d"%(iref, N_step*max_iter+Iter+1),' ','0')))
				else:
					sbo = Util.infomask(vol1, mask3D, False)
					sbe = Util.infomask(vol2, mask3D, False)
					fscc = fsc( (vol1 - sbo[0])*mask3D, (vol2 - sbe[0])*mask3D, 1.0, os.path.join(outdir, replace("resolution%2d%4d"%(iref, N_step*max_iter+Iter+1),' ','0')))

				del vol1
				del vol2
				"""
				print  ttime()
				# calculate new and improved 3D
				list_p = []
				for im in xrange(nima):
					if(data[im].get_attr('group') == iref):  list_p.append(im)
				print  iref,len(list_p)
				if(CTF): volref[iref] = recons3d_4nn_ctf(data, list_p, snr, 1, symmetry)
				else:	 volref[iref] = recons3d_4nn(data, list_p, symmetry)
				del list_p
				print ttime()
				dropImage(volref[iref],os.path.join(outdir, replace("vol%2d%4d.spi"%(iref, N_step*max_iter+Iter+1),' ','0')), "s")
				#if(fscc[1][0] < 0.5):  fscc[1][0] = 1.0
				#if(fscc[1][1] < 0.5):  fscc[1][1] = 1.0
				#fl, fh = filt_params(fscc)
				#filt = filt_from_fsc(fscc, 0.075)
				# here figure the filtration parameters and filter vol for the  next iteration
				#fl, fh = filt_params(res)
				#lk = 2
				#while(fscc[1][lk] >0.98 and fscc[0][lk]<0.25):
				#	lk+=1
				#fl = fscc[0][lk]
				#fh = min(fl+0.1,0.49)
				#print "fl, fh, iter",fl,fh,Iter
				#volref[iref] = filt_btwl(volref[iref], fl, fh)
				volref[iref] = filt_btwl(volref[iref], 0.2, 0.3)
				if(center == 1):
					cs   = volref[iref].phase_cog()
					volref[iref]  = fshift(volref[iref], -cs[0], -cs[1] -cs[2])
				dropImage(volref[iref],os.path.join(outdir, replace("volf%2d%4d.spi"%(iref, N_step*max_iter+Iter+1),' ','0')), "s")
	del  volref
	#  here we  write header info
	if(CTF and data_had_ctf == 0):
		for im in xrange(len(data)): data[im].set_attr('ctf_applied', 0)
	from utilities import write_headers
	write_headers( stack, data, range(nima))
	print_end_msg("ali3d_m")


def ali3d_m_MPI(stack, ref_vol, outdir, maskfile = None, ir=1, ou=-1, rs=1, 
            xr              ="4 2  2  1",      yr="-1",
            ts="1 1 0.5 0.25",   delta="10  6  4  4", an="-1",
	    center = 0, maxit= 5, CTF = False, snr = 1.0,  ref_a="S", symmetry="c1"):
	from utilities      import model_circle, reduce_EMData_to_root, bcast_EMData_to_all, dropImage
	from utilities      import bcast_list_to_all, bcast_string_to_all, getImage, get_input_from_string
	from utilities      import get_arb_params, set_arb_params, dropSpiderDoc, recv_attr_dict, send_attr_dict
	from filter	    import filt_params, filt_btwl, filt_from_fsc2, filt_table, fit_tanh, filt_tanl, filt_vols
	from alignment	    import proj_ali_incore_index
	from random         import randint
	from fundamentals   import fshift
	from utilities      import print_begin_msg, print_end_msg, print_msg
	import os
	import types
	from string         import replace
	from mpi 	    import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	    import mpi_barrier

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	if(myid == main_node):
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)


	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	if (an == "-1"):
		an = []
		for i in xrange(len(xrng)):   an.append(-1)
	else:
		from  alignment	    import proj_ali_incore_local
		an      = get_input_from_string(an)

	from alignment import proj_ali_incore

	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(maxit)
	center      = int(center)

	numref            = EMUtil.get_image_count(ref_vol)
	volref = []
	for  iref in xrange(numref):
		vol     = EMData()
		vol.read_image(ref_vol, iref)
		volref.append(vol)
	del vol
	nx      = volref[0].get_xsize()
	if last_ring < 0:	last_ring = int(nx/2) - 2
	if(maskfile):
		if(type(maskfile) is types.StringType):	 mask3D = getImage(maskfile)
		else: 	                                mask3D = maskfile
	else        :   mask3D = model_circle(last_ring, nx, nx, nx)
	mask = model_circle(last_ring, nx, nx)

	nima = EMUtil.get_image_count(stack)

	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)

	finfo = open( "progress%04d"%myid, "w" )
        finfo.write( "image_start, image_end: %6d %6d\n" %(image_start, image_end) )
        finfo.flush()

	#if(myid == main_node)       :
	#	for  iref in xrange(numref):
	#		volref[iref].write_image(os.path.join(outdir,replace("ref_vol%2d.hdf"%iref,' ','0')))
	if(CTF):
		ctf_dicts = ["defocus", "Cs", "voltage", "Pixel_size", "amp_contrast", "B_factor", "ctf_applied" ]
		#                0       1        2            3               4            5               6
		from reconstruction import rec3D_MPI
	else:
		from reconstruction import rec3D_MPI_noCTF

	data = EMData.read_images(stack, range(image_start, image_end))
	for im in xrange(image_start, image_end):
		data[im].set_attr('ID', im)
		group = data[im].get_attr_default('group', -1)
		if( group == -1 ):  data[im].set_attr_dict({'group':im%numref})
		if(CTF):
			ctf_params = get_arb_params(data[im], ctf_dicts)
			if(im == image_start): data_had_ctf = ctf_params[6]
			if(ctf_params[6] == 0):
				st = Util.infomask(data[im], mask, False)
				data[im] -= st[0]
				from filter import filt_ctf
				data[im] = filt_ctf(data[im], ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
				data[im].set_attr('ctf_applied', 1)
				#set_arb_params(data[im], ctf_dicts)  #I am not sure whether this is needed
	finfo.write( "data loaded\n" )
	finfo.flush()

	par_str = ['phi', 'theta', 'psi', 's2x', 's2y', 'group']
	# do the projection matching
	for N_step in xrange(len(xrng)):
 		for Iter in xrange(max_iter):
			if(myid == main_node):
				msg = " ITERATION #%2i  X range = %5.2f   Y range = %5.2f   Step = %4.2f  delta = %4.1f  angular range = %4.1f\n"%(N_step*max_iter + Iter+1, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], an[N_step])
				finfo.write( msg + "\n" )
				finfo.flush()
			old_params  = []
			best_params = []
			for im in xrange(image_start, image_end):
				data[im-image_start].set_attr_dict({'peak_b':-1.0e23})
				qtemp = get_arb_params(data[im-image_start], par_str)
				old_params.append( qtemp )
				best_params.append( qtemp )
			for iref in xrange(numref):
				if(myid == main_node):
					msg = "  Align particles for group #%d \n"%(iref)
					finfo.write(msg)
					finfo.flush()
				if(an[N_step] == -1): proj_ali_incore(volref[iref], mask3D, data, first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], ref_a, symmetry, CTF, finfo, MPI=True)
				else:                proj_ali_incore_local(volref[iref], mask3D, data, first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], an[N_step], ref_a, symmetry, CTF, finfo, MPI=True)
				for im in xrange(image_start, image_end):
					peaks = get_arb_params(data[im-image_start], ['peak', 'peak_b'] )
					finfo.write( "imn, grpid, peak, peak_best: %6d %d %f %f\n" % (im, iref, peaks[0], peaks[1]) )
					finfo.flush()
					if( peaks[0] > peaks[1] ):
						data[im-image_start].set_attr_dict({'peak_b': peaks[0]})
						data[im-image_start].set_attr_dict({'group': iref})
						best_params[im-image_start] = get_arb_params(data[im-image_start], par_str)
					set_arb_params(data[im-image_start], old_params[im-image_start], par_str)

			for im in xrange(image_start, image_end):
				set_arb_params(data[im-image_start], best_params[im-image_start], par_str)
				qtemp = data[im-image_start].get_attr('peak_b')
				data[im-image_start].set_attr_dict({'peak': qtemp})
				finfo.write( "image %6d assigned to group %d, peak value is %f\n" % (im, best_params[im-image_start][-1], qtemp) )
				finfo.flush()
			dropSpiderDoc(os.path.join(outdir, "new_params%03d_%03d"%(Iter,myid)), best_params, "phi, theta, psi, s2x, s2y, group")
			del old_params, best_params
			mpi_barrier(MPI_COMM_WORLD)

			fscc = [None]*numref
			for iref in xrange(numref):
				finfo.write( "reconstructing %d\n" % iref )
				finfo.flush()
				if(CTF): 
                                    volref[iref], fscc[iref] = rec3D_MPI(data, snr, symmetry, mask3D, os.path.join(outdir, "resolution%03d_%03d"%(iref, Iter)), myid, main_node, info=outf, index=iref)
				else:    
				    volref[iref], fscc[iref] = rec3D_MPI_noCTF(data, symmetry, mask3D, os.path.join(outdir, "resolution%03d_%03d"%(iref,Iter)), myid, main_node, 1.0, index = iref)
				if(myid == main_node):
					dropImage(volref[iref],os.path.join(outdir, "vol%03d_%03d.spi"%(iref, Iter)), 's')
					if(center == 1):
						cs   = volref[iref].phase_cog()
						volref[iref]  = fshift(volref[iref], -cs[0], -cs[1] -cs[2])
			mpi_barrier(MPI_COMM_WORLD)

			if(myid==main_node):
				volref = filt_vols( volref, fscc, mask3D )
				for iref in xrange(numref):
					dropImage( volref[iref], os.path.join(outdir, "volf%03d_%03d.spi"%(iref, Iter)), 's' )

			for iref in xrange(numref): bcast_EMData_to_all(volref[iref], myid, main_node)

		if(CTF and data_had_ctf == 0):
			for im in xrange(image_start, image_end): data[im-image_start].set_attr('ctf_applied', 0)
		if(myid == main_node):
			result_stack = 'ali3d_m_result.hdf'
			recv_attr_dict(main_node, result_stack, data, par_str, image_start, image_end, number_of_proc)
		else: send_attr_dict(main_node, data, par_str, image_start, image_end)


	mpi_barrier(MPI_COMM_WORLD)
	del  volref
	vol, fscc = rec3D_MPI_noCTF(data, symmetry, mask3D, os.path.join(outdir, "resolution_merged"), myid, main_node, info=myinfo)
	if(myid == main_node):  dropImage(vol,os.path.join(outdir, "vol_merged.spi"), "s")
	#  clean up
	del vol
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
	if(CTF and data_had_ctf == 0):
		for im in xrange(len(data)): data[im].set_attr('ctf_applied', 0)
	if(myid == main_node): recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	else: send_attr_dict(main_node, data, par_str, image_start, image_end)

def ali3d_em_MPI_origin(stack, ref_vol, outdir, maskfile, ou=-1,  delta=2, maxit=10, CTF = None, snr=1.0, sym="c1"):
	"""
		
	"""
	from alignment	    import eqprojEuler
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl, filt_tanl
	from fundamentals   import fshift, rot_avg_image
	from projection     import prep_vol, prgs
	from utilities      import amoeba, bcast_list_to_all, bcast_string_to_all, model_circle, get_arb_params, set_arb_params, dropSpiderDoc
	from utilities      import getImage, dropImage, bcast_EMData_to_all, send_attr_dict, recv_attr_dict
	from utilities      import readSpiderDoc, get_im
	from reconstruction import rec3D_MPI
	from statistics     import ccc
	from math           import pi, sqrt
	from string         import replace
	import os
	import sys
	from development    import ali_G3
	from mpi 	    import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	    import mpi_barrier

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	#  chose a random node as a main one...
	main_node = 0
	if(myid == main_node):
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)
	mpi_barrier(MPI_COMM_WORLD)

	nima = EMUtil.get_image_count(stack)

	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
	
	outf = file("progress%04d"%myid, "w")

        outf.write( "image_start, image_end: %6d %6d\n" %(image_start, image_end) )
        outf.flush()

	#  We heave Kref reference volumes
	Kref = EMUtil.get_image_count(ref_vol)
	vol = []
	for krf in xrange(Kref):
		vol.append(get_im(ref_vol, krf))
	nx  = vol[0].get_xsize()

	if(ou <= 0):  ou = nx//2-2

	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask3D = getImage(maskfile)
		else:                                   mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)
	mask2D = model_circle(ou, nx, nx)

	parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
	#                  0                 1            2            3            4                 5                 6
	#from utilities import readSpiderDoc, set_arb_params
	#prm = readSpiderDoc("params_new.doc")
	#prm_dict = ["phi", "theta", "psi", "s2x", "s2y"]
	dataim = EMData.read_images(stack, range(image_start, image_end))
	for im in xrange(image_start, image_end):
		dataim[im].set_attr('ID', im)
		#set_arb_params(dataim[im], prm[im], prm_dict)
		if(CTF):
			ctf_params = get_arb_params(dataim[im], parnames)
			if(im == image_start): data_had_ctf = ctf_params[6]
			if(ctf_params[6] == 0):
				st = Util.infomask(dataim[im], mask2D, False)
				dataim[im] -= st[0]
				from filter import filt_ctf
				dataim[im] = filt_ctf(dataim[im], ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
				dataim[im].set_attr('ctf_applied', 1)
		group = dataim[im].get_attr_default('group', 0)
		dataim[im].set_attr('group', group)
	outf.write("data read\n")
	outf.flush()

	from utilities      import bcast_number_to_all
	from morphology     import threshold, threshold_to_minval

	par_str=["phi", "theta", "psi", "s2x", "s2y"]

	from time import time
	ldef = -1
	for iteration in xrange(maxit):
		outf.write("iteration = "+str(iteration)+"   ")
		outf.write("\n")
		outf.flush()

		volft = [None]*Kref
		for krf in xrange(Kref):
			volft[krf],kb  = prep_vol(vol[krf])
		outf.write( "prepare volumes done\n" )
		outf.flush()

		iter_start_time = time()

		for imn in xrange(image_start, image_end):

			img_start_time = time()
			
			atparams = get_arb_params(dataim[imn-image_start], par_str)
			weight_phi = max(delta, delta*abs((atparams[1]-90.0)/180.0*pi))
			best = -1.0
			for krf in xrange(Kref):
				data = []
				#data.append(volft)
				#data.append(kb)
				#data.append(mask2D)
				# Only Euler angles
				data.append(volft[krf])
				data.append(kb)
				data.append(dataim[imn-image_start])

	
				#optm_params = ali_G3(data, atparams, dtheta)
				#  Align only Euler angles
				#  change signs of shifts for projections
				data.insert(3, -atparams[3])
				data.insert(4, -atparams[4])
				"""
				initial  = eqproj(atparams, data)  # this is if we need initial discrepancy
				outf.write("Image "+str(imn)+"\n")
				outf.write('Old %6.1f  %6.1f  %6.1f  %6.1f  %6.1f   %7.4f '%(atparams[0],atparams[1],atparams[2],-atparams[3],-atparams[4], initial))
				outf.write("\n")
				outf.flush()
				"""
				data.insert(5, mask2D)
				#from utilities import start_time, finish_time
				#t3=start_time()
				optm_params =  amoeba(atparams[0:3], [weight_phi, delta, weight_phi], eqprojEuler, 1.e-4,1.e-4,500, data)
				if(optm_params[1] > best):
					best = optm_params[1]
					ngroup = krf
					best_optm_params = optm_params
			best_optm_params[0].append(imn)																													     
			#new_params.append(optm_params[0])		

			'''														     
			outf.write('New %6.1f  %6.1f  %6.1f  %6.1f  %6.1f   %7.4f    %d4   %7.1f\n'%(optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4],optm_params[1], optm_params[2], ctf_params[1]))						     
			outf.write("\n")																															     
			outf.flush()				
			'''													     
																																				     
			del  data
			ima.set_attr('group', ngroup)
			set_arb_params(dataim[imn-image_start], [best_optm_params[0][0], best_optm_params[0][1], best_optm_params[0][2]], par_str[0:3])

                        end_time = time()

			outf.write( "imn, time, all time: %6d %10.3f %10.3f\n" % (imn, end_time-img_start_time, end_time-iter_start_time) )
			outf.flush()


			
			#  change signs of shifts for projections
			#atparams[3] *= -1
			#atparams[4] *= -1
			
			'''
			initial  = eqproj(atparams, data)  # this is if we need initial discrepancy
			outf.write("Image "+str(imn)+"\n")
			outf.write('Old %6.1f  %6.1f  %6.1f  %6.1f  %6.1f   %7.4f '%(atparams[0],atparams[1],atparams[2],-atparams[3],-atparams[4], initial))
			outf.write("\n")
			outf.flush()
			'''
			#weight_phi = max(delta, delta*abs((atparams[1]-90.0)/180.0*pi))
			#from utilities import start_time, finish_time
			#t3=start_time()
			#optm_params =  amoeba(atparams, [weight_phi, delta, weight_phi, 0.05, 0.05], eqproj, 1.e-4,1.e-4,500,data)
			#  change signs of shifts for header
			#optm_params[0][3] *= -1
			#optm_params[0][4] *= -1
			#optm_params[0].append(imn)
			#new_params.append(optm_params[0])
			
			'''
			outf.write('New %6.1f  %6.1f  %6.1f  %6.1f  %6.1f   %7.4f    %d4   %7.1f'%(optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4],optm_params[1], optm_params[2], ctf_params[1]))
			outf.write("\n")
			outf.flush()
			'''
			
			#del  data[2]
			#set_arb_params(dataim[imn-image_start], [optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4]], par_str)
			#t4 = finish_time(t3)
		soto = []
		for imn in xrange(image_start, image_end):
			s2x   = dataim[imn-image_start].get_attr('s2x')
			s2y   = dataim[imn-image_start].get_attr('s2y')
			phi   = dataim[imn-image_start].get_attr('phi')
			theta = dataim[imn-image_start].get_attr('theta')
			psi   = dataim[imn-image_start].get_attr('psi')
			group   = dataim[imn-image_start].get_attr('group')
			soto.append([phi,theta,psi,s2x,s2y,group,imn])
		dropSpiderDoc(os.path.join(outdir, replace("new_params%3d_%3d_%3d"%(iteration, ic, myid),' ','0')), soto," phi, theta, psi, s2x, s2y, group, image number")
		del soto

		fscc = [None]*Kref
		for krf in xrange(Kref):
 	    		# resolution
			outf.write("  begin reconstruction = "+str(image_start)+"   ")
			outf.write("\n")
			outf.flush()
			vol[krf], fscc[krf] = rec3D_MPI(dataim, snr, sym, mask3D, os.path.join(outdir, replace("resolution%3d_%3d_%3d"%(krf, iteration, ic),' ','0') ), myid, main_node, index = krf)
			outf.write("  done reconstruction = "+str(image_start)+"   ")
			outf.write("\n")
			outf.flush()
		mpi_barrier(MPI_COMM_WORLD)
		if(myid == main_node):
			flm = 1.0
			for krf in xrange(Kref):
				dropImage(vol[krf], os.path.join(outdir, replace("vol%3d_%3d_%3d.hdf"%(krf, iteration, ic),' ','0') ))
				stat = Util.infomask(vol[krf], mask3D, False)
				vol -= stat[0]
				vol /= stat[1]
				vol[krf] = threshold(vol[krf])
				fl, aa = fit_tanh(fscc[krf])
				if(fl < flm):
					flm = fl
					aam = aa
			outf.write('tanh params %8.4f  %8.4f '%(flm, aam))
			outf.write("\n")
			outf.flush()
			for krf in xrange(Kref):
				vol[krf] = filt_tanl(vol[krf], flm, aam)
				dropImage(vol[krf], os.path.join(outdir, replace("volf%3d_%3d_%3d.hdf"%(krf, iteration, ic),' ','0') ))
		for krf in xrange(Kref):
			bcast_EMData_to_all(vol[krf], myid, main_node)

		#  here we should write header info, just in case the program crashes...
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)



def ali3d_em_MPI(stack, ref_vol, outdir, maskfile, ou=-1,  delta=2, maxit=10, CTF = None, snr=1.0, sym="c1", chunk=0.101):
	"""
		
	"""
	from alignment	    import eqprojEuler
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl, filt_tanl, filt_vols
	from fundamentals   import fshift, rot_avg_image
	from projection     import prep_vol, prgs
	from utilities      import amoeba, bcast_list_to_all, bcast_string_to_all, model_circle, get_arb_params, set_arb_params, dropSpiderDoc
	from utilities      import getImage, dropImage, bcast_EMData_to_all, send_attr_dict, recv_attr_dict
	from utilities      import readSpiderDoc, get_im
	from reconstruction import rec3D_MPI
	from statistics     import ccc
	from math           import pi, sqrt
	from string         import replace
	import os
	import sys
	from development    import ali_G3
	from mpi 	    import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	    import mpi_barrier

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	#  chose a random node as a main one...
	main_node = 0
	if(myid == main_node):
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)
	mpi_barrier(MPI_COMM_WORLD)

	nima = EMUtil.get_image_count(stack)

	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
	
	outf = file( outdir+("/progress%04d"%myid), "w")
        outf.write( "image_start, image_end: %6d %6d\n" %(image_start, image_end) )
        outf.flush()

        # refine step define on which step refinement will be carried
        # if set to -1, new refinement only assignment 
        total_step=5
	refine_step=0
	nx  = getImage( ref_vol ).get_xsize()

	if maskfile:
		import  types
		if(type(maskfile) is types.StringType): 
			outf.write( "usage 3D mask " + maskfile  + "\n")
			outf.flush()                          
			mask3D = getImage(maskfile)
		else:   
			mask3D = maskfile
	else:
		outf.write( "using spherical 3D mask: rad=%d\n" % ou )
		outf.flush()
		mask3D = model_circle(ou, nx, nx, nx)
	mask2D = model_circle(ou, nx, nx)

	#  We heave Kref reference volumes
	Kref = EMUtil.get_image_count(ref_vol)
	vol = [None]*Kref
	for krf in xrange(Kref):
		vol[krf] = get_im(ref_vol, krf)
	if(ou <= 0):  ou = nx//2-2

	dataim = []
	parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
	#               0             1         2        3          4              5                 6
	#from utilities import readSpiderDoc, set_arb_params
	#prm = readSpiderDoc("params_new.doc")
	#prm_dict = ["phi", "theta", "psi", "s2x", "s2y"]

	for im in xrange(image_start, image_end):
		ima = EMData()
		ima.read_image(stack, im)
		if(CTF):
			ctfparms = get_arb_params(ima, parnames)
			if(ctfparms[6] == 1.0):
				print "Error: ali3d_em does not work on ctf_applied data"
				return

		group = ima.get_attr_default('group', 0)
		ima.set_attr('group', group)
		dataim.append(ima)


	outf.write("data read\n")
	outf.flush()

	from utilities      import bcast_number_to_all, info
	from morphology     import threshold, threshold_to_minval

	par_str=["phi", "theta", "psi", "s2x", "s2y"]

	from time import time
	ldef = -1
	for iteration in xrange(maxit):
		outf.write("iteration = "+str(iteration)+"   ")
		outf.write("\n")
		outf.flush()

		volft = [None]*Kref
		iter_start_time = time()

		if iteration%total_step==refine_step:
			mychunk = chunk
			nchunk = int(1.0/chunk) + 1
		else:
			mychunk = 1.0
			nchunk = 1

		outf.write( "N of chunks: %4d\n" % nchunk )
		outf.flush()

		curt_defocus = -1.0

		for ichunk in xrange(nchunk):
			image_chunk_start = image_start + int((image_end-image_start)*mychunk*ichunk)
			image_chunk_end   = image_start + int((image_end-image_start)*mychunk*(ichunk+1))
			if image_chunk_end > image_end:
				image_chunk_end = image_end

			for imn in xrange(image_chunk_start, image_chunk_end):
				ima = dataim[imn-image_start]

				if CTF:
					ctfparms = get_arb_params( ima, parnames )
					if ctfparms[1] != curt_defocus :
						curt_defocus = ctfparms[1]
						
						for krf in xrange(Kref):
							ctfvol = filt_ctf(vol[krf],ctfparms[1],ctfparms[3],ctfparms[2],ctfparms[0],ctfparms[4],ctfparms[5])
							volft[krf],kb  = prep_vol( vol[krf] )



				img_start_time = time()
				atparams = get_arb_params(dataim[imn-image_start], par_str)
				weight_phi = max(delta, delta*abs((atparams[1]-90.0)/180.0*pi))
				bestccc = -2.0
				bestgrp = None
				bestang = [None]*3
				for krf in xrange(Kref):

					if iteration%total_step==refine_step:
						# running refinement
						data = []
						data.append(volft[krf])
						data.append(kb)
						data.append(ima)     
						data.append(-atparams[3])
						data.append(-atparams[4])
						data.append(mask2D)
						[curtang,curtccc,niter] = amoeba(atparams[0:3], [weight_phi, delta, weight_phi], eqprojEuler, 1.e-4, 1.e-4, 500, data)
					else:
						# change assignment only
						prj = prgs( volft[krf], kb, [ atparams[0], atparams[1], atparams[2], -atparams[3], -atparams[4] ] )
						curtccc = prj.cmp( "ccc", ima, {"mask":mask2D, "negative":0} )
						outf.write( "        imn, krf, ccc: %d %d %10.5f\n" % (imn, krf, curtccc) )
						outf.flush()

					if(bestccc < curtccc):
						bestccc = curtccc
						bestgrp = krf
						if iteration%total_step==refine_step: bestang = curtang[0:3]																										     
				ima.set_attr('group', bestgrp)

				if iteration%total_step==refine_step:
					ima.set_attr( 'phi', bestang[0])
					ima.set_attr( 'theta', bestang[1] )
					ima.set_attr( 'psi', bestang[2] )
					outf.write( "imn,old: %6d %10.3f %10.3f %10.3f\n" % (imn, atparams[0], atparams[1], atparams[2]) )
					outf.write( "imn,new: %6d %10.3f %10.3f %10.3f %d iter %d\n" % (imn, bestang[0], bestang[1], bestang[2], bestgrp, iteration) )
					outf.flush()
				else:
					outf.write( "imn %6d assigned to group %d on iter %d\n" % (imn, bestgrp, iteration) )
					outf.flush()

			if iteration%total_step==refine_step:
				soto = []
				for imn in xrange(image_chunk_start, image_chunk_end):
					s2x   = dataim[imn-image_start].get_attr('s2x')
					s2y   = dataim[imn-image_start].get_attr('s2y')
					phi   = dataim[imn-image_start].get_attr('phi')
					theta = dataim[imn-image_start].get_attr('theta')
					psi   = dataim[imn-image_start].get_attr('psi')
					group = dataim[imn-image_start].get_attr('group')
					soto.append([phi,theta,psi,s2x,s2y,group,imn])
				dropSpiderDoc(os.path.join(outdir, "new_params%03d_%03d_%03d"%(iteration, ichunk, myid)), soto," phi, theta, psi, s2x, s2y, group, image number")


			fscc = [None]*Kref
			for krf in xrange(Kref):
 	    			# resolution
				outf.write("reconstructing volume: %d\n" % krf)
				outf.flush()
				vol[krf], fscc[krf] = rec3D_MPI(dataim, snr, sym, mask3D, os.path.join(outdir, "resolution%03d_%03d_%03d"%(krf, iteration, ichunk)), myid, main_node, index = krf)
				if(myid==main_node):
					dropImage( vol[krf], os.path.join(outdir,"vol%03d_%03d_%03d.spi"%(krf, iteration, ichunk)), 's')
			mpi_barrier(MPI_COMM_WORLD)

			if(myid == main_node):
				vol = filt_vols( vol, fscc, mask3D )
				for krf in xrange(Kref):
					dropImage( vol[krf], os.path.join(outdir,"volf%03d_%03d_%03d.spi"%(krf, iteration, ichunk)), 's')

			for krf in xrange(Kref):
				outf.write( "bcasting volume: %d\n" % krf )
				outf.flush()
				bcast_EMData_to_all(vol[krf], myid, main_node)

		#from sys import exit
		#exit()
		#  here we should write header info, just in case the program crashes...
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
	if(myid == main_node): recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	else: send_attr_dict(main_node, data, par_str, image_start, image_end)

def ali3d_en_MPI(stack, ref_vol, outdir, maskfile, ou=-1,  delta=2, maxit=10, CTF = None, snr=1.0, sym="c1", chunk=0.101):
	"""
		
	"""
	from alignment	    import eqprojEuler
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl, filt_tanl
	from fundamentals   import fshift, rot_avg_image
	from projection     import prep_vol, prgs
	from utilities      import amoeba, bcast_list_to_all, bcast_string_to_all, model_circle, get_arb_params, set_arb_params, dropSpiderDoc
	from utilities      import getImage, dropImage, bcast_EMData_to_all, send_attr_dict, recv_attr_dict
	from utilities      import readSpiderDoc, get_im
	from reconstruction import rec3D_MPI
	from statistics     import ccc
	from math           import pi, sqrt
	from string         import replace
	import os
	import sys
	from development    import ali_G3
	from mpi 	    import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	    import mpi_barrier

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	#  chose a random node as a main one...
	main_node = 0
	if(myid == main_node):
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)
	mpi_barrier(MPI_COMM_WORLD)

	nima = EMUtil.get_image_count(stack)

	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
	
	outf = file("progress%04d"%myid, "w")
        outf.write( "image_start, image_end: %6d %6d\n" %(image_start, image_end) )
        outf.flush()


	nx  = getImage( ref_vol ).get_xsize()

	if maskfile:
		import  types
		if(type(maskfile) is types.StringType): 
			outf.write( "usage 3D mask " + maskfile  + "\n")
			outf.flush()                          
			mask3D = getImage(maskfile)
		else:   
			mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)
	mask2D = model_circle(ou, nx, nx)

	#  We heave Kref reference volumes
	Kref = EMUtil.get_image_count(ref_vol)
	vol = [None]*Kref
	for krf in xrange(Kref):
		vol[krf] = get_im(ref_vol, krf)


	if(ou <= 0):  ou = nx//2-2

	dataim = []
	parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
	#               0             1         2        3          4              5                 6
	#from utilities import readSpiderDoc, set_arb_params
	#prm = readSpiderDoc("params_new.doc")
	#prm_dict = ["phi", "theta", "psi", "s2x", "s2y"]
	param_file = "first_ali3d_run4/new_params004_000_%03d" % myid 
	soto = readSpiderDoc( param_file )
	outf.write( "reading params from " + param_file + "\n")
	outf.flush()

	data = EMData.read_images(stack, range(image_start, image_end))
	for im in xrange(image_start, image_end):
		dataim[im].set_attr('ID', im)
		dataim[im].set_attr( 'phi', soto[im-image_start][0] )
		dataim[im].set_attr( 'theta', soto[im-image_start][1] )
		dataim[im].set_attr( 'psi', soto[im-image_start][2] )
		dataim[im].set_attr( 's2x', soto[im-image_start][3] )
		dataim[im].set_attr( 's2y', soto[im-image_start][4] )
		dataim[im].set_attr( 'group', int(soto[im-image_start][5]) )

		#set_arb_params(dataim[im], prm[im], prm_dict)
		if(CTF):
			ctf_params = get_arb_params(dataim[im], parnames)
			if(im == image_start): data_had_ctf = ctf_params[6]
			if(ctf_params[6] == 0):
				st = Util.infomask(dataim[im], mask2D, False)
				dataim[im] -= st[0]
				from filter import filt_ctf
				dataim[im] = filt_ctf(dataim[im], ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
				dataim[im].set_attr('ctf_applied', 1)

	outf.write("data read\n")
	outf.flush()

	from utilities      import bcast_number_to_all
	from morphology     import threshold, threshold_to_minval

	par_str=["phi", "theta", "psi", "s2x", "s2y"]

	from time import time
	ldef = -1
	for iteration in xrange(maxit):
		outf.write("iteration = "+str(iteration)+"   ")
		outf.write("\n")
		outf.flush()

		volft = [None]*Kref
		for krf in xrange(Kref):
			volft[krf],kb  = prep_vol( vol[krf] )
		outf.write( "prepare volumes done\n" )
		outf.flush()

		iter_start_time = time()
		mychunk = 1.0
		nchunk = 1

		outf.write( "N of chunks: %4d\n" % nchunk )
		outf.flush()

		for ichunk in xrange(nchunk):
			image_chunk_start = image_start + int((image_end-image_start)*mychunk*ichunk)
			image_chunk_end   = image_start + int((image_end-image_start)*mychunk*(ichunk+1))
			if image_chunk_end > image_end:
				image_chunk_end = image_end

			for imn in xrange(image_chunk_start, image_chunk_end):
				ima = dataim[imn-image_start]
				img_start_time = time()
				atparams = get_arb_params(dataim[imn-image_start], par_str)
				weight_phi = max(delta, delta*abs((atparams[1]-90.0)/180.0*pi))
				
				gid = ima.get_attr( "group" )
				data = []
				data.append(volft[gid])
				data.append(kb)
				data.append(ima)     
				data.append(-atparams[3])
				data.append(-atparams[4])
				data.append(mask2D)
				[bestang,curtccc,niter] = amoeba(atparams[0:3], [weight_phi, delta, weight_phi], eqprojEuler, 1.e-4, 1.e-4, 500, data)
			
				ima.set_attr( 'phi', bestang[0])
				ima.set_attr( 'theta', bestang[1] )
				ima.set_attr( 'psi', bestang[2] )
				outf.write( "imn,old: %6d %6.1f %6.1f %6.1f\n" % (imn, atparams[0], atparams[1], atparams[2]) )
				outf.write( "imn,new: %6d %6.1f %6.1f %6.1f\n" % (imn, bestang[0], bestang[1], bestang[2]) )
				outf.flush()

			soto = []
			for imn in xrange(image_chunk_start, image_chunk_end):
				s2x   = dataim[imn-image_start].get_attr('s2x')
				s2y   = dataim[imn-image_start].get_attr('s2y')
				phi   = dataim[imn-image_start].get_attr('phi')
				theta = dataim[imn-image_start].get_attr('theta')
				psi   = dataim[imn-image_start].get_attr('psi')
				group = dataim[imn-image_start].get_attr('group')
				soto.append([phi,theta,psi,s2x,s2y,group,imn])
			dropSpiderDoc(os.path.join(outdir, "new_params%03d_%03d_%03d"%(iteration, ichunk, myid)), soto," phi, theta, psi, s2x, s2y, group, image number")


			fscc = [None]*Kref
			for krf in xrange(Kref):
 	    			# resolution
				outf.write("reconstructing volume: %d\n" % krf)
				outf.flush()
				vol[krf], fscc[krf] = rec3D_MPI(dataim, snr, sym, mask3D, os.path.join(outdir, "resolution%03d_%03d_%03d"%(krf, iteration, ichunk)), myid, main_node, index = krf)
				if(myid==main_node):
					dropImage( vol[krf], os.path.join(outdir,"vol%03d_%03d_%03d.spi"%(krf, iteration, ichunk)), 's')
			mpi_barrier(MPI_COMM_WORLD)

			if(myid == main_node):
				from filter import fit_tanh, filt_tanl
				for krf in xrange(Kref):
					fl, aa = fit_tanh( fscc[krf] )
					vol[krf] = filt_tanl( vol[krf], fl, aa )
				
				for krf in xrange(Kref):
					dropImage( vol[krf], os.path.join(outdir,"volf%03d_%03d_%03d.spi"%(krf, iteration, ichunk)), 's')

			for krf in xrange(Kref):
				outf.write( "bcasting volume: %d\n" % krf )
				outf.flush()
				bcast_EMData_to_all(vol[krf], myid, main_node)

		#from sys import exit
		#exit()
		#  here we should write header info, just in case the program crashes...
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)

def ali3d_e(stack, ref_vol, outdir, maskfile = None, ou = -1,  delta = 2, maxit = 10, CTF = False, snr = 1.0, sym="c1", chunk = -1.0, user_func_name="ref_ali3d", info = True, MPI = False):
	"""
		
	"""
	if MPI:
		ali3d_e_MPI(stack, ref_vol, outdir, maskfile, ou, delta, maxit, CTF, snr, sym, chunk, user_func_name)
		return

	from alignment	    import eqproj
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl
	from fundamentals   import fshift
	from projection     import prep_vol
	from utilities      import amoeba, model_circle, get_arb_params, set_arb_params, dropSpiderDoc
	from utilities      import dropImage
	from math           import pi
	from string         import replace
	from statistics     import fsc_mask
	from fundamentals   import fshift
	import os 
	import sys
	
	from utilities  import print_begin_msg, print_end_msg, print_msg

	import user_functions
	user_func = user_functions.factory[user_func_name]

	print_begin_msg('ali3d_e')
	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Reference volume            : %s\n"%(ref_vol))
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Maskfile                    : %s\n"%(maskfile))

	vol = EMData()
	vol.read_image(ref_vol)
	nx  = vol.get_xsize()
	if (ou <= 0):  ou = nx//2-1

	print_msg("Outer radius                : %i\n"%(ou))
	print_msg("Angular bracket             : %f\n"%(delta))
	print_msg("Maximum iteration           : %i\n"%(maxit))
	print_msg("data with CTF               : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	print_msg("Symmetry group              : %s\n"%(sym))
	print_msg("Chunk of data used          : %-5.2f\n\n"%(chunk))
	
	if os.path.exists(outdir):  os.system('rm -rf '+outdir)
	os.mkdir(outdir)
	
	center = 0

	if CTF :
		parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
		#                  0                  1              2          3              4               5                   6
		#  ERROR if ctf applied
		ima = EMData()
		ima.read_image(stack)
		ctf_params = get_arb_params(ima, parnames)
		if(ctf_params[6] == 1):  ERROR("ali3d_e does not work for CTF-applied data","ali3d_e",1)
		from reconstruction import recons3d_4nn_ctf
	else   : from reconstruction import recons3d_4nn

	nima = EMUtil.get_image_count(stack)

	# figure the size of the chunk (3D is updated after each chunk).  Chunk should be given as 0.0< chunk <= 1.0.  1.0 means all projections
	if(chunk <= 0.0):  chunk = 1.0
	n_in_chunk  = max(int(chunk * nima), 1)
	n_of_chunks = nima//n_in_chunk + min(nima%n_in_chunk,1)
	image_start = 0
	
	if info:
		outf = file(os.path.join(outdir, "progress"), "w")
		outf.write("  chunk = "+str(chunk)+"   ")
		outf.write("\n")
		outf.flush()
		outf.write("  chunk = "+str(n_in_chunk)+"   ")
		outf.write("  chunk = "+str(n_of_chunks)+"   ")
		outf.write("\n")
		outf.flush()
	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask3D=getImage(maskfile)
		else: mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)
	mask2D = model_circle(ou, nx, nx)

	dataim = EMData.read_images(stack)

	# initialize data for the reference preparation function
	ref_data = []
	ref_data.append( mask3D )
	ref_data.append( center )
	ref_data.append( None )
	ref_data.append( None )

	jtep = 0
	par_str=["phi", "theta", "psi", "s2x", "s2y"]
	for iteration in xrange(maxit):
		msg = "ITERATION #%3d\n"%(iteration+1)
		print_msg(msg)
		for  ic  in xrange(n_of_chunks):
			image_start_in_chunk = ic*n_in_chunk
			image_end_in_chunk   = min(image_start_in_chunk + n_in_chunk, nima)
			if info:
				outf.write("image_start_in_chunk "+str(image_start_in_chunk)+"\n")
				outf.write("\n")
				outf.write("image_end_in_chunk "+str(image_end_in_chunk)+"\n")
				outf.write("\n")
				outf.flush()
			jtep += 1
			Util.mul_img(vol, mask3D)
			volft,kb  = prep_vol(vol)
			data = []
			data.append(volft)
			data.append(kb)
			data.append(None)
			data.append(mask2D)
			new_params = []
			for imn in xrange(image_start_in_chunk, image_end_in_chunk):
				"""
				if(imn%50 == 0):
					sys.stdout.write( "\n" )
					sys.stdout.write( " %6d " % imn )
					sys.stdout.flush()
				sys.stdout.write(".")
				sys.stdout.flush()
				"""
				data[2] = dataim[imn-image_start]
				atparams = get_arb_params(dataim[imn-image_start], par_str)
				#  change signs of shifts for projections
				atparams[3] *= -1
				atparams[4] *= -1
				if info:
					initial  = eqproj(atparams, data)  # this is if we need initial discrepancy
					outf.write("Image "+str(imn)+"\n")
					outf.write(' %6.1f  %6.1f  %6.1f  %6.1f  %6.1f   %7.4f '%(atparams[0],atparams[1],atparams[2],atparams[3],atparams[4], initial))
					outf.write("\n")
					outf.flush()
				weight_phi = max(delta, delta*abs((atparams[1]-90.0)/180.0*pi))
				#from utilities import start_time, finish_time
				#t3=start_time()
				optm_params =  amoeba(atparams, [weight_phi,delta,weight_phi,1.0,1.0], eqproj, 1.e-4,1.e-4,500,data)
				optm_params[0].append(imn)
				optm_params[0][3] *= -1
				optm_params[0][4] *= -1
				new_params.append(optm_params[0])
				if info:
					outf.write(' %6.1f  %6.1f  %6.1f  %6.1f  %6.1f   %7.4f    %d4'%(optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4],optm_params[1], optm_params[2]))
					outf.write("\n")
					outf.flush()
				set_arb_params(dataim[imn-image_start], [optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4]], par_str)
				#t4 = finish_time(t3)
			del data
			dropSpiderDoc(os.path.join(outdir, replace("new_params%6d"%(jtep),' ','0')), new_params," phi, theta, psi, s2x, s2y, image number")
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
			if(CTF): vol = recons3d_4nn_ctf(stack, range(nima), snr, 1, sym)
			else:	 vol = recons3d_4nn(stack, range(nima), sym)

			# store the reference volume
			dropImage(vol, os.path.join(outdir, "vol%04d.hdf"%(iteration*n_of_chunks+ic+1)))
			ref_data[2] = vol
			ref_data[3] = fscc
			#  call user-supplied function to prepare reference image, i.e., center and filter it
			vol, cs = user_func( ref_data )

			if center == 1:
				from utilities import rotate_3D_shift
				rotate_3D_shift(dataim, cs)
			dropImage(vol,os.path.join(outdir, "volf%04d.hdf"%(iteration*n_of_chunks+ic+1)))

	#sys.stdout.write( "\n\n" )
	#print  ttime()
	#  here we  write header info
	from utilities import write_headers
	write_headers( stack, data, range(nima))
	print_end_msg("ali3d_e")

def ali3d_e_MPI(stack, ref_vol, outdir, maskfile, ou=-1,  delta=2, maxit=10, CTF = None, snr=1.0, sym="c1", chunk = -1.0, user_func_name="ref_ali3d"):
	"""
		
	"""
	from alignment	    import eqproj
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl, filt_gaussl
	from fundamentals   import fshift, rot_avg_image
	from projection     import prep_vol
	from utilities      import amoeba, bcast_list_to_all, bcast_string_to_all, model_circle, get_arb_params, set_arb_params, dropSpiderDoc
	from utilities      import getImage, dropImage, bcast_EMData_to_all, send_attr_dict, recv_attr_dict
	from utilities      import readSpiderDoc, get_im
	from reconstruction import rec3D_MPI
	from math           import pi
	from string         import replace
	import os
	import sys
	from mpi 	    import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	    import mpi_barrier

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)

	if(CTF):
		from filter import filt_ctf
		parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
	          #             0             1          2          3         4               5              6
	main_node = 0
	if(myid == main_node):
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)
		import user_functions
		user_func = user_functions.factory[user_func_name]
		if  CTF:
			ima = EMData()
			ima.read_image(stack)
			ctf_params = get_arb_params(ima, parnames)
			if(ctf_params[6] == 1):  ERROR("ali3d_e does not work for CTF-applied data","ali3d_e_MPI",1)
	mpi_barrier(MPI_COMM_WORLD)

	nima = EMUtil.get_image_count(stack)

	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)

	# figure the size of the chunk (3D is updated after each chunk).  Chunk should be given as 0.0< chunk <= 1.0.  1.0 means all projections
	if(chunk <= 0.0):  chunk = 1.0
	n_in_chunk  = max(int(chunk * (image_end-image_start+1)), 1)
	n_of_chunks = (image_end-image_start+1)//n_in_chunk + min((image_end-image_start+1)%n_in_chunk,1)

	outf = file(os.path.join(outdir, replace("progress%4d"%myid,' ','0')), "w")
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
	if(ou <= 0):  ou = nx//2-2

	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask3D=getImage(maskfile)
		else:                                   mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)
	mask2D = model_circle(ou, nx, nx)

	#from utilities import readSpiderDoc, set_arb_params
	#prm = readSpiderDoc("params_new.doc")
	#prm_dict = ["phi", "theta", "psi", "s2x", "s2y"]
	dataim = EMData.read_images(stack)

	outf.write("  data read = "+str(image_start)+"   ")
	outf.write("\n")
	outf.flush()


	from  string        import replace

	jtep = 0
	par_str=["phi", "theta", "psi", "s2x", "s2y"]
	for iteration in xrange(maxit):
		outf.write("  iteration = "+str(iteration)+"   ")
		outf.write("\n")
		outf.flush()
		for  ic  in xrange(n_of_chunks):
			jtep += 1
			bcast_EMData_to_all(vol, myid, main_node)
			volft,kb  = prep_vol(vol)
			data = []
			data.append(volft)
			data.append(kb)
			data.append(None)
			data.append(mask2D)

			image_start_in_chunk = image_start + ic*n_in_chunk
			image_end_in_chunk   = min(image_start_in_chunk + n_in_chunk, image_end)
			outf.write("ic "+str(ic)+"   image_start "+str(image_start)+"   n_in_chunk "+str(n_in_chunk)+"   image_end "+str(image_end)+"\n")
			outf.write("image_start_in_chunk "+str(image_start_in_chunk)+"  image_end_in_chunk "+str(image_end_in_chunk)+"\n")
			outf.flush()
			if(CTF):  previous_defocus = -1.0
			for imn in xrange(image_start_in_chunk, image_end_in_chunk):
				if(CTF):
					ctf_params = get_arb_params(dataim[imn-image_start], parnames)
					if(ctf_params[1] != previous_defocus):
						previous_defocus = ctf_params[1]
						data[0],kb = prep_vol(filt_ctf(vol, ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5]))

				data[2] = dataim[imn-image_start]

				atparams = get_arb_params(dataim[imn-image_start], par_str)
				#  change signs of shifts for projections
				atparams[3] *= -1
				atparams[4] *= -1

				#optm_params = ali_G3(data, atparams, dtheta)
				#  Align only Euler angles
				#  change signs of shifts for projections
				#data.insert(3, -atparams[3])
				#data.insert(4, -atparams[4])

				#initial  = eqproj(atparams, data)  # this is if we need initial discrepancy
				outf.write("Image "+str(imn)+"\n")
				outf.write('Old  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f '%(atparams[0], atparams[1], atparams[2], -atparams[3], -atparams[4]))
				outf.write("\n")
				#outf.flush()

				weight_phi = max(delta, delta*abs((atparams[1]-90.0)/180.0*pi))
				#from utilities import start_time, finish_time
				#t3=start_time()
				optm_params =  amoeba(atparams, [weight_phi, delta, weight_phi, 1.0, 1.0], eqproj, 1.e-4, 1.e-4,500, data)
				optm_params[0].append(imn)
				optm_params[0][3] *= -1
				optm_params[0][4] *= -1

				outf.write('New  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f     %7.4f    %d4   %7.1f'%(optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4],optm_params[1], optm_params[2], ctf_params[1]))
				outf.write("\n")
				outf.flush()

				set_arb_params(dataim[imn-image_start], optm_params[0], par_str)

			del data
			soto = []
			for imn in xrange(image_start, image_end):
				s2x   = dataim[imn-image_start].get_attr('s2x')
				s2y   = dataim[imn-image_start].get_attr('s2y')
				phi   = dataim[imn-image_start].get_attr('phi')
				theta = dataim[imn-image_start].get_attr('theta')
				psi   = dataim[imn-image_start].get_attr('psi')
				soto.append([phi,theta,psi,s2x,s2y,imn])
			dropSpiderDoc(os.path.join(outdir, replace("new_params%3d_%3d"%(iteration, ic),' ','0')), soto," phi, theta, psi, s2x, s2y, image number")
			del soto

			# compute updated 3D after each chunk
 	    		# resolution
			outf.write("  begin reconstruction = "+str(image_start)+"   ")
			outf.write("\n")
			outf.flush()
			vol, fscc = rec3D_MPI(dataim, snr, sym, mask3D, os.path.join(outdir, replace("resolution%3d_%3d"%(iteration, ic),' ','0') ), myid, main_node)
			outf.write("  done reconstruction = "+str(image_start)+"   ")
			outf.write("\n")
			outf.flush()
			if(myid == main_node):
				dropImage(vol, os.path.join(outdir, replace("vol%3d_%3d.hdf"%(iteration, ic),' ','0') ))
				center = 1
				ref_data = []
				ref_data.append( mask3D )
				ref_data.append( center )
				ref_data.append( vol )
				ref_data.append( fscc )
				#  call user-supplied function to prepare reference image, i.e., filter it
				vol, cs = user_func( ref_data )
				if center == 1:
					from utilities import rotate_3D_shift
					rotate_3D_shift(dataim, cs)
				dropImage(vol, os.path.join(outdir, replace("volf%3d_%3d.hdf"%(iteration, ic),' ','0') ))

	del vol
	del volft
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
     
	if(myid == main_node): recv_attr_dict(main_node, stack, dataim, par_str, image_start, image_end, number_of_proc)
	else: send_attr_dict(main_node, dataim, par_str, image_start, image_end)

def eqprojG3(args, data):
	from utilities import peak_search, amoeba, dropImage, info, pad
	from fundamentals import fft, ccf
	from sys import exit
	from development import twoD_fine_search

	volft = data[0]
	kb = data[1]
	kb2 = data[6]
	params = args
	
	# This part is copied from prgs
	EULER_SPIDER = Transform3D.EulerType.SPIDER
	phi = params[0]
	theta = params[1]
	psi = params[2]
	R= Transform3D(EULER_SPIDER,phi,theta,psi)
	temp = volft.extractplane(R,kb)
	M = temp.get_ysize()	
	temp = temp.Four_shuf_ds_cen_us(M, M, 1, False)

	nx = M/2
	sx = (nx-data[7][0]*2)/2.0
	sy = (nx-data[7][1]*2)/2.0

	product = ccf(temp, data[5])
	data2 = [0]*2
	data2[0] = product
	data2[1] = kb2
	ps = amoeba([sx,sy],[1.0,1.0],twoD_fine_search,1.e-5,1.e-5,500,data2)
	
	s2x = (nx-ps[0][0]*2)/2
	s2y = (nx-ps[0][1]*2)/2
	
	params2 = {"filter_type" : Processor.fourier_filter_types.SHIFT, "x_shift" : s2x*2, "y_shift" : s2y*2, "z_shift" : 0.0}
	temp2 = Processor.EMFourierFilter(temp, params2)
	v = -temp2.cmp("SqEuclidean", data[4])

	return v, [s2x, s2y]

def prepij(image):
	M=image.get_ysize()
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

def ali3d_eB_ORIGINAL(stack, ref_vol, outdir, maskfile = None, ou = -1,  delta = 2, maxit = 10, CTF = False, snr = 1.0, sym="c1", chunk = -1.0, MPI=False):
	"""
		
	"""
	if MPI:
		ali3d_eB_MPI(stack, ref_vol, outdir, maskfile, ou, delta, maxit, CTF, snr, sym, chunk)
		return

	from alignment	    import eqproj
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl
	from fundamentals   import fshift
	from projection     import prep_vol
	from utilities      import amoeba2, model_circle, get_arb_params, set_arb_params, dropSpiderDoc
	from utilities      import dropImage
	from math           import pi,sin
	from string         import replace
	from statistics     import fsc_mask
	import os 
	import sys
	
	from utilities  import print_begin_msg, print_end_msg, print_msg
	print_begin_msg('ali3d_e')
	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Reference volume            : %s\n"%(ref_vol))
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Maskfile                    : %s\n"%(maskfile))

	vol = EMData()
	vol.read_image(ref_vol)
	nx  = vol.get_xsize()
	if (ou <= 0):  ou = nx//2-1

	print_msg("Outer radius                : %i\n"%(ou))
	print_msg("Angular bracket             : %f\n"%(delta))
	print_msg("Maximum iteration           : %i\n"%(maxit))
	print_msg("data with CTF               : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	print_msg("Symmetry group              : %s\n"%(sym))
	print_msg("Chunk of data used          : %-5.2f\n\n"%(chunk))
	
	if os.path.exists(outdir):  os.system('rm -rf '+outdir)
	os.mkdir(outdir)

	if CTF :from reconstruction import recons3d_4nn_ctf
	else   : from reconstruction import recons3d_4nn

	nima = EMUtil.get_image_count(stack)

	# figure the size of the chunk (3D is updated after each chunk).  Chunk should be given as 0.0< chunk <= 1.0.  1.0 means all projections
	if(chunk <= 0.0):  chunk = 1.0
	n_in_chunk  = max(int(chunk * nima), 1)
	n_of_chunks = nima//n_in_chunk + min(nima%n_in_chunk,1)
	image_start = 0
	
	outf = file(os.path.join(outdir, "progress"), "w")
	outf.write("  chunk = "+str(chunk)+"   ")
	outf.write("\n")
	outf.flush()
	outf.write("  chunk = "+str(n_in_chunk)+"   ")
	outf.write("  chunk = "+str(n_of_chunks)+"   ")
	outf.write("\n")
	outf.flush()
	vol.write_image(os.path.join(outdir,"ref_volf00.hdf"))
	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask3D=getImage(maskfile)
		else: mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)
	mask2D = model_circle(ou, nx, nx)

	parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
	#                  0                  1              2          3              4               5                   6
	dataim = EMData.read_images(stack)
	nima = len(dataim)
	for im in xrange(nima):
		dataim[im].set_attr('ID', im)
		if(CTF):
			ctf_params = get_arb_params(data[im], parnames)
			if(im == 0): data_had_ctf = ctf_params[6]
			if(ctf_params[6] == 0):
				st = Util.infomask(data[im], mask2D, False)
				data[im] -= st[0]
				from filter import filt_ctf
				data[im] = filt_ctf(data[im], ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
				data[im].set_attr('ctf_applied', 1)

	jtep = 0
	par_str=["phi", "theta", "psi", "s2x", "s2y"]
	for iteration in xrange(maxit):
		msg = "ITERATION #%3d\n"%(iteration+1)
		print_msg(msg)
		for  ic  in xrange(n_of_chunks):
			image_start_in_chunk = ic*n_in_chunk
			image_end_in_chunk   = min(image_start_in_chunk + n_in_chunk, nima)
			outf.write("image_start_in_chunk "+str(image_start_in_chunk)+"\n")
			outf.write("\n")
			outf.write("image_end_in_chunk "+str(image_end_in_chunk)+"\n")
			outf.write("\n")
			outf.flush()
			jtep += 1
			Util.mul_img(vol, mask3D)
			volft,kb  = prep_vol(vol)
			data = [0]*8
			data[0] = volft
			data[1] = kb
			data[3] = mask2D
			new_params = []
			for imn in xrange(image_start_in_chunk, image_end_in_chunk):
				"""
				if(imn%50 == 0):
					sys.stdout.write( "\n" )
					sys.stdout.write( " %6d " % imn )
					sys.stdout.flush()
				sys.stdout.write(".")
				sys.stdout.flush()
				"""
				data[2] = dataim[imn-image_start].copy()
			
				refi = dataim[imn-image_start].copy()
				oo, qq, kb2 = prepij(refi)
				data[4] = oo
				data[5] = qq
				data[6] = kb2 
				
				atparams = get_arb_params(dataim[imn-image_start], par_str)
				atparams[3] *= -1
				atparams[4] *= -1

				data[7] = [atparams[3], atparams[4]]
				del atparams[3]
				del atparams[3]
				weight_phi = max(delta, delta*abs((atparams[1]-90.0)/180.0*pi))			
			
				# For downhill simplex method 
				optm_params =  amoeba2(atparams, [weight_phi, delta, weight_phi], eqprojG3, 1.e-5,1.e-5,500,data)
				optm_params[3][0] *= -1
				optm_params[3][1] *= -1
				set_arb_params(dataim[imn-image_start], [optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[3][0], optm_params[3][1]], par_str)
			
			dropSpiderDoc(os.path.join(outdir, replace("new_params%6d"%(jtep),' ','0')), new_params," phi, theta, psi, s2x, s2y, image number")
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
			#dropImage(vol,os.path.join(outdir, replace("vol%4d.spi"%(N_step*max_iter+Iter+1),' ','0')), "s")
			dropImage(vol,os.path.join(outdir, replace("vol%4d.spi"%(iteration*n_of_chunks+ic+1),' ','0')), "s")
			#filt = filt_from_fsc(fscc, 0.05)
			#vol  = filt_table(vol, filt)
			# here figure the filtration parameters and filter vol for the  next iteration
			#fl, fh = filt_params(res)
			#vol    = filt_btwl(vol, fl, fh)
			# store the filtred reference volume
			lk = 0
			while(fscc[1][lk] >0.9 and fscc[0][lk]<0.25):
				lk+=1
			fl = fscc[0][lk]
			fh = min(fl+0.1,0.49)
			vol = filt_btwl(vol, fl, fh)
			cs   = vol.phase_cog()
			vol  = fshift(vol, -cs[0], -cs[1] -cs[2])
			#dropImage(vol,os.path.join(outdir, replace("volf%4d.spi"%(N_step*max_iter+Iter+1),' ','0')), "s")
			dropImage(vol,os.path.join(outdir, replace("volf%4d.spi"%(iteration*n_of_chunks+ic+1),' ','0')), "s")

		#sys.stdout.write( "\n\n" )
		#print  ttime()
		#  here we  write header info
	if(CTF and data_had_ctf == 0):
		for im in xrange(nima): data[im].set_attr('ctf_applied', 0)
	from utilities import write_headers
	write_headers( stack, data, range(nima))
	print_end_msg('ali3d_e')	    	    

def ali3d_eB(stack, ref_vol, outdir, maskfile = None, ou = -1,  delta = 2, maxit = 10, CTF = False, snr = 1.0, sym="c1", chunk = -1.0, user_func_name="ref_aliB_cone", MPI=False):
	"""
		This is modified MPI version to test CCC
		with local searches
	"""
	if MPI:
		ali3d_eB_MPI(stack, ref_vol, outdir, maskfile, ou, delta, maxit, CTF, snr, sym, chunk, user_func_name)
		return

	from alignment	    import eqproj, eqprojEuler
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl
	from fundamentals   import fshift
	from projection     import prep_vol
	from utilities      import amoeba, model_circle, get_arb_params, set_arb_params, dropSpiderDoc
	from utilities      import dropImage, get_im
	from math           import pi,sin
	from string         import replace
	from statistics     import fsc_mask
	import os
	import sys
	
	from utilities  import print_begin_msg, print_end_msg, print_msg
	import user_functions
	user_func = user_functions.factory[user_func_name]
	print_begin_msg('ali3d_e')
	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Reference volume            : %s\n"%(ref_vol))
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Maskfile                    : %s\n"%(maskfile))

	vol = EMData()
	vol.read_image(ref_vol)
	nx  = vol.get_xsize()
	if (ou <= 0):  ou = nx//2-1

	print_msg("Outer radius                : %i\n"%(ou))
	print_msg("Angular bracket             : %f\n"%(delta))
	print_msg("Maximum iteration           : %i\n"%(maxit))
	print_msg("data with CTF               : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	print_msg("Symmetry group              : %s\n"%(sym))
	print_msg("Chunk of data used          : %-5.2f\n"%(chunk))
	print_msg("User function               : %-s\n\n"%(user_func_name))
	
	if os.path.exists(outdir):  os.system('rm -rf '+outdir)
	os.mkdir(outdir)

	if CTF :from reconstruction import recons3d_4nn_ctf
	else   : from reconstruction import recons3d_4nn

	nima = EMUtil.get_image_count(stack)
	image_start = 0
	image_end   = nima

	# figure the size of the chunk (3D is updated after each chunk).  Chunk should be given as 0.0< chunk <= 1.0.  1.0 means all projections
	if(chunk <= 0.0):  chunk = 1.0
	n_in_chunk  = max(int(chunk * nima), 1)
	n_of_chunks = nima//n_in_chunk + min(nima%n_in_chunk,1)
	image_start = 0
	
	outf = file(os.path.join(outdir, "progress"), "w")
	outf.write("  chunk = "+str(chunk)+"   ")
	outf.write("\n")
	outf.flush()
	outf.write("  chunk = "+str(n_in_chunk)+"   ")
	outf.write("  chunk = "+str(n_of_chunks)+"   ")
	outf.write("\n")
	outf.flush()
	vol.write_image(os.path.join(outdir,"ref_volf00.hdf"))
	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask3D=get_im(maskfile)
		else: mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)
	mask2D = model_circle(ou, nx, nx)

	parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
	#                  0                  1              2          3              4               5                   6
	dataim = EMData.read_images(stack)
	nima = len(dataim)
	for im in xrange(nima):
		dataim[im].set_attr('ID', im)

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
	for iteration in xrange(maxit):
		msg = "ITERATION #%3d\n"%(iteration+1)
		print_msg(msg)
		for  ic  in xrange(n_of_chunks):
			jtep += 1
			dropImage(vol, os.path.join(outdir, replace("vol%3d_%3d.hdf"%(iteration, ic),' ','0') ))
			ref_data.append( vol )
			ref_data.append( fscc )
			#  call user-supplied function to prepare reference image, i.e., filter it
			vol = user_func( ref_data )
			#  HERE CS SHOULD BE USED TO MODIFY PROJECTIONS' PARAMETERS  !!!
			del ref_data[2]
			del ref_data[2]
			dropImage(vol, os.path.join(outdir, replace("volf%3d_%3d.hdf"%(iteration, ic),' ','0') ))

			volft,kb  = prep_vol(vol)
			data = []
			data.append(volft)
			data.append(kb)
			data.append(None)
			data.append(None)
			data.append(None)
			data.append(mask2D)

			image_start_in_chunk = image_start + ic*n_in_chunk
			image_end_in_chunk   = min(image_start_in_chunk + n_in_chunk, image_end)
			outf.write("ic "+str(ic)+"   image_start "+str(image_start)+"   n_in_chunk "+str(n_in_chunk)+"   image_end "+str(image_end)+"\n")
			outf.write("image_start_in_chunk "+str(image_start_in_chunk)+"  image_end_in_chunk "+str(image_end_in_chunk)+"\n")
			outf.flush()
			if(CTF):  previous_defocus = -1.0
			for imn in xrange(image_start_in_chunk, image_end_in_chunk):
			#for imn in xrange(1):
				if(CTF):
					ctf_params = get_arb_params(dataim[imn-image_start], parnames)
					if(ctf_params[6] == 0 and (ctf_params[1] != previous_defocus) ):
						previous_defocus = ctf_params[1]
						data[0],kb = prep_vol(filt_ctf(vol, ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5]))
				data[2] = dataim[imn-image_start]

				atparams = get_arb_params(dataim[imn-image_start], par_str)
				#atparams[0]=  240.083  #    112.392 #+3.5#110.560  #
				#atparams[1]=  66.828 #154.089 #154.728 #
				#atparams[2]=  308.566#278.59  #+3.5#  276.900#
				weight_phi = max(delta, delta*abs((atparams[1]-90.0)/180.0*pi))
			
				#optm_params = ali_G3(data, atparams, dtheta)
				#  Align only Euler angles
				#  change signs of shifts for projections
				data[3] = -atparams[3]
				data[4] = -atparams[4]
				#outf.write("Image "+str(imn)+"\n")
				#outf.write('Old %6.1f  %6.1f  %6.1f  %6.1f  %6.1f   %7.4f '%(atparams[0],atparams[1],atparams[2],-atparams[3],-atparams[4], initial))
				#outf.write('Old %8.3f  %8.3f  %8.3f  '%(atparams[0],atparams[1],atparams[2]))
				#outf.write("\n")
				initial  = eqprojEuler([atparams[0], atparams[1], atparams[2]], data)
				outf.write(' %11.7f    %8.3f  %8.3f  %8.3f\n'%(initial, atparams[0], atparams[1], atparams[2]))
				'''
				#  generate 3D ccf...
				ist = 25
				rng = 3.2  # +/- search range
				isr = float(ist//2)
				ict = isr+1
				angt = [0.0]*3
				for  iphi in xrange(ist):
					angt[0] = rng*(iphi - ict)/isr + atparams[0]
					for  itheta in xrange(ist):  # -0.9:0.9
						angt[1] = rng*(itheta - ict)/isr + atparams[1]
						for  ipsi in xrange(ist):
							angt[2] = rng*(ipsi - ict)/isr + atparams[2]
							initial  = eqprojEuler(angt, data)  # this is if we need initial discrepancy
							outf.write(' %11.7f    %8.3f  %8.3f  %8.3f\n'%(initial, angt[0],angt[1],angt[2]))
					outf.flush()
				'''

				#from utilities import start_time, finish_time
				#t3=start_time()
				optm_params =  amoeba(atparams[0:3], [weight_phi, delta, weight_phi], eqprojEuler, 1.e-4,1.e-4,500, data)
				#optm_params[0].append(imn)
				print  optm_params[2]
				#new_params.append(optm_params[0])
				#outf.write('New %6.1f  %6.1f  %6.1f  %6.1f  %6.1f   %7.4f    %d4   %7.1f'%(optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4],optm_params[1], optm_params[2], ctf_params[1]))
				#outf.write('New %8.3f  %8.3f  %8.3f  %11.7f  %6.1f   '%(optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[1], optm_params[2]))
				#outf.write("\n")
				initial  = eqprojEuler([optm_params[0][0], optm_params[0][1], optm_params[0][2]], data)

				outf.write(' %11.7f    %8.3f  %8.3f  %8.3f\n'%(initial, optm_params[0][0], optm_params[0][1], optm_params[0][2]))
				outf.flush()
				#set_arb_params(dataim[imn-image_start], [optm_params[0][0], optm_params[0][1], optm_params[0][2]], par_str[0:3])
			
			from sys import exit
			exit()
			del data
			soto = []
			for imn in xrange(image_start, image_end):
				s2x   = dataim[imn-image_start].get_attr('s2x')
				s2y   = dataim[imn-image_start].get_attr('s2y')
				phi   = dataim[imn-image_start].get_attr('phi')
				theta = dataim[imn-image_start].get_attr('theta')
				psi   = dataim[imn-image_start].get_attr('psi')
				soto.append([phi,theta,psi,s2x,s2y,imn])
			dropSpiderDoc(os.path.join(outdir, replace("new_params%3d_%3d"%(iteration, ic),' ','0')), soto," phi, theta, psi, s2x, s2y, image number")
			del soto
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
			#dropImage(vol,os.path.join(outdir, replace("vol%4d.spi"%(N_step*max_iter+Iter+1),' ','0')), "s")
			dropImage(vol,os.path.join(outdir, replace("vol%4d.hdf"%(iteration*n_of_chunks+ic+1),' ','0')), "s")
	'''
		#sys.stdout.write( "\n\n" )
		#print  ttime()
		#  here we  write header info
	if(CTF and data_had_ctf == 0):
		for im in xrange(nima): data[im].set_attr('ctf_applied', 0)
	for im in xrange(nima):
		dataim[im].write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)			    	    					    	    
	print_end_msg('ali3d_e')
	'''

def ali3d_eB_MPI_LAST_USED(stack, ref_vol, outdir, maskfile, ou=-1,  delta=2, maxit=10, CTF = None, snr=1.0, sym="c1", chunk = -1.0, user_func_name="ref_aliB_cone"):
	"""
		Cone
	"""
	from alignment	    import proj_ali_incore_cone
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl, filt_tanl
	from fundamentals   import fshift, rot_avg_image
	from projection     import prep_vol, prgs
	from utilities      import bcast_list_to_all, bcast_string_to_all, model_circle, get_arb_params, set_arb_params, dropSpiderDoc
	from utilities      import getImage, dropImage, bcast_EMData_to_all, send_attr_dict, recv_attr_dict
	from utilities      import readSpiderDoc, get_im
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
		if(type(maskfile) is types.StringType):  mask3D = getImage(maskfile)
		else:                                   mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)
	mask2D = model_circle(ou, nx, nx)

	parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
	#                  0                 1            2            3            4                 5                 6
	#from utilities import readSpiderDoc, set_arb_params
	#prm = readSpiderDoc("params_new.doc")
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
			ctf_params = get_arb_params(data[im], parnames)
			if(im == image_start): data_had_ctf = ctf_params[6]
			if(ctf_params[6] == 0):
				st = Util.infomask(data[im], mask2D, False)
				data[im] -= st[0]
				from filter import filt_ctf
				data[im] = filt_ctf(data[im], ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
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
					dropImage(vol, os.path.join(outdir, replace("vol%3d_%3d.hdf"%(iteration, ic),' ','0') ))
					
					ref_data.append( vol )
					ref_data.append( fscc )
					#  call user-supplied function to prepare reference image, i.e., filter it
					vol = user_func( ref_data )
					#  HERE CS SHOULD BE USED TO MODIFY PROJECTIONS' PARAMETERS  !!!
					del ref_data[2]
					del ref_data[2]
					dropImage(vol, os.path.join(outdir, replace("volf%3d_%3d.hdf"%(iteration, ic),' ','0') ))
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
				s2x   = dataim[imn-image_start].get_attr('s2x')
				s2y   = dataim[imn-image_start].get_attr('s2y')
				phi   = dataim[imn-image_start].get_attr('phi')
				theta = dataim[imn-image_start].get_attr('theta')
				psi   = dataim[imn-image_start].get_attr('psi')
				soto.append([phi,theta,psi,s2x,s2y,imn])
			dropSpiderDoc(os.path.join(outdir, replace("new_params%3d_%3d_%3d"%(iteration, ic, myid),' ','0')), soto," phi, theta, psi, s2x, s2y, image number")
			del soto

			#from sys import exit
			#exit()
		#  here we should write header info, just in case the program crashes...
	del vol
	del volft
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)

def ali3d_eB_CCC(stack, ref_vol, outdir, maskfile, ou=-1,  delta=2, maxit=10, CTF = None, snr=1.0, sym="c1", chunk = -1.0, user_func_name="ref_aliB_cone"):
	"""
		Cone, modified version to test CCC
		single processor version
	"""
	from utilities      import print_begin_msg, print_end_msg, print_msg

	from alignment	    import proj_ali_incore_cone
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl, filt_tanl
	from fundamentals   import fshift, rot_avg_image
	from projection     import prep_vol, prgs
	from utilities      import model_circle, get_arb_params, set_arb_params, dropSpiderDoc
	from utilities      import getImage, dropImage, send_attr_dict, recv_attr_dict
	from utilities      import readSpiderDoc, get_im
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
		if(type(maskfile) is types.StringType):  mask3D = getImage(maskfile)
		else:                                   mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)
	mask2D = model_circle(ou, nx, nx)

	dataim = []
	parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
	#                  0                 1            2            3            4                 5                 6
	#from utilities import readSpiderDoc, set_arb_params
	#prm = readSpiderDoc("params_new.doc")
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
			ctf_params = get_arb_params(ima, parnames)
			if(im == image_start): data_had_ctf = ctf_params[6]
			if(ctf_params[6] == 0):
				st = Util.infomask(ima, mask2D, False)
				ima -= st[0]
				from filter import filt_ctf
				ima = filt_ctf(ima, ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
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
			dropImage(vol, os.path.join(outdir, replace("vol%3d_%3d.hdf"%(iteration, ic),' ','0') ))
			ref_data.append( vol )
			ref_data.append( fscc )
			#  call user-supplied function to prepare reference image, i.e., filter it
			vol = user_func( ref_data )
			#  HERE CS SHOULD BE USED TO MODIFY PROJECTIONS' PARAMETERS  !!!
			del ref_data[2]
			del ref_data[2]
			dropImage(vol, os.path.join(outdir, replace("volf%3d_%3d.hdf"%(iteration, ic),' ','0') ))

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
				s2x   = dataim[imn-image_start].get_attr('s2x')
				s2y   = dataim[imn-image_start].get_attr('s2y')
				phi   = dataim[imn-image_start].get_attr('phi')
				theta = dataim[imn-image_start].get_attr('theta')
				psi   = dataim[imn-image_start].get_attr('psi')
				soto.append([phi,theta,psi,s2x,s2y,imn])
			dropSpiderDoc(os.path.join(outdir, replace("new_params%3d_%3d"%(iteration, ic),' ','0')), soto," phi, theta, psi, s2x, s2y, image number")
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
			#dropImage(vol,os.path.join(outdir, replace("vol%4d.spi"%(N_step*max_iter+Iter+1),' ','0')), "s")
			dropImage(vol,os.path.join(outdir, replace("vol%4d.hdf"%(iteration*n_of_chunks+ic+1),' ','0')), "s")

def ali3d_eB_MPI_conewithselect(stack, ref_vol, outdir, maskfile, ou=-1,  delta=2, maxit=10, CTF = None, snr=1.0, sym="c1", chunk = -1.0):
	"""
		Cone
	"""
	from alignment	    import proj_ali_incore_cone
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl, filt_tanl
	from fundamentals   import fshift, rot_avg_image
	from projection     import prep_vol, prgs
	from utilities      import bcast_list_to_all, bcast_string_to_all, model_circle, get_arb_params, set_arb_params, dropSpiderDoc
	from utilities      import getImage, dropImage, bcast_EMData_to_all, send_attr_dict, recv_attr_dict
	from utilities      import readSpiderDoc, get_im
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
		if(type(maskfile) is types.StringType):  mask3D = getImage(maskfile)
		else:                                   mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)
	mask2D = model_circle(ou, nx, nx)

	dataim = []
	parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
	#                  0                 1            2            3            4                 5                 6
	#from utilities import readSpiderDoc, set_arb_params
	#prm = readSpiderDoc("params_new.doc")
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
			ctf_params = get_arb_params(ima, parnames)
			if(im == image_start): data_had_ctf = ctf_params[6]
			if(ctf_params[6] == 0):
				st = Util.infomask(ima, mask2D, False)
				ima -= st[0]
				from filter import filt_ctf
				ima = filt_ctf(ima, ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
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
				"""
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
				"""
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
				dropImage(vol, os.path.join(outdir, replace("vol%3d_%3d.hdf"%(iteration, ic),' ','0') ))
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
				dropImage(vol, os.path.join(outdir, replace("volf%3d_%3d.hdf"%(iteration, ic),' ','0') ))
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
				s2x   = dataim[imn-image_start].get_attr('s2x')
				s2y   = dataim[imn-image_start].get_attr('s2y')
				phi   = dataim[imn-image_start].get_attr('phi')
				theta = dataim[imn-image_start].get_attr('theta')
				psi   = dataim[imn-image_start].get_attr('psi')
				soto.append([phi,theta,psi,s2x,s2y,imn])
			dropSpiderDoc(os.path.join(outdir, replace("new_params%3d_%3d_%3d"%(iteration, ic, myid),' ','0')), soto," phi, theta, psi, s2x, s2y, image number")
			del soto

			#from sys import exit
			#exit()
		#  here we should write header info, just in case the program crashes...
	del vol
	del volft
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)


def ali3d_eB_MPI(stack, ref_vol, outdir, maskfile, ou=-1,  delta=2, maxit=10, CTF = None, snr=1.0, sym="c1", chunk = -1.0, user_func_name="ref_aliB_cone"):
	"""
		Version 03/20/08, all particles used
	"""
	from alignment	    import eqprojEuler
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl, filt_tanl
	from fundamentals   import fshift, rot_avg_image
	from projection     import prep_vol, prgs
	from utilities      import amoeba, bcast_list_to_all, bcast_string_to_all, model_circle, get_arb_params, set_arb_params, dropSpiderDoc
	from utilities      import getImage, dropImage, bcast_EMData_to_all, send_attr_dict, recv_attr_dict
	from utilities      import readSpiderDoc, get_im
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

	parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
	#                  0                 1            2            3            4                 5                 6
	if(myid == main_node):
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)
		import user_functions
		user_func = user_functions.factory[user_func_name]
		#  ERROR if ctf applied
		ima = EMData()
		ima.read_image(stack)
		ctf_params = get_arb_params(ima, parnames)
		if(ctf_params[6] == 1):
			ERROR("ali3d_e does not work for CTF-applied data","ali3d_eB_MPI",1)
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
		if(type(maskfile) is types.StringType):  mask3D = getImage(maskfile)
		else:                                   mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)
	mask2D = model_circle(ou, nx, nx)

	dataim = []
	#from utilities import readSpiderDoc, set_arb_params
	#prm = readSpiderDoc("params_new.doc")
	#prm_dict = ["phi", "theta", "psi", "s2x", "s2y"]
	prm_dict = ["phi", "theta", "psi"]
	#from random import seed,gauss
	#seed()
	dataim = EMData.read_images(image_start, image_end)
	for im in xrange(image_start, image_end):
		dataim[im].set_attr('ID', im)
		#angn = get_arb_params(ima, prm_dict)
		#for ian in xrange(3):  angn[ian] += gauss(0.0, 1.0)
		#set_arb_params(ima, angn, prm_dict)
		#set_arb_params(ima, prm[im], prm_dict)
		'''
		if(CTF):
			ctf_params = get_arb_params(ima, parnames)
			if(im == image_start): data_had_ctf = ctf_params[6]
			if(ctf_params[6] == 0):
				st = Util.infomask(ima, mask2D, False)
				ima -= st[0]
				from filter import filt_ctf
				ima = filt_ctf(ima, ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
				ima.set_attr('ctf_applied', 1)
		'''
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
	'''
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(nima, number_of_proc, im)
		recvcount.append( ie - ib )
	'''
	#from filter import  filt_tophatb, filt_gaussl
	#volftb = vol.copy()

	for iteration in xrange(maxit):
		outf.write("  iteration = "+str(iteration)+"   ")
		outf.write("\n")
		outf.flush()
		for  ic  in xrange(n_of_chunks):
			# compute updated 3D after each chunk
			'''
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
				"""
				for im in xrange(len(templ)):
					outf.write('ccc, image %12.5f  %07d'%( templ[im][0], templ[im][1]  ))
					outf.write("\n")
				outf.flush()
				"""
				
				templ.sort()
				ilow = int(0.25*len(templ))  # reject 25% worst images.
				for im in xrange(ilow):  templ[im] = [templ[im][1], 0]
				for im in xrange(ilow, len(templ)):  templ[im] = [templ[im][1], 1]
				templ.sort()
				sendbuf = []
				for im in xrange(len(templ)):	sendbuf.append(templ[im][1])
				del templ
				"""
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
				"""
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
			"""
			nact = 0
			for imn in xrange(image_start, image_end):
				nact += dataim[imn-image_start].get_attr('active')
			nact = float(nact)
			tn = mpi_reduce(nact, 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
			if(myid == main_node):
				outf.write('total number of used images %12.2f  '%(float(tn)))
				outf.write("\n")
				outf.flush()
			"""
			# resolution
			outf.write("  begin reconstruction = "+str(image_start)+"   ")
			outf.write("\n")
			outf.flush()
			vol, fscc = rec3D_MPI(dataim, snr, sym, mask3D, os.path.join(outdir, replace("resolution%3d_%3d"%(iteration, ic),' ','0') ), myid, main_node)
			outf.write("  done reconstruction = "+str(image_start)+"   ")
			outf.write("\n")
			outf.flush()
			#volftb = vol.copy()

			#  restore original normalization
			#for imn in xrange(image_start, image_end):
			#	#dataim[imn-image_start].set_attr_dict({'active': int(recvbuf[imn-image_start])})
			#	dataim[imn-image_start] *= float(recvbuf[imn-image_start])
			#del recvbuf
			mpi_barrier(MPI_COMM_WORLD)
			if(myid == main_node):
				dropImage(vol, os.path.join(outdir, replace("vol%3d_%3d.hdf"%(iteration, ic),' ','0') ))
				ref_data.append( vol )
				ref_data.append( fscc )
				#  call user-supplied function to prepare reference image, i.e., filter it
				vol = user_func( ref_data )
				#  HERE CS SHOULD BE USED TO MODIFY PROJECTIONS' PARAMETERS  !!!
				del ref_data[2]
				del ref_data[2]
				dropImage(vol, os.path.join(outdir, replace("volf%3d_%3d.hdf"%(iteration, ic),' ','0') ))
			from sys import exit
			exit()

			#bcast_EMData_to_all(volftb, myid, main_node)
			bcast_EMData_to_all(vol, myid, main_node)
			#volft,kb  = prep_vol(vol)
			data = []
			# Only Euler angles
			data.append(None)
			data.append(None)
			data.append(None)
			data.append(None)
			data.append(None)
			data.append(mask2D)

			image_start_in_chunk = image_start + ic*n_in_chunk
			image_end_in_chunk   = min(image_start_in_chunk + n_in_chunk, image_end)
			outf.write("ic "+str(ic)+"   image_start "+str(image_start)+"   n_in_chunk "+str(n_in_chunk)+"   image_end "+str(image_end)+"\n")
			outf.write("image_start_in_chunk "+str(image_start_in_chunk)+"  image_end_in_chunk "+str(image_end_in_chunk)+"\n")
			outf.flush()
			if(CTF):  previous_defocus = -1.0
			for imn in xrange(image_start_in_chunk, image_end_in_chunk):
				if(CTF):
					ctf_params = get_arb_params(dataim[imn-image_start], parnames)
					if(ctf_params[1] != previous_defocus):
						previous_defocus = ctf_params[1]
						data[0],data[1] = prep_vol(filt_ctf(vol, ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5]))
				data[2] = dataim[imn-image_start]

				atparams = get_arb_params(dataim[imn-image_start], par_str)
				
				#optm_params = ali_G3(data, atparams, dtheta)
				#  Align only Euler angles
				#  change signs of shifts for projections
				data[3] = -atparams[3]
				data[4] = -atparams[4]

				#initial  = eqproj(atparams, data)  # this is if we need initial discrepancy
				outf.write("Image "+str(imn)+"\n")
				#outf.write('Old %6.1f  %6.1f  %6.1f  %6.1f  %6.1f   %7.4f '%(atparams[0],atparams[1],atparams[2],-atparams[3],-atparams[4], initial))
				outf.write('Old %8.3f  %8.3f  %8.3f  '%(atparams[0],atparams[1],atparams[2]))
				outf.write("\n")
				#outf.flush()

				weight_phi = max(delta, delta*abs((atparams[1]-90.0)/180.0*pi))
				#from utilities import start_time, finish_time
				#t3=start_time()
				optm_params =  amoeba(atparams[0:3], [weight_phi, delta, weight_phi], eqprojEuler, 1.e-4,1.e-4,500, data)
				optm_params[0].append(imn)
				#new_params.append(optm_params[0])

				#outf.write('New %6.1f  %6.1f  %6.1f  %6.1f  %6.1f   %7.4f    %d4   %7.1f'%(optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4],optm_params[1], optm_params[2], ctf_params[1]))
				outf.write('New %8.3f  %8.3f  %8.3f  %8.5f  %6.1f   '%(optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[1], optm_params[2]))
				outf.write("\n")
				outf.flush()

				set_arb_params(dataim[imn-image_start], [optm_params[0][0], optm_params[0][1], optm_params[0][2]], par_str[0:3])

				"""

				#  change signs of shifts for projections
				#atparams[3] *= -1
				#atparams[4] *= -1

				initial  = eqproj(atparams, data)  # this is if we need initial discrepancy
				outf.write("Image "+str(imn)+"\n")
				outf.write('Old %6.1f  %6.1f  %6.1f  %6.1f  %6.1f   %7.4f '%(atparams[0],atparams[1],atparams[2],-atparams[3],-atparams[4], initial))
				outf.write("\n")
				outf.flush()

				#weight_phi = max(delta, delta*abs((atparams[1]-90.0)/180.0*pi))
				#from utilities import start_time, finish_time
				#t3=start_time()
				#optm_params =  amoeba(atparams, [weight_phi, delta, weight_phi, 0.05, 0.05], eqproj, 1.e-4,1.e-4,500,data)
				#  change signs of shifts for header
				#optm_params[0][3] *= -1
				#optm_params[0][4] *= -1
				#optm_params[0].append(imn)
				#new_params.append(optm_params[0])

				outf.write('New %6.1f  %6.1f  %6.1f  %6.1f  %6.1f   %7.4f    %d4   %7.1f'%(optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4],optm_params[1], optm_params[2], ctf_params[1]))
				outf.write("\n")
				outf.flush()
				
				#del  data[2]
				#set_arb_params(dataim[imn-image_start], [optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4]], par_str)
				#t4 = finish_time(t3)
				"""
			del data
			soto = []
			for imn in xrange(image_start, image_end):
				s2x   = dataim[imn-image_start].get_attr('s2x')
				s2y   = dataim[imn-image_start].get_attr('s2y')
				phi   = dataim[imn-image_start].get_attr('phi')
				theta = dataim[imn-image_start].get_attr('theta')
				psi   = dataim[imn-image_start].get_attr('psi')
				soto.append([phi,theta,psi,s2x,s2y,imn])
			dropSpiderDoc(os.path.join(outdir, replace("new_params%3d_%3d_%3d"%(iteration, ic, myid),' ','0')), soto," phi, theta, psi, s2x, s2y, image number")
			del soto

			#from sys import exit
			#exit()
		#  here we should write header info, just in case the program crashes...
	del vol
	del volft
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)


def ali3d_eB_MPI_select(stack, ref_vol, outdir, maskfile, ou=-1,  delta=2, maxit=10, CTF = None, snr=1.0, sym="c1", chunk = -1.0, user_func_name="ref_aliB_cone"):
	"""
		Version 03/20/08, to test how many particles to use.
	"""
	from alignment	    import eqprojEuler
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl, filt_tanl
	from fundamentals   import fshift, rot_avg_image
	from projection     import prep_vol, prgs
	from utilities      import amoeba, bcast_list_to_all, bcast_string_to_all, model_circle, get_arb_params, set_arb_params, dropSpiderDoc
	from utilities      import getImage, dropImage, bcast_EMData_to_all, send_attr_dict, recv_attr_dict
	from utilities      import readSpiderDoc, get_im
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

	parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
	#                  0                 1            2            3            4                 5                 6
	if(myid == main_node):
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)
		import user_functions
		user_func = user_functions.factory[user_func_name]
		#  ERROR if ctf applied
		ima = EMData()
		ima.read_image(stack)
		ctf_params = get_arb_params(ima, parnames)
		if(ctf_params[6] == 1):
			ERROR("ali3d_e does not work for CTF-applied data","ali3d_eB_MPI",1)
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
		if(type(maskfile) is types.StringType):  mask3D = getImage(maskfile)
		else:                                   mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)
	mask2D = model_circle(ou, nx, nx)

	dataim = []
	#from utilities import readSpiderDoc, set_arb_params
	#prm = readSpiderDoc("params_new.doc")
	#prm_dict = ["phi", "theta", "psi", "s2x", "s2y"]
	prm_dict = ["phi", "theta", "psi"]
	#from random import seed,gauss
	#seed()
	dataim = EMData.read_images(image_start, image_end)
	for im in xrange(image_start, image_end):
		dataim[im].set_attr('ID', im)
		#angn = get_arb_params(ima, prm_dict)
		#for ian in xrange(3):  angn[ian] += gauss(0.0, 1.0)
		#set_arb_params(ima, angn, prm_dict)
		#set_arb_params(ima, prm[im], prm_dict)
		'''
		if(CTF):
			ctf_params = get_arb_params(ima, parnames)
			if(im == image_start): data_had_ctf = ctf_params[6]
			if(ctf_params[6] == 0):
				st = Util.infomask(ima, mask2D, False)
				ima -= st[0]
				from filter import filt_ctf
				ima = filt_ctf(ima, ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
				ima.set_attr('ctf_applied', 1)
		'''
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
				"""
				for im in xrange(len(templ)):
					outf.write('ccc, image %12.5f  %07d'%( templ[im][0], templ[im][1]  ))
					outf.write("\n")
				outf.flush()
				"""
				
				templ.sort()
				ilow = int(chunk*len(templ))  # reject 25% worst images.
				for im in xrange(ilow):  templ[im] = [templ[im][1], 0]
				for im in xrange(ilow, len(templ)):  templ[im] = [templ[im][1], 1]
				templ.sort()
				sendbuf = []
				for im in xrange(len(templ)):	sendbuf.append(templ[im][1])
				del templ
				"""
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
				"""
			else:
				sendbuf = []
			mpi_barrier(MPI_COMM_WORLD)
			#recvbuf = mpi_scatterv(sendbuf, recvcount, disps, MPI_FLOAT, recvcount[myid], MPI_FLOAT, main_node, MPI_COMM_WORLD)
			recvbuf = mpi_scatterv(sendbuf, recvcount, disps, MPI_INT, recvcount[myid], MPI_INT, main_node, MPI_COMM_WORLD)
			del sendbuf

			for imn in xrange(image_start, image_end):
				dataim[imn-image_start].set_attr_dict({'active': int(recvbuf[imn-image_start])})
				#dataim[imn-image_start] /= float(recvbuf[imn-image_start])

			"""
			nact = 0
			for imn in xrange(image_start, image_end):
				nact += dataim[imn-image_start].get_attr('active')
			nact = float(nact)
			tn = mpi_reduce(nact, 1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
			if(myid == main_node):
				outf.write('total number of used images %12.2f  '%(float(tn)))
				outf.write("\n")
				outf.flush()
			"""
			# resolution
			outf.write("  begin reconstruction = "+str(image_start)+"   ")
			outf.write("\n")
			outf.flush()
			vol, fscc = rec3D_MPI(dataim, snr, sym, mask3D, os.path.join(outdir, replace("resolution%3d_%3d"%(iteration, ic),' ','0') ), myid, main_node)
			outf.write("  done reconstruction = "+str(image_start)+"   ")
			outf.write("\n")
			outf.flush()
			#volftb = vol.copy()

			#  restore original normalization
			#for imn in xrange(image_start, image_end):
			#	#dataim[imn-image_start].set_attr_dict({'active': int(recvbuf[imn-image_start])})
			#	dataim[imn-image_start] *= float(recvbuf[imn-image_start])
			#del recvbuf
			mpi_barrier(MPI_COMM_WORLD)
			if(myid == main_node):
				dropImage(vol, os.path.join(outdir, replace("vol%3d_%3d.hdf"%(iteration, ic),' ','0') ))
				ref_data.append( vol )
				ref_data.append( fscc )
				#  call user-supplied function to prepare reference image, i.e., filter it
				vol = user_func( ref_data )
				#  HERE CS SHOULD BE USED TO MODIFY PROJECTIONS' PARAMETERS  !!!
				del ref_data[2]
				del ref_data[2]
				dropImage(vol, os.path.join(outdir, replace("volf%3d_%3d.hdf"%(iteration, ic),' ','0') ))
			from sys import exit
			exit()

			#bcast_EMData_to_all(volftb, myid, main_node)
			bcast_EMData_to_all(vol, myid, main_node)
			#volft,kb  = prep_vol(vol)
			data = []
			# Only Euler angles
			data.append(None)
			data.append(None)
			data.append(None)
			data.append(None)
			data.append(None)
			data.append(mask2D)

			image_start_in_chunk = image_start + ic*n_in_chunk
			image_end_in_chunk   = min(image_start_in_chunk + n_in_chunk, image_end)
			outf.write("ic "+str(ic)+"   image_start "+str(image_start)+"   n_in_chunk "+str(n_in_chunk)+"   image_end "+str(image_end)+"\n")
			outf.write("image_start_in_chunk "+str(image_start_in_chunk)+"  image_end_in_chunk "+str(image_end_in_chunk)+"\n")
			outf.flush()
			if(CTF):  previous_defocus = -1.0
			for imn in xrange(image_start_in_chunk, image_end_in_chunk):
				if(CTF):
					ctf_params = get_arb_params(dataim[imn-image_start], parnames)
					if(ctf_params[1] != previous_defocus):
						previous_defocus = ctf_params[1]
						data[0],data[1] = prep_vol(filt_ctf(vol, ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5]))
				data[2] = dataim[imn-image_start]

				atparams = get_arb_params(dataim[imn-image_start], par_str)
				
				#optm_params = ali_G3(data, atparams, dtheta)
				#  Align only Euler angles
				#  change signs of shifts for projections
				data[3] = -atparams[3]
				data[4] = -atparams[4]

				#initial  = eqproj(atparams, data)  # this is if we need initial discrepancy
				outf.write("Image "+str(imn)+"\n")
				#outf.write('Old %6.1f  %6.1f  %6.1f  %6.1f  %6.1f   %7.4f '%(atparams[0],atparams[1],atparams[2],-atparams[3],-atparams[4], initial))
				outf.write('Old %8.3f  %8.3f  %8.3f  '%(atparams[0],atparams[1],atparams[2]))
				outf.write("\n")
				#outf.flush()

				weight_phi = max(delta, delta*abs((atparams[1]-90.0)/180.0*pi))
				#from utilities import start_time, finish_time
				#t3=start_time()
				optm_params =  amoeba(atparams[0:3], [weight_phi, delta, weight_phi], eqprojEuler, 1.e-4,1.e-4,500, data)
				optm_params[0].append(imn)
				#new_params.append(optm_params[0])

				#outf.write('New %6.1f  %6.1f  %6.1f  %6.1f  %6.1f   %7.4f    %d4   %7.1f'%(optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4],optm_params[1], optm_params[2], ctf_params[1]))
				outf.write('New %8.3f  %8.3f  %8.3f  %8.5f  %6.1f   '%(optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[1], optm_params[2]))
				outf.write("\n")
				outf.flush()

				set_arb_params(dataim[imn-image_start], [optm_params[0][0], optm_params[0][1], optm_params[0][2]], par_str[0:3])

				"""

				#  change signs of shifts for projections
				#atparams[3] *= -1
				#atparams[4] *= -1

				initial  = eqproj(atparams, data)  # this is if we need initial discrepancy
				outf.write("Image "+str(imn)+"\n")
				outf.write('Old %6.1f  %6.1f  %6.1f  %6.1f  %6.1f   %7.4f '%(atparams[0],atparams[1],atparams[2],-atparams[3],-atparams[4], initial))
				outf.write("\n")
				outf.flush()

				#weight_phi = max(delta, delta*abs((atparams[1]-90.0)/180.0*pi))
				#from utilities import start_time, finish_time
				#t3=start_time()
				#optm_params =  amoeba(atparams, [weight_phi, delta, weight_phi, 0.05, 0.05], eqproj, 1.e-4,1.e-4,500,data)
				#  change signs of shifts for header
				#optm_params[0][3] *= -1
				#optm_params[0][4] *= -1
				#optm_params[0].append(imn)
				#new_params.append(optm_params[0])

				outf.write('New %6.1f  %6.1f  %6.1f  %6.1f  %6.1f   %7.4f    %d4   %7.1f'%(optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4],optm_params[1], optm_params[2], ctf_params[1]))
				outf.write("\n")
				outf.flush()
				
				#del  data[2]
				#set_arb_params(dataim[imn-image_start], [optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4]], par_str)
				#t4 = finish_time(t3)
				"""
			del data
			soto = []
			for imn in xrange(image_start, image_end):
				s2x   = dataim[imn-image_start].get_attr('s2x')
				s2y   = dataim[imn-image_start].get_attr('s2y')
				phi   = dataim[imn-image_start].get_attr('phi')
				theta = dataim[imn-image_start].get_attr('theta')
				psi   = dataim[imn-image_start].get_attr('psi')
				soto.append([phi,theta,psi,s2x,s2y,imn])
			dropSpiderDoc(os.path.join(outdir, replace("new_params%3d_%3d_%3d"%(iteration, ic, myid),' ','0')), soto," phi, theta, psi, s2x, s2y, image number")
			del soto

			#from sys import exit
			#exit()
		#  here we should write header info, just in case the program crashes...
	del vol
	del volft
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)

def ali3d_eB_MPI_(stack, ref_vol, outdir, maskfile, ou=-1,  delta=2, maxit=10, CTF = None, snr=1.0, sym="c1", chunk = -1.0):
	"""
		
	"""
	from alignment	    import eqprojEuler
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl, filt_gaussl
	from fundamentals   import fshift, rot_avg_image
	from projection     import prep_vol, prgs
	from utilities      import amoeba, bcast_list_to_all, bcast_string_to_all, model_circle, get_arb_params, set_arb_params, dropSpiderDoc
	from utilities      import getImage, dropImage, bcast_EMData_to_all, send_attr_dict, recv_attr_dict
	from utilities      import readSpiderDoc, get_im
	from reconstruction import rec3D_MPI
	from math           import pi
	from string         import replace
	import os
	import sys
	from mpi 	    import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	    import mpi_barrier, mpi_gatherv, mpi_scatterv
	from mpi 	    import MPI_FLOAT, MPI_INT, MPI_SUM

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	#  chose a random node as a main one...
	main_node = 0
	if(myid == main_node):
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)
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
	vol = EMData()
	vol.read_image(ref_vol)
	nx  = vol.get_xsize()
	if(ou <= 0):  ou = nx//2-2
	#if(myid == main_node):  vol.write_image(os.path.join(outdir,"ref_volf00.hdf"))
	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask3D = getImage(maskfile)
		else:                                   mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)
	mask2D = model_circle(ou, nx, nx)

	dataim = []

	for im in xrange(image_start, image_end):
		dataim.append(im - image_start)

	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                  disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(nima, number_of_proc, im)
		recvcount.append( ie - ib )
	recvbuf = mpi_gatherv(dataim, len(dataim), MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
	mpi_barrier(MPI_COMM_WORLD)
	if(myid == main_node):
		print  "recvcount", recvcount
		print  "disps", disps
		print  "len(recvbuf)",len(recvbuf)
		print  recvbuf
		templ = []
		for im in xrange(len(recvbuf)):	   templ.append([float(recvbuf[im]), im])
		del recvbuf
		print  " templ ",templ
		templ.sort()
		ilow = int(0.3*len(templ))  # skip 30%
		for im in xrange(ilow):  templ[im] = [templ[im][1], 0]
		for im in xrange(ilow, len(templ)):  templ[im] = [templ[im][1], 1]
		templ.sort()
		print  " templ 2",templ
		sendbuf = []
		for im in xrange(len(templ)):   sendbuf.append(templ[im][1])
		del templ
	else:
		sendbuf = []
	mpi_barrier(MPI_COMM_WORLD)
	recvbuf = mpi_scatterv(sendbuf, recvcount, disps, MPI_INT, recvcount[myid], MPI_INT, main_node, MPI_COMM_WORLD)
	del sendbuf
	print  " MY ID ", myid, image_start, image_end, len(dataim), len(recvbuf)
	print  " MY ID ", myid, recvbuf
	from sys import exit
	exit()
	
	parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
	#                  0                 1            2            3            4                 5                 6
	from utilities import readSpiderDoc, set_arb_params
	prm = readSpiderDoc("params_new.doc")
	prm_dict = ["phi", "theta", "psi", "s2x", "s2y"]
	dataim = EMData.read_images(image_start, image_end)
	for im in xrange(image_start, image_end):
		dataim[im].set_attr('ID', im)
		set_arb_params(ima, prm[im], prm_dict)
		if(CTF):
			ctf_params = get_arb_params(data[im], parnames)
			if(im == image_start): data_had_ctf = ctf_params[6]
			if(ctf_params[6] == 0):
				st = Util.infomask(data[im], mask2D, False)
				data[im] -= st[0]
				from filter import filt_ctf
				data[im] = filt_ctf(data[im], ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
				data[im].set_attr('ctf_applied', 1)
	outf.write("  data read = "+str(image_start)+"   ")
	outf.write("\n")
	outf.flush()

	from  string        import replace
	from utilities      import bcast_number_to_all
	from morphology     import adaptive_mask, threshold, threshold_to_minval
	from statistics     import histogram
	from reconstruction import recons3d_nn_SSNR_MPI

	par_str=["phi", "theta", "psi", "s2x", "s2y"]
	'''
	#  Berlin
	adam = adaptive_mask(vol, 3000, 2.22)
	vol = threshold( adam*vol)
	h = histogram( vol )
	vol *= get_im("rotav_tteftu_with_tRNA.spi") / threshold_to_minval( rot_avg_image(vol), h[len(h)//2])
	vol = filt_btwl(vol, 0.3, 0.4)
	'''

	ldef = -1
	for iteration in xrange(maxit):
		outf.write("  iteration = "+str(iteration)+"   ")
		outf.write("\n")
		outf.flush()
		for  ic  in xrange(n_of_chunks):
			# compute updated 3D after each chunk
 	    		# resolution
			"""
			outf.write("  begin ssnr = "+str(image_start)+"   ")
			outf.write("\n")
			outf.flush()
			if myid == 0:
				[ssnr1, vol_ssnr] = recons3d_nn_SSNR_MPI(myid, dataim, model_circle(ou, nx, nx), ctf = CTF)
				del ssnr1
				#  change volume to fsc
				vol_ssnr /= vol_ssnr + 1.0
				#  filter it somewhat
				vol_ssnr = threshold(filt_gaussl(vol_ssnr, 0.1))
			else:    recons3d_nn_SSNR_MPI(myid, dataim, model_circle(ou, nx, nx), ctf = CTF)
			outf.write("  begin reconstruction = "+str(image_start)+"   ")
			outf.write("\n")
			outf.flush()
			"""
			volo, fscc = rec3D_MPI(dataim, snr, sym, mask3D, os.path.join(outdir, replace("resolution%3d_%3d"%(iteration, ic),' ','0') ), myid, main_node)
			outf.write("  done reconstruction = "+str(image_start)+"   ")
			outf.write("\n")
			outf.flush()
			mpi_barrier(MPI_COMM_WORLD)
			if(myid == main_node):
				dropImage(volo, os.path.join(outdir, replace("vol%3d_%3d.hdf"%(iteration, ic),' ','0') ))
				from filter import fit_tanh, filt_tanl
				fl, aa = fit_tanh(fscc)
				vol = filt_tanl(vol, fl, aa)
				from math import exp, sqrt, tanh, pi
				#parameters of tanh filter are in absolute units.
				nct = int(1.8*nx)
				fl = 0.45
				aa = 0.1
				envt = []
				for i in xrange(nct):
					# e(x)*h(x)/(bckg(x)+e(x)**2*ctf(x)**2/20)
					xs = float(i)/2.22/nx
					et = exp(-100.0*xs**2)
					bckgt = exp(-0.8-100.0*xs**2)+0.01
					ht = 1.0-0.6*exp(-xs**2/2.0/0.012**2)
					H = 0.5*( tanh(pi*(xs*ctf_params[0]+fl)/2./aa/fl) - tanh(pi*(xs*ctf_params[0]-fl)/2./aa/fl) )
					fmt = H*ht/bckgt
					envt.append(fmt)
				from filter import filt_table
				vhlf = threshold(filt_table(vol, envt))
				Util.mul_img(vhlf, mask3D)
				dropImage(vol, os.path.join(outdir, replace("vhlf%3d_%3d.hdf"%(iteration, ic),' ','0') ))
				"""
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
					if(qtq > 0.0):
						lqtq = filt[irf]/qtq
						rfsc.set_value_at(irf, lqtq)
					else:            rfsc.set_value_at(irf, lqtq)
				vol_ssnr = vol_ssnr.mult_radial(rfsc)  # adjust radially fsc filter to fsc curve
				dropImage(vol_ssnr,  os.path.join(outdir, replace("filter%3d_%3d.hdf"%(iteration, ic),' ','0') ))
				#  Filter by volume
				vol = volo.filter_by_image(vol_ssnr)
				#vol  = filt_table(volo, filt)
				#if(center == 1):
				#	cs   = vol.phase_cog()
				#	vol  = fshift(vol, -cs[0], -cs[1] -cs[2])
				#	volo = fshift(volo, -cs[0], -cs[1] -cs[2])
				dropImage(vol, os.path.join(outdir, replace("volf%3d_%3d.hdf"%(iteration, ic),' ','0') ))
				adam = adaptive_mask(vol, 4000, 2.22)
				vol = threshold( adam*vol )
				h = histogram( vol )
				vol = threshold( adam*volo ) * get_im("rotav_tteftu_with_tRNA.spi") / threshold_to_minval( rot_avg_image(vol), h[len(h)//2])
				del volo
				#  Filter by volume
				vol = vol.filter_by_image(vol_ssnr)
				#vol  = filt_table(vol, filt)
				#vol = filt_btwl(vol, fl, fh)
				dropImage(vol, os.path.join(outdir, replace("vhlf%3d_%3d.hdf"%(iteration, ic),' ','0') ))
				del adam, vol_ssnr, rfsc, filt
				del h
				"""
			bcast_EMData_to_all(vol,  myid, main_node)
			bcast_EMData_to_all(vhlf, myid, main_node)
			from filter import  filt_tophatb
			volft,kb  = prep_vol( filt_tophatb(vhlf, 0.20, 0.40, False))
			from statistics import ccc
			qcc = []
			for imn in xrange(image_start, image_end):
				atparams = get_arb_params(dataim[imn-image_start], par_str)
				projt = prgs(volft, kb, [atparams[0], atparams[1], atparams[2], -atparams[3], atparams[4]])
				qcc.append(ccc(projt, dataim[imn-image_start], mask2D))
			del projt




			volft,kb  = prep_vol(vol)
			data = []
			data.append(volft)
			data.append(kb)
			data.append(mask2D)



			image_start_in_chunk = image_start + ic*n_in_chunk
			image_end_in_chunk   = min(image_start_in_chunk + n_in_chunk, image_end)
			outf.write("ic "+str(ic)+"   image_start "+str(image_start)+"   n_in_chunk "+str(n_in_chunk)+"   image_end "+str(image_end)+"\n")
			outf.write("image_start_in_chunk "+str(image_start_in_chunk)+"  image_end_in_chunk "+str(image_end_in_chunk)+"\n")
			outf.flush()
			#new_params = []

			for imn in xrange(image_start_in_chunk, image_end_in_chunk):
				ctf_params = get_arb_params(dataim[imn-image_start], parnames)
				from morphology import ctf_2
				ctf2 = ctf_2(nx, ctf_params[0], ctf_params[1])
				nct = len(ctf2)
				from math import exp, sqrt, tanh, pi
				#parameters of tanh filter are in absolute units.
				fl = 0.45
				aa = 0.1
				envt = []
				for i in xrange(nct):
					# e(x)*h(x)/(bckg(x)+e(x)**2*ctf(x)**2/20)
					xs = float(i)/2.22/nx
					et = exp(-100.0*xs**2)
					bckgt = exp(-0.8-100.0*xs**2)+0.01
					ht = 1.0-0.6*exp(-xs**2/2.0/0.012**2)
					H = 0.5*( tanh(pi*(xs*ctf_params[0]+fl)/2./aa/fl) - tanh(pi*(xs*ctf_params[0]-fl)/2./aa/fl) )
					fmt = H*ht/(bckgt + ctf2[i]*et**2/5.0)
					envt.append(fmt)
				from filter import filt_table
				ima = filt_table(dataim[imn-image_start], envt)
				data.insert(2, ima)

				atparams = get_arb_params(dataim[imn-image_start], par_str)
				#  Align only Euler angles
				#  change signs of shifts for projections
				data.insert(3, -atparams[3])
				data.insert(4, -atparams[4])
				"""
				initial  = eqproj(atparams, data)  # this is if we need initial discrepancy
				outf.write("Image "+str(imn)+"\n")
				outf.write('Old %6.1f  %6.1f  %6.1f  %6.1f  %6.1f   %7.4f '%(atparams[0],atparams[1],atparams[2],-atparams[3],-atparams[4], initial))
				outf.write("\n")
				outf.flush()
				"""
				weight_phi = max(delta, delta*abs((atparams[1]-90.0)/180.0*pi))
				#from utilities import start_time, finish_time
				#t3=start_time()
				optm_params =  amoeba(atparams[0:3], [weight_phi, delta, weight_phi], eqprojEuler, 1.e-4,1.e-4,500, data)
				optm_params[0].append(imn)
				#new_params.append(optm_params[0])
				"""
				outf.write('New %6.1f  %6.1f  %6.1f  %6.1f  %6.1f   %7.4f    %d4   %7.1f'%(optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4],optm_params[1], optm_params[2], ctf_params[1]))
				outf.write("\n")
				outf.flush()
				"""
				del  data[2]
				del  data[2]
				del  data[2]
				set_arb_params(dataim[imn-image_start], [optm_params[0][0], optm_params[0][1], optm_params[0][2]], par_str[0:3])
				#t4 = finish_time(t3)
			soto = []
			for imn in xrange(image_start, image_end):
				s2x   = dataim[imn-image_start].get_attr('s2x')
				s2y   = dataim[imn-image_start].get_attr('s2y')
				phi   = dataim[imn-image_start].get_attr('phi')
				theta = dataim[imn-image_start].get_attr('theta')
				psi   = dataim[imn-image_start].get_attr('psi')
				soto.append([phi,theta,psi,s2x,s2y,imn])
			dropSpiderDoc(os.path.join(outdir, replace("new_params%3d_%3d_%3d"%(iteration, ic, myid),' ','0')), soto," phi, theta, psi, s2x, s2y, image number")
			del soto

			'''
			#filt = filt_from_fsc(fscc, 0.05)
			#vol  = filt_table(vol, filt)
			lowfq,highfq = filt_params(fscc)
			vol  = filt_btwl(vol, lowfq, highfq)
			#  center!
			cs = vol.phase_cog()
			vol = fshift(vol, -cs[0], -cs[1] -cs[2])
			'''

			#from sys import exit
			#exit()
		#  here we should write header info, just in case the program crashes...
	del vol
	del volft
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
     
	#if(CTF and data_had_ctf == 0):
	#	for im in xrange(image_start, image_end): data[im-image_start].set_attr('ctf_applied', 0)
	#if(myid == main_node): recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	#else: send_attr_dict(main_node, data, par_str, image_start, image_end)


def ali3d_eB_MPI__(stack, ref_vol, outdir, maskfile, ou=-1,  delta=2, maxit=10, CTF = None, snr=1.0, sym="c1", chunk = -1.0):
	"""
		
	"""
	from alignment	    import eqproj
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc, filt_btwl, filt_gaussl
	from fundamentals   import fshift, rot_avg_image
	from projection     import prep_vol
	from utilities      import amoeba, bcast_list_to_all, bcast_string_to_all, model_circle, get_arb_params, set_arb_params, dropSpiderDoc
	from utilities      import getImage, dropImage, bcast_EMData_to_all, send_attr_dict, recv_attr_dict
	from utilities      import readSpiderDoc, get_im
	from reconstruction import rec3D_MPI
	from math           import pi
	from string         import replace
	import os
	import sys
	from mpi 	    import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	    import mpi_barrier, mpi_gatherv, mpi_scatterv
	from mpi 	    import MPI_FLOAT, MPI_INT, MPI_SUM

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	#  chose a random node as a main one...
	main_node = 0
	if(myid == main_node):
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)
	mpi_barrier(MPI_COMM_WORLD)

	nima = EMUtil.get_image_count(stack)

	#nimage_per_node = nima/number_of_proc
	#image_start = myid * nimage_per_node
	#if(myid == number_of_proc-1):  image_end = nima
	#else:                          image_end = image_start + nimage_per_node
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
	vol = EMData()
	vol.read_image(ref_vol)
	nx  = vol.get_xsize()
	if(ou <= 0):  ou = nx//2-2
	#if(myid == main_node):  vol.write_image(os.path.join(outdir,"ref_volf00.hdf"))
	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask3D=getImage(maskfile)
		else:                                   mask3D = maskfile
	else:
		mask3D = model_circle(ou, nx, nx, nx)
	mask2D = model_circle(ou, nx, nx)

	dataim = []
	parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
	#                  0                 1            2            3            4                 5                 6
	#from utilities import readSpiderDoc, set_arb_params
	#prm = readSpiderDoc("params_new.doc")
	#prm_dict = ["phi", "theta", "psi", "s2x", "s2y"]
	for im in xrange(image_start, image_end):
		ima = EMData()
		ima.read_image(stack, im)
		#set_arb_params(ima, prm[im], prm_dict)
		if(CTF):
			ctf_params = get_arb_params(ima, parnames)
			if(im == image_start): data_had_ctf = ctf_params[6]
			if(ctf_params[6] == 0):
				st = Util.infomask(ima, mask2D, False)
				ima -= st[0]
				from filter import filt_ctf
				ima = filt_ctf(ima, ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
				ima.set_attr('ctf_applied', 1)
		dataim.append(ima)
	outf.write("  data read = "+str(image_start)+"   ")
	outf.write("\n")
	outf.flush()

	from  string        import replace
	from utilities      import bcast_number_to_all
	from morphology     import adaptive_mask, threshold, threshold_to_minval
	from statistics     import histogram
	from reconstruction import recons3d_nn_SSNR_MPI

	par_str=["phi", "theta", "psi", "s2x", "s2y"]
	'''
	#  Berlin
	adam = adaptive_mask(vol, 3000, 2.22)
	vol = threshold( adam*vol)
	h = histogram( vol )
	vol *= get_im("rotav_tteftu_with_tRNA.spi") / threshold_to_minval( rot_avg_image(vol), h[len(h)//2])
	vol = filt_btwl(vol, 0.3, 0.4)
	'''

	ldef = -1
	for iteration in xrange(maxit):
		outf.write("  iteration = "+str(iteration)+"   ")
		outf.write("\n")
		outf.flush()
		for  ic  in xrange(n_of_chunks):
			# compute updated 3D after each chunk
 	    		# resolution
			outf.write("  begin ssnr = "+str(image_start)+"   ")
			outf.write("\n")
			outf.flush()
			if myid == 0:
				[ssnr1, vol_ssnr] = recons3d_nn_SSNR_MPI(myid, dataim, model_circle(ou, nx, nx), ctf = CTF)
				del ssnr1
				#  change volume to fsc
				vol_ssnr /= vol_ssnr + 1.0
				#  filter it somewhat
				vol_ssnr = threshold(filt_gaussl(vol_ssnr, 0.1))
			else:    recons3d_nn_SSNR_MPI(myid, dataim, model_circle(ou, nx, nx), ctf = CTF)
			outf.write("  begin reconstruction = "+str(image_start)+"   ")
			outf.write("\n")
			outf.flush()
			volo, fscc = rec3D_MPI(dataim, snr, sym, mask3D, os.path.join(outdir, replace("resolution%3d_%3d"%(iteration, ic),' ','0') ), myid, main_node)
			outf.write("  done reconstruction = "+str(image_start)+"   ")
			outf.write("\n")
			outf.flush()
			mpi_barrier(MPI_COMM_WORLD)
			if(myid == main_node):
				dropImage(volo, os.path.join(outdir, replace("vol%3d_%3d.hdf"%(iteration, ic),' ','0') ))
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
					if(qtq > 0.0):
						lqtq = filt[irf]/qtq
						rfsc.set_value_at(irf, lqtq)
					else:            rfsc.set_value_at(irf, lqtq)
				vol_ssnr = vol_ssnr.mult_radial(rfsc)  # adjust radially fsc filter to fsc curve
				dropImage(vol_ssnr,  os.path.join(outdir, replace("filter%3d_%3d.hdf"%(iteration, ic),' ','0') ))
				#  Filter by volume
				vol = volo.filter_by_image(vol_ssnr)
				#vol  = filt_table(volo, filt)
				#if(center == 1):
				#	cs   = vol.phase_cog()
				#	vol  = fshift(vol, -cs[0], -cs[1] -cs[2])
				#	volo = fshift(volo, -cs[0], -cs[1] -cs[2])
				dropImage(vol, os.path.join(outdir, replace("volf%3d_%3d.hdf"%(iteration, ic),' ','0') ))
				adam = adaptive_mask(vol, 4000, 2.22)
				vol = threshold( adam*vol )
				h = histogram( vol )
				vol = threshold( adam*volo ) * get_im("rotav_tteftu_with_tRNA.spi") / threshold_to_minval( rot_avg_image(vol), h[len(h)//2])
				del volo
				#  Filter by volume
				vol = vol.filter_by_image(vol_ssnr)
				#vol  = filt_table(vol, filt)
				#vol = filt_btwl(vol, fl, fh)
				dropImage(vol, os.path.join(outdir, replace("vhlf%3d_%3d.hdf"%(iteration, ic),' ','0') ))
				del adam, vol_ssnr, rfsc, filt
				del h
			bcast_EMData_to_all(vol, myid, main_node)

			volft,kb  = prep_vol(vol)
			data = []
			data.append(volft)
			data.append(kb)
			data.append(mask2D)

			image_start_in_chunk = image_start + ic*n_in_chunk
			image_end_in_chunk   = min(image_start_in_chunk + n_in_chunk, image_end)
			outf.write("ic "+str(ic)+"   image_start "+str(image_start)+"   n_in_chunk "+str(n_in_chunk)+"   image_end "+str(image_end)+"\n")
			outf.write("image_start_in_chunk "+str(image_start_in_chunk)+"  image_end_in_chunk "+str(image_end_in_chunk)+"\n")
			outf.flush()

			#new_params = []

			for imn in xrange(image_start_in_chunk, image_end_in_chunk):
				'''
				defocus = dataim[imn-image_start].get_attr('defocus')

				if(defocus != current_defocus):
					#maft = dataim[imn-image_start].get_attr('matched_filter')
					ldef += 1
					if(ldef == 221): ldef = 0
					c1 = mafi_params[ldef][0]
					c2 = mafi_params[ldef][1]
					c3 = mafi_params[ldef][2]
					c4 = mafi_params[ldef][3]
					maft = []
					for j in xrange(150):
						x = j/2.22/2.0/75.0
						B  = exp(c1*x**2+c2)
						E  = (1.-exp(-(x+c3)**2/2./.01**2))/(1.+c4*x**2.2)
						ffil =  E/B*.5*(tanh(pi*(x+fl)/2./fl/a)-tanh(pi*(x-fl)/2./fl/a))
						maft.append(ffil)

					from filter import filt_table
					vol_maft = filt_table(vol, maft)
					volft,kb  = prep_vol(vol_maft)
					#volft,kb  = prep_vol(vol)
					data = []
					data.append(volft)
					data.append(kb)
					data.append(mask2D)
					current_defocus = defocus
				'''
				ctf_params = get_arb_params(dataim[imn-image_start], parnames)
				from morphology import ctf_2
				ctf2 = ctf_2(nx, ctf_params[0], ctf_params[1])
				nct = len(ctf2)
				from math import exp, sqrt, tanh, pi
				#parameters of tanh filter are in absolute units.
				fl = 0.45
				aa = 0.1
				envt = []
				for i in xrange(nct):
					# e(x)*h(x)/(bckg(x)+e(x)**2*ctf(x)**2/20)
					xs = float(i)/2.22/nx
					et = exp(-100.0*xs**2)
					bckgt = exp(-0.8-100.0*xs**2)+0.01
					ht = 1.0-0.6*exp(-xs**2/2.0/0.012**2)
					H = 0.5*( tanh(pi*(xs*ctf_params[0]+fl)/2./aa/fl) - tanh(pi*(xs*ctf_params[0]-fl)/2./aa/fl) )
					fmt = H*ht/(bckgt + ctf2[i]*et**2/5.0)
					envt.append(fmt)
				from filter import filt_table
				ima = filt_table(dataim[imn-image_start], envt)
				data.insert(2, ima)

				atparams = get_arb_params(dataim[imn-image_start], par_str)
				#  change signs of shifts for projections
				atparams[3] *= -1
				atparams[4] *= -1
				"""
				initial  = eqproj(atparams, data)  # this is if we need initial discrepancy
				outf.write("Image "+str(imn)+"\n")
				outf.write('Old %6.1f  %6.1f  %6.1f  %6.1f  %6.1f   %7.4f '%(atparams[0],atparams[1],atparams[2],-atparams[3],-atparams[4], initial))
				outf.write("\n")
				outf.flush()
				"""
				weight_phi = max(delta, delta*abs((atparams[1]-90.0)/180.0*pi))
				#from utilities import start_time, finish_time
				#t3=start_time()
				optm_params =  amoeba(atparams, [weight_phi, delta, weight_phi, 0.05, 0.05], eqproj, 1.e-4,1.e-4,500,data)
				#  change signs of shifts for header
				optm_params[0][3] *= -1
				optm_params[0][4] *= -1
				optm_params[0].append(imn)
				#new_params.append(optm_params[0])
				"""
				outf.write('New %6.1f  %6.1f  %6.1f  %6.1f  %6.1f   %7.4f    %d4   %7.1f'%(optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4],optm_params[1], optm_params[2], ctf_params[1]))
				outf.write("\n")
				outf.flush()
				"""
				del  data[2]
				set_arb_params(dataim[imn-image_start], [optm_params[0][0], optm_params[0][1], optm_params[0][2], optm_params[0][3], optm_params[0][4]], par_str)
				#t4 = finish_time(t3)
			soto = []
			for imn in xrange(image_start, image_end):
				s2x   = dataim[imn-image_start].get_attr('s2x')
				s2y   = dataim[imn-image_start].get_attr('s2y')
				phi   = dataim[imn-image_start].get_attr('phi')
				theta = dataim[imn-image_start].get_attr('theta')
				psi   = dataim[imn-image_start].get_attr('psi')
				soto.append([phi,theta,psi,s2x,s2y,imn])
			dropSpiderDoc(os.path.join(outdir, replace("new_params%3d_%3d_%3d"%(iteration, ic, myid),' ','0')), soto," phi, theta, psi, s2x, s2y, image number")
			del soto

			'''
			#filt = filt_from_fsc(fscc, 0.05)
			#vol  = filt_table(vol, filt)
			lowfq,highfq = filt_params(fscc)
			vol  = filt_btwl(vol, lowfq, highfq)
			#  center!
			cs = vol.phase_cog()
			vol = fshift(vol, -cs[0], -cs[1] -cs[2])
			'''

			#from sys import exit
			#exit()
		#  here we should write header info, just in case the program crashes...
	del vol
	del volft
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
     
	#if(CTF and data_had_ctf == 0):
	#	for im in xrange(image_start, image_end): data[im-image_start].set_attr('ctf_applied', 0)
	#if(myid == main_node): recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	#else: send_attr_dict(main_node, data, par_str, image_start, image_end)

"""
def ssnr3d(stack, output_volume = None, ssnr_text_file = None, mask = None, ou = -1, rw = 1.0,  npad = 1, CTF = False, sign = 1, sym ="c1", random_angles = 0):
"""
"""
	Perform 3D reconstruction using selected particle images, 
	and calculate spectrum signal noise ratio (SSNR).
	1. The selection file is supposed to be in SPIDER text format.
	2. The 3D alignment parameters have been written in headers of the particle images.
""" 
""" 

	from global_def import MPI
	if MPI:
		ssnr3d_MPI(stack, output_volume, ssnr_text_file, mask, ou, rw, npad, CTF, sign, sym, random_angles)
		return

	from utilities import readSpiderDoc, dropImage, get_arb_params, set_arb_params, model_circle, get_im
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
	dropImage(vol_ssnr1, output_volume+"1.spi", "s")
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
	from utilities import dropSpiderDoc19289
	dropSpiderDoc(ssnr_text_file+".doc", ssnr1)
	dropImage(vol_ssnr2, output_volume+"2.spi", "s")
	'''

def ssnr3d_MPI(stack, output_volume = None, ssnr_text_file = None, mask = None, ou = -1, rw = 1.0, npad = 1, CTF = False, sign = 1, sym ="c1", random_angles = 0):
	from reconstruction import recons3d_nn_SSNR_MPI, recons3d_4nn_MPI, recons3d_4nn_ctf_MPI
	from string import replace
	from time import time
	from utilities import info, get_arb_params, set_arb_params, bcast_EMData_to_all, model_blank, model_circle, get_im, dropImage
	from filter import filt_ctf
	from projection import prep_vol, prgs

	if output_volume  is None: output_volume  = "SSNR.spi"
	if ssnr_text_file is None: ssnr_text_file = "ssnr"
	from utilities import readSpiderDoc

	#   TWO NEXT LINEs FOR EXTRA PROJECT
	#params_ref = readSpiderDoc("params_new.txt")
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
		dropImage(vol_ssnr1, output_volume+"1.spi", "s")
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
		from utilities import dropSpiderDoc
		dropSpiderDoc(ssnr_text_file+".doc", ssnr1)
		dropImage(vol_ssnr2, output_volume+"2.spi", "s")
	'''
"""
	
def ali3d_f(stack, ref_vol, outdir, maskfile, ali_maskfile, radius=-1, snr=1.0, dtheta=2, max_it=10, symmetry="c1", CTF = None, chunk = -1.0, MPI=False):
	if MPI:
		ali3d_f_MPI(stack, ref_vol, outdir, maskfile, ali_maskfile, radius, snr, dtheta, max_it,symmetry, CTF, chunk)
		return
	print 'non-MPI version not implemented!'
			
def ali3d_f_MPI(stack, ref_vol, outdir, maskfile, ali_maskfile, radius=-1, snr=1.0, dtheta=2, max_it=10, symmetry="c1", CTF = None, chunk = -1.0):
	"""
		
	"""
	from alignment	    import eqproj
	from filter         import filt_ctf, filt_params, filt_table, filt_from_fsc2, filt_btwl
	from fundamentals   import fshift
	from projection     import prep_vol
	from utilities      import amoeba, bcast_list_to_all, bcast_string_to_all, model_circle, get_arb_params, set_arb_params, dropSpiderDoc
	from utilities      import getImage, dropImage, bcast_EMData_to_all, send_attr_dict, recv_attr_dict,get_im
	from utilities      import readSpiderDoc
	from reconstruction import rec3D_MPI,rec3D_MPI_index
	from morphology     import refine_with_mask
	from math           import pi
	from random         import randint
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
	mpi_barrier(MPI_COMM_WORLD)

	nima = EMUtil.get_image_count(stack)

	#nimage_per_node = nima/number_of_proc
	#image_start = myid * nimage_per_node
	#if(myid == number_of_proc-1):             image_end = nima
	#else:                                     image_end = image_start + nimage_per_node
	image_start, image_end = MPI_start_end(nima, number_of_proc, myid)	

	# figure the size of the chunk (3D is updated after each chunk).  Chunk should be given as 0.0< chunk <= 1.0.  1.0 means all projections
	if(chunk <= 0.0):  chunk = 1.0
	
	outf = file(replace("progress%4d"%myid,' ','0'), "w")
	outf.write("  chunk = "+str(chunk)+"   ")
	outf.write("\n")
	outf.flush()
	refvols = []
	nrefvol = EMUtil.get_image_count(ref_vol)
	for iref in xrange(nrefvol):
	    vol = get_im(ref_vol, iref)
	    refvols.append(vol)

	nx  = refvols[0].get_xsize()
	if(radius <= 0):  radius = nx/2-1
	if(myid == main_node):
	    for i in xrange(nrefvol):
	        filename = "ref_volf%4d.hdf" % i
		filename = replace(filename, ' ', '0')
	        refvols[i].write_image(os.path.join(outdir,filename))

	if maskfile:
		import  types
		if(type(maskfile) is types.StringType):  mask3D=getImage(maskfile)
		else: mask3D = maskfile
	else:
		mask3D = model_circle(radius, nx, nx, nx)

        if ali_maskfile:
            ali_mask = getImage( ali_maskfile )
        else:
            print "Error: no alignment mask"
            exit(-1)


	mask2D = model_circle(radius, nx, nx)
	parnames = ["Pixel_size", "defocus", "voltage", "Cs", "amp_contrast", "B_factor",  "ctf_applied"]
	#                  0         1           2       3           4             5                 6
	dataim = EMData.read_images(image_start, image_end)
	for im in xrange(image_start, image_end):
		dataim[im].set_attr('ID', im)
		if(CTF):
			ctf_params = get_arb_params(data[im], parnames)
			if(im == image_start): data_had_ctf = ctf_params[6]
			if(ctf_params[6] == 0):
				st = Util.infomask(data[im], mask2D, False)
				data[im] -= st[0]
				from filter import filt_ctf
				data[im] = filt_ctf(data[im], ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4])
				data[im].set_attr('ctf_applied', 1)


	jtep = 0
	par_str=["phi", "theta", "psi", "s2x", "s2y"]
	from math import exp, tanh, pi
	a = 0.08
	ldef = -1
	for iteration in xrange(max_it):
                for imn in xrange(image_start, image_end):
                        igroup = randint(0, 1)
                        dataim[imn-image_start].set_attr( 'group', igroup )
                        if nrefvol==2:
                            dataim[imn-image_start].set_attr( 'refid', 0 )
                '''
                if iteration==0 :
                    n_in_chunk  = (image_end-image_start)
                    n_of_chunks = 1
                else :
                    n_in_chunk  = max( int( chunk*(image_end-image_start+1) ), 1 )
	            n_of_chunks = (image_end-image_start+1)//n_in_chunk + min((image_end-image_start+1)%n_in_chunk,1)
                '''

	        n_in_chunk  = max( int( chunk*(image_end-image_start+1) ), 1 )
	        n_of_chunks = (image_end-image_start+1)//n_in_chunk + min((image_end-image_start+1)%n_in_chunk,1)
        

                outf.write( "nchunk, chunksize: %d %d\n" % (n_of_chunks,n_in_chunk) ) 
                outf.flush()

		for  ic  in xrange(n_of_chunks):
			image_start_in_chunk = image_start + ic*n_in_chunk
			image_end_in_chunk   = min(image_start_in_chunk + n_in_chunk, image_end)
			outf.write("image_start_in_chunk "+str(image_start_in_chunk)+"\n")
			outf.write("\n")
			outf.write("image_end_in_chunk "+str(image_end_in_chunk)+"\n")
			outf.write("\n")
			outf.flush()
			jtep += 1

                        volfts = []
			kbs = []
                        for i in xrange(nrefvol):
			    vol = refvols[i]

			    Util.mul_img(vol, mask3D)
			    volft,kb  = prep_vol(vol)

                            volfts.append(volft)
			    kbs.append(kb)

			new_params = []
			current_defocus = [-1.0]*nrefvol
			for imn in xrange(image_start_in_chunk, image_end_in_chunk):
                                igroup = dataim[imn-image_start].get_attr('group')
               			data = [None]*4
			        data[0] = None
			        data[1] = None
				data[2] = dataim[imn-image_start]
			        data[3] = mask2D
				defocus = dataim[imn-image_start].get_attr('defocus')
        			atparams = get_arb_params(dataim[imn-image_start], par_str)
				atparams[3] *= -1
				atparams[4] *= -1


                                best_ref = None
                                best_score = None
				outf.write("Image "+str(imn)+"\n")
                                for iref in xrange(nrefvol/2):
                                    data[0] = volfts[2*iref+igroup]
                                    data[1] = kbs[2*iref+igroup]
           			    initial  = eqproj(atparams, data)  # this is if we need initial discrepancy
				    outf.write('Old %6.2f  %6.2f  %6.2f  %6.2f  %6.2f   %7.4f %d '%(atparams[0],atparams[1],atparams[2],-atparams[3],-atparams[4], initial,iref))
				    outf.write("\n")
				    outf.flush()
                                    if best_score is None or best_score < initial :
                                        best_ref = iref
                                        best_score = initial
                                
                                if iteration%4 == 5 :
                                    atparams.append( best_ref )
                                    new_params.append( atparams )
                                    dataim[imn-image_start].set_attr( 'refid', best_ref )
                                    outf.write( "Image %d assigned to vol %d on iter %d\n\n" % (imn, best_ref, iteration) )
                                    outf.flush()
                                    continue
                                 

                                weight_phi = max(dtheta, dtheta*abs((atparams[1]-90.0)/180.0*pi))
                               
                                best_ref = None
                                best_score = None
                                best_params = None
                                for iref in xrange(nrefvol/2):
                                    data[0] = volfts[2*iref+igroup]
                                    data[1] = kbs[2*iref+igroup]
	 			    optm_params,curt_score,niter =  amoeba(atparams, [weight_phi,dtheta,weight_phi,1.0,1.0], eqproj, 1.e-4,1.e-4,500,data)
				    #  change signs of shifts for header
				    optm_params[3] *= -1
				    optm_params[4] *= -1
				    optm_params.append(imn)
				    outf.write('New %6.2f  %6.2f  %6.2f  %6.2f  %6.2f   %7.4f   %7.1f\n'%(optm_params[0],optm_params[1],optm_params[2],optm_params[3],optm_params[4],curt_score,defocus))
				    outf.flush()
                                    if(best_score is None or best_score < curt_score):
                                        best_ref = iref
                                        best_score = curt_score
                                        best_params = optm_params

                                best_params.append(best_ref)

				new_params.append(best_params)

                                dataim[imn-image_start].set_attr( 'refid', best_ref )
                                outf.write( "Image %d assigned to vol %d on iter %d\n\n" % (imn, best_ref, iteration) )
                                outf.flush()

				set_arb_params(dataim[imn-image_start], [best_params[0], best_params[1], best_params[2], best_params[3], best_params[4]], par_str)

			dropSpiderDoc(os.path.join(outdir, replace("new_params%6d_%3d"%(jtep, myid),' ','0')), new_params," phi, theta, psi, s2x, s2y, refid, image number")
			# compute updated 3D after each chunk
 	    		# resolution
                        if(nrefvol==1):
 	                    refvols[0], fscc, oddvol, evevol = rec3D_MPI(dataim, snr, symmetry, mask3D, os.path.join(outdir, replace("resolution%6d"%jtep,' ','0')), myid, main_node)
			    if(myid == main_node):
				dropImage(refvols[0], os.path.join(outdir, replace("vol%6d.spi"%jtep,' ','0')), "s")
				#filt = filt_from_fsc(fscc, 0.05)
				#vol  = filt_table(vol, filt)
				lowfq,highfq = filt_params(fscc)
				refvols[0] = filt_btwl(refvols[0], lowfq, highfq)
				#  center!
				cs = vol.phase_cog()
				refvols[0] = fshift(refvols[0], -cs[0], -cs[1] -cs[2])
				dropImage(refvols[0],os.path.join(outdir, replace("volf%6d.spi"%jtep,' ','0')), "s")
			    bcast_EMData_to_all(refvols[0], myid, main_node)
                        elif nrefvol==2:
                            iref = 0
                            fscfile = os.path.join(outdir, replace("resolution%4d%2d"%(jtep,iref), ' ', '0'))
                            totvol, fscc, refvols[0], refvols[1] = rec3D_MPI_index( dataim, iref, snr, symmetry, mask3D, myid, main_node, fscfile )
  	                    if myid==main_node:
                                dropImage(totvol, os.path.join(outdir, replace("vol%4d%2d.spi"%(jtep,iref),' ','0')), "s")
                                #filt = filt_from_fsc2(fscc, 0.05)
                                #totvol = filt_table(totvol, filt)
                                if(fscc[1][0]<0.5) : fscc[1][0] = 1.0
                                lowfq,hghfq = filt_params(fscc)
                                totvol = filt_btwl(totvol, lowfq, hghfq)
                                cs = totvol.phase_cog()
                                vol = fshift(totvol, -cs[0], -cs[1], -cs[2])
                                refvols[0] = vol
                                refvols[1] = vol
                                dropImage(refvols[0],os.path.join(outdir, replace("volf%4d%2d.spi"%(jtep, 0),' ','0')), "s")
                                dropImage(refvols[1],os.path.join(outdir, replace("volf%4d%2d.spi"%(jtep, 1),' ','0')), "s")

                            bcast_EMData_to_all(refvols[0], myid, main_node)
                            bcast_EMData_to_all(refvols[1], myid, main_node)

                        else:
                            fscs =[None]*(nrefvol/2)
		 	    for iref in xrange(nrefvol/2):
                                fscfile = os.path.join(outdir, replace("resolution%4d%2d"%(jtep,iref), ' ', '0'))
                                totvol,fscs[iref],refvols[2*iref],refvols[2*iref+1]  = rec3D_MPI_index( dataim, iref, snr, symmetry, mask3D, myid, main_node, fscfile ) 
                                outf.write( '%d reconstructed\n' % iref )
                                outf.flush( )

                                if myid==main_node:
			   	    dropImage(totvol, os.path.join(outdir, replace("vol%4d%2d.spi"%(jtep,iref),' ','0')), "s")
                                    [mean,sigma,fmin,fmax] = Util.infomask( totvol, None, True )
                                    outf.write( 'vol after reconstruction, myid,iref: %d %d %10.3e %10.3e %10.3e %10.3e\n' % ( myid, iref, mean, sigma, fmin, fmax ) )
                                    outf.flush()


			    if(myid == main_node):
                                minlowfq = 1.0
                                minhghfq = 1.0
                                for iref in xrange(nrefvol/2):
                                    fscc = fscs[iref]  
                                    if(fscc[1][0] < 0.5) : fscc[1][0] = 1.0
			            lowfq, hghfq = filt_params(fscc)
                                    outf.write( "iref, lowfq, hghfq: %4d %10.3f %10.3f\n" % (iref, lowfq, hghfq) )
                                    if minlowfq > lowfq :
                                       minlowfq = lowfq

                                    if minhghfq > hghfq :
                                       minhghfq = hghfq

                                outf.write( "minlowfq, minhghfq: %10.3f %10.3f\n" %(minlowfq, minhghfq) )

                                for irefvol in xrange(nrefvol):
				    refvols[irefvol] = filt_btwl(refvols[irefvol], minlowfq, minhghfq)
			            cs = refvols[irefvol].phase_cog()
			            refvols[irefvol] = fshift(refvols[irefvol], -cs[0], -cs[1] -cs[2])
                                    refvols[irefvol] = refine_with_mask(refvols[irefvol])
 
                                    '''
                                    if iref > 0 :
                                        refvols[2*iref  ] = ali_with_mask(refvols[2*iref  ], refvols[0], ali_mask)
                                        outf.write( "%4d aligned\n" % (2*iref) )
                                        outf.flush()

                                        refvols[2*iref+1] = ali_with_mask(refvols[2*iref+1], refvols[1], ali_mask) 
                                        outf.write( "%4d aligned\n" % (2*iref+1) )
                                        outf.flush()
                                    '''

			    	    dropImage(refvols[irefvol],os.path.join(outdir, replace("volf%4d%2d.spi"%(jtep, irefvol),' ','0')), "s")


			    for iref in xrange(nrefvol):
                                bcast_EMData_to_all(refvols[iref], myid, main_node)
                                [mean,sigma,fmin,fmax] = Util.infomask( refvols[iref], None, True )
                                outf.write( 'after bcast, myid,iref: %d %d %10.3e %10.3e %10.3e %10.3e\n' % ( myid, iref, mean, sigma, fmin, fmax ) )

		#  here we should write header info, just in case the program crashes...
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
	if(CTF and data_had_ctf == 0):
		for im in xrange(image_start, image_end): data[im-image_start].set_attr('ctf_applied', 0)
	if(myid == main_node): recv_attr_dict(main_node, stack, dataim, par_str, image_start, image_end, number_of_proc)
	else:                  send_attr_dict(main_node, dataim, par_str, image_start, image_end)

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
	from utilities 		import getImage, model_circle, ce_fit, dropImage,info
	from utilities          import get_arb_params, set_arb_params
	from fundamentals 	import smallprime, window2d, ccf, ramp, fft
	from filter 		import filt_gaussh, filt_tanl
	from string 		import split
	from morphology 	import flcc
	import os
	if os.path.exists(indir)  is False: ERROR("micrograph directory does not exsit", "autowin.py",1)
	else                              : flist=os.listdir(indir)
	if os.path.exists(outdir)         : os.system("rm -rf "+outdir) 
	os.mkdir(outdir)
	t       = EMData()
	e_n     = EMData()
	e       = getImage(noisemic)
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
		if CTF : ctf_params    = get_arb_params(img1, ctf_dicts)
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
		if N_wi == 0 :	ERROR("Number of particles is zero", "autowin.py", 0)
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
			if CTF: set_arb_params(outlist[2], ctf_params, ctf_dicts)
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
	from utilities 		import getImage, model_circle,ce_fit,dropImage,info
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
	if os.path.exists(indir)  is False: ERROR("micrograph directory does not exsit", "autowin.py",1)	
	flist = os.listdir(indir)
	nima          = 0
	mic_name_list = []
	for i, v in enumerate(flist):
		micname                  = os.path.join(indir,v)
		(filename, filextension) = os.path.splitext(v  )
		if(filename[0:len(prefix_of_micrograph)] == prefix_of_micrograph):
			mic_name_list.append(micname)
			nima += 1
	if(myid == int(main_node)): # directory cleaning only performed by main node
		if os.path.exists(outdir):  os.system("rm -rf "+outdir)
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
	e       = getImage(noisemic)
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
		if N_wi == 0 :	ERROR("Number of particles is zero","autowin.py",0)
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


def ihrsr(stack, ref_vol, outdir, maskfile, ir, ou, rs, min_cc_peak, xr, max_x_shift, yr, 
          max_y_shift, max_tilt, ts, delta, an, maxit, CTF, snr, dp, dphi, pixel_size, 
	  rmin, rmax, fract, pol_ang_step, step_a, step_r, sym, user_func_name, datasym):
	if MPI:
		#ihrsr_MPI
		return

	from utilities      import model_circle, dropImage, readSpiderDoc
	from utilities      import getImage, get_input_from_string
	from utilities      import get_arb_params, set_arb_params
	#from filter	    import filt_params, filt_btwl, filt_from_fsc, filt_table, fit_tanh, filt_tanl
	from alignment	    import proj_ali_incore, proj_ali_incore_local, helios
	from statistics     import ccc
	from fundamentals   import cyclic_shift, rot_shift3D
	#from statistics    import fsc_mask
	import os
	import types
	from utilities      import print_begin_msg, print_end_msg, print_msg
	print_begin_msg("ihrsr")

	import user_functions
	user_func = user_functions.factory[user_func_name]

	if CTF :
		from reconstruction import recons3d_4nn_ctf
		ctf_dicts = ["defocus", "Cs", "voltage", "Pixel_size", "amp_contrast", "B_factor", "ctf_applied" ]
		#                     0          1         2             3                 4                   5               6
	else   : from reconstruction import recons3d_4nn

	if os.path.exists(outdir):  os.system('rm -rf '+outdir)
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

	print_msg("Input stack                 : %s\n"%(stack))
	print_msg("Reference volume            : %s\n"%(ref_vol))	
	print_msg("Output directory            : %s\n"%(outdir))
	print_msg("Maskfile                    : %s\n"%(maskfile))
	print_msg("Inner radius                : %i\n"%(first_ring))

	vol     = EMData()
	vol.read_image(ref_vol)
	nx      = vol.get_xsize()
	if (last_ring == -1):	last_ring = nx//2 - 2

	print_msg("Outer radius                              : %i\n"%(last_ring))
	print_msg("Ring step                                 : %i\n"%(rstep))
	print_msg("threshold for the CC peak                 : %f\n"%(min_cc_peak))
	print_msg("X search range                            : %s\n"%(xrng))
	print_msg("threshold for translation in X direction  : %f\n"%(max_x_shift))
	print_msg("Y search range                            : %s\n"%(yrng))
	print_msg("threshold for translation in Y direction  : %f\n"%(max_y_shift))
	print_msg("Translational step                        : %s\n"%(step))
	print_msg("Angular step                              : %s\n"%(delta))
	print_msg("Angular search range                      : %s\n"%(an))
	print_msg("pixel size in Angstoms                    : %f\n"%(pixel_size))
	print_msg("max radius for helical search (in Ang)    : %f\n"%(rmax))
	print_msg("search step for hsearch - angle           : %f\n"%(step_a))
	print_msg("search step for hsearch - rise            : %f\n"%(step_r))
	print_msg("fraction of volume used for helical search: %f\n"%(fract))
	print_msg("initial symmetry - angle                  : %f\n"%(dphi))
	print_msg("initial symmetry - axial rise             : %f\n"%(dp))
	print_msg("Maximum iteration                         : %i\n"%(max_iter))
	print_msg("data with CTF                             : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio                     : %f\n"%(snr))
	print_msg("Symmetry group                            : %s\n\n"%(sym))
	print_msg("symmetry doc file                         : %s\n\n"%(datasym))

	if (maskfile) :
		if (type(maskfile) is types.StringType): mask3D = getImage(maskfile)
		else                                  : mask3D = maskfile
	else          :   mask3D = model_circle(last_ring, nx, nx, nx)
	mask = model_circle(last_ring, nx, nx)
	#dropImage(vol, os.path.join(outdir,"ref_vol00.hdf"))
	ref_a = "s"
	names_params = ["psi", "s2x", "s2y", "peak"]
	data = EMData.read_images(stack)
	for im in xrange(nima):
		data[im].set_attr('ID', im)
		if(CTF):
			ctf_params = get_arb_params(data[im], ctf_dicts)
			if(im == 0): data_had_ctf = ctf_params[6]
			if(ctf_params[6] == 0):
				st = Util.infomask(data[im], mask, False)
				data[im] -= st[0]
				from filter import filt_ctf
				data[im] = filt_ctf(data[im], ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
				data[im].set_attr('ctf_applied', 1)

	finfo = None#open("desperado", 'w')
	# do the projection matching
	dropImage(vol, os.path.join(outdir, "aligned0000.spi"), "s")
	for N_step in xrange(lstp):
 		# max_iter=30
		for Iter in xrange(max_iter):
			print_msg("ITERATION #%3d\n"%(N_step*max_iter + Iter+1))
			#from filter import filt_gaussinv
			#vol = filt_gaussinv(vol, 0.175, True)
			#dropImage(vol, os.path.join(outdir, replace("vhl%4d.spi"%(N_step*max_iter+Iter+1),' ','0')), "s")

			if(an[N_step] == -1):  proj_ali_incore(vol, mask3D, data, first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], ref_a, sym, finfo = finfo)
			else:	               proj_ali_incore_local(vol, mask3D, data, first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], an[N_step], ref_a, sym)
			
			# exclusion code
			# psi, s2x s2y, peak
			for i in xrange(nima):
				paramali = get_arb_params(data[i], names_params)
				active = 1
				# peak
				if(paramali[3] < min_cc_peak): active = 0
				# s2x
				elif(abs(paramali[1]) > max_x_shift): active = 0
				# s2y
				elif(abs(paramali[2]) > max_y_shift): active = 0
				# psi correct value should be 90 +/- max_tilt, or 270 +/- max_tilt!
				elif(abs(paramali[0]-270.0) >max_tilt and paramali[0] >180.0): active = 0
				elif(abs(paramali[0]-90.0) >max_tilt and paramali[0] <180.0): active = 0
				data[i].set_attr_dict({'active':active})
				#print the alignment parameters into the LOG file!
				print_msg("image status                : %i\n"%(active))
				print_msg("image psi                   : %i\n"%(paramali[0]))
				print_msg("image s2x                   : %i\n"%(paramali[1]))
				print_msg("image s2y                   : %i\n"%(paramali[2]))
				print_msg("image peak                  : %i\n"%(paramali[3]))
			#  3D stuff
			#  I removed symmetry, by default the volume is not symmetrized
			#  calculate new and improved 3D
			if(CTF): vol = recons3d_4nn_ctf(data, range(nima), snr, npad = 2)
			else:	 vol = recons3d_4nn(data, range(nima), npad = 2)
			# store the reference volume
			vof  = os.path.join(outdir, "unsymmetrized%04d.spi"%(N_step*max_iter+Iter+1))
			dropImage(vol, vof, "s")
			'''
			vofs = os.path.join(outdir, "unsymmetrized%04d"%(N_step*max_iter+Iter+1))
			vofo = os.path.join(outdir, "unsymmetrized%04d_byte_swapped"%(N_step*max_iter+Iter+1))
			vofq = os.path.join(outdir, "symmetrized%04d.spi"%(N_step*max_iter+Iter+1))
			dropImage(vol, os.path.join(outdir, "raw%04d.spi"%(N_step*max_iter+Iter+1)), "s")
			dropImage(rot_shift3D(vol,90.,-90.,270.), vof, "s")

			# if center == 1:
			# from utilities import rotate_3D_shift
			# rotate_3D_shift(data, cs)
			# dropImage(vol, os.path.join(outdir, "volf%04d.hdf"%(N_step*max_iter+Iter+1)))
			#qtp = readSpiderDoc(datasym)
			outtemp = open("qtdwe", "w")
			outtemp.write("cp to opend\n")
			outtemp.write(vofs)
			outtemp.write("\n")
			outtemp.write(vofo)
			outtemp.write("\n")
			outtemp.write("en d\n")
			outtemp.close()
			os.system("ls -l    qtdwe")
			os.system("cat    qtdwe")
			os.system("/usr/local/bin/spider_linux_mpfftw_opt64   spi   <qtdwe")
			os.system("hsearch_lorentz  "+vofo+".spi  "+datasym+"   2.380000   0.00000  55.00000   0.10000   0.10000")
			os.system("himpose  "+vofo+".spi  "+datasym+"  "+vofq+" 2.380000   0.00000  55.00000 ")
			vol = getImage(vofq)
			'''
			vol, dp, dphi = helios(vol, pixel_size, dp, dphi, fract, rmax)
			print_msg("new delta z and delta phi      : %s,    %s\n\n"%(dp,dphi))
			previous_vol=getImage(os.path.join(outdir, "aligned%04d.spi"%(N_step*max_iter+Iter)))
			#Here we align the volume with the previous one to keep the polarity fixed.
			#peakmax - any large negative value for the CC peak to start from
			#360 - literaly 360 degrees
			peakmax=[-1000000.0]*3
			for i in xrange(0, 360, pol_ang_step):
				vtm = rot_shift3D(vol,float(i))
				s_r = int(dp/pixel_size)
				#s_r - search range, should be at least +/- 1 pixel, if s_r < 2 then we set it to 2
				if(int(s_r) < 2): s_r = 2 
				for j in xrange(s_r):
					zs = j-s_r//2
					peaks = ccc(cyclic_shift(vtm, 0, 0, zs), previous_vol)
					if(peaks>peakmax[0]):  peakmax=[peaks, i, zs]
			vol = cyclic_shift(rot_shift3D(vol, peakmax[1]),  0, 0, peakmax[2])
			dropImage(vol, os.path.join(outdir, "aligned%04d.spi"%(N_step*max_iter+Iter+1)), "s")
	#  here we  write header info
	if(CTF and data_had_ctf == 0):
		for im in xrange(nima): data[im].set_attr('ctf_applied', 0)
	from utilities import write_headers
	write_headers( stack, data, range(nima))
	print_end_msg("ihrsr")

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

	from utilities 		import getImage, dropImage
	from fundamentals 	import smallprime, window2d, resample, image_decimate
	from filter 		import filt_btwl
	import types
	import os
	if os.path.exists(indir) is False: ERROR("Input directory doesn't exist","copyfromtif",1)
	else                             : flist          = os.listdir(indir)
	if(type(outdir)          is types.StringType):
		if os.path.exists(outdir) is True:   os.system("rm -rf " +outdir)
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
			e = getImage(tifname)
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
			if output_extension == "spi": dropImage(e1,f_micname,"s")
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
	from utilities 		import getImage, dropImage
	from fundamentals 	import smallprime, window2d, resample, image_decimate
	from filter 		import filt_btwl
	from random     	import randint
	import types
	import os
	import sys
	from mpi 	    import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from mpi 	    import mpi_barrier, mpi_bcast
	from mpi 	    import MPI_INT

	if os.path.exists(indir) is False: ERROR("Input directory doesn't exist","copyfromtif",1)
	else                             : flist = os.listdir(indir)
	if(type(outdir)          is types.StringType):
		if os.path.exists(outdir):   
			os.system("rm -rf " +outdir)
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
	if nima < 1:    ERROR("No micrograph is found, check either directory or prefix of micrographs is correctly given","pw2sp",1)
       
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
		e = getImage(filename)
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
		if  gridding : e1    = resample(e, scaling_ratio, 1) # resample will pad image to four times 
		else         : e1    = image_decimate(e, scaling_ratio, 1)
		f_micrograph         = "micrograph_"+tifname[len(indir)+1:]+"."+ output_extension
		f_micname            = os.path.join(outdir, f_micrograph)
		e1.write_image(f_micname)

def cpy(ins, ous):
	nima = EMUtil.get_image_count(ins)
	data = EMData()
	from utilities import file_type
	iextension = file_type(ins)
	oextension = file_type(ous)
			
	if(nima == 1 and oextension == "spi"):
		data.read_image(ins)
		data.write_image(ous, 0, EMUtil.ImageType.IMAGE_SINGLE_SPIDER)

	elif(iextension == "bdb" and oextension == "bdb"):

		DB = db_open_dict(ins)
		OB = db_open_dict(ous)
		for i in range(nima):
			OB[i] = DB[i]
		DB.close()
		OB.close()

	elif(iextension == "bdb"):

		DB = db_open_dict(ins)
		for i in range(nima):
			a = DB[i]
			a.write_image(ous, i)
		DB.close()

	elif(oextension == "bdb"):

		DB = db_open_dict(ous)
		for i in range(nima):
			a = EMData()
			a.read_image(ins, i)
			DB[i] = a
		DB.close()
		
	else:
		for im in xrange(nima):
			data.read_image(ins, im)
			data.write_image(ous, im)

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

def project3d(volume, stack, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None ):
# 2D multi-reference alignment using rotational ccf in polar coords and quadratic interpolation
	from projection    import   prgs, prep_vol
	from utilities     import   even_angles, read_txt_col
	from string        import   split
	import os
	import types

	if listagls is None:
		angles = even_angles(delta, symmetry = symmetry, method = method, phiEqpsi = phiEqpsi)
	elif(type(listagls) is types.StringType):
		angles = read_txt_col(listagls, "", "")
	else:
		angles = listagls


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
		volft, kb = prep_vol(vol)
	else:
		if(mask):
			if(type(mask) is types.StringType):
				maski = EMData()
				maski.read_image(volume)
				Util.mul_img(vol, maski)
				del maski
			else:
				Util.mul_img(vol, mask)
		volft, kb = prep_vol(volume)


	if(type(stack) is types.StringType):
		Disk = True
		os.system("rm -f  "+stack)	
	else:
		out = []
		Disk = False
	
	s2x=0
	s2y=0
	for i in xrange(len(angles)):
		proj = prgs(volft, kb, [angles[i][0], angles[i][1], angles[i][2], 0.0, 0.0])
		proj.set_attr_dict({'phi':angles[i][0], 'theta':angles[i][1], 'psi':angles[i][2], 's2x':0.0, 's2y':0.0, 'active':1})
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

	from utilities    import getImage,dropImage, model_circle, info
	from fundamentals import welch_pw2, ro_textfile
	import sys
	import os
	import types
	if os.path.exists(indir) is False: ERROR("Input directory doesn't exist","pw2sp",1)
	else				 : flist    = os.listdir(indir)
	if(type(outdir)          is types.StringType):
		if os.path.exists(outdir) is True: os.system("rm -rf " +outdir)
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
			e         = getImage(micname)
 			pw2       = welch_pw2(e,int(w),int(xo),int(yo),int(xd),int(yd))
			pw2_mask  = pw2*mask
			pw2name   = os.path.join(outdir,pw2file)
			dropImage(pw2_mask, pw2name)
			rotxtname = os.path.join(outdir,roofile)
			ro_textfile(pw2, rotxtname)
	if ncount < 1: 	ERROR("No micrograph is found, check either directory or prefix of micrographs is correctly given","pw2sp",1)
 
def pw2sp_MPI(indir, outdir, w =256, xo =50, yo = 50, xd = 0, yd = 0, r = 0, prefix_of_micrograph="micrograph"):
	""" 
		Purpose: 
		Calculate power spectra of a list of micrographs in a given directory using Welch's periodogram
		The input options enable one selects area in micrographs to calculate overlapped periodogram.
	"""
	from utilities    	import getImage,dropImage, model_circle, info
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
	if (myid == int(main_node)):  # only main node do cleaning & creating jobs
		if os.path.exists(outdir): os.system("rm -rf " +outdir)
		os.mkdir(outdir)
	mpi_barrier(MPI_COMM_WORLD)
	# get the total micrograph number 
	if os.path.exists(indir) is False: ERROR("Input directory doesn't exist","pw2sp_MPI",1)
	else:	                           flist = os.listdir(indir)
	nima          = 0
	mic_name_list = []
	for i, v in enumerate(flist):
		micname                  = os.path.join(indir,v)
		(filename, filextension) = os.path.splitext(v)
		if(filename[0:len(prefix_of_micrograph)] == prefix_of_micrograph):
			mic_name_list.append(micname)
			nima += 1
	if nima < 1: 	ERROR("No micrograph is found, check spelling of either directory or micrographs","pw2sp_MPI",1)	

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
		e         = getImage(filename)
 		pw2       = welch_pw2(e,int(w),int(xo),int(yo),int(xd),int(yd))
		pw2_mask  = pw2*mask
		pw2name   = os.path.join(outdir,pw2file)
		dropImage(pw2_mask, pw2name)
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
	e_n    = getImage(noise)
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
	new_params = amoeba(new_params, scale, ali_vol_func, 1.e-1, 1.e-1, 500, data)
        vol = rot_shift3D(vol, new_params[0][0], new_params[0][1], new_params[0][2], new_params[0][3], new_params[0][4], new_params[0][5],1.0)
        return vol
        
def ali_vol_3(vol, refv, ang_scale, shift_scale, radius=None, discrepancy = "ccc"):
	#rotation and shift
	from alignment    import ali_vol_func
	from utilities    import model_circle

	nx = refv.get_xsize()
	ny = refv.get_ysize()
	nz = refv.get_zsize()
	if(radius != None): mask = model_circle(radius, nx, ny, nz)
	else:               mask = model_circle(float(min(nx, ny, nz)//2-2), nx, ny, nz)

	scale = [ang_scale, ang_scale, ang_scale, shift_scale, shift_scale, shift_scale]
	data=[vol, refv, mask]
	new_params = [0.0]*6
	opt_params,funval,niter = amoeba(new_params, scale, ali_vol_func, 1.e-1, 1.e-1, 500, data)
        return opt_params
 
def ali_vol(vol, refv, ang_scale, shift_scale, radius=None, discrepancy = "ccc"):
	#rotation and shift
	from alignment    import ali_vol_func
	from utilities    import get_im, model_circle, get_params3D, set_params3D
	from utilities    import amoeba, compose_transform3
	from fundamentals import rot_shift3D
	ref = get_im(refv)
	nx = ref.get_xsize()
	ny = ref.get_ysize()
	nz = ref.get_zsize()
	if(radius != None):    mask = model_circle(radius, nx, ny, nz)
	else:                  mask = model_circle(float(min(nx, ny, nz)//2-2), nx, ny, nz)

	#names_params = ["phi", "theta", "psi", "s3x", "s3y", "s3z", "scale"]
	params = get_params3D(ref)
	print  " params of the reference volume",params
	ref = rot_shift3D(ref, params[0], params[1], params[2], params[3], params[4], params[5], params[7])

	e = get_im(vol)
	params = get_params3D(e)
	print  params
	e = rot_shift3D(e, params[0], params[1], params[2], params[3], params[4], params[5], params[7])

	e = get_im(vol)
	params = get_arb_params(e, names_params)
	print  " input params ",params
	data=[e, ref, mask, params, discrepancy]
	new_params = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	new_params = amoeba(new_params, [ang_scale, ang_scale, ang_scale, shift_scale, shift_scale, shift_scale], ali_vol_func, 1.e-1, 1.e-1, 500, data)
	cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale= compose_transform3(params[0], params[1], params[2], params[3], params[4], params[5], params[7], new_params[0][0], new_params[0][1], new_params[0][2], new_params[0][3], new_params[0][4], new_params[0][5],1.0)
	print  " new params ", cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale, new_params[1]
	set_params3D(e, [cphi, ctheta, cpsi, cs2x, cs2y, cs2z, 0, cscale])
	from utilities import write_headers
	write_headers( vol, [e], [0])

def ali_vol_rotate(vol, refv, ang_scale, radius=None, discrepancy = "ccc"):
	#rotation 
	from alignment    import ali_vol_func_rotate
	from utilities    import get_im, model_circle, get_params3D, set_params3D
	from utilities    import amoeba, compose_transform3
	from fundamentals import rot_shift3D
	ref = get_im(refv)
	nx = ref.get_xsize()
	ny = ref.get_ysize()
	nz = ref.get_zsize()
	if(radius != None):    mask = model_circle(radius, nx, ny, nz)
	else:                  mask = model_circle(float(min(nx, ny, nz)//2-2), nx, ny, nz)

	#names_params = ["phi", "theta", "psi", "s3x", "s3y", "s3z", "scale"]
	params = get_params3D(ref)
	print  " params of the reference volume",params
	ref = rot_shift3D(ref, params[0], params[1], params[2], params[3], params[4], params[5], params[7])

	e = get_im(vol)
	params = get_params3D(e)
	print  " input params ",params
	data=[e, ref, mask, params, discrepancy]
	new_params = [0.0, 0.0, 0.0]
	new_params = amoeba(new_params, [ang_scale, ang_scale, ang_scale], ali_vol_func_rotate, 1.e-1, 1.e-1, 500, data)
	cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale= compose_transform3(params[0], params[1], params[2], params[3], params[4], params[5], params[7], new_params[0][0], new_params[0][1], new_params[0][2],0.0,0.0,0.0,1.0)
	print  " new params ", cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale, new_params[1]
	set_params3D(e, [cphi, ctheta, cpsi, cs2x, cs2y, cs2z, 0, cscale])
	from utilities import write_headers
	write_headers( vol, [e], [0])

def ali_vol_shift(vol, refv, shift_scale, radius=None, discrepancy = "ccc"):
	# shift
	from alignment    import ali_vol_func_shift
	from utilities    import get_im, model_circle, get_params3D, set_params3D
	from utilities    import amoeba, compose_transform3
	from fundamentals import rot_shift3D
	ref = get_im(refv)
	nx = ref.get_xsize()
	ny = ref.get_ysize()
	nz = ref.get_zsize()
	if(radius != None):    mask = model_circle(radius, nx, ny, nz)
	else:                  mask = model_circle(float(min(nx, ny, nz)//2-2), nx, ny, nz)

	#names_params = ["phi", "theta", "psi", "s3x", "s3y", "s3z", "scale"]
	params = get_params3D(ref)
	print  " params of the reference volume",params
	ref = rot_shift3D(ref, params[0], params[1], params[2], params[3], params[4], params[5], params[7])

	e = get_im(vol)
	params = get_params3D(e)
	print  " input params ",params
	data=[e, ref, mask, params, discrepancy]
	new_params = [0.0, 0.0, 0.0]
	new_params = amoeba(new_params, [shift_scale, shift_scale, shift_scale], ali_vol_func_shift, 1.e-1, 1.e-1, 500, data)
	cphi, ctheta, cpsi, cs3x, cs3y, cs3z, cscale= compose_transform3(params[0], params[1], params[2], params[3], params[4], params[5], params[7], 0.0,0.0,0.0, new_params[0][0], new_params[0][1], new_params[0][2],1.0)
	print  " new params ", cphi, ctheta, cpsi, cs3x, cs3y, cs3z, cscale, new_params[1]
	set_params3D(e, [cphi, ctheta, cpsi, cs3x, cs3y, cs3z, 0, cscale])
	from utilities import write_headers
	write_headers( vol, [e], [0])

def ali_vol_scale(vol, refv, ang_scale, shift_scale, mag_scale, radius=None, discrepancy = "ccc"):
	# rotation shift and scale
	from alignment    import ali_vol_func_scale
	from utilities    import get_im, model_circle, get_params3D, set_params3D
	from utilities    import amoeba, compose_transform3
	from fundamentals import rot_shift3D
	ref = get_im(refv)
	nx = ref.get_xsize()
	ny = ref.get_ysize()
	nz = ref.get_zsize()
	if(radius != None):    mask = model_circle(radius, nx, ny, nz)
	else:                  mask = model_circle(float(min(nx, ny, nz)//2-2), nx, ny, nz)

	#names_params = ["phi", "theta", "psi", "s3x", "s3y", "s3z", "scale"]
	params = get_params3D(ref)
	print  " params of the reference volume",params
	ref = rot_shift3D(ref, params[0], params[1], params[2], params[3], params[4], params[5], params[7])

	e = get_im(vol)
	params = get_params3D(e)
	print  " input params ",params
	data=[e, ref, mask, params, discrepancy]
	new_params = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
	new_params = amoeba(new_params, [ang_scale, ang_scale, ang_scale, shift_scale, shift_scale, shift_scale, mag_scale], ali_vol_func_scale, 1.e-1, 1.e-1, 500, data)
	cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale= compose_transform3(params[0], params[1], params[2], params[3], params[4], params[5], params[7], new_params[0][0], new_params[0][1], new_params[0][2], new_params[0][3], new_params[0][4], new_params[0][5], new_params[0][6])
	print  " new params ", cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale, new_params[1]
	set_params3D(e, [cphi, ctheta, cpsi, cs2x, cs2y, cs2z, 0, cscale])
	from utilities import write_headers
	write_headers( vol, [e], [0])

def ali_vol_only_scale(vol, refv, mag_scale, radius=None, discrepancy = "ccc"):
	# scale
	from alignment    import ali_vol_func_only_scale
	from utilities    import get_im, model_circle, get_params3D, set_params3D
	from utilities    import amoeba, compose_transform3
	from fundamentals import rot_shift3D
	ref = get_im(refv)
	nx = ref.get_xsize()
	ny = ref.get_ysize()
	nz = ref.get_zsize()
	if(radius != None):    mask = model_circle(radius, nx, ny, nz)
	else:                  mask = model_circle(float(min(nx, ny, nz)//2-2), nx, ny, nz)


	#names_params = ["phi", "theta", "psi", "s3x", "s3y", "s3z", "scale"]
	params = get_params3D(ref)
	print  " params of the reference volume",params
	ref = rot_shift3D(ref, params[0], params[1], params[2], params[3], params[4], params[5], params[6])

	e = get_im(vol)
	params = get_params3D(e)
	print  " input params ",params
	data=[e, ref, mask, params, discrepancy]
	new_params = [1.0]
	new_params = amoeba(new_params, [mag_scale], ali_vol_func_only_scale, 1.e-1, 1.e-1, 500, data)
	cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale= compose_transform3(params[0], params[1], params[2], params[3], params[4], params[5], params[7], 0.0,0.0,0.0,0.0,0.0,0.0, new_params[0][0])
	print  " new params ", cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale, new_params[1]
	set_params3D(e, [cphi, ctheta, cpsi, cs2x, cs2y, cs2z, 0, cscale])
	from utilities import write_headers
	write_headers( vol, [e], [0])

def rot_sym(infile, outfile, sym_gp="d4", \
			radius=None, phi=0, theta=0, psi=0, phirange=20, thetarange=20, psirange=20, ftolerance=1.e-4, xtolerance=1.e-4):
	
	from alignment import find_symm
	from utilities import dropImage, model_circle, sym_vol
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

	dropImage(sym, outfile)
	
def transform2d(stack_data, stack_data_ali):
# apply 2D alignment parameters stored in the header of the input stack file using gridding interpolation and create an output stack file
	from fundamentals	import rot_shift2D
	from utilities 		import set_arb_params, combine_params2
	from utilities      import print_begin_msg, print_end_msg, print_msg
	import os
	
	print_begin_msg("transform2d")	
	print_msg("Input stack                 : %s\n"%(stack_data))
	print_msg("Output stack                : %s\n\n"%(stack_data_ali))

	if os.path.exists(stack_data_ali): os.system("rm -f "+stack_data_ali)

	attributes = ['phi', 'theta', 'psi', 'alpha', 'sx', 'sy', 'mirror', 'nclass', 'assign', 's2x', 's2y', 'scale']
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
		# apply params to the image
		if(params[3] != 0.0 or params[4] != 0.0 or params[5] != 0.0):
			temp = rot_shift2D(data, params[3], params[4], params[5])
			params[3]=params[4]=params[5] = 0.0
		else:
			temp = data.copy()
		if  params[6]:
			temp.process_inplace("mirror", {"axis":'x'})
			params[6] = 0
		set_arb_params(temp, params, attributes)
		temp.write_image(stack_data_ali, im)
	print_end_msg("transform2d")

def recons3d_n(prj_stack, pid_list, vol_stack, CTF=False, snr=1.0, sign=1, npad=4, sym="c1", verbose=0, MPI=False):
	if MPI:
		recons3d_n_MPI(prj_stack, pid_list, vol_stack, CTF, snr, 1, npad, sym, verbose)
		return

	from reconstruction import recons3d_4nn_ctf, recons3d_4nn
	from utilities import dropImage
	from utilities import print_begin_msg, print_end_msg, print_msg

	print_begin_msg("recons3d_n")
	print_msg("Input stack                 : %s\n"%(prj_stack))
	print_msg("Output volume               : %s\n"%(vol_stack))
	print_msg("Padding factor              : %i\n"%(npad))
	print_msg("data with CTF               : %s\n"%(CTF))
	print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
	print_msg("CTF sign                    : %i\n"%(sign))
	print_msg("Symmetry group              : %s\n\n"%(sym))

	if CTF: vol = recons3d_4nn_ctf(prj_stack, pid_list, snr, 1, sym, verbose)
	else:   vol = recons3d_4nn(prj_stack,  pid_list, sym, npad)
	dropImage(vol, vol_stack)
	print_end_msg("recons3d_n")

def recons3d_n_MPI(prj_stack, pid_list, vol_stack, ctf, snr, sign, npad, sym, verbose):
	from reconstruction import recons3d_4nn_ctf_MPI, recons3d_4nn_MPI
	from utilities import get_im
	from string import replace
	from time import time
	from mpi 	    import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD

	myid  = mpi_comm_rank(MPI_COMM_WORLD)
	nproc = mpi_comm_size(MPI_COMM_WORLD)
	time_start = time()

	if verbose==0:
		info = None
	else:
		infofile = "progress%4d.txt" % (myid+1)
		infofile = replace(infofile, ' ', '0')
		info = open( infofile, 'w' )

	nimage = len(pid_list)

	#nimage_per_node = nimage/nproc
	#image_start = myid * nimage_per_node
	#if myid == nproc-1 : image_end = nimage
	#else:                image_end = image_start + nimage_per_node
	image_start, image_end = MPI_start_end(nimage, nproc, myid)

	# Berlin
	#from utilities import read_txt_col, set_arb_params
	#prm = read_txt_col("coco","u")
	#par_str=["phi", "theta", "psi", "s2x", "s2y"]

	prjlist = []
	for i in range(image_start,image_end):
		prj = get_im( prj_stack, pid_list[i] )

		#set_arb_params(prj, prm[pid_list[i]], par_str)

		prjlist.append( prj )
		if not(info is None): info.write( "%4d read\n" % i )
	#del prj
	#del prm

	if ctf: vol = recons3d_4nn_ctf_MPI(myid, prjlist, snr, sign, sym, info)
	else:	vol = recons3d_4nn_MPI(myid, prjlist, sym, info)

	if myid == 0 :
		vol.write_image(vol_stack)
		if not(info is None):
			info.write( "result wrote to " + vol_stack + "\n")
			info.write( "Total time: %10.3f\n" % (time()-time_start) )
			info.flush()

def recons3d_f(prj_stack, vol_stack, fsc_file, mask=None, CTF=True, snr=1.0, sym="c1", verbose=1, MPI=False):
	if MPI:
		recons3d_f_MPI(prj_stack, vol_stack, fsc_file, mask, CTF, snr, sym, verbose)
		return

	nimage = EMUtil.get_image_count( prj_stack )

	from reconstruction import recons3d_4nn_ctf, recons3d_4nn
	from statistics     import fsc_mask
	if CTF:
		volodd = recons3d_4nn_ctf(prj_stack, range(0, nimage, 2), snr, 1, sym, verbose)
		voleve = recons3d_4nn_ctf(prj_stack, range(1, nimage, 2), snr, 1, sym, verbose)
		volall = recons3d_4nn_ctf(prj_stack, range(nimage),       snr, 1, sym, verbose)
	else:
		volodd = recons3d_4nn(prj_stack, range(0, nimage, 2), sym)
		voleve = recons3d_4nn(prj_stack, range(1, nimage, 2), sym)
		volall = recons3d_4nn(prj_stack, range(nimage),       sym)
	volall.write_image(vol_stack)
	t = fsc_mask( volodd, voleve, mask, filename=fsc_file)

def recons3d_f_MPI(prj_stack, vol_stack, fsc_file, mask, CTF=True, snr=1.0, sym="c1", verbose=1):

	from mpi import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from utilities import get_im
	nproc = mpi_comm_size( MPI_COMM_WORLD )
	myid  = mpi_comm_rank( MPI_COMM_WORLD )

	img_number     = EMUtil.get_image_count( prj_stack )

	#img_per_node   = img_number/nproc
	#img_node_start = img_per_node*myid
	#if myid == nproc-1:  img_node_end = img_number
	#else:                 img_node_end = img_node_start + img_per_node
	img_node_start, img_node_end = MPI_start_end(img_number, nproc, myid)

	imgdata = []
	for i in xrange(img_node_start, img_node_end):
		img = get_im(prj_stack, i)
		imgdata.append(img)
	del img
	print  "  DATA  LOADED  ",myid
	if CTF:
		from reconstruction import rec3D_MPI
		odd_start = img_node_start % 2
		eve_start = (odd_start+1)%2

		if verbose==0:
			info = None
		else:
			from string import replace
			infofile = "progress%4d.txt" % (myid+1)
			infofile = replace(infofile, ' ', '0')
			info = open( infofile, 'w' )

		vol,fsc = rec3D_MPI(imgdata, snr, sym, mask, fsc_file, myid, 0, 1.0, odd_start, eve_start, info)
	else :
		from reconstruction import rec3D_MPI_noCTF
		vol,fsc = rec3D_MPI_noCTF(imgdata, sym, mask, fsc_file, myid)
	if myid == 0: vol.write_image(vol_stack)

def pca( input_stack, output_stack, imglist, nfile, subavg, mask_radius, nvec, type="out_of_core", maskfile="",verbose=False ) :
	from utilities import getImage, get_im, model_circle, model_blank
	from EMAN2 import Analyzers

	if mask_radius > 0 and maskfile !="":
		print "Error: mask radius and mask file cannot be used at the same time"
		return

	if mask_radius >0:
		if(verbose):
			print 'maskfile: ', maskfile
		assert maskfile==""
		if(verbose):
			print "Using spherical mask, rad=", mask_radius
		data = EMData()
		if nfile == 1:
			data.read_image( input_stack, 0, True)
		else:
			data.read_image( input_stack + "0000.hdf", 0, True)
		mask = model_circle(mask_radius, data.get_xsize(), data.get_ysize(), data.get_zsize())
	elif(maskfile!="") :
		if(verbose):
			print "Using mask: ", maskfile
		mask = getImage( maskfile )
	else:
		data = EMData()
		data.read_image( input_stack, 0, True)
		mask = model_blank(data.get_xsize(), data.get_ysize(), data.get_zsize(), bckg=1.0)

	if subavg != "":
		if(verbose):
			print "Subtracting ", subavg, " from each image"
		avg = getImage( subavg )
	else:
		avg = None
	if type == "out_of_core" : ana = Analyzers.get( "pca_large", {"mask":mask, "nvec":nvec} )
	else :                     ana = Analyzers.get( "pca", {"mask":mask, "nvec":nvec} )

	totimg = 0
	for j in xrange(nfile):
		
		if nfile==1:
			curt_input_stack = input_stack
		else:
			curt_input_stack = input_stack + ("%04d.hdf" % j)
			nimage = EMUtil.get_image_count( curt_input_stack )
			imglist = xrange(nimage)

		if(verbose):
			print "loading file ", curt_input_stack
		for i in imglist:
			data = get_im( curt_input_stack, i)
			if not(avg is None):
				# the folllowing is incerrect!
				#ERROR("DO NOT USE THIS OPTION","pca",1)
				[avg1,sigma,fmin,fmax]=Util.infomask(data, None,True)
				data -= avg
				[avg2,sigma,fmin,fmax]=Util.infomask(data, None,True)
			ana.insert_image( data )
			if(verbose):
				 print "Inserting image %4d" % totimg
			totimg += 1

	vecs = ana.analyze()
	iout = 0
	for vec in vecs:
		vec.write_image( output_stack, iout)
		iout = iout + 1

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
			ctf_params = get_arb_params(ima, parnames)
			ali = get_arb_params(ima, pali)
			ima    = rot_shift2D(ima, ali[0], ali[1], ali[2])
			if  ali[3]:  ima.process_inplace("mirror",{"axis":'x'})
			if avg:
				st = Util.infomask(ima, mask, False)
				ima -= st[0]
			oc = filt_ctf(fft(pad(ima, nx2, ny2, background = 0.0)), ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5])
			Util.add_img(ave, oc)
			Util.add_img2(ctf_2_sum, ctf_img(nx2, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5], ny = ny2, nz = 1))
		Util.div_filter(ave, ctf_2_sum)
		for i in xrange(n):
			ima = EMData()
			ima.read_image(input_stack, i)
			ctf_params = get_arb_params(ima, parnames)
			ali = get_arb_params(ima, pali)
			ima    = rot_shift2D(ima, ali[0], ali[1], ali[2])
			if  ali[3]:  ima.process_inplace("mirror",{"axis":'x'})
			if avg:
				st = Util.infomask(ima, mask, False)
				ima -= st[0]
			oc = filt_ctf(ave, ctf_params[1], ctf_params[3], ctf_params[2], ctf_params[0], ctf_params[4], ctf_params[5], pad= True)
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
			ima    = rot_shift2D(ima, ali[0], ali[1], ali[2])
			if  ali[3]:  ima.process_inplace("mirror",{"axis":'x'})
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
			ima    = rot_shift2D(ima, ali[0], ali[1], ali[2])
			if  ali[3]:  ima.process_inplace("mirror",{"axis":'x'})
			if avg:
				st = Util.infomask(ima, mask, False)
				ima -= st[0]
			Util.sub_img(ima, ave)
			set_arb_params(ima, [0.0,0.0,0.0,0], pali)
			ima.write_image(output_stack, i)

def varimax(input_stack, imglist, output_stack, mask_radius, verbose ) :
	from utilities import getImage, model_circle
	from EMAN2 import Analyzers

	data = getImage( input_stack )
	mask = model_circle( mask_radius, data.get_xsize(), data.get_ysize(), data.get_zsize() )

	ana = Analyzers.get( "varimax", {"mask":mask} )

	for i in imglist:
		data = EMData()
		data.read_image( input_stack, i)
		ana.insert_image( data )
		#print "Inserting image %4d" % i

	vecs = ana.analyze()

	iout = 0
	for vec in vecs:
		vec.write_image( output_stack, iout)
		iout = iout + 1

def bootstrap_genbuf(prj_stack, buf_prefix, verbose, MPI=False):
	from EMAN2 import file_store
	import string
	from mpi import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD, mpi_init

	npad = 4

	if(MPI):
		sys.argv = mpi_init(len(sys.argv),sys.argv)
		size = mpi_comm_size(MPI_COMM_WORLD)
		myid = mpi_comm_rank(MPI_COMM_WORLD)
	else:
		size = 1
		myid = 0

	store = file_store(buf_prefix, npad, 1)

	if verbose != 0 :
		mystatus = "genbuf%4d.txt" % ( myid )
		mystatus = string.replace( mystatus, ' ', '0' )
		output = open( mystatus, "w" )

	nimage = EMUtil.get_image_count( prj_stack )
	for i in xrange(nimage):
		proj = EMData()
		proj.read_image( prj_stack, i )
		store.add_image( proj )

		if( verbose !=0 and (i+1) % 100 == 0 ) :
			output.write( "proj %4d done\n" % (i+1) )
			output.flush()

	if verbose != 0:
		output.write( "proj %4d done\n" % nimage )
		output.flush()
 
def bootstrap_run(prj_stack, media, vol_prefix, nvol, snr, sym, verbose, MPI=False):

	import string
	from mpi import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD, mpi_init

	if MPI:
		sys.argv = mpi_init(len(sys.argv),sys.argv)
		size = mpi_comm_size(MPI_COMM_WORLD)
		myid = mpi_comm_rank(MPI_COMM_WORLD)
	else:
		size = 1
		myid = 0

	myvolume_file = "%s%4d.hdf" % ( vol_prefix, myid )
	myvolume_file = string.replace( myvolume_file, ' ', '0' )
	if verbose != 0 :
		mystatus_file = "status%4d.inf" % (myid)
		mystatus_file = string.replace( mystatus_file, ' ', '0' )
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
	from reconstruction import bootstrap_nnctf
	bootstrap_nnctf( prj_stack, myvolume_file, list_proj, mynvol, media, npad, sym, mystatus, snr, sign)
	
def params_2D_to_3D(stack):
	from utilities import params_2D_3D, print_begin_msg, print_end_msg, print_msg
	
	print_begin_msg("params_2D_to_3D")
	print_msg("Input stack                 : %s\n\n"%(stack))

	nima = EMUtil.get_image_count(stack)
	ima = EMData()
	for im in xrange(nima):
		ima.read_image(stack, im, True)
		alpha  = ima.get_attr('alpha')
		sx     = ima.get_attr('sx')
		sy     = ima.get_attr('sy')
		mirror = ima.get_attr('mirror')
		phi, theta, psi, s2x, s2y = params_2D_3D(alpha, sx, sy, mirror)
		ima.set_attr_dict({'phi':phi, 'theta':theta, 'psi':psi, 's2x':s2x, 's2y':s2y})
		ima.write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)	
	print_end_msg("params_2D_to_3D")
	
def params_3D_to_2D(stack):
	from utilities import params_3D_2D, print_begin_msg, print_end_msg, print_msg
	
	print_begin_msg("params_3D_to_2D")
	print_msg("Input stack                 : %s\n\n"%(stack))

	nima = EMUtil.get_image_count(stack)
	ima = EMData()
	for im in xrange(nima):
		ima.read_image(stack, im, True)
		s2x   = ima.get_attr('s2x')
		s2y   = ima.get_attr('s2y')
		phi   = ima.get_attr('phi')
		theta = ima.get_attr('theta')
		psi   = ima.get_attr('psi')
		alpha, sx, sy, mirror = params_3D_2D(phi, theta, psi, s2x, s2y)
		ima.set_attr_dict({'alpha':alpha, 'sx':sx, 'sy':sy, 'mirror':mirror})
		ima.write_image(stack, im, EMUtil.ImageType.IMAGE_HDF, True)
	print_end_msg("params_3D_to_2D")

def find_struct(prj_stack, outdir, delta, ir, ou, Iter, rand_seed, trials, refine, MPI = False, debug = False):
	#from projection import cml_voronoi
	#cml_voronoi(prj_stack, outdir, delta, rand_seed, refine, debug)

	from projection import cml_find_struc, cml_export_struc, cml_open_proj, cml_head_log, cml_end_log
	from projection import cml_init_rnd
	from utilities  import print_begin_msg, print_msg, print_end_msg, start_time, running_time
	from copy       import deepcopy
	import time
	import os

	if os.path.exists(outdir): os.system( 'rm -rf ' + outdir )
	os.mkdir(outdir)

	if MPI:
		#ERROR('No MPI version yet!', 'find_struct', 1)
		from projection import cml_init_MPI, cml_end_mpilog
		from mpi        import mpi_barrier, mpi_reduce, mpi_bcast
		from mpi        import MPI_COMM_WORLD, MPI_FLOAT, MPI_INT, MPI_SUM
				
		main_node, myid, ncpu, loop = cml_init_MPI(trials)

		if myid == main_node:
			#  what is that ?? PAP
			try: os.remove('.tmp_txt_1a32b4')
			except: pass
			t_start = start_time()
			print_begin_msg('find_struct_MPI')
			cml_head_log(prj_stack, outdir, delta, ir, ou, rand_seed, ncpu, refine, loop * ncpu)

		Rnd = cml_init_rnd(loop * ncpu, rand_seed)
		Prj = cml_open_proj(prj_stack, ir, ou)

		best_val = 1.0e10
		best_Prj = None

		for i in xrange(loop):
			newPrj = deepcopy(Prj)
			newPrj, disc = cml_find_struc(newPrj, delta, outdir, i * ncpu + myid, Iter, Rnd[i * ncpu + myid], refine, False)

			if disc < best_val:
				best_val = disc
				best_Prj = deepcopy(newPrj)

			mpi_barrier(MPI_COMM_WORLD)

			for n in xrange(ncpu):
				if n == myid:
					f = open('.tmp_txt_1a32b4', 'a')
					f.write('\ntrials %s______________________________rnd: %d\n' % (str(i * ncpu + n).rjust(3, '0'), Rnd[i * ncpu + n]))
					f.write('          \__discrepancy: %10.3f\n' % disc)
					f.write('           \_%s\n' % time.ctime())
					f.close()
				
				mpi_barrier(MPI_COMM_WORLD)

		mpi_barrier(MPI_COMM_WORLD)

		# Which node found the best result?
		RES       = [0.0] * ncpu
		RES[myid] = best_val
		RES       = mpi_reduce(RES, ncpu, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
		RES       = RES.tolist()
		best_node = -1
		if myid  == main_node: best_node = RES.index(min(RES))
		best_node = mpi_bcast(best_node, 1, MPI_INT, main_node, MPI_COMM_WORLD)
		best_node = best_node.tolist()[0]

		if myid == best_node:
			cml_export_struc(prj_stack, outdir, best_Prj)
			cml_end_mpilog(best_Prj, best_val)

		if myid == main_node:
			# add txt logfile from the others node
			txt = open('.tmp_txt_1a32b4', 'r').readlines()
			for line in txt: print_msg(line)
			running_time(t_start)
			print_end_msg('find_struct_MPI')
			os.remove('.tmp_txt_1a32b4')
		
	else:
		t_start = start_time()
		print_begin_msg('find_struct')
		cml_head_log(prj_stack, outdir, delta, ir, ou, rand_seed, 1, refine, trials)   # 1 is ncpu

		Rnd = cml_init_rnd(trials, rand_seed)
		Prj = cml_open_proj(prj_stack, ir, ou)

		best_val = 1.0e10
		best_Prj = None
		for trial in xrange(trials):
			newPrj = deepcopy(Prj)
			print_msg('\ntrials %s______________________________rnd: %d\n' % (str(trial).rjust(3, '0'), Rnd[trial]))
			newPrj, disc = cml_find_struc(newPrj, delta, outdir, trial, Iter, Rnd[trial], refine, debug)
			print_msg('          \__discrepancy: %10.3f\n' % disc)
			print_msg('           \_%s' % time.ctime())

			if disc < best_val:
				best_val = disc
				best_Prj = deepcopy(newPrj)
	
		cml_export_struc(prj_stack, outdir, best_Prj)

		cml_end_log(best_Prj, best_val)
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

def header( stack, params, zero, one, randomize, fimport, fexport, fprint ):
        from string import split

        params = split( params )

        if len(fimport)>0: fimp = open( fimport, 'r' )
	
	if len(fexport)>0: fexp = open( fexport, 'w' )

	from utilities import write_header, file_type
	ext = file_type(stack)
	if(ext == "bdb"): DB = db_open_dict(stack)
	nimage = EMUtil.get_image_count( stack )
	for i in xrange(nimage):
		img = EMData()
		img.read_image( stack, i, True)

		if len(fimport)>0:
			line = fimp.readline()
			if len(line)==0 :
				print "Error: file " + fimport + " has only " + str(i) + " lines, while there are " + str(nimage) + " images in the file."
				return

			parmvalues = split( line )
			if len(params)!=len(parmvalues):
				print "Error: %d params need to be set, while %d values are provided in line %d of file." % ( len(params), len(parmvalues), i )
				return

			for j in xrange(len(params)):
				img.set_attr( params[j], extract_value(parmvalues[j]) )

			write_header( stack, img, i)
			#img.write_image( stack, i, EMUtil.ImageType.IMAGE_HDF, True)

		else:
			for p in params:
				if zero:
					img.set_attr( p, 0.0 )
				elif one:
					img.set_attr( p, 1.0 )
				elif randomize:
					from random import random, randint
					if p == "alpha" or p == "phi" or p == "psi": 
						v = random()*360.0
					elif p == "theta":
						v = random()*180.0
					elif p == "sx" or p == "sy" or p == "s2x" or p == "s2y":
						v = random()*4.0-2.0
					elif p == "mirror":
						v = randint(0,1)
					else:
						v = 0.0
					img.set_attr(p, v)						
				elif len(fexport)>0:
					fexp.write( "%15s   " % str(img.get_attr(p)) )
				elif fprint:
					print ("%15s   " % str(img.get_attr(p))),
				else:
					print "Error: no operation selected"
					return

			if zero or one or randomize:
				write_header( stack, img, i)
				#img.write_image( stack, i, EMUtil.ImageType.IMAGE_HDF, True)
			elif( len(fexport)>0 ):
				fexp.write( "\n" )
			elif fprint:
				print " "
			else:
				assert False
	if(ext == "bdb"): DB.close(stack)


def imgstat_ccc( stacks, rad ):
	from EMAN2 import EMUtil
	from utilities import get_im, model_circle
	from statistics import ccc

	if len(stacks)>3: ERROR("Error: ccc should be run on two stacks","imgstat",1)

	nimg1 = EMUtil.get_image_count( stacks[0] )
	nimg2 = EMUtil.get_image_count( stacks[1] )

	
	nimg = max(nimg1,nimg2)
	if(nimg2<nimg1): nimg2 = 1

        imgtmp = get_im( stacks[0] )

	if rad==-1:
		if len(stacks) == 3:  mask = get_im(stacks[2])
		else:                 mask = None
	else:
		if len(stacks) == 3:    ERROR("Error: Mask radius and mask file canot be given simultaneously","imgstat",1)
		else:
			nx = imgtmp.get_xsize()
			ny = imgtmp.get_ysize()
			nz = imgtmp.get_zsize()
			mask = model_circle( rad, nx, ny, nz )

	for i in xrange(nimg):
		img1 = get_im( stacks[0], i )

		if nimg2==1:
			if i==0:
				img2 = get_im( stacks[1] )
		else:
			img2 = get_im( stacks[1], i )

		val = ccc(img1, img2, mask)

		print "%10.5f" % val

def imgstat_fsc( stacks, fscfile, rad ):
	from utilities import get_im, model_circle
	from statistics import fsc_mask

	if len(stacks)>3: ERROR("Error: fsc should be run on two images","imgstat",1)

	img1 = get_im( stacks[0] )
	img2 = get_im( stacks[1] )

	nx = img1.get_xsize()
	ny = img1.get_ysize()
	nz = img1.get_zsize()

	if  EMUtil.get_image_count(stacks[0])>1: ERROR("Error: %s is an stack, fsc should be run on images","imgstat",1)

	if img2.get_xsize() != nx or img2.get_ysize() != ny or img2.get_zsize() != nz: ERROR("Error: input images has different sizes","imgstat",1)

	if rad==-1:
		if len(stacks) == 3: mask = get_im(stacks[2])
		else:                mask = None
	else:
		if len(stacks) == 3:  ERROR("Error: Mask radius and mask file canot be given simultaneously","imgstat",1)
		else:    mask = model_circle( rad, nx, ny, nz )

	fsc_mask( img1, img2, mask, filename=fscfile )
	
def imgstat_inf( stacks, rad ):
	from EMAN2 import EMUtil
	from utilities import get_im, model_circle
	if len(stacks)>2: ERROR("Error: inf should be run on one file","imgstat",1)

	nimg = EMUtil.get_image_count( stacks[0] )
	img1 = get_im( stacks[0] )
	
	nx = img1.get_xsize()
	ny = img1.get_ysize()
	nz = img1.get_zsize()

	if rad==-1:
		if len(stacks) == 2:  mask = get_im(stacks[1])
		else:                 mask = None
	else:
		if len(stacks) == 2:    ERROR("Error: Mask radius and mask file canot be given simultaneously","imgstat",1)
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


def normal_prj( prj_stack, outdir, refvol, r, niter, snr, sym, MPI=False ):
	def peak_range( nx, defocus, Cs, voltage, pixel ):
		from morphology import ctf_1d
		ctf = ctf_1d( nx, pixel, defocus, voltage, Cs )
    
		for i in xrange( 1, len(ctf)-1 ):
			prev = ctf[i-1]
			curt = ctf[i]
			next = ctf[i+1]

			if curt > prev and curt > next:
				freq = float(i)/nx
				return [freq-0.03, freq+0.02]

		assert false

	from utilities import getImage, get_im, model_circle, dropSpiderDoc, bcast_EMData_to_all,dropImage
	from projection import prep_vol, prgs
	from filter import filt_ctf, filt_btwo,filt_tophatb
	from statistics import ccc
	import os

	if MPI:
		from mpi import mpi_comm_size, mpi_comm_rank, mpi_barrier
		sys.argv = mpi_init(len(sys.argv),sys.argv)
		nproc = mpi_comm_size( MPI_COMM_WORLD )
		myid  = mpi_comm_rank( MPI_COMM_WORLD )
	else:
		nproc = 1
		myid  = 0

	img1st = getImage( prj_stack )
	nx = img1st.get_xsize()
	ny = img1st.get_ysize()

	if r < 0:
	    r = nx/2-1

	img_number     = EMUtil.get_image_count( prj_stack )
	img_node_start, img_node_end = MPI_start_end(img_number, nproc, myid )

	if myid==0:
    		os.system( "rm -rf " + outdir )
    		os.system( "mkdir " + outdir )

	if MPI: 
		mpi_barrier( MPI_COMM_WORLD )


	infofile = outdir + ("/progress%04d.txt" % (myid+1))
	info = open( infofile, 'w' )


	imgdata = []
	for i in xrange(img_node_start, img_node_end):
		img = get_im(prj_stack, i)
		imgdata.append(img)
	info.write( ' all imgs loaden\n' )
	info.flush( )


	odd_start = img_node_start%2
	eve_start = 1 - odd_start

	newimg = [None]*len(imgdata)

	pred = [None]* len(imgdata)
	for i in xrange( len(imgdata) ):
		pred[i] = [1.0, 1.0]

	from reconstruction import rec3D_MPI,rec3D_MPI_noCTF,rec3D_MPI_noCTF
        if refvol is None:
		fsc_file = outdir + "/fsc_init.dat"
		vol_file = outdir + "/vol_init.hdf"
		refvol, fscc = rec3D_MPI( imgdata, snr, sym, None, fsc_file, myid, 0, 1.0, odd_start, eve_start, None )
		bcast_EMData_to_all( refvol, myid )
		if myid==0: 
			refvol.write_image( vol_file )
		info.write( "inital reconstructed volume wrote to " + vol_file + "\n" )
		info.flush()


	for iter in xrange(niter) :
		volft, kb = prep_vol( refvol )
    		mask = model_circle( r, nx, ny )

		scales = []
		for i in xrange( len(imgdata) ) :
			exp_prj = imgdata[i]
			nx = imgdata[i].get_xsize()
			phi = exp_prj.get_attr( 'phi' )
			theta = exp_prj.get_attr( 'theta' )
			psi = exp_prj.get_attr( 'psi' )
			s2x = exp_prj.get_attr( 's2x' )
			s2y = exp_prj.get_attr( 's2y' )

			ref_prj = prgs( volft, kb, [phi, theta, psi, -s2x, -s2y] )
			ref_prj = filt_btwo( ref_prj, 0.01,0.1,0.2)
 
			defocus = exp_prj.get_attr( 'defocus' )
			wgh = exp_prj.get_attr( 'amp_contrast' )
			Cs = exp_prj.get_attr( 'Cs' )
			voltage = exp_prj.get_attr( 'voltage' )
			pixel = exp_prj.get_attr( 'Pixel_size' )


			nx = exp_prj.get_xsize()
			frange = peak_range( nx, defocus, Cs, voltage, pixel )


			ref_ctfprj = filt_ctf( ref_prj, defocus, Cs, voltage, pixel, wgh )
			ref_ctfprj2 = filt_ctf( ref_ctfprj, defocus, Cs, voltage, pixel, wgh )


			if exp_prj.get_attr('ctf_applied')==0.0:
				exp_ctfprj2 = filt_ctf( exp_prj,    defocus, Cs, voltage, pixel, wgh )
			else:
				exp_ctfprj2 = exp_prj.copy()

			ex_refprj = ref_ctfprj2
			ex_expprj = exp_ctfprj2

			tophat_ref = filt_tophatb( ex_refprj, frange[0], frange[1], False )
			tophat_exp = filt_tophatb( ex_expprj, frange[0], frange[1], False )
         
			ref_mask = tophat_ref*mask
			exp_mask = tophat_exp*mask
			curtccc = ccc( ref_mask, exp_mask, mask )
       
			try:
				a = exp_mask.dot( ref_mask ) / exp_mask.dot(exp_mask)
			except:
				print 'exception at myid, i:', myid, i
				a = 1.0

		        scales.append( a )
			info.write( "i, a, b, ccc, defocus:  %4d %10.5f %10.5f %10.5f %10.5f\n" %(i, a, 0.0, curtccc, defocus) )
			info.flush()

 	   	sum_scale = sum( scales )

		total_sum_scale = mpi_reduce( sum_scale, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD )
		total_sum_scale = mpi_bcast( total_sum_scale, 1, MPI_FLOAT, 0, MPI_COMM_WORLD)

		avg_scale = total_sum_scale[0][0]/img_number

		assert( len(imgdata)==len(scales) )


		for i in xrange( len(imgdata) ):
			s = scales[i] / avg_scale
			imgdata[i] *= s
			pred[i][0] *= s
			pred[i][1] /= s

    		scale_file = outdir + ('/newscale%04d_%04d.txt' % (myid, iter))
    		dropSpiderDoc( scale_file, pred )

		fsc_file = outdir + "/" + ( "fsc_%04d.dat" % iter )
		vol_file = outdir + "/" + ( "vol_%04d.spi" % iter )
		info.write( 'running reconstruction\n' )
		info.flush()
		refvol,fscc = rec3D_MPI( imgdata, snr, sym, None, fsc_file, myid, 0, 1.0, odd_start, eve_start, None )
		info.write( 'reconstruction finished\n' )
		info.flush()



		if myid==0:
			dropImage( refvol, vol_file )
			info.write( "reconstructed volume wrote to " + vol_file  + "\n")
			info.flush()

		bcast_EMData_to_all( refvol, myid )
		[mean,sigma,fmin,fmax] = Util.infomask( refvol, None, True )
		info.write( 'vol all after reconstruction, myid: %d %10.3e %10.3e %10.3e %10.3e\n' % ( myid, mean, sigma, fmin, fmax ) )
		info.flush()

	prj_file = outdir + ("/prj%04d.hdf" % myid)
	info.write( "writing resulting projections to file " + prj_file + '\n' )
	info.flush()

	for iprj in xrange( len(imgdata) ):
		imgdata[iprj].write_image( prj_file, iprj )


	info.write( "output written to file " + prj_file + '\n' )
	info.flush()

def incvar(prefix, nfile, nprj, output, fl, fh, radccc, writelp, writestack):
	from statistics import variancer, ccc
	from string import atoi, replace, split, atof
	from EMAN2 import EMUtil
	from utilities import get_im, circumference, model_circle, dropImage
	from filter import filt_btwl
	from math import sqrt
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


class file_set :

	def __init__( self, prefix, nfile ):
		self.files = [None] * nfile
		self.fends = [None] * nfile

		totimg = 0
		for i in xrange(nfile):
			self.files[i] = prefix + ("%04d.hdf" % i )
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



def incvar_mpi(prefix, nfile, nprj, output, fl, fh, radccc, writelp, writestack):
	from statistics import variancer, ccc
	from string import atoi, replace, split, atof
	from EMAN2 import EMUtil
	from utilities import get_im, circumference, model_circle, dropImage, info
	from filter import filt_btwl, filt_gaussl, filt_tanl
	from math import sqrt
	import os
	from mpi import mpi_comm_rank, mpi_comm_size, MPI_COMM_WORLD

        myid = mpi_comm_rank( MPI_COMM_WORLD )
        ncpu = mpi_comm_size( MPI_COMM_WORLD )
        finf = open( "progress%04d.txt" % myid, "w" )
	
	all_varer = variancer()
	odd_varer = variancer()
	eve_varer = variancer()
 
	filname = prefix + "0000.hdf"
	n = get_im(filname, 0).get_xsize()
	
	if os.path.exists(output):		os.system("rm -f "+output)
	if os.path.exists('stack_'+output):	os.system("rm -f stack_"+output)
	if os.path.exists('odd_stack_'+output):	os.system("rm -f odd_stack_"+output)
	if os.path.exists('eve_stack_'+output):	os.system("rm -f eve_stack_"+output)
	if os.path.exists('avg_'+output):	os.system("rm -f avg_"+output)

	cccmask = model_circle(radccc, n, n, n)
	scale = sqrt( nprj )
	radcir = n/2

        mystack = file_set( prefix, nfile )
        nimage = mystack.nimg()
        mystart, myend = MPI_start_end( nimage, ncpu, myid ) 
        ndump = 4

        print 'myid, nimage, start, end: ', myid, nimage, mystart, myend

	lpstack = "btwl_cir_prj%04d.hdf" % myid
	iwrite = 0
	iprint = 0
	for i in xrange(mystart, myend):
		filename, imgid = mystack.get( i ) 	
		finf.write( "processing %4d (%s:%d)\n" % (i, filename, imgid) )
		finf.flush()

		img = get_im( filename, imgid )
		img *= scale
		img = circumference( img, radcir, radcir+1 )
		img = filt_tanl(img, fl, fh)

		if writelp:
			img.write_image(lpstack, iwrite)
			iwrite += 1

		if i%2==0: 
			odd_varer.insert(img)
		else: 
			eve_varer.insert(img)

		all_varer.insert(img)

		'''
		if (i-mystart)%ndump==(ndump-1):
			odd_var = odd_varer.mpi_getvar(myid, 0)
			eve_var = eve_varer.mpi_getvar(myid, 0)
			all_var = all_varer.mpi_getvar(myid, 0)
			

			if myid==0 :
				iprint += 1
				print 'ntot, ccc: %6d %10.3f' % (ncpu*ndump*iprint, ccc(odd_var, eve_var, cccmask))  
				if writestack:
					odd_var.write_image( 'odd_stack_' + output, iwrite )
					eve_var.write_image( 'eve_stack_' + output, iwrite )
					all_var.write_image( 'stack_' + output, iwrite )
					iwrite += 1
		'''

	all_var = all_varer.mpi_getvar(myid, 0)
	odd_var = odd_varer.mpi_getvar(myid, 0)
	eve_var = eve_varer.mpi_getvar(myid, 0)
	avg = all_varer.mpi_getavg(myid, 0)

	if myid==0:
		print 'ntot, ccc: %6d %10.3f' % (nimage, ccc(odd_var, eve_var, cccmask))  


	if myid==0 and writestack:
		all_var.write_image( 'stack_' + output, iwrite )
		odd_var.write_image( 'odd_stack_' + output, iwrite )
		eve_var.write_image( 'eve_stack_' + output, iwrite )

	if myid==0:
		avg.write_image( 'avg_' + output, 0 )
		#all_var = circumference( all_var, radcir, radcir+1 )
		all_var.write_image( output, 0 )






def factcoords2D( prj_stack, avgvol_stack = None, eigvol_stack = None, output = None, rad = -1, neigvol = -1, of = "txt"):
	from utilities import get_im, model_circle, model_blank

        if of=="txt":
		foutput = open( output, "w" );


	nx = get_im( prj_stack ).get_xsize()
	ny = nx

	eigvols = []
	if(neigvol < 0): neigvol = EMUtil.get_image_count( eigvol_stack )
	for i in xrange(neigvol):
		eigvols.append( get_im(eigvol_stack, i) )

	if( avgvol_stack != None):
		avgvol = get_im( avgvol_stack )

	m = model_circle( rad, nx, ny )
	nimage = EMUtil.get_image_count( prj_stack )

	for i in xrange( nimage ) :
		exp_prj = get_im( prj_stack, i )
		if(avgvol_stack != None):  Util.sub_img(exp_prj, avgvol)
		if of=="hdf": img = model_blank( len(eigvols) )

		for j in xrange( neigvol ) :

			d = exp_prj.cmp( "dot", eigvols[j], {"negative":0, "mask":m} )

			if of=="hdf":
				img.set_value_at( j, 0, 0, d )
			else:
				foutput.write( "    %e  " % d )
				foutput.flush()
                if of=="hdf":
			img.write_image( output, i )
		else:
			foutput.write( "\n" )

def factcoords3D( prj_stack, avgvol_stack, eigvol_stack, output, rad, neigvol, of):
	from utilities import get_im, getImage, model_circle, model_blank
	from projection import prgs, prep_vol
	from filter import filt_ctf
	from statistics import im_diff
	from utilities import memory_usage

        if of=="txt":
		foutput = open( output, "w" );


	nx = get_im( prj_stack ).get_xsize()
	ny = nx

	eigvols = []
	if(neigvol < 0): neigvol = EMUtil.get_image_count( eigvol_stack )
	for i in xrange(neigvol):
		eigvols.append( get_im(eigvol_stack, i) )

	avgvol = getImage( avgvol_stack )
	volft, kb = prep_vol( avgvol )
        print 'prepare vol done'
	eigvolfts = []
	for e in eigvols:
    		eigvolfts.append( prep_vol(e) )
        print 'prepare eigvol done '

	m = model_circle( rad, nx, ny )
	nimage = EMUtil.get_image_count( prj_stack )

	for i in xrange( nimage ) :
		exp_prj = get_im( prj_stack, i )
                
		phi = exp_prj.get_attr( 'phi' )
		theta = exp_prj.get_attr( 'theta' )
		psi = exp_prj.get_attr( 'psi' )
		s2x = exp_prj.get_attr( 's2x' )
		s2y = exp_prj.get_attr( 's2y' )
		defocus = exp_prj.get_attr( 'defocus' )
		wgh = exp_prj.get_attr( 'amp_contrast' )
		Cs = exp_prj.get_attr( 'Cs' )
		voltage = exp_prj.get_attr( 'voltage' )
		pixel = exp_prj.get_attr( 'Pixel_size' )

		assert exp_prj.get_attr( "ctf_applied" ) == 0.0

		shift_params = {"filter_type" : Processor.fourier_filter_types.SHIFT,
				"x_shift" : s2x, "y_shift" : s2y, "z_shift" : 0.0}
		exp_prj =  Processor.EMFourierFilter(exp_prj, shift_params)
		
		ref_prj = prgs( volft, kb, [phi, theta, psi, 0.0, 0.0] )
		ref_ctfprj = filt_ctf( ref_prj, defocus, Cs, voltage, pixel, wgh )
		
		diff,a,b = im_diff( ref_ctfprj, exp_prj, m)
		
		if of=="hdf":
			img = model_blank( len(eigvols) )


		for j in xrange( neigvol ) :
			eft = eigvolfts[j]

			ref_eigprj = prgs( eft[0], eft[1], [phi, theta, psi, 0.0, 0.0] )
			ref_ctfeigprj = filt_ctf( ref_eigprj, defocus, Cs, voltage, pixel, wgh )

			d = diff.cmp( "dot", ref_ctfeigprj, {"negative":0, "mask":m} )

			if of=="hdf":
				img.set_value_at( j, 0, 0, d )
			else:
				foutput.write( "    %e  " % d )
				foutput.flush()

                if of=="hdf":
			img.write_image( output, i )
		else:
			foutput.write( "\n" )

		print i, ' done'

def refvol( vollist, fsclist, output, mask ):
	from utilities     import getImage, read_fsc
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
		vols[i] = getImage( vollist[i] )
		print 'rawvol, resolution: ', vollist[i], fsclist[i]

	m    = getImage( mask )
	volfs = filt_vols( vols, fscs, m )

	for i in xrange(nvol):
		volfs[i].write_image( output, i )

# -- K-means main ---------------------------------------------------------------------------

# K-means main driver
def k_means_main(stack, out_dir, maskname, opt_method, K, rand_seed, maxit, trials, critname, CTF = False, F = 0, T0 = 0, MPI = False, SA2 = False, DEBUG = False):
	from utilities 	import print_begin_msg, print_end_msg, print_msg
	from statistics import k_means_criterion, k_means_export, k_means_open_im, k_means_headlog
	import sys

	#== NEW VERSION ==

	if stack.split(':')[0] == 'bdb': BDB = True
	else:                            BDB = False
	
	# check entry
	if maskname != None:
		if (maskname.split(':')[0] != 'bdb' and BDB) or (maskname.split(':')[0] == 'bdb' and not BDB):
		        ERROR('Both mask and stack must be on bdb format!', 'k-means groups', 1)
		
	N = EMUtil.get_image_count(stack)

	if MPI:
		from statistics import k_means_cla_MPI, k_means_SSE_MPI, k_means_init_MPI
		
		main_node, myid, ncpu, N_start, N_stop = k_means_init_MPI(N)

		if myid == main_node:
			print_begin_msg('k-means')
			k_means_headlog(stack, out_dir, opt_method, N, K, critname, maskname, trials, maxit, CTF, T0, F, SA2, rand_seed, ncpu)

		[im_M, mask, ctf, ctf2] = k_means_open_im(stack, maskname, N_start, N_stop, N, CTF, BDB)

		if   opt_method == 'cla':
			[Cls, assign] = k_means_cla_MPI(im_M, mask, K, rand_seed, maxit, trials, [CTF, ctf, ctf2], myid, main_node, N_start, N_stop, F, T0, SA2)
		elif opt_method == 'SSE':
			[Cls, assign] = k_means_SSE_MPI(im_M, mask, K, rand_seed, maxit, trials, [CTF, ctf, ctf2], myid, main_node, ncpu, N_start, N_stop, F, T0, SA2)
		else:
			ERROR('opt_method %s unknown!' % opt_method, 'k-means', 1)
			sys.exit()

		if myid == main_node:
			crit = k_means_criterion(Cls, critname)
			k_means_export(stack, Cls, crit, assign, out_dir, BDB)
			print_end_msg('k-means')

	else:
		from statistics import k_means_classical, k_means_SSE

		print_begin_msg('k-means')
		k_means_headlog(stack, out_dir, opt_method, N, K, critname, maskname, trials, maxit, CTF, T0, F, SA2, rand_seed, 1)
		[im_M, mask, ctf, ctf2] = k_means_open_im(stack, maskname, 0, N, N, CTF, BDB)
		
		if   opt_method == 'cla':
			[Cls, assign] = k_means_classical(im_M, mask, K, rand_seed, maxit, trials, [CTF, ctf, ctf2], F, T0, SA2, DEBUG)
		elif opt_method == 'SSE':
			[Cls, assign] = k_means_SSE(im_M, mask, K, rand_seed, maxit, trials, [CTF, ctf, ctf2], F, T0, SA2, DEBUG)
		else:
			ERROR('opt_method %s unknown!' % opt_method, 'k-means', 1)
			sys.exit()
			
		crit = k_means_criterion(Cls, critname)
		k_means_export(stack, Cls, crit, assign, out_dir, BDB)
		print_end_msg('k-means')
			
# K-means groups driver
def k_means_groups(stack, out_file, maskname, opt_method, K1, K2, rand_seed, maxit, trials, crit, CTF, F, T0, SA2, MPI=False, DEBUG=False):

	# check entry
	if stack.split(':')[0] == 'bdb': BDB = True
	else:                            BDB = False
	
	if maskname != None:
		if (maskname.split(':')[0] != 'bdb' and BDB) or (maskname.split(':')[0] == 'bdb' and not BDB):
		        ERROR('Both mask and stack must be on bdb format!', 'k-means groups', 1)
	
	if MPI:
		from statistics import k_means_groups_MPI
		[rescrit, KK] = k_means_groups_MPI(stack, out_file, maskname, opt_method, K1, K2, rand_seed, maxit, trials, crit, CTF, F, T0, SA2)
		
	else:
		from statistics import k_means_groups_serial
		[rescrit, KK] = k_means_groups_serial(stack, out_file, maskname, opt_method, K1, K2, rand_seed, maxit, trials, crit, CTF, F, T0, SA2, DEBUG)
		
		return rescrit, KK
	
# ----------------------------------------------------------------------------------------------
