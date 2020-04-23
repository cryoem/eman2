

































































































































































































































































"""0
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
		"""















































































































































































































































































































































































































































































































































































































































































































































































































































































































































































'''1
def ORGali2d_c(stack, outdir, maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", center=-1, maxit=0, \
		CTF=False, snr=1.0, Fourvar = False, user_func_name="ref_ali2d", CUDA=False, GPU=0, MPI=False):
	if MPI:
		ali2d_c_MPI(stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, center, maxit, CTF, snr, Fourvar, user_func_name, CUDA, GPU)
		return

	from sp_utilities    import model_circle, drop_image, get_image, get_input_from_string, get_params2D
	from sp_statistics   import fsc_mask, sum_oe, hist_list
	from sp_alignment    import Numrinit, ringwe, ali2d_single_iter
	from sp_pixel_error  import pixel_error_2D
	from sp_filter       import filt_ctf, filt_table, filt_tophatb
	from sp_fundamentals import fshift
	from sp_utilities    import print_begin_msg, print_end_msg, print_msg
	from sp_fundamentals import fft, rot_avg_table
	from sp_utilities    import write_text_file, file_type
	import os
		
	print_begin_msg("ali2d_c")

	if os.path.exists(outdir):   ERROR('Output directory exists, please change the name and restart the program', " ORGali2d_c", 1)
	os.mkdir(outdir)

	import sp_user_functions
	user_func = sp_user_functions.factory[user_func_name]

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
		if ima.get_attr_default('ctf_applied', 0) > 0:	ERROR("data cannot be ctf-applied", "ORGali2d_c", 1)
		from sp_morphology   import ctf_img
		ctf_2_sum = EMData(nx, nx, 1, False)
	else:
		ctf_2_sum = None
	if  Fourvar:
		from sp_statistics   import add_ave_varf

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
		        		this_error = pixel_error_2D(old_ali_params[im][0], old_ali_params[im][1], old_ali_params[im][2], alphan, sxn, syn, last_ring)
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
	from sp_utilities import write_headers
	write_headers(stack, data, list_of_particles)
	print_end_msg("ali2d_c")

def ORGali2d_c_MPI(stack, outdir, maskfile=None, ir=1, ou=-1, rs=1, xr="4 2 1 1", yr="-1", ts="2 1 0.5 0.25", center=-1, maxit=0, CTF=False, snr=1.0, \
			Fourvar = False, user_func_name="ref_ali2d", CUDA=False, GPU=0):

	from sp_utilities    import model_circle, model_blank, drop_image, get_image, get_input_from_string
	from sp_utilities    import reduce_EMData_to_root, bcast_EMData_to_all, send_attr_dict, file_type, bcast_number_to_all, bcast_list_to_all
	from sp_statistics   import fsc_mask, sum_oe, add_ave_varf_MPI, hist_list
	from sp_alignment    import Numrinit, ringwe, ali2d_single_iter
	from sp_pixel_error  import pixel_error_2D
	from sp_filter       import filt_table, filt_ctf, filt_tophatb
	from numpy        import reshape, shape
	from sp_fundamentals import fshift, fft, rot_avg_table
	from sp_utilities    import write_text_file, get_params2D, set_params2D
	from sp_utilities    import print_msg, print_begin_msg, print_end_msg
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
		import sp_global_def
		sp_global_def.LOGFILE =  os.path.join(outdir, sp_global_def.LOGFILE)

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
		import sp_user_functions
		user_func = sp_user_functions.factory[user_func_name]
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
	list_of_particles = bcast_list_to_all(list_of_particles, myid,  source_node = main_node)
	
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
		if ima.get_attr_default('ctf_applied', 0) > 0:	ERROR("data cannot be ctf-applied", "ORGali2d_c_MPI", 1,myid)
		from sp_filter import filt_ctf
		from sp_morphology   import ctf_img
		ctf_2_sum = EMData(nx, nx, 1, False)
	else:
		ctf_2_sum = None
	if  Fourvar:
		from sp_statistics   import add_ave_varf

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
					sx_sum, sy_sum, nope = ali2d_single_iter(data, numr, wr, cs, tavg, cnx, cny, xrng[N_step], yrng[N_step], step[N_step], mode, CTF)
					
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
		        			this_error = pixel_error_2D(old_ali_params[im*4:im*4+3], [alpha, sx, sy], last_ring)
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
		from sp_utilities import file_type
		if(file_type(stack) == "bdb"):
			from sp_utilities import recv_attr_dict_bdb
			recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
		else:
			from sp_utilities import recv_attr_dict
			recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	else:           send_attr_dict(main_node, data, par_str, image_start, image_end)
	if myid == main_node: print_end_msg("ali2d_c_MPI")
'''






































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































'''2
def ali3d_abandoned(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", apsi = "-1", deltapsi = "-1", startpsi = "-1",
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
	        	delta, an, apsi, deltapsi, startpsi, center, maxit, CTF, snr, ref_a, sym, user_func_name,
			fourvar, npad, debug, termprec)
		return

	from sp_alignment      import proj_ali_incore, proj_ali_incore_local
	from sp_utilities      import model_circle, drop_image, get_image, get_input_from_string
	from sp_utilities      import get_params_proj
	from sp_utilities      import estimate_3D_center, rotate_3D_shift
	from sp_filter         import filt_params, fit_tanh, filt_tanl, filt_ctf
	from sp_statistics     import fsc_mask
	from sp_utilities      import print_begin_msg, print_end_msg, print_msg
	from sp_alignment      import Numrinit, prepare_refrings
	from sp_projection     import prep_vol

	import sp_user_functions
	import os
	import types
	from math			import radians, sin, cos

	user_func = sp_user_functions.factory[user_func_name]

	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ali3d", 1)
	os.mkdir(outdir)
	import sp_global_def
	sp_global_def.LOGFILE =  os.path.join(outdir, sp_global_def.LOGFILE)
	print_begin_msg("ali3d")

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
	print_msg("User function               : %s\n"%(user_func_name))

	if maskfile :
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else                                  : mask3D = maskfile
	else          :   mask3D = model_circle(last_ring, nx, nx, nx)
	mask2D = model_circle(last_ring, nx, nx) - model_circle(first_ring, nx, nx)
	numr   = Numrinit(first_ring, last_ring, rstep, "F")

	if CTF:
		from sp_reconstruction import recons3d_4nn_ctf
		from sp_filter         import filt_ctf
	else: from sp_reconstruction import recons3d_4nn

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
			ref_angles = even_angles(delta[N_step], symmetry=sym, method = ref_a, phiEqpsi = "Minus")

			pixer  = [0.0]*nima
			neworient = [[0.0, 0.0, 0.0, 0.0, 0.0, -2.0e23] for i in xrange(nima)]
			Torg = []
			pikorg = [0.0]*nima
			for im in xrange( nima ):
				Torg.append(data[im].get_attr('xform.projection'))
				pikorg[im] = data[im].get_attr_default("previousmax",-1.0e23)

			for refang in ref_angles:
				n1 = sin(radians(refang[1]))*cos(radians(refang[0]))
				n2 = sin(radians(refang[1]))*sin(radians(refang[0]))
				n3 = cos(radians(refang[1]))

				for im in xrange(nima):
					volft, kb = prep_vol(vol)
					refrings = [None]
				
					if an[N_step] == -1:
						if(refrings[0] == None): refrings = refprojs( volft, kb, [refang], numr, mode, wr )
						peak, pixel_error = proj_ali_incore(data[im], refrings, numr, xrng[N_step], yrng[N_step], step[N_step], sym=sym)
						if(peak > neworient[im][-1]):
							# I stopped here realizing it would include conversion of data[im] to polar coords at the bottom of the loop,
							#    thus significantly slowing down the code.
							pass
					else:
						if(comparedirections):
							if(refrings[0] == None):
								refrings = refprojs( volft, kb, [refang], numr, mode, wr )
						peak, pixel_error = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step],sym=sym)
					data[im].set_attr("previousmax", peak)
			if center == -1 and sym[0] == 'c':
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center(data)
				msg = "Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
				print_msg(msg)
				if int(sym[1]) > 1:
					cs[0] = cs[1] = 0.0
					print_msg("For symmetry group cn (n>1), we only center the volume in z-direction\n")
				rotate_3D_shift(data, [-cs[0], -cs[1], -cs[2]])

			del volft, kb

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
			from sp_utilities import write_headers
			#from utilities import write_select_headers
			if CTF:
				for dat in data:  dat.set_attr('ctf_applied',0)
			write_headers(stack, data, list_of_particles)
			#list_params= ['ID','xform.projection']
			#write_select_headers(stack, data, list_params)
			if CTF:
				for dat in data:  dat.set_attr('ctf_applied', 1)
	print_end_msg("ali3d")
'''

'''3
def Xali3d_MPI_chunks(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", apsi = "-1", deltapsi = "-1", startpsi = "-1",
	    center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",  user_func_name = "ref_ali3d",
	    fourvar = True, npad = 4, debug = False, termprec = 0.0):

	from sp_alignment       import Numrinit, prepare_refrings, proj_ali_incore, proj_ali_incore_local, proj_ali_incore_local_psi
	from sp_utilities       import model_circle, get_image, drop_image, get_input_from_string
	from sp_utilities       import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all
	from sp_utilities       import send_attr_dict
	from sp_utilities       import get_params_proj, file_type
	from sp_fundamentals    import rot_avg_image
	import os
	import types
	from sp_utilities       import print_begin_msg, print_end_msg, print_msg
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi             import mpi_reduce, MPI_INT, MPI_SUM
	from sp_filter          import filt_ctf
	from sp_projection      import prep_vol, prgs
	from sp_statistics      import hist_list, varf3d_MPI
	from sp_applications    import MPI_start_end


	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ali3d_MPI", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		os.mkdir(outdir)
		import sp_global_def
		sp_global_def.LOGFILE =  os.path.join(outdir, sp_global_def.LOGFILE)
		print_begin_msg("ali3d_MPI")
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

	if apsi == "-1":
		apsi = [-1] * lstp
	else:
		apsi = get_input_from_string(apsi)

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
		import sp_user_functions
		user_func = sp_user_functions.factory[user_func_name]

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
		print_msg("Angular search range (phi and theta)       : %s\n"%(an))
		print_msg("Angular search range (psi)                 : %s\n"%(apsi))
		print_msg("Delta psi                   : %s\n"%(deltapsi))
		print_msg("Start psi                   : %s\n"%(startpsi))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Percentage of change for termination: %f\n"%(termprec))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Reference projection method : %s\n"%(ref_a))
		print_msg("Symmetry group              : %s\n\n"%(sym))
		print_msg("User function               : %s\n"%(user_func_name))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else: mask3D = model_circle(last_ring, nx, nx, nx)

	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	mask2D  = model_circle(last_ring,nx,nx) - model_circle(first_ring,nx,nx)

	fscmask = mask3D  #model_circle(last_ring,nx,nx,nx)  For a fancy mask circle would work better  PAP 7/21/11
	if CTF:
		from sp_reconstruction import rec3D_MPI
		from sp_filter         import filt_ctf
	else:	 from sp_reconstruction import rec3D_MPI_noCTF

	if myid == main_node:
       		if file_type(stack) == "bdb":
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
	list_of_particles = bcast_list_to_all(list_of_particles, myid,  source_node = main_node)

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
		if im == main_node :  disps.append(0)
		else:                 disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append(ie-ib)

	pixer = [0.0]*nima
	cs = [0.0]*3
	total_iter = 0
	# do the projection matching
	for N_step in xrange(lstp):
		terminate = 0
		Iter = -1
 		while Iter < max_iter-1 and terminate == 0:
			Iter += 1
			total_iter += 1
			if myid == main_node:
				start_time = time()
				print_msg("\nITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f, delta psi = %5.2f, start psi = %5.2f\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step],deltapsi[N_step],startpsi[N_step]))

			volft, kb = prep_vol(vol)
			refrings = prepare_refrings_chunks(volft, kb, nx, delta[N_step], ref_a, sym, numr, True, ant = max(an[N_step],0.0)*1.1)  # 1.1 is to have extra safety
			del volft, kb
			if myid == main_node:
				print_msg("Time to prepare rings: %d\n" % (time()-start_time))
				start_time = time()

			for im in xrange(nima):
				if deltapsi[N_step] > 0.0:
					from sp_alignment import proj_ali_incore_delta
					peak, pixer[im] = proj_ali_incore_delta(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],startpsi[N_step],deltapsi[N_step],finfo)						
				elif an[N_step] == -1:
					peak, pixer[im] = proj_ali_incore_chunks(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],finfo)
				else:
					if apsi[N_step] == -1:
						peak, pixer[im] = proj_ali_incore_local_chunks(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step],finfo, sym = sym)
					else:
						peak, pixer[im] = proj_ali_incore_local_psi(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step],apsi[N_step],finfo)
				data[im].set_attr("previousmax", peak)

			if myid == main_node:
				print_msg("Time of alignment = %d\n"%(time()-start_time))
				start_time = time()

			#output pixel errors
			from mpi import mpi_gatherv
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			terminate = 0
			if myid == main_node:
				recvbuf = map(float, recvbuf)
				from sp_statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				if region[0] < 0.0:  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				# Terminate if 95% within 1 pixel error
				im = 0
				for lhx in xrange(lhist):
					if region[lhx] > 1.0: break
					im += histo[lhx]
				precn = 100*float(total_nima-im)/float(total_nima)
				msg = " Number of particles that changed orientations %7d, percentage of total: %5.1f\n"%(total_nima-im, precn)
				print_msg(msg)
				if precn <= termprec:  terminate = 1
				del region, histo
			del recvbuf
			terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
			terminate = int(terminate[0])

			if center == -1 and sym[0] == 'c':
				from sp_utilities      import estimate_3D_center_MPI, rotate_3D_shift
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					print_msg(msg)
				if int(sym[1]) > 1:
					cs[0] = cs[1] = 0.0
					if myid == main_node:
						print_msg("For symmetry group cn (n>1), we only center the volume in z-direction\n")
				cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(data, cs)

			# write out headers, under MPI writing has to be done sequentially
			mpi_barrier(MPI_COMM_WORLD)
			par_str = ['xform.projection', 'previousmax', 'ID']
			if myid == main_node:
	   			if(file_type(stack) == "bdb"):
	        			from sp_utilities import recv_attr_dict_bdb
	        			recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        		else:
	        			from sp_utilities import recv_attr_dict
	        			recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
				print_msg("Time to write header information= %d\n"%(time()-start_time))
				start_time = time()
	        	else:	       send_attr_dict(main_node, data, par_str, image_start, image_end)

			if CTF: vol, fscc = rec3D_MPI(data, snr, sym, fscmask, os.path.join(outdir, "resolution%04d"%(total_iter)), myid, main_node, npad = npad)
			else:   vol, fscc = rec3D_MPI_noCTF(data, sym, fscmask, os.path.join(outdir, "resolution%04d"%(total_iter)), myid, main_node, npad = npad)

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


	if myid == main_node: print_end_msg("ali3d_MPI")
'''
























































































































































































































































































































































































































































































































































































































































"""4
		# finfo = None
		import os
		outdir = "./"
		# info_file = os.path.join(outdir, "progress%04d"%myid)
		
		# finfo = open(info_file, 'w')
		if "filename_first_time" not in sali3d_base.__dict__:
			import datetime
			sali3d_base.filename_first_time = os.path.join(outdir, "progress%04d_%s"%(myid, datetime.datetime.now().strftime('%Y-%m-%d--%H-%M-%S-%f')[:-3]))
			finfo = open(sali3d_base.filename_first_time, "w")
		else: 
			finfo = open(sali3d_base.filename_first_time, "a")
			finfo.write("=======================================================================================\n")
		# finfo = open(info_file, 'a')
		"""

























"""5
	user_func = Tracker["constants"]["user_func"]
	if ref_vol:
		#vol = do_volume_mrk01(ref_vol, Tracker, 0, mpi_comm)
		ref_data = [ref_vol, Tracker, 0, mpi_comm]
		vol = user_func(ref_data)
	else:
		#vol = do_volume_mrk01(data, Tracker, 0, mpi_comm)
		ref_data = [data, Tracker, 0, mpi_comm]
		vol = user_func(ref_data)
	"""
































































































































































"""  Have to think about it PAP6
				if(  changes > saturatecrit ):
					if( Iter == 1 ):
						log.add("Will continue even though %4.2f images did not find better orientations"%saturatecrit)
					else:
						terminate = 1
						log.add("...............")
						log.add(">>>>>>>>>>>>>>>   Will terminate as %4.2f images did not find better orientations"%saturatecrit)
				if( terminate == 0 ):
					historyofchanges.append(changes)
					historyofchanges = historyofchanges[:3]
					historyofchanges.sort()
					if( (historyofchanges[-1]-historyofchanges[0])/2/(historyofchanges[-1]+historyofchanges[0]) <0.05 ):
						terminate = 1
						log.add("...............")
						log.add(">>>>>>>>>>>>>>>   Will terminate as orientations do not improve anymore")
				"""














"""7
					if( lhx > saturatecrit):
						if( Iter == 1 ):
							log.add("Will continue even though %4.2f images had pixel error < %5.2f"%(saturatecrit,pixercutoff))
						else:
							terminate = 1
							log.add("...............")
							log.add(">>>>>>>>>>>>>>>   Will terminate as %4.2f images had pixel error < %5.2f"%(saturatecrit,pixercutoff))
					"""

















































"""8
			#=========================================================================
			if(False):  #total_iter%1 == 5 or terminate):
				# gather parameters
				params = []
				previousmax = []
				for im in data:
					t = get_params_proj(im)
					params.append( [t[0], t[1], t[2], t[3]/shrinkage, t[4]/shrinkage] )
					#if(t[3] >0.0 or t[4]>0.0):  print  "  ERRROR  ",t
					previousmax.append(im.get_attr("previousmax"))
				assert(nima == len(params))
				params = wrap_mpi_gatherv(params, 0, mpi_comm)
				if myid == 0:
					assert(total_nima == len(params))
				previousmax = wrap_mpi_gatherv(previousmax, 0, mpi_comm)
				if myid == main_node:
					from sp_utilities import write_text_row, write_text_file
					write_text_row(params, "soft/params%04d.txt"%total_iter)
					write_text_file(previousmax, "soft/previousmax%04d.txt"%total_iter)


				del previousmax, params
				i = 1
				while data[0].has_attr("xform.projection" + str(i)):
					params = []
					previousmax = []
					for im in data:

						try:
							#print  im.get_attr("xform.projection" + str(i))
							t = get_params_proj(im,"xform.projection" + str(i))
						except:
							print " NO XFORM  ",myid, i,im.get_attr('ID')
							from sys import exit
							exit()

						params.append( [t[0], t[1], t[2], t[3]/shrinkage, t[4]/shrinkage] )
					assert(nima == len(params))
					params = wrap_mpi_gatherv(params, 0, mpi_comm)
					if myid == 0:
						assert(total_nima == len(params))
					if myid == main_node:
						write_text_row(params, "soft/params-%04d-%04d.txt"%(i,total_iter))
					del previousmax, params
					i+=1


			if( ( terminate or (Iter == max_iter) ) and (myid == main_node) ):
				if( type(stack) is types.StringType ):
					from EMAN2 import Vec2f, Transform
					from EMAN2db import db_open_dict
					DB = db_open_dict(stack)
					for im in xrange(len(params)):
						t = Transform({"type":"spider","phi":params[im][0],"theta":params[im][1],"psi":params[im][2]})
						t.set_trans(Vec2f(-params[im][3], -params[im][4]))
						DB.set_attr(particle_ids[im], "xform.projection", t)
					DB.close()
				else:
					for im in xrange(len(params)): set_params_proj(stack[particle_ids[im]], params[im])
			"""







'''9
def sali3d_base_horatio_01(stack, ref_vol = None, Tracker = None, rangle = 0.0, rshift = 0.0, mpi_comm = None, log = None):	
	"""
		parameters: list of (all) projections | reference volume is optional, the data is shrank, 
		  the program does not know anything about shrinking| ...
		Data is assumed to be CTF multiplied and the ctf_applied flag to be set.
		The alignment done depends on nsoft:
					 nsoft = 0 & an = -1: exhaustive deterministic
					 nsoft = 0 & an > 0 : local deterministic
					 nsoft = 1 shc
					 nsoft >1  shc_multi
		
	"""

	from sp_alignment       import Numrinit, prepare_refrings, generate_indices_and_refrings
	from sp_alignment       import proj_ali_incore,  proj_ali_incore_zoom,  proj_ali_incore_local, proj_ali_incore_local_zoom
	from sp_alignment       import shc, center_projections_3D, ringwe
	from sp_utilities       import bcast_number_to_all, bcast_EMData_to_all, 	wrap_mpi_gatherv, wrap_mpi_bcast, model_blank, print_from_process
	from sp_utilities       import get_im, file_type, model_circle, get_input_from_string, get_params_proj, set_params_proj, pad
	from sp_utilities       import even_angles
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier, mpi_reduce, MPI_INT, MPI_SUM
	from sp_projection      import prep_vol
	from sp_statistics      import hist_list
	from sp_applications    import MPI_start_end
	from sp_filter          import filt_ctf, filt_table
	from sp_global_def      import Util
	from sp_fundamentals    import resample, fshift
	from sp_multi_shc       import shc_multi
	#from development     import do_volume_mrk01
	import sp_user_functions
	from EMAN2           import EMUtil, EMData
	import types
	from time            import time

	nsoft            = Tracker["nsoft"]
	saturatecrit     = Tracker["saturatecrit"]
	pixercutoff      = Tracker["pixercutoff"]
	zoom             = Tracker["zoom"]
	center           = Tracker["constants"]["center"]
	CTF              = Tracker["constants"]["CTF"]
	ref_a            = Tracker["constants"]["ref_a"]
	rstep            = Tracker["constants"]["rs"]
	sym              = Tracker["constants"]["sym"]
	first_ring       = 1
	last_ring        = Tracker["radius"]
	xr               = Tracker["xr"]
	yr               = Tracker["yr"]
	ts               = Tracker["ts"]
	an               = Tracker["an"]
	delta            = Tracker["delta"]
	max_iter         = Tracker["maxit"]

	if mpi_comm == None: mpi_comm = MPI_COMM_WORLD

	if log == None:
		from sp_logger import Logger
		log = Logger()

	number_of_proc = mpi_comm_size(mpi_comm)
	myid           = mpi_comm_rank(mpi_comm)
	main_node      = 0

	if myid == main_node:
		log.add("Start sali3d_base_h_01, nsoft = %1d"%nsoft)

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

	if( type(stack) is types.StringType ):
		if myid == main_node:
			total_nima = EMUtil.get_image_count( stack )
		else:
			total_nima = 0
		total_nima = wrap_mpi_bcast(total_nima, main_node, mpi_comm)
		list_of_particles = range(total_nima)
		image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
		# create a list of images for each node
		list_of_particles = list_of_particles[image_start: image_end]
		nima = len(list_of_particles)

	else:
		list_of_particles = range(len(stack))
		nima = len(list_of_particles)
		total_nima = len(list_of_particles)
		total_nima = mpi_reduce(total_nima, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
		total_nima = mpi_bcast(total_nima, 1, MPI_INT, 0, MPI_COMM_WORLD)
		total_nima = int(total_nima[0])


	if myid == 0 or myid == 5 or myid == 14:
		# finfo = None
		import os
		outdir = "./"
		# info_file = os.path.join(outdir, "progress%04d"%myid)
		
		# finfo = open(info_file, 'w')
		if "filename_first_time" not in sali3d_base_h_01.__dict__:
			import datetime
			sali3d_base_h_01.filename_first_time = os.path.join(outdir, "progress%04d_%s"%(myid, datetime.datetime.now().strftime('%Y-%m-%d--%H-%M-%S-%f')[:-3]))
			finfo = open(sali3d_base_h_01.filename_first_time, 'w')
		else: 
			finfo = open(sali3d_base_h_01.filename_first_time, 'a')
			finfo.write("=======================================================================================\n")
		# finfo = open(info_file, 'a')
	else:
		finfo = None

	if( myid == main_node):
		if( type(stack) is types.StringType ):  mask2D = get_im(stack, list_of_particles[0])
		else:                                   mask2D = stack[list_of_particles[0]]
		nx = mask2D.get_xsize()
	else:  nx = 0
	nx  = bcast_number_to_all(nx, source_node = main_node)
	mx = 2*nx
	if last_ring < 0:	last_ring = int(nx/2) - 2

	numr	= Numrinit(first_ring, last_ring, rstep, "F")

	data = [None]*nima
	for im in xrange(nima):
		if( type(stack) is types.StringType ):  data[im] = get_im(stack, list_of_particles[im])
		else:                                   data[im] = stack[list_of_particles[im]]
	mpi_barrier(mpi_comm)


	if myid == main_node:
		start_time = time()

	#  Read	template volume if provided or reconstruct it
	"""
	user_func = Tracker["constants"]["user_func"]
	if ref_vol:
		#vol = do_volume_mrk01(ref_vol, Tracker, 0, mpi_comm)
		ref_data = [ref_vol, Tracker, 0, mpi_comm]
		vol = user_func(ref_data)
	else:
		#vol = do_volume_mrk01(data, Tracker, 0, mpi_comm)
		ref_data = [data, Tracker, 0, mpi_comm]
		vol = user_func(ref_data)
	"""
	vol = ref_vol
	# log
	if myid == main_node:
		log.add("Setting of reference 3D reconstruction time = %10.1f\n"%(time()-start_time))
		start_time = time()


	pixer = [0.0]*nima
	#historyofchanges = [0.0, 0.5, 1.0]
	#par_r = [[] for im in list_of_particles ]
	cs = [0.0]*3
	total_iter = 0
	# do the projection matching,  it has a loop over iterations here, 
	#  but it can only do one iteration as many settings are done in meridien.  Perturbations are a good example, there is one per each iteration.
	if zoom: lstp = 1
	
	for N_step in xrange(lstp):
		# calculate_number_of_cones(volft, kb, delta, sym, cnx, cny, numr, mode, wr_four)
		cnx = cny = nx//2 + 1
		# numr = Numrinit(1,15)
		mode = "F"
		# wr_four = ringwe(numr, mode)
		# volft, kb = prep_vol(vol)

		terminate = 0
		Iter = 0
		while Iter < max_iter and terminate == 0:

			Iter += 1
			total_iter += 1

			mpi_barrier(mpi_comm)
			if myid == main_node:
				log.add("ITERATION #%3d,  inner iteration #%3d"%(total_iter, Iter))
				log.add("Delta = %7.4f, an = %7.4f, xrange = %7.4f, yrange = %7.4f, step = %7.4f   %7.4f  %7.4f\n"%\
							(delta[N_step], an[N_step], xrng[N_step], yrng[N_step], step[N_step], rshift, rangle))
				start_time = time()


			#=========================================================================
			# build references
			volft, kb = prep_vol(vol)
			projangles = [[] for i in xrange(nima)]
			for im in xrange(nima):
				projangles[im] = get_params_proj(data[im])[:3]
					

			cone_count = 0
			# start loop here ZCcw2oL8ZbcGbU		
			for image_indices, refrings, list_of_reference_angles in generate_indices_and_refrings(nima, projangles, volft, kb, nx, delta[N_step], an[N_step],	rangle, ref_a, sym, numr, MPI=mpi_comm, phiEqpsi = "Zero"):
				
				cone_count += 1
				print "cone_count", cone_count

				if myid == main_node:
					log.add("Time to prepare rings: %10.1f\n" % (time()-start_time))
					start_time = time()
	
				#=========================================================================
				#  there is no need for previousmax for deterministic searches
				if total_iter == 1 and nsoft > 0:
					if(an[N_step] < 0.0):
						# adjust params to references, calculate psi+shifts, calculate previousmax
						# for im in xrange(nima):
						for im in image_indices:
							previousmax = data[im].get_attr_default("previousmax", -1.0e23)
							if(previousmax == -1.0e23):
								peak, pixer[im] = proj_ali_incore_local(data[im], refrings, list_of_reference_angles, numr, \
										xrng[N_step], yrng[N_step], step[N_step], delta[N_step]*2.5, sym = sym)
								data[im].set_attr("previousmax", peak)
					else:
						#  Here it is supposed to be shake and bake for local SHC, but it would have to be signaled somehow
						for im in xrange(nima):
							data[im].set_attr("previousmax", -1.0e23)
					if myid == main_node:
						log.add("Time to calculate first psi+shifts+previousmax: %10.1f\n" % (time()-start_time))
						start_time = time()
				#=========================================================================
	
				# cannot have barriers in this loop because some cones might not have images assigned to them! In this case 'image_indices' is empty.
				# mpi_barrier(mpi_comm)
				if myid == main_node:  start_time = time()
				#=========================================================================
				# alignment
				#number_of_checked_refs = 0
				par_r = [0]*max(2,(nsoft+1))
				if(an[N_step] > 0):
					pass
					# these are already calculated by the generator at the top of the loop
					# generate list of angles
					# from alignment import generate_list_of_reference_angles_for_search
					# list_of_reference_angles = \
					# generate_list_of_reference_angles_for_search([[refrings[lr].get_attr("phi"), refrings[lr].get_attr("theta")] for lr in xrange(len(refrings))], sym=sym)			
				else:  list_of_reference_angles = [[1.0,1.0]]
				error_status = 0
				from sp_utilities import if_error_then_all_processes_exit_program

				# for im in xrange(nima):
				for im in image_indices:
					if Tracker["constants"]["pwsharpening"] :
						#  High-pass filtration of data[im]
						try:
							stmp = data[im].get_attr("ptcl_source_image")
						except:
							try:
								stmp = data[im].get_attr("ctf")
								stmp = round(stmp.defocus,4)
							except:
								ERROR("Either ptcl_source_image or ctf has to be present in the header.","meridien",1, myid)
						try:
							indx = Tracker["bckgnoise"][1].index(stmp)
						except:
							ERROR("Problem with indexing ptcl_source_image.","meridien",1, myid)
	
						tempdata = Util.window(pad(filt_table(data[im],[Tracker["bckgnoise"][0][i,indx] for i in xrange(nx)]), mx, mx,1,0.0), nx, nx)
					else:  tempdata = data[im].copy()
					if(nsoft == 0):
						if(an[N_step] == -1):
							#  In zoom option each projection goes through shift zoom alignment
							if  zoom: peak, pixer[im] = proj_ali_incore_zoom(tempdata, refrings, numr, \
															xrng, yrng, step, finfo = finfo, sym=sym)
							else:  peak, pixer[im] = proj_ali_incore(tempdata, refrings, numr, \
													xrng[N_step], yrng[N_step], step[N_step], finfo = finfo, sym=sym, delta_psi = delta[N_step], rshift = rshift*xrng[N_step])
						else:
							if  zoom: peak, pixer[im] = proj_ali_incore_local_zoom(tempdata, refrings, list_of_reference_angles, numr, \
										xrng, yrng, step, an, finfo = finfo, sym=sym)
							else:  
								
								peak, pixer[im] = proj_ali_incore_local(tempdata, refrings, list_of_reference_angles, numr, \
										xrng[N_step], yrng[N_step], step[N_step], an[N_step], finfo = finfo, sym=sym, delta_psi = delta[N_step], rshift = rshift)
						if(pixer[im] == 0.0):  par_r[0] += 1
					elif(nsoft == 1):
						tempdata.set_attr("previousmax", data[im].get_attr("previousmax"))
						peak, pixer[im], number_of_checked_refs, iref = \
							shc(tempdata, refrings, list_of_reference_angles, numr, xrng[N_step], yrng[N_step], step[N_step], an[N_step], sym, finfo = finfo)
						if(pixer[im] == 0.0):  par_r[0] += 1
						data[im].set_attr("previousmax", tempdata.get_attr("previousmax"))
					elif(nsoft > 1):
						#  This is not functional
						peak, pixer[im], checked_refs, number_of_peaks = shc_multi(data[im], refrings, numr, \
													xrng[N_step], yrng[N_step], step[N_step], an[N_step], nsoft, sym, finfo = finfo)
						par_r[number_of_peaks] += 1
						#number_of_checked_refs += checked_refs
					data[im].set_attr("xform.projection", tempdata.get_attr("xform.projection"))
				if len(image_indices)>0: del tempdata
				if(an[N_step] > 0):  del list_of_reference_angles
				#=========================================================================
				# if_error_then_all_processes_exit_program(error_status)
				# cannot have barriers in this loop because some cones might not have images assigned to them! In this case 'image_indices' is empty.
				# mpi_barrier(mpi_comm)
				if myid == main_node:
					#print  data[0].get_attr_dict()
					log.add("Time of alignment = %10.1f\n"%(time()-start_time))
					start_time = time()
				
				# end loop here ZCcw2oL8ZbcGbU
				
			del volft, kb
			mpi_barrier(mpi_comm)
			print_from_process(0, "passed")
			#=========================================================================
			#  Pixer errors available here are useless as they are done for shifts on the reduced image size.
			#output pixel errors, check stop criterion
			all_pixer = wrap_mpi_gatherv(pixer, 0, mpi_comm)
			par_r = mpi_reduce(par_r, len(par_r), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
			#total_checked_refs = wrap_mpi_gatherv([number_of_checked_refs], main_node, mpi_comm)
			terminate = 0
			if myid == main_node:
				#total_checked_refs = sum(total_checked_refs)
				if(nsoft < 2):  par_r[1] = total_nima - par_r[0]
				log.add("=========== Number of better orientations found ==============")
				for lhx in xrange(len(par_r)):
					msg = "            %5d     %7d"%(lhx, par_r[lhx])
					log.add(msg)
				log.add("_______________________________________________________")
				changes = par_r[0]/float(total_nima)
				"""  Have to think about it PAP
				if(  changes > saturatecrit ):
					if( Iter == 1 ):
						log.add("Will continue even though %4.2f images did not find better orientations"%saturatecrit)
					else:
						terminate = 1
						log.add("...............")
						log.add(">>>>>>>>>>>>>>>   Will terminate as %4.2f images did not find better orientations"%saturatecrit)
				if( terminate == 0 ):
					historyofchanges.append(changes)
					historyofchanges = historyofchanges[:3]
					historyofchanges.sort()
					if( (historyofchanges[-1]-historyofchanges[0])/2/(historyofchanges[-1]+historyofchanges[0]) <0.05 ):
						terminate = 1
						log.add("...............")
						log.add(">>>>>>>>>>>>>>>   Will terminate as orientations do not improve anymore")
				"""

				lhist = 20
				region, histo = hist_list(all_pixer, lhist)
				log.add("=========== Histogram of pixel errors ==============")
				for lhx in xrange(lhist):
					msg = "          %10.3f     %7d"%(region[lhx], histo[lhx])
					log.add(msg)
				log.add("____________________________________________________")
				if(nsoft<2 and terminate == 0):
					lhx = 0
					for msg in all_pixer:
						if(msg < pixercutoff): lhx += 1
					lhx = float(lhx)/float(total_nima)
					log.add(">>> %4.2f images had pixel error <%5.2f"%(lhx,pixercutoff))
					"""
					if( lhx > saturatecrit):
						if( Iter == 1 ):
							log.add("Will continue even though %4.2f images had pixel error < %5.2f"%(saturatecrit,pixercutoff))
						else:
							terminate = 1
							log.add("...............")
							log.add(">>>>>>>>>>>>>>>   Will terminate as %4.2f images had pixel error < %5.2f"%(saturatecrit,pixercutoff))
					"""
			terminate = True  #wrap_mpi_bcast(terminate, main_node, mpi_comm)
			#=========================================================================
			mpi_barrier(mpi_comm)
			if myid == main_node:
				#print  data[0].get_attr_dict()
				log.add("Time to compute histograms = %10.1f\n"%(time()-start_time))
				start_time = time()


			#=========================================================================
			mpi_barrier(mpi_comm)
			if( terminate or (Iter == max_iter) ):
				# gather parameters
				params = []
				for im in xrange(nima):
					t = get_params_proj(data[im])
					params.append( [t[0], t[1], t[2], t[3], t[4]] )
				params = wrap_mpi_gatherv(params, main_node, mpi_comm)
			# centering and volume reconstruction if not terminating
			else:
				#=========================================================================
				# centering
				if center == -1 and sym[0] == 'c':
					from sp_utilities      import estimate_3D_center_MPI, rotate_3D_shift
					cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node, mpi_comm=mpi_comm)
					if myid == main_node:
						msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
						log.add(msg)
					if int(sym[1]) > 1:
						cs[0] = cs[1] = 0.0
						if myid == main_node:
							log.add("For symmetry group cn (n>1), we only center the volume in z-direction\n")
					cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, mpi_comm)
					cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
					rotate_3D_shift(data, cs)
				#=========================================================================
				if myid == main_node:
					start_time = time()
				#vol = do_volume_mrk01(data, Tracker, total_iter, mpi_comm)
				ref_data = [data, Tracker, total_iter, mpi_comm]
				user_func = Tracker["constants"] ["user_func"]
				vol = user_func(ref_data)
				#if myid == main_node:  vol.write_image('soft/smvol%04d.hdf'%total_iter)
				# log
				if myid == main_node:
					log.add("3D reconstruction time = %10.1f\n"%(time()-start_time))
					start_time = time()
			#=========================================================================

			"""
			#=========================================================================
			if(False):  #total_iter%1 == 5 or terminate):
				# gather parameters
				params = []
				previousmax = []
				for im in data:
					t = get_params_proj(im)
					params.append( [t[0], t[1], t[2], t[3]/shrinkage, t[4]/shrinkage] )
					#if(t[3] >0.0 or t[4]>0.0):  print  "  ERRROR  ",t
					previousmax.append(im.get_attr("previousmax"))
				assert(nima == len(params))
				params = wrap_mpi_gatherv(params, 0, mpi_comm)
				if myid == 0:
					assert(total_nima == len(params))
				previousmax = wrap_mpi_gatherv(previousmax, 0, mpi_comm)
				if myid == main_node:
					from sp_utilities import write_text_row, write_text_file
					write_text_row(params, "soft/params%04d.txt"%total_iter)
					write_text_file(previousmax, "soft/previousmax%04d.txt"%total_iter)


				del previousmax, params
				i = 1
				while data[0].has_attr("xform.projection" + str(i)):
					params = []
					previousmax = []
					for im in data:

						try:
							#print  im.get_attr("xform.projection" + str(i))
							t = get_params_proj(im,"xform.projection" + str(i))
						except:
							print " NO XFORM  ",myid, i,im.get_attr('ID')
							from sys import exit
							exit()

						params.append( [t[0], t[1], t[2], t[3]/shrinkage, t[4]/shrinkage] )
					assert(nima == len(params))
					params = wrap_mpi_gatherv(params, 0, mpi_comm)
					if myid == 0:
						assert(total_nima == len(params))
					if myid == main_node:
						write_text_row(params, "soft/params-%04d-%04d.txt"%(i,total_iter))
					del previousmax, params
					i+=1


			if( ( terminate or (Iter == max_iter) ) and (myid == main_node) ):
				if( type(stack) is types.StringType ):
					from EMAN2 import Vec2f, Transform
					from EMAN2db import db_open_dict
					DB = db_open_dict(stack)
					for im in xrange(len(params)):
						t = Transform({"type":"spider","phi":params[im][0],"theta":params[im][1],"psi":params[im][2]})
						t.set_trans(Vec2f(-params[im][3], -params[im][4]))
						DB.set_attr(particle_ids[im], "xform.projection", t)
					DB.close()
				else:
					for im in xrange(len(params)): set_params_proj(stack[particle_ids[im]], params[im])
			"""


	if myid == main_node:
		log.add("Finish sali3d_base, nsoft = %1d"%nsoft)
	return params
'''
























































"""10
	if myid == main_node:
		os.mkdir(outdir)
		import sp_global_def
		sp_global_def.LOGFILE =  os.path.join(outdir, sp_global_def.LOGFILE)
		print_begin_msg("local_ali3d_MPI")
		import sp_user_functions
		user_func = sp_user_functions.factory[user_func_name]
		if CTF:
			ima = EMData()
			ima.read_image(stack, 0)
			ctf_applied = ima.get_attr_default("ctf_applied", 0)
			del ima
			if ctf_applied == 1:  ERROR("local_ali3d does not work for CTF-applied data", "local_ali3d_MPI", 1,myid)
	mpi_barrier(mpi_comm)
	"""












































































"""11
	if myid == main_node:
		import sp_user_functions
		user_func = sp_user_functions.factory[user_func_name]

		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Maskfile                    : %s\n"%(Tracker["constants"]["mask3D"]))
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Angular search range        : %s\n"%(delta))
		print_msg("Shift search range          : %f\n"%(ts))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Symmetry group              : %s\n"%(sym))
		print_msg("Chunk size                  : %f\n\n"%(chunk))
		print_msg("User function               : %s\n"%(user_func_name))
	"""


































































































































































































































































"""12
			if( lhx > saturatecrit):
				if( iteration == 1 ):
					log.add("First iteration, will continue even though %4.2f images did not find better orientations"%saturatecrit)
				else:
					terminate = 1
					log.add("..............................................................")
					log.add(">>>>>>>>>>>>>>>   Will terminate as %4.2f images did not find better orientations"%saturatecrit)
			"""
















"""13
		if( type(stack) is types.StringType ):
			from EMAN2 import Vec2f, Transform
			from EMAN2db import db_open_dict
			DB = db_open_dict(stack)
			for im in xrange(len(params)):
				t = Transform({"type":"spider","phi":params[im][0],"theta":params[im][1],"psi":params[im][2]})
				t.set_trans(Vec2f(-params[im][3], -params[im][4]))
				DB.set_attr(particle_ids[im], "xform.projection", t)
			DB.close()
		else:
			for im in xrange(len(params)): set_params_proj(stack[particle_ids[im]], params[im])

		log.add("Time to write header information= %d\n"%(time()-start_time))
		"""





'''14
#  This does not handle symmetries properly

def ali3dlocal_MPI(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", apsi = "-1", deltapsi = "-1", startpsi = "-1",
	    center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",  user_func_name = "ref_ali3d",
	    fourvar = True, npad = 4, debug = False, termprec = 0.0):

	from sp_alignment       import Numrinit, ringwe, proj_ali_incore, proj_ali_incore_local, proj_ali_incore_local_psi
	from sp_utilities       import model_circle, get_image, drop_image, get_input_from_string
	from sp_utilities       import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all
	from sp_utilities       import send_attr_dict
	from sp_utilities       import get_params_proj, file_type
	from sp_fundamentals    import rot_avg_image
	from sp_utilities       import print_begin_msg, print_end_msg, print_msg
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi             import mpi_reduce, MPI_INT, MPI_SUM
	from sp_filter          import filt_ctf
	from sp_projection      import prep_vol, prgs
	from sp_statistics      import hist_list, varf3d_MPI
	from sp_applications    import MPI_start_end
	import os
	import types


	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ali3dlocal_MPI", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		os.mkdir(outdir)
		import sp_global_def
		sp_global_def.LOGFILE =  os.path.join(outdir, sp_global_def.LOGFILE)
		print_begin_msg("ali3dlocal_MPI")
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

	if apsi == "-1":
		apsi = [-1] * lstp
	else:
		apsi = get_input_from_string(apsi)

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
		import sp_user_functions
		user_func = sp_user_functions.factory[user_func_name]

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
		print_msg("Angular search range (phi and theta)       : %s\n"%(an))
		#print_msg("Angular search range (psi)                 : %s\n"%(apsi))
		#print_msg("Delta psi                   : %s\n"%(deltapsi))
		#print_msg("Start psi                   : %s\n"%(startpsi))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Percentage of change for termination: %f\n"%(termprec))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Reference projection method : %s\n"%(ref_a))
		print_msg("Symmetry group              : %s\n\n"%(sym))
		print_msg("User function               : %s\n"%(user_func_name))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else: mask3D = model_circle(last_ring, nx, nx, nx)

	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	wr_four = ringwe(numr, "F")
	cx = cy = nx//2
	mask2D  = model_circle(last_ring,nx,nx) - model_circle(first_ring,nx,nx)

	fscmask = mask3D  #model_circle(last_ring,nx,nx,nx)  For a fancy mask circle would work better  PAP 7/21/11
	if CTF:
		from sp_reconstruction import rec3D_MPI
		from sp_filter         import filt_ctf
	else:	 from sp_reconstruction import rec3D_MPI_noCTF

	if myid == main_node:
		if file_type(stack) == "bdb":
			from EMAN2db import db_open_dict
			dummy = db_open_dict(stack, True)
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = EMUtil.get_all_attributes(stack, 'active')
		# list_of_particles = []
		# for im in xrange(len(active)):
		# 	if active[im]:  list_of_particles.append(im)
		# del active
		# nima = len(list_of_particles)

		nima = EMUtil.get_image_count(stack)
		list_of_particles = range(nima)
		
	else:
		nima = 0
	total_nima = bcast_number_to_all(nima, source_node = main_node)

	if myid != main_node:
		list_of_particles = [-1]*total_nima
	list_of_particles = bcast_list_to_all(list_of_particles, myid,  source_node = main_node)

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
		if im == main_node :  disps.append(0)
		else:                 disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append(ie-ib)

	#  Number of reference image that fit into the memory.  Will have to be hardwired.
	numberofrefs = 500


	pixer = [0.0]*nima
	cs = [0.0]*3
	total_iter = 0
	# do the projection matching
	for N_step in xrange(lstp):
		terminate = 0
		Iter = -1
 		while Iter < max_iter-1 and terminate == 0:
			Iter += 1
			total_iter += 1
			if myid == main_node:
				start_time = time()
				print_msg("\nITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f, delta psi = %5.2f, start psi = %5.2f\n"%\
						(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step],deltapsi[N_step],startpsi[N_step]))

			totrefs = len(even_angles(delta[N_step], method = ref_a, symmetry = sym))
			numberofcones = max(totrefs/numberofrefs,1)
			if myid == main_node:
				print_msg("\n   Number of references permitted in memory = %d  , total number of references = %d , number of cones = %d \n"%(numberofrefs, totrefs, numberofcones))


			volft, kb = prep_vol(vol)
			if( numberofcones == 1):
				# One cone good enough, use the original code
				refrings = prepare_refrings(volft, kb, nx, delta[N_step], ref_a, sym, numr, True)

				if myid == main_node:
					print_msg("Time to prepare rings: %d\n" % (time()-start_time))
					start_time = time()

				for im in xrange(nima):
					#  It needs list_of_reference_angles
					peak, pixer[im] = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step],finfo, sym = sym)
					data[im].set_attr("previousmax", peak)

			else:
				from sp_morphology import  bracket_def
				from sp_utilities  import  assign_projangles, cone_ang
				from sp_alignment  import  refprojs

				h = 1.0
				dat = [sym, numberofcones, ref_a]
				def1, def2 = bracket_def(computenumberofrefs,dat, rs, h)
				def1, val  = goldsearch_astigmatism(computenumberofrefs, dat, def1, def2, tol=1.0)
				coneangles = even_angles(def1, method = "S", symmetry = sym)
				if myid == main_node:
					print_msg("\n   Computed cone delta = %f  , and the number of cones = %d \n"%(def1, len(coneangles)))
				assignments = assign_projangles(projangles, coneangles)
				for k in xrange(len(coneangles)):
					if(len(assignements[k]) > 0):
						refsincone = even_angles(delta, method = ref_a, symmetry = sym)
						ant = 1.5*an[N_step]
						refsincone = cone_ang( refsincone, coneangles[k][0], coneangles[k][1], ant )
						refrings = refprojs( volf, kb, refsincone, cnx, cny, numr, "F", wr_four )
						#    match projections to its cone using an as a distance.
						for im in assignments[k]:
							#  It needs list_of_reference_angles
							peak, pixer[im] = proj_ali_incore_local(data[im], refrings, numr, xrng[N_step], yrng[N_step], step[N_step], an[N_step], finfo, sym = sym)
							data[im].set_attr("previousmax", peak)

			del volft, kb

			if myid == main_node:
				print_msg("Time of alignment = %d\n"%(time()-start_time))
				start_time = time()

			#output pixel errors
			from mpi import mpi_gatherv
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			terminate = 0
			if myid == main_node:
				recvbuf = map(float, recvbuf)
				from sp_statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				if region[0] < 0.0:  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				# Terminate if 95% within 1 pixel error
				im = 0
				for lhx in xrange(lhist):
					if region[lhx] > 1.0: break
					im += histo[lhx]
				precn = 100*float(total_nima-im)/float(total_nima)
				msg = " Number of particles that changed orientations %7d, percentage of total: %5.1f\n"%(total_nima-im, precn)
				print_msg(msg)
				if precn <= termprec:  terminate = 1
				del region, histo
			del recvbuf
			terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
			terminate = int(terminate[0])

			if center == -1 and sym[0] == 'c':
				from sp_utilities      import estimate_3D_center_MPI, rotate_3D_shift
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					print_msg(msg)
				if int(sym[1]) > 1:
					cs[0] = cs[1] = 0.0
					if myid == main_node:
						print_msg("For symmetry group cn (n>1), we only center the volume in z-direction\n")
				cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(data, cs)

			# write out headers, under MPI writing has to be done sequentially
			mpi_barrier(MPI_COMM_WORLD)
			par_str = ['xform.projection', 'previousmax', 'ID']
			if myid == main_node:
	   			if(file_type(stack) == "bdb"):
	        			from sp_utilities import recv_attr_dict_bdb
	        			recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        		else:
	        			from sp_utilities import recv_attr_dict
	        			recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
				print_msg("Time to write header information= %d\n"%(time()-start_time))
				start_time = time()
	        	else:	       send_attr_dict(main_node, data, par_str, image_start, image_end)

			if CTF: vol, fscc = rec3D_MPI(data, snr, sym, fscmask, os.path.join(outdir, "resolution%04d"%(total_iter)), myid, main_node, npad = npad)
			else:   vol, fscc = rec3D_MPI_noCTF(data, sym, fscmask, os.path.join(outdir, "resolution%04d"%(total_iter)), myid, main_node, npad = npad)

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


	if myid == main_node: print_end_msg("ali3dlocal_MPI")
'''









































































































































































































































































'''15
def Xali3d_shc0MPI(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an = "-1", apsi = "-1", deltapsi = "-1", startpsi = "-1",
	    center = -1, maxit = 5, CTF = False, snr = 1.0,  ref_a = "S", sym = "c1",  user_func_name = "ref_ali3d",
	    fourvar = True, npad = 4, debug = False, termprec = 0.0, gamma=-1):

	from sp_alignment       import Numrinit, prepare_refrings, proj_ali_incore, proj_ali_incore_local, proj_ali_incore_local_psi, shc
	from sp_utilities       import model_circle, get_image, drop_image, get_input_from_string
	from sp_utilities       import bcast_list_to_all, bcast_number_to_all, bcast_EMData_to_all
	from sp_utilities       import send_attr_dict, get_params_proj, file_type
	import os
	import types
	from sp_utilities       import print_begin_msg, print_end_msg, print_msg
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier, MPI_INT
	from sp_projection      import prep_vol, prgs
	from sp_statistics      import hist_list, varf3d_MPI
	from sp_applications    import MPI_start_end
	from math            import sqrt, acos, radians
	from random          import shuffle

	if gamma > 0:
		gamma = radians(gamma)

	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	
	if os.path.exists(outdir):  ERROR('Output directory exists, please change the name and restart the program', "ali3d_MPI", 1, myid)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		os.mkdir(outdir)
		import sp_global_def
		sp_global_def.LOGFILE =  os.path.join(outdir, sp_global_def.LOGFILE)
		print_begin_msg("ali3d_shcMPI")
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

	if apsi == "-1":
		apsi = [-1] * lstp
	else:
		apsi = get_input_from_string(apsi)

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
		import sp_user_functions
		user_func = sp_user_functions.factory[user_func_name]

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
		print_msg("Angular search range (phi and theta)       : %s\n"%(an))
		#print_msg("Angular search range (psi)                 : %s\n"%(apsi))
		#print_msg("Delta psi                   : %s\n"%(deltapsi))
		#print_msg("Start psi                   : %s\n"%(startpsi))
		print_msg("Maximum number of iterations : %i\n"%(max_iter))
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
		from sp_reconstruction import rec3D_MPI
		from sp_filter         import filt_ctf
	else:	 from sp_reconstruction import rec3D_MPI_noCTF

	if myid == main_node:
		if file_type(stack) == "bdb":
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
	list_of_particles = bcast_list_to_all(list_of_particles, myid,  source_node = main_node)

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
		if im == main_node :  disps.append(0)
		else:                 disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append(ie-ib)

	
	mode = "F"
	pixer = [0.0]*nima
	cs = [0.0]*3
	total_iter = 0
	final_params = None
	final_volume = None
	final_volume_filtered = None
	# do the projection matching
	for N_step in xrange(lstp):
		##convert data into polar coordinates.  @ming.
		ky = int(2*yrng[N_step]/step[N_step]+0.5)/2;
		kx = int(2*xrng[N_step]/step[N_step]+0.5)/2;
		cimages = []
		for im in xrange(nima):
			nx = data[im].get_xsize()
			ny = data[im].get_ysize()
			#  center is in SPIDER convention
			cnx  = nx//2 + 1
			cny  = ny//2 + 1
			ims = [None]*(2*ky+1)*(2*kx+1)
			
			nring = len(numr)/3
			inr = numr[3*(nring-1)]
			for i in xrange(-ky, ky+1):
				iy = i*step[N_step]
				if inr+int(cny+iy) <= ny-1 and -inr + int(cny+iy) >=1:
					for j in xrange(-kx, kx+1):
						ix = j*step[N_step]
						if inr+int(cnx+ix) <= nx-1 and -inr + int(cnx+ix) >=1:
							cimage = Util.Polar2Dm(data[im], cnx+ix, cny+iy, numr, mode)
							Util.Normalize_ring( cimage, numr, 0)
							Util.Frngs(cimage, numr)
							ims[(i+ky)*(2*kx+1)+j+kx] = cimage
			cimages.append(ims)
			
	
		terminate = 0
		Iter = 0
		while Iter < max_iter and terminate == 0:

			Iter += 1
			total_iter += 1

			if myid == main_node:
				start_time = time()
				print_msg("\nITERATION #%3d,  inner iteration #%3d\nDelta = %5.2f, an = %5.2f, xrange = %5.2f, yrange = %5.2f,translational step = %5.2f\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step], yrng[N_step], step[N_step]))
				#print_msg("\nITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f, delta psi = %5.2f, start psi = %5.2f\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step],deltapsi[N_step],startpsi[N_step]))

			#=========================================================================
			# build references
			volft, kb = prep_vol(vol)
			refrings = prepare_refrings(volft, kb, nx, delta[N_step], ref_a, sym, numr, True)
			lastdelta = delta[N_step]
			del volft, kb
			#=========================================================================

			if myid == main_node:
				print_msg("Time to prepare rings: %d\n" % (time()-start_time))
				start_time = time()

			#=========================================================================
			#  We assume previousmax exists
			"""
			if total_iter == 1:
				# adjust params to references, calculate psi+shifts, calculate previousmax
				for im in xrange(nima):
					stable = data[im].get_attr_default("stable", 0)
					if stable == 0:
						data[im].set_attr("previousmax", -1.0e23)
						data[im].set_attr("stable", 1)
					else:
						#print "  params  ",get_params_proj(data[im])
						peak, pixer[im] = proj_ali_incore_local(data[im],refrings,numr,0.0,0.0,1.0,delta[N_step]/4.0,finfo)
						data[im].set_attr("previousmax", peak)
						print "peak  ",im,peak,get_params_proj(data[im])
				if myid == main_node:
					print_msg("Time to calculate first psi+shifts+previousmax: %d\n" % (time()-start_time))
					start_time = time()
			"""
			#=========================================================================

			#=========================================================================
			# alignment
			iter_indexes = range(nima)
			shuffle(iter_indexes)
			for im in iter_indexes:
				from sp_utilities import get_params_proj
				#print "  IN  ",im,get_params_proj(data[im]),data[im].get_attr("previousmax")
				peak, pixer[im], number_of_checked_refs, iref = \
					shc0(data[im], cimages[im], refrings, numr, xrng[N_step], yrng[N_step], step[N_step], an[N_step], sym, finfo)
				#print "  OU  ",im,get_params_proj(data[im]),data[im].get_attr("previousmax")
				if gamma > 0:
					n1 = refrings[iref].get_attr("n1")
					n2 = refrings[iref].get_attr("n2")
					n3 = refrings[iref].get_attr("n3")
					to_be_deleted = []
					for irr in xrange(len(refrings)):
						nn1 = refrings[irr].get_attr("n1")
						nn2 = refrings[irr].get_attr("n2")
						nn3 = refrings[irr].get_attr("n3")
						if abs(acos(n1*nn1 + n2*nn2 * n3*nn3)) < gamma:
							to_be_deleted.append( irr )
					if len(to_be_deleted) > 0:
						to_be_deleted.sort(reverse=True)
						for irr in to_be_deleted:
							del refrings[irr]
				elif gamma == 0:
						del refrings[iref]
			#=========================================================================
			del refrings

			if myid == main_node:
				print_msg("Time of alignment = %d\n"%(time()-start_time))
				start_time = time()
			#=========================================================================
			#output pixel errors, check stop criterion
			from mpi import mpi_gatherv
			recvbuf = mpi_gatherv(pixer, nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, main_node, MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			terminate = 0
			if myid == main_node:
				recvbuf = map(float, recvbuf)
				from sp_statistics import hist_list
				lhist = 20
				region, histo = hist_list(recvbuf, lhist)
				if region[0] < 0.0:  region[0] = 0.0
				msg = "      Histogram of pixel errors\n      ERROR       number of particles\n"
				print_msg(msg)
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d\n"%(region[lhx], histo[lhx])
					print_msg(msg)
				# Terminate if 95% within 0.1 pixel error
				im = 0
				for lhx in xrange(lhist):
					if region[lhx] > 0.1: break
					im += histo[lhx]
				precn = 100*float(total_nima-im)/float(total_nima)
				msg = " Number of particles that changed orientations %7d, percentage of total: %5.1f\n"%(total_nima-im, precn)
				print_msg(msg)
				if precn <= termprec:  terminate = 1
				del region, histo
			del recvbuf
			terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
			terminate = int(terminate[0])
			#=========================================================================

			#=========================================================================
			# centering
			if center == -1 and sym[0] == 'c':
				from sp_utilities      import estimate_3D_center_MPI, rotate_3D_shift
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data, total_nima, myid, number_of_proc, main_node)
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					print_msg(msg)
				if int(sym[1]) > 1:
					cs[0] = cs[1] = 0.0
					if myid == main_node:
						print_msg("For symmetry group cn (n>1), we only center the volume in z-direction\n")
				cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, MPI_COMM_WORLD)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(data, cs)
			#=========================================================================

			#=========================================================================
			# write out headers, under MPI writing has to be done sequentially
			mpi_barrier(MPI_COMM_WORLD)
			# It takes plenty of times, so do it once in awhile...
			if( True ): #total_iter>3 and total_iter%5 == 0 ):
				par_str = ['xform.projection', 'previousmax', 'ID']
				if myid == main_node:
					if(file_type(stack) == "bdb"):
						from sp_utilities import recv_attr_dict_bdb
						recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
					else:
						from sp_utilities import recv_attr_dict
						recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
					"""
					# save parameters to file
					paro = [None]*total_nima
					projs_headers = EMData.read_images(stack, range(total_nima), True)
					for im in xrange(total_nima):
						a1,a2,a3,a4,a5 = get_params_proj(projs_headers[im])
						previousmax = projs_headers[im].get_attr("previousmax")
						paro[im] = [a1,a2,a3,a4,a5,previousmax]
					from sp_utilities import write_text_row
					write_text_row(paro,os.path.join(outdir, "params%04d.txt"%(total_iter)))
					final_params = paro
					del projs_headers
					del paro
					# ------- end of saving parameters to file
					print_msg("Time to write header information= %d\n"%(time()-start_time))
					start_time = time()
					"""
				else:
					send_attr_dict(main_node, data, par_str, image_start, image_end)
			#=========================================================================

			#=========================================================================
			# volume reconstruction
			#vol_previous = vol
			if CTF: vol, fscc = rec3D_MPI(data, snr, sym, fscmask, os.path.join(outdir, "resolution%04d"%(total_iter)), myid, main_node, npad = npad)
			else:   vol, fscc = rec3D_MPI_noCTF(data, sym, fscmask, os.path.join(outdir, "resolution%04d"%(total_iter)), myid, main_node, npad = npad)
			# log
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
			else:
				varf = None
			# user functions + save volume
			if myid == main_node:
				#drop_image(vol, os.path.join(outdir, "vol%04d.hdf"%(total_iter)))
				ref_data[2] = vol
				ref_data[3] = fscc
				ref_data[4] = varf
				#  call user-supplied function to prepare reference image, i.e., center and filter it
				vol, cs = user_func(ref_data)
				drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(total_iter)))
				#print_msg("Euclidean distance between the current and the previous volume: " + str(sqrt(vol.cmp("SqEuclidean",vol_previous,{"mask":mask3D,"zeromask":0,"normto":0}))) + "\n")
				print_msg("L2 norm of the volume: " + str(vol.cmp("dot", vol, {"negative":0, "mask":mask3D})) + "\n")
			del varf
			# broadcast volume
			bcast_EMData_to_all(vol, myid, main_node)
			#=========================================================================
		del cimages

	if myid == main_node: 
		print_end_msg("ali3d_shcMPI")
'''






































































































































































































































"""16
			if total_iter == 1:
				# adjust params to references, calculate psi+shifts, calculate previousmax
				for im in xrange(nima):
					stable = data[im].get_attr_default("stable", 0)
					if stable == 0:
						data[im].set_attr("previousmax", -1.0e23)
						data[im].set_attr("stable", 1)
					else:
						#print "  params  ",get_params_proj(data[im])
						peak, pixer[im] = proj_ali_incore_local(data[im],refrings,numr,0.0,0.0,1.0,delta[N_step]/4.0,finfo)
						data[im].set_attr("previousmax", peak)
						print "peak  ",im,peak,get_params_proj(data[im])
				if myid == main_node:
					print_msg("Time to calculate first psi+shifts+previousmax: %d\n" % (time()-start_time))
					start_time = time()
			"""


































































































'''17
					# save parameters to file
					paro = [None]*total_nima
					projs_headers = EMData.read_images(stack, range(total_nima), True)
					for im in xrange(total_nima):
						a1,a2,a3,a4,a5 = get_params_proj(projs_headers[im])
						previousmax = projs_headers[im].get_attr("previousmax")
						paro[im] = [a1,a2,a3,a4,a5,previousmax]
					from sp_utilities import write_text_row
					write_text_row(paro,os.path.join(outdir, "params%04d.txt"%(total_iter)))
					final_params = paro
					del projs_headers
					del paro
					# ------- end of saving parameters to file
					print_msg("Time to write header information= %d\n"%(time()-start_time))
					start_time = time()
					'''

























































































































































































































































































































































































































































'''18
	if(myid == main_node):	
		total_nima = EMUtil.get_image_count(stack)
		list_of_particles = range(total_nima)
	
	else:
		total_nima =0

	total_nima = bcast_number_to_all(total_nima, source_node = main_node)

	if(myid != main_node):
		list_of_particles = [-1]*total_nima

	list_of_particles = bcast_list_to_all(list_of_particles, myid,  source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	'''




































































































































































































































"""19
		#  Trying to use ISAC code for EQ-Kmeans  PAP 03/21/2015
		if myid == main_node:

			for imrefa in xrange(numrefang):
				from sp_utilities import findall
				N = findall(imrefa, assigntorefa)
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

			d = empty( (numref, total_nima), dtype = float32)
			for ima in xrange(total_nima):
				for iref in xrange(numref):  d[iref][ima] = dtot[iref][N[ima]]
			id_list_long = Util.assign_groups(str(d.__array_interface__['data'][0]), numref, nima) # string with memory address is passed as parameters
			del d
			id_list = [[] for i in xrange(numref)]
			maxasi = total_nima/numref
			for i in xrange(maxasi*numref):
				id_list[i/maxasi].append(id_list_long[i])
			for i in xrange(total_nima%maxasi):
				id_list[id_list_long[-1]].append(id_list_long[maxasi*numref+i])
			for iref in xrange(numref):
				id_list[iref].sort()

			assignment = [0]*total_nima
			for iref in xrange(numref):
				for im in id_list[iref]: assignment[im] = iref
		else:
			assignment = [0]*total_nima
		mpi_barrier(MPI_COMM_WORLD)
		#belongsto = mpi_bcast(belongsto, nima, MPI_INT, main_node, MPI_COMM_WORLD)
		#belongsto = map(int, belongsto)
		"""


















































































"""20
		if myid == main_node:
			assignment = [0]*total_nima
			for iref in xrange(numref):
				for im in xrange(len(asi[iref])):
					assignment[asi[iref][im]] = iref
			del asi
		"""

'''21
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































































































































































































































































































































































































































































































































































































"""22
	if runtype=="REFINEMENT":
		par_str = ['xform.projection', 'ID', 'group']
	else:
		par_str = ['group', 'ID' ]
	if myid == main_node:
		from sp_utilities import file_type
		if file_type(stack) == "bdb":
			from sp_utilities import recv_attr_dict_bdb
			recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
		else:
			from sp_utilities import recv_attr_dict
			recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	else:		send_attr_dict(main_node, data, par_str, image_start, image_end)
	"""









































































































































































































































































"""23
			if CTF:
				previous_defocus = -1.0
				if runtype=="REFINEMENT":
					start_time = time()
					prjref = prgq( volft, kb, nx, delta[N_step], ref_a, sym, MPI=True)
					if myid == main_node:
						log.add( "Calculation of projections: %d" % (time()-start_time) );start_time = time()
					del volft, kb
			else:
				if runtype=="REFINEMENT":
					start_time = time()
					refrings = prepare_refrings( volft, kb, nx, delta[N_step], ref_a, syms[0], numr)
					if myid == main_node:
						log.add( "Initial time to prepare rings: %d" % (time()-start_time) );start_time = time()
					del volft, kb
			"""


"""24
				if CTF:
					ctf = data[im].get_attr( "ctf" )
					if runtype=="REFINEMENT":
						if ctf.defocus != previous_defocus:
							previous_defocus = ctf.defocus
							rstart_time = time()
							refrings = gen_rings_ctf( prjref, nx, ctf, numr)
							if myid == main_node:
								log.add( "Repeated time to prepare rings: %d" % (time()-rstart_time) );rstart_time = time()
				"""






































































































































































































"""25
	if runtype=="REFINEMENT":
		par_str = ['xform.projection', 'ID', 'group']
	else:
		par_str = ['group', 'ID' ]
	if myid == main_node:
		from sp_utilities import file_type
		if file_type(stack) == "bdb":
			from sp_utilities import recv_attr_dict_bdb
			recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
		else:
			from sp_utilities import recv_attr_dict
			recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	else:		send_attr_dict(main_node, data, par_str, image_start, image_end)
	"""










































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































"""26
	if myid == main_node:
		os.mkdir(outdir)
		import sp_global_def
		sp_global_def.LOGFILE =  os.path.join(outdir, sp_global_def.LOGFILE)
		print_begin_msg("local_ali3d_MPI")
		import sp_user_functions
		user_func = sp_user_functions.factory[user_func_name]
		if CTF:
			ima = EMData()
			ima.read_image(stack, 0)
			ctf_applied = ima.get_attr_default("ctf_applied", 0)
			del ima
			if ctf_applied == 1:  ERROR("local_ali3d does not work for CTF-applied data", "local_ali3d_MPI", 1,myid)
	mpi_barrier(mpi_comm)
	"""








































































































"""27
	if myid == main_node:
		import sp_user_functions
		user_func = sp_user_functions.factory[user_func_name]

		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Maskfile                    : %s\n"%(ali3d_options.mask3D))
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Angular search range        : %s\n"%(delta))
		print_msg("Shift search range          : %f\n"%(ts))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Symmetry group              : %s\n"%(sym))
		print_msg("Chunk size                  : %f\n\n"%(chunk))
		print_msg("User function               : %s\n"%(user_func_name))
	"""








































































































































































































































































"""28
		if( type(stack) is types.StringType ):
			from EMAN2 import Vec2f, Transform
			from EMAN2db import db_open_dict
			DB = db_open_dict(stack)
			for im in xrange(len(params)):
				t = Transform({"type":"spider","phi":params[im][0],"theta":params[im][1],"psi":params[im][2]})
				t.set_trans(Vec2f(-params[im][3], -params[im][4]))
				DB.set_attr(particle_ids[im], "xform.projection", t)
			DB.close()
		else:
			for im in xrange(len(params)): set_params_proj(stack[particle_ids[im]], params[im])

		log.add("Time to write header information= %d\n"%(time()-start_time))
		"""















































































































































































































































































































































































































































































































































































'''29
					if myid == main_node:
						print_msg("\nphi = %5.2f, theta = %5.2f, psi=%5.2f\n"%( refrings[i].get_attr('phi'), refrings[i].get_attr('theta'), refrings[i].get_attr('psi') ) )
					if myid == main_node:
						print_msg("\nlen(ref1) = %4d, len(ref2) = %4d\n"%(len(refrings1), len(refrings2)) )
					'''		







"""30
					Logic of searches:
						refrings1 (even point-group symmetry AND theta = 90.), and refrings2 ( odd point-group symmetry AND any theta (including theta=90), OR even point-group symmetry AND theta <> 90.)
						an == -1 exhaustive search
						psi_max - how far rotation in plane can can deviate from 90 or 270 degrees
						psi_max - is used in both exhaustive and local searches
						:::
						>even point-group symmetry AND theta = 90. 
							exhaustive   proj_ali_helical_90   psi_max
											Util.multiref_polar_ali_helical_90  psi_max
											Crosrng_sm_psi  (in util_sparx.cpp)  psi_max - WILL SEARCH for BOTH PSI=0 AND 180 NO MATTER WHAT
																					flag - no mirror or mirror, DOES NOT CHECK MIRRORED
											
							local       reference projections phi= [0,180,delta], theta=90
										proj_ali_helical_90_local   an[N_step], psi_max,  (sp_alignment.py)
											Util.multiref_polar_ali_helical_90_local   psi_max
											Uses the following construct
											if ((psi-90.0f) < 90.0f) retvals = Crosrng_sm_psi(crefim[iref], cimage, numr,   0, 0, psi_max);
											else                     retvals = Crosrng_sm_psi(crefim[iref], cimage, numr, 180, 0, psi_max);
										where Crosrng_sm_psi will do search for psi around one angle and mirror is NOT checked.

						>odd point-group symmetry AND any theta (including theta=90), OR even point-group symmetry AND theta <> 90.
							exhaustive   proj_ali_helical           psi_max
											Util.multiref_polar_ali_helical  psi_max
												CALLS Crosrng_pi twice, for 0 and for 180, thus duplicates the work!!
											Crosrng_psi  checks MIRROR
											
							local        proj_ali_helical_local   psi_max
											Util.multiref_polar_ali_helical_local  psi_max
												Uses the following construct
												if (mirror_only == true) {
													if ((psi-90) < 90) retvals = Crosrng_sm_psi(crefim[iref], cimage, numr,   0, 1, psi_max);
													else               retvals = Crosrng_sm_psi(crefim[iref], cimage, numr, 180, 1, psi_max); 
												} else {
													if ((psi-90) < 90) retvals = Crosrng_sm_psi(crefim[iref], cimage, numr,   0, 0, psi_max);
													else               retvals = Crosrng_sm_psi(crefim[iref], cimage, numr, 180, 0, psi_max);
												}
											where Crosrng_sm_psi will do search for psi around one angle and do mirror or not.

				"""






























































































































































































































































































































"""31
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
			"""









































































































'''32
def gchelix_MPI(stack, ref_vol, outdir, maskfile, ir, ou, rs, xr, ynumber,\
	txs, delta, initial_theta, delta_theta, an, maxit, CTF, snr, dp, ndp, dp_step, dphi, ndphi, dphi_step, psi_max,\
	rmin, rmax, fract, nise, npad, sym, user_func_name, datasym,\
	pixel_size, debug, y_restrict, WRAP):

	from sp_alignment      import Numrinit, prepare_refrings, proj_ali_helical, proj_ali_helical_90, proj_ali_helical_local, proj_ali_helical_90_local, helios,helios7
	from sp_utilities      import model_circle, get_image, drop_image, get_input_from_string, pad, model_blank
	from sp_utilities      import bcast_list_to_all, bcast_number_to_all, reduce_EMData_to_root, bcast_EMData_to_all
	from sp_utilities      import send_attr_dict
	from sp_utilities      import get_params_proj, set_params_proj, file_type
	from sp_fundamentals   import rot_avg_image
	from sp_pixel_error    import max_3D_pixel_error
	import os
	import types
	from sp_utilities      import print_begin_msg, print_end_msg, print_msg
	from mpi            import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from mpi            import mpi_recv,  mpi_send
	from mpi            import mpi_reduce, MPI_INT, MPI_SUM
	from sp_filter         import filt_ctf
	from sp_projection     import prep_vol, prgs
	from sp_statistics     import hist_list, varf3d_MPI
	from sp_applications   import MPI_start_end
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
		import sp_global_def
		sp_global_def.LOGFILE =  os.path.join(outdir, sp_global_def.LOGFILE)
	mpi_barrier(MPI_COMM_WORLD)


	nlprms = (2*ndp+1)*(2*ndphi+1)
	if nlprms< number_of_proc:
		ERROR('number of CPUs is larger than the number of helical search, please reduce it or at this moment modify ndp,dphi in the program', "ihrsr_MPI", 1,myid)


	if debug:
		from time import sleep
		while not os.path.exists(outdir):
			print  "Node ",myid,"  waiting..."
			sleep(5)

		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None

	symref = "s"+sym

	ref_a= "P"
	symmetryLower = sym.lower()
	symmetry_string = split(symmetryLower)[0]

	xrng        = get_input_from_string(xr)
	y_restrict       = get_input_from_string(y_restrict)
	ynumber	    = get_input_from_string(ynumber)
	for i in xrange(len(ynumber)):
		if(ynumber[i]%2==1): ynumber[i]=ynumber[i]+1
	yrng =[]

	for i in xrange(len(xrng)): yrng.append(dp/2)
	
	stepx        = get_input_from_string(txs)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(stepx), len(delta))
	if an == "-1": an = [-1] * lstp
	else:          an = get_input_from_string(an)

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
		import sp_user_functions
		user_func = sp_user_functions.factory[user_func_name]

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
		print_msg("number of times to impose initial symmetry: %i\n"%(nise))
		print_msg("npad                                      : %i\n"%(npad))
		print_msg("User function                             : %s\n"%(user_func_name))

	if maskfile:
		if type(maskfile) is types.StringType: mask3D = get_image(maskfile)
		else:                                  mask3D = maskfile
	else: mask3D = None
	#else: mask3D = model_circle(last_ring, nx, nx, nx)

	numr	= Numrinit(first_ring, last_ring, rstep, "F")

	if CTF:
		from sp_reconstruction import recons3d_4nn_ctf_MPI
		from sp_filter         import filt_ctf
	else:	 from sp_reconstruction import recons3d_4nn_MPI

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
	list_of_particles = bcast_list_to_all(list_of_particles, myid,  source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	if debug:
		finfo.write("image_start, image_end: %d %d\n" %(image_start, image_end))
		finfo.flush()
	
	mask2D = pad( model_blank( int(nmax-20),nmax,1,bckg=1.0), nmax, nmax, 1,0.0)

	data = EMData.read_images(stack, list_of_particles)
	#if fourvar:  original_data = []
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		sttt = Util.infomask(data[im], mask2D, False)
		data[im] = data[im] - sttt[0]
		#if fourvar: original_data.append(data[im].copy())
		if CTF:
			st = data[im].get_attr_default("ctf_applied", 0)
			if(st == 0):
				ctf_params = data[im].get_attr("ctf")
				data[im] = filt_ctf(data[im], ctf_params)
				data[im].set_attr('ctf_applied', 1)
	del mask2D

	if debug:
		finfo.write( '%d loaded  \n' % nima )
		finfo.flush()
	
	for i in xrange(len(xrng)): yrng[i]=dp/(2*pixel_size)
	from math import sin, pi
	if ( ou > ( nmax/2.0)*sin( initial_theta*pi/180) - dp/2.0/pixel_size -1.0 ):
		ERROR('ou should be less than or equal to ----( nmax/2.0)*sin( initial_theta*pi/180) - dp/2.0/pixel_size -1.0 ', "ihrsr_MPI", 1,myid)

	if myid == main_node:
		print_msg("Pixel size in Angstroms                   : %5.4f\n\n"%(pixel_size))
		print_msg("Y search range (pix) initialized as       : %s\n\n"%(yrng))

	from time import time	

	#  this is needed for gathering of pixel errors
	disps = []
	recvcount = []
	for im in xrange(number_of_proc):
		if( im == main_node ):  disps.append(0)
		else:                   disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_nima, number_of_proc, im)
		recvcount.append( ie - ib )

	total_iter = 0
	# do the projection matching
	for N_step in xrange(lstp):
		terminate = 0
		Iter = 0
 		while(Iter < max_iter and terminate == 0):
			yrng[N_step]=float(dp)/(2*pixel_size) #will change it later according to dp
			if(ynumber[N_step]==0): stepy = 0.0
			else:                   stepy = (2*yrng[N_step]/ynumber[N_step])

			pixer  = [0.0]*nima
			modphi = [0.0]*nima
			Iter += 1
			total_iter += 1
			if myid == main_node:
				start_time = time()
				print_msg("\nITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.4f, xrange (Pixels) = %5.4f,stepx (Pixels) = %5.4f, yrng (Pixels) = %5.4f,  stepy (Pixels) = %5.4f, y_restrict (Pixels)=%5.4f, ynumber = %3d\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],stepx[N_step],yrng[N_step],stepy,y_restrict[N_step], ynumber[N_step]))
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
			#split refrings to two list: refrings1 (even point-group symmetry AND theta = 90. ), 
			#   or refrings2 ( odd point-group symmetry AND any theta (including theta=90), OR even point-group symmetry AND theta <> 90.
			refrings1= []
			refrings2= []
			sn = int(symmetry_string[1:])
			for i in xrange( len(refrings) ):
				if( sn%2 ==0 and abs( refrings[i].get_attr('n3') ) <1.0e-6 ):
					#  even point-group symmetry AND theta = 90. 
					refrings1.append( refrings[i] )
				else:
					# odd point-group symmetry AND any theta (including theta=90), OR even point-group symmetry AND theta <> 90.
					refrings2.append( refrings[i] )
					"""
					if myid == main_node:
						print_msg("\nphi = %5.2f, theta = %5.2f, psi=%5.2f\n"%( refrings[i].get_attr('phi'), refrings[i].get_attr('theta'), refrings[i].get_attr('psi') ) )
					if myid == main_node:
						print_msg("\nlen(ref1) = %4d, len(ref2) = %4d\n"%(len(refrings1), len(refrings2)) )
					"""		
			del refrings
			from numpy import float32
			dpp = float32(float(dp)/pixel_size)
			dpp = float( dpp )
			dpp_half = dpp/2.0

			for im in xrange( nima ):
				"""
					Logic of searches:
						refrings1 (even point-group symmetry AND theta = 90.), and refrings2 ( odd point-group symmetry AND any theta (including theta=90), OR even point-group symmetry AND theta <> 90.)
						an == -1 exhaustive search
						psi_max - how far rotation in plane can can deviate from 90 or 270 degrees
						psi_max -s used in both exhaustive and local searches
						:::
						>even point-group symmetry AND theta = 90. 
							exhaustive   proj_ali_helical_90   psi_max
											Util.multiref_polar_ali_helical_90  psi_max
											Crosrng_sm_psi  (in util_sparx.cpp)  psi_max - WILL SEARCH for BOTH PSI=0 AND 180 NO MATTER WHAT					
																					flag - no mirror or mirror, DOES NOT CHECK MIRRORED
											
							local       reference projections phi= [0,180,delta], theta=90 
										proj_ali_helical_90_local   an[N_step], psi_max,  (sp_alignment.py)
											Util.multiref_polar_ali_helical_90_local   psi_max
											Crosrng_sm_psi		(in util_sparx.cpp)   psi_max - WILL SEARCH AROUND BOTH PSI=0 AND 180 NO MATTER WHAT					

						>odd point-group symmetry AND any theta (including theta=90), OR even point-group symmetry AND theta <> 90.
							exhaustive   proj_ali_helical           psi_max
											Util.multiref_polar_ali_helical  psi_max
												CALLS Crosrng_pi twice, for 0 and for 180, thus duplicates the work!!
											Crosrng_psi  checks MIRROR
											
							local        proj_ali_helical_local   psi_max
											Util.multiref_polar_ali_helical_local  psi_max
												Uses the following construct
												if (mirror_only == true) {
													if ((psi-90) < 90) retvals = Crosrng_sm_psi(crefim[iref], cimage, numr,   0, 1, psi_max);
													else               retvals = Crosrng_sm_psi(crefim[iref], cimage, numr, 180, 1, psi_max); 
												} else {
													if ((psi-90) < 90) retvals = Crosrng_sm_psi(crefim[iref], cimage, numr,   0, 0, psi_max);
													else               retvals = Crosrng_sm_psi(crefim[iref], cimage, numr, 180, 0, psi_max);
												}
											where Crosrng_sm_psi will do search for psi around one angle and do mirror or not.

				"""
				peak1 = None
				peak2 = None
				#print im, get_params_proj(data[im])
				if( len(refrings1) > 0):
					if  an[N_step] == -1:
						peak1, phihi1, theta1, psi1, sxi1, syi1, t11 = \
						proj_ali_helical_90(data[im], refrings1, numr, xrng[N_step], yrng[N_step], stepx[N_step], ynumber[N_step], psi_max, finfo)
					else:
						peak1, phihi1, theta1, psi1, sxi1, syi1, t11 = \
						proj_ali_helical_90_local(data[im], refrings1, numr, xrng[N_step], yrng[N_step], stepx[N_step], ynumber[N_step], an[N_step], psi_max, finfo, yrnglocal=y_restrict[N_step])
					#print "  1  ",im, peak1, phihi1, theta1, psi1, sxi1, syi1
				if( len(refrings2) > 0):
					if  an[N_step] == -1:
						peak2, phihi2, theta2, psi2, sxi2, syi2, t12 = \
						proj_ali_helical(data[im], refrings2, numr, xrng[N_step], yrng[N_step], stepx[N_step], ynumber[N_step], psi_max, finfo)
					else:
						peak2, phihi2, theta2, psi2, sxi2, syi2, t12 = \
						proj_ali_helical_local(data[im], refrings2, numr, xrng[N_step], yrng[N_step], stepx[N_step], ynumber[N_step], an[N_step], psi_max, finfo, yrnglocal=y_restrict[N_step])
					#print "  2  ",im, peak2, phihi2, theta2, psi2, sxi2, syi2
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
				#print "  3  ",im, peak, phihi, theta, psi, sxi, syi
				if(peak > -1.0e22):
					if WRAP == 1:
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
					#print "  4  ",im, peak, phihi, theta, psi, sxi, syi

					# unique ranges of azimuthal angle for ortho-axial and non-ortho-axial projection directions are identified by [k0,k1) and [k2,k3), where k0, k1, k2, k3 are floats denoting azimuthal angles.
					# Eulerian angles whose azimuthal angles are mapped into [k2, k3) are related to Eulerian angles whose azimuthal angles are mapped into [k0, k1) by an in-plane mirror operaton along the x-axis.

					tp = Transform({"type":"spider","phi":phihi,"theta":theta,"psi":psi})
					tp.set_trans( Vec2f( -sxi, -syi ) )

					k0 =   0.0
					k2 = 180.0
					if( abs( tp.at(2,2) )<1.0e-6 ):
						if (symmetry_string[0] =="c"):
							if sn%2 == 0:  k1=360.0/sn
							else:          k1=360.0/2/sn
						elif (symmetry_string[0] =="d"):
							if sn%2 == 0:  k1=360.0/2/sn
							else:          k1=360.0/4/sn
					else:
						if (symmetry_string[0] =="c"):  k1=360.0/sn
						if (symmetry_string[0] =="d"):  k1=360.0/2/sn
					k3 = k1 +180.0

					from sp_utilities import get_sym
					T = get_sym(symmetry_string[0:])

					d1tp = tp.get_params('spider')
					sxnew    = -d1tp["tx"]
					synew    = -d1tp["ty"]
					phinew   =  d1tp['phi']
					thetanew =  d1tp["theta"]
					psinew   =  d1tp["psi"]
					del d1tp

					#print "  5  ",im, phinew, thetanew, psinew, sxnew, synew
					#print k0,k1,k2,k3

					for i in xrange( len(T) ):
						ttt = tp*Transform({"type":"spider","phi":T[i][0],"theta":T[i][1],"psi":T[i][2]})
						d1  = ttt.get_params("spider")

						if ( abs( tp.at(2,2) )<1.0e-6 ):
							if( sn%2==1 ): # theta=90 and n odd, only one of the two region match

								if( ( d1['phi'] >= k0 and d1['phi'] < k1 ) or ( d1['phi'] >= k2 and d1['phi'] < k3 )):

									sxnew    = -d1["tx"]
									synew    = -d1["ty"]
									phinew   =  d1['phi']
									thetanew =  d1["theta"]
									psinew   =  d1["psi"]

									# For boundary cases where phihi is exactly on the boundary of the unique range, there may be two symmetry related Eulerian angles which are both in the unique 
									# range but whose psi differ by 180. 
									# For example, (180,90,270) has two symmetry related angles in unique range: (180,90,270) and (180, 90, 90)
									# In local search, psi should stay within neighborhood of original value, so take the symmetry related
									# Eulerian angles in unique range which does not change psi by 180.
									if an[N_step] != -1:
										if abs(psinew - psi) < 90:
											break
							else: #for theta=90 and n even, there is no mirror version during aligment, so only consider region [k0,k1]

								if( d1['phi'] >= float(k0) and d1['phi'] < float(k1)  ) :

									sxnew    = - d1["tx"]
									synew    = - d1["ty"]
									phinew   = d1['phi']
									thetanew = d1["theta"]
									psinew   = d1["psi"]

						else: #theta !=90, # if theta >90, put the projection into [k2,k3]. Otherwise put it into the region [k0,k1]

							if( sn==1):
								sxnew    = sxi
								synew    = syi
								phinew   = phihi
								thetanew = theta
								psinew   = psi
							else:

								if (tp.at(2,2) >0.0): #theta <90
									
									if(  d1['phi'] >= float(k0) and d1['phi'] < float(k1)):
										if( cos( pi*float( d1['theta'] )/180.0 )>0.0 ):
			
											sxnew    = - d1["tx"]
											synew    = - d1["ty"]
											phinew   = d1['phi']
											thetanew = d1["theta"]
											psinew   = d1["psi"]
		
								else:
									if(  d1['phi'] >= float(k2) and d1['phi'] < float(k3)):
										if( cos( pi*float( d1['theta'] )/180.0 )<0.0 ):

											sxnew    = - d1["tx"]
											synew    = - d1["ty"]
											phinew   = d1['phi']
											thetanew = d1["theta"]
											psinew   = d1["psi"]

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
				#print "  6  ",im, phinew, thetanew, psinew, sxnew, synew

				#if(im==2):
				#	from sys import exit
				#	exit()

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
				from sp_utilities import write_text_file
				write_text_file([range(len(recvbuf)), recvbuf], os.path.join(outdir, "pixer_%04d_%04d.txt"%(N_step+1,Iter)) )
				from sp_statistics import hist_list
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
			from mpi import mpi_recv, mpi_send, MPI_COMM_WORLD, MPI_FLOAT
			if myid == main_node:
				
				fexp = open(os.path.join(outdir, "parameters_%04d_%04d.txt"%(N_step+1,Iter)),"w")
				for n in xrange(number_of_proc):
					if n!=main_node:
						import sp_global_def
						t = mpi_recv(recvcount[n]*m,MPI_FLOAT, n, sp_global_def.SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
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
				import sp_global_def
				mpi_send(nvalue, recvcount[myid]*m, MPI_FLOAT, main_node, sp_global_def.SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
				del nvalue
			if myid == main_node:
				print_msg("Time to write parameters = %d\n"%(time()-start_time))
				start_time = time()

			if CTF: vol = recons3d_4nn_ctf_MPI(myid, data, symmetry=sym, snr = snr, npad = npad, xysize = xysize, zsize = zsize)
			else:    vol = recons3d_4nn_MPI(myid, data, symmetry=sym, snr = snr, npad = npad, xysize = xysize, zsize = zsize)

			if myid == main_node:
				print_msg("\n3D reconstruction time = %d\n"%(time()-start_time))
				start_time = time()

			"""
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
			"""
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
						import sp_global_def
						if n!=main_node: mpi_send(lprms[2*recvpara[2*n]:2*recvpara[2*n+1]], 2*(recvpara[2*n+1]-recvpara[2*n]), MPI_FLOAT, n, sp_global_def.SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
						else:    list_dps = lprms[2*recvpara[2*0]:2*recvpara[2*0+1]]
				else:
					import sp_global_def
					list_dps = mpi_recv((para_end-para_start)*2, MPI_FLOAT, main_node, sp_global_def.SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)

				list_dps = map(float, list_dps)

				local_pos = [0.0, 0.0, -1.0e20]
				for i in xrange(para_end-para_start):
					fvalue = helios7(vol, pixel_size, list_dps[i*2], list_dps[i*2+1], fract, rmax, rmin)
					if(fvalue >= local_pos[2]):
						local_pos = [list_dps[i*2], list_dps[i*2+1], fvalue ]
				if myid == main_node:
					list_return = [0.0]*(3*number_of_proc)
					for n in xrange(number_of_proc):
						import sp_global_def
						if n != main_node: list_return[3*n:3*n+3]                 = mpi_recv(3,MPI_FLOAT, n, sp_global_def.SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)
 						else:              list_return[3*main_node:3*main_node+3]  = local_pos[:]
				else:
					import sp_global_def
					mpi_send(local_pos, 3, MPI_FLOAT, main_node, sp_global_def.SPARX_MPI_TAG_UNIVERSAL, MPI_COMM_WORLD)

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

					vol  = vol.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
					print_msg("New delta z and delta phi      : %s,    %s\n\n"%(dp,dphi))		
			
			else:
				if myid==main_node:
					#  in the first nise steps the symmetry is imposed
					vol = vol.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
					print_msg("Imposed delta z and delta phi      : %s,    %s\n\n"%(dp,dphi))
			if(myid==main_node):
				fofo = open(os.path.join(outdir,datasym),'a')
				fofo.write('  %12.4f   %12.4f\n'%(dp,dphi))
				fofo.close()
				ref_data = [vol, mask3D]
				#if  fourvar:  ref_data.append(varf)
				vol = user_func(ref_data)
				vol = vol.helicise(pixel_size, dp, dphi, fract, rmax, rmin)

				drop_image(vol, os.path.join(outdir, "volf%04d.hdf"%(total_iter)))
				print_msg("\nSymmetry search and user function time = %d\n"%(time()-start_time))
				start_time = time()

			bcast_EMData_to_all(vol, myid, main_node)
			dp   = bcast_number_to_all(dp,   source_node = main_node)
			dphi = bcast_number_to_all(dphi, source_node = main_node)
			# del varf
	par_str = ["xform.projection"]
	if myid == main_node:
	   	if(file_type(stack) == "bdb"):
	        	from sp_utilities import recv_attr_dict_bdb
	        	recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        else:
	        	from sp_utilities import recv_attr_dict
	        	recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
		print_msg("Time to write header information= %d\n"%(time()-start_time))
		start_time = time()
	else:	       send_attr_dict(main_node, data, par_str, image_start, image_end)
	if myid == main_node: print_end_msg("ihrsr_MPI")
'''


















































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































"""33

	if(myid == 0):
		pid_list = read_text_file("main000/chunk0.txt")
		nima = len(pid_list)
	else: nima = 0
	nima = bcast_number_to_all(nima, source_node = 0)
	if(myid != 0):
		pid = [-1]*nima
	pid_list = mpi_bcast(pid_list, nima, MPI_INT, 0, MPI_COMM_WORLD)
	pid_list = map(int, pid_list)
	
	image_start, image_end = MPI_start_end(nima, nproc, myid)
	prjlist = [EMData.read_images(prj_stack, pid_list[image_start:image_end])]
	"""
"""34
	if(myid == 0):
		pid_list = read_text_file("main000/chunk1.txt")
		nima = len(pid_list)
	else: nima = 0
	nima = bcast_number_to_all(nima, source_node = 0)
	if(myid != 0):
		pid = [-1]*nima
	"""











"""35
	from sp_fundamentals import fdecimate
	from sp_utilities import get_params_proj,set_params_proj
	scale = 384./54.
	for i in xrange(len(prjlist)):
		prjlist[k][i] = fdecimate(prjlist[i],54,54)
		ctf_params = prjlist[i].get_attr("ctf")
		ctf_params.apix *= scale
		prjlist[i].set_attr('ctf', ctf_params)
		phi,theta,psi,sx,sy = get_params_proj(prjlist[i])
		set_params_proj(prjlist[i],[phi,theta,psi,sx/scale,sy/scale])
	"""





"""36


	bckgnoise = get_im("bckgnoise.hdf")
	nlx = bckgnoise.get_xsize()
	datastamp = read_text_file("defgroup_stamp.txt")
	#nnnx = 200#prjlist[0].get_ysize()
	from sp_utilities import get_params_proj,set_params_proj
	for i in xrange(len(prjlist)):
		#phi,theta,psi,sxs,sys = get_params_proj(prjlist[i])
		try:
			stmp = prjlist[i].get_attr("ptcl_source_image")
		except:
			try:
				stmp = prjlist[i].get_attr("ctf")
				stmp = round(stmp.defocus,4)
			except:
				ERROR("Either ptcl_source_image or ctf has to be present in the header.","get_shrink_data",1, myid)
		try:
			indx = datastamp.index(stmp)
		except:
			ERROR("Problem with indexing ptcl_source_image.","get_shrink_data",1, myid)


		prjlist[i] = fft(prjlist[i])
		#prjlist[i] = fshift(prjlist[i],sxs,sys)
		prjlist[i].set_attr("padffted",1)
		prjlist[i].set_attr("npad",1)
		#set_params_proj(prjlist[i] ,[phi,theta,psi,0.0,0.0])
		prjlist[i].set_attr("bckgnoise", [bckgnoise[k,indx] for k in xrange(nlx)]) #  This is index on bckgnoise table



	"""









































"""37
	if(myid == 0):
		if(listfile):
			from sp_utilities import read_text_file
			pid_list = read_text_file(listfile, 0)
			pid_list = map(int, pid_list)
		elif(group > -1):
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
	"""























































































































































































































































'''38
	qt = 0.0
	for i in xrange(len(ssnr1)):
		tqt = ssnr1[i][1] - ssnr2[i][1]
		if( tqt<qt ): qt = tqt
	for i in xrange(len(ssnr1)): ssnr1[i][1] -= (ssnr2[i][1] + qt)
	from sp_utilities import dropSpiderDoc19289
	dropSpiderDoc(ssnr_text_file+".doc", ssnr1)
	dropImage(vol_ssnr2, output_volume+"2.spi", "s")
	'''










































































































"""39
		qt = 0.0
		for i in xrange(len(ssnr2)):
			tqt = ssnr1[i][1] - ssnr2[i][1]
			if( tqt<qt ): qt = tqt
		for i in xrange(len(ssnr1)): ssnr1[i][1] -= (ssnr2[i][1] + qt)
		
		dropSpiderDoc(ssnr_text_file+".doc", ssnr1)
		vol_ssnr2, output_volume+"2.spi", "s")
		"""



























































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































"""40
def incvar(prefix, nfile, nprj, output, fl, fh, radccc, writelp, writestack):
	from sp_statistics import variancer, ccc
	from string     import atoi, replace, split, atof
	from sp_utilities  import get_im, circumference, model_circle, drop_image
	from sp_filter     import filt_btwl
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




































































































"""41
	Util.mul_scalar(avg1, 1.0/float(total_img//2+total_img%2 - 1 ))
	avg1.write_image(avgfileE)
	Util.mul_scalar(avg2, 1.0/float(total_img//2 - 1) )
	avg2.write_image(avgfileO)
	"""



























"""42
	Util.mul_scalar(var1, 1.0/float(total_img//2+total_img%2 - 1 ))
	var1.write_image(varfileE)
	Util.mul_scalar(var2, 1.0/float(total_img//2 - 1) )
	var2.write_image(varfileO)
	"""





















"""43
	  writelp means overwrite original stacks with repaired ones
	  writestack means write new stacks with low-pass filtered data
	"""


































































































































































































































































































































































































































































































































































































































































































































































'''44
#  06-12-2014 code lifted
def within_group_refinement(data, maskfile, randomize, ir, ou, rs, xrng, yrng, step, dst, maxit, FH, FF):

	# Comment by Zhengfan Yang 03/11/11
	# This is a simple version of ali2d_data (down to the bone), no output dir, no logfile, no CTF, no MPI or CUDA, no Fourvar, no auto stop, no user function
	
	from sp_alignment    import Numrinit, ringwe, ali2d_single_iter
	from sp_filter	  import filt_tanl
	from sp_fundamentals import fshift
	from random	  import randint, random
	from sp_statistics   import ave_series
	from sp_utilities    import get_input_from_string, model_circle, set_params2D, get_params2D, combine_params2, inverse_transform2

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
			if( FH > 0.0):
				fl = 0.1+(FH-0.1)*Iter/float(max_iter-1)
				tavg = filt_tanl(tavg, fl, FF)
			if total_iter == len(xrng)*max_iter:  return tavg
			#if( xrng[0] > 0.0 ): cs[0] = sx_sum/float(nima)
			#if( yrng[0] > 0.0 ): cs[1] = sy_sum/float(nima)
			#tavg = fshift(tavg, -cs[0], -cs[1])
			if Iter%4 != 0 or total_iter > max_iter*len(xrng)-10: delta = 0.0
			else: delta = dst
			sx_sum, sy_sum = ali2d_single_iter(data, numr, wr, cs, tavg, cnx, cny, xrng[N_step], yrng[N_step], step[N_step], mode=mode, CTF=False, delta=delta)

'''
'''45
def Xwithin_group_refinement(data, maskfile, randomize, ir, ou, rs, xrng, yrng, step, dst, maxit, FH, FF, method = ""):

	# Comment by Zhengfan Yang 03/11/11
	# This is a simple version of ali2d_data (down to the bone), no output dir, no logfile, no CTF, no MPI or CUDA, no Fourvar, no auto stop, no user function

	from sp_alignment    import Numrinit, ringwe, ali2d_single_iter
	from sp_filter	      import filt_tanl
	from sp_fundamentals import fshift, fft
	from random	      import randint, random
	from sp_statistics   import ave_series
	from sp_utilities    import get_input_from_string, model_circle, center_2D
	from sp_utilities    import set_params2D, get_params2D, combine_params2, inverse_transform2

	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit);
	nima = len(data)
	nx = data[0].get_xsize()
	if last_ring == -1:  last_ring = nx/2-2
	if maskfile: mask = maskfile
	else:        mask = model_circle(last_ring, nx, nx)

	lx = [0]*nima
	ly = [0]*nima
	for im in range(nima):
		alpha, sx, sy, mirrorn, dummy = get_params2D(data[im])
		alphai, sxi, syi, dummy    = combine_params2(0.0, sx, sy, 0, -alpha, 0.,0.,0)
		lx[im] = int(round(sxi,0))
		ly[im] = int(round(syi,0))
		Util.cyclicshift(data[im] , {"dx":lx[im],"dy":ly[im]})
		sxi -= lx[im]
		syi -= ly[im]
		if randomize :
			alphan, sxn, syn, mirrorn    = combine_params2(0.0, sxi, syi, 0, random()*360.0, 0.0, 0.0, randint(0, 1))
			#alphan, sxn, syn, mirrorn = combine_params2(0.0, randint(-xrng[0],xrng[0]), randint(-xrng[0],xrng[0]), 0, random()*360.0, 0, 0, randint(0, 1))
		else:
			alphan, sxn, syn, dummy = combine_params2(0.0, sxi, syi, 0, alpha, 0.0, 0.0, 0)
		set_params2D(data[im], [alpha, sxn, syn, mirrorn, 1.0])
		

	cnx = nx/2+1
	cny = cnx
	mode = "F"
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr = ringwe(numr, mode)

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
		tavg = ave_series(data)
		for N_step in range(len(xrng)):
			nope = 0
			Iter = 0
			while(nope < len(data)//1 and Iter < max_iter ):
				total_iter += 1
				Iter += 1
				if( FH > 0.0):
					tavg = filt_tanl(fft(tavg), FH, FF)
					if( xrng[0] > 0.0 ): cs[0] = sx_sum/float(nima)
					if( yrng[0] > 0.0 ): cs[1] = sy_sum/float(nima)
					tavg = fft(fshift(tavg, -cs[0], -cs[1]))
				else:
					tavg = filt_tanl(tavg, FH, FF)
				sx_sum, sy_sum, nope = ali2d_single_iter(data, numr, wr, cs, tavg, cnx, cny, \
															xrng[N_step], yrng[N_step], step[N_step], \
															mode=mode, CTF=False, random_method = method)
				#print  "  iteration  shc   %03d   %03d   %7.2f    %7.2f  "%(total_iter,nope,cs[0],cs[1])
				#print total_iter,nope
				#for i in data:  print "  ",i.get_attr('previousmax'),
				#print "  "
				#tavg.write_image('tata.hdf',total_iter-1)
				tavg = ave_series(data)
		"""
		tavg.write_image('tata.hdf')
		for Iter in xrange(0):#max_iter):  # large number
			total_iter += 1
			tavg = ave_series(data)
			if( FH > 0.0):
				fl = FH
				tavg = filt_tanl(fft(tavg), fl, FF)
				if total_iter == len(xrng)*max_iter:  return fft(tavg)
				if( xrng[0] > 0.0 ): cs[0] = sx_sum/float(nima)
				if( yrng[0] > 0.0 ): cs[1] = sy_sum/float(nima)
				tavg = fft(fshift(tavg, -cs[0], -cs[1]))
			else:
				if total_iter == len(xrng)*max_iter:  return tavg
				if( xrng[0] > 0.0 ): cs[0] = sx_sum/float(nima)
				if( yrng[0] > 0.0 ): cs[1] = sy_sum/float(nima)
				tavg = fshift(tavg, -cs[0], -cs[1])
			
			print  "  iteration  ***   %03d   %7.2f    %7.2f  "%(total_iter,cs[0],cs[1])
			if Iter%4 != 0 or total_iter > max_iter*len(xrng)-10: delta = 0.0
			else:                                                 delta = dst
			sx_sum, sy_sum, nope = ali2d_single_iter(data, numr, wr, cs, tavg, cnx, cny, \
														xrng[N_step], yrng[N_step], step[N_step], \
														mode=mode, CTF=False, delta=delta)
			if( (abs(cs[0]) + abs(cs[1])) < 0.01 and Iter > 1):  break
		"""


	elif( method == "PCP"):
		from sp_isac import prepref
		from sp_utilities import model_circle
		stp = step[-1]
		rings = prepref(data, model_circle(nx//2-1,nx,nx), cnx, cnx, numr, mode, xrng[0], xrng[0], stp)
		sxprint(" rings  ",len(rings))
		for im in range(len(data)):
			rings[im][0][0].set_attr("sxi",0)
			rings[im][0][0].set_attr("syi",0)
			rings[im][0][0].set_attr("inx",nx)
		tavg = ave_series(data)
		for N_step in range(len(xrng)):
			sxprint(" xrng ",xrng[N_step])
			for Iter in range(max_iter):
				total_iter += 1
				if( FH > 0.0):
					fl = 0.1+(FH-0.1)*Iter/float(max_iter-1)
					tavg = filt_tanl(tavg, fl, FF)
					"""
					tavg = filt_tanl(fft(tavg), fl, FF)
					if total_iter == len(xrng)*max_iter:  return fft(tavg)
					if( xrng[0] > 0.0 ): cs[0] = sx_sum/float(nima)
					if( yrng[0] > 0.0 ): cs[1] = sy_sum/float(nima)
					tavg = fft(fshift(tavg, -cs[0], -cs[1]))
					"""
				else:
					"""
					if total_iter == len(xrng)*max_iter:  return tavg
					if( xrng[0] > 0.0 ): cs[0] = sx_sum/float(nima)
					if( yrng[0] > 0.0 ): cs[1] = sy_sum/float(nima)
					tavg = fshift(tavg, -cs[0], -cs[1])
					"""
				cs = [0,0]
				#print  "  iteration  std   %03d   %7.2f    %7.2f  "%(total_iter,cs[0],cs[1])
				if Iter%4 != 0 or total_iter > max_iter*len(xrng)-10: delta = 0.0
				else:                                                 delta = dst
				sx_sum, sy_sum, nope = ali2d_single_iter(rings, numr, wr, cs, tavg, cnx, cny, \
															xrng[N_step], yrng[N_step], step[N_step], \
															mode=mode, CTF=False, delta=delta, random_method = method)
				for im in range(len(data)):
					alpha, tx, ty, mir, scale = get_params2D(rings[im][0][0])
					set_params2D(data[im],[alpha, tx, ty, mir, scale])
				tavg = ave_series(data)
				#print  "tata ",total_iter-1,Util.infomask(tavg,None,True)
				#tavg.write_image('tata.hdf',total_iter-1)
	else:
		tavg = ave_series(data)
		for N_step in range(len(xrng)):
			for Iter in range(max_iter):
				total_iter += 1
				cs = Util.infomask(tavg, mask, False)
				tavg -= cs[0]
				if( FH > 0.0):
					fl = 0.1+(FH-0.1)*Iter/float(max_iter-1)
					tavg = filt_tanl(tavg, fl, FF)
				"""
					tavg = filt_tanl(fft(tavg), fl, FF)
					if total_iter == len(xrng)*max_iter:  return fft(tavg)
					if( xrng[0] > 0.0 ): cs[0] = sx_sum/float(nima)
					if( yrng[0] > 0.0 ): cs[1] = sy_sum/float(nima)
					tavg = fft(fshift(tavg, -cs[0], -cs[1]))
				else:
					if total_iter == len(xrng)*max_iter:  return tavg
					if( xrng[0] > 0.0 ): cs[0] = sx_sum/float(nima)
					if( yrng[0] > 0.0 ): cs[1] = sy_sum/float(nima)
					tavg = fshift(tavg, -cs[0], -cs[1])
				"""
				cs = [0,0]
				if Iter%4 != 0 or total_iter > max_iter*len(xrng)-10:
					delta = 0.0
				else:
					delta = dst
					asx = 0
					asy = 0
					#tavg, asx,asy = \
					#	center_2D(tavg, center_method = 7, searching_range = cnx//2, self_defined_reference = mask)
					if(asx != 0 or asy != 0):
						#  Shift images by this additional amount
						for im in range(nima):
							alpha, sx, sy, mir, scale = get_params2D(data[im])
							if mir == 0:  sxn = sx-asx
							else:  sxn = sx+asx
							syn = sy-asy
							alphai, sxn, syn, dummy  = combine_params2(0, sxn, syn, 0, -alpha, 0,0, 0)
							sxn += asx
							syn += asy
							Util.cyclicshift(data[im] , {"dx":-asx,"dy":-asy})
							lx[im] += asx
							ly[im] += asy
							alphai, sxn, syn, dummy  = combine_params2(0, sxn, syn, 0, alpha, 0,0, 0)
							set_params2D(data[im], [alpha, sxn, syn, mir, 1.0])

				sx_sum, sy_sum, nope = ali2d_single_iter(data, numr, wr, cs, tavg, cnx, cny, \
															xrng[N_step], yrng[N_step], step[N_step], \
															mode=mode, CTF=False, delta=delta)

				tavg = ave_series(data)
				#for im in data:  print get_params2D(im)
				#print  "tata ",total_iter-1,Util.infomask(tavg,None,True)
				#tavg.write_image('tata.hdf',total_iter-1)

		#  Shift data back and adjust parameters
		for im in range(nima):
			alpha, sx, sy, mir, scale = get_params2D(data[im])
			alphai, sxn, syn, dummy  = combine_params2(0, sx, sy, 0, -alpha, 0,0, 0)
			Util.cyclicshift(data[im] , {"dx":-lx[im],"dy":-ly[im]})
			alphai, sxn, syn, dummy  = combine_params2(0, sxn-lx[im], syn-ly[im], 0, alpha, 0,0, 0)
			set_params2D(data[im], [alpha, sxn, syn, mir, 1.0])

	return tavg

'''









































































"""46
		tavg.write_image('tata.hdf')
		for Iter in xrange(0):#max_iter):  # large number
			total_iter += 1
			tavg = ave_series(data)
			if( FH > 0.0):
				fl = FH
				tavg = filt_tanl(fft(tavg), fl, FF)
				if total_iter == len(xrng)*max_iter:  return fft(tavg)
				if( xrng[0] > 0.0 ): cs[0] = sx_sum/float(nima)
				if( yrng[0] > 0.0 ): cs[1] = sy_sum/float(nima)
				tavg = fft(fshift(tavg, -cs[0], -cs[1]))
			else:
				if total_iter == len(xrng)*max_iter:  return tavg
				if( xrng[0] > 0.0 ): cs[0] = sx_sum/float(nima)
				if( yrng[0] > 0.0 ): cs[1] = sy_sum/float(nima)
				tavg = fshift(tavg, -cs[0], -cs[1])
			
			print  "  iteration  ***   %03d   %7.2f    %7.2f  "%(total_iter,cs[0],cs[1])
			if Iter%4 != 0 or total_iter > max_iter*len(xrng)-10: delta = 0.0
			else:                                                 delta = dst
			sx_sum, sy_sum, nope = ali2d_single_iter(data, numr, wr, cs, tavg, cnx, cny, \
														xrng[N_step], yrng[N_step], step[N_step], \
														mode=mode, CTF=False, delta=delta)
			if( (abs(cs[0]) + abs(cs[1])) < 0.01 and Iter > 1):  break
		"""


















"""47
					tavg = filt_tanl(fft(tavg), fl, FF)
					if total_iter == len(xrng)*max_iter:  return fft(tavg)
					if( xrng[0] > 0.0 ): cs[0] = sx_sum/float(nima)
					if( yrng[0] > 0.0 ): cs[1] = sy_sum/float(nima)
					tavg = fft(fshift(tavg, -cs[0], -cs[1]))
					"""

"""48
					if total_iter == len(xrng)*max_iter:  return tavg
					if( xrng[0] > 0.0 ): cs[0] = sx_sum/float(nima)
					if( yrng[0] > 0.0 ): cs[1] = sy_sum/float(nima)
					tavg = fshift(tavg, -cs[0], -cs[1])
					"""




































"""49
					tavg = filt_tanl(fft(tavg), fl, FF)
					if total_iter == len(xrng)*max_iter:  return fft(tavg)
					if( xrng[0] > 0.0 ): cs[0] = sx_sum/float(nima)
					if( yrng[0] > 0.0 ): cs[1] = sy_sum/float(nima)
					tavg = fft(fshift(tavg, -cs[0], -cs[1]))
				else:
					if total_iter == len(xrng)*max_iter:  return tavg
					if( xrng[0] > 0.0 ): cs[0] = sx_sum/float(nima)
					if( yrng[0] > 0.0 ): cs[1] = sy_sum/float(nima)
					tavg = fshift(tavg, -cs[0], -cs[1])
				"""

















'''50


#  commented out to prevent problems 03/02/2015
def within_group_refinement_fast(data, dimage, maskfile, randomize, ir, ou, rs, xrng, yrng, step, maxrange, dst, maxit, FH, FF):
	#  It is not used anywhere, however, the check of boundaries has to be added or the code removed
	# Comment by Zhengfan Yang 03/11/11
	# This is a simple version of ali2d_data (down to the bone), no output dir, no logfile, no CTF, no MPI or CUDA, no Fourvar, no auto stop, no user function
	
	from sp_alignment    import Numrinit, ringwe, ali2d_single_iter_fast
	from sp_filter	      import filt_tanl
	from sp_fundamentals import fshift, rot_shift2D, cyclic_shift
	from random	      import randint, random
	from math         import cos, sin, radians
	from sp_statistics   import ave_series
	from sp_utilities    import get_input_from_string, model_circle, model_blank, set_params2D, get_params2D, combine_params2, inverse_transform2

	first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit);
	nima = len(data)
	nx = data[0].get_xsize()
	if last_ring == -1:  last_ring = nx/2-2
	if maskfile: mask = maskfile
	else: mask = model_circle(last_ring, nx, nx)

	params = [[0.,0.,0.,0] for im in xrange(nima) ]
	if randomize:
		for im in xrange(nima):
			#alpha, sx, sy, mirror, scale = get_params2D(data[im])
			#alphai, sxi, syi, mirrori = inverse_transform2(alpha, sx, sy)
			#alphan, sxn, syn, mirrorn = combine_params2(0.0, -sxi, -syi, 0, random()*360.0, 0,0, randint(0, 1))
			params[im] = [random()*360.0, 0, 0, randint(0, 1)]
	else:
		for im in xrange(nima):
			alpha, sx, sy, mirror, scale = get_params2D(data[im])
			params[im] = [alpha, 0, 0, mirror]

	tavg = model_blank(nx,nx)
	for im in xrange(nima):
		Util.add_img( tavg, rot_shift2D(data[im], params[im][0], params[im][1], params[im][2], params[im][3]) )
	tavg /= nima
	tavg = filt_tanl(tavg, 0.1, FF)
	#tavg = EMData('image.hdf')
	#for im in xrange(nima):  print  im,params[im]

	cnx = nx/2+1
	cny = cnx
	mode = "F"
	numr = Numrinit(first_ring, last_ring, rstep, mode)
	wr = ringwe(numr, mode)

	cs = [0.0]*2
	total_iter = 0
	for Iter in xrange(max_iter):
		for N_step in xrange(len(xrng)):
			total_iter += 1
			if Iter%4 != 0 or total_iter > max_iter*len(xrng)-10: delta = 0.0
			else:                                                 delta = dst
			#delta=0.0
			#for im in xrange(nima):		print  " sx, sy:  %2d   %2d"%(params[im][1], params[im][2]) 
			ali2d_single_iter_fast(data, dimage, params, numr, wr, cs, tavg, cnx, cny, \
					xrng[N_step], yrng[N_step], step[N_step], maxrange = maxrange, mode=mode, delta=delta)
			tavg = model_blank(nx,nx)
			sx_sum = 0.0
			sy_sum = 0.0
			for im in xrange(nima):
				sxi = -params[im][1]
				syi = -params[im][2]
				"""
				co =  cos(radians(alphai))
				so = -sin(radians(alphai))
				"""
				co =  cos(radians(params[im][0]))
				so = -sin(radians(params[im][0]))
				sxs = sxi*co - syi*so
				sys = sxi*so + syi*co
				Util.add_img( tavg, rot_shift2D(data[im], params[im][0], sxs, sys, params[im][3]) )

				if params[im][3] == 0: sx_sum += sxs
				else:                  sx_sum -= sxs
				sy_sum += sys

			tavg /= nima
			if( FH > 0.0):
				fl = 0.1+(FH-0.1)*Iter/float(max_iter-1)
			tavg = filt_tanl(tavg, fl, FF)
			"""
			if( xrng[0] > 0.0 ): sx_sum = int(sx_sum/float(nima)+0.5)
			if( yrng[0] > 0.0 ): sy_sum = int(sy_sum/float(nima)+0.5)
			#print ' ave shift',sx_sum, sy_sum
			tavg = cyclic_shift(tavg, -sx_sum, -sy_sum)
			"""
			if( xrng[0] > 0.0 ): sx_sum = sx_sum/float(nima)
			if( yrng[0] > 0.0 ): sy_sum = sy_sum/float(nima)
			tavg = fshift(tavg, -sx_sum, -sy_sum)
			#tavg.write_image('tavg.hdf',total_iter-1)
	for im in xrange(nima):
		sxi = -params[im][1]
		syi = -params[im][2]
		"""
		co =  cos(radians(alphai))
		so = -sin(radians(alphai))
		"""
		co =  cos(radians(params[im][0]))
		so = -sin(radians(params[im][0]))
		sxs = sxi*co - syi*so
		sys = sxi*so + syi*co
		params[im][1] = sxs
		params[im][2] = sys

	return tavg, params
'''














































































































































"""51
					sxprint("GBT: %6d    %12.1f    %12.1f   %6.2f   %6.2f   %6.2f"%(i,previous,mstack[i][6],mstack[i][3][0],mstack[i][3][1],mstack[i][3][2]))
				else:
					sxprint("GBC: %6d    %12.1f    %12.1f   %6.2f   %6.2f   %6.2f"%(i,previous,mstack[i][6],mstack[i][3][0],mstack[i][3][1],mstack[i][3][2]))
			else: sxprint("NBT: %6d    %12.1f    %12.1f"%(i,previous,mstack[i][6]))
			"""












































































































































































































































































































"""52
	ndat = 0
	for i in xrange(0, nfils, 2): ndat += len(filaments[i])
	td = [None]*ndat
	k = 0
	for i in xrange(0, nfils, 2):
		for j in xrange(indcs[i][0],indcs[i][1]):
			td[k] = data[j]
			k += 1

	if CTF:  vol_even = recons3d_4nn_ctf_MPI(myid, td, symmetry=sym, snr = snr, npad = npad, xysize = xysize, zsize = zsize)
	else:    vol_even = recons3d_4nn_MPI(myid, td, symmetry=sym, snr = snr, npad = npad, xysize = xysize, zsize = zsize)

	ndat = 0
	for i in xrange(1, nfils, 2): ndat += len(filaments[i])
	td = [None]*ndat
	k = 0
	for i in xrange(1, nfils, 2):
		for j in xrange(indcs[i][0],indcs[i][1]):
			td[k] = data[j]
			k += 1

	if CTF:  vol_odd = recons3d_4nn_ctf_MPI(myid, td, symmetry=sym, snr = snr, npad = npad, xysize = xysize, zsize = zsize)
	else:    vol_odd = recons3d_4nn_MPI(myid, td, symmetry=sym, snr = snr, npad = npad, xysize = xysize, zsize = zsize)

	del td

	if  myid == main_node  :
		vol_even = vol_even.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
		vol_even = sym_vol(vol_even, symmetry=sym)
		vol_odd  = vol_odd.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
		vol_odd = sym_vol(vol_odd, symmetry=sym)
	
		f = fsc_mask( vol_even, vol_odd, fscmask, 1.0, os.path.join(outdir, "hfsc.txt"))

	del vol_even, vol_odd
	"""






























































































































































"""53
		from sp_filter import filt_tanl
		fullvol0 = filt_tanl(fullvol0, 0.3, 0.2)
		fullvol0 = fullvol0.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
		fullvol0 = sym_vol(fullvol0, symmetry=sym)
		"""

































































































































"""54
	forg = []
	for ivol in xrange(nfils):
		for im in xrange(inds[ivol][0],inds[ivol][1]):
				T = data[im].get_attr("xform.projection")
				d = T.get_params('spider')
				forg.append([d["phi"],d["theta"],d["psi"],-d["tx"],-d["ty"],data[im].get_attr('ID')])
	write_text_row(forg,"pheader%03d.txt"%myid)
	"""















































































































































































































































































































































































































































'''55
	def rot2pad(imi, alpha=0.0, sx=0.0, sy=0.0):
		from sp_utilities    import pad
		from sp_fundamentals import rot_shift2D
		lnx = imi.get_xsize()
		lny = imi.get_ysize()
		ln = max(lnx,lny)
		if lnx == lny: return rot_shift2D(imi,alpha,sx,sy)
		else:          return Util.window(rot_shift2D(pad(imi,ln,ln,1,"circumference"), alpha,sx,sy), lnx, lny,1, 0,0,0)
	'''





































































































































































































'''56
		# check previous max and if does not exist set it to -1.e23
		p = data[im].get_attr_default('previousmax',-1.0e23)
		if( p == -1.0e23 ):
			resetatone = True
			data[im].set_attr('previousmax', p)
		'''
















































































'''57
			if doExhaustive:
				Util.constrained_helix_exhaustive(ldata, fdata[indcs[ifil][0]:indcs[ifil][1]], refproj, rotproj, [float(dp), float(dphi), float(rise), float(delta)], [int(nphi), int(phiwobble), int(rng), int(ywobble), int(Dsym), int(nwx), int(nwy), int(nwxc), int(nwyc)], FindPsi, float(psi_max), crefim, numr, int(maxrin), mode, int(cnx), int(cny))
				terminate = 0
			else:
				tempch = Util.constrained_helix(ldata, fdata[indcs[ifil][0]:indcs[ifil][1]], refproj, rotproj, [float(dp), float(dphi), float(rise), float(delta)], [int(nphi), int(phiwobble), int(rng), int(ywobble), int(Dsym), int(nwx), int(nwy), int(nwxc), int(nwyc)], FindPsi, float(psi_max), crefim, numr, int(maxrin), mode, int(cnx), int(cny))
				#print "tempch, Iter, myid: ", tempch, Iter, myid
				#tempch = constrained_helix_SHC(ldata, fdata[indcs[ifil][0]:indcs[ifil][1]], refproj, rotproj,  dp, dphi, rise, delta ,  nphi, phiwobble, rng, ywobble, Dsym, nwx, nwy, nwxc, nwyc , FindPsi, psi_max, crefim, numr, maxrin, mode, cnx, cny, myid, main_node)
				if tempch > -1:
					#if myid == main_node:
					#	print_msg("tempch %d\n"%tempch)
					terminate = 0
			'''




















"""58
		#  Should we continue??
		terminate = mpi_reduce(terminate, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
		terminate = mpi_bcast(terminate, 1, MPI_INT, 0, MPI_COMM_WORLD)
		terminate = int(terminate[0])
		if terminate == nproc:
			if myid == main_node: print_end_msg("helicon_MPI")
			return
		"""














































































































































































































































































































































































"""59
			# If the previous iteration did a reconstruction, then generate new refrings
			if ( (Iter - 1) % search_iter == 0):


				volft,kb = prep_vol( vol )
				#  What about cushion for a neighborhood?  PAP 06/04/2014
				refrings = prepare_refffts( volft, kb, data_nn,data_nn,nz, segmask, delta[N_step], \
					MPI=True, psimax=psi_max, psistep=psistep, initial_theta =initial_theta, delta_theta = delta_theta)
				#refrings = prepare_refrings2(  volft, kb, nmax, segmask, delta[N_step], ref_a, symref, numr, MPI = True, phiEqpsi = "Zero", initial_theta =initial_theta, delta_theta = delta_theta)
				del volft,kb

				if myid== main_node:
					print_msg( "Time to prepare rings: %d\n" % (time()-start_time) )
					start_time = time()
			"""



















































































































































































































































































































































































































































































































"""60
			# If the previous iteration did a reconstruction, then generate new refrings
			if ( (Iter - 1) % search_iter == 0):


				volft,kb = prep_vol( vol )
				#  What about cushion for a neighborhood?  PAP 06/04/2014
				refrings = prepare_refffts( volft, kb, data_nn,data_nn,nz, segmask, delta[N_step], \
					MPI=True, psimax=psi_max, psistep=psistep, initial_theta =initial_theta, delta_theta = delta_theta)
				#refrings = prepare_refrings2(  volft, kb, nmax, segmask, delta[N_step], ref_a, symref, numr, MPI = True, phiEqpsi = "Zero", initial_theta =initial_theta, delta_theta = delta_theta)
				del volft,kb

				if myid== main_node:
					print_msg( "Time to prepare rings: %d\n" % (time()-start_time) )
					start_time = time()
			"""





















































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































"""61
def setfilori_MA(fildata, pixel_size, dp, dphi):
	from sp_utilities		import get_params_proj, set_params_proj, get_dist
	from sp_applications	import filamentupdown
	from copy 			import copy
	from math 			import atan2, sin, cos, pi

	#if sym != 'c1':
	#	ERROR("does not handle any point-group symmetry other than c1 for the time being.", 'setfilori_MA')

	rise 	= dp/pixel_size
	ddphi   = pixel_size/dp*dphi
	ns 		= len(fildata)
	qv 		= pi/180.0

	phig 	= [0.0]*ns # given phi
	psig 	= [0.0]*ns # given psi
	yg 		= [0.0]*ns # given y
	xg 		= [0.0]*ns # given x
	thetag	= [0.0]*ns # given theta

	coords = [[] for im in xrange(ns)]
	for i in xrange(ns):
		coords[i] = fildata[i].get_attr('ptcl_source_coord')
		phig[i], thetag[i], psig[i] , xg[i], yg[i] = get_params_proj(fildata[i])
		#print "%3d  %5.1f   %5.1f   %5.1f   %5.1f   %5.1f"%(i,phig[i], thetag[i], psig[i] , xg[i], yg[i])
		if( abs(psig[i] - psig[0]) )> 90.0:
			ERROR('PSI should be pointing in the same direction for all segments belonging to same filament', 'setfilori_MA')

	# forward predicted psi of segment i is current psi of segment i + 1
	# backward predicted psi of segment i is current psi of segment i - 1
	# For now set the new psi to the psi closest to the two predicted psi
	conspsi = [0.0]*ns
	for i in xrange(1,ns-1):
		ff = psig[i+1]*qv
		bb = psig[i-1]*qv
		conspsi[i] = (atan2(  (sin(ff)+sin(bb)) , (cos(ff)+cos(bb)) )/qv)%360.0
	#  border psi values
	ff = psig[1]*qv
	bb = psig[0]*qv
	conspsi[0]    = (atan2(  (sin(ff)+sin(bb)) , (cos(ff)+cos(bb)) )/qv)%360.0
	ff = psig[ns-1]*qv
	bb = psig[ns-2]*qv
	conspsi[ns-1] = (atan2(  (sin(ff)+sin(bb)) , (cos(ff)+cos(bb)) )/qv)%360.0

	# predict theta in same way as psi. 
	constheta = [0.0]*ns
	for i in xrange(1,ns-1):
		ff = thetag[i+1]*qv
		bb = thetag[i-1]*qv
		constheta[i] = (atan2(  (sin(ff)+sin(bb)) , (cos(ff)+cos(bb)) )/qv)%360.0
	#  border theta values
	ff = thetag[1]*qv
	bb = thetag[0]*qv
	constheta[0]    = (atan2(  (sin(ff)+sin(bb)) , (cos(ff)+cos(bb)) )/qv)%360.0
	ff = thetag[ns-1]*qv
	bb = thetag[ns-2]*qv
	constheta[ns-1] = (atan2(  (sin(ff)+sin(bb)) , (cos(ff)+cos(bb)) )/qv)%360.0


	#  X-shift (straightforward)
	consx = [0.0]*ns
	for i in xrange(1,ns-1): consx[i] = (xg[i+1] + xg[i-1])/2.0
	consx[0]    = (xg[1] + xg[0])/2.0
	consx[ns-1] = (xg[ns-1] + xg[ns-2])/2.0

	#  Y-shift - include non-integer distance.

	# Predict phi angles based on approximate distances calculated using theta=90 between nearby segments.
	# (The other option is to calculate distance between segments using predicted theta. 
	#   But for now, just use theta=90 since for reasonable deviations of theta from 90 (say +/- 10 degrees), 
	#   and an intersegment distance of say >= 15 pixels, the difference theta makes is minimal.)

	sgn = (1 - fildata[0].get_attr("updown")*2)

	consy   = [0.0]*ns
	consphi = [0.0]*ns
	
	
	#per = 0.0
	#yer = 0.0
	
	for i in xrange(ns):
		if( i>0 ):
			# from i-1 predict i
			dist = get_dist(coords[i-1], coords[i])
			qd = round((yg[i-1] + dist)/rise)
			yb   = yg[i-1] + dist - rise*qd
			#print  " YB ",yg[i-1],dist,qd,yb
			phib = qv*((phig[i-1] + sgn*dphi*qd)%360.0)
		else:
			yb   = yg[0]
			phib = qv*phig[0]

		if( i < ns-1 ):
			# from i+1 predict i
			dist = get_dist(coords[i+1], coords[i])
			qd = round((yg[i+1] - dist)/rise)
			yf   = yg[i+1] - dist - rise*qd
			#print  " YF ",yg[i+1],dist,qd,yf
			phif = qv*((phig[i+1] + sgn*dphi*qd)%360.0)
		else:
			yf   = yg[-1]
			phif = qv*phig[-1]

		consy[i] = (yb + yf)/2.0
		#print "    YP",yg[max(0,i-1)], yg[min(i+1,ns-1)], yb,yf, yg[i],consy[i],"   >>>   ",yg[i]-consy[i]
		consphi[i] = (atan2(  (sin(phib)+sin(phif)) , (cos(phib)+cos(phif)) )/qv)%360.0
		#print " PHI ",i, round(phig[max(0,i-1)],1), round(phig[i],1), round(phig[min(i+1,ns-1)],1), \
		#	round(phib/qv,1), round(phif/qv,1), round(consphi[i],1),\
		#	"    >><>>>>> ",round((round((phig[i]-consphi[i])*100)/100.)%360.,1)
		#yer += abs(yg[i]-consy[i])
		#per += abs(phig[i]-consphi[i])
		set_params_proj(fildata[i], [consphi[i], constheta[i], conspsi[i], consx[i], consy[i]])
	#print yer, per
	return
"""
































































































































































































































































































































































































































"""62
		vol = sym_vol(vol, symmetry=sym)
		ref_data = [vol, mask3D]
		#if  fourvar:  ref_data.append(varf)
		vol = user_func(ref_data)
		vol = vol.helicise(pixel_size, dp, dphi, fract, rmax, rmin)
		vol = sym_vol(vol, symmetry=sym)
		drop_image(vol, os.path.join(outdir, "volf.hdf"))
		print_msg("\nSymmetry search and user function time = %d\n"%(time()-start_time))
		start_time = time()
		"""










































































































"""63
		import os
		outdir = "./"
		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
		"""















































































































































































"""  Have to think about it PAP64
					if( (historyofchanges[-1]-historyofchanges[0])/2/(historyofchanges[-1]+historyofchanges[0]) <0.05 ):
						terminate = 1
						log.add("...............")
						log.add(">>>>>>>>>>>>>>>   Will terminate as orientations do not improve anymore")
					"""






































































"""65
			#=========================================================================
			if(False):  #total_iter%1 == 5 or terminate):
				# gather parameters
				params = []
				previousmax = []
				for im in data:
					t = get_params_proj(im)
					params.append( [t[0], t[1], t[2], t[3]/shrinkage, t[4]/shrinkage] )
					#if(t[3] >0.0 or t[4]>0.0):  print  "  ERRROR  ",t
					previousmax.append(im.get_attr("previousmax"))
				assert(nima == len(params))
				params = wrap_mpi_gatherv(params, 0, mpi_comm)
				if myid == 0:
					assert(total_nima == len(params))
				previousmax = wrap_mpi_gatherv(previousmax, 0, mpi_comm)
				if myid == main_node:
					from sp_utilities import write_text_row, write_text_file
					write_text_row(params, "soft/params%04d.txt"%total_iter)
					write_text_file(previousmax, "soft/previousmax%04d.txt"%total_iter)


				del previousmax, params
				i = 1
				while data[0].has_attr("xform.projection" + str(i)):
					params = []
					previousmax = []
					for im in data:

						try:
							#print  im.get_attr("xform.projection" + str(i))
							t = get_params_proj(im,"xform.projection" + str(i))
						except:
							print " NO XFORM  ",myid, i,im.get_attr('ID')
							from sys import exit
							exit()

						params.append( [t[0], t[1], t[2], t[3]/shrinkage, t[4]/shrinkage] )
					assert(nima == len(params))
					params = wrap_mpi_gatherv(params, 0, mpi_comm)
					if myid == 0:
						assert(total_nima == len(params))
					if myid == main_node:
						write_text_row(params, "soft/params-%04d-%04d.txt"%(i,total_iter))
					del previousmax, params
					i+=1


			if( ( terminate or (Iter == max_iter) ) and (myid == main_node) ):
				if( type(stack) is types.StringType ):
					from EMAN2 import Vec2f, Transform
					from EMAN2db import db_open_dict
					DB = db_open_dict(stack)
					for im in xrange(len(params)):
						t = Transform({"type":"spider","phi":params[im][0],"theta":params[im][1],"psi":params[im][2]})
						t.set_trans(Vec2f(-params[im][3], -params[im][4]))
						DB.set_attr(particle_ids[im], "xform.projection", t)
					DB.close()
				else:
					for im in xrange(len(params)): set_params_proj(stack[particle_ids[im]], params[im])
			"""








































































"""66
	if myid == main_node:
		os.mkdir(outdir)
		import sp_global_def
		sp_global_def.LOGFILE =  os.path.join(outdir, sp_global_def.LOGFILE)
		print_begin_msg("local_ali3d_MPI")
		import sp_user_functions
		user_func = sp_user_functions.factory[user_func_name]
		if CTF:
			ima = EMData()
			ima.read_image(stack, 0)
			ctf_applied = ima.get_attr_default("ctf_applied", 0)
			del ima
			if ctf_applied == 1:  ERROR("local_ali3d does not work for CTF-applied data", "local_ali3d_MPI", 1,myid)
	mpi_barrier(mpi_comm)
	"""












































































"""67
	if myid == main_node:
		import sp_user_functions
		user_func = sp_user_functions.factory[user_func_name]

		print_msg("Input stack                 : %s\n"%(stack))
		print_msg("Output directory            : %s\n"%(outdir))
		print_msg("Maskfile                    : %s\n"%(Tracker["constants"]["mask3D"]))
		print_msg("Outer radius                : %i\n"%(last_ring))
		print_msg("Angular search range        : %s\n"%(delta))
		print_msg("Shift search range          : %f\n"%(ts))
		print_msg("Maximum iteration           : %i\n"%(max_iter))
		print_msg("Center type                 : %i\n"%(center))
		print_msg("CTF correction              : %s\n"%(CTF))
		print_msg("Signal-to-Noise Ratio       : %f\n"%(snr))
		print_msg("Symmetry group              : %s\n"%(sym))
		print_msg("Chunk size                  : %f\n\n"%(chunk))
		print_msg("User function               : %s\n"%(user_func_name))
	"""























































































































































































































































































"""68
		if( type(stack) is types.StringType ):
			from EMAN2 import Vec2f, Transform
			from EMAN2db import db_open_dict
			DB = db_open_dict(stack)
			for im in xrange(len(params)):
				t = Transform({"type":"spider","phi":params[im][0],"theta":params[im][1],"psi":params[im][2]})
				t.set_trans(Vec2f(-params[im][3], -params[im][4]))
				DB.set_attr(particle_ids[im], "xform.projection", t)
			DB.close()
		else:
			for im in xrange(len(params)): set_params_proj(stack[particle_ids[im]], params[im])

		log.add("Time to write header information= %d\n"%(time()-start_time))
		"""





































































































































































































'''69
	if fourvar:
		from sp_reconstruction import rec3D_MPI
		from sp_statistics     import varf3d_MPI
		#  Compute Fourier variance
		vol, fscc = rec3D_MPI(data,snr,sym,fscmask,os.path.join(outdir, "resolution0000"), myid, main_node, finfo=frec, npad=npad)
		varf = varf3d_MPI(data, os.path.join(outdir, "ssnr0000"), None, vol, last_ring, 1.0, 1, CTF, 1, sym, myid)
		if myid == main_node:   
			varf = 1.0/varf
			varf.write_image( os.path.join(outdir,"varf0000.hdf") )
	else:
		varf = None
	'''




















































































































































'''70
					if CTF:
						ref = fft(filt_ctf( prgl( volft, [phi,tht,psi,-s2x,-s2y],1,False), ctf ))
					else:
						ref = prgl( volft, [phi,tht,psi,-s2x,-s2y],1)
					peak = ref.cmp("ccc",data[im],{"mask":mask2D, "negative":0})
					'''
























'''71
					#  FSC distance
					#  Ref is in reciprocal space
					if CTF:
						ref = filt_ctf( prgl( volft, [phi,tht,psi,-s2x,-s2y], 1, False), ctf )
					else:
						ref = prgl( volft, [phi,tht,psi,-s2x,-s2y], 1, False)
					from sp_statistics import fsc
					if(focus):
						mask2D = binarize( prgl( focus, [phi,tht,psi,-s2x,-s2y]), 1)
						tempx = fsc(ref, fft(data[im]*mask2D))[1]
					else:
						tempx = fsc(ref, data[im])[1]
					peak = sum(tempx[1:highres[iref]])/highres[iref]

					if not(finfo is None):
						finfo.write( "ID,iref,peak: %6d %d %8.5f\n" % (list_of_particles[im],iref,peak) )
				else:
					if an[N_step] == -1:
						peak, pixel_error = proj_ali_incore(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step])
					else:
						peak, pixel_error = proj_ali_incore_local(data[im], refrings, list_of_reference_angles, numr,\
																	xrng[N_step], yrng[N_step], step[N_step], an[N_step])
					if not(finfo is None):
						phi,tht,psi,s2x,s2y = get_params_proj(data[im])
						finfo.write( "ID,iref,peak,trans: %6d %d %f %f %f %f %f %f\n"%(list_of_particles[im],iref,peak,phi,tht,psi,s2x,s2y) )
						finfo.flush()
					'''





























































































































































































"""72
		if myid == main_node:
			refdata    = [None]*7
			refdata[0] = numref
			refdata[1] = outdir
			refdata[2] = None
			refdata[3] = total_iter
			refdata[4] = varf
			refdata[5] = mask3D
			refdata[6] = Tracker["low_pass_filter"]#(runtype=="REFINEMENT") # whether align on 50S, this only happens at refinement step
			user_func( refdata )
		"""




































"""73
	if myid == main_node:
		from sp_utilities import file_type
		if file_type(stack) == "bdb":
			from sp_utilities import recv_attr_dict_bdb
			recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
		else:
			from sp_utilities import recv_attr_dict
			recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	else:		send_attr_dict(main_node, data, par_str, image_start, image_end)
	"""
































































































































































































'''74
	if(myid == main_node):	
		total_nima = EMUtil.get_image_count(stack)
		list_of_particles = range(total_nima)
	
	else:
		total_nima =0

	total_nima = bcast_number_to_all(total_nima, source_node = main_node)

	if(myid != main_node):
		list_of_particles = [-1]*total_nima

	list_of_particles = bcast_list_to_all(list_of_particles, myid,  source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	'''



























































































































































'''75
					if CTF:
						ref = fft(filt_ctf( prgl( volft, [phi,tht,psi,-s2x,-s2y],1,False), ctf ))
					else:
						ref = prgl( volft, [phi,tht,psi,-s2x,-s2y],1)
					if(focus != None):  mask2D = binarize( prgl( focus, [phi,tht,psi,-s2x,-s2y]),1)  #  Should be precalculated!!
					peak = ref.cmp("ccc",data[im],{"mask":mask2D, "negative":0})
					'''























'''76
					#  FSC distance
					#  Ref is in reciprocal space
					if CTF:
						ref = filt_ctf( prgl( volft, [phi,tht,psi,-s2x,-s2y], 1, False), ctf )
					else:
						ref = prgl( volft, [phi,tht,psi,-s2x,-s2y], 1, False)
					if(focus):
						mask2D = binarize( prgl( focus, [phi,tht,psi,-s2x,-s2y]), 1)
						tempx = fsc(ref, fft(data[im]*mask2D))[1]
					else:
						tempx = fsc(ref, data[im])[1]
					peak = sum(tempx[1:highres[iref]])/highres[iref]
					'''








































































"""77
		#  Trying to use ISAC code for EQ-Kmeans  PAP 03/21/2015
		if myid == main_node:

			for imrefa in xrange(numrefang):
				from sp_utilities import findall
				N = findall(imrefa, assigntorefa)
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

			e b = empty( (numref, total_nima), dtype = float32)
			for ima in xrange(total_nima):
				for iref in xrange(numref):  d[iref][ima] = dtot[iref][N[ima]]
			id_list_long = Util.assign_groups(str(d.__array_interface__['data'][0]), numref, nima) # string with memory address is passed as parameters
			del d
			id_list = [[] for i in xrange(numref)]
			maxasi = total_nima/numref
			for i in xrange(maxasi*numref):
				id_list[i/maxasi].append(id_list_long[i])
			for i in xrange(total_nima%maxasi):
				id_list[id_list_long[-1]].append(id_list_long[maxasi*numref+i])
			for iref in xrange(numref):
				id_list[iref].sort()

			assignment = [0]*total_nima
			for iref in xrange(numref):
				for im in id_list[iref]: assignment[im] = iref
		else:
			assignment = [0]*total_nima
		mpi_barrier(MPI_COMM_WORLD)
		#belongsto = mpi_bcast(belongsto, nima, MPI_INT, main_node, MPI_COMM_WORLD)
		#belongsto = map(int, belongsto)
		"""













































































"""78
		if myid == main_node:
			assignment = [0]*total_nima
			for iref in xrange(numref):
				for im in xrange(len(asi[iref])):
					assignment[asi[iref][im]] = iref
			del asi
		"""

'''79
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


























































































































































"""80
		if runtype=="REFINEMENT":
			if fourvar:
				varf = varf3d_MPI(data, os.path.join(outdir, "ssnr%04d"%total_iter), None,sumvol,last_ring, 1.0, 1, CTF, 1, sym, myid)
				if myid == main_node:   
					varf = 1.0/varf
					varf.write_image( os.path.join(outdir,"varf%04d.hdf"%total_iter) )
		"""

"""81
			frcs={}
			for iref in xrange(numref):
				frc=read_text_file(os.path.join(outdir, "resolution_%02d_%04d"%(iref, total_iter)),-1)
				frcs[iref]=frc
			refdata = [None]*7
			refdata[0] = numref
			refdata[1] = outdir
			refdata[2] = None #frequency_low_pass
			refdata[3] = total_iter
			refdata[4] = None
			refdata[5] = mask3D
			refdata[6] = Tracker["low_pass_filter"]#(runtype=="REFINEMENT") # whether to align on 50S, this only happens at refinement step
		"""





"""82
	        	if myid == main_node:
				from sp_utilities import file_type
	        		if(file_type(stack) == "bdb"):
	        			from sp_utilities import recv_attr_dict_bdb
	        			recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        		else:
	        			from sp_utilities import recv_attr_dict
	        			recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        	else:		send_attr_dict(main_node, data, par_str, image_start, image_end)
			"""














































































































































































































'''83
	if(myid == main_node):	
		total_nima = EMUtil.get_image_count(stack)
		list_of_particles = range(total_nima)
	
	else:
		total_nima =0

	total_nima = bcast_number_to_all(total_nima, source_node = main_node)

	if(myid != main_node):
		list_of_particles = [-1]*total_nima

	list_of_particles = bcast_list_to_all(list_of_particles, myid,  source_node = main_node)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)
	'''



















































































































































'''84
					if CTF:
						ref = fft(filt_ctf( prgl( volft, [phi,tht,psi,-s2x,-s2y],1,False), ctf ))
					else:
						ref = prgl( volft, [phi,tht,psi,-s2x,-s2y],1)
					if(focus != None):  mask2D = binarize( prgl( focus, [phi,tht,psi,-s2x,-s2y]),1)  #  Should be precalculated!!
					peak = ref.cmp("ccc",data[im],{"mask":mask2D, "negative":0})
					'''























'''85
					#  FSC distance
					#  Ref is in reciprocal space
					if CTF:
						ref = filt_ctf( prgl( volft, [phi,tht,psi,-s2x,-s2y], 1, False), ctf )
					else:
						ref = prgl( volft, [phi,tht,psi,-s2x,-s2y], 1, False)
					if(focus):
						mask2D = binarize( prgl( focus, [phi,tht,psi,-s2x,-s2y]), 1)
						tempx = fsc(ref, fft(data[im]*mask2D))[1]
					else:
						tempx = fsc(ref, data[im])[1]
					peak = sum(tempx[1:highres[iref]])/highres[iref]
					'''















































































"""86
		#  Trying to use ISAC code for EQ-Kmeans  PAP 03/21/2015
		if myid == main_node:

			for imrefa in xrange(numrefang):
				from sp_utilities import findall
				N = findall(imrefa, assigntorefa)
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

			e b = empty( (numref, total_nima), dtype = float32)
			for ima in xrange(total_nima):
				for iref in xrange(numref):  d[iref][ima] = dtot[iref][N[ima]]
			id_list_long = Util.assign_groups(str(d.__array_interface__['data'][0]), numref, nima) # string with memory address is passed as parameters
			del d
			id_list = [[] for i in xrange(numref)]
			maxasi = total_nima/numref
			for i in xrange(maxasi*numref):
				id_list[i/maxasi].append(id_list_long[i])
			for i in xrange(total_nima%maxasi):
				id_list[id_list_long[-1]].append(id_list_long[maxasi*numref+i])
			for iref in xrange(numref):
				id_list[iref].sort()

			assignment = [0]*total_nima
			for iref in xrange(numref):
				for im in id_list[iref]: assignment[im] = iref
		else:
			assignment = [0]*total_nima
		mpi_barrier(MPI_COMM_WORLD)
		#belongsto = mpi_bcast(belongsto, nima, MPI_INT, main_node, MPI_COMM_WORLD)
		#belongsto = map(int, belongsto)
		"""













































































"""87
		if myid == main_node:
			assignment = [0]*total_nima
			for iref in xrange(numref):
				for im in xrange(len(asi[iref])):
					assignment[asi[iref][im]] = iref
			del asi
		"""

'''88
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











































































































































"""89
		if runtype=="REFINEMENT":
			if fourvar:
				varf = varf3d_MPI(data, os.path.join(outdir, "ssnr%04d"%total_iter), None,sumvol,last_ring, 1.0, 1, CTF, 1, sym, myid)
				if myid == main_node:   
					varf = 1.0/varf
					varf.write_image( os.path.join(outdir,"varf%04d.hdf"%total_iter) )
		"""

"""90
			frcs={}
			for iref in xrange(numref):
				frc=read_text_file(os.path.join(outdir, "resolution_%02d_%04d"%(iref, total_iter)),-1)
				frcs[iref]=frc
			refdata = [None]*7
			refdata[0] = numref
			refdata[1] = outdir
			refdata[2] = None #frequency_low_pass
			refdata[3] = total_iter
			refdata[4] = None
			refdata[5] = mask3D
			refdata[6] = Tracker["low_pass_filter"]#(runtype=="REFINEMENT") # whether to align on 50S, this only happens at refinement step
		"""





"""91
	        	if myid == main_node:
				from sp_utilities import file_type
	        		if(file_type(stack) == "bdb"):
	        			from sp_utilities import recv_attr_dict_bdb
	        			recv_attr_dict_bdb(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        		else:
	        			from sp_utilities import recv_attr_dict
	        			recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	        	else:		send_attr_dict(main_node, data, par_str, image_start, image_end)
			"""




























