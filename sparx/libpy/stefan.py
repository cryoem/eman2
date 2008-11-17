from EMAN2_cppwrap import *
from global_def import *
"""
def ali3d_d_MPI_Stefan(stack, ref_vol, outdir, maskfile = None, ir = 1, ou = -1, rs = 1, 
            xr = "4 2 2 1", yr = "-1", ts = "1 1 0.5 0.25", delta = "10 6 4 4", an="-1",
	    center = 1.0, maxit = 5, CTF = False, snr = 1.0,  ref_a="S", sym="c1"):
	from utilities      import model_circle, reduce_EMData_to_root, bcast_EMData_to_all, drop_image
	from utilities      import bcast_list_to_all, bcast_string_to_all, get_image, get_input_from_string
	from utilities      import get_arb_params, set_arb_params, drop_spider_doc,recv_attr_dict, send_attr_dict
	from utilities      import drop_spider_doc,get_im
	from filter	    import filt_params,filt_btwl, filt_from_fsc, filt_table
	from alignment	    import proj_ali_incore
	from random	    import randint
	from statistics     import fsc
	from fundamentals   import fshift,rot_avg_image
	import os
	import types
	from string         import replace
	from reconstruction import rec3D_MPI, rec3D_MPI_noCTF
	
	# 2D alignment using rotational ccf in polar coords and linear
	# interpolation	
	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	if(myid == main_node):
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)
        #info_file = replace("progress%4d"%myid, ' ', '0') 
        #finfo = open(info_file, 'w')
	finfo = None

	max_iter    = int(maxit)
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
	nimage_per_node = nima/number_of_proc+1
	image_start     = myid * nimage_per_node
	if(myid == number_of_proc-1):  image_end = nima
	else                        :  image_end = image_start + nimage_per_node
	if(myid == main_node)       :  drop_image(vol, os.path.join(outdir,"ref_volf00.spi"), "s")
	if(CTF):
		ctf_dicts = ["defocus", "Cs", "voltage", "Pixel_size", "amp_contrast", "B_factor", "ctf_applied" ]
		#                 0      1         2            3               4                  5               6
		from reconstruction import rec3D_MPI
	else:
		from reconstruction import rec3D_MPI_noCTF
	data = []
	for im in xrange(image_start, image_end):
		ima = EMData()
		ima.read_image(stack, im)
		if(CTF):
			ctf_params = get_arb_params(ima, ctf_dicts)
			if(ctf_params[6] == 0):
				from filter import filt_ctf
				ima = filt_ctf(ima, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
				ima.set_attr('ctf_applied', 1)
		data.append(ima)
	# do the projection matching
	from  string  import replace

	for N_step in xrange(len(xrng)):
 		for Iter in xrange(max_iter):
			#from filter import filt_gaussinv
			#vol = filt_gaussinv(vol, 0.175, True)
			#from filter import filt_btwl,filt_btwo
			#vol = filt_btwl(filt_btwo(vol,0.25,0.35,0.1),0.35,0.45)
			if(myid == main_node):
				#drop_image(vol, os.path.join(outdir, replace("vhl%4d.spi"%(N_step*max_iter+Iter+1),' ','0')), "s")
				print " ali3d_d_MPI: ITERATION #",N_step*max_iter + Iter+1
			if(an[N_step] == -1): proj_ali_incore(vol, mask3D, data, first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], ref_a, sym, CTF, finfo)
			else:                proj_ali_incore_local(vol, mask3D, data, first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], an[N_step], ref_a, sym, CTF, finfo)

			if(CTF): 
				 vol, fscc_wos = rec3D_MPI(data, snr, sym, mask3D, os.path.join(outdir, replace("resolution%4d"%(N_step*max_iter+Iter+1),' ','0')), myid, main_node)
                    		 volws, fscc = rec3D_MPI(data, snr, "d4", mask3D, os.path.join(outdir, replace("resolution%4d"%(N_step*max_iter+Iter+1),' ','0')), myid, main_node)

			else:    
				 vol, fscc_wos = rec3D_MPI_noCTF(data, sym, mask3D, os.path.join(outdir, replace("resolution%4d"%(N_step*max_iter+Iter+1),' ','0')), myid, main_node)
				 volws, fscc = rec3D_MPI_noCTF(data, "d4", mask3D, os.path.join(outdir, replace("resolution%4d"%(N_step*max_iter+Iter+1),' ','0')), myid, main_node)

			mpi_barrier(MPI_COMM_WORLD)
			if(myid == main_node):
				drop_image(vol,  os.path.join(outdir, replace("vol%4d.spi"%(N_step*max_iter+Iter+1),' ','0')), "s")
				#filt = filt_from_fsc(fscc, 0.05)
				#vol  = filt_table(vol, filt)
				# here figure the filtration parameters and filter vol for the  next iteration
				#lk = 2
				#while(fscc[1][lk] >0.98 and fscc[0][lk]<0.25):
				#	lk+=1
				#fl = fscc[0][lk]
				#fh = min(fl+0.1,0.49)

				from utilities import get_im
				e = get_im("maskp_1RXOff.spi")
				volws *= e
				f = get_im("maskn_1RXOff.spi")
				vol *= f
				vol += volws

				fl, fh = filt_params(fscc,0.95,0.3)
				vol = filt_btwl(vol, fl, fh)
				if(center == 1):
					cs   = vol.phase_cog()
					vol  = fshift(vol, -cs[0], -cs[1] -cs[2])
				drop_image(vol,  os.path.join(outdir, replace("volf%4d.spi"%(N_step*max_iter+Iter+1),' ','0')), "s")
			bcast_EMData_to_all(vol, myid, main_node)
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
	#  At this moment parameters are not written to the file
	par_str = ['phi', 'theta', 'psi', 's2x', 's2y']
	if(myid == main_node): recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	else: send_attr_dict(main_node, data, par_str, image_start, image_end)
	
	
def ali3d_m_MPI_Stefan(stack, ref_vol, outdir, maskfile = None, ir=1, ou=-1, rs=1, 
            xr              ="4 2  2  1",      yr="-1",
            translation_step="1 1 0.5 0.25",   delta="10  6  4  4",
	    center = 0, maxit= 5, CTF = False, snr = 1.0,  ref_a="S", symmetry="c1"):
	from utilities      import model_circle, reduce_EMData_to_root, bcast_EMData_to_all, drop_image
	from utilities      import bcast_list_to_all, bcast_string_to_all, get_image, get_input_from_string
	from utilities      import get_arb_params, set_arb_params, drop_spider_doc
	from filter	    import filt_params,filt_btwl, filt_from_fsc2, filt_table
	from alignment	    import proj_ali_incore_index
	from random         import randint
	from statistics     import fsc
	from fundamentals   import fshift
	import os
	import types
	from string         import replace
	# 2D alignment using rotational ccf in polar coords and linear
	# interpolation	
	number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
	myid           = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0
	if(myid == main_node):
		if os.path.exists(outdir):  os.system('rm -rf '+outdir)
		os.mkdir(outdir)

	max_iter    = int(maxit);
	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(translation_step)
	delta       = get_input_from_string(delta)
	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	numref            = EMUtil.get_image_count(ref_vol)
	volref = []
	volrefws = []
	for  iref in xrange(numref):
		vol     = EMData()
		vol.read_image(ref_vol, iref)
		volref.append(vol)
	del vol
	nx      = volref[0].get_xsize()
	if last_ring < 0:	last_ring = int(nx/2) - 2
	if(maskfile):
		if(type(maskfile) is types.StringType):	 mask3D = get_image(maskfile)
		else: 	                                mask3D = maskfile
	else        :   mask3D = model_circle(last_ring, nx, nx, nx)
	nima            = EMUtil.get_image_count(stack)
	nimage_per_node = nima/number_of_proc
	image_start     = myid * nimage_per_node
	if(myid == number_of_proc-1):  image_end = nima
	else                        :  image_end = image_start + nimage_per_node	
	if(myid == main_node)       :
		for  iref in xrange(numref):
			volref[iref].write_image(os.path.join(outdir,replace("ref_vol%2d.spi"%iref,' ','0')))
	if(CTF):
		ctf_dicts = ["defocus", "Cs", "voltage", "Pixel_size", "amp_contrast", "B_factor", "ctf_applied" ]
		#                        0      1              2            3               4                  5               6
		from reconstruction import rec3D_MPI
	else:
		from reconstruction import rec3D_MPI_noCTF
	data = []
	for im in xrange(image_start, image_end):
		ima = EMData()
		ima.read_image(stack, im)
		ima.set_attr_dict({'group':im%numref})  # this pretends to be general.  However, for the time being no CTF
		if(CTF):
			ctf_params = get_arb_params(ima, ctf_dicts)
			if(im == image_start): data_had_ctf = ctf_params[6]
			if(ctf_params[6] == 0):
				from filter import filt_ctf
				ima = filt_ctf(ima, ctf_params[0], ctf_params[1], ctf_params[2], ctf_params[3], ctf_params[4], ctf_params[5])
				ima.set_attr('ctf_applied', 1)
				#set_arb_params(ima, ctf_dicts)  #I am not sure whether this is needed
		data.append(ima)
	# do the projection matching
	for N_step in xrange(len(xrng)):
 		for Iter in xrange(max_iter):
			for im in xrange(image_start, image_end):  data[im-image_start].set_attr_dict({'peak':-1.0e23})
			for iref in xrange(numref):
				#print " ali3d_d_MPI: ITERATION #",N_step*max_iter + Iter+1
				proj_ali_incore_index(volref[iref], iref, mask3D, data, first_ring, last_ring, rstep, xrng[N_step], yrng[N_step], step[N_step], delta[N_step], ref_a, symmetry, CTF)
			mpi_barrier(MPI_COMM_WORLD)
			for iref in xrange(numref):
				if(CTF): 
					volref[iref], fscc_wos = rec3D_MPI(data, snr, symmetry, mask3D, os.path.join(outdir, replace("resolution%4d"%(N_step*max_iter+Iter+1),' ','0')), myid, main_node)
					volrefws[iref], fscc   = rec3D_MPI(data, snr, "d4", mask3D, os.path.join(outdir, replace("resolution%4d"%(N_step*max_iter+Iter+1),' ','0')), myid, main_node)
				else:    
					volref[iref], fscc_wos = rec3D_MPI_noCTF(data, symmetry, mask3D, os.path.join(outdir, replace("resolution%7d"%((N_step*max_iter+Iter+1)*10 + iref),' ','0')), myid, main_node, 1.0, index = iref)
                                        volrefws[iref], fscc   = rec3D_MPI_noCTF(data, "d4", mask3D, os.path.join(outdir, replace("resolution%7d"%((N_step*max_iter+Iter+1)*10 + iref),' ','0')), myid, main_node, 1.0, index = iref)

				if(myid == main_node):
					drop_image(volref[iref],os.path.join(outdir, replace("vol%2d%4d.spi"%(iref, N_step*max_iter+Iter+1),' ','0')), "s")
					if(fscc[1][0] < 0.5):  fscc[1][0] = 1.0
					if(fscc[1][1] < 0.5):  fscc[1][1] = 1.0
					
					#filt = filt_from_fsc2(fscc, 0.075)
					# here figure the filtration parameters and filter vol for the  next iteration
					#fl, fh = filt_params(res)
					#lk = 2
					#while(fscc[1][lk] >0.98 and fscc[0][lk]<0.25):
					#	lk+=1
					#fl = fscc[0][lk]
					#fh = min(fl+0.1,0.49)
					#print "fl, fh, iter",fl,fh,Iter
					from utilities import get_im
					e = get_im("maskp_1RXOff.spi")
					volrefws[iref] *= e
					f = get_im("maskn_1RXOff.spi")
					volref[iref] *= f
					volref[iref] += volws[iref]
					
					fl, fh = filt_params(fscc,0.95,0.3)

					volref[iref] = filt_btwl(volref[iref], fl, fh)
					if(center == 1):
						cs   = volref[iref].phase_cog()
						volref[iref]  = fshift(volref[iref], -cs[0], -cs[1] -cs[2])
					drop_image(volref[iref],os.path.join(outdir, replace("volf%2d%4d.spi"%(iref, N_step*max_iter+Iter+1),' ','0')), "s")
			mpi_barrier(MPI_COMM_WORLD)
			for iref in xrange(numref): bcast_EMData_to_all(volref[iref], myid, main_node)
	mpi_barrier(MPI_COMM_WORLD)
	del  volref
	#vol, fscc = rec3D_MPI_noCTF(data, symmetry, mask3D, os.path.join(outdir, "resolution_merged"), myid, main_node, info=myinfo)
	#if(myid == main_node):  drop_image(vol,os.path.join(outdir, "vol_merged.spi"), "s")
	#  clean up
	#del vol
	# write out headers  and STOP, under MPI writing has to be done sequentially
	mpi_barrier(MPI_COMM_WORLD)
	if(CTF and data_had_ctf == 0):
		for im in xrange(image_start, image_end): data[im-image_start].set_attr('ctf_applied', 0)
	par_str = ['phi', 'theta', 'psi', 's2x', 's2y']
	if(myid == main_node): recv_attr_dict(main_node, stack, data, par_str, image_start, image_end, number_of_proc)
	else: send_attr_dict(main_node, data, par_str, image_start, image_end)
"""
