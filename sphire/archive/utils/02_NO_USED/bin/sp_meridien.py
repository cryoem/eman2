






























































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































from __future__ import print_function
def subdict(d,u):
	# substitute values in dictionary d by those given by dictionary u
	for q in u:  d[q] = u[q]
















































































def get_coarse_shifts():
	global Tracker, Blockdata

	k = int(numpy.ceil(Tracker["xr"]/Tracker["ts"]))
	radi = Tracker["xr"]*Tracker["xr"]

	indc = []
	for ix in range(-k,k+1,1):
		six = ix*Tracker["ts"]
		for iy in range(-k,k+1,1):
			siy = iy*Tracker["ts"]
			if(six*six+siy*siy <= radi):
				indc.append( [ix, iy] )

	cndc = []
	for ix in range(0,k+1,2):
		six = ix*Tracker["ts"]
		for iy in range(0,k+1,2):
			siy = iy*Tracker["ts"]
			if(six*six+siy*siy <= radi):
				cndc.append( [ix, iy] )

	for ix in range(1,len(cndc)):  cndc.append([-cndc[ix][0],-cndc[ix][1]])

	list_of_coarse_shifts = []
	for ix in range(len(cndc)):  list_of_coarse_shifts.append(indc.index(cndc[ix]))


	return list_of_coarse_shifts


























def prepdata_ali3d(projdata, rshifts, shrink, method = "DIRECT"):
	global Tracker, Blockdata
	pass#IMPORTIMPORTIMPORT from sp_fundamentals 	import prepi
	pass#IMPORTIMPORTIMPORT from sp_morphology 	import ctf_img_real
	#  Data is NOT CTF-applied.
	#  Data is shrank, in Fourier format
	data = [[] for i in range(len(projdata))]
	if Tracker["constants"]["CTF"]:
		nx = projdata[0].get_ysize()
		ctfs = [ sp_morphology.ctf_img_real(nx, q.get_attr('ctf')) for q in projdata ]
	else:  ctfs = None
	if Blockdata["bckgnoise"] :
		bckgnoise = [q.get_attr("bckgnoise") for q in projdata ]
	else:  bckgnoise = None
	for kl in range(len(projdata)-1,-1,-1):  #  Run backwards and keep deleting projdata, it is not needed anymore
		#  header shifts were shrank in get_shrink_data, shifts were already pre-applied, but we will leave the code as is.
		#phi, theta, psi, sxs, sys = get_params_proj(projdata[kl])
		particle_group = projdata[kl].get_attr("particle_group")
		ds = projdata[kl]
		for iq in rshifts:
			xx = iq[0]*shrink
			yy = iq[1]*shrink
			dss = sp_fundamentals.fshift(ds, xx, yy)
			dss.set_attr("is_complex",0)
			"""Multiline Comment9"""
			#MULTILINEMULTILINEMULTILINE 9
				#MULTILINEMULTILINEMULTILINE 9
				#MULTILINEMULTILINEMULTILINE 9
				#MULTILINEMULTILINEMULTILINE 9
			#MULTILINEMULTILINEMULTILINE 9
				#MULTILINEMULTILINEMULTILINE 9
				#MULTILINEMULTILINEMULTILINE 9
			#MULTILINEMULTILINEMULTILINE 9
			data[kl].append(dss)
		data[kl][0].set_attr("particle_group",particle_group)  #  Pass group number only in the first shifted copy.
		del projdata[kl]
	return data, ctfs, bckgnoise





















































































def steptwo(tvol, tweight, treg, cfsc = None, regularized = True):
	global Tracker, Blockdata
	nz = tweight.get_zsize()
	ny = tweight.get_ysize()
	nx = tweight.get_xsize()
	tvol.set_attr("is_complex",1)
	if regularized:
		nr = len(cfsc)
		limitres = 0
		for i in range(nr):
			cfsc[i] = min(max(cfsc[i], 0.0), 0.999)
			#print( i,cfsc[i] )
			if( cfsc[i] == 0.0 ):
				limitres = i-1
				break
		if( limitres == 0 ): limitres = nr-2;
		ovol = sp_utilities.reshape_1d(cfsc, nr, 2*nr)
		limitres = 2*min(limitres, Tracker["maxfrad"])  # 2 on account of padding, which is always on
		maxr2 = limitres**2
		for i in range(limitres+1, len(ovol), 1):   ovol[i] = 0.0
		ovol[0] = 1.0
		#print(" ovol  ", ovol)
		it = sp_utilities.model_blank(2*nr)
		for i in range(2*nr):  it[i] = ovol[i]
		del ovol
		#  Do not regularize first four
		for i in range(5):  treg[i] = 0.0
		EMAN2_cppwrap.Util.reg_weights(tweight, treg, it)
		del it
	else:
		limitres = 2*min(Tracker["constants"]["nnxo"]//2, Tracker["maxfrad"])
		maxr2 = limitres**2
	#  Iterative weights
	if( Tracker["constants"]["symmetry"] != "c1" ):
		tvol    = tvol.symfvol(Tracker["constants"]["symmetry"], limitres)
		tweight = tweight.symfvol(Tracker["constants"]["symmetry"], limitres)

	#  tvol is overwritten, meaning it is also an output
	EMAN2_cppwrap.Util.iterefa(tvol, tweight, maxr2, Tracker["constants"]["nnxo"])
	#  Either pad or window in F space to 2*nnxo
	nx = tvol.get_ysize()
	if( nx > 2*Tracker["constants"]["nnxo"] ):
		tvol = sp_fundamentals.fdecimate(tvol, 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], False, False)
	elif(nx < 2*Tracker["constants"]["nnxo"] ):
		tvol = sp_fundamentals.fpol(tvol, 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], RetReal = False, normalize = False)

	tvol = sp_fundamentals.fft(tvol)
	tvol = sp_fundamentals.cyclic_shift(tvol,Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	tvol.set_attr("npad",2)
	tvol.div_sinc(1)
	tvol.del_attr("npad")
	tvol = EMAN2_cppwrap.Util.window(tvol, Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	tvol = sp_morphology.cosinemask(tvol,Tracker["constants"]["nnxo"]//2-1-5,5, None) # clean artifacts in corners

	return tvol












































































































































































































































































































def Xali3D_direct_ccc(data, refang, shifts, ctfs = None, bckgnoise = None, kb3D = None):
	global Tracker, Blockdata
	pass#IMPORTIMPORTIMPORT from sp_projection 	import prgs,prgl
	pass#IMPORTIMPORTIMPORT from sp_fundamentals 	import fft
	pass#IMPORTIMPORTIMPORT from sp_utilities 		import wrap_mpi_gatherv
	pass#IMPORTIMPORTIMPORT from math 			import sqrt
	#  Input data has to be CTF-multiplied, preshifted
	#  Output - newpar, see structure
	#    newpar = [[i, [worst_similarity, sum_all_similarities], [[-1, -1.0e23] for j in range(Tracker["lentop"])]] for i in range(len(data))]
	#    newpar = [[params],[],... len(data)]
	#    params = [particleID, [worst_similarity, sum_all_similarities],[imageallparams]]]
	#    imageallparams = [[orientation, similarity],[],...  number of all orientations ]
	#  Coding of orientations:
	#    hash = ang*100000000 + lpsi*1000 + ishift
	#    ishift = hash%1000
	#    ipsi = (hash/1000)%100000
	#    iang  = hash/100000000
	#  To get best matching for particle #kl:
	#     hash_best = newpar[kl][-1][0][0]
	#     best_sim  = newpar[kl][-1][0][1]
	#  To sort:
	pass#IMPORTIMPORTIMPORT from operator 		import itemgetter#, attrgetter, methodcaller
	pass#IMPORTIMPORTIMPORT from math 			import exp
	
	#   params.sort(key=itemgetter(2))
	at = time.time()
	if(Blockdata["myid"] == 0):  sp_global_def.sxprint("  ENTERING Xali buffered exhaustive CCC  ")
	npsi = int(360./Tracker["delta"])
	nang = len(refang)
	ndat = len(data)

	ny = data[0][0].get_ysize()
	mask = EMAN2_cppwrap.Util.unrollmask(ny,ny)
	nxt = 2*(mask.get_xsize())

	"""Multiline Comment13"""
	#MULTILINEMULTILINEMULTILINE 13
		#MULTILINEMULTILINEMULTILINE 13
			#MULTILINEMULTILINEMULTILINE 13
	#MULTILINEMULTILINEMULTILINE 13

	if Tracker["mainiteration"]>1 :
		#first = True
		if Tracker["constants"]["CTF"] :
			for kl in range(ndat):
				for im in range(len(shifts)):
					EMAN2_cppwrap.Util.mulclreal(data[kl][im], ctfs[kl])
		del ctfs
		if bckgnoise:  #  This should be a flag to activate sharpening during refinement as bckgnoise is always present (for 3D later)
			for kl in range(ndat):
				temp = EMAN2_cppwrap.Util.unroll1dpw(ny,ny, bckgnoise[kl])
				for im in range(len(shifts)):
					EMAN2_cppwrap.Util.mulclreal(data[kl][im], temp)
			del bckgnoise
	#else:  first = False

	disp_unit = numpy.dtype("f4").itemsize

	#  REFVOL
	if( Blockdata["myid_on_node"] == 0 ):
		if( Blockdata["myid"] == Blockdata["main_node"] ):
			odo = get_refvol( Tracker["nxpolar"] )
			nxvol = odo.get_xsize()
			nyvol = odo.get_ysize()
			nzvol = odo.get_zsize()
		else:
			nxvol = 0
			nyvol = 0
			nzvol = 0

		nxvol = sp_utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nyvol = sp_utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nzvol = sp_utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

		if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = sp_utilities.model_blank( nxvol,nyvol, nzvol)

		sp_utilities.bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

		odo = sp_projection.prep_vol( odo, npad = 2, interpolation_method = 1 )
		ndo = EMAN2_cppwrap.EMNumPy.em2numpy(odo)
		nxvol = odo.get_xsize()
		nyvol = odo.get_ysize()
		nzvol = odo.get_zsize()
		orgsizevol = nxvol*nyvol*nzvol
		sizevol = orgsizevol
	else:
		orgsizevol = 0
		sizevol = 0
		nxvol = 0
		nyvol = 0
		nzvol = 0

	orgsizevol = sp_utilities.bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
	nxvol = sp_utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
	nyvol = sp_utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
	nzvol = sp_utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

	win_vol, base_vol  = mpi.mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
	sizevol = orgsizevol
	if( Blockdata["myid_on_node"] != 0 ):
		base_vol, = mpi.mpi_win_shared_query(win_vol, mpi.MPI_PROC_NULL)

	volbuf = numpy.frombuffer(numpy.core.multiarray.int_asbuffer(base_vol, sizevol*disp_unit), dtype = 'f4')
	volbuf = volbuf.reshape(nzvol, nyvol, nxvol)
	if( Blockdata["myid_on_node"] == 0 ):
		numpy.copyto(volbuf,ndo)
		del odo,ndo


	#volprep = EMNumPy.assign_numpy_to_emdata(volbuf)
	emnumpy1 = EMAN2_cppwrap.EMNumPy()
	volprep = emnumpy1.register_numpy_to_emdata(volbuf)
	volprep.set_attr_dict({'is_complex':1,  'is_complex_x': 0, 'is_fftodd': 0, 'is_fftpad': 1, 'is_shuffled': 1,'npad': 2})


	#  BIG BUFFER
	size_of_one_image = ny*nxt
	lenbigbuf = min(Blockdata["no_of_processes_per_group"],nang)*npsi
	orgsize = lenbigbuf*size_of_one_image #  This is number of projections to be computed simultaneously times their size

	if( Blockdata["myid_on_node"] == 0 ): size = orgsize
	else:  size = 0

	win_sm, base_ptr  = mpi.mpi_win_allocate_shared( size*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
	size = orgsize
	if( Blockdata["myid_on_node"] != 0 ):
		base_ptr, = mpi.mpi_win_shared_query(win_sm, mpi.MPI_PROC_NULL)

	buffer = numpy.frombuffer(numpy.core.multiarray.int_asbuffer(base_ptr, size*disp_unit), dtype = 'f4')
	buffer = buffer.reshape(lenbigbuf, ny, nxt)
	#ncbuf = lenbigbuf//2
	
	#bigbuffer = EMNumPy.assign_numpy_to_emdata(buffer)

	emnumpy2 = EMAN2_cppwrap.EMNumPy()
	bigbuffer = emnumpy2.register_numpy_to_emdata(buffer)

	emnumpy3 = EMAN2_cppwrap.EMNumPy()

	#  end of setup

	at = time.time()
	#  Here we simply search for a max
	newpar = [[i, [1.0], [[-1,-1.0e23]] ] for i in range(ndat)]

	for i in range(nang):
		if( ( Blockdata["myid"] == Blockdata["main_node"])  and  (i%(max(1,nang/5)) == 0) and (i>0)):
			sp_global_def.sxprint( "  Angle :%7d   %5d  %5.1f"%(i,ndat,float(i)/float(nang)*100.) + "%" +"   %10.1fmin"%((time.time()-at)/60.))

		if(i%Blockdata["no_of_processes_per_group"] == 0 ):  #  Time to fill up the buffer
			for itemp in range(i, min(i+Blockdata["no_of_processes_per_group"], nang)):
				if( itemp-i == Blockdata["myid_on_node"]):
					for j in range(npsi):
						psi = (refang[i][2] + j*Tracker["delta"])%360.0
						###if kb3D:  rtemp = fft(prgs(volprep, kb3D, [refang[i][0],refang[i][1],psi, 0.0,0.0]))
						###else:     
						temp = sp_projection.prgl(volprep,[ refang[itemp][0],refang[itemp][1],psi, 0.0,0.0], 1, False)
						EMAN2_cppwrap.Util.mulclreal(temp, mask)
						nrmref = numpy.sqrt(EMAN2_cppwrap.Util.innerproduct(temp, temp, None))
						temp.set_attr("is_complex",0)
						EMAN2_cppwrap.Util.mul_scalar(temp, 1.0/nrmref)
						bigbuffer.insert_clip(temp,(0,0,(itemp-i)*npsi+j))
	
			mpi.mpi_barrier(Blockdata["shared_comm"])

		iang = i*100000000
		for j in range(npsi):
			iangpsi = j*1000 + iang
			psi = (refang[i][2] + j*Tracker["delta"])%360.0
			#temp = Util.window(bigbuffer, nxt, ny, 1, 0, 0, -ncbuf + (i%Blockdata["no_of_processes_per_group"])*npsi + j)

			#  Here we get an image from a buffer by assigning an address instead of copy.
			pointer_location = base_ptr + ((i%Blockdata["no_of_processes_per_group"])*npsi + j)*size_of_one_image*disp_unit
			img_buffer = numpy.frombuffer(numpy.core.multiarray.int_asbuffer(pointer_location, size_of_one_image*disp_unit), dtype = 'f4')
			img_buffer = img_buffer.reshape(ny, nxt)
			#temp = EMNumPy.numpy2em(img_buffer)
			temp = emnumpy3.register_numpy_to_emdata(img_buffer)
			temp.set_attr("is_complex",1)

			#temp *= (1000.0/nrmref)
			#nrmref = 1000.
			for kl,emimage in enumerate(data):
				for im in range(len(shifts)):
					peak = EMAN2_cppwrap.Util.innerproduct(temp, emimage[im], None)
					if(peak>newpar[kl][2][0][1]):  newpar[kl][2] = [[im + iangpsi, peak]]

	#print  " >>>  %4d   %12.3e       %12.5f     %12.5f     %12.5f     %12.5f     %12.5f"%(best,simis[0],newpar[0][0],newpar[0][1],newpar[0][2],newpar[0][3],newpar[0][4])
	###if Blockdata["myid"] == Blockdata["main_node"]:  print "  Finished :",time()-at
	
	#print("  SEARCHES DONE  ",Blockdata["myid"])

	#mpi_barrier(MPI_COMM_WORLD)
	mpi.mpi_win_free(win_vol)
	mpi.mpi_win_free(win_sm)
	emnumpy1.unregister_numpy_from_emdata()
	emnumpy2.unregister_numpy_from_emdata()
	emnumpy3.unregister_numpy_from_emdata()
	del emnumpy1, emnumpy2, emnumpy3

	mpi.mpi_barrier(Blockdata["shared_comm"])

	#print("  NORMALIZATION DONE  ",Blockdata["myid"])
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	#print("  ALL DONE  ",Blockdata["myid"])
	return newpar

def XXali3D_direct_ccc(data, refang, shifts, coarse_angles, coarse_shifts, ctfs = None, bckgnoise = None, kb3D = None):
	global Tracker, Blockdata
	pass#IMPORTIMPORTIMPORT from sp_projection 	import prgs,prgl
	pass#IMPORTIMPORTIMPORT from sp_fundamentals 	import fft
	pass#IMPORTIMPORTIMPORT from sp_utilities 		import wrap_mpi_gatherv
	pass#IMPORTIMPORTIMPORT from math 			import sqrt
	#  Input data has to be CTF-multiplied, preshifted
	#  Output - newpar, see structure
	#    newpar = [[i, [worst_similarity, sum_all_similarities], [[-1, -1.0e23] for j in range(Tracker["lentop"])]] for i in range(len(data))]
	#    newpar = [[params],[],... len(data)]
	#    params = [particleID, [worst_similarity, sum_all_similarities],[imageallparams]]]
	#    imageallparams = [[orientation, similarity],[],...  number of all orientations ]
	#  Coding of orientations:
	#    hash = ang*100000000 + lpsi*1000 + ishift
	#    ishift = hash%1000
	#    ipsi = (hash/1000)%100000
	#    iang  = hash/100000000
	#  To get best matching for particle #kl:
	#     hash_best = newpar[kl][-1][0][0]
	#     best_sim  = newpar[kl][-1][0][1]
	#  To sort:
	pass#IMPORTIMPORTIMPORT from operator 		import itemgetter#, attrgetter, methodcaller
	pass#IMPORTIMPORTIMPORT from math 			import exp
	
	#   params.sort(key=itemgetter(2))
	at = time.time()
	if(Blockdata["myid"] == 0):
		print_dict(Tracker,"PROJECTION MATCHING parameters of buffered exhaustive CCC")
		sp_utilities.write_text_row(coarse_angles,"coarse_angles.txt")

	npsi = int(360./Tracker["delta"])
	nang = len(refang)
	ndat = len(data)

	#  FINE SEARCH CONSTANTS
	#  Fine grids are shifted by half-fine_step
	nang = len(refang)
	npsi = int(360./Tracker["delta"])
	nshifts = len(shifts)
	n_fine_shifts = 4

	#  COARSE SEARCH CONSTANTS
	n_coarse_ang = len(coarse_angles)
	coarse_delta = 2*Tracker["delta"]
	n_coarse_psi = int(360./coarse_delta)
	n_coarse_shifts = len(coarse_shifts)

	#  each data holds two lists: n_coarse_shifts of coarse shifted versions and then nhifts versions of Fine shifted images

	ny = data[0][0].get_ysize()
	mask = EMAN2_cppwrap.Util.unrollmask(ny,ny)
	nxt = 2*(mask.get_xsize())

	"""Multiline Comment14"""
	#MULTILINEMULTILINEMULTILINE 14
		#MULTILINEMULTILINEMULTILINE 14
			#MULTILINEMULTILINEMULTILINE 14
	#MULTILINEMULTILINEMULTILINE 14

	if Tracker["mainiteration"]>1 :
		#first = True
		if Tracker["constants"]["CTF"] :
			for kl in range(ndat):
				for im in range(n_coarse_shifts+nshifts):
					EMAN2_cppwrap.Util.mulclreal(data[kl][im], ctfs[kl])
		del ctfs
		if bckgnoise:  #  This should be a flag to activate sharpening during refinement as bckgnoise is always present (for 3D later)
			for kl in range(ndat):
				temp = EMAN2_cppwrap.Util.unroll1dpw(ny,ny, bckgnoise[kl])
				for im in range(n_coarse_shifts+nshifts):
					EMAN2_cppwrap.Util.mulclreal(data[kl][im], temp)
			del bckgnoise
	#else:  first = False

	disp_unit = numpy.dtype("f4").itemsize

	#  REFVOL
	if( Blockdata["myid_on_node"] == 0 ):
		if( Blockdata["myid"] == Blockdata["main_node"] ):
			odo = get_refvol( Tracker["nxpolar"] )
			nxvol = odo.get_xsize()
			nyvol = odo.get_ysize()
			nzvol = odo.get_zsize()
		else:
			nxvol = 0
			nyvol = 0
			nzvol = 0

		nxvol = sp_utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nyvol = sp_utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nzvol = sp_utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

		if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = sp_utilities.model_blank( nxvol,nyvol, nzvol)

		sp_utilities.bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

		odo = sp_projection.prep_vol( odo, npad = 2, interpolation_method = 1 )
		ndo = EMAN2_cppwrap.EMNumPy.em2numpy(odo)
		nxvol = odo.get_xsize()
		nyvol = odo.get_ysize()
		nzvol = odo.get_zsize()
		orgsizevol = nxvol*nyvol*nzvol
		sizevol = orgsizevol
	else:
		orgsizevol = 0
		sizevol = 0
		nxvol = 0
		nyvol = 0
		nzvol = 0

	orgsizevol = sp_utilities.bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
	nxvol = sp_utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
	nyvol = sp_utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
	nzvol = sp_utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

	win_vol, base_vol  = mpi.mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
	sizevol = orgsizevol
	if( Blockdata["myid_on_node"] != 0 ):
		base_vol, = mpi.mpi_win_shared_query(win_vol, mpi.MPI_PROC_NULL)

	volbuf = numpy.frombuffer(numpy.core.multiarray.int_asbuffer(base_vol, sizevol*disp_unit), dtype = 'f4')
	volbuf = volbuf.reshape(nzvol, nyvol, nxvol)
	if( Blockdata["myid_on_node"] == 0 ):
		numpy.copyto(volbuf,ndo)
		del odo,ndo


	#volprep = EMNumPy.assign_numpy_to_emdata(volbuf)
	emnumpy1 = EMAN2_cppwrap.EMNumPy()
	volprep = emnumpy1.register_numpy_to_emdata(volbuf)
	volprep.set_attr_dict({'is_complex':1,  'is_complex_x': 0, 'is_fftodd': 0, 'is_fftpad': 1, 'is_shuffled': 1,'npad': 2})


	#  BIG BUFFER
	size_of_one_image = ny*nxt
	#lenbigbuf = min(Blockdata["no_of_processes_per_group"],n_coarse_ang)*n_coarse_psi
	lenbigbuf = nang*npsi
	orgsize = lenbigbuf*size_of_one_image #  This is number of projections to be computed simultaneously times their size
	if( Blockdata["myid"] == 0 ):  sp_global_def.sxprint("  BIGBUFFER  ",float(orgsize)/1.e9,"GB")
	if( Blockdata["myid_on_node"] == 0 ): size = orgsize
	else:  size = 0

	win_sm, base_ptr  = mpi.mpi_win_allocate_shared( size*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
	size = orgsize
	if( Blockdata["myid_on_node"] != 0 ):
		base_ptr, = mpi.mpi_win_shared_query(win_sm, mpi.MPI_PROC_NULL)

	buffer = numpy.frombuffer(numpy.core.multiarray.int_asbuffer(base_ptr, size*disp_unit), dtype = 'f4')
	buffer = buffer.reshape(lenbigbuf, ny, nxt)
	#ncbuf = lenbigbuf//2
	
	#bigbuffer = EMNumPy.assign_numpy_to_emdata(buffer)

	emnumpy2 = EMAN2_cppwrap.EMNumPy()
	bigbuffer = emnumpy2.register_numpy_to_emdata(buffer)

	emnumpy3 = EMAN2_cppwrap.EMNumPy()

	#  end of setup

	at = time.time()
	#  Here we simply search for a max
	newpar = [[i, [1.0], [[-1,-1.0e23]] ] for i in range(ndat)]

	for i in range(n_coarse_ang):
		if( ( Blockdata["myid"] == Blockdata["main_node"])  and  (i%(max(1,n_coarse_ang/5)) == 0) and (i>0)):
			sp_global_def.sxprint( "  Angle :%7d   %5d  %5.1f"%(i,ndat,float(i)/float(n_coarse_ang)*100.) + "%" +"   %10.1fmin"%((time.time()-at)/60.))

		if(i%Blockdata["no_of_processes_per_group"] == 0 ):  #  Time to fill up the buffer
			for itemp in range(i, min(i+Blockdata["no_of_processes_per_group"], n_coarse_ang)):
				if( itemp-i == Blockdata["myid_on_node"]):
					for j in range(n_coarse_psi):
						psi = (coarse_angles[i][2] + j*coarse_delta)%360.0
						###if kb3D:  rtemp = fft(prgs(volprep, kb3D, [refang[i][0],refang[i][1],psi, 0.0,0.0]))
						###else:
						temp = sp_projection.prgl(volprep,[ coarse_angles[itemp][0],coarse_angles[itemp][1],psi, 0.0,0.0], 1, False)
						EMAN2_cppwrap.Util.mulclreal(temp, mask)
						nrmref = numpy.sqrt(EMAN2_cppwrap.Util.innerproduct(temp, temp, None))
						temp.set_attr("is_complex",0)
						EMAN2_cppwrap.Util.mul_scalar(temp, 1.0/nrmref)
						bigbuffer.insert_clip(temp,(0,0,(itemp-i)*npsi+j))

			mpi.mpi_barrier(Blockdata["shared_comm"])

		iang = i*100000000
		for j in range(n_coarse_psi):
			iangpsi = j*1000 + iang
			psi = (coarse_angles[i][2] + j*coarse_delta)%360.0
			#temp = Util.window(bigbuffer, nxt, ny, 1, 0, 0, -ncbuf + (i%Blockdata["no_of_processes_per_group"])*npsi + j)

			#  Here we get an image from a buffer by assigning an address instead of copy.
			pointer_location = base_ptr + ((i%Blockdata["no_of_processes_per_group"])*n_coarse_psi + j)*size_of_one_image*disp_unit
			img_buffer = numpy.frombuffer(numpy.core.multiarray.int_asbuffer(pointer_location, size_of_one_image*disp_unit), dtype = 'f4')
			img_buffer = img_buffer.reshape(ny, nxt)
			#temp = EMNumPy.numpy2em(img_buffer)
			temp = emnumpy3.register_numpy_to_emdata(img_buffer)
			temp.set_attr("is_complex",1)

			#temp *= (1000.0/nrmref)
			#nrmref = 1000.
			for kl,emimage in enumerate(data):
				for im in range(n_coarse_shifts):
					peak = EMAN2_cppwrap.Util.innerproduct(temp, emimage[im], None)
					if(peak>newpar[kl][2][0][1]):  newpar[kl][2] = [[im + iangpsi, peak]]

	#print  " >>>  %4d   %12.3e       %12.5f     %12.5f     %12.5f     %12.5f     %12.5f"%(best,simis[0],newpar[0][0],newpar[0][1],newpar[0][2],newpar[0][3],newpar[0][4])
	###if Blockdata["myid"] == Blockdata["main_node"]:  print "  Finished :",time()-at
	if Blockdata["myid"] == Blockdata["main_node"]:   sp_global_def.sxprint("  COARSE SEARCHES DONE  ","   %10.1fmin"%((time.time()-at)/60.))
	at = time.time()
	for itemp in range(0, min(Blockdata["no_of_processes_per_group"], nang)):
		if( itemp == Blockdata["myid_on_node"]):
			for j in range(npsi):
				psi = (refang[i][2] + j*Tracker["delta"])%360.0
				temp = sp_projection.prgl(volprep,[ refang[itemp][0],refang[itemp][1],psi, 0.0,0.0], 1, False)
				EMAN2_cppwrap.Util.mulclreal(temp, mask)
				nrmref = numpy.sqrt(EMAN2_cppwrap.Util.innerproduct(temp, temp, None))
				temp.set_attr("is_complex",0)
				EMAN2_cppwrap.Util.mul_scalar(temp, 1.0/nrmref)
				bigbuffer.insert_clip(temp,(0,0,itemp*npsi+j))

	mpi.mpi_barrier(Blockdata["shared_comm"])
	if( Blockdata["myid"] == Blockdata["main_node"] ):   sp_global_def.sxprint("  REFPROJ DONE  ","   %10.1fmin"%((time.time()-at)/60.))
	at = time.time()

	opar = []
	Blockdata["angle_set"] = refang
	#  Extract normals from rotation matrices
	refdirs = sp_utilities.angles_to_normals(refang)

	for kl,emimage in enumerate(data):
		hashparams = newpar[kl][2][0][0]
		ipsiandiang	= hashparams/1000
		oldiang = ipsiandiang/100000
		ipsi = ipsiandiang%100000
		ishift = hashparams%1000
		tshifts = get_shifts_neighbors(shifts, coarse_shifts[ishift])
		ltabang = find_nearest_k_refangles_to_many_angles(refdirs, [coarse_angles[oldiang]], Tracker["delta"], howmany = 4)
		opar.append([coarse_angles[oldiang][0],coarse_angles[oldiang][1],(coarse_angles[oldiang][2]+ipsi*coarse_delta)%360.0 , coarse_shifts[ishift][0],coarse_shifts[ishift][1]])
		for i2 in range(4):
			iang = ltabang[0][i2]
			for i3 in range(2):  # psi
				itpsi = int((coarse_angles[oldiang][2] + ipsi*coarse_delta - refang[iang][2]+360.0)/Tracker["delta"])
				itpsi = (itpsi + i3)%npsi
				pointer_location = base_ptr + (iang*n_coarse_psi + itpsi)*size_of_one_image*disp_unit
				img_buffer = numpy.frombuffer(numpy.core.multiarray.int_asbuffer(pointer_location, size_of_one_image*disp_unit), dtype = 'f4')
				img_buffer = img_buffer.reshape(ny, nxt)
				#temp = EMNumPy.numpy2em(img_buffer)
				temp = emnumpy3.register_numpy_to_emdata(img_buffer)
				temp.set_attr("is_complex",1)

				newpar[kl][2][0][1] = -1.0e23
				for i4 in range(len(tshifts)):
					peak = EMAN2_cppwrap.Util.innerproduct(temp, emimage[im+n_coarse_shifts], None) # take now fine images
					if(peak>newpar[kl][2][0][1]):  newpar[kl][2] = [[iang*100000000 + itpsi*1000 + tshifts[i4], peak]]

	#mpi_barrier(MPI_COMM_WORLD)
	mpi.mpi_win_free(win_vol)
	mpi.mpi_win_free(win_sm)
	emnumpy1.unregister_numpy_from_emdata()
	emnumpy2.unregister_numpy_from_emdata()
	emnumpy3.unregister_numpy_from_emdata()
	del emnumpy1, emnumpy2, emnumpy3

	mpi.mpi_barrier(Blockdata["shared_comm"])
	if Blockdata["myid"] == Blockdata["main_node"]:   sp_global_def.sxprint("  FINE SEARCH DONE  ","   %10.1fmin"%((time.time()-at)/60.))

	#print("  NORMALIZATION DONE  ",Blockdata["myid"])
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	opar = sp_utilities.wrap_mpi_gatherv(opar, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	if( Blockdata["myid"] == Blockdata["main_node"] ):   sp_utilities.write_text_row(opar,"opar.txt")
	#print("  ALL DONE  ",Blockdata["myid"])
	return newpar




























































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































def ali3D_local_polar_ccc(refang, shifts, coarse_angles, coarse_shifts, procid, original_data = None, oldparams = None, \
					preshift = False, apply_mask = True, nonorm = False, applyctf = True):
	global Tracker, Blockdata
	pass#IMPORTIMPORTIMPORT from sp_projection 	import prgs,prgl
	pass#IMPORTIMPORTIMPORT from sp_fundamentals 	import fft
	pass#IMPORTIMPORTIMPORT from sp_utilities 		import wrap_mpi_gatherv
	pass#IMPORTIMPORTIMPORT from math 			import sqrt
	#  Input data has to be CTF-multiplied, preshifted
	#  Output - newpar, see structure
	#    newpar = [[i, [worst_similarity, sum_all_similarities], [[-1, -1.0e23] for j in range(Tracker["lentop"])]] for i in range(len(data))]
	#    newpar = [[params],[],... len(data)]
	#    params = [particleID, [worst_similarity, sum_all_similarities],[imageallparams]]]
	#    imageallparams = [[orientation, similarity],[],...  number of all orientations ]
	#  Coding of orientations:
	#    hash = ang*100000000 + lpsi*1000 + ishift
	#    ishift = hash%1000
	#    ipsi = (hash/1000)%100000
	#    iang  = hash/100000000
	#  To get best matching for particle #kl:
	#     hash_best = newpar[kl][-1][0][0]
	#     best_sim  = newpar[kl][-1][0][1]
	#  To sort:
	#from operator 		import itemgetter#, attrgetter, methodcaller
	#   params.sort(key=itemgetter(2))
	#from fundamentals import resample
	pass#IMPORTIMPORTIMPORT from sp_utilities    import get_im, model_gauss_noise, set_params_proj, get_params_proj
	pass#IMPORTIMPORTIMPORT from sp_fundamentals import fdecimate, fshift, fft
	pass#IMPORTIMPORTIMPORT from sp_filter       import filt_ctf, filt_table
	pass#IMPORTIMPORTIMPORT from sp_applications import MPI_start_end
	pass#IMPORTIMPORTIMPORT from math         import sqrt

	if(Blockdata["myid"] == Blockdata["main_node"]):
		print_dict(Tracker,"PROJECTION MATCHING parameters of buffered local polar ccc")

	at = time.time()
	shrinkage = float(Tracker["nxpolar"])/float(Tracker["constants"]["nnxo"])
	shrink = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
	mode = "F"
	numr = Numrinit_local(1, radius, 1, mode)
	wr = sp_alignment.ringwe(numr, mode)
	cnx = float(Tracker["nxpolar"]//2 + 1)

	##if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  RADIUS IN POLAR  ",radius,numr[-1])
	#  FINE SEARCH CONSTANTS
	nang = len(refang)
	ac_fine = numpy.cos(numpy.radians(Tracker["an"]))
	npsi = int(360./Tracker["delta"])
	mpsi = 2
	c_fine_psi = mpsi//2
	nshifts = len(shifts)
	n_fine_shifts = 4

	nima = len(original_data)
	mask = EMAN2_cppwrap.Util.unrollmask(Tracker["nxinit"],Tracker["nxinit"])
	for j in range(Tracker["nxinit"]//2,Tracker["nxinit"]):  mask[0,j]=1.0
	mask2D = sp_utilities.model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	reachpw = mask.get_xsize()  # The last element of accumulated pw is zero so for the full size nothing is added.

	#  COARSE SEARCH CONSTANTS
	n_coarse_ang = len(coarse_angles)
	coarse_delta = 2*Tracker["delta"]
	n_coarse_psi = int(360./coarse_delta)
	###m_coarse_psi = int((2*2*Tracker["an"])/coarse_delta + 0.5) + 1
	m_coarse_psi = int(Tracker["an"]/Tracker["delta"] + 0.5)
	c_coarse_psi = m_coarse_psi//2
	n_coarse_shifts = len(coarse_shifts)

	coarse_shifts_shrank = [None]*n_coarse_shifts
	for ib in range(n_coarse_shifts):
		coarse_shifts_shrank[ib] = [coarse_shifts[ib][0]*shrinkage,coarse_shifts[ib][1]*shrinkage]

	###if(Blockdata["myid"] == Blockdata["main_node"]): sxprint("   TRETR  ",Tracker["constants"]["nnxo"],Tracker["nxinit"],reachpw,n_coarse_ang,coarse_delta,n_coarse_psi,m_coarse_psi,c_coarse_psi,n_coarse_shifts)

	#if(Blockdata["myid"] == Blockdata["main_node"]):
	#	sxprint( original_data[0].get_attr("identifier") )
	#	sxprint( original_data[1].get_attr("identifier") )
	
	disp_unit = numpy.dtype("f4").itemsize

	#  REFVOL
	if( Blockdata["myid_on_node"] == 0 ):
		if( Blockdata["myid"] == Blockdata["main_node"] ):
			odo = get_refvol( Tracker["nxpolar"] )
			nxvol = odo.get_xsize()
			nyvol = odo.get_ysize()
			nzvol = odo.get_zsize()
		else:
			nxvol = 0
			nyvol = 0
			nzvol = 0

		nxvol = sp_utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nyvol = sp_utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nzvol = sp_utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

		if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = sp_utilities.model_blank( nxvol,nyvol, nzvol)

		sp_utilities.bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

		odo = sp_projection.prep_vol( odo, npad = 2, interpolation_method = 1 )
		ndo = EMAN2_cppwrap.EMNumPy.em2numpy(odo)
		nxvol = odo.get_xsize()
		nyvol = odo.get_ysize()
		nzvol = odo.get_zsize()
		orgsizevol = nxvol*nyvol*nzvol
		sizevol = orgsizevol
	else:
		orgsizevol = 0
		sizevol = 0
		nxvol = 0
		nyvol = 0
		nzvol = 0

	orgsizevol = sp_utilities.bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
	nxvol = sp_utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
	nyvol = sp_utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
	nzvol = sp_utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

	win_vol, base_vol  = mpi.mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
	sizevol = orgsizevol
	if( Blockdata["myid_on_node"] != 0 ):
		base_vol, = mpi.mpi_win_shared_query(win_vol, mpi.MPI_PROC_NULL)

	volbuf = numpy.frombuffer(numpy.core.multiarray.int_asbuffer(base_vol, sizevol*disp_unit), dtype = 'f4')
	volbuf = volbuf.reshape(nzvol, nyvol, nxvol)
	if( Blockdata["myid_on_node"] == 0 ):
		numpy.copyto(volbuf,ndo)
		del odo,ndo

	#volprep = EMNumPy.assign_numpy_to_emdata(volbuf)
	emnumpy1 = EMAN2_cppwrap.EMNumPy()
	volprep = emnumpy1.register_numpy_to_emdata(volbuf)
	volprep.set_attr_dict({'is_complex':1,  'is_complex_x': 0, 'is_fftodd': 0, 'is_fftpad': 1, 'is_shuffled': 1,'npad': 2})

	if( Tracker["nxinit"] != Tracker["nxpolar"] ):
		mpi.mpi_win_free(win_vol)
		#  REFVOL FOR ML
		if( Blockdata["myid_on_node"] == 0 ):
			if( Blockdata["myid"] == Blockdata["main_node"] ):
				odo = get_refvol( Tracker["nxpolar"] )
				nxvol = odo.get_xsize()
				nyvol = odo.get_ysize()
				nzvol = odo.get_zsize()
			else:
				nxvol = 0
				nyvol = 0
				nzvol = 0

			nxvol = sp_utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
			nyvol = sp_utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
			nzvol = sp_utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

			if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = sp_utilities.model_blank( nxvol,nyvol, nzvol)

			sp_utilities.bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

			odo = sp_projection.prep_vol( odo, npad = 2, interpolation_method = 1 )
			ndo = EMAN2_cppwrap.EMNumPy.em2numpy(odo)
			nxvol = odo.get_xsize()
			nyvol = odo.get_ysize()
			nzvol = odo.get_zsize()
			orgsizevol = nxvol*nyvol*nzvol
			sizevol = orgsizevol
		else:
			orgsizevol = 0
			sizevol = 0
			nxvol = 0
			nyvol = 0
			nzvol = 0

		orgsizevol = sp_utilities.bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
		nxvol = sp_utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
		nyvol = sp_utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
		nzvol = sp_utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

		win_volinit, base_volinit  = mpi.mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
		sizevol = orgsizevol
		if( Blockdata["myid_on_node"] != 0 ):
			base_volinit, = mpi.mpi_win_shared_query(win_volinit, mpi.MPI_PROC_NULL)

		volbufinit = numpy.frombuffer(numpy.core.multiarray.int_asbuffer(base_volinit, sizevol*disp_unit), dtype = 'f4')
		volbufinit = volbufinit.reshape(nzvol, nyvol, nxvol)
		if( Blockdata["myid_on_node"] == 0 ):
			numpy.copyto(volbufinit,ndoinit)
			del odo,ndoinit

		#volinit = EMNumPy.assign_numpy_to_emdata(volbufinit)
		emnumpy4 = EMAN2_cppwrap.EMNumPy()
		volinit = emnumpy4.register_numpy_to_emdata(volbufinit)
		volinit.set_attr_dict({'is_complex':1,  'is_complex_x': 0, 'is_fftodd': 0, 'is_fftpad': 1, 'is_shuffled': 1,'npad': 2})
		if( Blockdata["myid_on_node"] == 0 ):  volinit.update()
		mpi.mpi_barrier(Blockdata["shared_comm"])
		###if( Blockdata["myid"] < 10 ):
		###	sxprint(" VOLPREP  ",volinit[0],volinit[1],info(volinit))
	else: volinit = volprep
	#  End of replaced volprep

    ### local search is plugged in from here
	#  START CONES
	#  This has to be systematically done per node
	#
	crefim = EMAN2_cppwrap.Util.Polar2Dm(sp_utilities.model_blank(Tracker["nxpolar"],Tracker["nxpolar"]), cnx, cnx, numr, mode)
	size_of_one_image = crefim.get_xsize()
	#  We will assume half of the memory is available.  We will do it betteer later.
	numberofrefs_inmem = int(Tracker["constants"]["memory_per_node"]/4/((size_of_one_image*disp_unit)/1.0e9))
	####if( Blockdata["myid_on_node"] == 0  ):  sxprint( " MEMEST ", n_coarse_ang,numberofrefs_inmem)
	#  number of references that will fit into one mode
	normals_set = sp_utilities.angles_to_normals(coarse_angles)
	Blockdata["angle_set"] = coarse_angles
	if( n_coarse_ang <= numberofrefs_inmem ):
		number_of_cones = 1
		numberofrefs_inmem = n_coarse_ang
		assignments_to_cones = [list(range(len(oldparams)))]

		assignments_of_refangles_to_cones = [[] for i in range(len(assignments_to_cones))]
		assignments_of_refangles_to_angles = [[] for i in range(nima)]  # for each myid separately, these are angles on this myid

		for i,q in enumerate(assignments_to_cones):
			#  find assignments of refdirs to this cone within an of each angledirs, so we need symmetry neighbors set of refdirs.
			###print( "in loop ", Blockdata["myid"],i,len(q))#,q

			for m in q:
				#print " m ",m,len(angles)

				assignments_of_refangles_to_angles[m] = find_assignments_of_refangles_to_angles(normals_set, oldparams[m], Tracker["an"])
				assignments_of_refangles_to_cones[i].extend(assignments_of_refangles_to_angles[m])

			assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i]))
			#print(  " assignments_of_refangles_to_cones on myid ",Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]) )
			assignments_of_refangles_to_cones[i] = sp_utilities.wrap_mpi_gatherv(assignments_of_refangles_to_cones[i], 0, Blockdata["shared_comm"])
			doit = 1
			if( Blockdata["myid_on_node"] == 0 ):
				assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i]))
				#print(  " assignments_of_refangles_to_cones on zero ",i,len(assignments_of_refangles_to_cones[i]) )#,assignments_of_refangles_to_cones[i]

		for i,q in enumerate(assignments_of_refangles_to_cones):
			###if( Blockdata["myid_on_node"] == 0 ): sxprint( " which refangles belong to which cone",Blockdata["color"],Blockdata["myid"],i,len(q) )#,q
			assignments_of_refangles_to_cones[i] = sp_utilities.wrap_mpi_bcast(q, 0, Blockdata["shared_comm"])

	else:

		angledirs = sp_utilities.angles_to_normals([u1[:3] for u1 in oldparams])

		number_of_cones = max(2, int(n_coarse_ang/numberofrefs_inmem*1.5 + 0.5))
		###if Blockdata["myid_on_node"] == 0:  sxprint( " LENGTHS  ",Blockdata["color"],nima,n_coarse_ang, numberofrefs_inmem,number_of_cones)
		cont = True
		while  cont :
			#  Translate number of cones to cone_delta
			cone_delta, number_of_cones = number_of_cones_to_delta(number_of_cones)
			#  Generate cone_angles
			##if Blockdata["myid"] == 0:  sxprint( "  WHILE  ",number_of_cones, cone_delta)
			if(Blockdata["symclass"].sym[0] == "c"):
				if( number_of_cones == 1 ):
					cone_delta  = 180.1
					cone_angles = [[0., 1.0, 0.]]
				else:
					cone_angles = Blockdata["symclass"].even_angles(cone_delta, theta1=1.0, theta2=89.0)
					cone_angles += [[(q[0]+90.0)%360., 180.-q[1],0] for q in cone_angles]
					cone_angles = Blockdata["symclass"].reduce_anglesets(cone_angles)
			else:
				cone_angles = Blockdata["symclass"].even_angles(cone_delta, theta1=1.0)

			#if Blockdata["myid"] == 0:  sxprint(  "  number of cones ",number_of_cones,cone_delta, len(cone_angles))
			assert(number_of_cones == len(cone_angles) )

			conedirs = sp_utilities.angles_to_normals(Blockdata["symclass"].symmetry_neighbors(cone_angles))
			neighbors = len(conedirs)/len(cone_angles)  #  Symmetry related neighbors
			#if Blockdata["myid"] == 0:  sxprint(  "  neighbors  ",Blockdata["myid"],neighbors, cone_angles)
			#  assign data directions to cone_angles
			assignments_to_cones = sp_utilities.assign_projdirs_f(angledirs, conedirs, neighbors)
			###print(  " assignments_to_cones ",Blockdata["myid"],len(assignments_to_cones),[len(q) for q in assignments_to_cones],assignments_to_cones[0])
			#  the above should have length of refdirs and each should have indexes of data that belong to this cone
			del conedirs
			#print "assignments_to_cones ",assignments_to_cones
			#  For each cone we have to find which refangles are needed to do the matching
			assignments_of_refangles_to_cones = [[] for i in range(len(assignments_to_cones))]
			assignments_of_refangles_to_angles = [[] for i in range(nima)]  # for each myid separately, these are angles on this myid


			for i,q in enumerate(assignments_to_cones):
				#  find assignments of refdirs to this cone within an of each angledirs, so we need symmetry neighbors set of refdirs.
				#if Blockdata["myid"] == 0:  sxprint( "in loop ", Blockdata["myid"],i,len(q),q)

				if(len(q) == 0):
					# empty assignment, on a given CPU there are no images assigned to a given cone
					assignments_of_refangles_to_cones[i] = [-1]
				else:
					for m in q:
						#print " m ",m,len(angles)

						assignments_of_refangles_to_angles[m] = find_assignments_of_refangles_to_angles(normals_set, oldparams[m], Tracker["an"])
						#if Blockdata["myid"] == 0:  sxprint( "assignments_of_refangles_to_angles[m] ", Blockdata["color"],i,m,assignments_of_refangles_to_angles[m])
						assignments_of_refangles_to_cones[i].extend(assignments_of_refangles_to_angles[m])

					#if( Blockdata["myid"] == 19 ):  sxprint(  " doit0 ",Blockdata["myid"], i,assignments_of_refangles_to_cones[i],q,assignments_to_cones)
					assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i]))
				###if Blockdata["myid"] == 19:  sxprint(  " assignments_of_refangles_to_cones on myid ",Blockdata["color"],Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]), assignments_of_refangles_to_cones[i] )
				assignments_of_refangles_to_cones[i] = sp_utilities.wrap_mpi_gatherv(assignments_of_refangles_to_cones[i], 0, Blockdata["shared_comm"])
				###if Blockdata["myid"] == 0:  sxprint(  " assignments_of_refangles_to_cones gatherv",Blockdata["color"],Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]) )
				doit = 1
				if( Blockdata["myid_on_node"] == 0 ):
					assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i])-set([-1]))
					###if( Blockdata["myid_on_node"] == 0 ):  sxprint(  " assignments_of_refangles_to_cones on zero ",i,len(assignments_of_refangles_to_cones[i]))
					#  POSSIBLE PROBLEM - IT IS POSSIBLE FOR A GIVEN CONE TO HAVE NO REFANGLES AND THUS BE EMPTY
					if( len(assignments_of_refangles_to_cones[i]) > numberofrefs_inmem ):
						number_of_cones = int(number_of_cones*1.25)
						#print(  " increased number_of_cones ",i,number_of_cones )
						doit = 0
				doit = sp_utilities.bcast_number_to_all(doit, source_node = 0)
				number_of_cones = sp_utilities.bcast_number_to_all(number_of_cones, source_node = 0, mpi_comm = Blockdata["shared_comm"] )
				###if( Blockdata["myid"] == 19 ):  sxprint(  " doit ",Blockdata["myid"], i,doit ,assignments_of_refangles_to_cones[i],assignments_to_cones)
				if( doit == 0 ):  break

			if( doit == 1 ):
				cont = False


		for i,q in enumerate(assignments_of_refangles_to_cones):
			###if( Blockdata["myid_on_node"] == 0 ): sxprint( " which refangles belong to which cone IOUT ",Blockdata["color"],Blockdata["myid"],i,len(q))#,q
			assignments_of_refangles_to_cones[i] = sp_utilities.wrap_mpi_bcast(q, 0, Blockdata["shared_comm"])
			###if( Blockdata["myid_on_node"] == 0 ): sxprint( " which refangles belong to which cone XOUT ",Blockdata["color"],Blockdata["myid"],i,len(q))#,q
			#if( myid == 1 ):
			#	print " which refangles belong to which cone",myid,i,len(assignments_of_refangles_to_cones[i])#,q

	if( Blockdata["myid"] == 0  ):  sp_global_def.sxprint( " number_of_cones: ",number_of_cones)
	
	#mpi_barrier(MPI_COMM_WORLD)
	#mpi_finalize()
	#exit()
	
	#  Maximum number of refangles assigned to angles (max number of references per image)
	nlocal_angles = max( [ len(q) for q in assignments_of_refangles_to_angles] )
	#  Most likely we have to delete some lists before proceeding
	del normals_set
	# We have to figure the maximum length of xod1, which is lang.  If there are no cones, it is an estimate.  If there are cones, we have list of assignments
	#  For number of cones I use refang and an.  This should give approximately the same number as coarse angles and 2*an, which is what is actually used in searches
	numberofrefs_inmem = max([len(q) for q in assignments_of_refangles_to_cones])

	###for i,q in enumerate(assignments_of_refangles_to_cones):
	###	sxprint( " which refangles belong to which cone OUT ",Blockdata["color"],Blockdata["myid_on_node"],Blockdata["myid"],i,numberofrefs_inmem,nlocal_angles,len(q))#,q

	#  BIG BUFFER
	lenbigbuf = numberofrefs_inmem  #MAXIMUM NUMBER OF REFERENCE IN ONE OF THE CONES
	orgsize = lenbigbuf*size_of_one_image #  This is number of projections to be computed simultaneously times their size

	if( Blockdata["myid_on_node"] == 0 ): size = orgsize
	else:  size = 0

	win_sm, base_ptr  = mpi.mpi_win_allocate_shared( size*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
	size = orgsize
	if( Blockdata["myid_on_node"] != 0 ):
		base_ptr, = mpi.mpi_win_shared_query(win_sm, mpi.MPI_PROC_NULL)

	buffer = numpy.frombuffer(numpy.core.multiarray.int_asbuffer(base_ptr, size*disp_unit), dtype = 'f4')
	buffer = buffer.reshape(lenbigbuf, size_of_one_image)
	#bigbuffer = EMNumPy.assign_numpy_to_emdata(buffer)

	emnumpy2 = EMAN2_cppwrap.EMNumPy()
	bigbuffer = emnumpy2.register_numpy_to_emdata(buffer)

	#  end of CONES setup

	if( Blockdata["myid"] == Blockdata["main_node"] ):
		sp_global_def.sxprint( "  " )
		line = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime()) + " =>"
		sp_global_def.sxprint(  line, "Processing data for polar onx: %3d, nx: %3d, nxpolar: %3d, CTF: %s, preshift: %s,  MEM: %6.2fGB."%(Tracker["constants"]["nnxo"], Tracker["nxinit"], Tracker["nxpolar"], Tracker["constants"]["CTF"], preshift,orgsize/1.0e9) )

	#  Note these are in Fortran notation for polar searches
	txm    	= float(Tracker["nxpolar"]-(Tracker["nxpolar"]//2+1) - radius)
	txl    	= float(radius - Tracker["nxpolar"]//2+1)

	if Blockdata["bckgnoise"] :
		oneover = []
		nnx = Blockdata["bckgnoise"][0].get_xsize()
		for i in range(len(Blockdata["bckgnoise"])):
			temp = [0.0]*nnx
			for k in range(nnx):
				if( Blockdata["bckgnoise"][i].get_value_at(k) > 0.0):  temp[k] = 1.0/numpy.sqrt(Blockdata["bckgnoise"][i].get_value_at(k))
			oneover.append(temp)
		del temp

	accumulatepw = [0.0]*nima
	norm_per_particle = [0.0]*nima

	##lxod1 = lang*len(list_of_coarse_shifts)*(int(2*Tracker["an"]/coarse_delta+0.5)+1)
	newpar = [[i, [0.0], []] for i in range(nima)]

	#  This is for auxiliary function searches.
	Blockdata["angle_set"] = refang
	#  Extract normals from rotation matrices
	refdirs = sp_utilities.angles_to_normals(refang)


	#  We have to make sure the number of cones is the same on all nodes, this is due to the strange MPI problem
	#   that forced me to insert overall barrier into iterations over cones
	max_number_of_cones = number_of_cones
	max_number_of_cones = mpi.mpi_reduce(max_number_of_cones, 1, mpi.MPI_INT, mpi.MPI_MAX, 0, mpi.MPI_COMM_WORLD)
	if Blockdata["myid"] == Blockdata["main_node"] :  max_number_of_cones = int(max_number_of_cones[0])
	max_number_of_cones = mpi.mpi_bcast(max_number_of_cones, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD)
	max_number_of_cones = int(max_number_of_cones[0])



	#if( Blockdata["myid_on_node"] == 0):
	#	sxprint("  FIRST  nima, nang, npsi, nshift, n_coarse_ang, nlocal_angles, n_coarse_psi, len(list_of_coarse_shifts), number_of_cones , max_number_of_cones ",nima, nang,npsi,nshifts,n_coarse_ang,nlocal_angles,n_coarse_psi, len(coarse_shifts),number_of_cones,max_number_of_cones)

	##firsti = True
	at = time.time()
	##eat = 0.0
	lima = 0  #  total counter of images
	#  PROCESSING OF CONES
	keepfirst_error = 0
	for icone in range(max_number_of_cones):
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		if( icone < number_of_cones ):  #  This is executed for individual number of cones, some nodes may have fewer.
			nang_start, nang_end = sp_applications.MPI_start_end(len(assignments_of_refangles_to_cones[icone]), Blockdata["no_of_processes_per_group"], Blockdata["myid_on_node"])
			#if(Blockdata["color"] == 1):
			###print( " ZXZ11  ",Blockdata["color"],Blockdata["myid_on_node"],Blockdata["myid"],len(assignments_of_refangles_to_cones[icone]),nang_start, nang_end)

			for i in range(nang_start, nang_end, 1):  # This will take care of no of process on a node less than nang.  Some loops will not be executed
				ic = assignments_of_refangles_to_cones[icone][i]
				temp = sp_projection.prgl(volprep,[ coarse_angles[ic][0],coarse_angles[ic][1],0.0, 0.0,0.0], 1, True)
				crefim = EMAN2_cppwrap.Util.Polar2Dm(temp, cnx, cnx, numr, mode)
				EMAN2_cppwrap.Util.Frngs(crefim, numr)
				EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)
				bigbuffer.insert_clip(crefim,(0,i) )

			mpi.mpi_barrier(Blockdata["shared_comm"])

			#if(Blockdata["myid"] == Blockdata["main_node"]):
			#	sxprint( "  Reference projections generated for cone %4d: %10.1fmin"%(icone,(time()-at)/60.))
			###print("  TOPOP    ",Blockdata["myid"],MPI_COMM_WORLD,len(assignments_to_cones[icone]))
			###mpi_finalize()
			###exit()
			###  <><><><><><><><>

			#  Preprocess the data

			#  We only process images in the current cone.  icnm is consecutive number, im the actual image number
			#for icnm,im in enumerate(assignments_to_cones[icone]):
			lenass = len(assignments_to_cones[icone])
			###print("   ENTERING  ",Blockdata["myid"],icone,lenass)
			for icnm in range(max(1,lenass)):  # I have to enter the loop even it there is no assignment
				if( lenass == 0 ):
					keepf = -1
					###print("  FOUNDEMPTY  ",Blockdata["myid"],icone,icnm,len(assignments_to_cones[icone]),assignments_to_cones)
				else:
					#  PREPARE ONE IMAGE
					im = assignments_to_cones[icone][icnm]

					phi,theta,psi,sx,sy, wnorm = oldparams[im][0], oldparams[im][1], oldparams[im][2], oldparams[im][3], oldparams[im][4], oldparams[im][7]

					#  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):  sxprint("\n\n  INPUT PARAMS  ",im,phi,theta,psi,sx,sy)
					if preshift:
						sx = int(round(sx))
						sy = int(round(sy))
						dataimage  = sp_fundamentals.cyclic_shift(original_data[im],sx,sy)
						#  Put rounded shifts on the list, note image has the original floats - check whether it may cause problems
						oldparams[im][3] = sx
						oldparams[im][4] = sy
						sx = 0.0
						sy = 0.0
					else:  dataimage = original_data[im].copy()


					st = EMAN2_cppwrap.Util.infomask(dataimage, mask2D, False)
					dataimage -= st[0]
					dataimage /= st[1]
					if dataimage.get_attr_default("bckgnoise", None) :  dataimage.delete_attr("bckgnoise")
					#  Do bckgnoise if exists
					if Blockdata["bckgnoise"]:
						if apply_mask:
							if Tracker["constants"]["hardmask"]:
								dataimage = sp_morphology.cosinemask(dataimage,radius = Tracker["constants"]["radius"])
							else:
								bckg = sp_utilities.model_gauss_noise(1.0,Tracker["constants"]["nnxo"]+2,Tracker["constants"]["nnxo"])
								bckg.set_attr("is_complex",1)
								bckg.set_attr("is_fftpad",1)
								bckg = sp_fundamentals.fft(sp_filter.filt_table(bckg, oneover[dataimage.get_attr("particle_group")]))
								#  Normalize bckg noise in real space, only region actually used.
								st = EMAN2_cppwrap.Util.infomask(bckg, mask2D, False)
								bckg -= st[0]
								bckg /= st[1]
								dataimage = sp_morphology.cosinemask(dataimage,radius = Tracker["constants"]["radius"], bckg = bckg)
					else:
						#  if no bckgnoise, do simple masking instead
						if apply_mask:  dataimage = sp_morphology.cosinemask(dataimage,radius = Tracker["constants"]["radius"] )

					#  Apply varadj
					if not nonorm:
						EMAN2_cppwrap.Util.mul_scalar(dataimage, Tracker["avgvaradj"][procid]/wnorm)

					###  FT
					dataimage = sp_fundamentals.fft(dataimage)
					sig = EMAN2_cppwrap.Util.rotavg_fourier( dataimage )
					accumulatepw[im] = sig[len(sig)//2:]+[0.0]

					#  We have to make sure the shifts are within correct range, shrinkage or not
					#set_params_proj(dataimage,[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
					
					if Blockdata["bckgnoise"]:
						temp = Blockdata["bckgnoise"][dataimage.get_attr("particle_group")]
						bckgn = EMAN2_cppwrap.Util.unroll1dpw(Tracker["nxinit"],Tracker["nxinit"], [temp[i] for i in range(temp.get_xsize())])
					else:
						bckgn = EMAN2_cppwrap.Util.unroll1dpw(Tracker["nxinit"],Tracker["nxinit"], [1.0]*600)
					bckgnoise = bckgn.copy()
					for j in range(Tracker["nxinit"]//2+1,Tracker["nxinit"]):  bckgn[0,j] = bckgn[0,Tracker["nxinit"]-j]

					if Tracker["constants"]["CTF"] :
						ctf_params = dataimage.get_attr("ctf")
						ctf_params.apix = ctf_params.apix/(float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"]))
						ctfa = sp_morphology.ctf_img_real(Tracker["nxinit"], ctf_params)
						ctfs = ctfa
					##if( ( Blockdata["myid"] == Blockdata["main_node"])   and firsti ):
					##	dataimage.set_attr("is_complex",0)
					##	dataimage.write_image("dataimagefft.hdf")
					##	dataimage.set_attr("is_complex",1)
					dataml = sp_fundamentals.fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False, False)
					data = []
					for iq in coarse_shifts:
						xx = iq[0]*shrink
						yy = iq[1]*shrink
						dss = sp_fundamentals.fshift(dataml, xx, yy)
						dss.set_attr("is_complex",0)
						data.append(dss)

					#  This will get it to real space
					#dataimage = fpol(Util.mulnclreal(Util.mulnclreal(fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False), Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)
					#  Here we can reuse dataml to save time.  It was not normalized, but it does not matter as dataimage is for POLAR.  PAP 08/06/2018
					dataimage = sp_fundamentals.fpol(EMAN2_cppwrap.Util.mulnclreal(EMAN2_cppwrap.Util.mulnclreal(dataml, EMAN2_cppwrap.Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)

					# Compute max number of angles on the fly
					lang = len(assignments_of_refangles_to_angles[im])
					###print("   BICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)

				if( ( Blockdata["myid"] == Blockdata["main_node"])  and  (lima%(max(1,nima/5)) == 0) and (lima>0)):
					sp_global_def.sxprint( "  Number of images :%7d   %5d  %5.1f"%(lima,nima,float(lima)/float(nima)*100.) + "%" +"   %10.1fmin"%((time.time()-at)/60.))
					##print( "  Number of images :%7d   %5d  %5.1f"%(lima,nima,float(lima)/float(nima)*100.) + "%" +"   %10.1fmin   %10.1fmin"%((time()-at)/60.,eat/60.0))
				lima += 1

				###print("  CONA1    ",Blockdata["myid"],lima)
				if( lima == 1 and procid == 0):
					###print("  CONA2    ",Blockdata["myid"])
					if( lenass > 0):
						###print("   CICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)
						keepfirst = lang*m_coarse_psi*n_coarse_shifts#150#lang*m_coarse_psi*len(coarse_shifts_shrank)  #500#min(200,nlocal_angles)#min(max(lxod1/25,200),lxod1)
						###if( Blockdata["myid"] == 18 and lima<5):  sxprint(" START   nlocal_angles* m_coarse_psi*len(coarse_shifts_shrank",nlocal_angles,coarse_delta,len(coarse_shifts_shrank),keepfirst)
						#if( Blockdata["myid"] == Blockdata["main_node"] ):  
						lxod1 = EMAN2_cppwrap.Util.multiref_Crosrng_msg_stack_stepsi_local(dataimage, bigbuffer, \
								coarse_shifts_shrank,\
								assignments_of_refangles_to_angles[im], assignments_of_refangles_to_cones[icone],\
								numr, [coarse_angles[i][2] for i in range(n_coarse_ang)], \
								oldparams[im][2], c_coarse_psi, coarse_delta, cnx, keepfirst)

						##'''
						assert(len(lxod1)/3 == keepfirst)

						xod1 = numpy.ndarray((keepfirst),dtype='f4',order="C")
						#xod1.fill(1.0)
						xod2 = numpy.ndarray((keepfirst),dtype='int',order="C")
						for iq in range(keepfirst):
							ioffset = 3*iq
							#          ishift         iang                      ipsi
							xod2[iq] = lxod1[ioffset] + lxod1[ioffset+1]*100000000 + lxod1[ioffset+2]*1000

						##'''
						#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

						#  Second step - find which coarse ones are significant

						# DO NOT order by angular directions to save time on reprojections.

						pre_ipsiandiang = -1
						for iln in range(keepfirst):
							#print("iln", iln)
							hashparams	= int(xod2[iln])
							ishift		= hashparams%1000
							ipsiandiang	= hashparams/1000
							if(ipsiandiang != pre_ipsiandiang):
								pre_ipsiandiang = ipsiandiang
								ipsi = ipsiandiang%100000
								iang = ipsiandiang/100000
								##junk = time()
								temp = sp_projection.prgl(volinit,[coarse_angles[iang][0],coarse_angles[iang][1],(coarse_angles[iang][2] + ipsi*coarse_delta)%360.0, 0.0,0.0], 1, False)
								##eat += time()-junk
								###   if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  SELECTEDSTEPTWO  ",iln,refang[iang][0],refang[iang][1],refang[iang][2]+ipsi*Tracker["delta"],ishift,xod1[iln])
								temp.set_attr("is_complex",0)
								nrmref = numpy.sqrt(EMAN2_cppwrap.Util.innerproduct(temp, temp, None))
							##junk = time()
							#xod1[iln] = -Util.sqed(data[ishift], temp, ctfa, bckgnoise)
							xod1[iln] = EMAN2_cppwrap.Util.innerproduct(data[ishift], temp, None)/nrmref
							##eat += time()-junk
							##xod2[iln] = hashparams

						xod1 -= numpy.max(xod1)
						lina = numpy.argwhere(xod1 > Tracker["constants"]["expthreshold"])
						#print("JJJ, expthreshold", Tracker["constants"]["expthreshold"])
						#print("LLL", xod1)
						#if( Blockdata["myid"] == 0 ):  
						###print("  STARTING1   ",Blockdata["myid"],np.max(xod1),np.min(xod1),len(lina),lina[-1])


						lina = lina.reshape(lina.size)
						
						keepf = int(lina[-1]) + 1
						#print("QQQQ, keepf", Blockdata["myid"], keepf)
						xod1 = xod1[lina]
						xod2 = xod2[lina]

						###print("  STARTING2    ",Blockdata["myid"])

						lina = numpy.argsort(xod1)
						xod1 = xod1[lina[::-1]]  # This sorts in reverse order
						xod2 = xod2[lina[::-1]]  # This sorts in reverse order
						numpy.exp(xod1, out=xod1)
						xod1 /= numpy.sum(xod1)
						cumprob = 0.0
						lit = len(xod1)
						###print("  STARTING3    ",Blockdata["myid"],lit)
						for j in range(len(xod1)):
							cumprob += xod1[j]
							if(cumprob > Tracker["constants"]["ccfpercentage"]):
								lit = j+1
								break


						###print("  STARTING4    ",Blockdata["myid"],lit,keepf)
					#  Turn into percentage of all possible
					keepf = [int(float(keepf*100)/float(keepfirst))]
					###mpi_barrier(MPI_COMM_WORLD)
					###print("  STARTING5    ",Blockdata["myid"],keepf,MPI_COMM_WORLD)
					mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
					keepf = sp_utilities.wrap_mpi_gatherv(keepf, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
					###print("  STARTING6    ",Blockdata["myid"],keepf)
					if( Blockdata["myid"] == 0 ):
						keepf = [junk for junk in keepf if junk >0]
						if( len(keepf) <2 ):
							keepf = 0
						else:
							keepf.sort()
							keepf = keepf[int(len(keepf)*Blockdata["rkeepf"])]
					else:  keepf = 0
					###print("  STARTING7    ",Blockdata["myid"],keepf)
					keepf = sp_utilities.wrap_mpi_bcast(keepf, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
					if(keepf == 0):
						keepf = 3
						if not keepfirst_error:
							sp_global_def.ERROR("\n\n###\n###\n###\nToo few images to estimate keepfirst.\nSet it to 3 for now.\nYou are entering unkown water!\nIn case the results do not look promising, try less processors.\n###\n###\n###\n\n","sxmeridien", 0, Blockdata["myid"])
							keepfirst_error = 1
						#ERROR("Too few images to estimate keepfirst","sxmeridien", 1, Blockdata["myid"])
						#mpi_finalize()
						#exit()
					###print("  STARTING8    ",Blockdata["myid"],keepf)
					Tracker["keepfirst"] = int(keepf)
					###if( Blockdata["myid"] == 0 ):  sxprint("  keepfirst first ",Tracker["keepfirst"])

				else:
					if(lenass>0):
						###print("   DICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)
						###if( Blockdata["myid"] == 18 and lima<5):  sxprint(" START   nlocal_angles* m_coarse_psi*len(coarse_shifts_shrank",nlocal_angles,coarse_delta,len(coarse_shifts_shrank),Tracker["keepfirst"])
						#if( Blockdata["myid"] == Blockdata["main_node"] ):
						keepfirst = lang*m_coarse_psi*n_coarse_shifts
						keepfirst = max(keepfirst*Tracker["keepfirst"]/100,min(keepfirst,3))
						lxod1 = EMAN2_cppwrap.Util.multiref_Crosrng_msg_stack_stepsi_local(dataimage, bigbuffer, \
								coarse_shifts_shrank,\
								assignments_of_refangles_to_angles[im], assignments_of_refangles_to_cones[icone],\
								numr, [coarse_angles[i][2] for i in range(n_coarse_ang)], \
								oldparams[im][2], c_coarse_psi, coarse_delta, cnx, keepfirst)
						#if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  ")
						##'''
						assert( keepfirst == len(lxod1)/3 )
						xod1 = numpy.ndarray((keepfirst),dtype='f4',order="C")
						#xod1.fill(1.0)
						xod2 = numpy.ndarray((keepfirst),dtype='int',order="C")
						for iq in range(keepfirst):
							ioffset = 3*iq
							#          ishift         iang                      ipsi
							xod2[iq] = lxod1[ioffset] + lxod1[ioffset+1]*100000000 + lxod1[ioffset+2]*1000

						##'''

						#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

						#  Second step - find which coarse ones are significant

						# order by angular directions to save time on reprojections.
						ipsiandiang = xod2/1000
						lina = numpy.argsort(ipsiandiang)
						xod2 = xod2[lina]  # order does not matter

						pre_ipsiandiang = -1
						for iln in range(keepfirst):
							hashparams	= int(xod2[iln])
							ishift		= hashparams%1000
							ipsiandiang	= hashparams/1000
							if(ipsiandiang != pre_ipsiandiang):
								pre_ipsiandiang = ipsiandiang
								ipsi = ipsiandiang%100000
								iang = ipsiandiang/100000
								##junk = time()
								temp = sp_projection.prgl(volinit,[coarse_angles[iang][0],coarse_angles[iang][1],(coarse_angles[iang][2] + ipsi*coarse_delta)%360.0, 0.0,0.0], 1, False)
								###   if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  SELECTEDSTEPTWO  ",iln,refang[iang][0],refang[iang][1],refang[iang][2]+ipsi*Tracker["delta"],ishift,xod1[iln])
								##eat += time()-junk
								temp.set_attr("is_complex",0)
							##junk = time()
							peak = -EMAN2_cppwrap.Util.sqed(data[ishift], temp, ctfa, bckgnoise)
							##eat += time()-junk
							#  Note I replace ccc by eqdist
							xod1[iln] = peak
							##xod2[iln] = hashparams

						lina = numpy.argsort(xod1)
						xod1 = xod1[lina[::-1]]  # This sorts in reverse order
						xod2 = xod2[lina[::-1]]  # This sorts in reverse order

						xod1 -= xod1[0]

						#if( Blockdata["myid"] == Blockdata["main_node"]):
						#	#print("  PROJECT   ",im,lit,johi)#,cod2)
						#	for iln in range(len(xod1)):  sxprint("  ROPE   ",iln,xod1[iln],xod2[iln])


						lina = numpy.argwhere(xod1 > Tracker["constants"]["expthreshold"])
						xod1 = xod1[lina]
						xod2 = xod2[lina]
						numpy.exp(xod1, out=xod1)
						xod1 /= numpy.sum(xod1)
						cumprob = 0.0
						lit = len(xod1)
						for j in range(len(xod1)):
							cumprob += xod1[j]
							if(cumprob > Tracker["constants"]["ccfpercentage"]):
								lit = j+1
								break

						#  To here
						###if( Blockdata["myid"] == 18 and lima<5):  sxprint("  SECOND KEPT  ",lit)
						#if( lima<5):  sxprint("  SECOND KEPT  ",lit)


				#mpi_barrier(MPI_COMM_WORLD)
				#mpi_finalize()
				#exit()
				#for j in range(lit):
				#	 newpar[kl][2].append([int(xod2[j]),float(xod1[j])])
				if( lenass > 0):
					###print("   EICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im)#,assignments_to_cones)

					firstdirections = [[0.0,0.0] for iln in range(lit)]
					firstshifts = [0]*lit
					for iln in range(lit):
						hashparams = int(xod2[iln])
						ishift = hashparams%1000
						ipsiandiang	= hashparams/1000
						#ipsi = ipsiandiang%100000
						iang = ipsiandiang/100000
						#try:
						firstdirections[iln] = [coarse_angles[iang][0], coarse_angles[iang][1], 0.0]
						#except:
						#	sxprint(" FAILED  ",Blockdata["myid"],Tracker["keepfirst"],iln,lima,lit,hashparams,iang,xod2[:lit],xod1[:lit])
						#	mpi_finalize()
						#	exit()
						firstshifts[iln] = ishift
						#  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):
						#	ipsi = ipsiandiang%100000
						#	ishift = hashparams%1000
						#	###print("  SECONDPAR%04d  "%im,iln,hashparams, ipsi,iang,ishift)
						#	sxprint(" SECONDPAR%04d  "%im,iln, refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, shifts[ishift])
					###del xod2
					###if( Blockdata["myid"] == 18 and lima<5):   sxprint("  SECOND KEPT  ",im,lit)#,len(xod2),xod2,firstdirections,firstshifts)
					#mpi_barrier(MPI_COMM_WORLD)
					#mpi_finalize()
					#exit()

					#if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  FIFI ",firstdirections)
					#if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  GUGU ",firstshifts)
					# Find neighbors, ltabang contains positions on refang list, no psis
					###ltabang = nearest_many_full_k_projangles(refang, firstdirections, howmany = 5, sym=Tracker["constants"]["symmetry"])
					ltabang = find_nearest_k_refangles_to_many_angles(refdirs, firstdirections, Tracker["delta"], howmany = 4)
					###if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  ltabang ",ltabang)
					##mpi_barrier(MPI_COMM_WORLD)
					##mpi_finalize()
					##exit()

					# ltabang has length lit, which is the same as length as firshifts.  However, original xod2 is still available,
					#   even though it is longer than lit.


					###ltabshi = [shiftneighbors[iln] for iln in firstshifts]
					#if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  HUHU ",ltabang)
					#if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  OUOU ",ltabshi)

					#  Prepare image for chi2.
					#  We have to repeat everything from get shrink data, including shifts
					#  The number of entries is lit=len(ltabang), each proj dir has 4 neighbors, so including itself there are 5,
					#    there are 3 psis, and at most n_fine_shifts. 
					#    Not not all shifts have to be populated. We may also have repeated triplets of angles, but each can have
					#      different shifts.  If so, we have to remove duplicates from the entire set.
					lcod1 = lit*4*2*n_fine_shifts

					#cod2 = np.ndarray((lit,5,3,n_fine_shifts),dtype=int,order="C")
					#cod2.fill(-1)  #  hashparams
					cod2 = []
					#lol = 0
					for i1 in range(lit):
						hashparams = int(xod2[i1])
						ipsiandiang	= hashparams/1000
						oldiang = ipsiandiang/100000
						ipsi = ipsiandiang%100000
						ishift = hashparams%1000
						tshifts = get_shifts_neighbors(shifts, coarse_shifts[ishift])
						#if( Blockdata["myid"] == Blockdata["main_node"]): sxprint(" tshifts  ",i1,len(tshifts))
						for i2 in range(4):
							iang = ltabang[i1][i2]
							for i3 in range(2):  # psi
								itpsi = int((coarse_angles[oldiang][2] + ipsi*coarse_delta - refang[iang][2]+360.0)/Tracker["delta"])
								itpsi = (itpsi + i3)%npsi
								for i4 in range(len(tshifts)):
									cod2.append(iang*100000000 + itpsi*1000 + tshifts[i4])
									#lol += 1
									#if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  zibzi  ",i1,i2,i3,i4, lol)


					del xod1, xod2

					###if( Blockdata["myid"] == 18 and lima<5):   sxprint("  THIRD   ",len(cod2))#,cod2)
					cod2 = list(set(cod2))
					cod1 = [[q/1000,i] for i,q in enumerate(cod2)]
					cod1.sort()

					lit = len(cod1)

					cod2 = numpy.asarray([cod2[cod1[i][1]] for i in range(lit)])

					cod1 = numpy.ndarray(lit,dtype='f4',order="C")
					#cod1.fill(np.finfo(dtype='f4').min)
					cod3 = numpy.ndarray(lit,dtype='f4',order="C")
					#cod3.fill(0.0)  #  varadj

					###if( Blockdata["myid"] == 18 and lima<5): sxprint("  THIRD   ",im,lit)#,cod2)

					#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					#  Third ML part.  Now the actual chi2 distances, we compute reprojections on the fly.
					#  Make sure volprep has nxinit size
					data = [None]*nshifts
					johi = 0
					iln = 0
					prevdir = -1
					while(iln<lit):
						hashparams = cod2[iln]
						#if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  COD2   ",im,lit,iln,cod2[iln])
						ipsiandiang	= hashparams/1000
						if(ipsiandiang != prevdir):
							prevdir = ipsiandiang
							ipsi = ipsiandiang%100000
							iang = ipsiandiang/100000
							##junk = time()
							temp = sp_projection.prgl(volinit,[ refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, 0.0,0.0], 1, False)
							##eat += time()-junk
							temp.set_attr("is_complex",0)
							johi += 1
						while( ipsiandiang == cod2[iln]/1000 ):
							hashparams = cod2[iln]
							ishift = hashparams%1000
							if( data[ishift] == None ):
								xx = shifts[ishift][0]*shrink
								yy = shifts[ishift][1]*shrink
								data[ishift] = sp_fundamentals.fshift(dataml, xx, yy)
								data[ishift].set_attr("is_complex",0)
							##junk = time()
							[peak,varadj] = EMAN2_cppwrap.Util.sqednorm(data[ishift], temp, ctfa, bckgnoise)
							##eat += time()-junk
							cod1[iln] = -peak
							cod3[iln] = varadj
							iln += 1
							if(iln == lit  ):  break
						#if( Blockdata["myid"] == Blockdata["main_node"]):
						#	temp.write_image("temp.hdf")
						#	data[iln].write_image("data.hdf")
						#	ctfa.write_image("ctfa.hdf")
						#	bckgnoise.write_image("bckgnoise.hdf")
						#	exit()
						###if( Blockdata["myid"] == Blockdata["main_node"] and iln%1000 ==0):  sxprint(" progress  ",iln,time()-at)
					#if( Blockdata["myid"] == Blockdata["main_node"]):
					#	sxprint("  PROJECT   ",im,lit,johi)#,cod2)
					#	#for iln in range(lit):  sxprint("  ROPE   ",iln,cod1[iln],cod2[iln],cod3[iln])
					del data
					del dataml

					lina = numpy.argsort(cod1)
					cod1 = cod1[lina[::-1]]  # This sorts in reverse order
					cod2 = cod2[lina[::-1]]  # This sorts in reverse order
					cod3 = cod3[lina[::-1]]  # This sorts in reverse order
					cod1 -= cod1[0]
					lina = numpy.argwhere(cod1 > Tracker["constants"]["expthreshold"])
					cod1 = cod1[lina]
					cod2 = cod2[lina]
					cod3 = cod3[lina]

					###if( Blockdata["myid"] == Blockdata["main_node"]):
					###for iui in range(len(lina)):
					###	for iui in range(len(cod1)):
					###		sxprint("  MLML  ",iui,cod1[iui],exp(cod1[iui]),cod2[iui],cod3[iui])

					numpy.exp(cod1, out=cod1)
					cod1 /= numpy.sum(cod1)
					cumprob = 0.0
					for j in range(len(cod1)):
						cumprob += cod1[j]
						if(cumprob > Tracker["constants"]["ccfpercentage"]):
							lit = j+1
							break

					#  New norm is a sum of eq distances multiplied by their probabilities augmented by PW.
					norm_per_particle[im] = numpy.sum(cod1[:lit]*cod3[:lit]) + accumulatepw[im][reachpw]
					###print("   CNORMPERPARTICLE  ",Blockdata["myid"],im,norm_per_particle[im])

					for iln in range(lit):
						newpar[im][2].append([int(cod2[iln]), float(cod1[iln])])
						#  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):
						#	#print("  NEWPAR%04d  "%im,iln,newpar[im][2][-1])
						#	hashparams = newpar[im][2][iln][0]
						#	ipsiandiang	= hashparams/1000
						#	ipsi = ipsiandiang%100000
						#	iang = ipsiandiang/100000
						#	ishift = hashparams%1000
						#	#print("  NEWPAR%04d  "%im,iln,newpar[im][2][-1],hashparams,ipsi,iang,ishift)
						#	sxprint(" NEWPAR%04d  "%im,iln,newpar[im][2][iln], refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, shifts[ishift])

					###if( Blockdata["myid"] == 18 and lima<5):   sxprint("  FINALLY  ",im,lit)
					del cod1, cod2, cod3, lina
					###mpi_barrier(MPI_COMM_WORLD)
					###mpi_finalize()
					###exit()

	"""Multiline Comment23"""
	#MULTILINEMULTILINEMULTILINE 23
	#MULTILINEMULTILINEMULTILINE 23
	#MULTILINEMULTILINEMULTILINE 23
	#MULTILINEMULTILINEMULTILINE 23
	#MULTILINEMULTILINEMULTILINE 23
	#MULTILINEMULTILINEMULTILINE 23
	#MULTILINEMULTILINEMULTILINE 23
	
	#  END OF CONES
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	if( Blockdata["myid"] == Blockdata["main_node"] ):
		sp_global_def.sxprint( "  Finished projection matching   %10.1fmin"%((time.time()-at)/60.))
	at = time.time()
	#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	#  All images were processed, now to the additional calculations

	###mpi_barrier(MPI_COMM_WORLD)
	###mpi_finalize()
	###exit()


	# norm correction ---- calc the norm correction per particle
	snormcorr = 0.0
	for kl in range(nima):
		###print("   NORMPERPARTICLE  ",Blockdata["myid"],kl,norm_per_particle[kl])
		norm_per_particle[kl]	= numpy.sqrt(norm_per_particle[kl]*2.0)*oldparams[kl][7]/Tracker["avgvaradj"][procid]
		snormcorr				+= norm_per_particle[kl]
	Tracker["avgvaradj"][procid] = snormcorr
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	#  Compute avgvaradj
	Tracker["avgvaradj"][procid] = mpi.mpi_reduce( Tracker["avgvaradj"][procid], 1, mpi.MPI_FLOAT, mpi.MPI_SUM, Blockdata["main_node"], mpi.MPI_COMM_WORLD )
	if(Blockdata["myid"] == Blockdata["main_node"]):
		Tracker["avgvaradj"][procid] = float(Tracker["avgvaradj"][procid])/Tracker["nima_per_chunk"][procid]
	else:  Tracker["avgvaradj"][procid] = 0.0
	Tracker["avgvaradj"][procid] = sp_utilities.bcast_number_to_all(Tracker["avgvaradj"][procid], Blockdata["main_node"])
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	
	at = time.time()
	mpi.mpi_barrier(Blockdata["shared_comm"])

	###if Blockdata["myid"] == Blockdata["main_node"]:  print "  Finished :",time()-at
	mpi.mpi_win_free(win_sm)
	if( Tracker["nxinit"] != Tracker["nxpolar"] ):
		mpi.mpi_win_free(win_volinit)
		emnumpy4.unregister_numpy_from_emdata()
		del emnumpy4
	else:   mpi.mpi_win_free(win_vol)

	mpi.mpi_barrier(Blockdata["shared_comm"])
	emnumpy1.unregister_numpy_from_emdata()
	emnumpy2.unregister_numpy_from_emdata()
	del emnumpy1, emnumpy2

	del volinit

	mpi.mpi_barrier(Blockdata["shared_comm"])

	#print("  NORMALIZATION DONE  ",Blockdata["myid"])
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	#if( Blockdata["myid"] == Blockdata["main_node"] ):
	#	sxprint( "  Projection matching finished : %10.1fmin"%((time()-at)/60.))
	return newpar, [1.0]*nima
























































































































































































































































































































































def XYXali3D_local_polar_ccc(refang, shifts, coarse_angles, coarse_shifts, procid, original_data = None, oldparams = None, \
					preshift = False, apply_mask = True, nonorm = False, applyctf = True):
	global Tracker, Blockdata
	global Tracker, Blockdata
	pass#IMPORTIMPORTIMPORT from sp_projection 	import prgs,prgl
	pass#IMPORTIMPORTIMPORT from sp_fundamentals 	import fft
	pass#IMPORTIMPORTIMPORT from sp_utilities 		import wrap_mpi_gatherv
	pass#IMPORTIMPORTIMPORT from math 			import sqrt
	#  Input data has to be CTF-multiplied, preshifted
	#  Output - newpar, see structure
	#    newpar = [[i, [worst_similarity, sum_all_similarities], [[-1, -1.0e23] for j in range(Tracker["lentop"])]] for i in range(len(data))]
	#    newpar = [[params],[],... len(data)]
	#    params = [particleID, [worst_similarity, sum_all_similarities],[imageallparams]]]
	#    imageallparams = [[orientation, similarity],[],...  number of all orientations ]
	#  Coding of orientations:
	#    hash = ang*100000000 + lpsi*1000 + ishift
	#    ishift = hash%1000
	#    ipsi = (hash/1000)%100000
	#    iang  = hash/100000000
	#  To get best matching for particle #kl:
	#     hash_best = newpar[kl][-1][0][0]
	#     best_sim  = newpar[kl][-1][0][1]
	#  To sort:
	pass#IMPORTIMPORTIMPORT from operator 		import itemgetter#, attrgetter, methodcaller
	#   params.sort(key=itemgetter(2))
	#from fundamentals import resample
	pass#IMPORTIMPORTIMPORT from sp_utilities    import get_im, model_gauss_noise, set_params_proj, get_params_proj
	pass#IMPORTIMPORTIMPORT from sp_fundamentals import fdecimate, fshift, fft
	pass#IMPORTIMPORTIMPORT from sp_filter       import filt_ctf, filt_table
	pass#IMPORTIMPORTIMPORT from sp_applications import MPI_start_end
	pass#IMPORTIMPORTIMPORT from math         import sqrt

	if(Blockdata["myid"] == Blockdata["main_node"]):
		print_dict(Tracker,"PROJECTION MATCHING parameters of buffered local polar")

	at = time.time()
	shrinkage = float(Tracker["nxpolar"])/float(Tracker["constants"]["nnxo"])
	shrink = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
	mode = "F"
	numr = Numrinit_local(1, radius, 1, mode)
	wr = sp_alignment.ringwe(numr, mode)
	cnx = float(Tracker["nxpolar"]//2 + 1)

	##if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  RADIUS IN POLAR  ",radius,numr[-1])
	#  FINE SEARCH CONSTANTS
	nang = len(refang)
	ac_fine = numpy.cos(numpy.radians(Tracker["an"]))
	npsi = int(360./Tracker["delta"])
	mpsi = 2
	c_fine_psi = mpsi//2
	nshifts = len(shifts)
	n_fine_shifts = 4

	nima = len(original_data)
	mask = EMAN2_cppwrap.Util.unrollmask(Tracker["nxinit"],Tracker["nxinit"])
	for j in range(Tracker["nxinit"]//2,Tracker["nxinit"]):  mask[0,j]=1.0
	mask2D = sp_utilities.model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	reachpw = mask.get_xsize()  # The last element of accumulated pw is zero so for the full size nothing is added.

	#  COARSE SEARCH CONSTANTS
	n_coarse_ang = len(coarse_angles)
	coarse_delta = 2*Tracker["delta"]
	n_coarse_psi = int(360./coarse_delta)
	###m_coarse_psi = int((2*2*Tracker["an"])/coarse_delta + 0.5) + 1
	m_coarse_psi = int(Tracker["an"]/Tracker["delta"] + 0.5)
	c_coarse_psi = m_coarse_psi//2
	n_coarse_shifts = len(coarse_shifts)

	coarse_shifts_shrank = [None]*n_coarse_shifts
	for ib in range(n_coarse_shifts):
		coarse_shifts_shrank[ib] = [coarse_shifts[ib][0]*shrinkage,coarse_shifts[ib][1]*shrinkage]

	###if(Blockdata["myid"] == Blockdata["main_node"]): sxprint("   TRETR  ",Tracker["constants"]["nnxo"],Tracker["nxinit"],reachpw,n_coarse_ang,coarse_delta,n_coarse_psi,m_coarse_psi,c_coarse_psi,n_coarse_shifts)

	#if(Blockdata["myid"] == Blockdata["main_node"]):
	#	sxprint( original_data[0].get_attr("identifier") )
	#	sxprint( original_data[1].get_attr("identifier") )

	disp_unit = numpy.dtype("f4").itemsize

	#  REFVOL
	if( Blockdata["myid_on_node"] == 0 ):
		if( Blockdata["myid"] == Blockdata["main_node"] ):
			odo = get_refvol( Tracker["nxpolar"] )
			nxvol = odo.get_xsize()
			nyvol = odo.get_ysize()
			nzvol = odo.get_zsize()
		else:
			nxvol = 0
			nyvol = 0
			nzvol = 0

		nxvol = sp_utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nyvol = sp_utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nzvol = sp_utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

		if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = sp_utilities.model_blank( nxvol,nyvol, nzvol)

		sp_utilities.bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

		odo = sp_projection.prep_vol( odo, npad = 2, interpolation_method = 1 )
		ndo = EMAN2_cppwrap.EMNumPy.em2numpy(odo)
		nxvol = odo.get_xsize()
		nyvol = odo.get_ysize()
		nzvol = odo.get_zsize()
		orgsizevol = nxvol*nyvol*nzvol
		sizevol = orgsizevol
	else:
		orgsizevol = 0
		sizevol = 0
		nxvol = 0
		nyvol = 0
		nzvol = 0

	orgsizevol = sp_utilities.bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
	nxvol = sp_utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
	nyvol = sp_utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
	nzvol = sp_utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

	win_vol, base_vol  = mpi.mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
	sizevol = orgsizevol
	if( Blockdata["myid_on_node"] != 0 ):
		base_vol, = mpi.mpi_win_shared_query(win_vol, mpi.MPI_PROC_NULL)

	volbuf = numpy.frombuffer(numpy.core.multiarray.int_asbuffer(base_vol, sizevol*disp_unit), dtype = 'f4')
	volbuf = volbuf.reshape(nzvol, nyvol, nxvol)
	if( Blockdata["myid_on_node"] == 0 ):
		numpy.copyto(volbuf,ndo)
		del odo,ndo

	#volprep = EMNumPy.assign_numpy_to_emdata(volbuf)
	emnumpy1 = EMAN2_cppwrap.EMNumPy()
	volprep = emnumpy1.register_numpy_to_emdata(volbuf)
	volprep.set_attr_dict({'is_complex':1,  'is_complex_x': 0, 'is_fftodd': 0, 'is_fftpad': 1, 'is_shuffled': 1,'npad': 2})

	if( Tracker["nxinit"] != Tracker["nxpolar"] ):
		mpi.mpi_win_free(win_vol)
		#  REFVOL FOR ML
		if( Blockdata["myid_on_node"] == 0 ):
			if( Blockdata["myid"] == Blockdata["main_node"] ):
				odo = get_refvol( Tracker["nxpolar"] )
				nxvol = odo.get_xsize()
				nyvol = odo.get_ysize()
				nzvol = odo.get_zsize()
			else:
				nxvol = 0
				nyvol = 0
				nzvol = 0

			nxvol = sp_utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
			nyvol = sp_utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
			nzvol = sp_utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

			if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = sp_utilities.model_blank( nxvol,nyvol, nzvol)

			sp_utilities.bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

			odo = sp_projection.prep_vol( odo, npad = 2, interpolation_method = 1 )
			ndo = EMAN2_cppwrap.EMNumPy.em2numpy(odo)
			nxvol = odo.get_xsize()
			nyvol = odo.get_ysize()
			nzvol = odo.get_zsize()
			orgsizevol = nxvol*nyvol*nzvol
			sizevol = orgsizevol
		else:
			orgsizevol = 0
			sizevol = 0
			nxvol = 0
			nyvol = 0
			nzvol = 0

		orgsizevol = sp_utilities.bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
		nxvol = sp_utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
		nyvol = sp_utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
		nzvol = sp_utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

		win_volinit, base_volinit  = mpi.mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
		sizevol = orgsizevol
		if( Blockdata["myid_on_node"] != 0 ):
			base_volinit, = mpi.mpi_win_shared_query(win_volinit, mpi.MPI_PROC_NULL)

		volbufinit = numpy.frombuffer(numpy.core.multiarray.int_asbuffer(base_volinit, sizevol*disp_unit), dtype = 'f4')
		volbufinit = volbufinit.reshape(nzvol, nyvol, nxvol)
		if( Blockdata["myid_on_node"] == 0 ):
			numpy.copyto(volbufinit,ndoinit)
			del odo,ndoinit

		#volinit = EMNumPy.assign_numpy_to_emdata(volbufinit)
		emnumpy4 = EMAN2_cppwrap.EMNumPy()
		volinit = emnumpy4.register_numpy_to_emdata(volbufinit)
		volinit.set_attr_dict({'is_complex':1,  'is_complex_x': 0, 'is_fftodd': 0, 'is_fftpad': 1, 'is_shuffled': 1,'npad': 2})
		if( Blockdata["myid_on_node"] == 0 ):  volinit.update()
		mpi.mpi_barrier(Blockdata["shared_comm"])
		###if( Blockdata["myid"] < 10 ):
		###	sxprint(" VOLPREP  ",volinit[0],volinit[1],info(volinit))
	else:
		volinit = volprep
	#  End of replaced volprep

    ### local search is plugged in from here
	#  START CONES
	#  This has to be systematically done per node
	#
	crefim = EMAN2_cppwrap.Util.Polar2Dm(sp_utilities.model_blank(Tracker["nxpolar"],Tracker["nxpolar"]), cnx, cnx, numr, mode)
	size_of_one_image = crefim.get_xsize()
	#  We will assume half of the memory is available.  We will do it betteer later.
	numberofrefs_inmem = int(Tracker["constants"]["memory_per_node"]/4/((size_of_one_image*disp_unit)/1.0e9))
	####if( Blockdata["myid_on_node"] == 0  ):  sxprint( " MEMEST ", n_coarse_ang,numberofrefs_inmem)
	#  number of references that will fit into one mode
	normals_set = sp_utilities.angles_to_normals(coarse_angles)
	Blockdata["angle_set"] = coarse_angles
	if( n_coarse_ang <= numberofrefs_inmem ):
		number_of_cones = 1
		numberofrefs_inmem = n_coarse_ang
		assignments_to_cones = [list(range(len(oldparams)))]

		assignments_of_refangles_to_cones = [[] for i in range(len(assignments_to_cones))]
		assignments_of_refangles_to_angles = [[] for i in range(nima)]  # for each myid separately, these are angles on this myid

		for i,q in enumerate(assignments_to_cones):
			#  find assignments of refdirs to this cone within an of each angledirs, so we need symmetry neighbors set of refdirs.
			###print( "in loop ", Blockdata["myid"],i,len(q))#,q

			for m in q:
				#print " m ",m,len(angles)

				assignments_of_refangles_to_angles[m] = find_assignments_of_refangles_to_angles(normals_set, oldparams[m], Tracker["an"])
				assignments_of_refangles_to_cones[i].extend(assignments_of_refangles_to_angles[m])

			assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i]))
			#print(  " assignments_of_refangles_to_cones on myid ",Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]) )
			assignments_of_refangles_to_cones[i] = sp_utilities.wrap_mpi_gatherv(assignments_of_refangles_to_cones[i], 0, Blockdata["shared_comm"])
			doit = 1
			if( Blockdata["myid_on_node"] == 0 ):
				assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i]))
				#print(  " assignments_of_refangles_to_cones on zero ",i,len(assignments_of_refangles_to_cones[i]) )#,assignments_of_refangles_to_cones[i]

		for i,q in enumerate(assignments_of_refangles_to_cones):
			###if( Blockdata["myid_on_node"] == 0 ): sxprint( " which refangles belong to which cone",Blockdata["color"],Blockdata["myid"],i,len(q) )#,q
			assignments_of_refangles_to_cones[i] = sp_utilities.wrap_mpi_bcast(q, 0, Blockdata["shared_comm"])

	else:

		angledirs = sp_utilities.angles_to_normals([u1[:3] for u1 in oldparams])

		number_of_cones = max(2, int(n_coarse_ang/numberofrefs_inmem*1.5 + 0.5))
		###if Blockdata["myid_on_node"] == 0:  sxprint( " LENGTHS  ",Blockdata["color"],nima,n_coarse_ang, numberofrefs_inmem,number_of_cones)
		cont = True
		while  cont :
			#  Translate number of cones to cone_delta
			cone_delta, number_of_cones = number_of_cones_to_delta(number_of_cones)
			#  Generate cone_angles
			##if Blockdata["myid"] == 0:  sxprint( "  WHILE  ",number_of_cones, cone_delta)
			if(Blockdata["symclass"].sym[0] == "c"):
				if( number_of_cones == 1 ):
					cone_delta  = 180.1
					cone_angles = [[0., 1.0, 0.]]
				else:
					cone_angles = Blockdata["symclass"].even_angles(cone_delta, theta1=1.0, theta2=89.0)
					cone_angles += [[(q[0]+90.0)%360., 180.-q[1],0] for q in cone_angles]
					cone_angles = Blockdata["symclass"].reduce_anglesets(cone_angles)
			else:
				cone_angles = Blockdata["symclass"].even_angles(cone_delta, theta1=1.0)

			#if Blockdata["myid"] == 0:  sxprint(  "  number of cones ",number_of_cones,cone_delta, len(cone_angles))
			assert(number_of_cones == len(cone_angles) )

			conedirs = sp_utilities.angles_to_normals(Blockdata["symclass"].symmetry_neighbors(cone_angles))
			neighbors = len(conedirs)/len(cone_angles)  #  Symmetry related neighbors
			#if Blockdata["myid"] == 0:  sxprint(  "  neighbors  ",Blockdata["myid"],neighbors, cone_angles)
			#  assign data directions to cone_angles
			assignments_to_cones = sp_utilities.assign_projdirs_f(angledirs, conedirs, neighbors)
			###print(  " assignments_to_cones ",Blockdata["myid"],len(assignments_to_cones),[len(q) for q in assignments_to_cones],assignments_to_cones[0])
			#  the above should have length of refdirs and each should have indexes of data that belong to this cone
			del conedirs
			#print "assignments_to_cones ",assignments_to_cones
			#  For each cone we have to find which refangles are needed to do the matching
			assignments_of_refangles_to_cones = [[] for i in range(len(assignments_to_cones))]
			assignments_of_refangles_to_angles = [[] for i in range(nima)]  # for each myid separately, these are angles on this myid


			for i,q in enumerate(assignments_to_cones):
				#  find assignments of refdirs to this cone within an of each angledirs, so we need symmetry neighbors set of refdirs.
				#if Blockdata["myid"] == 0:  sxprint( "in loop ", Blockdata["myid"],i,len(q),q)

				if(len(q) == 0):
					# empty assignment, on a given CPU there are no images assigned to a given cone
					assignments_of_refangles_to_cones[i] = [-1]
				else:
					for m in q:
						#print " m ",m,len(angles)

						assignments_of_refangles_to_angles[m] = find_assignments_of_refangles_to_angles(normals_set, oldparams[m], Tracker["an"])
						#if Blockdata["myid"] == 0:  sxprint( "assignments_of_refangles_to_angles[m] ", Blockdata["color"],i,m,assignments_of_refangles_to_angles[m])
						assignments_of_refangles_to_cones[i].extend(assignments_of_refangles_to_angles[m])

					#if( Blockdata["myid"] == 19 ):  sxprint(  " doit0 ",Blockdata["myid"], i,assignments_of_refangles_to_cones[i],q,assignments_to_cones)
					assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i]))
				###if Blockdata["myid"] == 19:  sxprint(  " assignments_of_refangles_to_cones on myid ",Blockdata["color"],Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]), assignments_of_refangles_to_cones[i] )
				assignments_of_refangles_to_cones[i] = sp_utilities.wrap_mpi_gatherv(assignments_of_refangles_to_cones[i], 0, Blockdata["shared_comm"])
				###if Blockdata["myid"] == 0:  sxprint(  " assignments_of_refangles_to_cones gatherv",Blockdata["color"],Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]) )
				doit = 1
				if( Blockdata["myid_on_node"] == 0 ):
					assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i])-set([-1]))
					###if( Blockdata["myid_on_node"] == 0 ):  sxprint(  " assignments_of_refangles_to_cones on zero ",i,len(assignments_of_refangles_to_cones[i]))
					#  POSSIBLE PROBLEM - IT IS POSSIBLE FOR A GIVEN CONE TO HAVE NO REFANGLES AND THUS BE EMPTY
					if( len(assignments_of_refangles_to_cones[i]) > numberofrefs_inmem ):
						number_of_cones = int(number_of_cones*1.25)
						#print(  " increased number_of_cones ",i,number_of_cones )
						doit = 0
				doit = sp_utilities.bcast_number_to_all(doit, source_node = 0)
				number_of_cones = sp_utilities.bcast_number_to_all(number_of_cones, source_node = 0, mpi_comm = Blockdata["shared_comm"] )
				###if( Blockdata["myid"] == 19 ):  sxprint(  " doit ",Blockdata["myid"], i,doit ,assignments_of_refangles_to_cones[i],assignments_to_cones)
				if( doit == 0 ):  break

			if( doit == 1 ):
				cont = False


		for i,q in enumerate(assignments_of_refangles_to_cones):
			###if( Blockdata["myid_on_node"] == 0 ): sxprint( " which refangles belong to which cone IOUT ",Blockdata["color"],Blockdata["myid"],i,len(q))#,q
			assignments_of_refangles_to_cones[i] = sp_utilities.wrap_mpi_bcast(q, 0, Blockdata["shared_comm"])
			###if( Blockdata["myid_on_node"] == 0 ): sxprint( " which refangles belong to which cone XOUT ",Blockdata["color"],Blockdata["myid"],i,len(q))#,q
			#if( myid == 1 ):
			#	print " which refangles belong to which cone",myid,i,len(assignments_of_refangles_to_cones[i])#,q

	if( Blockdata["myid"] == 0  ):  sp_global_def.sxprint( " number_of_cones: ",number_of_cones)
	
	#mpi_barrier(MPI_COMM_WORLD)
	#mpi_finalize()
	#exit()
	
	#  Maximum number of refangles assigned to angles (max number of references per image)
	nlocal_angles = max( [ len(q) for q in assignments_of_refangles_to_angles] )
	#  Most likely we have to delete some lists before proceeding
	del normals_set
	# We have to figure the maximum length of xod1, which is lang.  If there are no cones, it is an estimate.  If there are cones, we have list of assignments
	#  For number of cones I use refang and an.  This should give approximately the same number as coarse angles and 2*an, which is what is actually used in searches
	numberofrefs_inmem = max([len(q) for q in assignments_of_refangles_to_cones])

	###for i,q in enumerate(assignments_of_refangles_to_cones):
	###	sxprint( " which refangles belong to which cone OUT ",Blockdata["color"],Blockdata["myid_on_node"],Blockdata["myid"],i,numberofrefs_inmem,nlocal_angles,len(q))#,q

	#  BIG BUFFER
	lenbigbuf = numberofrefs_inmem  #MAXIMUM NUMBER OF REFERENCE IN ONE OF THE CONES
	orgsize = lenbigbuf*size_of_one_image #  This is number of projections to be computed simultaneously times their size

	if( Blockdata["myid_on_node"] == 0 ): size = orgsize
	else:  size = 0

	win_sm, base_ptr  = mpi.mpi_win_allocate_shared( size*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
	size = orgsize
	if( Blockdata["myid_on_node"] != 0 ):
		base_ptr, = mpi.mpi_win_shared_query(win_sm, mpi.MPI_PROC_NULL)

	buffer = numpy.frombuffer(numpy.core.multiarray.int_asbuffer(base_ptr, size*disp_unit), dtype = 'f4')
	buffer = buffer.reshape(lenbigbuf, size_of_one_image)
	#bigbuffer = EMNumPy.assign_numpy_to_emdata(buffer)

	emnumpy2 = EMAN2_cppwrap.EMNumPy()
	bigbuffer = emnumpy2.register_numpy_to_emdata(buffer)

	#  end of CONES setup

	if( Blockdata["myid"] == Blockdata["main_node"] ):
		sp_global_def.sxprint( "  " )
		line = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime()) + " =>"
		sp_global_def.sxprint(  line, "Processing data for polar onx: %3d, nx: %3d, nxpolar: %3d, CTF: %s, preshift: %s,  MEM: %6.2fGB."%(Tracker["constants"]["nnxo"], Tracker["nxinit"], Tracker["nxpolar"], Tracker["constants"]["CTF"], preshift,orgsize/1.0e9) )

	#  Note these are in Fortran notation for polar searches
	txm    	= float(Tracker["nxpolar"]-(Tracker["nxpolar"]//2+1) - radius)
	txl    	= float(radius - Tracker["nxpolar"]//2+1)

	if Blockdata["bckgnoise"] :
		oneover = []
		nnx = Blockdata["bckgnoise"][0].get_xsize()
		for i in range(len(Blockdata["bckgnoise"])):
			temp = [0.0]*nnx
			for k in range(nnx):
				if( Blockdata["bckgnoise"][i].get_value_at(k) > 0.0):  temp[k] = 1.0/numpy.sqrt(Blockdata["bckgnoise"][i].get_value_at(k))
			oneover.append(temp)
		del temp

	accumulatepw = [0.0]*nima
	norm_per_particle = [0.0]*nima

	##lxod1 = lang*len(list_of_coarse_shifts)*(int(2*Tracker["an"]/coarse_delta+0.5)+1)
	newpar = [[i, [0.0], []] for i in range(nima)]

	#  This is for auxiliary function searches.
	Blockdata["angle_set"] = refang
	#  Extract normals from rotation matrices
	refdirs = sp_utilities.angles_to_normals(refang)


	#  We have to make sure the number of cones is the same on all nodes, this is due to the strange MPI problem
	#   that forced me to insert overall barrier into iterations over cones
	max_number_of_cones = number_of_cones
	max_number_of_cones = mpi.mpi_reduce(max_number_of_cones, 1, mpi.MPI_INT, mpi.MPI_MAX, 0, mpi.MPI_COMM_WORLD)
	if Blockdata["myid"] == Blockdata["main_node"] :  max_number_of_cones = int(max_number_of_cones[0])
	max_number_of_cones = mpi.mpi_bcast(max_number_of_cones, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD)
	max_number_of_cones = int(max_number_of_cones[0])



	#if( Blockdata["myid_on_node"] == 0):
	#	sxprint("  FIRST  nima, nang, npsi, nshift, n_coarse_ang, nlocal_angles, n_coarse_psi, len(list_of_coarse_shifts), number_of_cones , max_number_of_cones ",nima, nang,npsi,nshifts,n_coarse_ang,nlocal_angles,n_coarse_psi, len(coarse_shifts),number_of_cones,max_number_of_cones)

	##firsti = True
	at = time.time()
	##eat = 0.0
	lima = 0  #  total counter of images
	#  PROCESSING OF CONES
	keepfirst_error = 0
	for icone in range(max_number_of_cones):
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		if( icone < number_of_cones ):  #  This is executed for individual number of cones, some nodes may have fewer.
			nang_start, nang_end = sp_applications.MPI_start_end(len(assignments_of_refangles_to_cones[icone]), Blockdata["no_of_processes_per_group"], Blockdata["myid_on_node"])
			#if(Blockdata["color"] == 1):
			###print( " ZXZ11  ",Blockdata["color"],Blockdata["myid_on_node"],Blockdata["myid"],len(assignments_of_refangles_to_cones[icone]),nang_start, nang_end)

			for i in range(nang_start, nang_end, 1):  # This will take care of no of process on a node less than nang.  Some loops will not be executed
				ic = assignments_of_refangles_to_cones[icone][i]
				temp = sp_projection.prgl(volprep,[ coarse_angles[ic][0],coarse_angles[ic][1],0.0, 0.0,0.0], 1, True)
				crefim = EMAN2_cppwrap.Util.Polar2Dm(temp, cnx, cnx, numr, mode)
				EMAN2_cppwrap.Util.Frngs(crefim, numr)
				EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)
				bigbuffer.insert_clip(crefim,(0,i) )

			mpi.mpi_barrier(Blockdata["shared_comm"])

			#if(Blockdata["myid"] == Blockdata["main_node"]):
			#	sxprint( "  Reference projections generated for cone %4d: %10.1fmin"%(icone,(time()-at)/60.))
			###print("  TOPOP    ",Blockdata["myid"],MPI_COMM_WORLD,len(assignments_to_cones[icone]))
			###mpi_finalize()
			###exit()
			###  <><><><><><><><>

			#  Preprocess the data

			#  We only process images in the current cone.  icnm is consecutive number, im the actual image number
			#for icnm,im in enumerate(assignments_to_cones[icone]):
			lenass = len(assignments_to_cones[icone])
			###print("   ENTERING  ",Blockdata["myid"],icone,lenass)
			for icnm in range(max(1,lenass)):  # I have to enter the loop even it there is no assignment
				if( lenass == 0 ):
					keepf = -1
					###print("  FOUNDEMPTY  ",Blockdata["myid"],icone,icnm,len(assignments_to_cones[icone]),assignments_to_cones)
				else:
					#  PREPARE ONE IMAGE
					im = assignments_to_cones[icone][icnm]

					phi,theta,psi,sx,sy, wnorm = oldparams[im][0], oldparams[im][1], oldparams[im][2], oldparams[im][3], oldparams[im][4], oldparams[im][7]

					#  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):  sxprint("\n\n  INPUT PARAMS  ",im,phi,theta,psi,sx,sy)
					if preshift:
						sx = int(round(sx))
						sy = int(round(sy))
						dataimage  = sp_fundamentals.cyclic_shift(original_data[im],sx,sy)
						#  Put rounded shifts on the list, note image has the original floats - check whether it may cause problems
						oldparams[im][3] = sx
						oldparams[im][4] = sy
						sx = 0.0
						sy = 0.0
					else:  dataimage = original_data[im].copy()


					st = EMAN2_cppwrap.Util.infomask(dataimage, mask2D, False)
					dataimage -= st[0]
					dataimage /= st[1]
					if dataimage.get_attr_default("bckgnoise", None) :  dataimage.delete_attr("bckgnoise")
					#  Do bckgnoise if exists
					if Blockdata["bckgnoise"]:
						if apply_mask:
							if Tracker["constants"]["hardmask"]:
								dataimage = sp_morphology.cosinemask(dataimage,radius = Tracker["constants"]["radius"])
							else:
								bckg = sp_utilities.model_gauss_noise(1.0,Tracker["constants"]["nnxo"]+2,Tracker["constants"]["nnxo"])
								bckg.set_attr("is_complex",1)
								bckg.set_attr("is_fftpad",1)
								bckg = sp_fundamentals.fft(sp_filter.filt_table(bckg, oneover[dataimage.get_attr("particle_group")]))
								#  Normalize bckg noise in real space, only region actually used.
								st = EMAN2_cppwrap.Util.infomask(bckg, mask2D, False)
								bckg -= st[0]
								bckg /= st[1]
								dataimage = sp_morphology.cosinemask(dataimage,radius = Tracker["constants"]["radius"], bckg = bckg)
					else:
						#  if no bckgnoise, do simple masking instead
						if apply_mask:  dataimage = sp_morphology.cosinemask(dataimage,radius = Tracker["constants"]["radius"] )

					#  Apply varadj
					if not nonorm:
						EMAN2_cppwrap.Util.mul_scalar(dataimage, Tracker["avgvaradj"][procid]/wnorm)

					###  FT
					dataimage = sp_fundamentals.fft(dataimage)
					sig = EMAN2_cppwrap.Util.rotavg_fourier( dataimage )
					accumulatepw[im] = sig[len(sig)//2:]+[0.0]

					#  We have to make sure the shifts are within correct range, shrinkage or not
					#set_params_proj(dataimage,[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
					if Blockdata["bckgnoise"]:
						temp = Blockdata["bckgnoise"][dataimage.get_attr("particle_group")]
						bckgn = EMAN2_cppwrap.Util.unroll1dpw(Tracker["nxinit"],Tracker["nxinit"], [temp[i] for i in range(temp.get_xsize())])
					else:
						bckgn = EMAN2_cppwrap.Util.unroll1dpw(Tracker["nxinit"],Tracker["nxinit"], [1.0]*600)
					bckgnoise = bckgn.copy()
					for j in range(Tracker["nxinit"]//2+1,Tracker["nxinit"]):  bckgn[0,j] = bckgn[0,Tracker["nxinit"]-j]

					if Tracker["constants"]["CTF"] :
						ctf_params = dataimage.get_attr("ctf")
						ctf_params.apix = ctf_params.apix/(float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"]))
						ctfa = sp_morphology.ctf_img_real(Tracker["nxinit"], ctf_params)
						ctfs = ctfa
					##if( ( Blockdata["myid"] == Blockdata["main_node"])   and firsti ):
					##	dataimage.set_attr("is_complex",0)
					##	dataimage.write_image("dataimagefft.hdf")
					##	dataimage.set_attr("is_complex",1)
					dataml = sp_fundamentals.fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False, False)
					data = []
					for iq in coarse_shifts:
						xx = iq[0]*shrink
						yy = iq[1]*shrink
						dss = sp_fundamentals.fshift(dataml, xx, yy)
						dss.set_attr("is_complex",0)
						data.append(dss)

					#  This will get it to real space
					#dataimage = fpol(Util.mulnclreal(Util.mulnclreal(fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False), Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)
					#  Here we can reuse dataml to save time.  It was not normalized, but it does not matter as dataimage is for POLAR.  PAP 08/06/2018
					dataimage = sp_fundamentals.fpol(EMAN2_cppwrap.Util.mulnclreal(EMAN2_cppwrap.Util.mulnclreal(dataml, EMAN2_cppwrap.Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)

					# Compute max number of angles on the fly
					lang = len(assignments_of_refangles_to_angles[im])
					###print("   BICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)

				if( ( Blockdata["myid"] == Blockdata["main_node"])  and  (lima%(max(1,nima/5)) == 0) and (lima>0)):
					sp_global_def.sxprint( "  Number of images :%7d   %5d  %5.1f"%(lima,nima,float(lima)/float(nima)*100.) + "%" +"   %10.1fmin"%((time.time()-at)/60.))
					##print( "  Number of images :%7d   %5d  %5.1f"%(lima,nima,float(lima)/float(nima)*100.) + "%" +"   %10.1fmin   %10.1fmin"%((time()-at)/60.,eat/60.0))
				lima += 1



				###print("  CONA1    ",Blockdata["myid"],lima)
				if( lima == 1 and procid == 0):
					###print("  CONA2    ",Blockdata["myid"])
					if( lenass > 0):
						###print("   CICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)
						keepfirst = lang*m_coarse_psi*n_coarse_shifts#150#lang*m_coarse_psi*len(coarse_shifts_shrank)  #500#min(200,nlocal_angles)#min(max(lxod1/25,200),lxod1)
						###if( Blockdata["myid"] == 18 and lima<5):  sxprint(" START   nlocal_angles* m_coarse_psi*len(coarse_shifts_shrank",nlocal_angles,coarse_delta,len(coarse_shifts_shrank),keepfirst)
						#if( Blockdata["myid"] == Blockdata["main_node"] ):  
						lxod1 = EMAN2_cppwrap.Util.multiref_Crosrng_msg_stack_stepsi_local(dataimage, bigbuffer, \
								coarse_shifts_shrank,\
								assignments_of_refangles_to_angles[im], assignments_of_refangles_to_cones[icone],\
								numr, [coarse_angles[i][2] for i in range(n_coarse_ang)], \
								oldparams[im][2], c_coarse_psi, coarse_delta, cnx, keepfirst)

						##'''
						assert(len(lxod1)/3 == keepfirst)

						xod1 = numpy.ndarray((keepfirst),dtype='f4',order="C")
						#xod1.fill(1.0)
						xod2 = numpy.ndarray((keepfirst),dtype='int',order="C")
						for iq in range(keepfirst):
							ioffset = 3*iq
							#          ishift         iang                      ipsi
							xod2[iq] = lxod1[ioffset] + lxod1[ioffset+1]*100000000 + lxod1[ioffset+2]*1000

						##'''
						#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

						#  Second step - find which coarse ones are significant

						# order by angular directions to save time on reprojections.

						pre_ipsiandiang = -1
						for iln in range(keepfirst):
							#print("iln", iln)
							hashparams	= int(xod2[iln])
							ishift		= hashparams%1000
							ipsiandiang	= hashparams/1000
							if(ipsiandiang != pre_ipsiandiang):
								pre_ipsiandiang = ipsiandiang
								ipsi = ipsiandiang%100000
								iang = ipsiandiang/100000
								##junk = time()
								temp = sp_projection.prgl(volinit,[coarse_angles[iang][0],coarse_angles[iang][1],(coarse_angles[iang][2] + ipsi*coarse_delta)%360.0, 0.0,0.0], 1, False)
								##eat += time()-junk
								###   if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  SELECTEDSTEPTWO  ",iln,refang[iang][0],refang[iang][1],refang[iang][2]+ipsi*Tracker["delta"],ishift,xod1[iln])
								temp.set_attr("is_complex",0)
								nrmref = numpy.sqrt(EMAN2_cppwrap.Util.innerproduct(temp, temp, None))
							##junk = time()
							#xod1[iln] = -Util.sqed(data[ishift], temp, ctfa, bckgnoise)
							xod1[iln] = EMAN2_cppwrap.Util.innerproduct(data[ishift], temp, None)/nrmref
							##eat += time()-junk
							##xod2[iln] = hashparams

						xod1 -= numpy.max(xod1)
						lina = numpy.argwhere(xod1 > Tracker["constants"]["expthreshold"])
						#print("JJJ, expthreshold", Tracker["constants"]["expthreshold"])
						#print("LLL", xod1)
						#if( Blockdata["myid"] == 0 ):  
						###print("  STARTING1   ",Blockdata["myid"],np.max(xod1),np.min(xod1),len(lina),lina[-1])


						lina = lina.reshape(lina.size)
						
						keepf = int(lina[-1]) + 1
						#print("QQQQ, keepf", Blockdata["myid"], keepf)
						xod1 = xod1[lina]
						xod2 = xod2[lina]

						###print("  STARTING2    ",Blockdata["myid"])

						lina = numpy.argsort(xod1)
						xod1 = xod1[lina[::-1]]  # This sorts in reverse order
						xod2 = xod2[lina[::-1]]  # This sorts in reverse order
						numpy.exp(xod1, out=xod1)
						xod1 /= numpy.sum(xod1)
						cumprob = 0.0
						lit = len(xod1)
						###print("  STARTING3    ",Blockdata["myid"],lit)
						for j in range(len(xod1)):
							cumprob += xod1[j]
							if(cumprob > Tracker["constants"]["ccfpercentage"]):
								lit = j+1
								break


						###print("  STARTING4    ",Blockdata["myid"],lit,keepf)
					keepf = [keepf]
					###mpi_barrier(MPI_COMM_WORLD)
					###print("  STARTING5    ",Blockdata["myid"],keepf,MPI_COMM_WORLD)
					mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
					keepf = sp_utilities.wrap_mpi_gatherv(keepf, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
					###print("  STARTING6    ",Blockdata["myid"],keepf)
					if( Blockdata["myid"] == 0 ):
						keepf = [junk for junk in keepf if junk >0]
						if( len(keepf) <2 ):
							keepf = 0
						else:
							keepf.sort()
							keepf = keepf[int(len(keepf)*Blockdata["rkeepf"])]
					else:  keepf = 0
					###print("  STARTING7    ",Blockdata["myid"],keepf)
					keepf = sp_utilities.wrap_mpi_bcast(keepf, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
					if(keepf == 0):
						keepf = 3
						if not keepfirst_error:
							sp_global_def.ERROR("\n\n###\n###\n###\nToo few images to estimate keepfirst.\nSet it to 3 for now.\nYou are entering unkown water!\nIn case the results do not look promising, try less processors.\n###\n###\n###\n\n","sxmeridien", 0, Blockdata["myid"])
							keepfirst_error = 1
						#ERROR("Too few images to estimate keepfirst","sxmeridien", 1, Blockdata["myid"])
						#mpi_finalize()
						#exit()
					###print("  STARTING8    ",Blockdata["myid"],keepf)
					Tracker["keepfirst"] = int(keepf)
					if( Blockdata["myid"] == 0 ):  sp_global_def.sxprint("  keepfirst first ",Tracker["keepfirst"])

				else:
					if(lenass>0):
						###print("   DICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)
						###if( Blockdata["myid"] == 18 and lima<5):  sxprint(" START   nlocal_angles* m_coarse_psi*len(coarse_shifts_shrank",nlocal_angles,coarse_delta,len(coarse_shifts_shrank),Tracker["keepfirst"])
						#if( Blockdata["myid"] == Blockdata["main_node"] ):  
						lxod1 = EMAN2_cppwrap.Util.multiref_Crosrng_msg_stack_stepsi_local(dataimage, bigbuffer, \
								coarse_shifts_shrank,\
								assignments_of_refangles_to_angles[im], assignments_of_refangles_to_cones[icone],\
								numr, [coarse_angles[i][2] for i in range(n_coarse_ang)], \
								oldparams[im][2], c_coarse_psi, coarse_delta, cnx, Tracker["keepfirst"])
						#if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  ")
						##'''

						xod1 = numpy.ndarray((Tracker["keepfirst"]),dtype='f4',order="C")
						#xod1.fill(1.0)
						xod2 = numpy.ndarray((Tracker["keepfirst"]),dtype='int',order="C")
						for iq in range(Tracker["keepfirst"]):
							ioffset = 3*iq
							#          ishift         iang                      ipsi
							xod2[iq] = lxod1[ioffset] + lxod1[ioffset+1]*100000000 + lxod1[ioffset+2]*1000

						##'''

						#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

						#  Second step - find which coarse ones are significant

						# order by angular directions to save time on reprojections.
						ipsiandiang = xod2/1000
						lina = numpy.argsort(ipsiandiang)
						xod2 = xod2[lina]  # order does not matter

						pre_ipsiandiang = -1
						for iln in range(Tracker["keepfirst"]):
							hashparams	= int(xod2[iln])
							ishift		= hashparams%1000
							ipsiandiang	= hashparams/1000
							if(ipsiandiang != pre_ipsiandiang):
								pre_ipsiandiang = ipsiandiang
								ipsi = ipsiandiang%100000
								iang = ipsiandiang/100000
								##junk = time()
								temp = sp_projection.prgl(volinit,[coarse_angles[iang][0],coarse_angles[iang][1],(coarse_angles[iang][2] + ipsi*coarse_delta)%360.0, 0.0,0.0], 1, False)
								###   if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  SELECTEDSTEPTWO  ",iln,refang[iang][0],refang[iang][1],refang[iang][2]+ipsi*Tracker["delta"],ishift,xod1[iln])
								##eat += time()-junk
								temp.set_attr("is_complex",0)
							##junk = time()
							peak = -EMAN2_cppwrap.Util.sqed(data[ishift], temp, ctfa, bckgnoise)
							##eat += time()-junk
							#  Note I replace ccc by eqdist
							xod1[iln] = peak
							##xod2[iln] = hashparams

						lina = numpy.argsort(xod1)
						xod1 = xod1[lina[::-1]]  # This sorts in reverse order
						xod2 = xod2[lina[::-1]]  # This sorts in reverse order

						xod1 -= xod1[0]

						#if( Blockdata["myid"] == Blockdata["main_node"]):
						#	#print("  PROJECT   ",im,lit,johi)#,cod2)
						#	for iln in range(len(xod1)):  sxprint("  ROPE   ",iln,xod1[iln],xod2[iln])


						lina = numpy.argwhere(xod1 > Tracker["constants"]["expthreshold"])
						xod1 = xod1[lina]
						xod2 = xod2[lina]
						numpy.exp(xod1, out=xod1)
						xod1 /= numpy.sum(xod1)
						cumprob = 0.0
						lit = len(xod1)
						for j in range(len(xod1)):
							cumprob += xod1[j]
							if(cumprob > Tracker["constants"]["ccfpercentage"]):
								lit = j+1
								break

						#  To here
						###if( Blockdata["myid"] == 18 and lima<5):  sxprint("  SECOND KEPT  ",lit)
						#if( lima<5):  sxprint("  SECOND KEPT  ",lit)


				#mpi_barrier(MPI_COMM_WORLD)
				#mpi_finalize()
				#exit()
				#for j in range(lit):
				#	 newpar[kl][2].append([int(xod2[j]),float(xod1[j])])
				if( lenass > 0):
					###print("   EICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im)#,assignments_to_cones)

					firstdirections = [[0.0,0.0] for iln in range(lit)]
					firstshifts = [0]*lit
					for iln in range(lit):
						hashparams = int(xod2[iln])
						ishift = hashparams%1000
						ipsiandiang	= hashparams/1000
						#ipsi = ipsiandiang%100000
						iang = ipsiandiang/100000
						#try:
						firstdirections[iln] = [coarse_angles[iang][0], coarse_angles[iang][1], 0.0]
						#except:
						#	sxprint(" FAILED  ",Blockdata["myid"],Tracker["keepfirst"],iln,lima,lit,hashparams,iang,xod2[:lit],xod1[:lit])
						#	mpi_finalize()
						#	exit()
						firstshifts[iln] = ishift
						#  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):
						#	ipsi = ipsiandiang%100000
						#	ishift = hashparams%1000
						#	###print("  SECONDPAR%04d  "%im,iln,hashparams, ipsi,iang,ishift)
						#	sxprint(" SECONDPAR%04d  "%im,iln, refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, shifts[ishift])
					###del xod2
					###if( Blockdata["myid"] == 18 and lima<5):   sxprint("  SECOND KEPT  ",im,lit)#,len(xod2),xod2,firstdirections,firstshifts)
					#mpi_barrier(MPI_COMM_WORLD)
					#mpi_finalize()
					#exit()

					#if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  FIFI ",firstdirections)
					#if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  GUGU ",firstshifts)
					# Find neighbors, ltabang contains positions on refang list, no psis
					###ltabang = nearest_many_full_k_projangles(refang, firstdirections, howmany = 5, sym=Tracker["constants"]["symmetry"])
					ltabang = find_nearest_k_refangles_to_many_angles(refdirs, firstdirections, Tracker["delta"], howmany = 4)
					###if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  ltabang ",ltabang)
					##mpi_barrier(MPI_COMM_WORLD)
					##mpi_finalize()
					##exit()

					# ltabang has length lit, which is the same as length as firshifts.  However, original xod2 is still available,
					#   even though it is longer than lit.


					###ltabshi = [shiftneighbors[iln] for iln in firstshifts]
					#if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  HUHU ",ltabang)
					#if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  OUOU ",ltabshi)

					#  Prepare image for chi2.
					#  We have to repeat everything from get shrink data, including shifts
					#  The number of entries is lit=len(ltabang), each proj dir has 4 neighbors, so including itself there are 5,
					#    there are 3 psis, and at most n_fine_shifts. 
					#    Not not all shifts have to be populated. We may also have repeated triplets of angles, but each can have
					#      different shifts.  If so, we have to remove duplicates from the entire set.
					lcod1 = lit*4*2*n_fine_shifts

					#cod2 = np.ndarray((lit,5,3,n_fine_shifts),dtype=int,order="C")
					#cod2.fill(-1)  #  hashparams
					cod2 = []
					#lol = 0
					for i1 in range(lit):
						hashparams = int(xod2[i1])
						ipsiandiang	= hashparams/1000
						oldiang = ipsiandiang/100000
						ipsi = ipsiandiang%100000
						ishift = hashparams%1000
						tshifts = get_shifts_neighbors(shifts, coarse_shifts[ishift])
						#if( Blockdata["myid"] == Blockdata["main_node"]): sxprint(" tshifts  ",i1,len(tshifts))
						for i2 in range(4):
							iang = ltabang[i1][i2]
							for i3 in range(2):  # psi
								itpsi = int((coarse_angles[oldiang][2] + ipsi*coarse_delta - refang[iang][2]+360.0)/Tracker["delta"])
								itpsi = (itpsi + i3)%npsi
								for i4 in range(len(tshifts)):
									cod2.append(iang*100000000 + itpsi*1000 + tshifts[i4])
									#lol += 1
									#if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  zibzi  ",i1,i2,i3,i4, lol)


					del xod1, xod2

					###if( Blockdata["myid"] == 18 and lima<5):   sxprint("  THIRD   ",len(cod2))#,cod2)
					cod2 = list(set(cod2))
					cod1 = [[q/1000,i] for i,q in enumerate(cod2)]
					cod1.sort()

					lit = len(cod1)

					cod2 = numpy.asarray([cod2[cod1[i][1]] for i in range(lit)])

					cod1 = numpy.ndarray(lit,dtype='f4',order="C")
					#cod1.fill(np.finfo(dtype='f4').min)
					cod3 = numpy.ndarray(lit,dtype='f4',order="C")
					#cod3.fill(0.0)  #  varadj

					###if( Blockdata["myid"] == 18 and lima<5): sxprint("  THIRD   ",im,lit)#,cod2)

					#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					#  Third ML part.  Now the actual chi2 distances, we compute reprojections on the fly.
					#  Make sure volprep has nxinit size
					data = [None]*nshifts
					johi = 0
					iln = 0
					prevdir = -1
					while(iln<lit):
						hashparams = cod2[iln]
						#if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  COD2   ",im,lit,iln,cod2[iln])
						ipsiandiang	= hashparams/1000
						if(ipsiandiang != prevdir):
							prevdir = ipsiandiang
							ipsi = ipsiandiang%100000
							iang = ipsiandiang/100000
							##junk = time()
							temp = sp_projection.prgl(volinit,[ refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, 0.0,0.0], 1, False)
							##eat += time()-junk
							temp.set_attr("is_complex",0)
							johi += 1
						while( ipsiandiang == cod2[iln]/1000 ):
							hashparams = cod2[iln]
							ishift = hashparams%1000
							if( data[ishift] == None ):
								xx = shifts[ishift][0]*shrink
								yy = shifts[ishift][1]*shrink
								data[ishift] = sp_fundamentals.fshift(dataml, xx, yy)
								data[ishift].set_attr("is_complex",0)
							##junk = time()
							[peak,varadj] = EMAN2_cppwrap.Util.sqednorm(data[ishift], temp, ctfa, bckgnoise)
							##eat += time()-junk
							cod1[iln] = -peak
							cod3[iln] = varadj
							iln += 1
							if(iln == lit  ):  break
						#if( Blockdata["myid"] == Blockdata["main_node"]):
						#	temp.write_image("temp.hdf")
						#	data[iln].write_image("data.hdf")
						#	ctfa.write_image("ctfa.hdf")
						#	bckgnoise.write_image("bckgnoise.hdf")
						#	exit()
						###if( Blockdata["myid"] == Blockdata["main_node"] and iln%1000 ==0):  sxprint(" progress  ",iln,time()-at)
					#if( Blockdata["myid"] == Blockdata["main_node"]):
					#	sxprint("  PROJECT   ",im,lit,johi)#,cod2)
					#	#for iln in range(lit):  sxprint("  ROPE   ",iln,cod1[iln],cod2[iln],cod3[iln])
					del data
					del dataml

					lina = numpy.argsort(cod1)
					cod1 = cod1[lina[::-1]]  # This sorts in reverse order
					cod2 = cod2[lina[::-1]]  # This sorts in reverse order
					cod3 = cod3[lina[::-1]]  # This sorts in reverse order
					cod1 -= cod1[0]
					lina = numpy.argwhere(cod1 > Tracker["constants"]["expthreshold"])
					cod1 = cod1[lina]
					cod2 = cod2[lina]
					cod3 = cod3[lina]

					###if( Blockdata["myid"] == Blockdata["main_node"]):
					###for iui in range(len(lina)):
					###	for iui in range(len(cod1)):
					###		sxprint("  MLML  ",iui,cod1[iui],exp(cod1[iui]),cod2[iui],cod3[iui])

					numpy.exp(cod1, out=cod1)
					cod1 /= numpy.sum(cod1)
					cumprob = 0.0
					for j in range(len(cod1)):
						cumprob += cod1[j]
						if(cumprob > Tracker["constants"]["ccfpercentage"]):
							lit = j+1
							break

					#  New norm is a sum of eq distances multiplied by their probabilities augmented by PW.
					norm_per_particle[im] = numpy.sum(cod1[:lit]*cod3[:lit]) + accumulatepw[im][reachpw]
					###print("   CNORMPERPARTICLE  ",Blockdata["myid"],im,norm_per_particle[im])

					for iln in range(lit):
						newpar[im][2].append([int(cod2[iln]), float(cod1[iln])])
						#  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):
						#	#print("  NEWPAR%04d  "%im,iln,newpar[im][2][-1])
						#	hashparams = newpar[im][2][iln][0]
						#	ipsiandiang	= hashparams/1000
						#	ipsi = ipsiandiang%100000
						#	iang = ipsiandiang/100000
						#	ishift = hashparams%1000
						#	#print("  NEWPAR%04d  "%im,iln,newpar[im][2][-1],hashparams,ipsi,iang,ishift)
						#	sxprint(" NEWPAR%04d  "%im,iln,newpar[im][2][iln], refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, shifts[ishift])

					###if( Blockdata["myid"] == 18 and lima<5):   sxprint("  FINALLY  ",im,lit)
					del cod1, cod2, cod3, lina
					###mpi_barrier(MPI_COMM_WORLD)
					###mpi_finalize()
					###exit()

	"""Multiline Comment24"""
	#MULTILINEMULTILINEMULTILINE 24
	#MULTILINEMULTILINEMULTILINE 24
	#MULTILINEMULTILINEMULTILINE 24
	#MULTILINEMULTILINEMULTILINE 24
	#MULTILINEMULTILINEMULTILINE 24
	#MULTILINEMULTILINEMULTILINE 24
	#MULTILINEMULTILINEMULTILINE 24
	
	
	
	#  END OF CONES
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	if( Blockdata["myid"] == Blockdata["main_node"] ):
		sp_global_def.sxprint( "  Finished projection matching   %10.1fmin"%((time.time()-at)/60.))
	at = time.time()
	#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	#  All images were processed, now to the additional calculations

	###mpi_barrier(MPI_COMM_WORLD)
	###mpi_finalize()
	###exit()


	# norm correction ---- calc the norm correction per particle
	snormcorr = 0.0
	for kl in range(nima):
		###print("   NORMPERPARTICLE  ",Blockdata["myid"],kl,norm_per_particle[kl])
		norm_per_particle[kl]	= numpy.sqrt(norm_per_particle[kl]*2.0)*oldparams[kl][7]/Tracker["avgvaradj"][procid]
		snormcorr				+= norm_per_particle[kl]
	Tracker["avgvaradj"][procid] = snormcorr
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	#  Compute avgvaradj
	Tracker["avgvaradj"][procid] = mpi.mpi_reduce( Tracker["avgvaradj"][procid], 1, mpi.MPI_FLOAT, mpi.MPI_SUM, Blockdata["main_node"], mpi.MPI_COMM_WORLD )
	if(Blockdata["myid"] == Blockdata["main_node"]):
		Tracker["avgvaradj"][procid] = float(Tracker["avgvaradj"][procid])/Tracker["nima_per_chunk"][procid]
	else:  Tracker["avgvaradj"][procid] = 0.0
	Tracker["avgvaradj"][procid] = sp_utilities.bcast_number_to_all(Tracker["avgvaradj"][procid], Blockdata["main_node"])
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	
	at = time.time()
	mpi.mpi_barrier(Blockdata["shared_comm"])

	###if Blockdata["myid"] == Blockdata["main_node"]:  print "  Finished :",time()-at
	mpi.mpi_win_free(win_sm)
	if( Tracker["nxinit"] != Tracker["nxpolar"] ):
		mpi.mpi_win_free(win_volinit)
		emnumpy4.unregister_numpy_from_emdata()
		del emnumpy4
	else:   mpi.mpi_win_free(win_vol)

	mpi.mpi_barrier(Blockdata["shared_comm"])
	emnumpy1.unregister_numpy_from_emdata()
	emnumpy2.unregister_numpy_from_emdata()
	del emnumpy1, emnumpy2

	del volinit

	mpi.mpi_barrier(Blockdata["shared_comm"])

	#print("  NORMALIZATION DONE  ",Blockdata["myid"])
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	#if( Blockdata["myid"] == Blockdata["main_node"] ):
	#	sxprint( "  Projection matching finished : %10.1fmin"%((time()-at)/60.))
	return newpar, [1.0]*nima































































































def init_Tracker_mpi(option_initvol = None):
	global Tracker, Blockdata
	# reset
	Tracker["mainiteration"]        = 0
	## nosmearing rec3d and evaluate resolution

	Tracker["previousoutputdir"]    =  None
	Blockdata["accumulatepw"]       = [None, None]
	Blockdata["bckgnoise"]          =  None
	Tracker["currentres"]		    = -1
	Tracker["fsc143"]			    = -1
	Tracker["maxfrad"]           	= -1
	Tracker["no_improvement"]    	= 0
	Tracker["no_params_changes"] 	= 0
	Tracker["large_at_Nyquist"]  	= False
	Tracker["anger"]             	= 1.e23
	Tracker["shifter"]           	= 1.e23
	Tracker["pixercutoff"]       	= 2.0
	Tracker["acc_rot"]           	= 0.0
	Tracker["acc_trans"]			= 0.0
	Tracker["avgvaradj"]			= [1.0,1.0]  # This has to be initialized to 1.0 !!
	Tracker["mainiteration"]     	= 0
	Tracker["lentop"]				= 2000
	Tracker["nxstep"]		        = 0
	try: nxinit = Tracker["nxinit"]
	except: Tracker["nxinit"] = Tracker["constants"]["nnxo"]
	Tracker["maxfrad"]              = Tracker["nxinit"]//2 #Tracker["constants"]["nnxo"]//2
	Tracker["refvol"]               = option_initvol
	Tracker["state"]                = "CONTINUATION_INITIAL"
	Tracker["constants"]["best"]    = 0
	Tracker["keepfirst"]            = -1
	Tracker["bestres"]          	= -1
	Tracker["bestres_143"]          = -1
	Tracker["directory"]         	= ""
	Tracker["previousoutputdir"] 	= ""
	Tracker["mainiteration"]     	= 0
	Tracker["nima_per_chunk"]    	= [0,0]
	try:    user_func = Tracker["constants"]["user_func"]
	except: Tracker["constants"]["user_func"] = "do_volume_mask"
	return





































































































































def compare_bckgnoise(bckgnoise1, bckgnoise2):
	ccsum = 0.0
	nx =  bckgnoise1.get_xsize()
	ny =  bckgnoise1.get_ysize()
	for i in range(ny):
		c1 = []
		c2 = []
		for j in range(1, nx):
			c1.append(bckgnoise1.get_value_at(j,i, 0))
			c2.append(bckgnoise2.get_value_at(j,i, 0))
		ccsum +=sp_statistics.pearson(c1, c2)
	return ccsum/float(ny)
	

































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































