'''1
		#  I am not sure whether what follows is correct.  This part should be recalculated upon restart
		Blockdata["accumulatepw"] = [[],[]]
		ndata = len(projdata)
		for i in range(ndata):
			if(i<first_procid):  iproc = 0 #  This points to correct procid
			else:                iproc = 1
			Blockdata["accumulatepw"][iproc].append([0.0]*200)
		'''
"""2
		#  Inverted Gaussian mask
		invg = model_gauss(Tracker["constants"]["radius"],nx,nx)
		invg /= invg[nx//2,nx//2]
		invg = model_blank(nx,nx,1,1.0) - invg
		"""
'''3
		if(myid == 0):  ndata = EMUtil.get_image_count(partstack)
		else:           ndata = 0
		ndata = bcast_number_to_all(ndata)
		if( ndata < Blockdata["nproc"]):
			if(myid<ndata):
				image_start = myid
				image_end   = myid+1
			else:
				image_start = 0
				image_end   = 1
		else:
			image_start, image_end = MPI_start_end(ndata, Blockdata["nproc"], myid)
		#data = EMData.read_images(stack, range(image_start, image_end))
		if(myid == 0):
			params = read_text_row( paramsname )
			params = [params[i][j]  for i in range(len(params))   for j in range(5)]
		else:           params = [0.0]*(5*ndata)
		params = bcast_list_to_all(params, myid, source_node=Blockdata["main_node"])
		params = [[params[i*5+j] for j in range(5)] for i in range(ndata)]
		'''
'''4
			if doac:
				if(i<first_procid):  iproc = 0 #  This points to correct procid
				else:                iproc = 1
				Blockdata["accumulatepw"][iproc].append(sig[nv:]+[0.0])  # add zero at the end so for the full size nothing is added.
			'''
"""5
				for k in range(6,nv):
					tsd.set_value_at(k,i,1.0/(tsd.get_value_at(k,i)/tocp[i]))  # Already inverted
				qt = tsd.get_value_at(6,i)
				for k in range(1,6):
					tsd.set_value_at(k,i,qt)
				"""
'''6
	particles_on_node = []
	parms_on_node     = []
	for i in range( group_start, group_end ):
		particles_on_node += lpartids[group_reference[i][2]:group_reference[i][3]+1]  #  +1 is on account of python idiosyncrasies
		parms_on_node     += partstack[group_reference[i][2]:group_reference[i][3]+1]


	Blockdata["nima_per_cpu"][procid] = len(particles_on_node)
	#ZZprint("groups_on_thread  ",Blockdata["myid"],procid, Tracker["groups_on_thread"][procid])
	#ZZprint("  particles  ",Blockdata["myid"],Blockdata["nima_per_cpu"][procid],len(parms_on_node))
	'''
"""7
            17            28            57            84    5
            18            14            85            98    6
            19            15            99           113    7
            25            20           114           133    8
            29             9           134           142    9

	"""
"""8
	This function will read from stack a subset of images specified in partids
	   and assign to them parameters from partstack with optional CTF application and shifting of the data.
	So, the lengths of partids and partstack are the same.
	  The read data is properly distributed among MPI threads.
	
	Flow of data:
	1. Read images, if there is enough memory, keep them as original_data.
	2. Read current params
	3.  Apply shift
	4.  Normalize outside of the radius
	5.  Do noise substitution and cosine mask.  (Optional?)
	6.  Shrink data.
	7.  Apply CTF.
	
	"""
'''9
		if preshift:
			data[im] = fshift(original_data[im], sx, sy)
			sx = 0.0
			sy = 0.0
		'''
"""10
	if(Tracker["delta"] == 15.0):  refang = read_text_row("refang15.txt")
	elif(Tracker["delta"] == 7.5):  refang = read_text_row("refang7p5.txt")
	elif(Tracker["delta"] == 3.75):  refang = read_text_row("refang3p75.txt")
	elif(Tracker["delta"] == 1.875):  refang = read_text_row("refang1p875.txt")
	elif(Tracker["delta"] == 0.9375):  refang = read_text_row("refang0p9375.txt")
	elif(Tracker["delta"] == 0.46875):  refang = read_text_row("refang0p46875.txt")
	"""
'''11
			if( method == "DIRECT" ):
				#dss = fshift(ds, xx+sxs, yy+sys)
				dss = fshift(ds, xx+sxs, yy+sys)
				dss.set_attr("is_complex",0)
			else:
				dss = fft(fshift(ds, x+sxs, yy+sys))
				dss,kb = prepi(dss)
			'''
"""12
	tvol, tweight, trol = recons3d_4nnstruct_MPI(myid = Blockdata["subgroup_myid"], main_node = Blockdata["nodes"][procid], prjlist = data, \
											paramstructure = newparams, refang = refang, delta = Tracker["delta"], CTF = Tracker["constants"]["CTF"],\
											upweighted = False, mpi_comm = mpi_comm, \
											target_size = (2*Tracker["nxinit"]+3), avgnorm = Tracker["avgvaradj"][procid], norm_per_particle = norm_per_particle)
	"""
'''13
def getangc5(p1,p2):
	pass#IMPORTIMPORTIMPORT from utilities import getfvec
	pass#IMPORTIMPORTIMPORT from math import acos, degrees
	n1 = getfvec(p1[0],p1[1])
	n2 = getfvec(p1[0]+72.,p1[1])
	n3 = getfvec(p1[0]-72.,p1[1])
	n4 = getfvec(p2[0],p2[1])
	return degrees(min(acos(max(min((n1[0]*n4[0]+n1[1]*n4[1]+n1[2]*n4[2]),1.0),-1.0)),\
	acos(max(min((n2[0]*n4[0]+n2[1]*n4[1]+n2[2]*n4[2]),1.0),-1.0)),acos(max(min((n3[0]*n4[0]+n3[1]*n4[1]+n3[2]*n4[2]),1.0),-1.0))))

def difangc5(n4,p1):
	pass#IMPORTIMPORTIMPORT from utilities import getfvec
	pass#IMPORTIMPORTIMPORT from math import acos, degrees
	n1 = getfvec(p1[0],p1[1])
	n2 = getfvec(p1[0]+72.,p1[1])
	n3 = getfvec(p1[0]-72.,p1[1])
	#n4 = getfvec(p2[0],p2[1])
	return degrees(min(acos(max(min((n1[0]*n4[0]+n1[1]*n4[1]+n1[2]*n4[2]),1.0),-1.0)),\
	acos(max(min((n2[0]*n4[0]+n2[1]*n4[1]+n2[2]*n4[2]),1.0),-1.0)),acos(max(min((n3[0]*n4[0]+n3[1]*n4[1]+n3[2]*n4[2]),1.0),-1.0))))
'''
'''14
def XNumrinit_local(first_ring, last_ring, skip=1, mode="F"):
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
	MAXFFT = 32768
	pass#IMPORTIMPORTIMPORT from math import pi
	skip = 1

	if (mode == 'f' or mode == 'F'): dpi = 2*pi
	else:                            dpi = pi
	numr = []
	lcirc = 1
	for k in range(first_ring, last_ring+1, skip):
		numr.append(k)
		jp = int(dpi * k+0.5)
		ip = 2**(log2(jp)+1)  # two times oversample each ring
		if (k+skip <= last_ring and jp > ip+ip//2): ip=min(MAXFFT,2*ip)
		if (k+skip  > last_ring and jp > ip+ip//5): ip=min(MAXFFT,2*ip)

		numr.append(lcirc)
		numr.append(ip)
		lcirc += ip

	return  numr
'''
'''15
	if(Blockdata["myid"] <3):
		for kl in range(0,ndat,ndat/2):
			for m in range(0,len(data[kl]),len(data[kl])/3):  print(" DNORM  ",Blockdata["myid"],kl,m, Util.innerproduct(data[kl][m],data[kl][m],mask))
	'''
'''16
	if(Blockdata["myid"] <3):
		for kl in range(0,ndat,ndat/2):
			for m in range(0,len(data[kl]),len(data[kl])/3):  print(" DNORM  ",Blockdata["myid"],kl,m, Util.innerproduct(data[kl][m],data[kl][m],mask))
	'''
"""17
	if Blockdata["bckgnoise"] :
		oneover = []
		nnx = Blockdata["bckgnoise"][0].get_xsize()
		for i in range(len(Blockdata["bckgnoise"])):
			temp = [0.0]*nnx
			for k in range(nnx):
				if( Blockdata["bckgnoise"][i].get_value_at(k) > 0.0):  temp[k] = 1.0/sqrt(Blockdata["bckgnoise"][i].get_value_at(k))
			oneover.append(temp)
		del temp

	accumulatepw = [None]*nima
	norm_per_particle = [None]*nima
	"""
"""18
		if Blockdata["bckgnoise"]:
			temp = Blockdata["bckgnoise"][dataimage.get_attr("particle_group")]
			bckgn = Util.unroll1dpw(Tracker["nxinit"], [temp[i] for i in range(temp.get_xsize())])
		else:
			bckgn = Util.unroll1dpw(Tracker["nxinit"], [1.0]*600)
		bckgnoise = bckgn.copy()
		for j in range(Tracker["nxinit"]//2+1,Tracker["nxinit"]):  bckgn[0,j] = bckgn[0,Tracker["nxinit"]-j]
		"""
"""19
			lina = np.argsort(xod1)
			xod1 = xod1[lina[::-1]]  # This sorts in reverse order
			xod2 = xod2[lina[::-1]]  # This sorts in reverse order
			np.exp(xod1, out=xod1)
			xod1 /= np.sum(xod1)
			cumprob = 0.0
			lit = len(xod1)
			for j in range(len(xod1)):
				cumprob += xod1[j]
				if(cumprob > Tracker["constants"]["ccfpercentage"]):
					lit = j+1
					break
			"""
"""20
			xod1 -= xod1[0]

			lina = np.argwhere(xod1 > Tracker["constants"]["expthreshold"])
			xod1 = xod1[lina]
			xod2 = xod2[lina]
			np.exp(xod1, out=xod1)
			xod1 /= np.sum(xod1)
			cumprob = 0.0
			lit = len(xod1)
			for j in range(len(xod1)):
				cumprob += xod1[j]
				if(cumprob > Tracker["constants"]["ccfpercentage"]):
					lit = j+1
					break
			"""
"""21
		np.exp(cod1, out=cod1)
		cod1 /= np.sum(cod1)
		cumprob = 0.0
		for j in range(len(cod1)):
			cumprob += cod1[j]
			if(cumprob > Tracker["constants"]["ccfpercentage"]):
				lit = j+1
				break
		"""
"""22
	# norm correction ---- calc the norm correction per particle
	snormcorr = 0.0
	for kl in range(nima):
		norm_per_particle[kl] = sqrt(norm_per_particle[kl]*2.0)*oldparams[kl][7]/Tracker["avgvaradj"][procid]
		snormcorr            += norm_per_particle[kl]
	Tracker["avgvaradj"][procid] = snormcorr
	mpi_barrier(MPI_COMM_WORLD)
	#  Compute avgvaradj
	Tracker["avgvaradj"][procid] = mpi_reduce( Tracker["avgvaradj"][procid], 1, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD )
	if(Blockdata["myid"] == Blockdata["main_node"]):
		Tracker["avgvaradj"][procid] = float(Tracker["avgvaradj"][procid])/Tracker["nima_per_chunk"][procid]
	else:  Tracker["avgvaradj"][procid] = 0.0
	Tracker["avgvaradj"][procid] = bcast_number_to_all(Tracker["avgvaradj"][procid], Blockdata["main_node"])
	mpi_barrier(MPI_COMM_WORLD)

	#  Compute statistics of smear -----------------
	smax = -1000000
	smin = 1000000
	sava = 0.0
	svar = 0.0
	snum = 0
	for kl in range(nima):
		j = len(newpar[kl][2])
		snum += 1
		sava += float(j)
		svar += j*float(j)
		smax = max(smax, j)
		smin = min(smin, j)
	snum = mpi_reduce(snum, 1, MPI_INT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	sava = mpi_reduce(sava, 1, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	svar = mpi_reduce(svar, 1, MPI_FLOAT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	smax = mpi_reduce(smax, 1, MPI_INT, MPI_MAX, Blockdata["main_node"], MPI_COMM_WORLD)
	smin = mpi_reduce(smin, 1, MPI_INT, MPI_MIN, Blockdata["main_node"], MPI_COMM_WORLD)
	if( Blockdata["myid"] == 0 ):
		pass#IMPORTIMPORTIMPORT from math import sqrt
		sava = float(sava)/snum
		svar = sqrt(max(0.0,(float(svar) - snum*sava**2)/(snum -1)))
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print(line, "Smear stat  (number of images, ave, sumsq, min, max)):  %7d    %12.3g   %12.3g  %7d  %7d"%(snum,sava,svar,smin,smax))
	"""
"""23
	Number of images :     82     410   20.0%          1.3min          0.8min
	Number of images :    164     410   40.0%          3.0min          1.8min
	Number of images :    246     410   60.0%          5.7min          3.1min
	Number of images :    328     410   80.0%          8.8min          4.4min
	#  Projection and equdist take 50% time, so on the whole run of the program one could
	#    reduce time from 311 to 233, (6h to 4h) if this time was totally eliminated.
	"""
"""24
	Number of images :     82     410   20.0%          1.3min          0.8min
	Number of images :    164     410   40.0%          3.0min          1.8min
	Number of images :    246     410   60.0%          5.7min          3.1min
	Number of images :    328     410   80.0%          8.8min          4.4min
	#  Projection and equdist take 50% time, so on the whole run of the program one could
	#    reduce time from 311 to 233, (6h to 4h) if this time was totally eliminated.
	"""
"""25
	Number of images :     82     410   20.0%          1.3min          0.8min
	Number of images :    164     410   40.0%          3.0min          1.8min
	Number of images :    246     410   60.0%          5.7min          3.1min
	Number of images :    328     410   80.0%          8.8min          4.4min
	#  Projection and equdist take 50% time, so on the whole run of the program one could
	#    reduce time from 311 to 233, (6h to 4h) if this time was totally eliminated.
	"""
"""26
	Number of images :     82     410   20.0%          1.3min          0.8min
	Number of images :    164     410   40.0%          3.0min          1.8min
	Number of images :    246     410   60.0%          5.7min          3.1min
	Number of images :    328     410   80.0%          8.8min          4.4min
	#  Projection and equdist take 50% time, so on the whole run of the program one could
	#    reduce time from 311 to 233, (6h to 4h) if this time was totally eliminated.
	"""
"""27
			if(Blockdata["myid"] == Blockdata["main_node"]):
				fff = read_text_file(os.path.join(initdir,"driver_%03d.txt"%(Tracker["mainiteration"])))
			else:
				fff = []
			fff = bcast_list_to_all(fff, Blockdata["myid"], source_node=Blockdata["main_node"])
			keepgoing = AI_continuation(fff)
			"""
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
	pass#IMPORTIMPORTIMPORT from fundamentals 	import prepi
	pass#IMPORTIMPORTIMPORT from morphology 	import ctf_img_real
	#  Data is NOT CTF-applied.
	#  Data is shrank, in Fourier format
	data = [[] for i in range(len(projdata))]
	if Tracker["constants"]["CTF"]:
		nx = projdata[0].get_ysize()
		ctfs = [ morphology.ctf_img_real(nx, q.get_attr('ctf')) for q in projdata ]
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
			dss = fundamentals.fshift(ds, xx, yy)
			dss.set_attr("is_complex",0)
			"""Multiline Comment10"""
			data[kl].append(dss)
		data[kl][0].set_attr("particle_group",particle_group)  #  Pass group number only in the first shifted copy.
		del projdata[kl]
	return data, ctfs, bckgnoise

def Xali3D_direct_ccc(data, refang, shifts, ctfs = None, bckgnoise = None, kb3D = None):
	global Tracker, Blockdata
	pass#IMPORTIMPORTIMPORT from projection 	import prgs,prgl
	pass#IMPORTIMPORTIMPORT from fundamentals 	import fft
	pass#IMPORTIMPORTIMPORT from utilities 		import wrap_mpi_gatherv
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
	if(Blockdata["myid"] == 0):  print("  ENTERING Xali buffered exhaustive CCC  ")
	npsi = int(360./Tracker["delta"])
	nang = len(refang)
	ndat = len(data)

	ny = data[0][0].get_ysize()
	mask = EMAN2_cppwrap.Util.unrollmask(ny)
	nxt = 2*(mask.get_xsize())

	"""Multiline Comment14"""

	if Tracker["mainiteration"]>1 :
		#first = True
		if Tracker["constants"]["CTF"] :
			for kl in range(ndat):
				for im in range(len(shifts)):
					EMAN2_cppwrap.Util.mulclreal(data[kl][im], ctfs[kl])
		del ctfs
		if bckgnoise:  #  This should be a flag to activate sharpening during refinement as bckgnoise is always present (for 3D later)
			for kl in range(ndat):
				temp = EMAN2_cppwrap.Util.unroll1dpw(ny, bckgnoise[kl])
				for im in range(len(shifts)):
					EMAN2_cppwrap.Util.mulclreal(data[kl][im], temp)
			del bckgnoise
	#else:  first = False

	disp_unit = np.dtype("f4").itemsize

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

		nxvol = utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nyvol = utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nzvol = utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

		if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = utilities.model_blank( nxvol,nyvol, nzvol)

		utilities.bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

		odo = projection.prep_vol( odo, npad = 2, interpolation_method = 1 )
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

	orgsizevol = utilities.bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
	nxvol = utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
	nyvol = utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
	nzvol = utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

	win_vol, base_vol  = mpi.mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
	sizevol = orgsizevol
	if( Blockdata["myid_on_node"] != 0 ):
		base_vol, = mpi.mpi_win_shared_query(win_vol, mpi.MPI_PROC_NULL)

	volbuf = np.frombuffer(np.core.multiarray.int_asbuffer(base_vol, sizevol*disp_unit), dtype = 'f4')
	volbuf = volbuf.reshape(nzvol, nyvol, nxvol)
	if( Blockdata["myid_on_node"] == 0 ):
		np.copyto(volbuf,ndo)
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

	buffer = np.frombuffer(np.core.multiarray.int_asbuffer(base_ptr, size*disp_unit), dtype = 'f4')
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
			print( "  Angle :%7d   %5d  %5.1f"%(i,ndat,float(i)/float(nang)*100.) + "%" +"   %10.1fmin"%((time.time()-at)/60.))

		if(i%Blockdata["no_of_processes_per_group"] == 0 ):  #  Time to fill up the buffer
			for itemp in range(i, min(i+Blockdata["no_of_processes_per_group"], nang)):
				if( itemp-i == Blockdata["myid_on_node"]):
					for j in range(npsi):
						psi = (refang[i][2] + j*Tracker["delta"])%360.0
						###if kb3D:  rtemp = fft(prgs(volprep, kb3D, [refang[i][0],refang[i][1],psi, 0.0,0.0]))
						###else:     
						temp = projection.prgl(volprep,[ refang[itemp][0],refang[itemp][1],psi, 0.0,0.0], 1, False)
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
			img_buffer = np.frombuffer(np.core.multiarray.int_asbuffer(pointer_location, size_of_one_image*disp_unit), dtype = 'f4')
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
	pass#IMPORTIMPORTIMPORT from projection 	import prgs,prgl
	pass#IMPORTIMPORTIMPORT from fundamentals 	import fft
	pass#IMPORTIMPORTIMPORT from utilities 		import wrap_mpi_gatherv
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
		utilities.write_text_row(coarse_angles,"coarse_angles.txt")

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
	mask = EMAN2_cppwrap.Util.unrollmask(ny)
	nxt = 2*(mask.get_xsize())

	"""Multiline Comment15"""

	if Tracker["mainiteration"]>1 :
		#first = True
		if Tracker["constants"]["CTF"] :
			for kl in range(ndat):
				for im in range(n_coarse_shifts+nshifts):
					EMAN2_cppwrap.Util.mulclreal(data[kl][im], ctfs[kl])
		del ctfs
		if bckgnoise:  #  This should be a flag to activate sharpening during refinement as bckgnoise is always present (for 3D later)
			for kl in range(ndat):
				temp = EMAN2_cppwrap.Util.unroll1dpw(ny, bckgnoise[kl])
				for im in range(n_coarse_shifts+nshifts):
					EMAN2_cppwrap.Util.mulclreal(data[kl][im], temp)
			del bckgnoise
	#else:  first = False

	disp_unit = np.dtype("f4").itemsize

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

		nxvol = utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nyvol = utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nzvol = utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

		if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = utilities.model_blank( nxvol,nyvol, nzvol)

		utilities.bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

		odo = projection.prep_vol( odo, npad = 2, interpolation_method = 1 )
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

	orgsizevol = utilities.bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
	nxvol = utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
	nyvol = utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
	nzvol = utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

	win_vol, base_vol  = mpi.mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
	sizevol = orgsizevol
	if( Blockdata["myid_on_node"] != 0 ):
		base_vol, = mpi.mpi_win_shared_query(win_vol, mpi.MPI_PROC_NULL)

	volbuf = np.frombuffer(np.core.multiarray.int_asbuffer(base_vol, sizevol*disp_unit), dtype = 'f4')
	volbuf = volbuf.reshape(nzvol, nyvol, nxvol)
	if( Blockdata["myid_on_node"] == 0 ):
		np.copyto(volbuf,ndo)
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
	if( Blockdata["myid"] == 0 ):  print("  BIGBUFFER  ",float(orgsize)/1.e9,"GB")
	if( Blockdata["myid_on_node"] == 0 ): size = orgsize
	else:  size = 0

	win_sm, base_ptr  = mpi.mpi_win_allocate_shared( size*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
	size = orgsize
	if( Blockdata["myid_on_node"] != 0 ):
		base_ptr, = mpi.mpi_win_shared_query(win_sm, mpi.MPI_PROC_NULL)

	buffer = np.frombuffer(np.core.multiarray.int_asbuffer(base_ptr, size*disp_unit), dtype = 'f4')
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
			print( "  Angle :%7d   %5d  %5.1f"%(i,ndat,float(i)/float(n_coarse_ang)*100.) + "%" +"   %10.1fmin"%((time.time()-at)/60.))

		if(i%Blockdata["no_of_processes_per_group"] == 0 ):  #  Time to fill up the buffer
			for itemp in range(i, min(i+Blockdata["no_of_processes_per_group"], n_coarse_ang)):
				if( itemp-i == Blockdata["myid_on_node"]):
					for j in range(n_coarse_psi):
						psi = (coarse_angles[i][2] + j*coarse_delta)%360.0
						###if kb3D:  rtemp = fft(prgs(volprep, kb3D, [refang[i][0],refang[i][1],psi, 0.0,0.0]))
						###else:
						temp = projection.prgl(volprep,[ coarse_angles[itemp][0],coarse_angles[itemp][1],psi, 0.0,0.0], 1, False)
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
			img_buffer = np.frombuffer(np.core.multiarray.int_asbuffer(pointer_location, size_of_one_image*disp_unit), dtype = 'f4')
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
	if Blockdata["myid"] == Blockdata["main_node"]:   print("  COARSE SEARCHES DONE  ","   %10.1fmin"%((time.time()-at)/60.))
	at = time.time()
	for itemp in range(0, min(Blockdata["no_of_processes_per_group"], nang)):
		if( itemp == Blockdata["myid_on_node"]):
			for j in range(npsi):
				psi = (refang[i][2] + j*Tracker["delta"])%360.0
				temp = projection.prgl(volprep,[ refang[itemp][0],refang[itemp][1],psi, 0.0,0.0], 1, False)
				EMAN2_cppwrap.Util.mulclreal(temp, mask)
				nrmref = numpy.sqrt(EMAN2_cppwrap.Util.innerproduct(temp, temp, None))
				temp.set_attr("is_complex",0)
				EMAN2_cppwrap.Util.mul_scalar(temp, 1.0/nrmref)
				bigbuffer.insert_clip(temp,(0,0,itemp*npsi+j))

	mpi.mpi_barrier(Blockdata["shared_comm"])
	if( Blockdata["myid"] == Blockdata["main_node"] ):   print("  REFPROJ DONE  ","   %10.1fmin"%((time.time()-at)/60.))
	at = time.time()

	opar = []
	Blockdata["angle_set"] = refang
	#  Extract normals from rotation matrices
	refdirs = utilities.angles_to_normals(refang)

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
				img_buffer = np.frombuffer(np.core.multiarray.int_asbuffer(pointer_location, size_of_one_image*disp_unit), dtype = 'f4')
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
	if Blockdata["myid"] == Blockdata["main_node"]:   print("  FINE SEARCH DONE  ","   %10.1fmin"%((time.time()-at)/60.))

	#print("  NORMALIZATION DONE  ",Blockdata["myid"])
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	opar = utilities.wrap_mpi_gatherv(opar, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	if( Blockdata["myid"] == Blockdata["main_node"] ):   utilities.write_text_row(opar,"opar.txt")
	#print("  ALL DONE  ",Blockdata["myid"])
	return newpar

def ali3D_local_polar_ccc(refang, shifts, coarse_angles, coarse_shifts, procid, original_data = None, oldparams = None, \
					preshift = False, apply_mask = True, nonorm = False, applyctf = True):
	global Tracker, Blockdata
	pass#IMPORTIMPORTIMPORT from projection 	import prgs,prgl
	pass#IMPORTIMPORTIMPORT from fundamentals 	import fft
	pass#IMPORTIMPORTIMPORT from utilities 		import wrap_mpi_gatherv
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
	pass#IMPORTIMPORTIMPORT from utilities    import get_im, model_gauss_noise, set_params_proj, get_params_proj
	pass#IMPORTIMPORTIMPORT from fundamentals import fdecimate, fshift, fft
	pass#IMPORTIMPORTIMPORT from filter       import filt_ctf, filt_table
	pass#IMPORTIMPORTIMPORT from applications import MPI_start_end
	pass#IMPORTIMPORTIMPORT from math         import sqrt

	if(Blockdata["myid"] == Blockdata["main_node"]):
		print_dict(Tracker,"PROJECTION MATCHING parameters of buffered local polar ccc")

	at = time.time()
	shrinkage = float(Tracker["nxpolar"])/float(Tracker["constants"]["nnxo"])
	shrink = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
	mode = "F"
	numr = Numrinit_local(1, radius, 1, mode)
	wr = alignment.ringwe(numr, mode)
	cnx = float(Tracker["nxpolar"]//2 + 1)

	##if(Blockdata["myid"] == Blockdata["main_node"]):  print("  RADIUS IN POLAR  ",radius,numr[-1])
	#  FINE SEARCH CONSTANTS
	nang = len(refang)
	ac_fine = numpy.cos(numpy.radians(Tracker["an"]))
	npsi = int(360./Tracker["delta"])
	mpsi = 2
	c_fine_psi = mpsi//2
	nshifts = len(shifts)
	n_fine_shifts = 4

	nima = len(original_data)
	mask = EMAN2_cppwrap.Util.unrollmask(Tracker["nxinit"])
	for j in range(Tracker["nxinit"]//2,Tracker["nxinit"]):  mask[0,j]=1.0
	mask2D = utilities.model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
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

	###if(Blockdata["myid"] == Blockdata["main_node"]): print("   TRETR  ",Tracker["constants"]["nnxo"],Tracker["nxinit"],reachpw,n_coarse_ang,coarse_delta,n_coarse_psi,m_coarse_psi,c_coarse_psi,n_coarse_shifts)

	#if(Blockdata["myid"] == Blockdata["main_node"]):
	#	print( original_data[0].get_attr("identifier") )
	#	print( original_data[1].get_attr("identifier") )
	
	disp_unit = np.dtype("f4").itemsize

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

		nxvol = utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nyvol = utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nzvol = utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

		if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = utilities.model_blank( nxvol,nyvol, nzvol)

		utilities.bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

		odo = projection.prep_vol( odo, npad = 2, interpolation_method = 1 )
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

	orgsizevol = utilities.bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
	nxvol = utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
	nyvol = utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
	nzvol = utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

	win_vol, base_vol  = mpi.mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
	sizevol = orgsizevol
	if( Blockdata["myid_on_node"] != 0 ):
		base_vol, = mpi.mpi_win_shared_query(win_vol, mpi.MPI_PROC_NULL)

	volbuf = np.frombuffer(np.core.multiarray.int_asbuffer(base_vol, sizevol*disp_unit), dtype = 'f4')
	volbuf = volbuf.reshape(nzvol, nyvol, nxvol)
	if( Blockdata["myid_on_node"] == 0 ):
		np.copyto(volbuf,ndo)
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

			nxvol = utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
			nyvol = utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
			nzvol = utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

			if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = utilities.model_blank( nxvol,nyvol, nzvol)

			utilities.bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

			odo = projection.prep_vol( odo, npad = 2, interpolation_method = 1 )
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

		orgsizevol = utilities.bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
		nxvol = utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
		nyvol = utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
		nzvol = utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

		win_volinit, base_volinit  = mpi.mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
		sizevol = orgsizevol
		if( Blockdata["myid_on_node"] != 0 ):
			base_volinit, = mpi.mpi_win_shared_query(win_volinit, mpi.MPI_PROC_NULL)

		volbufinit = np.frombuffer(np.core.multiarray.int_asbuffer(base_volinit, sizevol*disp_unit), dtype = 'f4')
		volbufinit = volbufinit.reshape(nzvol, nyvol, nxvol)
		if( Blockdata["myid_on_node"] == 0 ):
			np.copyto(volbufinit,ndoinit)
			del odo,ndoinit

		#volinit = EMNumPy.assign_numpy_to_emdata(volbufinit)
		emnumpy4 = EMAN2_cppwrap.EMNumPy()
		volinit = emnumpy4.register_numpy_to_emdata(volbufinit)
		volinit.set_attr_dict({'is_complex':1,  'is_complex_x': 0, 'is_fftodd': 0, 'is_fftpad': 1, 'is_shuffled': 1,'npad': 2})
		if( Blockdata["myid_on_node"] == 0 ):  volinit.update()
		mpi.mpi_barrier(Blockdata["shared_comm"])
		###if( Blockdata["myid"] < 10 ):
		###	print(" VOLPREP  ",volinit[0],volinit[1],info(volinit))
	else: volinit = volprep
	#  End of replaced volprep

    ### local search is plugged in from here
	#  START CONES
	#  This has to be systematically done per node
	#
	crefim = EMAN2_cppwrap.Util.Polar2Dm(utilities.model_blank(Tracker["nxpolar"],Tracker["nxpolar"]), cnx, cnx, numr, mode)
	size_of_one_image = crefim.get_xsize()
	#  We will assume half of the memory is available.  We will do it betteer later.
	numberofrefs_inmem = int(Tracker["constants"]["memory_per_node"]/4/((size_of_one_image*disp_unit)/1.0e9))
	####if( Blockdata["myid_on_node"] == 0  ):  print( " MEMEST ", n_coarse_ang,numberofrefs_inmem)
	#  number of references that will fit into one mode
	normals_set = utilities.angles_to_normals(coarse_angles)
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
			assignments_of_refangles_to_cones[i] = utilities.wrap_mpi_gatherv(assignments_of_refangles_to_cones[i], 0, Blockdata["shared_comm"])
			doit = 1
			if( Blockdata["myid_on_node"] == 0 ):
				assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i]))
				#print(  " assignments_of_refangles_to_cones on zero ",i,len(assignments_of_refangles_to_cones[i]) )#,assignments_of_refangles_to_cones[i]

		for i,q in enumerate(assignments_of_refangles_to_cones):
			###if( Blockdata["myid_on_node"] == 0 ): print( " which refangles belong to which cone",Blockdata["color"],Blockdata["myid"],i,len(q) )#,q
			assignments_of_refangles_to_cones[i] = utilities.wrap_mpi_bcast(q, 0, Blockdata["shared_comm"])

	else:

		angledirs = utilities.angles_to_normals([u1[:3] for u1 in oldparams])

		number_of_cones = max(2, int(n_coarse_ang/numberofrefs_inmem*1.5 + 0.5))
		###if Blockdata["myid_on_node"] == 0:  print( " LENGTHS  ",Blockdata["color"],nima,n_coarse_ang, numberofrefs_inmem,number_of_cones)
		cont = True
		while  cont :
			#  Translate number of cones to cone_delta
			cone_delta, number_of_cones = number_of_cones_to_delta(number_of_cones)
			#  Generate cone_angles
			##if Blockdata["myid"] == 0:  print( "  WHILE  ",number_of_cones, cone_delta)
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

			#if Blockdata["myid"] == 0:  print(  "  number of cones ",number_of_cones,cone_delta, len(cone_angles))
			assert(number_of_cones == len(cone_angles) )

			conedirs = utilities.angles_to_normals(Blockdata["symclass"].symmetry_neighbors(cone_angles))
			neighbors = len(conedirs)/len(cone_angles)  #  Symmetry related neighbors
			#if Blockdata["myid"] == 0:  print(  "  neighbors  ",Blockdata["myid"],neighbors, cone_angles)
			#  assign data directions to cone_angles
			assignments_to_cones = utilities.assign_projdirs_f(angledirs, conedirs, neighbors)
			###print(  " assignments_to_cones ",Blockdata["myid"],len(assignments_to_cones),[len(q) for q in assignments_to_cones],assignments_to_cones[0])
			#  the above should have length of refdirs and each should have indexes of data that belong to this cone
			del conedirs
			#print "assignments_to_cones ",assignments_to_cones
			#  For each cone we have to find which refangles are needed to do the matching
			assignments_of_refangles_to_cones = [[] for i in range(len(assignments_to_cones))]
			assignments_of_refangles_to_angles = [[] for i in range(nima)]  # for each myid separately, these are angles on this myid


			for i,q in enumerate(assignments_to_cones):
				#  find assignments of refdirs to this cone within an of each angledirs, so we need symmetry neighbors set of refdirs.
				#if Blockdata["myid"] == 0:  print( "in loop ", Blockdata["myid"],i,len(q),q)

				if(len(q) == 0):
					# empty assignment, on a given CPU there are no images assigned to a given cone
					assignments_of_refangles_to_cones[i] = [-1]
				else:
					for m in q:
						#print " m ",m,len(angles)

						assignments_of_refangles_to_angles[m] = find_assignments_of_refangles_to_angles(normals_set, oldparams[m], Tracker["an"])
						#if Blockdata["myid"] == 0:  print( "assignments_of_refangles_to_angles[m] ", Blockdata["color"],i,m,assignments_of_refangles_to_angles[m])
						assignments_of_refangles_to_cones[i].extend(assignments_of_refangles_to_angles[m])

					#if( Blockdata["myid"] == 19 ):  print(  " doit0 ",Blockdata["myid"], i,assignments_of_refangles_to_cones[i],q,assignments_to_cones)
					assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i]))
				###if Blockdata["myid"] == 19:  print(  " assignments_of_refangles_to_cones on myid ",Blockdata["color"],Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]), assignments_of_refangles_to_cones[i] )
				assignments_of_refangles_to_cones[i] = utilities.wrap_mpi_gatherv(assignments_of_refangles_to_cones[i], 0, Blockdata["shared_comm"])
				###if Blockdata["myid"] == 0:  print(  " assignments_of_refangles_to_cones gatherv",Blockdata["color"],Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]) )
				doit = 1
				if( Blockdata["myid_on_node"] == 0 ):
					assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i])-set([-1]))
					###if( Blockdata["myid_on_node"] == 0 ):  print(  " assignments_of_refangles_to_cones on zero ",i,len(assignments_of_refangles_to_cones[i]))
					#  POSSIBLE PROBLEM - IT IS POSSIBLE FOR A GIVEN CONE TO HAVE NO REFANGLES AND THUS BE EMPTY
					if( len(assignments_of_refangles_to_cones[i]) > numberofrefs_inmem ):
						number_of_cones = int(number_of_cones*1.25)
						#print(  " increased number_of_cones ",i,number_of_cones )
						doit = 0
				doit = utilities.bcast_number_to_all(doit, source_node = 0)
				number_of_cones = utilities.bcast_number_to_all(number_of_cones, source_node = 0, mpi_comm = Blockdata["shared_comm"] )
				###if( Blockdata["myid"] == 19 ):  print(  " doit ",Blockdata["myid"], i,doit ,assignments_of_refangles_to_cones[i],assignments_to_cones)
				if( doit == 0 ):  break

			if( doit == 1 ):
				cont = False


		for i,q in enumerate(assignments_of_refangles_to_cones):
			###if( Blockdata["myid_on_node"] == 0 ): print( " which refangles belong to which cone IOUT ",Blockdata["color"],Blockdata["myid"],i,len(q))#,q
			assignments_of_refangles_to_cones[i] = utilities.wrap_mpi_bcast(q, 0, Blockdata["shared_comm"])
			###if( Blockdata["myid_on_node"] == 0 ): print( " which refangles belong to which cone XOUT ",Blockdata["color"],Blockdata["myid"],i,len(q))#,q
			#if( myid == 1 ):
			#	print " which refangles belong to which cone",myid,i,len(assignments_of_refangles_to_cones[i])#,q

	if( Blockdata["myid"] == 0  ):  print( " number_of_cones: ",number_of_cones)
	
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
	###	print( " which refangles belong to which cone OUT ",Blockdata["color"],Blockdata["myid_on_node"],Blockdata["myid"],i,numberofrefs_inmem,nlocal_angles,len(q))#,q

	#  BIG BUFFER
	lenbigbuf = numberofrefs_inmem  #MAXIMUM NUMBER OF REFERENCE IN ONE OF THE CONES
	orgsize = lenbigbuf*size_of_one_image #  This is number of projections to be computed simultaneously times their size

	if( Blockdata["myid_on_node"] == 0 ): size = orgsize
	else:  size = 0

	win_sm, base_ptr  = mpi.mpi_win_allocate_shared( size*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
	size = orgsize
	if( Blockdata["myid_on_node"] != 0 ):
		base_ptr, = mpi.mpi_win_shared_query(win_sm, mpi.MPI_PROC_NULL)

	buffer = np.frombuffer(np.core.multiarray.int_asbuffer(base_ptr, size*disp_unit), dtype = 'f4')
	buffer = buffer.reshape(lenbigbuf, size_of_one_image)
	#bigbuffer = EMNumPy.assign_numpy_to_emdata(buffer)

	emnumpy2 = EMAN2_cppwrap.EMNumPy()
	bigbuffer = emnumpy2.register_numpy_to_emdata(buffer)

	#  end of CONES setup

	if( Blockdata["myid"] == Blockdata["main_node"] ):
		print( "  " )
		line = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime()) + " =>"
		print(  line, "Processing data for polar onx: %3d, nx: %3d, nxpolar: %3d, CTF: %s, preshift: %s,  MEM: %6.2fGB."%(Tracker["constants"]["nnxo"], Tracker["nxinit"], Tracker["nxpolar"], Tracker["constants"]["CTF"], preshift,orgsize/1.0e9) )

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
	refdirs = utilities.angles_to_normals(refang)


	#  We have to make sure the number of cones is the same on all nodes, this is due to the strange MPI problem
	#   that forced me to insert overall barrier into iterations over cones
	max_number_of_cones = number_of_cones
	max_number_of_cones = mpi.mpi_reduce(max_number_of_cones, 1, mpi.MPI_INT, mpi.MPI_MAX, 0, mpi.MPI_COMM_WORLD)
	if Blockdata["myid"] == Blockdata["main_node"] :  max_number_of_cones = int(max_number_of_cones[0])
	max_number_of_cones = mpi.mpi_bcast(max_number_of_cones, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD)
	max_number_of_cones = int(max_number_of_cones[0])



	#if( Blockdata["myid_on_node"] == 0):
	#	print("  FIRST  nima, nang, npsi, nshift, n_coarse_ang, nlocal_angles, n_coarse_psi, len(list_of_coarse_shifts), number_of_cones , max_number_of_cones ",nima, nang,npsi,nshifts,n_coarse_ang,nlocal_angles,n_coarse_psi, len(coarse_shifts),number_of_cones,max_number_of_cones)

	##firsti = True
	at = time.time()
	##eat = 0.0
	lima = 0  #  total counter of images
	#  PROCESSING OF CONES
	for icone in range(max_number_of_cones):
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		if( icone < number_of_cones ):  #  This is executed for individual number of cones, some nodes may have fewer.
			nang_start, nang_end = applications.MPI_start_end(len(assignments_of_refangles_to_cones[icone]), Blockdata["no_of_processes_per_group"], Blockdata["myid_on_node"])
			#if(Blockdata["color"] == 1):
			###print( " ZXZ11  ",Blockdata["color"],Blockdata["myid_on_node"],Blockdata["myid"],len(assignments_of_refangles_to_cones[icone]),nang_start, nang_end)

			for i in range(nang_start, nang_end, 1):  # This will take care of no of process on a node less than nang.  Some loops will not be executed
				ic = assignments_of_refangles_to_cones[icone][i]
				temp = projection.prgl(volprep,[ coarse_angles[ic][0],coarse_angles[ic][1],0.0, 0.0,0.0], 1, True)
				crefim = EMAN2_cppwrap.Util.Polar2Dm(temp, cnx, cnx, numr, mode)
				EMAN2_cppwrap.Util.Frngs(crefim, numr)
				EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)
				bigbuffer.insert_clip(crefim,(0,i) )

			mpi.mpi_barrier(Blockdata["shared_comm"])

			#if(Blockdata["myid"] == Blockdata["main_node"]):
			#	print( "  Reference projections generated for cone %4d: %10.1fmin"%(icone,(time()-at)/60.))
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

					#  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):  print("\n\n  INPUT PARAMS  ",im,phi,theta,psi,sx,sy)
					if preshift:
						sx = int(round(sx))
						sy = int(round(sy))
						dataimage  = fundamentals.cyclic_shift(original_data[im],sx,sy)
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
								dataimage = morphology.cosinemask(dataimage,radius = Tracker["constants"]["radius"])
							else:
								bckg = utilities.model_gauss_noise(1.0,Tracker["constants"]["nnxo"]+2,Tracker["constants"]["nnxo"])
								bckg.set_attr("is_complex",1)
								bckg.set_attr("is_fftpad",1)
								bckg = fundamentals.fft(filter.filt_table(bckg, oneover[dataimage.get_attr("particle_group")]))
								#  Normalize bckg noise in real space, only region actually used.
								st = EMAN2_cppwrap.Util.infomask(bckg, mask2D, False)
								bckg -= st[0]
								bckg /= st[1]
								dataimage = morphology.cosinemask(dataimage,radius = Tracker["constants"]["radius"], bckg = bckg)
					else:
						#  if no bckgnoise, do simple masking instead
						if apply_mask:  dataimage = morphology.cosinemask(dataimage,radius = Tracker["constants"]["radius"] )

					#  Apply varadj
					if not nonorm:
						EMAN2_cppwrap.Util.mul_scalar(dataimage, Tracker["avgvaradj"][procid]/wnorm)

					###  FT
					dataimage = fundamentals.fft(dataimage)
					sig = EMAN2_cppwrap.Util.rotavg_fourier( dataimage )
					accumulatepw[im] = sig[len(sig)//2:]+[0.0]

					#  We have to make sure the shifts are within correct range, shrinkage or not
					#set_params_proj(dataimage,[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
					
					if Blockdata["bckgnoise"]:
						temp = Blockdata["bckgnoise"][dataimage.get_attr("particle_group")]
						bckgn = EMAN2_cppwrap.Util.unroll1dpw(Tracker["nxinit"], [temp[i] for i in range(temp.get_xsize())])
					else:
						bckgn = EMAN2_cppwrap.Util.unroll1dpw(Tracker["nxinit"], [1.0]*600)
					bckgnoise = bckgn.copy()
					for j in range(Tracker["nxinit"]//2+1,Tracker["nxinit"]):  bckgn[0,j] = bckgn[0,Tracker["nxinit"]-j]

					if Tracker["constants"]["CTF"] :
						ctf_params = dataimage.get_attr("ctf")
						ctf_params.apix = ctf_params.apix/(float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"]))
						ctfa = morphology.ctf_img_real(Tracker["nxinit"], ctf_params)
						ctfs = ctfa
					##if( ( Blockdata["myid"] == Blockdata["main_node"])   and firsti ):
					##	dataimage.set_attr("is_complex",0)
					##	dataimage.write_image("dataimagefft.hdf")
					##	dataimage.set_attr("is_complex",1)
					dataml = fundamentals.fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False, False)
					data = []
					for iq in coarse_shifts:
						xx = iq[0]*shrink
						yy = iq[1]*shrink
						dss = fundamentals.fshift(dataml, xx, yy)
						dss.set_attr("is_complex",0)
						data.append(dss)

					#  This will get it to real space
					#dataimage = fpol(Util.mulnclreal(Util.mulnclreal(fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False), Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)
					#  Here we can reuse dataml to save time.  It was not normalized, but it does not matter as dataimage is for POLAR.  PAP 08/06/2018
					dataimage = fundamentals.fpol(EMAN2_cppwrap.Util.mulnclreal(EMAN2_cppwrap.Util.mulnclreal(dataml, EMAN2_cppwrap.Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)

					# Compute max number of angles on the fly
					lang = len(assignments_of_refangles_to_angles[im])
					###print("   BICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)

				if( ( Blockdata["myid"] == Blockdata["main_node"])  and  (lima%(max(1,nima/5)) == 0) and (lima>0)):
					print( "  Number of images :%7d   %5d  %5.1f"%(lima,nima,float(lima)/float(nima)*100.) + "%" +"   %10.1fmin"%((time.time()-at)/60.))
					##print( "  Number of images :%7d   %5d  %5.1f"%(lima,nima,float(lima)/float(nima)*100.) + "%" +"   %10.1fmin   %10.1fmin"%((time()-at)/60.,eat/60.0))
				lima += 1

				###print("  CONA1    ",Blockdata["myid"],lima)
				if( lima == 1 and procid == 0):
					###print("  CONA2    ",Blockdata["myid"])
					if( lenass > 0):
						###print("   CICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)
						keepfirst = lang*m_coarse_psi*n_coarse_shifts#150#lang*m_coarse_psi*len(coarse_shifts_shrank)  #500#min(200,nlocal_angles)#min(max(lxod1/25,200),lxod1)
						###if( Blockdata["myid"] == 18 and lima<5):  print(" START   nlocal_angles* m_coarse_psi*len(coarse_shifts_shrank",nlocal_angles,coarse_delta,len(coarse_shifts_shrank),keepfirst)
						#if( Blockdata["myid"] == Blockdata["main_node"] ):  
						lxod1 = EMAN2_cppwrap.Util.multiref_Crosrng_msg_stack_stepsi_local(dataimage, bigbuffer, \
								coarse_shifts_shrank,\
								assignments_of_refangles_to_angles[im], assignments_of_refangles_to_cones[icone],\
								numr, [coarse_angles[i][2] for i in range(n_coarse_ang)], \
								oldparams[im][2], c_coarse_psi, coarse_delta, cnx, keepfirst)

						##'''
						assert(len(lxod1)/3 == keepfirst)

						xod1 = np.ndarray((keepfirst),dtype='f4',order="C")
						#xod1.fill(1.0)
						xod2 = np.ndarray((keepfirst),dtype='int',order="C")
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
								temp = projection.prgl(volinit,[coarse_angles[iang][0],coarse_angles[iang][1],(coarse_angles[iang][2] + ipsi*coarse_delta)%360.0, 0.0,0.0], 1, False)
								##eat += time()-junk
								###   if( Blockdata["myid"] == Blockdata["main_node"] ):  print("  SELECTEDSTEPTWO  ",iln,refang[iang][0],refang[iang][1],refang[iang][2]+ipsi*Tracker["delta"],ishift,xod1[iln])
								temp.set_attr("is_complex",0)
								nrmref = numpy.sqrt(EMAN2_cppwrap.Util.innerproduct(temp, temp, None))
							##junk = time()
							#xod1[iln] = -Util.sqed(data[ishift], temp, ctfa, bckgnoise)
							xod1[iln] = EMAN2_cppwrap.Util.innerproduct(data[ishift], temp, None)/nrmref
							##eat += time()-junk
							##xod2[iln] = hashparams

						xod1 -= np.max(xod1)
						lina = np.argwhere(xod1 > Tracker["constants"]["expthreshold"])
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

						lina = np.argsort(xod1)
						xod1 = xod1[lina[::-1]]  # This sorts in reverse order
						xod2 = xod2[lina[::-1]]  # This sorts in reverse order
						np.exp(xod1, out=xod1)
						xod1 /= np.sum(xod1)
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
					keepf = utilities.wrap_mpi_gatherv(keepf, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
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
					keepf = utilities.wrap_mpi_bcast(keepf, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
					if(keepf == 0):
						global_def.ERROR("Too few images to estimate keepfirst","sxmeridien", 1, Blockdata["myid"])
						mpi.mpi_finalize()
						exit()
					###print("  STARTING8    ",Blockdata["myid"],keepf)
					Tracker["keepfirst"] = int(keepf)
					###if( Blockdata["myid"] == 0 ):  print("  keepfirst first ",Tracker["keepfirst"])

				else:
					if(lenass>0):
						###print("   DICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)
						###if( Blockdata["myid"] == 18 and lima<5):  print(" START   nlocal_angles* m_coarse_psi*len(coarse_shifts_shrank",nlocal_angles,coarse_delta,len(coarse_shifts_shrank),Tracker["keepfirst"])
						#if( Blockdata["myid"] == Blockdata["main_node"] ):
						keepfirst = lang*m_coarse_psi*n_coarse_shifts
						keepfirst = max(keepfirst*Tracker["keepfirst"]/100,min(keepfirst,3))
						lxod1 = EMAN2_cppwrap.Util.multiref_Crosrng_msg_stack_stepsi_local(dataimage, bigbuffer, \
								coarse_shifts_shrank,\
								assignments_of_refangles_to_angles[im], assignments_of_refangles_to_cones[icone],\
								numr, [coarse_angles[i][2] for i in range(n_coarse_ang)], \
								oldparams[im][2], c_coarse_psi, coarse_delta, cnx, keepfirst)
						#if( Blockdata["myid"] == Blockdata["main_node"] ):  print("  ")
						##'''
						assert( keepfirst == len(lxod1)/3 )
						xod1 = np.ndarray((keepfirst),dtype='f4',order="C")
						#xod1.fill(1.0)
						xod2 = np.ndarray((keepfirst),dtype='int',order="C")
						for iq in range(keepfirst):
							ioffset = 3*iq
							#          ishift         iang                      ipsi
							xod2[iq] = lxod1[ioffset] + lxod1[ioffset+1]*100000000 + lxod1[ioffset+2]*1000

						##'''

						#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

						#  Second step - find which coarse ones are significant

						# order by angular directions to save time on reprojections.
						ipsiandiang = xod2/1000
						lina = np.argsort(ipsiandiang)
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
								temp = projection.prgl(volinit,[coarse_angles[iang][0],coarse_angles[iang][1],(coarse_angles[iang][2] + ipsi*coarse_delta)%360.0, 0.0,0.0], 1, False)
								###   if( Blockdata["myid"] == Blockdata["main_node"] ):  print("  SELECTEDSTEPTWO  ",iln,refang[iang][0],refang[iang][1],refang[iang][2]+ipsi*Tracker["delta"],ishift,xod1[iln])
								##eat += time()-junk
								temp.set_attr("is_complex",0)
							##junk = time()
							peak = -EMAN2_cppwrap.Util.sqed(data[ishift], temp, ctfa, bckgnoise)
							##eat += time()-junk
							#  Note I replace ccc by eqdist
							xod1[iln] = peak
							##xod2[iln] = hashparams

						lina = np.argsort(xod1)
						xod1 = xod1[lina[::-1]]  # This sorts in reverse order
						xod2 = xod2[lina[::-1]]  # This sorts in reverse order

						xod1 -= xod1[0]

						#if( Blockdata["myid"] == Blockdata["main_node"]):
						#	#print("  PROJECT   ",im,lit,johi)#,cod2)
						#	for iln in range(len(xod1)):  print("  ROPE   ",iln,xod1[iln],xod2[iln])


						lina = np.argwhere(xod1 > Tracker["constants"]["expthreshold"])
						xod1 = xod1[lina]
						xod2 = xod2[lina]
						np.exp(xod1, out=xod1)
						xod1 /= np.sum(xod1)
						cumprob = 0.0
						lit = len(xod1)
						for j in range(len(xod1)):
							cumprob += xod1[j]
							if(cumprob > Tracker["constants"]["ccfpercentage"]):
								lit = j+1
								break

						#  To here
						###if( Blockdata["myid"] == 18 and lima<5):  print("  SECOND KEPT  ",lit)
						#if( lima<5):  print("  SECOND KEPT  ",lit)


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
						#	print(" FAILED  ",Blockdata["myid"],Tracker["keepfirst"],iln,lima,lit,hashparams,iang,xod2[:lit],xod1[:lit])
						#	mpi_finalize()
						#	exit()
						firstshifts[iln] = ishift
						#  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):
						#	ipsi = ipsiandiang%100000
						#	ishift = hashparams%1000
						#	###print("  SECONDPAR%04d  "%im,iln,hashparams, ipsi,iang,ishift)
						#	print(" SECONDPAR%04d  "%im,iln, refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, shifts[ishift])
					###del xod2
					###if( Blockdata["myid"] == 18 and lima<5):   print("  SECOND KEPT  ",im,lit)#,len(xod2),xod2,firstdirections,firstshifts)
					#mpi_barrier(MPI_COMM_WORLD)
					#mpi_finalize()
					#exit()

					#if(Blockdata["myid"] == Blockdata["main_node"]):  print("  FIFI ",firstdirections)
					#if(Blockdata["myid"] == Blockdata["main_node"]):  print("  GUGU ",firstshifts)
					# Find neighbors, ltabang contains positions on refang list, no psis
					###ltabang = nearest_many_full_k_projangles(refang, firstdirections, howmany = 5, sym=Tracker["constants"]["symmetry"])
					ltabang = find_nearest_k_refangles_to_many_angles(refdirs, firstdirections, Tracker["delta"], howmany = 4)
					###if( Blockdata["myid"] == Blockdata["main_node"]): print("  ltabang ",ltabang)
					##mpi_barrier(MPI_COMM_WORLD)
					##mpi_finalize()
					##exit()

					# ltabang has length lit, which is the same as length as firshifts.  However, original xod2 is still available,
					#   even though it is longer than lit.


					###ltabshi = [shiftneighbors[iln] for iln in firstshifts]
					#if(Blockdata["myid"] == Blockdata["main_node"]):  print("  HUHU ",ltabang)
					#if(Blockdata["myid"] == Blockdata["main_node"]):  print("  OUOU ",ltabshi)

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
						#if( Blockdata["myid"] == Blockdata["main_node"]): print(" tshifts  ",i1,len(tshifts))
						for i2 in range(4):
							iang = ltabang[i1][i2]
							for i3 in range(2):  # psi
								itpsi = int((coarse_angles[oldiang][2] + ipsi*coarse_delta - refang[iang][2]+360.0)/Tracker["delta"])
								itpsi = (itpsi + i3)%npsi
								for i4 in range(len(tshifts)):
									cod2.append(iang*100000000 + itpsi*1000 + tshifts[i4])
									#lol += 1
									#if( Blockdata["myid"] == Blockdata["main_node"] ):  print("  zibzi  ",i1,i2,i3,i4, lol)


					del xod1, xod2

					###if( Blockdata["myid"] == 18 and lima<5):   print("  THIRD   ",len(cod2))#,cod2)
					cod2 = list(set(cod2))
					cod1 = [[q/1000,i] for i,q in enumerate(cod2)]
					cod1.sort()

					lit = len(cod1)

					cod2 = np.asarray([cod2[cod1[i][1]] for i in range(lit)])

					cod1 = np.ndarray(lit,dtype='f4',order="C")
					#cod1.fill(np.finfo(dtype='f4').min)
					cod3 = np.ndarray(lit,dtype='f4',order="C")
					#cod3.fill(0.0)  #  varadj

					###if( Blockdata["myid"] == 18 and lima<5): print("  THIRD   ",im,lit)#,cod2)

					#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					#  Third ML part.  Now the actual chi2 distances, we compute reprojections on the fly.
					#  Make sure volprep has nxinit size
					data = [None]*nshifts
					johi = 0
					iln = 0
					prevdir = -1
					while(iln<lit):
						hashparams = cod2[iln]
						#if( Blockdata["myid"] == Blockdata["main_node"]): print("  COD2   ",im,lit,iln,cod2[iln])
						ipsiandiang	= hashparams/1000
						if(ipsiandiang != prevdir):
							prevdir = ipsiandiang
							ipsi = ipsiandiang%100000
							iang = ipsiandiang/100000
							##junk = time()
							temp = projection.prgl(volinit,[ refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, 0.0,0.0], 1, False)
							##eat += time()-junk
							temp.set_attr("is_complex",0)
							johi += 1
						while( ipsiandiang == cod2[iln]/1000 ):
							hashparams = cod2[iln]
							ishift = hashparams%1000
							if( data[ishift] == None ):
								xx = shifts[ishift][0]*shrink
								yy = shifts[ishift][1]*shrink
								data[ishift] = fundamentals.fshift(dataml, xx, yy)
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
						###if( Blockdata["myid"] == Blockdata["main_node"] and iln%1000 ==0):  print(" progress  ",iln,time()-at)
					#if( Blockdata["myid"] == Blockdata["main_node"]):
					#	print("  PROJECT   ",im,lit,johi)#,cod2)
					#	#for iln in range(lit):  print("  ROPE   ",iln,cod1[iln],cod2[iln],cod3[iln])
					del data
					del dataml

					lina = np.argsort(cod1)
					cod1 = cod1[lina[::-1]]  # This sorts in reverse order
					cod2 = cod2[lina[::-1]]  # This sorts in reverse order
					cod3 = cod3[lina[::-1]]  # This sorts in reverse order
					cod1 -= cod1[0]
					lina = np.argwhere(cod1 > Tracker["constants"]["expthreshold"])
					cod1 = cod1[lina]
					cod2 = cod2[lina]
					cod3 = cod3[lina]

					###if( Blockdata["myid"] == Blockdata["main_node"]):
					###for iui in range(len(lina)):
					###	for iui in range(len(cod1)):
					###		print("  MLML  ",iui,cod1[iui],exp(cod1[iui]),cod2[iui],cod3[iui])

					np.exp(cod1, out=cod1)
					cod1 /= np.sum(cod1)
					cumprob = 0.0
					for j in range(len(cod1)):
						cumprob += cod1[j]
						if(cumprob > Tracker["constants"]["ccfpercentage"]):
							lit = j+1
							break

					#  New norm is a sum of eq distances multiplied by their probabilities augmented by PW.
					norm_per_particle[im] = np.sum(cod1[:lit]*cod3[:lit]) + accumulatepw[im][reachpw]
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
						#	print(" NEWPAR%04d  "%im,iln,newpar[im][2][iln], refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, shifts[ishift])

					###if( Blockdata["myid"] == 18 and lima<5):   print("  FINALLY  ",im,lit)
					del cod1, cod2, cod3, lina
					###mpi_barrier(MPI_COMM_WORLD)
					###mpi_finalize()
					###exit()

	"""Multiline Comment24"""
	
	#  END OF CONES
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	if( Blockdata["myid"] == Blockdata["main_node"] ):
		print( "  Finished projection matching   %10.1fmin"%((time.time()-at)/60.))
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
	Tracker["avgvaradj"][procid] = utilities.bcast_number_to_all(Tracker["avgvaradj"][procid], Blockdata["main_node"])
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
	#	print( "  Projection matching finished : %10.1fmin"%((time()-at)/60.))
	return newpar, [1.0]*nima

def XYXali3D_local_polar_ccc(refang, shifts, coarse_angles, coarse_shifts, procid, original_data = None, oldparams = None, \
					preshift = False, apply_mask = True, nonorm = False, applyctf = True):
	global Tracker, Blockdata
	global Tracker, Blockdata
	pass#IMPORTIMPORTIMPORT from projection 	import prgs,prgl
	pass#IMPORTIMPORTIMPORT from fundamentals 	import fft
	pass#IMPORTIMPORTIMPORT from utilities 		import wrap_mpi_gatherv
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
	pass#IMPORTIMPORTIMPORT from utilities    import get_im, model_gauss_noise, set_params_proj, get_params_proj
	pass#IMPORTIMPORTIMPORT from fundamentals import fdecimate, fshift, fft
	pass#IMPORTIMPORTIMPORT from filter       import filt_ctf, filt_table
	pass#IMPORTIMPORTIMPORT from applications import MPI_start_end
	pass#IMPORTIMPORTIMPORT from math         import sqrt

	if(Blockdata["myid"] == Blockdata["main_node"]):
		print_dict(Tracker,"PROJECTION MATCHING parameters of buffered local polar")

	at = time.time()
	shrinkage = float(Tracker["nxpolar"])/float(Tracker["constants"]["nnxo"])
	shrink = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
	mode = "F"
	numr = Numrinit_local(1, radius, 1, mode)
	wr = alignment.ringwe(numr, mode)
	cnx = float(Tracker["nxpolar"]//2 + 1)

	##if(Blockdata["myid"] == Blockdata["main_node"]):  print("  RADIUS IN POLAR  ",radius,numr[-1])
	#  FINE SEARCH CONSTANTS
	nang = len(refang)
	ac_fine = numpy.cos(numpy.radians(Tracker["an"]))
	npsi = int(360./Tracker["delta"])
	mpsi = 2
	c_fine_psi = mpsi//2
	nshifts = len(shifts)
	n_fine_shifts = 4

	nima = len(original_data)
	mask = EMAN2_cppwrap.Util.unrollmask(Tracker["nxinit"])
	for j in range(Tracker["nxinit"]//2,Tracker["nxinit"]):  mask[0,j]=1.0
	mask2D = utilities.model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
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

	###if(Blockdata["myid"] == Blockdata["main_node"]): print("   TRETR  ",Tracker["constants"]["nnxo"],Tracker["nxinit"],reachpw,n_coarse_ang,coarse_delta,n_coarse_psi,m_coarse_psi,c_coarse_psi,n_coarse_shifts)

	#if(Blockdata["myid"] == Blockdata["main_node"]):
	#	print( original_data[0].get_attr("identifier") )
	#	print( original_data[1].get_attr("identifier") )

	disp_unit = np.dtype("f4").itemsize

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

		nxvol = utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nyvol = utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
		nzvol = utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

		if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = utilities.model_blank( nxvol,nyvol, nzvol)

		utilities.bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

		odo = projection.prep_vol( odo, npad = 2, interpolation_method = 1 )
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

	orgsizevol = utilities.bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
	nxvol = utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
	nyvol = utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
	nzvol = utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

	win_vol, base_vol  = mpi.mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
	sizevol = orgsizevol
	if( Blockdata["myid_on_node"] != 0 ):
		base_vol, = mpi.mpi_win_shared_query(win_vol, mpi.MPI_PROC_NULL)

	volbuf = np.frombuffer(np.core.multiarray.int_asbuffer(base_vol, sizevol*disp_unit), dtype = 'f4')
	volbuf = volbuf.reshape(nzvol, nyvol, nxvol)
	if( Blockdata["myid_on_node"] == 0 ):
		np.copyto(volbuf,ndo)
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

			nxvol = utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
			nyvol = utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])
			nzvol = utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"], mpi_comm = Blockdata["group_zero_comm"])

			if( Blockdata["myid"] != Blockdata["main_node"] ):  odo = utilities.model_blank( nxvol,nyvol, nzvol)

			utilities.bcast_EMData_to_all(odo, Blockdata["group_zero_myid"], source_node = Blockdata["main_node"], comm = Blockdata["group_zero_comm"])			

			odo = projection.prep_vol( odo, npad = 2, interpolation_method = 1 )
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

		orgsizevol = utilities.bcast_number_to_all(orgsizevol, source_node = Blockdata["main_node"])
		nxvol = utilities.bcast_number_to_all(nxvol, source_node = Blockdata["main_node"])
		nyvol = utilities.bcast_number_to_all(nyvol, source_node = Blockdata["main_node"])
		nzvol = utilities.bcast_number_to_all(nzvol, source_node = Blockdata["main_node"])

		win_volinit, base_volinit  = mpi.mpi_win_allocate_shared( sizevol*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
		sizevol = orgsizevol
		if( Blockdata["myid_on_node"] != 0 ):
			base_volinit, = mpi.mpi_win_shared_query(win_volinit, mpi.MPI_PROC_NULL)

		volbufinit = np.frombuffer(np.core.multiarray.int_asbuffer(base_volinit, sizevol*disp_unit), dtype = 'f4')
		volbufinit = volbufinit.reshape(nzvol, nyvol, nxvol)
		if( Blockdata["myid_on_node"] == 0 ):
			np.copyto(volbufinit,ndoinit)
			del odo,ndoinit

		#volinit = EMNumPy.assign_numpy_to_emdata(volbufinit)
		emnumpy4 = EMAN2_cppwrap.EMNumPy()
		volinit = emnumpy4.register_numpy_to_emdata(volbufinit)
		volinit.set_attr_dict({'is_complex':1,  'is_complex_x': 0, 'is_fftodd': 0, 'is_fftpad': 1, 'is_shuffled': 1,'npad': 2})
		if( Blockdata["myid_on_node"] == 0 ):  volinit.update()
		mpi.mpi_barrier(Blockdata["shared_comm"])
		###if( Blockdata["myid"] < 10 ):
		###	print(" VOLPREP  ",volinit[0],volinit[1],info(volinit))
	else:
		volinit = volprep
	#  End of replaced volprep

    ### local search is plugged in from here
	#  START CONES
	#  This has to be systematically done per node
	#
	crefim = EMAN2_cppwrap.Util.Polar2Dm(utilities.model_blank(Tracker["nxpolar"],Tracker["nxpolar"]), cnx, cnx, numr, mode)
	size_of_one_image = crefim.get_xsize()
	#  We will assume half of the memory is available.  We will do it betteer later.
	numberofrefs_inmem = int(Tracker["constants"]["memory_per_node"]/4/((size_of_one_image*disp_unit)/1.0e9))
	####if( Blockdata["myid_on_node"] == 0  ):  print( " MEMEST ", n_coarse_ang,numberofrefs_inmem)
	#  number of references that will fit into one mode
	normals_set = utilities.angles_to_normals(coarse_angles)
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
			assignments_of_refangles_to_cones[i] = utilities.wrap_mpi_gatherv(assignments_of_refangles_to_cones[i], 0, Blockdata["shared_comm"])
			doit = 1
			if( Blockdata["myid_on_node"] == 0 ):
				assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i]))
				#print(  " assignments_of_refangles_to_cones on zero ",i,len(assignments_of_refangles_to_cones[i]) )#,assignments_of_refangles_to_cones[i]

		for i,q in enumerate(assignments_of_refangles_to_cones):
			###if( Blockdata["myid_on_node"] == 0 ): print( " which refangles belong to which cone",Blockdata["color"],Blockdata["myid"],i,len(q) )#,q
			assignments_of_refangles_to_cones[i] = utilities.wrap_mpi_bcast(q, 0, Blockdata["shared_comm"])

	else:

		angledirs = utilities.angles_to_normals([u1[:3] for u1 in oldparams])

		number_of_cones = max(2, int(n_coarse_ang/numberofrefs_inmem*1.5 + 0.5))
		###if Blockdata["myid_on_node"] == 0:  print( " LENGTHS  ",Blockdata["color"],nima,n_coarse_ang, numberofrefs_inmem,number_of_cones)
		cont = True
		while  cont :
			#  Translate number of cones to cone_delta
			cone_delta, number_of_cones = number_of_cones_to_delta(number_of_cones)
			#  Generate cone_angles
			##if Blockdata["myid"] == 0:  print( "  WHILE  ",number_of_cones, cone_delta)
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

			#if Blockdata["myid"] == 0:  print(  "  number of cones ",number_of_cones,cone_delta, len(cone_angles))
			assert(number_of_cones == len(cone_angles) )

			conedirs = utilities.angles_to_normals(Blockdata["symclass"].symmetry_neighbors(cone_angles))
			neighbors = len(conedirs)/len(cone_angles)  #  Symmetry related neighbors
			#if Blockdata["myid"] == 0:  print(  "  neighbors  ",Blockdata["myid"],neighbors, cone_angles)
			#  assign data directions to cone_angles
			assignments_to_cones = utilities.assign_projdirs_f(angledirs, conedirs, neighbors)
			###print(  " assignments_to_cones ",Blockdata["myid"],len(assignments_to_cones),[len(q) for q in assignments_to_cones],assignments_to_cones[0])
			#  the above should have length of refdirs and each should have indexes of data that belong to this cone
			del conedirs
			#print "assignments_to_cones ",assignments_to_cones
			#  For each cone we have to find which refangles are needed to do the matching
			assignments_of_refangles_to_cones = [[] for i in range(len(assignments_to_cones))]
			assignments_of_refangles_to_angles = [[] for i in range(nima)]  # for each myid separately, these are angles on this myid


			for i,q in enumerate(assignments_to_cones):
				#  find assignments of refdirs to this cone within an of each angledirs, so we need symmetry neighbors set of refdirs.
				#if Blockdata["myid"] == 0:  print( "in loop ", Blockdata["myid"],i,len(q),q)

				if(len(q) == 0):
					# empty assignment, on a given CPU there are no images assigned to a given cone
					assignments_of_refangles_to_cones[i] = [-1]
				else:
					for m in q:
						#print " m ",m,len(angles)

						assignments_of_refangles_to_angles[m] = find_assignments_of_refangles_to_angles(normals_set, oldparams[m], Tracker["an"])
						#if Blockdata["myid"] == 0:  print( "assignments_of_refangles_to_angles[m] ", Blockdata["color"],i,m,assignments_of_refangles_to_angles[m])
						assignments_of_refangles_to_cones[i].extend(assignments_of_refangles_to_angles[m])

					#if( Blockdata["myid"] == 19 ):  print(  " doit0 ",Blockdata["myid"], i,assignments_of_refangles_to_cones[i],q,assignments_to_cones)
					assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i]))
				###if Blockdata["myid"] == 19:  print(  " assignments_of_refangles_to_cones on myid ",Blockdata["color"],Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]), assignments_of_refangles_to_cones[i] )
				assignments_of_refangles_to_cones[i] = utilities.wrap_mpi_gatherv(assignments_of_refangles_to_cones[i], 0, Blockdata["shared_comm"])
				###if Blockdata["myid"] == 0:  print(  " assignments_of_refangles_to_cones gatherv",Blockdata["color"],Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]) )
				doit = 1
				if( Blockdata["myid_on_node"] == 0 ):
					assignments_of_refangles_to_cones[i] = list(set(assignments_of_refangles_to_cones[i])-set([-1]))
					###if( Blockdata["myid_on_node"] == 0 ):  print(  " assignments_of_refangles_to_cones on zero ",i,len(assignments_of_refangles_to_cones[i]))
					#  POSSIBLE PROBLEM - IT IS POSSIBLE FOR A GIVEN CONE TO HAVE NO REFANGLES AND THUS BE EMPTY
					if( len(assignments_of_refangles_to_cones[i]) > numberofrefs_inmem ):
						number_of_cones = int(number_of_cones*1.25)
						#print(  " increased number_of_cones ",i,number_of_cones )
						doit = 0
				doit = utilities.bcast_number_to_all(doit, source_node = 0)
				number_of_cones = utilities.bcast_number_to_all(number_of_cones, source_node = 0, mpi_comm = Blockdata["shared_comm"] )
				###if( Blockdata["myid"] == 19 ):  print(  " doit ",Blockdata["myid"], i,doit ,assignments_of_refangles_to_cones[i],assignments_to_cones)
				if( doit == 0 ):  break

			if( doit == 1 ):
				cont = False


		for i,q in enumerate(assignments_of_refangles_to_cones):
			###if( Blockdata["myid_on_node"] == 0 ): print( " which refangles belong to which cone IOUT ",Blockdata["color"],Blockdata["myid"],i,len(q))#,q
			assignments_of_refangles_to_cones[i] = utilities.wrap_mpi_bcast(q, 0, Blockdata["shared_comm"])
			###if( Blockdata["myid_on_node"] == 0 ): print( " which refangles belong to which cone XOUT ",Blockdata["color"],Blockdata["myid"],i,len(q))#,q
			#if( myid == 1 ):
			#	print " which refangles belong to which cone",myid,i,len(assignments_of_refangles_to_cones[i])#,q

	if( Blockdata["myid"] == 0  ):  print( " number_of_cones: ",number_of_cones)
	
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
	###	print( " which refangles belong to which cone OUT ",Blockdata["color"],Blockdata["myid_on_node"],Blockdata["myid"],i,numberofrefs_inmem,nlocal_angles,len(q))#,q

	#  BIG BUFFER
	lenbigbuf = numberofrefs_inmem  #MAXIMUM NUMBER OF REFERENCE IN ONE OF THE CONES
	orgsize = lenbigbuf*size_of_one_image #  This is number of projections to be computed simultaneously times their size

	if( Blockdata["myid_on_node"] == 0 ): size = orgsize
	else:  size = 0

	win_sm, base_ptr  = mpi.mpi_win_allocate_shared( size*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
	size = orgsize
	if( Blockdata["myid_on_node"] != 0 ):
		base_ptr, = mpi.mpi_win_shared_query(win_sm, mpi.MPI_PROC_NULL)

	buffer = np.frombuffer(np.core.multiarray.int_asbuffer(base_ptr, size*disp_unit), dtype = 'f4')
	buffer = buffer.reshape(lenbigbuf, size_of_one_image)
	#bigbuffer = EMNumPy.assign_numpy_to_emdata(buffer)

	emnumpy2 = EMAN2_cppwrap.EMNumPy()
	bigbuffer = emnumpy2.register_numpy_to_emdata(buffer)

	#  end of CONES setup

	if( Blockdata["myid"] == Blockdata["main_node"] ):
		print( "  " )
		line = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime()) + " =>"
		print(  line, "Processing data for polar onx: %3d, nx: %3d, nxpolar: %3d, CTF: %s, preshift: %s,  MEM: %6.2fGB."%(Tracker["constants"]["nnxo"], Tracker["nxinit"], Tracker["nxpolar"], Tracker["constants"]["CTF"], preshift,orgsize/1.0e9) )

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
	refdirs = utilities.angles_to_normals(refang)


	#  We have to make sure the number of cones is the same on all nodes, this is due to the strange MPI problem
	#   that forced me to insert overall barrier into iterations over cones
	max_number_of_cones = number_of_cones
	max_number_of_cones = mpi.mpi_reduce(max_number_of_cones, 1, mpi.MPI_INT, mpi.MPI_MAX, 0, mpi.MPI_COMM_WORLD)
	if Blockdata["myid"] == Blockdata["main_node"] :  max_number_of_cones = int(max_number_of_cones[0])
	max_number_of_cones = mpi.mpi_bcast(max_number_of_cones, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD)
	max_number_of_cones = int(max_number_of_cones[0])



	#if( Blockdata["myid_on_node"] == 0):
	#	print("  FIRST  nima, nang, npsi, nshift, n_coarse_ang, nlocal_angles, n_coarse_psi, len(list_of_coarse_shifts), number_of_cones , max_number_of_cones ",nima, nang,npsi,nshifts,n_coarse_ang,nlocal_angles,n_coarse_psi, len(coarse_shifts),number_of_cones,max_number_of_cones)

	##firsti = True
	at = time.time()
	##eat = 0.0
	lima = 0  #  total counter of images
	#  PROCESSING OF CONES
	for icone in range(max_number_of_cones):
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		if( icone < number_of_cones ):  #  This is executed for individual number of cones, some nodes may have fewer.
			nang_start, nang_end = applications.MPI_start_end(len(assignments_of_refangles_to_cones[icone]), Blockdata["no_of_processes_per_group"], Blockdata["myid_on_node"])
			#if(Blockdata["color"] == 1):
			###print( " ZXZ11  ",Blockdata["color"],Blockdata["myid_on_node"],Blockdata["myid"],len(assignments_of_refangles_to_cones[icone]),nang_start, nang_end)

			for i in range(nang_start, nang_end, 1):  # This will take care of no of process on a node less than nang.  Some loops will not be executed
				ic = assignments_of_refangles_to_cones[icone][i]
				temp = projection.prgl(volprep,[ coarse_angles[ic][0],coarse_angles[ic][1],0.0, 0.0,0.0], 1, True)
				crefim = EMAN2_cppwrap.Util.Polar2Dm(temp, cnx, cnx, numr, mode)
				EMAN2_cppwrap.Util.Frngs(crefim, numr)
				EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)
				bigbuffer.insert_clip(crefim,(0,i) )

			mpi.mpi_barrier(Blockdata["shared_comm"])

			#if(Blockdata["myid"] == Blockdata["main_node"]):
			#	print( "  Reference projections generated for cone %4d: %10.1fmin"%(icone,(time()-at)/60.))
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

					#  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):  print("\n\n  INPUT PARAMS  ",im,phi,theta,psi,sx,sy)
					if preshift:
						sx = int(round(sx))
						sy = int(round(sy))
						dataimage  = fundamentals.cyclic_shift(original_data[im],sx,sy)
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
								dataimage = morphology.cosinemask(dataimage,radius = Tracker["constants"]["radius"])
							else:
								bckg = utilities.model_gauss_noise(1.0,Tracker["constants"]["nnxo"]+2,Tracker["constants"]["nnxo"])
								bckg.set_attr("is_complex",1)
								bckg.set_attr("is_fftpad",1)
								bckg = fundamentals.fft(filter.filt_table(bckg, oneover[dataimage.get_attr("particle_group")]))
								#  Normalize bckg noise in real space, only region actually used.
								st = EMAN2_cppwrap.Util.infomask(bckg, mask2D, False)
								bckg -= st[0]
								bckg /= st[1]
								dataimage = morphology.cosinemask(dataimage,radius = Tracker["constants"]["radius"], bckg = bckg)
					else:
						#  if no bckgnoise, do simple masking instead
						if apply_mask:  dataimage = morphology.cosinemask(dataimage,radius = Tracker["constants"]["radius"] )

					#  Apply varadj
					if not nonorm:
						EMAN2_cppwrap.Util.mul_scalar(dataimage, Tracker["avgvaradj"][procid]/wnorm)

					###  FT
					dataimage = fundamentals.fft(dataimage)
					sig = EMAN2_cppwrap.Util.rotavg_fourier( dataimage )
					accumulatepw[im] = sig[len(sig)//2:]+[0.0]

					#  We have to make sure the shifts are within correct range, shrinkage or not
					#set_params_proj(dataimage,[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
					if Blockdata["bckgnoise"]:
						temp = Blockdata["bckgnoise"][dataimage.get_attr("particle_group")]
						bckgn = EMAN2_cppwrap.Util.unroll1dpw(Tracker["nxinit"], [temp[i] for i in range(temp.get_xsize())])
					else:
						bckgn = EMAN2_cppwrap.Util.unroll1dpw(Tracker["nxinit"], [1.0]*600)
					bckgnoise = bckgn.copy()
					for j in range(Tracker["nxinit"]//2+1,Tracker["nxinit"]):  bckgn[0,j] = bckgn[0,Tracker["nxinit"]-j]

					if Tracker["constants"]["CTF"] :
						ctf_params = dataimage.get_attr("ctf")
						ctf_params.apix = ctf_params.apix/(float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"]))
						ctfa = morphology.ctf_img_real(Tracker["nxinit"], ctf_params)
						ctfs = ctfa
					##if( ( Blockdata["myid"] == Blockdata["main_node"])   and firsti ):
					##	dataimage.set_attr("is_complex",0)
					##	dataimage.write_image("dataimagefft.hdf")
					##	dataimage.set_attr("is_complex",1)
					dataml = fundamentals.fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False, False)
					data = []
					for iq in coarse_shifts:
						xx = iq[0]*shrink
						yy = iq[1]*shrink
						dss = fundamentals.fshift(dataml, xx, yy)
						dss.set_attr("is_complex",0)
						data.append(dss)

					#  This will get it to real space
					#dataimage = fpol(Util.mulnclreal(Util.mulnclreal(fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False), Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)
					#  Here we can reuse dataml to save time.  It was not normalized, but it does not matter as dataimage is for POLAR.  PAP 08/06/2018
					dataimage = fundamentals.fpol(EMAN2_cppwrap.Util.mulnclreal(EMAN2_cppwrap.Util.mulnclreal(dataml, EMAN2_cppwrap.Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)

					# Compute max number of angles on the fly
					lang = len(assignments_of_refangles_to_angles[im])
					###print("   BICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)

				if( ( Blockdata["myid"] == Blockdata["main_node"])  and  (lima%(max(1,nima/5)) == 0) and (lima>0)):
					print( "  Number of images :%7d   %5d  %5.1f"%(lima,nima,float(lima)/float(nima)*100.) + "%" +"   %10.1fmin"%((time.time()-at)/60.))
					##print( "  Number of images :%7d   %5d  %5.1f"%(lima,nima,float(lima)/float(nima)*100.) + "%" +"   %10.1fmin   %10.1fmin"%((time()-at)/60.,eat/60.0))
				lima += 1



				###print("  CONA1    ",Blockdata["myid"],lima)
				if( lima == 1 and procid == 0):
					###print("  CONA2    ",Blockdata["myid"])
					if( lenass > 0):
						###print("   CICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)
						keepfirst = lang*m_coarse_psi*n_coarse_shifts#150#lang*m_coarse_psi*len(coarse_shifts_shrank)  #500#min(200,nlocal_angles)#min(max(lxod1/25,200),lxod1)
						###if( Blockdata["myid"] == 18 and lima<5):  print(" START   nlocal_angles* m_coarse_psi*len(coarse_shifts_shrank",nlocal_angles,coarse_delta,len(coarse_shifts_shrank),keepfirst)
						#if( Blockdata["myid"] == Blockdata["main_node"] ):  
						lxod1 = EMAN2_cppwrap.Util.multiref_Crosrng_msg_stack_stepsi_local(dataimage, bigbuffer, \
								coarse_shifts_shrank,\
								assignments_of_refangles_to_angles[im], assignments_of_refangles_to_cones[icone],\
								numr, [coarse_angles[i][2] for i in range(n_coarse_ang)], \
								oldparams[im][2], c_coarse_psi, coarse_delta, cnx, keepfirst)

						##'''
						assert(len(lxod1)/3 == keepfirst)

						xod1 = np.ndarray((keepfirst),dtype='f4',order="C")
						#xod1.fill(1.0)
						xod2 = np.ndarray((keepfirst),dtype='int',order="C")
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
								temp = projection.prgl(volinit,[coarse_angles[iang][0],coarse_angles[iang][1],(coarse_angles[iang][2] + ipsi*coarse_delta)%360.0, 0.0,0.0], 1, False)
								##eat += time()-junk
								###   if( Blockdata["myid"] == Blockdata["main_node"] ):  print("  SELECTEDSTEPTWO  ",iln,refang[iang][0],refang[iang][1],refang[iang][2]+ipsi*Tracker["delta"],ishift,xod1[iln])
								temp.set_attr("is_complex",0)
								nrmref = numpy.sqrt(EMAN2_cppwrap.Util.innerproduct(temp, temp, None))
							##junk = time()
							#xod1[iln] = -Util.sqed(data[ishift], temp, ctfa, bckgnoise)
							xod1[iln] = EMAN2_cppwrap.Util.innerproduct(data[ishift], temp, None)/nrmref
							##eat += time()-junk
							##xod2[iln] = hashparams

						xod1 -= np.max(xod1)
						lina = np.argwhere(xod1 > Tracker["constants"]["expthreshold"])
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

						lina = np.argsort(xod1)
						xod1 = xod1[lina[::-1]]  # This sorts in reverse order
						xod2 = xod2[lina[::-1]]  # This sorts in reverse order
						np.exp(xod1, out=xod1)
						xod1 /= np.sum(xod1)
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
					keepf = utilities.wrap_mpi_gatherv(keepf, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
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
					keepf = utilities.wrap_mpi_bcast(keepf, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
					if(keepf == 0):
						global_def.ERROR("Too few images to estimate keepfirst","sxmeridien", 1, Blockdata["myid"])
						mpi.mpi_finalize()
						exit()
					###print("  STARTING8    ",Blockdata["myid"],keepf)
					Tracker["keepfirst"] = int(keepf)
					if( Blockdata["myid"] == 0 ):  print("  keepfirst first ",Tracker["keepfirst"])

				else:
					if(lenass>0):
						###print("   DICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)
						###if( Blockdata["myid"] == 18 and lima<5):  print(" START   nlocal_angles* m_coarse_psi*len(coarse_shifts_shrank",nlocal_angles,coarse_delta,len(coarse_shifts_shrank),Tracker["keepfirst"])
						#if( Blockdata["myid"] == Blockdata["main_node"] ):  
						lxod1 = EMAN2_cppwrap.Util.multiref_Crosrng_msg_stack_stepsi_local(dataimage, bigbuffer, \
								coarse_shifts_shrank,\
								assignments_of_refangles_to_angles[im], assignments_of_refangles_to_cones[icone],\
								numr, [coarse_angles[i][2] for i in range(n_coarse_ang)], \
								oldparams[im][2], c_coarse_psi, coarse_delta, cnx, Tracker["keepfirst"])
						#if( Blockdata["myid"] == Blockdata["main_node"] ):  print("  ")
						##'''

						xod1 = np.ndarray((Tracker["keepfirst"]),dtype='f4',order="C")
						#xod1.fill(1.0)
						xod2 = np.ndarray((Tracker["keepfirst"]),dtype='int',order="C")
						for iq in range(Tracker["keepfirst"]):
							ioffset = 3*iq
							#          ishift         iang                      ipsi
							xod2[iq] = lxod1[ioffset] + lxod1[ioffset+1]*100000000 + lxod1[ioffset+2]*1000

						##'''

						#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

						#  Second step - find which coarse ones are significant

						# order by angular directions to save time on reprojections.
						ipsiandiang = xod2/1000
						lina = np.argsort(ipsiandiang)
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
								temp = projection.prgl(volinit,[coarse_angles[iang][0],coarse_angles[iang][1],(coarse_angles[iang][2] + ipsi*coarse_delta)%360.0, 0.0,0.0], 1, False)
								###   if( Blockdata["myid"] == Blockdata["main_node"] ):  print("  SELECTEDSTEPTWO  ",iln,refang[iang][0],refang[iang][1],refang[iang][2]+ipsi*Tracker["delta"],ishift,xod1[iln])
								##eat += time()-junk
								temp.set_attr("is_complex",0)
							##junk = time()
							peak = -EMAN2_cppwrap.Util.sqed(data[ishift], temp, ctfa, bckgnoise)
							##eat += time()-junk
							#  Note I replace ccc by eqdist
							xod1[iln] = peak
							##xod2[iln] = hashparams

						lina = np.argsort(xod1)
						xod1 = xod1[lina[::-1]]  # This sorts in reverse order
						xod2 = xod2[lina[::-1]]  # This sorts in reverse order

						xod1 -= xod1[0]

						#if( Blockdata["myid"] == Blockdata["main_node"]):
						#	#print("  PROJECT   ",im,lit,johi)#,cod2)
						#	for iln in range(len(xod1)):  print("  ROPE   ",iln,xod1[iln],xod2[iln])


						lina = np.argwhere(xod1 > Tracker["constants"]["expthreshold"])
						xod1 = xod1[lina]
						xod2 = xod2[lina]
						np.exp(xod1, out=xod1)
						xod1 /= np.sum(xod1)
						cumprob = 0.0
						lit = len(xod1)
						for j in range(len(xod1)):
							cumprob += xod1[j]
							if(cumprob > Tracker["constants"]["ccfpercentage"]):
								lit = j+1
								break

						#  To here
						###if( Blockdata["myid"] == 18 and lima<5):  print("  SECOND KEPT  ",lit)
						#if( lima<5):  print("  SECOND KEPT  ",lit)


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
						#	print(" FAILED  ",Blockdata["myid"],Tracker["keepfirst"],iln,lima,lit,hashparams,iang,xod2[:lit],xod1[:lit])
						#	mpi_finalize()
						#	exit()
						firstshifts[iln] = ishift
						#  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):
						#	ipsi = ipsiandiang%100000
						#	ishift = hashparams%1000
						#	###print("  SECONDPAR%04d  "%im,iln,hashparams, ipsi,iang,ishift)
						#	print(" SECONDPAR%04d  "%im,iln, refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, shifts[ishift])
					###del xod2
					###if( Blockdata["myid"] == 18 and lima<5):   print("  SECOND KEPT  ",im,lit)#,len(xod2),xod2,firstdirections,firstshifts)
					#mpi_barrier(MPI_COMM_WORLD)
					#mpi_finalize()
					#exit()

					#if(Blockdata["myid"] == Blockdata["main_node"]):  print("  FIFI ",firstdirections)
					#if(Blockdata["myid"] == Blockdata["main_node"]):  print("  GUGU ",firstshifts)
					# Find neighbors, ltabang contains positions on refang list, no psis
					###ltabang = nearest_many_full_k_projangles(refang, firstdirections, howmany = 5, sym=Tracker["constants"]["symmetry"])
					ltabang = find_nearest_k_refangles_to_many_angles(refdirs, firstdirections, Tracker["delta"], howmany = 4)
					###if( Blockdata["myid"] == Blockdata["main_node"]): print("  ltabang ",ltabang)
					##mpi_barrier(MPI_COMM_WORLD)
					##mpi_finalize()
					##exit()

					# ltabang has length lit, which is the same as length as firshifts.  However, original xod2 is still available,
					#   even though it is longer than lit.


					###ltabshi = [shiftneighbors[iln] for iln in firstshifts]
					#if(Blockdata["myid"] == Blockdata["main_node"]):  print("  HUHU ",ltabang)
					#if(Blockdata["myid"] == Blockdata["main_node"]):  print("  OUOU ",ltabshi)

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
						#if( Blockdata["myid"] == Blockdata["main_node"]): print(" tshifts  ",i1,len(tshifts))
						for i2 in range(4):
							iang = ltabang[i1][i2]
							for i3 in range(2):  # psi
								itpsi = int((coarse_angles[oldiang][2] + ipsi*coarse_delta - refang[iang][2]+360.0)/Tracker["delta"])
								itpsi = (itpsi + i3)%npsi
								for i4 in range(len(tshifts)):
									cod2.append(iang*100000000 + itpsi*1000 + tshifts[i4])
									#lol += 1
									#if( Blockdata["myid"] == Blockdata["main_node"] ):  print("  zibzi  ",i1,i2,i3,i4, lol)


					del xod1, xod2

					###if( Blockdata["myid"] == 18 and lima<5):   print("  THIRD   ",len(cod2))#,cod2)
					cod2 = list(set(cod2))
					cod1 = [[q/1000,i] for i,q in enumerate(cod2)]
					cod1.sort()

					lit = len(cod1)

					cod2 = np.asarray([cod2[cod1[i][1]] for i in range(lit)])

					cod1 = np.ndarray(lit,dtype='f4',order="C")
					#cod1.fill(np.finfo(dtype='f4').min)
					cod3 = np.ndarray(lit,dtype='f4',order="C")
					#cod3.fill(0.0)  #  varadj

					###if( Blockdata["myid"] == 18 and lima<5): print("  THIRD   ",im,lit)#,cod2)

					#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					#  Third ML part.  Now the actual chi2 distances, we compute reprojections on the fly.
					#  Make sure volprep has nxinit size
					data = [None]*nshifts
					johi = 0
					iln = 0
					prevdir = -1
					while(iln<lit):
						hashparams = cod2[iln]
						#if( Blockdata["myid"] == Blockdata["main_node"]): print("  COD2   ",im,lit,iln,cod2[iln])
						ipsiandiang	= hashparams/1000
						if(ipsiandiang != prevdir):
							prevdir = ipsiandiang
							ipsi = ipsiandiang%100000
							iang = ipsiandiang/100000
							##junk = time()
							temp = projection.prgl(volinit,[ refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, 0.0,0.0], 1, False)
							##eat += time()-junk
							temp.set_attr("is_complex",0)
							johi += 1
						while( ipsiandiang == cod2[iln]/1000 ):
							hashparams = cod2[iln]
							ishift = hashparams%1000
							if( data[ishift] == None ):
								xx = shifts[ishift][0]*shrink
								yy = shifts[ishift][1]*shrink
								data[ishift] = fundamentals.fshift(dataml, xx, yy)
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
						###if( Blockdata["myid"] == Blockdata["main_node"] and iln%1000 ==0):  print(" progress  ",iln,time()-at)
					#if( Blockdata["myid"] == Blockdata["main_node"]):
					#	print("  PROJECT   ",im,lit,johi)#,cod2)
					#	#for iln in range(lit):  print("  ROPE   ",iln,cod1[iln],cod2[iln],cod3[iln])
					del data
					del dataml

					lina = np.argsort(cod1)
					cod1 = cod1[lina[::-1]]  # This sorts in reverse order
					cod2 = cod2[lina[::-1]]  # This sorts in reverse order
					cod3 = cod3[lina[::-1]]  # This sorts in reverse order
					cod1 -= cod1[0]
					lina = np.argwhere(cod1 > Tracker["constants"]["expthreshold"])
					cod1 = cod1[lina]
					cod2 = cod2[lina]
					cod3 = cod3[lina]

					###if( Blockdata["myid"] == Blockdata["main_node"]):
					###for iui in range(len(lina)):
					###	for iui in range(len(cod1)):
					###		print("  MLML  ",iui,cod1[iui],exp(cod1[iui]),cod2[iui],cod3[iui])

					np.exp(cod1, out=cod1)
					cod1 /= np.sum(cod1)
					cumprob = 0.0
					for j in range(len(cod1)):
						cumprob += cod1[j]
						if(cumprob > Tracker["constants"]["ccfpercentage"]):
							lit = j+1
							break

					#  New norm is a sum of eq distances multiplied by their probabilities augmented by PW.
					norm_per_particle[im] = np.sum(cod1[:lit]*cod3[:lit]) + accumulatepw[im][reachpw]
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
						#	print(" NEWPAR%04d  "%im,iln,newpar[im][2][iln], refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, shifts[ishift])

					###if( Blockdata["myid"] == 18 and lima<5):   print("  FINALLY  ",im,lit)
					del cod1, cod2, cod3, lina
					###mpi_barrier(MPI_COMM_WORLD)
					###mpi_finalize()
					###exit()

	"""Multiline Comment25"""
	
	
	
	#  END OF CONES
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	if( Blockdata["myid"] == Blockdata["main_node"] ):
		print( "  Finished projection matching   %10.1fmin"%((time.time()-at)/60.))
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
	Tracker["avgvaradj"][procid] = utilities.bcast_number_to_all(Tracker["avgvaradj"][procid], Blockdata["main_node"])
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
	#	print( "  Projection matching finished : %10.1fmin"%((time()-at)/60.))
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
		ccsum +=statistics.pearson(c1, c2)
	return ccsum/float(ny)
	
