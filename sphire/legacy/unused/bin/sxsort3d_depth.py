"""1
			if (do_freeze_groups == 0) and (nbox == 0) and (converged == 1) and \
			  (reset_number_of_groups == 0) and (Tracker["check_minimum_size"]):# Check it at the very last moment
				if (Blockdata["myid"] == Blockdata["main_node"]):
					Tracker["check_minimum_size"] = False
					stat_stable, Kmeans_size = AI_freeze_groups(minimum_grp_size, list_of_stable)
					log_main.add("--------- Checking minimum grp size ----------- ")
					for indx in range(len(list_of_stable)):
						log_main.add(" %3d   %8d   %8d"%(indx, Kmeans_size[indx], len(list_of_stable[indx])))
					log_main.add("  %12.1f    %8.3f   %12.1f    %12.1f "%(\
					   float(stat_stable[0]), sqrt(stat_stable[1]), float(stat_stable[2]), float(stat_stable[3])))
					new_stable = []
					# Remove small groups
					for indx in range(len(list_of_stable)):
						any = list_of_stable[indx]
						if Kmeans_size[indx]> sqrt(stat_stable[1])*.25:
							new_stable.append(any.tolist())
						else: unaccounted_list += any.tolist()
						
					if (len(list_of_stable) != len(new_stable)) and (len(new_stable)>=2):
						log_main.add("Minimum_grp_size checking removes %d groups "%(len(list_of_stable) -len(new_stable)))
						unaccounted_list  = sorted(unaccounted_list)
						converged = 0
						reset_number_of_groups = 1
						accounted_list, new_index = merge_classes_into_partition_list(new_stable)
						list_of_stable = new_stable
						write_text_row(new_index,         os.path.join(iter_dir, "Accounted.txt"))
						write_text_file(unaccounted_list, os.path.join(iter_dir, "Core_set.txt"))
					else:
						log_main.add("Minimum_grp_size checking removes no groups ") 
						converged = 1
				else:
					converged = 0
					Tracker  = 0
					list_of_stable = 0
				reset_number_of_groups = bcast_number_to_all(reset_number_of_groups, Blockdata["main_node"], MPI_COMM_WORLD)
				converged          = bcast_number_to_all(converged, Blockdata["main_node"], MPI_COMM_WORLD)
				Tracker            = wrap_mpi_bcast(Tracker,        Blockdata["main_node"], MPI_COMM_WORLD)
				list_of_stable     = wrap_mpi_bcast(list_of_stable, Blockdata["main_node"], MPI_COMM_WORLD)
			"""
"""2
		for iref in range(number_of_groups):
			if(Blockdata["myid"] == Blockdata["last_node"]):
				tag =7007
				refvol = get_im(os.path.join(Tracker["directory"],"vol_grp%03d_iter%03d.hdf"%(iref, total_iter)))
				nnn = refvol.get_xsize()
				if(Tracker["nxinit"] != nnn): refvol = fdecimate(refvol, Tracker["nxinit"], \
				    Tracker["nxinit"], Tracker["nxinit"], True, False)
				stat = Util.infomask(refvol, mask3D, False)
				refvol -= stat[0]
				if stat[1]!=0.0:Util.mul_scalar(refvol, 1.0/stat[1])
				refvol *=mask3D
				send_EMData(refvol, Blockdata["main_node"], tag, MPI_COMM_WORLD)
				del refvol
			if(Blockdata["myid"] == Blockdata["main_node"]):
				tag = 7007
				refvol = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			if(Blockdata["myid_on_node"] == 0):
				if(Blockdata["myid"] != Blockdata["main_node"]):
					refvol = model_blank(Tracker["nxinit"],Tracker["nxinit"], Tracker["nxinit"])
				bcast_EMData_to_all(refvol, Blockdata["group_zero_myid"], source_node = \
				     Blockdata["main_node"], comm = Blockdata["group_zero_comm"])	
				np.copyto(volbuf,EMNumPy.em2numpy(refvol))
				del refvol
			mpi_barrier(Blockdata["shared_comm"])
			ref_vol = emnumpy1.register_numpy_to_emdata(volbuf)
			if Tracker["constants"]["comparison_method"] =="cross": 
				ref_peaks = compare_two_images_cross(cdata, ref_vol, ctf_images)
			else: ref_peaks = compare_two_images_eucd(cdata, ref_vol, fdata, ctf_images)
			local_peaks[iref] = ref_peaks
			del ref_vol
			mpi_barrier(MPI_COMM_WORLD)
			"""
"""3
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
"""4
def assign_unaccounted_elements_mpi(glist_in, clusters_in, img_per_grp):
	# assign unaccounted images by group probabilities
	global Tracker, Blockdata
	pass#IMPORTIMPORTIMPORT import random
	pass#IMPORTIMPORTIMPORT import copy
	pass#IMPORTIMPORTIMPORT import numpy as np
	icut = int(1.5*img_per_grp)
	if Blockdata["myid"]== Blockdata["main_node"]:
		glist    = copy.copy(glist_in)
		clusters = copy.copy(clusters_in)
		for ic in range(len(clusters)):
			if len(clusters[ic])>(2*img_per_grp):
				shuffle(clusters[ic])
				glist.extend(clusters[ic][icut:])
				del clusters[ic][icut:]
		shuffle(glist)
		clusters = fill_clusters(clusters, glist)
	else: clusters = 0
	clusters = wrap_mpi_bcast(clusters, Blockdata["main_node"], MPI_COMM_WORLD)
	return clusters

def fill_clusters(clusters, ulist):
	pass#IMPORTIMPORTIMPORT from random import shuffle
	nuacc = len(ulist)
	shuffle(ulist)
	if nuacc>len(clusters):# deserves proceeding
		flist      = []
		clusters   = sorted(clusters, key=len, reverse = True)
		not_enough = True
		while not_enough and len(clusters)>1: 
			nl  = len(clusters[0])
			tot = 0
			for im in range(1, len(clusters)):
				tot += nl-len(clusters[im])
			if tot >nuacc:
				flist.append(clusters[0])
				del clusters[0]
			else: not_enough = False
		if len(clusters)==1:	
			flist.append(clusters[0].extend(ulist))
			new_clusters = []
			if len(flist)>=1:
				for a in flist:
					new_clusters.append(sorted(a))
			return new_clusters
		else:
			all = nuacc
			for im in range(len(clusters)):all +=len(clusters[im])
			group_size = all//len(clusters)
			nc  = 0
			for im in range(len(clusters)):
				nptl = group_size - len(clusters[im]) 
				clusters[im].extend(ulist[nc: nc+nptl])
				nc +=nptl
			if nc<nuacc:
				for im in range(nc, nuacc):
					clusters[im%len(clusters)].append(ulist[im])
			flist +=clusters
			new_clusters = []
			if len(flist)>=1:
				for a in flist:
					new_clusters.append(sorted(a))
			return new_clusters
	else:# only a few particles
		for im in range(nuacc):
			clusters[im/len(clusters)].append(ulist[im])
		new_clusters = []
		for a in clusters:
			new_clusters.append(sorted(a))
		return new_clusters
"""			
'''5
	vol_data = get_image_data(tvol)
	we_data =  get_image_data(tweight)
	#  tvol is overwritten, meaning it is also an output
	n_iter = 10

	if( Blockdata["myid_on_node"] == 0 ):  at = time()

	ifi = mpi_iterefa( vol_data.__array_interface__['data'][0] ,  we_data.__array_interface__['data'][0] , nx, ny, nz, maxr2, \
			Tracker["constants"]["nnxo"], Blockdata["myid_on_node"], color, Blockdata["no_of_processes_per_group"],  Blockdata["shared_comm"], n_iter)
	'''
'''6
	vol_data = get_image_data(tvol)
	we_data =  get_image_data(tweight)
	#  tvol is overwritten, meaning it is also an output
	n_iter = 10

	if( Blockdata["myid_on_node"] == 0 ):  at = time()
	ifi = mpi_iterefa( vol_data.__array_interface__['data'][0] ,  we_data.__array_interface__['data'][0] , nx, ny, nz, maxr2, \
			Tracker["constants"]["nnxo"], Blockdata["myid_on_node"], color, Blockdata["no_of_processes_per_group"],  Blockdata["shared_comm"], n_iter)	
	'''
"""7
		for ifreq in range(len(Tracker["constants"]["fsc_curve"])):
			Tracker["constants"]["fsc_curve"][ifreq] = \
		       Tracker["constants"]["fsc_curve"][ifreq]*2./(1.+Tracker["constants"]["fsc_curve"][ifreq])
		"""    
"""8
			for ifreq in range(len(Tracker["constants"]["fsc_curve"])):
				Tracker["constants"]["fsc_curve"][ifreq] = \
				  Tracker["constants"]["fsc_curve"][ifreq]*2./(1.+Tracker["constants"]["fsc_curve"][ifreq])
			"""
def copy_refinement_tracker(tracker_refinement):
	global Tracker, Blockdata
	for key, value in Tracker:
		try:
			value_refinement = tracker_refinement[key]
			if (value == None) and (value_refinement != None): 
				Tracker[key] = value_refinement
		except:
			if (Blockdata["myid"] == Blockdata["main_node"]): 
				print(key, " in sorting set as ", value, \
				     ", while in refinement, it is set as ", value_refinement)
	return

def fuzzy_to_hard(umat_arrays, cutoff= 0.8):
	pass#IMPORTIMPORTIMPORT import numpy as np
	(ndat, ngroups) = np.shape(umat_arrays)
	for i in range(ndat):
		if np.max(umat_arrays[i])>cutoff:
			for j in xrange(ngroups):
				if umat_arrays[i][j]>=cutoff:
					umat_arrays[i][j]   = 1.0
				else: umat_arrays[i][j] = 0.0
	return umat_arrays
	
def compute_umat_from_assignment(new_assign, old_assign, ngrps):
	pass#IMPORTIMPORTIMPORT import numpy as np
	pass#IMPORTIMPORTIMPORT from math import sqrt
	# 3:2
	umat = np.full((new_assign.shape[0],  ngrps), 0.0, dtype=np.float64)
	for im in range(new_assign.shape[0]):
		umat[im][new_assign[im]]  =  1.
		umat[im][old_assign[im]] +=  0.3
		norm = np.dot(umat[im], umat[im])
		umat[im] =np.multiply(umat[im], 1./numpy.sqrt(norm))
	return umat
####======================================================================================
def get_shrink_data_final(nxinit, procid, original_data = None, oldparams = None, \
		return_real = False, preshift = False, apply_mask = True, nonorm = False, npad = 1):
	global Tracker, Blockdata
	"""Multiline Comment2"""
	#from fundamentals import resample
	pass#IMPORTIMPORTIMPORT from utilities    import get_im, model_gauss_noise, set_params_proj, get_params_proj
	pass#IMPORTIMPORTIMPORT from fundamentals import fdecimate, fshift, fft
	pass#IMPORTIMPORTIMPORT from filter       import filt_ctf, filt_table
	pass#IMPORTIMPORTIMPORT from applications import MPI_start_end
	pass#IMPORTIMPORTIMPORT from math         import sqrt
	
	mask2D  	= utilities.model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	nima 		= len(original_data)
	shrinkage 	= nxinit/float(Tracker["constants"]["nnxo"])
	radius 	    = int(Tracker["constants"]["radius"]*shrinkage + 0.5)
	txm    	    = float(nxinit-(nxinit//2+1) - radius)
	txl    	    = float(radius - nxinit//2+1)

	if Blockdata["bckgnoise"] :
		oneover = []
		nnx     = Blockdata["bckgnoise"][0].get_xsize()
		for i in range(len(Blockdata["bckgnoise"])):
			temp = [0.0]*nnx
			for k in range(nnx):
				if( Blockdata["bckgnoise"][i].get_value_at(k) > 0.0):
					temp[k] = 1.0/numpy.sqrt(Blockdata["bckgnoise"][i].get_value_at(k))
			oneover.append(temp)
		del temp
	Blockdata["accumulatepw"][procid] = [None]*nima
	data = [None]*nima
	for im in range(nima):
		phi, theta, psi, sx, sy, wnorm = oldparams[im][0], oldparams[im][1], \
		    oldparams[im][2], oldparams[im][3], oldparams[im][4], oldparams[im][7]
		if preshift:
			sx = int(round(sx))
			sy = int(round(sy))
			data[im]  = fundamentals.cyclic_shift(original_data[im],sx,sy)
			#  Put rounded shifts on the list, note image has the original floats - check whether it may cause problems
			oldparams[im][3] = sx
			oldparams[im][4] = sy
			sx = 0.0
			sy = 0.0
		else:  data[im] = original_data[im].copy()
		st = EMAN2_cppwrap.Util.infomask(data[im], mask2D, False)
		data[im] -= st[0]
		data[im] /= st[1]
		if data[im].get_attr_default("bckgnoise", None): data[im].delete_attr("bckgnoise")
		#  Do bckgnoise if exists
		if Blockdata["bckgnoise"]:
			if apply_mask:
				if Tracker["constants"]["hardmask"]:
					data[im] = morphology.cosinemask(data[im],radius = Tracker["constants"]["radius"])
				else:
					bckg = utilities.model_gauss_noise(1.0,Tracker["constants"]["nnxo"]+2,Tracker["constants"]["nnxo"])
					bckg.set_attr("is_complex",1)
					bckg.set_attr("is_fftpad",1)
					bckg = fundamentals.fft(filter.filt_table(bckg, oneover[data[im].get_attr("particle_group")]))
					#  Normalize bckg noise in real space, only region actually used.
					st = EMAN2_cppwrap.Util.infomask(bckg, mask2D, False)
					bckg -= st[0]
					bckg /= st[1]
					data[im] = morphology.cosinemask(data[im],radius = Tracker["constants"]["radius"], bckg = bckg)
		else:
			#  if no bckgnoise, do simple masking instead
			if apply_mask:  data[im] = morphology.cosinemask(data[im],radius = Tracker["constants"]["radius"] )
		#  Apply varadj
		if not nonorm: EMAN2_cppwrap.Util.mul_scalar(data[im], Tracker["avgvaradj"][procid]/wnorm)
		#  FT
		data[im] = fundamentals.fft(data[im])
		sig = EMAN2_cppwrap.Util.rotavg_fourier( data[im] )
		Blockdata["accumulatepw"][procid][im] = sig[len(sig)//2:]+[0.0]
		if Tracker["constants"]["CTF"] :
			data[im]        = fundamentals.fdecimate(data[im], nxinit*npad, nxinit*npad, 1, False, False)
			ctf_params      = original_data[im].get_attr("ctf")
			ctf_params.apix = ctf_params.apix/shrinkage
			data[im].set_attr('ctf', ctf_params)
			data[im].set_attr('ctf_applied', 0)
			if return_real: data[im] = fundamentals.fft(data[im])
		else:
			ctf_params = original_data[im].get_attr_default("ctf", False)
			if ctf_params:
				ctf_params.apix = ctf_params.apix/shrinkage
				data[im].set_attr('ctf', ctf_params)
				data[im].set_attr('ctf_applied', 0)
			data[im] = fundamentals.fdecimate(data[im], nxinit*npad, nxinit*npad, 1, True, False)
			apix     = Tracker["constants"]["pixel_size"]
			data[im].set_attr('apix', apix/shrinkage)	
		#  We have to make sure the shifts are within correct range, shrinkage or not
		utilities.set_params_proj(data[im],[phi,theta,psi,max(min(sx*shrinkage,txm),txl),\
		    max(min(sy*shrinkage,txm),txl)])
		if not return_real: data[im].set_attr("padffted",1)
		data[im].set_attr("npad",npad)
		if Blockdata["bckgnoise"]:
			temp = Blockdata["bckgnoise"][data[im].get_attr("particle_group")]
			###  Do not adjust the values, we try to keep everything in the same Fourier values.
			data[im].set_attr("bckgnoise", [temp[i] for i in range(temp.get_xsize())])
	return data

###5
def find_smallest_group(clusters):
	min_size =[len(clusters[0]), [0]]
	for ic in range(1, len(clusters)):
		if len(cluster[ic]) < min_size[0]: min_size = [len(clusters[ic]), [ic]]
		elif len(cluster[ic]) == min_size[0]: min_size[1].append(ic)
	if len(min_size[1])>=1: random.shuffle(min_size[1])
	return min_size[1][0]

def even_assignment_alist_to_mclusters(glist, number_of_groups):
	# evenly assign glist to clusters
	pass#IMPORTIMPORTIMPORT import copy
	pass#IMPORTIMPORTIMPORT from random import shuffle
	if number_of_groups >0:
		clusters = [[] for i in range(number_of_groups)]
		ulist = copy.deepcopy(glist)
		nc = 0
		while len(ulist)>0:
			im =  nc%number_of_groups
			random.shuffle(ulist)
			clusters[im].append(ulist[0])
			del ulist[0]
			nc +=1
		return clusters
	else: return []
	
def MPI_volume_start_end(number_of_groups, ncolor, mycolor):
	igroup_start = int(round(float(number_of_groups)/ncolor*mycolor))
	igroup_end   = int(round(float(number_of_groups)/ncolor*(mycolor+1)))
	return igroup_start, igroup_end
####=================================================================	
###===================>>> refangles partition <<<====================	
def do3d(procid, data, newparams, refang, rshifts, norm_per_particle, myid, mpi_comm = -1):
	global Tracker, Blockdata
	#  Without filtration
	pass#IMPORTIMPORTIMPORT from reconstruction import recons3d_trl_struct_MPI
	if (mpi_comm < -1): mpi_comm = MPI_COMM_WORDLD
	if (Blockdata["subgroup_myid"]== Blockdata["main_node"]):
		if( procid == 0 ):
			if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
				os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
	shrinkage = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	tvol, tweight, trol = reconstruction.recons3d_trl_struct_MPI(myid = \
	    Blockdata["subgroup_myid"], main_node = Blockdata["main_node"], prjlist = data, \
		paramstructure = newparams, refang = refang, rshifts_shrank = [[q[0]*shrinkage,q[1]*shrinkage] for q in rshifts], \
		delta = Tracker["delta"], CTF = Tracker["constants"]["CTF"], upweighted = False, mpi_comm = mpi_comm, \
		target_size = (2*Tracker["nxinit"]+3), avgnorm = Tracker["avgvaradj"][procid], norm_per_particle = norm_per_particle)
	if Blockdata["subgroup_myid"] == Blockdata["main_node"]:
		tvol.set_attr("is_complex",0)
		tvol.write_image(   os.path.join(Tracker["directory"], "tempdir", "tvol_%01d_%03d.hdf"%(procid,Tracker["mainiteration"])))
		tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%01d_%03d.hdf"%(procid,Tracker["mainiteration"])))
		trol.write_image(   os.path.join(Tracker["directory"], "tempdir", "trol_%01d_%03d.hdf"%(procid,Tracker["mainiteration"])))
	mpi.mpi_barrier(mpi_comm)
	return
#######--------------------------------------------------------
def do3d_sorting_groups_rec3d(iteration, masterdir, log_main):
	global Tracker, Blockdata
	pass#IMPORTIMPORTIMPORT from utilities import get_im
	# reconstruct final unfiltered volumes from sorted clusters
	keepgoing = 1
	### ====
	Tracker["directory"]              = masterdir
	Tracker["constants"]["masterdir"] = masterdir
	Tracker["maxfrad"]                = Tracker["nxinit"]//2
	####
	if (Blockdata["no_of_groups"]>1): # multiple nodes
		sub_main_node_list = [-1 for i in range(Blockdata["no_of_groups"])]
		for index_of_colors in range(Blockdata["no_of_groups"]):
			for iproc in range(Blockdata["nproc"]-1):
				if (Blockdata["myid"]== iproc):
					if (Blockdata["color"] == index_of_colors) and (Blockdata["myid_on_node"] == 0):
						sub_main_node_list[index_of_colors] = Blockdata["myid"]
					utilities.wrap_mpi_send(sub_main_node_list, Blockdata["last_node"], mpi.MPI_COMM_WORLD)
				if (Blockdata["myid"] == Blockdata["last_node"]):
					dummy = utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
					for im in range(len(dummy)):
						if (dummy[im]>-1): sub_main_node_list[im] = dummy[im]
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		utilities.wrap_mpi_bcast(sub_main_node_list, Blockdata["last_node"], mpi.MPI_COMM_WORLD)
		####		
		if Tracker["number_of_groups"]%Blockdata["no_of_groups"]== 0: 
			nbig_loop = Tracker["number_of_groups"]//Blockdata["no_of_groups"]
		else: nbig_loop = Tracker["number_of_groups"]//Blockdata["no_of_groups"]+1
	
		big_loop_colors = [[] for i in range(nbig_loop)]
		big_loop_groups = [[] for i in range(nbig_loop)]
		nc = 0
		while nc <Tracker["number_of_groups"]:
			im =  nc//Blockdata["no_of_groups"]
			jm =  nc%Blockdata["no_of_groups"]
			big_loop_colors[im].append(jm)
			big_loop_groups[im].append(nc)
			nc +=1
		#####
		for iloop in range(nbig_loop):
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				Clusterdir = os.path.join(Tracker["directory"], "Cluster%d"%index_of_group, "main%03d"%iteration)
				if(Blockdata["myid"] == Blockdata["last_node"]):
					tvol2 		= utilities.get_im(os.path.join(Clusterdir, "tempdir", "tvol_0_%03d.hdf"%iteration))
					tweight2 	= utilities.get_im(os.path.join(Clusterdir, "tempdir", "tweight_0_%03d.hdf"%iteration))
					treg2 		= utilities.get_im(os.path.join(Clusterdir, "tempdir", "trol_0_%03d.hdf"%iteration))
					tag      = 7007
					utilities.send_EMData(tvol2,    sub_main_node_list[index_of_colors],  tag, mpi.MPI_COMM_WORLD)
					utilities.send_EMData(tweight2, sub_main_node_list[index_of_colors],  tag, mpi.MPI_COMM_WORLD)
					utilities.send_EMData(treg2,    sub_main_node_list[index_of_colors],  tag, mpi.MPI_COMM_WORLD)
				elif (Blockdata["myid"] == sub_main_node_list[index_of_colors]):
					tag      = 7007
					tvol2    = utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
					tweight2 = utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
					treg2    = utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				if Blockdata["color"] == index_of_colors:
					if( Blockdata["myid_on_node"] != 0):
						tvol2 		= utilities.model_blank(1)
						tweight2 	= utilities.model_blank(1)
						treg2		= utilities.model_blank(1)
					tvol2 = steptwo_mpi(tvol2, tweight2, treg2, None, False, color = index_of_colors) # has to be False!!!
					del tweight2, treg2
				mpi.mpi_barrier(Blockdata["shared_comm"])
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				if (Blockdata["color"] == index_of_colors) and (Blockdata["myid_on_node"] == 0):
					tag = 7007
					utilities.send_EMData(tvol2, Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
				elif(Blockdata["myid"] == Blockdata["last_node"]):
					tag = 7007
					tvol2 = utilities.recv_EMData(sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
					tvol2.write_image(os.path.join(Tracker["directory"], "vol_unfiltered_0_grp%03d.hdf"%index_of_group))
					del tvol2
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				Clusterdir = os.path.join(Tracker["directory"], "Cluster%d"%index_of_group, "main%03d"%iteration)
				if(Blockdata["myid"] == Blockdata["last_node"]):
					tvol2 		= utilities.get_im(os.path.join(Clusterdir, "tempdir", "tvol_1_%03d.hdf"%iteration))
					tweight2 	= utilities.get_im(os.path.join(Clusterdir, "tempdir", "tweight_1_%03d.hdf"%iteration))
					treg2 		= utilities.get_im(os.path.join(Clusterdir, "tempdir", "trol_1_%03d.hdf"%iteration))
					tag      = 7007
					utilities.send_EMData(tvol2,    sub_main_node_list[index_of_colors], tag,mpi.MPI_COMM_WORLD)
					utilities.send_EMData(tweight2, sub_main_node_list[index_of_colors], tag,mpi.MPI_COMM_WORLD)
					utilities.send_EMData(treg2,    sub_main_node_list[index_of_colors], tag,mpi.MPI_COMM_WORLD)
				
				elif (Blockdata["myid"] == sub_main_node_list[index_of_colors]):
					tag      = 7007
					tvol2       = utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
					tweight2    = utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
					treg2       = utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		
			for im in range(len(big_loop_colors[iloop])):
				index_of_group     = big_loop_groups[iloop][im]
				index_of_colors    = big_loop_colors[iloop][im]
				Tracker["maxfrad"] = Tracker["nxinit"]//2
				if (Blockdata["color"] == index_of_colors):
					if( Blockdata["myid_on_node"] != 0):
						tvol2 		= utilities.model_blank(1)
						tweight2 	= utilities.model_blank(1)
						treg2		= utilities.model_blank(1)
					tvol2 = steptwo_mpi(tvol2, tweight2, treg2, None, False, color = index_of_colors) # has to be False!!!
					del tweight2, treg2
				mpi.mpi_barrier(Blockdata["shared_comm"])
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				if (Blockdata["color"] == index_of_colors) and (Blockdata["myid_on_node"] == 0):
					tag = 7007
					utilities.send_EMData(tvol2, Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
				elif(Blockdata["myid"] == Blockdata["last_node"]):
					tag = 7007
					tvol2 = utilities.recv_EMData(sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
					tvol2.write_image(os.path.join(Tracker["directory"], "vol_unfiltered_1_grp%03d.hdf"%index_of_group))
					del tvol2
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			
	else: # Single workstation
		Tracker["maxfrad"] = Tracker["nxinit"]//2
		for index_of_group in range(Tracker["number_of_groups"]):
			Clusterdir = os.path.join(Tracker["directory"], "Cluster%d"%index_of_group, "main%03d"%iteration)
			
			if(Blockdata["myid"] == Blockdata["last_node"]):
				tvol2 		= utilities.get_im(os.path.join(Clusterdir, "tempdir", "tvol_0_%03d.hdf"%iteration))
				tweight2 	= utilities.get_im(os.path.join(Clusterdir, "tempdir", "tweight_0_%03d.hdf"%iteration))
				treg2 		= utilities.get_im(os.path.join(Clusterdir, "tempdir", "trol_0_%03d.hdf"%iteration))
				tag      = 7007
				utilities.send_EMData(tvol2, Blockdata["main_node"],    tag, mpi.MPI_COMM_WORLD)
				utilities.send_EMData(tweight2, Blockdata["main_node"], tag, mpi.MPI_COMM_WORLD)
				utilities.send_EMData(treg2, Blockdata["main_node"],    tag, mpi.MPI_COMM_WORLD)
			elif (Blockdata["myid"] == Blockdata["main_node"]):
				tag      = 7007
				tvol2    = utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
				tweight2 = utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
				treg2    = utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			if Blockdata["myid"] != Blockdata["main_node"]:
				tvol2 		= utilities.model_blank(1)
				tweight2 	= utilities.model_blank(1)
				treg2		= utilities.model_blank(1)
			tvol2 = steptwo_mpi(tvol2, tweight2, treg2, None, False, color = 0) # has to be False!!!
			del tweight2, treg2
			if( Blockdata["myid"] == Blockdata["main_node"]):
				tvol2.write_image(os.path.join(Tracker["directory"], "vol_unfiltered_0_grp%03d.hdf"%index_of_group))
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			
			if(Blockdata["myid"] == Blockdata["last_node"]):
				tvol2 		= utilities.get_im(os.path.join(Clusterdir, "tempdir", "tvol_1_%03d.hdf"%iteration))
				tweight2 	= utilities.get_im(os.path.join(Clusterdir, "tempdir", "tweight_1_%03d.hdf"%iteration))
				treg2 		= utilities.get_im(os.path.join(Clusterdir, "tempdir", "trol_1_%03d.hdf"%iteration))
				tag      = 7007
				utilities.send_EMData(tvol2,    Blockdata["main_node"],  tag, mpi.MPI_COMM_WORLD)
				utilities.send_EMData(tweight2, Blockdata["main_node"],  tag, mpi.MPI_COMM_WORLD)
				utilities.send_EMData(treg2,    Blockdata["main_node"],  tag, mpi.MPI_COMM_WORLD)
			elif (Blockdata["myid"] == Blockdata["main_node"]):
				tag      = 7007
				tvol2    = utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
				tweight2 = utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
				treg2    = utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			if Blockdata["myid"] != Blockdata["main_node"]:
				tvol2 		= utilities.model_blank(1)
				tweight2 	= utilities.model_blank(1)
				treg2		= utilities.model_blank(1)
			tvol2 = steptwo_mpi(tvol2, tweight2, treg2, None, False, color = 0) # has to be False!!!
			del tweight2, treg2
			if( Blockdata["myid"] == Blockdata["main_node"]):
				tvol2.write_image(os.path.join(Tracker["directory"], "vol_unfiltered_1_grp%03d.hdf"%index_of_group))
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			
	keepgoing = utilities.bcast_number_to_all(keepgoing, source_node = Blockdata["main_node"], mpi_comm = mpi.MPI_COMM_WORLD) # always check 
	Tracker = utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"])
	if keepgoing == 0: 
		global_def.ERROR("do3d_sorting_groups_trl_iter  %s"%os.path.join(Tracker["directory"], \
		   "tempdir"),"do3d_sorting_groups_trl_iter", 1, Blockdata["myid"]) 
	return
####=====<-----------------------------------------------------------------------------
### nofsc rec3d
def stacksize(since=0.0):
    '''Return stack size in bytes.
    '''
    return _VmB('VmStk:') - since
###=======================================================================================
####=============>>> Parsers, settings, and IO related utilities <<<======================
def AI_freeze_groups(min_size, groups):
	number_of_groups = len(groups)
	Kmeans_size      = [ None for im in range(number_of_groups)]
	if number_of_groups>1:
		for im in range(number_of_groups):	
			Kmeans_size[im] = max((len(groups[im]) - sum(min_size)//2), 1)
		return statistics.table_stat(Kmeans_size), Kmeans_size
	else:
		return [Kmeans_size[0], 0.0, Kmeans_size[0], Kmeans_size[0]], Kmeans_size
###---------------------------------------------------------------------------------------
