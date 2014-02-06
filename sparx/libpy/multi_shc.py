# For some reason these programs get stuck on MPI if I change the order of programs in the file.  Strange, PAP.
def generate_uneven_projections_directions(count, half_sphere=False, output_filename=None, 
										   density_ratio_on_poles_and_equator=1.0): 
	"""
	Generates set of uneven distributed projections directions with given size 
	Returns: list of projections directions
	"""
	from random import random, uniform
	from utilities import write_text_row
	from math import sin, pi
	max_theta = 180.0
	if half_sphere: max_theta = 90.0
	a1 = []
	while len(a1) < count:
		theta = uniform(0.0, max_theta)
		if random() < abs(sin(theta/180.0*pi)):
			if theta > 90.0:
				angle_from_pole = 180.0 - theta
			else:
				angle_from_pole = theta
			density_ratio = density_ratio_on_poles_and_equator + (angle_from_pole / 90.0) * (1.0 - density_ratio_on_poles_and_equator)
			if density_ratio_on_poles_and_equator > 1.0:
				density_ratio /= density_ratio_on_poles_and_equator
			if random() < density_ratio:
				a1.append( [uniform(0.0, 360.0), theta] )
	for pos in a1:
		for i in range(2,5):
			if len(pos) <= i:
				pos.append(0.0)
			else:
				pos[i] = 0.0
	if output_filename != None and output_filename != "":
		write_text_row(a1, output_filename)
	return a1


def mult_transform(v1, v2):
	from EMAN2 import Transform
	T1 = Transform({"type":"spider","phi":v1[0],"theta":v1[1],"psi":v1[2],"tx":v1[3],"ty":v1[4],"tz":0.0,"mirror":0,"scale":1.0})
	T2 = Transform({"type":"spider","phi":v2[0],"theta":v2[1],"psi":v2[2],"tx":v2[3],"ty":v2[4],"tz":0.0,"mirror":0,"scale":1.0})
	T = T1*T2
	return [ T.get_params("spider")["phi"], T.get_params("spider")["theta"], T.get_params("spider")["psi"], T.get_params("spider")["tx"], T.get_params("spider")["ty"]  ]


def orient_params(params, indexes=None):
	from utilities import rotation_between_anglesets
	from pixel_error import angle_diff

	m = len(params)
	n = len(params[0])

	if indexes == None:
		indexes = range(n)

	for i in xrange(1,m):
		cmp_par_i = []
		cmp_par_0 = []
		for j in indexes:
			cmp_par_i.append(params[i][j])
			cmp_par_0.append(params[0][j])
		t1,t2,t3 = rotation_between_anglesets(cmp_par_i, cmp_par_0)
		rot = [t1, t2, t3, 0.0, 0.0]
		for j in xrange(n):
			params[i][j] = mult_transform(params[i][j], rot)
		# mirror checking
		psi_diff = angle_diff( [params[i][j][2] for j in indexes], [params[0][j][2] for j in indexes] )
		if(abs(psi_diff-180.0) <90.0):
			#mirror
			for j in xrange(n):
				params[i][j][2] = (params[i][j][2] + 180.0) % 360.0


def shuffle_configurations(params):
	from random import shuffle

	m = len(params)
	n = len(params[0])
	new_params = [ ([0]*n) for i in xrange(m) ]

	for i in xrange(n):
		src = range(m)
		shuffle(src)
		for j in xrange(m):
			new_params[j][i] = params[src[j]][i]

	return new_params


def calculate_matrix_rot(projs):
	from utilities import rotation_between_anglesets
	sc = len(projs)
	matrix_rot  = [[[0.0,0.0,0.0,0.0,0.0] for i in xrange(sc)] for k in xrange(sc)]
	for i in xrange(sc):
		for j in xrange(i):
			t1, t2, t3 = rotation_between_anglesets(projs[i], projs[j])
			matrix_rot[i][j] = [t1, t2, t3, 0.0, 0.0]
	return matrix_rot


# returns subset_for_threshold, subset_for_minimal_size[, threshold_for_thr_subset, threshold_for_min_subset]
def find_common_subset_3(projs, target_threshold, minimal_subset_size=3, sym = "c1", thresholds=False):
	from global_def import Util

	n = len(projs[0])
	sc = len(projs)

	subset = range(n)

	minimal_subset_size = min( minimal_subset_size, n)

	res_thr_subset  = None
	res_size_subset = None
	error_thr_subset = -1
	error_size_subset = -1

	for iIter in xrange(n-2):
		projs2 = [0.0]*sc
		trans_projs = []
		for iConf in xrange(sc):
			projs2[iConf] = []
			for i in subset:
				projs2[iConf].append(projs[iConf][i][:])
				trans_projs.extend(projs[iConf][i][0:5])
		if( sym[0] == "d"):		matrix_rot  = [[[0.0,0.0,0.0,0.0,0.0] for i in xrange(sc)] for k in xrange(sc)]
		else:					matrix_rot = calculate_matrix_rot(projs2)

		trans_matrix = []
		for i in xrange(sc):
			for j in xrange(i):
				trans_matrix.extend(matrix_rot[i][j][0:3])
		avg_diff_per_image = Util.diff_between_matrix_of_3D_parameters_angles(trans_projs, trans_matrix)
		#print avg_diff_per_image
		max_error = -1.0
		the_worst_proj = -1
		for i in xrange(len(avg_diff_per_image)):
			if avg_diff_per_image[i] > max_error:
				max_error = avg_diff_per_image[i]
				the_worst_proj = subset[i]
		if max_error <= target_threshold:
			res_thr_subset = subset[:]
			error_thr_subset = max_error
			break
		if len(subset) == minimal_subset_size:
			res_size_subset = subset[:]
			error_size_subset = max_error
		subset.remove(the_worst_proj)

	if res_thr_subset == None:
		res_thr_subset = subset
		error_thr_subset = max_error

	if res_size_subset == None:
		res_size_subset = res_thr_subset
		error_size_subset = error_thr_subset

	if thresholds:
		return res_thr_subset, res_size_subset, error_thr_subset, error_size_subset
	return res_thr_subset, res_size_subset, avg_diff_per_image


# parameters: list of (all) projections | reference volume | ...
#  Genetic programming version
#  The data structure:
#  [[L2, [parameters row-wise]], [], []number_of_runs ]
#  It is kept on main proc
def ali3d_multishc(stack, ref_vol, ali3d_options, mpi_comm = None, log = None, number_of_runs=2 ):

	from alignment    import Numrinit, prepare_refrings, proj_ali_incore_local, shc
	from utilities    import model_circle, get_input_from_string, get_params_proj, set_params_proj, wrap_mpi_gatherv, wrap_mpi_bcast, wrap_mpi_send, wrap_mpi_recv
	from mpi          import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier, mpi_comm_split, mpi_comm_free
	from projection   import prep_vol
	from statistics   import hist_list
	from applications import MPI_start_end
	from filter       import filt_ctf
	from global_def   import Util, ERROR
	from time         import time
	from random       import shuffle

	ir     = ali3d_options.ir
	rs     = ali3d_options.rs
	ou     = ali3d_options.ou
	xr     = ali3d_options.xr
	yr     = ali3d_options.yr
	ts     = ali3d_options.ts
	an     = ali3d_options.an
	sym    = ali3d_options.sym
	sym = sym[0].lower() + sym[1:]
	delta  = ali3d_options.delta
	center = ali3d_options.center
	maxit  = ali3d_options.maxit
	CTF    = ali3d_options.CTF
	ref_a  = ali3d_options.ref_a

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	if log == None:
		from logger import Logger
		log = Logger()

	number_of_proc = mpi_comm_size(mpi_comm)
	myid           = mpi_comm_rank(mpi_comm)
	main_node = 0

	if myid == main_node:
		log.add("Start ali3d_multishc")

	if number_of_proc < number_of_runs:
		ERROR("number_of_proc < number_of_runs","ali3d_multishc")
	
	mpi_subcomm = mpi_comm_split(mpi_comm, myid % number_of_runs, myid / number_of_runs)
	mpi_subrank = mpi_comm_rank(mpi_subcomm)
	mpi_subsize = mpi_comm_size(mpi_subcomm)
	mpi_subroots = range(number_of_runs)

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

	vol = ref_vol
	nx      = vol.get_xsize()
	if last_ring < 0:
		last_ring = int(nx/2) - 2

	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	mask2D  = model_circle(last_ring,nx,nx) - model_circle(first_ring,nx,nx)

	if myid == main_node:
		list_of_particles = range(len(stack))
		total_nima = len(list_of_particles)
	else:
		list_of_particles = None
		total_nima = None
	total_nima = wrap_mpi_bcast(total_nima, main_node, mpi_comm)
	list_of_particles = wrap_mpi_bcast(list_of_particles, main_node, mpi_comm)
	nima = len(list_of_particles)

	image_start, image_end = MPI_start_end(total_nima, mpi_subsize, mpi_subrank)

	data = [ stack[im] for im in list_of_particles ]
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		ctf_applied = data[im].get_attr_default('ctf_applied', 0)
		if CTF and ctf_applied == 0:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)

	cs = [0.0]*3
	total_iter = 0

	# initialize GA data structure [ [L2, [params]], ...]
	if myid == main_node:
		GA = [ [0.0, [[0.0,0.,0.0,0.0,0.0] for j in xrange(total_nima)]] for i in xrange(number_of_runs)]

	orient_and_shuffle = False

	# do the projection matching
	for N_step in xrange(lstp):  # At this point there is just one value here, it cannot loop.
		
		terminate = 0
		Iter = 0
		while Iter < max_iter and terminate == 0:

			Iter += 1
			total_iter += 1

			mpi_barrier(mpi_comm)
			if myid == main_node:
				log.add("ITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step]))
				start_time = time()

			#=========================================================================
			# build references
			volft, kb = prep_vol(vol)
			refrings = prepare_refrings(volft, kb, nx, delta[N_step], ref_a, sym, numr, MPI=mpi_subcomm)
			del volft, kb
			all_ref_dirs = []
			for r in refrings:
				all_ref_dirs.append( [r.get_attr("phi"), r.get_attr("theta")] )
			#=========================================================================

			mpi_barrier(mpi_comm)
			if myid == main_node:
				log.add("Time to prepare rings: %f\n" % (time()-start_time))
				start_time = time()

			#=========================================================================
			if total_iter == 1 or orient_and_shuffle:
				# adjust params to references, calculate psi+shifts, calculate previousmax
				for im in xrange(nima):
					stable = data[im].get_attr_default("stable", 0)
					if stable == 0:
						data[im].set_attr("previousmax", -1.0e23)
						data[im].set_attr("stable", 1)
					else:
						peak, temp = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step], delta[N_step]*0.7 )
						data[im].set_attr("previousmax", peak)
				if myid == main_node:
					log.add("Time to calculate first psi+shifts+previousmax: %f\n" % (time()-start_time))
					start_time = time()
			#=========================================================================

			mpi_barrier(mpi_comm)
			if myid == main_node: start_time = time()
			#=========================================================================
			# alignment
			if mpi_subrank == 0:
				pixer = []
				number_of_checked_refs = 0
				proj_ids_to_process = range(nima)
				shuffle(proj_ids_to_process)
			while True:
				# --------- broadcast projections ids
				if mpi_subrank == 0:
					if len(proj_ids_to_process) >= mpi_subsize:
						proj_ids = proj_ids_to_process[:(mpi_subsize)]
					else:
						proj_ids = proj_ids_to_process[:]
				else:
					proj_ids = None
				proj_ids = wrap_mpi_bcast(proj_ids, 0, mpi_subcomm)
				if len(proj_ids) == 0:
					break
				if mpi_subrank < len(proj_ids):
					# -------- alignment
					im = proj_ids[mpi_subrank]
					peak, pixel_error, checked_refs, iref = shc(data[im], refrings, numr, xrng[N_step], yrng[N_step], step[N_step], an[N_step])
					# -------- gather results to root
					vector_assigned_refs = wrap_mpi_gatherv([iref], 0, mpi_subcomm)
					vector_previousmax   = wrap_mpi_gatherv([data[im].get_attr("previousmax")], 0, mpi_subcomm)
					vector_xformprojs    = wrap_mpi_gatherv([data[im].get_attr("xform.projection")], 0, mpi_subcomm)
					vector_pixel_error   = wrap_mpi_gatherv([pixel_error], 0, mpi_subcomm)
					vector_checked_ref   = wrap_mpi_gatherv([checked_refs], 0, mpi_subcomm)
				else:
					# -------- no projection assigned, send to root empty lists
					vector_assigned_refs = wrap_mpi_gatherv([], 0, mpi_subcomm)
					vector_previousmax   = wrap_mpi_gatherv([], 0, mpi_subcomm)
					vector_xformprojs    = wrap_mpi_gatherv([], 0, mpi_subcomm)
					vector_pixel_error   = wrap_mpi_gatherv([], 0, mpi_subcomm)
					vector_checked_ref   = wrap_mpi_gatherv([], 0, mpi_subcomm)
				# -------- merge results
				if mpi_subrank == 0:
					used_refs = set()
					for i in xrange(len(vector_assigned_refs)):
						ir = vector_assigned_refs[i]
						if ir in used_refs:
							# reference is already used - cancel all changes
							vector_previousmax[i] = data[proj_ids[i]].get_attr("previousmax")
							vector_xformprojs[i]  = data[proj_ids[i]].get_attr("xform.projection")
						else:
							used_refs.add(ir)
							proj_ids_to_process.remove(proj_ids[i])
							pixer.append(vector_pixel_error[i])
							number_of_checked_refs += vector_checked_ref[i]
					used_refs = list(used_refs)
					used_refs.sort(reverse=True)
				else:
					used_refs = None
				# ------- broadcast results
				used_refs = wrap_mpi_bcast(used_refs, 0, mpi_subcomm)
				vector_previousmax = wrap_mpi_bcast(vector_previousmax, 0, mpi_subcomm)
				vector_xformprojs  = wrap_mpi_bcast(vector_xformprojs, 0, mpi_subcomm)
				# ------- delete used references
				for ir in used_refs:  del refrings[ir]
				# ------- set projections parameters
				for i in xrange(len(vector_previousmax)):
					data[proj_ids[i]].set_attr("previousmax", vector_previousmax[i])
					data[proj_ids[i]].set_attr("xform.projection", vector_xformprojs[i])
			#=========================================================================
			mpi_barrier(mpi_comm)
			if myid == main_node:
				log.add("Time of alignment = %f\n"%(time()-start_time))

			#=========================================================================
			#output pixel errors, check stop criterion
			if mpi_subrank == 0:
				all_pixer          = wrap_mpi_gatherv(pixer, 0, mpi_comm)
				total_checked_refs = wrap_mpi_gatherv([number_of_checked_refs], main_node, mpi_comm)
			else:
				all_pixer          = wrap_mpi_gatherv([], 0, mpi_comm)
				total_checked_refs = wrap_mpi_gatherv([], main_node, mpi_comm)
			if myid == main_node:
				total_checked_refs = sum(total_checked_refs)
				lhist = 20
				region, histo = hist_list(all_pixer, lhist)
				log.add("=========================")
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d"%(region[lhx], histo[lhx])
					log.add(msg)
				temp = 0
				for i in all_pixer:
					if i < 1.0: temp += 1
				percent_of_pixerr_below_one = (temp * 1.0) / (total_nima * number_of_runs)
				orient_and_shuffle = ( percent_of_pixerr_below_one > 0.3 )  #  TODO - parameter ?
				#terminate          = ( percent_of_pixerr_below_one > 0.9 )  #  TODO - parameter ?
				log.add("=========================")
				log.add("Percent of positions with pixel error below 1.0 = ", (int(percent_of_pixerr_below_one*100)), "%","   Shuffling: ",orient_and_shuffle)
			orient_and_shuffle = wrap_mpi_bcast(orient_and_shuffle, 0, mpi_comm)
			#=========================================================================

			#=========================================================================
			# centering
			if center == -1 and sym[0] == 'c':
				from utilities      import estimate_3D_center_MPI, rotate_3D_shift
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data[image_start:image_end], total_nima, mpi_subrank, mpi_subsize, 0, mpi_comm=mpi_subcomm) #estimate_3D_center_MPI(data, number_of_runs*total_nima, myid, number_of_proc, main_node, mpi_comm=mpi_comm)
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					log.add(msg)
				if int(sym[1]) > 1:
					cs[0] = cs[1] = 0.0
					if myid == main_node:
						log.add("For symmetry group cn (n>1), we only center the volume in z-direction\n")
				cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, mpi_subcomm)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(data, cs)
			#=========================================================================

			mpi_barrier(mpi_comm)
			if myid == main_node:
				start_time = time()
			#  Temporary - write out params
			params = []
			for im in data:
				phi, theta, psi, sx, sy = get_params_proj(im)
				params.append([phi, theta, psi, sx, sy])
			params_0 = wrap_mpi_bcast(params, mpi_subroots[0], mpi_comm)
			if mpi_subrank == 0:
				from utilities import write_text_row
				write_text_row(params, "qparams%04d%04d.hdf"%(myid,total_iter) )


			#=========================================================================
			# volume reconstruction
			mpi_barrier(mpi_comm)
			if myid == main_node:
				start_time = time()
			vol = volume_reconstruction(data[image_start:image_end], ali3d_options, mpi_subcomm)
			if mpi_subrank == 0:  vol.write_image("qvolf%04d%04d.hdf"%(myid,total_iter))
			# log
			if myid == main_node:
				log.add("3D reconstruction time = %f\n"%(time()-start_time))
				start_time = time()

			#=========================================================================
			if orient_and_shuffle:
				terminate = 0
				params = []
				for im in data:
					phi, theta, psi, sx, sy = get_params_proj(im)
					params.append([phi, theta, psi, sx, sy])

				# ------ orientation - begin
				params_0 = wrap_mpi_bcast(params, mpi_subroots[0], mpi_comm)
				if mpi_subrank == 0:
					if(sym[0] == "d"):
						reduce_dsym_angles(params_0, sym)
						reduce_dsym_angles(params, sym)
					subset_thr, subset_min, avg_diff_per_image = find_common_subset_3([params_0, params], 2.0, len(params)/3, sym)
					if len(subset_thr) < len(subset_min):
						subset = subset_min
					else:
						subset = subset_thr
					if(sym[0] != "d"):  orient_params([params_0, params], subset)
				if(sym[0] != "d"):  params = wrap_mpi_bcast(params, 0, mpi_subcomm)
				# ------ orientation - end
				
				# ------ Compute L2s and gather to the root
				if mpi_subrank == 0:
					L2 = vol.cmp("dot", vol, dict(negative = 0, mask = model_circle(last_ring, nx, nx, nx)))
				if myid == 0:
					all_L2s = []
					for sr in mpi_subroots:
						if sr == myid:
							all_L2s.append(L2)
						else:
							all_L2s.append(wrap_mpi_recv(sr, mpi_comm))
				else:
					if mpi_subrank == 0:
						wrap_mpi_send(L2, 0, mpi_comm)				

				# ------ gather parameters to root
				if myid == 0:
					all_params = []
					for sr in mpi_subroots:
						if sr == myid:
							all_params.append(params)
						else:
							all_params.append(wrap_mpi_recv(sr, mpi_comm))
				else:
					if mpi_subrank == 0:
						wrap_mpi_send(params, 0, mpi_comm)
				
				# ---------------------------------

				#  Add params to GA, do mutations and send back
				if myid == 0:
					#all_params = shuffle_configurations(all_params)
					for i in xrange(number_of_runs):
						GA.append([all_L2s[i],params[i]])
					GA.sort(reverse=True)
					GA = GA[:number_of_runs]
					#  ---  Stopping criterion
					from statistics import table_stat
					from math import sqrt
					q1,q2,q3,q4 = table_stat([GA[i][0] for i in xrange(number_of_runs)])
					# Terminate if Vvariation of L2 norms less than 10% of their average
					terminate = sqrt(max(q2,0.0))/q1 <0.1

					if not terminate:
						#  Now do the mutation

						#  Put mutated params on one list
						all_params = []
						for i in xrange(number_of_runs):
							all_params.append(GA[i][1])

				terminate = wrap_mpi_bcast(terminate, main_node, mpi_comm)


				if not terminate:
					# Send params back
					if myid == 0:
						for i in xrange(number_of_runs):
							sr = mpi_subroots[i]
							if sr == myid:
								params = all_params[i]
							else:
								wrap_mpi_send(all_params[i], sr, mpi_comm)
					else:
						if mpi_subrank == 0:
							params = wrap_mpi_recv(0, mpi_comm)

					params = wrap_mpi_bcast(params, 0, mpi_subcomm)
					for i in xrange(nima):
						set_params_proj(data[i], params[i])

					#=========================================================================
					# volume reconstruction
					mpi_barrier(mpi_comm)
					if myid == main_node:
						start_time = time()
					vol = volume_reconstruction(data[image_start:image_end], ali3d_options, mpi_subcomm)
					if mpi_subrank == 0:  vol.write_image("qmutatedvolf%04d%04d.hdf"%(myid,total_iter))
					# log
					if myid == main_node:
						log.add("3D reconstruction time = %f\n"%(time()-start_time))
						start_time = time()
					#=========================================================================

				mpi_barrier(mpi_comm)
				if myid == main_node:
					log.add("Time of orientation and mutations = %f\n"%(time()-start_time))
					start_time = time()

	#=========================================================================
	# gather parameters to subroot-s
	params = []
	previousmax = []
	for im in data:
		t = get_params_proj(im)
		p = im.get_attr("previousmax")
		params.append( [t[0], t[1], t[2], t[3], t[4]] )
		previousmax.append(p)
	assert(len(params) == nima)

	# gather data to main root
	if mpi_subrank == 0:
		vol         = wrap_mpi_gatherv([vol], 0, mpi_comm)
		params      = wrap_mpi_gatherv([params], 0, mpi_comm)
		previousmax = wrap_mpi_gatherv([previousmax], 0, mpi_comm)
	else:
		vol         = wrap_mpi_gatherv([], 0, mpi_comm)
		params      = wrap_mpi_gatherv([], 0, mpi_comm)
		previousmax = wrap_mpi_gatherv([], 0, mpi_comm)

	mpi_comm_free(mpi_subcomm)
	
	
	if myid == main_node: 
		log.add("Finish ali3d_multishc")
		if(sym[0] == "d"):  reduce_dsym_angles(params, sym)
		return params, vol, previousmax
	else:
		return None, None, None  # results for the other processes

"""


# parameters: list of (all) projections | reference volume | ...
def ali3d_multishc(stack, ref_vol, ali3d_options, mpi_comm = None, log = None, number_of_runs=2 ):

	from alignment    import Numrinit, prepare_refrings, proj_ali_incore_local, shc
	from utilities    import model_circle, get_input_from_string, get_params_proj, set_params_proj, wrap_mpi_gatherv, wrap_mpi_bcast, wrap_mpi_send, wrap_mpi_recv
	from mpi          import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier, mpi_comm_split, mpi_comm_free
	from projection   import prep_vol
	from statistics   import hist_list
	from applications import MPI_start_end
	from filter       import filt_ctf
	from global_def   import Util, ERROR
	from time         import time
	from random       import shuffle

	ir     = ali3d_options.ir
	rs     = ali3d_options.rs
	ou     = ali3d_options.ou
	xr     = ali3d_options.xr
	yr     = ali3d_options.yr
	ts     = ali3d_options.ts
	an     = ali3d_options.an
	sym    = ali3d_options.sym
	sym = sym[0].lower() + sym[1:]
	delta  = ali3d_options.delta
	center = ali3d_options.center
	maxit  = ali3d_options.maxit
	CTF    = ali3d_options.CTF
	ref_a  = ali3d_options.ref_a

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	if log == None:
		from logger import Logger
		log = Logger()

	number_of_proc = mpi_comm_size(mpi_comm)
	myid           = mpi_comm_rank(mpi_comm)
	main_node = 0

	if myid == main_node:
		log.add("Start ali3d_multishc")

	if number_of_proc < number_of_runs:
		ERROR("number_of_proc < number_of_runs","ali3d_multishc")
	
	mpi_subcomm = mpi_comm_split(mpi_comm, myid % number_of_runs, myid / number_of_runs)
	mpi_subrank = mpi_comm_rank(mpi_subcomm)
	mpi_subsize = mpi_comm_size(mpi_subcomm)
	mpi_subroots = range(number_of_runs)

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

	vol = ref_vol
	nx      = vol.get_xsize()
	if last_ring < 0:
		last_ring = int(nx/2) - 2

	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	mask2D  = model_circle(last_ring,nx,nx) - model_circle(first_ring,nx,nx)

	if myid == main_node:
		list_of_particles = range(len(stack))
		total_nima = len(list_of_particles)
	else:
		list_of_particles = None
		total_nima = None
	total_nima = wrap_mpi_bcast(total_nima, main_node, mpi_comm)
	list_of_particles = wrap_mpi_bcast(list_of_particles, main_node, mpi_comm)
	nima = len(list_of_particles)

	image_start, image_end = MPI_start_end(total_nima, mpi_subsize, mpi_subrank)

	data = [ stack[im] for im in list_of_particles ]
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		ctf_applied = data[im].get_attr_default('ctf_applied', 0)
		if CTF and ctf_applied == 0:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)

	cs = [0.0]*3
	total_iter = 0

	orient_and_shuffle = False

	# do the projection matching
	for N_step in xrange(lstp):
		
		terminate = 0
		Iter = 0
		while Iter < max_iter and terminate == 0:

			Iter += 1
			total_iter += 1

			mpi_barrier(mpi_comm)
			if myid == main_node:
				log.add("ITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step]))
				start_time = time()

			#=========================================================================
			# build references
			volft, kb = prep_vol(vol)
			refrings = prepare_refrings(volft, kb, nx, delta[N_step], ref_a, sym, numr, MPI=mpi_subcomm)
			del volft, kb
			all_ref_dirs = []
			for r in refrings:
				all_ref_dirs.append( [r.get_attr("phi"), r.get_attr("theta")] )
			#=========================================================================

			mpi_barrier(mpi_comm)
			if myid == main_node:
				log.add("Time to prepare rings: %f\n" % (time()-start_time))
				start_time = time()

			#=========================================================================
			if total_iter == 1 or orient_and_shuffle:
				# adjust params to references, calculate psi+shifts, calculate previousmax
				for im in xrange(nima):
					stable = data[im].get_attr_default("stable", 0)
					if stable == 0:
						data[im].set_attr("previousmax", -1.0e23)
						data[im].set_attr("stable", 1)
					else:
						peak, temp = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step], delta[N_step]*0.7 )
						data[im].set_attr("previousmax", peak)
				if myid == main_node:
					log.add("Time to calculate first psi+shifts+previousmax: %f\n" % (time()-start_time))
					start_time = time()
			#=========================================================================

			mpi_barrier(mpi_comm)
			if myid == main_node: start_time = time()
			#=========================================================================
			# alignment
			if mpi_subrank == 0:
				pixer = []
				number_of_checked_refs = 0
				proj_ids_to_process = range(nima)
				shuffle(proj_ids_to_process)
			while True:
				# --------- broadcast projections ids
				if mpi_subrank == 0:
					if len(proj_ids_to_process) >= mpi_subsize:
						proj_ids = proj_ids_to_process[:(mpi_subsize)]
					else:
						proj_ids = proj_ids_to_process[:]
				else:
					proj_ids = None
				proj_ids = wrap_mpi_bcast(proj_ids, 0, mpi_subcomm)
				if len(proj_ids) == 0:
					break
				if mpi_subrank < len(proj_ids):
					# -------- alignment
					im = proj_ids[mpi_subrank]
					peak, pixel_error, checked_refs, iref = shc(data[im], refrings, numr, xrng[N_step], yrng[N_step], step[N_step], an[N_step])
					# -------- gather results to root
					vector_assigned_refs = wrap_mpi_gatherv([iref], 0, mpi_subcomm)
					vector_previousmax   = wrap_mpi_gatherv([data[im].get_attr("previousmax")], 0, mpi_subcomm)
					vector_xformprojs    = wrap_mpi_gatherv([data[im].get_attr("xform.projection")], 0, mpi_subcomm)
					vector_pixel_error   = wrap_mpi_gatherv([pixel_error], 0, mpi_subcomm)
					vector_checked_ref   = wrap_mpi_gatherv([checked_refs], 0, mpi_subcomm)
				else:
					# -------- no projection assigned, send to root empty lists
					vector_assigned_refs = wrap_mpi_gatherv([], 0, mpi_subcomm)
					vector_previousmax   = wrap_mpi_gatherv([], 0, mpi_subcomm)
					vector_xformprojs    = wrap_mpi_gatherv([], 0, mpi_subcomm)
					vector_pixel_error   = wrap_mpi_gatherv([], 0, mpi_subcomm)
					vector_checked_ref   = wrap_mpi_gatherv([], 0, mpi_subcomm)
				# -------- merge results
				if mpi_subrank == 0:
					used_refs = set()
					for i in xrange(len(vector_assigned_refs)):
						ir = vector_assigned_refs[i]
						if ir in used_refs:
							# reference is already used - cancel all changes
							vector_previousmax[i] = data[proj_ids[i]].get_attr("previousmax")
							vector_xformprojs[i]  = data[proj_ids[i]].get_attr("xform.projection")
						else:
							used_refs.add(ir)
							proj_ids_to_process.remove(proj_ids[i])
							pixer.append(vector_pixel_error[i])
							number_of_checked_refs += vector_checked_ref[i]
					used_refs = list(used_refs)
					used_refs.sort(reverse=True)
				else:
					used_refs = None
				# ------- broadcast results
				used_refs = wrap_mpi_bcast(used_refs, 0, mpi_subcomm)
				vector_previousmax = wrap_mpi_bcast(vector_previousmax, 0, mpi_subcomm)
				vector_xformprojs  = wrap_mpi_bcast(vector_xformprojs, 0, mpi_subcomm)
				# ------- delete used references
				for ir in used_refs:  del refrings[ir]
				# ------- set projections parameters
				for i in xrange(len(vector_previousmax)):
					data[proj_ids[i]].set_attr("previousmax", vector_previousmax[i])
					data[proj_ids[i]].set_attr("xform.projection", vector_xformprojs[i])
			#=========================================================================
			mpi_barrier(mpi_comm)
			if myid == main_node:
				log.add("Time of alignment = %f\n"%(time()-start_time))

			#=========================================================================
			#output pixel errors, check stop criterion
			if mpi_subrank == 0:
				all_pixer          = wrap_mpi_gatherv(pixer, 0, mpi_comm)
				total_checked_refs = wrap_mpi_gatherv([number_of_checked_refs], main_node, mpi_comm)
			else:
				all_pixer          = wrap_mpi_gatherv([], 0, mpi_comm)
				total_checked_refs = wrap_mpi_gatherv([], main_node, mpi_comm)
			terminate = 0
			if myid == main_node:
				total_checked_refs = sum(total_checked_refs)
				lhist = 20
				region, histo = hist_list(all_pixer, lhist)
				log.add("=========================")
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d"%(region[lhx], histo[lhx])
					log.add(msg)
				temp = 0
				for i in all_pixer:
					if i < 1.0: temp += 1
				percent_of_pixerr_below_one = (temp * 1.0) / (total_nima * number_of_runs)
				orient_and_shuffle = ( percent_of_pixerr_below_one > 0.3 )  #  TODO - parameter ?
				terminate          = ( percent_of_pixerr_below_one > 0.9 )  #  TODO - parameter ?
				log.add("=========================")
				log.add("Percent of positions with pixel error below 1.0 = ", (int(percent_of_pixerr_below_one*100)), "%","   Shuffling: ",orient_and_shuffle)
			terminate = wrap_mpi_bcast(terminate, main_node, mpi_comm)
			orient_and_shuffle = wrap_mpi_bcast(orient_and_shuffle, 0, mpi_comm)
			#=========================================================================

			#=========================================================================
			# centering
			if center == -1 and sym[0] == 'c':
				from utilities      import estimate_3D_center_MPI, rotate_3D_shift
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data[image_start:image_end], total_nima, mpi_subrank, mpi_subsize, 0, mpi_comm=mpi_subcomm) #estimate_3D_center_MPI(data, number_of_runs*total_nima, myid, number_of_proc, main_node, mpi_comm=mpi_comm)
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					log.add(msg)
				if int(sym[1]) > 1:
					cs[0] = cs[1] = 0.0
					if myid == main_node:
						log.add("For symmetry group cn (n>1), we only center the volume in z-direction\n")
				cs = mpi_bcast(cs, 3, MPI_FLOAT, main_node, mpi_subcomm)
				cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
				rotate_3D_shift(data, cs)
			#=========================================================================

			mpi_barrier(mpi_comm)
			if myid == main_node:
				start_time = time()
			#  Temporary - write out params
			params = []
			for im in data:
				phi, theta, psi, sx, sy = get_params_proj(im)
				params.append([phi, theta, psi, sx, sy])
			params_0 = wrap_mpi_bcast(params, mpi_subroots[0], mpi_comm)
			if mpi_subrank == 0:
				from utilities import write_text_row
				write_text_row(params, "qparams%04d%04d.hdf"%(myid,total_iter) )
			#=========================================================================
			if orient_and_shuffle and not terminate:
				params = []
				for im in data:
					phi, theta, psi, sx, sy = get_params_proj(im)
					params.append([phi, theta, psi, sx, sy])

				# ------ orientation - begin
				params_0 = wrap_mpi_bcast(params, mpi_subroots[0], mpi_comm)
				if mpi_subrank == 0:
					if(sym[0] == "d"):
						reduce_dsym_angles(params_0, sym)
						reduce_dsym_angles(params, sym)
					subset_thr, subset_min, avg_diff_per_image = find_common_subset_3([params_0, params], 2.0, len(params)/3, sym)
					if len(subset_thr) < len(subset_min):
						subset = subset_min
					else:
						subset = subset_thr
					if(sym[0] != "d"):  orient_params([params_0, params], subset)
				if(sym[0] != "d"):  params = wrap_mpi_bcast(params, 0, mpi_subcomm)
				# ------ orientation - end

				# ------ gather parameters to root
				if myid == 0:
					all_params = []
					for sr in mpi_subroots:
						if sr == myid:
							all_params.append(params)
						else:
							all_params.append(wrap_mpi_recv(sr, mpi_comm))
				else:
					if mpi_subrank == 0:
						wrap_mpi_send(params, 0, mpi_comm)
				# ---------------------------------

				if myid == 0:
					all_params = shuffle_configurations(all_params)

				if myid == 0:
					for i in xrange(number_of_runs):
						sr = mpi_subroots[i]
						if sr == myid:
							params = all_params[i]
						else:
							wrap_mpi_send(all_params[i], sr, mpi_comm)
				else:
					if mpi_subrank == 0:
						params = wrap_mpi_recv(0, mpi_comm)

				params = wrap_mpi_bcast(params, 0, mpi_subcomm)
				for i in xrange(nima):
					set_params_proj(data[i], params[i])

				mpi_barrier(mpi_comm)
				if myid == main_node:
					log.add("Time of orientation and shuffling = %f\n"%(time()-start_time))
					start_time = time()

			#=========================================================================
			# volume reconstruction
			mpi_barrier(mpi_comm)
			if myid == main_node:
				start_time = time()
			vol = volume_reconstruction(data[image_start:image_end], ali3d_options, mpi_subcomm)
			if mpi_subrank == 0:  vol.write_image("qvolf%04d%04d.hdf"%(myid,total_iter))
			# log
			if myid == main_node:
				log.add("3D reconstruction time = %f\n"%(time()-start_time))
				start_time = time()
			#=========================================================================

	#=========================================================================
	# gather parameters to subroot-s
	params = []
	previousmax = []
	for im in data:
		t = get_params_proj(im)
		p = im.get_attr("previousmax")
		params.append( [t[0], t[1], t[2], t[3], t[4]] )
		previousmax.append(p)
	assert(len(params) == nima)

	# gather data to main root
	if mpi_subrank == 0:
		vol         = wrap_mpi_gatherv([vol], 0, mpi_comm)
		params      = wrap_mpi_gatherv([params], 0, mpi_comm)
		previousmax = wrap_mpi_gatherv([previousmax], 0, mpi_comm)
	else:
		vol         = wrap_mpi_gatherv([], 0, mpi_comm)
		params      = wrap_mpi_gatherv([], 0, mpi_comm)
		previousmax = wrap_mpi_gatherv([], 0, mpi_comm)

	mpi_comm_free(mpi_subcomm)
	
	
	if myid == main_node: 
		log.add("Finish ali3d_multishc")
		if(sym[0] == "d"):  reduce_dsym_angles(params, sym)		
		return params, vol, previousmax
	else:
		return None, None, None  # results for the other processes



"""

def shc_multi(data, refrings, numr, xrng, yrng, step, an, number_of_runs, finfo=None):
	from utilities    import compose_transform2
	from math         import cos, pi
	from EMAN2 import Vec2f, Transform
	from global_def import Util

	ID = data.get_attr("ID")

	mode = "F"
	nx   = data.get_xsize()
	ny   = data.get_ysize()
	#  center is in SPIDER convention
	cnx  = nx//2 + 1
	cny  = ny//2 + 1

	ant = cos(an*pi/180.0)
	#phi, theta, psi, sxo, syo = get_params_proj(data)
	t1 = data.get_attr("xform.projection")
	dp = t1.get_params("spider")
	if finfo:
		finfo.write("Image id: %6d\n"%(ID))
		finfo.write("Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"]))
		finfo.flush()

	#[ang, sxs, sys, mirror, iref, peak, checked_refs] = Util.shc(data, refrings, xrng, yrng, step, ant, mode, numr, cnx+dp["tx"], cny+dp["ty"])
	peaks = Util.shc_multipeaks(data, refrings, xrng, yrng, step, ant, mode, numr, cnx+dp["tx"], cny+dp["ty"], number_of_runs)
	peaks_count = len(peaks) / 7
	pixel_error = 0.0
	number_of_checked_refs = 0
	peak = 0.0
	for i in xrange(peaks_count):
		ang    = peaks[i*7+0]
		sxs    = peaks[i*7+1]
		sys    = peaks[i*7+2]
		mirror = peaks[i*7+3]
		iref   = int(peaks[i*7+4])
		peak   = peaks[i*7+5]
		checked_refs = int(peaks[i*7+6])
		number_of_checked_refs += checked_refs
		#[ang,sxs,sys,mirror,peak,numref] = apmq_local(projdata[imn], ref_proj_rings, xrng, yrng, step, ant, mode, numr, cnx-sxo, cny-syo)
		#ang = (ang+360.0)%360.0

		# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
		# What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
		angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1, -ang, 0.0, 0.0, 1)
		if  mirror:
			phi   = (refrings[iref].get_attr("phi")+540.0)%360.0
			theta = 180.0-refrings[iref].get_attr("theta")
			psi   = (540.0-refrings[iref].get_attr("psi")+angb)%360.0
			s2x   = sxb - dp["tx"]
			s2y   = syb - dp["ty"]
		else:
			phi   = refrings[iref].get_attr("phi")
			theta = refrings[iref].get_attr("theta")
			psi   = (refrings[iref].get_attr("psi")+angb+360.0)%360.0
			s2x   = sxb - dp["tx"]
			s2y   = syb - dp["ty"]

		t2 = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
		t2.set_trans(Vec2f(-s2x, -s2y))
		if i == 0:
			data.set_attr("xform.projection", t2)
		else:
			data.set_attr("xform.projection" + str(i), t2)
		from pixel_error import max_3D_pixel_error
		pixel_error += max_3D_pixel_error(t1, t2, numr[-3])
		if finfo:
			finfo.write( "New parameters: %9.4f %9.4f %9.4f %9.4f %9.4f %10.5f  %11.3e\n\n" %(phi, theta, psi, s2x, s2y, peak, pixel_error))
			finfo.flush()

	# remove old xform.projection
	i = max(peaks_count, 1)
	while data.has_attr("xform.projection" + str(i)):
		data.del_attr("xform.projection" + str(i))
		i += 1

# -------- remove weights
# 	data.del_attr("weight")
# 	for i in xrange(1, 50):
# 		if data.has_attr("weight" + str(i)):
# 			data.del_attr("weight" + str(i))
	
	return peak, (pixel_error / 7), (number_of_checked_refs / 7), peaks_count


# parameters: list of (all) projections | reference volume | ...
def ali3d_multishc_2(stack, ref_vol, ali3d_options, mpi_comm = None, log = None ):

	from alignment       import Numrinit, prepare_refrings, proj_ali_incore_local, shc
	from utilities       import model_circle, get_input_from_string, get_params_proj, wrap_mpi_gatherv, wrap_mpi_bcast
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from projection      import prep_vol
	from statistics      import hist_list
	from applications    import MPI_start_end
	from filter         import filt_ctf
	from global_def import Util
	from time import time

	ir     = ali3d_options.ir
	rs     = ali3d_options.rs
	ou     = ali3d_options.ou
	xr     = ali3d_options.xr
	yr     = ali3d_options.yr
	ts     = ali3d_options.ts
	an     = ali3d_options.an
	sym    = ali3d_options.sym
	sym = sym[0].lower() + sym[1:]
	delta  = ali3d_options.delta
	center = ali3d_options.center
	maxit  = ali3d_options.maxit
	CTF    = ali3d_options.CTF
	ref_a  = ali3d_options.ref_a

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	if log == None:
		from logger import Logger
		log = Logger()

	number_of_proc = mpi_comm_size(mpi_comm)
	myid           = mpi_comm_rank(mpi_comm)
	main_node = 0

	if myid == main_node:
		log.add("Start ali3d_multishc_2")

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

	vol = ref_vol
	nx      = vol.get_xsize()
	if last_ring < 0:	last_ring = int(nx/2) - 2

	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	mask2D  = model_circle(last_ring,nx,nx) - model_circle(first_ring,nx,nx)

	if myid == main_node:
		list_of_particles = range(len(stack))
		total_nima = len(list_of_particles)
	else:
		list_of_particles = None
		total_nima = None
	total_nima = wrap_mpi_bcast(total_nima, main_node, mpi_comm)
	list_of_particles = wrap_mpi_bcast(list_of_particles, main_node, mpi_comm)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)

	data = [ stack[im] for im in list_of_particles ]
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		ctf_applied = data[im].get_attr_default('ctf_applied', 0)
		if CTF and ctf_applied == 0:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)

	pixer = [0.0]*nima
	par_r = [[] for im in list_of_particles ]
	cs = [0.0]*3
	total_iter = 0
	# do the projection matching
	for N_step in xrange(lstp):
		
		terminate = 0
		Iter = 0
		while Iter < max_iter and terminate == 0:

			Iter += 1
			total_iter += 1

			mpi_barrier(mpi_comm)
			if myid == main_node:
				log.add("ITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step]))
				start_time = time()

			#=========================================================================
			# build references
			volft, kb = prep_vol(vol)
			refrings = prepare_refrings(volft, kb, nx, delta[N_step], ref_a, sym, numr, MPI=mpi_comm)
			del volft, kb
			#=========================================================================

			if myid == main_node:
				log.add("Time to prepare rings: %f\n" % (time()-start_time))
				start_time = time()
			
			#=========================================================================
			#if total_iter == 1:
			# adjust params to references, calculate psi+shifts, calculate previousmax
			for im in xrange(nima):
				stable = data[im].get_attr_default("stable", 0)
				if stable == 0:
					data[im].set_attr("previousmax", -1.0e23)
					data[im].set_attr("stable", 1)
				else:
					peak, pixer[im] = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],1.0)
					data[im].set_attr("previousmax", peak)
			if myid == main_node:
				log.add("Time to calculate first psi+shifts+previousmax: %f\n" % (time()-start_time))
				start_time = time()
			#=========================================================================

			mpi_barrier(mpi_comm)
			if myid == main_node:
				start_time = time()
			#=========================================================================
			# alignment
			number_of_checked_refs = 0
			for im in xrange(nima):
				#peak, pixer[im], checked_refs, number_of_peaks = shc_multi(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step], number_of_runs=number_of_runs)
				peak, pixer[im], checked_refs, iref = shc(data[im], refrings, numr, xrng[N_step], yrng[N_step], step[N_step], an[N_step])
				number_of_checked_refs += checked_refs
			#=========================================================================
			mpi_barrier(mpi_comm)
			if myid == main_node:
				print  data[0].get_attr_dict()
				log.add("Time of alignment = %f\n"%(time()-start_time))
				start_time = time()

			#=========================================================================
			#output pixel errors, check stop criterion
			all_pixer = wrap_mpi_gatherv(pixer, 0, mpi_comm)
			total_checked_refs = wrap_mpi_gatherv([number_of_checked_refs], main_node, mpi_comm)
			terminate = 0
			if myid == main_node:
				total_checked_refs = sum(total_checked_refs)
				lhist = 20
				region, histo = hist_list(all_pixer, lhist)
				log.add("=========================")
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d"%(region[lhx], histo[lhx])
					log.add(msg)
				if (max(all_pixer) < 0.5) and (sum(all_pixer)/total_nima < 0.05):
					terminate = 1
			terminate = wrap_mpi_bcast(terminate, main_node, mpi_comm)
			#=========================================================================

			#=========================================================================
			# centering
			if center == -1 and sym[0] == 'c':
				from utilities      import estimate_3D_center_MPI, rotate_3D_shift
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

			#=========================================================================
			# volume reconstruction
			mpi_barrier(mpi_comm)
			if myid == main_node:
				start_time = time()
			vol = volume_reconstruction(data, ali3d_options, mpi_comm)
			# log
			if myid == main_node:
				log.add("3D reconstruction time = %f\n"%(time()-start_time))
				start_time = time()
			#=========================================================================

	#=========================================================================
	# gather parameters
	params = []
	previousmax = []
	for im in data:
		t = get_params_proj(im)
		p = im.get_attr("previousmax")
		params.append( [t[0], t[1], t[2], t[3], t[4]] )
		previousmax.append(p)
	assert(nima == len(params))
	params = wrap_mpi_gatherv(params, 0, mpi_comm)
	if myid == 0:
		assert(total_nima == len(params))
	previousmax = wrap_mpi_gatherv(previousmax, 0, mpi_comm)

	par_r = wrap_mpi_gatherv(par_r, 0, mpi_comm)

	if myid == main_node: 
		log.add("Finish ali3d_multishc_2")
		if(sym[0] == "d"):  reduce_dsym_angles(params, sym)
		return params, vol, previousmax, par_r
	else:
		return None, None, None, None  # results for the other processes

"""

# parameters: list of (all) projections | reference volume | ...
def ali3d_multishc_2(stack, ref_vol, ali3d_options, mpi_comm = None, log = None, number_of_runs=2 ):

	from alignment       import Numrinit, prepare_refrings, proj_ali_incore_local
	from utilities       import model_circle, get_input_from_string, get_params_proj, wrap_mpi_gatherv, wrap_mpi_bcast
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from projection      import prep_vol
	from statistics      import hist_list
	from applications    import MPI_start_end
	from filter         import filt_ctf
	from global_def import Util
	from time import time

	ir     = ali3d_options.ir
	rs     = ali3d_options.rs
	ou     = ali3d_options.ou
	xr     = ali3d_options.xr
	yr     = ali3d_options.yr
	ts     = ali3d_options.ts
	an     = ali3d_options.an
	sym    = ali3d_options.sym
	sym = sym[0].lower() + sym[1:]
	delta  = ali3d_options.delta
	center = ali3d_options.center
	maxit  = ali3d_options.maxit
	CTF    = ali3d_options.CTF
	ref_a  = ali3d_options.ref_a

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	if log == None:
		from logger import Logger
		log = Logger()

	number_of_proc = mpi_comm_size(mpi_comm)
	myid           = mpi_comm_rank(mpi_comm)
	main_node = 0

	if myid == main_node:
		log.add("Start ali3d_multishc_2")

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

	vol = ref_vol
	nx      = vol.get_xsize()
	if last_ring < 0:	last_ring = int(nx/2) - 2

	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	mask2D  = model_circle(last_ring,nx,nx) - model_circle(first_ring,nx,nx)

	if myid == main_node:
		list_of_particles = range(len(stack))
		total_nima = len(list_of_particles)
	else:
		list_of_particles = None
		total_nima = None
	total_nima = wrap_mpi_bcast(total_nima, main_node, mpi_comm)
	list_of_particles = wrap_mpi_bcast(list_of_particles, main_node, mpi_comm)

	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)

	data = [ stack[im] for im in list_of_particles ]
	for im in xrange(nima):
		data[im].set_attr('ID', list_of_particles[im])
		ctf_applied = data[im].get_attr_default('ctf_applied', 0)
		if CTF and ctf_applied == 0:
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)

	pixer = [0.0]*nima
	par_r = [[] for im in list_of_particles ]
	cs = [0.0]*3
	total_iter = 0
	# do the projection matching
	for N_step in xrange(lstp):
		
		terminate = 0
		Iter = 0
		while Iter < max_iter and terminate == 0:

			Iter += 1
			total_iter += 1

			mpi_barrier(mpi_comm)
			if myid == main_node:
				log.add("ITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f\n"%(total_iter, Iter, delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step]))
				start_time = time()

			#=========================================================================
			# build references
			volft, kb = prep_vol(vol)
			refrings = prepare_refrings(volft, kb, nx, delta[N_step], ref_a, sym, numr, MPI=mpi_comm)
			del volft, kb
			#=========================================================================

			if myid == main_node:
				log.add("Time to prepare rings: %f\n" % (time()-start_time))
				start_time = time()
			
			#=========================================================================
			#if total_iter == 1:
			# adjust params to references, calculate psi+shifts, calculate previousmax
			for im in xrange(nima):
				stable = data[im].get_attr_default("stable", 0)
				if stable == 0:
					data[im].set_attr("previousmax", -1.0e23)
					data[im].set_attr("stable", 1)
				else:
					peak, pixer[im] = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],1.0)
					data[im].set_attr("previousmax", peak)
			if myid == main_node:
				log.add("Time to calculate first psi+shifts+previousmax: %f\n" % (time()-start_time))
				start_time = time()
			#=========================================================================

			mpi_barrier(mpi_comm)
			if myid == main_node:
				start_time = time()
			#=========================================================================
			# alignment
			number_of_checked_refs = 0
			for im in xrange(nima):
				peak, pixer[im], checked_refs, number_of_peaks = shc_multi(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step], number_of_runs=number_of_runs)
				number_of_checked_refs += checked_refs
				par_r[im].append(number_of_peaks)
				print  myid,im,number_of_peaks
			#=========================================================================
			mpi_barrier(mpi_comm)
			if myid == main_node:
				print  data[0].get_attr_dict()
				log.add("Time of alignment = %f\n"%(time()-start_time))
				start_time = time()

			#=========================================================================
			#output pixel errors, check stop criterion
			all_pixer = wrap_mpi_gatherv(pixer, 0, mpi_comm)
			total_checked_refs = wrap_mpi_gatherv([number_of_checked_refs], main_node, mpi_comm)
			terminate = 0
			if myid == main_node:
				total_checked_refs = sum(total_checked_refs)
				lhist = 20
				region, histo = hist_list(all_pixer, lhist)
				log.add("=========================")
				for lhx in xrange(lhist):
					msg = " %10.3f     %7d"%(region[lhx], histo[lhx])
					log.add(msg)
				if (max(all_pixer) < 0.5) and (sum(all_pixer)/total_nima < 0.05):
					terminate = 1
			terminate = wrap_mpi_bcast(terminate, main_node, mpi_comm)
			#=========================================================================

			#=========================================================================
			# centering
			if center == -1 and sym[0] == 'c':
				from utilities      import estimate_3D_center_MPI, rotate_3D_shift
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

			#=========================================================================
			# volume reconstruction
			mpi_barrier(mpi_comm)
			if myid == main_node:
				start_time = time()
			vol = volume_reconstruction(data, ali3d_options, mpi_comm)
			# log
			if myid == main_node:
				log.add("3D reconstruction time = %f\n"%(time()-start_time))
				start_time = time()
			#=========================================================================

	#=========================================================================
	# gather parameters
	params = []
	previousmax = []
	for im in data:
		t = get_params_proj(im)
		p = im.get_attr("previousmax")
		params.append( [t[0], t[1], t[2], t[3], t[4]] )
		previousmax.append(p)
	assert(nima == len(params))
	params = wrap_mpi_gatherv(params, 0, mpi_comm)
	if myid == 0:
		assert(total_nima == len(params))
	previousmax = wrap_mpi_gatherv(previousmax, 0, mpi_comm)

	par_r = wrap_mpi_gatherv(par_r, 0, mpi_comm)

	if myid == main_node: 
		log.add("Finish ali3d_multishc_2")
		return params, vol, previousmax, par_r
	else:
		return None, None, None, None  # results for the other processes
"""

# data - projections (scattered between cpus)
# options - the same for all cpus
# return - volume the same for all cpus
def volume_reconstruction(data, options, mpi_comm):
	from mpi import mpi_comm_rank
	from reconstruction import recons3d_4nn_MPI, recons3d_4nn_ctf_MPI
	from utilities import bcast_EMData_to_all, model_circle
	
	myid = mpi_comm_rank(mpi_comm)
	sym  = options.sym
	sym = sym[0].lower() + sym[1:]
	npad      = options.npad
	user_func = options.user_func
	CTF       = options.CTF
	snr       = options.snr
	center    = options.center
	#=========================================================================
	# volume reconstruction
	if CTF: vol = recons3d_4nn_ctf_MPI(myid, data, snr, symmetry=sym, npad=npad, mpi_comm=mpi_comm)
	else:   vol = recons3d_4nn_MPI    (myid, data,      symmetry=sym, npad=npad, mpi_comm=mpi_comm)

	if myid == 0:
		nx = data[0].get_xsize()
		last_ring   = int(options.ou)
		mask3D = model_circle(last_ring, nx, nx, nx)
		ref_data = [ mask3D, max(center,0), None, None, None, None ]
		ref_data[2] = vol
		ref_data[3] = None #fscc
		ref_data[4] = None#varf
		#  call user-supplied function to prepare reference image, i.e., center and filter it
		vol, cs = user_func(ref_data)

	# broadcast volume
	bcast_EMData_to_all(vol, myid, 0, comm=mpi_comm)
	#=========================================================================
	return vol


# all_projs and subset must be set only for root (MPI rank == 0)
# remaining parameters must be set for all
# size of mpi_communicator must be >= runs_count
def multi_shc(all_projs, subset, runs_count, ali3d_options, mpi_comm, log=None, ref_vol=None):
	from applications import MPI_start_end
	from mpi import mpi_comm_rank, mpi_comm_size
	from utilities import set_params_proj, wrap_mpi_bcast, write_text_row, drop_image, write_text_file
	from random import random

	mpi_rank = mpi_comm_rank(mpi_comm)
	mpi_size = mpi_comm_size(mpi_comm)

	assert (mpi_size >= runs_count)

	if log == None:
		from logger import Logger
		log = Logger()
	
	projections = []
	if mpi_rank == 0:
		all_projs_params = generate_uneven_projections_directions(len(all_projs), half_sphere=False)
		for i in subset:
			all_projs_params[i][2] = random()*360.0
			set_params_proj(all_projs[i], all_projs_params[i])
			all_projs[i].set_attr("stable", 0)
			projections.append(all_projs[i])
			j = 1
			while all_projs[i].has_attr("xform.projection" + str(j)):
				all_projs[i].del_attr("xform.projection" + str(j))
				j += 1
	projections = wrap_mpi_bcast(projections, 0, mpi_comm)

	n_projs = len(projections)

	if ref_vol == None:
		proj_begin, proj_end = MPI_start_end(n_projs, mpi_size, mpi_rank)
		ref_vol = volume_reconstruction(projections[proj_begin:proj_end], ali3d_options, mpi_comm=mpi_comm)

	out_params, out_vol, out_peaks = ali3d_multishc(projections, ref_vol, ali3d_options, mpi_comm=mpi_comm, log=log, number_of_runs=runs_count)
	if mpi_rank == 0:
		assert(len(out_params) == runs_count)

	if mpi_rank == 0:
		write_text_file(subset, log.prefix + "indexes.txt")
		for i in xrange(len(out_params)):
			write_text_row(out_params[i], log.prefix + "part_" + str(i) + "_params.txt")
			drop_image(out_vol[i], log.prefix + "part_" + str(i) + "_volf.hdf")
			#write_text_row(out_peaks[i], log.prefix + "part_" + str(i) + "_peaks.txt")

		temp_projs = []
		for iP in xrange(len(out_params[0])):
			iBestPeak = 0
			for iC in xrange(len(out_params)):
				if out_peaks[iC][iP] > out_peaks[iBestPeak][iP]:  iBestPeak = iC
				temp_projs.append( projections[iP].copy() )
				set_params_proj( temp_projs[len(temp_projs)-1], out_params[iC][iP])
			set_params_proj( projections[iP], out_params[iBestPeak][iP] )
			projections[iP].set_attr("stable", 1)
	else:
		temp_projs = None

	temp_projs = wrap_mpi_bcast(temp_projs, 0, mpi_comm)
	proj_begin, proj_end  = MPI_start_end(3*n_projs, mpi_size, mpi_rank)
	ref_vol = volume_reconstruction(temp_projs[proj_begin:proj_end], ali3d_options, mpi_comm=mpi_comm)

	out_params, out_vol, out_peaks, out_r = ali3d_multishc_2(projections, ref_vol, ali3d_options, mpi_comm=mpi_comm, log=log)
	if mpi_rank == 0:
		assert(len(out_params) == n_projs)

	if mpi_rank == 0:
		write_text_row(out_params, log.prefix + "params.txt")
		drop_image(out_vol, log.prefix + "volf.hdf")

	return out_params, out_vol, out_peaks


def reduce_dsym_angles(p1, sym):
	#  works only for d symmetry
	from utilities import get_symt
	from EMAN2 import Vec2f, Transform
	t = get_symt(sym)
	ns = int(sym[1:])
	for i in xrange(len(t)):  t[i] = t[i].inverse()

	for i in xrange(len(p1)):
		 a = Transform({"type":"spider","phi":p1[i][0], "theta":p1[i][1], "psi":p1[i][2]})
		 a.set_trans(Vec2f(-p1[i][3], -p1[i][4]))
		 for l in xrange(len(t)):
			q = a*t[l]
			q = q.get_params("spider")
			if(q["phi"]<360./ns and q["theta"] <= 90.0): break
		 p1[i][0] = q["phi"]
		 p1[i][1] = q["theta"]
		 p1[i][2] = q["psi"]
		 p1[i][3] = -q["tx"]
		 p1[i][4] = -q["ty"]
