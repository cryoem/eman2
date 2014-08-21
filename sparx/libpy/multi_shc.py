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
	T = T1*v2
	return [ T.get_params("spider")["phi"], T.get_params("spider")["theta"], T.get_params("spider")["psi"], T.get_params("spider")["tx"], T.get_params("spider")["ty"]  ]


def orient_params(params, indexes=None, sym = "c1"):
	from utilities import rotation_between_anglesets
	from pixel_error import angle_diff

	m = len(params)
	n = len(params[0])

	if(sym[0] == "d"):
		# In this one the first params is taken as a reference
		mirror_and_reduce_dsym(params, sym)
		"""
		if indexes == None:
			mirror_and_reduce_dsym(params, sym)
		else:
			temp = []
			for i in xrange(m):
				temp.append([params[i][j] for j in indexes])
			mirror_and_reduce_dsym(temp, sym)
			for i in xrange(1,m):
				k = 0
				for j in indexes:
					params[i][j] = temp[i][k]
					k += 1
		"""
	else:
		from EMAN2 import Transform
		if indexes == None:   indexes = range(n)
		for i in xrange(1,m):
			cmp_par_i = []
			cmp_par_0 = []
			for j in indexes:
				cmp_par_i.append(params[i][j])
				cmp_par_0.append(params[0][j])
			t1,t2,t3 = rotation_between_anglesets(cmp_par_i, cmp_par_0)
			rot = Transform({"type":"spider","phi":t1,"theta":t2,"psi":t3})
			for j in xrange(n):
				params[i][j] = mult_transform(params[i][j], rot)
			# mirror checking
			psi_diff = angle_diff( [params[i][j][2] for j in indexes], [params[0][j][2] for j in indexes] )
			if(abs(psi_diff-180.0) <90.0):
				#mirror
				for j in indexes:
					params[i][j][2] = (params[i][j][2] + 180.0) % 360.0
"""

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
"""

def calculate_matrix_rot(projs):
	from utilities import rotation_between_anglesets
	sc = len(projs)
	matrix_rot  = [[[0.0,0.0,0.0,0.0,0.0] for i in xrange(sc)] for k in xrange(sc)]
	for i in xrange(sc-1):
		for j in xrange(i+1,sc):
			t1, t2, t3 = rotation_between_anglesets(projs[i], projs[j])
			matrix_rot[i][j] = [t1, t2, t3, 0.0, 0.0]
	return matrix_rot


# returns subset_for_threshold, subset_for_minimal_size[, threshold_for_thr_subset, threshold_for_min_subset]
def find_common_subset_3(projs, target_threshold, minimal_subset_size=3, sym = "c1", thresholds=False):
	from global_def import Util

	n  = len(projs[0])
	sc = len(projs)

	subset = range(n)

	minimal_subset_size = min( minimal_subset_size, n)

	res_thr_subset  = None
	res_size_subset = None
	error_thr_subset = -1
	error_size_subset = -1

	for iIter in xrange(n-2):
		projs2 = [0.0]*sc
		for iConf in xrange(sc):
			projs2[iConf] = []
			for i in subset:
				projs2[iConf].append(projs[iConf][i][:])
		if( sym[0] == "d"):
			# have to figure whether anything has to be mirrored and then reduce the angles.
			mirror_and_reduce_dsym(projs2, sym)
			from utilities import getvec
			from math import acos, degrees

			trans_vec = [0.0]*sc
			for iConf in xrange(sc):
				#temp = reduce_dsym_angles(projs2[iConf], sym)
				trans_vec[iConf] = []
				for i in xrange(len(projs2[0])):
					t1,t2,t3=getvec(projs2[iConf][i][0], projs2[iConf][i][1])
					trans_vec[iConf].append([t1,t2,t3])

			avg_diff_per_image = []
			for i in xrange(len(trans_vec[0])):
				qt = 0.0
				for k in xrange(sc-1):
					for l in xrange(k+1,sc):
						zt = 0.0
						for m in xrange(3):  zt += trans_vec[k][i][m]*trans_vec[l][i][m]
						qt += degrees(acos(min(1.0,max(-1.0,zt))))
				avg_diff_per_image.append(qt/sc/(sc-1)/2.0)

		else:
			trans_projs = []
			matrix_rot = calculate_matrix_rot(projs2)
			for iConf in xrange(sc):
				for i in subset:
					trans_projs.extend(projs[iConf][i][0:5])
			trans_matrix = []
			for i in xrange(sc):
				for j in xrange(i):
					trans_matrix.extend(matrix_rot[i][j][0:3])
			avg_diff_per_image = Util.diff_between_matrix_of_3D_parameters_angles(trans_projs, trans_matrix)
		#print  "  AAAA ",iIter
		#print avg_diff_per_image
		max_error = -1.0
		the_worst_proj = -1
		for i in xrange(len(avg_diff_per_image)):
			if avg_diff_per_image[i] > max_error:
				max_error = avg_diff_per_image[i]
				the_worst_proj = subset[i]
		if max_error <= target_threshold:
			res_thr_subset   = subset[:]
			error_thr_subset = max_error
			break
		if len(subset) == minimal_subset_size:
			res_size_subset   = subset[:]
			error_size_subset = max_error
		subset.remove(the_worst_proj)

	if res_thr_subset == None:
		res_thr_subset   = subset
		error_thr_subset = max_error

	if res_size_subset == None:
		res_size_subset   = res_thr_subset
		error_size_subset = error_thr_subset

	if thresholds:
		return res_thr_subset, res_size_subset, error_thr_subset, error_size_subset
	return res_thr_subset, res_size_subset, avg_diff_per_image


# parameters: list of (all) projections | reference volume | ...
#  Genetic programming version
#  The data structure:
#  [[L2, [parameters row-wise]], [], []...number_of_runs ]
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
	doga   = ali3d_options.doga
	center = ali3d_options.center
	CTF    = ali3d_options.CTF
	ref_a  = ali3d_options.ref_a
	L2threshold = ali3d_options.L2threshold

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	if log == None:
		from logger import Logger
		log = Logger()

	number_of_proc = mpi_comm_size(mpi_comm)
	myid           = mpi_comm_rank(mpi_comm)
	main_node = 0

	if myid == main_node:
		log.add("Start VIPER1")

	if number_of_proc < number_of_runs:
		ERROR("number_of_proc < number_of_runs","VIPER1")
	
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
	max_iter    = int(ali3d_options.maxit1)
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
	#print "  image_start, image_end  ", myid,image_start, image_end

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
	if myid == main_node:  noimprovement = 0
	# do the projection matching
	for N_step in xrange(lstp):  # At this point there is just one value here, it cannot loop.
		if myid == 0:  afterGAcounter = 0
		terminate = False
		Iter = 0
		while not terminate:

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
			#  It has to be here
			if orient_and_shuffle:
				# adjust params to references, calculate psi, calculate previousmax
				for im in xrange(nima):
					peak, temp = proj_ali_incore_local(data[im],refrings,numr,0.,0.,1., delta[N_step]*0.7 )
					data[im].set_attr("previousmax", peak)
				if myid == main_node:
					log.add("Time to calculate first psi+previousmax: %f\n" % (time()-start_time))
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

			storevol=False

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
				orient_and_shuffle = ( percent_of_pixerr_below_one > doga )  and ( afterGAcounter < 0 ) #  TODO - parameter ?
				afterGAcounter -= 1
				## if total_iter%3 == 0:  orient_and_shuffle = True
				## else:   orient_and_shuffle = False
				# terminate          = ( percent_of_pixerr_below_one > 0.9 )  #  TODO - parameter ?
				log.add("=========================")
				log.add("Percent of positions with pixel error below 1.0 = ", (int(percent_of_pixerr_below_one*100)), "%","   Mutations: ",orient_and_shuffle)
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
			#===================== CORRECT PARAMETERS ON DATA =======================

			mpi_barrier(mpi_comm)
			if myid == main_node:
				start_time = time()
			"""
			params = []
			for im in data:
				phi, theta, psi, sx, sy = get_params_proj(im)
				params.append([phi, theta, psi, sx, sy])
			params_0 = wrap_mpi_bcast(params, mpi_subroots[0], mpi_comm)

			if mpi_subrank == 0:
				from utilities import write_text_row
				write_text_row(params, "qparams%04d%04d.hdf"%(myid,total_iter) )
			"""

			#=========================================================================
			if orient_and_shuffle:   #  DO orient
				params = []
				for im in data:
					phi, theta, psi, sx, sy = get_params_proj(im)
					params.append([phi, theta, psi, sx, sy])
				# if myid == 2:  print  " initial params before orient  ",myid,[get_params_proj(data[i]) for i in xrange(4)]

				# ------ orientation - begin
				params_0 = wrap_mpi_bcast(params, mpi_subroots[0], mpi_comm)
				if mpi_subrank == 0:

					subset_thr, subset_min, avg_diff_per_image = find_common_subset_3([params_0, params], 2.0, len(params)/3, sym)
					if len(subset_thr) < len(subset_min):
						subset = subset_min
					else:
						subset = subset_thr
					# if myid == 2:  print  " params before orient  ",myid,params[:4],params[-4:]
					from utilities import write_text_row
					#write_text_row(params_0,"bparamszero%04d%04d.txt"%(myid,total_iter))
					#write_text_row(params,"bparams%04d%04d.txt"%(myid,total_iter))
					orient_params([params_0, params], subset, sym)
					"""
					if myid == 2:
						print  " subset  ",len(subset)#," ...  ",subset
					"""
					from utilities import write_text_row
					#write_text_row(params_0,"aparamszero%04d%04d.txt"%(myid,total_iter))
					#write_text_row(params,"aparams%04d%04d.txt"%(myid,total_iter))
					# if myid == 2:  print  " params after orient  ",myid,params[:4],params[-4:]
				params = wrap_mpi_bcast(params, 0, mpi_subcomm)
				# if myid == 2:  print  " params after wrap_mpi_bcast  ",myid,params[:4],params[-4:]
				# ------ orientation - end


				#=========================================================================
				# volume reconstruction
				mpi_barrier(mpi_comm)
				if myid == main_node:
					start_time = time()

				#temp = [None]*nima
				#for i in xrange(nima): temp[i] = data[i].get_attr("xform.projection")
				for i in xrange(nima):  set_params_proj(data[i], params[i])
				vol = volume_reconstruction(data[image_start:image_end], ali3d_options, mpi_subcomm)
				#for i in xrange(nima): data[i].set_attr("xform.projection",temp[i])
				#del temp


				if mpi_subrank == 0:
					L2 = vol.cmp("dot", vol, dict(negative = 0, mask = model_circle(last_ring, nx, nx, nx)))
					# if myid == 2:  print  " Right after reconstruction L2", myid, L2,[get_params_proj(data[i]) for i in xrange(4)]
					#print  " Right after reconstruction of oriented parameters L2", myid, total_iter,L2
					#vol.write_image("recvolf%04d%04d.hdf"%(myid,total_iter))
				# log
				if myid == main_node:
					log.add("3D reconstruction time = %f\n"%(time()-start_time))
					start_time = time()

				# ------ gather parameters to root
				if myid == 0:
					all_L2s = []
					all_params = []
					for sr in mpi_subroots:
						if sr == myid:
							all_L2s.append(L2)
							all_params.append(params)
						else:
							all_L2s.append(wrap_mpi_recv(sr, mpi_comm))
							all_params.append(wrap_mpi_recv(sr, mpi_comm))
				else:
					if mpi_subrank == 0:
						wrap_mpi_send(L2, 0, mpi_comm)
						wrap_mpi_send(params, 0, mpi_comm)

				# ---------------------------------

				#  Add params to GA, sort, check termination and if not terminate do mutations and send back
				if myid == 0:
					#  after GA move do 3 iterations to give the program a chance to improve mutated structures.
					#all_params = shuffle_configurations(all_params)
					for i in xrange(number_of_runs):
						GA.append([all_L2s[i],all_params[i]])
					#  check whether this move will improve anything
					all_L2s.sort(reverse=True)
					#print " sorted terminate  ",all_L2s
					#for i in xrange(number_of_runs): print GA[i][0]
					if(all_L2s[0]<GA[number_of_runs-1][0]):
						noimprovement += 1
						if(noimprovement == 2):  terminate = True
						GA = GA[:number_of_runs]
					else:
						noimprovement = 0
						GA.sort(reverse=True)
						GA = GA[:number_of_runs]
						if( sym[0] == "d"  and   GA[0][0]>0.0 ):
							for i in xrange(1,len(GA)):
								mirror_and_reduce_dsym([GA[0][1],GA[i][1]], sym)

						#  ---  Stopping criterion
						from statistics import table_stat
						from math import sqrt
						q1,q2,q3,q4 = table_stat([GA[i][0] for i in xrange(number_of_runs)])
						# Terminate if variation of L2 norms less than (L2threshold*100)% of their average
						crit = sqrt(max(q2,0.0))/q1
						terminate = Iter > max_iter or crit < L2threshold
						##  if  total_iter > 17: terminate = True
						##  else:  terminate=False
						for i in xrange(number_of_runs):
							log.add("L2 norm for volume %3d  = %f"%(i,GA[i][0]))
						log.add("L2 norm std dev %f\n"%crit)

					if not terminate and noimprovement == 0:
						afterGAcounter = 3
						#  Now do the mutations
						all_params = []

						from utilities import nearestk_projangles
						from random import random, randint, shuffle
						# select random pairs of solutions
						ipl = range(number_of_runs)
						shuffle(ipl)
						for ip in xrange(0,2*(len(ipl)/2)+len(ipl)%2,2):
							#  random reference projection:
							itmp = randint(0,total_nima-1)
							#  if(   )
							keepset = nearestk_projangles(GA[ipl[ip]][1], whichone = itmp, howmany = total_nima/2, sym=sym)
							keepset.append(itmp)
							otherset = set(range(total_nima)) - set(keepset)
							otherset = [i for i in otherset]
							keepset.sort()
							otherset.sort()
							newparms1 = [None]*total_nima
							newparms2 = [None]*total_nima
							for i in keepset:
								newparms1[i] = GA[ipl[ip]][1][i]
								newparms2[i] = GA[ipl[(ip+1)%number_of_runs]][1][i]
							for i in otherset:
								newparms1[i] = GA[ipl[(ip+1)%number_of_runs]][1][i]
								newparms2[i] = GA[ipl[ip]][1][i]
							#print "  PRINTOUT SHUFFLED   ",ipl[ip],ipl[ip+1]
							"""
							for i in xrange(total_nima):
								print  i,newparms1[i],GA[ipl[ip]][1][i]
							for i in xrange(total_nima):
								print  i,newparms2[i],GA[ipl[ip+1]][1][i]
							for i in xrange(total_nima):
								GA[ipl[ip]][1][i]   = newparms1[i]
								GA[ipl[ip+1]][1][i] = newparms2[i]
							"""

							#  Put mutated params on one list
							all_params.append(newparms1)
							all_params.append(newparms2)
						all_params = all_params[:number_of_runs]

				terminate = wrap_mpi_bcast(terminate, main_node, mpi_comm)
				if not terminate:

					storevol=True

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
					"""
					#=========================================================================
					# volume reconstruction
					mpi_barrier(mpi_comm)
					if myid == main_node:
						start_time = time()
					vol = volume_reconstruction(data[image_start:image_end], ali3d_options, mpi_subcomm)
					##if mpi_subrank == 0:  vol.write_image("mutatedvolf%04d%04d.hdf"%(myid,total_iter))
					if mpi_subrank == 0:
						L2 = vol.cmp("dot", vol, dict(negative = 0, mask = model_circle(last_ring, nx, nx, nx)))
						# if myid == 2:  print  " Right after reconstruction L2", myid, L2,[get_params_proj(data[i]) for i in xrange(4)]
						print  "Mutated L2", myid, total_iter,L2

					##if myid == main_node:
					##	from utilities import write_text_row
					##	write_text_row(GA[0][1], "qparams%04d%04d.txt"%(myid,total_iter) )
					"""
					#=========================================================================
					#
					mpi_barrier(mpi_comm)
					if myid == main_node:
						log.add("Time of orientation and mutations = %f\n"%(time()-start_time))
						start_time = time()
				#else:  continue
			if not terminate:
				#=========================================================================
				# volume reconstruction
				mpi_barrier(mpi_comm)
				if myid == main_node:
					start_time = time()
				vol = volume_reconstruction(data[image_start:image_end], ali3d_options, mpi_subcomm)

				if mpi_subrank == 0:
					L2 = vol.cmp("dot", vol, dict(negative = 0, mask = model_circle(last_ring, nx, nx, nx)))
					# if myid == 2:  print  " Right after reconstruction L2", myid, L2,[get_params_proj(data[i]) for i in xrange(4)]
					#print  " Right after reconstruction L2", myid, total_iter,L2
					#if storevol:   vol.write_image("mutated%04d%04d.hdf"%(myid,total_iter))

				# log
				if myid == main_node:
					log.add("3D reconstruction time = %f\n"%(time()-start_time))
					start_time = time()

			"""
			#VERIFY GA
			if myid == 0:
				temp = [None]*nima
				for i in xrange(nima): temp[i] = data[i].get_attr("xform.projection")
				for k in xrange(len(GA)):
					for i in xrange(nima):
						set_params_proj(data[i], GA[k][1][i])
					tvol = volume_reconstruction(data, ali3d_options, mpi_subcomm)
					LL2 = tvol.cmp("dot", tvol, dict(negative = 0, mask = model_circle(last_ring, nx, nx, nx)))
					print  "GA VERIFY  ",k,GA[k][0],LL2,GA[k][1][:4],GA[k][1][-5:]
				for i in xrange(nima): data[i].set_attr("xform.projection",temp[i])
			"""
			"""
			# Send params back
			if myid == 0:
				#print  all_params
				#print "GA  ",GA
				for i in xrange(number_of_runs):
					sr = mpi_subroots[i]
					if sr == myid:
						params = GA[i][1]
					else:
						wrap_mpi_send(GA[i][1], sr, mpi_comm)
			else:
				if mpi_subrank == 0:
					params = wrap_mpi_recv(0, mpi_comm)

			params = wrap_mpi_bcast(params, 0, mpi_subcomm)
			if myid == 0:
				print  "params  "
				print params
			for i in xrange(nima):
				set_params_proj(data[i], params[i])

			#=========================================================================
			# volume reconstruction
			mpi_barrier(mpi_comm)
			if myid == main_node:
				start_time = time()
			vol = volume_reconstruction(data[image_start:image_end], ali3d_options, mpi_subcomm)
			if mpi_subrank == 0:
				L2 = vol.cmp("dot", vol, dict(negative = 0, mask = model_circle(last_ring, nx, nx, nx)))
				print  "VERIFICATION  ",myid,L2
				#vol.write_image("qmutatedvolf%04d%04d.hdf"%(myid,total_iter))

			#if myid == main_node:
			#	from utilities import write_text_row
			#	write_text_row(GA[0][1], "qparams%04d%04d.txt"%(myid,total_iter) )

			#=========================================================================
			#
			mpi_barrier(mpi_comm)
			if myid == main_node:
				log.add("Time of verification = %f\n"%(time()-start_time))
				start_time = time()
			"""



	#=========================================================================
	mpi_comm_free(mpi_subcomm)
	
	
	if myid == main_node:
		log.add("Finish viper1")
		return GA[0][1]
	else:
		return None  # results for the other processes


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
		log.add("Start VIPER2")

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
	max_iter    = int(ali3d_options.maxit2)
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




	#=========================================================================
	# adjust params to references, calculate psi+shifts, calculate previousmax
	qvol = volume_reconstruction(data, ali3d_options, mpi_comm)
	# log
	"""
	if myid == main_node:
		start_time = time()
		#qvol.write_image("vitera.hdf")
		L2 = qvol.cmp("dot", qvol, dict(negative = 0, mask = model_circle(last_ring, nx, nx, nx)))
		log.add("3D reconstruction time = %f\n"%(time()-start_time)," START  L2 norm:  %f"%L2)
		start_time = time()
	del qvol
	"""

	from projection   import prep_vol, prgs
	from alignment import ringwe
	cnx = nx//2 + 1
	cny = nx//2 + 1
	wr_four  = ringwe(numr, "F")
	from math import pi, sin, cos
	qv = pi/180.
	volft, kb = prep_vol(ref_vol)
	from utilities import get_params_proj
	for im in xrange(nima):
		phi,theta,psi,tx,ty = get_params_proj(data[im])
		prjref = prgs(volft, kb, [phi,theta,psi, 0.0, 0.0])
		cimage = Util.Polar2Dm(prjref, cnx, cny, numr, "F")
		Util.Normalize_ring(cimage, numr)
		Util.Frngs(cimage, numr)
		Util.Applyws(cimage, numr, wr_four)
		refrings = [cimage]
		n1 = sin(theta*qv)*cos(phi*qv)
		n2 = sin(theta*qv)*sin(phi*qv)
		n3 = cos(theta*qv)
		refrings[0].set_attr_dict( {"n1":n1, "n2":n2, "n3":n3} )
		refrings[0].set_attr("phi",   phi)
		refrings[0].set_attr("theta", theta)
		refrings[0].set_attr("psi",   psi)

		#print "orin  ",data[im].get_attr("ID"),get_params_proj(data[im])
		peak, pixer = proj_ali_incore_local(data[im],refrings,numr,0.0,0.0,1.0,delta[0]/4)
		data[im].set_attr("previousmax", peak)
		#print  "peak ",data[im].get_attr("ID"), peak,pixer,get_params_proj(data[im])
	del volft
	# volume reconstruction
	mpi_barrier(mpi_comm)
	if myid == main_node:
		start_time = time()
	ref_vol = volume_reconstruction(data, ali3d_options, mpi_comm)
	# log
	if myid == main_node:
		##ref_vol.write_image("viterb.hdf")
		L2 = ref_vol.cmp("dot", ref_vol, dict(negative = 0, mask = model_circle(last_ring, nx, nx, nx)))
		log.add("3D reconstruction time = %f\n"%(time()-start_time),"   L2 norm:  %f"%L2)
		start_time = time()

	if myid == main_node:
		log.add("Time to calculate first psi+shifts+previousmax: %f\n" % (time()-start_time))
		start_time = time()

	#=========================================================================



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

			mpi_barrier(mpi_comm)
			if myid == main_node:
				start_time = time()
			#=========================================================================
			# alignment
			number_of_checked_refs = 0
			for im in xrange(nima):
				#peak, pixer[im], checked_refs, number_of_peaks = shc_multi(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step], number_of_runs=number_of_runs)
				# previousmax is set in shc
				peak, pixer[im], checked_refs, iref = shc(data[im], refrings, numr, xrng[N_step], yrng[N_step], step[N_step], an[N_step])
				number_of_checked_refs += checked_refs
			#=========================================================================
			mpi_barrier(mpi_comm)
			if myid == main_node:
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
				#vol.write_image("viter%03d.hdf"%total_iter)
				L2 = vol.cmp("dot", vol, dict(negative = 0, mask = model_circle(last_ring, nx, nx, nx)))
				log.add("3D reconstruction time = %f\n"%(time()-start_time),"   L2 norm:  %f"%L2)
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
	if myid == 0:
		assert(total_nima == len(previousmax))

	par_r = wrap_mpi_gatherv(par_r, 0, mpi_comm)

	if myid == main_node: 
		log.add("Finish VIPER2")
		return params, vol, previousmax, par_r
	else:
		return None, None, None, None  # results for the other processes

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


def volume_recsp(data, options):

	from reconstruction import recons3d_4nn, recons3d_4nn_ctf
	from utilities import bcast_EMData_to_all, model_circle
	
	sym  = options.sym
	sym = sym[0].lower() + sym[1:]
	npad      = options.npad
	user_func = options.user_func
	CTF       = options.CTF
	snr       = options.snr
	center    = options.center
	#=========================================================================
	# volume reconstruction
	if CTF: vol = recons3d_4nn_ctf(data, snr, symmetry=sym, npad=npad)
	else:   vol = recons3d_4nn(data,      symmetry=sym, npad=npad)

	nx = data[0].get_xsize()
	last_ring   = int(options.ou)
	mask3D = model_circle(last_ring, nx, nx, nx)
	ref_data = [ mask3D, max(center,0), None, None, None, None ]
	ref_data[2] = vol
	ref_data[3] = None #fscc
	ref_data[4] = None#varf
	#  call user-supplied function to prepare reference image, i.e., center and filter it
	vol, cs = user_func(ref_data)

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
			#all_projs[i].set_attr("stable", 0)
			all_projs[i].set_attr("previousmax", -1.e23)
			projections.append(all_projs[i])
			"""
			j = 1
			while all_projs[i].has_attr("xform.projection" + str(j)):
				all_projs[i].del_attr("xform.projection" + str(j))
				j += 1
			"""
	projections = wrap_mpi_bcast(projections, 0, mpi_comm)

	n_projs = len(projections)

	if ref_vol == None:
		proj_begin, proj_end = MPI_start_end(n_projs, mpi_size, mpi_rank)
		ref_vol = volume_reconstruction(projections[proj_begin:proj_end], ali3d_options, mpi_comm=mpi_comm)


	# Each node keeps all projection data, this would not work for large datasets
	out_params = ali3d_multishc(projections, ref_vol, ali3d_options, mpi_comm=mpi_comm, log=log, number_of_runs=runs_count)
	"""
	if mpi_rank == 0:
		assert(len(out_params) == runs_count)
	"""

	if mpi_rank == 0:
		"""
		write_text_file(subset, log.prefix + "indexes.txt")
		for i in xrange(len(out_params)):
			write_text_row(out_params[i], log.prefix + "run_" + str(i) + "_params.txt")
			drop_image(out_vol[i], log.prefix + "run_" + str(i) + "_volf.hdf")
			#write_text_row(out_peaks[i], log.prefix + "run_" + str(i) + "_peaks.txt")

		temp_projs = []
		for iP in xrange(len(out_params[0])):
			iBestPeak = 0
			for iC in xrange(len(out_params)):
				if out_peaks[iC][iP] > out_peaks[iBestPeak][iP]:  iBestPeak = iC
				temp_projs.append( projections[iP].copy() )
				set_params_proj( temp_projs[len(temp_projs)-1], out_params[iC][iP])
			set_params_proj( projections[iP], out_params[iBestPeak][iP] )
			projections[iP].set_attr("stable", 1)
		#  Use the best one to finish off
		for iP in xrange(len(out_params[0])):
			projections[iP].set_attr("stable", 1)
			set_params_proj( projections[iP], out_params[0][iP] )
		"""
		temp = []
		from utilities import get_params_proj
		for i in xrange(n_projs):
			set_params_proj( projections[i], out_params[i] )
		write_text_row(out_params, log.prefix + "refparams2.txt")
		"""
		log.add("  WILL RECONSTRUCT  ")
		from utilities import model_circle
		tvol = volume_recsp(projections, ali3d_options)
		LL2 = tvol.cmp("dot", tvol, dict(negative = 0, mask = model_circle(22, 48, 48, 48)))
		log.add(" LLLLLLL2 norm of reference volume:  %f"%LL2)
		for k in xrange(mpi_size):
			proj_begin, proj_end  = MPI_start_end(n_projs, mpi_size, k)
			print  " from within   ",k,proj_begin,get_params_proj(projections[proj_begin])
		"""
	"""
	else:
		temp_projs = None
	"""
	proj_begin, proj_end  = MPI_start_end(n_projs, mpi_size, mpi_rank)

	projections = wrap_mpi_bcast(projections, 0, mpi_comm)
	from utilities import get_params_proj
	#print  " from mpi   ",mpi_rank,proj_begin,get_params_proj(projections[proj_begin])
	ref_vol = volume_reconstruction(projections[proj_begin:proj_end], ali3d_options, mpi_comm=mpi_comm)
	if mpi_rank == 0:
		ref_vol.write_image(log.prefix + "refvol2.hdf")
		from utilities import model_circle
		nx = ref_vol.get_xsize()
		L2 = ref_vol.cmp("dot", ref_vol, dict(negative = 0, mask = model_circle(ali3d_options.ou, nx,nx,nx)))
		log.add(" L2 norm of reference volume:  %f"%L2)

	"""
	if mpi_rank == 17:
		temp = []
		from utilities import get_params_proj
		for i in xrange(n_projs):
			#projections[i].set_attr("stable", 1)
			t1,t2,t3,t4,t5 = get_params_proj( projections[i])
			t6 = projections[i].get_attr("previousmax")
			temp.append([t1,t2,t3,t4,t5,t6])
			#set_params_proj( projections[i], out_params[i] )
		write_text_row(temp, log.prefix + "refparams17.txt")
	"""

	out_params, out_vol, previousmax, out_r = ali3d_multishc_2(projections, ref_vol, ali3d_options, mpi_comm=mpi_comm, log=log)
	if mpi_rank == 0:
		write_text_file(out_params, log.prefix + "previousmax.txt")
		write_text_row(out_params, log.prefix + "params.txt")
		drop_image(out_vol, log.prefix + "volf.hdf")

	return out_params, out_vol, None#, out_peaks

def reduce_dsym_angles(p1, sym):
	#  WARNING - it returns incorrect parameters that are only suitable for calculation of angular distances
	#  works only for d symmetry
	pr = [[0.0 for i in xrange(5)] for q in xrange(len(p1))]
	for i in xrange(len(p1)):
		if( p1[i][1] >90.0):
			p1[i][1] = 180.0 - p1[i][1]
			p1[i][0] = (p1[i][0] +180.0)%360.0
	"""
	from utilities import get_symt
	from EMAN2 import Vec2f, Transform
	t = get_symt(sym)
	phir = 360.0/int(sym[1:])
	for i in xrange(len(t)):  t[i] = t[i].inverse()
	pr = [[0.0 for i in xrange(5)] for q in xrange(len(p1))]
	for i in xrange(len(p1)):
		 a = Transform({"type":"spider","phi":p1[i][0], "theta":p1[i][1], "psi":p1[i][2]})
		 a.set_trans(Vec2f(-p1[i][3], -p1[i][4]))
		 for l in xrange(len(t)):
			q = a*t[l]
			q = q.get_params("spider")
			if(q["phi"]<phir and q["theta"] <= 90.0): break
		 pr[i][0] = q["phi"]
		 pr[i][1] = q["theta"]
		 pr[i][2] = q["psi"]
		 pr[i][3] = -q["tx"]
		 pr[i][4] = -q["ty"]
	"""
	return pr

def mirror_and_reduce_dsym(params, sym):
	# for the time being do it in a silly way, i.e., select first as a reference and run with it.
	#
	
	from utilities import get_symt, get_sym, getfvec
	from EMAN2 import Vec2f, Transform, EMData
	from pixel_error import angle_diff

	def discangset(pari, para, sym = "d3"):
		from pixel_error import max_3D_pixel_error
		ts = get_symt(sym)
		ks = len(ts)
		for i in xrange(ks):  ts[i] = ts[i].inverse()
		per1 = 0.0
		for j in xrange(len(pari)):
			apixer = 1.e20
			qt = Transform({"type":"spider","phi":pari[j][0], "theta":pari[j][1], "psi":pari[j][2]})
			for k in xrange(ks):
				ut = qt*ts[k]
				tmp = max_3D_pixel_error(ut, para[j])
				if(tmp < apixer): apixer = tmp
			per1 += apixer
		return per1



	sc = len(params)
	ns = len(params[0])
	#  Convert 0 to transforms
	#t0 = [None]*ns
	vt0 = [None]*ns
	for j in xrange(ns):
		#vt0[j] = getfvec(params[0][j][0], params[0][j][1])
		vt0[j] = Transform({"type":"spider","phi":params[0][j][0], "theta":params[0][j][1], "psi":params[0][j][2]})
		#t0[j].set_trans(Vec2f(-params[0][j][3], -params[0][j][4]))
	#  get sym transforms
	ts = get_symt(sym)
	ks = len(ts)
	for i in xrange(ks):  ts[i] = ts[i].inverse()

	# set rotation for mirror
	mm = Transform({"type":"spider","phi":360./2/int(sym[1:]), "theta":0., "psi":0.})
	symphi = 360./int(sym[1:])

	#  Excluded values for phi
	badb = 0.0
	bade = 360.0/int(sym[1:])/4
	bbdb = 360.0/int(sym[1:])/2
	bbde = bbdb + 360.0/int(sym[1:])/4
	
	for i in xrange(1,sc):
		solvs = []
		for rphi in xrange(0,61,60):
			tpari = [None]*ns
			for j in xrange(ns):  tpari[j] = params[i][j]
			if(rphi == 60):
				for j in xrange(ns):  tpari[j][0] = (tpari[j][0]+60)%symphi

			# mirror checking
			psi_diff = angle_diff( [tpari[j][2] for j in xrange(ns)], [params[0][j][2] for j in xrange(ns)] )
			#print  psi_diff
			if(abs(psi_diff-180.0) <90.0):
				#mirror
				#print "  MIRROR "
				# For each projection direction from the reference set (here zero) find the nearest reduced from the other set
				# and replace it in the other set
				temp = [[0.,0.,0.,0.,0.] for j in xrange(ns)]
				for j in xrange(ns):
					apixer = -1.e20
					qt = Transform({"type":"spider","phi":tpari[j][0], "theta":tpari[j][1], "psi":180.0+tpari[j][2]})
					qt.set_trans(Vec2f(-tpari[j][3], -tpari[j][4]))
					k = 0
					while(k < ks):
						ut = qt*ts[k]
						ut = ut*mm
						bt = ut.get_params("spider")

						tp = bt["phi"]
						tt = bt["theta"]
						if(tt > 90.0):   mp = (tp+180.0)%360.0
						else:            mp = tp
						if( (mp>=badb and mp<bade) or (mp>=bbdb and mp<bbde) ): k = ks
						else: k += 1
					temp[j][0] = bt["phi"]
					temp[j][1] = bt["theta"]
					temp[j][2] = bt["psi"]
					temp[j][3] = -bt["tx"]
					temp[j][4] = -bt["ty"]
				solvs.append([discangset(temp, vt0, sym), temp, [rphi,"psidiff"]])

			#  check the other possibility of mirroring
			p2 = [None]*ns
			for j in xrange(ns):
				p2[j] = [ (-tpari[j][0])%360.0, (180-tpari[j][1])%360.0, tpari[j][2] ,tpari[j][3], tpari[j][4]]
			for j in xrange(ns):  p2[j] = mult_transform(p2[j], mm)
			#  Now check whether p2 is closer to params[0] than tpari is.
			per1 = discangset(tpari, vt0, sym)
			per2 = discangset(p2, vt0, sym)
			#print "  other mirror ",per1,per2
			#for j in xrange(ns):  print p2[j][:3]

			if(per2<per1):
				temp = [[0.,0.,0.,0.,0.] for j in xrange(ns)]
				for j in xrange(ns):
					qt = Transform({"type":"spider","phi":p2[j][0], "theta":p2[j][1], "psi":p2[j][2]})
					qt.set_trans(Vec2f(-p2[j][3], -p2[j][4]))
					k = 0
					while(k < ks):
						ut = qt*ts[k]
						bt = ut.get_params("spider")

						tp = bt["phi"]
						tt = bt["theta"]
						if(tt > 90.0):   mp = (tp+180.0)%360.0
						else:            mp = tp
						if( (mp>=badb and mp<bade) or (mp>=bbdb and mp<bbde) ): k = ks
						else: k += 1
					temp[j][0] = bt["phi"]
					temp[j][1] = bt["theta"]
					temp[j][2] = bt["psi"]
					temp[j][3] = -bt["tx"]
					temp[j][4] = -bt["ty"]
				solvs.append([discangset(temp, vt0, sym), temp, [rphi,"mirror"]])

			# For each projection direction from the reference set (here zero) find the nearest reduced from the other set,
			#    but it cannot be mirrored
			# and replace it in the other set
			temp = [[0.,0.,0.,0.,0.] for j in xrange(ns)]
			for j in xrange(ns):
				apixer = -1.e20
				qt = Transform({"type":"spider","phi":tpari[j][0], "theta":tpari[j][1], "psi":tpari[j][2]})
				qt.set_trans(Vec2f(-tpari[j][3], -tpari[j][4]))
				k = 0
				while(k < ks):
					ut = qt*ts[k]
					bt = ut.get_params("spider")

					tp = bt["phi"]
					tt = bt["theta"]
					if(tt > 90.0):   mp = (tp+180.0)%360.0
					else:            mp = tp
					if( (mp>=badb and mp<bade) or (mp>=bbdb and mp<bbde) ): k = ks
					else: k += 1
				temp[j][0] = bt["phi"]
				temp[j][1] = bt["theta"]
				temp[j][2] = bt["psi"]
				temp[j][3] = -bt["tx"]
				temp[j][4] = -bt["ty"]
			solvs.append([discangset(temp, vt0, sym), temp, [rphi,"straight"]])

		solvs.sort(reverse=False)
		#for k in xrange(len(solvs)):
		#	print  "  SOLVS  ",solvs[k][0],solvs[k][-1]
		#	#for j in xrange(ns):  print solvs[0][1][j][:3]
		psi_diff = angle_diff( [temp[j][2] for j in xrange(ns)], [params[0][j][2] for j in xrange(ns)] )
		if(abs(psi_diff-180.0) <90.0): temp[j][2] = (temp[j][2]+180.0)%360.0
		for j in xrange(ns):  params[i][j] = solvs[0][1][j]
	
def get_dsym_angles(p1, sym):
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
			print  " s ",q["phi"], q["theta"], q["psi"],-q["tx"],-q["tx"]


	mm = Transform({"type":"spider","phi":180., "theta":0., "psi":0.})

	for i in xrange(len(p1)):
		a = Transform({"type":"spider","phi":-p1[i][0], "theta":180-p1[i][1], "psi":p1[i][2]})
		a.set_trans(Vec2f(-p1[i][3], -p1[i][4]))
		a = mm*a
		for l in xrange(len(t)):
			q = a*t[l]
			q = q.get_params("spider")
			print  " m ",q["phi"], q["theta"], q["psi"],-q["tx"],-q["tx"]



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
	max_iter    = int(ali3d_options.maxit1)
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
					if(sym[0] != "d"):  orient_params([params_0, params], subset, sym)
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


def proj_ali_incore_multi(data, refrings, numr, xrng = 0.0, yrng = 0.0, step=1.0, an = 1.0, nsoft = -1, finfo=None):
	from utilities    import compose_transform2
	from math         import cos, pi, radians, degrees
	from EMAN2 import Vec2f, Transform
	from global_def import Util

	mode = "F"
	nx   = data.get_xsize()
	ny   = data.get_ysize()
	#  center is in SPIDER convention
	cnx  = nx//2 + 1
	cny  = ny//2 + 1
	ant = cos(radians(an))

	#phi, theta, psi, sxo, syo = get_params_proj(data)
	t1 = data.get_attr("xform.projection")
	dp = t1.get_params("spider")
	if finfo:
		ID = data.get_attr("ID")
		finfo.write("Image id: %6d\n"%(ID))
		finfo.write("Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"]))
		finfo.flush()
	#print "Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"])
	#[ang, sxs, sys, mirror, iref, peak, checked_refs] = Util.shc(data, refrings, xrng, yrng, step, ant, mode, numr, cnx+dp["tx"], cny+dp["ty"])
	peaks = Util.multiref_polar_ali_2d_peaklist_local(data, refrings, xrng, yrng, step, ant, mode, numr, cnx+dp["tx"], cny+dp["ty"])
	peaks_count = len(peaks) / 5
	#pixel_error = 0.0
	peak = 0.0
	if( peaks_count > 0 ):
		if( nsoft == -1 ):  nsoft = peaks_count
		params = [None]*peaks_count
		#                                              peak         iref      ang  sxs  sys 
		for i in xrange(peaks_count):  params[i] = [ peaks[i*5+0], int(peaks[i*5+4]), peaks[i*5+1], peaks[i*5+2], peaks[i*5+3]]
		params.sort(reverse=True)
		if(peaks_count < nsoft ):
			for i in xrange(peaks_count,nsoft,1): params.insert(0,params[0])
			peaks_count = nsoft
		elif( peaks_count > nsoft ):  peaks_count = min(peaks_count, nsoft)
		ws = sum([params[i][0] for i in xrange(peaks_count)])
		for i in xrange(peaks_count):
			iref   = params[i][1]
			ang    = params[i][2]
			sxs    = params[i][3]
			sys    = params[i][4]
			#mirror = 0
			peak   = params[i][0]/ws
			# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
			# What that means is that one has to change the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
			angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1, -ang, 0.0, 0.0, 1)
			"""
			if  mirror:
				phi   = (refrings[iref].get_attr("phi")+540.0)%360.0
				theta = 180.0-refrings[iref].get_attr("theta")
				psi   = (540.0-refrings[iref].get_attr("psi")+angb)%360.0
				s2x   = sxb - dp["tx"]
				s2y   = syb - dp["ty"]
			else:
			"""
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
			if i == 0:
				data.set_attr("weight", peak)
			else:
				data.set_attr("weight" + str(i), peak)
			#from pixel_error import max_3D_pixel_error
			#pixel_error += max_3D_pixel_error(t1, t2, numr[-3])
			if finfo:
				finfo.write( "New parameters: %9.4f %9.4f %9.4f %9.4f %9.4f %10.5f  %11.3e\n\n" %(phi, theta, psi, s2x, s2y, peak, pixel_error))
				finfo.flush()
			#print  "New parameters: %9.4f %9.4f %9.4f %9.4f %9.4f %10.5f  %11.3e\n\n" %(phi, theta, psi, s2x, s2y, peak, pixel_error)

		# remove old xform.projection
		i = max(peaks_count, 1)
		while data.has_attr("xform.projection" + str(i)):
			data.del_attr("xform.projection" + str(i))
			i += 1
		i = max(peaks_count, 1)
		while data.has_attr("weight" + str(i)):
			data.del_attr("weight" + str(i))
			i += 1
		#pixel_error /= peaks_count
	return ws
	"""
		peak = peaks[0]  # It is not used anywhere, but set it to the maximum.

	return peak, pixel_error, peaks_count, ws
	"""

def shc_multi(data, refrings, numr, xrng, yrng, step, an, nsoft, sym, finfo=None):
	from utilities    import compose_transform2
	from fundamentals import mirror
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

	if finfo:
		t1 = data.get_attr("xform.projection")
		dp = t1.get_params("spider")
		finfo.write("Image id: %6d\n"%(ID))
		finfo.write("Old parameters: %9.4f %9.4f %9.4f %5.2f %5.2f   %10.3f\n"%(dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"], \
					data.get_attr("previousmax")))
		
		i = 1
		while data.has_attr("xform.projection" + str(i)):
			t1 = data.get_attr("xform.projection" + str(i))
			dp = t1.get_params("spider")
			finfo.write("Add parameters: %9.4f %9.4f %9.4f %5.2f %5.2f   %10.3f\n"%(dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"]))
			i += 1

	t1 = data.get_attr("xform.projection")

	#[ang, sxs, sys, mir, iref, peak, checked_refs] = Util.shc(data, refrings, xrng, yrng, step, ant, mode, numr, cnx+dp["tx"], cny+dp["ty"])
	#peaks = Util.shc_multipeaks(data, refrings, xrng, yrng, step, ant, mode, numr, cnx+dp["tx"], cny+dp["ty"], nsoft)
	#  Do not shift the image to prevent sliding away
	peaks = Util.shc_multipeaks(data, refrings, xrng, yrng, step, ant, mode, numr, cnx, cny, nsoft)
	peaks_count = len(peaks) / 7
	pixel_error = 0.0
	number_of_checked_refs = 0
	peak = 0.0
	if( peaks_count > 0 ):
		params = [None]*peaks_count
		#                                              peak         iref                  ang        sxs           sys           mir           checked references
		for i in xrange(peaks_count):  params[i] = [ peaks[i*7+5], int(peaks[i*7+4]), peaks[i*7+0], peaks[i*7+1], peaks[i*7+2], int(peaks[i*7+3]), int(peaks[i*7+6])]
		#  Make sure nothing is repeated
		if(peaks_count>1):
			taken = [params[k][1] for k in xrange(peaks_count)]
			from utilities import findall
			i = 0
			while(i<peaks_count):
				ll = findall(taken[i], taken)
				if(len(ll) > 1):
					print  "  PROBLEM, found the same orientation more than once !  "
					for k in xrange(len(params)):  print  params[k]
					ll.sort(reverse=True)
					for k in xrange(0,len(ll)-1):
						del params[k]
						peaks_count -= 1
					taken = [params[k][1] for k in xrange(peaks_count)]
				i+=1
		params.sort(reverse=True)
		ws = sum([params[i][0] for i in xrange(peaks_count)])  # peaks could be stretched
		for i in xrange(peaks_count):
			ang    = params[i][2]
			sxs    = params[i][3]
			sys    = params[i][4]
			mir    = params[i][5]
			iref   = params[i][1]
			#peak   = peaks[i*7+5]
			#checked_refs = int(peaks[i*7+6])
			#number_of_checked_refs += checked_refs
			#if(sxs>0.0 or sys >0.0):  print  "  SERROR in shc_multi  ",i,params[i]

			# The ormqip returns parameters such that the transformation is applied first, the mir operation second.
			# What that means is that one has to change the the Eulerian angles so they point into mired direction: phi+180, 180-theta, 180-psi
			angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1, -ang, 0.0, 0.0, 1)
			if  mir:
				phi   = (refrings[iref].get_attr("phi")+540.0)%360.0
				theta = 180.0-refrings[iref].get_attr("theta")
				psi   = (540.0-refrings[iref].get_attr("psi")+angb)%360.0
				s2x   = sxb #- dp["tx"]
				s2y   = syb #- dp["ty"]
			else:
				phi   = refrings[iref].get_attr("phi")
				theta = refrings[iref].get_attr("theta")
				psi   = (refrings[iref].get_attr("psi")+angb+360.0)%360.0
				s2x   = sxb #- dp["tx"]
				s2y   = syb #- dp["ty"]
			#if(sxs>0.0 or sys >0.0):  print  "  SERROR2 in shc_multi  ",i,phi,theta,psi,s2x,s2y

			t2 = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
			t2.set_trans(Vec2f(-s2x, -s2y))
			#print i,phi,theta,psi
			if i == 0:
				data.set_attr("xform.projection", t2)
				data.set_attr("weight", params[i][0]/ws)
			else:
				data.set_attr("xform.projection" + str(i), t2)
				data.set_attr("weight" + str(i), params[i][0]/ws)
			from pixel_error import max_3D_pixel_error
			pixel_error += max_3D_pixel_error(t1, t2, numr[-3])
			#  preserve params, they might be needed if peaks_count<nsoft
			params[i] = [params[i][0], phi, theta, psi, s2x, s2y, iref]

		# Now set previousmax to a value halfway through
		data.set_attr("previousmax", params[peaks_count//2][0])

		#  Add orientations around the main peak with exclusion of those already taken
		#  Allow matching to bsoft>nsoft, but still keep nsoft.  This should allow elongated neighborhood
		bsoft = 2*nsoft
		if(peaks_count<nsoft):
			tempref = [refrings[i] for i in xrange(len(refrings))]
			taken   = [params[i][6] for i in xrange(peaks_count)]
			taken.sort(reverse=True)
			if(len(taken) > 1):
				for k in xrange(1,len(taken)):
					dod = []
					if( taken[k] == taken[k-1] ):
						print  "  PROBLEM 2, entries duplicated  ",taken
						dod.append(k)
				if(len(dod) >0):
					for k in dod:  del taken[k]
			#  delete taken
			try:
				for i in xrange(peaks_count):  del  tempref[taken[i]]
			except:
				print  "  failed deleting tempref "
				print i,peaks_count,nsoft
				print  " taken ",taken
				print  len(tempref), len(refrings)
				from sys import exit
				exit()

			from utilities import getfvec
			t1 = data.get_attr("xform.projection")
			dp = t1.get_params("spider")
			n1,n2,n3 = getfvec(dp["phi"],dp["theta"])
			datanvec = [n1,n2,n3]
			if(int(sym[1:]) >1):
				iq = len(tempref)
				iq3 = 3*iq
				iq6 = 6*iq
				tempref += (tempref+tempref)
				refvecs = [None]*3*iq3
				dphi = 360.0/int(sym[1:])
				for i in xrange(iq):
					phi   = tempref[i].get_attr("phi")
					theta = tempref[i].get_attr("theta")
					n1,n2,n3 = getfvec(phi-dphi,theta)
					refvecs[3*i+0] = n1
					refvecs[3*i+1] = n2
					refvecs[3*i+2] = n3
					n1,n2,n3 = getfvec(phi+dphi,theta)
					refvecs[3*i+0+iq6] = n1
					refvecs[3*i+1+iq6] = n2
					refvecs[3*i+2+iq6] = n3
				for i in xrange(iq,2*iq):
					n1 = tempref[i].get_attr("n1")
					n2 = tempref[i].get_attr("n2")
					n3 = tempref[i].get_attr("n3")
					refvecs[3*i+0] = n1
					refvecs[3*i+1] = n2
					refvecs[3*i+2] = n3
			else: 
				refvecs = [None]*3*len(tempref)
				for i in xrange(len(tempref)):
					n1 = tempref[i].get_attr("n1")
					n2 = tempref[i].get_attr("n2")
					n3 = tempref[i].get_attr("n3")
					refvecs[3*i+0] = n1
					refvecs[3*i+1] = n2
					refvecs[3*i+2] = n3
			from utilities import nearestk_to_refdir
			nrst = nearestk_to_refdir(refvecs, datanvec, howmany = bsoft-peaks_count)
			del refvecs
			#  it does not use mir, do it by hand
			if( dp["theta"] > 90.0 ):  tdata = mirror(data)
			else:                      tdata = data.copy()
			#  delete from tdata higher xform and weight and keep only base one as it will be used to do orientation search.
			#  In addition, zero shifts as here we always search around the origin to prevent sliding away.
			i = 1
			while tdata.has_attr("xform.projection" + str(i)):
				tdata.del_attr("xform.projection" + str(i))
				i += 1
			i = 1
			while tdata.has_attr("weight" + str(i)):
				tdata.del_attr("weight" + str(i))
				i += 1
			#  Search
			#if( dp["theta"] > 90.0 ): 
			#	print "  IS MIRRORED  ",dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"]
			pws = proj_ali_incore_multi(tdata, [tempref[k] for k in nrst], numr, xrng, yrng, step, 180.0, bsoft-peaks_count)
			#  Can there be a problem with (0,0) direction??  PAP  05/25/2014
			for i in xrange(bsoft-peaks_count):
				if i == 0:    t1 = tdata.get_attr("xform.projection")
				else:         t1 = tdata.get_attr("xform.projection" + str(i))
				d = t1.get_params("spider")
				phi   = d["phi"]
				theta = d["theta"]
				psi   = d["psi"]
				s2x   = d["tx"]
				s2y   = d["ty"]
				if( dp["theta"] > 90.0 ):
					#  Change parameters if mirrored
					#print "  BEFORE MIRRORED  ",i,phi, theta, psi, s2x, s2y
					phi   = (phi+540.0)%360.0
					theta = 180.0-theta
					psi   = (540.0-psi)%360.0
					#print "  AFTER MIRRORED  ",i,phi, theta, psi, s2x, s2y
				if i == 0 :   w = tdata.get_attr("weight")
				else:         w = tdata.get_attr("weight" + str(i))
				w *= pws  # remove normalization
				params.append([w,  phi, theta, psi, s2x, s2y])
			#  From now on process nsfot largest
			params.sort(reverse=True)
			ws = sum([params[i][0] for i in xrange(nsoft)])  # peaks could be stretched
			for i in xrange(nsoft):
				#print  "  ADDITIONAL SOFT ASSIGNMENT  ",i,peaks_count,params[i][1],params[i][1],params[i][2],params[i][0]/ws
				t2 = Transform({"type":"spider","phi":params[i][1],"theta":params[i][2],"psi":params[i][3]})
				t2.set_trans(Vec2f(-params[i][4], -params[i][5]))
				#print i,phi,theta,psi
				if i == 0:
					data.set_attr("xform.projection", t2)
					data.set_attr("weight", params[i][0]/ws)
				else:
					data.set_attr("xform.projection" + str(i), t2)
					data.set_attr("weight" + str(i), params[i][0]/ws)

		if finfo:
			t1 = data.get_attr("xform.projection")
			dp = t1.get_params("spider")
			#finfo.write("Image id: %6d\n"%(ID))
			finfo.write("New parameters: %9.4f %9.4f %9.4f %5.2f %5.2f  %10.3f  %5.3f\n"%(dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"], \
						data.get_attr("previousmax"), data.get_attr("weight")))
			i = 1
			while data.has_attr("xform.projection" + str(i)):
				t1 = data.get_attr("xform.projection" + str(i))
				dp = t1.get_params("spider")
				finfo.write("Add parameters: %9.4f %9.4f %9.4f %5.2f %5.2f         %5.3f\n"%(dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"], \
						data.get_attr("weight" + str(i)) ) )
				finfo.flush()
				i += 1
		"""
		# remove old xform.projection
		i = max(peaks_count, 1)
		while data.has_attr("xform.projection" + str(i)):
			data.del_attr("xform.projection" + str(i))
			i += 1
		i = max(peaks_count, 1)
		while data.has_attr("weight" + str(i)):
			data.del_attr("weight" + str(i))
			i += 1
		"""
		pixel_error /= peaks_count
		peak = params[0][0]  # It is not used anywhere, but set it to the maximum.
	
	#  if it did not find any higher peaks would do nothing and return peaks_count=0
	return peak, pixel_error, number_of_checked_refs, peaks_count


# parameters: list of (all) projections | reference volume | ...
def ali3d_multishc_soft(stack, ref_vol, ali3d_options, mpi_comm = None, log = None, nsoft=2 ):

	from alignment       import Numrinit, prepare_refrings, proj_ali_incore_local
	from utilities       import get_im, file_type, model_circle, get_input_from_string, get_params_proj, wrap_mpi_gatherv, wrap_mpi_bcast
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier, mpi_reduce, MPI_INT, MPI_SUM
	from projection      import prep_vol
	from statistics      import hist_list
	from applications    import MPI_start_end
	from filter          import filt_ctf
	from global_def      import Util
	from EMAN2           import EMUtil, EMData
	import types
	from time            import time

	ir     = ali3d_options.ir
	rs     = ali3d_options.rs
	ou     = ali3d_options.ou
	xr     = ali3d_options.xr
	yr     = ali3d_options.yr
	ts     = ali3d_options.ts
	an     = ali3d_options.an
	sym    = ali3d_options.sym
	sym    = sym[0].lower() + sym[1:]
	delta  = ali3d_options.delta
	center = ali3d_options.center
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
		log.add("Start ali3d_multishc_soft")

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
	max_iter    = int(ali3d_options.maxit)
	center      = int(center)

	if( type(ref_vol) is types.StringType ):  vol = get_im(ref_vol)
	else:	vol = ref_vol
	nx      = vol.get_xsize()
	if last_ring < 0:	last_ring = int(nx/2) - 2

	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	mask2D  = model_circle(last_ring,nx,nx) - model_circle(first_ring,nx,nx)

	if( type(stack) is types.StringType ):
		if myid == main_node:
			if file_type(stack) == "bdb":
				from EMAN2db import db_open_dict
				dummy = db_open_dict(stack, True)
			active = EMUtil.get_all_attributes(stack, 'active')
			list_of_particles = []
			for im in xrange(len(active)):
				if active[im]:  list_of_particles.append(im)
			del active
			total_nima = len(list_of_particles)
		else:
			list_of_particles = None
			total_nima = 0

	else:
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

	if( type(stack) is types.StringType ):  data = EMData.read_images(stack, list_of_particles)
	else:                                   data = [ stack[im] for im in list_of_particles ]
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
	#par_r = [[] for im in list_of_particles ]
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
				log.add("ITERATION #%3d,  inner iteration #%3d"%(total_iter, Iter))
				log.add("Delta = %4.1f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f\n"%(delta[N_step], an[N_step], xrng[N_step],yrng[N_step],step[N_step]))
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
			if total_iter == 1:
				# adjust params to references, calculate psi+shifts, calculate previousmax
				for im in xrange(nima):
					previousmax = data[im].get_attr_default("previousmax", -1.0e23)
					if(previousmax == -1.0e23):
						peak, pixer[im] = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],10.0)
						data[im].set_attr("previousmax", peak*0.9)
				if myid == main_node:
					log.add("Time to calculate first psi+shifts+previousmax: %f\n" % (time()-start_time))
					start_time = time()
			#=========================================================================

			mpi_barrier(mpi_comm)
			if myid == main_node:
				start_time = time()
			#=========================================================================
			# alignment
			#number_of_checked_refs = 0
			par_r = [0]*(nsoft+1)
			for im in xrange(nima):
				peak, pixer[im], checked_refs, number_of_peaks = shc_multi(data[im], refrings, numr, xrng[N_step], yrng[N_step], step[N_step],\
																			an[N_step], nsoft, sym)
				#number_of_checked_refs += checked_refs
				par_r[number_of_peaks] += 1
				#print  myid,im,number_of_peaks
				#t = get_params_proj(data[im])
				#if(t[3] >0.0 or t[4]>0.0):  print  "  MERRROR  ",t
				
			#=========================================================================
			mpi_barrier(mpi_comm)
			if myid == main_node:
				#print  data[0].get_attr_dict()
				log.add("Time of alignment = %f\n"%(time()-start_time))
				start_time = time()

			#=========================================================================
			#output pixel errors, check stop criterion
			all_pixer = wrap_mpi_gatherv(pixer, 0, mpi_comm)
			par_r = mpi_reduce(par_r, nsoft+1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
			#total_checked_refs = wrap_mpi_gatherv([number_of_checked_refs], main_node, mpi_comm)
			terminate = 0
			if myid == main_node:
				#total_checked_refs = sum(total_checked_refs)
				log.add("=========== Number of better peaks found ==============")
				for lhx in xrange(nsoft+1):
					msg = "            %5d     %7d"%(lhx, par_r[lhx])
					log.add(msg)
				log.add("_______________________________________________________")

				lhist = 20
				region, histo = hist_list(all_pixer, lhist)
				log.add("=========== Histogram of pixel errors ==============")
				for lhx in xrange(lhist):
					msg = "          %10.3f     %7d"%(region[lhx], histo[lhx])
					log.add(msg)
				log.add("____________________________________________________")
				if (max(all_pixer) < 0.5) and (sum(all_pixer)/total_nima < 0.05):
					terminate = 1
					log.add("...............")
					log.add(">>>>>>>>>>>>>>>   Will terminate due to small pixel errors")
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
			if myid == main_node:  vol.write_image('soft/smvol%04d.hdf'%total_iter)
			# log
			if myid == main_node:
				log.add("3D reconstruction time = %f\n"%(time()-start_time))
				start_time = time()
			#=========================================================================

			#=========================================================================
			if(total_iter%1 == 0 or terminate):
				# gather parameters
				params = []
				previousmax = []
				for im in data:
					t = get_params_proj(im)
					params.append( [t[0], t[1], t[2], t[3], t[4]] )
					#if(t[3] >0.0 or t[4]>0.0):  print  "  ERRROR  ",t
					previousmax.append(im.get_attr("previousmax"))
				assert(nima == len(params))
				params = wrap_mpi_gatherv(params, 0, mpi_comm)
				if myid == 0:
					assert(total_nima == len(params))
				previousmax = wrap_mpi_gatherv(previousmax, 0, mpi_comm)
				if myid == main_node:
					from utilities import write_text_row, write_text_file
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

						params.append( [t[0], t[1], t[2], t[3], t[4]] )
						#if(t[3] >0.0 or t[4]>0.0):  print  "  ERRROR  ",i,t
					assert(nima == len(params))
					params = wrap_mpi_gatherv(params, 0, mpi_comm)
					if myid == 0:
						assert(total_nima == len(params))
					if myid == main_node:
						write_text_row(params, "soft/params-%04d-%04d.txt"%(i,total_iter))
					del previousmax, params
					i+=1


	if myid == main_node:
		log.add("Finish ali3d_multishc_soft")
		return #params, vol, previousmax, par_r
	else:
		return #None, None, None, None  # results for the other processes

def get_softy(im):
	w = [im.get_attr('weight')]
	from utilities import get_params_proj
	p1,p2,p3,p4,p5 = get_params_proj(im)
	x = [[p1,p2,p3,p4,p5]]
	i = 1
	while im.has_attr("xform.projection" + str(i)):
		w.append(im.get_attr("weight" + str(i)))
		p1,p2,p3,p4,p5 = get_params_proj(im, "xform.projection" + str(i))
		x.append( [p1,p2,p3,p4,p5] )
		i += 1
	return w,x


# data - projections (scattered between cpus)
# options - the same for all cpus
# return - volume the same for all cpus
def do_volume(data, options, iter, mpi_comm):
	from EMAN2          import Util
	from mpi            import mpi_comm_rank
	from reconstruction import recons3d_4nn_MPI, recons3d_4nn_ctf_MPI
	from utilities      import bcast_EMData_to_all
	
	myid = mpi_comm_rank(mpi_comm)
	sym  = options.sym
	sym = sym[0].lower() + sym[1:]
	npad      = options.npad
	CTF       = options.CTF
	snr       = options.snr
	#=========================================================================
	# volume reconstruction
	if CTF: vol = recons3d_4nn_ctf_MPI(myid, data, snr, symmetry=sym, npad=npad, mpi_comm=mpi_comm)
	else:   vol = recons3d_4nn_MPI    (myid, data,      symmetry=sym, npad=npad, mpi_comm=mpi_comm)

	if myid == 0:
		from morphology import threshold
		from filter     import filt_tanl, filt_btwl
		from utilities  import model_circle, get_im
		nx = data[0].get_xsize()
		last_ring   = int(options.ou)
		if(options.mask3D == None):	mask3D = model_circle(last_ring, nx, nx, nx)
		elif(options.mask3D == "auto"):
			from utilities import adaptive_mask
			mask3D = adaptive_mask(vol)
		else:						mask3D = get_im(options.mask3D)
		stat = Util.infomask(vol, mask3D, True)
		vol -= stat[0]
		Util.mul_scalar(vol, 1.0/stat[1])
		vol = threshold(vol)
		Util.mul_img(vol, mask3D)
		if( options.pwreference ):
			from utilities    import read_text_file
			from fundamentals import rops_table, fftip, fft
			from filter       import filt_table
			rt = read_text_file( options.pwreference )
			fftip(vol)
			ro = rops_table(vol)
			#  Here unless I am mistaken it is enough to take the beginning of the reference pw.
			for i in xrange(1,len(ro)):  ro[i] = (rt[i]/ro[i])**0.5
			vol = fft( filt_table( filt_tanl(vol, options.fl, options.aa), ro) )
		else:  vol = filt_tanl(vol, options.fl, options.aa)
		stat = Util.infomask(vol, mask3D, True)
		vol -= stat[0]
		Util.mul_scalar(vol, 1.0/stat[1])
		vol = threshold(vol)
		volf = filt_btwl(volf, 0.38, 0.5)
		Util.mul_img(vol, mask3D)
		del mask3D
		vol.write_image('toto%03d.hdf'%iter)
	# broadcast volume
	bcast_EMData_to_all(vol, myid, 0, comm=mpi_comm)
	#=========================================================================
	return vol


# parameters: list of (all) projections | reference volume | ...
#  Add reduction
def ali3d_base(stack, ref_vol, ali3d_options, shrinkage = 1.0, mpi_comm = None, log = None, nsoft = 3 ):

	from alignment       import Numrinit, prepare_refrings, proj_ali_incore,  proj_ali_incore_local, shc
	from utilities       import bcast_number_to_all, bcast_EMData_to_all, 	wrap_mpi_gatherv, wrap_mpi_bcast, model_blank
	from utilities       import get_im, file_type, model_circle, get_input_from_string, get_params_proj, set_params_proj
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier, mpi_reduce, MPI_INT, MPI_SUM
	from projection      import prep_vol
	from statistics      import hist_list
	from applications    import MPI_start_end
	from filter          import filt_ctf
	from global_def      import Util
	from fundamentals    import resample
	from multi_shc       import do_volume, shc_multi
	from EMAN2           import EMUtil, EMData
	import types
	from time            import time

	ir     = ali3d_options.ir
	rs     = ali3d_options.rs
	ou     = ali3d_options.ou
	xr     = ali3d_options.xr
	yr     = ali3d_options.yr
	ts     = ali3d_options.ts
	an     = ali3d_options.an
	sym    = ali3d_options.sym
	sym    = sym[0].lower() + sym[1:]
	delta  = ali3d_options.delta
	center = ali3d_options.center
	CTF    = ali3d_options.CTF
	ref_a  = ali3d_options.ref_a

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	if log == None:
		from logger import Logger
		log = Logger()

	number_of_proc = mpi_comm_size(mpi_comm)
	myid           = mpi_comm_rank(mpi_comm)
	main_node      = 0

	if myid == main_node:
		log.add("Start ali3d_base")

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
	max_iter    = int(ali3d_options.maxit)
	center      = int(center)
		
	"""
	if( type(ref_vol) is types.StringType ):
		if(myid == main_node):
			vol = get_im(ref_vol)
	else:
		vol = ref_vol

	if(myid == main_node):
		onx      = vol.get_xsize()
		if(shrinkage == 1.0):  nx = onx
		else:		
			vol = resample(vol, shrinkage)
			nx = vol.get_xsize()
		from filter import filt_tanl
		vol = filt_tanl(vol, 0.45, 0.1)  #  PREFILTER!!
	else:
		nx = 0
		onx = 0
	nx  = bcast_number_to_all(nx, source_node = main_node)
	onx = bcast_number_to_all(onx, source_node = main_node)
	if last_ring < 0:	last_ring = int(onx/2) - 2
	if(shrinkage < 1.0):
		first_ring = max(1, int(first_ring*shrinkage))
		last_ring  = int(last_ring*shrinkage)
		ali3d_options.ou = last_ring
		ali3d_options.ir = first_ring
	if(myid != main_node):  vol = model_blank(nx, nx, nx)
	mpi_barrier(mpi_comm)
	bcast_EMData_to_all(vol, myid, main_node)

	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	mask2D  = model_circle(last_ring,nx,nx) - model_circle(first_ring,nx,nx)
	"""

	if( type(stack) is types.StringType ):
		if myid == main_node:
			if file_type(stack) == "bdb":
				from EMAN2db import db_open_dict
				dummy = db_open_dict(stack, True)
			active = EMUtil.get_all_attributes(stack, 'active')
			list_of_particles = []
			for im in xrange(len(active)):
				if active[im]:  list_of_particles.append(im)
			del active
			total_nima = len(list_of_particles)
		else:
			list_of_particles = None
			total_nima = 0

	else:
		if myid == main_node:
			list_of_particles = range(len(stack))
			total_nima = len(list_of_particles)
		else:
			list_of_particles = None
			total_nima = None
	
	total_nima = wrap_mpi_bcast(total_nima, main_node, mpi_comm)
	list_of_particles = wrap_mpi_bcast(list_of_particles, main_node, mpi_comm)
	if myid == main_node:
		particle_ids = [0]*total_nima
		for i in xrange(total_nima):  particle_ids[i] = list_of_particles[i]
	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)

	if myid == 7:
		import os
		outdir = "./"
		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
	else:
		finfo = None

	if(myid == main_node):
		if( type(stack) is types.StringType ):  data = get_im(stack, list_of_particles[0])
		else:                                   data = stack[list_of_particles[0]]
		onx      = data.get_xsize()
		if(shrinkage == 1.0):  nx = onx
		else:		
			st = resample(data, shrinkage)
			nx = st.get_xsize()
	else:
		nx = 0
		onx = 0
	nx  = bcast_number_to_all(nx, source_node = main_node)
	onx = bcast_number_to_all(onx, source_node = main_node)

	if last_ring < 0:	last_ring = int(onx/2) - 2
	mask2D  = model_circle(last_ring,onx,onx) - model_circle(first_ring,onx,onx)
	if(shrinkage < 1.0):
		first_ring = max(1, int(first_ring*shrinkage))
		last_ring  = int(last_ring*shrinkage)
		ali3d_options.ou = last_ring
		ali3d_options.ir = first_ring
	numr	= Numrinit(first_ring, last_ring, rstep, "F")


	data = [None]*nima
	for im in xrange(nima):
		if( type(stack) is types.StringType ):  data[im] = get_im(stack, list_of_particles[im])
		else:                                   data[im] = stack[list_of_particles[im]].copy()
		data[im].set_attr('ID', list_of_particles[im])
		ctf_applied = data[im].get_attr_default('ctf_applied', 0)
		if CTF :
			ctf_params = data[im].get_attr("ctf")
			if ctf_applied == 0:
				st = Util.infomask(data[im], mask2D, False)
				data[im] -= st[0]
				data[im] = filt_ctf(data[im], ctf_params)
				data[im].set_attr('ctf_applied', 1)
		if(shrinkage != 1.0):
			phi,theta,psi,sx,sy = get_params_proj(data[im])
			data[im] = resample(data[im], shrinkage)
			sx *= shrinkage
			sy *= shrinkage
			set_params_proj(data[im], [phi,theta,psi,sx,sy])
			if CTF :
				ctf_params.apix /= shrinkage
				data[im].set_attr('ctf', ctf_params)
	del mask2D

	# Reference volume reconstruction
	mpi_barrier(mpi_comm)
	if myid == main_node:
		start_time = time()
	vol = do_volume(data, ali3d_options, 0, mpi_comm)
	#if myid == main_node:  vol.write_image('soft/smvol%04d.hdf'%total_iter)
	# log
	if myid == main_node:
		log.add("Dimensions used (nx, onx, first_ring, last_ring, shrinkage)  %5d     %5d     %5d     %5d     %6.3f\n"%(nx, onx, first_ring, last_ring, shrinkage))
		log.add("Reference 3D reconstruction time = %f\n"%(time()-start_time))
		start_time = time()


	pixer = [0.0]*nima
	historyofchanges = [0.0, 0.5, 1.0]
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
				log.add("ITERATION #%3d,  inner iteration #%3d"%(total_iter, Iter))
				log.add("Delta = %5.2f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f\n"%(delta[N_step], an[N_step], xrng[N_step], yrng[N_step], step[N_step]))
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
			if total_iter == 1:
				# adjust params to references, calculate psi+shifts, calculate previousmax
				for im in xrange(nima):
					previousmax = data[im].get_attr_default("previousmax", -1.0e23)
					if(previousmax == -1.0e23):
						peak, pixer[im] = proj_ali_incore_local(data[im], refrings, numr, xrng[N_step], yrng[N_step], step[N_step], 10.0, sym = sym)
						data[im].set_attr("previousmax", peak)#*0.9)
				if myid == main_node:
					log.add("Time to calculate first psi+shifts+previousmax: %f\n" % (time()-start_time))
					start_time = time()
			#=========================================================================

			mpi_barrier(mpi_comm)
			if myid == main_node:
				start_time = time()
			#=========================================================================
			# alignment
			#number_of_checked_refs = 0
			par_r = [0]*max(2,(nsoft+1))
			for im in xrange(nima):
				if(nsoft == 0):
					if(an[N_step] == -1): peak, pixer[im] = proj_ali_incore(data[im], refrings, numr, \
														xrng[N_step], yrng[N_step], step[N_step])
					else:                 peak, pixer[im] = proj_ali_incore_local(data[im], refrings, numr, \
														xrng[N_step], yrng[N_step], step[N_step], an[N_step], sym = sym, finfo = finfo)
					if(pixer[im] == 0.0):  par_r[0] += 1
				elif(nsoft == 1):
					peak, pixer[im], number_of_checked_refs, iref = \
						shc(data[im], refrings, numr, xrng[N_step], yrng[N_step], step[N_step], an[N_step], sym)
					if(pixer[im] == 0.0):  par_r[0] += 1
				elif(nsoft > 1):
					peak, pixer[im], checked_refs, number_of_peaks = shc_multi(data[im], refrings, numr, \
												xrng[N_step], yrng[N_step], step[N_step], an[N_step], nsoft, sym, finfo)
					par_r[number_of_peaks] += 1
					#number_of_checked_refs += checked_refs

			#=========================================================================
			mpi_barrier(mpi_comm)
			if myid == main_node:
				#print  data[0].get_attr_dict()
				log.add("Time of alignment = %f\n"%(time()-start_time))
				start_time = time()

			#=========================================================================
			#output pixel errors, check stop criterion
			all_pixer = wrap_mpi_gatherv(pixer, 0, mpi_comm)
			par_r = mpi_reduce(par_r, len(par_r), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
			#total_checked_refs = wrap_mpi_gatherv([number_of_checked_refs], main_node, mpi_comm)
			terminate = 0
			if myid == main_node:
				#total_checked_refs = sum(total_checked_refs)
				if(nsoft < 2):  par_r[1] = total_nima - par_r[0]
				log.add("=========== Number of better peaks found ==============")
				for lhx in xrange(len(par_r)):
					msg = "            %5d     %7d"%(lhx, par_r[lhx])
					log.add(msg)
				log.add("_______________________________________________________")
				saturatecrit = 0.95
				changes = par_r[0]/float(total_nima)
				if(  changes > saturatecrit ):
					terminate = 1
					log.add("...............")
					log.add(">>>>>>>>>>>>>>>   Will terminate as %4.2f images did not find better orientations"%saturatecrit)
				if( terminate == 0 ):
					historyofchanges.append(changes)
					historyofchanges = historyofchanges[:3]
					historyofchanges.sort()
					"""  Have to think about it PAP
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
						if(msg < 2.0): lhx += 1
					if(float(lhx)/float(total_nima) > saturatecrit):
						terminate = 1
						log.add("...............")
						log.add(">>>>>>>>>>>>>>>   Will terminate as %4.2f images had pixel error <2.0"%saturatecrit)
					
				"""
				if (max(all_pixer) < 0.5) and (sum(all_pixer)/total_nima < 0.05):
					terminate = 1
					log.add("...............")
					log.add(">>>>>>>>>>>>>>>   Will terminate due to small pixel errors")
				"""
			terminate = wrap_mpi_bcast(terminate, main_node, mpi_comm)
			#=========================================================================
			mpi_barrier(mpi_comm)
			if myid == main_node:
				#print  data[0].get_attr_dict()
				log.add("Time to compute histograms = %f\n"%(time()-start_time))
				start_time = time()

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
			vol = do_volume(data, ali3d_options, total_iter, mpi_comm)
			#if myid == main_node:  vol.write_image('soft/smvol%04d.hdf'%total_iter)
			# log
			if myid == main_node:
				log.add("3D reconstruction time = %f\n"%(time()-start_time))
				start_time = time()
			#=========================================================================

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
					from utilities import write_text_row, write_text_file
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


			if terminate:
				# gather parameters
				params = []
				for im in data:
					t = get_params_proj(im)
					params.append( [t[0], t[1], t[2], t[3]/shrinkage, t[4]/shrinkage] )
				params = wrap_mpi_gatherv(params, main_node, mpi_comm)
			if( terminate and (myid == main_node) ):
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


	if myid == main_node:
		log.add("Finish ali3d_base")
	i = 1
	return i   #params, vol, previousmax, par_r
	#else:
	#	return #None, None, None, None  # results for the other processes
