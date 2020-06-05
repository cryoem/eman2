












































































"""0
	elif(symmetry_class.sym[0] == "d"):
		# In this one the first params is taken as a reference
		out = deepcopy(params)
		mirror_and_reduce_dsym([refparams,out], tindexes, symmmetry_class.sym)
	else: # o, t, i
	"""



























































































































"""1

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



























































































"""2
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)
	"""














































































































































































































































"""3
			params = []
			for im in data:
				phi, theta, psi, sx, sy = get_params_proj(im)
				params.append([phi, theta, psi, sx, sy])
			params_0 = wrap_mpi_bcast(params, mpi_subroots[0], mpi_comm)

			if mpi_subrank == 0:
				from sp_utilities import write_text_row
				write_text_row(params, "qparams%04d%04d.hdf"%(myid,total_iter) )
			"""























































































































































"""4
							for i in xrange(total_nima):
								print  i,newparms1[i],GA[ipl[ip]][1][i]
							for i in xrange(total_nima):
								print  i,newparms2[i],GA[ipl[ip+1]][1][i]
							for i in xrange(total_nima):
								GA[ipl[ip]][1][i]   = newparms1[i]
								GA[ipl[ip+1]][1][i] = newparms2[i]
							"""






































"""5
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






























"""6
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
"""7
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

















































































"""8
	if an == "-1":
		an = [-1] * lstp
	else:
		an = get_input_from_string(an)
	"""



























































































"""9
		refrings = [refrings]
		n1 = sin(theta*qv)*cos(phi*qv)
		n2 = sin(theta*qv)*sin(phi*qv)
		n3 = cos(theta*qv)
		refrings[0].set_attr_dict( {"n1":n1, "n2":n2, "n3":n3} )
		refrings[0].set_attr("phi",   phi)
		refrings[0].set_attr("theta", theta)
		refrings[0].set_attr("psi",   psi)
		"""










"""10
	from sp_projection   import prep_vol, prgs
	from sp_alignment import ringwe
	cnx = nx//2 + 1
	cny = nx//2 + 1
	wr_four  = ringwe(numr, "F")
	from math import pi, sin, cos
	qv = pi/180.
	#volft, kb = prep_vol(ref_vol)
	from sp_utilities import get_params_proj
	for im in xrange(nima):
		phi,theta,psi,tx,ty = get_params_proj(data[im])
		prjref = prgs(volft, kb, [phi,theta,psi, 0.0, 0.0])
		cimage = Util.Polar2Dm(prjref, cnx, cny, numr, "F")
		Util.Normalize_ring(cimage, numr, 0 )
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
		peak, pixer = proj_ali_incore_local(data[im], refrings, [[phi, theta], [phi, theta]], numr,0.0,0.0,1.0,delta[0]/4)
		data[im].set_attr("previousmax", peak)
		print  "peak ",myid,im,data[im].get_attr("ID"), peak,pixer,get_params_proj(data[im])
	"""



























































































































































































































"""11
Order of functions:
multi_shc - main
  compute initial volume
  refparams2 = ali3d_multishc() - refine on each node, outputs one set of params
  refvol2.hdf = do_volume()  store
  out_params, out_vol, previousmax, out_r = ali3d_multishc_2()  stores everything
  	computes new previousmax
  	stores volume viterb.hdf
  
"""






















































































































"""12
	if mpi_rank == 0:
		assert(len(out_params) == runs_count)
	"""


"""13
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





"""14
		log.add("  WILL RECONSTRUCT  ")
		from sp_utilities import model_circle
		tvol = volume_recsp(projections, ali3d_options)
		LL2 = tvol.cmp("dot", tvol, dict(negative = 0, mask = model_circle(22, 48, 48, 48)))
		log.add(" LLLLLLL2 norm of reference volume:  %f"%LL2)
		for k in xrange(mpi_size):
			proj_begin, proj_end  = MPI_start_end(n_projs, mpi_size, k)
			print  " from within   ",k,proj_begin,get_params_proj(projections[proj_begin])
		"""
"""15
	else:
		temp_projs = None
	"""


















































"""16
	if mpi_rank == 17:
		temp = []
		from sp_utilities import get_params_proj
		for i in xrange(n_projs):
			#projections[i].set_attr("stable", 1)
			t1,t2,t3,t4,t5 = get_params_proj( projections[i])
			t6 = projections[i].get_attr("previousmax")
			temp.append([t1,t2,t3,t4,t5,t6])
			#set_params_proj( projections[i], out_params[i] )
		write_text_row(temp, log.prefix + "refparams17.txt")
	"""




































































































































































"""17
			if  mirror:
				phi   = (refrings[iref].get_attr("phi")+540.0)%360.0
				theta = 180.0-refrings[iref].get_attr("theta")
				psi   = (540.0-refrings[iref].get_attr("psi")+angb)%360.0
				s2x   = sxb - dp["tx"]
				s2y   = syb - dp["ty"]
			else:
			"""


































"""18
		peak = peaks[0]  # It is not used anywhere, but set it to the maximum.

	return peak, pixel_error, peaks_count, ws
	"""






























































































































































































































































"""19
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
































































































































































































































































































































































































































"""20
# Not used anymore
def calculate_matrix_rot(projs):
	from sp_utilities import rotation_between_anglesets
	sc = len(projs)
	matrix_rot  = [[[0.0,0.0,0.0,0.0,0.0] for i in xrange(sc)] for k in xrange(sc)]
	for i in xrange(sc):
		for j in xrange(i):
			t1, t2, t3 = rotation_between_anglesets(projs[i], projs[j])
			matrix_rot[i][j] = [t1, t2, t3, 0.0, 0.0]
	return matrix_rot


def get_softy(im):
	w = [im.get_attr('weight')]
	from sp_utilities import get_params_proj
	p1,p2,p3,p4,p5 = get_params_proj(im)
	x = [[p1,p2,p3,p4,p5]]
	i = 1
	while im.has_attr("xform.projection" + str(i)):
		w.append(im.get_attr("weight" + str(i)))
		p1,p2,p3,p4,p5 = get_params_proj(im, "xform.projection" + str(i))
		x.append( [p1,p2,p3,p4,p5] )
		i += 1
	return w,x




# Not used anywhere
def get_dsym_angles(p1, sym):
	#  works only for d symmetry
	from sp_utilities import get_symt
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

def reduce_dsym_angles(p1, sym):
	#  WARNING - it returns incorrect parameters that are only suitable for calculation of angular distances
	#  works only for d symmetry
	pr = [[0.0 for i in xrange(5)] for q in xrange(len(p1))]
	for i in xrange(len(p1)):
		if( p1[i][1] >90.0):
			p1[i][1] = 180.0 - p1[i][1]
			p1[i][0] = (p1[i][0] +180.0)%360.0
	'''
	from sp_utilities import get_symt
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
	'''
	return pr



# parameters: list of (all) projections | reference volume | ...
def Xali3d_multishc(stack, ref_vol, ali3d_options, mpi_comm = None, log = None, number_of_runs=2 ):

	from sp_alignment    import Numrinit, prepare_refrings, proj_ali_incore_local, shc
	from sp_utilities    import model_circle, get_input_from_string, get_params_proj, set_params_proj, wrap_mpi_gatherv, wrap_mpi_bcast, wrap_mpi_send, wrap_mpi_recv
	from mpi          import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier, mpi_comm_split, mpi_comm_free
	from sp_projection   import prep_vol
	from sp_statistics   import hist_list
	from sp_applications import MPI_start_end
	from sp_filter       import filt_ctf
	from sp_global_def   import Util, ERROR
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
		from sp_logger import Logger
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
						peak, temp = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step], delta[N_step]*0.7, sym=sym )
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
				from sp_utilities      import estimate_3D_center_MPI, rotate_3D_shift
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
				from sp_utilities import write_text_row
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




'''21
# parameters: list of (all) projections | reference volume is optional, if provided might be shrank| ...
#  The alignment done depends on nsoft:
# 			 nsoft = 0 & an = -1: exhaustive deterministic
# 			 nsoft = 0 & an > 0 : local deterministic
# 			 nsoft = 1 shc
# 			 nsoft >1  shc_multi
def ali3d_base(stack, ref_vol = None, ali3d_options = None, shrinkage = 1.0, mpi_comm = None, log = None, nsoft = 3, \
		saturatecrit = 0.95, pixercutoff = 1.0, zoom = False):

	from sp_alignment       import Numrinit, prepare_refrings, proj_ali_incore,  proj_ali_incore_local, shc, center_projections_3D
	from sp_utilities       import bcast_number_to_all, bcast_EMData_to_all, 	wrap_mpi_gatherv, wrap_mpi_bcast, model_blank
	from sp_utilities       import get_im, file_type, model_circle, get_input_from_string, get_params_proj, set_params_proj
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier, mpi_reduce, MPI_INT, MPI_SUM
	from sp_projection      import prep_vol
	from sp_statistics      import hist_list
	from sp_applications    import MPI_start_end
	from sp_filter          import filt_ctf
	from sp_global_def      import Util
	from sp_fundamentals    import resample, fshift
	from sp_multi_shc       import do_volume, shc_multi
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
	#maskfile = ali3d_options.mask3D

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	if log == None:
		from sp_logger import Logger
		log = Logger()

	number_of_proc = mpi_comm_size(mpi_comm)
	myid           = mpi_comm_rank(mpi_comm)
	main_node      = 0

	if myid == main_node:
		log.add("Start ali3d_base, nsoft = %1d"%nsoft)

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

	if( type(stack) is types.StringType ):
		if myid == main_node:
			total_nima = EMUtil.get_image_count( stack )
		else:
			total_nima = 0
		total_nima = wrap_mpi_bcast(total_nima, main_node, mpi_comm)
		list_of_particles = range(total_nime)
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


	if myid == 0:
		finfo = None
		"""
		import os
		outdir = "./"
		info_file = os.path.join(outdir, "progress%04d"%myid)
		finfo = open(info_file, 'w')
		"""
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

	oldshifts = [None]*nima
	data = [None]*nima
	for im in xrange(nima):
		if( type(stack) is types.StringType ):  data[im] = get_im(stack, list_of_particles[im])
		else:                                   data[im] = stack[list_of_particles[im]]
		data[im].set_attr('ID', list_of_particles[im])
		ctf_applied = data[im].get_attr_default('ctf_applied', 0)
		phi,theta,psi,sx,sy = get_params_proj(data[im])
		data[im] = fshift(data[im], sx, sy)
		set_params_proj(data[im],[phi,theta,psi,0.0,0.0])
		#  For local SHC set anchor
		if(nsoft == 1 and an[0] > -1):
			set_params_proj(data[im],[phi,theta,psi,0.0,0.0], "xform.anchor")
		oldshifts[im] = [sx,sy]
		if CTF :
			ctf_params = data[im].get_attr("ctf")
			if ctf_applied == 0:
				st = Util.infomask(data[im], mask2D, False)
				data[im] -= st[0]
				data[im] = filt_ctf(data[im], ctf_params)
				data[im].set_attr('ctf_applied', 1)
		if(shrinkage < 1.0):
			#  resample will properly adjusts shifts and pixel size in ctf
			data[im] = resample(data[im], shrinkage)
	del mask2D
	mpi_barrier(mpi_comm)
	"""
	if maskfile:
		if type(maskfile) is types.StringType:
			if myid == main_node:
				mask3D = get_im(maskfile)
				i = mask3D.get_xsize()
				if( shrinkage != 1.0 ):
					if( i != nx ):
						mask3D = resample(mask3D, shrinkage)
			else:
				mask3D = model_blank(nx, nx, nx)
			bcast_EMData_to_all(mask3D, myid, main_node)
		else:
			mask3D = maskfile
	else:
		mask3D = model_circle(last_ring, nx, nx, nx)
	ali3d_options.mask3D = mask3D
	"""

	if myid == main_node:
		start_time = time()

	#  Read	template volume if provided or reconstruct it
	#  Apply initfl first, meaning true fl has to be preserved
	fl = ali3d_options.fl
	ali3d_options.fl = ali3d_options.initfl
	if ref_vol:
		if type(ref_vol) is types.StringType:
			if myid == main_node:
				vol = get_im(ref_vol)
				i = vol.get_xsize()
				if( shrinkage != 1.0 ):
					if( i != nx ):
						vol = resample(vol, shrinkage)
			else:
				vol = model_blank(nx, nx, nx)
		else:
			if myid == main_node:
				i = ref_vol.get_xsize()
				if( shrinkage != 1.0 ):
					if( i != nx ):
						vol = resample(ref_vol, shrinkage)
				else:
					vol = ref_vol.copy()
			else:
				vol = model_blank(nx, nx, nx)
		bcast_EMData_to_all(vol, myid, main_node)
		del ref_vol
		vol = do_volume(vol, ali3d_options, 0, mpi_comm)
	else:
		vol = do_volume(data, ali3d_options, 0, mpi_comm)
	#  Restore desired fl
	ali3d_options.fl = fl

	# log
	if myid == main_node:
		log.add("Dimensions used (nx, onx, first_ring, last_ring, shrinkage)  %5d     %5d     %5d     %5d     %6.3f\n"%(nx, onx, first_ring, last_ring, shrinkage))
		log.add("Reference 3D reconstruction time = %f\n"%(time()-start_time))
		start_time = time()


	pixer = [0.0]*nima
	historyofchanges = [0.0, 0.5, 1.0]
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
				log.add("Delta = %5.2f, an = %5.2f, xrange = %5.2f, yrange = %5.2f, step = %5.2f\n"%(delta[N_step], an[N_step], xrng[N_step], yrng[N_step], step[N_step]))
				start_time = time()

			#=========================================================================
			# build references
			volft, kb = prep_vol(vol)
			refrings = prepare_refrings(volft, kb, nx, delta[N_step], ref_a, sym, numr, MPI=mpi_comm, phiEqpsi = "Zero")
			del volft, kb
			#=========================================================================

			if myid == main_node:
				log.add("Time to prepare rings: %f\n" % (time()-start_time))
				start_time = time()
			
			#=========================================================================
			#  there is no need for previousmax for deterministic searches
			if total_iter == 1 and nsoft > 0:
				if(an[N_step] < 0.0):
					# adjust params to references, calculate psi+shifts, calculate previousmax
					for im in xrange(nima):
						previousmax = data[im].get_attr_default("previousmax", -1.0e23)
						if(previousmax == -1.0e23):
							peak, pixer[im] = proj_ali_incore_local(data[im], refrings, numr, xrng[N_step], yrng[N_step], step[N_step], delta[N_step]*2.5, sym = sym)
							data[im].set_attr("previousmax", peak)
				else:
					#  Here it is supposed to be shake and bake for local SHC, but it would have to be signaled somehow
					for im in xrange(nima):
						data[im].set_attr("previousmax", -1.0e23)
				if myid == main_node:
					log.add("Time to calculate first psi+shifts+previousmax: %f\n" % (time()-start_time))
					start_time = time()
			#=========================================================================

			mpi_barrier(mpi_comm)
			if myid == main_node:  start_time = time()
			#=========================================================================
			# alignment
			#number_of_checked_refs = 0
			par_r = [0]*max(2,(nsoft+1))
			for im in xrange(nima):
				if(nsoft == 0):
					if(an[N_step] == -1): peak, pixer[im] = proj_ali_incore(data[im], refrings, numr, \
														xrng[N_step], yrng[N_step], step[N_step], sym=sym)
					else:                 peak, pixer[im] = proj_ali_incore_local(data[im], refrings, numr, \
														xrng[N_step], yrng[N_step], step[N_step], an[N_step], finfo = finfo, sym=sym)
					if(pixer[im] == 0.0):  par_r[0] += 1
				elif(nsoft == 1):
					peak, pixer[im], number_of_checked_refs, iref = \
						shc(data[im], refrings, numr, xrng[N_step], yrng[N_step], step[N_step], an[N_step], finfo = finfo, sym=sym)
					if(pixer[im] == 0.0):  par_r[0] += 1
				elif(nsoft > 1):
					peak, pixer[im], checked_refs, number_of_peaks = shc_multi(data[im], refrings, numr, \
												xrng[N_step], yrng[N_step], step[N_step], an[N_step], nsoft, finfo = finfo, sym=sy,)
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
				changes = par_r[0]/float(total_nima)
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
						if(msg < pixercutoff): lhx += 1
					lhx = float(lhx)/float(total_nima)
					log.add(">>> %4.2f images had pixel error <%5.2f"%(lhx,pixercutoff))
					if( lhx > saturatecrit):
						if( Iter == 1 ):
							log.add("Will continue even though %4.2f images had pixel error < %5.2f"%(saturatecrit,pixercutoff))
						else:
							terminate = 1
							log.add("...............")
							log.add(">>>>>>>>>>>>>>>   Will terminate as %4.2f images had pixel error < %5.2f"%(saturatecrit,pixercutoff))
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


			if( terminate or (Iter == max_iter) ):
				# gather parameters
				params = []
				for im in xrange(nima):
					t = get_params_proj(data[im])
					params.append( [t[0], t[1], t[2], t[3]/shrinkage + oldshifts[im][0], t[4]/shrinkage+ oldshifts[im][1]] )
				params = wrap_mpi_gatherv(params, main_node, mpi_comm)
			"""
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
		log.add("Finish ali3d_base, nsoft = %1d"%nsoft)
	return params  #, vol, previousmax, par_r
	#else:
	#	return #None, None, None, None  # results for the other processes

'''
'''22
#  I do not know why this would make sense  PAP 04/20/2017
def generate_uneven_projections_directions(count, half_sphere=False, output_filename=None, 
										   density_ratio_on_poles_and_equator=1.0): 
	"""
	Generates set of uneven distributed projections directions with given size 
	Returns: list of projections directions
	"""
	from random import random, uniform
	from sp_utilities import write_text_row
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
'''
"""23
def mult_transform(v1, v2):
	from EMAN2 import Transform
	T1 = Transform({"type":"spider","phi":v1[0],"theta":v1[1],"psi":v1[2],"tx":v1[3],"ty":v1[4],"tz":0.0,"mirror":0,"scale":1.0})
	T = T1*v2
	return [ T.get_params("spider")["phi"], T.get_params("spider")["theta"], T.get_params("spider")["psi"], T.get_params("spider")["tx"], T.get_params("spider")["ty"]  ]
"""


