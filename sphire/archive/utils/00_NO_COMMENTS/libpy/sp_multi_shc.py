# For some reason these programs get stuck on MPI if I change the order of programs in the file.  Strange, PAP.
from __future__ import print_function
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2019 (toshio.moriya@kek.jp)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
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

def orient_params(params, refparams, indexes=None, symmetry_class = None):
	#
	#  The assumption here is that the angles are within the unique region
	#  Since they come from projection refinement, why would it be otherwise?
	#  Any problems would be due to how even_angles generates reference angles
	#
	#  The function returns rotation object and properly rotated/mirrored params
	from sp_utilities import rotation_between_anglesets
	from sp_fundamentals import rotate_params
	from sp_pixel_error import angle_diff, angle_diff_sym
	from EMAN2 import Transform

	n = len(params)
	if( indexes == None ): tindexes = list(range(n))
	else:  tindexes = indexes

	if(symmetry_class.sym[0] == "c" and symmetry_class.nsym>1):
		from copy import deepcopy
		divic = 360.0/symmetry_class.nsym
		phi = angle_diff_sym([params[j] for j in tindexes], [refparams[j] for j in tindexes], symmetry_class.nsym)
		out = deepcopy(params)
		for j in range(n):
			out[j][0] = (out[j][0]+phi)%divic
		# mirror checking
		psi_diff = angle_diff( [out[j][2] for j in tindexes], [refparams[j][2] for j in tindexes] )
		if(abs(psi_diff-180.0) <90.0):
			for j in range(n):
				# apply mirror
				out[j][2] = (out[j][2] + 180.0) % 360.0
	elif(symmetry_class.sym[0] == "c"):
		t1,t2,t3 = rotation_between_anglesets([params[j] for j in tindexes], [refparams[j] for j in tindexes])
		out = rotate_params([params[i][:3] for i in range(n)],[-t3,-t2,-t1])
		out = [out[i]+params[i][3:]  for i in range(n)]  # reattach shifts
		# mirror checking
		psi_diff = angle_diff( [out[j][2] for j in tindexes], [refparams[j][2] for j in tindexes] )
		if(abs(psi_diff-180.0) <90.0):
			for j in range(n):
				# apply mirror
				out[j][2] = (out[j][2] + 180.0) % 360.0
	"""Multiline Comment0"""
	#MULTILINEMULTILINEMULTILINE 0
		#MULTILINEMULTILINEMULTILINE 0
		#MULTILINEMULTILINEMULTILINE 0
		#MULTILINEMULTILINEMULTILINE 0
	#MULTILINEMULTILINEMULTILINE 0
	#MULTILINEMULTILINEMULTILINE 0
	return  out


def find_common_subset(projs, target_threshold=2.0, minimal_subset_size=3, symmetry_class = None):
	#  projs - [reference set of angles, set of angles1, set of angles2, ... ]
	#  the function is written for multiple sets of angles
	#  The transformation is found for a subset, but applied to the full set
	from sp_utilities import getvec, getfvec, getang3, lacos, angles_to_normals
	from math import acos, degrees
	from copy import deepcopy
	from EMAN2 import Vec2f, Transform

	n  = len(projs[0])
	sc = len(projs)

	subset = list(range(n))

	minimal_subset_size = min( minimal_subset_size, n)

	avg_diff_per_image = [0.0]*n

	if( symmetry_class.sym[:3] == "oct" or symmetry_class.sym[:3] == "tet" or symmetry_class.sym[:4] == "icos" ):
		# In these cases there is no rotation so distances between angular directions cannot change
		outp = deepcopy(projs)
		for i in subset:
			qt = 0.0
			for k in range(sc-1):
				for l in range(k+1,sc):
					neisym = symmetry_class.symmetry_neighbors([outp[l][i][:3]])
					dmin = 180.0
					for q in neisym: dmin = min(dmin, getang3(q, outp[k][i]))
					qt += dmin
			avg_diff_per_image[i] = (qt/sc/(sc-1)/2.0)


	#  Start from the entire set and the slowly decrease it by rejecting worst data one by one.
	#  It will stop either when the subset reaches the minimum subset size 
	#   or if there are no more angles with errors above target threshold
	#
	while(True):
		#  extract images in common subset
		if(symmetry_class.sym[0] == "c"):  #  c including cn symmetry
			divic = 360.0/symmetry_class.nsym
			for i in subset: avg_diff_per_image[i] = 0.0
			outp = [deepcopy(projs[0])]
			for i in range(sc-1,-1,-1):
				tv = angles_to_normals(projs[i])

				for j in range(i+1,sc):
					out = orient_params(projs[j], projs[i], subset, symmetry_class = symmetry_class)
					if(symmetry_class.nsym > 1):
						# for k in xrange(n):
						for k in subset:
							mind = 1.0e23
							for l in range(symmetry_class.nsym):
								u1,u2,u3 = getfvec(out[k][0] + l*divic, out[k][1])
								qt = lacos(tv[k][0]*u1+tv[k][1]*u2+tv[k][2]*u3)
								mind = min(mind, qt)
							avg_diff_per_image[k] += mind
							#print  "avg_diff_per_image  %3d  %8.2f="%(k,avg_diff_per_image[k]),\
							#"  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f"%( projs[i][k][0],projs[i][k][1],projs[i][k][2],out[k][0],out[k][1],out[k][2])
					else:
						# for k in xrange(n):
						for k in subset:
							u1,u2,u3 = getfvec(out[k][0], out[k][1])
							avg_diff_per_image[k] += lacos(tv[k][0]*u1+tv[k][1]*u2+tv[k][2]*u3)
							# print  "avg_diff_per_image  %3d  %8.2f %8.6f="%(k,avg_diff_per_image[k], tv[k][0]*u1+tv[k][1]*u2+tv[k][2]*u3),\
							# "  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f"%( projs[i][k][0],projs[i][k][1],projs[i][k][2],out[k][0],out[k][1],out[k][2])
					if(i == 0):
						outp.append(deepcopy(out))

			# for k in range(n):
			for k in subset:
				avg_diff_per_image[k] /= (sc*(sc-1)/2.0)
		elif( symmetry_class.sym[0] == "d" ):
			outp = deepcopy(projs)
			mirror_and_reduce_dsym(outp, subset, symmetry_class)

			for i in subset:
				qt = 0.0
				for k in range(sc-1):
					for l in range(k+1,sc):
						neisym = symmetry_class.symmetry_neighbors([outp[l][i][:3]])
						dmin = 180.0
						for q in neisym: dmin = min(dmin, getang3(q,outp[k][i]))
						qt += dmin
				avg_diff_per_image[i] = (qt/sc/(sc-1)/2.0)
				#k = subset[i]
				#lml = 1
				#print  "avg_diff_per_image  %3d  %8.2f="%(k,avg_diff_per_image[k]),\
				#"  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f"%( projs2[0][i][0],projs2[0][i][1],projs2[0][i][2],projs2[lml][i][0],projs2[lml][i][1],projs2[lml][i][2])

		else:  # o, t, i
			from sp_pixel_error import angle_diff
			outp = deepcopy(projs)
			for l in range(1,sc):
				psi_diff = angle_diff( [outp[l][j][2] for j in subset], [outp[0][j][2] for j in subset] )
				#  adjust psi if necessary
				if(abs(psi_diff-180.0) <90.0):
					for j in range(n):
						# apply mirror
						outp[l][j][2] = (outp[l][j][2] + 180.0)%360.0


		if(len(subset) == minimal_subset_size):
			break
		
		#  Remove element whose avg_diff_per_image is larger than max_error, if none, break
		max_error = -1.0
		the_worst_proj = -1
		for i in subset:
			if(avg_diff_per_image[i] > target_threshold):
				if(avg_diff_per_image[i]>max_error):
					max_error = avg_diff_per_image[i]
					the_worst_proj = i
		if( the_worst_proj > -1):	subset.remove(the_worst_proj)
		else:  break
		#print  "the_worst_proj",the_worst_proj
	#  End of pruning loop

	return subset, avg_diff_per_image, outp


"""Multiline Comment1"""

#MULTILINEMULTILINEMULTILINE 1
	#MULTILINEMULTILINEMULTILINE 1

	#MULTILINEMULTILINEMULTILINE 1
	#MULTILINEMULTILINEMULTILINE 1
	#MULTILINEMULTILINEMULTILINE 1

	#MULTILINEMULTILINEMULTILINE 1
		#MULTILINEMULTILINEMULTILINE 1
		#MULTILINEMULTILINEMULTILINE 1
		#MULTILINEMULTILINEMULTILINE 1
			#MULTILINEMULTILINEMULTILINE 1

	#MULTILINEMULTILINEMULTILINE 1
#MULTILINEMULTILINEMULTILINE 1


# parameters: list of (all) projections | reference volume | ...
#  Genetic programming version
#  The data structure:
#  [[L2, [parameters row-wise]], [], []...number_of_runs ]
#  It is kept on main proc
def ali3d_multishc(stack, ref_vol, ali3d_options, symmetry_class, mpi_comm = None, log = None, number_of_runs=2 ):

	from sp_alignment    import Numrinit, prepare_refrings, proj_ali_incore_local, shc
	from sp_utilities    import model_circle, get_input_from_string, get_params_proj, set_params_proj
	from sp_utilities    import write_text_row
	from sp_utilities    import wrap_mpi_gatherv, wrap_mpi_bcast, wrap_mpi_send, wrap_mpi_recv, wrap_mpi_split
	from mpi          import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD
	from mpi          import mpi_barrier, mpi_comm_split, mpi_comm_free, mpi_finalize
	from sp_projection   import prep_vol
	from sp_statistics   import hist_list
	from sp_applications import MPI_start_end
	from sp_filter       import filt_ctf
	from sp_global_def   import Util, ERROR
	from time         import time
	from random       import shuffle, random
	from copy         import deepcopy


	ir     = ali3d_options.ir
	rs     = ali3d_options.rs
	ou     = ali3d_options.ou
	xr     = ali3d_options.xr
	yr     = ali3d_options.yr
	ts     = ali3d_options.ts
	an     = ali3d_options.an
	delta  = ali3d_options.delta
	doga   = ali3d_options.doga
	center = ali3d_options.center
	CTF    = ali3d_options.CTF
	ref_a  = ali3d_options.ref_a
	L2threshold = ali3d_options.L2threshold

	# Optionally restrict out-of-plane angle
	if hasattr(ali3d_options, 'theta1'): theta1 = ali3d_options.theta1
	else: theta1 = -1.0
	
	if hasattr(ali3d_options, 'theta2'): theta2 = ali3d_options.theta2
	else: theta2 = -1.0
	
	if hasattr(ali3d_options, 'method'): method = ali3d_options.method
	else: method = "S"
	
	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	if log == None:
		from sp_logger import Logger
		log = Logger()

	number_of_proc = mpi_comm_size(mpi_comm)
	myid           = mpi_comm_rank(mpi_comm)
	main_node = 0

	if myid == main_node:
		log.add("Start VIPER1")

	if number_of_proc < number_of_runs:
		ERROR("number_of_proc < number_of_runs","VIPER1",1,myid)

	# if an != "-1":
	# 	ERROR("Option an not used","VIPER1",1,myid)

	# mpi_subcomm = mpi_comm_split(mpi_comm, myid % number_of_runs, myid / number_of_runs)
	# mpi_subrank = mpi_comm_rank(mpi_subcomm)
	# mpi_subsize = mpi_comm_size(mpi_subcomm)
	# mpi_subroots = range(number_of_runs)

	mpi_subcomm = wrap_mpi_split(mpi_comm, number_of_runs)
	mpi_subrank = mpi_comm_rank(mpi_subcomm)
	mpi_subsize = mpi_comm_size(mpi_subcomm)
	# do not make any assumptions about the subroots, collect the rank_id as they are already assigned
	if mpi_subrank == 0:
		mpi_subroots = wrap_mpi_gatherv([myid], 0, mpi_comm)
	else:
		mpi_subroots = wrap_mpi_gatherv([], 0, mpi_comm)
	mpi_subroots = wrap_mpi_bcast(mpi_subroots, main_node, mpi_comm)


	xrng        = get_input_from_string(xr)
	if  yr == "-1":  yrng = xrng
	else          :  yrng = get_input_from_string(yr)
	step        = get_input_from_string(ts)
	delta       = get_input_from_string(delta)
	lstp = min(len(xrng), len(yrng), len(step), len(delta))
	"""Multiline Comment2"""
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
		#MULTILINEMULTILINEMULTILINE 2
	#MULTILINEMULTILINEMULTILINE 2
	first_ring  = int(ir)
	rstep       = int(rs)
	last_ring   = int(ou)
	max_iter    = int(ali3d_options.maxit1)
	center      = int(center)

	vol = ref_vol
	nx      = vol.get_xsize()
	if last_ring < 0:
		last_ring = int(nx/2) - 2

	cnx = nx//2 + 1
	cny = cnx
	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	mask2D  = model_circle(last_ring,nx,nx) - model_circle(first_ring,nx,nx)

	if myid == main_node:
		list_of_particles = list(range(len(stack)))
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
	for im in range(nima):
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
		GA = [ [0.0, [[0.0,0.,0.0,0.0,0.0] for j in range(total_nima)]] for i in range(number_of_runs)]
		noimprovement = 0
		firstcheck    = True

	orient_and_shuffle = False
	# do the projection matching
	for N_step in range(lstp):  # At this point there is just one value here, it cannot loop.
		if myid == 0:  afterGAcounter = 0
		terminate = False
		Iter = 0
		while not terminate:

			Iter += 1
			total_iter += 1

			mpi_barrier(mpi_comm)
			if myid == main_node:
				log.add("ITERATION #%3d"%(total_iter))
				start_time = time()

			#=========================================================================
			# build references
			volft, kb = prep_vol(vol)
			#  We generate mirrored versions as well MAJOR CHANGE PAP 04/20/2017
			reference_angles = symmetry_class.even_angles(delta[N_step], theta1=theta1, theta2=theta2, method=method)
			refrings = prepare_refrings(volft, kb, nx, -1.0, reference_angles, "", numr, MPI=mpi_subcomm)
			del volft, kb
			#=========================================================================

			mpi_barrier(mpi_comm)
			if myid == main_node:
				log.add("Time to prepare rings: %f" % (time()-start_time))
				start_time = time()

			#=========================================================================
			#  It has to be here
			if orient_and_shuffle:
				# adjust params to references, calculate psi, calculate previousmax
				vecs = [[refrings[lr].get_attr("n1"), refrings[lr].get_attr("n2"), refrings[lr].get_attr("n3")] for lr in range(len(refrings))]
				for im in range(nima):
					#  For testing purposes make sure that directions did not change
					t1,t2,t3,t4,t5 = get_params_proj( data[im] )
					from sp_utilities import nearest_fang
					iqa = nearest_fang( vecs, t1, t2 ) # Here I could use more sophisticated distance for symmetries
					#if myid == 0 : 
					#print "  XXXX  ",myid,total_iter,im,iqa,t1,t2, reference_angles[iqa]
					#if total_iter>3:
					#	'''
					#	from mpi import mpi_finalize
					#	mpi_finalize()
					#	from sys import exit
					#	exit()
					#	'''

					#peak, temp = proj_ali_incore_local(data[im], [refrings[iqa]], [reference_angles[iqa]], numr, 0., 0., 1., 180.0 , sym="c1")
					cimage = Util.Polar2Dm(data[im], cnx, cny, numr, "F")
					Util.Normalize_ring(cimage, numr, 0 )
					Util.Frngs(cimage, numr)
					retvals = Util.Crosrng_e(refrings[iqa], cimage, numr, 0, 0.0)
					data[im].set_attr("previousmax", retvals["qn"])
					#u1,u2,t3,t4,t5 = get_params_proj( data[im])
					#if(abs(u1-t1)>1.0e-4 or  abs(u2-t2)>1.0e-4 ):  print "  PROBLEM IN  proj_ali_incore_local"
					#print  "peak1 ",myid,im,data[im].get_attr("ID"), retvals["qn"],[t1,t2,t3,t4,t5]
				if myid == main_node:
					log.add("Time to calculate first psi+previousmax: %f\n" % (time()-start_time))
					start_time = time()
				del vecs
			#=========================================================================

			mpi_barrier(mpi_comm)

			if myid == main_node: start_time = time()
			#=========================================================================
			# alignment
			if mpi_subrank == 0:
				pixer = []
				number_of_checked_refs = 0
				proj_ids_to_process = list(range(nima))
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
					peak, pixel_error, checked_refs, iref = shc(data[im], refrings, [[1.0,1.0]], numr, xrng[N_step], yrng[N_step], step[N_step], sym = "nomirror")
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
					for i in range(len(vector_assigned_refs)):
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
				for i in range(len(vector_previousmax)):
					data[proj_ids[i]].set_attr("previousmax", vector_previousmax[i])
					data[proj_ids[i]].set_attr("xform.projection", vector_xformprojs[i])
			#=========================================================================
			mpi_barrier(mpi_comm)
			if myid == main_node:
				log.add("Time of alignment = %f\n"%(time()-start_time))

			storevol = False

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
				log.add("= Pixel error        Number of images in all runs")
				for lhx in range(lhist):
					msg = " %10.3f                  %7d"%(region[lhx], histo[lhx])
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
				log.add("=================================================")
				log.add("Percent of positions with pixel error below 1.0 = ", (int(percent_of_pixerr_below_one*100)), "%","   Mutations: ",orient_and_shuffle)
			orient_and_shuffle = wrap_mpi_bcast(orient_and_shuffle, 0, mpi_comm)
			#=========================================================================

			#=========================================================================
			# centering, for d unnecessary, for cn, n>1 only z can move
			if center == -1 and symmetry_class.sym[0] == 'c':
				from sp_utilities      import estimate_3D_center_MPI, rotate_3D_shift
				cs[0], cs[1], cs[2], dummy, dummy = estimate_3D_center_MPI(data[image_start:image_end], total_nima, mpi_subrank, mpi_subsize, 0, mpi_comm=mpi_subcomm) #estimate_3D_center_MPI(data, number_of_runs*total_nima, myid, number_of_proc, main_node, mpi_comm=mpi_comm)
				if myid == main_node:
					msg = " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"%(cs[0], cs[1], cs[2])
					log.add(msg)
				if symmetry_class.nsym > 1 :
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
			"""Multiline Comment3"""
			#MULTILINEMULTILINEMULTILINE 3
			#MULTILINEMULTILINEMULTILINE 3
				#MULTILINEMULTILINEMULTILINE 3
				#MULTILINEMULTILINEMULTILINE 3
			#MULTILINEMULTILINEMULTILINE 3

			#MULTILINEMULTILINEMULTILINE 3
				#MULTILINEMULTILINEMULTILINE 3
				#MULTILINEMULTILINEMULTILINE 3
			#MULTILINEMULTILINEMULTILINE 3

			#=========================================================================
			if orient_and_shuffle:   #  DO orient
				params = []
				for im in data:
					phi, theta, psi, sx, sy = get_params_proj(im)
					params.append([phi, theta, psi, sx, sy])
				# if myid == 2:  print  " initial params before orient  ",myid,[get_params_proj(data[i]) for i in xrange(4)]

				# ------ orientation - begin
				#  Send solution from the main process of the first group to all processes in all groups
				params_0 = wrap_mpi_bcast(params, mpi_subroots[0], mpi_comm)
				if (mpi_subrank == 0) and (myid != 0):
					#  This is done on the main node of each group (master node for MPI_COMM_WORLD skips it)
					#  Minimal length of the subset is set to 1/3 of the number of parameters
					#  Error threshold is set somewhat arbitrarily to 1.5 angular step of reference projections
					#  params gets overwritten by rotated parameters,  subset is a list of indexes common
					subset, avg_diff_per_image, params = find_common_subset([params_0, params], delta[N_step]*1.5, len(params)/3, symmetry_class)
					params = params[1]
					#   neither is needed
					del subset, avg_diff_per_image
					# if myid == 2:  print  " params before orient  ",myid,params[:4],params[-4:]
					#write_text_row(params_0,"bparamszero%04d%04d.txt"%(myid,total_iter))
					#write_text_row(params,"bparams%04d%04d.txt"%(myid,total_iter))
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
				for i in range(nima):  set_params_proj(data[i], params[i])
				vol = do_volume(data[image_start:image_end], ali3d_options, 0, mpi_subcomm)
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
							all_params.append(deepcopy(params))
						else:
							all_L2s.append(wrap_mpi_recv(sr, mpi_comm))
							all_params.append(wrap_mpi_recv(sr, mpi_comm))
				else:
					if mpi_subrank == 0:
						wrap_mpi_send(L2, 0, mpi_comm)
						wrap_mpi_send(params, 0, mpi_comm)

				# ---------------------------------

				#  Add params to GA, sort, check termination and if not terminate do crossover and send back
				if myid == 0:
					#  after GA move do 3 iterations to give the program a chance to improve mutated structures.
					#all_params = shuffle_configurations(all_params)

					for i in range(number_of_runs):
						log.add("L2 incoming norm for volume %3d  = %f"%(i,all_L2s[i]))

					for i in range(number_of_runs):
						GA.append([all_L2s[i],deepcopy(all_params[i])])
					#  check whether this move will improve anything
					all_L2s.sort(reverse=True)
					#print " sorted terminate  ",all_L2s
					#for i in xrange(number_of_runs): print GA[i][0]

					noreseeding = True

					if(all_L2s[0]<GA[number_of_runs-1][0]):
						if firstcheck:  sp_global_def.sxprint("  SHOULD NOT BE HERE")
						noimprovement += 1
						if(noimprovement == 2):  terminate = True
						GA = GA[:number_of_runs]
						log.add("VIPER1 could not find better solutions, it will terminate now.")
					else:
						noimprovement = 0
						GA.sort(reverse=True)
						GA = GA[:number_of_runs]
						if( symmetry_class.sym[0] == "d"  and   GA[0][0]>0.0 ):
							for i in range(1,len(GA)):
								mirror_and_reduce_dsym([GA[0][1],GA[i][1]], list(range(len(GA[0][1]))), symmetry_class)

						#  ---  Stopping criterion
						from sp_statistics import table_stat
						from math import sqrt
						q1,q2,q3,q4 = table_stat([GA[i][0] for i in range(number_of_runs)])
						# Terminate if variation of L2 norms less than (L2threshold*100)% of their average
						crit = sqrt(max(q2,0.0))/q1
						for i in range(number_of_runs):
							log.add("L2 norm for volume %3d  = %f"%(i,GA[i][0]))
						log.add("L2 norm std dev %f\n"%crit)
						crit = crit < L2threshold
						if (Iter < max_iter) and (firstcheck and crit):
							noreseeding = False
							terminate   = False
							log.add("Insufficient initial variability, reseeding!\n")
							all_params = [[[random()*360.0,random()*180.0,random()*360.0,0.0,0.0]\
											 for j in range(total_nima)] for i in range(number_of_runs)]
						else:  	terminate = Iter > max_iter or crit

						firstcheck = False

					if not terminate and noimprovement == 0 and noreseeding:
						afterGAcounter = 3
						#  Now do the crossover
						all_params = []

						from sp_utilities import nearest_many_full_k_projangles, angles_to_normals
						from random import random, randint, shuffle
						# select random pairs of solutions
						ipl = list(range(number_of_runs))
						shuffle(ipl)
						for ip in range(0,2*(len(ipl)/2)+len(ipl)%2,2):
							#  random reference projection:
							itmp = randint(0,total_nima-1)
							#print  "  nearest_many_full_k_projangles  ",total_nima,itmp,ipl,ip,len(GA[ipl[ip]][1]),GA[ipl[ip]][1][itmp],GA[ipl[ip]][1]
							keepset = nearest_many_full_k_projangles(angles_to_normals(GA[ipl[ip]][1]), [GA[ipl[ip]][1][itmp]], howmany = total_nima/2, sym_class = symmetry_class)[0]
							#print  "  keepset  ",total_nima,len(keepset),itmp,keepset
							otherset = set(range(total_nima)) - set(keepset)
							otherset = [i for i in otherset]
							keepset.sort()
							otherset.sort()
							newparms1 = [None]*total_nima
							newparms2 = [None]*total_nima
							for i in keepset:
								newparms1[i] = deepcopy(GA[ipl[ip]][1][i])
								newparms2[i] = deepcopy(GA[ipl[(ip+1)%number_of_runs]][1][i])
							for i in otherset:
								newparms1[i] = deepcopy(GA[ipl[(ip+1)%number_of_runs]][1][i])
								newparms2[i] = deepcopy(GA[ipl[ip]][1][i])
							#print "  PRINTOUT SHUFFLED   ",ipl[ip],ipl[ip+1]
							"""Multiline Comment4"""
							#MULTILINEMULTILINEMULTILINE 4
								#MULTILINEMULTILINEMULTILINE 4
							#MULTILINEMULTILINEMULTILINE 4
								#MULTILINEMULTILINEMULTILINE 4
							#MULTILINEMULTILINEMULTILINE 4
								#MULTILINEMULTILINEMULTILINE 4
								#MULTILINEMULTILINEMULTILINE 4
							#MULTILINEMULTILINEMULTILINE 4

							#  Put mutated params on one list
							all_params.append(deepcopy(newparms1))
							all_params.append(deepcopy(newparms2))
						all_params = all_params[:number_of_runs]
						#  Try this 02/03/2015 PAP
						#  Always mutate the first ones
						#  for half of projections 'mirror' them by adding 180 to psi
						keepset = max(1,int(0.25*number_of_runs))
						#ipl = range(total_nima)
						#shuffle(ipl)
						#ipl = ipl[:total_nima//2]
						for i in range(keepset):
							all_params[0][i][2] += 180.0
							#  Always reseed the last ones
							all_params[-1-i] = [[random()*360.0,random()*180.0,random()*360.0,0.0,0.0]\
										 for j in range(total_nima)]

				terminate = wrap_mpi_bcast(terminate, main_node, mpi_comm)
				if not terminate:

					storevol = True

					# Send params back
					if myid == 0:
						for i in range(number_of_runs):
							sr = mpi_subroots[i]
							if sr == myid:
								params = all_params[i]
							else:
								wrap_mpi_send(all_params[i], sr, mpi_comm)
					else:
						if mpi_subrank == 0:
							params = wrap_mpi_recv(0, mpi_comm)

					params = wrap_mpi_bcast(params, 0, mpi_subcomm)
					for i in range(nima):
						set_params_proj(data[i], params[i])
					"""Multiline Comment5"""
					#MULTILINEMULTILINEMULTILINE 5
					#MULTILINEMULTILINEMULTILINE 5
					#MULTILINEMULTILINEMULTILINE 5
					#MULTILINEMULTILINEMULTILINE 5
						#MULTILINEMULTILINEMULTILINE 5
					#MULTILINEMULTILINEMULTILINE 5
					#MULTILINEMULTILINEMULTILINE 5
					#MULTILINEMULTILINEMULTILINE 5
						#MULTILINEMULTILINEMULTILINE 5
						#MULTILINEMULTILINEMULTILINE 5
						#MULTILINEMULTILINEMULTILINE 5

					#MULTILINEMULTILINEMULTILINE 5
					#MULTILINEMULTILINEMULTILINE 5
					#MULTILINEMULTILINEMULTILINE 5
					#MULTILINEMULTILINEMULTILINE 5
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
				vol = do_volume(data[image_start:image_end], ali3d_options, 0, mpi_subcomm)

				if len(ali3d_options.moon_elimination) > 0:
					from sp_utilities import eliminate_moons
					vol = eliminate_moons(vol, ali3d_options.moon_elimination)

				if mpi_subrank == 0:
					L2 = vol.cmp("dot", vol, dict(negative = 0, mask = model_circle(last_ring, nx, nx, nx)))
					# if myid == 2:  print  " Right after reconstruction L2", myid, L2,[get_params_proj(data[i]) for i in xrange(4)]
					#print  " Right after reconstruction L2", myid, total_iter,L2
					#if storevol:   vol.write_image("mutated%04d%04d.hdf"%(myid,total_iter))

				# log
				if myid == main_node:
					log.add("3D reconstruction time = %f\n"%(time()-start_time))
					start_time = time()

			"""Multiline Comment6"""
			#MULTILINEMULTILINEMULTILINE 6
			#MULTILINEMULTILINEMULTILINE 6
				#MULTILINEMULTILINEMULTILINE 6
				#MULTILINEMULTILINEMULTILINE 6
				#MULTILINEMULTILINEMULTILINE 6
					#MULTILINEMULTILINEMULTILINE 6
						#MULTILINEMULTILINEMULTILINE 6
					#MULTILINEMULTILINEMULTILINE 6
					#MULTILINEMULTILINEMULTILINE 6
					#MULTILINEMULTILINEMULTILINE 6
				#MULTILINEMULTILINEMULTILINE 6
			#MULTILINEMULTILINEMULTILINE 6
			"""Multiline Comment7"""
			#MULTILINEMULTILINEMULTILINE 7
			#MULTILINEMULTILINEMULTILINE 7
				#MULTILINEMULTILINEMULTILINE 7
				#MULTILINEMULTILINEMULTILINE 7
				#MULTILINEMULTILINEMULTILINE 7
					#MULTILINEMULTILINEMULTILINE 7
					#MULTILINEMULTILINEMULTILINE 7
						#MULTILINEMULTILINEMULTILINE 7
					#MULTILINEMULTILINEMULTILINE 7
						#MULTILINEMULTILINEMULTILINE 7
			#MULTILINEMULTILINEMULTILINE 7
				#MULTILINEMULTILINEMULTILINE 7
					#MULTILINEMULTILINEMULTILINE 7

			#MULTILINEMULTILINEMULTILINE 7
			#MULTILINEMULTILINEMULTILINE 7
				#MULTILINEMULTILINEMULTILINE 7
				#MULTILINEMULTILINEMULTILINE 7
			#MULTILINEMULTILINEMULTILINE 7
				#MULTILINEMULTILINEMULTILINE 7

			#MULTILINEMULTILINEMULTILINE 7
			#MULTILINEMULTILINEMULTILINE 7
			#MULTILINEMULTILINEMULTILINE 7
			#MULTILINEMULTILINEMULTILINE 7
				#MULTILINEMULTILINEMULTILINE 7
			#MULTILINEMULTILINEMULTILINE 7
			#MULTILINEMULTILINEMULTILINE 7
				#MULTILINEMULTILINEMULTILINE 7
				#MULTILINEMULTILINEMULTILINE 7
				#MULTILINEMULTILINEMULTILINE 7

			#MULTILINEMULTILINEMULTILINE 7
			#MULTILINEMULTILINEMULTILINE 7
			#MULTILINEMULTILINEMULTILINE 7

			#MULTILINEMULTILINEMULTILINE 7
			#MULTILINEMULTILINEMULTILINE 7
			#MULTILINEMULTILINEMULTILINE 7
			#MULTILINEMULTILINEMULTILINE 7
				#MULTILINEMULTILINEMULTILINE 7
				#MULTILINEMULTILINEMULTILINE 7
			#MULTILINEMULTILINEMULTILINE 7


	#if mpi_subrank == 0:
	#	for im in xrange(len(data)):
	#		print  "VIPER1 peak ",myid,im,data[im].get_attr("ID"), data[im].get_attr("previousmax"), get_params_proj(data[im])
	#=========================================================================
	mpi_comm_free(mpi_subcomm)
	
	
	if myid == main_node:
		log.add("Finished VIPER1")
		return GA[0][1]
	else:
		return None  # results for the other processes


# parameters: list of (all) projections | reference volume | ...
def ali3d_multishc_2(stack, ref_vol, ali3d_options, symmetry_class, mpi_comm = None, log = None ):

	from sp_alignment       import Numrinit, prepare_refrings, proj_ali_incore_local, shc
	from sp_utilities       import model_circle, get_input_from_string, get_params_proj, wrap_mpi_gatherv, wrap_mpi_bcast, wrap_mpi_split
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier
	from sp_projection      import prep_vol
	from sp_statistics      import hist_list
	from sp_applications    import MPI_start_end
	from sp_filter          import filt_ctf
	from sp_fundamentals    import symclass
	from sp_global_def import Util
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

	# Optionally restrict out-of-plane angle
	if hasattr(ali3d_options, 'theta1'): theta1 = ali3d_options.theta1
	else: theta1 = -1.0
	
	if hasattr(ali3d_options, 'theta2'): theta2 = ali3d_options.theta2
	else: theta2 = -1.0
	
	if hasattr(ali3d_options, 'method'): method = ali3d_options.method
	else: method = "S"
	
	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	if log == None:
		from sp_logger import Logger
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

	symmetry_class = symclass(sym)

	# if an != "-1":
	# 	ERROR("Option an not used","VIPER1",1,myid)
	"""Multiline Comment8"""
	#MULTILINEMULTILINEMULTILINE 8
		#MULTILINEMULTILINEMULTILINE 8
	#MULTILINEMULTILINEMULTILINE 8
		#MULTILINEMULTILINEMULTILINE 8
	#MULTILINEMULTILINEMULTILINE 8

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
		list_of_particles = list(range(len(stack)))
		total_nima = len(list_of_particles)
	else:
		list_of_particles = None
		total_nima = None
	total_nima = wrap_mpi_bcast(total_nima, main_node, mpi_comm)
	list_of_particles = wrap_mpi_bcast(list_of_particles, main_node, mpi_comm)

	#old_mpi_comm = mpi_comm
	#mpi_size = mpi_comm_size(mpi_comm)

	## if there are fewer images than processors then split processors in 2 groups
	## one in which each processor analyzes one image, and another in which
	## processors stay idle and wait for the other group to finish
	#if (mpi_size > total_nima):
		#if (myid < total_nima):
			#mpi_subcomm = mpi_comm_split(mpi_comm, 0, myid)
			#mpi_comm = mpi_subcomm
		#else:
			#mpi_subcomm = mpi_comm_split(mpi_comm, 1, myid - total_nima)
			#mpi_barrier(mpi_comm)
			#return None, None, None, None


	number_of_proc = mpi_comm_size(mpi_comm)
	myid           = mpi_comm_rank(mpi_comm)


	image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
	# create a list of images for each node
	list_of_particles = list_of_particles[image_start: image_end]
	nima = len(list_of_particles)

	data = [ stack[im] for im in list_of_particles ]
	for im in range(nima):
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
	#qvol = volume_reconstruction(data, ali3d_options, mpi_comm)
	qvol = do_volume(data, ali3d_options, 0, mpi_comm)
	# log

	if myid == main_node:
		start_time = time()
		#qvol.write_image("vitera.hdf")
		L2 = qvol.cmp("dot", qvol, dict(negative = 0, mask = model_circle(last_ring, nx, nx, nx)))
		log.add("3D reconstruction time = %f\n"%(time()-start_time)," START  L2 norm:  %f"%L2)
		start_time = time()
	del qvol


	from sp_projection   import prep_vol, prgs
	from sp_alignment import ringwe
	cnx = nx//2 + 1
	cny = nx//2 + 1
	wr_four  = ringwe(numr, "F")
	from math import pi, sin, cos
	qv = pi/180.
	volft, kb = prep_vol(ref_vol)
	from sp_utilities import get_params_proj
	for im in range(nima):
		phi,theta,psi,tx,ty = get_params_proj(data[im])
		refrings = prgs(volft, kb, [phi,theta,psi, 0.0, 0.0])
		refrings = Util.Polar2Dm(refrings, cnx, cny, numr, "F")
		Util.Normalize_ring(refrings, numr, 0 )
		Util.Frngs(refrings, numr)
		Util.Applyws(refrings, numr, wr_four)
		"""Multiline Comment9"""
		#MULTILINEMULTILINEMULTILINE 9
		#MULTILINEMULTILINEMULTILINE 9
		#MULTILINEMULTILINEMULTILINE 9
		#MULTILINEMULTILINEMULTILINE 9
		#MULTILINEMULTILINEMULTILINE 9
		#MULTILINEMULTILINEMULTILINE 9
		#MULTILINEMULTILINEMULTILINE 9
		#MULTILINEMULTILINEMULTILINE 9
		#MULTILINEMULTILINEMULTILINE 9

		#print "orin  ",data[im].get_attr("ID"),get_params_proj(data[im])
		#peak, pixer = proj_ali_incore_local(data[im],refrings, [[phi, theta], [phi, theta]], numr,0.0,0.0,1.0,delta[0]/4)
		cimage = Util.Polar2Dm(data[im], cnx, cny, numr, "F")
		Util.Normalize_ring(cimage, numr, 0 )
		Util.Frngs(cimage, numr)
		retvals = Util.Crosrng_e(refrings, cimage, numr, 0, 0.0)
		data[im].set_attr("previousmax", retvals["qn"])
		#print  " VIPER2  peak ",myid,im,data[im].get_attr("ID"), retvals["qn"],get_params_proj(data[im])

	"""Multiline Comment10"""
	#MULTILINEMULTILINEMULTILINE 10
	#MULTILINEMULTILINEMULTILINE 10
	#MULTILINEMULTILINEMULTILINE 10
	#MULTILINEMULTILINEMULTILINE 10
	#MULTILINEMULTILINEMULTILINE 10
	#MULTILINEMULTILINEMULTILINE 10
	#MULTILINEMULTILINEMULTILINE 10
	#MULTILINEMULTILINEMULTILINE 10
	#MULTILINEMULTILINEMULTILINE 10
	#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10

		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
		#MULTILINEMULTILINEMULTILINE 10
	#MULTILINEMULTILINEMULTILINE 10
	del volft
	# volume reconstruction
	mpi_barrier(mpi_comm)
	if myid == main_node:
		log.add("Time to calculate first psi+shifts+previousmax: %f\n" % (time()-start_time))
		start_time = time()
	ref_vol = do_volume(data, ali3d_options, 0, mpi_comm)

	if myid == main_node:
		#ref_vol.write_image("viterb.hdf")
		L2 = ref_vol.cmp("dot", ref_vol, dict(negative = 0, mask = model_circle(last_ring, nx, nx, nx)))
		log.add("3D reconstruction time = %f\n"%(time()-start_time),"   L2 norm:  %f"%L2)
		start_time = time()

	#=========================================================================



	pixer = [0.0]*nima
	par_r = [[] for im in list_of_particles ]
	cs = [0.0]*3
	total_iter = 0
	# do the projection matching
	for N_step in range(lstp):
		
		terminate = 0
		Iter = 0
		while Iter < max_iter and terminate == 0:

			Iter += 1
			total_iter += 1

			mpi_barrier(mpi_comm)
			if myid == main_node:
				log.add("ITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, xrange = %5.2f, yrange = %5.2f, step = %5.2f\n"%(total_iter, Iter, delta[N_step], xrng[N_step],yrng[N_step],step[N_step]))
				start_time = time()

			#=========================================================================
			# build references
			volft, kb = prep_vol(vol)
			#  For the local SHC it is essential reference projections have psi zero, as otherwise it will get messed up.
			reference_angles = symmetry_class.even_angles(delta[N_step], phiEqpsi = "Zero", theta1=theta1, theta2=theta2, method=method)
			refrings = prepare_refrings(volft, kb, nx, -1.0, reference_angles, "", numr, MPI=mpi_comm)
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
			for im in range(nima):
				#peak, pixer[im], checked_refs, number_of_peaks = shc_multi(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step], number_of_runs=number_of_runs)
				# previousmax is set in shc
				peak, pixer[im], checked_refs, iref = shc(data[im], refrings, [[1.0,1.0]], numr, xrng[N_step], yrng[N_step], step[N_step], sym = "nomirror") # cannot use 'an' here
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
				for lhx in range(lhist):
					msg = " %10.3f     %7d"%(region[lhx], histo[lhx])
					log.add(msg)
				if (max(all_pixer) < 0.5) and (sum(all_pixer)/total_nima < 0.05):
					terminate = 1
			terminate = wrap_mpi_bcast(terminate, main_node, mpi_comm)
			#=========================================================================

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
			# vol = volume_reconstruction(data, ali3d_options, mpi_comm)
			vol = do_volume(data, ali3d_options, 0, mpi_comm)
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

	## if there are fewer images than processors then synchronize
	## with the other group of processors that did not do any work
	#if (mpi_size > total_nima):
		#if (myid < no_of_images):
			#mpi_comm = old_mpi_comm
			#mpi_barrier(mpi_comm)

	if myid == main_node: 
		log.add("Finished VIPER2")
		return params, vol, previousmax, par_r
	else:
		return None, None, None, None  # results for the other processes

# data - projections (scattered between cpus)
# options - the same for all cpus
# return - volume the same for all cpus
def volume_reconstruction(data, options, mpi_comm):
	from mpi import mpi_comm_rank
	from sp_reconstruction import recons3d_4nn_MPI, recons3d_4nn_ctf_MPI
	from sp_utilities import bcast_EMData_to_all, model_circle
	
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
	else:   vol = recons3d_4nn_MPI    (myid, data,      symmetry=sym, snr=snr, npad=npad, mpi_comm=mpi_comm)

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

	from sp_reconstruction import recons3d_4nn, recons3d_4nn_ctf
	from sp_utilities import bcast_EMData_to_all, model_circle
	
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


"""Multiline Comment11"""
#MULTILINEMULTILINEMULTILINE 11
#MULTILINEMULTILINEMULTILINE 11
  #MULTILINEMULTILINEMULTILINE 11
  #MULTILINEMULTILINEMULTILINE 11
  #MULTILINEMULTILINEMULTILINE 11
  #MULTILINEMULTILINEMULTILINE 11
  	#MULTILINEMULTILINEMULTILINE 11
  	#MULTILINEMULTILINEMULTILINE 11

#MULTILINEMULTILINEMULTILINE 11

# all_projs and subset must be set only for root (MPI rank == 0)
# remaining parameters must be set for all
# size of mpi_communicator must be >= runs_count
def multi_shc(all_projs, subset, runs_count, ali3d_options, mpi_comm, log=None, ref_vol=None):
	"""
	Arguments:
		all_projs : Stack of projections
		subset : Selection of images
		runs_count : Number of quasi-independent shc runs
		ali3d_options : Command-line options (as a namespace)
			.sym : Symmetry
			.theta1 : Starting tilt angle (optional)
			.theta2 : Ending tilt angle (optional)
			.ref_a : Method for generating quasi-uniformly distributed projection directions
			.delta : Angular increment
			.ou : Outer radius
			.dpi : Resolution of BILD angular-distribution file (optional)
		mpi_comm : MPI communicator
		log : Log file instance (optional)
		ref_vol : Reference volume (optional)
	"""
	
	from sp_applications import MPI_start_end
	from sp_utilities import set_params_proj, wrap_mpi_bcast, write_text_row, drop_image, write_text_file
	from sp_utilities import bcast_EMData_to_all, bcast_list_to_all, bcast_number_to_all
	from sp_fundamentals import symclass
	from sp_global_def import ERROR
	import sp_global_def
	from random import random
	from mpi import mpi_comm_rank, mpi_comm_size, mpi_finalize, mpi_comm_split, mpi_barrier
	from argparse import Namespace
	from sp_utilities import angular_distribution, get_im
	import os
	
	mpi_rank = mpi_comm_rank(mpi_comm)
	mpi_size = mpi_comm_size(mpi_comm)

	global BATCH, MPI
	sp_global_def.BATCH = True
	sp_global_def.MPI   = True
	if(mpi_size < runs_count):  ERROR("multi_shc","mpi_size < runs_count",1,mpi_rank)

	if log == None:
		from sp_logger import Logger
		log = Logger()


	#  Initialize symmetries
	symmetry_class = symclass(ali3d_options.sym)

	error = 0
	projections = []
	if mpi_rank == 0:
		# Optionally restrict out-of-plane angle
		if hasattr(ali3d_options, 'theta1'): 
			theta1 = ali3d_options.theta1
		else: 
			theta1 = -1.0
		
		if hasattr(ali3d_options, 'theta2'):
			theta2 = ali3d_options.theta2
		else:
			theta2 = -1.0
		
		if hasattr(ali3d_options, 'method'):
			method = ali3d_options.ref_a
		else:
			method = "S"

		if not hasattr(ali3d_options, 'filament_width'):
			ali3d_options.filament_width = -1
		
		prms = symmetry_class.even_angles(float(ali3d_options.delta), theta1=theta1, theta2=theta2, method=method)
		write_text_row(prms, log.prefix + "initangles.txt")
		
		if(len(prms) < len(subset)): error = 1
		else:
			from random import shuffle
			shuffle(prms)
			for i in subset:
				prms[i][2] = random()*360.0
				set_params_proj(all_projs[i], prms[i]+[0.,0.])
				#all_projs[i].set_attr("stable", 0)
				all_projs[i].set_attr("previousmax", -1.e23)
				projections.append(all_projs[i])
		del prms

	error = bcast_number_to_all(error, source_node = 0)
	if(error == 1): ERROR("Angular step too large, decrease delta", 1, mpi_rank)

	###from sys import exit
	#if mpi_rank == 0:   print "  NEW   ",mpi_rank
	###exit()
	projections = wrap_mpi_bcast(projections, 0, mpi_comm)

	n_projs = len(projections)

	if ref_vol == None:
		if (mpi_size > n_projs):
			working = int(not(mpi_rank < n_projs))
			mpi_subcomm = mpi_comm_split(mpi_comm, working,  mpi_rank - working*n_projs)
			mpi_subsize = mpi_comm_size(mpi_subcomm)
			mpi_subrank = mpi_comm_rank(mpi_subcomm)
			if (mpi_rank < n_projs):
				proj_begin, proj_end = MPI_start_end(n_projs, mpi_subsize, mpi_subrank)
				ref_vol = do_volume(projections[proj_begin:proj_end], ali3d_options, 0, mpi_comm=mpi_subcomm)
			else:
				from sp_utilities import model_blank
				nx = projections[0].get_xsize()
				ref_vol = model_blank(nx,nx,nx)
			bcast_EMData_to_all(ref_vol, mpi_rank, 0, comm=mpi_comm)
		else:
			proj_begin, proj_end = MPI_start_end(n_projs, mpi_size, mpi_rank)
			ref_vol = do_volume(projections[proj_begin:proj_end], ali3d_options, 0, mpi_comm=mpi_comm)

	# Each node keeps all projection data, this would not work for large datasets
	out_params = ali3d_multishc(projections, ref_vol, ali3d_options, symmetry_class, mpi_comm=mpi_comm, log=log, number_of_runs=runs_count)
	"""Multiline Comment12"""
	#MULTILINEMULTILINEMULTILINE 12
		#MULTILINEMULTILINEMULTILINE 12
	#MULTILINEMULTILINEMULTILINE 12

	if mpi_rank == 0:
		"""Multiline Comment13"""
		#MULTILINEMULTILINEMULTILINE 13
		#MULTILINEMULTILINEMULTILINE 13
			#MULTILINEMULTILINEMULTILINE 13
			#MULTILINEMULTILINEMULTILINE 13
			#MULTILINEMULTILINEMULTILINE 13

		#MULTILINEMULTILINEMULTILINE 13
		#MULTILINEMULTILINEMULTILINE 13
			#MULTILINEMULTILINEMULTILINE 13
			#MULTILINEMULTILINEMULTILINE 13
				#MULTILINEMULTILINEMULTILINE 13
				#MULTILINEMULTILINEMULTILINE 13
				#MULTILINEMULTILINEMULTILINE 13
			#MULTILINEMULTILINEMULTILINE 13
			#MULTILINEMULTILINEMULTILINE 13
		#MULTILINEMULTILINEMULTILINE 13
		#MULTILINEMULTILINEMULTILINE 13
			#MULTILINEMULTILINEMULTILINE 13
			#MULTILINEMULTILINEMULTILINE 13
		#MULTILINEMULTILINEMULTILINE 13
		temp = []
		from sp_utilities import get_params_proj
		for i in range(n_projs):
			set_params_proj( projections[i], out_params[i] )
		write_text_row(out_params, log.prefix + "refparams2.txt")
		"""Multiline Comment14"""
		#MULTILINEMULTILINEMULTILINE 14
		#MULTILINEMULTILINEMULTILINE 14
		#MULTILINEMULTILINEMULTILINE 14
		#MULTILINEMULTILINEMULTILINE 14
		#MULTILINEMULTILINEMULTILINE 14
		#MULTILINEMULTILINEMULTILINE 14
			#MULTILINEMULTILINEMULTILINE 14
			#MULTILINEMULTILINEMULTILINE 14
		#MULTILINEMULTILINEMULTILINE 14
	"""Multiline Comment15"""
	#MULTILINEMULTILINEMULTILINE 15
		#MULTILINEMULTILINEMULTILINE 15
	#MULTILINEMULTILINEMULTILINE 15
	# proj_begin, proj_end  = MPI_start_end(n_projs, mpi_size, mpi_rank)

	projections = wrap_mpi_bcast(projections, 0, mpi_comm)
	from sp_utilities import get_params_proj

	if (mpi_size > n_projs):
		working = int(not(mpi_rank < n_projs))
		mpi_subcomm = mpi_comm_split(mpi_comm, working,  mpi_rank - working*n_projs)
		mpi_subsize = mpi_comm_size(mpi_subcomm)
		mpi_subrank = mpi_comm_rank(mpi_subcomm)
		if (mpi_rank < n_projs):
			proj_begin, proj_end = MPI_start_end(n_projs, mpi_subsize, mpi_subrank)
			ref_vol = do_volume(projections[proj_begin:proj_end], ali3d_options, 0, mpi_comm=mpi_subcomm)
		else:
			from sp_utilities import model_blank
			nx = projections[0].get_xsize()
			ref_vol = model_blank(nx,nx,nx)
		bcast_EMData_to_all(ref_vol, mpi_rank, 0, comm=mpi_comm)
	else:
		proj_begin, proj_end = MPI_start_end(n_projs, mpi_size, mpi_rank)
		ref_vol = do_volume(projections[proj_begin:proj_end], ali3d_options, 0, mpi_comm=mpi_comm)

	if mpi_rank == 0:
		ref_vol.write_image(log.prefix + "refvol2.hdf")
		from sp_utilities import model_circle
		nx = ref_vol.get_xsize()
		L2 = ref_vol.cmp("dot", ref_vol, dict(negative = 0, mask = model_circle(ali3d_options.ou, nx,nx,nx)))
		log.add(" L2 norm of reference volume:  %f"%L2)

		# Generate angular distribution
		if hasattr(ali3d_options, 'dpi'):
			dpi = ali3d_options.dpi
		else:
			dpi = 72
		pixel_size = 1  # Not going to upscale to the original dimensions, so in Chimera open reconstruction at 1 Angstrom/voxel, etc.
		
		angular_distribution(
			params_file=log.prefix+'refparams2.txt',
			output_folder=log.prefix,
			prefix='refvol2_angdist',
			method=method,
			pixel_size=pixel_size,
			delta=float(ali3d_options.delta),
			symmetry=ali3d_options.sym,
			box_size=nx,
			particle_radius=ali3d_options.ou,
			dpi=dpi,
			do_print=False,
			)
		
	"""Multiline Comment16"""
	#MULTILINEMULTILINEMULTILINE 16
		#MULTILINEMULTILINEMULTILINE 16
		#MULTILINEMULTILINEMULTILINE 16
		#MULTILINEMULTILINEMULTILINE 16
			#MULTILINEMULTILINEMULTILINE 16
			#MULTILINEMULTILINEMULTILINE 16
			#MULTILINEMULTILINEMULTILINE 16
			#MULTILINEMULTILINEMULTILINE 16
			#MULTILINEMULTILINEMULTILINE 16
		#MULTILINEMULTILINEMULTILINE 16
	#MULTILINEMULTILINEMULTILINE 16



	if (mpi_size > n_projs):
		working = int(not(mpi_rank < n_projs))
		mpi_subcomm = mpi_comm_split(mpi_comm, working,  mpi_rank - working*n_projs)
		mpi_subsize = mpi_comm_size(mpi_subcomm)
		mpi_subrank = mpi_comm_rank(mpi_subcomm)
		if (mpi_rank < n_projs):
			out_params, out_vol, previousmax, out_r = ali3d_multishc_2(projections, ref_vol, ali3d_options, symmetry_class, mpi_comm=mpi_subcomm, log=log)
		else:
			out_params = None
			out_vol    = None
		mpi_barrier(mpi_comm)
	else:
		out_params, out_vol, previousmax, out_r = ali3d_multishc_2(projections, ref_vol, ali3d_options, symmetry_class, mpi_comm=mpi_comm, log=log)


	if mpi_rank == 0:
		write_text_file(previousmax, log.prefix + "previousmax.txt")
		write_text_row(out_params, log.prefix + "params.txt")
		drop_image(out_vol, log.prefix + "volf.hdf")
		
		# Generate angular distribution
		independent_run_dir = log.prefix
		sp_global_def.sxprint('independent_run_dir', independent_run_dir)
		
		params_file = log.prefix + "params.txt"
		output_folder = independent_run_dir 
		prefix = 'volf_angdist'  # will overwrite input parameters file if blank
		delta = float(ali3d_options.delta)
		symmetry = ali3d_options.sym
		if hasattr(ali3d_options, 'dpi'):
			dpi = ali3d_options.dpi
		else:
			dpi = 72
		
		# Not going to upscale to the original dimensions, so in Chimera open reconstruction at 1 Angstrom/voxel, etc.
		pixel_size = 1
		box_size = get_im( os.path.join(log.prefix, 'volf.hdf') ).get_xsize()
		
		angular_distribution(
			params_file=params_file,
			output_folder=output_folder,
			prefix=prefix,
			method=method,
			pixel_size=pixel_size,
			delta=delta,
			symmetry=symmetry,
			box_size=box_size,
			particle_radius=ali3d_options.ou,
			dpi=dpi,
			do_print=False,
			)
		
	return out_params, out_vol, None#, out_peaks

def mirror_and_reduce_dsym(params, indexes, symmetry_class):
	#  Input params contains multiple datasets [ p0, p1, ...]
	#  We treat p0 as a reference
	# For D symmetry there are two equivalent positions that agree with given Dn symmetry
	#  The second is rotated by 360/n degrees.

	from sp_utilities import getang3
	from sp_pixel_error import angle_diff

	sc = len(params)
	ns = len(params[0])

	#  bbdb is 360.0/nsym and indicates position of the second symmetry
	bbdb = 360.0/symmetry_class.nsym
	symphi = 360.0/symmetry_class.nsym*2
	#  For each set we have four positions to consider: straight, straight psi mirrored, phi+bdb, phi+bdb and psi mirrored
	for i in range(1,sc):
		psi_diff = angle_diff( [params[i][j][2] for j in indexes], [params[0][j][2] for j in indexes] )
		#  adjust psi if necessary
		if(abs(psi_diff-180.0) <90.0):
			for j in range(ns):
				# apply mirror
				params[i][j][2] = (params[i][j][2] + 180.0)%360.0
		# Check which one of the two possible symmetry positions is closer
		#  compute angular errors including symmetry
		per1 = 0.0
		per2 = 0.0
		temp = [None]*ns
		for j in indexes:
			neisym = symmetry_class.symmetry_neighbors([params[i][j][:3]])
			dmin = 180.0
			for q in neisym: dmin = min(dmin, getang3(q,params[0][j]))
			per1 += dmin
			temp = symmetry_class.reduce_anglesets([params[i][j][0]+bbdb,params[i][j][1],params[i][j][2]] )
			neisym = symmetry_class.symmetry_neighbors([temp])
			dmin = 180.0
			for q in neisym: dmin = min(dmin, getang3(q,params[0][j]))
			per2 += dmin

		if(per2<per1):
			for j in range(ns):
				temp = symmetry_class.reduce_anglesets([params[i][j][0]+bbdb,params[i][j][1],params[i][j][2]] )
				params[i][j] = [temp[0],temp[1],temp[2],params[i][j][3],params[i][j][4]]


def proj_ali_incore_multi(data, refrings, numr, xrng = 0.0, yrng = 0.0, step=1.0, an = 1.0, nsoft = -1, finfo=None, sym="c1"):
	from sp_utilities    import compose_transform2
	from math         import cos, pi, radians, degrees
	from EMAN2 import Vec2f, Transform
	from sp_global_def import Util
	from sp_global_def import ERROR

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
		
	ou = numr[-3]
	sxi = dp["tx"]
	syi = dp["ty"]
	txrng = [0.0]*2 
	tyrng = [0.0]*2
	ERROR("proj_ali_incore_multi","Needs corrections",1)
	txrng[0] = max(0,min(cnx+sxi-ou, xrng+sxi))
	txrng[1] = max(0, min(nx-cnx-sxi-ou, xrng-sxi))
	tyrng[0] = max(0,min(cny+syi-ou, yrng+syi))
	tyrng[1] = max(0, min(ny-cny-syi-ou, yrng-syi))
		
	#print "Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"%(dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"])
	#[ang, sxs, sys, mirror, iref, peak, checked_refs] = Util.shc(data, refrings, xrng, yrng, step, ant, mode, numr, cnx+dp["tx"], cny+dp["ty"])
	peaks = Util.multiref_polar_ali_2d_peaklist_local(data, refrings, txrng, tyrng, step, ant, mode, numr, cnx+sxi, cny+syi)
	peaks_count = len(peaks) / 5
	#pixel_error = 0.0
	peak = 0.0
	if( peaks_count > 0 ):
		if( nsoft == -1 ):  nsoft = peaks_count
		params = [None]*peaks_count
		#                                              peak         iref      ang  sxs  sys 
		for i in range(peaks_count):  params[i] = [ peaks[i*5+0], int(peaks[i*5+4]), peaks[i*5+1], peaks[i*5+2], peaks[i*5+3]]
		params.sort(reverse=True)
		if(peaks_count < nsoft ):
			for i in range(peaks_count,nsoft,1): params.insert(0,params[0])
			peaks_count = nsoft
		elif( peaks_count > nsoft ):  peaks_count = min(peaks_count, nsoft)
		ws = sum([params[i][0] for i in range(peaks_count)])
		for i in range(peaks_count):
			iref   = params[i][1]
			ang    = params[i][2]
			sxs    = params[i][3]
			sys    = params[i][4]
			#mirror = 0
			peak   = params[i][0]/ws
			# The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
			# What that means is that one has to change the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
			angb, sxb, syb, ct = compose_transform2(0.0, sxs, sys, 1, -ang, 0.0, 0.0, 1)
			"""Multiline Comment17"""
			#MULTILINEMULTILINEMULTILINE 17
				#MULTILINEMULTILINEMULTILINE 17
				#MULTILINEMULTILINEMULTILINE 17
				#MULTILINEMULTILINEMULTILINE 17
				#MULTILINEMULTILINEMULTILINE 17
				#MULTILINEMULTILINEMULTILINE 17
			#MULTILINEMULTILINEMULTILINE 17
			#MULTILINEMULTILINEMULTILINE 17
			phi   = refrings[iref].get_attr("phi")
			theta = refrings[iref].get_attr("theta")
			psi   = (refrings[iref].get_attr("psi")+angb+360.0)%360.0
			s2x   = sxb - sxi
			s2y   = syb - syi

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
	"""Multiline Comment18"""
		#MULTILINEMULTILINEMULTILINE 18

	#MULTILINEMULTILINEMULTILINE 18
	#MULTILINEMULTILINEMULTILINE 18

def shc_multi(data, refrings, numr, xrng, yrng, step, an, nsoft, sym, finfo=None):
	from sp_utilities    import compose_transform2
	from sp_fundamentals import mirror
	from math         import cos, pi
	from EMAN2 import Vec2f, Transform
	from sp_global_def import Util

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

	ou = numr[-3]
	txrng = [0.0]*2 
	tyrng = [0.0]*2
	txrng[0] = max(0,min(cnx-ou, xrng))
	txrng[1] = max(0, min(nx-cnx-ou, xrng))
	tyrng[0] = max(0,min(cny-ou, yrng))
	tyrng[1] = max(0, min(ny-cny-ou, yrng))
		
	peaks = Util.shc_multipeaks(data, refrings, txrng, tyrng, step, ant, mode, numr, cnx, cny, nsoft)
	peaks_count = len(peaks) / 7
	pixel_error = 0.0
	number_of_checked_refs = 0
	peak = 0.0
	if( peaks_count > 0 ):
		params = [None]*peaks_count
		#                                              peak         iref                  ang        sxs           sys           mir           checked references
		for i in range(peaks_count):  params[i] = [ peaks[i*7+5], int(peaks[i*7+4]), peaks[i*7+0], peaks[i*7+1], peaks[i*7+2], int(peaks[i*7+3]), int(peaks[i*7+6])]
		#  Make sure nothing is repeated
		if(peaks_count>1):
			taken = [params[k][1] for k in range(peaks_count)]
			from sp_utilities import findall
			i = 0
			while(i<peaks_count):
				ll = findall(taken[i], taken)
				if(len(ll) > 1):
					sp_global_def.sxprint("  PROBLEM, found the same orientation more than once !  ")
					for k in range(len(params)):  sp_global_def.sxprint(params[k])
					ll.sort(reverse=True)
					for k in range(0,len(ll)-1):
						del params[k]
						peaks_count -= 1
					taken = [params[k][1] for k in range(peaks_count)]
				i+=1
		params.sort(reverse=True)
		ws = sum([params[i][0] for i in range(peaks_count)])  # peaks could be stretched
		for i in range(peaks_count):
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
			from sp_pixel_error import max_3D_pixel_error
			pixel_error += max_3D_pixel_error(t1, t2, numr[-3])
			#  preserve params, they might be needed if peaks_count<nsoft
			params[i] = [params[i][0], phi, theta, psi, s2x, s2y, iref]

		# Now set previousmax to a value halfway through
		data.set_attr("previousmax", params[peaks_count//2][0])

		#  Add orientations around the main peak with exclusion of those already taken
		#  Allow matching to bsoft>nsoft, but still keep nsoft.  This should allow elongated neighborhood
		bsoft = 2*nsoft
		if(peaks_count<nsoft):
			tempref = [refrings[i] for i in range(len(refrings))]
			taken   = [params[i][6] for i in range(peaks_count)]
			taken.sort(reverse=True)
			if(len(taken) > 1):
				for k in range(1,len(taken)):
					dod = []
					if( taken[k] == taken[k-1] ):
						sp_global_def.sxprint("  PROBLEM 2, entries duplicated  ",taken)
						dod.append(k)
				if(len(dod) >0):
					for k in dod:  del taken[k]
			#  delete taken
			try:
				for i in range(peaks_count):  del  tempref[taken[i]]
			except:
				sp_global_def.sxprint("  failed deleting tempref ")
				sp_global_def.sxprint(i,peaks_count,nsoft)
				sp_global_def.sxprint(" taken ",taken)
				sp_global_def.sxprint(len(tempref), len(refrings))
				from sys import exit
				exit()

			from sp_utilities import getfvec
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
				for i in range(iq):
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
				for i in range(iq,2*iq):
					n1 = tempref[i].get_attr("n1")
					n2 = tempref[i].get_attr("n2")
					n3 = tempref[i].get_attr("n3")
					refvecs[3*i+0] = n1
					refvecs[3*i+1] = n2
					refvecs[3*i+2] = n3
			else: 
				refvecs = [None]*3*len(tempref)
				for i in range(len(tempref)):
					n1 = tempref[i].get_attr("n1")
					n2 = tempref[i].get_attr("n2")
					n3 = tempref[i].get_attr("n3")
					refvecs[3*i+0] = n1
					refvecs[3*i+1] = n2
					refvecs[3*i+2] = n3
			from sp_utilities import nearestk_to_refdir
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
			pws = proj_ali_incore_multi(tdata, [tempref[k] for k in nrst], numr, xrng, yrng, step, 180.0, bsoft-peaks_count, sym=sym)
			#  Can there be a problem with (0,0) direction??  PAP  05/25/2014
			for i in range(bsoft-peaks_count):
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
			ws = sum([params[i][0] for i in range(nsoft)])  # peaks could be stretched
			for i in range(nsoft):
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
		"""Multiline Comment19"""
		#MULTILINEMULTILINEMULTILINE 19
		#MULTILINEMULTILINEMULTILINE 19
		#MULTILINEMULTILINEMULTILINE 19
			#MULTILINEMULTILINEMULTILINE 19
			#MULTILINEMULTILINEMULTILINE 19
		#MULTILINEMULTILINEMULTILINE 19
		#MULTILINEMULTILINEMULTILINE 19
			#MULTILINEMULTILINEMULTILINE 19
			#MULTILINEMULTILINEMULTILINE 19
		#MULTILINEMULTILINEMULTILINE 19
		pixel_error /= peaks_count
		peak = params[0][0]  # It is not used anywhere, but set it to the maximum.
	
	#  if it did not find any higher peaks would do nothing and return peaks_count=0
	return peak, pixel_error, number_of_checked_refs, peaks_count


# parameters: list of (all) projections | reference volume | ...
def ali3d_multishc_soft(stack, ref_vol, ali3d_options, mpi_comm = None, log = None, nsoft=2 ):

	from sp_alignment       import Numrinit, prepare_refrings, proj_ali_incore_local
	from sp_utilities       import get_im, file_type, model_circle, get_input_from_string, get_params_proj, wrap_mpi_gatherv, wrap_mpi_bcast
	from mpi             import mpi_bcast, mpi_comm_size, mpi_comm_rank, MPI_FLOAT, MPI_COMM_WORLD, mpi_barrier, mpi_reduce, MPI_INT, MPI_SUM
	from sp_projection      import prep_vol
	from sp_statistics      import hist_list
	from sp_applications    import MPI_start_end
	from sp_filter          import filt_ctf
	from sp_global_def      import Util
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
		from sp_logger import Logger
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

	if( type(ref_vol) is bytes ):  vol = get_im(ref_vol)
	else:	vol = ref_vol
	nx      = vol.get_xsize()
	if last_ring < 0:	last_ring = int(nx/2) - 2

	numr	= Numrinit(first_ring, last_ring, rstep, "F")
	mask2D  = model_circle(last_ring,nx,nx) - model_circle(first_ring,nx,nx)

	if( type(stack) is bytes ):
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
			nima = EMUtil.get_image_count(stack)
			list_of_particles = list(range(nima))
			
			total_nima = len(list_of_particles)
		else:
			list_of_particles = None
			total_nima = 0

	else:
		if myid == main_node:
			list_of_particles = list(range(len(stack)))
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

	if( type(stack) is bytes ):  data = EMData.read_images(stack, list_of_particles)
	else:                                   data = [ stack[im] for im in list_of_particles ]
	for im in range(nima):
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
	for N_step in range(lstp):

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
				for im in range(nima):
					previousmax = data[im].get_attr_default("previousmax", -1.0e23)
					if(previousmax == -1.0e23):
						peak, pixer[im] = proj_ali_incore_local(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],10.0, sym=sym)
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
			par_r = [0]*max(2,(nsoft+1))
			for im in range(nima):
				ERROR("shc_multi","Needs corrections")
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
			par_r = mpi_reduce(par_r, len(par_r), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
			#total_checked_refs = wrap_mpi_gatherv([number_of_checked_refs], main_node, mpi_comm)
			terminate = 0
			if myid == main_node:
				#total_checked_refs = sum(total_checked_refs)
				log.add("=========== Number of better peaks found ==============")
				for lhx in range(nsoft+1):
					msg = "            %5d     %7d"%(lhx, par_r[lhx])
					log.add(msg)
				log.add("_______________________________________________________")

				lhist = 20
				region, histo = hist_list(all_pixer, lhist)
				log.add("=========== Histogram of pixel errors ==============")
				for lhx in range(lhist):
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
							sp_global_def.sxprint(" NO XFORM  ",myid, i,im.get_attr('ID'))
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
		log.add("Finished ali3d_multishc_soft")
		return #params, vol, previousmax, par_r
	else:
		return #None, None, None, None  # results for the other processes


# data - projections (scattered between cpus) or the volume.  If volume, just do the volume processing
# options - the same for all cpus
# return - volume the same for all cpus
def do_volume(data, options, iter, mpi_comm):
	from EMAN2          import Util
	from mpi            import mpi_comm_rank
	from sp_filter       import filt_table
	from sp_reconstruction import recons3d_4nn_MPI, recons3d_4nn_ctf_MPI
	from sp_utilities      import bcast_EMData_to_all
	import types
	
	myid = mpi_comm_rank(mpi_comm)
	sym  = options.sym
	sym = sym[0].lower() + sym[1:]
	npad      = options.npad
	CTF       = options.CTF
	snr       = options.snr
	#=========================================================================
	# volume reconstruction
	if( type(data) == list ):
		if CTF: vol = recons3d_4nn_ctf_MPI(myid, data, snr, symmetry=sym, npad=npad, mpi_comm=mpi_comm)
		else:   vol = recons3d_4nn_MPI    (myid, data,      symmetry=sym, snr=snr, npad=npad, mpi_comm=mpi_comm)
	else:
		vol = data

	if myid == 0:
		from sp_morphology import threshold
		from sp_filter     import filt_tanl, filt_btwl
		from sp_utilities  import model_circle, model_cylinder, get_im
		import types
		nx = vol.get_xsize()
		if options.filament_width != -1:
			mask3D = model_cylinder(int(options.filament_width*options.resample_ratio+0.5)//2, nx, nx, nx)
			options.mask3D = mask3D
		elif(options.mask3D == None):
			last_ring   = int(options.ou)
			mask3D = model_circle(last_ring, nx, nx, nx)
		elif(options.mask3D == "auto"):
			from sp_utilities import adaptive_mask
			mask3D = adaptive_mask(vol)
		else:
			if( type(options.mask3D) == bytes ):  mask3D = get_im(options.mask3D)
			else:  mask3D = (options.mask3D).copy()
			nxm = mask3D.get_xsize()
			if( nx != nxm):
				from sp_fundamentals import rot_shift3D
				mask3D = Util.window(rot_shift3D(mask3D,scale=float(nx)/float(nxm)),nx,nx,nx)
				nxm = mask3D.get_xsize()
				assert(nx == nxm)

		stat = Util.infomask(vol, mask3D, False)
		vol -= stat[0]
		Util.mul_scalar(vol, 1.0/stat[1])
		vol = threshold(vol)
		#Util.mul_img(vol, mask3D)
		if( options.pwreference ):
			from sp_utilities    import read_text_file
			from sp_fundamentals import rops_table, fftip, fft
			rt = read_text_file( options.pwreference )
			fftip(vol)
			ro = rops_table(vol)
			#  Here unless I am mistaken it is enough to take the beginning of the reference pw.
			for i in range(1,len(ro)):  ro[i] = (rt[i]/ro[i])**0.5
			if( type(options.fl) == list ):
				vol = fft( filt_table( filt_table(vol, options.fl), ro) )
			else:
				vol = fft( filt_table( filt_tanl(vol, options.fl, options.aa), ro) )
		else:
			if( type(options.fl) == list ):
				vol = filt_table(vol, options.fl)
			else:
				vol = filt_tanl(vol, options.fl, options.aa)
		stat = Util.infomask(vol, mask3D, False)
		vol -= stat[0]
		Util.mul_scalar(vol, 1.0/stat[1])
		vol = threshold(vol)
		vol = filt_btwl(vol, 0.38, 0.5)
		Util.mul_img(vol, mask3D)
		del mask3D
		# vol.write_image('toto%03d.hdf'%iter)
	# broadcast volume
	bcast_EMData_to_all(vol, myid, 0, comm=mpi_comm)
	#=========================================================================
	return vol



def no_of_processors_restricted_by_data__do_volume(projections, ali3d_options, iter, mpi_comm):
	from mpi import mpi_comm_rank, mpi_comm_size, mpi_finalize, mpi_comm_split, mpi_barrier, MPI_COMM_WORLD
	from sp_utilities      import bcast_EMData_to_all
	from sp_applications import MPI_start_end

	mpi_size = mpi_comm_size(mpi_comm)
	n_projs = len(projections)
	mpi_rank = mpi_comm_rank(mpi_comm)

	if (mpi_size > n_projs):
		working = int(not(mpi_rank < n_projs))
		mpi_subcomm = mpi_comm_split(mpi_comm, working,  mpi_rank - working*n_projs)
		mpi_subsize = mpi_comm_size(mpi_subcomm)
		mpi_subrank = mpi_comm_rank(mpi_subcomm)
		if (mpi_rank < n_projs):
			proj_begin, proj_end = MPI_start_end(n_projs, mpi_subsize, mpi_subrank)
			ref_vol = do_volume(projections[proj_begin:proj_end], ali3d_options, 0, mpi_comm=mpi_subcomm)
		else:
			from sp_utilities import model_blank
			nx = projections[0].get_xsize()
			ref_vol = model_blank(nx,nx,nx)
		bcast_EMData_to_all(ref_vol, mpi_rank, 0, comm=mpi_comm)
	else:
		proj_begin, proj_end = MPI_start_end(n_projs, mpi_size, mpi_rank)
		ref_vol = do_volume(projections[proj_begin:proj_end], ali3d_options, 0, mpi_comm=mpi_comm)

	return ref_vol


"""Multiline Comment20"""
#MULTILINEMULTILINEMULTILINE 20
#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20


#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20




#MULTILINEMULTILINEMULTILINE 20
#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20
		 #MULTILINEMULTILINEMULTILINE 20
		 #MULTILINEMULTILINEMULTILINE 20
		 #MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20


	#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20

#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
		 #MULTILINEMULTILINEMULTILINE 20
		 #MULTILINEMULTILINEMULTILINE 20
		 #MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
		 #MULTILINEMULTILINEMULTILINE 20
		 #MULTILINEMULTILINEMULTILINE 20
		 #MULTILINEMULTILINEMULTILINE 20
		 #MULTILINEMULTILINEMULTILINE 20
		 #MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20



#MULTILINEMULTILINEMULTILINE 20
#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20

		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20

			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20

			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20

			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20

			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20

			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20

			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20
							#MULTILINEMULTILINEMULTILINE 20
							#MULTILINEMULTILINEMULTILINE 20
							#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20
							#MULTILINEMULTILINEMULTILINE 20
							#MULTILINEMULTILINEMULTILINE 20
							#MULTILINEMULTILINEMULTILINE 20
							#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20

			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20

			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20

			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20

				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20

				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20
							#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20
							#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20

				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20

				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20
							#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20
							#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
						#MULTILINEMULTILINEMULTILINE 20

				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20

				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20
					#MULTILINEMULTILINEMULTILINE 20

			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
				#MULTILINEMULTILINEMULTILINE 20
			#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20

	#MULTILINEMULTILINEMULTILINE 20


	#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20
	#MULTILINEMULTILINEMULTILINE 20
		#MULTILINEMULTILINEMULTILINE 20



#MULTILINEMULTILINEMULTILINE 20




"""Multiline Comment21"""
#MULTILINEMULTILINEMULTILINE 21
#MULTILINEMULTILINEMULTILINE 21
#MULTILINEMULTILINEMULTILINE 21
#MULTILINEMULTILINEMULTILINE 21
#MULTILINEMULTILINEMULTILINE 21
#MULTILINEMULTILINEMULTILINE 21
#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21

	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21

	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21

	#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21

	#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21

	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21

	#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21

	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21

	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21

	#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21

	#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21


	#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21

	#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21

	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21

	#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21

	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21

	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21


	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21

		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21

			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21

			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21

			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21

			#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21

			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
							#MULTILINEMULTILINEMULTILINE 21
							#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21

			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
														#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
														#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
												#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21

			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21

			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21

				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
							#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
							#MULTILINEMULTILINEMULTILINE 21
							#MULTILINEMULTILINEMULTILINE 21
							#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21

			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21

			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21

			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21


				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21

						#MULTILINEMULTILINEMULTILINE 21
							#MULTILINEMULTILINEMULTILINE 21
							#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
							#MULTILINEMULTILINEMULTILINE 21
							#MULTILINEMULTILINEMULTILINE 21
							#MULTILINEMULTILINEMULTILINE 21

						#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21


			#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
						#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
				#MULTILINEMULTILINEMULTILINE 21
					#MULTILINEMULTILINEMULTILINE 21
			#MULTILINEMULTILINEMULTILINE 21


	#MULTILINEMULTILINEMULTILINE 21
		#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21
	#MULTILINEMULTILINEMULTILINE 21

#MULTILINEMULTILINEMULTILINE 21
"""Multiline Comment22"""
#MULTILINEMULTILINEMULTILINE 22
#MULTILINEMULTILINEMULTILINE 22
										   #MULTILINEMULTILINEMULTILINE 22
	#MULTILINEMULTILINEMULTILINE 22
	#MULTILINEMULTILINEMULTILINE 22
	#MULTILINEMULTILINEMULTILINE 22
	#MULTILINEMULTILINEMULTILINE 22
	#MULTILINEMULTILINEMULTILINE 22
	#MULTILINEMULTILINEMULTILINE 22
	#MULTILINEMULTILINEMULTILINE 22
	#MULTILINEMULTILINEMULTILINE 22
	#MULTILINEMULTILINEMULTILINE 22
	#MULTILINEMULTILINEMULTILINE 22
	#MULTILINEMULTILINEMULTILINE 22
		#MULTILINEMULTILINEMULTILINE 22
		#MULTILINEMULTILINEMULTILINE 22
			#MULTILINEMULTILINEMULTILINE 22
				#MULTILINEMULTILINEMULTILINE 22
			#MULTILINEMULTILINEMULTILINE 22
				#MULTILINEMULTILINEMULTILINE 22
			#MULTILINEMULTILINEMULTILINE 22
			#MULTILINEMULTILINEMULTILINE 22
				#MULTILINEMULTILINEMULTILINE 22
			#MULTILINEMULTILINEMULTILINE 22
				#MULTILINEMULTILINEMULTILINE 22
	#MULTILINEMULTILINEMULTILINE 22
		#MULTILINEMULTILINEMULTILINE 22
			#MULTILINEMULTILINEMULTILINE 22
				#MULTILINEMULTILINEMULTILINE 22
			#MULTILINEMULTILINEMULTILINE 22
				#MULTILINEMULTILINEMULTILINE 22
	#MULTILINEMULTILINEMULTILINE 22
		#MULTILINEMULTILINEMULTILINE 22
	#MULTILINEMULTILINEMULTILINE 22
#MULTILINEMULTILINEMULTILINE 22
"""Multiline Comment23"""
#MULTILINEMULTILINEMULTILINE 23
	#MULTILINEMULTILINEMULTILINE 23
	#MULTILINEMULTILINEMULTILINE 23
	#MULTILINEMULTILINEMULTILINE 23
	#MULTILINEMULTILINEMULTILINE 23
#MULTILINEMULTILINEMULTILINE 23
from builtins import range
import sp_global_def
