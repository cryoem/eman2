
class mpi_env_type:
	main_comm = None
	main_rank = None
	sub_comm = None
	sub_rank = None
	subcomm_id = None
	subcomms_count = None
	subcomms_roots = None
	


def kernel(projections, stable_subset, target_threshold, options, minimal_subset_size, number_of_runs, number_of_winners, mpi_env, log, prefix=""):
	from multi_shc import multi_shc, find_common_subset_3
	from utilities import wrap_mpi_bcast, write_text_row, write_text_file, wrap_mpi_gatherv, average_angles
	from itertools import combinations
	import os
	
	if log == None:
		from logger import Logger
		log = Logger()
	
	stable_subset = wrap_mpi_bcast(stable_subset, 0, mpi_env.main_comm)
	
	if mpi_env.main_rank == 0:
		log.add("Start ", number_of_runs, "* 3SHC")
		for i in xrange(number_of_runs):
			log.add("3SHC --> " + log.prefix + prefix + "_" + str(i))

	completed_mshc = 0
	params = []
	while completed_mshc < number_of_runs:
		runs_to_do = min([ (number_of_runs - completed_mshc), mpi_env.subcomms_count ])
		if mpi_env.subcomm_id < runs_to_do:
			out_dir = prefix + "_" + str(completed_mshc + mpi_env.subcomm_id)
			if mpi_env.sub_rank == 0:
				os.mkdir(log.prefix + out_dir)
			out_params, out_vol, out_peaks = multi_shc(projections, stable_subset, 3, options, mpi_env.sub_comm, log=log.sublog(out_dir + "/"))
		else:
			out_params = None
		if mpi_env.main_rank in mpi_env.subcomms_roots and mpi_env.subcomm_id < runs_to_do:
			params_temp = wrap_mpi_gatherv([out_params], 0, mpi_env.main_comm)
		else:
			params_temp = wrap_mpi_gatherv([], 0, mpi_env.main_comm)
		if mpi_env.main_rank == 0:
			params.extend(params_temp)
		completed_mshc += runs_to_do

	# find common subset
	if mpi_env.main_rank == 0:
		log.add("Calculate common subset")
		best_confs = []
		largest_subset = []
		largest_subset_error = 999.0
		msg = ""
		for it in combinations(range(number_of_runs), number_of_winners):
			confs = list(it)
			input_params = []
			for c in confs:
				input_params.append(params[c])
			subset_thr, subset_size, err_thr, err_size = find_common_subset_3(input_params, target_threshold, minimal_subset_size, thresholds=True)
			msg += str(len(subset_size)) + "(" + str(err_size) + ") "
			if len(subset_size) > len(largest_subset) or ( len(subset_size) == len(largest_subset) and err_size < largest_subset_error):
				best_confs = confs
				largest_subset = subset_size
				largest_subset_error = err_size
		log.add(msg)
		subset = largest_subset
		threshold = largest_subset_error
		new_stable_subset = []
		for i in subset:
			new_stable_subset.append(stable_subset[i])
		log.add("Best solutions (winners): ", best_confs)
		average_params = []
		for i in subset:
			temp = [ params[j][i] for j in best_confs ]
			average_params.append(average_angles(temp))
		write_text_file(new_stable_subset, log.prefix + prefix + "_indexes.txt")
		write_text_row(average_params, log.prefix + prefix + "_params.txt")
	else:
		threshold = None
		new_stable_subset = None
	
	# broadcast threshold and new stable subset and exit
	threshold = wrap_mpi_bcast(threshold, 0, mpi_env.main_comm)
	new_stable_subset = wrap_mpi_bcast(new_stable_subset, 0, mpi_env.main_comm)
	return threshold, new_stable_subset


# INPUT:
# projections -> list of all projections (images) - the same for all processes in main communicator
# stable_subset -> list/set of indices of projections belonging to stable subset - the same for all processes in main communicator
# stable_threshold -> threshold error to use during recalculation of stable subset - the same for all processes in main communicator
# OUTPUT:
# new stable subset -> set of indices - the same for all processes in main communicator
def shrink_step(projections, stable_subset, target_threshold, options, min_stable_subset_size, number_of_runs, number_of_winners, mpi_env=None, log=None, iteration=-1):
	from utilities import wrap_mpi_bcast
	
	if log == None:
		from logger import Logger
		log = Logger()
	
	glob_iter = iteration
	
	if mpi_env.main_rank == 0:
		log.add("-------------> Contracting step - BEGIN")
	
	iteration = 0
	while True:
		iteration += 1
		if mpi_env.main_rank == 0:
			log.add("-----> Iteration", iteration)
		stable_threshold, out_subset = kernel(projections, stable_subset, target_threshold, options, min_stable_subset_size, number_of_runs, number_of_winners, mpi_env, log=log, prefix=(str(glob_iter) + "_contracting_" + str(iteration)))
		terminate = False
		if mpi_env.main_rank == 0:
			log.add("Subset", len(out_subset), out_subset)
			log.add("Threshold", stable_threshold)
			terminate = (out_subset == stable_subset)
			stable_subset = out_subset
		terminate = wrap_mpi_bcast(terminate, 0, mpi_env.main_comm)
		if terminate:
			break

	if mpi_env.main_rank == 0:
		log.add("-------------> Contracting step - END")

	return stable_subset, stable_threshold


#calculate t1*t2
def mult_transform(v1, v2):
	from EMAN2 import Transform
	T1 = Transform({"type":"spider","phi":v1[0],"theta":v1[1],"psi":v1[2],"tx":0.0,"ty":0.0,"tz":0.0,"mirror":0,"scale":1.0})
	T2 = Transform({"type":"spider","phi":v2[0],"theta":v2[1],"psi":v2[2],"tx":0.0,"ty":0.0,"tz":0.0,"mirror":0,"scale":1.0})
	T = T1*T2
	return [ T.get_params("spider")["phi"], T.get_params("spider")["theta"], T.get_params("spider")["psi"] ]


def calculate_diff(ang1, ang2):
	from utilities import angle_between_projections_directions, rotation_between_anglesets
	
	phi, theta, psi = rotation_between_anglesets(ang1, ang2)
	rot = [phi, theta, psi]
	n = len(ang1)
	diff = []
	for k in xrange(n):
		diff.append( angle_between_projections_directions( mult_transform(ang1[k],rot), ang2[k] ) )
	return diff


# edges - list of pairs (int, int)
# returns list of ints - max clique
def find_max_clique(edges):
	from global_def import Util
	
	edges_cpp = []
	for edg in edges:
		edges_cpp.append(edg[0])
		edges_cpp.append(edg[1])
	
	c = Util.max_clique(edges_cpp)
	c.sort()
	return c


# returns list of int - ids of correct configurations
def find_correct_configurations(params, subsets, threshold, log):
	
	log.add("Find aligned configurations")
	
	edges = []
	
	for iS1 in xrange(len(subsets)):
		for iS2 in xrange(iS1):
			s1 = subsets[iS1]
			s2 = subsets[iS2]
			indicies = []
			param1 = []
			param2 = []
			for iE1 in xrange(len(s1)):
				ind = s1[iE1]
				if ind not in set(s2):
					continue
				iE2 = s2.index(ind)
				indicies.append(ind)
#				print iS1, iE1, "  <--->  ", iS2, iE2
				param1.append( params[iS1][iE1] )
				param2.append( params[iS2][iE2] )
			if len(param1) > 2:
				diff = calculate_diff(param1, param2)
				if sum(diff) / len(diff) <= threshold:
					edges.append([iS1,iS2])
			else:
				edges.append([iS1,iS2])
	log.add("Number of edges (pairs of aligned configurations): ", len(edges))
	correct_configurations = find_max_clique(edges)
	log.add("Size of max-clique (number of correct configurations): ", len(correct_configurations))
	return correct_configurations


# returns list of int - ids of projections
def find_new_stable_projections(params, subsets, threshold, stable_projections, log):
	
	m = len(params)
	params_stable = [ [] for i in xrange(m) ]
	subset_stable = [ [] for i in xrange(m) ]
	for i in xrange(m):
		for j in xrange(len(params[i])):
			if subsets[i][j] in stable_projections:
				params_stable[i].append( params[i][j] )
				subset_stable[i].append( subsets[i][j] )
	correct_confs = find_correct_configurations(params_stable, subset_stable, threshold, log)

	new_params = []
	new_subsets = []
	for i in correct_confs:
		new_params.append(params[i])
		new_subsets.append(subsets[i])
	params = new_params
	subsets = new_subsets
	m = len(params)

	unstable_projections = set()
	candidates = set()

	for iS1 in xrange(m):
		for iS2 in xrange(iS1):
			s1 = subsets[iS1]
			s2 = subsets[iS2]
			indicies = []
			param1 = []
			param2 = []
			for iE1 in xrange(len(s1)):
				ind = s1[iE1]
				if ind not in set(s2):
					continue
				iE2 = s2.index(ind)
				indicies.append(ind)
#				print iS1, iE1, "  <--->  ", iS2, iE2
				param1.append( params[iS1][iE1] )
				param2.append( params[iS2][iE2] )
			stable_subset_size = len(param1)
			if stable_subset_size > 5:
				diff = calculate_diff(param1, param2)
				while max(diff) > threshold:
					to_delete = diff.index( max(diff) )
					if indicies[to_delete] not in stable_projections:
						unstable_projections.add(indicies[to_delete])
					del indicies[to_delete]
					del param1[to_delete]
					del param2[to_delete]
					diff = calculate_diff(param1, param2)
				for i in indicies:
					if i not in stable_projections:
						candidates.add(i)
	
	log.add("Number of candidates: ", len(candidates))
	log.add("Number of unstable_projections: ", len(unstable_projections))
	new_stable_projections = list(candidates - unstable_projections)
	log.add("New stable projections (set of candidates - set of unstable): ", len(new_stable_projections), new_stable_projections)
	return new_stable_projections


def recalculate_subset(threshold, params, subsets):
	
	if len(params) != len(subsets):
		print params
		print subsets
		assert(0 == 1)
	for i in xrange(len(params)):
		if len(params[i]) != len(subsets[i]):
			print params[i]
			print subsets[i]
			assert(0 == 2)
	
	# number of all projections
	N = 0
	for s in subsets:
		t = max(s) + 1
		if t > N:
			N = t
	
	errors = [ [] for i in xrange(N) ]
	
	subsets_score = [0]*len(subsets)
	matrix_of_diffs = [ [] for i in xrange(len(subsets)) ]
	matrix_of_indicies = [ [] for i in xrange(len(subsets)) ]
	for iS1 in xrange(len(subsets)):
		for iS2 in xrange(iS1):
			s1 = subsets[iS1]
			s2 = subsets[iS2]
			indicies = []
			param1 = []
			param2 = []
			for iE1 in xrange(len(s1)):
				ind = s1[iE1]
				if ind not in set(s2):
					continue
				iE2 = s2.index(ind)
				indicies.append(ind)
#				print iS1, iE1, "  <--->  ", iS2, iE2
				param1.append( params[iS1][iE1] )
				param2.append( params[iS2][iE2] )
			diff = calculate_diff(param1, param2)
			diff_per_subset = sum(diff) / len(diff)
			subsets_score[iS1] += diff_per_subset
			subsets_score[iS2] += diff_per_subset
			matrix_of_diffs[iS1].append( diff )
			matrix_of_indicies[iS1].append( indicies )
# 			print indicies, len(indicies)
# 			print diff, len(diff)
			#for i in xrange(len(diff)):
			#	errors[indicies[i]].append(diff[i])
	
	temp = subsets_score[:]
	temp.sort(reverse=True)
	temp_thr = temp[len(temp) / 3]
	
	for iS1 in xrange(len(subsets)):
		if subsets_score[iS1] > temp_thr:
			continue
		for iS2 in xrange(iS1):
			if subsets_score[iS2] > temp_thr:
				continue
			diff = matrix_of_diffs[iS1][iS2]
			indicies = matrix_of_indicies[iS1][iS2]
			for i in xrange(len(diff)):
				errors[indicies[i]].append(diff[i])

	new_subset = []
	for i in xrange(N):
		if len(errors[i]) > 1 and max(errors[i]) < threshold:
			new_subset.append(i)

	return new_subset


# From given list_of_indicies the ceil(len(list_of_indicies) * number_of_repetitions / trg_subset_size) subsets are generated.
# They meet the following conditions:
# - they contain X number of elements, where X is as close to the trg_subset_size as possible
# - each element from list_of_indicies is exactly in number_of_repetitions subsets
def generate_subsets(list_of_indicies, trg_subset_size, number_of_repetitions):
	from random import randrange

	N = len(list_of_indicies)
	list_of_indicies = list(list_of_indicies)
	ind = range(N)
	# TODO shuffle ind if necessary
	
	subsets_count = N * number_of_repetitions // trg_subset_size
	if N * number_of_repetitions % trg_subset_size > trg_subset_size / 2:
		subsets_count += 1
	
	subsets_avail = [ [] for x in xrange(subsets_count) ]
	subsets_frozen = []
	
	for i in ind:
		for r in xrange(number_of_repetitions):
			if len(subsets_avail) == 0:
				subsets_avail = subsets_frozen
				subsets_frozen = []
			while True:
				trg_ind = randrange(0, len(subsets_avail))
				if subsets_avail[trg_ind].count(list_of_indicies[i]) == 0:
					subset = subsets_avail.pop(trg_ind)
					break
			subset.append(list_of_indicies[i])
			subsets_frozen.append(subset)
	
	subsets_frozen.extend(subsets_avail)
	return subsets_frozen


# INPUT:
# projections -> list of all projections (images) - the same for all processes in main communicator
# stable_subset -> list/set of indices of projections belonging to stable subset - the same for all processes in main communicator
# stable_threshold -> threshold error to use during recalculation of stable subset - the same for all processes in main communicator
# OUTPUT:
# new stable subset -> set of indices - the same for all processes in main communicator
def expand_step(projections, stable_subset, stable_threshold, options, tries_per_unstable=4, mpi_env=None, log=None, iteration=-1):
	from applications import MPI_start_end
	from multi_shc import multi_shc
	from utilities import wrap_mpi_recv, wrap_mpi_send, wrap_mpi_bcast, wrap_mpi_gatherv
	import os
	
	if log == None:
		from logger import Logger
		logger = Logger()
	
	if mpi_env.main_rank == 0:
		log.add("-------------> Expanding step - BEGIN")
	
	# generate subsets (on main root)
	if mpi_env.main_rank == 0:
		nprojs = len(projections)
		stable_subset = set(stable_subset)
		unstable_subset = set(range(nprojs)) - stable_subset
		trg_subset_size = max( [len(unstable_subset) / 4, 1] )
		subsets = generate_subsets(unstable_subset, trg_subset_size, tries_per_unstable)
		for i in xrange(len(subsets)):
			subsets[i] = list( stable_subset | set(subsets[i]) )
		log.add("Count of subsets: ", len(subsets))

	# scatter subsets between sub-communicators (fragments of subsets are sent to roots of sub-communicators)
	assigned_subsets = None
	if mpi_env.main_rank == 0:
		for iSC in xrange( mpi_env.subcomms_count ):
			subset_begin, subset_end = MPI_start_end(len(subsets), mpi_env.subcomms_count, iSC)
			dest = mpi_env.subcomms_roots[iSC]
			if mpi_env.main_rank != dest:
				wrap_mpi_send( subsets[subset_begin:subset_end], dest, mpi_env.main_comm )
			else:
				assigned_subsets = subsets[subset_begin:subset_end]
	else:
		if mpi_env.sub_rank == 0:
			assigned_subsets = wrap_mpi_recv(0, mpi_env.main_comm)
	
	# broadcast subsets among sub-communicators (from roots to other processes)
	assigned_subsets = wrap_mpi_bcast(assigned_subsets, 0, mpi_env.sub_comm)

	# run SHC (in each sub-communicator separately)
	params = []
	for iAS in xrange(len(assigned_subsets)):
		out_dir = str(iteration) + "_expanding_" + str(mpi_env.subcomm_id) + "_" + str(iAS)
		if mpi_env.sub_rank == 0:
			os.mkdir(log.prefix + out_dir)
			log.add("3SHC --> " + log.prefix + out_dir)
		subset = assigned_subsets[iAS]
		out_p, out_vol, out_peaks = multi_shc(projections, subset, 3, options, mpi_comm=mpi_env.sub_comm, log=log.sublog(out_dir + "/"))
		if mpi_env.sub_rank == 0:
			assert( len(subset) == len(out_p) )
		if mpi_env.sub_rank == 0:
			params.append( out_p )

	# gather obtained parameters to main root
	if mpi_env.sub_rank == 0:
		params = wrap_mpi_gatherv(params, 0, mpi_env.main_comm)
	else:
		params = wrap_mpi_gatherv([], 0, mpi_env.main_comm)

	# calculate new stable subset
	if mpi_env.main_rank == 0:
		log.add("Calculate new common subset")
		#if prefix != None:
		#	for pi in xrange(len(params)):
		#		p = params[pi]
		#		s = subsets[pi]
		#		write_text_row(p, prefix + "_params_" + str(pi) + "_" + log.prefix + ".txt")
		#		write_text_file(s, prefix + "_subset_" + str(pi) + "_" + log.prefix + ".txt")
		#new_stable_subset = recalculate_subset(stable_threshold, params, subsets)
		new_stable_subset = find_new_stable_projections(params, subsets, stable_threshold, stable_subset, log=log)
		new_stable_subset = list( set(stable_subset) | set(new_stable_subset) )
		log.add("New subset:", len(new_stable_subset), new_stable_subset)
	else:
		new_stable_subset = None
	
	if mpi_env.main_rank == 0:
		log.add("-------------> Expanding step - END")
	
	return new_stable_subset


def air(projs, minimal_subset_size, target_threshold, options, number_of_runs, number_of_winners, mpi_env, log):
	from utilities import wrap_mpi_bcast
	
	subset = range(len(projs))
	new_subset = range(len(projs))
	threshold = target_threshold
	
	if mpi_env.main_rank == 0:
		log.add("Minimal subset size", minimal_subset_size)
		log.add("Input set", len(projs))
		log.add("Target threshold", target_threshold)
	
	iteration = 0;
	while True:
		iteration += 1
		if mpi_env.main_rank == 0:
			log.add("====================== ITERATION ", iteration, "=======================")
		new_subset, new_threshold = shrink_step(projs, new_subset, target_threshold, options, minimal_subset_size, number_of_runs, number_of_winners, mpi_env=mpi_env, log=log, iteration=iteration)
		terminate = False
		if mpi_env.main_rank == 0:
			log.add("Subset", len(new_subset), new_subset)
			log.add("Threshold", new_threshold)
			if set(subset) == set(new_subset) and threshold >= new_threshold:
				terminate = True
				log.add("New subset is the same as the previous one (and new_threshold <= prev_threshold) - END")
				log.add("=====================================================================")
		terminate = wrap_mpi_bcast(terminate, 0, mpi_env.main_comm)
		if terminate:
			break
		if mpi_env.main_rank == 0:
			subset = new_subset
			threshold = new_threshold
		new_subset = expand_step(projs, subset, threshold, options, mpi_env=mpi_env, log=log, iteration=iteration)
		if mpi_env.main_rank == 0:
			log.add("Subset", len(new_subset), new_subset)
			log.add("Threshold", new_threshold)

	return new_subset, new_threshold

