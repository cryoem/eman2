#!/usr/bin/env python

import global_def
from global_def import *
from mpi import MPI_SUM, mpi_reduce, mpi_init, mpi_finalize, MPI_COMM_WORLD, mpi_comm_rank, mpi_comm_size, mpi_barrier, \
	mpi_comm_split, mpi_bcast, MPI_INT, MPI_CHAR, MPI_FLOAT

from utilities import get_im, string_found_in_file, get_latest_directory_increment_value, store_value_of_simple_vars_in_json_file
from utilities import cmdexecute, if_error_all_processes_quit_program
from utilities import read_text_row, read_text_file, write_text_file, write_text_row, getindexdata, print_program_start_information
from multi_shc import find_common_subset, do_volume, multi_shc

import string
import os, sys
#from debug_mpi import mpi_barrier, mpi_bcast 


MAXIMUM_NO_OF_VIPER_RUNS_ANALYZED_TOGETHER = 10
# NORMALIZED_AREA_THRESHOLD_FOR_OUTLIER_DETECTION = 0.2
PERCENT_THRESHOLD_X = .8
PERCENT_THRESHOLD_Y = .2
ANGLE_ERROR_THRESHOLD = 24
TRIPLET_WITH_ANGLE_ERROR_LESS_THAN_THRESHOLD_HAS_BEEN_FOUND = -100
MUST_END_PROGRAM_THIS_ITERATION = -101
EMPTY_VIPER_RUN_INDICES_LIST = -102
DUMMY_INDEX_USED_AS_BUFFER = -103

def calculate_list_of_independent_viper_run_indices_used_for_outlier_elimination(no_of_viper_runs_analyzed_together, 
	no_of_viper_runs_analyzed_together_from_user_options, masterdir, rviper_iter, criterion_name):

	from utilities import combinations_of_n_taken_by_k

	# generate all possible combinations of (no_of_viper_runs_analyzed_together - 1) taken (3 - 1) at a time
	import itertools

	number_of_additional_combinations_for_this_viper_iteration = combinations_of_n_taken_by_k(no_of_viper_runs_analyzed_together - 1,
																		  no_of_viper_runs_analyzed_together_from_user_options - 1)

	criterion_measure = [0.0] * number_of_additional_combinations_for_this_viper_iteration
	all_n_minus_1_combinations_taken_k_minus_1_at_a_time = list(itertools.combinations(range(no_of_viper_runs_analyzed_together - 1),
																  no_of_viper_runs_analyzed_together_from_user_options - 1))

	no_of_processors = mpi_comm_size(MPI_COMM_WORLD)
	my_rank = mpi_comm_rank(MPI_COMM_WORLD)

	for idx, tuple_of_projection_indices in enumerate(all_n_minus_1_combinations_taken_k_minus_1_at_a_time):
		if (my_rank == idx % no_of_processors):
			list_of_viper_run_indices = list(tuple_of_projection_indices) + [no_of_viper_runs_analyzed_together - 1]
			criterion_measure[idx] = measure_for_outlier_criterion(criterion_name, masterdir, rviper_iter, list_of_viper_run_indices)
			plot_errors_between_any_number_of_projections(masterdir, rviper_iter, list_of_viper_run_indices, criterion_measure[idx])

	criterion_measure = mpi_reduce(criterion_measure, number_of_additional_combinations_for_this_viper_iteration, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD)

	if (my_rank == 0):
		index_of_sorted_criterion_measure_list = [i[0] for i in sorted(enumerate(criterion_measure), reverse=False, key=lambda x: x[1])]

		list_of_viper_run_indices_for_the_current_rrr_viper_iteration = list(all_n_minus_1_combinations_taken_k_minus_1_at_a_time[index_of_sorted_criterion_measure_list[0]]) + \
																		[no_of_viper_runs_analyzed_together - 1]

		mainoutputdir = masterdir + "/main%03d/" % (rviper_iter)

		if criterion_measure[index_of_sorted_criterion_measure_list[0]] == TRIPLET_WITH_ANGLE_ERROR_LESS_THAN_THRESHOLD_HAS_BEEN_FOUND:
			list_of_viper_run_indices_for_the_current_rrr_viper_iteration.insert(0,MUST_END_PROGRAM_THIS_ITERATION)
		else:
			list_of_viper_run_indices_for_the_current_rrr_viper_iteration.insert(0,DUMMY_INDEX_USED_AS_BUFFER)
			if criterion_name == "80th percentile":
				pass_criterion = criterion_measure[index_of_sorted_criterion_measure_list[0]] < PERCENT_THRESHOLD_Y
			elif criterion_name == "fastest increase in the last quartile":
				pass_criterion = criterion_measure[index_of_sorted_criterion_measure_list[-1]] > PERCENT_THRESHOLD_Y
			else:
				pass_criterion = False
	
			if not pass_criterion:
				list_of_viper_run_indices_for_the_current_rrr_viper_iteration = [EMPTY_VIPER_RUN_INDICES_LIST]

		import json; f = open(mainoutputdir + "list_of_viper_runs_included_in_outlier_elimination.json", 'w')
		json.dump(list_of_viper_run_indices_for_the_current_rrr_viper_iteration[1:],f); f.close()

		mpi_barrier(MPI_COMM_WORLD)
		return list_of_viper_run_indices_for_the_current_rrr_viper_iteration

	mpi_barrier(MPI_COMM_WORLD)

	return [EMPTY_VIPER_RUN_INDICES_LIST]


def identify_outliers(myid, main_node, rviper_iter, no_of_viper_runs_analyzed_together, 
	no_of_viper_runs_analyzed_together_from_user_options, masterdir, bdb_stack_location, outlier_percentile, 
	criterion_name, outlier_index_threshold_method):
	
	no_of_viper_runs_analyzed_together_must_be_incremented = 0
	do_calculation = 1

	if (myid == main_node):
		mainoutputdir = masterdir + "/main%03d/" % (rviper_iter)
		if(os.path.exists(mainoutputdir + "/list_of_viper_runs_included_in_outlier_elimination.json")):
			# list_of_independent_viper_run_indices_used_for_outlier_elimination = map(int, read_text_file(mainoutputdir + "/list_of_viper_runs_included_in_outlier_elimination.txt"))
			import json; f = open(mainoutputdir + "list_of_viper_runs_included_in_outlier_elimination.json", 'r')
			list_of_independent_viper_run_indices_used_for_outlier_elimination  = json.load(f); f.close()
			do_calculation = 0
		do_calculation = mpi_bcast(do_calculation, 1, MPI_INT, 0, MPI_COMM_WORLD)[0]
	else:
		do_calculation = mpi_bcast(do_calculation, 1, MPI_INT, 0, MPI_COMM_WORLD)[0]

	if do_calculation:
		list_of_independent_viper_run_indices_used_for_outlier_elimination = calculate_list_of_independent_viper_run_indices_used_for_outlier_elimination(no_of_viper_runs_analyzed_together,
			no_of_viper_runs_analyzed_together_from_user_options, masterdir, rviper_iter, criterion_name)

	# only master has the actual list: list_of_independent_viper_run_indices_used_for_outlier_elimination
	# only master has the actual list: list_of_independent_viper_run_indices_used_for_outlier_elimination
	# only master has the actual list: list_of_independent_viper_run_indices_used_for_outlier_elimination

	error_status = 0
	if (myid == main_node):
		# if len(list_of_independent_viper_run_indices_used_for_outlier_elimination) == 0:
		if list_of_independent_viper_run_indices_used_for_outlier_elimination[0] == EMPTY_VIPER_RUN_INDICES_LIST:
			if no_of_viper_runs_analyzed_together > MAXIMUM_NO_OF_VIPER_RUNS_ANALYZED_TOGETHER:
				error_status = 1
				cmd = "{} {}".format("mkdir ", masterdir + "MAXIMUM_NO_OF_VIPER_RUNS_ANALYZED_TOGETHER__Reached"); cmdexecute(cmd)
			else:
				# No set of solutions has been found to make a selection for outlier elimination.
				# A new independent viper run will be performed
				no_of_viper_runs_analyzed_together_must_be_incremented = 1
				cmd = "{} {}".format("rm ", mainoutputdir + "list_of_viper_runs_included_in_outlier_elimination.json"); cmdexecute(cmd)

		else:
			# Outliers are eliminated based on the viper runs contained in "list_of_independent_viper_run_indices_used_for_outlier_elimination"
			if list_of_independent_viper_run_indices_used_for_outlier_elimination[0] == MUST_END_PROGRAM_THIS_ITERATION:
				no_of_viper_runs_analyzed_together_must_be_incremented = MUST_END_PROGRAM_THIS_ITERATION
				found_outliers(list_of_independent_viper_run_indices_used_for_outlier_elimination[1:], outlier_percentile, 
					rviper_iter, masterdir, bdb_stack_location, "use all images")
			else:
				# still need to eliminate DUMMY_INDEX_USED_AS_BUFFER
				found_outliers(list_of_independent_viper_run_indices_used_for_outlier_elimination[1:], outlier_percentile, 
					rviper_iter, masterdir, bdb_stack_location, outlier_index_threshold_method)

	if_error_all_processes_quit_program(error_status)

	no_of_viper_runs_analyzed_together_must_be_incremented = mpi_bcast(no_of_viper_runs_analyzed_together_must_be_incremented, 1, MPI_INT, 0, MPI_COMM_WORLD)[0]

	return no_of_viper_runs_analyzed_together_must_be_incremented


def plot_errors_between_any_number_of_projections(masterdir, rviper_iter, list_of_projection_indices, error_value):
	import matplotlib.pyplot as plt

	# # for debugging purposes
	# if "counter" not in plot_errors_between_any_number_of_projections.__dict__:
	# 	plot_errors_between_any_number_of_projections.counter = 0
	# plot_errors_between_any_number_of_projections.counter += 1

	# main_iterations = ["main" + "%03d" % i for i in range(1, rviper_iter + 1)]
	main_iterations = ["main" + "%03d" % i for i in range(rviper_iter, rviper_iter + 1)]

	mainoutputdir = masterdir + "/" + main_iterations[0] + "/"
	p=[]
	for i1 in list_of_projection_indices:
		p.append(read_text_row(mainoutputdir + "run%03d"%(i1) + "/params.txt"))

	ti1, ti3, out = find_common_subset(p,0)
	u = []
	for i in xrange(len(ti3)):
		u.append([ti3[i],i])
	u.sort()
	# EMAN2.display([range(len(u)),[u[i][0] for i in xrange(len(u))]])
	plt.plot(range(len(u)),[u[i][0] for i in xrange(len(u))])

	# import json; f = open("error_curve%03d.json"%plot_errors_between_any_number_of_projections.counter, 'w')
	# json.dump([u[i][0] for i in xrange(len(u))],f); f.close()

	plt.ylabel('Error')
	plt.xlabel('Image index')
	plt.title("Sorted errors between projections")
	import StringIO
	which_projections = StringIO.StringIO()
	which_projections.write("_" + "%.6f"%error_value)
	for p_i in list_of_projection_indices: which_projections.write("_" + "%03d"%p_i)
	for p_i in list_of_projection_indices: which_projections.write("___" + "%03d"%get_already_processed_viper_runs.r_permutation[p_i])

	plt.savefig(mainoutputdir + '/sorted_errors_between_projections' + which_projections.getvalue() + '.png')
	which_projections.close()
	plt.close()


def find_index_of_discontinuity_in_derivative(error_curve_func, list_of_projection_indices, mainoutputdir,
	outlier_percentile):

	import numpy as np

	resolution = 100
	split_point_resolution = 29
	degree_of_the_fitting_polynomial = 1

	data_set_length = len(error_curve_func)

	# split_point_set = np.linspace(0.71,.99,split_point_resolution)
	# split_point_set = np.linspace(0.71,outlier_percentile/100.0,split_point_resolution)
	split_point_set = np.linspace(0.9,outlier_percentile/100.0,split_point_resolution)

	minimum_goodness_of_fit_for_both_lines = 1e20
	optimized_split_point = -1
	for split_point in split_point_set:
		first_line_x = map(int, np.linspace(0,split_point,resolution)*data_set_length)
		first_line_y = np.array([error_curve_func[x] for x in first_line_x])
		first_line_z = np.poly1d( np.polyfit(first_line_x, first_line_y, degree_of_the_fitting_polynomial) )

		second_line_x = map(int, np.linspace(split_point,1, resolution)*data_set_length)
		second_line_y = np.array([error_curve_func[x-1] for x in second_line_x])
		second_line_z = np.poly1d( np.polyfit(second_line_x, second_line_y, degree_of_the_fitting_polynomial) )

		goodness_of_fit_for_both_lines = sum((first_line_z(first_line_x) - first_line_y)**2)
		goodness_of_fit_for_both_lines += sum((second_line_z(second_line_x) - second_line_y)**2)
		# goodness_of_fit_for_both_lines = angle((1,first_line_z[1]), (1,second_line_z[1]))
		if goodness_of_fit_for_both_lines < minimum_goodness_of_fit_for_both_lines:
			minimum_goodness_of_fit_for_both_lines = goodness_of_fit_for_both_lines
			optimized_split_point = split_point

		import matplotlib.pyplot as plt
		# split_point = optimized_split_point
		plt.plot(range(len(error_curve_func)),error_curve_func)

		first_line_x = map(int, np.linspace(0,split_point,resolution)*data_set_length)
		first_line_y = np.array([error_curve_func[x] for x in first_line_x])
		first_line_z = np.poly1d( np.polyfit(first_line_x, first_line_y, degree_of_the_fitting_polynomial) )

		plt.plot(first_line_x,first_line_z(first_line_x))

		second_line_x = map(int, np.linspace(split_point,1, resolution)*data_set_length)
		second_line_y = np.array([error_curve_func[x-1] for x in second_line_x])
		second_line_z = np.poly1d( np.polyfit(second_line_x, second_line_y, degree_of_the_fitting_polynomial) )
		plt.plot(second_line_x,second_line_z(second_line_x))

		import StringIO
		which_projections = StringIO.StringIO()
		which_projections.write("_" + "%.03f__%.6f"%(split_point, goodness_of_fit_for_both_lines))
		for p_i in list_of_projection_indices: which_projections.write("_" + "%03d"%p_i)
		for p_i in list_of_projection_indices: which_projections.write("___" + "%03d"%get_already_processed_viper_runs.r_permutation[p_i])

		plt.title(mainoutputdir + '/sorted_errors' + which_projections.getvalue() + '.png')
		plt.savefig(mainoutputdir + '/sorted_errors' + which_projections.getvalue() + '.png')
		plt.close()

	import matplotlib.pyplot as plt
	split_point = optimized_split_point
	plt.plot(range(len(error_curve_func)),error_curve_func)

	first_line_x = map(int, np.linspace(0,split_point,resolution)*data_set_length)
	first_line_y = np.array([error_curve_func[x] for x in first_line_x])
	first_line_z = np.poly1d( np.polyfit(first_line_x, first_line_y, degree_of_the_fitting_polynomial) )

	plt.plot(first_line_x,first_line_z(first_line_x))

	second_line_x = map(int, np.linspace(split_point,1, resolution)*data_set_length)
	second_line_y = np.array([error_curve_func[x-1] for x in second_line_x])
	second_line_z = np.poly1d( np.polyfit(second_line_x, second_line_y, degree_of_the_fitting_polynomial) )
	plt.plot(second_line_x,second_line_z(second_line_x))

	import StringIO
	which_projections = StringIO.StringIO()
	which_projections.write("_" + "%.03f"%split_point)
	for p_i in list_of_projection_indices: which_projections.write("_" + "%03d"%p_i)
	for p_i in list_of_projection_indices: which_projections.write("___" + "%03d"%get_already_processed_viper_runs.r_permutation[p_i])

	plt.title(mainoutputdir + '/optimized_errors' + which_projections.getvalue() + '.png')
	plt.savefig(mainoutputdir + '/optimized_errors' + which_projections.getvalue() + '.png')
	plt.close()

	if optimized_split_point < 0:
		return -1

	return int(optimized_split_point*data_set_length)


def measure_for_outlier_criterion(criterion_name, masterdir, rviper_iter, list_of_viper_run_indices):

	# main_iterations = ["main" + "%03d" % i for i in range(1, rviper_iter + 1)]
	main_iterations = ["main" + "%03d" % i for i in range(rviper_iter, rviper_iter + 1)]
	mainoutputdir = masterdir + "/" + main_iterations[0] + "/"

	p = []
	for i1 in list_of_viper_run_indices:
		p.append(read_text_row(mainoutputdir + "run%03d" % (i1) + "/params.txt"))
	subset, avg_diff_per_image, outp = find_common_subset(p, 0)

	avg_diff_per_image.sort()
	x1 = len(avg_diff_per_image)
	y1 = avg_diff_per_image[-1]
	
	if y1 <= ANGLE_ERROR_THRESHOLD:
		return TRIPLET_WITH_ANGLE_ERROR_LESS_THAN_THRESHOLD_HAS_BEEN_FOUND

	if criterion_name == "80th percentile":
		return avg_diff_per_image[int(x1*PERCENT_THRESHOLD_X)]/y1
	elif criterion_name == "fastest increase in the last quartile":
		for k in range(5,6):
			avg_diff_per_image_diff = [x - avg_diff_per_image[i - k] for i, x in enumerate(avg_diff_per_image)][k:]
			
			avg_diff_per_image_diff_max = max(avg_diff_per_image_diff)
			avg_diff_per_image_diff_max_normalized = max(avg_diff_per_image_diff)/y1
			
			if avg_diff_per_image_diff.index(avg_diff_per_image_diff_max) >= int(x1*0.75):
				return avg_diff_per_image_diff_max_normalized
			return 0.0
	else:
		print "Error, no criterion name is specified!"
		mpi_finalize()
		sys.exit()



def found_outliers(list_of_projection_indices, outlier_percentile, rviper_iter, masterdir,  bdb_stack_location,
	outlier_index_threshold_method):
	
	# sxheader.py bdb:nj  --consecutive  --params=OID
	import numpy as np

	mainoutputdir = masterdir + "/main%03d/"%(rviper_iter)

	# if this data analysis step was already performed in the past then return
	for check_run in list_of_projection_indices:
		if not (os.path.exists(mainoutputdir + "/run%03d"%(check_run) + "/rotated_reduced_params.txt")):
			break
	else:
		return

	print "identify_outliers"
	projs = []
	for i1 in list_of_projection_indices:
		projs.append(read_text_row(mainoutputdir + "run%03d"%(i1) + "/params.txt"))

	# ti1, ti3, out = find_common_subset(projs, 1.0)
	subset, avg_diff_per_image, rotated_params = find_common_subset(projs, target_threshold = 0)
	# subset, avg_diff_per_image, rotated_params = find_common_subset(projs, target_threshold = 1.0)
	error_values_and_indices = []
	for i in xrange(len(avg_diff_per_image)):
		error_values_and_indices.append([avg_diff_per_image[i], i])
	del subset, avg_diff_per_image

	error_values_and_indices.sort()

	if outlier_index_threshold_method == "discontinuity_in_derivative":
		outlier_index_threshold = find_index_of_discontinuity_in_derivative([i[0] for i in error_values_and_indices],
		list_of_projection_indices, mainoutputdir, outlier_percentile)
	elif outlier_index_threshold_method == "percentile":
		outlier_index_threshold = outlier_percentile * (len(error_values_and_indices) - 1)/ 100.0
	elif outlier_index_threshold_method == "angle_measure":
		error_values = [i[0] for i in error_values_and_indices]
		outlier_index_threshold = min(range(len(error_values)), key=lambda i: abs(error_values[i]-outlier_percentile))
	elif outlier_index_threshold_method == "use all images":
		outlier_index_threshold = len(error_values_and_indices)


	
	index_keep_images = [i[1] for i in error_values_and_indices[:outlier_index_threshold]]
	index_outliers = [i[1] for i in error_values_and_indices[outlier_index_threshold:]]

	# print "error_values_and_indices: %f"%error_values_and_indices
	print "index_outliers: ", index_outliers

	import copy
	reversed_sorted_index_outliers = copy.deepcopy(index_outliers)
	reversed_sorted_index_outliers.sort(reverse=True)

	for k in xrange(len(projs)):
		for l in reversed_sorted_index_outliers:
			del rotated_params[k][l]

	index_outliers.sort()
	index_keep_images.sort()

	write_text_file(index_outliers, mainoutputdir + "this_iteration_index_outliers.txt")
	write_text_file(index_keep_images, mainoutputdir + "this_iteration_index_keep_images.txt")

	#if len(index_outliers) < 3:
		#return False

	if len(index_outliers) > 0:
		cmd = "{} {} {} {}".format("e2bdb.py ", bdb_stack_location + "_%03d"%(rviper_iter - 1), "--makevstack=" + bdb_stack_location + "_outliers_%03d"%(rviper_iter), "--list=" + mainoutputdir  +  "this_iteration_index_outliers.txt")
		cmdexecute(cmd)
	cmd = "{} {} {} {}".format("e2bdb.py ", bdb_stack_location + "_%03d"%(rviper_iter - 1), "--makevstack=" + bdb_stack_location + "_%03d"%(rviper_iter), "--list=" + mainoutputdir +  "this_iteration_index_keep_images.txt")
	cmdexecute(cmd)
	dat = EMData.read_images(bdb_stack_location + "_%03d"%(rviper_iter - 1))

	write_text_file([dat[i].get_attr("original_image_index")  for i in index_outliers],mainoutputdir + "index_outliers.txt")
	write_text_file([dat[i].get_attr("original_image_index")  for i in index_keep_images],mainoutputdir + "index_keep_images.txt")

	print "index_outliers:: " + str(index_outliers)

	# write rotated param files
	for i1 in range(len(list_of_projection_indices)):
		write_text_row(rotated_params[i1], mainoutputdir + "run%03d"%(list_of_projection_indices[i1]) + "/rotated_reduced_params.txt")

	return True


def calculate_volumes_after_rotation_and_save_them(ali3d_options, rviper_iter, masterdir, bdb_stack_location, mpi_rank, mpi_size,
												   no_of_viper_runs_analyzed_together, no_of_viper_runs_analyzed_together_from_user_options, mpi_comm = -1):
	
	# This function takes into account the case in which there are more processors than images

	if mpi_comm == -1:
		mpi_comm = MPI_COMM_WORLD

	# some arguments are for debugging purposes

	mainoutputdir = masterdir + "/main%03d/"%(rviper_iter)

	# list_of_projection_indices_used_for_outlier_elimination = map(int, read_text_file(mainoutputdir + "/list_of_viper_runs_included_in_outlier_elimination.txt"))
	import json; f = open(mainoutputdir + "list_of_viper_runs_included_in_outlier_elimination.json", 'r')
	list_of_independent_viper_run_indices_used_for_outlier_elimination  = json.load(f); f.close()

	if len(list_of_independent_viper_run_indices_used_for_outlier_elimination)==0:
		print "Error: len(list_of_independent_viper_run_indices_used_for_outlier_elimination)==0"
		mpi_finalize()
		sys.exit()

	# if this data analysis step was already performed in the past then return
	# for future changes make sure that the file checked is the last one to be processed !!!
	# if(os.path.exists(mainoutputdir + "/run%03d"%(no_of_viper_runs_analyzed_together - 1) + "/rotated_volume.hdf")):
	# check_last_run = max(get_latest_directory_increment_value(mainoutputdir, "run", start_value=0), no_of_viper_runs_analyzed_together_from_user_options)
	# if(os.path.exists(mainoutputdir + "/run%03d"%(check_last_run) + "/rotated_volume.hdf")):
	# 	return

	# if this data analysis step was already performed in the past then return
	for check_run in list_of_independent_viper_run_indices_used_for_outlier_elimination:
		if not (os.path.exists(mainoutputdir + "/run%03d"%(check_run) + "/rotated_volume.hdf")):
			break
	else:
		return

	partstack = []
	# for i1 in range(0,no_of_viper_runs_analyzed_together):
	for i1 in list_of_independent_viper_run_indices_used_for_outlier_elimination:
		partstack.append(mainoutputdir + "run%03d"%(i1) + "/rotated_reduced_params.txt")
	partids_file_name = mainoutputdir + "this_iteration_index_keep_images.txt"

	lpartids = map(int, read_text_file(partids_file_name) )
	n_projs = len(lpartids)


	if (mpi_size > n_projs):
		# if there are more processors than images
		working = int(not(mpi_rank < n_projs))
		mpi_subcomm = mpi_comm_split(mpi_comm, working,  mpi_rank - working*n_projs)
		mpi_subsize = mpi_comm_size(mpi_subcomm)
		mpi_subrank = mpi_comm_rank(mpi_subcomm)
		if (mpi_rank < n_projs):

			# for i in xrange(no_of_viper_runs_analyzed_together):
			for idx, i in enumerate(list_of_independent_viper_run_indices_used_for_outlier_elimination):
				projdata = getindexdata(bdb_stack_location + "_%03d"%(rviper_iter - 1), partids_file_name, partstack[idx], mpi_rank, mpi_subsize)
				vol = do_volume(projdata, ali3d_options, 0, mpi_comm = mpi_subcomm)
				del projdata
				if( mpi_rank == 0):
					vol.write_image(mainoutputdir + "/run%03d"%(i) + "/rotated_volume.hdf")
					line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " => "
					print line  + "Generated rec_ref_volume_run #%01d \n"%i
				del vol

		mpi_barrier(mpi_comm)
	else:
		for idx, i in enumerate(list_of_independent_viper_run_indices_used_for_outlier_elimination):
			projdata = getindexdata(bdb_stack_location + "_%03d"%(rviper_iter - 1), partids_file_name, partstack[idx], mpi_rank, mpi_size)
			vol = do_volume(projdata, ali3d_options, 0, mpi_comm = mpi_comm)
			del projdata
			if( mpi_rank == 0):
				vol.write_image(mainoutputdir + "/run%03d"%(i) + "/rotated_volume.hdf")
				line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " => "
				print line + "Generated rec_ref_volume_run #%01d"%i
			del vol

	if( mpi_rank == 0):
		# Align all rotated volumes, calculate their average and save as an overall result
		from utilities import get_params3D, set_params3D, get_im, model_circle
		from statistics import ave_var
		from applications import ali_vol
		# vls = [None]*no_of_viper_runs_analyzed_together
		vls = [None]*len(list_of_independent_viper_run_indices_used_for_outlier_elimination)
		# for i in xrange(no_of_viper_runs_analyzed_together):
		for idx, i in enumerate(list_of_independent_viper_run_indices_used_for_outlier_elimination):
			vls[idx] = get_im(mainoutputdir + "/run%03d"%(i) + "/rotated_volume.hdf")
			set_params3D(vls[idx],[0.,0.,0.,0.,0.,0.,0,1.0])
		asa,sas = ave_var(vls)
		# do the alignment
		nx = asa.get_xsize()
		radius = nx/2 - .5
		st = Util.infomask(asa*asa, model_circle(radius,nx,nx,nx), True)
		goal = st[0]
		going = True
		while(going):
			set_params3D(asa,[0.,0.,0.,0.,0.,0.,0,1.0])
			# for i in xrange(no_of_viper_runs_analyzed_together):
			for idx, i in enumerate(list_of_independent_viper_run_indices_used_for_outlier_elimination):
				o = ali_vol(vls[idx],asa,7.0,5.,radius)  # range of angles and shifts, maybe should be adjusted
				p = get_params3D(o)
				del o
				set_params3D(vls[idx],p)
			asa,sas = ave_var(vls)
			st = Util.infomask(asa*asa, model_circle(radius,nx,nx,nx), True)
			if(st[0] > goal):  goal = st[0]
			else:  going = False
		# over and out
		asa.write_image(mainoutputdir + "/average_volume.hdf")
		sas.write_image(mainoutputdir + "/variance_volume.hdf")
	return



def get_already_processed_viper_runs(run_get_already_processed_viper_runs):

	import random

	if run_get_already_processed_viper_runs:
		location_location = "/Users/hvoicu/Analysis/rrviper/particle__PIC_ISAC_g1_clean/0001__sim_r_viper_pool_001/"
		location_location = "/Users/hvoicu/Analysis/rrviper/particle__sp_MED_isac_clean_v1/0001__sim_r_viper_pool_001/"
		
		if "counter" not in get_already_processed_viper_runs.__dict__:
			# function needs to be called once before being used !
			get_already_processed_viper_runs.counter = -2
	
			import os
			path, dirs, files = os.walk(location_location).next()
			# dirs = filter(lambda x:'run' in x, dirs)
			import re
			dirs = filter(lambda x:re.search('run\d\d\d$', x), dirs)
			get_already_processed_viper_runs.r_permutation = range(len(dirs))
			random.shuffle(get_already_processed_viper_runs.r_permutation)
			print str(get_already_processed_viper_runs.r_permutation)
		get_already_processed_viper_runs.counter += 1
		print "get_already_processed_viper_runs.counter: " + str(get_already_processed_viper_runs.counter)
		# if get_already_processed_viper_runs.counter > 9:
		if get_already_processed_viper_runs.counter > (MAXIMUM_NO_OF_VIPER_RUNS_ANALYZED_TOGETHER - 1):
			print "get_already_processed_viper_runs.counter > 9"
			mpi_finalize()
			sys.exit()
	
		return location_location + "run%03d"%get_already_processed_viper_runs.r_permutation[get_already_processed_viper_runs.counter]
	else:
		get_already_processed_viper_runs.r_permutation = [0]*20


def main():

	from logger import Logger, BaseLogger_Files
	import user_functions
	from optparse import OptionParser
	from global_def import SPARXVERSION
	from EMAN2 import EMData



	main_node = 0
	mpi_init(0, [])
	mpi_comm = MPI_COMM_WORLD
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	mpi_size = mpi_comm_size(MPI_COMM_WORLD)	# Total number of processes, passed by --np option.

	# mpi_barrier(mpi_comm)
	# from mpi import mpi_finalize
	# mpi_finalize()
	# print "mpi finalize"
	# from sys import exit
	# exit()

	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack  [output_directory]  [initial_volume]  --ir=inner_radius --ou=outer_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translational_search_step  --delta=angular_step --an=angular_neighborhood  --center=center_type --maxit1=max_iter1 --maxit2=max_iter2 --L2threshold=0.1  --fl --aa --ref_a=S --sym=c1"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ir",		type= "int",   default= 1,                  help="inner radius for rotational correlation > 0 (set to 1)")
	parser.add_option("--ou",       type= "int",   default= -1,                 help="outer radius for rotational correlation < int(nx/2)-1 (set to the radius of the particle)")
	parser.add_option("--rs",       type= "int",   default= 1,                  help="step between rings in rotational correlation >0  (set to 1)" ) 
	parser.add_option("--xr",       type="string", default= "0",                help="range for translation search in x direction, search is +/xr (default 0)")
	parser.add_option("--yr",       type="string", default= "-1",               help="range for translation search in y direction, search is +/yr (default = same as xr)")
	parser.add_option("--ts",       type="string", default= "1",                help="step size of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional")
	parser.add_option("--delta",    type="string", default= "2",                help="angular step of reference projections (default 2)")
	#parser.add_option("--an",       type="string", default= "-1",              help="angular neighborhood for local searches (phi and theta)")
	parser.add_option("--center",   type="float",  default= -1,                 help="-1: average shift method; 0: no centering; 1: center of gravity (default=-1)")
	parser.add_option("--maxit1",    type="float",  default= 400,               help="maximum number of iterations performed for the GA part (set to 400) ")
	parser.add_option("--maxit2",    type="float",  default= 50,                help="maximum number of iterations performed for the finishing up part (set to 50) ")
	parser.add_option("--L2threshold", type="float",  default= 0.03,            help="Stopping criterion of GA given as a maximum relative dispersion of L2 norms (set to 0.03) ")
	parser.add_option("--doga",     type="float",  default= 0.1,                help="do GA when fraction of orientation changes less than 1.0 degrees is at least doga (default=0.1)")
	parser.add_option("--n_shc_runs",    type="int",    default= 3,            help="number of quasi-independent runs (shc) (default=3)")
	parser.add_option("--n_rv_runs",       type= "int",   default= 30,          help="number of r_viper runs")
	parser.add_option("--n_v_runs",       type= "int",   default= 3,            help="number of viper runs for each r_viper cycle")
	parser.add_option("--outlier_percentile",     type="float",    default= 95, help="percentile above which outliers are removed every iteration")
	parser.add_option("--iteration_start",     type="int",    default= 0,       help="starting iteration for rviper, 0 means go to the most recent one (default).")
	#parser.add_option("--CTF",      action="store_true", default=False,        help="NOT IMPLEMENTED Consider CTF correction during the alignment ")
	#parser.add_option("--snr",      type="float",  default= 1.0,               help="Signal-to-Noise Ratio of the data (default 1.0)")
	parser.add_option("--ref_a",    type="string", default= "S",                help="method for generating the quasi-uniformly distributed projection directions (default S)")
	parser.add_option("--sym",      type="string", default= "c1",               help="symmetry of the refined structure")
	parser.add_option("--function", type="string", default="ref_ali3d",         help="name of the reference preparation function (ref_ali3d by default)")
	parser.add_option("--npad",     type="int",    default= 2,                  help="padding size for 3D reconstruction (default=2)")

	#options introduced for the do_volume function
	parser.add_option("--fl",      type="float",  default=0.12,    help="cut-off frequency of hyperbolic tangent low-pass Fourier filte (default 0.12)")
	parser.add_option("--aa",      type="float",  default=0.1,    help="fall-off of hyperbolic tangent low-pass Fourier filter (default 0.1)")
	parser.add_option("--pwreference",      type="string",  default="",    help="text file with a reference power spectrum (default no power spectrum adjustment)")
	parser.add_option("--mask3D",      type="string",  default=None,    help="3D mask file (default a sphere  WHAT RADIUS??)")
	parser.add_option("--moon_elimination",      type="string",  default=None,    help="mass in KDa and resolution in px/A separated by comma, no space")


	parser.add_option("--my_random_seed",      type="int",  default=123,    help="random seed, default value: 123")
	parser.add_option("--criterion_name",      type="string",  default="80th percentile",    help="default: 80th percentile, other options:/fastest increase in the last quartile/")
	parser.add_option("--outlier_index_threshold_method",      type="string",  default="discontinuity_in_derivative",    help="default: discontinuity_in_derivative, other options:/percentile/angle_measure")
	
	parser.add_option("--run_get_already_processed_viper_runs", action="store_true", dest="run_get_already_processed_viper_runs", default=False)
	
	parser.add_option("--use_latest_master_directory", action="store_true", dest="use_latest_master_directory", default=False)
	
	(options, args) = parser.parse_args(sys.argv[1:])

	options.CTF = False
	options.snr = 1.0
	options.an = -1

	if options.moon_elimination==None:
		options.moon_elimination = []
	else:
		options.moon_elimination = map(float, options.moon_elimination.split(","))

	my_random_seed = options.my_random_seed
	criterion_name = options.criterion_name
	outlier_index_threshold_method = options.outlier_index_threshold_method
	use_latest_master_directory = options.use_latest_master_directory
	iteration_start_default = options.iteration_start
	number_of_rrr_viper_runs = options.n_rv_runs
	no_of_viper_runs_analyzed_together_from_user_options = options.n_v_runs
	no_of_shc_runs_analyzed_together = options.n_shc_runs 
	outlier_percentile = options.outlier_percentile 
	
	run_get_already_processed_viper_runs = options.run_get_already_processed_viper_runs
	get_already_processed_viper_runs(run_get_already_processed_viper_runs)

	import random
	random.seed(my_random_seed)

	if len(args) < 1 or len(args) > 3:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
		return 1

	if len(args) > 2:
		ref_vol = get_im(args[2])
	else:
		ref_vol = None

	masterdir = ""
	bdb_stack_location = ""
	if len(args) == 2:
		masterdir = args[1]
		if masterdir[-1] != "/":
			masterdir += "/"
	elif len(args) == 1:
		if use_latest_master_directory:
			all_dirs = [d for d in os.listdir(".") if os.path.isdir(d)]
			import re; r = re.compile("^master.*$")
			all_dirs = filter(r.match, all_dirs)
			if len(all_dirs)>0:
				# all_dirs = max(all_dirs, key=os.path.getctime)
				masterdir = max(all_dirs, key=os.path.getmtime)
				masterdir += "/"

	log = Logger(BaseLogger_Files())

	error_status = 0	
	if mpi_size % no_of_shc_runs_analyzed_together != 0:
		ERROR('Number of processes needs to be a multiple of total number of runs. '
		'Total quasi-independent runs by default are 3, you can change it by specifying '
		'--nruns option. Also, to improve communication time it is recommended that '
		'the number of processes divided by the number of quasi-independent runs is a power '
		'of 2 (e.g. 2, 4, 8 or 16 depending on how many physical cores each node has).', 'sxviper', 1)
		error_status = 1
	if_error_all_processes_quit_program(error_status)

	#Create folder for all results or check if there is one created already
	if(myid == main_node):
		#cmd = "{}".format("Rmycounter ccc")
		#cmdexecute(cmd)

		if( masterdir == ""):
			timestring = strftime("%Y_%m_%d__%H_%M_%S/", localtime())
			masterdir = "master"+timestring

			cmd = "{} {}".format("mkdir", masterdir)
			cmdexecute(cmd)
		if os.path.exists(masterdir):
			if ':' in args[0]:
				bdb_stack_location = args[0].split(":")[0] + ":" + masterdir + args[0].split(":")[1]
				org_stack_location = args[0]

				if(not os.path.exists(os.path.join(masterdir,"EMAN2DB/"))):
					# cmd = "{} {}".format("cp -rp EMAN2DB", masterdir, "EMAN2DB/")
					# cmdexecute(cmd)
					cmd = "{} {} {}".format("e2bdb.py", org_stack_location,"--makevstack=" + bdb_stack_location + "_000")
					cmdexecute(cmd)

					try:
						os_return_value = os.system("sxheader.py  " + bdb_stack_location + "_000 --print  --params=original_image_index")
					except:
						pass
					if os_return_value != 0:
						cmd = "{} {}".format("sxheader.py  ", bdb_stack_location + "_000 --consecutive  --params=original_image_index")
						cmdexecute(cmd)
			else:
				filename = os.path.basename(args[0])
				bdb_stack_location = "bdb:" + masterdir + os.path.splitext(filename)[0]
				if(not os.path.exists(os.path.join(masterdir,"EMAN2DB/"))):
					cmd = "{} {} {}".format("sxcpy.py  ", args[0], bdb_stack_location + "_000")
					cmdexecute(cmd)

					try:
						os_return_value = os.system("sxheader.py  " + bdb_stack_location + "_000 --print  --params=original_image_index")
					except:
						pass
					if os_return_value != 0:
						cmd = "{} {}".format("sxheader.py  ", bdb_stack_location + "_000 --consecutive  --params=original_image_index")
						cmdexecute(cmd)
				else:
					ERROR('Conflicting information: EMAN2DB exists, but provided *.hdf file', "sxrviper", 1)
					error_status = 1

		else:
			# os.path.exists(masterdir) does not exist
			ERROR('Output directory does not exist, please change the name and restart the program', "sxrviper", 1)
			error_status = 1

	if_error_all_processes_quit_program(error_status)

	# send masterdir to all processes
	dir_len  = len(masterdir)*int(myid == main_node)
	dir_len = mpi_bcast(dir_len,1,MPI_INT,0,MPI_COMM_WORLD)[0]
	masterdir = mpi_bcast(masterdir,dir_len,MPI_CHAR,main_node,MPI_COMM_WORLD)
	masterdir = string.join(masterdir,"")
	if masterdir[-1] != "/":
		masterdir += "/"
		
	global_def.LOGFILE =  os.path.join(masterdir, global_def.LOGFILE)
	print_program_start_information()
	

	# mpi_barrier(mpi_comm)
	# from mpi import mpi_finalize
	# mpi_finalize()
	# print "mpi finalize"
	# from sys import exit
	# exit()
		
	
	# send bdb_stack_location to all processes
	dir_len  = len(bdb_stack_location)*int(myid == main_node)
	dir_len = mpi_bcast(dir_len,1,MPI_INT,0,MPI_COMM_WORLD)[0]
	bdb_stack_location = mpi_bcast(bdb_stack_location,dir_len,MPI_CHAR,main_node,MPI_COMM_WORLD)
	bdb_stack_location = string.join(bdb_stack_location,"")

	iteration_start = get_latest_directory_increment_value(masterdir, "main")

	if (myid == main_node):
		if (iteration_start < iteration_start_default):
			ERROR('Starting iteration provided is greater than last iteration performed. Quiting program', 'sxviper', 1)
			error_status = 1
	if iteration_start_default!=0:
		iteration_start = iteration_start_default
	if (myid == main_node):
		if (number_of_rrr_viper_runs < iteration_start):
			ERROR('Please provide number of rviper runs (--n_rv_runs) greater than number of iterations already performed.', 'sxviper', 1)
			error_status = 1

	if_error_all_processes_quit_program(error_status)

	for rviper_iter in range(iteration_start, number_of_rrr_viper_runs + 1):
		if(myid == main_node):
			all_projs = EMData.read_images(bdb_stack_location + "_%03d"%(rviper_iter - 1))
			print "XXXXXXXXXXXXXXXXX"
			print "Number of projections (in loop): " + str(len(all_projs))
			print "XXXXXXXXXXXXXXXXX"
			subset = range(len(all_projs))
		else:
			all_projs = None
			subset = None

		runs_iter = get_latest_directory_increment_value(masterdir + "main%03d"%rviper_iter, "/run", start_value=0) - 1
		no_of_viper_runs_analyzed_together = max(runs_iter + 2, no_of_viper_runs_analyzed_together_from_user_options)

		first_time_entering_the_loop_need_to_do_full_check_up = True
		while True:
			runs_iter += 1

			if not first_time_entering_the_loop_need_to_do_full_check_up:
				if runs_iter >= no_of_viper_runs_analyzed_together:
					break
			first_time_entering_the_loop_need_to_do_full_check_up = False

			this_run_is_NOT_complete = 0
			if (myid == main_node):
				independent_run_dir = masterdir + '/main%03d/run%03d/'%(rviper_iter, runs_iter)
				if run_get_already_processed_viper_runs:
					cmd = "{} {}".format("mkdir -p", masterdir + '/main%03d/'%(rviper_iter)); cmdexecute(cmd)
					cmd = "{} {}".format("rm -rf", independent_run_dir); cmdexecute(cmd)
					cmd = "{} {}".format("cp -r", get_already_processed_viper_runs() + " " +  independent_run_dir); cmdexecute(cmd)
				
				if os.path.exists(independent_run_dir + "log.txt") and (string_found_in_file("Finish VIPER2", independent_run_dir + "log.txt")):
					this_run_is_NOT_complete = 0
				else:
					this_run_is_NOT_complete = 1
					cmd = "{} {}".format("rm -rf", independent_run_dir); cmdexecute(cmd)
					cmd = "{} {}".format("mkdir -p", independent_run_dir); cmdexecute(cmd)

				this_run_is_NOT_complete = mpi_bcast(this_run_is_NOT_complete,1,MPI_INT,main_node,MPI_COMM_WORLD)[0]
				dir_len = len(independent_run_dir)
				dir_len = mpi_bcast(dir_len,1,MPI_INT,main_node,MPI_COMM_WORLD)[0]
				independent_run_dir = mpi_bcast(independent_run_dir,dir_len,MPI_CHAR,main_node,MPI_COMM_WORLD)
				independent_run_dir = string.join(independent_run_dir,"")
			else:
				this_run_is_NOT_complete = mpi_bcast(this_run_is_NOT_complete,1,MPI_INT,main_node,MPI_COMM_WORLD)[0]
				dir_len = 0
				independent_run_dir = ""
				dir_len = mpi_bcast(dir_len,1,MPI_INT,main_node,MPI_COMM_WORLD)[0]
				independent_run_dir = mpi_bcast(independent_run_dir,dir_len,MPI_CHAR,main_node,MPI_COMM_WORLD)
				independent_run_dir = string.join(independent_run_dir,"")

			if this_run_is_NOT_complete:
				mpi_barrier(MPI_COMM_WORLD)

				if independent_run_dir[-1] != "/":
					independent_run_dir += "/"

				log.prefix = independent_run_dir

				options.user_func = user_functions.factory[options.function]

				# for debugging purposes
				#if (myid == main_node):
					#cmd = "{} {}".format("cp ~/log.txt ", independent_run_dir)
					#cmdexecute(cmd)
					#cmd = "{} {}{}".format("cp ~/paramdir/params$(mycounter ccc).txt ", independent_run_dir, "param%03d.txt"%runs_iter)
					#cmd = "{} {}{}".format("cp ~/paramdir/params$(mycounter ccc).txt ", independent_run_dir, "params.txt")
					#cmdexecute(cmd)

				if (myid == main_node):
					store_value_of_simple_vars_in_json_file(locals(), exclude_list_of_vars=["usage"], 
						vars_that_will_show_only_size = ["subset"])
					store_value_of_simple_vars_in_json_file(options.__dict__, write_or_append='a')
				
				# mpi_barrier(mpi_comm)
				# from mpi import mpi_finalize
				# mpi_finalize()
				# print "mpi finalize"
				# from sys import exit
				# exit()
				
				out_params, out_vol, out_peaks = multi_shc(all_projs, subset, no_of_shc_runs_analyzed_together, options,
				mpi_comm=mpi_comm, log=log, ref_vol=ref_vol)

				# end of: if this_run_is_NOT_complete:

			if runs_iter >= (no_of_viper_runs_analyzed_together_from_user_options - 1):
				increment_for_current_iteration = identify_outliers(myid, main_node, rviper_iter,
				no_of_viper_runs_analyzed_together, no_of_viper_runs_analyzed_together_from_user_options, masterdir,
				bdb_stack_location, outlier_percentile, criterion_name, outlier_index_threshold_method)
				
				if increment_for_current_iteration == MUST_END_PROGRAM_THIS_ITERATION:
					break
				
				no_of_viper_runs_analyzed_together += increment_for_current_iteration

		# end of independent viper loop

		calculate_volumes_after_rotation_and_save_them(options, rviper_iter, masterdir, bdb_stack_location, myid,
		mpi_size, no_of_viper_runs_analyzed_together, no_of_viper_runs_analyzed_together_from_user_options)
		
		if increment_for_current_iteration == MUST_END_PROGRAM_THIS_ITERATION:
			break

	# end of R viper loop

	#mpi_finalize()
	#sys.exit()

	
	mpi_finalize()


if __name__=="__main__":
	main()

