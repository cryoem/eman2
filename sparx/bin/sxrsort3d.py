#!/usr/bin/env python
#
#
#  08/26/2016
#  New version of sort3D.
#  
from __future__ import print_function
import  os
import  sys
import  types
import  global_def
from    global_def import *
from    optparse   import OptionParser
from    sparx      import *
from    EMAN2      import *
from    numpy      import array
from    logger     import Logger, BaseLogger_Files

from mpi   	import  *
from math  	import  *
from random import *

import os
import sys
import subprocess
import time
import string
import json
from   sys 	import exit
from   time import localtime, strftime, sleep

global Tracker, Blockdata

# ------------------------------------------------------------------------------------
mpi_init(0, [])
nproc    = mpi_comm_size(MPI_COMM_WORLD)
myid     = mpi_comm_rank(MPI_COMM_WORLD)


Blockdata = {}
#  MPI stuff
Blockdata["nproc"]              = nproc
Blockdata["myid"]               = myid
Blockdata["main_node"]          = 0

Blockdata["shared_comm"]                    = mpi_comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED,  0, MPI_INFO_NULL)
Blockdata["myid_on_node"]                   = mpi_comm_rank(Blockdata["shared_comm"])
Blockdata["no_of_processes_per_group"]      = mpi_comm_size(Blockdata["shared_comm"])
masters_from_groups_vs_everything_else_comm = mpi_comm_split(MPI_COMM_WORLD, Blockdata["main_node"] == Blockdata["myid_on_node"], Blockdata["myid_on_node"])
Blockdata["color"], Blockdata["no_of_groups"], balanced_processor_load_on_nodes = get_colors_and_subsets(Blockdata["main_node"], MPI_COMM_WORLD, Blockdata["myid"], \
         Blockdata["shared_comm"], Blockdata["myid_on_node"], masters_from_groups_vs_everything_else_comm)
#  We need two nodes for processing of volumes
Blockdata["node_volume"] = [Blockdata["no_of_groups"]-3, Blockdata["no_of_groups"]-2, Blockdata["no_of_groups"]-1]  # For 3D stuff take three last nodes
#  We need two CPUs for processing of volumes, they are taken to be main CPUs on each volume
#  We have to send the two myids to all nodes so we can identify main nodes on two selected groups.
Blockdata["nodes"] = [Blockdata["node_volume"][0]*Blockdata["no_of_processes_per_group"],Blockdata["node_volume"][1]*Blockdata["no_of_processes_per_group"], \
     Blockdata["node_volume"][2]*Blockdata["no_of_processes_per_group"]]
# End of Blockdata: sorting requires at least three nodes, and the used number of nodes be integer times of three

####--
def do_EQKmeans_nways_clustering(workdir, initial_partids, params, sort_res, log_main):
	global Tracker, Blockdata
	# Using EQK->two_way_comparision->Kmeans to split a dataset into clusters till the number of 
	# the unaccounted is less than the minimum size of a cluster
	# input:  initial_partids
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] == Blockdata["main_node"]:
		msg = "----->>>do_EQKmeans_nways_clustering<<<------"
		print(line, msg)
		log_main.add(msg)
		Tracker["unaccounted_list"] = read_text_file(initial_partids, -1) # read all columns
		if len(Tracker["unaccounted_list"])>1: # either one column or two columns
			Tracker["unaccounted_list"] = Tracker["unaccounted_list"][1] # two column entries
		else:
			Tracker["unaccounted_list"] = Tracker["unaccounted_list"][0] # only one column
	else:   Tracker["unaccounted_list"] = 0
	Tracker["unaccounted_list"] = wrap_mpi_bcast(Tracker["unaccounted_list"], Blockdata["main_node"], MPI_COMM_WORLD)
	
	generation		            = 0
	Tracker["total_stack"]      = len(Tracker["unaccounted_list"])
	Tracker["number_of_groups"] = get_number_of_groups(Tracker["total_stack"],Tracker["number_of_images_per_group"])
	partids                     = initial_partids
	
	if Tracker["number_of_groups"]>1: # In case the number of the input particles is small
	
		while Tracker["number_of_groups"] >1:
			if Blockdata["myid"] == Blockdata["main_node"]: Tracker["partition_list"] = []
			else:                                           Tracker["partition_list"] = 0
			
			if Blockdata["myid"] == Blockdata["main_node"]:
				Tracker["directory"]  = os.path.join(workdir, "generation%03d"%generation)
				cmd="{} {}".format("mkdir",Tracker["directory"])
				junk = cmdexecute(cmd)
				log_main.add("-------->>> generation         %5d"%generation)
				log_main.add("number of images per group:  %d"%Tracker["number_of_images_per_group"])
				log_main.add("the initial number of groups:  %d  number of independent runs:  %d"%(Tracker["number_of_groups"], Tracker["constants"]["indep_runs"]))
			else:  Tracker["directory"] = 0
			Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], MPI_COMM_WORLD)
			create_nrandom_lists(partids)
			for indep_run_iter in xrange(0, Tracker["constants"]["indep_runs"]): # N independent runs
				Tracker["indep_run_iter"] = indep_run_iter
				index_file = os.path.join(Tracker["directory"],"independent_index_%03d.txt"%indep_run_iter)
				line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
				if Blockdata["myid"] == Blockdata["main_node"]:
					Tracker["directory"]  = os.path.join(workdir,"generation%03d"%generation, "EQKmeans_%03d"%Tracker["indep_run_iter"])
					msg =  "indep_run_iter %d"%indep_run_iter
					print(line, msg)
					log_main.add(msg)
					cmd="{} {}".format("mkdir",Tracker["directory"])
					junk = cmdexecute(cmd)
				Tracker        = wrap_mpi_bcast(Tracker, Blockdata["main_node"], MPI_COMM_WORLD)


				tmp_final_list = mref_ali3d_EQ_Kmeans(index_file, params, Tracker["clean_volumes"])


				Tracker["directory"] =  os.path.join(workdir, "generation%03d"%generation)
				if Blockdata["myid"] == Blockdata["main_node"]: Tracker["partition_list"].append(Tracker["partition"])
			Tracker["partition_list"] = wrap_mpi_bcast(Tracker["partition_list"], Blockdata["main_node"], MPI_COMM_WORLD)
			
			do_two_way_comparison_over_nindepruns(log_main)
			Tracker["directory"]  = os.path.join(workdir,"generation%03d"%generation,"Kmeans")
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			if Blockdata["myid"] == Blockdata["main_node"]:
				msg ="------>>>Kmeans<<<-------"
				print(line, msg)
				log_main.add(msg)
				msg = "initial number of groups: %d"%Tracker["number_of_groups"] 
				print(line, msg)
				log_main.add(msg)
				cmd="{} {}".format("mkdir",Tracker["directory"])
				junk = cmdexecute(cmd)
			Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], MPI_COMM_WORLD)
			final_list, num_in_groups = mref_ali3d_Kmeans_remove_small_groups(Tracker["Accounted_on_disk"],  \
			os.path.join(Tracker["constants"]["masterdir"],"refinement_parameters.txt"), Tracker["clean_volumes"])	
			partids    = Tracker["Unaccounted_on_disk"]
		
			if Blockdata["myid"] == Blockdata["main_node"]:
				line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
				for iref in xrange(len(num_in_groups)):
					msg = "group  %d:    %d "%(iref, num_in_groups[iref])
					log_main.add(msg)
					print(line, msg)
				msg = "in total:    %d "%sum(num_in_groups)
				log_main.add(msg)
				print(line, msg)
				msg = "groups of size less than %d are absorbed during Kmeans"%Tracker["constants"]["smallest_group"]
				log_main.add(msg)
				print(line, msg)
				msg = "the final number of groups: %d"%Tracker["number_of_groups"] 
				Tracker["total_stack"]       = len(read_text_row(partids))
				Tracker["number_of_groups"]  = get_number_of_groups(Tracker["total_stack"], Tracker["number_of_images_per_group"])
	
			Tracker     = wrap_mpi_bcast(Tracker, Blockdata["main_node"], MPI_COMM_WORLD)
			generation +=1
			sort_res.append(Tracker["partition"])
	else:
		if Blockdata["myid"] == Blockdata["main_node"]:
			msg ="The total number of particles is less than number_of_particles_per_group"
			print(line, msg)
			log_main.add(msg)
			final_list = []
		else:
			final_list = 0
		final_list  = wrap_mpi_bcast(final_list, Blockdata["main_node"], MPI_COMM_WORLD)
	return final_list # one the result of the last iteration
	 
def create_nrandom_lists(partids):
	# the second column denotes orignal particle IDs
	# the first column is randomized group ID 
	global Tracker, Blockdata
	import copy
	import random
	from   utilities import wrap_mpi_bcast, read_text_file, write_text_file
	
	if Blockdata["myid"] == Blockdata["main_node"]:
		data_list = read_text_file(partids, -1)
		if len(data_list)==1: Tracker["sorting_data_list"]= data_list[0]
		else:                 Tracker["sorting_data_list"]= data_list[1]	
		
		if Tracker["constants"]["seed"] ==- 1: random.seed()
		else:                                  random.seed(Tracker["constants"]["seed"])
		Tracker["indep_runs_list"]  = []
		group_size = len(Tracker["sorting_data_list"])//Tracker["number_of_groups"]
		for index_of_random in xrange(Tracker["constants"]["indep_runs"]):
			particle_dict = {}
			ll = copy.copy(Tracker["sorting_data_list"])
			random.shuffle(ll)
			group_list = []
			for index_of_groups in xrange(Tracker["number_of_groups"]):
				if index_of_groups != Tracker["number_of_groups"]-1:
					for iparticle in ll[index_of_groups*group_size:(index_of_groups+1)*group_size]:
						particle_dict[iparticle] = index_of_groups
						group_list.append(index_of_groups)
				else:
					for iparticle in ll[index_of_groups*group_size:]:
						particle_dict[iparticle] = index_of_groups
						group_list.append(index_of_groups)
			assignment = []
			for im in xrange(len(Tracker["sorting_data_list"])):
				assignment.append([particle_dict[Tracker["sorting_data_list"][im]], Tracker["sorting_data_list"][im]])
			write_text_row(assignment, os.path.join(Tracker["directory"],"independent_index_%03d.txt"%index_of_random))
	else:	
		Tracker["indep_runs_list"]   = 0
		Tracker["sorting_data_list"] = 0
	Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"])
	return
	
def mref_ali3d_Kmeans_remove_small_groups(partids, partstack, clean_volumes = True):
	global Tracker, Blockdata
	# partids      		original particle indices in a text file for all particles for this run, can be a subset of the full dataset
	# parstack     		a text file contains ali3d parameters of particles for this run
	# check existence of small group in each iteration and assigne them to large ones
	from projection import prep_vol
	from logger		import Logger,BaseLogger_Files
	try:	termprec = Tracker["Kmeans_stopercent"]
	except: termprec = 3.0
	log        = Logger()
	log   	   = Logger(BaseLogger_Files())
	log.prefix = Tracker["directory"]+"/"
	try: bckgnoise = Blockdata["bckgnoise"]
	except: Blockdata["bckgnoise"] = None
	
	if Tracker["focus3D"]:
		if(Blockdata["myid"] == Blockdata["main_node"]):  focus = get_shrink_3dmask(Tracker["nxinit"], Tracker["focus3D"])
		else:	                                          focus = model_blank(Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"])
		bcast_EMData_to_all(focus, Blockdata["myid"], Blockdata["main_node"])
		st = Util.infomask(focus, None, True)
		if( st[0] == 0.0 ):  ERROR("Incorrect focus mask, after binarize all values zero", "mref_ali3d_Kmeans_remove_small_groups", 1, Blockdata["myid"])
		focus = prep_vol(focus, 1, 1)
	else: focus = None
	shrinkage  = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	if Tracker["constants"]["interpolation"]=="4nn": 
		 projdata = get_shrink_data_sorting(partids, partstack, return_real = True, preshift = True, apply_mask = True)
	else:
		if Tracker["focus3D"]:  projdata = get_shrink_data_sorting(partids, partstack, return_real = True, preshift = True, apply_mask = True)
		else:  	                projdata = get_shrink_data_sorting(partids, partstack, return_real = False, preshift = True, apply_mask = True)
	#oldangles = []
	npergroup  = [0]*Tracker["number_of_groups"]
	assignment = [-1]*len(projdata)
		
	if(Blockdata["myid"] == Blockdata["main_node"]):  #Let's define total_stack here 
		total_stack = len(read_text_file(partids,-1)[0])
	else: total_stack = 0
	total_stack = bcast_number_to_all(total_stack, Blockdata["main_node"], MPI_COMM_WORLD)
	assert(Tracker["total_stack"] == total_stack)
	
	
	if(Blockdata["myid"] == Blockdata["main_node"]):
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print(line, "mref_ali3d_Kmeans_remove_small_groups")
		log.add("mref_ali3d_Kmeans_remove_small_groups")
		msg = "total_stack:  %d"%Tracker["total_stack"]
		log.add(msg)
		print(line, msg)
		msg = "number_of_groups: %d"%Tracker["number_of_groups"]
		log.add(msg)
		print(line, msg)
		msg = "noctf:  %s"%Tracker["constants"]["noctf"]
		log.add(msg)
		print(line, msg)
		msg = "currrent directory:  %s"%Tracker["directory"]
		log.add(msg)
		print(line, msg)
		msg = "nxinit:   %d    nnxo:   %d"%(Tracker["nxinit"], Tracker["constants"]["nnxo"])
		log.add(msg)
		print(line, msg)
		msg = "symmetry group:  %s"%Tracker["constants"]["symmetry"]
		log.add(msg)
		print(line, msg)
		msg = "the Kmeans stop criterion:  %f"%termprec
		log.add(msg)
		print(line, msg)
		
	Tracker     = wrap_mpi_bcast(Tracker,Blockdata["main_node"])
	nima        = len(projdata)
	disps 		= []
	recvcount 	= []
	for im in xrange(Blockdata["nproc"]):
		if( im == Blockdata["main_node"]):  disps.append(0)
		else:   disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_stack, Blockdata["nproc"] , im)
		recvcount.append(ie - ib)
	
	# Initial reconstruction
	if    Tracker["constants"]["interpolation"]=="trl": do3d_sorting_groups_trl_iter(projdata, iteration = 0)
	elif  Tracker["constants"]["interpolation"]=="4nn": do3d_sorting_groups_4nn_iter(projdata, iteration = 0)
	else: ERROR("Wrong interpolation method for do3d", "mref_ali3d_Kmeans_remove_small_groups", 1, Blockdata["myid"]) 
	
	for im in xrange(len(projdata)): 
		phi, theta, psi, sxs, sys = get_params_proj(projdata[im])
		#oldangles.append([phi, theta, psi])
		nref = projdata[im].get_attr("group")
		if nref != -1:
			assignment[im]   = nref
			npergroup[nref] += 1
			
	npergroup  = mpi_reduce(npergroup, Tracker["number_of_groups"], MPI_INT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)	
	npergroup  = map(int, npergroup)
	
	if(Blockdata["myid"] == Blockdata["main_node"]):
		msg = "Initial assignment of particles:"
		log.add(msg)
		for iref in xrange(Tracker["number_of_groups"]):
			msg = " %5d          %7d              %7d        %d"%(iref+1, npergroup[iref], Tracker["fsc05"][iref], Tracker["fsc143"][iref])
			log.add(msg)
	mpi_barrier(MPI_COMM_WORLD)	
		
	small_group_list     = []
	reassign             = False
	do_Kmeans            = True
	saturated_clustering = 0
	Iter                 = 0
	small_group_found    = 0
	
	while Iter <Tracker["total_number_of_iterations"]:
		if Tracker["constants"]["interpolation"]=="4nn":
			if not  Tracker["focus3D"]:
				for im in xrange(len(projdata)): projdata[im] = fft(projdata[im])
		if Tracker["constants"]["interpolation"]=="trl":
			if Tracker["focus3D"]: 
				for im in xrange(len(projdata)): projdata[im] = fft(projdata[im])
			mpi_barrier(MPI_COMM_WORLD)	
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		peaks =  [-1.0e23 for im in xrange(nima)]	
		for iref in xrange(Tracker["number_of_groups"]):
			get_low_pass_filter_for_Kmeans(iref)
			if(Blockdata["myid"] == Blockdata["main_node"]):
				ref_vol = get_im(os.path.join(Tracker["directory"],"vol_grp%03d_iter%03d.hdf"%(iref, Iter)))
				if Tracker["mask3D"]:
					mask3D   = get_im( Tracker["mask3D"])
					if mask3D.get_xsize()!= ref_vol.get_xsize(): mask3D = get_shrink_3dmask(ref_vol.get_xsize(), Tracker["mask3D"])
					stat     = Util.infomask(ref_vol, mask3D, False)
					ref_vol -= stat[0]
					if stat[1] !=0.0: Util.mul_scalar(ref_vol, 1.0/stat[1])
					Util.mul_img(ref_vol, mask3D)
				nnn = ref_vol.get_xsize()
				if(Tracker["nxinit"] != nnn):
					ref_vol = fdecimate(ref_vol,Tracker["nxinit"],Tracker["nxinit"],Tracker["nxinit"], True, False)
				msg = "Iter  %d  low_pass filter  %f is applied to reference %d "%(Iter, Tracker["lowpass"], iref)
				log.add(msg)	
			else:	ref_vol = model_blank(Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"])
			bcast_EMData_to_all(ref_vol, Blockdata["myid"], Blockdata["main_node"])
			ref_peaks                       = compare_two_images(projdata, ref_vol, Tracker["lowpass"], focus, comparison_method = Tracker["constants"]["comparison_method"])
			for im in xrange(nima):	
				if peaks[im] < ref_peaks[im]:
					peaks[im] = ref_peaks[im]
					projdata[im].set_attr("group", iref)
			mpi_barrier(MPI_COMM_WORLD)
		if(Blockdata["myid"] == Blockdata["main_node"]): print(line, " peaks are computed!")
		####---------------------------------
		nchng     =  0
		npergroup = [0]*Tracker["number_of_groups"]
		for im in xrange(nima):
			iref = projdata[im].get_attr("group")
			npergroup[iref] += 1
			if(iref != assignment[im]): 
				nchng += 1
				assignment[im] = iref
		nchng                = mpi_reduce(nchng, 1, MPI_INT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
		npergroup            = mpi_reduce(npergroup, Tracker["number_of_groups"], MPI_INT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
		npergroup            = map(int, npergroup)
		small_group_found    = 0
		small_group_list     = []
		saturated_clustering = 0
	
		if(Blockdata["myid"] == Blockdata["main_node"]):
			ngroup=[]
			nchng = int(nchng[0])
			precn = 100*float(nchng)/float(total_stack)
			msg   = "Iteration number  %d"%Iter
			log.add(msg)
			msg  = "Number of particles that changed assignments %7d, percentage of total: %5.1f"%(nchng, precn)
			log.add(msg)
			msg  = "Group ID  number of particles   fsc05  fsc143 "
			log.add(msg)
			for iref in xrange(Tracker["number_of_groups"]):
				try: 		fsc143 = Tracker["fsc143"][iref]
				except:		fsc143 = 0.0
				try: 		fsc05  = Tracker["fsc05"][iref]
				except:		fsc05  = 0.0			
				msg = "%5d     %7d    %7d    %7d "%(iref+1, npergroup[iref], fsc05, fsc143)
				log.add(msg)
				ngroup.append(int(npergroup[iref]))
			if(precn <= termprec): saturated_clustering  = 1
			for iref in xrange(Tracker["number_of_groups"]):
				if npergroup[iref] < Tracker["constants"]["smallest_group"]: 
					small_group_found     = 1
					small_group_list.append(iref)
		else: 
			ngroup           = 0
			small_group_list = 0
			
		saturated_clustering        = bcast_number_to_all(saturated_clustering, Blockdata["main_node"], MPI_COMM_WORLD)
		small_group_found           = bcast_number_to_all(small_group_found, Blockdata["main_node"], MPI_COMM_WORLD)
		
		if  small_group_found:
			reassign = True
			small_group_list = wrap_mpi_bcast(small_group_list, Blockdata["main_node"]) 
			res_sort3d = get_sorting_all_params(projdata)
			if(Blockdata["myid"] == Blockdata["main_node"]):
				log.add("reset group IDs")
				group_list, ali3d_params_list   = parsing_sorting_params(partids, res_sort3d)
				new_group_list, Tracker["number_of_groups"] = reassign_groups(partids, group_list, small_group_list, Tracker["number_of_groups"])
				log.add("the updated number_of_groups:   %d"% Tracker["number_of_groups"])
			else:
				new_group_list = 0
			Tracker["number_of_groups"] = bcast_number_to_all(Tracker["number_of_groups"], Blockdata["main_node"], MPI_COMM_WORLD)
			new_group_list                  = wrap_mpi_bcast(new_group_list, Blockdata["main_node"])
			image_start, image_end = MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"])
			new_group_list = new_group_list[image_start : image_end]
			for im in xrange(len(projdata)): projdata[im].set_attr("group", new_group_list[im][0])
			mpi_barrier(MPI_COMM_WORLD)	
			if    Tracker["constants"]["interpolation"]=="trl": 
				if Tracker["focus3D"]: 
					for im in xrange(len(projdata)): projdata[im] = fft(projdata[im])
				mpi_barrier(MPI_COMM_WORLD)
				do3d_sorting_groups_trl_iter(projdata, iteration = Iter+1)
			elif  Tracker["constants"]["interpolation"]=="4nn":
				if not Tracker["focus3D"]:
					for im in xrange(len(projdata)): projdata[im] = fft(projdata[im])
				mpi_barrier(MPI_COMM_WORLD)  
				do3d_sorting_groups_4nn_iter(projdata, iteration = Iter+1)
			else: ERROR("Wrong interpolation method for do3d", "mref_ali3d_Kmeans_remove_small_groups", 1, Blockdata["myid"]) 
			Iter +=1
		else:
			if saturated_clustering:
				do_Kmeans = False 
				break
			else:
				reassign = False
				ngroup = wrap_mpi_bcast(ngroup,Blockdata["main_node"])
				if    Tracker["constants"]["interpolation"]=="trl": 
					if Tracker["focus3D"]: 
						for im in xrange(len(projdata)): projdata[im] = fft(projdata[im])
					mpi_barrier(MPI_COMM_WORLD)
					do3d_sorting_groups_trl_iter(projdata, iteration = Iter+1)	
				elif  Tracker["constants"]["interpolation"]=="4nn":
					if not Tracker["focus3D"]:
						for im in xrange(len(projdata)): projdata[im] = fft(projdata[im])
					mpi_barrier(MPI_COMM_WORLD) 
					do3d_sorting_groups_4nn_iter(projdata, iteration = Iter+1)
				else: ERROR("Wrong interpolation method for do3d", "mref_ali3d_Kmeans_remove_small_groups", 1, Blockdata["myid"]) 
				Iter +=1
		mpi_barrier(MPI_COMM_WORLD)
			
	res_sort3d       = get_sorting_all_params(projdata)		
	if(Blockdata["myid"] == Blockdata["main_node"]):
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print(line, "mref_ali3d_Kmeans_remove_small_groups done")
		log.add("mref_ali3d_Kmeans_remove_small_groups done")
		group_list, ali3d_params_list = parsing_sorting_params(partids, res_sort3d)
		write_text_row(group_list, os.path.join(Tracker["directory"],"list.txt"))
		clusters = extract_clusters_from_partition(group_list, Tracker["number_of_groups"])
		for icluster in xrange(len(clusters)):
			write_text_row(clusters[icluster], os.path.join(Tracker["directory"],"Class%d.txt"%icluster))
		Tracker["partition"] = group_list
	else:
		Tracker["partition"] = 0
		group_list           = 0
	Tracker["partition"] = wrap_mpi_bcast(Tracker["partition"], Blockdata["main_node"])
	group_list           = wrap_mpi_bcast(group_list, Blockdata["main_node"])
	if(Blockdata["myid"] == Blockdata["main_node"]):
		cmd = "{} {}".format("rm -rf", os.path.join(Tracker["directory"], "tempdir"))
		if os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
			junk = cmdexecute(cmd)
		if clean_volumes:
			 cmd = "{} {}".format("rm -rf", os.path.join(Tracker["directory"], "vol_*.hdf"))
			 junk = cmdexecute(cmd)
	mpi_barrier(MPI_COMM_WORLD)	
	return group_list, npergroup 	

def mref_ali3d_EQ_Kmeans(partids, partstack, clean_volumes = True):
	global Tracker, Blockdata
	from projection import prep_vol
	from logger		import Logger,BaseLogger_Files
	
	# at this stage, nxinit is already determined; the initial reference volumes are already reconstructed
	# partids      		original particle indices in a text file for all particles for this run, can be a subset of the full dataset
	# parstack     		a text file contains ali3d parameters of particles for this run
	# particle_groups   assign group IDs for particles in this run
	
	try:	termprec = Tracker["EQKmeans_stopercent"]
	except: termprec = 6.0
	log     = Logger()
	log   	= Logger(BaseLogger_Files())
	log.prefix = Tracker["directory"]+"/"
	try: 	bckgnoise = Blockdata["bckgnoise"]
	except: Blockdata["bckgnoise"] = None
	if Tracker["focus3D"]:
		if(Blockdata["myid"] == Blockdata["main_node"]): focus = get_shrink_3dmask(Tracker["nxinit"], Tracker["focus3D"])
		else:	                                         focus = model_blank(Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"])
		bcast_EMData_to_all(focus, Blockdata["myid"], Blockdata["main_node"])
		st = Util.infomask(focus, None, True)
		if( st[0] == 0.0 ):  ERROR("incorrect focus mask, after binarize all values zero", "mref_ali3d_EQ_Kmeans", 1, Blockdata["myid"])
		focus = prep_vol(focus, 1, 1)
	else:
		focus = None
	shrinkage = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	if Tracker["constants"]["interpolation"] =="4nn": 
		projdata = get_shrink_data_sorting(partids, partstack, return_real = True, preshift = True, apply_mask = True)
	elif Tracker["constants"]["interpolation"] =="trl":  #  WRONG  PAP
		if Tracker["focus3D"]: projdata = get_shrink_data_sorting(partids, partstack, return_real = True, preshift = True, apply_mask = True)
		else:  	               projdata = get_shrink_data_sorting(partids, partstack, return_real = False, preshift = True, apply_mask = True)
	else:  ERROR("Unknown method for 3-D reconstruction", "mref_ali3d_EQ_Kmeans", 1, Blockdata["myid"])

	#oldangles					= []
	if(Blockdata["myid"] == Blockdata["main_node"]):  #Let's define total_stack here
		total_stack = len(read_text_file(partids))
	else:  total_stack = 0
	total_stack = bcast_number_to_all(total_stack, Blockdata["main_node"], MPI_COMM_WORLD)
	assert(Tracker["total_stack"] == total_stack)	
	if(Blockdata["myid"] == Blockdata["main_node"]):
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print(line, "mref_ali3d_EQ_Kmeans")
		log.add("mref_ali3d_EQ_Kmeans")
		msg = "total_stack:  %d"%Tracker["total_stack"]
		log.add(msg)
		print(line, msg)
		msg = "number_of_groups:  %d"%Tracker["number_of_groups"]
		log.add(msg)
		print(line, msg)
		msg = "noctf:  %s"%Tracker["constants"]["noctf"]
		log.add(msg)
		print(line, msg)
		msg = "Currrent directory:  %s"%Tracker["directory"]
		log.add(msg)
		print(line, msg)
		msg = "nxinit:   %d"%Tracker["nxinit"]
		log.add(msg)
		print(line, msg)
		msg = "symmetry group:  %s"%Tracker["constants"]["symmetry"] 
		log.add(msg)
		print(line, msg)
		msg = "the total number of iterations:  %d"%Tracker["total_number_of_iterations"]
		log.add(msg)
		print(line, msg)
		msg = "the stop criterion:  %f"%termprec
		log.add(msg)
		print(line, msg)
		msg = "Iteration starts ...... "
		log.add(msg)
		print(line, msg)
	Tracker = wrap_mpi_bcast(Tracker,Blockdata["main_node"])
	nima = len(projdata)
	disps 		= []
	recvcount 	= []
	for im in xrange(Blockdata["nproc"]):
		if( im == Blockdata["main_node"]):  disps.append(0)
		else:                   disps.append(disps[im-1] + recvcount[im-1])
		ib, ie = MPI_start_end(total_stack, Blockdata["nproc"] , im)
		recvcount.append(ie - ib)
			
	# Initial reconstruction of random group partition
	#  WRONG
	if    Tracker["constants"]["interpolation"]=="trl": do3d_sorting_groups_trl_iter(projdata, iteration = 0)
	elif  Tracker["constants"]["interpolation"]=="4nn": do3d_sorting_groups_4nn_iter(projdata, iteration = 0)
	else: ERROR("Wrong interpolation method for do3d", "mref_ali3d_EQ_Kmeans", 1, Blockdata["myid"]) 
	
	for im in xrange(len(projdata)):
		phi, theta, psi, sxs, sys = get_params_proj(projdata[im])
		#oldangles.append([phi, theta, psi])
		projdata[im].set_attr("group", -1) # after reconstruction on initial partition, reset group id for all particles
	
	for Iter in xrange(Tracker["total_number_of_iterations"]):
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		if(Blockdata["myid"] == Blockdata["main_node"]): print(line, "Iteration    %d"%Iter)
		peaks =  [[ -1.0e23 for im in xrange(nima)] for iref in xrange(Tracker["number_of_groups"])]
		if Tracker["constants"]["low_pass_filter"] <0.0: Tracker["lowpass"] =  min(Tracker["fsc05"])/float(Tracker["nxinit"])
		else:                                            Tracker["lowpass"] =  Tracker["constants"]["low_pass_filter"]
		if Tracker["constants"]["interpolation"]=="4nn":
			if not Tracker["focus3D"]:
				for im in xrange(len(projdata)): projdata[im] = fft(projdata[im])
		if Tracker["constants"]["interpolation"]=="trl":
			if Tracker["focus3D"]: 
				for im in xrange(len(projdata)): projdata[im] = fft(projdata[im])
		mpi_barrier(MPI_COMM_WORLD)

		for iref in xrange(Tracker["number_of_groups"]):
			if(Blockdata["myid"] == Blockdata["main_node"]):
				#print("lowpass", Tracker["lowpass"])
				msg = "Iter  %d low_pass filter  %f is applied to reference %d "%(Iter, Tracker["lowpass"], iref)
				log.add(msg)
				ref_vol = get_im(os.path.join(Tracker["directory"],"vol_grp%03d_iter%03d.hdf"%(iref, Iter)))
				if Tracker["mask3D"]: 
					mask3D = get_im(Tracker["mask3D"])
					if mask3D.get_xsize()!= ref_vol.get_xsize(): mask3D = get_shrink_3dmask(ref_vol.get_xsize(), Tracker["mask3D"])
					stat     = Util.infomask(ref_vol, mask3D, False)
					ref_vol -= stat[0]
					if stat[1] !=0.0: Util.mul_scalar(ref_vol, 1.0/stat[1])
					Util.mul_img(ref_vol, mask3D)
				nnn = ref_vol.get_xsize()
				if(Tracker["nxinit"] != nnn): ref_vol = fdecimate(ref_vol, Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"], True, False)
			else:	ref_vol = model_blank(Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"])
			mpi_barrier(MPI_COMM_WORLD)
			bcast_EMData_to_all(ref_vol, Blockdata["myid"], Blockdata["main_node"])
			ref_peaks = compare_two_images(projdata, ref_vol, Tracker["lowpass"], focus, comparison_method = Tracker["constants"]["comparison_method"])
			for im in xrange(nima):	peaks[iref][im] = ref_peaks[im]
			mpi_barrier(MPI_COMM_WORLD)	
			
		if(Blockdata["myid"] == Blockdata["main_node"]): print(line, "Image and reference comparision is done  ")
		from numpy import float32, empty, inner, abs
		if( Blockdata["myid"] == Blockdata["main_node"]):
			dtot = empty( (Tracker["number_of_groups"], Tracker["total_stack"]), dtype = float32)
		for  iref in xrange(Tracker["number_of_groups"]):
			recvbuf = mpi_gatherv(peaks[iref], nima, MPI_FLOAT, recvcount, disps, MPI_FLOAT, Blockdata["main_node"], MPI_COMM_WORLD)
			if(Blockdata["myid"] == Blockdata["main_node"]): 
				dtot[iref] = recvbuf
		mpi_barrier(MPI_COMM_WORLD)	
		del recvbuf
		
		#  The while loop over even angles delta should start here.
		#  prepare reference directions
		from utilities import even_angles, getvec
		refa = even_angles(60.0)# globular proteins
		numrefang = len(refa)
		refanorm = empty( (numrefang, 3), dtype = float32)
		for i in xrange(numrefang):
			tmp = getvec(refa[i][0], refa[i][1])
			for j in xrange(3):	refanorm[i][j] = tmp[j]
		del  refa, tmp
		transv = empty((nima, 3), dtype = float32)
		
		for im in xrange(nima):
			trns = projdata[im].get_attr("xform.projection")
			for j in xrange(3): transv[im][j] = trns.at(2,j)
		#  We have all vectors, now create a list of assignments of images to references
		refassign = [-1]*nima
		for im in xrange(nima):
			refassign[im] = abs(inner(refanorm,transv[im])).argmax()
		assigntorefa = mpi_gatherv(refassign, nima, MPI_INT, recvcount, disps, MPI_INT, Blockdata["main_node"], MPI_COMM_WORLD)
		assigntorefa = map(int, assigntorefa)
		del refassign, refanorm, transv
		####
		if Blockdata["myid"] == Blockdata["main_node"]:
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			print(line, "move on one")
			SA  = False
			asi = [[] for iref in xrange(Tracker["number_of_groups"])]
			report_error = 0
			for imrefa in xrange(numrefang):
				from utilities import findall
				N = findall(imrefa, assigntorefa)
				current_nima = len(N)
				if( current_nima >= Tracker["number_of_groups"] and report_error == 0):
					tasi = [[] for iref in xrange(Tracker["number_of_groups"])]
					maxasi = current_nima//Tracker["number_of_groups"]
					nt = current_nima
					kt = Tracker["number_of_groups"]
					K = range(Tracker["number_of_groups"])

					d = empty( (Tracker["number_of_groups"], current_nima), dtype = float32)
					for ima in xrange(current_nima):
						for iref in xrange(Tracker["number_of_groups"]):  d[iref][ima] = dtot[iref][N[ima]]

					while nt > 0 and kt > 0:
						l = d.argmax()
						group = l//current_nima
						ima   = l-current_nima*group
						if SA:
							J = [0.0]*Tracker["number_of_groups"]
							sJ = 0
							Jc = [0.0]*Tracker["number_of_groups"]
							for iref in xrange(Tracker["number_of_groups"]):
								J[iref] = exp(d[iref][ima]/T)
								sJ += J[iref]
							for iref in xrange(Tracker["number_of_groups"]):
								J[iref] /= sJ
							Jc[0] = J[0]
							for iref in xrange(1, Tracker["number_of_groups"]):
								Jc[iref] = Jc[iref-1]+J[iref]
							sss = random()
							for group in xrange(Tracker["number_of_groups"]):
								if( sss <= Jc[group]): break
						tasi[group].append(N[ima])
						N[ima] = -1
						for iref in xrange(Tracker["number_of_groups"]):  d[iref][ima] = -1.e10
						nt -= 1
						masi = len(tasi[group])
						if masi == maxasi:
							for im in xrange(current_nima):  d[group][im] = -1.e10
							kt -= 1
					else:
						for ima in xrange(current_nima):
							if N[ima] > -1:
								qm = -1.e10
								for iref in xrange(Tracker["number_of_groups"]):
									qt = dtot[iref][N[ima]]
									if( qt > qm ):
										qm = qt
										group = iref
								tasi[group].append(N[ima])
					del d, N, K
					if  SA:  del J, Jc
					for iref in xrange(Tracker["number_of_groups"]):
						asi[iref] += tasi[iref]
					del tasi
				else:
					report_error = 1
			#  This should be deleted only once we know that the number of images is sufficiently large, see below.
			del dtot
		else:
			assignment = []
			report_error = 0
		report_error = bcast_number_to_all(report_error, source_node = Blockdata["main_node"])
		if report_error == 1:  ERROR('Number of images within a group too small', "mref_ali3d_MPI", 1, Blockdata["myid"])
		if Blockdata["myid"] == Blockdata["main_node"]:
			assignment = [0]*total_stack
			for iref in xrange(Tracker["number_of_groups"]):
				for im in xrange(len(asi[iref])): assignment[asi[iref][im]] = iref
			del asi
		#####
		assignment = mpi_scatterv(assignment, recvcount, disps, MPI_INT, recvcount[Blockdata["myid"]], MPI_INT, Blockdata["main_node"], MPI_COMM_WORLD)
		assignment = map(int, assignment)
		#  compute number of particles that changed assignment and how many are in which group
		nchng     =  0
		npergroup = [0]*Tracker["number_of_groups"]
		for im in xrange(nima):
			iref = projdata[im].get_attr('group')
			npergroup[assignment[im]] += 1
			if( iref != assignment[im]): nchng += 1
			projdata[im].set_attr('group', assignment[im]) # reset group ID
		nchng     = mpi_reduce(nchng, 1, MPI_INT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
		npergroup = mpi_reduce(npergroup, Tracker["number_of_groups"], MPI_INT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
		npergroup = map(int, npergroup)
		terminate = 0
		if( Blockdata["myid"] == Blockdata["main_node"] ):
			ngroup=[]
			nchng = int(nchng[0])
			precn = 100*float(nchng)/float(total_stack)
			msg = "Iteration number:   %5d"%Iter
			log.add(msg)
			msg = "Number of particles that changed assignments %7d, percentage of total: %5.1f"%(nchng, precn)
			log.add(msg)
			msg = "Group ID    number of particles   fsc05   fsc143 "
			log.add(msg)
			for iref in xrange(Tracker["number_of_groups"]):
				try: 		fsc143 = Tracker["fsc143"][iref]
				except:		fsc143 = 0.0
				try: 		fsc05  = Tracker["fsc05"][iref]
				except:		fsc05  = 0.0			
				msg = "%5d         %7d   %7d   %7d"%(iref+1, npergroup[iref], fsc05, fsc143)
				log.add(msg)
				ngroup.append(int(npergroup[iref]))
			if(precn <= termprec):  terminate = 1
		else:
			ngroup = 0
		terminate = bcast_number_to_all(terminate, Blockdata["main_node"], MPI_COMM_WORLD)
		if terminate: break
		else:
			if Tracker["constants"]["interpolation"]=="4nn": 
				if not Tracker["focus3D"]:
					for im in xrange(len(projdata)):  projdata[im] = fft(projdata[im])
			if Tracker["constants"]["interpolation"]=="trl":
				if Tracker["focus3D"]: 
					for im in xrange(len(projdata)): projdata[im] = fft(projdata[im])
			ngroup = wrap_mpi_bcast(ngroup,Blockdata["main_node"])
			if    Tracker["constants"]["interpolation"]=="trl": do3d_sorting_groups_trl_iter(projdata, iteration = Iter+1)
			elif  Tracker["constants"]["interpolation"]=="4nn": do3d_sorting_groups_4nn_iter(projdata, iteration = Iter+1)
			else: ERROR("Wrong interpolation method for do3d", "mref_ali3d_EQ_Kmeans", 1, Blockdata["myid"])
			mpi_barrier(MPI_COMM_WORLD)
			
	mpi_barrier(MPI_COMM_WORLD)
	res_sort3d = get_sorting_all_params(projdata)		
	if(Blockdata["myid"] == Blockdata["main_node"]):
		Tracker["partition"], ali3d_params_list = parsing_sorting_params(partids, res_sort3d)
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print(line, "mref_ali3d_EQ_Kmeans done")
		log.add("mref_ali3d_EQ_Kmeans done")
		write_text_row(Tracker["partition"], os.path.join(Tracker["directory"],"list.txt"))
		cmd = "{} {}".format("rm -rf", os.path.join(Tracker["directory"], "tempdir"))
		if os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
			junk = cmdexecute(cmd)
		if clean_volumes:
			cmd = "{} {}".format("rm -rf", os.path.join(Tracker["directory"], "vol_*.hdf"))
			junk = cmdexecute(cmd)
	else:   Tracker["partition"] = 0
	Tracker["partition"] = wrap_mpi_bcast(Tracker["partition"],Blockdata["main_node"])
	return Tracker["partition"]
	
### Two_way comparison
def do_two_way_comparison_over_nindepruns(log_main):#  multiple way comparison
	global Tracker, Blockdata
	from utilities  import read_text_file,write_text_file
	from statistics import k_means_match_clusters_asg_new
	from math       import sqrt
	import os
	## input Tracker["partition_list"]
	## output accounted_list, unaccounted_list
	if Blockdata["myid"] ==Blockdata["main_node"]:
		msg="------->>>Two_way comparisons analysis of %3d independent runs of equal Kmeans<<<-------"%Tracker["constants"]["indep_runs"]

	if(Tracker["constants"]["indep_runs"] < 2):
		ERROR("Intend to do two-way comparision while there is only one independent run", "do_two_way_comparison_over_nindepruns", 1, Blockdata["myid"])
		
	### Two-way comparision is carried out on all nodes
	
	ptp, Tracker["org_id"] = extract_groups_from_partitions(Tracker["partition_list"], Tracker["number_of_groups"])
	mpi_barrier(MPI_COMM_WORLD)
	nc = 0
	scores   = [0.0]*(Tracker["constants"]["indep_runs"]*(Tracker["constants"]["indep_runs"]-1)/2)
	res_list = [None]*(Tracker["constants"]["indep_runs"]*(Tracker["constants"]["indep_runs"]-1)/2)
	for iptp in xrange(len(ptp)):
		for jptp in xrange(len(ptp)):
			if iptp < jptp:
				if Blockdata["myid"] == nc%Blockdata["nproc"]: # small number of nproc while large independent runs
					newindeces, list_stable, nb_tot_objs = k_means_match_clusters_asg_new(ptp[iptp], ptp[jptp])
					tt = 0.0
					for m in xrange(len(list_stable)):
						tt +=len(list_stable[m])
						unaccounted        = Tracker["total_stack"]-tt
						ratio_unaccounted  = 100.-tt/Tracker["total_stack"]*100.
						ratio_accounted    = tt/Tracker["total_stack"]*100
					scores[nc] = tt/Tracker["total_stack"]*100.0
					new_list   = []
					for any in list_stable:
						any.tolist()
						new_list.append(any)
					#two_ways_stable_member_list[(iptp,jptp)] = new_list[:]
					res_list[nc] = new_list[:][:]
					del new_list
				nc +=1
	mpi_barrier(MPI_COMM_WORLD)
	
	for ires in xrange(len(res_list)): ## make sure it works
		if Blockdata["myid"] == ires%Blockdata["nproc"]:	tmplist = res_list[ires][:][:]
		else:							tmplist = 0
		tmplist = wrap_mpi_bcast(tmplist, ires, MPI_COMM_WORLD)
		if Blockdata["myid"] == Blockdata["main_node"]:	
			res_list[ires] = tmplist[:][:]
	res_list  = wrap_mpi_bcast(res_list, Blockdata["main_node"], MPI_COMM_WORLD)
	scores    = mpi_reduce(scores, len(res_list), MPI_INT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	scores    = map(int,scores )
	scores    = wrap_mpi_bcast(scores, Blockdata["main_node"], MPI_COMM_WORLD)
	#### Score each independent run by pairwise summation
	two_way_dict  = {}
	full_scores   = {}
	nc            = 0
	Tracker["stable_classes"] = {}
	for ipp in xrange(len(ptp)):
		for jpp in xrange(len(ptp)):
			if ipp<jpp:
				full_scores[(ipp, jpp)] = scores[nc]
				full_scores[(jpp, ipp)] = scores[nc]
				Tracker["stable_classes"][(ipp, jpp)] = res_list[nc]
				Tracker["stable_classes"][(jpp, ipp)] = res_list[nc]
				nc +=1
	summed_scores = [None]*len(ptp)
	duplicate     = {}
	for ipp in xrange(len(ptp)):
		avg_scores = 0.0
		for jpp in xrange(len(ptp)):
			if ipp != jpp:	avg_scores += full_scores[(ipp, jpp)]
		summed_scores[ipp]  = avg_scores/(len(ptp)-1.0)
		try:
			kpp = two_way_dict[summed_scores[ipp]]
			duplicate[summed_scores[ipp]]    = ipp
		except:
			two_way_dict[summed_scores[ipp]] = ipp
	#### Select two independent runs that have highest scores
	summed_scores.sort()
	rate1 = summed_scores[-1]
	rate2 = summed_scores[-2]
	if rate1 != rate2:
		tmp_run1= two_way_dict[rate1]
		tmp_run2= two_way_dict[rate2]
		run1 = min(tmp_run1,tmp_run2)
		run2 = max(tmp_run1,tmp_run2)
	else:
		tmp_run1 =  two_way_dict[rate1]
		tmp_run2 =  duplicate[rate2]
		run1 =  min(tmp_run1, tmp_run2)
		run2 =  max(tmp_run1, tmp_run2)
	try:	    Tracker["selected_stable_class"] = Tracker["stable_classes"][(run1, run2)]
	except:
		try:	Tracker["selected_stable_class"] = Tracker["stable_classes"][(run2, run1)]
		except:	ERROR("The selected two independent runs are not available", "do_two_way_comparison_over_nindepruns", 1, Blockdata["myid"])
		
	####  Save both accounted ones and unaccounted ones
	accounted_list, new_index = merge_classes_into_partition_list(Tracker["selected_stable_class"])
	#if Blockdata["myid"] == Blockdata["main_node"]:
	#	write_text_file(accounted_list, os.path.join(Tracker["directory"], "a.txt"))
	#	write_text_row(new_index, os.path.join(Tracker["directory"], "b.txt"))
	#   Compute unaccounted list:
	a =  set(range(Tracker["total_stack"]))
	b =  set(accounted_list)
	unaccounted_list = list(a.difference(b))
	Tracker["unaccounted_list"] = []
	for index in xrange(len(unaccounted_list)):Tracker["unaccounted_list"].append(Tracker["org_id"][unaccounted_list[index]]) # always check
	Tracker["unaccounted_list"] = sorted(Tracker["unaccounted_list"])
	Tracker["accounted_list"]   = []
	for index in xrange(len(accounted_list)):Tracker["accounted_list"].append([new_index[index][0], Tracker["org_id"][new_index[index][1]]]) # always check  
	if Blockdata["myid"] == Blockdata["main_node"]:
		for icluster in xrange(len(Tracker["selected_stable_class"])):
			write_text_file(Tracker["selected_stable_class"][icluster], os.path.join(Tracker["directory"], "Class%d.txt"%icluster))
		log_main.add("accounted:        %8d"%len(Tracker["accounted_list"]))
		log_main.add("unaccounted:      %8d"%len(Tracker["unaccounted_list"]))
		log_main.add("total:            %8d"%(len(Tracker["accounted_list"])+len(Tracker["unaccounted_list"])))
		write_text_file(Tracker["unaccounted_list"], os.path.join(Tracker["directory"], "Unaccounted.txt"))
		write_text_row(Tracker["accounted_list"], os.path.join(Tracker["directory"], "Accounted.txt"))
		Tracker["Accounted_on_disk"]   = os.path.join(Tracker["directory"], "Accounted.txt")
		#write_text_row(new_index, os.path.join(Tracker["directory"], "new_index.txt"))
		#write_text_file(Tracker["selected_stable_class"], os.path.join(Tracker["directory"], "secleted.txt"))
		Tracker["Unaccounted_on_disk"] = os.path.join(Tracker["directory"], "Unaccounted.txt")
		avg_two_ways        = sum(scores)
		avg_two_ways_square = 0.0
		for a in scores:	avg_two_ways_square +=a*a
		#two_ways_std        = sqrt(avg_two_ways_square/float(len(scores))-avg_two_ways**2)
		#Tracker["net_rate"] = avg_two_ways-1./Tracker["number_of_groups"]*100.
		#msg="average of two-way comparison  %5.3f"%avg_two_ways
		#Tracker["log_main"].add(msg)
		#msg="net rate of two-way comparison  %5.3f"%Tracker["net_rate"]
		#Tracker["log_main"].add(msg)
		#msg="std of two-way comparison %5.3f"%two_ways_std
		#Tracker["log_main"].add(msg)
		#msg ="Score table of two_way comparison when Kgroup =  %5d"%number_of_groups
		#Tracker["log_main"].add(msg)
		#print_upper_triangular_matrix(scores,Tracker["constants"]["indep_runs"],Tracker["log_main"])
	Tracker  = wrap_mpi_bcast(Tracker, Blockdata["main_node"], MPI_COMM_WORLD)
	return
	
def patch_to_do_k_means_match_clusters_asg_new(ptp1, ptp2):
	from statistics import k_means_match_clusters_asg_new
	# patch ad hoc elements to make equal number of classes for two partitions and thus two_way comparison becomes feasible
	patch_elements = []
	if len(ptp1) != len(ptp2):
		alist = []
		blist = []
		for a in ptp1:
			if len(a)>0: alist.append(max(a))
		for b in ptp2: 
			if len(b)>0: blist.append(max(b))
		if len(alist)>0 and len(blist)>0:
			max_number = max(max(alist), max(blist))
		else:  exit() # This would never happen
		if len(ptp1) > len(ptp2):
			ndiff = len(ptp1) - len(ptp2)
			for indiff in xrange(ndiff):
				l = []
				l.append(max_number+indiff+1)
				patch_elements.append(max_number+indiff+1)
				l = array(l,"int32")
				ptp2.append(l)
		else:
			ndiff = len(ptp2)-len(ptp1)
			for indiff in xrange(ndiff):
				l = []
				l.append(max_number+indiff+1)
				patch_elements.append(max_number + indiff + 1)
				l = array(l,"int32")
				ptp1.append(l)
	newindeces, list_stable, nb_tot_objs = k_means_match_clusters_asg_new(ptp1, ptp2)
	return newindeces, list_stable, nb_tot_objs, patch_elements
	
def do_two_way_comparison_single(ptp1, ptp2, total_stack):
	newindeces, list_stable, nb_tot_objs, patch_elements = patch_to_do_k_means_match_clusters_asg_new(ptp1, ptp2)
	tt = 0.0
	for m in xrange(len(list_stable)): tt +=len(list_stable[m])
	ratio_unaccounted  = 100.-tt/total_stack*100.
	ratio_accounted    = tt/total_stack*100.
	new_list           = []
	for any in list_stable:
		any.tolist()
		new_list.append(any)
	accounted_list, new_index = merge_classes_into_partition_list(new_list)
	a = set(range(total_stack))
	b = set(accounted_list)
	unaccounted_list = list(a.difference(b))
	unaccounted_list = sorted(unaccounted_list)
	return  accounted_list, unaccounted_list, new_index	
	
#### sorting utilities

def get_shrink_data_sorting(partids, partstack, return_real = False, preshift = True, apply_mask = True, npad = 1):
	# The function will read from stack a subset of images specified in partids
	#   and assign to them parameters from partstack with optional CTF application and shifting of the data.
	# So, the lengths of partids and partstack are the same.
	#  The read data is properly distributed among MPI threads.
	# 10142015 --- preshift is set to True when doing 3-D sorting.
	# chunk_id are set when data is read in
	global Tracker, Blockdata
	from fundamentals	import resample, fshift
	from filter			import filt_ctf
	from applications	import MPI_start_end
	from EMAN2          import Region
	
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if( Blockdata["myid"] == Blockdata["main_node"]):
		print(line,"get_shrink_data_sorting")
	'''
	if( myid == main_node ):
		print "  "
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print  line, "Reading data  onx: %3d, nx: %3d, noctf: %s, applyctf: %s, preshift: %s."%(Tracker["constants"]["nnxo"], nxinit, Tracker["constants"]["CTF"], Tracker["applyctf"], preshift)
		print  "                       stack:      %s\n                       partids:     %s\n                       partstack: %s\n"%(Tracker["constants"]["stack"], partids, partstack)
	'''
	mask2D		= model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	shrinkage 	= Tracker["nxinit"]/float(Tracker["constants"]["nnxo"])
	
	#print("shrinkage is !!!!!  %f   %d"%(shrinkage,  Tracker["nxinit"]))
	radius 		= int(Tracker["constants"]["radius"] * shrinkage +0.5)
	
	if( Blockdata["myid"] == Blockdata["main_node"]):
		lpartids = read_text_file(partids, -1)
		if len(lpartids) == 1:
			lpartids = lpartids[0]
			groupids = len(lpartids)*[-1]
		else:
			groupids = lpartids[0]
			lpartids = lpartids[1]
	else:  	
		lpartids   = 0
		groupids   = 0
	lpartids   = wrap_mpi_bcast(lpartids, Blockdata["main_node"])
	groupids   = wrap_mpi_bcast(groupids, Blockdata["main_node"])
	Tracker["total_stack"]  = len(lpartids)
	#Tracker["total_stack"] = len(lpartids)
	if(Blockdata["myid"] == Blockdata["main_node"]):  partstack = read_text_row(partstack)
	else:  partstack = 0
	partstack = wrap_mpi_bcast(partstack, Blockdata["main_node"])
	
	if(Tracker["total_stack"] < Blockdata["nproc"]):
		if(Blockdata["myid"] < Tracker["total_stack"]):
			image_start = Blockdata["myid"]
			image_end   = Blockdata["myid"]+1
		else:
			image_start = 0
			image_end   = 1
	else:
		image_start, image_end = MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"])
	lpartids  = lpartids[image_start:image_end]
	groupids  = groupids[image_start:image_end]
	#  Preprocess the data
	# mask2D    =	model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	nima      =	image_end - image_start
	oldshifts = [[0.0,0.0]]#*nima
	data      = [None]*nima

	#  Note these are in Fortran notation for polar searches
	#txm = float(nxinit-(nxinit//2+1) - radius -1)
	#txl = float(2 + radius - Tracker["nxinit"]//2+1)
	txm = float(Tracker["nxinit"]-(Tracker["nxinit"]/2+1) - radius)
	txl = float(radius - Tracker["nxinit"]//2+1)
	for im in xrange(nima):
		data[im] = get_im(Tracker["constants"]["orgstack"], lpartids[im])
		if im ==0:
			if data[im].get_xsize() > Tracker["constants"]["nnxo"]:
				window_particle = True
			else:
				window_particle =False
		phi, theta, psi, sx, sy, chunk_id  = partstack[lpartids[im]][0], partstack[lpartids[im]][1], partstack[lpartids[im]][2], partstack[lpartids[im]][3], partstack[lpartids[im]][4], partstack[lpartids[im]][5]
		"""
		if( Tracker["constants"]["CTF"] and Tracker["applyctf"] ):
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)
		"""
		if preshift:# always true
			sx = int(round(sx))
			sy = int(round(sy))
			data[im]  = cyclic_shift(data[im],sx,sy)
			#set_params_proj(data[im],[phi,theta,psi,0.0,0.0])
			sx = 0.0
			sy = 0.0
		"""
		if Tracker["constants"]["wn"] !=0:
			mx = data[im].get_xsize()//2-Tracker["constants"]["nnxo"]//2
			my = data[im].get_ysize()//2-Tracker["constants"]["nnxo"]//2
			data[im] = data[im].get_clip(Region(mx,my,Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"]))
			data[im].set_attr('ctf_applied', 0)
			set_params_proj(data[im],[phi,theta,psi,0.0,0.0])
		"""		
		st = Util.infomask(data[im], mask2D, False)
		data[im] -= st[0]
		data[im] /= st[1]
		
		if apply_mask:  data[im] = cosinemask(data[im],radius = Tracker["constants"]["radius"])
		# FT
		data[im] = fft(data[im])
		
		if not Tracker["constants"]["noctf"] :
			ctf_params = data[im].get_attr("ctf")
			data[im] = fdecimate(data[im], Tracker["nxinit"]*npad, Tracker["nxinit"]*npad, 1, False, False)
			ctf_params.apix = ctf_params.apix/shrinkage
			data[im].set_attr('ctf', ctf_params)
			#if Tracker["applyctf"] :  #  This should be always False
			#	data[im] = filt_ctf(data[im], ctf_params, dopad=False)
			#	data[im].set_attr('ctf_applied', 1)
			#else:
			data[im].set_attr('ctf_applied', 0)
			if return_real :  data[im] = fft(data[im])
		else:
			ctf_params = data[im].get_attr_default("ctf", False)
			if  ctf_params:
				ctf_params.apix = ctf_params.apix/shrinkage
				data[im].set_attr('ctf', ctf_params)
				data[im].set_attr('ctf_applied', 0)
			data[im] = fdecimate(data[im], nxinit*npad, nxinit*npad, 1, True, False)
			apix = Tracker["constants"]["pixel_size"]
			data[im].set_attr('apix', apix/shrinkage)
		if not return_real:	data[im].set_attr("padffted",1)
		data[im].set_attr("npad",npad)
		set_params_proj(data[im],[phi, theta, psi, 0.0, 0.0])
		data[im].set_attr("chunk_id",chunk_id)
		data[im].set_attr("group",groupids[im])
		#data[im].set_attr("is_complex",0)
	
		#sig = Util.rotavg_fourier( data[im])
		#Blockdata["accumulatepw"][procid][im] = sig[len(sig)//2:]+[0.0]
		#oldshifts[im] = [sx,sy]
		#  resample will properly adjusts shifts and pixel size in ctf
		#data[im] = resample(data[im], shrinkage)
		#  We have to make sure the shifts are within correct range, shrinkage or not
		#set_params_proj(data[im],[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
		#  For local SHC set anchor
		#if(nsoft == 1 and an[0] > -1):
		#  We will always set it to simplify the code
		
	return data     

def compare_two_images(data, ref_vol, lowpassfilter, focus, comparison_method = "cross"):
	from filter import filt_tophatl
	from math   import sqrt
	peaks = len(data)*[None]
	
	if comparison_method == "eucd":
	
		ref_vol = prep_vol(ref_vol, npad = 1, interpolation_method = 1)
		if focus:
			if data[0].get_attr("is_complex")==0: 
				image = fft(data[0])
			nx      = image.get_xsize()
			ny      = image.get_ysize()
		else:
			nx      = data[0].get_xsize()
			ny      = data[0].get_ysize()
		ctfs        = [ctf_img_real(ny, q.get_attr('ctf')) for q in data]
		bckgnoise   = model_blank(nx,ny)
		bckgnoise  +=1
		for im in xrange(len(data)):
			phi, theta, psi, s2x, s2y = get_params_proj(data[im], xform = "xform.projection")
			rtemp = prgl(ref_vol,[ phi, theta, psi, 0.0,0.0], 1, False)
			if lowpassfilter > 0.0: rtemp = filt_tophatl(rtemp, lowpassfilter)
			rtemp.set_attr("is_complex",0)
			if(focus):
				mask2D = binarize(prgl(focus, [phi,theta,psi,-s2x,-s2y]), 1)
				tempx  = fft(data[im]*mask2D)
				tempx.set_attr("is_complex",0)
				peaks[im] = -Util.sqed(tempx, rtemp, ctfs[im], bckgnoise)
			else:
				data[im].set_attr("is_complex",0)
				peaks[im] = -Util.sqed(data[im], rtemp, ctfs[im], bckgnoise)
				data[im].set_attr("is_complex",1)
				
	elif comparison_method == "cross":
	
		volft = prep_vol(ref_vol, 1, 1)
		#  Ref is in reciprocal space
		for im in xrange(len(data)):
			phi, theta, psi, s2x, s2y = get_params_proj(data[im], xform = "xform.projection")
			ref = filt_ctf( prgl( volft, [phi, theta, psi, -s2x, -s2y], 1, False), data[im].get_attr("ctf") )
			if lowpassfilter>0.0: ref = filt_tophatl(ref, lowpassfilter)
			ref.set_attr("is_complex",0)
			ref.set_value_at(0,0,0.0)
			nrmref = sqrt(Util.innerproduct(ref, ref, None))
			if(focus):
				mask2D = binarize( prgl(focus,[phi,theta,psi,-s2x,-s2y]), 1)
				tempx  = fft(data[im]*mask2D)
				tempx.set_attr("is_complex",0)
				peak = Util.innerproduct(ref, tempx, None)
			else:
				data[im].set_attr("is_complex",0)
				peak      = Util.innerproduct(ref, data[im], None)
				data[im].set_attr("is_complex",1)
			peaks[im]     = peak/nrmref
			
	else:  ERROR("Unknown comparison method", "compare_two_images", 1)
	return peaks
	
def split_partition_into_clusters(sorting_res):
	# split groupids from indexes of particles
	id_list        = []
	clusters       = []
	final_class_id = 0
	ptp            = []
	for igen in xrange(len(sorting_res)):
		cluster_id =[]
		for im in xrange(len(sorting_res[igen])):
			#id_list.append(sorting_res[igen][im][1])
			if  sorting_res[igen][im][0] not in cluster_id:
				cluster_id.append(sorting_res[igen][im][0])
		for b in cluster_id:
			one_cluster = []
			for a in sorting_res[igen]:
				if a[0]==b:	one_cluster.append(a[1])
			clusters.append(one_cluster)
			final_class_id +=1
	return clusters

def split_partition_into_ordered_clusters(partition):
	# split groupids from indexes of particles
	# reindex groups
	clusters       = []
	cluster_id     = []
	
	for im in xrange(len(partition)):
		if  partition[im][0] not in cluster_id:
			cluster_id.append(partition[im][0])
	####
	cluster_dict = {}
	for icluster in xrange(len(cluster_id)):
		one_cluster = []
		for a in partition:
			if a[0]==icluster: 
				one_cluster.append(a[1])
				cluster_dict[a[1]]=icluster
		clusters.append(one_cluster)
		
	# create a partition list:
	new_partition = [] 
	for iptl in xrange(len(partition)):
		new_partition.append([cluster_dict[partition[iptl][1]], partition[iptl][1]])
	return clusters, new_partition
	 
def prep_ptp_single(all_lists, full_list):
	# full_list contains the initial input indexes
	# the assignment is aligned to full_list
	# convert classes into a single list ptp denoted by group id
	ad_hoc_group_ID        = len(all_lists) + 1
	ad_hoc_particle_exists = False
	a = set([])
	for b in all_lists:
		a.union(b)
	c = set(full_list)
	if list(a.difference(c)) !=[]: ERROR("Accounted and unaccounted in total do not match the total number of particles", "prep_ptp_single", 1, Blockdata["myid"])
	else:
		pdict = {}
		for iclass in xrange(len(all_lists)):
			for iptl in xrange(len(all_lists[iclass])):
				pdict[all_lists[iclass][iptl]] = iclass
		assignment = []
		for im in xrange(len(full_list)):
			#pdict[full_list[im]]
			try: group_ID =  pdict[full_list[im]]
			except:
				group_ID = ad_hoc_group_ID
				ad_hoc_particle_exists = True
			assignment.append(group_ID)
		if ad_hoc_particle_exists: ptp = convertasi(assignment, ad_hoc_group_ID)
		else:                      ptp = convertasi(assignment, len(all_lists)+1)
	return ptp

def merge_original_id_lists(original_id_lists):
	# merge particles lists with original ID into one list while stamped by new group id 
	clusters     = []
	all_id_list  = []
	for index_of_list in xrange(len(original_id_lists)):
		cluster_dict ={}
		for index_of_particle in xrange(len(original_id_lists[index_of_list])):
			try:    cluster_dict [original_id_lists[index_of_list][index_of_particle][0]].append(original_id_lists[index_of_list][index_of_particle][1])
			except: cluster_dict [original_id_lists[index_of_list][index_of_particle][0]]= [original_id_lists[index_of_list][index_of_particle][1]]
			all_id_list.append(original_id_lists[index_of_list][index_of_particle][1])
		for a in cluster_dict: clusters.append(cluster_dict[a])
		del cluster_dict
	all_id_list = sorted(all_id_list)
	cluster_dict = {}
	for index_of_cluster in xrange(len(clusters)):
		for index_of_particle in xrange(len(clusters[index_of_cluster])):
			cluster_dict[clusters[index_of_cluster][index_of_particle]] = index_of_cluster
	final_list = []
	for index_of_particle in xrange(len(all_id_list)):
		final_list.append([cluster_dict[all_id_list[index_of_particle]], all_id_list[index_of_particle]])
	return final_list, len(clusters)

def merge_classes_into_partition_list(classes_list):
	group_dict = {}
	data_list  = []
	new_index  = []
	# write_text_file(classes_list, "classes_list.txt")
	for index_of_class in xrange(len(classes_list)):
		for index_of_particle in xrange(len(classes_list[index_of_class])):
			data_list.append(classes_list[index_of_class][index_of_particle])
			group_dict[classes_list[index_of_class][index_of_particle]] = index_of_class
	data_list = sorted(data_list)
	for index_of_particle in xrange(len(data_list)):
		new_index.append([group_dict[data_list[index_of_particle]], data_list[index_of_particle]])
	return data_list, new_index
		
def get_sorting_all_params(data):
	global Tracker, Blockdata
	from utilities    import wrap_mpi_bcast
	from applications import MPI_start_end
	if Blockdata["myid"] == Blockdata["main_node"]:	total_attr_value_list = [[]]*Tracker["total_stack"]
	else:											total_attr_value_list = 0
	for myproc in xrange(Blockdata["nproc"]):
		attr_value_list = 0
		if Blockdata["myid"] == myproc:	attr_value_list = get_sorting_attr_stack(data)
		attr_value_list = wrap_mpi_bcast(attr_value_list, myproc)
		if Blockdata["myid"] == Blockdata["main_node"]:
			image_start,image_end = MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], myproc)
			total_attr_value_list = fill_in_mpi_list(total_attr_value_list, attr_value_list, image_start,image_end)
		mpi_barrier(MPI_COMM_WORLD)
	total_attr_value_list = wrap_mpi_bcast(total_attr_value_list, Blockdata["main_node"])
	return total_attr_value_list
	
def get_sorting_attr_stack(data_in_core):
	# get partitioned group ID and xform.projection parameters
	from utilities import get_params_proj
	attr_value_list = []
	for idat in xrange(len(data_in_core)):
		attr_value_list.append([data_in_core[idat].get_attr("group"), get_params_proj(data_in_core[idat],xform = "xform.projection")])
	return attr_value_list
	
def fill_in_mpi_list(mpi_list, data_list, index_start, index_end):
	for index in xrange(index_start, index_end): mpi_list[index] = data_list[index - index_start]
	return mpi_list
	
def parsing_sorting_params(partid, sorting_params_list):
	from utilities import read_text_file
	group_list        = []
	ali3d_params_list = []
	partid_list       = read_text_file(partid, -1)
	if len(partid_list)==1:
		for ielement in xrange(len(sorting_params_list)):
			group_list.append([sorting_params_list[ielement][0], partid_list[0][ielement]])
			ali3d_params_list.append(sorting_params_list[ielement][1:])
	elif len(partid_list)==2:
		for ielement in xrange(len(sorting_params_list)):
			group_list.append([sorting_params_list[ielement][0], partid_list[1][ielement]])
			ali3d_params_list.append(sorting_params_list[ielement][1:])
	else:
		ERROR("Wrong columns", "parsing_sorting_params", 1, 0)
	return group_list, ali3d_params_list
	
def reassign_groups(partid, sorting_params_list, grpslist_to_be_assigned, number_of_groups):
	from utilities import read_text_file
	# group_list contains two columns:   group_id     original_id   
	group_list, ali3d_params_list = parsing_sorting_params(partid, sorting_params_list)
	new_number_of_groups          = number_of_groups - len(grpslist_to_be_assigned)
	if new_number_of_groups<2: ERROR("Wrong number_of_groups", "reassign_groups", 1, 0)
	else:
		old_to_new_group_dict = {}
		for small_group_id in grpslist_to_be_assigned: old_to_new_group_dict[small_group_id] = -1
		new_group_id = 0
		for igroup in xrange(number_of_groups):
			if igroup not in grpslist_to_be_assigned:
				old_to_new_group_dict[igroup] = new_group_id
				new_group_id +=1
		for ielement in xrange(len(group_list)): group_list[ielement][0] = old_to_new_group_dict[group_list[ielement][0]]
	return group_list, new_number_of_groups

def convertasi(asig, number_of_groups):
	from numpy import array
	p = []
	for k in xrange(number_of_groups):
		l = []
		for i in xrange(len(asig)):
			if( asig[i]== k ): l.append(i)
		l = array(l,"int32")
		l.sort()
		p.append(l)
	return p

def extract_groups_from_partitions(partition_list, number_of_groups):
	# Given multiple partitions in partition_list
	ptp=[None]*len(partition_list)
	for ipt in xrange(len(partition_list)):
		assignment  =[-1]*len(partition_list[ipt])
		for index_of_particle in xrange(len(partition_list[ipt])):
			assignment[index_of_particle] = partition_list[ipt][index_of_particle][0]
		ptp[ipt] = convertasi(assignment, number_of_groups)
	org_id = []
	for a in partition_list[0]:# extract org id from the first partition
		org_id.append(a[1])
	org_id = sorted(org_id)
	return ptp, org_id
	
def get_res(res_curve):
	fsc05  = 0
	fsc143 = 0 
	for ifreq in xrange(1, len(res_curve)):	
		if res_curve[ifreq] <0.5:	break
	res_05  = ifreq - 1
	for ifreq in xrange(1, len(res_curve)):
		if res_curve[ifreq]<0.143:	break
	res_143  = ifreq - 1
	return res_05, res_143
	
def get_low_pass_filter_for_Kmeans(iref):
	global Tracker, Blockdata
	if   Tracker["constants"]["Kmeans_lpf"] =="min":   Tracker["lowpass"] = float(min(Tracker["fsc05"]))/float(Tracker["nxinit"])
	elif Tracker["constants"]["Kmeans_lpf"] =="max":   Tracker["lowpass"] = float(max(Tracker["fsc05"]))/float(Tracker["nxinit"])
	elif Tracker["constants"]["Kmeans_lpf"] =="adhoc":
		if Tracker["constants"]["low_pass_filter"] < 0.0 or Tracker["constants"]["low_pass_filter"] > 0.5: 
			 ERROR("User provided low_pass filter is not set properly", "get_low_pass_filter_for_Kmeans", 1, Blockdata["myid"])
		else:   Tracker["lowpass"] = Tracker["constants"]["low_pass_filter"]
	elif Tracker["constants"]["Kmeans_lpf"] =="adaptive": Tracker["lowpass"] = float(Tracker["fsc05"][iref])/float(Tracker["nxinit"])
	elif Tracker["constants"]["Kmeans_lpf"] =="avg":      Tracker["lowpass"] = float(sum(Tracker["fsc05"]))/float(Tracker["nxinit"])/float(len(Tracker["fsc05"]))
	return
	
def extract_clusters_from_partition(partition_to_be_saved, number_of_cluster):
	clusters = []
	for i in xrange(number_of_cluster):
		clusters.append([])
	for ipar in xrange(len(partition_to_be_saved)):
		[cluster_ID, original_ID] = partition_to_be_saved[ipar]
		clusters[cluster_ID].append(original_ID)
	for icluster in xrange(len(clusters)): clusters[icluster] = sorted(clusters[icluster])
	return clusters
		
## rec3d for sorting	
def steptwo_mpi(tvol, tweight, treg, cfsc = None, regularized = True, color = 0):
	global Tracker, Blockdata
	if( Blockdata["color"] != color ):  return model_blank(1)  #  This should not be executed if called properly
	if( Blockdata["myid_on_node"] == 0 ):
		nz = tweight.get_zsize()
		ny = tweight.get_ysize()
		nx = tweight.get_xsize()
		tvol.set_attr("is_complex",1)
		if regularized:
			nr = len(cfsc)
			limitres = 0
			for i in xrange(nr):
				cfsc[i] = min(max(cfsc[i], 0.0), 0.999)
				#print( i,cfsc[i] )
				if( cfsc[i] == 0.0 ):
					limitres = i-1
					break
			if( limitres == 0 ): limitres = nr-2;
			ovol = reshape_1d(cfsc, nr, 2*nr)
			limitres = 2*min(limitres, Tracker["maxfrad"])  # 2 on account of padding, which is always on
			maxr2 = limitres**2
			for i in xrange(limitres+1, len(ovol), 1):   ovol[i] = 0.0
			ovol[0] = 1.0
			#print(" ovol  ", ovol)
			it = model_blank(2*nr)
			for i in xrange(2*nr):  it[i] = ovol[i]
			del ovol
			#  Do not regularize first four
			for i in xrange(5):  treg[i] = 0.0
			Util.reg_weights(tweight, treg, it)
			del it
		else:
			limitres = 2*min(Tracker["constants"]["nnxo"]//2, Tracker["maxfrad"])
			maxr2 = limitres**2
		#  Iterative weights
		if( Tracker["constants"]["symmetry"] != "c1" ):
			tvol    = tvol.symfvol(Tracker["constants"]["symmetry"], limitres)
			tweight = tweight.symfvol(Tracker["constants"]["symmetry"], limitres)

	else:
		tvol = model_blank(1)
		tweight = model_blank(1)
		nz    = 0
		ny    = 0
		nx    = 0
		maxr2 = 0

	nx    = bcast_number_to_all(nx, source_node = 0, mpi_comm = Blockdata["shared_comm"])
	ny    = bcast_number_to_all(ny, source_node = 0, mpi_comm = Blockdata["shared_comm"])
	nz    = bcast_number_to_all(nz, source_node = 0, mpi_comm = Blockdata["shared_comm"])
	maxr2 = bcast_number_to_all(maxr2, source_node = 0, mpi_comm = Blockdata["shared_comm"])

	vol_data = get_image_data(tvol)
	we_data = get_image_data(tweight)
	#  tvol is overwritten, meaning it is also an output
	ifi = mpi_iterefa( vol_data.__array_interface__['data'][0] ,  we_data.__array_interface__['data'][0] , nx, ny, nz, maxr2, \
			Tracker["constants"]["nnxo"], Blockdata["myid_on_node"], color, Blockdata["no_of_processes_per_group"],  Blockdata["shared_comm"])
	#Util.iterefa(tvol, tweight, maxr2, Tracker["constants"]["nnxo"])
	
	if( Blockdata["myid_on_node"] == 0 ):
		#  Either pad or window in F space to 2*nnxo
		nx = tvol.get_ysize()
		if( nx > 2*Tracker["constants"]["nnxo"]):
			tvol = fdecimate(tvol, 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], False, False)
		elif(nx < 2*Tracker["constants"]["nnxo"]):
			tvol = fpol(tvol, 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], RetReal = False, normalize = False)

		tvol = fft(tvol)
		tvol = cyclic_shift(tvol,Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
		tvol = Util.window(tvol, Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
		#tvol = cosinemask(tvol, Tracker["constants"]["nnxo"]//2-1,5, None)
		tvol.div_sinc(1)
		tvol = cosinemask(tvol, Tracker["constants"]["nnxo"]//2-1,5, None)
		return tvol
	else:  return None

def recons3d_4nnsorting_MPI(myid, main_node, prjlist, random_subset, CTF = True, upweighted = True, mpi_comm= None, target_size=-1):
	"""
		recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
		Input
			list_of_prjlist: list of lists of projections to be included in the reconstruction
	"""
	from utilities		import reduce_EMData_to_root, random_string, get_im, findall, model_blank
	from filter			import filt_table
	from reconstruction import insert_slices_pdf
	from fundamentals	import fft
	from statistics	    import fsc
	from EMAN2			import Reconstructors
	from mpi			import MPI_COMM_WORLD, mpi_barrier
	import types
	import datetime
	#if mpi_comm == None: 
	mpi_comm = MPI_COMM_WORLD

	imgsize = prjlist[0].get_ysize()  # It can be Fourier, so take y-size

	refvol = model_blank(target_size)
	refvol.set_attr("fudge", 1.0)


	if CTF: do_ctf = 1
	else:   do_ctf = 0

	fftvol = EMData()
	weight = EMData()

	from utilities import info
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = Reconstructors.get( "nn4_ctfw", params )
	r.setup()
	
	#if norm_per_particle == None: norm_per_particle = len(prjlist)*[1.0]
	
	for im in xrange(len(prjlist)):
		phi, theta, psi, s2x, s2y = get_params_proj(prjlist[im], xform = "xform.projection") # shifts are already applied 
		if random_subset == 2:
			try:	bckgn = prjlist[im].get_attr("bckgnoise")
			except:	bckgn = target_size*[1.]
			if prjlist[im].get_attr("is_complex")==0:	prjlist[im] = fft(prjlist[im]) 
			prjlist[im].set_attr_dict({"padffted":1, "is_complex":1})
			if not upweighted:  prjlist[im] = filt_table(prjlist[im], bckgn)
			prjlist[im].set_attr("bckgnoise", bckgn)
			r.insert_slice(prjlist[im], Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi}), 1.0)
		else:
			if prjlist[im].get_attr("chunk_id") == random_subset:
				try:	bckgn = prjlist[im].get_attr("bckgnoise")
				except:	bckgn = target_size*[1.]
				if prjlist[im].get_attr("is_complex")==0:	
					prjlist[im] = fft(prjlist[im]) 
					#prjlist[im].set_attr("is_complex",1)
				prjlist[im].set_attr_dict({"padffted":1, "is_complex":1})
				if not upweighted:  prjlist[im] = filt_table(prjlist[im], bckgn)
				prjlist[im].set_attr("bckgnoise", bckgn)
				r.insert_slice(prjlist[im], Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi}), 1.0)
						
	#  clean stuff
	reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
	reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)

	if myid == main_node: dummy = r.finish(True)
	mpi_barrier(mpi_comm)

	if myid == main_node: return fftvol, weight, refvol
	else: return None, None, None

def recons3d_4nnsorting_group_MPI(myid, main_node, prjlist, random_subset, group_ID, CTF = True, upweighted = True, mpi_comm= None, target_size=-1):
	"""
		recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
		Input
			list_of_prjlist: list of lists of projections to be included in the reconstruction
	"""
	from utilities      import reduce_EMData_to_root, random_string, get_im, findall
	from EMAN2          import Reconstructors
	from utilities      import model_blank, info
	from filter		    import filt_table
	from mpi            import MPI_COMM_WORLD, mpi_barrier
	from statistics     import fsc 
	from reconstruction import insert_slices_pdf
	from fundamentals   import fft
	import datetime, types
	if mpi_comm == None: mpi_comm = MPI_COMM_WORLD
	
	imgsize = prjlist[0].get_ysize()  # It can be Fourier, so take y-size

	refvol = model_blank(target_size)
	refvol.set_attr("fudge", 1.0)

	if CTF: do_ctf = 1
	else:   do_ctf = 0

	fftvol = EMData()
	weight = EMData()

	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = Reconstructors.get( "nn4_ctfw", params )
	r.setup()
	#if norm_per_particle == None: norm_per_particle = len(prjlist)*[1.0]
	for im in xrange(len(prjlist)):
		phi, theta, psi, s2x, s2y = get_params_proj(prjlist[im], xform = "xform.projection") # shifts are already applied
		if prjlist[im].get_attr("group") == group_ID:
			#print(im, Blockdata["myid"], group_ID)
			if random_subset == 2:
				try:	bckgn = prjlist[im].get_attr("bckgnoise")
				except:	bckgn = target_size*[1.]
				if prjlist[im].get_attr("is_complex") == 0:	prjlist[im] = fft(prjlist[im]) 
				prjlist[im].set_attr_dict({"padffted":1, "is_complex":1})
				if not upweighted:  prjlist[im] = filt_table(prjlist[im], bckgn)
				prjlist[im].set_attr("bckgnoise", bckgn)
				r.insert_slice(prjlist[im], Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi}), 1.0)
			else:
				if	prjlist[im].get_attr("chunk_id") == random_subset:
					try:	bckgn = prjlist[im].get_attr("bckgnoise")
					except:	bckgn = target_size*[1.]
					if prjlist[im].get_attr("is_complex")==0:	
						prjlist[im] = fft(prjlist[im]) 
						#prjlist[im].set_attr("is_complex",1)
					prjlist[im].set_attr_dict({"padffted":1, "is_complex":1})
					if not upweighted:  prjlist[im] = filt_table(prjlist[im], bckgn)
					prjlist[im].set_attr("bckgnoise", bckgn)
					r.insert_slice(prjlist[im], Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi}), 1.0)	
	reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
	reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)

	if myid == main_node: dummy = r.finish(True)
	mpi_barrier(mpi_comm)

	if myid == main_node: return fftvol, weight, refvol
	else: return None, None, None

def do3d_sorting(procid, data):
	global Tracker, Blockdata	
	if(Blockdata["myid"] == Blockdata["main_node"]):
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print(line, "do3d_sorting")
		
	tvol, tweight, trol = recons3d_4nnsorting_MPI(myid = Blockdata["myid"], main_node = Blockdata["nodes"][procid], prjlist = data,\
					random_subset = procid, CTF = (not Tracker["constants"]["noctf"]), upweighted = False, target_size = (2*Tracker["nxinit"]+3))
	if(Blockdata["myid"] == Blockdata["nodes"][procid]):
		if( procid == 0):
			cmd = "{} {}".format("mkdir", os.path.join(Tracker["directory"],"tempdir"))
			if os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
				print("tempdir exists")
			else:
				junk = cmdexecute(cmd)
		tvol.set_attr("is_complex",0)
		tvol.write_image(os.path.join(Tracker["directory"], "tempdir", "tvol_%01d.hdf"%procid))
		tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%01d.hdf"%procid))
		trol.write_image(os.path.join(Tracker["directory"], "tempdir", "trol_%01d.hdf"%procid))
	mpi_barrier(MPI_COMM_WORLD)
	return  

def stepone(tvol, tweight):
	global Tracker, Blockdata
	tvol.set_attr("is_complex",1)
	ovol = Util.shrinkfvol(tvol,2)
	owol = Util.shrinkfvol(tweight,2)
	if( Tracker["constants"]["symmetry"] != "c1" ):
		ovol = ovol.symfvol(Tracker["constants"]["symmetry"], -1)
		owol = owol.symfvol(Tracker["constants"]["symmetry"], -1)
	return Util.divn_cbyr(ovol,owol)

def steptwo(tvol, tweight, treg, cfsc = None, regularized = True):
	global Tracker, Blockdata
	nz = tweight.get_zsize()
	ny = tweight.get_ysize()
	nx = tweight.get_xsize()
	tvol.set_attr("is_complex", 1)
	if regularized:
		nr = len(cfsc)
		limitres = 0
		for i in xrange(nr):
			cfsc[i] = min(max(cfsc[i], 0.0), 0.999)
			#print( i,cfsc[i] )
			if( cfsc[i] == 0.0 ):
				limitres = i-1
				break
		if( limitres == 0 ): limitres = nr-2;
		ovol = reshape_1d(cfsc, nr, 2*nr)
		limitres = 2*min(limitres, Tracker["maxfrad"])  # 2 on account of padding, which is always on
		maxr2 = limitres**2
		for i in xrange(limitres+1, len(ovol), 1):   ovol[i] = 0.0
		ovol[0] = 1.0
		#print(" ovol  ", ovol)
		it = model_blank(2*nr)
		for i in xrange(2*nr):  it[i] = ovol[i]
		del ovol
		#  Do not regularize first four
		for i in xrange(5):  treg[i] = 0.0
		Util.reg_weights(tweight, treg, it)
		del it
	else:
		limitres = 2*min(Tracker["constants"]["nnxo"]//2, Tracker["maxfrad"])
		maxr2 = limitres**2
	#  Iterative weights
	if( Tracker["constants"]["symmetry"] != "c1" ):
		tvol    = tvol.symfvol(Tracker["constants"]["symmetry"], limitres)
		tweight = tweight.symfvol(Tracker["constants"]["symmetry"], limitres)

	#  tvol is overwritten, meaning it is also an output
	Util.iterefa(tvol, tweight, maxr2, Tracker["constants"]["nnxo"])
	#  Either pad or window in F space to 2*nnxo
	nx = tvol.get_ysize()
	if( nx > 2*Tracker["constants"]["nnxo"] ):
		tvol = fdecimate(tvol, 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], False, False)
	elif(nx < 2*Tracker["constants"]["nnxo"] ):
		tvol = fpol(tvol, 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], RetReal = False, normalize = False)

	tvol = fft(fshift(tvol,Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"]))
	tvol = Util.window(tvol, Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	tvol = cosinemask(tvol, Tracker["constants"]["nnxo"]//2-1,5, None)
	tvol.div_sinc(1)
	return tvol
	
def do3d_sorting_groups(particle_ID_index, partstack):
	global Tracker, Blockdata
	data = get_shrink_data_sorting(particle_ID_index, partstack)
	do3d_sorting_group_bp(data)
	mpi_barrier(MPI_COMM_WORLD)
	fsc143 = 0
	fsc05  = 0
	Tracker["fsc143"]				=	0
	Tracker["fsc05"]				=	0
	res_05 						    =	Blockdata["no_of_groups"]*[0]
	res_143 					    =	Blockdata["no_of_groups"]*[0]
	for index_of_colors in xrange(Blockdata["no_of_groups"]):
	
		group_start, group_end = MPI_volume_start_end(Tracker["number_of_groups"], Blockdata["no_of_groups"], index_of_colors)
		
		if Blockdata["color"] == index_of_colors:  # It has to be 1 to avoid problem with tvol1 not closed on the disk
			
			for index_of_group in xrange(group_start, group_end):
				cfsc = 0
				if Blockdata["myid_on_node"] == 0:					
					#print("   group    %d"%index_of_group)
					#--  memory_check(Blockdata["myid"],"first node, before stepone")
					#  read volumes, shrink
					tvol0 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0_%d.hdf")%index_of_group)
					tweight0 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_0_%d.hdf")%index_of_group)
					tvol1 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_1_%d.hdf")%index_of_group)
					tweight1 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_1_%d.hdf")%index_of_group)
					Util.fuse_low_freq(tvol0, tvol1, tweight0, tweight1, 2*Tracker["constants"]["fuse_freq"])
					#print("   fuse done ")
					tag = 7007
					send_EMData(tvol1, Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					send_EMData(tweight1, Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					shrank0 	= stepone(tvol0, tweight0)
				elif Blockdata["myid_on_node"] == Blockdata["no_of_processes_per_group"]-1:
					#print("   the last proc   ")
					#--  memory_check(Blockdata["myid"],"second node, before stepone")
					#  read volumes, shrink
					
					tag = 7007
					tvol1 		= recv_EMData(0, tag, Blockdata["shared_comm"])
					tweight1 	= recv_EMData(0, tag, Blockdata["shared_comm"])
					tvol1.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1} )
					shrank1 	= stepone(tvol1, tweight1)
					#print(" receive done %d"%Blockdata["myid"])
					#  Get shrank volume, do fsc, send it to all
				mpi_barrier(Blockdata["shared_comm"])
				if 	Blockdata["myid_on_node"] == 0:
					tag = 7007					
					send_EMData(shrank0, Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					del shrank0
					lcfsc = 0
					#print(" stepone done %d"%Blockdata["myid"])
				elif Blockdata["myid_on_node"] == Blockdata["no_of_processes_per_group"]-1:
					tag = 7007
					shrank0 	= recv_EMData(0, tag, Blockdata["shared_comm"])
					#  Note shrank volumes are Fourier uncentered.
					cfsc 		= fsc(shrank0, shrank1)[1]
					write_text_row(cfsc, os.path.join(Tracker["directory"], "fsc_driver_%d.txt")%index_of_group)
					del shrank0, shrank1
					if(Tracker["nxinit"]<Tracker["constants"]["nnxo"]):
						cfsc 	= cfsc[:Tracker["nxinit"]]
						for i in xrange(len(cfsc),Tracker["constants"]["nnxo"]//2+1):  cfsc.append(0.0)
					lcfsc  = len(cfsc)							
					fsc05  = 0
					fsc143 = 0 
					for ifreq in xrange(len(cfsc)):	
						if cfsc[ifreq] < 0.5: break
					fsc05  = ifreq - 1
					for ifreq in xrange(len(cfsc)):
						if cfsc[ifreq] < 0.143: break
					fsc143 = ifreq - 1
					Tracker["fsc143"] = fsc143
					Tracker["fsc05"]  = fsc05
					line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
					print(line, "do3d_sorting_groups is done")
					#--  memory_check(Blockdata["myid"],"second node, after stepone")
				Tracker = wrap_mpi_bcast(Tracker, Blockdata["no_of_processes_per_group"]-1, Blockdata["shared_comm"])
				cfsc = wrap_mpi_bcast(cfsc, Blockdata["no_of_processes_per_group"]-1, Blockdata["shared_comm"])
				Tracker["maxfrad"] = Tracker["nxinit"]//2
				#--  memory_check(Blockdata["myid"],"first node, before steptwo")
				#  compute filtered volume
				
				if( Blockdata["myid_on_node"] == 0):
					tvol2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")%index_of_group)
					tweight2 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")%index_of_group)
					treg2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_group))
				else:
					tvol2 		= model_blank(1)
					tweight2 	= model_blank(1)
					treg2		= model_blank(1)
				tvol2 = steptwo_mpi(tvol2, tweight2, treg2, cfsc, True, color = index_of_colors)
				del tweight2, treg2
				if( Blockdata["myid_on_node"] == 0):
					res_05[index_of_group]  = Tracker["fsc05"]
					res_143[index_of_group] = Tracker["fsc143"]
				#--  memory_check(Blockdata["myid"],"first node, before masking")
					if(Tracker["mask3D"] == None):  tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
					else:  Util.mul_img(tvol2, get_im(options.mask3D))
					#--  memory_check(Blockdata["myid"],"first node, after masking")
					tvol2.write_image(os.path.join(Tracker["directory"], "vol_grp%03d_iter000.hdf"%(index_of_group)))
					#--  memory_check(Blockdata["myid"],"first node, after 1 steptwo")
					del tvol2
				mpi_barrier(Blockdata["shared_comm"])
			mpi_barrier(Blockdata["shared_comm"])	
	mpi_barrier(MPI_COMM_WORLD)
	res_05         = mpi_reduce(res_05,  Tracker["number_of_groups"], MPI_INT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	res_143        = mpi_reduce(res_143, Tracker["number_of_groups"], MPI_INT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	res_05         = map(int, res_05)
	res_143        = map(int, res_143)
	if (Blockdata["myid"] == Blockdata["main_node"]):
		Tracker["fsc143"] = res_143
		Tracker["fsc05"]  = res_05
	Tracker = wrap_mpi_bcast(Tracker,Blockdata["main_node"])
	return 
			
def do3d_sorting_group_bp(data):
	global Tracker, Blockdata
	if(Blockdata["myid"] == Blockdata["nodes"][0]):
		cmd = "{} {}".format("mkdir", os.path.join(Tracker["directory"], "tempdir"))
		if os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
			print("tempdir exists")
		else:
			junk = cmdexecute(cmd)
		
	for index_of_groups in xrange(Tracker["number_of_groups"]):
		for procid in xrange(3):
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			#if Blockdata["myid"] == Blockdata["main_node"]: 
			#print(line, "Reconstruct volume  by particle from random subset %d  %d"%(procid, index_of_groups))	
			tvol, tweight, trol = recons3d_4nnsorting_group_MPI(myid = Blockdata["myid"], main_node = Blockdata["nodes"][procid], prjlist = data,  random_subset = procid, group_ID = index_of_groups, CTF = (not Tracker["constants"]["noctf"]),\
											upweighted = False, target_size = (2*Tracker["nxinit"]+3))
			if(Blockdata["myid"] == Blockdata["nodes"][procid]):
				tvol.set_attr("is_complex",0)
				tvol.write_image(os.path.join(Tracker["directory"], "tempdir", "tvol_%d_%d.hdf"%(procid, index_of_groups)))
				tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%d_%d.hdf"%(procid, index_of_groups)))
				trol.write_image(os.path.join(Tracker["directory"], "tempdir", "trol_%01d_%d.hdf"%(procid, index_of_groups)))
			mpi_barrier(MPI_COMM_WORLD)	
	mpi_barrier(MPI_COMM_WORLD)
	return
	
def do3d_sorting_groups_trl_iter(data, iteration):
	global Tracker, Blockdata
	keepgoing = 1
	if(Blockdata["myid"] == Blockdata["nodes"][0]):
		cmd = "{} {}".format("mkdir", os.path.join(Tracker["directory"], "tempdir"))
		if os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
			print("tempdir exists")
		else:
			junk = cmdexecute(cmd)
	do3d_sorting_group_bp(data)
	mpi_barrier(MPI_COMM_WORLD)
	
	fsc143 = 0
	fsc05  = 0
	Tracker["fsc143"]				=	0
	Tracker["fsc05"]				=	0
	res_05 						    =	Tracker["number_of_groups"]*[0]
	res_143 					    =	Tracker["number_of_groups"]*[0]
	
	if Blockdata["myid"] == Blockdata["main_node"]: print("number of groups is ", Tracker["number_of_groups"])
	
	for index_of_colors in xrange(Blockdata["no_of_groups"]):
	
		group_start, group_end = MPI_volume_start_end(Tracker["number_of_groups"], Blockdata["no_of_groups"], index_of_colors)
		
		if Blockdata["color"] == index_of_colors:  # It has to be 1 to avoid problem with tvol1 not closed on the disk
			
			for index_of_group in xrange(group_start, group_end):
				cfsc = 0
				if Blockdata["myid_on_node"] == 0:
					
					#print("   group    %d"%index_of_group)
				
					#--  memory_check(Blockdata["myid"],"first node, before stepone")
					#  read volumes, shrink
					tvol0 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0_%d.hdf")%index_of_group)
					tweight0 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_0_%d.hdf")%index_of_group)
					tvol1 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_1_%d.hdf")%index_of_group)
					tweight1 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_1_%d.hdf")%index_of_group)
					Util.fuse_low_freq(tvol0, tvol1, tweight0, tweight1, 2*Tracker["constants"]["fuse_freq"])
					tag = 7007
					send_EMData(tvol1, Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					send_EMData(tweight1, Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					shrank0 	= stepone(tvol0, tweight0)
					
				elif Blockdata["myid_on_node"] == Blockdata["no_of_processes_per_group"]-1:
					#print("   the last proc   ")
					#--  memory_check(Blockdata["myid"],"second node, before stepone")
					#  read volumes, shrink
					
					tag = 7007
					tvol1 		= recv_EMData(0, tag, Blockdata["shared_comm"])
					tweight1 	= recv_EMData(0, tag, Blockdata["shared_comm"])
					tvol1.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1} )
					shrank1 	= stepone(tvol1, tweight1)
					#print(" receive done %d"%Blockdata["myid"])
					#  Get shrank volume, do fsc, send it to all
				mpi_barrier(Blockdata["shared_comm"])
				if 	Blockdata["myid_on_node"] == 0:
					tag = 7007					
					send_EMData(shrank0, Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					del shrank0
					lcfsc = 0
					#print(" stepone done %d"%Blockdata["myid"])
					
				elif Blockdata["myid_on_node"] == Blockdata["no_of_processes_per_group"]-1:
					tag = 7007
					shrank0 	= recv_EMData(0, tag, Blockdata["shared_comm"])
					#  Note shrank volumes are Fourier uncentered.
					cfsc 		= fsc(shrank0, shrank1)[1]
					write_text_row(cfsc, os.path.join(Tracker["directory"], "fsc_driver_grp%03d_iter%03d.txt")%(index_of_group,iteration))
					del shrank0, shrank1
					if(Tracker["nxinit"]<Tracker["constants"]["nnxo"]):
						cfsc 	= cfsc[:Tracker["nxinit"]]
						for i in xrange(len(cfsc),Tracker["constants"]["nnxo"]//2+1):  cfsc.append(0.0)
					lcfsc = len(cfsc)							
					fsc05  = 0
					fsc143 = 0 
					for ifreq in xrange(len(cfsc)):	
						if cfsc[ifreq] <0.5:	break
					fsc05  = ifreq - 1
					for ifreq in xrange(len(cfsc)):
						if cfsc[ifreq]<0.143:	break
					fsc143 = ifreq - 1
					Tracker["fsc143"]				=	fsc143
					Tracker["fsc05"]				=	fsc05
					line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
					print(line, "group %d  of do3d_sorting_groups_trl_iter is done"%index_of_group)
				Tracker = wrap_mpi_bcast(Tracker, Blockdata["no_of_processes_per_group"]-1, Blockdata["shared_comm"])
				cfsc    = wrap_mpi_bcast(cfsc, Blockdata["no_of_processes_per_group"]-1, Blockdata["shared_comm"])
				#fsc143 = bcast_number_to_all(fsc143, Blockdata["main_node"], MPI_COMM_WORLD)
				#fsc05 = bcast_number_to_all(fsc05, Blockdata["main_node"], MPI_COMM_WORLD)
				Tracker["maxfrad"] = Tracker["nxinit"]//2
				#--  memory_check(Blockdata["myid"],"first node, before steptwo")
				#  compute filtered volume
				
				if( Blockdata["myid_on_node"] == 0):
					tvol2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")%index_of_group)
					tweight2 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")%index_of_group)
					treg2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_group))
				else:
					tvol2 		= model_blank(1)
					tweight2 	= model_blank(1)
					treg2		= model_blank(1)
				tvol2 = steptwo_mpi(tvol2, tweight2, treg2, cfsc, True, color = index_of_colors)
				del tweight2, treg2
				if( Blockdata["myid_on_node"] == 0):
				#--  memory_check(Blockdata["myid"],"first node, before masking")
					if(Tracker["mask3D"] == None):  tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
					else:  Util.mul_img(tvol2, get_im(Tracker["constants"]["mask3D"]))
					#--  memory_check(Blockdata["myid"],"first node, after masking")
					tvol2.write_image(os.path.join(Tracker["directory"], "vol_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
					#--  memory_check(Blockdata["myid"],"first node, after 1 steptwo")
					del tvol2
					#print(Blockdata["myid"], "fsc05", Tracker["fsc05"]) 
					res_05[index_of_group]  = Tracker["fsc05"]
					res_143[index_of_group] = Tracker["fsc143"]
				mpi_barrier(Blockdata["shared_comm"])
			mpi_barrier(Blockdata["shared_comm"])
	mpi_barrier(MPI_COMM_WORLD)
	res_05         = mpi_reduce(res_05,  Tracker["number_of_groups"], MPI_INT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	res_143        = mpi_reduce(res_143, Tracker["number_of_groups"], MPI_INT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	res_05         = map(int, res_05)
	res_143        = map(int, res_143)
	if (Blockdata["myid"] == Blockdata["main_node"]):
		Tracker["fsc143"] = res_143
		Tracker["fsc05"]  = res_05
	keepgoing = bcast_number_to_all(keepgoing, source_node = Blockdata["main_node"], mpi_comm = MPI_COMM_WORLD) # always check 
	Tracker   = wrap_mpi_bcast(Tracker, Blockdata["main_node"])
	if not keepgoing:	ERROR("do3d_sorting_groups_trl_iter  %s"%os.path.join(Tracker["directory"], "tempdir"),"do3d_sorting_groups_trl_iter", 1, Blockdata["myid"]) 
	return
	
def MPI_volume_start_end(number_of_groups, ncolor, mycolor):
	igroup_start = int(round(float(number_of_groups)/ncolor*mycolor))
	igroup_end   = int(round(float(number_of_groups)/ncolor*(mycolor+1)))
	return igroup_start, igroup_end				
## conversion

def copy_refinement_tracker(tracker_refinement):
	global Tracker, Blockdata
	for key, value in Tracker:
		try:
			value_refinement = tracker_refinement[key]
			#if value != value_refinement:
			#	if Blockdata["myid"] == Blockdata["main_node"]:
			#		print(key, " in sorting set as ", value, ", while in refinement, it is set as ", value_refinement)
			if value == None and value_refinement != None:
				Tracker[key] = value_refinement
		except:
			if Blockdata["myid"] == Blockdata["main_node"]: print(key, " in sorting set as ", value, ", while in refinement, it is set as ", value_refinement)
	return
	
def print_dict(dict,theme, exclude = "refinement"):
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	print(line,theme)
	spaces = "                    "
	if exclude =="refinement": exclude = ["constants", "nodes", "yr", "shared_comm", "bckgnoise", "myid", "myid_on_node", "accumulatepw", "chunk_dict", "PW_dict", "full_list"]
	else: exclude = ["constants", "chunk_dict", "PW_dict", "full_list"]
	for key, value in sorted( dict.items() ):
		pt = True
		for ll in exclude:
			if(key == ll):
				pt = False
				break
		if pt:  print("                    => ", key+spaces[len(key):],":  ",value)

# --------------------------------------------------------------------------		
# - "Tracker" (dictionary) object
#   Keeps the current state of option settings and dataset 
#   (i.e. particle stack, reference volume, reconstructed volume, and etc)
#   Each iteration is allowed to add new fields/keys
#   if necessary. This happes especially when type of 3D Refinement or metamove changes.
#   Conceptually, each iteration will be associated to a specific Tracker state.
#   Therefore, the list of Tracker state represents the history of process.
#
#   This can be used to restart process from an arbitrary iteration.
#   
#
def get_input_from_sparx_ref3d(log_main):# case one
	# import SPARX results
	global Tracker, Blockdata
	import json
	from string import split, atoi
	import_from_sparx_refinement		= 1
	selected_iter						= 0
	Tracker_refinement					= 0
	if not os.path.exists (Tracker["constants"]["refinement_dir"]): 
		ERROR("SPARX refinement dir does not exist", "get_input_from_sparx_ref3d", 1,  Blockdata["myid"])
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] == Blockdata["main_node"]:
		msg = "Import results from SPARX 3-D refinement"
		print(line, msg)
		log_main.add(msg)
		if Tracker["constants"]["niter_for_sorting"] == -1: # take the best solution to do sorting
			msg = "Search in the directory %s ......"%Tracker["constants"]["refinement_dir"]
			print(line, msg)
			log_main.add(msg)
			niter_refinement	= 0
			while os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%niter_refinement)) and os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"],"main%03d"%niter_refinement, "Tracker_%03d.json"%niter_refinement)):
				niter_refinement +=1
			niter_refinement -=1
			if niter_refinement !=0:
				fout				= open(os.path.join(Tracker["constants"]["refinement_dir"],"main%03d"%niter_refinement, \
				 "Tracker_%03d.json"%niter_refinement),'r')
				Tracker_refinement 	= convert_json_fromunicode(json.load(fout))
				fout.close()
				try:
					fout 	= open(os.path.join(Tracker["constants"]["refinement_dir"],"main%03d"%Tracker_refinement["constants"]["best"], \
					"Tracker_%03d.json"%Tracker_refinement["constants"]["best"]),'r')
					Tracker_refinement = convert_json_fromunicode(json.load(fout))
					fout.close()
					selected_iter = Tracker_refinement["constants"]["best"]
				except:		import_from_sparx_refinement = 0
			else:			import_from_sparx_refinement = 0	
		else:
			msg = "Try to load json file ...%s"%os.path.join(Tracker["constants"]["refinement_dir"],"main%03d"%Tracker["constants"]["niter_for_sorting"],\
			 "Tracker_%03d.json"%Tracker["constants"]["niter_for_sorting"])
			print(line, msg)
			log_main.add(msg)
			try:
				fout				= open(os.path.join(Tracker["constants"]["refinement_dir"],"main%03d"%Tracker["constants"]["niter_for_sorting"], \
				"Tracker_%03d.json"%Tracker["constants"]["niter_for_sorting"]),'r')
				Tracker_refinement	= convert_json_fromunicode(json.load(fout))
				fout.close()
				selected_iter 		= Tracker["constants"]["niter_for_sorting"]
			except:	import_from_sparx_refinement = 0
			
	import_from_sparx_refinement	= bcast_number_to_all(import_from_sparx_refinement, source_node = Blockdata["main_node"])
	if import_from_sparx_refinement == 0:	
		ERROR("The best solution is not found","get_input_from_sparx_ref3d", "get_input_from_sparx_ref3d", 1, Blockdata["myid"])
		from mpi import mpi_finalize
		mpi_finalize()
		exit()			
	Tracker_refinement = wrap_mpi_bcast(Tracker_refinement, Blockdata["main_node"], communicator = MPI_COMM_WORLD)

	# Check orgstack, set correct path
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] == Blockdata["main_node"]:
		refinement_dir_path, refinement_dir_name = os.path.split(Tracker["constants"]["refinement_dir"])
		Tracker["constants"]["refinement_dir"]   = Tracker["constants"]["refinement_dir"]
		if Tracker_refinement["constants"]["stack"][0:4]=="bdb:":
			refinement_stack = "bdb:"+refinement_dir_path+Tracker_refinement["constants"]["stack"][4:]
		else:
			refinement_stack = os.path.join(refinement_dir_path, Tracker_refinement["constants"]["stack"])
			
		if not Tracker["constants"]["orgstack"]: # Use refinement stack if instack is not provided
			msg = "refinement stack", refinement_stack
			print(line, msg)
			log_main.add(msg)
			Tracker["constants"]["orgstack"] = refinement_stack #Tracker_refinement["constants"]["stack"]
			print(line, "The refinement image stack is %s"%Tracker_refinement["constants"]["stack"])
			try:	image = get_im(Tracker["constants"]["orgstack"], 0)
			except:
				print(line, "Fail to read image stack")	
				import_from_sparx_refinement = 0
		else:
			if Tracker["constants"]["orgstack"] == Tracker_refinement["constants"]["stack"]: # instack and refinement data stack is the same
				msg = "The sorting instack is the same refinement instack: %s"%Tracker_refinement["constants"]["stack"]
				print(line, msg)
				log_main.add(msg)
				if not os.path.exists(Tracker["constants"]["orgstack"]): import_from_sparx_refinement = 0
			else: # complicated cases
				if (not os.path.exists(Tracker["constants"]["orgstack"])) and (not os.path.exists(Tracker_refinement["constants"]["stack"])): 
					import_from_sparx_refinement = 0
				elif (not os.path.exists(Tracker["constants"]["orgstack"])) and os.path.exists(Tracker_refinement["constants"]["stack"]):
					Tracker["constants"]["orgstack"] = refinement_stack
					msg = "Use refinement orgstack "
					print(line, msg)
					log_main.add(msg)
				else:
					msg = "Use orgstack provided by options"
					print(line, msg)
					log_main.add(msg)
		if import_from_sparx_refinement:
			msg =  "data stack for sorting is %s"%Tracker["constants"]["orgstack"]
			print(line, msg)
			log_main.add(msg)
		total_stack   = EMUtil.get_image_count(Tracker["constants"]["orgstack"])
	else: total_stack = 0
	import_from_sparx_refinement = bcast_number_to_all(import_from_sparx_refinement, source_node = Blockdata["main_node"])
	
	if import_from_sparx_refinement == 0:
		ERROR("The data stack is not accessible","get_input_from_sparx_ref3d",1)
		from mpi import mpi_finalize
		mpi_finalize()
		exit()
	
	total_stack = bcast_number_to_all(total_stack, source_node = Blockdata["main_node"])			
	Tracker["constants"]["total_stack"] = total_stack
	
	# Now copy relevant refinement files to sorting directory:
	if Blockdata["myid"] == Blockdata["main_node"]:
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iter, "params_%03d.txt"%selected_iter)):
			cmd = "{} {} {}".format("cp ",os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iter, \
			 "params_%03d.txt"%selected_iter), os.path.join(Tracker["constants"]["masterdir"], "sparx_refinement_params.txt"))
			junk = cmdexecute(cmd)
		else: import_from_sparx_refinement = 0
		Tracker["constants"]["selected_iter"] = selected_iter
	import_from_sparx_refinement			= bcast_number_to_all(import_from_sparx_refinement, source_node = Blockdata["main_node"])
	if import_from_sparx_refinement == 0:	
		ERROR("The parameter file of the best solution is not accessible", "get_input_from_sparx_ref3d", 1)
		from mpi import mpi_finalize
		mpi_finalize()
		exit()
			
	if Blockdata["myid"] == Blockdata["main_node"]:
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iter, "driver_%03d.txt"%selected_iter)):
			cmd = "{} {} {}".format("cp ",os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iter, \
			 "driver_%03d.txt"%selected_iter), os.path.join(Tracker["constants"]["masterdir"], "fsc_global.txt"))
			junk = cmdexecute(cmd)
		else: import_from_sparx_refinement = 0
		Tracker["constants"]["selected_iter"] = selected_iter
		if import_from_sparx_refinement: fsc_curve = read_text_row(os.path.join(Tracker["constants"]["masterdir"], "fsc_global.txt"))
	
		fsc143	= 0
		fsc05	= 0
		for ifreq in xrange(len(fsc_curve)): # drive has only one column
			if fsc_curve[ifreq][0] < 0.5: break
		fsc05  = ifreq - 1
		for ifreq in xrange(len(fsc_curve)):
			if fsc_curve[ifreq][0] < 0.143: break
		fsc143 = ifreq - 1
	
		Tracker["constants"]["fsc143"]  = fsc143
		Tracker["constants"]["fsc05"]   = fsc05	
	import_from_sparx_refinement = bcast_number_to_all(import_from_sparx_refinement, source_node = Blockdata["main_node"])
	if import_from_sparx_refinement == 0:	
		ERROR("The driver of the best solution is not accessible","get_input_from_sparx_ref3d", 1)
		from mpi import mpi_finalize
		mpi_finalize()
		exit()		

	if Blockdata["myid"] == Blockdata["main_node"]:
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main000/indexes_000.txt")):
			cmd = "{} {} {}".format("cp ", os.path.join(Tracker["constants"]["refinement_dir"], "main000/indexes_000.txt"), \
			os.path.join(Tracker["constants"]["masterdir"], "indexes.txt"))
			junk = cmdexecute(cmd)
		else:	import_from_sparx_refinement = 0	
	import_from_sparx_refinement			= bcast_number_to_all(import_from_sparx_refinement, source_node = Blockdata["main_node"])
	if import_from_sparx_refinement == 0:	
		ERROR("The index file of the best solution are not accessible","get_input_from_sparx_ref3d", 1)
		from mpi import mpi_finalize
		mpi_finalize()
		exit()

	if Blockdata["myid"] == Blockdata["main_node"]:
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main000/chunk_0_000.txt")):
			cmd = "{} {} {}".format("cp ", os.path.join(Tracker["constants"]["refinement_dir"], "main000/chunk_0_000.txt"), \
			os.path.join(Tracker["constants"]["masterdir"], "chunk_0.txt"))
			junk = cmdexecute(cmd)
		else:	
			import_from_sparx_refinement == 0
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main000/chunk_1_000.txt")):
			cmd = "{} {} {}".format("cp ", os.path.join(Tracker["constants"]["refinement_dir"], "main000/chunk_1_000.txt"), \
			os.path.join(Tracker["constants"]["masterdir"], "chunk_1.txt"))
			junk = cmdexecute(cmd)
		else: 	
			import_from_sparx_refinement == 0
			
	import_from_sparx_refinement = bcast_number_to_all(import_from_sparx_refinement, source_node = Blockdata["main_node"])
	if import_from_sparx_refinement == 0:	
		ERROR("The chunk files of the best solution are not accessible","get_input_from_sparx_ref3d",1)
		from mpi import mpi_finalize
		mpi_finalize()
		exit()
			
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	# copy all relavant parameters into sorting tracker
	if Blockdata["myid"] == Blockdata["main_node"]:
		if Tracker["constants"]["radius"] == -1: Tracker["constants"]["radius"] = Tracker_refinement["constants"]["radius"]
		Tracker["constants"]["nnxo"] = Tracker_refinement["constants"]["nnxo"]
		"""
		if Tracker["constants"]["wn"] !=0: 
			if    Tracker["constants"]["nnxo"] < Tracker["constants"]["wn"]: ERROR("window image size is not correctly chosen!", \
			"get_input_from_sparx_ref3d", 1, Blockdata["myid"])
			else:                                                            Tracker["constants"]["nnxo"] = Tracker["constants"]["wn"]
		"""
		Tracker["constants"]["orgres"]     = Tracker_refinement["bestres"]
		Tracker["delta"]                   = Tracker_refinement["delta"]
		Tracker["an"]                      = 1.5*Tracker["delta"]
		Tracker["ts"]                      = Tracker_refinement["ts"]
		Tracker["xr"]                      = Tracker_refinement["xr"]
		Tracker["constants"]["pixel_size"] = Tracker_refinement["constants"]["pixel_size"]
		try:     sym =  Tracker_refinement["constants"]["sym"]
		except:  sym =  Tracker_refinement["constants"]["symmetry"]
		Tracker["constants"]["symmetry"]                 =  sym
		print(line, "Parameters importing is done!")
		if not Tracker["constants"]["mask3D"] and Tracker_refinement["constants"]["mask3D"]:
			refinement_mask3D_path, refinement_mask3D_file = os.path.split(Tracker_refinement["constants"]["mask3D"])# MRK_DEBUG
			cmd = "{} {} {}".format("cp ", os.path.join(refinement_dir_path, Tracker_refinement["constants"]["mask3D"]), \
			os.path.join(Tracker["constants"]["masterdir"], refinement_mask3D_file))
			junk = cmdexecute(cmd)
			Tracker["constants"]["mask3D"] = os.path.join(Tracker["constants"]["masterdir"], refinement_mask3D_file)
	Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], communicator = MPI_COMM_WORLD)
	import_from_sparx_refinement = bcast_number_to_all(import_from_sparx_refinement, source_node = Blockdata["main_node"])
	if not import_from_sparx_refinement: 
		ERROR("Import parameters from SPARX refinement failed", "get_input_from_sparx_ref3d", 1,  Blockdata["myid"])
		from mpi import mpi_finalize
		mpi_finalize()
		exit()
	# Setting for margin error				
	chunk_dict = {}
	chunk_list = []
	if(Blockdata["myid"] == Blockdata["main_node"]):
		chunk_one = read_text_file(os.path.join(Tracker["constants"]["masterdir"],"chunk_0.txt"))
		chunk_two = read_text_file(os.path.join(Tracker["constants"]["masterdir"],"chunk_1.txt"))
	else:
		chunk_one = 0
		chunk_two = 0
	chunk_one = wrap_mpi_bcast(chunk_one, Blockdata["main_node"])
	chunk_two = wrap_mpi_bcast(chunk_two, Blockdata["main_node"])
	for element in chunk_one: chunk_dict[element] = 0
	for element in chunk_two: chunk_dict[element] = 1
	chunk_list 				= [chunk_one, chunk_two]
	Tracker["chunk_dict"] 	= chunk_dict
	Tracker["P_chunk_0"]   	= len(chunk_one)/float(total_stack)
	Tracker["P_chunk_1"]   	= len(chunk_two)/float(total_stack)
	if(Blockdata["myid"] == Blockdata["main_node"]):
		chunk_ids = []
		partids   = read_text_file(os.path.join(Tracker["constants"]["masterdir"], "indexes.txt"),-1)
		partids   = partids[0]
		Tracker["constants"]["total_stack"] = len(partids)
		params    = read_text_file(os.path.join(Tracker["constants"]["masterdir"], "sparx_refinement_params.txt"),-1)
		for index_of_particle in xrange(len(partids)): chunk_ids.append(chunk_dict[partids[index_of_particle]])
		refinement_params = [params[0], params[1], params[2],   params[3], params[4], chunk_ids]
		write_text_file(refinement_params, os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt"))		
		line      = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print(line, "Initialization of sorting from SPARX refinement is done")
	else:
		Tracker["constants"]["total_stack"] = 0
	Tracker["constants"]["total_stack"]     = bcast_number_to_all(Tracker["constants"]["total_stack"], Blockdata["main_node"], MPI_COMM_WORLD)
	Tracker["total_stack"]                  = Tracker["constants"]["total_stack"]
	Tracker["constants"]["partstack"]	    = os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt")
	total_stack                             = Tracker["constants"]["total_stack"]
	Tracker["currentres"]                   = float(Tracker["constants"]["fsc05"])/float(Tracker["constants"]["nxinit"])
	return import_from_sparx_refinement
	
def get_input_from_relion_ref3d(log_main):# case two
	global Tracker, Blockdata
	import json
	from   string import atoi, atof
	from   os     import walk
	import_from_relion_refinement = 1
	random_subsets_set	          = 1
	orgstack_exists		          = 1
	if not os.path.exists (Tracker["constants"]["refinement_dir"]): ERROR("relion refinement dir does not exist", \
	"get_input_from_relion_ref3d", 1, Blockdata["myid"])
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] == Blockdata["main_node"]: # main_node handles all relion files
		msg = "Import results from relion 3-D refinement"
		print(line, msg)
		log_main.add(msg)
		relion_files			= []
		data_star_for_sorting	= None
		for (dirpath, dirnames, filenames) in walk(Tracker["constants"]["refinement_dir"]):
			relion_files.extend(filenames)
			break
		if Tracker["constants"]["niter_for_sorting"] == -1: # Search for the final or the best for non-finished run
			data_star_dict	= {} 		  # sort out how many (including ct run) iterations in relion dir									
			for one_file in relion_files:
				if one_file [-9:] =="data.star":
					try:	data_star_dict[len(one_file)].append(one_file)
					except:	data_star_dict[len(one_file)]=[one_file]
			#print("all data star files are", data_star_dict)
			if len(data_star_dict[min(data_star_dict.keys())])== 1: # finished run, take the final star file
				relion_run_root_name	= data_star_dict[min(data_star_dict.keys())][0][min(data_star_dict.keys())-9]
				data_star_for_sorting 	= data_star_dict[min(data_star_dict.keys())][0]
				msg = "The final parameter file is %s"%data_star_for_sorting
				print(msg)
				log_main.add(msg)
			else:	#Unfinished run, search for the last iteration
				iter_num_list 	= []
				data_star_list 	= data_star_dict[max(data_star_dict.keys())]
				for one_file in data_star_list: iter_num_list.append(atoi(one_file[-13:-10]))
				data_star_for_sorting		= data_star_list[0][0:-13]+"%03d_data.star"%max(iter_num_list)
				optimiser_star_for_sorting	= data_star_list[0][0:-13]+"%03d_optimiser.star"%max(iter_num_list)
				sampling_star_for_sorting	= data_star_list[0][0:-13]+"%03d_sampling.star"%max(iter_num_list)
				model_star_for_sorting		= data_star_list[0][0:-13]+"%03d_half1_model.star"%max(iter_num_list)
				print("data for sorting", data_star_for_sorting)
		else: # start from the given iteration
			for one_file in relion_files:
				if one_file [:-15] =="it%03d_data.star"%Tracker["constants"]["niter_for_sorting"]:
					data_star_for_sorting = one_file
					optimiser_star_for_sorting	= one_file[0:-15]+"it%03d_optimiser.star"%Tracker["constants"]["niter_for_sorting"]
					sampling_star_for_sorting	= one_file[0:-15]+"it%03d_sampling.star"%Tracker["constants"]["niter_for_sorting"]
					model_star_for_sorting		= one_file[0:-15]+"it%03d_half1_model.star"%Tracker["constants"]["niter_for_sorting"]
					break
		if data_star_for_sorting is None:	import_from_relion_refinement = 0	
		else:								print("The selected parameter file is %s"%data_star_for_sorting)
	import_from_relion_refinement = bcast_number_to_all(import_from_relion_refinement, source_node = Blockdata["main_node"])
	if import_from_relion_refinement == 0: ERROR("relion refinement params file is not found","get_input_from_relion_ref3d", 1, Blockdata["myid"])

	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] == Blockdata["main_node"]: # 
		total_stack  = 0
		if ( not Tracker["constants"]["orgstack"]) or (not os.path.exists(Tracker["constants"]["orgstack"])):
			msg =  "Create orgstack from relion refinement data"
			print(line, msg)
			Tracker["constants"]["orgstack"]  = os.path.join(Tracker["constants"]["masterdir"],"stack_from_relion.hdf")
			orgstack_exists = 0
		else:
			msg = "orgstack is provided for sorting"
			print(line, msg)
			log_main.add(msg)
			total_stack = EMUtil.get_image_count(Tracker["constants"]["orgstack"])
	else: total_stack = 0
	total_stack = bcast_number_to_all(total_stack, source_node = Blockdata["main_node"])

	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] == Blockdata["main_node"]:# get xform.projection/ctf params
		msg = "input data file is %s"%os.path.join(Tracker["constants"]["refinement_dir"],data_star_for_sorting)
		print(lines, msg)
		log_main.add(msg)
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], data_star_for_sorting)):
			data_in_core  = read_text_row(os.path.join(Tracker["constants"]["refinement_dir"],data_star_for_sorting))
			#print("full size",len(data_in_core))
			entries = {}
			data 	= []
			for a  in data_in_core:
				if len(a)<5 and len(a) >1:# skip empty lines
					if a[0][0:1] =="_": entries[a[0]] = atoi(a[1][1:]) 
				elif len(a)>5: 			data.append(a)
			if total_stack !=0 and total_stack != len(data):  import_from_relion_refinement	 = 0
			relion_header = [] 
			for a in data_in_core:
				if len(a)>5:
					break
				relion_header.append(a)
			write_text_row(relion_header, os.path.join(Tracker["constants"]["masterdir"], "relion_header.star"))
		else:
			import_from_relion_refinement	 = 0
		msg = "Read in relion data star file"
		print(line, msg)
		log_main.add(msg)
		#print(entries, len(data))
		if total_stack ==0:	total_stack = len(data)
		if total_stack ==0: import_from_relion_refinement = 0
		if total_stack !=0:	
			write_text_row(data, os.path.join(Tracker["constants"]["masterdir"], "relion_data.star"))
			write_text_file(range(len(data)), os.path.join(Tracker["constants"]["masterdir"], "indexes.txt"))
	else:
		total_stack = 0
	import_from_relion_refinement = bcast_number_to_all(import_from_relion_refinement, source_node = Blockdata["main_node"])					
	if import_from_relion_refinement == 0:  ERROR("not able to read relion refinement params file","get_input_from_relion_ref3d", 1,  Blockdata["myid"])
	total_stack = bcast_number_to_all(total_stack, source_node = Blockdata["main_node"])
	Tracker["constants"]["total_stack"] = total_stack

	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] == Blockdata["main_node"]: # data_star 
		msg = "Create HDF stack file for sorting... "
		print(line, msg)
		log_main.add(msg)
		if orgstack_exists:	
			msg = "orgstack headers are to be altered! "
			print(line, msg)
			log_main.add(msg)
		chunk_0				= []
		chunk_1				= []
		data_xform_proj		= []
		data_ctf			= []
		for iparticle in xrange(len(data)): # conversion loop
			#<<<<<---------->>>xform.projection<<<--------------------
			try:
				sx 		= data[iparticle][entries["_rlnOriginX"]-1]
				sy 		= data[iparticle][entries["_rlnOriginY"]-1]
				phi 	= data[iparticle][entries["_rlnAngleRot"]-1]
				if phi <0.0: phi +=360.
				theta 	= data[iparticle][entries["_rlnAngleTilt"]-1]
				psi 	= data[iparticle][entries["_rlnAnglePsi"]-1]
				if psi<0.0: psi  +=360.
				norm	= data[iparticle][entries["_rlnNormCorrection"]-1]
				random_subset           = data[iparticle][entries["_rlnRandomSubset"]-1]-1
				xform_proj_p	= [phi,theta,psi,sx,sy, random_subset]
				data_xform_proj.append(xform_proj_p)
			except:	
				import_from_relion_refinement = 0
				break
			#<<<------------->>> CTF params <<<---------------------
			try:
				voltage 				= data[iparticle][entries["_rlnVoltage"]-1]
				defocus_u 				= data[iparticle][entries["_rlnDefocusU"]-1]
				defocus_v 				= data[iparticle][entries["_rlnDefocusV"]-1]
				defocus_angle 			= data[iparticle][entries["_rlnDefocusAngle"]-1]
				cs 						= data[iparticle][entries["_rlnSphericalAberration"]-1]
				dps						= data[iparticle][entries["_rlnDetectorPixelSize"]-1]
				mag						= data[iparticle][entries["_rlnMagnification"]-1]
				amplitude_contrast		= data[iparticle][entries["_rlnAmplitudeContrast"]-1]
				defocus					= ( defocus_u + defocus_v)/20000.
				astigmatism_amplitude	= (-defocus_u + defocus_v)/10000.
				astigmatism_angle		= 45. - defocus_angle
				if astigmatism_angle>180.:	astigmatism_angle -= 180.
				if astigmatism_angle<0.0:	astigmatism_angle += 180.
				bfactor    = 0.0
				apix       = 10000.*dps/mag
				ampcont    = amplitude_contrast*100.
				ctf_p      = [defocus, cs, voltage, apix, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle]
			except:
				Tracker["constants"]["noctf"] = True
			if not Tracker["constants"]["noctf"]:	data_ctf.append(ctf_p)
			if import_from_relion_refinement: 
				##<<<-----------  
				micrograph_name   = data[iparticle][entries["_rlnMicrographName"]-1]
				org_ptl_name      = data[iparticle][entries["_rlnOriginalParticleName"]-1]
				image_name        = data[iparticle][entries["_rlnImageName"]-1]
				mrcs_stack_index  = atoi(image_name[0:6])-1
				mrcs_stack_file   = image_name[7:]
				try:
					random_subset = data[iparticle][entries["_rlnRandomSubset"]-1] 
					if random_subset == 1:	chunk_0.append(iparticle)
					else:					chunk_1.append(iparticle)
				except:	random_subsets_set = 0
				image = EMData()
			
				if not orgstack_exists:
					image.read_image(mrcs_stack_file, mrcs_stack_index)
					set_params_proj(image, xform_proj_p)
					if not Tracker["constants"]["noctf"]: set_ctf(image, ctf_p)
					image.set_attr("ptl_source_image", micrograph_name)
					image.write_image(Tracker["constants"]["orgstack"], iparticle)
				else:							
					image.read_image(Tracker["constants"]["orgstack"], iparticle, True)
					set_params_proj(image, xform_proj_p)
					if not Tracker["constants"]["noctf"]: set_ctf(image, ctf_p)
					image_set_attr("ptl_source_image", micrograph_name)
					image.write_image(Tracker["constants"]["orgstack"], iparticle, EMUtil.ImageType.IMAGE_HDF, True)
				
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		if import_from_relion_refinement:
			msg = "Extract xform.projection/CTF parameters from relion refinement %s"%Tracker["constants"]["refinement_dir"]
			print(line, msg)
			log_main.add(msg)	
			write_text_row(data_xform_proj, os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt"))
			if not Tracker["constants"]["noctf"]: write_text_row(data_ctf, os.path.join(Tracker["constants"]["masterdir"], "ctf_params.txt"))
		if random_subsets_set:
			msg = "Extract random subsets from relion refinement %s"%Tracker["constants"]["refinement_dir"]
			print(line, msg)
			log_main.add(msg)			
			write_text_file(chunk_0, Tracker["constants"]["chunk_0"])
			write_text_file(chunk_1, Tracker["constants"]["chunk_1"])
	import_from_relion_refinement = bcast_number_to_all(import_from_relion_refinement, source_node = Blockdata["main_node"])
	if not import_from_relion_refinement: ERROR("relion fails in creating orgstack","get_input_from_relion_ref3d", 1, Blockdata["myid"])
	Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], communicator = MPI_COMM_WORLD)

	# read other parameters, optimiser, sampling, half1_model
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] == Blockdata["main_node"]:#
		msg ="Read relion refinement state parameters"
		print(line, msg)
		log_main.add(msg)
		relion_dict	= {}
		optimiser_data = read_text_row(os.path.join(Tracker["constants"]["refinement_dir"], optimiser_star_for_sorting))
		nc = 0
		for odata in optimiser_data:
			try:
				if odata[0][0:4]=="_rln": relion_dict[odata[0]] = odata[1]
			except:	nc +=1
		sampling_data = read_text_row(os.path.join(Tracker["constants"]["refinement_dir"], sampling_star_for_sorting))
		nc = 0
		for sdata in sampling_data:
			try:
				if sdata[0][0:4]=="_rln": relion_dict[sdata[0]] = sdata[1]
			except:	nc +=1
		#print(model_star_for_sorting)
		model_data = read_text_row(os.path.join(Tracker["constants"]["refinement_dir"], model_star_for_sorting))
		for idata in xrange(len(model_data)):
			if  len(model_data[idata])>=1: # always skip empty lines
				if model_data[idata][0] == "data_model_classes": break
				try:
					if model_data[idata][0][0:4]=="_rln": relion_dict[model_data[idata][0]] = model_data[idata][1]
				except:	nc +=1
		nstart = idata
		fsc_curve = []
		for idata in xrange(nstart, len(model_data)):
			if  len(model_data[idata])>=1:
				if model_data[idata][0] == "data_model_class_1":	break
		nstart = idata
		for idata in xrange(nstart, len(model_data)):
			if  len(model_data[idata])>=1:
				if model_data[idata][0] == "data_model_groups":	break
				if len(model_data[idata])>3:	fsc_curve.append(model_data[idata])
		fsc05  = 0
		fsc143 = 0 
		for ifreq in xrange(len(fsc_curve)):
			if fsc_curve[ifreq][4] <0.5:	break
		fsc05  = ifreq - 1
		for ifreq in xrange(len(fsc_curve)):
			if fsc_curve[ifreq][4] <0.143:	break
		fsc143 = ifreq - 1
		Tracker["constants"]["fsc143"]				=	fsc143
		Tracker["constants"]["fsc05"]				=	fsc05
		Tracker["constants"]["orgres"]				=	fsc_curve[fsc05][1]
		Tracker["constants"]["pixel_size"]			=	relion_dict["_rlnPixelSize"]
		Tracker["constants"]["nnxo"]				=	relion_dict["_rlnOriginalImageSize"]
		"""
		if Tracker["constants"]["wn"] != 0:
			if Tracker["constants"]["nnxo"] < Tracker["constants"]["wn"]: ERROR("sxsort3d","window size is not correctly chosen",\
			 "get_input_from_relion_ref3d", 1, Blockdata["myid"])
			else:                                                         Tracker["constants"]["nnxo"] = Tracker["constants"]["wn"]
		""" 
		Tracker["constants"]["radius"]              = int(relion_dict["_rlnParticleDiameter"]/(2.*Tracker["constants"]["pixel_size"]))
		Tracker["constants"]["refinement_an"]       = Tracker["constants"]["refinement_delta"]*6.0
		Tracker["constants"]["symmetry"]            = relion_dict["_rlnSymmetryGroup"]
		Tracker["constants"]["symmetry"]            = Tracker["constants"]["symmetry"][0].lower() + Tracker["constants"]["symmetry"][1:] # be sure it is low case
		Tracker["delta"]                            = relion_dict["_rlnPsiStep"]
		Tracker["an"]                               = Tracker["delta"]*1.5
		Tracker["ts"]                               = relion_dict["_rlnOffsetStep"]
		Tracker["xr"]                               = Tracker["ts"] #relion_dict["_rlnOffsetRange"]
	
		if not Tracker["constants"]["mask3D"]:
			refinement_path, refinement_dir  =  os.path.split(Tracker["constants"]["refinement_dir"])
			#relion_mask3D_file = os.path.join(Tracker["constants"]["refinement_dir"], relion_dict["_rlnSolventMaskName"])
			relion_mask3D_path, relion_mask3D_file = os.path.split(relion_dict["_rlnSolventMaskName"])
			cmd = "{} {} {}".format("cp ", os.path.join(refinement_path, relion_dict["_rlnSolventMaskName"]), \
			os.path.join(Tracker["constants"]["masterdir"], relion_mask3D_file))
			junk = cmdexecute(cmd)
			Tracker["constants"]["mask3D"] = os.path.join(Tracker["constants"]["masterdir"], relion_mask3D_file)
				
	Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], communicator = MPI_COMM_WORLD)
	import_from_relion_refinement = bcast_number_to_all(import_from_relion_refinement, source_node = Blockdata["main_node"])
	if not import_from_relion_refinement: ERROR("sort3d fail in importing relion refinement", "get_input_from_relion_ref3d", 1,  Blockdata["myid"])
	# Setting for margin error
	chunk_dict = {}
	chunk_list = []
	if(Blockdata["myid"] == Blockdata["main_node"]):
		chunk_one = read_text_file(os.path.join(Tracker["constants"]["masterdir"],"chunk_0.txt"))
		chunk_two = read_text_file(os.path.join(Tracker["constants"]["masterdir"],"chunk_1.txt"))
		for element in chunk_one: chunk_dict[element] = 0
		for element in chunk_two: chunk_dict[element] = 1
		chunk_list 				= [chunk_one, chunk_two]
		Tracker["chunk_dict"] 	= chunk_dict
		Tracker["P_chunk_0"]   	= len(chunk_one)/float(total_stack)
		Tracker["P_chunk_1"]   	= len(chunk_two)/float(total_stack)
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		msg = "Initialization of sorting from relion refinement is done"
		print(line, msg)
		log_main.add(msg)
	Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], communicator = MPI_COMM_WORLD)
	return import_from_relion_refinement
	
def get_input_from_datastack(log_main):# case three
	global Tracker, Blockdata
	import json
	from   string import split, atoi
	from   random import shuffle
	import_from_data_stack 	= 1
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if(Blockdata["myid"] == Blockdata["main_node"]):
		msg =  "Sorting starts from given data stack %s with both xform.projection and ctf parameters filled in the headers"%Tracker["constants"]["orgstack"]
		print(line, msg)
		log_main.add(msg)
		image = get_im(Tracker["constants"]["orgstack"])
		Tracker["constants"]["nnxo"] = image.get_xsize()
		"""
		if Tracker["constants"]["wn"] != 0:
			if Tracker["constants"]["nnxo"] <Tracker["constants"]["wn"]: 
				ERROR("window size is not correct chosen","get_input_from_datastack",1, Blockdata["main_node"])
			else: Tracker["constants"]["nnxo"] = Tracker["constants"]["wn"]
		"""
		if( Tracker["nxinit"] > Tracker["constants"]["nnxo"]):
			ERROR("Image size less than minimum permitted $d"%Tracker["nxinit"],"get_input_from_datastack",1, Blockdata["myid"])
			nnxo = -1
		else:
			if not Tracker["constants"]["noctf"]:
				ictf = image.get_attr('ctf')
				Tracker["constants"]["pixel_size"] = ictf.apix
				fq = Tracker["constants"]["pixel_size"]/Tracker["constants"]["fuse_freq"]
			else:
				Tracker["constants"]["pixel_size"] = 1.0
				del image
	else:
		Tracker["constants"]["nnxo"]       = 0
		Tracker["constants"]["pixel_size"] = 1.0
		fq                                 = 0.0
	Tracker["constants"]["nnxo"]		= bcast_number_to_all(Tracker["constants"]["nnxo"], source_node = Blockdata["main_node"])
	if( Tracker["constants"]["nnxo"] < 0): ERROR("Image size is negative", "get_input_from_datastack", 1, Blockdata["main_node"])
	Tracker["constants"]["pixel_size"]	= bcast_number_to_all(Tracker["constants"]["pixel_size"], source_node =  Blockdata["main_node"])
	fq									= bcast_number_to_all(fq, source_node =  Blockdata["main_node"])
	Tracker["fuse_freq"]				= fq
	if(Tracker["constants"]["radius"] < 1): Tracker["constants"]["radius"]  = Tracker["constants"]["nnxo"]//2-2
	elif((2*Tracker["constants"]["radius"] +2) > Tracker["constants"]["nnxo"]): ERROR("Particle radius set too large!", \
	"get_input_from_datastack",1, Blockdata["myid"])
	if Blockdata["myid"] == Blockdata["main_node"]:	total_stack = EMUtil.get_image_count(Tracker["constants"]["orgstack"])
	else:											total_stack = 0
	total_stack = bcast_number_to_all(total_stack, Blockdata["main_node"])
	# randomly assign two subsets
	Tracker["constants"]["total_stack"]	= total_stack
	
	#print(Tracker["constants"]["total_stack"], "TOTAL STACK")
	Tracker["constants"]["chunk_0"]		= os.path.join(Tracker["constants"]["masterdir"],"chunk_0.txt")
	Tracker["constants"]["chunk_1"]		= os.path.join(Tracker["constants"]["masterdir"],"chunk_1.txt")
	Tracker["constants"]["partstack"]	= os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt")
	chunk_dict = {}
	chunk_list = []

	if Blockdata["myid"] == Blockdata["main_node"]:
		chunk_dict  = {}
		tlist		= range(total_stack)
		write_text_file(tlist, os.path.join(Tracker["constants"]["masterdir"], "indexes.txt"))
		shuffle(tlist)
		chunk_one	= tlist[0:total_stack//2]
		chunk_two	= tlist[total_stack//2:]
		chunk_one	= sorted(chunk_one)
		chunk_two	= sorted(chunk_two)
		write_text_row(chunk_one,Tracker["constants"]["chunk_0"])
		write_text_row(chunk_two,Tracker["constants"]["chunk_1"])
		for particle in chunk_one: chunk_dict[particle] = 0
		for particle in chunk_two: chunk_dict[particle] = 1
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		msg =  "Save randomly assigned two subsets"
		print(line, msg)
		log_main.add(msg)
		xform_proj_list = EMUtil.get_all_attributes(Tracker["constants"]["orgstack"], "xform.projection")
		for index_of_particle in xrange(len(xform_proj_list)):
			dp = xform_proj_list[index_of_particle].get_params("spider")
			xform_proj_list[index_of_particle] = [dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"], chunk_dict[index_of_particle]]
		write_text_row(xform_proj_list, Tracker["constants"]["partstack"])
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		msg  = "Extract xform.project parameters from stack"
		print(line, msg)
		log_main.add(msg)
	else:
		chunk_one = 0
		chunk_two = 0
	chunk_one = wrap_mpi_bcast(chunk_one, Blockdata["main_node"])
	chunk_two = wrap_mpi_bcast(chunk_two, Blockdata["main_node"])
	for element in chunk_one: chunk_dict[element] = 0
	for element in chunk_two: chunk_dict[element] = 1
	chunk_list 				= [chunk_one, chunk_two]
	Tracker["chunk_dict"] 	= chunk_dict
	Tracker["P_chunk_0"]   	= len(chunk_one)/float(total_stack)
	Tracker["P_chunk_1"]   	= len(chunk_two)/float(total_stack)
	vols = []
	# Reconstruction to determine the resolution in orignal data size
	Tracker["nxinit"] 		= Tracker["constants"]["nnxo"]
	Tracker["shrinkage"]	= float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	Tracker["directory"]	= Tracker["constants"]["masterdir"]
	
	if Blockdata["myid"] == Blockdata["main_node"]:
		msg  = "reconstruct 3-D volumes from two random subsets to calculate FSC"
		print(line, msg)
		log_main.add(msg)
		
	for procid in xrange(2):
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		#if Blockdata["myid"] == Blockdata["main_node"]: print(line, "Reconstruct volume  by particle from random subset %d"%procid)
		data = get_shrink_data_sorting(os.path.join(Tracker["constants"]["masterdir"],"chunk_%01d.txt"%procid), Tracker["constants"]["partstack"])
		mpi_barrier(MPI_COMM_WORLD)
		do3d_sorting(procid, data)
	mpi_barrier(MPI_COMM_WORLD)

	if( Blockdata["myid"] == Blockdata["nodes"][1]):  # It has to be 1 to avoid problem with tvol1 not closed on the disk
		#--  memory_check(Blockdata["myid"],"first node, before stepone")
		#  read volumes, shrink
		tvol0 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0.hdf"))
		tweight0 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_0.hdf"))
		tvol1 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_1.hdf"))
		tweight1 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_1.hdf"))
		Util.fuse_low_freq(tvol0, tvol1, tweight0, tweight1, 2*Tracker["constants"]["fuse_freq"])
		tag = 7007
		send_EMData(tvol1, Blockdata["nodes"][0], tag, MPI_COMM_WORLD)
		send_EMData(tweight1, Blockdata["nodes"][0], tag, MPI_COMM_WORLD)
		shrank0 	= stepone(tvol0, tweight0)
		send_EMData(shrank0, Blockdata["nodes"][0], tag, MPI_COMM_WORLD)
		del shrank0
		lcfsc = 0
		#--  memory_check(Blockdata["myid"],"first node, after stepone")
	
	elif( Blockdata["myid"] == Blockdata["nodes"][0]):

		#--  memory_check(Blockdata["myid"],"second node, before stepone")
		#  read volumes, shrink
		tag = 7007
		tvol1 		= recv_EMData(Blockdata["nodes"][1], tag, MPI_COMM_WORLD)
		tweight1 	= recv_EMData(Blockdata["nodes"][1], tag, MPI_COMM_WORLD)
		tvol1.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1} )
		shrank1 	= stepone(tvol1, tweight1)
		#  Get shrank volume, do fsc, send it to all
		shrank0 	= recv_EMData(Blockdata["nodes"][1], tag, MPI_COMM_WORLD)
		#  Note shrank volumes are Fourier uncentered.
		cfsc 		= fsc(shrank0, shrank1)[1]
		write_text_row(cfsc, os.path.join(Tracker["directory"], "fsc_driver.txt"))
		del shrank0, shrank1
		if(Tracker["nxinit"]<Tracker["constants"]["nnxo"]):
			cfsc  = cfsc[:Tracker["nxinit"]]
			for i in xrange(len(cfsc),Tracker["constants"]["nnxo"]//2+1): cfsc.append(0.0)
		lcfsc  = len(cfsc)							
		fsc05  = 0
		fsc143 = 0 
		for ifreq in xrange(len(cfsc)):	
			if cfsc[ifreq] <0.5: break
		fsc05  = ifreq - 1
		for ifreq in xrange(len(cfsc)):
			if cfsc[ifreq]<0.143: break
		fsc143 = ifreq - 1
		Tracker["constants"]["fsc143"] = fsc143
		Tracker["constants"]["fsc05"]  = fsc05
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		msg ="Initialization of sorting from data stack is done"
		print(line, msg)
		log_main.add(msg)
		#--  memory_check(Blockdata["myid"],"second node, after stepone")
	else:
		#  receive fsc
		lcfsc = 0
	Tracker = wrap_mpi_bcast(Tracker, Blockdata["nodes"][0], communicator = MPI_COMM_WORLD)
	# initialize some parameters for this case
	Tracker["delta"] = 1.5 #Tracker["constants"]["delta"]
	Tracker["an"]    = 1.5 #Tracker["constants"]["an"]
	Tracker["ts"]    = 0.5 #	Tracker["constants"]["ts"]
	Tracker["xr"]    = 0.4 #	Tracker["constants"]["xr"]
	return import_from_data_stack
	
#### Old rec3ds
def do3d_sorting_groups_4nn_iter(projdata, iteration = 0):
	global Tracker, Blockdata
	from   reconstruction import rec3D_two_chunks_MPI
	# For both CTF and noCTF 
	snr                 = 1.0
	npad                = 2
	fscmask             = None #= model_circle(int(Tracker["constants"]["radius"]*Tracker["shrinkage"]), Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"])
	#print("NXINIT!!!!", Tracker["nxinit"])
	fscc                = Tracker["number_of_groups"]*[None]
	
	if Blockdata["myid"] == Blockdata["main_node"]: 
		res_05   = Tracker["number_of_groups"]*[0]
		res_143  = Tracker["number_of_groups"]*[0]
	else:
		res_05   = 0
		res_143  = 0
	
	for iref in xrange(Tracker["number_of_groups"]):
		if not Tracker["constants"]["noctf"]: volref, fscc[iref] = rec3D_two_chunks_MPI(projdata, snr, Tracker["constants"]["symmetry"],  \
		fscmask, os.path.join(Tracker["directory"], "resolution_grp%03d_%03d.txt"%(iref, iteration)), Blockdata["myid"], Blockdata["main_node"],\
		 index = iref, npad = npad, finfo= None)
		else: volref, fscc[iref] = rec3D_MPI_noCTF(projdata, Tracker["constants"]["symmetry"], fscmask, \
		os.path.join(Tracker["directory"], "resolution_grp%03d_%03d.txt"%(iref, iteration)), Blockdata["myid"], Blockdata["main_node"],\
		 index = iref, npad = npad, finfo= None)
		if Blockdata["myid"] == Blockdata["main_node"]: volref.write_image(os.path.join(Tracker["directory"],"vol_grp%03d_iter%03d.hdf"%(iref,iteration)))
		mpi_barrier(MPI_COMM_WORLD)
	
	if Blockdata["myid"] == Blockdata["main_node"]:	
		for iref in xrange(Tracker["number_of_groups"]): 
			res_05[iref], res_143[iref]= get_res(fscc[iref][1])
	res_143  = wrap_mpi_bcast(res_143, Blockdata["main_node"])
	res_05   = wrap_mpi_bcast(res_05, Blockdata["main_node"])
	Tracker["fsc143"] = res_143
	Tracker["fsc05"]  = res_05
	return

def main():
	from optparse   import OptionParser
	from global_def import SPARXVERSION
	from EMAN2      import EMData
	from logger     import Logger, BaseLogger_Files
	import sys, os, time
	global Tracker, Blockdata
	from global_def import ERROR
	       
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack  outdir --refinement_dir=masterdir_of_sxmeridien --mask3D=mask.hdf --focus=binarymask.hdf  --radius=outer_radius " +\
	"  --sym=c1  --independent=indenpendent_runs  --number_of_images_per_group=number_of_images_per_group  --low_pass_filter=.25 "
	parser = OptionParser(usage,version=SPARXVERSION)
	#(options, args) 				= parser.parse_args(sys.argv[1:])
	#parser = OptionParser(usage,version=SPARXVERSION)
	# sorting options
	parser.add_option("--refinement_dir",                  type   ="string",        default ='',                       help="3-D refinement directory, ussually the master directory of sxmeridien")
	parser.add_option("--masterdir",                       type   ="string",        default ='',					   help="masterdir for sorting")
	parser.add_option("--refinement_method",               type   ="string",        default ='',					   help="refinement method used for 3-D reconstruction, [SPARX], [relion], or others")
	parser.add_option("--niter_for_sorting",               type   ="int",    		default = -1,					   help="selected iteration of refinement used for sorting")
	parser.add_option("--focus",                           type   ="string",        default ='',                       help="file name, the bineary 3D mask for focused clustering ")
	parser.add_option("--mask3D",                          type   ="string",        default ='',                       help="file name, the 3-D global mask for clustering ")
	parser.add_option("--low_pass_filter",                 type   ="float",         default =-1.0,					   help="frequency, User provided low-pass filter in absolute frequencies on orignal image size" )
	parser.add_option("--Kmeans_lpf",                      type   ="string",        default ='adaptive',               help="frequency, low_pass filter options for Kmeans clustering. The options are [adaptive[, [max], [min], [adhoc], and [avg]." )
	parser.add_option("--instack",                         type   ="string",        default ='',					   help="file name, data stack for sorting provided by user")
	parser.add_option("--radius",                          type   ="int",           default =-1,	                   help="particle radius in pixel for rotational correlation <nx-1 (set to the radius of the particle)")
	parser.add_option("--nxinit",						   type   ="int",           default =-1,					   help="integer number, user provided image size for sorting. Otherwise, program determines it from resolution" )
	#parser.add_option("--wn",                              type   ="int",           default =0,					       help="optimal window size for data processing. Reduce image size of original data by chopping eduge off")
	parser.add_option("--noctf",	                       action ="store_true",    default =False,                    help="do no ctf correction during clustring")
	parser.add_option("--sym",                             type   ="string",        default ='c1',                     help="point group symmetry of the structure")
	parser.add_option("--nindependent",                    type   ="int",           default = 3,                       help="number of independent run for EQkmeans clustering, an odd number larger than 2")
	parser.add_option("--number_of_images_per_group",      type   ="int",           default =1000,                     help="number of images in a group")
	parser.add_option("--smallest_group",				   type   ="int",           default =500,					   help="minimum number of members for being identified as a group")
	parser.add_option("--interpolation",                   type   ="string",        default ='4nn',                    help="3-D reconstruction interpolation method, either [trl] or [4nn]")
	parser.add_option("--comparison_method",               type   ="string",        default ='cross',                  help="comparision method, either cross-correlaton coefficients [cross] or Euclidean distance [eucd] ")
	(options, args) 				= parser.parse_args(sys.argv[1:])

	## Sanity check
	from utilities import bcast_number_to_all
	if options.interpolation =='4nn' and options.comparison_method == 'eucd':
		ERROR("interpolation 4nn is compatable with comparison_method cross, not eucd", "sort3d", 1)
		exit()
		
	if  options.interpolation =='trl' and os.path.exists(options.focus) and options.comparison_method == 'eucd':
		ERROR("interpolation trl and comparison method eucd are incompatable with focus mask. Try either interpolation 4nn or comparison method cross to include focus mask in sorting", "sort3d", 1)
		exit()
			
	if options.nindependent<=2:
		ERROR("nindependent has to be an odd number larger than 2", "sort3d", 1)
		exit()
	
	if options.refinement_dir:
		if not os.path.exists(options.refinement_dir): ERROR("The specified refinement_dir does not exist", "sort3d", 1)
			
	if options.focus:
		if not os.path.exists(options.focus): 
			ERROR("The specified focus mask file does not exist", "sort3d", 1)
			exit()
			
	if options.mask3D:
		if not os.path.exists(options.mask3D): 
			ERROR("The specified mask3D file does not exist", "sort3d", 1)
			exit()
			
	if options.number_of_images_per_group <= options.smallest_group: 
		ERROR("number_of_images_per_group should be way larger than smallest_group", "sort3d", 1)
		exit()
	

					
	#--- Fill input parameters into dictionary named after Constants
	Constants		                         = {}
	Constants["orgstack"]                    = options.instack
	Constants["masterdir"]                   = options.masterdir
	Constants["refinement_dir"]              = options.refinement_dir
	Constants["niter_for_sorting"]           = options.niter_for_sorting
	Constants["mask3D"]                      = options.mask3D
	Constants["focus3Dmask"]                 = options.focus
	Constants["indep_runs"]                  = options.nindependent
	Constants["number_of_images_per_group"]  = options.number_of_images_per_group
	Constants["smallest_group"]      		 = options.smallest_group
	Constants["not_included_unacc"]          = False  #options.not_included_unacc
	Constants["noctf"]                 		 = options.noctf
	Constants["radius"]              		 = options.radius
	Constants["symmetry"]                    = options.sym
	Constants["low_pass_filter"]             = options.low_pass_filter # enforced low_pass_filter
	Constants["nxinit"]                      = options.nxinit
	Constants["seed"]                        = -1
	Constants["upscale"]                     = 0.5 #
	#Constants["wn"]                          = options.wn
	Constants["interpolation"]               = options.interpolation
	Constants["comparison_method"]           = options.comparison_method
	Constants["fuse_freq"]                   = 40  # Now in A, convert to absolute before using
	Constants["Kmeans_lpf"]                  = options.Kmeans_lpf
	
	# -------------------------------------------------------------
	#
	# Create and initialize Tracker dictionary with input options  # State Variables
	
	Tracker							               = {}
	Tracker["constants"]			               = Constants
	if Tracker["constants"]["mask3D"]:	           Tracker["mask3D"] = Tracker["constants"]["mask3D"]
	else:								           Tracker["mask3D"] = None
	Tracker["radius"]						       = Tracker["constants"]["radius"]
	Tracker["upscale"]						       = Tracker["constants"]["upscale"]
	Tracker["applyctf"]						       = False  #  Should the data be premultiplied by the CTF.  Set to False for local continuous.
	Tracker["nxinit"]						       = Tracker["constants"]["nxinit"]
	Tracker["lowpass"]						       = Tracker["constants"]["low_pass_filter"]
	
	##<<<--options for advanced users:  
	Tracker["total_sort3d_indepent_run"]	       = 2
	Tracker["total_number_of_iterations"] 	       = 35
	Tracker["EQKmeans_stopercent"]			       = 3.0 # change the converge of EQKmeans
	Tracker["Kmeans_stopercent"]                   = 1.0 # change the converge of Kmeans
	Tracker["total_iter_rsort"]                    = 2
	Tracker["falloff"]						       = 0.1
	Tracker["aa"]				                   = 0.1
	Tracker["clean_volumes"]                       = True
	Tracker["constants"]["total_sort3d_iteration"] = 2 
	
	if options.interpolation =='trl' and os.path.exists(options.focus) and options.comparison_method =='cross':
		Tracker["total_number_of_iterations"] +=10
	if options.interpolation =='4nn': 
		Tracker["Kmeans_stopercent"]   = 0.1
		Tracker["EQKmeans_stopercent"] = 1.0
	
	###--------------------------------------------------------------------------------------------
	#  three sorting scenarios 
	# 1. given data stack and xform.projection/ctf in header; 
	# 2. directly import from SPARX refinement;
	# 3. directly import from relion refinement;
	
	#<<<---------------------->>>imported functions<<<---------------------------------------------

	from statistics 	import k_means_match_clusters_asg_new,k_means_stab_bbenum
	from applications 	import recons3d_n_MPI, mref_ali3d_MPI, Kmref_ali3d_MPI
	from utilities 		import get_im,bcast_number_to_all, write_text_file,read_text_file,wrap_mpi_bcast, get_params_proj, write_text_row
	from utilities 		import get_initial_ID, print_upper_triangular_matrix, print_a_line_with_timestamp
	from utilities 		import get_resolution_mrk01,partition_to_groups,partition_independent_runs,get_outliers
	from utilities 		import merge_groups, save_alist, margin_of_error, get_margin_of_error, get_ali3d_params
	from utilities 		import counting_projections, unload_dict, load_dict, get_stat_proj, get_number_of_groups, recons_mref
	from utilities 		import apply_low_pass_filter, get_groups_from_partition, get_number_of_groups, get_complementary_elements_total, update_full_dict
	from utilities 		import count_chunk_members, set_filter_parameters_from_adjusted_fsc, adjust_fsc_down, get_two_chunks_from_stack
	from utilities 		import get_input_from_string, cmdexecute
	from filter			import filt_tanl
	from time           import sleep
	from logger         import Logger,BaseLogger_Files
	import user_functions
	import string
	from string         import split, atoi, atof
	import json
	
	####------------------------------------------------------------------
	restart		= 0
	keepgoing 	= 1
	# 	Create Master directory
	line      = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	masterdir = Tracker["constants"]["masterdir"]
	if(Blockdata["myid"] == Blockdata["main_node"]):
		print(line, "Sort3d starts")
		if not masterdir:
			timestring = strftime("_%d_%b_%Y_%H_%M_%S", localtime())
			masterdir ="sort3d"+timestring
			cmd="{} {}".format("mkdir", masterdir)
			junk = cmdexecute(cmd)
		else:
			if os.path.exists(masterdir): restart = 1
			else:
				cmd="{} {}".format("mkdir", masterdir)
				junk = cmdexecute(cmd)
		li =len(masterdir)
		#print_dict(Tracker["constants"],"Permanent sorting settings from input options")
	else:
		li = 0
	restart											= bcast_number_to_all(restart, Blockdata["main_node"])
	li                                              = mpi_bcast(li,1,MPI_INT,Blockdata["main_node"],MPI_COMM_WORLD)[0]
	masterdir										= mpi_bcast(masterdir,li,MPI_CHAR,Blockdata["main_node"],MPI_COMM_WORLD)
	masterdir                                       = string.join(masterdir,"")
	Tracker["constants"]["masterdir"]				= masterdir
	Tracker["constants"]["chunk_0"]					= os.path.join(Tracker["constants"]["masterdir"],"chunk_0.txt")
	Tracker["constants"]["chunk_1"]					= os.path.join(Tracker["constants"]["masterdir"],"chunk_1.txt")              

	log_main = Logger(BaseLogger_Files())
	log_main.prefix = Tracker["constants"]["masterdir"]+"/"
	
	import_from_relion_refinement 	= 0
	import_from_sparx_refinement 	= 0
	import_from_data_stack			= 0
	total_stack						= 0
	
	while not os.path.exists(Tracker["constants"]["masterdir"]):
		print("Node ", Blockdata["myid"], "  waiting...", Tracker["constants"]["masterdir"])
		sleep(1)
	mpi_barrier(MPI_COMM_WORLD)
	
	######### collect refinement information 
	try:    nxinit = Tracker["nxinit"] 
	except: Tracker["nxinit"] = -1				
	Tracker["constants"]["orgres"]				= 0.0
	Tracker["constants"]["refinement_delta"]	= 0.0
	Tracker["constants"]["refinement_ts"]		= 0.0
	Tracker["constants"]["refinement_xr"]		= 0.0
	Tracker["constants"]["refinement_an"]		= 0.0
	Tracker["constants"]["selected_iter"]		= 0

	if options.refinement_method =="SPARX": # Senario one
		import_from_sparx_refinement = get_input_from_sparx_ref3d(log_main)
		
	elif options.refinement_method =="relion":# Senario two
		import_from_relion_refinement = get_input_from_relion_ref3d(log_main)
		
	else:  # Senario three, sorting from a given data stack, general cases
		import_from_data_stack = get_input_from_datastack(log_main)
		 
	###<<<----------------------->>>>>  SORT3D MAIN PROGRAM <<<<<---------------------------------------------# For all cases
	
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if(Blockdata["myid"] == Blockdata["main_node"]):
		print(line, "Sort3d main program")
		log_main.add("---------------->>>SPARX sort3d<<<---------------------- ")
		log_main.add("The shell line command:")
		line=""
		for a in sys.argv:	line +=(a+" ")
		log_main.add(line)
		log_main.add("Sort3d master directory: %s"%Tracker["constants"]["masterdir"])
		print_dict(Tracker["constants"],"Permanent sorting settings after initialization")
	mpi_barrier(MPI_COMM_WORLD)
	
	####--->>>>>> mask preparation
	if Tracker["constants"]["mask3D"]: 		Tracker["mask3D"]  = os.path.join(Tracker["constants"]["masterdir"],"smask.hdf")
	else:									Tracker["mask3D"]  = None
	if Tracker["constants"]["focus3Dmask"]: Tracker["focus3D"] = Tracker["constants"]["focus3Dmask"]
	else:									Tracker["focus3D"] = None

	if(Blockdata["myid"] == Blockdata["main_node"]):
		bad_focus3Dmask = 0
		if Tracker["constants"]["focus3Dmask"]:
			try: 
				 focusmask = get_im(Tracker["constants"]["focus3Dmask"])
				 st = Util.infomask(binarize(focusmask), None, True)
				 if(st[0] == 0.0): bad_focus3Dmask = 1
				 else:             focusmask.write_image(Tracker["focus3D"])
			except:  bad_focus3Dmask = 1
	else: bad_focus3Dmask = 0
	bad_focus3Dmask = bcast_number_to_all(bad_focus3Dmask,	source_node =  Blockdata["main_node"])
	if bad_focus3Dmask:  ERROR("Incorrect focused mask, after binarize all values zero","sxsort3d.py", 1, Blockdata["myid"])
	
	if(Blockdata["myid"] == Blockdata["main_node"]):
		bad_3Dmask = 0
		if Tracker["constants"]["mask3D"]:
			try: 
				mask3D = get_im(Tracker["constants"]["mask3D"])
				st = Util.infomask(binarize(mask3D), None, True)
				if (st[0] ==0.0): bad_3Dmask = 1
				else: mask3D.write_image(Tracker["mask3D"]) 
			except: bad_3Dmask = 1
	else: bad_3Dmask = 0
	bad_3Dmask  = bcast_number_to_all(bad_focus3Dmask,	source_node =  Blockdata["main_node"])
	if bad_3Dmask:  ERROR("Incorrect 3D mask", "sxsort3d.py", 1, Blockdata["myid"])
	mpi_barrier(MPI_COMM_WORLD)

	############################################################################################
	###<<<------ Determine the image size ### reset nxinit
	if Tracker["nxinit"] < 0:
		try:     fsc143 = Tracker["constants"]["fsc143"]
		except:  fsc143 = 0
		try:     fsc05  = Tracker["constants"]["fsc05"]
		except:  fsc05  = 0
		if   fsc143  !=0: Tracker["nxinit"] = (Tracker["constants"]["fsc143"]+5)*2
		elif fsc05   !=0: Tracker["nxinit"] = (Tracker["constants"]["fsc05"]+12)*2
		else: ERROR("Incorrect nxinit, and check the input information", "sxsort3d.py", 1, Blockdata["myid"])
	if Tracker["constants"]["nxinit"] > 0: Tracker["nxinit"] = Tracker["constants"]["nxinit"] # User defined nxinit
	Tracker["currentres"] = float(Tracker["constants"]["fsc05"])/float(Tracker["nxinit"])
	
	##################---------------<<<shrinkage, current resolution  <<<<<<------------------------------------------
	Tracker["total_stack"] = Tracker["constants"]["total_stack"]	
	Tracker["shrinkage"]   = float(Tracker["nxinit"])/Tracker["constants"]["nnxo"]
	Tracker["radius"]      = Tracker["constants"]["radius"]*Tracker["shrinkage"]
	Tracker["fuse_freq"]   = Tracker["constants"]["pixel_size"]/Tracker["constants"]["fuse_freq"]
	
	###################--------------------------->>>>> 3-D masks <<<<<--------------------------------------------------------------------
	
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if(Blockdata["myid"] == Blockdata["main_node"]):
		print_dict(Tracker["constants"],"Permanent sorting settings from input options")
		fout = open(os.path.join(Tracker["constants"]["masterdir"], "Tracker.json"),'w')
		json.dump(Tracker, fout)
		fout.close()
		msg = "---------->>>The first phase <<<----------- "
		log_main.add(msg)
		print(line, msg)
		msg = "orgstack: %s"%Tracker["constants"]["orgstack"]
		print(line, msg)
		log_main.add(msg)
		msg ="image comparison method: %s"%Tracker["constants"]["comparison_method"]
		print(line, msg)
		log_main.add(msg)
		msg ="3-D reconstruction method: %s"%Tracker["constants"]["interpolation"]
		print(line, msg)
		log_main.add(msg)
		if Tracker ["constants"]["focus3Dmask"]:
			msg ="User provided focus mask file:     %s"%Tracker ["constants"]["focus3Dmask"]
			print(line, msg)
			log_main.add(msg)
		Tracker["full_list"] = read_text_file(os.path.join(Tracker["constants"]["masterdir"],"indexes.txt"), -1) # could have one or two columns
		if    len(Tracker["full_list"]) == 2: Tracker["full_list"] = Tracker["full_list"][1] # take just one column
		elif  len(Tracker["full_list"]) == 1: Tracker["full_list"] = Tracker["full_list"][0]
		else: ERROR("The original particle ID for sorting has wrong format", "sxsort3d.py", 1, Blockdata["main_node"])
	else: Tracker["full_list"] = 0
	Tracker["full_list"]       = wrap_mpi_bcast(Tracker["full_list"], Blockdata["main_node"], MPI_COMM_WORLD)	
	Tracker["shrinkage"]       = float(Tracker["nxinit"])/Tracker["constants"]["nnxo"]
	if(Blockdata["myid"] == Blockdata["main_node"]): print_dict(Tracker,"Current sorting settings")


	#  READ AND PREPREPARED DATA  PAP


	###<<<-------sort3d starts here
	for indep_sort3d in xrange(Tracker["total_sort3d_indepent_run"]):		
		sorting                = {}
		sorting["Unaccounted"] = None 
		sorting["Accounted"]   = []
		partid_file            = os.path.join(Tracker["constants"]["masterdir"],"indexes.txt")
		line                   = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		Tracker["indep_sort3d_dir"] = os.path.join(Tracker["constants"]["masterdir"], "sort3d_run%d"%indep_sort3d)
		if Blockdata["myid"] == Blockdata["main_node"]:
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"	
			junk = cmdexecute("{} {}".format("mkdir",Tracker["indep_sort3d_dir"]))
			msg = "---------->>>Independent sort3d  %d<<<----------- "%indep_sort3d
			print(line, msg)
			log_main.add(msg)
			log_main.add("nnxo : %d"%Tracker["constants"]["nnxo"])
			log_main.add("Current resolution: %6.3f  absolute unit(maximum is 0.5)    1./%6.2f Angstrom "%(Tracker["currentres"]*Tracker["shrinkage"], Tracker["constants"]["pixel_size"]/Tracker["currentres"]/Tracker["shrinkage"]))
			if Tracker["constants"]["low_pass_filter"]!=-1.0:  log_main.add("User provided enforced low_pass_filter: %f"%Tracker["constants"]["low_pass_filter"])
			else:                                              log_main.add("low_pass_filter is calculated adaptively during clustering ")
			if Tracker["mask3D"]:
				msg = "User provided 3-D mask:  %s"%Tracker["constants"]["mask3D"]
				log_main.add(msg)
				print(line, msg)
			
			if import_from_relion_refinement:
				msg =  "Sorting is initiated from relion refinement"
				print(line, msg)
				log_main.add(msg)
				
			elif import_from_sparx_refinement:
				msg = "Sorting is initiated from SPARX refinement"
				print(line, msg)
				log_main.add(msg)
				
			else:
				msg= "Sorting is initiated from  data stack %s"%Tracker["constants"]["orgstack"]                              
				print(line, msg)
				log_main.add(msg)
				
			sorting["total"] = read_text_file(partid_file, -1)
			if len(sorting["total"])>1: sorting["total"] = sorting["total"][1]
			else:                       sorting["total"] = sorting["total"][0]
		else:                           sorting["total"] = 0
		sorting["total"] = wrap_mpi_bcast(sorting["total"], Blockdata["main_node"]) # total number of records in indexes.txt file 
		
		if Blockdata["myid"] == Blockdata["main_node"]:
			Tracker["number_of_images_per_group"] = Tracker["constants"]["number_of_images_per_group"]
			Tracker["total_stack"] 				  = len(sorting["total"]) # start from beginning
			msg = "number_of_images_per_group:  %d"%Tracker["number_of_images_per_group"]
			print(line, msg)
			log_main.add(msg)
			Tracker["number_of_groups"]   =	get_number_of_groups(Tracker["total_stack"],Tracker["number_of_images_per_group"])
		Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], MPI_COMM_WORLD)
		final_list = do_EQKmeans_nways_clustering(Tracker["indep_sort3d_dir"], partid_file, \
		os.path.join(Tracker["constants"]["masterdir"],"refinement_parameters.txt"), sorting["Accounted"], log_main)
		sorting["Unaccounted"] = Tracker["unaccounted_list"]
		if Blockdata["myid"] == Blockdata["main_node"]:
			# sorting.json contains [Accounted] lists and one Unaccounted list
			fout = open(os.path.join(Tracker["indep_sort3d_dir"], "sorting.json"),'w')
			json.dump(sorting, fout)
			fout.close()
		mpi_barrier(MPI_COMM_WORLD)
		
	#########################################################<<<---rsort
	
	if( Blockdata["myid"] == Blockdata["main_node"] ):
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		msg = "---------->>>The second phase <<<----------- "
		log_main.add(msg)
		print(line, msg)
	iter_rsort    = 0
	while iter_rsort< Tracker["total_iter_rsort"]:
		ptp        = []
		for indep_sort3d in xrange(Tracker["total_sort3d_indepent_run"]):
			indep_sort3d_dir = os.path.join(Tracker["constants"]["masterdir"], "sort3d_run%d"%indep_sort3d)
			if os.path.exists(os.path.join(indep_sort3d_dir, "sorting.json")):				
				if(Blockdata["myid"] == Blockdata["main_node"]):
					fout       = open(os.path.join(indep_sort3d_dir, "sorting.json"), "r")
					res_sort3d = convert_json_fromunicode(json.load(fout))
					fout.close()
					merged_classes = split_partition_into_clusters(res_sort3d["Accounted"])
					if len(res_sort3d["Unaccounted"])>0: merged_classes.append(res_sort3d["Unaccounted"])
				else:
					res_sort3d     = 0
					merged_classes = 0
				res_sort3d     = wrap_mpi_bcast(res_sort3d,  Blockdata["main_node"], MPI_COMM_WORLD)
				merged_classes = wrap_mpi_bcast(merged_classes, Blockdata["main_node"], MPI_COMM_WORLD)
				sptp           = prep_ptp_single(merged_classes, res_sort3d["total"])
				ptp.append(sptp)
			else:	ERROR("sorting results do not exist", "sxsort3d.py", 1, Blockdata["myid"])
		mpi_barrier(MPI_COMM_WORLD)
		accounted_list, unaccounted_list, new_index = do_two_way_comparison_single(ptp[0], ptp[1], len(res_sort3d["total"]))
		Tracker["unaccounted_list"] = []
		for index in xrange(len(unaccounted_list)): Tracker["unaccounted_list"].append(res_sort3d["total"][unaccounted_list[index]]) # always check
		Tracker["accounted_list"]   = []
		for index in xrange(len(accounted_list)): Tracker["accounted_list"].append([new_index[index][0], res_sort3d["total"][new_index[index][1]]]) # always check
		Tracker["directory"] = os.path.join(Tracker["constants"]["masterdir"], "rsort%d"%iter_rsort)
	
		if Blockdata["myid"] == Blockdata["main_node"]:
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			msg = "---------->>>rsort%d<<<----------- "%iter_rsort
			log_main.add(msg)
			print(line, msg)
			msg = "Summary of two sort3d runs"
			log_main.add(msg)
			print(line, msg)
			msg = "Accounted:    %d    Unaccounted:    %d"%(len(Tracker["accounted_list"]), len(Tracker["unaccounted_list"]))
			log_main.add(msg)
			print(line, msg)
			cmd = "{} {}".format("mkdir",Tracker["directory"])
			junk = cmdexecute(cmd)
			write_text_file(Tracker["unaccounted_list"], os.path.join(Tracker["directory"], "Unaccounted.txt"))
			write_text_row(Tracker["accounted_list"], os.path.join(Tracker["directory"], "Accounted.txt"))
		mpi_barrier(MPI_COMM_WORLD)
			
		partid_file = os.path.join(os.path.join(Tracker["directory"], "Unaccounted.txt"))
		Tracker["total_stack"]	     = len(Tracker["unaccounted_list"])
		Tracker["number_of_groups"]  = get_number_of_groups(Tracker["total_stack"],Tracker["number_of_images_per_group"])
		
		if Tracker["number_of_groups"]>1: # continue N-ways_clustering on unaccounted particles 
			sorting                      = {}
			sorting["Accounted"]         = [Tracker["accounted_list"]] # keep the accounted
			final_sort                   = []
			final_list = do_EQKmeans_nways_clustering(os.path.join(Tracker["constants"]["masterdir"], "rsort%d"%iter_rsort), partid_file, \
			os.path.join(Tracker["constants"]["masterdir"],"refinement_parameters.txt"), sorting["Accounted"], log_main)
			for a in sorting["Accounted"]:
				final_sort.append(a)		
			if (Blockdata["myid"] == Blockdata["main_node"]):
				if len(Tracker["unaccounted_list"])> 0:
					for iparticle in xrange(len(Tracker["unaccounted_list"])):
						Tracker["unaccounted_list"][iparticle] = [0, Tracker["unaccounted_list"][iparticle]]
					final_sort.append(Tracker["unaccounted_list"])
				else:
					msg = "empty unaccounted_list is found"
					print(line, msg)
					log_main.add(msg)
				indexed_particle_list, Tracker["number_of_groups"] = merge_original_id_lists(final_sort)
				write_text_row(indexed_particle_list, os.path.join(Tracker["constants"]["masterdir"], "rsort%d"%iter_rsort, "index_for_Kmeans.txt"))
				line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
				msg  = "include %d unaccounted as one cluster  "%len(Tracker["unaccounted_list"])
				print(line, msg)
				log_main.add(msg)
	
		else: # carry out with unaccounted ones also included into the final exhaustive K-means
			if (Blockdata["myid"] == Blockdata["main_node"]):
				final_sort                   = []
				final_sort.append(Tracker["accounted_list"])
				# group-index unaccounted_list
				if len(Tracker["unaccounted_list"])>0:
					for iparticle in xrange(len(Tracker["unaccounted_list"])):
						Tracker["unaccounted_list"][iparticle] = [0, Tracker["unaccounted_list"][iparticle]]
					final_sort.append(Tracker["unaccounted_list"])
				else:
					msg = "empty unaccounted_list is found"
					print(line, msg)
					log_main.add(msg)
				indexed_particle_list, Tracker["number_of_groups"] = merge_original_id_lists(final_sort)
				write_text_row(indexed_particle_list, os.path.join(Tracker["constants"]["masterdir"], "rsort%d"%iter_rsort, "index_for_Kmeans.txt"))
				line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
				msg = "include %d unaccounted as one cluster   "%len(Tracker["unaccounted_list"])
				print(line, msg)
				log_main.add(msg)
	
		Tracker["number_of_groups"] = bcast_number_to_all(Tracker["number_of_groups"], Blockdata["main_node"], MPI_COMM_WORLD)				
		Tracker["directory"]        = os.path.join(Tracker["constants"]["masterdir"],"rsort%d"%iter_rsort,"Kmeans")
		if (Blockdata["myid"] == Blockdata["main_node"]):
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			msg = "Exhaustive Kmeans run on %d clusters "%Tracker["number_of_groups"]
			print(line, msg)
			log_main.add(msg)
			cmd = "{} {}".format("mkdir",Tracker["directory"])
			junk = cmdexecute(cmd)
		mpi_barrier(MPI_COMM_WORLD)
		
		final_list, num_in_groups = mref_ali3d_Kmeans_remove_small_groups(os.path.join(Tracker["constants"]["masterdir"], "rsort%d"%iter_rsort, "index_for_Kmeans.txt"), \
		os.path.join(Tracker["constants"]["masterdir"],"refinement_parameters.txt"), Tracker["clean_volumes"])
		if(Blockdata["myid"] == Blockdata["main_node"]):
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			for iref in xrange(len(num_in_groups)):
				msg = " group  %d:    %d "%(iref, num_in_groups[iref])
				log_main.add(msg)
				print(line, msg)
			log_main.add("The rsort%d final Kmeans on %d"%(iter_rsort, len(final_list)))
			log_main.add("The total number:   %d"%Tracker["constants"]["total_stack"])
			print("The final accounted: %d"%len(final_list))
			msg = "rsort%d finishes"%iter_rsort
			print(line, msg)
			log_main.add(msg)
			write_text_row(final_list, os.path.join(Tracker["constants"]["masterdir"], "rsort%d_results.txt"%iter_rsort))
			clusters = split_partition_into_clusters([final_list])
		else: clusters = 0
		clusters     = wrap_mpi_bcast(clusters, Blockdata["main_node"], MPI_COMM_WORLD)
		iter_rsort  +=1
	# rsort final comparison
	mpi_barrier(MPI_COMM_WORLD)
	if(Blockdata["myid"] == Blockdata["main_node"]):
		line    = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		tlist   = read_text_file(os.path.join(Tracker["constants"]["masterdir"],"indexes.txt"))
		ptp     = []
		for irsort in xrange(Tracker["total_iter_rsort"]):
			r1 = read_text_row(os.path.join(Tracker["constants"]["masterdir"], "rsort%d_results.txt"%irsort))
			merged_classes = split_partition_into_clusters([r1])
			sptp           = prep_ptp_single(merged_classes, tlist)
			ptp.append(sptp)
		accounted_list, unaccounted_list, new_index = do_two_way_comparison_single(ptp[0], ptp[1], len(tlist))
		# export_sorting_results(clusters)
		msg  = " total number of accounted for two rsort:  %d percentage:  %5.2f "%(len(accounted_list), float(len(accounted_list))/float(len(tlist))*100.0)
		print(line, msg)
		log_main.add(msg)
		Tracker["accounted_list"] = []
		for index in xrange(len(accounted_list)): 
			Tracker["accounted_list"].append([new_index[index][0], tlist[new_index[index][1]]]) # always check
		rsort_clusters, new_partition = split_partition_into_ordered_clusters(Tracker["accounted_list"])
		Tracker["accounted_list"] = new_partition
		write_text_row(Tracker["accounted_list"], os.path.join(Tracker["constants"]["masterdir"], "final_partition.txt"))
	else:
		rsort_clusters = 0
	rsort_clusters = wrap_mpi_bcast(rsort_clusters, Blockdata["main_node"], MPI_COMM_WORLD)
	
	if(Blockdata["myid"] == Blockdata["main_node"]):# save clusters with respect to the initial indexes
		for icluster in xrange(len(rsort_clusters)):
			write_text_file(rsort_clusters[icluster], os.path.join(Tracker["constants"]["masterdir"],"Cluster%d.txt"%icluster))
			Tracker["number_of_groups"] = len(rsort_clusters)
	else:   Tracker["number_of_groups"] = 0
	Tracker["number_of_groups"]         = bcast_number_to_all(Tracker["number_of_groups"], Blockdata["main_node"], MPI_COMM_WORLD)
	Tracker["directory"]                = Tracker["constants"]["masterdir"]

	if Tracker["constants"]["interpolation"]=="4nn":  
		data = get_shrink_data_sorting(os.path.join(Tracker["constants"]["masterdir"], "final_partition.txt"), \
		os.path.join(Tracker["constants"]["masterdir"],"refinement_parameters.txt"), \
		return_real = True, preshift = True, apply_mask = False)
		do3d_sorting_groups_4nn_iter(data, 0)
		del data
	else:
		data = get_shrink_data_sorting(os.path.join(Tracker["constants"]["masterdir"], "final_partition.txt"), \
		os.path.join(Tracker["constants"]["masterdir"],"refinement_parameters.txt"), \
		return_real = False, preshift = True, apply_mask = False)
		do3d_sorting_groups_trl_iter(data, 0)
		del data
				
	if(Blockdata["myid"] == Blockdata["main_node"]):
		line    = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		cmd = "{} {}".format("rm -rf", os.path.join(Tracker["directory"], "tempdir"))
		if os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
			junk = cmdexecute(cmd)
		msg = "Final sort3d summary"
		print(line, msg)
		log_main.add(msg)
		for icluster in xrange(Tracker["number_of_groups"]):
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			if Tracker["constants"]["interpolation"]=="4nn":  
				cfsc          = read_text_file(os.path.join(Tracker["constants"]["masterdir"],"resolution_grp%03d_000.txt"%icluster), -1)
				res05, res143 = get_res(cfsc[1])
			else:
				cfsc          = read_text_file(os.path.join(Tracker["constants"]["masterdir"],"fsc_driver_grp%03d_iter000.txt"%icluster), -1)
				res05, res143 = get_res(cfsc[0])    
			vol = get_im(os.path.join(Tracker["constants"]["masterdir"],"vol_grp%03d_iter000.hdf"%icluster))
			vol = filt_tophatl(vol, float(res05)/float(Tracker["nxinit"]))
			vol.write_image(os.path.join(Tracker["constants"]["masterdir"],"volf_grp%03d_iter000.hdf"%icluster))
			msg = "group %3d   number of images: %8d     FSC05    %f     FSC143   %f"%(icluster, len(rsort_clusters[icluster]), \
			float(res05)/float(Tracker["nxinit"]), float(res143)/float(Tracker["nxinit"]))
			print(line, msg)
			log_main.add(msg)
	mpi_barrier(MPI_COMM_WORLD)
	from mpi import mpi_finalize
	mpi_finalize()
	exit()
if __name__ == "__main__":
	main()
