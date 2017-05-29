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
nproc     = mpi_comm_size(MPI_COMM_WORLD)
myid      = mpi_comm_rank(MPI_COMM_WORLD)
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
if(Blockdata["no_of_groups"] > 1):
	Blockdata["node_volume"] = [Blockdata["no_of_groups"]-2, Blockdata["no_of_groups"]-1]
	#Blockdata["nodes"] = [Blockdata["no_of_groups"]-2, Blockdata["no_of_groups"]-1]  # For 3D stuff take last two nodes
else: 
	Blockdata["node_volume"] = [0, 0]
#  We need two CPUs for processing of volumes, they are taken to be main CPUs on each volume
#  We have to send the two myids to all nodes so we can identify main nodes on two selected groups.
if(Blockdata["no_of_groups"] > 1): Blockdata["main_shared_nodes"] = [Blockdata["node_volume"][0]*Blockdata["no_of_processes_per_group"],Blockdata["node_volume"][1]*Blockdata["no_of_processes_per_group"]]
else:  Blockdata["main_shared_nodes"] = [0, 1]
Blockdata["nproc_previous"]  = 0
# End of Blockdata: sorting requires at least three nodes, and the used number of nodes be integer times of three
global_def.BATCH = True
global_def.MPI   = True
global _proc_status, _scale, is_unix_cluster
try:			
	_proc_status = '/proc/%d/status' % os.getpid()
	_scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,'KB': 1024.0, 'MB': 1024.0*1024.0}
	is_unix_cluster = True
except:
	if Blockdata["myid"]==Blockdata["main_node"]:print("Not a unix machine")
	is_unix_cluster = False
	
def create_subgroup():
	# select a subset of myids to be in subdivision
	if( Blockdata["myid_on_node"] < Blockdata["ncpuspernode"] ): submyids = [Blockdata["myid"]]
	else:  submyids = []
	submyids = wrap_mpi_gatherv(submyids, Blockdata["main_node"], MPI_COMM_WORLD)
	submyids = wrap_mpi_bcast(submyids, Blockdata["main_node"], MPI_COMM_WORLD)
	#if( Blockdata["myid"] == Blockdata["main_node"] ): print(submyids)
	world_group = mpi_comm_group(MPI_COMM_WORLD)
	subgroup = mpi_group_incl(world_group,len(submyids),submyids)
	Blockdata["subgroup_comm"] = mpi_comm_create(MPI_COMM_WORLD, subgroup)
	mpi_barrier(MPI_COMM_WORLD)
	Blockdata["subgroup_size"] = -1
	Blockdata["subgroup_myid"] = -1
	if (MPI_COMM_NULL != Blockdata["subgroup_comm"]):
		Blockdata["subgroup_size"] = mpi_comm_size(Blockdata["subgroup_comm"])
		Blockdata["subgroup_myid"] = mpi_comm_rank(Blockdata["subgroup_comm"])
	#  "nodes" are zero nodes on subgroups on the two "node_volume" that compute backprojection
	Blockdata["nodes"] = [Blockdata["node_volume"][0]*Blockdata["ncpuspernode"], Blockdata["node_volume"][1]*Blockdata["ncpuspernode"]]
	mpi_barrier(MPI_COMM_WORLD)
	return
	
def check_mpi_settings(log):
	global Tracker, Blockdata
	from   utilities import wrap_mpi_bcast, read_text_file, bcast_number_to_all
	import os
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	current_mpi_settings_is_bad = 0
	if(Blockdata["myid"] == Blockdata["main_node"]):
		fsc_refinement = read_text_file(os.path.join(Tracker["constants"]["masterdir"], "fsc_global.txt"))
		q = float(Tracker["constants"]["img_per_grp"])/float(Tracker["constants"]["total_stack"])
		for ifreq in xrange(len(fsc_refinement)): fsc_refinement[ifreq] = fsc_refinement[ifreq]*q/(1.-fsc_refinement[ifreq]*(1.-q))
		res = 0.0
		for ifreq in xrange(len(fsc_refinement)):
			if fsc_refinement[ifreq]<0.143: break
		res = float(ifreq)/2./float(len(fsc_refinement))
		nxinit = int(2.*res*Tracker["constants"]["nnxo"])
		del fsc_refinement
	else: nxinit =0
	nxinit = bcast_number_to_all(nxinit, Blockdata["main_node"], MPI_COMM_WORLD)	
	sys_required_mem = 1.0*Blockdata["no_of_processes_per_group"]
	if( Blockdata["myid"] == Blockdata["main_node"]):
		msg = "------->>>>>>>Check memory and mpi settings<<<<-------------"
		log.add(msg)
		print(line, msg)
		msg ="cpus number: %5d  node number:  %5d  cpu number per group:  %5d"%(Blockdata["nproc"], Blockdata["no_of_groups"], Blockdata["no_of_processes_per_group"])
		log.add(msg)
		print(line, msg)
	try:
		image_org_size = Tracker["constants"]["nnxo"]
		image_in_core_size = nxinit
		ratio = float(nxinit)/float(image_org_size)
		raw_data_size = float(Tracker["constants"]["total_stack"]*image_org_size*image_org_size)*4.0/1.e9
		raw_data_size_per_node = float(Tracker["constants"]["total_stack"]*image_org_size*image_org_size)*4.0/1.e9/Blockdata["no_of_groups"]
		sorting_data_size_per_node = raw_data_size_per_node + 2.*raw_data_size_per_node*ratio**2
		volume_size_per_node = (4.*image_in_core_size**3*8.)*Blockdata["no_of_processes_per_group"]/1.e9
	except:  current_mpi_settings_is_bad = 1
	if current_mpi_settings_is_bad == 1:ERROR("initial info is not provided", "check_mpi_settings", 1, Blockdata["myid"])
	try:
		mem_bytes = os.sysconf('SC_PAGE_SIZE')*os.sysconf('SC_PHYS_PAGES')# e.g. 4015976448
		mem_gib = mem_bytes/(1024.**3) # e.g. 3.74
		if( Blockdata["myid"] == Blockdata["main_node"]):print(line, "system mem info: %5.1f  G"%mem_gib)
	except:
		mem_gib = None
		if( Blockdata["myid"] == Blockdata["main_node"]):print(line, "It is not an unix machine!")
		else: pass
	if Tracker["constants"]["memory_per_node"] == -1.:
		if mem_gib: total_memory = mem_gib
		else:
			total_memory =  Blockdata["no_of_processes_per_group"]*2.0 # assume each CPU has 2.0 G
			if( Blockdata["myid"] == Blockdata["main_node"]):
				msg ="memory per node is not provided, sort3d assumes 2G per node"
				log.add(msg)
				print(line, msg)
		Tracker["constants"]["memory_per_node"] = total_memory
	else:
		msg ="memory per node: %f"%Tracker["constants"]["memory_per_node"]
		total_memory =  Tracker["constants"]["memory_per_node"]
		if( Blockdata["myid"] == Blockdata["main_node"]):
			log.add(msg)
			print(line, msg)
	if( Blockdata["myid"] == Blockdata["main_node"]):
		msg = "total number of particles:  %d  number of particles per group:  %d"%(Tracker["constants"]["total_stack"], Tracker["constants"]["img_per_grp"])
		log.add(msg)
		print(line, msg)
	if(Blockdata["myid"] == Blockdata["main_node"]):
		msg = "the available memory:  %5.1f  G"%total_memory
		log.add(msg)
		print(line, msg)
		msg = "total raw data:  %5.1f G raw data per node: %5.1f G"%(raw_data_size, raw_data_size_per_node)
		log.add(msg)
		print(line, msg)
	if (total_memory - sys_required_mem - raw_data_size_per_node - volume_size_per_node - sorting_data_size_per_node - 5.0) <0.0: 
		current_mpi_settings_is_bad = 1
		new_nproc =  raw_data_size*(2.*ratio**2+1.)*Blockdata["no_of_processes_per_group"]/(total_memory - 5. - sys_required_mem - volume_size_per_node)
		new_nproc =  int(new_nproc)
		if( Blockdata["myid"] == Blockdata["main_node"]):
			msg ="Suggestion: use  number of processes %d"%int(new_nproc)
			print(line, msg)
			log.add(msg)
		ERROR("In sufficient memory", "check_mpi_settings", 1, Blockdata["myid"])
	images_per_cpu = float(Tracker["constants"]["total_stack"])/float(Blockdata["nproc"])
	images_per_cpu_for_unaccounted_data  = Tracker["constants"]["img_per_grp"]*1.5/float(Blockdata["nproc"])
	if( Blockdata["myid"] == Blockdata["main_node"]):
		msg="current images per cpu:  %d "%int(images_per_cpu)
		log.add(msg)
		print(line, msg)
	if images_per_cpu < 5.0: ERROR("image per cpu less than 5", "check_mpi_settings", 1, Blockdata["myid"])	
	return
def get_sorting_image_size(original_data, partids, sparamstructure, snorm_per_particle, log):
	global Tracker, Blockdata
	from utilities    import wrap_mpi_bcast, read_text_file, write_text_file
	from applications import MPI_start_end
	iter = 0
	if(Blockdata["myid"] == Blockdata["main_node"]):
		msg = "start reconstruction with refinement window_size  %d"%Tracker["nxinit_refinement"]
		print(msg)
		log.add(msg)
		lpartids = read_text_file(partids, -1)
		if len(lpartids) == 1:
			iter_assignment = []
			for im in xrange(len(lpartids[0])):
				iter_assignment.append(randint(0,number_of_groups-1))# simple version
		else:
			iter_assignment = lpartids[0]
	else:   iter_assignment = 0
	iter_assignment = wrap_mpi_bcast(iter_assignment, Blockdata["main_node"])
	proc_list = [[None, None] for iproc in xrange(Blockdata["nproc"])]
	for iproc in xrange(Blockdata["nproc"]):
		iproc_image_start, iproc_image_end = MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], iproc)
		proc_list[iproc] = [iproc_image_start, iproc_image_end]
	compute_noise(Tracker["nxinit_refinement"])
	rdata = downsize_data_for_rec3D(original_data, Tracker["nxinit_refinement"], return_real = False, npad = 1)
	update_rdata_assignment(iter_assignment, proc_list, Blockdata["myid"], rdata)
	Tracker["nxinit"] = Tracker["nxinit_refinement"]
	compute_noise(Tracker["nxinit"])
	do3d_sorting_groups_fsc_only_iter(rdata, sparamstructure, snorm_per_particle, iteration = iter)
	del rdata
	if( Blockdata["myid"] == Blockdata["main_node"]):
		msg = "reconstruction with refinement window_size %d finshes"%Tracker["nxinit_refinement"]
		print(msg)
		log.add(msg)
	if( Blockdata["myid"] == Blockdata["main_node"]):
		fsc_data = []
		for igroup in xrange(Tracker["number_of_groups"]):
			for ichunk in xrange(2):
				tmp_fsc_data = read_text_file(os.path.join(Tracker["directory"], "fsc_driver_chunk%d_grp%03d_iter%03d.txt"%(ichunk, igroup, iter)), -1)
				fsc_data.append(tmp_fsc_data[0])
	else: fsc_data = 0
	fsc_data = wrap_mpi_bcast(fsc_data, Blockdata["main_node"])
	avg_fsc = [0.0 for i in xrange(len(fsc_data[0]))]
	avg_fsc[0] = 1.0
	for igroup in xrange(1): # Use group zero first
		for ifreq in xrange(1, len(fsc_data[0])):avg_fsc[ifreq] += fsc_data[igroup][ifreq]
	fsc143 = len(fsc_data[0])
	for igroup in xrange(Tracker["number_of_groups"]*2):
		for ifreq in xrange(1, len(fsc_data[igroup])):
			fsc_data[igroup][ifreq] = 2.*fsc_data[igroup][ifreq]/(1.+fsc_data[igroup][ifreq])
			if fsc_data[igroup][ifreq] < 0.143: break
		fsc143 = min(fsc143, ifreq)
	if fsc143 !=0: nxinit = int(fsc143)*2
	else: ERROR("program obtains wrong image size", "EQKmeans_by_dmatrix_orien_groups", 1, Blockdata["myid"])
	if(Blockdata["myid"] == Blockdata["main_node"]): write_text_file(avg_fsc, os.path.join(Tracker["directory"], "fsc_image_size.txt"))
	del iter_assignment
	del proc_list
	del fsc_data
	del avg_fsc
	return nxinit
	
def compute_noise(image_size):
	global Tracker, Blockdata
	from utilities    import get_im, model_blank
	from fundamentals import fft
	if Tracker["applybckgnoise"]: # from SPARX refinement only
		tsd = get_im(Tracker["bckgnoise"]) # invert power spectrum 
		nnx = tsd.get_xsize()
		nny = tsd.get_ysize()
		temp_image = model_blank(image_size, image_size)
		temp_image = fft(temp_image)
		nx = temp_image.get_xsize()
		ny = temp_image.get_ysize()
		Blockdata["bckgnoise"]  = []
		Blockdata["unrolldata"] = []
		for i in xrange(nny):
			prj = nnx*[0.0]
			for k in xrange(nnx):
				if tsd.get_value_at(k,i)>0.0: prj[k] = tsd.get_value_at(k,i)
			Blockdata["bckgnoise"].append(prj)
		for i in xrange(len(Blockdata["bckgnoise"])): Blockdata["unrolldata"].append(Util.unroll1dpw(ny, Blockdata["bckgnoise"][i]))
		if Blockdata["myid"] == Blockdata["main_node"]: print("size of bckgnoise, %d   %d"%(len(Blockdata["bckgnoise"]), len(Blockdata["bckgnoise"][0])))
	else: # from datastack and relion
		temp_image = model_blank(image_size, image_size)
		temp_image = fft(temp_image)
		nx = temp_image.get_xsize()
		ny = temp_image.get_ysize()
		Blockdata["bckgnoise"] = [1.0]*nx
		Blockdata["unrolldata"] = Util.unroll1dpw(ny, nx*[1.0])
	return

####------major procedure
def do_EQKmeans_nways_clustering_stable_seeds(workdir, initial_partids, params, sort_res, log_main):
	global Tracker, Blockdata
	from utilities import read_text_file, wrap_mpi_bcast, write_text_file
	import copy
	less_than_random_assignment = 0
	# Iteratively using EQKmeans to split a dataset into clusters of equal size till the number of 
	# the unaccounted is less than the minimum size of a cluster
	# input:  initial_partids
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] == Blockdata["main_node"]:
		msg = "----->>>do_EQKmeans_nways_clustering<<<------"
		print(line, msg)
		log_main.add(msg)
		Tracker["unaccounted_list"] = read_text_file(initial_partids, -1) # read all columns
		if len(Tracker["unaccounted_list"])>1: Tracker["unaccounted_list"] = Tracker["unaccounted_list"][1] # two column entries
		else: Tracker["unaccounted_list"] = Tracker["unaccounted_list"][0] # only one column
	else: Tracker["unaccounted_list"] = 0
	Tracker["unaccounted_list"]  = wrap_mpi_bcast(Tracker["unaccounted_list"], Blockdata["main_node"], MPI_COMM_WORLD)
	generation		             = 0
	Tracker["total_stack"]       = len(Tracker["unaccounted_list"])
	Tracker["number_of_groups"]  = get_number_of_groups(Tracker["total_stack"],Tracker["img_per_grp"])
	partids                      = initial_partids
	angle_step                   = Tracker["angle_step"]
	if Tracker["number_of_groups"]>1: # In case the number of the input particles is small
		while Tracker["number_of_groups"] >1:
			Tracker["min_orien_group_size"] = Tracker["number_of_groups"]*Tracker["minimum_ptl_number"]
			if Blockdata["myid"] == Blockdata["main_node"]: Tracker["partition_list"] = []
			else: Tracker["partition_list"] =  0
			
			if Blockdata["myid"] == Blockdata["main_node"]:
				Tracker["directory"]  = os.path.join(workdir, "generation%03d"%generation)
				cmd="{} {}".format("mkdir",Tracker["directory"])
				cmdexecute(cmd)
				cmd="{} {}".format("mkdir",os.path.join(Tracker["directory"], "tempdir"))
				cmdexecute(cmd)
				log_main.add("-------->>> generation       %5d"%generation)
				log_main.add("number of images per group:  %d"%Tracker["img_per_grp"])
				log_main.add("the initial number of groups:  %d  number of independent runs:  %d"%(Tracker["number_of_groups"], Tracker["constants"]["indep_runs"]))
			else:  Tracker["directory"] = 0
			Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], MPI_COMM_WORLD)
			create_nrandom_lists(partids)
			
			if Blockdata["myid"] == Blockdata["main_node"]: print("before get_orien_assignment_mpi")
			if Tracker["constants"]["symmetry"][0:1]=="c" or Tracker["constants"]["symmetry"][0:1]=="d":
				ptls_in_orien_groups =  get_orien_assignment_mpi(angle_step, partids, params, log_main)
			else: ptls_in_orien_groups = [range(Tracker["total_stack"])]
			if Blockdata["myid"] == Blockdata["main_node"]:print("after get_orien_assignment_mpi")
			Tracker["nxinit"] = Tracker["nxinit_refinement"]
			previous_params   = Tracker["previous_parstack"]
			original_data, norm_per_particle  = read_data_for_sorting(partids, params, previous_params)
			if Tracker["nosmearing"]:
				parameterstructure             = None
				Tracker["paramstructure_dict"] = None
				Tracker["paramstructure_dir"]  = None
			else: parameterstructure = read_paramstructure_for_sorting(partids, Tracker["paramstructure_dict"], Tracker["paramstructure_dir"])			
			Tracker["nxinit"] = get_sorting_image_size(original_data, os.path.join(Tracker["directory"], "independent_index_000.txt"), parameterstructure, norm_per_particle, log_main)
			if Blockdata["myid"] == Blockdata["main_node"]:
				line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
				msg ="determined image size for sorting is %d"%Tracker["nxinit"]
				print(line, msg)
				log_main.add(msg)
				tmpdir = os.path.join(Tracker["directory"], "tempdir")
				if os.path.exists(tmpdir):
					cmd="{} {}".format("rm -rf", tmpdir)
					cmdexecute(cmd)
			npremature = 0
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
					cmdexecute(cmd)
					cmd="{} {}".format("mkdir",os.path.join(Tracker["directory"], "tempdir"))
					cmdexecute(cmd)
				Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], MPI_COMM_WORLD)
				tmp_final_list, premature = EQKmeans_by_dmatrix_orien_groups(original_data, index_file, ptls_in_orien_groups, parameterstructure, norm_per_particle, clean_volumes= True)
				Tracker["directory"] =  os.path.join(workdir, "generation%03d"%generation)
				if Blockdata["myid"] == Blockdata["main_node"]: Tracker["partition_list"].append(Tracker["partition"])
				if premature ==1: npremature +=1
			Tracker["partition_list"] = wrap_mpi_bcast(Tracker["partition_list"], Blockdata["main_node"], MPI_COMM_WORLD)
			do_two_way_comparison_over_nindepruns(log_main)
			Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], MPI_COMM_WORLD)		
			if Blockdata["myid"] == Blockdata["main_node"]:
				Tracker["directory"] =  os.path.join(workdir, "generation%03d"%generation, "Reassign_from_EQKmeans")
				cmd="{} {}".format("mkdir",Tracker["directory"])
				cmdexecute(cmd)
				msg = "accounted ratio  %f  for  %d independent EQKmeans "%((len(Tracker["accounted_list"])/float(Tracker["total_stack"])*100.), Tracker["constants"]["indep_runs"])
				log_main.add(msg)
				print(msg)				
				if npremature !=0: msg = " premature runs %d  out of %d independent runs"%(npremature, Tracker["constants"]["indep_runs"])
				else: msg =" %d EQKmeans runs all end up with mature stop"%Tracker["constants"]["indep_runs"]
				log_main.add(msg)
				print(msg)
				if len(Tracker["accounted_list"])/float(Tracker["total_stack"]) < 1./float(Tracker["number_of_groups"]): less_than_random_assignment =1
			less_than_random_assignment = bcast_number_to_all(less_than_random_assignment, Blockdata["main_node"], MPI_COMM_WORLD)
			if less_than_random_assignment ==1: ERROR("accounted ratio less than accounted ratio of random assignment ", "sxsort3d.py", 1, Blockdata["myid"])
			tmp_list   = []
			npremature = 0
			for riter in xrange(2): # Set it to two for the time being
				assign_list = resize_groups_from_stable_members_mpi(Tracker["Accounted_on_disk"], Tracker["Unaccounted_on_disk"])
				Tracker["directory"]  = os.path.join(workdir,"generation%03d"%generation,"Reassign_from_EQKmeans", "EQKmeans_random%d"%riter)
				line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
				if Blockdata["myid"] == Blockdata["main_node"]:
					msg ="------>>>Reassign unaccounted ones to accounted in groups %d <<<-------"%riter
					print(line, msg)
					log_main.add(msg)
					msg = "current number of groups:  %d   total data  %d  "%(Tracker["number_of_groups"], Tracker["total_stack"]) 
					print(line, msg)
					log_main.add(msg)
					if npremature>0: msg = " premature runs %d out of  two runs"%npremature
					else: msg = "mature stop"
					print(line, msg)
					log_main.add(msg)
					cmd="{} {}".format("mkdir",Tracker["directory"])
					cmdexecute(cmd)
					write_text_file(assign_list[0], os.path.join(Tracker["directory"], "random_list.txt"))
				index_file     = os.path.join(Tracker["directory"], "random_list.txt")
				tmp_final_list, premature = EQKmeans_by_dmatrix_orien_groups(original_data, index_file, ptls_in_orien_groups, parameterstructure, norm_per_particle, clean_volumes= True)
				if premature ==1: npremature +=1
				newlist = [[],[]]
				for im in xrange(len(tmp_final_list)):
					newlist[0].append(tmp_final_list[im][0])
					newlist[1].append(tmp_final_list[im][1])
				tmp_list.append(newlist)
				if Blockdata["myid"] == Blockdata["main_node"]: write_text_row(tmp_final_list, os.path.join(workdir, "generation%03d"%generation, "Reassign_from_EQKmeans", "res_list%d.txt"%riter))
			Tracker["directory"] = os.path.join(workdir,"generation%03d"%generation, "Reassign_from_EQKmeans")
			ratio, matched_pairs, clusters = compare_two_iterations(tmp_list[0][0], tmp_list[1][0], Tracker["number_of_groups"])
			if Blockdata["myid"] == Blockdata["main_node"]:
				line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
				final_list = [] 
				pall  = read_text_file(os.path.join(workdir, "generation%03d"%generation, "independent_index_000.txt"), -1)
				for ic in xrange(len(clusters)):
					any = clusters[ic]
					any.tolist()
					orig_any = []
					for a in any: orig_any.append(pall[1][a])
					final_list.append(orig_any)
					write_text_file(orig_any, os.path.join(Tracker["directory"], "ncluster%d.txt")%ic)
				plist, qlist, ulist = clusters_to_plist(clusters, pall[1])
				sort_res.append(qlist)
				msg = "the reproducible ratio after re-assignment  %f"%(ratio*100.)
				print(line, msg)
				log_main.add(msg)
				write_text_file( plist, os.path.join(Tracker["directory"], "plist.txt"))
				write_text_file( ulist, os.path.join(Tracker["directory"], "Unaccounted.txt"))
			else: 
				final_list = 0
				qlist = 0
			final_list = wrap_mpi_bcast(final_list, Blockdata["main_node"], MPI_COMM_WORLD)
			qlist = wrap_mpi_bcast(qlist, Blockdata["main_node"], MPI_COMM_WORLD)		
			Tracker["Unaccounted_on_disk"] = os.path.join(Tracker["directory"], "Unaccounted.txt") # the input for the next round iteration
			# resample to particle original data size from stable members
			partids = Tracker["Unaccounted_on_disk"]
			if Blockdata["myid"] == Blockdata["main_node"]:
				line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"		
				msg  = "Smallest group size %d "%Tracker["constants"]["minimum_grp_size"]
				log_main.add(msg)
				print(line, msg)
				msg = "the final number of groups: %d"%Tracker["number_of_groups"] 
				Tracker["total_stack"]       = len(read_text_row(partids))
				Tracker["number_of_groups"]  = get_number_of_groups(Tracker["total_stack"], Tracker["img_per_grp"])
				msg ="Intermediate results are reconstructed in directory %s"%Tracker["directory"]
				print(msg)
				log_main.add(msg)
				tmpdir = os.path.join(Tracker["directory"], "tempdir")
				if os.path.exists(tmpdir): 
					cmd="{} {}".format("rm -rf", tmpdir)
					cmdexecute(cmd)
			Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], MPI_COMM_WORLD)
			generation +=1
			del original_data
			del parameterstructure
			del norm_per_particle
	else:
		if Blockdata["myid"] == Blockdata["main_node"]:
			msg ="The total number of particles is less than number_of_particles_per_group"
			print(line, msg)
			log_main.add(msg)
			if os.path.exsists(os.path.join(Tracker["directory"], "tempdir")): cmd = "{} {}".format("rm -rf", os.path.join(Tracker["directory"], "tempdir"))
		else:  final_list = 0
		final_list  = wrap_mpi_bcast(final_list, Blockdata["main_node"], MPI_COMM_WORLD)
		del original_data
		del parameterstructure
		del norm_per_particle
	return final_list # the result of the last iteration
### EQKmeans
def EQKmeans_by_dmatrix_orien_groups(original_data, partids, ptls_in_orien_groups, paramstructure, norm_per_particle, clean_volumes = False):
	global Tracker, Blockdata
	#<<<<---------- >>>>EQKmeans starts<<<<------------ 
	log	                    = Logger()
	log                     = Logger(BaseLogger_Files())
	log.prefix              = Tracker["directory"]+"/"
	total_stack             = Tracker["total_stack"]
	premature               = 0
	changed_nptls           = 100.0
	number_of_groups        = Tracker["number_of_groups"]
	nima                    = len(original_data)
	image_start, image_end  = MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"])
	iter_max                = Tracker["total_number_of_iterations"]//3
	stopercnt               = 3.
	total_iter              = 0
	require_check_setting   = False
	partial_rec3d           = False
	if Tracker["mask3D"]:
		mask3D = get_im(Tracker["mask3D"])
		if mask3D.get_xsize() != Tracker["nxinit"]: mask3D = fdecimate(mask3D, Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"], True, False)
	else: 
		mask3D = model_circle(Tracker["constants"]["radius"], Tracker["constants"]["nnxo"], Tracker["constants"]["nnxo"], Tracker["constants"]["nnxo"])
		mask3D = fdecimate(mask3D, Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"], True, False)
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if( Blockdata["myid"] == Blockdata["main_node"]):
		msg = "------>>>>EQKmeans_by_dmatrix_orien_groups starts<<<<--------"
		log.add(msg)
		msg = "total_stack:  %d"%Tracker["total_stack"]
		log.add(msg)
		print(line, msg)
		msg = "number_of_groups:  %d"%Tracker["number_of_groups"]
		log.add(msg)
		print(line, msg)
		msg = "ctf:  %s"%Tracker["constants"]["CTF"]
		log.add(msg)
		print(line, msg)
		msg = "Currrent directory:  %s"%Tracker["directory"]
		log.add(msg)
		print(line, msg)
		msg = "nxinit:   %d"%Tracker["nxinit"]
		log.add(msg)
		print(line, msg)
		msg = "Symmetry group:  %s"%Tracker["constants"]["symmetry"] 
		log.add(msg)
		print(line, msg)
		msg = "The total number of iterations:  %d"%iter_max
		log.add(msg)
		print(line, msg)
		msg = "The stop criterion:  %f"%stopercnt
		log.add(msg)
		print(line, msg)
		msg = "Number of orien groups:  %d"%len(ptls_in_orien_groups)
		log.add(msg)
		print(line, msg)
		msg = "3-D mask:  %s"%Tracker["mask3D"]
		log.add(msg)
		print(line, msg)
		msg = "focus:  %s"%Tracker["focus3D"]
		log.add(msg)
		print(line, msg)
		msg = "Rec3d interpolation: %s"%Tracker["constants"]["interpolation"]
		log.add(msg)
		print(line, msg)
		msg = "Comparison method: %s"%Tracker["constants"]["comparison_method"]
		log.add(msg)
		print(line, msg)
		lpartids = read_text_file(partids, -1)
		if len(lpartids) == 1:
			iter_assignment = []
			for im in xrange(len(lpartids[0])):iter_assignment.append(randint(0,number_of_groups-1))# simple version
		else: iter_assignment = lpartids[0]
	else:   iter_assignment = 0
	iter_assignment = wrap_mpi_bcast(iter_assignment, Blockdata["main_node"])
	proc_list = [[None, None] for iproc in xrange(Blockdata["nproc"])]
	for iproc in xrange(Blockdata["nproc"]):
		iproc_image_start, iproc_image_end = MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], iproc)
		proc_list[iproc] = [iproc_image_start, iproc_image_end]
	compute_noise(Tracker["nxinit"])
	cdata, rdata = downsize_data_for_sorting(original_data, preshift = True, npad = 1)# pay attentions to shifts!
	srdata = precalculate_shifted_data_for_recons3D(rdata, paramstructure, Tracker["refang"], \
	   Tracker["rshifts"], Tracker["delta"], Tracker["avgnorm"], Tracker["nxinit"], Tracker["constants"]["nnxo"], Tracker["nosmearing"], norm_per_particle)
	del rdata
	while total_iter<Tracker["total_number_of_iterations"]:
		iter = 0
		while iter <iter_max:
			if(Blockdata["myid"] == Blockdata["main_node"]):
				msg = "Iteration   %d"%total_iter
				log.add(msg)
				msg = "particles changes assignment  %f, and current stop criterion is  %f, current image size %d"%(changed_nptls, stopercnt, Tracker["nxinit"])
				log.add(msg)
				msg = "assignment length  %d"%len(iter_assignment)
				log.add(msg)
				write_text_file(iter_assignment, os.path.join(Tracker["directory"], "assignment%03d.txt"%total_iter))
			if changed_nptls<50.0: partial_rec3d = True
			update_data_assignment(cdata, srdata, iter_assignment, proc_list, Tracker["nosmearing"], Blockdata["myid"])
			mpi_barrier(MPI_COMM_WORLD)
			do3d_sorting_groups_nofsc_smearing_iter(srdata, partial_rec3d, iteration = total_iter)
			mpi_barrier(MPI_COMM_WORLD)
			local_peaks = [0.0 for im in xrange(number_of_groups*nima)]
			total_im = 0
			for iref in xrange(number_of_groups):
				if(Blockdata["myid"] == Blockdata["main_node"]):
					try: fsc143 = Tracker["fsc143"][iref]
					except:	fsc143 = 0.0
					try: fsc05  = Tracker["fsc05"][iref]
					except:	fsc05  = 0.0			
					ref_vol = get_im(os.path.join(Tracker["directory"],"vol_grp%03d_iter%03d.hdf"%(iref, total_iter)))
					nnn = ref_vol.get_xsize()
					if(Tracker["nxinit"] != nnn): ref_vol = fdecimate(ref_vol, Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"], True, False)
					stat = Util.infomask(ref_vol, mask3D, False)
					ref_vol -= stat[0]
					if stat[1]!=0.0:Util.mul_scalar(ref_vol, 1.0/stat[1])
					if Tracker["mask3D"]:ref_vol *=mask3D
				else: ref_vol = model_blank(Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"])
				bcast_EMData_to_all(ref_vol, Blockdata["myid"], Blockdata["main_node"])		
				if Tracker["constants"]["comparison_method"] =="cross": ref_peaks = compare_two_images_cross(cdata, ref_vol)
				else: ref_peaks = compare_two_images_eucd(cdata, ref_vol)
				for im in xrange(nima):
					local_peaks[total_im] = ref_peaks[im]
					total_im +=1
				mpi_barrier(MPI_COMM_WORLD)
			del ref_vol
			# pass to main_node
			if Blockdata["myid"] == Blockdata["main_node"]:
				peaks =[[ 0.0 for im in xrange(Tracker["total_stack"])] for iref in xrange(number_of_groups)]
				for im in xrange(len(local_peaks)): peaks[im/nima][im%nima+image_start] = local_peaks[im]
			else: peaks = 0		
			if Blockdata["myid"] != Blockdata["main_node"]: wrap_mpi_send(local_peaks, Blockdata["main_node"], MPI_COMM_WORLD)
			else:
				for iproc in xrange(Blockdata["nproc"]):
					if iproc != Blockdata["main_node"]:
						local_peaks = wrap_mpi_recv(iproc, MPI_COMM_WORLD)
						iproc_nima = proc_list[iproc][1]-proc_list[iproc][0]
						for im in xrange(len(local_peaks)): peaks[im/iproc_nima][im%iproc_nima+proc_list[iproc][0]] = local_peaks[im]
			peaks = wrap_mpi_bcast(peaks, Blockdata["main_node"], MPI_COMM_WORLD)
			if Blockdata["myid"] == Blockdata["main_node"]:
				ss = 0
				for io in xrange(len(ptls_in_orien_groups)): ss  +=len(ptls_in_orien_groups[io])
				for ip in xrange(len(peaks)): print(ip, len(peaks[ip]))
			mpi_barrier(MPI_COMM_WORLD)
			last_iter_assignment = copy.copy(iter_assignment)
			iter_assignment = [-1 for iptl in xrange(Tracker["total_stack"])]
			for iorien in xrange(len(ptls_in_orien_groups)):
				if iorien%Blockdata["nproc"] == Blockdata["myid"]:
					local_assignment = do_assignment_by_dmatrix_orien_group(peaks, ptls_in_orien_groups[iorien], number_of_groups)
					for iptl in xrange(len(ptls_in_orien_groups[iorien])):
						iter_assignment[ptls_in_orien_groups[iorien][iptl]] = local_assignment[iptl]
			mpi_barrier(MPI_COMM_WORLD)		
			if Blockdata["myid"] != Blockdata["main_node"]: wrap_mpi_send(iter_assignment, Blockdata["main_node"], MPI_COMM_WORLD)
			else:
				for iproc in xrange(Blockdata["nproc"]):
					if iproc != Blockdata["main_node"]:
						dummy = wrap_mpi_recv(iproc, MPI_COMM_WORLD)
						for iptl in xrange(len(dummy)):
							if dummy[iptl] !=-1:iter_assignment[iptl] = dummy[iptl]
							else: pass							
			mpi_barrier(MPI_COMM_WORLD)
			iter_assignment = wrap_mpi_bcast(iter_assignment, Blockdata["main_node"], MPI_COMM_WORLD)
			ratio, newindices, stable_clusters = compare_two_iterations(iter_assignment, last_iter_assignment, number_of_groups)
			changed_nptls = 100.-ratio*100.
			if changed_nptls < stopercnt and total_iter <=10:stopercnt = changed_nptls/2. # reduce stop criterion to gain improvement in clustering
			stopercnt = max(stopercnt, Tracker["constants"]["stop_eqkmeans_percentage"]) # But not exceed the specified number
			iter +=1
			total_iter +=1
			if Blockdata["myid"] == Blockdata["main_node"]:
				msg ="compared to the last iteration"
				log.add(msg)
				for igrp in xrange(len(newindices)):
					msg= "iter  %d  group  %d  matches iter %d  group %d "%(total_iter, \
					newindices[igrp][0], total_iter-1, newindices[igrp][1])
					log.add(msg)
			mpi_barrier(MPI_COMM_WORLD)
			if changed_nptls < stopercnt: break
		if changed_nptls<stopercnt: break
	#Finalize
	update_data_assignment(cdata, srdata, iter_assignment, proc_list, Tracker["nosmearing"], Blockdata["myid"])
	res_sort3d = get_sorting_all_params(cdata)
	del cdata
	del srdata
	del iter_assignment
	del peaks
	del last_iter_assignment
	if mask3D: del mask3D
	if changed_nptls > 15.0: require_check_setting = True
	if(Blockdata["myid"] == Blockdata["main_node"]):
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		if changed_nptls > Tracker["constants"]["stop_eqkmeans_percentage"]: 
			msg ="EQKmeans premature stop with changed particles ratio %f and image size %d"%(changed_nptls,Tracker["nxinit"])
			premature  = 1
		else: msg = "EQKmeans mature stop with changed particles ratio %f within %d iterations and actually used stop percentage is %f"%(\
		        changed_nptls, total_iter, stopercnt)
		print(line, msg)
		log.add(msg)
		Tracker["partition"], ali3d_params_list = parsing_sorting_params(partids, res_sort3d)
		write_text_row(Tracker["partition"], os.path.join(Tracker["directory"],"list.txt"))
		cmd = "{} {}".format("rm -rf", os.path.join(Tracker["directory"], "tempdir"))
		if os.path.exists(os.path.join(Tracker["directory"], "tempdir")): cmdexecute(cmd)
		if clean_volumes:
			cmd = "{} {}".format("rm -rf", os.path.join(Tracker["directory"], "vol_*.hdf"))
			cmdexecute(cmd)
	else:Tracker["partition"] = 0
	Tracker["partition"] = wrap_mpi_bcast(Tracker["partition"], Blockdata["main_node"])
	premature  = wrap_mpi_bcast(premature, Blockdata["main_node"])
	if require_check_setting:
		if(Blockdata["myid"] == Blockdata["main_node"]): print("Too large changed particles, and the sorting settings, such as img_per_grp requires a check")
		#from mpi import mpi_finalize
		#mpi_finalize()
		#exit()
	return Tracker["partition"], premature
### various reading data
### 1
def get_shrink_data_sorting(partids, partstack, return_real = False, preshift = True, apply_mask = True, npad = 1):
	# The function will read from stack a subset of images specified in partids
	#   and assign to them parameters from partstack with optional CTF application and shifting of the data.
	# So, the lengths of partids and partstack are the same.
	#  The read data is properly distributed among MPI threads.
	# 10142015 --- preshift is set to True when doing 3-D sorting.
	# chunk_id are set when data is read in
	global Tracker, Blockdata
	from utilities      import wrap_mpi_bcast, read_text_row
	from fundamentals	import resample, fshift
	from filter			import filt_ctf
	from applications	import MPI_start_end
	from EMAN2          import Region
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if( Blockdata["myid"] == Blockdata["main_node"]): print(line,"get_shrink_data_sorting")
	mask2D		= model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	shrinkage 	= Tracker["nxinit"]/float(Tracker["constants"]["nnxo"])
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
	if(Blockdata["myid"] == Blockdata["main_node"]):  partstack = read_text_row(partstack)
	else:  partstack = 0
	partstack = wrap_mpi_bcast(partstack, Blockdata["main_node"])
	
	if(Tracker["total_stack"] < Blockdata["nproc"]): ERROR("Wrong MPI settings!", "get_shrink_data_sorting", 1, Blockdata["myid"])
	else: image_start, image_end = MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"])
	lpartids  = lpartids[image_start:image_end]
	groupids  = groupids[image_start:image_end]
	nima      =	image_end - image_start
	data      = [None]*nima
	for im in xrange(nima):
		data[im] = get_im(Tracker["constants"]["orgstack"], lpartids[im])	
		try: phi, theta, psi, sx, sy, chunk_id, particle_group_id = partstack[lpartids[im]][0], partstack[lpartids[im]][1],\
		 partstack[lpartids[im]][2], partstack[lpartids[im]][3], partstack[lpartids[im]][4], partstack[lpartids[im]][5], partstack[lpartids[im]][6]
		except: phi, theta, psi, sx, sy, chunk_id, particle_group_id = partstack[lpartids[im]][0], partstack[lpartids[im]][1],\
		 partstack[lpartids[im]][2], partstack[lpartids[im]][3], partstack[lpartids[im]][4], partstack[lpartids[im]][5], -1		 
		if preshift:# always true
			data[im]  = fshift(data[im],sx,sy)
			sx = 0.0
			sy = 0.0
		st = Util.infomask(data[im], mask2D, False)
		data[im] -= st[0]
		data[im] /= st[1]
		if apply_mask: data[im] = cosinemask(data[im],radius = Tracker["constants"]["radius"])
		# FT
		data[im] = fft(data[im])
		nny =  data[im].get_ysize()
		if Tracker["constants"]["CTF"]:
			ctf_params = data[im].get_attr("ctf")
			data[im]   = fdecimate(data[im], Tracker["nxinit"]*npad, Tracker["nxinit"]*npad, 1, False, False)
			ctf_params.apix = ctf_params.apix/shrinkage
			data[im].set_attr('ctf', ctf_params)
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
		data[im].set_attr("particle_group", particle_group_id)
		if Tracker["applybckgnoise"]:
			data[im].set_attr("bckgnoise", Blockdata["bckgnoise"][particle_group_id])
			data[im].set_attr("qt", float(Tracker["constants"]["nnxo"]*Tracker["constants"]["nnxo"]))
		else: data[im].set_attr("bckgnoise", Blockdata["bckgnoise"]) # constant list
	return data
###2
def get_shrink_data_sorting_smearing(partids, partstack, return_real = False, preshift = True, apply_mask = True, npad = 1):
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
	if( Blockdata["myid"] == Blockdata["main_node"]): print(line,"get_shrink_data_sorting")
	mask2D		= model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	shrinkage 	= Tracker["nxinit"]/float(Tracker["constants"]["nnxo"])
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
	if(Blockdata["myid"] == Blockdata["main_node"]): partstack = read_text_row(partstack)
	else:  partstack = 0
	partstack = wrap_mpi_bcast(partstack, Blockdata["main_node"])
	if(Tracker["total_stack"] < Blockdata["nproc"]): ERROR("Wrong MPI settings!", "get_shrink_data_sorting", 1, Blockdata["myid"])
	else:   image_start, image_end = MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"])
	lpartids  = lpartids[image_start:image_end]
	groupids  = groupids[image_start:image_end]
	nima      =	image_end - image_start
	data      = [None]*nima
	norm_per_particle = []
	for im in xrange(nima):
		data[im] = get_im(Tracker["constants"]["orgstack"], lpartids[im])	
		try: phi, theta, psi, sx, sy, chunk_id, particle_group_id, norm = partstack[lpartids[im]][0], partstack[lpartids[im]][1],\
		 partstack[lpartids[im]][2], partstack[lpartids[im]][3], partstack[lpartids[im]][4], partstack[lpartids[im]][5], partstack[lpartids[im]][6], partstack[lpartids[im]][7]
		except: phi, theta, psi, sx, sy, chunk_id, particle_group_id, norm = partstack[lpartids[im]][0], partstack[lpartids[im]][1],\
		 partstack[lpartids[im]][2], partstack[lpartids[im]][3], partstack[lpartids[im]][4], partstack[lpartids[im]][5], -1, 1
		if preshift:# always true
			data[im]  = fshift(data[im],sx,sy)
			sx = 0.0
			sy = 0.0
		st = Util.infomask(data[im], mask2D, False)
		data[im] -= st[0]
		data[im] /= st[1]	
		if apply_mask: data[im] = cosinemask(data[im],radius = Tracker["constants"]["radius"])
		# FT
		data[im] = fft(data[im])
		nny =  data[im].get_ysize()
		if Tracker["constants"]["CTF"] :
			ctf_params = data[im].get_attr("ctf")
			data[im]   = fdecimate(data[im], Tracker["nxinit"]*npad, Tracker["nxinit"]*npad, 1, False, False)
			ctf_params.apix = ctf_params.apix/shrinkage
			data[im].set_attr('ctf', ctf_params)
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
		data[im].set_attr("particle_group", particle_group_id)
		if Tracker["applybckgnoise"]:
			data[im].set_attr("bckgnoise", Blockdata["bckgnoise"][particle_group_id])
			data[im].set_attr("qt", float(Tracker["constants"]["nnxo"]*Tracker["constants"]["nnxo"]))
		else: data[im].set_attr("bckgnoise", Blockdata["bckgnoise"]) # constant list
		norm_per_particle.append(norm)
	return data, norm_per_particle
###3
def get_data_prep_compare_rec3d(partids, partstack, return_real = False, preshift = True, npad = 1):
	# The function will read from stack a subset of images specified in partids
	#   and assign to them parameters from partstack with optional CTF application and shifting of the data.
	# So, the lengths of partids and partstack are the same.
	
	global Tracker, Blockdata
	from fundamentals	import resample, fshift, fft
	from filter			import filt_ctf
	from applications	import MPI_start_end
	from EMAN2          import Region
	from utilities      import model_circle, wrap_mpi_bcast, get_im, model_blank, set_params_proj
	# functions:
	# read in data
	# apply mask, and prepare focus projection if focus3D is specified
	# return  1. cdata: data for image comparison, always in Fourier format
	#         2. rdata: data for reconstruction, 4nn return real image

	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if( Blockdata["myid"] == Blockdata["main_node"]): print(line,"read_data in ")		
	mask2D	  = model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	shrinkage = Tracker["nxinit"]/float(Tracker["constants"]["nnxo"])
	radius    = int(Tracker["constants"]["radius"] * shrinkage +0.5)
	if Tracker["applybckgnoise"]:
		oneover = []
		nnx = len(Blockdata["bckgnoise"][0])
		for i in xrange(len(Blockdata["bckgnoise"])):
			temp = [0.0]*nnx
			for k in xrange(nnx):
				if(Blockdata["bckgnoise"][i][k] > 0.0):  temp[k] = 1.0/sqrt(Blockdata["bckgnoise"][i][k])
			oneover.append(temp)
		del temp
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
	lpartids = wrap_mpi_bcast(lpartids, Blockdata["main_node"])
	groupids = wrap_mpi_bcast(groupids, Blockdata["main_node"])
	Tracker["total_stack"] = len(lpartids)
	if(Blockdata["myid"] == Blockdata["main_node"]):  partstack = read_text_row(partstack)
	else:  partstack = 0
	partstack = wrap_mpi_bcast(partstack, Blockdata["main_node"])
	if(Tracker["total_stack"] < Blockdata["nproc"]):
		ERROR("number of processors in use is larger than the total number of particles", \
		  "get_data_and_prep", 1, Blockdata["myid"])
	else: image_start, image_end = MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"])
	lpartids  = lpartids[image_start:image_end]
	groupids  = groupids[image_start:image_end]
	if Tracker["focus3D"]: # focus mask is applied
		if Blockdata["myid"] == Blockdata["main_node"]:
			focus3d     = get_im(Tracker["focus3D"])
			focus3d_nx  = focus3d.get_xsize()
			if focus3d_nx != Tracker["constants"]["nnxo"]: # So the decimated focus volume can be directly used
				focus3d = resample(focus3d, float(Tracker["constants"]["nnxo"])/float(focus3d_nx))
		else: focus3d = model_blank(Tracker["constants"]["nnxo"], Tracker["constants"]["nnxo"], Tracker["constants"]["nnxo"])
		bcast_EMData_to_all(focus3d, Blockdata["myid"], Blockdata["main_node"])
		focus3d = prep_vol(focus3d, 1, 1)
	#  Preprocess the data
	#  mask2D    =	model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	nima  = image_end - image_start
	cdata = [None]*nima
	rdata = [None]*nima
	for im in xrange(nima):
		image = get_im(Tracker["constants"]["orgstack"], lpartids[im])
		try: phi, theta, psi, sx, sy, chunk_id, particle_group_id  = partstack[lpartids[im]][0], partstack[lpartids[im]][1], partstack[lpartids[im]][2], \
			partstack[lpartids[im]][3], partstack[lpartids[im]][4], partstack[lpartids[im]][5], partstack[lpartids[im]][6]
		except: phi, theta, psi, sx, sy, chunk_id, particle_group_id  = partstack[lpartids[im]][0], partstack[lpartids[im]][1], partstack[lpartids[im]][2], \
		  partstack[lpartids[im]][3], partstack[lpartids[im]][4], partstack[lpartids[im]][5], -1 	  
		if preshift:# always true
			image = fshift(image,sx,sy)
			sx = 0.0
			sy = 0.0
		st = Util.infomask(image, mask2D, False)
		image -= st[0]
		image /= st[1]
		cimage = image.copy()
		if Tracker["applybckgnoise"]:
			if Tracker["applymask"]:
				if Tracker["constants"]["hardmask"]: cimage = cosinemask(cimage, radius = Tracker["constants"]["radius"])
				else:
					bckg = model_gauss_noise(1.0,Tracker["constants"]["nnxo"]+2,Tracker["constants"]["nnxo"])
					bckg.set_attr("is_complex",1)
					bckg.set_attr("is_fftpad",1)
					bckg = fft(filt_table(bckg, oneover[particle_group_id]))
					#  Normalize bckg noise in real space, only region actually used.
					st = Util.infomask(bckg, mask2D, False)
					bckg -= st[0]
					bckg /= st[1]
					cimage = cosinemask(cimage,radius = Tracker["constants"]["radius"], bckg = bckg)
		else:
			if Tracker["applymask"]: cimage  = cosinemask(cimage, radius = Tracker["constants"]["radius"])
		# FT
		image = fft(image)
		cimage = fft(cimage)		
		if Tracker["constants"]["CTF"] :
			ctf_params = image.get_attr("ctf")
			image = fdecimate(image, Tracker["nxinit"]*npad, Tracker["nxinit"]*npad, 1, False, False)
			cimage = fdecimate(cimage, Tracker["nxinit"]*npad, Tracker["nxinit"]*npad, 1, False, False)
			ctf_params.apix = ctf_params.apix/shrinkage
			image.set_attr('ctf', ctf_params)
			cimage.set_attr('ctf', ctf_params)
			image.set_attr('ctf_applied', 0)
			cimage.set_attr('ctf_applied', 0)
			if return_real:image = fft(image)
		else:
			ctf_params = image.get_attr_default("ctf", False)
			if  ctf_params:
				ctf_params.apix = ctf_params.apix/shrinkage
				image.set_attr('ctf', ctf_params)
				image.set_attr('ctf_applied', 0)
				cimage.set_attr('ctf', ctf_params)
				cimage.set_attr('ctf_applied', 0)
			image = fdecimate(image, nxinit*npad, nxinit*npad, 1, True, False)
			cimage = fdecimate(cimage, nxinit*npad, nxinit*npad, 1, True, False)
			apix = Tracker["constants"]["pixel_size"]
			image.set_attr('apix', apix/shrinkage)
			cimage.set_attr('apix', apix/shrinkage)
		cimage.set_attr("padffted",1)
		cimage.set_attr("npad", npad)
		if not return_real:
			image.set_attr("padffted",1)
			image.set_attr("npad", npad)
		set_params_proj(image,[phi, theta, psi, 0.0, 0.0])
		image.set_attr("chunk_id", chunk_id)
		image.set_attr("group", groupids[im])
		image.set_attr("particle_group", particle_group_id)		
		set_params_proj(cimage,[phi, theta, psi, 0.0, 0.0])
		cimage.set_attr("chunk_id", chunk_id)
		cimage.set_attr("group", groupids[im])
		cimage.set_attr("particle_group", particle_group_id)
		rdata[im] =  image
		cdata[im] =  cimage
		if Tracker["applybckgnoise"]: 
			rdata[im].set_attr("bckgnoise", Blockdata["bckgnoise"][particle_group_id])
			if Tracker["constants"]["comparison_method"] == "cross": Util.mulclreal(cdata[im], Blockdata["unrolldata"][particle_group_id])                                
		if Tracker["focus3D"]:
			cdata[im] = fft(binarize(prgl(focus3d, [phi, theta, psi, 0.0, 0.0], 1, True), 1)*fft(cdata[im]))
			if Tracker["constants"]["CTF"]: cdata[im].set_attr("ctf", rdata[im].get_attr("ctf"))
		cdata[im].set_attr("is_complex",0)	
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if( Blockdata["myid"] == Blockdata["main_node"]): print(line,"reading data finishes")
	return cdata, rdata
#####4
def get_shrink_data_final(nxinit, procid, original_data = None, oldparams = None, \
		return_real = False, preshift = False, apply_mask = True, nonorm = False, npad = 1):
	global Tracker, Blockdata
	"""
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
	#from fundamentals import resample
	from utilities    import get_im, model_gauss_noise, set_params_proj, get_params_proj
	from fundamentals import fdecimate, fshift, fft
	from filter       import filt_ctf, filt_table
	from applications import MPI_start_end
	from math         import sqrt
	
	if( Blockdata["myid"] == Blockdata["main_node"] ):
		print( "  " )
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print(  line, "Processing data  onx: %3d, nx: %3d, CTF: %s, applymask: %s, preshift: %s."%(Tracker["constants"]["nnxo"], nxinit, Tracker["constants"]["CTF"], apply_mask, preshift) )
	#  Preprocess the data
	mask2D  	= model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	nima 		= len(original_data)
	shrinkage 	= nxinit/float(Tracker["constants"]["nnxo"])
	#  Note these are in Fortran notation for polar searches
	#txm = float(nxinit-(nxinit//2+1) - radius -1)
	#txl = float(2 + radius - nxinit//2+1)
	radius 	= int(Tracker["constants"]["radius"]*shrinkage + 0.5)
	txm    	= float(nxinit-(nxinit//2+1) - radius)
	txl    	= float(radius - nxinit//2+1)

	if Blockdata["bckgnoise"] :
		oneover = []
		nnx = Blockdata["bckgnoise"][0].get_xsize()
		for i in xrange(len(Blockdata["bckgnoise"])):
			temp = [0.0]*nnx
			for k in xrange(nnx):
				if( Blockdata["bckgnoise"][i].get_value_at(k) > 0.0):  temp[k] = 1.0/sqrt(Blockdata["bckgnoise"][i].get_value_at(k))
			oneover.append(temp)
		del temp
	Blockdata["accumulatepw"][procid] = [None]*nima
	data = [None]*nima
	for im in xrange(nima):

		if Tracker["mainiteration"] ==1:
			phi, theta, psi, sx, sy, = oldparams[im][0], oldparams[im][1], oldparams[im][2], oldparams[im][3], oldparams[im][4]
			wnorm  = 1.0 
		else: phi, theta, psi, sx, sy, wnorm = oldparams[im][0], oldparams[im][1], oldparams[im][2], oldparams[im][3], oldparams[im][4], oldparams[im][7]
		if preshift:
			sx = int(round(sx))
			sy = int(round(sy))
			data[im]  = cyclic_shift(original_data[im],sx,sy)
			#  Put rounded shifts on the list, note image has the original floats - check whether it may cause problems
			oldparams[im][3] = sx
			oldparams[im][4] = sy
			sx = 0.0
			sy = 0.0
		else:  data[im] = original_data[im].copy()
		st = Util.infomask(data[im], mask2D, False)
		data[im] -= st[0]
		data[im] /= st[1]
		if data[im].get_attr_default("bckgnoise", None) :  data[im].delete_attr("bckgnoise")
		#  Do bckgnoise if exists
		if Blockdata["bckgnoise"]:
			if apply_mask:
				if Tracker["constants"]["hardmask"]:
					data[im] = cosinemask(data[im],radius = Tracker["constants"]["radius"])
				else:
					bckg = model_gauss_noise(1.0,Tracker["constants"]["nnxo"]+2,Tracker["constants"]["nnxo"])
					bckg.set_attr("is_complex",1)
					bckg.set_attr("is_fftpad",1)
					bckg = fft(filt_table(bckg, oneover[data[im].get_attr("particle_group")]))
					#  Normalize bckg noise in real space, only region actually used.
					st = Util.infomask(bckg, mask2D, False)
					bckg -= st[0]
					bckg /= st[1]
					data[im] = cosinemask(data[im],radius = Tracker["constants"]["radius"], bckg = bckg)
		else:
			#  if no bckgnoise, do simple masking instead
			if apply_mask:  data[im] = cosinemask(data[im],radius = Tracker["constants"]["radius"] )
		#  resample will properly adjusts shifts and pixel size in ctf
		#data[im] = resample(data[im], shrinkage)
		#  return Fourier image
		#if npad> 1:  data[im] = pad(data[im], Tracker["constants"]["nnxo"]*npad, Tracker["constants"]["nnxo"]*npad, 1, 0.0)

		#  Apply varadj
		if not nonorm: Util.mul_scalar(data[im], Tracker["avgvaradj"][procid]/wnorm)
		#  FT
		data[im] = fft(data[im])
		sig = Util.rotavg_fourier( data[im] )
		Blockdata["accumulatepw"][procid][im] = sig[len(sig)//2:]+[0.0]
		if Tracker["constants"]["CTF"] :
			data[im] = fdecimate(data[im], nxinit*npad, nxinit*npad, 1, False, False)
			ctf_params = original_data[im].get_attr("ctf")
			ctf_params.apix = ctf_params.apix/shrinkage
			data[im].set_attr('ctf', ctf_params)
			data[im].set_attr('ctf_applied', 0)
			if return_real: data[im] = fft(data[im])
		else:
			ctf_params = original_data[im].get_attr_default("ctf", False)
			if ctf_params:
				ctf_params.apix = ctf_params.apix/shrinkage
				data[im].set_attr('ctf', ctf_params)
				data[im].set_attr('ctf_applied', 0)
			data[im] = fdecimate(data[im], nxinit*npad, nxinit*npad, 1, True, False)
			apix = Tracker["constants"]["pixel_size"]
			data[im].set_attr('apix', apix/shrinkage)
			
		#  We have to make sure the shifts are within correct range, shrinkage or not
		set_params_proj(data[im],[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
		if not return_real: data[im].set_attr("padffted",1)
		data[im].set_attr("npad",npad)
		if Blockdata["bckgnoise"]:
			temp = Blockdata["bckgnoise"][data[im].get_attr("particle_group")]
			###  Do not adjust the values, we try to keep everything in the same Fourier values.
			data[im].set_attr("bckgnoise", [temp[i] for i in xrange(temp.get_xsize())])
	return data
###5
def read_data_for_sorting(partids, partstack, previous_partstack):
	# The function will read from stack a subset of images specified in partids
	# and assign to them parameters from partstack with optional CTF application and shifting of the data.
	# So, the lengths of partids and partstack are the same.
	global Tracker, Blockdata
	from fundamentals	import resample, fshift
	from filter			import filt_ctf
	from applications	import MPI_start_end
	from EMAN2          import Region
	from utilities      import wrap_mpi_bcast, read_text_row, get_im, set_params_proj
	# functions:
	# read in data
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if( Blockdata["myid"] == Blockdata["main_node"]):print(line, "read_data for sorting ")
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
	lpartids = wrap_mpi_bcast(lpartids, Blockdata["main_node"])
	groupids = wrap_mpi_bcast(groupids, Blockdata["main_node"])
	Tracker["total_stack"]  = len(lpartids)
	if(Blockdata["myid"] == Blockdata["main_node"]): partstack = read_text_row(partstack)
	else:  partstack = 0
	partstack = wrap_mpi_bcast(partstack, Blockdata["main_node"])
	if(Blockdata["myid"] == Blockdata["main_node"]): previous_partstack = read_text_row(previous_partstack)
	else:  previous_partstack = 0
	previous_partstack = wrap_mpi_bcast(previous_partstack, Blockdata["main_node"])
	if(Tracker["total_stack"] < Blockdata["nproc"]): ERROR("number of processors in use is larger than the total number of particles", \
		  "get_data_and_prep", 1, Blockdata["myid"])
	else: image_start, image_end = MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"])
	lpartids  = lpartids[image_start:image_end]
	groupids  = groupids[image_start:image_end]
	nima      =	image_end - image_start
	data      = [None]*nima
	norm_per_particle = [ None for im in xrange(nima)]
	#print(Tracker["constants"]["orgstack"])
	for im in xrange(nima):
		image = get_im(Tracker["constants"]["orgstack"], lpartids[im])
		#print(im, Blockdata["myid"])
		try: phi, theta, psi, sx, sy, chunk_id, particle_group_id, mnorm = partstack[lpartids[im]][0], \
		   partstack[lpartids[im]][1], partstack[lpartids[im]][2], \
			partstack[lpartids[im]][3], partstack[lpartids[im]][4], partstack[lpartids[im]][5], \
			   partstack[lpartids[im]][6], partstack[lpartids[im]][7]
		except:
			phi, theta, psi, sx, sy, chunk_id, particle_group_id, mnorm  = partstack[lpartids[im]][0], \
			    partstack[lpartids[im]][1], partstack[lpartids[im]][2], \
		  partstack[lpartids[im]][3], partstack[lpartids[im]][4], partstack[lpartids[im]][5], -1, 1.
		sx1, sy1 = previous_partstack[lpartids[im]][3], previous_partstack[lpartids[im]][4]
		set_params_proj(image,[phi, theta, psi, 0.0, 0.0])
		image.set_attr("chunk_id", chunk_id)
		image.set_attr("group", groupids[im])
		image.set_attr("particle_group", particle_group_id)
		#image.set_attr("mnorm", mnorm)
		image.set_attr("previous_shifts", [sx1, sy1])
		image.set_attr("current_shifts", [sx, sy])
		norm_per_particle[im] = mnorm
		data[im] = image
	return data, norm_per_particle
###6 read paramstructure	
def read_paramstructure_for_sorting(partids, paramstructure_dict_file, paramstructure_dir):
	global Tracker, Blockdata
	from utilities    import read_text_row, read_text_file, wrap_mpi_bcast
	from applications import MPI_start_end
	if( Blockdata["myid"] == Blockdata["main_node"]):lcore = read_text_file(partids, -1)
	else: lcore = 0
	lcore   = wrap_mpi_bcast(lcore, Blockdata["main_node"])
	if len(lcore) == 1: lcore = lcore[0]
	else: lcore = lcore[1]
	psize   = len(lcore)
	oldparamstructure  = []	
	im_start, im_end   = MPI_start_end(psize, Blockdata["nproc"], Blockdata["myid"])
	lcore              = lcore[im_start:im_end]
	nima               = len(lcore)
	if( Blockdata["myid"] == Blockdata["main_node"]): tmp_list = read_text_row(paramstructure_dict_file)
	else: tmp_list = 0
	tmp_list = wrap_mpi_bcast(tmp_list, Blockdata["main_node"])
	pdict    = {}
	for im in xrange(len(lcore)):pdict[im] = tmp_list[lcore[im]]
	oldparamstructure             =  []
	nptl                          =  0
	last_old_paramstructure_file  =  None
	while nptl<nima:
		[jason_of_cpu_id, chunk_id, iteration, ptl_id_on_cpu, global_index] = pdict[nptl]
		old_paramstructure_file = os.path.join(paramstructure_dir, "oldparamstructure_%d_%03d_%03d.json"%(chunk_id, jason_of_cpu_id, iteration))
		if old_paramstructure_file != last_old_paramstructure_file:
			fout = open(old_paramstructure_file,'r')
			paramstructure = convert_json_fromunicode(json.load(fout))
			fout.close()
		last_old_paramstructure_file = old_paramstructure_file
		oldparamstructure.append(paramstructure[ptl_id_on_cpu])	
		nptl +=1
	return oldparamstructure
###7 copy oldparamstructures from meridien	
def copy_oldparamstructure_from_meridien_MPI(selected_iteration, log):
	global Tracker, Blockdata
	from utilities    import read_text_row, cmdexecute, write_text_row, read_text_file,wrap_mpi_bcast
	from applications import MPI_start_end
	import json
	Tracker["directory"] = os.path.join(Tracker["constants"]["masterdir"], "main%03d"%selected_iteration)
	Tracker["paramstructure_dir"] = os.path.join(Tracker["directory"], "oldparamstructure")
	old_refinement_iter_directory = os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iteration)
	old_refinement_previous_iter_directory = os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%(selected_iteration-1))
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if( Blockdata["myid"] == Blockdata["main_node"]):
		if os.path.exists(Tracker["paramstructure_dir"]):
			msg = "%s does exist!"%Tracker["paramstructure_dir"]
			print(line, msg)
			log.add(msg)
		else:
			cmd="{} {}".format("mkdir",os.path.join(Tracker["constants"]["masterdir"], "main%03d"%selected_iteration))
			cmdexecute(cmd)  
			cmd="{} {}".format("mkdir",Tracker["paramstructure_dir"])
			cmdexecute(cmd)
			msg = "%s is created"%Tracker["paramstructure_dir"]
			print(line, msg)
			log.add(msg)
	Tracker["refang"] = read_text_row(os.path.join(old_refinement_iter_directory, "refang.txt"))
	if( Blockdata["myid"] == Blockdata["main_node"]): write_text_row(Tracker["refang"], os.path.join(Tracker["directory"], "refang.txt"))
	Tracker["rshifts"] = read_text_row(os.path.join(old_refinement_iter_directory, "rshifts.txt"))
	if( Blockdata["myid"] == Blockdata["main_node"]): write_text_row(Tracker["refang"], os.path.join(Tracker["directory"], "rshifts.txt"))
	my_last_params = read_text_file(os.path.join(old_refinement_previous_iter_directory, "params_%03d.txt"%(selected_iteration-1)), -1)
	my_parstack    = read_text_file(os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt"), -1)
	if( Blockdata["myid"] == Blockdata["main_node"]):
		my_parstack[3:5]= my_last_params[3:5]
		write_text_file(my_parstack, os.path.join(Tracker["constants"]["masterdir"], "previous_refinement_parameters.txt"))
	Tracker["previous_parstack"] = os.path.join(Tracker["constants"]["masterdir"], "previous_refinement_parameters.txt")
	nproc_previous = 0
	procid         = 0
	old_refinement_iter_dir = os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iteration)
	if Blockdata["myid"] == Blockdata["main_node"]:
		while os.path.exists(os.path.join(old_refinement_iter_dir,"oldparamstructure","oldparamstructure_%01d_%03d_%03d.json"%(procid, nproc_previous, selected_iteration))):
			nproc_previous += 1
	nproc_previous = bcast_number_to_all(nproc_previous, Blockdata["main_node"], MPI_COMM_WORLD)
	Blockdata["nproc_previous"] = nproc_previous
	oldparamstructure =[[], []]
	local_dict = {}
	for procid in xrange(2):
		if( Blockdata["myid"] == Blockdata["main_node"]): lcore = read_text_file(os.path.join(Tracker["constants"]["masterdir"], "chunk_%d.txt"%procid))
		else: lcore = 0
		lcore = wrap_mpi_bcast(lcore, Blockdata["main_node"], MPI_COMM_WORLD)	
		psize = len(lcore)
		oldparamstructure[procid] = []
		im_start, im_end   = MPI_start_end(psize, Blockdata["nproc"], Blockdata["myid"])
		local_lcore = lcore[im_start:im_end]
		istart_old_proc_id = -1
		iend_old_proc_id = -1
		plist = []
		for iproc_old in xrange(nproc_previous):
			im_start_old, im_end_old = MPI_start_end(psize, nproc_previous, iproc_old)
			if (im_start>= im_start_old) and im_start <=im_end_old: istart_old_proc_id = iproc_old
			if (im_end>= im_start_old) and im_end <=im_end_old: iend_old_proc_id = iproc_old
			plist.append([im_start_old, im_end_old])
		ptl_on_this_cpu = im_start
		nptl_total = 0
		for iproc_index_old in xrange(istart_old_proc_id, iend_old_proc_id+1):
			fout = open(os.path.join(Tracker["constants"]["refinement_dir"],"main%03d"%selected_iteration, "oldparamstructure", "oldparamstructure_%01d_%03d_%03d.json"%(procid, \
			 iproc_index_old, selected_iteration)),'r')
			oldparamstructure_on_old_cpu = convert_json_fromunicode(json.load(fout))
			fout.close()
			mlocal_id_on_old = ptl_on_this_cpu - plist[iproc_index_old][0]
			while (mlocal_id_on_old<len(oldparamstructure_on_old_cpu)) and (ptl_on_this_cpu<im_end):
				oldparamstructure[procid].append(oldparamstructure_on_old_cpu[mlocal_id_on_old])
				local_dict [local_lcore[nptl_total]] = [Blockdata["myid"], procid, selected_iteration, nptl_total, ptl_on_this_cpu]
				ptl_on_this_cpu  +=1
				mlocal_id_on_old +=1
				nptl_total       +=1
		del oldparamstructure_on_old_cpu
		mpi_barrier(MPI_COMM_WORLD)
		fout = open(os.path.join(Tracker["constants"]["masterdir"], "main%03d"%selected_iteration, "oldparamstructure", "oldparamstructure_%01d_%03d_%03d.json"%(procid, \
		    Blockdata["myid"], selected_iteration)),'w')
		json.dump(oldparamstructure[procid], fout)
		fout.close()
		mpi_barrier(MPI_COMM_WORLD)
		
	if Blockdata["myid"] == Blockdata["main_node"]:
		full_dict_list = [ None for im in xrange(Tracker["constants"]["total_stack"])]
		for key, value in local_dict.iteritems():full_dict_list[key] = value			 
	for icpu in xrange(Blockdata["nproc"]):
		if Blockdata["myid"] == icpu and Blockdata["myid"] != Blockdata["main_node"]: wrap_mpi_send(local_dict, Blockdata["main_node"], MPI_COMM_WORLD)
		elif Blockdata["myid"] != icpu and Blockdata["myid"] == Blockdata["main_node"]:
			local_dict = wrap_mpi_recv(icpu, MPI_COMM_WORLD)
			for key, value in local_dict.iteritems():
				full_dict_list[key] = value
		else: pass
		mpi_barrier(MPI_COMM_WORLD)
	Tracker["paramstructure_dict"] = os.path.join(Tracker["constants"]["masterdir"], "paramstructure_dict.txt")
	if Blockdata["myid"] == Blockdata["main_node"]: write_text_row(full_dict_list, Tracker["paramstructure_dict"])
	return
### 8
def precalculate_shifted_data_for_recons3D(prjlist, paramstructure, refang, rshifts, delta, avgnorms, nxinit, nnxo, nosmearing, norm_per_particle = None, upweighted=False):
	from utilities    import random_string, get_im, findall, info, model_blank
	from filter	      import filt_table
	from fundamentals import fshift
	import types
	import datetime
	import copy
	if norm_per_particle == None: norm_per_particle = len(prjlist)*[1.0]
	nnx = prjlist[0].get_xsize()
	nny = prjlist[0].get_ysize()
	if not nosmearing:
		recdata_list = [[] for im in xrange(len(prjlist))]	
		rshifts_shrank = copy.deepcopy(rshifts)
		for im in xrange(len(rshifts_shrank)):
			rshifts_shrank[im][0] *= float(nxinit)/float(nnxo)
			rshifts_shrank[im][1] *= float(nxinit)/float(nnxo)
		nshifts = len(rshifts_shrank)
	for im in xrange(len(prjlist)):
		bckgn = prjlist[im].get_attr("bckgnoise")
		ct = prjlist[im].get_attr("ctf")
		group_id = prjlist[im].get_attr("group")
		if nosmearing:
			phi,theta,psi,s2x,s2y = get_params_proj(prjlist[im], xform = "xform.projection")
			prjlist[im].set_attr("wprob", 1.0)
			prjlist[im].set_attr("group", group_id)
			prjlist[im].set_attr_dict( {"bckgnoise":bckgn, "ctf":ct})
			prjlist[im].set_attr_dict({"padffted":1, "is_fftpad":1,"is_fftodd":0, "is_complex_ri":1, "is_complex":1})
			if not upweighted:prjlist[im] = filt_table(prjlist[im], bckgn)
			set_params_proj(prjlist[im],[ phi, theta, psi, 0.0, 0.0], xform = "xform.projection")
		else:
			avgnorm = avgnorms[prjlist[im].get_attr("chunk_id")]
			numbor      = len(paramstructure[im][2])
			ipsiandiang = [ paramstructure[im][2][i][0]/1000  for i in xrange(numbor) ]
			allshifts   = [ paramstructure[im][2][i][0]%1000  for i in xrange(numbor) ]
			probs       = [ paramstructure[im][2][i][1] for i in xrange(numbor) ]
			tdir = list(set(ipsiandiang))
			data = [None]*nshifts
			for ii in xrange(len(tdir)):
				#  Find the number of times given projection direction appears on the list, it is the number of different shifts associated with it.
				lshifts = findall(tdir[ii], ipsiandiang)
				toprab  = 0.0
				for ki in xrange(len(lshifts)):toprab += probs[lshifts[ki]]
				recdata = EMData(nny,nny,1,False)
				recdata.set_attr("is_complex",0)
				for ki in xrange(len(lshifts)):
					lpt = allshifts[lshifts[ki]]
					if(data[lpt] == None):
						data[lpt] = fshift(prjlist[im], rshifts_shrank[lpt][0], rshifts_shrank[lpt][1])
						data[lpt].set_attr("is_complex",0)
					Util.add_img(recdata, Util.mult_scalar(data[lpt], probs[lshifts[ki]]/toprab))
				recdata.set_attr_dict({"padffted":1, "is_fftpad":1,"is_fftodd":0, "is_complex_ri":1, "is_complex":1}) # preset already
				if not upweighted:recdata = filt_table(recdata, bckgn)
				recdata.set_attr_dict( {"bckgnoise":bckgn, "ctf":ct})
				ipsi = tdir[ii]%100000
				iang = tdir[ii]/100000
				set_params_proj(recdata,[refang[iang][0],refang[iang][1], refang[iang][2]+ipsi*delta, 0.0, 0.0], xform = "xform.projection")
				recdata.set_attr("wprob", toprab*avgnorm/norm_per_particle[im])
				recdata.set_attr("group", group_id)
				recdata_list[im].append(recdata)
	if nosmearing:return prjlist
	else:
		del bckgn, recdata, tdir, ipsiandiang, allshifts, probs, data
		return recdata_list
##### read data/paramstructure ends
###<<<----downsize data---->>>>
def downsize_data_for_sorting(original_data, return_real = False, preshift = True, npad = 1):
	# The function will read from stack a subset of images specified in partids
	#   and assign to them parameters from partstack with optional CTF application and shifting of the data.
	# So, the lengths of partids and partstack are the same.
	global Tracker, Blockdata
	from fundamentals	import resample, fshift, cyclic_shift
	from filter			import filt_ctf
	from applications	import MPI_start_end
	from EMAN2          import Region
	# functions:
	# read in data
	# apply mask, and prepare focus projection if focus3D is specified
	# return  1. cdata: data for image comparison, always in Fourier format
	#         2. rdata: data for reconstruction, 4nn return real image
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if( Blockdata["myid"] == Blockdata["main_node"]): print(line,"--->>>downsize_data_for_sorting starts<<<---- ")		
	mask2D		= model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	shrinkage 	= Tracker["nxinit"]/float(Tracker["constants"]["nnxo"])
	radius      = int(Tracker["constants"]["radius"] * shrinkage +0.5)
	if Tracker["applybckgnoise"]:
		oneover = []
		nnx     = len(Blockdata["bckgnoise"][0])
		for i in xrange(len(Blockdata["bckgnoise"])):
			#print(i, Blockdata["myid"])
			temp = [0.0]*nnx
			for k in xrange(nnx):
				if(Blockdata["bckgnoise"][i][k] > 0.0): temp[k] = 1.0/sqrt(Blockdata["bckgnoise"][i][k])
			oneover.append(temp)
		del temp
	if Tracker["focus3D"]: # focus mask is applied
		if Blockdata["myid"] == Blockdata["main_node"]:
			focus3d    = get_im(Tracker["focus3D"])
			focus3d_nx = focus3d.get_xsize()
			if focus3d_nx != Tracker["nxinit"]: # So the decimated focus volume can be directly used
				focus3d = resample(focus3d, float(Tracker["nxinit"])/float(focus3d_nx))
		else: focus3d = model_blank(Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"])
		bcast_EMData_to_all(focus3d, Blockdata["myid"], Blockdata["main_node"])
		focus3d = prep_vol(focus3d, 1, 1)
	#  Preprocess the data
	nima   = len(original_data)
	cdata  = [None]*nima
	rdata  = [None]*nima	
	for im in xrange(nima):
		image = original_data[im].copy()
		chunk_id = image.get_attr("chunk_id")
		try:    group_id = image.set_attr("group", groupids[im])
		except: pass
		particle_group_id = image.get_attr("particle_group")
		phi,theta,psi,s2x,s2y = get_params_proj(image, xform = "xform.projection")
		[sx, sy]   = image.get_attr("previous_shifts")
		[sx1, sy1] = image.get_attr("current_shifts")
		rimage = cyclic_shift(image, int(round(sx)), int(round(sy)))
		cimage = fshift(image, sx1, sy1)		
		st = Util.infomask(rimage, mask2D, False)
		rimage -= st[0]
		rimage /= st[1]
		st = Util.infomask(cimage, mask2D, False)
		cimage -= st[0]
		cimage /= st[1]
		if Tracker["applybckgnoise"]:
			if Tracker["applymask"]:
				if Tracker["constants"]["hardmask"]: cimage = cosinemask(cimage, radius = Tracker["constants"]["radius"])
				else:
					bckg = model_gauss_noise(1.0,Tracker["constants"]["nnxo"]+2,Tracker["constants"]["nnxo"])
					bckg.set_attr("is_complex",1)
					bckg.set_attr("is_fftpad",1)
					bckg = fft(filt_table(bckg, oneover[particle_group_id]))
					#  Normalize bckg noise in real space, only region actually used.
					st = Util.infomask(bckg, mask2D, False)
					bckg -= st[0]
					bckg /= st[1]
					cimage = cosinemask(cimage,radius = Tracker["constants"]["radius"], bckg = bckg)
		else:
			if Tracker["applymask"]:cimage  = cosinemask(cimage, radius = Tracker["constants"]["radius"])
			else: pass
		# FT
		rimage  = fft(rimage)
		cimage  = fft(cimage)
		if Tracker["constants"]["CTF"] :
			ctf_params = rimage.get_attr("ctf")
			rimage      = fdecimate(rimage, Tracker["nxinit"]*npad, Tracker["nxinit"]*npad, 1, False, False)
			cimage     = fdecimate(cimage, Tracker["nxinit"]*npad, Tracker["nxinit"]*npad, 1, False, False)
			ctf_params.apix = ctf_params.apix/shrinkage
			rimage.set_attr('ctf', ctf_params)
			cimage.set_attr('ctf', ctf_params)
			rimage.set_attr('ctf_applied', 0)
			cimage.set_attr('ctf_applied', 0)
			if return_real :  rimage = fft(rimage)
		else:
			ctf_params = rimage.get_attr_default("ctf", False)
			if  ctf_params:
				ctf_params.apix = ctf_params.apix/shrinkage
				rimage.set_attr('ctf', ctf_params)
				rimage.set_attr('ctf_applied', 0)
				cimage.set_attr('ctf', ctf_params)
				cimage.set_attr('ctf_applied', 0)
				
			rimage  = fdecimate(rimage, nxinit*npad, nxinit*npad, 1, True, False)
			cimage = fdecimate(cimage, nxinit*npad, nxinit*npad, 1, True, False)
			apix   = Tracker["constants"]["pixel_size"]
			rimage.set_attr('apix', apix/shrinkage)
			cimage.set_attr('apix', apix/shrinkage)
		
		cimage.set_attr("padffted",1)
		cimage.set_attr("npad", npad)
		if not return_real:	
			rimage.set_attr("padffted",1)
			rimage.set_attr("npad", npad)
			
		set_params_proj(rimage,[phi, theta, psi, 0.0, 0.0])
		rimage.set_attr("chunk_id", chunk_id)
		#image.set_attr("group", groupids[im])
		rimage.set_attr("particle_group", particle_group_id)
		
		set_params_proj(cimage,[phi, theta, psi, 0.0, 0.0])
		cimage.set_attr("chunk_id", chunk_id)
		#cimage.set_attr("group", groupids[im])
		cimage.set_attr("particle_group", particle_group_id)
		rdata[im] =  rimage
		cdata[im] =  cimage		
		if Tracker["applybckgnoise"]: 
			rdata[im].set_attr("bckgnoise", Blockdata["bckgnoise"][particle_group_id])
			if Tracker["constants"]["comparison_method"] == "cross":Util.mulclreal(cdata[im], Blockdata["unrolldata"][particle_group_id])
		else:
			rdata[im].set_attr("bckgnoise",  Blockdata["bckgnoise"])
			cdata[im].set_attr("bckgnoise",  Blockdata["bckgnoise"])                    
		if Tracker["focus3D"]:
			cdata[im] = fft(binarize(prgl(focus3d, [phi, theta, psi, 0.0, 0.0], 1, True), 1)*fft(cdata[im]))
			if Tracker["constants"]["CTF"]: cdata[im].set_attr("ctf", rdata[im].get_attr("ctf"))
		cdata[im].set_attr("is_complex",0)	
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if( Blockdata["myid"] == Blockdata["main_node"]):print(line,"downsizing data for sorting finishes")
	return cdata, rdata
##<<<----for 3D----->>>>
def downsize_data_for_rec3D(original_data, particle_size, return_real = False, npad = 1):
	# The function will read from stack a subset of images specified in partids
	#   and assign to them parameters from partstack with optional CTF application and shifting of the data.
	# So, the lengths of partids and partstack are the same.	
	global Tracker, Blockdata
	from fundamentals	import resample, fshift
	from filter			import filt_ctf
	from applications	import MPI_start_end
	from EMAN2          import Region
	# functions:
	# read in data
	# apply mask, and prepare focus projection if focus3D is specified
	# return  1. cdata: data for image comparison, always in Fourier format
	#         2. rdata: data for reconstruction, 4nn return real image
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if( Blockdata["myid"] == Blockdata["main_node"]): print(line, "--->>>downsize_data_for_rec3D starts<<<---- ")		
	nima       = len(original_data)
	rdata      = [None]*nima
	mask2D     = model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	shrinkage  = particle_size/float(Tracker["constants"]["nnxo"])
	radius     = int(Tracker["constants"]["radius"] * shrinkage +0.5)
	for im in xrange(nima):
		image = original_data[im].copy()
		chunk_id = image.get_attr("chunk_id")
		try: particle_group_id = image.get_attr("particle_group")
		except: particle_group_id = -1
		phi,theta,psi,s2x,s2y = get_params_proj(image, xform = "xform.projection")
		[sx, sy] = image.get_attr("previous_shifts") # always for rec3D
		if Tracker["nosmearing"]: image = fshift(image, s2x, s2y)
		else: image = cyclic_shift(image, int(round(sx)), int(round(sy)))
		st = Util.infomask(image, mask2D, False)
		image -= st[0]
		image /= st[1]
		image  = fft(image)
		if Tracker["constants"]["CTF"]:
			ctf_params = image.get_attr("ctf")
			image = fdecimate(image, particle_size*npad, particle_size*npad, 1, False, False)
			ctf_params.apix = ctf_params.apix/shrinkage
			image.set_attr('ctf', ctf_params)
			image.set_attr('ctf_applied', 0)
		else:
			ctf_params = image.get_attr_default("ctf", False)
			if  ctf_params:
				ctf_params.apix = ctf_params.apix/shrinkage
				image.set_attr('ctf', ctf_params)
				image.set_attr('ctf_applied', 0)
			image = fdecimate(image, particle_size*npad, particle_size*npad, 1, True, False)
			apix  = Tracker["constants"]["pixel_size"]
			image.set_attr('apix', apix/shrinkage)		
		if not return_real:	
			image.set_attr("padffted",1)
			image.set_attr("npad", npad)
		image.set_attr("chunk_id", chunk_id)
		image.set_attr("particle_group", particle_group_id)
		set_params_proj(image,[phi, theta, psi, 0.0, 0.0])
		rdata[im] =  image
		if Tracker["applybckgnoise"]: rdata[im].set_attr("bckgnoise", Blockdata["bckgnoise"][rdata[im].get_attr("particle_group")])
		else: rdata[im].set_attr("bckgnoise", Blockdata["bckgnoise"])
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"                  
	if( Blockdata["myid"] == Blockdata["main_node"]):print(line,"--->>>downsize_data_for_rec3D finishes<<<---- ")
	return rdata
### end of downsize
###<<<--- comparison	    
def compare_two_images_eucd(data, ref_vol):
	global Tracker, Blockdata
	from filter import filt_tophatl
	from math   import sqrt
	peaks   = len(data)*[None]
	ny      = data[0].get_ysize()
	ref_vol = prep_vol(ref_vol, npad = 2, interpolation_method = 1)
	ctfs    = [ctf_img_real(ny, q.get_attr('ctf')) for q in data]
	qt = float(Tracker["constants"]["nnxo"]*Tracker["constants"]["nnxo"])
	for im in xrange(len(data)):
		phi, theta, psi, s2x, s2y = get_params_proj(data[im], xform = "xform.projection")
		rtemp = prgl(ref_vol,[phi, theta, psi, 0.0,0.0], 1, False)
		rtemp.set_attr("is_complex",0)
		if data[im].get_attr("is_complex") ==1: data[im].set_attr("is_complex",0)
		if Tracker["applybckgnoise"]:
			peaks[im] = -Util.sqed(data[im], rtemp, ctfs[im], Blockdata["unrolldata"][data[im].get_attr("particle_group")])/qt
		else: peaks[im] = -Util.sqed(data[im], rtemp, ctfs[im], Blockdata["unrolldata"])/qt
	return peaks
#
def compare_two_images_cross(data, ref_vol):
	global Tracker, Blockdata
	from filter import filt_tophatl
	from math   import sqrt
	ny = data[0].get_ysize()
	peaks   = len(data)*[None]
	volft   = prep_vol(ref_vol, 2, 1)
	ctfs  = [None for im in xrange(len(data))]
	for im in xrange(len(data)): ctfs[im]  = ctf_img_real(ny, data[im].get_attr('ctf'))
	#  Ref is in reciprocal space
	for im in xrange(len(data)):
		phi, theta, psi, s2x, s2y = get_params_proj(data[im], xform = "xform.projection")
		ref = prgl( volft, [phi, theta, psi, 0.0, 0.0], 1, False)
		Util.mulclreal(ref, ctfs[im])
		ref.set_attr("is_complex", 0)
		ref.set_value_at(0,0,0.0)
		nrmref = sqrt(Util.innerproduct(ref, ref, None))
		if data[im].get_attr("is_complex") ==1: data[im].set_attr("is_complex",0)
		if not Tracker["focus3D"]:
			if Tracker["applybckgnoise"]:  peak = Util.innerproduct(ref, data[im], Blockdata["unrolldata"][data[im].get_attr("particle_group")])
			else: peak = Util.innerproduct(ref, data[im], None)
			peaks[im] = peak/nrmref
		else:  peaks[im] = Util.innerproduct(ref, data[im], None)/nrmref
	return peaks
###<<<---various utilities
def clusters_to_plist(clusters, pall):
	# clusters contains the original ids
	pdict = {}
	plist = []
	qlist = []
	for igrp in xrange(len(clusters)):
		clusters[igrp].tolist()
		for a in clusters[igrp]: 
			pdict[pall[a]] = igrp
			plist.append(pall[a])
	plist      = sorted(plist)
	assignment = [ None for im in xrange(len(plist))]
	for im in xrange(len(plist)):
		assignment[im] = pdict[plist[im]]
		qlist.append([pdict[plist[im]], plist[im]])
	a =  set(pall)
	b =  set(plist)
	unaccounted = list(a.difference(b))
	return [assignment, plist], qlist, unaccounted

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
		else: Tracker["sorting_data_list"]= data_list[1]
		if Tracker["constants"]["seed"] == -1: random.seed()
		else: random.seed(Tracker["constants"]["seed"])
		Tracker["indep_runs_list"]  = []
		group_size = len(Tracker["sorting_data_list"])//Tracker["number_of_groups"]
		for index_of_random in xrange(Tracker["constants"]["indep_runs"]):
			particle_dict = {}
			ll = copy.deepcopy(Tracker["sorting_data_list"])
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
			del assignment
			del ll
	else:
		Tracker["indep_runs_list"]   = 0
		Tracker["sorting_data_list"] = 0
	Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"])
	return
	
def resize_groups_from_stable_members_mpi(Accounted_on_disk, Unaccounted_on_disk):
	global Tracker, Blockdata
	import random
	if Blockdata["myid"] == Blockdata["main_node"]:
		ptl_dict   = {}
		accounted  = read_text_file(Accounted_on_disk, -1)
		number_of_groups = max(accounted[0]) + 1
		groups = [[] for igrp in xrange(number_of_groups)]
		unaccounted = read_text_file(Unaccounted_on_disk, -1)
		for im in xrange(len(accounted[0])):
			groups[accounted[0][im]].append(accounted[1][im])
			ptl_dict[accounted[1][im]] = accounted[0][im]
		accounted_members = sorted(accounted[1])
		unaccounted = sorted(unaccounted[0])
		full_list = accounted_members + unaccounted
		full_list = sorted(full_list)
		total_stack = len(full_list)
		Tracker["number_of_groups"] = number_of_groups
		group_size = int(float(total_stack)/number_of_groups)
	else:
		number_of_groups = 0
		total_stack      = 0
	number_of_groups = bcast_number_to_all(number_of_groups, Blockdata["main_node"], MPI_COMM_WORLD)
	total_stack = bcast_number_to_all(total_stack, Blockdata["main_node"], MPI_COMM_WORLD)
	Tracker["total_stack"] = total_stack
	Tracker["number_of_groups"] = number_of_groups
	###-----
	if Blockdata["myid"] == Blockdata["main_node"]:
		assignment_list = []
		for indep in xrange(1):
			print("after iter %d"%indep)
			unaccounted_members = copy.deepcopy(unaccounted)
			new_groups = copy.deepcopy(groups)
			for im in xrange(len(unaccounted_members)):
				igroup = random.randrange(0, number_of_groups)	
				new_groups[igroup].append(unaccounted_members[im])
				ptl_dict[unaccounted_members[im]] = igroup
			"""
			shuffle(unaccounted_members)
			l1 = 0
			while l1 <len(unaccounted_members):
				for igrp in xrange(number_of_groups):
					if len(new_groups[igrp])< group_size + 1:
						new_groups[igrp].append(unaccounted_members[l1])
						ptl_dict[unaccounted_members[l1]] = igrp
						l1 +=1
						print(l1, len(unaccounted_members), len(new_groups[igrp]))
			"""
			print ("unaccounted, resize", len(unaccounted_members))
			assignment = [None for iptl in xrange(len(full_list))]
			for im in xrange(len(full_list)): assignment[im] = ptl_dict[full_list[im]]
			assignment_list.append([assignment,full_list])
			del unaccounted_members
	else: assignment_list = 0
	mpi_barrier(MPI_COMM_WORLD)
	assignment_list = wrap_mpi_bcast(assignment_list, Blockdata["main_node"], MPI_COMM_WORLD)
	if Blockdata["myid"] == Blockdata["main_node"]: del ptl_dict
	return assignment_list
			
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
		log_main.add(msg)
		print(msg)
	if(Tracker["constants"]["indep_runs"] < 2): ERROR("There is only one independent run!", "do_two_way_comparison_over_nindepruns", 1, Blockdata["myid"])
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
		if Blockdata["myid"] == ires%Blockdata["nproc"]: tmplist = res_list[ires][:][:]
		else: tmplist = 0
		tmplist = wrap_mpi_bcast(tmplist, ires, MPI_COMM_WORLD)
		if Blockdata["myid"] == Blockdata["main_node"]:res_list[ires] = tmplist[:][:]
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
	duplicate = {}
	for ipp in xrange(len(ptp)):
		avg_scores = 0.0
		for jpp in xrange(len(ptp)):
			if ipp != jpp:	avg_scores += full_scores[(ipp, jpp)]
		summed_scores[ipp]  = avg_scores/(len(ptp)-1.0)
		try:
			kpp = two_way_dict[summed_scores[ipp]]
			duplicate[summed_scores[ipp]] = ipp
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
		tmp_run1 = two_way_dict[rate1]
		tmp_run2 = duplicate[rate2]
		run1 = min(tmp_run1, tmp_run2)
		run2 = max(tmp_run1, tmp_run2)
	try: Tracker["selected_stable_class"] = Tracker["stable_classes"][(run1, run2)]
	except:
		try: Tracker["selected_stable_class"] = Tracker["stable_classes"][(run2, run1)]
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
		for icluster in xrange(len(Tracker["selected_stable_class"])):write_text_file(Tracker["selected_stable_class"][icluster], os.path.join(Tracker["directory"], "Class%d.txt"%icluster))
		log_main.add("accounted:        %8d"%len(Tracker["accounted_list"]))
		log_main.add("unaccounted:      %8d"%len(Tracker["unaccounted_list"]))
		log_main.add("total:            %8d"%(len(Tracker["accounted_list"])+len(Tracker["unaccounted_list"])))
		write_text_file(Tracker["unaccounted_list"], os.path.join(Tracker["directory"], "Unaccounted.txt"))
		write_text_row(Tracker["accounted_list"], os.path.join(Tracker["directory"], "Accounted.txt"))
		Tracker["Accounted_on_disk"]   = os.path.join(Tracker["directory"], "Accounted.txt")
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
#####
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
	global Tracker, Blockdata
	# split groupids from indexes of particles
	# reindex groups
	clusters          = []
	cluster_id        = []
	for im in xrange(len(partition)):
		if  partition[im][0] not in cluster_id:cluster_id.append(partition[im][0])
	####
	cluster_dict      = {}
	group_change_dict = {}
	new_group_id      = 0
	for icluster in xrange(len(cluster_id)):
		one_cluster = []
		for a in partition:
			if a[0]== icluster: 
				one_cluster.append(a[1])
				cluster_dict[a[1]] = icluster
		if len(one_cluster)>= Tracker["constants"]["minimum_grp_size"]: # clean small ones
			clusters.append(one_cluster)
			group_change_dict[icluster] = new_group_id
			new_group_id +=1
		else:	group_change_dict[icluster] = -1

	# create a partition list:
	new_partition = [] 
	for iptl in xrange(len(partition)):
		gid = group_change_dict[cluster_dict[partition[iptl][1]]]
		if gid >-1: new_partition.append([group_change_dict[cluster_dict[partition[iptl][1]]], partition[iptl][1]])
	return clusters, new_partition
	 
def prep_ptp_single(all_lists, full_list):
	# full_list contains the initial input indexes
	# the assignment is aligned to full_list
	# convert classes into a single list ptp denoted by group id
	ad_hoc_group_ID = len(all_lists)+1
	ad_hoc_particle_exists = False
	a = set([])
	for b in all_lists: a.union(b)
	c = set(full_list)
	if list(a.difference(c)) !=[]: ERROR("Accounted and unaccounted in total do not match the total number of particles", "prep_ptp_single", 1, Blockdata["myid"])
	else:
		pdict = {}
		for iclass in xrange(len(all_lists)):
			for iptl in xrange(len(all_lists[iclass])): pdict[all_lists[iclass][iptl]] = iclass
		assignment = []
		for im in xrange(len(full_list)):
			#pdict[full_list[im]]
			try: group_ID =  pdict[full_list[im]]
			except:
				group_ID = ad_hoc_group_ID
				ad_hoc_particle_exists = True
			assignment.append(group_ID)
		if ad_hoc_particle_exists: ptp = convertasi(assignment, ad_hoc_group_ID)
		else: ptp = convertasi(assignment, len(all_lists)+1)
		del pdict
	return ptp

def merge_original_id_lists(original_id_lists):
	# merge particles lists with original ID into one list while stamped by new group id 
	clusters     = []
	all_id_list  = []
	for index_of_list in xrange(len(original_id_lists)):
		cluster_dict ={}
		for index_of_particle in xrange(len(original_id_lists[index_of_list])):
			try: cluster_dict [original_id_lists[index_of_list][index_of_particle][0]].append(original_id_lists[index_of_list][index_of_particle][1])
			except: cluster_dict [original_id_lists[index_of_list][index_of_particle][0]]= [original_id_lists[index_of_list][index_of_particle][1]]
			all_id_list.append(original_id_lists[index_of_list][index_of_particle][1])
		for a in cluster_dict: clusters.append(cluster_dict[a])
		del cluster_dict
	all_id_list = sorted(all_id_list)
	cluster_dict = {}
	for index_of_cluster in xrange(len(clusters)):
		for index_of_particle in xrange(len(clusters[index_of_cluster])): cluster_dict[clusters[index_of_cluster][index_of_particle]] = index_of_cluster
	final_list = []
	for index_of_particle in xrange(len(all_id_list)):final_list.append([cluster_dict[all_id_list[index_of_particle]], all_id_list[index_of_particle]])
	del cluster_dict
	return final_list, len(clusters)

def merge_classes_into_partition_list(classes_list):
	group_dict = {}
	data_list  = []
	new_index  = []
	for index_of_class in xrange(len(classes_list)):
		for index_of_particle in xrange(len(classes_list[index_of_class])):
			data_list.append(classes_list[index_of_class][index_of_particle])
			group_dict[classes_list[index_of_class][index_of_particle]] = index_of_class
	data_list = sorted(data_list)
	for index_of_particle in xrange(len(data_list)):new_index.append([group_dict[data_list[index_of_particle]], data_list[index_of_particle]])
	del group_dict
	return data_list, new_index
		
def get_sorting_all_params(data):
	global Tracker, Blockdata
	from utilities    import wrap_mpi_bcast
	from applications import MPI_start_end
	if Blockdata["myid"] == Blockdata["main_node"]:	total_attr_value_list = [[]]*Tracker["total_stack"]
	else: total_attr_value_list = 0
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
	for idat in xrange(len(data_in_core)): attr_value_list.append([data_in_core[idat].get_attr("group"), get_params_proj(data_in_core[idat],xform = "xform.projection")])
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
	else: ERROR("Wrong columns", "parsing_sorting_params", 1, 0)
	return group_list, ali3d_params_list

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
		for index_of_particle in xrange(len(partition_list[ipt])): assignment[index_of_particle] = partition_list[ipt][index_of_particle][0]
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
		if res_curve[ifreq] <0.5: break
	res_05  = ifreq - 1
	for ifreq in xrange(1, len(res_curve)):
		if res_curve[ifreq]<0.143: break
	res_143  = ifreq - 1
	return res_05, res_143
	
def extract_clusters_from_partition(partition_to_be_saved, number_of_cluster):
	clusters = []
	for i in xrange(number_of_cluster):
		clusters.append([])
	for ipar in xrange(len(partition_to_be_saved)):
		[cluster_ID, original_ID] = partition_to_be_saved[ipar]
		clusters[cluster_ID].append(original_ID)
	for icluster in xrange(len(clusters)): clusters[icluster] = sorted(clusters[icluster])
	return clusters
		
def update_data_partition(cdata, rdata, partids):
	# update particle clustering partitions of independent EQKmeans run
	global Tracker, Blockdata
	from utilities import wrap_mpi_bcast
	import copy
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
	lpartids = wrap_mpi_bcast(lpartids, Blockdata["main_node"])
	groupids = wrap_mpi_bcast(groupids, Blockdata["main_node"])
	assignment = copy.copy(groupids)
	try: assert(Tracker["total_stack"] == len(groupids))
	except: ERROR("total stack in Tracker does not agree with the one is just read in", "update_data_partition", 1, Blockdata["myid"])
	image_start, image_end = MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"])	
	nima = image_end - image_start
	assert(nima == len(cdata))
	groupids  = groupids[image_start:image_end]
	for im in xrange(nima):
		cdata[im].set_attr("group",groupids[im])
		rdata[im].set_attr("group",groupids[im])
	return assignment
	
def partition_data_into_orientation_groups_nompi(refa_vecs, data_vecs):
	orien_assignment = [ None for im in xrange(len(data_vecs))]
	for im in xrange(len(data_vecs)):
		max_dist = -999.0
		for jm in xrange(len(refa_vecs)):
			this_dis = get_dist1(data_vecs[im], refa_vecs[jm])
			if this_dis > max_dist:
				max_dist = this_dis
				orien_assignment[im] = jm
	return orien_assignment
	
### dmatrix and refangles partition
def get_dist1(vec1, vec2):
	sum_dot = 0.0
	for icomp in xrange(len(vec1)): sum_dot +=vec1[icomp]*vec2[icomp]
	return sum_dot

def get_refa_vecs(angle_step, sym, tilt1, tilt2):
	from utilities import even_angles, getvec
	refa          = even_angles(angle_step, symmetry = sym, theta1= tilt1,  theta2 = tilt2, method='S', phiEqpsi="Zero")
	refa_vecs     = []
	for i in xrange(len(refa)):
		tmp = getvec(refa[i][0], refa[i][1])
		refa_vecs.append(tmp)
	return refa_vecs

def find_neighborhood(refa_vecs, minor_groups):
	matched_oriens  = [ [None, None] for i in xrange(len(minor_groups))]
	for iproj in xrange(len(minor_groups)):
		max_value = -999.0
		for jproj in xrange(len(refa_vecs)):
			if jproj not in minor_groups:
				this_dis = get_dist1(refa_vecs[minor_groups[iproj]], refa_vecs[jproj])
				if this_dis > max_value:
					max_value = this_dis
					matched_oriens[iproj] = [minor_groups[iproj], jproj]
	return matched_oriens
	
def reassign_ptls_in_orien_groups(assigned_ptls_in_groups, matched_pairs):
	tmplist = []
	for iorien in xrange(len(matched_pairs)):
		assigned_ptls_in_groups[matched_pairs[iorien][1]] +=assigned_ptls_in_groups[matched_pairs[iorien][0]]
		tmplist.append(matched_pairs[iorien][0])
	reassignment = []
	for iorien in xrange(len(assigned_ptls_in_groups)):
		if iorien not in tmplist: reassignment.append(sorted(assigned_ptls_in_groups[iorien]))
	return reassignment

def findall_dict(value, L, start=0):
	"""
	 return a list of all indices of a value on the list L beginning from position start
	"""
	positions = []
	lL = len(L) - 1
	i = start - 1
	while(i < lL):
		i +=1
		try:
			if value == L[i]: positions.append(i) 
		except: pass 
	return positions
	
def do_assignment_by_dmatrix_orien_group(peaks, orien_group_members, number_of_groups):
	import numpy as np
	import random
	#print(len(orien_group_members))
	nima    = len(orien_group_members)
	dmatrix = [None]*(nima*number_of_groups)
	for im in xrange(nima):
		for iref in xrange(number_of_groups):
			dmatrix[iref*nima+im] = peaks[iref][orien_group_members[im]]*(-1.)
	# do dmatrix
	dd = np.argsort(dmatrix)
	maxasi     = nima/number_of_groups
	ngs        = number_of_groups
	id_list    = [ [] for i in xrange(number_of_groups)]
	del_row    = [ False for i in xrange(number_of_groups)]
	del_column = [ False for i in xrange(nima)]
	walktrough = 0
	while ngs > 0:
		flag = True
		while flag:
			l =  dd[walktrough]
			igroup  =  l/nima
			iptl    =  l%nima
			if del_row[igroup] or del_column[iptl]: walktrough +=1
			else:  flag = False
		id_list[igroup].append(iptl)
		if ngs>1:
			if (len(id_list[igroup]) < maxasi): igroup = -1
			else: ngs -= 1
		else:
			if len(id_list[igroup]) < maxasi+nima%number_of_groups: igroup = -1
			else: ngs -= 1
		del_column[iptl] = True
		if (igroup != -1): del_row[igroup] = True
	
	id_list1 = []
	for iref in xrange(number_of_groups):
		for im in xrange(maxasi):
			id_list1.append(id_list[iref][im])
			
	if nima%number_of_groups !=0:
		for im in xrange(maxasi, maxasi+ nima%number_of_groups):
			id_list1.append(id_list[igroup][im])
		id_list1.append(igroup)
		
	id_list = [[] for i in xrange(number_of_groups)]
	maxasi = nima/number_of_groups
	for iref in xrange(maxasi*number_of_groups):
		id_list[iref/maxasi].append(id_list1[iref])
	if nima%number_of_groups !=0:
		for iptl in xrange(nima%maxasi):
			id_list[id_list1[-1]].append(id_list1[maxasi*number_of_groups+iptl])
	for iref in xrange(number_of_groups):
		id_list[iref].sort()
	del id_list1
	assignment = [None]*nima
	for iref in xrange(number_of_groups):
		for im in id_list[iref]: assignment[im] = iref
	del dmatrix
	del dd
	del id_list
	del del_column
	del del_row
	return assignment
	
def get_orien_assignment_mpi(angle_step, partids, params, log_main):
	global Tracker, Blockdata
	from applications import MPI_start_end
	from utilities    import even_angles, wrap_mpi_recv, wrap_mpi_bcast, wrap_mpi_send, read_text_row, read_text_file, getvec
	if Blockdata["myid"] == Blockdata["main_node"]:
		msg = " Generate sampling orientations for EQKmeans with step %f   theta1  %f  theta2  %f"%(Tracker["angle_step"], Tracker["tilt1"], Tracker["tilt2"])
		log_main.add(msg)
		print(msg)
	image_start, image_end = MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"])
	if Blockdata["myid"] == Blockdata["main_node"]:
		orien_group_assignment = [None for im in xrange(Tracker["total_stack"])]
	else:  orien_group_assignment = 0
	#refa = even_angles(angle_step, symmetry = Tracker["constants"]["symmetry"], \
	#     theta1 = Tracker["tilt1"], theta2 = Tracker["tilt2"], method='S', phiEqpsi="Zero")
	refa = Blockdata["symclass"].even_angles(angle_step, theta1 = Tracker["tilt1"], theta2 = Tracker["tilt2"])
	#print(refa)
	refa_vecs = []
	for i in xrange(len(refa)):
		tmp = getvec(refa[i][0], refa[i][1])
		refa_vecs.append(tmp)
	if Blockdata["myid"] == Blockdata["main_node"]: print("step1")	
	if Blockdata["main_node"] == Blockdata["myid"]:
		params  = read_text_row(params)
		partids = read_text_file(partids, -1)
		if len(partids) == 1: partids = partids[0]
		else: partids = partids[1]
		data_angles = [[None,None] for im in xrange(len(partids))]
		for im in xrange(len(partids)): 
			data_angles[im] = getvec(params[partids[im]][0], params[partids[im]][1])
		del params
		del partids
	else: data_angles = 0
	data_angles = wrap_mpi_bcast(data_angles, Blockdata["main_node"], MPI_COMM_WORLD)
	data_angles = data_angles[image_start: image_end]
	if Blockdata["myid"] == Blockdata["main_node"]: print("step2")
	local_orien_group_assignment = partition_data_into_orientation_groups_nompi(refa_vecs, data_angles)
	if Blockdata["myid"] == Blockdata["main_node"]: orien_group_assignment[image_start:image_end] = local_orien_group_assignment[:]
	else:  orien_group_assignment = 0	
	if Blockdata["main_node"] != Blockdata["myid"]: wrap_mpi_send(local_orien_group_assignment, Blockdata["main_node"], MPI_COMM_WORLD)
	else:
		for iproc in xrange(Blockdata["nproc"]):
			iproc_image_start, iproc_image_end = MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], iproc)
			if iproc != Blockdata["main_node"]:
				dummy = wrap_mpi_recv(iproc, MPI_COMM_WORLD)
				orien_group_assignment[iproc_image_start:iproc_image_end] = dummy[:]
				del dummy
	mpi_barrier(MPI_COMM_WORLD)
	orien_group_assignment = wrap_mpi_bcast(orien_group_assignment, Blockdata["main_node"], MPI_COMM_WORLD)
	if Blockdata["myid"] == Blockdata["main_node"]: print("step3")
	ptls_in_orien_groups = [ None for iref in xrange(len(refa_vecs))]
	for iorien in xrange(len(refa_vecs)):
		if iorien%Blockdata["nproc"] == Blockdata["myid"]: ptls_in_orien_groups[iorien] = findall_dict(iorien, orien_group_assignment)
	mpi_barrier(MPI_COMM_WORLD)
	if Blockdata["myid"] == Blockdata["main_node"]:print("step4")
	for iorien in xrange(len(refa_vecs)):
		if iorien%Blockdata["nproc"]!= Blockdata["main_node"]:
			if iorien%Blockdata["nproc"]==Blockdata["myid"]: wrap_mpi_send(ptls_in_orien_groups[iorien], Blockdata["main_node"], MPI_COMM_WORLD)
			if Blockdata["myid"] ==Blockdata["main_node"]: ptls_in_orien_groups[iorien] = wrap_mpi_recv(iorien%Blockdata["nproc"], MPI_COMM_WORLD)
		mpi_barrier(MPI_COMM_WORLD)	
	mpi_barrier(MPI_COMM_WORLD)
	if Blockdata["myid"] == Blockdata["main_node"]:print("step5")
	mpi_barrier(MPI_COMM_WORLD)
	zero_member_group_found = 0
	if Blockdata["myid"] == Blockdata["main_node"]:
		small_groups = []
		for iorien in xrange(len(refa_vecs)):
			if len(ptls_in_orien_groups[iorien]) <Tracker["min_orien_group_size"]:small_groups.append(iorien)
			if len(ptls_in_orien_groups[iorien]) == 0:
				zero_member_group_found += 1
		matched_pairs = find_neighborhood(refa_vecs, small_groups)
		ptls_in_orien_groups = reassign_ptls_in_orien_groups(ptls_in_orien_groups, matched_pairs)
	else: ptls_in_orien_groups = 0
	zero_member_group_found = bcast_number_to_all(zero_member_group_found, Blockdata["main_node"], MPI_COMM_WORLD)
	ptls_in_orien_groups = wrap_mpi_bcast(ptls_in_orien_groups, Blockdata["main_node"], MPI_COMM_WORLD)
	del refa_vecs, refa
	del local_orien_group_assignment
	del data_angles
	del orien_group_assignment
	return ptls_in_orien_groups

def compare_two_iterations(assignment1, assignment2, number_of_groups):
	# compare two assignments during clustering, either iteratively or independently
	import numpy as np
	assigned_groups1 =[[] for i in xrange(number_of_groups)]
	for im in xrange(len(assignment1)):assigned_groups1[assignment1[im]].append(im)
	res1 = []
	for iref in xrange(number_of_groups):
		a = np.array(assigned_groups1[iref],"int32")
		a.sort()
		res1.append(a)
	assigned_groups2 =[[] for i in xrange(number_of_groups)]
	for im in xrange(len(assignment2)): assigned_groups2[assignment2[im]].append(im)
	res2 = []
	for iref in xrange(number_of_groups):
		a = np.array(assigned_groups2[iref],"int32")
		a.sort()
		res2.append(a)
		del a
	newindeces, list_stable, nb_tot_objs = k_means_match_clusters_asg_new(res1, res2)
	del res1
	del res2
	return float(nb_tot_objs)/len(assignment1), newindeces, list_stable
	
def update_data_assignment(cdata, rdata, assignment, proc_list, nosmearing, myid):
	nima = len(cdata)
	groupids  = assignment[proc_list[myid][0]:proc_list[myid][1]]
	for im in xrange(nima):
		try: previous_group = cdata[im].get_attr("group")
		except: previous_group = -1
		cdata[im].set_attr("group", groupids[im])
		if nosmearing:
			rdata[im].set_attr("group", groupids[im])
			rdata[im].set_attr("previous_group", previous_group) 
		else:
			for jm in xrange(len(rdata[im])):
				rdata[im][jm].set_attr("previous_group", previous_group) 
				rdata[im][jm].set_attr("group", groupids[im])
	return
	
def update_rdata_assignment(assignment, proc_list, myid, rdata):
	nima = len(rdata)
	groupids = assignment[proc_list[myid][0]:proc_list[myid][1]]
	for im in xrange(nima): rdata[im].set_attr("group", groupids[im])
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
			if value == None and value_refinement != None: Tracker[key] = value_refinement
		except:
			if Blockdata["myid"] == Blockdata["main_node"]: print(key, " in sorting set as ", value, ", while in refinement, it is set as ", value_refinement)
	return
	
def print_dict(dict,theme, exclude = "refinement"):
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	print(line,theme)
	spaces = "                    "
	if exclude =="refinement": exclude = ["constants", "nodes", "yr", "shared_comm", "bckgnoise", "myid", "myid_on_node", "accumulatepw", "chunk_dict", "PW_dict", "full_list", "rshifts", "refang"]
	else: exclude = ["constants", "chunk_dict", "PW_dict", "full_list", "rshifts", "refang"]
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
##<<<-----------rec3d for sorting------->>>>>>>>>
def stepone(tvol, tweight):
	global Tracker, Blockdata
	tvol.set_attr("is_complex",1)
	ovol = Util.shrinkfvol(tvol,2)
	owol = Util.shrinkfvol(tweight,2)
	if( Tracker["constants"]["symmetry"] != "c1" ):
		ovol = ovol.symfvol(Tracker["constants"]["symmetry"], -1)
		owol = owol.symfvol(Tracker["constants"]["symmetry"], -1)
	return Util.divn_cbyr(ovol,owol)
	
def steptwo_mpi(tvol, tweight, treg, cfsc = None, regularized = True, color = 0):
	global Tracker, Blockdata
	if( Blockdata["color"] != color ):return model_blank(1)  #  This should not be executed if called properly
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
	we_data =  get_image_data(tweight)
	#  tvol is overwritten, meaning it is also an output
	
	ifi = mpi_iterefa( vol_data.__array_interface__['data'][0] ,  we_data.__array_interface__['data'][0] , nx, ny, nz, maxr2, \
			Tracker["constants"]["nnxo"], Blockdata["myid_on_node"], color, Blockdata["no_of_processes_per_group"],  Blockdata["shared_comm"])	
	if( Blockdata["myid_on_node"] == 0 ):
		#  Either pad or window in F space to 2*nnxo
		nx = tvol.get_ysize()
		if( nx > 2*Tracker["constants"]["nnxo"]): tvol = fdecimate(tvol, 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], False, False)
		elif(nx < 2*Tracker["constants"]["nnxo"]): tvol = fpol(tvol, 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], RetReal = False, normalize = False)
		tvol = fft(tvol)
		tvol = cyclic_shift(tvol,Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
		tvol = Util.window(tvol, Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
		tvol.div_sinc(1)
		tvol = cosinemask(tvol, Tracker["constants"]["nnxo"]//2-1,5, None)
		return tvol
	else:  return None
### non_mpi
def steptwo(tvol, tweight, treg, cfsc = None, regularized = True):
	global Tracker, Blockdata
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

	#  tvol is overwritten, meaning it is also an output
	Util.iterefa(tvol, tweight, maxr2, Tracker["constants"]["nnxo"])
	#  Either pad or window in F space to 2*nnxo
	nx = tvol.get_ysize()
	if( nx > 2*Tracker["constants"]["nnxo"] ):
		tvol = fdecimate(tvol, 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], False, False)
	elif(nx < 2*Tracker["constants"]["nnxo"]):
		tvol = fpol(tvol, 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], RetReal = False, normalize = False)
	tvol = fft(tvol)
	tvol = cyclic_shift(tvol,Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	tvol.set_attr("npad",2)
	tvol.div_sinc(1)
	tvol.del_attr("npad")
	tvol = Util.window(tvol, Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	tvol = cosinemask(tvol,Tracker["constants"]["nnxo"]//2-1,5, None)# clean artifacts in corners
	return tvol
####<<<<<-----------	
def recons3d_4nnsorting_MPI(myid, main_node, prjlist, random_subset, CTF = True, upweighted = True, mpi_comm= None, target_size=-1):
	"""
		recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
		Input
			list_of_prjlist: list of lists of projections to be included in the reconstruction
	"""
	from utilities		import reduce_EMData_to_root, random_string, get_im, findall, model_blank, info, get_params_proj
	from filter			import filt_table
	from reconstruction import insert_slices_pdf
	from fundamentals	import fft
	from statistics	    import fsc
	from EMAN2			import Reconstructors
	from mpi			import MPI_COMM_WORLD, mpi_barrier
	import types
	import datetime
	if mpi_comm == None: mpi_comm = MPI_COMM_WORLD
	imgsize = prjlist[0].get_ysize()  # It can be Fourier, so take y-size
	refvol = model_blank(target_size)
	refvol.set_attr("fudge", 1.0)
	if CTF: do_ctf = 1
	else: do_ctf = 0
	fftvol = EMData()
	weight = EMData()
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = Reconstructors.get( "nn4_ctfw", params )
	r.setup()
	#if norm_per_particle == None: norm_per_particle = len(prjlist)*[1.0]
	for im in xrange(len(prjlist)):
		phi, theta, psi, s2x, s2y = get_params_proj(prjlist[im], xform = "xform.projection") # shifts are already applied 
		if random_subset == 2:
			#try:	bckgn = prjlist[im].get_attr("bckgnoise")
			bckgn = target_size*[1.]
			if prjlist[im].get_attr("is_complex") == 0:	prjlist[im] = fft(prjlist[im]) 
			prjlist[im].set_attr_dict({"padffted":1, "is_complex":1})
			if not upweighted:  prjlist[im] = filt_table(prjlist[im], bckgn)
			prjlist[im].set_attr("bckgnoise", bckgn)
			r.insert_slice(prjlist[im], Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi}), 1.0)
		else:
			if prjlist[im].get_attr("chunk_id") == random_subset:
				#try:	bckgn = prjlist[im].get_attr("bckgnoise")
				bckgn = target_size*[1.]
				if prjlist[im].get_attr("is_complex")==0:	
					prjlist[im] = fft(prjlist[im])
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
	try:    qt = projlist[0].get_attr("qt")
	except: qt = 1.0
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = Reconstructors.get( "nn4_ctfw", params )
	r.setup()
	for im in xrange(len(prjlist)):
		phi, theta, psi, s2x, s2y = get_params_proj(prjlist[im], xform = "xform.projection") # shifts are already applied
		if prjlist[im].get_attr("group") == group_ID:
			if random_subset == 2:
				try:	bckgn = prjlist[im].get_attr("bckgnoise")
				except:	bckgn = target_size*[1.]
				if prjlist[im].get_attr("is_complex") == 0:	image = fft(prjlist[im])
				else: image =  prjlist[im].copy()
				image.set_attr_dict({"padffted":1, "is_complex":1})
				if not upweighted:  image = filt_table(image, bckgn)
				image.set_attr("bckgnoise", bckgn)
				r.insert_slice(image, Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi}), 1.0)
			else:
				if prjlist[im].get_attr("chunk_id") == random_subset:
					try:     bckgn = prjlist[im].get_attr("bckgnoise")
					except:	 bckgn = target_size*[1.]
					if prjlist[im].get_attr("is_complex")==0: image = fft(prjlist[im])
					else: image =  prjlist[im].copy()
					image.set_attr_dict({"padffted":1, "is_complex":1})
					if not upweighted:  image = filt_table(image, bckgn)              
					image.set_attr("bckgnoise", bckgn)
					r.insert_slice(image, Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi}), 1.0)	
	reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
	reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
	if myid == main_node: dummy = r.finish(True)
	mpi_barrier(mpi_comm)
	if myid == main_node: return fftvol, weight, refvol
	else: return None, None, None

def do3d_sorting(procid, data):
	global Tracker, Blockdata
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if(Blockdata["myid"] == Blockdata["main_node"]):print(line, "do3d_sorting")
	tvol, tweight, trol = recons3d_4nnsorting_MPI(myid = Blockdata["myid"], main_node = Blockdata["nodes"][procid], prjlist = data,\
		random_subset = procid, CTF = Tracker["constants"]["CTF"], upweighted = False, target_size = (2*Tracker["nxinit"]+3))
	if(Blockdata["myid"] == Blockdata["nodes"][procid]):
		if(procid == 0):
			cmd = "{} {}".format("mkdir", os.path.join(Tracker["directory"],"tempdir"))
			if os.path.exists(os.path.join(Tracker["directory"], "tempdir")):print("tempdir exists")
			else: cmdexecute(cmd)
		tvol.set_attr("is_complex",0)
		tvol.write_image(os.path.join(Tracker["directory"], "tempdir", "tvol_%01d.hdf"%procid))
		tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%01d.hdf"%procid))
		trol.write_image(os.path.join(Tracker["directory"], "tempdir", "trol_%01d.hdf"%procid))
	mpi_barrier(MPI_COMM_WORLD)
	return
	
def do3d_sorting_groups(particle_ID_index, partstack):
	global Tracker, Blockdata
	from utilities import get_im, wrap_mpi_bcast
	data = get_shrink_data_sorting(particle_ID_index, partstack)
	do3d_sorting_group_insertion(data)
	mpi_barrier(MPI_COMM_WORLD)
	fsc143                          =   0
	fsc05                           =   0
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
					tvol0 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0_%d.hdf")%index_of_group)
					tweight0 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_0_%d.hdf")%index_of_group)
					tvol1 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_1_%d.hdf")%index_of_group)
					tweight1 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_1_%d.hdf")%index_of_group)
					Util.fuse_low_freq(tvol0, tvol1, tweight0, tweight1, 2*Tracker["fuse_freq"])
					#print("   fuse done ")
					tag = 7007
					send_EMData(tvol1, Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					send_EMData(tweight1, Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					shrank0 	= stepone(tvol0, tweight0)
				elif Blockdata["myid_on_node"] == Blockdata["no_of_processes_per_group"]-1:				
					tag = 7007
					tvol1 		= recv_EMData(0, tag, Blockdata["shared_comm"])
					tweight1 	= recv_EMData(0, tag, Blockdata["shared_comm"])
					tvol1.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1} )
					shrank1 	= stepone(tvol1, tweight1)
				mpi_barrier(Blockdata["shared_comm"])
				if 	Blockdata["myid_on_node"] == 0:
					tag = 7007					
					send_EMData(shrank0, Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					del shrank0
					lcfsc = 0
				elif Blockdata["myid_on_node"] == Blockdata["no_of_processes_per_group"]-1:
					tag = 7007
					shrank0 = recv_EMData(0, tag, Blockdata["shared_comm"])
					cfsc 	= fsc(shrank0, shrank1)[1]
					write_text_row(cfsc, os.path.join(Tracker["directory"], "fsc_driver_%d.txt")%index_of_group)
					del shrank0, shrank1
					if(Tracker["nxinit"]<Tracker["constants"]["nnxo"]):
						cfsc = cfsc[:Tracker["nxinit"]]
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
				if( Blockdata["myid_on_node"] == 0): 
					res_05[index_of_group]  = Tracker["fsc05"]
					res_143[index_of_group] = Tracker["fsc143"]
				if Blockdata["fftwmpi"]: 
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
						if(Tracker["mask3D"] == None):  tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
						else:  Util.mul_img(tvol2, get_im(Tracker["mask3D"]))
						tvol2.write_image(os.path.join(Tracker["directory"], "vol_grp%03d_iter000.hdf"%(index_of_group)))
						del tvol2
				else:
					if( Blockdata["myid_on_node"] == 0):
						tvol2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")%index_of_group)
						tweight2 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")%index_of_group)
						treg2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_group))
						tvol2       = steptwo(tvol2, tweight2, treg2, cfsc, True)
						del tweight2, treg2
						if(Tracker["mask3D"] == None): tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
						else: Util.mul_img(tvol2, get_im(Tracker["mask3D"]))
						tvol2.write_image(os.path.join(Tracker["directory"], "vol_grp%03d_iter000.hdf"%(index_of_group)))
						del tvol2
				mpi_barrier(Blockdata["shared_comm"])
			mpi_barrier(Blockdata["shared_comm"])	
	mpi_barrier(MPI_COMM_WORLD)
	res_05  = mpi_reduce(res_05,  Tracker["number_of_groups"], MPI_INT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	res_143 = mpi_reduce(res_143, Tracker["number_of_groups"], MPI_INT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	res_05  = map(int, res_05)
	res_143 = map(int, res_143)
	if (Blockdata["myid"] == Blockdata["main_node"]):
		Tracker["fsc143"] = res_143
		Tracker["fsc05"]  = res_05
	Tracker = wrap_mpi_bcast(Tracker,Blockdata["main_node"])
	return 
			
def do3d_sorting_group_insertion(data, randomset=2):
	global Tracker, Blockdata
	if(Blockdata["myid"] == Blockdata["node_volume"][0]):
		cmd = "{} {}".format("mkdir", os.path.join(Tracker["directory"], "tempdir"))
		if os.path.exists(os.path.join(Tracker["directory"], "tempdir")):print("tempdir exists")
		else: cmdexecute(cmd)
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] ==  Blockdata["main_node"]:print(line, "start backprojection of %d volumes"%Tracker["number_of_groups"])
	if randomset ==1:
		for index_of_groups in xrange(Tracker["number_of_groups"]):
			for procid in xrange(2, 3):
				tvol, tweight, trol = recons3d_4nnsorting_group_MPI(myid = Blockdata["myid"], main_node = Blockdata["nodes"][0],\
				  prjlist = data,  random_subset = procid, group_ID = index_of_groups, CTF = Tracker["constants"]["CTF"],\
					upweighted = False, target_size = (2*Tracker["nxinit"]+3))
				
				if(Blockdata["myid"] == Blockdata["nodes"][procid]):
					tvol.set_attr("is_complex",0)
					tvol.write_image(os.path.join(Tracker["directory"], "tempdir", "tvol_%d_%d.hdf"%(procid, index_of_groups)))
					tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%d_%d.hdf"%(procid, index_of_groups)))
					trol.write_image(os.path.join(Tracker["directory"], "tempdir", "trol_%d_%d.hdf"%(procid, index_of_groups)))
				mpi_barrier(MPI_COMM_WORLD)
	else:		
		for index_of_groups in xrange(Tracker["number_of_groups"]):
			for procid in xrange(2):
				tvol, tweight, trol = recons3d_4nnsorting_group_MPI(myid = Blockdata["myid"], main_node = Blockdata["nodes"][procid],  \
				  prjlist = data, random_subset = procid, group_ID = index_of_groups, CTF = Tracker["constants"]["CTF"],\
					upweighted = False, target_size = (2*Tracker["nxinit"]+3))
				
				if(Blockdata["myid"] == Blockdata["nodes"][procid]):
					tvol.set_attr("is_complex",0)
					tvol.write_image(os.path.join(Tracker["directory"], "tempdir", "tvol_%d_%d.hdf"%(procid, index_of_groups)))
					tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%d_%d.hdf"%(procid, index_of_groups)))
					trol.write_image(os.path.join(Tracker["directory"], "tempdir", "trol_%d_%d.hdf"%(procid, index_of_groups)))
				mpi_barrier(MPI_COMM_WORLD)
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] ==  Blockdata["main_node"]:print(line, "end backprojection of %d volumes"%Tracker["number_of_groups"])		
	mpi_barrier(MPI_COMM_WORLD)
	return
	
def do3d_sorting_groups_trl_iter(data, iteration):
	global Tracker, Blockdata
	from utilities import get_im, write_text_row, bcast_number_to_all, wrap_mpi_bcast
	keepgoing = 1
	if(Blockdata["myid"] == Blockdata["nodes"][0]):
		cmd = "{} {}".format("mkdir", os.path.join(Tracker["directory"], "tempdir"))
		if os.path.exists(os.path.join(Tracker["directory"], "tempdir")): print("tempdir exists")
		else: cmdexecute(cmd)
	do3d_sorting_group_insertion(data)
	mpi_barrier(MPI_COMM_WORLD)
	fsc143                          =   0
	fsc05                           =   0
	Tracker["fsc143"]				=	0
	Tracker["fsc05"]				=	0
	res_05 						    =	Tracker["number_of_groups"]*[0]
	res_143 					    =	Tracker["number_of_groups"]*[0]
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] == Blockdata["main_node"]: print(line, "trl_iter starts to do volumes, number of groups is ", Tracker["number_of_groups"])
	for index_of_colors in xrange(Blockdata["no_of_groups"]):
		group_start, group_end = MPI_volume_start_end(Tracker["number_of_groups"], Blockdata["no_of_groups"], index_of_colors)
		if Blockdata["color"] == index_of_colors:  # It has to be 1 to avoid problem with tvol1 not closed on the disk
			for index_of_group in xrange(group_start, group_end):
				cfsc = 0
				if Blockdata["myid_on_node"] == 0:
					tvol0 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0_%d.hdf")%index_of_group)
					tweight0 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_0_%d.hdf")%index_of_group)
					tvol1 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_1_%d.hdf")%index_of_group)
					tweight1 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_1_%d.hdf")%index_of_group)
					Util.fuse_low_freq(tvol0, tvol1, tweight0, tweight1, 2*Tracker["fuse_freq"])
					tag = 7007
					send_EMData(tvol1, Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					send_EMData(tweight1, Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					shrank0 	= stepone(tvol0, tweight0)
				elif Blockdata["myid_on_node"] == Blockdata["no_of_processes_per_group"]-1:					
					tag = 7007
					tvol1 		= recv_EMData(0, tag, Blockdata["shared_comm"])
					tweight1 	= recv_EMData(0, tag, Blockdata["shared_comm"])
					tvol1.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1} )
					shrank1 	= stepone(tvol1, tweight1)
				mpi_barrier(Blockdata["shared_comm"])
				if 	Blockdata["myid_on_node"] == 0:
					tag = 7007					
					send_EMData(shrank0, Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					del shrank0
					lcfsc = 0
				elif Blockdata["myid_on_node"] == Blockdata["no_of_processes_per_group"]-1:
					tag = 7007
					shrank0 	= recv_EMData(0, tag, Blockdata["shared_comm"])
					#  Note shrank volumes are Fourier uncentered.
					cfsc 		= fsc(shrank0, shrank1)[1]
					write_text_row(cfsc, os.path.join(Tracker["directory"], "fsc_driver_grp%03d_iter%03d.txt")%(index_of_group,iteration))
					del shrank0, shrank1
					if(Tracker["nxinit"]<Tracker["constants"]["nnxo"]):
						cfsc 	= cfsc[:Tracker["nxinit"]]
						for i in xrange(len(cfsc),Tracker["constants"]["nnxo"]//2+1): cfsc.append(0.0)
					lcfsc = len(cfsc)							
					fsc05  = 0
					fsc143 = 0 
					for ifreq in xrange(len(cfsc)):	
						if cfsc[ifreq] <0.5: break
					fsc05  = ifreq - 1
					for ifreq in xrange(len(cfsc)):
						if cfsc[ifreq]<0.143: break
					fsc143 = ifreq - 1
					Tracker["fsc143"] = fsc143
					Tracker["fsc05"]  = fsc05
					line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
					print(line, "group %d  of do3d_sorting_groups_trl_iter is done"%index_of_group)
				Tracker = wrap_mpi_bcast(Tracker, Blockdata["no_of_processes_per_group"]-1, Blockdata["shared_comm"])
				cfsc = wrap_mpi_bcast(cfsc, Blockdata["no_of_processes_per_group"]-1, Blockdata["shared_comm"])
				######
				Tracker["maxfrad"] = Tracker["nxinit"]//2
				if( Blockdata["myid_on_node"] == 0):
					res_05[index_of_group]  = Tracker["fsc05"]
					res_143[index_of_group] = Tracker["fsc143"]
				if Blockdata["fftwmpi"]:
					if( Blockdata["myid_on_node"] == 0):
						tvol2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0_%d.hdf")%index_of_group)
						tweight2 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_0_%d.hdf")%index_of_group)
						treg2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "trol_0_%d.hdf"%index_of_group))
					else:
						tvol2 		= model_blank(1)
						tweight2 	= model_blank(1)
						treg2		= model_blank(1)
					tvol2 = steptwo_mpi(tvol2, tweight2, treg2, cfsc, False, color = index_of_colors) # has to be False!!!
					del tweight2, treg2
					if( Blockdata["myid_on_node"] == 0):
						if(Tracker["mask3D"] == None): tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
						else: Util.mul_img(tvol2, get_im(Tracker["constants"]["mask3D"]))
						tvol2.write_image(os.path.join(Tracker["directory"], "vol_unfiltered_0_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
						del tvol2
				else:
					if( Blockdata["myid_on_node"] == 0):
						tvol2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0_%d.hdf")%index_of_group)
						tweight2 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_0_%d.hdf")%index_of_group)
						treg2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "trol_0_%d.hdf"%index_of_group))
						tvol2       = steptwo(tvol2, tweight2, treg2, cfsc, False)
						del tweight2, treg2
						if(Tracker["mask3D"] == None):  tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
						else:  Util.mul_img(tvol2, get_im(Tracker["constants"]["mask3D"]))
						tvol2.write_image(os.path.join(Tracker["directory"], "vol_unfiltered_0_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
						del tvol2
				######
				mpi_barrier(Blockdata["shared_comm"])
				if( Blockdata["myid_on_node"] == 0):
					res_05[index_of_group]  = Tracker["fsc05"]
					res_143[index_of_group] = Tracker["fsc143"]
				if Blockdata["fftwmpi"]:
					if( Blockdata["myid_on_node"] == 0):
						tvol2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_1_%d.hdf")%index_of_group)
						tweight2 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_1_%d.hdf")%index_of_group)
						treg2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "trol_1_%d.hdf"%index_of_group))
					else:
						tvol2 		= model_blank(1)
						tweight2 	= model_blank(1)
						treg2		= model_blank(1)
					tvol2 = steptwo_mpi(tvol2, tweight2, treg2, cfsc, False, color = index_of_colors) # has to be False!!!
					del tweight2, treg2
					if( Blockdata["myid_on_node"] == 0):
						if(Tracker["mask3D"] == None): tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
						else: Util.mul_img(tvol2, get_im(Tracker["constants"]["mask3D"]))
						tvol2.write_image(os.path.join(Tracker["directory"], "vol_unfiltered_1_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
						del tvol2
				else:
					if( Blockdata["myid_on_node"] == 0):
						tvol2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_1_%d.hdf")%index_of_group)
						tweight2 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_1_%d.hdf")%index_of_group)
						treg2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "trol_1_%d.hdf"%index_of_group))
						tvol2 = steptwo(tvol2, tweight2, treg2, cfsc, False)
						del tweight2, treg2
						if(Tracker["mask3D"] == None):  tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
						else:  Util.mul_img(tvol2, get_im(Tracker["constants"]["mask3D"]))
						tvol2.write_image(os.path.join(Tracker["directory"], "vol_unfiltered_1_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
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
	keepgoing = bcast_number_to_all(keepgoing, source_node = Blockdata["main_node"], mpi_comm = MPI_COMM_WORLD) # always check 
	Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"])
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if (Blockdata["myid"] == Blockdata["main_node"]):print(line, " trl_iter do volumes end")
	if not keepgoing: ERROR("do3d_sorting_groups_trl_iter  %s"%os.path.join(Tracker["directory"], "tempdir"),"do3d_sorting_groups_trl_iter", 1, Blockdata["myid"]) 
	return
	
# Three ways of importing refinement results
def get_input_from_sparx_ref3d(log_main):# case one
	# import SPARX results
	global Tracker, Blockdata
	import json
	from   string import split, atoi
	import_from_sparx_refinement = 1
	selected_iter      = 0
	Tracker_refinement = 0
	if not os.path.exists (Tracker["constants"]["refinement_dir"]): ERROR("SPARX refinement dir does not exist", "get_input_from_sparx_ref3d", 1, Blockdata["myid"])
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] == Blockdata["main_node"]:
		msg = "Import results from SPARX 3-D refinement"
		print(line, msg)
		log_main.add(msg)
		if Tracker["constants"]["niter_for_sorting"] == -1: # take the best solution to do sorting
			msg = "Search in the directory %s ......"%Tracker["constants"]["refinement_dir"]
			print(line, msg)
			log_main.add(msg)
			niter_refinement = 0
			while os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%niter_refinement)) and os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"],"main%03d"%niter_refinement, "Tracker_%03d.json"%niter_refinement)):
				niter_refinement +=1
			niter_refinement -=1
			if niter_refinement !=0:
				fout = open(os.path.join(Tracker["constants"]["refinement_dir"],"main%03d"%niter_refinement, "Tracker_%03d.json"%niter_refinement),'r')
				Tracker_refinement 	= convert_json_fromunicode(json.load(fout))
				fout.close()
				selected_iter = Tracker_refinement["constants"]["best"]
			else:		import_from_sparx_refinement = 0
		else:
			msg = "Try to load json file ...%s"%os.path.join(Tracker["constants"]["refinement_dir"],"main%03d"%Tracker["constants"]["niter_for_sorting"],\
			 "Tracker_%03d.json"%Tracker["constants"]["niter_for_sorting"])
			print(line, msg)
			log_main.add(msg)
			try:
				fout = open(os.path.join(Tracker["constants"]["refinement_dir"],"main%03d"%Tracker["constants"]["niter_for_sorting"], \
				"Tracker_%03d.json"%Tracker["constants"]["niter_for_sorting"]),'r')
				Tracker_refinement	= convert_json_fromunicode(json.load(fout))
				fout.close()
				selected_iter = Tracker["constants"]["niter_for_sorting"]
			except:	import_from_sparx_refinement = 0
	else: selected_iter = -1	
	selected_iter = bcast_number_to_all(selected_iter, Blockdata["main_node"], MPI_COMM_WORLD)
	import_from_sparx_refinement = bcast_number_to_all(import_from_sparx_refinement, source_node = Blockdata["main_node"])
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
		if Tracker_refinement["constants"]["stack"][0:4]=="bdb:": refinement_stack = "bdb:"+os.path.join(refinement_dir_path, Tracker_refinement["constants"]["stack"][4:])
		else: refinement_stack = os.path.join(refinement_dir_path, Tracker_refinement["constants"]["stack"])
		if not Tracker["constants"]["orgstack"]: # Use refinement stack if instack is not provided
			msg = "refinement stack  %s"%refinement_stack
			print(line, msg)
			log_main.add(msg)
			Tracker["constants"]["orgstack"] = refinement_stack #Tracker_refinement["constants"]["stack"]
			print(line, "The refinement image stack is %s"%Tracker_refinement["constants"]["stack"])
			try: image = get_im(Tracker["constants"]["orgstack"], 0)
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
					old_stack = Tracker["constants"]["stack"]
					if old_stack[0:3] == "bdb":
						Tracker["constants"]["orgstack"] = "bdb:" + Tracker["constants"]["refinement_dir"]+"/../"+old_stack[4:]
					else: Tracker["constants"]["orgstack"] = os.path.join(option_old_refinement_dir, "../", old_stack)
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
	
	if import_from_sparx_refinement == 0:ERROR("The data stack is not accessible","get_input_from_sparx_ref3d",1, Blockdata["myid"])
	total_stack = bcast_number_to_all(total_stack, source_node = Blockdata["main_node"])			
	Tracker["constants"]["total_stack"] = total_stack
	# Now copy relevant refinement files to sorting directory:
	if Blockdata["myid"] == Blockdata["main_node"]:
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iter, "params_%03d.txt"%selected_iter)):
			cmd = "{} {} {}".format("cp ",os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iter, \
			 "params_%03d.txt"%selected_iter), os.path.join(Tracker["constants"]["masterdir"], "sparx_refinement_params.txt"))
			cmdexecute(cmd)
		else: import_from_sparx_refinement = 0
		Tracker["constants"]["selected_iter"] = selected_iter
	import_from_sparx_refinement = bcast_number_to_all(import_from_sparx_refinement, source_node = Blockdata["main_node"])
	if import_from_sparx_refinement == 0:ERROR("The parameter file of the best solution is not accessible", "get_input_from_sparx_ref3d", 1, Blockdata["myid"])
		
	if Blockdata["myid"] == Blockdata["main_node"]:
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iter, "bckgnoise.hdf")):
			cmd = "{} {} {}".format("cp ",os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iter, "bckgnoise.hdf"),\
			os.path.join(Tracker["constants"]["masterdir"], "bckgnoise.hdf"))
			cmdexecute(cmd)
		else:
			import_from_sparx_refinement == 0
			for search_iter in xrange(selected_iter-1, 0, -1):
				 if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%search_iter, "bckgnoise.hdf")):
					cmd = "{} {} {}".format("cp ",os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%search_iter, \
					"bckgnoise.hdf"), os.path.join(Tracker["constants"]["masterdir"], "bckgnoise.hdf"))
					cmdexecute(cmd)
					import_from_sparx_refinement = 1
					break
	import_from_sparx_refinement = bcast_number_to_all(import_from_sparx_refinement, source_node = Blockdata["main_node"])
	
	if import_from_sparx_refinement == 0:
		Tracker["bckgnoise"] = None
		if Blockdata["myid"] == Blockdata["main_node"]:	print("Noise file is not found. However we continue")
	else: Tracker["bckgnoise"] = os.path.join(Tracker["constants"]["masterdir"], "bckgnoise.hdf")
	
	import_from_sparx_refinement = 1	
	if Blockdata["myid"] == Blockdata["main_node"]:
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iter, "driver_%03d.txt"%selected_iter)):
			cmd = "{} {} {}".format("cp ",os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iter, \
			 "driver_%03d.txt"%selected_iter), os.path.join(Tracker["constants"]["masterdir"], "fsc_global.txt"))
			cmdexecute(cmd)
		else: import_from_sparx_refinement = 0
		#Tracker["constants"]["selected_iter"] = selected_iter
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
	if import_from_sparx_refinement == 0:ERROR("The driver of the best solution is not accessible","get_input_from_sparx_ref3d", 1, Blockdata["myid"])
	if Blockdata["myid"] == Blockdata["main_node"]:
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main000/indexes_000.txt")):
			cmd = "{} {} {}".format("cp ", os.path.join(Tracker["constants"]["refinement_dir"], "main000/indexes_000.txt"), \
			os.path.join(Tracker["constants"]["masterdir"], "indexes.txt"))
			cmdexecute(cmd)
		else:	import_from_sparx_refinement = 0	
	import_from_sparx_refinement = bcast_number_to_all(import_from_sparx_refinement, source_node = Blockdata["main_node"])
	if import_from_sparx_refinement == 0: ERROR("The index file of the best solution are not accessible","get_input_from_sparx_ref3d", 1, Blockdata["myid"])
	if Blockdata["myid"] == Blockdata["main_node"]:
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main000/chunk_0_000.txt")):
			cmd = "{} {} {}".format("cp ", os.path.join(Tracker["constants"]["refinement_dir"], "main000/chunk_0_000.txt"), \
			os.path.join(Tracker["constants"]["masterdir"], "chunk_0.txt"))
			cmdexecute(cmd)
		else: import_from_sparx_refinement == 0
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main000/chunk_1_000.txt")):
			cmd = "{} {} {}".format("cp ", os.path.join(Tracker["constants"]["refinement_dir"], "main000/chunk_1_000.txt"), \
			os.path.join(Tracker["constants"]["masterdir"], "chunk_1.txt"))
			cmdexecute(cmd)
		else: import_from_sparx_refinement == 0
			
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main000/particle_groups_0.txt")):
			cmd = "{} {} {}".format("cp ", os.path.join(Tracker["constants"]["refinement_dir"], "main000/particle_groups_0.txt"), \
			os.path.join(Tracker["constants"]["masterdir"], "particle_groups_0.txt"))
			cmdexecute(cmd)
		else: import_from_sparx_refinement == 0
			
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main000/particle_groups_1.txt")):
			cmd = "{} {} {}".format("cp ", os.path.join(Tracker["constants"]["refinement_dir"], "main000/particle_groups_1.txt"), \
			os.path.join(Tracker["constants"]["masterdir"], "particle_groups_1.txt"))
			cmdexecute(cmd)
		else:import_from_sparx_refinement == 0
	import_from_sparx_refinement = bcast_number_to_all(import_from_sparx_refinement, source_node = Blockdata["main_node"])
	if import_from_sparx_refinement == 0:ERROR("The chunk files and partice group files are not accessible","get_input_from_sparx_ref3d",1, Blockdata["myid"])
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	# copy all relavant parameters into sorting tracker
	if Blockdata["myid"] == Blockdata["main_node"]:
		if Tracker["constants"]["radius"] == -1: Tracker["constants"]["radius"] = Tracker_refinement["constants"]["radius"]
		Tracker["constants"]["nnxo"]       = Tracker_refinement["constants"]["nnxo"]
		Tracker["constants"]["orgres"]     = Tracker_refinement["bestres"]
		Tracker["delta"]                   = Tracker_refinement["delta"]
		Tracker["ts"]                      = Tracker_refinement["ts"]
		Tracker["xr"]                      = Tracker_refinement["xr"]
		Tracker["constants"]["pixel_size"] = Tracker_refinement["constants"]["pixel_size"]
		Tracker["avgnorm"]                 = Tracker_refinement["avgvaradj"]
		if Tracker["constants"]["nxinit"]<0: Tracker["nxinit_refinement"] = Tracker_refinement["nxinit"] #Sphire window size
		else:  Tracker["nxinit_refinement"] = Tracker["constants"]["nxinit"] #User defined window size
		 
		try:     sym =  Tracker_refinement["constants"]["sym"]
		except:  sym =  Tracker_refinement["constants"]["symmetry"]
		Tracker["constants"]["symmetry"] = sym
		print(line, "Parameters importing is done!")
		if not Tracker["constants"]["mask3D"] and Tracker_refinement["constants"]["mask3D"]:
			refinement_mask3D_path, refinement_mask3D_file = os.path.split(Tracker_refinement["constants"]["mask3D"])# MRK_DEBUG
			cmd = "{} {} {}".format("cp ", os.path.join(refinement_dir_path, Tracker_refinement["constants"]["mask3D"]), \
			os.path.join(Tracker["constants"]["masterdir"], refinement_mask3D_file))
			cmdexecute(cmd)
			Tracker["constants"]["mask3D"] = os.path.join(Tracker["constants"]["masterdir"], refinement_mask3D_file)
	Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], communicator = MPI_COMM_WORLD)
	import_from_sparx_refinement = bcast_number_to_all(import_from_sparx_refinement, source_node = Blockdata["main_node"])
	if not import_from_sparx_refinement:ERROR("Import parameters from SPARX refinement failed", "get_input_from_sparx_ref3d", 1,  Blockdata["myid"])
	
	# Setting for margin error				
	chunk_dict = {}
	group_dict = {}
	if(Blockdata["myid"] == Blockdata["main_node"]):
		chunk_one = read_text_file(os.path.join(Tracker["constants"]["masterdir"],"chunk_0.txt"))
		chunk_two = read_text_file(os.path.join(Tracker["constants"]["masterdir"],"chunk_1.txt"))
	else:
		chunk_one = 0
		chunk_two = 0
	chunk_one = wrap_mpi_bcast(chunk_one, Blockdata["main_node"])
	chunk_two = wrap_mpi_bcast(chunk_two, Blockdata["main_node"])
	#
	if(Blockdata["myid"] == Blockdata["main_node"]):
		chunk_one_group = read_text_file(os.path.join(Tracker["constants"]["masterdir"],"particle_groups_0.txt"))
		chunk_two_group = read_text_file(os.path.join(Tracker["constants"]["masterdir"],"particle_groups_1.txt"))
	else:
		chunk_one_group = 0
		chunk_two_group = 0
	chunk_one_group = wrap_mpi_bcast(chunk_one_group, Blockdata["main_node"])
	chunk_two_group = wrap_mpi_bcast(chunk_two_group, Blockdata["main_node"])
	for index_of_element in xrange(len(chunk_one)): 
		chunk_dict[chunk_one[index_of_element]] = 0
		group_dict[chunk_one[index_of_element]] = chunk_one_group[index_of_element]
	for index_of_element in xrange(len(chunk_two)): 
		chunk_dict[chunk_two[index_of_element]] = 1
		group_dict[chunk_two[index_of_element]] = chunk_two_group[index_of_element] 			
	Tracker["chunk_dict"] 	= chunk_dict
	Tracker["P_chunk_0"]   	= len(chunk_one)/float(total_stack)
	Tracker["P_chunk_1"]   	= len(chunk_two)/float(total_stack)
	if(Blockdata["myid"] == Blockdata["main_node"]):
		chunk_ids = []
		group_ids = []
		partids   = read_text_file(os.path.join(Tracker["constants"]["masterdir"], "indexes.txt"),-1)
		partids   = partids[0]
		Tracker["constants"]["total_stack"] = len(partids)
		params = read_text_file(os.path.join(Tracker["constants"]["masterdir"], "sparx_refinement_params.txt"),-1)
		for index_of_particle in xrange(len(partids)): 
			chunk_ids.append(chunk_dict[partids[index_of_particle]])
			group_ids.append(group_dict[partids[index_of_particle]])
		refinement_params = [ params[0], params[1], params[2], params[3], params[4], chunk_ids, group_ids, params[7]]
		write_text_file(refinement_params, os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt"))		
		line  = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print(line, "Initialization of sorting from SPARX refinement is done")
	else: Tracker["constants"]["total_stack"] = 0
	Tracker["constants"]["total_stack"]     = bcast_number_to_all(Tracker["constants"]["total_stack"], Blockdata["main_node"], MPI_COMM_WORLD)
	Tracker["total_stack"]                  = Tracker["constants"]["total_stack"]
	Tracker["constants"]["partstack"]	    = os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt")
	total_stack                             = Tracker["constants"]["total_stack"]
	Tracker["currentres"]                   = float(Tracker["constants"]["fsc05"])/float(Tracker["constants"]["nxinit"])
	Tracker["bckgnoise"]                    =  os.path.join(Tracker["constants"]["masterdir"], "bckgnoise.hdf")
	# Now copy oldparamstruture
	copy_oldparamstructure_from_meridien_MPI(selected_iter, log_main)
	return import_from_sparx_refinement
		
def get_input_from_datastack(log_main):# Case three
	global Tracker, Blockdata
	from utilities import write_text_file, write_text_row, wrap_mpi_bcast
	import json
	from   string import split, atoi
	from   random import shuffle
	import_from_data_stack = 1
	#print("fuse", Tracker["fuse_freq"])
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if(Blockdata["myid"] == Blockdata["main_node"]):
		msg =  "Sorting starts from given data stack %s with both xform.projection and ctf parameters filled in the headers"%Tracker["constants"]["orgstack"]
		print(line, msg)
		log_main.add(msg)
		image = get_im(Tracker["constants"]["orgstack"])
		Tracker["constants"]["nnxo"] = image.get_xsize()		
		if( Tracker["nxinit"] > Tracker["constants"]["nnxo"]):
				ERROR("Image size less than minimum permitted $d"%Tracker["nxinit"],"get_input_from_datastack",1, Blockdata["myid"])
				nnxo = -1
		else:
			if Tracker["constants"]["CTF"]:
				ictf = image.get_attr('ctf')
				Tracker["constants"]["pixel_size"] = ictf.apix
			else:
				Tracker["constants"]["pixel_size"] = 1.0
				del image
	else:
		Tracker["constants"]["nnxo"] = 0
		Tracker["constants"]["pixel_size"] = 1.0
	Tracker["constants"]["nnxo"] = bcast_number_to_all(Tracker["constants"]["nnxo"], source_node = Blockdata["main_node"])
	if( Tracker["constants"]["nnxo"] < 0): ERROR("Image size is negative", "get_input_from_datastack", 1, Blockdata["main_node"])
	Tracker["constants"]["pixel_size"]	= bcast_number_to_all(Tracker["constants"]["pixel_size"], source_node =  Blockdata["main_node"])
	if(Tracker["constants"]["radius"] < 1): Tracker["constants"]["radius"]  = Tracker["constants"]["nnxo"]//2-2
	elif((2*Tracker["constants"]["radius"] +2) > Tracker["constants"]["nnxo"]): ERROR("Particle radius set too large!", \
	"get_input_from_datastack",1, Blockdata["myid"])
	if Blockdata["myid"] == Blockdata["main_node"]:	total_stack = EMUtil.get_image_count(Tracker["constants"]["orgstack"])
	else: total_stack = 0
	total_stack = bcast_number_to_all(total_stack, Blockdata["main_node"])
	# randomly assign two subsets
	Tracker["constants"]["total_stack"]	= total_stack
	Tracker["constants"]["chunk_0"]		= os.path.join(Tracker["constants"]["masterdir"],"chunk_0.txt")
	Tracker["constants"]["chunk_1"]		= os.path.join(Tracker["constants"]["masterdir"],"chunk_1.txt")
	Tracker["constants"]["partstack"]	= os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt")
	Tracker["previous_parstack"]        = os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt")#
	###
	Tracker["refang"], Tracker["rshifts"], Tracker["delta"] = None, None, None
	Tracker["avgnorm"] =1.0
	chunk_dict = {}
	chunk_list = []
	if Blockdata["myid"] == Blockdata["main_node"]:
		chunk_dict  = {}
		tlist = range(total_stack)
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
	
	# Reconstruction to determine the resolution in orignal data size
	Tracker["nxinit"]     = Tracker["constants"]["nnxo"]
	Tracker["shrinkage"]  = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	Tracker["bckgnoise"]  =  None
	temp = model_blank(Tracker["constants"]["nnxo"], Tracker["constants"]["nnxo"])
	nny  =  temp.get_ysize()
	Blockdata["bckgnoise"] =  [1.0]*nny # set for initial recon3D of data from stack	
	Tracker["focus3D"]     =  None
	Tracker["fuse_freq"] = int(Tracker["constants"]["pixel_size"]*Tracker["constants"]["nnxo"]/Tracker["constants"]["fuse_freq"] +0.5)
	Tracker["directory"] = Tracker["constants"]["masterdir"]
	if Tracker["constants"]["nxinit"]<0: Tracker["nxinit_refinement"] = Tracker["constants"]["nnxo"]
	else: Tracker["nxinit_refinement"] = Tracker["constants"]["nxinit"]
	if Blockdata["myid"] == Blockdata["main_node"]:
		msg  = "reconstruct 3-D volumes from two random subsets to calculate FSC"
		print(line, msg)
		log_main.add(msg)
	for procid in xrange(2):
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		if Blockdata["myid"] == Blockdata["main_node"]: print(line, "Reconstruction of random subset %d"%procid)
		data = get_shrink_data_sorting(os.path.join(Tracker["constants"]["masterdir"],"chunk_%01d.txt"%procid), Tracker["constants"]["partstack"])
		mpi_barrier(MPI_COMM_WORLD)
		do3d_sorting(procid, data)
	mpi_barrier(MPI_COMM_WORLD)
	if( Blockdata["myid"] == Blockdata["main_shared_nodes"][1]):  # It has to be 1 to avoid problem with tvol1 not closed on the disk
		print("come to mainshared_nodes 1", Blockdata["myid"])
		print(Blockdata["main_shared_nodes"], Blockdata["no_of_processes_per_group"])
		tvol0 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0.hdf"))
		tweight0 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_0.hdf"))
		tvol1 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_1.hdf"))
		tweight1 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_1.hdf"))
		Util.fuse_low_freq(tvol0, tvol1, tweight0, tweight1, 2*Tracker["fuse_freq"])
		tag = 7007
		send_EMData(tvol1, Blockdata["main_shared_nodes"][0], tag, MPI_COMM_WORLD)
		send_EMData(tweight1, Blockdata["main_shared_nodes"][0], tag, MPI_COMM_WORLD)
		shrank0 	= stepone(tvol0, tweight0)
		send_EMData(shrank0, Blockdata["main_shared_nodes"][0], tag, MPI_COMM_WORLD)
		del shrank0
		lcfsc = 0
	elif( Blockdata["myid"] == Blockdata["main_shared_nodes"][0]):
		print("come to mainshared_nodes 1", Blockdata["myid"])
		tag = 7007
		tvol1 		= recv_EMData(Blockdata["main_shared_nodes"][1], tag, MPI_COMM_WORLD)
		tweight1 	= recv_EMData(Blockdata["main_shared_nodes"][1], tag, MPI_COMM_WORLD)
		tvol1.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1} )
		shrank1 	= stepone(tvol1, tweight1)
		#  Get shrank volume, do fsc, send it to all
		shrank0 	= recv_EMData(Blockdata["main_shared_nodes"][1], tag, MPI_COMM_WORLD)
		#  Note shrank volumes are Fourier uncentered.
		cfsc 		= fsc(shrank0, shrank1)[1]
		write_text_row(cfsc, os.path.join(Tracker["directory"], "fsc_global.txt"))
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
	else:lcfsc = 0
	Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_shared_nodes"][0], communicator = MPI_COMM_WORLD)
	#from mpi import mpi_finalize
	#mpi_finalize()
	#exit()
	return import_from_data_stack
	
#### Old rec3ds
# do sparx final reconstruction
def do_ctrefromsort3d_get_subset_data(masterdir, option_old_refinement_dir, option_selected_cluster, option_selected_iter, shell_line_command, comm=-1):
	global Tracker, Blockdata
	from utilities    import get_im, read_text_row, read_text_file, wrap_mpi_bcast, bcast_number_to_all, write_text_row, wrap_mpi_recv
	from applications import MPI_start_end
	import json
	if (Blockdata["subgroup_myid"] > -1):
		orgstack = Tracker["constants"]["orgstack"]
		selected_iter = option_selected_iter
		if Blockdata["subgroup_myid"] == Blockdata["main_node"]: cluster = sorted(read_text_file(option_selected_cluster))
		else: cluster = 0
		cluster = wrap_mpi_bcast(cluster, Blockdata["main_node"], comm) # balance processors	
		old_refinement_iter_dir    = os.path.join(option_old_refinement_dir, "main%03d"%selected_iter)
		old_oldparamstructure_dir  = os.path.join(old_refinement_iter_dir, "oldparamstructure")
		old_previousoutputdir      = os.path.join(option_old_refinement_dir, "main%03d"%(selected_iter-1))	
		if Blockdata["subgroup_myid"] == Blockdata["main_node"]: 
			nproc_old_ref3d = 0
			while os.path.exists(os.path.join(old_oldparamstructure_dir, "oldparamstructure_0_%03d_%03d.json"%(nproc_old_ref3d, selected_iter))):nproc_old_ref3d += 1
		else: nproc_old_ref3d = 0
		nproc_old_ref3d = bcast_number_to_all(nproc_old_ref3d, Blockdata["main_node"], comm)
	
		# read old refinement Tracker
		if Blockdata["subgroup_myid"] == Blockdata["main_node"]:
			fout = open(os.path.join(old_refinement_iter_dir, "Tracker_%03d.json"%selected_iter),"r")
			Tracker 	= convert_json_fromunicode(json.load(fout))
			fout.close()
		else: Tracker = 0
		Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], comm) # balance processors
		Tracker["constants"]["orgstack"] = orgstack
		print("myid", Blockdata["subgroup_myid"], Blockdata["myid"])
		
		if Blockdata["subgroup_myid"] == Blockdata["main_node"]:
			noiseimage        = get_im(os.path.join(old_previousoutputdir, "bckgnoise.hdf"))
			noiseimage1       = get_im(os.path.join(old_refinement_iter_dir, "bckgnoise.hdf"))
			params            = read_text_row(os.path.join(old_refinement_iter_dir, "params_%03d.txt"%selected_iter))
			params_last_iter  = read_text_row(os.path.join(old_previousoutputdir, "params_%03d.txt"%(selected_iter-1)))
			refang            = read_text_row(os.path.join(old_refinement_iter_dir, "refang.txt"))
			rshifts           = read_text_row(os.path.join(old_refinement_iter_dir, "rshifts.txt"))
			chunk_one         = read_text_file(os.path.join(old_refinement_iter_dir, "chunk_0_%03d.txt"%selected_iter))
			chunk_two         = read_text_file(os.path.join(old_refinement_iter_dir, "chunk_1_%03d.txt"%selected_iter))
			error_threshold   = read_text_row(os.path.join(old_refinement_iter_dir, "error_thresholds_%03d.txt"%selected_iter))
		else:
			params           = 0
			refang           = 0
			rshifts          = 0
			chunk_one        = 0
			chunk_two        = 0
			params_last_iter = 0
		params            = wrap_mpi_bcast(params, Blockdata["main_node"], comm)
		params_last_iter  = wrap_mpi_bcast(params_last_iter, Blockdata["main_node"], comm)
		refang            = wrap_mpi_bcast(refang,    Blockdata["main_node"], comm)
		rshifts           = wrap_mpi_bcast(rshifts,   Blockdata["main_node"], comm)
		chunk_one         = wrap_mpi_bcast(chunk_one, Blockdata["main_node"], comm)
		chunk_two         = wrap_mpi_bcast(chunk_two, Blockdata["main_node"], comm)
		chunk_dict = {}
		for a in chunk_one: chunk_dict[a] = 0
		for b in chunk_two: chunk_dict[b] = 1
		### handle the selected cluster
		# create directories
		main0_dir                 = os.path.join(masterdir, "main000")
		iter_dir                  = os.path.join(masterdir, "main%03d"%selected_iter)
		previousoutputdir         = os.path.join(masterdir, "main%03d"%(selected_iter-1))
		new_oldparamstructure_dir = os.path.join(iter_dir,"oldparamstructure")
		if Blockdata["subgroup_myid"] == Blockdata["main_node"]:
			if not os.path.exists(iter_dir):
				cmd         = "{} {}".format("mkdir", iter_dir)
				cmdexecute(cmd)
			if not os.path.exists(main0_dir):
				cmd         = "{} {}".format("mkdir", main0_dir)
				cmdexecute(cmd)
			if not os.path.exists(new_oldparamstructure_dir):	
				cmd = "{} {}".format("mkdir", new_oldparamstructure_dir)
				cmdexecute(cmd)
			if not os.path.exists(previousoutputdir):
				cmd = "{} {}".format("mkdir", previousoutputdir)
				cmdexecute(cmd)
		mpi_barrier(comm)
		# load selected iter
		new_chunk_one                  = []
		new_chunk_two                  = []
		new_params                     = []
		new_params_chunk_one           = []
		new_params_chunk_two           = []
		new_params_chunk_one_last_iter = []
		new_params_chunk_two_last_iter = []	
		Tracker["avgvaradj"] = [0.0, 0.0]		
		for index_of_particle in xrange(len(cluster)): 
			if chunk_dict[cluster[index_of_particle]] == 0:  
				new_chunk_one.append(cluster[index_of_particle])
				new_params_chunk_one.append(params[cluster[index_of_particle]])
				new_params_chunk_one_last_iter.append(params_last_iter[cluster[index_of_particle]])
				Tracker["avgvaradj"][0] += params[cluster[index_of_particle]][7]
			else:				                             
				new_chunk_two.append(cluster[index_of_particle])
				new_params_chunk_two.append(params[cluster[index_of_particle]])
				new_params_chunk_two_last_iter.append(params_last_iter[cluster[index_of_particle]]) 
				Tracker["avgvaradj"][1] += params[cluster[index_of_particle]][7] 
			new_params.append(params[cluster[index_of_particle]])		
		selected_new_params = new_params
		if Blockdata["subgroup_myid"] == Blockdata["main_node"]:# some numbers and path are required to be modified
			Tracker["constants"]["masterdir"] = masterdir
			Tracker["directory"]              = iter_dir
			try:    sym = Tracker["constants"]["sym"] # For those generated by old version meridians
			except: sym = Tracker["constants"]["symmetry"]
			Tracker["constants"]["symmetry"]  = sym
			Tracker["best"]                   = selected_iter +2 # reset the best to arbitrary iteration
			Tracker["bestres"]                = 0
			Tracker["no_improvement"]         = 0
			Tracker["no_params_changes"]      = 0
			Tracker["pixercutoff"]            = 0
			Tracker["saturated_sampling"]     = False
			Tracker["is_converged"]           = False
			Tracker["large_at_Nyquist"]       = False
			Tracker["previousoutputdir"]      = previousoutputdir
			Tracker["refvol"]                 = os.path.join(iter_dir, "vol_0_%03d.hdf"%selected_iter)
			Tracker["mainiteration"]          = selected_iter
			if shell_line_command: update_tracker(shell_line_command) # the updated could be any refinement parameters that user wish to make change
			error_angles, error_shifts = params_changes((new_params_chunk_one + new_params_chunk_two), (new_params_chunk_one_last_iter + new_params_chunk_two_last_iter))
			# varibles in Tracker to be updated
			if Tracker["constants"]["mask3D"]: 
				Tracker["constants"]["mask3D"] = os.path.join(option_old_refinement_dir, "../", Tracker["constants"]["mask3D"])
				if not os.path.exists(Tracker["constants"]["mask3D"]): Tracker["constants"]["mask3D"] =  None		
			noiseimage.write_image(os.path.join(Tracker["previousoutputdir"], "bckgnoise.hdf"))
			noiseimage1.write_image(os.path.join(iter_dir, "bckgnoise.hdf"))
			write_text_file(cluster, os.path.join(iter_dir, "indexes_%03d.txt"%selected_iter))
			write_text_row(refang, os.path.join(iter_dir, "refang.txt"))
			write_text_row(rshifts, os.path.join(iter_dir, "rshifts.txt"))
			write_text_row(new_params_chunk_one, os.path.join(iter_dir, "params-chunk_0_%03d.txt"%selected_iter))
			write_text_row(new_params_chunk_two, os.path.join(iter_dir, "params-chunk_1_%03d.txt"%selected_iter))
			write_text_row(new_params_chunk_one_last_iter, os.path.join(Tracker["previousoutputdir"], "params-chunk_0_%03d.txt"%(selected_iter -1)))
			write_text_row(new_params_chunk_two_last_iter, os.path.join(Tracker["previousoutputdir"], "params-chunk_1_%03d.txt"%(selected_iter -1)))
			write_text_file(new_chunk_one, os.path.join(iter_dir, "chunk_0_%03d.txt"%selected_iter))
			write_text_file(new_chunk_two, os.path.join(iter_dir, "chunk_1_%03d.txt"%selected_iter))
			write_text_row(new_params, os.path.join(iter_dir, "params_%03d.txt"%selected_iter))
			write_text_row([[error_angles, error_shifts]], os.path.join(iter_dir, "error_thresholds_%03d.txt"%selected_iter))		
			Tracker["nima_per_chunk"] = [len(new_chunk_one), len(new_chunk_two)]
			Tracker["avgvaradj"][0] /=float(len(new_chunk_one))
			Tracker["avgvaradj"][1] /=float(len(new_chunk_two))
			fout = open(os.path.join(iter_dir, "Tracker_%03d.json"%selected_iter),"w")
			json.dump(Tracker, fout)
			fout.close()
			# now partition new indexes into new oldparamstructure
			nproc_dict = {}
			for ichunk in xrange(2):
				if ichunk == 0: total_stack_on_chunk = len(chunk_one)
				else: 	        total_stack_on_chunk = len(chunk_two)
				for myproc in xrange(nproc_old_ref3d):
					image_start,image_end = MPI_start_end(total_stack_on_chunk, nproc_old_ref3d, myproc)
					for index_of_particle in xrange(image_start, image_end):
						if ichunk == 0: nproc_dict[chunk_one[index_of_particle]] = [ichunk, myproc, index_of_particle - image_start]
						else: nproc_dict[chunk_two[index_of_particle]] = [ichunk, myproc, index_of_particle - image_start]
		else: nproc_dict = 0
		nproc_dict = wrap_mpi_bcast(nproc_dict, Blockdata["main_node"], comm)
		### parse previous nproc in refinement to current nproc
		#proc_start, proc_end = MPI_start_end(Blockdata["nsubset"], Blockdata["nsubset"], Blockdata["subgroup_myid"])
		for ichunk in xrange(2):
			oldparams = []
			if ichunk == 0: total_stack_on_chunk = len(new_chunk_one)
			else: 	        total_stack_on_chunk = len(new_chunk_two)
			image_start,image_end = MPI_start_end(total_stack_on_chunk, Blockdata["nsubset"], Blockdata["subgroup_myid"])
			for index_of_particle in xrange(image_start,image_end):
				if ichunk == 0:   [old_chunk, old_proc, old_index_of_particle] = nproc_dict[new_chunk_one[index_of_particle]]
				else: 	          [old_chunk, old_proc, old_index_of_particle] = nproc_dict[new_chunk_two[index_of_particle]]
				fout = open(os.path.join(old_oldparamstructure_dir, "oldparamstructure_%d_%03d_%03d.json"%(old_chunk, old_proc, selected_iter)),"r")
				old_oldparams 	= convert_json_fromunicode(json.load(fout))
				fout.close()
				oldparams.append(old_oldparams[old_index_of_particle])
			if Blockdata["subgroup_myid"] == Blockdata["main_node"]:
				fout = open(os.path.join(new_oldparamstructure_dir, "oldparamstructure_%d_%03d_%03d.json"%(ichunk, Blockdata["subgroup_myid"], selected_iter)), "w")
				json.dump(oldparams, fout)
				fout.close()
				del oldparams
			mpi_barrier(comm)
			for iproc in xrange(1, Blockdata["nsubset"]): # always skip main node
				if iproc == Blockdata["subgroup_myid"]:wrap_mpi_send(oldparams, Blockdata["main_node"], comm)
				if Blockdata["subgroup_myid"] ==  Blockdata["main_node"]:
					dummy = wrap_mpi_recv(iproc, comm)
					fout = open(os.path.join(new_oldparamstructure_dir, "oldparamstructure_%d_%03d_%03d.json"%(ichunk, iproc, selected_iter)), "w")
					json.dump(dummy, fout)
					fout.close()
					del dummy
				mpi_barrier(comm)
			mpi_barrier(comm)
		mpi_barrier(comm)				
		### <<<-------load 0 iteration
		selected_iter              = 0
		old_refinement_iter_dir    = os.path.join(option_old_refinement_dir, "main%03d"%selected_iter)
		old_oldparamstructure_dir  = os.path.join(old_refinement_iter_dir, "oldparamstructure")
		iter_dir                   = os.path.join(masterdir, "main%03d"%selected_iter)
	
		if Blockdata["subgroup_myid"] == Blockdata["main_node"]:
			fout     = open(os.path.join(old_refinement_iter_dir, "Tracker_%03d.json"%selected_iter),"r")
			Tracker  = convert_json_fromunicode(json.load(fout))
			fout.close()
		else: Tracker = 0
		Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], comm) # balance processors
		Tracker["constants"]["orgstack"] = orgstack
		if Blockdata["subgroup_myid"] == Blockdata["main_node"]:
			if not os.path.exists(iter_dir):
				cmd         = "{} {}".format("mkdir", iter_dir)
				cmdexecute(cmd)
		mpi_barrier(comm)
	
		if Blockdata["subgroup_myid"] == Blockdata["main_node"]:
			params               = read_text_row(os.path.join(old_refinement_iter_dir,  "params_%03d.txt"%selected_iter))
			chunk_one            = read_text_file(os.path.join(old_refinement_iter_dir, "chunk_0_%03d.txt"%selected_iter))
			chunk_two            = read_text_file(os.path.join(old_refinement_iter_dir, "chunk_1_%03d.txt"%selected_iter))
			particle_group_one   = read_text_file(os.path.join(old_refinement_iter_dir, "particle_groups_0.txt"))
			particle_group_two   = read_text_file(os.path.join(old_refinement_iter_dir, "particle_groups_1.txt"))
			groupids             = read_text_file(os.path.join(old_refinement_iter_dir, "groupids.txt"))
	
		else:
			groupids   = 0
			params     = 0
			refang     = 0
			rshifts    = 0
			chunk_one  = 0
			chunk_two  = 0
			particle_group_one = 0
			particle_group_two = 0
		params              = wrap_mpi_bcast(params,    Blockdata["main_node"], comm)
		chunk_one           = wrap_mpi_bcast(chunk_one, Blockdata["main_node"], comm)
		chunk_two           = wrap_mpi_bcast(chunk_two, Blockdata["main_node"], comm)
		particle_group_one  = wrap_mpi_bcast(particle_group_one, Blockdata["main_node"], comm)
		particle_group_two  = wrap_mpi_bcast(particle_group_two, Blockdata["main_node"], comm)
		groupids            = wrap_mpi_bcast(groupids, Blockdata["main_node"], comm)
		group_ids_dict = {}
		for iptl in xrange(len(particle_group_one)): group_ids_dict[chunk_one[iptl]] =  particle_group_one[iptl]
		for iptl in xrange(len(particle_group_two)): group_ids_dict[chunk_two[iptl]] =  particle_group_two[iptl]
		chunk_dict    = {}
		for a in chunk_one: chunk_dict[a] = 0
		for b in chunk_two: chunk_dict[b] = 1
		### handle the selected cluster
		new_chunk_one          = []
		new_chunk_two          = []
		new_params             = []
		new_params_chunk_one   = []
		new_params_chunk_two   = []
		new_particle_group_one = []
		new_particle_group_two = []
		for index_of_particle in xrange(len(cluster)): 
			if chunk_dict[cluster[index_of_particle]] == 0:  
				new_chunk_one.append(cluster[index_of_particle])
				new_params_chunk_one.append(params[cluster[index_of_particle]])
			else:				                             
				new_chunk_two.append(cluster[index_of_particle])
				new_params_chunk_two.append(params[cluster[index_of_particle]])
			new_params.append(params[cluster[index_of_particle]])	
		for iptl in xrange(len(new_chunk_one)):new_particle_group_one.append(group_ids_dict[new_chunk_one[iptl]])
		for iptl in xrange(len(new_chunk_two)):new_particle_group_two.append(group_ids_dict[new_chunk_two[iptl]])
		if Blockdata["subgroup_myid"] == Blockdata["main_node"]:# some numbers and path are required to be modified
			# varibles in Tracker to be updated
			Tracker["constants"]["masterdir"] = masterdir
			Tracker["previousoutputdir"]      = Tracker["directory"]
			Tracker["refvol"]                 = os.path.join(iter_dir, "vol_0_%03d.hdf"%selected_iter)
			Tracker["mainiteration"]          = selected_iter
		
			if Tracker["constants"]["mask3D"]: 
				Tracker["constants"]["mask3D"]= os.path.join(option_old_refinement_dir, "../", Tracker["constants"]["mask3D"])
				if not os.path.exists(Tracker["constants"]["mask3D"]): Tracker["constants"]["mask3D"] =  None
			write_text_file(cluster, os.path.join(iter_dir, "indexes_%03d.txt"%selected_iter))
			write_text_file(groupids, os.path.join(iter_dir, "groupids.txt"))
			write_text_row(new_params, os.path.join(iter_dir, "params_%03d.txt"%selected_iter))
			write_text_row(new_params_chunk_one, os.path.join(iter_dir, "params-chunk_0_%03d.txt"%selected_iter))
			write_text_row(new_params_chunk_two, os.path.join(iter_dir, "params-chunk_1_%03d.txt"%selected_iter))
			write_text_file(new_chunk_one, os.path.join(iter_dir, "chunk_0_%03d.txt"%selected_iter))
			write_text_file(new_chunk_two, os.path.join(iter_dir, "chunk_1_%03d.txt"%selected_iter))
			write_text_file(new_particle_group_one, os.path.join(iter_dir, "particle_groups_0.txt"))
			write_text_file(new_particle_group_two, os.path.join(iter_dir, "particle_groups_1.txt"))
		Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], comm) # balance processors
		Tracker["constants"]["orgstack"] = orgstack
		if Blockdata["subgroup_myid"] == Blockdata["main_node"]:# some numbers and path are required to be modified
			fout = open(os.path.join(iter_dir, "Tracker_%03d.json"%selected_iter),"w")
			json.dump(Tracker, fout)
			fout.close()
		mpi_barrier(comm)
	return
	
#######functions for ctreffromsorting
def ctrefromsorting_rec3d_faked_iter(masterdir, selected_iter=-1, rec3d_image_size =-1, comm = -1):
	global Tracker, Blockdata
	import json
	#from mpi import mpi_barrier, MPI_COMM_WORLD
	if comm ==-1: comm =  MPI_COMM_WORLD
	if Blockdata["subgroup_myid"]>-1:
		Tracker["directory"]          = os.path.join(masterdir, "main%03d"%selected_iter)
		Tracker["previousoutputdir"]  = os.path.join(masterdir, "main%03d"%(selected_iter-1))
		oldparamstructure =[[],[]]
		newparamstructure =[[],[]]
		projdata          = [[model_blank(1,1)], [model_blank(1,1)]]
		original_data     = [None,None]
		oldparams         = [[],[]]
		partids           = [None, None]
		partstack         = [None, None]
		final_dir         = Tracker["directory"] 
		if Blockdata["subgroup_myid"] == Blockdata["main_node"]:
			fout = open(os.path.join(Tracker["directory"], "Tracker_%03d.json"%selected_iter),"r")
			Tracker = convert_json_fromunicode(json.load(fout))
			fout.close()
			if rec3d_image_size !=-1: Tracker["nxinit"] = rec3d_image_size
		else: Tracker = 0
		Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], comm) # balance processors
		Blockdata["accumulatepw"] = [[],[]]
		if selected_iter ==-1: ERROR("Iteration number has to be determined in advance.","ctrefromsorting_rec3d_faked_iter",1, Blockdata["subgroup_myid"])  
		carryon = 1
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		if(Blockdata["subgroup_myid"] == Blockdata["main_node"]):
			print(line, "ctrefromsorting_rec3d_faked_iter")
			print(line, "Reconstruction uses solution from  %d iteration"%selected_iter)
			print(line, "Reconstruction image size is:  %d"%(Tracker["nxinit"]))
			print(line, "Reconstruction directory is %s"%(Tracker["directory"]))
		if(Blockdata["subgroup_myid"] == Blockdata["main_node"]):
			try:
				refang  = read_text_row( os.path.join(Tracker["directory"], "refang.txt"))
				rshifts = read_text_row( os.path.join(Tracker["directory"], "rshifts.txt"))
			except:carryon =0
		else:
			refang  = 0
			rshifts = 0
		carryon = bcast_number_to_all(carryon, source_node = Blockdata["main_node"], mpi_comm = comm)
		if carryon == 0: ERROR("Failed to read refang and rshifts: %s %s "%(os.path.join(Tracker["directory"], "refang.txt"), os.path.join(Tracker["directory"], \
			"rshifts.txt")), "ctrefromsorting_rec3d_faked_iter", 1,Blockdata["subgroup_myid"])
		refang  = wrap_mpi_bcast(refang, Blockdata["main_node"], comm)
		rshifts = wrap_mpi_bcast(rshifts, Blockdata["main_node"], comm)
		partids =[None, None]
		if(Blockdata["subgroup_myid"] == Blockdata["main_node"]):
			cmd = "{} {} ".format("mkdir ",os.path.join(Tracker["directory"], "tempdir"))
			if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")): cmdexecute(cmd)
			l = 0
			for procid in xrange(2):
				partids[procid] = os.path.join(Tracker["directory"],"chunk_%01d_%03d.txt"%(procid,Tracker["mainiteration"]))
				l += len(read_text_file(partids[procid]))
		else:l  = 0
		l = bcast_number_to_all(l, source_node = Blockdata["main_node"], mpi_comm = comm)
		norm_per_particle = [[],[]]
		for procid in xrange(2):
			if procid ==0: original_data[1] = None	
			partids[procid]   = os.path.join(Tracker["directory"],"chunk_%01d_%03d.txt"%(procid,Tracker["mainiteration"]))
			partstack[procid] = os.path.join(Tracker["constants"]["masterdir"],"main%03d"%(Tracker["mainiteration"]-1),"params-chunk_%01d_%03d.txt"%(procid,(Tracker["mainiteration"]-1)))
			if Blockdata["subgroup_myid"] == Blockdata["main_node"]:
				nproc_previous = 0
				while os.path.exists(os.path.join(Tracker["directory"],"oldparamstructure","oldparamstructure_%01d_%03d_%03d.json"%(procid, nproc_previous, Tracker["mainiteration"]))):nproc_previous += 1
				psize = len(read_text_file(partids[procid]))
			for jproc in xrange(Blockdata["nsubset"]):
				if Blockdata["subgroup_myid"] == Blockdata["main_node"]:
					dummy = []
					im_start, im_end = MPI_start_end(psize, Blockdata["nsubset"], jproc)
					istart_old_proc_id = -1
					iend_old_proc_id   = -1
					plist = []
					for iproc_old in xrange(nproc_previous):
						im_start_old, im_end_old = MPI_start_end(psize, nproc_previous, iproc_old)
						if (im_start>= im_start_old) and im_start <=im_end_old: istart_old_proc_id = iproc_old
						if (im_end>= im_start_old) and im_end <=im_end_old: iend_old_proc_id = iproc_old
						plist.append([im_start_old, im_end_old])
					ptl_on_this_cpu = im_start
					for iproc_index_old in xrange(istart_old_proc_id, iend_old_proc_id+1):
						fout = open(os.path.join(final_dir,"oldparamstructure","oldparamstructure_%01d_%03d_%03d.json"%(procid,iproc_index_old,Tracker["mainiteration"])),'r')
						oldparamstructure_on_old_cpu = convert_json_fromunicode(json.load(fout))
						fout.close()
						mlocal_id_on_old = ptl_on_this_cpu - plist[iproc_index_old][0]
						while (mlocal_id_on_old<len(oldparamstructure_on_old_cpu)) and (ptl_on_this_cpu<im_end):
							dummy.append(oldparamstructure_on_old_cpu[mlocal_id_on_old])
							ptl_on_this_cpu  +=1
							mlocal_id_on_old +=1
					del oldparamstructure_on_old_cpu
					if jproc== Blockdata["main_node"]: oldparamstructure[procid] = dummy
					else: wrap_mpi_send(dummy, jproc, comm)
				else:
					if Blockdata["subgroup_myid"] ==jproc: oldparamstructure[procid] = wrap_mpi_recv(0, comm)
				mpi_barrier(comm)
			mpi_barrier(comm)
			#####
			original_data[procid], oldparams[procid] = getindexdata(partids[procid], partstack[procid], \
					os.path.join(Tracker["constants"]["masterdir"],"main000", "particle_groups_%01d.txt"%procid), \
					original_data[procid], small_memory = Tracker["constants"]["small_memory"], \
					nproc = Blockdata["nsubset"], myid = Blockdata["subgroup_myid"], mpi_comm = comm)										
			temp = Tracker["directory"]
			Tracker["directory"] = os.path.join(Tracker["directory"], "tempdir")
			mpi_barrier(comm)
			if procid == 0: compute_sigma([[]]*l, [[]]*l, len(oldparams[0]), True, myid = Blockdata["subgroup_myid"], mpi_comm = comm)
			Tracker["directory"] = temp
			mpi_barrier(comm)
			projdata[procid] = get_shrink_data_final(Tracker["nxinit"], procid, original_data[procid], oldparams[procid],\
			   return_real = False, preshift = True, apply_mask = False, nonorm = True)
			for ipar in xrange(len(oldparams[procid])):norm_per_particle[procid].append(oldparams[procid][ipar][7])
			oldparams[procid] = []
			original_data[procid] = None
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			if(Blockdata["subgroup_myid"] == Blockdata["nodes"][procid]): print(line, "3-D reconstruction of group %d"%procid)
			Tracker["maxfrad"] = Tracker["nxinit"] //2
			do3d(procid, projdata[procid], oldparamstructure[procid], refang, rshifts, norm_per_particle[procid], myid = Blockdata["subgroup_myid"], mpi_comm = comm)
			projdata[procid]          = []
			oldparamstructure[procid] = []
			norm_per_particle[procid] = []
			mpi_barrier(comm)
	return selected_iter
	
def out_fsc(f):
	global Tracker, Blockdata
	print(" ")
	print("  driver FSC  after  iteration#%3d"%Tracker["mainiteration"])
	print("  %4d        %7.2f         %5.3f"%(0,1000.00,f[0]))
	for i in xrange(1,len(f)): print("  %4d        %7.2f         %5.3f"%(i,Tracker["constants"]["pixel_size"]*Tracker["constants"]["nnxo"]/float(i),f[i]))
	print(" ")
	
### functions for faked rec3d from subsets
def params_changes( params, oldparams ):
	#  Indexes contain list of images processed - sorted integers, subset of the full range.
	#  params - contain parameters associated with these images
	#  Both lists can be of different sizes, so we have to find a common subset
	#  We do not compensate for random changes of grids.
	from utilities    	import getang3
	from utilities    	import rotate_shift_params
	from pixel_error  	import max_3D_pixel_error
	from EMAN2        	import Vec2f
	from math 			import sqrt
	import sets

	n = len(params)
	anger       = 0.0
	shifter     = 0.0
	#  The shifter is given in the full scale displacement
	for i in xrange(n):
		shifter     += (params[i][3] - oldparams[i][3] )**2 + (params[i][4] - oldparams[i][4] )**2
		anger += get_anger(params[i][0:3], oldparams[i][0:3],Tracker["constants"]["symmetry"])
	return round(anger/n,5), round(sqrt(shifter/n),5)
	
def get_anger(angle1, angle2, sym="c1"):
	from math import acos, pi
	R1               = Transform({"type":"spider","phi":  angle1[0], "theta":  angle1[1],  "psi": angle1[2]})
	R2               = Transform({"type":"spider","phi":  angle2[0], "theta":  angle2[1],  "psi": angle2[2]})
	R2               = R2.get_sym_proj(sym)
	axes_dis_min     = 1.0e23
	for isym in xrange(len(R2)):
		A1 		         = R1.get_matrix()
		A2 		         = R2[isym].get_matrix()
		X1               = A1[0]*A2[0] + A1[1]*A2[1] + A1[2]*A2[2] 
		X2               = A1[4]*A2[4] + A1[5]*A2[5] + A1[6]*A2[6]
		X3               = A1[8]*A2[8] + A1[9]*A2[9] + A1[10]*A2[10] 
		axes_dis         = acos(max(min(X1,1.),-1.0))*180./pi +acos(max(min(X2,1.),-1.0))*180./pi +acos(max(min(X3,1.),-1.0))*180./pi/3.0
		axes_dis_min     = min(axes_dis_min, axes_dis)
	return axes_dis_min

def compute_sigma(projdata, params, first_procid, dryrun = False, myid = -1, mpi_comm = -1):
	global Tracker, Blockdata
	# Input stack of particles with all params in header
	# Output: 1/sigma^2 and a dictionary
	#  It could be a parameter
	if( mpi_comm < 0 ): mpi_comm = MPI_COMM_WORLD
	npad = 1
	if  dryrun:
		#tsd = model_blank(nv + nv//2,len(sd), 1, 1.0)
		#tocp = model_blank(len(sd), 1, 1, 1.0)
		if( myid == Blockdata["main_node"] ):
			tsd = get_im(os.path.join(Tracker["previousoutputdir"],"bckgnoise.hdf"))
			tsd.write_image(os.path.join(Tracker["directory"],"bckgnoise.hdf"))
			nnx = tsd.get_xsize()
			nny = tsd.get_ysize()
		else:
			nnx = 0
			nny = 0
		nnx = bcast_number_to_all(nnx, source_node = Blockdata["main_node"], mpi_comm = mpi_comm)
		nny = bcast_number_to_all(nny, source_node = Blockdata["main_node"], mpi_comm = mpi_comm)
		if( myid != Blockdata["main_node"] ): tsd = model_blank(nnx,nny, 1, 1.0)
		bcast_EMData_to_all(tsd, myid, source_node = Blockdata["main_node"], comm = mpi_comm)
	else:
		if( myid == Blockdata["main_node"] ): ngroups = len(read_text_file(os.path.join(Tracker["constants"]["masterdir"],"main000", "groupids.txt")))
		else: ngroups = 0
		ngroups = bcast_number_to_all(ngroups, source_node = Blockdata["main_node"], mpi_comm = mpi_comm)
		ndata = len(projdata)
		nx = Tracker["constants"]["nnxo"]
		mx = npad*nx
		nv = mx//2+1
		"""
		#  Inverted Gaussian mask
		invg = model_gauss(Tracker["constants"]["radius"],nx,nx)
		invg /= invg[nx//2,nx//2]
		invg = model_blank(nx,nx,1,1.0) - invg
		"""

		mask = model_circle(Tracker["constants"]["radius"],nx,nx)
		tsd = model_blank(nv + nv//2, ngroups)

		#projdata, params = getalldata(partstack, params, myid, Blockdata["nproc"])
		'''
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
			params = [params[i][j]  for i in xrange(len(params))   for j in xrange(5)]
		else:           params = [0.0]*(5*ndata)
		params = bcast_list_to_all(params, myid, source_node=Blockdata["main_node"])
		params = [[params[i*5+j] for j in xrange(5)] for i in xrange(ndata)]
		'''
		if(Blockdata["accumulatepw"] == None):
			Blockdata["accumulatepw"] = [[],[]]
			doac = True
		else:  doac = False
		tocp = model_blank(ngroups)
		tavg = model_blank(nx,nx)
		for i in xrange(ndata):  # apply_shift; info_mask; norm consistent with get_shrink_data
			indx = projdata[i].get_attr("particle_group")
			phi,theta,psi,sx,sy = params[i][0],params[i][1],params[i][2],params[i][3],params[i][4]
			stmp = cyclic_shift( projdata[i], int(round(sx)), int(round(sy)))
			st = Util.infomask(stmp, mask, False)
			stmp -=st[0]
			stmp /=st[1]
			temp = cosinemask(stmp, radius = Tracker["constants"]["radius"], s = 0.0)
			Util.add_img(tavg, temp)
			sig = Util.rotavg_fourier( temp )
			#sig = rops(pad(((cyclic_shift( projdata[i], int(sx), int(round(sy)) ) - st[0])/st[1]), mx,mx,1,0.0))
			#sig = rops(pad(((cyclic_shift(projdata, int(round(params[i][-2])), int(round(params[i][-1])) ) - st[0])/st[1])*invg, mx,mx,1,0.0))
			for k in xrange(nv):tsd.set_value_at(k,indx,tsd.get_value_at(k,indx)+sig[k])
			tocp[indx] += 1
		####for lll in xrange(len(Blockdata["accumulatepw"])):  print(myid,ndata,lll,len(Blockdata["accumulatepw"][lll]))
		reduce_EMData_to_root(tsd,  myid, Blockdata["main_node"],  mpi_comm)
		reduce_EMData_to_root(tocp, myid, Blockdata["main_node"], mpi_comm)
		reduce_EMData_to_root(tavg, myid, Blockdata["main_node"], mpi_comm)
		if( myid == Blockdata["main_node"]):
			Util.mul_scalar(tavg, 1.0/float(sum(Tracker["nima_per_chunk"])))
			sig = Util.rotavg_fourier( tavg )
			#for k in xrange(1,nv):  print("  BACKG  ",k,tsd.get_value_at(k,0)/tocp[0] ,sig[k],tsd.get_value_at(k,0)/tocp[0] - sig[k])
			tmp1 = [0.0]*nv
			tmp2 = [0.0]*nv
			for i in xrange(ngroups):
				for k in xrange(1,nv):
					qt = tsd.get_value_at(k,i)/tocp[i] - sig[k]
					if( qt > 0.0 ):	tmp1[k] = 2.0/qt
				#smooth
				tmp1[0] = tmp1[1]
				tmp1[-1] = tmp1[-2]
				for ism in xrange(0):  #2
					for k in xrange(1,nv-1):  tmp2[k] = (tmp1[k-1]+tmp1[k]+tmp1[k+1])/3.0
					for k in xrange(1,nv-1):  tmp1[k] = tmp2[k]
				#  We will keep 0-element the same as first tsd.set_value_at(0,i,1.0)
				for k in xrange(1,nv):tsd.set_value_at(k,i,tmp1[k])
				tsd.set_value_at(0,i,1.0)
			tsd.write_image(os.path.join(Tracker["directory"],"bckgnoise.hdf"))
		bcast_EMData_to_all(tsd, myid, source_node = 0, comm = mpi_comm)
	nnx = tsd.get_xsize()
	nny = tsd.get_ysize()
	Blockdata["bckgnoise"] = []
	for i in xrange(nny):
		prj = model_blank(nnx)
		for k in xrange(nnx): prj[k] = tsd.get_value_at(k,i)
		Blockdata["bckgnoise"].append(prj)  #  1.0/sigma^2
	return
###
def do3d(procid, data, newparams, refang, rshifts, norm_per_particle, myid, mpi_comm = -1):
	global Tracker, Blockdata
	#  Without filtration
	from reconstruction import recons3d_trl_struct_MPI
	if( mpi_comm < -1 ): mpi_comm = MPI_COMM_WORDLD
	if Blockdata["subgroup_myid"]== Blockdata["main_node"]:
		if( procid == 0 ):
			cmd = "{} {}".format("mkdir", os.path.join(Tracker["directory"], "tempdir") )
			if os.path.exists(os.path.join(Tracker["directory"], "tempdir")): print("tempdir exists")
			else: cmdexecute(cmd)
		
	"""
	tvol, tweight, trol = recons3d_4nnstruct_MPI(myid = Blockdata["subgroup_myid"], main_node = Blockdata["nodes"][procid], prjlist = data, \
											paramstructure = newparams, refang = refang, delta = Tracker["delta"], CTF = Tracker["constants"]["CTF"],\
											upweighted = False, mpi_comm = mpi_comm, \
											target_size = (2*Tracker["nxinit"]+3), avgnorm = Tracker["avgvaradj"][procid], norm_per_particle = norm_per_particle)
	"""
	shrinkage = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	tvol, tweight, trol = recons3d_trl_struct_MPI(myid = Blockdata["subgroup_myid"], main_node = Blockdata["main_node"], prjlist = data, \
											paramstructure = newparams, refang = refang, rshifts_shrank = [[q[0]*shrinkage,q[1]*shrinkage] for q in rshifts], \
											delta = Tracker["delta"], CTF = Tracker["constants"]["CTF"], upweighted = False, mpi_comm = mpi_comm, \
											target_size = (2*Tracker["nxinit"]+3), avgnorm = Tracker["avgvaradj"][procid], norm_per_particle = norm_per_particle)
	if Blockdata["subgroup_myid"] == Blockdata["main_node"]:
		tvol.set_attr("is_complex",0)
		tvol.write_image(os.path.join(Tracker["directory"], "tempdir", "tvol_%01d_%03d.hdf"%(procid,Tracker["mainiteration"])))
		tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%01d_%03d.hdf"%(procid,Tracker["mainiteration"])))
		trol.write_image(os.path.join(Tracker["directory"], "tempdir", "trol_%01d_%03d.hdf"%(procid,Tracker["mainiteration"])))
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print(line,"Executed successfully backprojection for group ",procid)
	mpi_barrier(mpi_comm)
	return
##
def getindexdata(partids, partstack, particle_groups, original_data = None, small_memory= True, nproc =-1, myid = -1, mpi_comm = -1):
	global Tracker, Blockdata
	# The function will read from stack a subset of images specified in partids
	#   and assign to them parameters from partstack
	# So, the lengths of partids and partstack are the same.
	#  The read data is properly distributed among MPI threads.
	if( mpi_comm < 0 ):  mpi_comm = MPI_COMM_WORLD
	#print("partids", partids)
	#print("partstack", partstack)
	#print("particle_groups", particle_groups)
	#exit()
	from applications import MPI_start_end
	#  parameters
	if( myid == 0 ):  partstack = read_text_row(partstack)
	else:  			  partstack = 0
	partstack = wrap_mpi_bcast(partstack, 0, mpi_comm)
	#  particles IDs
	if( myid == 0 ):  partids = read_text_file(partids)
	else:          	  partids = 0
	partids = wrap_mpi_bcast(partids, 0, mpi_comm)
	#  Group assignments
	if( myid == 0 ):	group_reference = read_text_file(particle_groups)
	else:          		group_reference = 0
	group_reference = wrap_mpi_bcast(group_reference, 0, mpi_comm)

	im_start, im_end = MPI_start_end(len(partstack), nproc, myid)
	partstack = partstack[im_start:im_end]
	partids   = partids[im_start:im_end]
	group_reference = group_reference[im_start:im_end]
	'''
	particles_on_node = []
	parms_on_node     = []
	for i in xrange( group_start, group_end ):
		particles_on_node += lpartids[group_reference[i][2]:group_reference[i][3]+1]  #  +1 is on account of python idiosyncrasies
		parms_on_node     += partstack[group_reference[i][2]:group_reference[i][3]+1]


	Blockdata["nima_per_cpu"][procid] = len(particles_on_node)
	#ZZprint("groups_on_thread  ",Blockdata["myid"],procid, Tracker["groups_on_thread"][procid])
	#ZZprint("  particles  ",Blockdata["myid"],Blockdata["nima_per_cpu"][procid],len(parms_on_node))
	'''
	"""
            17            28            57            84    5
            18            14            85            98    6
            19            15            99           113    7
            25            20           114           133    8
            29             9           134           142    9

	"""
	#print("getindexdata", Tracker["constants"]["orgstack"])
	#print(len(partids), Blockdata["myid"])
	if( original_data == None or small_memory):
		original_data = EMData.read_images(Tracker["constants"]["orgstack"], partids)
		for im in xrange( len(original_data) ): 
			original_data[im].set_attr("particle_group", group_reference[im])
	return original_data, partstack
#######
def do3d_sorting_groups_rec3d(iteration, masterdir, log_main):
	global Tracker, Blockdata
	from utilities import get_im
	# reconstruct final unfiltered volumes from sorted clusters
	keepgoing = 1
	#if(Blockdata["myid"] == Blockdata["nodes"][0]):
	#	cmd = "{} {}".format("mkdir", os.path.join(Tracker["directory"], "tempdir"))
	#	if os.path.exists(os.path.join(Tracker["directory"], "tempdir")): print("tempdir exists")
	#	else:                                                             cmdexecute(cmd)
	### ====	
	fsc143                          =   0
	fsc05                           =   0
	Tracker["fsc143"]				=	0
	Tracker["fsc05"]				=	0
	res_05 						    =	Tracker["number_of_groups"]*[0]
	res_143 					    =	Tracker["number_of_groups"]*[0]
	Tracker["directory"]            =   masterdir
	Tracker["constants"]["masterdir"] = masterdir
	### ====
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] == Blockdata["main_node"]:
		msg = "trl_iter starts to do volumes, number of groups is  %d"%Tracker["number_of_groups"]
		print(line, msg)
		log_main.add(msg)
	for index_of_colors in xrange(Blockdata["no_of_groups"]):	
		group_start, group_end = MPI_volume_start_end(Tracker["number_of_groups"], Blockdata["no_of_groups"], index_of_colors)
		if Blockdata["color"] == index_of_colors:  # It has to be 1 to avoid problem with tvol1 not closed on the disk
			for index_of_group in xrange(group_start, group_end):
				Clusterdir = os.path.join(Tracker["directory"], "Cluster%d"%index_of_group, "main%03d"%iteration)
				cfsc       = 0
				if Blockdata["myid_on_node"] == 0:					
					tvol0 		= get_im(os.path.join(Clusterdir, "tempdir", "tvol_0_%03d.hdf")%iteration)
					tweight0 	= get_im(os.path.join(Clusterdir, "tempdir", "tweight_0_%03d.hdf")%iteration)
					tvol1 		= get_im(os.path.join(Clusterdir, "tempdir", "tvol_1_%03d.hdf")%iteration)
					tweight1 	= get_im(os.path.join(Clusterdir, "tempdir", "tweight_1_%03d.hdf")%iteration)
					Util.fuse_low_freq(tvol0, tvol1, tweight0, tweight1, 2*Tracker["fuse_freq"])
					tag = 7007
					send_EMData(tvol1,    Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					send_EMData(tweight1, Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					tvol0.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1} )
					shrank0 	= stepone(tvol0, tweight0)
				elif Blockdata["myid_on_node"] == Blockdata["no_of_processes_per_group"]-1:					
					tag = 7007
					tvol1 		= recv_EMData(0, tag, Blockdata["shared_comm"])
					tweight1 	= recv_EMData(0, tag, Blockdata["shared_comm"])
					tvol1.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1} )
					shrank1 	= stepone(tvol1, tweight1)
					treg1       = get_im(os.path.join(Clusterdir, "tempdir", "trol_1_%03d.hdf"%iteration))
				mpi_barrier(Blockdata["shared_comm"])
				if 	Blockdata["myid_on_node"] == 0:
					tag = 7007					
					send_EMData(shrank0, Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					del shrank0
					lcfsc = 0
				elif Blockdata["myid_on_node"] == Blockdata["no_of_processes_per_group"]-1:
					tag         = 7007
					shrank0 	= recv_EMData(0, tag, Blockdata["shared_comm"])
					#  Note shrank volumes are Fourier uncentered.
					cfsc 		= fsc(shrank0, shrank1)[1]
					write_text_row(cfsc, os.path.join(Tracker["directory"], "fsc_driver_grp%03d_iter%03d.txt")%(index_of_group,iteration))
					del shrank0, shrank1
					if(Tracker["nxinit"]<Tracker["constants"]["nnxo"]):
						cfsc 	= cfsc[:Tracker["nxinit"]]
						for i in xrange(len(cfsc),Tracker["constants"]["nnxo"]//2+1):  cfsc.append(0.0)
					lcfsc  = len(cfsc)							
					fsc05  = 0
					fsc143 = 0 		
					for ifreq in xrange(len(cfsc)):	
						if cfsc[ifreq] <0.5: break
					fsc05  = ifreq - 1
					for ifreq in xrange(len(cfsc)):
						if cfsc[ifreq]<0.143: break
					fsc143 = ifreq - 1
					Tracker["fsc143"] = fsc143
					Tracker["fsc05"]  = fsc05					
					line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
					msg = "group %d  of FSC is done"%index_of_group
					print(line, msg)
					log_main.add(msg)
				Tracker = wrap_mpi_bcast(Tracker, Blockdata["no_of_processes_per_group"]-1, Blockdata["shared_comm"])
				cfsc = wrap_mpi_bcast(cfsc, Blockdata["no_of_processes_per_group"]-1, Blockdata["shared_comm"])
				Tracker["maxfrad"] = Tracker["constants"]["nnxo"]//2
				if( Blockdata["myid_on_node"] == 0):
					res_05[index_of_group]  = Tracker["fsc05"]
					res_143[index_of_group] = Tracker["fsc143"]
				if Blockdata["fftwmpi"]:
					if( Blockdata["myid_on_node"] == 0):
						treg0 = get_im(os.path.join(Clusterdir, "tempdir", "trol_0_%03d.hdf"%iteration))
					else:
						treg0    = model_blank(1)
						tvol0    = model_blank(1)
						tweight0 = model_blank(1)
				 	tvol0 = steptwo_mpi(tvol0, tweight0, treg0, cfsc, False, color = index_of_colors)
				 	if( Blockdata["myid_on_node"] == 0): 
				 		tvol0.write_image(os.path.join(masterdir, "vol_unfiltered_0_grp%03d.hdf"%index_of_group))
				 	del tvol0, tweight0, treg0
				else:
					if( Blockdata["myid_on_node"] == 0):
						treg0 = get_im(os.path.join(Clusterdir, "tempdir", "trol_0_%03d.hdf"%iteration))
						tvol0 = steptwo(tvol0, tweight0, treg0, cfsc, False)
						tvol0.write_image(os.path.join(masterdir, "vol_unfiltered_0_grp%03d.hdf"%index_of_group))			
						del tvol0, tweight0, treg0
				mpi_barrier(Blockdata["shared_comm"])
				if Blockdata["fftwmpi"]:
					if( Blockdata["myid_on_node"] == 0):# has to be main cpu
						tvol1 		= get_im(os.path.join(Clusterdir, "tempdir", "tvol_1_%03d.hdf")%iteration)
						tweight1 	= get_im(os.path.join(Clusterdir, "tempdir", "tweight_1_%03d.hdf")%iteration)
						treg1 =       get_im(os.path.join(Clusterdir, "tempdir", "trol_1_%03d.hdf"%iteration))
						tvol0 		= get_im(os.path.join(Clusterdir, "tempdir", "tvol_0_%03d.hdf")%iteration)
						tweight0 	= get_im(os.path.join(Clusterdir, "tempdir", "tweight_0_%03d.hdf")%iteration)
						Util.fuse_low_freq(tvol0, tvol1, tweight0, tweight1, 2*Tracker["fuse_freq"])
					else:
						treg1    = model_blank(1)
						tvol1    = model_blank(1)
						tweight1 = model_blank(1)
					tvol1 = steptwo_mpi(tvol1, tweight1, treg1, cfsc, False, color = index_of_colors)
					if (Blockdata["myid_on_node"] == 0): 
						tvol1.write_image(os.path.join(masterdir, "vol_unfiltered_1_grp%03d.hdf"%index_of_group))
					del tvol1, tweight1, treg1
				else:
					if( Blockdata["myid_on_node"] ==0):
						tvol1    = get_im(os.path.join(Clusterdir, "tempdir", "tvol_1_%03d.hdf")%iteration)
						tweight1 = get_im(os.path.join(Clusterdir, "tempdir", "tweight_1_%03d.hdf")%iteration)
						treg1    = get_im(os.path.join(Clusterdir, "tempdir", "trol_1_%03d.hdf"%iteration))
						tvol1    = steptwo(tvol1, tweight1, treg1, cfsc, False)
						tvol1.write_image(os.path.join(masterdir, "vol_unfiltered_1_grp%03d.hdf"%index_of_group))					
						#--  memory_check(Blockdata["myid"],"first node, before steptwo")
						del tvol1, tweight1, treg1
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
		for icluster in xrange(Tracker["number_of_groups"]):
			res05  = Tracker["constants"]["pixel_size"]*Tracker["constants"]["nnxo"]/res_05[icluster]
			res143 = Tracker["constants"]["pixel_size"]*Tracker["constants"]["nnxo"]/res_143[icluster]
			msg = "cluster  %d   fsc143/fsc05   %d/%d, %f/%f"%(icluster, Tracker["fsc143"][icluster], Tracker["fsc05"][icluster], res143, res05)
			print(msg)
			log_main.add(msg)
	keepgoing = bcast_number_to_all(keepgoing, source_node = Blockdata["main_node"], mpi_comm = MPI_COMM_WORLD) # always check	
	Tracker   = wrap_mpi_bcast(Tracker, Blockdata["main_node"])
	if (Blockdata["myid"] == Blockdata["main_node"]):
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		msg = " trl_iter do volumes end"
		print(line, msg)
		log_main.add(msg)
		msg = "command for sharpening maps: sxprocess.py map0.hdf map1.hdf --postprocess  --mask=%s --fsc_adj pixel_size=%f"%(Tracker["constants"]["mask3D"], Tracker["constants"]["pixel_size"])
		print(line, msg)
		log_main.add(msg)
	if not keepgoing: ERROR("do3d_sorting_groups_trl_iter  %s"%os.path.join(Tracker["directory"], "tempdir"),"do3d_sorting_groups_trl_iter", 1, Blockdata["myid"]) 
	return
####<<<--------
### nofsc rec3d	
def do3d_sorting_groups_nofsc_smearing_iter(srdata, partial_rec3d, iteration):
	global Tracker, Blockdata
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] == Blockdata["main_node"]: print(line, "---->>>do3d_sorting_groups_nofsc_smearing_iter<<<-----")
	keepgoing = 1
	if(Blockdata["myid"] == Blockdata["main_node"]):
		cmd = "{} {}".format("mkdir", os.path.join(Tracker["directory"], "tempdir"))
		if os.path.exists(os.path.join(Tracker["directory"], "tempdir")):print("tempdir exists")
		else: cmdexecute(cmd)
	mpi_barrier(MPI_COMM_WORLD)
	for index_of_groups in xrange(Tracker["number_of_groups"]):
		if partial_rec3d:
			msg = "---->>>partial rec3d<<<---"
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			if Blockdata["myid"] == Blockdata["main_node"]:print(line, msg)
			tvol, tweight, trol = recons3d_trl_struct_group_nofsc_shifted_data_partial_MPI(Blockdata["myid"], Blockdata["main_node"], Blockdata["nproc"], srdata, index_of_groups, \
			os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_groups), \
			os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf"%index_of_groups), \
			os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf"%index_of_groups),\
      		None,  Tracker["constants"]["CTF"], (2*Tracker["nxinit"]+3), Tracker["nosmearing"])
		else:
			msg = "---->>>full rec3d<<<---"
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			if Blockdata["myid"] == Blockdata["main_node"]:print(line, msg)
			tvol, tweight, trol = recons3d_trl_struct_group_nofsc_shifted_data_MPI(Blockdata["myid"], Blockdata["main_node"], srdata,\
			     index_of_groups, None,  Tracker["constants"]["CTF"], (2*Tracker["nxinit"]+3), Tracker["nosmearing"])
		if(Blockdata["myid"] == Blockdata["main_node"]):
			tvol.set_attr("is_complex",0)
			tvol.write_image(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf"%index_of_groups))
			tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf"%index_of_groups))
			trol.write_image(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_groups))
			del tvol
			del tweight
			del trol
			"""
			if partial_rec3d:
				tvol3.set_attr("is_complex",0)
				tvol3.write_image(os.path.join(Tracker["directory"], "tempdir", "tvol_3_%d.hdf"%index_of_groups))
				tweight3.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_3_%d.hdf"%index_of_groups))
				trol3.write_image(os.path.join(Tracker["directory"], "tempdir", "trol_3_%d.hdf"%index_of_groups))
			"""
		mpi_barrier(MPI_COMM_WORLD)
	mpi_barrier(MPI_COMM_WORLD)
	#fsc143            = 0
	#fsc05             = 0
	Tracker["fsc143"] = 0
	Tracker["fsc05"]  = 0
	Tracker["maxfrad"]= Tracker["nxinit"]//2
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] == Blockdata["main_node"]: print(line, "trl_iter starts to do volumes, number of groups is ", Tracker["number_of_groups"])
	if Blockdata["no_of_groups"]>1:
		for index_of_colors in xrange(Blockdata["no_of_groups"]):
			group_start, group_end = MPI_volume_start_end(Tracker["number_of_groups"], Blockdata["no_of_groups"], index_of_colors)
			if Blockdata["color"] == index_of_colors:  # It has to be 1 to avoid problem with tvol1 not closed on the disk
				for index_of_group in xrange(group_start, group_end):
					Tracker["fsc143"] =	Tracker["nxinit"]//2
					Tracker["fsc05"]  =	Tracker["nxinit"]//2
					cfsc = [0.0 for i in xrange(Tracker["constants"]["nnxo"])]
					if Blockdata["fftwmpi"]:
						if(Blockdata["myid_on_node"] == 0):
							tvol2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")%index_of_group)
							tweight2 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")%index_of_group)
							treg2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_group))
						else:
							tvol2 		= model_blank(1)
							tweight2 	= model_blank(1)
							treg2		= model_blank(1)
						tvol2 = steptwo_mpi(tvol2, tweight2, treg2, cfsc, False, color = index_of_colors) # has to be False!!!
						del tweight2, treg2
						if(Blockdata["myid_on_node"] == 0):
							if(Tracker["mask3D"] == None):tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
							else: Util.mul_img(tvol2, get_im(Tracker["constants"]["mask3D"]))
							tvol2.write_image(os.path.join(Tracker["directory"], "vol_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
							del tvol2
					else:
						if(Blockdata["myid_on_node"] == 0):
							tvol2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")%index_of_group)
							tweight2 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")%index_of_group)
							treg2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_group))
							tvol2 = steptwo(tvol2, tweight2, treg2, cfsc, False)
							del tweight2, treg2
							if(Tracker["mask3D"] == None): tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
							else: Util.mul_img(tvol2, get_im(Tracker["constants"]["mask3D"]))
							tvol2.write_image(os.path.join(Tracker["directory"], "vol_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
							del tvol2
					mpi_barrier(Blockdata["shared_comm"])
					######
					"""
					if partial_rec3d:
						if Blockdata["fftwmpi"]:
							if(Blockdata["myid_on_node"] == 0):
								tvol2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_3_%d.hdf")%index_of_group)
								tweight2 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_3_%d.hdf")%index_of_group)
								treg2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "trol_3_%d.hdf"%index_of_group))
							else:
								tvol2 		= model_blank(1)
								tweight2 	= model_blank(1)
								treg2		= model_blank(1)
							tvol2 = steptwo_mpi(tvol2, tweight2, treg2, cfsc, False, color = index_of_colors) # has to be False!!!
							del tweight2, treg2
							if(Blockdata["myid_on_node"] == 0):
								if(Tracker["mask3D"] == None):tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
								else: Util.mul_img(tvol2, get_im(Tracker["constants"]["mask3D"]))
								tvol2.write_image(os.path.join(Tracker["directory"], "vol_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
								del tvol2
						else:
							if(Blockdata["myid_on_node"] == 0):
								tvol2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_3_%d.hdf")%index_of_group)
								tweight2 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_3_%d.hdf")%index_of_group)
								treg2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "trol_3_%d.hdf"%index_of_group))
								tvol2       = steptwo(tvol2, tweight2, treg2, cfsc, False)
								del tweight2, treg2
								if(Tracker["mask3D"] == None):tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
								else: Util.mul_img(tvol2, get_im(Tracker["constants"]["mask3D"]))
								tvol2.write_image(os.path.join(Tracker["directory"], "vol_3_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
								del tvol2
						"""
				mpi_barrier(Blockdata["shared_comm"])
			mpi_barrier(Blockdata["shared_comm"])
	else:# loop over all groups for single node
		for index_of_group in xrange(Tracker["number_of_groups"]):
				Tracker["fsc143"] =	Tracker["nxinit"]//2
				Tracker["fsc05"]  =	Tracker["nxinit"]//2
				cfsc = [0.0 for i in xrange(Tracker["constants"]["nnxo"])]
				if Blockdata["fftwmpi"]:
					if(Blockdata["myid_on_node"] == 0):
						tvol2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")%index_of_group)
						tweight2 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")%index_of_group)
						treg2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_group))
						tvol2.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1})
					else:
						tvol2 		= model_blank(1)
						tweight2 	= model_blank(1)
						treg2		= model_blank(1)
					tvol2 = steptwo_mpi(tvol2, tweight2, treg2, cfsc, False, color = Blockdata["node_volume"][0]) # has to be False!!!
					del tweight2, treg2
					if(Blockdata["myid_on_node"] == 0):
						if(Tracker["mask3D"] == None):tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
						else: Util.mul_img(tvol2, get_im(Tracker["constants"]["mask3D"]))
						tvol2.write_image(os.path.join(Tracker["directory"], "vol_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
						del tvol2
				else:# to be paralleled
					if(Blockdata["myid_on_node"] == 0):
						tvol2 	 = get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")%index_of_group)
						tweight2 = get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")%index_of_group)
						treg2 	 = get_im(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_group))
						tvol2.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1})
						tvol2    = steptwo(tvol2, tweight2, treg2, cfsc, False)
						del tweight2, treg2
						if(Tracker["mask3D"] == None):tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
						else: Util.mul_img(tvol2, get_im(Tracker["constants"]["mask3D"]))
						tvol2.write_image(os.path.join(Tracker["directory"], "vol_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
						del tvol2
				mpi_barrier(MPI_COMM_WORLD)
	mpi_barrier(MPI_COMM_WORLD)
	keepgoing = bcast_number_to_all(keepgoing, source_node = Blockdata["main_node"], mpi_comm = MPI_COMM_WORLD) # always check 
	Tracker   = wrap_mpi_bcast(Tracker, Blockdata["main_node"])
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if (Blockdata["myid"] == Blockdata["main_node"]): print(line, " trl_iter do volumes end")
	if not keepgoing: ERROR("do3d_sorting_groups_trl_iter  %s"%os.path.join(Tracker["directory"], "tempdir"),"do3d_sorting_groups_trl_iter", 1, Blockdata["myid"])
	return
### nofsc insertion #1
def recons3d_trl_struct_group_nofsc_shifted_data_partial_MPI(myid, main_node, nproc, prjlist, group_ID, refvol_file, fftvol_file, weight_file, mpi_comm= None, CTF = True, target_size=-1, nosmearing = False):
	"""
		partial rec3d for re-assigned particles
		reconstructor nn4_ctfws
	"""
	from utilities    import reduce_EMData_to_root, random_string, get_im, findall, info, model_blank
	from EMAN2        import Reconstructors
	from filter	      import filt_table
	from fundamentals import fshift
	from mpi          import MPI_COMM_WORLD, mpi_barrier
	import types
	import datetime
	if mpi_comm == None: mpi_comm = MPI_COMM_WORLD
	if CTF: do_ctf = 1
	else:   do_ctf = 0
	if not os.path.exists(refvol_file):ERROR("refvol does not exist", "recons3d_trl_struct_group_nofsc_shifted_data_partial_MPI", 1, myid)
	if not os.path.exists(fftvol_file):ERROR("fftvol does not exist", "recons3d_trl_struct_group_nofsc_shifted_data_partial_MPI", 1, myid)
	if not os.path.exists(weight_file):ERROR("weight does not exist", "recons3d_trl_struct_group_nofsc_shifted_data_partial_MPI", 1, myid)
	
	#refvol
	if myid == main_node: target_size = get_im(refvol_file).get_xsize()
	else: target_size = 0
	target_size = bcast_number_to_all(target_size, main_node, mpi_comm)
	refvol = model_blank(target_size)# set to zero
	
	# fftvol
	if (myid == main_node): 
		fftvol = get_im(fftvol_file)
		fftvol /=float(nproc)
		#Util.mult_scalar(fftvol, 1./float(Blockdata["nproc"]))
		nxfft  =fftvol.get_xsize()
		nyfft  =fftvol.get_ysize()
		nzfft  =fftvol.get_zsize()
	else: 
		nxfft = 0
		nyfft = 0
		nzfft = 0
	nxfft = bcast_number_to_all(nxfft, main_node, mpi_comm)
	nyfft = bcast_number_to_all(nyfft, main_node, mpi_comm)
	nzfft = bcast_number_to_all(nzfft, main_node, mpi_comm)
	if (myid!= main_node): fftvol = model_blank(nxfft, nyfft, nzfft)
	bcast_EMData_to_all(fftvol, myid, main_node)
	
	# weight
	if (myid == main_node): 
		weight = get_im(weight_file)
		weight /=float(nproc)
		#Util.mult_scalar(weight, 1./float(Blockdata["nproc"]))
		nxweight  = weight.get_xsize()
		nyweight  = weight.get_ysize()
		nzweight  = weight.get_zsize()
	else: 
		nxweight = 0
		nyweight = 0
		nzweight = 0
	nxweight = bcast_number_to_all(nxweight, main_node, mpi_comm)
	nyweight = bcast_number_to_all(nyweight, main_node, mpi_comm)
	nzweight = bcast_number_to_all(nzweight, main_node, mpi_comm)
	if(myid != main_node): weight = model_blank(nxweight, nyweight, nzweight)
	bcast_EMData_to_all(weight, myid, main_node)
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = Reconstructors.get("nn4_ctfws", params)
	r.setup()
	if nosmearing:
		nnx = prjlist[0].get_xsize()
		nny = prjlist[0].get_ysize()
	else:
		nnx = prjlist[0][0].get_xsize()
		nny = prjlist[0][0].get_ysize()
	for im in xrange(len(prjlist)):
		if nosmearing:
			current_group_ID  = prjlist[im].get_attr("group")
			previous_group_ID = prjlist[im].get_attr("previous_group")
			if current_group_ID !=previous_group_ID:
				if current_group_ID == group_ID:
					flag = 1.0
					[phi, theta, psi, s2x, s2y] = get_params_proj(prjlist[im], xform = "xform.projection")
					r.insert_slice(prjlist[im], Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi}), flag)
				if previous_group_ID == group_ID: 
					flag = -1.0
					[phi, theta, psi, s2x, s2y] = get_params_proj(prjlist[im], xform = "xform.projection")
					r.insert_slice(prjlist[im], Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi}), flag)
		else:
			current_group_ID  = prjlist[im][0].get_attr("group")
			previous_group_ID = prjlist[im][0].get_attr("previous_group")
			if current_group_ID !=previous_group_ID:
				if current_group_ID == group_ID:
					flag = 1.0
					for jm in xrange(len(prjlist[im])):
						[phi, theta, psi, s2x, s2y] = get_params_proj(prjlist[im][jm], xform = "xform.projection")
						r.insert_slice(prjlist[im][jm], Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi}), prjlist[im][jm].get_attr("wprob")*flag)
				if previous_group_ID == group_ID: 
					flag =-1.0
					for jm in xrange(len(prjlist[im])):
						[phi, theta, psi, s2x, s2y] = get_params_proj(prjlist[im][jm], xform = "xform.projection")
						r.insert_slice(prjlist[im][jm], Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi}), prjlist[im][jm].get_attr("wprob")*flag)
	reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
	reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
	if myid == main_node:dummy = r.finish(True)
	mpi_barrier(mpi_comm)
	if myid == main_node: return fftvol, weight, refvol
	else:
		del fftvol
		del weight
		del refvol 
		return None, None, None
### insertion 2
def recons3d_trl_struct_group_nofsc_shifted_data_MPI(myid, main_node, prjlist, group_ID, mpi_comm= None, CTF = True, target_size=-1, nosmearing = False):
	"""
	  rec3d for pre-shifted data list
	  reconstructor nn4_ctfw
	"""
	from utilities    import reduce_EMData_to_root, random_string, get_im, findall, info, model_blank
	from EMAN2        import Reconstructors
	from filter	      import filt_table
	from fundamentals import fshift
	from mpi          import MPI_COMM_WORLD, mpi_barrier
	import types
	import datetime
	if mpi_comm == None: mpi_comm = MPI_COMM_WORLD
	refvol = model_blank(target_size)
	refvol.set_attr("fudge", 1.0)
	if CTF: do_ctf = 1
	else:   do_ctf = 0
	fftvol = EMData()
	weight = EMData()
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = Reconstructors.get( "nn4_ctfw", params)
	r.setup()
	if nosmearing:
		nnx = prjlist[0].get_xsize()
		nny = prjlist[0].get_ysize()
	else:
		nnx = prjlist[0][0].get_xsize()
		nny = prjlist[0][0].get_ysize()
	for im in xrange(len(prjlist)):
		if nosmearing: 
			if prjlist[im].get_attr("group") == group_ID:
				[phi, theta, psi, s2x, s2y] = get_params_proj(prjlist[im], xform = "xform.projection")
				r.insert_slice(prjlist[im], Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi}), 1.0)
		else: 
			if prjlist[im][0].get_attr("group") == group_ID:
				for jm in xrange(len(prjlist[im])):
					[phi, theta, psi, s2x, s2y] = get_params_proj(prjlist[im][jm], xform = "xform.projection")
					r.insert_slice(prjlist[im][jm], Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi}), prjlist[im][jm].get_attr("wprob"))
	reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
	reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
	if myid == main_node:dummy = r.finish(True)
	mpi_barrier(mpi_comm)
	if myid == main_node: return fftvol, weight, refvol
	else:
		del fftvol
		del weight
		del refvol 
		return None, None, None
###end of nofsc
def recons3d_trl_struct_group_MPI(myid, main_node, prjlist, random_subset, group_ID, paramstructure, norm_per_particle = None, \
      upweighted = True, mpi_comm= None, CTF = True, target_size=-1, nosmearing = False):
	"""
		recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
		Input
			list_of_prjlist: list of lists of projections to be included in the reconstruction
	"""
	global Tracker, Blockdata
	from utilities    import reduce_EMData_to_root, random_string, get_im, findall, model_blank, info
	from EMAN2        import Reconstructors
	from filter	      import filt_table
	from fundamentals import fshift
	from mpi          import MPI_COMM_WORLD, mpi_barrier
	import types
	import datetime
	import copy
	if mpi_comm == None: mpi_comm = MPI_COMM_WORLD
	
	refvol = model_blank(target_size)
	refvol.set_attr("fudge", 1.0)
	if CTF: do_ctf = 1
	else:   do_ctf = 0
	fftvol = EMData()
	weight = EMData()
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = Reconstructors.get( "nn4_ctfw", params )
	r.setup()
	if norm_per_particle == None: norm_per_particle = len(prjlist)*[1.0]
	if not nosmearing:
		delta  = Tracker["delta"]
		refang = Tracker["refang"]
		rshifts_shrank = copy.deepcopy(Tracker["rshifts"])
		nshifts = len(rshifts_shrank)
		for im in xrange(len(rshifts_shrank)):
			rshifts_shrank[im][0] *= float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
			rshifts_shrank[im][1] *= float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	nnx = prjlist[0].get_xsize()
	nny = prjlist[0].get_ysize()
	for im in xrange(len(prjlist)):
		if not nosmearing: avgnorm = Tracker["avgnorm"][prjlist[im].get_attr("chunk_id")]
		#  parse projection structure, generate three lists:
		#  [ipsi+iang], [ishift], [probability]
		#  Number of orientations for a given image
		if prjlist[im].get_attr("group") == group_ID:
			if random_subset == 2:
				if nosmearing:
					bckgn = prjlist[im].get_attr("bckgnoise")
					ct = prjlist[im].get_attr("ctf")
					prjlist[im].set_attr_dict( {"bckgnoise":bckgn, "ctf":ct} )
					[phi, theta, psi, s2x, s2y] = get_params_proj(prjlist[im], xform = "xform.projection")
					r.insert_slice(prjlist[im], Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi}), 1.0)
				else:
					numbor = len(paramstructure[im][2])
					ipsiandiang = [ paramstructure[im][2][i][0]/1000  for i in xrange(numbor) ]
					allshifts   = [ paramstructure[im][2][i][0]%1000  for i in xrange(numbor) ]
					probs       = [ paramstructure[im][2][i][1] for i in xrange(numbor) ]
					#  Find unique projection directions
					tdir = list(set(ipsiandiang))
					bckgn = prjlist[im].get_attr("bckgnoise")
					ct = prjlist[im].get_attr("ctf")
					#  For each unique projection direction:
					data = [None]*nshifts
					for ii in xrange(len(tdir)):
						#  Find the number of times given projection direction appears on the list, it is the number of different shifts associated with it.
						lshifts = findall(tdir[ii], ipsiandiang)
						toprab  = 0.0
						for ki in xrange(len(lshifts)):  toprab += probs[lshifts[ki]]
						recdata = EMData(nny,nny,1,False)
						recdata.set_attr("is_complex",0)
						for ki in xrange(len(lshifts)):
							lpt = allshifts[lshifts[ki]]
							if( data[lpt] == None ):
								data[lpt] = fshift(prjlist[im], rshifts_shrank[lpt][0], rshifts_shrank[lpt][1])
								data[lpt].set_attr("is_complex",0)
							Util.add_img(recdata, Util.mult_scalar(data[lpt], probs[lshifts[ki]]/toprab))
						recdata.set_attr_dict({"padffted":1, "is_fftpad":1,"is_fftodd":0, "is_complex_ri":1, "is_complex":1})
						if not upweighted:  recdata = filt_table(recdata, bckgn )
						recdata.set_attr_dict( {"bckgnoise":bckgn, "ctf":ct} )
						ipsi = tdir[ii]%100000
						iang = tdir[ii]/100000
						#for iloop in xrange(10000000):
						#if iloop%1000==0:memory_check("before slice %d  myid  %d"%(iloop, Blockdata["myid"]))
						r.insert_slice( recdata, Transform({"type":"spider","phi":refang[iang][0],"theta":refang[iang][1],"psi":refang[iang][2]+ipsi*delta}), toprab*avgnorm/norm_per_particle[im])
						#if iloop%1000==0:memory_check("after slice %d  myid  %d"%(iloop, Blockdata["myid"]))
			else:
				if	prjlist[im].get_attr("chunk_id") == random_subset:
					if nosmearing:
						bckgn = prjlist[im].get_attr("bckgnoise")
						ct = prjlist[im].get_attr("ctf")
						prjlist[im].set_attr_dict({"bckgnoise":bckgn, "ctf":ct})
						[phi, theta, psi, s2x, s2y] = get_params_proj(prjlist[im], xform = "xform.projection")
						r.insert_slice(prjlist[im], Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi}), 1.0)
					else:
						numbor = len(paramstructure[im][2])
						ipsiandiang = [ paramstructure[im][2][i][0]/1000  for i in xrange(numbor) ]
						allshifts   = [ paramstructure[im][2][i][0]%1000  for i in xrange(numbor) ]
						probs       = [ paramstructure[im][2][i][1] for i in xrange(numbor) ]
						#  Find unique projection directions
						tdir = list(set(ipsiandiang))
						bckgn = prjlist[im].get_attr("bckgnoise")
						ct = prjlist[im].get_attr("ctf")
						#  For each unique projection direction:
						data = [None]*nshifts
						for ii in xrange(len(tdir)):
							#  Find the number of times given projection direction appears on the list, it is the number of different shifts associated with it.
							lshifts = findall(tdir[ii], ipsiandiang)
							toprab  = 0.0
							for ki in xrange(len(lshifts)):  toprab += probs[lshifts[ki]]
							recdata = EMData(nny,nny,1,False)
							recdata.set_attr("is_complex",0)
							for ki in xrange(len(lshifts)):
								lpt = allshifts[lshifts[ki]]
								if( data[lpt] == None ):
									data[lpt] = fshift(prjlist[im], rshifts_shrank[lpt][0], rshifts_shrank[lpt][1])
									data[lpt].set_attr("is_complex",0)
								Util.add_img(recdata, Util.mult_scalar(data[lpt], probs[lshifts[ki]]/toprab))
							recdata.set_attr_dict({"padffted":1, "is_fftpad":1,"is_fftodd":0, "is_complex_ri":1, "is_complex":1})
							if not upweighted:  recdata = filt_table(recdata, bckgn )
							recdata.set_attr_dict( {"bckgnoise":bckgn, "ctf":ct} )
							ipsi = tdir[ii]%100000
							iang = tdir[ii]/100000
							r.insert_slice(recdata, Transform({"type":"spider","phi":refang[iang][0],"theta":refang[iang][1],"psi":refang[iang][2]+ipsi*delta}), toprab*avgnorm/norm_per_particle[im])
	#  clean stuff
	#if not nosmearing: del recdata, tdir, ipsiandiang, allshifts, probs
	reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
	reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
	if not nosmearing: del rshifts_shrank
	if myid == main_node:dummy = r.finish(True)
	mpi_barrier(mpi_comm)
	if myid == main_node: return fftvol, weight, refvol
	else:
		del fftvol
		del weight
		del refvol 
		return None, None, None
####<<<<--------
#####  FSC rec3d
def do3d_sorting_groups_fsc_only_iter(data, paramstructure, norm_per_particle, iteration):
	global Tracker, Blockdata
	# do resolution each time
	keepgoing = 1
	if(Blockdata["myid"] == Blockdata["nodes"][0]):
		cmd = "{} {}".format("mkdir", os.path.join(Tracker["directory"], "tempdir"))
		if os.path.exists(os.path.join(Tracker["directory"], "tempdir")): print("tempdir exists")
		else: cmdexecute(cmd)
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] == Blockdata["main_node"]:print(line, "insertion starts")
	do3d_sorting_group_insertion_random_two_for_fsc(data, paramstructure, norm_per_particle)
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] == Blockdata["main_node"]:print(line, "insertion finishes")
	mpi_barrier(MPI_COMM_WORLD)
	fsc143             = 0
	fsc05              = 0
	Tracker["fsc143"]  = 0
	Tracker["fsc05"]   = 0
	res_05             = Tracker["number_of_groups"]*[0]
	res_143            = Tracker["number_of_groups"]*[0]
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] == Blockdata["main_node"]: print(line, "doing volumes starts, and the total number of groups is ",Tracker["number_of_groups"])
	for index_of_colors in xrange(Blockdata["no_of_groups"]):
		group_start, group_end = MPI_volume_start_end(Tracker["number_of_groups"], Blockdata["no_of_groups"], index_of_colors)
		if Blockdata["color"] == index_of_colors:  # It has to be 1 to avoid problem with tvol1 not closed on the disk
			for index_of_group in xrange(group_start, group_end):								
				if Blockdata["myid_on_node"] == 0:
					#print(" odd   group    %d"%index_of_group)	
					tvol0 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0_0_%d.hdf")%index_of_group)
					tweight0 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_0_0_%d.hdf")%index_of_group)
					tvol1 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_1_0_%d.hdf")%index_of_group)
					tweight1 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_1_0_%d.hdf")%index_of_group)
					Util.fuse_low_freq(tvol0, tvol1, tweight0, tweight1, 2*Tracker["fuse_freq"])
					tag = 7007
					send_EMData(tvol1, Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					send_EMData(tweight1, Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					shrank0 	= stepone(tvol0, tweight0)	
				elif Blockdata["myid_on_node"] == Blockdata["no_of_processes_per_group"]-1:
					#print(" odd   group    %d"%index_of_group)	
					tag = 7007
					tvol1 		= recv_EMData(0, tag, Blockdata["shared_comm"])
					tweight1 	= recv_EMData(0, tag, Blockdata["shared_comm"])
					tvol1.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1} )
					shrank1 	= stepone(tvol1, tweight1)
				if Blockdata["myid_on_node"] == 1:
					tvol0 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0_1_%d.hdf")%index_of_group)
					tweight0 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_0_1_%d.hdf")%index_of_group)
					tvol1 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_1_1_%d.hdf")%index_of_group)
					tweight1 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_1_1_%d.hdf")%index_of_group)
					Util.fuse_low_freq(tvol0, tvol1, tweight0, tweight1, 2*Tracker["fuse_freq"])
					tag = 7007
					send_EMData(tvol1, Blockdata["no_of_processes_per_group"]-2, tag, Blockdata["shared_comm"])
					send_EMData(tweight1, Blockdata["no_of_processes_per_group"]-2, tag, Blockdata["shared_comm"])
					shrank0 	= stepone(tvol0, tweight0)					
				elif Blockdata["myid_on_node"] == Blockdata["no_of_processes_per_group"]-2:
					#print(" even   group    %d"%index_of_group)
					tag = 7007
					tvol1 		= recv_EMData(1, tag, Blockdata["shared_comm"])
					tweight1 	= recv_EMData(1, tag, Blockdata["shared_comm"])
					tvol1.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1} )
					shrank1 	= stepone(tvol1, tweight1)					
				mpi_barrier(Blockdata["shared_comm"])				
				if 	Blockdata["myid_on_node"] == 0:
					tag = 7007					
					send_EMData(shrank0, Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					del shrank0
					lcfsc = 0
				elif Blockdata["myid_on_node"] == Blockdata["no_of_processes_per_group"]-1:
					#print(" now we do fsc  odd ")	
					tag = 7007
					shrank0 	= recv_EMData(0, tag, Blockdata["shared_comm"])
					cfsc 		= fsc(shrank0, shrank1)[1]
					write_text_row(cfsc, os.path.join(Tracker["directory"], "fsc_driver_chunk0_grp%03d_iter%03d.txt")%(index_of_group,iteration))
					del shrank0, shrank1
					if(Tracker["nxinit"]<Tracker["constants"]["nnxo"]):
						cfsc 	= cfsc[:Tracker["nxinit"]]
						for i in xrange(len(cfsc),Tracker["constants"]["nnxo"]//2+1):  cfsc.append(0.0)
					lcfsc = len(cfsc)							
					fsc05  = 0
					fsc143 = 0 
					for ifreq in xrange(len(cfsc)):	
						if cfsc[ifreq] <0.5: break
					fsc05  = ifreq - 1
					for ifreq in xrange(len(cfsc)):
						if cfsc[ifreq]<0.143: break
					fsc143 = ifreq - 1
					Tracker["fsc143"] = fsc143
					Tracker["fsc05"] = fsc05
					line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
					print(line, "group %d  of do3d_sorting_groups_trl_iter is done"%index_of_group)					
				if 	Blockdata["myid_on_node"] == 1:
					#print(" now we do step one  even")	
					tag = 7007					
					send_EMData(shrank0, Blockdata["no_of_processes_per_group"]-2, tag, Blockdata["shared_comm"])
					del shrank0
					lcfsc = 0					
				elif Blockdata["myid_on_node"] == Blockdata["no_of_processes_per_group"]-2:
					tag = 7007
					shrank0 = recv_EMData(1, tag, Blockdata["shared_comm"])
					cfsc = fsc(shrank0, shrank1)[1]
					write_text_row(cfsc, os.path.join(Tracker["directory"], "fsc_driver_chunk1_grp%03d_iter%03d.txt")%(index_of_group,iteration))
					del shrank0, shrank1
					if(Tracker["nxinit"]<Tracker["constants"]["nnxo"]):
						cfsc = cfsc[:Tracker["nxinit"]]
						for i in xrange(len(cfsc),Tracker["constants"]["nnxo"]//2+1):cfsc.append(0.0)
					lcfsc = len(cfsc)							
					fsc05  = 0
					fsc143 = 0
					for ifreq in xrange(len(cfsc)):	
						if cfsc[ifreq] <0.5: break
					fsc05  = ifreq - 1
					for ifreq in xrange(len(cfsc)):
						if cfsc[ifreq]<0.143: break
					fsc143 = ifreq - 1
					Tracker["fsc143"] = fsc143
					Tracker["fsc05"] = fsc05
					line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"	
				Tracker = wrap_mpi_bcast(Tracker, Blockdata["no_of_processes_per_group"]-1, Blockdata["shared_comm"])
				if( Blockdata["myid_on_node"] == 0):
					res_05[index_of_group]  = Tracker["fsc05"]
					res_143[index_of_group] = Tracker["fsc143"]				
				mpi_barrier(Blockdata["shared_comm"])
			mpi_barrier(Blockdata["shared_comm"])
	mpi_barrier(MPI_COMM_WORLD)
	keepgoing = bcast_number_to_all(keepgoing, source_node = Blockdata["main_node"], mpi_comm = MPI_COMM_WORLD) # always check 
	Tracker   = wrap_mpi_bcast(Tracker, Blockdata["main_node"])
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if (Blockdata["myid"] == Blockdata["main_node"]): print(line, " doing volumes finishes")
	if not keepgoing:ERROR("do3d_sorting_groups_trl_iter  %s"%os.path.join(Tracker["directory"], "tempdir"),"do3d_sorting_groups_trl_iter", 1, Blockdata["myid"]) 
	return
	
def do3d_sorting_group_insertion_random_two_for_fsc(data, sparamstructure, snorm_per_particle):
	global Tracker, Blockdata
	if(Blockdata["myid"] == Blockdata["nodes"][0]):
		cmd = "{} {}".format("mkdir", os.path.join(Tracker["directory"], "tempdir"))
		if os.path.exists(os.path.join(Tracker["directory"], "tempdir")):print("tempdir exists")
		else: cmdexecute(cmd)		
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] ==  Blockdata["main_node"]: print(line, "backprojection of total %d volumes starts"%Tracker["number_of_groups"])		
	for index_of_groups in xrange(Tracker["number_of_groups"]):	
		for procid in xrange(2):
			for ifsc in xrange(2):
				tvol, tweight, trol = recons3d_4nnsorting_group_fsc_MPI(myid = Blockdata["myid"], main_node = Blockdata["nodes"][procid],  \
				  prjlist = data, fsc_half = ifsc, random_subset = procid, group_ID = index_of_groups, paramstructure=sparamstructure, norm_per_particle=snorm_per_particle,\
				  CTF = Tracker["constants"]["CTF"], upweighted = False, target_size = (2*Tracker["nxinit"]+3))
				if(Blockdata["myid"] == Blockdata["nodes"][procid]):
					tvol.set_attr("is_complex",0)
					tvol.write_image(os.path.join(Tracker["directory"], "tempdir", "tvol_%d_%d_%d.hdf"%(ifsc, procid, index_of_groups)))
					tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%d_%d_%d.hdf"%(ifsc, procid, index_of_groups)))
					trol.write_image(os.path.join(Tracker["directory"], "tempdir", "trol_%d_%d_%d.hdf"%(ifsc, procid, index_of_groups)))
				mpi_barrier(MPI_COMM_WORLD)
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] ==  Blockdata["main_node"]: print(line, "backprojection of total %d  volumes finishes"%Tracker["number_of_groups"])		
	mpi_barrier(MPI_COMM_WORLD)
	return
		
def recons3d_4nnsorting_group_fsc_MPI(myid, main_node, prjlist, fsc_half, random_subset, group_ID, paramstructure, norm_per_particle, \
    CTF = True, upweighted = True, mpi_comm= None, target_size=-1):
	##      with smearing
	#####	recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
	####	Input
	####	list_of_prjlist: list of lists of projections to be included in the reconstruction
	global Tracker, Blockdata
	from utilities      import reduce_EMData_to_root, random_string, get_im, findall, model_blank, info, get_params_proj
	from EMAN2          import Reconstructors
	from filter		    import filt_table
	from mpi            import MPI_COMM_WORLD, mpi_barrier
	from statistics     import fsc 
	from reconstruction import insert_slices_pdf
	from fundamentals   import fft
	import datetime, types
	import copy
	if mpi_comm == None: mpi_comm = MPI_COMM_WORLD
	imgsize = prjlist[0].get_ysize()  # It can be Fourier, so take y-size
	refvol = model_blank(target_size)
	refvol.set_attr("fudge", 1.0)
	if CTF: do_ctf = 1
	else:   do_ctf = 0
	fftvol = EMData()
	weight = EMData()
	try:    qt = projlist[0].get_attr("qt")
	except: qt = 1.0
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r      = Reconstructors.get( "nn4_ctfw", params )
	r.setup()
	if norm_per_particle == None: norm_per_particle = len(prjlist)*[1.0]
	# Definitions for smearing ,all copied from refinement
	if not Tracker["nosmearing"]:
		delta         = Tracker["delta"]
		refang        = Tracker["refang"]
		rshifts_shrank  = copy.deepcopy(Tracker["rshifts"])
		nshifts = len(rshifts_shrank)
		for im in xrange(nshifts):
			rshifts_shrank[im][0] *= float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
			rshifts_shrank[im][1] *= float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	nnx = prjlist[0].get_xsize()
	nny = prjlist[0].get_ysize()
	nc  = 0
	for im in xrange(len(prjlist)):
		if prjlist[im].get_attr("group") == group_ID and prjlist[im].get_attr("chunk_id") == random_subset:
			if Tracker["nosmearing"]: avgnorm = 1.0
			else: avgnorm =  Tracker["avgnorm"][prjlist[im].get_attr("chunk_id")]#
			if nc %2 == fsc_half:
				if Tracker["nosmearing"]:
					ct = prjlist[im].get_attr("ctf")
					bckgn = prjlist[im].get_attr("bckgnoise")
					if not upweighted:  prjlist[im] = filt_table(prjlist[im], bckgn)
					prjlist[im].set_attr_dict( {"bckgnoise":bckgn, "ctf":ct})
					phi,theta,psi,s2x,s2y = get_params_proj(prjlist[im], xform = "xform.projection")
					r.insert_slice(prjlist[im], Transform({"type":"spider","phi":phi, "theta":theta, "psi":psi}), 1.0)
				else:
					numbor = len(paramstructure[im][2])
					ipsiandiang = [paramstructure[im][2][i][0]/1000  for i in xrange(numbor)]
					allshifts   = [paramstructure[im][2][i][0]%1000  for i in xrange(numbor)]
					probs       = [paramstructure[im][2][i][1] for i in xrange(numbor)]
					#  Find unique projection directions
					tdir = list(set(ipsiandiang))
					bckgn = prjlist[im].get_attr("bckgnoise")
					ct = prjlist[im].get_attr("ctf")
					#  For each unique projection direction:
					data = [None]*nshifts
					for ii in xrange(len(tdir)):
						#  Find the number of times given projection direction appears on the list, it is the number of different shifts associated with it.
						lshifts = findall(tdir[ii], ipsiandiang)
						toprab  = 0.0
						for ki in xrange(len(lshifts)):  toprab += probs[lshifts[ki]]
						recdata = EMData(nny,nny,1,False)
						recdata.set_attr("is_complex",0)
						for ki in xrange(len(lshifts)):
							lpt = allshifts[lshifts[ki]]
							if( data[lpt] == None ):
								data[lpt] = fshift(prjlist[im], rshifts_shrank[lpt][0], rshifts_shrank[lpt][1])
								data[lpt].set_attr("is_complex",0)
							Util.add_img(recdata, Util.mult_scalar(data[lpt], probs[lshifts[ki]]/toprab))
						recdata.set_attr_dict({"padffted":1, "is_fftpad":1,"is_fftodd":0, "is_complex_ri":1, "is_complex":1})
						if not upweighted:  recdata = filt_table(recdata, bckgn )
						recdata.set_attr_dict( {"bckgnoise":bckgn, "ctf":ct} )
						ipsi = tdir[ii]%100000
						iang = tdir[ii]/100000
						r.insert_slice( recdata, Transform({"type":"spider","phi":refang[iang][0],"theta":refang[iang][1],"psi":refang[iang][2]+ipsi*delta}), toprab*avgnorm/norm_per_particle[im])
			nc +=1
	reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
	reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
	if myid == main_node: dummy = r.finish(True)
	mpi_barrier(mpi_comm)
	if myid == main_node: return fftvol, weight, refvol
	else: return None, None, None
#####end of FSC
###<<<-----group rec3d
### insertion
def do3d_sorting_group_insertion_smearing(sdata, sparamstructure, snorm_per_particle, randomset=2):
	global Tracker, Blockdata
	if(Blockdata["myid"] == Blockdata["nodes"][0]):
		cmd = "{} {}".format("mkdir", os.path.join(Tracker["directory"], "tempdir"))
		if os.path.exists(os.path.join(Tracker["directory"], "tempdir")):print("tempdir exists")
		else: cmdexecute(cmd)
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] ==  Blockdata["main_node"]:print(line, "start backprojection of %d volumes"%Tracker["number_of_groups"])
	if randomset ==1:
		for index_of_groups in xrange(Tracker["number_of_groups"]):
			for procid in xrange(2, 3):
				tvol, tweight, trol = recons3d_trl_struct_group_MPI(myid = Blockdata["myid"], main_node = Blockdata["nodes"][procid//2],\
				  prjlist = sdata,  random_subset = procid, group_ID = index_of_groups, paramstructure = sparamstructure, \
				  norm_per_particle = snorm_per_particle, CTF = Tracker["constants"]["CTF"],\
					mpi_comm = None, upweighted = False, target_size = (2*Tracker["nxinit"]+3), nosmearing = Tracker["nosmearing"])			
				if(Blockdata["myid"] == Blockdata["nodes"][procid//2]):
					tvol.set_attr("is_complex",0)
					tvol.write_image(os.path.join(Tracker["directory"], "tempdir", "tvol_%d_%d.hdf"%(procid, index_of_groups)))
					tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%d_%d.hdf"%(procid, index_of_groups)))
					trol.write_image(os.path.join(Tracker["directory"], "tempdir", "trol_%d_%d.hdf"%(procid, index_of_groups)))
				del tvol
				del tweight
				del trol
				mpi_barrier(MPI_COMM_WORLD)
	else:		
		for index_of_groups in xrange(Tracker["number_of_groups"]):
			for procid in xrange(3):	
				tvol, tweight, trol = recons3d_trl_struct_group_MPI(myid = Blockdata["myid"], main_node = Blockdata["nodes"][procid//2],\
				  prjlist = sdata,  random_subset = procid, group_ID = index_of_groups, paramstructure = sparamstructure, \
				   norm_per_particle = snorm_per_particle, CTF = Tracker["constants"]["CTF"],\
					mpi_comm= None, upweighted = False, target_size = (2*Tracker["nxinit"]+3), nosmearing = Tracker["nosmearing"])
				if(Blockdata["myid"] == Blockdata["nodes"][procid//2]):
					tvol.set_attr("is_complex",0)
					tvol.write_image(os.path.join(Tracker["directory"], "tempdir", "tvol_%d_%d.hdf"%(procid, index_of_groups)))
					tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%d_%d.hdf"%(procid, index_of_groups)))
					trol.write_image(os.path.join(Tracker["directory"], "tempdir", "trol_%d_%d.hdf"%(procid, index_of_groups)))
				del tvol
				del tweight
				del trol
				mpi_barrier(MPI_COMM_WORLD)
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] ==  Blockdata["main_node"]:print(line, "end backprojection of %d volumes"%Tracker["number_of_groups"])		
	mpi_barrier(MPI_COMM_WORLD)
	return
### rec3d 
def do3d_sorting_groups_trl_smearing_iter(data, paramstructure, norm_per_particle, iteration, unfiltered = False):
	global Tracker, Blockdata
	keepgoing = 1
	if(Blockdata["myid"] == Blockdata["nodes"][0]):
		cmd = "{} {}".format("mkdir", os.path.join(Tracker["directory"], "tempdir"))
		if os.path.exists(os.path.join(Tracker["directory"], "tempdir")): print("tempdir exists")
		else: cmdexecute(cmd)
	do3d_sorting_group_insertion_smearing(data, paramstructure, norm_per_particle)
	mpi_barrier(MPI_COMM_WORLD)
	fsc143                = 0
	fsc05                 = 0
	Tracker["fsc143"]	  = 0
	Tracker["fsc05"]	  = 0
	res_05 				  = Tracker["number_of_groups"]*[0]
	res_143 			  = Tracker["number_of_groups"]*[0]	
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if Blockdata["myid"] == Blockdata["main_node"]: print(line, "trl_iter starts to do volumes, number of groups is ", Tracker["number_of_groups"])
	for index_of_colors in xrange(Blockdata["no_of_groups"]):
		group_start, group_end = MPI_volume_start_end(Tracker["number_of_groups"], Blockdata["no_of_groups"], index_of_colors)
		if Blockdata["color"] == index_of_colors:  # It has to be 1 to avoid problem with tvol1 not closed on the disk
			for index_of_group in xrange(group_start, group_end):
				cfsc = 0
				if Blockdata["myid_on_node"] == 0:
					tvol0 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0_%d.hdf")%index_of_group)
					tweight0 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_0_%d.hdf")%index_of_group)
					tvol1 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_1_%d.hdf")%index_of_group)
					tweight1 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_1_%d.hdf")%index_of_group)
					Util.fuse_low_freq(tvol0, tvol1, tweight0, tweight1, 2*Tracker["fuse_freq"])
					tag = 7007
					send_EMData(tvol1, Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					send_EMData(tweight1, Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					shrank0 	= stepone(tvol0, tweight0)
					
				elif Blockdata["myid_on_node"] == Blockdata["no_of_processes_per_group"]-1:
					tag = 7007
					tvol1 		= recv_EMData(0, tag, Blockdata["shared_comm"])
					tweight1 	= recv_EMData(0, tag, Blockdata["shared_comm"])
					tvol1.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1} )
					shrank1 	= stepone(tvol1, tweight1)
				mpi_barrier(Blockdata["shared_comm"])
				if 	Blockdata["myid_on_node"] == 0:
					tag = 7007					
					send_EMData(shrank0, Blockdata["no_of_processes_per_group"]-1, tag, Blockdata["shared_comm"])
					del shrank0
					lcfsc = 0
				elif Blockdata["myid_on_node"] == Blockdata["no_of_processes_per_group"]-1:
					tag = 7007
					shrank0 	= recv_EMData(0, tag, Blockdata["shared_comm"])
					cfsc 		= fsc(shrank0, shrank1)[1]
					write_text_row(cfsc, os.path.join(Tracker["directory"], "fsc_driver_grp%03d_iter%03d.txt")%(index_of_group,iteration))
					del shrank0, shrank1
					if(Tracker["nxinit"]<Tracker["constants"]["nnxo"]):
						cfsc = cfsc[:Tracker["nxinit"]]
						for i in xrange(len(cfsc),Tracker["constants"]["nnxo"]//2+1):  cfsc.append(0.0)
					lcfsc = len(cfsc)							
					fsc05  = 0
					fsc143 = 0 
					for ifreq in xrange(len(cfsc)):	
						if cfsc[ifreq] <0.5: break
					fsc05  = ifreq - 1
					for ifreq in xrange(len(cfsc)):
						if cfsc[ifreq]<0.143: break
					fsc143 = ifreq - 1
					Tracker["fsc143"] = fsc143
					Tracker["fsc05"] = fsc05
					line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
					print(line, "group %d  of do3d_sorting_groups_trl_iter is done"%index_of_group)
				Tracker = wrap_mpi_bcast(Tracker, Blockdata["no_of_processes_per_group"]-1, Blockdata["shared_comm"])
				cfsc = wrap_mpi_bcast(cfsc, Blockdata["no_of_processes_per_group"]-1, Blockdata["shared_comm"])
				Tracker["maxfrad"] = Tracker["nxinit"]//2
				res_05[index_of_group]  = Tracker["fsc05"]
				res_143[index_of_group] = Tracker["fsc143"]
				if Blockdata["fftwmpi"]:			
					if( Blockdata["myid_on_node"] == 0):
						tvol2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")%index_of_group)
						tweight2 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")%index_of_group)
						treg2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_group))
					else:
						tvol2 		= model_blank(1)
						tweight2 	= model_blank(1)
						treg2		= model_blank(1)
					tvol2 = steptwo_mpi(tvol2, tweight2, treg2, cfsc, unfiltered, color = index_of_colors) # has to be False!!!
					del tweight2, treg2
					if( Blockdata["myid_on_node"] == 0):
						if(Tracker["mask3D"] == None):  tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
						else: Util.mul_img(tvol2, get_im(Tracker["constants"]["mask3D"]))
						tvol2.write_image(os.path.join(Tracker["directory"], "vol_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
						del tvol2
				else:
					if( Blockdata["myid_on_node"] == 0):
						tvol2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")%index_of_group)
						tweight2 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")%index_of_group)
						treg2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_group))
						tvol2       = steptwo(tvol2, tweight2, treg2, cfsc, unfiltered)
						del tweight2, treg2
						if(Tracker["mask3D"] == None):  tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
						else: Util.mul_img(tvol2, get_im(Tracker["constants"]["mask3D"]))
						tvol2.write_image(os.path.join(Tracker["directory"], "vol_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
						del tvol2
				mpi_barrier(Blockdata["shared_comm"])
			mpi_barrier(Blockdata["shared_comm"])
	mpi_barrier(MPI_COMM_WORLD)	
	res_05  = mpi_reduce(res_05, Tracker["number_of_groups"], MPI_INT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	res_143 = mpi_reduce(res_143,Tracker["number_of_groups"], MPI_INT, MPI_SUM, Blockdata["main_node"], MPI_COMM_WORLD)
	res_05  = map(int, res_05)
	res_143 = map(int, res_143)	
	if (Blockdata["myid"] == Blockdata["main_node"]):
		Tracker["fsc143"] = res_143
		Tracker["fsc05"]  = res_05
	keepgoing = bcast_number_to_all(keepgoing, source_node = Blockdata["main_node"], mpi_comm = MPI_COMM_WORLD) # always check 
	Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"])
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if (Blockdata["myid"] == Blockdata["main_node"]):print(line, " trl_iter do volumes end")
	if not keepgoing:ERROR("do3d_sorting_groups_trl_iter  %s"%os.path.join(Tracker["directory"], "tempdir"),"do3d_sorting_groups_trl_iter", 1, Blockdata["myid"]) 
	return
####<<<-------MEM related functions
def _VmB(VmKey):
    global _proc_status, _scale
     # get pseudo file  /proc/<pid>/status
    try:
        t = open(_proc_status)
        v = t.read()
        t.close()
    except:
        return 0.0  # non-Linux?
     # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    i = v.index(VmKey)
    v = v[i:].split(None, 3)  # whitespace
    if len(v) < 3:
        return 0.0  # invalid format?
     # convert Vm value to bytes
    return float(v[1]) * _scale[v[2]]
def memory(since=0.0):
    '''Return memory usage in bytes.
    '''
    return _VmB('VmSize:') - since

def resident(since=0.0):
    '''Return resident memory usage in bytes.
    '''
    return _VmB('VmRSS:') - since
def stacksize(since=0.0):
    '''Return stack size in bytes.
    '''
    return _VmB('VmStk:') - since
def memory_check(s="check_memory"):
	import os
	print(s)
	print(s +"  memory ",  memory()/1.e9)
	print(s +" resident  ", resident()/1.e9)
	print(s +" stacksize ", stacksize()/1.e9)
####<<<----do final maps ---->>>
def do_final_maps(number_of_groups, minimum_size, selected_iter, refinement_dir, masterdir, rec3d_image_size, log_main):
	global Tracker, Blockdata
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if(Blockdata["myid"] == Blockdata["main_node"]):
		msg =" reconstruct non-filtered maps of original image size"
		print(line, msg)
		log_main.add(msg)
	# estimate memory
	if( Blockdata["myid"] == Blockdata["main_node"]):
		msg = "------->>>>>>>Check memory for do final maps<<<<----------"
		log_main.add(msg)
		print(line, msg)
	for icluster  in xrange(number_of_groups):
		clusterdir = os.path.join(masterdir, "Cluster%d"%icluster)
		if os.path.exists(clusterdir):
			if Blockdata["myid"] ==icluster:
				cmd="{} {}".format("rm -rf", clusterdir)
				cmdexecute(cmd)
	mpi_barrier(MPI_COMM_WORLD)
	basic_memory_per_cpu = 1.0
	total_data_in_mem = Tracker["constants"]["nnxo"]*Tracker["constants"]["nnxo"]*Tracker["constants"]["total_stack"]*4./1.e9
	one_volume_in_mem = Tracker["constants"]["nnxo"]*Tracker["constants"]["nnxo"]*Tracker["constants"]["nnxo"]*4.*8./1.e9
	nproc_do_final_per_node =(Tracker["constants"]["memory_per_node"] - total_data_in_mem -1.0)/(basic_memory_per_cpu + one_volume_in_mem)
	if( Blockdata["myid"] == Blockdata["main_node"]):
		msg = "total mem per node: %5.1f G"%Tracker["constants"]["memory_per_node"]
		log_main.add(msg)
		print(line, msg)
	nproc_do_final_per_node = int(nproc_do_final_per_node)
	if nproc_do_final_per_node > Blockdata["nproc"] //Blockdata["no_of_groups"]:
		nproc_do_final_per_node = Blockdata["nproc"] //Blockdata["no_of_groups"]
	if Blockdata["nproc_previous"] > 0: nproc_do_final_per_node = min(nproc_do_final_per_node, Blockdata["nproc_previous"]//Blockdata["no_of_groups"])
	ncpu_per_node = min(minimum_size//5//Blockdata["no_of_groups"]//2, nproc_do_final_per_node)
	ncpu_per_node = max(ncpu_per_node, 2)
	if( Blockdata["myid"] == Blockdata["main_node"]):
		msg = "determined CPU used per node: %d"%ncpu_per_node
		log_main.add(msg)
		print(line, msg)
	Blockdata["ncpuspernode"] = ncpu_per_node
	Blockdata["nsubset"] = Blockdata["ncpuspernode"]*Blockdata["no_of_groups"]
	#print(Blockdata["nsubset"])
	create_subgroup()
	#print(Blockdata["nsubset"], Blockdata["subgroup_myid"], Blockdata["myid"])
	fuse_freq = Tracker["fuse_freq"] # sort does it already
	mask3D    = Tracker["mask3D"]
	mtf       = Tracker["constants"]["mtf"]
	fsc_adj   = Tracker["constants"]["fsc_adj"]
	Bstart    = Tracker["constants"]["B_start"]
	Bstop     = Tracker["constants"]["B_stop"]
	aa        = Tracker["constants"]["aa"]
	postlowpassfilter = Tracker["constants"]["postlowpassfilter"]
	B_enhance = Tracker["constants"]["B_enhance"]
	Tracker["number_of_groups"] = number_of_groups
	if Tracker["nosmearing"]:
		if(Blockdata["myid"] == Blockdata["main_node"]):
			msg =" Reconstruct maps of clusters without smearing"
			log_main.add(msg)
			print(msg)
			map_dir = os.path.join(masterdir, "maps_dir")
			cmd="{} {}".format("mkdir",map_dir)
			cmdexecute(cmd)
		else:map_dir = 0
		map_dir = wrap_mpi_bcast(map_dir, Blockdata["main_node"], MPI_COMM_WORLD)
		Tracker["directory"] = map_dir
		Tracker["nxinit"] = Tracker["constants"]["nnxo"]
		compute_noise(Tracker["nxinit"])
		data = get_shrink_data_sorting(os.path.join(Tracker["constants"]["masterdir"], "final_partition.txt"), \
		os.path.join(Tracker["constants"]["masterdir"],"refinement_parameters.txt"), \
		return_real = False, preshift = True, apply_mask = False)
		do3d_sorting_groups_trl_iter(data, 0)
		del data
		if(Blockdata["myid"] == Blockdata["main_node"]):
			for icluster in xrange(number_of_groups):
				cmd = "{} {} {}".format("cp ", \
				os.path.join(Tracker["directory"], "vol_unfiltered_0_grp%03d_iter000.hdf"%icluster), \
				os.path.join(Tracker["constants"]["masterdir"], "vol_unfiltered_0_grp%03d.hdf"%icluster))
				cmdexecute(cmd)
				cmd = "{} {} {}".format("cp ", \
				os.path.join(Tracker["directory"], "vol_unfiltered_1_grp%03d_iter000.hdf"%icluster), \
				os.path.join(Tracker["constants"]["masterdir"], "vol_unfiltered_1_grp%03d.hdf"%icluster))
				cmdexecute(cmd)
		mpi_barrier(MPI_COMM_WORLD)
	else:
		if(Blockdata["myid"] == Blockdata["main_node"]):
			msg =" Reconstruct maps of clusters with smearing"
			log_main.add(msg)
			print(msg)
		for icluster in xrange(Tracker["number_of_groups"]):
			cluster_masterdir = os.path.join(masterdir,"Cluster%d"%icluster)
			if(Blockdata["myid"] == Blockdata["main_node"]):
				cmd = "{} {}".format("mkdir", os.path.join(cluster_masterdir))
				if os.path.exists(cluster_masterdir): print("directory already exists! ")
				else: cmdexecute(cmd)
			mpi_barrier(MPI_COMM_WORLD)
			do_ctrefromsort3d_get_subset_data(cluster_masterdir, refinement_dir, \
			  os.path.join(masterdir,"Cluster%d.txt"%icluster), selected_iter, None, Blockdata["subgroup_comm"])
			Tracker["constants"]["small_memory"] = False
			ctrefromsorting_rec3d_faked_iter(cluster_masterdir, selected_iter, rec3d_image_size, Blockdata["subgroup_comm"])
			mpi_barrier(MPI_COMM_WORLD)
		mpi_barrier(MPI_COMM_WORLD)
		Tracker["constants"]["B_enhance"] = B_enhance
		Tracker["constants"]["B_start"] = Bstart    
		Tracker["constants"]["B_stop"]  = Bstop    
		Tracker["constants"]["aa"]      = aa  
		Tracker["constants"]["postlowpassfilter"] = postlowpassfilter  
		Tracker["constants"]["fsc_adj"]=fsc_adj
		Tracker["constants"]["mtf"]    = mtf
		Tracker["mask3D"]              = mask3D
		Tracker["nxinit"]              = rec3d_image_size 
		Tracker["number_of_groups"]    = number_of_groups
		Tracker["fuse_freq"]           = fuse_freq # reset
		# Using all CPUS to do step two
		Blockdata["ncpuspernode"] = Blockdata["nproc"]//Blockdata["no_of_groups"]
		Blockdata["nsubset"]  = Blockdata["ncpuspernode"]*Blockdata["no_of_groups"]
		create_subgroup()
		do3d_sorting_groups_rec3d(selected_iter, masterdir, log_main)
		if(Blockdata["myid"] == Blockdata["main_node"]):
			for icluster in xrange(Tracker["number_of_groups"]):
				cluster_masterdir = os.path.join(masterdir,"Cluster%d"%icluster)
				cmd = "{} {}".format("rm -rf ", os.path.join(cluster_masterdir))
				if os.path.exists(cluster_masterdir): cmdexecute(cmd)
	return
#####<<<<--------------------------------------
def merge_two_unfiltered_maps(map1_file, map2_file, cluster_ID):
	global Tracker, Blockdata
	# single processor only
	from math import sqrt, log
	try: map1 = get_im(map1_file)
	except: ERROR("Sphire postprocess fails to read the first map "+map1_file, "--postprocess option for 3-D", 1, Blockdata["myid"])
	try: map2 = get_im(map2_file)
	except: ERROR("Sphire postprocess fails to read the second map "+map2_file, "--postprocess option for 3-D", 1, Blockdata["myid"])
	if (map2.get_xsize() != map1.get_xsize()) or (map2.get_ysize() != map1.get_ysize()) or (map2.get_zsize() != map1.get_zsize()):
		ERROR(" Two input maps have different image size", "--postprocess option for 3-D", 1, Blockdata["myid"])
	if Tracker["mask3D"]: 
		mask3D = get_im(Tracker["mask3D"])
		if mask3D.get_xsize() != Tracker["nxinit"]: 
			mask3D = fdecimate(mask3D, Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"], True, False)
	else: mask3D = None
	## prepare FSC
	resolution_FSC143   = 0.5 # for single volume, this is the default resolution
	resolution_FSChalf  = 0.5
	if mask3D: fsc_true = fsc(map1*mask3D, map2*mask3D, 1)
	else:      fsc_true = fsc(map1, map2, 1)
	resolution_in_angstrom = [None]*len(fsc_true[0])
	for ifreq in xrange(len(fsc_true[0])):
		if fsc_true[0][ifreq] !=0.0: resolution_in_angstrom [ifreq] = Tracker["constants"]["pixel_size"]/fsc_true[0][ifreq]
		else: resolution_in_angstrom [ifreq] = 0.0
	fsc_true[1][0] =1.0  # always reset fsc of zero frequency as 1.0
	for ifreq in xrange(len(fsc_true[0])): 
		fsc_true[1][ifreq] = max(fsc_true[1][ifreq], 0.0)*2./(1.+max(fsc_true[1][ifreq], 0.0))
	resolution_FSC143_right  = 0.0
	resolution_FSC143_left   = 0.0
	dip_at_fsc = False
	nfreq0 = 1
	for ifreq in xrange(1, len(fsc_true[1])):
		if fsc_true[1][ifreq] < 0.0:
			nfreq0  = ifreq - 1
			break
	if nfreq0 ==1: nfreq0= len(fsc_true[1]) - 1

	nfreq05 = len(fsc_true[1])-1 		
	for ifreq in xrange(1, len(fsc_true[1])):
		if fsc_true[1][ifreq] < 0.5:
			resolution_FSChalf = fsc_true[0][ifreq-1]
			nfreq05 = ifreq-1
			break

	resolution_FSC143_left = fsc_true[0][len(fsc_true[1])-1]
	for ifreq in xrange(nfreq05, len(fsc_true[1])):
		if fsc_true[1][ifreq] < 0.143:
			resolution_FSC143_left = fsc_true[0][ifreq-1]
			nfreq143 = ifreq - 1
			break
	
	#print(nfreq0, nfreq05, "check")
	nfreq143_right = nfreq0
	resolution_FSC143_right = fsc_true[0][nfreq05]
	for ifreq in xrange(nfreq0, nfreq05, -1):
		if fsc_true[1][ifreq] >= 0.143:
			resolution_FSC143_right = fsc_true[0][ifreq]
			nfreq143_right = ifreq
			break
	resolution_FSC143 = resolution_FSC143_right
	nfreq143 = nfreq143_right
	## smooth FSC after FSC143 and set other values to zero
	for ifreq in xrange(nfreq143+1, len(fsc_true[1])):
		if ifreq ==nfreq143+1: fsc_true[1][ifreq] = (fsc_true[1][nfreq143-2] + fsc_true[1][nfreq143-1])/5.
		elif ifreq ==nfreq143+2: fsc_true[1][ifreq] = (fsc_true[1][nfreq143-1])/5.
		else: fsc_true[1][ifreq] = 0.0
	fsc_out = []
	for ifreq in xrange(len(fsc_true[0])): fsc_out.append("%5d   %7.2f   %5.3f"%(ifreq, resolution_in_angstrom[ifreq],fsc_true[1][ifreq]))
	write_text_file(fsc_out, "fsc_%d.txt"%cluster_ID)														
	map1 +=map2 #(get_im(args[0])+get_im(args[1]))/2.0
	map1 /=2.0
	del map2
	outtext = [["Squaredfreq"],[ "LogOrig"]]
	guinierline = rot_avg_table(power(periodogram(map1),.5))
	for ig in xrange(len(guinierline)):
		x = ig*.5/float(len(guinierline))/Tracker["constants"]["pixel_size"]
		outtext[0].append("%10.6f"%(x*x))
		outtext[1].append("%10.6f"%log(guinierline[ig]))
	if Tracker["constants"]["mtf"]: # divided by the mtf #1
		log_main.add("MTF correction is applied")
		log_main.add("MTF file is %s"%Tracker["constants"]["mtf"])
		try: mtf_core  = read_text_file(Tracker["constants"]["mtf"], -1)
		except: ERROR("Sphire postprocess fails to read MTF file "+Tracker["constants"]["mtf"], "--postprocess option for 3-D", 1)
		map1 = fft(Util.divide_mtf(fft(map1), mtf_core[1], mtf_core[0]))
		outtext.append(["LogMTFdiv"])
		guinierline   = rot_avg_table(power(periodogram(map1),.5))
		for ig in xrange(len(guinierline)): outtext[-1].append("%10.6f"%log(guinierline[ig]))
	if Tracker["constants"]["fsc_adj"]: #2
		#### FSC adjustment ((2.*fsc)/(1+fsc)) to the powerspectrum;
		fil = len(fsc_true[1])*[None]
		for i in xrange(len(fil)): fil[i] = sqrt(fsc_true[1][i]) # fsc already matched to full dataset
		map1 = filt_table(map1,fil)
		guinierline = rot_avg_table(power(periodogram(map1),.5))
		outtext.append(["LogFSCadj"])
		for ig in xrange(len(guinierline)):outtext[-1].append("%10.6f"%log(guinierline[ig]))
	if Tracker["constants"]["B_enhance"] !=-1: #3 One specifies and then apply B-factor sharpen
		if Tracker["constants"]["B_enhance"] == 0.0: # auto mode
			cutoff_by_fsc = 0
			for ifreq in xrange(len(fsc_true[1])):
				if fsc_true[1][ifreq]<0.143: break
			cutoff_by_fsc = float(ifreq-1)
			freq_max = cutoff_by_fsc/(2.*len(fsc_true[0]))/Tracker["constants"]["pixel_size"]
			guinierline = rot_avg_table(power(periodogram(map1),.5))
			logguinierline = []
			for ig in xrange(len(guinierline)):logguinierline.append(log(guinierline[ig]))
			freq_min = 1./Tracker["constants"]["B_start"] # given frequencies in Angstrom unit, say, B_start is 10 Angstrom, or 15  Angstrom
			if Tracker["constants"]["B_stop"]!=0.0: 
				freq_max = 1./Tracker["constants"]["B_stop"]
				B_stop = Tracker["constants"]["B_stop"]
			else: B_stop = Tracker["constants"]["pixel_size"]/(float(nfreq143)/Tracker["constants"]["nnxo"])
			if freq_min>= freq_max:
				print("Your B_start is too high!")
				freq_min = 1./(B_stop + 8.)
			b, junk, ifreqmin, ifreqmax = compute_bfactor(guinierline, freq_min, freq_max, Tracker["constants"]["pixel_size"])
			global_b = min(4.*b, 400.) # B-factor should not be too large
			cc = pearson(junk[1],logguinierline)
			sigma_of_inverse = sqrt(2./(global_b/Tracker["constants"]["pixel_size"]**2))		
		else: # User provided value
			sigma_of_inverse = sqrt(2./((abs(Tracker["constants"]["B_enhance"]))/Tracker["constants"]["pixel_size"]**2))
			global_b = Tracker["constants"]["B_enhance"]
		print("B_factor", global_b, "Cluster ", cluster_ID)
		map1 = filt_gaussinv(map1, sigma_of_inverse)
		guinierline = rot_avg_table(power(periodogram(map1),.5))
		last_non_zero = -999.0
		for ig in xrange(len(guinierline)):
			if guinierline[ig]>0: 
				outtext[-1].append("%10.6f"%log(guinierline[ig]))
				last_non_zero = log(guinierline[ig])
			else: outtext[-1].append("%10.6f"%last_non_zero)
	if Tracker["constants"]["postlowpassfilter"] != 0.0: # User provided low-pass filter #4.
		if Tracker["constants"]["postlowpassfilter"] > 0.5: # Input is in Angstrom 
			map1 = filt_tanl(map1,Tracker["constants"]["pixel_size"]/Tracker["constants"]["postlowpassfilter"], min(Tracker["constants"]["aa"],.1))
		elif Tracker["constants"]["postlowpassfilter"]>0.0 and Tracker["constants"]["postlowpassfilter"]<0.5:  # input is in absolution frequency
			map1 = filt_tanl(map1,Tracker["constants"]["postlowpassfilter"], min(Tracker["constants"]["aa"],.1))
	else: 
		map1 = filt_tanl(map1,resolution_FSC143, Tracker["constants"]["aa"])
		print("low_pass_filter_to", resolution_FSC143, "Cluster ", cluster_ID)
	map1.write_image(os.path.join(Tracker["constants"]["masterdir"], "vol_final_nomask_cluster%d.hdf"%cluster_ID))
	if mask3D: map1 *=mask3D
	map1.write_image(os.path.join(Tracker["constants"]["masterdir"], "vol_final_cluster%d.hdf"%cluster_ID))
	if mask3D: del mask3D
	del map1
	return
	
def main():
	from optparse   import OptionParser
	from global_def import SPARXVERSION
	from EMAN2      import EMData
	from logger     import Logger, BaseLogger_Files
	from global_def import ERROR  
	import sys, os, time
	global Tracker, Blockdata 
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack  outdir --refinement_dir=masterdir_of_sxmeridien --mask3D=mask.hdf --focus=binarymask.hdf  --radius=outer_radius " +\
	"  --sym=c1  --nindependent=indenpendent_runs  --img_per_grp=img_per_grp  "
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--refinement_dir",                  type   ="string",        default ='',                       help="3-D refinement directory, the master directory of sxmeridien")
	parser.add_option("--output_dir",                      type   ="string",        default ='',					   help="name of the directory for sorting computing")
	parser.add_option("--niter_for_sorting",               type   ="int",           default =-1,					   help="selected number of iteration of meridien refinement for sorting, -1 implies program uses the best iteration to initiate sorting")
	parser.add_option("--focus",                           type   ="string",        default ='',                       help="file name, the bineary 3D mask for focused clustering ")
	parser.add_option("--mask3D",                          type   ="string",        default ='',                       help="file name, the 3-D global mask for clustering ")
	parser.add_option("--instack",                         type   ="string",        default ='',					   help="file name, data stack for sorting provided by user. It applies when sorting starts from a given data stack")
	parser.add_option("--radius",                          type   ="int",           default =-1,	                   help="particle radius in pixel for rotational correlation <nx-1 (set to the radius of the particle)")
	parser.add_option("--sym",                             type   ="string",        default ='c1',                     help="point group symmetry of macromolecular structure, can be inherited from refinement")
	parser.add_option("--nindependent",                    type   ="int",           default = 3,                       help="number of independent run for EQkmeans clustering, an odd number larger than 2")
	parser.add_option("--img_per_grp",                     type   ="int",           default =1000,                     help="number of images in a group")
	parser.add_option("--minimum_grp_size",				   type   ="int",           default =500,					   help="minimum number of members for being identified as a group")
	parser.add_option("--comparison_method",               type   ="string",        default ='cross',                  help="option for comparing two images, either using cross-correlaton coefficients [cross] or using Euclidean distance [eucd] ")
	parser.add_option("--memory_per_node",                 type   ="float",         default =-1.0,                     help="memory_per_node, the number used for evaluate the CPUs/NODE settings given by user")
	parser.add_option("--nofinal_sharpen",                 action ="store_true",    default =False,                    help="not reconstruct unfiltered final maps for post refinement process")
	parser.add_option("--eqkmeans_angle_step",             type   ="float",         default =15.,                      help="smapling anglular step used for EQKmeans orientation constraints")
	parser.add_option("--eqkmeans_tilt1",                  type   ="float",         default =0.,                       help="sampling starting theta angle used for EQKmeans orientation constraints")
	parser.add_option("--eqkmeans_tilt2",                  type   ="float",         default =180.,                     help="sampling end theta angle used for EQKmeans orientation constraints")
	parser.add_option("--post_sorting_sharpen",            action ="store_true",    default =False,                    help="make sharpen maps from sorted clusters ")
	parser.add_option("--stop_eqkmeans_percentage",        type   ="float",         default =2.0,                      help="particle change percentage for stopping equal size Kmeans")
	parser.add_option("--minimum_ptl_number",              type   ="int",           default =20,					   help="integer number, the smallest orien group size equals number_of_groups multiplies this number")
	parser.add_option("--notapplybckgnoise",               action ="store_true",    default =False,                    help="do not applynoise")
	#parser.add_option("--nofftwmpi",                       action ="store_true",    default =False,                    help="Use fftwmpi (default False)")
	# postprocessing options
	parser.add_option("--mtf",                             type   ="string",        default ='',                       help="mtf file")
	parser.add_option("--B_enhance",                       type   ="float",         default=0.0,                       help="apply Bfactor to enhance map or not")
	parser.add_option("--fl",                              type   ="float",         default=0.0,                       help="=0.0, low_pass filter to resolution limit; =some value, low_pass filter to some valume; =-1, not low_pass filter applied")
	parser.add_option("--aa",                              type   ="float",         default=.1,                        help="low pass filter falloff")
	parser.add_option("--B_start",                         type   ="float",         default=10.,                       help="starting frequency in Angstrom for B-factor estimation")
	parser.add_option("--B_stop",                          type   ="float",         default=0.0,                       help="cutoff frequency in Angstrom for B-factor estimation, cutoff is set to the frequency where fsc < 0.0")
	parser.add_option("--nofsc_adj",                       action ="store_true",    default =False,                    help="do not multiply sqrt(fsc)")
	(options, args) = parser.parse_args(sys.argv[1:])
	from utilities import bcast_number_to_all
	### Sanity check
	
	if options.refinement_dir !='':
		if not os.path.exists(options.refinement_dir): ERROR("The specified refinement_dir does not exist", "sort3d", 1, Blockdata["myid"])
	if options.focus !='':
		if not os.path.exists(options.focus): ERROR("The specified focus mask file does not exist", "sort3d", 1, Blockdata["myid"])
	if options.mask3D !='':
		if not os.path.exists(options.mask3D): ERROR("The specified mask3D file does not exist", "sort3d", 1, Blockdata["myid"])
	if options.img_per_grp <= options.minimum_grp_size:ERROR("img_per_grp should be way larger than minimum_grp_size", "sort3d", 1, Blockdata["myid"])
	
	#--- Fill input parameters into dictionary Constants
	Constants		                         = {}
	Constants["stop_eqkmeans_percentage"]    = options.stop_eqkmeans_percentage
	if options.nofinal_sharpen: Constants["final_sharpen"] = False
	else:  Constants["final_sharpen"] = True
	Constants["niter_for_sorting"]           = options.niter_for_sorting
	Constants["memory_per_node"]             = options.memory_per_node
	Constants["orgstack"]                    = options.instack
	Constants["masterdir"]                   = options.output_dir
	Constants["refinement_dir"]              = options.refinement_dir
	Constants["mask3D"]                      = options.mask3D
	Constants["focus3Dmask"]                 = options.focus	
	Constants["indep_runs"]                  = options.nindependent
	Constants["indep_runs"]                  = max(3, Constants["indep_runs"])
	if Constants["indep_runs"]%2 ==0: Constants["indep_runs"]+=1 # always an odd number 	
	Constants["img_per_grp"]  = options.img_per_grp
	Constants["minimum_grp_size"]      		 = options.minimum_grp_size
	Constants["CTF"]                		 = True
	Constants["radius"]              		 = options.radius
	Constants["symmetry"]                    = options.sym
	Constants["nxinit"]                      = -1 # retrive it under request, equals to the role of nxinit of meridien
	Constants["seed"]                        = -1
	Constants["upscale"]                     = 0.5 #
	Constants["interpolation"]               = "trl"
	Constants["comparison_method"]           = options.comparison_method # either cross or eucd
	### postprocessing
	if options.mtf =='':   Constants["mtf"]  = None
	else:                  Constants["mtf"]  = options.mtf
	Constants["B_enhance"]                   = options.B_enhance
	Constants["B_start"]                     = options.B_start
	Constants["B_stop"]                      = options.B_stop
	Constants["postlowpassfilter"]           = options.fl
	Constants["aa"]                          = options.aa
	if options.nofsc_adj: Constants["fsc_adj"] = False
	else:   Constants["fsc_adj"] = True	
	if options.focus:  Constants["comparison_method"] = "cross" # in case of focus3D, cross is used.
	Constants["fuse_freq"] = 45.  # Now in A, convert to pixels before being used
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
	if options.notapplybckgnoise: Tracker["applybckgnoise"] = False
	else: Tracker["applybckgnoise"] = True 	
	###<<<--options for advanced users:
	Tracker["total_sort3d_indepent_run"]	       = 2
	Tracker["total_number_of_iterations"] 	       = 25
	Tracker["total_iter_rsort"]                    = 2
	Tracker["clean_volumes"]                       = True
	Tracker["constants"]["total_sort3d_iteration"] = 2
	### -----------orientation constraints
	Tracker["angle_step"]            = options.eqkmeans_angle_step # orientation constrained angle step
	Tracker["tilt1"]                 = options.eqkmeans_tilt1
	Tracker["tilt2"]                 = options.eqkmeans_tilt2
	### ------------<<< option for proteins images that have preferred orientations
	Tracker["minimum_ptl_number"] = options.minimum_ptl_number  # for orientation groups	
	if Tracker["constants"]["memory_per_node"] ==-1 or Tracker["constants"]["memory_per_node"] <32.: Tracker["constants"]["small_memory"] = True
	else: Tracker["constants"]["small_memory"] = False
	## Additional check
	Tracker["constants"]["hardmask"] =True
	Tracker["applymask"]             =True
	if os.path.exists(options.refinement_dir):
		Tracker["constants"]["refinement_method"] ="SPARX"
		Tracker["nosmearing"]                     = False
	elif options.instack !='':
		Tracker["constants"]["refinement_method"] ="stack"
		Tracker["nosmearing"]                     = True
		Tracker["constants"]["refinement_dir"]    = None
		Tracker["paramstructure_dir"]             = None
		Tracker["refang"]                         = None
		Tracker["rshifts"]                        = None
		Tracker["paramstructure_dict"]            = None
		Tracker["constants"]["selected_iter"]     = -1
	else:  ERROR("Incorrect refinement_dir ", "sxsort3d_new.py", 1, Blockdata["myid"])
	if os.path.exists(options.refinement_dir) and options.instack !='': ERROR("contractdict refinement methods", "sxsort3d_smearing.py", 1, Blockdata["myid"])	
	Blockdata["fftwmpi"] = True
	
	Blockdata["symclass"] = symclass(Tracker["constants"]["symmetry"])
	
	#else:  Blockdata["fftwmpi"] = False
	Blockdata["ncpuspernode"] = Blockdata["no_of_processes_per_group"]
	Blockdata["nsubset"] = Blockdata["ncpuspernode"]*Blockdata["no_of_groups"]
	create_subgroup()
	###--------------------------------------------------------------------------------------------
	#    Two typical sorting scenarios  
	#
	# 1. import data and refinement parameters from meridien refinement;
	# 2. given data stack and xform.projection/ctf in header(For simulated test data);
	#<<<---------------------->>>imported functions<<<---------------------------------------------
	from statistics 	import k_means_match_clusters_asg_new,k_means_stab_bbenum
	from utilities 		import get_im,bcast_number_to_all,cmdexecute,write_text_file,read_text_file,wrap_mpi_bcast, get_params_proj, write_text_row
	from utilities 		import get_number_of_groups
	from filter			import filt_tanl
	from time           import sleep
	from logger         import Logger,BaseLogger_Files
	import string
	from string         import split, atoi, atof
	import json
	import user_functions
	####--------------------------------------------------------------
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
			cmdexecute(cmd)
		else:
			if os.path.exists(masterdir): restart = 1
			else:
				cmd="{} {}".format("mkdir", masterdir)
				cmdexecute(cmd)
		li =len(masterdir)
	else:li = 0
	restart						          = bcast_number_to_all(restart, Blockdata["main_node"])
	li                                    = mpi_bcast(li,1,MPI_INT,Blockdata["main_node"],MPI_COMM_WORLD)[0]
	masterdir                             = mpi_bcast(masterdir,li,MPI_CHAR,Blockdata["main_node"],MPI_COMM_WORLD)
	masterdir                             = string.join(masterdir,"")
	Tracker["constants"]["masterdir"]     = masterdir
	Tracker["constants"]["chunk_0"]       = os.path.join(Tracker["constants"]["masterdir"],"chunk_0.txt")
	Tracker["constants"]["chunk_1"]       = os.path.join(Tracker["constants"]["masterdir"],"chunk_1.txt")              
	log_main = Logger(BaseLogger_Files())
	log_main.prefix = Tracker["constants"]["masterdir"]+"/"
	import_from_relion_refinement = 0
	import_from_sparx_refinement  = 0
	import_from_data_stack		  = 0
	total_stack					  = 0
	while not os.path.exists(Tracker["constants"]["masterdir"]):
		print("Node ", Blockdata["myid"], "waiting...", Tracker["constants"]["masterdir"])
		sleep(1)
	mpi_barrier(MPI_COMM_WORLD)
	if(Blockdata["myid"] == Blockdata["main_node"]):
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		msg = "importing data ... "
		print(line, msg)
		log_main.add(msg)
	######### Extract refinement information
	if Tracker["constants"]["refinement_method"] =="SPARX": # Senario one
		import_from_sparx_refinement = get_input_from_sparx_ref3d(log_main)
		Tracker["smearing"] = True
	else:  # Senario three, sorting from a given data stack, general cases
		#if not Tracker["constants"]["nofinal_sharpen"]: Tracker["constants"]["nofinal_sharpen"] = True
		import_from_data_stack = get_input_from_datastack(log_main)
		Tracker["constants"]["hardmask"] = True
		Tracker["applybckgnoise"]        = False
		Tracker["applymask"]             = True
		Tracker["smearing"]              = False
	###<<<------------------------>>>>>>checks<<<<<-------------
	if Tracker["constants"]["symmetry"] != options.sym:
		msg = "input symmetry %s is altered to %s after reading refinement information! "%(options.sym, Tracker["constants"]["symmetry"])
		log_main.add(msg)
		print(msg)
	###<<<----------------------->>>>>  SORT3D MAIN PROGRAM <<<<<---------------------------------------------# For all cases
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if(Blockdata["myid"] == Blockdata["main_node"]):
		print(line, "Sort3d main program")
		log_main.add("---------->>>SPARX sort3d<<<--------------")
		log_main.add("The shell line command:")
		line = ""
		for a in sys.argv: line +=(a + " ")
		log_main.add(line)
		log_main.add("Sort3d master directory: %s"%Tracker["constants"]["masterdir"])
		print_dict(Tracker["constants"],"Permanent sorting settings after initialization")
	mpi_barrier(MPI_COMM_WORLD)
	####--->>>>>> Check 3-D mask<<<<----------
	if Tracker["constants"]["mask3D"]:Tracker["mask3D"] = os.path.join(Tracker["constants"]["masterdir"],"smask.hdf")
	else: Tracker["mask3D"] = None
	if Tracker["constants"]["focus3Dmask"]:Tracker["focus3D"] = Tracker["constants"]["focus3Dmask"]
	else: Tracker["focus3D"] = None
	if(Blockdata["myid"] == Blockdata["main_node"]):
		bad_focus3Dmask = 0
		if Tracker["constants"]["focus3Dmask"]:
			try:
				focusmask = get_im(Tracker["constants"]["focus3Dmask"])
				st = Util.infomask(binarize(focusmask), None, True)
				if(st[0] == 0.0): bad_focus3Dmask = 1
				else:
					focusmask.write_image(os.path.join(Tracker["constants"]["masterdir"], "focus3d.hdf"))
					Tracker["focus3D"] = os.path.join(Tracker["constants"]["masterdir"], "focus3d.hdf")
			except:  bad_focus3Dmask = 1
	else: bad_focus3Dmask = 0
	bad_focus3Dmask = bcast_number_to_all(bad_focus3Dmask,	source_node =  Blockdata["main_node"])
	if bad_focus3Dmask: ERROR("Incorrect focused mask, after binarize all values zero","sxsort3d.py", 1, Blockdata["myid"])	
	if(Blockdata["myid"] == Blockdata["main_node"]):
		bad_3Dmask = 0
		if Tracker["constants"]["mask3D"]:
			try: 
				mask3D = get_im(Tracker["constants"]["mask3D"])
				st = Util.infomask(binarize(mask3D), None, True)
				if (st[0] ==0.0): bad_3Dmask = 1
				else: 
					mask3D.write_image(os.path.join(Tracker["constants"]["masterdir"], "mask3D.hdf"))
					Tracker["mask3D"]= os.path.join(Tracker["constants"]["masterdir"], "mask3D.hdf")
			except: bad_3Dmask = 1
	else: bad_3Dmask = 0
	bad_3Dmask = bcast_number_to_all(bad_focus3Dmask,	source_node =  Blockdata["main_node"])
	if bad_3Dmask: ERROR("Incorrect 3D mask", "sxsort3d.py", 1, Blockdata["myid"])
	Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], MPI_COMM_WORLD)
	############################################################################################	
	###<<<------ Determine the image size ### reset nxinit simply for calculate currentres!
	Tracker["nxinit"] = Tracker["nxinit_refinement"]
	Tracker["currentres"] = float(Tracker["constants"]["fsc05"])/float(Tracker["nxinit"])
	##################---------------<<<shrinkage, current resolution, fuse_freq  <<<<<<------------------------------------------
	Tracker["total_stack"] = Tracker["constants"]["total_stack"]
	Tracker["shrinkage"] = float(Tracker["nxinit"])/Tracker["constants"]["nnxo"]
	Tracker["radius"]  = Tracker["constants"]["radius"]*Tracker["shrinkage"]
	try: fuse_freq = Tracker["fuse_freq"]
	except: Tracker["fuse_freq"] = int(Tracker["constants"]["pixel_size"]*Tracker["constants"]["nnxo"]/Tracker["constants"]["fuse_freq"]+0.5)
	###################------------------->>>>> 3-D masks <<<<<-----------------------------------------------------------
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if(Blockdata["myid"] == Blockdata["main_node"]):
		print_dict(Tracker["constants"],"Permanent sorting settings from input options")
		fout = open(os.path.join(Tracker["constants"]["masterdir"], "Tracker.json"),'w')
		json.dump(Tracker, fout)
		fout.close()
		msg = "---------->>>sort3d<<<-----------"
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
			msg ="User provided focus mask file:  %s"%Tracker ["constants"]["focus3Dmask"]
			print(line, msg)
			log_main.add(msg)
		Tracker["full_list"] = read_text_file(os.path.join(Tracker["constants"]["masterdir"],"indexes.txt"), -1) # could have one or two columns
		if len(Tracker["full_list"]) == 2: Tracker["full_list"] = Tracker["full_list"][1] # take just one column
		elif len(Tracker["full_list"]) == 1: Tracker["full_list"] = Tracker["full_list"][0]
		else: ERROR("The original particle ID for sorting has wrong format", "sxsort3d.py", 1, Blockdata["main_node"])
	else: Tracker["full_list"] = 0
	Tracker["full_list"] = wrap_mpi_bcast(Tracker["full_list"], Blockdata["main_node"], MPI_COMM_WORLD)
	Tracker["shrinkage"] = float(Tracker["nxinit"])/Tracker["constants"]["nnxo"]
	if(Blockdata["myid"] == Blockdata["main_node"]): print_dict(Tracker,"Current sorting settings")
### Check mpi setting
	check_mpi_settings(log_main)
	if options.post_sorting_sharpen: # post sorting option!
		try: nxinit = Tracker["nxinit"] 
		except: Tracker["nxinit"]                 = -1
		Tracker["constants"]["orgres"]            = 0.0
		Tracker["constants"]["refinement_delta"]  = 0.0
		Tracker["constants"]["refinement_ts"]     = 0.0
		Tracker["constants"]["refinement_xr"]     = 0.0
		Tracker["constants"]["refinement_an"]     = 0.0
		number_of_groups = 0
		minimum_size = Tracker["constants"]["img_per_grp"]
		if(Blockdata["myid"] == Blockdata["main_node"]):
			print(" reconstruct maps for sharpening ")
			while os.path.exists(os.path.join(Tracker["constants"]["masterdir"], "Cluster%d.txt"%number_of_groups)):
				class_in = read_text_file(os.path.join(Tracker["constants"]["masterdir"], "Cluster%d.txt"%number_of_groups))
				minimum_size = min(len(class_in), minimum_size)
				number_of_groups +=1
			print(" % clusters are found  "%number_of_groups)
		number_of_groups = bcast_number_to_all(number_of_groups, Blockdata["main_node"], MPI_COMM_WORLD)
		if number_of_groups == 0:ERROR("No cluster is found, and the program terminates. ", "option post_sorting_sharpen ", 1, Blockdata["myid"])
		minimum_size = bcast_number_to_all(minimum_size, Blockdata["main_node"], MPI_COMM_WORLD)
		compute_noise(Tracker["constants"]["nnxo"])
		do_final_maps(number_of_groups, minimum_size, Tracker["constants"]["selected_iter"], Tracker["constants"]["refinement_dir"], \
		   Tracker["constants"]["masterdir"], Tracker["constants"]["nnxo"], log_main)
		if(Blockdata["myid"] == Blockdata["main_node"]):
			for iproc in xrange(number_of_groups):
				merge_two_unfiltered_maps(os.path.join(Tracker["constants"]["masterdir"], "vol_unfiltered_0_grp%03d.hdf"%iproc), \
				os.path.join(Tracker["constants"]["masterdir"], "vol_unfiltered_1_grp%03d.hdf"%iproc), iproc)
		mpi_barrier(MPI_COMM_WORLD)
		from mpi import mpi_finalize
		mpi_finalize()
		exit()
	###<<<-------sort3d starts here
	for indep_sort3d in xrange(Tracker["total_sort3d_indepent_run"]):		
		sorting                = {}
		sorting["Unaccounted"] = None
		sorting["Accounted"]   = []
		partid_file  = os.path.join(Tracker["constants"]["masterdir"],"indexes.txt")
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		Tracker["indep_sort3d_dir"] = os.path.join(Tracker["constants"]["masterdir"], "sort3d_run%d"%indep_sort3d)
		if Blockdata["myid"] == Blockdata["main_node"]:
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			if os.path.exists(Tracker["indep_sort3d_dir"]):
				msg = "directory %s already exists "%Tracker["indep_sort3d_dir"]	
				log_main.add(msg)
			else:  cmdexecute("{} {}".format("mkdir",Tracker["indep_sort3d_dir"]))
			msg = "---------->>> Independent sort3d  %d<<<----------- "%indep_sort3d
			print(line, msg)
			log_main.add(msg)
			log_main.add("nnxo : %d"%Tracker["constants"]["nnxo"])
			log_main.add("Current resolution: %f  absolute unit(maximum is 0.5) 1./%7.2f Angstrom "%(Tracker["currentres"]*Tracker["shrinkage"],\
                Tracker["constants"]["pixel_size"]/Tracker["currentres"]/Tracker["shrinkage"]))
			if Tracker["mask3D"]:
				msg = "User provided 3-D mask:  %s"%Tracker["constants"]["mask3D"]
				log_main.add(msg)
				print(line, msg)	
			if import_from_sparx_refinement ==1:
				msg = "Sorting is initiated from meridien refinement"
				print(line, msg)
				log_main.add(msg)
			elif import_from_data_stack ==1:
				msg = "Sorting is initiated from data stack"
				print(line, msg)
				log_main.add(msg)
			else: ERROR("importing neither from a meridien refinement nor a given data stack", "sort3d", 1, Blockdata["myid"])
			sorting["total"] = read_text_file(partid_file, -1)
			if len(sorting["total"])>1: sorting["total"] = sorting["total"][1]
			else:   sorting["total"] = sorting["total"][0]
		else: sorting["total"] = 0
		sorting["total"] = wrap_mpi_bcast(sorting["total"], Blockdata["main_node"]) # total number of records in indexes.txt file
		if Blockdata["myid"] == Blockdata["main_node"]:
			Tracker["img_per_grp"] = Tracker["constants"]["img_per_grp"]
			Tracker["total_stack"] = len(sorting["total"]) # start from beginning
			msg = "img_per_grp:  %d"%Tracker["img_per_grp"]
			print(line, msg)
			log_main.add(msg)
			Tracker["number_of_groups"] = get_number_of_groups(Tracker["total_stack"],Tracker["img_per_grp"])
		Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], MPI_COMM_WORLD)
		final_list = do_EQKmeans_nways_clustering_stable_seeds(Tracker["indep_sort3d_dir"], partid_file, \
		os.path.join(Tracker["constants"]["masterdir"],"refinement_parameters.txt"), sorting["Accounted"], log_main)
		sorting["Unaccounted"] = Tracker["unaccounted_list"]
		if Blockdata["myid"] == Blockdata["main_node"]:
			# sorting.json contains [Accounted] lists and one Unaccounted list
			fout = open(os.path.join(Tracker["indep_sort3d_dir"], "sorting.json"),'w')
			json.dump(sorting, fout)
			fout.close()
			sort3d_clusters = split_partition_into_clusters(sorting["Accounted"])
			sort3d_account_list, sort3d_res_partition_list = merge_classes_into_partition_list(sort3d_clusters)	
			write_text_row(sort3d_res_partition_list, os.path.join(os.path.join(Tracker["indep_sort3d_dir"], "sort3d_partition.txt")))
			Tracker["number_of_groups"] = len(sort3d_clusters)
			for icluster in xrange(len(sort3d_clusters)): write_text_file(sort3d_clusters[icluster], os.path.join(Tracker["indep_sort3d_dir"],"Cluster%d.txt"%icluster))
		else:   Tracker["number_of_groups"] = 0	
		Tracker["number_of_groups"] = bcast_number_to_all(Tracker["number_of_groups"], Blockdata["main_node"], MPI_COMM_WORLD)		
		###>>>>>--------->>>>rec3D<<<---------
		Tracker["nxinit"] = Tracker["nxinit_refinement"] # report resolution and structure in refinement image size
		data, norm_per_particle = get_shrink_data_sorting_smearing(os.path.join(Tracker["indep_sort3d_dir"], "sort3d_partition.txt"), \
		os.path.join(Tracker["constants"]["masterdir"],"refinement_parameters.txt"), \
		return_real = False, preshift = True, apply_mask = False)
		if Tracker["nosmearing"]: parameterstructure  = None	
		else: parameterstructure = read_paramstructure_for_sorting(os.path.join(Tracker["indep_sort3d_dir"], "sort3d_partition.txt"),\
		   Tracker["paramstructure_dict"], Tracker["paramstructure_dir"])
		Tracker["directory"] = Tracker["indep_sort3d_dir"]
		do3d_sorting_groups_trl_smearing_iter(data, parameterstructure, norm_per_particle, 0, True)
		del data
		if not Tracker["nosmearing"]:
			del parameterstructure
			del norm_per_particle
		### Summary of results given by sort3d 
		if(Blockdata["myid"] == Blockdata["main_node"]):
			line    = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			cmd = "{} {}".format("rm -rf", os.path.join(Tracker["directory"], "tempdir"))
			if os.path.exists(os.path.join(Tracker["directory"], "tempdir")): cmdexecute(cmd)
			msg = "--->>>Summary of indepenent sort3d run %d<<<----"%indep_sort3d
			print(line, msg)
			log_main.add(msg)
			for icluster in xrange(Tracker["number_of_groups"]):
				line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
				cfsc          = read_text_file(os.path.join(Tracker["indep_sort3d_dir"],"fsc_driver_grp%03d_iter000.txt"%icluster), -1)
				res05, res143 = get_res(cfsc[0])    
				msg = "group %3d   number of images: %8d     FSC05    %f     FSC143   %f"%(icluster, len(sort3d_clusters[icluster]), \
				float(res05)/float(Tracker["nxinit"]), float(res143)/float(Tracker["nxinit"]))
				print(line, msg)
				log_main.add(msg)
		mpi_barrier(MPI_COMM_WORLD)
	#########################################################<<<---rsort
	iter_rsort = 0
	while iter_rsort< Tracker["total_iter_rsort"]:
		ptp = []
		for indep_sort3d in xrange(Tracker["total_sort3d_indepent_run"]):
			indep_sort3d_dir = os.path.join(Tracker["constants"]["masterdir"], "sort3d_run%d"%indep_sort3d)
			if os.path.exists(os.path.join(indep_sort3d_dir, "sorting.json")):				
				if(Blockdata["myid"] == Blockdata["main_node"]):
					fout = open(os.path.join(indep_sort3d_dir, "sorting.json"), "r")
					res_sort3d = convert_json_fromunicode(json.load(fout))
					fout.close()
					merged_classes = split_partition_into_clusters(res_sort3d["Accounted"])
					if len(res_sort3d["Unaccounted"]) > 0: merged_classes.append(res_sort3d["Unaccounted"])
				else:
					res_sort3d = 0
					merged_classes = 0
				res_sort3d     = wrap_mpi_bcast(res_sort3d,  Blockdata["main_node"], MPI_COMM_WORLD)
				merged_classes = wrap_mpi_bcast(merged_classes, Blockdata["main_node"], MPI_COMM_WORLD)
				sptp           = prep_ptp_single(merged_classes, res_sort3d["total"])
				ptp.append(sptp)
			else: ERROR("sorting results do not exist", "rsort", 1, Blockdata["myid"])
		mpi_barrier(MPI_COMM_WORLD)
		accounted_list, unaccounted_list, new_index = do_two_way_comparison_single(ptp[0], ptp[1], len(res_sort3d["total"]))
		Tracker["unaccounted_list"] = []
		for index in xrange(len(unaccounted_list)): Tracker["unaccounted_list"].append(res_sort3d["total"][unaccounted_list[index]]) # always check
		Tracker["accounted_list"] = []
		for index in xrange(len(accounted_list)): Tracker["accounted_list"].append([new_index[index][0], res_sort3d["total"][new_index[index][1]]]) # always check
		Tracker["directory"] = os.path.join(Tracker["constants"]["masterdir"], "rsort%d"%iter_rsort)	
		if Blockdata["myid"] == Blockdata["main_node"]:
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			msg = "------------->>>rsort%d<<<----------- "%iter_rsort
			log_main.add(msg)
			print(line, msg)
			msg = "Summary of two sort3d runs"
			log_main.add(msg)
			print(line, msg)
			msg = "Accounted:    %d    Unaccounted:    %d"%(len(Tracker["accounted_list"]), len(Tracker["unaccounted_list"]))
			log_main.add(msg)
			print(line, msg)
			cmd = "{} {}".format("mkdir",Tracker["directory"])
			cmdexecute(cmd)
			if len(Tracker["unaccounted_list"])>0: write_text_file(Tracker["unaccounted_list"], os.path.join(Tracker["directory"], "Unaccounted.txt"))
			write_text_row(Tracker["accounted_list"], os.path.join(Tracker["directory"], "Accounted.txt"))
		mpi_barrier(MPI_COMM_WORLD)
		partid_file                 = os.path.join(os.path.join(Tracker["directory"], "Unaccounted.txt"))
		Tracker["total_stack"]	    = len(Tracker["unaccounted_list"])
		Tracker["number_of_groups"] = get_number_of_groups(Tracker["total_stack"],Tracker["img_per_grp"])
		if Tracker["number_of_groups"]>1: #continue N-ways_clustering on unaccounted particles
			sorting                = {}
			sorting["Accounted"]   = [Tracker["accounted_list"]] # keep the accounted
			final_sort             = []
			final_list = do_EQKmeans_nways_clustering_stable_seeds(os.path.join(Tracker["constants"]["masterdir"], "rsort%d"%iter_rsort), partid_file,\
			os.path.join(Tracker["constants"]["masterdir"],"refinement_parameters.txt"), sorting["Accounted"], log_main)
			for a in sorting["Accounted"]: final_sort.append(a)
			if (Blockdata["myid"] == Blockdata["main_node"]):
				try: 
					Tracker["unaccounted_list"] = read_text_file(Tracker["Unaccounted_on_disk"])
					if len(Tracker["unaccounted_list"])> 0:
						for iparticle in xrange(len(Tracker["unaccounted_list"])):
							Tracker["unaccounted_list"][iparticle] = [0, Tracker["unaccounted_list"][iparticle]]
						final_sort.append(Tracker["unaccounted_list"])
					else:
						msg = "Empty unaccounted_list is found"
						print(line, msg)
						log_main.add(msg)
				except: print("no saved unaccounted text")								
				indexed_particle_list, Tracker["number_of_groups"] = merge_original_id_lists(final_sort)
				write_text_row(indexed_particle_list, os.path.join(Tracker["constants"]["masterdir"], "rsort%d"%iter_rsort, "index_for_Kmeans.txt"))
				line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
				msg  = "Include %d unaccounted as one cluster  "%len(Tracker["unaccounted_list"])
				print(line, msg)
				log_main.add(msg)
		else: 
			if (Blockdata["myid"] == Blockdata["main_node"]):
				final_sort  = []
				final_sort.append(Tracker["accounted_list"])
				# group-index unaccounted_list
				try: 
					Tracker["unaccounted_list"] = read_text_file(Tracker["Unaccounted_on_disk"])
					if len(Tracker["unaccounted_list"])> 0:
						for iparticle in xrange(len(Tracker["unaccounted_list"])):
							Tracker["unaccounted_list"][iparticle] = [0, Tracker["unaccounted_list"][iparticle]]
						final_sort.append(Tracker["unaccounted_list"])
					else:
						msg = "Empty unaccounted_list is found"
						print(line, msg)
						log_main.add(msg)
				except:
					print("no saved unaccounted text")
					Tracker["unaccounted_list"] = []
				indexed_particle_list, Tracker["number_of_groups"] = merge_original_id_lists(final_sort)
				write_text_row(indexed_particle_list, os.path.join(Tracker["constants"]["masterdir"], "rsort%d"%iter_rsort, "index_for_Kmeans.txt"))
				line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
				msg = "include %d unaccounted particles as one cluster to do two-way comparison "%len(Tracker["unaccounted_list"])
				print(line, msg)
				log_main.add(msg)
		Tracker["number_of_groups"] = bcast_number_to_all(Tracker["number_of_groups"], Blockdata["main_node"], MPI_COMM_WORLD)				
		Tracker["directory"] = os.path.join(Tracker["constants"]["masterdir"],"rsort%d"%iter_rsort,"Kmeans")
		iter_rsort  +=1	
	# rsort final comparison
	mpi_barrier(MPI_COMM_WORLD)
	if(Blockdata["myid"] == Blockdata["main_node"]):
		line    = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		tlist   = read_text_file(os.path.join(Tracker["constants"]["masterdir"],"indexes.txt"))
		ptp     = []
		for irsort in xrange(Tracker["total_iter_rsort"]):
			r1 = read_text_row(os.path.join(Tracker["constants"]["masterdir"], "rsort%d"%irsort, "index_for_Kmeans.txt"))
			merged_classes = split_partition_into_clusters([r1])
			sptp = prep_ptp_single(merged_classes, tlist)
			ptp.append(sptp)
		accounted_list, unaccounted_list, new_index = do_two_way_comparison_single(ptp[0], ptp[1], len(tlist))
		# export_sorting_results(clusters)
		msg  = " total number of accounted for two rsort:  %d percentage:  %5.2f "%(len(accounted_list), float(len(accounted_list))/float(len(tlist))*100.0)
		print(line, msg)
		log_main.add(msg)
		Tracker["accounted_list"] = []
		for index in xrange(len(accounted_list)): Tracker["accounted_list"].append([new_index[index][0], tlist[new_index[index][1]]]) # always check
		rsort_clusters, new_partition = split_partition_into_ordered_clusters(Tracker["accounted_list"])
		Tracker["accounted_list"] = new_partition
		write_text_row(Tracker["accounted_list"], os.path.join(Tracker["constants"]["masterdir"], "final_partition.txt"))
	else:  rsort_clusters = 0
	rsort_clusters = wrap_mpi_bcast(rsort_clusters, Blockdata["main_node"], MPI_COMM_WORLD)
	# save res
	if(Blockdata["myid"] == Blockdata["main_node"]):# save clusters with respect to the initial indexes
		for icluster in xrange(len(rsort_clusters)):
			write_text_file(rsort_clusters[icluster], os.path.join(Tracker["constants"]["masterdir"],"Cluster%d.txt"%icluster))
			Tracker["number_of_groups"] = len(rsort_clusters)
	else:   Tracker["number_of_groups"] = 0
	Tracker["number_of_groups"] = bcast_number_to_all(Tracker["number_of_groups"], Blockdata["main_node"], MPI_COMM_WORLD)
	### Final Rec3D unfiltered two halves, valid only in case of sorting initiated from sphire refinement 
	if Tracker["constants"]["final_sharpen"]:
		Tracker["constants"]["orgres"]				= 0.0
		Tracker["constants"]["refinement_delta"]	= 0.0
		Tracker["constants"]["refinement_ts"]		= 0.0
		Tracker["constants"]["refinement_xr"]		= 0.0
		Tracker["constants"]["refinement_an"]		= 0.0
		minimum_size = Tracker["constants"]["img_per_grp"]
		number_of_groups = 0
		if(Blockdata["myid"] == Blockdata["main_node"]):
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			msg  = "reconstruct unfitered cluster maps for sharpening"
			log_main.add(msg)
			print(line, msg)
			while os.path.exists(os.path.join(Tracker["constants"]["masterdir"], "Cluster%d.txt"%number_of_groups)):
				class_in = read_text_file(os.path.join(Tracker["constants"]["masterdir"], "Cluster%d.txt"%number_of_groups))
				minimum_size = min(len(class_in), minimum_size)
				number_of_groups +=1
				del class_in
			msg = " % clusters are found  "%number_of_groups
			log_main.add(msg)
			print(msg)
		number_of_groups = bcast_number_to_all(number_of_groups, Blockdata["main_node"], MPI_COMM_WORLD)
		if number_of_groups == 0:ERROR("No cluster is found, and the program terminates.", "do_final_maps", 1, Blockdata["myid"])
		minimum_size = bcast_number_to_all(minimum_size, Blockdata["main_node"], MPI_COMM_WORLD)
		compute_noise(Tracker["constants"]["nnxo"])
		do_final_maps(number_of_groups, minimum_size, Tracker["constants"]["selected_iter"], Tracker["constants"]["refinement_dir"], \
		  Tracker["constants"]["masterdir"], Tracker["constants"]["nnxo"], log_main)
	else:
		if(Blockdata["myid"] == Blockdata["main_node"]):
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			msg  = "sort3d finishes without do map sharpening!"
			print(line, msg)
			log_main.add(msg)
	mpi_barrier(MPI_COMM_WORLD)
	if(Blockdata["myid"] == Blockdata["main_node"]):
		for iproc in xrange(number_of_groups):
			merge_two_unfiltered_maps(os.path.join(Tracker["constants"]["masterdir"], "vol_unfiltered_0_grp%03d.hdf"%iproc), \
			  os.path.join(Tracker["constants"]["masterdir"], "vol_unfiltered_1_grp%03d.hdf"%iproc), iproc)
	mpi_barrier(MPI_COMM_WORLD)
	from mpi import mpi_finalize
	mpi_finalize()
	exit()
if __name__ == "__main__":
	main()
