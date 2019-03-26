#!/usr/bin/env python
#
#
#  08/26/2016
#  New version of sort3D.
#  
from __future__ import print_function
from sp_global_def import sxprint
import EMAN2_cppwrap
import sp_applications
import copy
import sp_filter
import sp_fundamentals
import sp_global_def
import json
import sp_logger
import sp_morphology
import mpi
import numpy
import numpy as np
import optparse
import os
import sp_projection
import random
import sp_reconstruction
import scipy
import scipy.stats
import shutil
import sp_statistics
import string
import sys
import time
import sp_utilities

from builtins import range

"""
There are two ways to run the program, and for each way there are three modes to do sorting:

I.a Import data from a meridien refinement. User defines group size, and small groups 
are not removed during within box sorting.

mpirun -np 64 --hostfile four_nodes.txt  sxsort3d_depth.py  --refinemet_dir=meridien_run --niter_for_sorting=28 \
      --memory_per_node=100. --img_per_grp=80000   --output_dir=SORT3D

I.b Import data from a meridien refinement. User defines group size, and small groups 
are removed during within box sorting.

mpirun -np 64 --hostfile four_nodes.txt  sxsort3d_depth.py  --refinemet_dir=meridien_run 
   --niter_for_sorting=28  --memory_per_node=100. --img_per_grp=80000   \
   --output_dir=SORT3D   --check_smearing  --compute_on_the_fly --not_freeze_groups
   
I.c Import data from a meridien refinement. Program searches for clusters without prior information. 
  mpirun -np 64 --hostfile four_nodes.txt  sxsort3d_depth.py  --refinemet_dir=meridien_run --niter_for_sorting=28 \
      --memory_per_node=100. --output_dir=SORT3D
  

II.a Import data from a given data stack. User defines group size, and small groups 
   are not removed during within box sorting. This way uses much less memory than the method I.

mpirun  -np 48  --hostfile ./node012.txt  sxsort3d_depth.py \
      --output_dir=sorting_bmask04 --sym=c1  \
     --radius=30  --mu=1   --img_per_grp=2800    --instack=bdb:data  >sorting_bmask04/printout &

II.b Import data from a given data stack. User defines group size, and small groups 
are removed during within box sorting.

mpirun  -np 48  --hostfile ./node012.txt  sxsort3d_depth.py \
    --output_dir=sorting_bmask04 --sym=c1 --radius=30   --img_per_grp=2800  --not_freeze_groups \
          --instack=bdb:data  >sorting_bmask04/printout &
          
II.c Import data from a given data stack. Program searches for clusters without prior information. 
mpirun  -np 48  --hostfile ./node012.txt  sxsort3d_depth.py \
      --output_dir=sorting_bmask04 --sym=c1 --radius=30 --instack=bdb:data  >sorting_bmask04/printout &

III. Run the program on a single node workstation. The minimum requirement for a successful run on 
     a single workstation is the 2D data size is less than the total memory of the workstation. 

mpirun  -np 8  --hostfile ./node0.txt  sxsort3d_depth.py --orientation_groups=40 \
     --output_dir=sorting_bmask04 --sym=c1  --radius=30  \
         --img_per_grp=2800    --instack=bdb:data  >sorting_bmask04/printout &
         

Frequently used options:
a. --compute_on_the_fly: (valid only for meridien refinement)It enables sorting done with iterations with large 
     number of smearing. The shifted data are computed on the fly.
b. --nstep: Number of steps for sorting shrinks the minimum group size from high bound to low bound. Defalut value is 3.
c. --use_umat: using fuzzy membership to stabilize sorting. 
d. --check_smearing: (valid only for meridien refinement) program prints out the averaged number of smearings in each 
     iteration of meridien refinement


Nomenclatures in sorting intermediate results:
NACC:  total number of accounted
NUACC: number of unaccounted
K:     number of groups
mu:    maximum tolerable resolution discrepancy in Fourier pixels between the largest group and the smallest group
"""

global  Tracker, Blockdata

mpi.mpi_init(0, [])
nproc     = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
myid      = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
Blockdata = {}
#  MPI stuff
Blockdata["nproc"]              = nproc
Blockdata["myid"]               = myid
Blockdata["main_node"]          = 0
Blockdata["last_node"]          = nproc -1

Blockdata["shared_comm"]                    = mpi.mpi_comm_split_type(mpi.MPI_COMM_WORLD, mpi.MPI_COMM_TYPE_SHARED,  0, mpi.MPI_INFO_NULL)
Blockdata["myid_on_node"]                   = mpi.mpi_comm_rank(Blockdata["shared_comm"])
Blockdata["no_of_processes_per_group"]      = mpi.mpi_comm_size(Blockdata["shared_comm"])
masters_from_groups_vs_everything_else_comm = mpi.mpi_comm_split(mpi.MPI_COMM_WORLD, Blockdata["main_node"] == Blockdata["myid_on_node"], Blockdata["myid_on_node"])
Blockdata["color"], Blockdata["no_of_groups"], balanced_processor_load_on_nodes = sp_utilities.get_colors_and_subsets(Blockdata["main_node"], mpi.MPI_COMM_WORLD, Blockdata["myid"], \
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
sp_global_def.BATCH = True
sp_global_def.MPI   = True
##========================================================================================
global _proc_status, _scale, is_unix_cluster
try:			
	_proc_status = '/proc/%d/status' % os.getpid()
	_scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,'KB': 1024.0, 'MB': 1024.0*1024.0}
	is_unix_cluster = True
except:
	if Blockdata["myid"]== Blockdata["main_node"]: sxprint("Not a unix machine")
	is_unix_cluster = False
##========================================================================================	
def create_subgroup():
	# select a subset of myids to be in subdivision
	if( Blockdata["myid_on_node"] < Blockdata["ncpuspernode"] ): submyids = [Blockdata["myid"]]
	else:  submyids = []
	submyids = sp_utilities.wrap_mpi_gatherv(submyids, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	submyids = sp_utilities.wrap_mpi_bcast(submyids,   Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	world_group = mpi.mpi_comm_group(mpi.MPI_COMM_WORLD)
	subgroup = mpi.mpi_group_incl(world_group,len(submyids),submyids)
	Blockdata["subgroup_comm"] = mpi.mpi_comm_create(mpi.MPI_COMM_WORLD, subgroup)
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	Blockdata["subgroup_size"] = -1
	Blockdata["subgroup_myid"] = -1
	if (mpi.MPI_COMM_NULL != Blockdata["subgroup_comm"]):
		Blockdata["subgroup_size"] = mpi.mpi_comm_size(Blockdata["subgroup_comm"])
		Blockdata["subgroup_myid"] = mpi.mpi_comm_rank(Blockdata["subgroup_comm"])
	#  "nodes" are zero nodes on subgroups on the two "node_volume" that compute backprojection
	Blockdata["nodes"] = [Blockdata["node_volume"][0]*Blockdata["ncpuspernode"], Blockdata["node_volume"][1]*Blockdata["ncpuspernode"]]
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	return
	
def create_zero_group():
	# select a subset of myids to be in subdivision, This is a group of all zero IDs on nodes, taken from isac2
	if( Blockdata["myid_on_node"] == 0 ): submyids = [Blockdata["myid"]]
	else:  submyids = []

	submyids = sp_utilities.wrap_mpi_gatherv(submyids,  Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	submyids = sp_utilities.wrap_mpi_bcast(submyids,    Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	world_group = mpi.mpi_comm_group(mpi.MPI_COMM_WORLD)
	subgroup    = mpi.mpi_group_incl(world_group,len(submyids),submyids)
	Blockdata["group_zero_comm"] = mpi.mpi_comm_create(mpi.MPI_COMM_WORLD, subgroup)
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

	Blockdata["group_zero_size"] = -1
	Blockdata["group_zero_myid"] = -1
	if (mpi.MPI_COMM_NULL != Blockdata["group_zero_comm"]):
		Blockdata["group_zero_size"] = mpi.mpi_comm_size(Blockdata["group_zero_comm"])
		Blockdata["group_zero_myid"] = mpi.mpi_comm_rank(Blockdata["group_zero_comm"])
	#  "nodes" are zero nodes on subgroups on the two "node_volume" that compute backprojection
	#Blockdata["nodes"] = [Blockdata["node_volume"][0]*Blockdata["ncpuspernode"], Blockdata["node_volume"][1]*Blockdata["ncpuspernode"]]
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	return
###============================================================================================
######### ================>>>  depth clustering functions <<<==================================
def depth_clustering(work_dir, depth_order, initial_id_file, params, previous_params, log_main):
	global Tracker, Blockdata
	##-----------------------------------------------------------------------------------------
	keys_to_be_saved = [ "current_img_per_grp", "mu", "number_of_groups", "rand_index", "nruns"]
	###----------------
	if(Blockdata["myid"] == Blockdata["main_node"]):
		log_main.add('                        Total data : %d   number of groups: %d                       '%\
		   (Tracker["total_stack"], Tracker["number_of_groups"]))
		if not os.path.exists(os.path.join(work_dir, "layer0")):
			os.mkdir(os.path.join(work_dir, "layer0"))
		partition_per_box_per_layer_list = []
		for iparti in range(0, 2**(depth_order+1), 2):
			partition_per_box_per_layer_list.append([sp_utilities.read_text_file(initial_id_file), None])
	else: partition_per_box_per_layer_list = 0
	partition_per_box_per_layer_list = sp_utilities.wrap_mpi_bcast(partition_per_box_per_layer_list, Blockdata["main_node"])
	keepchecking     = 1
	bad_clustering   = 0
	stat_list        = []
	Tracker["depth"] = 0
	for depth in range(depth_order): #  layers, depth_order = 1 means one layer and two boxes.
		Tracker["no_cluster_generation"]  = False
		time_layer_start = time.time()
		n_cluster_boxes  = 2**(depth_order - depth)
		depth_dir        = os.path.join(work_dir, "layer%d"%depth)
		Tracker["depth"] = depth
		if(Blockdata["myid"] == Blockdata["main_node"]):
			log_main.add('                        Layer %d has %d pairs of independent sorting runs'%(depth,n_cluster_boxes) )
			log_main.add(' ')
			log_main.add(' ')
			if not os.path.exists(depth_dir): 
				os.mkdir(depth_dir)
				keepchecking = 0
				mark_sorting_state(depth_dir, False)
			else: 
				keepchecking = check_sorting_state(depth_dir, keepchecking, log_main)
				if keepchecking == 0: mark_sorting_state(depth_dir, False)
		else: keepchecking = 0
		keepchecking = sp_utilities.bcast_number_to_all(keepchecking, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		## box loop
		if keepchecking == 0:
			checkingbox = 1
			for nbox in range(n_cluster_boxes):
				time_box_start = time.time()
				input_accounted_file   = partition_per_box_per_layer_list[nbox][0]
				input_unaccounted_file = partition_per_box_per_layer_list[nbox][1]
				nbox_dir = os.path.join(depth_dir, "nbox%d"%nbox)
				if(Blockdata["myid"] == Blockdata["main_node"]):
					if not os.path.exists(nbox_dir):
						os.mkdir(nbox_dir)
						checkingbox = 0
						mark_sorting_state(nbox_dir, False)
					else:
						checkingbox = check_sorting_state(nbox_dir, checkingbox, log_main)
						if checkingbox == 0:# found not finished box
							msg ="Pair %d is not finished. Remove it and recompute..."%nbox
							log_main.add(msg)
							shutil.rmtree(nbox_dir)
							os.mkdir(nbox_dir)
							mark_sorting_state(nbox_dir, False)
						else:
							msg ="Pair %d is completed"%nbox
							log_main.add(msg)
				else: checkingbox = 0
				checkingbox = sp_utilities.bcast_number_to_all(checkingbox, Blockdata["main_node"], mpi.MPI_COMM_WORLD) 
				Tracker     = sp_utilities.wrap_mpi_bcast(Tracker,          Blockdata["main_node"], mpi.MPI_COMM_WORLD)
				# Code structure of the box
				if checkingbox == 0:
					bad_box = depth_clustering_box(nbox_dir, input_accounted_file, input_unaccounted_file, \
					    params, previous_params, nbox, log_main)
					if(Blockdata["myid"] == Blockdata["main_node"]):# save parameters per box
						mark_sorting_state(nbox_dir, True, dict = copy_state_variables_from_tracker(keys_to_be_saved))
				else: continue
				if bad_box == 1: bad_clustering = 1
				time_box_finish = (time.time() - time_box_start)//60.
				if(Blockdata["myid"] == Blockdata["main_node"]):
					log_main.add('Box %d of layer %d is done and it costs %f minutes'%(nbox, depth, time_box_finish))
			if bad_clustering !=1:
				partition_per_box_per_layer_list = []
				stop_generation_list             = []
				stat_lists                       = []
				for nbox in range(0,n_cluster_boxes, 2):
					input_box_parti1 = os.path.join(depth_dir, "nbox%d"%nbox,     "partition.txt")
					input_box_parti2 = os.path.join(depth_dir, "nbox%d"%(nbox+1), "partition.txt")
					accounted_list, unaccounted_list, bad_clustering, stop_generation, stat_list = \
					do_boxes_two_way_comparison_mpi(nbox, input_box_parti1, input_box_parti2, depth_order - depth, log_main)
					stop_generation_list.append(stop_generation)
					partition_per_box_per_layer_list.append([accounted_list, unaccounted_list])
					stat_lists.append(stat_list)
					
				if sum(stop_generation_list)>=1:# one cluster case: output the one with largest accounted ptls
					if len(partition_per_box_per_layer_list)>1:
						output_box = 0
						output_acc = len(partition_per_box_per_layer_list[0][0])
						for im in range(len(partition_per_box_per_layer_list)):
							if len(partition_per_box_per_layer_list[im][0]) >output_acc:
								output_box = im
						partition_per_box_per_layer_list = [partition_per_box_per_layer_list[output_box]]
						stat_list = stat_lists[output_box]
					break
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
				
			if(Blockdata["myid"] != Blockdata["main_node"]):Tracker = 0
			Tracker = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			if(Blockdata["myid"] == Blockdata["main_node"]):
				mark_sorting_state(depth_dir, True)
			if(bad_clustering  == 1):  break
			if(stop_generation == 1):  break ### only one cluster survives
		else:
			if(Blockdata["myid"] == Blockdata["main_node"]):
				log_main.add('Layer %d sorting has completed. Recompute two_way comparison'%depth)
			partition_per_box_per_layer_list = []
			stop_generation_list = []
			stat_lists           = []
			for nbox in range(0, n_cluster_boxes,2):
				input_box_parti1 = os.path.join(depth_dir, "nbox%d"%nbox,     "partition.txt")
				input_box_parti2 = os.path.join(depth_dir, "nbox%d"%(nbox+1), "partition.txt")
				accounted_list, unaccounted_list, bad_clustering, stop_generation, stat_list = \
				do_boxes_two_way_comparison_mpi(nbox, input_box_parti1, input_box_parti2, depth_order - depth, log_main)
				stop_generation_list.append(stop_generation)
				partition_per_box_per_layer_list.append([accounted_list, unaccounted_list])
				stat_lists.append(stat_list)
				
			if sum(stop_generation_list)>=1:# one cluster case: output the one with largest accounted ptls
				if len(partition_per_box_per_layer_list)>1:
					output_box = 0
					output_acc = len(partition_per_box_per_layer_list[0][0])
					for im in range(len(partition_per_box_per_layer_list)):
						if len(partition_per_box_per_layer_list[im][0]) >output_acc:output_box = im
					partition_per_box_per_layer_list = [partition_per_box_per_layer_list[output_box]]
					stat_list = stat_lists[output_box]
				break
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			
		time_of_sorting_h,  time_of_sorting_m = get_time(time_layer_start)
		if Blockdata["myid"] == Blockdata["main_node"]:
			log_main.add('            Execution of layer %d took %d hours %d minutes.'%(depth, time_of_sorting_h, time_of_sorting_m))
			log_main.add(' ')
			log_main.add(' ')
			log_main.add('================================================================================================================')
	return partition_per_box_per_layer_list, bad_clustering, stat_list
###---------------------------------------------------------------------------------------
def depth_clustering_box(work_dir, input_accounted_file, \
     input_unaccounted_file, params, previous_params, nbox, log_main):
	global Tracker, Blockdata
	# parameters  When nbox = 0, program does K training;
	#---------------------------------------------
	iter_converge_ratio      = 90.  # percentage
	iter_converge_criterion  = 1.0  
	quick_converge           = 3
	#---------------------------------------------
	box_start        = time.time()
	acc_time         = 0.0
	box_niter        = Tracker["constants"]["box_niter"]
	no_groups_runs   = 0
	nruns_rand_index = 0
	
	if(Blockdata["myid"] == Blockdata["main_node"]):
		if not os.path.exists(work_dir): os.mkdir(work_dir)
		freq_cutoff_dict = {}
		fout = open(os.path.join(work_dir, "freq_cutoff.json"),'w')
		json.dump(freq_cutoff_dict, fout)
		fout.close()
		log_main.add('----------------------------------------------------------------------------------------------------------------' )
		log_main.add('                        Executing pair of quasi-independent sortings, pair number %d'%nbox)
		log_main.add('----------------------------------------------------------------------------------------------------------------' )
	else:Tracker = 0
	Tracker = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD)	
	### -------------->>>>> Initialization<<<<<------------------------------
	ncluster         = 0
	nruns            = 0
	keepgoing        = 1
	converged        = 0
	do_freeze_groups = Tracker["freeze_groups"]
	
	if nbox == 0: iter_mstep  = max(Tracker["nstep"], 5)
	else:         iter_mstep  = max(min(Tracker["nstep"], 5), 3)
	
	####--------------------------------------------------------------------
	number_of_groups_init, total_stack_init, new_assignment, do_recompute = \
	   depth_box_initialization(work_dir, input_accounted_file, \
	       input_unaccounted_file, nbox, log_main)
	       
	NUACC                    = total_stack_init
	NACC                     = 0
	assignment_list          = new_assignment[:]
	total_stack              = total_stack_init
	current_number_of_groups = number_of_groups_init
	unaccounted_list         = new_assignment[:]
	###-----------------------------------------------------------------------------------
	### check recompute:		
	while (keepgoing == 1):### start a run
		Tracker["nruns"] = nruns
		if(Blockdata["myid"] == Blockdata["main_node"]):
			log_main.add('Box %d has been under processing for %f minutes...'%(nbox, (time.time()-box_start)/60.))	
		within_box_run_dir = os.path.join(work_dir, "run%d"%nruns)
		unaccounted_file   = os.path.join(within_box_run_dir, "Unaccounted_from_previous_run.txt")
		if(Blockdata["myid"] == Blockdata["main_node"]):
			time_box_start = time.time()
			os.mkdir(within_box_run_dir)
			sp_utilities.write_text_file(unaccounted_list, unaccounted_file)# new starting point
		if nruns >0:
			assignment_list = create_nrandom_lists(unaccounted_file, current_number_of_groups, 2)
		if(Blockdata["myid"] == Blockdata["main_node"]):
			for indep in range(2):
				sp_utilities.write_text_row(assignment_list[indep], os.path.join(within_box_run_dir, \
				   "independent_index_%03d.txt"%indep))
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		
		#===============>>>>  Iter initialization  <<<====================================
		iter                = 0
		previous_iter_ratio = 0.0
		current_iter_ratio  = 0.0
		converged           = 0
		iter_dir = os.path.join(within_box_run_dir, "iter%d"%iter)
		if Blockdata["myid"] == Blockdata["main_node"]:
			os.mkdir(iter_dir)
			for indep in range(2):
				sp_utilities.write_text_row(assignment_list[indep], \
				    os.path.join(iter_dir, "random_assignment_%03d.txt"%indep))
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		
		iter_id_init_file = os.path.join(iter_dir, "random_assignment_000.txt")
		if Blockdata["myid"] == Blockdata["main_node"]:
			id_list          = sp_utilities.read_text_file(os.path.join(within_box_run_dir, "independent_index_000.txt") , -1)
			total_stack, current_number_of_groups = len(id_list[0]), int((np.unique(np.array(id_list[0]))).shape[0])
		else:
			total_stack              = 0
			current_number_of_groups = 0
		total_stack                  = sp_utilities.bcast_number_to_all(total_stack,              Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		current_number_of_groups     = sp_utilities.bcast_number_to_all(current_number_of_groups, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		# Read raw data
		original_data, norm_per_particle  = read_data_for_sorting(iter_id_init_file, params, previous_params)
		if Tracker["nosmearing"]: parameterstructure  = None
		else: parameterstructure  = read_paramstructure_for_sorting(\
		   iter_id_init_file,Tracker["paramstructure_dict"], Tracker["paramstructure_dir"])
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		
		# ---------------->>>  Estimate image size  <<<-----------------------------------
		Tracker["directory"] = within_box_run_dir
		Tracker["nxinit"], Tracker["freq_fsc143_cutoff"] = get_sorting_image_size(original_data, \
		      iter_id_init_file, current_number_of_groups, log_main)
		      
		if Blockdata["myid"] == Blockdata["main_node"]:
			log_main.add('Sorting cutoff frequency:  %f'%(round(Tracker["freq_fsc143_cutoff"], 4)) +\
			  '  Sorting image size:  %d'%(Tracker["nxinit"]))
		Tracker = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		
		### ----->>> preset variables in trackers <<<-------------------------------------
		Tracker["total_stack"]         = total_stack
		Tracker["number_of_groups"]    = current_number_of_groups
		Tracker["current_img_per_grp"] = total_stack//Tracker["number_of_groups"]
		
		if nruns >0:# reset 
			if Blockdata["myid"] == Blockdata["main_node"]:
				calculate_grp_size_variation(log_main, state="run_init") 
				set_minimum_group_size(log_main, printing_dres_table = True)# When total number of images changes
			else: Tracker = 0
			Tracker = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		###-------------------------------------------------------------------------------	
		iter_previous_iter_ratio = 0.0
		iter_current_iter_ratio  = 0.0
		if nbox == 0: swap_ratio = 2*Tracker["swap_ratio"]
		else:         swap_ratio =   Tracker["swap_ratio"]
		### prepare cdata, fdata, srdata, ctf_images
		cdata, rdata, fdata, ctf_images = downsize_data_for_sorting(original_data, preshift = True, \
		    npad = 1, norms = norm_per_particle)# pay attentions to shifts!
		del original_data
		####------------------------------------------------------------------------------
		# Compute data statistics and estimate permitted number of smearing images computed on fly
		
		if Tracker["constants"]["compute_on_the_fly"]:
			if Blockdata["myid"] == Blockdata["main_node"]:
				precalculated_images_per_cpu = mem_calc_and_output_info(Tracker["constants"]["smearing_file"], \
				 Tracker["nxinit"], iter_id_init_file,log_main)
			else: precalculated_images_per_cpu = 0
			precalculated_images_per_cpu = sp_utilities.wrap_mpi_bcast(precalculated_images_per_cpu, \
				 Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			Tracker["num_on_the_fly"] = precalculated_images_per_cpu
				
		srdata = precalculate_shifted_data_for_recons3D(rdata, parameterstructure, Tracker["refang"], \
	      Tracker["rshifts"], Tracker["delta"], Tracker["avgnorm"], Tracker["nxinit"], \
	      Tracker["constants"]["nnxo"], Tracker["nosmearing"], norm_per_particle, Tracker["constants"]["nsmear"])
		del rdata
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD) # uneven smearing requires
		
		#---------->>>  Initial values per run  <<<<<-------------------------------------
		max_iter                      = Tracker["total_number_of_iterations"]
		minimum_grp_size              = [0, 0] # adjustable
		do_freeze                     = [0, 0]
		iter_init_K                   = Tracker["number_of_groups"]
		last_recompute                = 0
		quick_converge_counter        = 0
		Tracker["check_minimum_size"] = True
		### ------------------------------------------------------------------------------
		
		while((iter < box_niter) and (converged == 0)):# Within a pair of Kmeans, start an iter
			info_table             = None
			reset_number_of_groups = 0
			for indep_run_iter in range(2):
				Tracker["directory"] = os.path.join(iter_dir, "MGSKmeans_%03d"%indep_run_iter)
				MGSKmeans_index_file = os.path.join(iter_dir, "random_assignment_%03d.txt"%indep_run_iter)
				if Blockdata["myid"] == Blockdata["main_node"]:
					os.mkdir(Tracker["directory"])
					os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
				###---------------->>>     Kmeans clustering      <<<----------------------
				tmp_final_list, kmeans_iterations, minimum_grp_size[indep_run_iter], do_freeze[indep_run_iter]=\
				    Kmeans_minimum_group_size_orien_groups(nbox, iter_mstep, iter, cdata, fdata, srdata, ctf_images, \
					parameterstructure, norm_per_particle, MGSKmeans_index_file, params, \
					 minimum_grp_size[indep_run_iter], do_freeze[indep_run_iter], indep_run_iter, \
					    clean_volumes = Tracker["clean_volumes"], global_log=log_main)
				         
				max_iter = min(max_iter, kmeans_iterations)	   
				if Blockdata["myid"] == Blockdata["main_node"]:
					sp_utilities.write_text_row(tmp_final_list, os.path.join(iter_dir,\
					    "partition_%03d.txt"%indep_run_iter))
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			
			if Blockdata["myid"] == Blockdata["main_node"]:
				smallest_size, largest_size, list_of_stable, unaccounted_list, \
				      iter_current_iter_ratio, cngroups, info_table = \
				do_withinbox_two_way_comparison(iter_dir, nbox, nruns, iter) # two partitions are written in partition_dir as partition_%03d.txt
			else: 
				unaccounted_list   =  0
				list_of_stable     =  0
				smallest_size      =  0
				Tracker            =  0
				largest_size       =  0
				cngroups           =  0
			smallest_size       = sp_utilities.bcast_number_to_all(smallest_size, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			largest_size        = sp_utilities.bcast_number_to_all(largest_size,  Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			cngroups            = sp_utilities.bcast_number_to_all(cngroups,      Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			unaccounted_list    = sp_utilities.wrap_mpi_bcast(unaccounted_list,   Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			list_of_stable      = sp_utilities.wrap_mpi_bcast(list_of_stable,     Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			Tracker             = sp_utilities.wrap_mpi_bcast(Tracker,            Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			
			accounted_file      = os.path.join(iter_dir, "Accounted.txt")
			unaccounted_file    = os.path.join(iter_dir, "Core_set.txt")
			
			###### --------->>> Compute rand_index for two independent runs <<<--------------------
			if Blockdata["myid"] == Blockdata["main_node"]:
				core1      = sp_utilities.read_text_file(os.path.join(iter_dir, "partition_000.txt"),  -1)
				core2      = sp_utilities.read_text_file(os.path.join(iter_dir, "partition_001.txt"),  -1)
				core1      = core1[0]
				core2      = core2[0]
			else:
				core1 = 0
				core2 = 0
			core1       = sp_utilities.wrap_mpi_bcast(core1, Blockdata["main_node"])
			core2       = sp_utilities.wrap_mpi_bcast(core2, Blockdata["main_node"])
			rand_index  = compute_rand_index_mpi(core1, core2)
			#### ---------->>>>AI makes decisions for iter run <<<<-----------------------
			if Blockdata["myid"] == Blockdata["main_node"]:
				score_diff = abs(rand_index*100. - iter_current_iter_ratio)
				info_table.append('Rand index of two paritions is %5.4f '%rand_index)
				info_table.append('Difference between two indeces is %5.1f percent'%score_diff)
				######------ Check converge
				if (do_freeze_groups == 0):
					if (iter>=2) and (nbox == 0):#####
						if (abs(iter_current_iter_ratio - iter_previous_iter_ratio < iter_converge_criterion) or \
						   iter_current_iter_ratio > 90.) or (rand_index>0.90):#
							converged = 1 # no much improved
							
						if (max_iter <= quick_converge) and (iter>0):
							quick_converge_counter +=1
							if  quick_converge_counter>=2: 
								converged = 1 # Kmeans goes for only one step.
						####------------- adjust number of groups		
						if (nbox == 0) and (nruns == 0): #
							if (len(list_of_stable) < iter_init_K ) and len(list_of_stable)>2: 
								if (converged == 1) or (iter-last_recompute)>=2:
									reset_number_of_groups = 1
									acc  = 0
									for a in list_of_stable: acc +=len(a)
									unacc          = len(unaccounted_list)
									gsize          = acc//len(list_of_stable)
									funacc         = int(Tracker["total_stack"]*0.01)
									ngrp           = (unacc-funacc)//gsize
									iter_init_K    = len(list_of_stable) + ngrp
									ngrp = (unacc-funacc)//gsize
									if (ngrp >=1) and (iter_init_K < Tracker["number_of_groups"]):
										random.shuffle(unaccounted_list)
										last_recompute = iter
										info_table.append('unacc %d  funacc %d  gsize  %d'%(unacc, funacc, gsize))
										info_table.append('Reset number of groups to %d'%iter_init_K)
										for im in range(ngrp):
											list_of_stable.append(unaccounted_list[0:gsize])
											del unaccounted_list[0:gsize]
										accounted_list, new_index = merge_classes_into_partition_list(list_of_stable)
										sp_utilities.write_text_row(new_index, accounted_file)
										sp_utilities.write_text_file(unaccounted_list, unaccounted_file)
									else:
										iter_init_K = len(list_of_stable) 
										info_table.append('Reset number of groups to %d'%len(list_of_stable))
								else: pass
					elif nbox > 0:
						if (abs(iter_current_iter_ratio - iter_previous_iter_ratio < iter_converge_criterion) or \
						   iter_current_iter_ratio > 90.) or (rand_index>0.90):#
						   if (iter>=2): converged = 1 # n
						if (max_iter <= quick_converge) and (iter>=2):converged = 1
				else:
					# -------->>> freeze groups <<<---------
					if (abs(iter_current_iter_ratio - iter_previous_iter_ratio < iter_converge_criterion) \
					   or iter_current_iter_ratio > 90.):#
					   if iter>=2: converged = 1
					   	
					if (max_iter <= quick_converge) and (iter>=2):quick_converge_counter +=1
						
					if (quick_converge_counter>=1) and iter>=2:
						converged = 1 # Kmeans goes for only one step.
					
				###------------------print info	------------------------------------------	 
				with open(os.path.join(iter_dir, "info.txt"),'w') as fout:
					for line in info_table:
						fout.write(line +'\n')
				fout.close()
				###-----------------------------------------------------------------------
			else: 
				converged              = 0
				reset_number_of_groups = 0
				iter_init_K            = 0
				last_recompute         = 0
				quick_converge_counter = 0
				Tracker                = 0
				
			reset_number_of_groups = sp_utilities.bcast_number_to_all(reset_number_of_groups, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			iter_init_K            = sp_utilities.bcast_number_to_all(iter_init_K,            Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			last_recompute         = sp_utilities.bcast_number_to_all(last_recompute,         Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			converged              = sp_utilities.bcast_number_to_all(converged,              Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			quick_converge_counter = sp_utilities.bcast_number_to_all(quick_converge_counter, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			Tracker                = sp_utilities.wrap_mpi_bcast(Tracker,                     Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			######---------------------- Minimum group size checking ---------------------
			"""Multiline Comment0"""
			####------------------------------------------------------------>>> Case 1
			if (converged == 0) and (iter < box_niter+1) and (reset_number_of_groups == 0):
				new_assignment_list = []
				for indep in range(2):
					tmp_list = swap_accounted_with_unaccounted_elements_mpi(accounted_file, \
					 unaccounted_file, log_main, current_number_of_groups, swap_ratio)
					new_assignment_list.append(tmp_list)
				iter_previous_iter_ratio = iter_current_iter_ratio
				iter                    +=1
				iter_dir                 = os.path.join(within_box_run_dir, "iter%d"%iter)
				###++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				if Blockdata["myid"] == Blockdata["main_node"]:
					os.mkdir(iter_dir)
					for indep in range(2):
						sp_utilities.write_text_file(new_assignment_list[indep], \
						os.path.join(iter_dir, "random_assignment_%03d.txt"%indep))
						
			if (reset_number_of_groups == 1): #----------------------------->>> Case 2
				minimum_grp_size         = [0, 0] # adjustable
				do_freeze                = [0, 0]
				box_niter               += 1
				iter_previous_iter_ratio = iter_current_iter_ratio
				
				iter +=1
				iter_dir = os.path.join(within_box_run_dir, "iter%d"%iter)
				if converged ==1: converged =  0
				reset_number_of_groups = 0
				quick_converge_counter = 0
				##+++++++++++++++++++++++++++++++++++++++++++++
				if Blockdata["myid"] == Blockdata["main_node"]:
					log_main.add('Decrease number of groups within box run to %d'%len(list_of_stable))
					Tracker["number_of_groups"]    = len(list_of_stable) # Excessively large number of groups
					Tracker["current_img_per_grp"] = Tracker["total_stack"]//Tracker["number_of_groups"]
					calculate_grp_size_variation(log_main, state = "number_of_groups")
					set_minimum_group_size(log_main, printing_dres_table = True)
					os.mkdir(iter_dir)
				else: Tracker = 0
				Tracker             = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
				new_assignment_list = []
				for indep in range(2):
					tmp_list = swap_accounted_with_unaccounted_elements_mpi(accounted_file, unaccounted_file, \
					      log_main, Tracker["number_of_groups"], swap_ratio)
					new_assignment_list.append(tmp_list)
				if Blockdata["myid"] == Blockdata["main_node"]:
					for indep in range(2):
						sp_utilities.write_text_file(new_assignment_list[indep], \
						  os.path.join(iter_dir, "random_assignment_%03d.txt"%indep))     
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			if (converged ==1) and (reset_number_of_groups == 0):# Case 3 --------------->>>>  converged and finished
				break
		del cdata
		del srdata
		if Tracker["constants"]["focus3D"]: del fdata
		del ctf_images
		if not Tracker["nosmearing"]: del parameterstructure
		del assignment_list
		nruns_rand_index +=rand_index
		####-------------------->>>>>>   Run run part  <<<<<<<----------------------------
		if Blockdata["myid"] == Blockdata["main_node"]:
			ncluster, NACC, NUACC, unaccounted_list, new_clusters = \
			  output_iter_results(work_dir, ncluster, NACC, NUACC, Tracker["minimum_grp_size"],\
			    list_of_stable, unaccounted_list, log_main)
			if os.path.exists(os.path.join(work_dir, 'run%d'%nruns, 'tempdir')):
				shutil.rmtree(os.path.join(work_dir, 'run%d'%nruns, 'tempdir'))
		else:
			ncluster = 0
			NACC     = 0
			NUACC    = 0
			unaccounted_list = 0
			new_clusters     = 0
			Tracker          = 0
		new_clusters     = sp_utilities.bcast_number_to_all(new_clusters,  Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		ncluster         = sp_utilities.bcast_number_to_all(ncluster,      Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		NACC             = sp_utilities.bcast_number_to_all(NACC,          Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		NUACC            = sp_utilities.bcast_number_to_all(NUACC,         Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		unaccounted_list = sp_utilities.wrap_mpi_bcast(unaccounted_list,   Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		Tracker          = sp_utilities.wrap_mpi_bcast(Tracker,            Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		###-------------------------------------------------------------------------------------------
		if new_clusters >0: 
			no_cluster     = False
			no_groups_runs = 0
		else:               
			no_cluster      = True
			no_groups_runs +=1
		###----------------------->>> Output the final info table to log file <<<<<------------------------------ 		
		if Blockdata["myid"] == Blockdata["main_node"]:# report current state
			if new_clusters > 0:
				log_main.add(' ')
				log_main.add('In phase %d, the program found %d groups.'%(nruns, new_clusters))
				for itable in range(len(info_table)): 
					log_main.add(info_table[itable])
				log_main.add("Check the info.txt file in each iter directory for details of within box comparison results. \n")
			keepgoing, nruns, total_stack, current_number_of_groups = \
			   check_state_within_box_run(keepgoing, nruns, unaccounted_list, no_cluster, log_main)
		else:
			Tracker     = 0
			keepgoing   = 0
			nruns       = 0
			total_stack = 0
			current_number_of_groups = 0
		Tracker     = sp_utilities.wrap_mpi_bcast(Tracker,          Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		nruns       = sp_utilities.bcast_number_to_all(nruns,       Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		total_stack = sp_utilities.bcast_number_to_all(total_stack, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		keepgoing   = sp_utilities.bcast_number_to_all(keepgoing,   Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		current_number_of_groups = \
	     sp_utilities.bcast_number_to_all(current_number_of_groups, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		if no_groups_runs >= 2: 
			bad_clustering = 1
			break
	########----------------->>>> Box part <<<<<------------------------------------------
	if Blockdata["myid"] == Blockdata["main_node"]:
		partition = get_box_partition(work_dir, ncluster, unaccounted_list)
	else: partition = 0
	partition = sp_utilities.wrap_mpi_bcast(partition, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	
	if(Blockdata["myid"] == Blockdata["main_node"]):
		if ncluster > 0:
			log_main.add('In independent box run  %d, the program found %d group.'%(nbox, ncluster))
			bad_clustering = 0
		else:
			log_main.add('In independent box run  %d, no reproducible groups were found '%nbox)
			bad_clustering = 1
			partition      =['']
		sp_utilities.write_text_row(partition, os.path.join(work_dir, "partition.txt"))
		Tracker["number_of_groups"]  =  ncluster
		Tracker["rand_index"]        = nruns_rand_index/(nruns+1) # nruns start from 0
	else: 
		bad_clustering = 0
		Tracker        = 0
	bad_clustering = sp_utilities.bcast_number_to_all(bad_clustering, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	Tracker        = sp_utilities.wrap_mpi_bcast(Tracker,             Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	return bad_clustering
###----------------------------------------------end of box clustering
def depth_box_initialization(box_dir, input_list1, input_list2, nbox, log_file):
	# input_list1:  accounted partition
	# input_list2:  unaccounted
	global Tracker, Blockdata
	#---------->>>> parameters<<<<--------------------------------------------------------
	learning_start_rindex  = 0.85 # 
	do_recompute           = 0
	Tracker["rand_index"]  = 0
	if not Tracker["search_mode"]:
		Tracker["current_img_per_grp"] = Tracker["constants"]["img_per_grp"]
	else:
		try: current_img_per_grp = Tracker["current_img_per_grp"]
		except: 
			Tracker["current_img_per_grp"] = \
			   Tracker["constants"]["total_stack"]//Tracker["constants"]["init_K"]
	#####---------------------------------------------------------------------------------
	swap_ratio = Tracker["swap_ratio"]
	if input_list2 is not None: #Tracker 2; layer 1
		total_stack = len(input_list1)+ len(input_list2)
		Tracker["total_stack"] = total_stack
		groups                 = np.unique((np.array(input_list1, dtype=np.int32)).transpose()[0])
		number_of_groups       = groups.shape[0]
		img_per_grp            = total_stack//number_of_groups
		Tracker["current_img_per_grp"] = img_per_grp
		if Blockdata["myid"] == Blockdata["main_node"]:
			msg = "Number of groups found in initialization: %d, the largest possible number of groups: %d"%\
			  (number_of_groups, total_stack//img_per_grp)
			log_file.add(msg)
			sp_utilities.write_text_row(input_list1,  os.path.join(box_dir, "previous_NACC.txt"))
			sp_utilities.write_text_file(input_list2, os.path.join(box_dir, "previous_NUACC.txt"))
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		new_assignment = []
		for indep in range(2):
			tmp_assignment = swap_accounted_with_unaccounted_elements_mpi(os.path.join(box_dir,\
			    "previous_NACC.txt"),os.path.join(box_dir, \
			      "previous_NUACC.txt"), log_file, number_of_groups, swap_ratio)
			new_assignment.append(tmp_assignment)
			
		if Blockdata["myid"] == Blockdata["main_node"]:
			for indep in range(2):
				sp_utilities.write_text_file(new_assignment[indep], \
				  os.path.join(box_dir, "independent_index_%03d.txt"%indep))
			new_assignment = []
			for indep in range(2):
				new_assignment.append(sp_utilities.read_text_row(os.path.join(box_dir,"independent_index_%03d.txt"%indep)))
			id_list = sp_utilities.read_text_file( os.path.join(box_dir, "independent_index_000.txt"), -1)
			if len(id_list)>1: id_list=id_list[0]
			total_stack = len(id_list)
			group_id = np.sort(np.unique(np.array(id_list, dtype=np.int32)))
			number_of_groups = group_id.shape[0] # assume K be 0, ...,number_of_groups-1
			Tracker["total_stack"]         = total_stack
			Tracker["number_of_groups"]    = number_of_groups
			Tracker["current_img_per_grp"] = total_stack//number_of_groups
			Tracker["img_per_grp"]         = [0, 0, 0, 0, 0]
			calculate_grp_size_variation(log_file, state="box_init")
			set_minimum_group_size(log_file, printing_dres_table = False)
		else:
			number_of_groups = 0
			total_stack      = 0
			new_assignment   = 0
			Tracker = 0
		new_assignment   = sp_utilities.wrap_mpi_bcast(new_assignment,        Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		number_of_groups = sp_utilities.bcast_number_to_all(number_of_groups, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		total_stack      = sp_utilities.bcast_number_to_all(total_stack,      Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		Tracker          = sp_utilities.wrap_mpi_bcast(Tracker,               Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		
	else: #Tracker 1; layer 0;
		if(Tracker["box_learning"]==1):# applicable only for not_freeze_groups
			total_stack            = len(input_list1)
			Tracker["total_stack"] = total_stack
			number_of_groups       = total_stack//Tracker["current_img_per_grp"]
			if number_of_groups <=1:
				sp_global_def.ERROR("Number_of_groups should be at least larger than 1", \
				        "depth_box_initialization", 1, Blockdata["myid"])
				
			Tracker["current_img_per_grp"] = total_stack//number_of_groups
			if (nbox >= 1):
				if(Blockdata["myid"] == Blockdata["main_node"]):
					with open(os.path.join(box_dir, "../", "nbox%d"%(nbox-1), "state.json"),'r') as fout:
						dict = sp_utilities.convert_json_fromunicode(json.load(fout))
					fout.close()
					del dict["done"]
					if (dict["rand_index"] > learning_start_rindex):
						Tracker.update(dict)# speed up computation
						total_stack            = len(input_list1)
						Tracker["total_stack"] = total_stack
						number_of_groups       = Tracker["number_of_groups"]
						Tracker["current_img_per_grp"] = Tracker["total_stack"]//Tracker["number_of_groups"]
						log_file.add("In the last box clustering, determined number_of_groups = %d ; current_group_size = %d"%\
						    (Tracker["number_of_groups"], Tracker["current_img_per_grp"]))
						calculate_grp_size_variation(log_file, state="number_of_groups")
						set_minimum_group_size(log_file, printing_dres_table = False)
					else: 
						do_recompute = 1
					sp_utilities.write_text_file(input_list1, os.path.join(box_dir, "previous_all_indexes.txt"))
					log_file.add("Number of groups (K) is inhereted from the last box clustering")
				else: Tracker = 0
				Tracker          = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
				do_recompute     = sp_utilities.bcast_number_to_all(do_recompute, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
				number_of_groups = Tracker["number_of_groups"]
				total_stack      = Tracker["total_stack"]
				new_assignment   = create_nrandom_lists_from_given_pids(box_dir, os.path.join(box_dir, \
					 "previous_all_indexes.txt"), number_of_groups, 2)
			
			if ((nbox == 0) or (do_recompute ==1)):
				if number_of_groups<2: 
					sp_global_def.ERROR("Number_of_groups should be at least larger than 1", "depth_box_initialization", 1, Blockdata["myid"])
				Tracker["total_stack"]         = total_stack
				Tracker["number_of_groups"]    = number_of_groups
				Tracker["current_img_per_grp"] = total_stack//number_of_groups
				Tracker["img_per_grp"]         = [0, 0, 0, 0, 0]
				if (Blockdata["myid"] == Blockdata["main_node"]):	
					calculate_grp_size_variation(log_file, state="box_init")
					set_minimum_group_size(log_file, printing_dres_table = False)
					sp_utilities.write_text_file(input_list1, os.path.join(box_dir, "previous_all_indexes.txt"))
				else: Tracker = 0
				Tracker        = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
				new_assignment = create_nrandom_lists_from_given_pids(box_dir, os.path.join(box_dir, \
					 "previous_all_indexes.txt"), number_of_groups, 2)
		else:
			total_stack            = len(input_list1)
			number_of_groups       = total_stack//Tracker["current_img_per_grp"]
			Tracker["total_stack"] = total_stack
			if number_of_groups<2: 
				sp_global_def.ERROR("Number_of_groups should be at least larger than 1",\
				     "depth_box_initialization", 1, Blockdata["myid"])
			Tracker["number_of_groups"]    = number_of_groups
			Tracker["current_img_per_grp"] = total_stack//number_of_groups
			Tracker["img_per_grp"]         = [0, 0, 0, 0, 0]
			if (Blockdata["myid"] == Blockdata["main_node"]):	
				calculate_grp_size_variation(log_file, state="box_init")
				set_minimum_group_size(log_file, printing_dres_table = False)
				sp_utilities.write_text_file(input_list1, os.path.join(box_dir, "previous_all_indexes.txt"))
			else: Tracker = 0
			Tracker = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			new_assignment = create_nrandom_lists_from_given_pids(box_dir, os.path.join(box_dir, \
		     	 "previous_all_indexes.txt"), number_of_groups, 2)
		      
	if Blockdata["myid"] == Blockdata["main_node"]: # Final printout 
		log_file.add('---->>>>Parameters used in the box sorting<<<----')  
		log_file.add('Number of images:       %d'%total_stack)
		log_file.add('Number of groups:       %d'%number_of_groups)
	return number_of_groups, total_stack, new_assignment, do_recompute
####--------------------------------------- end of box initialization
def sort3d_init(to_be_decided, log_main):
	global Tracker, Blockdata
	keepsorting = 1
	if (Tracker["constants"]["img_per_grp"]<= 2) and (not Tracker["search_mode"]):
		log_main.add("Number of images per group is too small.")
		keepsorting = 0
	if Tracker["total_stack"] <= Blockdata["nproc"]*Tracker["minimum_img_per_cpu"]:
		log_main.add("Too few images are assigned to one single CPU.")
		keepsorting = 0
	return keepsorting
			
def check_state_within_box_run(keepgoing, nruns, unaccounted_list, \
        no_cluster_last_run, log_file):
	global Tracker, Blockdata
	### old AI of box clustering 
	total_stack      = len(unaccounted_list)
	number_of_groups = total_stack//max(Tracker["constants"]["img_per_grp"],\
			Tracker["current_img_per_grp"])
	if number_of_groups>=2: 
		keepgoing = 1
		Tracker["current_img_per_grp"] = total_stack//number_of_groups
		nruns +=1
	else:  
		keepgoing        = 0
		total_stack      = 0
		number_of_groups = 0
	if no_cluster_last_run:  number_of_groups -=1 # Old rules. Now it rarely goes to the next run.
	if number_of_groups <=1: keepgoing = 0
	if not Tracker["search_mode"]:
		if nruns> min(Tracker["constants"]["total_stack"]\
			 //Tracker["constants"]["img_per_grp"], 5):
				keepgoing = 0
				log_file.add("Sort3d stops because of number of runs exceeds desired number of iterations")
	if number_of_groups >2:
		log_file.add("Current parameters for the next run: %d, %d, %d"%\
		   (number_of_groups, total_stack, Tracker["current_img_per_grp"]))
	else:
		log_file.add("Box clustering converges at %d th run. "%nruns) 
	return keepgoing, nruns, total_stack, number_of_groups
	
def get_box_partition(box_dir, ncluster, unaccounted_list):
	unaccounted_list = sorted(unaccounted_list)
	if ncluster >=1:
		clusters_in_box = [sp_utilities.read_text_file(os.path.join(box_dir, \
		      "Cluster_%03d.txt"%ic)) for ic in range(ncluster)]
		if len(unaccounted_list)>0: clusters_in_box.append(unaccounted_list)
		alist, plist = merge_classes_into_partition_list(clusters_in_box)
	else:  plist = []
	return plist
						
def check_sorting_state(current_dir, keepchecking, log_file):
	try:
		with open(os.path.join(current_dir, "state.json"),'r') as fout:
			current_state = sp_utilities.convert_json_fromunicode(json.load(fout))
		fout.close()
		if current_state["done"]:keepchecking = 1
		else: keepchecking = 0
	except:   keepchecking = 0
	return    keepchecking

def mark_sorting_state(current_dir, sorting_done, dict = None):
	# single processor job
	with open(os.path.join(current_dir, "state.json"),'w') as fout:
		if dict is None: dict = {}
		if sorting_done: dict["done"] = True
		else:            dict["done"] = False
		json.dump(dict, fout)
	fout.close()
	return
	
def dump_tracker(path_of_the_tracker):
	global Tracker, Blockdata
	if(Blockdata["myid"] == Blockdata["main_node"]):
		with open(os.path.join(path_of_the_tracker, "Tracker.json"),'w') as fout:
			fout.seek(0)
			json.dump(Tracker, fout)
			fout.truncate()
		fout.close()
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	return
	
def copy_state_variables_from_tracker(state_keys):
	global Tracker, Blockdata
	dict = {}
	for akey in state_keys: 
		dict[akey] = Tracker[akey] 
	return dict

def read_tracker_mpi(current_dir):
	global Tracker, Blockdata
	open_tracker = 1
	if(Blockdata["myid"] == Blockdata["main_node"]):
		try:
			with open(os.path.join(current_dir, "Tracker.json"),'r') as fout:
				Tracker = sp_utilities.convert_json_fromunicode(json.load(fout))
			fout.close()
		except: open_tracker = 0
	else: open_tracker = 0
	open_tracker = sp_utilities.bcast_number_to_all(open_tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	###-----------------------------------------------------------------------------------
	if open_tracker ==1:
		if(Blockdata["myid"] != Blockdata["main_node"]): Tracker = 0
		Tracker = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	else: sp_global_def.ERROR("Fail to load tracker", "%s"%current_dir, 1, Blockdata["myid"])
	return
###---------------------------------------------------------------------------------------
def check_sorting(total_data, keepsorting, log_file):
	global Tracker, Blockdata
	########## add some code to check current settings
	if Blockdata["myid"] == Blockdata["main_node"]:
		if not Tracker["search_mode"]:
			number_of_groups = int(round(total_data/Tracker["constants"]["img_per_grp"]))
			log_file.add("Check_sorting: In the next generation, total =  %d   K = %d  img_per_grp   %d"%\
			     (total_data, number_of_groups, Tracker["constants"]["img_per_grp"]))
		else:
			current_img_per_grp = Tracker["constants"]["total_stack"]//Tracker["constants"]["init_K"]
			number_of_groups    = int(round(total_data/current_img_per_grp))
			
			log_file.add("Check_sorting: In the next generation, total = %8d   K = %3d  img_per_grp = %8d"%\
			     (total_data, number_of_groups, current_img_per_grp))
		
		if number_of_groups>=2:
			Tracker["current_img_per_grp"] = total_data//number_of_groups
			Tracker["total_stack"]         = total_data
			Tracker["number_of_groups"]    = number_of_groups
			keepsorting                    = sort3d_init("Initialization", log_file)
			log_file.add("Check_sorting: SORT3D continues to the next generation. ") 
		else:
			log_file.add("Check_sorting: SORT3D is completed. ") 
			keepsorting = 0
	keepsorting = sp_utilities.bcast_number_to_all(keepsorting, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	Tracker     = sp_utilities.wrap_mpi_bcast(Tracker,          Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	return keepsorting
###---------------------------------------------------------------------------------------			
def print_matching_pairs(pair_list, log_file):
	log_file.add(' ')
	log_file.add('                        Two-way matching of sorting results                ')
	log_file.add('M indicates that respective group of P0 sorting (row number) matches respective group of P1 sorting (column number)')
	msg ='   '
	for i in range(len(pair_list)):
		msg += '{:^5d}'.format(i)
	log_file.add(msg)
	for im in range(len(pair_list)):
		msg ='{:^3d}'.format(im)
		for jm in range(len(pair_list)):
			not_found = True
			for km in range(len(pair_list)):
				if pair_list[km][0] == im and pair_list[km][1] == jm:
					msg +='{:^5s}'.format('M')
					not_found = False
			if not_found: msg +='{:^5s}'.format(' ')
		log_file.add(msg)
	return

def create_masterdir():
	global Tracker, Blockdata
	masterdir = Tracker["constants"]["masterdir"]
	restart   = 0
	if(Blockdata["myid"] == Blockdata["main_node"]):
		if not os.path.exists(os.path.join(Tracker["constants"]["masterdir"])):
			if not masterdir:
				timestring = time.strftime("_%d_%b_%Y_%H_%M_%S", time.localtime())
				masterdir  ="sort3d"+timestring
				os.makedirs(masterdir)
			else:
				if not os.path.exists(masterdir):
					os.makedirs(masterdir)
			li =len(masterdir)
		else:
			li =len(masterdir)
			if os.path.exists(os.path.join(Tracker["constants"]["masterdir"], "generation_000")): restart = 1
			if os.path.exists(os.path.join(Tracker["constants"]["masterdir"], "Tracker.json")):   restart = 1
		sp_global_def.write_command(masterdir)
	else: 
		restart = 0
		li = 0
	restart   = sp_utilities.bcast_number_to_all(restart,       Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	li        = mpi.mpi_bcast(li,       1,   mpi.MPI_INT,  Blockdata["main_node"], mpi.MPI_COMM_WORLD)[0]
	masterdir = mpi.mpi_bcast(masterdir,li,  mpi.MPI_CHAR, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	masterdir = string.join(masterdir,"")
	if (not Tracker["constants"]["masterdir"]): Tracker["constants"]["masterdir"]  = masterdir
	Tracker["constants"]["chunk_0"]  = os.path.join(Tracker["constants"]["masterdir"],"chunk_0.txt")
	Tracker["constants"]["chunk_1"]  = os.path.join(Tracker["constants"]["masterdir"],"chunk_1.txt")
	return restart
	
def print_shell_command(args_list, log_main):
	global Tracker, Blockdata
	if(Blockdata["myid"] == Blockdata["main_node"]):
		log_main.add("The shell line command:")
		line = ""
		for a in args_list: line +=(a + " ")
		log_main.add(line)
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	return
	
def get_time(time_start):
	current_time   = time.time() - time_start
	current_time_h = current_time //3600
	current_time_m = (current_time - current_time_h*3600)//60
	return int(current_time_h), int(current_time_m)
###=======================================================================================
###===================>>> Image processing utilities <<<<<<===============================	
def get_sorting_image_size(original_data, partids, number_of_groups, log_main):
	global Tracker, Blockdata
	
	if Tracker["fixed_sorting_size"]:
		if Tracker["search_mode"]:
			if (Tracker["constants"]["img_per_grp"] <1):
				img_per_grp  = Tracker["constants"]["total_stack"]//Tracker["constants"]["init_K"]
			else:
				img_per_grp  = Tracker["constants"]["img_per_grp"]
			avg_fsc = sp_statistics.scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], \
		 		float(Tracker["constants"]["total_stack"]), img_per_grp)
	else:
		avg_fsc = sp_statistics.scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], \
		 	float(Tracker["constants"]["total_stack"]), Tracker["total_stack"]//number_of_groups)
		 	
	fsc143 = get_res143(avg_fsc)
	if (fsc143 !=0): 
		nxinit = sp_fundamentals.smallprime(min((int(fsc143)+ max(int(Tracker["constants"]["nnxo"]*0.03), 5))*2, \
		     Tracker["constants"]["nnxo"]))
		nxinit = sp_fundamentals.smallprime(nxinit-nxinit%2)
		nxinit = nxinit - nxinit%2
	else:
		sp_global_def.ERROR("Program obtains wrong image size", "get_sorting_image_size", 1, Blockdata["myid"])
		
	if (Tracker["nosmearing"]) and (Tracker["constants"]["nxinit"] !=-1):
		nxinit = Tracker["constants"]["nxinit"]
		if (nxinit >Tracker["constants"]["nnxo"]):
			sp_global_def.ERROR("User provides wrong nxinit", "get_sorting_image_size", 1, Blockdata["myid"])
	Tracker["nxinit"] = nxinit
	compute_noise(Tracker["nxinit"])
	freq_fsc143_cutoff = float(fsc143)/float(nxinit)
	if(Blockdata["myid"] == Blockdata["main_node"]): 
		sp_utilities.write_text_file(avg_fsc, os.path.join(Tracker["directory"], "fsc_image_size.txt"))
	del avg_fsc
	return nxinit, freq_fsc143_cutoff
	
def compute_noise(image_size):
	global Tracker, Blockdata
	if Tracker["applybckgnoise"]: # from SPARX refinement only
		if(Blockdata["myid"] == Blockdata["main_node"]):
			tsd = sp_utilities.get_im(Tracker["bckgnoise"]) # inverted power spectrum
			nnx = tsd.get_xsize()
			nny = tsd.get_ysize()
		else:
			nnx = 0
			nny = 0
		nnx = sp_utilities.bcast_number_to_all(nnx, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		nny = sp_utilities.bcast_number_to_all(nny, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		if( Blockdata["myid"] != Blockdata["main_node"]):tsd = sp_utilities.model_blank(nnx, nny)
		sp_utilities.bcast_EMData_to_all(tsd, Blockdata["myid"], Blockdata["main_node"])
		temp_image = sp_utilities.model_blank(image_size, image_size)
		temp_image = sp_fundamentals.fft(temp_image)
		nx = temp_image.get_xsize()
		ny = temp_image.get_ysize()
		Blockdata["bckgnoise"]  = []
		Blockdata["unrolldata"] = []
		for im in range(nny):
			prj = nnx*[0.0]
			for k in range(nnx):
				if tsd.get_value_at(k,im)>0.0: prj[k] = tsd.get_value_at(k,im)
			Blockdata["bckgnoise"].append(prj)
		for im in range(len(Blockdata["bckgnoise"])): 
			Blockdata["unrolldata"].append(EMAN2_cppwrap.Util.unroll1dpw(ny, ny, Blockdata["bckgnoise"][im]))
	else: # from datastack and relion
		temp_image = sp_utilities.model_blank(image_size, image_size)
		temp_image = sp_fundamentals.fft(temp_image)
		nx = temp_image.get_xsize()
		ny = temp_image.get_ysize()
		Blockdata["bckgnoise"]  = [1.0]*nx
		Blockdata["unrolldata"] = EMAN2_cppwrap.Util.unroll1dpw(ny, ny, nx*[1.0])
	return
	
def check_3dmask(log_main):
	global Tracker, Blockdata
	Tracker["nxinit"]      = Tracker["nxinit_refinement"]
	#   shrinkage, current resolution, fuse_freq
	Tracker["total_stack"] = Tracker["constants"]["total_stack"]
	Tracker["shrinkage"]   = float(Tracker["nxinit"])/Tracker["constants"]["nnxo"]
	Tracker["radius"]      = Tracker["constants"]["radius"]*Tracker["shrinkage"]
	try: fuse_freq = Tracker["fuse_freq"]
	except: Tracker["fuse_freq"] = int(Tracker["constants"]["pixel_size"]\
	     *Tracker["constants"]["nnxo"]/Tracker["constants"]["fuse_freq"]+0.5)
	Tracker = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	dump_tracker(Tracker["constants"]["masterdir"])
	Tracker["shrinkage"] = float(Tracker["nxinit"])/Tracker["constants"]["nnxo"]
	
	if Tracker["constants"]["focus3D"]:
		if(Blockdata["myid"] == Blockdata["main_node"]):
			log_main.add(' ')
			log_main.add('Sort3d is run in focus mode, and focus mask is %s'%\
			   Tracker["constants"]["focus3D"]+'\n')
	return
	
def estimate_tanhl_params(cutoff, taa, image_size):
	def tanhfl(x, cutoff, taa):
		omega = cutoff
		cnst  = numpy.pi/(2.0*omega*taa)
		v1    = (cnst*(x + omega))
		v2    = (cnst*(x - omega))
		return 0.5*(numpy.tanh(v1) - numpy.tanh(v2))
	
	def get_filter(cutoff1, taa, image_size):
		values = []
		N = image_size//2
		for im in range(N):
			x = float(im)/float(image_size)
			values.append(tanhfl(x, cutoff1, taa))
		return values

	values  = get_filter(cutoff, taa, image_size)
	icutoff = image_size
	init    = int(cutoff*image_size)
	while icutoff>=cutoff*image_size:
		cutoff1 = float(init)/image_size
		values  = get_filter(cutoff1, taa, image_size)
		for im in range(len(values)):
			if values[im] <=0.0:
				icutoff = im
				break
		init -=1
	return cutoff1, taa
#####=====================================================================================
###### ==================>>> Statistical analysis <<<=====================================
def do_analysis_on_identified_clusters(clusters, log_main):
	global Tracker, Blockdata
	if Tracker["nosmearing"]:
		ds, ss, norms = get_params_for_analysis(Tracker["constants"]["orgstack"], \
				os.path.join(Tracker["constants"]["masterdir"],"refinement_parameters.txt"),\
				None, None)
	else:
		ds, ss, norms = get_params_for_analysis(Tracker["constants"]["orgstack"], \
			os.path.join(Tracker["constants"]["masterdir"],"refinement_parameters.txt"),\
			os.path.join(Tracker["constants"]["masterdir"], "all_smearing.txt"), \
			       Tracker["constants"]["nsmear"])

	tmpres1, tmpres2 = do_one_way_anova_scipy(clusters, ds,                       name_of_variable="defocus",  log_main = log_main)
	if ss is not None: tmpres1, tmpres2 = do_one_way_anova_scipy(clusters, ss,    name_of_variable="smearing", log_main = log_main)
	if norms:          tmpres1, tmpres2 = do_one_way_anova_scipy(clusters, norms, name_of_variable="norm",     log_main = log_main)
	return
			
def get_params_for_analysis(orgstack, ali3d_params, smearing_file, smearing_number):
	if ali3d_params is not None:
		ali3d_params   = sp_utilities.read_text_file(ali3d_params, -1)
		try: norm_list = ali3d_params[7]
		except: norm_list = None
	if( orgstack is not None ):
		defo_list = []
		ctfs      = EMAN2_cppwrap.EMUtil.get_all_attributes(orgstack, "ctf")
		for im in range(len(ctfs)): defo_list.append(ctfs[im].defocus)
	else: defo_list = None
	if smearing_file is not None:
		try: smearing_list = sp_utilities.read_text_file(smearing_file)
		except: smearing_list = None
	else: smearing_list = None
	return defo_list, smearing_list, norm_list

def do_one_way_anova_scipy(clusters, value_list, name_of_variable="variable", log_main = None):
	# single cpu program
	NMAX = 30
	log_main.add('================================================================================================================')
	log_main.add('                                       ANOVA analysis')
	log_main.add('----------------------------------------------------------------------------------------------------------------')
	if len(clusters)<=1:
		log_main.add("There is only one group. No %s ANOVA analysis"%name_of_variable)
		return None, None, None
	K = min(NMAX, len(clusters))
	replicas = []
	for ic in range(K):
		ll = copy.deepcopy(clusters[ic])
		ll1 = [None for i in range(len(ll))]
		for ie in range(len(ll)): 
			ll1[ie] = value_list[ll[ie]]
		if name_of_variable == "defocus": ll1 = list(set(ll1))
		replicas.append(ll1)
	x0 = replicas[0]
	x1 = replicas[1]
	try: x2 = replicas[2]
	except: pass
	try: x3 = replicas[3]
	except: pass
	try: x4 = replicas[4]
	except: pass
	try: x5 = replicas[5]
	except: pass
	try: x6 = replicas[6]
	except: pass
	try: x7 = replicas[7]
	except: pass
	try: x8 = replicas[8]
	except: pass
	try: x9 = replicas[9]
	except: pass
	try: x10 = replicas[10]
	except: pass
	try: x11 = replicas[11]
	except: pass
	try: x12 = replicas[12]
	except: pass
	try: x13 = replicas[13]
	except: pass
	try: x14 = replicas[14]
	except: pass
	try: x15 = replicas[15]
	except: pass
	try: x16 = replicas[16]
	except: pass
	try: x17 = replicas[17]
	except: pass
	try: x18 = replicas[18]
	except: pass
	try: x19 = replicas[19]
	except: pass
	try: x20 = replicas[20]
	except: pass
	try: x21 = replicas[21]
	except: pass
	try: x22 = replicas[22]
	except: pass
	try: x23 = replicas[23]
	except: pass
	try: x24 = replicas[24]
	except: pass
	try: x25 = replicas[25]
	except: pass
	try: x26 = replicas[26]
	except: pass
	try: x27 = replicas[27]
	except: pass
	try: x28 = replicas[28]
	except: pass
	try: x29 = replicas[29]
	except: pass
	
	if   K==2: res = scipy.stats.f_oneway(x0, x1)
	elif K==3: res = scipy.stats.f_oneway(x0, x1, x2)
	elif K==4: res = scipy.stats.f_oneway(x0, x1, x2, x3)
	elif K==5: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4)
	elif K==6: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5)
	elif K==7: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6)
	elif K==8: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7)
	elif K==9: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8)
	elif K==10: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9)
	elif K==11: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
	elif K==12: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11)
	elif K==13: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12)
	elif K==14: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13)
	elif K==15: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14)
	elif K==16: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15)
	elif K==17: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16)
	elif K==18: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17)
	elif K==19: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18)
	elif K==20: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19)
	elif K==21: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19, x20)
	elif K==22: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19, x20, x21)
	elif K==23: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19, x20, x21, x22)
	elif K==24: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19, x20, x21, x22, x23)
	elif K==25: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19, x20, x21, x22, x23, x24)
	elif K==26: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19, x20, x21, x22, x23, x24, x25)
	elif K==27: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19, x20, x21, x22, x23, x24, x25, x26)
	elif K==28: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19, x20, x21, x22, x23, x24, x25, x26, x27)
	elif K==29: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28)
	elif K==30: res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29)
	else:
		return None, None, None
	res_table_stat = []
	for im in range(K):
		alist = sp_statistics.table_stat(replicas[im])
		res_table_stat.append(alist)
	avgs =[]
	global_mean = 0.0
	for ir in range(K): 
		global_mean +=sum(replicas[ir])
		avgs.append(sum(replicas[ir])/float(len(replicas[ir])))
	
	summed_squared_elements = 0.0
	summed_squared_elements_within_groups = [None for i in range(K)]
	std_list =[]
	sst = 0.0
	nsamples = 0.0
	for i in range(K): nsamples +=len(replicas[i])
	for i in range(K):
		ssa_per_group = 0.0
		std_per_group = 0.0 
		for j in range(len(replicas[i])):
			sst +=replicas[i][j]**2
			summed_squared_elements +=replicas[i][j]*replicas[i][j]
			ssa_per_group +=replicas[i][j]
			std_per_group += (replicas[i][j] -avgs[i])**2
		std_list.append(numpy.sqrt(std_per_group/float(len(replicas[i]))))
		summed_squared_elements_within_groups[i] = ssa_per_group**2/float(len(replicas[i]))
		
	sst -=global_mean**2/nsamples
	ssa  = sum(summed_squared_elements_within_groups) - global_mean**2/nsamples
	sse  = sst - ssa
	n1   = 0
	for i in range(K): n1 +=len(replicas[i])-1
	msa = ssa/(K-1.0)
	mse = sse/float(n1)
	mst = sst/float(n1)
	try:    f_ratio = msa/mse
	except: f_ratio = 999999.
	log_main.add('                              ANOVA of %s'%name_of_variable)
	log_main.add('{:5} {:^12} {:^12} '.format('ANOVA', 'F-value',  'Significance'))
	log_main.add('{:5} {:12.2f} {:12.2f}'.format('ANOVA', res[0], res[1]*100.))
	log_main.add(' ')

	log_main.add('ANOVA:  %s mean of all clusters: %f'%(name_of_variable, round(global_mean/(float(nsamples)), 4)))
	log_main.add('ANOVA:  Group averages')
	log_main.add('{:5} {:^7} {:^8} {:^12} {:^12} '.format('ANOVA', 'GID', 'N',  'mean',   'std'))
	for i in range(K):
		log_main.add('{:5} {:^7d} {:^8d} {:12.4f} {:12.4f}'.format('ANOVA', i, len(replicas[i]), res_table_stat[i][0], numpy.sqrt(res_table_stat[i][1])))
	log_main.add(' ')
	log_main.add('ANOVA  Pair-wise tests')
	log_main.add('{:5} {:^3} {:^3} {:^12} {:^12} {:^12} {:^12}'.format('ANOVA', 'A', 'B', 'avgA','avgB', 'F_value', 'Significance')) 
	for ires in range(K-1):
		for jres in range(ires+1, K):
			cres = scipy.stats.f_oneway(replicas[ires], replicas[jres])
			log_main.add('{:5} {:^3d} {:^3d} {:12.4f} {:12.4f} {:12.3f} {:12.4f} '.format('ANOVA', ires, jres, avgs[ires], avgs[jres], cres[0], cres[1]*100.))
	log_main.add(' ')
	log_main.add('================================================================================================================\n')
	return res[0], res[1]
####============================== FCM related functions =================================
def assignment_to_umat(iter_assignment, number_of_groups):
	umat = np.full((len(iter_assignment), number_of_groups), 0.0, dtype = np.float64)
	for im in range(len(iter_assignment)): 
		umat[im][iter_assignment[im]] = 1.0
	return umat
	
def mix_assignment(umat, ngroups, iter_assignment, scale = 1.0, shake_rate = 0.1):
	rlist = np.random.random(len(iter_assignment))
	for im in range(len(iter_assignment)):
		tmp = np.multiply(umat[im], scale)
		if rlist[im] < shake_rate:
			tmp = np.random.random(ngroups)
		else:
			tmp[iter_assignment[im]] +=1.0
		norm = np.dot(tmp, tmp)
		umat[im] = np.multiply(tmp, 1./numpy.sqrt(norm))
	return umat
	
def compute_umat_cross(dmat, m =2.):
	(ndata, ngroups) = dmat.shape
	umat = np.zeros((ndata, ngroups))
	for k in xrange(ndata):
		scale1 = np.amin(dmat[k])
		scale2 = np.amax(dmat[k])
		for im in range(ngroups):
			u = 0.0
			tmp1 = (scale2-scale1)/(dmat[k][im] - scale1+     (scale2-scale1)/ngroups)
			for jm in range(ngroups):
				tmp2 = (scale2-scale1)/(dmat[k][jm] - scale1+ (scale2-scale1)/ngroups)
				u +=(tmp1/tmp2)**(2./(m-1))
			umat[k][im] = 1./u
	return umat

def dmat_to_umat(dmat, m =1):
	(ndat, ngrp)= dmat.shape
	umat = np.full((ndat,ngrp), 0.0, dtype = np.float64)
	for im in range(ndat):
		assi = np.argsort(np.multiply(dmat[im], -1.0))
		for jm in range(len(assi)):
			umat[im][assi[jm]]= (ngrp - jm)**m
		norm = np.dot(umat[im], umat[im])
		umat[im] = np.multiply(umat[im], 1./numpy.sqrt(norm))
	return umat

def get_PE(umat):
	(ndata, ngroups) = np.shape(umat)
	pe = 0.0
	for i in xrange(ndata):
		for j in xrange(ngroups):
			if umat[i][j] >0.0:
				pe +=umat[i][j]*umat[i][j]*numpy.log(umat[i][j])
	return pe

def do3d_fcm_groups_nofsc_smearing_iter(srdata, paramstructure, norm_per_particle,\
           umat, current_group_sizes, iteration):
	global Tracker, Blockdata
	at = time.time()
	keepgoing = 1
	if(Blockdata["myid"] == Blockdata["last_node"]):
		if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
			os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
		try:
			with open(os.path.join(Tracker["directory"],"freq_cutoff.json"),'r') as fout:
				freq_cutoff_dict = sp_utilities.convert_json_fromunicode(json.load(fout))
			fout.close()
		except: freq_cutoff_dict = 0
	else: freq_cutoff_dict = 0
	freq_cutoff_dict = sp_utilities.wrap_mpi_bcast(freq_cutoff_dict, Blockdata["last_node"], mpi.MPI_COMM_WORLD)
	rest_time  = time.time()
	
	for index_of_groups in range(Tracker["number_of_groups"]):
		tvol, tweight, trol = recons3d_trl_struct_group_nofsc_shifted_data_fcm_MPI(\
		 Blockdata["myid"], Blockdata["last_node"], srdata, paramstructure, norm_per_particle,\
		   index_of_groups, umat, None,  Tracker["constants"]["CTF"],\
		        (2*Tracker["nxinit"]+3), Tracker["nosmearing"])
			     
		if(Blockdata["myid"] == Blockdata["last_node"]):
			tvol.set_attr("is_complex",0)
			tvol.write_image(os.path.join(Tracker["directory"],    "tempdir", "tvol_2_%d.hdf"%index_of_groups))
			tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf"%index_of_groups))
			trol.write_image(os.path.join(Tracker["directory"],    "tempdir", "trol_2_%d.hdf"%index_of_groups))
			del tvol
			del tweight
			del trol
			if Tracker["do_timing"]:
				sxprint("volumes written:",(time.time()-at)/60.)
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	acc_rest = time.time() - rest_time
	if Blockdata["myid"] == Blockdata["main_node"]:
		if Tracker["do_timing"]:
			sxprint("insertion  take  %f minutes "%(acc_rest/60.))
	Tracker["maxfrad"] = Tracker["nxinit"]//2
	fsc_groups = [None for im in range(len(current_group_sizes))]
	total_size = sum(current_group_sizes)
	for igrp in range(len(current_group_sizes)):
		if Tracker["even_scale_fsc"]: 
			fsc_groups[igrp] = sp_statistics.scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], \
		     float(Tracker["constants"]["total_stack"]), total_size//len(current_group_sizes))
		else:
			fsc_groups[igrp] = sp_statistics.scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], \
		     float(Tracker["constants"]["total_stack"]), current_group_sizes[igrp])
	if Blockdata["no_of_groups"]>1:
	 	# new starts
	 	rest_time  = time.time()
		sub_main_node_list = [ -1 for i in range(Blockdata["no_of_groups"])]
		for index_of_colors in range(Blockdata["no_of_groups"]):
			for iproc in range(Blockdata["nproc"]-1):
				if Blockdata["myid"]== iproc:
					if Blockdata["color"] == index_of_colors and Blockdata["myid_on_node"] == 0:
						sub_main_node_list[index_of_colors] = Blockdata["myid"]
					sp_utilities.wrap_mpi_send(sub_main_node_list, Blockdata["last_node"], mpi.MPI_COMM_WORLD)
				if Blockdata["myid"] == Blockdata["last_node"]:
					dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
					for im in range(len(dummy)):
						if dummy[im]>-1: sub_main_node_list[im] = dummy[im]
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		sp_utilities.wrap_mpi_bcast(sub_main_node_list, Blockdata["last_node"], mpi.MPI_COMM_WORLD)
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
		for iloop in range(nbig_loop):
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
			
				if(Blockdata["myid"] == Blockdata["last_node"]):
					tvol2    = sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")%index_of_group)
					tweight2 = sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")%index_of_group)
					treg2 	 = sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_group))
					tvol2.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1})
					tag = 7007
					sp_utilities.send_EMData(tvol2,    sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
					sp_utilities.send_EMData(tweight2, sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
					sp_utilities.send_EMData(treg2,    sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
				
				elif (Blockdata["myid"] == sub_main_node_list[index_of_colors]):
					tag      = 7007
					tvol2    = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
					tweight2 = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
					treg2    = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				try: 
					Tracker["freq_fsc143_cutoff"] = freq_cutoff_dict["Cluster_%03d.txt"%index_of_group]
				except: pass	
				if Blockdata["color"] == index_of_colors:  # It has to be 1 to avoid problem with tvol1 not closed on the disk 							
					if(Blockdata["myid_on_node"] != 0): 
						tvol2     = sp_utilities.model_blank(1)
						tweight2  = sp_utilities.model_blank(1)
						treg2     = sp_utilities.model_blank(1)
					tvol2 = steptwo_mpi_reg(tvol2, tweight2, treg2,  fsc_groups[big_loop_groups[iloop][im]], \
					         Tracker["freq_fsc143_cutoff"], 0.01, True, color = index_of_colors) # has to be False!!!
					del tweight2, treg2
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				if (Blockdata["color"] == index_of_colors) and (Blockdata["myid_on_node"] == 0):
					tag = 7007
					sp_utilities.send_EMData(tvol2, Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
				elif(Blockdata["myid"] == Blockdata["last_node"]):
					tag = 7007
					tvol2    = sp_utilities.recv_EMData(sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
					tvol2.write_image(os.path.join(Tracker["directory"], "vol_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
					del tvol2
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		acc_rest = time.time() - rest_time
		if Blockdata["myid"] == Blockdata["main_node"]:
			if Tracker["do_timing"]: sxprint("rec3d steptwo  take  %f minutes "%(acc_rest/60.))
		
	else:# loop over all groups for single node
		for index_of_group in range(Tracker["number_of_groups"]):
			if(Blockdata["myid_on_node"] == 0):
				tvol2 		= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")%index_of_group)
				tweight2 	= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")%index_of_group)
				treg2 		= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_group))
				tvol2.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1})
			else:
				tvol2 		= sp_utilities.model_blank(1)
				tweight2 	= sp_utilities.model_blank(1)
				treg2		= sp_utilities.model_blank(1)
			try: 
				Tracker["freq_fsc143_cutoff"] = freq_cutoff_dict["Cluster_%03d.txt"%index_of_group]
			except: pass
			
			tvol2 = steptwo_mpi_reg(tvol2, tweight2, treg2, fsc_groups[index_of_group],\
			      Tracker["freq_fsc143_cutoff"], 0.01, True) # has to be False!!!
			del tweight2, treg2
			if(Blockdata["myid_on_node"] == 0):
				tvol2.write_image(os.path.join(Tracker["directory"], "vol_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
				del tvol2
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	keepgoing = sp_utilities.bcast_number_to_all(keepgoing, source_node = Blockdata["main_node"], mpi_comm = mpi.MPI_COMM_WORLD) # always check 
	Tracker   = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"])
	if not keepgoing: sp_global_def.ERROR("do3d_sorting_groups_trl_iter  %s"%os.path.join(Tracker["directory"], "tempdir"),"do3d_sorting_groups_trl_iter", 1, Blockdata["myid"])
	if(Blockdata["myid"] == 0) and Tracker["do_timing"]:  sxprint("Reconstructions done  ",(time.time()-at)/60.)
	return
####--------------------------------------------------------------------------------------
def recons3d_trl_struct_group_nofsc_shifted_data_fcm_MPI(myid, main_node,\
       prjlist, paramstructure, norm_per_particle, group_ID, umat, \
       mpi_comm = None, CTF = True, target_size=-1, nosmearing = False):
	"""
	  rec3d for pre-shifted data list
	  reconstructor nn4_ctfw
	"""
	if mpi_comm == None: mpi_comm = mpi.MPI_COMM_WORLD
	refvol = sp_utilities.model_blank(target_size)
	refvol.set_attr("fudge", 1.0)
	if CTF: do_ctf = 1
	else:   do_ctf = 0
	fftvol = EMAN2_cppwrap.EMData()
	weight = EMAN2_cppwrap.EMData()
	
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1,\
	    "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = EMAN2_cppwrap.Reconstructors.get( "nn4_ctfw", params)
	r.setup()
	if nosmearing:
		for im in range(len(prjlist)):
			r.insert_slice(prjlist[im], prjlist[im].get_attr( "xform.projection" ), umat[im][group_ID])
	else:
		nny = prjlist[0][0].get_ysize()
		rshifts_shrank = copy.deepcopy(Tracker["rshifts"])
		for im in range(len(rshifts_shrank)):
			rshifts_shrank[im][0] *= float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
			rshifts_shrank[im][1] *= float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
		refang =  Tracker["refang"]
		delta  =  Tracker["delta"]
		if not Tracker["constants"]["compute_on_the_fly"]: num_on_the_fly = len(prjlist)
		else: num_on_the_fly = Tracker["num_on_the_fly"][Blockdata["myid"]]
		for im in range(num_on_the_fly):
			if prjlist[im][0].get_attr("group") == group_ID:
				for jm in range(len(prjlist[im])):
					r.insert_slice(prjlist[im][jm], prjlist[im][jm].get_attr( "xform.projection" ), \
					  prjlist[im][jm].get_attr("wprob")*umat[im][group_ID])
		if (num_on_the_fly < len(prjlist)) and Tracker["constants"]["compute_on_the_fly"]:
			for im in range(num_on_the_fly, len(prjlist)):
				if prjlist[im][0].get_attr("group") == group_ID:
					avgnorm = Tracker["avgnorm"][prjlist[im][0].get_attr("chunk_id")]
					ct      = prjlist[im][0].get_attr("ctf")
					bckgn   = prjlist[im][0].get_attr("bckgnoise")
					numbor      = len(paramstructure[im][2])
					ipsiandiang = [ paramstructure[im][2][i][0]/1000  for i in range(numbor) ]
					allshifts   = [ paramstructure[im][2][i][0]%1000  for i in range(numbor) ]
					probs       = [ paramstructure[im][2][i][1] for i in range(numbor) ]
					tdir = list(set(ipsiandiang))
					data = [None]*len(Tracker["rshifts"])
					for ii in range(len(tdir)):
						#  Find the number of times given projection direction appears on the list, it is the number of different shifts associated with it.
						lshifts = sp_utilities.findall(tdir[ii], ipsiandiang)
						toprab  = 0.0
						for ki in range(len(lshifts)):toprab += probs[lshifts[ki]]
						recdata = EMAN2_cppwrap.EMData(nny,nny,1,False)
						recdata.set_attr("is_complex",0)
						for ki in range(len(lshifts)):
							lpt = allshifts[lshifts[ki]]
							if(data[lpt] == None):
								data[lpt] = sp_fundamentals.fshift(prjlist[im][0], rshifts_shrank[lpt][0], rshifts_shrank[lpt][1])
								data[lpt].set_attr("is_complex",0)
							EMAN2_cppwrap.Util.add_img(recdata, EMAN2_cppwrap.Util.mult_scalar(data[lpt], probs[lshifts[ki]]/toprab))
						recdata.set_attr_dict({"padffted":1, "is_fftpad":1,"is_fftodd":0, "is_complex_ri":1, "is_complex":1}) # preset already
						recdata = sp_filter.filt_table(recdata, bckgn)
						ipsi = tdir[ii]%100000
						iang = tdir[ii]/100000
						recdata.set_attr_dict({"bckgnoise":bckgn, "ctf":ct})
						r.insert_slice(recdata, EMAN2_cppwrap.Transform({"type":"spider", "phi":refang[iang][0], "theta":refang[iang][1], \
						  "psi":refang[iang][2]+ipsi*delta}), toprab*avgnorm/norm_per_particle[im]*umat[im][group_ID])		
	sp_utilities.reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
	sp_utilities.reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
	if myid == main_node:dummy = r.finish(True)
	mpi.mpi_barrier(mpi_comm)
	if myid == main_node: return fftvol, weight, refvol
	else:
		del fftvol
		del weight
		del refvol 
		return None, None, None
#####=====================================================================================
#####========>>> Constrained Kmeans clustering and respective utilities <<<<==============
def Kmeans_minimum_group_size_orien_groups(nbox, iter_mstep, run_iter, cdata, fdata, srdata,\
       ctf_images,paramstructure, norm_per_particle, partids, params, \
         minimum_group_size_init, freeze_changes, indep_iter, global_log, clean_volumes = True):
	global Tracker, Blockdata
	##-----------------------------------------------
	### 1. assignment np.int32; 
	### 2. particle ids np.int32; 
	### 3. EMAN2 functions python native int;
	
	#==========---------- >>>>EQKmeans initialization ==========------------ 
	#log_main	             = sp_logger.Logger(do_print=False)
	log_main                 = sp_logger.Logger(sp_logger.BaseLogger_Files(), do_print=False)
	log_main.prefix          = Tracker["directory"]+"/"
	number_of_groups         = Tracker["number_of_groups"]
	stopercnt                = Tracker["stop_mgskmeans_percentage"]
	total_iter               = 0
	partial_rec3d            = False
	changed_nptls            = 100.0
	best_score               = 100.0
	last_score               = 100.0
	best_assignment          = []
	max_iter                 = Tracker["total_number_of_iterations"]
	do_move_one_up           = 0
	move_up_counter          = 0
	changed_nptls_list       = []
	####------------->>>>>>> Determine minimum_group_size <<<<----------------------------
	if Tracker["freeze_groups"] == 0:
	
		if(Blockdata["myid"] == Blockdata["main_node"]):
			mlist, mdict = compute_nstep(Tracker["minimum_grp_size"], iter_mstep + indep_iter)
			if indep_iter ==1:
				tmp_mlist = [None for im in range(iter_mstep)]
				for im in range(0, iter_mstep):
					tmp_mlist[im] = (mlist[im] + mlist[im+1])//2
				mlist = tmp_mlist
			if minimum_group_size_init == 0: 
				minimum_group_size = mlist[0] # always start from high bound
			else: minimum_group_size = minimum_group_size_init
			log_main.add("Minimum_grp size step used in the current run:")
			log_main.add("Step            Minimum_grp size")
			global_log.add("Minimum_grp size step used in the current run:")
			global_log.add("Step            Minimum_grp size")
			for im in range(iter_mstep): # nstep controls the quality and stability of converge.
				log_main.add("%4d                %6d"%(im+1, mlist[im]))
				global_log.add("%4d                %6d"%(im+1, mlist[im]))
		else:
			mlist = 0
			mdict = 0
			minimum_group_size = 0
		mlist = sp_utilities.wrap_mpi_bcast(mlist, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		mdict = sp_utilities.wrap_mpi_bcast(mdict, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		minimum_group_size = sp_utilities.bcast_number_to_all(minimum_group_size, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		###-------------------------------------------------------------------------------
	else: # freeze small groups
		if( Blockdata["myid"] == Blockdata["main_node"]):
			mlist, mdict = compute_nstep(Tracker["minimum_grp_size"], 2)
			if minimum_group_size_init == 0: minimum_group_size = (mlist[0]+mlist[1])//2
			else: minimum_group_size =  minimum_group_size_init
		else:
			mlist = 0
			mdict = 0
			minimum_group_size = 0
		mlist = sp_utilities.wrap_mpi_bcast(mlist, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		mdict = sp_utilities.wrap_mpi_bcast(mdict, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		minimum_group_size = sp_utilities.bcast_number_to_all(minimum_group_size, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	###-----------------------------------------------------------------------------------------------------
	###-------------------------------->>>>> mask3D <<<<<<--------------------------------------------------
	if( Blockdata["myid"] == Blockdata["main_node"]):
		try:
			if os.path.exists(Tracker["constants"]["mask3D"]): # prepare mask
				mask3D = sp_utilities.get_im(Tracker["constants"]["mask3D"])
				if mask3D.get_xsize() != Tracker["nxinit"]: mask3D = sp_fundamentals.fdecimate(mask3D, \
				    Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"], True, False)
		except:
			mask3D = sp_utilities.model_circle(Tracker["constants"]["radius"], Tracker["constants"]["nnxo"], \
			    Tracker["constants"]["nnxo"], Tracker["constants"]["nnxo"])
			mask3D = sp_fundamentals.fdecimate(mask3D, Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"], True, False)
	else:
		mask3D = sp_utilities.model_blank(Tracker["nxinit"],Tracker["nxinit"],Tracker["nxinit"])
	sp_utilities.bcast_EMData_to_all(mask3D, Blockdata["myid"], Blockdata["main_node"])
	###-----------------------------------------------------------------------------------------------------
	
	### ------------->>> get ptl_ids and intial assignment <<<-------------- 
	if( Blockdata["myid"] == Blockdata["main_node"]):
		lpartids = sp_utilities.read_text_file(partids, -1)
		if len(lpartids) == 1: iter_assignment = np.random.randint(0, number_of_groups, \
		         size=len(lpartids[0]), dtype=np.int32)
		else:                  iter_assignment = np.array(lpartids[0], dtype=np.int32)
	else: iter_assignment = 0
	iter_assignment                 = sp_utilities.wrap_mpi_bcast(iter_assignment, Blockdata["main_node"]) # initial assignment
	Tracker["total_stack"]          = iter_assignment.shape[0]
	nima                            = len(cdata) # local size
	image_start, image_end          = sp_applications.MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"])
	Tracker["min_orien_group_size"] = Tracker["number_of_groups"]*Tracker["minimum_ptl_number"]
	rest_time  = time.time()
	angle_step = get_angle_step_from_number_of_orien_groups(Tracker["orientation_groups"])
	acc_rest   = time.time() - rest_time
	if Blockdata["myid"] == Blockdata["main_node"]:
		if Tracker["do_timing"]: sxprint("angle_step  take  %f seconds "%(acc_rest/1.))
	rest_time  = time.time()
	ptls_in_orien_groups = get_angle_step_and_orien_groups_mpi(params, partids, angle_step)
	acc_rest             = time.time() - rest_time
	if Blockdata["myid"] == Blockdata["main_node"]:
		if Tracker["do_timing"]: sxprint("orien search take  %f seconds "%(acc_rest/1.))
	rest_time  = time.time()
	current_mu = Tracker["mu"]
	minimum_group_size       = max(minimum_group_size, len(ptls_in_orien_groups)) # At least one particle in one orien group
	minimum_group_size_ratio = min((minimum_group_size*Tracker["number_of_groups"])/float(Tracker["total_stack"]), 0.95)	
	### -------  printed info----------------------------------------------------------------------------------------------------
	if( Blockdata["myid"] == Blockdata["main_node"]):
		log_main.add('----------------------------------------------------------------------------------------------------------------' )
		log_main.add('      ==================> MGSKmeans clustering of run_iter %d<========== '%run_iter)
		log_main.add('----------------------------------------------------------------------------------------------------------------' )
		msg = "Total_stack:  %d K = : %d  nxinit: %d  CTF:  %s  Symmetry:  %s  stop percentage: %f  3-D mask: %s focus mask: %s  Comparison method: %s  minimum_group_size: %d orien  %d"% \
		   (Tracker["total_stack"], Tracker["number_of_groups"], Tracker["nxinit"],  Tracker["constants"]["CTF"], Tracker["constants"]["symmetry"], \
		   stopercnt, Tracker["constants"]["mask3D"], Tracker["constants"]["focus3D"], Tracker["constants"]["comparison_method"], \
		       minimum_group_size, len(ptls_in_orien_groups))
		log_main.add(msg)
		log_main.add("Minimum_grp_size   low  high", minimum_group_size, Tracker["minimum_grp_size"][0],Tracker["minimum_grp_size"][1])
		log_main.add("Img_per_grp  low  high", Tracker["constants"]["img_per_grp"], Tracker["img_per_grp"][0],Tracker["img_per_grp"][1])
		global_log.add('----------------------------------------------------------------------------------------------------------------' )
		global_log.add('      ==================> MGSKmeans clustering of run_iter %d<========== '%run_iter)
		global_log.add('----------------------------------------------------------------------------------------------------------------' )
		msg = "Total_stack:  %d K = : %d  nxinit: %d  CTF:  %s  Symmetry:  %s  stop percentage: %f  3-D mask: %s focus mask: %s  Comparison method: %s  minimum_group_size: %d orien  %d"% \
		   (Tracker["total_stack"], Tracker["number_of_groups"], Tracker["nxinit"],  Tracker["constants"]["CTF"], Tracker["constants"]["symmetry"], \
		   stopercnt, Tracker["constants"]["mask3D"], Tracker["constants"]["focus3D"], Tracker["constants"]["comparison_method"], \
		       minimum_group_size, len(ptls_in_orien_groups))
		global_log.add(msg)
		global_log.add("Minimum_grp_size   low  high", minimum_group_size, Tracker["minimum_grp_size"][0],Tracker["minimum_grp_size"][1])
		global_log.add("Img_per_grp  low  high", Tracker["constants"]["img_per_grp"], Tracker["img_per_grp"][0],Tracker["img_per_grp"][1])
	###-------------------------------------------------------------------------------------------------------------------------------	
	proc_list = [[None, None] for iproc in range(Blockdata["nproc"])]
	for iproc in range(Blockdata["nproc"]):
		iproc_image_start, iproc_image_end = sp_applications.MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], iproc)
		proc_list[iproc] = [iproc_image_start, iproc_image_end]
		
	compute_noise(Tracker["nxinit"])
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	last_iter_assignment = np.copy(iter_assignment)
	best_assignment      = np.copy(iter_assignment)
	total_iter           = 0
	keepgoing            = 1
	do_partial_rec3d     = 0
	partial_rec3d        = False
	decrease_last_iter   = 0
	current_fsc_curve    = Tracker["constants"]["fsc_curve"]
	
	#-------------------------------------------------------------------------
	### share memory preparation
	#disp_unit  = np.dtype("f4").itemsize
	#refvolsize = Tracker["nxinit"]*Tracker["nxinit"]*Tracker["nxinit"]
	#emnumpy1   = EMNumPy()
	#win_vol, base_vol = mpi_win_allocate_shared(refvolsize*disp_unit, disp_unit, \
	#   MPI_INFO_NULL, Blockdata["shared_comm"])
	#if(Blockdata["myid_on_node"]!= 0): base_vol, = mpi_win_shared_query(win_vol, MPI_PROC_NULL)
	#volbuf = np.frombuffer(np.core.multiarray.int_asbuffer(base_vol, refvolsize*disp_unit), dtype = 'f4')
	#volbuf = volbuf.reshape(Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"])
	#### ---end of shared memory----------------------------------------------
	umat = assignment_to_umat(iter_assignment[image_start:image_end], Tracker['number_of_groups'])
	while total_iter < max_iter: #### Kmeans
		rest_time  = time.time()
		if(Blockdata["myid"] == Blockdata["main_node"]):
			msg = "======================================================================="
			log_main.add(msg)
			global_log.add(msg)
			msg = "Iteration %3d: particle assignment changed ratio  %f "%(total_iter, changed_nptls)
			log_main.add(msg)
			global_log.add(msg)
			sp_utilities.write_text_file(np.copy(iter_assignment).tolist(), os.path.join(Tracker["directory"],\
			     "assignment%03d.txt"%total_iter))
			if changed_nptls< 50.0: do_partial_rec3d = 1
			else:                   do_partial_rec3d = 0
		else:do_partial_rec3d = 0
		do_partial_rec3d = sp_utilities.bcast_number_to_all(do_partial_rec3d, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		update_data_assignment(cdata, srdata, iter_assignment, proc_list, Tracker["nosmearing"], Blockdata["myid"])
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)####----------
		acc_rest = time.time() - rest_time
		if Blockdata["myid"] == Blockdata["main_node"]:
			if Tracker["do_timing"]: sxprint("orien_groups take  %f minutes"%(acc_rest/60.))
		rest_time           = time.time()
		current_group_sizes = get_group_size_from_iter_assign(iter_assignment)
		
		if Tracker["use_umat"]: # fuzzy membership reconstruction
			do3d_fcm_groups_nofsc_smearing_iter(srdata, paramstructure, norm_per_particle,\
		  umat, current_group_sizes, iteration = total_iter)
		  
		else:
			do3d_sorting_groups_nofsc_smearing_iter(srdata, paramstructure, norm_per_particle,\
           partial_rec3d, current_group_sizes, iteration = total_iter)
           
		acc_rest = time.time() - rest_time
		if Blockdata["myid"] == Blockdata["main_node"]:
			if Tracker["do_timing"]: sxprint("The entire recon3d takes  %f minutes"%(acc_rest/60.))
		rest_time   = time.time()
		local_peaks = np.full((number_of_groups, nima), -1.0, dtype = np.float32)
		## compute peaks and save them in 1D list--------------------------------------------------------------
		### shared memory
		"""Multiline Comment1"""
		###----------------------------------------------------------------------------------
		if do_move_one_up == 1:
			if(Blockdata["myid"] == Blockdata["main_node"]):
				do_move_one_up  = 0
				move_up_counter+= 1
				max_iter       += 1
				log_main.add('Minimum_grp_size is updated to %d'%minimum_group_size) 
				global_log.add('Minimum_grp_size is updated to %d'%minimum_group_size) 
			else: move_up_counter = 0
			minimum_group_size = sp_utilities.bcast_number_to_all(minimum_group_size, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			do_move_one_up     = sp_utilities.bcast_number_to_all(do_move_one_up,     Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			move_up_counter    = sp_utilities.bcast_number_to_all(move_up_counter,    Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			max_iter           = sp_utilities.bcast_number_to_all(max_iter,           Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			minimum_group_size_ratio = min((minimum_group_size*Tracker["number_of_groups"])/float(Tracker["total_stack"]),0.95)		
			
		#### Preprocess reference volumes--------- Compute distance	----------------------
		dmat = [ None for im in range(number_of_groups)]	
		for iref in range(number_of_groups):
			if(Blockdata["myid"] == Blockdata["last_node"]):
				ref_vol = sp_utilities.get_im(os.path.join(Tracker["directory"],"vol_grp%03d_iter%03d.hdf"%(iref, total_iter)))
				nnn = ref_vol.get_xsize()
				if(Tracker["nxinit"] != nnn): ref_vol = sp_fundamentals.fdecimate(ref_vol, Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"], True, False)
				stat = EMAN2_cppwrap.Util.infomask(ref_vol, mask3D, False)
				ref_vol -= stat[0]
				if stat[1]!= 0.0: EMAN2_cppwrap.Util.mul_scalar(ref_vol, 1.0/stat[1])
				ref_vol *=mask3D
			else: ref_vol = sp_utilities.model_blank(Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"])
			sp_utilities.bcast_EMData_to_all(ref_vol, Blockdata["myid"], Blockdata["last_node"])
			## Image comparison	
			if Tracker["constants"]["comparison_method"] =="cross": ref_peaks = compare_two_images_cross(cdata, ref_vol, ctf_images)
			else:                                                   ref_peaks = compare_two_images_eucd(cdata, ref_vol, fdata, ctf_images)
			local_peaks[iref] = ref_peaks
			dmat[iref]        = ref_peaks
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		###------------------------------------------------------- -----------------------
		## Compute umat
		dmat = np.array(dmat).transpose()
		#fcm_umat = compute_umat_cross(dmat, 9.0)
		umat = dmat_to_umat(dmat, m = 1)# rank group assignment per particle
		###-------------------Compute dmatrix---------------------------------------------
		local_peaks = local_peaks.reshape(number_of_groups*nima)
		acc_rest    = time.time() - rest_time
		if Blockdata["myid"] == Blockdata["main_node"]:
			if Tracker["do_timing"]: sxprint("comparison and reference preparation take  %f minutes"%(acc_rest/60.))
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		rest_time  = time.time()
		tdmat      = time.time()
		# pass to main_node
		if Blockdata["myid"] == Blockdata["main_node"]:
			dmatrix = np.full((number_of_groups, Tracker["total_stack"]), 0.0, dtype=np.float32)
			for im in range(local_peaks.shape[0]): dmatrix[im//nima][im%nima + image_start] = local_peaks[im]
		else: dmatrix = 0
		if Blockdata["myid"] != Blockdata["main_node"]:
			sp_utilities.wrap_mpi_send(local_peaks, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		else:
			for iproc in range(Blockdata["nproc"]):
				if iproc != Blockdata["main_node"]:
					local_peaks = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
					iproc_nima  = proc_list[iproc][1] - proc_list[iproc][0]
					for im in range(len(local_peaks)): 
						dmatrix[im/iproc_nima][im%iproc_nima + proc_list[iproc][0]] = local_peaks[im]
		dmatrix = sp_utilities.wrap_mpi_bcast(dmatrix, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		###-------------------------------------------------------------------------------
		acc_rest = time.time() - rest_time
		if Blockdata["myid"] == Blockdata["main_node"]:
			if Tracker["do_timing"]: sxprint("compute dmatrix of various orien groups step1 take  %f minutes"%(acc_rest/60.))
		rest_time = time.time()
		###------------Ranking and Shortest distance assignment in orien groups-----------
		last_iter_assignment = np.copy(iter_assignment)
		iter_assignment      = np.full(Tracker["total_stack"], -1, dtype=np.int32)
	 	for iorien in range(len(ptls_in_orien_groups)):
			if iorien%Blockdata["nproc"] == Blockdata["myid"]:
				iter_assignment[ptls_in_orien_groups[iorien]] = \
				do_assignment_by_dmatrix_orien_group_minimum_group_size(\
				 dmatrix, ptls_in_orien_groups[iorien], Tracker["number_of_groups"], \
				   minimum_group_size_ratio)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		acc_rest = time.time() - rest_time
		if Blockdata["myid"] == Blockdata["main_node"]:
			if Tracker["do_timing"]:print("compute dmatrix of various orien groups step2 take  %f minutes"%(acc_rest/60.))
		rest_time = time.time()
		if Blockdata["myid"] != Blockdata["main_node"]: 
			sp_utilities.wrap_mpi_send(iter_assignment, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		else:
			for iproc in range(Blockdata["nproc"]):
				if iproc != Blockdata["main_node"]:
					dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
					iter_assignment[np.where(dummy>-1)] = dummy[np.where(dummy>-1)]
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		###-------------------------------------------------------------------------------
		acc_rest = time.time() - rest_time
		td = time.time() - tdmat
		if Blockdata["myid"] == Blockdata["main_node"]:
			if Tracker["do_timing"]:
				sxprint("compute dmatrix of various orien groups step3 take  %f minutes total %f minutes"%(acc_rest/60., td/60.))
		rest_time = time.time()
		###-------------------------------------->>>>>  AI make decisions <<<<<<<<--------
		if Blockdata["myid"] == Blockdata["main_node"]:
			last_score = changed_nptls
			best_score, changed_nptls, keepgoing, best_assignment, iter_assignment,\
			    do_move_one_up, minimum_group_size, freeze_changes, decrease_last_iter = \
			     AI_MGSKmeans(total_iter, iter_assignment, last_iter_assignment, best_assignment, \
			       keepgoing, best_score, stopercnt, last_score, minimum_group_size, \
			           current_mu, Tracker["constants"]["fsc_curve"], Tracker["constants"]["total_stack"],\
			               do_move_one_up, mdict, mlist, freeze_changes, Tracker["search_mode"], \
			                  Tracker["freeze_groups"], decrease_last_iter, log_main, global_log)
			changed_nptls_list.append(changed_nptls)
			tmp_list =[best_score,changed_nptls]
		else:
			iter_assignment    = 0
			best_assignment    = 0
			keepgoing          = 1
			do_move_one_up     = 0
			changed_nptls_list = 0
			tmp_list           = 0
		iter_assignment      = sp_utilities.wrap_mpi_bcast(iter_assignment,          Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		best_assignment      = sp_utilities.wrap_mpi_bcast(best_assignment,          Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		changed_nptls_list   = sp_utilities.wrap_mpi_bcast(changed_nptls_list,       Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		tmp_list             = sp_utilities.wrap_mpi_bcast(tmp_list,                 Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		keepgoing            = sp_utilities.bcast_number_to_all(keepgoing,           Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		do_move_one_up       = sp_utilities.bcast_number_to_all(do_move_one_up,      Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		freeze_changes       = sp_utilities.bcast_number_to_all(freeze_changes,      Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		decrease_last_iter   = sp_utilities.bcast_number_to_all(decrease_last_iter,  Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		####======================== Fuzzy membership =====================================
		#shake = tmp_list[-1]/100.* 0.25
		#shake = 0.1*float(minimum_group_size)/max(mlist)*tmp_list[-1]/100.
		shake = 0.0
		if nbox == 0: scale = (tmp_list[-1]/100.)**2*1.0*float(minimum_group_size)/max(mlist)
		else: 		  scale = (tmp_list[-1]/100.)**2*0.5*float(minimum_group_size)/max(mlist)
		pe1   = get_PE(umat)
		umat  = mix_assignment(umat, number_of_groups, \
		     iter_assignment[image_start:image_end], scale,  shake)
		pe2   = get_PE(umat)
		score_sub_list = np.array([pe1, pe2], dtype=np.float64)
		if Blockdata["myid"] != Blockdata["main_node"]:
			sp_utilities.wrap_mpi_send(score_sub_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		else:
			for iproc in range(Blockdata["nproc"]):
				if iproc !=Blockdata["main_node"]:
					dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
					for im in range(2):score_sub_list[im] += dummy[im] 
			score_sub_list = np.multiply(score_sub_list, 1./Tracker["total_stack"])
			log_main.add("Simple PE: %f  Mixed PE: %f "%(score_sub_list[0], score_sub_list[1]))
			if Tracker["use_umat"]:log_main.add("Use_umat is turned on.")
			else:                  log_main.add("Use_umat is turned off.")
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		total_iter +=1
		acc_rest    = time.time() - rest_time
		if Blockdata["myid"] == Blockdata["main_node"]:
			if Tracker["do_timing"]: sxprint("Compute analysis takes  %f minutes"%(acc_rest/60.))
		###-------------------------------------------------------------------------------
		if (Tracker["freeze_groups"] == 0): stop_pstd = 1.0
		else:                               stop_pstd = 2.0
		###-------------------------------------------------------------------------------
		if keepgoing == 0:
			if Blockdata["myid"] == Blockdata["main_node"]:
				log_main.add("Converge criterion reaches. Stop MGSKmeans.")
		else:
			if (total_iter >=iter_mstep + 3) and (100. not in changed_nptls_list[-3:]):
				pave, pstd, pmin, pmax = sp_statistics.table_stat(changed_nptls_list[-3:])
				if Blockdata["myid"] == Blockdata["main_node"]:
					log_main.add("Stat of changed_nptls within the last 3 iters: ")
					log_main.add("%f  %f  %f  %f "%(pave, numpy.sqrt(pstd), pmin, pmax))
				if (numpy.sqrt(pstd) < stop_pstd):
					keepgoing = 0
					if Blockdata["myid"] == Blockdata["main_node"]:
						log_main.add("No improvment in clustering exceeds limit. Stop MGSKmeans. ")
		if (keepgoing == 0):break	
	#-------------->>>>>>>> Finalize  <<<<<<<<<-------------------------------------------
	update_data_assignment(cdata, srdata, best_assignment, proc_list, Tracker["nosmearing"], Blockdata["myid"])
	partition = get_sorting_parti_from_stack(cdata)
	del iter_assignment
	del last_iter_assignment
	del best_assignment
	#mpi_win_free(win_vol)
	#emnumpy1.unregister_numpy_from_emdata()
	#del emnumpy1
	if mask3D: del mask3D
	###-----------------------------------------------------------------------------------
	if(Blockdata["myid"] == Blockdata["main_node"]):
		lpartids = sp_utilities.read_text_file(partids, -1)
		if len(lpartids) ==1: lpartids =lpartids[0]
		else: lpartids =lpartids[1]
		sp_utilities.write_text_file([partition.tolist(), lpartids], os.path.join(Tracker["directory"],"list.txt"))
		shutil.rmtree(os.path.join(Tracker["directory"], "tempdir"))
		fplist  = np.array([partition, np.array(lpartids)], dtype = np.int32)
	else: fplist = 0
	fplist  = sp_utilities.wrap_mpi_bcast(fplist, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	#######-------------------------------------------------------------------------------
	if(Blockdata["myid"] == Blockdata["last_node"]):
		if clean_volumes:# always true unless in debug mode
			for jter in range(total_iter):
				for igroup in range(Tracker["number_of_groups"]):
					os.remove(os.path.join(Tracker["directory"], "vol_grp%03d_iter%03d.hdf"%(igroup,jter)))
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	return (fplist).transpose().tolist(), total_iter, minimum_group_size, freeze_changes
###---------------------------------------------------------------------------------------
def AI_MGSKmeans(total_iters, iter_assignment, last_iter_assignment, best_assignment, \
        keepgoing, best_score, stopercnt, last_score, minimum_grp_size, current_mu, global_fsc, \
             global_total, do_move_one_up, mdict, mlist, freeze_changes, search_mode, \
               freeze_groups, decrease_last_iter, log_file, global_log):
	#------------ Single cpu function ----------------------------------------------------
	#---------->>> Empirical parameters <<<-----------
	large_cluster_scale    = 2.
	stopercnt_scale        = 2.0
	cutoff_two_rate        = 0.3
	####------->>>> Assignment analysis <<<<----------
	group_dict       = np.sort(np.unique(iter_assignment))
	number_of_groups = group_dict.shape[0]
	clusters         = []
	groups           = [] 
	for igrp in range(group_dict.shape[0]):
		one_group = np.where(iter_assignment == group_dict[igrp])
		clusters.append(one_group[0].shape[0])
		groups.append(one_group[0].tolist())
		
	### ------------------>>>  Check unbalanced groups  <<<-------------------------------
	current_img_per_grp = len(iter_assignment)//number_of_groups
	ulist               = []
	unbalanced_groups   = False
	glist, coid         = AI_split_groups(clusters)
	ratio               = coid[1]/float(coid[0])*100.
	sublist             = []
	sublist1            = []
	for im in range(len(clusters)):
		if glist[im] == 0: sublist.append(clusters[im])
		else:              sublist1.append(clusters[im])
	if len(sublist)>1: sstat = sp_statistics.table_stat(sublist)
	else:              sstat = [sublist[0], 0.0, sublist[0], sublist[0]]
	log_file.add('Avg %f  1*Std %f  Min %f  Max  %f'%(sstat[0], numpy.sqrt(sstat[1]), sstat[2], sstat[3]))
	global_log.add('Avg %f  1*Std %f  Min %f  Max  %f'%(sstat[0], numpy.sqrt(sstat[1]), sstat[2], sstat[3]))
	if len(sublist1)>1: sstat1 = sp_statistics.table_stat(sublist1)
	else:               sstat1 = [sublist1[0], 0.0, sublist1[0], sublist1[0]]
	log_file.add('Avg %f  1*Std %f  Min %f  Max  %f'%(sstat1[0], numpy.sqrt(sstat1[1]), sstat1[2], sstat1[3]))
	global_log.add('Avg %f  1*Std %f  Min %f  Max  %f'%(sstat1[0], numpy.sqrt(sstat1[1]), sstat1[2], sstat1[3]))
	###-----------------------------------------------------------------------------------
	if search_mode: # Search AI
		if len(clusters)>2:
			if sum(glist) == len(glist)-1:cutoff = coid[1]
			else:
				if ratio >33.:
					if numpy.sqrt(sstat[1])/sstat[0] > 0.50: cutoff   = int(coid[0] + numpy.sqrt(sstat[1]))  #int(2*current_img_per_grp)
					else: cutoff = int(coid[0] + 2.*numpy.sqrt(sstat[1]))
				else: cutoff = int(2.*current_img_per_grp)
		else:
			if ratio < 33.: cutoff = (coid[0]+coid[1])//2  #always balance the minor size group
			else:           cutoff = int(2.*current_img_per_grp)
		for igrp in range(group_dict.shape[0]):
			if len(groups[igrp])>cutoff:
				random.shuffle(groups[igrp])
				ulist += groups[igrp][cutoff:]
				del groups[igrp][cutoff:]
				unbalanced_groups = True
		
		if unbalanced_groups:
			random.shuffle(ulist)
			clusters = []
			for ptl in range(len(ulist)):
				igrp = ptl%number_of_groups
				groups[igrp] +=[ulist[ptl]]
			
			new_assi = np.full(len(iter_assignment), -1, dtype = np.int32)
			for igrp in range(group_dict.shape[0]):
				for im in range(len(groups[igrp])):
					new_assi[groups[igrp][im]] = group_dict[igrp]
				clusters.append(len(groups[igrp]))
			iter_assignment[:] = new_assi[:]
	####==================================================================================
	minres143 = get_res143(sp_statistics.scale_fsc_datasetsize(global_fsc, global_total, minimum_grp_size))
	ares143   = get_res143(sp_statistics.scale_fsc_datasetsize(global_fsc, global_total, \
	    len(iter_assignment)//number_of_groups))
	#####---------------------------------------------------------------------------------
	log_file.add('Group ID  group size fsc143  move up/down  size over min     min')
	log_file.add('--------->>> Current resolutions versus sizes <<<------------------')
	#log_file.add('  AVG   %8d        %3d       %3d  '%(len(iter_assignment)/number_of_groups, ares143, 0))
	log_file.add('  Min   %8d        %3d       %3d  '%(minimum_grp_size, minres143, (minres143-ares143)))
	log_file.add('------------------                            ---------------------')
	global_log.add('Group ID  group size fsc143  move up/down  size over min     min')
	global_log.add('--------->>> Current resolutions versus sizes <<<------------------')
	#global_loglog_file.add('  AVG   %8d        %3d       %3d  '%(len(iter_assignment)/number_of_groups, ares143, 0))
	global_log.add('  Min   %8d        %3d       %3d  '%(minimum_grp_size, minres143, (minres143-ares143)))
	global_log.add('------------------                            ---------------------')
	res_list = [ 0 for jl in range(number_of_groups)]
	for il in range(len(clusters)):
		cfsc = sp_statistics.scale_fsc_datasetsize(global_fsc, len(iter_assignment), clusters[il])
		res_list[il] = get_res143(cfsc)
		log_file.add('%5d   %8d        %3d       %3d     %8d    %8d'%(il, clusters[il], res_list[il],\
		    (res_list[il] - ares143), (clusters[il]-minimum_grp_size),  minimum_grp_size))
	log_file.add('------------------------------------------------------------------------')
	global_log.add('------------------------------------------------------------------------')
	ratio, newindices, stable_clusters = \
	   compare_two_iterations(iter_assignment, last_iter_assignment)
	########------------------------------------------------------------------------------
	changed_nptls = 100.- ratio*100.
	if best_score >= changed_nptls:
		best_score = changed_nptls
		best_assignment = copy.copy(iter_assignment)
	### ----------------------------------------------------------------------------------
	res_smallest  = min(res_list)
	max_res       = max(res_list)
	size_smallest = min(clusters)
	size_largest  = max(clusters)
	
	do_move_one_up =  0 # normal case
	min_index      = -1
	make_change    = False
	###-------------------------- K=2 checking -------------------------------------------
	if (number_of_groups == 2) and (float(size_smallest)/float(size_largest)<cutoff_two_rate)\
	     and (freeze_changes == 0):
		minimum_grp_size = int((len(iter_assignment)//2)*0.75) # always use large minimum_grp_size
		changed_nptls    = 100
		keepgoing        = 1
		do_move_one_up   = 1
		freeze_changes   = 1
		iter_assignment  = np.random.randint(0, number_of_groups, size = len(iter_assignment))
		log_file.add("Freeze minimum_grp_size as %d and shuffle assignment."%minimum_grp_size)	
		global_log.add("Freeze minimum_grp_size as %d and shuffle assignment."%minimum_grp_size)	
	log_file.add("freeze_changes  %d;  best score %6.2f"%(freeze_changes, best_score))
	global_log.add("freeze_changes  %d;  best score %6.2f"%(freeze_changes, best_score))
	##------------------------------------------------------------------------------------
	# 1. K is too large: Some groups have minimum_grp_size elements
	# 2. K is too small
	# Uneven clusters
	# Case 1: Normal case; decrease minimum group size; conservative mode
	####----------------------------------------------------------------------------------			
	if (freeze_groups == 0): # adjust min group size
		if (changed_nptls < 50.) and (freeze_changes == 0) and \
		       (last_score > changed_nptls) and (decrease_last_iter == 0):
			dis = float('inf')
			for im in range(len(mlist)):
				if abs(minimum_grp_size - mlist[im]) < dis:
					min_index = im
					dis = abs(minimum_grp_size - mlist[im])
			if (min_index>-1):
				if (min_index+1 <= len(mlist) - 1):
					new_size = mlist [min_index+1]
					msg ="Shrink minimum_grp_size to %d from %d."%\
					   (new_size, minimum_grp_size)
					minimum_grp_size = new_size
					log_file.add(msg)
					log_file.add("Before switching minimum_grp_size, the changed_nptls is %f"%\
					   changed_nptls)
					global_log.add(msg)
					global_log.add("Before switching minimum_grp_size, the changed_nptls is %f"%\
					   changed_nptls)
					keepgoing          = 1
					do_move_one_up     = 1
					decrease_last_iter = 1
				else:
					decrease_last_iter =1 
					freeze_changes     =1
			else: pass
		else:
			decrease_last_iter = 0
			log_file.add("Normal progress. No minimum_grp_size update.")
			global_log.add("Normal progress. No minimum_grp_size update.")
	####----------------------------------------------------------------------------------
		if (changed_nptls < stopercnt) and (minimum_grp_size == mlist[-1]): keepgoing = 0
	else:
		if (changed_nptls < stopercnt): keepgoing = 0
		
	return best_score, changed_nptls, keepgoing, best_assignment, \
	    iter_assignment, do_move_one_up, minimum_grp_size, \
	       freeze_changes, decrease_last_iter
	    
###------------------------------------------------------------------------ end of AI	
def do_assignment_by_dmatrix_orien_group_minimum_group_size(dmatrix, \
      orien_group_members, number_of_groups, minimum_group_size_ratio):
	###---------------------------------------------------------------------
	results = [[] for i in range(number_of_groups)]
	nima               = len(orien_group_members)
	minimum_group_size = int(minimum_group_size_ratio*nima/number_of_groups)
	auxmatrix          = np.full(nima, -1.0, dtype=np.float32)
	submatrix          = np.full((number_of_groups, nima), 0.0, dtype=np.float32)
	####----------------------------------------------------------------------
	for im in range(number_of_groups):
		submatrix[im] = np.multiply(auxmatrix, dmatrix[im][orien_group_members])# sort in descending order
	tmp_array = np.argsort(submatrix, axis = 1)
	####------------------------- Ranking ------------------------------------
	rmatrix   = []
	for i in range(number_of_groups): rmatrix.append(tmp_array[i].tolist())
	del tmp_array
	while len(rmatrix[0])> nima - minimum_group_size*number_of_groups:
		tarray = []
		for i in range(number_of_groups): tarray.append(rmatrix[i][0])
		value_list, index_list = np.unique(np.array(tarray, dtype = np.int32), return_index= True)
		duplicate_list = (np.setdiff1d(np.arange(number_of_groups), index_list)).tolist()
		index_list     = index_list.tolist()
		value_list     = value_list.tolist()
		if len(value_list)<  number_of_groups:
			for i in range(len(index_list)):
				if tarray[index_list[i]] ==  tarray[duplicate_list[0]]: 
					duplicate_list.append(index_list[i])# find all duplicated ones
			random.shuffle(duplicate_list)
			duplicate_list.remove(duplicate_list[0])
			for i in range(len(duplicate_list)):# swap the first row with the next row
				index_column = 1
				while rmatrix[duplicate_list[i]][index_column] in value_list: index_column +=1 # search along column non-equal ones
				value_list.append(rmatrix[duplicate_list[i]][index_column])
				rmatrix[duplicate_list[i]][0], rmatrix[duplicate_list[i]][index_column] = \
				    rmatrix[duplicate_list[i]][index_column], rmatrix[duplicate_list[i]][0]
		for i in range(number_of_groups):
			results[i].append(rmatrix[i][0])
			for j in range(number_of_groups): rmatrix[i].remove(value_list[j]) # remove K elements from each column
	kmeans_ptl_list = np.delete(np.array(range(nima), dtype = np.int32), np.array(results).ravel()) # ravel works only for even size
	del rmatrix
	####----------------------Assign images to K groups by distance--------------------------
	for iptl in range(len(kmeans_ptl_list)):
		max_indexes = np.argwhere(submatrix[:, kmeans_ptl_list[iptl]]<=submatrix[:, \
		    kmeans_ptl_list[iptl]][submatrix[:, kmeans_ptl_list[iptl]].argmin()])
		if len(max_indexes) >1:
			t = list(range(len(max_indexes)))
			random.shuffle(t)
			results[max_indexes[t[0]][0]].append(kmeans_ptl_list[iptl])
		else: results[max_indexes[0][0]].append(kmeans_ptl_list[iptl])
	iter_assignment = np.full(nima, -1, dtype=np.int32)
	for im in range(number_of_groups): iter_assignment[results[im]] = im
	del results
	del submatrix
	return iter_assignment
###=======================================================================================
### =============>>>>  various reading data <<<===========================================
### 1
def get_shrink_data_sorting(partids, partstack, return_real = False, preshift = True, apply_mask = True, npad = 1):
	# The function will read from stack a subset of images specified in partids
	#   and assign to them parameters from partstack with optional CTF application and shifting of the data.
	# So, the lengths of partids and partstack are the same.
	#  The read data is properly distributed among MPI threads.
	# 10142015 --- preshift is set to True when doing 3-D sorting.
	# chunk_id are set when data is read in
	global Tracker, Blockdata
	mask2D		= sp_utilities.model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	shrinkage 	= Tracker["nxinit"]/float(Tracker["constants"]["nnxo"])
	radius 		= int(Tracker["constants"]["radius"] * shrinkage +0.5)	
	if( Blockdata["myid"] == Blockdata["main_node"]):
		lpartids = sp_utilities.read_text_file(partids, -1)
		if len(lpartids) == 1:
			lpartids = lpartids[0]
			groupids = len(lpartids)*[-1]
		else:
			groupids = lpartids[0]
			lpartids = lpartids[1]
	else:  	
		lpartids   = 0
		groupids   = 0
	lpartids   = sp_utilities.wrap_mpi_bcast(lpartids, Blockdata["main_node"])
	groupids   = sp_utilities.wrap_mpi_bcast(groupids, Blockdata["main_node"])
	Tracker["total_stack"]  = len(lpartids)
	if(Blockdata["myid"] == Blockdata["main_node"]): partstack = sp_utilities.read_text_row(partstack)
	else:  partstack = 0
	partstack = sp_utilities.wrap_mpi_bcast(partstack, Blockdata["main_node"])
	
	if(Tracker["total_stack"] < Blockdata["nproc"]):
		sp_global_def.ERROR("Wrong MPI settings!", "get_shrink_data_sorting", 1, Blockdata["myid"])
	else: image_start, image_end = sp_applications.MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"])
	lpartids  = lpartids[image_start:image_end]
	groupids  = groupids[image_start:image_end]
	nima      =	image_end - image_start
	data      = [None]*nima
	for im in range(nima):
		data[im] = sp_utilities.get_im(Tracker["constants"]["orgstack"], lpartids[im])	
		try: phi, theta, psi, sx, sy, chunk_id, particle_group_id = \
		      partstack[lpartids[im]][0], partstack[lpartids[im]][1],\
		     partstack[lpartids[im]][2], partstack[lpartids[im]][3], partstack[lpartids[im]][4], \
		     partstack[lpartids[im]][5], partstack[lpartids[im]][6]
		except: phi, theta, psi, sx, sy, chunk_id, particle_group_id = \
		         partstack[lpartids[im]][0], partstack[lpartids[im]][1],\
		       partstack[lpartids[im]][2], partstack[lpartids[im]][3], partstack[lpartids[im]][4], \
		       partstack[lpartids[im]][5], -1		 
		if preshift:# always true
			data[im]  = sp_fundamentals.fshift(data[im],sx,sy)
			sx = 0.0
			sy = 0.0
		st = EMAN2_cppwrap.Util.infomask(data[im], mask2D, False)
		data[im] -= st[0]
		data[im] /= st[1]
		if apply_mask: data[im] = sp_morphology.cosinemask(data[im],radius = Tracker["constants"]["radius"])
		# FT
		data[im] = sp_fundamentals.fft(data[im])
		nny      = data[im].get_ysize()
		if Tracker["constants"]["CTF"]:
			ctf_params = data[im].get_attr("ctf")
			data[im]   = sp_fundamentals.fdecimate(data[im], Tracker["nxinit"]*npad, Tracker["nxinit"]*npad, 1, False, False)
			ctf_params.apix = ctf_params.apix/shrinkage
			data[im].set_attr('ctf', ctf_params)
			data[im].set_attr('ctf_applied', 0)
			if return_real :  data[im] = sp_fundamentals.fft(data[im])
		else:
			ctf_params = data[im].get_attr_default("ctf", False)
			if  ctf_params:
				ctf_params.apix = ctf_params.apix/shrinkage
				data[im].set_attr('ctf', ctf_params)
				data[im].set_attr('ctf_applied', 0)
			data[im] = sp_fundamentals.fdecimate(data[im], Tracker["nxinit"]*npad, Tracker["nxinit"]*npad, 1, True, False)
			apix     = Tracker["constants"]["pixel_size"]
			data[im].set_attr('apix', apix/shrinkage)
		if not return_real:	data[im].set_attr("padffted",1)
		t = EMAN2_cppwrap.Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
		t.set_trans(EMAN2_cppwrap.Vec2f(0.0, 0.0))
		data[im].set_attr_dict({"chunk_id":chunk_id, "group": groupids[im], "particle_group": particle_group_id,\
		    "xform.projection": t, "npad": npad})
		if Tracker["applybckgnoise"]:
			data[im].set_attr("bckgnoise", Blockdata["bckgnoise"][particle_group_id])
			data[im].set_attr("qt", float(Tracker["constants"]["nnxo"]*Tracker["constants"]["nnxo"]))
		else: data[im].set_attr("bckgnoise", Blockdata["bckgnoise"]) # constant list
	return data
#####4
def read_data_for_sorting(partids, partstack, previous_partstack):
	# The function will read from stack a subset of images specified in partids
	# and assign to them parameters from partstack with optional CTF application and shifting of the data.
	# So, the lengths of partids and partstack are the same.
	global Tracker, Blockdata
	# functions:
	# read in data
	if( Blockdata["myid"] == Blockdata["main_node"]):
		lpartids = sp_utilities.read_text_file(partids, -1)
		if len(lpartids) == 1: 
			lpartids = lpartids[0]
			groupids = [-1 for im in range(len(lpartids))]
		else:
			groupids = lpartids[0]                  
			lpartids = lpartids[1]
	else: 
		lpartids   = 0
		groupids   = 0
	lpartids = sp_utilities.wrap_mpi_bcast(lpartids, Blockdata["main_node"])
	groupids = sp_utilities.wrap_mpi_bcast(groupids, Blockdata["main_node"])
	
	Tracker["total_stack"] = len(lpartids)
	if(Blockdata["myid"] == Blockdata["main_node"]):partstack = sp_utilities.read_text_row(partstack)
	else:  partstack = 0
	partstack = sp_utilities.wrap_mpi_bcast(partstack, Blockdata["main_node"])
	if(Blockdata["myid"] == Blockdata["main_node"]): previous_partstack = sp_utilities.read_text_row(previous_partstack)
	else:  previous_partstack = 0
	previous_partstack = sp_utilities.wrap_mpi_bcast(previous_partstack, Blockdata["main_node"])
	if(Tracker["total_stack"] < Blockdata["nproc"]):
		sp_global_def.ERROR("Number of processors in use is larger than the total number of images", \
		  "get_data_and_prep", 1, Blockdata["myid"])
	else: image_start, image_end = sp_applications.MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"])
	lpartids          = lpartids[image_start:image_end]
	nima              = image_end - image_start
	data              = [ None for im in range(nima)]
	norm_per_particle = [ None for im in range(nima)]
	groupids = groupids[image_start: image_end]
	for im in range(nima):
		image = sp_utilities.get_im(Tracker["constants"]["orgstack"], lpartids[im])
		try: phi, theta, psi, sx, sy, chunk_id, particle_group_id, mnorm = partstack[lpartids[im]][0], \
		   partstack[lpartids[im]][1], partstack[lpartids[im]][2], \
			partstack[lpartids[im]][3], partstack[lpartids[im]][4], partstack[lpartids[im]][5], \
			   partstack[lpartids[im]][6], partstack[lpartids[im]][7]
		except:
			phi, theta, psi, sx, sy, chunk_id, particle_group_id, mnorm  = partstack[lpartids[im]][0], \
			    partstack[lpartids[im]][1], partstack[lpartids[im]][2], \
		  partstack[lpartids[im]][3], partstack[lpartids[im]][4], partstack[lpartids[im]][5], -1, 1.
		sx1, sy1 = previous_partstack[lpartids[im]][3], previous_partstack[lpartids[im]][4]
		t = EMAN2_cppwrap.Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
		t.set_trans(EMAN2_cppwrap.Vec2f(0.0, 0.0))
		image.set_attr_dict( {"chunk_id":chunk_id, "group":groupids[im], "particle_group":\
		    particle_group_id, "previous_shifts": [sx1, sy1], "current_shifts": [sx, sy], "xform.projection": t})
		norm_per_particle[im] = mnorm
		data[im]              = image
	return data, norm_per_particle

###6 read paramstructure
def read_paramstructure_for_sorting(partids, paramstructure_dict_file, paramstructure_dir):
	global Tracker, Blockdata
	if( Blockdata["myid"] == Blockdata["main_node"]):lcore = sp_utilities.read_text_file(partids, -1)
	else: lcore = 0
	lcore   = sp_utilities.wrap_mpi_bcast(lcore, Blockdata["main_node"])
	if len(lcore) == 1: lcore = lcore[0]
	else: lcore = lcore[1]
	psize   = len(lcore)
	oldparamstructure  = []	
	im_start, im_end   = sp_applications.MPI_start_end(psize, Blockdata["nproc"], Blockdata["myid"])
	lcore              = lcore[im_start:im_end]
	nima               = len(lcore)
	if( Blockdata["myid"] == Blockdata["main_node"]):
		tmp_list = sp_utilities.read_text_row(paramstructure_dict_file)
	else: tmp_list = 0
	tmp_list = sp_utilities.wrap_mpi_bcast(tmp_list, Blockdata["main_node"])
	pdict    = {}
	for im in range(len(lcore)):pdict[im] = tmp_list[lcore[im]]
	oldparamstructure             =  []
	nptl                          =  0
	last_old_paramstructure_file  =  None
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	for iproc in range(Blockdata["nproc"]):
		if (Blockdata["myid"] == iproc): #always read oldparamstructure sequentially
			while nptl < nima:
				[jason_of_cpu_id, chunk_id, iteration, ptl_id_on_cpu, global_index] = pdict[nptl]
				old_paramstructure_file = os.path.join(paramstructure_dir, "oldparamstructure_%d_%03d_%03d.json"%\
				  (chunk_id, jason_of_cpu_id, iteration))
				if old_paramstructure_file != last_old_paramstructure_file:
					fout = open(old_paramstructure_file,'r')
					paramstructure = sp_utilities.convert_json_fromunicode(json.load(fout))
					fout.close()
				last_old_paramstructure_file = old_paramstructure_file
				oldparamstructure.append(paramstructure[ptl_id_on_cpu])	
				nptl +=1
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	return oldparamstructure
	
###7 copy oldparamstructures from meridien regardless the previous number of cpus 
def copy_oldparamstructure_from_meridien_MPI(selected_iteration):
	global Tracker, Blockdata
	Tracker["directory"]          = os.path.join(Tracker["constants"]["masterdir"], "main%03d"%selected_iteration)
	Tracker["paramstructure_dir"] = os.path.join(Tracker["directory"], "oldparamstructure")
	old_refinement_iter_directory = os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iteration)
	old_refinement_previous_iter_directory = os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%(selected_iteration-1))
	nproc_previous = 0
	if( Blockdata["myid"] == Blockdata["main_node"]):
		if not os.path.exists(Tracker["paramstructure_dir"]):
			os.mkdir(os.path.join(Tracker["constants"]["masterdir"], "main%03d"%selected_iteration))
			os.mkdir(Tracker["paramstructure_dir"])
		Tracker["refang"]  = sp_utilities.read_text_row(os.path.join(old_refinement_iter_directory, "refang.txt"))
		sp_utilities.write_text_row(Tracker["refang"], os.path.join(Tracker["directory"], "refang.txt"))
		Tracker["rshifts"] = sp_utilities.read_text_row(os.path.join(old_refinement_iter_directory, "rshifts.txt"))
		sp_utilities.write_text_row(Tracker["refang"], os.path.join(Tracker["directory"], "rshifts.txt"))
		my_last_params = sp_utilities.read_text_file(os.path.join(old_refinement_previous_iter_directory, "params_%03d.txt"%(selected_iteration-1)), -1)
		my_parstack    = sp_utilities.read_text_file(os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt"), -1)
		my_parstack[3:5]= my_last_params[3:5]
		sp_utilities.write_text_file( my_parstack, os.path.join(Tracker["constants"]["masterdir"], "previous_refinement_parameters.txt"))
		Tracker["previous_parstack"] = os.path.join(Tracker["constants"]["masterdir"], "previous_refinement_parameters.txt")
		nproc_previous = 0
		procid         = 0
		old_refinement_iter_dir = os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iteration)
		if Blockdata["myid"] == Blockdata["main_node"]:
			while os.path.exists(os.path.join(old_refinement_iter_dir,\
		  	 "oldparamstructure","oldparamstructure_%01d_%03d_%03d.json"%(procid, \
		        nproc_previous, selected_iteration))):
				nproc_previous += 1
	else: 
		Tracker        = 0
		nproc_previous = 0
	Tracker        = sp_utilities.wrap_mpi_bcast(Tracker,             Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	nproc_previous = sp_utilities.bcast_number_to_all(nproc_previous, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	Blockdata["nproc_previous"] = nproc_previous
	oldparamstructure           =[[], []]
	local_dict                  = {}
	for procid in range(2):
		smearing_list = []
		if( Blockdata["myid"] == Blockdata["main_node"]):
			lcore = sp_utilities.read_text_file(os.path.join(Tracker["constants"]["masterdir"], \
			     "chunk_%d.txt"%procid))
		else: lcore = 0
		lcore = sp_utilities.wrap_mpi_bcast(lcore, Blockdata["main_node"], mpi.MPI_COMM_WORLD)	
		psize = len(lcore)
		oldparamstructure[procid] = []
		im_start, im_end   = sp_applications.MPI_start_end(psize, Blockdata["nproc"], Blockdata["myid"])
		local_lcore = lcore[im_start:im_end]
		istart_old_proc_id = -1
		iend_old_proc_id   = -1
		plist              = []	
		for iproc_old in range(nproc_previous):
			im_start_old, im_end_old = sp_applications.MPI_start_end(psize, Blockdata["nproc_previous"], iproc_old)
			if (im_start>= im_start_old) and im_start <=im_end_old: 
				istart_old_proc_id = iproc_old
			if (im_end  >= im_start_old) and im_end <=im_end_old: 
				iend_old_proc_id   = iproc_old
			plist.append([im_start_old, im_end_old])
		ptl_on_this_cpu = im_start
		nptl_total      = 0
		for iproc_index_old in range(istart_old_proc_id, iend_old_proc_id+1):
			with open(os.path.join(Tracker["constants"]["refinement_dir"],\
			   "main%03d"%selected_iteration, "oldparamstructure", \
			      "oldparamstructure_%01d_%03d_%03d.json"%(procid, \
			 iproc_index_old, selected_iteration)),'r') as fout:
				oldparamstructure_on_old_cpu = sp_utilities.convert_json_fromunicode(json.load(fout))
			fout.close()
			mlocal_id_on_old = ptl_on_this_cpu - plist[iproc_index_old][0]
			while (mlocal_id_on_old<len(oldparamstructure_on_old_cpu)) and (ptl_on_this_cpu<im_end):
				oldparamstructure[procid].append(oldparamstructure_on_old_cpu[mlocal_id_on_old])
				local_dict [local_lcore[nptl_total]] = [Blockdata["myid"], procid, \
				   selected_iteration, nptl_total, ptl_on_this_cpu]
				ptl_on_this_cpu  +=1
				mlocal_id_on_old +=1
				nptl_total       +=1
		del oldparamstructure_on_old_cpu
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	
		for icpu in range(Blockdata["nproc"]):#dump to disk one by one
			if  Blockdata["myid"] == icpu:
				with open(os.path.join(Tracker["constants"]["masterdir"], "main%03d"%selected_iteration, \
				   "oldparamstructure", "oldparamstructure_%01d_%03d_%03d.json"%(procid, \
					Blockdata["myid"], selected_iteration)),'w') as fout:
					json.dump(oldparamstructure[procid], fout)
				fout.close()
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	# output number of smearing
	smearing_dict = {}
	tchunk        = []
	for procid in range(2):
		if Blockdata["myid"] == Blockdata["main_node"]:
			chunk = sp_utilities.read_text_file(os.path.join(Tracker["constants"]["masterdir"], "chunk_%d.txt"%procid))
			chunk_size = len(chunk)
			smearing_list =[ None for i in range(chunk_size) ]
		else: chunk_size  = 0
		chunk_size = sp_utilities.bcast_number_to_all(chunk_size, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		local_smearing_list = []
		for im in range(len(oldparamstructure[procid])):
			local_smearing_list.append(len(oldparamstructure[procid][im][2]))
			
		if Blockdata["myid"] == Blockdata["main_node"]:
			im_start_old, im_end_old = sp_applications.MPI_start_end(chunk_size, Blockdata["nproc"], Blockdata["main_node"])
			for im in range(len(local_smearing_list)): 
				smearing_list[im_start_old+im] = local_smearing_list[im]
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		if  Blockdata["myid"] != Blockdata["main_node"]:
			sp_utilities.wrap_mpi_send(local_smearing_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		else:
			for iproc in range(Blockdata["nproc"]):
				if iproc != Blockdata["main_node"]:
					im_start_old, im_end_old = sp_applications.MPI_start_end(chunk_size, Blockdata["nproc"], iproc)
					dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
					for idum in range(len(dummy)): smearing_list[idum + im_start_old] = dummy[idum]
				else: pass
			sp_utilities.write_text_file(smearing_list, os.path.join(Tracker["constants"]["masterdir"], "smearing_%d.txt"%procid))
			for im in range(len(chunk)):
				smearing_dict[chunk[im]] =  smearing_list[im]
			tchunk +=chunk
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		
	if Blockdata["myid"] == Blockdata["main_node"]:
		tchunk.sort()
		all_smearing = [None]*len(tchunk)
		for im in range(len(tchunk)): 
			all_smearing[im] = smearing_dict[tchunk[im]]
		sp_utilities.write_text_file(all_smearing, os.path.join(Tracker["constants"]["masterdir"], "all_smearing.txt"))
		full_dict_list = [ None for im in range(Tracker["constants"]["total_stack"])]
		for key, value in list(local_dict.items()):
			full_dict_list[key] = value
		Tracker["constants"]["smearing_file"] = os.path.join(Tracker["constants"]["masterdir"], "all_smearing.txt")
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	
	for icpu in range(Blockdata["nproc"]):
		if Blockdata["myid"] == icpu and Blockdata["myid"] != Blockdata["main_node"]:
			sp_utilities.wrap_mpi_send(local_dict, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		elif Blockdata["myid"] != icpu and Blockdata["myid"] == Blockdata["main_node"]:
			local_dict = sp_utilities.wrap_mpi_recv(icpu, mpi.MPI_COMM_WORLD)
			for key, value in list(local_dict.items()):
				full_dict_list[key] = value
		else: pass
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	Tracker["paramstructure_dict"] = \
	   os.path.join(Tracker["constants"]["masterdir"], "paramstructure_dict.txt")
	if Blockdata["myid"] == Blockdata["main_node"]: 
		sp_utilities.write_text_row(full_dict_list, Tracker["paramstructure_dict"])
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	return

def get_smearing_info(nproc_previous, selected_iteration, total_stack, my_dir, refinement_dir):
	global Tracker, Blockdata
	# Do not alter any values in Tracker and Blockdata
	oldparamstructure =[[], []]
	local_dict        = {}
	for procid in range(2):
		smearing_list = []
		if( Blockdata["myid"] == Blockdata["main_node"]): lcore = sp_utilities.read_text_file(\
		   os.path.join(my_dir, "chunk_%d.txt"%procid))
		else: lcore = 0
		lcore = sp_utilities.wrap_mpi_bcast(lcore, Blockdata["main_node"], mpi.MPI_COMM_WORLD)	
		psize = len(lcore)
		oldparamstructure[procid] = []
		im_start, im_end   = sp_applications.MPI_start_end(psize, Blockdata["nproc"], Blockdata["myid"])
		local_lcore        = lcore[im_start:im_end]
		istart_old_proc_id = -1
		iend_old_proc_id   = -1
		plist              = []	
		for iproc_old in range(nproc_previous):
			im_start_old, im_end_old = sp_applications.MPI_start_end(psize, nproc_previous, iproc_old)
			if (im_start>= im_start_old) and im_start <=im_end_old: istart_old_proc_id = iproc_old
			if (im_end  >= im_start_old) and im_end <=im_end_old:   iend_old_proc_id   = iproc_old
			plist.append([im_start_old, im_end_old])
		ptl_on_this_cpu = im_start
		nptl_total      = 0
		for iproc_index_old in range(istart_old_proc_id, iend_old_proc_id+1):
			with open(os.path.join(refinement_dir,\
			   "main%03d"%selected_iteration, "oldparamstructure", \
			      "oldparamstructure_%01d_%03d_%03d.json"%(procid, \
			 iproc_index_old, selected_iteration)),'r') as fout:
				oldparamstructure_on_old_cpu = sp_utilities.convert_json_fromunicode(json.load(fout))
			fout.close()
			mlocal_id_on_old = ptl_on_this_cpu - plist[iproc_index_old][0]
			while (mlocal_id_on_old<len(oldparamstructure_on_old_cpu)) and (ptl_on_this_cpu<im_end):
				oldparamstructure[procid].append(oldparamstructure_on_old_cpu[mlocal_id_on_old])
				local_dict [local_lcore[nptl_total]] = [Blockdata["myid"], procid, \
				   selected_iteration, nptl_total, ptl_on_this_cpu]
				ptl_on_this_cpu  +=1
				mlocal_id_on_old +=1
				nptl_total       +=1
		del oldparamstructure_on_old_cpu
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		
	# output number of smearing
	smearing_dict = {}
	tchunk        = []
	for procid in range(2):
		if Blockdata["myid"] == Blockdata["main_node"]:
			chunk = sp_utilities.read_text_file(os.path.join(my_dir, "chunk_%d.txt"%procid))
			chunk_size = len(chunk)
			smearing_list =[ None for i in range(chunk_size) ]
		else: chunk_size  = 0
		chunk_size = sp_utilities.bcast_number_to_all(chunk_size, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		local_smearing_list = []
		for im in range(len(oldparamstructure[procid])):
			local_smearing_list.append(len(oldparamstructure[procid][im][2]))
		if Blockdata["myid"] == Blockdata["main_node"]:
			im_start_old, im_end_old = sp_applications.MPI_start_end(chunk_size, Blockdata["nproc"], Blockdata["main_node"])
			for im in range(len(local_smearing_list)): smearing_list[im_start_old+im] = local_smearing_list[im]
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		if  Blockdata["myid"] != Blockdata["main_node"]:
			sp_utilities.wrap_mpi_send(local_smearing_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		else:
			for iproc in range(Blockdata["nproc"]):
				if iproc != Blockdata["main_node"]:
					im_start_old, im_end_old = sp_applications.MPI_start_end(chunk_size, Blockdata["nproc"], iproc)
					dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
					for idum in range(len(dummy)): smearing_list[idum + im_start_old] = dummy[idum]
				else: pass
			sp_utilities.write_text_file(smearing_list, os.path.join(my_dir, "smearing_%d.txt"%procid))
			for im in range(len(chunk)): smearing_dict[chunk[im]] =  smearing_list[im]
			tchunk +=chunk
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	
	if Blockdata["myid"] == Blockdata["main_node"]:
		tchunk.sort()
		all_smearing = np.full(len(tchunk), 0.0, dtype=np.float32)
		for im in range(len(tchunk)): all_smearing[im] = smearing_dict[tchunk[im]]
	else: all_smearing = 0
	all_smearing = sp_utilities.wrap_mpi_bcast(all_smearing, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	return all_smearing
### 8
def precalculate_shifted_data_for_recons3D(prjlist, paramstructure, refang, rshifts, delta, avgnorms, \
       nxinit, nnxo, nosmearing, norm_per_particle = None, upweighted=False, nsmear =-1):
	nny = prjlist[0].get_ysize()
	if nosmearing:
		for im in range(len(prjlist)):
			prjlist[im].set_attr_dict({"padffted":1, "is_fftpad":1,"is_fftodd":0, "is_complex_ri":1, "is_complex":1})
			if not upweighted:prjlist[im] = sp_filter.filt_table(prjlist[im], prjlist[im].get_attr("bckgnoise"))
			prjlist[im].set_attr("wprob", 1.0)
		return prjlist
	else:
		if norm_per_particle == None: norm_per_particle = len(prjlist)*[1.0]
		recdata_list = [[] for im in range(len(prjlist))]
		rshifts_shrank = copy.deepcopy(rshifts)
		for im in range(len(rshifts_shrank)):
			rshifts_shrank[im][0] *= float(nxinit)/float(nnxo)
			rshifts_shrank[im][1] *= float(nxinit)/float(nnxo)
		if not Tracker["constants"]["compute_on_the_fly"]:num_on_the_fly = len(prjlist)
		else:
			num_on_the_fly = Tracker["num_on_the_fly"][Blockdata["myid"]] 
		for im in range(num_on_the_fly):
			bckgn       = prjlist[im].get_attr("bckgnoise")
			ct          = prjlist[im].get_attr("ctf")
			avgnorm     = Tracker["avgnorm"][prjlist[im].get_attr("chunk_id")]
			numbor      = len(paramstructure[im][2])
			ipsiandiang = [ paramstructure[im][2][i][0]/1000  for i in range(numbor) ]
			allshifts   = [ paramstructure[im][2][i][0]%1000  for i in range(numbor) ]
			probs       = [ paramstructure[im][2][i][1] for i in range(numbor) ]
			tdir = list(set(ipsiandiang))
			data = [None]*len(rshifts_shrank)
			for ii in range(len(tdir)):
				#  Find the number of times given projection direction appears on the list, it is the number of different shifts associated with it.
				lshifts = sp_utilities.findall(tdir[ii], ipsiandiang)
				toprab  = 0.0
				for ki in range(len(lshifts)):toprab += probs[lshifts[ki]]
				recdata = EMAN2_cppwrap.EMData(nny,nny, 1, False)
				recdata.set_attr("is_complex",0)
				for ki in range(len(lshifts)):
					lpt = allshifts[lshifts[ki]]
					if(data[lpt] == None):
						data[lpt] = sp_fundamentals.fshift(prjlist[im], rshifts_shrank[lpt][0], rshifts_shrank[lpt][1])
						data[lpt].set_attr("is_complex",0)
					EMAN2_cppwrap.Util.add_img(recdata, EMAN2_cppwrap.Util.mult_scalar(data[lpt], probs[lshifts[ki]]/toprab))
				recdata.set_attr_dict({"padffted":1, "is_fftpad":1,"is_fftodd":0, "is_complex_ri":1, "is_complex":1}) # preset already
				if not upweighted:recdata = sp_filter.filt_table(recdata, bckgn)
				ipsi = tdir[ii]%100000
				iang = tdir[ii]/100000
				recdata.set_attr_dict({"wprob": toprab*avgnorm/norm_per_particle[im],\
				   "xform.projection":EMAN2_cppwrap.Transform({"type":"spider", "phi":refang[iang][0],\
					 "theta":refang[iang][1], "psi":refang[iang][2]+ipsi*delta}), "bckgnoise":bckgn, "ctf":ct})
				recdata_list[im].append(recdata)
			del recdata, tdir, ipsiandiang, allshifts, probs, data
		if num_on_the_fly<len(prjlist):
			for im in range(num_on_the_fly, len(prjlist)):recdata_list[im].append(prjlist[im])
		return recdata_list
##### read data/paramstructure ends

###=====<----downsize data---->>>>
def downsize_data_for_sorting(original_data, return_real = False, preshift = True, npad = 1, norms = None):
	# The function will read from stack a subset of images specified in partids
	#   and assign to them parameters from partstack with optional CTF application and shifting of the data.
	# So, the lengths of partids and partstack are the same.
	global Tracker, Blockdata
	# functions:
	# read in data
	# apply mask, and prepare focus projection if focus3D is specified
	# return  1. cdata: data for image comparison, always in Fourier format
	#         2. rdata: data for reconstruction, 4nn return real image
	mask2D		= sp_utilities.model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	shrinkage 	= Tracker["nxinit"]/float(Tracker["constants"]["nnxo"])
	radius      = int(Tracker["constants"]["radius"] * shrinkage +0.5)
	compute_noise(Tracker["nxinit"])
	if Tracker["applybckgnoise"]:
		oneover = []
		nnx     = len(Blockdata["bckgnoise"][0])
		for i in range(len(Blockdata["bckgnoise"])):
			temp = [0.0]*nnx
			for k in range(nnx):
				if(Blockdata["bckgnoise"][i][k] > 0.0): temp[k] = 1.0/numpy.sqrt(Blockdata["bckgnoise"][i][k])
			oneover.append(temp)
		del temp
	if Tracker["constants"]["focus3D"]: # focus mask is applied
		if Blockdata["myid"] == Blockdata["main_node"]:
			focus3d    = sp_utilities.get_im(Tracker["constants"]["focus3D"])
			focus3d_nx = focus3d.get_xsize()
			if focus3d_nx != Tracker["nxinit"]: # So the decimated focus volume can be directly used
				focus3d = sp_fundamentals.resample(focus3d, float(Tracker["nxinit"])/float(focus3d_nx))
		else: focus3d = sp_utilities.model_blank(Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"])
		sp_utilities.bcast_EMData_to_all(focus3d, Blockdata["myid"], Blockdata["main_node"])
		focus3d = sp_projection.prep_vol(focus3d, 1, 1)
	#  Preprocess the data
	nima   = len(original_data)
	cdata  = [None for im in range(nima)]
	rdata  = [None for im in range(nima)]
	fdata  = [None for im in range(nima)]
	ctfs   = [None for im in range(nima)]
	ctfa = EMAN2_cppwrap.EMAN2Ctf()
	for im in range(nima):
		image = original_data[im].copy()
		chunk_id          = image.get_attr("chunk_id")
		particle_group_id = image.get_attr("particle_group")
		phi,theta,psi,s2x,s2y = sp_utilities.get_params_proj(image, xform = "xform.projection")
		[sx, sy]   = image.get_attr("previous_shifts")
		[sx1, sy1] = image.get_attr("current_shifts")
		rimage = sp_fundamentals.cyclic_shift(image, int(round(sx)), int(round(sy)))
		cimage = sp_fundamentals.fshift(image, sx1, sy1)		
		st = EMAN2_cppwrap.Util.infomask(rimage, mask2D, False)
		rimage -= st[0]
		rimage /= st[1]
		st = EMAN2_cppwrap.Util.infomask(cimage, mask2D, False)
		cimage -= st[0]
		cimage /= st[1]
		if not Tracker["nosmearing"] and norms: cimage *= Tracker["avgnorm"][chunk_id]/norms[im] #norm correction
		if Tracker["applybckgnoise"]:
			if Tracker["applymask"]:
				if Tracker["constants"]["hardmask"]: 
					cimage = sp_morphology.cosinemask(cimage, radius = Tracker["constants"]["radius"])
				else:
					bckg = sp_utilities.model_gauss_noise(1.0,Tracker["constants"]["nnxo"]+2,Tracker["constants"]["nnxo"])
					bckg.set_attr("is_complex",1)
					bckg.set_attr("is_fftpad",1)
					bckg = sp_fundamentals.fft(sp_filter.filt_table(bckg, oneover[particle_group_id]))
					#  Normalize bckg noise in real space, only region actually used.
					st = EMAN2_cppwrap.Util.infomask(bckg, mask2D, False)
					bckg -= st[0]
					bckg /= st[1]
					cimage = sp_morphology.cosinemask(cimage,radius = Tracker["constants"]["radius"], bckg = bckg)
		else:
			if Tracker["applymask"]:cimage  = sp_morphology.cosinemask(cimage, radius = Tracker["constants"]["radius"])
			else: pass
		# FT
		rimage  = sp_fundamentals.fft(rimage)
		cimage  = sp_fundamentals.fft(cimage)
		if Tracker["constants"]["CTF"] :
			ctf_params = rimage.get_attr('ctf')
			rimage      = sp_fundamentals.fdecimate(rimage, Tracker["nxinit"]*npad, Tracker["nxinit"]*npad, 1, False, False)
			cimage      = sp_fundamentals.fdecimate(cimage, Tracker["nxinit"]*npad, Tracker["nxinit"]*npad, 1, False, False)
			ctf_params.apix = ctf_params.apix/shrinkage
			rimage.set_attr('ctf', ctf_params)
			rimage.set_attr('ctf_applied', 0)
			if return_real :  rimage = sp_fundamentals.fft(rimage)
		else:
			ctf_params = rimage.get_attr_default("ctf", False)
			if ctf_params:
				ctf_params.apix = ctf_params.apix/shrinkage
				rimage.set_attr('ctf', ctf_params)
				rimage.set_attr('ctf_applied', 0)
			rimage  = sp_fundamentals.fdecimate(rimage, Tracker["nxinit"]*npad, Tracker["nxinit"]*npad, 1, True, False)
			cimage  = sp_fundamentals.fdecimate(cimage, Tracker["nxinit"]*npad, Tracker["nxinit"]*npad, 1, True, False)
			apix   = Tracker["constants"]["pixel_size"]
			rimage.set_attr('apix', apix/shrinkage)
		cimage.set_attr("padffted",1)
		cimage.set_attr("npad", npad)
		if not return_real:	
			rimage.set_attr("padffted",1)
			rimage.set_attr("npad", npad)	
		t = EMAN2_cppwrap.Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
		t.set_trans(EMAN2_cppwrap.Vec2f(0.0, 0.0))
		rimage.set_attr_dict({"particle_group":particle_group_id,"chunk_id":chunk_id, "xform.projection": t})
		cimage.set_attr_dict({"particle_group":particle_group_id,"chunk_id":chunk_id, "xform.projection": t})
		rdata[im] =  rimage
		cdata[im] =  cimage
		if Tracker["applybckgnoise"]:
			rdata[im].set_attr("bckgnoise",  Blockdata["bckgnoise"][particle_group_id])
		else:
			rdata[im].set_attr("bckgnoise",  Blockdata["bckgnoise"])
			cdata[im].set_attr("bckgnoise",  Blockdata["bckgnoise"])                    
		if Tracker["constants"]["focus3D"]:
			focusmask = sp_morphology.binarize(sp_projection.prgl(focus3d, [phi, theta, psi, 0.0, 0.0], 1, True), 1)
			cdata[im] = sp_fundamentals.fft(focusmask*sp_fundamentals.fft(cdata[im]))
			fdata[im] = focusmask
		cdata[im].set_attr("is_complex",0)
		cdata[im].set_attr("particle_group", particle_group_id)
		if Tracker["constants"]["CTF"]:
			if im ==0: ny = cdata[im].get_ysize()
			if not sp_utilities.same_ctf(rdata[im].get_attr("ctf"), ctfa):
				ctfs[im]   = sp_morphology.ctf_img_real(ny, rdata[im].get_attr("ctf"))
			else: ctfs[im] = ctfs[im-1].copy()
	return cdata, rdata, fdata, ctfs
##=====<----for 3D----->>>>
def downsize_data_for_rec3D(original_data, particle_size, return_real = False, npad = 1):
	# The function will read from stack a subset of images specified in partids
	#   and assign to them parameters from partstack with optional CTF application and shifting of the data.
	# So, the lengths of partids and partstack are the same.	
	global Tracker, Blockdata
	# functions:
	# read in data
	# apply mask, and prepare focus projection if focus3D is specified
	# return  1. cdata: data for image comparison, always in Fourier format
	#         2. rdata: data for reconstruction, 4nn return real image
	nima       = len(original_data)
	rdata      = [None for im in range(nima)]
	mask2D     = sp_utilities.model_circle(Tracker["constants"]["radius"],\
	   Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	shrinkage  = particle_size/float(Tracker["constants"]["nnxo"])
	radius     = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
	for im in range(nima):
		image = original_data[im].copy()
		chunk_id = image.get_attr("chunk_id")
		try: particle_group_id = image.get_attr("particle_group")
		except: particle_group_id = -1
		phi,theta,psi,s2x,s2y = sp_utilities.get_params_proj(image, xform = "xform.projection")
		t = EMAN2_cppwrap.Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
		t.set_trans(EMAN2_cppwrap.Vec2f(0.0, 0.0))
		[sx, sy] = image.get_attr("previous_shifts") # always for rec3D
		if Tracker["nosmearing"]: image = sp_fundamentals.fshift(image, s2x, s2y)
		else: image = sp_fundamentals.cyclic_shift(image, int(round(sx)), int(round(sy)))
		st = EMAN2_cppwrap.Util.infomask(image, mask2D, False)
		image -= st[0]
		image /= st[1]
		image  = sp_fundamentals.fft(image)
		if Tracker["constants"]["CTF"]:
			ctf_params = image.get_attr("ctf")
			image = sp_fundamentals.fdecimate(image, particle_size*npad, particle_size*npad, 1, False, False)
			ctf_params.apix = ctf_params.apix/shrinkage
			image.set_attr('ctf', ctf_params)
			image.set_attr('ctf_applied', 0)
		else:
			ctf_params = image.get_attr_default("ctf", False)
			if  ctf_params:
				ctf_params.apix = ctf_params.apix/shrinkage
				image.set_attr('ctf', ctf_params)
				image.set_attr('ctf_applied', 0)
			image = sp_fundamentals.fdecimate(image, particle_size*npad, particle_size*npad, 1, True, False)
			apix  = Tracker["constants"]["pixel_size"]
			image.set_attr('apix', apix/shrinkage)		
		if not return_real:	
			image.set_attr("padffted",1)
			image.set_attr("npad", npad)
		image.set_attr_dict({"chunk_id": chunk_id, "particle_group": particle_group_id, \
		   "xform.projection": t})
		rdata[im] =  image
		if Tracker["applybckgnoise"]:
			rdata[im].set_attr("bckgnoise", Blockdata["bckgnoise"][rdata[im].get_attr("particle_group")])
		else: rdata[im].set_attr("bckgnoise", Blockdata["bckgnoise"])
	return rdata
### end of downsize

###==========================>>> Distance comparison <<<==============================================
def compare_two_images_eucd(data, ref_vol, fdata, ctfimgs):
	global Tracker, Blockdata
	peaks = np.zeros(shape = len(data), dtype=float)
	ny      = data[0].get_ysize()
	ref_vol = sp_projection.prep_vol(ref_vol, npad = 2, interpolation_method = 1)
	#ctfs    = EMAN2Ctf()
	qt = float(Tracker["constants"]["nnxo"]*Tracker["constants"]["nnxo"])
	m = EMAN2_cppwrap.Util.unrollmask(ny, ny)
	for im in range(len(data)):
		#current_ctf = data[im].get_attr('ctf')
		#if( not same_ctf(current_ctf,ctfs) ):
		#	ctfs = ctf_img_real(ny, current_ctf)
		phi, theta, psi, s2x, s2y = sp_utilities.get_params_proj(data[im], xform = "xform.projection")
		rtemp = sp_projection.prgl(ref_vol,[phi, theta, psi, 0.0,0.0], 1, False)
		if Tracker["constants"]["focus3D"]: rtemp = sp_fundamentals.fft(sp_fundamentals.fft(rtemp)*fdata[im])
		if Tracker["applybckgnoise"]:
			peaks[im] = -EMAN2_cppwrap.Util.sqed(data[im], rtemp, ctfimgs[im], \
			    Blockdata["unrolldata"][data[im].get_attr("particle_group")])/qt
		else:
			peaks[im] = -EMAN2_cppwrap.Util.sqed(data[im], rtemp, ctfimgs[im], m)/qt
	return peaks

def compare_two_images_cross(data, ref_vol, ctfimgs):
	global Tracker, Blockdata
	ny    = data[0].get_ysize()
	m     = EMAN2_cppwrap.Util.unrollmask(ny, ny)
	peaks = np.zeros(shape = len(data), dtype=float)
	volft = sp_projection.prep_vol(ref_vol, 2, 1)
	at = time.time()
	#  Ref is in reciprocal space
	#ctfa = EMAN2Ctf()		
	for im in range(len(data)):
		phi, theta, psi, s2x, s2y = sp_utilities.get_params_proj(data[im], xform = "xform.projection")
		ref = sp_projection.prgl( volft, [phi, theta, psi, 0.0, 0.0], 1, False)
		#current_ctf = data[im].get_attr('ctf')
		#if not same_ctf(current_ctf, ctfa):
		#	ctfa = current_ctf
		#	ctfimg  = ctf_img_real(ny, ctfa)
		EMAN2_cppwrap.Util.mulclreal(ref, ctfimgs[im])
		nrmref = numpy.sqrt(EMAN2_cppwrap.Util.innerproduct(ref, ref, m))
		if Tracker["applybckgnoise"]:
			peaks[im] = EMAN2_cppwrap.Util.innerproduct(ref, data[im], \
			   Blockdata["unrolldata"][data[im].get_attr("particle_group")])/nrmref
		else:         
			peaks[im] = EMAN2_cppwrap.Util.innerproduct(ref, data[im], m)/nrmref
	if Blockdata["myid"] == Blockdata["main_node"]:
		if Tracker["do_timing"]:
			sxprint("computing distances    ",(time.time()-at)/60.)
	return peaks		
#####==========>>>  clustering assignment utilities  <<<=====================
def isin(element, test_elements, assume_unique=False, invert=False):
	element = np.asarray(element)
	return np.in1d(element, test_elements, assume_unique=assume_unique, \
		invert=invert).reshape(element.shape)

def split_partition_into_ordered_clusters_split_ucluster(partition_in, input_row_wise = True):
	# group particles  by their cluster ids; take the last one as unaccounted group
	clusters   = []
	partition = np.array(partition, dtype=np.int32)
	if input_row_wise:
		partition = partition.transpose()
	group_id = np.sort(np.unique(partition[0]))
	if group_id.shape[0] > 1: 
		for icluster in range(group_id.shape[0]):
			cluster = np.sort(
					partition[1][
						isin(
							partition[0],
							group_id[icluster]
							)
						]
					).tolist()
			if icluster < group_id.shape[0] - 1:
				clusters.append(cluster)
			else:
				ucluster = cluster
	else:
		clusters = [partition[1].tolist()]
		ucluster = []
	return clusters, ucluster

def split_partition_into_ordered_clusters(partition_in, input_is_row_wise = True):
	# partition column 0 cluster  IDs
	# partition column 1 particle IDs
	# re-index partition as consecutive groups if their cluster IDs are not
	# Extract clusters; sort clusters in size-ascendent order; transpose partitions
	if input_is_row_wise:
		partition   = (np.array(partition_in).transpose()).tolist()
	else: partition = np.array(partition_in, dtype=np.int32)
	np_cluster_ids  = np.array(partition[0], dtype=np.int32)
	np_particle_ids = np.array(partition[1], dtype=np.int32)
	group_id        = np.unique(np_cluster_ids)
	cluster_list    = []
	mask_list       = []
	N = len(partition[0])
	llist        = np.full(group_id.shape[0], 0, dtype=np.int32)
	cluster_list = [None for im in range(group_id.shape[0])]
	for ic in range(group_id.shape[0]):
		m  = isin(np_cluster_ids, group_id[ic])
		cluster_list[ic] = np_particle_ids[m].tolist()
		mask_list.append(m)
		llist[ic]    = -len(np_particle_ids[m])
	sort_indx        = np.argsort(llist)
	new_clusters_ids = np.full(N, -1, dtype=np.int32)
	new_clusters     = []
	for ic in range(group_id.shape[0]):
		new_clusters.append(cluster_list[sort_indx[ic]])
		np.place(new_clusters_ids, mask_list[sort_indx[ic]], ic)
	partition[0] = new_clusters_ids
	return new_clusters, (np.array(partition, dtype=np.int32).transpose()).tolist()

def merge_classes_into_partition_list(classes_list, do_sort=True):
	if type(classes_list) is np.ndarray: classes_list.tolist()
	if len(classes_list) == 0: return [], [[]] # Will do nothing if empty
	if len(classes_list) == 1: # rare case, however providing solution here
		parti_list = sorted(classes_list[0])
		new_index  = [[0 for im in range(len(parti_list))], parti_list]
		return parti_list, (np.array(new_index, dtype=np.int32)).transpose().tolist()
	else:# normal case
		size_list  = [ None for im in range(len(classes_list))]
		parti_list = []
		for ic  in range(len(classes_list)):
			size_list[ic] = -len(classes_list[ic])
			for im in range(len(classes_list[ic])):
				parti_list.append(classes_list[ic][im])
		parti_list = sorted(parti_list) # ptl IDs in descendant order
		# mask by group ID and replace by group ID
		parti_list = np.array(parti_list, dtype=np.int32)
		if do_sort:
			indx       = np.argsort(np.array(size_list)) # GID in ascendent order
			indx.tolist()
		else:
			indx = list(range(len(size_list)))
		cluster_ids = np.full(parti_list.shape[0], -1, dtype=np.int32)
		for im in range(len(classes_list)):
			np.place(cluster_ids, isin(parti_list, \
			np.array(classes_list[indx[im]])), np.int32(im))
		return parti_list.tolist(), (np.array([cluster_ids, \
		     parti_list],dtype=np.int32 ).transpose()).tolist()
		
def get_sorting_parti_from_stack(data):
	global Tracker, Blockdata
	#data: a subset on a CPU
	if Blockdata["myid"] == Blockdata["main_node"]:
		total_plist = np.full(Tracker["total_stack"], -1, dtype=np.int32)
	else: total_plist = 0
	for myproc in range(Blockdata["nproc"]):
		image_start,image_end = sp_applications.MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], myproc)
		plist = 0
		if Blockdata["myid"] == myproc:
			plist = np.full(len(data), -1, dtype=np.int32)
			for im in range(len(data)):
				plist[im] = data[im].get_attr("group")
		plist = sp_utilities.wrap_mpi_bcast(plist, myproc, mpi.MPI_COMM_WORLD)
		if Blockdata["myid"] == Blockdata["main_node"]:
			total_plist[image_start:image_end] = plist
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	total_plist = sp_utilities.wrap_mpi_bcast(total_plist, Blockdata["main_node"])
	return total_plist

def create_nrandom_lists(partids, number_of_groups, number_of_runs):
	# the second column denotes orignal particle IDs
	# the first column is randomized group ID 
	global Tracker, Blockdata
	
	if Blockdata["myid"] == Blockdata["main_node"]:
		random_assignment = []
		data_list = sp_utilities.read_text_file(partids, -1)
		if len(data_list)==1: sorting_data_list = data_list[0]
		else:                 sorting_data_list = data_list[1]
		group_size = len(sorting_data_list)//number_of_groups
		for index_of_random in range(number_of_runs):
			particle_dict = {}
			ll = copy.deepcopy(sorting_data_list)
			random.shuffle(ll)
			group_list = []
			for index_of_groups in range(number_of_groups):
				if index_of_groups != number_of_groups-1:
					for iparticle in ll[index_of_groups*group_size:(index_of_groups+1)*group_size]:
						particle_dict[iparticle] = index_of_groups
						group_list.append(index_of_groups)
				else:
					for iparticle in ll[index_of_groups*group_size:]:
						particle_dict[iparticle] = index_of_groups
						group_list.append(index_of_groups)
			assignment = []
			for im in range(len(sorting_data_list)):
				assignment.append([particle_dict[sorting_data_list[im]], sorting_data_list[im]])
			random_assignment.append(assignment)
			del assignment
			del ll
	else: random_assignment = 0
	random_assignment = sp_utilities.wrap_mpi_bcast(random_assignment , Blockdata["main_node"])
	return random_assignment

def create_nrandom_lists_from_given_pids(work_dir, partids, number_of_groups, number_of_runs):
	# the second column denotes orignal particle IDs
	# the first column is randomized group ID 
	global Tracker, Blockdata
	
	if Blockdata["myid"] == Blockdata["main_node"]:
		random_assignment = []
		data_list = sp_utilities.read_text_file(partids, -1)
		if len(data_list)==1: sorting_data_list= data_list[0]
		else: sorting_data_list = data_list[1]
		group_size = len(sorting_data_list)//number_of_groups
		for index_of_random in range(number_of_runs):
			particle_dict = {}
			ll = copy.deepcopy(sorting_data_list)
			random.shuffle(ll)
			group_list = []
			for index_of_groups in range(number_of_groups):
				if index_of_groups != number_of_groups -1:
					for iparticle in ll[index_of_groups*group_size:(index_of_groups+1)*group_size]:
						particle_dict[iparticle] = index_of_groups
						group_list.append(index_of_groups)
				else:
					for iparticle in ll[index_of_groups*group_size:]:
						particle_dict[iparticle] = index_of_groups
						group_list.append(index_of_groups)
			assignment = []
			for im in range(len(sorting_data_list)):
				assignment.append([particle_dict[sorting_data_list[im]], sorting_data_list[im]])
			sp_utilities.write_text_row(assignment, os.path.join(work_dir,"independent_index_%03d.txt"%index_of_random))
			random_assignment.append(assignment)
			del assignment
			del ll
	else: random_assignment = 0
	random_assignment = sp_utilities.wrap_mpi_bcast(random_assignment, Blockdata["main_node"])
	return  random_assignment
	
def assign_unaccounted_elements_mpi(glist, clusters, img_per_grp):
	# assign unaccounted images by group probabilities
	# It has better behavior in randomizing ulist than the one below.
	global Tracker, Blockdata
	icut = int(1.5*img_per_grp)
	if Blockdata["myid"]== Blockdata["main_node"]:
		for ic in range(len(clusters)):
			if len(clusters)>icut:
				random.shuffle(clusters[ic])
				glist +=clusters[ic][icut:]
				del clusters[ic][icut:]
	else:
		glist    = 0
		clusters = 0
	clusters = sp_utilities.wrap_mpi_bcast(clusters, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	glist    = sp_utilities.wrap_mpi_bcast(glist,    Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	
	if len(glist) <= Blockdata["nproc"]*30:
		if Blockdata["myid"]== Blockdata["main_node"]:
			clusters = assign_unaccounted_inverse_proportion_to_size(glist, clusters, img_per_grp)
		else: clusters = 0
		clusters = sp_utilities.wrap_mpi_bcast(clusters, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	else:
		slist = []
		clist = []
		for ic in range(len(clusters)):
			if len(clusters[ic])<=img_per_grp:
				slist.append(max(1.- float(len(clusters[ic]))/float(img_per_grp), 0.05))
			else: slist.append(0.05)
		uplist = copy.deepcopy(glist)
		image_start, image_end = sp_applications.MPI_start_end(len(uplist), Blockdata["nproc"], Blockdata["myid"])
		uplist = uplist[image_start:image_end]
		nsize  = len(uplist)//3
		for ichunk in range(3):
			if ichunk !=2: ulist = uplist[ichunk*nsize:(ichunk+1)*nsize]
			else: ulist = uplist[ichunk*nsize:]
			while len(ulist)>0:
				im =  random.randint(0, len(clusters)-1)
				random.shuffle(ulist)
				if slist[im] > random.random():
					clusters[im].append(ulist[0])
					if len(clusters[ic])<=img_per_grp:
						slist[im] = max(1.- float(len(clusters[im]))/float(img_per_grp), 0.05)
					else: slist[im] = 0.05
					del ulist[0]
					if len(ulist)== 0: break
				else: continue
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			if Blockdata["myid"]!= Blockdata["main_node"]:
				 sp_utilities.wrap_mpi_send(clusters, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			else:
				for iproc in range(Blockdata["nproc"]):
					if iproc != Blockdata["main_node"]:
						dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
						for ic in range(len(clusters)):clusters[ic]+=dummy[ic]
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			if Blockdata["myid"]== Blockdata["main_node"]:
				for ic in  range(len(clusters)): 
					clusters[ic] = list(set(clusters[ic]))
					if len(clusters[ic])<=img_per_grp:
						slist[ic] = max(1.- float(len(clusters[ic]))/float(img_per_grp), 0.05)
					else: slist[ic] = 0.05
			else:
				slist    = 0 
				clusters = 0
			slist    = sp_utilities.wrap_mpi_bcast(slist,    Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			clusters = sp_utilities.wrap_mpi_bcast(clusters, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	return clusters
"""Multiline Comment3"""
def refilling_global_scheme_mpi(clusters, unaccounted_list, number_of_clusters, log_file, swap_ratio):
	global Tracker, Blockdata
	m     = 0
	NACC  = 0
	NUACC = len(unaccounted_list)
	for cluster in clusters: NACC += len(cluster)	
	swap_ratio /=100.
	N  = NUACC + NACC
	m  = number_of_clusters - len(clusters)
	if swap_ratio > 0.0:
		if Blockdata["myid"] == Blockdata["main_node"]:
			if int(swap_ratio*NACC) > NUACC:
				unaccounted_list, clusters = swap_clusters_small_NUACC(\
						 unaccounted_list, clusters, swap_ratio)
			else:
				unaccounted_list, clusters = swap_clusters_large_NUACC(\
						 unaccounted_list, clusters, swap_ratio)
		else:
			unaccounted_list = 0
			clusters         = 0
		clusters         = sp_utilities.wrap_mpi_bcast(clusters,         Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		unaccounted_list = sp_utilities.wrap_mpi_bcast(unaccounted_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	
	avg_size = N//number_of_clusters
	m = number_of_clusters - len(clusters)
	if (Blockdata["myid"] == Blockdata["main_node"]):
		large_clusters = []
		for ic in range(len(clusters)):
			if len(clusters[ic]) > 2*avg_size:
				large_clusters.append(clusters[ic])
	else: large_clusters = 0
	large_clusters = sp_utilities.wrap_mpi_bcast(large_clusters, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	L = len(large_clusters)
	if (m == 0) and (L == 0):
		out_clusters = assign_unaccounted_elements_mpi(unaccounted_list, clusters, avg_size)
	else:
		if (m !=0): empty_clusters =[[] for ie in range(m)]
		else: empty_clusters = []
		out_clusters = fill_no_large_groups_and_unaccounted_to_m_and_rcluster_mpi(\
				 unaccounted_list, empty_clusters, clusters, NUACC, NACC)
	sorted_out_clusters = []
	for im in range(len(out_clusters)):
		sorted_out_clusters.append(sorted(out_clusters[im]))
	return sorted_out_clusters

def select_fixed_size_cluster_from_alist(ulist, img_per_grp):
	cluster = []
	random.shuffle(ulist)
	cluster += ulist[0:img_per_grp]
	del ulist[0:img_per_grp]
	return cluster, ulist
		
def swap_clusters_small_NUACC(glist, clusters, swap_ratio):
	slist     = [None for im in range(len(clusters))]
	temp_list = []
	ulist = copy.deepcopy(glist)
	for ic in range(len(clusters)):
		slist[ic] = int(swap_ratio*len(clusters[ic]))
		random.shuffle(clusters[ic])
		temp_list.extend(clusters[ic][0:slist[ic]])
		del clusters[ic][0:slist[ic]]
	ulist +=temp_list
	for ic in range(len(clusters)):
		random.shuffle(ulist)
		clusters[ic].extend(ulist[0:slist[ic]])
		del ulist[0:slist[ic]]
	return ulist, clusters
		
def swap_clusters_large_NUACC(glist, clusters, swap_ratio):
	slist     = [None for im in range(len(clusters))]
	temp_list = []
	ulist = copy.deepcopy(glist)
	for ic in range(len(clusters)):
		slist[ic] = int(swap_ratio*len(clusters[ic]))
		random.shuffle(clusters[ic])
		temp_list.extend(clusters[ic][0:slist[ic]])
		del clusters[ic][0:slist[ic]]
	for ic in range(len(clusters)):
		random.shuffle(ulist)
		clusters[ic].extend(ulist[0:slist[ic]])
		del ulist[0:slist[ic]]
	ulist +=temp_list
	random.shuffle(ulist)
	return ulist, clusters

def fill_no_large_groups_and_unaccounted_to_m_and_rcluster_mpi(\
		unaccounted_list, empty_clusters, clusters, NUACC, NACC):
	global Tracker, Blockdata
	N = NUACC + NACC
	m = len(empty_clusters)
	number_of_groups = m + len(clusters)
	avg_size = N//number_of_groups
	if m*avg_size < NUACC - avg_size:
		if Blockdata["myid"] == Blockdata["main_node"]:
			for ic in range(len(clusters)):
				random.shuffle(clusters[ic])
				if len(clusters[ic])> avg_size:
					unaccounted_list +=clusters[ic][avg_size:]
					del clusters[ic][avg_size:]
			random.shuffle(unaccounted_list)
		else: unaccounted_list = 0
		unaccounted_list = sp_utilities.wrap_mpi_bcast(unaccounted_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		clusters         = sp_utilities.wrap_mpi_bcast(clusters,         Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		tmp_clusters = []
		if m > 0:
			if Blockdata["myid"] == Blockdata["main_node"]:
				for im in range(m):
					cluster, unaccounted_list = \
					   select_fixed_size_cluster_from_alist(unaccounted_list, avg_size//2)
					tmp_clusters.append(cluster)
			else: tmp_clusters = 0
			tmp_clusters     = sp_utilities.wrap_mpi_bcast(tmp_clusters,     Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			unaccounted_list = sp_utilities.wrap_mpi_bcast(unaccounted_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD)	 
		if len(tmp_clusters)>0:
			for cluster in tmp_clusters: clusters.append(cluster)
		clusters = assign_unaccounted_elements_mpi(unaccounted_list, clusters, avg_size)
		del unaccounted_list
	else:
		for a in empty_clusters: clusters.append(a)
		clusters = assign_unaccounted_elements_mpi(unaccounted_list, clusters, avg_size)
		del unaccounted_list
	return clusters

def assign_unaccounted_inverse_proportion_to_size(glist, clusters, img_per_grp):
	# assign unaccounted images by group probabilities, single processor version
	ulist = copy.deepcopy(glist)
	number_of_groups = len(clusters)
	slist = []
	for ic in range(len(clusters)):
		if len(clusters[ic])<=img_per_grp:
			slist.append(max(1.- float(len(clusters[ic]))/float(img_per_grp), 0.05))
		else: slist.append(0.05)
	nc = 0
	while len(ulist)>0:
		im = random.randint(0, number_of_groups - 1)
		random.shuffle(ulist)
		r = random.uniform(0.0, 1.0)
		if r<slist[im]:
			clusters[im].append(ulist[0])
			if len(clusters[im])<=img_per_grp:
				slist[im] = max(1.- float(len(clusters[im])/float(img_per_grp)), 0.05)
			else: slist[im] = 0.05
			del ulist[0]
			if len(ulist)== 0: break
		else: continue
	del ulist
	return clusters

def swap_accounted_with_unaccounted_elements_mpi(accounted_file, \
       unaccounted_file, log_file, number_of_groups, swap_ratio):
	global Tracker, Blockdata
	checking_flag = 0
	if Blockdata["myid"] == Blockdata["main_node"]:
		if len(sp_utilities.read_text_row(accounted_file)) <= 1:checking_flag = 1
	checking_flag = sp_utilities.bcast_number_to_all(checking_flag, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	if checking_flag == 0:
		if Blockdata["myid"] == Blockdata["main_node"]:
			clusters, npart  = split_partition_into_ordered_clusters(sp_utilities.read_text_file(accounted_file, -1), False)
			unaccounted_list = np.array(sp_utilities.read_text_file(unaccounted_file), dtype = np.int32)
		else: 
			clusters = 0
			unaccounted_list = 0
		clusters          = sp_utilities.wrap_mpi_bcast(clusters,         Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		unaccounted_list  = sp_utilities.wrap_mpi_bcast(unaccounted_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		clusters = refilling_global_scheme_mpi(clusters, unaccounted_list.tolist(), \
		   number_of_groups, log_file, swap_ratio)	
		if Blockdata["myid"] == Blockdata["main_node"]:
			dlist, assignment_list = merge_classes_into_partition_list(clusters)
			converted_assignment_list = [ ]
			assignment_list = np.array(assignment_list).transpose()
			for jm in range(2):converted_assignment_list.append(assignment_list[jm].tolist())
		else: converted_assignment_list = 0
		converted_assignment_list = sp_utilities.wrap_mpi_bcast(converted_assignment_list, \
		    Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	else: # there exits unaccounted
		assignment_list = create_nrandom_lists(unaccounted_file, number_of_groups, 1)#MPI
		converted_assignment_list = [ ]
		assignment_list = np.array(assignment_list[0]).transpose()
		for jm in range(2):converted_assignment_list.append(assignment_list[jm].tolist())
	return converted_assignment_list
	
def patch_to_do_k_means_match_clusters_asg_new(ptp1, ptp2):
	# input are clusters
	# patch ad hoc elements to make equal number of classes for two partitions and thus two_way comparison becomes feasible
	patch_elements = []
	if len(ptp1) != len(ptp2):
		max_number = -1
		for im in range(len(ptp1)): 
			ptp1[im]   = np.array(ptp1[im], dtype = np.int32)
			max_number = max(max_number, np.max(ptp1[im]))
		for im in range(len(ptp2)):
		    ptp2[im]   = np.array(ptp2[im], dtype = np.int32)
		    max_number = max(max_number, np.max(ptp2[im]))
		
		if len(ptp1) > len(ptp2):
			ndiff = len(ptp1) - len(ptp2)
			for indiff in range(ndiff):
				l = []
				l.append(max_number+indiff+1)
				patch_elements.append(max_number+indiff+1)
				l = numpy.array(l,"int32")
				ptp2.append(l)
		else:
			ndiff = len(ptp2)-len(ptp1)
			for indiff in range(ndiff):
				l = []
				l.append(max_number+indiff+1)
				patch_elements.append(max_number + indiff + 1)
				l = numpy.array(l,"int32")
				ptp1.append(l)
	else:
		for im in range(len(ptp1)): ptp1[im] = np.array(ptp1[im],dtype = np.int32)
		for im in range(len(ptp2)): ptp2[im] = np.array(ptp2[im],dtype = np.int32)
	newindeces, list_stable, nb_tot_objs = sp_statistics.k_means_match_clusters_asg_new(ptp1, ptp2)
	new_list_stable = []
	for a in list_stable:
		a.tolist()
		if len(a)>0: new_list_stable.append(a) # remove empty ones
	return newindeces, new_list_stable, nb_tot_objs, patch_elements
	
##### ===================>>> Assignment comparison <<<=====================================
def do_boxes_two_way_comparison_mpi(nbox, input_box_parti1,\
          input_box_parti2, depth, log_main):
	global Tracker, Blockdata
	NT = 1000
	if (Blockdata["myid"] == Blockdata["main_node"]):
		stop_generation  = 0
		log_main.add('================================================================================================================')
		log_main.add(' Two-way comparison of generation %d and layer %d computed between two pairs of independent runs: %d and %d.'%(\
		      Tracker["current_generation"],Tracker["depth"],nbox, nbox+1))
		log_main.add('----------------------------------------------------------------------------------------------------------------')
		bad_clustering =  0
		ipair      =      0
		core1      = sp_utilities.read_text_file(input_box_parti1, -1)
		ptp1, tmp1 = split_partition_into_ordered_clusters(core1, False)
		core2      = sp_utilities.read_text_file(input_box_parti2, -1)
		ptp2, tmp2 = split_partition_into_ordered_clusters(core2, False)
	else:
		ptp1  = 0
		ptp2  = 0
		core1 = 0
		core2 = 0
	ptp1  = sp_utilities.wrap_mpi_bcast(ptp1,  Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	ptp2  = sp_utilities.wrap_mpi_bcast(ptp2,  Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	core1 = sp_utilities.wrap_mpi_bcast(core1, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	core2 = sp_utilities.wrap_mpi_bcast(core2, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	rand_index = compute_rand_index_mpi(core1[0], core2[0])
	if Blockdata["myid"]==Blockdata["main_node"]:
		log_main.add('Rand index of two pairs of independent runs is  %5.4f'%rand_index)
	#### before comparison we do a simulation --------------------------
	gave, gvar = do_random_groups_simulation_mpi(ptp1, ptp2)
	if Blockdata["myid"] == Blockdata["main_node"]:# sorting data processing
		#####
		msg  ='P0      '
		msg1 ='Group ID'
		length = max(len(ptp1), len(ptp2))
		for im in range(length):
			try:    msg +='{:8d} '.format(len(ptp1[im]))
			except: pass
			msg1 +='{:8d} '.format(im)
		log_main.add(msg1)	
		log_main.add(msg)
		msg ='P1      '
		for im in range(len(ptp2)): msg +='{:8d} '.format(len(ptp2[im]))
		log_main.add(msg)
		full_list  = np.sort(np.array(core1[1], dtype = np.int32))
		total_data = full_list.shape[0]
		newindeces, list_stable, nb_tot_objs, patch_elements = \
		   patch_to_do_k_means_match_clusters_asg_new(ptp1, ptp2)
		ratio_unaccounted  = 100. - nb_tot_objs/float(total_data)*100.
		current_iter_ratio = nb_tot_objs/float(total_data)*100.
		new_list   = []
		print_matching_pairs(newindeces, log_main)
		nclass     = 0
		stat_list  = []
		tmp_list   = []
		log_main.add('               Post-matching results.')
		log_main.add('{:>5} {:>8}  {:^8}   {:>15} {:>22}  {:>5}'.format(\
		  'Group', '    size',  ' status ',   'reproducibility', 'random reproducibility', ' std '))
		log_main.add("-------------------------------------------------------------")
		
		number_of_groups = Tracker["number_of_groups"]
		ares143          = get_res143(sp_statistics.scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], \
		   Tracker["constants"]["total_stack"], Tracker["total_stack"]//number_of_groups))
		minres143        = get_res143(sp_statistics.scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], \
		   Tracker["constants"]["total_stack"], (Tracker["minimum_grp_size"][0]+\
		        Tracker["minimum_grp_size"][1])//2))
		####
		largest_size  = -float("inf")
		smallest_size =  float("inf")
		
		res_list = [ 0 for jl in range(len(list_stable))   ]
		groups   = [ None for jl in range(len(list_stable))]
		
		for il in range(len(list_stable)):
			cfsc = sp_statistics.scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], \
			   Tracker["constants"]["total_stack"], len(list_stable[il]))
			res_list[il]  = get_res143(cfsc)
			largest_size  = max(largest_size,  len(list_stable[il]))
			smallest_size = min(smallest_size, len(list_stable[il]))
			groups[il]    = len(list_stable[il])
			
		dtable, coid = AI_split_groups(groups)
		binary_ratio = coid[0]/float(coid[1])
		#print('dtable', dtable, coid)
		log_main.add('------------>>> Bineary decision on groups <<<--------------------')
		for im in range(len(list_stable)):
			log_main.add('%3d    %3d '%(im, dtable[im]))
		log_main.add('-------- Two centroids of binary decision -----------')
		log_main.add('%8d   %8d   %8.3f'%(coid[0], coid[1], binary_ratio*100.))
			 
		nclass         = 0
		decision_table = [[0 for im in range(len(list_stable))] for jm in range(4)]
		decision_stat  = []
		for index_of_any in range(len(list_stable)):
			any = np.sort(list_stable[index_of_any])
			any.tolist()
			score3 = float((np.intersect1d(ptp1[newindeces[index_of_any][0]], ptp2[newindeces[index_of_any][1]])).size)\
				  /float((np.union1d(ptp1[newindeces[index_of_any][0]], ptp2[newindeces[index_of_any][1]])).size)*100.
			cres = get_res143(sp_statistics.scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], \
			   Tracker["constants"]["total_stack"], len(any)))  
			decision_table[0][index_of_any] = len(any)
			decision_table[1][index_of_any] = score3
			decision_table[2][index_of_any] = cres
			decision_table[3][index_of_any] = largest_size/float(len(any))#
		for im in range(4):
			decision_stat.append(sp_statistics.table_stat(decision_table[im]))
				####-----------------------------------------------------------------------------------------------------------		  
		for indx in range(len(list_stable)):
			any = np.sort(list_stable[indx])
			any.tolist()
			if (decision_table[1][indx] > gave[indx]+Tracker["random_group_elimination_threshold"]*gvar[indx])\
			    and (max(res_list) - res_list[indx]<=Tracker["mu"]+1) and \
			     float(largest_size)/len(any)< 3.0:# Use simple criteria for the time being.
				if ((binary_ratio > 2.0) and (dtable[indx] == 0)) or (binary_ratio <2.0):
					new_list.append(any)
					nclass +=1
					log_main.add('{:>5} {:>8d}   {:^8}      {:>7.1f}      {:>15.1f}           {:>5.1f}'.format(\
						  indx, len(any),'accepted', decision_table[1][indx], \
							gave[indx], gvar[indx]))
					stat_list.append([decision_table[1][indx],  gave[indx], gvar[indx]])
					tmp_list.append(len(any)*(-1))
				else:
					log_main.add('{:>5} {:>8d}   {:^8}      {:>7.1f}      {:>15.1f}           {:>5.1f}'.format(\
					   indx, len(any), 'rejected', decision_table[1][indx], gave[indx], gvar[indx]))
			else:
				log_main.add('{:>5} {:>8d}   {:^8}      {:>7.1f}      {:>15.1f}           {:>5.1f}'.format(\
				    indx, len(any), 'rejected', decision_table[1][indx], gave[indx], gvar[indx]))
		###---------------------------------------------------------------------------------------------------------------
		log_main.add('-------------------------------------------------------------------------')
		log_main.add('----->>>Current FSC0.143 versus group sizes<<<-----')
		for il in range(len(list_stable)):
			log_main.add('%5d   %8d        %3d       %3d'%(il,len(list_stable[il]), res_list[il], (res_list[il] - ares143)))
		log_main.add('-------------------------------------------------------------------------')
		if nclass>1: log_main.add("nclass is %d"%nclass)
			
		### ----------------------------printout -------------------------
		log_main.add("----------- Decision making table ------------- ")
		log_main.add("----------  Stat of clusters ( mean, std, min, max)  ------------")
		log_main.add("------>>>Row 0: group size; 1: reproducibility; 2: mu; 3: Max(group size)/group size  <<<------")
		for indx in range(4):
			log_main.add("%3d    %12.1f    %8.3f   %12.1f    %12.1f "%(indx, float(decision_stat[indx][0]),\
			   numpy.sqrt(decision_stat[indx][1]), float(decision_stat[indx][2]), float(decision_stat[indx][3])))
		log_main.add("-----------------------------------------------------------------\n")
		###
		if len(tmp_list)>1:
			tmp_new_list  = []
			tmp_stat_list = []
			tmp_list = np.array(tmp_list, "int32")
			tmp_list = np.argsort(tmp_list)
			for ik in range(len(tmp_list)): 
				tmp_stat_list.append(stat_list[tmp_list[ik]])
				tmp_new_list.append(new_list[tmp_list[ik]])
			new_list[:]   = tmp_new_list[:]
			stat_list[:]  = tmp_stat_list[:]
	else: nclass = 0
	nclass = sp_utilities.bcast_number_to_all(nclass,  Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	############### parse output
	if nclass == 0:
		### 
		if (Blockdata["myid"] == Blockdata["main_node"]):
			log_main.add('There are no clusters that can pass criterion checking')
			log_main.add('Sorting eliminates the smallest group, and continues.')
			Tracker["orientation_groups"]     = max(Tracker["orientation_groups"]//2, 4)
			Tracker["no_cluster_generation"]  = True
			if number_of_groups>=3: assi = np.random.randint(0, number_of_groups, size=len(full_list))
			else:                   assi = np.random.randint(0, 3, size=len(full_list))
			ptp, unaccounted_list = split_partition_into_ordered_clusters_split_ucluster([assi,\
			      np.array(full_list)], False)
			accounted_list, new_index = merge_classes_into_partition_list(ptp)
			stat_list = []
			for im in range(number_of_groups):stat_list.append([0.,  0.,  0.])
			log_main.add('================================================================================================================\n')
		else: new_index, unaccounted_list, bad_clustering, stop_generation, stat_list = 0, 0, 0, 0, 0
		new_index          = sp_utilities.wrap_mpi_bcast(new_index,               Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		unaccounted_list   = sp_utilities.wrap_mpi_bcast(unaccounted_list,        Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		stat_list          = sp_utilities.wrap_mpi_bcast(stat_list,               Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		bad_clustering     = sp_utilities.bcast_number_to_all(bad_clustering,     Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		stop_generation    = sp_utilities.bcast_number_to_all(stop_generation,    Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		Tracker            = sp_utilities.wrap_mpi_bcast(Tracker,                 Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		return new_index, unaccounted_list, bad_clustering, stop_generation, stat_list
			
	elif nclass == 1: # Force to stop this generation, and output the cluster; do not do any other box comparison
		if Blockdata["myid"]==Blockdata["main_node"]:
			stop_generation  = 1
			accounted_list, new_index = merge_classes_into_partition_list(new_list)
			utmp = np.setdiff1d(np.array(full_list, dtype=np.int32), np.array(accounted_list, dtype=np.int32))
			if utmp.shape[0]>1: unaccounted_list = np.sort(utmp)
			else:               unaccounted_list = copy.copy(utmp)
			log_main.add('Only one group found. The program will output it and stop executing the current generation.')

			box1_dir =  os.path.join(Tracker["constants"]["masterdir"], "generation_%03d"%Tracker["current_generation"], "layer%d"%Tracker["depth"], "nbox%d"%nbox)
			box2_dir =  os.path.join(Tracker["constants"]["masterdir"], "generation_%03d"%Tracker["current_generation"], "layer%d"%Tracker["depth"], "nbox%d"%(nbox+1))
			gendir   =  os.path.join(Tracker["constants"]["masterdir"], "generation_%03d"%Tracker["current_generation"])
			
			with open(os.path.join(box1_dir,"freq_cutoff.json"),'r') as fout:
				freq_cutoff_dict1 = sp_utilities.convert_json_fromunicode(json.load(fout))
			fout.close()
			
			with open(os.path.join(box2_dir,"freq_cutoff.json"),'r') as fout:
				freq_cutoff_dict2 = sp_utilities.convert_json_fromunicode(json.load(fout))
			fout.close()
			
			try:
				with open(os.path.join(gendir,"freq_cutoff.json"),'r') as fout:
					freq_cutoff_dict3 = sp_utilities.convert_json_fromunicode(json.load(fout))
				fout.close()
				freq_cutoff_dict3 = {}
			except: freq_cutoff_dict3 = {}
			ncluster = 0
			for im in range(len(newindeces)):
				try:
					f1 = freq_cutoff_dict1["Cluster_%03d.txt"%newindeces[im][0]] 
					f2 = freq_cutoff_dict2["Cluster_%03d.txt"%newindeces[im][1]]
					freq_cutoff_dict3["Cluster_%03d.txt"%ncluster] = min(f1, f2)
					ncluster +=1
				except: pass
			with open(os.path.join(gendir,"freq_cutoff.json"),'w') as fout:
				json.dump(freq_cutoff_dict3, fout)
			fout.close()
			log_main.add('================================================================================================================\n')
		else: new_index, unaccounted_list, bad_clustering, stop_generation, stat_list = 0, 0, 0, 0, 0
		unaccounted_list   = sp_utilities.wrap_mpi_bcast(unaccounted_list,        Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		new_index          = sp_utilities.wrap_mpi_bcast(new_index,               Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		stat_list          = sp_utilities.wrap_mpi_bcast(stat_list,               Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		bad_clustering     = sp_utilities.bcast_number_to_all(bad_clustering,     Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		stop_generation    = sp_utilities.bcast_number_to_all(stop_generation,    Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		return new_index, unaccounted_list.tolist(), bad_clustering, stop_generation, stat_list
	else:
		if Blockdata["myid"]==Blockdata["main_node"]:
			if not Tracker["search_mode"]: ngroups  = total_data//Tracker["constants"]["img_per_grp"]
			else:                          ngroups  = total_data//Tracker["current_img_per_grp"]
			new_list = sorted(new_list, key=len, reverse = True)
			if len(new_list)>ngroups: del new_list[-1]
			accounted_list, new_index = merge_classes_into_partition_list(new_list)
			utmp = np.setdiff1d(np.array(full_list, dtype=np.int32), np.array(accounted_list, dtype=np.int32))
			if utmp.shape[0]>1: unaccounted_list = np.sort(utmp)
			else:               unaccounted_list = copy.copy(utmp)
			log_main.add('{} {} {} {}'.format('The number of accounted for images:', \
			    len(accounted_list),'The number of unaccounted for images:', len(unaccounted_list)))
			log_main.add('The current smallest group size: %d and the largest group size: %d'%(smallest_size, largest_size))
		
			if depth <= 1: # the last layer
				box1_dir =  os.path.join(Tracker["constants"]["masterdir"], "generation_%03d"%Tracker["current_generation"], "layer%d"%Tracker["depth"], "nbox%d"%nbox)
				box2_dir =  os.path.join(Tracker["constants"]["masterdir"], "generation_%03d"%Tracker["current_generation"], "layer%d"%Tracker["depth"], "nbox%d"%(nbox+1))
				gendir   =  os.path.join(Tracker["constants"]["masterdir"], "generation_%03d"%Tracker["current_generation"])
				with open(os.path.join(box1_dir,"freq_cutoff.json"),'r') as fout:
					freq_cutoff_dict1 = sp_utilities.convert_json_fromunicode(json.load(fout))
				fout.close()
				with open(os.path.join(box2_dir,"freq_cutoff.json"),'r') as fout:
					freq_cutoff_dict2 = sp_utilities.convert_json_fromunicode(json.load(fout))
				fout.close()
				try:
					with open(os.path.join(gendir,"freq_cutoff.json"),'r') as fout:
						freq_cutoff_dict3 = sp_utilities.convert_json_fromunicode(json.load(fout))
					fout.close()
					freq_cutoff_dict3     = {}
				except: freq_cutoff_dict3 = {}
				ncluster = 0
				for im in range(len(newindeces)):
					try:
						f1 = freq_cutoff_dict1["Cluster_%03d.txt"%newindeces[im][0]] 
						f2 = freq_cutoff_dict2["Cluster_%03d.txt"%newindeces[im][1]]
						freq_cutoff_dict3["Cluster_%03d.txt"%ncluster] = min(f1, f2)
						ncluster +=1
					except: pass
				with open(os.path.join(gendir,"freq_cutoff.json"),'w') as fout:
					json.dump(freq_cutoff_dict3, fout)
				fout.close()
				log_main.add('================================================================================================================\n')
			else:
				log_main.add('================================================================================================================\n')
		else: new_index, unaccounted_list, bad_clustering, stop_generation, stat_list = 0, 0, 0, 0, 0
		unaccounted_list   = sp_utilities.wrap_mpi_bcast(unaccounted_list,        Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		new_index          = sp_utilities.wrap_mpi_bcast(new_index,               Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		stat_list          = sp_utilities.wrap_mpi_bcast(stat_list,               Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		bad_clustering     = sp_utilities.bcast_number_to_all(bad_clustering,     Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		stop_generation    = sp_utilities.bcast_number_to_all(stop_generation,    Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		return new_index, unaccounted_list.tolist(), bad_clustering, stop_generation, stat_list
		
#### ------------------- >>>   within-box comparison <<< ---------------------------------
def do_withinbox_two_way_comparison(partition_dir, nbox, nrun, niter):
	global Tracker, Blockdata
	###------------------ adjustable local parameters --------------
	iter_start_removing = 2
	cutoff_small        = 0.333
	min_size_relax_rate = 1.25
	###--------------------------------------------------------------
	log_list = []
	## For single node only
	log_list.append(' ')
	log_list.append('----------------------------------------------------------------------------------------------------------------')	      
	ipair      = 0
	core1      = sp_utilities.read_text_file(os.path.join(partition_dir, "partition_%03d.txt"%(2*ipair)),  -1)
	ptp1, tmp1 = split_partition_into_ordered_clusters(core1, False)
	core2      = sp_utilities.read_text_file(os.path.join(partition_dir, "partition_%03d.txt"%(2*ipair+1)),-1)
	ptp2, tmp2 = split_partition_into_ordered_clusters(core2, False)
	log_list.append('       Matching of sorting results of two quasi-independent runs')
	# before comparison
	group_size = [[],[]]
	msg        = 'P0      '
	msg1       ='Group ID'
	for im in range(len(ptp1)):
		msg  +='{:8d} '.format(len(ptp1[im]))
		msg1 +='{:8d} '.format(im)
		group_size[0].append(len(ptp1[im]))
	log_list.append(msg1)
	log_list.append(msg)
	msg = 'P1      '
	for im in range(len(ptp2)): 
		msg +='{:8d} '.format(len(ptp2[im]))
		group_size[1].append(len(ptp2[im]))
	log_list.append(msg)
	full_list     = np.sort(np.array(core1[1], dtype=np.int32))
	total_data    = full_list.shape[0]
	smallest_size = total_data
	largest_size  = 0
	current_number_of_groups =  Tracker["number_of_groups"]
	newindeces, list_stable, nb_tot_objs, patch_elements = \
	          patch_to_do_k_means_match_clusters_asg_new(ptp1, ptp2)
	ratio_unaccounted = 100.-nb_tot_objs/float(total_data)*100.
	current_iter_ratio = nb_tot_objs/float(total_data)*100.
	log_list.append(' ------------------------------------------------------------------------------------ ')
	log_list.append('         Two-way matching output between a MGSKmeans partition pair P0 and P1          ')
	log_list.append('M indicates that respective group of P0 sorting (row number) matches respective group of P1 sorting (column number)')
	msg ='   '
	for i in range(len(newindeces)): msg += '{:^3d}'.format(i)
	log_list.append(msg)
	for im in range(len(newindeces)):
		msg ='{:^3d}'.format(im)
		for jm in range(len(newindeces)):
			not_found = True
			for km in range(len(newindeces)):
				if newindeces[km][0] == im and newindeces[km][1] == jm:
					msg +='{:^3s}'.format(' M' )
					not_found = False
			if not_found: msg +='{:^3s}'.format('   ')
		log_list.append(msg)
	###-----------------------------------------------------------------------------------
	
	stable_clusters   = []
	selected_clusters = []
	smallest = float('inf')
	largest = -float('inf')
	
	# col 1. size;  2. reproducibility; 3. mu diff;  4.size ratio
	decision_table = [[None for im in range(len(list_stable))] for jm in range(4)]
	decision_stat  = []
	groups         = [None for im in range(len(list_stable))] 
	for im in range(len(list_stable)):
		smallest_size = min(smallest_size, len(list_stable[im]))
		largest_size  = max(largest_size,  len(list_stable[im]))
		groups[im]    = len(list_stable[im])
	maxres143 = get_res143(sp_statistics.scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"],\
	   Tracker["constants"]["total_stack"], largest_size))
	minres143 = get_res143(sp_statistics.scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"],\
	   Tracker["constants"]["total_stack"], smallest_size))
	
	dtable, coid = AI_split_groups(groups)
	ratio = coid[1]/float(coid[0])*100.
	log_list.append('------------>>> Bineary decision on groups <<<--------------------')
	for im in range(len(list_stable)): log_list.append('%3d    %3d '%(im, dtable[im]))
	log_list.append('------ Two centroids of binary decision --------------')
	log_list.append('%8d   %8d    %8.3f '%(coid[0], coid[1], ratio))
	###-----------------------------------------------------------------------------------
	score_list  = [ ]
	current_MGR = get_MGR_from_two_way_comparison(newindeces, ptp1, ptp2, total_data)
	###----------------------------------------AI-----------------------------------------
	for index_of_any in range(len(list_stable)):
		any = np.sort(list_stable[index_of_any])
		any.tolist()
		score3 = float((np.intersect1d(ptp1[newindeces[index_of_any][0]], \
		  ptp2[newindeces[index_of_any][1]])).size)\
			  /float((np.union1d(ptp1[newindeces[index_of_any][0]], \
			     ptp2[newindeces[index_of_any][1]])).size)*100.
		cres = get_res143(sp_statistics.scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], \
		  Tracker["constants"]["total_stack"], len(any)))
		decision_table[0][index_of_any] = len(any)
		decision_table[1][index_of_any] = score3
		decision_table[2][index_of_any] = cres
		decision_table[3][index_of_any] = largest_size/float(len(any))
	for im in range(4): decision_stat.append(sp_statistics.table_stat(decision_table[im]))
	### ----------------------------------------------------------------------------------
	log_list.append('--------------------------------------------------------')
	log_list.append('               Post-matching results.                   ')
	log_list.append('{:>8} {:>10} {:>17} {:>8} {:>15}'.format('Group ID', '   size', 'min random size', ' status ',   'reproducibility'))
	for indx in range(len(list_stable)):
		any = np.sort(list_stable[indx])
		any.tolist()
		if (Tracker["freeze_groups"]== 1):
			log_list.append('{:>8} {:>10d} {:>10d}        {:>8} {:>8.1f}   {:>8}  '.format(indx, \
			  len(any), current_MGR[indx],'accepted', decision_table[1][indx], decision_table[2][indx]))
			selected_clusters.append(any)
		else:
			if (nbox == 0):
				if (niter >= iter_start_removing) and (len(list_stable)> 2): # Will not disturb sorting before it is stabilized.
					if ((ratio<50. and dtable[indx] == 0) or (ratio>50.) or (sum(dtable)== (len(dtable)-1))): # size checking
					
						if (len(any) > (Tracker["minimum_grp_size"][0]*min_size_relax_rate)) and \
						   (maxres143- decision_table[2][indx] <= Tracker["mu"]+1) and \
							  (float(len(any))/largest_size >cutoff_small) and (len(any)>numpy.sqrt(decision_stat[0][1])):
							log_list.append('{:>8} {:>10d} {:>10d}        {:>8} {:>8.1f}   {:>8}  '.format(indx, \
							decision_table[0][indx], current_MGR[indx],'accepted', decision_table[1][indx], decision_table[2][indx]))
							selected_clusters.append(any)
							
						else:
							log_list.append('{:>8} {:>10d} {:>10d}        {:>8} {:>8.1f}   {:>8}   '.format(indx, \
								decision_table[0][indx], current_MGR[indx], 'rejected', score3, cres))
					else:
						log_list.append('{:>8} {:>10d} {:>10d}        {:>8} {:>8.1f}   {:>8}   '.format(indx,\
						    decision_table[0][indx], current_MGR[indx], 'rejected', score3, cres))
				else:
					log_list.append('{:>8} {:>10d} {:>10d}        {:>8} {:>8.1f}   {:>8}  '.format(indx, \
					   len(any), current_MGR[indx],'accepted', decision_table[1][indx], decision_table[2][indx]))
					selected_clusters.append(any)
			else:
				log_list.append('{:>8} {:>10d} {:>10d}        {:>8} {:>8.1f}   {:>8}  '.format(indx, \
				   len(any), current_MGR[indx],'accepted', decision_table[1][indx], decision_table[2][indx]))
				selected_clusters.append(any)
	if (len(selected_clusters) == 0): selected_clusters[:] = list_stable[:]
	###-----------------------------------------------------------------------------------	     		
	accounted_list, new_index = merge_classes_into_partition_list(selected_clusters)
	utmp = np.setdiff1d(full_list, np.array(accounted_list, dtype=np.int32))
	if utmp.shape[0] >1: unaccounted_list = (np.sort(utmp)).tolist()
	else:                unaccounted_list = utmp.tolist()
	sp_utilities.write_text_row(new_index,         os.path.join(partition_dir, "Accounted.txt"))
	sp_utilities.write_text_file(unaccounted_list, os.path.join(partition_dir, "Core_set.txt"))
	##------------------------------------------------------------------------------------
	log_list.append('The overall reproducibility is %5.1f%%.'%current_iter_ratio)
	log_list.append('The number of accounted for images: %d.  The number of core set images: %d.'%(len(accounted_list), len(unaccounted_list)))
	log_list.append('The smallest group size: %d the largest group size: %d.'%(smallest_size, largest_size))
	log_list.append('-------------------------------------------------------------------------')
	### printout
	log_list.append("------  Stat of clusters (mean, std, min, max) ------")
	log_list.append("------>>>Row 0: group size; 1: reproducibility; 2: mu; 3: Max(group size)/group size 4.Kmeans_size <<<------")
	for indx in range(4):
		log_list.append("%3d    %12.1f    %8.3f   %12.1f    %12.1f "%(indx, float(decision_stat[indx][0]),\
		    numpy.sqrt(decision_stat[indx][1]), float(decision_stat[indx][2]), float(decision_stat[indx][3])))
	log_list.append("----------------------------------------------------------\n")
	### ----------------------------------------------------------------------------------
	b1, b2 = AI_within_box(decision_table, cutoff_ratio = 4.0, \
	     cutoff_size =  (Tracker["minimum_grp_size"][0]+Tracker["minimum_grp_size"][0])//2,\
	       cutoff_mu = Tracker["mu"]+1)
	
	log_list.append("-------- Decision making table ---------- ")
	for im in range(len(list_stable)):
		if b1[im]: c1 = 1
		else:      c1 = 0
		if b2[im]: c2 = 1
		else:      c2 = 0
		log_list.append(" %3d   %8d  c1  %1d  c2  %1d"%(im, len(list_stable[im]), c1, c2))
	log_list.append("-------- Cumulative statistics on group sizes ---------- ")
	
	nlist  = [ None for im in range(len(list_stable)) ]
	for im in range(len(list_stable)): nlist[im] = len(list_stable[im])
		
	nlist = sorted(nlist, reverse = True)
	log_list.append("%3d    %12.3f    %12.3f   %12.3f    %12.3f "%(1,\
		  float(nlist[0]), 0.0, float(nlist[0]), float(nlist[0])))
		  
	for im in range(1, len(list_stable)): 
		astat = sp_statistics.table_stat(nlist[0:im+1])
		log_list.append("%3d    %12.3f    %12.3f   %12.3f    %12.3f "%(im+1,\
		  astat[0], numpy.sqrt(astat[1]), astat[2], astat[3]))
		  
	log_list.append("----------------------------------------------------------\n")
	return smallest_size, largest_size, selected_clusters, unaccounted_list, \
	      current_iter_ratio, current_number_of_groups, log_list	      
####------Rand index----------------------------------------------------------------------	
def compute_rand_index_mpi(inassign_in1, inassign_in2):
	# Two partitions always have the same length
	# rand_index = (a + b)/[(n-1)*n/2] ; n is the full length of both partitions
	# a: pairs in the same class in both partitions
	# b: pairs not in the same class in both partitions
	# Journal of the American Statistical Association.
	# American Statistical Association. 66(336): 846-850
	# rand_index should be between [0, 1]; 
	# =0; two partitions different; =1; two partitions identical
	global Tracker, Blockdata
	assign1 = np.array(inassign_in1, dtype=np.int8)
	assign2 = np.array(inassign_in2, dtype=np.int8)
	num_in_both     = 0
	num_in_neither  = 0
	ntot            = len(assign1)
	num_all_pairs   = (ntot-1)/2.
	parti_list      = []
	for im in range(ntot-1):
		if im%Blockdata["nproc"] == Blockdata["myid"]:
			parti_list.append(im)
	for lm in range(len(parti_list)):
		for im in range(parti_list[lm], parti_list[lm]+1):
			for jm in range(im +1, ntot):
				if (assign1[im] == assign1[jm]) and (assign2[im] == assign2[jm]): num_in_both    +=1
				if (assign1[im] != assign1[jm]) and (assign2[im] != assign2[jm]): num_in_neither +=1
	nomin = (num_in_both+num_in_neither)/float(ntot)
	del assign1, assign2
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	if (Blockdata["myid"] != Blockdata["main_node"]):
		sp_utilities.wrap_mpi_send(nomin, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	else:
		for iproc in range(Blockdata["nproc"]):
			if (iproc !=Blockdata["main_node"]):
				nomin += sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
	nomin = sp_utilities.bcast_number_to_all(nomin,  Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	return float(nomin)/num_all_pairs
### reproducibility 	
def get_MGR_from_two_way_comparison(newindeces, clusters1, clusters2, N):
	rnd_grp_sizes = {}
	K = len(newindeces)
	reordered_cluster2 =[ None for i in range(len(clusters2))]
	for ij in range(len(newindeces)):
		reordered_cluster2[newindeces[ij][0]] = clusters2[newindeces[ij][1]]
	table_k_k =[[ None for i in range(K)] for j in range(K)]
	for i in range(K):
		for j in range(K):
			if not (clusters1[i] is None) and not (reordered_cluster2[j] is None):
				table_k_k[i][j] = len(set(clusters1[i]).intersection(set(reordered_cluster2[j])))
			else: table_k_k[i][j] = 0
	sum_rows = [ 0 for i in range(K)]
	sum_cols = [ 0 for i in range(K)]
	for i in range(K):
		for j in range(K):
			sum_rows[i] +=table_k_k[i][j]
			sum_cols[j] +=table_k_k[i][j]
	diagonal_k = [None for i in range(K)]
	for i in range(K): diagonal_k[i] = (sum_rows[i] + sum_cols[i])//2 # nk
	min_sizes = [None for i in range(K)]
	for i in range(K): min_sizes[i] = diagonal_k[i]**2/N
	return min_sizes
### -------->>>  between boxes  <<<-------------------------
def do_random_groups_simulation_mpi(ptp1, ptp2):
	global Tracker, Blockdata
	# return two lists: group avgs and group stds. The last one of two lists are the total avg and std.
	if (len(ptp1)>= 50) or (len(ptp2)>= 50):
		if(Blockdata["myid"] == Blockdata["main_node"]):
			sxprint('Warning: too many simulaton groups')
	Nloop = max(1000//Blockdata["nproc"], 5)
	NT    = 1000
	a     = []
	b     = []
	for ic in range(len(ptp1)): a+=ptp1[ic]
	for ic in range(len(ptp2)): b+=ptp2[ic]
	tsize = float(len(set(a+b)))
	
	nsize1   = 0
	plist1   = []
	for i1 in range(len(ptp1)):
		plist1.append([nsize1, nsize1+ max(int(float(len(ptp1[i1]))/tsize*100.), 1)])
		nsize1 += max(int(float(len(ptp1[i1]))/tsize*100.), 1)
	nsize2   = 0
	plist2   = []
	for i1 in range(len(ptp2)):
		plist2.append([nsize2, nsize2 + max(int(float(len(ptp2[i1]))/tsize*100.), 1)])
		nsize2 += max(int(float(len(ptp2[i1]))/tsize*100.), 1)
		
	if len(ptp1)>len(ptp2):
		for j in range(len(ptp1)-len(ptp2)):
			plist2.append([nsize2, nsize2+1])
			nsize2 += 1
			
	elif len(ptp2)>len(ptp1):
		for j in range(len(ptp2)-len(ptp1)):
			plist1.append([nsize1, nsize1+1])
			nsize1 += 1
	if (len(ptp1)>=50) or (len(ptp2)>=50):# will under-estimate random reproducibility
		alist = list(range(3*max(len(ptp1), len(ptp2))))
		blist = list(range(3*max(len(ptp1), len(ptp2))))
	else:
		alist = list(range(100))
		blist = list(range(100))
	k     = len(plist1)
	gave  = [ 0.0 for i in range(k)]
	gvar  = [ 0.0 for i in range(k)]
	svar  = 0.0
	save  = 0.0
	alist = np.array(alist, "int32")
	blist = np.array(blist, "int32")	
	for iloop in range(Nloop):
		tlist = []
		clist = [[] for i in range(k)]
		for i in range(NT):
			new_clusters1 = []
			new_clusters2 = []
			np.random.shuffle(alist)
			np.random.shuffle(blist)	
			for j in range(k):new_clusters1.append(alist[plist1[j][0]:plist1[j][1]])
			for j in range(k):new_clusters2.append(blist[plist2[j][0]:plist2[j][1]])
			for j in range(k): new_clusters1[j] = np.sort(new_clusters1[j])
			for j in range(k): new_clusters2[j] = np.sort(new_clusters2[j])
			newindeces, list_stable, nb_tot_objs = sp_statistics.k_means_match_clusters_asg_new(new_clusters1,new_clusters2)
			ts = []
			for ii in range(len(newindeces)):
				ts += list(set(new_clusters1[newindeces[ii][0]].tolist() + new_clusters2[newindeces[ii][1]].tolist()))
			tlist.append(nb_tot_objs/float(len(set(ts)))*100.)
			for j in range(k):
				try:
					clist[j].append(float((np.intersect1d(new_clusters1[newindeces[j][0]], \
					      new_clusters2[newindeces[j][1]])).size)/float((np.union1d(new_clusters1[newindeces[j][0]], \
					      new_clusters2[newindeces[j][1]])).size)*100.)
				except: pass		
		svar +=sp_statistics.table_stat(tlist)[0]
		save +=sp_statistics.table_stat(tlist)[1]
		for j in range(k):
			try: gvar[j] +=sp_statistics.table_stat(clist[j])[1]
			except: gvar[j] +=0.0
			try: gave[j] +=sp_statistics.table_stat(clist[j])[0]
			except: gave[j] +=0.0
			
	save = mpi.mpi_reduce(save, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, 0, mpi.MPI_COMM_WORLD)
	svar = mpi.mpi_reduce(svar, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, 0, mpi.MPI_COMM_WORLD)
	
	for iproc in range(1, Blockdata["nproc"]):
		if Blockdata["myid"] == iproc:
			sp_utilities.wrap_mpi_send(gvar,  Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			
		elif Blockdata["myid"] == Blockdata["main_node"]:
			dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
			for im in range(len(dummy)): gvar[im]+=dummy[im]
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

	for iproc in range(1, Blockdata["nproc"]):
		if Blockdata["myid"] == iproc:
			sp_utilities.wrap_mpi_send(gave,  Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			
		elif Blockdata["myid"] == Blockdata["main_node"]:
			dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
			for im  in  range(len(dummy)): gave[im]+=dummy[im]
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

	if(Blockdata["myid"] == Blockdata["main_node"]):
		svar = numpy.sqrt(max(svar, 0.0)/float(Blockdata["nproc"]*Nloop))
		save = save/float(Blockdata["nproc"]*Nloop)
		for i in range(k):
			gave[i] = gave[i]/float(Blockdata["nproc"]*Nloop)
			gvar[i] = numpy.sqrt(max(gvar[i], 0.0)/float(Blockdata["nproc"]*Nloop))
		gave.append(save)
		gvar.append(svar)
	else:
		gave = 0
		gvar = 0
	gave = sp_utilities.wrap_mpi_bcast(gave, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	gvar = sp_utilities.wrap_mpi_bcast(gvar, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	return gave, gvar
###------------------------------------------------
def compare_two_iterations(assignment1, assignment2):
	# compare two assignments during clustering, either iteratively or independently
	if len(assignment1) !=len(assignment2):
		sp_global_def.ERROR("Two assignments don't have the same length", "compare_two_iterations", 0, 0)
	else:
		N          = len(assignment1)
		full_list  = np.array(range(N), dtype = np.int32)
	res1 = []
	assignment1 = np.array(assignment1,  dtype = np.int32)
	group1_ids  = np.unique(assignment1)
	for im in range(group1_ids.shape[0]):
		res1.append(np.sort(full_list[isin(assignment1, group1_ids[im])]))
	res2        = []
	assignment2 = np.array(assignment2, dtype = np.int32)
	group2_ids  = np.unique(assignment2)
	for im in range(group2_ids.shape[0]):
		res2.append(np.sort(full_list[isin(assignment2, group2_ids[im])]))
	newindeces, list_stable, nb_tot_objs = sp_statistics.k_means_match_clusters_asg_new(res1, res2)
	del res1
	del res2
	return float(nb_tot_objs)/N, newindeces, list_stable
	
def update_data_assignment(cdata, rdata, assignment, proc_list, nosmearing, myid):
	groupids = np.array(assignment[proc_list[myid][0]:proc_list[myid][1]], dtype=np.int64)
	if nosmearing:
		for im in range(len(cdata)):
			try:    previous_group = cdata[im].get_attr("group")
			except: previous_group = -1
			cdata[im].set_attr("group", groupids[im])
			rdata[im].set_attr("previous_group", previous_group)
			rdata[im].set_attr("group", groupids[im]) 
	else:
		for im in range(len(cdata)):
			try:    previous_group = cdata[im].get_attr("group")
			except: previous_group = -1
			cdata[im].set_attr("group", groupids[im])
			rdata[im][0].set_attr("previous_group", previous_group)
			rdata[im][0].set_attr("group", groupids[im])
	return
	
def get_angle_step_from_number_of_orien_groups(orien_groups):
	global Tracker, Blockdata
	N = orien_groups
	angle_step = 60.
	while len(Blockdata["symclass"].even_angles(angle_step))< N: angle_step /=2.
	while len(Blockdata["symclass"].even_angles(angle_step))> N: angle_step +=0.1
	return angle_step
	
def parti_oriens(params, angstep, smc):
	ntot = len(params)
	eah  = smc.even_angles(angstep, inc_mirror=0)
	leah = len(eah)
	u    = []
	for q in eah:
		m = smc.symmetry_related([(180.0+q[0])%360.0,180.0-q[1],0.0])
		itst = len(u)
		for c in m:
			if smc.is_in_subunit(c[0],c[1],1) :
				if not smc.is_in_subunit(c[0],c[1],0) :
					u.append(c)
					break
		if(len(u) != itst+1):
			u.append(q)  #  This is for exceptions that cannot be easily handled
	seaf = []
	for q in eah+u:seaf += smc.symmetry_related(q)
	lseaf = 2*leah
	seaf = sp_utilities.angles_to_normals(seaf)
	occupancy = [[] for i in range(leah)]
	for i,q in enumerate(params):
		l = sp_utilities.nearest_fang(seaf,q[0],q[1])
		l = l%lseaf
		if(l>=leah):  l = l-leah
		occupancy[l].append(i)
	return occupancy, eah
	
def get_angle_step_and_orien_groups_mpi(params_in, partids_in, angstep):
	global Tracker, Blockdata
	if Blockdata["main_node"] == Blockdata["myid"]:
		params  = sp_utilities.read_text_row(params_in)
		partids = sp_utilities.read_text_file(partids_in, -1)
		if len(partids) == 1: partids = partids[0]
		else:                 partids = partids[1]
		subparams = []
		for im in range(len(partids)): 
			subparams.append(params[partids[im]])
		occu, eah = parti_oriens(subparams, angstep, Blockdata["symclass"])
		ptls_in_orien_groups = []
		reassign_list        = []
		for l,q in enumerate(occu):
			if len(q)<Tracker["min_orien_group_size"]: reassign_list +=q
			else:      ptls_in_orien_groups.append(q)
		random.shuffle(reassign_list)
		for a in reassign_list:
			img = random.randint(0, len(ptls_in_orien_groups) - 1)
			ptls_in_orien_groups[img].append(a)
		del reassign_list
		for img in range(len(ptls_in_orien_groups)):
			tmp = sorted(ptls_in_orien_groups[img])
			ptls_in_orien_groups[img][:] = tmp[:]
	else: ptls_in_orien_groups = 0
	ptls_in_orien_groups = sp_utilities.wrap_mpi_bcast( \
	  ptls_in_orien_groups,Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	  
	return ptls_in_orien_groups
##### ======================================	
## conversion
###=====================>>>    rec3d for sorting    <<<=====================
def stepone(tvol, tweight):
	global Tracker, Blockdata
	tvol.set_attr("is_complex",1)
	ovol = EMAN2_cppwrap.Util.shrinkfvol(tvol,2)
	owol = EMAN2_cppwrap.Util.shrinkfvol(tweight,2)
	if( Tracker["constants"]["symmetry"] != "c1" ):
		ovol = ovol.symfvol(Tracker["constants"]["symmetry"], -1)
		owol = owol.symfvol(Tracker["constants"]["symmetry"], -1)
	return EMAN2_cppwrap.Util.divn_cbyr(ovol,owol)
	
def steptwo_mpi(tvol, tweight, treg, cfsc = None, regularized = True, color = 0):
	global Tracker, Blockdata
	if( Blockdata["color"] != color ):return sp_utilities.model_blank(1)  #  This should not be executed if called properly
	if( Blockdata["myid_on_node"] == 0 ):
		nz = tweight.get_zsize()
		ny = tweight.get_ysize()
		nx = tweight.get_xsize()
		tvol.set_attr_dict({"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1})
		if regularized:
			nr = len(cfsc)
			limitres = 0
			for i in range(nr):
				cfsc[i] = min(max(cfsc[i], 0.0), 0.999)
				if( cfsc[i] == 0.0 ):
					limitres = i-1
					break
			if( limitres == 0 ): limitres = nr-2;
			ovol = sp_utilities.reshape_1d(cfsc, nr, 2*nr)
			limitres = 2*min(limitres, Tracker["maxfrad"])  # 2 on account of padding, which is always on
			maxr2 = limitres**2
			for i in range(limitres+1, len(ovol), 1):   ovol[i] = 0.0
			ovol[0] = 1.0
			it = sp_utilities.model_blank(2*nr)
			for i in range(2*nr):  it[i] = ovol[i]
			del ovol
			#  Do not regularize first four
			for i in range(5):  treg[i] = 0.0
			EMAN2_cppwrap.Util.reg_weights(tweight, treg, it)
			del it
		else:
			limitres = 2*min(Tracker["constants"]["nnxo"]//2, Tracker["maxfrad"])
			maxr2 = limitres**2
		#  Iterative weights
		if( Tracker["constants"]["symmetry"] != "c1" ):
			tvol    = tvol.symfvol(Tracker["constants"]["symmetry"], limitres)
			tweight = tweight.symfvol(Tracker["constants"]["symmetry"], limitres)
		EMAN2_cppwrap.Util.divclreal(tvol, tweight, 1.0e-5)

	else:
		tvol = sp_utilities.model_blank(1)
		tweight = sp_utilities.model_blank(1)
		nz    = 0
		ny    = 0
		nx    = 0
		maxr2 = 0
	nx    = sp_utilities.bcast_number_to_all(nx, source_node = 0,    mpi_comm = Blockdata["shared_comm"])
	ny    = sp_utilities.bcast_number_to_all(ny, source_node = 0,    mpi_comm = Blockdata["shared_comm"])
	nz    = sp_utilities.bcast_number_to_all(nz, source_node = 0,    mpi_comm = Blockdata["shared_comm"])
	maxr2 = sp_utilities.bcast_number_to_all(maxr2, source_node = 0, mpi_comm = Blockdata["shared_comm"])
	"""Multiline Comment4"""
	if( Blockdata["myid_on_node"] == 0 ):
		sxprint(" iterefa  ",Blockdata["myid"],"   ",(time.time())/60.0)
		#  Either pad or window in F space to 2*nnxo
		nx = tvol.get_ysize()
		if( nx > 2*Tracker["constants"]["nnxo"]):
			tvol = sp_fundamentals.fdecimate(tvol, 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], \
			  2*Tracker["constants"]["nnxo"], False, False)
		elif(nx < 2*Tracker["constants"]["nnxo"]):
			tvol = sp_fundamentals.fpol(tvol, 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], \
			   2*Tracker["constants"]["nnxo"], RetReal = False, normalize = False)
		tvol = sp_fundamentals.fft(tvol)
		tvol = sp_fundamentals.cyclic_shift(tvol,Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
		tvol = EMAN2_cppwrap.Util.window(tvol, Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
		tvol.div_sinc(1)
		tvol = sp_morphology.cosinemask(tvol, Tracker["constants"]["nnxo"]//2-1,5, None)
		return tvol
	else:  return None
	
def steptwo_mpi_reg(tvol, tweight, treg, cfsc = None, cutoff_freq = 0.45, aa = 0.01, regularized = True, color = 0):
	global Tracker, Blockdata
	cutoff_freq2, aa = estimate_tanhl_params(cutoff_freq, aa, 2*Tracker["constants"]["nnxo"])
	if( Blockdata["color"] != color ):return sp_utilities.model_blank(1)  #  This should not be executed if called properly
	if( Blockdata["myid_on_node"] == 0 ):
		nz = tweight.get_zsize()
		ny = tweight.get_ysize()
		nx = tweight.get_xsize()
		tvol.set_attr_dict({"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1})
		if regularized:
			nr = len(cfsc)
			limitres = 0
			for i in range(nr):
				cfsc[i] = min(max(cfsc[i], 0.0), 0.999)
				if( cfsc[i] == 0.0 ):
					limitres = i-1
					break
			if( limitres == 0 ): limitres = nr-2;
			ovol = sp_utilities.reshape_1d(cfsc, nr, 2*nr)
			limitres = 2*min(limitres, Tracker["maxfrad"])  # 2 on account of padding, which is always on
			maxr2 = limitres**2
			for i in range(limitres+1, len(ovol), 1):   ovol[i] = 0.0
			ovol[0] = 1.0
			it = sp_utilities.model_blank(2*nr)
			for i in range(2*nr):  it[i] = ovol[i]
			del ovol
			#  Do not regularize first four
			for i in range(5):  treg[i] = 0.0
			EMAN2_cppwrap.Util.reg_weights(tweight, treg, it)
			del it
		else:
			limitres = 2*min(Tracker["constants"]["nnxo"]//2, Tracker["maxfrad"])
			maxr2 = limitres**2
		#  Iterative weights
		if( Tracker["constants"]["symmetry"] != "c1" ):
			tvol    = tvol.symfvol(Tracker["constants"]["symmetry"], limitres)
			tweight = tweight.symfvol(Tracker["constants"]["symmetry"], limitres)
		EMAN2_cppwrap.Util.divclreal(tvol, tweight, 1.0e-5)
	else:
		tvol = sp_utilities.model_blank(1)
		tweight = sp_utilities.model_blank(1)
		nz    = 0
		ny    = 0
		nx    = 0
		maxr2 = 0
	nx    = sp_utilities.bcast_number_to_all(nx, source_node = 0, mpi_comm = Blockdata["shared_comm"])
	ny    = sp_utilities.bcast_number_to_all(ny, source_node = 0, mpi_comm = Blockdata["shared_comm"])
	nz    = sp_utilities.bcast_number_to_all(nz, source_node = 0, mpi_comm = Blockdata["shared_comm"])
	maxr2 = sp_utilities.bcast_number_to_all(maxr2, source_node = 0, mpi_comm = Blockdata["shared_comm"])
	"""Multiline Comment5"""
	if( Blockdata["myid_on_node"] == 0 ):
		at = time.time()
		#print(" iterefa  ",Blockdata["myid"],"   ",strftime("%a, %d %b %Y %H:%M:%S", localtime()),"   ",(time()-at)/60.0)
		#  Either pad or window in F space to 2*nnxo
		if not cfsc: tvol = sp_filter.filt_tanl(tvol, cutoff_freq2, aa)
		nx = tvol.get_ysize()
		if( nx > 2*Tracker["constants"]["nnxo"]):
			tvol = sp_fundamentals.fdecimate(tvol, 2*Tracker["constants"]["nnxo"], \
			    2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], False, False)
		elif(nx < 2*Tracker["constants"]["nnxo"]):
			tvol = sp_fundamentals.fpol(tvol, 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"],\
			  2*Tracker["constants"]["nnxo"], RetReal = False, normalize = False)
		if Tracker["do_timing"]:
			sxprint(" pre fft  ",Blockdata["myid"],"   ",(time.time()-at)/60.0)
		fat = time.time()
		tvol = sp_fundamentals.fft(tvol)
		if Tracker["do_timing"]: 
			sxprint(" fft  ",Blockdata["myid"],"   ",(time.time()-fat)/60.0, tvol.get_xsize(), tvol.get_ysize(), tvol.get_zsize())
		tvol = sp_fundamentals.cyclic_shift(tvol,Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
		tvol = EMAN2_cppwrap.Util.window(tvol, Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
		tvol.div_sinc(1)
		tvol = sp_morphology.cosinemask(tvol, Tracker["constants"]["nnxo"]//2-1,5, None)
		if Tracker["do_timing"]:
			sxprint(" rec3d_rest  ",Blockdata["myid"],"   ",(time.time()-at)/60.0)
		return tvol
	else:  return None
####-----------------------------------------------
def recons3d_4nnsorting_MPI(myid, main_node, prjlist, random_subset, CTF = True, upweighted = True, mpi_comm= None, target_size=-1):
	"""
		recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
		Input
			list_of_prjlist: list of lists of projections to be included in the reconstruction
	"""
	if mpi_comm == None: mpi_comm = mpi.MPI_COMM_WORLD
	imgsize = prjlist[0].get_ysize()  # It can be Fourier, so take y-size
	refvol  = sp_utilities.model_blank(target_size)
	refvol.set_attr("fudge", 1.0)
	if CTF: do_ctf = 1
	else: do_ctf = 0
	fftvol = EMAN2_cppwrap.EMData()
	weight = EMAN2_cppwrap.EMData()
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, \
	   "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = EMAN2_cppwrap.Reconstructors.get( "nn4_ctfw", params )
	r.setup()
	#if norm_per_particle == None: norm_per_particle = len(prjlist)*[1.0]
	for im in range(len(prjlist)):
		phi, theta, psi, s2x, s2y = sp_utilities.get_params_proj(prjlist[im], xform = "xform.projection") # shifts are already applied 
		if random_subset == 2:
			bckgn = target_size*[1.]
			if prjlist[im].get_attr("is_complex") == 0:	prjlist[im] = sp_fundamentals.fft(prjlist[im]) 
			prjlist[im].set_attr_dict({"padffted":1, "is_complex":1})
			if not upweighted:  prjlist[im] = sp_filter.filt_table(prjlist[im], bckgn)
			prjlist[im].set_attr("bckgnoise", bckgn)
			r.insert_slice(prjlist[im], EMAN2_cppwrap.Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi}), 1.0)
		else:
			if prjlist[im].get_attr("chunk_id") == random_subset:
				#try:	bckgn = prjlist[im].get_attr("bckgnoise")
				bckgn = target_size*[1.]
				if prjlist[im].get_attr("is_complex")==0: prjlist[im] = sp_fundamentals.fft(prjlist[im])
				prjlist[im].set_attr_dict({"padffted":1, "is_complex":1})
				if not upweighted:  prjlist[im] = sp_filter.filt_table(prjlist[im], bckgn)
				prjlist[im].set_attr("bckgnoise", bckgn)
				r.insert_slice(prjlist[im], EMAN2_cppwrap.Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi}), 1.0)		
	#  clean stuff
	sp_utilities.reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
	sp_utilities.reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
	if myid == main_node: dummy = r.finish(True)
	mpi.mpi_barrier(mpi_comm)
	if myid == main_node: return fftvol, weight, refvol
	else: return None, None, None

def recons3d_4nnsorting_group_MPI(myid, main_node, prjlist, random_subset, group_ID, \
       CTF = True, upweighted = True, mpi_comm= None, target_size=-1):
	"""
		recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
		Input
			list_of_prjlist: list of lists of projections to be included in the reconstruction
	"""
	if mpi_comm == None: mpi_comm = mpi.MPI_COMM_WORLD
	imgsize = prjlist[0].get_ysize()  # It can be Fourier, so take y-size
	refvol = sp_utilities.model_blank(target_size)
	refvol.set_attr("fudge", 1.0)
	if CTF: do_ctf = 1
	else:   do_ctf = 0
	fftvol = EMAN2_cppwrap.EMData()
	weight = EMAN2_cppwrap.EMData()
	try:    qt = prjlist[0].get_attr("qt")
	except: qt = 1.0
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, "symmetry":"c1", \
	   "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = EMAN2_cppwrap.Reconstructors.get( "nn4_ctfw", params )
	r.setup()
	for im in range(len(prjlist)):
		phi, theta, psi, s2x, s2y = sp_utilities.get_params_proj(prjlist[im], xform = "xform.projection") # shifts are already applied
		if prjlist[im].get_attr("group") == group_ID:
			if random_subset == 2:
				try:	bckgn = prjlist[im].get_attr("bckgnoise")
				except:	bckgn = target_size*[1.]
				if prjlist[im].get_attr("is_complex") == 0:	image = sp_fundamentals.fft(prjlist[im])
				else: image =  prjlist[im].copy()
				image.set_attr_dict({"padffted":1, "is_complex":1})
				if not upweighted:  image = sp_filter.filt_table(image, bckgn)
				image.set_attr("bckgnoise", bckgn)
				r.insert_slice(image, EMAN2_cppwrap.Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi}), 1.0)
			else:
				if prjlist[im].get_attr("chunk_id") == random_subset:
					try:     bckgn = prjlist[im].get_attr("bckgnoise")
					except:	 bckgn = target_size*[1.]
					if prjlist[im].get_attr("is_complex")==0: image = sp_fundamentals.fft(prjlist[im])
					else: image =  prjlist[im].copy()
					image.set_attr_dict({"padffted":1, "is_complex":1})
					if not upweighted:  image = sp_filter.filt_table(image, bckgn)              
					image.set_attr("bckgnoise", bckgn)
					r.insert_slice(image, EMAN2_cppwrap.Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi}), 1.0)	
	sp_utilities.reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
	sp_utilities.reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
	if myid == main_node: dummy = r.finish(True)
	mpi.mpi_barrier(mpi_comm)
	if myid == main_node: return fftvol, weight, refvol
	else: return None, None, None

def do3d_sorting(procid, data, myid, mpi_comm = -1):
	global Tracker, Blockdata
	if (mpi_comm == -1): mpi_comm = mpi.MPI_COMM_WORLD
	if(procid == 0):
		if(Blockdata["no_of_groups"] >1):
			if(Blockdata["myid"] == Blockdata["nodes"][procid]):
				if os.path.exists(os.path.join(Tracker["directory"], "tempdir")): sxprint("tempdir exists")
				else:
					try: os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
					except:  sxprint("tempdir exists")
		else:
			if myid == Blockdata["main_node"]:
				if not os.path.exists(os.path.join(Tracker["directory"],"tempdir")): 
					try: os.mkdir(os.path.join(Tracker["directory"],    "tempdir"))
					except: sxprint("tempdir exists")
				else: sxprint("tempdir exists")
	mpi.mpi_barrier(mpi_comm)
	
	tvol, tweight, trol = recons3d_4nnsorting_MPI(myid = Blockdata["myid"], \
	   main_node = Blockdata["nodes"][procid], prjlist = data, random_subset = procid, \
	     CTF = Tracker["constants"]["CTF"], upweighted = False, target_size = (2*Tracker["nxinit"]+3), \
	         mpi_comm  = mpi_comm)
		
	if(Blockdata["no_of_groups"] >1):
		if(Blockdata["myid"] == Blockdata["nodes"][procid]):
			if(procid == 0):
				if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
					os.mkdir(os.path.join(Tracker["directory"],"tempdir"))
			tvol.set_attr("is_complex",0)
			tvol.write_image(os.path.join(Tracker["directory"],    "tempdir", "tvol_%01d.hdf"%procid))
			tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%01d.hdf"%procid))
			trol.write_image(os.path.join(Tracker["directory"],    "tempdir", "trol_%01d.hdf"%procid))
	else:
		if myid == Blockdata["main_node"]:
			tvol.set_attr("is_complex",0)
			tvol.write_image(os.path.join(Tracker["directory"],    "tempdir", "tvol_%01d.hdf"%procid))
			tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%01d.hdf"%procid))
			trol.write_image(os.path.join(Tracker["directory"],    "tempdir", "trol_%01d.hdf"%procid))
	mpi.mpi_barrier(mpi_comm)
	return
			
def do3d_sorting_group_insertion(data, randomset=2):
	global Tracker, Blockdata
	if(Blockdata["myid"] == Blockdata["last_node"]):
		if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
			os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
	if randomset ==1:
		for index_of_groups in range(Tracker["number_of_groups"]):
			for procid in range(2, 3):
				tvol, tweight, trol = recons3d_4nnsorting_group_MPI(myid = Blockdata["myid"], main_node = Blockdata["nodes"][0],\
				  prjlist = data,  random_subset = procid, group_ID = index_of_groups, CTF = Tracker["constants"]["CTF"],\
					upweighted = False, target_size = (2*Tracker["nxinit"]+3))
				
				if(Blockdata["myid"] == Blockdata["nodes"][procid]):
					tvol.set_attr("is_complex",0)
					tvol.write_image(os.path.join(Tracker["directory"],    "tempdir", "tvol_%d_%d.hdf"%(procid,    index_of_groups)))
					tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%d_%d.hdf"%(procid, index_of_groups)))
					trol.write_image(os.path.join(Tracker["directory"],    "tempdir", "trol_%d_%d.hdf"%(procid,    index_of_groups)))
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	else:
		for index_of_groups in range(Tracker["number_of_groups"]):
			for procid in range(2):
				tvol, tweight, trol = recons3d_4nnsorting_group_MPI(myid = Blockdata["myid"], \
				 main_node = Blockdata["nodes"][procid], prjlist = data, random_subset = procid, \
				     group_ID = index_of_groups, CTF = Tracker["constants"]["CTF"], upweighted = False, \
				         target_size = (2*Tracker["nxinit"]+3))
				
				if(Blockdata["myid"] == Blockdata["nodes"][procid]):
					tvol.set_attr("is_complex",0)
					tag =7007
					sp_utilities.send_EMData(tvol,    Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
					sp_utilities.send_EMData(tweight, Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
					sp_utilities.send_EMData(trol,    Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
					
				elif Blockdata["myid"] == Blockdata["last_node"]:
					tag =7007
					tvol    = sp_utilities.recv_EMData(Blockdata["nodes"][procid], tag, mpi.MPI_COMM_WORLD)
					tweight = sp_utilities.recv_EMData(Blockdata["nodes"][procid], tag, mpi.MPI_COMM_WORLD)
					trol    = sp_utilities.recv_EMData(Blockdata["nodes"][procid], tag, mpi.MPI_COMM_WORLD)
					tvol.write_image(os.path.join(Tracker["directory"],    "tempdir", "tvol_%d_%d.hdf"%(procid,    index_of_groups)))
					tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%d_%d.hdf"%(procid, index_of_groups)))
					trol.write_image(os.path.join(Tracker["directory"],    "tempdir", "trol_%d_%d.hdf"%(procid,    index_of_groups)))
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	return
	
def do3d_sorting_groups_trl_iter(data, iteration):
	global Tracker, Blockdata
	keepgoing = 1
	if(Blockdata["myid"] == Blockdata["last_node"]):
		if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
			os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	do3d_sorting_group_insertion(data)
	#####
	if Blockdata["no_of_groups"]>1:
		sub_main_node_list = [-1 for i in range(Blockdata["no_of_groups"])]
		for index_of_colors in range(Blockdata["no_of_groups"]):
			for iproc in range(Blockdata["nproc"]-1):
				if Blockdata["myid"]== iproc:
					if Blockdata["color"] == index_of_colors and Blockdata["myid_on_node"] == 0:
						sub_main_node_list[index_of_colors] = Blockdata["myid"]
					sp_utilities.wrap_mpi_send(sub_main_node_list, Blockdata["last_node"], mpi.MPI_COMM_WORLD)
				if Blockdata["myid"] == Blockdata["last_node"]:
					dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
					for im in range(len(dummy)):
						if dummy[im]>-1: sub_main_node_list[im] = dummy[im]
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		sp_utilities.wrap_mpi_bcast(sub_main_node_list, Blockdata["last_node"], mpi.MPI_COMM_WORLD)
		####		
		if Tracker["number_of_groups"]%Blockdata["no_of_groups"] == 0: 
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
		for iloop in range(nbig_loop):
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				if(Blockdata["myid"] == Blockdata["last_node"]):
					tvol2 		= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0_%d.hdf")%index_of_group)
					tweight2 	= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_0_%d.hdf")%index_of_group)
					treg2 		= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0_%d.hdf")%index_of_group)
					tag         = 7007
					sp_utilities.send_EMData(tvol2, sub_main_node_list[index_of_colors],    tag, mpi.MPI_COMM_WORLD)
					sp_utilities.send_EMData(tweight2, sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
					sp_utilities.send_EMData(treg2, sub_main_node_list[index_of_colors],    tag, mpi.MPI_COMM_WORLD)
				elif (Blockdata["myid"] == sub_main_node_list[index_of_colors]):
					tag      = 7007
					tvol2    = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
					tweight2 = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
					treg2    = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				Tracker["maxfrad"] = Tracker["nxinit"]//2
				if Blockdata["color"] == index_of_colors:
					if( Blockdata["myid_on_node"] != 0):
						tvol2 		= sp_utilities.model_blank(1)
						tweight2 	= sp_utilities.model_blank(1)
						treg2		= sp_utilities.model_blank(1)
					tvol2 = steptwo_mpi(tvol2, tweight2, treg2, Tracker["constants"]["fsc_curve"], \
					   True, color = index_of_colors) # has to be False!!!
					del tweight2, treg2
				mpi.mpi_barrier(Blockdata["shared_comm"])
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				if (Blockdata["color"] == index_of_colors) and (Blockdata["myid_on_node"] == 0):
					tag = 7007
					sp_utilities.send_EMData(tvol2, Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
				elif(Blockdata["myid"] == Blockdata["last_node"]):
					tag = 7007
					tvol2 = sp_utilities.recv_EMData(sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
					tvol2.write_image(os.path.join(Tracker["directory"], \
					  "vol_unfiltered_0_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
					del tvol2
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				if(Blockdata["myid"] == Blockdata["last_node"]):
					tvol2 		= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_1_%d.hdf")%index_of_group)
					tweight2 	= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_1_%d.hdf")%index_of_group)
					treg2 		= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "trol_1_%d.hdf"%index_of_group))
					tag      = 7007
					sp_utilities.send_EMData(tvol2,    sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
					sp_utilities.send_EMData(tweight2, sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
					sp_utilities.send_EMData(treg2,    sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
				elif (Blockdata["myid"] == sub_main_node_list[index_of_colors]):
					tag      = 7007
					tvol2    = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
					tweight2 = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
					treg2    = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				Tracker["maxfrad"] = Tracker["nxinit"]//2
				if Blockdata["color"] == index_of_colors:
					if( Blockdata["myid_on_node"] != 0):
						tvol2 		= sp_utilities.model_blank(1)
						tweight2 	= sp_utilities.model_blank(1)
						treg2		= sp_utilities.model_blank(1)
					tvol2 = steptwo_mpi(tvol2, tweight2, treg2, \
					  Tracker["constants"]["fsc_curve"], True, color = index_of_colors) # has to be False!!!
					del tweight2, treg2
					#if( Blockdata["myid_on_node"] == 0):
					#	tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
				mpi.mpi_barrier(Blockdata["shared_comm"])
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				if (Blockdata["color"] == index_of_colors) and (Blockdata["myid_on_node"] == 0):
					tag = 7007
					sp_utilities.send_EMData(tvol2, Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
				elif(Blockdata["myid"] == Blockdata["last_node"]):
					tag = 7007
					tvol2 = sp_utilities.recv_EMData(sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
					tvol2.write_image(os.path.join(Tracker["directory"], \
					  "vol_unfiltered_1_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
					del tvol2
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	else:
		Tracker["maxfrad"] = Tracker["nxinit"]//2
		for index_of_group in range(Tracker["number_of_groups"]):
			for iprocid in range(2):
				if(Blockdata["myid"] == Blockdata["last_node"]):
					tvol2 		= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_%d_%d.hdf")%(iprocid, index_of_group))
					tweight2 	= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_%d_%d.hdf")%(iprocid, index_of_group))
					treg2 		= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "trol_%d_%d.hdf"%(iprocid, index_of_group)))
					tag      = 7007
					sp_utilities.send_EMData(tvol2,    Blockdata["main_node"], tag, mpi.MPI_COMM_WORLD)
					sp_utilities.send_EMData(tweight2, Blockdata["main_node"], tag, mpi.MPI_COMM_WORLD)
					sp_utilities.send_EMData(treg2,    Blockdata["main_node"], tag, mpi.MPI_COMM_WORLD)
				elif (Blockdata["myid"] == Blockdata["main_node"]):
					tag      = 7007
					tvol2    = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
					tweight2 = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
					treg2    = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
				if( Blockdata["myid"] != Blockdata["main_node"]):
					tvol2 		= sp_utilities.model_blank(1)
					tweight2 	= sp_utilities.model_blank(1)
					treg2		= sp_utilities.model_blank(1)
				tvol2 = steptwo_mpi(tvol2, tweight2, treg2, None, False, color = 0) # has to be False!!!
				del tweight2, treg2
				if( Blockdata["myid"] == Blockdata["main_node"]):
					#tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
					tvol2.write_image(os.path.join(Tracker["directory"], \
					  "vol_unfiltered_%d_grp%03d_iter%03d.hdf"%(iprocid, index_of_group,iteration)))
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	keepgoing = sp_utilities.bcast_number_to_all(keepgoing, source_node = Blockdata["main_node"], mpi_comm = mpi.MPI_COMM_WORLD) # always check 
	Tracker = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"])
	if keepgoing == 0: sp_global_def.ERROR("do3d_sorting_groups_trl_iter  %s"%os.path.join(Tracker["directory"], \
	   "tempdir"),"do3d_sorting_groups_trl_iter", 1, Blockdata["myid"]) 
	return
###=======================================================================================
# Three ways of importing refinement results
def import_data(log_main):
	global Tracker, Blockdata
	# Two typical sorting scenarios
	# 1. import data and refinement parameters from meridien refinement;
	# 2. given data stack and xform.projection/ctf in header(For simulated test data);
	import_from_relion_refinement = 0
	import_from_sp_refinement  = 0
	import_from_data_stack		  = 0
	total_stack					  = 0
	####-----------------------------
	if Tracker["constants"]["refinement_method"] =="SPARX": # Senario one
		import_from_sp_refinement = get_input_from_sp_ref3d(log_main)
		Tracker["smearing"]          = True
	else:  # Senario three, sorting from a given data stack, general cases
		import_from_data_stack = get_input_from_datastack(log_main)
		Tracker["constants"]["hardmask"] = True
		Tracker["applybckgnoise"]        = False
		Tracker["applymask"]             = True
		Tracker["smearing"]              = False
	Tracker["total_stack"] = Tracker["constants"]["total_stack"]
	###=====<----------------====>checks=====---------------------------------------------
	if Tracker["constants"]["symmetry"] != Tracker["constants"]["sym"]:
		if(Blockdata["myid"] == Blockdata["main_node"]):
			msg = "Input symmetry %s is altered to %s after reading refinement information! "%\
			       (Tracker["constants"]["sym"], Tracker["constants"]["symmetry"])
			log_main.add(msg)

	## ------- checking settings!--------
	if not Tracker["search_mode"]:
		if Tracker["constants"]["total_stack"]//Tracker["constants"]["img_per_grp"]<=1: 
			sp_global_def.ERROR("Your input parameter img_per_grp is too large", "sxsort3d_depth.py", 1,  Blockdata["myid"])
	return
	
def get_input_from_sp_ref3d(log_main):# case one
	# import SPARX results
	global Tracker, Blockdata
	###-------------------------------
	import_from_sp_refinement = 1
	selected_iter                = 0
	Tracker_refinement           = 0
	checking_flag                = 0
	###-------------------------------
	if Blockdata["myid"] == Blockdata["main_node"]:
		if not os.path.exists (Tracker["constants"]["refinement_dir"]): checking_flag = 1
	checking_flag = sp_utilities.bcast_number_to_all(checking_flag, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	if checking_flag ==1:
		sp_global_def.ERROR("SPARX refinement directory does not exist", "get_input_from_sp_ref3d", 1, Blockdata["myid"])
	if Blockdata["myid"] == Blockdata["main_node"]:
		if Tracker["constants"]["niter_for_sorting"] == -1: # take the best solution to do sorting
			niter_refinement = 0
			while os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], \
			   "main%03d"%niter_refinement)) and os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"],\
			     "main%03d"%niter_refinement, "Tracker_%03d.json"%niter_refinement)): niter_refinement +=1
			niter_refinement -=1
			if niter_refinement !=0:
				with open(os.path.join(Tracker["constants"]["refinement_dir"], \
				   "main%03d"%niter_refinement, "Tracker_%03d.json"%niter_refinement),'r') as fout:
					Tracker_refinement 	= sp_utilities.convert_json_fromunicode(json.load(fout))
				fout.close()
				selected_iter = Tracker_refinement["constants"]["best"]
			else: import_from_sp_refinement = 0
		else:
			try:
				with open(os.path.join(Tracker["constants"]["refinement_dir"],\
				  "main%03d"%Tracker["constants"]["niter_for_sorting"], \
				"Tracker_%03d.json"%Tracker["constants"]["niter_for_sorting"]),'r') as fout:
					Tracker_refinement	= sp_utilities.convert_json_fromunicode(json.load(fout))
				fout.close()
				selected_iter = Tracker["constants"]["niter_for_sorting"]
			except:	import_from_sp_refinement = 0
	else: selected_iter = -1
	selected_iter = sp_utilities.bcast_number_to_all(selected_iter, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	import_from_sp_refinement = sp_utilities.bcast_number_to_all(import_from_sp_refinement, source_node = Blockdata["main_node"])
	if import_from_sp_refinement == 0:
		sp_global_def.ERROR("The best solution is not found","get_input_from_sp_ref3d", 1, Blockdata["myid"])
		mpi.mpi_finalize()
		exit()
	Tracker_refinement = sp_utilities.wrap_mpi_bcast(Tracker_refinement, Blockdata["main_node"], communicator = mpi.MPI_COMM_WORLD)
	# Check orgstack, set correct path
	if Blockdata["myid"] == Blockdata["main_node"]:
		Tracker["constants"]["orgstack"] = Tracker_refinement["constants"]["stack"]
		try:
			image = sp_utilities.get_im(Tracker["constants"]["orgstack"], 0)
		except Exception as e: 
			refinement_dir_path, refinement_dir_name = os.path.split(Tracker["constants"]["refinement_dir"])
			
			if Tracker_refinement["constants"]["stack"][0:4]=="bdb:":
				refinement_stack = "bdb:" + os.path.join(refinement_dir_path, Tracker_refinement["constants"]["stack"][4:])
			else:
				refinement_stack = os.path.join(refinement_dir_path, Tracker_refinement["constants"]["stack"]) # very rare case		
			
			if not Tracker["constants"]["orgstack"]:
				Tracker["constants"]["orgstack"] = refinement_stack
			
			try:    
				image = sp_utilities.get_im(Tracker["constants"]["orgstack"], 0)
			except: 
				import_from_sp_refinement = 0

		try:
			total_stack = EMUtil.get_image_count(Tracker["constants"]["orgstack"])
		except:
			total_stack = 0
	else: total_stack = 0
	import_from_sp_refinement = sp_utilities.bcast_number_to_all(import_from_sp_refinement, source_node = Blockdata["main_node"])
	if import_from_sp_refinement == 0:
		sp_global_def.ERROR("The data stack is not accessible","get_input_from_sp_ref3d",1, Blockdata["myid"])
	total_stack = sp_utilities.bcast_number_to_all(total_stack, source_node = Blockdata["main_node"])			
	Tracker["constants"]["total_stack"] = total_stack
	
	# Now copy relevant refinement files to sorting directory:
	if Blockdata["myid"] == Blockdata["main_node"]:
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], \
		    "main%03d"%selected_iter, "params_%03d.txt"%selected_iter)):
			shutil.copyfile( os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iter, \
			 "params_%03d.txt"%selected_iter), os.path.join(Tracker["constants"]["masterdir"], "sp_refinement_params.txt"))
		else: import_from_sp_refinement = 0
		Tracker["constants"]["selected_iter"] = selected_iter
	import_from_sp_refinement = sp_utilities.bcast_number_to_all(import_from_sp_refinement, source_node = Blockdata["main_node"])
	if import_from_sp_refinement == 0:
		sp_global_def.ERROR("The parameter file of the best solution is not accessible", "get_input_from_sp_ref3d", 1, Blockdata["myid"])
		
	if Blockdata["myid"] == Blockdata["main_node"]:
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iter, "bckgnoise.hdf")):
			shutil.copyfile(os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iter, "bckgnoise.hdf"),\
			os.path.join(Tracker["constants"]["masterdir"], "bckgnoise.hdf"))
		else:
			import_from_sp_refinement == 0
			for search_iter in range(selected_iter-1, 0, -1):
				if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%search_iter, "bckgnoise.hdf")):
					shutil.copyfile(os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%search_iter, \
					"bckgnoise.hdf"), os.path.join(Tracker["constants"]["masterdir"], "bckgnoise.hdf"))
					import_from_sp_refinement = 1
					break
	import_from_sp_refinement = sp_utilities.bcast_number_to_all(import_from_sp_refinement, source_node = Blockdata["main_node"])
	
	if import_from_sp_refinement == 0:
		Tracker["bckgnoise"] = None
		if Blockdata["myid"] == Blockdata["main_node"]:	
			if Tracker["do_timing"]: sxprint("Noise file is not found. However we continue")
	else: Tracker["bckgnoise"] = os.path.join(Tracker["constants"]["masterdir"], "bckgnoise.hdf")
	
	import_from_sp_refinement = 1	
	if Blockdata["myid"] == Blockdata["main_node"]:
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], \
		  "main%03d"%selected_iter, "driver_%03d.txt"%selected_iter)):
			shutil.copyfile(os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iter, \
			 "driver_%03d.txt"%selected_iter), os.path.join(Tracker["constants"]["masterdir"], "fsc_global.txt"))
		else: import_from_sp_refinement = 0
		
		if import_from_sp_refinement:
			fsc_curve = sp_utilities.read_text_row(os.path.join(Tracker["constants"]["masterdir"], "fsc_global.txt"))
		Tracker["constants"]["fsc_curve"] = \
		  sp_utilities.read_text_file(os.path.join(Tracker["constants"]["masterdir"], "fsc_global.txt"))
		# scale 2. * x / (1 + x)
		"""Multiline Comment6"""
	import_from_sp_refinement = sp_utilities.bcast_number_to_all(import_from_sp_refinement, source_node = Blockdata["main_node"])
	if import_from_sp_refinement == 0:
		sp_global_def.ERROR("The driver of the best solution is not accessible","get_input_from_sp_ref3d", 1, Blockdata["myid"])
	if Blockdata["myid"] == Blockdata["main_node"]:
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main000/indexes_000.txt")):
			shutil.copyfile(os.path.join(Tracker["constants"]["refinement_dir"], "main000/indexes_000.txt"), \
			os.path.join(Tracker["constants"]["masterdir"], "indexes.txt"))
		else: import_from_sp_refinement = 0
	import_from_sp_refinement = sp_utilities.bcast_number_to_all(import_from_sp_refinement, source_node = Blockdata["main_node"])
	
	if import_from_sp_refinement == 0:
		sp_global_def.ERROR("The index file of the best solution are not accessible","get_input_from_sp_ref3d", 1, Blockdata["myid"])
	if Blockdata["myid"] == Blockdata["main_node"]:
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main000/chunk_0_000.txt")):
			shutil.copyfile( os.path.join(Tracker["constants"]["refinement_dir"], "main000/chunk_0_000.txt"), \
			os.path.join(Tracker["constants"]["masterdir"], "chunk_0.txt"))
		else: import_from_sp_refinement == 0
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main000/chunk_1_000.txt")):
			shutil.copyfile(os.path.join(Tracker["constants"]["refinement_dir"], "main000/chunk_1_000.txt"), \
			os.path.join(Tracker["constants"]["masterdir"], "chunk_1.txt"))
		else: import_from_sp_refinement == 0
			
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main000/particle_groups_0.txt")):
			shutil.copyfile(os.path.join(Tracker["constants"]["refinement_dir"], "main000/particle_groups_0.txt"), \
			os.path.join(Tracker["constants"]["masterdir"], "particle_groups_0.txt"))
		else: import_from_sp_refinement == 0
			
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main000/particle_groups_1.txt")):
			shutil.copyfile( os.path.join(Tracker["constants"]["refinement_dir"], "main000/particle_groups_1.txt"), \
			os.path.join(Tracker["constants"]["masterdir"], "particle_groups_1.txt"))
		else:import_from_sp_refinement == 0
	import_from_sp_refinement = sp_utilities.bcast_number_to_all(import_from_sp_refinement, source_node = Blockdata["main_node"])
	if import_from_sp_refinement == 0:
		sp_global_def.ERROR("The chunk files and partice group files are not accessible","get_input_from_sp_ref3d",1, Blockdata["myid"])
	
	# copy all relavant parameters into sorting tracker
	if Blockdata["myid"] == Blockdata["main_node"]:
		if Tracker["constants"]["radius"] == -1: 
			Tracker["constants"]["radius"] = Tracker_refinement["constants"]["radius"]
		Tracker["constants"]["nnxo"]       = Tracker_refinement["constants"]["nnxo"]
		Tracker["constants"]["orgres"]     = Tracker_refinement["bestres"]
		Tracker["delta"]                   = Tracker_refinement["delta"]
		Tracker["ts"]                      = Tracker_refinement["ts"]
		Tracker["xr"]                      = Tracker_refinement["xr"]
		Tracker["constants"]["pixel_size"] = Tracker_refinement["constants"]["pixel_size"]
		Tracker["avgnorm"]                 = Tracker_refinement["avgvaradj"]
		if Tracker["constants"]["nxinit"]<0: 
			Tracker["nxinit_refinement"]   = Tracker_refinement["nxinit"] #Sphire window size
		else: Tracker["nxinit_refinement"] = Tracker["constants"]["nxinit"] #User defined window size
		
		try:     sym =  Tracker_refinement["constants"]["sym"]
		except:  sym =  Tracker_refinement["constants"]["symmetry"]
		if sym !='c1' and Tracker["constants"]["symmetry"] =='c1':
			Tracker["constants"]["symmetry"] = sym
			update_sym = 1
		else: update_sym = 0
		
		if not Tracker["constants"]["mask3D"]:
			if Tracker_refinement["constants"]["mask3D"] and (not Tracker["constants"]["do_not_use_3dmask"]):
				refinement_mask3D_path, refinement_mask3D_file = os.path.split(Tracker_refinement["constants"]["mask3D"])# MRK_DEBUG
				shutil.copyfile( os.path.join(refinement_dir_path, Tracker_refinement["constants"]["mask3D"]), \
				os.path.join(Tracker["constants"]["masterdir"], refinement_mask3D_file))
				Tracker["constants"]["mask3D"] = os.path.join(Tracker["constants"]["masterdir"], refinement_mask3D_file)
	else: 
		update_sym  = 0
		Tracker     = 0 
	Tracker    = sp_utilities.wrap_mpi_bcast(Tracker,         Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	update_sym = sp_utilities.bcast_number_to_all(update_sym, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	
	if update_sym ==1:
		Blockdata["symclass"] = sp_fundamentals.symclass(Tracker["constants"]["symmetry"])
		Tracker["orientation_groups"] = max(4, Tracker["constants"]["orientation_groups"]//Blockdata["symclass"].nsym)
	
	# Setting for margin error				
	chunk_dict = {}
	group_dict = {}
	if(Blockdata["myid"] == Blockdata["main_node"]):
		chunk_one = sp_utilities.read_text_file(os.path.join(Tracker["constants"]["masterdir"],"chunk_0.txt"))
		chunk_two = sp_utilities.read_text_file(os.path.join(Tracker["constants"]["masterdir"],"chunk_1.txt"))
	else:
		chunk_one = 0
		chunk_two = 0
	chunk_one = sp_utilities.wrap_mpi_bcast(chunk_one, Blockdata["main_node"])
	chunk_two = sp_utilities.wrap_mpi_bcast(chunk_two, Blockdata["main_node"])
	#
	if(Blockdata["myid"] == Blockdata["main_node"]):
		chunk_one_group = sp_utilities.read_text_file(os.path.join(Tracker["constants"]["masterdir"],"particle_groups_0.txt"))
		chunk_two_group = sp_utilities.read_text_file(os.path.join(Tracker["constants"]["masterdir"],"particle_groups_1.txt"))
	else:
		chunk_one_group = 0
		chunk_two_group = 0
	chunk_one_group = sp_utilities.wrap_mpi_bcast(chunk_one_group, Blockdata["main_node"])
	chunk_two_group = sp_utilities.wrap_mpi_bcast(chunk_two_group, Blockdata["main_node"])
	for index_of_element in range(len(chunk_one)): 
		chunk_dict[chunk_one[index_of_element]] = 0
		group_dict[chunk_one[index_of_element]] = chunk_one_group[index_of_element]
	for index_of_element in range(len(chunk_two)): 
		chunk_dict[chunk_two[index_of_element]] = 1
		group_dict[chunk_two[index_of_element]] = chunk_two_group[index_of_element] 			
	if(Blockdata["myid"] == Blockdata["main_node"]):
		chunk_ids = []
		group_ids = []
		partids   = sp_utilities.read_text_file(os.path.join(Tracker["constants"]["masterdir"], "indexes.txt"),-1)
		partids   = partids[0]
		Tracker["constants"]["total_stack"] = len(partids)
		params = sp_utilities.read_text_file(os.path.join(Tracker["constants"]["masterdir"], "sp_refinement_params.txt"),-1)
		for index_of_particle in range(len(partids)): 
			chunk_ids.append(chunk_dict[partids[index_of_particle]])
			group_ids.append(group_dict[partids[index_of_particle]])
		refinement_params = [params[0], params[1], params[2], params[3], params[4], chunk_ids, group_ids, params[7]]
		sp_utilities.write_text_file(refinement_params, os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt"))
	else: Tracker["constants"]["total_stack"] = 0
	Tracker["constants"]["total_stack"] = sp_utilities.bcast_number_to_all(Tracker["constants"]["total_stack"], \
	    Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	Tracker["total_stack"]            = Tracker["constants"]["total_stack"]
	Tracker["constants"]["partstack"] = os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt")
	total_stack          = Tracker["constants"]["total_stack"]
	Tracker["bckgnoise"] =  os.path.join(Tracker["constants"]["masterdir"], "bckgnoise.hdf")
	###
	# ------- Now copy oldparamstruture -------------------
	copy_oldparamstructure_from_meridien_MPI(selected_iter)
	###
	if Tracker["constants"]["check_smearing"]:
		if not Tracker["search_mode"]:
			number_of_groups  = Tracker["constants"]["total_stack"]//Tracker["constants"]["img_per_grp"]
		else: number_of_groups = Tracker["constants"]["init_K"]
			 
		avg_fsc = sp_statistics.scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], \
			 float(Tracker["constants"]["total_stack"]), \
			 float(Tracker["constants"]["total_stack"])//number_of_groups)
		fsc143 = get_res143(avg_fsc)
		if fsc143 !=0: 
			nxinit = min((int(fsc143)+ max(int(Tracker["constants"]["nnxo"]*0.03), 5))*2, Tracker["constants"]["nnxo"])
		else:
			sp_global_def.ERROR("Error in compute res143", "get_input_from_sp_ref3d", 1, Blockdata["myid"])
		nxinit = sp_fundamentals.smallprime(nxinit)
		nxinit = nxinit - nxinit%2 
		cdata_in_core =(Tracker["constants"]["total_stack"]*nxinit*nxinit*4.0)/1.e9/Blockdata["no_of_groups"]
		if not Tracker["constants"]["focus3D"]:fdata_in_core = 0.0
		else: fdata_in_core = cdata_in_core
		ctfdata = cdata_in_core
		refvol_size = (Tracker["nxinit"]*Tracker["nxinit"]*Tracker["nxinit"]*4.0*2)/1.e9 # including the 3D mask
		iterations  = 2
		if(Blockdata["myid"] == Blockdata["main_node"]):
			while os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%iterations)):
				iterations +=1
			log_main.add("---------->>> Smearing summary of %s  <<<----------"%Tracker["constants"]["refinement_dir"])
			log_main.add("Iter   average smearing    precalculated data/node     Percents of memory/node")
		else: iterations = 0
		iterations = sp_utilities.bcast_number_to_all(iterations, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		if Tracker["constants"]["memory_per_node"] ==-1.:
			try:
				mem_bytes = os.sysconf('SC_PAGE_SIZE')*os.sysconf('SC_PHYS_PAGES')# e.g. 4015976448
				memory_per_node = mem_bytes/(1024.**3) # e.g. 3.74
			except: memory_per_node = 32.
		else: memory_per_node = Tracker["constants"]["memory_per_node"]
		smearing_all_iters    = []
		for iteration in range(2, iterations):
			iter_smearings = get_smearing_info(Blockdata["nproc_previous"], \
			   iteration, Tracker["constants"]["total_stack"], \
			   Tracker["constants"]["masterdir"], Tracker["constants"]["refinement_dir"])
			savg           = np.sum(iter_smearings)/Tracker["constants"]["total_stack"]
			srdata_in_core = (nxinit*nxinit*np.sum(iter_smearings)*4.)/1.e9/Blockdata["no_of_groups"]
			tdata          = cdata_in_core+srdata_in_core+ctfdata+refvol_size
			percnt         = tdata/memory_per_node*100.
			smearing_all_iters.append([iteration, savg, tdata, percnt])
			if(Blockdata["myid"] == Blockdata["main_node"]):
				log_main.add("%5d        %5.1f             %7.1f GB                  %7.1f%% "%(iteration, savg, tdata, percnt))
		if(Blockdata["myid"] == Blockdata["main_node"]):
			log_main.add('================================================================================================================')
			sp_utilities.write_text_row(smearing_all_iters, os.path.join(Tracker["constants"]["masterdir"], "smearing_all.txt"))
	return import_from_sp_refinement
		
def get_input_from_datastack(log_main):# Case three
	global Tracker, Blockdata
	import_from_data_stack = 1
	if(Blockdata["myid"] == Blockdata["main_node"]):
		image = sp_utilities.get_im(Tracker["constants"]["orgstack"], 0)
		Tracker["constants"]["nnxo"] = image.get_xsize()		
		if( Tracker["nxinit"] > Tracker["constants"]["nnxo"]):
				sp_global_def.ERROR("Image size less than minimum permitted %d"%Tracker["nxinit"], \
				  "get_input_from_datastack",1, Blockdata["myid"])
				nnxo = -1
		else:
			if Tracker["constants"]["CTF"]:
				ictf = image.get_attr('ctf')
				Tracker["constants"]["pixel_size"] = ictf.apix
			else:
				Tracker["constants"]["pixel_size"] = 1.0
				del image
	else:
		Tracker["constants"]["nnxo"]       = 0
		Tracker["constants"]["pixel_size"] = 1.0
	Tracker["constants"]["nnxo"] = \
	    sp_utilities.bcast_number_to_all(Tracker["constants"]["nnxo"], source_node = Blockdata["main_node"])
	if( Tracker["constants"]["nnxo"] < 0):
		sp_global_def.ERROR("Image size is negative", "get_input_from_datastack", 1, Blockdata["main_node"])
	Tracker["constants"]["pixel_size"]	= \
	    sp_utilities.bcast_number_to_all(Tracker["constants"]["pixel_size"], source_node =  Blockdata["main_node"])
	if(Tracker["constants"]["radius"] < 1):
		Tracker["constants"]["radius"]  = Tracker["constants"]["nnxo"]//2-2
	elif((2*Tracker["constants"]["radius"] +2) > Tracker["constants"]["nnxo"]):
		sp_global_def.ERROR("Particle radius set too large!", "get_input_from_datastack",1, Blockdata["myid"])
	if Blockdata["myid"] == Blockdata["main_node"]:
		total_stack = EMAN2_cppwrap.EMUtil.get_image_count(Tracker["constants"]["orgstack"])
	else: total_stack = 0
	total_stack = sp_utilities.bcast_number_to_all(total_stack, Blockdata["main_node"])
	# randomly assign two subsets
	Tracker["constants"]["total_stack"]	  = total_stack
	Tracker["constants"]["chunk_0"]		  = os.path.join(Tracker["constants"]["masterdir"],"chunk_0.txt")
	Tracker["constants"]["chunk_1"]		  = os.path.join(Tracker["constants"]["masterdir"],"chunk_1.txt")
	Tracker["constants"]["partstack"]	  = os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt")
	Tracker["previous_parstack"]          = os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt")#
	Tracker["constants"]["smearing_file"] = os.path.join(os.path.join(Tracker["constants"]["masterdir"], "all_smearing.txt"))
	
	if Blockdata["myid"] == Blockdata["main_node"]:
		sp_utilities.write_text_file(np.full(Tracker["constants"]["total_stack"], 1, dtype = np.int16).tolist(),\
		   Tracker["constants"]["smearing_file"])
		   
	Tracker["refang"], Tracker["rshifts"], Tracker["delta"] = None, None, None
	Tracker["avgnorm"] = 1.0
	chunk_dict         = {}
	chunk_list         = []
	if Blockdata["myid"] == Blockdata["main_node"]:
		chunk_dict  = {}
		tlist = list(range(total_stack))
		sp_utilities.write_text_file(tlist, os.path.join(Tracker["constants"]["masterdir"], "indexes.txt"))
		random.shuffle(tlist)
		chunk_one	= tlist[0:total_stack//2]
		chunk_two	= tlist[total_stack//2:]
		chunk_one	= sorted(chunk_one)
		chunk_two	= sorted(chunk_two)
		sp_utilities.write_text_row(chunk_one,Tracker["constants"]["chunk_0"])
		sp_utilities.write_text_row(chunk_two,Tracker["constants"]["chunk_1"])
		for particle in chunk_one: chunk_dict[particle] = 0
		for particle in chunk_two: chunk_dict[particle] = 1
		xform_proj_list = EMAN2_cppwrap.EMUtil.get_all_attributes(Tracker["constants"]["orgstack"], "xform.projection")
		for index_of_particle in range(len(xform_proj_list)):
			dp = xform_proj_list[index_of_particle].get_params("spider")
			xform_proj_list[index_of_particle] = \
			    [dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"], chunk_dict[index_of_particle]]
		sp_utilities.write_text_row(xform_proj_list, Tracker["constants"]["partstack"])
	else:
		chunk_one = 0
		chunk_two = 0
	chunk_one = sp_utilities.wrap_mpi_bcast(chunk_one, Blockdata["main_node"])
	chunk_two = sp_utilities.wrap_mpi_bcast(chunk_two, Blockdata["main_node"])
	for element in chunk_one: chunk_dict[element] = 0
	for element in chunk_two: chunk_dict[element] = 1
	chunk_list 				= [chunk_one, chunk_two]	
	# Reconstruction to determine the resolution in orignal data size
	Tracker["nxinit"]     = Tracker["constants"]["nnxo"]
	Tracker["shrinkage"]  = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	Tracker["bckgnoise"]  =  None
	temp = sp_utilities.model_blank(Tracker["constants"]["nnxo"], Tracker["constants"]["nnxo"])
	nny  = temp.get_ysize()
	Blockdata["bckgnoise"] = [1.0]*nny # set for initial recon3D of data from stack	
	Tracker["fuse_freq"]   = int(Tracker["constants"]["pixel_size"]*Tracker["constants"]["nnxo"]/Tracker["constants"]["fuse_freq"] +0.5)
	Tracker["directory"]   = Tracker["constants"]["masterdir"]
	
	if Tracker["constants"]["nxinit"]< 0: Tracker["nxinit_refinement"] = Tracker["constants"]["nnxo"]
	else:                                 Tracker["nxinit_refinement"] = Tracker["constants"]["nxinit"]
	
	for procid in range(2):
		data = get_shrink_data_sorting(os.path.join(Tracker["constants"]["masterdir"], \
		   "chunk_%01d.txt"%procid), Tracker["constants"]["partstack"])
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		do3d_sorting(procid, data, myid = Blockdata["myid"],  mpi_comm = mpi.MPI_COMM_WORLD)# 1
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	
	if(Blockdata["no_of_groups"] == 1):
		if( Blockdata["myid"] == Blockdata["main_node"]):
			tvol0 		= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0.hdf"))
			tweight0 	= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_0.hdf"))
			tvol1 		= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_1.hdf"))
			tweight1 	= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_1.hdf"))
			EMAN2_cppwrap.Util.fuse_low_freq(tvol0, tvol1, tweight0, tweight1, 2*Tracker["fuse_freq"])
			shrank0 	= stepone(tvol0, tweight0)
			shrank1 	= stepone(tvol1, tweight1)
		    #  Note shrank volumes are Fourier uncentered.
			cfsc 		= sp_statistics.fsc(shrank0, shrank1)[1]
			del shrank0, shrank1
			if(Tracker["nxinit"]<Tracker["constants"]["nnxo"]): cfsc  = cfsc[:Tracker["nxinit"]]
			for i in range(len(cfsc),Tracker["constants"]["nnxo"]//2+1): cfsc.append(0.0)
			sp_utilities.write_text_row(cfsc, os.path.join(Tracker["directory"], "fsc_global.txt"))
			lcfsc  = len(cfsc)							
			Tracker["constants"]["fsc_curve"]  = copy.copy(cfsc)
			# scale 2. * x / (1 + x)
			"""Multiline Comment7"""
		else: Tracker = 0
		Tracker = sp_utilities.wrap_mpi_bcast(Tracker,Blockdata["nodes"][0], communicator = mpi.MPI_COMM_WORLD)
	else:
		if(Blockdata["myid"] == Blockdata["nodes"][1]):  # It has to be 1 to avoid problem with tvol1 not closed on the disk
			tvol0 		= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0.hdf"))
			tweight0 	= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_0.hdf"))
			tvol1 		= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_1.hdf"))
			tweight1 	= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_1.hdf"))
			EMAN2_cppwrap.Util.fuse_low_freq(tvol0, tvol1, tweight0, tweight1, 2*Tracker["fuse_freq"])
			tag = 7007
			sp_utilities.send_EMData(tvol1, Blockdata["nodes"][0],    tag, mpi.MPI_COMM_WORLD)
			sp_utilities.send_EMData(tweight1, Blockdata["nodes"][0], tag, mpi.MPI_COMM_WORLD)
			shrank0 	= stepone(tvol0, tweight0)
			sp_utilities.send_EMData(shrank0, Blockdata["nodes"][0], tag, mpi.MPI_COMM_WORLD)
			del shrank0
			lcfsc = 0
		
		elif( Blockdata["myid"] == Blockdata["nodes"][0]):
			tag = 7007
			tvol1 		= sp_utilities.recv_EMData(Blockdata["nodes"][1], tag, mpi.MPI_COMM_WORLD)
			tweight1 	= sp_utilities.recv_EMData(Blockdata["nodes"][1], tag, mpi.MPI_COMM_WORLD)
			tvol1.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1} )
			shrank1 	= stepone(tvol1, tweight1)
			#  Get shrank volume, do fsc, send it to all
			shrank0 	= sp_utilities.recv_EMData(Blockdata["nodes"][1], tag, mpi.MPI_COMM_WORLD)
			#  Note shrank volumes are Fourier uncentered.
			cfsc 		= sp_statistics.fsc(shrank0, shrank1)[1]
			sp_utilities.write_text_row(cfsc, os.path.join(Tracker["directory"], "fsc_global.txt"))
			del shrank0, shrank1
			if(Tracker["nxinit"]<Tracker["constants"]["nnxo"]):
				cfsc  = cfsc[:Tracker["nxinit"]]
				for i in range(len(cfsc),Tracker["constants"]["nnxo"]//2+1): cfsc.append(0.0)							
			Tracker["constants"]["fsc_curve"]  = copy.copy(cfsc)
			for ifreq in range(len(Tracker["constants"]["fsc_curve"])):
				Tracker["constants"]["fsc_curve"][ifreq] = \
				 Tracker["constants"]["fsc_curve"][ifreq]*2./(1.+Tracker["constants"]["fsc_curve"][ifreq])
		Tracker = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["nodes"][0], communicator = mpi.MPI_COMM_WORLD)
	return import_from_data_stack
####=======================================================================================================	
### functions for faked rec3d from subsets
###
def do3d_sorting_groups_nofsc_smearing_iter(srdata, paramstructure, norm_per_particle,\
           partial_rec3d, current_group_sizes, iteration):
	global Tracker, Blockdata
	at = time.time()
	keepgoing = 1
	if(Blockdata["myid"] == Blockdata["last_node"]):
		if (not os.path.exists(os.path.join(Tracker["directory"], "tempdir"))):
			os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
		try:
			with open(os.path.join(Tracker["directory"],"freq_cutoff.json"),'r') as fout:
				freq_cutoff_dict = sp_utilities.convert_json_fromunicode(json.load(fout))
			fout.close()
		except: freq_cutoff_dict = 0
	else: freq_cutoff_dict = 0
	freq_cutoff_dict = sp_utilities.wrap_mpi_bcast(freq_cutoff_dict, Blockdata["last_node"], mpi.MPI_COMM_WORLD)
	rest_time  = time.time()
	for index_of_groups in range(Tracker["number_of_groups"]):
		if partial_rec3d:
			tvol, tweight, trol = recons3d_trl_struct_group_nofsc_shifted_data_partial_MPI(\
			 Blockdata["myid"], Blockdata["last_node"], Blockdata["nproc"], \
			 srdata, paramstructure, norm_per_particle, index_of_groups, \
			os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_groups), \
			os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf"%index_of_groups), \
			os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf"%index_of_groups),\
      		None,  Tracker["constants"]["CTF"], (2*Tracker["nxinit"]+3), Tracker["nosmearing"])
		else:
			tvol, tweight, trol = recons3d_trl_struct_group_nofsc_shifted_data_MPI(\
			    Blockdata["myid"], Blockdata["last_node"],\
			     srdata, paramstructure, norm_per_particle,\
			     index_of_groups, None,  Tracker["constants"]["CTF"], \
			     (2*Tracker["nxinit"]+3), Tracker["nosmearing"])
			     
		if(Blockdata["myid"] == Blockdata["last_node"]):
			tvol.set_attr("is_complex",0)
			tvol.write_image(os.path.join(Tracker["directory"],    "tempdir", "tvol_2_%d.hdf"%index_of_groups))
			tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf"%index_of_groups))
			trol.write_image(os.path.join(Tracker["directory"],    "tempdir", "trol_2_%d.hdf"%index_of_groups))
			del tvol
			del tweight
			del trol
			if Tracker["do_timing"]:
				sxprint("volumes written    ",(time.time()-at)/60.)
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	acc_rest = time.time() - rest_time
	if Blockdata["myid"] == Blockdata["main_node"]:
		if Tracker["do_timing"]:
			sxprint("insertion  take  %f minutes "%(acc_rest/60.))
	Tracker["maxfrad"] = Tracker["nxinit"]//2
	fsc_groups = [None for im in range(len(current_group_sizes))]
	total_size = sum(current_group_sizes)
	for igrp in range(len(current_group_sizes)):
		if Tracker["even_scale_fsc"]: 
			fsc_groups[igrp] = sp_statistics.scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], \
		     float(Tracker["constants"]["total_stack"]), total_size//len(current_group_sizes))
		else:
			fsc_groups[igrp] = sp_statistics.scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], \
		     float(Tracker["constants"]["total_stack"]), current_group_sizes[igrp])
	if (Blockdata["no_of_groups"]>1):
	 	# new starts
	 	rest_time  = time.time()
		sub_main_node_list = [ -1 for i in range(Blockdata["no_of_groups"])]
		for index_of_colors in range(Blockdata["no_of_groups"]):
			for iproc in range(Blockdata["nproc"]-1):
				if (Blockdata["myid"]== iproc):
					if (Blockdata["color"] == index_of_colors) and (Blockdata["myid_on_node"] == 0):
						sub_main_node_list[index_of_colors] = Blockdata["myid"]
					sp_utilities.wrap_mpi_send(sub_main_node_list, Blockdata["last_node"], mpi.MPI_COMM_WORLD)
				if (Blockdata["myid"] == Blockdata["last_node"]):
					dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
					for im in range(len(dummy)):
						if (dummy[im]>-1): sub_main_node_list[im] = dummy[im]
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		sp_utilities.wrap_mpi_bcast(sub_main_node_list, Blockdata["last_node"], mpi.MPI_COMM_WORLD)
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
		for iloop in range(nbig_loop):
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
			
				if(Blockdata["myid"] == Blockdata["last_node"]):
					tvol2    = sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")%index_of_group)
					tweight2 = sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")%index_of_group)
					treg2 	 = sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_group))
					tvol2.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1})
					tag = 7007
					sp_utilities.send_EMData(tvol2,    sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
					sp_utilities.send_EMData(tweight2, sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
					sp_utilities.send_EMData(treg2,    sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
				
				elif (Blockdata["myid"] == sub_main_node_list[index_of_colors]):
					tag      = 7007
					tvol2    = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
					tweight2 = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
					treg2    = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				try: 
					Tracker["freq_fsc143_cutoff"] = freq_cutoff_dict["Cluster_%03d.txt"%index_of_group]
				except: pass	
				if Blockdata["color"] == index_of_colors:  # It has to be 1 to avoid problem with tvol1 not closed on the disk 							
					if(Blockdata["myid_on_node"] != 0): 
						tvol2     = sp_utilities.model_blank(1)
						tweight2  = sp_utilities.model_blank(1)
						treg2     = sp_utilities.model_blank(1)
					tvol2 = steptwo_mpi_reg(tvol2, tweight2, treg2,  fsc_groups[big_loop_groups[iloop][im]], \
					         Tracker["freq_fsc143_cutoff"], 0.01, True, color = index_of_colors) # has to be False!!!
					del tweight2, treg2
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				if (Blockdata["color"] == index_of_colors) and (Blockdata["myid_on_node"] == 0):
					tag = 7007
					sp_utilities.send_EMData(tvol2, Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
					
				elif(Blockdata["myid"] == Blockdata["last_node"]):
					tag = 7007
					tvol2    = sp_utilities.recv_EMData(sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
					tvol2.write_image(os.path.join(Tracker["directory"], "vol_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
					del tvol2
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		acc_rest = time.time() - rest_time
		if Blockdata["myid"] == Blockdata["main_node"]:
			if Tracker["do_timing"]: sxprint("rec3d steptwo  take  %f minutes "%(acc_rest/60.))
		
	else:# loop over all groups for single node
		for index_of_group in range(Tracker["number_of_groups"]):
			if(Blockdata["myid_on_node"] == 0):
				tvol2 		= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")%index_of_group)
				tweight2 	= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")%index_of_group)
				treg2 		= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_group))
				tvol2.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1})
			else:
				tvol2 		= sp_utilities.model_blank(1)
				tweight2 	= sp_utilities.model_blank(1)
				treg2		= sp_utilities.model_blank(1)
			try: 
				Tracker["freq_fsc143_cutoff"] = freq_cutoff_dict["Cluster_%03d.txt"%index_of_group]
			except: pass
			
			tvol2 = steptwo_mpi_reg(tvol2, tweight2, treg2, fsc_groups[index_of_group],\
			      Tracker["freq_fsc143_cutoff"], 0.01, True) # has to be False!!!
			del tweight2, treg2
			if(Blockdata["myid_on_node"] == 0):
				tvol2.write_image(os.path.join(Tracker["directory"], "vol_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
				del tvol2
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	keepgoing = sp_utilities.bcast_number_to_all(keepgoing, source_node = Blockdata["main_node"], mpi_comm = mpi.MPI_COMM_WORLD) # always check 
	Tracker   = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"])
	if (not keepgoing):
		sp_global_def.ERROR("do3d_sorting_groups_trl_iter  %s"%os.path.join(Tracker["directory"], "tempdir"),"do3d_sorting_groups_trl_iter", 1, Blockdata["myid"])
	if(Blockdata["myid"] == 0) and Tracker["do_timing"]:
		sxprint("Reconstructions done    ",(time.time()-at)/60.)
	return
	
### nofsc insertion #1
def recons3d_trl_struct_group_nofsc_shifted_data_partial_MPI(myid, main_node, nproc,\
         prjlist, paramstructure, norm_per_particle, group_ID,\
         refvol_file, fftvol_file, weight_file,\
         mpi_comm = None, CTF = True, target_size=-1, nosmearing = False):
	"""
		partial rec3d for re-assigned images
		reconstructor nn4_ctfws
	"""
	if mpi_comm == None: mpi_comm = mpi.MPI_COMM_WORLD
	if CTF: do_ctf = 1
	else:   do_ctf = 0
	if not os.path.exists(refvol_file):
		sp_global_def.ERROR("refvol does not exist", "recons3d_trl_struct_group_nofsc_shifted_data_partial_MPI", 1, myid)
	if not os.path.exists(fftvol_file):
		sp_global_def.ERROR("fftvol does not exist", "recons3d_trl_struct_group_nofsc_shifted_data_partial_MPI", 1, myid)
	if not os.path.exists(weight_file):
		sp_global_def.ERROR("weight does not exist", "recons3d_trl_struct_group_nofsc_shifted_data_partial_MPI", 1, myid)
	
	#refvol
	if myid == main_node: target_size = sp_utilities.get_im(refvol_file).get_xsize()
	else: target_size = 0
	target_size = sp_utilities.bcast_number_to_all(target_size, main_node, mpi_comm)
	refvol = sp_utilities.model_blank(target_size)# set to zero
	
	# fftvol
	if (myid == main_node): 
		fftvol = sp_utilities.get_im(fftvol_file)
		fftvol.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1} )
		fftvol /=float(nproc)
		#Util.mult_scalar(fftvol, 1./float(Blockdata["nproc"]))
		nxfft  =fftvol.get_xsize()
		nyfft  =fftvol.get_ysize()
		nzfft  =fftvol.get_zsize()
	else: 
		nxfft = 0
		nyfft = 0
		nzfft = 0
	nxfft = sp_utilities.bcast_number_to_all(nxfft, main_node, mpi_comm)
	nyfft = sp_utilities.bcast_number_to_all(nyfft, main_node, mpi_comm)
	nzfft = sp_utilities.bcast_number_to_all(nzfft, main_node, mpi_comm)
	if (myid!= main_node): fftvol = sp_utilities.model_blank(nxfft, nyfft, nzfft)
	sp_utilities.bcast_EMData_to_all(fftvol, myid, main_node)
	
	# weight
	if (myid == main_node): 
		weight = sp_utilities.get_im(weight_file)
		weight /=float(nproc)
		#Util.mult_scalar(weight, 1./float(Blockdata["nproc"]))
		nxweight  = weight.get_xsize()
		nyweight  = weight.get_ysize()
		nzweight  = weight.get_zsize()
	else: 
		nxweight = 0
		nyweight = 0
		nzweight = 0
	nxweight = sp_utilities.bcast_number_to_all(nxweight, main_node, mpi_comm)
	nyweight = sp_utilities.bcast_number_to_all(nyweight, main_node, mpi_comm)
	nzweight = sp_utilities.bcast_number_to_all(nzweight, main_node, mpi_comm)
	if(myid != main_node): weight = sp_utilities.model_blank(nxweight, nyweight, nzweight)
	sp_utilities.bcast_EMData_to_all(weight, myid, main_node)
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, \
	  "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = EMAN2_cppwrap.Reconstructors.get("nn4_ctfws", params)
	r.setup()
	rest_time  = time.time()
	if nosmearing:
		for im in range(len(prjlist)):
			current_group_ID  = prjlist[im].get_attr("group")
			previous_group_ID = prjlist[im].get_attr("previous_group")
			if current_group_ID !=previous_group_ID:
				if previous_group_ID == group_ID:
					r.insert_slice(prjlist[im], prjlist[im].get_attr("xform.projection"), -1)
				if current_group_ID  == group_ID:
					r.insert_slice(prjlist[im], prjlist[im].get_attr("xform.projection"), 1)
	else:
		nny = prjlist[0][0].get_ysize()
		if not Tracker["constants"]["compute_on_the_fly"]:
			num_on_the_fly = len(prjlist)
		else:
			num_on_the_fly = Tracker["num_on_the_fly"][Blockdata["myid"]]
		for im in range(num_on_the_fly):
			current_group_ID  = prjlist[im][0].get_attr("group")
			previous_group_ID = prjlist[im][0].get_attr("previous_group")
			if current_group_ID != previous_group_ID:
				if current_group_ID  == group_ID:
					for jm in range(len(prjlist[im])):
						r.insert_slice(prjlist[im][jm], prjlist[im][jm].get_attr("xform.projection"), \
							  prjlist[im][jm].get_attr("wprob"))
				if previous_group_ID == group_ID:
					for jm in range(len(prjlist[im])):
						r.insert_slice(prjlist[im][jm], prjlist[im][jm].get_attr("xform.projection"), \
							  prjlist[im][jm].get_attr("wprob")*(-1.))
		if (Tracker["num_on_the_fly"][Blockdata["myid"]]<len(prjlist)) and \
		         Tracker["constants"]["compute_on_the_fly"]:
			rshifts_shrank = copy.deepcopy(Tracker["rshifts"])
			for lm in range(len(rshifts_shrank)):
				rshifts_shrank[lm][0] *= float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
				rshifts_shrank[lm][1] *= float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
			refang =  Tracker["refang"]
			delta  =  Tracker["delta"]
			for im in range(num_on_the_fly, len(prjlist)):
				current_group_ID  = prjlist[im][0].get_attr("group")
				previous_group_ID = prjlist[im][0].get_attr("previous_group")
				if (current_group_ID != previous_group_ID):
					ct      = prjlist[im][0].get_attr("ctf")
					bckgn   = prjlist[im][0].get_attr("bckgnoise")
					avgnorm = Tracker["avgnorm"][prjlist[im][0].get_attr("chunk_id")]
					numbor      = len(paramstructure[im][2])
					ipsiandiang = [ paramstructure[im][2][i][0]/1000  for i in range(numbor) ]
					allshifts   = [ paramstructure[im][2][i][0]%1000  for i in range(numbor) ]
					probs       = [ paramstructure[im][2][i][1] for i in range(numbor) ]
					tdir = list(set(ipsiandiang))
					data = [None]*len(Tracker["rshifts"])
					for ii in range(len(tdir)):
						#  Find the number of times given projection direction appears on the list, it is the number of different shifts associated with it.
						lshifts = sp_utilities.findall(tdir[ii], ipsiandiang)
						toprab  = 0.0
						for ki in range(len(lshifts)):toprab += probs[lshifts[ki]]
						recdata = EMAN2_cppwrap.EMData(nny,nny,1,False)
						recdata.set_attr("is_complex",0)
						for ki in range(len(lshifts)):
							lpt = allshifts[lshifts[ki]]
							if(data[lpt] == None):
								data[lpt] = sp_fundamentals.fshift(prjlist[im][0], rshifts_shrank[lpt][0], rshifts_shrank[lpt][1])
								data[lpt].set_attr("is_complex",0)
							EMAN2_cppwrap.Util.add_img(recdata, EMAN2_cppwrap.Util.mult_scalar(data[lpt], probs[lshifts[ki]]/toprab))
						recdata.set_attr_dict({"padffted":1, "is_fftpad":1,"is_fftodd":0, "is_complex_ri":1, "is_complex":1}) # preset already
						recdata = sp_filter.filt_table(recdata, bckgn)
						ipsi = tdir[ii]%100000
						iang = tdir[ii]/100000
						recdata.set_attr_dict({"bckgnoise":bckgn, "ctf":ct})
						if (current_group_ID  == group_ID):
							r.insert_slice(recdata, EMAN2_cppwrap.Transform({"type":"spider", "phi":refang[iang][0], "theta":refang[iang][1], \
								"psi":refang[iang][2]+ipsi*delta}), toprab*avgnorm/norm_per_particle[im])
						if (previous_group_ID == group_ID):
							r.insert_slice(recdata, EMAN2_cppwrap.Transform({"type":"spider", "phi":refang[iang][0], \
							  "theta":refang[iang][1], "psi":refang[iang][2]+ipsi*delta}), \
							     toprab*avgnorm/norm_per_particle[im]*(-1.))
	sp_utilities.reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
	sp_utilities.reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
	if myid == main_node:dummy = r.finish(True)
	mpi.mpi_barrier(mpi_comm)
	if myid == main_node: return fftvol, weight, refvol
	else:
		del fftvol
		del weight
		del refvol 
		return None, None, None
### insertion 2
def recons3d_trl_struct_group_nofsc_shifted_data_MPI(myid, main_node,\
       prjlist, paramstructure, norm_per_particle, group_ID,\
       mpi_comm = None, CTF = True, target_size=-1, nosmearing = False):
	"""
	  rec3d for pre-shifted data list
	  reconstructor nn4_ctfw
	"""
	if mpi_comm == None: mpi_comm = mpi.MPI_COMM_WORLD
	refvol = sp_utilities.model_blank(target_size)
	refvol.set_attr("fudge", 1.0)
	if CTF: do_ctf = 1
	else:   do_ctf = 0
	fftvol = EMAN2_cppwrap.EMData()
	weight = EMAN2_cppwrap.EMData()
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1,\
	    "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = EMAN2_cppwrap.Reconstructors.get( "nn4_ctfw", params)
	r.setup()
	if nosmearing:
		for im in range(len(prjlist)):
			if prjlist[im].get_attr("group") == group_ID:
				r.insert_slice(prjlist[im], prjlist[im].get_attr( "xform.projection" ), 1.0)
	else:
		nny = prjlist[0][0].get_ysize()
		rshifts_shrank = copy.deepcopy(Tracker["rshifts"])
		for im in range(len(rshifts_shrank)):
			rshifts_shrank[im][0] *= float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
			rshifts_shrank[im][1] *= float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
		refang =  Tracker["refang"]
		delta  =  Tracker["delta"]
		if not Tracker["constants"]["compute_on_the_fly"]: num_on_the_fly = len(prjlist)
		else: num_on_the_fly = Tracker["num_on_the_fly"][Blockdata["myid"]]
		for im in range(num_on_the_fly):
			if prjlist[im][0].get_attr("group") == group_ID:
				for jm in range(len(prjlist[im])):
					r.insert_slice(prjlist[im][jm], prjlist[im][jm].get_attr( "xform.projection" ), \
					  prjlist[im][jm].get_attr("wprob"))
		if (num_on_the_fly < len(prjlist)) and Tracker["constants"]["compute_on_the_fly"]:
			for im in range(num_on_the_fly, len(prjlist)):
				if prjlist[im][0].get_attr("group") == group_ID:
					avgnorm = Tracker["avgnorm"][prjlist[im][0].get_attr("chunk_id")]
					ct      = prjlist[im][0].get_attr("ctf")
					bckgn   = prjlist[im][0].get_attr("bckgnoise")
					numbor      = len(paramstructure[im][2])
					ipsiandiang = [ paramstructure[im][2][i][0]/1000  for i in range(numbor) ]
					allshifts   = [ paramstructure[im][2][i][0]%1000  for i in range(numbor) ]
					probs       = [ paramstructure[im][2][i][1] for i in range(numbor) ]
					tdir = list(set(ipsiandiang))
					data = [None]*len(Tracker["rshifts"])
					for ii in range(len(tdir)):
						#  Find the number of times given projection direction appears on the list, it is the number of different shifts associated with it.
						lshifts = sp_utilities.findall(tdir[ii], ipsiandiang)
						toprab  = 0.0
						for ki in range(len(lshifts)):toprab += probs[lshifts[ki]]
						recdata = EMAN2_cppwrap.EMData(nny,nny,1,False)
						recdata.set_attr("is_complex",0)
						for ki in range(len(lshifts)):
							lpt = allshifts[lshifts[ki]]
							if(data[lpt] == None):
								data[lpt] = sp_fundamentals.fshift(prjlist[im][0], rshifts_shrank[lpt][0], rshifts_shrank[lpt][1])
								data[lpt].set_attr("is_complex",0)
							EMAN2_cppwrap.Util.add_img(recdata, EMAN2_cppwrap.Util.mult_scalar(data[lpt], probs[lshifts[ki]]/toprab))
						recdata.set_attr_dict({"padffted":1, "is_fftpad":1,"is_fftodd":0, "is_complex_ri":1, "is_complex":1}) # preset already
						recdata = sp_filter.filt_table(recdata, bckgn)
						ipsi = tdir[ii]%100000
						iang = tdir[ii]/100000
						recdata.set_attr_dict({"bckgnoise":bckgn, "ctf":ct})
						r.insert_slice(recdata, EMAN2_cppwrap.Transform({"type":"spider", "phi":refang[iang][0], "theta":refang[iang][1], \
						  "psi":refang[iang][2]+ipsi*delta}), toprab*avgnorm/norm_per_particle[im])		
	sp_utilities.reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
	sp_utilities.reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
	if myid == main_node:dummy = r.finish(True)
	mpi.mpi_barrier(mpi_comm)
	if myid == main_node: return fftvol, weight, refvol
	else:
		del fftvol
		del weight
		del refvol 
		return None, None, None
###end of nofsc
def recons3d_trl_struct_group_MPI(myid, main_node, prjlist, random_subset, group_ID, \
      paramstructure, norm_per_particle = None, upweighted = True, mpi_comm= None, \
          CTF = True, target_size=-1, nosmearing = False):
	"""
		recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
		Input
			list_of_prjlist: list of lists of projections to be included in the reconstruction
	"""
	global Tracker, Blockdata
	if mpi_comm == None: mpi_comm = mpi.MPI_COMM_WORLD
	
	refvol = sp_utilities.model_blank(target_size)
	refvol.set_attr("fudge", 1.0)
	if CTF: do_ctf = 1
	else:   do_ctf = 0
	fftvol = EMAN2_cppwrap.EMData()
	weight = EMAN2_cppwrap.EMData()
	fftodd = 0
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, \
	  "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = EMAN2_cppwrap.Reconstructors.get( "nn4_ctfw", params )
	r.setup()
	if norm_per_particle == None: norm_per_particle = len(prjlist)*[1.0]
	if not nosmearing:
		delta  = Tracker["delta"]
		refang = Tracker["refang"]
		rshifts_shrank = copy.deepcopy(Tracker["rshifts"])
		for im in range(len(rshifts_shrank)):
			rshifts_shrank[im][0] *= float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
			rshifts_shrank[im][1] *= float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
		nny = prjlist[0].get_ysize()
	for im in range(len(prjlist)):
		#  parse projection structure, generate three lists:
		#  [ipsi+iang], [ishift], [probability]
		#  Number of orientations for a given image
		if prjlist[im].get_attr("group") == group_ID:
			if random_subset == 2:
				if nosmearing:
					r.insert_slice(prjlist[im], prjlist[im].get_attr("xform.projection"), 1.0)
				else:
					avgnorm = Tracker["avgnorm"][prjlist[im].get_attr("chunk_id")]
					numbor  = len(paramstructure[im][2])
					ipsiandiang = [ paramstructure[im][2][i][0]/1000  for i in range(numbor) ]
					allshifts   = [ paramstructure[im][2][i][0]%1000  for i in range(numbor) ]
					probs       = [ paramstructure[im][2][i][1] for i in range(numbor) ]
					#  Find unique projection directions
					tdir = list(set(ipsiandiang))
					bckgn = prjlist[im].get_attr("bckgnoise")
					ct = prjlist[im].get_attr("ctf")
					#  For each unique projection direction:
					data = [None]*len(Tracker["rshifts"])
					for ii in range(len(tdir)):
						#  Find the number of times given projection direction appears on the list, it is the number of different shifts associated with it.
						lshifts = sp_utilities.findall(tdir[ii], ipsiandiang)
						toprab  = 0.0
						for ki in range(len(lshifts)):  toprab += probs[lshifts[ki]]
						recdata = EMAN2_cppwrap.EMData(nny,nny,1,False)
						recdata.set_attr("is_complex",0)
						for ki in range(len(lshifts)):
							lpt = allshifts[lshifts[ki]]
							if( data[lpt] == None ):
								data[lpt] = sp_fundamentals.fshift(prjlist[im], rshifts_shrank[lpt][0], rshifts_shrank[lpt][1])
								data[lpt].set_attr("is_complex",0)
							EMAN2_cppwrap.Util.add_img(recdata, EMAN2_cppwrap.Util.mult_scalar(data[lpt], probs[lshifts[ki]]/toprab))
						recdata.set_attr_dict({"padffted":1, "is_fftpad":1,"is_fftodd": fftodd, "is_complex_ri":1, "is_complex":1})
						if not upweighted:  recdata = sp_filter.filt_table(recdata, bckgn )
						recdata.set_attr_dict( {"bckgnoise":bckgn, "ctf":ct} )
						ipsi = tdir[ii]%100000
						iang = tdir[ii]/100000
						r.insert_slice( recdata, \
						   EMAN2_cppwrap.Transform({"type":"spider","phi":refang[iang][0],"theta":refang[iang][1],"psi":refang[iang][2]+ipsi*delta}), \
						     toprab*avgnorm/norm_per_particle[im])
			else:
				if	prjlist[im].get_attr("chunk_id") == random_subset:
					if nosmearing:
						r.insert_slice(prjlist[im], prjlist[im].get_attr("xform.projection"), 1.0)
					else:
						if Tracker["constants"]["nsmear"]<=0.0: numbor = len(paramstructure[im][2])
						else:  numbor = 1
						ipsiandiang = [ paramstructure[im][2][i][0]/1000  for i in range(numbor) ]
						allshifts   = [ paramstructure[im][2][i][0]%1000  for i in range(numbor) ]
						probs       = [ paramstructure[im][2][i][1] for i in range(numbor) ]
						#  Find unique projection directions
						tdir  = list(set(ipsiandiang))
						bckgn = prjlist[im].get_attr("bckgnoise")
						ct    = prjlist[im].get_attr("ctf")
						#  For each unique projection direction:
						data = [None]*len(rshifts_shrank)
						for ii in range(len(tdir)):
							#  Find the number of times given projection direction appears on the list, it is the number of different shifts associated with it.
							lshifts = sp_utilities.findall(tdir[ii], ipsiandiang)
							toprab  = 0.0
							for ki in range(len(lshifts)):  toprab += probs[lshifts[ki]]
							recdata = EMAN2_cppwrap.EMData(nny,nny,1,False)
							recdata.set_attr("is_complex",0)
							for ki in range(len(lshifts)):
								lpt = allshifts[lshifts[ki]]
								if( data[lpt] == None ):
									data[lpt] = sp_fundamentals.fshift(prjlist[im], rshifts_shrank[lpt][0], rshifts_shrank[lpt][1])
									data[lpt].set_attr("is_complex",0)
								EMAN2_cppwrap.Util.add_img(recdata, EMAN2_cppwrap.Util.mult_scalar(data[lpt], probs[lshifts[ki]]/toprab))
							recdata.set_attr_dict({"padffted":1, "is_fftpad":1,"is_fftodd":fftodd, "is_complex_ri":1, "is_complex":1})
							if not upweighted:  recdata = sp_filter.filt_table(recdata, bckgn )
							recdata.set_attr_dict( {"bckgnoise":bckgn, "ctf":ct} )
							ipsi = tdir[ii]%100000
							iang = tdir[ii]/100000
							r.insert_slice(recdata, \
							 EMAN2_cppwrap.Transform({"type":"spider","phi":refang[iang][0],"theta":refang[iang][1],"psi":refang[iang][2]+ipsi*delta}), \
							  toprab*avgnorm/norm_per_particle[im])
	#  clean stuff
	#if not nosmearing: del recdata, tdir, ipsiandiang, allshifts, probs
	sp_utilities.reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
	sp_utilities.reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
	if not nosmearing: del rshifts_shrank
	if myid == main_node:dummy = r.finish(True)
	mpi.mpi_barrier(mpi_comm)
	if myid == main_node: return fftvol, weight, refvol
	else:
		del fftvol
		del weight
		del refvol 
		return None, None, None
###===============>>>  group rec3d  <<<===========================================================
### final rec3d
def do3d_sorting_groups_nofsc_final(rdata, parameterstructure, norm_per_particle):
	global Tracker, Blockdata
	keepgoing = 1
	if(Blockdata["myid"] == Blockdata["last_node"]):
		if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
			os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
		try:
			with open(os.path.join(Tracker["directory"],"freq_cutoff.json"),'r') as fout:
				freq_cutoff_dict = sp_utilities.convert_json_fromunicode(json.load(fout))
			fout.close()
		except: freq_cutoff_dict = 0
	else: freq_cutoff_dict = 0
	freq_cutoff_dict = sp_utilities.wrap_mpi_bcast(freq_cutoff_dict, Blockdata["last_node"], mpi.MPI_COMM_WORLD)
	for index_of_groups in range(Tracker["number_of_groups"]):
		tvol, tweight, trol = recons3d_trl_struct_group_MPI(Blockdata["myid"], \
		  Blockdata["last_node"], rdata, 2, index_of_groups, parameterstructure, norm_per_particle, \
      True, None, Tracker["constants"]["CTF"], (2*Tracker["nxinit"]+3), Tracker["nosmearing"])
		if(Blockdata["myid"] == Blockdata["last_node"]):
			tvol.set_attr("is_complex",0)
			tvol.write_image(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf"%index_of_groups))
			tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf"%index_of_groups))
			trol.write_image(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_groups))
			del tvol
			del tweight
			del trol
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	Tracker["maxfrad"] = Tracker["nxinit"]//2
	if Blockdata["no_of_groups"]>1:
	 	# new starts
		sub_main_node_list = [ -1 for i in range(Blockdata["no_of_groups"])]
		for index_of_colors in range(Blockdata["no_of_groups"]):
			for iproc in range(Blockdata["nproc"]-1):
				if Blockdata["myid"]== iproc:
					if Blockdata["color"] == index_of_colors and Blockdata["myid_on_node"] == 0:
						sub_main_node_list[index_of_colors] = Blockdata["myid"]
					sp_utilities.wrap_mpi_send(sub_main_node_list, Blockdata["last_node"], mpi.MPI_COMM_WORLD)
				if Blockdata["myid"] == Blockdata["last_node"]:
					dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
					for im in range(len(dummy)):
						if dummy[im]>-1: sub_main_node_list[im] = dummy[im]
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		sp_utilities.wrap_mpi_bcast(sub_main_node_list, Blockdata["last_node"], mpi.MPI_COMM_WORLD)
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
		for iloop in range(nbig_loop):
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
			
				if(Blockdata["myid"] == Blockdata["last_node"]):
					tvol2    = sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")%index_of_group)
					tweight2 = sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")%index_of_group)
					treg2 	 = sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_group))
					tvol2.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1})
					tag = 7007
					sp_utilities.send_EMData(tvol2,    sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
					sp_utilities.send_EMData(tweight2, sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
					sp_utilities.send_EMData(treg2,    sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
				
				elif (Blockdata["myid"] == sub_main_node_list[index_of_colors]):
					tag      = 7007
					tvol2    = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
					tweight2 = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
					treg2    = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
				mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				try: 
					Tracker["freq_fsc143_cutoff"] = freq_cutoff_dict["Cluster_%03d.txt"%index_of_group]
				except: pass	
				if Blockdata["color"] == index_of_colors:  # It has to be 1 to avoid problem with tvol1 not closed on the disk 							
					if(Blockdata["myid_on_node"] != 0): 
						tvol2     = sp_utilities.model_blank(1)
						tweight2  = sp_utilities.model_blank(1)
						treg2     = sp_utilities.model_blank(1)
					tvol2 = steptwo_mpi_reg(tvol2, tweight2, treg2, None, \
					 Tracker["freq_fsc143_cutoff"], 0.01, False, color = index_of_colors) # has to be False!!!
					del tweight2, treg2
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				if (Blockdata["color"] == index_of_colors) and (Blockdata["myid_on_node"] == 0):
					tag = 7007
					sp_utilities.send_EMData(tvol2, Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
				elif(Blockdata["myid"] == Blockdata["last_node"]):
					tag = 7007
					tvol2    = sp_utilities.recv_EMData(sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
					if Tracker["Core_set"]:
						if index_of_group !=Tracker["number_of_groups"]-1:
							tvol2.write_image(os.path.join(Tracker["directory"], "vol_cluster%03d.hdf"%index_of_group))
						else:
							tvol2.write_image(os.path.join(Tracker["directory"], "volume_core.hdf"))
					else: tvol2.write_image(os.path.join(Tracker["directory"], "vol_cluster%03d.hdf"%index_of_group))
					del tvol2
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
		
	else:# loop over all groups for single node
		for index_of_group in range(Tracker["number_of_groups"]):
			if(Blockdata["myid_on_node"] == 0):
				tvol2 		= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")%index_of_group)
				tweight2 	= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")%index_of_group)
				treg2 		= sp_utilities.get_im(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_group))
				tvol2.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1})
			else:
				tvol2 		= sp_utilities.model_blank(1)
				tweight2 	= sp_utilities.model_blank(1)
				treg2		= sp_utilities.model_blank(1)
			try: 
				Tracker["freq_fsc143_cutoff"] = freq_cutoff_dict["Cluster_%03d.txt"%index_of_group]
			except: pass
			tvol2 = steptwo_mpi_reg(tvol2, tweight2, treg2, None, Tracker["freq_fsc143_cutoff"], 0.01, False) # has to be False!!!
			del tweight2, treg2
			if(Blockdata["myid_on_node"] == 0):
				if Tracker["Core_set"]:
					if index_of_group !=Tracker["number_of_groups"]-1:
						tvol2.write_image(os.path.join(Tracker["directory"], "vol_cluster%03d.hdf"%index_of_group))
					else:
						tvol2.write_image(os.path.join(Tracker["directory"], "volume_core.hdf"))
				else: tvol2.write_image(os.path.join(Tracker["directory"], "vol_cluster%03d.hdf"%index_of_group))
				del tvol2
			mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	keepgoing = sp_utilities.bcast_number_to_all(keepgoing, source_node = Blockdata["main_node"], mpi_comm = mpi.MPI_COMM_WORLD) # always check 
	Tracker   = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"])
	if not keepgoing: sp_global_def.ERROR("do3d_sorting_groups_nofsc_final  %s"%os.path.join(Tracker["directory"], \
	     "tempdir"),"do3d_sorting_groups_nofsc_final", 1, Blockdata["myid"])
	return
#####=============================================================================
####------->>> UNIX MEM related functions <<<--------------
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
def sorting_main_mpi(log_main, not_include_unaccounted):
	global Tracker, Blockdata
	time_sorting_start = time.time()
	read_tracker_mpi(Tracker["constants"]["masterdir"])
	depth_order = Tracker["depth_order"]
	Tracker["generation"]         =  {}
	Tracker["current_generation"] =  0
	keepsorting                   =  1
	keepchecking                  =  1
	Tracker["current_generation"] = -1
	igen              = -1
	my_pids           = os.path.join(Tracker["constants"]["masterdir"], "indexes.txt")
	params            = os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt")
	previous_params   = Tracker["previous_parstack"]
	all_gen_stat_list = []
	bad_clustering    = 0
	while (keepsorting == 1) and (bad_clustering== 0): ### generation 0 sort out the groups with the most strong features 
		Tracker["current_generation"] +=1
		igen +=1
		work_dir = os.path.join(Tracker["constants"]["masterdir"], "generation_%03d"%igen)
		if Blockdata["myid"] == Blockdata["main_node"]:
			set_minimum_group_size(log_main, printing_dres_table = False)
			keepchecking = check_sorting_state(work_dir, keepchecking, log_main)
		else: 
			keepchecking = 0
			Tracker      = 0
		keepchecking = sp_utilities.bcast_number_to_all(keepchecking, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		Tracker      = sp_utilities.wrap_mpi_bcast(Tracker,           Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		if keepchecking == 0: # new, do it
			if Blockdata["myid"] == Blockdata["main_node"]:
				time_generation_start = time.time()
				if not os.path.exists(work_dir):
					os.mkdir(work_dir)# need check each box
					within_generation_restart   = 0
				else: within_generation_restart = 1
				freq_cutoff_dict = {}
				with open(os.path.join(work_dir, "freq_cutoff.json"),'w') as fout:
					json.dump(freq_cutoff_dict, fout)
				fout.close()
				log_main.add('================================================================================================================')
				log_main.add('                                    SORT3D MULTI-LAYER   generation %d'%igen)
				log_main.add('----------------------------------------------------------------------------------------------------------------')
			else: within_generation_restart = 0
			within_generation_restart = sp_utilities.bcast_number_to_all(within_generation_restart, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			if within_generation_restart == 1: read_tracker_mpi(work_dir)
			else: dump_tracker(work_dir)
			output_list, bad_clustering, stat_list = depth_clustering(work_dir, depth_order, \
			    my_pids, params, previous_params, log_main)
			all_gen_stat_list.append(stat_list)
			
			if bad_clustering !=1:
				if Blockdata["myid"] == Blockdata["main_node"]:
					clusters, nclusters, nuacc  = output_clusters(work_dir, output_list[0][0], \
					    output_list[0][1], not_include_unaccounted, log_main)
					try: del Tracker["generation"][str(igen)]
					except: pass
					 
					if not Tracker["no_cluster_generation"]: 
						Tracker["generation"][str(igen)] = len(clusters)
					else:  Tracker["generation"][str(igen)]  = 0
				else: 
					Tracker   = 0
					nclusters = 0
					nuacc     = 0
					
				Tracker   = sp_utilities.wrap_mpi_bcast(Tracker,        Blockdata["main_node"], mpi.MPI_COMM_WORLD)
				nclusters = sp_utilities.bcast_number_to_all(nclusters, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
				nuacc     = sp_utilities.bcast_number_to_all(nuacc,     Blockdata["main_node"], mpi.MPI_COMM_WORLD)
				dump_tracker(work_dir)
				if nclusters == 0:
					if Blockdata["myid"] == Blockdata["main_node"]:
						log_main.add('No cluster is found in generation %d'%igen)
					break
				else:
					if Blockdata["myid"] == Blockdata["main_node"]:
						time_of_generation_h,  time_of_generation_m = get_time(time_generation_start)
						log_main.add('SORT3D generation%d time: %d hours %d minutes.'%(igen, \
						   time_of_generation_h, time_of_generation_m))
					my_pids = os.path.join( work_dir, 'indexes_next_generation.txt')
					if Blockdata["myid"] == Blockdata["main_node"]:
						sp_utilities.write_text_file(output_list[0][1], my_pids)
						sp_utilities.write_text_row(stat_list, os.path.join(work_dir, 'gen_rep.txt'))
						mark_sorting_state(work_dir, True)
					mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
				keepsorting = check_sorting(nuacc, keepsorting, log_main)
				if keepsorting != 1: break #
		else: # restart run
			read_tracker_mpi(work_dir)
			my_pids  = os.path.join( work_dir, 'indexes_next_generation.txt') 
			if Blockdata["myid"] == Blockdata["main_node"]: 
				stat_list = sp_utilities.read_text_row( os.path.join(work_dir, 'gen_rep.txt'))
			else: stat_list = 0
			stat_list = sp_utilities.wrap_mpi_bcast(stat_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			all_gen_stat_list.append(stat_list)
			if Blockdata["myid"] == Blockdata["main_node"]:			
				log_main.add('================================================================================================================')
				log_main.add('                                    SORT3D MULTI-LAYER   generation %d completed'%igen)
				log_main.add('----------------------------------------------------------------------------------------------------------------')
	copy_results(log_main, all_gen_stat_list)# all nodes function
	compute_final_map(Tracker["constants"]["masterdir"], log_main)
	if Blockdata["myid"] == Blockdata["main_node"]:
		time_of_sorting_h,  time_of_sorting_m = get_time(time_sorting_start)
		log_main.add('SORT3D execution time: %d hours %d minutes.'%(time_of_sorting_h, time_of_sorting_m))
	return
###---------------------------------------------------------------------------------------
def check_mpi_settings(log_main):
	global Tracker, Blockdata
	####-------------------->>>>>Parameters <<<------------
	largest_K = 30
	#####---------------------------------------------------------------------------------
	if not Tracker["search_mode"]:
		number_of_groups = Tracker["constants"]["total_stack"]//Tracker["constants"]["img_per_grp"]
	else:
		number_of_groups =  Tracker["constants"]["init_K"]
	if number_of_groups>=largest_K:
		sp_global_def.ERROR("User-provided img_per_grp is too small", "check_mpi_settings", 1, Blockdata["myid"])	
	avg_fsc = sp_statistics.scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], \
	   float(Tracker["constants"]["total_stack"]), \
	   float(Tracker["constants"]["total_stack"])//number_of_groups)
	fsc143 = get_res143(avg_fsc)
	if fsc143 !=0: # Empirical image size estimation
		nxinit = sp_fundamentals.smallprime(min((int(fsc143)+ max(int(Tracker["constants"]["nnxo"]*0.03), 5))*2, \
		     Tracker["constants"]["nnxo"]))
		nxinit = sp_fundamentals.smallprime(nxinit-nxinit%2)
		nxinit = nxinit - nxinit%2
	else: sp_global_def.ERROR("Incorrect FSC0.143", "check_mpi_settings", 1, Blockdata["myid"])
	### pre-calculated 2D data: ----------------------------------------------------------
	#1. Cdata(data for comparison)
	#2. Rdata(data for rec3D) # smearing could dramatically increase data size
	#3. fdata (focus projection); 
	#4. 2D CTFs 
	# other large data
	#1. reference volumes
	#2. rec3D volumes 2 times padded  
	##------------------------------------------------------------------------------------
	sys_required_mem = Tracker["constants"]["overhead"]
	if( Blockdata["myid"] == Blockdata["main_node"]):
		log_main.add('\n')
		log_main.add('----------------------------------------------------------------------------------------------------------------' )
		log_main.add('                 =======>     Number of input images and memory information     <=====')
		log_main.add('----------------------------------------------------------------------------------------------------------------' )
		log_main.add("Number of processes:                                                %8d"%Blockdata["nproc"])
		log_main.add("Number of nodes:                                                    %8d"%Blockdata["no_of_groups"])
		log_main.add("Number of processes per node:                                       %8d"%Blockdata["no_of_processes_per_group"])
		
	try:
		image_org_size             = Tracker["constants"]["nnxo"]
		image_in_core_size         = nxinit
		ratio                      = float(nxinit)/float(image_org_size)
		raw_data_size              = float(Tracker["constants"]["total_stack"]*image_org_size*image_org_size)*4.0/1.e9
		raw_data_size_per_node     = float(Tracker["constants"]["total_stack"]*image_org_size*image_org_size)*4.0/1.e9/Blockdata["no_of_groups"]
		data_for_comparison        = raw_data_size_per_node*ratio**2
		data_for_rec3d             = raw_data_size_per_node*ratio**2
		data_ctf                   = raw_data_size_per_node*ratio**2
		if Tracker["constants"]["focus3D"]:focusdata = raw_data_size_per_node*ratio**2
		else:                              focusdata = 0.0
		sorting_data_size_per_node = data_for_comparison + data_for_rec3d + data_ctf + focusdata
		volume_size_per_node = (4.*image_in_core_size**3*8.)*Blockdata["no_of_processes_per_group"]/1.e9
		sorting_data_size_per_node += volume_size_per_node # memory peak 
	except:
		sp_global_def.ERROR("Initial info is not provided", "check_mpi_settings", 1, Blockdata["myid"])
	try:
		mem_bytes = os.sysconf('SC_PAGE_SIZE')*os.sysconf('SC_PHYS_PAGES')# e.g. 4015976448
		mem_gib   = mem_bytes/(1024.**3) # e.g. 3.74
		if( Blockdata["myid"] == Blockdata["main_node"]):
			log_main.add("Available memory information provided by the operating system(GB):  %8.3f "%mem_gib)
	except: mem_gib = None
	
	if Tracker["constants"]["memory_per_node"] == -1.:
		if mem_gib: total_memory = mem_gib
		else:
			total_memory =  Blockdata["no_of_processes_per_group"]*2.0 # assume each CPU has 2.0 G
			if( Blockdata["myid"] == Blockdata["main_node"]):
				log_main.add("Memory per node is not provided, sort3d assumes 2GB per node")
		Tracker["constants"]["memory_per_node"] = total_memory
	else:
		total_memory =  Tracker["constants"]["memory_per_node"]
		if( Blockdata["myid"] == Blockdata["main_node"]):
			log_main.add("Memory per node:                                            %f"%Tracker["constants"]["memory_per_node"])
			
	if(Blockdata["myid"] == Blockdata["main_node"]):
		log_main.add("Total number of images:                                             %8d     "%Tracker["constants"]["total_stack"])
		if not Tracker["search_mode"]:  
			log_main.add("Number of images per group:                                         %8d     "%Tracker["constants"]["img_per_grp"])
		else:
			log_main.add("Number of images per group:                                         %8d     "%(Tracker["constants"]["total_stack"]//Tracker["constants"]["init_K"]))
		log_main.add("The total available memory per node(in GB):                         %8.3f   "%total_memory)
		log_main.add("The size of total input 2D stack(in GB):                            %8.3f   "%(raw_data_size))
		log_main.add("The per-node amount of memory 2D data will occupy(in GB):           %8.3f   "%(raw_data_size_per_node))
		log_main.add("Precalculated 2D data (no smearing) for sorting(in GB)              %8.3f   "%(sorting_data_size_per_node))
		log_main.add("Reserved memory for overhead loading per node(in GB):               %8.3f   "%sys_required_mem)
		
	if (total_memory - sys_required_mem - sorting_data_size_per_node ) <0.0:
		sp_global_def.ERROR("Insufficient memory to pass the sorting processs. Increase a node", "check_mpi_settings", 1, Blockdata["myid"])
		
	if (total_memory - sys_required_mem - raw_data_size_per_node - 2*raw_data_size_per_node*ratio**2 ) <0.0: 
		sp_global_def.ERROR("Insufficient memory to read in the raw data to produce downsized data. Increase a node", \
		   "check_mpi_settings", 1, Blockdata["myid"])
			
	images_per_cpu = int(float(Tracker["constants"]["total_stack"])/Blockdata["nproc"])
	if not Tracker["search_mode"]:
		images_per_cpu_for_unaccounted_data  = Tracker["constants"]["img_per_grp"]*1.5/float(Blockdata["nproc"])
	else:
		images_per_cpu_for_unaccounted_data  = \
		   Tracker["constants"]["total_stack"]//Tracker["constants"]["init_K"]*1.5/float(Blockdata["nproc"])
	if( Blockdata["myid"] == Blockdata["main_node"]):
		log_main.add("Number of images per processor:                                     %8d "%images_per_cpu)
		
	if( images_per_cpu < 5):
		sp_global_def.ERROR("Number of images per processor is less than 5", \
		   "one may want to consider decreasing the number of processors used in MPI setting", \
		   0, Blockdata["myid"])
		   
	if(Blockdata["myid"] == Blockdata["main_node"]):
		avg_smearing_on_node = sum(sp_utilities.read_text_file(Tracker["constants"]["smearing_file"]))//Blockdata["no_of_groups"]
	else: avg_smearing_on_node = 0
	avg_smearing_on_node  = sp_utilities.bcast_number_to_all(avg_smearing_on_node, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	avg_images_on_node    = Tracker["constants"]["total_stack"]//Blockdata["no_of_groups"]
	precalculated_2D_data = raw_data_size_per_node*ratio**2*(max(avg_smearing_on_node/avg_images_on_node,1.) + 3.)
	
	if ((precalculated_2D_data - sys_required_mem)>total_memory) and (not Tracker["constants"]["compute_on_the_fly"]):
		sp_global_def.ERROR("The size of precalculated shifted 2D data is too large. Turn on the option compute_on_the_fly", \
		  "check_mpi_settings", 1, Blockdata["myid"])
	else:
		if(Blockdata["myid"] == Blockdata["main_node"]):
			log_main.add("Percentage of memory occupied by precalculated shifted 2D data:     %8.3f%% "%\
			   (precalculated_2D_data/total_memory*100.))
			if Tracker["constants"]["compute_on_the_fly"]:
				log_main.add("The compute_on_the_fly option is on. The multiplely assiged 2D data is precacluated till 90 percents of the memory is filled up")
	if(Blockdata["myid"] == Blockdata["main_node"]):
		log_main.add('----------------------------------------------------------------------------------------------------------------\n')
	return
#####-------------------------------------------------------------------------------------
def set_sorting_global_variables_mpi(log_file):
	global Tracker, Blockdata 
	if Blockdata["myid"] == Blockdata["main_node"]:
		###=====<--options for advanced users:
		log_file.add("---------------------------------------------------------------------------------------------------")
		log_file.add(" Global adjustable varibales (to be modified for advanced users to fit specific clustering requirements)")
		log_file.add("---------------------------------------------------------------------------------------------------")
		Tracker["total_number_of_iterations"]         = 15
		Tracker["clean_volumes"]                      = True  # always true
		Tracker["swap_ratio"]                         = 5.
		Tracker["stop_mgskmeans_percentage"]          = 10.
		Tracker["depth_order"]                        = 2
		Tracker["random_group_elimination_threshold"] = 0.
		if (Tracker["freeze_groups"] ==1): 
			Tracker["nstep"]         = 1
			Tracker["box_learning"]  = 0
		### -----------Orientation constraints--------------------------------------------
		Tracker["tilt1"]                      =  0.0
		Tracker["tilt2"]                      =  180.0
		Tracker["minimum_ptl_number"]         =  20
        ### ------------------------------------------------------------------------------
		Tracker["minimum_img_per_cpu"]        =  5 # User can adjust it to other number
		Tracker["minimum_grp_size"]           = [0, 0, 0, 0] # low, high, res, total
		Tracker["num_on_the_fly"]             = [ -1  for im in range(Blockdata["nproc"])]
		if not Tracker["search_mode"]:          Tracker["current_img_per_grp"] = Tracker["constants"]["img_per_grp"]
		else:  Tracker["current_img_per_grp"] = Tracker["constants"]["total_stack"]//Tracker["constants"]["init_K"]
		Tracker["constants"]["box_niter"]     = 5
		Tracker["even_scale_fsc"]             = True
		Tracker["rand_index"]                 = 0.0
		Tracker["nruns"]                      = -1
		Tracker["minimum_grp_size"]           = [0, 0, 0, 0]
		Tracker["img_per_grp"]                = [0, 0, 0, 0]
		Tracker["do_timing"]                  = False
		Tracker["no_cluster_generation"]      = False
		Tracker["fixed_sorting_size"]         = False
		### ---------------------------------------------------------------------------------------------
		if Tracker["search_mode"]:
			log_file.add("Sorting is set on search_mode and starts from K = %d"%Tracker["constants"]["init_K"])
			
		elif (Tracker["freeze_groups"]== 1) and (not Tracker["search_mode"]):
			log_file.add("Sorting is set on freeze_groups mode (prior to v1.2)")
			
		else: log_file.add("Sorting is set on general mode")
		
		if Tracker["use_umat"]: log_file.add("Fuzzy membership is applied to stabilize sorting")
			
		log_file.add("Sorting depth order:                                   %8d"%Tracker["depth_order"])
		log_file.add("Maximum box iteration:                                 %8d"%Tracker["constants"]["box_niter"])
		log_file.add("Minimum_img_per_cpu:                                   %8d"%Tracker["minimum_img_per_cpu"])
		log_file.add("MGSKmeans maximum iteration:                           %8d"%Tracker["total_number_of_iterations"])
		log_file.add("Minimum number of images per orientation angle:        %8d"%Tracker["minimum_ptl_number"])
		log_file.add("Minimum group shrink size steps:                       %8d"%Tracker["nstep"])
		log_file.add("Box_learning:                                          %8d"%Tracker["box_learning"])
		log_file.add("Ratio for swapping accounted for and unaccounted:      %8.2f"%Tracker["swap_ratio"])
		log_file.add("Freeze small clusters during within box comparison:    %8d"%Tracker["freeze_groups"])
		log_file.add("Particles change perecent for stopping MGSKmeans:      %8.2f"%Tracker["stop_mgskmeans_percentage"])
		log_file.add("Random group elimination threshold:                    %8.2f"%Tracker["random_group_elimination_threshold"])
		log_file.add("Num_core_set:                                          %8d"%Tracker["constants"]["num_core_set"])
			
		###---------------------------------------------------------------------------------------------------
		calculate_grp_size_variation(log_file, state = "sort_init")
		####------------------------------------------------------------------------------
		log_file.add("-----------------------------------------------------------------------------------")
		log_file.add(" Compute scaling relationship (dres_table) between FSC0.143 and image group size   ")
		log_file.add(" FSC0.143 (in pixels)  size  low bound     high bound    width of pixel bin        ")
		log_file.add("-----------------------------------------------------------------------------------")
		dict         = {}
		previous_res = 0
		res_list     = []
		growing_size = Tracker["constants"]["img_per_grp"]//100
		step = max(growing_size//1000, 20)
		while (growing_size + step<Tracker["total_stack"]):
			growing_size +=step
			res = get_res143(sp_statistics.scale_fsc_datasetsize(\
			   Tracker["constants"]["fsc_curve"], Tracker["constants"]["total_stack"], growing_size))
			if res == previous_res: dict[str(res)][1] = growing_size
			else: 
				res_list.append(str(res))
				dict[str(res)] = [growing_size, growing_size]
			previous_res       = res
		for il in range(len(res_list)):
			log_file.add("%5s               %10d           %10d        %10d"%(res_list[il], dict[res_list[il]][0], \
			   dict[res_list[il]][1],  (dict[res_list[il]][1]- dict[res_list[il]][0])))
		log_file.add("-----------------------------------------------------------------------------------")
		###-----------------------------------------------------------
		set_minimum_group_size(log_file, printing_dres_table = True)
		###-----------------------------------------------------------
	else: Tracker = 0
	Tracker = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	return
###---------------------------------------------------------------------------------------
def mem_calc_and_output_info(smearing_file, nxinit, iter_id_init_file, log_main):
	global Blockdata, Tracker
	# single CPU ; no info update
	log_main.add( "-------------------------------------------------------------")
	log_main.add( "              Precalculated images info                      ")
	log_main.add( "-------------------------------------------------------------")
	
	precalculated_images_per_cpu = np.full(Blockdata["nproc"], 0, dtype=np.int32)
	smearing_list = np.array(sp_utilities.read_text_file(smearing_file),       dtype=np.int32)
	indx_list     = sp_utilities.read_text_file(iter_id_init_file, -1)
	if len(indx_list) ==1: indx_list= indx_list[0]
	else:                  indx_list= indx_list[1]
	indx_list = np.sort(np.array(indx_list, dtype = np.int32))
	if (smearing_list.shape[0] !=indx_list.shape[0]): smearing_list = smearing_list[indx_list]
	avg_smear = np.sum(smearing_list)/smearing_list.shape[0]
	cdata_in_core  = (Tracker["total_stack"]*nxinit*nxinit*4.0)/1.e9/Blockdata["no_of_groups"]
	srdata_in_core = (nxinit*nxinit*np.sum(smearing_list)*4.)/1.e9/Blockdata["no_of_groups"]
	
	if (not Tracker["constants"]["focus3D"]): fdata_in_core = 0.0
	else: fdata_in_core = cdata_in_core
	ctfdata = cdata_in_core
	refvol_size = (nxinit*nxinit*nxinit*4.0*2)/1.e9*Blockdata["no_of_processes_per_group"]# including the 3D mask
	log_main.add( "Precalculated data (GB) in core per node (available memory per node: %6.2f):"%Tracker["constants"]["memory_per_node"])
	log_main.add( "Images for comparison: %6.2f GB; shifted images: %6.2f GB; focus images: %6.2f GB; ctfs: %6.2f GB"%\
			(cdata_in_core, srdata_in_core, fdata_in_core, ctfdata))
	tdata = cdata_in_core+srdata_in_core+ctfdata+refvol_size+fdata_in_core
	
	if (tdata/Tracker["constants"]["memory_per_node"]*100.> 90.):
		if (not Tracker["constants"]["compute_on_the_fly"]):
			sp_global_def.ERROR("More than 90% memory is used. Turn on compute_on_the_fly and rerun the program", \
			   "mem_calc_and_output_info", 1, Blockdata["myid"])
			   
	log_main.add("The data to be in core for sorting occupies %7.3f percents of memory;  average smearing is %5.1f"%\
		(min(tdata/Tracker["constants"]["memory_per_node"]*100., 90.0), avg_smear))

	log_main.add("Estimated required memory for sorting on each node:\n")
	smearings_on_nodes = np.full(Blockdata["no_of_groups"], 0.0, dtype=np.float32)
	smearings_per_cpu  = [None for im in range(Blockdata["nproc"])]
	
	dict = [None for i in range(Blockdata["nproc"])]
	for iproc in range(Blockdata["nproc"]):
		image_start, image_end = sp_applications.MPI_start_end(smearing_list.shape[0], Blockdata["nproc"],iproc)
		smearings_on_nodes[iproc//Blockdata["no_of_processes_per_group"]] += \
			np.sum(smearing_list[image_start:image_end])*(nxinit*nxinit*4.)/1.e9
		smearings_per_cpu[iproc] = smearing_list[image_start:image_end]
	### output info ===================================
	msg = ""
	for icolor in range(Blockdata["no_of_groups"]):
		tdata = cdata_in_core + ctfdata + refvol_size + smearings_on_nodes[icolor]+fdata_in_core
		if Tracker["constants"]["focus3D"]: tdata +=fdata_in_core
		msg += "Node %5d :  Mem %7.2f GB  "%(icolor, min(tdata, Tracker["constants"]["memory_per_node"]*0.9))+"; "
		if (icolor%3==2):
			log_main.add(msg)
			msg =""
	if Blockdata["no_of_groups"]%3!=0:log_main.add(msg)
	# Estimate number of precomputed images
	overhead_mem = Tracker["constants"]["overhead"] #
	mem_leftover = max(Tracker["constants"]["memory_per_node"], 32.) - \
	  overhead_mem - 2.*cdata_in_core - fdata_in_core - refvol_size
	  
	mem_leftover = mem_leftover/float(Blockdata["no_of_processes_per_group"])
	log_main.add("Number of particles precomputed on each CPU:\n")
	log_main.add("Node       Number   Pecentage")
	msg = ""
	for im in range(Blockdata["nproc"]):
		image_start, image_end = sp_applications.MPI_start_end(smearing_list.shape[0], Blockdata["nproc"],im)
		size       = image_end - image_start 
		mem_on_cpu = mem_leftover
		jm = 0
		while (mem_on_cpu>0.0) and (jm<size):
			mem_on_cpu -=smearings_per_cpu[im][jm]*(nxinit*nxinit*4.)/1.e9
			jm +=1
		dict[im] = jm
		msg += "%5d   %8d     %6.1f%% "%(im, jm, jm/float(size)*100.)+"; "
		if (im%3==2):
			log_main.add(msg)
			msg =""
	if Blockdata["nproc"]%3!=0:log_main.add(msg+"\n")
	log_main.add( "-------------------------------------------------------------")
	return dict
##-------------------------------------- Final -------------------------------------------
def copy_results(log_file, all_gen_stat_list):
	global Tracker, Blockdata
	nsgen = 0
	for i in range(len(all_gen_stat_list)): nsgen += len(all_gen_stat_list[i])
	if nsgen >0:	
		if Blockdata["myid"] == Blockdata["main_node"]:
			log_file.add('================================================================================================================')
			log_file.add('                     Final results saved in %s'%Tracker["constants"]["masterdir"])
			log_file.add('----------------------------------------------------------------------------------------------------------------' )
			nclusters = 0
			log_file.add( '{:^8} {:>8}   {:^24}  {:>15} {:^22} {:^5} {:^20} {:^20} '.format('Group ID', '    size','determined in generation', 'reproducibility', 'random reproducibility', ' std ', ' selection file', '       map file     '))
			clusters  = []
			NACC      = 0           
			for ig1, value in list(Tracker["generation"].items()):
				ig = string.atoi('%s'%ig1)
				for ic in range(value):
					cluster_file = os.path.join(Tracker["constants"]["masterdir"], "generation_%03d"%ig, "Cluster_%03d.txt"%ic)
					shutil.copyfile(cluster_file, os.path.join(Tracker["constants"]["masterdir"], "Cluster_%03d.txt"%nclusters))
					clusters.append(sp_utilities.read_text_file(cluster_file))
					cluster      = sp_utilities.read_text_file(os.path.join(Tracker["constants"]["masterdir"], "generation_%03d"%ig, "Cluster_%03d.txt"%ic))
					cluster_file = "Cluster_%03d.txt"%nclusters
					vol_file     = "vol_cluster%03d.hdf"%nclusters
					msg          = '{:>8} {:>8}   {:^24}        {:^6}          {:^6}          {:>5} {:^20} {:^20} '.format(nclusters, \
					  len(cluster), ig, round(all_gen_stat_list[ig][ic][0], 1), round(all_gen_stat_list[ig][ic][1], 1), \
					    round(all_gen_stat_list[ig][ic][2],1), cluster_file,  vol_file)
					nclusters   +=1
					NACC        +=len(cluster)
					log_file.add(msg)
					
			if nclusters>=1:			
				try:
					Unaccounted_file = os.path.join(Tracker["constants"]["masterdir"], \
					   "generation_%03d"%Tracker["current_generation"], "Core_set.txt")
					cluster = sp_utilities.read_text_file(Unaccounted_file)
					shutil.copyfile(Unaccounted_file, os.path.join(Tracker["constants"]["masterdir"], "Core_set.txt"))
					cluster_file = "Core_set.txt"
					if (Tracker["constants"]["num_core_set"]>0) and (len(cluster) > Tracker["constants"]["num_core_set"]): 
						vol_file = "volume_core.hdf"
					elif (Tracker["constants"]["num_core_set"]==-1) and (len(cluster) > max(Tracker["minimum_grp_size"][0], 100)): 
						vol_file = "volume_core.hdf"
					else: vol_file = "               "
					msg = '{:>8} {:>8}   {:^24}        {:^6}          {:^6}          {:>5} {:^20} {:^20} '.format(nclusters, \
					  len(cluster), Tracker["current_generation"], 0.0, 0.0, 0.0, cluster_file,  vol_file)
					log_file.add(msg)
				except: log_file.add("Core_set.txt does not exist")
			
				NUACC = Tracker["constants"]["total_stack"] - NACC
				log_file.add('{:^7} {:^8} {:^22} {:^8} {:^24} {:^8} '.format(' Images', \
				   Tracker["constants"]["total_stack"], 'accounted for images: ', NACC, 'unaccounted for images: ', NUACC))
				log_file.add('The unaccounted for images are saved in Core_set.txt')
				if (len(clusters) >=2):do_analysis_on_identified_clusters(clusters, log_file)
				else: log_file.add(' ANOVA analysis is skipped ')
				with open(os.path.join(Tracker["constants"]["masterdir"], "Tracker.json"), 'w') as fout:
					json.dump(Tracker, fout)
				fout.close()
			else: log_file.add('No groups are found. Increase your img_per_grp and rerun program. \n')
		mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	return
####----------------------------------------------------------	
def output_clusters(output_dir, partition, unaccounted_list, \
    not_include_unaccounted, log_main):
	global Tracker, Blockdata
	### Single cpu function per generation  -------------------
	nclasses, npart = split_partition_into_ordered_clusters(partition)
	if (not Tracker["no_cluster_generation"]):
		nc = 0
		identified_clusters = []
		largest  = -float("inf")
		smallest =  float("inf") 
		for any in nclasses:
			largest  = max(largest, len(any))
			smallest = min(smallest, len(any))
		for ic in range(len(nclasses)):
			if (largest/float(len(nclasses[ic]))< 3.0):# Can be adaptive when criterion is established.
				sp_utilities.write_text_file(nclasses[ic], os.path.join(output_dir,"Cluster_%03d.txt"%nc))
				nc +=1
				identified_clusters.append(nclasses[ic])
			else: unaccounted_list +=nclasses[ic]
		log_main.add("Output  %d clusters. "%nc)
		###------------------------------------------------------
		if (len(unaccounted_list) > 1):
			unaccounted_list = sorted(unaccounted_list)
			sp_utilities.write_text_file(unaccounted_list, os.path.join(output_dir, "Core_set.txt"))
		nclasses = copy.deepcopy(identified_clusters)
		del identified_clusters ###-------------------------------
		if (len(unaccounted_list)>1): # output unaccounted as the last cluster
			if (not not_include_unaccounted):
				sp_utilities.write_text_file(unaccounted_list, os.path.join(output_dir,"Cluster_%03d.txt"%nc))
		if not not_include_unaccounted: 
			unclasses = copy.deepcopy(nclasses)
			unclasses.append(unaccounted_list)
			alist, partition = merge_classes_into_partition_list(unclasses)
			del unclasses
		else: alist, partition = merge_classes_into_partition_list(nclasses)
		sp_utilities.write_text_row(partition, os.path.join(output_dir, "final_partition.txt"))
	else:
		nc = 0
		for ic in range(len(nclasses)): unaccounted_list +=nclasses[ic]
		unaccounted_list = sorted(unaccounted_list)
		sp_utilities.write_text_file(unaccounted_list, os.path.join(output_dir, "Core_set.txt"))
		nclasses = [[]]
		log_main.add("Output  no clusters. ")
	return nclasses, nc, len(unaccounted_list)
###-----------------------------------------------------
def compute_final_map(work_dir, log_main):
	global Tracker, Blockdata
	###
	Tracker["constants"]["orgres"]			 = 0.0
	Tracker["constants"]["refinement_delta"] = 0.0
	Tracker["constants"]["refinement_ts"]	 = 0.0
	Tracker["constants"]["refinement_xr"]	 = 0.0
	Tracker["constants"]["refinement_an"]	 = 0.0
	number_of_groups                         = 0
	Tracker["Core_set"]                      = False
	###--------------------------------------------------
	if(Blockdata["myid"] == Blockdata["main_node"]):
		final_accounted_ptl = 0
		with open(os.path.join(work_dir, "Tracker.json"),'w') as fout:
			json.dump(Tracker, fout)
		fout.close()
		clusters = []
		while os.path.exists(os.path.join(work_dir, "Cluster_%03d.txt"%number_of_groups)):
			class_in = sp_utilities.read_text_file(os.path.join(work_dir, "Cluster_%03d.txt"%number_of_groups))
			number_of_groups += 1
			final_accounted_ptl +=len(class_in)
			clusters.append(class_in)
			del class_in
		try: # check the unaccounted for particles
			class_in = sp_utilities.read_text_file(os.path.join(work_dir, "Core_set.txt"))
			log_main.add('Core_set contains  %d images'%len(class_in))
			log_main.add('Current minimum_grp_size is [%d,%d]'%(Tracker["minimum_grp_size"][0],\
			    Tracker["minimum_grp_size"][1]))
			
			if (Tracker["constants"]["num_core_set"] ==-1):
				if len(class_in)> max(100, Tracker["minimum_grp_size"][0]):
					Tracker["Core_set"] = True
					clusters.append(class_in)
					number_of_groups += 1
					final_accounted_ptl  += len(class_in)
				else: Tracker["Core_set"] = False
			else: 
				if len(class_in) >=Tracker["constants"]["num_core_set"]:
					Tracker["Core_set"] = True
					clusters.append(class_in)
					number_of_groups    += 1
					final_accounted_ptl += len(class_in)
				else:Tracker["Core_set"]= False
			del class_in
		except:
			Tracker["Core_set"] = False
			log_main.add("Core_set.txt file does not exist.")
		Tracker["total_stack"]      = final_accounted_ptl
		Tracker["number_of_groups"] = number_of_groups
		Tracker["nxinit"]           = Tracker["nxinit_refinement"]
		log_main.add("Reconstruct %d maps in the final reconstruction."%number_of_groups)
		if Tracker["Core_set"]:log_main.add("Core_set is reconstructed.")
		else:                  log_main.add("Core_set is not reconstructed.")
	else: Tracker = 0
	number_of_groups = sp_utilities.bcast_number_to_all(number_of_groups, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	Tracker          = sp_utilities.wrap_mpi_bcast(Tracker,               Blockdata["main_node"], mpi.MPI_COMM_WORLD)
	if(number_of_groups == 0):
		sp_global_def.ERROR("No clusters  found, the program terminates.", "compute_final_map", 1, Blockdata["myid"])
	compute_noise(Tracker["nxinit"])
	if(Blockdata["myid"] == Blockdata["main_node"]):
		alist, partition = merge_classes_into_partition_list(clusters, do_sort=False)
		sp_utilities.write_text_row(partition, os.path.join(work_dir, "generation_partition.txt"))
	parti_file      = os.path.join(work_dir, "generation_partition.txt")
	params          = os.path.join(Tracker["constants"]["masterdir"],"refinement_parameters.txt")
	previous_params = Tracker["previous_parstack"]
	original_data, norm_per_particle = read_data_for_sorting(parti_file, params, previous_params)
	
	if Tracker["nosmearing"]: parameterstructure  = None
	else: parameterstructure  = read_paramstructure_for_sorting(parti_file, \
	    Tracker["paramstructure_dict"], Tracker["paramstructure_dir"])
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	
	Tracker["directory"] = work_dir
	compute_noise(Tracker["nxinit"])
	rdata = downsize_data_for_rec3D(original_data, Tracker["nxinit"], False, 1)# pay attentions to shifts!
	for im in range(len(original_data)):rdata[im].set_attr('group', original_data[im].get_attr('group'))
	del original_data
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	do3d_sorting_groups_nofsc_final(rdata, parameterstructure, norm_per_particle)
	del rdata
	if Blockdata["myid"] == Blockdata["main_node"]:
		if os.path.exists(os.path.join(Tracker["constants"]["masterdir"], 'tempdir')):
			shutil.rmtree(os.path.join(Tracker["constants"]["masterdir"], 'tempdir'))
	mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
	return
####--------------------------------------------------------------------------------------
def output_iter_results(box_dir, ncluster, NACC, NUACC, \
    minimum_grp_size, list_of_stable, unaccounted_list, log_main):
	### Single node
	new_list    = []
	nc          = 0
	NACC        = 0
	try:
		with open(os.path.join(box_dir,"freq_cutoff.json"),'r') as fout:
			freq_cutoff_dict = sp_utilities.convert_json_fromunicode(json.load(fout))
		fout.close()
	except: freq_cutoff_dict = {}
	
	for index_of_any in range(len(list_of_stable)):
		any = np.sort(list_of_stable[index_of_any])
		any.tolist()
		new_list.append(any)
		sp_utilities.write_text_file(any, os.path.join(box_dir, "Cluster_%03d.txt"%ncluster))
		freq_cutoff_dict["Cluster_%03d.txt"%ncluster] = Tracker["freq_fsc143_cutoff"]
		ncluster += 1
		nc       += 1
		NACC +=len(any)
	unaccounted_list = sorted(unaccounted_list)
	NUACC = len(unaccounted_list)
	with open(os.path.join(box_dir, "freq_cutoff.json"),'w') as fout:
		json.dump(freq_cutoff_dict, fout)
	fout.close()
	return ncluster, NACC, NUACC, unaccounted_list, nc	

#######================>>> FSC and group size related functions <<<=========================	
def get_group_size_from_iter_assign(iter_assignment):
	group_sizes  = []
	my_assign    = np.array(iter_assignment, dtype=np.int32)
	group_ids    = np.sort(np.unique(my_assign))
	for im in range(group_ids.shape[0]):
		group_sizes.append((np.sort(my_assign[isin(my_assign, group_ids[im])])).shape[0])
	return group_sizes

def set_minimum_group_size(log_main, printing_dres_table = True):
	global Tracker, Blockdata
	#------------>>> Empirical parameters <<<---------------
	# Guiding rule: high should be large, low should be small
	min_size_low_bound1  = 0.2
	min_size_low_bound2  = 0.3
	min_size_low_bound   = 0.25
	min_size_high_bound1 = 0.75
	min_size_high_bound2 = 0.90
	#-------------------------------------------------------
	total_stack      = Tracker["total_stack"]
	img_per_grp      = Tracker["current_img_per_grp"]
	number_of_groups = Tracker["number_of_groups"]
	try:    mu = Tracker["mu"]
	except: mu = Tracker["constants"]["mu"]
	current_fsc_curve = Tracker["constants"]["fsc_curve"]
	dres_table        = get_current_max_diff(current_fsc_curve, \
	    Tracker["constants"]["total_stack"],total_stack, number_of_groups)
	######----------------------------------------------------------------
	if mu<0:
		min_dis  = float("inf")
		selected = 0
		for il in range(len(dres_table)):
			diff = abs(dres_table[il][3]/float(dres_table[il][4]) - 1./(number_of_groups+1.))
			if (diff < min_dis):
				min_dis  = diff
				selected = il  # Always move forward one to get sufficiently small minimum grp size
		mu       = dres_table[min(selected, len(dres_table)-1)][0]
		min_size = dres_table[selected][3]
		min_res, mlow, mhigh = find_bounds(current_fsc_curve, total_stack, min_size)
		
	else:###
		mlow, mhigh, min_res = 0, 0, 0
		min_dis  = float("inf")
		selected = 0	
		for il in range(len(dres_table)):
			diff = abs(dres_table[il][0] - mu)
			if (diff < min_dis):
				min_dis  = diff
				selected = il # Always move forward one to get sufficiently small minimum grp size
		mu = dres_table[min(selected, len(dres_table)-1)][0]
		min_size = dres_table[selected][3]
		min_res, mlow, mhigh = find_bounds(current_fsc_curve, total_stack, min_size)

	Tracker["mu"]= max(1, mu)# mu at least be 1. 
	
	###--AI to determine mlow and mhigh	--------------
	img_per_grp = total_stack//number_of_groups
	l1          =  mlow/float(img_per_grp)*100.
	l2          = mhigh/float(img_per_grp)*100.
	###-----------------------------------------------
	log_main.add("Minimum grp size search gives low: %8d  high:  %8d  min_res: %3d  low:  %6.1f   high: %6.1f  "%\
	   (mlow, mhigh, min_res, l1, l2))
	   
	if(mlow >=mhigh): mhigh = mlow + int(0.1*img_per_grp) # Use img_per_grp as a ruler
	
	if (mhigh - mlow)/float(img_per_grp) > 0.25: # large gap between two bounds; tolerate more
		mhigh = int(min_size_high_bound2*img_per_grp)
		mlow  = min(img_per_grp//2, max(mlow, img_per_grp//3))
		mhigh = max(int(0.90*img_per_grp), mhigh)
				
		l1 = mlow/float(img_per_grp)*100.
		l2 = mhigh/float(img_per_grp)*100.
		log_main.add("Case1 low: %8d  high:  %8d  min_res: %3d   low:  %6.1f   high: %6.1f "%\
	       (mlow, mhigh, min_res, l1, l2 ))
		Tracker["large_gap"] = True
	   
	else:# high resolution, low gap between two bounds
		if (Tracker["freeze_groups"] == 0):
			mhigh = max(int(img_per_grp*.90), mhigh)# need a large gap between two bounds
			mlow  = min(img_per_grp//3, mlow)
		l1 = mlow/float(img_per_grp)*100.
		l2 = mhigh/float(img_per_grp)*100.
		log_main.add("Case2 low: %8d  high:  %8d  min_res: %3d  low:  %6.1f   high: %6.1f "%\
	      (mlow, mhigh, min_res, l1, l2))
	   	Tracker["large_gap"]    = False
	Tracker["minimum_grp_size"] = [mlow, mhigh, min_res, total_stack]
	####---------------------------------------------------------------------------------------------------
	if printing_dres_table:
		mlow, mhigh, min_res = Tracker["minimum_grp_size"][0], Tracker["minimum_grp_size"][1], Tracker["minimum_grp_size"][2]
		log_main.add("--------------------------------------------------------------------------------------")
		log_main.add("                           Dres_table                                                 ")
		log_main.add("FSC0.143 differences (mu) between the (K-1) small groups and one large group ")
		log_main.add("Total number of images: %8d   number of groups: %3d  "%(total_stack, number_of_groups))
		log_main.add("Mu   down   up  (to total/K)  small size   large size  ratio")
		log_main.add("--------------------------------------------------------------------------------------")
		for il in range(len(dres_table)):
			log_main.add("%3d   %3d     %3d        %10d      %10d    %5.4f"%(dres_table[il][0], dres_table[il][1], \
			  dres_table[il][2], dres_table[il][3], dres_table[il][4], dres_table[il][3]/float(dres_table[il][4])))
	if mu<Tracker["constants"]["mu"]: log_main.add("The user-provided mu is too large. It is changed to %d"%mu)
	log_main.add("Determined minimum_grp_size in range [%d,%d] at FSC0.143 %d  with given mu=%d, and total= %d   K=%d\n"%\
	(mlow, mhigh, min_res, mu, Tracker["total_stack"], Tracker["number_of_groups"]))
	return
####--------------------------------------------------------------------------------------
def get_res143(fsc_curve):
	fsc143 = len(fsc_curve)
	for ifreq in range(len(fsc_curve)):
		if (fsc_curve[ifreq] < 0.143):
			fsc143 = ifreq -1
			break
	return fsc143

def find_bounds(global_fsc, total_num, group_size):
	step = max(group_size//1000, 5)
	low_bound  = group_size
	high_bound = group_size
	cfsc       = sp_statistics.scale_fsc_datasetsize(global_fsc, total_num, group_size)
	res0       = get_res143(cfsc)
	res_low    = res0
	res_high   = res0
	while (res_low == res0):
		low_bound -=step
		res_low = get_res143(sp_statistics.scale_fsc_datasetsize(global_fsc, total_num, low_bound))
	while (res_high == res0):
		high_bound +=step
		res_high = get_res143(sp_statistics.scale_fsc_datasetsize(global_fsc, total_num, high_bound))
	return res0, int(low_bound), min(int(high_bound), total_num)

def get_current_max_diff(fsc_curve, org_total, NT, K):
	min_size   = NT//K
	max_size   = NT//K
	step       = max(min_size//1000, 5)
	dres_table = []
	pdres = -1 
	res0 = get_res143(sp_statistics.scale_fsc_datasetsize(fsc_curve, org_total, NT//K))
	while min_size >step:
		min_size -=step
		max_size +=step*(K-1)
		res1 = get_res143(sp_statistics.scale_fsc_datasetsize(fsc_curve, org_total, min_size))
		res2 = get_res143(sp_statistics.scale_fsc_datasetsize(fsc_curve, org_total, max_size))
		dres = res2-res1
		if dres == pdres:
			del dres_table[len(dres_table)-1]
		dres_table.append((dres, res0 -res1, res2 -res0, min_size, \
		    max_size, (max_size-min_size)))
		pdres = dres
	return dres_table
	
def compute_nstep(input_mdict, nstep=2):
	mlist = [0 for im in range(nstep)]
	mdict = {}
	
	if nstep == 1:
		mlist[0] = input_mdict[1]
		mdict[str(mlist[0])] = 0
	elif nstep ==2:
		mlist[0] = input_mdict[1]
		mdict[str(mlist[0])] = 0
		mlist[1] = input_mdict[0]
		mdict[str(mlist[1])] = 1
	elif nstep ==3:
		mlist[0] = input_mdict[1]
		mdict[str(mlist[0])] = 0
		mlist[2] = input_mdict[0]
		mdict[str(mlist[2])] = 2
		mlist[1] = (input_mdict[1]+input_mdict[0])//2
		mdict[str(mlist[1])] = 1
	else:
		mlist[0]  = input_mdict[1]
		mlist[-1] = input_mdict[0]
		kstep = 0
		t0 = mlist[0]
		t1 = mlist[-1]
		while kstep< nstep-1:
			kstep +=1
			mlist[kstep] = (t0+t1)//2
			t0 = mlist[kstep]
	return mlist, mdict

def calculate_grp_size_variation(log_file, state = "box_init"):
	global Tracker, Blockdata 
	# Single CPU
	####
	total_stack = Tracker["total_stack"]
	if state == "box_init":
		number_of_groups    = Tracker["number_of_groups"]
		current_img_per_grp = Tracker["current_img_per_grp"]
		
	elif state == "sort_init":
		current_img_per_grp = Tracker["current_img_per_grp"]
		number_of_groups    = Tracker["total_stack"]//current_img_per_grp
		current_img_per_grp = Tracker["total_stack"]//number_of_groups
		Tracker["number_of_groups"] = number_of_groups
		
	elif state == "run_init":
		number_of_groups    = Tracker["number_of_groups"]
		current_img_per_grp = Tracker["total_stack"]//number_of_groups
		Tracker["current_img_per_grp"] = current_img_per_grp
		
	elif state == "number_of_groups":
		number_of_groups    = Tracker["number_of_groups"]
		current_img_per_grp = Tracker["total_stack"]//number_of_groups
		Tracker["current_img_per_grp"] = current_img_per_grp
	#####
	min_res, mlow, mhigh = find_bounds(\
	   sp_statistics.scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], Tracker["constants"]["total_stack"], \
	   total_stack), total_stack, current_img_per_grp)
	   
	Tracker["img_per_grp"] = [mlow, mhigh, min_res, total_stack, number_of_groups]
	mlow1  = total_stack//(Tracker["number_of_groups"]+1) -1
	if number_of_groups>1: mhigh1 = min(total_stack//(Tracker["number_of_groups"]-1)+1, total_stack)
	else:                  mhigh1 = total_stack//Tracker["number_of_groups"]
	low_res1  = get_res143(sp_statistics.scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], Tracker["constants"]["total_stack"], mlow1))
	high_res1 = get_res143(sp_statistics.scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], Tracker["constants"]["total_stack"], mhigh1))
	log_file.add("Images group size range [%d, %d]  (Fall into the same Fourier pixel) at FSC0.143  %5d   "%(mlow1, mhigh1, high_res1))
	Tracker["img_per_grp"] = [max(mlow, mlow1),  mhigh1,  min_res, total_stack, Tracker["number_of_groups"]]
	log_file.add("                                                                   ")
	log_file.add("Current_img_per_grp = total/K=%d; total=%d; K=%d "%(Tracker["current_img_per_grp"], total_stack, number_of_groups))
	return
#  End of various utilities
###---------------------------------------------------------------------------------------
def AI_within_box(dtable, cutoff_ratio, cutoff_size, cutoff_mu):
	dtable = np.array(dtable, dtype=np.float64)
	(ntp, nsize) = dtable.shape
	# size, rep, cres, size ratio
	bratio  = dtable[3] < cutoff_ratio 
	bsize   = dtable[0] > cutoff_size
	cresmax = np.amax(dtable[2])
	#cresmin = np.minimum(dtable[2])
	bmu     = np.ones(nsize, dtype=bool)
	for im in range(nsize):
		 if cresmax - dtable[2][im]>cutoff_mu:
		 	bmu[im]= False
	bratio_mu   = np.multiply(bratio, bmu)
	bratio_size = np.multiply(bratio, bsize)
	return bratio_mu, bratio_size
	
def AI_split_groups(group_size):
	if len(group_size)>2:
		coid  = [max(group_size),  min(group_size)]
		glist = [None for im in range(len(group_size))]
		iter  = 0
		while iter<=3:
			for im in range(len(group_size)):
				if abs(float(coid[0])/group_size[im])<abs( group_size[im]/float(coid[1])):
					glist[im] = 0
				else: glist[im] = 1
			nc   = [0,   0]
			coid = [0,   0]	
			for ig in range(len(group_size)):
				if glist[ig]== 0:
					coid[0] +=group_size[ig]
					nc [0]  +=1
				else:
					coid[1] +=group_size[ig]
					nc [1]  +=1
			for ic in range(2):
				coid[ic] = 	coid[ic]//nc[ic]
			iter +=1
	else: 
		glist = [0, 1]
		coid  = [max(group_size),  min(group_size)]
	return glist, coid
###======================================================================================= 
def main():
	global Tracker, Blockdata
	
	progname = os.path.basename(sys.argv[0])
	usage = progname + " --refinement_dir=masterdir_of_sxmeridien   --output_dir=sort3d_output --mask3D=mask.hdf --focus=binarymask.hdf  --radius=outer_radius " +\
	"  --sym=c1  --img_per_grp=img_per_grp  --minimum_grp_size=minimum_grp_size "
	parser = optparse.OptionParser(usage,version=sp_global_def.SPARXVERSION)
	parser.add_option("--refinement_dir", type ="string",  default ='',  help="sxmeridien 3-D refinement directory")
	parser.add_option("--instack",        type ="string",  default ='',  help="file name, data stack for sorting provided by user. It applies when sorting starts from a given data stack")

	initiate_from_meridien_mode = False
	for q in sys.argv[1:]:
		if( q[:16] == "--refinement_dir" ):
			initiate_from_meridien_mode = True
			break
	
	initiate_from_data_stack_mode = False
	for q in sys.argv[1:]:
		if( q[:9] == "--instack" ):
			initiate_from_data_stack_mode = True
			break
	
  	# priority
	if initiate_from_data_stack_mode and initiate_from_meridien_mode:
		initiate_from_data_stack_mode = False
		
	if (not initiate_from_data_stack_mode) and (not initiate_from_meridien_mode):
		if Blockdata["myid"] == Blockdata["main_node"]:
			sxprint("Specify one of the two options to start the program: --refinement_dir, --instack")
	
	if initiate_from_meridien_mode:
		parser.add_option("--output_dir",                  type   ="string",        default ='',	    help="Sort3d output directory name")
		parser.add_option("--niter_for_sorting",           type   ="int",           default =-1,	    help="User specified iteration number of 3D refinement for sorting")
		parser.add_option("--focus",                       type   ="string",        default ='',        help="Focus 3D mask. File path of a binary 3D mask for focused clustering ")
		parser.add_option("--mask3D",                      type   ="string",        default ='',        help="3D mask. File path of the global 3D mask for clustering")
		parser.add_option("--radius",                      type   ="int",           default =-1,	    help="Estimated protein radius in pixels")
		parser.add_option("--sym",                         type   ="string",        default ='c1',      help="Point-group symmetry")
		parser.add_option("--img_per_grp",                 type   ="int",           default =-1,        help="Number of images per group, default value will activate automated group search")
		parser.add_option("--nsmear",                      type   ="float",         default =-1.,       help="Number of smears used in sorting. Fill it with 1 if user does not want to use all smears")
		parser.add_option("--memory_per_node",             type   ="float",         default =-1.0,      help="Memory_per_node, the number used for computing the CPUs/NODE settings given by user")
		parser.add_option("--orientation_groups",          type   ="int",           default = 20,       help="Number of orientation groups in the asymmetric unit")
		parser.add_option("--not_include_unaccounted",     action ="store_true",    default =False,     help="Do not reconstruct unaccounted elements in each generation")
		parser.add_option("--notapplybckgnoise",           action ="store_true",    default =False,     help="Do not applynoise")
		parser.add_option("--num_core_set",                type   ="int",           default =-1,	    help="Number of images for reconstructing core set images. Will not reconstruct core set images if the total number of core set images is less than this")
		parser.add_option("--check_smearing",              action ="store_true",    default =False,     help="Check the smearing per iteration and estimate the precalculated data")
		parser.add_option("--compute_on_the_fly",          action ="store_true",    default =False,     help="Pre-compute part of multiple assignment images for reconstruction till memory is filled up and then the other part on the fly")
		parser.add_option("--nstep",                       type   ="int",           default =5,		    help="Number of steps to decrease minimum group size from high bound to low bound")
		parser.add_option("--overhead",                    type   ="float",         default =5.0,       help="Estimated python overhead per node in GB")
		parser.add_option("--not_freeze_groups",           action ="store_true",    default =False,     help="Do not remove small groups during within-box clustering")
		parser.add_option("--use_umat",                    action ="store_true",    default =False,     help="Use fuzzy membership to stabalize sorting")

		(options, args) = parser.parse_args(sys.argv[1:])
		### Sanity check
	
		checking_flag = 0
		if Blockdata["myid"] == Blockdata["main_node"]:
			if options.refinement_dir !='':
				if not os.path.exists(options.refinement_dir): checking_flag = 1
				
		checking_flag = sp_utilities.bcast_number_to_all(checking_flag, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
		if checking_flag ==1:
			sp_global_def.ERROR("The specified refinement_dir does not exist", "sort3d", 1, Blockdata["myid"])
	
		if options.focus !='':
			if Blockdata["myid"] == Blockdata["main_node"]:
				if not os.path.exists(options.focus):  checking_flag = 1
			checking_flag = sp_utilities.bcast_number_to_all(checking_flag, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			if checking_flag ==1:
				sp_global_def.ERROR("The specified focus mask file does not exist", "sort3d", 1, Blockdata["myid"])
		
		if options.mask3D !='':
			if Blockdata["myid"] == Blockdata["main_node"]:
				if not os.path.exists(options.mask3D): checking_flag = 1
			checking_flag = sp_utilities.bcast_number_to_all(checking_flag, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			if checking_flag ==1:
				sp_global_def.ERROR("The specified mask3D file does not exist", "sort3d", 1, Blockdata["myid"])
				
	
		#--- Fill input parameters into dictionary Constants
		Constants		                         = {}
		Constants["niter_for_sorting"]           = options.niter_for_sorting
		Constants["memory_per_node"]             = options.memory_per_node
		Constants["orgstack"]                    = options.instack
		Constants["masterdir"]                   = options.output_dir
		Constants["refinement_dir"]              = options.refinement_dir
	
		if options.mask3D == '': Constants["mask3D"] = False
		else:                    Constants["mask3D"] = options.mask3D
		if options.focus!='':   Constants["focus3D"] = options.focus
		else:                   Constants["focus3D"] = False
	
		Constants["img_per_grp"]                 = options.img_per_grp
		Constants["mu"]      		             = -1
		Constants["radius"]              		 = options.radius
		Constants["sym"]                         = options.sym
		Constants["nsmear"]                      = options.nsmear

		Constants["restart_from_nbox"]           =  0 #options.restart_from_nbox
		Constants["restart_from_depth_order"]    = -1 #options.restart_from_depth_order
		Constants["restart_from_generation"]     = -1 #options.restart_from_generation
	
		#### options for advanced users
		Constants["not_include_unaccounted"]     = False
		Constants["nxinit"]                      = -1
		### Frozen options
		Constants["upscale"]                     = 0.5 #
		Constants["interpolation"]               = "trl"
		Constants["comparison_method"]           = "cross" #options.comparison_method # either cross or eucd
		Constants["symmetry"]                    = Constants["sym"]
		Constants["CTF"]                		 = True
		Constants["do_not_use_3dmask"]           = False 
	
		if options.focus:  Constants["comparison_method"] = "cross" # in case of focus3D, cross is used.
		Constants["fuse_freq"]           = 45.  # Now in A, convert to pixels before being used
		Constants["orientation_groups"]  = options.orientation_groups # orientation constrained angle step
		Constants["num_core_set"]        = options.num_core_set
		Constants["check_smearing"]      = options.check_smearing
		Constants["compute_on_the_fly"]  = options.compute_on_the_fly
		Constants["overhead"]            = max(options.overhead, 1.)
		#
		#
		# Create and initialize Tracker dictionary with input options  # State Variables	
		Tracker                     = {}
		Tracker["constants"]	    = Constants
		Tracker["radius"]           = Tracker["constants"]["radius"]
		Tracker["upscale"]          = Tracker["constants"]["upscale"]
		Tracker["applyctf"]         = False  # Should the data be premultiplied by the CTF.  Set to False for local continuous.
		Tracker["nxinit"]           = Tracker["constants"]["nxinit"]
		if options.notapplybckgnoise: Tracker["applybckgnoise"] = False
		else:                         Tracker["applybckgnoise"] = True
		### ------------------------------------------------------------------------------
		if Tracker["constants"]["memory_per_node"] == -1 or Tracker["constants"]["memory_per_node"] <32.:
			Tracker["constants"]["small_memory"]   = True
		else: Tracker["constants"]["small_memory"] = False
		
		Tracker["constants"]["hardmask"]           = True
		Tracker["applymask"]                       = True
		Tracker["constants"]["refinement_method"]  = "SPARX" 
		Tracker["nosmearing"]                      = False
		checking_flag                              = 0 # reset
		Blockdata["fftwmpi"]                       = True
		
		###----------------- Set sorting mode and nstep ----------------------------------
		# Three modes: 1. freeze small groups; 2. not freeze; 3. search mode
		Tracker["box_learning"]  = 1
		Tracker["nstep"]         = max(options.nstep, 1)
		Tracker["search_mode"]   = True		
		if options.not_freeze_groups: Tracker["freeze_groups"] = 0
		else:                         Tracker["freeze_groups"] = 1
		if Tracker["constants"]["img_per_grp"]>1: Tracker["search_mode"] = False
		if Tracker["search_mode"]:
			Tracker["constants"]["init_K"] = 5
			Tracker["freeze_groups"]       = 0
		else:Tracker["constants"]["init_K"] = -1
		Tracker["use_umat"]      =  options.use_umat
		###-------------------------------------------------------------------------------
				
		try : 
			Blockdata["symclass"] = sp_fundamentals.symclass(Tracker["constants"]["symmetry"])
			Tracker["orientation_groups"] = \
			    max(4, Tracker["constants"]["orientation_groups"]//Blockdata["symclass"].nsym)
		except: pass
		Blockdata["ncpuspernode"] = Blockdata["no_of_processes_per_group"]
		Blockdata["nsubset"]      = Blockdata["ncpuspernode"]*Blockdata["no_of_groups"]
		
		create_subgroup()
		create_zero_group()		
		#
		#
		continue_from_interuption = 0
		continue_from_interuption = create_masterdir()
		log_main                  = sp_logger.Logger(sp_logger.BaseLogger_Files())
		log_main.prefix           = Tracker["constants"]["masterdir"]+"/"
		if continue_from_interuption == 0:
			if Blockdata["myid"] == Blockdata["main_node"]:
				log_main.add('================================================================================================================')
				log_main.add('                                 SORT3D MULTI-LAYER v2.0')
				log_main.add('================================================================================================================')
			import_data(log_main)
			print_shell_command(sys.argv, log_main)
			set_sorting_global_variables_mpi(log_main)
			check_3dmask(log_main)
			check_mpi_settings(log_main)
			keepsorting = sort3d_init("Initialization", log_main)
			dump_tracker(Tracker["constants"]["masterdir"])
			if keepsorting == 0:
				mpi.mpi_finalize()
				exit()				 
		sorting_main_mpi(log_main, options.not_include_unaccounted)
		if Blockdata["myid"] == Blockdata["main_node"]:
			log_main.add('----------------------------------------------------------------------------------------------------------------' )
			log_main.add('                                 SORT3D MULTI-LAYER finished')
			log_main.add('----------------------------------------------------------------------------------------------------------------' )
		mpi.mpi_finalize()
		exit()
			
	elif initiate_from_data_stack_mode:
		parser.add_option("--nxinit",                            type   ="int",           default =-1,         help="User provided image size")
		parser.add_option("--output_dir",                        type   ="string",        default ='',		   help="Name of the sort3d directory")
		parser.add_option("--focus",                             type   ="string",        default ='',         help="Focus 3D mask. File path of a binary 3D mask for focused clustering ")
		parser.add_option("--mask3D",                            type   ="string",        default ='',         help="3D mask. File path of the global 3D mask for clustering")
		parser.add_option("--radius",                            type   ="int",           default =-1,	       help="Estimated protein radius in pixels")
		parser.add_option("--sym",                               type   ="string",        default ='c1',       help="Point-group symmetry")
		parser.add_option("--img_per_grp",                       type   ="int",           default =-1,         help="Number of images per group, default value will activate automated group search")
		parser.add_option("--nsmear",                            type   ="float",         default =10.,        help="Number of smears used in sorting. Fill it with 1 if user does not want to use all smears")
		parser.add_option("--memory_per_node",                   type   ="float",         default =-1.0,       help="Memory_per_node, the number used for computing the CPUs/NODE settings given by user")
		parser.add_option("--orientation_groups",                type   ="int",           default = 20,        help="Number of orientation groups in the asymmetric unit")
		parser.add_option("--not_include_unaccounted",           action ="store_true",    default =False,      help="Do not reconstruct unaccounted elements in each generation")
		parser.add_option("--notapplybckgnoise",                 action ="store_true",    default =False,      help="Flag to turn off background noise")
		parser.add_option("--num_core_set",                      type   ="int",           default =-1,		   help="Number of images for reconstructing core set images. Will not reconstruct core set images if the total number of core set images is less than this")
		parser.add_option("--nstep",                             type   ="int",           default =7,		   help="number of steps to decrease minimum group size from high bound to low bound")
		parser.add_option("--overhead",                          type   ="float",         default =5.0,        help="Estimated python overhead per node in GB")
		parser.add_option("--not_freeze_groups",                 action ="store_true",    default =False,      help="Do not remove small groups during within-box clustering")
		parser.add_option("--use_umat",                          action ="store_true",    default =False,      help="Use fuzzy membership to stabilize sorting")
		

		(options, args) = parser.parse_args(sys.argv[1:])
		### Sanity check
		
		checking_flag = 0
		if options.focus !='':
			if Blockdata["myid"] == Blockdata["main_node"]:
				if not os.path.exists(options.focus):  checking_flag = 1
			checking_flag = sp_utilities.bcast_number_to_all(checking_flag, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			if checking_flag ==1:
				sp_global_def.ERROR("The specified focus mask file does not exist", "sort3d", 1, Blockdata["myid"])
		
		if options.mask3D !='':
			if Blockdata["myid"] == Blockdata["main_node"]:
				if not os.path.exists(options.mask3D): checking_flag = 1
			checking_flag = sp_utilities.bcast_number_to_all(checking_flag, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
			if checking_flag ==1:
				sp_global_def.ERROR("The specified mask3D file does not exist", "sort3d", 1, Blockdata["myid"])
					
				
		#--- Fill input parameters into dictionary Constants------------------------------
		Constants		                             = {}
		Constants["memory_per_node"]                 = options.memory_per_node
		Constants["orgstack"]                        = options.instack
		Constants["masterdir"]                       = options.output_dir
		
		if options.mask3D == '': Constants["mask3D"] = False
		else:                    Constants["mask3D"] = options.mask3D
		if options.focus!='':   Constants["focus3D"] = options.focus
		else:                   Constants["focus3D"] = False
		
		Constants["nsmear"]                      = 1
		Constants["img_per_grp"]                 = options.img_per_grp
		Constants["mu"]                          = -1
		Constants["radius"]              		 = options.radius
		Constants["sym"]                         = options.sym
	
		Constants["restart_from_nbox"]           =  0 #options.restart_from_nbox
		Constants["restart_from_depth_order"]    = -1 #options.restart_from_depth_order
		Constants["restart_from_generation"]     = -1 #options.restart_from_generation
	
		#### options for advanced users
		Constants["not_include_unaccounted"]     = False
		Constants["nxinit"]                      = options.nxinit
		
		### Frozen options
		Constants["upscale"]                     = 0.5 #
		Constants["interpolation"]               = "trl"
		Constants["comparison_method"]           = "cross" #options.comparison_method # either cross or eucd
		Constants["symmetry"]                    = Constants["sym"]
		Constants["CTF"]                		 = True
		Constants["do_not_use_3dmask"]           = False
	
		if options.focus:  Constants["comparison_method"] = "cross" # in case of focus3D, cross is used.
		Constants["fuse_freq"] = 45.  # Now in A, convert to pixels before being used
		Constants["orientation_groups"]  = options.orientation_groups # orientation constrained angle step
		Constants["num_core_set"]        = options.num_core_set
		Constants["compute_on_the_fly"]  = False
		Constants["overhead"]            = max(options.overhead, 1.)
		#
		#
		# Create and initialize Tracker dictionary with input options  # State Variables<<<------------	
		Tracker                     = {}
		Tracker["constants"]	    = Constants
		Tracker["radius"]           = Tracker["constants"]["radius"]
		Tracker["upscale"]          = Tracker["constants"]["upscale"]
		Tracker["applyctf"]         = False  # Should the data be premultiplied by the CTF.  Set to False for local continuous.
		Tracker["nxinit"]           = Tracker["constants"]["nxinit"]
		if options.notapplybckgnoise: Tracker["applybckgnoise"] = False
		else:                         Tracker["applybckgnoise"] = True
	
		### ------------=====< option for proteins images that have preferred orientations<<<------------
		 # for orientation groups
		if (Tracker["constants"]["memory_per_node"] == -1) or (Tracker["constants"]["memory_per_node"] <32.):
			Tracker["constants"]["small_memory"]   = True
		else: Tracker["constants"]["small_memory"] = False
		## additional check
		
		###----------------- Set sorting mode and nstep ----------------------------------
		# Three modes: 1. freeze small groups; 2. not freeze; 3. search mode
		Tracker["box_learning"]  = 1
		Tracker["nstep"]         = max(options.nstep, 1)
		Tracker["search_mode"]   = True		
		if options.not_freeze_groups: Tracker["freeze_groups"] = 0
		else:                         Tracker["freeze_groups"] = 1
		
		if (Tracker["constants"]["img_per_grp"]>1):Tracker["search_mode"] = False
		if Tracker["search_mode"]: 
			Tracker["constants"]["init_K"]   = 7
			Tracker["freeze_groups"]         = 0
		else: Tracker["constants"]["init_K"] = -1
		Tracker["use_umat"]                  =  options.use_umat
		###-------------------------------------------------------------------------------
		
		Tracker["constants"]["hardmask"]           = True
		Tracker["applymask"]                       = True
		Tracker["constants"]["refinement_method"]  ="stack"
		Tracker["constants"]["refinement_dir"]     = None
		Tracker["paramstructure_dir"]              = None
		Tracker["refang"]                          = None
		Tracker["rshifts"]                         = None
		Tracker["paramstructure_dict"]             = None
		Tracker["constants"]["selected_iter"]      = -1
		Tracker["nosmearing"]                      = True
		checking_flag                              =  0 # reset
		Blockdata["fftwmpi"]                       = True
		Blockdata["symclass"]                      = sp_fundamentals.symclass(Tracker["constants"]["symmetry"])
		Tracker["orientation_groups"] = max(4, Tracker["constants"]["orientation_groups"]//Blockdata["symclass"].nsym)		
		Blockdata["ncpuspernode"]     = Blockdata["no_of_processes_per_group"]
		Blockdata["nsubset"]          = Blockdata["ncpuspernode"]*Blockdata["no_of_groups"]
		create_subgroup()
		create_zero_group()
		# 
		
		continue_from_interuption = 0
		continue_from_interuption = create_masterdir()
		log_main = sp_logger.Logger(sp_logger.BaseLogger_Files())
		log_main.prefix = Tracker["constants"]["masterdir"]+"/"
		if continue_from_interuption == 0:# Fresh run
			if Blockdata["myid"] == Blockdata["main_node"]:
				log_main.add('================================================================================================================')
				log_main.add('                                  SORT3D MULTI-LAYER v2.0')
				log_main.add('================================================================================================================\n')	
			import_data(log_main)
			print_shell_command(sys.argv, log_main)
			set_sorting_global_variables_mpi(log_main)
			check_3dmask(log_main)
			check_mpi_settings(log_main)
			keepsorting = sort3d_init("Initialization", log_main)
			dump_tracker(Tracker["constants"]["masterdir"])
			if (keepsorting ==0):
				mpi.mpi_finalize()
				exit()
		sorting_main_mpi(log_main, options.not_include_unaccounted)
		if Blockdata["myid"] == Blockdata["main_node"]:
			log_main.add('----------------------------------------------------------------------------------------------------------------' )
			log_main.add('                                 SORT3D MULTI-LAYER finished')
			log_main.add('----------------------------------------------------------------------------------------------------------------' )
if __name__ == "__main__":
	sp_global_def.print_timestamp( "Start" )
	main()
	sp_global_def.print_timestamp( "Finish" )
	mpi.mpi_finalize()
