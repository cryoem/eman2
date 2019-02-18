#!/usr/bin/env python
#
#
#  08/26/2016
#  New version of sort3D.
#  
from __future__ import print_function

from builtins import range

"""
There are two ways to run the program:

1.a Import data from a meridien refinement. Data are precomputed and save in the memory.
mpirun -np 64 --hostfile four_nodes.txt  sxsort3d_depth.py  --refinemet_dir=meridien_run --niter_for_sorting=28 \
      --memory_per_node=100. --img_per_grp=80000  --minimum_grp_size=15000   \
      --stop_mgskmeans_percentage=8. --output_dir=SORT3D

1.b Import data from a meridien refinement and handles large number of smearing.

mpirun -np 64 --hostfile four_nodes.txt  sxsort3d_depth.py  --refinemet_dir=meridien_run 
   --niter_for_sorting=28  --memory_per_node=100. --img_per_grp=80000  --minimum_grp_size=15000  \
   --stop_mgskmeans_percentage=8. --output_dir=SORT3D   --check_smearing  --compute_on_the_fly

2. Import data from a given data stack.
mpirun  -np 48  --hostfile ./node012.txt  sxsort3d_depth.py --stop_mgskmeans_percentage=15. \
     --orientation_cones=40  --do_swap_au  --swap_ratio=5. --output_dir=sorting_bmask04 --sym=c1  \
     --radius=30  --minimum_grp_size=2000   --img_per_grp=2800    --instack=bdb:data  >sorting_bmask04/printout &

3. Run the program on a single node workstation. The minimum requirement for a successful run is
    the 2D data size is less than the memory of the workstation. 

mpirun  -np 8  --hostfile ./node0.txt  sxsort3d_depth.py --stop_mgskmeans_percentage=15. \
    --orientation_cones=40  --do_swap_au  --swap_ratio=5. --output_dir=sorting_bmask04 --sym=c1  \
    --radius=30  --minimum_grp_size=2000   --img_per_grp=2800    --instack=bdb:data  >sorting_bmask04/printout &

Notices on options:

a. --do_swap_au  --swap_ratio=5.   : ratio of the elements of determined clusters that are exchanged with the unaccounted elements.
b. --stop_mgskmeans_percentage=15. : criterion to stop Kmeans run.
c. --check_smearing : program prints out the averaged number of smearings in each iteration of meridien refinement
d. --compute_on_the_fly: it enables sorting done with iterations with large number of smearing. The shifted data are computed on the fly.


Nomenclatures in sorting intermediate results:

NACC:  total number of accounted
NUACC: number of unaccounted
MGU:   user defined minimum_grp_size  
MGR:   random group size
K:     number of groups
"""

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
from    mpi        import  *
from    math  	   import  *
from    random     import  *
import  shutil
import  os
import  sys
import  subprocess
import  string
import  json
from    sys 	import exit
from    time    import time, localtime, strftime, sleep
global  Tracker, Blockdata

mpi_init(0, [])
nproc     = mpi_comm_size(MPI_COMM_WORLD)
myid      = mpi_comm_rank(MPI_COMM_WORLD)
Blockdata = {}
#  MPI stuff
Blockdata["nproc"]              = nproc
Blockdata["myid"]               = myid
Blockdata["main_node"]          = 0
Blockdata["last_node"]          = nproc -1

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
	if Blockdata["myid"]== Blockdata["main_node"]: print("Not a unix machine")
	is_unix_cluster = False
	
def create_subgroup():
	# select a subset of myids to be in subdivision
	if( Blockdata["myid_on_node"] < Blockdata["ncpuspernode"] ): submyids = [Blockdata["myid"]]
	else:  submyids = []
	submyids = wrap_mpi_gatherv(submyids, Blockdata["main_node"], MPI_COMM_WORLD)
	submyids = wrap_mpi_bcast(submyids,   Blockdata["main_node"], MPI_COMM_WORLD)
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
	
def create_zero_group():
	# select a subset of myids to be in subdivision, This is a group of all zero IDs on nodes, taken from isac2
	if( Blockdata["myid_on_node"] == 0 ): submyids = [Blockdata["myid"]]
	else:  submyids = []

	submyids = wrap_mpi_gatherv(submyids,  Blockdata["main_node"], MPI_COMM_WORLD)
	submyids = wrap_mpi_bcast(submyids,    Blockdata["main_node"], MPI_COMM_WORLD)
	world_group = mpi_comm_group(MPI_COMM_WORLD)
	subgroup    = mpi_group_incl(world_group,len(submyids),submyids)
	Blockdata["group_zero_comm"] = mpi_comm_create(MPI_COMM_WORLD, subgroup)
	mpi_barrier(MPI_COMM_WORLD)

	Blockdata["group_zero_size"] = -1
	Blockdata["group_zero_myid"] = -1
	if (MPI_COMM_NULL != Blockdata["group_zero_comm"]):
		Blockdata["group_zero_size"] = mpi_comm_size(Blockdata["group_zero_comm"])
		Blockdata["group_zero_myid"] = mpi_comm_rank(Blockdata["group_zero_comm"])
	#  "nodes" are zero nodes on subgroups on the two "node_volume" that compute backprojection
	#Blockdata["nodes"] = [Blockdata["node_volume"][0]*Blockdata["ncpuspernode"], Blockdata["node_volume"][1]*Blockdata["ncpuspernode"]]
	mpi_barrier(MPI_COMM_WORLD)
	return
### utilities 
	
def dump_tracker(path_of_the_tracker):
	global Tracker, Blockdata
	if(Blockdata["myid"] == Blockdata["main_node"]):
		fout = open(os.path.join(path_of_the_tracker, "Tracker.json"),'w')
		fout.seek(0)
		json.dump(Tracker, fout)
		fout.truncate()
		fout.close()
	mpi_barrier(MPI_COMM_WORLD)
	return
				
######### depth clustering functions
def depth_clustering(work_dir, depth_order, initial_id_file, params, previous_params, log_main):
	global Tracker, Blockdata
	if(Blockdata["myid"] == Blockdata["main_node"]):
		if not os.path.exists(os.path.join(work_dir, "layer0")):
			os.mkdir(os.path.join(work_dir, "layer0"))
		partition_per_box_per_layer_list = []
		for iparti in range(0, 2**(depth_order+1), 2):
			partition_per_box_per_layer_list.append([read_text_file(initial_id_file), None])
	else: partition_per_box_per_layer_list = 0
	partition_per_box_per_layer_list = wrap_mpi_bcast(partition_per_box_per_layer_list, Blockdata["main_node"])
	keepchecking     = 1
	bad_clustering   = 0
	stat_list        = []
	Tracker["depth"] = 0
	for depth in range(depth_order): #  layers, depth_order = 1 means one layer and two boxes.
		time_layer_start = time()
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
				mark_sorting_state(depth_dir, False, log_main)
			else: 
				keepchecking = check_sorting_state(depth_dir, keepchecking, log_main)
				if keepchecking == 0: mark_sorting_state(depth_dir, False,  log_main)
		else: keepchecking = 0
		keepchecking = bcast_number_to_all(keepchecking, Blockdata["main_node"], MPI_COMM_WORLD)
		## box loop
		if keepchecking == 0:
			checkingbox = 1
			for nbox in range(n_cluster_boxes):
				time_box_start = time()
				input_accounted_file   = partition_per_box_per_layer_list[nbox][0]
				input_unaccounted_file = partition_per_box_per_layer_list[nbox][1]
				nbox_dir = os.path.join(depth_dir, "nbox%d"%nbox)
				if(Blockdata["myid"] == Blockdata["main_node"]):
					if not os.path.exists(nbox_dir):
						os.mkdir(nbox_dir)
						checkingbox = 0
						mark_sorting_state(nbox_dir, False, log_main)
					else:
						checkingbox = check_sorting_state(nbox_dir, checkingbox, log_main)
						if checkingbox == 0:# found not finished box
							msg ="Pair %d is not finished. Remove it and recompute..."%nbox
							log_main.add(msg)
							shutil.rmtree(nbox_dir)
							os.mkdir(nbox_dir)
							mark_sorting_state(nbox_dir, False, log_main)
						else:
							msg ="Pair %d is completed"%nbox
							log_main.add(msg)
				else: checkingbox = 0
				checkingbox = bcast_number_to_all(checkingbox, Blockdata["main_node"], MPI_COMM_WORLD)
				# Code structure of the box
				if checkingbox == 0:
					bad_box = depth_clustering_box(nbox_dir, input_accounted_file, input_unaccounted_file, \
					    params, previous_params, nbox, log_main)
					if(Blockdata["myid"] == Blockdata["main_node"]):
						mark_sorting_state(nbox_dir, True, log_main)
				else: continue
				if bad_box == 1: bad_clustering = 1
				time_box_finish = (time() - time_box_start)//60.
				if(Blockdata["myid"] == Blockdata["main_node"]):
					log_main.add(' Box %d of layer %d is done, calculations took %10.1f minutes'%(nbox, depth, time_box_finish))
					log_main.add('----------------------------------------------------------------------------------------------------------------\n\n')
		if( bad_clustering !=1 ):
			partition_per_box_per_layer_list = []
			stop_generation = 0
			for nbox in range(0,n_cluster_boxes, 2):
				input_box_parti1 = os.path.join(depth_dir, "nbox%d"%nbox,     "partition.txt")
				input_box_parti2 = os.path.join(depth_dir, "nbox%d"%(nbox+1), "partition.txt")
				minimum_grp_size, maximum_grp_size, accounted_list, unaccounted_list, bad_clustering, stop_generation, stat_list = \
				do_boxes_two_way_comparison_mpi(nbox, input_box_parti1, input_box_parti2, depth_order - depth, log_main)
				
				if(stop_generation == 1):
					partition_per_box_per_layer_list = []
					partition_per_box_per_layer_list.append([accounted_list, unaccounted_list])
					break
				else:
					partition_per_box_per_layer_list.append([accounted_list, unaccounted_list])
			mpi_barrier(MPI_COMM_WORLD)
			if(Blockdata["myid"] != Blockdata["main_node"]): Tracker = 0
			Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], MPI_COMM_WORLD)
			
			if(Blockdata["myid"] == Blockdata["main_node"]): mark_sorting_state(depth_dir, True, log_main)
			if( bad_clustering == 1):   break
			if( stop_generation == 1 ): break ### only one cluster survives
		else:
			if(Blockdata["myid"] == Blockdata["main_node"]):
				log_main.add('Layer %d sorting has completed. Recompute two_way comparison'%depth)
			partition_per_box_per_layer_list = []
			for nbox in range(0, n_cluster_boxes,2):
				input_box_parti1 = os.path.join(depth_dir, "nbox%d"%nbox,     "partition.txt")
				input_box_parti2 = os.path.join(depth_dir, "nbox%d"%(nbox+1), "partition.txt")
				minimum_grp_size, maximum_grp_size, accounted_list, unaccounted_list, bad_clustering, stop_generation, stat_list = \
				do_boxes_two_way_comparison_mpi(nbox, input_box_parti1, input_box_parti2, depth_order - depth, log_main)
				
				if(stop_generation == 1):
					partition_per_box_per_layer_list = []
					partition_per_box_per_layer_list.append([accounted_list, unaccounted_list])
					break
				else:
					partition_per_box_per_layer_list.append([accounted_list, unaccounted_list])
			mpi_barrier(MPI_COMM_WORLD)
			
		time_of_sorting_h,  time_of_sorting_m = get_time(time_layer_start)
		if Blockdata["myid"] == Blockdata["main_node"]:
			log_main.add('            Execution of layer %d took %5d hours %3d minutes.'%(depth, time_of_sorting_h, time_of_sorting_m))
			log_main.add(' ')
			log_main.add(' ')
			log_main.add('================================================================================================================')
	return partition_per_box_per_layer_list, bad_clustering, stat_list

def depth_box_initialization(box_dir, input_list1, input_list2, log_file):
	# input_list1:  accounted partition
	# input_list2:  unacccounted
	global Tracker, Blockdata
	import numpy as np
	img_per_grp      = Tracker["current_img_per_grp"]
	minimum_grp_size = Tracker["constants"]["minimum_grp_size"]
	if Tracker["constants"]["do_swap_au"]: swap_ratio = Tracker["constants"]["swap_ratio"]
	else: swap_ratio = 0.0
	if input_list2 is not None: #Track 2
		total_stack = len(input_list1)+ len(input_list2)
		groups = np.unique((np.array(input_list1, dtype=np.int32)).transpose()[0])
		number_of_groups = groups.shape[0]
		if Blockdata["myid"] == Blockdata["main_node"]:
			msg = "Number of groups found in initialization: %d, the largest possible number of groups: %d"%\
			  (number_of_groups, total_stack//img_per_grp)
			log_file.add(msg)
			write_text_row(input_list1,  os.path.join(box_dir, "previous_NACC.txt"))
			write_text_file(input_list2, os.path.join(box_dir, "previous_NUACC.txt"))
		mpi_barrier(MPI_COMM_WORLD)
		new_assignment   = []
		for indep in range(2):
			tmp_assignment = swap_accounted_with_unaccounted_elements_mpi(os.path.join(box_dir, "previous_NACC.txt"), \
				os.path.join(box_dir, "previous_NUACC.txt"), log_file, number_of_groups, swap_ratio)
			new_assignment.append(tmp_assignment)
			
		if Blockdata["myid"] == Blockdata["main_node"]:
			for indep in range(2):
				write_text_file(new_assignment[indep], os.path.join(box_dir, "independent_index_%03d.txt"%indep))
			new_assignment = []
			for indep in range(2):
				new_assignment.append(read_text_row(os.path.join(box_dir,"independent_index_%03d.txt"%indep)))
			id_list = read_text_file( os.path.join(box_dir, "independent_index_000.txt"), -1)
			if len(id_list)>1: id_list=id_list[0]
			total_stack = len(id_list)
			group_id = np.sort(np.unique(np.array(id_list, dtype=np.int32)))
			number_of_groups = group_id.shape[0] # assume K be 0, ...,number_of_groups-1
		else:
			number_of_groups = 0
			total_stack      = 0
			new_assignment   = 0
		new_assignment   = wrap_mpi_bcast(new_assignment,        Blockdata["main_node"], MPI_COMM_WORLD)
		number_of_groups = bcast_number_to_all(number_of_groups, Blockdata["main_node"], MPI_COMM_WORLD)
		total_stack      = bcast_number_to_all(total_stack,      Blockdata["main_node"], MPI_COMM_WORLD)
		
	else: #Tracker 1
		total_stack = len(input_list1)
		if Blockdata["myid"] == Blockdata["main_node"]:
			write_text_file(input_list1, os.path.join(box_dir, "previous_all_indexes.txt"))
		
		number_of_groups = total_stack//Tracker["current_img_per_grp"]				
		if number_of_groups<2: ERROR("The number_of_groups should be at least two", "depth_box_initialization", 1, Blockdata["myid"])
			
		new_assignment = create_nrandom_lists_from_given_pids(box_dir, os.path.join(box_dir, \
		      "previous_all_indexes.txt"), number_of_groups, 2)
	
	if Blockdata["myid"] == Blockdata["main_node"]:
		log_file.add('Sorting settings:  Number of images: %d.  Number of groups: %d.  Minimum group size: %d.'%(\
		      total_stack, number_of_groups, Tracker["constants"]["minimum_grp_size"]))
	else: Tracker = 0
	Tracker   = wrap_mpi_bcast(Tracker, Blockdata["main_node"], MPI_COMM_WORLD)
	#mpi_barrier(MPI_COMM_WORLD)
	return img_per_grp, number_of_groups, total_stack, minimum_grp_size, new_assignment
		
def output_iter_results(box_dir, ncluster, NACC, NUACC, minimum_grp_size, list_of_stable, unaccounted_list, img_per_grp, log_main):
	import numpy as np
	### single node
	iter_cutoff = min(minimum_grp_size, int(img_per_grp*0.75))
	new_list    = []
	nc          = 0
	NACC        = 0
	try:
		fout             = open(os.path.join(box_dir,"freq_cutoff.json"),'r')
		freq_cutoff_dict = convert_json_fromunicode(json.load(fout))
		fout.close()
	except: freq_cutoff_dict = {}
	
	for index_of_any in range(len(list_of_stable)):
		any = np.sort(list_of_stable[index_of_any])
		any.tolist()
		if len(any) >= iter_cutoff:
			new_list.append(any)
			write_text_file(any, os.path.join(box_dir, "Cluster_%03d.txt"%ncluster))
			freq_cutoff_dict["Cluster_%03d.txt"%ncluster] = Tracker["freq_fsc143_cutoff"]
			ncluster += 1
			nc       += 1
			NACC +=len(any)
		else:
			for im in range(len(any)):unaccounted_list.append(any[im])
	unaccounted_list = sorted(unaccounted_list)
	NUACC = len(unaccounted_list)
	fout  = open(os.path.join(box_dir, "freq_cutoff.json"),'w')
	json.dump(freq_cutoff_dict, fout)
	fout.close()
	return ncluster, NACC, NUACC, unaccounted_list, nc
	
def check_state_within_box_run(keepgoing, nruns, unaccounted_list, no_cluster_last_run, log_main):
	global Tracker, Blockdata
	total_stack             = len(unaccounted_list)
	number_of_groups        = total_stack//Tracker["current_img_per_grp"]
	if number_of_groups>=2: keepgoing = 1
	else:                   keepgoing = 0
	if keepgoing ==1: nruns +=1
	else:
		total_stack      = 0
		number_of_groups = 0
	if no_cluster_last_run:  number_of_groups -=1
	if number_of_groups<=1:  keepgoing =0
	if nruns> min(Tracker["constants"]["total_stack"]\
	     //Tracker["constants"]["img_per_grp"], 5):
	     	keepgoing = 0
	     	log_main.add("Sort3d stops because the number of runs exceeds desired number of iterations")
	return keepgoing, nruns, total_stack, number_of_groups
	
def get_box_partition(box_dir, ncluster, unaccounted_list):
	unaccounted_list = sorted(unaccounted_list)
	if ncluster >=1:
		clusters_in_box = [read_text_file(os.path.join(box_dir, \
		      "Cluster_%03d.txt"%ic)) for ic in range(ncluster)]
		if len(unaccounted_list)>0: clusters_in_box.append(unaccounted_list)
		alist, plist = merge_classes_into_partition_list(clusters_in_box)
	else:  plist = []
	return plist
				
def output_clusters(output_dir, partition, unaccounted_list, not_include_unaccounted, log_main):
	global Tracker, Blockdata
	import copy
	### single cpu function
	nclasses, npart = split_partition_into_ordered_clusters(partition)
	nc              = 0
	identified_clusters = []
	for ic in range(len(nclasses)):
		if len(nclasses[ic])>= max(Tracker["constants"]["img_per_grp"]\
		    //2**max(Tracker["constants"]["img_per_grp_split_rate"],1), \
		      min(Tracker["constants"]["minimum_grp_size"], Tracker["constants"]["img_per_grp"]*0.75)):
			write_text_file(nclasses[ic], os.path.join(output_dir,"Cluster_%03d.txt"%nc))
			nc +=1
			identified_clusters.append(nclasses[ic])
		else: unaccounted_list +=nclasses[ic]
	if len(unaccounted_list)>1:
		unaccounted_list = sorted(unaccounted_list)
		write_text_file(unaccounted_list, os.path.join(output_dir, "Core_set.txt"))
	nclasses = copy.deepcopy(identified_clusters)
	del identified_clusters
	if len(unaccounted_list)>1: # output unaccounted as the last cluster
		if not not_include_unaccounted:
			write_text_file(unaccounted_list, os.path.join(output_dir,"Cluster_%03d.txt"%nc))
	if not not_include_unaccounted: 
		import copy
		unclasses = copy.deepcopy(nclasses)
		unclasses.append(unaccounted_list)
		alist, partition = merge_classes_into_partition_list(unclasses)
		del unclasses
	else: alist, partition = merge_classes_into_partition_list(nclasses)
	write_text_row(partition, os.path.join(output_dir, "final_partition.txt"))
	return nclasses, nc, len(unaccounted_list)
	
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
	
def check_sorting_state(current_dir, keepchecking, log_file):
	import json
	try:
		fout = open(os.path.join(current_dir, "state.json"),'r')
		current_state = convert_json_fromunicode(json.load(fout))
		fout.close()
		if current_state["done"]:keepchecking = 1
		else: keepchecking = 0
	except:   keepchecking = 0
	return    keepchecking
	
def read_tracker_mpi(current_dir):
	global Tracker, Blockdata
	open_tracker = 1
	if(Blockdata["myid"] == Blockdata["main_node"]):
		try:
			fout    = open(os.path.join(current_dir, "Tracker.json"),'r')
			Tracker = convert_json_fromunicode(json.load(fout))
			fout.close()
		except:
			open_tracker = 0
	else: open_tracker = 0
	open_tracker =  bcast_number_to_all(open_tracker, Blockdata["main_node"], MPI_COMM_WORLD)
	if open_tracker ==1:
		if(Blockdata["myid"] != Blockdata["main_node"]): Tracker = 0
		Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], MPI_COMM_WORLD)
	else: ERROR("Failed to load Tracker", "%s"%current_dir, 1, Blockdata["myid"])
	return

def mark_sorting_state(current_dir, sorting_done, log_file):
	# single processor job
	import json
	current_state = {}
	fout = open(os.path.join(current_dir, "state.json"),'w')
	if sorting_done: current_state["done"] = True
	else:            current_state["done"] = False
	json.dump(current_state, fout)
	fout.close()
	return
	
def mem_calc_and_output_info(smearing_file, nxinit, iter_id_init_file, log_main):
	global Blockdata, Tracker
	from utilities    import read_text_file
	from applications import MPI_start_end
	import numpy as np
	# single CPU ; no info update
	precalculated_images_per_cpu = np.full(Blockdata["nproc"], 0, dtype=np.int32)
	smearing_list = np.array(read_text_file(smearing_file), dtype=np.int32)
	indx_list     = read_text_file(iter_id_init_file, -1)
	if len(indx_list) ==1: indx_list= indx_list[0]
	else:                  indx_list= indx_list[1]
	indx_list = np.sort(np.array(indx_list, dtype = np.int32))
	smearing_list = smearing_list[indx_list]
	avg_smear = np.sum(smearing_list)/smearing_list.shape[0]
	cdata_in_core  = (Tracker["total_stack"]*nxinit*nxinit*4.0)/1.e9/Blockdata["no_of_groups"]
	srdata_in_core = (nxinit*nxinit*np.sum(smearing_list)*4.)/1.e9/Blockdata["no_of_groups"]
	if not Tracker["constants"]["focus3D"]:	fdata_in_core = 0.0
	else: fdata_in_core = cdata_in_core
	ctfdata     = cdata_in_core
	refvol_size = (nxinit*nxinit*nxinit*4.0*2)/1.e9*Blockdata["no_of_processes_per_group"]# including the 3D mask
	log_main.add( "Precalculated data (GB) in core per node (available memory per node: %6.1f):"%Tracker["constants"]["memory_per_node"])
	log_main.add( "Images for comparison: %6.1f GB; shifted images: %6.1f GB; focus images: %6.1f GB; ctfs: %6.1f GB"%\
			(cdata_in_core, srdata_in_core, fdata_in_core, ctfdata))
	tdata = cdata_in_core+srdata_in_core+ctfdata+refvol_size+fdata_in_core
	if tdata/Tracker["constants"]["memory_per_node"]*100.> 90.:
		if not Tracker["constants"]["compute_on_the_fly"]: 
			ERROR("More than 90% memory is used. Please turn on compute_on_the_fly and rerun the program", \
			   "mem_calc_and_output_info", 1, Blockdata["myid"])
	log_main.add("The data kept in core for sorting will occupy %7.1f percent of memory;  average smearing is %5.1f"%\
		(min(tdata/Tracker["constants"]["memory_per_node"]*100., 90.0), avg_smear))

	log_main.add("Estimated per-node memory required for sorting:")
	smearings_on_nodes = np.full(Blockdata["no_of_groups"], 0.0, dtype=np.float32)
	smearings_per_cpu  = [None for im in range(Blockdata["nproc"])]
	dict = [None for i in range(Blockdata["nproc"])]
	for iproc in range(Blockdata["nproc"]):
		image_start, image_end = MPI_start_end(smearing_list.shape[0], Blockdata["nproc"],iproc)
		smearings_on_nodes[iproc//Blockdata["no_of_processes_per_group"]] += \
			np.sum(smearing_list[image_start:image_end])*(nxinit*nxinit*4.)/1.e9
		smearings_per_cpu[iproc] = smearing_list[image_start:image_end]
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
	if Tracker["constants"]["compute_on_the_fly"]:
		overhead_mem = 5. #
		mem_leftover = max(Tracker["constants"]["memory_per_node"], 32.) - \
		  overhead_mem - 2.*cdata_in_core - fdata_in_core - refvol_size
		mem_leftover = mem_leftover/float(Blockdata["no_of_processes_per_group"])
		log_main.add("Number of particles precomputed on each CPU:")
		log_main.add("Node       Number          Pecentage")
		msg = ""
		for im in range(Blockdata["nproc"]):
			image_start, image_end = MPI_start_end(smearing_list.shape[0], Blockdata["nproc"],im)
			size  = image_end - image_start 
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
		if Blockdata["nproc"]%3!=0:log_main.add(msg)
	else:
		for im in range(Blockdata["nproc"]):dict[im] = -1
	return dict
	
def depth_clustering_box(work_dir, input_accounted_file, input_unaccounted_file, \
         params, previous_params, nbox, log_main):
	global Tracker, Blockdata
	from utilities import read_text_file, wrap_mpi_bcast, write_text_file, write_text_row
	import copy, shutil
	from   shutil import copyfile
	from   math   import sqrt
	import numpy as np
	box_start      = time()
	acc_time       = 0.0
	box_niter      = 5
	no_groups_runs = 0
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
	Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], MPI_COMM_WORLD)	
	### ------- Initialization
	ncluster  = 0
	nruns     = 0
	keepgoing = 1
	converged = 0
	####
	img_per_grp, number_of_groups_init, total_stack_init, minimum_grp_size_init, new_assignment = \
	   depth_box_initialization(work_dir, input_accounted_file, input_unaccounted_file, log_main)
	NUACC = total_stack_init
	NACC  = 0
	####==========--------------------------
	assignment_list          = new_assignment[:]
	total_stack              = total_stack_init
	current_number_of_groups = number_of_groups_init
	unaccounted_list         = new_assignment[:]
	
	while (keepgoing == 1):
		if(Blockdata["myid"] == Blockdata["main_node"]):
			log_main.add(' Box %d is running for %10.1f minutes...'%(nbox, (time()-box_start)/60.))
		within_box_run_dir = os.path.join(work_dir, "run%d"%nruns)
		unaccounted_file   = os.path.join(within_box_run_dir, "Unaccounted_from_previous_run.txt")
		if(Blockdata["myid"] == Blockdata["main_node"]):
			time_box_start = time()
			os.mkdir(within_box_run_dir)
			write_text_file(unaccounted_list, unaccounted_file)# new starting point
		if nruns >0:
			assignment_list = create_nrandom_lists(unaccounted_file, current_number_of_groups, 2)
		if(Blockdata["myid"] == Blockdata["main_node"]):
			for indep in range(2):
					write_text_row(assignment_list[indep], os.path.join(within_box_run_dir,\
					 "independent_index_%03d.txt"%indep))
		#====>>>>Iter initialization<<<======
		iter                = 0
		previous_iter_ratio = 0.0
		current_iter_ratio  = 0.0
		iter_dir = os.path.join(within_box_run_dir, "iter%d"%iter)
		if Blockdata["myid"] == Blockdata["main_node"]:
			os.mkdir(iter_dir)
			for indep in range(2):
				write_text_row(assignment_list[indep], \
				    os.path.join(iter_dir, "random_assignment_%03d.txt"%indep))
		mpi_barrier(MPI_COMM_WORLD)
		
		iter_id_init_file = os.path.join(iter_dir, "random_assignment_000.txt")
		if Blockdata["myid"] == Blockdata["main_node"]:
			id_list          = read_text_file(os.path.join(within_box_run_dir, "independent_index_000.txt") , -1)
			total_stack, current_number_of_groups = len(id_list[0]), int((np.unique(np.array(id_list[0]))).shape[0])
		else:
			total_stack              = 0
			current_number_of_groups = 0
		total_stack                  = bcast_number_to_all(total_stack, Blockdata["main_node"], MPI_COMM_WORLD)
		current_number_of_groups     = bcast_number_to_all(current_number_of_groups, Blockdata["main_node"], MPI_COMM_WORLD)
		# read raw data
		original_data, norm_per_particle  = read_data_for_sorting(iter_id_init_file, params, previous_params)
		if Tracker["nosmearing"]: parameterstructure  = None
		else: parameterstructure  = read_paramstructure_for_sorting(\
		   iter_id_init_file,Tracker["paramstructure_dict"], Tracker["paramstructure_dir"])
		mpi_barrier(MPI_COMM_WORLD)
		
		# Estimate image size		
		Tracker["directory"] = within_box_run_dir
		Tracker["nxinit"], Tracker["freq_fsc143_cutoff"] = get_sorting_image_size(original_data, \
		      iter_id_init_file, current_number_of_groups, parameterstructure, norm_per_particle, log_main)
							
		if Blockdata["myid"] == Blockdata["main_node"]:
			log_main.add('Sorting cutoff frequency:  %5.2f'%(Tracker["freq_fsc143_cutoff"]) +\
			  '  The actual image size used in sorting:  %d'%(Tracker["nxinit"]))
		Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], MPI_COMM_WORLD)
			
		### preset variables in trackers
		Tracker["total_stack"]      = total_stack
		Tracker["number_of_groups"] = current_number_of_groups
		minimum_grp_size            = minimum_grp_size_init
		iter_previous_iter_ratio    = 0.0
		iter_current_iter_ratio     = 0.0
		if Tracker["constants"]["do_swap_au"]: swap_ratio = Tracker["constants"]["swap_ratio"]
		else:                                  swap_ratio = 0.0
		### prepare data for sorting: three datasets are generated cdata, fdata, srdata 
		cdata, rdata, fdata, ctf_images = downsize_data_for_sorting(original_data, preshift = True, \
		    npad = 1, norms = norm_per_particle)# pay attentions to shifts!
		del original_data
		# Compute data statistics
		
		if Blockdata["myid"] == Blockdata["main_node"]:
			precalculated_images_per_cpu = mem_calc_and_output_info(Tracker["constants"]["smearing_file"], \
			 Tracker["nxinit"], iter_id_init_file,log_main)
		else: precalculated_images_per_cpu = 0
		precalculated_images_per_cpu = wrap_mpi_bcast(precalculated_images_per_cpu, \
		     Blockdata["main_node"], MPI_COMM_WORLD)
		if Tracker["constants"]["compute_on_the_fly"]:
			Tracker["num_on_the_fly"] = precalculated_images_per_cpu
		srdata = precalculate_shifted_data_for_recons3D(rdata, parameterstructure, Tracker["refang"], \
	      Tracker["rshifts"], Tracker["delta"], Tracker["avgnorm"], Tracker["nxinit"], \
	      Tracker["constants"]["nnxo"], Tracker["nosmearing"], norm_per_particle, Tracker["constants"]["nsmear"])
		del rdata
		mpi_barrier(MPI_COMM_WORLD) # uneven smearing requires
		while((iter < box_niter) and (converged == 0)):
			max_iter = 0
			for indep_run_iter in range(2):
				Tracker["directory"] = os.path.join(iter_dir, "MGSKmeans_%03d"%indep_run_iter)
				MGSKmeans_index_file = os.path.join(iter_dir, "random_assignment_%03d.txt"%indep_run_iter)
				
				if Blockdata["myid"] == Blockdata["main_node"]:
					os.mkdir(Tracker["directory"])
					os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
				mpi_barrier(MPI_COMM_WORLD)
				tmp_final_list, premature, kmeans_iterations = \
				    Kmeans_minimum_group_size_orien_cones(cdata, fdata, srdata, ctf_images, \
					parameterstructure, norm_per_particle, \
				     MGSKmeans_index_file, params, minimum_grp_size, clean_volumes= Tracker["clean_volumes"] )
				max_iter = max(max_iter, kmeans_iterations)	   
				if Blockdata["myid"] == Blockdata["main_node"]:
					write_text_row(tmp_final_list, os.path.join(iter_dir,\
					    "partition_%03d.txt"%indep_run_iter))
			mpi_barrier(MPI_COMM_WORLD)
			if Blockdata["myid"] == Blockdata["main_node"]:
				minimum_grp_size1, maximum_grp_size1, list_of_stable, unaccounted_list, \
				      iter_current_iter_ratio, selected_number_of_groups, info_table = \
				do_withinbox_two_way_comparison(iter_dir, nbox, nruns, iter) # two partitions are written in partition_dir as partition_%03d.txt
			else: 
				unaccounted_list  =  0
				list_of_stable    =  0
				minimum_grp_size1 =  0
				Tracker           =  0
			minimum_grp_size1 = bcast_number_to_all(minimum_grp_size1, Blockdata["main_node"], MPI_COMM_WORLD)
			unaccounted_list  = wrap_mpi_bcast(unaccounted_list,       Blockdata["main_node"], MPI_COMM_WORLD)
			list_of_stable    = wrap_mpi_bcast(list_of_stable,         Blockdata["main_node"], MPI_COMM_WORLD)
			Tracker           = wrap_mpi_bcast(Tracker,                Blockdata["main_node"], MPI_COMM_WORLD)
			accounted_file    = os.path.join(iter_dir, "Accounted.txt")
			unaccounted_file  = os.path.join(iter_dir, "Core_set.txt")
			###### Compute rand_index for two independent runs
			if Blockdata["myid"] == Blockdata["main_node"]:
				core1      = read_text_file(os.path.join(iter_dir, "partition_000.txt"),  -1)
				core2      = read_text_file(os.path.join(iter_dir, "partition_001.txt"),  -1)
				core1      = core1[0]
				core2      = core2[0]
			else:
				core1 = 0
				core2 = 0
			core1  = wrap_mpi_bcast(core1, Blockdata["main_node"])
			core2  = wrap_mpi_bcast(core2, Blockdata["main_node"])
			rand_index  = compute_rand_index_mpi(core1, core2)
			if Blockdata["myid"] == Blockdata["main_node"]:
				info_table.append('  Rand index of two paritions is %5.4f '%rand_index)
				info_table.append('----------------------------------------------------------------------------------------------------------------')
				if abs(iter_current_iter_ratio - iter_previous_iter_ratio< 1.0) and \
				   iter_current_iter_ratio > 90.:converged = 1
				else: converged = 0
				if max_iter <=2: converged = 1 # Kmeans goes for only one step.
			else: converged = 0
			converged = bcast_number_to_all(converged, Blockdata["main_node"], MPI_COMM_WORLD)
			if (converged == 0) and (iter < box_niter):
				new_assignment_list = []
				for indep in range(2):
					tmp_list = swap_accounted_with_unaccounted_elements_mpi(accounted_file, unaccounted_file, \
						log_main, current_number_of_groups, swap_ratio)
					new_assignment_list.append(tmp_list)
				iter_previous_iter_ratio = iter_current_iter_ratio
				iter    +=1
				iter_dir = os.path.join(within_box_run_dir, "iter%d"%iter)
				if Blockdata["myid"] == Blockdata["main_node"]:
					os.mkdir(iter_dir)
					for indep in range(2):
						write_text_file(new_assignment_list[indep], \
							os.path.join(iter_dir, "random_assignment_%03d.txt"%indep))
			else:break
		del cdata
		del srdata
		if Tracker["constants"]["focus3D"]: del fdata
		del ctf_images
		if not Tracker["nosmearing"]: del parameterstructure
		del assignment_list
		if Blockdata["myid"] == Blockdata["main_node"]:
			ncluster, NACC, NUACC, unaccounted_list, new_clusters = \
			  output_iter_results(work_dir, ncluster, NACC, NUACC, Tracker["current_img_per_grp"],\
			    list_of_stable, unaccounted_list, Tracker["constants"]["img_per_grp"], log_main)
			if os.path.exists(os.path.join(work_dir, 'run%d'%nruns, 'tempdir')):
				shutil.rmtree(os.path.join(work_dir, 'run%d'%nruns, 'tempdir'))
		else:
			ncluster = 0
			NACC     = 0
			NUACC    = 0
			unaccounted_list = 0
			new_clusters     = 0
			Tracker          = 0
		new_clusters     = bcast_number_to_all(new_clusters,  Blockdata["main_node"], MPI_COMM_WORLD)
		ncluster         = bcast_number_to_all(ncluster,      Blockdata["main_node"], MPI_COMM_WORLD)
		NACC             = bcast_number_to_all(NACC,          Blockdata["main_node"], MPI_COMM_WORLD)
		NUACC            = bcast_number_to_all(NUACC,         Blockdata["main_node"], MPI_COMM_WORLD)
		unaccounted_list = wrap_mpi_bcast(unaccounted_list,   Blockdata["main_node"], MPI_COMM_WORLD)
		Tracker          = wrap_mpi_bcast(Tracker,            Blockdata["main_node"], MPI_COMM_WORLD)
		if new_clusters >0: 
			no_cluster     = False
			no_groups_runs = 0
		else:               
			no_cluster      = True
			no_groups_runs +=1
		keepgoing, nruns, total_stack, current_number_of_groups = \
		   check_state_within_box_run(keepgoing, nruns, unaccounted_list, no_cluster, log_main)
		   
		if Blockdata["myid"] == Blockdata["main_node"]:# report current state
			if new_clusters > 0:
				log_main.add(' ')
				log_main.add('In phase %d, the program found %d groups.'%(nruns, new_clusters))
				for itable in range(len(info_table)): log_main.add(info_table[itable])
		if no_groups_runs >= 2: 
			bad_clustering = 1
			break
	
	if Blockdata["myid"] == Blockdata["main_node"]:
		partition = get_box_partition(work_dir, ncluster, unaccounted_list)
	else: partition = 0
	partition = wrap_mpi_bcast(partition, Blockdata["main_node"], MPI_COMM_WORLD)
	
	if(Blockdata["myid"] == Blockdata["main_node"]):
		if ncluster > 0:
			if(ncluster == 1):
				log_main.add('In independent run  %d, the program found %d group.'%(nbox, ncluster))
			else: log_main.add('In independent run  %d, the program found %d groups'%(nbox, ncluster))
			bad_clustering = 0
		else:
			log_main.add('In independent run  %d, there were no reproducible groups found,'%nbox)
			bad_clustering = 1
			partition      =['']
		log_main.add('----------------------------------------------------------------------------------------------------------------')
		write_text_row(partition, os.path.join(work_dir, "partition.txt"))
	else: 
		bad_clustering = 0
		Tracker        = 0
	bad_clustering = bcast_number_to_all(bad_clustering, Blockdata["main_node"], MPI_COMM_WORLD)
	Tracker        = wrap_mpi_bcast(Tracker,             Blockdata["main_node"], MPI_COMM_WORLD)
	return bad_clustering

def check_mpi_settings(log_main):
	global Tracker, Blockdata
	from utilities    import wrap_mpi_bcast, read_text_file, bcast_number_to_all
	from fundamentals import smallprime
	from statistics   import scale_fsc_datasetsize
	import os
	import numpy as np
	number_of_groups = float(Tracker["constants"]["total_stack"])/float(Tracker["constants"]["img_per_grp"])
	avg_fsc = scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], \
	   float(Tracker["constants"]["total_stack"]), \
	   float(Tracker["constants"]["total_stack"])//number_of_groups)
	fsc143 = len(avg_fsc)
	for ifreq in range(len(avg_fsc)):
		if avg_fsc[ifreq] < 0.143:
			fsc143 = ifreq -1
			break
	if fsc143 !=0: 
		nxinit = smallprime(min((int(fsc143)+ max(int(Tracker["constants"]["nnxo"]*0.03), 5))*2, \
		     Tracker["constants"]["nnxo"]))
		nxinit = smallprime(nxinit-nxinit%2)
		nxinit = nxinit - nxinit%2
	else: ERROR("Incorrect fsc143", "check_mpi_settings", 1, Blockdata["myid"])
	### pre-calculated 2D data: 
	#1. Cdata(data for comparison)
	#2. Rdata(data for rec3D) # smearing could dramatically increase data size
	#3. fdata (focus projection); 
	#4. 2D CTFs 
	# other large data
	#1. reference volumes
	#2. rec3D volumes 2 times padded  
	sys_required_mem = 0.3*Blockdata["no_of_processes_per_group"]# Overhead
	if( Blockdata["myid"] == Blockdata["main_node"]):
		log_main.add('\n')
		log_main.add('----------------------------------------------------------------------------------------------------------------' )
		log_main.add('                 =======>     Number of input images and memory information     <=====')
		log_main.add('----------------------------------------------------------------------------------------------------------------' )
		log_main.add("Number of processes: %d number of nodes:  %d.  Number of processes per node:  %d."%(Blockdata["nproc"], Blockdata["no_of_groups"], Blockdata["no_of_processes_per_group"]))
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
		ERROR("Initial info is not provided", "check_mpi_settings", 1, Blockdata["myid"])
	try:
		mem_bytes = os.sysconf('SC_PAGE_SIZE')*os.sysconf('SC_PHYS_PAGES')# e.g. 4015976448
		mem_gib   = mem_bytes/(1024.**3) # e.g. 3.74
		if( Blockdata["myid"] == Blockdata["main_node"]):
			log_main.add("Available memory information provided by the operating system: %5.1f GB"%mem_gib)
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
			log_main.add("Memory per node: %f"%Tracker["constants"]["memory_per_node"])
			
	if(Blockdata["myid"] == Blockdata["main_node"]):
		log_main.add("Total number of images: %d.  Number of images per group: %d."%\
		   (Tracker["constants"]["total_stack"], Tracker["constants"]["img_per_grp"]))
		
	if(Blockdata["myid"] == Blockdata["main_node"]):
		log_main.add("The total available memory per node:               %5.1f GB"%total_memory)
		log_main.add("The size of total input 2D stack:                  %5.1f GB"%raw_data_size)
		log_main.add("The per-node amount of memory 2D data will occupy: %5.1f GB"%raw_data_size_per_node)
		log_main.add("Precalculated 2D data (no smearing) for sorting:   %5.1f GB"%sorting_data_size_per_node)
		log_main.add("Reserved overhead memory per node:                 %5.1f GB"%sys_required_mem)
		
	if (total_memory - sys_required_mem - sorting_data_size_per_node ) <0.0:
		ERROR("Insufficient memory to pass the sorting processs. Increase a node", "check_mpi_settings", 1, Blockdata["myid"])
		
	if (total_memory - sys_required_mem - raw_data_size_per_node - 2*raw_data_size_per_node*ratio**2 ) <0.0: 
		ERROR("Insufficient memory to read in the raw data to produce downsized data. Increase a node", \
		   "check_mpi_settings", 1, Blockdata["myid"])
			
	images_per_cpu = int(float(Tracker["constants"]["total_stack"])/Blockdata["nproc"])
	images_per_cpu_for_unaccounted_data  = Tracker["constants"]["img_per_grp"]*1.5/float(Blockdata["nproc"])
	if( Blockdata["myid"] == Blockdata["main_node"]):
		log_main.add("Number of images per processor:                    %5d"%images_per_cpu )
	if( images_per_cpu < 5):
		ERROR("Number of images per processor is less than 5", \
		   "one may want to consider decreasing the number of processors used in MPI setting", \
		   0, Blockdata["myid"])
		   
	if(Blockdata["myid"] == Blockdata["main_node"]):
		avg_smearing_on_node = sum(read_text_file(Tracker["constants"]["smearing_file"]))//Blockdata["no_of_groups"]
	else: avg_smearing_on_node = 0
	avg_smearing_on_node  = bcast_number_to_all(avg_smearing_on_node, Blockdata["main_node"], MPI_COMM_WORLD)
	avg_images_on_node    = Tracker["constants"]["total_stack"]//Blockdata["no_of_groups"]
	precalculated_2D_data = raw_data_size_per_node*ratio**2*(max(avg_smearing_on_node/avg_images_on_node,1.) + 3.)
	if ((precalculated_2D_data - sys_required_mem)>total_memory) and (not Tracker["constants"]["compute_on_the_fly"]):
		ERROR("The size of precalculated shifted 2D data is too large. Turn on option compute_on_the_fly", \
		  "check_mpi_settings", 1, Blockdata["myid"])
	else:
		if(Blockdata["myid"] == Blockdata["main_node"]):
			log_main.add("Percentage of memory occupied by precalculated shifted 2D data: %6.1f%%"%\
			  (precalculated_2D_data/total_memory*100.))
			if Tracker["constants"]["compute_on_the_fly"]:
				log_main.add("The compute_on_the_fly option is on and 2D data is precacluated until 90 percent of memory is filled up")
	if(Blockdata["myid"] == Blockdata["main_node"]):
		log_main.add('----------------------------------------------------------------------------------------------------------------\n')
	return
	
def get_sorting_image_size(original_data, partids, number_of_groups, sparamstructure, snorm_per_particle, log_main):
	global Tracker, Blockdata
	from utilities    import wrap_mpi_bcast, read_text_file, write_text_file
	from fundamentals import smallprime
	from applications import MPI_start_end
	import numpy as np		
	if Tracker["compute_fsc"]:
		iter = 0
		Tracker["number_of_groups"] = number_of_groups
		if(Blockdata["myid"] == Blockdata["main_node"]):
			log_main.add("3D reconstruction is computed using window size:  %d"%Tracker["nxinit_refinement"])
			lpartids = read_text_file(partids, -1)
			if len(lpartids) == 1:
				iter_assignment = []
				for im in range(len(lpartids[0])):
					iter_assignment.append(randint(0,number_of_groups - 1))# simple version
			else:
				iter_assignment = lpartids[0]
		else:   iter_assignment = 0
		iter_assignment = wrap_mpi_bcast(iter_assignment, Blockdata["main_node"])
	
		Tracker["total_stack"] = len(iter_assignment)
		proc_list = [[None, None] for iproc in range(Blockdata["nproc"])]
		for iproc in range(Blockdata["nproc"]):
			iproc_image_start, iproc_image_end = MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], iproc)
			proc_list[iproc] = [iproc_image_start, iproc_image_end]
	
		compute_noise(Tracker["nxinit_refinement"])
		rdata = downsize_data_for_rec3D(original_data, Tracker["nxinit_refinement"], return_real = False, npad = 1)
		
		for im in range(len(rdata)):# update rdata
			rdata[im].set_attr("group", np.array(iter_assignment[proc_list[Blockdata["myid"]][0]:\
			    proc_list[Blockdata["myid"]][1]], dtype=np.int64)[im])

		Tracker["nxinit"] = Tracker["nxinit_refinement"]
		compute_noise(Tracker["nxinit"])
		do3d_sorting_groups_fsc_only_iter(rdata, sparamstructure, snorm_per_particle, iteration = iter)
		del rdata
		if( Blockdata["myid"] == Blockdata["main_node"]):
			fsc_data = []
			for igroup in range(Tracker["number_of_groups"]):
				for ichunk in range(2):
					tmp_fsc_data = read_text_file(os.path.join(Tracker["directory"], \
					     "fsc_driver_chunk%d_grp%03d_iter%03d.txt"%(ichunk, igroup, iter)), -1)
					fsc_data.append(tmp_fsc_data[0])
		else: fsc_data = 0
		fsc_data = wrap_mpi_bcast(fsc_data, Blockdata["main_node"])
		avg_fsc = [ 0.0 for i in range(len(fsc_data[0]))]
		avg_fsc[0] = 1.0
		for igroup in range(1): # Use group zero first
			for ifreq in range(1, len(fsc_data[0])):avg_fsc[ifreq] += fsc_data[igroup][ifreq]
		fsc143 = len(fsc_data[0])
	else:
		from statistics import scale_fsc_datasetsize
		avg_fsc = scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], \
		     float(Tracker["constants"]["total_stack"]), Tracker["total_stack"]//number_of_groups)
		fsc143 = len(avg_fsc)
		
	for ifreq in range(len(avg_fsc)):
		if avg_fsc[ifreq] < 0.143:
			fsc143 = ifreq -1
			break
	if fsc143 !=0: 
		nxinit = smallprime(min((int(fsc143)+ max(int(Tracker["constants"]["nnxo"]*0.03), 5))*2, \
		     Tracker["constants"]["nnxo"]))
		nxinit = smallprime(nxinit-nxinit%2)
		nxinit = nxinit - nxinit%2
	else: ERROR("Incorrect image size encountered", "get_sorting_image_size", 1, Blockdata["myid"])
	
	if Tracker["nosmearing"] and Tracker["constants"]["nxinit"] !=-1:
		nxinit = Tracker["constants"]["nxinit"]
		if nxinit >Tracker["constants"]["nnxo"]:
			ERROR("User-provided nxinit is incorrect, it is larger than the original image size", "get_sorting_image_size", 1, Blockdata["myid"])
			
	Tracker["nxinit"] = nxinit
	compute_noise(Tracker["nxinit"])
	freq_fsc143_cutoff = float(fsc143)/float(nxinit)
	if(Blockdata["myid"] == Blockdata["main_node"]): 
		write_text_file(avg_fsc, os.path.join(Tracker["directory"], "fsc_image_size.txt"))
	if Tracker["compute_fsc"]:
		del iter_assignment
		del proc_list
		del fsc_data
	del avg_fsc
	return nxinit, freq_fsc143_cutoff
	
def compute_noise(image_size):
	global Tracker, Blockdata
	from utilities    import get_im, model_blank
	from fundamentals import fft
	if Tracker["applybckgnoise"]: # from SPARX refinement only
		if(Blockdata["myid"] == Blockdata["main_node"]):
			tsd = get_im(Tracker["bckgnoise"]) # inverted power spectrum
			nnx = tsd.get_xsize()
			nny = tsd.get_ysize()
		else:
			nnx = 0
			nny = 0
		nnx = bcast_number_to_all(nnx, Blockdata["main_node"], MPI_COMM_WORLD)
		nny = bcast_number_to_all(nny, Blockdata["main_node"], MPI_COMM_WORLD)
		if( Blockdata["myid"] != Blockdata["main_node"]):tsd = model_blank(nnx, nny)
		bcast_EMData_to_all(tsd, Blockdata["myid"], Blockdata["main_node"])
		temp_image = model_blank(image_size, image_size)
		temp_image = fft(temp_image)
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
			Blockdata["unrolldata"].append(Util.unroll1dpw(ny, ny, Blockdata["bckgnoise"][im]))
	else: # from datastack and relion
		temp_image = model_blank(image_size, image_size)
		temp_image = fft(temp_image)
		nx = temp_image.get_xsize()
		ny = temp_image.get_ysize()
		Blockdata["bckgnoise"]  = [1.0]*nx
		Blockdata["unrolldata"] = Util.unroll1dpw(ny, ny, nx*[1.0])
	return
			
def get_params_for_analysis(orgstack, ali3d_params, smearing_file, smearing_number):
	if ali3d_params is not None:
		ali3d_params   = read_text_file(ali3d_params, -1)
		try: norm_list = ali3d_params[7]
		except: norm_list = None
	if( orgstack is not None ):
		defo_list = []
		ctfs      = EMUtil.get_all_attributes(orgstack, "ctf")
		for im in range(len(ctfs)): defo_list.append(ctfs[im].defocus)
	else: defo_list = None
	if smearing_file is not None:
		try: smearing_list = read_text_file(smearing_file)
		except: smearing_list = None
	else: smearing_list = None
	return defo_list, smearing_list, norm_list

def do_one_way_anova_scipy(clusters, value_list, name_of_variable="variable", log_main = None):
	# single cpu program
	import copy
	from math   import sqrt
	import scipy
	from scipy import stats
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
	
	if   K==2: res = stats.f_oneway(x0, x1)
	elif K==3: res = stats.f_oneway(x0, x1, x2)
	elif K==4: res = stats.f_oneway(x0, x1, x2, x3)
	elif K==5: res = stats.f_oneway(x0, x1, x2, x3, x4)
	elif K==6: res = stats.f_oneway(x0, x1, x2, x3, x4, x5)
	elif K==7: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6)
	elif K==8: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7)
	elif K==9: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8)
	elif K==10: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9)
	elif K==11: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
	elif K==12: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11)
	elif K==13: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12)
	elif K==14: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13)
	elif K==15: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14)
	elif K==16: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15)
	elif K==17: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16)
	elif K==18: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17)
	elif K==19: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18)
	elif K==20: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19)
	elif K==21: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19, x20)
	elif K==22: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19, x20, x21)
	elif K==23: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19, x20, x21, x22)
	elif K==24: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19, x20, x21, x22, x23)
	elif K==25: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19, x20, x21, x22, x23, x24)
	elif K==26: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19, x20, x21, x22, x23, x24, x25)
	elif K==27: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19, x20, x21, x22, x23, x24, x25, x26)
	elif K==28: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19, x20, x21, x22, x23, x24, x25, x26, x27)
	elif K==29: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28)
	elif K==30: res = stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29)
	else:
		return None, None, None
	res_table_stat = []
	for im in range(K):
		alist = table_stat(replicas[im])
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
		std_list.append(sqrt(std_per_group/float(len(replicas[i]))))
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
	from math import sqrt
	for i in range(K):
		log_main.add('{:5} {:^7d} {:^8d} {:12.4f} {:12.4f}'.format('ANOVA', i, len(replicas[i]), res_table_stat[i][0], sqrt(res_table_stat[i][1])))
	log_main.add(' ')
	log_main.add('ANOVA  Pair-wise tests')
	log_main.add('{:5} {:^3} {:^3} {:^12} {:^12} {:^12} {:^12}'.format('ANOVA', 'A', 'B', 'avgA','avgB', 'F_value', 'Significance')) 
	for ires in range(K-1):
		for jres in range(ires+1, K):
			cres = stats.f_oneway(replicas[ires], replicas[jres])
			log_main.add('{:5} {:^3d} {:^3d} {:12.4f} {:12.4f} {:12.3f} {:12.4f} '.format('ANOVA', ires, jres, avgs[ires], avgs[jres], cres[0], cres[1]*100.))
	log_main.add(' ')
	log_main.add('================================================================================================================\n')
	return res[0], res[1]

def check_3dmask(log_main):
	global Tracker, Blockdata
	###########################################################################	
	Tracker["nxinit"]     = Tracker["nxinit_refinement"]
	#   shrinkage, current resolution, fuse_freq
	Tracker["total_stack"] = Tracker["constants"]["total_stack"]
	Tracker["shrinkage"]   = float(Tracker["nxinit"])/Tracker["constants"]["nnxo"]
	Tracker["radius"]      = Tracker["constants"]["radius"]*Tracker["shrinkage"]
	try: fuse_freq = Tracker["fuse_freq"]
	except: Tracker["fuse_freq"] = int(Tracker["constants"]["pixel_size"]\
	     *Tracker["constants"]["nnxo"]/Tracker["constants"]["fuse_freq"]+0.5)
	Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"], MPI_COMM_WORLD)
	dump_tracker(Tracker["constants"]["masterdir"])
	Tracker["shrinkage"] = float(Tracker["nxinit"])/Tracker["constants"]["nnxo"]
	if Tracker["constants"]["focus3D"]:
		if(Blockdata["myid"] == Blockdata["main_node"]):
			log_main.add(' ')
			log_main.add('Sort3d is run in focus mode, and focus mask is %s'%Tracker["constants"]["focus3D"]+'\n')
	return

def import_data(log_main, override):
	global Tracker, Blockdata
	# Two sorting scenarios
	# 1. import data and refinement parameters from meridien refinement;
	# 2. given data stack and xform.projection/ctf in header(For simulated test data);
	import_from_relion_refinement = 0
	import_from_sparx_refinement  = 0
	import_from_data_stack		  = 0
	total_stack					  = 0
	if Tracker["constants"]["refinement_method"] =="SPARX": # Senario one
		import_from_sparx_refinement = get_input_from_sparx_ref3d(log_main)
		Tracker["smearing"] = True
	else:  # Scenario three, sorting from a given data stack, general cases
		import_from_data_stack = get_input_from_datastack(log_main)
		Tracker["constants"]["hardmask"] = True
		Tracker["applybckgnoise"]        = False
		Tracker["applymask"]             = True
		Tracker["smearing"]              = False
	Tracker["total_stack"] = Tracker["constants"]["total_stack"]
	###=====<----------------====>checks=====-------------
	if Tracker["constants"]["symmetry"] != Tracker["constants"]["sym"]:
		if(Blockdata["myid"] == Blockdata["main_node"]):
			msg = "Input symmetry %s is altered to %s after reading refinement information! "%\
			       (Tracker["constants"]["sym"], Tracker["constants"]["symmetry"])
			log_main.add(msg)

	## checking settings!
	if Tracker["constants"]["total_stack"]//Tracker["constants"]["img_per_grp"]<=1: 
		ERROR("Requested number of images per group is too large, corresponds to requesting one group", "sxsort3d_depth.py", 1,  Blockdata["myid"])
	if Tracker["constants"]["minimum_grp_size"] > Tracker["constants"]["img_per_grp"]:
		ERROR("Minimum group size cannot be larger that the requested number of images per group", "sxsort3d_depth.py", 1,  Blockdata["myid"])
	if( not override):
		# Check sanity of parameters
		nnew = int(round(Tracker["constants"]["img_per_grp"]/4.0))
		if(Blockdata["myid"] == Blockdata["main_node"]):
			log_main.add("Suggested minimum number of images per group: %d"%nnew)
		if(Tracker["constants"]["minimum_grp_size"]<nnew): ERROR("Minimum number of images per group is too low", "import_data", 1, Blockdata["myid"])

	return
	
def create_masterdir():
	global Tracker, Blockdata
	masterdir = Tracker["constants"]["masterdir"]
	restart   = 0
	if(Blockdata["myid"] == Blockdata["main_node"]):
		if not os.path.exists(os.path.join(Tracker["constants"]["masterdir"])):
			if not masterdir:
				timestring = strftime("_%d_%b_%Y_%H_%M_%S", localtime())
				masterdir  ="sort3d"+timestring
				os.mkdir(masterdir)
			else:
				if not os.path.exists(masterdir): os.mkdir(masterdir)
			li =len(masterdir)
		else:
			li =len(masterdir)
			if os.path.exists(os.path.join(Tracker["constants"]["masterdir"], "generation_000")): restart = 1
			if os.path.exists(os.path.join(Tracker["constants"]["masterdir"], "Tracker.json")):   restart = 1
	else: 
		restart = 0
		li = 0
	restart   = bcast_number_to_all(restart,       Blockdata["main_node"], MPI_COMM_WORLD)
	li        = mpi_bcast(li,       1,   MPI_INT,  Blockdata["main_node"],MPI_COMM_WORLD)[0]
	masterdir = mpi_bcast(masterdir,li,  MPI_CHAR, Blockdata["main_node"],MPI_COMM_WORLD)
	masterdir = string.join(masterdir,"")
	if not Tracker["constants"]["masterdir"]: Tracker["constants"]["masterdir"]  = masterdir
	Tracker["constants"]["chunk_0"]  = os.path.join(Tracker["constants"]["masterdir"],"chunk_0.txt")
	Tracker["constants"]["chunk_1"]  = os.path.join(Tracker["constants"]["masterdir"],"chunk_1.txt")
	return restart

def sort3d_init(to_be_decided, log_main):
	global Tracker, Blockdata
	keepsorting = 1
	if Tracker["constants"]["img_per_grp"]<= 2:
		log_main.add("Number of images per group is too small.")
		keepsorting = 0
	if Tracker["total_stack"] <= Blockdata["nproc"]*Tracker["minimum_img_per_cpu"]:
		log_main.add("Either user requires too many processors, or number of images is too small.")
		keepsorting = 0
	return keepsorting

def print_shell_command(args_list, log_main):
	global Tracker, Blockdata
	if(Blockdata["myid"] == Blockdata["main_node"]):
		log_main.add("The shell line command:")
		line = ""
		for a in args_list: line +=(a + " ")
		log_main.add(line)
	mpi_barrier(MPI_COMM_WORLD)
	return
	
def AI_MGSKmeans(iter_assignment, last_iter_assignment, best_assignment, \
        keepgoing, best_score, stopercnt, minimum_grp_size, log_file):
	# single cpu function
	import numpy as np
	group_dict       = np.sort(np.unique(iter_assignment))
	number_of_groups = group_dict.shape[0]
	msg = "group id     group size"
	log_file.add(msg)
	clusters = []
	for igrp in range(group_dict.shape[0]):
		one_group = np.sort(iter_assignment[isin(iter_assignment, group_dict[igrp])])
		msg ="%5d    %10d"%(igrp, one_group.shape[0])
		log_file.add(msg)
		clusters.append(one_group.shape[0])
	is_unicorn_cluster = 0
	nc  = 0
	tot = 0
	for ic in range(len(clusters)):
		tot += clusters[ic]
		if clusters[ic] < minimum_grp_size + len(clusters): nc +=1
	if tot//minimum_grp_size>2*len(clusters) and nc+1==len(clusters):
		is_unicorn_cluster =1 # Unicorn condition
	if is_unicorn_cluster == 0:
		sum_newindices1 = 0
		sum_newindices2 = 0
		ratio, newindices, stable_clusters = \
		   compare_two_iterations(iter_assignment, last_iter_assignment)
		for idx in range(len(newindices)):
			sum_newindices1 += newindices[idx][0]
			sum_newindices2 += newindices[idx][1]
			if newindices[idx][0] != newindices[idx][1]:
				msg ="Group %d  swaps with group %d "%(newindices[idx][0], newindices[idx][1])
				log_file.add(msg)
		changed_nptls = 100.- ratio*100.
		if best_score >= changed_nptls:
			best_score = changed_nptls
			best_assignment = copy.copy(iter_assignment)
		if changed_nptls < stopercnt: keepgoing = 0
	else:
		msg ="Unicorn cluster is found. Shuffle assignment."
		log_file.add(msg)
		iter_assignment = np.random.randint(0,  number_of_groups, size=len(iter_assignment))
		best_score    = 100.
		changed_nptls = 100.
		keepgoing     = 1
		best_assignment = copy.copy(iter_assignment)
	return best_score, changed_nptls, keepgoing, best_assignment, iter_assignment
#####
def Kmeans_minimum_group_size_orien_cones(cdata, fdata, srdata, ctf_images,\
	paramstructure, norm_per_particle,\
    partids, params, minimum_group_size_init, clean_volumes = True):
	global Tracker, Blockdata
	import shutil
	import numpy as np
	### data type: 1. assignment np.int32; particle ids np.int32; EMAN2 functions python native int;
	
	#==========---------- >>>>EQKmeans initialization ==========------------ 
	log_main	             = Logger()
	log_main                 = Logger(BaseLogger_Files())
	log_main.prefix          = Tracker["directory"]+"/"
	premature                = 0
	changed_nptls            = 100.0
	number_of_groups         = Tracker["number_of_groups"]
	stopercnt                = Tracker["constants"]["stop_mgskmeans_percentage"]
	total_iter               = 0
	partial_rec3d            = False
	best_score               = 100.0
	best_assignment          = []
	max_iter                 = Tracker["total_number_of_iterations"]
	last_score               = 100.0
	fixed_value              = 100.0
	has_converged            = 0
	times_around_fixed_value = 0  
	###=====<------------mask3D
	if( Blockdata["myid"] == Blockdata["main_node"]):
		try:
			if os.path.exists(Tracker["constants"]["mask3D"]): # prepare mask
				mask3D = get_im(Tracker["constants"]["mask3D"])
				if mask3D.get_xsize() != Tracker["nxinit"]: mask3D = fdecimate(mask3D, \
				    Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"], True, False)
		except:
			mask3D = model_circle(Tracker["constants"]["radius"], Tracker["constants"]["nnxo"], \
			    Tracker["constants"]["nnxo"], Tracker["constants"]["nnxo"])
			mask3D = fdecimate(mask3D, Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"], True, False)
	else:
		mask3D = model_blank(Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"])
	bcast_EMData_to_all(mask3D, Blockdata["myid"], Blockdata["main_node"])
	### get ptl_ids and intial assignment
	if( Blockdata["myid"] == Blockdata["main_node"]):
		lpartids = read_text_file(partids, -1)
		if len(lpartids) == 1: iter_assignment = np.random.randint(0, number_of_groups, size=len(lpartids[0]), dtype=np.int32)
		else:                  iter_assignment = np.array(lpartids[0], dtype=np.int32)
	else: iter_assignment = 0
	iter_assignment                 = wrap_mpi_bcast(iter_assignment, Blockdata["main_node"]) # initial assignment
	Tracker["total_stack"]          = iter_assignment.shape[0]
	nima                            = len(cdata) # local size
	image_start, image_end          = MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"])
	Tracker["min_orien_cone_size"] = Tracker["number_of_groups"]*Tracker["minimum_ptl_number"]
	angle_step                      = number_of_cones_to_delta(Tracker["constants"]["orientation_cones"])
	ptls_in_orien_cones				= get_angle_step_and_orien_cones_mpi(params, partids, angle_step)
	minimum_group_size              = max(minimum_group_size_init, Tracker["number_of_groups"]*len(ptls_in_orien_cones))
	minimum_group_size_ratio        = min((minimum_group_size*Tracker["number_of_groups"])/float(Tracker["total_stack"]), 0.95)

	### printed info
	if( Blockdata["myid"] == Blockdata["main_node"]):
		log_main.add('----------------------------------------------------------------------------------------------------------------' )
		log_main.add(' ==========> MGSKmeans clustering <========== ')
		log_main.add('----------------------------------------------------------------------------------------------------------------' )
		msg = "Total_stack:  %d K = : %d  nxinit: %d  CTF:  %s  Symmetry:  %s  stop percentage: %f  3-D mask: %s focus mask: %s  Comparison method: %s  minimum_group_size: %d orien  %d"% \
		   (Tracker["total_stack"], Tracker["number_of_groups"], Tracker["nxinit"],  Tracker["constants"]["CTF"], Tracker["constants"]["symmetry"], \
		   stopercnt, Tracker["constants"]["mask3D"], Tracker["constants"]["focus3D"], Tracker["constants"]["comparison_method"], \
		       minimum_group_size, len(ptls_in_orien_cones))
		log_main.add(msg)
			
	proc_list = [[None, None] for iproc in range(Blockdata["nproc"])]
	for iproc in range(Blockdata["nproc"]):
		iproc_image_start, iproc_image_end = MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], iproc)
		proc_list[iproc] = [iproc_image_start, iproc_image_end]
		
	compute_noise(Tracker["nxinit"])
	mpi_barrier(MPI_COMM_WORLD)	
	last_iter_assignment = np.copy(iter_assignment)
	best_assignment      = np.copy(iter_assignment)
	total_iter           = 0
	keepgoing            = 1
	do_partial_rec3d     = 0
	partial_rec3d        = False
	### share memory preparation
	#disp_unit  = np.dtype("f4").itemsize
	#refvolsize = Tracker["nxinit"]*Tracker["nxinit"]*Tracker["nxinit"]
	#emnumpy1   = EMNumPy()
	#win_vol, base_vol = mpi_win_allocate_shared(refvolsize*disp_unit, disp_unit, \
	#   MPI_INFO_NULL, Blockdata["shared_comm"])
	#if(Blockdata["myid_on_node"]!= 0): base_vol, = mpi_win_shared_query(win_vol, MPI_PROC_NULL)
	#volbuf = np.frombuffer(np.core.multiarray.int_asbuffer(base_vol, refvolsize*disp_unit), dtype = 'f4')
	#volbuf = volbuf.reshape(Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"])
	#### ---end of shared memory
	while total_iter < max_iter:
		if(Blockdata["myid"] == Blockdata["main_node"]):
			msg = "Iteration %d particle assignment changed ratio  %f "%(total_iter, changed_nptls)
			log_main.add(msg)
			write_text_file(np.copy(iter_assignment).tolist(), os.path.join(Tracker["directory"],\
			     "assignment%03d.txt"%total_iter))
			if changed_nptls< 50.0: do_partial_rec3d = 1
			else:                   do_partial_rec3d = 0
		else: do_partial_rec3d = 0
		do_partial_rec3d       = bcast_number_to_all(do_partial_rec3d, Blockdata["main_node"], MPI_COMM_WORLD)
		if do_partial_rec3d ==1: partial_rec3d = True
		else:                    partial_rec3d = False
		update_data_assignment(cdata, srdata, iter_assignment, proc_list, Tracker["nosmearing"], Blockdata["myid"])
		mpi_barrier(MPI_COMM_WORLD)
		current_group_sizes = get_group_size_from_iter_assign(iter_assignment)
		do3d_sorting_groups_nofsc_smearing_iter(srdata, paramstructure, norm_per_particle,\
		  partial_rec3d, current_group_sizes, iteration = total_iter)
		local_peaks = np.full((number_of_groups, nima), -1.0, dtype = np.float32)
		## compute peaks and save them in 1D list
		### shared memor
		"""
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
		for iref in range(number_of_groups):
			if(Blockdata["myid"] == Blockdata["last_node"]):
				ref_vol = get_im(os.path.join(Tracker["directory"],"vol_grp%03d_iter%03d.hdf"%(iref, total_iter)))
				nnn = ref_vol.get_xsize()
				if(Tracker["nxinit"] != nnn): ref_vol = fdecimate(ref_vol, Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"], True, False)
				stat = Util.infomask(ref_vol, mask3D, False)
				ref_vol -= stat[0]
				if stat[1]!=0.0:Util.mul_scalar(ref_vol, 1.0/stat[1])
				ref_vol *=mask3D
			else: ref_vol = model_blank(Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"])
			bcast_EMData_to_all(ref_vol, Blockdata["myid"], Blockdata["last_node"])
			## Image comparison optimal solution is the larger one	
			if Tracker["constants"]["comparison_method"] =="cross": ref_peaks = compare_two_images_cross(cdata, ref_vol, ctf_images)
			else:                                                   ref_peaks = compare_two_images_eucd(cdata, ref_vol, fdata, ctf_images)
			local_peaks[iref] = ref_peaks
			mpi_barrier(MPI_COMM_WORLD)
		local_peaks = local_peaks.reshape(number_of_groups*nima)
		mpi_barrier(MPI_COMM_WORLD)
		# pass to main_node
		if Blockdata["myid"] == Blockdata["main_node"]:
			dmatrix = np.full((number_of_groups, Tracker["total_stack"]), 0.0, dtype=np.float32)
			for im in range(local_peaks.shape[0]): dmatrix[im//nima][im%nima + image_start] = local_peaks[im]
		else: dmatrix = 0
		if Blockdata["myid"] != Blockdata["main_node"]:
			wrap_mpi_send(local_peaks, Blockdata["main_node"], MPI_COMM_WORLD)
		else:
			for iproc in range(Blockdata["nproc"]):
				if iproc != Blockdata["main_node"]:
					local_peaks = wrap_mpi_recv(iproc, MPI_COMM_WORLD)
					iproc_nima  = proc_list[iproc][1] - proc_list[iproc][0]
					for im in range(len(local_peaks)): 
						dmatrix[im/iproc_nima][im%iproc_nima + proc_list[iproc][0]] = local_peaks[im]
		dmatrix = wrap_mpi_bcast(dmatrix, Blockdata["main_node"], MPI_COMM_WORLD)
		last_iter_assignment = np.copy(iter_assignment)
		iter_assignment      = np.full(Tracker["total_stack"], -1, dtype=np.int32)
	 	for iorien in range(len(ptls_in_orien_cones)):
			if iorien%Blockdata["nproc"] == Blockdata["myid"]:
				iter_assignment[ptls_in_orien_cones[iorien]] = \
				do_assignment_by_dmatrix_orien_group_minimum_group_size(\
				 dmatrix, ptls_in_orien_cones[iorien], Tracker["number_of_groups"], \
				   minimum_group_size_ratio)
		mpi_barrier(MPI_COMM_WORLD)
		if Blockdata["myid"] != Blockdata["main_node"]: 
			wrap_mpi_send(iter_assignment, Blockdata["main_node"], MPI_COMM_WORLD)
		else:
			for iproc in range(Blockdata["nproc"]):
				if iproc != Blockdata["main_node"]:
					dummy = wrap_mpi_recv(iproc, MPI_COMM_WORLD)
					iter_assignment[np.where(dummy>-1)] = dummy[np.where(dummy>-1)]
		mpi_barrier(MPI_COMM_WORLD)
		iter_assignment = wrap_mpi_bcast(iter_assignment, Blockdata["main_node"], MPI_COMM_WORLD)
		if Blockdata["myid"] == Blockdata["main_node"]:
			last_score =  changed_nptls
			best_score, changed_nptls, keepgoing, best_assignmen, iter_assignment = \
			     AI_MGSKmeans(iter_assignment, last_iter_assignment, best_assignment, \
			       keepgoing, best_score, stopercnt, minimum_group_size, log_main)
			if abs(last_score - changed_nptls)<1.0:
				if  times_around_fixed_value == 0:
					fixed_value = changed_nptls
					times_around_fixed_value +=1
				else:
					if abs(changed_nptls - fixed_value)<1.0: 
						times_around_fixed_value +=1
					else:
						times_around_fixed_value = 0
						fixed_value = changed_nptls
		else:
			iter_assignment = 0
			best_assignment = 0
			keepgoing       = 1
		iter_assignment          = wrap_mpi_bcast(iter_assignment,               Blockdata["main_node"], MPI_COMM_WORLD)
		best_assignment          = wrap_mpi_bcast(best_assignment,               Blockdata["main_node"], MPI_COMM_WORLD)
		times_around_fixed_value = bcast_number_to_all(times_around_fixed_value, Blockdata["main_node"], MPI_COMM_WORLD)
		keepgoing                = bcast_number_to_all(keepgoing,                Blockdata["main_node"], MPI_COMM_WORLD)
		total_iter +=1
		if times_around_fixed_value>=3: keepgoing = 0
		if keepgoing == 0: break
	# Finalize
	update_data_assignment(cdata, srdata, iter_assignment, proc_list, Tracker["nosmearing"], Blockdata["myid"])
	partition = get_sorting_parti_from_stack(cdata)
	del iter_assignment
	del last_iter_assignment
	del best_assignment
	#mpi_win_free(win_vol)
	#emnumpy1.unregister_numpy_from_emdata()
	#del emnumpy1
	if mask3D: del mask3D
	if(Blockdata["myid"] == Blockdata["main_node"]):
		lpartids = read_text_file(partids, -1)
		if len(lpartids) ==1: lpartids =lpartids[0]
		else: lpartids =lpartids[1]
		if best_score > Tracker["constants"]["stop_mgskmeans_percentage"]: premature  = 1
		write_text_file([partition.tolist(), lpartids], os.path.join(Tracker["directory"],"list.txt"))
		shutil.rmtree(os.path.join(Tracker["directory"], "tempdir"))
		fplist  = np.array([partition, np.array(lpartids)], dtype = np.int32)
	else: fplist = 0
	fplist     = wrap_mpi_bcast(fplist,         Blockdata["main_node"], MPI_COMM_WORLD)
	premature  = bcast_number_to_all(premature, Blockdata["main_node"], MPI_COMM_WORLD)
	if(Blockdata["myid"] == Blockdata["last_node"]):
		if clean_volumes:# always true unless in debug mode
			for jter in range(total_iter):
				for igroup in range(Tracker["number_of_groups"]): 
					os.remove(os.path.join(Tracker["directory"], "vol_grp%03d_iter%03d.hdf"%(igroup,jter)))	
	return (fplist).transpose().tolist(), premature, total_iter
	
def do_assignment_by_dmatrix_orien_group_minimum_group_size(dmatrix, \
      orien_group_members, number_of_groups, minimum_group_size_ratio):
	import numpy as np
	import random
	results = [[] for i in range(number_of_groups)]
	nima               = len(orien_group_members)
	minimum_group_size = int(minimum_group_size_ratio*nima/number_of_groups)
	auxmatrix          = np.full(nima, -1.0, dtype=np.float32)
	submatrix          = np.full((number_of_groups, nima), 0.0, dtype=np.float32)
	for im in range(number_of_groups):
		submatrix[im] = np.multiply(auxmatrix, dmatrix[im][orien_group_members])# sort in descending order
	tmp_array = np.argsort(submatrix, axis = 1)
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
			shuffle(duplicate_list)
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
	for iptl in range(len(kmeans_ptl_list)):
		max_indexes = np.argwhere(submatrix[:, kmeans_ptl_list[iptl]]<=submatrix[:, \
		    kmeans_ptl_list[iptl]][submatrix[:, kmeans_ptl_list[iptl]].argmin()])
		if len(max_indexes) >1:
			t = list(range(len(max_indexes)))
			shuffle(t)
			results[max_indexes[t[0]][0]].append(kmeans_ptl_list[iptl])
		else: results[max_indexes[0][0]].append(kmeans_ptl_list[iptl])
	iter_assignment = np.full(nima, -1, dtype=np.int32)
	for im in range(number_of_groups): iter_assignment[results[im]] = im
	del results
	del submatrix
	return iter_assignment
	
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
	from EMAN2          import Region, Vec2f
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
	
	if(Tracker["total_stack"] < Blockdata["nproc"]):
		ERROR("Wrong MPI settings!", "get_shrink_data_sorting", 1, Blockdata["myid"])
	else: image_start, image_end = MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"])
	lpartids  = lpartids[image_start:image_end]
	groupids  = groupids[image_start:image_end]
	nima      =	image_end - image_start
	data      = [None]*nima
	for im in range(nima):
		data[im] = get_im(Tracker["constants"]["orgstack"], lpartids[im])	
		try: phi, theta, psi, sx, sy, chunk_id, particle_group_id = \
		      partstack[lpartids[im]][0], partstack[lpartids[im]][1],\
		     partstack[lpartids[im]][2], partstack[lpartids[im]][3], partstack[lpartids[im]][4], \
		     partstack[lpartids[im]][5], partstack[lpartids[im]][6]
		except: phi, theta, psi, sx, sy, chunk_id, particle_group_id = \
		         partstack[lpartids[im]][0], partstack[lpartids[im]][1],\
		       partstack[lpartids[im]][2], partstack[lpartids[im]][3], partstack[lpartids[im]][4], \
		       partstack[lpartids[im]][5], -1		 
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
		t = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
		t.set_trans(Vec2f(0.0, 0.0))
		data[im].set_attr_dict({"chunk_id":chunk_id, "group": groupids[im], "particle_group": particle_group_id,\
		    "xform.projection": t, "npad": npad})
		if Tracker["applybckgnoise"]:
			data[im].set_attr("bckgnoise", Blockdata["bckgnoise"][particle_group_id])
			data[im].set_attr("qt", float(Tracker["constants"]["nnxo"]*Tracker["constants"]["nnxo"]))
		else: data[im].set_attr("bckgnoise", Blockdata["bckgnoise"]) # constant list
	return data
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
	
	mask2D  	= model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	nima 		= len(original_data)
	shrinkage 	= nxinit/float(Tracker["constants"]["nnxo"])
	radius 	    = int(Tracker["constants"]["radius"]*shrinkage + 0.5)
	txm    	    = float(nxinit-(nxinit//2+1) - radius)
	txl    	    = float(radius - nxinit//2+1)

	if Blockdata["bckgnoise"] :
		oneover = []
		nnx = Blockdata["bckgnoise"][0].get_xsize()
		for i in range(len(Blockdata["bckgnoise"])):
			temp = [0.0]*nnx
			for k in range(nnx):
				if( Blockdata["bckgnoise"][i].get_value_at(k) > 0.0):
					temp[k] = 1.0/sqrt(Blockdata["bckgnoise"][i].get_value_at(k))
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
		set_params_proj(data[im],[phi,theta,psi,max(min(sx*shrinkage,txm),txl),\
		    max(min(sy*shrinkage,txm),txl)])
		if not return_real: data[im].set_attr("padffted",1)
		data[im].set_attr("npad",npad)
		if Blockdata["bckgnoise"]:
			temp = Blockdata["bckgnoise"][data[im].get_attr("particle_group")]
			###  Do not adjust the values, we try to keep everything in the same Fourier values.
			data[im].set_attr("bckgnoise", [temp[i] for i in range(temp.get_xsize())])
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
	from EMAN2          import Region, Vec2f
	from utilities      import wrap_mpi_bcast, read_text_row, get_im, set_params_proj
	import numpy as np
	# functions:
	# read in data
	if( Blockdata["myid"] == Blockdata["main_node"]):
		lpartids = read_text_file(partids, -1)
		if len(lpartids) == 1: 
			lpartids = lpartids[0]
			groupids = [-1 for im in range(len(lpartids))]
		else:
			groupids = lpartids[0]                  
			lpartids = lpartids[1]
	else: 
		lpartids   = 0
		groupids   = 0
	lpartids = wrap_mpi_bcast(lpartids, Blockdata["main_node"])
	groupids = wrap_mpi_bcast(groupids, Blockdata["main_node"])
	
	Tracker["total_stack"] = len(lpartids)
	if(Blockdata["myid"] == Blockdata["main_node"]): partstack = read_text_row(partstack)
	else:  partstack = 0
	partstack = wrap_mpi_bcast(partstack, Blockdata["main_node"])
	if(Blockdata["myid"] == Blockdata["main_node"]): previous_partstack = read_text_row(previous_partstack)
	else:  previous_partstack = 0
	previous_partstack = wrap_mpi_bcast(previous_partstack, Blockdata["main_node"])
	if(Tracker["total_stack"] < Blockdata["nproc"]):
		ERROR("Number of processors in use is larger than the total number of images", \
		  "get_data_and_prep", 1, Blockdata["myid"])
	else: image_start, image_end = MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"])
	lpartids          = lpartids[image_start:image_end]
	nima              = image_end - image_start
	data              = [ None for im in range(nima)]
	norm_per_particle = [ None for im in range(nima)]
	groupids = groupids[image_start: image_end]
	for im in range(nima):
		image = get_im(Tracker["constants"]["orgstack"], lpartids[im])
		try: phi, theta, psi, sx, sy, chunk_id, particle_group_id, mnorm = partstack[lpartids[im]][0], \
		   partstack[lpartids[im]][1], partstack[lpartids[im]][2], \
			partstack[lpartids[im]][3], partstack[lpartids[im]][4], partstack[lpartids[im]][5], \
			   partstack[lpartids[im]][6], partstack[lpartids[im]][7]
		except:
			phi, theta, psi, sx, sy, chunk_id, particle_group_id, mnorm  = partstack[lpartids[im]][0], \
			    partstack[lpartids[im]][1], partstack[lpartids[im]][2], \
		  partstack[lpartids[im]][3], partstack[lpartids[im]][4], partstack[lpartids[im]][5], -1, 1.
		sx1, sy1 = previous_partstack[lpartids[im]][3], previous_partstack[lpartids[im]][4]
		t = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
		t.set_trans(Vec2f(0.0, 0.0))
		image.set_attr_dict( {"chunk_id":chunk_id, "group":groupids[im], "particle_group":\
		    particle_group_id, "previous_shifts": [sx1, sy1], "current_shifts": [sx, sy], "xform.projection": t})
		norm_per_particle[im] = mnorm
		data[im]              = image
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
	if( Blockdata["myid"] == Blockdata["main_node"]):
		tmp_list = read_text_row(paramstructure_dict_file)
	else: tmp_list = 0
	tmp_list = wrap_mpi_bcast(tmp_list, Blockdata["main_node"])
	pdict    = {}
	for im in range(len(lcore)):pdict[im] = tmp_list[lcore[im]]
	oldparamstructure             =  []
	nptl                          =  0
	last_old_paramstructure_file  =  None
	mpi_barrier(MPI_COMM_WORLD)
	for iproc in range(Blockdata["nproc"]):
		if (Blockdata["myid"] == iproc): #always read oldparamstructure sequentially
			while nptl < nima:
				[jason_of_cpu_id, chunk_id, iteration, ptl_id_on_cpu, global_index] = pdict[nptl]
				old_paramstructure_file = os.path.join(paramstructure_dir, "oldparamstructure_%d_%03d_%03d.json"%\
				  (chunk_id, jason_of_cpu_id, iteration))
				if old_paramstructure_file != last_old_paramstructure_file:
					fout = open(old_paramstructure_file,'r')
					paramstructure = convert_json_fromunicode(json.load(fout))
					fout.close()
				last_old_paramstructure_file = old_paramstructure_file
				oldparamstructure.append(paramstructure[ptl_id_on_cpu])	
				nptl +=1
			mpi_barrier(MPI_COMM_WORLD)
	return oldparamstructure
	
###7 copy oldparamstructures from meridien regardless the previous number of cpus 
def copy_oldparamstructure_from_meridien_MPI(selected_iteration, log_main):
	global Tracker, Blockdata
	from utilities    import read_text_row, write_text_row, read_text_file,wrap_mpi_bcast
	from applications import MPI_start_end
	import json
	Tracker["directory"]          = os.path.join(Tracker["constants"]["masterdir"], "main%03d"%selected_iteration)
	Tracker["paramstructure_dir"] = os.path.join(Tracker["directory"], "oldparamstructure")
	old_refinement_iter_directory = os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iteration)
	old_refinement_previous_iter_directory = os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%(selected_iteration-1))
	if( Blockdata["myid"] == Blockdata["main_node"]):
		if not os.path.exists(Tracker["paramstructure_dir"]):
			os.mkdir(os.path.join(Tracker["constants"]["masterdir"], "main%03d"%selected_iteration))
			os.mkdir(Tracker["paramstructure_dir"])
	Tracker["refang"]     = read_text_row(os.path.join(old_refinement_iter_directory, "refang.txt"))
	if( Blockdata["myid"] == Blockdata["main_node"]): write_text_row(Tracker["refang"], os.path.join(Tracker["directory"], "refang.txt"))
	Tracker["rshifts"] = read_text_row(os.path.join(old_refinement_iter_directory, "rshifts.txt"))
	if( Blockdata["myid"] == Blockdata["main_node"]): write_text_row(Tracker["refang"], os.path.join(Tracker["directory"], "rshifts.txt"))
	my_last_params = read_text_file(os.path.join(old_refinement_previous_iter_directory, "params_%03d.txt"%(selected_iteration-1)), -1)
	my_parstack    = read_text_file(os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt"), -1)
	if( Blockdata["myid"] == Blockdata["main_node"]):
		my_parstack[3:5]= my_last_params[3:5]
		write_text_file(my_parstack, os.path.join(Tracker["constants"]["masterdir"], "previous_refinement_parameters.txt"))
	Tracker["previous_parstack"]   = os.path.join(Tracker["constants"]["masterdir"], "previous_refinement_parameters.txt")
	nproc_previous = 0
	procid         = 0
	old_refinement_iter_dir = os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iteration)
	if Blockdata["myid"] == Blockdata["main_node"]:
		while os.path.exists(os.path.join(old_refinement_iter_dir,\
		   "oldparamstructure","oldparamstructure_%01d_%03d_%03d.json"%(procid, \
		      nproc_previous, selected_iteration))):
			nproc_previous += 1
	nproc_previous = bcast_number_to_all(nproc_previous, Blockdata["main_node"], MPI_COMM_WORLD)
	Blockdata["nproc_previous"] = nproc_previous
	oldparamstructure           =[[], []]
	local_dict                  = {}
	for procid in range(2):
		smearing_list = []
		if( Blockdata["myid"] == Blockdata["main_node"]): lcore = read_text_file(\
		   os.path.join(Tracker["constants"]["masterdir"], "chunk_%d.txt"%procid))
		else: lcore = 0
		lcore = wrap_mpi_bcast(lcore, Blockdata["main_node"], MPI_COMM_WORLD)	
		psize = len(lcore)
		oldparamstructure[procid] = []
		im_start, im_end   = MPI_start_end(psize, Blockdata["nproc"], Blockdata["myid"])
		local_lcore = lcore[im_start:im_end]
		istart_old_proc_id = -1
		iend_old_proc_id   = -1
		plist              = []	
		for iproc_old in range(nproc_previous):
			im_start_old, im_end_old = MPI_start_end(psize, nproc_previous, iproc_old)
			if (im_start>= im_start_old) and im_start <=im_end_old: istart_old_proc_id = iproc_old
			if (im_end  >= im_start_old) and im_end <=im_end_old:   iend_old_proc_id   = iproc_old
			plist.append([im_start_old, im_end_old])
		ptl_on_this_cpu = im_start
		nptl_total      = 0
		for iproc_index_old in range(istart_old_proc_id, iend_old_proc_id+1):
			fout = open(os.path.join(Tracker["constants"]["refinement_dir"],\
			   "main%03d"%selected_iteration, "oldparamstructure", \
			      "oldparamstructure_%01d_%03d_%03d.json"%(procid, \
			 iproc_index_old, selected_iteration)),'r')
			oldparamstructure_on_old_cpu = convert_json_fromunicode(json.load(fout))
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
		mpi_barrier(MPI_COMM_WORLD)
	
		for icpu in range(Blockdata["nproc"]):#dump to disk one by one
			if  Blockdata["myid"] == icpu:
				fout = open(os.path.join(Tracker["constants"]["masterdir"], "main%03d"%selected_iteration, \
				   "oldparamstructure", "oldparamstructure_%01d_%03d_%03d.json"%(procid, \
					Blockdata["myid"], selected_iteration)),'w')
				json.dump(oldparamstructure[procid], fout)
				fout.close()
				mpi_barrier(MPI_COMM_WORLD)
		mpi_barrier(MPI_COMM_WORLD)
	# output number of smearing
	smearing_dict = {}
	tchunk        = []
	for procid in range(2):
		if Blockdata["myid"] == Blockdata["main_node"]:
			chunk = read_text_file(os.path.join(Tracker["constants"]["masterdir"], "chunk_%d.txt"%procid))
			chunk_size = len(chunk)
			smearing_list =[ None for i in range(chunk_size) ]
		else: chunk_size  = 0
		chunk_size = bcast_number_to_all(chunk_size, Blockdata["main_node"], MPI_COMM_WORLD)
		local_smearing_list = []
		for im in range(len(oldparamstructure[procid])):local_smearing_list.append(len(oldparamstructure[procid][im][2]))
			
		if Blockdata["myid"] == Blockdata["main_node"]:
			im_start_old, im_end_old = MPI_start_end(chunk_size, Blockdata["nproc"], Blockdata["main_node"])
			for im in range(len(local_smearing_list)): smearing_list[im_start_old+im] = local_smearing_list[im]
		mpi_barrier(MPI_COMM_WORLD)
		
		if  Blockdata["myid"] != Blockdata["main_node"]:
			wrap_mpi_send(local_smearing_list, Blockdata["main_node"], MPI_COMM_WORLD)
		else:
			for iproc in range(Blockdata["nproc"]):
				if iproc != Blockdata["main_node"]:
					im_start_old, im_end_old = MPI_start_end(chunk_size, Blockdata["nproc"], iproc)
					dummy = wrap_mpi_recv(iproc, MPI_COMM_WORLD)
					for idum in range(len(dummy)): smearing_list[idum + im_start_old] = dummy[idum]
				else: pass
			write_text_file(smearing_list, os.path.join(Tracker["constants"]["masterdir"], "smearing_%d.txt"%procid))
			for im in range(len(chunk)): smearing_dict[chunk[im]] =  smearing_list[im]
			tchunk +=chunk
		mpi_barrier(MPI_COMM_WORLD)
		
	if Blockdata["myid"] == Blockdata["main_node"]:
		tchunk.sort()
		all_smearing = [None]*len(tchunk)
		for im in range(len(tchunk)): all_smearing[im] = smearing_dict[tchunk[im]]
		write_text_file(all_smearing, os.path.join(Tracker["constants"]["masterdir"], "all_smearing.txt"))
		full_dict_list = [ None for im in range(Tracker["constants"]["total_stack"])]
		for key, value in list(local_dict.items()):full_dict_list[key] = value
		Tracker["constants"]["smearing_file"] = os.path.join(Tracker["constants"]["masterdir"], "all_smearing.txt")
	mpi_barrier(MPI_COMM_WORLD)
	
	for icpu in range(Blockdata["nproc"]):
		if Blockdata["myid"] == icpu and Blockdata["myid"] != Blockdata["main_node"]:
			wrap_mpi_send(local_dict, Blockdata["main_node"], MPI_COMM_WORLD)
		elif Blockdata["myid"] != icpu and Blockdata["myid"] == Blockdata["main_node"]:
			local_dict = wrap_mpi_recv(icpu, MPI_COMM_WORLD)
			for key, value in list(local_dict.items()):
				full_dict_list[key] = value
		else: pass
		mpi_barrier(MPI_COMM_WORLD)
	Tracker["paramstructure_dict"] = \
	   os.path.join(Tracker["constants"]["masterdir"], "paramstructure_dict.txt")
	if Blockdata["myid"] == Blockdata["main_node"]: 
		write_text_row(full_dict_list, Tracker["paramstructure_dict"])
	mpi_barrier(MPI_COMM_WORLD)
	return

def get_smearing_info(nproc_previous, selected_iteration, total_stack, my_dir, refinement_dir):
	global Tracker, Blockdata
	import numpy as np
	# Do not alter any values in Tracker and Blockdata
	oldparamstructure    =[[], []]
	local_dict           = {}
	for procid in range(2):
		smearing_list = []
		if( Blockdata["myid"] == Blockdata["main_node"]): lcore = read_text_file(\
		   os.path.join(my_dir, "chunk_%d.txt"%procid))
		else: lcore = 0
		lcore = wrap_mpi_bcast(lcore, Blockdata["main_node"], MPI_COMM_WORLD)	
		psize = len(lcore)
		oldparamstructure[procid] = []
		im_start, im_end   = MPI_start_end(psize, Blockdata["nproc"], Blockdata["myid"])
		local_lcore        = lcore[im_start:im_end]
		istart_old_proc_id = -1
		iend_old_proc_id   = -1
		plist              = []	
		for iproc_old in range(nproc_previous):
			im_start_old, im_end_old = MPI_start_end(psize, nproc_previous, iproc_old)
			if (im_start>= im_start_old) and im_start <=im_end_old: istart_old_proc_id = iproc_old
			if (im_end  >= im_start_old) and im_end <=im_end_old:   iend_old_proc_id   = iproc_old
			plist.append([im_start_old, im_end_old])
		ptl_on_this_cpu = im_start
		nptl_total      = 0
		for iproc_index_old in range(istart_old_proc_id, iend_old_proc_id+1):
			fout = open(os.path.join(refinement_dir,\
			   "main%03d"%selected_iteration, "oldparamstructure", \
			      "oldparamstructure_%01d_%03d_%03d.json"%(procid, \
			 iproc_index_old, selected_iteration)),'r')
			oldparamstructure_on_old_cpu = convert_json_fromunicode(json.load(fout))
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
		mpi_barrier(MPI_COMM_WORLD)
		
	# output number of smearing
	smearing_dict = {}
	tchunk        = []
	for procid in range(2):
		if Blockdata["myid"] == Blockdata["main_node"]:
			chunk = read_text_file(os.path.join(my_dir, "chunk_%d.txt"%procid))
			chunk_size = len(chunk)
			smearing_list =[ None for i in range(chunk_size) ]
		else: chunk_size  = 0
		chunk_size = bcast_number_to_all(chunk_size, Blockdata["main_node"], MPI_COMM_WORLD)
		local_smearing_list = []
		for im in range(len(oldparamstructure[procid])):
			local_smearing_list.append(len(oldparamstructure[procid][im][2]))
		if Blockdata["myid"] == Blockdata["main_node"]:
			im_start_old, im_end_old = MPI_start_end(chunk_size, Blockdata["nproc"], Blockdata["main_node"])
			for im in range(len(local_smearing_list)): smearing_list[im_start_old+im] = local_smearing_list[im]
		mpi_barrier(MPI_COMM_WORLD)
		if  Blockdata["myid"] != Blockdata["main_node"]:
			wrap_mpi_send(local_smearing_list, Blockdata["main_node"], MPI_COMM_WORLD)
		else:
			for iproc in range(Blockdata["nproc"]):
				if iproc != Blockdata["main_node"]:
					im_start_old, im_end_old = MPI_start_end(chunk_size, Blockdata["nproc"], iproc)
					dummy = wrap_mpi_recv(iproc, MPI_COMM_WORLD)
					for idum in range(len(dummy)): smearing_list[idum + im_start_old] = dummy[idum]
				else: pass
			write_text_file(smearing_list, os.path.join(my_dir, "smearing_%d.txt"%procid))
			for im in range(len(chunk)): smearing_dict[chunk[im]] =  smearing_list[im]
			tchunk +=chunk
		mpi_barrier(MPI_COMM_WORLD)
	
	if Blockdata["myid"] == Blockdata["main_node"]:
		tchunk.sort()
		all_smearing = np.full(len(tchunk), 0.0, dtype=np.float32)
		for im in range(len(tchunk)): 
			all_smearing[im] = smearing_dict[tchunk[im]]
	else: all_smearing = 0
	all_smearing = wrap_mpi_bcast(all_smearing, Blockdata["main_node"], MPI_COMM_WORLD)
	return all_smearing
### 8
def precalculate_shifted_data_for_recons3D(prjlist, paramstructure, refang, rshifts, delta, avgnorms, \
       nxinit, nnxo, nosmearing, norm_per_particle = None, upweighted=False, nsmear =-1):
	from utilities    import random_string, get_im, findall, info, model_blank
	from filter	      import filt_table
	from fundamentals import fshift
	import types, datetime, copy
	nny = prjlist[0].get_ysize()
	if nosmearing:
		for im in range(len(prjlist)):
			prjlist[im].set_attr_dict({"padffted":1, "is_fftpad":1,"is_fftodd":0, "is_complex_ri":1, "is_complex":1})
			if not upweighted:prjlist[im] = filt_table(prjlist[im], prjlist[im].get_attr("bckgnoise"))
			prjlist[im].set_attr("wprob", 1.0)
		return prjlist
	else:
		if norm_per_particle == None: norm_per_particle = len(prjlist)*[1.0]
		recdata_list = [[] for im in range(len(prjlist))]
		rshifts_shrank = copy.deepcopy(rshifts)
		for im in range(len(rshifts_shrank)):
			rshifts_shrank[im][0] *= float(nxinit)/float(nnxo)
			rshifts_shrank[im][1] *= float(nxinit)/float(nnxo)
		if not Tracker["constants"]["compute_on_the_fly"]:
			num_on_the_fly = len(prjlist)
		else:
			num_on_the_fly =  Tracker["num_on_the_fly"][Blockdata["myid"]] 
		for im in range(num_on_the_fly):
			bckgn = prjlist[im].get_attr("bckgnoise")
			ct = prjlist[im].get_attr("ctf")
			avgnorm     = Tracker["avgnorm"][prjlist[im].get_attr("chunk_id")]
			numbor      = len(paramstructure[im][2])
			ipsiandiang = [ paramstructure[im][2][i][0]/1000  for i in range(numbor) ]
			allshifts   = [ paramstructure[im][2][i][0]%1000  for i in range(numbor) ]
			probs       = [ paramstructure[im][2][i][1] for i in range(numbor) ]
			tdir = list(set(ipsiandiang))
			data = [None]*len(rshifts_shrank)
			for ii in range(len(tdir)):
				#  Find the number of times given projection direction appears on the list, it is the number of different shifts associated with it.
				lshifts = findall(tdir[ii], ipsiandiang)
				toprab  = 0.0
				for ki in range(len(lshifts)):toprab += probs[lshifts[ki]]
				recdata = EMData(nny,nny, 1, False)
				recdata.set_attr("is_complex",0)
				for ki in range(len(lshifts)):
					lpt = allshifts[lshifts[ki]]
					if(data[lpt] == None):
						data[lpt] = fshift(prjlist[im], rshifts_shrank[lpt][0], rshifts_shrank[lpt][1])
						data[lpt].set_attr("is_complex",0)
					Util.add_img(recdata, Util.mult_scalar(data[lpt], probs[lshifts[ki]]/toprab))
				recdata.set_attr_dict({"padffted":1, "is_fftpad":1,"is_fftodd":0, "is_complex_ri":1, "is_complex":1}) # preset already
				if not upweighted:recdata = filt_table(recdata, bckgn)
				ipsi = tdir[ii]%100000
				iang = tdir[ii]/100000
				recdata.set_attr_dict({"wprob": toprab*avgnorm/norm_per_particle[im],\
				   "xform.projection":Transform({"type":"spider", "phi":refang[iang][0],\
					 "theta":refang[iang][1], "psi":refang[iang][2]+ipsi*delta}), "bckgnoise":bckgn, "ctf":ct})
				recdata_list[im].append(recdata)
			del recdata, tdir, ipsiandiang, allshifts, probs, data
		if num_on_the_fly<len(prjlist):
			for im in range(num_on_the_fly, len(prjlist)):
				recdata_list[im].append(prjlist[im])
		return recdata_list
##### read data/paramstructure ends

###=====<----downsize data---->>>>
def downsize_data_for_sorting(original_data, return_real = False, preshift = True, npad = 1, norms = None):
	# The function will read from stack a subset of images specified in partids
	#   and assign to them parameters from partstack with optional CTF application and shifting of the data.
	# So, the lengths of partids and partstack are the same.
	global Tracker, Blockdata
	from fundamentals	import resample, fshift, cyclic_shift
	from filter			import filt_ctf
	from applications	import MPI_start_end
	from EMAN2          import Region, Vec2f
	from utilities import same_ctf
	# functions:
	# read in data
	# apply mask, and prepare focus projection if focus3D is specified
	# return  1. cdata: data for image comparison, always in Fourier format
	#         2. rdata: data for reconstruction, 4nn return real image
	mask2D		= model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	shrinkage 	= Tracker["nxinit"]/float(Tracker["constants"]["nnxo"])
	radius      = int(Tracker["constants"]["radius"] * shrinkage +0.5)
	compute_noise(Tracker["nxinit"])
	if Tracker["applybckgnoise"]:
		oneover = []
		nnx     = len(Blockdata["bckgnoise"][0])
		for i in range(len(Blockdata["bckgnoise"])):
			temp = [0.0]*nnx
			for k in range(nnx):
				if(Blockdata["bckgnoise"][i][k] > 0.0): temp[k] = 1.0/sqrt(Blockdata["bckgnoise"][i][k])
			oneover.append(temp)
		del temp
	if Tracker["constants"]["focus3D"]: # focus mask is applied
		if Blockdata["myid"] == Blockdata["main_node"]:
			focus3d    = get_im(Tracker["constants"]["focus3D"])
			focus3d_nx = focus3d.get_xsize()
			if focus3d_nx != Tracker["nxinit"]: # So the decimated focus volume can be directly used
				focus3d = resample(focus3d, float(Tracker["nxinit"])/float(focus3d_nx))
		else: focus3d = model_blank(Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"])
		bcast_EMData_to_all(focus3d, Blockdata["myid"], Blockdata["main_node"])
		focus3d = prep_vol(focus3d, 1, 1)
	#  Preprocess the data
	nima   = len(original_data)
	cdata  = [None for im in range(nima)]
	rdata  = [None for im in range(nima)]
	fdata  = [None for im in range(nima)]
	ctfs   = [None for im in range(nima)]
	ctfa = EMAN2Ctf()
	for im in range(nima):
		image = original_data[im].copy()
		chunk_id          = image.get_attr("chunk_id")
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
		if not Tracker["nosmearing"] and norms: cimage *= Tracker["avgnorm"][chunk_id]/norms[im] #norm correction
		if Tracker["applybckgnoise"]:
			if Tracker["applymask"]:
				if Tracker["constants"]["hardmask"]: 
					cimage = cosinemask(cimage, radius = Tracker["constants"]["radius"])
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
			ctf_params = rimage.get_attr('ctf')
			rimage      = fdecimate(rimage, Tracker["nxinit"]*npad, Tracker["nxinit"]*npad, 1, False, False)
			cimage      = fdecimate(cimage, Tracker["nxinit"]*npad, Tracker["nxinit"]*npad, 1, False, False)
			ctf_params.apix = ctf_params.apix/shrinkage
			rimage.set_attr('ctf', ctf_params)
			rimage.set_attr('ctf_applied', 0)
			if return_real :  rimage = fft(rimage)
		else:
			ctf_params = rimage.get_attr_default("ctf", False)
			if ctf_params:
				ctf_params.apix = ctf_params.apix/shrinkage
				rimage.set_attr('ctf', ctf_params)
				rimage.set_attr('ctf_applied', 0)
			rimage  = fdecimate(rimage, nxinit*npad, nxinit*npad, 1, True, False)
			cimage  = fdecimate(cimage, nxinit*npad, nxinit*npad, 1, True, False)
			apix   = Tracker["constants"]["pixel_size"]
			rimage.set_attr('apix', apix/shrinkage)
		cimage.set_attr("padffted",1)
		cimage.set_attr("npad", npad)
		if not return_real:	
			rimage.set_attr("padffted",1)
			rimage.set_attr("npad", npad)	
		t = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
		t.set_trans(Vec2f(0.0, 0.0))
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
			focusmask = binarize(prgl(focus3d, [phi, theta, psi, 0.0, 0.0], 1, True), 1)
			cdata[im] = fft(focusmask*fft(cdata[im]))
			fdata[im] = focusmask
		cdata[im].set_attr("is_complex",0)
		cdata[im].set_attr("particle_group", particle_group_id)
		if Tracker["constants"]["CTF"]:
			if im ==0: ny = cdata[im].get_ysize()
			if not same_ctf(rdata[im].get_attr("ctf"), ctfa):
				ctfs[im]   = ctf_img_real(ny, rdata[im].get_attr("ctf"))
			else: ctfs[im] = ctfs[im-1].copy()
	return cdata, rdata, fdata, ctfs
##=====<----for 3D----->>>>
def downsize_data_for_rec3D(original_data, particle_size, return_real = False, npad = 1):
	# The function will read from stack a subset of images specified in partids
	#   and assign to them parameters from partstack with optional CTF application and shifting of the data.
	# So, the lengths of partids and partstack are the same.	
	global Tracker, Blockdata
	from fundamentals	import resample, fshift
	from filter			import filt_ctf
	from applications	import MPI_start_end
	from EMAN2          import Region, Vec2f
	# functions:
	# read in data
	# apply mask, and prepare focus projection if focus3D is specified
	# return  1. cdata: data for image comparison, always in Fourier format
	#         2. rdata: data for reconstruction, 4nn return real image
	nima       = len(original_data)
	rdata      = [None for im in range(nima)]
	mask2D     = model_circle(Tracker["constants"]["radius"],\
	   Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	shrinkage  = particle_size/float(Tracker["constants"]["nnxo"])
	radius     = int(Tracker["constants"]["radius"] * shrinkage +0.5)
	for im in range(nima):
		image = original_data[im].copy()
		chunk_id = image.get_attr("chunk_id")
		try: particle_group_id = image.get_attr("particle_group")
		except: particle_group_id = -1
		phi,theta,psi,s2x,s2y = get_params_proj(image, xform = "xform.projection")
		t = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
		t.set_trans(Vec2f(0.0, 0.0))
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
		image.set_attr_dict({"chunk_id": chunk_id, "particle_group": particle_group_id, \
		   "xform.projection": t})
		rdata[im] =  image
		if Tracker["applybckgnoise"]:
			rdata[im].set_attr("bckgnoise", Blockdata["bckgnoise"][rdata[im].get_attr("particle_group")])
		else: rdata[im].set_attr("bckgnoise", Blockdata["bckgnoise"])
	return rdata
### end of downsize

###=====<--- comparison
def compare_two_images_eucd(data, ref_vol, fdata, ctfimgs):
	global Tracker, Blockdata
	import numpy as np
	peaks = np.zeros(shape = len(data), dtype=float)
	ny      = data[0].get_ysize()
	ref_vol = prep_vol(ref_vol, npad = 2, interpolation_method = 1)
	#ctfs    = EMAN2Ctf()
	qt = float(Tracker["constants"]["nnxo"]*Tracker["constants"]["nnxo"])
	m = Util.unrollmask(ny,ny)
	for im in range(len(data)):
		#current_ctf = data[im].get_attr('ctf')
		#if( not same_ctf(current_ctf,ctfs) ):
		#	ctfs = ctf_img_real(ny, current_ctf)
		phi, theta, psi, s2x, s2y = get_params_proj(data[im], xform = "xform.projection")
		rtemp = prgl(ref_vol,[phi, theta, psi, 0.0,0.0], 1, False)
		if Tracker["constants"]["focus3D"]: rtemp = fft(fft(rtemp)*fdata[im])
		if Tracker["applybckgnoise"]:
			peaks[im] = -Util.sqed(data[im], rtemp, ctfimgs[im], \
			    Blockdata["unrolldata"][data[im].get_attr("particle_group")])/qt
		else:
			peaks[im] = -Util.sqed(data[im], rtemp, ctfimgs[im], m)/qt
	return peaks

def compare_two_images_cross(data, ref_vol, ctfimgs):
	global Tracker, Blockdata
	import numpy as np
	from utilities import same_ctf
	ny    = data[0].get_ysize()
	m     = Util.unrollmask(ny,ny)
	peaks = np.zeros(shape = len(data), dtype=float)
	volft = prep_vol(ref_vol, 2, 1)
	#  Ref is in reciprocal space
	#ctfa = EMAN2Ctf()		
	for im in range(len(data)):
		phi, theta, psi, s2x, s2y = get_params_proj(data[im], xform = "xform.projection")
		ref = prgl( volft, [phi, theta, psi, 0.0, 0.0], 1, False)
		#current_ctf = data[im].get_attr('ctf')
		#if not same_ctf(current_ctf, ctfa):
		#	ctfa = current_ctf
		#	ctfimg  = ctf_img_real(ny, ctfa)
		Util.mulclreal(ref, ctfimgs[im])
		nrmref = sqrt(Util.innerproduct(ref, ref, m))
		if Tracker["applybckgnoise"]:
			peaks[im] = Util.innerproduct(ref, data[im], \
			   Blockdata["unrolldata"][data[im].get_attr("particle_group")])/nrmref
		else:         
			peaks[im] = Util.innerproduct(ref, data[im], m)/nrmref
	return peaks		
#####==========--------------------utilities of creating random assignments
def create_nrandom_lists(partids, number_of_groups, number_of_runs):
	# the second column denotes orignal particle IDs
	# the first column is randomized group ID 
	global Tracker, Blockdata
	import copy
	import random
	from   utilities import wrap_mpi_bcast, read_text_file, write_text_file
	if Blockdata["myid"] == Blockdata["main_node"]:
		random_assignment = []
		data_list = read_text_file(partids, -1)
		if len(data_list)==1: sorting_data_list = data_list[0]
		else:                 sorting_data_list = data_list[1]
		random.seed()
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
	random_assignment = wrap_mpi_bcast(random_assignment , Blockdata["main_node"])
	return random_assignment

def create_nrandom_lists_from_given_pids(work_dir, partids, number_of_groups, number_of_runs):
	# the second column denotes orignal particle IDs
	# the first column is randomized group ID 
	global Tracker, Blockdata
	import copy
	import random
	from   utilities import wrap_mpi_bcast, read_text_file, write_text_file
	if Blockdata["myid"] == Blockdata["main_node"]:
		random_assignment = []
		data_list = read_text_file(partids, -1)
		if len(data_list)==1: sorting_data_list= data_list[0]
		else: sorting_data_list = data_list[1]
		random.seed()
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
			write_text_row(assignment, os.path.join(work_dir,"independent_index_%03d.txt"%index_of_random))
			random_assignment.append(assignment)
			del assignment
			del ll
	else: random_assignment = 0
	random_assignment = wrap_mpi_bcast(random_assignment, Blockdata["main_node"])
	return  random_assignment
	
def find_smallest_group(clusters):
	min_size =[len(clusters[0]), [0]]
	for ic in range(1, len(clusters)):
		if len(cluster[ic]) < min_size[0]: min_size = [len(clusters[ic]), [ic]]
		elif len(cluster[ic]) == min_size[0]: min_size[1].append(ic)
	if len(min_size[1])>=1: shuffle(min_size[1])
	return min_size[1][0]
	
def assign_unaccounted_elements_mpi(glist, clusters, img_per_grp):
	# assign unaccounted images by group probabilities
	global Tracker, Blockdata
	import random
	import copy
	icut = 3*img_per_grp//2
	if Blockdata["myid"]== Blockdata["main_node"]:
		for ic in range(len(clusters)):
			if len(clusters)>(3*img_per_grp//2):
				shuffle(clusters[ic])
				glist +=clusters[ic][icut:]
	else:
		glist    = 0
		clusters = 0
	clusters = wrap_mpi_bcast(clusters, Blockdata["main_node"], MPI_COMM_WORLD)
	glist    = wrap_mpi_bcast(glist,    Blockdata["main_node"], MPI_COMM_WORLD)
	
	if len(glist) <= Blockdata["nproc"]*30:
		if Blockdata["myid"]== Blockdata["main_node"]:
			clusters = assign_unaccounted_inverse_proportion_to_size(glist, clusters, img_per_grp)
		else: clusters = 0
		clusters = wrap_mpi_bcast(clusters, Blockdata["main_node"], MPI_COMM_WORLD)
	else:
		slist = []
		clist = []
		for ic in range(len(clusters)):
			if len(clusters[ic])<=img_per_grp:
				slist.append(max(1.- float(len(clusters[ic]))/float(img_per_grp), 0.05))
			else: slist.append(0.05)
		uplist = copy.deepcopy(glist)
		image_start, image_end = MPI_start_end(len(uplist), Blockdata["nproc"], Blockdata["myid"])
		uplist = uplist[image_start:image_end]
		nsize  = len(uplist)//3
		for ichunk in range(3):
			if ichunk !=2: ulist = uplist[ichunk*nsize:(ichunk+1)*nsize]
			else: ulist = uplist[ichunk*nsize:]
			while len(ulist)>0:
				im =  random.randint(0, len(clusters)-1)
				shuffle(ulist)
				if slist[im] > random.random():
					clusters[im].append(ulist[0])
					if len(clusters[ic])<=img_per_grp:
						slist[im] = max(1.- float(len(clusters[im]))/float(img_per_grp), 0.05)
					else: slist[im] = 0.05
					del ulist[0]
					if len(ulist)== 0: break
				else: continue
			mpi_barrier(MPI_COMM_WORLD)
			if Blockdata["myid"]!= Blockdata["main_node"]:
				 wrap_mpi_send(clusters, Blockdata["main_node"], MPI_COMM_WORLD)
			else:
				for iproc in range(Blockdata["nproc"]):
					if iproc != Blockdata["main_node"]:
						dummy = wrap_mpi_recv(iproc, MPI_COMM_WORLD)
						for ic in range(len(clusters)):clusters[ic]+=dummy[ic]
			mpi_barrier(MPI_COMM_WORLD)
			if Blockdata["myid"]== Blockdata["main_node"]:
				for ic in  range(len(clusters)): 
					clusters[ic] = list(set(clusters[ic]))
					if len(clusters[ic])<=img_per_grp:
						slist[ic] = max(1.- float(len(clusters[ic]))/float(img_per_grp), 0.05)
					else: slist[ic] = 0.05
			else:
				slist    = 0 
				clusters = 0
			slist    = wrap_mpi_bcast(slist,    Blockdata["main_node"], MPI_COMM_WORLD)
			clusters = wrap_mpi_bcast(clusters, Blockdata["main_node"], MPI_COMM_WORLD)
	return clusters
	
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
		clusters         = wrap_mpi_bcast(clusters,         Blockdata["main_node"], MPI_COMM_WORLD)
		unaccounted_list = wrap_mpi_bcast(unaccounted_list, Blockdata["main_node"], MPI_COMM_WORLD)
	
	avg_size = N//number_of_clusters
	m = number_of_clusters - len(clusters)
	if Blockdata["myid"] == Blockdata["main_node"]:
		large_clusters = []
		for ic in range(len(clusters)):
			if len(clusters[ic]) > 2*avg_size:
				large_clusters.append(clusters[ic])
	else: large_clusters = 0
	large_clusters = wrap_mpi_bcast(large_clusters, Blockdata["main_node"], MPI_COMM_WORLD)
	L = len(large_clusters)
	if m == 0 and L == 0:
		out_clusters = assign_unaccounted_elements_mpi(unaccounted_list, clusters, avg_size)
	else:
		if m !=0: empty_clusters =[[] for ie in range(m)]
		else: empty_clusters = []
		out_clusters = fill_no_large_groups_and_unaccounted_to_m_and_rcluster_mpi(\
				 unaccounted_list, empty_clusters, clusters, NUACC, NACC)
	sorted_out_clusters = []
	for im in range(len(out_clusters)):sorted_out_clusters.append(sorted(out_clusters[im]))
	return sorted_out_clusters

def select_fixed_size_cluster_from_alist(ulist, img_per_grp):
	cluster = []
	shuffle(ulist)
	cluster += ulist[0:img_per_grp]
	del ulist[0:img_per_grp]
	return cluster, ulist
		
def swap_clusters_small_NUACC(glist, clusters, swap_ratio):
	from random import shuffle
	slist     = [None for im in range(len(clusters))]
	temp_list = []
	ulist = copy.deepcopy(glist)
	for ic in range(len(clusters)):
		slist[ic] = int(swap_ratio*len(clusters[ic]))
		shuffle(clusters[ic])
		temp_list +=clusters[ic][0:slist[ic]]
		del clusters[ic][0:slist[ic]]
	ulist +=temp_list
	for ic in range(len(clusters)):
		shuffle(ulist)
		clusters[ic] +=ulist[0:slist[ic]]
		del ulist[0:slist[ic]]
	return ulist, clusters
		
def swap_clusters_large_NUACC(glist, clusters, swap_ratio):
	from random import shuffle
	import copy
	slist     = [None for im in range(len(clusters))]
	temp_list = []
	ulist = copy.deepcopy(glist)
	for ic in range(len(clusters)):
		slist[ic] = int(swap_ratio*len(clusters[ic]))
		shuffle(clusters[ic])
		temp_list +=clusters[ic][0:slist[ic]]
		del clusters[ic][0:slist[ic]]
	for ic in range(len(clusters)):
		shuffle(ulist)
		clusters[ic] +=ulist[0:slist[ic]]
		del ulist[0:slist[ic]]
	ulist +=temp_list
	shuffle(ulist)
	return ulist, clusters

def even_assignment_alist_to_mclusters(glist, number_of_groups):
	# evenly assign glist to clusters
	import copy
	from random import shuffle
	if number_of_groups >0:
		clusters = [[] for i in range(number_of_groups)]
		ulist = copy.deepcopy(glist)
		nc = 0
		while len(ulist)>0:
			im =  nc%number_of_groups
			shuffle(ulist)
			clusters[im].append(ulist[0])
			del ulist[0]
			nc +=1
		return clusters
	else: return []
	
def fill_no_large_groups_and_unaccounted_to_m_and_rcluster_mpi(\
		unaccounted_list, empty_clusters, clusters, NUACC, NACC):
	global Tracker, Blockdata
	import numpy as np
	N = NUACC + NACC
	m = len(empty_clusters)
	number_of_groups = m + len(clusters)
	avg_size = N//number_of_groups
	if m*avg_size < NUACC - avg_size:
		if Blockdata["myid"] == Blockdata["main_node"]:
			for ic in range(len(clusters)):
				shuffle(clusters[ic])
				if len(clusters[ic])> avg_size:
					unaccounted_list +=clusters[ic][avg_size:]
					del clusters[ic][avg_size:]
			shuffle(unaccounted_list)
		else: unaccounted_list = 0
		unaccounted_list = wrap_mpi_bcast(unaccounted_list, Blockdata["main_node"], MPI_COMM_WORLD)
		clusters         = wrap_mpi_bcast(clusters,         Blockdata["main_node"], MPI_COMM_WORLD)
		tmp_clusters = []
		if m > 0:
			if Blockdata["myid"] == Blockdata["main_node"]:
				for im in range(m):
					cluster, unaccounted_list = \
					   select_fixed_size_cluster_from_alist(unaccounted_list, avg_size//2)
					tmp_clusters.append(cluster)
			else: tmp_clusters = 0
			tmp_clusters     = wrap_mpi_bcast(tmp_clusters,     Blockdata["main_node"], MPI_COMM_WORLD)
			unaccounted_list = wrap_mpi_bcast(unaccounted_list, Blockdata["main_node"], MPI_COMM_WORLD)	 
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
	import random
	import copy
	from random import shuffle
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
		shuffle(ulist)
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
	import numpy as np
	import random
	import copy
	checking_flag = 0
	if Blockdata["myid"] == Blockdata["main_node"]:
		if len(read_text_row(accounted_file)) <= 1:checking_flag = 1
	checking_flag = bcast_number_to_all(checking_flag, Blockdata["main_node"], MPI_COMM_WORLD)
	if checking_flag == 0:
		if Blockdata["myid"] == Blockdata["main_node"]:
			clusters, npart  = split_partition_into_ordered_clusters(read_text_file(accounted_file, -1), False)
			unaccounted_list = np.array(read_text_file(unaccounted_file), dtype = np.int32)
		else: 
			clusters = 0
			unaccounted_list = 0
		clusters          = wrap_mpi_bcast(clusters,         Blockdata["main_node"], MPI_COMM_WORLD)
		unaccounted_list  = wrap_mpi_bcast(unaccounted_list, Blockdata["main_node"], MPI_COMM_WORLD)
		clusters = refilling_global_scheme_mpi(clusters, unaccounted_list.tolist(), \
		   number_of_groups, log_file, swap_ratio)	
		if Blockdata["myid"] == Blockdata["main_node"]:
			dlist, assignment_list = merge_classes_into_partition_list(clusters)
			converted_assignment_list = [ ]
			assignment_list = np.array(assignment_list).transpose()
			for jm in range(2):converted_assignment_list.append(assignment_list[jm].tolist())
		else: converted_assignment_list = 0
		converted_assignment_list = wrap_mpi_bcast(converted_assignment_list, \
		    Blockdata["main_node"], MPI_COMM_WORLD)
	else: # there exits unaccounted
		assignment_list = create_nrandom_lists(unaccounted_file, number_of_groups, 1)#MPI
		converted_assignment_list = [ ]
		assignment_list = np.array(assignment_list[0]).transpose()
		for jm in range(2):converted_assignment_list.append(assignment_list[jm].tolist())
	return converted_assignment_list
	
def patch_to_do_k_means_match_clusters_asg_new(ptp1, ptp2):
	from statistics import k_means_match_clusters_asg_new
	import numpy as np
	# input are clusters
	# patch ad hoc elements to make equal number of classes for two partitions and thus two_way comparison becomes feasible
	patch_elements = []
	if len(ptp1) != len(ptp2):
		max_number = -1
		for im in range(len(ptp1)): 
			ptp1[im]   = np.array(ptp1[im], dtype = np.int32)
			max_number = max(max_number, np.max(ptp1[im]))
		for im in range(len(ptp2)):
		    ptp2[im]  = np.array(ptp2[im], dtype = np.int32)
		    max_number = max(max_number, np.max(ptp2[im]))
		
		if len(ptp1) > len(ptp2):
			ndiff = len(ptp1) - len(ptp2)
			for indiff in range(ndiff):
				l = []
				l.append(max_number+indiff+1)
				patch_elements.append(max_number+indiff+1)
				l = array(l,"int32")
				ptp2.append(l)
		else:
			ndiff = len(ptp2)-len(ptp1)
			for indiff in range(ndiff):
				l = []
				l.append(max_number+indiff+1)
				patch_elements.append(max_number + indiff + 1)
				l = array(l,"int32")
				ptp1.append(l)
	else:
		for im in range(len(ptp1)): 
			ptp1[im]   = np.array(ptp1[im],dtype = np.int32)
		for im in range(len(ptp2)): 
			ptp2[im]   = np.array(ptp2[im],dtype = np.int32)
	newindeces, list_stable, nb_tot_objs = k_means_match_clusters_asg_new(ptp1, ptp2)
	new_list_stable = []
	for a in list_stable:
		a.tolist()
		if len(a)>0: new_list_stable.append(a) # remove empty ones
	return newindeces, new_list_stable, nb_tot_objs, patch_elements

def do_boxes_two_way_comparison_mpi(nbox, input_box_parti1, input_box_parti2, depth, log_main):
	global Tracker, Blockdata
	import json
	import numpy as np
	NT = 1000
	if Blockdata["myid"]==Blockdata["main_node"]:
		stop_generation  = 0
		log_main.add('================================================================================================================')
		log_main.add(' Two-way comparison of generation %d and layer %d computed between two pairs of independent runs: %d and %d.'%(\
		      Tracker["current_generation"],Tracker["depth"],nbox, nbox+1))
		log_main.add('----------------------------------------------------------------------------------------------------------------')
		bad_clustering =  0
		ipair = 0
		core1 = read_text_file(input_box_parti1, -1)
		ptp1, tmp1 = split_partition_into_ordered_clusters(core1, False)
		core2 = read_text_file(input_box_parti2, -1)
		ptp2, tmp2 = split_partition_into_ordered_clusters(core2, False)
	else:
		ptp1  = 0
		ptp2  = 0
		core1 = 0
		core2 = 0
	ptp1 = wrap_mpi_bcast(ptp1, Blockdata["main_node"], MPI_COMM_WORLD)
	ptp2 = wrap_mpi_bcast(ptp2, Blockdata["main_node"], MPI_COMM_WORLD)
	core1 = wrap_mpi_bcast(core1, Blockdata["main_node"], MPI_COMM_WORLD)
	core2 = wrap_mpi_bcast(core2, Blockdata["main_node"], MPI_COMM_WORLD)
	rand_index = compute_rand_index_mpi(core1[0], core2[0])
	if Blockdata["myid"]==Blockdata["main_node"]:
		log_main.add('Rand index of two pairs of independent runs is  %5.4f'%rand_index)
	#### before comparison we do a simulation
	gave, gvar = do_random_groups_simulation_mpi(ptp1, ptp2)
	if Blockdata["myid"]==Blockdata["main_node"]:# sorting data processing
		#####
		msg = 'P0      '
		msg1 ='Group ID'
		length = max(len(ptp1), len(ptp2))
		for im in range(length):
			try:    msg +='{:8d} '.format(len(ptp1[im]))
			except: pass
			msg1 +='{:8d} '.format(im)
		log_main.add(msg1)	
		log_main.add(msg)
		msg = 'P1      '
		for im in range(len(ptp2)): msg +='{:8d} '.format(len(ptp2[im]))
		log_main.add(msg)
		full_list  = np.sort(np.array(core1[1], dtype = np.int32))
		total_data = full_list.shape[0]
		minimum_group_size = total_data
		maximum_group_size = 0
		newindeces, list_stable, nb_tot_objs, patch_elements = \
		   patch_to_do_k_means_match_clusters_asg_new(ptp1, ptp2)
		ratio_unaccounted  = 100. - nb_tot_objs/float(total_data)*100.
		Tracker["current_iter_ratio"] = nb_tot_objs/float(total_data)*100.
		new_list = []
		print_matching_pairs(newindeces, log_main)
		score_list = [ ]
		nclass     = 0
		stat_list  = []
		tmp_list   = []
		log_main.add(' ')
		log_main.add('               Post-matching results.')
		log_main.add('{:>5} {:>8}  {:^8}   {:>15} {:>22}  {:>5}'.format(\
		  'Group', '    size',  ' status ',   'reproducibility', 'random reproducibility', ' std '))
		from math import sqrt
		for index_of_any in range(len(list_stable)):
			any = np.sort(list_stable[index_of_any])
			any.tolist()
			score1 = float(len(any))*100./float(len(ptp1[newindeces[index_of_any][0]]))
			score2 = float(len(any))*100./float(len(ptp2[newindeces[index_of_any][1]]))
			score3 = float((np.intersect1d(ptp1[newindeces[index_of_any][0]], ptp2[newindeces[index_of_any][1]])).size)\
				  /float((np.union1d(ptp1[newindeces[index_of_any][0]], ptp2[newindeces[index_of_any][1]])).size)*100.
			if (score3 > gave[index_of_any]+Tracker["constants"]["random_group_elimination_threshold"]*gvar[index_of_any]):
				score_list.append([score1, score2])
				minimum_group_size = min(minimum_group_size, len(any))
				maximum_group_size = max(maximum_group_size, len(any))
				new_list.append(any)
				nclass +=1
				log_main.add('{:>5} {:>8d}   {:^8}      {:>7.1f}      {:>15.1f}           {:>5.1f}'.format(\
				      index_of_any, len(any),'accepted', score3, \
					    gave[index_of_any], gvar[index_of_any]))
				stat_list.append([score3,  gave[index_of_any], gvar[index_of_any]])
				tmp_list.append(len(any)*(-1))
			else:
				log_main.add('{:>5} {:>8d}   {:^8}      {:>7.1f}      {:>15.1f}           {:>5.1f}'.format(\
				    index_of_any, len(any), 'rejected', score3, \
					   gave[index_of_any], gvar[index_of_any]))
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
	nclass = bcast_number_to_all(nclass, Blockdata["main_node"], MPI_COMM_WORLD)
		
	if nclass == 0:
		### redo two way comparison
		if depth >1:
			if Blockdata["myid"]==Blockdata["main_node"]:
				log_main.add('There are no clusters larger than the user provided minimum group size %d.'%\
				  Tracker["constants"]["minimum_grp_size"])
				log_main.add('Sorting eliminates the smallest group, and continues')
				ptp1, ucluster1 = split_partition_into_ordered_clusters_split_ucluster(core1, False)
				ptp2, ucluster2 = split_partition_into_ordered_clusters_split_ucluster(core2, False)
				newindeces, list_stable, nb_tot_objs, patch_elements = \
				    patch_to_do_k_means_match_clusters_asg_new(ptp1, ptp2)
				if len(list_stable)>3:# No change for two groups
					fake_list = sorted(list_stable, key=len)
					list_stable.remove(fake_list[0])
				list_stable = sorted(list_stable, key=len, reverse = True)
				accounted_list, new_index = merge_classes_into_partition_list(list_stable)
				unaccounted_list = np.sort(np.setdiff1d(full_list, np.array(accounted_list, dtype = np.int32)))
				log_main.add('================================================================================================================\n')
			else: minimum_group_size, maximum_group_size, new_index, unaccounted_list, bad_clustering, stop_generation, stat_list = 0, 0, 0, 0, 0, 0, 0
			new_index          = wrap_mpi_bcast(new_index,               Blockdata["main_node"], MPI_COMM_WORLD)
			unaccounted_list   = wrap_mpi_bcast(unaccounted_list,        Blockdata["main_node"], MPI_COMM_WORLD)
			stat_list          = wrap_mpi_bcast(stat_list,               Blockdata["main_node"], MPI_COMM_WORLD)
			bad_clustering     = bcast_number_to_all(bad_clustering,     Blockdata["main_node"], MPI_COMM_WORLD)
			stop_generation    = bcast_number_to_all(stop_generation,    Blockdata["main_node"], MPI_COMM_WORLD)
			maximum_group_size = bcast_number_to_all(maximum_group_size, Blockdata["main_node"], MPI_COMM_WORLD)
			minimum_group_size = bcast_number_to_all(minimum_group_size, Blockdata["main_node"], MPI_COMM_WORLD)
			return minimum_group_size, maximum_group_size, new_index, unaccounted_list, bad_clustering, stop_generation, stat_list
		else:
			if Blockdata["myid"]==Blockdata["main_node"]:
				bad_clustering = 1
				log_main.add('There are no clusters larger than the user provided minimum group size %d.'%Tracker["constants"]["minimum_grp_size"])
				log_main.add('The reason can be: (1) There are no groups in the data set. 2. Minimum group size set is too large. (3) Desired number of groups K is too large.')
				log_main.add('================================================================================================================\n')
			else: 
				minimum_group_size, maximum_group_size, full_list, bad_clustering, stop_generation, stat_list = 0, 0, 0, 0, 0, 0
			full_list          = wrap_mpi_bcast(full_list,               Blockdata["main_node"], MPI_COMM_WORLD)
			stat_list          = wrap_mpi_bcast(stat_list,               Blockdata["main_node"], MPI_COMM_WORLD)
			bad_clustering     = bcast_number_to_all(bad_clustering,     Blockdata["main_node"], MPI_COMM_WORLD)
			stop_generation    = bcast_number_to_all(stop_generation,    Blockdata["main_node"], MPI_COMM_WORLD)
			maximum_group_size = bcast_number_to_all(maximum_group_size, Blockdata["main_node"], MPI_COMM_WORLD)
			minimum_group_size = bcast_number_to_all(minimum_group_size, Blockdata["main_node"], MPI_COMM_WORLD)
			return minimum_group_size, maximum_group_size, [ ], full_list.tolist(), bad_clustering, stop_generation, stat_list
			
	elif nclass == 1: # Force to stop this generation, and output the cluster; do not do any other box comparison
		if Blockdata["myid"]==Blockdata["main_node"]:
			stop_generation  = 1
			accounted_list, new_index = merge_classes_into_partition_list(new_list)
			unaccounted_list = np.sort(np.setdiff1d(np.array(full_list), np.array(accounted_list, dtype=np.int32)))
			log_main.add('Only one group found. The program will output it and stop executing the current generation.')

			box1_dir =  os.path.join(Tracker["constants"]["masterdir"], "generation_%03d"%Tracker["current_generation"], "layer%d"%Tracker["depth"], "nbox%d"%nbox)
			box2_dir =  os.path.join(Tracker["constants"]["masterdir"], "generation_%03d"%Tracker["current_generation"], "layer%d"%Tracker["depth"], "nbox%d"%(nbox+1))
			gendir   =  os.path.join(Tracker["constants"]["masterdir"], "generation_%03d"%Tracker["current_generation"])
			fout     = open(os.path.join(box1_dir,"freq_cutoff.json"),'r')
			freq_cutoff_dict1 = convert_json_fromunicode(json.load(fout))
			fout.close()
			fout = open(os.path.join(box2_dir,"freq_cutoff.json"),'r')
			freq_cutoff_dict2 = convert_json_fromunicode(json.load(fout))
			fout.close()
			try:
				fout = open(os.path.join(gendir,"freq_cutoff.json"),'r')
				freq_cutoff_dict3 = convert_json_fromunicode(json.load(fout))
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
			fout = open(os.path.join(gendir,"freq_cutoff.json"),'w')
			json.dump(freq_cutoff_dict3, fout)
			fout.close()
			log_main.add('================================================================================================================\n')
		else:
			minimum_group_size, maximum_group_size, new_index, unaccounted_list, \
			   bad_clustering, stop_generation, stat_list = 0, 0, 0, 0, 0, 0, 0
		unaccounted_list   = wrap_mpi_bcast(unaccounted_list,        Blockdata["main_node"], MPI_COMM_WORLD)
		new_index          = wrap_mpi_bcast(new_index,               Blockdata["main_node"], MPI_COMM_WORLD)
		stat_list          = wrap_mpi_bcast(stat_list,               Blockdata["main_node"], MPI_COMM_WORLD)
		bad_clustering     = bcast_number_to_all(bad_clustering,     Blockdata["main_node"], MPI_COMM_WORLD)
		stop_generation    = bcast_number_to_all(stop_generation,    Blockdata["main_node"], MPI_COMM_WORLD)
		maximum_group_size = bcast_number_to_all(maximum_group_size, Blockdata["main_node"], MPI_COMM_WORLD)
		minimum_group_size = bcast_number_to_all(minimum_group_size, Blockdata["main_node"], MPI_COMM_WORLD)
		return minimum_group_size, maximum_group_size, new_index, \
		   unaccounted_list.tolist(), bad_clustering, stop_generation, stat_list
	else:
		if Blockdata["myid"]==Blockdata["main_node"]:
			accounted_list, new_index = merge_classes_into_partition_list(new_list)
			unaccounted_list = np.sort(np.setdiff1d(np.array(full_list, dtype=np.int32),\
			    np.array(accounted_list, dtype=np.int32)))
			log_main.add(' {} {} {} {}'.format('  The number of accounted for images:', \
			    len(accounted_list),'  The number of unaccounted for images:', len(unaccounted_list)))
			log_main.add('  The current minimum group size: %6d\n  The current maximum group size: %6d'%(minimum_group_size, maximum_group_size))
		
			if depth <= 1: # the last layer
				box1_dir =  os.path.join(Tracker["constants"]["masterdir"], "generation_%03d"%Tracker["current_generation"], "layer%d"%Tracker["depth"], "nbox%d"%nbox)
				box2_dir =  os.path.join(Tracker["constants"]["masterdir"], "generation_%03d"%Tracker["current_generation"], "layer%d"%Tracker["depth"], "nbox%d"%(nbox+1))
				gendir   =  os.path.join(Tracker["constants"]["masterdir"], "generation_%03d"%Tracker["current_generation"])
				fout     = open(os.path.join(box1_dir,"freq_cutoff.json"),'r')
				freq_cutoff_dict1 = convert_json_fromunicode(json.load(fout))
				fout.close()
				fout = open(os.path.join(box2_dir,"freq_cutoff.json"),'r')
				freq_cutoff_dict2 = convert_json_fromunicode(json.load(fout))
				fout.close()
				try:
					fout = open(os.path.join(gendir,"freq_cutoff.json"),'r')
					freq_cutoff_dict3 = convert_json_fromunicode(json.load(fout))
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
				fout = open(os.path.join(gendir,"freq_cutoff.json"),'w')
				json.dump(freq_cutoff_dict3, fout)
				fout.close()
			log_main.add('================================================================================================================\n')
		else:
			minimum_group_size, maximum_group_size, new_index, unaccounted_list, \
			   bad_clustering, stop_generation, stat_list = 0, 0, 0, 0, 0, 0, 0
		unaccounted_list   = wrap_mpi_bcast(unaccounted_list,        Blockdata["main_node"], MPI_COMM_WORLD)
		new_index          = wrap_mpi_bcast(new_index,               Blockdata["main_node"], MPI_COMM_WORLD)
		stat_list          = wrap_mpi_bcast(stat_list,               Blockdata["main_node"], MPI_COMM_WORLD)
		bad_clustering     = bcast_number_to_all(bad_clustering,     Blockdata["main_node"], MPI_COMM_WORLD)
		stop_generation    = bcast_number_to_all(stop_generation,    Blockdata["main_node"], MPI_COMM_WORLD)
		maximum_group_size = bcast_number_to_all(maximum_group_size, Blockdata["main_node"], MPI_COMM_WORLD)
		minimum_group_size = bcast_number_to_all(minimum_group_size, Blockdata["main_node"], MPI_COMM_WORLD)
		return minimum_group_size, maximum_group_size, new_index, unaccounted_list.tolist(), \
		   bad_clustering, stop_generation, stat_list

def do_withinbox_two_way_comparison(partition_dir, nbox, nrun, niter):
	global Tracker, Blockdata
	import numpy as np
	log_list = []
	## for single node only
	log_list.append(' ')
	log_list.append('----------------------------------------------------------------------------------------------------------------')	      
	ipair      = 0
	core1      = read_text_file(os.path.join(partition_dir, "partition_%03d.txt"%(2*ipair)),  -1)
	ptp1, tmp1 = split_partition_into_ordered_clusters(core1, False)
	core2      = read_text_file(os.path.join(partition_dir, "partition_%03d.txt"%(2*ipair+1)),-1)
	ptp2, tmp2 = split_partition_into_ordered_clusters(core2, False)
	log_list.append('       Matching of sorting results of two quasi-independent runs')
	# before comparison
	msg = 'P0      '
	msg1 ='Group ID'
	for im in range(len(ptp1)):
		msg  +='{:8d} '.format(len(ptp1[im]))
		msg1 +='{:8d} '.format(im)
	log_list.append(msg1)
	log_list.append(msg)
	msg = 'P1      '
	for im in range(len(ptp2)): msg +='{:8d} '.format(len(ptp2[im]))
	log_list.append(msg)
	full_list          = np.sort(np.array(core1[1], dtype=np.int32))
	total_data         = full_list.shape[0]
	minimum_group_size = total_data
	maximum_group_size = 0
	newindeces, list_stable, nb_tot_objs, patch_elements = \
	          patch_to_do_k_means_match_clusters_asg_new(ptp1, ptp2)
	ratio_unaccounted = 100.-nb_tot_objs/float(total_data)*100.
	Tracker["current_iter_ratio"] = nb_tot_objs/float(total_data)*100.

	log_list.append(' ')
	log_list.append('                        Two-way matching of sorting results.')
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
	score_list = [ ]
	nclass     = 0
	log_list.append('               Post-matching results.')
	log_list.append('{:>8} {:>10} {:>17} {:>8} {:>15}'.format('Group ID', '   size', 'min random size', ' status ',   'reproducibility'))
	current_MGR = get_MGR_from_two_way_comparison(newindeces, ptp1, ptp2, total_data)
	stable_clusters   = []
	selected_clusters = []
	for index_of_any in range(len(list_stable)):
		any = np.sort(list_stable[index_of_any])
		any.tolist()
		score1 = float(len(any))*100./float(len(ptp1[newindeces[index_of_any][0]]))
		score2 = float(len(any))*100./float(len(ptp2[newindeces[index_of_any][1]]))
		score3 = float((np.intersect1d(ptp1[newindeces[index_of_any][0]], ptp2[newindeces[index_of_any][1]])).size)\
			  /float((np.union1d(ptp1[newindeces[index_of_any][0]], ptp2[newindeces[index_of_any][1]])).size)*100.
		if len(any) > current_MGR[index_of_any]:
			score_list.append([score1, score2])
			minimum_group_size = min(minimum_group_size, len(any))
			maximum_group_size = max(maximum_group_size, len(any))
			nclass +=1
			log_list.append('{:>8} {:>10d} {:>10d}        {:>8} {:>8.1f}       '.format(index_of_any, \
			     len(any), current_MGR[index_of_any],'accepted', score3))
			selected_clusters.append(any)
		else:
			log_list.append('{:>8} {:>10d} {:>10d}        {:>8}  {:>8.1f}       '.format(index_of_any, \
			     len(any), current_MGR[index_of_any], 'rejected', score3))		
	accounted_list, new_index = merge_classes_into_partition_list(selected_clusters)
	unaccounted_list = (np.sort(np.setdiff1d(full_list, np.array(accounted_list, dtype=np.int32)))).tolist()
	write_text_row(new_index, os.path.join(partition_dir, "Accounted.txt"))
	write_text_file(unaccounted_list, os.path.join(partition_dir, "Core_set.txt"))
	log_list.append('  The overall reproducibility is %5.1f%%.'%Tracker["current_iter_ratio"])
	log_list.append('  The number of accounted for images: %d.  The number of core set images: %d.'%(len(accounted_list), len(unaccounted_list)))
	log_list.append('  The current minimum group size: %d and the maximum group size: %d.'%(minimum_group_size, maximum_group_size))
	return minimum_group_size, maximum_group_size, selected_clusters, unaccounted_list, Tracker["current_iter_ratio"], len(list_stable), log_list

def compute_rand_index_mpi(inassign1, inassign2):
	# Two partitions always have the same length
	# rand_index = (a + b)/[(n-1)*n/2] ; n is the full length of both partitions
	# a: pairs in the same class in both partitions
	# b: pairs not in the same class in both partitions
	# Journal of the American Statistical Association.
	# American Statistical Association. 66(336): 846-850
	# rand_index should be between [0, 1]; 
	# =0; two partitions different; =1; two partitions identical
	global Tracker, Blockdata
	import numpy as np
	assign1 = np.array(inassign1, dtype=np.int8)
	assign2 = np.array(inassign2, dtype=np.int8)
	num_in_both     = 0
	num_in_neither  = 0
	ntot = len(assign1)
	num_all_pairs = (ntot-1)/2.
	parti_list = []
	for im in range(ntot-1):
		if im%Blockdata["nproc"] == Blockdata["myid"]:
			parti_list.append(im)
	for lm in range(len(parti_list)):
		for im in range(parti_list[lm], parti_list[lm]+1):
			for jm in range(im +1, ntot):
				if (assign1[im] == assign1[jm]) and (assign2[im] == assign2[jm]):
					num_in_both    +=1
				if (assign1[im] != assign1[jm]) and (assign2[im] != assign2[jm]):
					num_in_neither +=1
	nomin = (num_in_both+num_in_neither)/float(ntot)
	del assign1, assign2
	mpi_barrier(MPI_COMM_WORLD)
	if (Blockdata["myid"] != Blockdata["main_node"]):
		wrap_mpi_send(nomin, Blockdata["main_node"], MPI_COMM_WORLD)
	else:
		for iproc in range(Blockdata["nproc"]):
			if (iproc !=Blockdata["main_node"]):
				nomin += wrap_mpi_recv(iproc, MPI_COMM_WORLD)
	nomin = bcast_number_to_all(nomin,  Blockdata["main_node"], MPI_COMM_WORLD)
	return float(nomin)/num_all_pairs
	
def isin(element, test_elements, assume_unique=False, invert=False):
    import numpy as np
    element = np.asarray(element)
    return np.in1d(element, test_elements, assume_unique=assume_unique, \
        invert=invert).reshape(element.shape)

def split_partition_into_ordered_clusters_split_ucluster(partition, input_row_wise = True):
	# group particles  by their cluster ids; take the last one as unaccounted group
	import numpy as np
	clusters   = []
	if input_row_wise: partition = np.array(partition, dtype=np.int32).transpose()
	group_id = np.sort(np.unique(partition[0]))
	if group_id.shape[0] >1: 
		for icluster in range(group_id.shape[0]):
			if icluster <group_id.shape[0] -1: clusters.append((np.sort(partition[1][isin(partition[0],\
			     group_id[icluster])])).tolist())
			else:  ucluster =(np.sort(partition[1][isin(partition[0], group_id[icluster])])).tolist()
		return clusters, ucluster
	else: return [partition[1].tolist()], []
                           
def split_partition_into_ordered_clusters(partition, input_is_row_wise = True):
	# partition column 0 cluster  IDs
	# partition column 1 particle IDs
	# re-index partition as consecutive groups if their cluster IDs are not
	# Extract clusters; sort clusters in ascendent order; transpose partitions
	import numpy as np
	if input_is_row_wise:partition = (np.array(partition).transpose()).tolist()
	np_cluster_ids  = np.array(partition[0], dtype=np.int32)
	np_particle_ids = np.array(partition[1], dtype=np.int32)
	group_id        = np.unique(np_cluster_ids)
	cluster_list    = []
	mask_list       = []
	N = len(partition[0])
	llist = np.full(group_id.shape[0], 0, dtype=np.int32)
	cluster_list = [None for im in range(group_id.shape[0])]
	for ic in np.nditer(group_id, order='C'):
		m  = isin(np_cluster_ids, ic)
		l  = np_particle_ids[m].tolist()
		cluster_list[ic] = l
		mask_list.append(m)
		llist[ic] = -len(l)
	sort_indx = np.argsort(llist)
	new_clusters_ids = np.full(N, -1, dtype=np.int32)
	new_clusters     = []
	for ic in range(group_id.shape[0]):
		new_clusters.append(cluster_list[sort_indx[ic]])
		np.place(new_clusters_ids, mask_list[sort_indx[ic]], ic)
	partition[0] = new_clusters_ids
	return new_clusters, (np.array(partition, dtype=np.int32).transpose()).tolist()

def merge_classes_into_partition_list(classes_list):
	import numpy as np
	import copy
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
		indx = np.argsort(np.array(size_list)) # GID in ascendent order
		indx.tolist()
		cluster_ids = np.full(parti_list.shape[0], -1, dtype=np.int32)
		for im in range(len(classes_list)):
			np.place(cluster_ids, isin(parti_list, \
			np.array(classes_list[indx[im]])), np.int32(im))
		return parti_list.tolist(), (np.array([cluster_ids, \
		     parti_list],dtype=np.int32 ).transpose()).tolist()
		
def get_sorting_parti_from_stack(data):
	global Tracker, Blockdata
	from utilities    import wrap_mpi_bcast
	from applications import MPI_start_end
	import numpy as np
	#data: a subset on a CPU
	if Blockdata["myid"] == Blockdata["main_node"]:
		total_plist = np.full(Tracker["total_stack"], -1, dtype=np.int32)
	else: total_plist = 0
	for myproc in range(Blockdata["nproc"]):
		image_start,image_end = MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], myproc)
		plist = 0
		if Blockdata["myid"] == myproc:
			plist = np.full(len(data), -1, dtype=np.int32)
			for im in range(len(data)):
				plist[im] = data[im].get_attr("group")
		plist = wrap_mpi_bcast(plist, myproc, MPI_COMM_WORLD)
		if Blockdata["myid"] == Blockdata["main_node"]:
			total_plist[image_start:image_end] = plist
		mpi_barrier(MPI_COMM_WORLD)
	total_plist = wrap_mpi_bcast(total_plist, Blockdata["main_node"])
	return total_plist

### dmatrix and refangles partition			
def number_of_cones_to_delta(number_of_cones):

	if( number_of_cones == 1):  return Blockdata["symclass"][1][3] + 0.5, 1
	else:
		if( Blockdata["symclass"].sym[0] == "c"):
			t2 = 89.0
			for i in range(92,0,-1):
				a = Blockdata["symclass"].even_angles(i, theta1=1.0, theta2=t2)
				a += [[(q[0]+90.0)%360., 180.-q[1],0] for q in a]
				nren = len(a)
				if(nren>number_of_cones): return float(i)#, nren
			while(True):
				i /= 2.0
				a = Blockdata["symclass"].even_angles(i, theta1=1.0, theta2=t2)
				a += [[(q[0]+90.0)%360., 180.-q[1],0] for q in a]
				nren = len(a)
				if(nren>number_of_cones): return float(i)#, nren

		else:
			#t2 = Blockdata["symclass"].brackets[1][3] + 0.5
			for i in range(int(t2+1),0,-1):
				nren = len(Blockdata["symclass"].even_angles(i, theta1=1.0))
				if(nren>number_of_cones): return float(i)#, nren
			while(True):
				i /= 2.0
				nren = len(Blockdata["symclass"].even_angles(i, theta1=1.0))
				if(nren>number_of_cones): return float(i)#, nren

		ERROR("number_of_cones_to_delta","should not be here",0)
		return -1, 0

def get_angle_step_and_orien_cones_mpi(params_in, partids_in, angstep):
	global Tracker, Blockdata
	from applications import MPI_start_end
	from utilities    import wrap_mpi_recv, wrap_mpi_bcast, wrap_mpi_send, bcast_number_to_all
	from utilities    import read_text_row, read_text_file
	from fundamentals import symclass
	if Blockdata["main_node"] == Blockdata["myid"]:
		params  = read_text_row(params_in)
		partids = read_text_file(partids_in, -1)
		if len(partids) == 1: partids = partids[0]
		else:                 partids = partids[1]
		subparams = []
		for im in range(len(partids)): 
			subparams.append(params[partids[im]])
		from utilities import angular_occupancy
		occu, eah = angular_occupancy(subparams, angstep, Blockdata["symclass"].sym)
		ptls_in_orien_cones = []
		reassign_list        = []
		for l,q in enumerate(occu):
			if len(q)<Tracker["min_orien_cone_size"]: reassign_list +=q
			else:      ptls_in_orien_cones.append(q)
		import random
		from   random import shuffle
		shuffle(reassign_list)
		for a in reassign_list:
			img = random.randint(0, len(ptls_in_orien_cones) - 1)
			ptls_in_orien_cones[img].append(a)
		del reassign_list
		for img in range(len(ptls_in_orien_cones)):
			tmp = sorted(ptls_in_orien_cones[img])
			ptls_in_orien_cones[img][:] = tmp[:]
	else:
		ptls_in_orien_cones = 0
	ptls_in_orien_cones = wrap_mpi_bcast( ptls_in_orien_cones, Blockdata["main_node"], MPI_COMM_WORLD)
	return ptls_in_orien_cones

##### ================= orientation groups	
def compare_two_iterations(assignment1, assignment2):
	# compare two assignments during clustering, either iteratively or independently
	import numpy as np
	if len(assignment1) !=len(assignment2):
		ERROR("Two assignments do not have the same length", "compare_two_iterations", 0, 0)
	else:
		N = len(assignment1)
		full_list  = np.array(range(N), dtype = np.int32)
	res1 = []
	assignment1 = np.array(assignment1,  dtype = np.int32)
	group1_ids  = np.unique(assignment1)
	for im in range(group1_ids.shape[0]):
		res1.append(np.sort(full_list[isin(assignment1, group1_ids[im])]))
	res2 = []
	assignment2 = np.array(assignment2, dtype = np.int32)
	group2_ids  = np.unique(assignment2)
	for im in range(group2_ids.shape[0]):
		res2.append(np.sort(full_list[isin(assignment2, group2_ids[im])]))
	newindeces, list_stable, nb_tot_objs = k_means_match_clusters_asg_new(res1, res2)
	del res1
	del res2
	return float(nb_tot_objs)/N, newindeces, list_stable
	
def update_data_assignment(cdata, rdata, assignment, proc_list, nosmearing, myid):
	import numpy as np
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
			if value == None and value_refinement != None: Tracker[key] = value_refinement
		except:
			if Blockdata["myid"] == Blockdata["main_node"]: print(key, \
			   " in sorting set as ", value, ", while in refinement, it is set as ", value_refinement)
	return
'''
def print_dict(dict, theme):
	spaces = "                    "
	exclude = ["nodes", "yr", "output", "shared_comm", "bckgnoise", \
	   "myid", "myid_on_node", "accumulatepw", \
	     "chunk_dict", "PW_dict", "full_list", "rshifts", "refang",\
	         "chunk_dict", "PW_dict", "full_list", "rshifts", \
	        "refang", "sorting_data_list", "partition", \
	          "constants", "random_assignment"]
	for key, value in sorted( dict.items() ):
		pt = True
		for ll in exclude:
			if(key == ll):
				pt = False
				break
		if pt:  print("                    => ", key+spaces[len(key):],":  ",value)
'''
# 
# - "Tracker" (dictionary) object
#   Keeps the current state of option settings and dataset 
#   (i.e. particle stack, reference volume, reconstructed volume, and etc)
#   Each iteration is allowed to add new fields/keys
#   if necessary. This happes especially when type of 3D Refinement or metamove changes.
#   Conceptually, each iteration will be associated to a specific Tracker state.
#   Therefore, the list of Tracker state represents the history of process.
#
#   This can be used to restart process from an arbitrary iteration.
#                     rec3d for sorting
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
			ovol = reshape_1d(cfsc, nr, 2*nr)
			limitres = 2*min(limitres, Tracker["maxfrad"])  # 2 on account of padding, which is always on
			maxr2 = limitres**2
			for i in range(limitres+1, len(ovol), 1):   ovol[i] = 0.0
			ovol[0] = 1.0
			it = model_blank(2*nr)
			for i in range(2*nr):  it[i] = ovol[i]
			del ovol
			#  Do not regularize first four
			for i in range(5):  treg[i] = 0.0
			Util.reg_weights(tweight, treg, it)
			del it
		else:
			limitres = 2*min(Tracker["constants"]["nnxo"]//2, Tracker["maxfrad"])
			maxr2 = limitres**2
		#  Iterative weights
		if( Tracker["constants"]["symmetry"] != "c1" ):
			tvol    = tvol.symfvol(Tracker["constants"]["symmetry"], limitres)
			tweight = tweight.symfvol(Tracker["constants"]["symmetry"], limitres)
		Util.divclreal(tvol, tweight, 1.0e-5)

	else:
		tvol = model_blank(1)
		tweight = model_blank(1)
		nz    = 0
		ny    = 0
		nx    = 0
		maxr2 = 0
	nx    = bcast_number_to_all(nx, source_node = 0,    mpi_comm = Blockdata["shared_comm"])
	ny    = bcast_number_to_all(ny, source_node = 0,    mpi_comm = Blockdata["shared_comm"])
	nz    = bcast_number_to_all(nz, source_node = 0,    mpi_comm = Blockdata["shared_comm"])
	maxr2 = bcast_number_to_all(maxr2, source_node = 0, mpi_comm = Blockdata["shared_comm"])
	'''
	vol_data = get_image_data(tvol)
	we_data =  get_image_data(tweight)
	#  tvol is overwritten, meaning it is also an output
	n_iter = 10


	ifi = mpi_iterefa( vol_data.__array_interface__['data'][0] ,  we_data.__array_interface__['data'][0] , nx, ny, nz, maxr2, \
			Tracker["constants"]["nnxo"], Blockdata["myid_on_node"], color, Blockdata["no_of_processes_per_group"],  Blockdata["shared_comm"], n_iter)
	'''
	if( Blockdata["myid_on_node"] == 0 ):
		#  Either pad or window in F space to 2*nnxo
		nx = tvol.get_ysize()
		if( nx > 2*Tracker["constants"]["nnxo"]):
			tvol = fdecimate(tvol, 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], \
			  2*Tracker["constants"]["nnxo"], False, False)
		elif(nx < 2*Tracker["constants"]["nnxo"]):
			tvol = fpol(tvol, 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], \
			   2*Tracker["constants"]["nnxo"], RetReal = False, normalize = False)
		tvol = fft(tvol)
		tvol = cyclic_shift(tvol,Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
		tvol = Util.window(tvol, Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
		tvol.div_sinc(1)
		tvol = cosinemask(tvol, Tracker["constants"]["nnxo"]//2-1,5, None)
		return tvol
	else:  return None
	
def steptwo_mpi_reg(tvol, tweight, treg, cfsc = None, cutoff_freq = 0.45, aa = 0.01, regularized = True, color = 0):
	global Tracker, Blockdata
	cutoff_freq2, aa = estimate_tanhl_params(cutoff_freq, aa, 2*Tracker["constants"]["nnxo"])
	if( Blockdata["color"] != color ):return model_blank(1)  #  This should not be executed if called properly
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
			ovol = reshape_1d(cfsc, nr, 2*nr)
			limitres = 2*min(limitres, Tracker["maxfrad"])  # 2 on account of padding, which is always on
			maxr2 = limitres**2
			for i in range(limitres+1, len(ovol), 1):   ovol[i] = 0.0
			ovol[0] = 1.0
			it = model_blank(2*nr)
			for i in range(2*nr):  it[i] = ovol[i]
			del ovol
			#  Do not regularize first four
			for i in range(5):  treg[i] = 0.0
			Util.reg_weights(tweight, treg, it)
			del it
		else:
			limitres = 2*min(Tracker["constants"]["nnxo"]//2, Tracker["maxfrad"])
			maxr2 = limitres**2
		#  Iterative weights
		if( Tracker["constants"]["symmetry"] != "c1" ):
			tvol    = tvol.symfvol(Tracker["constants"]["symmetry"], limitres)
			tweight = tweight.symfvol(Tracker["constants"]["symmetry"], limitres)
		Util.divclreal(tvol, tweight, 1.0e-5)
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
	'''
	vol_data = get_image_data(tvol)
	we_data =  get_image_data(tweight)
	#  tvol is overwritten, meaning it is also an output
	n_iter = 10

	ifi = mpi_iterefa( vol_data.__array_interface__['data'][0] ,  we_data.__array_interface__['data'][0] , nx, ny, nz, maxr2, \
			Tracker["constants"]["nnxo"], Blockdata["myid_on_node"], color, Blockdata["no_of_processes_per_group"],  Blockdata["shared_comm"], n_iter)	
	'''
	if( Blockdata["myid_on_node"] == 0 ):
		from filter import  filt_tanl
		#  Either pad or window in F space to 2*nnxo
		if not cfsc: tvol = filt_tanl(tvol, cutoff_freq2, aa)
		nx = tvol.get_ysize()
		if( nx > 2*Tracker["constants"]["nnxo"]):
			tvol = fdecimate(tvol, 2*Tracker["constants"]["nnxo"], \
			    2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"], False, False)
		elif(nx < 2*Tracker["constants"]["nnxo"]):
			tvol = fpol(tvol, 2*Tracker["constants"]["nnxo"], 2*Tracker["constants"]["nnxo"],\
			  2*Tracker["constants"]["nnxo"], RetReal = False, normalize = False)
		tvol = fft(tvol)
		tvol = cyclic_shift(tvol,Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
		tvol = Util.window(tvol, Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
		tvol.div_sinc(1)
		tvol = cosinemask(tvol, Tracker["constants"]["nnxo"]//2-1,5, None)
		return tvol
	else:  return None
####=====-----------	
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
	refvol  = model_blank(target_size)
	refvol.set_attr("fudge", 1.0)
	if CTF: do_ctf = 1
	else: do_ctf = 0
	fftvol = EMData()
	weight = EMData()
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, \
	   "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = Reconstructors.get( "nn4_ctfw", params )
	r.setup()
	#if norm_per_particle == None: norm_per_particle = len(prjlist)*[1.0]
	for im in range(len(prjlist)):
		phi, theta, psi, s2x, s2y = get_params_proj(prjlist[im], xform = "xform.projection") # shifts are already applied 
		if random_subset == 2:
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
				if prjlist[im].get_attr("is_complex")==0: prjlist[im] = fft(prjlist[im])
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
	try:    qt = prjlist[0].get_attr("qt")
	except: qt = 1.0
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, "symmetry":"c1", \
	   "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = Reconstructors.get( "nn4_ctfw", params )
	r.setup()
	for im in range(len(prjlist)):
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

def do3d_sorting(procid, data, myid, mpi_comm = -1, log_main = None):
	global Tracker, Blockdata
	if (mpi_comm == -1): mpi_comm = MPI_COMM_WORLD
	if(procid == 0):
		if(Blockdata["no_of_groups"] >1):
			if(Blockdata["myid"] == Blockdata["nodes"][procid]):
				if os.path.exists(os.path.join(Tracker["directory"], "tempdir")): print("tempdir exists")
				else:
					try: os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
					except:  log_main.add("tempdir exists")
		else:
			if myid == Blockdata["main_node"]:
				if not os.path.exists(os.path.join(Tracker["directory"],"tempdir")): 
					try: os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
					except: log_main.add("tempdir exists")
				else: log_main.add("tempdir exists")
	mpi_barrier(mpi_comm)
	
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
	mpi_barrier(mpi_comm)
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
					tvol.write_image(os.path.join(Tracker["directory"], "tempdir", "tvol_%d_%d.hdf"%(procid, index_of_groups)))
					tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%d_%d.hdf"%(procid, index_of_groups)))
					trol.write_image(os.path.join(Tracker["directory"], "tempdir", "trol_%d_%d.hdf"%(procid, index_of_groups)))
				mpi_barrier(MPI_COMM_WORLD)
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
					send_EMData(tvol,    Blockdata["last_node"], tag, MPI_COMM_WORLD)
					send_EMData(tweight, Blockdata["last_node"], tag, MPI_COMM_WORLD)
					send_EMData(trol,    Blockdata["last_node"], tag, MPI_COMM_WORLD)
					
				elif Blockdata["myid"] == Blockdata["last_node"]:
					tag =7007
					tvol    = recv_EMData(Blockdata["nodes"][procid], tag, MPI_COMM_WORLD)
					tweight = recv_EMData(Blockdata["nodes"][procid], tag, MPI_COMM_WORLD)
					trol    = recv_EMData(Blockdata["nodes"][procid], tag, MPI_COMM_WORLD)
					tvol.write_image(os.path.join(Tracker["directory"],    "tempdir", "tvol_%d_%d.hdf"%(procid,    index_of_groups)))
					tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%d_%d.hdf"%(procid, index_of_groups)))
					trol.write_image(os.path.join(Tracker["directory"],    "tempdir", "trol_%d_%d.hdf"%(procid,    index_of_groups)))
				mpi_barrier(MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
		mpi_barrier(MPI_COMM_WORLD)
	return
	
def do3d_sorting_groups_trl_iter(data, iteration):
	global Tracker, Blockdata
	from utilities import get_im, write_text_row, bcast_number_to_all, wrap_mpi_bcast
	keepgoing = 1
	if(Blockdata["myid"] == Blockdata["last_node"]):
		if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
			os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
	mpi_barrier(MPI_COMM_WORLD)
	do3d_sorting_group_insertion(data)
	#####
	if Blockdata["no_of_groups"]>1:
		sub_main_node_list = [-1 for i in range(Blockdata["no_of_groups"])]
		for index_of_colors in range(Blockdata["no_of_groups"]):
			for iproc in range(Blockdata["nproc"]-1):
				if Blockdata["myid"]== iproc:
					if Blockdata["color"] == index_of_colors and Blockdata["myid_on_node"] == 0:
						sub_main_node_list[index_of_colors] = Blockdata["myid"]
					wrap_mpi_send(sub_main_node_list, Blockdata["last_node"], MPI_COMM_WORLD)
				if Blockdata["myid"] == Blockdata["last_node"]:
					dummy = wrap_mpi_recv(iproc, MPI_COMM_WORLD)
					for im in range(len(dummy)):
						if dummy[im]>-1: sub_main_node_list[im] = dummy[im]
				mpi_barrier(MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
		mpi_barrier(MPI_COMM_WORLD)
		wrap_mpi_bcast(sub_main_node_list, Blockdata["last_node"], MPI_COMM_WORLD)
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
					tvol2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0_%d.hdf")%index_of_group)
					tweight2 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_0_%d.hdf")%index_of_group)
					treg2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0_%d.hdf")%index_of_group)
					tag         = 7007
					send_EMData(tvol2, sub_main_node_list[index_of_colors],    tag, MPI_COMM_WORLD)
					send_EMData(tweight2, sub_main_node_list[index_of_colors], tag, MPI_COMM_WORLD)
					send_EMData(treg2, sub_main_node_list[index_of_colors],    tag, MPI_COMM_WORLD)
				elif (Blockdata["myid"] == sub_main_node_list[index_of_colors]):
					tag      = 7007
					tvol2    = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
					tweight2 = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
					treg2    = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
				mpi_barrier(MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
		
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				Tracker["maxfrad"] = Tracker["nxinit"]//2
				if Blockdata["color"] == index_of_colors:
					if( Blockdata["myid_on_node"] != 0):
						tvol2 		= model_blank(1)
						tweight2 	= model_blank(1)
						treg2		= model_blank(1)
					tvol2 = steptwo_mpi(tvol2, tweight2, treg2, Tracker["constants"]["fsc_curve"], \
					   True, color = index_of_colors) # has to be False!!!
					del tweight2, treg2
				mpi_barrier(Blockdata["shared_comm"])
			mpi_barrier(MPI_COMM_WORLD)
		
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				if (Blockdata["color"] == index_of_colors) and (Blockdata["myid_on_node"] == 0):
					tag = 7007
					send_EMData(tvol2, Blockdata["last_node"], tag, MPI_COMM_WORLD)
				elif(Blockdata["myid"] == Blockdata["last_node"]):
					tag = 7007
					tvol2 = recv_EMData(sub_main_node_list[index_of_colors], tag, MPI_COMM_WORLD)
					tvol2.write_image(os.path.join(Tracker["directory"], \
					  "vol_unfiltered_0_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
					del tvol2
				mpi_barrier(MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
		
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				if(Blockdata["myid"] == Blockdata["last_node"]):
					tvol2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_1_%d.hdf")%index_of_group)
					tweight2 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_1_%d.hdf")%index_of_group)
					treg2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "trol_1_%d.hdf"%index_of_group))
					tag      = 7007
					send_EMData(tvol2,    sub_main_node_list[index_of_colors], tag, MPI_COMM_WORLD)
					send_EMData(tweight2, sub_main_node_list[index_of_colors], tag, MPI_COMM_WORLD)
					send_EMData(treg2,    sub_main_node_list[index_of_colors], tag, MPI_COMM_WORLD)
				elif (Blockdata["myid"] == sub_main_node_list[index_of_colors]):
					tag      = 7007
					tvol2    = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
					tweight2 = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
					treg2    = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
				mpi_barrier(MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
		
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				Tracker["maxfrad"] = Tracker["nxinit"]//2
				if Blockdata["color"] == index_of_colors:
					if( Blockdata["myid_on_node"] != 0):
						tvol2 		= model_blank(1)
						tweight2 	= model_blank(1)
						treg2		= model_blank(1)
					tvol2 = steptwo_mpi(tvol2, tweight2, treg2, \
					  Tracker["constants"]["fsc_curve"], True, color = index_of_colors) # has to be False!!!
					del tweight2, treg2
					#if( Blockdata["myid_on_node"] == 0):
					#	tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
				mpi_barrier(Blockdata["shared_comm"])
			mpi_barrier(MPI_COMM_WORLD)
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				if (Blockdata["color"] == index_of_colors) and (Blockdata["myid_on_node"] == 0):
					tag = 7007
					send_EMData(tvol2, Blockdata["last_node"], tag, MPI_COMM_WORLD)
				elif(Blockdata["myid"] == Blockdata["last_node"]):
					tag = 7007
					tvol2 = recv_EMData(sub_main_node_list[index_of_colors], tag, MPI_COMM_WORLD)
					tvol2.write_image(os.path.join(Tracker["directory"], \
					  "vol_unfiltered_1_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
					del tvol2
				mpi_barrier(MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
	else:
		Tracker["maxfrad"] = Tracker["nxinit"]//2
		for index_of_group in range(Tracker["number_of_groups"]):
			for iprocid in range(2):
				if(Blockdata["myid"] == Blockdata["last_node"]):
					tvol2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_%d_%d.hdf")%(iprocid, index_of_group))
					tweight2 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_%d_%d.hdf")%(iprocid, index_of_group))
					treg2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "trol_%d_%d.hdf"%(iprocid, index_of_group)))
					tag      = 7007
					send_EMData(tvol2, Blockdata["main_node"],    tag, MPI_COMM_WORLD)
					send_EMData(tweight2, Blockdata["main_node"], tag, MPI_COMM_WORLD)
					send_EMData(treg2, Blockdata["main_node"],    tag, MPI_COMM_WORLD)
				elif (Blockdata["myid"] == Blockdata["main_node"]):
					tag      = 7007
					tvol2    = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
					tweight2 = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
					treg2    = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
				mpi_barrier(MPI_COMM_WORLD)
				if( Blockdata["myid"] != Blockdata["main_node"]):
					tvol2 		= model_blank(1)
					tweight2 	= model_blank(1)
					treg2		= model_blank(1)
				tvol2 = steptwo_mpi(tvol2, tweight2, treg2, None, False, color = 0) # has to be False!!!
				del tweight2, treg2
				if( Blockdata["myid"] == Blockdata["main_node"]):
					#tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
					tvol2.write_image(os.path.join(Tracker["directory"], \
					  "vol_unfiltered_%d_grp%03d_iter%03d.hdf"%(iprocid, index_of_group,iteration)))
				mpi_barrier(MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
	keepgoing = bcast_number_to_all(keepgoing, source_node = Blockdata["main_node"], mpi_comm = MPI_COMM_WORLD) # always check 
	Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"])
	if keepgoing == 0: ERROR("do3d_sorting_groups_trl_iter  %s"%os.path.join(Tracker["directory"], \
	   "tempdir"),"do3d_sorting_groups_trl_iter", 1, Blockdata["myid"]) 
	return
	
# Three ways of importing refinement results
def get_input_from_sparx_ref3d(log_main):# case one
	# import SPARX results
	global Tracker, Blockdata
	from utilities import write_text_file, write_text_row, bcast_number_to_all, get_im, \
	      bcast_number_to_all
	import numpy as np
	import json
	from  shutil import copyfile
	from  string import split, atoi
	import_from_sparx_refinement = 1
	selected_iter                = 0
	Tracker_refinement           = 0
	checking_flag                = 0
	if Blockdata["myid"] == Blockdata["main_node"]:
		if not os.path.exists (Tracker["constants"]["refinement_dir"]): checking_flag = 1
	checking_flag = bcast_number_to_all(checking_flag, Blockdata["main_node"], MPI_COMM_WORLD)
	if checking_flag ==1:
		ERROR("MERIDIEN refinement directory does not exist", "get_input_from_sparx_ref3d", 1, Blockdata["myid"])
	if Blockdata["myid"] == Blockdata["main_node"]:
		if Tracker["constants"]["niter_for_sorting"] == -1: # take the best solution to do sorting
			niter_refinement = 0
			while os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], \
			   "main%03d"%niter_refinement)) and os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"],\
			     "main%03d"%niter_refinement, "Tracker_%03d.json"%niter_refinement)): niter_refinement +=1
			niter_refinement -=1
			if niter_refinement !=0:
				fout = open(os.path.join(Tracker["constants"]["refinement_dir"], \
				   "main%03d"%niter_refinement, "Tracker_%03d.json"%niter_refinement),'r')
				Tracker_refinement 	= convert_json_fromunicode(json.load(fout))
				fout.close()
				selected_iter = Tracker_refinement["constants"]["best"]
			else: import_from_sparx_refinement = 0
		else:
			try:
				fout = open(os.path.join(Tracker["constants"]["refinement_dir"],\
				  "main%03d"%Tracker["constants"]["niter_for_sorting"], \
				"Tracker_%03d.json"%Tracker["constants"]["niter_for_sorting"]),'r')
				Tracker_refinement	= convert_json_fromunicode(json.load(fout))
				fout.close()
				selected_iter = Tracker["constants"]["niter_for_sorting"]
			except:	import_from_sparx_refinement = 0
	else: selected_iter = -1
	selected_iter = bcast_number_to_all(selected_iter, Blockdata["main_node"], MPI_COMM_WORLD)
	import_from_sparx_refinement = bcast_number_to_all(import_from_sparx_refinement, source_node = Blockdata["main_node"])
	if import_from_sparx_refinement == 0:
		ERROR("The best solution not found","get_input_from_sparx_ref3d", 1, Blockdata["myid"])
		from mpi import mpi_finalize
		mpi_finalize()
		exit()
	Tracker_refinement = wrap_mpi_bcast(Tracker_refinement, Blockdata["main_node"], communicator = MPI_COMM_WORLD)
	# Check orgstack, set correct path
	if Blockdata["myid"] == Blockdata["main_node"]:
		refinement_dir_path, refinement_dir_name = os.path.split(Tracker["constants"]["refinement_dir"])	
		if Tracker_refinement["constants"]["stack"][0:4]=="bdb:":
			refinement_stack = "bdb:" + os.path.join(refinement_dir_path, Tracker_refinement["constants"]["stack"][4:])
		else:
			refinement_stack = os.path.join(refinement_dir_path, Tracker_refinement["constants"]["stack"]) # very rare case		
		if not Tracker["constants"]["orgstack"]: Tracker["constants"]["orgstack"] = refinement_stack
		try:    image = get_im(Tracker["constants"]["orgstack"], 0)
		except: import_from_sparx_refinement = 0
		try:    total_stack = EMUtil.get_image_count(Tracker["constants"]["orgstack"])
		except: total_stack = 0
	else: total_stack = 0
	import_from_sparx_refinement = bcast_number_to_all(import_from_sparx_refinement, source_node = Blockdata["main_node"])
	if import_from_sparx_refinement == 0:
		ERROR("The data stack is not accessible","get_input_from_sparx_ref3d",1, Blockdata["myid"])
	total_stack = bcast_number_to_all(total_stack, source_node = Blockdata["main_node"])			
	Tracker["constants"]["total_stack"] = total_stack
	
	# Now copy relevant refinement files to sorting directory:
	if Blockdata["myid"] == Blockdata["main_node"]:
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], \
		    "main%03d"%selected_iter, "params_%03d.txt"%selected_iter)):
			copyfile( os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iter, \
			 "params_%03d.txt"%selected_iter), os.path.join(Tracker["constants"]["masterdir"], "sparx_refinement_params.txt"))
		else: import_from_sparx_refinement = 0
		Tracker["constants"]["selected_iter"] = selected_iter
	import_from_sparx_refinement = bcast_number_to_all(import_from_sparx_refinement, source_node = Blockdata["main_node"])
	if import_from_sparx_refinement == 0:
		ERROR("The parameter file of the best solution is not accessible", "get_input_from_sparx_ref3d", 1, Blockdata["myid"])
		
	if Blockdata["myid"] == Blockdata["main_node"]:
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iter, "bckgnoise.hdf")):
			copyfile(os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iter, "bckgnoise.hdf"),\
			os.path.join(Tracker["constants"]["masterdir"], "bckgnoise.hdf"))
		else:
			import_from_sparx_refinement == 0
			for search_iter in range(selected_iter-1, 0, -1):
				if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%search_iter, "bckgnoise.hdf")):
					copyfile(os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%search_iter, \
					"bckgnoise.hdf"), os.path.join(Tracker["constants"]["masterdir"], "bckgnoise.hdf"))
					import_from_sparx_refinement = 1
					break
	import_from_sparx_refinement = bcast_number_to_all(import_from_sparx_refinement, source_node = Blockdata["main_node"])
	
	if import_from_sparx_refinement == 0:
		Tracker["bckgnoise"] = None
		if Blockdata["myid"] == Blockdata["main_node"]:	print("Noise file is not found. However, we continue")
	else: Tracker["bckgnoise"] = os.path.join(Tracker["constants"]["masterdir"], "bckgnoise.hdf")
	
	import_from_sparx_refinement = 1	
	if Blockdata["myid"] == Blockdata["main_node"]:
		fifi = os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%selected_iter, "driver_%03d.txt"%selected_iter)
		if os.path.exists(fifi):
			copyfile(fifi, os.path.join(Tracker["constants"]["masterdir"], "fsc_global.txt"))
		else: import_from_sparx_refinement = 0
		
		if import_from_sparx_refinement:
			fsc_curve = read_text_row(os.path.join(Tracker["constants"]["masterdir"], "fsc_global.txt"))
		Tracker["constants"]["fsc_curve"] = \
		  read_text_file(os.path.join(Tracker["constants"]["masterdir"], "fsc_global.txt"))
		
	import_from_sparx_refinement = bcast_number_to_all(import_from_sparx_refinement, source_node = Blockdata["main_node"])
	if import_from_sparx_refinement == 0:
		ERROR("The driver FSC of the best solution is not accessible","get_input_from_sparx_ref3d", 1, Blockdata["myid"])
	if Blockdata["myid"] == Blockdata["main_node"]:
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main000/indexes_000.txt")):
			copyfile(os.path.join(Tracker["constants"]["refinement_dir"], "main000/indexes_000.txt"), \
			os.path.join(Tracker["constants"]["masterdir"], "indexes.txt"))
		else: import_from_sparx_refinement = 0
	import_from_sparx_refinement = bcast_number_to_all(import_from_sparx_refinement, source_node = Blockdata["main_node"])
	
	if import_from_sparx_refinement == 0:
		ERROR("The index file of the best solution is not accessible","get_input_from_sparx_ref3d", 1, Blockdata["myid"])
	if Blockdata["myid"] == Blockdata["main_node"]:
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main000/chunk_0_000.txt")):
			copyfile( os.path.join(Tracker["constants"]["refinement_dir"], "main000/chunk_0_000.txt"), \
			os.path.join(Tracker["constants"]["masterdir"], "chunk_0.txt"))
		else: import_from_sparx_refinement == 0
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main000/chunk_1_000.txt")):
			copyfile(os.path.join(Tracker["constants"]["refinement_dir"], "main000/chunk_1_000.txt"), \
			os.path.join(Tracker["constants"]["masterdir"], "chunk_1.txt"))
		else: import_from_sparx_refinement == 0
			
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main000/particle_groups_0.txt")):
			copyfile(os.path.join(Tracker["constants"]["refinement_dir"], "main000/particle_groups_0.txt"), \
			os.path.join(Tracker["constants"]["masterdir"], "particle_groups_0.txt"))
		else: import_from_sparx_refinement == 0
			
		if os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main000/particle_groups_1.txt")):
			copyfile( os.path.join(Tracker["constants"]["refinement_dir"], "main000/particle_groups_1.txt"), \
			os.path.join(Tracker["constants"]["masterdir"], "particle_groups_1.txt"))
		else:import_from_sparx_refinement == 0
	import_from_sparx_refinement = bcast_number_to_all(import_from_sparx_refinement, source_node = Blockdata["main_node"])
	if import_from_sparx_refinement == 0:
		ERROR("The chunk files and partice group files are not accessible","get_input_from_sparx_ref3d",1, Blockdata["myid"])
	
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
			Tracker["nxinit_refinement"]    = Tracker_refinement["nxinit"] #Sphire window size
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
				copyfile( os.path.join(refinement_dir_path, Tracker_refinement["constants"]["mask3D"]), \
				os.path.join(Tracker["constants"]["masterdir"], refinement_mask3D_file))
				Tracker["constants"]["mask3D"] = os.path.join(Tracker["constants"]["masterdir"], refinement_mask3D_file)
	else: 
		update_sym  = 0
		Tracker     = 0 
	Tracker    = wrap_mpi_bcast(Tracker,         Blockdata["main_node"], MPI_COMM_WORLD)
	update_sym = bcast_number_to_all(update_sym, Blockdata["main_node"], MPI_COMM_WORLD)
	
	if( update_sym == 1):
		Blockdata["symclass"] = symclass(Tracker["constants"]["symmetry"])
		#  Default was set to -1, in which case adjust fo symmetry, else use input directly
		if(Tracker["constants"]["orientation_cones"] < 3):
			Tracker["constants"]["orientation_cones"] = \
			    max(3, int(round(20.0/sqrt(Blockdata["symclass"].nsym))))
	
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
	for index_of_element in range(len(chunk_one)): 
		chunk_dict[chunk_one[index_of_element]] = 0
		group_dict[chunk_one[index_of_element]] = chunk_one_group[index_of_element]
	for index_of_element in range(len(chunk_two)): 
		chunk_dict[chunk_two[index_of_element]] = 1
		group_dict[chunk_two[index_of_element]] = chunk_two_group[index_of_element] 			
	if(Blockdata["myid"] == Blockdata["main_node"]):
		chunk_ids = []
		group_ids = []
		partids   = read_text_file(os.path.join(Tracker["constants"]["masterdir"], "indexes.txt"),-1)
		partids   = partids[0]
		Tracker["constants"]["total_stack"] = len(partids)
		params = read_text_file(os.path.join(Tracker["constants"]["masterdir"], "sparx_refinement_params.txt"),-1)
		for index_of_particle in range(len(partids)): 
			chunk_ids.append(chunk_dict[partids[index_of_particle]])
			group_ids.append(group_dict[partids[index_of_particle]])
		refinement_params = [params[0], params[1], params[2], params[3], params[4], chunk_ids, group_ids, params[7]]
		write_text_file(refinement_params, os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt"))
	else: Tracker["constants"]["total_stack"] = 0
	Tracker["constants"]["total_stack"] = bcast_number_to_all(Tracker["constants"]["total_stack"], \
	    Blockdata["main_node"], MPI_COMM_WORLD)
	Tracker["total_stack"] = Tracker["constants"]["total_stack"]
	Tracker["constants"]["partstack"] = os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt")
	total_stack = Tracker["constants"]["total_stack"]
	Tracker["bckgnoise"] =  os.path.join(Tracker["constants"]["masterdir"], "bckgnoise.hdf")
	###
	from string import atoi
	if( Tracker["constants"]["minimum_grp_size"] ==-1 ):
		Tracker["constants"]["minimum_grp_size"] = Tracker["constants"]["total_stack"]\
		    //Tracker["constants"]["img_per_grp"]*(100//Blockdata["symclass"].nsym)
	if Tracker["constants"]["minimum_grp_size"] > Tracker["constants"]["img_per_grp"]:
		ERROR("User provided img_per_grp is smaller than minimum_grp_size", "get_input_from_sparx_ref3d", 1, Blockdata["myid"])
	# Now copy oldparamstruture
	copy_oldparamstructure_from_meridien_MPI(selected_iter, log_main)
	###
	if Tracker["constants"]["check_smearing"]:
		from statistics import scale_fsc_datasetsize
		number_of_groups = float(Tracker["constants"]["total_stack"])/Tracker["constants"]["img_per_grp"] 
		avg_fsc = scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], \
			 float(Tracker["constants"]["total_stack"]), \
			 float(Tracker["constants"]["total_stack"])//number_of_groups)
		fsc143 = len(avg_fsc)
		for ifreq in range(len(avg_fsc)):
			if avg_fsc[ifreq] < 0.143:
				fsc143 = ifreq -1
				break
		if fsc143 !=0: 
			nxinit = min((int(fsc143)+ max(int(Tracker["constants"]["nnxo"]*0.03), 5))*2, Tracker["constants"]["nnxo"])
		nxinit = smallprime(nxinit)
		nxinit = nxinit - nxinit%2 
		cdata_in_core =(Tracker["constants"]["total_stack"]*nxinit*nxinit*4.0)/1.e9/Blockdata["no_of_groups"]
		if not Tracker["constants"]["focus3D"]:fdata_in_core = 0.0
		else: fdata_in_core = cdata_in_core
		ctfdata = cdata_in_core
		refvol_size = (Tracker["nxinit"]*Tracker["nxinit"]*Tracker["nxinit"]*4.0*2)/1.e9 # including the 3D mask
		iterations = 2
		if(Blockdata["myid"] == Blockdata["main_node"]):
			while os.path.exists(os.path.join(Tracker["constants"]["refinement_dir"], "main%03d"%iterations)):
				iterations +=1
			log_main.add("---------->>> Smearing summary of %s  <<<----------"%Tracker["constants"]["refinement_dir"])
			log_main.add("Iter   average smearing    precalculated data/node     Percents of memory/node")
		else: iterations = 0
		iterations = bcast_number_to_all(iterations, Blockdata["main_node"], MPI_COMM_WORLD)
		if Tracker["constants"]["memory_per_node"] ==-1.:
			try:
				mem_bytes = os.sysconf('SC_PAGE_SIZE')*os.sysconf('SC_PHYS_PAGES')# e.g. 4015976448
				memory_per_node = mem_bytes/(1024.**3) # e.g. 3.74
			except: memory_per_node = 32.
		else: memory_per_node = Tracker["constants"]["memory_per_node"]
		smearing_all_iters = []
		for iteration in range(2, iterations):
			iter_smearings = get_smearing_info(Blockdata["nproc_previous"], iteration, Tracker["constants"]["total_stack"], \
			   Tracker["constants"]["masterdir"], Tracker["constants"]["refinement_dir"])
			savg = np.sum(iter_smearings)/Tracker["constants"]["total_stack"]
			srdata_in_core = (nxinit*nxinit*np.sum(iter_smearings)*4.)/1.e9/Blockdata["no_of_groups"]
			tdata = cdata_in_core+srdata_in_core+ctfdata+refvol_size
			percnt = tdata/memory_per_node*100.
			smearing_all_iters.append([iteration, savg, tdata, percnt])
			if(Blockdata["myid"] == Blockdata["main_node"]):
				log_main.add("%5d        %5.1f             %7.1f GB                  %7.1f%% "%(iteration, savg, tdata, percnt))
		if(Blockdata["myid"] == Blockdata["main_node"]):
			log_main.add('================================================================================================================')
			write_text_row(smearing_all_iters, os.path.join(Tracker["constants"]["masterdir"], "smearing_all.txt"))
	return import_from_sparx_refinement
		
def get_input_from_datastack(log_main):# Case three
	global Tracker, Blockdata
	import numpy as np
	from utilities import write_text_file, write_text_row, wrap_mpi_bcast
	import json
	from   string import split, atoi
	from   random import shuffle
	import_from_data_stack = 1
	import copy
	if(Blockdata["myid"] == Blockdata["main_node"]):
		image = get_im(Tracker["constants"]["orgstack"], 0)
		Tracker["constants"]["nnxo"] = image.get_xsize()		
		if( Tracker["nxinit"] > Tracker["constants"]["nnxo"]):
				ERROR("Image size less than minimum permitted %d"%Tracker["nxinit"], \
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
	    bcast_number_to_all(Tracker["constants"]["nnxo"], source_node = Blockdata["main_node"])
	if( Tracker["constants"]["nnxo"] < 0):
		ERROR("Image size is negative", "get_input_from_datastack", 1, Blockdata["main_node"])
	Tracker["constants"]["pixel_size"]	= \
	    bcast_number_to_all(Tracker["constants"]["pixel_size"], source_node =  Blockdata["main_node"])
	if(Tracker["constants"]["radius"] < 1):
		Tracker["constants"]["radius"]  = Tracker["constants"]["nnxo"]//2-2
	elif((2*Tracker["constants"]["radius"] +2) > Tracker["constants"]["nnxo"]):
		ERROR("Particle radius set too large!", "get_input_from_datastack",1, Blockdata["myid"])
	if Blockdata["myid"] == Blockdata["main_node"]:
		total_stack = EMUtil.get_image_count(Tracker["constants"]["orgstack"])
	else: total_stack = 0
	total_stack = bcast_number_to_all(total_stack, Blockdata["main_node"])
	# randomly assign two subsets
	Tracker["constants"]["total_stack"]	  = total_stack
	Tracker["constants"]["chunk_0"]		  = os.path.join(Tracker["constants"]["masterdir"],"chunk_0.txt")
	Tracker["constants"]["chunk_1"]		  = os.path.join(Tracker["constants"]["masterdir"],"chunk_1.txt")
	Tracker["constants"]["partstack"]	  = os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt")
	Tracker["previous_parstack"]          = os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt")#
	Tracker["constants"]["smearing_file"] = os.path.join(os.path.join(Tracker["constants"]["masterdir"], "all_smearing.txt"))
	if Blockdata["myid"] == Blockdata["main_node"]:
		write_text_file(np.full(Tracker["constants"]["total_stack"], 1, dtype = np.int16).tolist(),\
		   Tracker["constants"]["smearing_file"])
	
	if( Tracker["constants"]["minimum_grp_size"] == -1 ):
		Tracker["constants"]["minimum_grp_size"] = Tracker["constants"]["total_stack"]\
		     //Tracker["constants"]["img_per_grp"]*(100//Blockdata["symclass"].nsym)
	if( Tracker["constants"]["minimum_grp_size"] > Tracker["constants"]["img_per_grp"] ):
		ERROR("User provided img_per_grp is smaller than minimum_grp_size", \
		     "get_input_from_sparx_ref3d", 1, Blockdata["myid"])

	###
	Tracker["refang"], Tracker["rshifts"], Tracker["delta"] = None, None, None
	Tracker["avgnorm"] = 1.0
	chunk_dict         = {}
	chunk_list         = []
	if Blockdata["myid"] == Blockdata["main_node"]:
		chunk_dict  = {}
		tlist = list(range(total_stack))
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
		xform_proj_list = EMUtil.get_all_attributes(Tracker["constants"]["orgstack"], "xform.projection")
		for index_of_particle in range(len(xform_proj_list)):
			dp = xform_proj_list[index_of_particle].get_params("spider")
			xform_proj_list[index_of_particle] = \
			    [dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"], chunk_dict[index_of_particle]]
		write_text_row(xform_proj_list, Tracker["constants"]["partstack"])
	else:
		chunk_one = 0
		chunk_two = 0
	chunk_one = wrap_mpi_bcast(chunk_one, Blockdata["main_node"])
	chunk_two = wrap_mpi_bcast(chunk_two, Blockdata["main_node"])
	for element in chunk_one: chunk_dict[element] = 0
	for element in chunk_two: chunk_dict[element] = 1
	chunk_list 				= [chunk_one, chunk_two]	
	# Reconstruction to determine the resolution in orignal data size
	Tracker["nxinit"]     = Tracker["constants"]["nnxo"]
	Tracker["shrinkage"]  = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	Tracker["bckgnoise"]  =  None
	temp = model_blank(Tracker["constants"]["nnxo"], Tracker["constants"]["nnxo"])
	nny  =  temp.get_ysize()
	Blockdata["bckgnoise"] = [1.0]*nny # set for initial recon3D of data from stack	
	Tracker["fuse_freq"] = int(Tracker["constants"]["pixel_size"]*Tracker["constants"]["nnxo"]/Tracker["constants"]["fuse_freq"] +0.5)
	Tracker["directory"] = Tracker["constants"]["masterdir"]
	if Tracker["constants"]["nxinit"]< 0: Tracker["nxinit_refinement"] = Tracker["constants"]["nnxo"]
	else: Tracker["nxinit_refinement"] =  Tracker["constants"]["nxinit"]
	
	for procid in range(2):
		data = get_shrink_data_sorting(os.path.join(Tracker["constants"]["masterdir"],"chunk_%01d.txt"%procid), Tracker["constants"]["partstack"])
		mpi_barrier(MPI_COMM_WORLD)
		do3d_sorting(procid, data, myid = Blockdata["myid"],  mpi_comm = MPI_COMM_WORLD, log_main = log_main)# 1
	mpi_barrier(MPI_COMM_WORLD)
	
	if(Blockdata["no_of_groups"] == 1):
		if( Blockdata["myid"] == Blockdata["main_node"]):
			tvol0 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0.hdf"))
			tweight0 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_0.hdf"))
			tvol1 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_1.hdf"))
			tweight1 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_1.hdf"))
			Util.fuse_low_freq(tvol0, tvol1, tweight0, tweight1, 2*Tracker["fuse_freq"])
			shrank0 	= stepone(tvol0, tweight0)
			shrank1 	= stepone(tvol1, tweight1)
		    #  Note shrank volumes are Fourier uncentered.
			cfsc 		= fsc(shrank0, shrank1)[1]
			del shrank0, shrank1
			if(Tracker["nxinit"]<Tracker["constants"]["nnxo"]): cfsc  = cfsc[:Tracker["nxinit"]]
			for i in range(len(cfsc),Tracker["constants"]["nnxo"]//2+1): cfsc.append(0.0)
			write_text_row(cfsc, os.path.join(Tracker["directory"], "fsc_global.txt"))
			lcfsc  = len(cfsc)							
			Tracker["constants"]["fsc_curve"]  = copy.copy(cfsc)
			# scale 2. * x / (1 + x)
			for ifreq in range(len(Tracker["constants"]["fsc_curve"])):
				Tracker["constants"]["fsc_curve"][ifreq] = \
				  Tracker["constants"]["fsc_curve"][ifreq]*2./(1.+Tracker["constants"]["fsc_curve"][ifreq])
		else: Tracker = 0
		Tracker = wrap_mpi_bcast(Tracker,Blockdata["nodes"][0], communicator = MPI_COMM_WORLD)
	else:
		if(Blockdata["myid"] == Blockdata["nodes"][1]):  # It has to be 1 to avoid problem with tvol1 not closed on the disk
			tvol0 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_0.hdf"))
			tweight0 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_0.hdf"))
			tvol1 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_1.hdf"))
			tweight1 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_1.hdf"))
			Util.fuse_low_freq(tvol0, tvol1, tweight0, tweight1, 2*Tracker["fuse_freq"])
			tag = 7007
			send_EMData(tvol1, Blockdata["nodes"][0],    tag, MPI_COMM_WORLD)
			send_EMData(tweight1, Blockdata["nodes"][0], tag, MPI_COMM_WORLD)
			shrank0 	= stepone(tvol0, tweight0)
			send_EMData(shrank0, Blockdata["nodes"][0], tag, MPI_COMM_WORLD)
			del shrank0
			lcfsc = 0
		
		elif( Blockdata["myid"] == Blockdata["nodes"][0]):
			tag = 7007
			tvol1 		= recv_EMData(Blockdata["nodes"][1], tag, MPI_COMM_WORLD)
			tweight1 	= recv_EMData(Blockdata["nodes"][1], tag, MPI_COMM_WORLD)
			tvol1.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1} )
			shrank1 	= stepone(tvol1, tweight1)
			#  Get shrank volume, do fsc, send it to all
			shrank0 	= recv_EMData(Blockdata["nodes"][1], tag, MPI_COMM_WORLD)
			#  Note shrank volumes are Fourier uncentered.
			cfsc 		= fsc(shrank0, shrank1)[1]
			write_text_row(cfsc, os.path.join(Tracker["directory"], "fsc_global.txt"))
			del shrank0, shrank1
			if(Tracker["nxinit"]<Tracker["constants"]["nnxo"]):
				cfsc  = cfsc[:Tracker["nxinit"]]
				for i in range(len(cfsc),Tracker["constants"]["nnxo"]//2+1): cfsc.append(0.0)							
			Tracker["constants"]["fsc_curve"]  = copy.copy(cfsc)
			for ifreq in range(len(Tracker["constants"]["fsc_curve"])):
				Tracker["constants"]["fsc_curve"][ifreq] = \
				 Tracker["constants"]["fsc_curve"][ifreq]*2./(1.+Tracker["constants"]["fsc_curve"][ifreq])
		Tracker = wrap_mpi_bcast(Tracker, Blockdata["nodes"][0], communicator = MPI_COMM_WORLD)
	return import_from_data_stack
	
####	
### functions for fake rec3d from subsets
###
def do3d(procid, data, newparams, refang, rshifts, norm_per_particle, myid, mpi_comm = -1):
	global Tracker, Blockdata
	#  Without filtration
	from reconstruction import recons3d_trl_struct_MPI
	if( mpi_comm < -1 ): mpi_comm = MPI_COMM_WORDLD
	if Blockdata["subgroup_myid"]== Blockdata["main_node"]:
		if( procid == 0 ):
			if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
				os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
	shrinkage = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	tvol, tweight, trol = recons3d_trl_struct_MPI(myid = \
	    Blockdata["subgroup_myid"], main_node = Blockdata["main_node"], prjlist = data, \
		paramstructure = newparams, refang = refang, rshifts_shrank = [[q[0]*shrinkage,q[1]*shrinkage] for q in rshifts], \
		delta = Tracker["delta"], CTF = Tracker["constants"]["CTF"], upweighted = False, mpi_comm = mpi_comm, \
		target_size = (2*Tracker["nxinit"]+3), avgnorm = Tracker["avgvaradj"][procid], norm_per_particle = norm_per_particle)
	if Blockdata["subgroup_myid"] == Blockdata["main_node"]:
		tvol.set_attr("is_complex",0)
		tvol.write_image(os.path.join(Tracker["directory"], "tempdir", "tvol_%01d_%03d.hdf"%(procid,Tracker["mainiteration"])))
		tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%01d_%03d.hdf"%(procid,Tracker["mainiteration"])))
		trol.write_image(os.path.join(Tracker["directory"], "tempdir", "trol_%01d_%03d.hdf"%(procid,Tracker["mainiteration"])))
	mpi_barrier(mpi_comm)
	return
#######
def do3d_sorting_groups_rec3d(iteration, masterdir, log_main):
	global Tracker, Blockdata
	from utilities import get_im
	# reconstruct final unfiltered volumes from sorted clusters
	keepgoing = 1
	### ====
	Tracker["directory"]              = masterdir
	Tracker["constants"]["masterdir"] = masterdir
	Tracker["maxfrad"] = Tracker["nxinit"]//2
	####
	if Blockdata["no_of_groups"]>1:
		sub_main_node_list = [-1 for i in range(Blockdata["no_of_groups"])]
		for index_of_colors in range(Blockdata["no_of_groups"]):
			for iproc in range(Blockdata["nproc"]-1):
				if Blockdata["myid"]== iproc:
					if Blockdata["color"] == index_of_colors and Blockdata["myid_on_node"] == 0:
						sub_main_node_list[index_of_colors] = Blockdata["myid"]
					wrap_mpi_send(sub_main_node_list, Blockdata["last_node"], MPI_COMM_WORLD)
				if Blockdata["myid"] == Blockdata["last_node"]:
					dummy = wrap_mpi_recv(iproc, MPI_COMM_WORLD)
					for im in range(len(dummy)):
						if dummy[im]>-1: sub_main_node_list[im] = dummy[im]
				mpi_barrier(MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
		mpi_barrier(MPI_COMM_WORLD)
		wrap_mpi_bcast(sub_main_node_list, Blockdata["last_node"], MPI_COMM_WORLD)
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
					tvol2 		= get_im(os.path.join(Clusterdir, "tempdir", "tvol_0_%03d.hdf"%iteration))
					tweight2 	= get_im(os.path.join(Clusterdir, "tempdir", "tweight_0_%03d.hdf"%iteration))
					treg2 		= get_im(os.path.join(Clusterdir, "tempdir", "trol_0_%03d.hdf"%iteration))
					tag      = 7007
					send_EMData(tvol2, sub_main_node_list[index_of_colors],    tag,MPI_COMM_WORLD)
					send_EMData(tweight2, sub_main_node_list[index_of_colors], tag,MPI_COMM_WORLD)
					send_EMData(treg2, sub_main_node_list[index_of_colors],    tag,MPI_COMM_WORLD)
				elif (Blockdata["myid"] == sub_main_node_list[index_of_colors]):
					tag      = 7007
					tvol2    = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
					tweight2 = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
					treg2    = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
				mpi_barrier(MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
		
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				if Blockdata["color"] == index_of_colors:
					if( Blockdata["myid_on_node"] != 0):
						tvol2 		= model_blank(1)
						tweight2 	= model_blank(1)
						treg2		= model_blank(1)
					tvol2 = steptwo_mpi(tvol2, tweight2, treg2, None, False, color = index_of_colors) # has to be False!!!
					del tweight2, treg2
				mpi_barrier(Blockdata["shared_comm"])
			mpi_barrier(MPI_COMM_WORLD)
		
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				if (Blockdata["color"] == index_of_colors) and (Blockdata["myid_on_node"] == 0):
					tag = 7007
					send_EMData(tvol2, Blockdata["last_node"], tag, MPI_COMM_WORLD)
				elif(Blockdata["myid"] == Blockdata["last_node"]):
					tag = 7007
					tvol2 = recv_EMData(sub_main_node_list[index_of_colors], tag, MPI_COMM_WORLD)
					tvol2.write_image(os.path.join(Tracker["directory"], "vol_unfiltered_0_grp%03d.hdf"%index_of_group))
					del tvol2
				mpi_barrier(MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
		
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				Clusterdir = os.path.join(Tracker["directory"], "Cluster%d"%index_of_group, "main%03d"%iteration)
				if(Blockdata["myid"] == Blockdata["last_node"]):
					tvol2 		= get_im(os.path.join(Clusterdir, "tempdir", "tvol_1_%03d.hdf"%iteration))
					tweight2 	= get_im(os.path.join(Clusterdir, "tempdir", "tweight_1_%03d.hdf"%iteration))
					treg2 		= get_im(os.path.join(Clusterdir, "tempdir", "trol_1_%03d.hdf"%iteration))
					tag      = 7007
					send_EMData(tvol2,    sub_main_node_list[index_of_colors], tag,MPI_COMM_WORLD)
					send_EMData(tweight2, sub_main_node_list[index_of_colors], tag,MPI_COMM_WORLD)
					send_EMData(treg2,    sub_main_node_list[index_of_colors], tag,MPI_COMM_WORLD)
				
				elif (Blockdata["myid"] == sub_main_node_list[index_of_colors]):
					tag      = 7007
					tvol2       = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
					tweight2    = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
					treg2       = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
				mpi_barrier(MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
		
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				Tracker["maxfrad"] = Tracker["nxinit"]//2
				if Blockdata["color"] == index_of_colors:
					if( Blockdata["myid_on_node"] != 0):
						tvol2 		= model_blank(1)
						tweight2 	= model_blank(1)
						treg2		= model_blank(1)
					tvol2 = steptwo_mpi(tvol2, tweight2, treg2, None, False, color = index_of_colors) # has to be False!!!
					del tweight2, treg2
				mpi_barrier(Blockdata["shared_comm"])
			mpi_barrier(MPI_COMM_WORLD)
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				if (Blockdata["color"] == index_of_colors) and (Blockdata["myid_on_node"] == 0):
					tag = 7007
					send_EMData(tvol2, Blockdata["last_node"], tag, MPI_COMM_WORLD)
				elif(Blockdata["myid"] == Blockdata["last_node"]):
					tag = 7007
					tvol2 = recv_EMData(sub_main_node_list[index_of_colors], tag, MPI_COMM_WORLD)
					tvol2.write_image(os.path.join(Tracker["directory"], "vol_unfiltered_1_grp%03d.hdf"%index_of_group))
					del tvol2
				mpi_barrier(MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
	else:
		Tracker["maxfrad"] = Tracker["nxinit"]//2
		for index_of_group in range(Tracker["number_of_groups"]):
			Clusterdir = os.path.join(Tracker["directory"], "Cluster%d"%index_of_group, "main%03d"%iteration)
			
			if(Blockdata["myid"] == Blockdata["last_node"]):
				tvol2 		= get_im(os.path.join(Clusterdir, "tempdir", "tvol_0_%03d.hdf"%iteration))
				tweight2 	= get_im(os.path.join(Clusterdir, "tempdir", "tweight_0_%03d.hdf"%iteration))
				treg2 		= get_im(os.path.join(Clusterdir, "tempdir", "trol_0_%03d.hdf"%iteration))
				tag      = 7007
				send_EMData(tvol2, Blockdata["main_node"],    tag, MPI_COMM_WORLD)
				send_EMData(tweight2, Blockdata["main_node"], tag, MPI_COMM_WORLD)
				send_EMData(treg2, Blockdata["main_node"],    tag, MPI_COMM_WORLD)
			elif (Blockdata["myid"] == Blockdata["main_node"]):
				tag      = 7007
				tvol2    = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
				tweight2 = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
				treg2    = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			if Blockdata["myid"] != Blockdata["main_node"]:
				tvol2 		= model_blank(1)
				tweight2 	= model_blank(1)
				treg2		= model_blank(1)
			tvol2 = steptwo_mpi(tvol2, tweight2, treg2, None, False, color = 0) # has to be False!!!
			del tweight2, treg2
			if( Blockdata["myid"] == Blockdata["main_node"]):
				tvol2.write_image(os.path.join(Tracker["directory"], "vol_unfiltered_0_grp%03d.hdf"%index_of_group))
			mpi_barrier(MPI_COMM_WORLD)
			
			if(Blockdata["myid"] == Blockdata["last_node"]):
				tvol2 		= get_im(os.path.join(Clusterdir, "tempdir", "tvol_1_%03d.hdf"%iteration))
				tweight2 	= get_im(os.path.join(Clusterdir, "tempdir", "tweight_1_%03d.hdf"%iteration))
				treg2 		= get_im(os.path.join(Clusterdir, "tempdir", "trol_1_%03d.hdf"%iteration))
				tag      = 7007
				send_EMData(tvol2, Blockdata["main_node"],    tag, MPI_COMM_WORLD)
				send_EMData(tweight2, Blockdata["main_node"], tag, MPI_COMM_WORLD)
				send_EMData(treg2, Blockdata["main_node"],    tag, MPI_COMM_WORLD)
			elif (Blockdata["myid"] == Blockdata["main_node"]):
				tag      = 7007
				tvol2    = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
				tweight2 = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
				treg2    = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			if Blockdata["myid"] != Blockdata["main_node"]:
				tvol2 		= model_blank(1)
				tweight2 	= model_blank(1)
				treg2		= model_blank(1)
			tvol2 = steptwo_mpi(tvol2, tweight2, treg2, None, False, color = 0) # has to be False!!!
			del tweight2, treg2
			if( Blockdata["myid"] == Blockdata["main_node"]):
				tvol2.write_image(os.path.join(Tracker["directory"], "vol_unfiltered_1_grp%03d.hdf"%index_of_group))
			mpi_barrier(MPI_COMM_WORLD)
			
	keepgoing = bcast_number_to_all(keepgoing, source_node = Blockdata["main_node"], mpi_comm = MPI_COMM_WORLD) # always check 
	Tracker = wrap_mpi_bcast(Tracker, Blockdata["main_node"])
	if keepgoing == 0: 
		ERROR("do3d_sorting_groups_trl_iter  %s"%os.path.join(Tracker["directory"], \
		   "tempdir"),"do3d_sorting_groups_trl_iter", 1, Blockdata["myid"]) 
	return
####=====<--------
### nofsc rec3d
def do3d_sorting_groups_nofsc_smearing_iter(srdata, paramstructure, norm_per_particle,\
    partial_rec3d, current_group_sizes, iteration):
	global Tracker, Blockdata
	keepgoing = 1
	if(Blockdata["myid"] == Blockdata["last_node"]):
		if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
			os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
		try:
			fout = open(os.path.join(Tracker["directory"],"freq_cutoff.json"),'r')
			freq_cutoff_dict = convert_json_fromunicode(json.load(fout))
			fout.close()
		except: freq_cutoff_dict = 0
	else: freq_cutoff_dict = 0
	freq_cutoff_dict = wrap_mpi_bcast(freq_cutoff_dict, Blockdata["last_node"], MPI_COMM_WORLD)
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
			tvol.write_image(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf"%index_of_groups))
			tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf"%index_of_groups))
			trol.write_image(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_groups))
			del tvol
			del tweight
			del trol
	mpi_barrier(MPI_COMM_WORLD)
	Tracker["maxfrad"] = Tracker["nxinit"]//2
	fsc_groups = [None for im in range(len(current_group_sizes))]
	from statistics import scale_fsc_datasetsize
	total_size = sum(current_group_sizes)
	for igrp in range(len(current_group_sizes)):
		if Tracker["even_scale_fsc"]: 
			fsc_groups[igrp] = scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], \
		     float(Tracker["constants"]["total_stack"]), total_size//len(current_group_sizes))
		else:
			fsc_groups[igrp] = scale_fsc_datasetsize(Tracker["constants"]["fsc_curve"], \
		     float(Tracker["constants"]["total_stack"]), current_group_sizes[igrp])
	if Blockdata["no_of_groups"]>1:
	 	# new starts
		sub_main_node_list = [ -1 for i in range(Blockdata["no_of_groups"])]
		for index_of_colors in range(Blockdata["no_of_groups"]):
			for iproc in range(Blockdata["nproc"]-1):
				if Blockdata["myid"]== iproc:
					if Blockdata["color"] == index_of_colors and Blockdata["myid_on_node"] == 0:
						sub_main_node_list[index_of_colors] = Blockdata["myid"]
					wrap_mpi_send(sub_main_node_list, Blockdata["last_node"], MPI_COMM_WORLD)
				if Blockdata["myid"] == Blockdata["last_node"]:
					dummy = wrap_mpi_recv(iproc, MPI_COMM_WORLD)
					for im in range(len(dummy)):
						if dummy[im]>-1: sub_main_node_list[im] = dummy[im]
				mpi_barrier(MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
		mpi_barrier(MPI_COMM_WORLD)
		wrap_mpi_bcast(sub_main_node_list, Blockdata["last_node"], MPI_COMM_WORLD)
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
					tvol2    = get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")%index_of_group)
					tweight2 = get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")%index_of_group)
					treg2 	 = get_im(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_group))
					tvol2.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1})
					tag = 7007
					send_EMData(tvol2,    sub_main_node_list[index_of_colors], tag, MPI_COMM_WORLD)
					send_EMData(tweight2, sub_main_node_list[index_of_colors], tag, MPI_COMM_WORLD)
					send_EMData(treg2,    sub_main_node_list[index_of_colors], tag, MPI_COMM_WORLD)
				
				elif (Blockdata["myid"] == sub_main_node_list[index_of_colors]):
					tag      = 7007
					tvol2    = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
					tweight2 = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
					treg2    = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
				mpi_barrier(MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				try: 
					Tracker["freq_fsc143_cutoff"] = freq_cutoff_dict["Cluster_%03d.txt"%index_of_group]
				except: pass	
				if Blockdata["color"] == index_of_colors:  # It has to be 1 to avoid problem with tvol1 not closed on the disk 							
					if(Blockdata["myid_on_node"] != 0): 
						tvol2     = model_blank(1)
						tweight2  = model_blank(1)
						treg2     = model_blank(1)
					tvol2 = steptwo_mpi_reg(tvol2, tweight2, treg2,  fsc_groups[big_loop_groups[iloop][im]], \
					         Tracker["freq_fsc143_cutoff"], 0.01, True, color = index_of_colors) # has to be False!!!
					del tweight2, treg2
			mpi_barrier(MPI_COMM_WORLD)
		
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				if (Blockdata["color"] == index_of_colors) and (Blockdata["myid_on_node"] == 0):
					tag = 7007
					send_EMData(tvol2, Blockdata["last_node"], tag, MPI_COMM_WORLD)
				elif(Blockdata["myid"] == Blockdata["last_node"]):
					tag = 7007
					tvol2    = recv_EMData(sub_main_node_list[index_of_colors], tag, MPI_COMM_WORLD)
					tvol2.write_image(os.path.join(Tracker["directory"], "vol_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
					del tvol2
			mpi_barrier(MPI_COMM_WORLD)
		mpi_barrier(MPI_COMM_WORLD)
		
	else:# loop over all groups for single node
		for index_of_group in range(Tracker["number_of_groups"]):
			if(Blockdata["myid_on_node"] == 0):
				tvol2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")%index_of_group)
				tweight2 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")%index_of_group)
				treg2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_group))
				tvol2.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1})
			else:
				tvol2 		= model_blank(1)
				tweight2 	= model_blank(1)
				treg2		= model_blank(1)
			try: 
				Tracker["freq_fsc143_cutoff"] = freq_cutoff_dict["Cluster_%03d.txt"%index_of_group]
			except: pass
			
			tvol2 = steptwo_mpi_reg(tvol2, tweight2, treg2, fsc_groups[index_of_group],\
			      Tracker["freq_fsc143_cutoff"], 0.01, True) # has to be False!!!
			del tweight2, treg2
			if(Blockdata["myid_on_node"] == 0):
				tvol2.write_image(os.path.join(Tracker["directory"], "vol_grp%03d_iter%03d.hdf"%(index_of_group,iteration)))
				del tvol2
			mpi_barrier(MPI_COMM_WORLD)
	mpi_barrier(MPI_COMM_WORLD)
	keepgoing = bcast_number_to_all(keepgoing, source_node = Blockdata["main_node"], mpi_comm = MPI_COMM_WORLD) # always check 
	Tracker   = wrap_mpi_bcast(Tracker, Blockdata["main_node"])
	if not keepgoing: ERROR("do3d_sorting_groups_trl_iter  %s"%os.path.join(Tracker["directory"], "tempdir"),"do3d_sorting_groups_trl_iter", 1, Blockdata["myid"])
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
	if not os.path.exists(refvol_file):
		ERROR("refvol does not exist", "recons3d_trl_struct_group_nofsc_shifted_data_partial_MPI", 1, myid)
	if not os.path.exists(fftvol_file):
		ERROR("fftvol does not exist", "recons3d_trl_struct_group_nofsc_shifted_data_partial_MPI", 1, myid)
	if not os.path.exists(weight_file):
		ERROR("weight does not exist", "recons3d_trl_struct_group_nofsc_shifted_data_partial_MPI", 1, myid)
	
	#refvol
	if myid == main_node: target_size = get_im(refvol_file).get_xsize()
	else: target_size = 0
	target_size = bcast_number_to_all(target_size, main_node, mpi_comm)
	refvol = model_blank(target_size)# set to zero
	
	# fftvol
	if (myid == main_node): 
		fftvol = get_im(fftvol_file)
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
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, \
	  "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = Reconstructors.get("nn4_ctfws", params)
	r.setup()
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
						lshifts = findall(tdir[ii], ipsiandiang)
						toprab  = 0.0
						for ki in range(len(lshifts)):toprab += probs[lshifts[ki]]
						recdata = EMData(nny,nny,1,False)
						recdata.set_attr("is_complex",0)
						for ki in range(len(lshifts)):
							lpt = allshifts[lshifts[ki]]
							if(data[lpt] == None):
								data[lpt] = fshift(prjlist[im][0], rshifts_shrank[lpt][0], rshifts_shrank[lpt][1])
								data[lpt].set_attr("is_complex",0)
							Util.add_img(recdata, Util.mult_scalar(data[lpt], probs[lshifts[ki]]/toprab))
						recdata.set_attr_dict({"padffted":1, "is_fftpad":1,"is_fftodd":0, "is_complex_ri":1, "is_complex":1}) # preset already
						recdata = filt_table(recdata, bckgn)
						ipsi = tdir[ii]%100000
						iang = tdir[ii]/100000
						recdata.set_attr_dict({"bckgnoise":bckgn, "ctf":ct})
						if (current_group_ID  == group_ID):
							r.insert_slice(recdata, Transform({"type":"spider", "phi":refang[iang][0], "theta":refang[iang][1], \
								"psi":refang[iang][2]+ipsi*delta}), toprab*avgnorm/norm_per_particle[im])
						if (previous_group_ID == group_ID):
							r.insert_slice(recdata, Transform({"type":"spider", "phi":refang[iang][0], \
							  "theta":refang[iang][1], "psi":refang[iang][2]+ipsi*delta}), \
							     toprab*avgnorm/norm_per_particle[im]*(-1.))
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
def recons3d_trl_struct_group_nofsc_shifted_data_MPI(myid, main_node,\
       prjlist, paramstructure, norm_per_particle, group_ID,\
       mpi_comm = None, CTF = True, target_size=-1, nosmearing = False):
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
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1,\
	    "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = Reconstructors.get( "nn4_ctfw", params)
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
						lshifts = findall(tdir[ii], ipsiandiang)
						toprab  = 0.0
						for ki in range(len(lshifts)):toprab += probs[lshifts[ki]]
						recdata = EMData(nny,nny,1,False)
						recdata.set_attr("is_complex",0)
						for ki in range(len(lshifts)):
							lpt = allshifts[lshifts[ki]]
							if(data[lpt] == None):
								data[lpt] = fshift(prjlist[im][0], rshifts_shrank[lpt][0], rshifts_shrank[lpt][1])
								data[lpt].set_attr("is_complex",0)
							Util.add_img(recdata, Util.mult_scalar(data[lpt], probs[lshifts[ki]]/toprab))
						recdata.set_attr_dict({"padffted":1, "is_fftpad":1,"is_fftodd":0, "is_complex_ri":1, "is_complex":1}) # preset already
						recdata = filt_table(recdata, bckgn)
						ipsi = tdir[ii]%100000
						iang = tdir[ii]/100000
						recdata.set_attr_dict({"bckgnoise":bckgn, "ctf":ct})
						r.insert_slice(recdata, Transform({"type":"spider", "phi":refang[iang][0], "theta":refang[iang][1], \
						  "psi":refang[iang][2]+ipsi*delta}), toprab*avgnorm/norm_per_particle[im])		
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
	fftodd = 0
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, \
	  "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = Reconstructors.get( "nn4_ctfw", params )
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
						lshifts = findall(tdir[ii], ipsiandiang)
						toprab  = 0.0
						for ki in range(len(lshifts)):  toprab += probs[lshifts[ki]]
						recdata = EMData(nny,nny,1,False)
						recdata.set_attr("is_complex",0)
						for ki in range(len(lshifts)):
							lpt = allshifts[lshifts[ki]]
							if( data[lpt] == None ):
								data[lpt] = fshift(prjlist[im], rshifts_shrank[lpt][0], rshifts_shrank[lpt][1])
								data[lpt].set_attr("is_complex",0)
							Util.add_img(recdata, Util.mult_scalar(data[lpt], probs[lshifts[ki]]/toprab))
						recdata.set_attr_dict({"padffted":1, "is_fftpad":1,"is_fftodd": fftodd, "is_complex_ri":1, "is_complex":1})
						if not upweighted:  recdata = filt_table(recdata, bckgn )
						recdata.set_attr_dict( {"bckgnoise":bckgn, "ctf":ct} )
						ipsi = tdir[ii]%100000
						iang = tdir[ii]/100000
						r.insert_slice( recdata, \
						   Transform({"type":"spider","phi":refang[iang][0],"theta":refang[iang][1],"psi":refang[iang][2]+ipsi*delta}), \
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
							lshifts = findall(tdir[ii], ipsiandiang)
							toprab  = 0.0
							for ki in range(len(lshifts)):  toprab += probs[lshifts[ki]]
							recdata = EMData(nny,nny,1,False)
							recdata.set_attr("is_complex",0)
							for ki in range(len(lshifts)):
								lpt = allshifts[lshifts[ki]]
								if( data[lpt] == None ):
									data[lpt] = fshift(prjlist[im], rshifts_shrank[lpt][0], rshifts_shrank[lpt][1])
									data[lpt].set_attr("is_complex",0)
								Util.add_img(recdata, Util.mult_scalar(data[lpt], probs[lshifts[ki]]/toprab))
							recdata.set_attr_dict({"padffted":1, "is_fftpad":1,"is_fftodd":fftodd, "is_complex_ri":1, "is_complex":1})
							if not upweighted:  recdata = filt_table(recdata, bckgn )
							recdata.set_attr_dict( {"bckgnoise":bckgn, "ctf":ct} )
							ipsi = tdir[ii]%100000
							iang = tdir[ii]/100000
							r.insert_slice(recdata, \
							 Transform({"type":"spider","phi":refang[iang][0],"theta":refang[iang][1],"psi":refang[iang][2]+ipsi*delta}), \
							  toprab*avgnorm/norm_per_particle[im])
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
####==========--------
#####  FSC rec3d
def do3d_sorting_groups_fsc_only_iter(data, paramstructure, norm_per_particle, iteration):
	global Tracker, Blockdata
	# do resolution each time
	keepgoing = 1
	if(Blockdata["myid"] == Blockdata["nodes"][0]):
		if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
			os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
	do3d_sorting_group_insertion_random_two_for_fsc(data, paramstructure, norm_per_particle)
	for index_of_colors in range(Blockdata["no_of_groups"]):
		group_start, group_end = MPI_volume_start_end(Tracker["number_of_groups"], Blockdata["no_of_groups"], index_of_colors)
		if Blockdata["color"] == index_of_colors:  # It has to be 1 to avoid problem with tvol1 not closed on the disk
			for index_of_group in range(group_start, group_end):								
				if Blockdata["myid_on_node"] == 0:
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
					send_EMData(tvol1,    Blockdata["no_of_processes_per_group"]-2, tag, Blockdata["shared_comm"])
					send_EMData(tweight1, Blockdata["no_of_processes_per_group"]-2, tag, Blockdata["shared_comm"])
					shrank0 	= stepone(tvol0, tweight0)					
				elif Blockdata["myid_on_node"] == Blockdata["no_of_processes_per_group"]-2:
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
					tag = 7007
					shrank0 	= recv_EMData(0, tag, Blockdata["shared_comm"])
					cfsc 		= fsc(shrank0, shrank1)[1]
					write_text_row(cfsc, os.path.join(Tracker["directory"], "fsc_driver_chunk0_grp%03d_iter%03d.txt")%(index_of_group,iteration))
					del shrank0, shrank1
					if(Tracker["nxinit"]<Tracker["constants"]["nnxo"]):
						cfsc 	= cfsc[:Tracker["nxinit"]]
						for i in range(len(cfsc),Tracker["constants"]["nnxo"]//2+1):  cfsc.append(0.0)
					lcfsc = len(cfsc)
				if 	Blockdata["myid_on_node"] == 1:
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
						for i in range(len(cfsc),Tracker["constants"]["nnxo"]//2+1):cfsc.append(0.0)
					lcfsc = len(cfsc)							
				Tracker = wrap_mpi_bcast(Tracker, Blockdata["no_of_processes_per_group"]-1, Blockdata["shared_comm"])			
				mpi_barrier(Blockdata["shared_comm"])
			mpi_barrier(Blockdata["shared_comm"])
	mpi_barrier(MPI_COMM_WORLD)
	keepgoing = bcast_number_to_all(keepgoing, source_node = Blockdata["main_node"], mpi_comm = MPI_COMM_WORLD) # always check 
	Tracker   = wrap_mpi_bcast(Tracker, Blockdata["main_node"])
	if not keepgoing:ERROR("do3d_sorting_groups_trl_iter  %s"%os.path.join(Tracker["directory"], "tempdir"),"do3d_sorting_groups_trl_iter", 1, Blockdata["myid"]) 
	return
	
def do3d_sorting_group_insertion_random_two_for_fsc(data, sparamstructure, snorm_per_particle):
	global Tracker, Blockdata
	if(Blockdata["myid"] == Blockdata["nodes"][0]):
		if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
			os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
	for index_of_groups in range(Tracker["number_of_groups"]):
		for procid in range(2):
			for ifsc in range(2):
				tvol, tweight, trol = recons3d_4nnsorting_group_fsc_MPI(myid = Blockdata["myid"], main_node = Blockdata["nodes"][procid],  \
				  prjlist = data, fsc_half = ifsc, random_subset = procid, group_ID = index_of_groups, \
				   paramstructure=sparamstructure, norm_per_particle=snorm_per_particle,\
				  CTF = Tracker["constants"]["CTF"], upweighted = False, target_size = (2*Tracker["nxinit"]+3))
				if(Blockdata["myid"] == Blockdata["nodes"][procid]):
					tvol.set_attr("is_complex",0)
					tvol.write_image(os.path.join(   Tracker["directory"], "tempdir", "tvol_%d_%d_%d.hdf"%(ifsc,    procid, index_of_groups)))
					tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%d_%d_%d.hdf"%(ifsc, procid, index_of_groups)))
					trol.write_image(os.path.join(   Tracker["directory"], "tempdir", "trol_%d_%d_%d.hdf"%(ifsc,    procid, index_of_groups)))
				mpi_barrier(MPI_COMM_WORLD)
		mpi_barrier(MPI_COMM_WORLD)
	mpi_barrier(MPI_COMM_WORLD)
	return
#### Never read volumes
		
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
	refvol  = model_blank(target_size)
	refvol.set_attr("fudge", 1.0)
	if CTF: do_ctf = 1
	else:   do_ctf = 0
	fftvol = EMData()
	weight = EMData()
	try:    qt = prjlist[0].get_attr("qt")
	except: qt = 1.0
	params = {"size":target_size, "npad":2, "snr":1.0, "sign":1, \
	   "symmetry":"c1", "refvol":refvol, "fftvol":fftvol, "weight":weight, "do_ctf": do_ctf}
	r = Reconstructors.get( "nn4_ctfw", params )
	r.setup()
	if norm_per_particle == None: norm_per_particle = [1.0 for im in range(len(prjlist))]
	# Definitions for smearing ,all copied from refinement
	if not Tracker["nosmearing"]:
		delta           = Tracker["delta"]
		refang          = Tracker["refang"]
		rshifts_shrank  = copy.deepcopy(Tracker["rshifts"])
		for im in range(len(rshifts_shrank)):
			rshifts_shrank[im][0] *= float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
			rshifts_shrank[im][1] *= float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	nnx = prjlist[0].get_xsize()
	nny = prjlist[0].get_ysize()
	nc  = 0
	for im in range(len(prjlist)):
		if prjlist[im].get_attr("group") == group_ID and prjlist[im].get_attr("chunk_id") == random_subset:
			if Tracker["nosmearing"]: avgnorm = 1.0
			else: avgnorm = Tracker["avgnorm"][prjlist[im].get_attr("chunk_id")]#
			if nc %2 == fsc_half:
				if Tracker["nosmearing"]:
					if not upweighted: prjlist[im] = filt_table(prjlist[im], prjlist[im].get_attr("bckgnoise"))
					r.insert_slice(prjlist[im], prjlist[im].get_attr("xform.projection"), 1.0)
				else:
					if Tracker["constants"]["nsmear"]<=0.0: numbor = len(paramstructure[im][2])
					else:  numbor =1
					ipsiandiang = [paramstructure[im][2][i][0]/1000  for i in range(numbor)]
					allshifts   = [paramstructure[im][2][i][0]%1000  for i in range(numbor)]
					probs       = [paramstructure[im][2][i][1] for i in range(numbor)]
					#  Find unique projection directions
					tdir  = list(set(ipsiandiang))
					bckgn = prjlist[im].get_attr("bckgnoise")
					ct    = prjlist[im].get_attr("ctf")
					#  For each unique projection direction:
					data = [None]*len(rshifts_shrank)
					for ii in range(len(tdir)):
						#  Find the number of times given projection direction appears on the list, it is the number of different shifts associated with it.
						lshifts = findall(tdir[ii], ipsiandiang)
						toprab  = 0.0
						for ki in range(len(lshifts)):  toprab += probs[lshifts[ki]]
						recdata = EMData(nny,nny,1,False)
						recdata.set_attr("is_complex",0)
						for ki in range(len(lshifts)):
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
						r.insert_slice( recdata, \
						 Transform({"type":"spider","phi":refang[iang][0],"theta":refang[iang][1],"psi":refang[iang][2]+ipsi*delta}), \
						    toprab*avgnorm/norm_per_particle[im])
			nc +=1
	reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
	reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
	if myid == main_node: dummy = r.finish(True)
	mpi_barrier(mpi_comm)
	if myid == main_node: return fftvol, weight, refvol
	else: return None, None, None
#####end of FSC
###=====<-----group rec3d
### insertion
def do3d_sorting_group_insertion_smearing(sdata, sparamstructure, snorm_per_particle, randomset=2):
	global Tracker, Blockdata
	if(Blockdata["myid"] == Blockdata["nodes"][0]):
		if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
			os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
	if randomset ==1:
		for index_of_groups in range(Tracker["number_of_groups"]):
			for procid in range(2, 3):
				tvol, tweight, trol = recons3d_trl_struct_group_MPI(myid = \
				  Blockdata["myid"], main_node = Blockdata["nodes"][procid//2],prjlist = sdata, \
				  random_subset = procid, group_ID = index_of_groups, paramstructure = sparamstructure, \
				  norm_per_particle = snorm_per_particle, CTF = Tracker["constants"]["CTF"],\
				  mpi_comm = None, upweighted = False, target_size = (2*Tracker["nxinit"]+3), \
				  nosmearing = Tracker["nosmearing"])
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
		for index_of_groups in range(Tracker["number_of_groups"]):
			for procid in range(3):
				tvol, tweight, trol = recons3d_trl_struct_group_MPI(myid = \
				    Blockdata["myid"], main_node = Blockdata["nodes"][procid//2],\
				    prjlist = sdata,  random_subset = procid, group_ID = index_of_groups, \
				    paramstructure = sparamstructure, norm_per_particle = snorm_per_particle, \
				    CTF = Tracker["constants"]["CTF"],mpi_comm= None, upweighted = False, \
				    target_size = (2*Tracker["nxinit"]+3), nosmearing = Tracker["nosmearing"])
				if(Blockdata["myid"] == Blockdata["nodes"][procid//2]):
					tvol.set_attr("is_complex",0)
					tvol.write_image(os.path.join(   Tracker["directory"], "tempdir", "tvol_%d_%d.hdf"%(procid,    index_of_groups)))
					tweight.write_image(os.path.join(Tracker["directory"], "tempdir", "tweight_%d_%d.hdf"%(procid, index_of_groups)))
					trol.write_image(os.path.join(   Tracker["directory"], "tempdir", "trol_%d_%d.hdf"%(procid,    index_of_groups)))
				del tvol
				del tweight
				del trol
				mpi_barrier(MPI_COMM_WORLD)
	mpi_barrier(MPI_COMM_WORLD)
	return
### rec3d 
####=====<-------MEM related functions
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
    
#####==========-------------------------Functions for post processing
def compute_final_map(work_dir, log_main):
	global Tracker, Blockdata
	Tracker["constants"]["orgres"]			 = 0.0
	Tracker["constants"]["refinement_delta"] = 0.0
	Tracker["constants"]["refinement_ts"]	 = 0.0
	Tracker["constants"]["refinement_xr"]	 = 0.0
	Tracker["constants"]["refinement_an"]	 = 0.0
	minimum_size     = Tracker["constants"]["img_per_grp"]
	number_of_groups = 0
	
	if(Blockdata["myid"] == Blockdata["main_node"]):
		final_accounted_ptl = 0
		fout = open(os.path.join(work_dir, "Tracker.json"),'w')
		json.dump(Tracker, fout)
		fout.close()
		clusters = []
		while os.path.exists(os.path.join(work_dir, "Cluster_%03d.txt"%number_of_groups)):
			class_in          = read_text_file(os.path.join(work_dir, "Cluster_%03d.txt"%number_of_groups))
			minimum_size      = min(len(class_in), minimum_size)
			number_of_groups += 1
			final_accounted_ptl +=len(class_in)
			clusters.append(class_in)
			del class_in
		try: # check the unaccounted for particles
			class_in = read_text_file(os.path.join(work_dir, "Core_set.txt"))
			if Tracker["constants"]["num_core_set"] ==-1:
				if len(class_in)> max(100, Tracker["constants"]["minimum_grp_size"]): 
					Tracker["Core_set"]   = True
					clusters.append(class_in)
					minimum_size      = min(len(class_in), minimum_size)
					number_of_groups += 1
					final_accounted_ptl +=len(class_in)
				else: Tracker["Core_set"] = False
				
			elif len(class_in) >Tracker["constants"]["num_core_set"]:
				Tracker["Core_set"] = True
				clusters.append(class_in)
				minimum_size = min(len(class_in), minimum_size)
				number_of_groups += 1
				final_accounted_ptl +=len(class_in)
				
			else: Tracker["Core_set"] = False
			del class_in
		except:
			log_main.add(' Core_set.txt does not exist   ')
			Tracker["Core_set"] = False
		Tracker["total_stack"]      = final_accounted_ptl
		Tracker["number_of_groups"] = number_of_groups
		Tracker["nxinit"]           = Tracker["nxinit_refinement"]
	else: Tracker = 0
	number_of_groups = bcast_number_to_all(number_of_groups, Blockdata["main_node"], MPI_COMM_WORLD)
	Tracker          = wrap_mpi_bcast(Tracker,               Blockdata["main_node"], MPI_COMM_WORLD)
	if(number_of_groups == 0): ERROR("No clusters  found, the program terminates.", "compute_final_map", 1, Blockdata["myid"])
	compute_noise( Tracker["nxinit"])
	if(Blockdata["myid"] == Blockdata["main_node"]):
		alist, partition = merge_classes_into_partition_list(clusters)
		write_text_row(partition, os.path.join(work_dir, "generation_partition.txt"))
	parti_file      = os.path.join(work_dir, "generation_partition.txt")
	params          = os.path.join(Tracker["constants"]["masterdir"],"refinement_parameters.txt")
	previous_params = Tracker["previous_parstack"]
	original_data, norm_per_particle  = read_data_for_sorting(parti_file, params, previous_params)
	
	if Tracker["nosmearing"]: parameterstructure  = None
	else:
		parameterstructure  = read_paramstructure_for_sorting(parti_file,  Tracker["paramstructure_dict"], Tracker["paramstructure_dir"])
	mpi_barrier(MPI_COMM_WORLD)
	
	Tracker["directory"] = work_dir
	compute_noise(Tracker["nxinit"])
	rdata = downsize_data_for_rec3D(original_data, Tracker["nxinit"], False, 1)# pay attentions to shifts!
	for im in range(len(original_data)): rdata[im].set_attr('group', original_data[im].get_attr('group'))
	del original_data
	mpi_barrier(MPI_COMM_WORLD)
	do3d_sorting_groups_nofsc_final(rdata, parameterstructure, norm_per_particle)
	del rdata
	import  shutil
	if Blockdata["myid"] == Blockdata["main_node"]:
		if os.path.exists(os.path.join(Tracker["constants"]["masterdir"], 'tempdir')):
			shutil.rmtree(os.path.join(Tracker["constants"]["masterdir"], 'tempdir'))
	mpi_barrier(MPI_COMM_WORLD)
	return
### final rec3d
def do3d_sorting_groups_nofsc_final(rdata, parameterstructure, norm_per_particle):
	global Tracker, Blockdata
	keepgoing = 1
	if(Blockdata["myid"] == Blockdata["last_node"]):
		if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
			os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
		try:
			fout = open(os.path.join(Tracker["directory"],"freq_cutoff.json"),'r')
			freq_cutoff_dict = convert_json_fromunicode(json.load(fout))
			fout.close()
		except: freq_cutoff_dict = 0
	else: freq_cutoff_dict = 0
	freq_cutoff_dict = wrap_mpi_bcast(freq_cutoff_dict, Blockdata["last_node"], MPI_COMM_WORLD)
		
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
	mpi_barrier(MPI_COMM_WORLD)
	Tracker["maxfrad"] = Tracker["nxinit"]//2
	if Blockdata["no_of_groups"]>1:
	 	# new starts
		sub_main_node_list = [ -1 for i in range(Blockdata["no_of_groups"])]
		for index_of_colors in range(Blockdata["no_of_groups"]):
			for iproc in range(Blockdata["nproc"]-1):
				if Blockdata["myid"]== iproc:
					if Blockdata["color"] == index_of_colors and Blockdata["myid_on_node"] == 0:
						sub_main_node_list[index_of_colors] = Blockdata["myid"]
					wrap_mpi_send(sub_main_node_list, Blockdata["last_node"], MPI_COMM_WORLD)
				if Blockdata["myid"] == Blockdata["last_node"]:
					dummy = wrap_mpi_recv(iproc, MPI_COMM_WORLD)
					for im in range(len(dummy)):
						if dummy[im]>-1: sub_main_node_list[im] = dummy[im]
				mpi_barrier(MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
		mpi_barrier(MPI_COMM_WORLD)
		wrap_mpi_bcast(sub_main_node_list, Blockdata["last_node"], MPI_COMM_WORLD)
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
					tvol2    = get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")%index_of_group)
					tweight2 = get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")%index_of_group)
					treg2 	 = get_im(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_group))
					tvol2.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1})
					tag = 7007
					send_EMData(tvol2,    sub_main_node_list[index_of_colors], tag, MPI_COMM_WORLD)
					send_EMData(tweight2, sub_main_node_list[index_of_colors], tag, MPI_COMM_WORLD)
					send_EMData(treg2,    sub_main_node_list[index_of_colors], tag, MPI_COMM_WORLD)
				
				elif (Blockdata["myid"] == sub_main_node_list[index_of_colors]):
					tag      = 7007
					tvol2    = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
					tweight2 = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
					treg2    = recv_EMData(Blockdata["last_node"], tag, MPI_COMM_WORLD)
				mpi_barrier(MPI_COMM_WORLD)
			mpi_barrier(MPI_COMM_WORLD)
			
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				try: 
					Tracker["freq_fsc143_cutoff"] = freq_cutoff_dict["Cluster_%03d.txt"%index_of_group]
				except: pass	
				if Blockdata["color"] == index_of_colors:  # It has to be 1 to avoid problem with tvol1 not closed on the disk 							
					if(Blockdata["myid_on_node"] != 0): 
						tvol2     = model_blank(1)
						tweight2  = model_blank(1)
						treg2     = model_blank(1)
					tvol2 = steptwo_mpi_reg(tvol2, tweight2, treg2, None, \
					 Tracker["freq_fsc143_cutoff"], 0.01, False, color = index_of_colors) # has to be False!!!
					del tweight2, treg2
			mpi_barrier(MPI_COMM_WORLD)
			for im in range(len(big_loop_colors[iloop])):
				index_of_group  = big_loop_groups[iloop][im]
				index_of_colors = big_loop_colors[iloop][im]
				if (Blockdata["color"] == index_of_colors) and (Blockdata["myid_on_node"] == 0):
					tag = 7007
					send_EMData(tvol2, Blockdata["last_node"], tag, MPI_COMM_WORLD)
				elif(Blockdata["myid"] == Blockdata["last_node"]):
					tag = 7007
					tvol2    = recv_EMData(sub_main_node_list[index_of_colors], tag, MPI_COMM_WORLD)
					if Tracker["Core_set"]:
						if index_of_group !=Tracker["number_of_groups"]-1:
							tvol2.write_image(os.path.join(Tracker["directory"], "vol_cluster%03d.hdf"%index_of_group))
						else:
							tvol2.write_image(os.path.join(Tracker["directory"], "volume_core.hdf"))
					else: tvol2.write_image(os.path.join(Tracker["directory"], "vol_cluster%03d.hdf"%index_of_group))
					del tvol2
			mpi_barrier(MPI_COMM_WORLD)
		mpi_barrier(MPI_COMM_WORLD)
		
	else:# loop over all groups for single node
		for index_of_group in range(Tracker["number_of_groups"]):
			if(Blockdata["myid_on_node"] == 0):
				tvol2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")%index_of_group)
				tweight2 	= get_im(os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")%index_of_group)
				treg2 		= get_im(os.path.join(Tracker["directory"], "tempdir", "trol_2_%d.hdf"%index_of_group))
				tvol2.set_attr_dict( {"is_complex":1, "is_fftodd":1, 'is_complex_ri': 1, 'is_fftpad': 1})
			else:
				tvol2 		= model_blank(1)
				tweight2 	= model_blank(1)
				treg2		= model_blank(1)
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
			mpi_barrier(MPI_COMM_WORLD)
	mpi_barrier(MPI_COMM_WORLD)
	keepgoing = bcast_number_to_all(keepgoing, source_node = Blockdata["main_node"], mpi_comm = MPI_COMM_WORLD) # always check 
	Tracker   = wrap_mpi_bcast(Tracker, Blockdata["main_node"])
	if not keepgoing: ERROR("do3d_sorting_groups_nofsc_final  %s"%os.path.join(Tracker["directory"], \
	     "tempdir"),"do3d_sorting_groups_nofsc_final", 1, Blockdata["myid"])
	return
#####==========----various utilities
def get_time(time_start):
	current_time   = time() - time_start
	current_time_h = current_time // 3600
	current_time_m = (current_time - current_time_h*3600)// 60
	return int(current_time_h), int(current_time_m)
	
def check_sorting(total_data, keepsorting, log_file):
	global Tracker, Blockdata
	import json
	if Blockdata["myid"] == Blockdata["main_node"]:
		fout         = open(os.path.join(Tracker["constants"]["masterdir"],"Tracker.json"),'r')
		Tracker_main = convert_json_fromunicode(json.load(fout))
		fout.close()
	else: Tracker_main = 0
	Tracker_main = wrap_mpi_bcast(Tracker_main, Blockdata["main_node"])
	n            = 0
	keepsorting  = 0
	while n <= Tracker_main["constants"]["img_per_grp_split_rate"]:
		if Tracker_main["constants"]["img_per_grp"]//2**n> Tracker_main["constants"]["minimum_grp_size"]:
			Tracker["number_of_groups"] = (2**n*total_data)//Tracker_main["constants"]["img_per_grp"]
			if Tracker["number_of_groups"] >=2:
				Tracker["current_img_per_grp"] = Tracker_main["constants"]["img_per_grp"]//2**n
				keepsorting = 1
		else:
			Tracker["number_of_groups"] = total_data//Tracker_main["constants"]["minimum_grp_size"]
			if Tracker["number_of_groups"]>=2:
				Tracker["current_img_per_grp"] = Tracker_main["constants"]["minimum_grp_size"]
				keepsorting = 1
			else: keepsorting = 0
			
		if Tracker["number_of_groups"] >=2:
			keepsorting = 1
			break
		else: n+=1
	if keepsorting ==1:
		Tracker["total_stack"] = total_data
		keepsorting            = sort3d_init("initialization", log_file)
	return keepsorting

def copy_results(log_file, all_gen_stat_list):
	global Tracker, Blockdata
	import json
	from   shutil import copyfile
	from   string import atoi
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
				ig = atoi('%s'%ig1)
				for ic in range(value):
					cluster_file = os.path.join(Tracker["constants"]["masterdir"], "generation_%03d"%ig, "Cluster_%03d.txt"%ic)
					copyfile(cluster_file, os.path.join(Tracker["constants"]["masterdir"], "Cluster_%03d.txt"%nclusters))
					clusters.append(read_text_file(cluster_file))
					cluster      = read_text_file(os.path.join(Tracker["constants"]["masterdir"], "generation_%03d"%ig, "Cluster_%03d.txt"%ic))
					cluster_file = "Cluster_%03d.txt"%nclusters
					vol_file     = "vol_cluster%03d.hdf"%nclusters
					msg          = '{:>8} {:>8}   {:^24}        {:^6}          {:^6}          {:>5} {:^20} {:^20} '.format(nclusters, \
					  len(cluster), ig, round(all_gen_stat_list[ig][ic][0],1), round(all_gen_stat_list[ig][ic][1],1), \
					    round(all_gen_stat_list[ig][ic][2],1), cluster_file,  vol_file)
					nclusters   +=1
					NACC        +=len(cluster)
					log_file.add(msg)			
			try:
				Unaccounted_file = os.path.join(Tracker["constants"]["masterdir"], "generation_%03d"%Tracker["current_generation"], "Core_set.txt")
				cluster = read_text_file(Unaccounted_file)
				copyfile(Unaccounted_file, os.path.join(Tracker["constants"]["masterdir"], "Core_set.txt"))
				cluster_file = "Core_set.txt"
				if (Tracker["constants"]["num_core_set"]>0) and (len(cluster) > Tracker["constants"]["num_core_set"]): 
					vol_file = "volume_core.hdf"
				elif (Tracker["constants"]["num_core_set"]==-1) and (len(cluster) > max(Tracker["constants"]["minimum_grp_size"], 100)): 
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
			if    len(clusters) >=2: do_analysis_on_identified_clusters(clusters, log_file)
			else: log_file.add(' ANOVA analysis is skipped ')
			fout = open(os.path.join(Tracker["constants"]["masterdir"], "Tracker.json"), 'w')
			json.dump(Tracker, fout)
			fout.close()
		mpi_barrier(MPI_COMM_WORLD)
	else:
		if Blockdata["myid"] == Blockdata["main_node"]: log_file.add('No groups are found.\n')
	return

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
			else:
				table_k_k[i][j] = 0
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

def estimate_tanhl_params(cutoff, taa, image_size):
	from math import sqrt
	def tanhfl(x, cutoff, taa):
		from math import pi, tanh
		omega = cutoff
		cnst  = pi/(2.0*omega*taa)
		v1    = (cnst*(x + omega))
		v2    = (cnst*(x - omega))
		return 0.5*(tanh(v1) - tanh(v2))
	
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

def print_matching_pairs(pair_list, log_file):
	log_file.add(' ')
	log_file.add('                        Two-way matching of sorting results.')
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
	
def do_random_groups_simulation_mpi(ptp1, ptp2):
	global Tracker, Blockdata
	import numpy as np
	# return two lists: group avgs and group stds. The last one of two lists are the total avg and std.
	if (len(ptp1)>=50) or (len(ptp2)>=50):
		if(Blockdata["myid"] == Blockdata["main_node"]): print('Warning: too many simulation groups')
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
			newindeces, list_stable, nb_tot_objs = k_means_match_clusters_asg_new(new_clusters1,new_clusters2)
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
		svar +=table_stat(tlist)[0]
		save +=table_stat(tlist)[1]
		for j in range(k):
			try: gvar[j] +=table_stat(clist[j])[1]
			except: gvar[j] +=0.0
			try: gave[j] +=table_stat(clist[j])[0]
			except: gave[j] +=0.0
			
	save = mpi_reduce(save, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD)
	svar = mpi_reduce(svar, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD)
	
	for iproc in range(1, Blockdata["nproc"]):
		if Blockdata["myid"] == iproc:
			wrap_mpi_send(gvar,  Blockdata["main_node"], MPI_COMM_WORLD)
			
		elif Blockdata["myid"] == Blockdata["main_node"]:
			dummy = wrap_mpi_recv(iproc, MPI_COMM_WORLD)
			for im  in  range(len(dummy)): gvar[im]+=dummy[im]
			
		mpi_barrier(MPI_COMM_WORLD)

	for iproc in range(1, Blockdata["nproc"]):
		if Blockdata["myid"] == iproc:
			wrap_mpi_send(gave,  Blockdata["main_node"], MPI_COMM_WORLD)
			
		elif Blockdata["myid"] == Blockdata["main_node"]:
			dummy = wrap_mpi_recv(iproc, MPI_COMM_WORLD)
			for im  in  range(len(dummy)): gave[im]+=dummy[im]
		mpi_barrier(MPI_COMM_WORLD)

	from math import sqrt
	if(Blockdata["myid"] == Blockdata["main_node"]):
		svar = sqrt(svar/float(Blockdata["nproc"]*Nloop))
		save = save/float(Blockdata["nproc"]*Nloop)
		for i in range(k):
			gave[i] = gave[i]/float(Blockdata["nproc"]*Nloop)
			gvar[i] = sqrt(gvar[i]/float(Blockdata["nproc"]*Nloop))
		gave.append(save)
		gvar.append(svar)
	else:
		gave = 0
		gvar = 0
	gave = wrap_mpi_bcast(gave, Blockdata["main_node"], MPI_COMM_WORLD)
	gvar = wrap_mpi_bcast(gvar, Blockdata["main_node"], MPI_COMM_WORLD)
	return gave, gvar
	
def get_group_size_from_iter_assign(iter_assignment):
	import numpy as np
	group_sizes  = []
	my_assign    = np.array(iter_assignment, dtype=np.int32)
	group_ids    = np.sort(np.unique(my_assign))
	for im in range(group_ids.shape[0]):
		group_sizes.append((np.sort(my_assign[isin(my_assign, group_ids[im])])).shape[0])
	return group_sizes

def get_res143(fsc_curve):
	fsc143 = len(fsc_curve)
	for ifreq in range(len(fsc_curve)):
		if (fsc_curve[ifreq] < 0.143):
			fsc143 = ifreq -1
			break
	return fsc143
		
def sorting_main_mpi(log_main, depth_order, not_include_unaccounted):
	global Tracker, Blockdata
	time_sorting_start = time()
	read_tracker_mpi(Tracker["constants"]["masterdir"])
	Tracker["generation"]         =  {}
	Tracker["current_generation"] =  0
	keepsorting                   =  1
	keepchecking                  =  1
	Tracker["current_generation"] = -1
	igen              = -1
	my_pids           = os.path.join(Tracker["constants"]["masterdir"], "indexes.txt")
	params            = os.path.join(Tracker["constants"]["masterdir"],"refinement_parameters.txt")
	previous_params   = Tracker["previous_parstack"]

	all_gen_stat_list = []
	bad_clustering    = 0
	while (keepsorting == 1) and (bad_clustering== 0):
		Tracker["current_generation"] +=1
		igen +=1
		work_dir = os.path.join(Tracker["constants"]["masterdir"], "generation_%03d"%igen)
		if Blockdata["myid"] == Blockdata["main_node"]:
			keepchecking = check_sorting_state(work_dir, keepchecking, log_main)
		else: keepchecking = 0
		keepchecking = bcast_number_to_all(keepchecking, Blockdata["main_node"], MPI_COMM_WORLD)
		if keepchecking == 0: # new, do it
			if Blockdata["myid"] == Blockdata["main_node"]:
				time_generation_start = time()
				if not os.path.exists(work_dir):
					os.mkdir(work_dir)# need check each box
					within_generation_restart =   0
				else: within_generation_restart = 1
				freq_cutoff_dict = {}
				fout = open(os.path.join(work_dir, "freq_cutoff.json"),'w')
				json.dump(freq_cutoff_dict, fout)
				fout.close()
				log_main.add('================================================================================================================')
				log_main.add('                                    SORT3D MULTI-LAYER   generation %d'%igen)
				log_main.add('----------------------------------------------------------------------------------------------------------------')
			else: within_generation_restart = 0
			within_generation_restart       = bcast_number_to_all(within_generation_restart, Blockdata["main_node"], MPI_COMM_WORLD)
			if within_generation_restart == 1: read_tracker_mpi(work_dir)
			else: dump_tracker(work_dir)
			output_list, bad_clustering, stat_list = depth_clustering(work_dir, depth_order, my_pids, params, previous_params, log_main)
			all_gen_stat_list.append(stat_list)
			if bad_clustering !=1:
				if Blockdata["myid"] == Blockdata["main_node"]:
					clusters, nclusters, nuacc  = output_clusters(work_dir, output_list[0][0], output_list[0][1], not_include_unaccounted, log_main)
					try:  del Tracker["generation"][str(igen)]
					except: pass 
					Tracker["generation"][igen] = len(clusters)
				else: 
					Tracker   = 0
					nclusters = 0
					nuacc     = 0
				Tracker   = wrap_mpi_bcast(Tracker,        Blockdata["main_node"], MPI_COMM_WORLD)
				nclusters = bcast_number_to_all(nclusters, Blockdata["main_node"], MPI_COMM_WORLD)
				nuacc     = bcast_number_to_all(nuacc,     Blockdata["main_node"], MPI_COMM_WORLD)
				dump_tracker(work_dir)
				if nclusters == 0:
					if Blockdata["myid"] == Blockdata["main_node"]:
						log_main.add('No cluster is found in generation %d'%igen)
					break
				else:
					if Blockdata["myid"] == Blockdata["main_node"]:
						time_of_generation_h,  time_of_generation_m = get_time(time_generation_start)
						log_main.add('SORT3D generation%d time: %d hours %d minutes.'%(igen, time_of_generation_h, time_of_generation_m))
					my_pids = os.path.join( work_dir, 'indexes_next_generation.txt')
					if Blockdata["myid"] == Blockdata["main_node"]:
						write_text_file(output_list[0][1], my_pids)
						write_text_row(stat_list, os.path.join(work_dir, 'gen_rep.txt'))
						mark_sorting_state(work_dir, True, log_main)
					mpi_barrier(MPI_COMM_WORLD)
				keepsorting = check_sorting(nuacc, keepsorting, log_main)
				if keepsorting != 1: break #
		else: # restart run
			read_tracker_mpi(work_dir)
			my_pids  = os.path.join( work_dir, 'indexes_next_generation.txt') 
			if Blockdata["myid"] == Blockdata["main_node"]: 
				stat_list = read_text_row( os.path.join(work_dir, 'gen_rep.txt'))
			else: stat_list = 0
			stat_list = wrap_mpi_bcast(stat_list, Blockdata["main_node"], MPI_COMM_WORLD)
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

#  End of various utilities	
def main():
	from optparse   import OptionParser
	from global_def import SPARXVERSION
	from EMAN2      import EMData
	from logger     import Logger, BaseLogger_Files
	from global_def import ERROR
	import sys, os, time, shutil
	global Tracker, Blockdata
	progname = os.path.basename(sys.argv[0])
	usage = progname + " --refinement_dir=masterdir_of_sxmeridien   --output_dir=sort3d_output --mask3D=mask.hdf --focus=binarymask.hdf  --radius=outer_radius " +\
	"  --sym=c1  --img_per_grp=img_per_grp  --minimum_grp_size=minimum_grp_size "
	parser = OptionParser(usage,version=SPARXVERSION)
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
			print("Specify one of the two options to start the program: --refinement_dir, --instack")
	
	if initiate_from_meridien_mode:
		parser.add_option("--output_dir",                        type   ="string",        default ='',					   help="Sort3d output directory name")
		parser.add_option("--niter_for_sorting",                 type   ="int",           default =-1,					   help="User specified iteration number of 3D refinement for sorting")
		parser.add_option("--focus",                             type   ="string",        default ='',                     help="Focus 3D mask. File path of a binary 3D mask for focused clustering ")
		parser.add_option("--mask3D",                            type   ="string",        default ='',                     help="3D mask. File path of the global 3D mask for clustering")
		parser.add_option("--radius",                            type   ="int",           default =-1,	                   help="Estimated protein radius in pixels")
		parser.add_option("--sym",                               type   ="string",        default ='c1',                   help="Point-group symmetry")
		parser.add_option("--img_per_grp",                       type   ="int",           default =1000,                   help="Number of images per group")
		parser.add_option("--img_per_grp_split_rate",            type   ="int",           default =1,                      help="rate for splitting img_per_grp")
		parser.add_option("--nsmear",                            type   ="float",         default =-1.,                    help="Number of smears used in sorting. Fill it with 1 if user does not want to use all smears")
		parser.add_option("--minimum_grp_size",				     type   ="int",           default =-1,					   help="Cluster selection size")
		parser.add_option("--depth_order",				         type   ="int",           default =2,					   help="Depth order. A number defines the number of initial independent MGSKmeans runs (2^depth_order)")
		parser.add_option("--memory_per_node",						type   ="float",		default =-1.0,		help="Memory_per_node, the number used for computing the CPUs/NODE settings given by user")
		parser.add_option("--orientation_cones",					type   ="int",			default =-1,		help="Number of orientation cones in the asymmetric unit (default 20, adjusted for symmetry)")
		parser.add_option("--not_include_unaccounted",				action ="store_true",	default =False,		help="Do not reconstruct unaccounted elements in each generation")
		parser.add_option("--stop_mgskmeans_percentage",			type   ="float",		default =10.0,		help="Percentage of particles with changed group assignment. Converge criteria for Kmeans clustering with constrained group size")
		parser.add_option("--swap_ratio",							type   ="float",		default =5.0,		help="Randomness ratio of swapping accounted elements with unaccounted for elemetns per cluster. A float number between 0.0 and 50.0")
		parser.add_option("--notapplybckgnoise",					action ="store_true",	default =False,		help="Do not applynoise")
		parser.add_option("--do_swap_au",							action ="store_true",	default =False,		help="Flag to turn on swapping the accounted for images with the unaccounted for images")
		parser.add_option("--random_group_elimination_threshold",	type   ="float",		default =2.0,		help="Number of random group reproducibility standard deviation for eliminating random groups")
		parser.add_option("--num_core_set",                        type   ="int",         default =-1,					   help="Number of images for reconstructing core set images. Will not reconstruct core set images if their number if less")
		parser.add_option("--check_smearing",						action ="store_true",	default =False,		help="Check the smearing per iteration and estimate the precalculated data")
		parser.add_option("--compute_on_the_fly",					action ="store_true",	default =False,		help="Pre-compute part of multiple assignment images for reconstruction till memory is filled up and then the other part on the fly")
		parser.add_option("--override",								action ="store_true",	default =False,		help="Override sanity checks of input parameters.  To be used at user's discretion.")

		(options, args) = parser.parse_args(sys.argv[1:])
		from utilities import bcast_number_to_all
		### Sanity check
	
		checking_flag = 0
		if Blockdata["myid"] == Blockdata["main_node"]:
			if options.refinement_dir !='':
				if not os.path.exists(options.refinement_dir): checking_flag = 1
		checking_flag = bcast_number_to_all(checking_flag, Blockdata["main_node"], MPI_COMM_WORLD)
		if checking_flag ==1:
			ERROR("The specified refinement_dir does not exist", "sort3d", 1, Blockdata["myid"])
	
		if options.focus !='':
			if Blockdata["myid"] == Blockdata["main_node"]:
				if not os.path.exists(options.focus):  checking_flag = 1
			checking_flag = bcast_number_to_all(checking_flag, Blockdata["main_node"], MPI_COMM_WORLD)
			if checking_flag ==1:
				ERROR("The specified focus mask file does not exist", "sort3d", 1, Blockdata["myid"])
		
		if options.mask3D !='':
			if Blockdata["myid"] == Blockdata["main_node"]:
				if not os.path.exists(options.mask3D): checking_flag = 1
			checking_flag = bcast_number_to_all(checking_flag, Blockdata["main_node"], MPI_COMM_WORLD)
			if checking_flag ==1:
				ERROR("The specified mask3D file does not exist", "sort3d", 1, Blockdata["myid"])
		
		if options.img_per_grp <=1: ERROR("Incorrect input parameter for img_per_grp", "sort3d", 1, Blockdata["myid"])
		elif options.img_per_grp < options.minimum_grp_size: 
			ERROR("Img_per_grp should be always larger than minimum_grp_size", "sort3d", 1, Blockdata["myid"])
	
		#--- Fill input parameters into dictionary Constants
		Constants		                         = {}
		Constants["stop_mgskmeans_percentage"]   = options.stop_mgskmeans_percentage
		Constants["niter_for_sorting"]           = options.niter_for_sorting
		Constants["memory_per_node"]             = options.memory_per_node
		Constants["orgstack"]                    = options.instack
		Constants["masterdir"]                   = options.output_dir
		Constants["refinement_dir"]              = options.refinement_dir
		Constants["random_group_elimination_threshold"]  = options.random_group_elimination_threshold
	
		if options.mask3D == '': Constants["mask3D"] = False
		else:                    Constants["mask3D"] = options.mask3D
		if options.focus!='':   Constants["focus3D"] = options.focus
		else:                   Constants["focus3D"] = False
	
		Constants["depth_order"]                 = options.depth_order
		Constants["img_per_grp"]                 = options.img_per_grp
		Constants["img_per_grp_split_rate"]      = options.img_per_grp_split_rate
		Constants["minimum_grp_size"]      		 = options.minimum_grp_size
		Constants["radius"]              		 = options.radius
		Constants["sym"]                         = options.sym
		Constants["nsmear"]                      = options.nsmear

		Constants["restart_from_nbox"]           =  0 #options.restart_from_nbox
		Constants["restart_from_depth_order"]    = -1 #options.restart_from_depth_order
		Constants["restart_from_generation"]     = -1 #options.restart_from_generation
	
		#### options for advanced users
		Constants["do_swap_au"]                  = options.do_swap_au
		Constants["swap_ratio"]                  = options.swap_ratio
		Constants["not_include_unaccounted"]     = False
		Constants["nxinit"]                      = -1
		Constants["box_niter"]                   = 5 
	
		### Frozen options
		Constants["upscale"]                     = 0.5 #
		Constants["interpolation"]               = "trl"
		Constants["comparison_method"]           = "cross" #options.comparison_method # either cross or eucd
		Constants["symmetry"]                    = Constants["sym"]
		Constants["CTF"]                		 = True
		Constants["do_not_use_3dmask"]           = False 
	
		if options.focus:  Constants["comparison_method"] = "cross" # in case of focus3D, cross is used.
		Constants["fuse_freq"]			= 45.  # Now in A, convert to pixels before being used
		Constants["orientation_cones"]  = options.orientation_cones # orientation constrained angle step
		Constants["num_core_set"]        = options.num_core_set
		Constants["check_smearing"]      = options.check_smearing
		Constants["compute_on_the_fly"]  = options.compute_on_the_fly
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
	
		###=====<--options for advanced users:
		Tracker["total_number_of_iterations"] = 15
		Tracker["clean_volumes"]              = True # always true
		Tracker["compute_fsc"]                = False
		### -----------Orientation constraints
		Tracker["tilt1"]                =  0.0
		Tracker["tilt2"]                = 180.0
		Tracker["grp_size_relx_ratio"]  = 0.98
		Tracker["minimum_ptl_number"]   = 20
		### ------------=====< option for proteins images that have preferred orientations
		 # for orientation groups
		if    Tracker["constants"]["memory_per_node"] == -1 or Tracker["constants"]["memory_per_node"] <32.:
			Tracker["constants"]["small_memory"] = True
		else: Tracker["constants"]["small_memory"] = False
	
		## additional check
		Tracker["constants"]["hardmask"]          = True
		Tracker["applymask"]                      = True
		Tracker["constants"]["refinement_method"] = "SPARX" 
		Tracker["nosmearing"]                     = False
		checking_flag                             = 0 # reset
		Blockdata["fftwmpi"]                      = True
		Tracker["even_scale_fsc"]                 = True
		Tracker["minimum_img_per_cpu"]            =   5 # User can adjust it to other number
		
		Tracker["num_on_the_fly"] = [ -1  for im in range(Blockdata["nproc"])]
		Tracker["current_img_per_grp"] = Tracker["constants"]["img_per_grp"]
		Blockdata["symclass"] = symclass(Tracker["constants"]["symmetry"])
		#  Default was set to -1, in which case adjust fo symmetry, else use input directly
		if(Tracker["constants"]["orientation_cones"] < 3):
			Tracker["constants"]["orientation_cones"] = \
			    max(3, int(round(20.0/sqrt(Blockdata["symclass"].nsym))))
		
		ast = number_of_cones_to_delta(Tracker["constants"]["orientation_cones"])
		Blockdata["ncpuspernode"] = Blockdata["no_of_processes_per_group"]
		Blockdata["nsubset"]      = Blockdata["ncpuspernode"]*Blockdata["no_of_groups"]
		create_subgroup()
		create_zero_group()		
		#
		from statistics 	import k_means_match_clusters_asg_new,k_means_stab_bbenum
		from utilities 		import get_im,bcast_number_to_all, write_text_file,read_text_file,wrap_mpi_bcast, get_params_proj, write_text_row
		from filter			import filt_tanl
		from logger         import Logger,BaseLogger_Files
		import string
		import json
		import user_functions
		from string         import split, atoi, atof
		#
		continue_from_interuption = 0
		continue_from_interuption = create_masterdir()
		log_main                  = Logger(BaseLogger_Files())
		log_main.prefix           = Tracker["constants"]["masterdir"]+"/"
		if continue_from_interuption == 0:
			if Blockdata["myid"] == Blockdata["main_node"]:
				log_main.add('================================================================================================================')
				log_main.add('                                 SORT3D MULTI-LAYER v1.3')
				log_main.add('================================================================================================================')
			import_data(log_main, options.override)
			print_shell_command(sys.argv, log_main)
			check_3dmask(log_main)
			check_mpi_settings(log_main)
			keepsorting = sort3d_init("initialization", log_main)
			dump_tracker(Tracker["constants"]["masterdir"])
			if keepsorting == 0:
				from mpi import mpi_finalize
				mpi_finalize()
				exit()				 

		sorting_main_mpi(log_main, options.depth_order, options.not_include_unaccounted)
		if Blockdata["myid"] == Blockdata["main_node"]:
			log_main.add('----------------------------------------------------------------------------------------------------------------')
			log_main.add('                                 SORT3D MULTI-LAYER finished')
			log_main.add('================================================================================================================\n')
		from mpi import mpi_finalize
		mpi_finalize()
		exit()
			
	elif initiate_from_data_stack_mode:
		parser.add_option("--nxinit",                            type   ="int",           default =-1,                     help="User provided image size")
		parser.add_option("--output_dir",                        type   ="string",        default ='',					   help="Name of the sort3d directory")
		parser.add_option("--focus",                             type   ="string",        default ='',                     help="Focus 3D mask. File path of a binary 3D mask for focused clustering ")
		parser.add_option("--mask3D",                            type   ="string",        default ='',                     help="3D mask. File path of the global 3D mask for clustering")
		parser.add_option("--radius",                            type   ="int",           default =-1,	                   help="Estimated protein radius in pixels")
		parser.add_option("--sym",                               type   ="string",        default ='c1',                   help="Point-group symmetry")
		parser.add_option("--img_per_grp",                       type   ="int",           default =1000,                   help="Number of images per group")
		parser.add_option("--img_per_grp_split_rate",            type   ="int",           default =1,                      help="rate for splitting img_per_grp")
		parser.add_option("--nsmear",                            type   ="float",         default =-1.,                    help="Number of smears used in sorting. Fill it with 1 if user does not want to use all smears")
		parser.add_option("--minimum_grp_size",				     type   ="int",           default =-1,					   help="Cluster selection size")
		parser.add_option("--depth_order",				         type   ="int",           default =2,					   help="Depth order. A number defines the number of initial independent MGSKmeans runs (2^depth_order)")
		parser.add_option("--memory_per_node",                   type   ="float",         default =-1.0,                   help="Memory_per_node, the number used for computing the CPUs/NODE settings given by user")
		parser.add_option("--orientation_cones",			type   ="int",           default =20,		help="Number of orientation cones in the asymmetric unit")
		parser.add_option("--not_include_unaccounted",		action ="store_true",    default =False,    help="Do not reconstruct unaccounted elements in each generation")
		parser.add_option("--stop_mgskmeans_percentage",	type   ="float",         default =10.0,		help="Percentage of particles with changed group assignment. Converge criteria for Kmeans clustering with constrained group size")
		parser.add_option("--swap_ratio",                        type   ="float",         default =5.0,                    help="Randomness ratio of swapping accounted elements with unaccounted for elemetns per cluster. A float number between 0.0 and 50.0")
		parser.add_option("--notapplybckgnoise",                 action ="store_true",    default =False,                  help="Flag to turn off background noise")
		parser.add_option("--do_swap_au",                        action ="store_true",    default =False,                  help="Flag to turn on swapping the accounted for images with the unaccounted for images")
		parser.add_option("--random_group_elimination_threshold",  type   ="float",       default =2.0,                    help="Number of random group reproducibility standard deviation for eliminating random groups")
		parser.add_option("--num_core_set",                        type   ="int",         default =-1,					   help="Number of images for reconstructing core set images. Will not reconstruct core set images if their number if less")
		parser.add_option("--override",								action ="store_true",	default =False,		help="Override sanity checks of input parameters.  To be used at user's discretion.")

		(options, args) = parser.parse_args(sys.argv[1:])
		from utilities import bcast_number_to_all
		### Sanity check
		
		checking_flag = 0
		if options.focus !='':
			if Blockdata["myid"] == Blockdata["main_node"]:
				if not os.path.exists(options.focus):  checking_flag = 1
			checking_flag = bcast_number_to_all(checking_flag, Blockdata["main_node"], MPI_COMM_WORLD)
			if checking_flag ==1:
				ERROR("The specified focus mask file does not exist", "sort3d", 1, Blockdata["myid"])
		
		if options.mask3D !='':
			if Blockdata["myid"] == Blockdata["main_node"]:
				if not os.path.exists(options.mask3D): checking_flag = 1
			checking_flag = bcast_number_to_all(checking_flag, Blockdata["main_node"], MPI_COMM_WORLD)
			if checking_flag ==1:
				ERROR("The specified mask3D file does not exist", "sort3d", 1, Blockdata["myid"])
		
		if( options.img_per_grp <=1 ):
			ERROR("Incorrect number of images per group", "sort3d", 1, Blockdata["myid"])
		elif( options.img_per_grp < options.minimum_grp_size ):
			ERROR("Parameter img_per_grp should be always larger than parameter minimum_grp_size", "sort3d", 1, Blockdata["myid"])
	
		#--- Fill input parameters into dictionary Constants
		Constants		                         = {}
		Constants["stop_mgskmeans_percentage"]   = options.stop_mgskmeans_percentage
		Constants["memory_per_node"]             = options.memory_per_node
		Constants["orgstack"]                    = options.instack
		Constants["masterdir"]                   = options.output_dir
		Constants["random_group_elimination_threshold"]  = options.random_group_elimination_threshold
		
		if options.mask3D == '': Constants["mask3D"] = False
		else:                    Constants["mask3D"] = options.mask3D
		if options.focus!='':   Constants["focus3D"] = options.focus
		else:                   Constants["focus3D"] = False
		
		Constants["nsmear"]                      = 1
		Constants["depth_order"]                 = options.depth_order
		Constants["img_per_grp"]                 = options.img_per_grp
		Constants["img_per_grp_split_rate"]      = options.img_per_grp_split_rate
		Constants["minimum_grp_size"]      		 = options.minimum_grp_size
		Constants["radius"]              		 = options.radius
		Constants["sym"]                         = options.sym
	
		Constants["restart_from_nbox"]           =  0 #options.restart_from_nbox
		Constants["restart_from_depth_order"]    = -1 #options.restart_from_depth_order
		Constants["restart_from_generation"]     = -1 #options.restart_from_generation
	
		#### options for advanced users
		Constants["do_swap_au"]                  = options.do_swap_au
		Constants["swap_ratio"]                  = options.swap_ratio
		Constants["not_include_unaccounted"]     = False
		Constants["nxinit"]                      = options.nxinit
		Constants["box_niter"]                   = 5 
		
		### Frozen options
		Constants["upscale"]                      = 0.5 #
		Constants["interpolation"]                = "trl"
		Constants["comparison_method"]            = "cross" #options.comparison_method # either cross or eucd
		Constants["symmetry"]                     = Constants["sym"]
		Constants["CTF"]                		  = True
		Constants["do_not_use_3dmask"]            = False 
	
		if options.focus:  Constants["comparison_method"] = "cross" # in case of focus3D, cross is used.
		Constants["fuse_freq"] = 45.  # Now in A, convert to pixels before being used
		Constants["orientation_cones"]  = options.orientation_cones # orientation constrained angle step
		Constants["num_core_set"]        = options.num_core_set
		Constants["compute_on_the_fly"]  = False
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
	
		###=====<--options for advanced users:
		Tracker["total_number_of_iterations"] = 15
		Tracker["clean_volumes"]              = True  # always true
		Tracker["compute_fsc"]                = False
		### -----------Orientation constraints
		Tracker["tilt1"]                =  0.0
		Tracker["tilt2"]                = 180.0
		Tracker["grp_size_relx_ratio"]  = 0.98
		Tracker["minimum_ptl_number"]   = 20
		### ------------=====< option for proteins images that have preferred orientations
		 # for orientation groups
		if    Tracker["constants"]["memory_per_node"] == -1 or Tracker["constants"]["memory_per_node"] <32.:
			Tracker["constants"]["small_memory"]   = True
		else: Tracker["constants"]["small_memory"] = False
	
		## additional check
		Tracker["constants"]["hardmask"]          = True
		Tracker["applymask"]                      = True
		Tracker["constants"]["refinement_method"] ="stack"
		Tracker["constants"]["refinement_dir"]    = None
		Tracker["paramstructure_dir"]             = None
		Tracker["refang"]                         = None
		Tracker["rshifts"]                        = None
		Tracker["paramstructure_dict"]            = None
		Tracker["constants"]["selected_iter"]     = -1
		Tracker["nosmearing"]                     = True
		Tracker["even_scale_fsc"]                 = True
		Tracker["minimum_img_per_cpu"]            =   5 # User can adjust it to other number
		checking_flag                              =  0 # reset
		Blockdata["fftwmpi"]                       = True
		Blockdata["symclass"]                      = symclass(Tracker["constants"]["symmetry"])
		#  Default was set to -1, in which case adjust fo symmetry, else use input directly
		if(Tracker["constants"]["orientation_cones"] < 3):
			Tracker["constants"]["orientation_cones"] = \
			    max(3, int(round(20.0/sqrt(Blockdata["symclass"].nsym))))
		
		Tracker["num_on_the_fly"]      = [ -1 for im in range(Blockdata["nproc"])]
		Tracker["current_img_per_grp"] = Tracker["constants"]["img_per_grp"]
		
		ast = number_of_cones_to_delta(Tracker["constants"]["orientation_cones"])
		Blockdata["ncpuspernode"] = Blockdata["no_of_processes_per_group"]
		Blockdata["nsubset"]      = Blockdata["ncpuspernode"]*Blockdata["no_of_groups"]
		create_subgroup()
		create_zero_group()
		# 
		from statistics 	import k_means_match_clusters_asg_new,k_means_stab_bbenum
		from utilities 		import get_im,bcast_number_to_all, write_text_file,read_text_file,wrap_mpi_bcast, get_params_proj, write_text_row
		from filter			import filt_tanl
		from logger         import Logger,BaseLogger_Files
		import string
		import json
		import user_functions
		from string         import split, atoi, atof
	
		continue_from_interuption = 0
		continue_from_interuption = create_masterdir()
		log_main = Logger(BaseLogger_Files())
		log_main.prefix = Tracker["constants"]["masterdir"]+"/"
		if continue_from_interuption == 0:
			if Blockdata["myid"] == Blockdata["main_node"]:
				log_main.add('================================================================================================================')
				log_main.add('                                  SORT3D MULTI-LAYER v1.3')
				log_main.add('================================================================================================================\n')	
			import_data(log_main, options.override)
			print_shell_command(sys.argv, log_main)
			check_3dmask(log_main)
			check_mpi_settings(log_main)
			keepsorting = sort3d_init("Initialization", log_main)
			dump_tracker(Tracker["constants"]["masterdir"])
			if keepsorting ==0:
				from mpi import mpi_finalize
				mpi_finalize()
				exit()
		sorting_main_mpi(log_main, options.depth_order, options.not_include_unaccounted)
		if Blockdata["myid"] == Blockdata["main_node"]:
			log_main.add('----------------------------------------------------------------------------------------------------------------')
			log_main.add('                                 SORT3D MULTI-LAYER finished')
			log_main.add('================================================================================================================\n')
		from mpi import mpi_finalize
		mpi_finalize()
		exit()
if __name__ == "__main__":
	main()
