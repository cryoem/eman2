#!/usr/bin/env python
#
#
#  08/26/2016
#  New version of sort3D.
#
from __future__ import print_function
from __future__ import division
from past.utils import old_div
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
#
#

import EMAN2_cppwrap
import copy
import json
import mpi
import numpy
import optparse
import os
import random
import scipy
import shutil
from ..libpy import sp_applications
from ..libpy import sp_filter
from ..libpy import sp_fundamentals
from ..libpy import sp_global_def
from ..libpy import sp_logger
from ..libpy import sp_morphology
from ..libpy import sp_projection
from ..libpy import sp_statistics
from ..libpy import sp_utilities
import string
import sys
import time
from builtins import range

"""
There are two ways to run the program:

1.Import data from a meridien refinement.
mpirun -np 64 --hostfile four_nodes.txt  sxsort3d_depth.py  --refinemet_dir=meridien_run --niter_for_sorting=28  --memory_per_node=100. --img_per_grp=80000  --minimum_grp_size=15000   --stop_mgskmeans_percentage=8. --output_dir=SORT3D

2.Import data from a given data stack.
mpirun  -np 48  --hostfile ./node012.txt  sxsort3d_depth.py --stop_mgskmeans_percentage=15. --orientation_groups=40  --do_swap_au  --swap_ratio=5. --output_dir=sorting_bmask04 --sym=c1   --radius=30  --minimum_grp_size=2000   --img_per_grp=2800    --instack=bdb:data  >sorting_bmask04/printout &

Notices on options:

a. --do_swap_au  --swap_ratio=5.   : Ratio of the elements of determined clusters that are exchanged with the unaccounted elements.
b. --stop_mgskmeans_percentage=15. : criterion to stop Kmeans run.

Nomenclatures in sorting intermediate results:

NACC:  total number of accounted
NUACC: number of unaccounted
MGU:   user defined minimum_grp_size  
MGR:   random group siz
K:     number of groups

"""


global Tracker, Blockdata

mpi.mpi_init(0, [])
nproc = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
Blockdata = {}
#  MPI stuff
Blockdata["nproc"] = nproc
Blockdata["myid"] = myid
Blockdata["main_node"] = 0
Blockdata["last_node"] = nproc - 1

Blockdata["shared_comm"] = mpi.mpi_comm_split_type(
    mpi.MPI_COMM_WORLD, mpi.MPI_COMM_TYPE_SHARED, 0, mpi.MPI_INFO_NULL
)
Blockdata["myid_on_node"] = mpi.mpi_comm_rank(Blockdata["shared_comm"])
Blockdata["no_of_processes_per_group"] = mpi.mpi_comm_size(Blockdata["shared_comm"])
masters_from_groups_vs_everything_else_comm = mpi.mpi_comm_split(
    mpi.MPI_COMM_WORLD,
    Blockdata["main_node"] == Blockdata["myid_on_node"],
    Blockdata["myid_on_node"],
)
Blockdata["color"], Blockdata[
    "no_of_groups"
], balanced_processor_load_on_nodes = sp_utilities.get_colors_and_subsets(
    Blockdata["main_node"],
    mpi.MPI_COMM_WORLD,
    Blockdata["myid"],
    Blockdata["shared_comm"],
    Blockdata["myid_on_node"],
    masters_from_groups_vs_everything_else_comm,
)

#  We need two nodes for processing of volumes
if Blockdata["no_of_groups"] > 1:
    Blockdata["node_volume"] = [
        Blockdata["no_of_groups"] - 2,
        Blockdata["no_of_groups"] - 1,
    ]
# Blockdata["nodes"] = [Blockdata["no_of_groups"]-2, Blockdata["no_of_groups"]-1]  # For 3D stuff take last two nodes
else:
    Blockdata["node_volume"] = [0, 0]
#  We need two CPUs for processing of volumes, they are taken to be main CPUs on each volume
#  We have to send the two myids to all nodes so we can identify main nodes on two selected groups.
if Blockdata["no_of_groups"] > 1:
    Blockdata["main_shared_nodes"] = [
        Blockdata["node_volume"][0] * Blockdata["no_of_processes_per_group"],
        Blockdata["node_volume"][1] * Blockdata["no_of_processes_per_group"],
    ]
else:
    Blockdata["main_shared_nodes"] = [0, 1]
Blockdata["nproc_previous"] = 0
# End of Blockdata: sorting requires at least three nodes, and the used number of nodes be integer times of three
sp_global_def.BATCH = True
sp_global_def.MPI = True
global _proc_status, _scale, is_unix_cluster
try:
    _proc_status = "/proc/%d/status" % os.getpid()
    _scale = {"kB": 1024.0, "mB": 1024.0 * 1024.0, "KB": 1024.0, "MB": 1024.0 * 1024.0}
    is_unix_cluster = True
except:
    if Blockdata["myid"] == Blockdata["main_node"]:
        print("Not a unix machine")
    is_unix_cluster = False


def create_subgroup():
    # select a subset of myids to be in subdivision
    if Blockdata["myid_on_node"] < Blockdata["ncpuspernode"]:
        submyids = [Blockdata["myid"]]
    else:
        submyids = []
    submyids = sp_utilities.wrap_mpi_gatherv(
        submyids, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    submyids = sp_utilities.wrap_mpi_bcast(
        submyids, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    world_group = mpi.mpi_comm_group(mpi.MPI_COMM_WORLD)
    subgroup = mpi.mpi_group_incl(world_group, len(submyids), submyids)
    Blockdata["subgroup_comm"] = mpi.mpi_comm_create(mpi.MPI_COMM_WORLD, subgroup)
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    Blockdata["subgroup_size"] = -1
    Blockdata["subgroup_myid"] = -1
    if mpi.MPI_COMM_NULL != Blockdata["subgroup_comm"]:
        Blockdata["subgroup_size"] = mpi.mpi_comm_size(Blockdata["subgroup_comm"])
        Blockdata["subgroup_myid"] = mpi.mpi_comm_rank(Blockdata["subgroup_comm"])
    #  "nodes" are zero nodes on subgroups on the two "node_volume" that compute backprojection
    Blockdata["nodes"] = [
        Blockdata["node_volume"][0] * Blockdata["ncpuspernode"],
        Blockdata["node_volume"][1] * Blockdata["ncpuspernode"],
    ]
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    return


def create_zero_group():
    # select a subset of myids to be in subdivision, This is a group of all zero IDs on nodes, taken from isac2
    if Blockdata["myid_on_node"] == 0:
        submyids = [Blockdata["myid"]]
    else:
        submyids = []

    submyids = sp_utilities.wrap_mpi_gatherv(
        submyids, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    submyids = sp_utilities.wrap_mpi_bcast(
        submyids, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    world_group = mpi.mpi_comm_group(mpi.MPI_COMM_WORLD)
    subgroup = mpi.mpi_group_incl(world_group, len(submyids), submyids)
    Blockdata["group_zero_comm"] = mpi.mpi_comm_create(mpi.MPI_COMM_WORLD, subgroup)
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    Blockdata["group_zero_size"] = -1
    Blockdata["group_zero_myid"] = -1
    if mpi.MPI_COMM_NULL != Blockdata["group_zero_comm"]:
        Blockdata["group_zero_size"] = mpi.mpi_comm_size(Blockdata["group_zero_comm"])
        Blockdata["group_zero_myid"] = mpi.mpi_comm_rank(Blockdata["group_zero_comm"])
    #  "nodes" are zero nodes on subgroups on the two "node_volume" that compute backprojection
    # Blockdata["nodes"] = [Blockdata["node_volume"][0]*Blockdata["ncpuspernode"], Blockdata["node_volume"][1]*Blockdata["ncpuspernode"]]
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    return


### utilities


def dump_tracker(path_of_the_tracker):
    global Tracker, Blockdata
    if Blockdata["myid"] == Blockdata["main_node"]:
        fout = open(os.path.join(path_of_the_tracker, "Tracker.json"), "w")
        fout.seek(0)
        json.dump(Tracker, fout)
        fout.truncate()
        fout.close()
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    return


######### depth clustering functions


def depth_clustering(
    work_dir, depth_order, initial_id_file, params, previous_params, log_main
):
    global Tracker, Blockdata
    keepchecking = 1
    bad_clustering = 0
    stat_list = []
    init_layer_dir = os.path.join(work_dir, "layer0")
    if Blockdata["myid"] == Blockdata["main_node"]:
        if not os.path.exists(init_layer_dir):
            os.mkdir(init_layer_dir)

    if Blockdata["myid"] == Blockdata["main_node"]:
        partition_per_box_per_layer_list = []
        initial_id_list = sp_utilities.read_text_file(initial_id_file)
        for iparti in range(0, 2 ** (depth_order + 1), 2):
            reassign = [initial_id_list, None]
            partition_per_box_per_layer_list.append(reassign)
    else:
        partition_per_box_per_layer_list = 0
    partition_per_box_per_layer_list = sp_utilities.wrap_mpi_bcast(
        partition_per_box_per_layer_list, Blockdata["main_node"]
    )

    Tracker["depth"] = 0
    for depth in range(
        depth_order
    ):  # layers, depth_order = 1 means one layer and two boxes.
        time_layer_start = time.time()
        n_cluster_boxes = 2 ** (depth_order - depth)
        depth_dir = os.path.join(work_dir, "layer%d" % depth)
        Tracker["depth"] = depth
        if Blockdata["myid"] == Blockdata["main_node"]:
            log_main.add(
                "                        Layer %d has %d pairs of independent sorting runs"
                % (depth, n_cluster_boxes)
            )
            log_main.add(" ")
            log_main.add(" ")
            if not os.path.exists(depth_dir):
                os.mkdir(depth_dir)
                keepchecking = 0
                mark_sorting_state(depth_dir, False, log_main)
            else:
                keepchecking = check_sorting_state(depth_dir, keepchecking, log_main)
                if keepchecking == 0:
                    mark_sorting_state(depth_dir, False, log_main)
        else:
            keepchecking = 0
        keepchecking = sp_utilities.bcast_number_to_all(
            keepchecking, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        ## box loop
        if keepchecking == 0:
            checkingbox = 1
            Tracker["box_nxinit"] = -1
            Tracker["box_nxinit_freq"] = -1.0
            for nbox in range(n_cluster_boxes):
                input_accounted_file = partition_per_box_per_layer_list[nbox][0]
                input_unaccounted_file = partition_per_box_per_layer_list[nbox][1]
                nbox_dir = os.path.join(depth_dir, "nbox%d" % nbox)
                if Blockdata["myid"] == Blockdata["main_node"]:
                    if not os.path.exists(nbox_dir):
                        os.mkdir(nbox_dir)
                        checkingbox = 0
                        mark_sorting_state(nbox_dir, False, log_main)
                    else:
                        checkingbox = check_sorting_state(
                            nbox_dir, checkingbox, log_main
                        )
                        if checkingbox == 0:  # found not finished box
                            msg = (
                                "Pair %d is not finished. Remove it and recompute..."
                                % nbox
                            )
                            log_main.add(msg)
                            shutil.rmtree(nbox_dir)
                            os.mkdir(nbox_dir)
                            mark_sorting_state(nbox_dir, False, log_main)
                        else:
                            msg = "Pair %d is completed" % nbox
                            log_main.add(msg)
                else:
                    checkingbox = 0
                checkingbox = sp_utilities.bcast_number_to_all(
                    checkingbox, Blockdata["main_node"], mpi.MPI_COMM_WORLD
                )
                # Code structure of the box
                if checkingbox == 0:
                    bad_box = depth_clustering_box(
                        nbox_dir,
                        input_accounted_file,
                        input_unaccounted_file,
                        params,
                        previous_params,
                        nbox,
                        log_main,
                    )
                    if Blockdata["myid"] == Blockdata["main_node"]:
                        mark_sorting_state(nbox_dir, True, log_main)
                else:
                    continue
                if bad_box == 1:
                    bad_clustering = 1
            if bad_clustering != 1:
                partition_per_box_per_layer_list = []
                stop_generation = 0
                for nbox in range(0, n_cluster_boxes, 2):
                    input_box_parti1 = os.path.join(
                        depth_dir, "nbox%d" % nbox, "partition.txt"
                    )
                    input_box_parti2 = os.path.join(
                        depth_dir, "nbox%d" % (nbox + 1), "partition.txt"
                    )
                    minimum_grp_size, maximum_grp_size, accounted_list, unaccounted_list, bad_clustering, stop_generation, stat_list = do_boxes_two_way_comparison_mpi(
                        nbox,
                        input_box_parti1,
                        input_box_parti2,
                        depth_order - depth,
                        log_main,
                    )

                    if stop_generation == 1:
                        partition_per_box_per_layer_list = []
                        partition_per_box_per_layer_list.append(
                            [accounted_list, unaccounted_list]
                        )
                        break
                    else:
                        partition_per_box_per_layer_list.append(
                            [accounted_list, unaccounted_list]
                        )
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            if Blockdata["myid"] != Blockdata["main_node"]:
                Tracker = 0
            Tracker = sp_utilities.wrap_mpi_bcast(
                Tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )

            if Blockdata["myid"] == Blockdata["main_node"]:
                mark_sorting_state(depth_dir, True, log_main)
            if bad_clustering == 1:
                break
            if stop_generation == 1:
                break  ### only one cluster survives
        else:
            if Blockdata["myid"] == Blockdata["main_node"]:
                log_main.add(
                    "Layer %d sorting has completed. Recompute two_way comparison"
                    % depth
                )
            partition_per_box_per_layer_list = []
            for nbox in range(0, n_cluster_boxes, 2):
                input_box_parti1 = os.path.join(
                    depth_dir, "nbox%d" % nbox, "partition.txt"
                )
                input_box_parti2 = os.path.join(
                    depth_dir, "nbox%d" % (nbox + 1), "partition.txt"
                )
                minimum_grp_size, maximum_grp_size, accounted_list, unaccounted_list, bad_clustering, stop_generation, stat_list = do_boxes_two_way_comparison_mpi(
                    nbox,
                    input_box_parti1,
                    input_box_parti2,
                    depth_order - depth,
                    log_main,
                )

                if stop_generation == 1:
                    partition_per_box_per_layer_list = []
                    partition_per_box_per_layer_list.append(
                        [accounted_list, unaccounted_list]
                    )
                    break
                else:
                    partition_per_box_per_layer_list.append(
                        [accounted_list, unaccounted_list]
                    )
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

        time_of_sorting_h, time_of_sorting_m = get_time(time_layer_start)
        if Blockdata["myid"] == Blockdata["main_node"]:
            log_main.add(
                "            Execution of layer %d took %d hours %d minutes."
                % (depth, time_of_sorting_h, time_of_sorting_m)
            )
            log_main.add(" ")
            log_main.add(" ")
            log_main.add(
                "================================================================================================================"
            )
    return partition_per_box_per_layer_list, bad_clustering, stat_list


def depth_box_initialization(box_dir, input_list1, input_list2, log_file):
    global Tracker, Blockdata
    img_per_grp = Tracker["current_img_per_grp"]
    minimum_grp_size = Tracker["constants"]["minimum_grp_size"]

    if input_list2 is not None:  # Track 2
        total_stack = len(input_list1) + len(input_list2)
        groups = []
        for one in input_list1:
            if one[0] not in groups:
                groups.append(
                    one[0]
                )  # safe in speed when the number of groups is not large.
        number_of_groups = len(groups)
        if Blockdata["myid"] == Blockdata["main_node"]:
            msg = (
                "Number of groups found in initialization: %d, the largest possible number of groups: %d"
                % (number_of_groups, old_div(total_stack, img_per_grp))
            )
            log_file.add(msg)
            sp_utilities.write_text_row(
                input_list1, os.path.join(box_dir, "previous_NACC.txt")
            )
            sp_utilities.write_text_file(
                input_list2, os.path.join(box_dir, "previous_NUACC.txt")
            )
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

        if Tracker["constants"]["do_swap_au"]:
            swap_ratio = Tracker["constants"]["swap_ratio"]
        else:
            swap_ratio = 0.0
        new_assignment = []

        for indep in range(2):
            tmp_assignment = swap_accounted_with_unaccounted_elements_mpi(
                os.path.join(box_dir, "previous_NACC.txt"),
                os.path.join(box_dir, "previous_NUACC.txt"),
                log_file,
                number_of_groups,
                swap_ratio,
            )
            new_assignment.append(tmp_assignment)

        if Blockdata["myid"] == Blockdata["main_node"]:
            for indep in range(2):
                sp_utilities.write_text_file(
                    new_assignment[indep],
                    os.path.join(box_dir, "independent_index_%03d.txt" % indep),
                )
            new_assignment = []
            for indep in range(2):
                new_assignment.append(
                    sp_utilities.read_text_row(
                        os.path.join(box_dir, "independent_index_%03d.txt" % indep)
                    )
                )
            id_list = sp_utilities.read_text_file(
                os.path.join(box_dir, "independent_index_000.txt"), -1
            )
            if len(id_list) > 1:
                id_list = id_list[0]
            total_stack = len(id_list)
            number_of_groups = max(id_list) + 1  # assume K be 0, ...,number_of_groups-1
        else:
            number_of_groups = 0
            total_stack = 0
            new_assignment = 0
        new_assignment = sp_utilities.wrap_mpi_bcast(
            new_assignment, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        number_of_groups = sp_utilities.bcast_number_to_all(
            number_of_groups, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        total_stack = sp_utilities.bcast_number_to_all(
            total_stack, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )

    else:  # Tracker 1
        total_stack = len(input_list1)
        if Blockdata["myid"] == Blockdata["main_node"]:
            sp_utilities.write_text_file(
                input_list1, os.path.join(box_dir, "previous_all_indexes.txt")
            )

        number_of_groups = old_div(total_stack, Tracker["current_img_per_grp"])
        if number_of_groups < 2:
            sp_global_def.ERROR(
                "number_of_groups should be at least larger than 1",
                "depth_box_initialization",
                1,
                Blockdata["myid"],
            )

        new_assignment = create_nrandom_lists_from_given_pids(
            box_dir,
            os.path.join(box_dir, "previous_all_indexes.txt"),
            number_of_groups,
            2,
        )

    if Blockdata["myid"] == Blockdata["main_node"]:
        log_file.add(
            "Sorting settings:  Number of images: %d.  Number of groups: %d.  Minimum group size: %d."
            % (total_stack, number_of_groups, Tracker["constants"]["minimum_grp_size"])
        )
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    return img_per_grp, number_of_groups, total_stack, minimum_grp_size, new_assignment


def depth_iter_initialization(run_id_file):
    id_list = sp_utilities.read_text_file(run_id_file, -1)
    number_of_groups = max(id_list[0]) + 1
    total_stack = len(id_list[0])
    return total_stack, number_of_groups


def output_iter_results(
    box_dir,
    ncluster,
    NACC,
    NUACC,
    minimum_grp_size,
    list_of_stable,
    unaccounted_list,
    img_per_grp,
    log_main,
):
    ### single node
    iter_cutoff = min(minimum_grp_size, int(img_per_grp * 0.75))
    new_list = []
    nc = 0
    NACC = 0
    try:
        fout = open(os.path.join(current_dir, "freq_cutoff.json"), "r")
        freq_cutoff_dict = sp_utilities.convert_json_fromunicode(json.load(fout))
        fout.close()
    except:
        freq_cutoff_dict = {}

    for index_of_any in range(len(list_of_stable)):
        any = list_of_stable[index_of_any]
        any.tolist()
        if len(any) >= iter_cutoff:
            any.sort()
            new_list.append(any)
            sp_utilities.write_text_file(
                any, os.path.join(box_dir, "Cluster_%03d.txt" % ncluster)
            )
            freq_cutoff_dict["Cluster_%03d.txt" % ncluster] = Tracker[
                "freq_fsc143_cutoff"
            ]
            ncluster += 1
            nc += 1
            NACC += len(any)
        else:
            for element in any:
                unaccounted_list.append(element)
    unaccounted_list.sort()
    NUACC = len(unaccounted_list)
    fout = open(os.path.join(box_dir, "freq_cutoff.json"), "w")
    json.dump(freq_cutoff_dict, fout)
    fout.close()
    return ncluster, NACC, NUACC, unaccounted_list, nc


def check_state_within_box_run(keepgoing, nruns, unaccounted_list, no_cluster_last_run):
    global Tracker, Blockdata
    total_stack = len(unaccounted_list)
    number_of_groups = old_div(total_stack, Tracker["current_img_per_grp"])
    if number_of_groups >= 2:
        keepgoing = 1
    else:
        keepgoing = 0
    if keepgoing == 1:
        nruns += 1
    else:
        total_stack = 0
        number_of_groups = 0
    if no_cluster_last_run:
        number_of_groups -= 1
    if number_of_groups <= 1:
        keepgoing = 0
    return keepgoing, nruns, total_stack, number_of_groups


def get_box_partition(box_dir, ncluster, unaccounted_list):
    if ncluster >= 1:
        unaccounted_list.sort()
        clusters_in_box = []
        for ic in range(ncluster):
            one_cluster = sp_utilities.read_text_file(
                os.path.join(box_dir, "Cluster_%03d.txt" % ic)
            )
            clusters_in_box.append(one_cluster)
        if len(unaccounted_list) > 0:
            clusters_in_box.append(unaccounted_list)
        alist, plist = merge_classes_into_partition_list(clusters_in_box)
    else:
        plist = []
    return plist


def output_clusters(
    output_dir, partition, unaccounted_list, not_include_unaccounted, log_main
):
    global Tracker, Blockdata
    ### single cpu function
    nclasses, npart = split_partition_into_ordered_clusters(partition)
    nc = 0
    identified_clusters = []
    for ic in range(len(nclasses)):
        if len(nclasses[ic]) >= max(
            old_div(Tracker["constants"]["img_per_grp"], 2)
            ** max(Tracker["constants"]["img_per_grp_split_rate"], 1),
            min(
                Tracker["constants"]["minimum_grp_size"],
                Tracker["constants"]["img_per_grp"] * 0.75,
            ),
        ):
            sp_utilities.write_text_file(
                nclasses[ic], os.path.join(output_dir, "Cluster_%03d.txt" % nc)
            )
            nc += 1
            identified_clusters.append(nclasses[ic])
        else:
            unaccounted_list += nclasses[ic]
    if len(unaccounted_list) > 1:
        unaccounted_list.sort()
        sp_utilities.write_text_file(
            unaccounted_list, os.path.join(output_dir, "Unaccounted.txt")
        )
    nclasses = copy.deepcopy(identified_clusters)
    del identified_clusters

    if len(unaccounted_list) > 1:  # output unaccounted as the last cluster
        if not not_include_unaccounted:
            sp_utilities.write_text_file(
                unaccounted_list, os.path.join(output_dir, "Cluster_%03d.txt" % nc)
            )
    if not not_include_unaccounted:
        unclasses = copy.deepcopy(nclasses)
        unclasses.append(unaccounted_list)
        alist, partition = merge_classes_into_partition_list(unclasses)
        del unclasses
    else:
        alist, partition = merge_classes_into_partition_list(nclasses)
    sp_utilities.write_text_row(
        partition, os.path.join(output_dir, "final_partition.txt")
    )
    return nclasses, nc, len(unaccounted_list)


def do_analysis_on_identified_clusters(clusters, log_main):
    global Tracker, Blockdata
    if Tracker["nosmearing"]:
        vs, ds, ss, norms = get_params_for_analysis(
            Tracker["constants"]["orgstack"],
            os.path.join(
                Tracker["constants"]["masterdir"], "refinement_parameters.txt"
            ),
            None,
            None,
        )
    else:
        vs, ds, ss, norms = get_params_for_analysis(
            Tracker["constants"]["orgstack"],
            os.path.join(
                Tracker["constants"]["masterdir"], "refinement_parameters.txt"
            ),
            os.path.join(Tracker["constants"]["masterdir"], "all_smearing.txt"),
            Tracker["constants"]["nsmear"],
        )

    tmpres1, tmpres2 = do_one_way_anova_scipy(
        clusters, ds, name_of_variable="defocus", log_main=log_main
    )
    if ss is not None:
        tmpres1, tmpres2 = do_one_way_anova_scipy(
            clusters, ss, name_of_variable="smearing", log_main=log_main
        )
    if norms:
        tmpres1, tmpres2 = do_one_way_anova_scipy(
            clusters, norms, name_of_variable="norm", log_main=log_main
        )
    return


def check_sorting_state(current_dir, keepchecking, log_file):
    try:
        fout = open(os.path.join(current_dir, "state.json"), "r")
        current_state = sp_utilities.convert_json_fromunicode(json.load(fout))
        fout.close()
        if current_state["done"]:
            keepchecking = 1
        else:
            keepchecking = 0
    except:
        keepchecking = 0
    return keepchecking


def read_tracker_mpi(current_dir):
    global Tracker, Blockdata
    open_tracker = 1
    if Blockdata["myid"] == Blockdata["main_node"]:
        try:
            fout = open(os.path.join(current_dir, "Tracker.json"), "r")
            Tracker = sp_utilities.convert_json_fromunicode(json.load(fout))
            fout.close()
        except:
            open_tracker = 0
    else:
        open_tracker = 0
    open_tracker = sp_utilities.bcast_number_to_all(
        open_tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    if open_tracker == 1:
        if Blockdata["myid"] != Blockdata["main_node"]:
            Tracker = 0
        Tracker = sp_utilities.wrap_mpi_bcast(
            Tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
    else:
        sp_global_def.ERROR(
            "Fail to load tracker", "%s" % current_dir, 1, Blockdata["myid"]
        )
    return


def mark_sorting_state(current_dir, sorting_done, log_file):
    # single processor job
    current_state = {}
    fout = open(os.path.join(current_dir, "state.json"), "w")
    if sorting_done:
        current_state["done"] = True
    else:
        current_state["done"] = False
    json.dump(current_state, fout)
    fout.close()
    return


def depth_clustering_box(
    work_dir,
    input_accounted_file,
    input_unaccounted_file,
    params,
    previous_params,
    nbox,
    log_main,
):
    global Tracker, Blockdata
    box_niter = 5
    no_groups_runs = 0
    if Blockdata["myid"] == Blockdata["main_node"]:
        if not os.path.exists(work_dir):
            os.mkdir(work_dir)
        freq_cutoff_dict = {}
        fout = open(os.path.join(work_dir, "freq_cutoff.json"), "w")
        json.dump(freq_cutoff_dict, fout)
        fout.close()
        log_main.add(
            "----------------------------------------------------------------------------------------------------------------"
        )
        log_main.add(
            "                        Executing pair of quasi-independent sortings, pair number %d"
            % nbox
        )
        log_main.add(
            "----------------------------------------------------------------------------------------------------------------"
        )

    ### ------- Initialization
    ncluster = 0
    nruns = 0
    keepgoing = 1
    converged = 0
    ####
    img_per_grp, number_of_groups_init, total_stack_init, minimum_grp_size_init, new_assignment = depth_box_initialization(
        work_dir, input_accounted_file, input_unaccounted_file, log_main
    )
    NUACC = total_stack_init
    NACC = 0
    ####==========--------------------------
    assignment_list = new_assignment[:]
    total_stack = total_stack_init
    current_number_of_groups = number_of_groups_init
    unaccounted_list = new_assignment[:]

    while keepgoing == 1:
        within_box_run_dir = os.path.join(work_dir, "run%d" % nruns)
        unaccounted_file = os.path.join(
            within_box_run_dir, "Unaccounted_from_previous_run.txt"
        )
        if Blockdata["myid"] == Blockdata["main_node"]:
            time_box_start = time.time()
            os.mkdir(within_box_run_dir)
            sp_utilities.write_text_file(
                unaccounted_list, unaccounted_file
            )  # new starting point
        nreassign_list = []
        if nruns > 0:
            assignment_list = create_nrandom_lists(
                unaccounted_file, current_number_of_groups, 2
            )
        if Blockdata["myid"] == Blockdata["main_node"]:
            for indep in range(2):
                sp_utilities.write_text_row(
                    assignment_list[indep],
                    os.path.join(
                        within_box_run_dir, "independent_index_%03d.txt" % indep
                    ),
                )
        run_id_file = os.path.join(within_box_run_dir, "independent_index_000.txt")
        #   iter initialization
        iter = 0
        previous_iter_ratio = 0.0
        current_iter_ratio = 0.0
        iter_dir = os.path.join(within_box_run_dir, "iter%d" % iter)
        if Blockdata["myid"] == Blockdata["main_node"]:
            os.mkdir(iter_dir)
            for indep in range(2):
                sp_utilities.write_text_row(
                    assignment_list[indep],
                    os.path.join(iter_dir, "random_assignment_%03d.txt" % indep),
                )
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

        iter_id_init_file = os.path.join(iter_dir, "random_assignment_000.txt")
        if Blockdata["myid"] == Blockdata["main_node"]:
            iter_total_stack, iter_number_of_groups = depth_iter_initialization(
                iter_id_init_file
            )
        else:
            iter_total_stack = 0
            iter_number_of_groups = 0

        iter_total_stack = sp_utilities.bcast_number_to_all(
            iter_total_stack, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        iter_number_of_groups = sp_utilities.bcast_number_to_all(
            iter_number_of_groups, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        total_stack = iter_total_stack
        current_number_of_groups = iter_number_of_groups

        #  computation starts,    prepare data
        original_data, norm_per_particle = read_data_for_sorting(
            iter_id_init_file, params, previous_params
        )

        if Tracker["nosmearing"]:
            parameterstructure = None
            paramstructure_dict = None
            paramstructure_dir = None
        else:
            paramstructure_dict = Tracker["paramstructure_dict"]
            paramstructure_dir = Tracker["paramstructure_dir"]
            parameterstructure = read_paramstructure_for_sorting(
                iter_id_init_file, paramstructure_dict, paramstructure_dir
            )
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

        Tracker["directory"] = within_box_run_dir
        if nruns == 0:  # only do it in the first box
            if Tracker["box_nxinit"] == -1:
                Tracker["nxinit"], Tracker[
                    "freq_fsc143_cutoff"
                ] = get_sorting_image_size(
                    original_data,
                    iter_id_init_file,
                    current_number_of_groups,
                    parameterstructure,
                    norm_per_particle,
                    log_main,
                )
                Tracker["box_nxinit"] = Tracker["nxinit"]
                Tracker["box_nxinit_freq"] = Tracker["freq_fsc143_cutoff"]
            else:
                Tracker["nxinit"] = Tracker["box_nxinit"]
                Tracker["freq_fsc143_cutoff"] = Tracker["box_nxinit_freq"]
        else:  # nruns>1 always estimate image size
            Tracker["nxinit"], Tracker["freq_fsc143_cutoff"] = get_sorting_image_size(
                original_data,
                iter_id_init_file,
                current_number_of_groups,
                parameterstructure,
                norm_per_particle,
                log_main,
            )

        if Blockdata["myid"] == Blockdata["main_node"]:
            log_main.add(
                "Sorting cutoff frequency is %f" % (Tracker["freq_fsc143_cutoff"])
            )

        Tracker = sp_utilities.wrap_mpi_bcast(
            Tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )

        ### preset variables in trackers
        Tracker["total_stack"] = total_stack
        Tracker["number_of_groups"] = current_number_of_groups
        minimum_grp_size = minimum_grp_size_init
        iter_previous_iter_ratio = 0.0
        iter_current_iter_ratio = 0.0
        ### prepare data for sorting: three datasets are generated cdata, fdata, srdata
        cdata, rdata, fdata = downsize_data_for_sorting(
            original_data, preshift=True, npad=1, norms=norm_per_particle
        )  # pay attentions to shifts!
        del original_data
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        srdata = precalculate_shifted_data_for_recons3D(
            rdata,
            parameterstructure,
            Tracker["refang"],
            Tracker["rshifts"],
            Tracker["delta"],
            Tracker["avgnorm"],
            Tracker["nxinit"],
            Tracker["constants"]["nnxo"],
            Tracker["nosmearing"],
            norm_per_particle,
            Tracker["constants"]["nsmear"],
        )
        del rdata
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

        while (iter < box_niter) and (converged == 0):
            for indep_run_iter in range(2):
                Tracker["directory"] = os.path.join(
                    iter_dir, "MGSKmeans_%03d" % indep_run_iter
                )
                MGSKmeans_index_file = os.path.join(
                    iter_dir, "random_assignment_%03d.txt" % indep_run_iter
                )
                if Blockdata["myid"] == Blockdata["main_node"]:
                    os.mkdir(Tracker["directory"])
                    os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

                tmp_final_list, premature = Kmeans_minimum_group_size_orien_groups(
                    cdata,
                    fdata,
                    srdata,
                    MGSKmeans_index_file,
                    params,
                    minimum_grp_size,
                    clean_volumes=True,
                )

                if Blockdata["myid"] == Blockdata["main_node"]:
                    sp_utilities.write_text_row(
                        tmp_final_list,
                        os.path.join(iter_dir, "partition_%03d.txt" % indep_run_iter),
                    )
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            if Blockdata["myid"] == Blockdata["main_node"]:
                minimum_grp_size1, maximum_grp_size1, list_of_stable, unaccounted_list, iter_current_iter_ratio, selected_number_of_groups, info_table = do_withinbox_two_way_comparison(
                    iter_dir, nbox, nruns, iter
                )  # two partitions are written in partition_dir as partition_%03d.txt
            else:
                unaccounted_list = 0
                list_of_stable = 0
                minimum_grp_size1 = 0
                Tracker = 0
            minimum_grp_size1 = sp_utilities.bcast_number_to_all(
                minimum_grp_size1, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            unaccounted_list = sp_utilities.wrap_mpi_bcast(
                unaccounted_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            list_of_stable = sp_utilities.wrap_mpi_bcast(
                list_of_stable, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            Tracker = sp_utilities.wrap_mpi_bcast(
                Tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            accounted_file = os.path.join(iter_dir, "Accounted.txt")
            unaccounted_file = os.path.join(iter_dir, "Unaccounted.txt")

            if Tracker["constants"]["do_swap_au"]:
                swap_ratio = Tracker["constants"]["swap_ratio"]
            else:
                swap_ratio = 0.0
            if Blockdata["myid"] == Blockdata["main_node"]:
                if (
                    abs(iter_current_iter_ratio - iter_previous_iter_ratio < 1.0)
                    and iter_current_iter_ratio > 90.0
                ):
                    converged = 1
            else:
                converged = 0
            converged = sp_utilities.bcast_number_to_all(
                converged, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            if (converged == 0) and (iter < box_niter):
                new_assignment_list = []
                for indep in range(2):
                    tmp_list = swap_accounted_with_unaccounted_elements_mpi(
                        accounted_file,
                        unaccounted_file,
                        log_main,
                        current_number_of_groups,
                        swap_ratio,
                    )
                    new_assignment_list.append(tmp_list)
                iter_previous_iter_ratio = iter_current_iter_ratio
                iter += 1
                iter_dir = os.path.join(within_box_run_dir, "iter%d" % iter)
                if Blockdata["myid"] == Blockdata["main_node"]:
                    os.mkdir(iter_dir)
                    for indep in range(2):
                        sp_utilities.write_text_file(
                            new_assignment_list[indep],
                            os.path.join(
                                iter_dir, "random_assignment_%03d.txt" % indep
                            ),
                        )
            else:
                break

        if Blockdata["myid"] == Blockdata["main_node"]:
            ncluster, NACC, NUACC, unaccounted_list, new_clusters = output_iter_results(
                work_dir,
                ncluster,
                NACC,
                NUACC,
                Tracker["current_img_per_grp"],
                list_of_stable,
                unaccounted_list,
                Tracker["constants"]["img_per_grp"],
                log_main,
            )
            # remove tempdir
            if os.path.exists(os.path.join(work_dir, "run%d" % nruns, "tempdir")):
                shutil.rmtree(os.path.join(work_dir, "run%d" % nruns, "tempdir"))
        else:
            ncluster = 0
            NACC = 0
            NUACC = 0
            unaccounted_list = 0
            new_clusters = 0

        new_clusters = sp_utilities.bcast_number_to_all(
            new_clusters, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        ncluster = sp_utilities.bcast_number_to_all(
            ncluster, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        NACC = sp_utilities.bcast_number_to_all(
            NACC, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        NUACC = sp_utilities.bcast_number_to_all(
            NUACC, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        unaccounted_list = sp_utilities.wrap_mpi_bcast(
            unaccounted_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )

        if new_clusters > 0:
            no_cluster = False
            no_groups_runs = 0
        else:
            no_cluster = True
            no_groups_runs += 1

        keepgoing, nruns, total_stack, current_number_of_groups = check_state_within_box_run(
            keepgoing, nruns, unaccounted_list, no_cluster
        )

        if Blockdata["myid"] == Blockdata["main_node"]:  # report current state
            if new_clusters > 0:
                log_main.add(" ")
                log_main.add(
                    "In phase %d, the program found %d groups." % (nruns, new_clusters)
                )
                for itable in range(len(info_table)):
                    log_main.add(info_table[itable])
        if no_groups_runs >= 2:
            bad_clustering = 1
            break

    if Blockdata["myid"] == Blockdata["main_node"]:
        partition = get_box_partition(work_dir, ncluster, unaccounted_list)
    else:
        partition = 0
    partition = sp_utilities.wrap_mpi_bcast(
        partition, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )

    if Blockdata["myid"] == Blockdata["main_node"]:
        if ncluster > 0:
            if ncluster == 1:
                log_main.add(
                    "In independent run  %d, the program found %d group."
                    % (nbox, ncluster)
                )
            else:
                log_main.add(
                    "In independent run  %d, the program found %d groups"
                    % (nbox, ncluster)
                )
            bad_clustering = 0
        else:
            log_main.add(
                "In independent run  %d, no reproducible groups were found " % nbox
            )
            bad_clustering = 1
            partition = [""]
        log_main.add(
            "----------------------------------------------------------------------------------------------------------------\n"
        )
        sp_utilities.write_text_row(partition, os.path.join(work_dir, "partition.txt"))
    else:
        bad_clustering = 0
    bad_clustering = sp_utilities.bcast_number_to_all(
        bad_clustering, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    return bad_clustering


def check_mpi_settings(log_main):
    global Tracker, Blockdata
    current_mpi_settings_is_bad = 0

    if Blockdata["myid"] == Blockdata["main_node"]:
        fsc_refinement = sp_utilities.read_text_file(
            os.path.join(Tracker["constants"]["masterdir"], "fsc_global.txt")
        )
        q = old_div(
            float(Tracker["constants"]["img_per_grp"]),
            float(Tracker["constants"]["total_stack"]),
        )
        for ifreq in range(len(fsc_refinement)):
            fsc_refinement[ifreq] = old_div(
                fsc_refinement[ifreq] * q, (1.0 - fsc_refinement[ifreq] * (1.0 - q))
            )
        res = 0.0
        for ifreq in range(len(fsc_refinement)):
            if fsc_refinement[ifreq] < 0.143:
                break
        res = old_div(old_div(float(ifreq), 2.0), float(len(fsc_refinement)))
        nxinit = int(2.0 * res * Tracker["constants"]["nnxo"])
        del fsc_refinement
    else:
        nxinit = 0
    nxinit = sp_utilities.bcast_number_to_all(
        nxinit, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    sys_required_mem = 1.0 * Blockdata["no_of_processes_per_group"]

    if Blockdata["myid"] == Blockdata["main_node"]:
        log_main.add("\n")
        log_main.add(
            "----------------------------------------------------------------------------------------------------------------"
        )
        log_main.add(
            "                 =======>     Number of input images and memory information     <====="
        )
        log_main.add(
            "----------------------------------------------------------------------------------------------------------------"
        )
        log_main.add(
            "Number of processes: %d number of nodes:  %d.  Number of processes per node:  %d."
            % (
                Blockdata["nproc"],
                Blockdata["no_of_groups"],
                Blockdata["no_of_processes_per_group"],
            )
        )
    try:
        image_org_size = Tracker["constants"]["nnxo"]
        image_in_core_size = nxinit
        ratio = old_div(float(nxinit), float(image_org_size))
        raw_data_size = old_div(
            float(Tracker["constants"]["total_stack"] * image_org_size * image_org_size)
            * 4.0,
            1.0e9,
        )
        raw_data_size_per_node = old_div(
            old_div(
                float(
                    Tracker["constants"]["total_stack"]
                    * image_org_size
                    * image_org_size
                )
                * 4.0,
                1.0e9,
            ),
            Blockdata["no_of_groups"],
        )
        sorting_data_size_per_node = (
            raw_data_size_per_node + 2.0 * raw_data_size_per_node * ratio ** 2
        )
        volume_size_per_node = old_div(
            (4.0 * image_in_core_size ** 3 * 8.0)
            * Blockdata["no_of_processes_per_group"],
            1.0e9,
        )
    except:
        current_mpi_settings_is_bad = 1
    if current_mpi_settings_is_bad == 1:
        sp_global_def.ERROR(
            "Initial info is not provided", "check_mpi_settings", 1, Blockdata["myid"]
        )
    try:
        mem_bytes = os.sysconf("SC_PAGE_SIZE") * os.sysconf(
            "SC_PHYS_PAGES"
        )  # e.g. 4015976448
        mem_gib = old_div(mem_bytes, (1024.0 ** 3))  # e.g. 3.74
        if Blockdata["myid"] == Blockdata["main_node"]:
            log_main.add(
                "Available memory information provided by the operating system: %5.1f GB"
                % mem_gib
            )
    except:
        mem_gib = None
    if Tracker["constants"]["memory_per_node"] == -1.0:
        if mem_gib:
            total_memory = mem_gib
        else:
            total_memory = (
                Blockdata["no_of_processes_per_group"] * 2.0
            )  # assume each CPU has 2.0 G
            if Blockdata["myid"] == Blockdata["main_node"]:
                log_main.add(
                    "Memory per node is not provided, sort3d assumes 2GB per node"
                )
        Tracker["constants"]["memory_per_node"] = total_memory
    else:
        total_memory = Tracker["constants"]["memory_per_node"]
        if Blockdata["myid"] == Blockdata["main_node"]:
            log_main.add(
                "Memory per node: %f" % Tracker["constants"]["memory_per_node"]
            )
    if Blockdata["myid"] == Blockdata["main_node"]:
        log_main.add(
            "Total number of images: %d.  Number of images per group: %d."
            % (Tracker["constants"]["total_stack"], Tracker["constants"]["img_per_grp"])
        )

    if Blockdata["myid"] == Blockdata["main_node"]:
        log_main.add("The total available memory:  %5.1f GB" % total_memory)
        log_main.add("The size of input 2D stack: %5.1f GB" % (raw_data_size))
        log_main.add(
            "The per-node amount of memory 2D data will occupy: %5.1f GB"
            % (raw_data_size_per_node)
        )

    if (
        total_memory
        - sys_required_mem
        - raw_data_size_per_node
        - volume_size_per_node
        - sorting_data_size_per_node
        - 5.0
    ) < 0.0:
        current_mpi_settings_is_bad = 1
        new_nproc = old_div(
            raw_data_size
            * (2.0 * ratio ** 2 + 1.0)
            * Blockdata["no_of_processes_per_group"],
            (total_memory - 5.0 - sys_required_mem - volume_size_per_node),
        )
        new_nproc = int(new_nproc)
        sp_global_def.ERROR(
            "Insufficient memory",
            "Suggestion: set number of processes to: %d" % new_nproc,
            1,
            Blockdata["myid"],
        )

    images_per_cpu = int(
        old_div(float(Tracker["constants"]["total_stack"]), Blockdata["nproc"])
    )
    images_per_cpu_for_unaccounted_data = old_div(
        Tracker["constants"]["img_per_grp"] * 1.5, float(Blockdata["nproc"])
    )
    if Blockdata["myid"] == Blockdata["main_node"]:
        log_main.add("Number of images per processor: %d " % images_per_cpu)
    if images_per_cpu < 5:
        sp_global_def.ERROR(
            "Number of images per processor is less than 5",
            "one may want to consider decreasing the number of processors used in MPI setting",
            0,
            Blockdata["myid"],
        )
    if Blockdata["myid"] == Blockdata["main_node"]:
        log_main.add(
            "----------------------------------------------------------------------------------------------------------------\n"
        )
    return


def get_sorting_image_size(
    original_data,
    partids,
    number_of_groups,
    sparamstructure,
    snorm_per_particle,
    log_main,
):
    global Tracker, Blockdata
    iter = 0
    Tracker["number_of_groups"] = number_of_groups
    if Blockdata["myid"] == Blockdata["main_node"]:
        log_main.add(
            "3D reconstruction is computed using window size:  %d"
            % Tracker["nxinit_refinement"]
        )
        lpartids = sp_utilities.read_text_file(partids, -1)
        if len(lpartids) == 1:
            iter_assignment = []
            for im in range(len(lpartids[0])):
                iter_assignment.append(
                    random.randint(0, number_of_groups - 1)
                )  # simple version
        else:
            iter_assignment = lpartids[0]
    else:
        iter_assignment = 0
    iter_assignment = sp_utilities.wrap_mpi_bcast(
        iter_assignment, Blockdata["main_node"]
    )

    Tracker["total_stack"] = len(iter_assignment)
    proc_list = [[None, None] for iproc in range(Blockdata["nproc"])]
    for iproc in range(Blockdata["nproc"]):
        iproc_image_start, iproc_image_end = sp_applications.MPI_start_end(
            Tracker["total_stack"], Blockdata["nproc"], iproc
        )
        proc_list[iproc] = [iproc_image_start, iproc_image_end]

    compute_noise(Tracker["nxinit_refinement"])
    rdata = downsize_data_for_rec3D(
        original_data, Tracker["nxinit_refinement"], return_real=False, npad=1
    )
    update_rdata_assignment(iter_assignment, proc_list, Blockdata["myid"], rdata)
    Tracker["nxinit"] = Tracker["nxinit_refinement"]
    compute_noise(Tracker["nxinit"])
    do3d_sorting_groups_fsc_only_iter(
        rdata, sparamstructure, snorm_per_particle, iteration=iter
    )
    del rdata

    if Blockdata["myid"] == Blockdata["main_node"]:
        fsc_data = []
        for igroup in range(Tracker["number_of_groups"]):
            for ichunk in range(2):
                tmp_fsc_data = sp_utilities.read_text_file(
                    os.path.join(
                        Tracker["directory"],
                        "fsc_driver_chunk%d_grp%03d_iter%03d.txt"
                        % (ichunk, igroup, iter),
                    ),
                    -1,
                )
                fsc_data.append(tmp_fsc_data[0])
    else:
        fsc_data = 0
    fsc_data = sp_utilities.wrap_mpi_bcast(fsc_data, Blockdata["main_node"])
    avg_fsc = [0.0 for i in range(len(fsc_data[0]))]
    avg_fsc[0] = 1.0
    for igroup in range(1):  # Use group zero first
        for ifreq in range(1, len(fsc_data[0])):
            avg_fsc[ifreq] += fsc_data[igroup][ifreq]
    fsc143 = len(fsc_data[0])
    for ifreq in range(len(avg_fsc)):
        if avg_fsc[ifreq] < 0.143:
            fsc143 = ifreq - 1
            break
    if fsc143 != 0:
        nxinit = min(
            (int(fsc143) + max(int(Tracker["constants"]["nnxo"] * 0.03), 5)) * 2,
            Tracker["constants"]["nnxo"],
        )
    else:
        sp_global_def.ERROR(
            "Program obtains wrong image size",
            "get_sorting_image_size",
            1,
            Blockdata["myid"],
        )
    freq_fsc143_cutoff = old_div(float(fsc143), float(nxinit))
    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_utilities.write_text_file(
            avg_fsc, os.path.join(Tracker["directory"], "fsc_image_size.txt")
        )
    del iter_assignment
    del proc_list
    del fsc_data
    del avg_fsc
    return nxinit, freq_fsc143_cutoff


def compute_noise(image_size):
    global Tracker, Blockdata
    if Tracker["applybckgnoise"]:  # from SPARX refinement only
        if Blockdata["myid"] == Blockdata["main_node"]:
            tsd = sp_utilities.get_im(Tracker["bckgnoise"])  # inverted power spectrum
            nnx = tsd.get_xsize()
            nny = tsd.get_ysize()
        else:
            nnx = 0
            nny = 0
        nnx = sp_utilities.bcast_number_to_all(
            nnx, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        nny = sp_utilities.bcast_number_to_all(
            nny, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        if Blockdata["myid"] != Blockdata["main_node"]:
            tsd = sp_utilities.model_blank(nnx, nny)
        sp_utilities.bcast_EMData_to_all(tsd, Blockdata["myid"], Blockdata["main_node"])
        temp_image = sp_utilities.model_blank(image_size, image_size)
        temp_image = sp_fundamentals.fft(temp_image)
        nx = temp_image.get_xsize()
        ny = temp_image.get_ysize()
        Blockdata["bckgnoise"] = []
        Blockdata["unrolldata"] = []
        for i in range(nny):
            prj = nnx * [0.0]
            for k in range(nnx):
                if tsd.get_value_at(k, i) > 0.0:
                    prj[k] = tsd.get_value_at(k, i)
            Blockdata["bckgnoise"].append(prj)
        for i in range(len(Blockdata["bckgnoise"])):
            Blockdata["unrolldata"].append(
                EMAN2_cppwrap.Util.unroll1dpw(ny, ny, Blockdata["bckgnoise"][i])
            )
    else:  # from datastack and relion
        temp_image = sp_utilities.model_blank(image_size, image_size)
        temp_image = sp_fundamentals.fft(temp_image)
        nx = temp_image.get_xsize()
        ny = temp_image.get_ysize()
        Blockdata["bckgnoise"] = [1.0] * nx
        Blockdata["unrolldata"] = EMAN2_cppwrap.Util.unroll1dpw(ny, ny, nx * [1.0])
    return


def get_params_for_analysis(orgstack, ali3d_params, smearing_file, smearing_number):
    if ali3d_params is not None:
        vecs_list = []
        norm_list = []
        ali3d_params = sp_utilities.read_text_row(ali3d_params)
        for im in range(len(ali3d_params)):
            vecs_list.append(
                sp_utilities.getvec(ali3d_params[im][0], ali3d_params[im][1])
            )
            try:
                norm_list.append(ali3d_params[im][7])
            except:
                norm_list = None
    ##
    if orgstack is not None:
        defo_list = []
        ctfs = EMAN2_cppwrap.EMUtil.get_all_attributes(orgstack, "ctf")
        for im in range(len(ctfs)):
            defo_list.append(ctfs[im].defocus)
    else:
        defo_list = None
    ##
    if smearing_file is not None:
        try:
            smearing_list = sp_utilities.read_text_file(smearing_file)
        except:
            smearing_list = None
    else:
        smearing_list = None
    ##
    return vecs_list, defo_list, smearing_list, norm_list


def do_one_way_anova_scipy(
    clusters, value_list, name_of_variable="variable", log_main=None
):
    # single cpu program
    NMAX = 30
    log_main.add(
        "================================================================================================================"
    )
    log_main.add("                                       ANOVA analysis")
    log_main.add(
        "----------------------------------------------------------------------------------------------------------------"
    )
    if len(clusters) <= 1:
        return None, None, None
    K = min(NMAX, len(clusters))
    replicas = []
    for ic in range(K):
        ll = copy.deepcopy(clusters[ic])
        ll1 = [None for i in range(len(ll))]
        for ie in range(len(ll)):
            ll1[ie] = value_list[ll[ie]]
        if name_of_variable == "defocus":
            ll1 = list(set(ll1))
        replicas.append(ll1)
    x0 = replicas[0]
    x1 = replicas[1]
    try:
        x2 = replicas[2]
    except:
        pass
    try:
        x3 = replicas[3]
    except:
        pass
    try:
        x4 = replicas[4]
    except:
        pass
    try:
        x5 = replicas[5]
    except:
        pass
    try:
        x6 = replicas[6]
    except:
        pass
    try:
        x7 = replicas[7]
    except:
        pass
    try:
        x8 = replicas[8]
    except:
        pass
    try:
        x9 = replicas[9]
    except:
        pass
    try:
        x10 = replicas[10]
    except:
        pass
    try:
        x11 = replicas[11]
    except:
        pass
    try:
        x12 = replicas[12]
    except:
        pass
    try:
        x13 = replicas[13]
    except:
        pass
    try:
        x14 = replicas[14]
    except:
        pass
    try:
        x15 = replicas[15]
    except:
        pass
    try:
        x16 = replicas[16]
    except:
        pass
    try:
        x17 = replicas[17]
    except:
        pass
    try:
        x18 = replicas[18]
    except:
        pass
    try:
        x19 = replicas[19]
    except:
        pass
    try:
        x20 = replicas[20]
    except:
        pass
    try:
        x21 = replicas[21]
    except:
        pass
    try:
        x22 = replicas[22]
    except:
        pass
    try:
        x23 = replicas[23]
    except:
        pass
    try:
        x24 = replicas[24]
    except:
        pass
    try:
        x25 = replicas[25]
    except:
        pass
    try:
        x26 = replicas[26]
    except:
        pass
    try:
        x27 = replicas[27]
    except:
        pass
    try:
        x28 = replicas[28]
    except:
        pass
    try:
        x29 = replicas[29]
    except:
        pass

    if K == 2:
        res = scipy.stats.f_oneway(x0, x1)
    elif K == 3:
        res = scipy.stats.f_oneway(x0, x1, x2)
    elif K == 4:
        res = scipy.stats.f_oneway(x0, x1, x2, x3)
    elif K == 5:
        res = scipy.stats.f_oneway(x0, x1, x2, x3, x4)
    elif K == 6:
        res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5)
    elif K == 7:
        res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6)
    elif K == 8:
        res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7)
    elif K == 9:
        res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8)
    elif K == 10:
        res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9)
    elif K == 11:
        res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
    elif K == 12:
        res = scipy.stats.f_oneway(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11)
    elif K == 13:
        res = scipy.stats.f_oneway(
            x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12
        )
    elif K == 14:
        res = scipy.stats.f_oneway(
            x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13
        )
    elif K == 15:
        res = scipy.stats.f_oneway(
            x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14
        )
    elif K == 16:
        res = scipy.stats.f_oneway(
            x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15
        )
    elif K == 17:
        res = scipy.stats.f_oneway(
            x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16
        )
    elif K == 18:
        res = scipy.stats.f_oneway(
            x0,
            x1,
            x2,
            x3,
            x4,
            x5,
            x6,
            x7,
            x8,
            x9,
            x10,
            x11,
            x12,
            x13,
            x14,
            x15,
            x16,
            x17,
        )
    elif K == 19:
        res = scipy.stats.f_oneway(
            x0,
            x1,
            x2,
            x3,
            x4,
            x5,
            x6,
            x7,
            x8,
            x9,
            x10,
            x11,
            x12,
            x13,
            x14,
            x15,
            x16,
            x17,
            x18,
        )
    elif K == 20:
        res = scipy.stats.f_oneway(
            x0,
            x1,
            x2,
            x3,
            x4,
            x5,
            x6,
            x7,
            x8,
            x9,
            x10,
            x11,
            x12,
            x13,
            x14,
            x15,
            x16,
            x17,
            x18,
            x19,
        )
    elif K == 21:
        res = scipy.stats.f_oneway(
            x0,
            x1,
            x2,
            x3,
            x4,
            x5,
            x6,
            x7,
            x8,
            x9,
            x10,
            x11,
            x12,
            x13,
            x14,
            x15,
            x16,
            x17,
            x18,
            x19,
            x20,
        )
    elif K == 22:
        res = scipy.stats.f_oneway(
            x0,
            x1,
            x2,
            x3,
            x4,
            x5,
            x6,
            x7,
            x8,
            x9,
            x10,
            x11,
            x12,
            x13,
            x14,
            x15,
            x16,
            x17,
            x18,
            x19,
            x20,
            x21,
        )
    elif K == 23:
        res = scipy.stats.f_oneway(
            x0,
            x1,
            x2,
            x3,
            x4,
            x5,
            x6,
            x7,
            x8,
            x9,
            x10,
            x11,
            x12,
            x13,
            x14,
            x15,
            x16,
            x17,
            x18,
            x19,
            x20,
            x21,
            x22,
        )
    elif K == 24:
        res = scipy.stats.f_oneway(
            x0,
            x1,
            x2,
            x3,
            x4,
            x5,
            x6,
            x7,
            x8,
            x9,
            x10,
            x11,
            x12,
            x13,
            x14,
            x15,
            x16,
            x17,
            x18,
            x19,
            x20,
            x21,
            x22,
            x23,
        )
    elif K == 25:
        res = scipy.stats.f_oneway(
            x0,
            x1,
            x2,
            x3,
            x4,
            x5,
            x6,
            x7,
            x8,
            x9,
            x10,
            x11,
            x12,
            x13,
            x14,
            x15,
            x16,
            x17,
            x18,
            x19,
            x20,
            x21,
            x22,
            x23,
            x24,
        )
    elif K == 26:
        res = scipy.stats.f_oneway(
            x0,
            x1,
            x2,
            x3,
            x4,
            x5,
            x6,
            x7,
            x8,
            x9,
            x10,
            x11,
            x12,
            x13,
            x14,
            x15,
            x16,
            x17,
            x18,
            x19,
            x20,
            x21,
            x22,
            x23,
            x24,
            x25,
        )
    elif K == 27:
        res = scipy.stats.f_oneway(
            x0,
            x1,
            x2,
            x3,
            x4,
            x5,
            x6,
            x7,
            x8,
            x9,
            x10,
            x11,
            x12,
            x13,
            x14,
            x15,
            x16,
            x17,
            x18,
            x19,
            x20,
            x21,
            x22,
            x23,
            x24,
            x25,
            x26,
        )
    elif K == 28:
        res = scipy.stats.f_oneway(
            x0,
            x1,
            x2,
            x3,
            x4,
            x5,
            x6,
            x7,
            x8,
            x9,
            x10,
            x11,
            x12,
            x13,
            x14,
            x15,
            x16,
            x17,
            x18,
            x19,
            x20,
            x21,
            x22,
            x23,
            x24,
            x25,
            x26,
            x27,
        )
    elif K == 29:
        res = scipy.stats.f_oneway(
            x0,
            x1,
            x2,
            x3,
            x4,
            x5,
            x6,
            x7,
            x8,
            x9,
            x10,
            x11,
            x12,
            x13,
            x14,
            x15,
            x16,
            x17,
            x18,
            x19,
            x20,
            x21,
            x22,
            x23,
            x24,
            x25,
            x26,
            x27,
            x28,
        )
    elif K == 30:
        res = scipy.stats.f_oneway(
            x0,
            x1,
            x2,
            x3,
            x4,
            x5,
            x6,
            x7,
            x8,
            x9,
            x10,
            x11,
            x12,
            x13,
            x14,
            x15,
            x16,
            x17,
            x18,
            x19,
            x20,
            x21,
            x22,
            x23,
            x24,
            x25,
            x26,
            x27,
            x28,
            x29,
        )
    else:
        return None, None, None
    res_table_stat = []
    for im in range(K):
        alist = sp_statistics.table_stat(replicas[im])
        res_table_stat.append(alist)
    avgs = []
    global_mean = 0.0
    for ir in range(K):
        global_mean += sum(replicas[ir])
        avgs.append(old_div(sum(replicas[ir]), float(len(replicas[ir]))))

    summed_squared_elements = 0.0
    summed_squared_elements_within_groups = [None for i in range(K)]
    std_list = []
    sst = 0.0
    nsamples = 0.0
    for i in range(K):
        nsamples += len(replicas[i])
    for i in range(K):
        ssa_per_group = 0.0
        std_per_group = 0.0
        for j in range(len(replicas[i])):
            sst += replicas[i][j] ** 2
            summed_squared_elements += replicas[i][j] * replicas[i][j]
            ssa_per_group += replicas[i][j]
            std_per_group += (replicas[i][j] - avgs[i]) ** 2
        std_list.append(numpy.sqrt(old_div(std_per_group, float(len(replicas[i])))))
        summed_squared_elements_within_groups[i] = old_div(
            ssa_per_group ** 2, float(len(replicas[i]))
        )

    sst -= old_div(global_mean ** 2, nsamples)
    ssa = sum(summed_squared_elements_within_groups) - old_div(
        global_mean ** 2, nsamples
    )
    sse = sst - ssa
    n1 = 0
    for i in range(K):
        n1 += len(replicas[i]) - 1
    msa = old_div(ssa, (K - 1.0))
    mse = old_div(sse, float(n1))
    mst = old_div(sst, float(n1))
    f_ratio = old_div(msa, mse)
    log_main.add("                              ANOVA of %s" % name_of_variable)
    log_main.add("{:5} {:^12} {:^12} ".format("ANOVA", "F-value", "Significance"))
    log_main.add("{:5} {:12.2f} {:12.2f}".format("ANOVA", res[0], res[1] * 100.0))
    log_main.add(" ")

    log_main.add(
        "ANOVA:  %s mean of all clusters: %f"
        % (name_of_variable, round(old_div(global_mean, (float(nsamples))), 4))
    )
    log_main.add("ANOVA:  Group averages")
    log_main.add(
        "{:5} {:^7} {:^8} {:^12} {:^12} ".format("ANOVA", "GID", "N", "mean", "std")
    )
    for i in range(K):
        log_main.add(
            "{:5} {:^7d} {:^8d} {:12.4f} {:12.4f}".format(
                "ANOVA",
                i,
                len(replicas[i]),
                res_table_stat[i][0],
                numpy.sqrt(res_table_stat[i][1]),
            )
        )
    log_main.add(" ")
    log_main.add("ANOVA  Pair-wise tests")
    log_main.add(
        "{:5} {:^3} {:^3} {:^12} {:^12} {:^12} {:^12}".format(
            "ANOVA", "A", "B", "avgA", "avgB", "F_value", "Significance"
        )
    )
    for ires in range(K - 1):
        for jres in range(ires + 1, K):
            cres = scipy.stats.f_oneway(replicas[ires], replicas[jres])
            log_main.add(
                "{:5} {:^3d} {:^3d} {:12.4f} {:12.4f} {:12.3f} {:12.4f} ".format(
                    "ANOVA",
                    ires,
                    jres,
                    avgs[ires],
                    avgs[jres],
                    cres[0],
                    cres[1] * 100.0,
                )
            )
    log_main.add(" ")
    log_main.add(
        "================================================================================================================\n"
    )
    return res[0], res[1]


def check_3dmask(log_main):
    global Tracker, Blockdata
    ###########################################################################
    Tracker["nxinit"] = Tracker["nxinit_refinement"]
    Tracker["currentres"] = old_div(
        float(Tracker["constants"]["fsc05"]), float(Tracker["nxinit"])
    )
    #   shrinkage, current resolution, fuse_freq
    Tracker["total_stack"] = Tracker["constants"]["total_stack"]
    Tracker["shrinkage"] = old_div(
        float(Tracker["nxinit"]), Tracker["constants"]["nnxo"]
    )
    Tracker["radius"] = Tracker["constants"]["radius"] * Tracker["shrinkage"]
    try:
        fuse_freq = Tracker["fuse_freq"]
    except:
        Tracker["fuse_freq"] = int(
            old_div(
                Tracker["constants"]["pixel_size"] * Tracker["constants"]["nnxo"],
                Tracker["constants"]["fuse_freq"],
            )
            + 0.5
        )
    Tracker = sp_utilities.wrap_mpi_bcast(
        Tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    dump_tracker(Tracker["constants"]["masterdir"])
    Tracker["shrinkage"] = old_div(
        float(Tracker["nxinit"]), Tracker["constants"]["nnxo"]
    )
    if Tracker["constants"]["focus3D"]:
        if Blockdata["myid"] == Blockdata["main_node"]:
            log_main.add(" ")
            log_main.add(
                "Sort3d is run in focus mode, and focus mask is %s"
                % Tracker["constants"]["focus3D"]
                + "\n"
            )
    return


def import_data(log_main):
    global Tracker, Blockdata
    # Two typical sorting scenarios
    # 1. import data and refinement parameters from meridien refinement;
    # 2. given data stack and xform.projection/ctf in header(For simulated test data);
    import_from_relion_refinement = 0
    import_from_sparx_refinement = 0
    import_from_data_stack = 0
    total_stack = 0
    if Tracker["constants"]["refinement_method"] == "SPARX":  # Senario one
        import_from_sparx_refinement = get_input_from_sparx_ref3d(log_main)
        Tracker["smearing"] = True
    else:  # Senario three, sorting from a given data stack, general cases
        import_from_data_stack = get_input_from_datastack(log_main)
        Tracker["constants"]["hardmask"] = True
        Tracker["applybckgnoise"] = False
        Tracker["applymask"] = True
        Tracker["smearing"] = False
    Tracker["total_stack"] = Tracker["constants"]["total_stack"]
    ###=====<----------------====>checks=====-------------
    if Tracker["constants"]["symmetry"] != Tracker["constants"]["sym"]:
        if Blockdata["myid"] == Blockdata["main_node"]:
            msg = (
                "Input symmetry %s is altered to %s after reading refinement information! "
                % (Tracker["constants"]["sym"], Tracker["constants"]["symmetry"])
            )
            log_main.add(msg)

    ## checking settings!
    if (
        old_div(
            Tracker["constants"]["total_stack"], Tracker["constants"]["img_per_grp"]
        )
        <= 1
    ):
        sp_global_def.ERROR(
            "Your img_per_grp is too large", "sxsort3d_depth.py", 1, Blockdata["myid"]
        )
    if Tracker["constants"]["minimum_grp_size"] > Tracker["constants"]["img_per_grp"]:
        sp_global_def.ERROR(
            "Minimum_grp_size is too large", "sxsort3d_depth.py", 1, Blockdata["myid"]
        )
    return


def create_masterdir():
    global Tracker, Blockdata
    masterdir = Tracker["constants"]["masterdir"]
    restart = 0
    if Blockdata["myid"] == Blockdata["main_node"]:
        if not os.path.exists(os.path.join(Tracker["constants"]["masterdir"])):
            if not masterdir:
                timestring = time.strftime("_%d_%b_%Y_%H_%M_%S", time.localtime())
                masterdir = "sort3d" + timestring
                os.makedirs(masterdir)
            else:
                if not os.path.exists(masterdir):
                    os.makedirs(masterdir)
            li = len(masterdir)
        else:
            li = len(masterdir)
            if os.path.exists(
                os.path.join(Tracker["constants"]["masterdir"], "generation_000")
            ):
                restart = 1
            if os.path.exists(
                os.path.join(Tracker["constants"]["masterdir"], "Tracker.json")
            ):
                restart = 1
        sp_global_def.write_command(masterdir)
    else:
        restart = 0
        li = 0
    restart = sp_utilities.bcast_number_to_all(
        restart, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    li = mpi.mpi_bcast(li, 1, mpi.MPI_INT, Blockdata["main_node"], mpi.MPI_COMM_WORLD)[
        0
    ]
    masterdir = sp_utilities.wrap_mpi_bcast(
        masterdir, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )

    masterdir = "".join(masterdir)
    if not Tracker["constants"]["masterdir"]:
        Tracker["constants"]["masterdir"] = masterdir
    Tracker["constants"]["chunk_0"] = os.path.join(
        Tracker["constants"]["masterdir"], "chunk_0.txt"
    )
    Tracker["constants"]["chunk_1"] = os.path.join(
        Tracker["constants"]["masterdir"], "chunk_1.txt"
    )
    return restart


def sort3d_init(to_be_decided, log_main):
    global Tracker, Blockdata
    keepsorting = 1
    if Tracker["constants"]["img_per_grp"] <= 2:
        log_main.add("Number of images per group is too small.")
        keepsorting = 0
    if Tracker["total_stack"] <= Blockdata["nproc"] * 2:
        log_main.add(
            "Either user requires too many processors, or number of images is too small."
        )
        keepsorting = 0
    return keepsorting


def print_shell_command(args_list, log_main):
    global Tracker, Blockdata
    if Blockdata["myid"] == Blockdata["main_node"]:
        log_main.add("The shell line command:")
        line = ""
        for a in args_list:
            line += a + " "
        log_main.add(line)
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    return


def AI_MGSKmeans(
    iter_assignment,
    last_iter_assignment,
    best_assignment,
    keepgoing,
    best_score,
    stopercnt,
    minimum_grp_size,
    log_file,
):
    # single cpu function
    group_dict = shutil.collections.Counter(iter_assignment)
    number_of_groups = len(group_dict)
    msg = "group id     group size"
    log_file.add(msg)
    clusters = []
    for igrp in range(len(group_dict)):
        msg = "%5d    %10d" % (igrp, group_dict[igrp])
        log_file.add(msg)
        clusters.append(group_dict[igrp])
    is_unicorn_cluster = 0
    nc = 0
    tot = 0
    for ic in range(len(clusters)):
        tot += clusters[ic]
        if clusters[ic] < minimum_grp_size + len(clusters):
            nc += 1
    if old_div(tot, minimum_grp_size) > 2 * len(clusters) and nc + 1 == len(clusters):
        is_unicorn_cluster = 1
    if is_unicorn_cluster == 0:
        sum_newindices1 = 0
        sum_newindices2 = 0
        ratio, newindices, stable_clusters = compare_two_iterations(
            iter_assignment, last_iter_assignment, number_of_groups
        )
        for idx in range(len(newindices)):
            sum_newindices1 += newindices[idx][0]
            sum_newindices2 += newindices[idx][1]
            if newindices[idx][0] != newindices[idx][1]:
                msg = "Group %d  swaps with group %d " % (
                    newindices[idx][0],
                    newindices[idx][1],
                )
                log_file.add(msg)
        changed_nptls = 100.0 - ratio * 100.0
        if best_score >= changed_nptls:
            best_score = changed_nptls
            best_assignment = copy.copy(iter_assignment)
        if changed_nptls < stopercnt:
            keepgoing = 0
    else:
        msg = "Unicorn cluster is found. shuffle assignment."
        log_file.add(msg)

        iter_assignment = shuffle_assignment(iter_assignment, number_of_groups)
        best_score = 100.0
        changed_nptls = 100.0
        keepgoing = 1
        best_assignment = copy.copy(iter_assignment)
    return best_score, changed_nptls, keepgoing, best_assignment, iter_assignment


def shuffle_assignment(iter_assignment, number_of_groups):
    new_assignment = list(range(0, len(iter_assignment)))
    tmp_assignment = list(range(0, len(iter_assignment)))
    ngroup = 0
    while len(tmp_assignment) >= 1:
        random.shuffle(tmp_assignment)
        im = tmp_assignment[0]
        new_assignment[im] = ngroup % number_of_groups
        del tmp_assignment[0]
        ngroup += 1
    return new_assignment


#####


def Kmeans_minimum_group_size_orien_groups(
    cdata, fdata, srdata, partids, params, minimum_group_size_init, clean_volumes=False
):
    global Tracker, Blockdata

    # ==========---------- >>>>EQKmeans initialization ==========------------
    log_main = sp_logger.Logger()
    log_main = sp_logger.Logger(sp_logger.BaseLogger_Files())
    log_main.prefix = Tracker["directory"] + "/"
    premature = 0
    changed_nptls = 100.0
    number_of_groups = Tracker["number_of_groups"]
    stopercnt = Tracker["constants"]["stop_mgskmeans_percentage"]
    total_iter = 0
    require_check_setting = False
    partial_rec3d = False
    best_score = 100.0
    best_assignment = []
    max_iter = Tracker["total_number_of_iterations"]
    last_score = 100.0
    fixed_value = 100.0
    has_converged = 0
    times_around_fixed_value = 0
    ###=====<------------

    if Blockdata["myid"] == Blockdata["main_node"]:
        try:
            if os.path.exists(Tracker["constants"]["mask3D"]):  # prepare mask
                mask3D = sp_utilities.get_im(Tracker["constants"]["mask3D"])
                if mask3D.get_xsize() != Tracker["nxinit"]:
                    mask3D = sp_fundamentals.fdecimate(
                        mask3D,
                        Tracker["nxinit"],
                        Tracker["nxinit"],
                        Tracker["nxinit"],
                        True,
                        False,
                    )
        except:
            mask3D = sp_utilities.model_circle(
                Tracker["constants"]["radius"],
                Tracker["constants"]["nnxo"],
                Tracker["constants"]["nnxo"],
                Tracker["constants"]["nnxo"],
            )
            mask3D = sp_fundamentals.fdecimate(
                mask3D,
                Tracker["nxinit"],
                Tracker["nxinit"],
                Tracker["nxinit"],
                True,
                False,
            )
    else:
        mask3D = sp_utilities.model_blank(
            Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"]
        )
    sp_utilities.bcast_EMData_to_all(mask3D, Blockdata["myid"], Blockdata["main_node"])

    if Blockdata["myid"] == Blockdata["main_node"]:
        lpartids = sp_utilities.read_text_file(partids, -1)
        if len(lpartids) == 1:  # Not true for sorting
            iter_assignment = []
            for im in range(len(lpartids[0])):
                iter_assignment.append(
                    random.randint(0, number_of_groups - 1)
                )  # simple version
        else:
            iter_assignment = lpartids[0]
    else:
        iter_assignment = 0
    iter_assignment = sp_utilities.wrap_mpi_bcast(
        iter_assignment, Blockdata["main_node"]
    )  # initial assignment
    total_stack = len(iter_assignment)
    Tracker["total_stack"] = total_stack
    nima = len(cdata)
    image_start, image_end = sp_applications.MPI_start_end(
        Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"]
    )

    Tracker["min_orien_group_size"] = (
        Tracker["number_of_groups"] * Tracker["minimum_ptl_number"]
    )
    angle_step = get_angle_step_from_number_of_orien_groups(
        Tracker["constants"]["orientation_groups"]
    )
    ptls_in_orien_groups = get_angle_step_and_orien_groups_mpi(
        params, partids, angle_step
    )
    minimum_group_size = max(
        minimum_group_size_init, Tracker["number_of_groups"] * len(ptls_in_orien_groups)
    )
    minimum_group_size_ratio = min(
        old_div(
            (minimum_group_size * Tracker["number_of_groups"]),
            float(Tracker["total_stack"]),
        ),
        0.95,
    )

    ### printed info
    if Blockdata["myid"] == Blockdata["main_node"]:
        log_main.add(
            "----------------------------------------------------------------------------------------------------------------"
        )
        log_main.add(" ==========> MGSKmeans clustering <========== ")
        log_main.add(
            "----------------------------------------------------------------------------------------------------------------"
        )
        msg = (
            "Total_stack:  %d K = : %d  nxinit: %d  CTF:  %s  Symmetry:  %s  stop percentage: %f  3-D mask: %s focus mask: %s  Comparison method: %s  minimum_group_size: %d orien  %d"
            % (
                Tracker["total_stack"],
                Tracker["number_of_groups"],
                Tracker["nxinit"],
                Tracker["constants"]["CTF"],
                Tracker["constants"]["symmetry"],
                stopercnt,
                Tracker["constants"]["mask3D"],
                Tracker["constants"]["focus3D"],
                Tracker["constants"]["comparison_method"],
                minimum_group_size,
                len(ptls_in_orien_groups),
            )
        )
        log_main.add(msg)

    proc_list = [[None, None] for iproc in range(Blockdata["nproc"])]

    for iproc in range(Blockdata["nproc"]):
        iproc_image_start, iproc_image_end = sp_applications.MPI_start_end(
            Tracker["total_stack"], Blockdata["nproc"], iproc
        )
        proc_list[iproc] = [iproc_image_start, iproc_image_end]

    compute_noise(Tracker["nxinit"])
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    # cdata, rdata, fdata = downsize_data_for_sorting(original_data, preshift = True, npad = 1, norms =norm_per_particle)# pay attentions to shifts!
    # mpi_barrier(MPI_COMM_WORLD)

    # srdata = precalculate_shifted_data_for_recons3D(rdata, paramstructure, Tracker["refang"], \
    #   Tracker["rshifts"], Tracker["delta"], Tracker["avgnorm"], Tracker["nxinit"], \
    #     Tracker["constants"]["nnxo"], Tracker["nosmearing"], norm_per_particle, Tracker["constants"]["nsmear"])
    # del rdata
    # mpi_barrier(MPI_COMM_WORLD)

    last_iter_assignment = copy.copy(iter_assignment)
    best_assignment = copy.copy(iter_assignment)
    total_iter = 0
    keepgoing = 1
    do_partial_rec3d = 0
    partial_rec3d = False

    while total_iter < max_iter:
        ptls_in_orien_groups = get_angle_step_and_orien_groups_mpi(
            params, partids, angle_step
        )
        if Blockdata["myid"] == Blockdata["main_node"]:
            msg = "Iteration %d particle assignment changed ratio  %f " % (
                total_iter,
                changed_nptls,
            )
            log_main.add(msg)
            sp_utilities.write_text_file(
                iter_assignment,
                os.path.join(Tracker["directory"], "assignment%03d.txt" % total_iter),
            )
            if changed_nptls < 50.0:
                do_partial_rec3d = 1
            else:
                do_partial_rec3d = 0
        else:
            do_partial_rec3d = 0
        do_partial_rec3d = sp_utilities.bcast_number_to_all(
            do_partial_rec3d, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        if do_partial_rec3d == 1:
            partial_rec3d = True
        else:
            partial_rec3d = False

        update_data_assignment(
            cdata,
            srdata,
            iter_assignment,
            proc_list,
            Tracker["nosmearing"],
            Blockdata["myid"],
        )
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        do3d_sorting_groups_nofsc_smearing_iter(
            srdata, partial_rec3d, iteration=total_iter
        )

        local_peaks = [0.0 for im in range(number_of_groups * nima)]
        total_im = 0
        local_kmeans_peaks = [-1.0e23 for im in range(nima)]
        ## compute peaks and save them in 1D list
        for iref in range(number_of_groups):
            if Blockdata["myid"] == Blockdata["last_node"]:
                try:
                    fsc143 = Tracker["fsc143"][iref]
                except:
                    fsc143 = 0.0
                try:
                    fsc05 = Tracker["fsc05"][iref]
                except:
                    fsc05 = 0.0
                ref_vol = sp_utilities.get_im(
                    os.path.join(
                        Tracker["directory"],
                        "vol_grp%03d_iter%03d.hdf" % (iref, total_iter),
                    )
                )
                nnn = ref_vol.get_xsize()
                if Tracker["nxinit"] != nnn:
                    ref_vol = sp_fundamentals.fdecimate(
                        ref_vol,
                        Tracker["nxinit"],
                        Tracker["nxinit"],
                        Tracker["nxinit"],
                        True,
                        False,
                    )
                stat = EMAN2_cppwrap.Util.infomask(ref_vol, mask3D, False)
                ref_vol -= stat[0]
                if stat[1] != 0.0:
                    EMAN2_cppwrap.Util.mul_scalar(ref_vol, old_div(1.0, stat[1]))
                ref_vol *= mask3D
            else:
                ref_vol = sp_utilities.model_blank(
                    Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"]
                )
            sp_utilities.bcast_EMData_to_all(
                ref_vol, Blockdata["myid"], Blockdata["last_node"]
            )
            ## Image comparison optimal solution is the larger one
            if Tracker["constants"]["comparison_method"] == "cross":
                ref_peaks = compare_two_images_cross(cdata, ref_vol)
            else:
                ref_peaks = compare_two_images_eucd(cdata, ref_vol, fdata)
            for im in range(nima):
                local_peaks[total_im] = ref_peaks[im]
                total_im += 1
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        del ref_vol
        # pass to main_node
        if Blockdata["myid"] == Blockdata["main_node"]:
            dmatrix = [
                [0.0 for im in range(Tracker["total_stack"])]
                for iref in range(number_of_groups)
            ]
            for im in range(len(local_peaks)):
                dmatrix[old_div(im, nima)][im % nima + image_start] = local_peaks[im]
        else:
            dmatrix = 0
        if Blockdata["myid"] != Blockdata["main_node"]:
            sp_utilities.wrap_mpi_send(
                local_peaks, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
        else:
            for iproc in range(Blockdata["nproc"]):
                if iproc != Blockdata["main_node"]:
                    local_peaks = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
                    iproc_nima = proc_list[iproc][1] - proc_list[iproc][0]
                    for im in range(len(local_peaks)):
                        dmatrix[old_div(im, iproc_nima)][
                            im % iproc_nima + proc_list[iproc][0]
                        ] = local_peaks[im]
        dmatrix = sp_utilities.wrap_mpi_bcast(
            dmatrix, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        last_iter_assignment = copy.copy(iter_assignment)
        iter_assignment = [-1 for iptl in range(Tracker["total_stack"])]
        for iorien in range(len(ptls_in_orien_groups)):
            if iorien % Blockdata["nproc"] == Blockdata["myid"]:
                local_assignment = do_assignment_by_dmatrix_orien_group_minimum_group_size(
                    dmatrix,
                    ptls_in_orien_groups[iorien],
                    Tracker["number_of_groups"],
                    minimum_group_size_ratio,
                )
                for iptl in range(len(ptls_in_orien_groups[iorien])):
                    iter_assignment[
                        ptls_in_orien_groups[iorien][iptl]
                    ] = local_assignment[iptl]
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        if Blockdata["myid"] != Blockdata["main_node"]:
            sp_utilities.wrap_mpi_send(
                iter_assignment, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
        else:
            for iproc in range(Blockdata["nproc"]):
                if iproc != Blockdata["main_node"]:
                    dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
                    for iptl in range(len(dummy)):
                        if dummy[iptl] != -1:
                            iter_assignment[iptl] = dummy[iptl]
                        else:
                            pass
        iter_assignment = sp_utilities.wrap_mpi_bcast(
            iter_assignment, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        if Blockdata["myid"] == Blockdata["main_node"]:
            last_score = changed_nptls
            best_score, changed_nptls, keepgoing, best_assignmen, iter_assignment = AI_MGSKmeans(
                iter_assignment,
                last_iter_assignment,
                best_assignment,
                keepgoing,
                best_score,
                stopercnt,
                minimum_group_size,
                log_main,
            )
            if abs(last_score - changed_nptls) < 1.0:
                if times_around_fixed_value == 0:
                    fixed_value = changed_nptls
                    times_around_fixed_value += 1
                else:
                    if abs(changed_nptls - fixed_value) < 1.0:
                        times_around_fixed_value += 1
                    else:
                        times_around_fixed_value = 0
                        fixed_value = changed_nptls
        else:
            iter_assignment = 0
            best_assignment = 0
            keepgoing = 1
        iter_assignment = sp_utilities.wrap_mpi_bcast(
            iter_assignment, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        best_assignment = sp_utilities.wrap_mpi_bcast(
            best_assignment, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        times_around_fixed_value = sp_utilities.bcast_number_to_all(
            times_around_fixed_value, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        keepgoing = sp_utilities.bcast_number_to_all(
            keepgoing, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        total_iter += 1
        last_iter_assignment = copy.copy(iter_assignment)
        if times_around_fixed_value >= 3:
            keepgoing = 0
        if keepgoing == 0:
            break
    # Finalize
    update_data_assignment(
        cdata,
        srdata,
        iter_assignment,
        proc_list,
        Tracker["nosmearing"],
        Blockdata["myid"],
    )
    res_sort3d = get_sorting_all_params(cdata)
    del cdata
    del srdata
    del iter_assignment
    del last_iter_assignment
    del best_assignment
    if mask3D:
        del mask3D
    if Blockdata["myid"] == Blockdata["main_node"]:
        if best_score > Tracker["constants"]["stop_mgskmeans_percentage"]:
            premature = 1
        partition, ali3d_params_list = parsing_sorting_params(partids, res_sort3d)
        sp_utilities.write_text_row(
            partition, os.path.join(Tracker["directory"], "list.txt")
        )
        shutil.rmtree(os.path.join(Tracker["directory"], "tempdir"))
    else:
        partition = 0
    partition = sp_utilities.wrap_mpi_bcast(
        partition, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    premature = sp_utilities.bcast_number_to_all(
        premature, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    if Blockdata["myid"] == Blockdata["last_node"]:
        if clean_volumes:
            for jter in range(total_iter):
                for igroup in range(Tracker["number_of_groups"]):
                    os.remove(
                        os.path.join(
                            Tracker["directory"],
                            "vol_grp%03d_iter%03d.hdf" % (igroup, jter),
                        )
                    )
    if require_check_setting:
        if Blockdata["myid"] == Blockdata["main_node"]:
            print("Warning: the sorting settings, such as img_per_grp requires a check")
    return partition, premature


def do_assignment_by_dmatrix_orien_group_minimum_group_size(
    dmatrix, orien_group_members, number_of_groups, minimum_group_size_ratio
):
    results = [[] for i in range(number_of_groups)]
    nima = len(orien_group_members)
    minimum_group_size = int(old_div(minimum_group_size_ratio * nima, number_of_groups))
    submatrix = numpy.zeros((number_of_groups, nima))
    for i in range(number_of_groups):
        for j in range(len(orien_group_members)):
            submatrix[i][j] = dmatrix[i][orien_group_members[j]] * (
                -1.0
            )  # sort in descending order
    tmp_array = numpy.argsort(submatrix, axis=1)
    rmatrix = []
    for i in range(number_of_groups):
        rmatrix.append(tmp_array[i].tolist())
    del tmp_array
    while len(rmatrix[0]) > nima - minimum_group_size * number_of_groups:
        tarray = []
        for i in range(number_of_groups):
            tarray.append(rmatrix[i][0])
        value_list, index_list = numpy.unique(numpy.array(tarray), return_index=True)
        duplicate_list = (
            numpy.setdiff1d(numpy.arange(number_of_groups), index_list)
        ).tolist()
        index_list = index_list.tolist()
        value_list = value_list.tolist()
        if len(value_list) < number_of_groups:
            for i in range(len(index_list)):
                if tarray[index_list[i]] == tarray[duplicate_list[0]]:
                    duplicate_list.append(index_list[i])  # find all duplicated ones
            random.shuffle(duplicate_list)
            duplicate_list.remove(duplicate_list[0])
            for i in range(len(duplicate_list)):  # swap the first row with the next row
                index_column = 1
                while rmatrix[duplicate_list[i]][index_column] in value_list:
                    index_column += 1  # search along column non-equal ones
                value_list.append(rmatrix[duplicate_list[i]][index_column])
                rmatrix[duplicate_list[i]][0], rmatrix[duplicate_list[i]][
                    index_column
                ] = (
                    rmatrix[duplicate_list[i]][index_column],
                    rmatrix[duplicate_list[i]][0],
                )
        for i in range(number_of_groups):
            results[i].append(rmatrix[i][0])
            for j in range(number_of_groups):
                rmatrix[i].remove(value_list[j])  # remove K elements from each column
    kmeans_ptl_list = (
        numpy.delete(numpy.array(list(range(nima))), numpy.array(results).ravel())
    ).tolist()  # ravel works only for even size
    del rmatrix
    for iptl in range(len(kmeans_ptl_list)):
        max_indexes = numpy.argwhere(
            submatrix[:, kmeans_ptl_list[iptl]]
            <= submatrix[:, kmeans_ptl_list[iptl]][
                submatrix[:, kmeans_ptl_list[iptl]].argmin()
            ]
        )
        if len(max_indexes) > 1:
            t = list(range(len(max_indexes)))
            random.shuffle(t)
            results[max_indexes[t[0]][0]].append(kmeans_ptl_list[iptl])
        else:
            results[max_indexes[0][0]].append(kmeans_ptl_list[iptl])
    iter_assignment = [-1 for i in range(nima)]
    for i in range(number_of_groups):
        results[i].sort()
        for j in range(len(results[i])):
            iter_assignment[results[i][j]] = i
    del results
    del submatrix
    return iter_assignment


### various reading data
### 1


def get_shrink_data_sorting(
    partids, partstack, return_real=False, preshift=True, apply_mask=True, npad=1
):
    # The function will read from stack a subset of images specified in partids
    #   and assign to them parameters from partstack with optional CTF application and shifting of the data.
    # So, the lengths of partids and partstack are the same.
    #  The read data is properly distributed among MPI threads.
    # 10142015 --- preshift is set to True when doing 3-D sorting.
    # chunk_id are set when data is read in
    global Tracker, Blockdata
    mask2D = sp_utilities.model_circle(
        Tracker["constants"]["radius"],
        Tracker["constants"]["nnxo"],
        Tracker["constants"]["nnxo"],
    )
    shrinkage = old_div(Tracker["nxinit"], float(Tracker["constants"]["nnxo"]))
    radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
    if Blockdata["myid"] == Blockdata["main_node"]:
        lpartids = sp_utilities.read_text_file(partids, -1)
        if len(lpartids) == 1:
            lpartids = lpartids[0]
            groupids = len(lpartids) * [-1]
        else:
            groupids = lpartids[0]
            lpartids = lpartids[1]
    else:
        lpartids = 0
        groupids = 0
    lpartids = sp_utilities.wrap_mpi_bcast(lpartids, Blockdata["main_node"])
    groupids = sp_utilities.wrap_mpi_bcast(groupids, Blockdata["main_node"])
    Tracker["total_stack"] = len(lpartids)
    if Blockdata["myid"] == Blockdata["main_node"]:
        partstack = sp_utilities.read_text_row(partstack)
    else:
        partstack = 0
    partstack = sp_utilities.wrap_mpi_bcast(partstack, Blockdata["main_node"])

    if Tracker["total_stack"] < Blockdata["nproc"]:
        sp_global_def.ERROR(
            "Wrong MPI settings!", "get_shrink_data_sorting", 1, Blockdata["myid"]
        )
    else:
        image_start, image_end = sp_applications.MPI_start_end(
            Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"]
        )
    lpartids = lpartids[image_start:image_end]
    groupids = groupids[image_start:image_end]
    nima = image_end - image_start
    data = [None] * nima
    for im in range(nima):
        data[im] = sp_utilities.get_im(Tracker["constants"]["orgstack"], lpartids[im])
        try:
            phi, theta, psi, sx, sy, chunk_id, particle_group_id = (
                partstack[lpartids[im]][0],
                partstack[lpartids[im]][1],
                partstack[lpartids[im]][2],
                partstack[lpartids[im]][3],
                partstack[lpartids[im]][4],
                partstack[lpartids[im]][5],
                partstack[lpartids[im]][6],
            )
        except:
            phi, theta, psi, sx, sy, chunk_id, particle_group_id = (
                partstack[lpartids[im]][0],
                partstack[lpartids[im]][1],
                partstack[lpartids[im]][2],
                partstack[lpartids[im]][3],
                partstack[lpartids[im]][4],
                partstack[lpartids[im]][5],
                -1,
            )
        if preshift:  # always true
            data[im] = sp_fundamentals.fshift(data[im], sx, sy)
            sx = 0.0
            sy = 0.0
        st = EMAN2_cppwrap.Util.infomask(data[im], mask2D, False)
        data[im] -= st[0]
        data[im] = old_div(data[im], st[1])
        if apply_mask:
            data[im] = sp_morphology.cosinemask(
                data[im], radius=Tracker["constants"]["radius"]
            )
        # FT
        data[im] = sp_fundamentals.fft(data[im])
        nny = data[im].get_ysize()
        if Tracker["constants"]["CTF"]:
            ctf_params = data[im].get_attr("ctf")
            data[im] = sp_fundamentals.fdecimate(
                data[im],
                Tracker["nxinit"] * npad,
                Tracker["nxinit"] * npad,
                1,
                False,
                False,
            )
            ctf_params.apix = old_div(ctf_params.apix, shrinkage)
            data[im].set_attr("ctf", ctf_params)
            data[im].set_attr("ctf_applied", 0)
            if return_real:
                data[im] = sp_fundamentals.fft(data[im])
        else:
            ctf_params = data[im].get_attr_default("ctf", False)
            if ctf_params:
                ctf_params.apix = old_div(ctf_params.apix, shrinkage)
                data[im].set_attr("ctf", ctf_params)
                data[im].set_attr("ctf_applied", 0)
            data[im] = sp_fundamentals.fdecimate(
                data[im], nxinit * npad, nxinit * npad, 1, True, False
            )
            apix = Tracker["constants"]["pixel_size"]
            data[im].set_attr("apix", old_div(apix, shrinkage))
        if not return_real:
            data[im].set_attr("padffted", 1)
        data[im].set_attr("npad", npad)
        sp_utilities.set_params_proj(data[im], [phi, theta, psi, 0.0, 0.0])
        data[im].set_attr("chunk_id", chunk_id)
        data[im].set_attr("group", groupids[im])
        data[im].set_attr("particle_group", particle_group_id)
        if Tracker["applybckgnoise"]:
            data[im].set_attr("bckgnoise", Blockdata["bckgnoise"][particle_group_id])
            data[im].set_attr(
                "qt", float(Tracker["constants"]["nnxo"] * Tracker["constants"]["nnxo"])
            )
        else:
            data[im].set_attr("bckgnoise", Blockdata["bckgnoise"])  # constant list
    return data


###2


def read_data_for_sorting(partids, partstack, previous_partstack):
    # The function will read from stack a subset of images specified in partids
    # and assign to them parameters from partstack with optional CTF application and shifting of the data.
    # So, the lengths of partids and partstack are the same.
    global Tracker, Blockdata
    # functions:
    # read in data
    if Blockdata["myid"] == Blockdata["main_node"]:
        lpartids = sp_utilities.read_text_file(partids, -1)
        if len(lpartids) == 1:
            lpartids = lpartids[0]
            groupids = len(lpartids) * [-1]
        else:
            groupids = lpartids[0]
            lpartids = lpartids[1]
    else:
        lpartids = 0
        groupids = 0
    lpartids = sp_utilities.wrap_mpi_bcast(lpartids, Blockdata["main_node"])
    groupids = sp_utilities.wrap_mpi_bcast(groupids, Blockdata["main_node"])
    Tracker["total_stack"] = len(lpartids)
    if Blockdata["myid"] == Blockdata["main_node"]:
        partstack = sp_utilities.read_text_row(partstack)
    else:
        partstack = 0
    partstack = sp_utilities.wrap_mpi_bcast(partstack, Blockdata["main_node"])
    if Blockdata["myid"] == Blockdata["main_node"]:
        previous_partstack = sp_utilities.read_text_row(previous_partstack)
    else:
        previous_partstack = 0
    previous_partstack = sp_utilities.wrap_mpi_bcast(
        previous_partstack, Blockdata["main_node"]
    )
    if Tracker["total_stack"] < Blockdata["nproc"]:
        sp_global_def.ERROR(
            "Number of processors in use is larger than the total number of images",
            "get_data_and_prep",
            1,
            Blockdata["myid"],
        )
    else:
        image_start, image_end = sp_applications.MPI_start_end(
            Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"]
        )
    lpartids = lpartids[image_start:image_end]
    groupids = groupids[image_start:image_end]
    nima = image_end - image_start
    data = [None] * nima
    norm_per_particle = [None for im in range(nima)]
    for im in range(nima):
        image = sp_utilities.get_im(Tracker["constants"]["orgstack"], lpartids[im])
        try:
            phi, theta, psi, sx, sy, chunk_id, particle_group_id, mnorm = (
                partstack[lpartids[im]][0],
                partstack[lpartids[im]][1],
                partstack[lpartids[im]][2],
                partstack[lpartids[im]][3],
                partstack[lpartids[im]][4],
                partstack[lpartids[im]][5],
                partstack[lpartids[im]][6],
                partstack[lpartids[im]][7],
            )
        except:
            phi, theta, psi, sx, sy, chunk_id, particle_group_id, mnorm = (
                partstack[lpartids[im]][0],
                partstack[lpartids[im]][1],
                partstack[lpartids[im]][2],
                partstack[lpartids[im]][3],
                partstack[lpartids[im]][4],
                partstack[lpartids[im]][5],
                -1,
                1.0,
            )
        sx1, sy1 = (
            previous_partstack[lpartids[im]][3],
            previous_partstack[lpartids[im]][4],
        )
        sp_utilities.set_params_proj(image, [phi, theta, psi, 0.0, 0.0])
        image.set_attr("chunk_id", chunk_id)
        image.set_attr("group", groupids[im])
        image.set_attr("particle_group", particle_group_id)
        image.set_attr("previous_shifts", [sx1, sy1])
        image.set_attr("current_shifts", [sx, sy])
        norm_per_particle[im] = mnorm
        data[im] = image
    return data, norm_per_particle


###6 read paramstructure


def read_paramstructure_for_sorting(
    partids, paramstructure_dict_file, paramstructure_dir
):
    global Tracker, Blockdata
    if Blockdata["myid"] == Blockdata["main_node"]:
        lcore = sp_utilities.read_text_file(partids, -1)
    else:
        lcore = 0
    lcore = sp_utilities.wrap_mpi_bcast(lcore, Blockdata["main_node"])
    if len(lcore) == 1:
        lcore = lcore[0]
    else:
        lcore = lcore[1]
    psize = len(lcore)
    oldparamstructure = []
    im_start, im_end = sp_applications.MPI_start_end(
        psize, Blockdata["nproc"], Blockdata["myid"]
    )
    lcore = lcore[im_start:im_end]
    nima = len(lcore)
    if Blockdata["myid"] == Blockdata["main_node"]:
        tmp_list = sp_utilities.read_text_row(paramstructure_dict_file)
    else:
        tmp_list = 0
    tmp_list = sp_utilities.wrap_mpi_bcast(tmp_list, Blockdata["main_node"])
    pdict = {}
    for im in range(len(lcore)):
        pdict[im] = tmp_list[lcore[im]]
    oldparamstructure = []
    nptl = 0
    last_old_paramstructure_file = None
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    for iproc in range(Blockdata["nproc"]):
        if Blockdata["myid"] == iproc:  # always read oldparamstructure sequentially
            while nptl < nima:
                [
                    jason_of_cpu_id,
                    chunk_id,
                    iteration,
                    ptl_id_on_cpu,
                    global_index,
                ] = pdict[nptl]
                old_paramstructure_file = os.path.join(
                    paramstructure_dir,
                    "oldparamstructure_%d_%03d_%03d.json"
                    % (chunk_id, jason_of_cpu_id, iteration),
                )
                if old_paramstructure_file != last_old_paramstructure_file:
                    fout = open(old_paramstructure_file, "r")
                    paramstructure = sp_utilities.convert_json_fromunicode(
                        json.load(fout)
                    )
                    fout.close()
                last_old_paramstructure_file = old_paramstructure_file
                oldparamstructure.append(paramstructure[ptl_id_on_cpu])
                nptl += 1
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    return oldparamstructure


###7 copy oldparamstructures from meridien


def copy_oldparamstructure_from_meridien_MPI(selected_iteration, log_main):
    global Tracker, Blockdata
    Tracker["directory"] = os.path.join(
        Tracker["constants"]["masterdir"], "main%03d" % selected_iteration
    )
    Tracker["paramstructure_dir"] = os.path.join(
        Tracker["directory"], "oldparamstructure"
    )
    old_refinement_iter_directory = os.path.join(
        Tracker["constants"]["refinement_dir"], "main%03d" % selected_iteration
    )
    old_refinement_previous_iter_directory = os.path.join(
        Tracker["constants"]["refinement_dir"], "main%03d" % (selected_iteration - 1)
    )
    if Blockdata["myid"] == Blockdata["main_node"]:
        if not os.path.exists(Tracker["paramstructure_dir"]):
            os.mkdir(
                os.path.join(
                    Tracker["constants"]["masterdir"], "main%03d" % selected_iteration
                )
            )
            os.mkdir(Tracker["paramstructure_dir"])
    Tracker["refang"] = sp_utilities.read_text_row(
        os.path.join(old_refinement_iter_directory, "refang.txt")
    )
    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_utilities.write_text_row(
            Tracker["refang"], os.path.join(Tracker["directory"], "refang.txt")
        )
    Tracker["rshifts"] = sp_utilities.read_text_row(
        os.path.join(old_refinement_iter_directory, "rshifts.txt")
    )
    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_utilities.write_text_row(
            Tracker["refang"], os.path.join(Tracker["directory"], "rshifts.txt")
        )
    my_last_params = sp_utilities.read_text_file(
        os.path.join(
            old_refinement_previous_iter_directory,
            "params_%03d.txt" % (selected_iteration - 1),
        ),
        -1,
    )
    my_parstack = sp_utilities.read_text_file(
        os.path.join(Tracker["constants"]["masterdir"], "refinement_parameters.txt"), -1
    )
    if Blockdata["myid"] == Blockdata["main_node"]:
        my_parstack[3:5] = my_last_params[3:5]
        sp_utilities.write_text_file(
            my_parstack,
            os.path.join(
                Tracker["constants"]["masterdir"], "previous_refinement_parameters.txt"
            ),
        )
    Tracker["previous_parstack"] = os.path.join(
        Tracker["constants"]["masterdir"], "previous_refinement_parameters.txt"
    )
    nproc_previous = 0
    procid = 0
    old_refinement_iter_dir = os.path.join(
        Tracker["constants"]["refinement_dir"], "main%03d" % selected_iteration
    )
    if Blockdata["myid"] == Blockdata["main_node"]:
        while os.path.exists(
            os.path.join(
                old_refinement_iter_dir,
                "oldparamstructure",
                "oldparamstructure_%01d_%03d_%03d.json"
                % (procid, nproc_previous, selected_iteration),
            )
        ):
            nproc_previous += 1
    nproc_previous = sp_utilities.bcast_number_to_all(
        nproc_previous, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    Blockdata["nproc_previous"] = nproc_previous
    oldparamstructure = [[], []]
    local_dict = {}
    for procid in range(2):
        smearing_list = []
        if Blockdata["myid"] == Blockdata["main_node"]:
            lcore = sp_utilities.read_text_file(
                os.path.join(Tracker["constants"]["masterdir"], "chunk_%d.txt" % procid)
            )
        else:
            lcore = 0
        lcore = sp_utilities.wrap_mpi_bcast(
            lcore, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        psize = len(lcore)
        oldparamstructure[procid] = []
        im_start, im_end = sp_applications.MPI_start_end(
            psize, Blockdata["nproc"], Blockdata["myid"]
        )
        local_lcore = lcore[im_start:im_end]
        istart_old_proc_id = -1
        iend_old_proc_id = -1
        plist = []

        for iproc_old in range(nproc_previous):
            im_start_old, im_end_old = sp_applications.MPI_start_end(
                psize, nproc_previous, iproc_old
            )
            if (im_start >= im_start_old) and im_start <= im_end_old:
                istart_old_proc_id = iproc_old
            if (im_end >= im_start_old) and im_end <= im_end_old:
                iend_old_proc_id = iproc_old
            plist.append([im_start_old, im_end_old])
        ptl_on_this_cpu = im_start
        nptl_total = 0

        for iproc_index_old in range(istart_old_proc_id, iend_old_proc_id + 1):
            fout = open(
                os.path.join(
                    Tracker["constants"]["refinement_dir"],
                    "main%03d" % selected_iteration,
                    "oldparamstructure",
                    "oldparamstructure_%01d_%03d_%03d.json"
                    % (procid, iproc_index_old, selected_iteration),
                ),
                "r",
            )
            oldparamstructure_on_old_cpu = sp_utilities.convert_json_fromunicode(
                json.load(fout)
            )
            fout.close()
            mlocal_id_on_old = ptl_on_this_cpu - plist[iproc_index_old][0]
            while (mlocal_id_on_old < len(oldparamstructure_on_old_cpu)) and (
                ptl_on_this_cpu < im_end
            ):
                oldparamstructure[procid].append(
                    oldparamstructure_on_old_cpu[mlocal_id_on_old]
                )
                local_dict[local_lcore[nptl_total]] = [
                    Blockdata["myid"],
                    procid,
                    selected_iteration,
                    nptl_total,
                    ptl_on_this_cpu,
                ]
                ptl_on_this_cpu += 1
                mlocal_id_on_old += 1
                nptl_total += 1

        del oldparamstructure_on_old_cpu
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

        for icpu in range(Blockdata["nproc"]):  # dump to disk one by one
            if Blockdata["myid"] == icpu:
                fout = open(
                    os.path.join(
                        Tracker["constants"]["masterdir"],
                        "main%03d" % selected_iteration,
                        "oldparamstructure",
                        "oldparamstructure_%01d_%03d_%03d.json"
                        % (procid, Blockdata["myid"], selected_iteration),
                    ),
                    "w",
                )
                json.dump(oldparamstructure[procid], fout)
                fout.close()
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    # output number of smearing
    smearing_dict = {}
    tchunk = []
    for procid in range(2):
        if Blockdata["myid"] == Blockdata["main_node"]:
            chunk = sp_utilities.read_text_file(
                os.path.join(Tracker["constants"]["masterdir"], "chunk_%d.txt" % procid)
            )
            chunk_size = len(chunk)
            smearing_list = [None for i in range(chunk_size)]
        else:
            chunk_size = 0
        chunk_size = sp_utilities.bcast_number_to_all(
            chunk_size, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        local_smearing_list = []
        for im in range(len(oldparamstructure[procid])):
            local_smearing_list.append(len(oldparamstructure[procid][im][2]))

        if Blockdata["myid"] == Blockdata["main_node"]:
            im_start_old, im_end_old = sp_applications.MPI_start_end(
                chunk_size, Blockdata["nproc"], Blockdata["main_node"]
            )
            for im in range(len(local_smearing_list)):
                smearing_list[im_start_old + im] = local_smearing_list[im]
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

        if Blockdata["myid"] != Blockdata["main_node"]:
            sp_utilities.wrap_mpi_send(
                local_smearing_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
        else:
            for iproc in range(Blockdata["nproc"]):
                if iproc != Blockdata["main_node"]:
                    im_start_old, im_end_old = sp_applications.MPI_start_end(
                        chunk_size, Blockdata["nproc"], iproc
                    )
                    dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
                    for idum in range(len(dummy)):
                        smearing_list[idum + im_start_old] = dummy[idum]
                else:
                    pass

            sp_utilities.write_text_file(
                smearing_list,
                os.path.join(
                    Tracker["constants"]["masterdir"], "smearing_%d.txt" % procid
                ),
            )
            for im in range(len(chunk)):
                smearing_dict[chunk[im]] = smearing_list[im]
            tchunk += chunk
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    if Blockdata["myid"] == Blockdata["main_node"]:
        tchunk.sort()
        all_smearing = [None] * len(tchunk)
        for im in range(len(tchunk)):
            all_smearing[im] = smearing_dict[tchunk[im]]
        sp_utilities.write_text_file(
            all_smearing,
            os.path.join(Tracker["constants"]["masterdir"], "all_smearing.txt"),
        )
        full_dict_list = [None for im in range(Tracker["constants"]["total_stack"])]
        for key, value in list(local_dict.items()):
            full_dict_list[key] = value
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    for icpu in range(Blockdata["nproc"]):
        if Blockdata["myid"] == icpu and Blockdata["myid"] != Blockdata["main_node"]:
            sp_utilities.wrap_mpi_send(
                local_dict, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
        elif Blockdata["myid"] != icpu and Blockdata["myid"] == Blockdata["main_node"]:
            local_dict = sp_utilities.wrap_mpi_recv(icpu, mpi.MPI_COMM_WORLD)
            for key, value in list(local_dict.items()):
                full_dict_list[key] = value
        else:
            pass
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    Tracker["paramstructure_dict"] = os.path.join(
        Tracker["constants"]["masterdir"], "paramstructure_dict.txt"
    )
    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_utilities.write_text_row(full_dict_list, Tracker["paramstructure_dict"])
    return


### 8


def precalculate_shifted_data_for_recons3D(
    prjlist,
    paramstructure,
    refang,
    rshifts,
    delta,
    avgnorms,
    nxinit,
    nnxo,
    nosmearing,
    norm_per_particle=None,
    upweighted=False,
    nsmear=-1,
):
    if norm_per_particle == None:
        norm_per_particle = len(prjlist) * [1.0]
    nnx = prjlist[0].get_xsize()
    nny = prjlist[0].get_ysize()
    if not nosmearing:
        recdata_list = [[] for im in range(len(prjlist))]
        rshifts_shrank = copy.deepcopy(rshifts)
        for im in range(len(rshifts_shrank)):
            rshifts_shrank[im][0] *= old_div(float(nxinit), float(nnxo))
            rshifts_shrank[im][1] *= old_div(float(nxinit), float(nnxo))
        nshifts = len(rshifts_shrank)
    for im in range(len(prjlist)):
        bckgn = prjlist[im].get_attr("bckgnoise")
        ct = prjlist[im].get_attr("ctf")
        group_id = prjlist[im].get_attr("group")

        if nosmearing:
            phi, theta, psi, s2x, s2y = sp_utilities.get_params_proj(
                prjlist[im], xform="xform.projection"
            )
            prjlist[im].set_attr("wprob", 1.0)
            prjlist[im].set_attr("group", group_id)
            prjlist[im].set_attr_dict({"bckgnoise": bckgn, "ctf": ct})
            prjlist[im].set_attr_dict(
                {
                    "padffted": 1,
                    "is_fftpad": 1,
                    "is_fftodd": 0,
                    "is_complex_ri": 1,
                    "is_complex": 1,
                }
            )
            if not upweighted:
                prjlist[im] = sp_filter.filt_table(prjlist[im], bckgn)
            sp_utilities.set_params_proj(
                prjlist[im], [phi, theta, psi, 0.0, 0.0], xform="xform.projection"
            )
        else:
            avgnorm = avgnorms[prjlist[im].get_attr("chunk_id")]
            # if nsmear <=0.0: numbor = len(paramstructure[im][2])
            # else:         numbor = 1
            numbor = len(paramstructure[im][2])
            ipsiandiang = [
                old_div(paramstructure[im][2][i][0], 1000) for i in range(numbor)
            ]
            allshifts = [paramstructure[im][2][i][0] % 1000 for i in range(numbor)]
            probs = [paramstructure[im][2][i][1] for i in range(numbor)]
            tdir = list(set(ipsiandiang))
            data = [None] * nshifts
            for ii in range(len(tdir)):
                #  Find the number of times given projection direction appears on the list, it is the number of different shifts associated with it.
                lshifts = sp_utilities.findall(tdir[ii], ipsiandiang)
                toprab = 0.0
                for ki in range(len(lshifts)):
                    toprab += probs[lshifts[ki]]
                recdata = EMAN2_cppwrap.EMData(nny, nny, 1, False)
                recdata.set_attr("is_complex", 0)
                for ki in range(len(lshifts)):
                    lpt = allshifts[lshifts[ki]]
                    if data[lpt] == None:
                        data[lpt] = sp_fundamentals.fshift(
                            prjlist[im], rshifts_shrank[lpt][0], rshifts_shrank[lpt][1]
                        )
                        data[lpt].set_attr("is_complex", 0)
                    EMAN2_cppwrap.Util.add_img(
                        recdata,
                        EMAN2_cppwrap.Util.mult_scalar(
                            data[lpt], old_div(probs[lshifts[ki]], toprab)
                        ),
                    )
                recdata.set_attr_dict(
                    {
                        "padffted": 1,
                        "is_fftpad": 1,
                        "is_fftodd": 0,
                        "is_complex_ri": 1,
                        "is_complex": 1,
                    }
                )  # preset already
                if not upweighted:
                    recdata = sp_filter.filt_table(recdata, bckgn)
                recdata.set_attr_dict({"bckgnoise": bckgn, "ctf": ct})
                ipsi = tdir[ii] % 100000
                iang = old_div(tdir[ii], 100000)
                sp_utilities.set_params_proj(
                    recdata,
                    [
                        refang[iang][0],
                        refang[iang][1],
                        refang[iang][2] + ipsi * delta,
                        0.0,
                        0.0,
                    ],
                    xform="xform.projection",
                )
                recdata.set_attr(
                    "wprob", old_div(toprab * avgnorm, norm_per_particle[im])
                )
                recdata.set_attr("group", group_id)
                recdata_list[im].append(recdata)
    if nosmearing:
        return prjlist
    else:
        del bckgn, recdata, tdir, ipsiandiang, allshifts, probs, data
        return recdata_list


##### read data/paramstructure ends
###=====<----downsize data---->>>>


def downsize_data_for_sorting(
    original_data, return_real=False, preshift=True, npad=1, norms=None
):
    # The function will read from stack a subset of images specified in partids
    #   and assign to them parameters from partstack with optional CTF application and shifting of the data.
    # So, the lengths of partids and partstack are the same.
    global Tracker, Blockdata
    # functions:
    # read in data
    # apply mask, and prepare focus projection if focus3D is specified
    # return  1. cdata: data for image comparison, always in Fourier format
    #         2. rdata: data for reconstruction, 4nn return real image
    mask2D = sp_utilities.model_circle(
        Tracker["constants"]["radius"],
        Tracker["constants"]["nnxo"],
        Tracker["constants"]["nnxo"],
    )
    shrinkage = old_div(Tracker["nxinit"], float(Tracker["constants"]["nnxo"]))
    radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
    if Tracker["applybckgnoise"]:
        oneover = []
        nnx = len(Blockdata["bckgnoise"][0])
        for i in range(len(Blockdata["bckgnoise"])):
            temp = [0.0] * nnx
            for k in range(nnx):
                if Blockdata["bckgnoise"][i][k] > 0.0:
                    temp[k] = old_div(1.0, numpy.sqrt(Blockdata["bckgnoise"][i][k]))
            oneover.append(temp)
        del temp
    if Tracker["constants"]["focus3D"]:  # focus mask is applied
        if Blockdata["myid"] == Blockdata["main_node"]:
            focus3d = sp_utilities.get_im(Tracker["constants"]["focus3D"])
            focus3d_nx = focus3d.get_xsize()
            if (
                focus3d_nx != Tracker["nxinit"]
            ):  # So the decimated focus volume can be directly used
                focus3d = sp_fundamentals.resample(
                    focus3d, old_div(float(Tracker["nxinit"]), float(focus3d_nx))
                )
        else:
            focus3d = sp_utilities.model_blank(
                Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"]
            )
        sp_utilities.bcast_EMData_to_all(
            focus3d, Blockdata["myid"], Blockdata["main_node"]
        )
        focus3d = sp_projection.prep_vol(focus3d, 1, 1)
    #  Preprocess the data
    nima = len(original_data)
    cdata = [None] * nima
    rdata = [None] * nima
    fdata = [None] * nima  # focusmask projections
    for im in range(nima):
        image = original_data[im].copy()
        chunk_id = image.get_attr("chunk_id")
        try:
            group_id = image.set_attr("group", groupids[im])
        except:
            pass
        particle_group_id = image.get_attr("particle_group")
        phi, theta, psi, s2x, s2y = sp_utilities.get_params_proj(
            image, xform="xform.projection"
        )
        [sx, sy] = image.get_attr("previous_shifts")
        [sx1, sy1] = image.get_attr("current_shifts")
        rimage = sp_fundamentals.cyclic_shift(image, int(round(sx)), int(round(sy)))
        cimage = sp_fundamentals.fshift(image, sx1, sy1)
        st = EMAN2_cppwrap.Util.infomask(rimage, mask2D, False)
        rimage -= st[0]
        rimage = old_div(rimage, st[1])
        st = EMAN2_cppwrap.Util.infomask(cimage, mask2D, False)
        cimage -= st[0]
        cimage = old_div(cimage, st[1])

        if not Tracker["nosmearing"] and norms:
            cimage *= old_div(
                Tracker["avgnorm"][chunk_id], norms[im]
            )  # norm correction

        if Tracker["applybckgnoise"]:
            if Tracker["applymask"]:
                if Tracker["constants"]["hardmask"]:
                    cimage = sp_morphology.cosinemask(
                        cimage, radius=Tracker["constants"]["radius"]
                    )
                else:
                    bckg = sp_utilities.model_gauss_noise(
                        1.0,
                        Tracker["constants"]["nnxo"] + 2,
                        Tracker["constants"]["nnxo"],
                    )
                    bckg.set_attr("is_complex", 1)
                    bckg.set_attr("is_fftpad", 1)
                    bckg = sp_fundamentals.fft(
                        sp_filter.filt_table(bckg, oneover[particle_group_id])
                    )
                    #  Normalize bckg noise in real space, only region actually used.
                    st = EMAN2_cppwrap.Util.infomask(bckg, mask2D, False)
                    bckg -= st[0]
                    bckg = old_div(bckg, st[1])
                    cimage = sp_morphology.cosinemask(
                        cimage, radius=Tracker["constants"]["radius"], bckg=bckg
                    )
        else:
            if Tracker["applymask"]:
                cimage = sp_morphology.cosinemask(
                    cimage, radius=Tracker["constants"]["radius"]
                )
            else:
                pass
        # FT
        rimage = sp_fundamentals.fft(rimage)
        cimage = sp_fundamentals.fft(cimage)
        if Tracker["constants"]["CTF"]:
            ctf_params = rimage.get_attr("ctf")
            rimage = sp_fundamentals.fdecimate(
                rimage,
                Tracker["nxinit"] * npad,
                Tracker["nxinit"] * npad,
                1,
                False,
                False,
            )
            cimage = sp_fundamentals.fdecimate(
                cimage,
                Tracker["nxinit"] * npad,
                Tracker["nxinit"] * npad,
                1,
                False,
                False,
            )
            ctf_params.apix = old_div(ctf_params.apix, shrinkage)
            rimage.set_attr("ctf", ctf_params)
            cimage.set_attr("ctf", ctf_params)
            rimage.set_attr("ctf_applied", 0)
            cimage.set_attr("ctf_applied", 0)
            if return_real:
                rimage = sp_fundamentals.fft(rimage)
        else:
            ctf_params = rimage.get_attr_default("ctf", False)
            if ctf_params:
                ctf_params.apix = old_div(ctf_params.apix, shrinkage)
                rimage.set_attr("ctf", ctf_params)
                rimage.set_attr("ctf_applied", 0)
                cimage.set_attr("ctf", ctf_params)
                cimage.set_attr("ctf_applied", 0)

            rimage = sp_fundamentals.fdecimate(
                rimage, nxinit * npad, nxinit * npad, 1, True, False
            )
            cimage = sp_fundamentals.fdecimate(
                cimage, nxinit * npad, nxinit * npad, 1, True, False
            )
            apix = Tracker["constants"]["pixel_size"]
            rimage.set_attr("apix", old_div(apix, shrinkage))
            cimage.set_attr("apix", old_div(apix, shrinkage))

        cimage.set_attr("padffted", 1)
        cimage.set_attr("npad", npad)
        if not return_real:
            rimage.set_attr("padffted", 1)
            rimage.set_attr("npad", npad)

        sp_utilities.set_params_proj(rimage, [phi, theta, psi, 0.0, 0.0])
        rimage.set_attr("chunk_id", chunk_id)
        # image.set_attr("group", groupids[im])
        rimage.set_attr("particle_group", particle_group_id)

        sp_utilities.set_params_proj(cimage, [phi, theta, psi, 0.0, 0.0])
        cimage.set_attr("chunk_id", chunk_id)
        # cimage.set_attr("group", groupids[im])
        cimage.set_attr("particle_group", particle_group_id)
        rdata[im] = rimage
        cdata[im] = cimage
        if Tracker["applybckgnoise"]:
            rdata[im].set_attr("bckgnoise", Blockdata["bckgnoise"][particle_group_id])
        else:
            rdata[im].set_attr("bckgnoise", Blockdata["bckgnoise"])
            cdata[im].set_attr("bckgnoise", Blockdata["bckgnoise"])
        if Tracker["constants"]["focus3D"]:
            focusmask = sp_morphology.binarize(
                sp_projection.prgl(focus3d, [phi, theta, psi, 0.0, 0.0], 1, True), 1
            )
            cdata[im] = sp_fundamentals.fft(focusmask * sp_fundamentals.fft(cdata[im]))
            if Tracker["constants"]["CTF"]:
                cdata[im].set_attr("ctf", rdata[im].get_attr("ctf"))
            fdata[im] = focusmask
        cdata[im].set_attr("is_complex", 0)
    return cdata, rdata, fdata


##=====<----for 3D----->>>>


def downsize_data_for_rec3D(original_data, particle_size, return_real=False, npad=1):
    # The function will read from stack a subset of images specified in partids
    #   and assign to them parameters from partstack with optional CTF application and shifting of the data.
    # So, the lengths of partids and partstack are the same.
    global Tracker, Blockdata
    # functions:
    # read in data
    # apply mask, and prepare focus projection if focus3D is specified
    # return  1. cdata: data for image comparison, always in Fourier format
    #         2. rdata: data for reconstruction, 4nn return real image
    nima = len(original_data)
    rdata = [None] * nima
    mask2D = sp_utilities.model_circle(
        Tracker["constants"]["radius"],
        Tracker["constants"]["nnxo"],
        Tracker["constants"]["nnxo"],
    )
    shrinkage = old_div(particle_size, float(Tracker["constants"]["nnxo"]))
    radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
    for im in range(nima):
        image = original_data[im].copy()
        chunk_id = image.get_attr("chunk_id")
        try:
            particle_group_id = image.get_attr("particle_group")
        except:
            particle_group_id = -1
        phi, theta, psi, s2x, s2y = sp_utilities.get_params_proj(
            image, xform="xform.projection"
        )
        [sx, sy] = image.get_attr("previous_shifts")  # always for rec3D
        if Tracker["nosmearing"]:
            image = sp_fundamentals.fshift(image, s2x, s2y)
        else:
            image = sp_fundamentals.cyclic_shift(image, int(round(sx)), int(round(sy)))
        st = EMAN2_cppwrap.Util.infomask(image, mask2D, False)
        image -= st[0]
        image = old_div(image, st[1])
        image = sp_fundamentals.fft(image)
        if Tracker["constants"]["CTF"]:
            ctf_params = image.get_attr("ctf")
            image = sp_fundamentals.fdecimate(
                image, particle_size * npad, particle_size * npad, 1, False, False
            )
            ctf_params.apix = old_div(ctf_params.apix, shrinkage)
            image.set_attr("ctf", ctf_params)
            image.set_attr("ctf_applied", 0)
        else:
            ctf_params = image.get_attr_default("ctf", False)
            if ctf_params:
                ctf_params.apix = old_div(ctf_params.apix, shrinkage)
                image.set_attr("ctf", ctf_params)
                image.set_attr("ctf_applied", 0)
            image = sp_fundamentals.fdecimate(
                image, particle_size * npad, particle_size * npad, 1, True, False
            )
            apix = Tracker["constants"]["pixel_size"]
            image.set_attr("apix", old_div(apix, shrinkage))
        if not return_real:
            image.set_attr("padffted", 1)
            image.set_attr("npad", npad)
        image.set_attr("chunk_id", chunk_id)
        image.set_attr("particle_group", particle_group_id)
        sp_utilities.set_params_proj(image, [phi, theta, psi, 0.0, 0.0])
        rdata[im] = image
        if Tracker["applybckgnoise"]:
            rdata[im].set_attr(
                "bckgnoise",
                Blockdata["bckgnoise"][rdata[im].get_attr("particle_group")],
            )
        else:
            rdata[im].set_attr("bckgnoise", Blockdata["bckgnoise"])
    return rdata


### end of downsize

###=====<--- comparison


def compare_two_images_eucd(data, ref_vol, fdata):
    global Tracker, Blockdata
    peaks = len(data) * [None]
    ny = data[0].get_ysize()
    ref_vol = sp_projection.prep_vol(ref_vol, npad=2, interpolation_method=1)
    ctfs = [sp_morphology.ctf_img_real(ny, q.get_attr("ctf")) for q in data]
    qt = float(Tracker["constants"]["nnxo"] * Tracker["constants"]["nnxo"])
    for im in range(len(data)):
        phi, theta, psi, s2x, s2y = sp_utilities.get_params_proj(
            data[im], xform="xform.projection"
        )
        if Tracker["constants"]["focus3D"]:
            rtemp = sp_projection.prgl(ref_vol, [phi, theta, psi, 0.0, 0.0], 1, True)
            rtemp = sp_fundamentals.fft(rtemp * fdata[im])
        else:
            rtemp = sp_projection.prgl(ref_vol, [phi, theta, psi, 0.0, 0.0], 1, False)
        rtemp.set_attr("is_complex", 0)
        if data[im].get_attr("is_complex") == 1:
            data[im].set_attr("is_complex", 0)

        if Tracker["applybckgnoise"]:
            peaks[im] = old_div(
                -EMAN2_cppwrap.Util.sqed(
                    data[im],
                    rtemp,
                    ctfs[im],
                    Blockdata["unrolldata"][data[im].get_attr("particle_group")],
                ),
                qt,
            )
        else:
            peaks[im] = old_div(
                -EMAN2_cppwrap.Util.sqed(
                    data[im], rtemp, ctfs[im], Blockdata["unrolldata"]
                ),
                qt,
            )
    return peaks


#


def compare_two_images_cross(data, ref_vol):
    global Tracker, Blockdata
    ny = data[0].get_ysize()
    peaks = len(data) * [None]
    volft = sp_projection.prep_vol(ref_vol, 2, 1)
    ctfs = [None for im in range(len(data))]
    for im in range(len(data)):
        if im == 0:
            current_ctf = data[im].get_attr("ctf")
            ctfimg = sp_morphology.ctf_img_real(ny, current_ctf)
        else:
            if not sp_utilities.same_ctf(current_ctf, data[im].get_attr("ctf")):
                current_ctf = data[im].get_attr("ctf")
                ctfimg = sp_morphology.ctf_img_real(ny, current_ctf)
        ctfs[im] = ctfimg
    #  Ref is in reciprocal space
    for im in range(len(data)):
        phi, theta, psi, s2x, s2y = sp_utilities.get_params_proj(
            data[im], xform="xform.projection"
        )
        ref = sp_projection.prgl(volft, [phi, theta, psi, 0.0, 0.0], 1, False)
        EMAN2_cppwrap.Util.mulclreal(ref, ctfs[im])
        ref.set_attr("is_complex", 0)
        ref.set_value_at(0, 0, 0.0)
        nrmref = numpy.sqrt(EMAN2_cppwrap.Util.innerproduct(ref, ref, None))
        if data[im].get_attr("is_complex") == 1:
            data[im].set_attr("is_complex", 0)
        if Tracker["constants"]["focus3D"]:
            peaks[im] = old_div(
                EMAN2_cppwrap.Util.innerproduct(ref, data[im], None), nrmref
            )
        else:
            if Tracker["applybckgnoise"]:
                peaks[im] = old_div(
                    EMAN2_cppwrap.Util.innerproduct(
                        ref,
                        data[im],
                        Blockdata["unrolldata"][data[im].get_attr("particle_group")],
                    ),
                    nrmref,
                )
            else:
                peaks[im] = old_div(
                    EMAN2_cppwrap.Util.innerproduct(ref, data[im], None), nrmref
                )
    return peaks


#####==========--------------------utilities of creating random assignments


def create_nrandom_lists(partids, number_of_groups, number_of_runs):
    # the second column denotes orignal particle IDs
    # the first column is randomized group ID
    global Tracker, Blockdata
    if Blockdata["myid"] == Blockdata["main_node"]:
        random_assignment = []
        data_list = sp_utilities.read_text_file(partids, -1)
        if len(data_list) == 1:
            sorting_data_list = data_list[0]
        else:
            sorting_data_list = data_list[1]
        random.seed()
        group_size = old_div(len(sorting_data_list), number_of_groups)
        for index_of_random in range(number_of_runs):
            particle_dict = {}
            ll = copy.deepcopy(sorting_data_list)
            random.shuffle(ll)
            group_list = []
            for index_of_groups in range(number_of_groups):
                if index_of_groups != number_of_groups - 1:
                    for iparticle in ll[
                        index_of_groups
                        * group_size : (index_of_groups + 1)
                        * group_size
                    ]:
                        particle_dict[iparticle] = index_of_groups
                        group_list.append(index_of_groups)
                else:
                    for iparticle in ll[index_of_groups * group_size :]:
                        particle_dict[iparticle] = index_of_groups
                        group_list.append(index_of_groups)
            assignment = []
            for im in range(len(sorting_data_list)):
                assignment.append(
                    [particle_dict[sorting_data_list[im]], sorting_data_list[im]]
                )
            random_assignment.append(assignment)
            del assignment
            del ll
    else:
        random_assignment = 0
    random_assignment = sp_utilities.wrap_mpi_bcast(
        random_assignment, Blockdata["main_node"]
    )
    return random_assignment


def create_nrandom_lists_from_given_pids(
    work_dir, partids, number_of_groups, number_of_runs
):
    # the second column denotes orignal particle IDs
    # the first column is randomized group ID
    global Tracker, Blockdata
    if Blockdata["myid"] == Blockdata["main_node"]:
        random_assignment = []
        data_list = sp_utilities.read_text_file(partids, -1)
        if len(data_list) == 1:
            sorting_data_list = data_list[0]
        else:
            sorting_data_list = data_list[1]
        random.seed()
        group_size = old_div(len(sorting_data_list), number_of_groups)
        for index_of_random in range(number_of_runs):
            particle_dict = {}
            ll = copy.deepcopy(sorting_data_list)
            random.shuffle(ll)
            group_list = []
            for index_of_groups in range(number_of_groups):
                if index_of_groups != number_of_groups - 1:
                    for iparticle in ll[
                        index_of_groups
                        * group_size : (index_of_groups + 1)
                        * group_size
                    ]:
                        particle_dict[iparticle] = index_of_groups
                        group_list.append(index_of_groups)
                else:
                    for iparticle in ll[index_of_groups * group_size :]:
                        particle_dict[iparticle] = index_of_groups
                        group_list.append(index_of_groups)
            assignment = []
            for im in range(len(sorting_data_list)):
                assignment.append(
                    [particle_dict[sorting_data_list[im]], sorting_data_list[im]]
                )
            sp_utilities.write_text_row(
                assignment,
                os.path.join(work_dir, "independent_index_%03d.txt" % index_of_random),
            )
            random_assignment.append(assignment)
            del assignment
            del ll
    else:
        random_assignment = 0
    random_assignment = sp_utilities.wrap_mpi_bcast(
        random_assignment, Blockdata["main_node"]
    )
    return random_assignment


def assign_unaccounted_elements_mpi(glist, clusters, img_per_grp):
    # assign unaccounted images by group probabilities
    global Tracker, Blockdata
    icut = old_div(3 * img_per_grp, 2)
    if Blockdata["myid"] == Blockdata["main_node"]:
        for ic in range(len(clusters)):
            if len(clusters) > (old_div(3 * img_per_grp, 2)):
                random.shuffle(clusters[ic])
                glist += clusters[ic][icut:]
    else:
        glist = 0
        clusters = 0
    clusters = sp_utilities.wrap_mpi_bcast(
        clusters, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    glist = sp_utilities.wrap_mpi_bcast(
        glist, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )

    if len(glist) <= Blockdata["nproc"] * 30:
        if Blockdata["myid"] == Blockdata["main_node"]:
            clusters = assign_unaccounted_inverse_proportion_to_size(
                glist, clusters, img_per_grp
            )
        else:
            clusters = 0
        clusters = sp_utilities.wrap_mpi_bcast(
            clusters, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
    else:
        slist = []
        clist = []
        for ic in range(len(clusters)):
            if len(clusters[ic]) <= img_per_grp:
                slist.append(
                    max(
                        1.0 - old_div(float(len(clusters[ic])), float(img_per_grp)),
                        0.05,
                    )
                )
            else:
                slist.append(0.05)
        uplist = copy.deepcopy(glist)
        image_start, image_end = sp_applications.MPI_start_end(
            len(uplist), Blockdata["nproc"], Blockdata["myid"]
        )
        uplist = uplist[image_start:image_end]
        nsize = old_div(len(uplist), 3)
        for ichunk in range(3):
            if ichunk != 2:
                ulist = uplist[ichunk * nsize : (ichunk + 1) * nsize]
            else:
                ulist = uplist[ichunk * nsize :]
            while len(ulist) > 0:
                im = random.randint(0, len(clusters) - 1)
                random.shuffle(ulist)
                if slist[im] > random.random():
                    clusters[im].append(ulist[0])
                    if len(clusters[ic]) <= img_per_grp:
                        slist[im] = max(
                            1.0 - old_div(float(len(clusters[im])), float(img_per_grp)),
                            0.05,
                        )
                    else:
                        slist[im] = 0.05
                    del ulist[0]
                    if len(ulist) == 0:
                        break
                else:
                    continue
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            if Blockdata["myid"] != Blockdata["main_node"]:
                sp_utilities.wrap_mpi_send(
                    clusters, Blockdata["main_node"], mpi.MPI_COMM_WORLD
                )
            else:
                for iproc in range(Blockdata["nproc"]):
                    if iproc != Blockdata["main_node"]:
                        dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
                        for ic in range(len(clusters)):
                            clusters[ic] += dummy[ic]
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            if Blockdata["myid"] == Blockdata["main_node"]:
                for ic in range(len(clusters)):
                    clusters[ic] = list(set(clusters[ic]))
                    if len(clusters[ic]) <= img_per_grp:
                        slist[ic] = max(
                            1.0 - old_div(float(len(clusters[ic])), float(img_per_grp)),
                            0.05,
                        )
                    else:
                        slist[ic] = 0.05
            else:
                slist = 0
                clusters = 0
            slist = sp_utilities.wrap_mpi_bcast(
                slist, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            clusters = sp_utilities.wrap_mpi_bcast(
                clusters, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
    return clusters


def refilling_global_scheme_mpi(
    clusters, unaccounted_list, number_of_clusters, log_file, swap_ratio
):
    global Tracker, Blockdata
    m = 0
    NACC = 0
    NUACC = len(unaccounted_list)
    for cluster in clusters:
        NACC += len(cluster)
    swap_ratio = old_div(swap_ratio, 100.0)
    N = NUACC + NACC
    m = number_of_clusters - len(clusters)
    if swap_ratio > 0.0:
        if Blockdata["myid"] == Blockdata["main_node"]:
            if int(swap_ratio * NACC) > NUACC:
                unaccounted_list, clusters = swap_clusters_small_NUACC(
                    unaccounted_list, clusters, swap_ratio
                )
            else:
                unaccounted_list, clusters = swap_clusters_large_NUACC(
                    unaccounted_list, clusters, swap_ratio
                )
        else:
            unaccounted_list = 0
            clusters = 0
        clusters = sp_utilities.wrap_mpi_bcast(
            clusters, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        unaccounted_list = sp_utilities.wrap_mpi_bcast(
            unaccounted_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )

    avg_size = old_div(N, number_of_clusters)
    m = number_of_clusters - len(clusters)
    if Blockdata["myid"] == Blockdata["main_node"]:
        large_clusters = []
        for ic in range(len(clusters)):
            if len(clusters[ic]) > 2 * avg_size:
                large_clusters.append(clusters[ic])
    else:
        large_clusters = 0
    large_clusters = sp_utilities.wrap_mpi_bcast(
        large_clusters, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    L = len(large_clusters)

    if m == 0 and L == 0:
        out_clusters = assign_unaccounted_elements_mpi(
            unaccounted_list, clusters, avg_size
        )
    else:
        if m != 0:
            empty_clusters = [[] for ie in range(m)]
        else:
            empty_clusters = []
        out_clusters = fill_no_large_groups_and_unaccounted_to_m_and_rcluster_mpi(
            unaccounted_list, empty_clusters, clusters, NUACC, NACC
        )
    for i in range(len(out_clusters)):
        out_clusters[i].sort()
    return out_clusters


def select_fixed_size_cluster_from_alist(ulist, img_per_grp):
    cluster = []
    random.shuffle(ulist)
    cluster += ulist[0:img_per_grp]
    del ulist[0:img_per_grp]
    return cluster, ulist


def swap_clusters_small_NUACC(glist, clusters, swap_ratio):
    slist = [None for i in range(len(clusters))]
    temp_list = []
    for ic in range(len(clusters)):
        slist[ic] = int(swap_ratio * len(clusters[ic]))
        random.shuffle(clusters[ic])
        temp_list += clusters[ic][0 : slist[ic]]
        del clusters[ic][0 : slist[ic]]
    temp_list += glist
    for ic in range(len(clusters)):
        random.shuffle(temp_list)
        clusters[ic] += temp_list[0 : slist[ic]]
        del temp_list[0 : slist[ic]]
    return temp_list, clusters


def swap_clusters_large_NUACC(glist, clusters, swap_ratio):
    slist = [None for i in range(len(clusters))]
    temp_list = []
    ulist = copy.deepcopy(glist)
    for ic in range(len(clusters)):
        slist[ic] = int(swap_ratio * len(clusters[ic]))
        random.shuffle(clusters[ic])
        temp_list += clusters[ic][0 : slist[ic]]
        del clusters[ic][0 : slist[ic]]
    for ic in range(len(clusters)):
        random.shuffle(ulist)
        clusters[ic] += ulist[0 : slist[ic]]
        del ulist[0 : slist[ic]]
    ulist += temp_list
    random.shuffle(ulist)
    return ulist, clusters


def fill_no_large_groups_and_unaccounted_to_m_and_rcluster_mpi(
    unaccounted_list, empty_clusters, clusters, NUACC, NACC
):
    global Tracker, Blockdata

    N = NUACC + NACC
    m = len(empty_clusters)
    number_of_groups = m + len(clusters)
    avg_size = old_div(N, number_of_groups)
    if m * avg_size < NUACC - avg_size:
        if Blockdata["myid"] == Blockdata["main_node"]:
            for ic in range(len(clusters)):
                random.shuffle(clusters[ic])
                if len(clusters[ic]) > avg_size:
                    unaccounted_list += clusters[ic][avg_size:]
                    del clusters[ic][avg_size:]
            random.shuffle(unaccounted_list)
        else:
            unaccounted_list = 0
        unaccounted_list = sp_utilities.wrap_mpi_bcast(
            unaccounted_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        clusters = sp_utilities.wrap_mpi_bcast(
            clusters, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        tmp_clusters = []
        if m > 0:
            if Blockdata["myid"] == Blockdata["main_node"]:
                for im in range(m):
                    cluster, unaccounted_list = select_fixed_size_cluster_from_alist(
                        unaccounted_list, old_div(avg_size, 2)
                    )
                    tmp_clusters.append(cluster)
            else:
                tmp_clusters = 0
            tmp_clusters = sp_utilities.wrap_mpi_bcast(
                tmp_clusters, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            unaccounted_list = sp_utilities.wrap_mpi_bcast(
                unaccounted_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
        if len(tmp_clusters) > 0:
            for cluster in tmp_clusters:
                clusters.append(cluster)
        clusters = assign_unaccounted_elements_mpi(unaccounted_list, clusters, avg_size)
        del unaccounted_list
    else:
        for a in empty_clusters:
            clusters.append(a)
        clusters = assign_unaccounted_elements_mpi(unaccounted_list, clusters, avg_size)
        del unaccounted_list
    return clusters


def assign_unaccounted_inverse_proportion_to_size(glist, clusters, img_per_grp):
    # assign unaccounted images by group probabilities, single processor version
    ulist = copy.deepcopy(glist)
    number_of_groups = len(clusters)
    slist = []
    for ic in range(len(clusters)):
        if len(clusters[ic]) <= img_per_grp:
            slist.append(
                max(1.0 - old_div(float(len(clusters[ic])), float(img_per_grp)), 0.05)
            )
        else:
            slist.append(0.05)
    nc = 0
    while len(ulist) > 0:
        im = random.randint(0, number_of_groups - 1)
        random.shuffle(ulist)
        r = random.uniform(0.0, 1.0)
        if r < slist[im]:
            clusters[im].append(ulist[0])
            if len(clusters[im]) <= img_per_grp:
                slist[im] = max(
                    1.0 - float(old_div(len(clusters[im]), float(img_per_grp))), 0.05
                )
            else:
                slist[im] = 0.05
            del ulist[0]
            if len(ulist) == 0:
                break
        else:
            continue
    del ulist
    return clusters


def swap_accounted_with_unaccounted_elements_mpi(
    accounted_file, unaccounted_file, log_file, number_of_groups, swap_ratio
):
    global Tracker, Blockdata
    checking_flag = 0

    if Blockdata["myid"] == Blockdata["main_node"]:
        p1 = sp_utilities.read_text_row(accounted_file)
        if len(p1) <= 1:
            checking_flag = 1
    checking_flag = sp_utilities.bcast_number_to_all(
        checking_flag, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )

    if checking_flag == 0:
        tot = 0
        if Blockdata["myid"] == Blockdata["main_node"]:
            p1 = sp_utilities.read_text_row(accounted_file)
            clusters, npart = split_partition_into_ordered_clusters(p1)
            unaccounted_list = sp_utilities.read_text_file(unaccounted_file)
            total_stack = len(unaccounted_list)
            for any in clusters:
                total_stack += len(any)
        else:
            clusters = 0
            unaccounted_list = 0
            total_stack = 0
        clusters = sp_utilities.wrap_mpi_bcast(
            clusters, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        total_stack = sp_utilities.bcast_number_to_all(
            total_stack, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        unaccounted_list = sp_utilities.wrap_mpi_bcast(
            unaccounted_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        clusters = refilling_global_scheme_mpi(
            clusters, unaccounted_list, number_of_groups, log_file, swap_ratio
        )

        if Blockdata["myid"] == Blockdata["main_node"]:
            dlist, assignment_list = merge_classes_into_partition_list(clusters)
            converted_assignment_list = [[], []]
            for im in range(len(assignment_list)):
                for jm in range(2):
                    converted_assignment_list[jm].append(assignment_list[im][jm])
        else:
            converted_assignment_list = 0
        converted_assignment_list = sp_utilities.wrap_mpi_bcast(
            converted_assignment_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
    else:
        assignment_list = create_nrandom_lists(unaccounted_file, number_of_groups, 1)
        assignment_list = assignment_list[0]
        converted_assignment_list = [[], []]
        for jm in range(2):
            for im in range(len(assignment_list)):
                converted_assignment_list[jm].append(assignment_list[im][jm])
    return converted_assignment_list


def patch_to_do_k_means_match_clusters_asg_new(ptp1, ptp2):
    # patch ad hoc elements to make equal number of classes for two partitions and thus two_way comparison becomes feasible
    patch_elements = []
    if len(ptp1) != len(ptp2):
        for i in range(len(ptp1)):
            ptp1[i] = numpy.array(ptp1[i], "int32")
        for i in range(len(ptp2)):
            ptp2[i] = numpy.array(ptp2[i], "int32")
        alist = []
        blist = []
        for a in ptp1:
            if len(a) > 0:
                alist.append(max(a))
        for b in ptp2:
            if len(b) > 0:
                blist.append(max(b))
        if len(alist) > 0 and len(blist) > 0:
            max_number = max(max(alist), max(blist))
        else:
            exit()  # This would never happen
        if len(ptp1) > len(ptp2):
            ndiff = len(ptp1) - len(ptp2)
            for indiff in range(ndiff):
                l = []
                l.append(max_number + indiff + 1)
                patch_elements.append(max_number + indiff + 1)
                l = numpy.array(l, "int32")
                ptp2.append(l)
        else:
            ndiff = len(ptp2) - len(ptp1)
            for indiff in range(ndiff):
                l = []
                l.append(max_number + indiff + 1)
                patch_elements.append(max_number + indiff + 1)
                l = numpy.array(l, "int32")
                ptp1.append(l)
    else:
        for i in range(len(ptp1)):
            ptp1[i] = numpy.array(ptp1[i], "int32")
            ptp2[i] = numpy.array(ptp2[i], "int32")
    newindeces, list_stable, nb_tot_objs = sp_statistics.k_means_match_clusters_asg_new(
        ptp1, ptp2
    )
    new_list_stable = []
    for a in list_stable:
        a.tolist()
        if len(a) > 0:
            new_list_stable.append(a)  # remove empty ones
    return newindeces, new_list_stable, nb_tot_objs, patch_elements


def do_boxes_two_way_comparison_mpi(
    nbox, input_box_parti1, input_box_parti2, depth, log_main
):
    global Tracker, Blockdata
    NT = 1000
    if Blockdata["myid"] == Blockdata["main_node"]:
        stop_generation = 0
        log_main.add(
            "================================================================================================================"
        )
        log_main.add(
            " Two-way comparison of generation %d and layer %d computed between two pairs of independent runs: %d and %d."
            % (Tracker["current_generation"], Tracker["depth"], nbox, nbox + 1)
        )
        log_main.add(
            "----------------------------------------------------------------------------------------------------------------"
        )
        bad_clustering = 0
        ipair = 0
        core1 = sp_utilities.read_text_row(input_box_parti1)
        ptp1, tmp1 = split_partition_into_ordered_clusters(core1)
        core2 = sp_utilities.read_text_row(input_box_parti2)
        ptp2, tmp2 = split_partition_into_ordered_clusters(core2)
    else:
        ptp1 = 0
        ptp2 = 0
    ptp1 = sp_utilities.wrap_mpi_bcast(ptp1, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
    ptp2 = sp_utilities.wrap_mpi_bcast(ptp2, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
    #### before comparison we do a simulation
    gave, gvar = do_random_groups_simulation_mpi(ptp1, ptp2)
    if Blockdata["myid"] == Blockdata["main_node"]:
        #####
        msg = "P0      "
        msg1 = "Group ID"
        length = max(len(ptp1), len(ptp2))
        for im in range(length):
            try:
                msg += "{:8d} ".format(len(ptp1[im]))
            except:
                pass
            msg1 += "{:8d} ".format(im)
        log_main.add(msg1)
        log_main.add(msg)
        msg = "P1      "
        for im in range(len(ptp2)):
            msg += "{:8d} ".format(len(ptp2[im]))
        log_main.add(msg)

        if len(core1) != len(core2):
            sp_global_def.ERROR(
                "Two partitions have different lengths",
                "do_boxes_two_way_comparison",
                1,
                0,
            )

        full_list = []
        for a in core1:
            full_list.append(a[1])
        full_list.sort()
        total_data = len(full_list)
        minimum_group_size = total_data
        maximum_group_size = 0

        newindeces, list_stable, nb_tot_objs, patch_elements = patch_to_do_k_means_match_clusters_asg_new(
            ptp1, ptp2
        )
        ratio_unaccounted = 100.0 - old_div(nb_tot_objs, float(total_data)) * 100.0
        ratio_accounted = old_div(nb_tot_objs, float(total_data)) * 100.0
        new_list = []
        print_matching_pairs(newindeces, log_main)

        Tracker["current_iter_ratio"] = ratio_accounted
        score_list = []
        nclass = 0
        stat_list = []
        tmp_list = []
        log_main.add("               Post-matching results.")
        log_main.add(
            "{:>5} {:>8}  {:^8}   {:>15} {:>22}  {:>5}".format(
                "Group",
                "    size",
                " status ",
                "reproducibility",
                "random reproducibility",
                " std ",
            )
        )
        for index_of_any in range(len(list_stable)):
            any = list_stable[index_of_any]
            any.tolist()
            any.sort()
            score1 = old_div(
                float(len(any)) * 100.0, float(len(ptp1[newindeces[index_of_any][0]]))
            )
            score2 = old_div(
                float(len(any)) * 100.0, float(len(ptp2[newindeces[index_of_any][1]]))
            )
            score3 = (
                old_div(
                    float(
                        (
                            numpy.intersect1d(
                                ptp1[newindeces[index_of_any][0]],
                                ptp2[newindeces[index_of_any][1]],
                            )
                        ).size
                    ),
                    float(
                        (
                            numpy.union1d(
                                ptp1[newindeces[index_of_any][0]],
                                ptp2[newindeces[index_of_any][1]],
                            )
                        ).size
                    ),
                )
                * 100.0
            )
            if (
                score3
                > gave[index_of_any]
                + Tracker["constants"]["random_group_elimination_threshold"]
                * gvar[index_of_any]
            ):
                score_list.append([score1, score2])
                minimum_group_size = min(minimum_group_size, len(any))
                maximum_group_size = max(maximum_group_size, len(any))
                new_list.append(any)
                nclass += 1
                log_main.add(
                    "{:>5} {:>8d}   {:^8}      {:>7.1f}      {:>15.1f}           {:>5.1f}".format(
                        index_of_any,
                        len(any),
                        "accepted",
                        score3,
                        gave[index_of_any],
                        gvar[index_of_any],
                    )
                )
                stat_list.append([score3, gave[index_of_any], gvar[index_of_any]])
                tmp_list.append(len(any) * (-1))
            else:
                log_main.add(
                    "{:>5} {:>8d}   {:^8}      {:>7.1f}      {:>15.1f}           {:>5.1f}".format(
                        index_of_any,
                        len(any),
                        "rejected",
                        score3,
                        gave[index_of_any],
                        gvar[index_of_any],
                    )
                )
        ###
        if len(tmp_list) > 1:
            tmp_new_list = []
            tmp_stat_list = []
            tmp_list = numpy.array(tmp_list, "int32")
            tmp_list = numpy.argsort(tmp_list)
            for ik in range(len(tmp_list)):
                tmp_stat_list.append(stat_list[tmp_list[ik]])
                tmp_new_list.append(new_list[tmp_list[ik]])
            new_list[:] = tmp_new_list[:]
            stat_list[:] = tmp_stat_list[:]
    else:
        nclass = 0
    nclass = sp_utilities.bcast_number_to_all(
        nclass, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )

    if nclass == 0:
        ### redo two way comparison
        if depth > 1:
            if Blockdata["myid"] == Blockdata["main_node"]:
                log_main.add(
                    "There are no clusters larger than the user provided minimum group size %d."
                    % Tracker["constants"]["minimum_grp_size"]
                )
                log_main.add("Sorting eliminates the smallest group, and continues")

                ptp1, ucluster1 = split_partition_into_ordered_clusters_split_ucluster(
                    core1
                )
                ptp2, ucluster2 = split_partition_into_ordered_clusters_split_ucluster(
                    core2
                )
                newindeces, list_stable, nb_tot_objs, patch_elements = patch_to_do_k_means_match_clusters_asg_new(
                    ptp1, ptp2
                )
                if len(list_stable) > 3:  # No change for two groups
                    fake_list = sorted(list_stable, key=len)
                    list_stable.remove(fake_list[0])
                list_stable = sorted(list_stable, key=len, reverse=True)
                accounted_list, new_index = merge_classes_into_partition_list(
                    list_stable
                )
                a = set(full_list)
                b = set(accounted_list)
                unaccounted_list = sorted(list(a.difference(b)))
                log_main.add(
                    "================================================================================================================\n"
                )
            else:
                minimum_group_size, maximum_group_size, new_index, unaccounted_list, bad_clustering, stop_generation, stat_list = (
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                )
            new_index = sp_utilities.wrap_mpi_bcast(
                new_index, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            unaccounted_list = sp_utilities.wrap_mpi_bcast(
                unaccounted_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            stat_list = sp_utilities.wrap_mpi_bcast(
                stat_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            bad_clustering = sp_utilities.bcast_number_to_all(
                bad_clustering, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            stop_generation = sp_utilities.bcast_number_to_all(
                stop_generation, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            maximum_group_size = sp_utilities.bcast_number_to_all(
                maximum_group_size, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            minimum_group_size = sp_utilities.bcast_number_to_all(
                minimum_group_size, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            return (
                minimum_group_size,
                maximum_group_size,
                new_index,
                unaccounted_list,
                bad_clustering,
                stop_generation,
                stat_list,
            )
        else:
            if Blockdata["myid"] == Blockdata["main_node"]:
                bad_clustering = 1
                log_main.add(
                    "There are no clusters larger than the user provided minimum group size %d."
                    % Tracker["constants"]["minimum_grp_size"]
                )
                log_main.add(
                    "The reason can be: (1) There are no groups in the data set. 2. Minimum group size set is too large. (3) Desired number of groups K is too large."
                )
                log_main.add(
                    "================================================================================================================\n"
                )
            else:
                minimum_group_size, maximum_group_size, full_list, bad_clustering, stop_generation, stat_list = (
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                )
            full_list = sp_utilities.wrap_mpi_bcast(
                full_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            stat_list = sp_utilities.wrap_mpi_bcast(
                stat_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            bad_clustering = sp_utilities.bcast_number_to_all(
                bad_clustering, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            stop_generation = sp_utilities.bcast_number_to_all(
                stop_generation, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            maximum_group_size = sp_utilities.bcast_number_to_all(
                maximum_group_size, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            minimum_group_size = sp_utilities.bcast_number_to_all(
                minimum_group_size, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            return (
                minimum_group_size,
                maximum_group_size,
                [],
                full_list,
                bad_clustering,
                stop_generation,
                stat_list,
            )

    elif (
        nclass == 1
    ):  # Force to stop this generation, and output the cluster; do not do any other box comparison
        if Blockdata["myid"] == Blockdata["main_node"]:
            stop_generation = 1
            accounted_list, new_index = merge_classes_into_partition_list(new_list)
            a = set(full_list)
            b = set(accounted_list)
            unaccounted_list = sorted(list(a.difference(b)))
            log_main.add(
                "Only one group found. The program will output it and stop executing the current generation."
            )

            box1_dir = os.path.join(
                Tracker["constants"]["masterdir"],
                "generation_%03d" % Tracker["current_generation"],
                "layer%d" % Tracker["depth"],
                "nbox%d" % nbox,
            )
            box2_dir = os.path.join(
                Tracker["constants"]["masterdir"],
                "generation_%03d" % Tracker["current_generation"],
                "layer%d" % Tracker["depth"],
                "nbox%d" % (nbox + 1),
            )
            gendir = os.path.join(
                Tracker["constants"]["masterdir"],
                "generation_%03d" % Tracker["current_generation"],
            )
            fout = open(os.path.join(box1_dir, "freq_cutoff.json"), "r")
            freq_cutoff_dict1 = sp_utilities.convert_json_fromunicode(json.load(fout))
            fout.close()
            fout = open(os.path.join(box2_dir, "freq_cutoff.json"), "r")
            freq_cutoff_dict2 = sp_utilities.convert_json_fromunicode(json.load(fout))
            fout.close()
            try:
                fout = open(os.path.join(gendir, "freq_cutoff.json"), "r")
                freq_cutoff_dict3 = sp_utilities.convert_json_fromunicode(
                    json.load(fout)
                )
                fout.close()
                freq_cutoff_dict3 = {}
            except:
                freq_cutoff_dict3 = {}
            ncluster = 0
            for im in range(len(newindeces)):
                try:
                    f1 = freq_cutoff_dict1["Cluster_%03d.txt" % newindeces[im][0]]
                    f2 = freq_cutoff_dict2["Cluster_%03d.txt" % newindeces[im][1]]
                    freq_cutoff_dict3["Cluster_%03d.txt" % ncluster] = min(f1, f2)
                    ncluster += 1
                except:
                    pass
            fout = open(os.path.join(gendir, "freq_cutoff.json"), "w")
            json.dump(freq_cutoff_dict3, fout)
            fout.close()
            log_main.add(
                "================================================================================================================\n"
            )
        else:
            minimum_group_size, maximum_group_size, new_index, unaccounted_list, bad_clustering, stop_generation, stat_list = (
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            )
        unaccounted_list = sp_utilities.wrap_mpi_bcast(
            unaccounted_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        new_index = sp_utilities.wrap_mpi_bcast(
            new_index, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        stat_list = sp_utilities.wrap_mpi_bcast(
            stat_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        bad_clustering = sp_utilities.bcast_number_to_all(
            bad_clustering, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        stop_generation = sp_utilities.bcast_number_to_all(
            stop_generation, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        maximum_group_size = sp_utilities.bcast_number_to_all(
            maximum_group_size, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        minimum_group_size = sp_utilities.bcast_number_to_all(
            minimum_group_size, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        return (
            minimum_group_size,
            maximum_group_size,
            new_index,
            unaccounted_list,
            bad_clustering,
            stop_generation,
            stat_list,
        )
    else:
        if Blockdata["myid"] == Blockdata["main_node"]:
            accounted_list, new_index = merge_classes_into_partition_list(new_list)
            a = set(full_list)
            b = set(accounted_list)
            unaccounted_list = sorted(list(a.difference(b)))
            log_main.add(
                " {} {} {} {}".format(
                    "  The number of accounted for images:",
                    len(accounted_list),
                    "  The number of unaccounted for images:",
                    len(unaccounted_list),
                )
            )
            log_main.add(
                "  The current minimum group size: %d and the maximum group size: %d"
                % (minimum_group_size, maximum_group_size)
            )

            if depth <= 1:  # the last layer
                box1_dir = os.path.join(
                    Tracker["constants"]["masterdir"],
                    "generation_%03d" % Tracker["current_generation"],
                    "layer%d" % Tracker["depth"],
                    "nbox%d" % nbox,
                )
                box2_dir = os.path.join(
                    Tracker["constants"]["masterdir"],
                    "generation_%03d" % Tracker["current_generation"],
                    "layer%d" % Tracker["depth"],
                    "nbox%d" % (nbox + 1),
                )
                gendir = os.path.join(
                    Tracker["constants"]["masterdir"],
                    "generation_%03d" % Tracker["current_generation"],
                )
                fout = open(os.path.join(box1_dir, "freq_cutoff.json"), "r")
                freq_cutoff_dict1 = sp_utilities.convert_json_fromunicode(
                    json.load(fout)
                )
                fout.close()
                fout = open(os.path.join(box2_dir, "freq_cutoff.json"), "r")
                freq_cutoff_dict2 = sp_utilities.convert_json_fromunicode(
                    json.load(fout)
                )
                fout.close()
                try:
                    fout = open(os.path.join(gendir, "freq_cutoff.json"), "r")
                    freq_cutoff_dict3 = sp_utilities.convert_json_fromunicode(
                        json.load(fout)
                    )
                    fout.close()
                    freq_cutoff_dict3 = {}
                except:
                    freq_cutoff_dict3 = {}
                ncluster = 0
                for im in range(len(newindeces)):
                    try:
                        f1 = freq_cutoff_dict1["Cluster_%03d.txt" % newindeces[im][0]]
                        f2 = freq_cutoff_dict2["Cluster_%03d.txt" % newindeces[im][1]]
                        freq_cutoff_dict3["Cluster_%03d.txt" % ncluster] = min(f1, f2)
                        ncluster += 1
                    except:
                        pass
                fout = open(os.path.join(gendir, "freq_cutoff.json"), "w")
                json.dump(freq_cutoff_dict3, fout)
                fout.close()
                log_main.add(
                    "================================================================================================================\n"
                )
            else:
                log_main.add(
                    "================================================================================================================\n"
                )
        else:
            minimum_group_size, maximum_group_size, new_index, unaccounted_list, bad_clustering, stop_generation, stat_list = (
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            )
        unaccounted_list = sp_utilities.wrap_mpi_bcast(
            unaccounted_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        new_index = sp_utilities.wrap_mpi_bcast(
            new_index, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        stat_list = sp_utilities.wrap_mpi_bcast(
            stat_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        bad_clustering = sp_utilities.bcast_number_to_all(
            bad_clustering, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        stop_generation = sp_utilities.bcast_number_to_all(
            stop_generation, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        maximum_group_size = sp_utilities.bcast_number_to_all(
            maximum_group_size, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        minimum_group_size = sp_utilities.bcast_number_to_all(
            minimum_group_size, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        return (
            minimum_group_size,
            maximum_group_size,
            new_index,
            unaccounted_list,
            bad_clustering,
            stop_generation,
            stat_list,
        )


def split_partition_into_ordered_clusters_split_ucluster(partition):
    # split groupids from indexes of images
    # reindex groups
    ucluster = []
    clusters = []
    cluster_id = []
    for im in range(len(partition)):
        if partition[im][0] not in cluster_id:
            cluster_id.append(partition[im][0])
    ####
    cluster_dict = {}
    group_change_dict = {}
    new_group_id = 0
    if len(cluster_id) > 1:
        cluster_id.sort()
    for icluster in range(len(cluster_id)):
        one_cluster = []
        for a in partition:
            if a[0] == icluster:
                one_cluster.append(a[1])
                cluster_dict[a[1]] = icluster
        one_cluster.sort()
        if icluster < len(cluster_id) - 1:
            clusters.append(one_cluster)
        else:
            ucluster.append(one_cluster)
        group_change_dict[icluster] = new_group_id
        new_group_id += 1
    # create a partition list:
    return clusters, ucluster[0]


def do_withinbox_two_way_comparison(partition_dir, nbox, nrun, niter):
    global Tracker, Blockdata
    log_list = []
    ## for single node only
    log_list.append(" ")
    log_list.append(
        "----------------------------------------------------------------------------------------------------------------"
    )
    ipair = 0
    core1 = sp_utilities.read_text_row(
        os.path.join(partition_dir, "partition_%03d.txt" % (2 * ipair))
    )
    ptp1, tmp1 = split_partition_into_ordered_clusters(core1)
    core2 = sp_utilities.read_text_row(
        os.path.join(partition_dir, "partition_%03d.txt" % (2 * ipair + 1))
    )
    ptp2, tmp2 = split_partition_into_ordered_clusters(core2)

    log_list.append("       Matching of sorting results of two quasi-independent runs")
    # before comparison
    msg = "P0      "
    msg1 = "Group ID"
    for im in range(len(ptp1)):
        msg += "{:8d} ".format(len(ptp1[im]))
        msg1 += "{:8d} ".format(im)
    log_list.append(msg1)
    log_list.append(msg)
    msg = "P1      "
    for im in range(len(ptp2)):
        msg += "{:8d} ".format(len(ptp2[im]))
    log_list.append(msg)
    if len(core1) != len(core2):
        sp_global_def.ERROR(
            "The two partitions have different lengths",
            "do_withinbox_two_way_comparison",
            1,
            0,
        )
    full_list = []
    for a in core1:
        full_list.append(a[1])
    full_list.sort()
    total_data = len(full_list)
    minimum_group_size = total_data
    maximum_group_size = 0

    newindeces, list_stable, nb_tot_objs, patch_elements = patch_to_do_k_means_match_clusters_asg_new(
        ptp1, ptp2
    )
    ratio_unaccounted = 100.0 - old_div(nb_tot_objs, float(total_data)) * 100.0
    ratio_accounted = old_div(nb_tot_objs, float(total_data)) * 100.0

    log_list.append(" ")
    log_list.append("                        Two-way matching of sorting results.")
    log_list.append(
        "M indicates that respective group of P0 sorting (row number) matches respective group of P1 sorting (column number)"
    )
    msg = "   "
    for i in range(len(newindeces)):
        msg += "{:^3d}".format(i)
    log_list.append(msg)
    for im in range(len(newindeces)):
        msg = "{:^3d}".format(im)
        for jm in range(len(newindeces)):
            not_found = True
            for km in range(len(newindeces)):
                if newindeces[km][0] == im and newindeces[km][1] == jm:
                    msg += "{:^3s}".format(" M")
                    not_found = False
            if not_found:
                msg += "{:^3s}".format("   ")
        log_list.append(msg)

    Tracker["current_iter_ratio"] = ratio_accounted
    score_list = []
    nclass = 0
    log_list.append("               Post-matching results.")
    log_list.append(
        "{:>8} {:>10} {:>17} {:>8} {:>15}".format(
            "Group ID", "   size", "min random size", " status ", "reproducibility"
        )
    )
    current_MGR = get_MGR_from_two_way_comparison(newindeces, ptp1, ptp2, total_data)
    stable_clusters = []
    selected_clusters = []

    for index_of_any in range(len(list_stable)):
        any = list_stable[index_of_any]
        any.tolist()
        any.sort()
        score1 = old_div(
            float(len(any)) * 100.0, float(len(ptp1[newindeces[index_of_any][0]]))
        )
        score2 = old_div(
            float(len(any)) * 100.0, float(len(ptp2[newindeces[index_of_any][1]]))
        )
        score3 = (
            old_div(
                float(
                    (
                        numpy.intersect1d(
                            ptp1[newindeces[index_of_any][0]],
                            ptp2[newindeces[index_of_any][1]],
                        )
                    ).size
                ),
                float(
                    (
                        numpy.union1d(
                            ptp1[newindeces[index_of_any][0]],
                            ptp2[newindeces[index_of_any][1]],
                        )
                    ).size
                ),
            )
            * 100.0
        )
        if len(any) > current_MGR[index_of_any]:
            score_list.append([score1, score2])
            minimum_group_size = min(minimum_group_size, len(any))
            maximum_group_size = max(maximum_group_size, len(any))
            nclass += 1
            log_list.append(
                "{:>8} {:>10d} {:>10d}        {:>8} {:>8.1f}       ".format(
                    index_of_any,
                    len(any),
                    current_MGR[index_of_any],
                    "accepted",
                    score3,
                )
            )
            selected_clusters.append(any)
        else:
            log_list.append(
                "{:>8} {:>10d} {:>10d}        {:>8}  {:>8.1f}       ".format(
                    index_of_any,
                    len(any),
                    current_MGR[index_of_any],
                    "rejected",
                    score3,
                )
            )

    accounted_list, new_index = merge_classes_into_partition_list(selected_clusters)
    a = set(full_list)
    b = set(accounted_list)
    unaccounted_list = sorted(list(a.difference(b)))
    sp_utilities.write_text_row(new_index, os.path.join(partition_dir, "Accounted.txt"))
    sp_utilities.write_text_file(
        unaccounted_list, os.path.join(partition_dir, "Unaccounted.txt")
    )
    log_list.append("  The overall reproducibility is %5.1f%%." % ratio_accounted)
    log_list.append(
        "  The number of accounted for images: %d.  The number of unaccounted for images: %d."
        % (len(accounted_list), len(unaccounted_list))
    )
    log_list.append(
        "  The current minimum group size: %d and the maximum group size: %d."
        % (minimum_group_size, maximum_group_size)
    )
    log_list.append(
        "----------------------------------------------------------------------------------------------------------------"
    )
    return (
        minimum_group_size,
        maximum_group_size,
        selected_clusters,
        unaccounted_list,
        ratio_accounted,
        len(list_stable),
        log_list,
    )


def split_partition_into_ordered_clusters(partition):
    # split groupids from indexes of images
    # reindex groups
    clusters = []
    cluster_id = []
    for im in range(len(partition)):
        if partition[im][0] not in cluster_id:
            cluster_id.append(partition[im][0])
    ####
    cluster_dict = {}
    group_change_dict = {}
    new_group_id = 0
    for icluster in range(len(cluster_id)):
        one_cluster = []
        for a in partition:
            if a[0] == icluster:
                one_cluster.append(a[1])
                cluster_dict[a[1]] = icluster
        # if len(one_cluster)>= Tracker["constants"]["minimum_grp_size"]: # clean small ones
        one_cluster.sort()
        clusters.append(one_cluster)
        group_change_dict[icluster] = new_group_id
        new_group_id += 1

    # create a partition list:
    new_partition = []
    for iptl in range(len(partition)):
        gid = group_change_dict[cluster_dict[partition[iptl][1]]]
        if gid > -1:
            new_partition.append(
                [
                    group_change_dict[cluster_dict[partition[iptl][1]]],
                    partition[iptl][1],
                ]
            )
    return clusters, new_partition


def merge_classes_into_partition_list(classes_list):
    # keep the order of classes
    group_dict = {}
    data_list = []
    new_index = []
    if len(classes_list) > 0:
        for index_of_class in range(len(classes_list)):
            for index_of_particle in range(len(classes_list[index_of_class])):
                data_list.append(classes_list[index_of_class][index_of_particle])
                group_dict[
                    classes_list[index_of_class][index_of_particle]
                ] = index_of_class
        data_list = sorted(data_list)
        for index_of_particle in range(len(data_list)):
            new_index.append(
                [group_dict[data_list[index_of_particle]], data_list[index_of_particle]]
            )
        del group_dict
    else:
        data_list = []
        new_index = [[]]
    return data_list, new_index


def get_sorting_all_params(data):
    global Tracker, Blockdata
    if Blockdata["myid"] == Blockdata["main_node"]:
        total_attr_value_list = [[]] * Tracker["total_stack"]
    else:
        total_attr_value_list = 0
    for myproc in range(Blockdata["nproc"]):
        attr_value_list = 0
        if Blockdata["myid"] == myproc:
            attr_value_list = get_sorting_attr_stack(data)
        attr_value_list = sp_utilities.wrap_mpi_bcast(attr_value_list, myproc)
        if Blockdata["myid"] == Blockdata["main_node"]:
            image_start, image_end = sp_applications.MPI_start_end(
                Tracker["total_stack"], Blockdata["nproc"], myproc
            )
            total_attr_value_list = fill_in_mpi_list(
                total_attr_value_list, attr_value_list, image_start, image_end
            )
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    total_attr_value_list = sp_utilities.wrap_mpi_bcast(
        total_attr_value_list, Blockdata["main_node"]
    )
    return total_attr_value_list


def get_sorting_attr_stack(data_in_core):
    # get partitioned group ID and xform.projection parameters
    attr_value_list = []
    for idat in range(len(data_in_core)):
        attr_value_list.append(
            [
                data_in_core[idat].get_attr("group"),
                sp_utilities.get_params_proj(
                    data_in_core[idat], xform="xform.projection"
                ),
            ]
        )
    return attr_value_list


def fill_in_mpi_list(mpi_list, data_list, index_start, index_end):
    for index in range(index_start, index_end):
        mpi_list[index] = data_list[index - index_start]
    return mpi_list


def parsing_sorting_params(partid, sorting_params_list):
    group_list = []
    ali3d_params_list = []
    partid_list = sp_utilities.read_text_file(partid, -1)
    if len(partid_list) == 1:
        for ielement in range(len(sorting_params_list)):
            group_list.append(
                [sorting_params_list[ielement][0], partid_list[0][ielement]]
            )
            ali3d_params_list.append(sorting_params_list[ielement][1:])
    elif len(partid_list) == 2:
        for ielement in range(len(sorting_params_list)):
            group_list.append(
                [sorting_params_list[ielement][0], partid_list[1][ielement]]
            )
            ali3d_params_list.append(sorting_params_list[ielement][1:])
    else:
        sp_global_def.ERROR("Wrong columns", "parsing_sorting_params", 1, 0)
    return group_list, ali3d_params_list


def get_angle_step_from_number_of_orien_groups(orien_groups):
    global Tracker, Blockdata
    sym_class = Blockdata["symclass"]
    N = orien_groups
    angle_step = 180.0
    while len(sym_class.even_angles(angle_step)) < N:
        angle_step = old_div(angle_step, 2.0)
    while len(sym_class.even_angles(angle_step)) > N:
        angle_step += 0.1
    del sym_class
    return angle_step


def parti_oriens(params, angstep, smc):
    ntot = len(params)
    eah = smc.even_angles(angstep, inc_mirror=0)
    leah = len(eah)
    u = []
    for q in eah:
        m = smc.symmetry_related([(180.0 + q[0]) % 360.0, 180.0 - q[1], 0.0])
        itst = len(u)
        for c in m:
            if smc.is_in_subunit(c[0], c[1], 1):
                if not smc.is_in_subunit(c[0], c[1], 0):
                    u.append(c)
                    break
        if len(u) != itst + 1:
            u.append(q)  # This is for exceptions that cannot be easily handled
    seaf = []
    for q in eah + u:
        seaf += smc.symmetry_related(q)
    lseaf = 2 * leah
    seaf = sp_utilities.angles_to_normals(seaf)
    occupancy = [[] for i in range(leah)]
    for i, q in enumerate(params):
        l = sp_utilities.nearest_fang(seaf, q[0], q[1])
        l = l % lseaf
        if l >= leah:
            l = l - leah
        occupancy[l].append(i)
    return occupancy, eah


def get_angle_step_and_orien_groups_mpi(params_in, partids_in, angstep):
    global Tracker, Blockdata
    if Blockdata["main_node"] == Blockdata["myid"]:
        params = sp_utilities.read_text_row(params_in)
        partids = sp_utilities.read_text_file(partids_in, -1)
        if len(partids) == 1:
            partids = partids[0]
        else:
            partids = partids[1]
        subparams = []
        for im in range(len(partids)):
            subparams.append(params[partids[im]])
        occu, eah = parti_oriens(subparams, angstep, Blockdata["symclass"])
        ptls_in_orien_groups = []
        reassign_list = []
        for l, q in enumerate(occu):
            if len(q) < Tracker["min_orien_group_size"]:
                reassign_list += q
            else:
                ptls_in_orien_groups.append(q)
        random.shuffle(reassign_list)
        for a in reassign_list:
            img = random.randint(0, len(ptls_in_orien_groups) - 1)
            ptls_in_orien_groups[img].append(a)
        del reassign_list
        for img in range(len(ptls_in_orien_groups)):
            tmp = sorted(ptls_in_orien_groups[img])
            ptls_in_orien_groups[img][:] = tmp[:]
    else:
        ptls_in_orien_groups = 0
    ptls_in_orien_groups = sp_utilities.wrap_mpi_bcast(
        ptls_in_orien_groups, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    return ptls_in_orien_groups


##### ================= orientation groups


def compare_two_iterations(assignment1, assignment2, number_of_groups):
    # compare two assignments during clustering, either iteratively or independently
    assigned_groups1 = [[] for i in range(number_of_groups)]
    for im in range(len(assignment1)):
        assigned_groups1[assignment1[im]].append(im)
    res1 = []
    for iref in range(number_of_groups):
        a = numpy.array(assigned_groups1[iref], "int32")
        a.sort()
        res1.append(a)
    assigned_groups2 = [[] for i in range(number_of_groups)]
    for im in range(len(assignment2)):
        assigned_groups2[assignment2[im]].append(im)
    res2 = []
    for iref in range(number_of_groups):
        a = numpy.array(assigned_groups2[iref], "int32")
        a.sort()
        res2.append(a)
        del a
    newindeces, list_stable, nb_tot_objs = sp_statistics.k_means_match_clusters_asg_new(
        res1, res2
    )
    del res1
    del res2
    return old_div(float(nb_tot_objs), len(assignment1)), newindeces, list_stable


def update_data_assignment(cdata, rdata, assignment, proc_list, nosmearing, myid):
    nima = len(cdata)
    groupids = assignment[proc_list[myid][0] : proc_list[myid][1]]
    for im in range(nima):
        try:
            previous_group = cdata[im].get_attr("group")
        except:
            previous_group = -1
        cdata[im].set_attr("group", groupids[im])
        if nosmearing:
            rdata[im].set_attr("group", groupids[im])
            rdata[im].set_attr("previous_group", previous_group)
        else:
            for jm in range(len(rdata[im])):
                rdata[im][jm].set_attr("previous_group", previous_group)
                rdata[im][jm].set_attr("group", groupids[im])
    return


def update_rdata_assignment(assignment, proc_list, myid, rdata):
    nima = len(rdata)
    groupids = assignment[proc_list[myid][0] : proc_list[myid][1]]
    for im in range(nima):
        rdata[im].set_attr("group", groupids[im])
    return


def MPI_volume_start_end(number_of_groups, ncolor, mycolor):
    igroup_start = int(round(old_div(float(number_of_groups), ncolor) * mycolor))
    igroup_end = int(round(old_div(float(number_of_groups), ncolor) * (mycolor + 1)))
    return igroup_start, igroup_end


## conversion


def stepone(tvol, tweight):
    global Tracker, Blockdata
    tvol.set_attr("is_complex", 1)
    ovol = EMAN2_cppwrap.Util.shrinkfvol(tvol, 2)
    owol = EMAN2_cppwrap.Util.shrinkfvol(tweight, 2)
    if Tracker["constants"]["symmetry"] != "c1":
        ovol = ovol.symfvol(Tracker["constants"]["symmetry"], -1)
        owol = owol.symfvol(Tracker["constants"]["symmetry"], -1)
    return EMAN2_cppwrap.Util.divn_cbyr(ovol, owol)


def steptwo_mpi(tvol, tweight, treg, cfsc=None, regularized=True, color=0):
    global Tracker, Blockdata
    if Blockdata["color"] != color:
        return sp_utilities.model_blank(
            1
        )  # This should not be executed if called properly
    if Blockdata["myid_on_node"] == 0:
        nz = tweight.get_zsize()
        ny = tweight.get_ysize()
        nx = tweight.get_xsize()
        tvol.set_attr_dict(
            {"is_complex": 1, "is_fftodd": 1, "is_complex_ri": 1, "is_fftpad": 1}
        )
        if regularized:
            nr = len(cfsc)
            limitres = 0
            for i in range(nr):
                cfsc[i] = min(max(cfsc[i], 0.0), 0.999)
                if cfsc[i] == 0.0:
                    limitres = i - 1
                    break
            if limitres == 0:
                limitres = nr - 2
            ovol = sp_utilities.reshape_1d(cfsc, nr, 2 * nr)
            limitres = 2 * min(
                limitres, Tracker["maxfrad"]
            )  # 2 on account of padding, which is always on
            maxr2 = limitres ** 2
            for i in range(limitres + 1, len(ovol), 1):
                ovol[i] = 0.0
            ovol[0] = 1.0
            it = sp_utilities.model_blank(2 * nr)
            for i in range(2 * nr):
                it[i] = ovol[i]
            del ovol
            #  Do not regularize first four
            for i in range(5):
                treg[i] = 0.0
            EMAN2_cppwrap.Util.reg_weights(tweight, treg, it)
            del it
        else:
            limitres = 2 * min(
                old_div(Tracker["constants"]["nnxo"], 2), Tracker["maxfrad"]
            )
            maxr2 = limitres ** 2
        #  Iterative weights
        if Tracker["constants"]["symmetry"] != "c1":
            tvol = tvol.symfvol(Tracker["constants"]["symmetry"], limitres)
            tweight = tweight.symfvol(Tracker["constants"]["symmetry"], limitres)

    else:
        tvol = sp_utilities.model_blank(1)
        tweight = sp_utilities.model_blank(1)
        nz = 0
        ny = 0
        nx = 0
        maxr2 = 0
    nx = sp_utilities.bcast_number_to_all(
        nx, source_node=0, mpi_comm=Blockdata["shared_comm"]
    )
    ny = sp_utilities.bcast_number_to_all(
        ny, source_node=0, mpi_comm=Blockdata["shared_comm"]
    )
    nz = sp_utilities.bcast_number_to_all(
        nz, source_node=0, mpi_comm=Blockdata["shared_comm"]
    )
    maxr2 = sp_utilities.bcast_number_to_all(
        maxr2, source_node=0, mpi_comm=Blockdata["shared_comm"]
    )

    vol_data = sp_utilities.get_image_data(tvol)
    we_data = sp_utilities.get_image_data(tweight)
    #  tvol is overwritten, meaning it is also an output
    n_iter = 10
    ifi = mpi.mpi_iterefa(
        vol_data.__array_interface__["data"][0],
        we_data.__array_interface__["data"][0],
        nx,
        ny,
        nz,
        maxr2,
        Tracker["constants"]["nnxo"],
        Blockdata["myid_on_node"],
        color,
        Blockdata["no_of_processes_per_group"],
        Blockdata["shared_comm"],
        n_iter,
    )
    if Blockdata["myid_on_node"] == 0:
        #  Either pad or window in F space to 2*nnxo
        nx = tvol.get_ysize()
        if nx > 2 * Tracker["constants"]["nnxo"]:
            tvol = sp_fundamentals.fdecimate(
                tvol,
                2 * Tracker["constants"]["nnxo"],
                2 * Tracker["constants"]["nnxo"],
                2 * Tracker["constants"]["nnxo"],
                False,
                False,
            )
        elif nx < 2 * Tracker["constants"]["nnxo"]:
            tvol = sp_fundamentals.fpol(
                tvol,
                2 * Tracker["constants"]["nnxo"],
                2 * Tracker["constants"]["nnxo"],
                2 * Tracker["constants"]["nnxo"],
                RetReal=False,
                normalize=False,
            )
        tvol = sp_fundamentals.fft(tvol)
        tvol = sp_fundamentals.cyclic_shift(
            tvol,
            Tracker["constants"]["nnxo"],
            Tracker["constants"]["nnxo"],
            Tracker["constants"]["nnxo"],
        )
        tvol = EMAN2_cppwrap.Util.window(
            tvol,
            Tracker["constants"]["nnxo"],
            Tracker["constants"]["nnxo"],
            Tracker["constants"]["nnxo"],
        )
        tvol.div_sinc(1)
        tvol = sp_morphology.cosinemask(
            tvol, old_div(Tracker["constants"]["nnxo"], 2) - 1, 5, None
        )
        return tvol
    else:
        return None


def steptwo_mpi_filter(
    tvol, tweight, treg, cfsc=None, cutoff_freq=0.45, aa=0.01, regularized=True, color=0
):
    global Tracker, Blockdata
    cutoff_freq2, aa = estimate_tanhl_params(
        cutoff_freq, aa, 2 * Tracker["constants"]["nnxo"]
    )
    if Blockdata["color"] != color:
        return sp_utilities.model_blank(
            1
        )  # This should not be executed if called properly
    if Blockdata["myid_on_node"] == 0:
        nz = tweight.get_zsize()
        ny = tweight.get_ysize()
        nx = tweight.get_xsize()
        tvol.set_attr_dict(
            {"is_complex": 1, "is_fftodd": 1, "is_complex_ri": 1, "is_fftpad": 1}
        )
        if regularized:
            nr = len(cfsc)
            limitres = 0
            for i in range(nr):
                cfsc[i] = min(max(cfsc[i], 0.0), 0.999)
                if cfsc[i] == 0.0:
                    limitres = i - 1
                    break
            if limitres == 0:
                limitres = nr - 2
            ovol = sp_utilities.reshape_1d(cfsc, nr, 2 * nr)
            limitres = 2 * min(
                limitres, Tracker["maxfrad"]
            )  # 2 on account of padding, which is always on
            maxr2 = limitres ** 2
            for i in range(limitres + 1, len(ovol), 1):
                ovol[i] = 0.0
            ovol[0] = 1.0
            it = sp_utilities.model_blank(2 * nr)
            for i in range(2 * nr):
                it[i] = ovol[i]
            del ovol
            #  Do not regularize first four
            for i in range(5):
                treg[i] = 0.0
            EMAN2_cppwrap.Util.reg_weights(tweight, treg, it)
            del it
        else:
            limitres = 2 * min(
                old_div(Tracker["constants"]["nnxo"], 2), Tracker["maxfrad"]
            )
            maxr2 = limitres ** 2
        #  Iterative weights
        if Tracker["constants"]["symmetry"] != "c1":
            tvol = tvol.symfvol(Tracker["constants"]["symmetry"], limitres)
            tweight = tweight.symfvol(Tracker["constants"]["symmetry"], limitres)

    else:
        tvol = sp_utilities.model_blank(1)
        tweight = sp_utilities.model_blank(1)
        nz = 0
        ny = 0
        nx = 0
        maxr2 = 0
    nx = sp_utilities.bcast_number_to_all(
        nx, source_node=0, mpi_comm=Blockdata["shared_comm"]
    )
    ny = sp_utilities.bcast_number_to_all(
        ny, source_node=0, mpi_comm=Blockdata["shared_comm"]
    )
    nz = sp_utilities.bcast_number_to_all(
        nz, source_node=0, mpi_comm=Blockdata["shared_comm"]
    )
    maxr2 = sp_utilities.bcast_number_to_all(
        maxr2, source_node=0, mpi_comm=Blockdata["shared_comm"]
    )

    vol_data = sp_utilities.get_image_data(tvol)
    we_data = sp_utilities.get_image_data(tweight)
    #  tvol is overwritten, meaning it is also an output
    n_iter = 10
    ifi = mpi.mpi_iterefa(
        vol_data.__array_interface__["data"][0],
        we_data.__array_interface__["data"][0],
        nx,
        ny,
        nz,
        maxr2,
        Tracker["constants"]["nnxo"],
        Blockdata["myid_on_node"],
        color,
        Blockdata["no_of_processes_per_group"],
        Blockdata["shared_comm"],
        n_iter,
    )
    if Blockdata["myid_on_node"] == 0:
        #  Either pad or window in F space to 2*nnxo
        tvol = sp_filter.filt_tanl(tvol, cutoff_freq2, aa)
        nx = tvol.get_ysize()
        if nx > 2 * Tracker["constants"]["nnxo"]:
            tvol = sp_fundamentals.fdecimate(
                tvol,
                2 * Tracker["constants"]["nnxo"],
                2 * Tracker["constants"]["nnxo"],
                2 * Tracker["constants"]["nnxo"],
                False,
                False,
            )
        elif nx < 2 * Tracker["constants"]["nnxo"]:
            tvol = sp_fundamentals.fpol(
                tvol,
                2 * Tracker["constants"]["nnxo"],
                2 * Tracker["constants"]["nnxo"],
                2 * Tracker["constants"]["nnxo"],
                RetReal=False,
                normalize=False,
            )
        tvol = sp_fundamentals.fft(tvol)
        tvol = sp_fundamentals.cyclic_shift(
            tvol,
            Tracker["constants"]["nnxo"],
            Tracker["constants"]["nnxo"],
            Tracker["constants"]["nnxo"],
        )
        tvol = EMAN2_cppwrap.Util.window(
            tvol,
            Tracker["constants"]["nnxo"],
            Tracker["constants"]["nnxo"],
            Tracker["constants"]["nnxo"],
        )
        tvol.div_sinc(1)
        tvol = sp_morphology.cosinemask(
            tvol, old_div(Tracker["constants"]["nnxo"], 2) - 1, 5, None
        )
        return tvol
    else:
        return None


####=====-----------


def recons3d_4nnsorting_MPI(
    myid,
    main_node,
    prjlist,
    random_subset,
    CTF=True,
    upweighted=True,
    mpi_comm=None,
    target_size=-1,
):
    """
        recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
        Input
            list_of_prjlist: list of lists of projections to be included in the reconstruction
    """
    if mpi_comm == None:
        mpi_comm = mpi.MPI_COMM_WORLD
    imgsize = prjlist[0].get_ysize()  # It can be Fourier, so take y-size
    refvol = sp_utilities.model_blank(target_size)
    refvol.set_attr("fudge", 1.0)
    if CTF:
        do_ctf = 1
    else:
        do_ctf = 0
    fftvol = EMAN2_cppwrap.EMData()
    weight = EMAN2_cppwrap.EMData()
    params = {
        "size": target_size,
        "npad": 2,
        "snr": 1.0,
        "sign": 1,
        "symmetry": "c1",
        "refvol": refvol,
        "fftvol": fftvol,
        "weight": weight,
        "do_ctf": do_ctf,
    }
    r = EMAN2_cppwrap.Reconstructors.get("nn4_ctfw", params)
    r.setup()
    # if norm_per_particle == None: norm_per_particle = len(prjlist)*[1.0]
    for im in range(len(prjlist)):
        phi, theta, psi, s2x, s2y = sp_utilities.get_params_proj(
            prjlist[im], xform="xform.projection"
        )  # shifts are already applied
        if random_subset == 2:
            bckgn = target_size * [1.0]
            if prjlist[im].get_attr("is_complex") == 0:
                prjlist[im] = sp_fundamentals.fft(prjlist[im])
            prjlist[im].set_attr_dict({"padffted": 1, "is_complex": 1})
            if not upweighted:
                prjlist[im] = sp_filter.filt_table(prjlist[im], bckgn)
            prjlist[im].set_attr("bckgnoise", bckgn)
            r.insert_slice(
                prjlist[im],
                EMAN2_cppwrap.Transform(
                    {"type": "spider", "phi": phi, "theta": theta, "psi": psi}
                ),
                1.0,
            )
        else:
            if prjlist[im].get_attr("chunk_id") == random_subset:
                # try:    bckgn = prjlist[im].get_attr("bckgnoise")
                bckgn = target_size * [1.0]
                if prjlist[im].get_attr("is_complex") == 0:
                    prjlist[im] = sp_fundamentals.fft(prjlist[im])
                prjlist[im].set_attr_dict({"padffted": 1, "is_complex": 1})
                if not upweighted:
                    prjlist[im] = sp_filter.filt_table(prjlist[im], bckgn)
                prjlist[im].set_attr("bckgnoise", bckgn)
                r.insert_slice(
                    prjlist[im],
                    EMAN2_cppwrap.Transform(
                        {"type": "spider", "phi": phi, "theta": theta, "psi": psi}
                    ),
                    1.0,
                )
    #  clean stuff
    sp_utilities.reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
    sp_utilities.reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
    if myid == main_node:
        dummy = r.finish(True)
    mpi.mpi_barrier(mpi_comm)
    if myid == main_node:
        return fftvol, weight, refvol
    else:
        return None, None, None


def recons3d_4nnsorting_group_MPI(
    myid,
    main_node,
    prjlist,
    random_subset,
    group_ID,
    CTF=True,
    upweighted=True,
    mpi_comm=None,
    target_size=-1,
):
    """
        recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
        Input
            list_of_prjlist: list of lists of projections to be included in the reconstruction
    """
    if mpi_comm == None:
        mpi_comm = mpi.MPI_COMM_WORLD
    imgsize = prjlist[0].get_ysize()  # It can be Fourier, so take y-size
    refvol = sp_utilities.model_blank(target_size)
    refvol.set_attr("fudge", 1.0)
    if CTF:
        do_ctf = 1
    else:
        do_ctf = 0
    fftvol = EMAN2_cppwrap.EMData()
    weight = EMAN2_cppwrap.EMData()
    try:
        qt = prjlist[0].get_attr("qt")
    except:
        qt = 1.0
    params = {
        "size": target_size,
        "npad": 2,
        "snr": 1.0,
        "sign": 1,
        "symmetry": "c1",
        "refvol": refvol,
        "fftvol": fftvol,
        "weight": weight,
        "do_ctf": do_ctf,
    }
    r = EMAN2_cppwrap.Reconstructors.get("nn4_ctfw", params)
    r.setup()
    for im in range(len(prjlist)):
        phi, theta, psi, s2x, s2y = sp_utilities.get_params_proj(
            prjlist[im], xform="xform.projection"
        )  # shifts are already applied
        if prjlist[im].get_attr("group") == group_ID:
            if random_subset == 2:
                try:
                    bckgn = prjlist[im].get_attr("bckgnoise")
                except:
                    bckgn = target_size * [1.0]
                if prjlist[im].get_attr("is_complex") == 0:
                    image = sp_fundamentals.fft(prjlist[im])
                else:
                    image = prjlist[im].copy()
                image.set_attr_dict({"padffted": 1, "is_complex": 1})
                if not upweighted:
                    image = sp_filter.filt_table(image, bckgn)
                image.set_attr("bckgnoise", bckgn)
                r.insert_slice(
                    image,
                    EMAN2_cppwrap.Transform(
                        {"type": "spider", "phi": phi, "theta": theta, "psi": psi}
                    ),
                    1.0,
                )
            else:
                if prjlist[im].get_attr("chunk_id") == random_subset:
                    try:
                        bckgn = prjlist[im].get_attr("bckgnoise")
                    except:
                        bckgn = target_size * [1.0]
                    if prjlist[im].get_attr("is_complex") == 0:
                        image = sp_fundamentals.fft(prjlist[im])
                    else:
                        image = prjlist[im].copy()
                    image.set_attr_dict({"padffted": 1, "is_complex": 1})
                    if not upweighted:
                        image = sp_filter.filt_table(image, bckgn)
                    image.set_attr("bckgnoise", bckgn)
                    r.insert_slice(
                        image,
                        EMAN2_cppwrap.Transform(
                            {"type": "spider", "phi": phi, "theta": theta, "psi": psi}
                        ),
                        1.0,
                    )
    sp_utilities.reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
    sp_utilities.reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
    if myid == main_node:
        dummy = r.finish(True)
    mpi.mpi_barrier(mpi_comm)
    if myid == main_node:
        return fftvol, weight, refvol
    else:
        return None, None, None


def do3d_sorting(procid, data, myid, mpi_comm=-1):
    global Tracker, Blockdata
    if mpi_comm == -1:
        mpi_comm = mpi.MPI_COMM_WORLD
    if procid == 0:
        if Blockdata["no_of_groups"] > 1:
            if Blockdata["myid"] == Blockdata["nodes"][procid]:
                if os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
                    print("tempdir exists")
                else:
                    try:
                        os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
                    except:
                        print("tempdir exists")
        else:
            if myid == Blockdata["main_node"]:
                if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
                    try:
                        os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
                    except:
                        print("tempdir exists")
                else:
                    print("tempdir exists")
    mpi.mpi_barrier(mpi_comm)

    tvol, tweight, trol = recons3d_4nnsorting_MPI(
        myid=Blockdata["myid"],
        main_node=Blockdata["nodes"][procid],
        prjlist=data,
        random_subset=procid,
        CTF=Tracker["constants"]["CTF"],
        upweighted=False,
        target_size=(2 * Tracker["nxinit"] + 3),
        mpi_comm=mpi_comm,
    )

    if Blockdata["no_of_groups"] > 1:
        if Blockdata["myid"] == Blockdata["nodes"][procid]:
            if procid == 0:
                if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
                    os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
            tvol.set_attr("is_complex", 0)
            tvol.write_image(
                os.path.join(Tracker["directory"], "tempdir", "tvol_%01d.hdf" % procid)
            )
            tweight.write_image(
                os.path.join(
                    Tracker["directory"], "tempdir", "tweight_%01d.hdf" % procid
                )
            )
            trol.write_image(
                os.path.join(Tracker["directory"], "tempdir", "trol_%01d.hdf" % procid)
            )
    else:
        if myid == Blockdata["main_node"]:
            tvol.set_attr("is_complex", 0)
            tvol.write_image(
                os.path.join(Tracker["directory"], "tempdir", "tvol_%01d.hdf" % procid)
            )
            tweight.write_image(
                os.path.join(
                    Tracker["directory"], "tempdir", "tweight_%01d.hdf" % procid
                )
            )
            trol.write_image(
                os.path.join(Tracker["directory"], "tempdir", "trol_%01d.hdf" % procid)
            )
    mpi.mpi_barrier(mpi_comm)
    return


def do3d_sorting_group_insertion(data, randomset=2):
    global Tracker, Blockdata
    if Blockdata["myid"] == Blockdata["last_node"]:
        if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
            os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
    if randomset == 1:
        for index_of_groups in range(Tracker["number_of_groups"]):
            for procid in range(2, 3):
                tvol, tweight, trol = recons3d_4nnsorting_group_MPI(
                    myid=Blockdata["myid"],
                    main_node=Blockdata["nodes"][0],
                    prjlist=data,
                    random_subset=procid,
                    group_ID=index_of_groups,
                    CTF=Tracker["constants"]["CTF"],
                    upweighted=False,
                    target_size=(2 * Tracker["nxinit"] + 3),
                )

                if Blockdata["myid"] == Blockdata["nodes"][procid]:
                    tvol.set_attr("is_complex", 0)
                    tvol.write_image(
                        os.path.join(
                            Tracker["directory"],
                            "tempdir",
                            "tvol_%d_%d.hdf" % (procid, index_of_groups),
                        )
                    )
                    tweight.write_image(
                        os.path.join(
                            Tracker["directory"],
                            "tempdir",
                            "tweight_%d_%d.hdf" % (procid, index_of_groups),
                        )
                    )
                    trol.write_image(
                        os.path.join(
                            Tracker["directory"],
                            "tempdir",
                            "trol_%d_%d.hdf" % (procid, index_of_groups),
                        )
                    )
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    else:
        for index_of_groups in range(Tracker["number_of_groups"]):
            for procid in range(2):
                tvol, tweight, trol = recons3d_4nnsorting_group_MPI(
                    myid=Blockdata["myid"],
                    main_node=Blockdata["nodes"][procid],
                    prjlist=data,
                    random_subset=procid,
                    group_ID=index_of_groups,
                    CTF=Tracker["constants"]["CTF"],
                    upweighted=False,
                    target_size=(2 * Tracker["nxinit"] + 3),
                )

                if Blockdata["myid"] == Blockdata["nodes"][procid]:
                    tvol.set_attr("is_complex", 0)
                    tag = 7007
                    sp_utilities.send_EMData(
                        tvol, Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                    sp_utilities.send_EMData(
                        tweight, Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                    sp_utilities.send_EMData(
                        trol, Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )

                elif Blockdata["myid"] == Blockdata["last_node"]:
                    tag = 7007
                    tvol = sp_utilities.recv_EMData(
                        Blockdata["nodes"][procid], tag, mpi.MPI_COMM_WORLD
                    )
                    tweight = sp_utilities.recv_EMData(
                        Blockdata["nodes"][procid], tag, mpi.MPI_COMM_WORLD
                    )
                    trol = sp_utilities.recv_EMData(
                        Blockdata["nodes"][procid], tag, mpi.MPI_COMM_WORLD
                    )
                    tvol.write_image(
                        os.path.join(
                            Tracker["directory"],
                            "tempdir",
                            "tvol_%d_%d.hdf" % (procid, index_of_groups),
                        )
                    )
                    tweight.write_image(
                        os.path.join(
                            Tracker["directory"],
                            "tempdir",
                            "tweight_%d_%d.hdf" % (procid, index_of_groups),
                        )
                    )
                    trol.write_image(
                        os.path.join(
                            Tracker["directory"],
                            "tempdir",
                            "trol_%d_%d.hdf" % (procid, index_of_groups),
                        )
                    )
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    return


def do3d_sorting_groups_trl_iter(data, iteration):
    global Tracker, Blockdata
    keepgoing = 1
    if Blockdata["myid"] == Blockdata["last_node"]:
        if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
            os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    do3d_sorting_group_insertion(data)
    fsc143 = 0
    fsc05 = 0
    Tracker["fsc143"] = 0
    Tracker["fsc05"] = 0
    res_05 = Tracker["number_of_groups"] * [0]
    res_143 = Tracker["number_of_groups"] * [0]
    #####
    if Blockdata["no_of_groups"] > 1:
        sub_main_node_list = [-1 for i in range(Blockdata["no_of_groups"])]
        for index_of_colors in range(Blockdata["no_of_groups"]):
            for iproc in range(Blockdata["nproc"] - 1):
                if Blockdata["myid"] == iproc:
                    if (
                        Blockdata["color"] == index_of_colors
                        and Blockdata["myid_on_node"] == 0
                    ):
                        sub_main_node_list[index_of_colors] = Blockdata["myid"]
                    sp_utilities.wrap_mpi_send(
                        sub_main_node_list, Blockdata["last_node"], mpi.MPI_COMM_WORLD
                    )
                if Blockdata["myid"] == Blockdata["last_node"]:
                    dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
                    for im in range(len(dummy)):
                        if dummy[im] > -1:
                            sub_main_node_list[im] = dummy[im]
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        sp_utilities.wrap_mpi_bcast(
            sub_main_node_list, Blockdata["last_node"], mpi.MPI_COMM_WORLD
        )
        ####
        if Tracker["number_of_groups"] % Blockdata["no_of_groups"] == 0:
            nbig_loop = old_div(Tracker["number_of_groups"], Blockdata["no_of_groups"])
        else:
            nbig_loop = (
                old_div(Tracker["number_of_groups"], Blockdata["no_of_groups"]) + 1
            )

        big_loop_colors = [[] for i in range(nbig_loop)]
        big_loop_groups = [[] for i in range(nbig_loop)]
        nc = 0
        while nc < Tracker["number_of_groups"]:
            im = old_div(nc, Blockdata["no_of_groups"])
            jm = nc % Blockdata["no_of_groups"]
            big_loop_colors[im].append(jm)
            big_loop_groups[im].append(nc)
            nc += 1
        for iloop in range(nbig_loop):
            for im in range(len(big_loop_colors[iloop])):
                index_of_group = big_loop_groups[iloop][im]
                index_of_colors = big_loop_colors[iloop][im]
                if Blockdata["myid"] == Blockdata["last_node"]:
                    tvol2 = sp_utilities.get_im(
                        os.path.join(Tracker["directory"], "tempdir", "tvol_0_%d.hdf")
                        % index_of_group
                    )
                    tweight2 = sp_utilities.get_im(
                        os.path.join(
                            Tracker["directory"], "tempdir", "tweight_0_%d.hdf"
                        )
                        % index_of_group
                    )
                    treg2 = sp_utilities.get_im(
                        os.path.join(Tracker["directory"], "tempdir", "tvol_0_%d.hdf")
                        % index_of_group
                    )
                    tag = 7007
                    sp_utilities.send_EMData(
                        tvol2,
                        sub_main_node_list[index_of_colors],
                        tag,
                        mpi.MPI_COMM_WORLD,
                    )
                    sp_utilities.send_EMData(
                        tweight2,
                        sub_main_node_list[index_of_colors],
                        tag,
                        mpi.MPI_COMM_WORLD,
                    )
                    sp_utilities.send_EMData(
                        treg2,
                        sub_main_node_list[index_of_colors],
                        tag,
                        mpi.MPI_COMM_WORLD,
                    )
                elif Blockdata["myid"] == sub_main_node_list[index_of_colors]:
                    tag = 7007
                    tvol2 = sp_utilities.recv_EMData(
                        Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                    tweight2 = sp_utilities.recv_EMData(
                        Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                    treg2 = sp_utilities.recv_EMData(
                        Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            for im in range(len(big_loop_colors[iloop])):
                index_of_group = big_loop_groups[iloop][im]
                index_of_colors = big_loop_colors[iloop][im]
                Tracker["maxfrad"] = old_div(Tracker["nxinit"], 2)
                if Blockdata["color"] == index_of_colors:
                    if Blockdata["myid_on_node"] != 0:
                        tvol2 = sp_utilities.model_blank(1)
                        tweight2 = sp_utilities.model_blank(1)
                        treg2 = sp_utilities.model_blank(1)
                    tvol2 = steptwo_mpi(
                        tvol2, tweight2, treg2, None, False, color=index_of_colors
                    )  # has to be False!!!
                    del tweight2, treg2
                mpi.mpi_barrier(Blockdata["shared_comm"])
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            for im in range(len(big_loop_colors[iloop])):
                index_of_group = big_loop_groups[iloop][im]
                index_of_colors = big_loop_colors[iloop][im]
                if (Blockdata["color"] == index_of_colors) and (
                    Blockdata["myid_on_node"] == 0
                ):
                    tag = 7007
                    sp_utilities.send_EMData(
                        tvol2, Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                elif Blockdata["myid"] == Blockdata["last_node"]:
                    tag = 7007
                    tvol2 = sp_utilities.recv_EMData(
                        sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD
                    )
                    tvol2.write_image(
                        os.path.join(
                            Tracker["directory"],
                            "vol_unfiltered_0_grp%03d_iter%03d.hdf"
                            % (index_of_group, iteration),
                        )
                    )
                    del tvol2
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            for im in range(len(big_loop_colors[iloop])):
                index_of_group = big_loop_groups[iloop][im]
                index_of_colors = big_loop_colors[iloop][im]
                if Blockdata["myid"] == Blockdata["last_node"]:
                    tvol2 = sp_utilities.get_im(
                        os.path.join(Tracker["directory"], "tempdir", "tvol_1_%d.hdf")
                        % index_of_group
                    )
                    tweight2 = sp_utilities.get_im(
                        os.path.join(
                            Tracker["directory"], "tempdir", "tweight_1_%d.hdf"
                        )
                        % index_of_group
                    )
                    treg2 = sp_utilities.get_im(
                        os.path.join(
                            Tracker["directory"],
                            "tempdir",
                            "trol_1_%d.hdf" % index_of_group,
                        )
                    )
                    tag = 7007
                    sp_utilities.send_EMData(
                        tvol2,
                        sub_main_node_list[index_of_colors],
                        tag,
                        mpi.MPI_COMM_WORLD,
                    )
                    sp_utilities.send_EMData(
                        tweight2,
                        sub_main_node_list[index_of_colors],
                        tag,
                        mpi.MPI_COMM_WORLD,
                    )
                    sp_utilities.send_EMData(
                        treg2,
                        sub_main_node_list[index_of_colors],
                        tag,
                        mpi.MPI_COMM_WORLD,
                    )
                elif Blockdata["myid"] == sub_main_node_list[index_of_colors]:
                    tag = 7007
                    tvol2 = sp_utilities.recv_EMData(
                        Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                    tweight2 = sp_utilities.recv_EMData(
                        Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                    treg2 = sp_utilities.recv_EMData(
                        Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            for im in range(len(big_loop_colors[iloop])):
                index_of_group = big_loop_groups[iloop][im]
                index_of_colors = big_loop_colors[iloop][im]
                Tracker["maxfrad"] = old_div(Tracker["nxinit"], 2)
                if Blockdata["color"] == index_of_colors:
                    if Blockdata["myid_on_node"] != 0:
                        tvol2 = sp_utilities.model_blank(1)
                        tweight2 = sp_utilities.model_blank(1)
                        treg2 = sp_utilities.model_blank(1)
                    tvol2 = steptwo_mpi(
                        tvol2, tweight2, treg2, None, False, color=index_of_colors
                    )  # has to be False!!!
                    del tweight2, treg2
                # if( Blockdata["myid_on_node"] == 0):
                #    tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
                mpi.mpi_barrier(Blockdata["shared_comm"])
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            for im in range(len(big_loop_colors[iloop])):
                index_of_group = big_loop_groups[iloop][im]
                index_of_colors = big_loop_colors[iloop][im]
                if (Blockdata["color"] == index_of_colors) and (
                    Blockdata["myid_on_node"] == 0
                ):
                    tag = 7007
                    sp_utilities.send_EMData(
                        tvol2, Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                elif Blockdata["myid"] == Blockdata["last_node"]:
                    tag = 7007
                    tvol2 = sp_utilities.recv_EMData(
                        sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD
                    )
                    tvol2.write_image(
                        os.path.join(
                            Tracker["directory"],
                            "vol_unfiltered_1_grp%03d_iter%03d.hdf"
                            % (index_of_group, iteration),
                        )
                    )
                    del tvol2
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    else:
        Tracker["maxfrad"] = old_div(Tracker["nxinit"], 2)
        for index_of_group in range(Tracker["number_of_groups"]):
            for iprocid in range(2):
                if Blockdata["myid"] == Blockdata["last_node"]:
                    tvol2 = sp_utilities.get_im(
                        os.path.join(Tracker["directory"], "tempdir", "tvol_%d_%d.hdf")
                        % (iprocid, index_of_group)
                    )
                    tweight2 = sp_utilities.get_im(
                        os.path.join(
                            Tracker["directory"], "tempdir", "tweight_%d_%d.hdf"
                        )
                        % (iprocid, index_of_group)
                    )
                    treg2 = sp_utilities.get_im(
                        os.path.join(
                            Tracker["directory"],
                            "tempdir",
                            "trol_%d_%d.hdf" % (iprocid, index_of_group),
                        )
                    )
                    tag = 7007
                    sp_utilities.send_EMData(
                        tvol2, Blockdata["main_node"], tag, mpi.MPI_COMM_WORLD
                    )
                    sp_utilities.send_EMData(
                        tweight2, Blockdata["main_node"], tag, mpi.MPI_COMM_WORLD
                    )
                    sp_utilities.send_EMData(
                        treg2, Blockdata["main_node"], tag, mpi.MPI_COMM_WORLD
                    )
                elif Blockdata["myid"] == Blockdata["main_node"]:
                    tag = 7007
                    tvol2 = sp_utilities.recv_EMData(
                        Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                    tweight2 = sp_utilities.recv_EMData(
                        Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                    treg2 = sp_utilities.recv_EMData(
                        Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
                if Blockdata["myid"] != Blockdata["main_node"]:
                    tvol2 = sp_utilities.model_blank(1)
                    tweight2 = sp_utilities.model_blank(1)
                    treg2 = sp_utilities.model_blank(1)
                tvol2 = steptwo_mpi(
                    tvol2, tweight2, treg2, None, False, color=0
                )  # has to be False!!!
                del tweight2, treg2
                if Blockdata["myid"] == Blockdata["main_node"]:
                    # tvol2 = cosinemask(tvol2, radius = Tracker["constants"]["radius"])
                    tvol2.write_image(
                        os.path.join(
                            Tracker["directory"],
                            "vol_unfiltered_%d_grp%03d_iter%03d.hdf"
                            % (iprocid, index_of_group, iteration),
                        )
                    )
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    keepgoing = sp_utilities.bcast_number_to_all(
        keepgoing, source_node=Blockdata["main_node"], mpi_comm=mpi.MPI_COMM_WORLD
    )  # always check
    Tracker = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"])
    if keepgoing == 0:
        sp_global_def.ERROR(
            "do3d_sorting_groups_trl_iter  %s"
            % os.path.join(Tracker["directory"], "tempdir"),
            "do3d_sorting_groups_trl_iter",
            1,
            Blockdata["myid"],
        )
    return


# Three ways of importing refinement results


def get_input_from_sparx_ref3d(log_main):  # case one
    # import SPARX results
    global Tracker, Blockdata
    import_from_sparx_refinement = 1
    selected_iter = 0
    Tracker_refinement = 0
    checking_flag = 0
    if Blockdata["myid"] == Blockdata["main_node"]:
        if not os.path.exists(Tracker["constants"]["refinement_dir"]):
            checking_flag = 1
    checking_flag = sp_utilities.bcast_number_to_all(
        checking_flag, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    if checking_flag == 1:
        sp_global_def.ERROR(
            "SPARX refinement directory does not exist",
            "get_input_from_sparx_ref3d",
            1,
            Blockdata["myid"],
        )
    if Blockdata["myid"] == Blockdata["main_node"]:
        if (
            Tracker["constants"]["niter_for_sorting"] == -1
        ):  # take the best solution to do sorting
            niter_refinement = 0
            while os.path.exists(
                os.path.join(
                    Tracker["constants"]["refinement_dir"],
                    "main%03d" % niter_refinement,
                )
            ) and os.path.exists(
                os.path.join(
                    Tracker["constants"]["refinement_dir"],
                    "main%03d" % niter_refinement,
                    "Tracker_%03d.json" % niter_refinement,
                )
            ):
                niter_refinement += 1
            niter_refinement -= 1
            if niter_refinement != 0:
                fout = open(
                    os.path.join(
                        Tracker["constants"]["refinement_dir"],
                        "main%03d" % niter_refinement,
                        "Tracker_%03d.json" % niter_refinement,
                    ),
                    "r",
                )
                Tracker_refinement = sp_utilities.convert_json_fromunicode(
                    json.load(fout)
                )
                fout.close()
                selected_iter = Tracker_refinement["constants"]["best"]
            else:
                import_from_sparx_refinement = 0
        else:
            try:
                fout = open(
                    os.path.join(
                        Tracker["constants"]["refinement_dir"],
                        "main%03d" % Tracker["constants"]["niter_for_sorting"],
                        "Tracker_%03d.json" % Tracker["constants"]["niter_for_sorting"],
                    ),
                    "r",
                )
                Tracker_refinement = sp_utilities.convert_json_fromunicode(
                    json.load(fout)
                )
                fout.close()
                selected_iter = Tracker["constants"]["niter_for_sorting"]
            except:
                import_from_sparx_refinement = 0
    else:
        selected_iter = -1
    selected_iter = sp_utilities.bcast_number_to_all(
        selected_iter, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    import_from_sparx_refinement = sp_utilities.bcast_number_to_all(
        import_from_sparx_refinement, source_node=Blockdata["main_node"]
    )
    if import_from_sparx_refinement == 0:
        sp_global_def.ERROR(
            "The best solution is not found",
            "get_input_from_sparx_ref3d",
            1,
            Blockdata["myid"],
        )
        mpi.mpi_finalize()
        exit()
    Tracker_refinement = sp_utilities.wrap_mpi_bcast(
        Tracker_refinement, Blockdata["main_node"], communicator=mpi.MPI_COMM_WORLD
    )

    # Check orgstack, set correct path

    if Blockdata["myid"] == Blockdata["main_node"]:
        Tracker["constants"]["orgstack"] = Tracker_refinement["constants"]["stack"]
        try:
            image = sp_utilities.get_im(Tracker["constants"]["orgstack"], 0)
        except:
            refinement_dir_path, refinement_dir_name = os.path.split(
                Tracker["constants"]["refinement_dir"]
            )
            if Tracker_refinement["constants"]["stack"][0:4] == "bdb:":
                refinement_stack = "bdb:" + os.path.join(
                    refinement_dir_path, Tracker_refinement["constants"]["stack"][4:]
                )
            else:
                refinement_stack = os.path.join(
                    refinement_dir_path, Tracker_refinement["constants"]["stack"]
                )  # very rare case

            if not Tracker["constants"]["orgstack"]:
                Tracker["constants"]["orgstack"] = refinement_stack
            try:
                image = sp_utilities.get_im(Tracker["constants"]["orgstack"], 0)
            except:
                import_from_sparx_refinement = 0
        try:
            total_stack = EMAN2_cppwrap.EMUtil.get_image_count(
                Tracker["constants"]["orgstack"]
            )
        except:
            total_stack = 0
    else:
        total_stack = 0

    import_from_sparx_refinement = sp_utilities.bcast_number_to_all(
        import_from_sparx_refinement, source_node=Blockdata["main_node"]
    )
    if import_from_sparx_refinement == 0:
        sp_global_def.ERROR(
            "The data stack is not accessible",
            "get_input_from_sparx_ref3d",
            1,
            Blockdata["myid"],
        )
    total_stack = sp_utilities.bcast_number_to_all(
        total_stack, source_node=Blockdata["main_node"]
    )
    Tracker["constants"]["total_stack"] = total_stack

    # Now copy relevant refinement files to sorting directory:
    if Blockdata["myid"] == Blockdata["main_node"]:
        if os.path.exists(
            os.path.join(
                Tracker["constants"]["refinement_dir"],
                "main%03d" % selected_iter,
                "params_%03d.txt" % selected_iter,
            )
        ):
            shutil.copyfile(
                os.path.join(
                    Tracker["constants"]["refinement_dir"],
                    "main%03d" % selected_iter,
                    "params_%03d.txt" % selected_iter,
                ),
                os.path.join(
                    Tracker["constants"]["masterdir"], "sparx_refinement_params.txt"
                ),
            )
        else:
            import_from_sparx_refinement = 0
        Tracker["constants"]["selected_iter"] = selected_iter
    import_from_sparx_refinement = sp_utilities.bcast_number_to_all(
        import_from_sparx_refinement, source_node=Blockdata["main_node"]
    )
    if import_from_sparx_refinement == 0:
        sp_global_def.ERROR(
            "The parameter file of the best solution is not accessible",
            "get_input_from_sparx_ref3d",
            1,
            Blockdata["myid"],
        )

    if Blockdata["myid"] == Blockdata["main_node"]:
        if os.path.exists(
            os.path.join(
                Tracker["constants"]["refinement_dir"],
                "main%03d" % selected_iter,
                "bckgnoise.hdf",
            )
        ):
            shutil.copyfile(
                os.path.join(
                    Tracker["constants"]["refinement_dir"],
                    "main%03d" % selected_iter,
                    "bckgnoise.hdf",
                ),
                os.path.join(Tracker["constants"]["masterdir"], "bckgnoise.hdf"),
            )
        else:
            import_from_sparx_refinement == 0
            for search_iter in range(selected_iter - 1, 0, -1):
                if os.path.exists(
                    os.path.join(
                        Tracker["constants"]["refinement_dir"],
                        "main%03d" % search_iter,
                        "bckgnoise.hdf",
                    )
                ):
                    shutil.copyfile(
                        os.path.join(
                            Tracker["constants"]["refinement_dir"],
                            "main%03d" % search_iter,
                            "bckgnoise.hdf",
                        ),
                        os.path.join(
                            Tracker["constants"]["masterdir"], "bckgnoise.hdf"
                        ),
                    )
                    import_from_sparx_refinement = 1
                    break
    import_from_sparx_refinement = sp_utilities.bcast_number_to_all(
        import_from_sparx_refinement, source_node=Blockdata["main_node"]
    )

    if import_from_sparx_refinement == 0:
        Tracker["bckgnoise"] = None
        if Blockdata["myid"] == Blockdata["main_node"]:
            print("Noise file is not found. However we continue")
    else:
        Tracker["bckgnoise"] = os.path.join(
            Tracker["constants"]["masterdir"], "bckgnoise.hdf"
        )

    import_from_sparx_refinement = 1
    if Blockdata["myid"] == Blockdata["main_node"]:
        if os.path.exists(
            os.path.join(
                Tracker["constants"]["refinement_dir"],
                "main%03d" % selected_iter,
                "driver_%03d.txt" % selected_iter,
            )
        ):
            shutil.copyfile(
                os.path.join(
                    Tracker["constants"]["refinement_dir"],
                    "main%03d" % selected_iter,
                    "driver_%03d.txt" % selected_iter,
                ),
                os.path.join(Tracker["constants"]["masterdir"], "fsc_global.txt"),
            )
        else:
            import_from_sparx_refinement = 0

        if import_from_sparx_refinement:
            fsc_curve = sp_utilities.read_text_row(
                os.path.join(Tracker["constants"]["masterdir"], "fsc_global.txt")
            )
        fsc143 = 0
        fsc05 = 0
        for ifreq in range(len(fsc_curve)):  # drive has only one column
            if fsc_curve[ifreq][0] < 0.5:
                break
        fsc05 = ifreq - 1
        for ifreq in range(len(fsc_curve)):
            if fsc_curve[ifreq][0] < 0.143:
                break
        fsc143 = ifreq - 1
        Tracker["constants"]["fsc143"] = fsc143
        Tracker["constants"]["fsc05"] = fsc05

    import_from_sparx_refinement = sp_utilities.bcast_number_to_all(
        import_from_sparx_refinement, source_node=Blockdata["main_node"]
    )
    if import_from_sparx_refinement == 0:
        sp_global_def.ERROR(
            "The driver of the best solution is not accessible",
            "get_input_from_sparx_ref3d",
            1,
            Blockdata["myid"],
        )
    if Blockdata["myid"] == Blockdata["main_node"]:
        if os.path.exists(
            os.path.join(
                Tracker["constants"]["refinement_dir"], "main000/indexes_000.txt"
            )
        ):
            shutil.copyfile(
                os.path.join(
                    Tracker["constants"]["refinement_dir"], "main000/indexes_000.txt"
                ),
                os.path.join(Tracker["constants"]["masterdir"], "indexes.txt"),
            )
        else:
            import_from_sparx_refinement = 0
    import_from_sparx_refinement = sp_utilities.bcast_number_to_all(
        import_from_sparx_refinement, source_node=Blockdata["main_node"]
    )

    if import_from_sparx_refinement == 0:
        sp_global_def.ERROR(
            "The index file of the best solution are not accessible",
            "get_input_from_sparx_ref3d",
            1,
            Blockdata["myid"],
        )
    if Blockdata["myid"] == Blockdata["main_node"]:
        if os.path.exists(
            os.path.join(
                Tracker["constants"]["refinement_dir"], "main000/chunk_0_000.txt"
            )
        ):
            shutil.copyfile(
                os.path.join(
                    Tracker["constants"]["refinement_dir"], "main000/chunk_0_000.txt"
                ),
                os.path.join(Tracker["constants"]["masterdir"], "chunk_0.txt"),
            )
        else:
            import_from_sparx_refinement == 0
        if os.path.exists(
            os.path.join(
                Tracker["constants"]["refinement_dir"], "main000/chunk_1_000.txt"
            )
        ):
            shutil.copyfile(
                os.path.join(
                    Tracker["constants"]["refinement_dir"], "main000/chunk_1_000.txt"
                ),
                os.path.join(Tracker["constants"]["masterdir"], "chunk_1.txt"),
            )
        else:
            import_from_sparx_refinement == 0

        if os.path.exists(
            os.path.join(
                Tracker["constants"]["refinement_dir"], "main000/particle_groups_0.txt"
            )
        ):
            shutil.copyfile(
                os.path.join(
                    Tracker["constants"]["refinement_dir"],
                    "main000/particle_groups_0.txt",
                ),
                os.path.join(
                    Tracker["constants"]["masterdir"], "particle_groups_0.txt"
                ),
            )
        else:
            import_from_sparx_refinement == 0

        if os.path.exists(
            os.path.join(
                Tracker["constants"]["refinement_dir"], "main000/particle_groups_1.txt"
            )
        ):
            shutil.copyfile(
                os.path.join(
                    Tracker["constants"]["refinement_dir"],
                    "main000/particle_groups_1.txt",
                ),
                os.path.join(
                    Tracker["constants"]["masterdir"], "particle_groups_1.txt"
                ),
            )
        else:
            import_from_sparx_refinement == 0
    import_from_sparx_refinement = sp_utilities.bcast_number_to_all(
        import_from_sparx_refinement, source_node=Blockdata["main_node"]
    )
    if import_from_sparx_refinement == 0:
        sp_global_def.ERROR(
            "The chunk files and partice group files are not accessible",
            "get_input_from_sparx_ref3d",
            1,
            Blockdata["myid"],
        )

    # copy all relavant parameters into sorting tracker
    if Blockdata["myid"] == Blockdata["main_node"]:
        if Tracker["constants"]["radius"] == -1:
            Tracker["constants"]["radius"] = Tracker_refinement["constants"]["radius"]
        Tracker["constants"]["nnxo"] = Tracker_refinement["constants"]["nnxo"]
        Tracker["constants"]["orgres"] = Tracker_refinement["bestres"]
        Tracker["delta"] = Tracker_refinement["delta"]
        Tracker["ts"] = Tracker_refinement["ts"]
        Tracker["xr"] = Tracker_refinement["xr"]
        Tracker["constants"]["pixel_size"] = Tracker_refinement["constants"][
            "pixel_size"
        ]
        Tracker["avgnorm"] = Tracker_refinement["avgvaradj"]
        if Tracker["constants"]["nxinit"] < 0:
            Tracker["nxinit_refinement"] = Tracker_refinement[
                "nxinit"
            ]  # Sphire window size
        else:
            Tracker["nxinit_refinement"] = Tracker["constants"][
                "nxinit"
            ]  # User defined window size

        try:
            sym = Tracker_refinement["constants"]["sym"]
        except:
            sym = Tracker_refinement["constants"]["symmetry"]
        if sym != "c1" and Tracker["constants"]["symmetry"] == "c1":
            Tracker["constants"]["symmetry"] = sym
            update_sym = 1
        else:
            update_sym = 0

        if not Tracker["constants"]["mask3D"]:
            if Tracker_refinement["constants"]["mask3D"] and (
                not Tracker["constants"]["do_not_use_3dmask"]
            ):
                refinement_mask3D_path, refinement_mask3D_file = os.path.split(
                    Tracker_refinement["constants"]["mask3D"]
                )  # MRK_DEBUG
                shutil.copyfile(
                    os.path.join(
                        refinement_dir_path, Tracker_refinement["constants"]["mask3D"]
                    ),
                    os.path.join(
                        Tracker["constants"]["masterdir"], refinement_mask3D_file
                    ),
                )
                Tracker["constants"]["mask3D"] = os.path.join(
                    Tracker["constants"]["masterdir"], refinement_mask3D_file
                )
    else:
        update_sym = 0
        Tracker = 0
    Tracker = sp_utilities.wrap_mpi_bcast(
        Tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    update_sym = sp_utilities.bcast_number_to_all(
        update_sym, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )

    if update_sym == 1:
        Blockdata["symclass"] = sp_fundamentals.symclass(
            Tracker["constants"]["symmetry"]
        )
        Tracker["constants"]["orientation_groups"] = max(
            4,
            old_div(
                Tracker["constants"]["orientation_groups"], Blockdata["symclass"].nsym
            ),
        )

    # Setting for margin error
    chunk_dict = {}
    group_dict = {}
    if Blockdata["myid"] == Blockdata["main_node"]:
        chunk_one = sp_utilities.read_text_file(
            os.path.join(Tracker["constants"]["masterdir"], "chunk_0.txt")
        )
        chunk_two = sp_utilities.read_text_file(
            os.path.join(Tracker["constants"]["masterdir"], "chunk_1.txt")
        )
    else:
        chunk_one = 0
        chunk_two = 0
    chunk_one = sp_utilities.wrap_mpi_bcast(chunk_one, Blockdata["main_node"])
    chunk_two = sp_utilities.wrap_mpi_bcast(chunk_two, Blockdata["main_node"])
    #
    if Blockdata["myid"] == Blockdata["main_node"]:
        chunk_one_group = sp_utilities.read_text_file(
            os.path.join(Tracker["constants"]["masterdir"], "particle_groups_0.txt")
        )
        chunk_two_group = sp_utilities.read_text_file(
            os.path.join(Tracker["constants"]["masterdir"], "particle_groups_1.txt")
        )
    else:
        chunk_one_group = 0
        chunk_two_group = 0
    chunk_one_group = sp_utilities.wrap_mpi_bcast(
        chunk_one_group, Blockdata["main_node"]
    )
    chunk_two_group = sp_utilities.wrap_mpi_bcast(
        chunk_two_group, Blockdata["main_node"]
    )
    for index_of_element in range(len(chunk_one)):
        chunk_dict[chunk_one[index_of_element]] = 0
        group_dict[chunk_one[index_of_element]] = chunk_one_group[index_of_element]
    for index_of_element in range(len(chunk_two)):
        chunk_dict[chunk_two[index_of_element]] = 1
        group_dict[chunk_two[index_of_element]] = chunk_two_group[index_of_element]
    Tracker["chunk_dict"] = chunk_dict
    Tracker["P_chunk_0"] = old_div(len(chunk_one), float(total_stack))
    Tracker["P_chunk_1"] = old_div(len(chunk_two), float(total_stack))
    if Blockdata["myid"] == Blockdata["main_node"]:
        chunk_ids = []
        group_ids = []
        partids = sp_utilities.read_text_file(
            os.path.join(Tracker["constants"]["masterdir"], "indexes.txt"), -1
        )
        partids = partids[0]
        Tracker["constants"]["total_stack"] = len(partids)
        params = sp_utilities.read_text_file(
            os.path.join(
                Tracker["constants"]["masterdir"], "sparx_refinement_params.txt"
            ),
            -1,
        )
        for index_of_particle in range(len(partids)):
            chunk_ids.append(chunk_dict[partids[index_of_particle]])
            group_ids.append(group_dict[partids[index_of_particle]])
        refinement_params = [
            params[0],
            params[1],
            params[2],
            params[3],
            params[4],
            chunk_ids,
            group_ids,
            params[7],
        ]
        sp_utilities.write_text_file(
            refinement_params,
            os.path.join(
                Tracker["constants"]["masterdir"], "refinement_parameters.txt"
            ),
        )
    else:
        Tracker["constants"]["total_stack"] = 0
    Tracker["constants"]["total_stack"] = sp_utilities.bcast_number_to_all(
        Tracker["constants"]["total_stack"], Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    Tracker["total_stack"] = Tracker["constants"]["total_stack"]
    Tracker["constants"]["partstack"] = os.path.join(
        Tracker["constants"]["masterdir"], "refinement_parameters.txt"
    )
    total_stack = Tracker["constants"]["total_stack"]
    Tracker["currentres"] = old_div(
        float(Tracker["constants"]["fsc05"]), float(Tracker["constants"]["nxinit"])
    )
    Tracker["bckgnoise"] = os.path.join(
        Tracker["constants"]["masterdir"], "bckgnoise.hdf"
    )
    ###
    if Tracker["constants"]["minimum_grp_size"] == -1:
        Tracker["constants"]["minimum_grp_size"] = old_div(
            Tracker["constants"]["total_stack"], Tracker["constants"]["img_per_grp"]
        ) * (old_div(100, Blockdata["symclass"].nsym))
    # else:
    #    if Tracker["constants"]["minimum_grp_size"] < Tracker["constants"]["total_stack"]//Tracker["constants"]["img_per_grp"]*(100//Blockdata["symclass"].nsym):
    #        Tracker["constants"]["minimum_grp_size"] = Tracker["constants"]["total_stack"]//Tracker["constants"]["img_per_grp"]*(100//Blockdata["symclass"].nsym)
    #        if(Blockdata["myid"] == Blockdata["main_node"]):
    #            log_main.add("User provided minimum_grp_size is replaced by %d"%Tracker["constants"]["minimum_grp_size"])
    if Tracker["constants"]["minimum_grp_size"] > Tracker["constants"]["img_per_grp"]:
        sp_global_def.ERROR(
            "User provided img_per_grp is smaller than minimum_grp_size",
            "get_input_from_sparx_ref3d",
            1,
            Blockdata["myid"],
        )
    # Now copy oldparamstruture
    copy_oldparamstructure_from_meridien_MPI(selected_iter, log_main)
    return import_from_sparx_refinement


def get_input_from_datastack(log_main):  # Case three
    global Tracker, Blockdata
    import_from_data_stack = 1

    if Blockdata["myid"] == Blockdata["main_node"]:
        image = sp_utilities.get_im(Tracker["constants"]["orgstack"], 0)
        Tracker["constants"]["nnxo"] = image.get_xsize()
        if Tracker["nxinit"] > Tracker["constants"]["nnxo"]:
            sp_global_def.ERROR(
                "Image size less than minimum permitted $d" % Tracker["nxinit"],
                "get_input_from_datastack",
                1,
                Blockdata["myid"],
            )
            nnxo = -1
        else:
            if Tracker["constants"]["CTF"]:
                ictf = image.get_attr("ctf")
                Tracker["constants"]["pixel_size"] = ictf.apix
            else:
                Tracker["constants"]["pixel_size"] = 1.0
                del image
    else:
        Tracker["constants"]["nnxo"] = 0
        Tracker["constants"]["pixel_size"] = 1.0
    Tracker["constants"]["nnxo"] = sp_utilities.bcast_number_to_all(
        Tracker["constants"]["nnxo"], source_node=Blockdata["main_node"]
    )
    if Tracker["constants"]["nnxo"] < 0:
        sp_global_def.ERROR(
            "Image size is negative",
            "get_input_from_datastack",
            1,
            Blockdata["main_node"],
        )
    Tracker["constants"]["pixel_size"] = sp_utilities.bcast_number_to_all(
        Tracker["constants"]["pixel_size"], source_node=Blockdata["main_node"]
    )
    if Tracker["constants"]["radius"] < 1:
        Tracker["constants"]["radius"] = old_div(Tracker["constants"]["nnxo"], 2) - 2
    elif (2 * Tracker["constants"]["radius"] + 2) > Tracker["constants"]["nnxo"]:
        sp_global_def.ERROR(
            "Particle radius set too large!",
            "get_input_from_datastack",
            1,
            Blockdata["myid"],
        )
    if Blockdata["myid"] == Blockdata["main_node"]:
        total_stack = EMAN2_cppwrap.EMUtil.get_image_count(
            Tracker["constants"]["orgstack"]
        )
    else:
        total_stack = 0
    total_stack = sp_utilities.bcast_number_to_all(total_stack, Blockdata["main_node"])
    # randomly assign two subsets
    Tracker["constants"]["total_stack"] = total_stack
    Tracker["constants"]["chunk_0"] = os.path.join(
        Tracker["constants"]["masterdir"], "chunk_0.txt"
    )
    Tracker["constants"]["chunk_1"] = os.path.join(
        Tracker["constants"]["masterdir"], "chunk_1.txt"
    )
    Tracker["constants"]["partstack"] = os.path.join(
        Tracker["constants"]["masterdir"], "refinement_parameters.txt"
    )
    Tracker["previous_parstack"] = os.path.join(
        Tracker["constants"]["masterdir"], "refinement_parameters.txt"
    )  #

    if Tracker["constants"]["minimum_grp_size"] == -1:
        Tracker["constants"]["minimum_grp_size"] = old_div(
            Tracker["constants"]["total_stack"], Tracker["constants"]["img_per_grp"]
        ) * (old_div(100, Blockdata["symclass"].nsym))
    """Multiline Comment0"""
    if Tracker["constants"]["minimum_grp_size"] > Tracker["constants"]["img_per_grp"]:
        sp_global_def.ERROR(
            "User provided img_per_grp is smaller than minimum_grp_size",
            "get_input_from_sparx_ref3d",
            1,
            Blockdata["myid"],
        )

    ###
    Tracker["refang"], Tracker["rshifts"], Tracker["delta"] = None, None, None
    Tracker["avgnorm"] = 1.0
    chunk_dict = {}
    chunk_list = []
    if Blockdata["myid"] == Blockdata["main_node"]:
        chunk_dict = {}
        tlist = list(range(total_stack))
        sp_utilities.write_text_file(
            tlist, os.path.join(Tracker["constants"]["masterdir"], "indexes.txt")
        )
        random.shuffle(tlist)
        chunk_one = tlist[0 : old_div(total_stack, 2)]
        chunk_two = tlist[old_div(total_stack, 2) :]
        chunk_one = sorted(chunk_one)
        chunk_two = sorted(chunk_two)
        sp_utilities.write_text_row(chunk_one, Tracker["constants"]["chunk_0"])
        sp_utilities.write_text_row(chunk_two, Tracker["constants"]["chunk_1"])
        for particle in chunk_one:
            chunk_dict[particle] = 0
        for particle in chunk_two:
            chunk_dict[particle] = 1
        xform_proj_list = EMAN2_cppwrap.EMUtil.get_all_attributes(
            Tracker["constants"]["orgstack"], "xform.projection"
        )
        for index_of_particle in range(len(xform_proj_list)):
            dp = xform_proj_list[index_of_particle].get_params("spider")
            xform_proj_list[index_of_particle] = [
                dp["phi"],
                dp["theta"],
                dp["psi"],
                -dp["tx"],
                -dp["ty"],
                chunk_dict[index_of_particle],
            ]
        sp_utilities.write_text_row(xform_proj_list, Tracker["constants"]["partstack"])
    else:
        chunk_one = 0
        chunk_two = 0
    chunk_one = sp_utilities.wrap_mpi_bcast(chunk_one, Blockdata["main_node"])
    chunk_two = sp_utilities.wrap_mpi_bcast(chunk_two, Blockdata["main_node"])
    for element in chunk_one:
        chunk_dict[element] = 0
    for element in chunk_two:
        chunk_dict[element] = 1
    chunk_list = [chunk_one, chunk_two]
    Tracker["chunk_dict"] = chunk_dict
    Tracker["P_chunk_0"] = old_div(len(chunk_one), float(total_stack))
    Tracker["P_chunk_1"] = old_div(len(chunk_two), float(total_stack))

    # Reconstruction to determine the resolution in orignal data size
    Tracker["nxinit"] = Tracker["constants"]["nnxo"]
    Tracker["shrinkage"] = old_div(
        float(Tracker["nxinit"]), float(Tracker["constants"]["nnxo"])
    )
    Tracker["bckgnoise"] = None
    temp = sp_utilities.model_blank(
        Tracker["constants"]["nnxo"], Tracker["constants"]["nnxo"]
    )
    nny = temp.get_ysize()
    Blockdata["bckgnoise"] = [1.0] * nny  # set for initial recon3D of data from stack
    Tracker["fuse_freq"] = int(
        old_div(
            Tracker["constants"]["pixel_size"] * Tracker["constants"]["nnxo"],
            Tracker["constants"]["fuse_freq"],
        )
        + 0.5
    )
    Tracker["directory"] = Tracker["constants"]["masterdir"]
    if Tracker["constants"]["nxinit"] < 0:
        Tracker["nxinit_refinement"] = Tracker["constants"]["nnxo"]
    else:
        Tracker["nxinit_refinement"] = Tracker["constants"]["nxinit"]

    for procid in range(2):
        data = get_shrink_data_sorting(
            os.path.join(Tracker["constants"]["masterdir"], "chunk_%01d.txt" % procid),
            Tracker["constants"]["partstack"],
        )
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        do3d_sorting(
            procid, data, myid=Blockdata["myid"], mpi_comm=mpi.MPI_COMM_WORLD
        )  # 1
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    if Blockdata["no_of_groups"] == 1:
        if Blockdata["myid"] == Blockdata["main_node"]:
            tvol0 = sp_utilities.get_im(
                os.path.join(Tracker["directory"], "tempdir", "tvol_0.hdf")
            )
            tweight0 = sp_utilities.get_im(
                os.path.join(Tracker["directory"], "tempdir", "tweight_0.hdf")
            )
            tvol1 = sp_utilities.get_im(
                os.path.join(Tracker["directory"], "tempdir", "tvol_1.hdf")
            )
            tweight1 = sp_utilities.get_im(
                os.path.join(Tracker["directory"], "tempdir", "tweight_1.hdf")
            )
            EMAN2_cppwrap.Util.fuse_low_freq(
                tvol0, tvol1, tweight0, tweight1, 2 * Tracker["fuse_freq"]
            )
            shrank0 = stepone(tvol0, tweight0)
            shrank1 = stepone(tvol1, tweight1)
            #  Note shrank volumes are Fourier uncentered.
            cfsc = sp_statistics.fsc(shrank0, shrank1)[1]
            del shrank0, shrank1
            if Tracker["nxinit"] < Tracker["constants"]["nnxo"]:
                cfsc = cfsc[: Tracker["nxinit"]]
            for i in range(len(cfsc), old_div(Tracker["constants"]["nnxo"], 2) + 1):
                cfsc.append(0.0)
            sp_utilities.write_text_row(
                cfsc, os.path.join(Tracker["directory"], "fsc_global.txt")
            )
            lcfsc = len(cfsc)
            fsc05 = 0
            fsc143 = 0
            for ifreq in range(len(cfsc)):
                if cfsc[ifreq] < 0.5:
                    break
            fsc05 = ifreq - 1
            for ifreq in range(len(cfsc)):
                if cfsc[ifreq] < 0.143:
                    break
            fsc143 = ifreq - 1
            Tracker["constants"]["fsc143"] = fsc143
            Tracker["constants"]["fsc05"] = fsc05
        Tracker = sp_utilities.wrap_mpi_bcast(
            Tracker, Blockdata["nodes"][0], communicator=mpi.MPI_COMM_WORLD
        )
    else:
        if (
            Blockdata["myid"] == Blockdata["nodes"][1]
        ):  # It has to be 1 to avoid problem with tvol1 not closed on the disk
            tvol0 = sp_utilities.get_im(
                os.path.join(Tracker["directory"], "tempdir", "tvol_0.hdf")
            )
            tweight0 = sp_utilities.get_im(
                os.path.join(Tracker["directory"], "tempdir", "tweight_0.hdf")
            )
            tvol1 = sp_utilities.get_im(
                os.path.join(Tracker["directory"], "tempdir", "tvol_1.hdf")
            )
            tweight1 = sp_utilities.get_im(
                os.path.join(Tracker["directory"], "tempdir", "tweight_1.hdf")
            )
            EMAN2_cppwrap.Util.fuse_low_freq(
                tvol0, tvol1, tweight0, tweight1, 2 * Tracker["fuse_freq"]
            )
            tag = 7007
            sp_utilities.send_EMData(
                tvol1, Blockdata["nodes"][0], tag, mpi.MPI_COMM_WORLD
            )
            sp_utilities.send_EMData(
                tweight1, Blockdata["nodes"][0], tag, mpi.MPI_COMM_WORLD
            )
            shrank0 = stepone(tvol0, tweight0)
            sp_utilities.send_EMData(
                shrank0, Blockdata["nodes"][0], tag, mpi.MPI_COMM_WORLD
            )
            del shrank0
            lcfsc = 0

        elif Blockdata["myid"] == Blockdata["nodes"][0]:
            tag = 7007
            tvol1 = sp_utilities.recv_EMData(
                Blockdata["nodes"][1], tag, mpi.MPI_COMM_WORLD
            )
            tweight1 = sp_utilities.recv_EMData(
                Blockdata["nodes"][1], tag, mpi.MPI_COMM_WORLD
            )
            tvol1.set_attr_dict(
                {"is_complex": 1, "is_fftodd": 1, "is_complex_ri": 1, "is_fftpad": 1}
            )
            shrank1 = stepone(tvol1, tweight1)
            #  Get shrank volume, do fsc, send it to all
            shrank0 = sp_utilities.recv_EMData(
                Blockdata["nodes"][1], tag, mpi.MPI_COMM_WORLD
            )
            #  Note shrank volumes are Fourier uncentered.
            cfsc = sp_statistics.fsc(shrank0, shrank1)[1]
            sp_utilities.write_text_row(
                cfsc, os.path.join(Tracker["directory"], "fsc_global.txt")
            )
            del shrank0, shrank1
            if Tracker["nxinit"] < Tracker["constants"]["nnxo"]:
                cfsc = cfsc[: Tracker["nxinit"]]
                for i in range(len(cfsc), old_div(Tracker["constants"]["nnxo"], 2) + 1):
                    cfsc.append(0.0)
            lcfsc = len(cfsc)
            fsc05 = 0
            fsc143 = 0
            for ifreq in range(len(cfsc)):
                if cfsc[ifreq] < 0.5:
                    break
            fsc05 = ifreq - 1
            for ifreq in range(len(cfsc)):
                if cfsc[ifreq] < 0.143:
                    break
            fsc143 = ifreq - 1
            Tracker["constants"]["fsc143"] = fsc143
            Tracker["constants"]["fsc05"] = fsc05
        Tracker = sp_utilities.wrap_mpi_bcast(
            Tracker, Blockdata["nodes"][0], communicator=mpi.MPI_COMM_WORLD
        )
    return import_from_data_stack


####
### functions for faked rec3d from subsets
###


def do3d_sorting_groups_nofsc_smearing_iter(srdata, partial_rec3d, iteration):
    global Tracker, Blockdata
    keepgoing = 1
    if Blockdata["myid"] == Blockdata["last_node"]:
        if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
            os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
        try:
            fout = open(os.path.join(Tracker["directory"], "freq_cutoff.json"), "r")
            freq_cutoff_dict = sp_utilities.convert_json_fromunicode(json.load(fout))
            fout.close()
        except:
            freq_cutoff_dict = 0
    else:
        freq_cutoff_dict = 0
    freq_cutoff_dict = sp_utilities.wrap_mpi_bcast(
        freq_cutoff_dict, Blockdata["last_node"], mpi.MPI_COMM_WORLD
    )

    for index_of_groups in range(Tracker["number_of_groups"]):
        if partial_rec3d:
            tvol, tweight, trol = recons3d_trl_struct_group_nofsc_shifted_data_partial_MPI(
                Blockdata["myid"],
                Blockdata["last_node"],
                Blockdata["nproc"],
                srdata,
                index_of_groups,
                os.path.join(
                    Tracker["directory"], "tempdir", "trol_2_%d.hdf" % index_of_groups
                ),
                os.path.join(
                    Tracker["directory"], "tempdir", "tvol_2_%d.hdf" % index_of_groups
                ),
                os.path.join(
                    Tracker["directory"],
                    "tempdir",
                    "tweight_2_%d.hdf" % index_of_groups,
                ),
                None,
                Tracker["constants"]["CTF"],
                (2 * Tracker["nxinit"] + 3),
                Tracker["nosmearing"],
            )
        else:
            tvol, tweight, trol = recons3d_trl_struct_group_nofsc_shifted_data_MPI(
                Blockdata["myid"],
                Blockdata["last_node"],
                srdata,
                index_of_groups,
                None,
                Tracker["constants"]["CTF"],
                (2 * Tracker["nxinit"] + 3),
                Tracker["nosmearing"],
            )

        if Blockdata["myid"] == Blockdata["last_node"]:
            tvol.set_attr("is_complex", 0)
            tvol.write_image(
                os.path.join(
                    Tracker["directory"], "tempdir", "tvol_2_%d.hdf" % index_of_groups
                )
            )
            tweight.write_image(
                os.path.join(
                    Tracker["directory"],
                    "tempdir",
                    "tweight_2_%d.hdf" % index_of_groups,
                )
            )
            trol.write_image(
                os.path.join(
                    Tracker["directory"], "tempdir", "trol_2_%d.hdf" % index_of_groups
                )
            )
            del tvol
            del tweight
            del trol
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    Tracker["fsc143"] = 0
    Tracker["fsc05"] = 0
    Tracker["maxfrad"] = old_div(Tracker["nxinit"], 2)
    if Blockdata["no_of_groups"] > 1:
        # new starts
        sub_main_node_list = [-1 for i in range(Blockdata["no_of_groups"])]
        for index_of_colors in range(Blockdata["no_of_groups"]):
            for iproc in range(Blockdata["nproc"] - 1):
                if Blockdata["myid"] == iproc:
                    if (
                        Blockdata["color"] == index_of_colors
                        and Blockdata["myid_on_node"] == 0
                    ):
                        sub_main_node_list[index_of_colors] = Blockdata["myid"]
                    sp_utilities.wrap_mpi_send(
                        sub_main_node_list, Blockdata["last_node"], mpi.MPI_COMM_WORLD
                    )
                if Blockdata["myid"] == Blockdata["last_node"]:
                    dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
                    for im in range(len(dummy)):
                        if dummy[im] > -1:
                            sub_main_node_list[im] = dummy[im]
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        sp_utilities.wrap_mpi_bcast(
            sub_main_node_list, Blockdata["last_node"], mpi.MPI_COMM_WORLD
        )
        if Tracker["number_of_groups"] % Blockdata["no_of_groups"] == 0:
            nbig_loop = old_div(Tracker["number_of_groups"], Blockdata["no_of_groups"])
        else:
            nbig_loop = (
                old_div(Tracker["number_of_groups"], Blockdata["no_of_groups"]) + 1
            )

        big_loop_colors = [[] for i in range(nbig_loop)]
        big_loop_groups = [[] for i in range(nbig_loop)]
        nc = 0
        while nc < Tracker["number_of_groups"]:
            im = old_div(nc, Blockdata["no_of_groups"])
            jm = nc % Blockdata["no_of_groups"]
            big_loop_colors[im].append(jm)
            big_loop_groups[im].append(nc)
            nc += 1
        for iloop in range(nbig_loop):
            for im in range(len(big_loop_colors[iloop])):
                index_of_group = big_loop_groups[iloop][im]
                index_of_colors = big_loop_colors[iloop][im]

                if Blockdata["myid"] == Blockdata["last_node"]:
                    tvol2 = sp_utilities.get_im(
                        os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")
                        % index_of_group
                    )
                    tweight2 = sp_utilities.get_im(
                        os.path.join(
                            Tracker["directory"], "tempdir", "tweight_2_%d.hdf"
                        )
                        % index_of_group
                    )
                    treg2 = sp_utilities.get_im(
                        os.path.join(
                            Tracker["directory"],
                            "tempdir",
                            "trol_2_%d.hdf" % index_of_group,
                        )
                    )
                    tvol2.set_attr_dict(
                        {
                            "is_complex": 1,
                            "is_fftodd": 1,
                            "is_complex_ri": 1,
                            "is_fftpad": 1,
                        }
                    )
                    tag = 7007
                    sp_utilities.send_EMData(
                        tvol2,
                        sub_main_node_list[index_of_colors],
                        tag,
                        mpi.MPI_COMM_WORLD,
                    )
                    sp_utilities.send_EMData(
                        tweight2,
                        sub_main_node_list[index_of_colors],
                        tag,
                        mpi.MPI_COMM_WORLD,
                    )
                    sp_utilities.send_EMData(
                        treg2,
                        sub_main_node_list[index_of_colors],
                        tag,
                        mpi.MPI_COMM_WORLD,
                    )

                elif Blockdata["myid"] == sub_main_node_list[index_of_colors]:
                    tag = 7007
                    tvol2 = sp_utilities.recv_EMData(
                        Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                    tweight2 = sp_utilities.recv_EMData(
                        Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                    treg2 = sp_utilities.recv_EMData(
                        Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            for im in range(len(big_loop_colors[iloop])):
                index_of_group = big_loop_groups[iloop][im]
                index_of_colors = big_loop_colors[iloop][im]
                try:
                    Tracker["freq_fsc143_cutoff"] = freq_cutoff_dict[
                        "Cluster_%03d.txt" % index_of_group
                    ]
                except:
                    pass
                if (
                    Blockdata["color"] == index_of_colors
                ):  # It has to be 1 to avoid problem with tvol1 not closed on the disk
                    if Blockdata["myid_on_node"] != 0:
                        tvol2 = sp_utilities.model_blank(1)
                        tweight2 = sp_utilities.model_blank(1)
                        treg2 = sp_utilities.model_blank(1)
                    tvol2 = steptwo_mpi_filter(
                        tvol2,
                        tweight2,
                        treg2,
                        None,
                        Tracker["freq_fsc143_cutoff"],
                        0.01,
                        False,
                        color=index_of_colors,
                    )  # has to be False!!!
                    del tweight2, treg2
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            for im in range(len(big_loop_colors[iloop])):
                index_of_group = big_loop_groups[iloop][im]
                index_of_colors = big_loop_colors[iloop][im]
                if (Blockdata["color"] == index_of_colors) and (
                    Blockdata["myid_on_node"] == 0
                ):
                    tag = 7007
                    sp_utilities.send_EMData(
                        tvol2, Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                elif Blockdata["myid"] == Blockdata["last_node"]:
                    tag = 7007
                    tvol2 = sp_utilities.recv_EMData(
                        sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD
                    )
                    tvol2.write_image(
                        os.path.join(
                            Tracker["directory"],
                            "vol_grp%03d_iter%03d.hdf" % (index_of_group, iteration),
                        )
                    )
                    del tvol2
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    else:  # loop over all groups for single node
        for index_of_group in range(Tracker["number_of_groups"]):
            if Blockdata["myid_on_node"] == 0:
                tvol2 = sp_utilities.get_im(
                    os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")
                    % index_of_group
                )
                tweight2 = sp_utilities.get_im(
                    os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")
                    % index_of_group
                )
                treg2 = sp_utilities.get_im(
                    os.path.join(
                        Tracker["directory"],
                        "tempdir",
                        "trol_2_%d.hdf" % index_of_group,
                    )
                )
                tvol2.set_attr_dict(
                    {
                        "is_complex": 1,
                        "is_fftodd": 1,
                        "is_complex_ri": 1,
                        "is_fftpad": 1,
                    }
                )
            else:
                tvol2 = sp_utilities.model_blank(1)
                tweight2 = sp_utilities.model_blank(1)
                treg2 = sp_utilities.model_blank(1)
            try:
                Tracker["freq_fsc143_cutoff"] = freq_cutoff_dict[
                    "Cluster_%03d.txt" % index_of_group
                ]
            except:
                pass
            tvol2 = steptwo_mpi_filter(
                tvol2, tweight2, treg2, None, Tracker["freq_fsc143_cutoff"], 0.01, False
            )  # has to be False!!!
            del tweight2, treg2
            if Blockdata["myid_on_node"] == 0:
                tvol2.write_image(
                    os.path.join(
                        Tracker["directory"],
                        "vol_grp%03d_iter%03d.hdf" % (index_of_group, iteration),
                    )
                )
                del tvol2
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    keepgoing = sp_utilities.bcast_number_to_all(
        keepgoing, source_node=Blockdata["main_node"], mpi_comm=mpi.MPI_COMM_WORLD
    )  # always check
    Tracker = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"])
    if not keepgoing:
        sp_global_def.ERROR(
            "do3d_sorting_groups_trl_iter  %s"
            % os.path.join(Tracker["directory"], "tempdir"),
            "do3d_sorting_groups_trl_iter",
            1,
            Blockdata["myid"],
        )
    return


### nofsc insertion #1


def recons3d_trl_struct_group_nofsc_shifted_data_partial_MPI(
    myid,
    main_node,
    nproc,
    prjlist,
    group_ID,
    refvol_file,
    fftvol_file,
    weight_file,
    mpi_comm=None,
    CTF=True,
    target_size=-1,
    nosmearing=False,
):
    """
        partial rec3d for re-assigned images
        reconstructor nn4_ctfws
    """
    if mpi_comm == None:
        mpi_comm = mpi.MPI_COMM_WORLD
    if CTF:
        do_ctf = 1
    else:
        do_ctf = 0
    if not os.path.exists(refvol_file):
        sp_global_def.ERROR(
            "refvol does not exist",
            "recons3d_trl_struct_group_nofsc_shifted_data_partial_MPI",
            1,
            myid,
        )
    if not os.path.exists(fftvol_file):
        sp_global_def.ERROR(
            "fftvol does not exist",
            "recons3d_trl_struct_group_nofsc_shifted_data_partial_MPI",
            1,
            myid,
        )
    if not os.path.exists(weight_file):
        sp_global_def.ERROR(
            "weight does not exist",
            "recons3d_trl_struct_group_nofsc_shifted_data_partial_MPI",
            1,
            myid,
        )

    # refvol
    if myid == main_node:
        target_size = sp_utilities.get_im(refvol_file).get_xsize()
    else:
        target_size = 0
    target_size = sp_utilities.bcast_number_to_all(target_size, main_node, mpi_comm)
    refvol = sp_utilities.model_blank(target_size)  # set to zero

    # fftvol
    if myid == main_node:
        fftvol = sp_utilities.get_im(fftvol_file)
        fftvol.set_attr_dict(
            {"is_complex": 1, "is_fftodd": 1, "is_complex_ri": 1, "is_fftpad": 1}
        )
        fftvol = old_div(fftvol, float(nproc))
        # Util.mult_scalar(fftvol, 1./float(Blockdata["nproc"]))
        nxfft = fftvol.get_xsize()
        nyfft = fftvol.get_ysize()
        nzfft = fftvol.get_zsize()
    else:
        nxfft = 0
        nyfft = 0
        nzfft = 0
    nxfft = sp_utilities.bcast_number_to_all(nxfft, main_node, mpi_comm)
    nyfft = sp_utilities.bcast_number_to_all(nyfft, main_node, mpi_comm)
    nzfft = sp_utilities.bcast_number_to_all(nzfft, main_node, mpi_comm)
    if myid != main_node:
        fftvol = sp_utilities.model_blank(nxfft, nyfft, nzfft)
    sp_utilities.bcast_EMData_to_all(fftvol, myid, main_node)

    # weight
    if myid == main_node:
        weight = sp_utilities.get_im(weight_file)
        weight = old_div(weight, float(nproc))
        # Util.mult_scalar(weight, 1./float(Blockdata["nproc"]))
        nxweight = weight.get_xsize()
        nyweight = weight.get_ysize()
        nzweight = weight.get_zsize()
    else:
        nxweight = 0
        nyweight = 0
        nzweight = 0
    nxweight = sp_utilities.bcast_number_to_all(nxweight, main_node, mpi_comm)
    nyweight = sp_utilities.bcast_number_to_all(nyweight, main_node, mpi_comm)
    nzweight = sp_utilities.bcast_number_to_all(nzweight, main_node, mpi_comm)
    if myid != main_node:
        weight = sp_utilities.model_blank(nxweight, nyweight, nzweight)
    sp_utilities.bcast_EMData_to_all(weight, myid, main_node)
    params = {
        "size": target_size,
        "npad": 2,
        "snr": 1.0,
        "sign": 1,
        "symmetry": "c1",
        "refvol": refvol,
        "fftvol": fftvol,
        "weight": weight,
        "do_ctf": do_ctf,
    }
    r = EMAN2_cppwrap.Reconstructors.get("nn4_ctfws", params)
    r.setup()
    if nosmearing:
        nnx = prjlist[0].get_xsize()
        nny = prjlist[0].get_ysize()
    else:
        nnx = prjlist[0][0].get_xsize()
        nny = prjlist[0][0].get_ysize()
    for im in range(len(prjlist)):
        if nosmearing:
            current_group_ID = prjlist[im].get_attr("group")
            previous_group_ID = prjlist[im].get_attr("previous_group")
            if current_group_ID != previous_group_ID:
                if current_group_ID == group_ID:
                    flag = 1.0
                    [phi, theta, psi, s2x, s2y] = sp_utilities.get_params_proj(
                        prjlist[im], xform="xform.projection"
                    )
                    r.insert_slice(
                        prjlist[im],
                        EMAN2_cppwrap.Transform(
                            {"type": "spider", "phi": phi, "theta": theta, "psi": psi}
                        ),
                        flag,
                    )
                if previous_group_ID == group_ID:
                    flag = -1.0
                    [phi, theta, psi, s2x, s2y] = sp_utilities.get_params_proj(
                        prjlist[im], xform="xform.projection"
                    )
                    r.insert_slice(
                        prjlist[im],
                        EMAN2_cppwrap.Transform(
                            {"type": "spider", "phi": phi, "theta": theta, "psi": psi}
                        ),
                        flag,
                    )
        else:
            current_group_ID = prjlist[im][0].get_attr("group")
            previous_group_ID = prjlist[im][0].get_attr("previous_group")
            if current_group_ID != previous_group_ID:
                if current_group_ID == group_ID:
                    flag = 1.0
                    for jm in range(len(prjlist[im])):
                        [phi, theta, psi, s2x, s2y] = sp_utilities.get_params_proj(
                            prjlist[im][jm], xform="xform.projection"
                        )
                        r.insert_slice(
                            prjlist[im][jm],
                            EMAN2_cppwrap.Transform(
                                {
                                    "type": "spider",
                                    "phi": phi,
                                    "theta": theta,
                                    "psi": psi,
                                }
                            ),
                            prjlist[im][jm].get_attr("wprob") * flag,
                        )
                if previous_group_ID == group_ID:
                    flag = -1.0
                    for jm in range(len(prjlist[im])):
                        [phi, theta, psi, s2x, s2y] = sp_utilities.get_params_proj(
                            prjlist[im][jm], xform="xform.projection"
                        )
                        r.insert_slice(
                            prjlist[im][jm],
                            EMAN2_cppwrap.Transform(
                                {
                                    "type": "spider",
                                    "phi": phi,
                                    "theta": theta,
                                    "psi": psi,
                                }
                            ),
                            prjlist[im][jm].get_attr("wprob") * flag,
                        )
    sp_utilities.reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
    sp_utilities.reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
    if myid == main_node:
        dummy = r.finish(True)
    mpi.mpi_barrier(mpi_comm)
    if myid == main_node:
        return fftvol, weight, refvol
    else:
        del fftvol
        del weight
        del refvol
        return None, None, None


### insertion 2


def recons3d_trl_struct_group_nofsc_shifted_data_MPI(
    myid,
    main_node,
    prjlist,
    group_ID,
    mpi_comm=None,
    CTF=True,
    target_size=-1,
    nosmearing=False,
):
    """
      rec3d for pre-shifted data list
      reconstructor nn4_ctfw
    """
    if mpi_comm == None:
        mpi_comm = mpi.MPI_COMM_WORLD
    refvol = sp_utilities.model_blank(target_size)
    refvol.set_attr("fudge", 1.0)
    if CTF:
        do_ctf = 1
    else:
        do_ctf = 0
    fftvol = EMAN2_cppwrap.EMData()
    weight = EMAN2_cppwrap.EMData()
    params = {
        "size": target_size,
        "npad": 2,
        "snr": 1.0,
        "sign": 1,
        "symmetry": "c1",
        "refvol": refvol,
        "fftvol": fftvol,
        "weight": weight,
        "do_ctf": do_ctf,
    }
    r = EMAN2_cppwrap.Reconstructors.get("nn4_ctfw", params)
    r.setup()
    if nosmearing:
        nnx = prjlist[0].get_xsize()
        nny = prjlist[0].get_ysize()
    else:
        nnx = prjlist[0][0].get_xsize()
        nny = prjlist[0][0].get_ysize()
    for im in range(len(prjlist)):
        if nosmearing:
            if prjlist[im].get_attr("group") == group_ID:
                [phi, theta, psi, s2x, s2y] = sp_utilities.get_params_proj(
                    prjlist[im], xform="xform.projection"
                )
                r.insert_slice(
                    prjlist[im],
                    EMAN2_cppwrap.Transform(
                        {"type": "spider", "phi": phi, "theta": theta, "psi": psi}
                    ),
                    1.0,
                )
        else:
            if prjlist[im][0].get_attr("group") == group_ID:
                for jm in range(len(prjlist[im])):
                    [phi, theta, psi, s2x, s2y] = sp_utilities.get_params_proj(
                        prjlist[im][jm], xform="xform.projection"
                    )
                    r.insert_slice(
                        prjlist[im][jm],
                        EMAN2_cppwrap.Transform(
                            {"type": "spider", "phi": phi, "theta": theta, "psi": psi}
                        ),
                        prjlist[im][jm].get_attr("wprob"),
                    )
    sp_utilities.reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
    sp_utilities.reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
    if myid == main_node:
        dummy = r.finish(True)
    mpi.mpi_barrier(mpi_comm)
    if myid == main_node:
        return fftvol, weight, refvol
    else:
        del fftvol
        del weight
        del refvol
        return None, None, None


###end of nofsc


def recons3d_trl_struct_group_MPI(
    myid,
    main_node,
    prjlist,
    random_subset,
    group_ID,
    paramstructure,
    norm_per_particle=None,
    upweighted=True,
    mpi_comm=None,
    CTF=True,
    target_size=-1,
    nosmearing=False,
):
    """
        recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
        Input
            list_of_prjlist: list of lists of projections to be included in the reconstruction
    """
    global Tracker, Blockdata
    if mpi_comm == None:
        mpi_comm = mpi.MPI_COMM_WORLD

    refvol = sp_utilities.model_blank(target_size)
    refvol.set_attr("fudge", 1.0)
    if CTF:
        do_ctf = 1
    else:
        do_ctf = 0
    fftvol = EMAN2_cppwrap.EMData()
    weight = EMAN2_cppwrap.EMData()
    params = {
        "size": target_size,
        "npad": 2,
        "snr": 1.0,
        "sign": 1,
        "symmetry": "c1",
        "refvol": refvol,
        "fftvol": fftvol,
        "weight": weight,
        "do_ctf": do_ctf,
    }
    r = EMAN2_cppwrap.Reconstructors.get("nn4_ctfw", params)
    r.setup()
    if norm_per_particle == None:
        norm_per_particle = len(prjlist) * [1.0]
    if not nosmearing:
        delta = Tracker["delta"]
        refang = Tracker["refang"]
        rshifts_shrank = copy.deepcopy(Tracker["rshifts"])
        nshifts = len(rshifts_shrank)
        for im in range(len(rshifts_shrank)):
            rshifts_shrank[im][0] *= old_div(
                float(Tracker["nxinit"]), float(Tracker["constants"]["nnxo"])
            )
            rshifts_shrank[im][1] *= old_div(
                float(Tracker["nxinit"]), float(Tracker["constants"]["nnxo"])
            )
    nnx = prjlist[0].get_xsize()
    nny = prjlist[0].get_ysize()
    for im in range(len(prjlist)):
        if not nosmearing:
            avgnorm = Tracker["avgnorm"][prjlist[im].get_attr("chunk_id")]
        #  parse projection structure, generate three lists:
        #  [ipsi+iang], [ishift], [probability]
        #  Number of orientations for a given image
        if prjlist[im].get_attr("group") == group_ID:
            if random_subset == 2:
                if nosmearing:
                    bckgn = prjlist[im].get_attr("bckgnoise")
                    ct = prjlist[im].get_attr("ctf")
                    prjlist[im].set_attr_dict({"bckgnoise": bckgn, "ctf": ct})
                    [phi, theta, psi, s2x, s2y] = sp_utilities.get_params_proj(
                        prjlist[im], xform="xform.projection"
                    )
                    r.insert_slice(
                        prjlist[im],
                        EMAN2_cppwrap.Transform(
                            {"type": "spider", "phi": phi, "theta": theta, "psi": psi}
                        ),
                        1.0,
                    )
                else:
                    if Tracker["constants"]["nsmear"] <= 0.0:
                        numbor = len(paramstructure[im][2])
                    else:
                        numbor = 1
                    ipsiandiang = [
                        old_div(paramstructure[im][2][i][0], 1000)
                        for i in range(numbor)
                    ]
                    allshifts = [
                        paramstructure[im][2][i][0] % 1000 for i in range(numbor)
                    ]
                    probs = [paramstructure[im][2][i][1] for i in range(numbor)]
                    #  Find unique projection directions
                    tdir = list(set(ipsiandiang))
                    bckgn = prjlist[im].get_attr("bckgnoise")
                    ct = prjlist[im].get_attr("ctf")
                    #  For each unique projection direction:
                    data = [None] * nshifts
                    for ii in range(len(tdir)):
                        #  Find the number of times given projection direction appears on the list, it is the number of different shifts associated with it.
                        lshifts = sp_utilities.findall(tdir[ii], ipsiandiang)
                        toprab = 0.0
                        for ki in range(len(lshifts)):
                            toprab += probs[lshifts[ki]]
                        recdata = EMAN2_cppwrap.EMData(nny, nny, 1, False)
                        recdata.set_attr("is_complex", 0)
                        for ki in range(len(lshifts)):
                            lpt = allshifts[lshifts[ki]]
                            if data[lpt] == None:
                                data[lpt] = sp_fundamentals.fshift(
                                    prjlist[im],
                                    rshifts_shrank[lpt][0],
                                    rshifts_shrank[lpt][1],
                                )
                                data[lpt].set_attr("is_complex", 0)
                            EMAN2_cppwrap.Util.add_img(
                                recdata,
                                EMAN2_cppwrap.Util.mult_scalar(
                                    data[lpt], old_div(probs[lshifts[ki]], toprab)
                                ),
                            )
                        recdata.set_attr_dict(
                            {
                                "padffted": 1,
                                "is_fftpad": 1,
                                "is_fftodd": 0,
                                "is_complex_ri": 1,
                                "is_complex": 1,
                            }
                        )
                        if not upweighted:
                            recdata = sp_filter.filt_table(recdata, bckgn)
                        recdata.set_attr_dict({"bckgnoise": bckgn, "ctf": ct})
                        ipsi = tdir[ii] % 100000
                        iang = old_div(tdir[ii], 100000)
                        r.insert_slice(
                            recdata,
                            EMAN2_cppwrap.Transform(
                                {
                                    "type": "spider",
                                    "phi": refang[iang][0],
                                    "theta": refang[iang][1],
                                    "psi": refang[iang][2] + ipsi * delta,
                                }
                            ),
                            old_div(toprab * avgnorm, norm_per_particle[im]),
                        )
            else:
                if prjlist[im].get_attr("chunk_id") == random_subset:
                    if nosmearing:
                        bckgn = prjlist[im].get_attr("bckgnoise")
                        ct = prjlist[im].get_attr("ctf")
                        prjlist[im].set_attr_dict({"bckgnoise": bckgn, "ctf": ct})
                        [phi, theta, psi, s2x, s2y] = sp_utilities.get_params_proj(
                            prjlist[im], xform="xform.projection"
                        )
                        r.insert_slice(
                            prjlist[im],
                            EMAN2_cppwrap.Transform(
                                {
                                    "type": "spider",
                                    "phi": phi,
                                    "theta": theta,
                                    "psi": psi,
                                }
                            ),
                            1.0,
                        )
                    else:
                        if Tracker["constants"]["nsmear"] <= 0.0:
                            numbor = len(paramstructure[im][2])
                        else:
                            numbor = 1
                        ipsiandiang = [
                            old_div(paramstructure[im][2][i][0], 1000)
                            for i in range(numbor)
                        ]
                        allshifts = [
                            paramstructure[im][2][i][0] % 1000 for i in range(numbor)
                        ]
                        probs = [paramstructure[im][2][i][1] for i in range(numbor)]
                        #  Find unique projection directions
                        tdir = list(set(ipsiandiang))
                        bckgn = prjlist[im].get_attr("bckgnoise")
                        ct = prjlist[im].get_attr("ctf")
                        #  For each unique projection direction:
                        data = [None] * nshifts
                        for ii in range(len(tdir)):
                            #  Find the number of times given projection direction appears on the list, it is the number of different shifts associated with it.
                            lshifts = sp_utilities.findall(tdir[ii], ipsiandiang)
                            toprab = 0.0
                            for ki in range(len(lshifts)):
                                toprab += probs[lshifts[ki]]
                            recdata = EMAN2_cppwrap.EMData(nny, nny, 1, False)
                            recdata.set_attr("is_complex", 0)
                            for ki in range(len(lshifts)):
                                lpt = allshifts[lshifts[ki]]
                                if data[lpt] == None:
                                    data[lpt] = sp_fundamentals.fshift(
                                        prjlist[im],
                                        rshifts_shrank[lpt][0],
                                        rshifts_shrank[lpt][1],
                                    )
                                    data[lpt].set_attr("is_complex", 0)
                                EMAN2_cppwrap.Util.add_img(
                                    recdata,
                                    EMAN2_cppwrap.Util.mult_scalar(
                                        data[lpt], old_div(probs[lshifts[ki]], toprab)
                                    ),
                                )
                            recdata.set_attr_dict(
                                {
                                    "padffted": 1,
                                    "is_fftpad": 1,
                                    "is_fftodd": 0,
                                    "is_complex_ri": 1,
                                    "is_complex": 1,
                                }
                            )
                            if not upweighted:
                                recdata = sp_filter.filt_table(recdata, bckgn)
                            recdata.set_attr_dict({"bckgnoise": bckgn, "ctf": ct})
                            ipsi = tdir[ii] % 100000
                            iang = old_div(tdir[ii], 100000)
                            r.insert_slice(
                                recdata,
                                EMAN2_cppwrap.Transform(
                                    {
                                        "type": "spider",
                                        "phi": refang[iang][0],
                                        "theta": refang[iang][1],
                                        "psi": refang[iang][2] + ipsi * delta,
                                    }
                                ),
                                old_div(toprab * avgnorm, norm_per_particle[im]),
                            )
    #  clean stuff
    # if not nosmearing: del recdata, tdir, ipsiandiang, allshifts, probs
    sp_utilities.reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
    sp_utilities.reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
    if not nosmearing:
        del rshifts_shrank
    if myid == main_node:
        dummy = r.finish(True)
    mpi.mpi_barrier(mpi_comm)
    if myid == main_node:
        return fftvol, weight, refvol
    else:
        del fftvol
        del weight
        del refvol
        return None, None, None


####==========--------
#####  FSC rec3d


def do3d_sorting_groups_fsc_only_iter(
    data, paramstructure, norm_per_particle, iteration
):
    global Tracker, Blockdata
    # do resolution each time
    keepgoing = 1
    if Blockdata["myid"] == Blockdata["nodes"][0]:
        if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
            os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
    do3d_sorting_group_insertion_random_two_for_fsc(
        data, paramstructure, norm_per_particle
    )
    fsc143 = 0
    fsc05 = 0
    Tracker["fsc143"] = 0
    Tracker["fsc05"] = 0
    res_05 = Tracker["number_of_groups"] * [0]
    res_143 = Tracker["number_of_groups"] * [0]
    for index_of_colors in range(Blockdata["no_of_groups"]):
        group_start, group_end = MPI_volume_start_end(
            Tracker["number_of_groups"], Blockdata["no_of_groups"], index_of_colors
        )
        if (
            Blockdata["color"] == index_of_colors
        ):  # It has to be 1 to avoid problem with tvol1 not closed on the disk
            for index_of_group in range(group_start, group_end):
                if Blockdata["myid_on_node"] == 0:
                    tvol0 = sp_utilities.get_im(
                        os.path.join(Tracker["directory"], "tempdir", "tvol_0_0_%d.hdf")
                        % index_of_group
                    )
                    tweight0 = sp_utilities.get_im(
                        os.path.join(
                            Tracker["directory"], "tempdir", "tweight_0_0_%d.hdf"
                        )
                        % index_of_group
                    )
                    tvol1 = sp_utilities.get_im(
                        os.path.join(Tracker["directory"], "tempdir", "tvol_1_0_%d.hdf")
                        % index_of_group
                    )
                    tweight1 = sp_utilities.get_im(
                        os.path.join(
                            Tracker["directory"], "tempdir", "tweight_1_0_%d.hdf"
                        )
                        % index_of_group
                    )
                    EMAN2_cppwrap.Util.fuse_low_freq(
                        tvol0, tvol1, tweight0, tweight1, 2 * Tracker["fuse_freq"]
                    )
                    tag = 7007
                    sp_utilities.send_EMData(
                        tvol1,
                        Blockdata["no_of_processes_per_group"] - 1,
                        tag,
                        Blockdata["shared_comm"],
                    )
                    sp_utilities.send_EMData(
                        tweight1,
                        Blockdata["no_of_processes_per_group"] - 1,
                        tag,
                        Blockdata["shared_comm"],
                    )
                    shrank0 = stepone(tvol0, tweight0)
                elif (
                    Blockdata["myid_on_node"]
                    == Blockdata["no_of_processes_per_group"] - 1
                ):
                    tag = 7007
                    tvol1 = sp_utilities.recv_EMData(0, tag, Blockdata["shared_comm"])
                    tweight1 = sp_utilities.recv_EMData(
                        0, tag, Blockdata["shared_comm"]
                    )
                    tvol1.set_attr_dict(
                        {
                            "is_complex": 1,
                            "is_fftodd": 1,
                            "is_complex_ri": 1,
                            "is_fftpad": 1,
                        }
                    )
                    shrank1 = stepone(tvol1, tweight1)

                if Blockdata["myid_on_node"] == 1:
                    tvol0 = sp_utilities.get_im(
                        os.path.join(Tracker["directory"], "tempdir", "tvol_0_1_%d.hdf")
                        % index_of_group
                    )
                    tweight0 = sp_utilities.get_im(
                        os.path.join(
                            Tracker["directory"], "tempdir", "tweight_0_1_%d.hdf"
                        )
                        % index_of_group
                    )
                    tvol1 = sp_utilities.get_im(
                        os.path.join(Tracker["directory"], "tempdir", "tvol_1_1_%d.hdf")
                        % index_of_group
                    )
                    tweight1 = sp_utilities.get_im(
                        os.path.join(
                            Tracker["directory"], "tempdir", "tweight_1_1_%d.hdf"
                        )
                        % index_of_group
                    )
                    EMAN2_cppwrap.Util.fuse_low_freq(
                        tvol0, tvol1, tweight0, tweight1, 2 * Tracker["fuse_freq"]
                    )
                    tag = 7007
                    sp_utilities.send_EMData(
                        tvol1,
                        Blockdata["no_of_processes_per_group"] - 2,
                        tag,
                        Blockdata["shared_comm"],
                    )
                    sp_utilities.send_EMData(
                        tweight1,
                        Blockdata["no_of_processes_per_group"] - 2,
                        tag,
                        Blockdata["shared_comm"],
                    )
                    shrank0 = stepone(tvol0, tweight0)
                elif (
                    Blockdata["myid_on_node"]
                    == Blockdata["no_of_processes_per_group"] - 2
                ):
                    tag = 7007
                    tvol1 = sp_utilities.recv_EMData(1, tag, Blockdata["shared_comm"])
                    tweight1 = sp_utilities.recv_EMData(
                        1, tag, Blockdata["shared_comm"]
                    )
                    tvol1.set_attr_dict(
                        {
                            "is_complex": 1,
                            "is_fftodd": 1,
                            "is_complex_ri": 1,
                            "is_fftpad": 1,
                        }
                    )
                    shrank1 = stepone(tvol1, tweight1)

                mpi.mpi_barrier(Blockdata["shared_comm"])
                if Blockdata["myid_on_node"] == 0:
                    tag = 7007
                    sp_utilities.send_EMData(
                        shrank0,
                        Blockdata["no_of_processes_per_group"] - 1,
                        tag,
                        Blockdata["shared_comm"],
                    )
                    del shrank0
                    lcfsc = 0
                elif (
                    Blockdata["myid_on_node"]
                    == Blockdata["no_of_processes_per_group"] - 1
                ):
                    tag = 7007
                    shrank0 = sp_utilities.recv_EMData(0, tag, Blockdata["shared_comm"])
                    cfsc = sp_statistics.fsc(shrank0, shrank1)[1]
                    sp_utilities.write_text_row(
                        cfsc,
                        os.path.join(
                            Tracker["directory"],
                            "fsc_driver_chunk0_grp%03d_iter%03d.txt",
                        )
                        % (index_of_group, iteration),
                    )
                    del shrank0, shrank1
                    if Tracker["nxinit"] < Tracker["constants"]["nnxo"]:
                        cfsc = cfsc[: Tracker["nxinit"]]
                        for i in range(
                            len(cfsc), old_div(Tracker["constants"]["nnxo"], 2) + 1
                        ):
                            cfsc.append(0.0)
                    lcfsc = len(cfsc)
                    fsc05 = 0
                    fsc143 = 0
                    for ifreq in range(len(cfsc)):
                        if cfsc[ifreq] < 0.5:
                            break
                    fsc05 = ifreq - 1
                    for ifreq in range(len(cfsc)):
                        if cfsc[ifreq] < 0.143:
                            break
                    fsc143 = ifreq - 1
                    Tracker["fsc143"] = fsc143
                    Tracker["fsc05"] = fsc05
                if Blockdata["myid_on_node"] == 1:
                    tag = 7007
                    sp_utilities.send_EMData(
                        shrank0,
                        Blockdata["no_of_processes_per_group"] - 2,
                        tag,
                        Blockdata["shared_comm"],
                    )
                    del shrank0
                    lcfsc = 0
                elif (
                    Blockdata["myid_on_node"]
                    == Blockdata["no_of_processes_per_group"] - 2
                ):
                    tag = 7007
                    shrank0 = sp_utilities.recv_EMData(1, tag, Blockdata["shared_comm"])
                    cfsc = sp_statistics.fsc(shrank0, shrank1)[1]
                    sp_utilities.write_text_row(
                        cfsc,
                        os.path.join(
                            Tracker["directory"],
                            "fsc_driver_chunk1_grp%03d_iter%03d.txt",
                        )
                        % (index_of_group, iteration),
                    )
                    del shrank0, shrank1
                    if Tracker["nxinit"] < Tracker["constants"]["nnxo"]:
                        cfsc = cfsc[: Tracker["nxinit"]]
                        for i in range(
                            len(cfsc), old_div(Tracker["constants"]["nnxo"], 2) + 1
                        ):
                            cfsc.append(0.0)
                    lcfsc = len(cfsc)
                    fsc05 = 0
                    fsc143 = 0
                    for ifreq in range(len(cfsc)):
                        if cfsc[ifreq] < 0.5:
                            break
                    fsc05 = ifreq - 1
                    for ifreq in range(len(cfsc)):
                        if cfsc[ifreq] < 0.143:
                            break
                    fsc143 = ifreq - 1
                    Tracker["fsc143"] = fsc143
                    Tracker["fsc05"] = fsc05
                Tracker = sp_utilities.wrap_mpi_bcast(
                    Tracker,
                    Blockdata["no_of_processes_per_group"] - 1,
                    Blockdata["shared_comm"],
                )
                if Blockdata["myid_on_node"] == 0:
                    res_05[index_of_group] = Tracker["fsc05"]
                    res_143[index_of_group] = Tracker["fsc143"]
                mpi.mpi_barrier(Blockdata["shared_comm"])
            mpi.mpi_barrier(Blockdata["shared_comm"])
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    keepgoing = sp_utilities.bcast_number_to_all(
        keepgoing, source_node=Blockdata["main_node"], mpi_comm=mpi.MPI_COMM_WORLD
    )  # always check
    Tracker = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"])
    if not keepgoing:
        sp_global_def.ERROR(
            "do3d_sorting_groups_trl_iter  %s"
            % os.path.join(Tracker["directory"], "tempdir"),
            "do3d_sorting_groups_trl_iter",
            1,
            Blockdata["myid"],
        )
    return


def do3d_sorting_group_insertion_random_two_for_fsc(
    data, sparamstructure, snorm_per_particle
):
    global Tracker, Blockdata
    if Blockdata["myid"] == Blockdata["nodes"][0]:
        if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
            os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
    for index_of_groups in range(Tracker["number_of_groups"]):
        for procid in range(2):
            for ifsc in range(2):
                tvol, tweight, trol = recons3d_4nnsorting_group_fsc_MPI(
                    myid=Blockdata["myid"],
                    main_node=Blockdata["nodes"][procid],
                    prjlist=data,
                    fsc_half=ifsc,
                    random_subset=procid,
                    group_ID=index_of_groups,
                    paramstructure=sparamstructure,
                    norm_per_particle=snorm_per_particle,
                    CTF=Tracker["constants"]["CTF"],
                    upweighted=False,
                    target_size=(2 * Tracker["nxinit"] + 3),
                )
                if Blockdata["myid"] == Blockdata["nodes"][procid]:
                    tvol.set_attr("is_complex", 0)
                    tvol.write_image(
                        os.path.join(
                            Tracker["directory"],
                            "tempdir",
                            "tvol_%d_%d_%d.hdf" % (ifsc, procid, index_of_groups),
                        )
                    )
                    tweight.write_image(
                        os.path.join(
                            Tracker["directory"],
                            "tempdir",
                            "tweight_%d_%d_%d.hdf" % (ifsc, procid, index_of_groups),
                        )
                    )
                    trol.write_image(
                        os.path.join(
                            Tracker["directory"],
                            "tempdir",
                            "trol_%d_%d_%d.hdf" % (ifsc, procid, index_of_groups),
                        )
                    )
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    return


#### Never read volumes


def recons3d_4nnsorting_group_fsc_MPI(
    myid,
    main_node,
    prjlist,
    fsc_half,
    random_subset,
    group_ID,
    paramstructure,
    norm_per_particle,
    CTF=True,
    upweighted=True,
    mpi_comm=None,
    target_size=-1,
):
    ##      with smearing
    #####    recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
    ####    Input
    ####    list_of_prjlist: list of lists of projections to be included in the reconstruction
    global Tracker, Blockdata
    if mpi_comm == None:
        mpi_comm = mpi.MPI_COMM_WORLD
    imgsize = prjlist[0].get_ysize()  # It can be Fourier, so take y-size
    refvol = sp_utilities.model_blank(target_size)
    refvol.set_attr("fudge", 1.0)
    if CTF:
        do_ctf = 1
    else:
        do_ctf = 0
    fftvol = EMAN2_cppwrap.EMData()
    weight = EMAN2_cppwrap.EMData()
    try:
        qt = prjlist[0].get_attr("qt")
    except:
        qt = 1.0
    params = {
        "size": target_size,
        "npad": 2,
        "snr": 1.0,
        "sign": 1,
        "symmetry": "c1",
        "refvol": refvol,
        "fftvol": fftvol,
        "weight": weight,
        "do_ctf": do_ctf,
    }
    r = EMAN2_cppwrap.Reconstructors.get("nn4_ctfw", params)
    r.setup()
    if norm_per_particle == None:
        norm_per_particle = len(prjlist) * [1.0]
    # Definitions for smearing ,all copied from refinement
    if not Tracker["nosmearing"]:
        delta = Tracker["delta"]
        refang = Tracker["refang"]
        rshifts_shrank = copy.deepcopy(Tracker["rshifts"])
        nshifts = len(rshifts_shrank)
        for im in range(nshifts):
            rshifts_shrank[im][0] *= old_div(
                float(Tracker["nxinit"]), float(Tracker["constants"]["nnxo"])
            )
            rshifts_shrank[im][1] *= old_div(
                float(Tracker["nxinit"]), float(Tracker["constants"]["nnxo"])
            )
    nnx = prjlist[0].get_xsize()
    nny = prjlist[0].get_ysize()
    nc = 0
    for im in range(len(prjlist)):
        if (
            prjlist[im].get_attr("group") == group_ID
            and prjlist[im].get_attr("chunk_id") == random_subset
        ):
            if Tracker["nosmearing"]:
                avgnorm = 1.0
            else:
                avgnorm = Tracker["avgnorm"][prjlist[im].get_attr("chunk_id")]  #
            if nc % 2 == fsc_half:
                if Tracker["nosmearing"]:
                    ct = prjlist[im].get_attr("ctf")
                    bckgn = prjlist[im].get_attr("bckgnoise")
                    if not upweighted:
                        prjlist[im] = sp_filter.filt_table(prjlist[im], bckgn)
                    prjlist[im].set_attr_dict({"bckgnoise": bckgn, "ctf": ct})
                    phi, theta, psi, s2x, s2y = sp_utilities.get_params_proj(
                        prjlist[im], xform="xform.projection"
                    )
                    r.insert_slice(
                        prjlist[im],
                        EMAN2_cppwrap.Transform(
                            {"type": "spider", "phi": phi, "theta": theta, "psi": psi}
                        ),
                        1.0,
                    )
                else:
                    if Tracker["constants"]["nsmear"] <= 0.0:
                        numbor = len(paramstructure[im][2])
                    else:
                        numbor = 1
                    ipsiandiang = [
                        old_div(paramstructure[im][2][i][0], 1000)
                        for i in range(numbor)
                    ]
                    allshifts = [
                        paramstructure[im][2][i][0] % 1000 for i in range(numbor)
                    ]
                    probs = [paramstructure[im][2][i][1] for i in range(numbor)]
                    #  Find unique projection directions
                    tdir = list(set(ipsiandiang))
                    bckgn = prjlist[im].get_attr("bckgnoise")
                    ct = prjlist[im].get_attr("ctf")
                    #  For each unique projection direction:
                    data = [None] * nshifts
                    for ii in range(len(tdir)):
                        #  Find the number of times given projection direction appears on the list, it is the number of different shifts associated with it.
                        lshifts = sp_utilities.findall(tdir[ii], ipsiandiang)
                        toprab = 0.0
                        for ki in range(len(lshifts)):
                            toprab += probs[lshifts[ki]]
                        recdata = EMAN2_cppwrap.EMData(nny, nny, 1, False)
                        recdata.set_attr("is_complex", 0)
                        for ki in range(len(lshifts)):
                            lpt = allshifts[lshifts[ki]]
                            if data[lpt] == None:
                                data[lpt] = sp_fundamentals.fshift(
                                    prjlist[im],
                                    rshifts_shrank[lpt][0],
                                    rshifts_shrank[lpt][1],
                                )
                                data[lpt].set_attr("is_complex", 0)
                            EMAN2_cppwrap.Util.add_img(
                                recdata,
                                EMAN2_cppwrap.Util.mult_scalar(
                                    data[lpt], old_div(probs[lshifts[ki]], toprab)
                                ),
                            )
                        recdata.set_attr_dict(
                            {
                                "padffted": 1,
                                "is_fftpad": 1,
                                "is_fftodd": 0,
                                "is_complex_ri": 1,
                                "is_complex": 1,
                            }
                        )
                        if not upweighted:
                            recdata = sp_filter.filt_table(recdata, bckgn)
                        recdata.set_attr_dict({"bckgnoise": bckgn, "ctf": ct})
                        ipsi = tdir[ii] % 100000
                        iang = old_div(tdir[ii], 100000)
                        r.insert_slice(
                            recdata,
                            EMAN2_cppwrap.Transform(
                                {
                                    "type": "spider",
                                    "phi": refang[iang][0],
                                    "theta": refang[iang][1],
                                    "psi": refang[iang][2] + ipsi * delta,
                                }
                            ),
                            old_div(toprab * avgnorm, norm_per_particle[im]),
                        )
            nc += 1
    sp_utilities.reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
    sp_utilities.reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)
    if myid == main_node:
        dummy = r.finish(True)
    mpi.mpi_barrier(mpi_comm)
    if myid == main_node:
        return fftvol, weight, refvol
    else:
        return None, None, None


#####end of FSC
###=====<-----group rec3d
### insertion


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
    """Return memory usage in bytes.
    """
    return _VmB("VmSize:") - since


def compute_final_map(work_dir):
    global Tracker, Blockdata
    Tracker["constants"]["orgres"] = 0.0
    Tracker["constants"]["refinement_delta"] = 0.0
    Tracker["constants"]["refinement_ts"] = 0.0
    Tracker["constants"]["refinement_xr"] = 0.0
    Tracker["constants"]["refinement_an"] = 0.0
    minimum_size = Tracker["constants"]["img_per_grp"]
    number_of_groups = 0

    if Blockdata["myid"] == Blockdata["main_node"]:
        final_accounted_ptl = 0
        fout = open(os.path.join(work_dir, "Tracker.json"), "w")
        json.dump(Tracker, fout)
        fout.close()
        clusters = []
        while os.path.exists(
            os.path.join(work_dir, "Cluster_%03d.txt" % number_of_groups)
        ):
            class_in = sp_utilities.read_text_file(
                os.path.join(work_dir, "Cluster_%03d.txt" % number_of_groups)
            )
            minimum_size = min(len(class_in), minimum_size)
            number_of_groups += 1
            final_accounted_ptl += len(class_in)
            clusters.append(class_in)
            del class_in
        Tracker["total_stack"] = final_accounted_ptl
        Tracker["number_of_groups"] = number_of_groups
        Tracker["nxinit"] = Tracker["nxinit_refinement"]
    else:
        Tracker = 0
    number_of_groups = sp_utilities.bcast_number_to_all(
        number_of_groups, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    Tracker = sp_utilities.wrap_mpi_bcast(
        Tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    if number_of_groups == 0:
        sp_global_def.ERROR(
            "No clusters  found, the program terminates.",
            "compute_final_map",
            1,
            Blockdata["myid"],
        )
    compute_noise(Tracker["nxinit"])
    if Blockdata["myid"] == Blockdata["main_node"]:
        alist, partition = merge_classes_into_partition_list(clusters)
        sp_utilities.write_text_row(
            partition, os.path.join(work_dir, "generation_partition.txt")
        )
    parti_file = os.path.join(work_dir, "generation_partition.txt")
    params = os.path.join(
        Tracker["constants"]["masterdir"], "refinement_parameters.txt"
    )
    previous_params = Tracker["previous_parstack"]
    original_data, norm_per_particle = read_data_for_sorting(
        parti_file, params, previous_params
    )

    if Tracker["nosmearing"]:
        parameterstructure = None
        paramstructure_dict = None
        paramstructure_dir = None
    else:
        paramstructure_dict = Tracker["paramstructure_dict"]
        paramstructure_dir = Tracker["paramstructure_dir"]
        parameterstructure = read_paramstructure_for_sorting(
            parti_file, paramstructure_dict, paramstructure_dir
        )
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    Tracker["directory"] = work_dir
    compute_noise(Tracker["nxinit"])
    rdata = downsize_data_for_rec3D(
        original_data, Tracker["nxinit"], False, 1
    )  # pay attentions to shifts!
    for im in range(len(original_data)):
        rdata[im].set_attr("group", original_data[im].get_attr("group"))
    del original_data
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    # srdata = precalculate_shifted_data_for_recons3D(rdata, parameterstructure, Tracker["refang"], Tracker["rshifts"], \
    #  Tracker["delta"], Tracker["avgnorm"], Tracker["nxinit"], Tracker["constants"]["nnxo"], Tracker["nosmearing"], \
    #      norm_per_particle,  Tracker["constants"]["nsmear"])
    # del rdata, parameterstructure, norm_per_particle
    # mpi_barrier(MPI_COMM_WORLD)
    do3d_sorting_groups_nofsc_final(rdata, parameterstructure, norm_per_particle)
    del rdata
    if Blockdata["myid"] == Blockdata["main_node"]:
        if os.path.exists(os.path.join(Tracker["constants"]["masterdir"], "tempdir")):
            shutil.rmtree(os.path.join(Tracker["constants"]["masterdir"], "tempdir"))
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    return


### final rec3d


def do3d_sorting_groups_nofsc_final(rdata, parameterstructure, norm_per_particle):
    global Tracker, Blockdata
    keepgoing = 1
    if Blockdata["myid"] == Blockdata["last_node"]:
        if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
            os.mkdir(os.path.join(Tracker["directory"], "tempdir"))
        try:
            fout = open(os.path.join(Tracker["directory"], "freq_cutoff.json"), "r")
            freq_cutoff_dict = sp_utilities.convert_json_fromunicode(json.load(fout))
            fout.close()
        except:
            freq_cutoff_dict = 0
    else:
        freq_cutoff_dict = 0
    freq_cutoff_dict = sp_utilities.wrap_mpi_bcast(
        freq_cutoff_dict, Blockdata["last_node"], mpi.MPI_COMM_WORLD
    )

    for index_of_groups in range(Tracker["number_of_groups"]):
        tvol, tweight, trol = recons3d_trl_struct_group_MPI(
            Blockdata["myid"],
            Blockdata["last_node"],
            rdata,
            2,
            index_of_groups,
            parameterstructure,
            norm_per_particle,
            True,
            None,
            Tracker["constants"]["CTF"],
            (2 * Tracker["nxinit"] + 3),
            Tracker["nosmearing"],
        )
        if Blockdata["myid"] == Blockdata["last_node"]:
            tvol.set_attr("is_complex", 0)
            tvol.write_image(
                os.path.join(
                    Tracker["directory"], "tempdir", "tvol_2_%d.hdf" % index_of_groups
                )
            )
            tweight.write_image(
                os.path.join(
                    Tracker["directory"],
                    "tempdir",
                    "tweight_2_%d.hdf" % index_of_groups,
                )
            )
            trol.write_image(
                os.path.join(
                    Tracker["directory"], "tempdir", "trol_2_%d.hdf" % index_of_groups
                )
            )
            del tvol
            del tweight
            del trol
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    Tracker["fsc143"] = 0
    Tracker["fsc05"] = 0
    Tracker["maxfrad"] = old_div(Tracker["nxinit"], 2)
    if Blockdata["no_of_groups"] > 1:
        # new starts
        sub_main_node_list = [-1 for i in range(Blockdata["no_of_groups"])]
        for index_of_colors in range(Blockdata["no_of_groups"]):
            for iproc in range(Blockdata["nproc"] - 1):
                if Blockdata["myid"] == iproc:
                    if (
                        Blockdata["color"] == index_of_colors
                        and Blockdata["myid_on_node"] == 0
                    ):
                        sub_main_node_list[index_of_colors] = Blockdata["myid"]
                    sp_utilities.wrap_mpi_send(
                        sub_main_node_list, Blockdata["last_node"], mpi.MPI_COMM_WORLD
                    )
                if Blockdata["myid"] == Blockdata["last_node"]:
                    dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
                    for im in range(len(dummy)):
                        if dummy[im] > -1:
                            sub_main_node_list[im] = dummy[im]
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        sp_utilities.wrap_mpi_bcast(
            sub_main_node_list, Blockdata["last_node"], mpi.MPI_COMM_WORLD
        )
        if Tracker["number_of_groups"] % Blockdata["no_of_groups"] == 0:
            nbig_loop = old_div(Tracker["number_of_groups"], Blockdata["no_of_groups"])
        else:
            nbig_loop = (
                old_div(Tracker["number_of_groups"], Blockdata["no_of_groups"]) + 1
            )

        big_loop_colors = [[] for i in range(nbig_loop)]
        big_loop_groups = [[] for i in range(nbig_loop)]
        nc = 0
        while nc < Tracker["number_of_groups"]:
            im = old_div(nc, Blockdata["no_of_groups"])
            jm = nc % Blockdata["no_of_groups"]
            big_loop_colors[im].append(jm)
            big_loop_groups[im].append(nc)
            nc += 1
        for iloop in range(nbig_loop):
            for im in range(len(big_loop_colors[iloop])):
                index_of_group = big_loop_groups[iloop][im]
                index_of_colors = big_loop_colors[iloop][im]

                if Blockdata["myid"] == Blockdata["last_node"]:
                    tvol2 = sp_utilities.get_im(
                        os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")
                        % index_of_group
                    )
                    tweight2 = sp_utilities.get_im(
                        os.path.join(
                            Tracker["directory"], "tempdir", "tweight_2_%d.hdf"
                        )
                        % index_of_group
                    )
                    treg2 = sp_utilities.get_im(
                        os.path.join(
                            Tracker["directory"],
                            "tempdir",
                            "trol_2_%d.hdf" % index_of_group,
                        )
                    )
                    tvol2.set_attr_dict(
                        {
                            "is_complex": 1,
                            "is_fftodd": 1,
                            "is_complex_ri": 1,
                            "is_fftpad": 1,
                        }
                    )
                    tag = 7007
                    sp_utilities.send_EMData(
                        tvol2,
                        sub_main_node_list[index_of_colors],
                        tag,
                        mpi.MPI_COMM_WORLD,
                    )
                    sp_utilities.send_EMData(
                        tweight2,
                        sub_main_node_list[index_of_colors],
                        tag,
                        mpi.MPI_COMM_WORLD,
                    )
                    sp_utilities.send_EMData(
                        treg2,
                        sub_main_node_list[index_of_colors],
                        tag,
                        mpi.MPI_COMM_WORLD,
                    )

                elif Blockdata["myid"] == sub_main_node_list[index_of_colors]:
                    tag = 7007
                    tvol2 = sp_utilities.recv_EMData(
                        Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                    tweight2 = sp_utilities.recv_EMData(
                        Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                    treg2 = sp_utilities.recv_EMData(
                        Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            for im in range(len(big_loop_colors[iloop])):
                index_of_group = big_loop_groups[iloop][im]
                index_of_colors = big_loop_colors[iloop][im]
                try:
                    Tracker["freq_fsc143_cutoff"] = freq_cutoff_dict[
                        "Cluster_%03d.txt" % index_of_group
                    ]
                except:
                    pass
                if (
                    Blockdata["color"] == index_of_colors
                ):  # It has to be 1 to avoid problem with tvol1 not closed on the disk
                    if Blockdata["myid_on_node"] != 0:
                        tvol2 = sp_utilities.model_blank(1)
                        tweight2 = sp_utilities.model_blank(1)
                        treg2 = sp_utilities.model_blank(1)
                    tvol2 = steptwo_mpi_filter(
                        tvol2,
                        tweight2,
                        treg2,
                        None,
                        Tracker["freq_fsc143_cutoff"],
                        0.01,
                        False,
                        color=index_of_colors,
                    )  # has to be False!!!
                    del tweight2, treg2
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            for im in range(len(big_loop_colors[iloop])):
                index_of_group = big_loop_groups[iloop][im]
                index_of_colors = big_loop_colors[iloop][im]
                if (Blockdata["color"] == index_of_colors) and (
                    Blockdata["myid_on_node"] == 0
                ):
                    tag = 7007
                    sp_utilities.send_EMData(
                        tvol2, Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD
                    )
                elif Blockdata["myid"] == Blockdata["last_node"]:
                    tag = 7007
                    tvol2 = sp_utilities.recv_EMData(
                        sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD
                    )
                    tvol2.write_image(
                        os.path.join(
                            Tracker["directory"], "vol_cluster%03d.hdf" % index_of_group
                        )
                    )
                    del tvol2
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    else:  # loop over all groups for single node
        for index_of_group in range(Tracker["number_of_groups"]):
            if Blockdata["myid_on_node"] == 0:
                tvol2 = sp_utilities.get_im(
                    os.path.join(Tracker["directory"], "tempdir", "tvol_2_%d.hdf")
                    % index_of_group
                )
                tweight2 = sp_utilities.get_im(
                    os.path.join(Tracker["directory"], "tempdir", "tweight_2_%d.hdf")
                    % index_of_group
                )
                treg2 = sp_utilities.get_im(
                    os.path.join(
                        Tracker["directory"],
                        "tempdir",
                        "trol_2_%d.hdf" % index_of_group,
                    )
                )
                tvol2.set_attr_dict(
                    {
                        "is_complex": 1,
                        "is_fftodd": 1,
                        "is_complex_ri": 1,
                        "is_fftpad": 1,
                    }
                )
            else:
                tvol2 = sp_utilities.model_blank(1)
                tweight2 = sp_utilities.model_blank(1)
                treg2 = sp_utilities.model_blank(1)
            try:
                Tracker["freq_fsc143_cutoff"] = freq_cutoff_dict[
                    "Cluster_%03d.txt" % index_of_group
                ]
            except:
                pass
            tvol2 = steptwo_mpi_filter(
                tvol2, tweight2, treg2, None, Tracker["freq_fsc143_cutoff"], 0.01, False
            )  # has to be False!!!
            del tweight2, treg2
            if Blockdata["myid_on_node"] == 0:
                tvol2.write_image(
                    os.path.join(
                        Tracker["directory"], "vol_cluster%03d.hdf" % index_of_group
                    )
                )
                del tvol2
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    keepgoing = sp_utilities.bcast_number_to_all(
        keepgoing, source_node=Blockdata["main_node"], mpi_comm=mpi.MPI_COMM_WORLD
    )  # always check
    Tracker = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"])
    if not keepgoing:
        sp_global_def.ERROR(
            "do3d_sorting_groups_nofsc_final  %s"
            % os.path.join(Tracker["directory"], "tempdir"),
            "do3d_sorting_groups_nofsc_final",
            1,
            Blockdata["myid"],
        )
    return


#####==========----various utilities


def get_time(time_start):
    current_time = time.time() - time_start
    current_time_h = old_div(current_time, 3600)
    current_time_m = old_div((current_time - current_time_h * 3600), 60)
    return int(current_time_h), int(current_time_m)


def check_sorting(total_data, keepsorting, log_file):
    global Tracker, Blockdata
    if Blockdata["myid"] == Blockdata["main_node"]:
        fout = open(
            os.path.join(Tracker["constants"]["masterdir"], "Tracker.json"), "r"
        )
        Tracker_main = sp_utilities.convert_json_fromunicode(json.load(fout))
        fout.close()
    else:
        Tracker_main = 0
    Tracker_main = sp_utilities.wrap_mpi_bcast(Tracker_main, Blockdata["main_node"])
    n = 0
    keepsorting = 0
    while n <= Tracker_main["constants"]["img_per_grp_split_rate"]:
        if (
            old_div(Tracker_main["constants"]["img_per_grp"], 2) ** n
            > Tracker_main["constants"]["minimum_grp_size"]
        ):
            Tracker["number_of_groups"] = old_div(
                (2 ** n * total_data), Tracker_main["constants"]["img_per_grp"]
            )
            if Tracker["number_of_groups"] >= 2:
                Tracker["current_img_per_grp"] = (
                    old_div(Tracker_main["constants"]["img_per_grp"], 2) ** n
                )
                keepsorting = 1
        else:
            Tracker["number_of_groups"] = old_div(
                total_data, Tracker_main["constants"]["minimum_grp_size"]
            )
            if Tracker["number_of_groups"] >= 2:
                Tracker["current_img_per_grp"] = Tracker_main["constants"][
                    "minimum_grp_size"
                ]
                keepsorting = 1
            else:
                keepsorting = 0

        if Tracker["number_of_groups"] >= 2:
            keepsorting = 1
            break
        else:
            n += 1
    if keepsorting == 1:
        Tracker["total_stack"] = total_data
        keepsorting = sort3d_init("initialization", log_file)
    return keepsorting


def copy_results(log_file, all_gen_stat_list):
    global Tracker, Blockdata
    nsgen = 0
    for i in range(len(all_gen_stat_list)):
        nsgen += len(all_gen_stat_list[i])
    if nsgen > 0:
        if Blockdata["myid"] == Blockdata["main_node"]:
            log_file.add(
                "================================================================================================================"
            )
            log_file.add(
                "                     Final results saved in %s"
                % Tracker["constants"]["masterdir"]
            )
            log_file.add(
                "----------------------------------------------------------------------------------------------------------------"
            )
            nclusters = 0
            log_file.add(
                "{:^8} {:>8}   {:^24}  {:>15} {:^22} {:^5} {:^15} {:^20} ".format(
                    "Group ID",
                    "    size",
                    "determined in generation",
                    "reproducibility",
                    "random reproducibility",
                    " std ",
                    " selection file",
                    "       map file     ",
                )
            )
            clusters = []
            NACC = 0
            for ig1, value in list(Tracker["generation"].items()):
                ig = int(ig1)
                for ic in range(value):
                    cluster_file = os.path.join(
                        Tracker["constants"]["masterdir"],
                        "generation_%03d" % ig,
                        "Cluster_%03d.txt" % ic,
                    )
                    shutil.copyfile(
                        cluster_file,
                        os.path.join(
                            Tracker["constants"]["masterdir"],
                            "Cluster_%03d.txt" % nclusters,
                        ),
                    )
                    clusters.append(sp_utilities.read_text_file(cluster_file))
                    # copyfile(os.path.join(Tracker["constants"]["masterdir"], "generation_%03d"%ig, "vol_grp%03d_iter000.hdf"%ic), os.path.join(Tracker["constants"]["masterdir"], "vol_cluster%03d.hdf"%nclusters))
                    cluster = sp_utilities.read_text_file(
                        os.path.join(
                            Tracker["constants"]["masterdir"],
                            "generation_%03d" % ig,
                            "Cluster_%03d.txt" % ic,
                        )
                    )
                    cluster_file = "Cluster_%03d.txt" % nclusters
                    vol_file = "vol_cluster%03d.hdf" % nclusters
                    msg = "{:>8} {:>8}   {:^24}        {:^6}          {:^6}          {:>5}  {:^15} {:^20} ".format(
                        nclusters,
                        len(cluster),
                        ig,
                        round(all_gen_stat_list[ig][ic][0], 1),
                        round(all_gen_stat_list[ig][ic][1], 1),
                        round(all_gen_stat_list[ig][ic][2], 1),
                        cluster_file,
                        vol_file,
                    )
                    nclusters += 1
                    NACC += len(cluster)
                    log_file.add(msg)
            Unaccounted_file = os.path.join(
                Tracker["constants"]["masterdir"],
                "generation_%03d" % ig,
                "Unaccounted.txt",
            )
            shutil.copyfile(
                Unaccounted_file,
                os.path.join(Tracker["constants"]["masterdir"], "Unaccounted.txt"),
            )
            NUACC = Tracker["constants"]["total_stack"] - NACC
            log_file.add(
                "{:^7} {:^8} {:^22} {:^8} {:^24} {:^8} ".format(
                    " Images",
                    Tracker["constants"]["total_stack"],
                    "accounted for images: ",
                    NACC,
                    "unaccounted for images: ",
                    NUACC,
                )
            )
            log_file.add("Unaccounted images saved in Unaccounted.txt")
            if len(clusters) >= 2:
                do_analysis_on_identified_clusters(clusters, log_file)
            else:
                log_file.add(" ANOVA analysis is skipped ")
            fout = open(
                os.path.join(Tracker["constants"]["masterdir"], "Tracker.json"), "w"
            )
            json.dump(Tracker, fout)
            fout.close()
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    else:
        if Blockdata["myid"] == Blockdata["main_node"]:
            log_file.add("No groups are found.\n")
    return


def get_MGR_from_two_way_comparison(newindeces, clusters1, clusters2, N):
    rnd_grp_sizes = {}
    K = len(newindeces)
    reordered_cluster2 = [None for i in range(len(clusters2))]
    for ij in range(len(newindeces)):
        reordered_cluster2[newindeces[ij][0]] = clusters2[newindeces[ij][1]]
    table_k_k = [[None for i in range(K)] for j in range(K)]
    for i in range(K):
        for j in range(K):
            if not (clusters1[i] is None) and not (reordered_cluster2[j] is None):
                table_k_k[i][j] = len(
                    set(clusters1[i]).intersection(set(reordered_cluster2[j]))
                )
            else:
                table_k_k[i][j] = 0
    sum_rows = [0 for i in range(K)]
    sum_cols = [0 for i in range(K)]
    for i in range(K):
        for j in range(K):
            sum_rows[i] += table_k_k[i][j]
            sum_cols[j] += table_k_k[i][j]
    diagonal_k = [None for i in range(K)]
    for i in range(K):
        diagonal_k[i] = old_div((sum_rows[i] + sum_cols[i]), 2)  # nk
    min_sizes = [None for i in range(K)]
    for i in range(K):
        min_sizes[i] = old_div(diagonal_k[i] ** 2, N)
    return min_sizes


def estimate_tanhl_params(cutoff, taa, image_size):
    def tanhfl(x, cutoff, taa):
        omega = cutoff
        cnst = old_div(numpy.pi, (2.0 * omega * taa))
        v1 = cnst * (x + omega)
        v2 = cnst * (x - omega)
        return 0.5 * (numpy.tanh(v1) - numpy.tanh(v2))

    def get_filter(cutoff1, taa, image_size):
        values = []
        N = old_div(image_size, 2)
        for im in range(N):
            x = old_div(float(im), float(image_size))
            values.append(tanhfl(x, cutoff1, taa))
        return values

    values = get_filter(cutoff, taa, image_size)
    icutoff = image_size
    init = int(cutoff * image_size)
    while icutoff >= cutoff * image_size:
        cutoff1 = old_div(float(init), image_size)
        values = get_filter(cutoff1, taa, image_size)
        for im in range(len(values)):
            if values[im] <= 0.0:
                icutoff = im
                break
        init -= 1
    return cutoff1, taa


def print_matching_pairs(pair_list, log_file):
    log_file.add(" ")
    log_file.add("                        Two-way matching of sorting results.")
    log_file.add(
        "M indicates that respective group of P0 sorting (row number) matches respective group of P1 sorting (column number)"
    )

    msg = "   "
    for i in range(len(pair_list)):
        msg += "{:^5d}".format(i)
    log_file.add(msg)
    for im in range(len(pair_list)):
        msg = "{:^3d}".format(im)
        for jm in range(len(pair_list)):
            not_found = True
            for km in range(len(pair_list)):
                if pair_list[km][0] == im and pair_list[km][1] == jm:
                    msg += "{:^5s}".format("M")
                    not_found = False
            if not_found:
                msg += "{:^5s}".format(" ")
        log_file.add(msg)
    return


def do_random_groups_simulation_mpi(ptp1, ptp2):
    global Tracker, Blockdata
    # return two lists: group avgs and group stds. The last one of two lists are the total avg and std.
    if (len(ptp1) >= 50) or (len(ptp2) >= 50):
        if Blockdata["myid"] == Blockdata["main_node"]:
            print("Warning: too many simulaton groups")
    Nloop = max(old_div(1000, Blockdata["nproc"]), 1)
    NT = 1000
    a = []
    b = []
    for ic in range(len(ptp1)):
        a += ptp1[ic]
    for ic in range(len(ptp2)):
        b += ptp2[ic]
    tsize = float(len(set(a + b)))

    nsize1 = 0
    plist1 = []
    for i1 in range(len(ptp1)):
        plist1.append(
            [nsize1, nsize1 + max(int(old_div(float(len(ptp1[i1])), tsize) * 100.0), 1)]
        )
        nsize1 += max(int(old_div(float(len(ptp1[i1])), tsize) * 100.0), 1)

    nsize2 = 0
    plist2 = []

    for i1 in range(len(ptp2)):
        plist2.append(
            [nsize2, nsize2 + max(int(old_div(float(len(ptp2[i1])), tsize) * 100.0), 1)]
        )
        nsize2 += max(int(old_div(float(len(ptp2[i1])), tsize) * 100.0), 1)

    if len(ptp1) > len(ptp2):
        for j in range(len(ptp1) - len(ptp2)):
            plist2.append([nsize2, nsize2 + 1])
            nsize2 += 1

    elif len(ptp2) > len(ptp1):
        for j in range(len(ptp2) - len(ptp1)):
            plist1.append([nsize1, nsize1 + 1])
            nsize1 += 1
    if (len(ptp1) >= 50) or (
        len(ptp2) >= 50
    ):  # will under-estimate random reproducibility
        alist = list(range(3 * max(len(ptp1), len(ptp2))))
        blist = list(range(3 * max(len(ptp1), len(ptp2))))
    else:
        alist = list(range(100))
        blist = list(range(100))
    k = len(plist1)
    gave = [0.0 for i in range(k)]
    gvar = [0.0 for i in range(k)]
    svar = 0.0
    save = 0.0
    alist = numpy.array(alist, "int32")
    blist = numpy.array(blist, "int32")
    for iloop in range(Nloop):
        tlist = []
        clist = [[] for i in range(k)]
        for i in range(NT):
            new_clusters1 = []
            new_clusters2 = []
            numpy.random.shuffle(alist)
            numpy.random.shuffle(blist)
            for j in range(k):
                new_clusters1.append(alist[plist1[j][0] : plist1[j][1]])
            for j in range(k):
                new_clusters2.append(blist[plist2[j][0] : plist2[j][1]])
            for j in range(k):
                new_clusters1[j] = numpy.sort(new_clusters1[j])
            for j in range(k):
                new_clusters2[j] = numpy.sort(new_clusters2[j])
            newindeces, list_stable, nb_tot_objs = sp_statistics.k_means_match_clusters_asg_new(
                new_clusters1, new_clusters2
            )
            ts = []
            for ii in range(len(newindeces)):
                ts += list(
                    set(
                        new_clusters1[newindeces[ii][0]].tolist()
                        + new_clusters2[newindeces[ii][1]].tolist()
                    )
                )
            tlist.append(old_div(nb_tot_objs, float(len(set(ts)))) * 100.0)
            for j in range(k):
                try:
                    clist[j].append(
                        old_div(
                            float(
                                (
                                    numpy.intersect1d(
                                        new_clusters1[newindeces[j][0]],
                                        new_clusters2[newindeces[j][1]],
                                    )
                                ).size
                            ),
                            float(
                                (
                                    numpy.union1d(
                                        new_clusters1[newindeces[j][0]],
                                        new_clusters2[newindeces[j][1]],
                                    )
                                ).size
                            ),
                        )
                        * 100.0
                    )
                except:
                    pass
        svar += sp_statistics.table_stat(tlist)[0]
        save += sp_statistics.table_stat(tlist)[1]
        for j in range(k):
            try:
                gvar[j] += sp_statistics.table_stat(clist[j])[1]
            except:
                gvar[j] += 0.0
            try:
                gave[j] += sp_statistics.table_stat(clist[j])[0]
            except:
                gave[j] += 0.0

    save = mpi.mpi_reduce(save, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, 0, mpi.MPI_COMM_WORLD)
    svar = mpi.mpi_reduce(svar, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, 0, mpi.MPI_COMM_WORLD)

    for iproc in range(1, Blockdata["nproc"]):
        if Blockdata["myid"] == iproc:
            sp_utilities.wrap_mpi_send(gvar, Blockdata["main_node"], mpi.MPI_COMM_WORLD)

        elif Blockdata["myid"] == Blockdata["main_node"]:
            dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
            for im in range(len(dummy)):
                gvar[im] += dummy[im]

        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    for iproc in range(1, Blockdata["nproc"]):
        if Blockdata["myid"] == iproc:
            sp_utilities.wrap_mpi_send(gave, Blockdata["main_node"], mpi.MPI_COMM_WORLD)

        elif Blockdata["myid"] == Blockdata["main_node"]:
            dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
            for im in range(len(dummy)):
                gave[im] += dummy[im]
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    if Blockdata["myid"] == Blockdata["main_node"]:
        svar = numpy.sqrt(old_div(svar, float(Blockdata["nproc"] * Nloop)))
        save = old_div(save, float(Blockdata["nproc"] * Nloop))
        for i in range(k):
            gave[i] = old_div(gave[i], float(Blockdata["nproc"] * Nloop))
            gvar[i] = numpy.sqrt(old_div(gvar[i], float(Blockdata["nproc"] * Nloop)))
        gave.append(save)
        gvar.append(svar)
    else:
        gave = 0
        gvar = 0
    gave = sp_utilities.wrap_mpi_bcast(gave, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
    gvar = sp_utilities.wrap_mpi_bcast(gvar, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
    return gave, gvar


def sorting_main_mpi(log_main, depth_order, not_include_unaccounted):
    global Tracker, Blockdata

    time_sorting_start = time.time()
    read_tracker_mpi(Tracker["constants"]["masterdir"])
    Tracker["generation"] = {}
    Tracker["current_generation"] = 0
    keepsorting = 1
    keepchecking = 1
    Tracker["current_generation"] = -1
    igen = -1
    my_pids = os.path.join(Tracker["constants"]["masterdir"], "indexes.txt")
    params = os.path.join(
        Tracker["constants"]["masterdir"], "refinement_parameters.txt"
    )
    previous_params = Tracker["previous_parstack"]
    all_gen_stat_list = []
    bad_clustering = 0
    Tracker["current_img_per_grp"] = Tracker["constants"]["img_per_grp"]
    while (keepsorting == 1) and (bad_clustering == 0):
        Tracker["current_generation"] += 1
        igen += 1
        work_dir = os.path.join(
            Tracker["constants"]["masterdir"], "generation_%03d" % igen
        )
        if Blockdata["myid"] == Blockdata["main_node"]:
            keepchecking = check_sorting_state(work_dir, keepchecking, log_main)
        else:
            keepchecking = 0
        keepchecking = sp_utilities.bcast_number_to_all(
            keepchecking, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        if keepchecking == 0:  # new, do it
            if Blockdata["myid"] == Blockdata["main_node"]:
                time_generation_start = time.time()
                if not os.path.exists(work_dir):
                    os.mkdir(work_dir)  # need check each box
                    within_generation_restart = 0
                else:
                    within_generation_restart = 1
                freq_cutoff_dict = {}
                fout = open(os.path.join(work_dir, "freq_cutoff.json"), "w")
                json.dump(freq_cutoff_dict, fout)
                fout.close()
                log_main.add(
                    "================================================================================================================"
                )
                log_main.add(
                    "                                    SORT3D IN-DEPTH   generation %d"
                    % igen
                )
                log_main.add(
                    "----------------------------------------------------------------------------------------------------------------"
                )
            else:
                within_generation_restart = 0
            within_generation_restart = sp_utilities.bcast_number_to_all(
                within_generation_restart, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            if within_generation_restart == 1:
                read_tracker_mpi(work_dir)
            else:
                dump_tracker(work_dir)
            output_list, bad_clustering, stat_list = depth_clustering(
                work_dir, depth_order, my_pids, params, previous_params, log_main
            )
            all_gen_stat_list.append(stat_list)

            if bad_clustering != 1:
                if Blockdata["myid"] == Blockdata["main_node"]:
                    clusters, nclusters, nuacc = output_clusters(
                        work_dir,
                        output_list[0][0],
                        output_list[0][1],
                        not_include_unaccounted,
                        log_main,
                    )
                    try:
                        del Tracker["generation"][str(igen)]
                    except:
                        pass
                    Tracker["generation"][igen] = len(clusters)
                else:
                    Tracker = 0
                    nclusters = 0
                    nuacc = 0
                Tracker = sp_utilities.wrap_mpi_bcast(
                    Tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD
                )
                nclusters = sp_utilities.bcast_number_to_all(
                    nclusters, Blockdata["main_node"], mpi.MPI_COMM_WORLD
                )
                nuacc = sp_utilities.bcast_number_to_all(
                    nuacc, Blockdata["main_node"], mpi.MPI_COMM_WORLD
                )
                dump_tracker(work_dir)
                if nclusters == 0:
                    if Blockdata["myid"] == Blockdata["main_node"]:
                        log_main.add("No cluster is found in generation %d" % igen)
                    break
                else:
                    if Blockdata["myid"] == Blockdata["main_node"]:
                        time_of_generation_h, time_of_generation_m = get_time(
                            time_generation_start
                        )
                        log_main.add(
                            "SORT3D generation%d time: %d hours %d minutes."
                            % (igen, time_of_generation_h, time_of_generation_m)
                        )
                    my_pids = os.path.join(work_dir, "indexes_next_generation.txt")
                    if Blockdata["myid"] == Blockdata["main_node"]:
                        sp_utilities.write_text_file(output_list[0][1], my_pids)
                        if stat_list:
                            sp_utilities.write_text_row(
                                stat_list, os.path.join(work_dir, "gen_rep.txt")
                            )
                        else:
                            with open(os.path.join(work_dir, "gen_rep.txt"), "w"):
                                pass
                        mark_sorting_state(work_dir, True, log_main)
                    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
                keepsorting = check_sorting(nuacc, keepsorting, log_main)
                if keepsorting != 1:
                    break  #
        else:  # restart run
            read_tracker_mpi(work_dir)
            my_pids = os.path.join(work_dir, "indexes_next_generation.txt")
            if Blockdata["myid"] == Blockdata["main_node"]:
                stat_list = sp_utilities.read_text_row(
                    os.path.join(work_dir, "gen_rep.txt")
                )
            else:
                stat_list = 0
            stat_list = sp_utilities.wrap_mpi_bcast(
                stat_list, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            all_gen_stat_list.append(stat_list)
            if Blockdata["myid"] == Blockdata["main_node"]:
                log_main.add(
                    "================================================================================================================"
                )
                log_main.add(
                    "                                    SORT3D IN-DEPTH   generation %d completed"
                    % igen
                )
                log_main.add(
                    "----------------------------------------------------------------------------------------------------------------"
                )

    copy_results(log_main, all_gen_stat_list)  # all nodes function
    compute_final_map(Tracker["constants"]["masterdir"])
    if Blockdata["myid"] == Blockdata["main_node"]:
        time_of_sorting_h, time_of_sorting_m = get_time(time_sorting_start)
        log_main.add(
            "SORT3D execution time: %d hours %d minutes."
            % (time_of_sorting_h, time_of_sorting_m)
        )
    return


#  End of various utilities


def main():
    global Tracker, Blockdata
    progname = os.path.basename(sys.argv[0])
    usage = (
        progname
        + " --refinement_dir=masterdir_of_sxmeridien   --output_dir=sort3d_output --mask3D=mask.hdf --focus=binarymask.hdf  --radius=outer_radius "
        + "  --sym=c1  --img_per_grp=img_per_grp  --minimum_grp_size=minimum_grp_size "
    )
    parser = optparse.OptionParser(usage, version=sp_global_def.SPARXVERSION)
    parser.add_option(
        "--refinement_dir",
        type="string",
        default="",
        help="sxmeridien 3-D refinement directory",
    )
    parser.add_option(
        "--instack",
        type="string",
        default="",
        help="file name, data stack for sorting provided by user. It applies when sorting starts from a given data stack",
    )

    initiate_from_meridien_mode = False
    for q in sys.argv[1:]:
        if q[:16] == "--refinement_dir":
            initiate_from_meridien_mode = True
            break

    initiate_from_data_stack_mode = False
    for q in sys.argv[1:]:
        if q[:9] == "--instack":
            initiate_from_data_stack_mode = True
            break

    # priority
    if initiate_from_data_stack_mode and initiate_from_meridien_mode:
        initiate_from_data_stack_mode = False

    if (not initiate_from_data_stack_mode) and (not initiate_from_meridien_mode):
        if Blockdata["myid"] == Blockdata["main_node"]:
            print(
                "Specify one of the two options to start the program: --refinement_dir, --instack"
            )

    if initiate_from_meridien_mode:
        parser.add_option(
            "--output_dir",
            type="string",
            default="",
            help="Sort3d output directory name",
        )
        parser.add_option(
            "--niter_for_sorting",
            type="int",
            default=-1,
            help="User specified iteration number of 3D refinement for sorting",
        )
        parser.add_option(
            "--focus",
            type="string",
            default="",
            help="Focus 3D mask. File path of a binary 3D mask for focused clustering ",
        )
        parser.add_option(
            "--mask3D",
            type="string",
            default="",
            help="3D mask. File path of the global 3D mask for clustering",
        )
        parser.add_option(
            "--radius",
            type="int",
            default=-1,
            help="Estimated protein radius in pixels",
        )
        parser.add_option(
            "--sym", type="string", default="c1", help="Point-group symmetry"
        )
        parser.add_option(
            "--img_per_grp", type="int", default=1000, help="Number of images per group"
        )
        parser.add_option(
            "--img_per_grp_split_rate",
            type="int",
            default=1,
            help="rate for splitting img_per_grp",
        )
        parser.add_option(
            "--nsmear",
            type="float",
            default=-1.0,
            help="Number of smears used in sorting. Fill it with 1 if user does not want to use all smears",
        )
        parser.add_option(
            "--minimum_grp_size", type="int", default=-1, help="Cluster selection size"
        )
        parser.add_option(
            "--depth_order",
            type="int",
            default=2,
            help="Depth order. A number defines the number of initial independent MGSKmeans runs (2^depth_order)",
        )
        parser.add_option(
            "--memory_per_node",
            type="float",
            default=-1.0,
            help="Memory_per_node, the number used for computing the CPUs/NODE settings given by user",
        )
        parser.add_option(
            "--orientation_groups",
            type="int",
            default=100,
            help="Number of orientation groups in the asymmetric unit",
        )
        parser.add_option(
            "--not_include_unaccounted",
            action="store_true",
            default=False,
            help="Do not reconstruct unaccounted elements in each generation",
        )
        parser.add_option(
            "--stop_mgskmeans_percentage",
            type="float",
            default=10.0,
            help="Swap ratio. A float number between 0.0 and 50",
        )
        parser.add_option(
            "--swap_ratio",
            type="float",
            default=1.0,
            help="Randomness ratio of swapping accounted elements with unaccounted elemetns per cluster",
        )
        parser.add_option(
            "--notapplybckgnoise",
            action="store_true",
            default=False,
            help="Do not applynoise",
        )
        parser.add_option(
            "--do_swap_au",
            action="store_true",
            default=False,
            help="Flag to turn on swapping the accounted for images with the unaccounted for images",
        )
        parser.add_option(
            "--random_group_elimination_threshold",
            type="float",
            default=2.0,
            help="Number of random group reproducibility standard deviation for eliminating random groups",
        )
        (options, args) = parser.parse_args(sys.argv[1:])
        ### Sanity check

        checking_flag = 0
        if Blockdata["myid"] == Blockdata["main_node"]:
            if options.refinement_dir != "":
                if not os.path.exists(options.refinement_dir):
                    checking_flag = 1
        checking_flag = sp_utilities.bcast_number_to_all(
            checking_flag, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        if checking_flag == 1:
            sp_global_def.ERROR(
                "The specified refinement_dir does not exist",
                "sort3d",
                1,
                Blockdata["myid"],
            )

        if options.focus != "":
            if Blockdata["myid"] == Blockdata["main_node"]:
                if not os.path.exists(options.focus):
                    checking_flag = 1
            checking_flag = sp_utilities.bcast_number_to_all(
                checking_flag, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            if checking_flag == 1:
                sp_global_def.ERROR(
                    "The specified focus mask file does not exist",
                    "sort3d",
                    1,
                    Blockdata["myid"],
                )

        if options.mask3D != "":
            if Blockdata["myid"] == Blockdata["main_node"]:
                if not os.path.exists(options.mask3D):
                    checking_flag = 1
            checking_flag = sp_utilities.bcast_number_to_all(
                checking_flag, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            if checking_flag == 1:
                sp_global_def.ERROR(
                    "The specified mask3D file does not exist",
                    "sort3d",
                    1,
                    Blockdata["myid"],
                )

        if options.img_per_grp <= 1:
            sp_global_def.ERROR(
                "Improperiate input paramter for img_per_grp",
                "sort3d",
                1,
                Blockdata["myid"],
            )
        elif options.img_per_grp < options.minimum_grp_size:
            sp_global_def.ERROR(
                "Img_per_grp should be always larger than minimum_grp_size",
                "sort3d",
                1,
                Blockdata["myid"],
            )

        # --- Fill input parameters into dictionary Constants
        Constants = {}
        Constants["stop_mgskmeans_percentage"] = options.stop_mgskmeans_percentage
        Constants["niter_for_sorting"] = options.niter_for_sorting
        Constants["memory_per_node"] = options.memory_per_node
        Constants["orgstack"] = options.instack
        Constants["masterdir"] = options.output_dir
        Constants["refinement_dir"] = options.refinement_dir
        Constants[
            "random_group_elimination_threshold"
        ] = options.random_group_elimination_threshold

        if options.mask3D == "":
            Constants["mask3D"] = False
        else:
            Constants["mask3D"] = options.mask3D
        if options.focus != "":
            Constants["focus3D"] = options.focus
        else:
            Constants["focus3D"] = False

        Constants["depth_order"] = options.depth_order
        Constants["img_per_grp"] = options.img_per_grp
        Constants["img_per_grp_split_rate"] = options.img_per_grp_split_rate
        Constants["minimum_grp_size"] = options.minimum_grp_size
        Constants["radius"] = options.radius
        Constants["sym"] = options.sym
        Constants["nsmear"] = options.nsmear

        Constants["restart_from_nbox"] = 0  # options.restart_from_nbox
        Constants["restart_from_depth_order"] = -1  # options.restart_from_depth_order
        Constants["restart_from_generation"] = -1  # options.restart_from_generation

        #### options for advanced users
        Constants["relax_oriens"] = False
        Constants["do_swap_au"] = options.do_swap_au
        Constants["swap_ratio"] = options.swap_ratio
        Constants["not_include_unaccounted"] = False
        Constants["final_sharpen"] = True  # options.do_not_combinemaps
        Constants["nxinit"] = -1
        Constants["box_niter"] = 5

        ### Frozen options
        Constants["upscale"] = 0.5  #
        Constants["interpolation"] = "trl"
        Constants[
            "comparison_method"
        ] = "cross"  # options.comparison_method # either cross or eucd
        Constants["symmetry"] = Constants["sym"]
        Constants["CTF"] = True
        Constants["do_not_use_3dmask"] = False

        if options.focus:
            Constants[
                "comparison_method"
            ] = "cross"  # in case of focus3D, cross is used.
        Constants["fuse_freq"] = 45.0  # Now in A, convert to pixels before being used
        Constants[
            "orientation_groups"
        ] = options.orientation_groups  # orientation constrained angle step
        #
        #
        # Create and initialize Tracker dictionary with input options  # State Variables
        Tracker = {}
        Tracker["constants"] = Constants
        Tracker["radius"] = Tracker["constants"]["radius"]
        Tracker["upscale"] = Tracker["constants"]["upscale"]
        Tracker[
            "applyctf"
        ] = (
            False
        )  # Should the data be premultiplied by the CTF.  Set to False for local continuous.
        Tracker["nxinit"] = Tracker["constants"]["nxinit"]
        if options.notapplybckgnoise:
            Tracker["applybckgnoise"] = False
        else:
            Tracker["applybckgnoise"] = True

        ###=====<--options for advanced users:
        Tracker["total_number_of_iterations"] = 25
        Tracker["clean_volumes"] = True  # always true

        ### -----------Orientation constraints
        Tracker["tilt1"] = 0.0
        Tracker["tilt2"] = 180.0
        Tracker["grp_size_relx_ratio"] = 0.98
        Tracker["minimum_ptl_number"] = 20
        ### ------------=====< option for proteins images that have preferred orientations
        # for orientation groups
        if (
            Tracker["constants"]["memory_per_node"] == -1
            or Tracker["constants"]["memory_per_node"] < 32.0
        ):
            Tracker["constants"]["small_memory"] = True
        else:
            Tracker["constants"]["small_memory"] = False

        ## additional check
        Tracker["constants"]["hardmask"] = True
        Tracker["applymask"] = True
        Tracker["constants"]["refinement_method"] = "SPARX"
        Tracker["nosmearing"] = False
        checking_flag = 0  # reset
        Blockdata["fftwmpi"] = True

        try:
            Blockdata["symclass"] = sp_fundamentals.symclass(
                Tracker["constants"]["symmetry"]
            )
            Tracker["constants"]["orientation_groups"] = max(
                4,
                old_div(
                    Tracker["constants"]["orientation_groups"],
                    Blockdata["symclass"].nsym,
                ),
            )
        except:
            pass

        ast = get_angle_step_from_number_of_orien_groups(
            Tracker["constants"]["orientation_groups"]
        )
        Blockdata["ncpuspernode"] = Blockdata["no_of_processes_per_group"]
        Blockdata["nsubset"] = Blockdata["ncpuspernode"] * Blockdata["no_of_groups"]
        create_subgroup()
        create_zero_group()

        #

        #

        continue_from_interuption = 0
        continue_from_interuption = create_masterdir()
        log_main = sp_logger.Logger(sp_logger.BaseLogger_Files())
        log_main.prefix = Tracker["constants"]["masterdir"] + "/"
        if continue_from_interuption == 0:
            if Blockdata["myid"] == Blockdata["main_node"]:
                log_main.add(
                    "================================================================================================================"
                )
                log_main.add("                                 SORT3D IN-DEPTH v1.1")
                log_main.add(
                    "================================================================================================================"
                )
            import_data(log_main)
            print_shell_command(sys.argv, log_main)
            check_3dmask(log_main)
            check_mpi_settings(log_main)
            keepsorting = sort3d_init("initialization", log_main)
            dump_tracker(Tracker["constants"]["masterdir"])
            if keepsorting == 0:
                mpi.mpi_finalize()
                exit()
        sorting_main_mpi(log_main, options.depth_order, options.not_include_unaccounted)
        if Blockdata["myid"] == Blockdata["main_node"]:
            log_main.add(
                "----------------------------------------------------------------------------------------------------------------"
            )
            log_main.add("                                 SORT3D IN-DEPTH finished")
            log_main.add(
                "----------------------------------------------------------------------------------------------------------------"
            )
        mpi.mpi_finalize()
        exit()

    elif initiate_from_data_stack_mode:
        parser.add_option(
            "--nxinit", type="int", default=-1, help="User provided image size"
        )
        parser.add_option(
            "--output_dir",
            type="string",
            default="",
            help="Name of the sort3d directory",
        )
        parser.add_option(
            "--focus",
            type="string",
            default="",
            help="Focus 3D mask. File path of a binary 3D mask for focused clustering ",
        )
        parser.add_option(
            "--mask3D",
            type="string",
            default="",
            help="3D mask. File path of the global 3D mask for clustering",
        )
        parser.add_option(
            "--radius",
            type="int",
            default=-1,
            help="Estimated protein radius in pixels",
        )
        parser.add_option(
            "--sym", type="string", default="c1", help="Point-group symmetry"
        )
        parser.add_option(
            "--img_per_grp", type="int", default=1000, help="Number of images per group"
        )
        parser.add_option(
            "--img_per_grp_split_rate",
            type="int",
            default=1,
            help="rate for splitting img_per_grp",
        )
        parser.add_option(
            "--nsmear",
            type="float",
            default=-1.0,
            help="Number of smears used in sorting. Fill it with 1 if user does not want to use all smears",
        )
        parser.add_option(
            "--minimum_grp_size", type="int", default=-1, help="Cluster selection size"
        )
        parser.add_option(
            "--depth_order",
            type="int",
            default=2,
            help="Depth order. A number defines the number of initial independent MGSKmeans runs (2^depth_order)",
        )
        parser.add_option(
            "--memory_per_node",
            type="float",
            default=-1.0,
            help="Memory_per_node, the number used for computing the CPUs/NODE settings given by user",
        )
        parser.add_option(
            "--orientation_groups",
            type="int",
            default=100,
            help="Number of orientation groups in the asymmetric unit",
        )
        parser.add_option(
            "--not_include_unaccounted",
            action="store_true",
            default=False,
            help="Do not reconstruct unaccounted elements in each generation",
        )
        parser.add_option(
            "--stop_mgskmeans_percentage",
            type="float",
            default=10.0,
            help="Swap ratio. A float number between 0.0 and 50.0",
        )
        parser.add_option(
            "--swap_ratio",
            type="float",
            default=1.0,
            help="Randomness ratio of swapping accounted elements with unaccounted elemetns per cluster",
        )
        parser.add_option(
            "--notapplybckgnoise",
            action="store_true",
            default=False,
            help="Flag to turn off background noise",
        )
        parser.add_option(
            "--do_swap_au",
            action="store_true",
            default=False,
            help="Flag to turn on swapping the accounted for images with the unaccounted for images",
        )
        parser.add_option(
            "--random_group_elimination_threshold",
            type="float",
            default=2.0,
            help="Number of random group reproducibility standard deviation for eliminating random groups",
        )
        (options, args) = parser.parse_args(sys.argv[1:])
        ### Sanity check

        checking_flag = 0
        if options.focus != "":
            if Blockdata["myid"] == Blockdata["main_node"]:
                if not os.path.exists(options.focus):
                    checking_flag = 1
            checking_flag = sp_utilities.bcast_number_to_all(
                checking_flag, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            if checking_flag == 1:
                sp_global_def.ERROR(
                    "The specified focus mask file does not exist",
                    "sort3d",
                    1,
                    Blockdata["myid"],
                )

        if options.mask3D != "":
            if Blockdata["myid"] == Blockdata["main_node"]:
                if not os.path.exists(options.mask3D):
                    checking_flag = 1
            checking_flag = sp_utilities.bcast_number_to_all(
                checking_flag, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            if checking_flag == 1:
                sp_global_def.ERROR(
                    "The specified mask3D file does not exist",
                    "sort3d",
                    1,
                    Blockdata["myid"],
                )

        if options.img_per_grp <= 1:
            sp_global_def.ERROR(
                "Improperiate number for img_per_grp", "sort3d", 1, Blockdata["myid"]
            )
        elif options.img_per_grp < options.minimum_grp_size:
            sp_global_def.ERROR(
                "Parameter img_per_grp should be always larger than parameter minimum_grp_size",
                "sort3d",
                1,
                Blockdata["myid"],
            )

        # --- Fill input parameters into dictionary Constants
        Constants = {}
        Constants["stop_mgskmeans_percentage"] = options.stop_mgskmeans_percentage
        Constants["memory_per_node"] = options.memory_per_node
        Constants["orgstack"] = options.instack
        Constants["masterdir"] = options.output_dir
        Constants[
            "random_group_elimination_threshold"
        ] = options.random_group_elimination_threshold

        if options.mask3D == "":
            Constants["mask3D"] = False
        else:
            Constants["mask3D"] = options.mask3D
        if options.focus != "":
            Constants["focus3D"] = options.focus
        else:
            Constants["focus3D"] = False

        Constants["nsmear"] = 1
        Constants["depth_order"] = options.depth_order
        Constants["img_per_grp"] = options.img_per_grp
        Constants["img_per_grp_split_rate"] = options.img_per_grp_split_rate
        Constants["minimum_grp_size"] = options.minimum_grp_size
        Constants["radius"] = options.radius
        Constants["sym"] = options.sym

        Constants["restart_from_nbox"] = 0  # options.restart_from_nbox
        Constants["restart_from_depth_order"] = -1  # options.restart_from_depth_order
        Constants["restart_from_generation"] = -1  # options.restart_from_generation

        #### options for advanced users
        Constants["relax_oriens"] = False
        Constants["do_swap_au"] = options.do_swap_au
        Constants["swap_ratio"] = options.swap_ratio
        Constants["not_include_unaccounted"] = False
        Constants["final_sharpen"] = True  # options.do_not_combinemaps
        Constants["nxinit"] = options.nxinit
        Constants["box_niter"] = 5

        ### Frozen options
        Constants["upscale"] = 0.5  #
        Constants["interpolation"] = "trl"
        Constants[
            "comparison_method"
        ] = "cross"  # options.comparison_method # either cross or eucd
        Constants["symmetry"] = Constants["sym"]
        Constants["CTF"] = True
        Constants["do_not_use_3dmask"] = False

        if options.focus:
            Constants[
                "comparison_method"
            ] = "cross"  # in case of focus3D, cross is used.
        Constants["fuse_freq"] = 45.0  # Now in A, convert to pixels before being used
        Constants[
            "orientation_groups"
        ] = options.orientation_groups  # orientation constrained angle step
        #
        #
        # Create and initialize Tracker dictionary with input options  # State Variables
        Tracker = {}
        Tracker["constants"] = Constants
        Tracker["radius"] = Tracker["constants"]["radius"]
        Tracker["upscale"] = Tracker["constants"]["upscale"]
        Tracker[
            "applyctf"
        ] = (
            False
        )  # Should the data be premultiplied by the CTF.  Set to False for local continuous.
        Tracker["nxinit"] = Tracker["constants"]["nxinit"]
        if options.notapplybckgnoise:
            Tracker["applybckgnoise"] = False
        else:
            Tracker["applybckgnoise"] = True

        ###=====<--options for advanced users:
        Tracker["total_number_of_iterations"] = 25
        Tracker["clean_volumes"] = True  # always true

        ### -----------Orientation constraints
        Tracker["tilt1"] = 0.0
        Tracker["tilt2"] = 180.0
        Tracker["grp_size_relx_ratio"] = 0.98
        Tracker["minimum_ptl_number"] = 20
        ### ------------=====< option for proteins images that have preferred orientations
        # for orientation groups
        if (
            Tracker["constants"]["memory_per_node"] == -1
            or Tracker["constants"]["memory_per_node"] < 32.0
        ):
            Tracker["constants"]["small_memory"] = True
        else:
            Tracker["constants"]["small_memory"] = False

        ## additional check
        Tracker["constants"]["hardmask"] = True
        Tracker["applymask"] = True
        Tracker["constants"]["refinement_method"] = "stack"
        Tracker["constants"]["refinement_dir"] = None
        Tracker["paramstructure_dir"] = None
        Tracker["refang"] = None
        Tracker["rshifts"] = None
        Tracker["paramstructure_dict"] = None
        Tracker["constants"]["selected_iter"] = -1
        Tracker["nosmearing"] = True

        checking_flag = 0  # reset
        Blockdata["fftwmpi"] = True
        Blockdata["symclass"] = sp_fundamentals.symclass(
            Tracker["constants"]["symmetry"]
        )
        Tracker["constants"]["orientation_groups"] = max(
            4,
            old_div(
                Tracker["constants"]["orientation_groups"], Blockdata["symclass"].nsym
            ),
        )

        ast = get_angle_step_from_number_of_orien_groups(
            Tracker["constants"]["orientation_groups"]
        )
        Blockdata["ncpuspernode"] = Blockdata["no_of_processes_per_group"]
        Blockdata["nsubset"] = Blockdata["ncpuspernode"] * Blockdata["no_of_groups"]
        create_subgroup()

        #

        continue_from_interuption = 0
        continue_from_interuption = create_masterdir()
        log_main = sp_logger.Logger(sp_logger.BaseLogger_Files())
        log_main.prefix = Tracker["constants"]["masterdir"] + "/"

        if continue_from_interuption == 0:
            if Blockdata["myid"] == Blockdata["main_node"]:
                log_main.add(
                    "================================================================================================================"
                )
                log_main.add("                                  SORT3D IN-DEPTH v1.1")
                log_main.add(
                    "================================================================================================================\n"
                )
            import_data(log_main)
            print_shell_command(sys.argv, log_main)
            check_3dmask(log_main)
            check_mpi_settings(log_main)
            keepsorting = sort3d_init("Initialization", log_main)
            dump_tracker(Tracker["constants"]["masterdir"])
            if keepsorting == 0:
                mpi.mpi_finalize()
                exit()
        sorting_main_mpi(log_main, options.depth_order, options.not_include_unaccounted)
        if Blockdata["myid"] == Blockdata["main_node"]:
            log_main.add(
                "----------------------------------------------------------------------------------------------------------------"
            )
            log_main.add("                                 SORT3D IN-DEPTH finished")
            log_main.add(
                "----------------------------------------------------------------------------------------------------------------"
            )
        mpi.mpi_finalize()
        exit()


if __name__ == "__main__":
    main()
