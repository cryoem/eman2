#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
from past.utils import old_div
#
# Author: Pawel A.Penczek 05/27/2009 (Pawel.A.Penczek@uth.tmc.edu)
# Please do not copy or modify this file without written consent of the author.
# Copyright (c) 2000-2019 The University of Texas - Houston Medical School
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
#  09/09/2016
#
#  CPU subgroup
#  10/27/2016  Added sigma2 updating in the first phased called PRIMARY
#  11/07       Shared refvol
#  10/28/2016 - Polar
#  11/18/2016 change in strategy
#  04/10/2017 - Enabled for one node
#  04/18/2017 - Introduce symclass to handle angles in a unified manner
#  01/21/2018 - Rationalize the code, particularly restart
#  07/21/2018 - Full size reconstruction after delta change
import EMAN2_cppwrap
import json
import math
import mpi
import numpy
import operator
import optparse
import os
import random
from ..libpy import sp_alignment
from ..libpy import sp_applications
from ..libpy import sp_filter
from ..libpy import sp_fundamentals
from ..libpy.prior_calculation import sp_helix_fundamentals
from ..libpy import sp_global_def
from ..libpy.prior_calculation import sp_helix_sphire
from ..libpy import sp_logger
from ..libpy import sp_morphology
from ..libpy import sp_projection
from ..libpy import sp_reconstruction
from ..libpy import sp_statistics
from ..libpy import sp_user_functions
from ..libpy import sp_utilities
import string
import subprocess
import sys
import time
from builtins import range
import ctypes


#todo: BUG TO FIX 29 april 2020
"""
t2 unresolved reference :  https://github.com/cryoem/eman2/blob/python3_transition/sphire/bin/sp_meridien_alpha.py#L1059
ndoinit unresolved reference :  https://github.com/cryoem/eman2/blob/python3_transition/sphire/bin/sp_meridien_alpha.py#L2543
ndoinit unresolved reference :  https://github.com/cryoem/eman2/blob/python3_transition/sphire/bin/sp_meridien_alpha.py#L3423
final_iter unresolved reference :  https://github.com/cryoem/eman2/blob/python3_transition/sphire/bin/sp_meridien_alpha.py#L9475
final_iter unresolved reference :  https://github.com/cryoem/eman2/blob/python3_transition/sphire/bin/sp_meridien_alpha.py#L9482
"""


"""
There are four ways to run the program:

1. Standard default run, starts from exhaustive searches, uses initial reference structure
mpirun -np 64 --hostfile four_nodes.txt  sxmeridien.py  bdb:sparx_stack vton1 mask15.hdf --sym=c5  --initialshifts  --radius=120  --mask3D=mask15.hdf    >1ovotn &

2. Restart after the last fully finished iteration, one can change some parameters (MPI settings have to be the same)
mpirun -np 64 --hostfile four_nodes.txt  sxmeridien.py  vton1 --radius=100 >2ovotn &

3. Local refinement, starts from user-provided orientation parameters, delta has to be <= 3.75
mpirun -np 64 --hostfile four_nodes.txt sxmeridien.py --local_refinement bdb:sparx_stack   vton3 --delta=1.875 --xr=2.0  --inires=5.5  --sym=c5  --radius=120  --mask3D=mask15.hdf >5ovotn &

4. Restart of local refinement after the last fully finished iteration, one can change some parameters (MPI settings have to be the same)
mpirun -np 64 --hostfile four_nodes.txt  sxmeridien.py --local_refinement  vton3  --xr=0.6 >6ovotn &


"""


"""
Exhaustive:
          delta        projdir         projdir*npsi
         15.00000        189            4536  
          7.50000        756           36288  
          3.75000       3025          290400  
          1.87500      12101         2323392  
          0.93750      48405        18587520  
          0.46875     193622       148701696  
          0.23438     774489      1189615104  
          0.11719    3097958      9516926976  

Assume image 192x192
         15.00000        189       0.028 GB            4536          0.669 GB 
          7.50000        756       0.111 GB           36288          5.351 GB 
          3.75000       3025       0.446 GB          290400         42.821 GB 
          1.87500      12101       1.784 GB         2323392        342.598 GB 
          0.93750      48405       7.138 GB        18587520       2740.841 GB 
          0.46875     193622      28.551 GB       148701696      21926.957 GB 
          0.23438     774489     114.203 GB      1189615104     175415.885 GB 
          0.11719    3097958     456.812 GB      9516926976    1403327.984 GB 


Local with an = 6*delta
         15.00000         96            1248  
          7.50000        113            1469  
          3.75000        115            1495  
          1.87500        115            1495  
          0.93750        117            1521  
          0.46875        117            1521  
          0.23438        117            1521  
          0.11719        115            1495  

Local with an = 12*delta
         15.00000        189            4725  
          7.50000        377            9425  
          3.75000        446           11150  
          1.87500        461           11525  
          0.93750        470           11750  
          0.46875        464           11600  
          0.23438        463           11575  
          0.11719        463           11575  


"""
"""
08/14/2018
Normalization issues:
#This is for real space filter
nx = 1024
p1 = model_gauss_noise(1.0,nx,nx)
mask = Util.unrollmask(nx,nx)
for j in range(nx//2,nx):  mask[0,j]=1.0
#This is for Fourier valid region
m = Util.unrollmask(nx,nx)
p2 = fft(p1)
fp1=fft(Util.mulnclreal(p2,mask))
Util.innerproduct(p2,p2,m)/(nx*nx/2) = Util.innerproduct(fp1,fp1,None)
"""

global Tracker, Blockdata
global target_theta, refang

mpi.mpi_init(0, [])
Tracker = {}
Blockdata = {}
#  MPI stuff
Blockdata["nproc"] = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
Blockdata["myid"] = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
Blockdata["main_node"] = 0
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
    ]  # For 3D stuff take two last nodes
else:
    Blockdata["node_volume"] = [0, 0]
#  We need two CPUs for processing of volumes, they are taken to be main CPUs on each volume
#  We have to send the two myids to all nodes so we can identify main nodes on two selected groups.
Blockdata["main_shared_nodes"] = [
    Blockdata["node_volume"][0] * Blockdata["no_of_processes_per_group"],
    Blockdata["node_volume"][1] * Blockdata["no_of_processes_per_group"],
]
# end of Blockdata
sp_global_def.BATCH = True
sp_global_def.MPI = True


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
    # if( Blockdata["myid"] == Blockdata["main_node"] ): sxprint(submyids)
    world_group = mpi.mpi_comm_group(mpi.MPI_COMM_WORLD)
    subgroup = mpi.mpi_group_incl(world_group, len(submyids), submyids)
    # print(" XXX world group  ",Blockdata["myid"],world_group,subgroup)
    Blockdata["subgroup_comm"] = mpi.mpi_comm_create(mpi.MPI_COMM_WORLD, subgroup)
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    # print(" ZZZ subgroup  ",Blockdata["myid"],world_group,subgroup,subgroup_comm)

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
    # if( Blockdata["myid"] == Blockdata["main_node"] ): sxprint(submyids)
    world_group = mpi.mpi_comm_group(mpi.MPI_COMM_WORLD)
    subgroup = mpi.mpi_group_incl(world_group, len(submyids), submyids)
    # print(" XXX world group  ",Blockdata["myid"],world_group,subgroup)
    Blockdata["group_zero_comm"] = mpi.mpi_comm_create(mpi.MPI_COMM_WORLD, subgroup)
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    # print(" ZZZ subgroup  ",Blockdata["myid"],world_group,subgroup,subgroup_comm)

    Blockdata["group_zero_size"] = -1
    Blockdata["group_zero_myid"] = -1
    if mpi.MPI_COMM_NULL != Blockdata["group_zero_comm"]:
        Blockdata["group_zero_size"] = mpi.mpi_comm_size(Blockdata["group_zero_comm"])
        Blockdata["group_zero_myid"] = mpi.mpi_comm_rank(Blockdata["group_zero_comm"])
    #  "nodes" are zero nodes on subgroups on the two "node_volume" that compute backprojection
    # Blockdata["nodes"] = [Blockdata["node_volume"][0]*Blockdata["ncpuspernode"], Blockdata["node_volume"][1]*Blockdata["ncpuspernode"]]
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    return


# if( Blockdata["subgroup_myid"] > -1 ):
# 	dudu = [Blockdata["subgroup_myid"]]
# 	dudu = wrap_mpi_gatherv(dudu, 0, Blockdata["subgroup_comm"])
# 	if Blockdata["subgroup_myid"] == 0 :  sxprint("  HERE  ",dudu)

# we may want to free it in order to use different number of CPUs
#  create_subgroup()
# if( Blockdata["subgroup_myid"] > -1 ): mpi_comm_free(Blockdata["subgroup_comm"])


def AI(fff, anger, shifter, chout=False):
    global Tracker, Blockdata
    #  chout - if true, one can print, call the program with, chout = (Blockdata["myid"] == Blockdata["main_node"])
    #  fff (fsc), anger, shifter are coming from the previous iteration
    #
    #  Possibilities we will consider:
    #    1.  resolution improved: keep going with current settings.
    #    2.  resolution stalled and no pwadjust: turn on pwadjust
    #    3.  resolution stalled and pwadjust: move to the next phase
    #    4.  resolution decreased: back off and move to the next phase
    #    5.  All phases tried and nxinit < nnxo: set nxinit == nnxo and run local searches.

    ###  checkconvergence  merged in AI  02/16/2017
    # when the following conditions are all true
    # 1. has_fine_enough_angular_sampling  True  #   Current sampling are fine enough
    # 2. nr_iter_wo_resol_gain >= MAX_NR_ITER_WO_RESOL_GAIN #
    # 3. nr_iter_wo_large_hidden_variable_changes >= MAX_NR_ITER_WO_LARGE_HIDDEN_VARIABLE_CHANGES

    keepgoing = 1
    if Blockdata["myid"] == Blockdata["main_node"]:
        line = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime()) + " =>"

    if Tracker["mainiteration"] == 1:
        Tracker["state"] = "INITIAL"

        inc = Tracker["currentres"]
        if Tracker["large_at_Nyquist"]:
            inc += int(old_div(0.25 * Tracker["constants"]["nnxo"], 2) + 0.5)
        else:
            inc += Tracker["nxstep"]
        Tracker["nxinit"] = min(
            2 * inc, Tracker["constants"]["nnxo"]
        )  # Cannot exceed image size
        Tracker["local"] = False
        Tracker["changed_delta"] = False
        #  Do not use CTF during first iteration
        # Tracker["applyctf"]    = False
        # Tracker["constants"]["best"] = Tracker["mainiteration"]
    else:
        if Tracker["mainiteration"] == 2:
            Tracker["state"] = "PRIMARY"
        l05 = -1
        l01 = -1
        for i in range(len(fff)):
            if fff[i] < 0.5:
                l05 = i - 1
                break
        for i in range(l05 + 1, len(fff)):
            if fff[i] < 0.143:
                l01 = i - 1
                break
        l01 = max(l01, -1)

        if chout:
            sp_global_def.sxprint(
                "  AI: TR[nxstep], TR[currentres], TR[fsc143], l05, l01, fff[Tracker[nxinit]//2-1]:",
                Tracker["nxstep"],
                Tracker["currentres"],
                Tracker["fsc143"],
                l05,
                l01,
                fff[old_div(Tracker["nxinit"], 2) - 1],
            )
        Tracker["nxstep"] = max(Tracker["nxstep"], l01 - l05 + 5)
        if Tracker["state"] == "FINAL" or Tracker["state"] == "RESTRICTED":
            Tracker["large_at_Nyquist"] = (
                fff[old_div(Tracker["nxinit"], 2)] > 0.1
                or fff[old_div(Tracker["nxinit"], 2) - 1] > 0.2
            )
        else:
            Tracker["large_at_Nyquist"] = fff[old_div(Tracker["nxinit"], 2) - 1] > 0.2

        if Tracker["mainiteration"] == 2:
            maxres = Tracker["constants"]["inires"]
            maxres_143 = l01
        else:
            maxres = max(
                l05, 5
            )  # 5 is minimum resolution of the map, could be set by the user
            maxres_143 = l01
        try:
            bestres_143 = Tracker["bestres_143"]
        except:
            Tracker["bestres_143"] = maxres_143
        if (maxres >= Tracker["bestres"]) and (maxres_143 >= Tracker["bestres_143"]):
            Tracker["bestres"] = maxres
            Tracker["bestres_143"] = maxres_143
            Tracker["constants"]["best"] = max(Tracker["mainiteration"] - 1, 3)  #

        if maxres > Tracker["currentres"]:
            Tracker["no_improvement"] = 0
            Tracker["no_params_changes"] = 0
        else:
            Tracker["no_improvement"] += 1

        Tracker["currentres"] = maxres
        Tracker["fsc143"] = maxres_143

        params_changes = (
            anger >= 1.03 * Tracker["anger"] and shifter >= 1.03 * Tracker["shifter"]
        )

        #  figure changes in params
        if chout:
            sp_global_def.sxprint(
                "  Incoming  parameters  %10.3f  %10.3f  %10.3f  %10.3f   %s"
                % (Tracker["anger"], anger, Tracker["shifter"], shifter, params_changes)
            )
        if params_changes:
            Tracker["no_params_changes"] += 1
        else:
            Tracker["no_params_changes"] = 0

        if anger < Tracker["anger"]:
            Tracker["anger"] = anger
        if shifter < Tracker["shifter"]:
            Tracker["shifter"] = shifter

        inc = Tracker["currentres"]
        if Tracker["large_at_Nyquist"]:
            inc += int(old_div(0.25 * Tracker["constants"]["nnxo"], 2) + 0.5)
            slim = int(Tracker["nxinit"] * 1.09)
            tmp = min(max(2 * inc, slim + slim % 2), Tracker["constants"]["nnxo"])
        else:
            inc += Tracker["nxstep"]
            tmp = min(2 * inc, Tracker["constants"]["nnxo"])  # Cannot exceed image size

        if chout:
            sp_global_def.sxprint(
                "  IN AI: nxstep, large at Nyq, outcoming current res, adjusted current, inc, estimated image size",
                Tracker["nxstep"],
                Tracker["large_at_Nyquist"],
                Tracker["currentres"],
                inc,
                tmp,
            )

        Tracker["nxinit"] = tmp
        Tracker["changed_delta"] = False
        #  decide angular step and translations
        if (
            (Tracker["no_improvement"] >= Tracker["constants"]["limit_improvement"])
            and (Tracker["no_params_changes"] >= Tracker["constants"]["limit_changes"])
            and (not Tracker["large_at_Nyquist"])
        ):
            if (
                Tracker["delta"] < 0.75 * Tracker["acc_rot"]
            ):  # <<<----it might cause converge issues when shake is 0.0
                keepgoing = 0
                if Blockdata["myid"] == Blockdata["main_node"]:
                    sp_global_def.sxprint(
                        line,
                        "Convergence criterion A is reached (angular step delta smaller than 3/4 changes in angles))",
                    )
            else:
                step_range, step = compute_search_params(
                    Tracker["acc_trans"], Tracker["shifter"], Tracker["xr"]
                )
                if chout:
                    sp_global_def.sxprint(
                        "  Computed  pares  ",
                        Tracker["anger"],
                        anger,
                        Tracker["shifter"],
                        shifter,
                        Tracker["xr"],
                        step_range,
                        step,
                    )
                Tracker["xr"] = step_range
                Tracker["ts"] = step
                Tracker["delta"] = old_div(Tracker["delta"], 2.0)
                Tracker["changed_delta"] = True
                if (
                    (Tracker["state"] == "PRIMARY")
                    and (Tracker["delta"] <= old_div(3.75, 2.0))
                    and Tracker["constants"]["do_exhaustive"]
                ):
                    Tracker["state"] = "EXHAUSTIVE"
                    Tracker["delta"] *= 2.0
                    Tracker["changed_delta"] = False
                elif Tracker["delta"] <= old_div(3.75, 2.0):  # MOVE DOWN TO RESTRICTED
                    Tracker["an"] = 6 * Tracker["delta"]
                    Tracker["howmany"] = 4
                    Tracker["theta_min"] = -1
                    Tracker["theta_max"] = -1
                    if Tracker["delta"] <= numpy.degrees(
                        math.atan(old_div(0.25, Tracker["constants"]["radius"]))
                    ):
                        Tracker["state"] = "FINAL"
                    else:
                        Tracker["state"] = "RESTRICTED"
                else:
                    Tracker["an"] = -1
                    if Tracker["state"] == "PRIMARY":
                        Tracker["state"] = "EXHAUSTIVE"
                if chout:
                    sp_global_def.sxprint(
                        "  IN AI there was reset due to no changes, adjust stuff  ",
                        Tracker["no_improvement"],
                        Tracker["no_params_changes"],
                        Tracker["delta"],
                        Tracker["xr"],
                        Tracker["ts"],
                        Tracker["state"],
                    )
                # check convergence before reset
                if (Tracker["state"] == "FINAL") and (
                    Tracker["no_improvement"]
                    >= Tracker["constants"]["limit_improvement"]
                ):
                    keepgoing = 0
                    if Blockdata["myid"] == Blockdata["main_node"]:
                        sp_global_def.sxprint(
                            line,
                            "Convergence criterion B is reached (angular step delta smaller than the limit imposed by the structure radius)",
                        )
                Tracker["no_improvement"] = 0
                Tracker["no_params_changes"] = 0
                Tracker["anger"] = 1.0e23
                Tracker["shifter"] = 1.0e23
    Tracker["keepfirst"] = -1
    return keepgoing


def AI_continuation(fff, anger=-1.0, shifter=-1.0, chout=False):
    global Tracker, Blockdata
    #  chout - if true, one can print, call the program with, chout = (Blockdata["myid"] == Blockdata["main_node"])
    #  fff (fsc), anger, shifter are coming from the previous iteration
    #
    #  Possibilities we will consider:
    #    1.  resolution improved: keep going with current settings.
    #    2.  resolution stalled and no pwadjust: turn on pwadjust
    #    3.  resolution stalled and pwadjust: move to the next phase
    #    4.  resolution decreased: back off and move to the next phase
    #    5.  All phases tried and nxinit < nnxo: set nxinit == nnxo and run local searches.

    ###  checkconvergence  merged in AI  02/16/2017
    # when the following conditions are all true
    # 1. has_fine_enough_angular_sampling  True  #   Current sampling are fine enough
    # 2. nr_iter_wo_resol_gain >= MAX_NR_ITER_WO_RESOL_GAIN #
    # 3. nr_iter_wo_large_hidden_variable_changes >= MAX_NR_ITER_WO_LARGE_HIDDEN_VARIABLE_CHANGES

    keepgoing = 1

    l05 = -1
    l01 = -1
    for i in range(len(fff)):
        if fff[i] < 0.5:
            l05 = i - 1
            break
    for i in range(l05 + 1, len(fff)):
        if fff[i] < 0.143:
            l01 = i - 1
            break
    l01 = max(l01, -1)

    if chout:
        sp_global_def.sxprint(
            "  AI: Tracker[nxstep], TR[currentres], Tracker[fsc143], l05, l01, fff[Tracker[nxinit]//2-1]:",
            Tracker["nxstep"],
            Tracker["currentres"],
            Tracker["fsc143"],
            l05,
            l01,
            fff[old_div(Tracker["nxinit"], 2) - 1],
        )
    if Blockdata["myid"] == Blockdata["main_node"]:
        line = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime()) + " =>"

    if Tracker["mainiteration"] == 1:
        # Tracker["state"]		= "PRIMARY"
        Tracker["currentres"] = l05
        Tracker["fsc143"] = l01
        Tracker["large_at_Nyquist"] = (
            fff[old_div(Tracker["nxinit"], 2)] > 0.1
            or fff[old_div(Tracker["nxinit"], 2) - 1] > 0.2
        )
        Tracker["nxinit"] = min(
            2 * Tracker["fsc143"], Tracker["constants"]["nnxo"]
        )  # Cannot exceed image size
        Tracker["local"] = True
        Tracker["an"] = 6 * Tracker["delta"]
        Tracker["no_improvement"] = 0
        Tracker["no_params_changes"] = 0
        Tracker["anger"] = 1.0e23
        Tracker["shifter"] = 1.0e23
        Tracker["constants"]["best"] = Tracker["mainiteration"]
        if Blockdata["myid"] == Blockdata["main_node"]:
            sp_global_def.sxprint(
                line,
                "ITERATION  #%2d. Resolution achieved       : %3d/%3d pixels, %5.2fA/%5.2fA."
                % (
                    Tracker["mainiteration"],
                    Tracker["currentres"],
                    Tracker["fsc143"],
                    old_div(
                        Tracker["constants"]["pixel_size"]
                        * Tracker["constants"]["nnxo"],
                        float(Tracker["currentres"]),
                    ),
                    old_div(
                        Tracker["constants"]["pixel_size"]
                        * Tracker["constants"]["nnxo"],
                        float(Tracker["fsc143"]),
                    ),
                ),
            )
    else:
        if Tracker["mainiteration"] > 3:
            Tracker["nxstep"] = max(Tracker["nxstep"], l01 - l05 + 5)
        if Tracker["state"] == "FINAL" or Tracker["state"] == "RESTRICTED":
            Tracker["large_at_Nyquist"] = (
                fff[old_div(Tracker["nxinit"], 2)] > 0.1
                or fff[old_div(Tracker["nxinit"], 2) - 1] > 0.2
            )
        else:
            Tracker["large_at_Nyquist"] = fff[old_div(Tracker["nxinit"], 2) - 1] > 0.2

        maxres = max(l05, 5)
        maxres_143 = l01
        try:
            bestres_143 = Tracker["bestres_143"]
        except:
            Tracker["bestres_143"] = maxres_143
        if (maxres >= Tracker["bestres"]) and (maxres_143 >= Tracker["bestres_143"]):
            Tracker["bestres"] = maxres
            Tracker["bestres_143"] = maxres_143
            Tracker["constants"]["best"] = max(Tracker["mainiteration"] - 1, 3)

        if maxres > Tracker["currentres"]:
            Tracker["no_improvement"] = 0
            Tracker["no_params_changes"] = 0
        else:
            Tracker["no_improvement"] += 1

        Tracker["currentres"] = maxres
        Tracker["fsc143"] = maxres_143

        params_changes = (
            anger >= 1.03 * Tracker["anger"] and shifter >= 1.03 * Tracker["shifter"]
        )

        #  figure changes in params
        if chout:
            sp_global_def.sxprint(
                "  Incoming  parameters  %10.3f  %10.3f  %10.3f  %10.3f   %s"
                % (Tracker["anger"], anger, Tracker["shifter"], shifter, params_changes)
            )
        if params_changes:
            Tracker["no_params_changes"] += 1
        else:
            Tracker["no_params_changes"] = 0

        if anger < Tracker["anger"]:
            Tracker["anger"] = anger
        if shifter < Tracker["shifter"]:
            Tracker["shifter"] = shifter

        inc = Tracker["currentres"]
        if Tracker["large_at_Nyquist"]:
            inc += int(old_div(0.25 * Tracker["constants"]["nnxo"], 2) + 0.5)
            slim = int(Tracker["nxinit"] * 1.09)
            tmp = min(max(2 * inc, slim + slim % 2), Tracker["constants"]["nnxo"])
        else:
            inc += Tracker["nxstep"]
            tmp = min(2 * inc, Tracker["constants"]["nnxo"])  # Cannot exceed image size

        if chout:
            sp_global_def.sxprint(
                "  IN AI: nxstep, large at Nyq, outcoming current res, adjusted current, inc, estimated image size",
                Tracker["nxstep"],
                Tracker["large_at_Nyquist"],
                Tracker["currentres"],
                inc,
                tmp,
            )

        Tracker["nxinit"] = tmp
        Tracker["changed_delta"] = False
        #  decide angular step and translations
        if (
            (Tracker["no_improvement"] >= Tracker["constants"]["limit_improvement"])
            and (Tracker["no_params_changes"] >= Tracker["constants"]["limit_changes"])
            and (not Tracker["large_at_Nyquist"])
        ):
            if (
                Tracker["delta"] < 0.75 * Tracker["acc_rot"]
            ):  # <<<----it might cause converge issues when shake is 0.0
                keepgoing = 0
                if Blockdata["myid"] == Blockdata["main_node"]:
                    sp_global_def.sxprint(
                        line,
                        "Convergence criterion A is reached (angular step delta smaller than 3/4 changes in angles))",
                    )
            else:
                step_range, step = compute_search_params(
                    Tracker["acc_trans"], Tracker["shifter"], Tracker["xr"]
                )
                if chout:
                    sp_global_def.sxprint(
                        "  Computed  pares  ",
                        Tracker["anger"],
                        anger,
                        Tracker["shifter"],
                        shifter,
                        Tracker["xr"],
                        step_range,
                        step,
                    )
                Tracker["xr"] = step_range
                Tracker["ts"] = step
                Tracker["delta"] = old_div(Tracker["delta"], 2.0)
                Tracker["changed_delta"] = True
                if (
                    Tracker["delta"] <= old_div(3.75, 2.0) or True
                ):  # MOVE DOWN TO RESTRICTED
                    Tracker["an"] = 6 * Tracker["delta"]
                    Tracker["howmany"] = 4
                    Tracker["theta_min"] = -1
                    Tracker["theta_max"] = -1
                    if Tracker["delta"] <= numpy.degrees(
                        math.atan(old_div(0.25, Tracker["constants"]["radius"]))
                    ):
                        Tracker["state"] = "FINAL"
                    else:
                        Tracker["state"] = "RESTRICTED"
                else:
                    Tracker["an"] = -1
                    if Tracker["state"] == "PRIMARY":
                        Tracker["state"] = "EXHAUSTIVE"
                if chout:
                    sp_global_def.sxprint(
                        "  IN AI there was reset due to no changes, adjust stuff  ",
                        Tracker["no_improvement"],
                        Tracker["no_params_changes"],
                        Tracker["delta"],
                        Tracker["xr"],
                        Tracker["ts"],
                        Tracker["state"],
                    )
                # check convergence before reset
                if (Tracker["state"] == "FINAL") and (
                    Tracker["no_improvement"]
                    >= Tracker["constants"]["limit_improvement"]
                ):
                    keepgoing = 0
                    if Blockdata["myid"] == Blockdata["main_node"]:
                        sp_global_def.sxprint(
                            line,
                            "Convergence criterion B is reached (angular step delta smaller than the limit imposed by the structure radius)",
                        )
                Tracker["no_improvement"] = 0
                Tracker["no_params_changes"] = 0
                Tracker["anger"] = 1.0e23
                Tracker["shifter"] = 1.0e23
    Tracker["keepfirst"] = -1
    return keepgoing


def params_changes(params, oldparams):
    #  Indexes contain list of images processed - sorted integers, subset of the full range.
    #  params - contain parameters associated with these images
    #  Both lists can be of different sizes, so we have to find a common subset
    #  We do not compensate for random changes of grids.

    n = len(params)
    anger = 0.0
    shifter = 0.0
    #  The shifter is given in the full scale displacement
    for i in range(n):
        shifter += (params[i][3] - oldparams[i][3]) ** 2 + (
            params[i][4] - oldparams[i][4]
        ) ** 2
        anger += get_anger(
            params[i][0:3], oldparams[i][0:3]
        )  # Symmetry is in Blockdata

    return (
        round(old_div(anger, n), 5),
        round(numpy.sqrt(old_div(old_div(shifter, 2), n)), 5),
    )


def compute_search_params(acc_trans, shifter, old_range):
    # step refer to the fine sampled step; while the range remains
    if old_range == 0.0 and shifter != 0.0:
        old_range = acc_trans
    step = min(1.5, 0.75 * acc_trans)  # remove 2
    new_range = min(1.3 * old_range, 5.0 * shifter)
    new_range = max(new_range, 3.0 * step)  # change 1.5 to 3.0
    if new_range > 8.0 * step:
        new_range = old_div(new_range, 2.0)  # change 4 to 8
    if new_range > 8.0 * step:
        step = old_div(new_range, 8.0)  # change 4 to 8
    #  change 4. to 8. on account of the fact that we return actual shift which is then doubled for coarse search.
    if new_range == 0.0:
        step = 0.5  # change 1.0 to 0.5
    return new_range, step


def assign_particles_to_groups(
    minimum_group_size=10, asubset=None, name_tag="ptcl_source_image"
):
    global Tracker, Blockdata
    #  Input data does not have to be consecutive in terms of ptcl_source_image/filament_id or defocus
    #
    if not asubset:
        try:
            stmp = EMAN2_cppwrap.EMUtil.get_all_attributes(
                Tracker["constants"]["stack"], name_tag
            )
            if Tracker["constants"]["CTF"]:
                defstmp = EMAN2_cppwrap.EMUtil.get_all_attributes(
                    Tracker["constants"]["stack"], "ctf"
                )
            else:
                defstmp = [-1.0] * len(stmp)
            for i in range(len(defstmp)):
                defstmp[i] = round(defstmp[i].defocus, 4)
        except:
            if Tracker["constants"]["CTF"]:
                stmp = EMAN2_cppwrap.EMUtil.get_all_attributes(
                    Tracker["constants"]["stack"], "ctf"
                )
                for i in range(len(stmp)):
                    stmp[i] = round(stmp[i].defocus, 4)
                defstmp = stmp[:]
            else:
                sp_global_def.ERROR(
                    "Either ptcl_source_image/filament/filament_id or ctf has to be present in the header."
                )
    else:
        try:
            stmp_junk = EMAN2_cppwrap.EMUtil.get_all_attributes(
                Tracker["constants"]["stack"], name_tag
            )
            stmp = [None] * len(asubset)
            for isub in range(len(asubset)):
                stmp[isub] = stmp_junk[asubset[isub]]
            if Tracker["constants"]["CTF"]:
                defstmp_junk = EMAN2_cppwrap.EMUtil.get_all_attributes(
                    Tracker["constants"]["stack"], "ctf"
                )
                defstmp = [None] * len(asubset)
                for isub in range(len(asubset)):
                    defstmp[isub] = round(defstmp_junk[asubset[isub]].defocus, 4)
            else:
                defstmp = [-1.0] * len(asubset)
        except:
            if Tracker["constants"]["CTF"]:
                stmp_junk = EMAN2_cppwrap.EMUtil.get_all_attributes(
                    Tracker["constants"]["stack"], "ctf"
                )
                stmp = [None] * len(asubset)
                defstmp = [-1.0] * len(asubset)
                for isub in range(len(asubset)):
                    stmp[isub] = stmp_junk[asubset[isub]]
                    stmp[isub] = round(stmp[isub].defocus, 4)
                defstmp[:] = stmp[:]
            else:
                sp_global_def.ERROR(
                    "Either ptcl_source_image/filament/filament_id or ctf has to be present in the header."
                )
    tt = [[stmp[i], i] for i in range(len(stmp))]
    tt.sort()
    tt.append([-1, -1])
    st = tt[0][0]
    sd = []
    occup = []
    groups = []
    ig = 0
    ib = 0
    for i in range(len(tt)):
        if st != tt[i][0]:
            # create a group
            groups.append([tt[k][1] for k in range(ib, i)])
            sd.append([st, defstmp[tt[ib][1]]])
            occup.append(len(groups[ig]))
            groups[ig].sort()
            ib = i
            st = tt[i][0]
            ig += 1
    del tt, stmp, defstmp
    # print(" UUU  ", sd)
    #  [0]ID, [1]stamp, [2]defocus, [3]occupancy, [4]groups
    cross_reference_txt = [
        [[i] for i in range(len(sd))],
        [sd[i][0] for i in range(len(sd))],
        [sd[i][1] for i in range(len(sd))],
        [occup[i] for i in range(len(sd))],
        [groups[i] for i in range(len(sd))],
    ]
    del occup, groups

    #  Remove small groups
    while min(cross_reference_txt[3]) < minimum_group_size:
        # print("  minimum occupancy ",min(cross_reference_txt[3]),len(cross_reference_txt[3]))
        #  Find smallest group
        lax = minimum_group_size
        for i in range(len(cross_reference_txt[3])):
            if lax > cross_reference_txt[3][i]:
                lax = cross_reference_txt[3][i]
                togo = i
        if Tracker["constants"]["CTF"]:
            # find nearest group by defocus
            sdef = 1.0e23
            for i in range(len(cross_reference_txt[3])):
                if i != togo:
                    qt = abs(cross_reference_txt[2][i] - cross_reference_txt[2][togo])
                    if qt < sdef:
                        target = i
                        sdef = qt
        else:
            # find the next smallest
            lax = minimum_group_size
            for i in range(len(cross_reference_txt[3])):
                if i != togo:
                    if lax > cross_reference_txt[3][i]:
                        lax = cross_reference_txt[3][i]
                        target = i

        # print("  merging groups  ",target,togo,cross_reference_txt[2][target],cross_reference_txt[2][togo],cross_reference_txt[3][target],cross_reference_txt[3][togo],len(cross_reference_txt[4][target]),len(cross_reference_txt[4][togo]))
        cross_reference_txt[2][target] = cross_reference_txt[2][target] * sum(
            cross_reference_txt[0][target]
        ) + cross_reference_txt[2][togo] * sum(cross_reference_txt[0][togo])
        cross_reference_txt[0][target] += cross_reference_txt[0][togo]
        cross_reference_txt[2][target] = old_div(
            cross_reference_txt[2][target], sum(cross_reference_txt[0][target])
        )
        cross_reference_txt[3][target] += cross_reference_txt[3][togo]
        cross_reference_txt[4][target] += cross_reference_txt[4][togo]
        # print("  merged  ",cross_reference_txt[0][target],cross_reference_txt[3][target],len(cross_reference_txt[4][target]))

        #  remove the group
        for i in range(len(cross_reference_txt)):
            del cross_reference_txt[i][togo]

    #  Sort as much as possible by the original particle number
    for i in range(len(cross_reference_txt[4])):
        cross_reference_txt[4][i].sort()

    temp = [
        [i, cross_reference_txt[4][i][0]] for i in range(len(cross_reference_txt[0]))
    ]

    temp.sort(key=operator.itemgetter(1))

    cross_reference_txt = [
        [cross_reference_txt[j][temp[i][0]] for i in range(len(cross_reference_txt[0]))]
        for j in range(5)
    ]

    sp_utilities.write_text_row(
        cross_reference_txt[0],
        os.path.join(Tracker["constants"]["masterdir"], "main000", "groupids.txt"),
    )
    sp_utilities.write_text_row(
        [
            [
                sd[cross_reference_txt[0][i][j]][0]
                for j in range(len(cross_reference_txt[0][i]))
            ]
            for i in range(len(cross_reference_txt[0]))
        ],
        os.path.join(Tracker["constants"]["masterdir"], "main000", "micids.txt"),
    )

    Tracker["constants"]["number_of_groups"] = len(cross_reference_txt[0])
    #  split into two chunks by groups
    lili = [[], list(range(Tracker["constants"]["number_of_groups"]))]
    random.shuffle(lili[1])
    lili[0] = lili[1][: old_div(len(lili[1]), 2)]
    lili[1] = lili[1][old_div(len(lili[1]), 2) :]
    lili[0].sort()
    lili[1].sort()

    #  Create output tables
    for iproc in range(2):
        sp_utilities.write_text_row(
            [cross_reference_txt[0][i] for i in lili[iproc]],
            os.path.join(
                Tracker["constants"]["masterdir"],
                "main000",
                "groupids_%03d.txt" % iproc,
            ),
        )
        sp_utilities.write_text_row(
            [
                [
                    sd[cross_reference_txt[0][i][j]][0]
                    for j in range(len(cross_reference_txt[0][i]))
                ]
                for i in lili[iproc]
            ],
            os.path.join(
                Tracker["constants"]["masterdir"], "main000", "micids_%03d.txt" % iproc
            ),
        )
    del sd

    sp_utilities.write_text_file(
        [len(cross_reference_txt[4][i]) for i in range(len(cross_reference_txt[4]))],
        os.path.join(
            Tracker["constants"]["masterdir"],
            "main000",
            "number_of_particles_per_group.txt",
        ),
    )

    q0 = []
    g0 = []
    q1 = []
    g1 = []
    for i in lili[0]:
        g0 += [i] * len(cross_reference_txt[4][i])
        q0 += cross_reference_txt[4][i]
    for i in lili[1]:
        g1 += [i] * len(cross_reference_txt[4][i])
        q1 += cross_reference_txt[4][i]
    Tracker["nima_per_chunk"] = [len(q0), len(q1)]

    ### conversion

    # for iproc in range(2):
    # 	if( Tracker["nima_per_chunk"][iproc] < Blockdata["nproc"] ):  ERROR("Number of particles per chunk smaller than the number of CPUs","assign_particles_to_groups",1,Blockdata["myid"])
    # write_text_file(q0, os.path.join(Tracker["constants"]["masterdir"],"main000","tchunk_0.txt") )
    sp_utilities.write_text_file(
        g0,
        os.path.join(
            Tracker["constants"]["masterdir"], "main000", "particle_groups_0.txt"
        ),
    )
    # write_text_file(q1, os.path.join(Tracker["constants"]["masterdir"],"main000","tchunk_1.txt") )
    sp_utilities.write_text_file(
        g1,
        os.path.join(
            Tracker["constants"]["masterdir"], "main000", "particle_groups_1.txt"
        ),
    )
    if asubset:
        sq0 = [None] * len(q0)
        sq1 = [None] * len(q1)
        for iptl in range(len(sq0)):
            sq0[iptl] = asubset[q0[iptl]]
        for iptl in range(len(sq1)):
            sq1[iptl] = asubset[q1[iptl]]
        return sq0, sq1
    else:
        return q0, q1


#  CONES support functions


def number_of_cones_to_delta(number_of_cones):
    if number_of_cones == 1:
        return Blockdata["symclass"][1][3] + 0.5, 1
    else:
        if Blockdata["symclass"].sym[0] == "c":
            t2 = 89.0
            for i in range(92, 0, -1):
                a = Blockdata["symclass"].even_angles(i, theta1=1.0, theta2=t2)
                a += [[(q[0] + 90.0) % 360.0, 180.0 - q[1], 0] for q in a]
                nren = len(a)
                if nren > number_of_cones:
                    return float(i), nren
            while True:
                i = old_div(i, 2.0)
                a = Blockdata["symclass"].even_angles(i, theta1=1.0, theta2=t2)
                a += [[(q[0] + 90.0) % 360.0, 180.0 - q[1], 0] for q in a]
                nren = len(a)
                if nren > number_of_cones:
                    return float(i), nren

        else:
            # t2 = Blockdata["symclass"].brackets[1][3] + 0.5
            for i in range(int(t2 + 1), 0, -1):
                nren = len(Blockdata["symclass"].even_angles(i, theta1=1.0))
                if nren > number_of_cones:
                    return float(i), nren
            while True:
                i = old_div(i, 2.0)
                nren = len(Blockdata["symclass"].even_angles(i, theta1=1.0))
                if nren > number_of_cones:
                    return float(i), nren

        sp_global_def.ERROR("number_of_cones_to_delta", "should not be here", 0)
        return -1, 0


def find_assignments_of_refangles_to_angles(normals_set, ancor_angle, an):
    #  returns list of angles within normals_set that are within an angles from ancor_angle
    global Blockdata

    nang1 = len(normals_set) - 1
    Blockdata["target_theta"] = ancor_angle[1]
    u1, u2 = sp_fundamentals.goldsearch(auxiliary_funcdef, 0, nang1, tol=1.0e-2)
    u1 = int(u1 + 0.5)

    Blockdata["target_theta"] = max(0.0, ancor_angle[1] - 1.2 * an)
    u3, u4 = sp_fundamentals.goldsearch(auxiliary_funcdef, 0, nang1, tol=1.0e-2)
    u3 = int(u3 + 0.5)
    if u3 < 10:
        u3 = 0

    Blockdata["target_theta"] = min(
        Blockdata["symclass"].brackets[1][3], ancor_angle[1] + 1.2 * an
    )
    u5, u6 = sp_fundamentals.goldsearch(auxiliary_funcdef, 0, nang1, tol=1.0e-2)
    u5 = int(u5 + 0.5)
    if u5 > nang1 - 10:
        u5 = nang1 + 1

    ancordir = sp_utilities.angles_to_normals(
        Blockdata["symclass"].symmetry_neighbors([ancor_angle[:3]])
    )
    ltemp = EMAN2_cppwrap.Util.cone_dirs_f(normals_set[u3:u5], ancordir, an)
    ###ltemp = cone_dirs_f( normals_set[u3:u5], ancordir, an )#Util.cone_dirs_f( normals_set[u3:u5], ancordir, an )
    # print(" us ",Blockdata["myid"],u1,u3,u5,m,"  ltemp ",ltemp)
    return [qtemp + u3 for qtemp in ltemp]


def find_nearest_k_refangles_to_many_angles(normals_set, angles, delta, howmany):
    assignments = [-1] * len(angles)
    for i, ancor_angle in enumerate(angles):
        assignments[i] = find_nearest_k_refangles_to_angles(
            normals_set, ancor_angle, delta, howmany
        )
    return assignments


def find_nearest_k_refangles_to_angles(normals_set, ancor_angle, delta, howmany):
    #  returns list of angles within normals_set that are within an angles from ancor_angle
    global Blockdata
    nang1 = len(normals_set) - 1
    qtl = 0
    bigger = 1.0
    while qtl < howmany:
        an = bigger * max(delta * (old_div(howmany, 2)), delta)
        Blockdata["target_theta"] = ancor_angle[1]
        u1, u2 = sp_fundamentals.goldsearch(auxiliary_funcdef, 0, nang1, tol=1.0e-2)
        u1 = int(u1 + 0.5)

        Blockdata["target_theta"] = max(0.0, ancor_angle[1] - 1.1 * an)
        u3, u4 = sp_fundamentals.goldsearch(auxiliary_funcdef, 0, nang1, tol=1.0e-2)
        u3 = int(u3 + 0.5)
        if u3 < 10:
            u3 = 0

        Blockdata["target_theta"] = min(
            Blockdata["symclass"].brackets[1][3], ancor_angle[1] + 1.1 * an
        )
        u5, u6 = sp_fundamentals.goldsearch(auxiliary_funcdef, 0, nang1, tol=1.0e-2)
        u5 = int(u5 + 0.5)
        if u5 > nang1 - 10:
            u5 = nang1 + 1
        qtl = len(normals_set[u3:u5])
        bigger *= 1.25

    if Blockdata["symclass"].sym == "c1":
        ancordir = sp_utilities.getfvec(ancor_angle[0], ancor_angle[1])
        ltemp = EMAN2_cppwrap.Util.nearest_fang_select(
            normals_set[u3:u5], ancordir[0], ancordir[1], ancordir[2], howmany
        )
    else:
        ancordir = sp_utilities.angles_to_normals(
            Blockdata["symclass"].symmetry_neighbors([ancor_angle])
        )
        ltemp = EMAN2_cppwrap.Util.nearest_fang_sym(
            ancordir, normals_set[u3:u5], len(ancordir), howmany
        )
    return [qtemp + u3 for qtemp in ltemp]


def auxiliary_funcdef(xxx):
    global Blockdata
    l = min(max(int(xxx + 0.5), 0), len(Blockdata["angle_set"]) - 1)
    # print l,abs(target_theta - refang[l][1])
    return abs(Blockdata["target_theta"] - Blockdata["angle_set"][l][1])


#  ------


def compute_sigma(projdata, params, first_procid, dryrun=False, myid=-1, mpi_comm=-1):
    global Tracker, Blockdata
    # Input stack of particles with all params in header
    # Output: 1/sigma^2 and a dictionary
    #  It could be a parameter
    if mpi_comm < 0:
        mpi_comm = mpi.MPI_COMM_WORLD
    npad = 1

    if dryrun:
        # tsd = model_blank(nv + nv//2,len(sd), 1, 1.0)
        # tocp = model_blank(len(sd), 1, 1, 1.0)
        if myid == Blockdata["main_node"]:
            tsd = sp_utilities.get_im(
                os.path.join(Tracker["previousoutputdir"], "bckgnoise.hdf")
            )
            tsd.write_image(os.path.join(Tracker["directory"], "bckgnoise.hdf"))
            nnx = tsd.get_xsize()
            nny = tsd.get_ysize()
        else:
            nnx = 0
            nny = 0
        nnx = sp_utilities.bcast_number_to_all(
            nnx, source_node=Blockdata["main_node"], mpi_comm=mpi_comm
        )
        nny = sp_utilities.bcast_number_to_all(
            nny, source_node=Blockdata["main_node"], mpi_comm=mpi_comm
        )
        if myid != Blockdata["main_node"]:
            tsd = sp_utilities.model_blank(nnx, nny, 1, 1.0)
        sp_utilities.bcast_EMData_to_all(
            tsd, myid, source_node=Blockdata["main_node"], comm=mpi_comm
        )
        """Multiline Comment0"""

    else:

        if myid == Blockdata["main_node"]:
            ngroups = len(
                sp_utilities.read_text_file(
                    os.path.join(
                        Tracker["constants"]["masterdir"], "main000", "groupids.txt"
                    )
                )
            )
        else:
            ngroups = 0
        ngroups = sp_utilities.bcast_number_to_all(
            ngroups, source_node=Blockdata["main_node"], mpi_comm=mpi_comm
        )

        ndata = len(projdata)
        nx = Tracker["constants"]["nnxo"]
        mx = npad * nx
        nv = old_div(mx, 2) + 1
        """Multiline Comment1"""

        mask = sp_utilities.model_circle(Tracker["constants"]["radius"], nx, nx)
        tsd = sp_utilities.model_blank(nv + old_div(nv, 2), ngroups)

        # projdata, params = getalldata(partstack, params, myid, Blockdata["nproc"])
        """Multiline Comment2"""
        if Blockdata["accumulatepw"] == None:
            Blockdata["accumulatepw"] = [[], []]
            doac = True
        else:
            doac = False
        tocp = sp_utilities.model_blank(ngroups)
        tavg = sp_utilities.model_blank(nx, nx)
        for i in range(
            ndata
        ):  # apply_shift; info_mask; norm consistent with get_shrink_data
            indx = projdata[i].get_attr("particle_group")
            phi, theta, psi, sx, sy = (
                params[i][0],
                params[i][1],
                params[i][2],
                params[i][3],
                params[i][4],
            )
            stmp = sp_fundamentals.cyclic_shift(
                projdata[i], int(round(sx)), int(round(sy))
            )
            st = get_image_statistics(stmp, mask, False)
            stmp -= st[0]
            stmp = old_div(stmp, st[1])
            temp = sp_morphology.cosinemask(
                stmp, radius=Tracker["constants"]["radius"], s=0.0
            )
            EMAN2_cppwrap.Util.add_img(tavg, temp)
            sig = EMAN2_cppwrap.Util.rotavg_fourier(temp)
            # sig = rops(pad(((cyclic_shift( projdata[i], int(sx), int(round(sy)) ) - st[0])/st[1]), mx,mx,1,0.0))
            # sig = rops(pad(((cyclic_shift(projdata, int(round(params[i][-2])), int(round(params[i][-1])) ) - st[0])/st[1])*invg, mx,mx,1,0.0))
            for k in range(nv):
                tsd.set_value_at(k, indx, tsd.get_value_at(k, indx) + sig[k])
            """Multiline Comment3"""
            tocp[indx] += 1

        ####for lll in range(len(Blockdata["accumulatepw"])):  sxprint(myid,ndata,lll,len(Blockdata["accumulatepw"][lll]))
        sp_utilities.reduce_EMData_to_root(tsd, myid, Blockdata["main_node"], mpi_comm)
        sp_utilities.reduce_EMData_to_root(tocp, myid, Blockdata["main_node"], mpi_comm)
        sp_utilities.reduce_EMData_to_root(tavg, myid, Blockdata["main_node"], mpi_comm)
        if myid == Blockdata["main_node"]:
            EMAN2_cppwrap.Util.mul_scalar(
                tavg, old_div(1.0, float(sum(Tracker["nima_per_chunk"])))
            )
            sig = EMAN2_cppwrap.Util.rotavg_fourier(tavg)
            # for k in range(1,nv):  sxprint("  BACKG  ",k,tsd.get_value_at(k,0)/tocp[0] ,sig[k],tsd.get_value_at(k,0)/tocp[0] - sig[k])
            tmp1 = [0.0] * nv
            tmp2 = [0.0] * nv
            for i in range(ngroups):
                for k in range(1, nv):
                    qt = old_div(tsd.get_value_at(k, i), tocp[i]) - sig[k]
                    if qt > 0.0:
                        tmp1[k] = old_div(2.0, qt)
                #  smooth
                tmp1[0] = tmp1[1]
                tmp1[-1] = tmp1[-2]
                for ism in range(0):  # 2
                    for k in range(1, nv - 1):
                        tmp2[k] = old_div((tmp1[k - 1] + tmp1[k] + tmp1[k + 1]), 3.0)
                    for k in range(1, nv - 1):
                        tmp1[k] = tmp2[k]
                """Multiline Comment4"""
                #  We will keep 0-element the same as first tsd.set_value_at(0,i,1.0)
                for k in range(1, nv):
                    tsd.set_value_at(k, i, tmp1[k])
                tsd.set_value_at(0, i, 1.0)
            tsd.write_image(os.path.join(Tracker["directory"], "bckgnoise.hdf"))
        sp_utilities.bcast_EMData_to_all(tsd, myid, source_node=0, comm=mpi_comm)
    nnx = tsd.get_xsize()
    nny = tsd.get_ysize()
    Blockdata["bckgnoise"] = []
    for i in range(nny):
        prj = sp_utilities.model_blank(nnx)
        for k in range(nnx):
            prj[k] = tsd.get_value_at(k, i)
        Blockdata["bckgnoise"].append(prj)  # 1.0/sigma^2
    return
    # return Blockdata["bckgnoise"]#tsd, sd#, [int(tocp[i]) for i in range(len(sd))]


def getindexdata(
    partids,
    partstack,
    particle_groups,
    original_data=None,
    small_memory=True,
    nproc=-1,
    myid=-1,
    mpi_comm=-1,
):
    global Tracker, Blockdata
    # The function will read from stack a subset of images specified in partids
    #   and assign to them parameters from partstack
    # So, the lengths of partids and partstack are the same.
    #  The read data is properly distributed among MPI threads.
    if mpi_comm < 0:
        mpi_comm = mpi.MPI_COMM_WORLD
    #  parameters
    if myid == 0:
        partstack = sp_utilities.read_text_row(partstack)
    else:
        partstack = 0
    partstack = sp_utilities.wrap_mpi_bcast(partstack, 0, mpi_comm)
    #  particles IDs
    if myid == 0:
        partids = sp_utilities.read_text_file(partids)
    else:
        partids = 0
    partids = sp_utilities.wrap_mpi_bcast(partids, 0, mpi_comm)
    #  Group assignments
    if myid == 0:
        group_reference = sp_utilities.read_text_file(particle_groups)
    else:
        group_reference = 0
    group_reference = sp_utilities.wrap_mpi_bcast(group_reference, 0, mpi_comm)

    im_start, im_end = sp_applications.MPI_start_end(len(partstack), nproc, myid)
    partstack = partstack[im_start:im_end]
    partids = partids[im_start:im_end]
    group_reference = group_reference[im_start:im_end]
    """Multiline Comment5"""

    """Multiline Comment6"""

    if original_data == None or small_memory:
        original_data = EMAN2_cppwrap.EMData.read_images(
            Tracker["constants"]["stack"], partids
        )
        for im in range(len(original_data)):
            original_data[im].set_attr("particle_group", group_reference[im])
    return original_data, partstack, [im_start, im_end]


def get_shrink_data(
    nxinit,
    procid,
    original_data=None,
    oldparams=None,
    return_real=False,
    preshift=False,
    apply_mask=True,
    nonorm=False,
    nosmearing=False,
    npad=1,
):
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
    # from fundamentals import resample

    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_global_def.sxprint("  ")
        line = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime()) + " =>"
        sp_global_def.sxprint(
            line,
            "Processing data  onx: %3d, nx: %3d, CTF: %s, applymask: %s, preshift: %s."
            % (
                Tracker["constants"]["nnxo"],
                nxinit,
                Tracker["constants"]["CTF"],
                apply_mask,
                preshift,
            ),
        )
    #  Preprocess the data
    mask2D = sp_utilities.model_circle(
        Tracker["constants"]["radius"],
        Tracker["constants"]["nnxo"],
        Tracker["constants"]["nnxo"],
    )
    nima = len(original_data)
    shrinkage = old_div(nxinit, float(Tracker["constants"]["nnxo"]))

    #  Note these are in Fortran notation for polar searches
    # txm = float(nxinit-(nxinit//2+1) - radius -1)
    # txl = float(2 + radius - nxinit//2+1)
    radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
    txm = float(nxinit - (old_div(nxinit, 2) + 1) - radius)
    txl = float(radius - old_div(nxinit, 2) + 1)

    if Blockdata["bckgnoise"]:
        oneover = []
        nnx = Blockdata["bckgnoise"][0].get_xsize()
        for i in range(len(Blockdata["bckgnoise"])):
            temp = [0.0] * nnx
            for k in range(nnx):
                if Blockdata["bckgnoise"][i].get_value_at(k) > 0.0:
                    temp[k] = old_div(
                        1.0, numpy.sqrt(Blockdata["bckgnoise"][i].get_value_at(k))
                    )
            oneover.append(temp)
        del temp

    Blockdata["accumulatepw"][procid] = [None] * nima
    data = [None] * nima
    for im in range(nima):

        if Tracker["mainiteration"] <= 1:
            phi, theta, psi, sx, sy, = (
                oldparams[im][0],
                oldparams[im][1],
                oldparams[im][2],
                oldparams[im][3],
                oldparams[im][4],
            )
            wnorm = 1.0
        else:
            phi, theta, psi, sx, sy, wnorm = (
                oldparams[im][0],
                oldparams[im][1],
                oldparams[im][2],
                oldparams[im][3],
                oldparams[im][4],
                oldparams[im][7],
            )

        """Multiline Comment7"""

        if preshift:
            # data[im] = fshift(original_data[im], sx, sy)
            if nosmearing:
                data[im] = sp_fundamentals.fshift(original_data[im], sx, sy)
                oldparams[im][3] = 0.0
                oldparams[im][4] = 0.0
            else:
                sx = int(round(sx))
                sy = int(round(sy))
                data[im] = sp_fundamentals.cyclic_shift(original_data[im], sx, sy)
                #  Put rounded shifts on the list, note image has the original floats - check whether it may cause problems
                oldparams[im][3] = sx
                oldparams[im][4] = sy
            sx = 0.0
            sy = 0.0
        else:
            data[im] = original_data[im].copy()

        st = get_image_statistics(data[im], mask2D, False)
        data[im] -= st[0]
        data[im] = old_div(data[im], st[1])
        if data[im].get_attr_default("bckgnoise", None):
            data[im].delete_attr("bckgnoise")
        #  Do bckgnoise if exists
        if Blockdata["bckgnoise"]:
            if apply_mask:
                if Tracker["constants"]["hardmask"]:
                    data[im] = sp_morphology.cosinemask(
                        data[im], radius=Tracker["constants"]["radius"]
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
                        sp_filter.filt_table(
                            bckg, oneover[data[im].get_attr("particle_group")]
                        )
                    )
                    #  Normalize bckg noise in real space, only region actually used.
                    st = get_image_statistics(bckg, mask2D, False)
                    bckg -= st[0]
                    bckg = old_div(bckg, st[1])
                    data[im] = sp_morphology.cosinemask(
                        data[im], radius=Tracker["constants"]["radius"], bckg=bckg
                    )
        else:
            #  if no bckgnoise, do simple masking instead
            if apply_mask:
                data[im] = sp_morphology.cosinemask(
                    data[im], radius=Tracker["constants"]["radius"]
                )
        #  resample will properly adjusts shifts and pixel size in ctf
        # data[im] = resample(data[im], shrinkage)
        #  return Fourier image
        # if npad> 1:  data[im] = pad(data[im], Tracker["constants"]["nnxo"]*npad, Tracker["constants"]["nnxo"]*npad, 1, 0.0)

        #  Apply varadj
        if not nonorm:
            EMAN2_cppwrap.Util.mul_scalar(
                data[im], old_div(Tracker["avgvaradj"][procid], wnorm)
            )
            # print(Tracker["avgvaradj"][procid]/wnorm)

        #  FT
        data[im] = sp_fundamentals.fft(data[im])
        sig = EMAN2_cppwrap.Util.rotavg_fourier(data[im])
        Blockdata["accumulatepw"][procid][im] = sig[old_div(len(sig), 2) :] + [0.0]

        if Tracker["constants"]["CTF"]:
            data[im] = sp_fundamentals.fdecimate(
                data[im], nxinit * npad, nxinit * npad, 1, False, False
            )
            ctf_params = original_data[im].get_attr("ctf")
            ctf_params.apix = old_div(ctf_params.apix, shrinkage)
            data[im].set_attr("ctf", ctf_params)
            # if Tracker["applyctf"] :  #  This should be always False
            # 	data[im] = filt_ctf(data[im], ctf_params, dopad=False)
            # 	data[im].set_attr('ctf_applied', 1)
            # else:
            data[im].set_attr("ctf_applied", 0)
            if return_real:
                data[im] = sp_fundamentals.fft(data[im])
        else:
            ctf_params = original_data[im].get_attr_default("ctf", False)
            if ctf_params:
                ctf_params.apix = old_div(ctf_params.apix, shrinkage)
                data[im].set_attr("ctf", ctf_params)
                data[im].set_attr("ctf_applied", 0)
            data[im] = sp_fundamentals.fdecimate(
                data[im], nxinit * npad, nxinit * npad, 1, True, False
            )
            apix = Tracker["constants"]["pixel_size"]
            data[im].set_attr("apix", old_div(apix, shrinkage))

        #  We have to make sure the shifts are within correct range, shrinkage or not
        sp_utilities.set_params_proj(
            data[im],
            [
                phi,
                theta,
                psi,
                max(min(sx * shrinkage, txm), txl),
                max(min(sy * shrinkage, txm), txl),
            ],
        )
        if not return_real:
            data[im].set_attr("padffted", 1)
        data[im].set_attr("npad", npad)
        if Blockdata["bckgnoise"]:
            temp = Blockdata["bckgnoise"][data[im].get_attr("particle_group")]
            ###  Do not adjust the values, we try to keep everything in the same Fourier values.
            data[im].set_attr("bckgnoise", [temp[i] for i in range(temp.get_xsize())])
    return data


def get_anger(angle1, angle2):
    A1 = sp_fundamentals.rotmatrix(angle1[0], angle1[1], angle1[2])
    ar = Blockdata["symclass"].symmetry_related(angle2)
    axes_dis_min = 1.0e23
    for q in ar:
        A2 = sp_fundamentals.rotmatrix(q[0], q[1], q[2])
        axes_dis = 0.0
        for i in range(3):
            axes_dis += sp_utilities.lacos(
                A1[i][0] * A2[i][0] + A1[i][1] * A2[i][1] + A1[i][2] * A2[i][2]
            )
        axes_dis_min = min(axes_dis_min, old_div(axes_dis, 3.0))
    return axes_dis_min


def checkstep(item, keepchecking):
    global Tracker, Blockdata
    if Blockdata["myid"] == Blockdata["main_node"]:
        if keepchecking:
            if os.path.exists(item):
                doit = 0
            else:
                doit = 1
                keepchecking = False
        else:
            doit = 1
    else:
        doit = 1
    doit = sp_utilities.bcast_number_to_all(doit, source_node=Blockdata["main_node"])
    return doit, keepchecking


def out_fsc(f):
    global Tracker, Blockdata
    sp_global_def.sxprint(" ")
    sp_global_def.sxprint(
        "  driver FSC  after  iteration#%3d" % Tracker["mainiteration"]
    )
    sp_global_def.sxprint("  %4d        %7.2f         %7.3f" % (0, 1000.00, f[0]))
    for i in range(1, len(f)):
        sp_global_def.sxprint(
            "  %4d        %7.2f         %7.3f"
            % (
                i,
                old_div(
                    Tracker["constants"]["pixel_size"] * Tracker["constants"]["nnxo"],
                    float(i),
                ),
                f[i],
            )
        )
    sp_global_def.sxprint(" ")


def get_refangs_and_shifts():
    global Tracker, Blockdata

    refang = Blockdata["symclass"].even_angles(Tracker["delta"])
    coarse = Blockdata["symclass"].even_angles(
        2 * Tracker["delta"],
        theta1=Tracker["theta_min"],
        theta2=Tracker["theta_max"],
        method=Tracker["constants"]["angle_method"],
    )
    refang = Blockdata["symclass"].reduce_anglesets(
        sp_fundamentals.rotate_params(
            refang,
            [-0.5 * Tracker["delta"], -0.5 * Tracker["delta"], -0.5 * Tracker["delta"]],
        )
    )

    """Multiline Comment8"""

    k = int(numpy.ceil(old_div(Tracker["xr"], Tracker["ts"])))
    radi = (Tracker["xr"] + Tracker["ts"]) ** 2
    rshifts = []
    for ix in range(-k, k, 1):
        six = ix * Tracker["ts"] + old_div(Tracker["ts"], 2)
        for iy in range(-k, k, 1):
            siy = iy * Tracker["ts"] + old_div(Tracker["ts"], 2)
            if six * six + siy * siy <= radi:
                rshifts.append([six, siy])

    ts_coarse = 2 * Tracker["ts"]
    k = int(numpy.ceil(old_div(Tracker["xr"], ts_coarse)))
    radi = Tracker["xr"] * Tracker["xr"]
    coarse_shifts = []
    for ix in range(-k, k + 1, 1):
        six = ix * ts_coarse
        for iy in range(-k, k + 1, 1):
            siy = iy * ts_coarse
            if six * six + siy * siy <= radi:
                coarse_shifts.append([six, siy])

    return refang, rshifts, coarse, coarse_shifts


def get_shifts_neighbors(rshifts, cs):
    """
	  	output is a list of shift neighbors to cs within ts*1.5 distance\
	  	This yields at most five shifts
	"""
    shiftneighbors = []
    rad = Tracker["ts"] * 1.5
    for l, q in enumerate(rshifts):
        if sp_utilities.get_dist(q, cs) <= rad:
            shiftneighbors.append(l)
    return shiftneighbors


def shakegrid(rshifts, qt):
    for i in range(len(rshifts)):
        rshifts[i][0] += qt
        rshifts[i][1] += qt


###----------------


def get_refvol(nxinit):
    ref_vol = sp_utilities.get_im(Tracker["refvol"])
    nnn = ref_vol.get_xsize()
    if nxinit != nnn:
        ref_vol = sp_fundamentals.fdecimate(
            ref_vol, nxinit, nxinit, nxinit, True, False
        )
    return ref_vol


def do3d(
    procid,
    data,
    newparams,
    refang,
    rshifts,
    norm_per_particle,
    myid,
    smearing=True,
    mpi_comm=-1,
):
    global Tracker, Blockdata

    #  Without filtration
    if mpi_comm == -1:
        mpi_comm = mpi.MPI_COMM_WORLD

    if procid == 0:
        if Blockdata["no_of_groups"] > 1:
            if myid == Blockdata["nodes"][0]:
                if os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
                    sp_global_def.sxprint("tempdir exists")
                else:
                    try:
                        os.makedirs(os.path.join(Tracker["directory"], "tempdir"))
                    except:
                        sp_global_def.sxprint("tempdir exists")
        else:
            if myid == Blockdata["main_node"]:
                if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
                    try:
                        os.makedirs(os.path.join(Tracker["directory"], "tempdir"))
                    except:
                        sp_global_def.sxprint("tempdir exists")
                else:
                    sp_global_def.sxprint("tempdir exists")
    mpi.mpi_barrier(mpi_comm)
    """Multiline Comment10"""
    shrinkage = old_div(float(Tracker["nxinit"]), float(Tracker["constants"]["nnxo"]))
    if smearing:
        tvol, tweight, trol = sp_reconstruction.recons3d_trl_struct_MPI(
            myid=myid,
            main_node=Blockdata["nodes"][procid],
            prjlist=data,
            paramstructure=newparams,
            refang=refang,
            rshifts_shrank=[[q[0] * shrinkage, q[1] * shrinkage] for q in rshifts],
            delta=Tracker["delta"],
            CTF=Tracker["constants"]["CTF"],
            upweighted=False,
            mpi_comm=mpi_comm,
            target_size=(2 * Tracker["nxinit"] + 3),
            avgnorm=Tracker["avgvaradj"][procid],
            norm_per_particle=norm_per_particle,
        )
    else:
        tvol, tweight, trol = recons3d_trl_struct_MPI_nosmearing(
            myid=myid,
            main_node=Blockdata["nodes"][procid],
            prjlist=data,
            parameters=newparams,
            CTF=Tracker["constants"]["CTF"],
            upweighted=False,
            mpi_comm=mpi_comm,
            target_size=(2 * Tracker["nxinit"] + 3),
        )

    if Blockdata["no_of_groups"] > 1:
        if myid == Blockdata["nodes"][procid]:
            while not os.path.exists(os.path.join(Tracker["directory"], "tempdir")):
                time.sleep(5)
            tvol.set_attr("is_complex", 0)
            tvol.write_image(
                os.path.join(
                    Tracker["directory"],
                    "tempdir",
                    "tvol_%01d_%03d.hdf" % (procid, Tracker["mainiteration"]),
                )
            )
            tweight.write_image(
                os.path.join(
                    Tracker["directory"],
                    "tempdir",
                    "tweight_%01d_%03d.hdf" % (procid, Tracker["mainiteration"]),
                )
            )
            trol.write_image(
                os.path.join(
                    Tracker["directory"],
                    "tempdir",
                    "trol_%01d_%03d.hdf" % (procid, Tracker["mainiteration"]),
                )
            )
            line = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime()) + " =>"
            sp_global_def.sxprint(
                line, "Executed successfully backprojection for group ", procid
            )
    else:
        if myid == Blockdata["main_node"]:
            tvol.set_attr("is_complex", 0)
            tvol.write_image(
                os.path.join(
                    Tracker["directory"],
                    "tempdir",
                    "tvol_%01d_%03d.hdf" % (procid, Tracker["mainiteration"]),
                )
            )
            tweight.write_image(
                os.path.join(
                    Tracker["directory"],
                    "tempdir",
                    "tweight_%01d_%03d.hdf" % (procid, Tracker["mainiteration"]),
                )
            )
            trol.write_image(
                os.path.join(
                    Tracker["directory"],
                    "tempdir",
                    "trol_%01d_%03d.hdf" % (procid, Tracker["mainiteration"]),
                )
            )
            line = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime()) + " =>"
            sp_global_def.sxprint(
                line, "Executed successfully backprojection for group ", procid
            )
    mpi.mpi_barrier(mpi_comm)
    return


def print_dict(dict, theme):
    line = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime()) + " =>"
    sp_global_def.sxprint(line, theme)
    spaces = "                    "
    exclude = [
        "constants",
        "maxit",
        "nodes",
        "yr",
        "shared_comm",
        "bckgnoise",
        "myid",
        "myid_on_node",
        "accumulatepw",
        "changed_delta",
    ]
    for key, value in sorted(dict.items()):
        pt = True
        for ll in exclude:
            if key == ll:
                pt = False
                break
        if pt:
            sp_global_def.sxprint(
                "                    => ", key + spaces[len(key) :], ":  ", value
            )


def stepone(tvol, tweight):
    global Tracker, Blockdata
    tvol.set_attr("is_complex", 1)
    ovol = EMAN2_cppwrap.Util.shrinkfvol(tvol, 2)
    owol = EMAN2_cppwrap.Util.shrinkfvol(tweight, 2)
    if Tracker["constants"]["symmetry"] != "c1":
        ovol = ovol.symfvol(Tracker["constants"]["symmetry"], -1)
        owol = owol.symfvol(Tracker["constants"]["symmetry"], -1)
    # print(info(ovol,Comment = " shrank ovol"))
    return EMAN2_cppwrap.Util.divn_cbyr(ovol, owol)


def steptwo_mpi(tvol, tweight, treg, cfsc=None, regularized=True, color=0):
    global Tracker, Blockdata

    n_iter = 10  # number of iterations for iterative process for doing map
    if Blockdata["color"] != color:
        return sp_utilities.model_blank(
            1
        )  # This should not be executed if called properly
    if Blockdata["myid_on_node"] == 0:
        nz = tweight.get_zsize()
        ny = tweight.get_ysize()
        nx = tweight.get_xsize()
        tvol.set_attr("is_complex", 1)
        if regularized:
            nr = len(cfsc)
            ovol = [0.0] * nr
            limitres = 0
            for i in range(nr):
                ovol[i] = min(max(cfsc[i], 0.0), 0.999)
                # print( i,cfsc[i] )
                if ovol[i] == 0.0:
                    limitres = i - 1
                    break
            if limitres == 0:
                limitres = nr - 2
            ovol = sp_utilities.reshape_1d(ovol, nr, 2 * nr)
            limitres = 2 * min(
                limitres, Tracker["maxfrad"]
            )  # 2 on account of padding, which is always on
            maxr2 = limitres ** 2
            for i in range(limitres + 1, len(ovol), 1):
                ovol[i] = 0.0
            ovol[0] = 1.0
            # print(" ovol  ", ovol)
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
    # Util.iterefa(tvol, tweight, maxr2, Tracker["constants"]["nnxo"])

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
        tvol.set_attr("npad", 2)
        tvol.div_sinc(1)
        tvol.del_attr("npad")
        tvol = EMAN2_cppwrap.Util.window(
            tvol,
            Tracker["constants"]["nnxo"],
            Tracker["constants"]["nnxo"],
            Tracker["constants"]["nnxo"],
        )
        tvol = sp_morphology.cosinemask(
            tvol, old_div(Tracker["constants"]["nnxo"], 2) - 1 - 5, 5, None
        )  # clean artifacts in corners
        return tvol
    else:
        return None


def calculate_2d_params_for_centering(kwargs):
    #################################################################################################################################################################
    # get parameters from the dictionary
    init2dir = kwargs["init2dir"]
    myid = kwargs["myid"]
    main_node = kwargs["main_node"]
    number_of_images_in_stack = kwargs["number_of_images_in_stack"]
    nproc = kwargs["nproc"]

    target_radius = kwargs["target_radius"]
    # target_nx = kwargs["target_nx"]
    radi = kwargs["radi"]

    center_method = kwargs["center_method"]

    nxrsteps = kwargs["nxrsteps"]

    # stack_processed_by_ali2d_base__filename = kwargs["stack_processed_by_ali2d_base__filename"]
    command_line_provided_stack_filename = kwargs[
        "command_line_provided_stack_filename"
    ]

    # masterdir = kwargs["masterdir"]

    options_skip_prealignment = kwargs["options_skip_prealignment"]
    options_CTF = kwargs["options_CTF"]

    mpi_comm = kwargs["mpi_comm"]
    #################################################################################################################################################################

    if options_skip_prealignment:
        if Blockdata["myid"] == 0:
            sp_global_def.sxprint("=========================================")
            sp_global_def.sxprint(" >>> There is no pre-alignment step.")
            sp_global_def.sxprint("=========================================")
            return [[0, 0, 0, 0, 0] for i in range(number_of_images_in_stack)]
        else:
            return [0.0]
    Finished_initial_2d_alignment = 1
    if Blockdata["myid"] == Blockdata["main_node"]:
        if os.path.exists(os.path.join(init2dir, "Finished_initial_2d_alignment.txt")):
            Finished_initial_2d_alignment = 0
    Finished_initial_2d_alignment = sp_utilities.bcast_number_to_all(
        Finished_initial_2d_alignment, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    if Finished_initial_2d_alignment == 1:

        if Blockdata["myid"] == 0:
            #  Create output directory
            log2d = sp_logger.Logger(sp_logger.BaseLogger_Files())
            log2d.prefix = os.path.join(init2dir)
            cmd = "mkdir -p " + log2d.prefix
            outcome = subprocess.call(cmd, shell=True)
            log2d.prefix += "/"
            # outcome = subprocess.call("sxheader.py  "+command_line_provided_stack_filename+"   --params=xform.align2d  --zero", shell=True)
        else:
            outcome = 0
            log2d = None

        if Blockdata["myid"] == Blockdata["main_node"]:
            a = sp_utilities.get_im(command_line_provided_stack_filename)
            nnxo = a.get_xsize()
        else:
            nnxo = 0
        nnxo = sp_utilities.bcast_number_to_all(
            nnxo, source_node=Blockdata["main_node"]
        )

        image_start, image_end = sp_applications.MPI_start_end(
            number_of_images_in_stack, Blockdata["nproc"], Blockdata["myid"]
        )

        original_images = EMAN2_cppwrap.EMData.read_images(
            command_line_provided_stack_filename, list(range(image_start, image_end))
        )
        #  We assume the target radius will be 29, and xr = 1.
        shrink_ratio = old_div(float(target_radius), float(radi))

        for im in range(len(original_images)):
            if shrink_ratio != 1.0:
                original_images[im] = sp_fundamentals.resample(
                    original_images[im], shrink_ratio
                )

        nx = original_images[0].get_xsize()
        # nx = int(nx*shrink_ratio + 0.5)

        txrm = old_div((nx - 2 * (target_radius + 1)), 2)
        if txrm < 0:
            sp_global_def.ERROR(
                "ERROR!!   Radius of the structure larger than the window data size permits   %d"
                % (radi),
                "2d prealignment",
                1,
                Blockdata["myid"],
            )
        if old_div(txrm, nxrsteps) > 0:
            tss = ""
            txr = ""
            while old_div(txrm, nxrsteps) > 0:
                tts = old_div(txrm, nxrsteps)
                tss += "  %d" % tts
                txr += "  %d" % (tts * nxrsteps)
                txrm = old_div(txrm, 2)
        else:
            tss = "1"
            txr = "%d" % txrm

        # section ali2d_base

        params2d = sp_applications.ali2d_base(
            original_images,
            init2dir,
            None,
            1,
            target_radius,
            1,
            txr,
            txr,
            tss,
            False,
            90.0,
            center_method,
            14,
            options_CTF,
            1.0,
            False,
            "ref_ali2d",
            "",
            log2d,
            Blockdata["nproc"],
            Blockdata["myid"],
            Blockdata["main_node"],
            mpi_comm,
            write_headers=False,
        )

        del original_images

        for i in range(len(params2d)):
            alpha, sx, sy, mirror = sp_utilities.combine_params2(
                0, params2d[i][1], params2d[i][2], 0, -params2d[i][0], 0, 0, 0
            )
            sx = old_div(sx, shrink_ratio)
            sy = old_div(sy, shrink_ratio)
            params2d[i][0] = 0.0
            params2d[i][1] = sx
            params2d[i][2] = sy
            params2d[i][3] = 0

        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

        params2d = sp_utilities.wrap_mpi_gatherv(
            params2d, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        if Blockdata["myid"] == Blockdata["main_node"]:
            sp_utilities.write_text_row(
                params2d, os.path.join(init2dir, "initial2Dparams.txt")
            )
            return params2d
        else:
            return [0.0]
    else:
        if Blockdata["myid"] == Blockdata["main_node"]:
            params2d = sp_utilities.read_text_row(
                os.path.join(init2dir, "initial2Dparams.txt")
            )
            sp_global_def.sxprint("Skipping 2d alignment since it was already done!")
            return params2d
        else:
            return [0.0]


"""Multiline Comment11"""


def Numrinit_local(first_ring, last_ring, skip=1, mode="F"):
    """This function calculates the necessary information for the 2D
	   polar interpolation. For each ring, three elements are recorded:
	   numr[i*3]:  Radius of this ring
	   numr[i*3+1]: Total number of samples of all inner rings+1
	   		(Or, the beginning point of this ring)
	   numr[i*3+2]: Number of samples of this ring. This number is an
	   		FFT-friendly power of the 2.

	   "F" means a full circle interpolation
	   "H" means a half circle interpolation
	"""
    MAXFFT = 32768
    skip = 1
    if mode == "f" or mode == "F":
        dpi = 2 * numpy.pi
    else:
        dpi = numpy.pi
    numr = []
    lcirc = 1
    for k in range(first_ring, last_ring + 1, skip):
        numr.append(k)
        jp = int(dpi * 1.5 * k)
        ip = 2 ** (sp_alignment.log2(jp))  # two times oversample each ring
        ip = min(MAXFFT, ip)
        # if ip >16 : ip/=2  #  subsample rings
        # if (k+skip <= last_ring and jp > ip+ip//2): ip=min(MAXFFT,2*ip)
        # if (k+skip  > last_ring and jp > ip+ip//5): ip=min(MAXFFT,2*ip)

        numr.append(lcirc)
        numr.append(ip)
        lcirc += ip

    return numr


"""Multiline Comment12"""


def ali3D_polar_ccc(
    refang,
    shifts,
    coarse_angles,
    coarse_shifts,
    procid,
    original_data=None,
    oldparams=None,
    preshift=False,
    apply_mask=True,
    nonorm=False,
    applyctf=True,
):
    global Tracker, Blockdata
    #  Input data has to be CTF-multiplied, preshifted
    #  Output - newpar, see structure
    #    newpar = [[i, [worst_similarity, sum_all_similarities], [[-1, -1.0e23] for j in range(Tracker["lentop"])]] for i in range(len(data))]
    #    newpar = [[params],[],... len(data)]
    #    params = [particleID, [worst_similarity, sum_all_similarities],[imageallparams]]]
    #    imageallparams = [[orientation, similarity],[],...  number of all orientations ]
    #  Coding of orientations:
    #    hash = ang*100000000 + lpsi*1000 + ishift
    #    ishift = hash%1000
    #    ipsi = (hash/1000)%100000
    #    iang  = hash/100000000
    #  To get best matching for particle #kl:
    #     hash_best = newpar[kl][-1][0][0]
    #     best_sim  = newpar[kl][-1][0][1]
    #  To sort:
    # from operator 		import itemgetter#, attrgetter, methodcaller
    #   params.sort(key=itemgetter(2))
    # from fundamentals import resample

    if Blockdata["myid"] == Blockdata["main_node"]:
        print_dict(Tracker, "PROJECTION MATCHING parameters of polar CCC")

    at = time.time()
    shrinkage = old_div(float(Tracker["nxpolar"]), float(Tracker["constants"]["nnxo"]))
    shrink = old_div(float(Tracker["nxinit"]), float(Tracker["constants"]["nnxo"]))
    radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
    mode = "F"
    numr = Numrinit_local(1, radius, 1, mode)
    wr = sp_alignment.ringwe(numr, mode)
    cnx = float(old_div(Tracker["nxpolar"], 2) + 1)

    ##if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  RADIUS IN POLAR  ",radius,numr[-1])
    #  FINE SEARCH CONSTANTS
    #  Fine grids are shifted by half-fine_step
    nang = len(refang)
    npsi = int(old_div(360.0, Tracker["delta"]))
    nshifts = len(shifts)
    n_fine_shifts = 4

    nima = len(original_data)
    mask = EMAN2_cppwrap.Util.unrollmask(Tracker["nxinit"], Tracker["nxinit"])
    for j in range(old_div(Tracker["nxinit"], 2), Tracker["nxinit"]):
        mask[0, j] = 1.0
    mask2D = sp_utilities.model_circle(
        Tracker["constants"]["radius"],
        Tracker["constants"]["nnxo"],
        Tracker["constants"]["nnxo"],
    )
    reachpw = (
        mask.get_xsize()
    )  # The last element of accumulated pw is zero so for the full size nothing is added.

    #  COARSE SEARCH CONSTANTS
    n_coarse_ang = len(coarse_angles)
    coarse_delta = 2 * Tracker["delta"]
    n_coarse_psi = int(old_div(360.0, coarse_delta))
    n_coarse_shifts = len(coarse_shifts)

    coarse_shifts_shrank = [None] * n_coarse_shifts
    for ib in range(n_coarse_shifts):
        coarse_shifts_shrank[ib] = [
            coarse_shifts[ib][0] * shrinkage,
            coarse_shifts[ib][1] * shrinkage,
        ]

    ###if(Blockdata["myid"] == Blockdata["main_node"]): sxprint("   TRETR  ",Tracker["constants"]["nnxo"],Tracker["nxinit"],reachpw)

    disp_unit = sp_helix_sphire.np.dtype("f4").itemsize

    #  REFVOL
    if Blockdata["myid_on_node"] == 0:
        if Blockdata["myid"] == Blockdata["main_node"]:
            odo = get_refvol(Tracker["nxpolar"])
            nxvol = odo.get_xsize()
            nyvol = odo.get_ysize()
            nzvol = odo.get_zsize()
        else:
            nxvol = 0
            nyvol = 0
            nzvol = 0

        nxvol = sp_utilities.bcast_number_to_all(
            nxvol,
            source_node=Blockdata["main_node"],
            mpi_comm=Blockdata["group_zero_comm"],
        )
        nyvol = sp_utilities.bcast_number_to_all(
            nyvol,
            source_node=Blockdata["main_node"],
            mpi_comm=Blockdata["group_zero_comm"],
        )
        nzvol = sp_utilities.bcast_number_to_all(
            nzvol,
            source_node=Blockdata["main_node"],
            mpi_comm=Blockdata["group_zero_comm"],
        )

        if Blockdata["myid"] != Blockdata["main_node"]:
            odo = sp_utilities.model_blank(nxvol, nyvol, nzvol)

        sp_utilities.bcast_EMData_to_all(
            odo,
            Blockdata["group_zero_myid"],
            source_node=Blockdata["main_node"],
            comm=Blockdata["group_zero_comm"],
        )

        odo = sp_projection.prep_vol(odo, npad=2, interpolation_method=1)
        ndo = EMAN2_cppwrap.EMNumPy.em2numpy(odo)
        nxvol = odo.get_xsize()
        nyvol = odo.get_ysize()
        nzvol = odo.get_zsize()
        orgsizevol = nxvol * nyvol * nzvol
        sizevol = orgsizevol
    else:
        orgsizevol = 0
        sizevol = 0
        nxvol = 0
        nyvol = 0
        nzvol = 0

    orgsizevol = sp_utilities.bcast_number_to_all(
        orgsizevol, source_node=Blockdata["main_node"]
    )
    nxvol = sp_utilities.bcast_number_to_all(nxvol, source_node=Blockdata["main_node"])
    nyvol = sp_utilities.bcast_number_to_all(nyvol, source_node=Blockdata["main_node"])
    nzvol = sp_utilities.bcast_number_to_all(nzvol, source_node=Blockdata["main_node"])

    win_vol, base_vol = mpi.mpi_win_allocate_shared(
        sizevol * disp_unit, disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"]
    )
    sizevol = orgsizevol
    if Blockdata["myid_on_node"] != 0:
        base_vol, = mpi.mpi_win_shared_query(win_vol, mpi.MPI_PROC_NULL)

    """
    Comments from Adnan:
    numpy.core.mutliarray.int_asbuffer function is not available anymore . That is why a different alternative is tried here.
    """
    # volbuf = sp_helix_sphire.np.frombuffer(
    #     sp_helix_sphire.np.core.multiarray.int_asbuffer(base_vol, sizevol * disp_unit),
    #     dtype="f4",
    # )
    ptr_n = ctypes.cast(base_vol, ctypes.POINTER(ctypes.c_int * sizevol))
    volbuf = numpy.frombuffer(ptr_n.contents, dtype="f4")
    volbuf = volbuf.reshape(nzvol, nyvol, nxvol)
    if Blockdata["myid_on_node"] == 0:
        sp_helix_sphire.np.copyto(volbuf, ndo)
        del odo, ndo

    # volprep = EMNumPy.assign_numpy_to_emdata(volbuf)
    emnumpy1 = EMAN2_cppwrap.EMNumPy()
    volprep = emnumpy1.register_numpy_to_emdata(volbuf)
    volprep.set_attr_dict(
        {
            "is_complex": 1,
            "is_complex_x": 0,
            "is_fftodd": 0,
            "is_fftpad": 1,
            "is_shuffled": 1,
            "npad": 2,
        }
    )
    crefim = EMAN2_cppwrap.Util.Polar2Dm(
        sp_utilities.model_blank(Tracker["nxpolar"], Tracker["nxpolar"]),
        cnx,
        cnx,
        numr,
        mode,
    )

    #  BIG BUFFER (for n_coarse_ang polar arrays)
    size_of_one_image = crefim.get_xsize()
    lenbigbuf = n_coarse_ang
    orgsize = (
        lenbigbuf * size_of_one_image
    )  # This is number of projections to be computed simultaneously times their size

    if Blockdata["myid_on_node"] == 0:
        size = orgsize
    else:
        size = 0

    win_sm, base_ptr = mpi.mpi_win_allocate_shared(
        size * disp_unit, disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"]
    )
    size = orgsize
    if Blockdata["myid_on_node"] != 0:
        base_ptr, = mpi.mpi_win_shared_query(win_sm, mpi.MPI_PROC_NULL)

    """
    Comments from Adnan:
    numpy.core.mutliarray.int_asbuffer function is not available anymore . That is why a different alternative is tried here.
    """
    # buffer = sp_helix_sphire.np.frombuffer(
    #     sp_helix_sphire.np.core.multiarray.int_asbuffer(base_ptr, size * disp_unit),
    #     dtype="f4",
    # )
    ptr_n = ctypes.cast(base_ptr, ctypes.POINTER(ctypes.c_int * size))
    buffer = numpy.frombuffer(ptr_n.contents, dtype="f4")

    buffer = buffer.reshape(lenbigbuf, size_of_one_image)
    # bigbuffer = EMNumPy.assign_numpy_to_emdata(buffer)

    emnumpy2 = EMAN2_cppwrap.EMNumPy()
    bigbuffer = emnumpy2.register_numpy_to_emdata(buffer)

    #  end of setup

    if Tracker["nxinit"] != Tracker["nxpolar"]:
        mpi.mpi_win_free(win_vol)
        #  REFVOL FOR ML
        if Blockdata["myid_on_node"] == 0:
            if Blockdata["myid"] == Blockdata["main_node"]:
                odo = get_refvol(Tracker["nxpolar"])
                nxvol = odo.get_xsize()
                nyvol = odo.get_ysize()
                nzvol = odo.get_zsize()
            else:
                nxvol = 0
                nyvol = 0
                nzvol = 0

            nxvol = sp_utilities.bcast_number_to_all(
                nxvol,
                source_node=Blockdata["main_node"],
                mpi_comm=Blockdata["group_zero_comm"],
            )
            nyvol = sp_utilities.bcast_number_to_all(
                nyvol,
                source_node=Blockdata["main_node"],
                mpi_comm=Blockdata["group_zero_comm"],
            )
            nzvol = sp_utilities.bcast_number_to_all(
                nzvol,
                source_node=Blockdata["main_node"],
                mpi_comm=Blockdata["group_zero_comm"],
            )

            if Blockdata["myid"] != Blockdata["main_node"]:
                odo = sp_utilities.model_blank(nxvol, nyvol, nzvol)

            sp_utilities.bcast_EMData_to_all(
                odo,
                Blockdata["group_zero_myid"],
                source_node=Blockdata["main_node"],
                comm=Blockdata["group_zero_comm"],
            )

            odo = sp_projection.prep_vol(odo, npad=2, interpolation_method=1)
            ndo = EMAN2_cppwrap.EMNumPy.em2numpy(odo)
            nxvol = odo.get_xsize()
            nyvol = odo.get_ysize()
            nzvol = odo.get_zsize()
            orgsizevol = nxvol * nyvol * nzvol
            sizevol = orgsizevol
        else:
            orgsizevol = 0
            sizevol = 0
            nxvol = 0
            nyvol = 0
            nzvol = 0

        orgsizevol = sp_utilities.bcast_number_to_all(
            orgsizevol, source_node=Blockdata["main_node"]
        )
        nxvol = sp_utilities.bcast_number_to_all(
            nxvol, source_node=Blockdata["main_node"]
        )
        nyvol = sp_utilities.bcast_number_to_all(
            nyvol, source_node=Blockdata["main_node"]
        )
        nzvol = sp_utilities.bcast_number_to_all(
            nzvol, source_node=Blockdata["main_node"]
        )

        win_volinit, base_volinit = mpi.mpi_win_allocate_shared(
            sizevol * disp_unit, disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"]
        )
        sizevol = orgsizevol
        if Blockdata["myid_on_node"] != 0:
            base_volinit, = mpi.mpi_win_shared_query(win_volinit, mpi.MPI_PROC_NULL)
        """
        Comments from Adnan:
        numpy.core.mutliarray.int_asbuffer function is not available anymore . That is why a different alternative is tried here.
        """
        # volbufinit = sp_helix_sphire.np.frombuffer(
        #     sp_helix_sphire.np.core.multiarray.int_asbuffer(
        #         base_volinit, sizevol * disp_unit
        #     ),
        #     dtype="f4",
        # )
        ptr = ctypes.cast(base_volinit, ctypes.POINTER(ctypes.c_int * sizevol))
        volbufinit = numpy.frombuffer(ptr.contents, dtype="f4")
        volbufinit = volbufinit.reshape(nzvol, nyvol, nxvol)
        if Blockdata["myid_on_node"] == 0:
            sp_helix_sphire.np.copyto(volbufinit, ndoinit)
            del odo, ndoinit

        # volinit = EMNumPy.assign_numpy_to_emdata(volbufinit)
        emnumpy4 = EMAN2_cppwrap.EMNumPy()
        volinit = emnumpy4.register_numpy_to_emdata(volbufinit)
        volinit.set_attr_dict(
            {
                "is_complex": 1,
                "is_complex_x": 0,
                "is_fftodd": 0,
                "is_fftpad": 1,
                "is_shuffled": 1,
                "npad": 2,
            }
        )
        if Blockdata["myid_on_node"] == 0:
            volinit.update()
        mpi.mpi_barrier(Blockdata["shared_comm"])
        ###if( Blockdata["myid"] < 10 ):
        ###	sxprint(" VOLPREP  ",volinit[0],volinit[1],info(volinit))
    else:
        volinit = volprep
    #  End of replaced volprep

    at = time.time()

    nang_start, nang_end = sp_applications.MPI_start_end(
        n_coarse_ang, Blockdata["no_of_processes_per_group"], Blockdata["myid_on_node"]
    )

    for i in range(
        nang_start, nang_end, 1
    ):  # This will take care of process on a node less than nang.  Some loops will not be executed
        temp = sp_projection.prgl(
            volprep, [coarse_angles[i][0], coarse_angles[i][1], 0.0, 0.0, 0.0], 1, True
        )
        crefim = EMAN2_cppwrap.Util.Polar2Dm(temp, cnx, cnx, numr, mode)
        EMAN2_cppwrap.Util.Normalize_ring(crefim, numr, 0)
        EMAN2_cppwrap.Util.Frngs(crefim, numr)
        EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)
        bigbuffer.insert_clip(crefim, (0, i))

    mpi.mpi_barrier(Blockdata["shared_comm"])
    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_global_def.sxprint(
            "  Reference projections generated : %10.1fmin"
            % (old_div((time.time() - at), 60.0))
        )

    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_global_def.sxprint("  ")
        # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
        line = ""
        sp_global_def.sxprint(
            line,
            "Processing data for polar onx: %3d, nx: %3d, nxpolar: %3d, CTF: %s, preshift: %s."
            % (
                Tracker["constants"]["nnxo"],
                Tracker["nxinit"],
                Tracker["nxpolar"],
                Tracker["constants"]["CTF"],
                preshift,
            ),
        )

    #  Preprocess the data

    #  Note these are in Fortran notation for polar searches
    txm = float(Tracker["nxpolar"] - (old_div(Tracker["nxpolar"], 2) + 1) - radius)
    txl = float(radius - old_div(Tracker["nxpolar"], 2) + 1)
    """Multiline Comment15"""

    lxod1 = n_coarse_ang * len(coarse_shifts) * n_coarse_psi
    newpar = [[i, [0.0], []] for i in range(nima)]

    #  This is for auxiliary function searches.
    Blockdata["angle_set"] = refang
    #  Extract normals from rotation matrices
    refdirs = sp_utilities.angles_to_normals(refang)

    # if( Blockdata["myid"] == Blockdata["main_node"]):
    # 	sxprint("  FIRST  nima, nang, npsi, nshift, n_coarse_ang, n_coarse_psi, len(list_of_coarse_shifts), lxod1 ",nima, nang,npsi,nshifts, n_coarse_ang,n_coarse_psi, len(coarse_shifts), lxod1)

    at = time.time()

    for im in range(nima):

        # phi,theta,psi,sx,sy, wnorm = oldparams[im][0], oldparams[im][1], oldparams[im][2], oldparams[im][3], oldparams[im][4], oldparams[im][7]
        phi, theta, psi, sx, sy = (
            oldparams[im][0],
            oldparams[im][1],
            oldparams[im][2],
            oldparams[im][3],
            oldparams[im][4],
        )  # First ITER

        if preshift:
            sx, sy = reduce_shifts(sx, sy, original_data[im])
            dataimage = sp_fundamentals.cyclic_shift(original_data[im], sx, sy)
            #  Put rounded shifts on the list, note image has the original floats - check whether it may cause problems
            oldparams[im][3] = sx
            oldparams[im][4] = sy
            sx = 0.0
            sy = 0.0
        else:
            dataimage = original_data[im].copy()

        st = get_image_statistics(dataimage, mask2D, False)
        dataimage -= st[0]
        dataimage = old_div(dataimage, st[1])
        ##if dataimage.get_attr_default("bckgnoise", None) :  dataimage.delete_attr("bckgnoise")
        #  Do bckgnoise if exists
        if apply_mask:
            if Tracker["constants"]["hardmask"]:
                dataimage = sp_morphology.cosinemask(
                    dataimage, radius=Tracker["constants"]["radius"]
                )
            else:
                bckg = sp_utilities.model_gauss_noise(
                    1.0, Tracker["constants"]["nnxo"] + 2, Tracker["constants"]["nnxo"]
                )
                bckg.set_attr("is_complex", 1)
                bckg.set_attr("is_fftpad", 1)
                bckg = sp_fundamentals.fft(
                    sp_filter.filt_table(
                        bckg, oneover[dataimage.get_attr("particle_group")]
                    )
                )
                #  Normalize bckg noise in real space, only region actually used.
                st = get_image_statistics(bckg, mask2D, False)
                bckg -= st[0]
                bckg = old_div(bckg, st[1])
                dataimage = sp_morphology.cosinemask(
                    dataimage, radius=Tracker["constants"]["radius"], bckg=bckg
                )
        # else:
        # 	#  if no bckgnoise, do simple masking instead
        # 	if apply_mask:  dataimage = cosinemask(dataimage,radius = Tracker["constants"]["radius"] )

        #  Apply varadj, FIRST ITERATION, THERE IS NONE
        # if not nonorm:
        # 	Util.mul_scalar(dataimage, Tracker["avgvaradj"][procid]/wnorm)

        ###  FT
        dataimage = sp_fundamentals.fft(dataimage)
        ##sig = Util.rotavg_fourier( dataimage )
        ##accumulatepw[im] = sig[len(sig)//2:]+[0.0]

        #  We have to make sure the shifts are within correct range, shrinkage or not
        # set_params_proj(dataimage,[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
        """Multiline Comment16"""
        if Tracker["constants"]["CTF"]:
            ctf_params = dataimage.get_attr("ctf")
            ctf_params.apix = old_div(
                ctf_params.apix,
                (
                    old_div(
                        float(Tracker["nxinit"]), float(Tracker["constants"]["nnxo"])
                    )
                ),
            )
            ctfa = sp_morphology.ctf_img_real(Tracker["nxinit"], ctf_params)
            ctfs = ctfa
        dataml = sp_fundamentals.fdecimate(
            dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False, False
        )
        data = []
        for iq in coarse_shifts:
            xx = iq[0] * shrink
            yy = iq[1] * shrink
            dss = sp_fundamentals.fshift(dataml, xx, yy)
            dss.set_attr("is_complex", 0)
            data.append(dss)

        #  This will get it to real space
        # dataimage = fpol(Util.mulnclreal(Util.mulnclreal(fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False), Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)
        dataimage = sp_fundamentals.fpol(
            EMAN2_cppwrap.Util.mulnclreal(
                sp_fundamentals.fdecimate(
                    dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False
                ),
                mask,
            ),
            Tracker["nxpolar"],
            Tracker["nxpolar"],
            1,
            True,
        )

        if (
            (Blockdata["myid"] == Blockdata["main_node"])
            and (im % (max(1, old_div(nima, 5))) == 0)
            and (im > 0)
        ):
            sp_global_def.sxprint(
                "  Number of images :%7d   %5d  %5.1f"
                % (im, nima, old_div(float(im), float(nima)) * 100.0)
                + "%"
                + "   %10.1fmin" % (old_div((time.time() - at), 60.0))
            )

        if im < min(nima, 1) and procid == 0:
            #   Search all and compare with direct to figure what keepfirst might be
            keepfirst = old_div(
                (n_coarse_ang * n_coarse_psi), 10
            )  # keepfirst = (n_coarse_ang *  n_coarse_psi * n_coarse_shifts)/10

            xod2 = sp_helix_sphire.np.asarray(
                EMAN2_cppwrap.Util.multiref_Crosrng_msg_stack_stepsi(
                    dataimage,
                    bigbuffer,
                    coarse_shifts_shrank,
                    numr,
                    [coarse_angles[ic][2] for ic in range(n_coarse_ang)],
                    coarse_delta,
                    cnx,
                    keepfirst,
                )
            )

            assert len(xod2) == keepfirst

            xod1 = sp_helix_sphire.np.ndarray((keepfirst), dtype="f4", order="C")

            for iln in range(keepfirst):
                m = xod2[iln]
                j = m % n_coarse_psi
                ic = (old_div(m, n_coarse_psi)) % n_coarse_ang
                ib = old_div(m, (n_coarse_ang * n_coarse_psi))
                xod2[iln] = j * 1000 + ic * 100000000 + ib  # hashparams
            # DO NOT order by angular directions to save time on reprojections.
            pre_ipsiandiang = -1
            for iln in range(keepfirst):
                hashparams = int(xod2[iln])
                ishift = hashparams % 1000
                ipsiandiang = old_div(hashparams, 1000)
                if ipsiandiang != pre_ipsiandiang:
                    pre_ipsiandiang = ipsiandiang
                    ipsi = ipsiandiang % 100000
                    iang = old_div(ipsiandiang, 100000)
                    temp = sp_projection.prgl(
                        volinit,
                        [
                            coarse_angles[iang][0],
                            coarse_angles[iang][1],
                            (coarse_angles[iang][2] + ipsi * coarse_delta) % 360.0,
                            0.0,
                            0.0,
                        ],
                        1,
                        False,
                    )
                    nrmref = numpy.sqrt(
                        EMAN2_cppwrap.Util.innerproduct(temp, temp, None)
                    )
                    # Util.mul_scalar(temp, 1.0/nrmref)

                xod1[iln] = old_div(
                    EMAN2_cppwrap.Util.innerproduct(data[ishift], temp, None), nrmref
                )
                # peak
                ##xod2[iln] = hashparams

            z = sp_helix_sphire.np.max(xod1)
            lina = sp_helix_sphire.np.argwhere(xod1 == z)[0]
            # if( Blockdata["myid"] == 0 ):
            # print("  STARTING ",Blockdata["myid"],z,lina)
            keepf = (
                int(lina[0]) + 1
            )  # if position of the highest is on lina[0] it means we have to keep lina[0] + 1 elements
            xod1 = xod1[lina]
            xod2 = xod2[lina]

            """Multiline Comment17"""
            lit = 1
            keepf = [keepf]
            keepf = sp_utilities.wrap_mpi_gatherv(
                keepf, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            if Blockdata["myid"] == 0:
                keepf.sort()
                keepf = keepf[int(len(keepf) * Blockdata["rkeepf"])]
            else:
                keepf = 0
            keepf = sp_utilities.wrap_mpi_bcast(
                keepf, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            Tracker["keepfirst"] = int(keepf)
            ###if( Blockdata["myid"] == 0 ):  sxprint("  keepfirst first ",Tracker["keepfirst"])

        else:
            # Tracker["keepfirst"] = min(200,nang)#min(max(lxod1/25,200),lxod1)

            xod2 = sp_helix_sphire.np.asarray(
                EMAN2_cppwrap.Util.multiref_Crosrng_msg_stack_stepsi(
                    dataimage,
                    bigbuffer,
                    coarse_shifts_shrank,
                    numr,
                    [coarse_angles[ic][2] for ic in range(n_coarse_ang)],
                    coarse_delta,
                    cnx,
                    Tracker["keepfirst"],
                )
            )

            xod1 = sp_helix_sphire.np.ndarray(
                (Tracker["keepfirst"]), dtype="f4", order="C"
            )

            for iln in range(Tracker["keepfirst"]):
                m = xod2[iln]
                j = m % n_coarse_psi
                ic = (old_div(m, n_coarse_psi)) % n_coarse_ang
                ib = old_div(m, (n_coarse_ang * n_coarse_psi))
                xod2[iln] = j * 1000 + ic * 100000000 + ib  # hashparams
            # order by angular directions to save time on reprojections.
            ipsiandiang = old_div(xod2, 1000)
            lina = sp_helix_sphire.np.argsort(ipsiandiang)
            xod2 = xod2[lina]  # order does not matter
            pre_ipsiandiang = -1
            for iln in range(Tracker["keepfirst"]):
                hashparams = int(xod2[iln])
                ishift = hashparams % 1000
                ipsiandiang = old_div(hashparams, 1000)
                if ipsiandiang != pre_ipsiandiang:
                    pre_ipsiandiang = ipsiandiang
                    ipsi = ipsiandiang % 100000
                    iang = old_div(ipsiandiang, 100000)
                    temp = sp_projection.prgl(
                        volinit,
                        [
                            coarse_angles[iang][0],
                            coarse_angles[iang][1],
                            (coarse_angles[iang][2] + ipsi * coarse_delta) % 360.0,
                            0.0,
                            0.0,
                        ],
                        1,
                        False,
                    )
                    temp.set_attr("is_complex", 0)
                    nrmref = numpy.sqrt(
                        EMAN2_cppwrap.Util.innerproduct(temp, temp, None)
                    )
                peak = old_div(
                    EMAN2_cppwrap.Util.innerproduct(data[ishift], temp, None), nrmref
                )
                xod1[iln] = peak
                ##xod2[iln] = hashparams

            lina = sp_helix_sphire.np.argsort(xod1)
            xod1 = xod1[lina[::-1]]  # This sorts in reverse order
            xod2 = xod2[lina[::-1]]  # This sorts in reverse order
            lit = Tracker["keepfirst"]
            """Multiline Comment18"""

        del data

        # if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  SECOND KEPT  ",lit)

        firstdirections = [[0.0, 0.0] for iln in range(lit)]
        firstshifts = [0] * lit
        for iln in range(lit):
            hashparams = int(xod2[iln])
            ishift = hashparams % 1000
            ipsiandiang = old_div(hashparams, 1000)
            # ipsi = ipsiandiang%100000
            iang = old_div(ipsiandiang, 100000)
            firstdirections[iln] = [coarse_angles[iang][0], coarse_angles[iang][1], 0.0]
            firstshifts[iln] = ishift
        ###del xod2
        ###if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  SECOND KEPT  ",im,lit)#,len(xod2),xod2,firstdirections,firstshifts)
        # mpi_barrier(MPI_COMM_WORLD)
        # mpi_finalize()
        # exit()

        # Find neighbors, ltabang contains positions on refang list, no psis
        ltabang = find_nearest_k_refangles_to_many_angles(
            refdirs, firstdirections, Tracker["delta"], howmany=Tracker["howmany"]
        )

        # ltabang has length lit, which is the same as length as firshifts.  However, original xod2 is still available,
        #   even though it is longer than lit.

        #  Prepare image for chi2.
        #  We have to repeat everything from get shrink data, including shifts
        #  The number of entries is lit=len(ltabang), each proj dir has 4 neighbors
        #    there are 2 psis, and at most n_fine_shifts. which should be 4.
        #    Not not all shifts have to be populated. We may also have repeated triplets of angles, but each can have
        #      different shifts.  If so, we have to remove duplicates from the entire set.
        lcod1 = lit * 4 * 2 * n_fine_shifts
        cod2 = []
        # lol = 0
        for i1 in range(lit):
            hashparams = int(xod2[i1])
            ipsiandiang = old_div(hashparams, 1000)
            oldiang = old_div(ipsiandiang, 100000)
            ipsi = ipsiandiang % 100000
            ishift = hashparams % 1000
            tshifts = get_shifts_neighbors(shifts, coarse_shifts[ishift])
            for i2 in range(Tracker["howmany"]):
                iang = ltabang[i1][i2]
                for i3 in range(2):  # psi
                    itpsi = int(
                        old_div(
                            (
                                coarse_angles[oldiang][2]
                                + ipsi * coarse_delta
                                - refang[iang][2]
                                + 360.0
                            ),
                            Tracker["delta"],
                        )
                    )
                    itpsi = (itpsi + i3) % npsi
                    for i4 in range(len(tshifts)):
                        cod2.append(iang * 100000000 + itpsi * 1000 + tshifts[i4])

        del xod1, xod2

        # if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  THIRD1   ",len(cod2),cod2)
        cod2 = list(set(cod2))
        cod1 = [[old_div(q, 1000), i] for i, q in enumerate(cod2)]
        cod1.sort()

        lit = len(cod1)

        cod2 = sp_helix_sphire.np.asarray([cod2[cod1[i][1]] for i in range(lit)])

        cod1 = sp_helix_sphire.np.ndarray(lit, dtype="f4", order="C")
        # cod1.fill(np.finfo(dtype='f4').min)
        # cod3 = np.ndarray(lit,dtype='f4',order="C")
        # cod3.fill(0.0)  #  varadj

        # if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  THIRD2   ",im,lit,cod2)

        # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        #  Third ML part.  Now the actual chi2 distances, we compute reprojections on the fly.

        data = [None] * nshifts
        johi = 0
        iln = 0
        prevdir = -1
        while iln < lit:
            hashparams = cod2[iln]
            ipsiandiang = old_div(hashparams, 1000)
            if ipsiandiang != prevdir:
                prevdir = ipsiandiang
                ipsi = ipsiandiang % 100000
                iang = old_div(ipsiandiang, 100000)
                temp = sp_projection.prgl(
                    volinit,
                    [
                        refang[iang][0],
                        refang[iang][1],
                        (refang[iang][2] + ipsi * Tracker["delta"]) % 360.0,
                        0.0,
                        0.0,
                    ],
                    1,
                    False,
                )
                temp.set_attr("is_complex", 0)
                nrmref = numpy.sqrt(EMAN2_cppwrap.Util.innerproduct(temp, temp, None))
                johi += 1
            while ipsiandiang == old_div(cod2[iln], 1000):
                hashparams = cod2[iln]
                ishift = hashparams % 1000
                if data[ishift] == None:
                    xx = shifts[ishift][0] * shrink
                    yy = shifts[ishift][1] * shrink
                    data[ishift] = sp_fundamentals.fshift(dataml, xx, yy)
                    data[ishift].set_attr("is_complex", 0)

                cod1[iln] = old_div(
                    EMAN2_cppwrap.Util.innerproduct(data[ishift], temp, None), nrmref
                )
                iln += 1
                if iln == lit:
                    break

        del data
        del dataml

        lina = sp_helix_sphire.np.argsort(cod1)
        cod1 = cod1[lina[::-1]]  # This sorts in reverse order
        cod2 = cod2[lina[::-1]]  # This sorts in reverse order
        ##cod3 = cod3[lina[::-1]]  # This sorts in reverse order
        ##cod1 -= cod1[0]
        ##lina = np.argwhere(cod1 > Tracker["constants"]["expthreshold"])
        cod1 = cod1[0]
        cod2 = cod2[0]
        # = cod3[lina]
        """Multiline Comment19"""
        #  New norm is a sum of eq distances multiplied by their probabilities augmented by PW.
        ###norm_per_particle[im] = np.sum(cod1[:lit]*cod3[:lit]) + accumulatepw[im][reachpw]

        ###for iln in range(lit):
        newpar[im][2] = [[int(cod2), float(cod1)]]

        del cod1, cod2, lina
        ###mpi_barrier(MPI_COMM_WORLD)
        ###mpi_finalize()
        ###exit()
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_global_def.sxprint(
            "  Finished projection matching   %10.1fmin"
            % (old_div((time.time() - at), 60.0))
        )
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    #  All images were processed, now to the additional calculations
    ###mpi_barrier(MPI_COMM_WORLD)
    ###mpi_finalize()
    ###exit()

    """Multiline Comment20"""

    at = time.time()
    mpi.mpi_barrier(Blockdata["shared_comm"])

    ###if Blockdata["myid"] == Blockdata["main_node"]:  sxprint "  Finished :",time()-at
    mpi.mpi_win_free(win_sm)
    if Tracker["nxinit"] != Tracker["nxpolar"]:
        mpi.mpi_win_free(win_volinit)
        emnumpy4.unregister_numpy_from_emdata()
        del emnumpy4
    else:
        mpi.mpi_win_free(win_vol)

    mpi.mpi_barrier(Blockdata["shared_comm"])
    emnumpy1.unregister_numpy_from_emdata()
    emnumpy2.unregister_numpy_from_emdata()
    del emnumpy1, emnumpy2

    del volinit

    # print("  NORMALIZATION DONE  ",Blockdata["myid"])
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    # if( Blockdata["myid"] == Blockdata["main_node"] ):
    # 	#write_text_row([[newpar[0][2][j][0],newpar[0][2][j][1]] for j in range(len(newpar[0][2]))],os.path.join(Tracker["directory"], "polar%1d.txt"%procid))
    # 	sxprint( "  Statistics finished : %10.1fmin"%((time()-at)/60.))
    return newpar, [1.0] * nima  # norm_per_particle


def ali3D_primary_polar(
    refang,
    shifts,
    coarse_angles,
    coarse_shifts,
    procid,
    original_data=None,
    oldparams=None,
    preshift=False,
    apply_mask=True,
    nonorm=False,
    applyctf=True,
):
    global Tracker, Blockdata
    #  Input data has to be CTF-multiplied, preshifted
    #  Output - newpar, see structure
    #    newpar = [[i, [worst_similarity, sum_all_similarities], [[-1, -1.0e23] for j in range(Tracker["lentop"])]] for i in range(len(data))]
    #    newpar = [[params],[],... len(data)]
    #    params = [particleID, [worst_similarity, sum_all_similarities],[imageallparams]]]
    #    imageallparams = [[orientation, similarity],[],...  number of all orientations ]
    #  Coding of orientations:
    #    hash = ang*100000000 + lpsi*1000 + ishift
    #    ishift = hash%1000
    #    ipsi = (hash/1000)%100000
    #    iang  = hash/100000000
    #  To get best matching for particle #kl:
    #     hash_best = newpar[kl][-1][0][0]
    #     best_sim  = newpar[kl][-1][0][1]
    #  To sort:
    # from operator 		import itemgetter, attrgetter, methodcaller
    #   params.sort(key=itemgetter(2))
    # from fundamentals import resample

    if Blockdata["myid"] == Blockdata["main_node"]:
        print_dict(Tracker, "PROJECTION MATCHING parameters of primary polar CCC")

    at = time.time()
    shrinkage = old_div(float(Tracker["nxpolar"]), float(Tracker["constants"]["nnxo"]))
    shrink = old_div(float(Tracker["nxinit"]), float(Tracker["constants"]["nnxo"]))
    radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
    mode = "F"
    numr = Numrinit_local(1, radius, 1, mode)
    wr = sp_alignment.ringwe(numr, mode)
    cnx = float(old_div(Tracker["nxpolar"], 2) + 1)

    ##if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  RADIUS IN POLAR  ",radius,numr[-1])
    #  FINE SEARCH CONSTANTS
    #  Fine grids are shifted by half-fine_step
    nang = len(refang)
    npsi = int(old_div(360.0, Tracker["delta"]))
    nshifts = len(shifts)
    n_fine_shifts = 4

    nima = len(original_data)
    mask = EMAN2_cppwrap.Util.unrollmask(Tracker["nxinit"], Tracker["nxinit"])
    for j in range(old_div(Tracker["nxinit"], 2), Tracker["nxinit"]):
        mask[0, j] = 1.0
    mask2D = sp_utilities.model_circle(
        Tracker["constants"]["radius"],
        Tracker["constants"]["nnxo"],
        Tracker["constants"]["nnxo"],
    )
    reachpw = (
        mask.get_xsize()
    )  # The last element of accumulated pw is zero so for the full size nothing is added.

    #  COARSE SEARCH CONSTANTS
    n_coarse_ang = len(coarse_angles)
    coarse_delta = 2 * Tracker["delta"]
    n_coarse_psi = int(old_div(360.0, coarse_delta))
    n_coarse_shifts = len(coarse_shifts)

    coarse_shifts_shrank = [None] * n_coarse_shifts
    for ib in range(n_coarse_shifts):
        coarse_shifts_shrank[ib] = [
            coarse_shifts[ib][0] * shrinkage,
            coarse_shifts[ib][1] * shrinkage,
        ]

    ny = Tracker["nxinit"]
    nyp2 = old_div(ny, 2)
    nxth = old_div((Tracker["nxinit"] + 2), 2)
    indx = sp_utilities.model_blank(nxth, Tracker["nxinit"], 1, -1)
    tfrac = sp_utilities.model_blank(nxth, Tracker["nxinit"])
    tcount = sp_utilities.model_blank(nxth)
    for iy in range(1, ny + 1):
        jy = iy - 1
        if jy > nyp2:
            jy = jy - ny
        argy = float(jy * jy)
        for ix in range(1, nxth + 1):
            jx = ix - 1
            roff = jx + (iy - 1) * nxth
            if mask[ix - 1, iy - 1] > 0.0:
                rf = numpy.sqrt(argy + float(jx * jx))
                ir = int(rf)
                # print  ix-1,iy-1,roff,mask[ix-1,iy-1],rf,ir

                if ir < nxth - 1:
                    frac = rf - float(ir)
                    qres = 1.0 - frac
                    tfrac[ix - 1, iy - 1] = frac
                    # ioff = 2*roff
                    tcount[ir] += qres
                    tcount[ir + 1] += frac
                    indx[ix - 1, iy - 1] = ir

    disp_unit = sp_helix_sphire.np.dtype("f4").itemsize

    #  REFVOL
    if Blockdata["myid_on_node"] == 0:
        if Blockdata["myid"] == Blockdata["main_node"]:
            odo = get_refvol(Tracker["nxpolar"])
            nxvol = odo.get_xsize()
            nyvol = odo.get_ysize()
            nzvol = odo.get_zsize()
        else:
            nxvol = 0
            nyvol = 0
            nzvol = 0

        nxvol = sp_utilities.bcast_number_to_all(
            nxvol,
            source_node=Blockdata["main_node"],
            mpi_comm=Blockdata["group_zero_comm"],
        )
        nyvol = sp_utilities.bcast_number_to_all(
            nyvol,
            source_node=Blockdata["main_node"],
            mpi_comm=Blockdata["group_zero_comm"],
        )
        nzvol = sp_utilities.bcast_number_to_all(
            nzvol,
            source_node=Blockdata["main_node"],
            mpi_comm=Blockdata["group_zero_comm"],
        )

        if Blockdata["myid"] != Blockdata["main_node"]:
            odo = sp_utilities.model_blank(nxvol, nyvol, nzvol)

        sp_utilities.bcast_EMData_to_all(
            odo,
            Blockdata["group_zero_myid"],
            source_node=Blockdata["main_node"],
            comm=Blockdata["group_zero_comm"],
        )

        odo = sp_projection.prep_vol(odo, npad=2, interpolation_method=1)
        ndo = EMAN2_cppwrap.EMNumPy.em2numpy(odo)
        nxvol = odo.get_xsize()
        nyvol = odo.get_ysize()
        nzvol = odo.get_zsize()
        orgsizevol = nxvol * nyvol * nzvol
        sizevol = orgsizevol
    else:
        orgsizevol = 0
        sizevol = 0
        nxvol = 0
        nyvol = 0
        nzvol = 0

    orgsizevol = sp_utilities.bcast_number_to_all(
        orgsizevol, source_node=Blockdata["main_node"]
    )
    nxvol = sp_utilities.bcast_number_to_all(nxvol, source_node=Blockdata["main_node"])
    nyvol = sp_utilities.bcast_number_to_all(nyvol, source_node=Blockdata["main_node"])
    nzvol = sp_utilities.bcast_number_to_all(nzvol, source_node=Blockdata["main_node"])

    win_vol, base_vol = mpi.mpi_win_allocate_shared(
        sizevol * disp_unit, disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"]
    )
    sizevol = orgsizevol
    if Blockdata["myid_on_node"] != 0:
        base_vol, = mpi.mpi_win_shared_query(win_vol, mpi.MPI_PROC_NULL)
    """
    Comments from Adnan:
    numpy.core.mutliarray.int_asbuffer function is not available anymore . That is why a different alternative is tried here.
    """
    # volbuf = sp_helix_sphire.np.frombuffer(
    #     sp_helix_sphire.np.core.multiarray.int_asbuffer(base_vol, sizevol * disp_unit),
    #     dtype="f4",
    # )
    ptr = ctypes.cast(base_vol, ctypes.POINTER(ctypes.c_int * sizevol))
    volbuf = numpy.frombuffer(ptr.contents, dtype="f4")
    volbuf = volbuf.reshape(nzvol, nyvol, nxvol)
    if Blockdata["myid_on_node"] == 0:
        sp_helix_sphire.np.copyto(volbuf, ndo)
        del odo, ndo

    # volprep = EMNumPy.assign_numpy_to_emdata(volbuf)
    emnumpy1 = EMAN2_cppwrap.EMNumPy()
    volprep = emnumpy1.register_numpy_to_emdata(volbuf)
    volprep.set_attr_dict(
        {
            "is_complex": 1,
            "is_complex_x": 0,
            "is_fftodd": 0,
            "is_fftpad": 1,
            "is_shuffled": 1,
            "npad": 2,
        }
    )
    crefim = EMAN2_cppwrap.Util.Polar2Dm(
        sp_utilities.model_blank(Tracker["nxpolar"], Tracker["nxpolar"]),
        cnx,
        cnx,
        numr,
        mode,
    )

    #  BIG BUFFER (for n_coarse_ang polar arrays)
    size_of_one_image = crefim.get_xsize()
    lenbigbuf = n_coarse_ang
    orgsize = (
        lenbigbuf * size_of_one_image
    )  # This is number of projections to be computed simultaneously times their size

    if Blockdata["myid_on_node"] == 0:
        size = orgsize
    else:
        size = 0

    win_sm, base_ptr = mpi.mpi_win_allocate_shared(
        size * disp_unit, disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"]
    )
    size = orgsize
    if Blockdata["myid_on_node"] != 0:
        base_ptr, = mpi.mpi_win_shared_query(win_sm, mpi.MPI_PROC_NULL)
    """
    Comments from Adnan:
    numpy.core.mutliarray.int_asbuffer function is not available anymore . That is why a different alternative is tried here.
    """
    # buffer = sp_helix_sphire.np.frombuffer(
    #     sp_helix_sphire.np.core.multiarray.int_asbuffer(base_ptr, size * disp_unit),
    #     dtype="f4",
    # )
    ptr = ctypes.cast(base_ptr, ctypes.POINTER(ctypes.c_int * size))
    buffer = numpy.frombuffer(ptr.contents, dtype="f4")
    buffer = buffer.reshape(lenbigbuf, size_of_one_image)
    # bigbuffer = EMNumPy.assign_numpy_to_emdata(buffer)

    emnumpy2 = EMAN2_cppwrap.EMNumPy()
    bigbuffer = emnumpy2.register_numpy_to_emdata(buffer)

    #  end of setup

    if Tracker["nxinit"] != Tracker["nxpolar"]:
        mpi.mpi_win_free(win_vol)
        #  REFVOL FOR ML
        if Blockdata["myid_on_node"] == 0:
            if Blockdata["myid"] == Blockdata["main_node"]:
                odo = get_refvol(Tracker["nxpolar"])
                nxvol = odo.get_xsize()
                nyvol = odo.get_ysize()
                nzvol = odo.get_zsize()
            else:
                nxvol = 0
                nyvol = 0
                nzvol = 0

            nxvol = sp_utilities.bcast_number_to_all(
                nxvol,
                source_node=Blockdata["main_node"],
                mpi_comm=Blockdata["group_zero_comm"],
            )
            nyvol = sp_utilities.bcast_number_to_all(
                nyvol,
                source_node=Blockdata["main_node"],
                mpi_comm=Blockdata["group_zero_comm"],
            )
            nzvol = sp_utilities.bcast_number_to_all(
                nzvol,
                source_node=Blockdata["main_node"],
                mpi_comm=Blockdata["group_zero_comm"],
            )

            if Blockdata["myid"] != Blockdata["main_node"]:
                odo = sp_utilities.model_blank(nxvol, nyvol, nzvol)

            sp_utilities.bcast_EMData_to_all(
                odo,
                Blockdata["group_zero_myid"],
                source_node=Blockdata["main_node"],
                comm=Blockdata["group_zero_comm"],
            )

            odo = sp_projection.prep_vol(odo, npad=2, interpolation_method=1)
            ndo = EMAN2_cppwrap.EMNumPy.em2numpy(odo)
            nxvol = odo.get_xsize()
            nyvol = odo.get_ysize()
            nzvol = odo.get_zsize()
            orgsizevol = nxvol * nyvol * nzvol
            sizevol = orgsizevol
        else:
            orgsizevol = 0
            sizevol = 0
            nxvol = 0
            nyvol = 0
            nzvol = 0

        orgsizevol = sp_utilities.bcast_number_to_all(
            orgsizevol, source_node=Blockdata["main_node"]
        )
        nxvol = sp_utilities.bcast_number_to_all(
            nxvol, source_node=Blockdata["main_node"]
        )
        nyvol = sp_utilities.bcast_number_to_all(
            nyvol, source_node=Blockdata["main_node"]
        )
        nzvol = sp_utilities.bcast_number_to_all(
            nzvol, source_node=Blockdata["main_node"]
        )

        win_volinit, base_volinit = mpi.mpi_win_allocate_shared(
            sizevol * disp_unit, disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"]
        )
        sizevol = orgsizevol
        if Blockdata["myid_on_node"] != 0:
            base_volinit, = mpi.mpi_win_shared_query(win_volinit, mpi.MPI_PROC_NULL)
        """
        Comments from Adnan:
        numpy.core.mutliarray.int_asbuffer function is not available anymore . That is why a different alternative is tried here.
        """
        # volbufinit = sp_helix_sphire.np.frombuffer(
        #     sp_helix_sphire.np.core.multiarray.int_asbuffer(
        #         base_volinit, sizevol * disp_unit
        #     ),
        #     dtype="f4",
        # )
        ptr = ctypes.cast(base_volinit, ctypes.POINTER(ctypes.c_int * sizevol))
        volbufinit = numpy.frombuffer(ptr.contents, dtype="f4")
        volbufinit = volbufinit.reshape(nzvol, nyvol, nxvol)
        if Blockdata["myid_on_node"] == 0:
            sp_helix_sphire.np.copyto(volbufinit, ndoinit)
            del odo, ndoinit

        # volinit = EMNumPy.assign_numpy_to_emdata(volbufinit)
        emnumpy4 = EMAN2_cppwrap.EMNumPy()
        volinit = emnumpy4.register_numpy_to_emdata(volbufinit)
        volinit.set_attr_dict(
            {
                "is_complex": 1,
                "is_complex_x": 0,
                "is_fftodd": 0,
                "is_fftpad": 1,
                "is_shuffled": 1,
                "npad": 2,
            }
        )
        if Blockdata["myid_on_node"] == 0:
            volinit.update()
        mpi.mpi_barrier(Blockdata["shared_comm"])
    else:
        volinit = volprep
    #  End of replaced volprep

    at = time.time()

    nang_start, nang_end = sp_applications.MPI_start_end(
        n_coarse_ang, Blockdata["no_of_processes_per_group"], Blockdata["myid_on_node"]
    )

    for i in range(
        nang_start, nang_end, 1
    ):  # This will take care of no of process on a node less than nang.  Some loops will not be executed
        temp = sp_projection.prgl(
            volprep, [coarse_angles[i][0], coarse_angles[i][1], 0.0, 0.0, 0.0], 1, True
        )
        crefim = EMAN2_cppwrap.Util.Polar2Dm(temp, cnx, cnx, numr, mode)
        EMAN2_cppwrap.Util.Frngs(crefim, numr)
        EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)
        bigbuffer.insert_clip(crefim, (0, i))

    mpi.mpi_barrier(Blockdata["shared_comm"])
    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_global_def.sxprint(
            "  Reference projections generated : %10.1fmin"
            % (old_div((time.time() - at), 60.0))
        )

    ###mpi_finalize()
    ###exit()
    ###  <><><><><><><><>

    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_global_def.sxprint("  ")
        # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
        line = ""
        sp_global_def.sxprint(
            line,
            "Processing data for polar onx: %3d, nx: %3d, nxpolar: %3d, CTF: %s, preshift: %s."
            % (
                Tracker["constants"]["nnxo"],
                Tracker["nxinit"],
                Tracker["nxpolar"],
                Tracker["constants"]["CTF"],
                preshift,
            ),
        )

    #  Preprocess the data

    #  Note these are in Fortran notation for polar searches
    txm = float(Tracker["nxpolar"] - (old_div(Tracker["nxpolar"], 2) + 1) - radius)
    txl = float(radius - old_div(Tracker["nxpolar"], 2) + 1)

    if Blockdata["bckgnoise"]:
        oneover = []
        nxb = Blockdata["bckgnoise"][0].get_xsize()
        nyb = len(Blockdata["bckgnoise"])
        for i in range(nyb):
            temp = [0.0] * nxb
            for k in range(nxb):
                if Blockdata["bckgnoise"][i].get_value_at(k) > 0.0:
                    temp[k] = old_div(
                        1.0, numpy.sqrt(Blockdata["bckgnoise"][i].get_value_at(k))
                    )
            oneover.append(temp)
        del temp

        if procid == 0:
            Blockdata["totprob"] = [0.0] * nyb
            Blockdata["newbckgnoise"] = sp_utilities.model_blank(nxb, nyb)

    accumulatepw = [None] * nima
    norm_per_particle = [None] * nima

    lxod1 = n_coarse_ang * len(coarse_shifts) * n_coarse_psi
    newpar = [[i, [0.0], []] for i in range(nima)]

    #  This is for auxiliary function searches.
    Blockdata["angle_set"] = refang
    #  Extract normals from rotation matrices
    refdirs = sp_utilities.angles_to_normals(refang)

    # if( Blockdata["myid"] == Blockdata["main_node"]):
    # 	sxprint("  FIRST  nima, nang, npsi, nshift, n_coarse_ang, n_coarse_psi, len(list_of_coarse_shifts), lxod1 ",nima, nang,npsi,nshifts, n_coarse_ang,n_coarse_psi, len(coarse_shifts), lxod1)

    at = time.time()

    for im in range(nima):
        particle_group = original_data[im].get_attr("particle_group")

        phi, theta, psi, sx, sy, wnorm = (
            oldparams[im][0],
            oldparams[im][1],
            oldparams[im][2],
            oldparams[im][3],
            oldparams[im][4],
            oldparams[im][7],
        )

        if preshift:
            sx, sy = reduce_shifts(sx, sy, original_data[im])
            dataimage = sp_fundamentals.cyclic_shift(original_data[im], sx, sy)
            #  Put rounded shifts on the list, note image has the original floats - check whether it may cause problems
            oldparams[im][3] = sx
            oldparams[im][4] = sy
            sx = 0.0
            sy = 0.0
        else:
            dataimage = original_data[im].copy()

        st = get_image_statistics(dataimage, mask2D, False)
        dataimage -= st[0]
        dataimage = old_div(dataimage, st[1])
        if dataimage.get_attr_default("bckgnoise", None):
            dataimage.delete_attr("bckgnoise")
        #  Do bckgnoise if exists
        if apply_mask:
            if Tracker["constants"]["hardmask"]:
                dataimage = sp_morphology.cosinemask(
                    dataimage, radius=Tracker["constants"]["radius"]
                )
            else:
                bckg = sp_utilities.model_gauss_noise(
                    1.0, Tracker["constants"]["nnxo"] + 2, Tracker["constants"]["nnxo"]
                )
                bckg.set_attr("is_complex", 1)
                bckg.set_attr("is_fftpad", 1)
                bckg = sp_fundamentals.fft(
                    sp_filter.filt_table(bckg, oneover[particle_group])
                )
                #  Normalize bckg noise in real space, only region actually used.
                st = get_image_statistics(bckg, mask2D, False)
                bckg -= st[0]
                bckg = old_div(bckg, st[1])
                dataimage = sp_morphology.cosinemask(
                    dataimage, radius=Tracker["constants"]["radius"], bckg=bckg
                )
        # else:
        # 	#  if no bckgnoise, do simple masking instead
        # 	if apply_mask:  dataimage = cosinemask(dataimage,radius = Tracker["constants"]["radius"] )

        #  Apply varadj
        if not nonorm:
            EMAN2_cppwrap.Util.mul_scalar(
                dataimage, old_div(Tracker["avgvaradj"][procid], wnorm)
            )

        ###  FT
        dataimage = sp_fundamentals.fft(dataimage)
        sig = EMAN2_cppwrap.Util.rotavg_fourier(dataimage)
        accumulatepw[im] = sig[old_div(len(sig), 2) :] + [0.0]

        #  We have to make sure the shifts are within correct range, shrinkage or not
        # set_params_proj(dataimage,[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
        if Blockdata["bckgnoise"]:
            temp = Blockdata["bckgnoise"][particle_group]
            bckgn = EMAN2_cppwrap.Util.unroll1dpw(
                Tracker["nxinit"],
                Tracker["nxinit"],
                [temp[i] for i in range(temp.get_xsize())],
            )
        else:
            bckgn = EMAN2_cppwrap.Util.unroll1dpw(
                Tracker["nxinit"], Tracker["nxinit"], [1.0] * 600
            )
        bckgnoise = bckgn.copy()
        for j in range(old_div(Tracker["nxinit"], 2) + 1, Tracker["nxinit"]):
            bckgn[0, j] = bckgn[0, Tracker["nxinit"] - j]

        if Tracker["constants"]["CTF"]:
            ctf_params = dataimage.get_attr("ctf")
            ctf_params.apix = old_div(
                ctf_params.apix,
                (
                    old_div(
                        float(Tracker["nxinit"]), float(Tracker["constants"]["nnxo"])
                    )
                ),
            )
            ctfa = sp_morphology.ctf_img_real(Tracker["nxinit"], ctf_params)
            ctfs = ctfa
        dataml = sp_fundamentals.fdecimate(
            dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False, False
        )
        data = []
        for iq in coarse_shifts:
            xx = iq[0] * shrink
            yy = iq[1] * shrink
            dss = sp_fundamentals.fshift(dataml, xx, yy)
            dss.set_attr("is_complex", 0)
            data.append(dss)

        #  This will get it to real space
        # dataimage = fpol(Util.mulnclreal(Util.mulnclreal(fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False), Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)
        #  Here we can reuse dataml to save time.  It was not normalized, but it does not matter as dataimage is for POLAR.  PAP 08/06/2018
        dataimage = sp_fundamentals.fpol(
            EMAN2_cppwrap.Util.mulnclreal(
                EMAN2_cppwrap.Util.mulnclreal(
                    dataml, EMAN2_cppwrap.Util.muln_img(bckgn, ctfs)
                ),
                mask,
            ),
            Tracker["nxpolar"],
            Tracker["nxpolar"],
            1,
            True,
        )

        if (
            (Blockdata["myid"] == Blockdata["main_node"])
            and (im % (max(1, old_div(nima, 5))) == 0)
            and (im > 0)
        ):
            sp_global_def.sxprint(
                "  Number of images :%7d   %5d  %5.1f"
                % (im, nima, old_div(float(im), float(nima)) * 100.0)
                + "%"
                + "   %10.1fmin" % (old_div((time.time() - at), 60.0))
            )

        if im < min(nima, 1) and procid == 0:
            #   Search all and compare with direct to figure what keepfirst might be
            keepfirst = old_div(
                (n_coarse_ang * n_coarse_psi), 10
            )  # keepfirst = (n_coarse_ang *  n_coarse_psi * n_coarse_shifts)/10

            xod2 = sp_helix_sphire.np.asarray(
                EMAN2_cppwrap.Util.multiref_Crosrng_msg_stack_stepsi(
                    dataimage,
                    bigbuffer,
                    coarse_shifts_shrank,
                    numr,
                    [coarse_angles[ic][2] for ic in range(n_coarse_ang)],
                    coarse_delta,
                    cnx,
                    keepfirst,
                )
            )

            assert len(xod2) == keepfirst

            xod1 = sp_helix_sphire.np.ndarray((keepfirst), dtype="f4", order="C")

            for iln in range(keepfirst):
                m = xod2[iln]
                j = m % n_coarse_psi
                ic = (old_div(m, n_coarse_psi)) % n_coarse_ang
                ib = old_div(m, (n_coarse_ang * n_coarse_psi))
                xod2[iln] = j * 1000 + ic * 100000000 + ib  # hashparams
            # DO NOT order by angular directions to save time on reprojections.
            pre_ipsiandiang = -1
            for iln in range(keepfirst):
                hashparams = int(xod2[iln])
                ishift = hashparams % 1000
                ipsiandiang = old_div(hashparams, 1000)
                if ipsiandiang != pre_ipsiandiang:
                    pre_ipsiandiang = ipsiandiang
                    ipsi = ipsiandiang % 100000
                    iang = old_div(ipsiandiang, 100000)
                    temp = sp_projection.prgl(
                        volinit,
                        [
                            coarse_angles[iang][0],
                            coarse_angles[iang][1],
                            (coarse_angles[iang][2] + ipsi * coarse_delta) % 360.0,
                            0.0,
                            0.0,
                        ],
                        1,
                        False,
                    )
                    temp.set_attr("is_complex", 0)
                xod1[iln] = -EMAN2_cppwrap.Util.sqed(
                    data[ishift], temp, ctfa, bckgnoise
                )  # peak
                ##xod2[iln] = hashparams

            xod1 -= sp_helix_sphire.np.max(xod1)
            lina = sp_helix_sphire.np.argwhere(
                xod1 > Tracker["constants"]["expthreshold"]
            )
            # if( Blockdata["myid"] == 0 ):  sxprint("  STARTING ",np.max(xod1),np.min(xod1),len(lina),lina[-1])
            lina = lina.reshape(lina.size)
            keepf = int(lina[-1]) + 1

            xod1 = xod1[lina]
            xod2 = xod2[lina]

            lina = sp_helix_sphire.np.argsort(xod1)
            xod1 = xod1[lina[::-1]]  # This sorts in reverse order
            xod2 = xod2[lina[::-1]]  # This sorts in reverse order
            sp_helix_sphire.np.exp(xod1, out=xod1)
            xod1 = old_div(xod1, sp_helix_sphire.np.sum(xod1))
            cumprob = 0.0
            lit = len(xod1)
            for j in range(len(xod1)):
                cumprob += xod1[j]
                if cumprob > Tracker["constants"]["ccfpercentage"]:
                    lit = j + 1
                    break
            # keepf = mpi_reduce(keepf, 1, MPI_INT, MPI_MAX, Blockdata["main_node"], MPI_COMM_WORLD)
            # if( Blockdata["myid"] == 0 ):
            # 	keepf = max(int(float(keepf)*0.9),1)
            keepf = [keepf]
            keepf = sp_utilities.wrap_mpi_gatherv(
                keepf, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            if Blockdata["myid"] == 0:
                keepf.sort()
                keepf = keepf[int(len(keepf) * Blockdata["rkeepf"])]
            else:
                keepf = 0
            keepf = sp_utilities.wrap_mpi_bcast(
                keepf, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            Tracker["keepfirst"] = int(keepf)
            ###if( Blockdata["myid"] == 0 ):  sxprint("  keepfirst first ",Tracker["keepfirst"])

        else:
            # Tracker["keepfirst"] = min(200,nang)#min(max(lxod1/25,200),lxod1)

            xod2 = sp_helix_sphire.np.asarray(
                EMAN2_cppwrap.Util.multiref_Crosrng_msg_stack_stepsi(
                    dataimage,
                    bigbuffer,
                    coarse_shifts_shrank,
                    numr,
                    [coarse_angles[ic][2] for ic in range(n_coarse_ang)],
                    coarse_delta,
                    cnx,
                    Tracker["keepfirst"],
                )
            )

            xod1 = sp_helix_sphire.np.ndarray(
                (Tracker["keepfirst"]), dtype="f4", order="C"
            )

            for iln in range(Tracker["keepfirst"]):
                m = xod2[iln]
                j = m % n_coarse_psi
                ic = (old_div(m, n_coarse_psi)) % n_coarse_ang
                ib = old_div(m, (n_coarse_ang * n_coarse_psi))
                xod2[iln] = j * 1000 + ic * 100000000 + ib  # hashparams
            # order by angular directions to save time on reprojections.
            ipsiandiang = old_div(xod2, 1000)
            lina = sp_helix_sphire.np.argsort(ipsiandiang)
            xod2 = xod2[lina]  # order does not matter
            pre_ipsiandiang = -1
            for iln in range(Tracker["keepfirst"]):
                hashparams = int(xod2[iln])
                ishift = hashparams % 1000
                ipsiandiang = old_div(hashparams, 1000)
                if ipsiandiang != pre_ipsiandiang:
                    pre_ipsiandiang = ipsiandiang
                    ipsi = ipsiandiang % 100000
                    iang = old_div(ipsiandiang, 100000)
                    temp = sp_projection.prgl(
                        volinit,
                        [
                            coarse_angles[iang][0],
                            coarse_angles[iang][1],
                            (coarse_angles[iang][2] + ipsi * coarse_delta) % 360.0,
                            0.0,
                            0.0,
                        ],
                        1,
                        False,
                    )
                    temp.set_attr("is_complex", 0)
                peak = -EMAN2_cppwrap.Util.sqed(data[ishift], temp, ctfa, bckgnoise)
                xod1[iln] = peak
                ##xod2[iln] = hashparams

            lina = sp_helix_sphire.np.argsort(xod1)
            xod1 = xod1[lina[::-1]]  # This sorts in reverse order
            xod2 = xod2[lina[::-1]]  # This sorts in reverse order

            xod1 -= xod1[0]

            lina = sp_helix_sphire.np.argwhere(
                xod1 > Tracker["constants"]["expthreshold"]
            )
            xod1 = xod1[lina]
            xod2 = xod2[lina]
            sp_helix_sphire.np.exp(xod1, out=xod1)
            xod1 = old_div(xod1, sp_helix_sphire.np.sum(xod1))
            cumprob = 0.0
            lit = len(xod1)
            for j in range(len(xod1)):
                cumprob += xod1[j]
                if cumprob > Tracker["constants"]["ccfpercentage"]:
                    lit = j + 1
                    break

        del data

        # if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  SECOND KEPT  ",lit)

        firstdirections = [[0.0, 0.0] for iln in range(lit)]
        firstshifts = [0] * lit
        for iln in range(lit):
            hashparams = int(xod2[iln])
            ishift = hashparams % 1000
            ipsiandiang = old_div(hashparams, 1000)
            # ipsi = ipsiandiang%100000
            iang = old_div(ipsiandiang, 100000)
            firstdirections[iln] = [coarse_angles[iang][0], coarse_angles[iang][1], 0.0]
            firstshifts[iln] = ishift

        # Find neighbors, ltabang contains positions on refang list, no psis
        ltabang = find_nearest_k_refangles_to_many_angles(
            refdirs, firstdirections, Tracker["delta"], howmany=Tracker["howmany"]
        )

        # ltabang has length lit, which is the same as length as firshifts.  However, original xod2 is still available,
        #   even though it is longer than lit.

        #  Prepare image for chi2.
        #  We have to repeat everything from get shrink data, including shifts
        #  The number of entries is lit=len(ltabang), each proj dir has 4 neighbors
        #    there are 2 psis, and at most n_fine_shifts. which should be 4.
        #    Not not all shifts have to be populated. We may also have repeated triplets of angles, but each can have
        #      different shifts.  If so, we have to remove duplicates from the entire set.
        lcod1 = lit * 4 * 2 * n_fine_shifts
        cod2 = []
        for i1 in range(lit):
            hashparams = int(xod2[i1])
            ipsiandiang = old_div(hashparams, 1000)
            oldiang = old_div(ipsiandiang, 100000)
            ipsi = ipsiandiang % 100000
            ishift = hashparams % 1000
            tshifts = get_shifts_neighbors(shifts, coarse_shifts[ishift])
            for i2 in range(Tracker["howmany"]):
                iang = ltabang[i1][i2]
                for i3 in range(2):  # psi
                    itpsi = int(
                        old_div(
                            (
                                coarse_angles[oldiang][2]
                                + ipsi * coarse_delta
                                - refang[iang][2]
                                + 360.0
                            ),
                            Tracker["delta"],
                        )
                    )
                    itpsi = (itpsi + i3) % npsi
                    for i4 in range(len(tshifts)):
                        cod2.append(iang * 100000000 + itpsi * 1000 + tshifts[i4])

        del xod1, xod2

        # if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  THIRD1   ",len(cod2),cod2)
        cod2 = list(set(cod2))
        cod1 = [[old_div(q, 1000), i] for i, q in enumerate(cod2)]
        cod1.sort()

        lit = len(cod1)

        cod2 = sp_helix_sphire.np.asarray([cod2[cod1[i][1]] for i in range(lit)])

        cod1 = sp_helix_sphire.np.ndarray(lit, dtype="f4", order="C")
        # cod1.fill(np.finfo(dtype='f4').min)
        cod3 = sp_helix_sphire.np.ndarray(lit, dtype="f4", order="C")
        # cod3.fill(0.0)  #  varadj

        # if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  THIRD2   ",im,lit,cod2)

        # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        #  Third ML part.  Now the actual chi2 distances, we compute reprojections on the fly.
        #  Make sure volprep has nxinit size
        tbckg = []
        data = [None] * nshifts
        johi = 0
        iln = 0
        prevdir = -1
        while iln < lit:
            hashparams = cod2[iln]
            ipsiandiang = old_div(hashparams, 1000)
            if ipsiandiang != prevdir:
                prevdir = ipsiandiang
                ipsi = ipsiandiang % 100000
                iang = old_div(ipsiandiang, 100000)
                temp = sp_projection.prgl(
                    volinit,
                    [
                        refang[iang][0],
                        refang[iang][1],
                        (refang[iang][2] + ipsi * Tracker["delta"]) % 360.0,
                        0.0,
                        0.0,
                    ],
                    1,
                    False,
                )
                temp.set_attr("is_complex", 0)
                johi += 1
            while ipsiandiang == old_div(cod2[iln], 1000):
                hashparams = cod2[iln]
                ishift = hashparams % 1000
                if data[ishift] == None:
                    xx = shifts[ishift][0] * shrink
                    yy = shifts[ishift][1] * shrink
                    data[ishift] = sp_fundamentals.fshift(dataml, xx, yy)
                    data[ishift].set_attr("is_complex", 0)

                ##[peak,varadj] = Util.sqednorm(data[ishift], temp, ctfa, bckgnoise)
                fofa = EMAN2_cppwrap.Util.sqednormbckg(
                    data[ishift], temp, ctfa, bckgnoise, indx, tfrac, tcount
                )
                cod1[iln] = -fofa[-2]  # -peak
                cod3[iln] = fofa[-1]  # varadj
                tbckg.append(fofa[:-2])
                iln += 1
                if iln == lit:
                    break

        del data
        del dataml

        lina = sp_helix_sphire.np.argsort(cod1)
        lina = lina[::-1]  # This sorts in reverse order
        cod1 = cod1[lina]
        cod2 = cod2[lina]
        cod3 = cod3[lina]
        cod1 -= cod1[0]
        tbckg = [tbckg[int(q)] for q in lina]
        lina = sp_helix_sphire.np.argwhere(cod1 > Tracker["constants"]["expthreshold"])
        cod1 = cod1[lina]
        cod2 = cod2[lina]
        cod3 = cod3[lina]
        tbckg = [tbckg[int(q)] for q in lina]

        sp_helix_sphire.np.exp(cod1, out=cod1)
        cod1 = old_div(cod1, sp_helix_sphire.np.sum(cod1))
        cumprob = 0.0
        for j in range(len(cod1)):
            cumprob += cod1[j]
            if cumprob > Tracker["constants"]["ccfpercentage"]:
                lit = j + 1
                break

        #  New norm is a sum of eq distances multiplied by their probabilities augmented by PW.
        norm_per_particle[im] = (
            sp_helix_sphire.np.sum(cod1[:lit] * cod3[:lit]) + accumulatepw[im][reachpw]
        )
        atbckg = [0.0] * len(tbckg[0])
        for iln in range(lit):
            prob = float(cod1[iln])
            Blockdata["totprob"][particle_group] += prob
            for iq in range(len(tbckg[0])):
                atbckg[iq] += tbckg[iln][iq] * prob

        del tbckg
        for iq in range(nxth):
            Blockdata["newbckgnoise"][iq, particle_group] += atbckg[iq]
        del atbckg
        for iln in range(lit):
            newpar[im][2].append([int(cod2[iln]), float(cod1[iln])])

        del cod1, cod2, cod3, lina
        ###mpi_barrier(MPI_COMM_WORLD)
        ###mpi_finalize()
        ###exit()
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_global_def.sxprint(
            "  Finished projection matching   %10.1fmin"
            % (old_div((time.time() - at), 60.0))
        )
    at = time.time()
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    #  All images were processed, now to the additional calculations
    ###mpi_barrier(MPI_COMM_WORLD)
    ###mpi_finalize()
    ###exit()

    # norm correction ---- calc the norm correction per particle
    snormcorr = 0.0
    for kl in range(nima):
        norm_per_particle[kl] = old_div(
            numpy.sqrt(norm_per_particle[kl] * 2.0) * oldparams[kl][7],
            Tracker["avgvaradj"][procid],
        )
        snormcorr += norm_per_particle[kl]
    Tracker["avgvaradj"][procid] = snormcorr
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    #  Compute avgvaradj
    Tracker["avgvaradj"][procid] = mpi.mpi_reduce(
        Tracker["avgvaradj"][procid],
        1,
        mpi.MPI_FLOAT,
        mpi.MPI_SUM,
        Blockdata["main_node"],
        mpi.MPI_COMM_WORLD,
    )
    if Blockdata["myid"] == Blockdata["main_node"]:
        Tracker["avgvaradj"][procid] = old_div(
            float(Tracker["avgvaradj"][procid]), Tracker["nima_per_chunk"][procid]
        )
    else:
        Tracker["avgvaradj"][procid] = 0.0
    Tracker["avgvaradj"][procid] = sp_utilities.bcast_number_to_all(
        Tracker["avgvaradj"][procid], Blockdata["main_node"]
    )
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    #  Compute statistics of smear -----------------
    smax = -1000000
    smin = 1000000
    sava = 0.0
    svar = 0.0
    snum = 0
    for kl in range(nima):
        j = len(newpar[kl][2])
        snum += 1
        sava += float(j)
        svar += j * float(j)
        smax = max(smax, j)
        smin = min(smin, j)
    snum = mpi.mpi_reduce(
        snum, 1, mpi.MPI_INT, mpi.MPI_SUM, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    sava = mpi.mpi_reduce(
        sava, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    svar = mpi.mpi_reduce(
        svar, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    smax = mpi.mpi_reduce(
        smax, 1, mpi.MPI_INT, mpi.MPI_MAX, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    smin = mpi.mpi_reduce(
        smin, 1, mpi.MPI_INT, mpi.MPI_MIN, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    if Blockdata["myid"] == 0:
        sava = old_div(float(sava), snum)
        svar = numpy.sqrt(
            max(0.0, old_div((float(svar) - snum * sava ** 2), (snum - 1)))
        )
        # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
        line = ""
        sp_global_def.sxprint(
            line,
            "Smear stat  (number of images, ave, sumsq, min, max)):  %7d    %12.3g   %12.3g  %7d  %7d"
            % (snum, sava, svar, smin, smax),
        )

    at = time.time()
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    mpi.mpi_win_free(win_sm)
    if Tracker["nxinit"] != Tracker["nxpolar"]:
        mpi.mpi_win_free(win_volinit)
        emnumpy4.unregister_numpy_from_emdata()
        del emnumpy4
    else:
        mpi.mpi_win_free(win_vol)

    mpi.mpi_barrier(Blockdata["shared_comm"])
    emnumpy1.unregister_numpy_from_emdata()
    emnumpy2.unregister_numpy_from_emdata()
    del emnumpy1, emnumpy2

    del volinit

    # Compute new background noise
    # Reduce stuff
    if procid == 1:
        Blockdata["totprob"] = mpi.mpi_reduce(
            Blockdata["totprob"],
            nyb,
            mpi.MPI_FLOAT,
            mpi.MPI_SUM,
            Blockdata["main_node"],
            mpi.MPI_COMM_WORLD,
        )
        sp_utilities.reduce_EMData_to_root(
            Blockdata["newbckgnoise"], Blockdata["myid"], Blockdata["main_node"]
        )
        if Blockdata["myid"] == 0:
            for igrp in range(nyb):
                Blockdata["newbckgnoise"][0, igrp] = 1.0
                for i in range(1, nxth):
                    if Blockdata["newbckgnoise"][i, igrp] > 0.0:
                        Blockdata["newbckgnoise"][i, igrp] = old_div(
                            2.0 * Blockdata["totprob"][igrp],
                            Blockdata["newbckgnoise"][i, igrp],
                        )  # normalize and invert
                for i in range(nxth, nxb):
                    Blockdata["newbckgnoise"][i, igrp] = Blockdata["bckgnoise"][igrp][i]
            Blockdata["newbckgnoise"].write_image(
                os.path.join(Tracker["directory"], "bckgnoise.hdf")
            )  # Write updated bckgnoise to current directory

        sp_utilities.bcast_EMData_to_all(
            Blockdata["newbckgnoise"],
            Blockdata["myid"],
            source_node=Blockdata["main_node"],
            comm=mpi.MPI_COMM_WORLD,
        )
        for igrp in range(nyb):
            for i in range(nxb):
                Blockdata["bckgnoise"][igrp][i] = Blockdata["newbckgnoise"][i, igrp]
        del Blockdata["newbckgnoise"]
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_global_def.sxprint(
            "  Finished sigma2   %10.1fmin" % (old_div((time.time() - at), 60.0))
        )
    return newpar, norm_per_particle


def ali3D_polar(
    refang,
    shifts,
    coarse_angles,
    coarse_shifts,
    procid,
    original_data=None,
    oldparams=None,
    preshift=False,
    apply_mask=True,
    nonorm=False,
    applyctf=True,
):
    global Tracker, Blockdata
    #  Input data has to be CTF-multiplied, preshifted
    #  Output - newpar, see structure
    #    newpar = [[i, [worst_similarity, sum_all_similarities], [[-1, -1.0e23] for j in range(Tracker["lentop"])]] for i in range(len(data))]
    #    newpar = [[params],[],... len(data)]
    #    params = [particleID, [worst_similarity, sum_all_similarities],[imageallparams]]]
    #    imageallparams = [[orientation, similarity],[],...  number of all orientations ]
    #  Coding of orientations:
    #    hash = ang*100000000 + lpsi*1000 + ishift
    #    ishift = hash%1000
    #    ipsi = (hash/1000)%100000
    #    iang  = hash/100000000
    #  To get best matching for particle #kl:
    #     hash_best = newpar[kl][-1][0][0]
    #     best_sim  = newpar[kl][-1][0][1]
    #  To sort:
    #   params.sort(key=itemgetter(2))
    # from fundamentals import resample

    if Blockdata["myid"] == Blockdata["main_node"]:
        print_dict(Tracker, "PROJECTION MATCHING parameters of buffered exhaustive CCC")

    at = time.time()
    shrinkage = old_div(float(Tracker["nxpolar"]), float(Tracker["constants"]["nnxo"]))
    shrink = old_div(float(Tracker["nxinit"]), float(Tracker["constants"]["nnxo"]))
    radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
    mode = "F"
    numr = Numrinit_local(1, radius, 1, mode)
    wr = sp_alignment.ringwe(numr, mode)
    cnx = float(old_div(Tracker["nxpolar"], 2) + 1)

    ##if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  RADIUS IN POLAR  ",radius,numr[-1])
    #  FINE SEARCH CONSTANTS
    #  Fine grids are shifted by half-fine_step
    nang = len(refang)
    npsi = int(old_div(360.0, Tracker["delta"]))
    nshifts = len(shifts)
    n_fine_shifts = 4

    nima = len(original_data)
    mask = EMAN2_cppwrap.Util.unrollmask(Tracker["nxinit"], Tracker["nxinit"])
    for j in range(old_div(Tracker["nxinit"], 2), Tracker["nxinit"]):
        mask[0, j] = 1.0
    mask2D = sp_utilities.model_circle(
        Tracker["constants"]["radius"],
        Tracker["constants"]["nnxo"],
        Tracker["constants"]["nnxo"],
    )
    reachpw = (
        mask.get_xsize()
    )  # The last element of accumulated pw is zero so for the full size nothing is added.

    #  COARSE SEARCH CONSTANTS
    n_coarse_ang = len(coarse_angles)
    coarse_delta = 2 * Tracker["delta"]
    n_coarse_psi = int(old_div(360.0, coarse_delta))
    n_coarse_shifts = len(coarse_shifts)

    coarse_shifts_shrank = [None] * n_coarse_shifts
    for ib in range(n_coarse_shifts):
        coarse_shifts_shrank[ib] = [
            coarse_shifts[ib][0] * shrinkage,
            coarse_shifts[ib][1] * shrinkage,
        ]

    ###if(Blockdata["myid"] == Blockdata["main_node"]): sxprint("   TRETR  ",Tracker["constants"]["nnxo"],Tracker["nxinit"],reachpw)

    disp_unit = sp_helix_sphire.np.dtype("f4").itemsize

    #  REFVOL
    if Blockdata["myid_on_node"] == 0:
        if Blockdata["myid"] == Blockdata["main_node"]:
            odo = get_refvol(Tracker["nxpolar"])
            nxvol = odo.get_xsize()
            nyvol = odo.get_ysize()
            nzvol = odo.get_zsize()
        else:
            nxvol = 0
            nyvol = 0
            nzvol = 0

        nxvol = sp_utilities.bcast_number_to_all(
            nxvol,
            source_node=Blockdata["main_node"],
            mpi_comm=Blockdata["group_zero_comm"],
        )
        nyvol = sp_utilities.bcast_number_to_all(
            nyvol,
            source_node=Blockdata["main_node"],
            mpi_comm=Blockdata["group_zero_comm"],
        )
        nzvol = sp_utilities.bcast_number_to_all(
            nzvol,
            source_node=Blockdata["main_node"],
            mpi_comm=Blockdata["group_zero_comm"],
        )

        if Blockdata["myid"] != Blockdata["main_node"]:
            odo = sp_utilities.model_blank(nxvol, nyvol, nzvol)

        sp_utilities.bcast_EMData_to_all(
            odo,
            Blockdata["group_zero_myid"],
            source_node=Blockdata["main_node"],
            comm=Blockdata["group_zero_comm"],
        )

        odo = sp_projection.prep_vol(odo, npad=2, interpolation_method=1)
        ndo = EMAN2_cppwrap.EMNumPy.em2numpy(odo)
        nxvol = odo.get_xsize()
        nyvol = odo.get_ysize()
        nzvol = odo.get_zsize()
        orgsizevol = nxvol * nyvol * nzvol
        sizevol = orgsizevol
    else:
        orgsizevol = 0
        sizevol = 0
        nxvol = 0
        nyvol = 0
        nzvol = 0

    orgsizevol = sp_utilities.bcast_number_to_all(
        orgsizevol, source_node=Blockdata["main_node"]
    )
    nxvol = sp_utilities.bcast_number_to_all(nxvol, source_node=Blockdata["main_node"])
    nyvol = sp_utilities.bcast_number_to_all(nyvol, source_node=Blockdata["main_node"])
    nzvol = sp_utilities.bcast_number_to_all(nzvol, source_node=Blockdata["main_node"])

    win_vol, base_vol = mpi.mpi_win_allocate_shared(
        sizevol * disp_unit, disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"]
    )
    sizevol = orgsizevol
    if Blockdata["myid_on_node"] != 0:
        base_vol, = mpi.mpi_win_shared_query(win_vol, mpi.MPI_PROC_NULL)
    """
    Comments from Adnan:
    numpy.core.mutliarray.int_asbuffer function is not available anymore . That is why a different alternative is tried here.
    """
    # volbuf = sp_helix_sphire.np.frombuffer(
    #     sp_helix_sphire.np.core.multiarray.int_asbuffer(base_vol, sizevol * disp_unit),
    #     dtype="f4",
    # )
    ptr = ctypes.cast(base_vol, ctypes.POINTER(ctypes.c_int * sizevol))
    volbuf = numpy.frombuffer(ptr.contents, dtype="f4")
    volbuf = volbuf.reshape(nzvol, nyvol, nxvol)
    if Blockdata["myid_on_node"] == 0:
        sp_helix_sphire.np.copyto(volbuf, ndo)
        del odo, ndo

    # volprep = EMNumPy.assign_numpy_to_emdata(volbuf)
    emnumpy1 = EMAN2_cppwrap.EMNumPy()
    volprep = emnumpy1.register_numpy_to_emdata(volbuf)
    volprep.set_attr_dict(
        {
            "is_complex": 1,
            "is_complex_x": 0,
            "is_fftodd": 0,
            "is_fftpad": 1,
            "is_shuffled": 1,
            "npad": 2,
        }
    )

    crefim = EMAN2_cppwrap.Util.Polar2Dm(
        sp_utilities.model_blank(Tracker["nxpolar"], Tracker["nxpolar"]),
        cnx,
        cnx,
        numr,
        mode,
    )

    #  BIG BUFFER (for n_coarse_ang polar arrays)
    size_of_one_image = crefim.get_xsize()
    lenbigbuf = n_coarse_ang
    orgsize = (
        lenbigbuf * size_of_one_image
    )  # This is number of projections to be computed simultaneously times their size

    if Blockdata["myid_on_node"] == 0:
        size = orgsize
    else:
        size = 0

    win_sm, base_ptr = mpi.mpi_win_allocate_shared(
        size * disp_unit, disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"]
    )
    size = orgsize
    if Blockdata["myid_on_node"] != 0:
        base_ptr, = mpi.mpi_win_shared_query(win_sm, mpi.MPI_PROC_NULL)
    """
    Comments from Adnan:
    numpy.core.mutliarray.int_asbuffer function is not available anymore . That is why a different alternative is tried here.
    """
    # buffer = sp_helix_sphire.np.frombuffer(
    #     sp_helix_sphire.np.core.multiarray.int_asbuffer(base_ptr, size * disp_unit),
    #     dtype="f4",
    # )
    ptr = ctypes.cast(base_ptr, ctypes.POINTER(ctypes.c_int * size))
    buffer = numpy.frombuffer(ptr.contents, dtype="f4")
    buffer = buffer.reshape(lenbigbuf, size_of_one_image)
    # bigbuffer = EMNumPy.assign_numpy_to_emdata(buffer)

    emnumpy2 = EMAN2_cppwrap.EMNumPy()
    bigbuffer = emnumpy2.register_numpy_to_emdata(buffer)

    #  end of setup

    if Tracker["nxinit"] != Tracker["nxpolar"]:
        mpi.mpi_win_free(win_vol)
        #  REFVOL FOR ML
        if Blockdata["myid_on_node"] == 0:
            if Blockdata["myid"] == Blockdata["main_node"]:
                odo = get_refvol(Tracker["nxpolar"])
                nxvol = odo.get_xsize()
                nyvol = odo.get_ysize()
                nzvol = odo.get_zsize()
            else:
                nxvol = 0
                nyvol = 0
                nzvol = 0

            nxvol = sp_utilities.bcast_number_to_all(
                nxvol,
                source_node=Blockdata["main_node"],
                mpi_comm=Blockdata["group_zero_comm"],
            )
            nyvol = sp_utilities.bcast_number_to_all(
                nyvol,
                source_node=Blockdata["main_node"],
                mpi_comm=Blockdata["group_zero_comm"],
            )
            nzvol = sp_utilities.bcast_number_to_all(
                nzvol,
                source_node=Blockdata["main_node"],
                mpi_comm=Blockdata["group_zero_comm"],
            )

            if Blockdata["myid"] != Blockdata["main_node"]:
                odo = sp_utilities.model_blank(nxvol, nyvol, nzvol)

            sp_utilities.bcast_EMData_to_all(
                odo,
                Blockdata["group_zero_myid"],
                source_node=Blockdata["main_node"],
                comm=Blockdata["group_zero_comm"],
            )

            odo = sp_projection.prep_vol(odo, npad=2, interpolation_method=1)
            ndo = EMAN2_cppwrap.EMNumPy.em2numpy(odo)
            nxvol = odo.get_xsize()
            nyvol = odo.get_ysize()
            nzvol = odo.get_zsize()
            orgsizevol = nxvol * nyvol * nzvol
            sizevol = orgsizevol
        else:
            orgsizevol = 0
            sizevol = 0
            nxvol = 0
            nyvol = 0
            nzvol = 0

        orgsizevol = sp_utilities.bcast_number_to_all(
            orgsizevol, source_node=Blockdata["main_node"]
        )
        nxvol = sp_utilities.bcast_number_to_all(
            nxvol, source_node=Blockdata["main_node"]
        )
        nyvol = sp_utilities.bcast_number_to_all(
            nyvol, source_node=Blockdata["main_node"]
        )
        nzvol = sp_utilities.bcast_number_to_all(
            nzvol, source_node=Blockdata["main_node"]
        )

        win_volinit, base_volinit = mpi.mpi_win_allocate_shared(
            sizevol * disp_unit, disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"]
        )
        sizevol = orgsizevol
        if Blockdata["myid_on_node"] != 0:
            base_volinit, = mpi.mpi_win_shared_query(win_volinit, mpi.MPI_PROC_NULL)
        """
        Comments from Adnan:
        numpy.core.mutliarray.int_asbuffer function is not available anymore . That is why a different alternative is tried here.
        """
        # volbufinit = sp_helix_sphire.np.frombuffer(
        #     sp_helix_sphire.np.core.multiarray.int_asbuffer(
        #         base_volinit, sizevol * disp_unit
        #     ),
        #     dtype="f4",
        # )
        ptr = ctypes.cast(base_volinit, ctypes.POINTER(ctypes.c_int * sizevol))
        volbufinit = numpy.frombuffer(ptr.contents, dtype="f4")
        volbufinit = volbufinit.reshape(nzvol, nyvol, nxvol)
        if Blockdata["myid_on_node"] == 0:
            sp_helix_sphire.np.copyto(volbufinit, ndoinit)
            del odo, ndoinit

        # volinit = EMNumPy.assign_numpy_to_emdata(volbufinit)
        emnumpy4 = EMAN2_cppwrap.EMNumPy()
        volinit = emnumpy4.register_numpy_to_emdata(volbufinit)
        volinit.set_attr_dict(
            {
                "is_complex": 1,
                "is_complex_x": 0,
                "is_fftodd": 0,
                "is_fftpad": 1,
                "is_shuffled": 1,
                "npad": 2,
            }
        )
        if Blockdata["myid_on_node"] == 0:
            volinit.update()
        mpi.mpi_barrier(Blockdata["shared_comm"])
        ###if( Blockdata["myid"] < 10 ):
        ###	sxprint(" VOLPREP  ",volinit[0],volinit[1],info(volinit))
    else:
        volinit = volprep
    #  End of replaced volprep

    at = time.time()

    nang_start, nang_end = sp_applications.MPI_start_end(
        n_coarse_ang, Blockdata["no_of_processes_per_group"], Blockdata["myid_on_node"]
    )

    for i in range(
        nang_start, nang_end, 1
    ):  # This will take care of process on a node less than nang.  Some loops will not be executed
        temp = sp_projection.prgl(
            volprep, [coarse_angles[i][0], coarse_angles[i][1], 0.0, 0.0, 0.0], 1, True
        )
        crefim = EMAN2_cppwrap.Util.Polar2Dm(temp, cnx, cnx, numr, mode)
        EMAN2_cppwrap.Util.Frngs(crefim, numr)
        EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)
        bigbuffer.insert_clip(crefim, (0, i))

    mpi.mpi_barrier(Blockdata["shared_comm"])
    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_global_def.sxprint(
            "  Reference projections generated : %10.1fmin"
            % (old_div((time.time() - at), 60.0))
        )

    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_global_def.sxprint("  ")
        # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
        line = ""
        sp_global_def.sxprint(
            line,
            "Processing data for polar onx: %3d, nx: %3d, nxpolar: %3d, CTF: %s, preshift: %s."
            % (
                Tracker["constants"]["nnxo"],
                Tracker["nxinit"],
                Tracker["nxpolar"],
                Tracker["constants"]["CTF"],
                preshift,
            ),
        )

    #  Preprocess the data

    #  Note these are in Fortran notation for polar searches
    txm = float(Tracker["nxpolar"] - (old_div(Tracker["nxpolar"], 2) + 1) - radius)
    txl = float(radius - old_div(Tracker["nxpolar"], 2) + 1)

    if Blockdata["bckgnoise"]:
        oneover = []
        nnx = Blockdata["bckgnoise"][0].get_xsize()
        for i in range(len(Blockdata["bckgnoise"])):
            temp = [0.0] * nnx
            for k in range(nnx):
                if Blockdata["bckgnoise"][i].get_value_at(k) > 0.0:
                    temp[k] = old_div(
                        1.0, numpy.sqrt(Blockdata["bckgnoise"][i].get_value_at(k))
                    )
            oneover.append(temp)
        del temp

    accumulatepw = [None] * nima
    norm_per_particle = [None] * nima

    lxod1 = n_coarse_ang * len(coarse_shifts) * n_coarse_psi
    newpar = [[i, [0.0], []] for i in range(nima)]

    #  This is for auxiliary function searches.
    Blockdata["angle_set"] = refang
    #  Extract normals from rotation matrices
    refdirs = sp_utilities.angles_to_normals(refang)

    # if( Blockdata["myid"] == Blockdata["main_node"]):
    # 	sxprint("  FIRST  nima, nang, npsi, nshift, n_coarse_ang, n_coarse_psi, len(list_of_coarse_shifts), lxod1 ",nima, nang,npsi,nshifts, n_coarse_ang,n_coarse_psi, len(coarse_shifts), lxod1)

    ###firsti = True

    at = time.time()

    for im in range(nima):

        phi, theta, psi, sx, sy, wnorm = (
            oldparams[im][0],
            oldparams[im][1],
            oldparams[im][2],
            oldparams[im][3],
            oldparams[im][4],
            oldparams[im][7],
        )

        if preshift:
            sx, sy = reduce_shifts(sx, sy, original_data[im])
            dataimage = sp_fundamentals.cyclic_shift(original_data[im], sx, sy)
            #  Put rounded shifts on the list, note image has the original floats - check whether it may cause problems
            oldparams[im][3] = sx
            oldparams[im][4] = sy
            sx = 0.0
            sy = 0.0
        else:
            dataimage = original_data[im].copy()

        st = get_image_statistics(dataimage, mask2D, False)
        dataimage -= st[0]
        dataimage = old_div(dataimage, st[1])
        if dataimage.get_attr_default("bckgnoise", None):
            dataimage.delete_attr("bckgnoise")
        #  Do bckgnoise if exists
        # if Blockdata["bckgnoise"]:  # I think we should assume it always exists
        if apply_mask:
            if Tracker["constants"]["hardmask"]:
                dataimage = sp_morphology.cosinemask(
                    dataimage, radius=Tracker["constants"]["radius"]
                )
            else:
                bckg = sp_utilities.model_gauss_noise(
                    1.0, Tracker["constants"]["nnxo"] + 2, Tracker["constants"]["nnxo"]
                )
                bckg.set_attr("is_complex", 1)
                bckg.set_attr("is_fftpad", 1)
                bckg = sp_fundamentals.fft(
                    sp_filter.filt_table(
                        bckg, oneover[dataimage.get_attr("particle_group")]
                    )
                )
                #  Normalize bckg noise in real space, only region actually used.
                st = get_image_statistics(bckg, mask2D, False)
                bckg -= st[0]
                bckg = old_div(bckg, st[1])
                dataimage = sp_morphology.cosinemask(
                    dataimage, radius=Tracker["constants"]["radius"], bckg=bckg
                )
        # else:
        # 	#  if no bckgnoise, do simple masking instead
        # 	if apply_mask:  dataimage = cosinemask(dataimage,radius = Tracker["constants"]["radius"] )

        #  Apply varadj
        if not nonorm:
            EMAN2_cppwrap.Util.mul_scalar(
                dataimage, old_div(Tracker["avgvaradj"][procid], wnorm)
            )

        ###  FT
        dataimage = sp_fundamentals.fft(dataimage)
        sig = EMAN2_cppwrap.Util.rotavg_fourier(dataimage)
        accumulatepw[im] = sig[old_div(len(sig), 2) :] + [0.0]

        #  We have to make sure the shifts are within correct range, shrinkage or not
        # set_params_proj(dataimage,[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
        if Blockdata["bckgnoise"]:
            temp = Blockdata["bckgnoise"][dataimage.get_attr("particle_group")]
            bckgn = EMAN2_cppwrap.Util.unroll1dpw(
                Tracker["nxinit"],
                Tracker["nxinit"],
                [temp[i] for i in range(temp.get_xsize())],
            )
        else:
            bckgn = EMAN2_cppwrap.Util.unroll1dpw(
                Tracker["nxinit"], Tracker["nxinit"], [1.0] * 600
            )
        bckgnoise = bckgn.copy()
        for j in range(old_div(Tracker["nxinit"], 2) + 1, Tracker["nxinit"]):
            bckgn[0, j] = bckgn[0, Tracker["nxinit"] - j]

        if Tracker["constants"]["CTF"]:
            ctf_params = dataimage.get_attr("ctf")
            ctf_params.apix = old_div(
                ctf_params.apix,
                (
                    old_div(
                        float(Tracker["nxinit"]), float(Tracker["constants"]["nnxo"])
                    )
                ),
            )
            ctfa = sp_morphology.ctf_img_real(Tracker["nxinit"], ctf_params)
            ctfs = ctfa
        dataml = sp_fundamentals.fdecimate(
            dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False, False
        )
        data = []
        for iq in coarse_shifts:
            xx = iq[0] * shrink
            yy = iq[1] * shrink
            dss = sp_fundamentals.fshift(dataml, xx, yy)
            dss.set_attr("is_complex", 0)
            data.append(dss)

        #  This will get it to real space
        # dataimage = fpol(Util.mulnclreal(Util.mulnclreal(fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False), Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)
        #  Here we can reuse dataml to save time.  It was not normalized, but it does not matter as dataimage is for POLAR.  PAP 08/06/2018
        dataimage = sp_fundamentals.fpol(
            EMAN2_cppwrap.Util.mulnclreal(
                EMAN2_cppwrap.Util.mulnclreal(
                    dataml, EMAN2_cppwrap.Util.muln_img(bckgn, ctfs)
                ),
                mask,
            ),
            Tracker["nxpolar"],
            Tracker["nxpolar"],
            1,
            True,
        )

        if (
            (Blockdata["myid"] == Blockdata["main_node"])
            and (im % (max(1, old_div(nima, 5))) == 0)
            and (im > 0)
        ):
            sp_global_def.sxprint(
                "  Number of images :%7d   %5d  %5.1f"
                % (im, nima, old_div(float(im), float(nima)) * 100.0)
                + "%"
                + "   %10.1fmin" % (old_div((time.time() - at), 60.0))
            )

        if im < min(nima, 1) and procid == 0:
            #   Search all and compare with direct to figure what keepfirst might be
            keepfirst = old_div(
                (n_coarse_ang * n_coarse_psi), 10
            )  # keepfirst = (n_coarse_ang *  n_coarse_psi * n_coarse_shifts)/10

            xod2 = sp_helix_sphire.np.asarray(
                EMAN2_cppwrap.Util.multiref_Crosrng_msg_stack_stepsi(
                    dataimage,
                    bigbuffer,
                    coarse_shifts_shrank,
                    numr,
                    [coarse_angles[ic][2] for ic in range(n_coarse_ang)],
                    coarse_delta,
                    cnx,
                    keepfirst,
                )
            )

            assert len(xod2) == keepfirst

            xod1 = sp_helix_sphire.np.ndarray((keepfirst), dtype="f4", order="C")

            for iln in range(keepfirst):
                m = xod2[iln]
                j = m % n_coarse_psi
                ic = (old_div(m, n_coarse_psi)) % n_coarse_ang
                ib = old_div(m, (n_coarse_ang * n_coarse_psi))
                xod2[iln] = j * 1000 + ic * 100000000 + ib  # hashparams
            # DO NOT order by angular directions to save time on reprojections.
            pre_ipsiandiang = -1
            for iln in range(keepfirst):
                hashparams = int(xod2[iln])
                ishift = hashparams % 1000
                ipsiandiang = old_div(hashparams, 1000)
                if ipsiandiang != pre_ipsiandiang:
                    pre_ipsiandiang = ipsiandiang
                    ipsi = ipsiandiang % 100000
                    iang = old_div(ipsiandiang, 100000)
                    temp = sp_projection.prgl(
                        volinit,
                        [
                            coarse_angles[iang][0],
                            coarse_angles[iang][1],
                            (coarse_angles[iang][2] + ipsi * coarse_delta) % 360.0,
                            0.0,
                            0.0,
                        ],
                        1,
                        False,
                    )
                    temp.set_attr("is_complex", 0)
                xod1[iln] = -EMAN2_cppwrap.Util.sqed(
                    data[ishift], temp, ctfa, bckgnoise
                )  # peak
                ##xod2[iln] = hashparams

            xod1 -= sp_helix_sphire.np.max(xod1)
            lina = sp_helix_sphire.np.argwhere(
                xod1 > Tracker["constants"]["expthreshold"]
            )
            # if( Blockdata["myid"] == 0 ):
            # print("  STARTING ",Blockdata["myid"],np.max(xod1),np.min(xod1),len(lina),lina[-1])
            lina = lina.reshape(lina.size)
            keepf = int(lina[-1]) + 1

            xod1 = xod1[lina]
            xod2 = xod2[lina]

            lina = sp_helix_sphire.np.argsort(xod1)
            xod1 = xod1[lina[::-1]]  # This sorts in reverse order
            xod2 = xod2[lina[::-1]]  # This sorts in reverse order
            sp_helix_sphire.np.exp(xod1, out=xod1)
            xod1 = old_div(xod1, sp_helix_sphire.np.sum(xod1))
            cumprob = 0.0
            lit = len(xod1)
            for j in range(len(xod1)):
                cumprob += xod1[j]
                if cumprob > Tracker["constants"]["ccfpercentage"]:
                    lit = j + 1
                    break
            # keepf = mpi_reduce(keepf, 1, MPI_INT, MPI_MAX, Blockdata["main_node"], MPI_COMM_WORLD)
            # if( Blockdata["myid"] == 0 ):
            # 	keepf = max(int(float(keepf)*0.9),1)
            keepf = [keepf]
            keepf = sp_utilities.wrap_mpi_gatherv(
                keepf, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            if Blockdata["myid"] == 0:
                keepf.sort()
                keepf = keepf[int(len(keepf) * Blockdata["rkeepf"])]
            else:
                keepf = 0
            keepf = sp_utilities.wrap_mpi_bcast(
                keepf, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )
            Tracker["keepfirst"] = int(keepf)
            ###if( Blockdata["myid"] == 0 ):  sxprint("  keepfirst first ",Tracker["keepfirst"])

        else:
            # Tracker["keepfirst"] = min(200,nang)#min(max(lxod1/25,200),lxod1)

            xod2 = sp_helix_sphire.np.asarray(
                EMAN2_cppwrap.Util.multiref_Crosrng_msg_stack_stepsi(
                    dataimage,
                    bigbuffer,
                    coarse_shifts_shrank,
                    numr,
                    [coarse_angles[ic][2] for ic in range(n_coarse_ang)],
                    coarse_delta,
                    cnx,
                    Tracker["keepfirst"],
                )
            )

            xod1 = sp_helix_sphire.np.ndarray(
                (Tracker["keepfirst"]), dtype="f4", order="C"
            )

            for iln in range(Tracker["keepfirst"]):
                m = xod2[iln]
                j = m % n_coarse_psi
                ic = (old_div(m, n_coarse_psi)) % n_coarse_ang
                ib = old_div(m, (n_coarse_ang * n_coarse_psi))
                xod2[iln] = j * 1000 + ic * 100000000 + ib  # hashparams
            # order by angular directions to save time on reprojections.
            ipsiandiang = old_div(xod2, 1000)
            lina = sp_helix_sphire.np.argsort(ipsiandiang)
            xod2 = xod2[lina]  # order does not matter
            pre_ipsiandiang = -1
            for iln in range(Tracker["keepfirst"]):
                hashparams = int(xod2[iln])
                ishift = hashparams % 1000
                ipsiandiang = old_div(hashparams, 1000)
                if ipsiandiang != pre_ipsiandiang:
                    pre_ipsiandiang = ipsiandiang
                    ipsi = ipsiandiang % 100000
                    iang = old_div(ipsiandiang, 100000)
                    temp = sp_projection.prgl(
                        volinit,
                        [
                            coarse_angles[iang][0],
                            coarse_angles[iang][1],
                            (coarse_angles[iang][2] + ipsi * coarse_delta) % 360.0,
                            0.0,
                            0.0,
                        ],
                        1,
                        False,
                    )
                    temp.set_attr("is_complex", 0)
                peak = -EMAN2_cppwrap.Util.sqed(data[ishift], temp, ctfa, bckgnoise)
                xod1[iln] = peak
                ##xod2[iln] = hashparams

            lina = sp_helix_sphire.np.argsort(xod1)
            xod1 = xod1[lina[::-1]]  # This sorts in reverse order
            xod2 = xod2[lina[::-1]]  # This sorts in reverse order

            xod1 -= xod1[0]

            lina = sp_helix_sphire.np.argwhere(
                xod1 > Tracker["constants"]["expthreshold"]
            )
            xod1 = xod1[lina]
            xod2 = xod2[lina]
            sp_helix_sphire.np.exp(xod1, out=xod1)
            xod1 = old_div(xod1, sp_helix_sphire.np.sum(xod1))
            cumprob = 0.0
            lit = len(xod1)
            for j in range(len(xod1)):
                cumprob += xod1[j]
                if cumprob > Tracker["constants"]["ccfpercentage"]:
                    lit = j + 1
                    break

        del data

        # if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  SECOND KEPT  ",lit)

        firstdirections = [[0.0, 0.0] for iln in range(lit)]
        firstshifts = [0] * lit
        for iln in range(lit):
            hashparams = int(xod2[iln])
            ishift = hashparams % 1000
            ipsiandiang = old_div(hashparams, 1000)
            # ipsi = ipsiandiang%100000
            iang = old_div(ipsiandiang, 100000)
            firstdirections[iln] = [coarse_angles[iang][0], coarse_angles[iang][1], 0.0]
            firstshifts[iln] = ishift
        ###del xod2
        ###if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  SECOND KEPT  ",im,lit)#,len(xod2),xod2,firstdirections,firstshifts)
        # mpi_barrier(MPI_COMM_WORLD)
        # mpi_finalize()
        # exit()

        # Find neighbors, ltabang contains positions on refang list, no psis
        ltabang = find_nearest_k_refangles_to_many_angles(
            refdirs, firstdirections, Tracker["delta"], howmany=Tracker["howmany"]
        )

        # ltabang has length lit, which is the same as length as firshifts.  However, original xod2 is still available,
        #   even though it is longer than lit.

        #  Prepare image for chi2.
        #  We have to repeat everything from get shrink data, including shifts
        #  The number of entries is lit=len(ltabang), each proj dir has 4 neighbors
        #    there are 2 psis, and at most n_fine_shifts. which should be 4.
        #    Not not all shifts have to be populated. We may also have repeated triplets of angles, but each can have
        #      different shifts.  If so, we have to remove duplicates from the entire set.
        lcod1 = lit * 4 * 2 * n_fine_shifts
        cod2 = []
        # lol = 0
        for i1 in range(lit):
            hashparams = int(xod2[i1])
            ipsiandiang = old_div(hashparams, 1000)
            oldiang = old_div(ipsiandiang, 100000)
            ipsi = ipsiandiang % 100000
            ishift = hashparams % 1000
            tshifts = get_shifts_neighbors(shifts, coarse_shifts[ishift])
            for i2 in range(Tracker["howmany"]):
                iang = ltabang[i1][i2]
                for i3 in range(2):  # psi
                    itpsi = int(
                        old_div(
                            (
                                coarse_angles[oldiang][2]
                                + ipsi * coarse_delta
                                - refang[iang][2]
                                + 360.0
                            ),
                            Tracker["delta"],
                        )
                    )
                    itpsi = (itpsi + i3) % npsi
                    for i4 in range(len(tshifts)):
                        cod2.append(iang * 100000000 + itpsi * 1000 + tshifts[i4])

        del xod1, xod2

        # if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  THIRD1   ",len(cod2),cod2)
        cod2 = list(set(cod2))
        cod1 = [[old_div(q, 1000), i] for i, q in enumerate(cod2)]
        cod1.sort()

        lit = len(cod1)

        cod2 = sp_helix_sphire.np.asarray([cod2[cod1[i][1]] for i in range(lit)])

        cod1 = sp_helix_sphire.np.ndarray(lit, dtype="f4", order="C")
        # cod1.fill(np.finfo(dtype='f4').min)
        cod3 = sp_helix_sphire.np.ndarray(lit, dtype="f4", order="C")
        # cod3.fill(0.0)  #  varadj

        # if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  THIRD2   ",im,lit,cod2)

        # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        #  Third ML part.  Now the actual chi2 distances, we compute reprojections on the fly.

        data = [None] * nshifts
        johi = 0
        iln = 0
        prevdir = -1
        while iln < lit:
            hashparams = cod2[iln]
            ipsiandiang = old_div(hashparams, 1000)
            if ipsiandiang != prevdir:
                prevdir = ipsiandiang
                ipsi = ipsiandiang % 100000
                iang = old_div(ipsiandiang, 100000)
                temp = sp_projection.prgl(
                    volinit,
                    [
                        refang[iang][0],
                        refang[iang][1],
                        (refang[iang][2] + ipsi * Tracker["delta"]) % 360.0,
                        0.0,
                        0.0,
                    ],
                    1,
                    False,
                )
                temp.set_attr("is_complex", 0)
                johi += 1
            while ipsiandiang == old_div(cod2[iln], 1000):
                hashparams = cod2[iln]
                ishift = hashparams % 1000
                if data[ishift] == None:
                    xx = shifts[ishift][0] * shrink
                    yy = shifts[ishift][1] * shrink
                    data[ishift] = sp_fundamentals.fshift(dataml, xx, yy)
                    data[ishift].set_attr("is_complex", 0)

                [peak, varadj] = EMAN2_cppwrap.Util.sqednorm(
                    data[ishift], temp, ctfa, bckgnoise
                )
                cod1[iln] = -peak
                cod3[iln] = varadj
                iln += 1
                if iln == lit:
                    break

        del data
        del dataml

        lina = sp_helix_sphire.np.argsort(cod1)
        cod1 = cod1[lina[::-1]]  # This sorts in reverse order
        cod2 = cod2[lina[::-1]]  # This sorts in reverse order
        cod3 = cod3[lina[::-1]]  # This sorts in reverse order
        cod1 -= cod1[0]
        lina = sp_helix_sphire.np.argwhere(cod1 > Tracker["constants"]["expthreshold"])
        cod1 = cod1[lina]
        cod2 = cod2[lina]
        cod3 = cod3[lina]

        sp_helix_sphire.np.exp(cod1, out=cod1)
        cod1 = old_div(cod1, sp_helix_sphire.np.sum(cod1))
        cumprob = 0.0
        for j in range(len(cod1)):
            cumprob += cod1[j]
            if cumprob > Tracker["constants"]["ccfpercentage"]:
                lit = j + 1
                break

        #  New norm is a sum of eq distances multiplied by their probabilities augmented by PW.
        norm_per_particle[im] = (
            sp_helix_sphire.np.sum(cod1[:lit] * cod3[:lit]) + accumulatepw[im][reachpw]
        )

        for iln in range(lit):
            newpar[im][2].append([int(cod2[iln]), float(cod1[iln])])

        del cod1, cod2, cod3, lina
        ###mpi_barrier(MPI_COMM_WORLD)
        ###mpi_finalize()
        ###exit()
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_global_def.sxprint(
            "  Finished projection matching   %10.1fmin"
            % (old_div((time.time() - at), 60.0))
        )
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    #  All images were processed, now to the additional calculations
    ###mpi_barrier(MPI_COMM_WORLD)
    ###mpi_finalize()
    ###exit()

    # norm correction ---- calc the norm correction per particle
    snormcorr = 0.0
    for kl in range(nima):
        norm_per_particle[kl] = old_div(
            numpy.sqrt(norm_per_particle[kl] * 2.0) * oldparams[kl][7],
            Tracker["avgvaradj"][procid],
        )
        snormcorr += norm_per_particle[kl]
    Tracker["avgvaradj"][procid] = snormcorr
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    #  Compute avgvaradj
    Tracker["avgvaradj"][procid] = mpi.mpi_reduce(
        Tracker["avgvaradj"][procid],
        1,
        mpi.MPI_FLOAT,
        mpi.MPI_SUM,
        Blockdata["main_node"],
        mpi.MPI_COMM_WORLD,
    )
    if Blockdata["myid"] == Blockdata["main_node"]:
        Tracker["avgvaradj"][procid] = old_div(
            float(Tracker["avgvaradj"][procid]), Tracker["nima_per_chunk"][procid]
        )
    else:
        Tracker["avgvaradj"][procid] = 0.0
    Tracker["avgvaradj"][procid] = sp_utilities.bcast_number_to_all(
        Tracker["avgvaradj"][procid], Blockdata["main_node"]
    )
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    #  Compute statistics of smear -----------------
    smax = -1000000
    smin = 1000000
    sava = 0.0
    svar = 0.0
    snum = 0
    for kl in range(nima):
        j = len(newpar[kl][2])
        snum += 1
        sava += float(j)
        svar += j * float(j)
        smax = max(smax, j)
        smin = min(smin, j)
    snum = mpi.mpi_reduce(
        snum, 1, mpi.MPI_INT, mpi.MPI_SUM, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    sava = mpi.mpi_reduce(
        sava, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    svar = mpi.mpi_reduce(
        svar, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    smax = mpi.mpi_reduce(
        smax, 1, mpi.MPI_INT, mpi.MPI_MAX, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    smin = mpi.mpi_reduce(
        smin, 1, mpi.MPI_INT, mpi.MPI_MIN, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    if Blockdata["myid"] == 0:
        sava = old_div(float(sava), snum)
        svar = numpy.sqrt(
            max(0.0, old_div((float(svar) - snum * sava ** 2), (snum - 1)))
        )
        # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
        line = ""
        sp_global_def.sxprint(
            line,
            "Smear stat  (number of images, ave, sumsq, min, max)):  %7d    %12.3g   %12.3g  %7d  %7d"
            % (snum, sava, svar, smin, smax),
        )

    at = time.time()
    mpi.mpi_barrier(Blockdata["shared_comm"])

    ###if Blockdata["myid"] == Blockdata["main_node"]:  sxprint "  Finished :",time()-at
    mpi.mpi_win_free(win_sm)
    if Tracker["nxinit"] != Tracker["nxpolar"]:
        mpi.mpi_win_free(win_volinit)
        emnumpy4.unregister_numpy_from_emdata()
        del emnumpy4
    else:
        mpi.mpi_win_free(win_vol)

    mpi.mpi_barrier(Blockdata["shared_comm"])
    emnumpy1.unregister_numpy_from_emdata()
    emnumpy2.unregister_numpy_from_emdata()
    del emnumpy1, emnumpy2

    del volinit

    # print("  NORMALIZATION DONE  ",Blockdata["myid"])
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    if Blockdata["myid"] == Blockdata["main_node"]:
        # write_text_row([[newpar[0][2][j][0],newpar[0][2][j][1]] for j in range(len(newpar[0][2]))],os.path.join(Tracker["directory"], "polar%1d.txt"%procid))
        sp_global_def.sxprint(
            "  Statistics finished : %10.1fmin" % (old_div((time.time() - at), 60.0))
        )
    return newpar, norm_per_particle


# MODIFIED FROM TRUE POLAR  12/06/2017  PAP


def ali3D_primary_local_polar(
    refang,
    shifts,
    coarse_angles,
    coarse_shifts,
    procid,
    original_data=None,
    oldparams=None,
    preshift=False,
    apply_mask=True,
    nonorm=False,
    applyctf=True,
):
    """
	12/06/2017
	"""
    global Tracker, Blockdata
    #  Input data has to be CTF-multiplied, preshifted
    #  Output - newpar, see structure
    #    newpar = [[i, [worst_similarity, sum_all_similarities], [[-1, -1.0e23] for j in range(Tracker["lentop"])]] for i in range(len(data))]
    #    newpar = [[params],[],... len(data)]
    #    params = [particleID, [worst_similarity, sum_all_similarities],[imageallparams]]]
    #    imageallparams = [[orientation, similarity],[],...  number of all orientations ]
    #  Coding of orientations:
    #    hash = ang*100000000 + lpsi*1000 + ishift
    #    ishift = hash%1000
    #    ipsi = (hash/1000)%100000
    #    iang  = hash/100000000
    #  To get best matching for particle #kl:
    #     hash_best = newpar[kl][-1][0][0]
    #     best_sim  = newpar[kl][-1][0][1]
    #  To sort:
    # from operator 		import itemgetter, attrgetter, methodcaller
    #   params.sort(key=itemgetter(2))
    # from fundamentals import resample

    if Blockdata["myid"] == Blockdata["main_node"]:
        print_dict(
            Tracker, "PROJECTION MATCHING parameters of buffered primary local polar"
        )

    at = time.time()
    shrinkage = old_div(float(Tracker["nxpolar"]), float(Tracker["constants"]["nnxo"]))
    shrink = old_div(float(Tracker["nxinit"]), float(Tracker["constants"]["nnxo"]))
    radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
    mode = "F"
    numr = Numrinit_local(1, radius, 1, mode)
    wr = sp_alignment.ringwe(numr, mode)
    cnx = float(old_div(Tracker["nxpolar"], 2) + 1)

    ##if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  RADIUS IN POLAR  ",radius,numr[-1])
    #  FINE SEARCH CONSTANTS
    nang = len(refang)
    ac_fine = numpy.cos(numpy.radians(Tracker["an"]))
    npsi = int(old_div(360.0, Tracker["delta"]))
    mpsi = 2
    c_fine_psi = old_div(mpsi, 2)
    nshifts = len(shifts)
    n_fine_shifts = 4

    nima = len(original_data)
    mask = EMAN2_cppwrap.Util.unrollmask(Tracker["nxinit"], Tracker["nxinit"])
    for j in range(old_div(Tracker["nxinit"], 2), Tracker["nxinit"]):
        mask[0, j] = 1.0
    mask2D = sp_utilities.model_circle(
        Tracker["constants"]["radius"],
        Tracker["constants"]["nnxo"],
        Tracker["constants"]["nnxo"],
    )
    reachpw = (
        mask.get_xsize()
    )  # The last element of accumulated pw is zero so for the full size nothing is added.

    #  COARSE SEARCH CONSTANTS
    n_coarse_ang = len(coarse_angles)
    coarse_delta = 2 * Tracker["delta"]
    n_coarse_psi = int(old_div(360.0, coarse_delta))
    ###m_coarse_psi = int((2*2*Tracker["an"])/coarse_delta + 0.5) + 1
    m_coarse_psi = int(old_div(Tracker["an"], Tracker["delta"]) + 0.5)
    c_coarse_psi = old_div(m_coarse_psi, 2)
    n_coarse_shifts = len(coarse_shifts)

    coarse_shifts_shrank = [None] * n_coarse_shifts
    for ib in range(n_coarse_shifts):
        coarse_shifts_shrank[ib] = [
            coarse_shifts[ib][0] * shrinkage,
            coarse_shifts[ib][1] * shrinkage,
        ]

    ny = Tracker["nxinit"]
    nyp2 = old_div(ny, 2)
    nxth = old_div((Tracker["nxinit"] + 2), 2)
    indx = sp_utilities.model_blank(nxth, Tracker["nxinit"], 1, -1)
    tfrac = sp_utilities.model_blank(nxth, Tracker["nxinit"])
    tcount = sp_utilities.model_blank(nxth)
    for iy in range(1, ny + 1):
        jy = iy - 1
        if jy > nyp2:
            jy = jy - ny
        argy = float(jy * jy)
        for ix in range(1, nxth + 1):
            jx = ix - 1
            roff = jx + (iy - 1) * nxth
            if mask[ix - 1, iy - 1] > 0.0:
                rf = numpy.sqrt(argy + float(jx * jx))
                ir = int(rf)
                # print  ix-1,iy-1,roff,mask[ix-1,iy-1],rf,ir

                if ir < nxth - 1:
                    frac = rf - float(ir)
                    qres = 1.0 - frac
                    tfrac[ix - 1, iy - 1] = frac
                    # ioff = 2*roff
                    tcount[ir] += qres
                    tcount[ir + 1] += frac
                    indx[ix - 1, iy - 1] = ir

    disp_unit = sp_helix_sphire.np.dtype("f4").itemsize

    #  REFVOL
    if Blockdata["myid_on_node"] == 0:
        if Blockdata["myid"] == Blockdata["main_node"]:
            odo = get_refvol(Tracker["nxpolar"])
            nxvol = odo.get_xsize()
            nyvol = odo.get_ysize()
            nzvol = odo.get_zsize()
        else:
            nxvol = 0
            nyvol = 0
            nzvol = 0

        nxvol = sp_utilities.bcast_number_to_all(
            nxvol,
            source_node=Blockdata["main_node"],
            mpi_comm=Blockdata["group_zero_comm"],
        )
        nyvol = sp_utilities.bcast_number_to_all(
            nyvol,
            source_node=Blockdata["main_node"],
            mpi_comm=Blockdata["group_zero_comm"],
        )
        nzvol = sp_utilities.bcast_number_to_all(
            nzvol,
            source_node=Blockdata["main_node"],
            mpi_comm=Blockdata["group_zero_comm"],
        )

        if Blockdata["myid"] != Blockdata["main_node"]:
            odo = sp_utilities.model_blank(nxvol, nyvol, nzvol)

        sp_utilities.bcast_EMData_to_all(
            odo,
            Blockdata["group_zero_myid"],
            source_node=Blockdata["main_node"],
            comm=Blockdata["group_zero_comm"],
        )

        odo = sp_projection.prep_vol(odo, npad=2, interpolation_method=1)
        ndo = EMAN2_cppwrap.EMNumPy.em2numpy(odo)
        nxvol = odo.get_xsize()
        nyvol = odo.get_ysize()
        nzvol = odo.get_zsize()
        orgsizevol = nxvol * nyvol * nzvol
        sizevol = orgsizevol
    else:
        orgsizevol = 0
        sizevol = 0
        nxvol = 0
        nyvol = 0
        nzvol = 0

    orgsizevol = sp_utilities.bcast_number_to_all(
        orgsizevol, source_node=Blockdata["main_node"]
    )
    nxvol = sp_utilities.bcast_number_to_all(nxvol, source_node=Blockdata["main_node"])
    nyvol = sp_utilities.bcast_number_to_all(nyvol, source_node=Blockdata["main_node"])
    nzvol = sp_utilities.bcast_number_to_all(nzvol, source_node=Blockdata["main_node"])

    win_vol, base_vol = mpi.mpi_win_allocate_shared(
        sizevol * disp_unit, disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"]
    )
    sizevol = orgsizevol
    if Blockdata["myid_on_node"] != 0:
        base_vol, = mpi.mpi_win_shared_query(win_vol, mpi.MPI_PROC_NULL)
    """
    Comments from Adnan:
    numpy.core.mutliarray.int_asbuffer function is not available anymore . That is why a different alternative is tried here.
    """
    # volbuf = sp_helix_sphire.np.frombuffer(
    #     sp_helix_sphire.np.core.multiarray.int_asbuffer(base_vol, sizevol * disp_unit),
    #     dtype="f4",
    # )
    ptr = ctypes.cast(base_vol, ctypes.POINTER(ctypes.c_int * sizevol))
    volbuf = numpy.frombuffer(ptr.contents, dtype="f4")
    volbuf = volbuf.reshape(nzvol, nyvol, nxvol)
    if Blockdata["myid_on_node"] == 0:
        sp_helix_sphire.np.copyto(volbuf, ndo)
        del odo, ndo

    # volprep = EMNumPy.assign_numpy_to_emdata(volbuf)
    emnumpy1 = EMAN2_cppwrap.EMNumPy()
    volprep = emnumpy1.register_numpy_to_emdata(volbuf)
    volprep.set_attr_dict(
        {
            "is_complex": 1,
            "is_complex_x": 0,
            "is_fftodd": 0,
            "is_fftpad": 1,
            "is_shuffled": 1,
            "npad": 2,
        }
    )

    if Tracker["nxinit"] != Tracker["nxpolar"]:
        mpi.mpi_win_free(win_vol)
        #  REFVOL FOR ML
        if Blockdata["myid_on_node"] == 0:
            if Blockdata["myid"] == Blockdata["main_node"]:
                odo = get_refvol(Tracker["nxpolar"])
                nxvol = odo.get_xsize()
                nyvol = odo.get_ysize()
                nzvol = odo.get_zsize()
            else:
                nxvol = 0
                nyvol = 0
                nzvol = 0

            nxvol = sp_utilities.bcast_number_to_all(
                nxvol,
                source_node=Blockdata["main_node"],
                mpi_comm=Blockdata["group_zero_comm"],
            )
            nyvol = sp_utilities.bcast_number_to_all(
                nyvol,
                source_node=Blockdata["main_node"],
                mpi_comm=Blockdata["group_zero_comm"],
            )
            nzvol = sp_utilities.bcast_number_to_all(
                nzvol,
                source_node=Blockdata["main_node"],
                mpi_comm=Blockdata["group_zero_comm"],
            )

            if Blockdata["myid"] != Blockdata["main_node"]:
                odo = sp_utilities.model_blank(nxvol, nyvol, nzvol)

            sp_utilities.bcast_EMData_to_all(
                odo,
                Blockdata["group_zero_myid"],
                source_node=Blockdata["main_node"],
                comm=Blockdata["group_zero_comm"],
            )

            odo = sp_projection.prep_vol(odo, npad=2, interpolation_method=1)
            ndo = EMAN2_cppwrap.EMNumPy.em2numpy(odo)
            nxvol = odo.get_xsize()
            nyvol = odo.get_ysize()
            nzvol = odo.get_zsize()
            orgsizevol = nxvol * nyvol * nzvol
            sizevol = orgsizevol
        else:
            orgsizevol = 0
            sizevol = 0
            nxvol = 0
            nyvol = 0
            nzvol = 0

        orgsizevol = sp_utilities.bcast_number_to_all(
            orgsizevol, source_node=Blockdata["main_node"]
        )
        nxvol = sp_utilities.bcast_number_to_all(
            nxvol, source_node=Blockdata["main_node"]
        )
        nyvol = sp_utilities.bcast_number_to_all(
            nyvol, source_node=Blockdata["main_node"]
        )
        nzvol = sp_utilities.bcast_number_to_all(
            nzvol, source_node=Blockdata["main_node"]
        )

        win_volinit, base_volinit = mpi.mpi_win_allocate_shared(
            sizevol * disp_unit, disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"]
        )
        sizevol = orgsizevol
        if Blockdata["myid_on_node"] != 0:
            base_volinit, = mpi.mpi_win_shared_query(win_volinit, mpi.MPI_PROC_NULL)
        """
        Comments from Adnan:
        numpy.core.mutliarray.int_asbuffer function is not available anymore . That is why a different alternative is tried here.
        """
        # volbufinit = sp_helix_sphire.np.frombuffer(
        #     sp_helix_sphire.np.core.multiarray.int_asbuffer(
        #         base_volinit, sizevol * disp_unit
        #     ),
        #     dtype="f4",
        # )
        ptr = ctypes.cast(base_volinit, ctypes.POINTER(ctypes.c_int * sizevol))
        volbufinit = numpy.frombuffer(ptr.contents, dtype="f4")
        volbufinit = volbufinit.reshape(nzvol, nyvol, nxvol)
        if Blockdata["myid_on_node"] == 0:
            sp_helix_sphire.np.copyto(volbufinit, ndoinit)
            del odo, ndoinit

        # volinit = EMNumPy.assign_numpy_to_emdata(volbufinit)
        emnumpy4 = EMAN2_cppwrap.EMNumPy()
        volinit = emnumpy4.register_numpy_to_emdata(volbufinit)
        volinit.set_attr_dict(
            {
                "is_complex": 1,
                "is_complex_x": 0,
                "is_fftodd": 0,
                "is_fftpad": 1,
                "is_shuffled": 1,
                "npad": 2,
            }
        )
        if Blockdata["myid_on_node"] == 0:
            volinit.update()
        mpi.mpi_barrier(Blockdata["shared_comm"])
    else:
        volinit = volprep
    #  End of replaced volprep

    #  START CONES
    #  This has to be systematically done per node
    #
    crefim = EMAN2_cppwrap.Util.Polar2Dm(
        sp_utilities.model_blank(Tracker["nxpolar"], Tracker["nxpolar"]),
        cnx,
        cnx,
        numr,
        mode,
    )
    size_of_one_image = crefim.get_xsize()
    #  We will assume half of the memory is available.  We will do it betteer later.
    numberofrefs_inmem = int(
        old_div(
            old_div(Tracker["constants"]["memory_per_node"], 4),
            (old_div((size_of_one_image * disp_unit), 1.0e9)),
        )
    )
    ####if( Blockdata["myid_on_node"] == 0  ):  sxprint( " MEMEST ", n_coarse_ang,numberofrefs_inmem)
    #  number of references that will fit into one mode
    normals_set = sp_utilities.angles_to_normals(coarse_angles)
    Blockdata["angle_set"] = coarse_angles
    if n_coarse_ang <= numberofrefs_inmem:
        number_of_cones = 1
        numberofrefs_inmem = n_coarse_ang
        assignments_to_cones = [list(range(len(oldparams)))]

        assignments_of_refangles_to_cones = [
            [] for i in range(len(assignments_to_cones))
        ]
        assignments_of_refangles_to_angles = [
            [] for i in range(nima)
        ]  # for each myid separately, these are angles on this myid

        for i, q in enumerate(assignments_to_cones):
            #  find assignments of refdirs to this cone within an of each angledirs, so we need symmetry neighbors set of refdirs.
            ###print( "in loop ", Blockdata["myid"],i,len(q))#,q

            for m in q:
                # print " m ",m,len(angles)

                assignments_of_refangles_to_angles[
                    m
                ] = find_assignments_of_refangles_to_angles(
                    normals_set, oldparams[m], Tracker["an"]
                )
                assignments_of_refangles_to_cones[i].extend(
                    assignments_of_refangles_to_angles[m]
                )

            assignments_of_refangles_to_cones[i] = list(
                set(assignments_of_refangles_to_cones[i])
            )
            # print(  " assignments_of_refangles_to_cones on myid ",Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]) )
            assignments_of_refangles_to_cones[i] = sp_utilities.wrap_mpi_gatherv(
                assignments_of_refangles_to_cones[i], 0, Blockdata["shared_comm"]
            )
            doit = 1
            if Blockdata["myid_on_node"] == 0:
                assignments_of_refangles_to_cones[i] = list(
                    set(assignments_of_refangles_to_cones[i])
                )
                # print(  " assignments_of_refangles_to_cones on zero ",i,len(assignments_of_refangles_to_cones[i]) )#,assignments_of_refangles_to_cones[i]

        for i, q in enumerate(assignments_of_refangles_to_cones):
            ###if( Blockdata["myid_on_node"] == 0 ): sxprint( " which refangles belong to which cone",Blockdata["color"],Blockdata["myid"],i,len(q) )#,q
            assignments_of_refangles_to_cones[i] = sp_utilities.wrap_mpi_bcast(
                q, 0, Blockdata["shared_comm"]
            )

    else:

        angledirs = sp_utilities.angles_to_normals([u1[:3] for u1 in oldparams])

        number_of_cones = max(
            2, int(old_div(n_coarse_ang, numberofrefs_inmem) * 1.5 + 0.5)
        )
        ###if Blockdata["myid_on_node"] == 0:  sxprint( " LENGTHS  ",Blockdata["color"],nima,n_coarse_ang, numberofrefs_inmem,number_of_cones)
        cont = True
        while cont:
            #  Translate number of cones to cone_delta
            cone_delta, number_of_cones = number_of_cones_to_delta(number_of_cones)
            #  Generate cone_angles
            ##if Blockdata["myid"] == 0:  sxprint( "  WHILE  ",number_of_cones, cone_delta)
            if Blockdata["symclass"].sym[0] == "c":
                if number_of_cones == 1:
                    cone_delta = 180.1
                    cone_angles = [[0.0, 1.0, 0.0]]
                else:
                    cone_angles = Blockdata["symclass"].even_angles(
                        cone_delta, theta1=1.0, theta2=89.0
                    )
                    cone_angles += [
                        [(q[0] + 90.0) % 360.0, 180.0 - q[1], 0] for q in cone_angles
                    ]
                    cone_angles = Blockdata["symclass"].reduce_anglesets(cone_angles)
            else:
                cone_angles = Blockdata["symclass"].even_angles(cone_delta, theta1=1.0)

            # if Blockdata["myid"] == 0:  sxprint(  "  number of cones ",number_of_cones,cone_delta, len(cone_angles))
            assert number_of_cones == len(cone_angles)

            conedirs = sp_utilities.angles_to_normals(
                Blockdata["symclass"].symmetry_neighbors(cone_angles)
            )
            neighbors = old_div(
                len(conedirs), len(cone_angles)
            )  # Symmetry related neighbors
            # if Blockdata["myid"] == 0:  sxprint(  "  neighbors  ",Blockdata["myid"],neighbors, cone_angles)
            #  assign data directions to cone_angles
            assignments_to_cones = sp_utilities.assign_projdirs_f(
                angledirs, conedirs, neighbors
            )
            ###print(  " assignments_to_cones ",Blockdata["myid"],len(assignments_to_cones),[len(q) for q in assignments_to_cones],assignments_to_cones[0])
            #  the above should have length of refdirs and each should have indexes of data that belong to this cone
            del conedirs
            # print "assignments_to_cones ",assignments_to_cones
            #  For each cone we have to find which refangles are needed to do the matching
            assignments_of_refangles_to_cones = [
                [] for i in range(len(assignments_to_cones))
            ]
            assignments_of_refangles_to_angles = [
                [] for i in range(nima)
            ]  # for each myid separately, these are angles on this myid

            for i, q in enumerate(assignments_to_cones):
                #  find assignments of refdirs to this cone within an of each angledirs, so we need symmetry neighbors set of refdirs.
                # if Blockdata["myid"] == 0:  sxprint( "in loop ", Blockdata["myid"],i,len(q),q)

                if len(q) == 0:
                    # empty assignment, on a given CPU there are no images assigned to a given cone
                    assignments_of_refangles_to_cones[i] = [-1]
                else:
                    for m in q:
                        # print " m ",m,len(angles)

                        assignments_of_refangles_to_angles[
                            m
                        ] = find_assignments_of_refangles_to_angles(
                            normals_set, oldparams[m], Tracker["an"]
                        )
                        # if Blockdata["myid"] == 0:  sxprint( "assignments_of_refangles_to_angles[m] ", Blockdata["color"],i,m,assignments_of_refangles_to_angles[m])
                        assignments_of_refangles_to_cones[i].extend(
                            assignments_of_refangles_to_angles[m]
                        )

                    # if( Blockdata["myid"] == 19 ):  sxprint(  " doit0 ",Blockdata["myid"], i,assignments_of_refangles_to_cones[i],q,assignments_to_cones)
                    assignments_of_refangles_to_cones[i] = list(
                        set(assignments_of_refangles_to_cones[i])
                    )
                ###if Blockdata["myid"] == 19:  sxprint(  " assignments_of_refangles_to_cones on myid ",Blockdata["color"],Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]), assignments_of_refangles_to_cones[i] )
                assignments_of_refangles_to_cones[i] = sp_utilities.wrap_mpi_gatherv(
                    assignments_of_refangles_to_cones[i], 0, Blockdata["shared_comm"]
                )
                ###if Blockdata["myid"] == 0:  sxprint(  " assignments_of_refangles_to_cones gatherv",Blockdata["color"],Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]) )
                doit = 1
                if Blockdata["myid_on_node"] == 0:
                    assignments_of_refangles_to_cones[i] = list(
                        set(assignments_of_refangles_to_cones[i]) - set([-1])
                    )
                    ###if( Blockdata["myid_on_node"] == 0 ):  sxprint(  " assignments_of_refangles_to_cones on zero ",i,len(assignments_of_refangles_to_cones[i]))
                    #  POSSIBLE PROBLEM - IT IS POSSIBLE FOR A GIVEN CONE TO HAVE NO REFANGLES AND THUS BE EMPTY
                    if len(assignments_of_refangles_to_cones[i]) > numberofrefs_inmem:
                        number_of_cones = int(number_of_cones * 1.25)
                        # print(  " increased number_of_cones ",i,number_of_cones )
                        doit = 0
                doit = sp_utilities.bcast_number_to_all(doit, source_node=0)
                number_of_cones = sp_utilities.bcast_number_to_all(
                    number_of_cones, source_node=0, mpi_comm=Blockdata["shared_comm"]
                )
                ###if( Blockdata["myid"] == 19 ):  sxprint(  " doit ",Blockdata["myid"], i,doit ,assignments_of_refangles_to_cones[i],assignments_to_cones)
                if doit == 0:
                    break

            if doit == 1:
                cont = False

        for i, q in enumerate(assignments_of_refangles_to_cones):
            ###if( Blockdata["myid_on_node"] == 0 ): sxprint( " which refangles belong to which cone IOUT ",Blockdata["color"],Blockdata["myid"],i,len(q))#,q
            assignments_of_refangles_to_cones[i] = sp_utilities.wrap_mpi_bcast(
                q, 0, Blockdata["shared_comm"]
            )
            ###if( Blockdata["myid_on_node"] == 0 ): sxprint( " which refangles belong to which cone XOUT ",Blockdata["color"],Blockdata["myid"],i,len(q))#,q
            # if( myid == 1 ):
            # 	sxprint " which refangles belong to which cone",myid,i,len(assignments_of_refangles_to_cones[i])#,q

    if Blockdata["myid"] == 0:
        sp_global_def.sxprint(" number_of_cones: ", number_of_cones)
    #  Maximum number of refangles assigned to angles (max number of references per image)
    nlocal_angles = max([len(q) for q in assignments_of_refangles_to_angles])
    #  Most likely we have to delete some lists before proceeding
    del normals_set
    # We have to figure the maximum length of xod1, which is lang.  If there are no cones, it is an estimate.  If there are cones, we have list of assignments
    #  For number of cones I use refang and an.  This should give approximately the same number as coarse angles and 2*an, which is what is actually used in searches
    numberofrefs_inmem = max([len(q) for q in assignments_of_refangles_to_cones])

    ###for i,q in enumerate(assignments_of_refangles_to_cones):
    ###	sxprint( " which refangles belong to which cone OUT ",Blockdata["color"],Blockdata["myid_on_node"],Blockdata["myid"],i,numberofrefs_inmem,nlocal_angles,len(q))#,q

    #  BIG BUFFER
    lenbigbuf = numberofrefs_inmem  # MAXIMUM NUMBER OF REFERENCE IN ONE OF THE CONES
    orgsize = (
        lenbigbuf * size_of_one_image
    )  # This is number of projections to be computed simultaneously times their size

    if Blockdata["myid_on_node"] == 0:
        size = orgsize
    else:
        size = 0

    win_sm, base_ptr = mpi.mpi_win_allocate_shared(
        size * disp_unit, disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"]
    )
    size = orgsize
    if Blockdata["myid_on_node"] != 0:
        base_ptr, = mpi.mpi_win_shared_query(win_sm, mpi.MPI_PROC_NULL)
    """
    Comments from Adnan:
    numpy.core.mutliarray.int_asbuffer function is not available anymore . That is why a different alternative is tried here.
    """
    # buffer = sp_helix_sphire.np.frombuffer(
    #     sp_helix_sphire.np.core.multiarray.int_asbuffer(base_ptr, size * disp_unit),
    #     dtype="f4",
    # )
    ptr = ctypes.cast(base_ptr, ctypes.POINTER(ctypes.c_int * size))
    buffer = numpy.frombuffer(ptr.contents, dtype="f4")
    buffer = buffer.reshape(lenbigbuf, size_of_one_image)
    # bigbuffer = EMNumPy.assign_numpy_to_emdata(buffer)

    emnumpy2 = EMAN2_cppwrap.EMNumPy()
    bigbuffer = emnumpy2.register_numpy_to_emdata(buffer)

    #  end of CONES setup

    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_global_def.sxprint("  ")
        # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
        line = ""
        sp_global_def.sxprint(
            line,
            "Processing data for polar onx: %3d, nx: %3d, nxpolar: %3d, CTF: %s, preshift: %s,  MEM: %6.2fGB."
            % (
                Tracker["constants"]["nnxo"],
                Tracker["nxinit"],
                Tracker["nxpolar"],
                Tracker["constants"]["CTF"],
                preshift,
                old_div(orgsize, 1.0e9),
            ),
        )

    #  Note these are in Fortran notation for polar searches
    txm = float(Tracker["nxpolar"] - (old_div(Tracker["nxpolar"], 2) + 1) - radius)
    txl = float(radius - old_div(Tracker["nxpolar"], 2) + 1)

    if Blockdata["bckgnoise"]:
        oneover = []
        nxb = Blockdata["bckgnoise"][0].get_xsize()
        nyb = len(Blockdata["bckgnoise"])
        for i in range(nyb):
            temp = [0.0] * nxb
            for k in range(nxb):
                if Blockdata["bckgnoise"][i].get_value_at(k) > 0.0:
                    temp[k] = old_div(
                        1.0, numpy.sqrt(Blockdata["bckgnoise"][i].get_value_at(k))
                    )
            oneover.append(temp)
        del temp

        if procid == 0:
            Blockdata["totprob"] = [0.0] * nyb
            Blockdata["newbckgnoise"] = sp_utilities.model_blank(nxb, nyb)

    accumulatepw = [0.0] * nima
    norm_per_particle = [0.0] * nima

    ##lxod1 = lang*len(list_of_coarse_shifts)*(int(2*Tracker["an"]/coarse_delta+0.5)+1)
    newpar = [[i, [0.0], []] for i in range(nima)]

    #  This is for auxiliary function searches.
    Blockdata["angle_set"] = refang
    #  Extract normals from rotation matrices
    refdirs = sp_utilities.angles_to_normals(refang)

    #  We have to make sure the number of cones is the same on all nodes, this is due to the strange MPI problem
    #   that forced me to insert overall barrier into iterations over cones
    max_number_of_cones = number_of_cones
    max_number_of_cones = mpi.mpi_reduce(
        max_number_of_cones, 1, mpi.MPI_INT, mpi.MPI_MAX, 0, mpi.MPI_COMM_WORLD
    )
    if Blockdata["myid"] == Blockdata["main_node"]:
        max_number_of_cones = int(max_number_of_cones[0])
    max_number_of_cones = mpi.mpi_bcast(
        max_number_of_cones, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD
    )
    max_number_of_cones = int(max_number_of_cones[0])

    # if( Blockdata["myid_on_node"] == 0):
    # 	sxprint("  FIRST  nima, nang, npsi, nshift, n_coarse_ang, nlocal_angles, n_coarse_psi, len(list_of_coarse_shifts), number_of_cones , max_number_of_cones ",nima, nang,npsi,nshifts,n_coarse_ang,nlocal_angles,n_coarse_psi, len(coarse_shifts),number_of_cones,max_number_of_cones)

    ##firsti = True
    at = time.time()
    ##eat = 0.0
    lima = 0  # total counter of images
    #  PROCESSING OF CONES
    for icone in range(max_number_of_cones):
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        if (
            icone < number_of_cones
        ):  # This is executed for individual number of cones, some nodes may have fewer.
            nang_start, nang_end = sp_applications.MPI_start_end(
                len(assignments_of_refangles_to_cones[icone]),
                Blockdata["no_of_processes_per_group"],
                Blockdata["myid_on_node"],
            )
            # if(Blockdata["color"] == 1):
            ###print( " ZXZ11  ",Blockdata["color"],Blockdata["myid_on_node"],Blockdata["myid"],len(assignments_of_refangles_to_cones[icone]),nang_start, nang_end)

            for i in range(
                nang_start, nang_end, 1
            ):  # This will take care of no of process on a node less than nang.  Some loops will not be executed
                ic = assignments_of_refangles_to_cones[icone][i]
                temp = sp_projection.prgl(
                    volprep,
                    [coarse_angles[ic][0], coarse_angles[ic][1], 0.0, 0.0, 0.0],
                    1,
                    True,
                )
                crefim = EMAN2_cppwrap.Util.Polar2Dm(temp, cnx, cnx, numr, mode)
                EMAN2_cppwrap.Util.Frngs(crefim, numr)
                EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)
                bigbuffer.insert_clip(crefim, (0, i))

            mpi.mpi_barrier(Blockdata["shared_comm"])

            # if(Blockdata["myid"] == Blockdata["main_node"]):
            # 	sxprint( "  Reference projections generated for cone %4d: %10.1fmin"%(icone,(time()-at)/60.))
            ###print("  TOPOP    ",Blockdata["myid"],MPI_COMM_WORLD,len(assignments_to_cones[icone]))
            ###mpi_finalize()
            ###exit()
            ###  <><><><><><><><>

            #  Preprocess the data

            #  We only process images in the current cone.  icnm is consecutive number, im the actual image number
            # for icnm,im in enumerate(assignments_to_cones[icone]):
            lenass = len(assignments_to_cones[icone])
            ###print("   ENTERING  ",Blockdata["myid"],icone,lenass)
            for icnm in range(
                max(1, lenass)
            ):  # I have to enter the loop even it there is no assignment
                if lenass == 0:
                    keepf = -1
                    ###print("  FOUNDEMPTY  ",Blockdata["myid"],icone,icnm,len(assignments_to_cones[icone]),assignments_to_cones)
                else:
                    #  PREPARE ONE IMAGE
                    im = assignments_to_cones[icone][icnm]
                    particle_group = original_data[im].get_attr("particle_group")

                    phi, theta, psi, sx, sy, wnorm = (
                        oldparams[im][0],
                        oldparams[im][1],
                        oldparams[im][2],
                        oldparams[im][3],
                        oldparams[im][4],
                        oldparams[im][7],
                    )

                    #  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):  sxprint("\n\n  INPUT PARAMS  ",im,phi,theta,psi,sx,sy)
                    if preshift:
                        sx = int(round(sx))
                        sy = int(round(sy))
                        dataimage = sp_fundamentals.cyclic_shift(
                            original_data[im], sx, sy
                        )
                        #  Put rounded shifts on the list, note image has the original floats - check whether it may cause problems
                        oldparams[im][3] = sx
                        oldparams[im][4] = sy
                        sx = 0.0
                        sy = 0.0
                    else:
                        dataimage = original_data[im].copy()

                    st = get_image_statistics(dataimage, mask2D, False)
                    dataimage -= st[0]
                    dataimage = old_div(dataimage, st[1])
                    if dataimage.get_attr_default("bckgnoise", None):
                        dataimage.delete_attr("bckgnoise")
                    #  Do bckgnoise if exists
                    if Blockdata["bckgnoise"]:
                        if apply_mask:
                            if Tracker["constants"]["hardmask"]:
                                dataimage = sp_morphology.cosinemask(
                                    dataimage, radius=Tracker["constants"]["radius"]
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
                                    sp_filter.filt_table(bckg, oneover[particle_group])
                                )
                                #  Normalize bckg noise in real space, only region actually used.
                                st = get_image_statistics(bckg, mask2D, False)
                                bckg -= st[0]
                                bckg = old_div(bckg, st[1])
                                dataimage = sp_morphology.cosinemask(
                                    dataimage,
                                    radius=Tracker["constants"]["radius"],
                                    bckg=bckg,
                                )
                    else:
                        #  if no bckgnoise, do simple masking instead
                        if apply_mask:
                            dataimage = sp_morphology.cosinemask(
                                dataimage, radius=Tracker["constants"]["radius"]
                            )

                    #  Apply varadj
                    if not nonorm:
                        EMAN2_cppwrap.Util.mul_scalar(
                            dataimage, old_div(Tracker["avgvaradj"][procid], wnorm)
                        )

                    ###  FT
                    dataimage = sp_fundamentals.fft(dataimage)
                    sig = EMAN2_cppwrap.Util.rotavg_fourier(dataimage)
                    accumulatepw[im] = sig[old_div(len(sig), 2) :] + [0.0]

                    #  We have to make sure the shifts are within correct range, shrinkage or not
                    # set_params_proj(dataimage,[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
                    if Blockdata["bckgnoise"]:
                        temp = Blockdata["bckgnoise"][particle_group]
                        bckgn = EMAN2_cppwrap.Util.unroll1dpw(
                            Tracker["nxinit"],
                            Tracker["nxinit"],
                            [temp[i] for i in range(temp.get_xsize())],
                        )
                    else:
                        bckgn = EMAN2_cppwrap.Util.unroll1dpw(
                            Tracker["nxinit"], Tracker["nxinit"], [1.0] * 600
                        )
                    bckgnoise = bckgn.copy()
                    for j in range(
                        old_div(Tracker["nxinit"], 2) + 1, Tracker["nxinit"]
                    ):
                        bckgn[0, j] = bckgn[0, Tracker["nxinit"] - j]

                    if Tracker["constants"]["CTF"]:
                        ctf_params = dataimage.get_attr("ctf")
                        ctf_params.apix = old_div(
                            ctf_params.apix,
                            (
                                old_div(
                                    float(Tracker["nxinit"]),
                                    float(Tracker["constants"]["nnxo"]),
                                )
                            ),
                        )
                        ctfa = sp_morphology.ctf_img_real(Tracker["nxinit"], ctf_params)
                        ctfs = ctfa
                    ##if( ( Blockdata["myid"] == Blockdata["main_node"])   and firsti ):
                    ##	dataimage.set_attr("is_complex",0)
                    ##	dataimage.write_image("dataimagefft.hdf")
                    ##	dataimage.set_attr("is_complex",1)
                    dataml = sp_fundamentals.fdecimate(
                        dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False, False
                    )
                    data = []
                    for iq in coarse_shifts:
                        xx = iq[0] * shrink
                        yy = iq[1] * shrink
                        dss = sp_fundamentals.fshift(dataml, xx, yy)
                        dss.set_attr("is_complex", 0)
                        data.append(dss)

                    #  This will get it to real space
                    # dataimage = fpol(Util.mulnclreal(Util.mulnclreal(fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False), Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)
                    #  Here we can reuse dataml to save time.  It was not normalized, but it does not matter as dataimage is for POLAR.  PAP 08/06/2018
                    dataimage = sp_fundamentals.fpol(
                        EMAN2_cppwrap.Util.mulnclreal(
                            EMAN2_cppwrap.Util.mulnclreal(
                                dataml, EMAN2_cppwrap.Util.muln_img(bckgn, ctfs)
                            ),
                            mask,
                        ),
                        Tracker["nxpolar"],
                        Tracker["nxpolar"],
                        1,
                        True,
                    )

                    # Compute max number of angles on the fly
                    lang = len(assignments_of_refangles_to_angles[im])
                    ###print("   BICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)

                if (
                    (Blockdata["myid"] == Blockdata["main_node"])
                    and (lima % (max(1, old_div(nima, 5))) == 0)
                    and (lima > 0)
                ):
                    sp_global_def.sxprint(
                        "  Number of images :%7d   %5d  %5.1f"
                        % (lima, nima, old_div(float(lima), float(nima)) * 100.0)
                        + "%"
                        + "   %10.1fmin" % (old_div((time.time() - at), 60.0))
                    )
                    ##print( "  Number of images :%7d   %5d  %5.1f"%(lima,nima,float(lima)/float(nima)*100.) + "%" +"   %10.1fmin   %10.1fmin"%((time()-at)/60.,eat/60.0))
                lima += 1

                ###print("  CONA1    ",Blockdata["myid"],lima)
                if lima == 1 and procid == 0:
                    ###print("  CONA2    ",Blockdata["myid"])
                    if lenass > 0:
                        ###print("   CICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)
                        keepfirst = (
                            lang * m_coarse_psi * n_coarse_shifts
                        )  # 150#lang*m_coarse_psi*len(coarse_shifts_shrank)  #500#min(200,nlocal_angles)#min(max(lxod1/25,200),lxod1)
                        ###if( Blockdata["myid"] == 18 and lima<5):  sxprint(" START   nlocal_angles* m_coarse_psi*len(coarse_shifts_shrank",nlocal_angles,coarse_delta,len(coarse_shifts_shrank),keepfirst)
                        # if( Blockdata["myid"] == Blockdata["main_node"] ):
                        lxod1 = EMAN2_cppwrap.Util.multiref_Crosrng_msg_stack_stepsi_local(
                            dataimage,
                            bigbuffer,
                            coarse_shifts_shrank,
                            assignments_of_refangles_to_angles[im],
                            assignments_of_refangles_to_cones[icone],
                            numr,
                            [coarse_angles[i][2] for i in range(n_coarse_ang)],
                            oldparams[im][2],
                            c_coarse_psi,
                            coarse_delta,
                            cnx,
                            keepfirst,
                        )

                        ##'''
                        assert old_div(len(lxod1), 3) == keepfirst

                        xod1 = sp_helix_sphire.np.ndarray(
                            (keepfirst), dtype="f4", order="C"
                        )
                        # xod1.fill(1.0)
                        xod2 = sp_helix_sphire.np.ndarray(
                            (keepfirst), dtype="int", order="C"
                        )
                        for iq in range(keepfirst):
                            ioffset = 3 * iq
                            #          ishift         iang                      ipsi
                            xod2[iq] = (
                                lxod1[ioffset]
                                + lxod1[ioffset + 1] * 100000000
                                + lxod1[ioffset + 2] * 1000
                            )

                        ##'''
                        # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

                        #  Second step - find which coarse ones are significant

                        # DO NOT order by angular directions to save time on reprojections.

                        pre_ipsiandiang = -1
                        for iln in range(keepfirst):
                            hashparams = int(xod2[iln])
                            ishift = hashparams % 1000
                            ipsiandiang = old_div(hashparams, 1000)
                            if ipsiandiang != pre_ipsiandiang:
                                pre_ipsiandiang = ipsiandiang
                                ipsi = ipsiandiang % 100000
                                iang = old_div(ipsiandiang, 100000)
                                ##junk = time()
                                temp = sp_projection.prgl(
                                    volinit,
                                    [
                                        coarse_angles[iang][0],
                                        coarse_angles[iang][1],
                                        (coarse_angles[iang][2] + ipsi * coarse_delta)
                                        % 360.0,
                                        0.0,
                                        0.0,
                                    ],
                                    1,
                                    False,
                                )
                                ##eat += time()-junk
                                ###   if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  SELECTEDSTEPTWO  ",iln,refang[iang][0],refang[iang][1],refang[iang][2]+ipsi*Tracker["delta"],ishift,xod1[iln])
                                temp.set_attr("is_complex", 0)
                            ##junk = time()
                            xod1[iln] = -EMAN2_cppwrap.Util.sqed(
                                data[ishift], temp, ctfa, bckgnoise
                            )
                            ##eat += time()-junk
                            ##xod2[iln] = hashparams

                        xod1 -= sp_helix_sphire.np.max(xod1)
                        lina = sp_helix_sphire.np.argwhere(
                            xod1 > Tracker["constants"]["expthreshold"]
                        )
                        # if( Blockdata["myid"] == 0 ):
                        ###print("  STARTING1   ",Blockdata["myid"],np.max(xod1),np.min(xod1),len(lina),lina[-1])

                        lina = lina.reshape(lina.size)
                        keepf = int(lina[-1]) + 1

                        xod1 = xod1[lina]
                        xod2 = xod2[lina]

                        ###print("  STARTING2    ",Blockdata["myid"])

                        lina = sp_helix_sphire.np.argsort(xod1)
                        xod1 = xod1[lina[::-1]]  # This sorts in reverse order
                        xod2 = xod2[lina[::-1]]  # This sorts in reverse order
                        sp_helix_sphire.np.exp(xod1, out=xod1)
                        xod1 = old_div(xod1, sp_helix_sphire.np.sum(xod1))
                        cumprob = 0.0
                        lit = len(xod1)
                        ###print("  STARTING3    ",Blockdata["myid"],lit)
                        for j in range(len(xod1)):
                            cumprob += xod1[j]
                            if cumprob > Tracker["constants"]["ccfpercentage"]:
                                lit = j + 1
                                break

                        ###print("  STARTING4    ",Blockdata["myid"],lit,keepf)
                    #  Turn into percentage of all possible
                    keepf = [int(old_div(float(keepf * 100), float(keepfirst)))]
                    ###mpi_barrier(MPI_COMM_WORLD)
                    ###print("  STARTING5    ",Blockdata["myid"],keepf,MPI_COMM_WORLD)
                    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
                    keepf = sp_utilities.wrap_mpi_gatherv(
                        keepf, Blockdata["main_node"], mpi.MPI_COMM_WORLD
                    )
                    ###print("  STARTING6    ",Blockdata["myid"],keepf)
                    if Blockdata["myid"] == 0:
                        keepf = [junk for junk in keepf if junk > 0]
                        if len(keepf) < 2:
                            keepf = 3
                        else:
                            keepf.sort()
                            keepf = keepf[int(len(keepf) * Blockdata["rkeepf"])]
                    else:
                        keepf = 0
                    ###print("  STARTING7    ",Blockdata["myid"],keepf)
                    keepf = sp_utilities.wrap_mpi_bcast(
                        keepf, Blockdata["main_node"], mpi.MPI_COMM_WORLD
                    )
                    # if(keepf == 0):
                    # 	ERROR("Too few images to estimate keepfirst","sxmeridien", 1, Blockdata["myid"])
                    # 	mpi_finalize()
                    # 	exit()
                    ###print("  STARTING8    ",Blockdata["myid"],keepf)
                    Tracker["keepfirst"] = int(keepf)
                    ###if( Blockdata["myid"] == 0 ):  sxprint("  keepfirst first ",Tracker["keepfirst"])

                else:
                    if lenass > 0:
                        ###print("   DICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)
                        ###if( Blockdata["myid"] == 18 and lima<5):  sxprint(" START   nlocal_angles* m_coarse_psi*len(coarse_shifts_shrank",nlocal_angles,coarse_delta,len(coarse_shifts_shrank),Tracker["keepfirst"])
                        # if( Blockdata["myid"] == Blockdata["main_node"] ):
                        keepfirst = lang * m_coarse_psi * n_coarse_shifts
                        keepfirst = max(
                            old_div(keepfirst * Tracker["keepfirst"], 100),
                            min(keepfirst, 3),
                        )
                        lxod1 = EMAN2_cppwrap.Util.multiref_Crosrng_msg_stack_stepsi_local(
                            dataimage,
                            bigbuffer,
                            coarse_shifts_shrank,
                            assignments_of_refangles_to_angles[im],
                            assignments_of_refangles_to_cones[icone],
                            numr,
                            [coarse_angles[i][2] for i in range(n_coarse_ang)],
                            oldparams[im][2],
                            c_coarse_psi,
                            coarse_delta,
                            cnx,
                            keepfirst,
                        )
                        # if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  ")
                        ##'''
                        assert keepfirst == old_div(len(lxod1), 3)
                        xod1 = sp_helix_sphire.np.ndarray(
                            (keepfirst), dtype="f4", order="C"
                        )
                        # xod1.fill(1.0)
                        xod2 = sp_helix_sphire.np.ndarray(
                            (keepfirst), dtype="int", order="C"
                        )
                        for iq in range(keepfirst):
                            ioffset = 3 * iq
                            #          ishift         iang                      ipsi
                            xod2[iq] = (
                                lxod1[ioffset]
                                + lxod1[ioffset + 1] * 100000000
                                + lxod1[ioffset + 2] * 1000
                            )

                        ##'''

                        # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

                        #  Second step - find which coarse ones are significant

                        # order by angular directions to save time on reprojections.
                        ipsiandiang = old_div(xod2, 1000)
                        lina = sp_helix_sphire.np.argsort(ipsiandiang)
                        xod2 = xod2[lina]  # order does not matter

                        pre_ipsiandiang = -1
                        for iln in range(keepfirst):
                            hashparams = int(xod2[iln])
                            ishift = hashparams % 1000
                            ipsiandiang = old_div(hashparams, 1000)
                            if ipsiandiang != pre_ipsiandiang:
                                pre_ipsiandiang = ipsiandiang
                                ipsi = ipsiandiang % 100000
                                iang = old_div(ipsiandiang, 100000)
                                ##junk = time()
                                temp = sp_projection.prgl(
                                    volinit,
                                    [
                                        coarse_angles[iang][0],
                                        coarse_angles[iang][1],
                                        (coarse_angles[iang][2] + ipsi * coarse_delta)
                                        % 360.0,
                                        0.0,
                                        0.0,
                                    ],
                                    1,
                                    False,
                                )
                                ###   if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  SELECTEDSTEPTWO  ",iln,refang[iang][0],refang[iang][1],refang[iang][2]+ipsi*Tracker["delta"],ishift,xod1[iln])
                                ##eat += time()-junk
                                temp.set_attr("is_complex", 0)
                            ##junk = time()
                            peak = -EMAN2_cppwrap.Util.sqed(
                                data[ishift], temp, ctfa, bckgnoise
                            )
                            ##eat += time()-junk
                            #  Note I replace ccc by eqdist
                            xod1[iln] = peak
                            ##xod2[iln] = hashparams

                        lina = sp_helix_sphire.np.argsort(xod1)
                        xod1 = xod1[lina[::-1]]  # This sorts in reverse order
                        xod2 = xod2[lina[::-1]]  # This sorts in reverse order

                        xod1 -= xod1[0]

                        # if( Blockdata["myid"] == Blockdata["main_node"]):
                        # 	#print("  PROJECT   ",im,lit,johi)#,cod2)
                        # 	for iln in range(len(xod1)):  sxprint("  ROPE   ",iln,xod1[iln],xod2[iln])

                        lina = sp_helix_sphire.np.argwhere(
                            xod1 > Tracker["constants"]["expthreshold"]
                        )
                        xod1 = xod1[lina]
                        xod2 = xod2[lina]
                        sp_helix_sphire.np.exp(xod1, out=xod1)
                        xod1 = old_div(xod1, sp_helix_sphire.np.sum(xod1))
                        cumprob = 0.0
                        lit = len(xod1)
                        for j in range(len(xod1)):
                            cumprob += xod1[j]
                            if cumprob > Tracker["constants"]["ccfpercentage"]:
                                lit = j + 1
                                break

                        #  To here
                        ###if( Blockdata["myid"] == 18 and lima<5):  sxprint("  SECOND KEPT  ",lit)
                        # if( lima<5):  sxprint("  SECOND KEPT  ",lit)

                # mpi_barrier(MPI_COMM_WORLD)
                # mpi_finalize()
                # exit()
                # for j in range(lit):
                # 	 newpar[kl][2].append([int(xod2[j]),float(xod1[j])])
                if lenass > 0:
                    ###print("   EICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im)#,assignments_to_cones)

                    firstdirections = [[0.0, 0.0] for iln in range(lit)]
                    firstshifts = [0] * lit
                    for iln in range(lit):
                        hashparams = int(xod2[iln])
                        ishift = hashparams % 1000
                        ipsiandiang = old_div(hashparams, 1000)
                        # ipsi = ipsiandiang%100000
                        iang = old_div(ipsiandiang, 100000)
                        # try:
                        firstdirections[iln] = [
                            coarse_angles[iang][0],
                            coarse_angles[iang][1],
                            0.0,
                        ]
                        # except:
                        # 	sxprint(" FAILED  ",Blockdata["myid"],Tracker["keepfirst"],iln,lima,lit,hashparams,iang,xod2[:lit],xod1[:lit])
                        # 	mpi_finalize()
                        # 	exit()
                        firstshifts[iln] = ishift
                        #  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):
                        # 	ipsi = ipsiandiang%100000
                        # 	ishift = hashparams%1000
                        # 	###print("  SECONDPAR%04d  "%im,iln,hashparams, ipsi,iang,ishift)
                        # 	sxprint(" SECONDPAR%04d  "%im,iln, refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, shifts[ishift])
                    ###del xod2
                    ###if( Blockdata["myid"] == 18 and lima<5):   sxprint("  SECOND KEPT  ",im,lit)#,len(xod2),xod2,firstdirections,firstshifts)
                    # mpi_barrier(MPI_COMM_WORLD)
                    # mpi_finalize()
                    # exit()

                    # if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  FIFI ",firstdirections)
                    # if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  GUGU ",firstshifts)
                    # Find neighbors, ltabang contains positions on refang list, no psis
                    ###ltabang = nearest_many_full_k_projangles(refang, firstdirections, howmany = 5, sym=Tracker["constants"]["symmetry"])
                    ltabang = find_nearest_k_refangles_to_many_angles(
                        refdirs,
                        firstdirections,
                        Tracker["delta"],
                        howmany=Tracker["howmany"],
                    )
                    ###if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  ltabang ",ltabang)
                    ##mpi_barrier(MPI_COMM_WORLD)
                    ##mpi_finalize()
                    ##exit()

                    # ltabang has length lit, which is the same as length as firshifts.  However, original xod2 is still available,
                    #   even though it is longer than lit.

                    ###ltabshi = [shiftneighbors[iln] for iln in firstshifts]
                    # if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  HUHU ",ltabang)
                    # if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  OUOU ",ltabshi)

                    #  Prepare image for chi2.
                    #  We have to repeat everything from get shrink data, including shifts
                    #  The number of entries is lit=len(ltabang), each proj dir has 4 neighbors, so including itself there are 5,
                    #    there are 3 psis, and at most n_fine_shifts.
                    #    Not not all shifts have to be populated. We may also have repeated triplets of angles, but each can have
                    #      different shifts.  If so, we have to remove duplicates from the entire set.
                    lcod1 = lit * 4 * 2 * n_fine_shifts

                    # cod2 = np.ndarray((lit,5,3,n_fine_shifts),dtype=int,order="C")
                    # cod2.fill(-1)  #  hashparams
                    cod2 = []
                    # lol = 0
                    for i1 in range(lit):
                        hashparams = int(xod2[i1])
                        ipsiandiang = old_div(hashparams, 1000)
                        oldiang = old_div(ipsiandiang, 100000)
                        ipsi = ipsiandiang % 100000
                        ishift = hashparams % 1000
                        tshifts = get_shifts_neighbors(shifts, coarse_shifts[ishift])
                        # if( Blockdata["myid"] == Blockdata["main_node"]): sxprint(" tshifts  ",i1,len(tshifts))
                        for i2 in range(Tracker["howmany"]):
                            iang = ltabang[i1][i2]
                            for i3 in range(2):  # psi
                                itpsi = int(
                                    old_div(
                                        (
                                            coarse_angles[oldiang][2]
                                            + ipsi * coarse_delta
                                            - refang[iang][2]
                                            + 360.0
                                        ),
                                        Tracker["delta"],
                                    )
                                )
                                itpsi = (itpsi + i3) % npsi
                                for i4 in range(len(tshifts)):
                                    cod2.append(
                                        iang * 100000000 + itpsi * 1000 + tshifts[i4]
                                    )
                                    # lol += 1
                                    # if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  zibzi  ",i1,i2,i3,i4, lol)

                    del xod1, xod2

                    ###if( Blockdata["myid"] == 18 and lima<5):   sxprint("  THIRD   ",len(cod2))#,cod2)
                    cod2 = list(set(cod2))
                    cod1 = [[old_div(q, 1000), i] for i, q in enumerate(cod2)]
                    cod1.sort()

                    lit = len(cod1)

                    cod2 = sp_helix_sphire.np.asarray(
                        [cod2[cod1[i][1]] for i in range(lit)]
                    )

                    cod1 = sp_helix_sphire.np.ndarray(lit, dtype="f4", order="C")
                    # cod1.fill(np.finfo(dtype='f4').min)
                    cod3 = sp_helix_sphire.np.ndarray(lit, dtype="f4", order="C")
                    # cod3.fill(0.0)  #  varadj

                    ###if( Blockdata["myid"] == 18 and lima<5): sxprint("  THIRD   ",im,lit)#,cod2)

                    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    #  Third ML part.  Now the actual chi2 distances, we compute reprojections on the fly.
                    #  Make sure volprep has nxinit size
                    tbckg = []
                    data = [None] * nshifts
                    johi = 0
                    iln = 0
                    prevdir = -1
                    while iln < lit:
                        hashparams = cod2[iln]
                        # if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  COD2   ",im,lit,iln,cod2[iln])
                        ipsiandiang = old_div(hashparams, 1000)
                        if ipsiandiang != prevdir:
                            prevdir = ipsiandiang
                            ipsi = ipsiandiang % 100000
                            iang = old_div(ipsiandiang, 100000)
                            ##junk = time()
                            temp = sp_projection.prgl(
                                volinit,
                                [
                                    refang[iang][0],
                                    refang[iang][1],
                                    (refang[iang][2] + ipsi * Tracker["delta"]) % 360.0,
                                    0.0,
                                    0.0,
                                ],
                                1,
                                False,
                            )
                            ##eat += time()-junk
                            temp.set_attr("is_complex", 0)
                            johi += 1
                        while ipsiandiang == old_div(cod2[iln], 1000):
                            hashparams = cod2[iln]
                            ishift = hashparams % 1000
                            if data[ishift] == None:
                                xx = shifts[ishift][0] * shrink
                                yy = shifts[ishift][1] * shrink
                                data[ishift] = sp_fundamentals.fshift(dataml, xx, yy)
                                data[ishift].set_attr("is_complex", 0)
                            ##junk = time()
                            # [peak,varadj] = Util.sqednorm(data[ishift], temp, ctfa, bckgnoise)
                            fofa = EMAN2_cppwrap.Util.sqednormbckg(
                                data[ishift], temp, ctfa, bckgnoise, indx, tfrac, tcount
                            )
                            cod1[iln] = -fofa[-2]  # -peak
                            cod3[iln] = fofa[-1]  # varadj
                            tbckg.append(fofa[:-2])
                            iln += 1
                            if iln == lit:
                                break
                        # if( Blockdata["myid"] == Blockdata["main_node"]):
                        # 	temp.write_image("temp.hdf")
                        # 	data[iln].write_image("data.hdf")
                        # 	ctfa.write_image("ctfa.hdf")
                        # 	bckgnoise.write_image("bckgnoise.hdf")
                        # 	exit()
                        ###if( Blockdata["myid"] == Blockdata["main_node"] and iln%1000 ==0):  sxprint(" progress  ",iln,time()-at)
                    # if( Blockdata["myid"] == Blockdata["main_node"]):
                    # 	sxprint("  PROJECT   ",im,lit,johi)#,cod2)
                    # 	#for iln in range(lit):  sxprint("  ROPE   ",iln,cod1[iln],cod2[iln],cod3[iln])
                    del data
                    del dataml

                    lina = sp_helix_sphire.np.argsort(cod1)
                    cod1 = cod1[lina[::-1]]  # This sorts in reverse order
                    cod2 = cod2[lina[::-1]]  # This sorts in reverse order
                    cod3 = cod3[lina[::-1]]  # This sorts in reverse order
                    cod1 -= cod1[0]
                    tbckg = [tbckg[int(q)] for q in lina]
                    lina = sp_helix_sphire.np.argwhere(
                        cod1 > Tracker["constants"]["expthreshold"]
                    )
                    cod1 = cod1[lina]
                    cod2 = cod2[lina]
                    cod3 = cod3[lina]
                    tbckg = [tbckg[int(q)] for q in lina]

                    ###if( Blockdata["myid"] == Blockdata["main_node"]):
                    ###for iui in range(len(lina)):
                    ###	for iui in range(len(cod1)):
                    ###		sxprint("  MLML  ",iui,cod1[iui],exp(cod1[iui]),cod2[iui],cod3[iui])

                    sp_helix_sphire.np.exp(cod1, out=cod1)
                    cod1 = old_div(cod1, sp_helix_sphire.np.sum(cod1))
                    cumprob = 0.0
                    for j in range(len(cod1)):
                        cumprob += cod1[j]
                        if cumprob > Tracker["constants"]["ccfpercentage"]:
                            lit = j + 1
                            break

                    #  New norm is a sum of eq distances multiplied by their probabilities augmented by PW.
                    norm_per_particle[im] = (
                        sp_helix_sphire.np.sum(cod1[:lit] * cod3[:lit])
                        + accumulatepw[im][reachpw]
                    )
                    atbckg = [0.0] * len(tbckg[0])
                    for iln in range(lit):
                        prob = float(cod1[iln])
                        Blockdata["totprob"][particle_group] += prob
                        for iq in range(len(tbckg[0])):
                            atbckg[iq] += tbckg[iln][iq] * prob

                    del tbckg
                    for iq in range(nxth):
                        Blockdata["newbckgnoise"][iq, particle_group] += atbckg[iq]
                    del atbckg

                    for iln in range(lit):
                        newpar[im][2].append([int(cod2[iln]), float(cod1[iln])])
                        #  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):
                        # 	#print("  NEWPAR%04d  "%im,iln,newpar[im][2][-1])
                        # 	hashparams = newpar[im][2][iln][0]
                        # 	ipsiandiang	= hashparams/1000
                        # 	ipsi = ipsiandiang%100000
                        # 	iang = ipsiandiang/100000
                        # 	ishift = hashparams%1000
                        # 	#print("  NEWPAR%04d  "%im,iln,newpar[im][2][-1],hashparams,ipsi,iang,ishift)
                        # 	sxprint(" NEWPAR%04d  "%im,iln,newpar[im][2][iln], refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, shifts[ishift])

                    ###if( Blockdata["myid"] == 18 and lima<5):   sxprint("  FINALLY  ",im,lit)
                    del cod1, cod2, cod3, lina
                    ###mpi_barrier(MPI_COMM_WORLD)
                    ###mpi_finalize()
                    ###exit()

    """Multiline Comment21"""

    #  END OF CONES
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_global_def.sxprint(
            "  Finished projection matching   %10.1fmin"
            % (old_div((time.time() - at), 60.0))
        )
    at = time.time()
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    #  All images were processed, now to the additional calculations

    ###mpi_barrier(MPI_COMM_WORLD)
    ###mpi_finalize()
    ###exit()

    # norm correction ---- calc the norm correction per particle
    snormcorr = 0.0
    for kl in range(nima):
        ###print("   NORMPERPARTICLE  ",Blockdata["myid"],kl,norm_per_particle[kl])
        norm_per_particle[kl] = old_div(
            numpy.sqrt(norm_per_particle[kl] * 2.0) * oldparams[kl][7],
            Tracker["avgvaradj"][procid],
        )
        snormcorr += norm_per_particle[kl]
    Tracker["avgvaradj"][procid] = snormcorr
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    #  Compute avgvaradj
    Tracker["avgvaradj"][procid] = mpi.mpi_reduce(
        Tracker["avgvaradj"][procid],
        1,
        mpi.MPI_FLOAT,
        mpi.MPI_SUM,
        Blockdata["main_node"],
        mpi.MPI_COMM_WORLD,
    )
    if Blockdata["myid"] == Blockdata["main_node"]:
        Tracker["avgvaradj"][procid] = old_div(
            float(Tracker["avgvaradj"][procid]), Tracker["nima_per_chunk"][procid]
        )
    else:
        Tracker["avgvaradj"][procid] = 0.0
    Tracker["avgvaradj"][procid] = sp_utilities.bcast_number_to_all(
        Tracker["avgvaradj"][procid], Blockdata["main_node"]
    )
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    #  Compute statistics of smear -----------------
    smax = -1000000
    smin = 1000000
    sava = 0.0
    svar = 0.0
    snum = 0
    for kl in range(nima):
        j = len(newpar[kl][2])
        snum += 1
        sava += float(j)
        svar += j * float(j)
        smax = max(smax, j)
        smin = min(smin, j)
    snum = mpi.mpi_reduce(
        snum, 1, mpi.MPI_INT, mpi.MPI_SUM, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    sava = mpi.mpi_reduce(
        sava, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    svar = mpi.mpi_reduce(
        svar, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    smax = mpi.mpi_reduce(
        smax, 1, mpi.MPI_INT, mpi.MPI_MAX, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    smin = mpi.mpi_reduce(
        smin, 1, mpi.MPI_INT, mpi.MPI_MIN, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    if Blockdata["myid"] == 0:
        sava = old_div(float(sava), snum)
        svar = numpy.sqrt(
            max(0.0, old_div((float(svar) - snum * sava ** 2), (snum - 1)))
        )
        # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
        line = ""
        sp_global_def.sxprint(
            line,
            "Smear stat  (number of images, ave, sumsq, min, max)):  %7d    %12.3g   %12.3g  %7d  %7d"
            % (snum, sava, svar, smin, smax),
        )

    at = time.time()
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    ###if Blockdata["myid"] == Blockdata["main_node"]:  print "  Finished :",time()-at
    mpi.mpi_win_free(win_sm)
    if Tracker["nxinit"] != Tracker["nxpolar"]:
        mpi.mpi_win_free(win_volinit)
        emnumpy4.unregister_numpy_from_emdata()
        del emnumpy4
    else:
        mpi.mpi_win_free(win_vol)

    mpi.mpi_barrier(Blockdata["shared_comm"])
    emnumpy1.unregister_numpy_from_emdata()
    emnumpy2.unregister_numpy_from_emdata()
    del emnumpy1, emnumpy2

    del volinit

    mpi.mpi_barrier(Blockdata["shared_comm"])

    # Compute new background noise
    # Reduce stuff
    if procid == 1:
        Blockdata["totprob"] = mpi.mpi_reduce(
            Blockdata["totprob"],
            nyb,
            mpi.MPI_FLOAT,
            mpi.MPI_SUM,
            Blockdata["main_node"],
            mpi.MPI_COMM_WORLD,
        )
        sp_utilities.reduce_EMData_to_root(
            Blockdata["newbckgnoise"], Blockdata["myid"], Blockdata["main_node"]
        )
        if Blockdata["myid"] == 0:
            for igrp in range(nyb):
                Blockdata["newbckgnoise"][0, igrp] = 1.0
                for i in range(1, nxth):
                    if Blockdata["newbckgnoise"][i, igrp] > 0.0:
                        Blockdata["newbckgnoise"][i, igrp] = old_div(
                            2.0 * Blockdata["totprob"][igrp],
                            Blockdata["newbckgnoise"][i, igrp],
                        )  # normalize and invert
                for i in range(nxth, nxb):
                    Blockdata["newbckgnoise"][i, igrp] = Blockdata["bckgnoise"][igrp][i]
            Blockdata["newbckgnoise"].write_image(
                os.path.join(Tracker["directory"], "bckgnoise.hdf")
            )  # Write updated bckgnoise to current directory

        sp_utilities.bcast_EMData_to_all(
            Blockdata["newbckgnoise"],
            Blockdata["myid"],
            source_node=Blockdata["main_node"],
            comm=mpi.MPI_COMM_WORLD,
        )
        for igrp in range(nyb):
            for i in range(nxb):
                Blockdata["bckgnoise"][igrp][i] = Blockdata["newbckgnoise"][i, igrp]
        del Blockdata["newbckgnoise"]

    # print("  NORMALIZATION DONE  ",Blockdata["myid"])
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_global_def.sxprint(
            "  Finished sigma2   %10.1fmin" % (old_div((time.time() - at), 60.0))
        )
    return newpar, norm_per_particle


# TRUE POLAR


def ali3D_local_polar(
    refang,
    shifts,
    coarse_angles,
    coarse_shifts,
    procid,
    original_data=None,
    oldparams=None,
    preshift=False,
    apply_mask=True,
    nonorm=False,
    applyctf=True,
):
    """
	02/07/2017
	"""
    global Tracker, Blockdata
    #  Input data has to be CTF-multiplied, preshifted
    #  Output - newpar, see structure
    #    newpar = [[i, [worst_similarity, sum_all_similarities], [[-1, -1.0e23] for j in range(Tracker["lentop"])]] for i in range(len(data))]
    #    newpar = [[params],[],... len(data)]
    #    params = [particleID, [worst_similarity, sum_all_similarities],[imageallparams]]]
    #    imageallparams = [[orientation, similarity],[],...  number of all orientations ]
    #  Coding of orientations:
    #    hash = ang*100000000 + lpsi*1000 + ishift
    #    ishift = hash%1000
    #    ipsi = (hash/1000)%100000
    #    iang  = hash/100000000
    #  To get best matching for particle #kl:
    #     hash_best = newpar[kl][-1][0][0]
    #     best_sim  = newpar[kl][-1][0][1]
    #  To sort:
    # from operator 		import itemgetter#, attrgetter, methodcaller
    #   params.sort(key=itemgetter(2))
    # from fundamentals import resample

    if Blockdata["myid"] == Blockdata["main_node"]:
        print_dict(Tracker, "PROJECTION MATCHING parameters of buffered local polar")

    at = time.time()
    shrinkage = old_div(float(Tracker["nxpolar"]), float(Tracker["constants"]["nnxo"]))
    shrink = old_div(float(Tracker["nxinit"]), float(Tracker["constants"]["nnxo"]))
    radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
    mode = "F"
    numr = Numrinit_local(1, radius, 1, mode)
    wr = sp_alignment.ringwe(numr, mode)
    cnx = float(old_div(Tracker["nxpolar"], 2) + 1)

    ##if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  RADIUS IN POLAR  ",radius,numr[-1])
    #  FINE SEARCH CONSTANTS
    nang = len(refang)
    ac_fine = numpy.cos(numpy.radians(Tracker["an"]))
    npsi = int(old_div(360.0, Tracker["delta"]))
    mpsi = 2
    c_fine_psi = old_div(mpsi, 2)
    nshifts = len(shifts)
    n_fine_shifts = 4

    nima = len(original_data)
    mask = EMAN2_cppwrap.Util.unrollmask(Tracker["nxinit"], Tracker["nxinit"])
    for j in range(old_div(Tracker["nxinit"], 2), Tracker["nxinit"]):
        mask[0, j] = 1.0
    mask2D = sp_utilities.model_circle(
        Tracker["constants"]["radius"],
        Tracker["constants"]["nnxo"],
        Tracker["constants"]["nnxo"],
    )
    reachpw = (
        mask.get_xsize()
    )  # The last element of accumulated pw is zero so for the full size nothing is added.

    #  COARSE SEARCH CONSTANTS
    n_coarse_ang = len(coarse_angles)
    coarse_delta = 2 * Tracker["delta"]
    n_coarse_psi = int(old_div(360.0, coarse_delta))
    ###m_coarse_psi = int((2*2*Tracker["an"])/coarse_delta + 0.5) + 1
    m_coarse_psi = int(old_div(Tracker["an"], Tracker["delta"]) + 0.5)
    c_coarse_psi = old_div(m_coarse_psi, 2)
    n_coarse_shifts = len(coarse_shifts)

    coarse_shifts_shrank = [None] * n_coarse_shifts
    for ib in range(n_coarse_shifts):
        coarse_shifts_shrank[ib] = [
            coarse_shifts[ib][0] * shrinkage,
            coarse_shifts[ib][1] * shrinkage,
        ]

    ###if(Blockdata["myid"] == Blockdata["main_node"]): sxprint("   TRETR  ",Tracker["constants"]["nnxo"],Tracker["nxinit"],reachpw,n_coarse_ang,coarse_delta,n_coarse_psi,m_coarse_psi,c_coarse_psi,n_coarse_shifts)

    # if(Blockdata["myid"] == Blockdata["main_node"]):
    # 	sxprint( original_data[0].get_attr("identifier") )
    # 	sxprint( original_data[1].get_attr("identifier") )

    disp_unit = sp_helix_sphire.np.dtype("f4").itemsize

    #  REFVOL
    if Blockdata["myid_on_node"] == 0:
        if Blockdata["myid"] == Blockdata["main_node"]:
            odo = get_refvol(Tracker["nxpolar"])
            nxvol = odo.get_xsize()
            nyvol = odo.get_ysize()
            nzvol = odo.get_zsize()
        else:
            nxvol = 0
            nyvol = 0
            nzvol = 0

        nxvol = sp_utilities.bcast_number_to_all(
            nxvol,
            source_node=Blockdata["main_node"],
            mpi_comm=Blockdata["group_zero_comm"],
        )
        nyvol = sp_utilities.bcast_number_to_all(
            nyvol,
            source_node=Blockdata["main_node"],
            mpi_comm=Blockdata["group_zero_comm"],
        )
        nzvol = sp_utilities.bcast_number_to_all(
            nzvol,
            source_node=Blockdata["main_node"],
            mpi_comm=Blockdata["group_zero_comm"],
        )

        if Blockdata["myid"] != Blockdata["main_node"]:
            odo = sp_utilities.model_blank(nxvol, nyvol, nzvol)

        sp_utilities.bcast_EMData_to_all(
            odo,
            Blockdata["group_zero_myid"],
            source_node=Blockdata["main_node"],
            comm=Blockdata["group_zero_comm"],
        )

        odo = sp_projection.prep_vol(odo, npad=2, interpolation_method=1)
        ndo = EMAN2_cppwrap.EMNumPy.em2numpy(odo)
        nxvol = odo.get_xsize()
        nyvol = odo.get_ysize()
        nzvol = odo.get_zsize()
        orgsizevol = nxvol * nyvol * nzvol
        sizevol = orgsizevol
    else:
        orgsizevol = 0
        sizevol = 0
        nxvol = 0
        nyvol = 0
        nzvol = 0

    orgsizevol = sp_utilities.bcast_number_to_all(
        orgsizevol, source_node=Blockdata["main_node"]
    )
    nxvol = sp_utilities.bcast_number_to_all(nxvol, source_node=Blockdata["main_node"])
    nyvol = sp_utilities.bcast_number_to_all(nyvol, source_node=Blockdata["main_node"])
    nzvol = sp_utilities.bcast_number_to_all(nzvol, source_node=Blockdata["main_node"])

    win_vol, base_vol = mpi.mpi_win_allocate_shared(
        sizevol * disp_unit, disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"]
    )
    sizevol = orgsizevol
    if Blockdata["myid_on_node"] != 0:
        base_vol, = mpi.mpi_win_shared_query(win_vol, mpi.MPI_PROC_NULL)
    """
    Comments from Adnan:
    numpy.core.mutliarray.int_asbuffer function is not available anymore . That is why a different alternative is tried here.
    """
    # volbuf = sp_helix_sphire.np.frombuffer(
    #     sp_helix_sphire.np.core.multiarray.int_asbuffer(base_vol, sizevol * disp_unit),
    #     dtype="f4",
    # )
    ptr = ctypes.cast(base_vol, ctypes.POINTER(ctypes.c_int * sizevol))
    volbuf = numpy.frombuffer(ptr.contents, dtype="f4")
    volbuf = volbuf.reshape(nzvol, nyvol, nxvol)
    if Blockdata["myid_on_node"] == 0:
        sp_helix_sphire.np.copyto(volbuf, ndo)
        del odo, ndo

    # volprep = EMNumPy.assign_numpy_to_emdata(volbuf)
    emnumpy1 = EMAN2_cppwrap.EMNumPy()
    volprep = emnumpy1.register_numpy_to_emdata(volbuf)
    volprep.set_attr_dict(
        {
            "is_complex": 1,
            "is_complex_x": 0,
            "is_fftodd": 0,
            "is_fftpad": 1,
            "is_shuffled": 1,
            "npad": 2,
        }
    )

    if Tracker["nxinit"] != Tracker["nxpolar"]:
        mpi.mpi_win_free(win_vol)
        #  REFVOL FOR ML
        if Blockdata["myid_on_node"] == 0:
            if Blockdata["myid"] == Blockdata["main_node"]:
                odo = get_refvol(Tracker["nxpolar"])
                nxvol = odo.get_xsize()
                nyvol = odo.get_ysize()
                nzvol = odo.get_zsize()
            else:
                nxvol = 0
                nyvol = 0
                nzvol = 0

            nxvol = sp_utilities.bcast_number_to_all(
                nxvol,
                source_node=Blockdata["main_node"],
                mpi_comm=Blockdata["group_zero_comm"],
            )
            nyvol = sp_utilities.bcast_number_to_all(
                nyvol,
                source_node=Blockdata["main_node"],
                mpi_comm=Blockdata["group_zero_comm"],
            )
            nzvol = sp_utilities.bcast_number_to_all(
                nzvol,
                source_node=Blockdata["main_node"],
                mpi_comm=Blockdata["group_zero_comm"],
            )

            if Blockdata["myid"] != Blockdata["main_node"]:
                odo = sp_utilities.model_blank(nxvol, nyvol, nzvol)

            sp_utilities.bcast_EMData_to_all(
                odo,
                Blockdata["group_zero_myid"],
                source_node=Blockdata["main_node"],
                comm=Blockdata["group_zero_comm"],
            )

            odo = sp_projection.prep_vol(odo, npad=2, interpolation_method=1)
            ndo = EMAN2_cppwrap.EMNumPy.em2numpy(odo)
            nxvol = odo.get_xsize()
            nyvol = odo.get_ysize()
            nzvol = odo.get_zsize()
            orgsizevol = nxvol * nyvol * nzvol
            sizevol = orgsizevol
        else:
            orgsizevol = 0
            sizevol = 0
            nxvol = 0
            nyvol = 0
            nzvol = 0

        orgsizevol = sp_utilities.bcast_number_to_all(
            orgsizevol, source_node=Blockdata["main_node"]
        )
        nxvol = sp_utilities.bcast_number_to_all(
            nxvol, source_node=Blockdata["main_node"]
        )
        nyvol = sp_utilities.bcast_number_to_all(
            nyvol, source_node=Blockdata["main_node"]
        )
        nzvol = sp_utilities.bcast_number_to_all(
            nzvol, source_node=Blockdata["main_node"]
        )

        win_volinit, base_volinit = mpi.mpi_win_allocate_shared(
            sizevol * disp_unit, disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"]
        )
        sizevol = orgsizevol
        if Blockdata["myid_on_node"] != 0:
            base_volinit, = mpi.mpi_win_shared_query(win_volinit, mpi.MPI_PROC_NULL)
        """
        Comments from Adnan:
        numpy.core.mutliarray.int_asbuffer function is not available anymore . That is why a different alternative is tried here.
        """
        # volbufinit = sp_helix_sphire.np.frombuffer(
        #     sp_helix_sphire.np.core.multiarray.int_asbuffer(
        #         base_volinit, sizevol * disp_unit
        #     ),
        #     dtype="f4",
        # )
        ptr = ctypes.cast(base_volinit, ctypes.POINTER(ctypes.c_int * sizevol))
        volbufinit = numpy.frombuffer(ptr.contents, dtype="f4")
        volbufinit = volbufinit.reshape(nzvol, nyvol, nxvol)
        if Blockdata["myid_on_node"] == 0:
            sp_helix_sphire.np.copyto(volbufinit, ndoinit)
            del odo, ndoinit

        # volinit = EMNumPy.assign_numpy_to_emdata(volbufinit)
        emnumpy4 = EMAN2_cppwrap.EMNumPy()
        volinit = emnumpy4.register_numpy_to_emdata(volbufinit)
        volinit.set_attr_dict(
            {
                "is_complex": 1,
                "is_complex_x": 0,
                "is_fftodd": 0,
                "is_fftpad": 1,
                "is_shuffled": 1,
                "npad": 2,
            }
        )
        if Blockdata["myid_on_node"] == 0:
            volinit.update()
        mpi.mpi_barrier(Blockdata["shared_comm"])
        ###if( Blockdata["myid"] < 10 ):
        ###	sxprint(" VOLPREP  ",volinit[0],volinit[1],info(volinit))
    else:
        volinit = volprep
    #  End of replaced volprep

    #  START CONES
    #  This has to be systematically done per node
    #
    crefim = EMAN2_cppwrap.Util.Polar2Dm(
        sp_utilities.model_blank(Tracker["nxpolar"], Tracker["nxpolar"]),
        cnx,
        cnx,
        numr,
        mode,
    )
    size_of_one_image = crefim.get_xsize()
    #  We will assume half of the memory is available.  We will do it betteer later.
    numberofrefs_inmem = int(
        old_div(
            old_div(Tracker["constants"]["memory_per_node"], 4),
            (old_div((size_of_one_image * disp_unit), 1.0e9)),
        )
    )
    ####if( Blockdata["myid_on_node"] == 0  ):  sxprint( " MEMEST ", n_coarse_ang,numberofrefs_inmem)
    #  number of references that will fit into one mode
    normals_set = sp_utilities.angles_to_normals(coarse_angles)
    Blockdata["angle_set"] = coarse_angles
    if n_coarse_ang <= numberofrefs_inmem:
        number_of_cones = 1
        numberofrefs_inmem = n_coarse_ang
        assignments_to_cones = [list(range(len(oldparams)))]

        assignments_of_refangles_to_cones = [
            [] for i in range(len(assignments_to_cones))
        ]
        assignments_of_refangles_to_angles = [
            [] for i in range(nima)
        ]  # for each myid separately, these are angles on this myid

        for i, q in enumerate(assignments_to_cones):
            #  find assignments of refdirs to this cone within an of each angledirs, so we need symmetry neighbors set of refdirs.
            ###print( "in loop ", Blockdata["myid"],i,len(q))#,q

            for m in q:
                # print " m ",m,len(angles)

                assignments_of_refangles_to_angles[
                    m
                ] = find_assignments_of_refangles_to_angles(
                    normals_set, oldparams[m], Tracker["an"]
                )
                assignments_of_refangles_to_cones[i].extend(
                    assignments_of_refangles_to_angles[m]
                )

            assignments_of_refangles_to_cones[i] = list(
                set(assignments_of_refangles_to_cones[i])
            )
            # print(  " assignments_of_refangles_to_cones on myid ",Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]) )
            assignments_of_refangles_to_cones[i] = sp_utilities.wrap_mpi_gatherv(
                assignments_of_refangles_to_cones[i], 0, Blockdata["shared_comm"]
            )
            doit = 1
            if Blockdata["myid_on_node"] == 0:
                assignments_of_refangles_to_cones[i] = list(
                    set(assignments_of_refangles_to_cones[i])
                )
                # print(  " assignments_of_refangles_to_cones on zero ",i,len(assignments_of_refangles_to_cones[i]) )#,assignments_of_refangles_to_cones[i]

        for i, q in enumerate(assignments_of_refangles_to_cones):
            ###if( Blockdata["myid_on_node"] == 0 ): sxprint( " which refangles belong to which cone",Blockdata["color"],Blockdata["myid"],i,len(q) )#,q
            assignments_of_refangles_to_cones[i] = sp_utilities.wrap_mpi_bcast(
                q, 0, Blockdata["shared_comm"]
            )

    else:

        angledirs = sp_utilities.angles_to_normals([u1[:3] for u1 in oldparams])

        number_of_cones = max(
            2, int(old_div(n_coarse_ang, numberofrefs_inmem) * 1.5 + 0.5)
        )
        ###if Blockdata["myid_on_node"] == 0:  sxprint( " LENGTHS  ",Blockdata["color"],nima,n_coarse_ang, numberofrefs_inmem,number_of_cones)
        cont = True
        while cont:
            #  Translate number of cones to cone_delta
            cone_delta, number_of_cones = number_of_cones_to_delta(number_of_cones)
            #  Generate cone_angles
            ##if Blockdata["myid"] == 0:  sxprint( "  WHILE  ",number_of_cones, cone_delta)
            if Blockdata["symclass"].sym[0] == "c":
                if number_of_cones == 1:
                    cone_delta = 180.1
                    cone_angles = [[0.0, 1.0, 0.0]]
                else:
                    cone_angles = Blockdata["symclass"].even_angles(
                        cone_delta, theta1=1.0, theta2=89.0
                    )
                    cone_angles += [
                        [(q[0] + 90.0) % 360.0, 180.0 - q[1], 0] for q in cone_angles
                    ]
                    cone_angles = Blockdata["symclass"].reduce_anglesets(cone_angles)
            else:
                cone_angles = Blockdata["symclass"].even_angles(cone_delta, theta1=1.0)

            # if Blockdata["myid"] == 0:  sxprint(  "  number of cones ",number_of_cones,cone_delta, len(cone_angles))
            assert number_of_cones == len(cone_angles)

            conedirs = sp_utilities.angles_to_normals(
                Blockdata["symclass"].symmetry_neighbors(cone_angles)
            )
            neighbors = old_div(
                len(conedirs), len(cone_angles)
            )  # Symmetry related neighbors
            # if Blockdata["myid"] == 0:  sxprint(  "  neighbors  ",Blockdata["myid"],neighbors, cone_angles)
            #  assign data directions to cone_angles
            assignments_to_cones = sp_utilities.assign_projdirs_f(
                angledirs, conedirs, neighbors
            )
            ###print(  " assignments_to_cones ",Blockdata["myid"],len(assignments_to_cones),[len(q) for q in assignments_to_cones],assignments_to_cones[0])
            #  the above should have length of refdirs and each should have indexes of data that belong to this cone
            del conedirs
            # print "assignments_to_cones ",assignments_to_cones
            #  For each cone we have to find which refangles are needed to do the matching
            assignments_of_refangles_to_cones = [
                [] for i in range(len(assignments_to_cones))
            ]
            assignments_of_refangles_to_angles = [
                [] for i in range(nima)
            ]  # for each myid separately, these are angles on this myid

            for i, q in enumerate(assignments_to_cones):
                #  find assignments of refdirs to this cone within an of each angledirs, so we need symmetry neighbors set of refdirs.
                # if Blockdata["myid"] == 0:  sxprint( "in loop ", Blockdata["myid"],i,len(q),q)

                if len(q) == 0:
                    # empty assignment, on a given CPU there are no images assigned to a given cone
                    assignments_of_refangles_to_cones[i] = [-1]
                else:
                    for m in q:
                        # print " m ",m,len(angles)

                        assignments_of_refangles_to_angles[
                            m
                        ] = find_assignments_of_refangles_to_angles(
                            normals_set, oldparams[m], Tracker["an"]
                        )
                        # if Blockdata["myid"] == 0:  sxprint( "assignments_of_refangles_to_angles[m] ", Blockdata["color"],i,m,assignments_of_refangles_to_angles[m])
                        assignments_of_refangles_to_cones[i].extend(
                            assignments_of_refangles_to_angles[m]
                        )

                    # if( Blockdata["myid"] == 19 ):  sxprint(  " doit0 ",Blockdata["myid"], i,assignments_of_refangles_to_cones[i],q,assignments_to_cones)
                    assignments_of_refangles_to_cones[i] = list(
                        set(assignments_of_refangles_to_cones[i])
                    )
                ###if Blockdata["myid"] == 19:  sxprint(  " assignments_of_refangles_to_cones on myid ",Blockdata["color"],Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]), assignments_of_refangles_to_cones[i] )
                assignments_of_refangles_to_cones[i] = sp_utilities.wrap_mpi_gatherv(
                    assignments_of_refangles_to_cones[i], 0, Blockdata["shared_comm"]
                )
                ###if Blockdata["myid"] == 0:  sxprint(  " assignments_of_refangles_to_cones gatherv",Blockdata["color"],Blockdata["myid"],i,len(assignments_of_refangles_to_cones[i]) )
                doit = 1
                if Blockdata["myid_on_node"] == 0:
                    assignments_of_refangles_to_cones[i] = list(
                        set(assignments_of_refangles_to_cones[i]) - set([-1])
                    )
                    ###if( Blockdata["myid_on_node"] == 0 ):  sxprint(  " assignments_of_refangles_to_cones on zero ",i,len(assignments_of_refangles_to_cones[i]))
                    #  POSSIBLE PROBLEM - IT IS POSSIBLE FOR A GIVEN CONE TO HAVE NO REFANGLES AND THUS BE EMPTY
                    if len(assignments_of_refangles_to_cones[i]) > numberofrefs_inmem:
                        number_of_cones = int(number_of_cones * 1.25)
                        # print(  " increased number_of_cones ",i,number_of_cones )
                        doit = 0
                doit = sp_utilities.bcast_number_to_all(doit, source_node=0)
                number_of_cones = sp_utilities.bcast_number_to_all(
                    number_of_cones, source_node=0, mpi_comm=Blockdata["shared_comm"]
                )
                ###if( Blockdata["myid"] == 19 ):  sxprint(  " doit ",Blockdata["myid"], i,doit ,assignments_of_refangles_to_cones[i],assignments_to_cones)
                if doit == 0:
                    break

            if doit == 1:
                cont = False

        for i, q in enumerate(assignments_of_refangles_to_cones):
            ###if( Blockdata["myid_on_node"] == 0 ): sxprint( " which refangles belong to which cone IOUT ",Blockdata["color"],Blockdata["myid"],i,len(q))#,q
            assignments_of_refangles_to_cones[i] = sp_utilities.wrap_mpi_bcast(
                q, 0, Blockdata["shared_comm"]
            )
            ###if( Blockdata["myid_on_node"] == 0 ): sxprint( " which refangles belong to which cone XOUT ",Blockdata["color"],Blockdata["myid"],i,len(q))#,q
            # if( myid == 1 ):
            # 	print " which refangles belong to which cone",myid,i,len(assignments_of_refangles_to_cones[i])#,q

    if Blockdata["myid"] == 0:
        sp_global_def.sxprint(" number_of_cones: ", number_of_cones)
    #  Maximum number of refangles assigned to angles (max number of references per image)
    nlocal_angles = max([len(q) for q in assignments_of_refangles_to_angles])
    #  Most likely we have to delete some lists before proceeding
    del normals_set
    # We have to figure the maximum length of xod1, which is lang.  If there are no cones, it is an estimate.  If there are cones, we have list of assignments
    #  For number of cones I use refang and an.  This should give approximately the same number as coarse angles and 2*an, which is what is actually used in searches
    numberofrefs_inmem = max([len(q) for q in assignments_of_refangles_to_cones])

    ###for i,q in enumerate(assignments_of_refangles_to_cones):
    ###	sxprint( " which refangles belong to which cone OUT ",Blockdata["color"],Blockdata["myid_on_node"],Blockdata["myid"],i,numberofrefs_inmem,nlocal_angles,len(q))#,q

    #  BIG BUFFER
    lenbigbuf = numberofrefs_inmem  # MAXIMUM NUMBER OF REFERENCE IN ONE OF THE CONES
    orgsize = (
        lenbigbuf * size_of_one_image
    )  # This is number of projections to be computed simultaneously times their size

    if Blockdata["myid_on_node"] == 0:
        size = orgsize
    else:
        size = 0

    win_sm, base_ptr = mpi.mpi_win_allocate_shared(
        size * disp_unit, disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"]
    )
    size = orgsize
    if Blockdata["myid_on_node"] != 0:
        base_ptr, = mpi.mpi_win_shared_query(win_sm, mpi.MPI_PROC_NULL)
    """
    Comments from Adnan:
    numpy.core.mutliarray.int_asbuffer function is not available anymore . That is why a different alternative is tried here.
    """
    # buffer = sp_helix_sphire.np.frombuffer(
    #     sp_helix_sphire.np.core.multiarray.int_asbuffer(base_ptr, size * disp_unit),
    #     dtype="f4",
    # )

    ptr = ctypes.cast(base_ptr, ctypes.POINTER(ctypes.c_int * size))
    buffer = numpy.frombuffer(ptr.contents, dtype="f4")
    buffer = buffer.reshape(lenbigbuf, size_of_one_image)
    # bigbuffer = EMNumPy.assign_numpy_to_emdata(buffer)

    emnumpy2 = EMAN2_cppwrap.EMNumPy()
    bigbuffer = emnumpy2.register_numpy_to_emdata(buffer)

    #  end of CONES setup

    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_global_def.sxprint("  ")
        # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
        line = ""
        sp_global_def.sxprint(
            line,
            "Processing data for polar onx: %3d, nx: %3d, nxpolar: %3d, CTF: %s, preshift: %s,  MEM: %6.2fGB."
            % (
                Tracker["constants"]["nnxo"],
                Tracker["nxinit"],
                Tracker["nxpolar"],
                Tracker["constants"]["CTF"],
                preshift,
                old_div(orgsize, 1.0e9),
            ),
        )

    #  Note these are in Fortran notation for polar searches
    txm = float(Tracker["nxpolar"] - (old_div(Tracker["nxpolar"], 2) + 1) - radius)
    txl = float(radius - old_div(Tracker["nxpolar"], 2) + 1)

    if Blockdata["bckgnoise"]:
        oneover = []
        nnx = Blockdata["bckgnoise"][0].get_xsize()
        for i in range(len(Blockdata["bckgnoise"])):
            temp = [0.0] * nnx
            for k in range(nnx):
                if Blockdata["bckgnoise"][i].get_value_at(k) > 0.0:
                    temp[k] = old_div(
                        1.0, numpy.sqrt(Blockdata["bckgnoise"][i].get_value_at(k))
                    )
            oneover.append(temp)
        del temp

    accumulatepw = [0.0] * nima
    norm_per_particle = [0.0] * nima

    ##lxod1 = lang*len(list_of_coarse_shifts)*(int(2*Tracker["an"]/coarse_delta+0.5)+1)
    newpar = [[i, [0.0], []] for i in range(nima)]

    #  This is for auxiliary function searches.
    Blockdata["angle_set"] = refang
    #  Extract normals from rotation matrices
    refdirs = sp_utilities.angles_to_normals(refang)

    #  We have to make sure the number of cones is the same on all nodes, this is due to the strange MPI problem
    #   that forced me to insert overall barrier into iterations over cones
    max_number_of_cones = number_of_cones
    max_number_of_cones = mpi.mpi_reduce(
        max_number_of_cones, 1, mpi.MPI_INT, mpi.MPI_MAX, 0, mpi.MPI_COMM_WORLD
    )
    if Blockdata["myid"] == Blockdata["main_node"]:
        max_number_of_cones = int(max_number_of_cones[0])
    max_number_of_cones = mpi.mpi_bcast(
        max_number_of_cones, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD
    )
    max_number_of_cones = int(max_number_of_cones[0])

    # if( Blockdata["myid_on_node"] == 0):
    # 	sxprint("  FIRST  nima, nang, npsi, nshift, n_coarse_ang, nlocal_angles, n_coarse_psi, len(list_of_coarse_shifts), number_of_cones , max_number_of_cones ",nima, nang,npsi,nshifts,n_coarse_ang,nlocal_angles,n_coarse_psi, len(coarse_shifts),number_of_cones,max_number_of_cones)

    ##firsti = True
    at = time.time()
    ##eat = 0.0
    lima = 0  # total counter of images
    #  PROCESSING OF CONES
    for icone in range(max_number_of_cones):
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        if (
            icone < number_of_cones
        ):  # This is executed for individual number of cones, some nodes may have fewer.
            nang_start, nang_end = sp_applications.MPI_start_end(
                len(assignments_of_refangles_to_cones[icone]),
                Blockdata["no_of_processes_per_group"],
                Blockdata["myid_on_node"],
            )
            # if(Blockdata["color"] == 1):
            ###print( " ZXZ11  ",Blockdata["color"],Blockdata["myid_on_node"],Blockdata["myid"],len(assignments_of_refangles_to_cones[icone]),nang_start, nang_end)

            for i in range(
                nang_start, nang_end, 1
            ):  # This will take care of no of process on a node less than nang.  Some loops will not be executed
                ic = assignments_of_refangles_to_cones[icone][i]
                temp = sp_projection.prgl(
                    volprep,
                    [coarse_angles[ic][0], coarse_angles[ic][1], 0.0, 0.0, 0.0],
                    1,
                    True,
                )
                crefim = EMAN2_cppwrap.Util.Polar2Dm(temp, cnx, cnx, numr, mode)
                EMAN2_cppwrap.Util.Frngs(crefim, numr)
                EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)
                bigbuffer.insert_clip(crefim, (0, i))

            mpi.mpi_barrier(Blockdata["shared_comm"])

            # if(Blockdata["myid"] == Blockdata["main_node"]):
            # 	sxprint( "  Reference projections generated for cone %4d: %10.1fmin"%(icone,(time()-at)/60.))
            ###print("  TOPOP    ",Blockdata["myid"],MPI_COMM_WORLD,len(assignments_to_cones[icone]))
            ###mpi_finalize()
            ###exit()
            ###  <><><><><><><><>

            #  Preprocess the data

            #  We only process images in the current cone.  icnm is consecutive number, im the actual image number
            # for icnm,im in enumerate(assignments_to_cones[icone]):
            lenass = len(assignments_to_cones[icone])
            ###print("   ENTERING  ",Blockdata["myid"],icone,lenass)
            for icnm in range(
                max(1, lenass)
            ):  # I have to enter the loop even it there is no assignment
                if lenass == 0:
                    keepf = -1
                    ###print("  FOUNDEMPTY  ",Blockdata["myid"],icone,icnm,len(assignments_to_cones[icone]),assignments_to_cones)
                else:
                    #  PREPARE ONE IMAGE
                    im = assignments_to_cones[icone][icnm]

                    phi, theta, psi, sx, sy, wnorm = (
                        oldparams[im][0],
                        oldparams[im][1],
                        oldparams[im][2],
                        oldparams[im][3],
                        oldparams[im][4],
                        oldparams[im][7],
                    )

                    #  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):  sxprint("\n\n  INPUT PARAMS  ",im,phi,theta,psi,sx,sy)
                    if preshift:
                        sx = int(round(sx))
                        sy = int(round(sy))
                        dataimage = sp_fundamentals.cyclic_shift(
                            original_data[im], sx, sy
                        )
                        #  Put rounded shifts on the list, note image has the original floats - check whether it may cause problems
                        oldparams[im][3] = sx
                        oldparams[im][4] = sy
                        sx = 0.0
                        sy = 0.0
                    else:
                        dataimage = original_data[im].copy()

                    st = get_image_statistics(dataimage, mask2D, False)
                    dataimage -= st[0]
                    dataimage = old_div(dataimage, st[1])
                    if dataimage.get_attr_default("bckgnoise", None):
                        dataimage.delete_attr("bckgnoise")
                    #  Do bckgnoise if exists
                    if Blockdata["bckgnoise"]:
                        if apply_mask:
                            if Tracker["constants"]["hardmask"]:
                                dataimage = sp_morphology.cosinemask(
                                    dataimage, radius=Tracker["constants"]["radius"]
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
                                    sp_filter.filt_table(
                                        bckg,
                                        oneover[dataimage.get_attr("particle_group")],
                                    )
                                )
                                #  Normalize bckg noise in real space, only region actually used.
                                st = get_image_statistics(bckg, mask2D, False)
                                bckg -= st[0]
                                bckg = old_div(bckg, st[1])
                                dataimage = sp_morphology.cosinemask(
                                    dataimage,
                                    radius=Tracker["constants"]["radius"],
                                    bckg=bckg,
                                )
                    else:
                        #  if no bckgnoise, do simple masking instead
                        if apply_mask:
                            dataimage = sp_morphology.cosinemask(
                                dataimage, radius=Tracker["constants"]["radius"]
                            )

                    #  Apply varadj
                    if not nonorm:
                        EMAN2_cppwrap.Util.mul_scalar(
                            dataimage, old_div(Tracker["avgvaradj"][procid], wnorm)
                        )

                    ###  FT
                    dataimage = sp_fundamentals.fft(dataimage)
                    sig = EMAN2_cppwrap.Util.rotavg_fourier(dataimage)
                    accumulatepw[im] = sig[old_div(len(sig), 2) :] + [0.0]

                    #  We have to make sure the shifts are within correct range, shrinkage or not
                    # set_params_proj(dataimage,[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
                    if Blockdata["bckgnoise"]:
                        temp = Blockdata["bckgnoise"][
                            dataimage.get_attr("particle_group")
                        ]
                        bckgn = EMAN2_cppwrap.Util.unroll1dpw(
                            Tracker["nxinit"],
                            Tracker["nxinit"],
                            [temp[i] for i in range(temp.get_xsize())],
                        )
                    else:
                        bckgn = EMAN2_cppwrap.Util.unroll1dpw(
                            Tracker["nxinit"], Tracker["nxinit"], [1.0] * 600
                        )
                    bckgnoise = bckgn.copy()
                    for j in range(
                        old_div(Tracker["nxinit"], 2) + 1, Tracker["nxinit"]
                    ):
                        bckgn[0, j] = bckgn[0, Tracker["nxinit"] - j]

                    if Tracker["constants"]["CTF"]:
                        ctf_params = dataimage.get_attr("ctf")
                        ctf_params.apix = old_div(
                            ctf_params.apix,
                            (
                                old_div(
                                    float(Tracker["nxinit"]),
                                    float(Tracker["constants"]["nnxo"]),
                                )
                            ),
                        )
                        ctfa = sp_morphology.ctf_img_real(Tracker["nxinit"], ctf_params)
                        ctfs = ctfa
                    ##if( ( Blockdata["myid"] == Blockdata["main_node"])   and firsti ):
                    ##	dataimage.set_attr("is_complex",0)
                    ##	dataimage.write_image("dataimagefft.hdf")
                    ##	dataimage.set_attr("is_complex",1)
                    dataml = sp_fundamentals.fdecimate(
                        dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False, False
                    )
                    data = []
                    for iq in coarse_shifts:
                        xx = iq[0] * shrink
                        yy = iq[1] * shrink
                        dss = sp_fundamentals.fshift(dataml, xx, yy)
                        dss.set_attr("is_complex", 0)
                        data.append(dss)

                    #  This will get it to real space
                    # dataimage = fpol(Util.mulnclreal(Util.mulnclreal(fdecimate(dataimage, Tracker["nxinit"], Tracker["nxinit"], 1, False), Util.muln_img(bckgn, ctfs)), mask ), Tracker["nxpolar"], Tracker["nxpolar"],1, True)
                    #  Here we can reuse dataml to save time.  It was not normalized, but it does not matter as dataimage is for POLAR.  PAP 08/06/2018
                    dataimage = sp_fundamentals.fpol(
                        EMAN2_cppwrap.Util.mulnclreal(
                            EMAN2_cppwrap.Util.mulnclreal(
                                dataml, EMAN2_cppwrap.Util.muln_img(bckgn, ctfs)
                            ),
                            mask,
                        ),
                        Tracker["nxpolar"],
                        Tracker["nxpolar"],
                        1,
                        True,
                    )
                    # Compute max number of angles on the fly
                    lang = len(assignments_of_refangles_to_angles[im])
                    ###print("   BICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)

                if (
                    (Blockdata["myid"] == Blockdata["main_node"])
                    and (lima % (max(1, old_div(nima, 5))) == 0)
                    and (lima > 0)
                ):
                    sp_global_def.sxprint(
                        "  Number of images :%7d   %5d  %5.1f"
                        % (lima, nima, old_div(float(lima), float(nima)) * 100.0)
                        + "%"
                        + "   %10.1fmin" % (old_div((time.time() - at), 60.0))
                    )
                    ##print( "  Number of images :%7d   %5d  %5.1f"%(lima,nima,float(lima)/float(nima)*100.) + "%" +"   %10.1fmin   %10.1fmin"%((time()-at)/60.,eat/60.0))
                lima += 1

                ###print("  CONA1    ",Blockdata["myid"],lima)
                if lima == 1 and procid == 0:
                    ###print("  CONA2    ",Blockdata["myid"])
                    if lenass > 0:
                        ###print("   CICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)
                        keepfirst = (
                            lang * m_coarse_psi * n_coarse_shifts
                        )  # 150#lang*m_coarse_psi*len(coarse_shifts_shrank)  #500#min(200,nlocal_angles)#min(max(lxod1/25,200),lxod1)
                        ###if( Blockdata["myid"] == 18 and lima<5):  sxprint(" START   nlocal_angles* m_coarse_psi*len(coarse_shifts_shrank",nlocal_angles,coarse_delta,len(coarse_shifts_shrank),keepfirst)
                        # if( Blockdata["myid"] == Blockdata["main_node"] ):
                        lxod1 = EMAN2_cppwrap.Util.multiref_Crosrng_msg_stack_stepsi_local(
                            dataimage,
                            bigbuffer,
                            coarse_shifts_shrank,
                            assignments_of_refangles_to_angles[im],
                            assignments_of_refangles_to_cones[icone],
                            numr,
                            [coarse_angles[i][2] for i in range(n_coarse_ang)],
                            oldparams[im][2],
                            c_coarse_psi,
                            coarse_delta,
                            cnx,
                            keepfirst,
                        )

                        ##'''
                        assert old_div(len(lxod1), 3) == keepfirst

                        xod1 = sp_helix_sphire.np.ndarray(
                            (keepfirst), dtype="f4", order="C"
                        )
                        # xod1.fill(1.0)
                        xod2 = sp_helix_sphire.np.ndarray(
                            (keepfirst), dtype="int", order="C"
                        )
                        for iq in range(keepfirst):
                            ioffset = 3 * iq
                            #          ishift         iang                      ipsi
                            xod2[iq] = (
                                lxod1[ioffset]
                                + lxod1[ioffset + 1] * 100000000
                                + lxod1[ioffset + 2] * 1000
                            )

                        ##'''
                        # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

                        #  Second step - find which coarse ones are significant

                        # DO NOT order by angular directions to save time on reprojections.

                        pre_ipsiandiang = -1
                        for iln in range(keepfirst):
                            hashparams = int(xod2[iln])
                            ishift = hashparams % 1000
                            ipsiandiang = old_div(hashparams, 1000)
                            if ipsiandiang != pre_ipsiandiang:
                                pre_ipsiandiang = ipsiandiang
                                ipsi = ipsiandiang % 100000
                                iang = old_div(ipsiandiang, 100000)
                                ##junk = time()
                                temp = sp_projection.prgl(
                                    volinit,
                                    [
                                        coarse_angles[iang][0],
                                        coarse_angles[iang][1],
                                        (coarse_angles[iang][2] + ipsi * coarse_delta)
                                        % 360.0,
                                        0.0,
                                        0.0,
                                    ],
                                    1,
                                    False,
                                )
                                ##eat += time()-junk
                                ###   if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  SELECTEDSTEPTWO  ",iln,refang[iang][0],refang[iang][1],refang[iang][2]+ipsi*Tracker["delta"],ishift,xod1[iln])
                                temp.set_attr("is_complex", 0)
                            ##junk = time()
                            xod1[iln] = -EMAN2_cppwrap.Util.sqed(
                                data[ishift], temp, ctfa, bckgnoise
                            )
                            ##eat += time()-junk
                            ##xod2[iln] = hashparams

                        xod1 -= sp_helix_sphire.np.max(xod1)
                        lina = sp_helix_sphire.np.argwhere(
                            xod1 > Tracker["constants"]["expthreshold"]
                        )
                        # if( Blockdata["myid"] == 0 ):
                        ###print("  STARTING1   ",Blockdata["myid"],np.max(xod1),np.min(xod1),len(lina),lina[-1])

                        lina = lina.reshape(lina.size)
                        keepf = int(lina[-1]) + 1

                        xod1 = xod1[lina]
                        xod2 = xod2[lina]

                        ###print("  STARTING2    ",Blockdata["myid"])

                        lina = sp_helix_sphire.np.argsort(xod1)
                        xod1 = xod1[lina[::-1]]  # This sorts in reverse order
                        xod2 = xod2[lina[::-1]]  # This sorts in reverse order
                        sp_helix_sphire.np.exp(xod1, out=xod1)
                        xod1 = old_div(xod1, sp_helix_sphire.np.sum(xod1))
                        cumprob = 0.0
                        lit = len(xod1)
                        ###print("  STARTING3    ",Blockdata["myid"],lit)
                        for j in range(len(xod1)):
                            cumprob += xod1[j]
                            if cumprob > Tracker["constants"]["ccfpercentage"]:
                                lit = j + 1
                                break

                        ###print("  STARTING4    ",Blockdata["myid"],lit,keepf)
                    #  Turn into percentage of all possible
                    keepf = [int(old_div(float(keepf * 100), float(keepfirst)))]
                    ###mpi_barrier(MPI_COMM_WORLD)
                    ###print("  STARTING5    ",Blockdata["myid"],keepf,MPI_COMM_WORLD)
                    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
                    keepf = sp_utilities.wrap_mpi_gatherv(
                        keepf, Blockdata["main_node"], mpi.MPI_COMM_WORLD
                    )
                    ###print("  STARTING6    ",Blockdata["myid"],keepf)
                    if Blockdata["myid"] == 0:
                        keepf = [junk for junk in keepf if junk > 0]
                        if len(keepf) < 2:
                            keepf = 0
                        else:
                            keepf.sort()
                            keepf = keepf[int(len(keepf) * Blockdata["rkeepf"])]
                    else:
                        keepf = 0
                    ###print("  STARTING7    ",Blockdata["myid"],keepf)
                    keepf = sp_utilities.wrap_mpi_bcast(
                        keepf, Blockdata["main_node"], mpi.MPI_COMM_WORLD
                    )
                    if keepf == 0:
                        keepf = 2
                    ###print("  STARTING8    ",Blockdata["myid"],keepf)
                    Tracker["keepfirst"] = int(keepf)
                    ###if( Blockdata["myid"] == 0 ):  sxprint("  keepfirst first ",Tracker["keepfirst"])

                else:
                    if lenass > 0:
                        ###print("   DICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im,lang)#,assignments_to_cones)
                        ###if( Blockdata["myid"] == 18 and lima<5):  sxprint(" START   nlocal_angles* m_coarse_psi*len(coarse_shifts_shrank",nlocal_angles,coarse_delta,len(coarse_shifts_shrank),Tracker["keepfirst"])
                        # if( Blockdata["myid"] == Blockdata["main_node"] ):
                        keepfirst = lang * m_coarse_psi * n_coarse_shifts
                        keepfirst = max(
                            old_div(keepfirst * Tracker["keepfirst"], 100),
                            min(keepfirst, 3),
                        )
                        lxod1 = EMAN2_cppwrap.Util.multiref_Crosrng_msg_stack_stepsi_local(
                            dataimage,
                            bigbuffer,
                            coarse_shifts_shrank,
                            assignments_of_refangles_to_angles[im],
                            assignments_of_refangles_to_cones[icone],
                            numr,
                            [coarse_angles[i][2] for i in range(n_coarse_ang)],
                            oldparams[im][2],
                            c_coarse_psi,
                            coarse_delta,
                            cnx,
                            keepfirst,
                        )
                        # if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  ")
                        ##'''
                        assert keepfirst == old_div(len(lxod1), 3)
                        xod1 = sp_helix_sphire.np.ndarray(
                            (keepfirst), dtype="f4", order="C"
                        )
                        # xod1.fill(1.0)
                        xod2 = sp_helix_sphire.np.ndarray(
                            (keepfirst), dtype="int", order="C"
                        )
                        for iq in range(keepfirst):
                            ioffset = 3 * iq
                            #          ishift         iang                      ipsi
                            xod2[iq] = (
                                lxod1[ioffset]
                                + lxod1[ioffset + 1] * 100000000
                                + lxod1[ioffset + 2] * 1000
                            )

                        ##'''

                        # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

                        #  Second step - find which coarse ones are significant

                        # order by angular directions to save time on reprojections.
                        ipsiandiang = old_div(xod2, 1000)
                        lina = sp_helix_sphire.np.argsort(ipsiandiang)
                        xod2 = xod2[lina]  # order does not matter

                        pre_ipsiandiang = -1
                        for iln in range(keepfirst):
                            hashparams = int(xod2[iln])
                            ishift = hashparams % 1000
                            ipsiandiang = old_div(hashparams, 1000)
                            if ipsiandiang != pre_ipsiandiang:
                                pre_ipsiandiang = ipsiandiang
                                ipsi = ipsiandiang % 100000
                                iang = old_div(ipsiandiang, 100000)
                                ##junk = time()
                                temp = sp_projection.prgl(
                                    volinit,
                                    [
                                        coarse_angles[iang][0],
                                        coarse_angles[iang][1],
                                        (coarse_angles[iang][2] + ipsi * coarse_delta)
                                        % 360.0,
                                        0.0,
                                        0.0,
                                    ],
                                    1,
                                    False,
                                )
                                ###   if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  SELECTEDSTEPTWO  ",iln,refang[iang][0],refang[iang][1],refang[iang][2]+ipsi*Tracker["delta"],ishift,xod1[iln])
                                ##eat += time()-junk
                                temp.set_attr("is_complex", 0)
                            ##junk = time()
                            peak = -EMAN2_cppwrap.Util.sqed(
                                data[ishift], temp, ctfa, bckgnoise
                            )
                            ##eat += time()-junk
                            #  Note I replace ccc by eqdist
                            xod1[iln] = peak
                            ##xod2[iln] = hashparams

                        lina = sp_helix_sphire.np.argsort(xod1)
                        xod1 = xod1[lina[::-1]]  # This sorts in reverse order
                        xod2 = xod2[lina[::-1]]  # This sorts in reverse order

                        xod1 -= xod1[0]

                        # if( Blockdata["myid"] == Blockdata["main_node"]):
                        # 	#print("  PROJECT   ",im,lit,johi)#,cod2)
                        # 	for iln in range(len(xod1)):  sxprint("  ROPE   ",iln,xod1[iln],xod2[iln])

                        lina = sp_helix_sphire.np.argwhere(
                            xod1 > Tracker["constants"]["expthreshold"]
                        )
                        xod1 = xod1[lina]
                        xod2 = xod2[lina]
                        sp_helix_sphire.np.exp(xod1, out=xod1)
                        xod1 = old_div(xod1, sp_helix_sphire.np.sum(xod1))
                        cumprob = 0.0
                        lit = len(xod1)
                        for j in range(len(xod1)):
                            cumprob += xod1[j]
                            if cumprob > Tracker["constants"]["ccfpercentage"]:
                                lit = j + 1
                                break

                        #  To here
                        ###if( Blockdata["myid"] == 18 and lima<5):  sxprint("  SECOND KEPT  ",lit)
                        # if( lima<5):  sxprint("  SECOND KEPT  ",lit)

                # mpi_barrier(MPI_COMM_WORLD)
                # mpi_finalize()
                # exit()
                # for j in range(lit):
                # 	 newpar[kl][2].append([int(xod2[j]),float(xod1[j])])
                if lenass > 0:
                    ###print("   EICONE icnm,im in enumerateassignments_to_cones[icone]  ",Blockdata["myid"],icone,icnm,im)#,assignments_to_cones)

                    firstdirections = [[0.0, 0.0] for iln in range(lit)]
                    firstshifts = [0] * lit
                    for iln in range(lit):
                        hashparams = int(xod2[iln])
                        ishift = hashparams % 1000
                        ipsiandiang = old_div(hashparams, 1000)
                        # ipsi = ipsiandiang%100000
                        iang = old_div(ipsiandiang, 100000)
                        # try:
                        firstdirections[iln] = [
                            coarse_angles[iang][0],
                            coarse_angles[iang][1],
                            0.0,
                        ]
                        # except:
                        # 	sxprint(" FAILED  ",Blockdata["myid"],Tracker["keepfirst"],iln,lima,lit,hashparams,iang,xod2[:lit],xod1[:lit])
                        # 	mpi_finalize()
                        # 	exit()
                        firstshifts[iln] = ishift
                        #  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):
                        # 	ipsi = ipsiandiang%100000
                        # 	ishift = hashparams%1000
                        # 	###print("  SECONDPAR%04d  "%im,iln,hashparams, ipsi,iang,ishift)
                        # 	sxprint(" SECONDPAR%04d  "%im,iln, refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, shifts[ishift])
                    ###del xod2
                    ###if( Blockdata["myid"] == 18 and lima<5):   sxprint("  SECOND KEPT  ",im,lit)#,len(xod2),xod2,firstdirections,firstshifts)
                    # mpi_barrier(MPI_COMM_WORLD)
                    # mpi_finalize()
                    # exit()

                    # if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  FIFI ",firstdirections)
                    # if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  GUGU ",firstshifts)
                    # Find neighbors, ltabang contains positions on refang list, no psis
                    ###ltabang = nearest_many_full_k_projangles(refang, firstdirections, howmany = 5, sym=Tracker["constants"]["symmetry"])
                    ltabang = find_nearest_k_refangles_to_many_angles(
                        refdirs,
                        firstdirections,
                        Tracker["delta"],
                        howmany=Tracker["howmany"],
                    )
                    ###if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  ltabang ",ltabang)
                    ##mpi_barrier(MPI_COMM_WORLD)
                    ##mpi_finalize()
                    ##exit()

                    # ltabang has length lit, which is the same as length as firshifts.  However, original xod2 is still available,
                    #   even though it is longer than lit.

                    ###ltabshi = [shiftneighbors[iln] for iln in firstshifts]
                    # if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  HUHU ",ltabang)
                    # if(Blockdata["myid"] == Blockdata["main_node"]):  sxprint("  OUOU ",ltabshi)

                    #  Prepare image for chi2.
                    #  We have to repeat everything from get shrink data, including shifts
                    #  The number of entries is lit=len(ltabang), each proj dir has 4 neighbors, so including itself there are 5,
                    #    there are 3 psis, and at most n_fine_shifts.
                    #    Not not all shifts have to be populated. We may also have repeated triplets of angles, but each can have
                    #      different shifts.  If so, we have to remove duplicates from the entire set.
                    lcod1 = lit * 4 * 2 * n_fine_shifts

                    # cod2 = np.ndarray((lit,5,3,n_fine_shifts),dtype=int,order="C")
                    # cod2.fill(-1)  #  hashparams
                    cod2 = []
                    # lol = 0
                    for i1 in range(lit):
                        hashparams = int(xod2[i1])
                        ipsiandiang = old_div(hashparams, 1000)
                        oldiang = old_div(ipsiandiang, 100000)
                        ipsi = ipsiandiang % 100000
                        ishift = hashparams % 1000
                        tshifts = get_shifts_neighbors(shifts, coarse_shifts[ishift])
                        # if( Blockdata["myid"] == Blockdata["main_node"]): sxprint(" tshifts  ",i1,len(tshifts))
                        for i2 in range(Tracker["howmany"]):
                            iang = ltabang[i1][i2]
                            for i3 in range(2):  # psi
                                itpsi = int(
                                    old_div(
                                        (
                                            coarse_angles[oldiang][2]
                                            + ipsi * coarse_delta
                                            - refang[iang][2]
                                            + 360.0
                                        ),
                                        Tracker["delta"],
                                    )
                                )
                                itpsi = (itpsi + i3) % npsi
                                for i4 in range(len(tshifts)):
                                    cod2.append(
                                        iang * 100000000 + itpsi * 1000 + tshifts[i4]
                                    )
                                    # lol += 1
                                    # if( Blockdata["myid"] == Blockdata["main_node"] ):  sxprint("  zibzi  ",i1,i2,i3,i4, lol)

                    del xod1, xod2

                    ###if( Blockdata["myid"] == 18 and lima<5):   sxprint("  THIRD   ",len(cod2))#,cod2)
                    cod2 = list(set(cod2))
                    cod1 = [[old_div(q, 1000), i] for i, q in enumerate(cod2)]
                    cod1.sort()

                    lit = len(cod1)

                    cod2 = sp_helix_sphire.np.asarray(
                        [cod2[cod1[i][1]] for i in range(lit)]
                    )

                    cod1 = sp_helix_sphire.np.ndarray(lit, dtype="f4", order="C")
                    # cod1.fill(np.finfo(dtype='f4').min)
                    cod3 = sp_helix_sphire.np.ndarray(lit, dtype="f4", order="C")
                    # cod3.fill(0.0)  #  varadj

                    ###if( Blockdata["myid"] == 18 and lima<5): sxprint("  THIRD   ",im,lit)#,cod2)

                    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    #  Third ML part.  Now the actual chi2 distances, we compute reprojections on the fly.
                    #  Make sure volprep has nxinit size
                    data = [None] * nshifts
                    johi = 0
                    iln = 0
                    prevdir = -1
                    while iln < lit:
                        hashparams = cod2[iln]
                        # if( Blockdata["myid"] == Blockdata["main_node"]): sxprint("  COD2   ",im,lit,iln,cod2[iln])
                        ipsiandiang = old_div(hashparams, 1000)
                        if ipsiandiang != prevdir:
                            prevdir = ipsiandiang
                            ipsi = ipsiandiang % 100000
                            iang = old_div(ipsiandiang, 100000)
                            ##junk = time()
                            temp = sp_projection.prgl(
                                volinit,
                                [
                                    refang[iang][0],
                                    refang[iang][1],
                                    (refang[iang][2] + ipsi * Tracker["delta"]) % 360.0,
                                    0.0,
                                    0.0,
                                ],
                                1,
                                False,
                            )
                            ##eat += time()-junk
                            temp.set_attr("is_complex", 0)
                            johi += 1
                        while ipsiandiang == old_div(cod2[iln], 1000):
                            hashparams = cod2[iln]
                            ishift = hashparams % 1000
                            if data[ishift] == None:
                                xx = shifts[ishift][0] * shrink
                                yy = shifts[ishift][1] * shrink
                                data[ishift] = sp_fundamentals.fshift(dataml, xx, yy)
                                data[ishift].set_attr("is_complex", 0)
                            ##junk = time()
                            [peak, varadj] = EMAN2_cppwrap.Util.sqednorm(
                                data[ishift], temp, ctfa, bckgnoise
                            )
                            ##eat += time()-junk
                            cod1[iln] = -peak
                            cod3[iln] = varadj
                            iln += 1
                            if iln == lit:
                                break
                        # if( Blockdata["myid"] == Blockdata["main_node"]):
                        # 	temp.write_image("temp.hdf")
                        # 	data[iln].write_image("data.hdf")
                        # 	ctfa.write_image("ctfa.hdf")
                        # 	bckgnoise.write_image("bckgnoise.hdf")
                        # 	exit()
                        ###if( Blockdata["myid"] == Blockdata["main_node"] and iln%1000 ==0):  sxprint(" progress  ",iln,time()-at)
                    # if( Blockdata["myid"] == Blockdata["main_node"]):
                    # 	sxprint("  PROJECT   ",im,lit,johi)#,cod2)
                    # 	#for iln in range(lit):  sxprint("  ROPE   ",iln,cod1[iln],cod2[iln],cod3[iln])
                    del data
                    del dataml

                    lina = sp_helix_sphire.np.argsort(cod1)
                    cod1 = cod1[lina[::-1]]  # This sorts in reverse order
                    cod2 = cod2[lina[::-1]]  # This sorts in reverse order
                    cod3 = cod3[lina[::-1]]  # This sorts in reverse order
                    cod1 -= cod1[0]
                    lina = sp_helix_sphire.np.argwhere(
                        cod1 > Tracker["constants"]["expthreshold"]
                    )
                    cod1 = cod1[lina]
                    cod2 = cod2[lina]
                    cod3 = cod3[lina]

                    ###if( Blockdata["myid"] == Blockdata["main_node"]):
                    ###for iui in range(len(lina)):
                    ###	for iui in range(len(cod1)):
                    ###		sxprint("  MLML  ",iui,cod1[iui],exp(cod1[iui]),cod2[iui],cod3[iui])

                    sp_helix_sphire.np.exp(cod1, out=cod1)
                    cod1 = old_div(cod1, sp_helix_sphire.np.sum(cod1))
                    cumprob = 0.0
                    for j in range(len(cod1)):
                        cumprob += cod1[j]
                        if cumprob > Tracker["constants"]["ccfpercentage"]:
                            lit = j + 1
                            break

                    #  New norm is a sum of eq distances multiplied by their probabilities augmented by PW.
                    norm_per_particle[im] = (
                        sp_helix_sphire.np.sum(cod1[:lit] * cod3[:lit])
                        + accumulatepw[im][reachpw]
                    )
                    ###print("   CNORMPERPARTICLE  ",Blockdata["myid"],im,norm_per_particle[im])

                    for iln in range(lit):
                        newpar[im][2].append([int(cod2[iln]), float(cod1[iln])])
                        #  CONTROLPRINTOUT   if( Blockdata["myid"] == Blockdata["main_node"] and icnm<5):
                        # 	#print("  NEWPAR%04d  "%im,iln,newpar[im][2][-1])
                        # 	hashparams = newpar[im][2][iln][0]
                        # 	ipsiandiang	= hashparams/1000
                        # 	ipsi = ipsiandiang%100000
                        # 	iang = ipsiandiang/100000
                        # 	ishift = hashparams%1000
                        # 	#print("  NEWPAR%04d  "%im,iln,newpar[im][2][-1],hashparams,ipsi,iang,ishift)
                        # 	sxprint(" NEWPAR%04d  "%im,iln,newpar[im][2][iln], refang[iang][0],refang[iang][1],(refang[iang][2] + ipsi*Tracker["delta"])%360.0, shifts[ishift])

                    ###if( Blockdata["myid"] == 18 and lima<5):   sxprint("  FINALLY  ",im,lit)
                    del cod1, cod2, cod3, lina
                    ###mpi_barrier(MPI_COMM_WORLD)
                    ###mpi_finalize()
                    ###exit()

    """Multiline Comment22"""

    #  END OF CONES
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_global_def.sxprint(
            "  Finished projection matching   %10.1fmin"
            % (old_div((time.time() - at), 60.0))
        )
    at = time.time()
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    #  All images were processed, now to the additional calculations

    ###mpi_barrier(MPI_COMM_WORLD)
    ###mpi_finalize()
    ###exit()

    # norm correction ---- calc the norm correction per particle
    snormcorr = 0.0
    for kl in range(nima):
        ###print("   NORMPERPARTICLE  ",Blockdata["myid"],kl,norm_per_particle[kl])
        norm_per_particle[kl] = old_div(
            numpy.sqrt(norm_per_particle[kl] * 2.0) * oldparams[kl][7],
            Tracker["avgvaradj"][procid],
        )
        snormcorr += norm_per_particle[kl]
    Tracker["avgvaradj"][procid] = snormcorr
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    #  Compute avgvaradj
    Tracker["avgvaradj"][procid] = mpi.mpi_reduce(
        Tracker["avgvaradj"][procid],
        1,
        mpi.MPI_FLOAT,
        mpi.MPI_SUM,
        Blockdata["main_node"],
        mpi.MPI_COMM_WORLD,
    )
    if Blockdata["myid"] == Blockdata["main_node"]:
        Tracker["avgvaradj"][procid] = old_div(
            float(Tracker["avgvaradj"][procid]), Tracker["nima_per_chunk"][procid]
        )
    else:
        Tracker["avgvaradj"][procid] = 0.0
    Tracker["avgvaradj"][procid] = sp_utilities.bcast_number_to_all(
        Tracker["avgvaradj"][procid], Blockdata["main_node"]
    )
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    #  Compute statistics of smear -----------------
    smax = -1000000
    smin = 1000000
    sava = 0.0
    svar = 0.0
    snum = 0
    for kl in range(nima):
        j = len(newpar[kl][2])
        snum += 1
        sava += float(j)
        svar += j * float(j)
        smax = max(smax, j)
        smin = min(smin, j)
    snum = mpi.mpi_reduce(
        snum, 1, mpi.MPI_INT, mpi.MPI_SUM, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    sava = mpi.mpi_reduce(
        sava, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    svar = mpi.mpi_reduce(
        svar, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    smax = mpi.mpi_reduce(
        smax, 1, mpi.MPI_INT, mpi.MPI_MAX, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    smin = mpi.mpi_reduce(
        smin, 1, mpi.MPI_INT, mpi.MPI_MIN, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    if Blockdata["myid"] == 0:
        sava = old_div(float(sava), snum)
        svar = numpy.sqrt(
            max(0.0, old_div((float(svar) - snum * sava ** 2), (snum - 1)))
        )
        # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
        line = ""
        sp_global_def.sxprint(
            line,
            "Smear stat  (number of images, ave, sumsq, min, max)):  %7d    %12.3g   %12.3g  %7d  %7d"
            % (snum, sava, svar, smin, smax),
        )

    at = time.time()
    mpi.mpi_barrier(Blockdata["shared_comm"])

    ###if Blockdata["myid"] == Blockdata["main_node"]:  print "  Finished :",time()-at
    mpi.mpi_win_free(win_sm)
    if Tracker["nxinit"] != Tracker["nxpolar"]:
        mpi.mpi_win_free(win_volinit)
        emnumpy4.unregister_numpy_from_emdata()
        del emnumpy4
    else:
        mpi.mpi_win_free(win_vol)

    mpi.mpi_barrier(Blockdata["shared_comm"])
    emnumpy1.unregister_numpy_from_emdata()
    emnumpy2.unregister_numpy_from_emdata()
    del emnumpy1, emnumpy2

    del volinit

    mpi.mpi_barrier(Blockdata["shared_comm"])

    # print("  NORMALIZATION DONE  ",Blockdata["myid"])
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    # if( Blockdata["myid"] == Blockdata["main_node"] ):
    # 	sxprint( "  Projection matching finished : %10.1fmin"%((time()-at)/60.))
    return newpar, norm_per_particle


def cerrs(params, ctfs, particle_groups):
    global Tracker, Blockdata

    shrinkage = old_div(float(Tracker["nxinit"]), float(Tracker["constants"]["nnxo"]))
    procid = 0
    if Blockdata["myid"] == Blockdata["nodes"][procid]:
        ref_vol = sp_utilities.get_im(Tracker["refvol"])
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
    else:
        # log = None
        ref_vol = sp_utilities.model_blank(
            Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"]
        )
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    sp_utilities.bcast_EMData_to_all(
        ref_vol, Blockdata["myid"], Blockdata["nodes"][procid]
    )
    interpolation_method = 1
    ref_vol = sp_projection.prep_vol(
        ref_vol, npad=2, interpolation_method=interpolation_method
    )

    lb = Blockdata["bckgnoise"].get_xsize()
    acc_rot = 0.0
    acc_trans = 0.0

    # // P(X | X_1) / P(X | X_2) = exp ( |F_1 - F_2|^2 / (-2 sigma2) )
    # // exp(-4.60517) = 0.01
    pvalue = 4.60517

    for itry in range(len(params)):

        # // Get orientations (angles1) for this particle
        phi1 = params[itry][0]
        theta1 = params[itry][1]
        psi1 = params[itry][2]
        # // Get CTF for this particle
        #   Get F1 = Proj(refvol; angles1, shifts=0)
        F1 = sp_projection.prgl(
            ref_vol,
            [phi1, theta1, psi1, 0.0, 0.0],
            interpolation_method=1,
            return_real=False,
        )
        ctfs[itry].apix = old_div(ctfs[itry].apix, shrinkage)
        ct = sp_morphology.ctf_img_real(Tracker["nxinit"], ctfs[itry])
        EMAN2_cppwrap.Util.mul_img(ct, ct)
        ctfsbckgnoise = EMAN2_cppwrap.Util.muln_img(
            EMAN2_cppwrap.Util.unroll1dpw(
                Tracker["nxinit"],
                Tracker["nxinit"],
                [Blockdata["bckgnoise"][i, particle_groups[itry]] for i in range(lb)],
            ),
            ct,
        )

        # // Search 2 times: angles and shifts
        for imode in range(2):
            ang_error = 0.0
            sh_error = 0.0
            peak = 0.0

            # // Search for ang_error and sh_error where there are at least 3-sigma differences!
            while peak <= pvalue:
                # // Graduallly increase the step size
                if ang_error < 0.2:
                    ang_step = 0.05
                elif ang_error < 1.0:
                    ang_step = 0.1
                elif ang_error < 2.0:
                    ang_step = 0.2
                elif ang_error < 5.0:
                    ang_step = 0.5
                elif ang_error < 10.0:
                    ang_step = 1.0
                elif ang_error < 20.0:
                    ang_step = 2.0
                else:
                    ang_step = 5.0

                if sh_error < 0.2:
                    sh_step = 0.05
                elif sh_error < 1.0:
                    sh_step = 0.1
                elif sh_error < 2.0:
                    sh_step = 0.2
                elif sh_error < 5.0:
                    sh_step = 0.5
                elif sh_error < 10.0:
                    sh_step = 1.0
                else:
                    sh_step = 2.0

                ang_error += ang_step
                sh_error += sh_step

                # // Prevent an endless while by putting boundaries on ang_error and sh_error
                if (imode == 0 and ang_error > 30.0) or (
                    imode == 1 and sh_error > 10.0
                ):
                    break

                phi2 = phi1
                theta2 = theta1
                psi2 = psi1
                xoff1 = yoff1 = 0.0
                xshift = yshift = 0.0

                # // Perturb angle or shift , depending on the mode
                ran = random.random()
                if imode == 0:
                    if ran < 0.3333:
                        phi2 = phi1 + ang_error
                    elif ran < 0.6667:
                        theta2 = theta1 + ang_error
                    else:
                        psi2 = psi1 + ang_error
                else:
                    if ran < 0.5:
                        xshift = xoff1 + sh_error
                        yshift = 0.0
                    else:
                        xshift = 0.0
                        yshift = yoff1 + sh_error

                if imode == 0:
                    F2 = sp_projection.prgl(
                        ref_vol, [phi2, theta2, psi2, 0.0, 0.0], 1, False
                    )
                else:
                    F2 = sp_fundamentals.fshift(
                        F1, xshift * shrinkage, yshift * shrinkage
                    )

                peak = EMAN2_cppwrap.Util.sqedac(F1, F2, ctfsbckgnoise)

            if imode == 0:
                acc_rot += ang_error
            elif imode == 1:
                acc_trans += sh_error

    acc_rot = mpi.mpi_reduce(
        acc_rot,
        1,
        mpi.MPI_FLOAT,
        mpi.MPI_SUM,
        Blockdata["main_node"],
        mpi.MPI_COMM_WORLD,
    )
    acc_trans = mpi.mpi_reduce(
        acc_trans,
        1,
        mpi.MPI_FLOAT,
        mpi.MPI_SUM,
        Blockdata["main_node"],
        mpi.MPI_COMM_WORLD,
    )
    acc_rot = mpi.mpi_bcast(
        acc_rot, 1, mpi.MPI_FLOAT, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    acc_trans = mpi.mpi_bcast(
        acc_trans, 1, mpi.MPI_FLOAT, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )

    acc_rot = float(acc_rot[0])
    acc_trans = float(acc_trans[0])
    n_trials = Blockdata["nproc"] * len(params)

    acc_rot = old_div(acc_rot, n_trials)
    acc_trans = old_div(acc_trans, n_trials)

    if Blockdata["myid"] == Blockdata["main_node"]:
        # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
        line = ""
        sp_global_def.sxprint(
            line,
            "Estimated accuracy of angles = %6.2f degrees; and shifts = %5.1f pixels"
            % (acc_rot, acc_trans),
        )

    Tracker["acc_rot"] = acc_rot
    Tracker["acc_trans"] = acc_trans


def do3d_final(
    partids,
    partstack,
    original_data,
    oldparams,
    oldparamstructure,
    projdata,
    final_iter=-1,
    comm=-1,
):
    global Tracker, Blockdata

    final_dir = Tracker["directory"]
    if Blockdata["subgroup_myid"] > -1:
        this_color = 1
    else:
        this_color = 0
    subgroup_comm = mpi.mpi_comm_split(
        mpi.MPI_COMM_WORLD, this_color, Blockdata["myid"]
    )
    if Blockdata["subgroup_myid"] > -1:
        # load datastructure, read data, do two reconstructions(stepone, steptwo)
        if final_iter == -1:
            final_iter = Tracker["constants"]["best"]
        carryon = 1
        # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
        line = ""
        if Blockdata["subgroup_myid"] == Blockdata["main_node"]:
            sp_global_def.sxprint(line, "do_final_rec3d")
            sp_global_def.sxprint(
                "Reconstruction uses solution of %d iteration" % final_iter
            )
            sp_global_def.sxprint(
                "Final reconstruction image size is:  %d"
                % (Tracker["constants"]["nnxo"])
            )
            sp_global_def.sxprint("Final directory is %s" % (Tracker["directory"]))
        final_dir = Tracker["directory"]
        if Blockdata["subgroup_myid"] == Blockdata["main_node"]:
            try:
                refang = sp_utilities.read_text_row(
                    os.path.join(final_dir, "refang.txt")
                )
                rshifts = sp_utilities.read_text_row(
                    os.path.join(final_dir, "rshifts.txt")
                )
            except:
                carryon = 0
        else:
            refang = 0
            rshifts = 0
        carryon = sp_utilities.bcast_number_to_all(
            carryon, source_node=Blockdata["main_node"], mpi_comm=comm
        )
        if carryon == 0:
            sp_global_def.ERROR(
                "Failed to read refang and rshifts: %s %s "
                % (
                    os.path.join(final_dir, "refang.txt"),
                    os.path.join(final_dir, "rshifts.txt"),
                ),
                "do_final_rec3d",
                1,
                Blockdata["subgroup_myid"],
            )
        refang = sp_utilities.wrap_mpi_bcast(refang, Blockdata["main_node"], comm)
        rshifts = sp_utilities.wrap_mpi_bcast(rshifts, Blockdata["main_node"], comm)

        partids = [None, None]
        if Blockdata["subgroup_myid"] == Blockdata["main_node"]:
            if not os.path.exists(
                os.path.join(Tracker["constants"]["masterdir"], "tempdir")
            ):
                os.makedirs(os.path.join(Tracker["constants"]["masterdir"], "tempdir"))
            l = 0
            for procid in range(2):
                partids[procid] = os.path.join(
                    final_dir,
                    "chunk_%01d_%03d.txt" % (procid, Tracker["mainiteration"]),
                )
                l += len(sp_utilities.read_text_file(partids[procid]))
        else:
            l = 0
        l = sp_utilities.bcast_number_to_all(
            l, source_node=Blockdata["main_node"], mpi_comm=comm
        )

        norm_per_particle = [[], []]
        # get the previous number of CPUs
        nproc_previous = 0
        if Blockdata["subgroup_myid"] == 0:
            while os.path.exists(
                os.path.join(
                    final_dir,
                    "oldparamstructure",
                    "oldparamstructure_%01d_%03d_%03d.json"
                    % (procid, nproc_previous, Tracker["mainiteration"]),
                )
            ):
                nproc_previous += 1
        nproc_previous = sp_utilities.bcast_number_to_all(
            nproc_previous, source_node=Blockdata["main_node"], mpi_comm=comm
        )

        for procid in range(2):
            if procid == 0:
                original_data[1] = None
            partids[procid] = os.path.join(
                final_dir, "chunk_%01d_%03d.txt" % (procid, Tracker["mainiteration"])
            )
            partstack[procid] = os.path.join(
                final_dir,
                "params-chunk_%01d_%03d.txt" % (procid, (Tracker["mainiteration"])),
            )
            ###
            psize = len(sp_utilities.read_text_file(partids[procid]))
            oldparamstructure[procid] = []
            im_start, im_end = sp_applications.MPI_start_end(
                psize, Blockdata["subgroup_size"], Blockdata["subgroup_myid"]
            )
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
            for iproc_index_old in range(istart_old_proc_id, iend_old_proc_id + 1):
                oldparamstructure_on_old_cpu = load_object_from_json(
                    os.path.join(
                        final_dir,
                        "oldparamstructure",
                        "oldparamstructure_%01d_%03d_%03d.json"
                        % (procid, iproc_index_old, Tracker["mainiteration"]),
                    )
                )
                mlocal_id_on_old = ptl_on_this_cpu - plist[iproc_index_old][0]
                while (mlocal_id_on_old < len(oldparamstructure_on_old_cpu)) and (
                    ptl_on_this_cpu < im_end
                ):
                    oldparamstructure[procid].append(
                        oldparamstructure_on_old_cpu[mlocal_id_on_old]
                    )
                    ptl_on_this_cpu += 1
                    mlocal_id_on_old += 1
            del oldparamstructure_on_old_cpu
            mpi.mpi_barrier(Blockdata["subgroup_comm"])
            #####
            original_data[procid], oldparams[procid], start_end = getindexdata(
                partids[procid],
                partstack[procid],
                os.path.join(
                    Tracker["constants"]["masterdir"],
                    "main000",
                    "particle_groups_%01d.txt" % procid,
                ),
                original_data[procid],
                small_memory=Tracker["constants"]["small_memory"],
                nproc=Blockdata["subgroup_size"],
                myid=Blockdata["subgroup_myid"],
                mpi_comm=comm,
            )
            temp = Tracker["directory"]
            Tracker["directory"] = os.path.join(
                Tracker["constants"]["masterdir"], "tempdir"
            )
            mpi.mpi_barrier(Blockdata["subgroup_comm"])
            if procid == 0:
                compute_sigma(
                    [[]] * l,
                    [[]] * l,
                    len(oldparams[0]),
                    True,
                    myid=Blockdata["subgroup_myid"],
                    mpi_comm=comm,
                )
            Tracker["directory"] = temp
            mpi.mpi_barrier(Blockdata["subgroup_comm"])
            projdata[procid] = get_shrink_data(
                Tracker["constants"]["nnxo"],
                procid,
                original_data[procid],
                oldparams[procid],
                return_real=False,
                preshift=True,
                apply_mask=False,
                nonorm=True,
            )
            for ipar in range(len(oldparams[procid])):
                norm_per_particle[procid].append(oldparams[procid][ipar][7])
            oldparams[procid] = []
            original_data[procid] = None
            # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
            line = ""
            if Blockdata["subgroup_myid"] == Blockdata["nodes"][procid]:
                sp_global_def.sxprint(line, "3-D reconstruction of group %d" % procid)
            ###--------------------------------------------------------- Force
            Tracker["directory"] = Tracker["constants"]["masterdir"]
            Tracker["nxinit"] = Tracker["constants"]["nnxo"]
            Tracker["maxfrad"] = old_div(Tracker["constants"]["nnxo"], 2)
            ###---------------------------------------------------------

            outlier_file = os.path.join(
                final_dir,
                "outlier-params-chunk_%01d_%03d.txt"
                % (procid, (Tracker["mainiteration"])),
            )

            if Tracker["prior"]["apply_prior"]:
                outlier_list = sp_utilities.read_text_file(outlier_file)[
                    start_end[0] : start_end[1]
                ]
            else:
                outlier_list = [0] * len(projdata[procid])
            norm_per_particle_outlier = []
            oldparamstructure_outlier = []
            projdata_outlier = []
            for idx, entry in enumerate(outlier_list):
                if int(entry) == 0:
                    projdata_outlier.append(projdata[procid][idx])
                    norm_per_particle_outlier.append(norm_per_particle[procid][idx])
                    oldparamstructure_outlier.append(oldparamstructure[procid][idx])
            total_left_particles = sp_utilities.wrap_mpi_gatherv(
                norm_per_particle_outlier, Blockdata["main_node"], subgroup_comm
            )
            if Blockdata["myid"] == Blockdata["main_node"]:
                sp_global_def.sxprint(
                    "Use {0} particles for final reconstruction!".format(
                        len(total_left_particles)
                    )
                )

            do3d(
                procid,
                projdata_outlier,
                oldparamstructure_outlier,
                refang,
                rshifts,
                norm_per_particle_outlier,
                myid=Blockdata["subgroup_myid"],
                smearing=True,
                mpi_comm=comm,
            )
            projdata[procid] = []
            oldparamstructure[procid] = []
            norm_per_particle[procid] = []
            mpi.mpi_barrier(Blockdata["subgroup_comm"])
        mpi.mpi_barrier(Blockdata["subgroup_comm"])
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
    line = ""
    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_global_def.sxprint(line, "final rec3d_make_maps")
    rec3d_make_maps(compute_fsc=False, regularized=False)

    # also copy params to masterdir as final params
    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_helix_sphire.shutil.copyfile(
            os.path.join(
                Tracker["constants"]["masterdir"],
                "main%03d" % Tracker["mainiteration"],
                "params_%03d.txt" % Tracker["mainiteration"],
            ),
            os.path.join(
                Tracker["constants"]["masterdir"],
                "final_params_%03d.txt" % Tracker["mainiteration"],
            ),
        )
        sp_helix_sphire.shutil.rmtree(
            os.path.join(Tracker["constants"]["masterdir"], "tempdir")
        )
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    return


def recons3d_final(masterdir, do_final_iter_init, memory_per_node, orgstack=None):
    global Tracker, Blockdata
    # search for best solution, and load respective tracker
    carryon = 1
    # line     = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
    line = ""
    if Blockdata["myid"] == Blockdata["main_node"]:
        print(line, "recons3d_final")
    do_final_iter = 3
    if do_final_iter_init == 0:
        if Blockdata["myid"] == Blockdata["main_node"]:
            try:
                Tracker = load_tracker_from_json(
                    os.path.join(masterdir, "Tracker_final.json")
                )
                sp_global_def.sxprint(
                    "The best solution is %d  " % Tracker["constants"]["best"]
                )
                do_final_iter = Tracker["constants"][
                    "best"
                ]  # set the best as do_final iteration
            except:
                carryon = 0
        carryon = sp_utilities.bcast_number_to_all(carryon)
        if carryon == 0:
            sp_global_def.ERROR(
                "Best resolution is not found, do_final will not be computed",
                "recons3d_final",
                1,
                Blockdata["myid"],
            )  # Now work on selected directory
        do_final_iter = sp_utilities.bcast_number_to_all(do_final_iter)
    elif do_final_iter_init == -1:
        do_final_iter = Tracker["constants"]["best"]
    else:
        do_final_iter = do_final_iter_init
        if Blockdata["myid"] == Blockdata["main_node"]:
            sp_global_def.sxprint(
                "User selected %d iteration to compute the 3D reconstruction "
                % do_final_iter
            )
        if do_final_iter <= 2:
            sp_global_def.ERROR(
                "The selected iteration should be larger than 2",
                "recons3d_final",
                1,
                Blockdata["myid"],
            )

    final_dir = os.path.join(masterdir, "main%03d" % do_final_iter)
    if Blockdata["myid"] == Blockdata["main_node"]:  # check json file and load tracker
        try:
            Tracker = load_tracker_from_json(
                os.path.join(final_dir, "Tracker_%03d.json" % do_final_iter)
            )
            stack_prior = Tracker["constants"]["stack_prior"]
            Tracker["constants"]["stack_prior"] = None
        except:
            carryon = 0
        if orgstack:
            Tracker["constants"]["stack"] = orgstack
    else:
        Tracker = 0
        stack_prior = None
    carryon = sp_utilities.bcast_number_to_all(
        carryon, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    if carryon == 0:
        sp_global_def.ERROR(
            "Failed to load Tracker file %s, program terminates "
            % os.path.join(final_dir, "Tracker_%03d.json" % do_final_iter),
            "recons3d_final",
            1,
            Blockdata["myid"],
        )
    Tracker = sp_utilities.wrap_mpi_bcast(
        Tracker, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    Tracker["constants"]["stack_prior"] = stack_prior
    if Blockdata["myid"] == Blockdata["main_node"]:  # check stack
        try:
            image = sp_utilities.get_im(Tracker["constants"]["stack"], 0)
        except:
            carryon = 0
    carryon = sp_utilities.bcast_number_to_all(
        carryon, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    if carryon == 0:
        sp_global_def.ERROR(
            "The orignal data stack for reconstruction %s does not exist, final reconstruction terminates"
            % Tracker["constants"]["stack"],
            "recons3d_final",
            1,
            Blockdata["myid"],
        )

    if Blockdata["myid"] == Blockdata["main_node"]:
        #  Estimated volume size
        volume_size = old_div(
            (1.5 * 4 * (2.0 * Tracker["constants"]["nnxo"] + 3.0) ** 3), 1.0e9
        )
        #  Estimated data size
        data_size = old_div(
            old_div(
                max(Tracker["nima_per_chunk"])
                * 4
                * float(Tracker["constants"]["nnxo"] ** 2),
                float(Blockdata["no_of_groups"]),
            ),
            1.0e9,
        )
        nnprocs = min(
            Blockdata["no_of_processes_per_group"],
            int((old_div((memory_per_node - data_size * 1.2), volume_size))),
        )
        sp_global_def.sxprint(
            "  MEMORY ESTIMATION.  memory per node = %6.1fGB,  volume size = %6.2fGB, data size per node = %6.2fGB, estimated number of CPUs = %d"
            % (memory_per_node, volume_size, data_size, nnprocs)
        )
        if (memory_per_node - data_size * 1.2 - volume_size) < 0 or (nnprocs == 0):
            nogo = 1
        else:
            nogo = 0
    else:
        nnprocs = 0
        nogo = 0

    nogo = sp_utilities.bcast_number_to_all(
        nogo, source_node=Blockdata["main_node"], mpi_comm=mpi.MPI_COMM_WORLD
    )
    if nogo == 1:
        sp_global_def.ERROR(
            "Insufficient memory to compute final reconstruction",
            "recons3d_final",
            1,
            Blockdata["myid"],
        )
    nnprocs = sp_utilities.bcast_number_to_all(
        nnprocs, source_node=Blockdata["main_node"], mpi_comm=mpi.MPI_COMM_WORLD
    )
    Blockdata["ncpuspernode"] = nnprocs
    Blockdata["nsubset"] = Blockdata["ncpuspernode"] * Blockdata["no_of_groups"]
    create_subgroup()

    oldparamstructure = [[], []]
    newparamstructure = [[], []]
    projdata = [[sp_utilities.model_blank(1, 1)], [sp_utilities.model_blank(1, 1)]]
    original_data = [None, None]
    oldparams = [[], []]
    partids = [None, None]
    partstack = [None, None]

    do3d_final(
        partids,
        partstack,
        original_data,
        oldparams,
        oldparamstructure,
        projdata,
        do_final_iter,
        Blockdata["subgroup_comm"],
    )
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
    line = ""
    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_global_def.sxprint(line, "Final reconstruction is successfully done")
    return


def recons3d_trl_struct_MPI_nosmearing(
    myid, main_node, prjlist, parameters, CTF, upweighted, mpi_comm, target_size
):
    global Tracker, Blockdata
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
    shrink = old_div(float(Tracker["nxinit"]), float(Tracker["constants"]["nnxo"]))
    for im in range(len(prjlist)):
        ct = prjlist[im].get_attr("ctf")
        try:
            bckgn = prjlist[im].get_attr("bckgnoise")
        except:
            bckgn = [1.0] * (old_div(Tracker["constants"]["nnxo"], 2))
        if not upweighted:
            prjlist[im] = sp_filter.filt_table(prjlist[im], bckgn)
        prjlist[im].set_attr_dict({"bckgnoise": bckgn, "ctf": ct})
        phi, theta, psi, s2x, s2y = sp_utilities.get_params_proj(
            prjlist[im], xform="xform.projection"
        )
        junkphi, junktheta, junkpsi, js2x, js2y = (
            parameters[im][0],
            parameters[im][1],
            parameters[im][2],
            parameters[im][3],
            parameters[im][4],
        )
        s2x = js2x - round(js2x)
        s2y = js2y - round(js2y)
        prjlist[im] = sp_fundamentals.fshift(prjlist[im], s2x, s2y)
        r.insert_slice(
            prjlist[im],
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


def rec3d_continuation_nosmearing(original_data, mpi_comm):
    global Tracker, Blockdata

    original_data = [None, None]
    oldparams = [None, None]
    projdata = [None, None]
    partstack = [None, None]
    partids = [None, None]

    temp = Tracker["directory"]
    Tracker["directory"] = os.path.join(
        Tracker["constants"]["masterdir"], "main%03d" % Tracker["mainiteration"]
    )

    for procid in range(2):
        partids[procid] = os.path.join(
            Tracker["constants"]["masterdir"],
            "main%03d" % Tracker["mainiteration"],
            "chunk_%01d_%03d.txt" % (procid, Tracker["mainiteration"]),
        )
        partstack[procid] = os.path.join(
            Tracker["constants"]["masterdir"],
            "main%03d" % Tracker["mainiteration"],
            "params-chunk_%01d_%03d.txt" % (procid, Tracker["mainiteration"]),
        )

        original_data[procid], oldparams[procid], start_end = getindexdata(
            partids[procid],
            partstack[procid],
            os.path.join(
                Tracker["constants"]["masterdir"],
                "main000",
                "particle_groups_%01d.txt" % procid,
            ),
            original_data[procid],
            small_memory=Tracker["constants"]["small_memory"],
            nproc=Blockdata["nproc"],
            myid=Blockdata["myid"],
            mpi_comm=mpi_comm,
        )
        mpi.mpi_barrier(mpi_comm)

        projdata[procid] = get_shrink_data(
            Tracker["nxinit"],
            procid,
            original_data[procid],
            oldparams[procid],
            return_real=False,
            preshift=True,
            apply_mask=False,
            nonorm=True,
            nosmearing=True,
        )
        mpi.mpi_barrier(mpi_comm)

    compute_sigma(
        original_data[0] + original_data[1],
        oldparams[0] + oldparams[1],
        len(oldparams[0]),
        False,
        Blockdata["myid"],
        mpi_comm=mpi_comm,
    )

    # Estimate initial resolution/image size
    if Blockdata["myid"] == Blockdata["nodes"][0]:
        # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
        line = ""
        sp_global_def.sxprint(line, "continuation: reconstruct initial reference")

    for procid in range(2):
        # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
        line = ""
        if Blockdata["myid"] == Blockdata["nodes"][procid]:
            sp_global_def.sxprint(line, "3-D reconstruction of group %d" % procid)

        do3d(
            procid,
            projdata[procid],
            oldparams[procid],
            None,
            None,
            None,
            Blockdata["myid"],
            smearing=False,
            mpi_comm=mpi_comm,
        )

        projdata[procid] = []
        mpi.mpi_barrier(mpi_comm)

    Tracker["maxfrad"] = old_div(Tracker["nxinit"], 2)
    # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
    line = ""
    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_global_def.sxprint(line, "do3d_continuation_get_maps_mpi")
    rec3d_make_maps(compute_fsc=True, regularized=True)
    # if(Blockdata["myid"] == Blockdata["nodes"][0]):  shutil.rmtree(os.path.join(Tracker["directory"], "tempdir"))
    Tracker["directory"] = temp
    return


def update_memory_estimation():
    global Tracker, Blockdata
    if Blockdata["myid"] == Blockdata["main_node"]:
        if Tracker["constants"]["memory_per_node"] == -1.0:
            Tracker["constants"]["memory_per_node"] = (
                Blockdata["no_of_processes_per_group"] * 2.0
            )  # reasonable approximation
        try:
            total_stack = EMAN2_cppwrap.EMUtil.get_image_count(
                Tracker["constants"]["stack"]
            )
        except:
            sp_global_def.sxprint(Tracker["constants"]["stack"], "does not exist")
        image_size = max(
            Tracker["nxinit"], old_div(Tracker["constants"]["nnxo"] * 1.0, 2.0)
        )
        data_size = old_div(
            old_div(
                total_stack * 4 * float(image_size ** 2),
                float(Blockdata["no_of_groups"]),
            ),
            1.0e9,
        )
        volume_size = old_div((1.5 * 4 * (2.0 * image_size + 3.0) ** 3), 1.0e9)
        # nnprocs = Blockdata["no_of_processes_per_group"]
        nnprocs = min(
            Blockdata["no_of_processes_per_group"],
            int(
                old_div(
                    (Tracker["constants"]["memory_per_node"] - data_size * 1.2),
                    volume_size,
                )
            ),
        )
        sp_global_def.sxprint(
            "  MEMORY ESTIMATION.  memory per node = %6.1fGB,  volume size = %6.2fGB, data size per node = %6.2fGB, estimated number of CPUs = %d"
            % (Tracker["constants"]["memory_per_node"], volume_size, data_size, nnprocs)
        )
        memory_per_cpu_3d = (
            old_div(data_size, Blockdata["no_of_processes_per_group"]) * 1.2
            + volume_size
        )
        sp_global_def.sxprint(
            "  Estimated memory consumption per CPU in reconstruction %6.2f"
            % memory_per_cpu_3d
        )
        if (
            Tracker["constants"]["memory_per_node"] - data_size * 1.2 - volume_size
        ) < 0 or (nnprocs == 0):
            nogo = 1
        else:
            nogo = 0
    else:
        nnprocs = 0
        nogo = 0
    nogo = sp_utilities.bcast_number_to_all(
        nogo, source_node=Blockdata["main_node"], mpi_comm=mpi.MPI_COMM_WORLD
    )
    if nogo == 1:
        sp_global_def.ERROR(
            "Insufficient memory to continue refinement from subset",
            "continue_from_subset",
            1,
            Blockdata["myid"],
        )
    nnprocs = sp_utilities.bcast_number_to_all(
        nnprocs, source_node=Blockdata["main_node"], mpi_comm=mpi.MPI_COMM_WORLD
    )
    Blockdata["ncpuspernode"] = nnprocs
    Blockdata["nsubset"] = Blockdata["ncpuspernode"] * Blockdata["no_of_groups"]
    create_subgroup()


def update_tracker(shell_line_command):
    global Tracker, Blockdata
    # reset parameters for a restart run; update only those specified options in restart
    # 1. maxit is not included.
    # 2. those sigmas for local search can be considered included
    parser_no_default = optparse.OptionParser()
    parser_no_default.add_option("--radius", type="int")
    parser_no_default.add_option("--xr", type="float")
    parser_no_default.add_option("--ts", type="float")
    parser_no_default.add_option("--inires", type="float")
    parser_no_default.add_option("--mask3D", type="string")
    parser_no_default.add_option("--function", type="string")
    parser_no_default.add_option(
        "--symmetry", type="string"
    )  # rare to change sym; however, keep it an option.
    parser_no_default.add_option("--delta", type="float")
    parser_no_default.add_option("--shake", type="float")
    parser_no_default.add_option("--small_memory", action="store_true")
    parser_no_default.add_option("--ccfpercentage", type="float")
    parser_no_default.add_option("--nonorm", action="store_true")
    parser_no_default.add_option("--memory_per_node", type="float")
    parser_no_default.add_option("--an", type="float")

    parser_no_default.add_option("--do_final", type="int")
    parser_no_default.add_option("--local_refinement", action="store_true")
    parser_no_default.add_option("--continuation_orgstack", type="string")
    parser_no_default.add_option("--continuation_initvol", type="string")
    parser_no_default.add_option("--subset", type="string")
    parser_no_default.add_option("--oldrefdir", type="string")
    parser_no_default.add_option("--continuation_iter", type="int")
    parser_no_default.add_option("--continuation_smearing", type="int")
    parser_no_default.add_option("--keep_groups", action="store_true")

    (options_no_default_value, args) = parser_no_default.parse_args(shell_line_command)

    #  This is for printout only
    tempdict = {}

    if options_no_default_value.radius != None:
        Tracker["constants"]["radius"] = options_no_default_value.radius
        tempdict["radius"] = Tracker["constants"]["radius"]

    if options_no_default_value.xr != None:
        Tracker["xr"] = options_no_default_value.xr
        tempdict["xr"] = Tracker["xr"]

    if options_no_default_value.ts != None:
        Tracker["ts"] = options_no_default_value.ts
        tempdict["ts"] = Tracker["ts"]

    if options_no_default_value.inires != None:
        Tracker["constants"]["inires"] = options_no_default_value.inires
        Tracker["constants"]["inires"] = int(
            old_div(
                Tracker["constants"]["nnxo"] * Tracker["constants"]["pixel_size"],
                Tracker["constants"]["inires"],
            )
            + 0.5
        )
        Tracker["currentres"] = Tracker["constants"]["inires"]
        tempdict["currentres"] = Tracker["currentres"]

    if options_no_default_value.delta != None:
        Tracker["delta"] = options_no_default_value.delta
        tempdict["delta"] = Tracker["delta"]
    if options_no_default_value.shake != None:
        Tracker["constants"]["shake"] = options_no_default_value.shake
        tempdict["shake"] = Tracker["constants"]["shake"]
    if (
        options_no_default_value.symmetry != None
    ):  # this rarely happens. However, keep it an option.
        sym = options_no_default_value.symmetry
        Tracker["constants"]["symmetry"] = sym[0].lower() + sym[1:]
        tempdict["symmetry"] = Tracker["constants"]["symmetry"]
    if options_no_default_value.mask3D != None:
        Tracker["constants"]["mask3D"] = options_no_default_value.mask3D
        tempdict["mask3D"] = Tracker["constants"]["mask3D"]
    if options_no_default_value.ccfpercentage != None:
        Tracker["constants"]["ccfpercentage"] = old_div(
            options_no_default_value.ccfpercentage, 100.0
        )
        tempdict["ccfpercentage"] = Tracker["constants"]["ccfpercentage"]
    if options_no_default_value.nonorm != None:
        Tracker["constants"]["nonorm"] = options_no_default_value.nonorm
        tempdict["nonorm"] = Tracker["constants"]["nonorm"]
    if options_no_default_value.small_memory != None:
        Tracker["constants"]["small_memory"] = options_no_default_value.small_memory
        tempdict["small_memory"] = Tracker["constants"]["small_memory"]
    if options_no_default_value.memory_per_node != None:
        Tracker["constants"][
            "memory_per_node"
        ] = options_no_default_value.memory_per_node
        tempdict["memory_per_node"] = Tracker["constants"]["memory_per_node"]
    ### continuation
    if options_no_default_value.continuation_orgstack != None:
        Tracker["constants"][
            "continuation_orgstack"
        ] = options_no_default_value.continuation_orgstack
    if options_no_default_value.keep_groups != None:
        Tracker["constants"]["keep_groups"] = options_no_default_value.keep_groups
    if options_no_default_value.subset != None:
        Tracker["constants"]["subset"] = options_no_default_value.subset
    if options_no_default_value.oldrefdir != None:
        Tracker["constants"]["oldrefdir"] = options_no_default_value.oldrefdir
    if options_no_default_value.continuation_initvol != None:
        Tracker["constants"][
            "continuation_initvol"
        ] = options_no_default_value.continuation_initvol
    if options_no_default_value.continuation_iter != None:
        Tracker["constants"][
            "continuation_iter"
        ] = options_no_default_value.continuation_iter
    if options_no_default_value.function != None:
        Tracker["constants"]["function"] = options_no_default_value.function
    if options_no_default_value.an != None:
        Tracker["constants"]["an"] = options_no_default_value.an
        tempdict["an"] = Tracker["constants"]["an"]

    if (Blockdata["myid"] == Blockdata["main_node"]) and (len(tempdict) > 0):
        print_dict(tempdict, "Updated settings")

    Blockdata["symclass"] = sp_fundamentals.symclass(Tracker["constants"]["symmetry"])
    return


def rec3d_make_maps(compute_fsc=True, regularized=True):
    global Tracker, Blockdata

    # final reconstruction: compute_fsc = False; regularized = False
    # tempdir is removed in the end of the function
    if compute_fsc:
        if Blockdata["no_of_groups"] == 1:
            if Blockdata["myid"] == Blockdata["nodes"][0]:
                tvol0 = sp_utilities.get_im(
                    os.path.join(
                        Tracker["directory"],
                        "tempdir",
                        "tvol_0_%03d.hdf" % (Tracker["mainiteration"]),
                    )
                )
                tweight0 = sp_utilities.get_im(
                    os.path.join(
                        Tracker["directory"],
                        "tempdir",
                        "tweight_0_%03d.hdf" % (Tracker["mainiteration"]),
                    )
                )
                tvol1 = sp_utilities.get_im(
                    os.path.join(
                        Tracker["directory"],
                        "tempdir",
                        "tvol_1_%03d.hdf" % (Tracker["mainiteration"]),
                    )
                )
                tweight1 = sp_utilities.get_im(
                    os.path.join(
                        Tracker["directory"],
                        "tempdir",
                        "tweight_1_%03d.hdf" % (Tracker["mainiteration"]),
                    )
                )
                EMAN2_cppwrap.Util.fuse_low_freq(
                    tvol0,
                    tvol1,
                    tweight0,
                    tweight1,
                    2 * Tracker["constants"]["fuse_freq"],
                )
                shrank0 = stepone(tvol0, tweight0)
                shrank1 = stepone(tvol1, tweight1)
                #  Note shrank volumes are Fourier uncentered.
                cfsc = sp_statistics.fsc(shrank0, shrank1)[1]
                del shrank0, shrank1
                if Tracker["nxinit"] < Tracker["constants"]["nnxo"]:
                    cfsc = cfsc[: old_div(Tracker["nxinit"], 2) + 1]
                    for i in range(
                        len(cfsc), old_div(Tracker["constants"]["nnxo"], 2) + 1
                    ):
                        cfsc.append(0.0)
                lcfsc = len(cfsc)
                # --  memory_check(Blockdata["myid"],"second node, after stepone")
            else:
                #  receive fsc
                lcfsc = 0
        else:
            if (
                Blockdata["myid"] == Blockdata["nodes"][1]
            ):  # It has to be 1 to avoid problem with tvol1 not closed on the disk
                # --  memory_check(Blockdata["myid"],"first node, before stepone")
                #  read volumes, shrink
                tvol0 = sp_utilities.get_im(
                    os.path.join(
                        Tracker["directory"],
                        "tempdir",
                        "tvol_0_%03d.hdf" % (Tracker["mainiteration"]),
                    )
                )
                tweight0 = sp_utilities.get_im(
                    os.path.join(
                        Tracker["directory"],
                        "tempdir",
                        "tweight_0_%03d.hdf" % (Tracker["mainiteration"]),
                    )
                )
                tvol1 = sp_utilities.get_im(
                    os.path.join(
                        Tracker["directory"],
                        "tempdir",
                        "tvol_1_%03d.hdf" % (Tracker["mainiteration"]),
                    )
                )
                tweight1 = sp_utilities.get_im(
                    os.path.join(
                        Tracker["directory"],
                        "tempdir",
                        "tweight_1_%03d.hdf" % (Tracker["mainiteration"]),
                    )
                )
                EMAN2_cppwrap.Util.fuse_low_freq(
                    tvol0,
                    tvol1,
                    tweight0,
                    tweight1,
                    2 * Tracker["constants"]["fuse_freq"],
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
                # --  memory_check(Blockdata["myid"],"first node, after stepone")
            elif Blockdata["myid"] == Blockdata["nodes"][0]:
                # --  memory_check(Blockdata["myid"],"second node, before stepone")
                #  read volumes, shrink
                tag = 7007
                tvol1 = sp_utilities.recv_EMData(
                    Blockdata["nodes"][1], tag, mpi.MPI_COMM_WORLD
                )
                tweight1 = sp_utilities.recv_EMData(
                    Blockdata["nodes"][1], tag, mpi.MPI_COMM_WORLD
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
                #  Get shrank volume, do fsc, send it to all
                shrank0 = sp_utilities.recv_EMData(
                    Blockdata["nodes"][1], tag, mpi.MPI_COMM_WORLD
                )
                #  Note shrank volumes are Fourier uncentered.
                cfsc = sp_statistics.fsc(shrank0, shrank1)[1]
                del shrank0, shrank1
                if Tracker["nxinit"] < Tracker["constants"]["nnxo"]:
                    cfsc = cfsc[: old_div(Tracker["nxinit"], 2) + 1]
                    for i in range(
                        len(cfsc), old_div(Tracker["constants"]["nnxo"], 2) + 1
                    ):
                        cfsc.append(0.0)
                lcfsc = len(cfsc)
                # --  memory_check(Blockdata["myid"],"second node, after stepone")
            else:
                #  receive fsc
                lcfsc = 0
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        lcfsc = sp_utilities.bcast_number_to_all(lcfsc)
        if Blockdata["myid"] != Blockdata["nodes"][0]:
            cfsc = [0.0] * lcfsc
        cfsc = sp_utilities.bcast_list_to_all(
            cfsc, Blockdata["myid"], Blockdata["nodes"][0]
        )
        if Blockdata["myid"] == Blockdata["main_node"]:
            sp_utilities.write_text_file(
                cfsc,
                os.path.join(
                    Tracker["directory"], "driver_%03d.txt" % (Tracker["mainiteration"])
                ),
            )
            out_fsc(cfsc)

    # Now that we have the curve, do the reconstruction
    Tracker["maxfrad"] = old_div(Tracker["nxinit"], 2)
    if Blockdata["no_of_groups"] > 1:
        lorder = [0, 0]  # Two blocks in parallel
    elif Blockdata["no_of_groups"] == 1:
        lorder = [0, 1]  # One after another

    if regularized:
        for iorder in range(2):
            if iorder == lorder[0]:
                if Blockdata["color"] == Blockdata["node_volume"][1]:
                    # --  memory_check(Blockdata["myid"],"first node, before steptwo")
                    #  compute filtered volume
                    if Blockdata["myid_on_node"] == 0:
                        treg0 = sp_utilities.get_im(
                            os.path.join(
                                Tracker["directory"],
                                "tempdir",
                                "trol_0_%03d.hdf" % (Tracker["mainiteration"]),
                            )
                        )
                    else:
                        tvol0 = sp_utilities.model_blank(1)
                        tweight0 = sp_utilities.model_blank(1)
                        treg0 = sp_utilities.model_blank(1)
                    tvol0 = steptwo_mpi(
                        tvol0,
                        tweight0,
                        treg0,
                        cfsc,
                        True,
                        color=Blockdata["node_volume"][1],
                    )
                    del tweight0, treg0
                    if Blockdata["myid_on_node"] == 0:
                        # --  memory_check(Blockdata["myid"],"first node, before masking")
                        if Tracker["mainiteration"] == 1:
                            # At a first iteration truncate resolution at the initial resolution set by the user
                            for i in range(len(cfsc)):
                                if i < Tracker["constants"]["inires"] + 1:
                                    cfsc[i] = 1.0
                                if i == Tracker["constants"]["inires"] + 1:
                                    cfsc[i] = 0.5
                                elif i > Tracker["constants"]["inires"] + 1:
                                    cfsc[i] = 0.0
                            tvol0 = sp_filter.filt_table(tvol0, cfsc)
                            if Blockdata["no_of_groups"] > 1:
                                del cfsc

                        user_func = sp_user_functions.factory[
                            Tracker["constants"]["user_func"]
                        ]
                        # ref_data = [tvol0, Tracker, mainiteration]
                        ref_data = [tvol0, Tracker, Tracker["mainiteration"]]
                        # --  #--  memory_check(Blockdata["myid"],"first node, after masking")
                        user_func(ref_data).write_image(
                            os.path.join(
                                Tracker["directory"],
                                "vol_0_%03d.hdf" % (Tracker["mainiteration"]),
                            )
                        )
                    del tvol0
                    # --  memory_check(Blockdata["myid"],"first node, after 2 steptwo")
            if iorder == lorder[1]:
                if Blockdata["color"] == Blockdata["node_volume"][0]:
                    # --  memory_check(Blockdata["myid"],"second node, before steptwo")
                    #  compute filtered volume
                    if Blockdata["myid_on_node"] == 0:
                        treg1 = sp_utilities.get_im(
                            os.path.join(
                                Tracker["directory"],
                                "tempdir",
                                "trol_1_%03d.hdf" % (Tracker["mainiteration"]),
                            )
                        )
                    else:
                        tvol1 = sp_utilities.model_blank(1)
                        tweight1 = sp_utilities.model_blank(1)
                        treg1 = sp_utilities.model_blank(1)

                    tvol1 = steptwo_mpi(
                        tvol1,
                        tweight1,
                        treg1,
                        cfsc,
                        True,
                        color=Blockdata["node_volume"][0],
                    )
                    del tweight1, treg1
                    if Blockdata["myid_on_node"] == 0:
                        # --  memory_check(Blockdata["myid"],"second node, before masking")
                        if Tracker["mainiteration"] == 1:
                            # At a first iteration truncate resolution at the initial resolution set by the user
                            for i in range(len(cfsc)):
                                if i < Tracker["constants"]["inires"] + 1:
                                    cfsc[i] = 1.0
                                if i == Tracker["constants"]["inires"] + 1:
                                    cfsc[i] = 0.5
                                elif i > Tracker["constants"]["inires"] + 1:
                                    cfsc[i] = 0.0
                            tvol1 = sp_filter.filt_table(tvol1, cfsc)
                            del cfsc
                        user_func = sp_user_functions.factory[
                            Tracker["constants"]["user_func"]
                        ]
                        # ref_data = [tvol1, Tracker, mainiteration]
                        ref_data = [tvol1, Tracker, Tracker["mainiteration"]]
                        # --  #--  memory_check(Blockdata["myid"],"first node, after masking")
                        user_func(ref_data).write_image(
                            os.path.join(
                                Tracker["directory"],
                                "vol_1_%03d.hdf" % (Tracker["mainiteration"]),
                            )
                        )
                    del tvol1
                    # --  memory_check(Blockdata["myid"],"second node, after 2 steptwo")
            #  Here end per node execution.
        if Blockdata["myid"] == Blockdata["nodes"][0]:
            sp_helix_sphire.shutil.rmtree(os.path.join(Tracker["directory"], "tempdir"))
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    else:
        if Blockdata["no_of_groups"] == 1:
            for iproc in range(2):
                if Blockdata["myid_on_node"] == 0:
                    tvol0 = sp_utilities.get_im(
                        os.path.join(
                            Tracker["directory"],
                            os.path.join(
                                "tempdir",
                                "tvol_0_%03d.hdf" % (Tracker["mainiteration"]),
                            ),
                        )
                    )
                    tweight0 = sp_utilities.get_im(
                        os.path.join(
                            Tracker["directory"],
                            os.path.join(
                                "tempdir",
                                "tweight_0_%03d.hdf" % (Tracker["mainiteration"]),
                            ),
                        )
                    )
                    tvol1 = sp_utilities.get_im(
                        os.path.join(
                            Tracker["directory"],
                            os.path.join(
                                "tempdir",
                                "tvol_1_%03d.hdf" % (Tracker["mainiteration"]),
                            ),
                        )
                    )
                    tweight1 = sp_utilities.get_im(
                        os.path.join(
                            Tracker["directory"],
                            os.path.join(
                                "tempdir",
                                "tweight_1_%03d.hdf" % (Tracker["mainiteration"]),
                            ),
                        )
                    )
                    EMAN2_cppwrap.Util.fuse_low_freq(
                        tvol0,
                        tvol1,
                        tweight0,
                        tweight1,
                        2 * Tracker["constants"]["fuse_freq"],
                    )
                    treg = sp_utilities.get_im(
                        os.path.join(
                            Tracker["directory"],
                            "tempdir",
                            "trol_%d_%03d.hdf" % ((iproc, Tracker["mainiteration"])),
                        )
                    )
                else:
                    treg = sp_utilities.model_blank(1)
                    if iproc == 0:
                        tvol0 = sp_utilities.model_blank(1)
                        tweight0 = sp_utilities.model_blank(1)
                    else:
                        tvol1 = sp_utilities.model_blank(1)
                        tweight1 = sp_utilities.model_blank(1)
                if iproc == 0:
                    tvol0 = steptwo_mpi(
                        tvol0,
                        tweight0,
                        treg,
                        None,
                        False,
                        color=Blockdata["node_volume"][0],
                    )
                else:
                    tvol1 = steptwo_mpi(
                        tvol1,
                        tweight1,
                        treg,
                        None,
                        False,
                        color=Blockdata["node_volume"][0],
                    )
                if Blockdata["myid_on_node"] == 0:
                    if iproc == 0:
                        tvol0.write_image(
                            os.path.join(
                                Tracker["constants"]["masterdir"],
                                "vol_%d_unfil_%03d.hdf" % (iproc, final_iter),
                            )
                        )
                    else:
                        tvol1.write_image(
                            os.path.join(
                                Tracker["constants"]["masterdir"],
                                "vol_%d_unfil_%03d.hdf" % (iproc, final_iter),
                            )
                        )
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        else:
            if Blockdata["myid"] == Blockdata["main_shared_nodes"][1]:
                # post-insertion operations, done only in main_node
                tvol0 = sp_utilities.get_im(
                    os.path.join(
                        Tracker["directory"],
                        os.path.join(
                            "tempdir", "tvol_0_%03d.hdf" % Tracker["mainiteration"]
                        ),
                    )
                )
                tweight0 = sp_utilities.get_im(
                    os.path.join(
                        Tracker["directory"],
                        os.path.join(
                            "tempdir", "tweight_0_%03d.hdf" % Tracker["mainiteration"]
                        ),
                    )
                )
                tvol1 = sp_utilities.get_im(
                    os.path.join(
                        Tracker["directory"],
                        os.path.join(
                            "tempdir", "tvol_1_%03d.hdf" % Tracker["mainiteration"]
                        ),
                    )
                )
                tweight1 = sp_utilities.get_im(
                    os.path.join(
                        Tracker["directory"],
                        os.path.join(
                            "tempdir", "tweight_1_%03d.hdf" % Tracker["mainiteration"]
                        ),
                    )
                )
                EMAN2_cppwrap.Util.fuse_low_freq(
                    tvol0,
                    tvol1,
                    tweight0,
                    tweight1,
                    2 * Tracker["constants"]["fuse_freq"],
                )
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            if Blockdata["myid"] == Blockdata["main_shared_nodes"][1]:
                tag = 7007
                sp_utilities.send_EMData(
                    tvol1, Blockdata["main_shared_nodes"][0], tag, mpi.MPI_COMM_WORLD
                )
                sp_utilities.send_EMData(
                    tweight1, Blockdata["main_shared_nodes"][0], tag, mpi.MPI_COMM_WORLD
                )
                tvol0.set_attr_dict(
                    {
                        "is_complex": 1,
                        "is_fftodd": 1,
                        "is_complex_ri": 1,
                        "is_fftpad": 1,
                    }
                )

            elif Blockdata["myid"] == Blockdata["main_shared_nodes"][0]:
                tag = 7007
                tvol1 = sp_utilities.recv_EMData(
                    Blockdata["main_shared_nodes"][1], tag, mpi.MPI_COMM_WORLD
                )
                tweight1 = sp_utilities.recv_EMData(
                    Blockdata["main_shared_nodes"][1], tag, mpi.MPI_COMM_WORLD
                )
                tvol1.set_attr_dict(
                    {
                        "is_complex": 1,
                        "is_fftodd": 1,
                        "is_complex_ri": 1,
                        "is_fftpad": 1,
                    }
                )
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            if Blockdata["color"] == Blockdata["node_volume"][1]:
                if Blockdata["myid"] == Blockdata["main_shared_nodes"][1]:
                    treg0 = sp_utilities.get_im(
                        os.path.join(
                            Tracker["directory"],
                            "tempdir",
                            "trol_0_%03d.hdf" % (Tracker["mainiteration"]),
                        )
                    )
                else:
                    tvol0 = sp_utilities.model_blank(1)
                    tweight0 = sp_utilities.model_blank(1)
                    treg0 = sp_utilities.model_blank(1)
                tvol0 = steptwo_mpi(
                    tvol0,
                    tweight0,
                    treg0,
                    None,
                    False,
                    color=Blockdata["node_volume"][1],
                )
                del tweight0, treg0
                if Blockdata["myid_on_node"] == 0:
                    tvol0.write_image(
                        os.path.join(
                            Tracker["constants"]["masterdir"],
                            "vol_0_unfil_%03d.hdf" % Tracker["mainiteration"],
                        )
                    )

            elif Blockdata["color"] == Blockdata["node_volume"][0]:
                if Blockdata["myid"] == Blockdata["main_shared_nodes"][0]:
                    treg1 = sp_utilities.get_im(
                        os.path.join(
                            Tracker["directory"],
                            "tempdir",
                            "trol_1_%03d.hdf" % (Tracker["mainiteration"]),
                        )
                    )
                else:
                    tvol1 = sp_utilities.model_blank(1)
                    tweight1 = sp_utilities.model_blank(1)
                    treg1 = sp_utilities.model_blank(1)
                tvol1 = steptwo_mpi(
                    tvol1,
                    tweight1,
                    treg1,
                    None,
                    False,
                    color=Blockdata["node_volume"][0],
                )
                del tweight1, treg1
                if Blockdata["myid_on_node"] == 0:
                    tvol1.write_image(
                        os.path.join(
                            Tracker["constants"]["masterdir"],
                            "vol_1_unfil_%03d.hdf" % Tracker["mainiteration"],
                        )
                    )
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    return


def refinement_one_iteration(
    partids,
    partstack,
    original_data,
    oldparams,
    projdata,
    general_mode=True,
    continuation_mode=False,
):
    global Tracker, Blockdata
    #  READ DATA AND COMPUTE SIGMA2   ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
    start_end = [[], []]
    for procid in range(2):
        original_data[procid], oldparams[procid], start_end[procid] = getindexdata(
            partids[procid],
            partstack[procid],
            os.path.join(
                Tracker["constants"]["masterdir"],
                "main000",
                "particle_groups_%01d.txt" % procid,
            ),
            original_data[procid],
            small_memory=Tracker["constants"]["small_memory"],
            nproc=Blockdata["nproc"],
            myid=Blockdata["myid"],
            mpi_comm=mpi.MPI_COMM_WORLD,
        )

    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    if general_mode:
        if Tracker["mainiteration"] == 1:
            dryrun = False
        else:
            dryrun = True
    elif continuation_mode:
        dryrun = True
    else:
        pass

    compute_sigma(
        original_data[0] + original_data[1],
        oldparams[0] + oldparams[1],
        len(oldparams[0]),
        dryrun,
        Blockdata["myid"],
    )

    #  REFINEMENT   ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    refang, rshifts, coarse_angles, coarse_shifts = get_refangs_and_shifts()
    if Tracker["constants"]["shake"] > 0.0:
        if Blockdata["myid"] == Blockdata["main_node"]:
            shakenumber = random.uniform(
                -Tracker["constants"]["shake"], Tracker["constants"]["shake"]
            )
        else:
            shakenumber = 0.0
        shakenumber = sp_utilities.bcast_number_to_all(
            shakenumber, source_node=Blockdata["main_node"]
        )
        # it has to be rounded as the number written to the disk is rounded,
        #  so if there is discrepancy one cannot reproduce iteration.
        shakenumber = round(shakenumber, 5)

        rangle = shakenumber * Tracker["delta"]
        rshift = shakenumber * Tracker["ts"]
        refang = Blockdata["symclass"].reduce_anglesets(
            sp_fundamentals.rotate_params(refang, [-rangle, -rangle, -rangle])
        )
        coarse_angles = Blockdata["symclass"].reduce_anglesets(
            sp_fundamentals.rotate_params(coarse_angles, [-rangle, -rangle, -rangle])
        )
        shakegrid(rshifts, rshift)
        shakegrid(coarse_shifts, rshift)

        if Blockdata["myid"] == Blockdata["main_node"]:
            sp_utilities.write_text_row(
                [[shakenumber, rangle, rshift]],
                os.path.join(Tracker["directory"], "randomize_search.txt"),
            )
    else:
        rangle = 0.0
        rshift = 0.0

    if Blockdata["myid"] == Blockdata["main_node"]:
        sp_utilities.write_text_row(
            refang, os.path.join(Tracker["directory"], "refang.txt")
        )
        sp_utilities.write_text_row(
            rshifts, os.path.join(Tracker["directory"], "rshifts.txt")
        )
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    newparamstructure = [[], []]
    raw_vol = [[], []]
    norm_per_particle = [[], []]

    for procid in range(2):
        Tracker["refvol"] = os.path.join(
            Tracker["previousoutputdir"],
            "vol_%01d_%03d.hdf" % (procid, Tracker["mainiteration"] - 1),
        )

        Tracker["nxpolar"] = Tracker[
            "nxinit"
        ]  # min( 3*Tracker["nxinit"], Tracker["constants"]["nnxo"] )
        # Tracker["nxpolar"] = min( 2*Tracker["nxinit"], Tracker["constants"]["nnxo"] )
        if general_mode:
            if Tracker["state"] == "INITIAL":
                newparamstructure[procid], norm_per_particle[procid] = ali3D_polar_ccc(
                    refang,
                    rshifts,
                    coarse_angles,
                    coarse_shifts,
                    procid,
                    original_data[procid],
                    oldparams[procid],
                    preshift=True,
                    apply_mask=True,
                    nonorm=Tracker["constants"]["nonorm"],
                    applyctf=False,
                )
            elif Tracker["state"] == "PRIMARY":
                newparamstructure[procid], norm_per_particle[
                    procid
                ] = ali3D_primary_polar(
                    refang,
                    rshifts,
                    coarse_angles,
                    coarse_shifts,
                    procid,
                    original_data[procid],
                    oldparams[procid],
                    preshift=True,
                    apply_mask=True,
                    nonorm=Tracker["constants"]["nonorm"],
                    applyctf=True,
                )
            elif Tracker["state"] == "EXHAUSTIVE":
                newparamstructure[procid], norm_per_particle[procid] = ali3D_polar(
                    refang,
                    rshifts,
                    coarse_angles,
                    coarse_shifts,
                    procid,
                    original_data[procid],
                    oldparams[procid],
                    preshift=True,
                    apply_mask=True,
                    nonorm=Tracker["constants"]["nonorm"],
                    applyctf=True,
                )
            elif (Tracker["state"] == "RESTRICTED") or (Tracker["state"] == "FINAL"):
                newparamstructure[procid], norm_per_particle[
                    procid
                ] = ali3D_local_polar(
                    refang,
                    rshifts,
                    coarse_angles,
                    coarse_shifts,
                    procid,
                    original_data[procid],
                    oldparams[procid],
                    preshift=True,
                    apply_mask=True,
                    nonorm=Tracker["constants"]["nonorm"],
                    applyctf=True,
                )
            else:
                sp_global_def.ERROR(
                    "sxmeridien",
                    "Incorrect state  %s" % Tracker["state"],
                    1,
                    Blockdata["myid"],
                )

        elif continuation_mode:
            if Tracker["state"] == "PRIMARY":
                ###print("   ",Blockdata["myid"],len(refang),len(rshifts),len(coarse_angles),len(coarse_shifts),len(original_data[procid]), len(oldparams[procid]))
                ###print("   ",Blockdata["myid"],original_data[0][0], oldparams[0][0])
                newparamstructure[procid], norm_per_particle[
                    procid
                ] = ali3D_primary_local_polar(
                    refang,
                    rshifts,
                    coarse_angles,
                    coarse_shifts,
                    procid,
                    original_data[procid],
                    oldparams[procid],
                    preshift=True,
                    apply_mask=True,
                    nonorm=Tracker["constants"]["nonorm"],
                    applyctf=True,
                )
            elif (Tracker["state"] == "RESTRICTED") or (Tracker["state"] == "FINAL"):
                ###print("   ",Blockdata["myid"],len(refang),len(rshifts),len(coarse_angles),len(coarse_shifts),len(original_data[procid]), len(oldparams[procid]))
                ###print("   ",Blockdata["myid"],original_data[0][0], oldparams[0][0])
                newparamstructure[procid], norm_per_particle[
                    procid
                ] = ali3D_local_polar(
                    refang,
                    rshifts,
                    coarse_angles,
                    coarse_shifts,
                    procid,
                    original_data[procid],
                    oldparams[procid],
                    preshift=True,
                    apply_mask=True,
                    nonorm=Tracker["constants"]["nonorm"],
                    applyctf=True,
                )
            else:
                sp_global_def.ERROR(
                    "sxmeridien",
                    "Incorrect state  %s" % Tracker["state"],
                    1,
                    Blockdata["myid"],
                )
        else:
            pass

        qt = old_div(
            old_div(1.0, Tracker["constants"]["nnxo"]), Tracker["constants"]["nnxo"]
        )
        params = []
        for im in range(len(newparamstructure[procid])):
            #  Select only one best
            hash = newparamstructure[procid][im][2][0][0]
            ishift = hash % 1000
            ipsi = (old_div(hash, 1000)) % 100000
            iang = old_div(hash, 100000000)
            params.append(
                [
                    refang[iang][0],
                    refang[iang][1],
                    (refang[iang][2] + ipsi * Tracker["delta"]) % 360.0,
                    rshifts[ishift][0] + oldparams[procid][im][3],
                    rshifts[ishift][1] + oldparams[procid][im][4],
                    newparamstructure[procid][im][-1][0][1],
                    norm_per_particle[procid][im] * qt,
                    norm_per_particle[procid][im],
                ]
            )

        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

        params = sp_utilities.wrap_mpi_gatherv(
            params, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        #  store params
        if Blockdata["myid"] == Blockdata["main_node"]:
            # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
            line = ""
            sp_global_def.sxprint(
                line,
                "Executed successfully: ",
                "Projection matching, state: %s, number of images:%7d"
                % (Tracker["state"], len(params)),
            )
            sp_utilities.write_text_row(
                params,
                os.path.join(
                    Tracker["directory"],
                    "params-chunk_%01d_%03d.txt" % (procid, Tracker["mainiteration"]),
                ),
            )
        del params

        projdata[procid] = []
        if Tracker["constants"]["small_memory"]:
            original_data[procid], oldparams[procid], start_end[procid] = getindexdata(
                partids[procid],
                partstack[procid],
                os.path.join(
                    Tracker["constants"]["masterdir"],
                    "main000",
                    "particle_groups_%01d.txt" % procid,
                ),
                original_data[procid],
                small_memory=Tracker["constants"]["small_memory"],
                nproc=Blockdata["nproc"],
                myid=Blockdata["myid"],
                mpi_comm=mpi.MPI_COMM_WORLD,
            )

        if Tracker["changed_delta"]:
            org_nxinit = Tracker["nxinit"]
            Tracker["nxinit"] = Tracker["constants"]["nnxo"]

        projdata[procid] = get_shrink_data(
            Tracker["nxinit"],
            procid,
            original_data[procid],
            oldparams[procid],
            return_real=False,
            preshift=True,
            apply_mask=False,
            nonorm=True,
        )

        oldparams[procid] = []
        if Tracker["constants"]["small_memory"]:
            original_data[procid] = []

        if Blockdata["myid"] == Blockdata["main_node"]:
            # Carry over chunk information
            cmd = "{} {} {}".format(
                "cp -p",
                os.path.join(
                    Tracker["previousoutputdir"],
                    "chunk_%01d_%03d.txt" % (procid, Tracker["mainiteration"] - 1),
                ),
                os.path.join(
                    Tracker["directory"],
                    "chunk_%01d_%03d.txt" % (procid, Tracker["mainiteration"]),
                ),
            )
            junk = sp_utilities.cmdexecute(cmd)

        chunk_file = os.path.join(
            Tracker["directory"],
            "chunk_%01d_%03d.txt" % (procid, Tracker["mainiteration"]),
        )
        params_file = os.path.join(
            Tracker["directory"],
            "params-chunk_%01d_%03d.txt" % (procid, Tracker["mainiteration"]),
        )
        outlier_file = os.path.join(
            Tracker["directory"],
            "outlier-params-chunk_%01d_%03d.txt" % (procid, Tracker["mainiteration"]),
        )

        outlier_list, outlier_full = calculate_prior_values(
            tracker=Tracker,
            blockdata=Blockdata,
            outlier_file=outlier_file,
            params_file=params_file,
            chunk_file=chunk_file,
            im_start=start_end[procid][0],
            im_end=start_end[procid][1],
            procid=procid,
        )

        if Tracker["prior"]["apply_prior"]:
            pass
        else:
            outlier_list = [0] * len(projdata[procid])

        assert len(newparamstructure[procid]) == len(projdata[procid])
        assert len(newparamstructure[procid]) == len(norm_per_particle[procid])
        assert len(newparamstructure[procid]) == len(outlier_list), (
            len(newparamstructure[procid]),
            len(outlier_list),
        )
        assert len(norm_per_particle[procid]) == len(outlier_list), (
            len(norm_per_particle[procid]),
            len(outlier_list),
        )
        assert len(projdata[procid]) == len(outlier_list), (
            len(projdata[procid]),
            len(outlier_list),
        )

        norm_per_particle_outlier = []
        newparamstructure_outlier = []
        projdata_outlier = []
        for idx, entry in enumerate(outlier_list):
            if int(entry) == 0:
                projdata_outlier.append(projdata[procid][idx])
                norm_per_particle_outlier.append(norm_per_particle[procid][idx])
                newparamstructure_outlier.append(newparamstructure[procid][idx])
        projdata[procid] = []

        assert len(projdata_outlier) == len(
            [entry for entry in outlier_list if int(entry) == 0]
        ), (
            len(projdata_outlier),
            len([entry for entry in outlier_list if int(entry) == 0]),
        )
        assert len(norm_per_particle_outlier) == len(newparamstructure_outlier)
        assert len(norm_per_particle_outlier) == len(projdata_outlier)

        total_left_particles = sp_utilities.wrap_mpi_gatherv(
            norm_per_particle_outlier, Blockdata["main_node"], mpi.MPI_COMM_WORLD
        )
        if Blockdata["myid"] == Blockdata["main_node"]:
            sp_global_def.sxprint(
                "Use {0} particles for final reconstruction!".format(
                    len(total_left_particles)
                )
            )

        do3d(
            procid,
            projdata_outlier,
            newparamstructure_outlier,
            refang,
            rshifts,
            norm_per_particle_outlier,
            Blockdata["myid"],
            smearing=True,
            mpi_comm=mpi.MPI_COMM_WORLD,
        )
        projdata_outlier = []
        newparamstructure_outlier = []
        norm_per_particle_outlier = []
        if Tracker["changed_delta"]:
            Tracker["nxinit"] = org_nxinit

        if Blockdata["myid_on_node"] == 0:
            for kproc in range(Blockdata["no_of_processes_per_group"]):
                if kproc == 0:
                    dump_object_to_json(
                        os.path.join(
                            Tracker["constants"]["masterdir"],
                            "main%03d" % Tracker["mainiteration"],
                            "oldparamstructure",
                            "oldparamstructure_%01d_%03d_%03d.json"
                            % (procid, Blockdata["myid"], Tracker["mainiteration"]),
                        ),
                        newparamstructure[procid],
                    )
                else:
                    dummy = sp_utilities.wrap_mpi_recv(kproc, Blockdata["shared_comm"])
                    dump_object_to_json(
                        os.path.join(
                            Tracker["constants"]["masterdir"],
                            "main%03d" % Tracker["mainiteration"],
                            "oldparamstructure",
                            "oldparamstructure_%01d_%03d_%03d.json"
                            % (
                                procid,
                                (
                                    Blockdata["color"]
                                    * Blockdata["no_of_processes_per_group"]
                                    + kproc
                                ),
                                Tracker["mainiteration"],
                            ),
                        ),
                        dummy,
                    )
                    del dummy
                try:
                    os.remove(
                        os.path.join(
                            Tracker["constants"]["masterdir"],
                            "main%03d" % Tracker["mainiteration"],
                            "oldparamstructure",
                            "oldparamstructure_%01d_%03d_%03d.json"
                            % (
                                procid,
                                (
                                    Blockdata["color"]
                                    * Blockdata["no_of_processes_per_group"]
                                    + kproc
                                ),
                                Tracker["mainiteration"] - 1,
                            ),
                        )
                    )
                except OSError:
                    pass
        else:
            sp_utilities.wrap_mpi_send(
                newparamstructure[procid], 0, Blockdata["shared_comm"]
            )

        ###fout = open(os.path.join(Tracker["constants"]["masterdir"],"main%03d"%Tracker["mainiteration"],"oldparamstructure","oldparamstructure_%01d_%03d_%03d.json"%(procid,Blockdata["myid"],Tracker["mainiteration"])),'w')
        ###json.dump(newparamstructure[procid], fout)
        ###fout.close()
        newparamstructure[procid] = []
        norm_per_particle[procid] = []
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    del refang, rshifts

    #  DRIVER RESOLUTION ASSESSMENT and RECONSTRUCTION <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    if Tracker["changed_delta"]:
        org_nxinit = Tracker["nxinit"]
        Tracker["nxinit"] = Tracker["constants"]["nnxo"]

    rec3d_make_maps(compute_fsc=True, regularized=True)

    if Tracker["changed_delta"]:
        Tracker["nxinit"] = org_nxinit

    # from sys import exit
    # mpi_finalize()
    # exit()
    #
    #  Change to current params
    partids = [None] * 2
    for procid in range(2):
        partids[procid] = os.path.join(
            Tracker["directory"],
            "chunk_%01d_%03d.txt" % (procid, Tracker["mainiteration"]),
        )
    partstack = [None] * 2
    vol = [None] * 2
    for procid in range(2):
        partstack[procid] = os.path.join(
            Tracker["directory"],
            "params-chunk_%01d_%03d.txt" % (procid, Tracker["mainiteration"]),
        )
    if Blockdata["myid"] == Blockdata["main_node"]:
        pinids = sp_utilities.read_text_file(partids[0]) + sp_utilities.read_text_file(
            partids[1]
        )
        params = sp_utilities.read_text_row(partstack[0]) + sp_utilities.read_text_row(
            partstack[1]
        )

        assert len(pinids) == len(params)

        for i in range(len(pinids)):
            pinids[i] = [pinids[i], params[i]]
        del params
        pinids.sort()

        sp_utilities.write_text_file(
            [pinids[i][0] for i in range(len(pinids))],
            os.path.join(
                Tracker["directory"], "indexes_%03d.txt" % (Tracker["mainiteration"])
            ),
        )
        sp_utilities.write_text_row(
            [pinids[i][1] for i in range(len(pinids))],
            os.path.join(
                Tracker["directory"], "params_%03d.txt" % (Tracker["mainiteration"])
            ),
        )
        del pinids

    for procid in range(2):
        if (
            Blockdata["myid"] == Blockdata["main_node"]
            and Tracker["constants"]["plot_ang_dist"]
        ):
            ### NEEDS TO BE REACTIVATED AFTER THE SYMCLASS CHANGE
            outlier_full = []
            sp_global_def.sxprint(
                "Create angular distribution plot for chunk {0}".format(procid)
            )
            delta = sp_helix_sphire.np.maximum(Tracker["delta"], 3.75)
            exclude = []
            exclude.append(
                [
                    None,
                    os.path.join(Tracker["directory"], "ang_dist_{0}".format(procid)),
                    "",
                ]
            )
            if sp_helix_sphire.np.array(outlier_full).any():
                exclude.append(
                    [
                        outlier_full,
                        os.path.join(
                            Tracker["directory"], "ang_dist_{0}_outlier".format(procid)
                        ),
                        "_outlier",
                    ]
                )

            for exclude_list, dir_name, suffix in exclude:
                sp_utilities.angular_distribution(
                    params_file=partstack[procid],
                    output_folder=dir_name,
                    prefix="ang_dist{0}".format(suffix),
                    method=Tracker["constants"]["angle_method"],
                    pixel_size=1,
                    delta=delta,
                    symmetry=Tracker["constants"]["symmetry"],
                    box_size=Tracker["constants"]["nnxo"],
                    particle_radius=Tracker["constants"]["radius"],
                    dpi=72,
                    do_print=False,
                    exclude=exclude_list,
                )
                sp_utilities.angular_distribution(
                    params_file=partstack[procid],
                    output_folder=dir_name,
                    prefix="full_ang_dist{0}".format(suffix),
                    method=Tracker["constants"]["angle_method"],
                    pixel_size=1,
                    delta=delta,
                    symmetry=Tracker["constants"]["symmetry"] + "_full",
                    box_size=Tracker["constants"]["nnxo"],
                    particle_radius=Tracker["constants"]["radius"],
                    dpi=72,
                    do_print=False,
                    exclude=exclude_list,
                )
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    if Tracker["mainiteration"] == 1:
        acc_rot = acc_trans = 1.0e23
    else:
        if Blockdata["myid"] == Blockdata["main_node"]:
            Blockdata["bckgnoise"] = sp_utilities.get_im(
                os.path.join(Tracker["directory"], "bckgnoise.hdf")
            )
            nnx = Blockdata["bckgnoise"].get_xsize()
            nny = Blockdata["bckgnoise"].get_ysize()
        else:
            nnx = 0
            nny = 0
        nnx = sp_utilities.bcast_number_to_all(nnx)
        nny = sp_utilities.bcast_number_to_all(nny)
        if Blockdata["myid"] != Blockdata["main_node"]:
            Blockdata["bckgnoise"] = sp_utilities.model_blank(nnx, nny, 1, 1.0)
        sp_utilities.bcast_EMData_to_all(
            Blockdata["bckgnoise"],
            Blockdata["myid"],
            source_node=Blockdata["main_node"],
        )

        if Blockdata["myid"] == Blockdata["main_node"]:
            params = sp_utilities.read_text_row(
                os.path.join(
                    Tracker["directory"],
                    "params-chunk_0_%03d.txt" % (Tracker["mainiteration"]),
                )
            ) + sp_utilities.read_text_row(
                os.path.join(
                    Tracker["directory"],
                    "params-chunk_1_%03d.txt" % (Tracker["mainiteration"]),
                )
            )
            li = sp_utilities.read_text_file(
                os.path.join(
                    Tracker["directory"],
                    "chunk_0_%03d.txt" % (Tracker["mainiteration"]),
                )
            ) + sp_utilities.read_text_file(
                os.path.join(
                    Tracker["directory"],
                    "chunk_1_%03d.txt" % (Tracker["mainiteration"]),
                )
            )
            ctfs = EMAN2_cppwrap.EMUtil.get_all_attributes(
                Tracker["constants"]["stack"], "ctf"
            )
            ctfs = [ctfs[i] for i in li]
            particle_groups = sp_utilities.read_text_file(
                os.path.join(
                    Tracker["constants"]["masterdir"],
                    "main000",
                    "particle_groups_0.txt",
                )
            ) + sp_utilities.read_text_file(
                os.path.join(
                    Tracker["constants"]["masterdir"],
                    "main000",
                    "particle_groups_1.txt",
                )
            )
            npart = old_div(500, Blockdata["nproc"]) + 1
            li = list(range(len(ctfs)))
            random.shuffle(li)
            li = li[: npart * Blockdata["nproc"]]
            params = [params[i] for i in li]
            ctfs = [
                [
                    ctfs[i].defocus,
                    ctfs[i].cs,
                    ctfs[i].voltage,
                    ctfs[i].apix,
                    ctfs[i].bfactor,
                    ctfs[i].ampcont,
                    ctfs[i].dfdiff,
                    ctfs[i].dfang,
                ]
                for i in li
            ]
            particle_groups = [particle_groups[i] for i in li]
        else:
            params = 0
            ctfs = 0
            particle_groups = 0
        params = sp_utilities.wrap_mpi_bcast(params, Blockdata["main_node"])
        ctfs = sp_utilities.wrap_mpi_bcast(ctfs, Blockdata["main_node"])
        particle_groups = sp_utilities.wrap_mpi_bcast(
            particle_groups, Blockdata["main_node"]
        )
        # print(" A ",Blockdata["myid"] ,len(params),len(ctfs),len(particle_groups),len(params)/Blockdata["nproc"])
        npart = old_div(len(params), Blockdata["nproc"])
        params = params[Blockdata["myid"] * npart : (Blockdata["myid"] + 1) * npart]
        ctfs = [
            sp_utilities.generate_ctf(ctfs[i])
            for i in range(Blockdata["myid"] * npart, (Blockdata["myid"] + 1) * npart)
        ]
        particle_groups = particle_groups[
            Blockdata["myid"] * npart : (Blockdata["myid"] + 1) * npart
        ]
        Tracker["refvol"] = os.path.join(
            Tracker["directory"], "vol_0_%03d.hdf" % (Tracker["mainiteration"])
        )
        # print(" B ",Blockdata["myid"] ,len(params),len(ctfs),len(particle_groups),npart)
        cerrs(params, ctfs, particle_groups)
        del params, ctfs, particle_groups
        if Blockdata["myid"] == Blockdata["main_node"]:
            sp_utilities.write_text_row(
                [[Tracker["acc_rot"], Tracker["acc_trans"]]],
                os.path.join(
                    Tracker["directory"],
                    "accuracy_%03d.txt" % (Tracker["mainiteration"]),
                ),
            )

    if Blockdata["myid"] == Blockdata["main_node"]:
        anger, shifter = params_changes(
            sp_utilities.read_text_row(
                os.path.join(
                    Tracker["directory"], "params_%03d.txt" % (Tracker["mainiteration"])
                )
            ),
            sp_utilities.read_text_row(
                os.path.join(
                    Tracker["previousoutputdir"],
                    "params_%03d.txt" % (Tracker["mainiteration"] - 1),
                )
            ),
        )
        sp_utilities.write_text_row(
            [[anger, shifter]],
            os.path.join(
                Tracker["directory"],
                "error_thresholds_%03d.txt" % (Tracker["mainiteration"]),
            ),
        )

        # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
        line = ""
        sp_global_def.sxprint(
            line,
            "Average displacements for angular directions = %6.2f degrees; and shifts = %5.1f pixels"
            % (anger, shifter),
        )

        #  Write current Trucker

        if Blockdata["bckgnoise"]:
            Blockdata["bckgnoise"] = "computed"
        dump_tracker_to_json(
            os.path.join(
                Tracker["constants"]["masterdir"],
                "main%03d" % Tracker["mainiteration"],
                "Tracker_%03d.json" % Tracker["mainiteration"],
            ),
            Tracker,
        )
    Tracker["previousoutputdir"] = Tracker["directory"]
    return  # parameters are all passed by Tracker


#
# - "Tracker" (dictionary) object
#   Keeps the current state of option settings and dataset
#   (i.e. particle stack, reference structure, reconstructed volume, and etc)
#   Each iteration is allowed to add new fields/keys
#   if necessary. This happes especially when type of 3D Refinement or metamove changes.
#   Conceptually, each iteration will be associated to a specific Tracker state.
#   Therefore, the list of Tracker state represents the history of process.
#
#   This can be used to restart process from an arbitrary iteration.
#
#


def reduce_shifts(sx, sy, img):
    def rot_matrix(angle):
        angle = sp_helix_sphire.np.radians(angle)
        matrix = sp_helix_sphire.np.array(
            [
                [sp_helix_sphire.np.cos(angle), -sp_helix_sphire.np.sin(angle)],
                [sp_helix_sphire.np.sin(angle), sp_helix_sphire.np.cos(angle)],
            ]
        )
        return matrix

    if Tracker["constants"]["helical_rise"] is not None:
        try:
            rotation_angle = img.get_attr("segment_angle")
        except AttributeError:
            pass
        else:
            rise = old_div(
                Tracker["constants"]["helical_rise"],
                float(Tracker["constants"]["pixel_size"]),
            )
            rise_half = old_div(rise, 2.0)
            point = sp_helix_sphire.np.array([sx, sy])
            rot_point = sp_helix_sphire.np.dot(rot_matrix(rotation_angle), point.T)
            rot_point[0] = ((rot_point[0] + rise_half) % rise) - rise_half
            sx, sy = sp_helix_sphire.np.dot(rot_matrix(rotation_angle).T, rot_point.T)

    return int(round(sx)), int(round(sy))


def get_image_statistics(image, mask, invert):
    if Tracker["constants"]["filament_width"] is None:
        mask2d = mask
    else:
        mask2d = sp_utilities.model_rotated_rectangle2D(
            radius_long=int(
                old_div(sp_helix_sphire.np.sqrt(2 * image.get_xsize() ** 2), 2)
            ),
            radius_short=old_div(
                int(
                    old_div(
                        Tracker["constants"]["filament_width"] * image.get_xsize(),
                        float(Tracker["constants"]["nnxo"]),
                    )
                    + 0.5
                ),
                2,
            ),
            nx=image.get_xsize(),
            ny=image.get_ysize(),
            angle=image.get_attr("segment_angle"),
        )

    return EMAN2_cppwrap.Util.infomask(image, mask2d, invert)


def calculate_prior_values(
    tracker, blockdata, outlier_file, chunk_file, params_file, im_start, im_end, procid
):
    """Calculate the prior values and identify outliers"""

    if not tracker["constants"]["outlier_by"]:
        return [0] * (im_end - im_start), None

    # Print to screen
    if blockdata["myid"] == blockdata["main_node"]:
        sp_global_def.sxprint("Executed successfully: ", "Prior calculation")

    # Calculate outliers
    if blockdata["myid"] == blockdata["main_node"]:
        # Calculate priors
        outliers, new_params, new_index = sp_helix_fundamentals.calculate_priors(
            tracker=Tracker,
            params_file=params_file,
            index_file=chunk_file,
            group_id=Tracker["constants"]["outlier_by"],
            typ="sphire",
            tol_psi=Tracker["prior"]["tol_psi"],
            tol_theta=Tracker["prior"]["tol_theta"],
            tol_filament=Tracker["prior"]["tol_filament"],
            tol_std=Tracker["prior"]["tol_std"],
            tol_mean=Tracker["prior"]["tol_mean"],
            outlier_method=Tracker["prior"]["outlier_method"],
            prior_method=Tracker["prior"]["prior_method"],
            force_outlier=Tracker["prior"]["force_outlier"],
            window_size=Tracker["prior"]["window_size"],
            remove_outlier=Tracker["prior"]["remove_outlier"],
            symclass=Blockdata["symclass"],
            plot=Tracker["prior"]["plot"],
        )

        # Print to screen
        nr_outliers = len(outliers[outliers == 1])
        len_data = len(outliers)
        no_outliers = len_data - nr_outliers
        sp_global_def.sxprint(
            "Chunk {0}: Discarded {1}|{2:.1f}%; Kept {3}|{4:.1f}%; Nr. particles {5}".format(
                procid,
                nr_outliers,
                old_div(100 * nr_outliers, float(len_data)),
                no_outliers,
                old_div(100 * no_outliers, float(len_data)),
                len_data,
            )
        )
        if Tracker["prior"]["force_outlier"]:
            sp_helix_sphire.shutil.copy(params_file, "{0}_old".format(params_file))
            sp_helix_sphire.shutil.copy(new_params, params_file)

        if Tracker["prior"]["apply_prior"]:
            outliers = outliers.tolist()
        else:
            outliers = [0] * len_data

        if (
            old_div(
                Tracker["constants"]["pixel_size"] * Tracker["constants"]["nnxo"],
                float(Tracker["fsc143"]),
            )
            > 10
        ):
            sp_global_def.sxprint(
                "Resolution not sufficient to remove outliers! Do not discard outlier!"
            )
            outliers = [0] * len_data
        elif old_div(100 * no_outliers, float(len_data)) < 50:
            sp_global_def.sxprint(
                "Number of outliers too large! Do not discard outlier!"
            )
            outliers = [0] * len_data

        sp_helix_sphire.np.savetxt(outlier_file, outliers)
    else:
        # Dummy variable
        outliers = 0

    # Distribute outlier list to all processes
    outliers = sp_utilities.bcast_list_to_all(
        outliers, blockdata["myid"], blockdata["main_node"]
    )

    # Get the node specific outlier information
    outliers_node = outliers[im_start:im_end]

    return outliers_node, outliers


def load_tracker_from_json(file_name):
    # Expects this only beeing executed on the main node
    tracker = load_object_from_json(file_name)
    try:
        if tracker["constants"]["stack_prior"] is not None:
            tracker["constants"]["stack_prior_dtype"] = [
                tuple(entry) for entry in tracker["constants"]["stack_prior_dtype"]
            ]
            tracker["constants"]["stack_prior"] = sp_helix_sphire.np.genfromtxt(
                tracker["constants"]["stack_prior"],
                dtype=tracker["constants"]["stack_prior_dtype"],
            )
    except KeyError:
        tracker["constants"]["stack_prior"] = None
    return tracker


def load_object_from_json(file_name):
    with open(file_name, "r") as fin:
        json_object = sp_utilities.convert_json_fromunicode(json.load(fin))
    return json_object


def dump_tracker_to_json(file_name, tracker):
    if tracker["constants"]["stack_prior"] is not None:
        numpy.savetxt(
            os.path.join(tracker["constants"]["masterdir"], "stack_prior.txt"),
            Tracker["constants"]["stack_prior"],
            fmt=Tracker["constants"]["stack_prior_fmt"],
        )
        tracker["constants"]["stack_prior"] = os.path.join(
            tracker["constants"]["masterdir"], "stack_prior.txt"
        )
    dump_object_to_json(file_name, tracker, indent=4)


def dump_object_to_json(file_name, data_object, indent=None):
    with open(file_name, "w") as fout:
        json.dump(data_object, fout, indent=indent)


def prior_stack_fmt(sphire_prior_stack):
    fmt = []
    for entry in sphire_prior_stack.dtype.names:
        if isinstance(sphire_prior_stack[entry][0], (str, numpy.character)):
            fmt.append(
                "% {0}s".format(
                    sp_helix_sphire.np.max(
                        [len(str_entry) for str_entry in sphire_prior_stack[entry]]
                    )
                    + 3
                )
            )

        elif isinstance(sphire_prior_stack[entry][0], (int, numpy.integer)):
            fmt.append(
                "% {0}d".format(
                    len(str(sp_helix_sphire.np.max(sphire_prior_stack[entry]))) + 3
                )
            )

        elif isinstance(sphire_prior_stack[entry][0], (float, numpy.floating)):
            fmt.append(
                "% {0}.6f".format(
                    len(str(sp_helix_sphire.np.max(sphire_prior_stack[entry]))) + 3
                )
            )
        else:
            sp_global_def.sxprint("UNKNOWN", entry, type(entry), sphire_prior_stack[entry][0], type(sphire_prior_stack[entry][0]))
    return " ".join(fmt)


def run():
    global Tracker, Blockdata

    # ------------------------------------------------------------------------------------
    # PARSE COMMAND OPTIONS
    progname = os.path.basename(sys.argv[0])
    usage = (
        progname
        + """ stack  [output_directory]  initial_volume  --radius=particle_radius --symmetry=c1 --initialshifts --inires=25  --mask3D=surface_mask.hdf --function=user_function


	There are five ways to run the program:

1. Standard default run, starts from exhaustive searches, uses initial reference structure
mpirun -np 64 --hostfile four_nodes.txt  sxmeridien.py  bdb:sparx_stack vton1 mask15.hdf --sym=c5  --initialshifts  --radius=120  --mask3D=mask15.hdf    >1ovotn &

2. Restart after the last fully finished iteration, one can change some parameters (MPI settings have to be the same)
mpirun -np 64 --hostfile four_nodes.txt  sxmeridien.py  vton1 --radius=100 >2ovotn &

3. Local refinement, starts from user-provided orientation parameters, delta has to be <= 3.75
mpirun -np 64 --hostfile four_nodes.txt sxmeridien.py --local_refinement bdb:sparx_stack   vton3 --delta=1.875 --xr=2.0  --inires=5.5  --sym=c5  --radius=120  --mask3D=mask15.hdf >5ovotn &

4. Restart of local refinement after the last fully finished iteration, one can change some parameters (MPI settings have to be the same)
mpirun -np 64 --hostfile four_nodes.txt  sxmeridien.py --local_refinement  vton3  --xr=0.6 >6ovotn &

5.  mpirun -np 64 --hostfile four_nodes.txt  sxmeridien.py vton3 --do_final=21

	"""
    )
    parser = optparse.OptionParser(usage, version=sp_global_def.SPARXVERSION)
    parser.add_option(
        "--do_final",
        type="int",
        default=-1,
        help="Do unfiltered odd and even volume 3-D reconstruction from an existing meridien refinement with optional specified iteration",
    )
    parser.add_option(
        "--local_refinement",
        action="store_true",
        default=False,
        help="Perform local refinement starting from user-provided orientation parameters",
    )

    do_final_mode = False
    for q in sys.argv[1:]:
        if q[:10] == "--do_final":
            do_final_mode = True
            break

    do_continuation_mode = False
    for q in sys.argv[1:]:
        if q[:18] == "--local_refinement":
            do_continuation_mode = True
            break

    if (not do_final_mode) and (not do_continuation_mode):
        # case1: standard meridien run
        # case2: restart mode of standard meridien run. Parameters can be altered in the restart run.
        parser.add_option(
            "--radius",
            type="int",
            default=-1,
            help="Outer radius [in pixels] of particles < int(nx/2)-1",
        )
        parser.add_option(
            "--xr",
            type="float",
            default=5.0,
            help="Range for translation search in both directions, search is +/xr (default 5), can be fractional",
        )
        parser.add_option(
            "--ts",
            type="float",
            default=1.0,
            help="Step size of the translation search in both directions, search is within a circle of radius xr on a grid with steps ts, (default 1), can be fractional",
        )
        parser.add_option(
            "--inires",
            type="float",
            default=25.0,
            help="Resolution of the initial_volume volume (default 25A)",
        )
        parser.add_option(
            "--mask3D",
            type="string",
            default=None,
            help="3D mask file (default a sphere with radius (nx/2)-1)",
        )
        parser.add_option(
            "--function",
            type="string",
            default="do_volume_mask",
            help="Vame of the reference preparation function (default do_volume_mask)",
        )
        parser.add_option(
            "--symmetry",
            type="string",
            default="c1",
            help="Point-group symmetry of the refined structure (default c1)",
        )
        parser.add_option(
            "--skip_prealignment",
            action="store_true",
            default=False,
            help="Skip 2-D pre-alignment step: to be used if images are already centered. (default False)",
        )
        parser.add_option(
            "--initialshifts",
            action="store_true",
            default=False,
            help="Use orientation parameters in the input file header to jumpstart the procedure. (default False)",
        )
        parser.add_option(
            "--center_method",
            type="int",
            default=-1,
            help="Method for centering: of average during initial 2D prealignment of data (0 : no centering; -1 : average shift  method;  please see center_2D in utilities.py for methods 1-7) (default -1)",
        )
        parser.add_option(
            "--target_radius",
            type="int",
            default=29,
            help="Target particle radius for 2D prealignment. Images will be shrank/enlarged to this radius (default 29)",
        )
        parser.add_option(
            "--delta",
            type="float",
            default=7.5,
            help="Initial angular sampling step (default 7.5)",
        )
        parser.add_option(
            "--an",
            type="float",
            default=-1.0,
            help="Angular neighborhood for local search",
        )
        parser.add_option("--shake", type="float", default=0.5, help="Shake (0.5)")
        parser.add_option(
            "--small_memory",
            action="store_true",
            default=False,
            help="Data will not be kept in memory if small_memory is true. (default False)",
        )
        parser.add_option(
            "--ccfpercentage",
            type="float",
            default=99.9,
            help="Percentage of the correlation peak area to be included, 0.0 corresponds to hard matching (default 99.9%)",
        )
        parser.add_option(
            "--memory_per_node",
            type="float",
            default=-1.0,
            help="User provided information about memory per node (NOT per CPU) [in GB] (default 2GB*(number of CPUs per node))",
        )
        parser.add_option(
            "--nonorm",
            action="store_true",
            default=False,
            help="Do not apply image norm correction. (default False)",
        )

        parser.add_option(
            "--plot_ang_dist",
            action="store_true",
            default=False,
            help="Plot the angular distribution plot for every iteration. This will take some time for high symmetries. (Default False)",
        )
        parser.add_option(
            "--theta_min",
            type="float",
            default=-1,
            help="Lower limit for the out-of-plane rotation angle. Default is the full range based on the symmetry. (Default -1)",
        )
        parser.add_option(
            "--theta_max",
            type="float",
            default=-1,
            help="Upper limit for the out-of-plane rotation angle.  Default is the full range based on the symmetry. (Default -1)",
        )
        parser.add_option(
            "--howmany",
            type="int",
            default=4,
            help="Upper limit for the out-of-plane rotation angle.  Default is the full range based on the symmetry. (Default -1)",
        )
        parser.add_option(
            "--angle_method",
            type="str",
            default="S",
            help="Even angle creation strategy. Choices: S, P, M. (Default S)",
        )
        parser.add_option(
            "--helical_rise",
            type="float",
            default=None,
            help="Helical rise in angstrom. This is used to limit the shift along the helical axis. (Default None)",
        )
        parser.add_option(
            "--filament_width",
            type="int",
            default=None,
            help="Filament width used to normalize the particles. (Default None)",
        )
        parser.add_option(
            "--chunk_by",
            type="str",
            default="ptcl_source_image",
            help="Group particles by header information. For helical refinement use filament or filament_id if present. (Default ptcl_source_image)",
        )
        parser.add_option(
            "--outlier_by",
            type="str",
            default=None,
            help="Group particles by header information. For helical refinement use filament or filament_id if present. (Default ptcl_source_image)",
        )
        parser.add_option(
            "--outlier_tracker",
            type="str",
            default=None,
            help="Skip stack file creation and load the stack from an existing stack.",
        )
        (options, args) = parser.parse_args(sys.argv[1:])

        if len(args) == 3:
            volinit = args[2]
            masterdir = args[1]
            orgstack = args[0]

        elif len(args) == 2:
            orgstack = args[0]
            volinit = args[1]
            masterdir = ""

        elif len(args) == 1:
            masterdir = args[0]
        else:
            sp_global_def.sxprint("usage: " + usage)
            sp_global_def.sxprint(
                "Please run '" + progname + " -h' for detailed options"
            )
            return 1

        #  Check whether we are restarting the program, in the least main000 should exist, otherwise there is nothing to restart
        keepgoing1 = 1
        keepgoing2 = 1
        restart_flag = 0
        if Blockdata["myid"] == Blockdata["main_node"]:
            if os.path.exists(os.path.join(masterdir, "main000", "Tracker_000.json")):
                if len(args) > 1:
                    keepgoing1 = 0
                restart_flag = 1
            else:
                if len(args) == 1:
                    keepgoing2 = 0
                restart_flag = 0
        restart_flag = sp_utilities.bcast_number_to_all(
            restart_flag,
            source_node=Blockdata["main_node"],
            mpi_comm=mpi.MPI_COMM_WORLD,
        )
        keepgoing1 = sp_utilities.bcast_number_to_all(
            keepgoing1, source_node=Blockdata["main_node"], mpi_comm=mpi.MPI_COMM_WORLD
        )
        keepgoing2 = sp_utilities.bcast_number_to_all(
            keepgoing2, source_node=Blockdata["main_node"], mpi_comm=mpi.MPI_COMM_WORLD
        )
        if keepgoing1 == 0:
            sp_global_def.ERROR(
                "To restart, meridien requires only the name of existing refinement directory.",
                "meridien",
                1,
                Blockdata["myid"],
            )
        if keepgoing2 == 0:
            sp_global_def.ERROR(
                "To start, meridien requires at least the stack name and the name of reference structure",
                "meridien",
                1,
                Blockdata["myid"],
            )
        if restart_flag == 1:
            restart_mode = True
        else:
            restart_mode = False

        # ------------------------------------------------------------------------------------
        # Initialize MPI related variables
        ###  MPI SANITY CHECKES
        if not balanced_processor_load_on_nodes:
            sp_global_def.ERROR(
                "Nodes do not have the same number of CPUs, please check configuration of the cluster.",
                "meridien",
                1,
                Blockdata["myid"],
            )
        if Blockdata["myid"] == Blockdata["main_node"]:
            line = ""
            for a in sys.argv:
                line += a + "  "
            sp_global_def.sxprint(" shell line command ")
            sp_global_def.sxprint(line)
        # ------------------------------------------------------------------------------------
        #  INPUT PARAMETERS
        sp_global_def.BATCH = True
        sp_global_def.MPI = True
        ###  VARIOUS SANITY CHECKES <-----------------------
        if options.memory_per_node < 0.0:
            options.memory_per_node = 2.0 * Blockdata["no_of_processes_per_group"]
        #  For the time being we use all CPUs during refinement
        Blockdata["ncpuspernode"] = Blockdata["no_of_processes_per_group"]
        Blockdata["nsubset"] = Blockdata["ncpuspernode"] * Blockdata["no_of_groups"]
        create_subgroup()
        create_zero_group()

        Blockdata["rkeepf"] = 0.90

        if not restart_mode:  # <<<-------Fresh run
            Prior = {}
            Prior["tol_psi"] = 30
            Prior["tol_theta"] = 15
            Prior["tol_filament"] = 0.2
            Prior["tol_std"] = 1
            Prior["tol_mean"] = 30
            Prior["outlier_method"] = "deg"
            Prior["prior_method"] = "running"
            Prior["force_outlier"] = False
            Prior["remove_outlier"] = False
            Prior["apply_prior"] = False
            Prior["window_size"] = 3
            Prior["plot"] = True

            #  Constant settings of the project
            Constants = {}
            Constants["stack"] = args[0]
            Constants["rs"] = 1
            Constants["radius"] = options.radius
            Constants["an"] = "-1"
            Constants["maxit"] = 1
            Constants["fuse_freq"] = 45  # Now in A, convert to absolute before using
            sym = options.symmetry
            Constants["symmetry"] = sym[0].lower() + sym[1:]
            Constants["npad"] = 1
            Constants["center"] = 0
            Constants["shake"] = options.shake  # move params every iteration
            Constants["CTF"] = True  # internally set
            Constants["mask3D"] = options.mask3D
            Constants["nnxo"] = -1
            Constants["pixel_size"] = None  # read from data
            Constants[
                "inires"
            ] = options.inires  # Now in A, convert to absolute before using
            Constants["refvol"] = volinit
            Constants["masterdir"] = masterdir
            Constants["best"] = 3
            Constants["limit_improvement"] = 1
            Constants[
                "limit_changes"
            ] = 1  # reduce delta by half if both limits are reached simultaneously
            Constants["states"] = [
                "INITIAL",
                "PRIMARY",
                "EXHAUSTIVE",
                "RESTRICTED",
                "LOCAL",
                "FINAL",
            ]  # will add two states, CONINUATION_INITIAL, CONINUATION_PRIMARY
            Constants["user_func"] = options.function
            Constants["hardmask"] = True  # options.hardmask
            Constants["ccfpercentage"] = old_div(options.ccfpercentage, 100.0)
            Constants["expthreshold"] = -10
            Constants[
                "number_of_groups"
            ] = -1  # number of defocus groups, to be set by assign_particles_to_groups
            Constants["nonorm"] = options.nonorm
            Constants["small_memory"] = options.small_memory
            Constants["initialshifts"] = options.initialshifts
            Constants["memory_per_node"] = options.memory_per_node
            Constants["plot_ang_dist"] = options.plot_ang_dist
            Constants["angle_method"] = options.angle_method
            Constants["helical_rise"] = options.helical_rise
            Constants["filament_width"] = options.filament_width
            Constants["outlier_by"] = options.outlier_by
            Constants["do_exhaustive"] = False

            if options.outlier_by is None:
                Constants["stack_prior"] = None
                Constants["stack_prior_fmt"] = None
                Constants["stack_prior_dtype"] = None
            else:
                if Blockdata["myid"] == Blockdata["main_node"]:
                    sp_global_def.sxprint("Import outlier information")
                    if options.outlier_tracker:
                        tmp_tracker = load_tracker_from_json(options.outlier_tracker)
                        Constants["stack_prior"] = tmp_tracker["constants"][
                            "stack_prior"
                        ]
                        del tmp_tracker
                    else:
                        Constants["stack_prior"] = sp_helix_sphire.import_sphire_stack(
                            args[0], options.outlier_by
                        )
                    Constants["stack_prior_fmt"] = prior_stack_fmt(
                        Constants["stack_prior"]
                    )
                    Constants["stack_prior_dtype"] = Constants[
                        "stack_prior"
                    ].dtype.descr
                else:
                    Constants["stack_prior"] = None
                    Constants["stack_prior_fmt"] = None
                    Constants["stack_prior_dtype"] = None

            #
            #  The program will use three different meanings of x-size
            #  nnxo         - original nx of the data, will not be changed
            #  nxinit       - window size used by the program during given iteration,
            #                 will be increased in steps of 32 with the resolution
            #
            #  nxstep       - step by wich window size increases
            #
            # Initialize Tracker Dictionary with input options
            Tracker = {}
            Tracker["constants"] = Constants
            Tracker["prior"] = Prior
            Tracker["maxit"] = Tracker["constants"]["maxit"]

            Tracker["xr"] = options.xr
            Tracker[
                "yr"
            ] = options.xr  # Do not change!  I do not think it is used anywhere
            Tracker["ts"] = options.ts
            Tracker["an"] = "-1"
            Tracker["delta"] = options.delta  # How to decide it
            Tracker["refvol"] = None
            Tracker["nxinit"] = -1  # will be figured in first AI.
            Tracker["nxstep"] = 10
            Tracker["theta_min"] = options.theta_min
            Tracker["theta_max"] = options.theta_max
            Tracker["howmany"] = options.howmany
            #  Resolution in pixels at 0.5 cutoff
            Tracker["currentres"] = -1
            Tracker["fsc143"] = -1
            Tracker["maxfrad"] = -1
            Tracker["no_improvement"] = 0
            Tracker["no_params_changes"] = 0
            Tracker["large_at_Nyquist"] = False
            Tracker["anger"] = 1.0e23
            Tracker["shifter"] = 1.0e23
            Tracker["pixercutoff"] = 2.0
            Tracker["directory"] = ""
            Tracker["previousoutputdir"] = ""
            Tracker["acc_rot"] = 0.0
            Tracker["acc_trans"] = 0.0
            Tracker["avgvaradj"] = [1.0, 1.0]  # This has to be initialized to 1.0 !!
            Tracker["mainiteration"] = 0
            Tracker["lentop"] = 2000
            Tracker["state"] = Tracker["constants"]["states"][0]
            Tracker["nima_per_chunk"] = [0, 0]
            ###<<<----state
            Tracker["bestres"] = 0
            Tracker["bestres_143"] = 0
            Blockdata["bckgnoise"] = None
            Blockdata["accumulatepw"] = [[], []]

            # ------------------------------------------------------------------------------------
            # Get the pixel size; if none, set to 1.0, and the original image size
            if Blockdata["myid"] == Blockdata["main_node"]:
                # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
                line = ""
                sp_global_def.sxprint(line, "INITIALIZATION OF MERIDIEN")
                a = sp_utilities.get_im(orgstack)
                nnxo = a.get_xsize()
                if Tracker["constants"]["CTF"]:
                    i = a.get_attr("ctf")
                    pixel_size = i.apix
                    fq = int(
                        old_div(pixel_size * nnxo, Tracker["constants"]["fuse_freq"])
                        + 0.5
                    )
                else:
                    pixel_size = Tracker["constants"]["pixel_size"]
                    #  No pixel size, fusing computed as 5 Fourier pixels
                    fq = 5
                del a
            else:
                nnxo = 0
                fq = 0
                pixel_size = 1.0

            #  Object to handle symmetries, for now only used by oct
            # Blockdata["parsesym"] = parsesym(Tracker["constants"]["symmetry"])
            #  Initialize symmetry
            Blockdata["symclass"] = sp_fundamentals.symclass(
                Tracker["constants"]["symmetry"]
            )

            nnxo = sp_utilities.bcast_number_to_all(
                nnxo, source_node=Blockdata["main_node"]
            )
            if nnxo < 0:
                sp_global_def.ERROR(
                    "Incorrect image size  ", "meridien", 1, Blockdata["myid"]
                )
            pixel_size = sp_utilities.bcast_number_to_all(
                pixel_size, source_node=Blockdata["main_node"]
            )
            fq = sp_utilities.bcast_number_to_all(
                fq, source_node=Blockdata["main_node"]
            )
            Tracker["constants"]["nnxo"] = nnxo
            Tracker["constants"]["pixel_size"] = pixel_size
            Tracker["constants"]["fuse_freq"] = fq
            del fq, nnxo, pixel_size
            # Resolution is always in full size image pixel units.
            Tracker["constants"]["inires"] = int(
                old_div(
                    Tracker["constants"]["nnxo"] * Tracker["constants"]["pixel_size"],
                    Tracker["constants"]["inires"],
                )
                + 0.5
            )
            Tracker["currentres"] = Tracker["constants"]["inires"]
            Tracker["fsc143"] = Tracker["constants"]["inires"]

            checking_flag = 1
            ###  VARIOUS SANITY CHECKES
            if options.initialshifts:
                options.skip_prealignment = (
                    True
                )  # No prealignment if initial shifts are set
            if Blockdata["myid"] == Blockdata["main_node"]:
                if Tracker["constants"]["mask3D"] and (
                    not os.path.exists(Tracker["constants"]["mask3D"])
                ):
                    checking_flag = 0
            checking_flag = sp_utilities.bcast_number_to_all(
                checking_flag,
                source_node=Blockdata["main_node"],
                mpi_comm=mpi.MPI_COMM_WORLD,
            )
            if checking_flag == 0:
                sp_global_def.ERROR(
                    "mask3D file does  not exists ", "meridien", 1, Blockdata["myid"]
                )

            if old_div(options.xr, options.ts) < 1.0:
                sp_global_def.ERROR(
                    "Incorrect translational searching settings, search range cannot be smaller than translation step ",
                    "meridien",
                    1,
                    Blockdata["myid"],
                )
            if (
                2 * (Tracker["currentres"] + Tracker["nxstep"])
                > Tracker["constants"]["nnxo"]
            ):
                sp_global_def.ERROR(
                    "Image size less than what would follow from the initial resolution provided %d  %d  %d"
                    % (
                        Tracker["currentres"],
                        Tracker["nxstep"],
                        2 * (Tracker["currentres"] + Tracker["nxstep"]),
                    ),
                    "sxmeridien",
                    1,
                    Blockdata["myid"],
                )

            if Tracker["constants"]["radius"] < 1:
                Tracker["constants"]["radius"] = (
                    old_div(Tracker["constants"]["nnxo"], 2) - 2
                )
            elif (2 * Tracker["constants"]["radius"] + 2) > Tracker["constants"][
                "nnxo"
            ]:
                sp_global_def.ERROR(
                    "Particle radius set too large!", "sxmeridien", 1, Blockdata["myid"]
                )
            ###<-----end of sanity check <----------------------
            ###<<<----------------------------- parse program

            # ------------------------------------------------------------------------------------
            #  MASTER DIRECTORY
            if Blockdata["myid"] == Blockdata["main_node"]:
                if masterdir == "":
                    timestring = time.strftime("_%d_%b_%Y_%H_%M_%S", time.localtime())
                    masterdir = "master" + timestring
                    li = len(masterdir)
                    os.makedirs(masterdir)
                    keepchecking = 0
                else:
                    if not os.path.exists(masterdir):
                        os.makedirs(masterdir)
                    li = 0
                    keepchecking = 1
                sp_global_def.write_command(masterdir)
            else:
                li = 0
                keepchecking = 1

            li = mpi.mpi_bcast(
                li, 1, mpi.MPI_INT, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )[0]

            if li > 0:
                masterdir = mpi.mpi_bcast(
                    masterdir,
                    li,
                    mpi.MPI_CHAR,
                    Blockdata["main_node"],
                    mpi.MPI_COMM_WORLD,
                )
                masterdir = string.join(masterdir, "")

            Tracker["constants"]["masterdir"] = masterdir
            if Blockdata["myid"] == Blockdata["main_node"]:
                print_dict(Tracker["constants"], "Permanent settings of meridien")
                print_dict(Blockdata, "MPI settings of meridien")

            # Initialization of orgstack
            Tracker["constants"]["stack"] = orgstack
            if Blockdata["myid"] == Blockdata["main_node"]:
                total_stack = EMAN2_cppwrap.EMUtil.get_image_count(
                    Tracker["constants"]["stack"]
                )
            else:
                total_stack = 0
            total_stack = sp_utilities.bcast_number_to_all(
                total_stack, source_node=Blockdata["main_node"]
            )
            # ------------------------------------------------------------------------------------
            #  	Fresh start INITIALIZATION
            initdir = os.path.join(Tracker["constants"]["masterdir"], "main000")

            partids = os.path.join(initdir, "indexes_000.txt")
            #### add prealignment like in isac

            #########################################################################################################################
            # prepare parameters to call calculate_2d_params_for_centering

            radi = options.radius
            target_radius = options.target_radius
            # target_nx = options.target_nx
            center_method = options.center_method
            if radi < 1:
                sp_global_def.ERROR(
                    "Particle radius has to be provided!",
                    "2d prealignment",
                    1,
                    Blockdata["myid"],
                )

            nxrsteps = 4

            init2dir = os.path.join(masterdir, "2dalignment")

            ##########################################################################################################################
            # put all parameters in a dictionary
            kwargs = dict()

            kwargs["init2dir"] = init2dir
            kwargs["myid"] = Blockdata["myid"]
            kwargs["main_node"] = Blockdata["main_node"]
            kwargs["number_of_images_in_stack"] = total_stack
            kwargs["nproc"] = Blockdata["nproc"]

            kwargs["target_radius"] = target_radius
            # kwargs["target_nx"] = target_nx
            kwargs["radi"] = radi

            kwargs["center_method"] = center_method

            kwargs["nxrsteps"] = nxrsteps

            kwargs["command_line_provided_stack_filename"] = Tracker["constants"][
                "stack"
            ]

            # kwargs["masterdir"] = masterdir

            kwargs["options_skip_prealignment"] = options.skip_prealignment
            kwargs["options_CTF"] = True

            kwargs["mpi_comm"] = mpi.MPI_COMM_WORLD
            #################################################################################################################################################################
            if (Blockdata["myid"] == Blockdata["main_node"]) and (
                not options.skip_prealignment
            ):
                # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
                line = ""
                sp_global_def.sxprint(line, "2D pre-alignment step")
            ## only the master has the parameters
            params2d = calculate_2d_params_for_centering(kwargs)
            del kwargs
            if (Blockdata["myid"] == Blockdata["main_node"]) and (
                not options.skip_prealignment
            ):
                # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
                line = ""
                sp_global_def.sxprint(line, "2D pre-alignment completed")

            if Blockdata["myid"] == Blockdata["main_node"]:
                os.makedirs(initdir)
                sp_utilities.write_text_file(list(range(total_stack)), partids)
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            #  store params
            partids = [None] * 2
            for procid in range(2):
                partids[procid] = os.path.join(initdir, "chunk_%01d_000.txt" % procid)
            partstack = [None] * 2
            for procid in range(2):
                partstack[procid] = os.path.join(
                    initdir, "params-chunk_%01d_000.txt" % procid
                )
            if Blockdata["myid"] == Blockdata["main_node"]:
                l1, l2 = assign_particles_to_groups(
                    minimum_group_size=10, name_tag=options.chunk_by
                )
                sp_utilities.write_text_file(l1, partids[0])
                sp_utilities.write_text_file(l2, partids[1])
                if options.initialshifts:
                    tp_list = EMAN2_cppwrap.EMUtil.get_all_attributes(
                        Tracker["constants"]["stack"], "xform.projection"
                    )
                    for i in range(len(tp_list)):
                        dp = tp_list[i].get_params("spider")
                        tp_list[i] = [
                            dp["phi"],
                            dp["theta"],
                            dp["psi"],
                            -dp["tx"],
                            -dp["ty"],
                            0.0,
                            1.0,
                        ]
                    sp_utilities.write_text_row(
                        tp_list, os.path.join(initdir, "params_000.txt")
                    )
                    sp_utilities.write_text_row([tp_list[i] for i in l1], partstack[0])
                    sp_utilities.write_text_row([tp_list[i] for i in l2], partstack[1])
                    # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
                    line = ""
                    sp_global_def.sxprint(
                        line,
                        "Executed successfully: Imported initial parameters from the input stack",
                    )
                    del tp_list

                else:
                    sp_utilities.write_text_row(
                        [
                            [0, 0, 0, params2d[i][1], params2d[i][2], 0.0, 1.0]
                            for i in l1
                        ],
                        partstack[0],
                    )
                    sp_utilities.write_text_row(
                        [
                            [0, 0, 0, params2d[i][1], params2d[i][2], 0.0, 1.0]
                            for i in l2
                        ],
                        partstack[1],
                    )
                    sp_utilities.write_text_row(
                        [
                            [0, 0, 0, params2d[i][1], params2d[i][2], 0.0, 1.0]
                            for i in range(len(l1) + len(l2))
                        ],
                        os.path.join(initdir, "params_000.txt"),
                    )

                del l1, l2

                # Create reference models for each particle group
                if Tracker["constants"]["mask3D"] == None:
                    viv = sp_filter.filt_table(
                        sp_morphology.cosinemask(
                            sp_utilities.get_im(volinit),
                            radius=Tracker["constants"]["radius"],
                        ),
                        [1.0] * Tracker["constants"]["inires"]
                        + [0.5]
                        + [0.0] * Tracker["constants"]["nnxo"],
                    )
                else:
                    viv = sp_filter.filt_table(
                        sp_utilities.get_im(volinit)
                        * sp_utilities.get_im(Tracker["constants"]["mask3D"]),
                        [1.0] * Tracker["constants"]["inires"]
                        + [0.5]
                        + [0.0] * Tracker["constants"]["nnxo"],
                    )
                # make a copy of original reference model for this particle group (procid)
                for procid in range(2):
                    viv.write_image(
                        os.path.join(
                            initdir,
                            "vol_%01d_%03d.hdf" % (procid, Tracker["mainiteration"]),
                        )
                    )
                del viv
            else:
                Tracker["nima_per_chunk"] = [0, 0]
            Tracker["nima_per_chunk"][0] = sp_utilities.bcast_number_to_all(
                Tracker["nima_per_chunk"][0], Blockdata["main_node"]
            )
            Tracker["nima_per_chunk"][1] = sp_utilities.bcast_number_to_all(
                Tracker["nima_per_chunk"][1], Blockdata["main_node"]
            )
            Tracker["constants"]["number_of_groups"] = sp_utilities.bcast_number_to_all(
                Tracker["constants"]["number_of_groups"], Blockdata["main_node"]
            )
            del params2d

            mainiteration = 0
            Tracker["mainiteration"] = mainiteration
            if Blockdata["myid"] == Blockdata["main_node"]:
                dump_tracker_to_json(
                    os.path.join(
                        Tracker["constants"]["masterdir"],
                        "main%03d" % Tracker["mainiteration"],
                        "Tracker_%03d.json" % Tracker["mainiteration"],
                    ),
                    Tracker,
                )

        else:  # simple restart, at least main000 is completed. Otherwise no need restart

            Blockdata["bckgnoise"] = None
            Blockdata["accumulatepw"] = [[], []]
            initdir = os.path.join(masterdir, "main000")
            keepchecking = 1
            if Blockdata["myid"] == Blockdata["main_node"]:
                Tracker = load_tracker_from_json(
                    os.path.join(initdir, "Tracker_000.json")
                )
                print_dict(
                    Tracker["constants"],
                    "Permanent settings of the original run recovered from main000",
                )
                stack_prior = Tracker["constants"]["stack_prior"]
                Tracker["constants"]["stack_prior"] = None
                sp_global_def.write_command(Tracker["constants"]["masterdir"])
            else:
                Tracker = None
                stack_prior = None
            Tracker = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"])
            Tracker["constants"]["stack_prior"] = stack_prior
            mainiteration = 0
            Tracker["mainiteration"] = mainiteration

        Tracker["previousoutputdir"] = initdir

        # ------------------------------------------------------------------------------------
        #  MAIN ITERATION
        projdata = [[sp_utilities.model_blank(1, 1)], [sp_utilities.model_blank(1, 1)]]
        oldparams = [[], []]
        currentparams = [[], []]
        original_data = [None, None]
        keepgoing = 1
        while keepgoing:
            Tracker["previousoutputdir"] = os.path.join(
                Tracker["constants"]["masterdir"], "main%03d" % Tracker["mainiteration"]
            )
            mainiteration += 1
            Tracker["mainiteration"] = mainiteration
            Tracker["directory"] = os.path.join(
                Tracker["constants"]["masterdir"], "main%03d" % Tracker["mainiteration"]
            )
            doit, keepchecking = checkstep(Tracker["directory"], keepchecking)
            if not doit:
                li = True
                doit2, keepchecking2 = checkstep(
                    os.path.join(
                        Tracker["directory"],
                        "Tracker_%03d.json" % Tracker["mainiteration"],
                    ),
                    li,
                )
                if doit2:
                    doit = True
                    if Blockdata["myid"] == Blockdata["main_node"]:
                        sp_helix_sphire.shutil.rmtree(Tracker["directory"])
                    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            if doit:
                if Blockdata["myid"] == Blockdata["main_node"]:
                    Tracker = load_tracker_from_json(
                        os.path.join(
                            Tracker["previousoutputdir"],
                            "Tracker_%03d.json" % (Tracker["mainiteration"] - 1),
                        )
                    )
                    stack_prior = Tracker["constants"]["stack_prior"]
                    Tracker["constants"]["stack_prior"] = None
                    #  It has to be repeated here as Tracker is from previous iteration, I see no other way.
                    Tracker["previousoutputdir"] = os.path.join(
                        Tracker["constants"]["masterdir"],
                        "main%03d" % Tracker["mainiteration"],
                    )
                    Tracker["mainiteration"] = mainiteration
                    Tracker["directory"] = os.path.join(
                        Tracker["constants"]["masterdir"],
                        "main%03d" % Tracker["mainiteration"],
                    )
                else:
                    Tracker = None
                    stack_prior = None
                Tracker = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"])
                Tracker["constants"]["stack_prior"] = stack_prior
                ###
                if restart_mode:
                    update_tracker(sys.argv[1:])
                    update_memory_estimation()
                    restart_mode = False

                # prepare names of input file names, they are in main directory,
                #   log subdirectories contain outputs from specific refinements
                partids = [None] * 2
                for procid in range(2):
                    partids[procid] = os.path.join(
                        Tracker["previousoutputdir"],
                        "chunk_%01d_%03d.txt" % (procid, Tracker["mainiteration"] - 1),
                    )
                partstack = [None] * 2
                for procid in range(2):
                    partstack[procid] = os.path.join(
                        Tracker["previousoutputdir"],
                        "params-chunk_%01d_%03d.txt"
                        % (procid, Tracker["mainiteration"] - 1),
                    )

                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

                if Tracker["mainiteration"] == 1:
                    fff = None
                    anger = 1.0e9
                    shifter = 1.0e9
                else:
                    if Blockdata["myid"] == Blockdata["main_node"]:
                        fff = sp_utilities.read_text_file(
                            os.path.join(
                                Tracker["previousoutputdir"],
                                "driver_%03d.txt" % (Tracker["mainiteration"] - 1),
                            )
                        )
                        [anger, shifter] = sp_utilities.read_text_row(
                            os.path.join(
                                Tracker["previousoutputdir"],
                                "error_thresholds_%03d.txt"
                                % (Tracker["mainiteration"] - 1),
                            )
                        )[0]
                    else:
                        fff = []
                        anger = 0.0
                        shifter = 0.0
                    fff = sp_utilities.bcast_list_to_all(
                        fff, Blockdata["myid"], source_node=Blockdata["main_node"]
                    )
                    anger = sp_utilities.bcast_number_to_all(
                        anger, source_node=Blockdata["main_node"]
                    )
                    shifter = sp_utilities.bcast_number_to_all(
                        shifter, source_node=Blockdata["main_node"]
                    )

                keepgoing = AI(
                    fff, anger, shifter, Blockdata["myid"] == Blockdata["main_node"]
                )

                if keepgoing == 1:  # not converged
                    if Blockdata["myid"] == Blockdata["main_node"]:
                        if Tracker["mainiteration"] > 1:
                            # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
                            line = ""
                            sp_global_def.sxprint(
                                line,
                                "Resolution achieved in ITERATION  #%2d: %3d/%3d pixels, %5.2fA/%5.2fA."
                                % (
                                    Tracker["mainiteration"] - 1,
                                    Tracker["currentres"],
                                    Tracker["fsc143"],
                                    old_div(
                                        Tracker["constants"]["pixel_size"]
                                        * Tracker["constants"]["nnxo"],
                                        float(Tracker["currentres"]),
                                    ),
                                    old_div(
                                        Tracker["constants"]["pixel_size"]
                                        * Tracker["constants"]["nnxo"],
                                        float(Tracker["fsc143"]),
                                    ),
                                ),
                            )
                        sp_global_def.sxprint("\n\n\n\n")
                        # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
                        line = ""
                        sp_global_def.sxprint(
                            line,
                            "ITERATION  #%2d. Current state: %14s, nxinit: %3d, delta: %9.4f, xr: %9.4f, ts: %9.4f"
                            % (
                                Tracker["mainiteration"],
                                Tracker["state"],
                                Tracker["nxinit"],
                                Tracker["delta"],
                                Tracker["xr"],
                                Tracker["ts"],
                            ),
                        )
                    # print("RACING  A ",Blockdata["myid"])
                    li = True
                    doit2, keepchecking2 = checkstep(Tracker["directory"], li)
                    if Blockdata["myid"] == Blockdata["main_node"] and doit2:
                        cmd = "{} {}".format("mkdir -p", Tracker["directory"])
                        junk = sp_utilities.cmdexecute(cmd)
                        cmd = "{} {}".format(
                            "mkdir -p",
                            os.path.join(Tracker["directory"], "oldparamstructure"),
                        )
                        junk = sp_utilities.cmdexecute(cmd)
                    if not doit2:
                        sp_global_def.ERROR(
                            "There was a gap in main directories, program cannot proceed",
                            "sxmeridien",
                            1,
                            Blockdata["myid"],
                        )

                    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

                    refinement_one_iteration(
                        partids,
                        partstack,
                        original_data,
                        oldparams,
                        projdata,
                        general_mode=True,
                        continuation_mode=False,
                    )

                    # 	sxprint("  MOVING  ON --------------------------------------------------------------------")
                else:  # converged, do final
                    if Blockdata["subgroup_myid"] > -1:
                        mpi.mpi_comm_free(Blockdata["subgroup_comm"])
                    if Blockdata["myid"] == Blockdata["main_node"]:
                        sp_global_def.sxprint(
                            line,
                            "Resolution achieved in ITERATION  #%2d: %3d/%3d pixels, %5.2fA/%5.2fA."
                            % (
                                Tracker["mainiteration"] - 1,
                                Tracker["currentres"],
                                Tracker["fsc143"],
                                old_div(
                                    Tracker["constants"]["pixel_size"]
                                    * Tracker["constants"]["nnxo"],
                                    float(Tracker["currentres"]),
                                ),
                                old_div(
                                    Tracker["constants"]["pixel_size"]
                                    * Tracker["constants"]["nnxo"],
                                    float(Tracker["fsc143"]),
                                ),
                            ),
                        )
                        sp_global_def.sxprint("\n\n\n\n")

                        sp_global_def.sxprint(
                            "The iteration that contains the best resolution is %d"
                            % Tracker["constants"]["best"]
                        )
                        if Tracker["constants"]["best"] == 2:
                            sp_global_def.sxprint(
                                "No resolution improvement in refinement "
                            )
                        dump_tracker_to_json(
                            os.path.join(masterdir, "Tracker_final.json"), Tracker
                        )
                    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

                    newparamstructure = [[], []]
                    projdata = [
                        [sp_utilities.model_blank(1, 1)],
                        [sp_utilities.model_blank(1, 1)],
                    ]
                    original_data = [None, None]
                    oldparams = [[], []]
                    Blockdata["accumulatepw"] = [None, None]
                    # if Tracker["constants"]["memory_per_node"] <0.0: Tracker["constants"]["memory_per_node"] = 2.0*Blockdata["no_of_processes_per_group"]
                    recons3d_final(
                        Tracker["constants"]["masterdir"],
                        Tracker["constants"]["best"],
                        Tracker["constants"]["memory_per_node"],
                    )

            #  End of if doit
        #   end of while
        mpi.mpi_finalize()
        exit()

    elif do_continuation_mode:
        # case1: local meridien run using parameters stored in headers
        # case2: restart mode of standard meridien run. Parameters can be altered in the restart run.
        parser.add_option(
            "--radius",
            type="int",
            default=-1,
            help="Outer radius [in pixels] of particles < int(nx/2)-1",
        )
        parser.add_option(
            "--xr",
            type="float",
            default=5.0,
            help="Range for translation search in both directions, search is +/xr (default 5), can be fractional",
        )
        parser.add_option(
            "--ts",
            type="float",
            default=1.0,
            help="Step size of the translation search in both directions, search is within a circle of radius xr on a grid with steps ts, (default 1), can be fractional",
        )
        parser.add_option(
            "--inires",
            type="float",
            default=-1.0,
            help="Resolution of the initial_volume volume (default 25A)",
        )
        parser.add_option(
            "--mask3D",
            type="string",
            default=None,
            help="3D mask file (default a sphere with radius (nx/2)-1)",
        )
        parser.add_option(
            "--function",
            type="string",
            default="do_volume_mask",
            help="Vame of the reference preparation function (default do_volume_mask)",
        )
        parser.add_option(
            "--symmetry",
            type="string",
            default="c1",
            help="Point-group symmetry of the refined structure (default c1)",
        )
        parser.add_option(
            "--delta",
            type="float",
            default=3.75,
            help="Initial angular sampling step (default 7.5)",
        )
        parser.add_option(
            "--an",
            type="float",
            default=-1.0,
            help="Angular neighborhood for local search",
        )
        parser.add_option("--shake", type="float", default=0.5, help="Shake (0.5)")
        parser.add_option(
            "--small_memory",
            action="store_true",
            default=False,
            help="Data will not be kept in memory if small_memory is true. (default False)",
        )
        parser.add_option(
            "--ccfpercentage",
            type="float",
            default=99.9,
            help="Percentage of the correlation peak area to be included, 0.0 corresponds to hard matching (default 99.9%)",
        )
        parser.add_option(
            "--nonorm",
            action="store_true",
            default=False,
            help="Do not apply image norm correction. (default False)",
        )
        parser.add_option(
            "--memory_per_node",
            type="float",
            default=-1.0,
            help="User provided information about memory per node (NOT per CPU) [in GB] (default 2GB*(number of CPUs per node))",
        )
        parser.add_option(
            "--skip_primary",
            action="store_true",
            default=False,
            help="Skip the primary step and go directly to RESTRICTED.",
        )
        (options, args) = parser.parse_args(sys.argv[1:])

        if len(args) == 2:
            masterdir = args[1]
            orgstack = args[0]

        elif len(args) == 1:
            masterdir = args[0]

        else:
            sp_global_def.sxprint("usage: " + usage)
            sp_global_def.sxprint(
                "Please run '" + progname + " -h' for detailed options"
            )
            return 1
        #  Check whether we are restarting the program, in the least main000 should exist, otherwise there is nothing to restart
        keepgoing1 = 1
        keepgoing2 = 1
        restart_flag = 0
        if Blockdata["myid"] == Blockdata["main_node"]:
            if os.path.exists(os.path.join(masterdir, "main000", "Tracker_000.json")):
                if len(args) > 1:
                    keepgoing1 = 0
                restart_flag = 1
            else:
                if len(args) == 1:
                    keepgoing2 = 0
                restart_flag = 0
        restart_flag = sp_utilities.bcast_number_to_all(
            restart_flag,
            source_node=Blockdata["main_node"],
            mpi_comm=mpi.MPI_COMM_WORLD,
        )
        keepgoing1 = sp_utilities.bcast_number_to_all(
            keepgoing1, source_node=Blockdata["main_node"], mpi_comm=mpi.MPI_COMM_WORLD
        )
        keepgoing2 = sp_utilities.bcast_number_to_all(
            keepgoing2, source_node=Blockdata["main_node"], mpi_comm=mpi.MPI_COMM_WORLD
        )

        if keepgoing1 == 0:
            sp_global_def.ERROR(
                "To restart, meridien requires only the name of existing refinement directory.",
                "meridien local",
                1,
                Blockdata["myid"],
            )

        if keepgoing2 == 0:
            sp_global_def.ERROR(
                "To start, meridien requires at least the stack name and the name of reference structure",
                "meridien local",
                1,
                Blockdata["myid"],
            )

        if restart_flag == 1:
            restart_mode = True
        else:
            restart_mode = False

        # ------------------------------------------------------------------------------------
        # Initialize MPI related variables
        ###  MPI SANITY CHECKES
        if not balanced_processor_load_on_nodes:
            sp_global_def.ERROR(
                "Nodes do not have the same number of CPUs, please check configuration of the cluster.",
                "meridien",
                1,
                Blockdata["myid"],
            )
        if Blockdata["myid"] == Blockdata["main_node"]:
            line = ""
            for a in sys.argv:
                line += a + "  "
            sp_global_def.sxprint(" shell line command ")
            sp_global_def.sxprint(line)
        # ------------------------------------------------------------------------------------
        #  INPUT PARAMETERS
        sp_global_def.BATCH = True
        sp_global_def.MPI = True
        ###  VARIOUS SANITY CHECKES <-----------------------
        if options.delta > 3.75:
            sp_global_def.ERROR(
                "Local searches requested, delta cannot be larger than 3.73.",
                "meridien",
                1,
                Blockdata["myid"],
            )
        if options.memory_per_node < 0.0:
            options.memory_per_node = 2.0 * Blockdata["no_of_processes_per_group"]
        #  For the time being we use all CPUs during refinement
        Blockdata["ncpuspernode"] = Blockdata["no_of_processes_per_group"]
        Blockdata["nsubset"] = Blockdata["ncpuspernode"] * Blockdata["no_of_groups"]

        create_subgroup()
        create_zero_group()

        Blockdata["rkeepf"] = 0.90

        if not restart_mode:  # <<<-------Fresh run
            #  Constant settings of the project
            Constants = {}
            Constants["stack"] = args[0]
            Constants["rs"] = 1
            Constants["radius"] = options.radius
            Constants["an"] = "-1"
            Constants["maxit"] = 1
            Constants["fuse_freq"] = 45  # Now in A, convert to absolute before using
            sym = options.symmetry
            Constants["symmetry"] = sym[0].lower() + sym[1:]
            Constants["npad"] = 1
            Constants["center"] = 0
            Constants["shake"] = options.shake  # move params every iteration
            Constants["CTF"] = True  # internally set
            Constants["mask3D"] = options.mask3D
            Constants["nnxo"] = -1
            Constants["pixel_size"] = None  # read from data
            Constants[
                "inires"
            ] = options.inires  # Now in A, convert to absolute before using
            # Constants["refvol"]            			= volinit
            Constants["masterdir"] = masterdir
            Constants["best"] = 3
            Constants["limit_improvement"] = 1
            Constants[
                "limit_changes"
            ] = 1  # reduce delta by half if both limits are reached simultaneously
            Constants["states"] = ["PRIMARY", "RESTRICTED", "FINAL"]
            Constants["user_func"] = options.function
            Constants["hardmask"] = True  # options.hardmask
            Constants["ccfpercentage"] = old_div(options.ccfpercentage, 100.0)
            Constants["expthreshold"] = -10
            Constants[
                "number_of_groups"
            ] = -1  # number of defocus groups, to be set by assign_particles_to_groups
            Constants["nonorm"] = options.nonorm
            Constants["small_memory"] = options.small_memory
            Constants["memory_per_node"] = options.memory_per_node

            #
            #  The program will use three different meanings of x-size
            #  nnxo         - original nx of the data, will not be changed
            #  nxinit       - window size used by the program during given iteration,
            #                 will be increased in steps of 32 with the resolution
            #
            #  nxstep       - step by wich window size increases
            #
            # Initialize Tracker Dictionary with input options
            Tracker = {}
            Tracker["constants"] = Constants
            Tracker["maxit"] = Tracker["constants"]["maxit"]

            Tracker["xr"] = options.xr
            Tracker[
                "yr"
            ] = options.xr  # Do not change!  I do not think it is used anywhere
            Tracker["ts"] = options.ts
            Tracker["an"] = "-1"
            Tracker["delta"] = options.delta  # How to decide it
            Tracker["refvol"] = None
            Tracker["nxinit"] = -1  # will be figured in first AI.
            Tracker["nxstep"] = 10
            #  Resolution in pixels at 0.5 cutoff
            Tracker["currentres"] = -1
            Tracker["fsc143"] = -1
            Tracker["maxfrad"] = -1
            Tracker["no_improvement"] = 0
            Tracker["no_params_changes"] = 0
            Tracker["large_at_Nyquist"] = False
            Tracker["anger"] = 1.0e23
            Tracker["shifter"] = 1.0e23
            Tracker["pixercutoff"] = 2.0
            Tracker["directory"] = ""
            Tracker["previousoutputdir"] = ""
            Tracker["acc_rot"] = 0.0
            Tracker["acc_trans"] = 0.0
            Tracker["avgvaradj"] = [1.0, 1.0]  # This has to be initialized to 1.0 !!
            Tracker["mainiteration"] = 0
            Tracker["lentop"] = 2000
            Tracker["state"] = Tracker["constants"]["states"][options.skip_primary]
            Tracker["nima_per_chunk"] = [0, 0]
            ###<<<----state
            Tracker["bestres"] = 0
            Tracker["bestres_143"] = 0
            Blockdata["bckgnoise"] = None
            Blockdata["accumulatepw"] = [[], []]
            Tracker["changed_delta"] = False

            # ------------------------------------------------------------------------------------
            # Get the pixel size; if none, set to 1.0, and the original image size
            if Blockdata["myid"] == Blockdata["main_node"]:
                # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
                line = ""
                sp_global_def.sxprint(line, "INITIALIZATION OF LOCAL MERIDIEN")
                a = sp_utilities.get_im(orgstack)
                nnxo = a.get_xsize()
                if Tracker["constants"]["CTF"]:
                    i = a.get_attr("ctf")
                    pixel_size = i.apix
                    fq = int(
                        old_div(pixel_size * nnxo, Tracker["constants"]["fuse_freq"])
                        + 0.5
                    )
                else:
                    pixel_size = Tracker["constants"]["pixel_size"]
                    #  No pixel size, fusing computed as 5 Fourier pixels
                    fq = 5
                del a
            else:
                nnxo = 0
                fq = 0
                pixel_size = 1.0

            #  Object to handle symmetries, for now only used by oct
            # Blockdata["parsesym"] = parsesym(Tracker["constants"]["symmetry"])
            #  Initialize symmetry
            Blockdata["symclass"] = sp_fundamentals.symclass(
                Tracker["constants"]["symmetry"]
            )

            nnxo = sp_utilities.bcast_number_to_all(
                nnxo, source_node=Blockdata["main_node"]
            )
            if nnxo < 0:
                sp_global_def.ERROR(
                    "Incorrect image size  ", "meridien", 1, Blockdata["myid"]
                )
            pixel_size = sp_utilities.bcast_number_to_all(
                pixel_size, source_node=Blockdata["main_node"]
            )
            fq = sp_utilities.bcast_number_to_all(
                fq, source_node=Blockdata["main_node"]
            )
            Tracker["constants"]["nnxo"] = nnxo
            Tracker["constants"]["pixel_size"] = pixel_size
            Tracker["constants"]["fuse_freq"] = fq
            del fq, nnxo, pixel_size
            # Resolution is always in full size image pixel units.
            # HOHO
            if Tracker["constants"]["inires"] > 0.0:
                Tracker["constants"]["inires"] = int(
                    old_div(
                        Tracker["constants"]["nnxo"]
                        * Tracker["constants"]["pixel_size"],
                        Tracker["constants"]["inires"],
                    )
                    + 0.5
                )
            Tracker["currentres"] = Tracker["constants"]["inires"]
            Tracker["fsc143"] = Tracker["constants"]["inires"]

            checking_flag = 1
            ###  VARIOUS SANITY CHECKES
            if Blockdata["myid"] == Blockdata["main_node"]:
                if Tracker["constants"]["mask3D"] and (
                    not os.path.exists(Tracker["constants"]["mask3D"])
                ):
                    checking_flag = 0
            checking_flag = sp_utilities.bcast_number_to_all(
                checking_flag, Blockdata["main_node"], mpi_comm=mpi.MPI_COMM_WORLD
            )
            if checking_flag == 0:
                sp_global_def.ERROR(
                    "mask3D file does  not exists ", "meridien", 1, Blockdata["myid"]
                )

            if old_div(options.xr, options.ts) < 1.0:
                sp_global_def.ERROR(
                    "Incorrect translational searching settings, search range cannot be smaller than translation step ",
                    "meridien",
                    1,
                    Blockdata["myid"],
                )
            # HOHO
            if (
                2 * (Tracker["currentres"] + Tracker["nxstep"])
                > Tracker["constants"]["nnxo"]
            ):
                sp_global_def.ERROR(
                    "Image size less than what would follow from the initial resolution provided %d  %d  %d"
                    % (
                        Tracker["currentres"],
                        Tracker["nxstep"],
                        2 * (Tracker["currentres"] + Tracker["nxstep"]),
                    ),
                    "sxmeridien",
                    1,
                    Blockdata["myid"],
                )

            if Tracker["constants"]["radius"] < 1:
                Tracker["constants"]["radius"] = (
                    old_div(Tracker["constants"]["nnxo"], 2) - 2
                )
            elif (2 * Tracker["constants"]["radius"] + 2) > Tracker["constants"][
                "nnxo"
            ]:
                sp_global_def.ERROR(
                    "Particle radius set too large!", "sxmeridien", 1, Blockdata["myid"]
                )
            ###<-----end of sanity check <----------------------
            ###<<<----------------------------- parse program

            #  	Fresh start INITIALIZATION
            # ------------------------------------------------------------------------------------
            mainiteration = 0
            Tracker["mainiteration"] = mainiteration

            #  MASTER DIRECTORY
            if Blockdata["myid"] == Blockdata["main_node"]:
                if masterdir == "":
                    timestring = time.strftime("_%d_%b_%Y_%H_%M_%S", time.localtime())
                    masterdir = "master" + timestring
                    li = len(masterdir)
                    os.makedirs(masterdir)
                    keepchecking = 0
                else:
                    if not os.path.exists(masterdir):
                        os.makedirs(masterdir)
                    li = 0
                    keepchecking = 1
                sp_global_def.write_command(masterdir)
            else:
                li = 0
                keepchecking = 1

            li = mpi.mpi_bcast(
                li, 1, mpi.MPI_INT, Blockdata["main_node"], mpi.MPI_COMM_WORLD
            )[0]

            if li > 0:
                masterdir = mpi.mpi_bcast(
                    masterdir,
                    li,
                    mpi.MPI_CHAR,
                    Blockdata["main_node"],
                    mpi.MPI_COMM_WORLD,
                )
                masterdir = string.join(masterdir, "")

            Tracker["constants"]["masterdir"] = masterdir
            initdir = os.path.join(Tracker["constants"]["masterdir"], "main000")
            if Blockdata["myid"] == Blockdata["main_node"]:
                if os.path.exists(initdir):
                    sp_helix_sphire.shutil.rmtree(initdir)
                os.makedirs(initdir)

            if Blockdata["myid"] == Blockdata["main_node"]:
                print_dict(Tracker["constants"], "Permanent settings of meridien")
                print_dict(Blockdata, "MPI settings of meridien")

            # Initialization of orgstack
            Tracker["constants"]["stack"] = orgstack
            if Blockdata["myid"] == Blockdata["main_node"]:
                total_stack = EMAN2_cppwrap.EMUtil.get_image_count(
                    Tracker["constants"]["stack"]
                )
            else:
                total_stack = 0
            total_stack = sp_utilities.bcast_number_to_all(
                total_stack, source_node=Blockdata["main_node"]
            )

            # ------------------------------------------------------------------------------------
            partids = os.path.join(initdir, "indexes_000.txt")
            if Blockdata["myid"] == Blockdata["main_node"]:
                sp_utilities.write_text_file(list(range(total_stack)), partids)
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            #  store params
            partids = [None] * 2
            for procid in range(2):
                partids[procid] = os.path.join(initdir, "chunk_%01d_000.txt" % procid)
            partstack = [None] * 2
            for procid in range(2):
                partstack[procid] = os.path.join(
                    initdir, "params-chunk_%01d_000.txt" % procid
                )

            if Blockdata["myid"] == Blockdata["main_node"]:
                l1, l2 = assign_particles_to_groups(minimum_group_size=10)
                sp_utilities.write_text_file(l1, partids[0])
                sp_utilities.write_text_file(l2, partids[1])
                tp_list = EMAN2_cppwrap.EMUtil.get_all_attributes(
                    Tracker["constants"]["stack"], "xform.projection"
                )
                for i in range(len(tp_list)):
                    dp = tp_list[i].get_params("spider")
                    tp_list[i] = [
                        dp["phi"],
                        dp["theta"],
                        dp["psi"],
                        -dp["tx"],
                        -dp["ty"],
                        0.0,
                        1.0,
                        1.0,
                    ]

                sp_utilities.write_text_row(
                    tp_list, os.path.join(initdir, "params_000.txt")
                )
                sp_utilities.write_text_row([tp_list[i] for i in l1], partstack[0])
                sp_utilities.write_text_row([tp_list[i] for i in l2], partstack[1])

                del tp_list
                del l1, l2
            else:
                Tracker["nima_per_chunk"] = [0, 0]
            Tracker["nima_per_chunk"][0] = sp_utilities.bcast_number_to_all(
                Tracker["nima_per_chunk"][0], Blockdata["main_node"]
            )
            Tracker["nima_per_chunk"][1] = sp_utilities.bcast_number_to_all(
                Tracker["nima_per_chunk"][1], Blockdata["main_node"]
            )
            Tracker["constants"]["number_of_groups"] = sp_utilities.bcast_number_to_all(
                Tracker["constants"]["number_of_groups"], Blockdata["main_node"]
            )

            projdata = [
                [sp_utilities.model_blank(1, 1)],
                [sp_utilities.model_blank(1, 1)],
            ]
            oldparams = [[], []]
            currentparams = [[], []]
            original_data = [None, None]
            # HOHO
            if Tracker["constants"]["inires"] > 0:
                Tracker["nxinit"] = min(
                    2 * Tracker["constants"]["inires"], Tracker["constants"]["nnxo"]
                )
            else:
                Tracker["nxinit"] = Tracker["constants"]["nnxo"]

            rec3d_continuation_nosmearing(original_data, mpi.MPI_COMM_WORLD)

            """Multiline Comment25"""

            if Blockdata["myid"] == Blockdata["main_node"]:
                dump_tracker_to_json(
                    os.path.join(
                        initdir, "Tracker_%03d.json" % Tracker["mainiteration"]
                    ),
                    Tracker,
                )

        else:  # simple restart, at least main000 is completed. Otherwise no need restart

            Blockdata["bckgnoise"] = None
            Blockdata["accumulatepw"] = [[], []]
            projdata = [
                [sp_utilities.model_blank(1, 1)],
                [sp_utilities.model_blank(1, 1)],
            ]
            oldparams = [[], []]
            currentparams = [[], []]
            original_data = [None, None]
            initdir = os.path.join(masterdir, "main000")
            keepchecking = 1
            if Blockdata["myid"] == Blockdata["main_node"]:
                Tracker = load_tracker_from_json(
                    os.path.join(initdir, "Tracker_000.json")
                )
                print_dict(
                    Tracker["constants"],
                    "Permanent settings of the original run recovered from main000",
                )
                stack_prior = Tracker["constants"]["stack_prior"]
                Tracker["constants"]["stack_prior"] = None
                sp_global_def.write_command(Tracker["constants"]["masterdir"])
            else:
                Tracker = None
                stack_prior = None
            Tracker = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"])
            Tracker["constants"]["stack_prior"] = stack_prior
            mainiteration = 0
            Tracker["mainiteration"] = mainiteration
        # ------------------------------------------------------------------------------------
        #  MAIN ITERATION
        keepgoing = 1
        while keepgoing:
            Tracker["previousoutputdir"] = os.path.join(
                Tracker["constants"]["masterdir"], "main%03d" % Tracker["mainiteration"]
            )
            mainiteration += 1
            Tracker["mainiteration"] = mainiteration
            Tracker["directory"] = os.path.join(
                Tracker["constants"]["masterdir"], "main%03d" % Tracker["mainiteration"]
            )
            doit, keepchecking = checkstep(Tracker["directory"], keepchecking)
            if not doit:
                li = True
                doit2, keepchecking2 = checkstep(
                    os.path.join(
                        Tracker["directory"],
                        "Tracker_%03d.json" % Tracker["mainiteration"],
                    ),
                    li,
                )
                if doit2:
                    doit = True
                    if Blockdata["myid"] == Blockdata["main_node"]:
                        sp_helix_sphire.shutil.rmtree(Tracker["directory"])
                    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            if doit:
                if Blockdata["myid"] == Blockdata["main_node"]:
                    Tracker = load_tracker_from_json(
                        os.path.join(
                            Tracker["previousoutputdir"],
                            "Tracker_%03d.json" % (Tracker["mainiteration"] - 1),
                        )
                    )
                    stack_prior = Tracker["constants"]["stack_prior"]
                    Tracker["constants"]["stack_prior"] = None
                    #  It has to be repeated here as Tracker is from previous iteration, I see no other way.
                    Tracker["previousoutputdir"] = os.path.join(
                        Tracker["constants"]["masterdir"],
                        "main%03d" % Tracker["mainiteration"],
                    )
                    Tracker["mainiteration"] = mainiteration
                    Tracker["directory"] = os.path.join(
                        Tracker["constants"]["masterdir"],
                        "main%03d" % Tracker["mainiteration"],
                    )
                else:
                    Tracker = None
                    stack_prior = None
                Tracker = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"])
                Tracker["constants"]["stack_prior"] = stack_prior
                ###
                if restart_mode:
                    update_tracker(sys.argv[1:])
                    update_memory_estimation()
                    restart_mode = False

                # prepare names of input file names, they are in main directory,
                #   log subdirectories contain outputs from specific refinements
                partids = [None] * 2
                for procid in range(2):
                    partids[procid] = os.path.join(
                        Tracker["previousoutputdir"],
                        "chunk_%01d_%03d.txt" % (procid, Tracker["mainiteration"] - 1),
                    )
                partstack = [None] * 2
                for procid in range(2):
                    partstack[procid] = os.path.join(
                        Tracker["previousoutputdir"],
                        "params-chunk_%01d_%03d.txt"
                        % (procid, Tracker["mainiteration"] - 1),
                    )

                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

                if Blockdata["myid"] == Blockdata["main_node"]:
                    fff = sp_utilities.read_text_file(
                        os.path.join(
                            Tracker["previousoutputdir"],
                            "driver_%03d.txt" % (Tracker["mainiteration"] - 1),
                        )
                    )
                    if Tracker["mainiteration"] == 1:
                        anger = 1.0e9
                        shifter = 1.0e9
                    else:
                        [anger, shifter] = sp_utilities.read_text_row(
                            os.path.join(
                                Tracker["previousoutputdir"],
                                "error_thresholds_%03d.txt"
                                % (Tracker["mainiteration"] - 1),
                            )
                        )[0]
                else:
                    fff = []
                    anger = 0.0
                    shifter = 0.0
                fff = sp_utilities.bcast_list_to_all(
                    fff, Blockdata["myid"], source_node=Blockdata["main_node"]
                )
                anger = sp_utilities.bcast_number_to_all(
                    anger, source_node=Blockdata["main_node"]
                )
                shifter = sp_utilities.bcast_number_to_all(
                    shifter, source_node=Blockdata["main_node"]
                )

                keepgoing = AI_continuation(
                    fff, anger, shifter, Blockdata["myid"] == Blockdata["main_node"]
                )

                if keepgoing == 1:  # not converged
                    if Blockdata["myid"] == Blockdata["main_node"]:
                        if Tracker["mainiteration"] > 1:
                            # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
                            line = ""
                            sp_global_def.sxprint(
                                line,
                                "Resolution achieved in ITERATION  #%2d: %3d/%3d pixels, %5.2fA/%5.2fA."
                                % (
                                    Tracker["mainiteration"] - 1,
                                    Tracker["currentres"],
                                    Tracker["fsc143"],
                                    old_div(
                                        Tracker["constants"]["pixel_size"]
                                        * Tracker["constants"]["nnxo"],
                                        float(Tracker["currentres"]),
                                    ),
                                    old_div(
                                        Tracker["constants"]["pixel_size"]
                                        * Tracker["constants"]["nnxo"],
                                        float(Tracker["fsc143"]),
                                    ),
                                ),
                            )
                        sp_global_def.sxprint("\n\n\n\n")
                        # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
                        line = ""
                        sp_global_def.sxprint(
                            line,
                            "ITERATION  #%2d. Current state: %14s, nxinit: %3d, delta: %9.4f, xr: %9.4f, ts: %9.4f"
                            % (
                                Tracker["mainiteration"],
                                Tracker["state"],
                                Tracker["nxinit"],
                                Tracker["delta"],
                                Tracker["xr"],
                                Tracker["ts"],
                            ),
                        )
                    # print("RACING  A ",Blockdata["myid"])
                    li = True
                    doit2, keepchecking2 = checkstep(Tracker["directory"], li)
                    if Blockdata["myid"] == Blockdata["main_node"] and doit2:
                        cmd = "{} {}".format("mkdir -p", Tracker["directory"])
                        junk = sp_utilities.cmdexecute(cmd)
                        cmd = "{} {}".format(
                            "mkdir -p",
                            os.path.join(Tracker["directory"], "oldparamstructure"),
                        )
                        junk = sp_utilities.cmdexecute(cmd)
                    if not doit2:
                        sp_global_def.ERROR(
                            "There was a gap in main directories, program cannot proceed",
                            "sxmeridien",
                            1,
                            Blockdata["myid"],
                        )
                    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

                    #  READ DATA AND COMPUTE SIGMA2   ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
                    refinement_one_iteration(
                        partids,
                        partstack,
                        original_data,
                        oldparams,
                        projdata,
                        general_mode=False,
                        continuation_mode=True,
                    )
                    # 	sxprint("  MOVING  ON --------------------------------------------------------------------")
                else:  # converged, do final
                    if Blockdata["subgroup_myid"] > -1:
                        mpi.mpi_comm_free(Blockdata["subgroup_comm"])

                    if Blockdata["myid"] == Blockdata["main_node"]:
                        sp_global_def.sxprint(
                            line,
                            "Resolution achieved in ITERATION  #%2d: %3d/%3d pixels, %5.2fA/%5.2fA."
                            % (
                                Tracker["mainiteration"] - 1,
                                Tracker["currentres"],
                                Tracker["fsc143"],
                                old_div(
                                    Tracker["constants"]["pixel_size"]
                                    * Tracker["constants"]["nnxo"],
                                    float(Tracker["currentres"]),
                                ),
                                old_div(
                                    Tracker["constants"]["pixel_size"]
                                    * Tracker["constants"]["nnxo"],
                                    float(Tracker["fsc143"]),
                                ),
                            ),
                        )
                        sp_global_def.sxprint("\n\n\n\n")

                        sp_global_def.sxprint(
                            "Iteration that contains the best resolution is %d"
                            % Tracker["constants"]["best"]
                        )
                        if Tracker["constants"]["best"] == 2:
                            sp_global_def.sxprint(
                                "No resolution improvement in refinement "
                            )

                        dump_tracker_to_json(
                            os.path.join(masterdir, "Tracker_final.json"), Tracker
                        )
                    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

                    newparamstructure = [[], []]
                    projdata = [
                        [sp_utilities.model_blank(1, 1)],
                        [sp_utilities.model_blank(1, 1)],
                    ]
                    original_data = [None, None]
                    oldparams = [[], []]
                    Blockdata["accumulatepw"] = [None, None]
                    # if Tracker["constants"]["memory_per_node"] <0.0: Tracker["constants"]["memory_per_node"] = 2.0*Blockdata["no_of_processes_per_group"]
                    recons3d_final(
                        Tracker["constants"]["masterdir"],
                        Tracker["constants"]["best"],
                        Tracker["constants"]["memory_per_node"],
                    )

            #  End of if doit
        #   end of while
        mpi.mpi_finalize()
        exit()

    elif do_final_mode:  # DO FINAL
        parser.add_option(
            "--memory_per_node",
            type="float",
            default=-1.0,
            help="User provided information about memory per node (NOT per CPU) [in GB] (default 2GB*(number of CPUs per node))",
        )
        (options, args) = parser.parse_args(sys.argv[1:])
        # global Tracker, Blockdata
        # print( "  args  ",args)
        checking_flag = 1
        orgstack = None
        if len(args) == 3:
            sp_global_def.ERROR(
                "do_final option requires only one or two arguments ",
                "meridien",
                1,
                Blockdata["myid"],
            )

        elif len(args) == 2:  # option for signal subtraction
            masterdir = args[1]
            orgstack = args[0]
            if Blockdata["myid"] == Blockdata["main_node"]:
                if not os.path.exists(masterdir):
                    checking_flag = 0
            checking_flag = sp_utilities.bcast_number_to_all(
                checking_flag, source_node=Blockdata["main_node"]
            )
            if checking_flag == 0:
                sp_global_def.ERROR(
                    "do_final: refinement directory for final reconstruction does not exist ",
                    "meridien",
                    1,
                    Blockdata["myid"],
                )

        elif len(args) == 1:
            masterdir = args[0]
            if Blockdata["myid"] == Blockdata["main_node"]:
                if not os.path.exists(masterdir):
                    checking_flag = 0
            checking_flag = sp_utilities.bcast_number_to_all(
                checking_flag, source_node=Blockdata["main_node"]
            )
            if checking_flag == 0:
                sp_global_def.ERROR(
                    "do_final: refinement directory for final reconstruction does not exist ",
                    "meridien",
                    1,
                    Blockdata["myid"],
                )

        else:
            sp_global_def.sxprint("usage: " + usage)
            sp_global_def.sxprint(
                "Please run '" + progname + " -h' for detailed options"
            )
            return 1

        if options.do_final < 0:
            sp_global_def.ERROR(
                "Incorrect iteration number in do_final  %d" % options.do_final,
                "meridien",
                1,
                Blockdata["myid"],
            )
        # print(  orgstack,masterdir,volinit )
        # ------------------------------------------------------------------------------------
        # Initialize MPI related variables

        ###print("  MPIINFO  ",Blockdata)
        ###  MPI SANITY CHECKES
        if not balanced_processor_load_on_nodes:
            sp_global_def.ERROR(
                "Nodes do not have the same number of CPUs, please check configuration of the cluster.",
                "meridien",
                1,
                Blockdata["myid"],
            )
        # if( Blockdata["no_of_groups"] < 2 ):  ERROR("To run, program requires cluster with at least two nodes.","meridien",1,Blockdata["myid"])
        ###
        if Blockdata["myid"] == Blockdata["main_node"]:
            line = ""
            for a in sys.argv:
                line += a + "  "
            sp_global_def.sxprint(" shell line command ")
            sp_global_def.sxprint(line)

        if Blockdata["myid"] == Blockdata["main_node"]:
            sp_global_def.write_command(masterdir)

        # ------------------------------------------------------------------------------------
        #  INPUT PARAMETERS
        sp_global_def.BATCH = True
        sp_global_def.MPI = True

        ###  VARIOUS SANITY CHECKES <-----------------------
        if options.memory_per_node < 0.0:
            options.memory_per_node = 2.0 * Blockdata["no_of_processes_per_group"]

        Blockdata["accumulatepw"] = [[], []]
        recons3d_final(masterdir, options.do_final, options.memory_per_node, orgstack)
        mpi.mpi_finalize()
        exit()
    else:
        sp_global_def.ERROR("Incorrect input options", "meridien", 1, Blockdata["myid"])


def main():
    run()

if __name__ == "__main__":
    main()
