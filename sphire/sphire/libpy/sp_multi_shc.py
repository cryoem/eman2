# For some reason these programs get stuck on MPI if I change the order of programs in the file.  Strange, PAP.
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

import EMAN2_cppwrap
import copy
from copy import deepcopy
import math
import mpi
import os
import random
from . import sp_alignment
from . import sp_applications
from . import sp_filter
from . import sp_fundamentals
from . import sp_global_def
from . import sp_logger
from . import sp_morphology
from . import sp_pixel_error
from . import sp_projection
from . import sp_reconstruction
from . import sp_statistics
from . import sp_utilities
from ..libpy.sp_utilities import getvec, getfvec, getang3, lacos, angles_to_normals
import time


def orient_params(params, refparams, indexes=None, symmetry_class=None):
    #
    #  The assumption here is that the angles are within the unique region
    #  Since they come from projection refinement, why would it be otherwise?
    #  Any problems would be due to how even_angles generates reference angles
    #
    #  The function returns rotation object and properly rotated/mirrored params

    n = len(params)
    if indexes == None:
        tindexes = list(range(n))
    else:
        tindexes = indexes

    if symmetry_class.sym[0] == "c" and symmetry_class.nsym > 1:
        divic = old_div(360.0, symmetry_class.nsym)
        phi = sp_pixel_error.angle_diff_sym(
            [params[j] for j in tindexes],
            [refparams[j] for j in tindexes],
            symmetry_class.nsym,
        )
        out = copy.deepcopy(params)
        for j in range(n):
            out[j][0] = (out[j][0] + phi) % divic
        # mirror checking
        psi_diff = sp_pixel_error.angle_diff(
            [out[j][2] for j in tindexes], [refparams[j][2] for j in tindexes]
        )
        if abs(psi_diff - 180.0) < 90.0:
            for j in range(n):
                # apply mirror
                out[j][2] = (out[j][2] + 180.0) % 360.0
    elif symmetry_class.sym[0] == "c":
        t1, t2, t3 = sp_utilities.rotation_between_anglesets(
            [params[j] for j in tindexes], [refparams[j] for j in tindexes]
        )
        out = sp_fundamentals.rotate_params(
            [params[i][:3] for i in range(n)], [-t3, -t2, -t1]
        )
        out = [out[i] + params[i][3:] for i in range(n)]  # reattach shifts
        # mirror checking
        psi_diff = sp_pixel_error.angle_diff(
            [out[j][2] for j in tindexes], [refparams[j][2] for j in tindexes]
        )
        if abs(psi_diff - 180.0) < 90.0:
            for j in range(n):
                # apply mirror
                out[j][2] = (out[j][2] + 180.0) % 360.0
    """Multiline Comment0"""
    return out


def find_common_subset(
    projs, target_threshold=2.0, minimal_subset_size=3, symmetry_class=None
):
    #  projs - [reference set of angles, set of angles1, set of angles2, ... ]
    #  the function is written for multiple sets of angles
    #  The transformation is found for a subset, but applied to the full set

    n = len(projs[0])
    sc = len(projs)

    subset = list(range(n))

    minimal_subset_size = min(minimal_subset_size, n)

    avg_diff_per_image = [0.0] * n

    if (
        symmetry_class.sym[:3] == "oct"
        or symmetry_class.sym[:3] == "tet"
        or symmetry_class.sym[:4] == "icos"
    ):
        # In these cases there is no rotation so distances between angular directions cannot change
        outp = copy.deepcopy(projs)
        for i in subset:
            qt = 0.0
            for k in range(sc - 1):
                for l in range(k + 1, sc):
                    neisym = symmetry_class.symmetry_neighbors([outp[l][i][:3]])
                    dmin = 180.0
                    for q in neisym:
                        dmin = min(dmin, getang3(q, outp[k][i]))
                    qt += dmin
            avg_diff_per_image[i] = old_div(old_div(old_div(qt, sc), (sc - 1)), 2.0)

    #  Start from the entire set and the slowly decrease it by rejecting worst data one by one.
    #  It will stop either when the subset reaches the minimum subset size
    #   or if there are no more angles with errors above target threshold
    #
    while True:
        #  extract images in common subset
        if symmetry_class.sym[0] == "c":  #  c including cn symmetry
            divic = old_div(360.0, symmetry_class.nsym)
            for i in subset:
                avg_diff_per_image[i] = 0.0
            outp = [copy.deepcopy(projs[0])]
            for i in range(sc - 1, -1, -1):
                tv = sp_utilities.angles_to_normals(projs[i])

                for j in range(i + 1, sc):
                    out = orient_params(
                        projs[j], projs[i], subset, symmetry_class=symmetry_class
                    )
                    if symmetry_class.nsym > 1:
                        # for k in xrange(n):
                        for k in subset:
                            mind = 1.0e23
                            for l in range(symmetry_class.nsym):
                                u1, u2, u3 = sp_utilities.getfvec(
                                    out[k][0] + l * divic, out[k][1]
                                )
                                qt = sp_utilities.lacos(
                                    tv[k][0] * u1 + tv[k][1] * u2 + tv[k][2] * u3
                                )
                                mind = min(mind, qt)
                            avg_diff_per_image[k] += mind
                            # print  "avg_diff_per_image  %3d  %8.2f="%(k,avg_diff_per_image[k]),\
                            # "  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f"%( projs[i][k][0],projs[i][k][1],projs[i][k][2],out[k][0],out[k][1],out[k][2])
                    else:
                        # for k in xrange(n):
                        for k in subset:
                            u1, u2, u3 = sp_utilities.getfvec(out[k][0], out[k][1])
                            avg_diff_per_image[k] += sp_utilities.lacos(
                                tv[k][0] * u1 + tv[k][1] * u2 + tv[k][2] * u3
                            )
                            # print  "avg_diff_per_image  %3d  %8.2f %8.6f="%(k,avg_diff_per_image[k], tv[k][0]*u1+tv[k][1]*u2+tv[k][2]*u3),\
                            # "  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f"%( projs[i][k][0],projs[i][k][1],projs[i][k][2],out[k][0],out[k][1],out[k][2])
                    if i == 0:
                        outp.append(copy.deepcopy(out))

            # for k in range(n):
            for k in subset:
                avg_diff_per_image[k] = old_div(
                    avg_diff_per_image[k], old_div(sc * (sc - 1), 2.0)
                )
        elif symmetry_class.sym[0] == "d":
            outp = copy.deepcopy(projs)
            mirror_and_reduce_dsym(outp, subset, symmetry_class)

            for i in subset:
                qt = 0.0
                for k in range(sc - 1):
                    for l in range(k + 1, sc):
                        neisym = symmetry_class.symmetry_neighbors([outp[l][i][:3]])
                        dmin = 180.0
                        for q in neisym:
                            dmin = min(dmin, getang3(q, outp[k][i]))
                        qt += dmin
                avg_diff_per_image[i] = old_div(old_div(old_div(qt, sc), (sc - 1)), 2.0)
                # k = subset[i]
                # lml = 1
                # print  "avg_diff_per_image  %3d  %8.2f="%(k,avg_diff_per_image[k]),\
                # "  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f"%( projs2[0][i][0],projs2[0][i][1],projs2[0][i][2],projs2[lml][i][0],projs2[lml][i][1],projs2[lml][i][2])

        else:  # o, t, i
            outp = copy.deepcopy(projs)
            for l in range(1, sc):
                psi_diff = sp_pixel_error.angle_diff(
                    [outp[l][j][2] for j in subset], [outp[0][j][2] for j in subset]
                )
                #  adjust psi if necessary
                if abs(psi_diff - 180.0) < 90.0:
                    for j in range(n):
                        # apply mirror
                        outp[l][j][2] = (outp[l][j][2] + 180.0) % 360.0

        if len(subset) == minimal_subset_size:
            break

        #  Remove element whose avg_diff_per_image is larger than max_error, if none, break
        max_error = -1.0
        the_worst_proj = -1
        for i in subset:
            if avg_diff_per_image[i] > target_threshold:
                if avg_diff_per_image[i] > max_error:
                    max_error = avg_diff_per_image[i]
                    the_worst_proj = i
        if the_worst_proj > -1:
            subset.remove(the_worst_proj)
        else:
            break
        # print  "the_worst_proj",the_worst_proj
    #  End of pruning loop

    return subset, avg_diff_per_image, outp


"""Multiline Comment1"""


# parameters: list of (all) projections | reference volume | ...
#  Genetic programming version
#  The data structure:
#  [[L2, [parameters row-wise]], [], []...number_of_runs ]
#  It is kept on main proc


def ali3d_multishc(
    stack,
    ref_vol,
    ali3d_options,
    symmetry_class,
    mpi_comm=None,
    log=None,
    number_of_runs=2,
):

    ir = ali3d_options.ir
    rs = ali3d_options.rs
    ou = ali3d_options.ou
    xr = ali3d_options.xr
    yr = ali3d_options.yr
    ts = ali3d_options.ts
    an = ali3d_options.an
    delta = ali3d_options.delta
    doga = ali3d_options.doga
    center = ali3d_options.center
    CTF = ali3d_options.CTF
    ref_a = ali3d_options.ref_a
    L2threshold = ali3d_options.L2threshold

    # Optionally restrict out-of-plane angle
    if hasattr(ali3d_options, "theta1"):
        theta1 = ali3d_options.theta1
    else:
        theta1 = -1.0

    if hasattr(ali3d_options, "theta2"):
        theta2 = ali3d_options.theta2
    else:
        theta2 = -1.0

    if hasattr(ali3d_options, "method"):
        method = ali3d_options.method
    else:
        method = "S"

    if mpi_comm == None:
        mpi_comm = mpi.MPI_COMM_WORLD

    if log == None:
        log = sp_logger.Logger()

    number_of_proc = mpi.mpi_comm_size(mpi_comm)
    myid = mpi.mpi_comm_rank(mpi_comm)
    main_node = 0

    if myid == main_node:
        log.add("Start VIPER1")

    if number_of_proc < number_of_runs:
        sp_global_def.ERROR("number_of_proc < number_of_runs", "VIPER1", 1, myid)

    # if an != "-1":
    # 	ERROR("Option an not used","VIPER1",1,myid)

    # mpi_subcomm = mpi_comm_split(mpi_comm, myid % number_of_runs, myid / number_of_runs)
    # mpi_subrank = mpi_comm_rank(mpi_subcomm)
    # mpi_subsize = mpi_comm_size(mpi_subcomm)
    # mpi_subroots = range(number_of_runs)

    mpi_subcomm = sp_utilities.wrap_mpi_split(mpi_comm, number_of_runs)
    mpi_subrank = mpi.mpi_comm_rank(mpi_subcomm)
    mpi_subsize = mpi.mpi_comm_size(mpi_subcomm)
    # do not make any assumptions about the subroots, collect the rank_id as they are already assigned
    if mpi_subrank == 0:
        mpi_subroots = sp_utilities.wrap_mpi_gatherv([myid], 0, mpi_comm)
    else:
        mpi_subroots = sp_utilities.wrap_mpi_gatherv([], 0, mpi_comm)
    mpi_subroots = sp_utilities.wrap_mpi_bcast(mpi_subroots, main_node, mpi_comm)

    xrng = sp_utilities.get_input_from_string(xr)
    if yr == "-1":
        yrng = xrng
    else:
        yrng = sp_utilities.get_input_from_string(yr)
    step = sp_utilities.get_input_from_string(ts)
    delta = sp_utilities.get_input_from_string(delta)
    lstp = min(len(xrng), len(yrng), len(step), len(delta))
    """Multiline Comment2"""
    first_ring = int(ir)
    rstep = int(rs)
    last_ring = int(ou)
    max_iter = int(ali3d_options.maxit1)
    center = int(center)

    vol = ref_vol
    nx = vol.get_xsize()
    if last_ring < 0:
        last_ring = int(old_div(nx, 2)) - 2

    cnx = old_div(nx, 2) + 1
    cny = cnx
    numr = sp_alignment.Numrinit(first_ring, last_ring, rstep, "F")
    mask2D = sp_utilities.model_circle(last_ring, nx, nx) - sp_utilities.model_circle(
        first_ring, nx, nx
    )

    if myid == main_node:
        list_of_particles = list(range(len(stack)))
        total_nima = len(list_of_particles)
    else:
        list_of_particles = None
        total_nima = None
    total_nima = sp_utilities.wrap_mpi_bcast(total_nima, main_node, mpi_comm)
    list_of_particles = sp_utilities.wrap_mpi_bcast(
        list_of_particles, main_node, mpi_comm
    )
    nima = len(list_of_particles)

    image_start, image_end = sp_applications.MPI_start_end(
        total_nima, mpi_subsize, mpi_subrank
    )
    # print "  image_start, image_end  ", myid,image_start, image_end

    data = [stack[im] for im in list_of_particles]
    for im in range(nima):
        data[im].set_attr("ID", list_of_particles[im])
        ctf_applied = data[im].get_attr_default("ctf_applied", 0)
        if CTF and ctf_applied == 0:
            ctf_params = data[im].get_attr("ctf")
            st = EMAN2_cppwrap.Util.infomask(data[im], mask2D, False)
            data[im] -= st[0]
            data[im] = sp_filter.filt_ctf(data[im], ctf_params)
            data[im].set_attr("ctf_applied", 1)

    cs = [0.0] * 3
    total_iter = 0

    # initialize GA data structure [ [L2, [params]], ...]
    if myid == main_node:
        GA = [
            [0.0, [[0.0, 0.0, 0.0, 0.0, 0.0] for j in range(total_nima)]]
            for i in range(number_of_runs)
        ]
        noimprovement = 0
        firstcheck = True

    orient_and_shuffle = False
    # do the projection matching
    for N_step in range(
        lstp
    ):  # At this point there is just one value here, it cannot loop.
        if myid == 0:
            afterGAcounter = 0
        terminate = False
        Iter = 0
        while not terminate:

            Iter += 1
            total_iter += 1

            mpi.mpi_barrier(mpi_comm)
            if myid == main_node:
                log.add("ITERATION #%3d" % (total_iter))
                start_time = time.time()

            # =========================================================================
            # build references
            volft, kb = sp_projection.prep_vol(vol)
            #  We generate mirrored versions as well MAJOR CHANGE PAP 04/20/2017
            reference_angles = symmetry_class.even_angles(
                delta[N_step], theta1=theta1, theta2=theta2, method=method
            )
            refrings = sp_alignment.prepare_refrings(
                volft, kb, nx, -1.0, reference_angles, "", numr, MPI=mpi_subcomm
            )
            del volft, kb
            # =========================================================================

            mpi.mpi_barrier(mpi_comm)
            if myid == main_node:
                log.add("Time to prepare rings: %f" % (time.time() - start_time))
                start_time = time.time()

            # =========================================================================
            #  It has to be here
            if orient_and_shuffle:
                # adjust params to references, calculate psi, calculate previousmax
                vecs = [
                    [
                        refrings[lr].get_attr("n1"),
                        refrings[lr].get_attr("n2"),
                        refrings[lr].get_attr("n3"),
                    ]
                    for lr in range(len(refrings))
                ]
                for im in range(nima):
                    #  For testing purposes make sure that directions did not change
                    t1, t2, t3, t4, t5 = sp_utilities.get_params_proj(data[im])
                    iqa = sp_utilities.nearest_fang(
                        vecs, t1, t2
                    )  # Here I could use more sophisticated distance for symmetries
                    # if myid == 0 :
                    # print "  XXXX  ",myid,total_iter,im,iqa,t1,t2, reference_angles[iqa]
                    # if total_iter>3:
                    # 	'''
                    # 	from mpi import mpi_finalize
                    # 	mpi_finalize()
                    # 	from sys import exit
                    # 	exit()
                    # 	'''

                    # peak, temp = proj_ali_incore_local(data[im], [refrings[iqa]], [reference_angles[iqa]], numr, 0., 0., 1., 180.0 , sym="c1")
                    cimage = EMAN2_cppwrap.Util.Polar2Dm(data[im], cnx, cny, numr, "F")
                    EMAN2_cppwrap.Util.Normalize_ring(cimage, numr, 0)
                    EMAN2_cppwrap.Util.Frngs(cimage, numr)
                    retvals = EMAN2_cppwrap.Util.Crosrng_e(
                        refrings[iqa], cimage, numr, 0, 0.0
                    )
                    data[im].set_attr("previousmax", retvals["qn"])
                    # u1,u2,t3,t4,t5 = get_params_proj( data[im])
                    # if(abs(u1-t1)>1.0e-4 or  abs(u2-t2)>1.0e-4 ):  print "  PROBLEM IN  proj_ali_incore_local"
                    # print  "peak1 ",myid,im,data[im].get_attr("ID"), retvals["qn"],[t1,t2,t3,t4,t5]
                if myid == main_node:
                    log.add(
                        "Time to calculate first psi+previousmax: %f\n"
                        % (time.time() - start_time)
                    )
                    start_time = time.time()
                del vecs
            # =========================================================================

            mpi.mpi_barrier(mpi_comm)

            if myid == main_node:
                start_time = time.time()
            # =========================================================================
            # alignment
            if mpi_subrank == 0:
                pixer = []
                number_of_checked_refs = 0
                proj_ids_to_process = list(range(nima))
                random.shuffle(proj_ids_to_process)
            while True:
                # --------- broadcast projections ids
                if mpi_subrank == 0:
                    if len(proj_ids_to_process) >= mpi_subsize:
                        proj_ids = proj_ids_to_process[:(mpi_subsize)]
                    else:
                        proj_ids = proj_ids_to_process[:]
                else:
                    proj_ids = None
                proj_ids = sp_utilities.wrap_mpi_bcast(proj_ids, 0, mpi_subcomm)
                if len(proj_ids) == 0:
                    break
                if mpi_subrank < len(proj_ids):
                    # -------- alignment
                    im = proj_ids[mpi_subrank]
                    peak, pixel_error, checked_refs, iref = sp_alignment.shc(
                        data[im],
                        refrings,
                        [[1.0, 1.0]],
                        numr,
                        xrng[N_step],
                        yrng[N_step],
                        step[N_step],
                        sym="nomirror",
                    )
                    # -------- gather results to root
                    vector_assigned_refs = sp_utilities.wrap_mpi_gatherv(
                        [iref], 0, mpi_subcomm
                    )
                    vector_previousmax = sp_utilities.wrap_mpi_gatherv(
                        [data[im].get_attr("previousmax")], 0, mpi_subcomm
                    )
                    vector_xformprojs = sp_utilities.wrap_mpi_gatherv(
                        [data[im].get_attr("xform.projection")], 0, mpi_subcomm
                    )
                    vector_pixel_error = sp_utilities.wrap_mpi_gatherv(
                        [pixel_error], 0, mpi_subcomm
                    )
                    vector_checked_ref = sp_utilities.wrap_mpi_gatherv(
                        [checked_refs], 0, mpi_subcomm
                    )
                else:
                    # -------- no projection assigned, send to root empty lists
                    vector_assigned_refs = sp_utilities.wrap_mpi_gatherv(
                        [], 0, mpi_subcomm
                    )
                    vector_previousmax = sp_utilities.wrap_mpi_gatherv(
                        [], 0, mpi_subcomm
                    )
                    vector_xformprojs = sp_utilities.wrap_mpi_gatherv(
                        [], 0, mpi_subcomm
                    )
                    vector_pixel_error = sp_utilities.wrap_mpi_gatherv(
                        [], 0, mpi_subcomm
                    )
                    vector_checked_ref = sp_utilities.wrap_mpi_gatherv(
                        [], 0, mpi_subcomm
                    )
                # -------- merge results
                if mpi_subrank == 0:
                    used_refs = set()
                    for i in range(len(vector_assigned_refs)):
                        ir = vector_assigned_refs[i]
                        if ir in used_refs:
                            # reference is already used - cancel all changes
                            vector_previousmax[i] = data[proj_ids[i]].get_attr(
                                "previousmax"
                            )
                            vector_xformprojs[i] = data[proj_ids[i]].get_attr(
                                "xform.projection"
                            )
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
                used_refs = sp_utilities.wrap_mpi_bcast(used_refs, 0, mpi_subcomm)
                vector_previousmax = sp_utilities.wrap_mpi_bcast(
                    vector_previousmax, 0, mpi_subcomm
                )
                vector_xformprojs = sp_utilities.wrap_mpi_bcast(
                    vector_xformprojs, 0, mpi_subcomm
                )
                # ------- delete used references
                for ir in used_refs:
                    del refrings[ir]
                # ------- set projections parameters
                for i in range(len(vector_previousmax)):
                    data[proj_ids[i]].set_attr("previousmax", vector_previousmax[i])
                    data[proj_ids[i]].set_attr("xform.projection", vector_xformprojs[i])
            # =========================================================================
            mpi.mpi_barrier(mpi_comm)
            if myid == main_node:
                log.add("Time of alignment = %f\n" % (time.time() - start_time))

            storevol = False

            # =========================================================================
            # output pixel errors, check stop criterion
            if mpi_subrank == 0:
                all_pixer = sp_utilities.wrap_mpi_gatherv(pixer, 0, mpi_comm)
                total_checked_refs = sp_utilities.wrap_mpi_gatherv(
                    [number_of_checked_refs], main_node, mpi_comm
                )
            else:
                all_pixer = sp_utilities.wrap_mpi_gatherv([], 0, mpi_comm)
                total_checked_refs = sp_utilities.wrap_mpi_gatherv(
                    [], main_node, mpi_comm
                )
            if myid == main_node:
                total_checked_refs = sum(total_checked_refs)
                lhist = 20
                region, histo = sp_statistics.hist_list(all_pixer, lhist)
                log.add("= Pixel error        Number of images in all runs")
                for lhx in range(lhist):
                    msg = " %10.3f                  %7d" % (region[lhx], histo[lhx])
                    log.add(msg)
                temp = 0
                for i in all_pixer:
                    if i < 1.0:
                        temp += 1
                percent_of_pixerr_below_one = old_div(
                    (temp * 1.0), (total_nima * number_of_runs)
                )
                orient_and_shuffle = (percent_of_pixerr_below_one > doga) and (
                    afterGAcounter < 0
                )  #  TODO - parameter ?
                afterGAcounter -= 1
                ## if total_iter%3 == 0:  orient_and_shuffle = True
                ## else:   orient_and_shuffle = False
                # terminate          = ( percent_of_pixerr_below_one > 0.9 )  #  TODO - parameter ?
                log.add("=================================================")
                log.add(
                    "Percent of positions with pixel error below 1.0 = ",
                    (int(percent_of_pixerr_below_one * 100)),
                    "%",
                    "   Mutations: ",
                    orient_and_shuffle,
                )
            orient_and_shuffle = sp_utilities.wrap_mpi_bcast(
                orient_and_shuffle, 0, mpi_comm
            )
            # =========================================================================

            # =========================================================================
            # centering, for d unnecessary, for cn, n>1 only z can move
            if center == -1 and symmetry_class.sym[0] == "c":
                cs[0], cs[1], cs[2], dummy, dummy = sp_utilities.estimate_3D_center_MPI(
                    data[image_start:image_end],
                    total_nima,
                    mpi_subrank,
                    mpi_subsize,
                    0,
                    mpi_comm=mpi_subcomm,
                )  # estimate_3D_center_MPI(data, number_of_runs*total_nima, myid, number_of_proc, main_node, mpi_comm=mpi_comm)
                if myid == main_node:
                    msg = (
                        " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"
                        % (cs[0], cs[1], cs[2])
                    )
                    log.add(msg)
                if symmetry_class.nsym > 1:
                    cs[0] = cs[1] = 0.0
                    if myid == main_node:
                        log.add(
                            "For symmetry group cn (n>1), we only center the volume in z-direction\n"
                        )
                cs = mpi.mpi_bcast(cs, 3, mpi.MPI_FLOAT, main_node, mpi_subcomm)
                cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
                sp_utilities.rotate_3D_shift(data, cs)
            # ===================== CORRECT PARAMETERS ON DATA =======================

            mpi.mpi_barrier(mpi_comm)
            if myid == main_node:
                start_time = time.time()
            """Multiline Comment3"""

            # =========================================================================
            if orient_and_shuffle:  #  DO orient
                params = []
                for im in data:
                    phi, theta, psi, sx, sy = sp_utilities.get_params_proj(im)
                    params.append([phi, theta, psi, sx, sy])
                # if myid == 2:  print  " initial params before orient  ",myid,[get_params_proj(data[i]) for i in xrange(4)]

                # ------ orientation - begin
                #  Send solution from the main process of the first group to all processes in all groups
                params_0 = sp_utilities.wrap_mpi_bcast(
                    params, mpi_subroots[0], mpi_comm
                )
                if (mpi_subrank == 0) and (myid != 0):
                    #  This is done on the main node of each group (master node for MPI_COMM_WORLD skips it)
                    #  Minimal length of the subset is set to 1/3 of the number of parameters
                    #  Error threshold is set somewhat arbitrarily to 1.5 angular step of reference projections
                    #  params gets overwritten by rotated parameters,  subset is a list of indexes common
                    subset, avg_diff_per_image, params = find_common_subset(
                        [params_0, params],
                        delta[N_step] * 1.5,
                        old_div(len(params), 3),
                        symmetry_class,
                    )
                    params = params[1]
                    #   neither is needed
                    del subset, avg_diff_per_image
                    # if myid == 2:  print  " params before orient  ",myid,params[:4],params[-4:]
                    # write_text_row(params_0,"bparamszero%04d%04d.txt"%(myid,total_iter))
                    # write_text_row(params,"bparams%04d%04d.txt"%(myid,total_iter))
                params = sp_utilities.wrap_mpi_bcast(params, 0, mpi_subcomm)
                # if myid == 2:  print  " params after wrap_mpi_bcast  ",myid,params[:4],params[-4:]
                # ------ orientation - end

                # =========================================================================
                # volume reconstruction
                mpi.mpi_barrier(mpi_comm)
                if myid == main_node:
                    start_time = time.time()

                # temp = [None]*nima
                # for i in xrange(nima): temp[i] = data[i].get_attr("xform.projection")
                for i in range(nima):
                    sp_utilities.set_params_proj(data[i], params[i])
                vol = do_volume(
                    data[image_start:image_end], ali3d_options, 0, mpi_subcomm
                )
                # for i in xrange(nima): data[i].set_attr("xform.projection",temp[i])
                # del temp

                if mpi_subrank == 0:
                    L2 = vol.cmp(
                        "dot",
                        vol,
                        dict(
                            negative=0,
                            mask=sp_utilities.model_circle(last_ring, nx, nx, nx),
                        ),
                    )
                    # if myid == 2:  print  " Right after reconstruction L2", myid, L2,[get_params_proj(data[i]) for i in xrange(4)]
                    # print  " Right after reconstruction of oriented parameters L2", myid, total_iter,L2
                    # vol.write_image("recvolf%04d%04d.hdf"%(myid,total_iter))
                # log
                if myid == main_node:
                    log.add(
                        "3D reconstruction time = %f\n" % (time.time() - start_time)
                    )
                    start_time = time.time()

                # ------ gather parameters to root
                if myid == 0:
                    all_L2s = []
                    all_params = []
                    for sr in mpi_subroots:
                        if sr == myid:
                            all_L2s.append(L2)
                            all_params.append(copy.deepcopy(params))
                        else:
                            all_L2s.append(sp_utilities.wrap_mpi_recv(sr, mpi_comm))
                            all_params.append(sp_utilities.wrap_mpi_recv(sr, mpi_comm))
                else:
                    if mpi_subrank == 0:
                        sp_utilities.wrap_mpi_send(L2, 0, mpi_comm)
                        sp_utilities.wrap_mpi_send(params, 0, mpi_comm)

                # ---------------------------------

                #  Add params to GA, sort, check termination and if not terminate do crossover and send back
                if myid == 0:
                    #  after GA move do 3 iterations to give the program a chance to improve mutated structures.
                    # all_params = shuffle_configurations(all_params)

                    for i in range(number_of_runs):
                        log.add(
                            "L2 incoming norm for volume %3d  = %f" % (i, all_L2s[i])
                        )

                    for i in range(number_of_runs):
                        GA.append([all_L2s[i], copy.deepcopy(all_params[i])])
                    #  check whether this move will improve anything
                    all_L2s.sort(reverse=True)
                    # print " sorted terminate  ",all_L2s
                    # for i in xrange(number_of_runs): print GA[i][0]

                    noreseeding = True

                    if all_L2s[0] < GA[number_of_runs - 1][0]:
                        if firstcheck:
                            sp_global_def.sxprint("  SHOULD NOT BE HERE")
                        noimprovement += 1
                        if noimprovement == 2:
                            terminate = True
                        GA = GA[:number_of_runs]
                        log.add(
                            "VIPER1 could not find better solutions, it will terminate now."
                        )
                    else:
                        noimprovement = 0
                        GA.sort(reverse=True)
                        GA = GA[:number_of_runs]
                        if symmetry_class.sym[0] == "d" and GA[0][0] > 0.0:
                            for i in range(1, len(GA)):
                                mirror_and_reduce_dsym(
                                    [GA[0][1], GA[i][1]],
                                    list(range(len(GA[0][1]))),
                                    symmetry_class,
                                )

                        #  ---  Stopping criterion
                        q1, q2, q3, q4 = sp_statistics.table_stat(
                            [GA[i][0] for i in range(number_of_runs)]
                        )
                        # Terminate if variation of L2 norms less than (L2threshold*100)% of their average
                        crit = old_div(math.sqrt(max(q2, 0.0)), q1)
                        for i in range(number_of_runs):
                            log.add("L2 norm for volume %3d  = %f" % (i, GA[i][0]))
                        log.add("L2 norm std dev %f\n" % crit)
                        crit = crit < L2threshold
                        if (Iter < max_iter) and (firstcheck and crit):
                            noreseeding = False
                            terminate = False
                            log.add("Insufficient initial variability, reseeding!\n")
                            all_params = [
                                [
                                    [
                                        random.random() * 360.0,
                                        random.random() * 180.0,
                                        random.random() * 360.0,
                                        0.0,
                                        0.0,
                                    ]
                                    for j in range(total_nima)
                                ]
                                for i in range(number_of_runs)
                            ]
                        else:
                            terminate = Iter > max_iter or crit

                        firstcheck = False

                    if not terminate and noimprovement == 0 and noreseeding:
                        afterGAcounter = 3
                        #  Now do the crossover
                        all_params = []

                        # select random pairs of solutions
                        ipl = list(range(number_of_runs))
                        random.shuffle(ipl)
                        for ip in range(
                            0, 2 * (old_div(len(ipl), 2)) + len(ipl) % 2, 2
                        ):
                            #  random reference projection:
                            itmp = random.randint(0, total_nima - 1)
                            # print  "  nearest_many_full_k_projangles  ",total_nima,itmp,ipl,ip,len(GA[ipl[ip]][1]),GA[ipl[ip]][1][itmp],GA[ipl[ip]][1]
                            keepset = sp_utilities.nearest_many_full_k_projangles(
                                sp_utilities.angles_to_normals(GA[ipl[ip]][1]),
                                [GA[ipl[ip]][1][itmp]],
                                howmany=old_div(total_nima, 2),
                                sym_class=symmetry_class,
                            )[0]
                            # print  "  keepset  ",total_nima,len(keepset),itmp,keepset
                            otherset = set(range(total_nima)) - set(keepset)
                            otherset = [i for i in otherset]
                            keepset.sort()
                            otherset.sort()
                            newparms1 = [None] * total_nima
                            newparms2 = [None] * total_nima
                            for i in keepset:
                                newparms1[i] = copy.deepcopy(GA[ipl[ip]][1][i])
                                newparms2[i] = copy.deepcopy(
                                    GA[ipl[(ip + 1) % number_of_runs]][1][i]
                                )
                            for i in otherset:
                                newparms1[i] = copy.deepcopy(
                                    GA[ipl[(ip + 1) % number_of_runs]][1][i]
                                )
                                newparms2[i] = copy.deepcopy(GA[ipl[ip]][1][i])
                            # print "  PRINTOUT SHUFFLED   ",ipl[ip],ipl[ip+1]
                            """Multiline Comment4"""

                            #  Put mutated params on one list
                            all_params.append(copy.deepcopy(newparms1))
                            all_params.append(copy.deepcopy(newparms2))
                        all_params = all_params[:number_of_runs]
                        #  Try this 02/03/2015 PAP
                        #  Always mutate the first ones
                        #  for half of projections 'mirror' them by adding 180 to psi
                        keepset = max(1, int(0.25 * number_of_runs))
                        # ipl = range(total_nima)
                        # shuffle(ipl)
                        # ipl = ipl[:total_nima//2]
                        for i in range(keepset):
                            all_params[0][i][2] += 180.0
                            #  Always reseed the last ones
                            all_params[-1 - i] = [
                                [
                                    random.random() * 360.0,
                                    random.random() * 180.0,
                                    random.random() * 360.0,
                                    0.0,
                                    0.0,
                                ]
                                for j in range(total_nima)
                            ]

                terminate = sp_utilities.wrap_mpi_bcast(terminate, main_node, mpi_comm)
                if not terminate:

                    storevol = True

                    # Send params back
                    if myid == 0:
                        for i in range(number_of_runs):
                            sr = mpi_subroots[i]
                            if sr == myid:
                                params = all_params[i]
                            else:
                                sp_utilities.wrap_mpi_send(all_params[i], sr, mpi_comm)
                    else:
                        if mpi_subrank == 0:
                            params = sp_utilities.wrap_mpi_recv(0, mpi_comm)

                    params = sp_utilities.wrap_mpi_bcast(params, 0, mpi_subcomm)
                    for i in range(nima):
                        sp_utilities.set_params_proj(data[i], params[i])
                    """Multiline Comment5"""

                    # =========================================================================
                    #
                    mpi.mpi_barrier(mpi_comm)
                    if myid == main_node:
                        log.add(
                            "Time of orientation and mutations = %f\n"
                            % (time.time() - start_time)
                        )
                        start_time = time.time()
                # else:  continue
            if not terminate:
                # =========================================================================
                # volume reconstruction
                mpi.mpi_barrier(mpi_comm)
                if myid == main_node:
                    start_time = time.time()
                vol = do_volume(
                    data[image_start:image_end], ali3d_options, 0, mpi_subcomm
                )

                if len(ali3d_options.moon_elimination) > 0:
                    vol = sp_utilities.eliminate_moons(
                        vol, ali3d_options.moon_elimination
                    )

                if mpi_subrank == 0:
                    L2 = vol.cmp(
                        "dot",
                        vol,
                        dict(
                            negative=0,
                            mask=sp_utilities.model_circle(last_ring, nx, nx, nx),
                        ),
                    )
                    # if myid == 2:  print  " Right after reconstruction L2", myid, L2,[get_params_proj(data[i]) for i in xrange(4)]
                    # print  " Right after reconstruction L2", myid, total_iter,L2
                    # if storevol:   vol.write_image("mutated%04d%04d.hdf"%(myid,total_iter))

                # log
                if myid == main_node:
                    log.add(
                        "3D reconstruction time = %f\n" % (time.time() - start_time)
                    )
                    start_time = time.time()

            """Multiline Comment6"""
            """Multiline Comment7"""

    # if mpi_subrank == 0:
    # 	for im in xrange(len(data)):
    # 		print  "VIPER1 peak ",myid,im,data[im].get_attr("ID"), data[im].get_attr("previousmax"), get_params_proj(data[im])
    # =========================================================================
    mpi.mpi_comm_free(mpi_subcomm)

    if myid == main_node:
        log.add("Finished VIPER1")
        return GA[0][1]
    else:
        return None  # results for the other processes


# parameters: list of (all) projections | reference volume | ...


def ali3d_multishc_2(
    stack, ref_vol, ali3d_options, symmetry_class, mpi_comm=None, log=None
):

    ir = ali3d_options.ir
    rs = ali3d_options.rs
    ou = ali3d_options.ou
    xr = ali3d_options.xr
    yr = ali3d_options.yr
    ts = ali3d_options.ts
    an = ali3d_options.an
    sym = ali3d_options.sym
    sym = sym[0].lower() + sym[1:]
    delta = ali3d_options.delta
    center = ali3d_options.center
    CTF = ali3d_options.CTF
    ref_a = ali3d_options.ref_a

    # Optionally restrict out-of-plane angle
    if hasattr(ali3d_options, "theta1"):
        theta1 = ali3d_options.theta1
    else:
        theta1 = -1.0

    if hasattr(ali3d_options, "theta2"):
        theta2 = ali3d_options.theta2
    else:
        theta2 = -1.0

    if hasattr(ali3d_options, "method"):
        method = ali3d_options.method
    else:
        method = "S"

    if mpi_comm == None:
        mpi_comm = mpi.MPI_COMM_WORLD

    if log == None:
        log = sp_logger.Logger()

    number_of_proc = mpi.mpi_comm_size(mpi_comm)
    myid = mpi.mpi_comm_rank(mpi_comm)
    main_node = 0

    if myid == main_node:
        log.add("Start VIPER2")

    xrng = sp_utilities.get_input_from_string(xr)
    if yr == "-1":
        yrng = xrng
    else:
        yrng = sp_utilities.get_input_from_string(yr)
    step = sp_utilities.get_input_from_string(ts)
    delta = sp_utilities.get_input_from_string(delta)
    lstp = min(len(xrng), len(yrng), len(step), len(delta))

    symmetry_class = sp_fundamentals.symclass(sym)

    # if an != "-1":
    # 	ERROR("Option an not used","VIPER1",1,myid)
    """Multiline Comment8"""

    first_ring = int(ir)
    rstep = int(rs)
    last_ring = int(ou)
    max_iter = int(ali3d_options.maxit2)
    center = int(center)

    vol = ref_vol
    nx = vol.get_xsize()
    if last_ring < 0:
        last_ring = int(old_div(nx, 2)) - 2

    numr = sp_alignment.Numrinit(first_ring, last_ring, rstep, "F")
    mask2D = sp_utilities.model_circle(last_ring, nx, nx) - sp_utilities.model_circle(
        first_ring, nx, nx
    )

    if myid == main_node:
        list_of_particles = list(range(len(stack)))
        total_nima = len(list_of_particles)
    else:
        list_of_particles = None
        total_nima = None
    total_nima = sp_utilities.wrap_mpi_bcast(total_nima, main_node, mpi_comm)
    list_of_particles = sp_utilities.wrap_mpi_bcast(
        list_of_particles, main_node, mpi_comm
    )

    # old_mpi_comm = mpi_comm
    # mpi_size = mpi_comm_size(mpi_comm)

    ## if there are fewer images than processors then split processors in 2 groups
    ## one in which each processor analyzes one image, and another in which
    ## processors stay idle and wait for the other group to finish
    # if (mpi_size > total_nima):
    # if (myid < total_nima):
    # mpi_subcomm = mpi_comm_split(mpi_comm, 0, myid)
    # mpi_comm = mpi_subcomm
    # else:
    # mpi_subcomm = mpi_comm_split(mpi_comm, 1, myid - total_nima)
    # mpi_barrier(mpi_comm)
    # return None, None, None, None

    number_of_proc = mpi.mpi_comm_size(mpi_comm)
    myid = mpi.mpi_comm_rank(mpi_comm)

    image_start, image_end = sp_applications.MPI_start_end(
        total_nima, number_of_proc, myid
    )
    # create a list of images for each node
    list_of_particles = list_of_particles[image_start:image_end]
    nima = len(list_of_particles)

    data = [stack[im] for im in list_of_particles]
    for im in range(nima):
        data[im].set_attr("ID", list_of_particles[im])
        ctf_applied = data[im].get_attr_default("ctf_applied", 0)
        if CTF and ctf_applied == 0:
            ctf_params = data[im].get_attr("ctf")
            st = EMAN2_cppwrap.Util.infomask(data[im], mask2D, False)
            data[im] -= st[0]
            data[im] = sp_filter.filt_ctf(data[im], ctf_params)
            data[im].set_attr("ctf_applied", 1)

    # =========================================================================
    # adjust params to references, calculate psi+shifts, calculate previousmax
    # qvol = volume_reconstruction(data, ali3d_options, mpi_comm)
    qvol = do_volume(data, ali3d_options, 0, mpi_comm)
    # log

    if myid == main_node:
        start_time = time.time()
        # qvol.write_image("vitera.hdf")
        L2 = qvol.cmp(
            "dot",
            qvol,
            dict(negative=0, mask=sp_utilities.model_circle(last_ring, nx, nx, nx)),
        )
        log.add(
            "3D reconstruction time = %f\n" % (time.time() - start_time),
            " START  L2 norm:  %f" % L2,
        )
        start_time = time.time()
    del qvol

    cnx = old_div(nx, 2) + 1
    cny = old_div(nx, 2) + 1
    wr_four = sp_alignment.ringwe(numr, "F")
    qv = old_div(math.pi, 180.0)
    volft, kb = sp_projection.prep_vol(ref_vol)
    for im in range(nima):
        phi, theta, psi, tx, ty = sp_utilities.get_params_proj(data[im])
        refrings = sp_projection.prgs(volft, kb, [phi, theta, psi, 0.0, 0.0])
        refrings = EMAN2_cppwrap.Util.Polar2Dm(refrings, cnx, cny, numr, "F")
        EMAN2_cppwrap.Util.Normalize_ring(refrings, numr, 0)
        EMAN2_cppwrap.Util.Frngs(refrings, numr)
        EMAN2_cppwrap.Util.Applyws(refrings, numr, wr_four)
        """Multiline Comment9"""

        # print "orin  ",data[im].get_attr("ID"),get_params_proj(data[im])
        # peak, pixer = proj_ali_incore_local(data[im],refrings, [[phi, theta], [phi, theta]], numr,0.0,0.0,1.0,delta[0]/4)
        cimage = EMAN2_cppwrap.Util.Polar2Dm(data[im], cnx, cny, numr, "F")
        EMAN2_cppwrap.Util.Normalize_ring(cimage, numr, 0)
        EMAN2_cppwrap.Util.Frngs(cimage, numr)
        retvals = EMAN2_cppwrap.Util.Crosrng_e(refrings, cimage, numr, 0, 0.0)
        data[im].set_attr("previousmax", retvals["qn"])
        # print  " VIPER2  peak ",myid,im,data[im].get_attr("ID"), retvals["qn"],get_params_proj(data[im])

    """Multiline Comment10"""

    del volft
    # volume reconstruction
    mpi.mpi_barrier(mpi_comm)
    if myid == main_node:
        log.add(
            "Time to calculate first psi+shifts+previousmax: %f\n"
            % (time.time() - start_time)
        )
        start_time = time.time()
    ref_vol = do_volume(data, ali3d_options, 0, mpi_comm)

    if myid == main_node:
        # ref_vol.write_image("viterb.hdf")
        L2 = ref_vol.cmp(
            "dot",
            ref_vol,
            dict(negative=0, mask=sp_utilities.model_circle(last_ring, nx, nx, nx)),
        )
        log.add(
            "3D reconstruction time = %f\n" % (time.time() - start_time),
            "   L2 norm:  %f" % L2,
        )
        start_time = time.time()

    # =========================================================================

    pixer = [0.0] * nima
    par_r = [[] for im in list_of_particles]
    cs = [0.0] * 3
    total_iter = 0
    # do the projection matching
    for N_step in range(lstp):

        terminate = 0
        Iter = 0
        while Iter < max_iter and terminate == 0:

            Iter += 1
            total_iter += 1

            mpi.mpi_barrier(mpi_comm)
            if myid == main_node:
                log.add(
                    "ITERATION #%3d,  inner iteration #%3d\nDelta = %4.1f, xrange = %5.2f, yrange = %5.2f, step = %5.2f\n"
                    % (
                        total_iter,
                        Iter,
                        delta[N_step],
                        xrng[N_step],
                        yrng[N_step],
                        step[N_step],
                    )
                )
                start_time = time.time()

            # =========================================================================
            # build references
            volft, kb = sp_projection.prep_vol(vol)
            #  For the local SHC it is essential reference projections have psi zero, as otherwise it will get messed up.
            reference_angles = symmetry_class.even_angles(
                delta[N_step],
                phiEqpsi="Zero",
                theta1=theta1,
                theta2=theta2,
                method=method,
            )
            refrings = sp_alignment.prepare_refrings(
                volft, kb, nx, -1.0, reference_angles, "", numr, MPI=mpi_comm
            )
            del volft, kb
            # =========================================================================

            if myid == main_node:
                log.add("Time to prepare rings: %f\n" % (time.time() - start_time))
                start_time = time.time()

            mpi.mpi_barrier(mpi_comm)
            if myid == main_node:
                start_time = time.time()
            # =========================================================================
            # alignment
            number_of_checked_refs = 0
            for im in range(nima):
                # peak, pixer[im], checked_refs, number_of_peaks = shc_multi(data[im],refrings,numr,xrng[N_step],yrng[N_step],step[N_step],an[N_step], number_of_runs=number_of_runs)
                # previousmax is set in shc
                peak, pixer[im], checked_refs, iref = sp_alignment.shc(
                    data[im],
                    refrings,
                    [[1.0, 1.0]],
                    numr,
                    xrng[N_step],
                    yrng[N_step],
                    step[N_step],
                    sym="nomirror",
                )  # cannot use 'an' here
                number_of_checked_refs += checked_refs
            # =========================================================================
            mpi.mpi_barrier(mpi_comm)
            if myid == main_node:
                log.add("Time of alignment = %f\n" % (time.time() - start_time))
                start_time = time.time()

            # =========================================================================
            # output pixel errors, check stop criterion
            all_pixer = sp_utilities.wrap_mpi_gatherv(pixer, 0, mpi_comm)
            total_checked_refs = sp_utilities.wrap_mpi_gatherv(
                [number_of_checked_refs], main_node, mpi_comm
            )
            terminate = 0
            if myid == main_node:
                total_checked_refs = sum(total_checked_refs)
                lhist = 20
                region, histo = sp_statistics.hist_list(all_pixer, lhist)
                log.add("=========================")
                for lhx in range(lhist):
                    msg = " %10.3f     %7d" % (region[lhx], histo[lhx])
                    log.add(msg)
                if (max(all_pixer) < 0.5) and (
                    old_div(sum(all_pixer), total_nima) < 0.05
                ):
                    terminate = 1
            terminate = sp_utilities.wrap_mpi_bcast(terminate, main_node, mpi_comm)
            # =========================================================================

            # =========================================================================
            # centering
            if center == -1 and sym[0] == "c":
                cs[0], cs[1], cs[2], dummy, dummy = sp_utilities.estimate_3D_center_MPI(
                    data, total_nima, myid, number_of_proc, main_node, mpi_comm=mpi_comm
                )
                if myid == main_node:
                    msg = (
                        " Average center x = %10.3f        Center y = %10.3f        Center z = %10.3f\n"
                        % (cs[0], cs[1], cs[2])
                    )
                    log.add(msg)
                if int(sym[1]) > 1:
                    cs[0] = cs[1] = 0.0
                    if myid == main_node:
                        log.add(
                            "For symmetry group cn (n>1), we only center the volume in z-direction\n"
                        )
                cs = mpi.mpi_bcast(cs, 3, mpi.MPI_FLOAT, main_node, mpi_comm)
                cs = [-float(cs[0]), -float(cs[1]), -float(cs[2])]
                sp_utilities.rotate_3D_shift(data, cs)
            # =========================================================================

            # =========================================================================
            # volume reconstruction
            mpi.mpi_barrier(mpi_comm)
            if myid == main_node:
                start_time = time.time()
            # vol = volume_reconstruction(data, ali3d_options, mpi_comm)
            vol = do_volume(data, ali3d_options, 0, mpi_comm)
            # log
            if myid == main_node:
                # vol.write_image("viter%03d.hdf"%total_iter)
                L2 = vol.cmp(
                    "dot",
                    vol,
                    dict(
                        negative=0,
                        mask=sp_utilities.model_circle(last_ring, nx, nx, nx),
                    ),
                )
                log.add(
                    "3D reconstruction time = %f\n" % (time.time() - start_time),
                    "   L2 norm:  %f" % L2,
                )
                start_time = time.time()
            # =========================================================================

    # =========================================================================
    # gather parameters
    params = []
    previousmax = []
    for im in data:
        t = sp_utilities.get_params_proj(im)
        p = im.get_attr("previousmax")
        params.append([t[0], t[1], t[2], t[3], t[4]])
        previousmax.append(p)
    assert nima == len(params)
    params = sp_utilities.wrap_mpi_gatherv(params, 0, mpi_comm)
    if myid == 0:
        assert total_nima == len(params)
    previousmax = sp_utilities.wrap_mpi_gatherv(previousmax, 0, mpi_comm)
    if myid == 0:
        assert total_nima == len(previousmax)

    par_r = sp_utilities.wrap_mpi_gatherv(par_r, 0, mpi_comm)

    ## if there are fewer images than processors then synchronize
    ## with the other group of processors that did not do any work
    # if (mpi_size > total_nima):
    # if (myid < no_of_images):
    # mpi_comm = old_mpi_comm
    # mpi_barrier(mpi_comm)

    if myid == main_node:
        log.add("Finished VIPER2")
        return params, vol, previousmax, par_r
    else:
        return None, None, None, None  # results for the other processes


# data - projections (scattered between cpus)
# options - the same for all cpus
# return - volume the same for all cpus


"""Multiline Comment11"""


# all_projs and subset must be set only for root (MPI rank == 0)
# remaining parameters must be set for all
# size of mpi_communicator must be >= runs_count


def multi_shc(
    all_projs, subset, runs_count, ali3d_options, mpi_comm, log=None, ref_vol=None
):
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

    mpi_rank = mpi.mpi_comm_rank(mpi_comm)
    mpi_size = mpi.mpi_comm_size(mpi_comm)

    global BATCH, MPI
    sp_global_def.BATCH = True
    sp_global_def.MPI = True
    if mpi_size < runs_count:
        sp_global_def.ERROR("multi_shc", "mpi_size < runs_count", 1, mpi_rank)

    if log == None:
        log = sp_logger.Logger()

    #  Initialize symmetries
    symmetry_class = sp_fundamentals.symclass(ali3d_options.sym)

    error = 0
    projections = []
    if mpi_rank == 0:
        # Optionally restrict out-of-plane angle
        if hasattr(ali3d_options, "theta1"):
            theta1 = ali3d_options.theta1
        else:
            theta1 = -1.0

        if hasattr(ali3d_options, "theta2"):
            theta2 = ali3d_options.theta2
        else:
            theta2 = -1.0

        if hasattr(ali3d_options, "method"):
            method = ali3d_options.ref_a
        else:
            method = "S"

        if not hasattr(ali3d_options, "filament_width"):
            ali3d_options.filament_width = -1

        prms = symmetry_class.even_angles(
            float(ali3d_options.delta), theta1=theta1, theta2=theta2, method=method
        )
        sp_utilities.write_text_row(prms, log.prefix + "initangles.txt")

        if len(prms) < len(subset):
            error = 1
        else:
            random.shuffle(prms)
            for i in subset:
                prms[i][2] = random.random() * 360.0
                sp_utilities.set_params_proj(all_projs[i], prms[i] + [0.0, 0.0])
                # all_projs[i].set_attr("stable", 0)
                all_projs[i].set_attr("previousmax", -1.0e23)
                projections.append(all_projs[i])
        del prms

    error = sp_utilities.bcast_number_to_all(error, source_node=0)
    if error == 1:
        sp_global_def.ERROR("Angular step too large, decrease delta", 1, mpi_rank)

    ###from sys import exit
    # if mpi_rank == 0:   print "  NEW   ",mpi_rank
    ###exit()
    projections = sp_utilities.wrap_mpi_bcast(projections, 0, mpi_comm)

    n_projs = len(projections)

    if ref_vol == None:
        if mpi_size > n_projs:
            working = int(not (mpi_rank < n_projs))
            mpi_subcomm = mpi.mpi_comm_split(
                mpi_comm, working, mpi_rank - working * n_projs
            )
            mpi_subsize = mpi.mpi_comm_size(mpi_subcomm)
            mpi_subrank = mpi.mpi_comm_rank(mpi_subcomm)
            if mpi_rank < n_projs:
                proj_begin, proj_end = sp_applications.MPI_start_end(
                    n_projs, mpi_subsize, mpi_subrank
                )
                ref_vol = do_volume(
                    projections[proj_begin:proj_end],
                    ali3d_options,
                    0,
                    mpi_comm=mpi_subcomm,
                )
            else:
                nx = projections[0].get_xsize()
                ref_vol = sp_utilities.model_blank(nx, nx, nx)
            sp_utilities.bcast_EMData_to_all(ref_vol, mpi_rank, 0, comm=mpi_comm)
        else:
            proj_begin, proj_end = sp_applications.MPI_start_end(
                n_projs, mpi_size, mpi_rank
            )
            ref_vol = do_volume(
                projections[proj_begin:proj_end], ali3d_options, 0, mpi_comm=mpi_comm
            )

    # Each node keeps all projection data, this would not work for large datasets
    out_params = ali3d_multishc(
        projections,
        ref_vol,
        ali3d_options,
        symmetry_class,
        mpi_comm=mpi_comm,
        log=log,
        number_of_runs=runs_count,
    )
    """Multiline Comment12"""

    if mpi_rank == 0:
        """Multiline Comment13"""

        temp = []
        for i in range(n_projs):
            sp_utilities.set_params_proj(projections[i], out_params[i])
        sp_utilities.write_text_row(out_params, log.prefix + "refparams2.txt")
        """Multiline Comment14"""
    """Multiline Comment15"""
    # proj_begin, proj_end  = MPI_start_end(n_projs, mpi_size, mpi_rank)

    projections = sp_utilities.wrap_mpi_bcast(projections, 0, mpi_comm)

    if mpi_size > n_projs:
        working = int(not (mpi_rank < n_projs))
        mpi_subcomm = mpi.mpi_comm_split(
            mpi_comm, working, mpi_rank - working * n_projs
        )
        mpi_subsize = mpi.mpi_comm_size(mpi_subcomm)
        mpi_subrank = mpi.mpi_comm_rank(mpi_subcomm)
        if mpi_rank < n_projs:
            proj_begin, proj_end = sp_applications.MPI_start_end(
                n_projs, mpi_subsize, mpi_subrank
            )
            ref_vol = do_volume(
                projections[proj_begin:proj_end], ali3d_options, 0, mpi_comm=mpi_subcomm
            )
        else:
            nx = projections[0].get_xsize()
            ref_vol = sp_utilities.model_blank(nx, nx, nx)
        sp_utilities.bcast_EMData_to_all(ref_vol, mpi_rank, 0, comm=mpi_comm)
    else:
        proj_begin, proj_end = sp_applications.MPI_start_end(
            n_projs, mpi_size, mpi_rank
        )
        ref_vol = do_volume(
            projections[proj_begin:proj_end], ali3d_options, 0, mpi_comm=mpi_comm
        )

    if mpi_rank == 0:
        ref_vol.write_image(log.prefix + "refvol2.hdf")
        nx = ref_vol.get_xsize()
        L2 = ref_vol.cmp(
            "dot",
            ref_vol,
            dict(
                negative=0, mask=sp_utilities.model_circle(ali3d_options.ou, nx, nx, nx)
            ),
        )
        log.add(" L2 norm of reference volume:  %f" % L2)

        # Generate angular distribution
        if hasattr(ali3d_options, "dpi"):
            dpi = ali3d_options.dpi
        else:
            dpi = 72
        pixel_size = (
            1
        )  # Not going to upscale to the original dimensions, so in Chimera open reconstruction at 1 Angstrom/voxel, etc.

        sp_utilities.angular_distribution(
            params_file=log.prefix + "refparams2.txt",
            output_folder=log.prefix,
            prefix="refvol2_angdist",
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

    if mpi_size > n_projs:
        working = int(not (mpi_rank < n_projs))
        mpi_subcomm = mpi.mpi_comm_split(
            mpi_comm, working, mpi_rank - working * n_projs
        )
        mpi_subsize = mpi.mpi_comm_size(mpi_subcomm)
        mpi_subrank = mpi.mpi_comm_rank(mpi_subcomm)
        if mpi_rank < n_projs:
            out_params, out_vol, previousmax, out_r = ali3d_multishc_2(
                projections,
                ref_vol,
                ali3d_options,
                symmetry_class,
                mpi_comm=mpi_subcomm,
                log=log,
            )
        else:
            out_params = None
            out_vol = None
        mpi.mpi_barrier(mpi_comm)
    else:
        out_params, out_vol, previousmax, out_r = ali3d_multishc_2(
            projections,
            ref_vol,
            ali3d_options,
            symmetry_class,
            mpi_comm=mpi_comm,
            log=log,
        )

    if mpi_rank == 0:
        sp_utilities.write_text_file(previousmax, log.prefix + "previousmax.txt")
        sp_utilities.write_text_row(out_params, log.prefix + "params.txt")
        sp_utilities.drop_image(out_vol, log.prefix + "volf.hdf")

        # Generate angular distribution
        independent_run_dir = log.prefix
        sp_global_def.sxprint("independent_run_dir", independent_run_dir)

        params_file = log.prefix + "params.txt"
        output_folder = independent_run_dir
        prefix = "volf_angdist"  # will overwrite input parameters file if blank
        delta = float(ali3d_options.delta)
        symmetry = ali3d_options.sym
        if hasattr(ali3d_options, "dpi"):
            dpi = ali3d_options.dpi
        else:
            dpi = 72

        # Not going to upscale to the original dimensions, so in Chimera open reconstruction at 1 Angstrom/voxel, etc.
        pixel_size = 1
        box_size = sp_utilities.get_im(os.path.join(log.prefix, "volf.hdf")).get_xsize()

        sp_utilities.angular_distribution(
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

    return out_params, out_vol, None  # , out_peaks


def mirror_and_reduce_dsym(params, indexes, symmetry_class):
    #  Input params contains multiple datasets [ p0, p1, ...]
    #  We treat p0 as a reference
    # For D symmetry there are two equivalent positions that agree with given Dn symmetry
    #  The second is rotated by 360/n degrees.

    sc = len(params)
    ns = len(params[0])

    #  bbdb is 360.0/nsym and indicates position of the second symmetry
    bbdb = old_div(360.0, symmetry_class.nsym)
    symphi = old_div(360.0, symmetry_class.nsym * 2)
    #  For each set we have four positions to consider: straight, straight psi mirrored, phi+bdb, phi+bdb and psi mirrored
    for i in range(1, sc):
        psi_diff = sp_pixel_error.angle_diff(
            [params[i][j][2] for j in indexes], [params[0][j][2] for j in indexes]
        )
        #  adjust psi if necessary
        if abs(psi_diff - 180.0) < 90.0:
            for j in range(ns):
                # apply mirror
                params[i][j][2] = (params[i][j][2] + 180.0) % 360.0
        # Check which one of the two possible symmetry positions is closer
        #  compute angular errors including symmetry
        per1 = 0.0
        per2 = 0.0
        temp = [None] * ns
        for j in indexes:
            neisym = symmetry_class.symmetry_neighbors([params[i][j][:3]])
            dmin = 180.0
            for q in neisym:
                dmin = min(dmin, getang3(q, params[0][j]))
            per1 += dmin
            temp = symmetry_class.reduce_anglesets(
                [params[i][j][0] + bbdb, params[i][j][1], params[i][j][2]]
            )
            neisym = symmetry_class.symmetry_neighbors([temp])
            dmin = 180.0
            for q in neisym:
                dmin = min(dmin, getang3(q, params[0][j]))
            per2 += dmin

        if per2 < per1:
            for j in range(ns):
                temp = symmetry_class.reduce_anglesets(
                    [params[i][j][0] + bbdb, params[i][j][1], params[i][j][2]]
                )
                params[i][j] = [
                    temp[0],
                    temp[1],
                    temp[2],
                    params[i][j][3],
                    params[i][j][4],
                ]


def do_volume(data, options, iter, mpi_comm):

    myid = mpi.mpi_comm_rank(mpi_comm)
    sym = options.sym
    sym = sym[0].lower() + sym[1:]
    npad = options.npad
    CTF = options.CTF
    snr = options.snr
    # =========================================================================
    # volume reconstruction
    if type(data) == list:
        if CTF:
            vol = sp_reconstruction.recons3d_4nn_ctf_MPI(
                myid, data, snr, symmetry=sym, npad=npad, mpi_comm=mpi_comm
            )
        else:
            vol = sp_reconstruction.recons3d_4nn_MPI(
                myid, data, symmetry=sym, snr=snr, npad=npad, mpi_comm=mpi_comm
            )
    else:
        vol = data

    if myid == 0:
        nx = vol.get_xsize()
        if options.filament_width != -1:
            mask3D = sp_utilities.model_cylinder(
                old_div(int(options.filament_width * options.resample_ratio + 0.5), 2),
                nx,
                nx,
                nx,
            )
            options.mask3D = mask3D
        elif options.mask3D == None:
            last_ring = int(options.ou)
            mask3D = sp_utilities.model_circle(last_ring, nx, nx, nx)
        elif options.mask3D == "auto":
            mask3D = sp_morphology.adaptive_mask(vol)
        else:
            if isinstance(options.mask3D, (bytes, str)):
                mask3D = sp_utilities.get_im(options.mask3D)
            else:
                mask3D = (options.mask3D).copy()
            nxm = mask3D.get_xsize()
            if nx != nxm:
                mask3D = EMAN2_cppwrap.Util.window(
                    sp_fundamentals.rot_shift3D(
                        mask3D, scale=old_div(float(nx), float(nxm))
                    ),
                    nx,
                    nx,
                    nx,
                )
                nxm = mask3D.get_xsize()
                assert nx == nxm

        stat = EMAN2_cppwrap.Util.infomask(vol, mask3D, False)
        vol -= stat[0]
        EMAN2_cppwrap.Util.mul_scalar(vol, old_div(1.0, stat[1]))
        vol = sp_morphology.threshold(vol)
        # Util.mul_img(vol, mask3D)
        if options.pwreference:
            rt = sp_utilities.read_text_file(options.pwreference)
            sp_fundamentals.fftip(vol)
            ro = sp_fundamentals.rops_table(vol)
            #  Here unless I am mistaken it is enough to take the beginning of the reference pw.
            for i in range(1, len(ro)):
                ro[i] = (old_div(rt[i], ro[i])) ** 0.5
            if type(options.fl) == list:
                vol = sp_fundamentals.fft(
                    sp_filter.filt_table(sp_filter.filt_table(vol, options.fl), ro)
                )
            else:
                vol = sp_fundamentals.fft(
                    sp_filter.filt_table(
                        sp_filter.filt_tanl(vol, options.fl, options.aa), ro
                    )
                )
        else:
            if type(options.fl) == list:
                vol = sp_filter.filt_table(vol, options.fl)
            else:
                vol = sp_filter.filt_tanl(vol, options.fl, options.aa)
        stat = EMAN2_cppwrap.Util.infomask(vol, mask3D, False)
        vol -= stat[0]
        EMAN2_cppwrap.Util.mul_scalar(vol, old_div(1.0, stat[1]))
        vol = sp_morphology.threshold(vol)
        vol = sp_filter.filt_btwl(vol, 0.38, 0.5)
        EMAN2_cppwrap.Util.mul_img(vol, mask3D)
        del mask3D
        # vol.write_image('toto%03d.hdf'%iter)
    # broadcast volume
    sp_utilities.bcast_EMData_to_all(vol, myid, 0, comm=mpi_comm)
    # =========================================================================
    return vol


"""Multiline Comment20"""


"""Multiline Comment21"""


"""Multiline Comment22"""
"""Multiline Comment23"""
from builtins import range
