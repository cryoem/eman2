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
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holfds
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
import EMAN2db
import EMAN2
import mpi
import numpy
import os
import random
from . import sp_alignment
from . import sp_filter
from . import sp_fundamentals
from . import sp_global_def
from . import sp_logger
from . import sp_morphology
from . import sp_pixel_error
from . import sp_projection
from . import sp_reconstruction
from . import sp_statistics
from . import sp_user_functions
from . import sp_utilities
import sys
import time


def ali2d_MPI(
    stack,
    outdir,
    maskfile=None,
    ir=1,
    ou=-1,
    rs=1,
    xr="4 2 1 1",
    yr="-1",
    ts="2 1 0.5 0.25",
    nomirror=False,
    dst=0.0,
    center=-1,
    maxit=0,
    CTF=False,
    snr=1.0,
    Fourvar=False,
    Ng=-1,
    user_func_name="ref_ali2d",
    CUDA=False,
    GPUID="",
    random_method="",
):

    number_of_proc = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
    myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
    main_node = 0

    ftp = sp_utilities.file_type(stack)

    if outdir:
        if os.path.exists(outdir):
            sp_global_def.ERROR(
                "Output directory exists, please change the name and restart the program",
                "ali2d_MPI",
                1,
                myid,
            )
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    if myid == main_node:
        if outdir:
            os.mkdir(outdir)
        sp_global_def.LOGFILE = os.path.join(outdir, sp_global_def.LOGFILE)
        sp_utilities.print_begin_msg("ali2d_MPI")

    xrng = sp_utilities.get_input_from_string(xr)
    if yr == "-1":
        yrng = xrng
    else:
        yrng = sp_utilities.get_input_from_string(yr)
    step = sp_utilities.get_input_from_string(ts)

    first_ring = int(ir)
    last_ring = int(ou)
    rstep = int(rs)
    max_iter = int(maxit)

    if max_iter == 0:
        max_iter = 10
        auto_stop = True
    else:
        auto_stop = False

    if myid == main_node:
        if ftp == "bdb":
            dummy = EMAN2db.db_open_dict(stack, True)
        # horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
        # active = EMUtil.get_all_attributes(stack, 'active')
        # list_of_particles = []
        # for im in xrange(len(active)):
        # 	if active[im]:  list_of_particles.append(im)
        # del active
        # nima = len(list_of_particles)

        nima = EMAN2_cppwrap.EMUtil.get_image_count(stack)
        list_of_particles = list(range(nima))

    else:
        nima = 0
    nima = sp_utilities.bcast_number_to_all(nima, source_node=main_node)

    if myid != main_node:
        list_of_particles = [-1] * nima
    list_of_particles = sp_utilities.bcast_list_to_all(
        list_of_particles, myid, source_node=main_node
    )

    image_start, image_end = MPI_start_end(nima, number_of_proc, myid)
    list_of_particles = list_of_particles[image_start:image_end]

    if Ng == -1:
        Ng = nima
    elif Ng == -2:
        Ng = int(0.98 * nima)

    # read nx and ctf_app (if CTF) and broadcast to all nodes
    if myid == main_node:
        ima = EMAN2_cppwrap.EMData()
        ima.read_image(stack, list_of_particles[0], True)
        nx = ima.get_xsize()
        if CTF:
            ctf_app = ima.get_attr_default("ctf_applied", 0)
        del ima
    else:
        nx = 0
        if CTF:
            ctf_app = 0
    nx = sp_utilities.bcast_number_to_all(nx, source_node=main_node)
    if CTF:
        ctf_app = sp_utilities.bcast_number_to_all(ctf_app, source_node=main_node)
        if ctf_app > 0:
            sp_global_def.ERROR("data cannot be ctf-applied", "ali2d_MPI", 1, myid)
        phase_flip = True
    else:
        phase_flip = False
    CTF = False

    # default value for the last ring
    if last_ring == -1:
        last_ring = old_div(nx, 2) - 2

    if last_ring + max([max(xrng), max(yrng)]) > old_div((nx - 1), 2):
        sp_global_def.ERROR(
            "Shift or radius is too large - particle crosses image boundary",
            "ali2d_MPI",
            1,
        )

    if CUDA:
        GPUID = sp_utilities.get_input_from_string(GPUID)
        GPUID = list(map(int, GPUID))

    if myid == main_node:
        sp_utilities.print_msg("Input stack                 : %s\n" % (stack))
        # horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
        # print_msg("Number of active images     : %d\n"%(nima))
        sp_utilities.print_msg("Number of images            : %d\n" % (nima))
        sp_utilities.print_msg("Output directory            : %s\n" % (outdir))
        sp_utilities.print_msg("Inner radius                : %i\n" % (first_ring))
        sp_utilities.print_msg("Outer radius                : %i\n" % (last_ring))
        sp_utilities.print_msg("Ring step                   : %i\n" % (rstep))
        sp_utilities.print_msg("X search range              : %s\n" % (xrng))
        sp_utilities.print_msg("Y search range              : %s\n" % (yrng))
        sp_utilities.print_msg("Translational step          : %s\n" % (step))
        sp_utilities.print_msg("Disable checking mirror     : %s\n" % (nomirror))
        sp_utilities.print_msg("Discrete angle used         : %d\n" % (dst))
        sp_utilities.print_msg("Center type                 : %i\n" % (center))
        sp_utilities.print_msg("Maximum iteration           : %i\n" % (max_iter))
        sp_utilities.print_msg("Use Fourier variance        : %s\n" % (Fourvar))
        # print_msg("Number of groups            : %d\n"%(Ng))
        sp_utilities.print_msg("CTF correction              : %s\n" % (CTF))
        sp_utilities.print_msg("Phase flip                  : %s\n" % (phase_flip))
        sp_utilities.print_msg("Signal-to-Noise Ratio       : %f\n" % (snr))
        if auto_stop:
            sp_utilities.print_msg("Stop iteration with         : criterion\n")
        else:
            sp_utilities.print_msg("Stop iteration with         : maxit\n")

        user_func = sp_user_functions.factory[user_func_name]

        sp_utilities.print_msg("User function               : %s\n" % (user_func_name))
        sp_utilities.print_msg("Number of processors used   : %d\n" % (number_of_proc))
        sp_utilities.print_msg("Using CUDA                  : %s\n" % (CUDA))
        if CUDA:
            sp_utilities.print_msg("GPU IDs                     : %s\n" % (GPUID))

    if maskfile:
        if isinstance(maskfile, (bytes, str)):
            if myid == main_node:
                sp_utilities.print_msg(
                    "Maskfile                    : %s\n\n" % (maskfile)
                )
            mask = sp_utilities.get_image(maskfile)
        else:
            if myid == main_node:
                sp_utilities.print_msg(
                    "Maskfile                    : user provided in-core mask\n\n"
                )
            mask = maskfile
    else:
        if myid == main_node:
            sp_utilities.print_msg(
                "Maskfile                    : default, a circle with radius %i\n\n"
                % (last_ring)
            )
        mask = sp_utilities.model_circle(last_ring, nx, nx)

    cnx = old_div(nx, 2) + 1
    cny = cnx
    if random_method == "SCF":
        mode = "H"
    else:
        mode = "F"
    data = []
    if CTF:
        ctf_abs_sum = EMAN2_cppwrap.EMData(nx, nx, 1, False)
        ctf_2_sum = EMAN2_cppwrap.EMData(nx, nx, 1, False)
    else:
        ctf_2_sum = None

    if sp_global_def.CACHE_DISABLE:
        data = EMAN2_cppwrap.EMData.read_images(stack, list_of_particles)
    else:
        for i in range(number_of_proc):
            if myid == i:
                data = EMAN2_cppwrap.EMData.read_images(stack, list_of_particles)
            if ftp == "bdb":
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    if CUDA:
        nGPU = len(GPUID)
        GPUID = GPUID[myid % nGPU]
        R = CUDA_Aligner(GPUID)
        all_ali_params = []
        all_ctf_params = []

    for im in range(len(data)):
        data[im].set_attr("ID", list_of_particles[im])
        st = EMAN2_cppwrap.Util.infomask(data[im], mask, False)
        data[im] -= st[0]
        if CTF:
            ctf_params = data[im].get_attr("ctf")
            ctfimg = sp_morphology.ctf_img(nx, ctf_params)
            if CUDA:
                all_ctf_params.extend(
                    [
                        ctf_params.defocus,
                        ctf_params.cs,
                        ctf_params.voltage,
                        ctf_params.apix,
                        ctf_params.bfactor,
                        ctf_params.ampcont,
                    ]
                )
            EMAN2_cppwrap.Util.add_img2(ctf_2_sum, ctfimg)
            EMAN2_cppwrap.Util.add_img_abs(ctf_abs_sum, ctfimg)
        if CUDA:
            alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(data[im])
            all_ali_params.extend([alpha, sx, sy, mirror])
        if random_method == "SHC":
            data[im].set_attr("previousmax", 1.0e-23)
        if phase_flip:
            data[im] = sp_filter.filt_ctf(
                data[im], data[im].get_attr("ctf"), binary=True
            )

    if CTF:
        sp_utilities.reduce_EMData_to_root(ctf_2_sum, myid, main_node)
        sp_utilities.reduce_EMData_to_root(ctf_abs_sum, myid, main_node)
        if myid == main_node:
            adw_img = EMAN2_cppwrap.Util.mult_scalar(ctf_2_sum, snr)
            EMAN2_cppwrap.Util.div_filter(adw_img, ctf_abs_sum)
            EMAN2_cppwrap.Util.mul_scalar(adw_img, old_div(float(Ng - 1), (nima - 1)))
            adw_img += old_div(float(nima - Ng), (nima - 1))
    else:
        ctf_2_sum = None
    # startup
    numr = sp_alignment.Numrinit(
        first_ring, last_ring, rstep, mode
    )  # precalculate rings
    wr = sp_alignment.ringwe(numr, mode)

    if myid == main_node:
        # initialize data for the reference preparation function
        ref_data = [mask, center, None, None]
        sx_sum = 0.0
        sy_sum = 0.0
        a0 = -1.0e22

    recvcount = []
    disp = []
    for i in range(number_of_proc):
        ib, ie = MPI_start_end(nima, number_of_proc, i)
        recvcount.append(ie - ib)
        if i == 0:
            disp.append(0)
        else:
            disp.append(disp[i - 1] + recvcount[i - 1])

    again = 1
    total_iter = 0
    cs = [0.0] * 2

    if CUDA:
        RING_LENGTH = 2 ** (
            int(old_div(numpy.log(2 * numpy.pi * last_ring), numpy.log(2))) + 1
        )
        NRING = 2 ** (int(old_div(numpy.log(last_ring), numpy.log(2))) + 1)

    for N_step in range(len(xrng)):

        if CUDA:
            R.setup(
                len(data),
                nx,
                nx,
                RING_LENGTH,
                NRING,
                last_ring,
                step[N_step],
                int(old_div(xrng[N_step], step[N_step]) + 0.5),
                int(old_div(yrng[N_step], step[N_step]) + 0.5),
                CTF,
            )
            for im in range(len(data)):
                R.insert_image(data[im], im)
            if CTF:
                R.filter_stack(all_ctf_params)

        msg = "\nX range = %5.2f   Y range = %5.2f   Step = %5.2f\n" % (
            xrng[N_step],
            yrng[N_step],
            step[N_step],
        )
        if myid == main_node:
            sp_utilities.print_msg(msg)
        for Iter in range(max_iter):
            total_iter += 1
            if CUDA:
                ave1 = sp_utilities.model_blank(nx, nx)
                ave2 = sp_utilities.model_blank(nx, nx)
                R.sum_oe(all_ctf_params, all_ali_params, ave1, ave2)
                # Comment by Zhengfan Yang on 02/01/10
                # The reason for this step is that in CUDA 2-D FFT, the image is multipled by NX*NY times after
                # FFT and IFFT, so we want to decrease it such that the criterion is in line with non-CUDA version
                # However, this step is not mandatory.
                if CTF:
                    ave1 = old_div(ave1, (nx * 2) ** 2)
                    ave2 = old_div(ave2, (nx * 2) ** 2)
            else:
                ave1, ave2 = sp_statistics.sum_oe(
                    data, "a", CTF, EMAN2_cppwrap.EMData()
                )  # pass empty object to prevent calculation of ctf^2
            sp_utilities.reduce_EMData_to_root(ave1, myid, main_node)
            sp_utilities.reduce_EMData_to_root(ave2, myid, main_node)
            if myid == main_node:
                sp_utilities.print_msg("Iteration #%4d\n" % (total_iter))
                if CTF:
                    tavg_Ng = sp_fundamentals.fft(
                        EMAN2_cppwrap.Util.divn_filter(
                            EMAN2_cppwrap.Util.muln_img(
                                sp_fundamentals.fft(
                                    EMAN2_cppwrap.Util.addn_img(ave1, ave2)
                                ),
                                adw_img,
                            ),
                            ctf_2_sum,
                        )
                    )
                    tavg = sp_fundamentals.fft(
                        EMAN2_cppwrap.Util.divn_filter(
                            sp_fundamentals.fft(
                                EMAN2_cppwrap.Util.addn_img(ave1, ave2)
                            ),
                            ctf_2_sum,
                        )
                    )
                else:
                    tavg = old_div((ave1 + ave2), nima)
                if outdir:
                    tavg.write_image(os.path.join(outdir, "aqc.hdf"), total_iter - 1)
                    if CTF:
                        tavg_Ng.write_image(
                            os.path.join(outdir, "aqc_view.hdf"), total_iter - 1
                        )
                    frsc = sp_statistics.fsc_mask(
                        ave1,
                        ave2,
                        mask,
                        1.0,
                        os.path.join(outdir, "resolution%03d" % (total_iter)),
                    )
                else:
                    frsc = sp_statistics.fsc_mask(ave1, ave2, mask, 1.0)
            else:
                tavg = sp_utilities.model_blank(nx, nx)
            del ave1, ave2

            if Fourvar:
                sp_utilities.bcast_EMData_to_all(tavg, myid, main_node)
                vav, rvar = sp_statistics.varf2d_MPI(myid, data, tavg, mask, "a", CTF)

            if myid == main_node:
                if Fourvar:
                    tavg = sp_fundamentals.fft(
                        EMAN2_cppwrap.Util.divn_img(sp_fundamentals.fft(tavg), vav)
                    )
                    vav_r = EMAN2_cppwrap.Util.pack_complex_to_real(vav)
                    if outdir:
                        vav_r.write_image(
                            os.path.join(outdir, "varf.hdf"), total_iter - 1
                        )

                # a0 should increase; stop algorithm when it decreases.
                #     However, the result will depend on filtration, so it is not quite right.
                #  moved it here, so it is for unfiltered average and thus hopefully makes more sense
                a1 = tavg.cmp("dot", tavg, dict(negative=0, mask=ref_data[0]))
                msg = "Criterion %d = %15.8e" % (total_iter, a1)
                log.add(msg)

                ref_data[2] = tavg
                ref_data[3] = frsc

                #  call user-supplied function to prepare reference image, i.e., center and filter it
                if center == -1:
                    # When center = -1, which is by default, we use the average center method
                    ref_data[1] = 0
                    tavg, cs = user_func(ref_data)
                    cs[0] = old_div(float(sx_sum), nima)
                    cs[1] = old_div(float(sy_sum), nima)
                    tavg = sp_fundamentals.fshift(tavg, -cs[0], -cs[1])
                    msg = (
                        "Average center x =      %10.3f        Center y       = %10.3f\n"
                        % (cs[0], cs[1])
                    )
                    sp_utilities.print_msg(msg)
                else:
                    tavg, cs = user_func(ref_data)

                # write the current filtered average
                if outdir:
                    tavg.write_image(os.path.join(outdir, "aqf.hdf"), total_iter - 1)

                if a1 < a0:
                    if auto_stop:
                        again = 0
                else:
                    a0 = a1
            else:
                tavg = sp_utilities.model_blank(nx, nx)
                cs = [0.0] * 2

            if auto_stop:
                again = mpi.mpi_bcast(
                    again, 1, mpi.MPI_INT, main_node, mpi.MPI_COMM_WORLD
                )
                if int(again[0]) == 0:
                    break

            if Fourvar:
                del vav
            sp_utilities.bcast_EMData_to_all(tavg, myid, main_node)
            cs = mpi.mpi_bcast(cs, 2, mpi.MPI_FLOAT, main_node, mpi.MPI_COMM_WORLD)
            cs = list(map(float, cs))
            if total_iter != max_iter * len(xrng):
                if CUDA:
                    old_ali_params = all_ali_params[:]
                else:
                    old_ali_params = []
                    for im in range(len(data)):
                        alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(
                            data[im]
                        )
                        old_ali_params.extend([alpha, sx, sy, mirror])

                if Iter % 4 != 0 or total_iter > max_iter * len(xrng) - 10:
                    delta = 0.0
                else:
                    delta = dst
                if CUDA:
                    all_ali_params = R.ali2d_single_iter(
                        tavg, all_ali_params, cs[0], cs[1], 1, delta
                    )
                    sx_sum = all_ali_params[-2]
                    sy_sum = all_ali_params[-1]
                    for im in range(len(data)):
                        all_ali_params[im * 4 + 3] = int(all_ali_params[im * 4 + 3])
                else:
                    sx_sum, sy_sum, nope = sp_alignment.ali2d_single_iter(
                        data,
                        numr,
                        wr,
                        cs,
                        tavg,
                        cnx,
                        cny,
                        xrng[N_step],
                        yrng[N_step],
                        step[N_step],
                        nomirror=nomirror,
                        mode=mode,
                        CTF=CTF,
                        delta=delta,
                        random_method=random_method,
                    )

                sx_sum = mpi.mpi_reduce(
                    sx_sum, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, main_node, mpi.MPI_COMM_WORLD
                )
                sy_sum = mpi.mpi_reduce(
                    sy_sum, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, main_node, mpi.MPI_COMM_WORLD
                )
                #  for SHC
                if random_method == "SHC":
                    nope = mpi.mpi_reduce(
                        nope, 1, mpi.MPI_INT, mpi.MPI_SUM, main_node, mpi.MPI_COMM_WORLD
                    )
                    nope = mpi.mpi_bcast(
                        nope, 1, mpi.MPI_INT, main_node, mpi.MPI_COMM_WORLD
                    )
                    if int(nope[0]) == nima:
                        break

                pixel_error = 0.0
                mirror_consistent = 0
                pixel_error_list = []
                for im in range(len(data)):
                    if CUDA:
                        alpha = all_ali_params[im * 4]
                        sx = all_ali_params[im * 4 + 1]
                        sy = all_ali_params[im * 4 + 2]
                        mirror = all_ali_params[im * 4 + 3]
                    else:
                        alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(
                            data[im]
                        )
                        if old_ali_params[im * 4 + 3] == mirror:
                            this_error = sp_pixel_error.pixel_error_2D(
                                old_ali_params[im * 4 : im * 4 + 3],
                                [alpha, sx, sy],
                                last_ring,
                            )
                            pixel_error += this_error
                            pixel_error_list.append(this_error)
                            mirror_consistent += 1
                        else:
                            pixel_error_list.append(-1)
                mirror_consistent = mpi.mpi_reduce(
                    mirror_consistent,
                    1,
                    mpi.MPI_INT,
                    mpi.MPI_SUM,
                    main_node,
                    mpi.MPI_COMM_WORLD,
                )
                pixel_error = mpi.mpi_reduce(
                    pixel_error,
                    1,
                    mpi.MPI_FLOAT,
                    mpi.MPI_SUM,
                    main_node,
                    mpi.MPI_COMM_WORLD,
                )
                pixel_error_list = mpi.mpi_gatherv(
                    pixel_error_list,
                    len(data),
                    mpi.MPI_FLOAT,
                    recvcount,
                    disp,
                    mpi.MPI_FLOAT,
                    main_node,
                    mpi.MPI_COMM_WORLD,
                )
                if myid == main_node:
                    sp_utilities.print_msg(
                        "Mirror consistency rate = %8.4f%%\n"
                        % (old_div(float(mirror_consistent), nima) * 100)
                    )
                    if mirror_consistent != 0:
                        sp_utilities.print_msg(
                            "Among the mirror-consistent images, average of pixel errors is %0.4f, and their distribution is:\n"
                            % (old_div(float(pixel_error), float(mirror_consistent)))
                        )
                        pixel_error_list = list(map(float, pixel_error_list))
                        for i in range(nima - 1, -1, -1):
                            if pixel_error_list[i] < 0:
                                del pixel_error_list[i]
                        region, hist = sp_statistics.hist_list(pixel_error_list, 20)
                        for p in range(20):
                            sp_utilities.print_msg(
                                "      %10.6f: %5d\n" % (region[p], hist[p])
                            )
                    sp_utilities.print_msg("\n\n\n")
        if CUDA:
            R.finish()

    if CUDA:
        for im in range(len(data)):
            sp_utilities.set_params2D(
                data[im],
                [
                    all_ali_params[im * 4],
                    all_ali_params[im * 4 + 1],
                    all_ali_params[im * 4 + 2],
                    all_ali_params[im * 4 + 3],
                    1.0,
                ],
            )

    if myid == main_node and outdir:
        sp_utilities.drop_image(tavg, os.path.join(outdir, "aqfinal.hdf"))
    # write out headers and STOP, under MPI writing has to be done sequentially
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    par_str = ["xform.align2d", "ID"]
    if myid == main_node:
        if sp_utilities.file_type(stack) == "bdb":
            sp_utilities.recv_attr_dict_bdb(
                main_node, stack, data, par_str, image_start, image_end, number_of_proc
            )
        else:
            sp_utilities.recv_attr_dict(
                main_node, stack, data, par_str, image_start, image_end, number_of_proc
            )
    else:
        sp_utilities.send_attr_dict(main_node, data, par_str, image_start, image_end)
    if myid == main_node:
        sp_utilities.print_end_msg("ali2d_MPI")


def ali2d_base(
    stack,
    outdir,
    maskfile=None,
    ir=1,
    ou=-1,
    rs=1,
    xr="4 2 1 1",
    yr="-1",
    ts="2 1 0.5 0.25",
    nomirror=False,
    dst=0.0,
    center=-1,
    maxit=0,
    CTF=False,
    snr=1.0,
    Fourvar=False,
    user_func_name="ref_ali2d",
    random_method="",
    log=None,
    number_of_proc=1,
    myid=0,
    main_node=0,
    mpi_comm=None,
    write_headers=False,
):

    if log == None:
        log = sp_logger.Logger()

    if mpi_comm == None:
        mpi_comm = mpi.MPI_COMM_WORLD

    # ftp = file_type(stack)

    if myid == main_node:
        sp_global_def.LOGFILE = os.path.join(outdir, sp_global_def.LOGFILE)
        log.add("Start  ali2d_MPI")

    xrng = sp_utilities.get_input_from_string(xr)
    if yr == "-1":
        yrng = xrng
    else:
        yrng = sp_utilities.get_input_from_string(yr)
    step = sp_utilities.get_input_from_string(ts)

    first_ring = int(ir)
    last_ring = int(ou)
    rstep = int(rs)
    max_iter = int(maxit)

    if max_iter == 0:
        max_iter = 10
        auto_stop = True
    else:
        auto_stop = False

    if isinstance(stack, (bytes, str)):
        if myid == main_node:
            total_nima = EMAN2_cppwrap.EMUtil.get_image_count(stack)
        else:
            total_nima = 0
        total_nima = sp_utilities.bcast_number_to_all(total_nima)
        list_of_particles = list(range(total_nima))

        image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
        list_of_particles = list_of_particles[image_start:image_end]
        nima = len(list_of_particles)
        data = EMAN2_cppwrap.EMData.read_images(stack, list_of_particles)

    else:
        data = stack
        total_nima = len(data)
        total_nima = mpi.mpi_reduce(
            total_nima, 1, mpi.MPI_INT, mpi.MPI_SUM, 0, mpi.MPI_COMM_WORLD
        )
        total_nima = mpi.mpi_bcast(
            total_nima, 1, mpi.MPI_INT, main_node, mpi.MPI_COMM_WORLD
        )[0]
        list_of_particles = list(range(total_nima))
        image_start, image_end = MPI_start_end(total_nima, number_of_proc, myid)
        list_of_particles = list_of_particles[image_start:image_end]
        nima = len(list_of_particles)
        assert nima == len(data)

    # read nx and ctf_app (if CTF) and broadcast to all nodes
    if myid == main_node:
        nx = data[0].get_xsize()
        if CTF:
            ctf_app = data[0].get_attr_default("ctf_applied", 0)
        # del ima
    else:
        nx = 0
        if CTF:
            ctf_app = 0
    nx = sp_utilities.bcast_number_to_all(nx, source_node=main_node)
    if CTF:
        ctf_app = sp_utilities.bcast_number_to_all(ctf_app, source_node=main_node)
        if ctf_app > 0:
            sp_global_def.ERROR("data cannot be ctf-applied", "ali2d_MPI", 1, myid)
        phase_flip = True
    else:
        phase_flip = False
    CTF = False

    # default value for the last ring
    if last_ring == -1:
        last_ring = old_div(nx, 2) - 2

    if last_ring + max([max(xrng), max(yrng)]) > old_div((nx - 1), 2):
        sp_global_def.ERROR(
            "Shift or radius is too large - particle crosses image boundary",
            "ali2d_MPI",
            1,
        )

    if myid == main_node:
        # log.add("Input stack                 : %s"%(stack))
        log.add("Number of images            : %d" % (total_nima))
        log.add("Output directory            : %s" % (outdir))
        log.add("Inner radius                : %i" % (first_ring))
        log.add("Outer radius                : %i" % (last_ring))
        log.add("Ring step                   : %i" % (rstep))
        log.add("X search range              : %s" % (xrng))
        log.add("Y search range              : %s" % (yrng))
        log.add("Translational step          : %s" % (step))
        log.add("Disable checking mirror     : %s" % (nomirror))
        log.add("Discrete angle used         : %d" % (dst))
        log.add("Center type                 : %i" % (center))
        log.add("Maximum iteration           : %i" % (max_iter))
        # log.add("Use Fourier variance        : %s\n"%(Fourvar))
        log.add("CTF correction              : %s" % (CTF))
        log.add("Phase flip                  : %s" % (phase_flip))
        # log.add("Signal-to-Noise Ratio       : %f\n"%(snr))
        if auto_stop:
            log.add("Stop iteration with         : criterion")
        else:
            log.add("Stop iteration with         : maxit")

        user_func = sp_user_functions.factory[user_func_name]

        log.add("User function               : %s" % (user_func_name))
        log.add("Number of processors used   : %d" % (number_of_proc))

    if maskfile:
        if isinstance(maskfile, (bytes, str)):
            if myid == main_node:
                log.add("Maskfile                    : %s" % (maskfile))
            mask = sp_utilities.get_image(maskfile)
        else:
            if myid == main_node:
                log.add("Maskfile                    : user provided in-core mask")
            mask = maskfile
    else:
        if myid == main_node:
            log.add(
                "Maskfile                    : default, a circle with radius %i"
                % (last_ring)
            )
        mask = sp_utilities.model_circle(last_ring, nx, nx)

    cnx = old_div(nx, 2) + 1
    cny = cnx
    if random_method == "SCF":
        mode = "H"
    else:
        mode = "F"

    if CTF:
        ctf_abs_sum = EMAN2_cppwrap.EMData(nx, nx, 1, False)
        ctf_2_sum = EMAN2_cppwrap.EMData(nx, nx, 1, False)
    else:
        ctf_2_sum = None

    for im in range(nima):
        data[im].set_attr("ID", list_of_particles[im])
        sp_utilities.set_params2D(data[im], [0.0, 0.0, 0.0, 0, 1.0], "xform.align2d")
        st = EMAN2_cppwrap.Util.infomask(data[im], mask, False)
        data[im] -= st[0]
        if CTF:
            ctf_params = data[im].get_attr("ctf")
            ctfimg = sp_morphology.ctf_img(nx, ctf_params)
            EMAN2_cppwrap.Util.add_img2(ctf_2_sum, ctfimg)
            EMAN2_cppwrap.Util.add_img_abs(ctf_abs_sum, ctfimg)
        if random_method == "SHC":
            data[im].set_attr("previousmax", 1.0e-23)
        if phase_flip:
            data[im] = sp_filter.filt_ctf(
                data[im], data[im].get_attr("ctf"), binary=True
            )

    if CTF:
        sp_utilities.reduce_EMData_to_root(ctf_2_sum, myid, main_node)
        sp_utilities.reduce_EMData_to_root(ctf_abs_sum, myid, main_node)
        if myid == main_node:
            adw_img = EMAN2_cppwrap.Util.mult_scalar(ctf_2_sum, snr)
            EMAN2_cppwrap.Util.div_filter(adw_img, ctf_abs_sum)
            EMAN2_cppwrap.Util.mul_scalar(adw_img, old_div(float(Ng - 1), (nima - 1)))
            adw_img += old_div(float(nima - Ng), (nima - 1))
    else:
        ctf_2_sum = None

    # startup
    numr = sp_alignment.Numrinit(
        first_ring, last_ring, rstep, mode
    )  # precalculate rings
    wr = sp_alignment.ringwe(numr, mode)

    if myid == main_node:
        # initialize data for the reference preparation function
        ref_data = [mask, center, None, None]
        sx_sum = 0.0
        sy_sum = 0.0
        a0 = -1.0e22

    recvcount = []
    disp = []
    for i in range(number_of_proc):
        ib, ie = MPI_start_end(total_nima, number_of_proc, i)
        recvcount.append(ie - ib)
        if i == 0:
            disp.append(0)
        else:
            disp.append(disp[i - 1] + recvcount[i - 1])

    again = 1
    total_iter = 0
    cs = [0.0] * 2
    delta = 0.0
    for N_step in range(len(xrng)):

        for Iter in range(max_iter):
            total_iter += 1
            ave1, ave2 = sp_statistics.sum_oe(
                data, "a", CTF, EMAN2_cppwrap.EMData()
            )  # pass empty object to prevent calculation of ctf^2
            sp_utilities.reduce_EMData_to_root(ave1, myid, main_node)
            sp_utilities.reduce_EMData_to_root(ave2, myid, main_node)
            sys.stdout.flush()
            if myid == main_node:
                log.add("Iteration #%4d" % (total_iter))
                msg = "X range = %5.2f   Y range = %5.2f   Step = %5.2f" % (
                    xrng[N_step],
                    yrng[N_step],
                    step[N_step],
                )
                log.add(msg)
                if CTF:
                    tavg_Ng = sp_fundamentals.fft(
                        EMAN2_cppwrap.Util.divn_filter(
                            EMAN2_cppwrap.Util.muln_img(
                                sp_fundamentals.fft(
                                    EMAN2_cppwrap.Util.addn_img(ave1, ave2)
                                ),
                                adw_img,
                            ),
                            ctf_2_sum,
                        )
                    )
                    tavg = sp_fundamentals.fft(
                        EMAN2_cppwrap.Util.divn_filter(
                            sp_fundamentals.fft(
                                EMAN2_cppwrap.Util.addn_img(ave1, ave2)
                            ),
                            ctf_2_sum,
                        )
                    )
                else:
                    tavg = old_div((ave1 + ave2), total_nima)
                if outdir:
                    tavg.write_image(os.path.join(outdir, "aqc.hdf"), total_iter - 1)

                    if CTF:
                        tavg_Ng.write_image(
                            os.path.join(outdir, "aqc_view.hdf"), total_iter - 1
                        )
                    frsc = sp_statistics.fsc_mask(
                        ave1,
                        ave2,
                        mask,
                        1.0,
                        os.path.join(outdir, "resolution%03d" % (total_iter)),
                    )
                else:
                    frsc = sp_statistics.fsc_mask(ave1, ave2, mask, 1.0)
            else:
                tavg = sp_utilities.model_blank(nx, nx)
            del ave1, ave2

            if Fourvar:
                sp_utilities.bcast_EMData_to_all(tavg, myid, main_node)
                vav, rvar = sp_statistics.varf2d_MPI(myid, data, tavg, mask, "a", CTF)

            if myid == main_node:
                if Fourvar:
                    tavg = sp_fundamentals.fft(
                        EMAN2_cppwrap.Util.divn_img(sp_fundamentals.fft(tavg), vav)
                    )
                    vav_r = EMAN2_cppwrap.Util.pack_complex_to_real(vav)
                    if outdir:
                        vav_r.write_image(
                            os.path.join(outdir, "varf.hdf"), total_iter - 1
                        )

                # a0 should increase; stop algorithm when it decreases.
                #     However, the result will depend on filtration, so it is not quite right.
                #  moved it here, so it is for unfiltered average and thus hopefully makes more sense
                a1 = tavg.cmp("dot", tavg, dict(negative=0, mask=ref_data[0]))
                msg = "Criterion %d = %15.8e" % (total_iter, a1)
                log.add(msg)

                ref_data[2] = tavg
                ref_data[3] = frsc

                #  call user-supplied function to prepare reference image, i.e., center and filter it
                if center == -1:
                    # When center = -1, which is by default, we use the average center method
                    ref_data[1] = 0
                    tavg, cs = user_func(ref_data)
                    cs[0] = old_div(float(sx_sum), total_nima)
                    cs[1] = old_div(float(sy_sum), total_nima)
                    tavg = sp_fundamentals.fshift(tavg, -cs[0], -cs[1])
                    msg = (
                        "Average center x =      %10.3f        Center y       = %10.3f"
                        % (cs[0], cs[1])
                    )
                    log.add(msg)
                else:
                    if delta != 0.0:
                        cnt = ref_data[1]
                        ref_data[1] = 0
                    tavg, cs = user_func(ref_data)
                    if delta != 0.0:
                        ref_data[1] = cnt
                # write the current filtered average
                if outdir:
                    tavg.write_image(os.path.join(outdir, "aqf.hdf"), total_iter - 1)

                if a1 < a0:
                    if auto_stop:
                        again = 0
                else:
                    a0 = a1
            else:
                tavg = sp_utilities.model_blank(nx, nx)
                cs = [0.0] * 2

            if auto_stop:
                again = mpi.mpi_bcast(again, 1, mpi.MPI_INT, main_node, mpi_comm)
                if int(again[0]) == 0:
                    break

            if Fourvar:
                del vav
            sp_utilities.bcast_EMData_to_all(tavg, myid, main_node)
            cs = mpi.mpi_bcast(cs, 2, mpi.MPI_FLOAT, main_node, mpi_comm)
            cs = list(map(float, cs))
            if total_iter != max_iter * len(xrng):
                old_ali_params = []
                for im in range(nima):
                    alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(data[im])
                    old_ali_params.extend([alpha, sx, sy, mirror])

                if Iter % 4 != 0 or total_iter > max_iter * len(xrng) - 10:
                    delta = 0.0
                else:
                    delta = dst
                sx_sum, sy_sum, nope = sp_alignment.ali2d_single_iter(
                    data,
                    numr,
                    wr,
                    cs,
                    tavg,
                    cnx,
                    cny,
                    xrng[N_step],
                    yrng[N_step],
                    step[N_step],
                    nomirror=nomirror,
                    mode=mode,
                    CTF=CTF,
                    delta=delta,
                    random_method=random_method,
                )

                sx_sum = mpi.mpi_reduce(
                    sx_sum, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, main_node, mpi_comm
                )
                sy_sum = mpi.mpi_reduce(
                    sy_sum, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, main_node, mpi_comm
                )
                #  for SHC
                if random_method == "SHC":
                    nope = mpi.mpi_reduce(
                        nope, 1, mpi.MPI_INT, mpi.MPI_SUM, main_node, mpi_comm
                    )
                    nope = mpi.mpi_bcast(nope, 1, mpi.MPI_INT, main_node, mpi_comm)
                    if int(nope[0]) == total_nima:
                        break

                pixel_error = 0.0
                mirror_consistent = 0
                pixel_error_list = [-1.0] * nima
                for im in range(nima):
                    alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(data[im])
                    if old_ali_params[im * 4 + 3] == mirror:
                        this_error = sp_pixel_error.pixel_error_2D(
                            old_ali_params[im * 4 : im * 4 + 3],
                            [alpha, sx, sy],
                            last_ring,
                        )
                        pixel_error += this_error
                        pixel_error_list[im] = this_error
                        mirror_consistent += 1
                del old_ali_params
                mirror_consistent = mpi.mpi_reduce(
                    mirror_consistent, 1, mpi.MPI_INT, mpi.MPI_SUM, main_node, mpi_comm
                )
                pixel_error = mpi.mpi_reduce(
                    pixel_error, 1, mpi.MPI_FLOAT, mpi.MPI_SUM, main_node, mpi_comm
                )
                pixel_error_list = mpi.mpi_gatherv(
                    pixel_error_list,
                    nima,
                    mpi.MPI_FLOAT,
                    recvcount,
                    disp,
                    mpi.MPI_FLOAT,
                    main_node,
                    mpi_comm,
                )
                if myid == main_node:
                    log.add(
                        "Mirror consistency rate = %8.4f%%"
                        % (old_div(float(mirror_consistent), total_nima) * 100)
                    )
                    if mirror_consistent != 0:
                        log.add(
                            "Among the mirror-consistent images, average of pixel errors is %0.4f, and their distribution is:"
                            % (old_div(float(pixel_error), float(mirror_consistent)))
                        )
                        pixel_error_list = list(map(float, pixel_error_list))
                        for i in range(total_nima - 1, -1, -1):
                            if pixel_error_list[i] < 0:
                                del pixel_error_list[i]
                        region, hist = sp_statistics.hist_list(pixel_error_list, 20)
                        for p in range(20):
                            log.add("      %14.2f: %6d" % (region[p], hist[p]))
                    log.add("\n\n")

    if myid == main_node and outdir:
        sp_utilities.drop_image(tavg, os.path.join(outdir, "aqfinal.hdf"))
    # write out headers and STOP, under MPI writing has to be done sequentially
    mpi.mpi_barrier(mpi_comm)
    if write_headers:
        par_str = ["xform.align2d", "ID"]
        if myid == main_node:
            if sp_utilities.file_type(stack) == "bdb":
                sp_utilities.recv_attr_dict_bdb(
                    main_node,
                    stack,
                    data,
                    par_str,
                    image_start,
                    image_end,
                    number_of_proc,
                )
            else:
                sp_utilities.recv_attr_dict(
                    main_node,
                    stack,
                    data,
                    par_str,
                    image_start,
                    image_end,
                    number_of_proc,
                )
        else:
            sp_utilities.send_attr_dict(
                main_node, data, par_str, image_start, image_end
            )
    params = []
    for im in range(nima):
        alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(data[im])
        params.append([alpha, sx, sy, mirror])
    # params = wrap_mpi_gatherv(params, main_node, mpi_comm)

    if myid == main_node:
        log.add("Finished ali2d_base")

    return params  # , data


"""Multiline Comment1"""


def mref_ali2d(
    stack,
    refim,
    outdir,
    maskfile=None,
    ir=1,
    ou=-1,
    rs=1,
    xrng=0,
    yrng=0,
    step=1,
    center=1,
    maxit=0,
    CTF=False,
    snr=1.0,
    user_func_name="ref_ali2d",
    rand_seed=1000,
    MPI=False,
):
    """
		Name
			mref_ali2d - Perform 2-D multi-reference alignment of an image series
		Input
			stack: set of 2-D images in a stack file, images have to be squares
			refim: set of initial reference 2-D images in a stack file
			maskfile: optional maskfile to be used in the alignment
			inner_radius: inner radius for rotational correlation > 0
			outer_radius: outer radius for rotational correlation < nx/2-1
			ring_step: step between rings in rotational correlation >0
			x_range: range for translation search in x direction, search is +/xr
			y_range: range for translation search in y direction, search is +/yr
			translation_step: step of translation search in both directions
			center: center the average
			max_iter: maximum number of iterations the program will perform
			CTF: if this flag is set, the program will use CTF information provided in file headers
			snr: signal-to-noise ratio of the data
			rand_seed: the seed used for generating random numbers
			MPI: whether to use MPI version
		Output
			output_directory: directory name into which the output files will be written.
			header: the alignment parameters are stored in the headers of input files as 'xform.align2d'.
	"""
    # 2D multi-reference alignment using rotational ccf in polar coordinates and quadratic interpolation
    if MPI:
        mref_ali2d_MPI(
            stack,
            refim,
            outdir,
            maskfile,
            ir,
            ou,
            rs,
            xrng,
            yrng,
            step,
            center,
            maxit,
            CTF,
            snr,
            user_func_name,
            rand_seed,
        )
        return

    # create the output directory, if it does not exist
    if os.path.exists(outdir):
        sp_global_def.ERROR(
            "Output directory exists, please change the name and restart the program",
            "mref_ali2d",
            1,
        )
    os.mkdir(outdir)
    sp_global_def.LOGFILE = os.path.join(outdir, sp_global_def.LOGFILE)

    first_ring = int(ir)
    last_ring = int(ou)
    rstep = int(rs)
    max_iter = int(maxit)
    if max_iter == 0:
        max_iter = 10
        auto_stop = True
    else:
        auto_stop = False

    sp_utilities.print_begin_msg("mref_ali2d")

    sp_utilities.print_msg("Input stack                 : %s\n" % (stack))
    sp_utilities.print_msg("Reference stack             : %s\n" % (refim))
    sp_utilities.print_msg("Output directory            : %s\n" % (outdir))
    sp_utilities.print_msg("Maskfile                    : %s\n" % (maskfile))
    sp_utilities.print_msg("Inner radius                : %i\n" % (first_ring))

    ima = EMAN2_cppwrap.EMData()
    ima.read_image(stack, 0)
    nx = ima.get_xsize()
    # default value for the last ring
    if last_ring == -1:
        last_ring = old_div(nx, 2) - 2

    sp_utilities.print_msg("Outer radius                : %i\n" % (last_ring))
    sp_utilities.print_msg("Ring step                   : %i\n" % (rstep))
    sp_utilities.print_msg("X search range              : %i\n" % (xrng))
    sp_utilities.print_msg("Y search range              : %i\n" % (yrng))
    sp_utilities.print_msg("Translational step          : %i\n" % (step))
    sp_utilities.print_msg("Center type                 : %i\n" % (center))
    sp_utilities.print_msg("Maximum iteration           : %i\n" % (max_iter))
    sp_utilities.print_msg("CTF correction              : %s\n" % (CTF))
    sp_utilities.print_msg("Signal-to-Noise Ratio       : %f\n" % (snr))
    sp_utilities.print_msg("Random seed                 : %i\n\n" % (rand_seed))
    sp_utilities.print_msg("User function               : %s\n" % (user_func_name))
    output = sys.stdout

    user_func = sp_user_functions.factory[user_func_name]

    if maskfile:
        if isinstance(maskfile, (bytes, str)):
            mask = sp_utilities.get_image(maskfile)
        else:
            mask = maskfile
    else:
        mask = sp_utilities.model_circle(last_ring, nx, nx)
    #  references
    refi = []
    numref = EMAN2_cppwrap.EMUtil.get_image_count(refim)
    #  CTF stuff
    if CTF:
        ctf_params = ima.get_attr("ctf")
        data_had_ctf = ima.get_attr("ctf_applied")
        ctm = sp_morphology.ctf_2(nx, ctf_params)
        lctf = len(ctm)
        ctf2 = [[[0.0] * lctf for k in range(2)] for j in range(numref)]

    # IMAGES ARE SQUARES! center is in SPIDER convention
    cnx = old_div(nx, 2) + 1
    cny = cnx

    mode = "F"
    # precalculate rings
    numr = sp_alignment.Numrinit(first_ring, last_ring, rstep, mode)
    wr = sp_alignment.ringwe(numr, mode)
    # reference images
    params = []
    # read all data
    data = EMAN2_cppwrap.EMData.read_images(stack)
    nima = len(data)
    # prepare the reference
    ima.to_zero()
    for j in range(numref):
        temp = EMAN2_cppwrap.EMData()
        temp.read_image(refim, j)
        #  eve, odd, numer of even, number of images.  After frc, totav
        refi.append([temp, ima.copy(), 0])
    random.seed(rand_seed)
    a0 = -1.0
    again = True
    Iter = 0

    ref_data = [mask, center, None, None]

    while Iter < max_iter and again:
        # again = False
        ringref = []
        # print "numref",numref
        mashi = cnx - last_ring - 2
        for j in range(numref):
            refi[j][0].process_inplace("normalize.mask", {"mask": mask, "no_sigma": 1})
            cimage = EMAN2_cppwrap.Util.Polar2Dm(refi[j][0], cnx, cny, numr, mode)
            EMAN2_cppwrap.Util.Frngs(cimage, numr)
            EMAN2_cppwrap.Util.Applyws(cimage, numr, wr)
            ringref.append(cimage)
            # zero refi
            refi[j][0].to_zero()
            refi[j][1].to_zero()
            refi[j][2] = 0
            if CTF:
                for i in range(lctf):
                    ctf2[j][0][i] = 0.0
                    ctf2[j][1][i] = 0.0
        assign = [[] for i in range(numref)]
        sx_sum = [0.0] * numref
        sy_sum = [0.0] * numref
        for im in range(nima):
            if CTF:
                ctf_params = data[im].get_attr("ctf")
                if data[im].get_attr("ctf_applied") == 0:
                    st = EMAN2_cppwrap.Util.infomask(data[im], mask, False)
                    data[im] -= st[0]
                    data[im] = sp_filter.filt_ctf(data[im], ctf_params)
                    data[im].set_attr("ctf_applied", 1)
            alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(data[im])
            #  Why inverse?  07/11/2015  PAP
            alphai, sxi, syi, scalei = sp_utilities.inverse_transform2(alpha, sx, sy)
            # normalize
            data[im].process_inplace("normalize.mask", {"mask": mask, "no_sigma": 0})
            # If shifts are outside of the permissible range, reset them
            if abs(sxi) > mashi or abs(syi) > mashi:
                sxi = 0.0
                syi = 0.0
                sp_utilities.set_params2D(data[im], [0.0, 0.0, 0.0, 0, 1.0])
            ny = nx
            txrng = sp_alignment.search_range(nx, last_ring, sxi, xrng, "mref_ali2d")
            txrng = [txrng[1], txrng[0]]
            tyrng = sp_alignment.search_range(ny, last_ring, syi, yrng, "mref_ali2d")
            tyrng = [tyrng[1], tyrng[0]]
            # align current image to the reference
            [
                angt,
                sxst,
                syst,
                mirrort,
                xiref,
                peakt,
            ] = EMAN2_cppwrap.Util.multiref_polar_ali_2d(
                data[im], ringref, txrng, tyrng, step, mode, numr, cnx + sxi, cny + syi
            )
            iref = int(xiref)
            # combine parameters and set them to the header, ignore previous angle and mirror
            [alphan, sxn, syn, mn] = sp_utilities.combine_params2(
                0.0, -sxi, -syi, 0, angt, sxst, syst, int(mirrort)
            )
            sp_utilities.set_params2D(data[im], [alphan, sxn, syn, int(mn), scale])
            if mn == 0:
                sx_sum[iref] += sxn
            else:
                sx_sum[iref] -= sxn
            sy_sum[iref] += syn
            data[im].set_attr("assign", iref)
            # apply current parameters and add to the average
            temp = sp_fundamentals.rot_shift2D(data[im], alphan, sxn, syn, mn)
            it = im % 2
            EMAN2_cppwrap.Util.add_img(refi[iref][it], temp)
            if CTF:
                ctm = sp_morphology.ctf_2(nx, ctf_params)
                for i in range(lctf):
                    ctf2[iref][it][i] += ctm[i]
            assign[iref].append(im)
            refi[iref][2] += 1
        del ringref
        if again:
            a1 = 0.0
            for j in range(numref):
                msg = "   group #%3d   number of particles = %7d\n" % (j, refi[j][2])
                sp_utilities.print_msg(msg)
                if refi[j][2] < 4:
                    # ERROR("One of the references vanished","mref_ali2d",1)
                    #  if vanished, put a random image there
                    assign[j] = []
                    assign[j].append(random.randint(0, nima - 1))
                    refi[j][0] = data[assign[j][0]].copy()
                else:
                    max_inter = 0  # switch off fine refi.
                    br = 1.75
                    #  the loop has to
                    for INter in range(max_inter + 1):
                        # Calculate averages at least ones, meaning even if no within group refinement was requested
                        if CTF:
                            for i in range(lctf):
                                ctm[i] = old_div(
                                    1.0, (ctf2[j][0][i] + old_div(1.0, snr))
                                )
                            av1 = sp_filter.filt_table(refi[j][0], ctm)
                            for i in range(lctf):
                                ctm[i] = old_div(
                                    1.0, (ctf2[j][1][i] + old_div(1.0, snr))
                                )
                            av2 = sp_filter.filt_table(refi[j][1], ctm)
                            frsc = sp_statistics.fsc(
                                av1,
                                av2,
                                1.0,
                                os.path.join(outdir, "drm_%03d_%04d.txt" % (Iter, j)),
                            )
                            # Now the total average
                            for i in range(lctf):
                                ctm[i] = old_div(
                                    1.0,
                                    (ctf2[j][0][i] + ctf2[j][1][i] + old_div(1.0, snr)),
                                )
                            refi[j][0] = sp_filter.filt_table(
                                EMAN2_cppwrap.Util.addn_img(refi[j][0], refi[j][1]), ctm
                            )
                        else:
                            frsc = sp_statistics.fsc(
                                refi[j][0],
                                refi[j][1],
                                1.0,
                                os.path.join(outdir, "drm_%03d_%04d.txt" % (Iter, j)),
                            )
                            EMAN2_cppwrap.Util.add_img(refi[j][0], refi[j][1])
                            EMAN2_cppwrap.Util.mul_scalar(
                                refi[j][0], old_div(1.0, float(refi[j][2]))
                            )

                        ref_data[2] = refi[j][0]
                        ref_data[3] = frsc
                        refi[j][0], cs = user_func(ref_data)
                        if center == -1:
                            cs[0] = old_div(sx_sum[j], len(assign[j]))
                            cs[1] = old_div(sy_sum[j], len(assign[j]))
                            refi[j][0] = sp_fundamentals.fshift(
                                refi[j][0], -cs[0], -cs[1]
                            )
                        for i in range(len(assign[j])):
                            im = assign[j][i]
                            alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(
                                data[im]
                            )
                            alphan, sxn, syn, mirrorn = sp_utilities.combine_params2(
                                alpha, sx, sy, mirror, 0.0, -cs[0], -cs[1], 0
                            )
                            sp_utilities.set_params2D(
                                data[im], [alphan, sxn, syn, int(mirrorn), scale]
                            )
                        # refine images within the group
                        #  Do the refinement only if max_inter>0, but skip it for the last iteration.
                        if INter < max_inter:
                            sp_alignment.fine_2D_refinement(
                                data, br, mask, refi[j][0], j
                            )
                            #  Calculate updated average
                            refi[j][0].to_zero()
                            refi[j][1].to_zero()
                            for i in range(len(assign[j])):
                                im = assign[j][i]
                                alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(
                                    data[im]
                                )
                                # apply current parameters and add to the average
                                temp = sp_fundamentals.rot_shift2D(
                                    data[im], alpha, sx, sy, mn
                                )
                                it = im % 2
                                EMAN2_cppwrap.Util.add_img(refi[j][it], temp)

                # write the current average
                TMP = []
                for i_tmp in range(len(assign[j])):
                    TMP.append(float(assign[j][i_tmp]))
                TMP.sort()
                refi[j][0].set_attr_dict({"ave_n": refi[j][2], "members": TMP})
                del TMP
                # replace the name of the stack with reference with the current one
                newrefim = os.path.join(outdir, "aqm%03d.hdf" % Iter)
                refi[j][0].write_image(newrefim, j)
                a1 += refi[j][0].cmp("dot", refi[j][0], {"negative": 0, "mask": mask})
            Iter += 1
            msg = "ITERATION #%3d        criterion = %20.7e\n" % (Iter, a1)
            sp_utilities.print_msg(msg)
            if a1 < a0:
                if auto_stop == True:
                    break
            else:
                a0 = a1

    newrefim = os.path.join(outdir, "multi_ref.hdf")
    for j in range(numref):
        refi[j][0].write_image(newrefim, j)
    if CTF:
        if data_had_ctf == 0:
            for im in range(nima):
                data[im].set_attr("ctf_applied", 0)
    sp_utilities.write_headers(stack, data, list(range(nima)))
    sp_utilities.print_end_msg("mref_ali2d")


def mref_ali2d_MPI(
    stack,
    refim,
    outdir,
    maskfile=None,
    ir=1,
    ou=-1,
    rs=1,
    xrng=0,
    yrng=0,
    step=1,
    center=1,
    maxit=10,
    CTF=False,
    snr=1.0,
    user_func_name="ref_ali2d",
    rand_seed=1000,
):
    # 2D multi-reference alignment using rotational ccf in polar coordinates and quadratic interpolation

    number_of_proc = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
    myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
    main_node = 0

    # create the output directory, if it does not exist
    if os.path.exists(outdir):
        sp_global_def.ERROR(
            "Output directory exists, please change the name and restart the program",
            "mref_ali2d_MPI ",
            1,
            myid,
        )
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    if myid == main_node:
        os.mkdir(outdir)
        sp_global_def.LOGFILE = os.path.join(outdir, sp_global_def.LOGFILE)
        sp_utilities.print_begin_msg("mref_ali2d_MPI")

    nima = EMAN2_cppwrap.EMUtil.get_image_count(stack)

    image_start, image_end = MPI_start_end(nima, number_of_proc, myid)

    nima = EMAN2_cppwrap.EMUtil.get_image_count(stack)
    ima = EMAN2_cppwrap.EMData()
    ima.read_image(stack, image_start)

    first_ring = int(ir)
    last_ring = int(ou)
    rstep = int(rs)
    max_iter = int(maxit)

    if max_iter == 0:
        max_iter = 10
        auto_stop = True
    else:
        auto_stop = False

    if myid == main_node:
        sp_utilities.print_msg("Input stack                 : %s\n" % (stack))
        sp_utilities.print_msg("Reference stack             : %s\n" % (refim))
        sp_utilities.print_msg("Output directory            : %s\n" % (outdir))
        sp_utilities.print_msg("Maskfile                    : %s\n" % (maskfile))
        sp_utilities.print_msg("Inner radius                : %i\n" % (first_ring))

    nx = ima.get_xsize()
    # default value for the last ring
    if last_ring == -1:
        last_ring = old_div(nx, 2) - 2

    if myid == main_node:
        sp_utilities.print_msg("Outer radius                : %i\n" % (last_ring))
        sp_utilities.print_msg("Ring step                   : %i\n" % (rstep))
        sp_utilities.print_msg("X search range              : %f\n" % (xrng))
        sp_utilities.print_msg("Y search range              : %f\n" % (yrng))
        sp_utilities.print_msg("Translational step          : %f\n" % (step))
        sp_utilities.print_msg("Center type                 : %i\n" % (center))
        sp_utilities.print_msg("Maximum iteration           : %i\n" % (max_iter))
        sp_utilities.print_msg("CTF correction              : %s\n" % (CTF))
        sp_utilities.print_msg("Signal-to-Noise Ratio       : %f\n" % (snr))
        sp_utilities.print_msg("Random seed                 : %i\n\n" % (rand_seed))
        sp_utilities.print_msg("User function               : %s\n" % (user_func_name))
    user_func = sp_user_functions.factory[user_func_name]

    if maskfile:
        if isinstance(maskfile, (bytes, str)):
            mask = sp_utilities.get_image(maskfile)
        else:
            mask = maskfile
    else:
        mask = sp_utilities.model_circle(last_ring, nx, nx)
    #  references, do them on all processors...
    refi = []
    numref = EMAN2_cppwrap.EMUtil.get_image_count(refim)
    #  CTF stuff
    if CTF:
        ctf_params = ima.get_attr("ctf")
        data_had_ctf = ima.get_attr("ctf_applied")
        ctm = sp_morphology.ctf_2(nx, ctf_params)
        lctf = len(ctm)

    # IMAGES ARE SQUARES! center is in SPIDER convention
    cnx = old_div(nx, 2) + 1
    cny = cnx

    mode = "F"
    # precalculate rings
    numr = sp_alignment.Numrinit(first_ring, last_ring, rstep, mode)
    wr = sp_alignment.ringwe(numr, mode)
    # reference images
    again = True
    params = []
    # prepare reference images on all nodes
    ima.to_zero()
    for j in range(numref):
        #  even, odd, numer of even, number of images.  After frc, totav
        refi.append([sp_utilities.get_im(refim, j), ima.copy(), 0])
    #  for each node read its share of data
    data = EMAN2_cppwrap.EMData.read_images(stack, list(range(image_start, image_end)))
    for im in range(image_start, image_end):
        data[im - image_start].set_attr("ID", im)
        if CTF:
            ctf_params = data[im - image_start].get_attr("ctf")
            if data[im - image_start].get_attr("ctf_applied") == 0:
                st = EMAN2_cppwrap.Util.infomask(data[im - image_start], mask, False)
                data[im - image_start] -= st[0]
                data[im - image_start] = sp_filter.filt_ctf(
                    data[im - image_start], ctf_params
                )
                data[im - image_start].set_attr("ctf_applied", 1)
    if myid == main_node:
        random.seed(rand_seed)

    a0 = -1.0
    again = True
    Iter = 0

    ref_data = [mask, center, None, None]

    while Iter < max_iter and again:
        ringref = []
        mashi = cnx - last_ring - 2
        for j in range(numref):
            refi[j][0].process_inplace(
                "normalize.mask", {"mask": mask, "no_sigma": 1}
            )  # normalize reference images to N(0,1)
            cimage = EMAN2_cppwrap.Util.Polar2Dm(refi[j][0], cnx, cny, numr, mode)
            EMAN2_cppwrap.Util.Frngs(cimage, numr)
            EMAN2_cppwrap.Util.Applyws(cimage, numr, wr)
            ringref.append(cimage)
            # zero refi
            refi[j][0].to_zero()
            refi[j][1].to_zero()
            refi[j][2] = 0
        if CTF:
            ctf2 = [[[0.0] * lctf for k in range(2)] for j in range(numref)]
        assign = [[] for i in range(numref)]
        # begin MPI section
        for im in range(image_start, image_end):
            alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(
                data[im - image_start]
            )
            #  Why inverse?  07/11/2015 PAP
            alphai, sxi, syi, scalei = sp_utilities.inverse_transform2(alpha, sx, sy)
            # normalize
            data[im - image_start].process_inplace(
                "normalize.mask", {"mask": mask, "no_sigma": 0}
            )  # subtract average under the mask
            # If shifts are outside of the permissible range, reset them
            if abs(sxi) > mashi or abs(syi) > mashi:
                sxi = 0.0
                syi = 0.0
                sp_utilities.set_params2D(
                    data[im - image_start], [0.0, 0.0, 0.0, 0, 1.0]
                )
            ny = nx
            txrng = sp_alignment.search_range(
                nx, last_ring, sxi, xrng, "mref_ali2d_MPI"
            )
            txrng = [txrng[1], txrng[0]]
            tyrng = sp_alignment.search_range(
                ny, last_ring, syi, yrng, "mref_ali2d_MPI"
            )
            tyrng = [tyrng[1], tyrng[0]]
            # align current image to the reference
            [
                angt,
                sxst,
                syst,
                mirrort,
                xiref,
                peakt,
            ] = EMAN2_cppwrap.Util.multiref_polar_ali_2d(
                data[im - image_start],
                ringref,
                txrng,
                tyrng,
                step,
                mode,
                numr,
                cnx + sxi,
                cny + syi,
            )
            iref = int(xiref)
            # combine parameters and set them to the header, ignore previous angle and mirror
            [alphan, sxn, syn, mn] = sp_utilities.combine_params2(
                0.0, -sxi, -syi, 0, angt, sxst, syst, (int)(mirrort)
            )
            sp_utilities.set_params2D(
                data[im - image_start], [alphan, sxn, syn, int(mn), scale]
            )
            data[im - image_start].set_attr("assign", iref)
            # apply current parameters and add to the average
            temp = sp_fundamentals.rot_shift2D(
                data[im - image_start], alphan, sxn, syn, mn
            )
            it = im % 2
            EMAN2_cppwrap.Util.add_img(refi[iref][it], temp)
            assign[iref].append(im)
            if CTF:
                #  I wonder whether params are still there....
                ctf_params = data[im - image_start].get_attr("ctf")
                ctm = sp_morphology.ctf_2(nx, ctf_params)
                for i in range(lctf):
                    ctf2[iref][it][i] += ctm[i]
            # assign[im] = iref
            refi[iref][2] += 1.0
        del ringref
        # end MPI section, bring partial things together, calculate new reference images, broadcast them back
        if CTF:
            # bring ctf2 together on main node
            s = numpy.shape(ctf2)
            ctf2 = mpi.mpi_reduce(
                ctf2,
                2 * lctf * numref,
                mpi.MPI_FLOAT,
                mpi.MPI_SUM,
                main_node,
                mpi.MPI_COMM_WORLD,
            )
            if myid == main_node:
                ctf2 = numpy.reshape(ctf2, s)
        for j in range(numref):
            sp_utilities.reduce_EMData_to_root(refi[j][0], myid, main_node)
            sp_utilities.reduce_EMData_to_root(refi[j][1], myid, main_node)
            refi[j][2] = mpi.mpi_reduce(
                refi[j][2], 1, mpi.MPI_FLOAT, mpi.MPI_SUM, main_node, mpi.MPI_COMM_WORLD
            )
            if myid == main_node:
                refi[j][2] = int(refi[j][2][0])
        # gather assignements
        for j in range(numref):
            if myid == main_node:
                for n in range(number_of_proc):
                    if n != main_node:
                        ln = mpi.mpi_recv(
                            1,
                            mpi.MPI_INT,
                            n,
                            sp_global_def.SPARX_MPI_TAG_UNIVERSAL,
                            mpi.MPI_COMM_WORLD,
                        )
                        lis = mpi.mpi_recv(
                            ln[0],
                            mpi.MPI_INT,
                            n,
                            sp_global_def.SPARX_MPI_TAG_UNIVERSAL,
                            mpi.MPI_COMM_WORLD,
                        )
                        for l in range(ln[0]):
                            assign[j].append(int(lis[l]))
            else:
                mpi.mpi_send(
                    len(assign[j]),
                    1,
                    mpi.MPI_INT,
                    main_node,
                    sp_global_def.SPARX_MPI_TAG_UNIVERSAL,
                    mpi.MPI_COMM_WORLD,
                )
                mpi.mpi_send(
                    assign[j],
                    len(assign[j]),
                    mpi.MPI_INT,
                    main_node,
                    sp_global_def.SPARX_MPI_TAG_UNIVERSAL,
                    mpi.MPI_COMM_WORLD,
                )

        if myid == main_node:
            # replace the name of the stack with reference with the current one
            refim = os.path.join(outdir, "aqm%03d.hdf" % Iter)
            a1 = 0.0
            ave_fsc = []
            for j in range(numref):
                if refi[j][2] < 4:
                    # ERROR("One of the references vanished","mref_ali2d_MPI",1)
                    #  if vanished, put a random image (only from main node!) there
                    assign[j] = []
                    assign[j].append(
                        random.randint(image_start, image_end - 1) - image_start
                    )
                    refi[j][0] = data[assign[j][0]].copy()
                    # print 'ERROR', j
                else:
                    if CTF:
                        for i in range(lctf):
                            ctm[i] = old_div(1.0, (ctf2[j][0][i] + old_div(1.0, snr)))
                        av1 = sp_filter.filt_table(refi[j][0], ctm)
                        for i in range(lctf):
                            ctm[i] = old_div(1.0, (ctf2[j][1][i] + old_div(1.0, snr)))
                        av2 = sp_filter.filt_table(refi[j][1], ctm)
                        # frsc = fsc_mask(av1, av2, mask, 1.0, os.path.join(outdir,"drm%03d%04d"%(Iter, j)))
                        frsc = sp_statistics.fsc(
                            av1,
                            av2,
                            1.0,
                            os.path.join(outdir, "drm%03d%04d.txt" % (Iter, j)),
                        )
                        # Now the total average
                        for i in range(lctf):
                            ctm[i] = old_div(
                                1.0, (ctf2[j][0][i] + ctf2[j][1][i] + old_div(1.0, snr))
                            )
                        refi[j][0] = sp_filter.filt_table(
                            EMAN2_cppwrap.Util.addn_img(refi[j][0], refi[j][1]), ctm
                        )
                    else:
                        # frsc = fsc_mask(refi[j][0], refi[j][1], mask, 1.0, os.path.join(outdir,"drm%03d%04d"%(Iter, j)))
                        frsc = sp_statistics.fsc(
                            refi[j][0],
                            refi[j][1],
                            1.0,
                            os.path.join(outdir, "drm%03d%04d.txt" % (Iter, j)),
                        )
                        EMAN2_cppwrap.Util.add_img(refi[j][0], refi[j][1])
                        EMAN2_cppwrap.Util.mul_scalar(
                            refi[j][0], old_div(1.0, float(refi[j][2]))
                        )

                    if ave_fsc == []:
                        for i in range(len(frsc[1])):
                            ave_fsc.append(frsc[1][i])
                        c_fsc = 1
                    else:
                        for i in range(len(frsc[1])):
                            ave_fsc[i] += frsc[1][i]
                        c_fsc += 1
                    # print 'OK', j, len(frsc[1]), frsc[1][0:5], ave_fsc[0:5]

            # print 'sum', sum(ave_fsc)
            if sum(ave_fsc) != 0:
                for i in range(len(ave_fsc)):
                    ave_fsc[i] = old_div(ave_fsc[i], float(c_fsc))
                    frsc[1][i] = ave_fsc[i]

            for j in range(numref):
                ref_data[2] = refi[j][0]
                ref_data[3] = frsc
                refi[j][0], cs = user_func(ref_data)

                # write the current average
                TMP = []
                for i_tmp in range(len(assign[j])):
                    TMP.append(float(assign[j][i_tmp]))
                TMP.sort()
                refi[j][0].set_attr_dict({"ave_n": refi[j][2], "members": TMP})
                del TMP
                refi[j][0].process_inplace(
                    "normalize.mask", {"mask": mask, "no_sigma": 1}
                )
                refi[j][0].write_image(refim, j)
                a1 += refi[j][0].cmp("dot", refi[j][0], {"negative": 0, "mask": mask})
                # send refi[j][0]  back!

            Iter += 1
            msg = "ITERATION #%3d        criterion = %20.7e\n\n" % (Iter, a1)
            sp_utilities.print_msg(msg)
            for j in range(numref):
                msg = "   group #%3d   number of particles = %7d\n" % (j, refi[j][2])
                sp_utilities.print_msg(msg)

            if a1 < a0:
                if auto_stop == True:
                    again = False
            else:
                a0 = a1
        # again = mpi_bcast(again, 1, MPI_INT, main_node, MPI_COMM_WORLD)
        Iter = sp_utilities.bcast_number_to_all(Iter, main_node)
        if CTF:
            del ctf2
        if again:
            for j in range(numref):
                sp_utilities.bcast_EMData_to_all(refi[j][0], myid, main_node)

    #  clean up
    del assign
    # write out headers  and STOP, under MPI writing has to be done sequentially
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    if CTF and data_had_ctf == 0:
        for im in range(len(data)):
            data[im].set_attr("ctf_applied", 0)
    par_str = ["xform.align2d", "assign", "ID"]
    if myid == main_node:
        if sp_utilities.file_type(stack) == "bdb":
            sp_utilities.recv_attr_dict_bdb(
                main_node, stack, data, par_str, image_start, image_end, number_of_proc
            )
        else:
            sp_utilities.recv_attr_dict(
                main_node, stack, data, par_str, image_start, image_end, number_of_proc
            )
    else:
        sp_utilities.send_attr_dict(main_node, data, par_str, image_start, image_end)

    if myid == main_node:
        newrefim = os.path.join(outdir, "multi_ref.hdf")
        for j in range(numref):
            refi[j][0].write_image(newrefim, j)
        sp_utilities.print_end_msg("mref_ali2d_MPI")


"""Multiline Comment2"""


"""Multiline Comment3"""


"""Multiline Comment9"""


# from development import sali3d_base_h_01
# sali3d_base = sali3d_base_horatio_01


"""Multiline Comment14"""


# Auxiliary function to compute number of cones in ali3dlocal


"""Multiline Comment15"""


"""Multiline Comment32"""

def transform2d(stack_data, stack_data_ali, shift = False, ignore_mirror = False, method = "quadratic"):
# apply 2D alignment parameters stored in the header of the input stack file using gridding interpolation and create an output stack file
    from ..libpy.sp_fundamentals   import rot_shift2D
    from ..libpy.sp_utilities 	    import set_params2D, get_params2D, get_im
    import os
    if  shift:
        from ..libpy.sp_utilities     import compose_transform2m
        from ..libpy.sp_fundamentals  import fshift, mirror

    t = EMAN2_cppwrap.Transform({"type":"2D"})# Transform({"type":"2D"})
    nima = EMAN2_cppwrap.EMUtil.get_image_count(stack_data)
    for im in range(nima):
        data = get_im(stack_data, im)
        al2d = get_params2D(data)
        if(shift):
            angb, sxb, syb, nm, ct = compose_transform2m(0.0, al2d[1], al2d[2], 0, 1.0, -al2d[0], 0.0, 0.0, al2d[3], 1.0)
            data = fshift(data, sxb, syb)
            if ignore_mirror: nm = 0
            if(nm == 1):  data = mirror(data)
        else:
            if ignore_mirror: al2d[3] = 0
            data = rot_shift2D(data, al2d[0], al2d[1], al2d[2], al2d[3], al2d[4], interpolation_method = method)
        data.set_attr("xform.align2d", t)
        data.write_image(stack_data_ali, im)


def cpy(ins_list, ous):
    # reworked to include lists, since we want to be able to copy lists of images
    #    into one file, concatenating.
    if isinstance(ins_list, list):
        # got a list of input files
        image_list = ins_list
    else:
        # got a single input file
        image_list = [ins_list]

    gl_index = 0

    oextension = sp_utilities.file_type(ous)

    if oextension == "bdb":
        DB = EMAN2db.db_open_dict(ous)

    # iterate over all images in the list, even if it's only one...
    for ins in image_list:

        # print ins
        nima = EMAN2_cppwrap.EMUtil.get_image_count(ins)
        iextension = sp_utilities.file_type(ins)

        if nima == 1 and oextension == "spi":
            sp_utilities.get_im(ins).write_image(
                ous, 0, EMAN2_cppwrap.EMUtil.ImageType.IMAGE_SINGLE_SPIDER
            )

        elif iextension == "bdb" and oextension == "bdb":

            OB = EMAN2db.db_open_dict(ins)
            for i in range(nima):
                DB[gl_index] = OB[i]
                gl_index += 1
            OB.close()

        elif iextension == "bdb":

            DB = EMAN2db.db_open_dict(ins)
            for i in range(nima):
                a = DB[i]
                a.write_image(ous, gl_index)
                gl_index += 1
            DB.close()

        elif oextension == "bdb":

            for i in range(nima):
                DB[gl_index] = sp_utilities.get_im(ins, i)
                gl_index += 1

        else:
            for im in range(nima):
                sp_utilities.get_im(ins, im).write_image(ous, gl_index)
                gl_index += 1

    if oextension == "bdb":
        DB.close()


def project3d(
    volume,
    stack=None,
    mask=None,
    delta=5,
    method="S",
    phiEqpsi="Minus",
    symmetry="c1",
    listagls=None,
    listctfs=None,
    noise=None,
    realsp=False,
    trillinear=False,
):
    if trillinear and realsp:
        sp_global_def.ERROR(
            "Both trilinear mode and realsp mode are specified", "project3d", 1
        )

    if listagls is None:
        angles = sp_utilities.even_angles(
            delta, symmetry=symmetry, method=method, phiEqpsi=phiEqpsi
        )
    elif isinstance(listagls, (bytes, str)):
        angles = sp_utilities.read_text_row(listagls, "", "")
    else:
        angles = listagls

    # try to parse the CTFs list. this is either not set (None), a filename or a list of values
    if listctfs is None:
        # not set, so simply ignore it
        ctfs = None
    elif isinstance(listctfs, (bytes, str)):
        # a string, so assume this is a filename and try to open the file
        try:
            ctfs = sp_utilities.read_text_row(listctfs, "", "")
        except:
            ctfs = [None for ii in range(len(angles))]
    else:
        # assume this a list of len(angles)
        ctfs = listctfs

    if not noise is None:
        # try to convert noise string to float. ignore noise if this fails
        try:
            noise_level = float(noise)
        except:
            noise_level = None
    # ignore noise, since it was not requested
    else:
        noise_level = None

    if isinstance(volume, (bytes, str)):
        vol = EMAN2_cppwrap.EMData()
        vol.read_image(volume)
        if mask:
            if isinstance(mask, (bytes, str)):
                maski = EMAN2_cppwrap.EMData()
                maski.read_image(volume)
                EMAN2_cppwrap.Util.mul_img(vol, maski)
                del maski
            else:
                EMAN2_cppwrap.Util.mul_img(vol, mask)
        nx = vol.get_xsize()
        ny = vol.get_ysize()
        nz = vol.get_zsize()

        if realsp:
            volft = vol
        elif trillinear:
            volft = sp_projection.prep_vol(vol, npad=2, interpolation_method=1)
        else:
            if nx == nz & ny == nz:
                volft, kb = sp_projection.prep_vol(vol)
            else:
                volft, kbx, kby, kbz = sp_projection.prep_vol(vol)
    else:
        vol = volume
        if mask:
            vol = vol.copy()
            if isinstance(mask, (bytes, str)):
                maski = EMAN2_cppwrap.EMData()
                maski.read_image(volume)
                EMAN2_cppwrap.Util.mul_img(vol, maski)
                del maski
            else:
                EMAN2_cppwrap.Util.mul_img(vol, mask)
        nx = vol.get_xsize()
        ny = vol.get_ysize()
        nz = vol.get_zsize()

        if realsp:
            volft = vol
        elif trillinear:
            volft = sp_projection.prep_vol(vol, npad=2, interpolation_method=1)
        else:
            if nx == nz & ny == nz:
                volft, kb = sp_projection.prep_vol(vol)
            else:
                volft, kbx, kby, kbz = sp_projection.prep_vol(vol)

    if isinstance(stack, (bytes, str)):
        Disk = True
        os.system("rm -f  " + stack)
    else:
        out = []
        Disk = False

    s2x = 0
    s2y = 0

    for i in range(len(angles)):
        if len(angles[i]) == 3:
            if realsp:
                proj = sp_projection.project(
                    volft, [angles[i][0], angles[i][1], angles[i][2], 0.0, 0.0], 10 * nx
                )
            elif trillinear:
                if ctfs is not None:
                    proj = sp_projection.prgl(
                        volft,
                        [angles[i][0], angles[i][1], angles[i][2], 0.0, 0.0],
                        1,
                        False,
                    )
                else:
                    proj = sp_projection.prgl(
                        volft,
                        [angles[i][0], angles[i][1], angles[i][2], 0.0, 0.0],
                        1,
                        True,
                    )
            else:
                if nx == nz & ny == nz:
                    proj = sp_projection.prgs(
                        volft, kb, [angles[i][0], angles[i][1], angles[i][2], 0.0, 0.0]
                    )
                else:
                    proj = sp_projection.prgs(
                        volft,
                        kbz,
                        [angles[i][0], angles[i][1], angles[i][2], 0.0, 0.0],
                        kbx,
                        kby,
                    )
            sp_utilities.set_params_proj(
                proj, [angles[i][0], angles[i][1], angles[i][2], 0.0, 0.0]
            )
        else:
            if realsp:
                proj = sp_projection.project(
                    volft,
                    [
                        angles[i][0],
                        angles[i][1],
                        angles[i][2],
                        -angles[i][3],
                        -angles[i][4],
                    ],
                    10 * nx,
                )
            elif trillinear:
                if ctfs is not None:
                    proj = sp_projection.prgl(
                        volft,
                        [
                            angles[i][0],
                            angles[i][1],
                            angles[i][2],
                            -angles[i][3],
                            -angles[i][4],
                        ],
                        1,
                        False,
                    )
                else:
                    proj = sp_projection.prgl(
                        volft,
                        [
                            angles[i][0],
                            angles[i][1],
                            angles[i][2],
                            -angles[i][3],
                            -angles[i][4],
                        ],
                        1,
                        True,
                    )
            else:
                if nx == nz & ny == nz:
                    proj = sp_projection.prgs(
                        volft,
                        kb,
                        [
                            angles[i][0],
                            angles[i][1],
                            angles[i][2],
                            -angles[i][3],
                            -angles[i][4],
                        ],
                    )
                else:
                    proj = sp_projection.prgs(
                        volft,
                        kbz,
                        [
                            angles[i][0],
                            angles[i][1],
                            angles[i][2],
                            -angles[i][3],
                            -angles[i][4],
                        ],
                        kbx,
                        kby,
                    )
            if trillinear:
                sp_utilities.set_params_proj(proj, angles[i][0:5])
            else:
                sp_utilities.set_params_proj(proj, angles[i])
        # horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
        # proj.set_attr_dict({'active':1})

        # add noise, if noise is set. this is two-fold: application of noise before
        #    ctf filtering and after it.
        if noise is not None:
            try:
                # no mask, so call w/ false
                noise_ima = sp_utilities.model_gauss_noise(
                    noise_level, proj.get_xsize(), proj.get_ysize()
                )
            except:
                pass
            else:
                proj += noise_ima

        # apply ctf, if ctf option is set and if we can create a valid CTF object
        if ctfs is not None:
            try:
                if len(ctfs[i]) == 6:
                    ctf = sp_utilities.generate_ctf(
                        [
                            ctfs[i][0],
                            ctfs[i][1],
                            ctfs[i][2],
                            ctfs[i][3],
                            ctfs[i][4],
                            ctfs[i][5],
                        ]
                    )
                elif len(ctfs[i]) == 8:
                    ctf = sp_utilities.generate_ctf(
                        [
                            ctfs[i][0],
                            ctfs[i][1],
                            ctfs[i][2],
                            ctfs[i][3],
                            ctfs[i][4],
                            ctfs[i][5],
                            ctfs[i][6],
                            ctfs[i][7],
                        ]
                    )
                else:
                    old_div(1.0, 0.0)
            except:
                # there are no ctf values, so ignore this and set no values
                sp_global_def.ERROR("Incorrect ctf values", "project3d", 1)
            # setting of values worked, so apply ctf and set the header info correctly
            if trillinear:
                if ctfs is not None:
                    proj.set_attr_dict({"is_complex": 0})
                    EMAN2_cppwrap.Util.mulclreal(
                        proj, sp_morphology.ctf_img_real(proj.get_ysize(), ctf)
                    )
                proj.set_attr_dict({"padffted": 1, "is_complex": 1})
                proj = sp_fundamentals.fft(proj)
            else:
                proj = sp_filter.filt_ctf(proj, ctf)
            proj.set_attr("ctf", ctf)
            proj.set_attr("ctf_applied", 0)

        # add second noise level that is not affected by CTF
        if noise is not None:
            try:
                noise_ima = sp_utilities.model_gauss_noise(
                    noise_level, proj.get_xsize(), proj.get_ysize()
                )
            except:
                pass
            else:
                proj += noise_ima

        if Disk:
            proj.write_image(stack, i)
        else:
            out.append(proj)
    if not Disk:
        return out


def ali_vol(vol, refv, ang_scale, shift_scale, radius=None, discrepancy="ccc"):
    """
		Name
			sxali_vol - align a 3D structure with respect of a 3D reference structure
		Input
			aligned_volume.hdf: 3D structure to be aligned.
		reference_volume.hdf
			3D reference structure.
		ang_scale
			correct angles are expected to be within +/-ang_scale of the values preset in the header of the structure to be aligned
		shift_scale
			correct shifts are expected to be within +/-shift_scale of the values preset in the header of the structure to be aligned
		mag_scale
			correct magnification is expected to be within +/-mag_scale of the value preset in the header of the structure to be aligned
		r
			radius of a spherical mask centered at nx/2, ny/2, nz/2
		Note - there are no defaults for three scale parameters. At least one has to appear.
	"""

    # rotation and shift

    ref = sp_utilities.get_image(refv)
    nx = ref.get_xsize()
    ny = ref.get_ysize()
    nz = ref.get_zsize()
    if radius != None:
        mask = sp_utilities.model_circle(radius, nx, ny, nz)
    else:
        mask = sp_utilities.model_circle(
            float(old_div(min(nx, ny, nz), 2) - 2), nx, ny, nz
        )

    # names_params = ["phi", "theta", "psi", "s3x", "s3y", "s3z", "scale"]
    phi, theta, psi, s3x, s3y, s3z, mirror, scale = sp_utilities.get_params3D(ref)
    paramsr = [phi, theta, psi, s3x, s3y, s3z, mirror, scale]
    # print  " params of the reference volume", paramsr
    ref = sp_fundamentals.rot_shift3D(
        ref,
        paramsr[0],
        paramsr[1],
        paramsr[2],
        paramsr[3],
        paramsr[4],
        paramsr[5],
        paramsr[7],
    )

    e = sp_utilities.get_image(vol)
    phi, theta, psi, s3x, s3y, s3z, mirror, scale = sp_utilities.get_params3D(e)
    paramsv = [phi, theta, psi, s3x, s3y, s3z, mirror, scale]
    e = sp_fundamentals.rot_shift3D(e, phi, theta, psi, s3x, s3y, s3z, scale)
    # print  " input params ", paramsv
    params = [phi, theta, psi, s3x, s3y, s3z]
    data = [e, ref, mask, params, discrepancy]
    new_params = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    new_params = sp_utilities.amoeba(
        new_params,
        [ang_scale, ang_scale, ang_scale, shift_scale, shift_scale, shift_scale],
        sp_alignment.ali_vol_func,
        1.0e-4,
        1.0e-4,
        500,
        data,
    )
    sp_global_def.sxprint(
        "amoeba: func_value =", new_params[1], "iter =", new_params[2]
    )

    cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale = sp_utilities.compose_transform3(
        paramsv[0],
        paramsv[1],
        paramsv[2],
        paramsv[3],
        paramsv[4],
        paramsv[5],
        paramsv[7],
        new_params[0][0],
        new_params[0][1],
        new_params[0][2],
        new_params[0][3],
        new_params[0][4],
        new_params[0][5],
        1.0,
    )
    # print  " new params ", cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale, new_params[1]
    sp_utilities.set_params3D(e, [cphi, ctheta, cpsi, cs2x, cs2y, cs2z, 0, cscale])
    if type(vol) == type(""):
        sp_utilities.write_headers(vol, [e], [0])
    else:
        return e


def recons3d_n(
    prj_stack,
    pid_list,
    vol_stack,
    CTF=False,
    snr=1.0,
    sign=1,
    npad=4,
    sym="c1",
    listfile="",
    group=-1,
    verbose=0,
    MPI=False,
    xysize=-1,
    zsize=-1,
    smearstep=0.0,
    upweighted=False,
    compensate=False,
    chunk_id=-1,
):
    if MPI:
        recons3d_n_MPI(
            prj_stack,
            pid_list,
            vol_stack,
            CTF,
            snr,
            1,
            npad,
            sym,
            listfile,
            group,
            verbose,
            xysize,
            zsize,
            smearstep,
        )
        ##newrecons3d_n_MPI(prj_stack, pid_list, vol_stack, CTF, snr, 1, npad, sym, listfile, group, verbose,xysize, zsize)
        return

    if listfile:
        pid_list = sp_utilities.read_text_file(listfile, 0)
        pid_list = list(map(int, pid_list))
    elif group > -1:
        tmp_list = EMAN2_cppwrap.EMUtil.get_all_attributes(prj_stack, "group")
        pid_list = []
        for i in range(len(tmp_list)):
            if tmp_list[i] == group:
                pid_list.append(i)
        del tmp_list

    if CTF:
        vol = sp_reconstruction.recons3d_4nn_ctf(
            prj_stack, pid_list, snr, 1, sym, verbose, npad, xysize=xysize, zsize=zsize
        )
    else:
        vol = sp_reconstruction.recons3d_4nn(
            prj_stack, pid_list, sym, npad, snr=snr, xysize=xysize, zsize=zsize
        )
    if vol_stack[-3:] == "spi":
        sp_utilities.drop_image(vol, vol_stack, "s")
    else:
        sp_utilities.drop_image(vol, vol_stack)


def recons3d_n_MPI(
    prj_stack,
    pid_list,
    vol_stack,
    CTF=False,
    snr=1.0,
    sign=1,
    npad=2,
    sym="c1",
    listfile="",
    group=-1,
    verbose=0,
    xysize=-1,
    zsize=-1,
    smearstep=0.0,
):

    myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
    nproc = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
    time_start = time.time()

    if myid == 0:
        if listfile:
            pid_list = sp_utilities.read_text_file(listfile, 0)
        elif group > -1:
            tmp_list = EMAN2_cppwrap.EMUtil.get_all_attributes(prj_stack, "group")
            pid_list = []
            for i in range(len(tmp_list)):
                if tmp_list[i] == group:
                    pid_list.append(i)
            del tmp_list
        nima = len(pid_list)
    else:
        nima = 0

    nima = sp_utilities.bcast_number_to_all(nima, source_node=0)

    if listfile or group > -1:
        if myid != 0:
            pid_list = [-1] * nima
        pid_list = mpi.mpi_bcast(pid_list, nima, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD)
        pid_list = list(map(int, pid_list))
    else:
        if not pid_list:
            pid_list = list(range(nima))

    if verbose == 0:
        finfo = None
    else:
        infofile = "progress%04d.txt" % (myid + 1)
        finfo = open(infofile, "w")

    image_start, image_end = MPI_start_end(nima, nproc, myid)

    prjlist = sp_utilities.iterImagesStack(prj_stack, pid_list[image_start:image_end])
    del pid_list

    if CTF:
        vol = sp_reconstruction.recons3d_4nn_ctf_MPI(
            myid,
            prjlist,
            snr,
            sign,
            sym,
            finfo,
            npad,
            xysize,
            zsize,
            smearstep=smearstep,
        )
    else:
        vol = sp_reconstruction.recons3d_4nn_MPI(
            myid, prjlist, sym, finfo, snr, npad, xysize, zsize
        )
    if myid == 0:
        if vol_stack[-3:] == "spi":
            sp_utilities.drop_image(vol, vol_stack, "s")
        else:
            sp_utilities.drop_image(vol, vol_stack)
        if not (finfo is None):
            finfo.write("result written to " + vol_stack + "\n")
            finfo.write("Total time: %10.3f\n" % (time.time() - time_start))
            finfo.flush()


def recons3d_trl_MPI(
    prj_stack,
    pid_list,
    vol_stack,
    CTF,
    snr,
    sign,
    npad,
    sym,
    verbose=None,
    niter=10,
    compensate=False,
    target_window_size=-1,
):
    # unregularized reconstruction  flags reconstruct(Iunreg(), gridding_nr_iter,false, 1., dummy, dummy, dummy, dummy, 1., false, true, nr_threads, -1
    myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
    nproc = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
    mpi_comm = mpi.MPI_COMM_WORLD
    time_start = time.time()
    nnxo = 0
    if myid == 0:
        nnxo = sp_utilities.get_im(prj_stack).get_xsize()
        sp_global_def.sxprint("  trilinear interpolation used in reconstruction")
        nima = len(pid_list)
    else:
        nima = 0
    nima = sp_utilities.bcast_number_to_all(nima, source_node=0)
    nnxo = sp_utilities.bcast_number_to_all(nnxo, source_node=0)
    if target_window_size == -1:
        target_size = 2 * nnxo + 3
    else:
        target_size = target_window_size * 2 + 3

    if verbose == 0:
        finfo = None
    else:
        infofile = "progress%04d.txt" % (myid + 1)
        finfo = open(infofile, "w")

    image_start, image_end = MPI_start_end(nima, nproc, myid)
    prjlist = EMAN2_cppwrap.EMData.read_images(
        prj_stack, pid_list[image_start:image_end]
    )

    # nnnx = ((prjlist[0].get_ysize())*2+3)

    # reconstruction step
    refvol = sp_utilities.model_blank(target_size)
    refvol.set_attr("fudge", 1.0)
    if CTF:
        do_ctf = 1
    else:
        do_ctf = 0
    if not (finfo is None):
        nimg = 0

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
    m = [1.0] * target_size
    is_complex = prjlist[0].get_attr("is_complex")
    for image in prjlist:
        if not is_complex:
            image = sp_fundamentals.fft(image)
        if CTF:
            image = sp_filter.filt_ctf(image, image.get_attr("ctf"))
        image.set_attr("padffted", 1)
        image.set_attr("npad", 1)
        image.set_attr("bckgnoise", m)
        r.insert_slice(image, image.get_attr_default("xform.projection", None), 1.0)
    if not (finfo is None):
        finfo.write("begin reduce\n")
        finfo.flush()

    sp_utilities.reduce_EMData_to_root(fftvol, myid, 0, comm=mpi_comm)
    sp_utilities.reduce_EMData_to_root(weight, myid, 0, comm=mpi_comm)

    if not (finfo is None):
        finfo.write("after reduce\n")
        finfo.flush()

    if myid == 0:
        dummy = r.finish(compensate)
    mpi.mpi_barrier(mpi_comm)

    if myid == 0:  # post-insertion operations, done only in main_node
        fftvol.set_attr("is_complex", 1)
        # no regularization here.
        if sym != "c1":
            fftvol = fftvol.symfvol(sym, -1)
            weight = weight.symfvol(sym, -1)  # symmetrize if not asymmetric
        maxr2 = ((old_div(target_window_size, -1)) * 2) ** 2
        EMAN2_cppwrap.Util.iterefa(fftvol, weight, maxr2, target_window_size)
        fftvol = sp_fundamentals.fft(
            sp_fundamentals.fshift(
                fftvol, target_window_size, target_window_size, target_window_size
            )
        )
        fftvol = EMAN2_cppwrap.Util.window(
            fftvol, target_window_size, target_window_size, target_window_size
        )
        fftvol = sp_morphology.cosinemask(
            fftvol, old_div(target_window_size, 2) - 1, 5, None
        )
        fftvol.div_sinc(1)
        fftvol.write_image(vol_stack)


def extract_value(s):

    try:
        i = int(s)
        return i
    except:
        pass

    try:
        f = float(s)
        return f
    except:
        pass

    return s


def header(
    stack,
    params,
    zero=False,
    one=False,
    set=0.0,
    randomize=False,
    rand_alpha=False,
    fimport=None,
    fexport=None,
    fprint=False,
    backup=False,
    suffix="_backup",
    restore=False,
    delete=False,
    consecutive=False,
    idx_list=None,
):

    if set == 0.0:
        doset = False
    else:
        doset = True

    op = (
        zero
        + one
        + +consecutive
        + randomize
        + rand_alpha
        + (fimport != None)
        + (fexport != None)
        + fprint
        + backup
        + restore
        + delete
        + doset
    )
    if op == 0:
        sp_global_def.sxprint("Error: no operation selected!")
        return
    elif op > 1:
        sp_global_def.sxprint("Error: more than one operation at the same time!")
        return

    params = params.split()

    if idx_list != None:
        fidx = sp_utilities.read_text_file(idx_list)
    else:
        fidx = None
    if fimport != None:
        fimp = open(fimport, "r")
    if fexport != None:
        fexp = open(fexport, "w")
    if idx_list and not fimport and not fexport:
        sp_global_def.sxprint(
            "Error: list option can only be used in combination with an import/export file"
        )
        return

    nimage = EMAN2_cppwrap.EMUtil.get_image_count(stack)
    ext = sp_utilities.file_type(stack)
    if ext == "bdb":
        DB = EMAN2db.db_open_dict(stack)
    for i in range(nimage):
        if fimport != None:
            if fidx:
                try:
                    i = fidx[i]
                except IndexError:
                    continue
            line = fimp.readline()
            if len(line) == 0:
                sp_global_def.sxprint(
                    "Error: file "
                    + fimport
                    + " has only "
                    + str(i)
                    + " lines, while there are "
                    + str(nimage)
                    + " images in the file."
                )
                return
            parmvalues = line.split()
            il = 0
            for p in params:
                if p[:13] == "xform.align2d":
                    if len(parmvalues) < il + 3:
                        sp_global_def.sxprint("Not enough parameters!")
                        return
                    alpha = extract_value(parmvalues[il])
                    sx = extract_value(parmvalues[il + 1])
                    sy = extract_value(parmvalues[il + 2])
                    if len(parmvalues) > (il + 3):
                        mirror = int(extract_value(parmvalues[il + 3]))
                    else:
                        mirror = 0
                    if len(parmvalues) > (il + 4):
                        scale = extract_value(parmvalues[il + 4])
                    else:
                        scale = 1.0
                    # set_params2D(img, [alpha, sx, sy, mirror, scale], params[0])
                    t = EMAN2_cppwrap.Transform(
                        {
                            "type": "2D",
                            "alpha": alpha,
                            "tx": sx,
                            "ty": sy,
                            "mirror": mirror,
                            "scale": scale,
                        }
                    )
                    if ext == "bdb":
                        DB.set_attr(i, "xform.align2d", t)
                    elif ext == "hdf":
                        EMAN2_cppwrap.EMUtil.write_hdf_attribute(
                            stack, "xform.align2d", t, i
                        )
                    il += 5

                elif p[:16] == "xform.projection":
                    if len(parmvalues) < il + 3:
                        sp_global_def.sxprint("Not enough parameters!")
                        return
                    phi = extract_value(parmvalues[il])
                    theta = extract_value(parmvalues[il + 1])
                    psi = extract_value(parmvalues[il + 2])
                    if len(parmvalues) > il + 3:
                        s2x = extract_value(parmvalues[il + 3])
                    else:
                        s2x = 0.0
                    if len(parmvalues) > il + 4:
                        s2y = extract_value(parmvalues[il + 4])
                    else:
                        s2y = 0.0
                    # set_params_proj(img, [phi, theta, psi, s2x, s2y], params[0])
                    t = EMAN2_cppwrap.Transform(
                        {"type": "spider", "phi": phi, "theta": theta, "psi": psi}
                    )
                    t.set_trans(EMAN2_cppwrap.Vec2f(-s2x, -s2y))
                    if ext == "bdb":
                        DB.set_attr(i, "xform.projection", t)
                    elif ext == "hdf":
                        EMAN2_cppwrap.EMUtil.write_hdf_attribute(
                            stack, "xform.projection", t, i
                        )
                    il = min(il + 5, len(parmvalues))
                elif p[:13] == "xform.align3d":
                    if len(parmvalues) < il + 8:
                        sp_global_def.sxprint("Not enough parameters!")
                        return
                    phi = extract_value(parmvalues[il])
                    theta = extract_value(parmvalues[il + 1])
                    psi = extract_value(parmvalues[il + 2])
                    s3x = extract_value(parmvalues[il + 3])
                    s3y = extract_value(parmvalues[il + 4])
                    s3z = extract_value(parmvalues[il + 5])
                    mirror = int(extract_value(parmvalues[il + 6]))
                    scale = extract_value(parmvalues[il + 7])
                    # set_params3D(img, [phi, theta, psi, s3x, s3y, s3z, mirror, scale], params[0])
                    t = EMAN2_cppwrap.Transform(
                        {
                            "type": "spider",
                            "phi": phi,
                            "theta": theta,
                            "psi": psi,
                            "tx": s3x,
                            "ty": s3y,
                            "tz": s3z,
                            "mirror": mirror,
                            "scale": scale,
                        }
                    )
                    if ext == "bdb":
                        DB.set_attr(i, "xform.align3d", t)
                    elif ext == "hdf":
                        EMAN2_cppwrap.EMUtil.write_hdf_attribute(
                            stack, "xform.align3d", t, i
                        )
                    il += 8
                elif p[: len("members")] == "members":
                    members = [int(entry) for entry in line.strip("[]\n \t").split(",")]
                    if ext == "bdb":
                        DB.set_attr(i, "members", members)
                    elif ext == "hdf":
                        EMAN2_cppwrap.EMUtil.write_hdf_attribute(
                            stack, "members", members, i
                        )
                elif p.startswith("ISAC_SPLIT_"):
                    if ext == "bdb":
                        DB.set_attr(i, p.strip(), line.strip())
                    elif ext == "hdf":
                        EMAN2_cppwrap.EMUtil.write_hdf_attribute(
                            stack, p.strip(), line.strip(), i
                        )
                elif p == "ctf":
                    if len(parmvalues) < il + 6:
                        sp_global_def.sxprint("Not enough parameters!")
                        return
                    defocus = extract_value(parmvalues[il])
                    cs = extract_value(parmvalues[il + 1])
                    voltage = extract_value(parmvalues[il + 2])
                    apix = extract_value(parmvalues[il + 3])
                    bfactor = extract_value(parmvalues[il + 4])
                    ampcont = extract_value(parmvalues[il + 5])
                    dfdiff = extract_value(parmvalues[il + 6])
                    dfang = extract_value(parmvalues[il + 7])
                    # set_ctf(img, [defocus, cs, voltage, apix, bfactor, ampcont])
                    ctf = sp_utilities.generate_ctf(
                        [defocus, cs, voltage, apix, bfactor, ampcont, dfdiff, dfang]
                    )
                    if ext == "bdb":
                        DB.set_attr(i, "ctf", ctf)
                    elif ext == "hdf":
                        EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, "ctf", ctf, i)
                    il += 6
                else:
                    # if len(params)!=len(parmvalues):
                    # print "Error: %d params need to be set, while %d values are provided in line %d of file." % ( len(params), len(parmvalues), i )
                    # return
                    if ext == "bdb":
                        DB.set_attr(i, p, extract_value(parmvalues[il]))
                    elif ext == "hdf":
                        EMAN2_cppwrap.EMUtil.write_hdf_attribute(
                            stack, p, extract_value(parmvalues[il]), i
                        )
                    il += 1

        else:
            if fidx:
                try:
                    i = fidx[i]
                except IndexError:
                    continue
            for p in params:

                if zero:
                    if p[:13] == "xform.align2d":
                        # set_params2D(img, [0.0, 0.0, 0.0, 0, 1.0], p)
                        t = EMAN2_cppwrap.Transform(
                            {
                                "type": "2D",
                                "alpha": 0.0,
                                "tx": 0.0,
                                "ty": 0.0,
                                "mirror": 0,
                                "scale": 1.0,
                            }
                        )
                        if ext == "bdb":
                            DB.set_attr(i, "xform.align2d", t)
                        elif ext == "hdf":
                            EMAN2_cppwrap.EMUtil.write_hdf_attribute(
                                stack, "xform.align2d", t, i
                            )
                    elif p[:16] == "xform.projection":
                        # set_params_proj(img, [0.0, 0.0, 0.0, 0.0, 0.0], p)
                        t = EMAN2_cppwrap.Transform({"type": "spider"})
                        if ext == "bdb":
                            DB.set_attr(i, "xform.projection", t)
                        elif ext == "hdf":
                            EMAN2_cppwrap.EMUtil.write_hdf_attribute(
                                stack, "xform.projection", t, i
                            )
                    elif p[:13] == "xform.align3d":
                        # set_params3D(img, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 1.0], p)
                        t = EMAN2_cppwrap.Transform({"type": "spider"})
                        if ext == "bdb":
                            DB.set_attr(i, "xform.align3d", t)
                        elif ext == "hdf":
                            EMAN2_cppwrap.EMUtil.write_hdf_attribute(
                                stack, "xform.align3d", t, i
                            )
                    elif p == "ctf":
                        sp_global_def.sxprint("Invalid operation!")
                        return
                    else:
                        # img.set_attr(p, 0)
                        if ext == "bdb":
                            DB.set_attr(i, p, 0)
                        elif ext == "hdf":
                            EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, p, 0, i)
                elif one:
                    if p[:6] == "xform." or p == "ctf":
                        sp_global_def.sxprint("Invalid operation!")
                        return
                    else:
                        # img.set_attr(p, 1)
                        if ext == "bdb":
                            DB.set_attr(i, p, 1)
                        elif ext == "hdf":
                            EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, p, 1, i)
                elif doset:
                    if p[:6] == "xform." or p == "ctf":
                        sp_global_def.sxprint("Invalid operation!")
                        return
                    else:
                        # img.set_attr(p, 1)
                        if ext == "bdb":
                            DB.set_attr(i, p, set)
                        elif ext == "hdf":
                            EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, p, set, i)
                elif consecutive:
                    if ext == "bdb":
                        DB.set_attr(i, p, i)
                    elif ext == "hdf":
                        EMAN2_cppwrap.EMUtil.write_hdf_attribute(stack, p, i, i)
                elif randomize:
                    if p[:13] == "xform.align2d":
                        alpha = random.random() * 360.0
                        sx = random.random() * 2.0 - 1.0
                        sy = random.random() * 2.0 - 1.0
                        mirror = random.randint(0, 1)
                        scale = 1.0
                        # set_params2D(img, [alpha, sx, sy, mirror, scale], p)
                        t = EMAN2_cppwrap.Transform(
                            {
                                "type": "2D",
                                "alpha": alpha,
                                "tx": sx,
                                "ty": sy,
                                "mirror": mirror,
                                "scale": scale,
                            }
                        )
                        if ext == "bdb":
                            DB.set_attr(i, "xform.align2d", t)
                        elif ext == "hdf":
                            EMAN2_cppwrap.EMUtil.write_hdf_attribute(
                                stack, "xform.align2d", t, i
                            )
                    elif p[:16] == "xform.projection":
                        phi = random.random() * 360.0
                        theta = random.random() * 180.0
                        psi = random.random() * 360.0
                        s2x = random.random() * 4.0 - 2.0
                        s2y = random.random() * 4.0 - 2.0
                        # set_params_proj(img, [phi, theta, psi, s2x, s2y], p)
                        t = EMAN2_cppwrap.Transform(
                            {"type": "spider", "phi": phi, "theta": theta, "psi": psi}
                        )
                        t.set_trans(EMAN2_cppwrap.Vec2f(-s2x, -s2y))
                        if ext == "bdb":
                            DB.set_attr(i, "xform.projection", t)
                        elif ext == "hdf":
                            EMAN2_cppwrap.EMUtil.write_hdf_attribute(
                                stack, "xform.projection", t, i
                            )
                    elif p[:13] == "xform.align3d":
                        phi = random.random() * 360.0
                        theta = random.random() * 180.0
                        psi = random.random() * 360.0
                        s3x = random.random() * 4.0 - 2.0
                        s3y = random.random() * 4.0 - 2.0
                        s3z = random.random() * 4.0 - 2.0
                        mirror = random.randint(0, 1)
                        scale = 1.0
                        # set_params3D(img, [phi, theta, psi, s3x, s3y, s3z, mirror, scale], p)
                        t = EMAN2_cppwrap.Transform(
                            {
                                "type": "spider",
                                "phi": phi,
                                "theta": theta,
                                "psi": psi,
                                "tx": s3x,
                                "ty": s3y,
                                "tz": s3z,
                                "mirror": mirror,
                                "scale": scale,
                            }
                        )
                        if ext == "bdb":
                            DB.set_attr(i, "xform.align3d", t)
                        elif ext == "hdf":
                            EMAN2_cppwrap.EMUtil.write_hdf_attribute(
                                stack, "xform.align3d", t, i
                            )
                    else:
                        sp_global_def.sxprint("Invalid operation!")
                        return
                elif rand_alpha:
                    if p[:13] == "xform.align2d":
                        alpha = random.random() * 360.0
                        sx = 0.0
                        sy = 0.0
                        mirror = random.randint(0, 1)
                        scale = 1.0
                        # set_params2D(img, [alpha, sx, sy, mirror, scale], p)
                        t = EMAN2_cppwrap.Transform(
                            {
                                "type": "2D",
                                "alpha": alpha,
                                "tx": sx,
                                "ty": sy,
                                "mirror": mirror,
                                "scale": scale,
                            }
                        )
                        if ext == "bdb":
                            DB.set_attr(i, "xform.align2d", t)
                        elif ext == "hdf":
                            EMAN2_cppwrap.EMUtil.write_hdf_attribute(
                                stack, "xform.align2d", t, i
                            )
                    elif p[:16] == "xform.projection":
                        phi = random.random() * 360.0
                        theta = random.random() * 180.0
                        psi = random.random() * 360.0
                        s2x = 0.0
                        s2y = 0.0
                        # set_params_proj(img, [phi, theta, psi, s2x, s2y], p)

                        t = EMAN2_cppwrap.Transform(
                            {"type": "spider", "phi": phi, "theta": theta, "psi": psi}
                        )
                        t.set_trans(EMAN2_cppwrap.Vec2f(-s2x, -s2y))
                        if ext == "bdb":
                            DB.set_attr(i, "xform.projection", t)
                        elif ext == "hdf":
                            EMAN2_cppwrap.EMUtil.write_hdf_attribute(
                                stack, "xform.projection", t, i
                            )
                    elif p[:13] == "xform.align3d":
                        phi = random.random() * 360.0
                        theta = random.random() * 180.0
                        psi = random.random() * 360.0
                        s3x = 0.0
                        s3y = 0.0
                        s3z = 0.0
                        mirror = random.randint(0, 1)
                        scale = 1.0
                        # set_params3D(img, [phi, theta, psi, s3x, s3y, s3z, mirror, scale], p)
                        t = EMAN2_cppwrap.Transform(
                            {
                                "type": "spider",
                                "phi": phi,
                                "theta": theta,
                                "psi": psi,
                                "tx": s3x,
                                "ty": s3y,
                                "tz": s3z,
                                "mirror": mirror,
                                "scale": scale,
                            }
                        )
                        if ext == "bdb":
                            DB.set_attr(i, "xform.align3d", t)
                        elif ext == "hdf":
                            EMAN2_cppwrap.EMUtil.write_hdf_attribute(
                                stack, "xform.align3d", t, i
                            )
                    else:
                        sp_global_def.sxprint("Invalid operation!")
                        return

                elif fexport != None:
                    if p[:13] == "xform.align2d":
                        # alpha, sx, sy, mirror, scale = get_params2D(img, p)
                        if ext == "bdb":
                            t = DB.get_attr(i, "xform.align2d")
                            d = t.get_params("2D")
                            fexp.write(
                                "%16.6f %16.6f %16.6f %10d %10.3f "
                                % (
                                    d["alpha"],
                                    d["tx"],
                                    d["ty"],
                                    d["mirror"],
                                    d["scale"],
                                )
                            )

                        elif ext == "hdf":
                            t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(
                                stack, "xform.align2d", i
                            )
                            d = t.get_params("2D")
                            fexp.write(
                                "%16.6f %16.6f %16.6f %10d %10.3f "
                                % (
                                    d["alpha"],
                                    d["tx"],
                                    d["ty"],
                                    d["mirror"],
                                    d["scale"],
                                )
                            )

                    elif p[:16] == "xform.projection":
                        # phi, theta, psi, s2x, s2y = get_params_proj(img, p)
                        if ext == "bdb":
                            t = DB.get_attr(i, "xform.projection")
                            d = t.get_params("spider")
                            fexp.write(
                                "%16.6f %16.6f %16.6f %16.6f %16.6f "
                                % (d["phi"], d["theta"], d["psi"], -d["tx"], -d["ty"])
                            )

                        elif ext == "hdf":
                            t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(
                                stack, "xform.projection", i
                            )
                            d = t.get_params("spider")
                            fexp.write(
                                "%16.6f %16.6f %16.6f %16.6f %16.6f "
                                % (d["phi"], d["theta"], d["psi"], -d["tx"], -d["ty"])
                            )

                    elif p[:13] == "xform.align3d":
                        # phi, theta, psi, s3x, s3y, s3z, mirror, scale = get_params3D(img, p)
                        if ext == "bdb":
                            t = DB.get_attr(i, "xform.align3d")
                            d = t.get_params("spider")
                            fexp.write(
                                "%16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %10d %10.3f "
                                % (
                                    d["phi"],
                                    d["theta"],
                                    d["psi"],
                                    d["tx"],
                                    d["ty"],
                                    d["tz"],
                                    d["mirror"],
                                    d["scale"],
                                )
                            )

                        elif ext == "hdf":
                            t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(
                                stack, "xform.align3d", i
                            )
                            d = t.get_params("spider")
                            fexp.write(
                                "%16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %10d %10.3f "
                                % (
                                    d["phi"],
                                    d["theta"],
                                    d["psi"],
                                    d["tx"],
                                    d["ty"],
                                    d["tz"],
                                    d["mirror"],
                                    d["scale"],
                                )
                            )

                    elif p == "ctf":
                        # defocus, cs, voltage, apix, bfactor, ampcont = get_ctf(img)
                        if ext == "bdb":
                            t = DB.get_attr(i, "ctf")
                        elif ext == "hdf":
                            t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, "ctf", i)
                        fexp.write(
                            "%16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f"
                            % (
                                t.defocus,
                                t.cs,
                                t.voltage,
                                t.apix,
                                t.bfactor,
                                t.ampcont,
                                t.dfdiff,
                                t.dfang,
                            )
                        )

                    else:
                        if ext == "bdb":
                            fexp.write("%15s " % str(DB.get_attr(i, p)))

                        elif ext == "hdf":
                            fexp.write(
                                "%15s "
                                % str(
                                    EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, p, i)
                                )
                            )
                elif fprint:
                    if p[:13] == "xform.align2d":
                        # alpha, sx, sy, mirror, scale = get_params2D(img, p)
                        if ext == "bdb":
                            t = DB.get_attr(i, "xform.align2d")
                            d = t.get_params("2D")
                            print(
                                "%16.6f %16.6f %16.6f %10d %10.3f"
                                % (
                                    d["alpha"],
                                    d["tx"],
                                    d["ty"],
                                    d["mirror"],
                                    d["scale"],
                                ),
                                end=" ",
                            )

                        elif ext == "hdf":
                            t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(
                                stack, "xform.align2d", i
                            )
                            d = t.get_params("2D")
                            print(
                                "%16.6f %16.6f %16.6f %10d %10.3f"
                                % (
                                    d["alpha"],
                                    d["tx"],
                                    d["ty"],
                                    d["mirror"],
                                    d["scale"],
                                ),
                                end=" ",
                            )
                    elif p[:16] == "xform.projection":
                        # phi, theta, psi, s2x, s2y = get_params_proj(img, p)
                        if ext == "bdb":
                            t = DB.get_attr(i, "xform.projection")
                            d = t.get_params("spider")
                            print(
                                "%16.6f %16.6f %16.6f %16.6f %16.6f"
                                % (d["phi"], d["theta"], d["psi"], -d["tx"], -d["ty"]),
                                end=" ",
                            )
                        elif ext == "hdf":
                            t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(
                                stack, "xform.projection", i
                            )
                            d = t.get_params("spider")
                            print(
                                "%16.6f %16.6f %16.6f %16.6f %16.6f"
                                % (d["phi"], d["theta"], d["psi"], -d["tx"], -d["ty"]),
                                end=" ",
                            )
                    elif p[:13] == "xform.align3d":
                        # phi, theta, psi, s3x, s3y, s3z, mirror, scale = get_params3D(img, p)
                        if ext == "bdb":
                            t = DB.get_attr(i, "xform.align3d")
                            d = t.get_params("spider")
                            print(
                                "%16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %10d %10.3f"
                                % (
                                    d["phi"],
                                    d["theta"],
                                    d["psi"],
                                    d["tx"],
                                    d["ty"],
                                    d["tz"],
                                    d["mirror"],
                                    d["scale"],
                                ),
                                end=" ",
                            )
                        elif ext == "hdf":
                            t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(
                                stack, "xform.align3d", i
                            )
                            d = t.get_params("spider")
                            print(
                                "%16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %10d %10.3f"
                                % (
                                    d["phi"],
                                    d["theta"],
                                    d["psi"],
                                    d["tx"],
                                    d["ty"],
                                    d["tz"],
                                    d["mirror"],
                                    d["scale"],
                                ),
                                end=" ",
                            )
                    elif p == "ctf":
                        # defocus, cs, voltage, apix, bfactor, ampcont = get_ctf(img)
                        if ext == "bdb":
                            t = DB.get_attr(i, "ctf")
                        elif ext == "hdf":
                            t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, "ctf", i)
                        print(
                            "%16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f"
                            % (
                                t.defocus,
                                t.cs,
                                t.voltage,
                                t.apix,
                                t.bfactor,
                                t.ampcont,
                                t.dfdiff,
                                t.dfang,
                            ),
                            end=" ",
                        )

                    else:
                        if ext == "bdb":
                            print("%15s" % str(DB.get_attr(i, p)), end=" ")
                        elif ext == "hdf":
                            print(
                                "%15s"
                                % str(
                                    EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, p, i)
                                ),
                                end=" ",
                            )
                elif backup:
                    # t = img.get_attr(p)
                    # img.set_attr(p+suffix, t)
                    if ext == "bdb":
                        t = DB.get_attr(i, p)
                        DB.set_attr(i, p + suffix, t)
                    elif ext == "hdf":
                        t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, p, i)
                        EMAN2_cppwrap.EMUtil.write_hdf_attribute(
                            stack, p + suffix, t, i
                        )

                elif restore:
                    if (
                        p == "xform.align2d"
                        or p == "xform.align3d"
                        or p == "xform.projection"
                    ):
                        sp_global_def.sxprint("ERROR, no suffix in xform!")
                        return
                    # t = img.get_attr(p)
                    # if ext == "bdb":
                    # for i in xrange(nimage):
                    # t= DB.get_attr(i,p)
                    # elif ext == "hdf":
                    # for i in xrange(nimage):
                    # t= EMUtil.read_hdf_attribute(stack,p,i)
                    if p[:13] == "xform.align2d":
                        # img.set_attr(p[:13], t)
                        if ext == "bdb":
                            t = DB.get_attr(i, p)
                            DB.set_attr(i, "xform.align2d", t)
                        elif ext == "hdf":
                            t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, p, i)
                            EMAN2_cppwrap.EMUtil.write_hdf_attribute(
                                stack, "xform.align2d", t, i
                            )
                    elif p[:16] == "xform.projection":
                        # img.set_attr(p[:10], t)
                        if ext == "bdb":
                            t = DB.get_attr(i, p)
                            DB.set_attr(i, "xform.projection", t)
                        elif ext == "hdf":
                            t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, p, i)
                            EMAN2_cppwrap.EMUtil.write_hdf_attribute(
                                stack, "xform.projection", t, i
                            )
                    elif p[:13] == "xform.align3d":
                        # img.set_attr(p[:13], t)
                        if ext == "bdb":
                            t = DB.get_attr(i, p)
                            DB.set_attr(i, "xform.align3d", t)
                    elif ext == "hdf":
                        for i in range(nimage):
                            t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, p, i)
                            EMAN2_cppwrap.EMUtil.write_hdf_attribute(
                                stack, "xform.align3d", t, i
                            )
                    else:
                        # img.set_attr(p[:-len(suffix)], t)
                        if ext == "bdb":
                            t = DB.get_attr(i, p)
                            DB.set_attr(i, p[: -len(suffix)], t)
                        elif ext == "hdf":
                            t = EMAN2_cppwrap.EMUtil.read_hdf_attribute(stack, p, i)
                            EMAN2_cppwrap.EMUtil.write_hdf_attribute(
                                stack, p[: -len(suffix)], t, i
                            )
                elif delete:
                    img = EMAN2_cppwrap.EMData()
                    img.read_image(stack, i, True)
                    img.del_attr(p)
                    sp_utilities.write_header(stack, img, i)
            # if zero or one or randomize or rand_alpha or backup or restore or delete:
            # write_header(stack, img, i)
            if fexport != None:
                fexp.write("\n")
            if fprint:
                print(" ")
    if ext == "bdb":
        DB.close()


def MPI_start_end(nima, nproc, myid):
    image_start = int(round(old_div(float(nima), nproc) * myid))
    image_end = int(round(old_div(float(nima), nproc) * (myid + 1)))
    return image_start, image_end


"""Multiline Comment40"""


def refvol(vollist, fsclist, output, mask):

    nvol = len(vollist)
    assert len(fsclist) == nvol

    fscs = [None] * nvol
    vols = [None] * nvol
    for i in range(nvol):
        fscs[i] = sp_utilities.read_fsc(fsclist[i])
        vols[i] = sp_utilities.get_image(vollist[i])
        sp_global_def.sxprint("rawvol, resolution: ", vollist[i], fsclist[i])

    m = sp_utilities.get_image(mask)
    volfs = sp_filter.filt_vols(vols, fscs, m)

    for i in range(nvol):
        volfs[i].write_image(output, i)


# -- K-means main ---------------------------------------------------------------------------
# K-means main driver


"""Multiline Comment44"""


"""Multiline Comment45"""


#  Version before I started messing with centering of averages  PAP 07/10/2015


def within_group_refinement(
    data,
    maskfile,
    randomize,
    ir,
    ou,
    rs,
    xrng,
    yrng,
    step,
    dst,
    maxit,
    FH,
    FF,
    method="",
    CTF=False,
):

    # Comment by Zhengfan Yang 03/11/11
    # This is a simple version of ali2d_data (down to the bone), no output dir, no logfile, no MPI or CUDA, no Fourvar, no auto stop, no user function
    #  Added CTF fot simple version PAP 03/30/2017

    first_ring = int(ir)
    last_ring = int(ou)
    rstep = int(rs)
    max_iter = int(maxit)
    nima = len(data)
    nx = data[0].get_xsize()
    if last_ring == -1:
        last_ring = old_div(nx, 2) - 2
    if maskfile:
        mask = maskfile
    else:
        mask = sp_utilities.model_circle(last_ring, nx, nx)

    if randomize:
        for im in data:
            alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(im)
            alphai, sxi, syi, mirrori = sp_utilities.inverse_transform2(alpha, sx, sy)
            alphan, sxn, syn, mirrorn = sp_utilities.combine_params2(
                0.0,
                -sxi,
                -syi,
                0,
                random.random() * 360.0,
                0.0,
                0.0,
                random.randint(0, 1),
            )
            # alphan, sxn, syn, mirrorn    = combine_params2(0.0, -sxi+randint(-xrng[0],xrng[0]), -syi+randint(-xrng[0],xrng[0]), 0, random()*360.0, 0, 0, randint(0, 1))
            sp_utilities.set_params2D(im, [alphan, sxn, syn, mirrorn, 1.0])

    cnx = old_div(nx, 2) + 1
    cny = cnx
    mode = "F"
    numr = sp_alignment.Numrinit(first_ring, last_ring, rstep, mode)
    wr = sp_alignment.ringwe(numr, mode)

    sx_sum = 0.0
    sy_sum = 0.0
    cs = [0.0] * 2
    total_iter = 0
    if method == "SHC":
        #  This is my failed attempt to use SHC for 2D alignment.
        #    Inexplicably, it did not do all that well.  While initially it converges fast
        #     and generally yields a very good solution, it converges to cluster of images scattered
        #     around the 'best' solution, i.e., the shifts are within a fraction of a pixel of what they
        #     should be and, as a result, some are in wrong positions and overall pixel error is large.
        #     Overall, Yang's method works much better, so I am leaving it at that.  PAP 01/22/2015
        for im in data:
            im.set_attr("previousmax", -1.0e23)
        tavg = sp_statistics.ave_series(data)
        for N_step in range(len(xrng)):
            nope = 0
            Iter = 0
            while nope < old_div(len(data), 1) and Iter < max_iter:
                total_iter += 1
                Iter += 1
                if FH > 0.0:
                    tavg = sp_filter.filt_tanl(sp_fundamentals.fft(tavg), FH, FF)
                    if xrng[0] > 0.0:
                        cs[0] = old_div(sx_sum, float(nima))
                    if yrng[0] > 0.0:
                        cs[1] = old_div(sy_sum, float(nima))
                    tavg = sp_fundamentals.fft(
                        sp_fundamentals.fshift(tavg, -cs[0], -cs[1])
                    )
                else:
                    tavg = sp_filter.filt_tanl(tavg, FH, FF)
                sx_sum, sy_sum, nope = sp_alignment.ali2d_single_iter(
                    data,
                    numr,
                    wr,
                    cs,
                    tavg,
                    cnx,
                    cny,
                    xrng[N_step],
                    yrng[N_step],
                    step[N_step],
                    mode=mode,
                    CTF=False,
                    random_method=method,
                )
                # print  "  iteration  shc   %03d   %03d   %7.2f    %7.2f  "%(total_iter,nope,cs[0],cs[1])
                # print total_iter,nope
                # for i in data:  print "  ",i.get_attr('previousmax'),
                # print "  "
                # tavg.write_image('tata.hdf',total_iter-1)
                tavg = sp_statistics.ave_series(data)
        """Multiline Comment46"""

    elif method == "PCP":
        stp = step[-1]
        rings = sp_alignment.prepref(
            data,
            sp_utilities.model_circle(old_div(nx, 2) - 1, nx, nx),
            cnx,
            cnx,
            numr,
            mode,
            xrng[0],
            xrng[0],
            stp,
        )
        sp_global_def.sxprint(" rings  ", len(rings))
        for im in range(len(data)):
            rings[im][0][0].set_attr("inx", nx)
        tavg = sp_statistics.ave_series(data)
        for N_step in range(len(xrng)):
            sp_global_def.sxprint(" xrng ", xrng[N_step])
            for Iter in range(max_iter):
                total_iter += 1
                if FH > 0.0:
                    fl = 0.1 + old_div((FH - 0.1) * Iter, float(max_iter - 1))
                    tavg = sp_filter.filt_tanl(tavg, fl, FF)
                    """Multiline Comment47"""
                else:
                    """Multiline Comment48"""
                cs = [0, 0]
                # print  "  iteration  std   %03d   %7.2f    %7.2f  "%(total_iter,cs[0],cs[1])
                if Iter % 4 != 0 or total_iter > max_iter * len(xrng) - 10:
                    delta = 0.0
                else:
                    delta = dst
                sx_sum, sy_sum, nope = sp_alignment.ali2d_single_iter(
                    rings,
                    numr,
                    wr,
                    cs,
                    tavg,
                    cnx,
                    cny,
                    xrng[N_step],
                    yrng[N_step],
                    step[N_step],
                    mode=mode,
                    CTF=False,
                    delta=delta,
                    random_method=method,
                )
                for im in range(len(data)):
                    alpha, tx, ty, mirror, scale = sp_utilities.get_params2D(
                        rings[im][0][0]
                    )
                    sp_utilities.set_params2D(data[im], [alpha, tx, ty, mirror, scale])
                tavg = sp_statistics.ave_series(data)
                # print  "tata ",total_iter-1,Util.infomask(tavg,None,True)
                # tavg.write_image('tata.hdf',total_iter-1)
    else:
        if CTF:
            ctf2 = EMAN2_cppwrap.EMData(nx, nx, 1, False)
            cdata = []
            for im in data:
                ctt = sp_morphology.ctf_img(nx, im.get_attr("ctf"))
                EMAN2_cppwrap.Util.add_img2(ctf2, ctt)
                cdata.append(
                    sp_fundamentals.fft(
                        EMAN2_cppwrap.Util.muln_img(sp_fundamentals.fft(im), ctt)
                    )
                )
        else:
            cdata = [None] * len(data)
            for i in range(len(data)):
                cdata[i] = data[i]

        for N_step in range(len(xrng)):
            for Iter in range(max_iter):
                total_iter += 1
                tavg = sp_statistics.ave_series(cdata)
                if CTF:
                    tavg = sp_fundamentals.fft(
                        EMAN2_cppwrap.Util.divn_img(sp_fundamentals.fft(tavg), ctf2)
                    )
                if FH > 0.0:
                    fl = 0.1 + old_div((FH - 0.1) * Iter, float(max_iter - 1))
                    tavg = sp_filter.filt_tanl(tavg, fl, FF)
                """Multiline Comment49"""
                if total_iter == len(xrng) * max_iter:
                    if CTF:
                        for i in range(len(data)):
                            data[i].set_attr(
                                "xform.align2d", cdata[i].get_attr("xform.align2d")
                            )
                    return tavg
                cs = [0, 0]
                # print  "  iteration  std   %03d   %7.2f    %7.2f  "%(total_iter,cs[0],cs[1])
                if Iter % 4 != 0 or total_iter > max_iter * len(xrng) - 10:
                    delta = 0.0
                else:
                    delta = dst
                sx_sum, sy_sum, nope = sp_alignment.ali2d_single_iter(
                    cdata,
                    numr,
                    wr,
                    cs,
                    tavg,
                    cnx,
                    cny,
                    xrng[N_step],
                    yrng[N_step],
                    step[N_step],
                    mode=mode,
                    CTF=False,
                    delta=delta,
                )
                # for im in data:  print get_params2D(im)
                # print  "tata ",total_iter-1,Util.infomask(tavg,None,True)
                # tavg.write_image('tata.hdf',total_iter-1)

    return tavg


"""Multiline Comment50"""


##############################################################################################


def refinement_2d_local(data, ou, arange, xrng, yrng, CTF=True, SNR=1.0e10):
    """
		Note applicability of this code is limited to certain class of images.
		It performs well for single-particle "averages, i.e., centered images
		surrounded by zero background to ~20% of window size.
	"""

    #  There are two functions internal to the alignment
    def current_ave(mstack, SNR, fl=0.47, aa=0.01):
        """
		[
		0	original prepf of CTF-multiplied image, shift-eliminated, mirrored, but no (1-2*mir)*alpha applied  CONST
		1	[original params alpha_o, sx_o, sy_o, mir_o]  CONST
		2	CTF object  CONST
		3	[current params alpha_c, sx_c, sy_c, 0]  # we do not change mir so we keep it at zero
		4	FT of CTF-multiplied image, shift-eliminated, mirrored, with (1-2*mir)*alpha applied  CONST
		5	CTF^2 in the current combined orientation
		6	current ccc
		7	FT of current CTF-multiplied image, shift-eliminated, mirrored, in the current orientation, to be used for AVE
		]


		global settings
		current sum of images
		current sum of CTF^2
		"""
        nx = mstack[0][7].get_ysize()  # Input is Fourier file
        nima = len(mstack)
        tima = list(range(nima))
        random.shuffle(tima)
        mask = EMAN2_cppwrap.Util.unrollmask(nx, nx)

        sqct2 = []
        sqave = []
        parts = []
        for isub in range(2):
            ib = old_div(isub * nima, 2)
            ie = old_div((isub + 1) * nima, 2)
            qave = mstack[ib][7].copy()
            qct2 = mstack[ib][5].copy()
            for i in range(ib + 1, ie):
                EMAN2_cppwrap.Util.add_img(qave, mstack[i][7])
                EMAN2_cppwrap.Util.add_img(qct2, mstack[i][5])

            qct2 += old_div(1.0, SNR)
            sqave.append(qave)
            sqct2.append(qct2)
            diva = sp_utilities.model_blank(qct2.get_xsize(), qct2.get_ysize(), 1, 1.0)
            EMAN2_cppwrap.Util.div_img(diva, qct2)
            parts.append(EMAN2_cppwrap.Util.mulnclreal(qave, diva))

        # fff = fsc(parts[0],parts[1])
        ffm = sp_statistics.fsc(
            sp_morphology.adaptive_mask(sp_fundamentals.fft(parts[0])),
            sp_morphology.adaptive_mask(sp_fundamentals.fft(parts[1])),
        )
        # fl1 = -1.0
        fl2 = -1.0
        for i in range(len(ffm[1][1:])):
            # print("UUU:  %3d   %6.4f     %6.2f   %6.2f"%(i,fff[0][i],q,ffm[1][i]))
            # if(fl1<0.0):
            # 	if(q<0.1):  fl1 = fff[0][i-1]
            # if(fl2<0.0):
            if ffm[1][i] < 0.1:
                fl2 = ffm[0][i - 1]
                break

        if fl2 == -1.0:
            fl2 = fl
        # print("  FSC     %6.4f   %6.4f"%(fl1,fl2))

        qave = EMAN2_cppwrap.Util.addn_img(sqave[0], sqave[1])
        qct2 = EMAN2_cppwrap.Util.addn_img(sqct2[0], sqct2[1])
        diva = sp_utilities.model_blank(qct2.get_xsize(), qct2.get_ysize(), 1, 1.0)
        EMAN2_cppwrap.Util.div_img(diva, qct2)
        EMAN2_cppwrap.Util.mulclreal(qave, diva)
        if fl > 0.0:
            qave = sp_filter.filt_tanl(qave, fl2, aa)

        L2 = EMAN2_cppwrap.Util.innerproduct(qave, qave, None)
        return qave, L2, fl2

    def iteration_2Dshc(mstack, ave, shifts, angles):
        nx = mstack[0][4].get_ysize()  # Input is Fourier file

        nima = len(mstack)
        kl = 2 * shifts + 1
        # prepare rotated average in Fourier space
        fave = sp_fundamentals.prepf(ave)
        sfave = []
        for a in angles:
            sfave.append(
                sp_fundamentals.rot_shift2D(fave, a, interpolation_method="fourier")
            )
        del fave
        window_indexes = []
        for i in range(kl):
            for j in range(kl):
                window_indexes.append([i, j])

        npeak = 0
        got_better = False
        for i in range(nima):
            previous = mstack[i][6]
            tangle = list(range(len(angles)))
            random.shuffle(tangle)
            image_got_better = False
            for ia in tangle:
                if image_got_better:
                    break
                temp = sp_fundamentals.window2d(
                    sp_fundamentals.ccf(sfave[ia], mstack[i][4]), kl, kl
                )
                #  If no candidate here, do nothing
                if temp["maximum"] > mstack[i][6]:
                    random.shuffle(window_indexes)
                    for ti in window_indexes:
                        # print("SHC:  ",i,ia,-angles[ia],ix-shifts,iy-shifts,temp[ti[0],ti[1]])
                        if temp[ti[0], ti[1]] > mstack[i][6]:
                            image_got_better = True
                            mstack[i][6] = temp[ti[0], ti[1]]
                            six = ti[0] - shifts
                            siy = ti[1] - shifts
                            alphaf = -angles[ia]
                            break
                        npeak += 1

            if (
                image_got_better
            ):  # The reason to have this if here is that it makes it easy to switch between exhaust and SHC
                alphan, six, siy, _ = sp_utilities.combine_params2(
                    (1 - 2 * mstack[i][1][3]) * mstack[i][1][0],
                    0,
                    0,
                    0,
                    alphaf,
                    six,
                    siy,
                    0,
                )
                if (
                    (mstack[i][3][0] != alphan)
                    or (mstack[i][3][1] != six)
                    or (mstack[i][3][2] != siy)
                ):
                    mstack[i][3][0] = alphan
                    mstack[i][3][1] = six
                    mstack[i][3][2] = siy
                    mstack[i][7] = sp_fundamentals.rot_shift2D(
                        mstack[i][0],
                        mstack[i][3][0],
                        mstack[i][3][1],
                        mstack[i][3][2],
                        interpolation_method="fourier",
                    )
                    ctf_rot = EMAN2_cppwrap.EMAN2Ctf()
                    ctf_rot.copy_from(mstack[i][2])
                    if mir == 0:
                        ctf_rot.dfang = ctf_rot.dfang + mstack[i][1][0] + alphaf
                    else:
                        ctf_rot.dfang = 270.0 - ctf_rot.dfang - mstack[i][1][0] + alphaf
                    mstack[i][5] = sp_morphology.square(
                        sp_morphology.ctf_img_real(nx, ctf_rot)
                    )
                    got_better = True
            """Multiline Comment51"""

        return got_better

    nx = data[0].get_xsize()
    nima = len(data)
    cosine_mask = sp_morphology.cosinemask(
        sp_utilities.model_blank(nx, nx, 1, 1.0), ou, s=0.0
    )
    outside_mask = sp_utilities.model_blank(nx, nx, 1, 1.0) - cosine_mask
    angle_step = numpy.degrees(numpy.tan(old_div(1.0, ou)))

    mstack = []
    pwrot = []
    for i in range(nima):
        ima = data[i].copy()
        if CTF:
            ctf_params = ima.get_attr("ctf")
        else:
            ctf_params = sp_utilities.generate_ctf(
                [1.0e-5, 0.0, 300.0, 1.0, 0.0, 100.0, 0.0, 0]
            )  # fake, values approximately one.
        ctf_rot = EMAN2_cppwrap.EMAN2Ctf()
        ctf_rot.copy_from(ctf_params)
        ctc = sp_morphology.ctf_img_real(nx, ctf_params)
        alpha, sx, sy, mir, _ = sp_utilities.get_params2D(ima)
        # Eliminate shifts
        _, sx_c, sy_c, _ = sp_utilities.combine_params2(
            0.0, sx, sy, 0, -alpha, 0.0, 0.0, 0
        )
        ima = sp_fundamentals.fshift(ima, sx_c, sy_c)
        st = EMAN2_cppwrap.Util.infomask(ima, cosine_mask, True)
        ima -= st[0]
        ima = old_div(ima, st[1])
        pwrot.append(
            sp_fundamentals.rops_table(EMAN2_cppwrap.Util.muln_img(ima, outside_mask))
        )
        #  Now that the image is centered, normalize and apply mask
        EMAN2_cppwrap.Util.mul_img(ima, cosine_mask)
        ima = sp_fundamentals.fft(
            EMAN2_cppwrap.Util.mulnclreal(sp_fundamentals.fft(ima), ctc)
        )
        if mir:
            ima.process_inplace("xform.mirror", {"axis": "x"})

        qt = sp_fundamentals.prepf(ima)
        #  Now deal with parameters
        alpha_c = sx_c = sy_c = 0.0
        if True:
            # this randomizes input parameters to get bland starting point, not really necessary since we are using SHC.
            alpha_i = random.randint(-arange, arange) * angle_step
            sx_i = random.randint(-xrng, xrng)
            sy_i = random.randint(-yrng, yrng)
            _, sx_c, sy_c, _ = sp_utilities.combine_params2(
                0.0, sx_i, sy_i, 0, -alpha_i, 0, 0, 0
            )
            alpha_c = alpha_i
        alpha_e, sx_e, sy_e, _ = sp_utilities.combine_params2(
            (1 - 2 * mir) * alpha, 0, 0, 0, alpha_c, sx_c, sy_c, 0
        )
        if mir == 0:
            ctf_rot.dfang = ctf_rot.dfang + alpha + alpha_c
        else:
            ctf_rot.dfang = 270.0 - ctf_rot.dfang - alpha + alpha_c
        ctc = sp_morphology.square(sp_morphology.ctf_img_real(nx, ctf_rot))
        mstack.append(
            [
                qt,
                [alpha, sx, sy, mir],
                ctf_params,
                [alpha_e, sx_e, sy_e, 0],
                sp_fundamentals.rot_shift2D(
                    qt, (1 - 2 * mir) * alpha, interpolation_method="fourier"
                ),
                ctc,
                -1.0e23,
                sp_fundamentals.rot_shift2D(
                    qt, alpha_e, sx_e, sy_e, interpolation_method="fourier"
                ),
            ]
        )

    avp = [0.0] * (len(pwrot[0]))
    # svp = [0.0]*(len(pwrot[0]))
    for i in range(len(pwrot)):
        for j, q in enumerate(pwrot[i]):
            avp[j] += q

    avp[0] = avp[1]
    for i in range(1, len(avp)):
        avp[i] = numpy.sqrt(old_div(avp[0], avp[i]))
    avp[0] = avp[1]
    #  We may want to output the noise shaping function
    # write_text_file([avp],"asvp.txt")

    del pwrot

    qave, L2, fl2 = current_ave(mstack, SNR, 0.1)
    got_better = True

    shifts = max(xrng, yrng)
    angles = list(range(-arange, arange + 1, 1))
    for i in range(len(angles)):
        angles[i] *= angle_step
    aa = 0.02
    iter = 0
    while got_better:
        avem = EMAN2_cppwrap.Util.muln_img(
            sp_fundamentals.fft(sp_filter.filt_table(qave, avp)), cosine_mask
        )
        got_better = iteration_2Dshc(mstack, avem, shifts, angles)
        if got_better:
            qave, L2, q = current_ave(mstack, SNR, fl2, aa)
            fl2 = max(q, fl2)
            iter += 1
            # print("L2  ITERATION  %4d:"%iter,"   FL: %6.4f    %12.1f"%(fl2,L2))
            # fft(qave).write_image("qave%03d.hdf"%iter)
    fipar = []
    for i in range(nima):
        # get the shift applied before anything happened
        _, sx_c, sy_c, _ = sp_utilities.combine_params2(
            0.0, mstack[i][1][1], mstack[i][1][2], 0, -mstack[i][1][0], 0.0, 0.0, 0
        )
        #  combine them with mirror and subsequent alignment parameters
        alpha_e, sx_e, sy_e, mir_e = sp_utilities.combine_params2(
            0.0,
            sx_c,
            sy_c,
            mstack[i][1][3],
            mstack[i][3][0],
            mstack[i][3][1],
            mstack[i][3][2],
            0,
        )
        fipar.append([alpha_e, sx_e, sy_e, mir_e])
    fl = -1.0
    qave, _, _ = current_ave(mstack, SNR, fl, aa)
    return sp_fundamentals.fft(qave), fipar, fl2


def ali_vol_3(
    vol, refv, ang_scale, shift_scale, radius=None, discrepancy="ccc", mask=None
):
    # rotation and shift
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_alignment    import ali_vol_func
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_utilities    import model_circle, amoeba

    nx = refv.get_xsize()
    ny = refv.get_ysize()
    nz = refv.get_zsize()
    if mask is None:
        if radius != None:
            mask = sp_utilities.model_circle(radius, nx, ny, nz)
        else:
            mask = sp_utilities.model_circle(
                float(old_div(min(nx, ny, nz), 2) - 2), nx, ny, nz
            )

    scale = [ang_scale, ang_scale, ang_scale, shift_scale, shift_scale, shift_scale]
    data = [vol, refv, mask, discrepancy]
    new_params = [0.0] * 6
    opt_params, funval, niter = sp_utilities.amoeba(
        new_params, scale, sp_alignment.ali_vol_func, 1.0e-4, 1.0e-4, 500, data
    )
    return opt_params


def ali_vol_rotate(vol, refv, ang_scale, radius=None, discrepancy="ccc"):
    # rotation
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_alignment    import ali_vol_func_rotate
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_utilities    import get_image, model_circle, get_params3D, set_params3D
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_utilities    import amoeba, compose_transform3
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_fundamentals import rot_shift3D

    ref = sp_utilities.get_image(refv)
    nx = ref.get_xsize()
    ny = ref.get_ysize()
    nz = ref.get_zsize()
    if radius != None:
        mask = sp_utilities.model_circle(radius, nx, ny, nz)
    else:
        mask = sp_utilities.model_circle(
            float(old_div(min(nx, ny, nz), 2) - 2), nx, ny, nz
        )

    # names_params = ["phi", "theta", "psi", "s3x", "s3y", "s3z", "scale"]
    params = sp_utilities.get_params3D(ref)
    # print  " params of the reference volume",params
    ref = sp_fundamentals.rot_shift3D(
        ref, params[0], params[1], params[2], params[3], params[4], params[5], params[7]
    )

    e = sp_utilities.get_image(vol)
    params = sp_utilities.get_params3D(e)
    # e = rot_shift3D(e, params[0], params[1], params[2], params[3], params[4], params[5], params[7])
    # print  " input params ", params
    data = [e, ref, mask, params, discrepancy]
    new_params = [0.0, 0.0, 0.0]
    new_params = sp_utilities.amoeba(
        new_params,
        [ang_scale, ang_scale, ang_scale],
        sp_alignment.ali_vol_func_rotate,
        1.0e-4,
        1.0e-4,
        500,
        data,
    )
    cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale = sp_utilities.compose_transform3(
        params[0],
        params[1],
        params[2],
        params[3],
        params[4],
        params[5],
        params[7],
        new_params[0][0],
        new_params[0][1],
        new_params[0][2],
        0.0,
        0.0,
        0.0,
        1.0,
    )
    # print  " new params ", cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale, new_params[1]
    sp_utilities.set_params3D(e, [cphi, ctheta, cpsi, cs2x, cs2y, cs2z, 0, cscale])
    if type(vol) == type(""):
        pass  # IMPORTIMPORTIMPORT from ..libpy.sp_utilities import write_headers
        sp_utilities.write_headers(vol, [e], [0])
    else:
        return e


def ali_vol_shift(vol, refv, shift_scale, radius=None, discrepancy="ccc"):
    # shift
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_alignment    import ali_vol_func_shift
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_utilities    import get_image, model_circle, get_params3D, set_params3D
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_utilities    import amoeba, compose_transform3
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_fundamentals import rot_shift3D

    ref = sp_utilities.get_image(refv)
    nx = ref.get_xsize()
    ny = ref.get_ysize()
    nz = ref.get_zsize()
    if radius != None:
        mask = sp_utilities.model_circle(radius, nx, ny, nz)
    else:
        mask = sp_utilities.model_circle(
            float(old_div(min(nx, ny, nz), 2) - 2), nx, ny, nz
        )

    # names_params = ["phi", "theta", "psi", "s3x", "s3y", "s3z", "scale"]
    params = sp_utilities.get_params3D(ref)
    # print  " params of the reference volume",params
    ref = sp_fundamentals.rot_shift3D(
        ref, params[0], params[1], params[2], params[3], params[4], params[5], params[7]
    )

    e = sp_utilities.get_image(vol)
    params = sp_utilities.get_params3D(e)
    # e = rot_shift3D(e, params[0], params[1], params[2], params[3], params[4], params[5], params[7])
    # print  " input params ",params
    data = [e, ref, mask, params, discrepancy]
    new_params = [0.0, 0.0, 0.0]
    new_params = sp_utilities.amoeba(
        new_params,
        [shift_scale, shift_scale, shift_scale],
        sp_alignment.ali_vol_func_shift,
        1.0e-4,
        1.0e-4,
        500,
        data,
    )
    cphi, ctheta, cpsi, cs3x, cs3y, cs3z, cscale = sp_utilities.compose_transform3(
        params[0],
        params[1],
        params[2],
        params[3],
        params[4],
        params[5],
        params[7],
        0.0,
        0.0,
        0.0,
        new_params[0][0],
        new_params[0][1],
        new_params[0][2],
        1.0,
    )
    # print  " new params ", cphi, ctheta, cpsi, cs3x, cs3y, cs3z, cscale, new_params[1]
    sp_utilities.set_params3D(e, [cphi, ctheta, cpsi, cs3x, cs3y, cs3z, 0, cscale])
    if type(vol) == type(""):
        pass  # IMPORTIMPORTIMPORT from ..libpy.sp_utilities import write_headers
        sp_utilities.write_headers(vol, [e], [0])
    else:
        return e


def wrapper_params_2D_to_3D(stack):
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_utilities import params_2D_3D, print_begin_msg, print_end_msg, print_msg, get_params2D, set_params_proj, write_header

    # print_begin_msg("params_2D_to_3D")
    # print_msg("Input stack                 : %s\n\n"%(stack))

    nima = EMAN2.EMUtil.get_image_count(stack)
    ima = EMAN2.EMData()
    for im in range(nima):
        ima.read_image(stack, im, True)
        p = sp_utilities.get_params2D(ima)
        p = sp_utilities.params_2D_3D(p[0], p[1], p[2], int(p[3]))
        sp_utilities.set_params_proj(ima, p)
        sp_utilities.write_header(stack, ima, im)


# print_end_msg("params_2D_to_3D")

"""Multiline Comment61"""


from builtins import range
