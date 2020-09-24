
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
import mpi
import numpy
from . import sp_filter
from . import sp_fundamentals
from . import sp_global_def
from . import sp_morphology
from . import sp_utilities
import sys
from itertools import zip_longest
from future import standard_library

standard_library.install_aliases()


def add_ave_varf_MPI(
    myid,
    data,
    mask=None,
    mode="a",
    CTF=False,
    ctf_2_sum=None,
    ali_params="xform.align2d",
    main_node=0,
    comm=-1,
):
    """
		Calculate sum of an image series and sum of squares in Fourier space
		Since this is the MPI version, we need to reduce sum and sum of squares
		on the main node and calculate variance there.
		mode - "a": use current alignment parameters
		CTF  - if True, use CTF for calculations of the sum.
	"""

    if comm == -1:
        comm = mpi.MPI_COMM_WORLD

    n = len(data)
    nx = data[0].get_xsize()
    ny = data[0].get_ysize()
    ave1 = EMAN2_cppwrap.EMData(nx, ny, 1, False)
    ave2 = EMAN2_cppwrap.EMData(nx, ny, 1, False)
    var = EMAN2_cppwrap.EMData(nx, ny, 1, False)

    if CTF:
        if data[0].get_attr_default("ctf_applied", 1) == 1:
            sp_global_def.ERROR("data cannot be ctf-applied", "add_ave_varf_MPI", 1)
        if ctf_2_sum:
            get_ctf2 = False
        else:
            get_ctf2 = True
        if get_ctf2:
            ctf_2_sum = EMAN2_cppwrap.EMData(nx, ny, 1, False)
        ctf_params = data[i].get_attr("ctf")
        for i in range(n):
            if mode == "a":
                alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(
                    data[i], ali_params
                )
                ima = sp_fundamentals.rot_shift2D(
                    data[i], alpha, sx, sy, mirror, scale, "quadratic"
                )
                if mask:
                    EMAN2_cppwrap.Util.mul_img(ima, mask)
                sp_fundamentals.fftip(ima)
                ctf_params.dfang += alpha
                if mirror == 1:
                    ctf_params.dfang = 270.0 - ctf_params.dfang
            else:
                if mask:
                    ima = sp_fundamentals.fft(
                        EMAN2_cppwrap.Util.muln_img(data[i], mask)
                    )
                else:
                    ima = sp_fundamentals.fft(data[i])
            ima_filt = sp_filter.filt_ctf(ima, ctf_params, dopad=False)
            if i % 2 == 0:
                EMAN2_cppwrap.Util.add_img(ave1, ima_filt)
            else:
                EMAN2_cppwrap.Util.add_img(ave2, ima_filt)
            EMAN2_cppwrap.Util.add_img2(var, ima)
            if get_ctf2:
                EMAN2_cppwrap.Util.add_img2(
                    ctf_2_sum, sp_morphology.ctf_img(nx, ctf_params)
                )
    else:
        get_ctf2 = False
        for i in range(n):
            if mode == "a":
                alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(
                    data[i], ali_params
                )
                ima = sp_fundamentals.rot_shift2D(
                    data[i], alpha, sx, sy, mirror, scale, "quadratic"
                )
                if mask:
                    EMAN2_cppwrap.Util.mul_img(ima, mask)
                sp_fundamentals.fftip(ima)
            else:
                if mask:
                    ima = sp_fundamentals.fft(
                        EMAN2_cppwrap.Util.muln_img(data[i], mask)
                    )
                else:
                    ima = sp_fundamentals.fft(data[i])
            if i % 2 == 0:
                EMAN2_cppwrap.Util.add_img(ave1, ima)
            else:
                EMAN2_cppwrap.Util.add_img(ave2, ima)
            EMAN2_cppwrap.Util.add_img2(var, ima)
    sp_utilities.reduce_EMData_to_root(ave1, myid, main_node, comm)
    sp_utilities.reduce_EMData_to_root(ave2, myid, main_node, comm)
    sp_utilities.reduce_EMData_to_root(var, myid, main_node, comm)
    if get_ctf2:
        sp_utilities.reduce_EMData_to_root(ctf_2_sum, myid, main_node, comm)
    nima = n
    nima = mpi.mpi_reduce(nima, 1, mpi.MPI_INT, mpi.MPI_SUM, main_node, comm)
    if myid == main_node:
        nima = int(nima)
        sumsq = EMAN2_cppwrap.Util.addn_img(ave1, ave2)
        if CTF:
            tavg = EMAN2_cppwrap.Util.divn_img(sumsq, ctf_2_sum)
            EMAN2_cppwrap.Util.mul_img(sumsq, sumsq)
            EMAN2_cppwrap.Util.div_img(sumsq, ctf_2_sum)
        else:
            tavg = EMAN2_cppwrap.Util.mult_scalar(sumsq, old_div(1.0, float(nima)))
            EMAN2_cppwrap.Util.mul_img(sumsq, sumsq)
            EMAN2_cppwrap.Util.mul_scalar(sumsq, old_div(1.0, float(nima)))
        EMAN2_cppwrap.Util.sub_img(var, sumsq)
        EMAN2_cppwrap.Util.mul_scalar(var, old_div(1.0, float(nima - 1)))
        var.set_value_at(0, 0, 1.0)
        st = EMAN2_cppwrap.Util.infomask(var, None, True)
        if st[2] < 0.0:
            sp_global_def.ERROR("Negative variance!", "add_ave_varf_MPI", 1)
    else:
        tavg = EMAN2_cppwrap.EMData()
        sumsq = EMAN2_cppwrap.EMData()
    return tavg, ave1, ave2, var, sumsq


def sum_oe(
    data, mode="a", CTF=False, ctf_2_sum=None, ctf_eo_sum=False, return_params=False
):
    """
		Calculate average of an image series
		mode - "a": use current alignment parameters
		CTF  - if True, use CTF for calculations of the average.
		In addition, calculate odd and even sums, these are not divided by the ctf^2
		If ctf_2_sum not provided, sum of ctf^2 will be calculated and returned
		If ctf_eo_sum is True, then compute ctf^2 in odd and even form
		If return_params is True, then return ali2d.xform parameters
	"""
    if CTF:
        sp_global_def.ERROR("DEBUG_TEST", "sum_oe", 0)
    n = len(data)
    if return_params:
        params_list = [None] * n
    if CTF:
        origin_size = data[0].get_xsize()
        ave1 = EMAN2_cppwrap.EMData(origin_size, origin_size, 1, False)
        ave2 = EMAN2_cppwrap.EMData(origin_size, origin_size, 1, False)
        ctf_2_sumo = EMAN2_cppwrap.EMData(origin_size, origin_size, 1, False)
        ctf_2_sume = EMAN2_cppwrap.EMData(origin_size, origin_size, 1, False)

        if data[0].get_attr_default("ctf_applied", 1) == 1:
            sp_global_def.ERROR("data cannot be ctf-applied", "sum_oe", 1)
        if ctf_2_sum:
            get_ctf2 = False
        else:
            get_ctf2 = True
        for im in range(n):
            current_ctf = data[im].get_attr("ctf")
            if im == 0:
                myctf = copy.deepcopy(current_ctf)
                ctt = sp_morphology.ctf_img(origin_size, myctf)
            else:
                if not sp_utilities.same_ctf(current_ctf, myctf):
                    myctf = copy.deepcopy(current_ctf)
                    ctt = sp_morphology.ctf_img(origin_size, myctf)
            if mode == "a":
                alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(data[im])
                ima = sp_fundamentals.rot_shift2D(
                    data[im], alpha, sx, sy, mirror, scale, "quadratic"
                )
                if return_params:
                    params_list[im] = [alpha, sx, sy, mirror, scale]
            else:
                ima = data[im]
                if return_params:
                    alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(data[im])
                    params_list[im] = [alpha, sx, sy, mirror, scale]
            ima = sp_fundamentals.fft(ima)
            EMAN2_cppwrap.Util.mul_img(ima, ctt)
            if im % 2 == 0:
                EMAN2_cppwrap.Util.add_img(ave1, ima)
                EMAN2_cppwrap.Util.add_img2(ctf_2_sume, ctt)
            else:
                EMAN2_cppwrap.Util.add_img(ave2, ima)
                EMAN2_cppwrap.Util.add_img2(ctf_2_sumo, ctt)
    else:
        nx = data[0].get_xsize()
        ny = data[0].get_ysize()
        ave1 = sp_utilities.model_blank(nx, ny, 1)
        ave2 = sp_utilities.model_blank(nx, ny, 1)
        for im in range(n):
            if mode == "a":
                alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(data[im])
                ima = sp_fundamentals.rot_shift2D(
                    data[im], alpha, sx, sy, mirror, scale, "quadratic"
                )
                if return_params:
                    params_list[im] = [alpha, sx, sy, mirror, scale]
            else:
                ima = data[im]
                if return_params:
                    alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(data[im])
                    params_list[im] = [alpha, sx, sy, mirror, scale]
            if im % 2 == 0:
                EMAN2_cppwrap.Util.add_img(ave1, ima)
            else:
                EMAN2_cppwrap.Util.add_img(ave2, ima)

    if CTF:
        if get_ctf2:
            if not ctf_eo_sum:  # Old usage
                ctf_2_sum = EMAN2_cppwrap.Util.addn_img(ctf_2_sume, ctf_2_sumo)
                return sp_fundamentals.fft(ave1), sp_fundamentals.fft(ave2), ctf_2_sum
            else:  # return Fourier images
                if return_params:
                    return ave1, ave2, ctf_2_sume, ctf_2_sumo, params_list
                else:
                    return ave1, ave2, ctf_2_sume, ctf_2_sumo
        else:
            if return_params:
                return ave1, ave2, params_list
            else:
                return ave1, ave2
    else:
        if not return_params:
            return ave1, ave2
        else:
            return ave1, ave2, params_list


def ave_var(data, mode="a", listID=None):
    """
		Calculate average and variance of a 2D or 3D image series
		with optional application of orientation parameters
		data can be either in-core stack or a disk file
	"""
    if type(data) == type(""):
        n = EMAN2_cppwrap.EMUtil.get_image_count(data)
    else:
        n = len(data)
    if listID == None:
        listID = list(range(n))
    img = sp_utilities.get_im(data, 0)
    nx = img.get_xsize()
    ny = img.get_ysize()
    nz = img.get_zsize()
    if mode == "a":
        if nz > 1:
            ali_params = "xform.align3d"
        else:
            ali_params = "xform.align2d"

    ave = sp_utilities.model_blank(nx, ny, nz)
    var = sp_utilities.model_blank(nx, ny, nz)
    nlistID = len(listID)
    for i in range(nlistID):
        img = sp_utilities.get_im(data, listID[i])
        if mode == "a":
            if nz > 1:
                phi, theta, psi, s3x, s3y, s3z, mirror, scale = sp_utilities.get_params3D(
                    img
                )
                img = sp_fundamentals.rot_shift3D(
                    img, phi, theta, psi, s3x, s3y, s3z, scale
                )
            else:
                angle, sx, sy, mirror, scale = sp_utilities.get_params2D(img)
                img = sp_fundamentals.rot_shift2D(img, angle, sx, sy, mirror, scale)
        EMAN2_cppwrap.Util.add_img(ave, img)
        EMAN2_cppwrap.Util.add_img2(var, img)
    EMAN2_cppwrap.Util.mul_scalar(ave, old_div(1.0, float(nlistID)))

    return ave, old_div((var - ave * ave * nlistID), (nlistID - 1))


def ave_series(data, pave=True, mask=None):
    """
		Calculate average of a image series using current alignment parameters
		data - real space image series
	"""
    n = len(data)
    nx = data[0].get_xsize()
    ny = data[0].get_ysize()
    ave = sp_utilities.model_blank(nx, ny)
    for i in range(n):
        alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(data[i])
        temp = sp_fundamentals.rot_shift2D(data[i], alpha, sx, sy, mirror)
        EMAN2_cppwrap.Util.add_img(ave, temp)
    if mask:
        EMAN2_cppwrap.Util.mul_img(ave, mask)
    if pave:
        EMAN2_cppwrap.Util.mul_scalar(ave, old_div(1.0, float(n)))
    return ave


"""Multiline Comment0"""

"""Multiline Comment1"""


def aves_wiener(input_stack, mode="a", SNR=1.0, interpolation_method="linear"):
    """
		Apply alignment parameters, and calculate Wiener average using CTF info
		mode="a" will apply alignment parameters to the input image.
	"""

    if type(input_stack) == type(""):
        n = EMAN2_cppwrap.EMUtil.get_image_count(input_stack)
    else:
        n = len(input_stack)
    ima = sp_utilities.get_im(input_stack, 0)
    nx = ima.get_xsize()
    ny = ima.get_xsize()

    if ima.get_attr_default("ctf_applied", -2) > 0:
        sp_global_def.ERROR("data cannot be ctf-applied", "aves_wiener", 1)

    ave = sp_utilities.model_blank(nx, ny)
    ctf_2_sum = EMAN2_cppwrap.EMData(nx, ny, 1, False)
    snrsqrt = numpy.sqrt(SNR)

    for i in range(n):
        ima = sp_utilities.get_im(input_stack, i)
        ctf_params = ima.get_attr("ctf")
        oc = sp_filter.filt_ctf(ima, ctf_params, dopad=False)
        if mode == "a":
            alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(ima)
            oc = sp_fundamentals.rot_shift2D(
                oc, alpha, sx, sy, mirror, interpolation_method=interpolation_method
            )
            ctf_params.dfang += alpha
            if mirror == 1:
                ctf_params.dfang = 270.0 - ctf_params.dfang
        EMAN2_cppwrap.Util.mul_scalar(oc, SNR)
        EMAN2_cppwrap.Util.add_img(ave, oc)
        EMAN2_cppwrap.Util.add_img2(
            ctf_2_sum, snrsqrt * sp_morphology.ctf_img(nx, ctf_params, ny=ny, nz=1)
        )
    ctf_2_sum += 1.0
    ave = EMAN2_cppwrap.Util.divn_filter(
        sp_fundamentals.fft(ave), ctf_2_sum
    )  # result in Fourier space
    # variance
    var = EMAN2_cppwrap.EMData(nx, ny)
    var.to_zero()
    for i in range(n):
        ima = sp_utilities.get_im(input_stack, i)
        ctf_params = ima.get_attr("ctf")
        if mode == "a":
            alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(ima)
            ima = sp_fundamentals.rot_shift2D(
                ima, alpha, sx, sy, mirror, interpolation_method=interpolation_method
            )
            ctf_params.dfang += alpha
            if mirror == 1:
                ctf_params.dfang = 270.0 - ctf_params.dfang
        """Multiline Comment2"""
        #  Subtract in Fourier space, multiply again by the ctf, divide by the ctf^2, fft, square in real space
        ima = EMAN2_cppwrap.Util.divn_filter(
            sp_filter.filt_ctf(
                EMAN2_cppwrap.Util.subn_img(
                    sp_fundamentals.fft(ima),
                    sp_filter.filt_ctf(ave, ctf_params, dopad=False),
                ),
                ctf_params,
                dopad=False,
            ),
            ctf_2_sum,
        )
        EMAN2_cppwrap.Util.mul_scalar(ima, SNR)
        EMAN2_cppwrap.Util.add_img2(var, sp_fundamentals.fft(ima))
    # ave = Util.window(fft(ave),nx,ny,1,0,0,0)
    # Util.mul_scalar(var, 1.0/(n-1))
    return sp_fundamentals.fft(ave), var


"""Multiline Comment5"""


def varf2d_MPI(myid, data, ave, mask=None, mode="a", CTF=False, main_node=0, comm=-1):
    """
	Calculate variance in Fourier space for 2D images optionally including ctf correction
	ave is the Wiener average of data
	If mode = "a" apply alignment parameters
	This command is for Wiener average, i.e., A = sum_k (CTF_k SNR_k F_k) / [sum_k ( CTF_k^2 SNR_k) + 1]
	"""

    if comm == -1:
        comm = mpi.MPI_COMM_WORLD

    n = len(data)
    nx = data[0].get_xsize()
    ny = data[0].get_ysize()
    var = EMAN2_cppwrap.EMData(nx, ny, 1, False)

    if CTF:
        if data[0].get_attr_default("ctf_applied", 1) == 1:
            sp_global_def.ERROR("data cannot be ctf-applied", "add_ave_varf_MPI", 1)
        for i in range(n):
            ctf_params = data[i].get_attr("ctf")
            if mode == "a":
                alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(data[i])
                ima = sp_fundamentals.rot_shift2D(data[i], alpha, sx, sy, mirror)
                ctf_params.dfang += alpha
                if mirror == 1:
                    ctf_params.dfang = 270.0 - ctf_params.dfang
            else:
                ima = data[i].copy()
            if mask:
                EMAN2_cppwrap.Util.mul_img(ima, mask)
            oc = sp_filter.filt_ctf(ave, ctf_params, dopad=True)
            EMAN2_cppwrap.Util.add_img2(
                var, sp_fundamentals.fft(EMAN2_cppwrap.Util.subn_img(ima, oc))
            )
    else:
        for i in range(n):
            if mode == "a":
                alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(data[i])
                ima = sp_fundamentals.rot_shift2D(
                    data[i], alpha, sx, sy, mirror, scale, "quadratic"
                )
                if mask:
                    EMAN2_cppwrap.Util.mul_img(ima, mask)
                sp_fundamentals.fftip(ima)
            else:
                if mask:
                    ima = sp_fundamentals.fft(
                        EMAN2_cppwrap.Util.muln_img(data[i], mask)
                    )
                else:
                    ima = sp_fundamentals.fft(data[i])
            EMAN2_cppwrap.Util.add_img2(var, ima)
    sp_utilities.reduce_EMData_to_root(var, myid, main_node, comm)
    nima = n
    nima = mpi.mpi_reduce(nima, 1, mpi.MPI_INT, mpi.MPI_SUM, main_node, comm)
    if myid == main_node:
        EMAN2_cppwrap.Util.mul_scalar(var, old_div(1.0, float(nima - 1)))
        if var.get_value_at(0, 0) < 0.0:
            var.set_value_at(0, 0, 0.0)
        st = EMAN2_cppwrap.Util.infomask(var, None, True)
        if st[2] < 0.0:
            sp_global_def.ERROR("Negative variance!", "varf2_MPI", 1)
        return (
            var,
            sp_fundamentals.rot_avg_table(EMAN2_cppwrap.Util.pack_complex_to_real(var)),
        )
    else:
        return (
            EMAN2_cppwrap.EMData(),
            [0],
        )  # return minimum what has to be returned, but no meaning.


def ccc(img1, img2, mask=None):
    """Cross-correlation coefficient.
	   Usage: result = ccc(image1, image2 [, mask])
	"""
    return img1.cmp("ccc", img2, {"mask": mask, "negative": 0})


def fsc(img1, img2, w=1.0, filename=None):
    """Fourier Shell (or Ring) Correlation.

	   Usage: [frsc =] fsc(image1, image2 [, w, filename])

	   Computes the Fourier Shell (3d) or Fourier Ring (2d) correlation
	   function of two images.  If a filename is provided, then the
	   result is saved using that filename.
	"""
    result = img1.calc_fourier_shell_correlation(img2, w)
    # repack results as a list of (freq, fsc, n) triplets
    size = old_div(len(result), 3)
    frsc = []
    for i in range(3):
        frsc.append(result[i * size : (i + 1) * size])
    if filename:
        outf = open(filename, "w")
        for i in range(size):
            datstrings = []
            datstrings.append("  %12f" % (frsc[0][i]))
            datstrings.append("  %12f" % (frsc[1][i]))
            datstrings.append("  %12f" % (frsc[2][i]))
            datstrings.append("\n")
            outf.write("".join(datstrings))
        outf.close()
    return frsc


def fsc_mask(img1, img2, mask=None, w=1.0, filename=None):
    """
	        Compute fsc using mask and normalization.

		Usage: [frsc =] fsc(image1, image2 [, mask, w, filename])

		If no mask is provided, using circular mask with R=nx//2

	"""
    nx = img1.get_xsize()
    ny = img1.get_ysize()
    nz = img1.get_zsize()
    if mask == None:
        mask = sp_utilities.model_circle(old_div(nx, 2), nx, ny, nz)
    m = sp_morphology.binarize(mask, 0.5)
    s1 = EMAN2_cppwrap.Util.infomask(img1, m, True)
    s2 = EMAN2_cppwrap.Util.infomask(img2, m, True)
    return fsc((img1 - s1[0]) * mask, (img2 - s2[0]) * mask, w, filename)


def locres(vi, ui, m, nk, cutoff, step, myid, main_node, number_of_proc):
    nx = m.get_xsize()
    ny = m.get_ysize()
    nz = m.get_zsize()

    mc = sp_utilities.model_blank(nx, ny, nz, 1.0) - m

    if myid == main_node:
        st = EMAN2_cppwrap.Util.infomask(vi, m, True)
        vi -= st[0]

        st = EMAN2_cppwrap.Util.infomask(ui, m, True)
        ui -= st[1]

    sp_utilities.bcast_EMData_to_all(vi, myid, main_node)
    sp_utilities.bcast_EMData_to_all(ui, myid, main_node)

    vf = sp_fundamentals.fft(vi)
    uf = sp_fundamentals.fft(ui)

    if myid == 0:
        freqvol = sp_utilities.model_blank(nx, ny, nz)
        resolut = []
    lp = int(old_div(old_div(max(nx, ny, nz), 2), step) + 0.5)
    step = old_div(0.5, lp)
    lt = old_div(lp, number_of_proc)
    lp = (lt + 1) * number_of_proc
    bailout = 0
    for i in range(myid, lp, number_of_proc):
        fl = step * i
        fh = fl + step
        freq = old_div((fl + fh), 2.0)
        # print number_of_proc,myid,lp,i,step,fl,fh,freq

        if i > 0:
            v = sp_fundamentals.fft(sp_filter.filt_tophatb(vf, fl, fh))
            u = sp_fundamentals.fft(sp_filter.filt_tophatb(uf, fl, fh))
            tmp1 = EMAN2_cppwrap.Util.muln_img(v, v)
            tmp2 = EMAN2_cppwrap.Util.muln_img(u, u)
            tmp3 = EMAN2_cppwrap.Util.muln_img(u, v)
            do = EMAN2_cppwrap.Util.infomask(
                sp_morphology.square_root(
                    sp_morphology.threshold(EMAN2_cppwrap.Util.muln_img(tmp1, tmp2))
                ),
                m,
                True,
            )[0]
            dp = EMAN2_cppwrap.Util.infomask(tmp3, m, True)[0]
            # print "dpdo   ",myid,dp,do
            if do == 0.0:
                dis = [freq, 0.0]
            else:
                dis = [freq, old_div(dp, do)]
        else:
            tmp1 = sp_utilities.model_blank(nx, ny, nz, 1.0)
            tmp2 = sp_utilities.model_blank(nx, ny, nz, 1.0)
            tmp3 = sp_utilities.model_blank(nx, ny, nz, 1.0)
            dis = [freq, 1.0]

        tmp1 = EMAN2_cppwrap.Util.box_convolution(tmp1, nk)
        tmp2 = EMAN2_cppwrap.Util.box_convolution(tmp2, nk)
        tmp3 = EMAN2_cppwrap.Util.box_convolution(tmp3, nk)

        EMAN2_cppwrap.Util.mul_img(tmp1, tmp2)

        tmp1 = sp_morphology.square_root(sp_morphology.threshold(tmp1))

        EMAN2_cppwrap.Util.mul_img(tmp1, m)
        EMAN2_cppwrap.Util.add_img(tmp1, mc)

        EMAN2_cppwrap.Util.mul_img(tmp3, m)
        EMAN2_cppwrap.Util.add_img(tmp3, mc)

        EMAN2_cppwrap.Util.div_img(tmp3, tmp1)

        EMAN2_cppwrap.Util.mul_img(tmp3, m)

        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

        if myid == main_node:
            for k in range(number_of_proc):
                if k != main_node:
                    # print " start receiving",myid,i
                    tag_node = k + 1001
                    dis = mpi.mpi_recv(
                        2,
                        mpi.MPI_FLOAT,
                        k,
                        sp_global_def.SPARX_MPI_TAG_UNIVERSAL,
                        mpi.MPI_COMM_WORLD,
                    )
                    # print  "received ",myid,dis
                    tmp3 = sp_utilities.recv_EMData(k, tag_node)
                    # print  "received ",myid
                if dis[0] <= 0.5:
                    resolut.append(dis)
                fl = step * (i + k)
                fh = fl + step
                freq = old_div((fl + fh), 2.0)
                # print k,dis,Util.infomask(tmp3,m,True)
                # if(k == number_of_proc-1):  bailout = 1
                bailout = 0
                # print  "setting freqvol  ",k
                EMAN2_cppwrap.Util.set_freq_sphire(freqvol, tmp3, m, cutoff, freq)
                """Multiline Comment6"""
            """Multiline Comment7"""

        else:
            tag_node = myid + 1001
            # print   "sent from", myid,dis
            mpi.mpi_send(
                dis,
                2,
                mpi.MPI_FLOAT,
                main_node,
                sp_global_def.SPARX_MPI_TAG_UNIVERSAL,
                mpi.MPI_COMM_WORLD,
            )
            # print   "sending EMD from", myid
            sp_utilities.send_EMData(tmp3, main_node, tag_node)
            # print   "sent EMD from",myid

        bailout = sp_utilities.bcast_number_to_all(bailout, main_node)
        if bailout == 1:
            break

    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    if myid == main_node:
        return freqvol, resolut
    else:
        return None, None


"""Multiline Comment11"""


class Munkres(object):
    """
Calculate the Munkres solution to the classical assignment problem.
See the module documentation for usage.
    """

    def __init__(self):
        """Create a new instance"""
        self.C = None
        self.row_covered = []
        self.col_covered = []
        self.n = 0
        self.Z0_r = 0
        self.Z0_c = 0
        self.marked = None
        self.path = None

    def make_cost_matrix(profit_matrix, inversion_function):
        """
        Create a cost matrix from a profit matrix by calling
        'inversion_function' to invert each value. The inversion
        function must take one numeric argument (of any type) and return
        another numeric argument which is presumed to be the cost inverse
        of the original profit.

        This is a static method. Call it like this:

        cost_matrix = Munkres.make_cost_matrix(matrix, inversion_func)

        For example:

        cost_matrix = Munkres.make_cost_matrix(matrix, lambda x : sys.maxint - x)
        """
        cost_matrix = []
        for row in profit_matrix:
            cost_row = []
            for value in row:
                cost_row += [inversion_function(value)]
            cost_matrix += [cost_row]
        return cost_matrix

    make_cost_matrix = staticmethod(make_cost_matrix)

    def compute(self, cost_matrix):
        """
        Compute the indexes for the lowest-cost pairings between rows and
        columns in the database. Returns a list of (row, column) tuples
        that can be used to traverse the matrix.

        The matrix must be square.
        """
        self.C = self.__copy_matrix(cost_matrix)
        self.n = len(cost_matrix)
        self.row_covered = [False for i in range(self.n)]
        self.col_covered = [False for i in range(self.n)]
        self.Z0_r = 0
        self.Z0_c = 0
        self.path = self.__make_matrix(self.n * 2, 0)
        self.marked = self.__make_matrix(self.n, 0)

        done = False
        step = 1

        steps = {
            1: self.__step1,
            2: self.__step2,
            3: self.__step3,
            4: self.__step4,
            5: self.__step5,
            6: self.__step6,
        }

        while not done:
            try:
                func = steps[step]
                # print 'calling ' + str(func)
                step = func()
            except KeyError:
                done = True

        # Look for the starred columns
        results = []
        for i in range(self.n):
            for j in range(self.n):
                if self.marked[i][j] == 1:
                    results += [(i, j)]
        assert len(results) == self.n

        return results

    def __copy_matrix(self, matrix):
        """Return an exact copy of the supplied matrix"""
        copy = []
        for row in matrix:
            new_row = []
            for item in row:
                new_row += [item]
            copy += [new_row]
        return copy

    def __make_matrix(self, n, val):
        """Create an NxN matrix, populating it with the specific value."""
        matrix = []
        for i in range(n):
            matrix += [[val for j in range(n)]]
        return matrix

    def __step1(self):
        """
        For each row of the matrix, find the smallest element and
        subtract it from every element in its row. Go to Step 2.
        """
        C = self.C
        n = self.n
        for i in range(n):
            minval = self.C[i][0]
            # Find the minimum value for this row
            for j in range(n):
                if minval > self.C[i][j]:
                    minval = self.C[i][j]

            # Subtract that minimum from every element in the row.
            for j in range(n):
                self.C[i][j] -= minval

        return 2

    def __step2(self):
        """
        Find a zero (Z) in the resulting matrix. If there is no starred
        zero in its row or column, star Z. Repeat for each element in the
        matrix. Go to Step 3.
        """
        n = self.n
        for i in range(n):
            for j in range(n):
                if (
                    (self.C[i][j] == 0)
                    and (not self.col_covered[j])
                    and (not self.row_covered[i])
                ):
                    self.marked[i][j] = 1
                    self.col_covered[j] = True
                    self.row_covered[i] = True

        self.__clear_covers()
        return 3

    def __step3(self):
        """
        Cover each column containing a starred zero. If K columns are
        covered, the starred zeros describe a complete set of unique
        assignments. In this case, Go to DONE, otherwise, Go to Step 4.
        """
        n = self.n
        count = 0
        for i in range(n):
            for j in range(n):
                if self.marked[i][j] == 1:
                    self.col_covered[j] = True
                    count += 1

        if count >= n:
            step = 7  # done
        else:
            step = 4

        return step

    def __step4(self):
        """
        Find a noncovered zero and prime it. If there is no starred zero
        in the row containing this primed zero, Go to Step 5. Otherwise,
        cover this row and uncover the column containing the starred
        zero. Continue in this manner until there are no uncovered zeros
        left. Save the smallest uncovered value and Go to Step 6.
        """
        step = 0
        done = False
        row = -1
        col = -1
        star_col = -1
        while not done:
            (row, col) = self.__find_a_zero()
            if row < 0:
                done = True
                step = 6
            else:
                self.marked[row][col] = 2
                star_col = self.__find_star_in_row(row)
                if star_col >= 0:
                    col = star_col
                    self.row_covered[row] = True
                    self.col_covered[col] = False
                else:
                    done = True
                    self.Z0_r = row
                    self.Z0_c = col
                    step = 5

        return step

    def __step5(self):
        """
        Construct a series of alternating primed and starred zeros as
        follows. Let Z0 represent the uncovered primed zero found in Step 4.
        Let Z1 denote the starred zero in the column of Z0 (if any).
        Let Z2 denote the primed zero in the row of Z1 (there will always
        be one). Continue until the series terminates at a primed zero
        that has no starred zero in its column. Unstar each starred zero
        of the series, star each primed zero of the series, erase all
        primes and uncover every line in the matrix. Return to Step 3
        """
        count = 0
        path = self.path
        path[count][0] = self.Z0_r
        path[count][1] = self.Z0_c
        done = False
        while not done:
            row = self.__find_star_in_col(path[count][1])
            if row >= 0:
                count += 1
                path[count][0] = row
                path[count][1] = path[count - 1][1]
            else:
                done = True

            if not done:
                col = self.__find_prime_in_row(path[count][0])
                count += 1
                path[count][0] = path[count - 1][0]
                path[count][1] = col

        self.__convert_path(path, count)
        self.__clear_covers()
        self.__erase_primes()
        return 3

    def __step6(self):
        """
        Add the value found in Step 4 to every element of each covered
        row, and subtract it from every element of each uncovered column.
        Return to Step 4 without altering any stars, primes, or covered
        lines.
        """
        minval = self.__find_smallest()
        for i in range(self.n):
            for j in range(self.n):
                if self.row_covered[i]:
                    self.C[i][j] += minval
                if not self.col_covered[j]:
                    self.C[i][j] -= minval
        return 4

    def __find_smallest(self):
        """Find the smallest uncovered value in the matrix."""
        minval = sys.maxsize
        for i in range(self.n):
            for j in range(self.n):
                if (not self.row_covered[i]) and (not self.col_covered[j]):
                    if minval > self.C[i][j]:
                        minval = self.C[i][j]
        return minval

    def __find_a_zero(self):
        """Find the first uncovered element with value 0"""
        row = -1
        col = -1
        i = 0
        n = self.n
        done = False

        while not done:
            j = 0
            while True:
                if (
                    (self.C[i][j] == 0)
                    and (not self.row_covered[i])
                    and (not self.col_covered[j])
                ):
                    row = i
                    col = j
                    done = True
                j += 1
                if j >= n:
                    break
            i += 1
            if i >= n:
                done = True

        return (row, col)

    def __find_star_in_row(self, row):
        """
        Find the first starred element in the specified row. Returns
        the column index, or -1 if no starred element was found.
        """
        col = -1
        for j in range(self.n):
            if self.marked[row][j] == 1:
                col = j
                break

        return col

    def __find_star_in_col(self, col):
        """
        Find the first starred element in the specified row. Returns
        the row index, or -1 if no starred element was found.
        """
        row = -1
        for i in range(self.n):
            if self.marked[i][col] == 1:
                row = i
                break

        return row

    def __find_prime_in_row(self, row):
        """
        Find the first prime element in the specified row. Returns
        the column index, or -1 if no starred element was found.
        """
        col = -1
        for j in range(self.n):
            if self.marked[row][j] == 2:
                col = j
                break

        return col

    def __convert_path(self, path, count):
        for i in range(count + 1):
            if self.marked[path[i][0]][path[i][1]] == 1:
                self.marked[path[i][0]][path[i][1]] = 0
            else:
                self.marked[path[i][0]][path[i][1]] = 1

    def __clear_covers(self):
        """Clear all covered matrix cells"""
        for i in range(self.n):
            self.row_covered[i] = False
            self.col_covered[i] = False

    def __erase_primes(self):
        """Erase all prime markings"""
        for i in range(self.n):
            for j in range(self.n):
                if self.marked[i][j] == 2:
                    self.marked[i][j] = 0


# NEEDS TO BE FIX
"""Multiline Comment12"""


# K-means main stability stream command line


def k_means_match_clusters_asg_new(asg1, asg2, T=0):
    # asg1 and asg2 are numpy array

    K = len(asg1)
    MAT = [[0] * K for i in range(K)]
    cost_MAT = [[0] * K for i in range(K)]
    dummy = numpy.array([0], "int32")
    for k1 in range(K):
        for k2 in range(K):
            MAT[k1][k2] = EMAN2_cppwrap.Util.k_means_cont_table(
                asg1[k1], asg2[k2], dummy, asg1[k1].size, asg2[k2].size, 0
            )
            if MAT[k1][k2] <= T:
                MAT[k1][k2] = 0
    for i in range(K):
        for j in range(K):
            cost_MAT[i][j] = sys.maxsize - MAT[i][j]
    m = Munkres()
    indexes = m.compute(cost_MAT)
    newindexes = []
    list_stable = []
    nb_tot_objs = 0
    for r, c in indexes:
        cont = MAT[r][c]
        if cont <= T:
            list_stable.append(numpy.array([], "int32"))
            continue
        nb_tot_objs += cont
        objs = numpy.zeros(cont, "int32")
        dummy = EMAN2_cppwrap.Util.k_means_cont_table(
            asg1[r], asg2[c], objs, asg1[r].size, asg2[c].size, 1
        )
        list_stable.append(objs)
        newindexes.append([r, c])

    return newindexes, list_stable, nb_tot_objs


# Hierarchical stability between partitions given by k-means


"""Multiline Comment13"""

"""Multiline Comment14"""


def hist_list(data, nbin=-1, fminiu=None, fmaxiu=None):
    """
	  Calculate histogram of the list elements
	  nbin - number of bins, if not provided it will be set such that in average there is 10 elements per bin
	  fminiu - user provided minimum value for the histogram, it has to be smaller than the smallest element in data
	  fmaxiu - user provided maximum value for the histogram, it has to be larger than the largest element in data
	"""
    if nbin < 0:
        nbin = old_div(len(data), 10)
    fmaxi = max(data)
    fmini = min(data)

    if fmaxi == fmini:
        hist = [0] * nbin
        hist[0] = len(data)
        return [fmaxi] * nbin, hist
    if fminiu != None:
        if fminiu < fmini:
            fmini = fminiu
    if fmaxiu != None:
        if fmaxiu > fmaxi:
            fmaxi = fmaxiu

    binsize_i = old_div((fmaxi - fmini), float(nbin))
    start_i = fmini

    region = [None] * nbin
    hist = [None] * nbin
    for i in range(nbin):
        region[i] = start_i + i * binsize_i
        hist[i] = 0

    for d in data:
        i = min(int(old_div((d - start_i), binsize_i)), nbin - 1)
        hist[i] += 1

    return region, hist


def linreg(X, Y):
    """
	  Linear regression y=ax+b
	"""

    Sx = Sy = Sxx = Syy = Sxy = 0.0
    N = len(X)
    for x, y in zip_longest(X, Y):
        Sx += x
        Sy += y
        Sxx += x * x
        Syy += y * y
        Sxy += x * y
    det = Sxx * N - Sx * Sx
    a, b = old_div(Sxy * N - Sy * Sx, det), old_div(Sxx * Sy - Sx * Sxy, det)
    """Multiline Comment15"""
    return a, b


def pearson(X, Y):
    """
	  Pearson correlation coefficient between two lists
	"""
    Sx = Sy = Sxx = Syy = Sxy = 0.0
    N = len(X)
    # for x, y in map(None, X, Y):
    for x,y in zip_longest(X,Y):
        Sx += x
        Sy += y
        Sxx += x * x
        Syy += y * y
        Sxy += x * y
    return old_div(
        (Sxy - old_div(Sx * Sy, N)),
        numpy.sqrt((Sxx - old_div(Sx * Sx, N)) * (Syy - old_div(Sy * Sy, N))),
    )


def table_stat(X):
    """
	  Basic statistics of numbers stored in a list: average, variance, minimum, maximum
	"""
    N = len(X)
    if N == 1:
        return X[0], 0.0, X[0], X[0]
    av = X[0]
    va = X[0] * X[0]
    mi = X[0]
    ma = X[0]
    for i in range(1, N):
        av += X[i]
        va += X[i] * X[i]
        mi = min(mi, X[i])
        ma = max(ma, X[i])
    return old_div(av, N), old_div((va - old_div(av * av, N)), float(N - 1)), mi, ma


def k_means_stab_bbenum(
    PART,
    T=10,
    nguesses=5,
    J=50,
    max_branching=40,
    stmult=0.25,
    branchfunc=2,
    LIM=-1,
    doMPI_init=False,
    njobs=-1,
    do_mpi=False,
    K=-1,
    cdim=[],
    nstart=-1,
    nstop=-1,
    top_Matches=[],
):
    """

		Input:
			PART:   list of list of arrays. So PART[0] corresponds to the first partition, PART[1] the second partition and so on.
				Each partition is a list of arrays of distinct integers sorted in ascending order. The arrays corresponds
				to classes, so PART[i][j] is the j-th class of the i-th partition.

				Example of how to construct a partition:
					K    = EMUtil.get_image_count(stack_name)
					part0 = []
					for k in xrange(K):
						im  = get_im(stack_name, k)
						lid = im.get_attr('members')
						lid = array(lid, 'int32')
						lid.sort()
						part0.append(lid.copy())

			K:	Number of classes in the first partition.

			T:      An integer. It is the threshold such that the stable sets corresponding to the matches in the output have size larger than T.
				Specifically, if there are say four partitions, and the i-th match in the output, i.e., MATCH[i], is [2,12,1,5], then
				the 2nd class of the first partition, 12th class of the second, first class of the third, and fifth class of the
				fourth have at least T elements in common.

			nguesses: Not used anymore. I'm going to remove it.

			J:	An integer. In the branching matching algorithm we use, each step corresponds to choosing a match to add to the
				collection of matches we will eventually output.
				Since at each step there are many different possibilities and
				it is generally not feasible to branch on them all, we consider only J matches with the largest weights in choosing which matches to branch on.
				See branchfunc below for more details.
				Intuitively, a larger J should give better results.
				If J=1, the algorithm defaults to the greedy algorithm, i.e, at each step, we choose the match which has the largest cost and add it to the collection of matches to output.


			branchfunc: An integer. Determines which branching function should be used.

				    Roughly speaking, the algorithm we use builds up the collection of matches iteratively. During each step, we determine
				    the next match that should be added to the collection. Since there are many possibilities, the role of the branching function
				    is to determine both how many and which possibilities to explore.

				    The branching functions chooses the possibilities to explore from the
				    J matches with the largest weights.

				    There are two possible choices for branchfunc: 2 or 4

				    branchfunc = 2:  The J matches are considered in order according to weight, beginning with the match with the largest weight.
				    	             Always branches on the current match with the largest weight.
						     Branch on the match with the second largest weight only if it is infeasible with the match with the largest weight.
				    	             Similarly, for each subsequent match, only branch on it if it is infeasible with at least LIM of the matches which have already been chosen to branch on.


				    branchfunc = 4:  The J matches are considered in order according to weight, beginning with the match with the largest weight.
				    		     We branch based on distribution of the cost of the J largest matches.
				    		     Let stdev be the standard deviation of the weights of the J matches.
						     We always branch on the match with the largest weight.
						     For the remaining J-1 matches, we branch on those which are within stmult*stdev of the weight of match with the largest weight.


			LIM: 	An integer smaller than J. See explanation for branchfunc above.

			max_branching: This is an upper bound on the maximum number of branches to explore. See explanation for branchfunc above. This is to ensure we get a result in reasonable time. Intuitively, the larger max_branching is, the likelihood of getting a better result (i.e, a collection of matches with a large total weight) is increased.

			stmult: An integer. See explanation for branchfunc above.

		Output: MATCH, STB_PART, CT_s, CT_t, ST, st

			If there are exactly two partitions, i.e., len(PART) == 2, then the matching will be computed using Hungarian algorithm, and those matches for which the corresponding stable set has size greater than T will be output. If T==0, then the output matches will be optimal.

			(EDIT 10/11/11: If there are exactly two partitions, i.e., len(PART) == 2, then the matching will be computed using Hungarian algorithm such that ONLY matches with weights greater than T are considered during the course of the algorithm. This means the matching computed using the Hungarian algorithm will contain only matches with weight greater than T. See k_means_match_clusters_asg_new.)

			MATCH: A list of arrays. Each array has len(PART) elements, and the i-th element corresponds to a class in the i-th partition.
			       So MATCH[0][0] is an index into PART[0], MATCH[0][1] is an index into PART[1], etc.

			STB_PART: A list of lists. The stable set corresponding to the i-th MATCH, i.e., MATCH[i], is stored in STB_PART[MATCH[i][0]].

			CT_s:   A list of integers. The size of the stable set corresponding to  the i-th MATCH, i.e., MATCH[i], is stored in CT_s[MATCH[i][0]].
				The quality of the output collection of matches, i.e., MATCH, can be determined by summing CT_s. The larger the better.

			CT_t:	A list of integers. The number of elements in the union of the classes in the i-th match is stored in CT_t[MATCH[i][0]].

			st:     st = 100.0 * sum(CT_s) / float(sum(CT_t))

			ST:	ST[k] = 100.0 * CT_s[k] / float(CT_t[k])

	"""

    MATCH = []
    np = len(PART)
    # do MPI_init: compute pruned partitions and Njobs top matches and return
    # if doMPI_init:
    # 	newParts,topMatches,class_dim, num_found=k_means_match_bbenum(PART,T=T, J=J, nguesses=nguesses, DoMPI_init=True,Njobs=njobs)
    # 	return newParts, topMatches, class_dim,num_found

    # if do_mpi:
    # 	MATCH,cost= k_means_match_bbenum(PART, T=T, J=J, nguesses=nguesses, DoMPI=True, K=K, np=np,c_dim=cdim, N_start=nstart,N_stop=nstop,topMatches=top_Matches)
    # 	return MATCH,cost

    if np > 2:
        MATCH = k_means_match_bbenum(
            PART,
            T=T,
            J=J,
            max_branching=max_branching,
            stmult=stmult,
            branchfunc=branchfunc,
            LIM=LIM,
            nguesses=nguesses,
            DoMPI=False,
        )

    list_stb = []
    tot_n = 0
    if np == 2:
        # First make sure the two partitions have the same number of classes. If not, pad the one with less with junk.
        K1 = len(PART[0])
        K2 = len(PART[1])
        maxK = max(K1, K2)
        if K1 < maxK or K2 < maxK:
            # ffind garbage value to pad empty classes with
            garbage_value = 923456
            garbage_incr = 1
            for i in range(np):
                K = len(PART[i])
                for j in range(K):
                    pSize = PART[i][j].size
                    if pSize >= 1:
                        for p in range(pSize):
                            if PART[i][j][p] >= garbage_value:
                                garbage_value = PART[i][j][p] + garbage_incr
            for i in range(np):
                if len(PART[i]) < maxK:
                    # pad with empty arrays
                    df = maxK - len(PART[i])
                    for pd in range(df):
                        PART[i].append(
                            numpy.array([garbage_value, garbage_value + 1], "int32")
                        )
                        garbage_value = garbage_value + 2

        # now call
        MATCH, list_stb, tot_n = k_means_match_clusters_asg_new(PART[0], PART[1], T)

    K = len(PART[0])

    STB_PART = [[] for i in range(K)]
    nm = len(MATCH)
    CT_t = [0] * K
    CT_s = [0] * K
    ST = [0] * K
    ct_t = 0
    ct_s = 0
    st = 0

    for k in range(nm):
        kk = int(MATCH[k][0])  # due to numpy obj
        vmax = [0] * np
        vmin = [0] * np
        for i in range(np):
            vmax[i] = max(PART[i][int(MATCH[k][i])])
            vmin[i] = min(PART[i][int(MATCH[k][i])])

        vmax = int(max(vmax))
        vmin = int(min(vmin))
        vd = vmax - vmin + 1

        asg = [0] * vd
        for i in range(np):
            for item in PART[i][int(MATCH[k][i])]:
                asg[int(item) - vmin] += 1

        stb = []
        for i in range(vd):
            if asg[i] != 0:
                CT_t[kk] += 1
                if asg[i] == np:
                    CT_s[kk] += 1
                    stb.append(i + vmin)

        STB_PART[kk] = copy.deepcopy(stb)

    for k in range(K):
        if CT_t[k] == 0:
            continue
        ST[k] = old_div(100.0 * CT_s[k], float(CT_t[k]))

    if sum(CT_t) == 0:
        st = 0
    else:
        st = old_div(100.0 * sum(CT_s), float(sum(CT_t)))

    if np == 2:
        if sum(CT_s) != tot_n:
            sp_global_def.sxprint(sum(CT_s), tot_n)
            sp_global_def.sxprint("something wrong!!")
            sys.exit()

    return MATCH, STB_PART, CT_s, CT_t, ST, st


# DO NOT copy memory - could lead to crashes
# This is the wrapper function for bb_enumerateMPI. It packages the arguments and formats the output....


def k_means_match_bbenum(
    PART,
    T=10,
    J=1,
    max_branching=40,
    stmult=0.25,
    nguesses=5,
    branchfunc=2,
    LIM=-1,
    DoMPI_init=False,
    Njobs=-1,
    DoMPI=False,
    K=-1,
    np=-1,
    c_dim=[],
    N_start=-1,
    N_stop=-1,
    topMatches=[],
):
    MATCH = []
    output = []

    if DoMPI == True and DoMPI_init == False:
        sp_global_def.sxprint("Not supporting MPI currently")
        # if len(levels)<1:
        # 	for i in xrange(K):
        # 		levels.append(1)

        # ar_levels = array(levels, 'int32')
        # ar_class_dim	= array(c_dim,'int32')
        # ar_newParts = array(PART, 'int32')
        # firstmatches = topMatches[N_start*(np+1):(N_stop+1)*(np+1)]
        # output=Util.branchMPIpy(ar_newParts,ar_class_dim,np,K,T,ar_levels,K,nguesses,N_stop-N_start+1, array(firstmatches, 'int32'))
    else:
        LARGEST_CLASS = 0
        np = len(PART)

        # ar_argParts = array([],'int32')
        onedParts = []
        class_dim = []

        # figure out the garbage value to pad empty classes with

        garbage_value = 923456
        garbage_incr = 10
        for i in range(np):
            K = len(PART[i])
            for j in range(K):
                pSize = PART[i][j].size
                if pSize > LARGEST_CLASS:
                    LARGEST_CLASS = pSize
                if pSize >= 1:
                    for p in range(pSize):
                        if PART[i][j][p] >= garbage_value:
                            garbage_value = PART[i][j][p] + garbage_incr

        # deal with the case where not all the partitions have the same number of classes
        max_K = len(PART[0])
        for i in range(np):
            if len(PART[i]) > max_K:
                max_K = len(PART[i])
        for i in range(np):
            if not (len(PART[i]) == max_K):
                # pad with empty arrays
                df = max_K - len(PART[i])
                for j in range(df):
                    PART[i].append(
                        numpy.array([garbage_value, garbage_value + 1], "int32")
                    )
                    garbage_value = garbage_value + 2

        K = len(PART[0])
        # serialize all the arguments in preparation for passing into c++ function

        for i in range(np):
            for j in range(K):
                pSize = PART[i][j].size

                onedParts.append(j)
                onedParts.append(0)
                zero_pad = 0
                if pSize > 0:
                    for p in range(pSize):
                        onedParts.append(PART[i][j][p])
                else:
                    # pad with garbage value
                    zero_pad = 2
                    onedParts.append(garbage_value)
                    onedParts.append(garbage_value + 1)
                    garbage_value = garbage_value + zero_pad

                class_dim.append(pSize + zero_pad + 2)

                # ar_argParts = append(ar_argParts,[j,0])
                # ar_argParts=append(ar_argParts,PART[i][j])
        ar_argParts = numpy.array(onedParts, "int32")
        ar_class_dim = numpy.array(class_dim, "int32")

        # Single processor version
        output = EMAN2_cppwrap.Util.bb_enumerateMPI(
            ar_argParts,
            ar_class_dim,
            np,
            K,
            T,
            nguesses,
            LARGEST_CLASS,
            J,
            max_branching,
            stmult,
            branchfunc,
            LIM,
        )

    # first element of output is the total  cost of the solution, second element is the number of matches
    # in the output solution, and then follows the list of matches.

    # convert each match into an array and insert into MATCH
    num_matches = output[1]

    for j in range(num_matches):
        # get the j-th match
        ar_match = numpy.array(output[j * np + 2 : j * np + 2 + np], "int32")
        MATCH.append(ar_match)

    # order Matches in Match by group in first partition
    outMATCH = []
    for j in range(num_matches):
        amatch = []
        for js in range(len(MATCH[j])):
            amatch.append(MATCH[j][js])
        outMATCH.append(amatch)
    outMATCH.sort()

    if DoMPI:
        sp_global_def.sxprint("Not supporting MPI currently")
        return outMATCH, output[0]
    else:
        return outMATCH

def mono(k1,k2):
	"""
	get index of a square nxn matrix stored in a triangular form
	for i in xrange(1,n):
	    for j in xrange(i):
		print  i,j,mono(i,j)

	"""
	mk = max(k1,k2)
	return  min(k1,k2) + old_div(mk*(mk-1),2)
# match is a list, where every five tuple corresponds to a match


from builtins import range
from builtins import object
