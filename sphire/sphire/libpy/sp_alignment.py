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
import mpi
import numpy
from . import sp_applications
from . import sp_filter
from . import sp_fundamentals
from . import sp_global_def
from . import sp_pixel_error
from . import sp_projection
from . import sp_utilities
from . import sp_statistics
import types


def ali2d_single_iter(
    data,
    numr,
    wr,
    cs,
    tavg,
    cnx,
    cny,
    xrng,
    yrng,
    step,
    nomirror=False,
    mode="F",
    CTF=False,
    random_method="",
    T=1.0,
    ali_params="xform.align2d",
    delta=0.0,
):
    """
		single iteration of 2D alignment using ormq
		if CTF = True, apply CTF to data (not to reference!)
	"""

    maxrin = numr[-1]  #  length
    ou = numr[-3]  #  maximum radius
    if random_method == "SCF":
        frotim = [sp_fundamentals.fft(tavg)]
        xrng = int(xrng + 0.5)
        yrng = int(yrng + 0.5)
        cimage = EMAN2_cppwrap.Util.Polar2Dm(
            sp_fundamentals.scf(tavg), cnx, cny, numr, mode
        )

        EMAN2_cppwrap.Util.Frngs(cimage, numr)
        EMAN2_cppwrap.Util.Applyws(cimage, numr, wr)
    else:
        # 2D alignment using rotational ccf in polar coords and quadratic interpolation
        cimage = EMAN2_cppwrap.Util.Polar2Dm(tavg, cnx, cny, numr, mode)
        EMAN2_cppwrap.Util.Frngs(cimage, numr)
        EMAN2_cppwrap.Util.Applyws(cimage, numr, wr)

    sx_sum = 0.0
    sy_sum = 0.0
    sxn = 0.0
    syn = 0.0
    mn = 0
    nope = 0
    mashi = cnx - ou - 2
    for im in range(len(data)):
        if CTF:
            # Apply CTF to image
            ctf_params = data[im].get_attr("ctf")
            ima = sp_filter.filt_ctf(data[im], ctf_params, True)
        else:
            ima = data[im]

        if random_method == "PCP":
            sxi = data[im][0][0].get_attr("sxi")
            syi = data[im][0][0].get_attr("syi")
            nx = ny = data[im][0][0].get_attr("inx")
        else:
            nx = ima.get_xsize()
            ny = ima.get_ysize()
            alpha, sx, sy, mirror, dummy = sp_utilities.get_params2D(
                data[im], ali_params
            )
            alpha, sx, sy, dummy = sp_utilities.combine_params2(
                alpha, sx, sy, mirror, 0.0, -cs[0], -cs[1], 0
            )
            alphai, sxi, syi, scalei = sp_utilities.inverse_transform2(alpha, sx, sy)
            #  introduce constraints on parameters to accomodate use of cs centering
            sxi = min(max(sxi, -mashi), mashi)
            syi = min(max(syi, -mashi), mashi)

        #  The search range procedure was adjusted for 3D searches, so since in 2D the order of operations is inverted, we have to invert ranges
        txrng = search_range(nx, ou, sxi, xrng, "ali2d_single_iter")
        txrng = [txrng[1], txrng[0]]
        tyrng = search_range(ny, ou, syi, yrng, "ali2d_single_iter")
        tyrng = [tyrng[1], tyrng[0]]
        # print im, "B",cnx,sxi,syi,txrng, tyrng
        # align current image to the reference
        if random_method == "SHC":
            """Multiline Comment0"""
            #  For shc combining of shifts is problematic as the image may randomly slide away and never come back.
            #  A possibility would be to reject moves that results in too large departure from the center.
            #  On the other hand, one cannot simply do searches around the proper center all the time,
            #    as if xr is decreased, the image cannot be brought back if the established shifts are further than new range
            olo = EMAN2_cppwrap.Util.shc(
                ima,
                [cimage],
                txrng,
                tyrng,
                step,
                -1.0,
                mode,
                numr,
                cnx + sxi,
                cny + syi,
                "c1",
            )
            ##olo = Util.shc(ima, [cimage], xrng, yrng, step, -1.0, mode, numr, cnx, cny, "c1")
            if data[im].get_attr("previousmax") < olo[5]:
                # [angt, sxst, syst, mirrort, peakt] = ormq(ima, cimage, xrng, yrng, step, mode, numr, cnx+sxi, cny+syi, delta)
                # print  angt, sxst, syst, mirrort, peakt,olo
                angt = olo[0]
                sxst = olo[1]
                syst = olo[2]
                mirrort = int(olo[3])
                # combine parameters and set them to the header, ignore previous angle and mirror
                [alphan, sxn, syn, mn] = sp_utilities.combine_params2(
                    0.0, -sxi, -syi, 0, angt, sxst, syst, mirrort
                )
                sp_utilities.set_params2D(
                    data[im], [alphan, sxn, syn, mn, 1.0], ali_params
                )
                ##set_params2D(data[im], [angt, sxst, syst, mirrort, 1.0], ali_params)
                data[im].set_attr("previousmax", olo[5])
            else:
                # Did not find a better peak, but we have to set shifted parameters, as the average shifted
                sp_utilities.set_params2D(
                    data[im], [alpha, sx, sy, mirror, 1.0], ali_params
                )
                nope += 1
                mn = 0
                sxn = 0.0
                syn = 0.0
        elif random_method == "SCF":
            sxst, syst, iref, angt, mirrort, totpeak = multalign2d_scf(
                data[im], [cimage], frotim, numr, xrng, yrng, ou=ou
            )
            [alphan, sxn, syn, mn] = sp_utilities.combine_params2(
                0.0, -sxi, -syi, 0, angt, sxst, syst, mirrort
            )
            sp_utilities.set_params2D(data[im], [alphan, sxn, syn, mn, 1.0], ali_params)
        elif random_method == "PCP":
            [angt, sxst, syst, mirrort, peakt] = ormq_fast(
                data[im], cimage, txrng, tyrng, step, numr, mode, delta
            )
            sxst = rings[0][0][0].get_attr("sxi")
            syst = rings[0][0][0].get_attr("syi")
            sp_global_def.sxprint(sxst, syst, sx, sy)
            dummy, sxs, sys, dummy = sp_utilities.inverse_transform2(
                -angt, sx + sxst, sy + syst
            )
            sp_utilities.set_params2D(
                data[im][0][0], [angt, sxs, sys, mirrort, 1.0], ali_params
            )
        else:
            if nomirror:
                [angt, sxst, syst, mirrort, peakt] = ornq(
                    ima, cimage, txrng, tyrng, step, mode, numr, cnx + sxi, cny + syi
                )
            else:
                [angt, sxst, syst, mirrort, peakt] = ormq(
                    ima,
                    cimage,
                    txrng,
                    tyrng,
                    step,
                    mode,
                    numr,
                    cnx + sxi,
                    cny + syi,
                    delta,
                )
            # combine parameters and set them to the header, ignore previous angle and mirror
            [alphan, sxn, syn, mn] = sp_utilities.combine_params2(
                0.0, -sxi, -syi, 0, angt, sxst, syst, mirrort
            )
            sp_utilities.set_params2D(data[im], [alphan, sxn, syn, mn, 1.0], ali_params)

        if mn == 0:
            sx_sum += sxn
        else:
            sx_sum -= sxn
        sy_sum += syn

    return sx_sum, sy_sum, nope


"""Multiline Comment1"""


def ang_n(tot, mode, maxrin):
    """
	  Calculate angle based on the position of the peak
	"""
    if mode == "f" or mode == "F":
        return numpy.fmod((old_div((tot - 1.0), maxrin) + 1.0) * 360.0, 360.0)
    else:
        return numpy.fmod((old_div((tot - 1.0), maxrin) + 1.0) * 180.0, 180.0)


# Copy of this function is implemented in C++ in Util (Util.Applyws). It works much faster than this one.
"""Multiline Comment2"""


def crit2d(args, data):
    # print  " AMOEBA ",args
    #  data: 0 - kb,  1 - mask, 2 - nima,  3 - current ave, 4 - current image in the gridding format
    # from utilities import info
    mn = data[4].get_attr("mirror")
    temp = sp_fundamentals.rtshgkb(data[4], args[0], args[1], args[2], data[0])
    if mn:
        temp.process_inplace("xform.mirror", {"axis": "x"})
    # temp2 = data[3] + temp/data[2]
    temp2 = EMAN2_cppwrap.Util.madn_scalar(data[3], temp, old_div(1.0, data[2]))
    v = temp2.cmp("dot", temp2, {"negative": 0, "mask": data[1]})
    # print  " AMOEBA ",args,mn,v
    return v


def ringwe(numr, mode="F"):
    """
	   Calculate ring weights for rotational alignment
	   The weights are r*delta(r)*delta(phi).
	"""
    if mode == "f" or mode == "F":
        dpi = 2 * numpy.pi
    else:
        dpi = numpy.pi
    nring = old_div(len(numr), 3)
    wr = [0.0] * nring
    maxrin = float(numr[len(numr) - 1])
    for i in range(0, nring):
        wr[i] = old_div(numr[i * 3] * dpi, float(numr[2 + i * 3])) * old_div(
            maxrin, float(numr[2 + i * 3])
        )
    return wr


def ornq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi=0.0):
    """Determine shift and rotation between image and reference image (refim)
	   no mirror
		quadratic interpolation
		cnx, cny in FORTRAN convention
	"""
    # from utilities import info
    # print "ORNQ"
    peak = -1.0e23

    lkx = int(old_div(xrng[0], step))
    rkx = int(old_div(xrng[-1], step))

    lky = int(old_div(yrng[0], step))
    rky = int(old_div(yrng[-1], step))

    for i in range(-lky, rky + 1):
        iy = i * step
        for j in range(-lkx, rkx + 1):
            ix = j * step
            cimage = EMAN2_cppwrap.Util.Polar2Dm(image, cnx + ix, cny + iy, numr, mode)
            EMAN2_cppwrap.Util.Frngs(cimage, numr)
            retvals = EMAN2_cppwrap.Util.Crosrng_e(crefim, cimage, numr, 0, deltapsi)
            qn = retvals["qn"]
            if qn >= peak:
                sx = -ix
                sy = -iy
                ang = ang_n(retvals["tot"], mode, numr[-1])
                peak = qn
    # mirror is returned as zero for consistency
    mirror = 0
    co = numpy.cos(numpy.radians(ang))
    so = -numpy.sin(numpy.radians(ang))
    sxs = sx * co - sy * so
    sys = sx * so + sy * co
    return ang, sxs, sys, mirror, peak


def ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, delta=0.0):
    """Determine shift and rotation between image and reference image (crefim)
		crefim should be as FT of polar coords with applied weights
	        consider mirror
		quadratic interpolation
		cnx, cny in FORTRAN convention
	"""
    # print "ORMQ"
    peak = -1.0e23

    lkx = int(old_div(xrng[0], step))
    rkx = int(old_div(xrng[-1], step))

    lky = int(old_div(yrng[0], step))
    rky = int(old_div(yrng[-1], step))

    for i in range(-lky, rky + 1):
        iy = i * step
        for j in range(-lkx, rkx + 1):
            ix = j * step
            cimage = EMAN2_cppwrap.Util.Polar2Dm(image, cnx + ix, cny + iy, numr, mode)
            EMAN2_cppwrap.Util.Frngs(cimage, numr)
            # The following code it used when mirror is considered
            if delta == 0.0:
                retvals = EMAN2_cppwrap.Util.Crosrng_ms(crefim, cimage, numr, 0.0)
            else:
                retvals = EMAN2_cppwrap.Util.Crosrng_ms_delta(
                    crefim, cimage, numr, 0.0, delta
                )
            qn = retvals["qn"]
            qm = retvals["qm"]
            if qn >= peak or qm >= peak:
                sx = -ix
                sy = -iy
                if qn >= qm:
                    ang = ang_n(retvals["tot"], mode, numr[-1])
                    peak = qn
                    mirror = 0
                else:
                    ang = ang_n(retvals["tmt"], mode, numr[-1])
                    peak = qm
                    mirror = 1
            """Multiline Comment4"""
    co = numpy.cos(numpy.radians(ang))
    so = -numpy.sin(numpy.radians(ang))
    sxs = sx * co - sy * so
    sys = sx * so + sy * co
    return ang, sxs, sys, mirror, peak


def ormq_fast(dimage, crefim, xrng, yrng, step, numr, mode, delta=0.0):
    """Determine shift and rotation between image and reference image (crefim)
		crefim should be as FT of polar coords with applied weights
	        consider mirror
		cnx, cny in FORTRAN convention
	"""
    # from math import pi, cos, sin, radians
    # print "ORMQ_FAST"
    maxrange = old_div(len(dimage), 2)
    # istep = int(2*step)
    istep = int(step)
    """Multiline Comment5"""

    lkx = rkx = int(xrng * istep)

    lky = rky = int(yrng * istep)

    peak = -1.0e23
    for j in range(-lky, rky + 1, istep):
        for i in range(-lkx, rkx + 1, istep):
            if delta == 0.0:
                retvals = EMAN2_cppwrap.Util.Crosrng_ms(
                    crefim, dimage[i + maxrange][j + maxrange], numr, 0.0
                )
            else:
                retvals = EMAN2_cppwrap.Util.Crosrng_ms_delta(
                    crefim, dimage[i + maxrange][j + maxrange], numr, delta
                )
            qn = retvals["qn"]
            qm = retvals["qm"]
            if qn >= peak or qm >= peak:
                sx = i
                sy = j
                if qn >= qm:
                    ang = ang_n(retvals["tot"], mode, numr[-1])
                    peak = qn
                    mirror = 0
                else:
                    ang = ang_n(retvals["tmt"], mode, numr[-1])
                    peak = qm
                    mirror = 1
    """Multiline Comment6"""
    if peak < -1.0e20:
        sp_global_def.ERROR("ormq_fast", "failed, most likely due to search ranges", 1)
    # return  ang, sx/2.0, sy/2.0, mirror, peak
    return ang, sx, sy, mirror, peak


def prepref(data, maskfile, cnx, cny, numr, mode, maxrangex, maxrangey, step):
    # step = 1
    mashi = cnx - numr[-3] - 2
    nima = len(data)
    istep = int(old_div(1.0, step))
    dimage = [
        [
            [None for j in range(2 * maxrangey * istep + 1)]
            for i in range(2 * maxrangex * istep + 1)
        ]
        for im in range(nima)
    ]
    for im in range(nima):
        sts = EMAN2_cppwrap.Util.infomask(data[im], maskfile, False)
        data[im] -= sts[0]
        data[im] = old_div(data[im], sts[1])
        alpha, sx, sy, mirror, dummy = sp_utilities.get_params2D(data[im])
        # alpha, sx, sy, dummy         = combine_params2(alpha, sx, sy, mirror, 0.0, -cs[0], -cs[1], 0)
        alphai, sxi, syi, dummy = sp_utilities.combine_params2(
            0.0, sx, sy, 0, -alpha, 0, 0, 0
        )
        #  introduce constraints on parameters to accomodate use of cs centering
        sxi = min(max(sxi, -mashi), mashi)
        syi = min(max(syi, -mashi), mashi)
        for j in range(-maxrangey * istep, maxrangey * istep + 1):
            iy = j * step
            for i in range(-maxrangex * istep, maxrangex * istep + 1):
                ix = i * step
                dimage[im][i + maxrangex][j + maxrangey] = EMAN2_cppwrap.Util.Polar2Dm(
                    data[im], cnx + sxi + ix, cny + syi + iy, numr, mode
                )
                # print ' prepref  ',j,i,j+maxrangey,i+maxrangex
                EMAN2_cppwrap.Util.Frngs(dimage[im][i + maxrangex][j + maxrangey], numr)
        dimage[im][0][0].set_attr("sxi", sxi)
        dimage[im][0][0].set_attr("syi", syi)

    return dimage


"""Multiline Comment7"""


def prepare_refrings(
    volft,
    kb,
    nz=-1,
    delta=2.0,
    ref_a="P",
    sym="c1",
    numr=None,
    MPI=False,
    phiEqpsi="Zero",
    kbx=None,
    kby=None,
    initial_theta=None,
    delta_theta=None,
    initial_phi=None,
):
    """
		Generate quasi-evenly distributed reference projections converted to rings
		ref_a can be a list of angles, in which case it is used instead of being generated
	"""
    # mpi communicator can be sent by the MPI parameter
    if type(MPI) is bool:
        if MPI:
            mpi_comm = mpi.MPI_COMM_WORLD
    else:
        mpi_comm = MPI
        MPI = True

    mode = "F"

    if type(ref_a) is list:
        # if ref_a is  list, it has to be a list of projection directions, use it
        ref_angles = ref_a
    else:
        # generate list of Eulerian angles for reference projections
        #  phi, theta, psi
        if initial_theta and initial_phi:
            ref_angles = sp_utilities.even_angles(
                delta,
                theta1=initial_theta,
                phi1=initial_phi,
                symmetry=sym,
                method=ref_a,
                phiEqpsi=phiEqpsi,
            )
        else:
            if initial_theta is None:
                if sym[:1] == "c" or sym[:1] == "d":
                    ref_angles = sp_utilities.even_angles(
                        delta, symmetry=sym, method=ref_a, phiEqpsi=phiEqpsi
                    )
                else:
                    psp = sp_fundamentals.symclass(sym)
                    ref_angles = psp.even_angles(delta)
                    del psp
            else:
                if delta_theta is None:
                    delta_theta = 1.0
                ref_angles = sp_utilities.even_angles(
                    delta,
                    theta1=initial_theta,
                    theta2=delta_theta,
                    symmetry=sym,
                    method=ref_a,
                    phiEqpsi=phiEqpsi,
                )

    wr_four = ringwe(numr, mode)
    cnx = old_div(nz, 2) + 1
    cny = old_div(nz, 2) + 1
    num_ref = len(ref_angles)

    if MPI:
        myid = mpi.mpi_comm_rank(mpi_comm)
        ncpu = mpi.mpi_comm_size(mpi_comm)
    else:
        ncpu = 1
        myid = 0

    if nz < 1:
        sp_global_def.ERROR(
            "Data size has to be given (nz)", "prepare_refrings", 1, myid
        )

    ref_start, ref_end = sp_applications.MPI_start_end(num_ref, ncpu, myid)

    refrings = (
        []
    )  # list of (image objects) reference projections in Fourier representation

    sizex = numr[len(numr) - 2] + numr[len(numr) - 1] - 1

    for i in range(num_ref):
        prjref = EMAN2_cppwrap.EMData()
        prjref.set_size(sizex, 1, 1)
        refrings.append(prjref)

    if kbx is None:
        for i in range(ref_start, ref_end):
            prjref = sp_projection.prgs(
                volft,
                kb,
                [ref_angles[i][0], ref_angles[i][1], ref_angles[i][2], 0.0, 0.0],
            )
            cimage = EMAN2_cppwrap.Util.Polar2Dm(
                prjref, cnx, cny, numr, mode
            )  # currently set to quadratic....
            EMAN2_cppwrap.Util.Normalize_ring(cimage, numr, 0)
            EMAN2_cppwrap.Util.Frngs(cimage, numr)
            EMAN2_cppwrap.Util.Applyws(cimage, numr, wr_four)
            refrings[i] = cimage
    else:
        for i in range(ref_start, ref_end):
            prjref = sp_projection.prgs(
                volft,
                kb,
                [ref_angles[i][0], ref_angles[i][1], ref_angles[i][2], 0.0, 0.0],
                kbx,
                kby,
            )
            cimage = EMAN2_cppwrap.Util.Polar2Dm(
                prjref, cnx, cny, numr, mode
            )  # currently set to quadratic....
            EMAN2_cppwrap.Util.Normalize_ring(cimage, numr, 0)
            EMAN2_cppwrap.Util.Frngs(cimage, numr)
            EMAN2_cppwrap.Util.Applyws(cimage, numr, wr_four)
            refrings[i] = cimage

    if MPI:
        sp_utilities.bcast_compacted_EMData_all_to_all(refrings, myid, comm=mpi_comm)

    for i in range(len(ref_angles)):
        n1, n2, n3 = sp_utilities.getfvec(ref_angles[i][0], ref_angles[i][1])
        refrings[i].set_attr_dict(
            {
                "phi": ref_angles[i][0],
                "theta": ref_angles[i][1],
                "psi": ref_angles[i][2],
                "n1": n1,
                "n2": n2,
                "n3": n3,
            }
        )

    return refrings


#   Implemented in c in utli_sparx
#  Helper functions for ali2d_r
def kbt(nx, npad=2):
    # padd two times
    N = nx * npad
    # support of the window
    K = 6
    alpha = 1.75
    r = old_div(nx, 2)
    v = old_div(old_div(K, 2.0), N)
    return EMAN2_cppwrap.Util.KaiserBessel(alpha, K, r, old_div(K, (2.0 * N)), N)


def log2(n):
    """ Returns the smallest power by which 2 has to be raised to obtain
	    an integer less equal n
	"""
    m = 1
    k = -1
    while m <= n:
        i = m
        k += 1
        m = 2 * i
    return k


def Numrinit(first_ring, last_ring, skip=1, mode="F"):
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
    pass  # IMPORTIMPORTIMPORT from math import pi

    if mode == "f" or mode == "F":
        dpi = 2 * numpy.pi
    else:
        dpi = numpy.pi
    numr = []
    lcirc = 1
    for k in range(first_ring, last_ring + 1, skip):
        numr.append(k)
        jp = int(dpi * k + 0.5)
        ip = 2 ** (log2(jp) + 1)  # two times oversample each ring
        if k + skip <= last_ring and jp > ip + old_div(ip, 2):
            ip = min(MAXFFT, 2 * ip)
        if k + skip > last_ring and jp > ip + old_div(ip, 5):
            ip = min(MAXFFT, 2 * ip)

        numr.append(lcirc)
        numr.append(ip)
        lcirc += ip

    return numr


"""Multiline Comment9"""


def ali_vol_func(params, data):
    # print  params
    # print  data[3]
    # cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale= compose_transform3(data[3][0], data[3][1], data[3][2], data[3][3], data[3][4], data[3][5], data[3][6], params[0], params[1], params[2],params[3], params[4], params[5],1.0)
    # print  cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale
    x = sp_fundamentals.rot_shift3D(
        data[0], params[0], params[1], params[2], params[3], params[4], params[5], 1.0
    )

    res = -x.cmp("ccc", data[1], {"mask": data[2]})
    # print  " %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f  %10.5f" %(params[0], params[1], params[2],params[3], params[4], params[5], -res)
    return res


def fine_2D_refinement(data, br, mask, tavg, group=-1):

    # IMAGES ARE SQUARES!
    nx = data[0].get_xsize()
    #  center is in SPIDER convention
    cnx = int(old_div(nx, 2)) + 1
    cny = cnx

    if group > -1:
        nima = 0
        for im in range(len(data)):
            if data[im].get_attr("ref_num") == group:
                nima += 1
    else:
        nima = len(data)

    # prepare KB interpolants
    kb = kbt(nx)
    # load stuff for amoeba
    stuff = []
    stuff.insert(0, kb)
    stuff.insert(1, mask)
    stuff.insert(2, nima)
    # stuff.insert(3,tave)  # current average
    # stuff.insert(4,data)  # current image in the gridding format
    weights = [
        br
    ] * 3  # weights define initial bracketing, so one would have to figure how to set them correctly

    for im in range(len(data)):
        if group > -1:
            if data[im].get_attr("ref_num") != group:
                continue
        # subtract current image from the average
        alpha = data[im].get_attr("alpha")
        sx = data[im].get_attr("sx")
        sy = data[im].get_attr("sy")
        mirror = data[im].get_attr("mirror")
        ddata = sp_fundamentals.prepg(data[im], kb)
        ddata.set_attr_dict({"alpha": alpha, "sx": sx, "sy": sy, "mirror": mirror})
        temp = sp_fundamentals.rtshgkb(ddata, alpha, sx, sy, kb)
        if mirror:
            temp.process_inplace("xform.mirror", {"axis": "x"})
        #  Subtract current image from the average
        refim = EMAN2_cppwrap.Util.madn_scalar(tavg, temp, old_div(-1.0, float(nima)))
        stuff.append(refim)  # curent ave-1
        stuff.append(ddata)  # curent image
        # perform amoeba alignment
        params = [alpha, sx, sy]
        outparams = sp_utilities.amoeba(
            params, weights, crit2d, 1.0e-4, 1.0e-4, 500, stuff
        )
        del stuff[3]
        del stuff[3]
        # set parameters to the header
        data[im].set_attr_dict(
            {
                "alpha": outparams[0][0],
                "sx": outparams[0][1],
                "sy": outparams[0][2],
                "mirror": mirror,
            }
        )
        # update the average
        temp = sp_fundamentals.rtshgkb(
            ddata, outparams[0][0], outparams[0][1], outparams[0][2], kb
        )
        if mirror:
            temp.process_inplace("xform.mirror", {"axis": "x"})
        # check whether the criterion actually increased
        # add current image to the average
        tavg = EMAN2_cppwrap.Util.madn_scalar(refim, temp, old_div(1.0, float(nima)))
        # print  im,tave.cmp("dot", tave, {"negative":0,"mask":mask}),params,outparams[0],outparams[2]
        # tave,tvar = ave_var_series_g(data,kb)
        # print  " Criterium on the fly ", tave.cmp("dot", tave, {"negative":0,"mask":mask})


def align2d(
    image,
    refim,
    xrng=[0, 0],
    yrng=[0, 0],
    step=1.0,
    first_ring=1,
    last_ring=0,
    rstep=1,
    mode="F",
):
    """  Determine shift and rotation between image and reference image
	     quadratic interpolation
	     xrng[k,m] - translation search will be performed in range [-k,m] in steps equal step, which can non-integer
	     Output: ang, sxs, sys, mirror, peak
	"""
    # from utilities import print_col

    step = float(step)
    nx = refim.get_xsize()
    ny = refim.get_ysize()
    if last_ring == 0:
        last_ring = old_div(nx, 2) - 2 - int(max(max(xrng), max(yrng)))
    # center in SPIDER convention
    cnx = old_div(nx, 2) + 1
    cny = old_div(ny, 2) + 1
    # precalculate rings
    numr = Numrinit(first_ring, last_ring, rstep, mode)
    wr = ringwe(numr, mode)
    # cimage=Util.Polar2Dmi(refim, cnx, cny, numr, mode, kb)
    crefim = EMAN2_cppwrap.Util.Polar2Dm(refim, cnx, cny, numr, mode)
    # crefim = Util.Polar2D(refim, numr, mode)
    # print_col(crefim)
    EMAN2_cppwrap.Util.Frngs(crefim, numr)
    EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)
    return ormq(image, crefim, xrng, yrng, step, mode, numr, cnx, cny)


"""Multiline Comment14"""


def align2d_scf(image, refim, xrng=-1, yrng=-1, ou=-1):
    nx = image.get_xsize()
    ny = image.get_xsize()
    if ou < 0:
        ou = min(old_div(nx, 2) - 1, old_div(ny, 2) - 1)
    if yrng < 0:
        yrng = xrng
    if ou < 2:
        sp_global_def.ERROR(
            "Radius of the object (ou) has to be given", "align2d_scf", 1
        )
    sci = sp_fundamentals.scf(image)
    scr = sp_fundamentals.scf(refim)
    first_ring = 1

    # alpha1, sxs, sys, mirr, peak1 = align2d_no_mirror(scf(image), scr, last_ring=ou, mode="H")
    # alpha2, sxs, sys, mirr, peak2 = align2d_no_mirror(scf(mirror(image)), scr, last_ring=ou, mode="H")
    # alpha1, sxs, sys, mirr, peak1 = align2d_no_mirror(sci, scr, first_ring = 1, last_ring=ou, mode="H")
    # alpha2, sxs, sys, mirr, peak2 = align2d_no_mirror(mirror(sci), scr,  first_ring = 1, last_ring=ou, mode="H")

    # center in SPIDER convention
    cnx = old_div(nx, 2) + 1
    cny = old_div(ny, 2) + 1
    # precalculate rings
    numr = Numrinit(first_ring, ou, 1, "H")
    wr = ringwe(numr, "H")
    crefim = EMAN2_cppwrap.Util.Polar2Dm(scr, cnx, cny, numr, "H")
    EMAN2_cppwrap.Util.Frngs(crefim, numr)
    EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)
    alpha1, sxs, sys, mirr, peak1 = ornq(
        sci, crefim, [0.0], [0.0], 1.0, "H", numr, cnx, cny
    )
    alpha2, sxs, sys, mirr, peak2 = ornq(
        sp_fundamentals.mirror(sci), crefim, [0.0], [0.0], 1.0, "H", numr, cnx, cny
    )

    if peak1 > peak2:
        mirr = 0
        alpha = alpha1
    else:
        mirr = 1
        alpha = -alpha2
    nrx = min(2 * (xrng + 1) + 1, ((old_div((nx - 2), 2)) * 2 + 1))
    nry = min(2 * (yrng + 1) + 1, ((old_div((ny - 2), 2)) * 2 + 1))
    frotim = sp_fundamentals.fft(refim)
    ccf1 = EMAN2_cppwrap.Util.window(
        sp_fundamentals.ccf(
            sp_fundamentals.rot_shift2D(image, alpha, 0.0, 0.0, mirr), frotim
        ),
        nrx,
        nry,
    )
    p1 = sp_utilities.peak_search(ccf1)

    ccf2 = EMAN2_cppwrap.Util.window(
        sp_fundamentals.ccf(
            sp_fundamentals.rot_shift2D(image, alpha + 180.0, 0.0, 0.0, mirr), frotim
        ),
        nrx,
        nry,
    )
    p2 = sp_utilities.peak_search(ccf2)
    # print p1
    # print p2

    peak_val1 = p1[0][0]
    peak_val2 = p2[0][0]

    if peak_val1 > peak_val2:
        sxs = -p1[0][4]
        sys = -p1[0][5]
        cx = int(p1[0][1])
        cy = int(p1[0][2])
        peak = peak_val1
    else:
        alpha += 180.0
        sxs = -p2[0][4]
        sys = -p2[0][5]
        peak = peak_val2
        cx = int(p2[0][1])
        cy = int(p2[0][2])
        ccf1 = ccf2
    # print cx,cy
    z = sp_utilities.model_blank(3, 3)
    for i in range(3):
        for j in range(3):
            z[i, j] = ccf1[i + cx - 1, j + cy - 1]
    # print  ccf1[cx,cy],z[1,1]
    XSH, YSH, PEAKV = parabl(z)
    # print sxs, sys, XSH, YSH, PEAKV, peak
    if mirr == 1:
        sx = -sxs + XSH
    else:
        sx = sxs - XSH
    return alpha, sx, sys - YSH, mirr, PEAKV


def multalign2d_scf(image, refrings, frotim, numr, xrng=-1, yrng=-1, ou=-1):

    nx = image.get_xsize()
    ny = image.get_xsize()
    if ou < 0:
        ou = min(old_div(nx, 2) - 1, old_div(ny, 2) - 1)
    if yrng < 0:
        yrng = xrng
    if ou < 2:
        sp_global_def.ERROR(
            "Radius of the object (ou) has to be given", "align2d_scf", 1
        )
    sci = sp_fundamentals.scf(image)
    first_ring = 1
    # center in SPIDER convention
    cnx = old_div(nx, 2) + 1
    cny = old_div(ny, 2) + 1

    cimage = EMAN2_cppwrap.Util.Polar2Dm(sci, cnx, cny, numr, "H")
    EMAN2_cppwrap.Util.Frngs(cimage, numr)
    mimage = EMAN2_cppwrap.Util.Polar2Dm(
        sp_fundamentals.mirror(sci), cnx, cny, numr, "H"
    )
    EMAN2_cppwrap.Util.Frngs(mimage, numr)

    nrx = min(2 * (xrng + 1) + 1, (old_div(nx - 2, 2) * 2 + 1))
    nry = min(2 * (yrng + 1) + 1, (old_div(ny - 2, 2) * 2 + 1))

    totpeak = -1.0e23

    for iki in range(len(refrings)):
        # print  "TEMPLATE  ",iki
        #  Find angle
        retvals = EMAN2_cppwrap.Util.Crosrng_e(refrings[iki], cimage, numr, 0, 0.0)
        alpha1 = ang_n(retvals["tot"], "H", numr[-1])
        peak1 = retvals["qn"]
        retvals = EMAN2_cppwrap.Util.Crosrng_e(refrings[iki], mimage, numr, 0, 0.0)
        alpha2 = ang_n(retvals["tot"], "H", numr[-1])
        peak2 = retvals["qn"]
        # print  alpha1, peak1
        # print  alpha2, peak2

        if peak1 > peak2:
            mirr = 0
            alpha = alpha1
        else:
            mirr = 1
            alpha = -alpha2

        ccf1 = EMAN2_cppwrap.Util.window(
            sp_fundamentals.ccf(
                sp_fundamentals.rot_shift2D(image, alpha, 0.0, 0.0, mirr), frotim[iki]
            ),
            nrx,
            nry,
        )
        p1 = sp_utilities.peak_search(ccf1)

        ccf2 = EMAN2_cppwrap.Util.window(
            sp_fundamentals.ccf(
                sp_fundamentals.rot_shift2D(image, alpha + 180.0, 0.0, 0.0, mirr),
                frotim[iki],
            ),
            nrx,
            nry,
        )
        p2 = sp_utilities.peak_search(ccf2)



        peak_val1 = p1[0][0]
        peak_val2 = p2[0][0]

        if peak_val1 > peak_val2:
            sxs = -p1[0][4]
            sys = -p1[0][5]
            cx = int(p1[0][1])
            cy = int(p1[0][2])
            peak = peak_val1
        else:
            alpha += 180.0
            sxs = -p2[0][4]
            sys = -p2[0][5]
            peak = peak_val2
            cx = int(p2[0][1])
            cy = int(p2[0][2])
            ccf1 = ccf2
        # print cx,cy
        z = sp_utilities.model_blank(3, 3)
        for i in range(3):
            for j in range(3):
                z[i, j] = ccf1[i + cx - 1, j + cy - 1]
        XSH, YSH, PEAKV = parabl(z)
        # print  PEAKV
        if PEAKV > totpeak:
            totpeak = PEAKV
            iref = iki
            if mirr == 1:
                sx = -sxs + XSH
            else:
                sx = sxs - XSH
            sy = sys - YSH
            talpha = alpha
            tmirr = mirr
            # print "BETTER",sx,sy,iref,talpha,tmirr,totpeak
            # return alpha, sx, sys-YSH, mirr, PEAKV
    return sx, sy, iref, talpha, tmirr, totpeak


def parabl(Z):
    #  parabolic fit to a peak, C indexing
    C1 = old_div(
        (
            26.0 * Z[0, 0]
            - Z[0, 1]
            + 2 * Z[0, 2]
            - Z[1, 0]
            - 19.0 * Z[1, 1]
            - 7.0 * Z[1, 2]
            + 2.0 * Z[2, 0]
            - 7.0 * Z[2, 1]
            + 14.0 * Z[2, 2]
        ),
        9.0,
    )

    C2 = old_div(
        (
            8.0 * Z[0, 0]
            - 8.0 * Z[0, 1]
            + 5.0 * Z[1, 0]
            - 8.0 * Z[1, 1]
            + 3.0 * Z[1, 2]
            + 2.0 * Z[2, 0]
            - 8.0 * Z[2, 1]
            + 6.0 * Z[2, 2]
        ),
        (-6.0),
    )

    C3 = old_div(
        (
            Z[0, 0]
            - 2.0 * Z[0, 1]
            + Z[0, 2]
            + Z[1, 0]
            - 2.0 * Z[1, 1]
            + Z[1, 2]
            + Z[2, 0]
            - 2.0 * Z[2, 1]
            + Z[2, 2]
        ),
        6.0,
    )

    C4 = old_div(
        (
            8.0 * Z[0, 0]
            + 5.0 * Z[0, 1]
            + 2.0 * Z[0, 2]
            - 8.0 * Z[1, 0]
            - 8.0 * Z[1, 1]
            - 8.0 * Z[1, 2]
            + 3.0 * Z[2, 1]
            + 6.0 * Z[2, 2]
        ),
        (-6.0),
    )

    C5 = old_div((Z[0, 0] - Z[0, 2] - Z[2, 0] + Z[2, 2]), 4.0)

    C6 = old_div(
        (
            Z[0, 0]
            + Z[0, 1]
            + Z[0, 2]
            - 2.0 * Z[1, 0]
            - 2.0 * Z[1, 1]
            - 2.0 * Z[1, 2]
            + Z[2, 0]
            + Z[2, 1]
            + Z[2, 2]
        ),
        6.0,
    )

    DENOM = 4.0 * C3 * C6 - C5 * C5
    if DENOM == 0.0:
        return 0.0, 0.0, 0.0

    YSH = old_div((C4 * C5 - 2.0 * C2 * C6), DENOM ) - 2.0
    XSH = old_div((C2 * C5 - 2.0 * C4 * C3), DENOM ) - 2.0

    PEAKV = (
        4.0 * C1 * C3 * C6 - C1 * C5 * C5 - C2 * C2 * C6 + C2 * C4 * C5 - C4 * C4 * C3
    )
    PEAKV = old_div(PEAKV, DENOM)
    # print  "  in PARABL  ",XSH,YSH,Z[1,1],PEAKV

    XSH = min(max(XSH, -1.0), 1.0)
    YSH = min(max(YSH, -1.0), 1.0)

    return XSH, YSH, PEAKV


"""Multiline Comment15"""


def align2d_direct3(
    input_images, refim, xrng=1, yrng=1, psimax=180, psistep=1, ou=-1, CTF=None
):

    nx = input_images[0].get_xsize()
    if ou < 0:
        ou = old_div(nx, 2) - 1
    mask = sp_utilities.model_circle(ou, nx, nx)
    nk = int(old_div(psimax, psistep))
    nm = 2 * nk + 1
    nc = nk + 1
    refs = [None] * nm * 2
    for i in range(nm):
        temp = sp_fundamentals.rot_shift2D(refim, (i - nc) * psistep) * mask
        refs[2 * i] = [
            sp_fundamentals.fft(temp),
            sp_fundamentals.fft(sp_fundamentals.mirror(temp)),
        ]
        temp = sp_fundamentals.rot_shift2D(refim, (i - nc) * psistep + 180.0) * mask
        refs[2 * i + 1] = [
            sp_fundamentals.fft(temp),
            sp_fundamentals.fft(sp_fundamentals.mirror(temp)),
        ]
    del temp

    results = []
    mir = 0
    for image in input_images:
        if CTF:
            ims = sp_filter.filt_ctf(sp_fundamentals.fft(image), image.get_attr("ctf"))
        else:
            ims = sp_fundamentals.fft(image)
        ama = -1.0e23
        bang = 0.0
        bsx = 0.0
        bsy = 0.0
        for i in range(nm * 2):
            for mirror_flag in [0, 1]:
                c = sp_fundamentals.ccf(ims, refs[i][mirror_flag])
                w = EMAN2_cppwrap.Util.window(c, 2 * xrng + 1, 2 * yrng + 1)
                pp = sp_utilities.peak_search(w)[0]
                px = int(pp[4])
                py = int(pp[5])
                if pp[0] == 1.0 and px == 0 and py == 0:
                    pass  # XSH, YSH, PEAKV = 0.,0.,0.
                else:
                    ww = sp_utilities.model_blank(3, 3)
                    ux = int(pp[1])
                    uy = int(pp[2])
                    for k in range(3):
                        for l in range(3):
                            ww[k, l] = w[k + ux - 1, l + uy - 1]
                    XSH, YSH, PEAKV = parabl(ww)
                    # print i,pp[-1],XSH, YSH,px+XSH, py+YSH, PEAKV
                    if PEAKV > ama:
                        ama = PEAKV
                        bsx = px + round(XSH, 2)
                        bsy = py + round(YSH, 2)
                        bang = i
                        mir = mirror_flag
        # returned parameters have to be inverted
        bang = (old_div(bang, 2) - nc) * psistep + 180.0 * (bang % 2)
        bang, bsx, bsy, _ = sp_utilities.inverse_transform2(
            bang, (1 - 2 * mir) * bsx, bsy, mir
        )
        results.append([bang, bsx, bsy, mir, ama])
    return results


"""Multiline Comment49"""


def shc(
    data,
    refrings,
    list_of_reference_angles,
    numr,
    xrng,
    yrng,
    step,
    an=-1.0,
    sym="c1",
    finfo=None,
):

    number_of_checked_refs = 0

    mode = "F"
    nx = data.get_xsize()
    ny = data.get_ysize()
    #  center is in SPIDER convention
    cnx = old_div(nx, 2) + 1
    cny = old_div(ny, 2) + 1

    if an >= 0.0:
        ant = numpy.cos(numpy.radians(an))
    else:
        ant = -1.0
    # phi, theta, psi, sxo, syo = get_params_proj(data)
    t1 = data.get_attr("xform.projection")
    dp = t1.get_params("spider")
    ou = numr[-3]
    sxi = round(-dp["tx"], 2)
    syi = round(-dp["ty"], 2)
    txrng = search_range(nx, ou, sxi, xrng)
    tyrng = search_range(ny, ou, syi, yrng)

    if finfo:
        finfo.write(
            "Old parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n"
            % (dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"])
        )
        finfo.flush()
        z1, z2, z3, z4, z5 = sp_utilities.get_params_proj(data, "xform.anchor")
        finfo.write(
            "Anc parameters: %9.4f %9.4f %9.4f %9.4f %9.4f\n" % (z1, z2, z3, -z4, -z5)
        )
        finfo.flush()

    previousmax = data.get_attr("previousmax")
    [ang, sxs, sys, mirror, iref, peak, checked_refs] = EMAN2_cppwrap.Util.shc(
        data,
        refrings,
        list_of_reference_angles,
        txrng,
        tyrng,
        step,
        ant,
        mode,
        numr,
        cnx - sxi,
        cny - syi,
        sym,
    )
    iref = int(iref)
    number_of_checked_refs += int(checked_refs)
    if peak <= previousmax:
        return -1.0e23, 0.0, number_of_checked_refs, -1
        """Multiline Comment50"""
    else:
        # The ormqip returns parameters such that the transformation is applied first, the mirror operation second.
        # What that means is that one has to change the the Eulerian angles so they point into mirrored direction: phi+180, 180-theta, 180-psi
        if mirror:
            phi = (refrings[iref].get_attr("phi") + 540.0) % 360.0
            theta = 180.0 - refrings[iref].get_attr("theta")
            psi = (540.0 - refrings[iref].get_attr("psi") - ang) % 360.0
        else:
            phi = refrings[iref].get_attr("phi")
            theta = refrings[iref].get_attr("theta")
            psi = (360.0 + refrings[iref].get_attr("psi") - ang) % 360.0
        s2x = sxs + sxi
        s2y = sys + syi

        # set_params_proj(data, [phi, theta, psi, s2x, s2y])
        t2 = EMAN2_cppwrap.Transform(
            {"type": "spider", "phi": phi, "theta": theta, "psi": psi}
        )
        t2.set_trans(EMAN2_cppwrap.Vec2f(-s2x, -s2y))
        data.set_attr("xform.projection", t2)
        data.set_attr("previousmax", peak)
        #  Find the pixel error that is minimum over symmetry transformations
        if sym == "nomirror" or sym == "c1":
            pixel_error = sp_pixel_error.max_3D_pixel_error(t1, t2, numr[-3])
        else:
            ts = t2.get_sym_proj(sym)
            # only do it if it is not c1
            pixel_error = +1.0e23
            for ut in ts:
                # we do not care which position minimizes the error
                pixel_error = min(
                    sp_pixel_error.max_3D_pixel_error(t1, ut, numr[-3]), pixel_error
                )
        if finfo:
            finfo.write(
                "New parameters: %9.4f %9.4f %9.4f %9.4f %9.4f %10.5f  %11.3e\n\n"
                % (phi, theta, psi, s2x, s2y, peak, pixel_error)
            )
            finfo.flush()
        return peak, pixel_error, number_of_checked_refs, iref


# parameters: list of (all) projections | reference volume is optional, if provided might be shrank| ...
#  This functions centers projections using an self-correlation-based exhaustive search
#  It only returns shifts
#  Data is assumed to be shrunk and CTF-applied
#  The input volume is assumed to be shrunk but not filtered, if not provided, it will be reconstructed and shrunk
#  We apply ali3d_options.fl


def search_range(n, radius, shift, range, location=""):
    """
		Find permissible ranges for translational searches by resampling into polar coordinates
		n - image size; radius - particle radius, the circle has to fit into the square image;
		shift - current particle shift; range - desired maximum range search
		Output: a list of two elements:
		  left range (positive)
		  right range
		NOTE - ranges are with respect to the point n//2+1-shift within image (in 3D)
	"""
    cn = old_div(n, 2) + 1
    ql = cn + shift - radius - 2  # lower end is positive
    qe = n - cn - shift - radius  # upper end
    if ql < 0 or qe < 0:
        sp_global_def.ERROR(
            "Shift of particle too large, results may be incorrect:  %4d   %3d   %f  %f  %f  %f  %f"
            % (n, cn, radius, shift, range, ql, qe),
            "search_range  " + location,
            0,
        )
        ql = max(ql, 0)
        qe = max(qe, 0)
    # ???for mysterious reasons it has to be this way as C code changes the order of searches.
    return [min(qe, range), min(ql, range)]


def alivol_mask_getref(v, mask):
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_utilities import set_params3D
    v50S_ref = v.copy()
    v50S_ref *= mask
    cnt = v50S_ref.phase_cog()
    sp_utilities.set_params3D(
        v50S_ref, (0.0, 0.0, 0.0, -cnt[0], -cnt[1], -cnt[2], 0, 1.0)
    )
    return v50S_ref


def alivol_mask(v, vref, mask):
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_utilities    import set_params3D, get_params3D,compose_transform3
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_applications import ali_vol_shift, ali_vol_rotate
    v50S_i = v.copy()
    v50S_i *= mask
    cnt = v50S_i.phase_cog()
    sp_utilities.set_params3D(
        v50S_i, (0.0, 0.0, 0.0, -cnt[0], -cnt[1], -cnt[2], 0, 1.0)
    )

    v50S_i = sp_applications.ali_vol_shift(v50S_i, vref, 1.0)
    v50S_i = sp_applications.ali_vol_rotate(v50S_i, vref, 5.0)
    v50S_i = sp_applications.ali_vol_shift(v50S_i, vref, 0.5)
    v50S_i = sp_applications.ali_vol_rotate(v50S_i, vref, 1.0)
    phi, tht, psi, s3x, s3y, s3z, mirror, scale = sp_utilities.get_params3D(v50S_i)
    dun, dum, dum, cnx, cny, cnz, mirror, scale = sp_utilities.get_params3D(vref)
    phi, tht, psi, s3x, s3y, s3z, scale = sp_utilities.compose_transform3(
        phi, tht, psi, s3x, s3y, s3z, 1.0, 0.0, 0.0, 0.0, -cnx, -cny, -cnz, 1.0
    )
    return phi, tht, psi, s3x, s3y, s3z


def ali_nvol(v, mask):
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_alignment    import alivol_mask_getref, alivol_mask
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_statistics   import ave_var
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_utilities    import set_params3D, get_params3D ,compose_transform3

    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_fundamentals import rot_shift3D
    ocrit = 1.0e20
    gogo = True
    niter = 0
    for l in range(len(v)):
        sp_utilities.set_params3D(v[l], (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 1.0))
    while gogo:
        ave, var = sp_statistics.ave_var(v)
        p = EMAN2_cppwrap.Util.infomask(var, mask, True)
        crit = p[1]
        if (
            old_div(old_div((crit - ocrit), (crit + ocrit)), 2.0) > -1.0e-2
            or niter > 10
        ):
            gogo = False
        niter += 1
        ocrit = crit
        ref = alivol_mask_getref(ave, mask)
        for l in range(len(v)):
            ophi, otht, opsi, os3x, os3y, os3z, dum, dum = sp_utilities.get_params3D(
                v[l]
            )
            vor = sp_fundamentals.rot_shift3D(v[l], ophi, otht, opsi, os3x, os3y, os3z)
            phi, tht, psi, s3x, s3y, s3z = alivol_mask(vor, ref, mask)
            phi, tht, psi, s3x, s3y, s3z, scale = sp_utilities.compose_transform3(
                phi,
                tht,
                psi,
                s3x,
                s3y,
                s3z,
                1.0,
                ophi,
                otht,
                opsi,
                os3x,
                os3y,
                os3z,
                1.0,
            )
            sp_utilities.set_params3D(v[l], (phi, tht, psi, s3x, s3y, s3z, 0, 1.0))
            # print "final align3d params: %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f" % (phi,tht,psi,s3x,s3y,s3z)
    for l in range(len(v)):
        ophi, otht, opsi, os3x, os3y, os3z, dum, dum = sp_utilities.get_params3D(v[l])
        sp_global_def.sxprint(l, ophi, otht, opsi, os3x, os3y, os3z)
        v[l] = sp_fundamentals.rot_shift3D(v[l], ophi, otht, opsi, os3x, os3y, os3z)
        v[l].del_attr("xform.align3d")
    return v


def ali_vol_func_rotate(params, data):
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_utilities    import compose_transform3
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_fundamentals import rot_shift3D
    cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale = sp_utilities.compose_transform3(
        data[3][0],
        data[3][1],
        data[3][2],
        data[3][3],
        data[3][4],
        data[3][5],
        data[3][7],
        params[0],
        params[1],
        params[2],
        0.0,
        0.0,
        0.0,
        1.0,
    )
    x = sp_fundamentals.rot_shift3D(
        data[0], cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale
    )
    res = -x.cmp(data[4], data[1], {"mask": data[2]})
    # print  " %9.3f %9.3f %9.3f  %12.3e" %(params[0], params[1], params[2], -res)
    return res


def ali_vol_func_shift(params, data):
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_utilities    import compose_transform3
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_fundamentals import rot_shift3D
    cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale = sp_utilities.compose_transform3(
        data[3][0],
        data[3][1],
        data[3][2],
        data[3][3],
        data[3][4],
        data[3][5],
        data[3][7],
        0.0,
        0.0,
        0.0,
        params[0],
        params[1],
        params[2],
        1.0,
    )
    x = sp_fundamentals.rot_shift3D(
        data[0], cphi, ctheta, cpsi, cs2x, cs2y, cs2z, cscale
    )
    res = -x.cmp(data[4], data[1], {"mask": data[2]})
    # print  " %9.3f %9.3f %9.3f  %12.3e" %(params[0], params[1], params[2], -res)
    return res


"""Multiline Comment53"""


#  06-12-14 code lifted
"""Multiline Comment54"""


from builtins import range
