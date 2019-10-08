from __future__ import print_function
from past.utils import old_div
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
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
import sparx_alignment
import sparx_applications
import copy
import sparx_filter
import sparx_fundamentals
import sparx_global_def
import sparx_morphology
import mpi
import numpy
import numpy as np
import random
import time
import sparx_utilities

pass  # IMPORTIMPORTIMPORT import EMAN2
pass  # IMPORTIMPORTIMPORT import EMAN2_cppwrap
pass  # IMPORTIMPORTIMPORT import alignment
pass  # IMPORTIMPORTIMPORT import applications
pass  # IMPORTIMPORTIMPORT import copy
pass  # IMPORTIMPORTIMPORT import filter
pass  # IMPORTIMPORTIMPORT import fundamentals
pass  # IMPORTIMPORTIMPORT import global_def
pass  # IMPORTIMPORTIMPORT import math
pass  # IMPORTIMPORTIMPORT import morphology
pass  # IMPORTIMPORTIMPORT import mpi
pass  # IMPORTIMPORTIMPORT import numpy
pass  # IMPORTIMPORTIMPORT import numpy as np
pass  # IMPORTIMPORTIMPORT import projection
pass  # IMPORTIMPORTIMPORT import random
pass  # IMPORTIMPORTIMPORT import sys
pass  # IMPORTIMPORTIMPORT import time
pass  # IMPORTIMPORTIMPORT import utilities
from builtins import range

pass  # IMPORTIMPORTIMPORT from global_def import *


def project(volume, params, radius=-1):
    """
        Name
            project - calculate 2-D projection of a 3-D volume using trilinear interpolation
        Input
            vol: input volume, all dimensions have to be the same
            params: input parameters given as a list [phi, theta, psi, s2x, s2y], projection in calculated using the three Eulerian angles and then shifted by s2x,s2y
        radius: radius of a sphere within which the projection of the volume will be calculated
        Output
        proj: generated 2-D projection
    """
    # angles phi, theta, psi
    pass  # IMPORTIMPORTIMPORT from fundamentals import rot_shift2D
    pass  # IMPORTIMPORTIMPORT from utilities import set_params_proj
    pass  # IMPORTIMPORTIMPORT from EMAN2 import Processor

    if (radius > 0):
        myparams = {"transform": EMAN2_cppwrap.Transform(
            {"type": "spider", "phi": params[0], "theta": params[1], "psi": params[2]}), "radius": radius}
    else:
        myparams = {"transform": EMAN2_cppwrap.Transform(
            {"type": "spider", "phi": params[0], "theta": params[1], "psi": params[2]})}
    proj = volume.project("pawel", myparams)
    if (params[3] != 0. or params[4] != 0.):
        params2 = {"filter_type": EMAN2_cppwrap.Processor.fourier_filter_types.SHIFT, "x_shift": params[3],
                   "y_shift": params[4], "z_shift": 0.0}
        proj = EMAN2_cppwrap.Processor.EMFourierFilter(proj, params2)
    # proj = rot_shift2D(proj, sx = params[3], sy = params[4], interpolation_method = "linear")
    sparx_utilities.set_params_proj(proj, [params[0], params[1], params[2], -params[3], -params[4]])
    proj.set_attr_dict({'ctf_applied': 0})
    return proj


"""Multiline Comment0"""


def prgs(volft, kb, params, kbx=None, kby=None):
    """
        Name
            prg - calculate 2-D projection of a 3-D volume
        Input
            vol: input volume, the volume can be either cubic or rectangular
            kb: interpolants generated using prep_vol (tabulated Kaiser-Bessel function). If the volume is cubic, kb is the only interpolant.
                Otherwise, kb is the for caculating weigthing along z direction.
            kbx,kby: interpolants generated using prep_vol used to calculae weighting aling x and y directin. Default is none when the volume is cubic.
                 If the volume is rectangular, kbx and kby must be given.
            params: input parameters given as a list [phi, theta, psi, s2x, s2y], projection in calculated using the three Eulerian angles and then shifted by sx,sy
        Output
            proj: generated 2-D projection
    """
    #  params:  phi, theta, psi, sx, sy
    pass  # IMPORTIMPORTIMPORT from fundamentals import fft
    pass  # IMPORTIMPORTIMPORT from utilities import set_params_proj
    pass  # IMPORTIMPORTIMPORT from EMAN2 import Processor

    R = EMAN2_cppwrap.Transform({"type": "spider", "phi": params[0], "theta": params[1], "psi": params[2]})
    if kbx is None:
        temp = volft.extract_plane(R, kb)
    else:
        temp = volft.extract_plane_rect(R, kbx, kby, kb)

    temp.fft_shuffle()
    temp.center_origin_fft()

    if (params[3] != 0. or params[4] != 0.):
        filt_params = {"filter_type": EMAN2_cppwrap.Processor.fourier_filter_types.SHIFT,
                       "x_shift": params[3], "y_shift": params[4], "z_shift": 0.0}
        temp = EMAN2_cppwrap.Processor.EMFourierFilter(temp, filt_params)
    temp.do_ift_inplace()
    sparx_utilities.set_params_proj(temp, [params[0], params[1], params[2], -params[3], -params[4]])
    temp.set_attr_dict({'ctf_applied': 0, 'npad': 2})
    temp.depad()
    return temp


def prgl(volft, params, interpolation_method=0, return_real=True):
    """
        Name
            prgl - calculate 2-D projection of a 3-D volume using either NN Fourier or or trilinear Fourier
        Input
            vol: input volume, the volume has to be cubic
            params: input parameters given as a list [phi, theta, psi, s2x, s2y], projection in calculated using the three Eulerian angles and then shifted by sx,sy
            interpolation_method = 0  NN
            interpolation_method = 1  trilinear
            return_real:  True - return real; False - return FT of a projection.
        Output
            proj: generated 2-D projection
    """
    #  params:  phi, theta, psi, sx, sy
    pass  # IMPORTIMPORTIMPORT from fundamentals import fft
    pass  # IMPORTIMPORTIMPORT from utilities import set_params_proj, info
    pass  # IMPORTIMPORTIMPORT from EMAN2 import Processor
    if (interpolation_method < 0 or interpolation_method > 1):  sparx_global_def.ERROR(
        'Unsupported interpolation method', "interpolation_method", 1, 0)
    npad = volft.get_attr_default("npad", 1)
    R = EMAN2_cppwrap.Transform({"type": "spider", "phi": params[0], "theta": params[1], "psi": params[2]})
    if (npad == 1):
        temp = volft.extract_section(R, interpolation_method)
    elif (npad == 2):
        temp = volft.extract_section2(R, interpolation_method)
    temp.fft_shuffle()
    temp.center_origin_fft()

    if (params[3] != 0. or params[4] != 0.):
        filt_params = {"filter_type": EMAN2_cppwrap.Processor.fourier_filter_types.SHIFT,
                       "x_shift": params[3], "y_shift": params[4], "z_shift": 0.0}
        temp = EMAN2_cppwrap.Processor.EMFourierFilter(temp, filt_params)
    if return_real:
        temp.do_ift_inplace()
        temp.set_attr_dict({'ctf_applied': 0, 'npad': 1})
        temp.depad()
    else:
        temp.set_attr_dict({'ctf_applied': 0, 'npad': 1})
    sparx_utilities.set_params_proj(temp, [params[0], params[1], params[2], -params[3], -params[4]])
    return temp


def prgq(volft, kb, nx, delta, ref_a, sym, MPI=False):
    """
      Generate set of projections based on even angles
      The command returns list of ffts of projections
    """
    pass  # IMPORTIMPORTIMPORT from projection   import prep_vol, prgs
    pass  # IMPORTIMPORTIMPORT from applications import MPI_start_end
    pass  # IMPORTIMPORTIMPORT from utilities    import even_angles, model_blank
    pass  # IMPORTIMPORTIMPORT from fundamentals import fft
    # generate list of Eulerian angles for reference projections
    #  phi, theta, psi
    mode = "F"
    ref_angles = sparx_utilities.even_angles(delta, symmetry=sym, method=ref_a, phiEqpsi="Minus")
    cnx = old_div(nx, 2) + 1
    cny = old_div(nx, 2) + 1
    num_ref = len(ref_angles)

    if MPI:
        pass  # IMPORTIMPORTIMPORT from mpi import mpi_comm_rank, mpi_comm_size, MPI_COMM_WORLD
        myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
        ncpu = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
    else:
        ncpu = 1
        myid = 0
    pass  # IMPORTIMPORTIMPORT from applications import MPI_start_end
    ref_start, ref_end = sparx_applications.MPI_start_end(num_ref, ncpu, myid)

    prjref = []  # list of (image objects) reference projections in Fourier representation

    for i in range(num_ref):
        prjref.append(sparx_utilities.model_blank(nx, nx))  # I am not sure why is that necessary, why not put None's??

    for i in range(ref_start, ref_end):
        prjref[i] = prgs(volft, kb, [ref_angles[i][0], ref_angles[i][1], ref_angles[i][2], 0.0, 0.0])

    if MPI:
        pass  # IMPORTIMPORTIMPORT from utilities import bcast_EMData_to_all
        for i in range(num_ref):
            for j in range(ncpu):
                ref_start, ref_end = sparx_applications.MPI_start_end(num_ref, ncpu, j)
                if i >= ref_start and i < ref_end: rootid = j
            sparx_utilities.bcast_EMData_to_all(prjref[i], myid, rootid)

    for i in range(len(ref_angles)):
        prjref[i].set_attr_dict({"phi": ref_angles[i][0], "theta": ref_angles[i][1], "psi": ref_angles[i][2]})

    return prjref


def prg(volume, params):
    """Given a volume, a set of projection angles, and Kaiser-Bessel
       window parameters, use gridding to generate projection
    """
    Mx = volume.get_xsize()
    My = volume.get_ysize()
    Mz = volume.get_zsize()
    if (Mx == Mz & My == Mz):
        volft, kb = prep_vol(volume)
        return prgs(volft, kb, params)
    else:
        volft, kbx, kby, kbz = prep_vol(volume)
        return prgs(volft, kbz, params, kbx, kby)


def prep_vol(vol, npad=2, interpolation_method=-1):
    """
        Name
            prep_vol - prepare the volume for calculation of gridding projections and generate the interpolants.
        Input
            vol: input volume for which projections will be calculated using prgs (interpolation_method=-1) or prgl (interpolation_method>0)
            interpolation_method = -1  gridding
            interpolation_method =  0  NN
            interpolation_method =  1  trilinear
        Output
            volft: volume prepared for gridding projections using prgs
            kb: interpolants (tabulated Kaiser-Bessel function) when the volume is cubic.
            kbx,kby: interpolants along x, y and z direction (tabulated Kaiser-Bessel function) when the volume is rectangular
    """
    # prepare the volume
    Mx = vol.get_xsize()
    My = vol.get_ysize()
    Mz = vol.get_zsize()
    #  gridding
    if interpolation_method == -1:
        K = 6
        alpha = 1.75
        assert npad == 2
        if (Mx == Mz & My == Mz):
            M = vol.get_xsize()
            # padd two times
            N = M * npad
            # support of the window
            kb = EMAN2_cppwrap.Util.KaiserBessel(alpha, K, old_div(M, 2), old_div(K, (2. * N)), N)
            volft = vol.copy()
            volft.divkbsinh(kb)
            volft = volft.norm_pad(False, npad)
            volft.do_fft_inplace()
            volft.center_origin_fft()
            volft.fft_shuffle()
            return volft, kb
        else:
            Nx = Mx * npad
            Ny = My * npad
            Nz = Mz * npad
            # support of the window
            kbx = EMAN2_cppwrap.Util.KaiserBessel(alpha, K, old_div(Mx, 2), old_div(K, (2. * Nx)), Nx)
            kby = EMAN2_cppwrap.Util.KaiserBessel(alpha, K, old_div(My, 2), old_div(K, (2. * Ny)), Ny)
            kbz = EMAN2_cppwrap.Util.KaiserBessel(alpha, K, old_div(Mz, 2), old_div(K, (2. * Nz)), Nz)
            volft = vol.copy()
            volft.divkbsinh_rect(kbx, kby, kbz)
            volft = volft.norm_pad(False, npad)
            volft.do_fft_inplace()
            volft.center_origin_fft()
            volft.fft_shuffle()
            return volft, kbx, kby, kbz
    else:
        # NN and trilinear
        assert interpolation_method >= 0
        pass  # IMPORTIMPORTIMPORT from utilities import pad
        volft = sparx_utilities.pad(vol, Mx * npad, My * npad, My * npad, 0.0)
        volft.set_attr("npad", npad)
        volft.div_sinc(interpolation_method)
        volft = volft.norm_pad(False, 1)
        volft.do_fft_inplace()
        volft.center_origin_fft()
        volft.fft_shuffle()
        volft.set_attr("npad", npad)
        return volft


def gen_rings_ctf(prjref, nx, ctf, numr):
    """
      Convert set of ffts of projections to Fourier rings with additional multiplication by a ctf
      The command returns list of rings
    """
    pass  # IMPORTIMPORTIMPORT from math         import sin, cos, pi
    pass  # IMPORTIMPORTIMPORT from fundamentals import fft
    pass  # IMPORTIMPORTIMPORT from alignment    import ringwe
    pass  # IMPORTIMPORTIMPORT from filter       import filt_ctf
    mode = "F"
    wr_four = sparx_alignment.ringwe(numr, "F")
    cnx = old_div(nx, 2) + 1
    cny = old_div(nx, 2) + 1
    qv = old_div(numpy.pi, 180.0)

    refrings = []  # list of (image objects) reference projections in Fourier representation

    for i in range(len(prjref)):
        cimage = EMAN2_cppwrap.Util.Polar2Dm(sparx_filter.filt_ctf(prjref[i], ctf, True), cnx, cny, numr,
                                             mode)  # currently set to quadratic....
        EMAN2_cppwrap.Util.Normalize_ring(cimage, numr, 0)

        EMAN2_cppwrap.Util.Frngs(cimage, numr)
        EMAN2_cppwrap.Util.Applyws(cimage, numr, wr_four)
        refrings.append(cimage)
        phi = prjref[i].get_attr('phi')
        theta = prjref[i].get_attr('theta')
        psi = prjref[i].get_attr('psi')
        n1 = numpy.sin(theta * qv) * numpy.cos(phi * qv)
        n2 = numpy.sin(theta * qv) * numpy.sin(phi * qv)
        n3 = numpy.cos(theta * qv)
        refrings[i].set_attr_dict({"n1": n1, "n2": n2, "n3": n3, "phi": phi, "theta": theta, "psi": psi})

    return refrings

###############################################################################################
## COMMON LINES NEW VERSION ###################################################################

# plot angles, map on half-sphere
# agls: [[phi0, theta0, psi0], [phi1, theta1, psi1], ..., [phin, thetan, psin]]
