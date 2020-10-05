
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
from . import sp_global_def
from . import sp_utilities


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

    if radius > 0:
        myparams = {
            "transform": EMAN2_cppwrap.Transform(
                {
                    "type": "spider",
                    "phi": params[0],
                    "theta": params[1],
                    "psi": params[2],
                }
            ),
            "radius": radius,
        }
    else:
        myparams = {
            "transform": EMAN2_cppwrap.Transform(
                {
                    "type": "spider",
                    "phi": params[0],
                    "theta": params[1],
                    "psi": params[2],
                }
            )
        }
    proj = volume.project("pawel", myparams)
    if params[3] != 0.0 or params[4] != 0.0:
        params2 = {
            "filter_type": EMAN2_cppwrap.Processor.fourier_filter_types.SHIFT,
            "x_shift": params[3],
            "y_shift": params[4],
            "z_shift": 0.0,
        }
        proj = EMAN2_cppwrap.Processor.EMFourierFilter(proj, params2)
        # proj = rot_shift2D(proj, sx = params[3], sy = params[4], interpolation_method = "linear")
    sp_utilities.set_params_proj(
        proj, [params[0], params[1], params[2], -params[3], -params[4]]
    )
    proj.set_attr_dict({"ctf_applied": 0})
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

    R = EMAN2_cppwrap.Transform(
        {"type": "spider", "phi": params[0], "theta": params[1], "psi": params[2]}
    )
    if kbx is None:
        temp = volft.extract_plane(R, kb)
    else:
        temp = volft.extract_plane_rect(R, kbx, kby, kb)

    temp.fft_shuffle()
    temp.center_origin_fft()

    if params[3] != 0.0 or params[4] != 0.0:
        filt_params = {
            "filter_type": EMAN2_cppwrap.Processor.fourier_filter_types.SHIFT,
            "x_shift": params[3],
            "y_shift": params[4],
            "z_shift": 0.0,
        }
        temp = EMAN2_cppwrap.Processor.EMFourierFilter(temp, filt_params)
    temp.do_ift_inplace()
    sp_utilities.set_params_proj(
        temp, [params[0], params[1], params[2], -params[3], -params[4]]
    )
    temp.set_attr_dict({"ctf_applied": 0, "npad": 2})
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
    if interpolation_method < 0 or interpolation_method > 1:
        sp_global_def.ERROR(
            "Unsupported interpolation method", "interpolation_method", 1, 0
        )
    npad = volft.get_attr_default("npad", 1)
    R = EMAN2_cppwrap.Transform(
        {"type": "spider", "phi": params[0], "theta": params[1], "psi": params[2]}
    )
    if npad == 1:
        temp = volft.extract_section(R, interpolation_method)
    elif npad == 2:
        temp = volft.extract_section2(R, interpolation_method)
    temp.fft_shuffle()
    temp.center_origin_fft()

    if params[3] != 0.0 or params[4] != 0.0:
        filt_params = {
            "filter_type": EMAN2_cppwrap.Processor.fourier_filter_types.SHIFT,
            "x_shift": params[3],
            "y_shift": params[4],
            "z_shift": 0.0,
        }
        temp = EMAN2_cppwrap.Processor.EMFourierFilter(temp, filt_params)
    if return_real:
        temp.do_ift_inplace()
        temp.set_attr_dict({"ctf_applied": 0, "npad": 1})
        temp.depad()
    else:
        temp.set_attr_dict({"ctf_applied": 0, "npad": 1})
    sp_utilities.set_params_proj(
        temp, [params[0], params[1], params[2], -params[3], -params[4]]
    )
    return temp


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
        if Mx == Mz & My == Mz:
            M = vol.get_xsize()
            # padd two times
            N = M * npad
            # support of the window
            kb = EMAN2_cppwrap.Util.KaiserBessel(
                alpha, K, old_div(M, 2), old_div(K, (2.0 * N)), N
            )
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
            kbx = EMAN2_cppwrap.Util.KaiserBessel(
                alpha, K, old_div(Mx, 2), old_div(K, (2.0 * Nx)), Nx
            )
            kby = EMAN2_cppwrap.Util.KaiserBessel(
                alpha, K, old_div(My, 2), old_div(K, (2.0 * Ny)), Ny
            )
            kbz = EMAN2_cppwrap.Util.KaiserBessel(
                alpha, K, old_div(Mz, 2), old_div(K, (2.0 * Nz)), Nz
            )
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
        volft = sp_utilities.pad(vol, Mx * npad, My * npad, My * npad, 0.0)
        volft.set_attr("npad", npad)
        volft.div_sinc(interpolation_method)
        volft = volft.norm_pad(False, 1)
        volft.do_fft_inplace()
        volft.center_origin_fft()
        volft.fft_shuffle()
        volft.set_attr("npad", npad)
        return volft
