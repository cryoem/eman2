
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
import math
import mpi
from . import sp_fundamentals
from . import sp_global_def
from . import sp_morphology
from . import sp_utilities


def filt_tophatb(e, freql, freqh, pad=False):
    """
		Name
			filt_tophatb - top-hat band-pass Fourier filter (truncation of a Fourier series)
		Input
			e: input image (can be either real or Fourier)
			freql: low-end frequency of the filter pass-band
			freqh: high-end frequency of the filter pass-band
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			All frequencies are in absolute frequency units fa and their valid range is [0:0.5].
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
    params = {
        "filter_type": EMAN2_cppwrap.Processor.fourier_filter_types.TOP_HAT_BAND_PASS,
        "low_cutoff_frequency": freql,
        "high_cutoff_frequency": freqh,
        "dopad": pad,
    }
    return EMAN2_cppwrap.Processor.EMFourierFilter(e, params)


def filt_gaussl(e, sigma, pad=False):
    """
		Name
			filt_gaussl - Gaussian low-pass Fourier filter
		Input
			e: input image (can be either real or Fourier)
			sigma: standard deviation of the Gaussian function in absolute frequency units fa.
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			Note: for sigma = 0.5 (the Nyquist frequency) the value of the filter at the maximum frequency is G(fN)=1e=0.61.
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
    params = {
        "filter_type": EMAN2_cppwrap.Processor.fourier_filter_types.GAUSS_LOW_PASS,
        "cutoff_abs": sigma,
        "dopad": pad,
    }
    return EMAN2_cppwrap.Processor.EMFourierFilter(e, params)


def filt_gaussinv(e, sigma, pad=False):
    """
		Name
			filt_gaussinv - inverse Gaussian (high-pass) Fourier filter (division by a Gaussian function in Fourier space)
		Input
			e: input image (can be either real or Fourier)
			sigma: standard deviation of the Gaussian function in absolute frequency units fa.
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			Note: for sigma = 0.5 (the Nyquist frequency) the value of the filter at the maximum frequency is G(fN)=e=1.65.
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
    params = {
        "filter_type": EMAN2_cppwrap.Processor.fourier_filter_types.GAUSS_INVERSE,
        "cutoff_abs": sigma,
        "dopad": pad,
    }
    return EMAN2_cppwrap.Processor.EMFourierFilter(e, params)


def filt_gaussh(e, sigma, pad=False):
    """
		Name
			filt_gaussh - Gaussian high-pass Fourier filter
		Input
			e: input image (can be either real or Fourier)
			sigma: standard deviation of the Gaussian function in absolute frequency units fa.
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			Note: for sigma = 0.5 (the Nyquist frequency) the value of the filter at the maximum frequency is 1-G(0)=1-1e=0.39
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
    params = {
        "filter_type": EMAN2_cppwrap.Processor.fourier_filter_types.GAUSS_HIGH_PASS,
        "cutoff_abs": sigma,
        "dopad": pad,
    }
    return EMAN2_cppwrap.Processor.EMFourierFilter(e, params)


def filt_btwl(e, freql, freqh, pad=False):
    """
		Name
			filt_btwl - Butterworth low-pass Fourier filter
		Input
			e: input image (can be either real or Fourier)
			freql: low - pass-band frequency
			freqh: high - stop-band frequency
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			All frequencies are in absolute frequency units fa and their valid range is [0:0.5].
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
    params = {
        "filter_type": EMAN2_cppwrap.Processor.fourier_filter_types.BUTTERWORTH_LOW_PASS,
        "low_cutoff_frequency": freql,
        "high_cutoff_frequency": freqh,
        "dopad": pad,
    }
    return EMAN2_cppwrap.Processor.EMFourierFilter(e, params)


def filt_tanl(e, freq, fall_off, pad=False):
    """
		Name
			filt_tanl - hyperbolic tangent low-pass Fourier filter
		Input
			e: input image (can be either real or Fourier)
			freq: stop-band frequency
			fall_off: fall off of the filter
			pad: logical flag specifying whether before filtering the image should be padded with zeroes in real space to twice the size (this helps avoiding aliasing artifacts). (Default pad = False).
			All frequencies are in absolute frequency units fa and their valid range is [0:0.5].
		Output
			filtered image. Output image is real when input image is real or Fourier when input image is Fourier
	"""
    params = {
        "filter_type": EMAN2_cppwrap.Processor.fourier_filter_types.TANH_LOW_PASS,
        "cutoff_abs": freq,
        "fall_off": fall_off,
        "dopad": pad,
    }
    return EMAN2_cppwrap.Processor.EMFourierFilter(e, params)


def filt_table(e, table):
    """
		Name
			filt_table - filter image using a user-provided (as a list) of Fourier filter values
		Input
			e: input image (real or Fourier)
			table: user-provided 1-D table of filter values.
		Output
			image Fourier-filtered using coefficients in table (real or Fourier, according the input image)
		Options
			fast: use the fast method; may combust certain computers.
			huge: gobble memory; there is plenty of it, anyway.
	"""
    params = {
        "filter_type": EMAN2_cppwrap.Processor.fourier_filter_types.RADIAL_TABLE,
        "table": table,
    }
    return EMAN2_cppwrap.Processor.EMFourierFilter(e, params)


def filt_ctf(img, ctf, dopad=True, sign=1, binary=0):
    """
		Name
			filt_ctf - apply Contrast Transfer Function (CTF) to an image in Fourier space
		Input
			image: input image, it can be real or complex
			ctf: an CTF object, please see CTF_info for description.
			pad: apply padding with zeroes in real space to twice the size before CTF application (Default is True, if set to False, no padding, if input image is Fourier, the flag has no effect).
			sign: sign of the CTF. If cryo data had inverted contrast, i.e., images were multiplied by -1 and particle projections appear bright on dark background, it has to be set to -1). (Default is 1).
			binary: phase flipping if set to 1 (default is 0).
		Output
			image multiplied in Fourier space by the CTF, the output image has the same format as the input image.
	"""
    assert img.get_ysize() > 1
    dict = ctf.to_dict()
    dz = dict["defocus"]
    cs = dict["cs"]
    voltage = dict["voltage"]
    pixel_size = dict["apix"]
    b_factor = dict["bfactor"]
    ampcont = dict["ampcont"]
    dza = dict["dfdiff"]
    azz = dict["dfang"]

    if dopad and not img.is_complex():
        ip = 1
    else:
        ip = 0

    params = {
        "filter_type": EMAN2_cppwrap.Processor.fourier_filter_types.CTF_,
        "defocus": dz,
        "Cs": cs,
        "voltage": voltage,
        "Pixel_size": pixel_size,
        "B_factor": b_factor,
        "amp_contrast": ampcont,
        "dopad": ip,
        "binary": binary,
        "sign": sign,
        "dza": dza,
        "azz": azz,
    }
    tmp = EMAN2_cppwrap.Processor.EMFourierFilter(img, params)
    tmp.set_attr_dict({"ctf": ctf})
    return tmp


def fit_tanh(dres, low=0.1):
    """
		dres - list produced by the fsc funcion
		dres[0] - absolute frequencies
		dres[1] - fsc, because it was calculated from the dataset split into halves, convert it to full using rn = 2r/(1+r)
		dres[2] - number of points use to calculate fsc coeff
		low cutoff of the fsc curve
		return parameters of the tanh filter: freq - cutoff frequency at which filter value is 0.5, and fall_off, the 'width' of the filter
	"""

    def fit_tanh_func(args, data):
        v = 0.0

        if data[1][0] < 0.0:
            data[1][0] *= -1.0

        for i in range(len(data[0])):
            fsc = old_div(2 * data[1][i], (1.0 + data[1][i]))
            if args[0] == 0 or args[1] == 0:
                qt = 0
            else:
                qt = fsc - 0.5 * (
                    math.tanh(
                        old_div(
                            old_div(
                                old_div(math.pi * (data[0][i] + args[0]), 2.0), args[1]
                            ),
                            args[0],
                        )
                    )
                    - math.tanh(
                        old_div(
                            old_div(
                                old_div(math.pi * (data[0][i] - args[0]), 2.0), args[1]
                            ),
                            args[0],
                        )
                    )
                )
            v -= qt * qt
        # print args,v
        return v

    setzero = False
    for i in range(1, len(dres[0])):
        if not setzero:
            if old_div(2 * dres[1][i], (1.0 + dres[1][i])) < low:
                setzero = True
        if setzero:
            dres[1][i] = 0.0

    freq = -1.0
    for i in range(1, len(dres[0]) - 1):
        if (old_div(2 * dres[1][i], (1.0 + dres[1][i]))) < 0.5:
            freq = dres[0][i - 1]
            break
    if freq < 0.0:
        # the curve never falls below 0.5, most likely something's wrong; however, return reasonable values
        freq = 0.4
        fall_off = 0.2
        return freq, fall_off

    args = [freq, 0.1]
    scale = [0.05, 0.05]
    result = sp_utilities.amoeba(args, scale, fit_tanh_func, data=dres)

    """Multiline Comment0"""
    return result[0][0], result[0][1]


def filt_vols(vols, fscs, mask3D):

    flmin = 1.0
    flmax = -1.0
    nvol = len(vols)
    for i in range(nvol):
        fl, aa = fit_tanh(fscs[i])
        if fl < flmin:
            flmin = fl
            aamin = aa
        if fl > flmax:
            flmax = fl
            idmax = i
    sp_global_def.sxprint(" Filter tanl, parameters: ", flmin - 0.05, "  ", aamin)
    volmax = vols[idmax]
    volmax = filt_tanl(volmax, flmin - 0.05, aamin)
    pmax = sp_fundamentals.rops_table(volmax)

    for i in range(nvol):
        ptab = sp_fundamentals.rops_table(vols[i])
        for j in range(len(ptab)):
            ptab[j] = math.sqrt(old_div(pmax[j], ptab[j]))

        vols[i] = filt_table(vols[i], ptab)
        # stat = Util.infomask( vols[i], mask3D, False )
        # volf -= stat[0]
        EMAN2_cppwrap.Util.mul_img(vols[i], mask3D)
        # volf = threshold( volf )

    return vols


def filterlocal(ui, vi, m, falloff, myid, main_node, number_of_proc):

    if myid == main_node:

        nx = vi.get_xsize()
        ny = vi.get_ysize()
        nz = vi.get_zsize()
        #  Round all resolution numbers to two digits
        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    ui.set_value_at_fast(x, y, z, round(ui.get_value_at(x, y, z), 2))
        dis = [nx, ny, nz]
    else:
        falloff = 0.0
        radius = 0
        dis = [0, 0, 0]
    falloff = sp_utilities.bcast_number_to_all(falloff, main_node)
    dis = sp_utilities.bcast_list_to_all(dis, myid, source_node=main_node)

    if myid != main_node:
        nx = int(dis[0])
        ny = int(dis[1])
        nz = int(dis[2])

        vi = sp_utilities.model_blank(nx, ny, nz)
        ui = sp_utilities.model_blank(nx, ny, nz)

    sp_utilities.bcast_EMData_to_all(vi, myid, main_node)
    sp_utilities.bcast_EMData_to_all(ui, myid, main_node)

    sp_fundamentals.fftip(vi)  #  volume to be filtered

    st = EMAN2_cppwrap.Util.infomask(ui, m, True)

    filteredvol = sp_utilities.model_blank(nx, ny, nz)
    cutoff = max(st[2] - 0.01, 0.0)
    while cutoff < st[3]:
        cutoff = round(cutoff + 0.01, 2)
        # if(myid == main_node):  print  cutoff,st
        pt = EMAN2_cppwrap.Util.infomask(
            sp_morphology.threshold_outside(ui, cutoff - 0.00501, cutoff + 0.005),
            m,
            True,
        )  # Ideally, one would want to check only slices in question...
        if pt[0] != 0.0:
            # print cutoff,pt[0]
            vovo = sp_fundamentals.fft(filt_tanl(vi, cutoff, falloff))
            for z in range(myid, nz, number_of_proc):
                for x in range(nx):
                    for y in range(ny):
                        if m.get_value_at(x, y, z) > 0.5:
                            if round(ui.get_value_at(x, y, z), 2) == cutoff:
                                filteredvol.set_value_at_fast(
                                    x, y, z, vovo.get_value_at(x, y, z)
                                )

    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    sp_utilities.reduce_EMData_to_root(filteredvol, myid, main_node, mpi.MPI_COMM_WORLD)
    return filteredvol


from builtins import range
