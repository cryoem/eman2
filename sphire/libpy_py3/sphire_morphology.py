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
import EMAN2db
import sparx_alignment
import sparx_applications
import sparx_filter
import sparx_fundamentals
import glob
import sparx_global_def
import inspect
import math
import mpi
import numpy
import numpy as np
import operator
import os
import sparx_pixel_error
import sparx_projection
import random
import sparx_statistics
import sys
import sparx_utilities

pass  # IMPORTIMPORTIMPORT import EMAN2
pass  # IMPORTIMPORTIMPORT import EMAN2_cppwrap
pass  # IMPORTIMPORTIMPORT import EMAN2db
pass  # IMPORTIMPORTIMPORT import alignment
pass  # IMPORTIMPORTIMPORT import applications
pass  # IMPORTIMPORTIMPORT import filter
pass  # IMPORTIMPORTIMPORT import fundamentals
pass  # IMPORTIMPORTIMPORT import glob
pass  # IMPORTIMPORTIMPORT import global_def
pass  # IMPORTIMPORTIMPORT import inspect
pass  # IMPORTIMPORTIMPORT import math
pass  # IMPORTIMPORTIMPORT import morphology
pass  # IMPORTIMPORTIMPORT import mpi
pass  # IMPORTIMPORTIMPORT import numpy
pass  # IMPORTIMPORTIMPORT import numpy as np
pass  # IMPORTIMPORTIMPORT import operator
pass  # IMPORTIMPORTIMPORT import os
pass  # IMPORTIMPORTIMPORT import pixel_error
pass  # IMPORTIMPORTIMPORT import projection
pass  # IMPORTIMPORTIMPORT import random
pass  # IMPORTIMPORTIMPORT import scipy.optimize
pass  # IMPORTIMPORTIMPORT import statistics
pass  # IMPORTIMPORTIMPORT import sys
pass  # IMPORTIMPORTIMPORT import time
pass  # IMPORTIMPORTIMPORT import types
pass  # IMPORTIMPORTIMPORT import utilities
from builtins import range

pass  # IMPORTIMPORTIMPORT from global_def import *


def binarize(img, minval=0.0):
    """
		Name
			binarize - create a binary image from an input image
		Input
			img: input image
			minval: value below which image pixels will be set to zero.
		Output
			binarized image.
	"""
    return img.process("threshold.binary", {"value": minval})


def collapse(img, minval=-1.0, maxval=1.0):
    """
		Name
			collapse - binarize image by setting to zero values within predefined range, and to one outside of this range
		Input
			input: input image
			minval: minimum bracket value (default is -1.0).
			maxval: maximum bracket value (default is 1.0).
		Output
			collapsed image.
	"""
    # for values between minval and maxval set to one, to zero outside
    return img.process("threshold.binaryrange", {"low": minval, "high": maxval})


def dilation(f, mask=None, morphtype="BINARY"):
    """
		Name
			dilation - Calculate the dilated image.
		Input
			The first input image
			mask: The second input image used as the mask.
			The size of the mask has to be odd so that the center of mask can be well defined.
			The size of the mask should be smaller than the size of the first input image.
			morph_type: Type of the dilation
			BINARY is for binary dilation;
			GRAYLEVEL is for graylevel dilation.
		Output
			dilated image
	"""
    pass  # IMPORTIMPORTIMPORT from EMAN2 import morph_type, filt_dilation_

    if not mask:
        pass  # IMPORTIMPORTIMPORT from utilities import model_circle
        nx = f.get_xsize()
        ny = f.get_ysize()
        nz = f.get_zsize()
        if (nz == 1):
            mask = sparx_utilities.model_circle(2, 5, 5)
        elif (nz > 1):
            mask = sparx_utilities.model_circle(2, 5, 5)
        else:
            sparx_global_def.ERROR("Command does not work for 1D images", "dilation", 1)

    if morphtype == "BINARY":
        return EMAN2_cppwrap.filt_dilation_(f, mask, EMAN2_cppwrap.morph_type.BINARY)
    elif morphtype == "GRAYLEVEL":
        return EMAN2_cppwrap.filt_dilation_(f, mask, EMAN2_cppwrap.morph_type.GRAYLEVEL)
    else:
        sparx_global_def.ERROR("Unknown dilation type", "dilation", 1)


def erosion(f, mask=None, morphtype="BINARY"):
    """
		Name
			erosion - Calculate the eroded image.
		Input
			The first input image
			mask: The second input image used as the mask.
			The size of the mask has to be odd so that the center of mask can be well defined.
			The size of the mask should be smaller than the size of the first input image.
			morph_type: Type of the erosion
			BINARY is for binary erosion (DEFAULT);
			GRAYLEVEL is for graylevel erosion.
		Output
			eroded image
	"""
    pass  # IMPORTIMPORTIMPORT from EMAN2 import morph_type, filt_erosion_

    if not mask:
        pass  # IMPORTIMPORTIMPORT from utilities import model_blank
        nx = f.get_xsize()
        ny = f.get_ysize()
        nz = f.get_zsize()
        pass  # IMPORTIMPORTIMPORT from utilities import model_circle
        nx = f.get_xsize()
        ny = f.get_ysize()
        nz = f.get_zsize()
        if (nz == 1):
            mask = sparx_utilities.model_circle(2, 5, 5)
        elif (nz > 1):
            mask = sparx_utilities.model_circle(2, 5, 5)
        else:
            sparx_global_def.ERROR("Command does not work for 1D images", "erosion", 1)

    if morphtype == "BINARY":
        return EMAN2_cppwrap.filt_erosion_(f, mask, EMAN2_cppwrap.morph_type.BINARY)
    elif morphtype == "GRAYLEVEL":
        return EMAN2_cppwrap.filt_erosion_(f, mask, EMAN2_cppwrap.morph_type.GRAYLEVEL)
    else:
        sparx_global_def.ERROR("Unknown erosion type", "erosion", 1)


def power(img, x=3.0):
    """
		Name
			power - generate image whose pixels are generated of raising to a given power pixels of the input image
		Input
			image: input real image
		Output
			the output image whose pixels are given by o=ir
			x: power
	"""
    return img.process("math.pow", {"pow": x})


def square_root(img):
    """
		Name
			square_root - create image whose pixels will be square roots of pixels of the input image
		Input
			input image
		Output
			output image.
	"""
    [a, b, c, d] = EMAN2_cppwrap.Util.infomask(img, None, False)
    if (c < 0.0):  sparx_global_def.ERROR("Cannot calculate square root of negative pixels", "square_root", 1)
    return img.process("math.sqrt")


def square(img):
    """
		Name
			square - create image whose pixels will be squared values of pixels of the input image
		Input
			input image
		Output
			output image.
	"""
    return img.process("math.squared")


def threshold(img, minval=0.0):
    """
		Name
			threshold - replace values below given threshold by zero
		Input
			img: input image
			minval: value below which image pixels will be set to zero.
		Output
			thresholded image.
	"""
    return img.process("threshold.belowtozero", {"minval": minval})


def threshold_outside(img, minval, maxval):
    """
		Name
			threshold_outside - replace values outside given thresholds by respective threshold values
		Input
			img: input image
			minval: value below which image pixels will be set to this value.
			maxval: value above which image pixels will be set to this value.
	"""
    return img.process("threshold.clampminmax", {"minval": minval, "maxval": maxval})


def notzero(img):
    """
		Name
			notzero - replace values that are not zero by 1.0
		Input
			img: input image
		Output
			binary image.
	"""
    return img.process("threshold.notzero")


def rotavg_ctf(img, defocus, Cs, voltage, Pixel_size, amp=0.0, ang=0.0):
    """1D rotational average of a 2D power spectrum (img)
	   based on estimated CTF parameters, including astigmatism amplitude and angle
	"""
    pass  # IMPORTIMPORTIMPORT from math import sqrt,atan2,pi,sin,radians
    defc = defocus * 10000
    astigmamp = amp * 10000
    lam = old_div(12.398, numpy.sqrt(voltage * (1022.0 + voltage)))
    angrad = numpy.radians(ang)
    nx = img.get_xsize()
    lr = [0.0] * 2 * ( old_div(nx, 2)  + 1)
    cnt = [0.0] * 2 * ( old_div(nx, 2)  + 1)
    nc = old_div(nx, 2)
    nr2 = nc * nc + 1
    if (Cs == 0.0):
        for ix in range(nx):
            x = ix - nc
            for iy in range(nx):
                y = iy - nc
                r2 = x * x + y * y
                if (r2 < nr2):
                    dfa = defc - old_div(astigmamp, 2) * numpy.sin(2 * (math.atan2(y, x) + angrad))
                    try:
                        u = numpy.sqrt( old_div(dfa, defc)) * numpy.sqrt(r2)
                        iu = int(u)
                        du = u - iu
                        lr[iu] += (1.0 - du) * img.get_value_at(ix, iy)
                        lr[iu + 1] += du * img.get_value_at(ix, iy)
                        cnt[iu] += 1.0 - du
                        cnt[iu + 1] += du
                    except:
                        pass
    else:
        Cst = Cs * 1.e7
        for ix in range(nx):
            x = ix - nc
            for iy in range(nx):
                y = iy - nc
                r2 = x * x + y * y
                if (r2 < nr2):
                    s = old_div(numpy.sqrt(r2), (nc * 2 * Pixel_size))
                    dfa = defc - old_div(astigmamp, 2) * numpy.sin(2 * (math.atan2(y, x) + angrad))
                    # u = sqrt(r2)*sqrt(1.0 -  astigmamp/2./defc*sin(2*(atan2(y,x) - angrad)))
                    # print ix,iy,sqrt(r2),defc,dfa,lam,s,u
                    # print  ix,iy,sqrt(r2),defc**2 + Cst**2*lam**4*s**4 - 2*Cst*lam**2*s**2*dfa
                    # print  defc
                    # print  defc - sqrt( defc**2 + Cst**2*lam**4*s**4 - 2*Cst*lam**2*s**2*dfa)
                    try:
                        u = old_div(numpy.sqrt(Cst * (defc - numpy.sqrt(
                            defc ** 2 + Cst ** 2 * lam ** 4 * s ** 4 - 2 * Cst * lam ** 2 * s ** 2 * dfa))), \
                            (Cst * lam)) * nc * 2 * Pixel_size
                        iu = int(u)
                        du = u - iu
                        lr[iu] += (1.0 - du) * img.get_value_at(ix, iy)
                        lr[iu + 1] += du * img.get_value_at(ix, iy)
                        cnt[iu] += 1.0 - du
                        cnt[iu + 1] += du
                    except:
                        pass
    for ix in range(nc):  lr[ix] = old_div(lr[ix], max(cnt[ix], 1.0))
    return lr[:nc]


def ctf_1d(nx, ctf, sign=1, doabs=False):
    """
		Generate a list of 1D CTF values
		Input
			nx: image size to which CTF will be applied.
			ctf: CTF object created using generate_ctf
			sign: sign of the CTF.
		Output
			a list of CTF values.
	"""
    dict = ctf.to_dict()
    dz = dict["defocus"]
    cs = dict["cs"]
    voltage = dict["voltage"]
    pixel_size = dict["apix"]
    bfactor = dict["bfactor"]
    ampcont = dict["ampcont"]

    ctf_1 = []
    scl =  old_div(old_div(1., pixel_size) , nx)
    length = int(1.41 * float( old_div(nx, 2))) + 1
    ctf_1 = [0.0] * length
    if doabs:
        for i in range(length): ctf_1[i] = abs(EMAN2_cppwrap.Util.tf(dz, i * scl, voltage, cs, ampcont, bfactor, sign))
    else:
        for i in range(length): ctf_1[i] = EMAN2_cppwrap.Util.tf(dz, i * scl, voltage, cs, ampcont, bfactor, sign)
    return ctf_1


def ctf_2(nx, ctf):
    """
		Generate a list of 1D CTF^2 values
		Input
			nx: image size to which CTF will be applied.
			ctf: ctf object created using generate_ctf
		Output
			a list of CTF2 values.
	"""
    dict = ctf.to_dict()
    dz = dict["defocus"]
    cs = dict["cs"]
    voltage = dict["voltage"]
    pixel_size = dict["apix"]
    b_factor = dict["bfactor"]
    ampcont = dict["ampcont"]

    ctf_2 = []
    scl = old_div(old_div(1.0, pixel_size) , nx)
    length = int(1.7321 * float(old_div(nx, 2))) + 2
    ctf_2 = [0.0] * length
    for i in range(length):
        ctf_val = EMAN2_cppwrap.Util.tf(dz, i * scl, voltage, cs, ampcont, b_factor)
        ctf_2[i] = ctf_val * ctf_val
    return ctf_2


def ctf_img(nx, ctf, sign=1, ny=0, nz=1):
    """
		Generate a 1-2-3-D complex image containing the CTF.
	 	Default is 2D output.
	  	Input
			nx: x image size.
			ctf: ctf object, see CTF_info for description.
			sign: sign of the CTF.
			ny: y image size
			nz: z image size
		Output
			ctfimg: complex image containing CTF.
	"""
    dict = ctf.to_dict()
    dz = dict["defocus"]
    cs = dict["cs"]
    voltage = dict["voltage"]
    pixel_size = dict["apix"]
    b_factor = dict["bfactor"]
    ampcont = dict["ampcont"]
    dza = dict["dfdiff"]
    azz = dict["dfang"]

    if (ny < 1):  ny = nx
    return EMAN2_cppwrap.Util.ctf_img(nx, ny, nz, dz, pixel_size, voltage, cs, ampcont, b_factor, dza, azz, sign)


def ctf_img_real(nx, ctf, sign=1, ny=0, nz=1):
    """
		Generate a 1-2-3-D real image containing the CTF.
	 	Default is 2D output.
	  	Input
			nx: x image size.
			ctf: ctf object, see CTF_info for description.
			sign: sign of the CTF.
			ny: y image size
			nz: z image size
		Output
			ctfimg: real image containing CTF, x-size half of the complex
	"""
    dict = ctf.to_dict()
    dz = dict["defocus"]
    cs = dict["cs"]
    voltage = dict["voltage"]
    pixel_size = dict["apix"]
    b_factor = dict["bfactor"]
    ampcont = dict["ampcont"]
    dza = dict["dfdiff"]
    azz = dict["dfang"]

    if (ny < 1):  ny = nx
    return EMAN2_cppwrap.Util.ctf_img_real(nx, ny, nz, dz, pixel_size, voltage, cs, ampcont, b_factor, dza, azz, sign)


def ctf_rimg(nx, ctf, sign=1, ny=0, nz=1):
    """
		Generate a 1-2-3-D real image containing the CTF.
	 	Default is 2D output.
	  	Input
			nx: x image size.
			ctf: ctf object, see CTF_info for description.
			sign: sign of the CTF, if sign=0 compute |CTF|
			ny: y image size
			nz: z image size
		Output
			ctfimg: image containing CTF^2.
	"""
    dict = ctf.to_dict()
    dz = dict["defocus"]
    cs = dict["cs"]
    voltage = dict["voltage"]
    pixel_size = dict["apix"]
    b_factor = dict["bfactor"]
    ampcont = dict["ampcont"]
    dza = dict["dfdiff"]
    azz = dict["dfang"]

    if (ny < 1):  ny = nx
    return EMAN2_cppwrap.Util.ctf_rimg(nx, ny, nz, dz, pixel_size, voltage, cs, ampcont, b_factor, dza, azz, sign)


def ctf2_rimg(nx, ctf, sign=1, ny=0, nz=1):
    """
		Generate a 1-2-3-D real image containing the CTF^2.
	 	Default is 2D output.
	  	Input
			nx: x image size.
			ctf: ctf object, see CTF_info for description.
			sign: sign of the CTF.
			ny: y image size
			nz: z image size
		Output
			ctfimg: image containing CTF^2.
	"""
    dict = ctf.to_dict()
    dz = dict["defocus"]
    cs = dict["cs"]
    voltage = dict["voltage"]
    pixel_size = dict["apix"]
    b_factor = dict["bfactor"]
    ampcont = dict["ampcont"]
    dza = dict["dfdiff"]
    azz = dict["dfang"]

    if (ny < 1):  ny = nx
    return EMAN2_cppwrap.Util.ctf2_rimg(nx, ny, nz, dz, pixel_size, voltage, cs, ampcont, b_factor, dza, azz, sign)


def ctflimit(nx, defocus, cs, voltage, pix):
    pass  # IMPORTIMPORTIMPORT import numpy as np

    def ctfperiod(defocus, Cs, lam, freq):
        # find local "period" T by solving fourth order polynomial resulting from equation:
        #  sin( 2pi (gamma(freq) + 1) ) = sin( 2pi (gamma(freq+T) )
        cis = Cs * 1.e7
        A = 0.5 * defocus * 10000.0 * lam
        B = 0.25 * cis * lam * lam * lam
        f2 = freq * freq
        """Multiline Comment1"""
        rot = np.roots([B, 4 * B * freq, 6 * B * f2 - A, 4 * B * f2 * freq - 2 * A * freq, -1.0])
        # print np.roots([A,2*A*freq,1.0]),-freq-np.sqrt(f2/2-1./A),-freq+np.sqrt(f2-1./A)
        # print  Util.tf(defocus, freq, voltage, Cs, 10., 0.0, 1.0),Util.tf(defocus, freq+min(rot), voltage, Cs, 10., 0.0, 1.0)
        return min(abs(rot))

    n = old_div(nx, 2) + 1
    #  Width of Fourier pixel
    fwpix = old_div(old_div(1., (2 * pix)), n)

    # Fourier cycle
    fcycle = old_div(1., (2 * fwpix))
    # Fourier period
    fper = old_div(1.0 , fcycle)
    # print "Image size %6d,   pixel size  %7.4f  Width of Fourier pixel %7.5f   Fourier period  %8.5f "%(nx,pix,fwpix,fper)

    # CTF
    lam = old_div(12.398 , np.sqrt(voltage * (1022.0 + voltage)))  # All units in A
    z1 = defocus * 10000.0
    ctfper = ctfperiod(defocus, cs, lam, old_div(1., (2 * pix)))
    # print " At Nyquist, the CTF period is ",ctfper
    for ii in range(n - 1, 1, -1):
        # print ii
        xr = old_div(old_div(ii ,float(n - 1)) , (2 * pix))
        ctfper = ctfperiod(defocus, cs, lam, xr)
        # print ii,xr,ctfper
        if (ctfper > fper):
            # print  " Limiting frequency is:",xr,"  limiting resolution is:",1.0/xr
            return int(old_div(xr , fwpix) + 0.5), xr
    return old_div(nx, 2), old_div(1.0, (2 * pix))


def imf_params_cl1(pw, n=2, iswi=3, Pixel_size=1):
    """
		Extract image formation parameters using contrained simplex method
		The output is a list of list, which contains the following four elements:
		1. frequencies in 1/Angstrom
		2. fitted curve, either background noise or envelope function
		3. original power spectrum to be fitted
		4. The parameters
		Attention:
		    iswi= 2 using polynomial n rank to fit no-Gaussian envelope function
			iswi =3 using polynomial n rank to fit background
			n = the polynomial rank +1
			The optimization tend to fail when the polynomial rank is higher than 6
	"""
    feq = []
    cur = []
    parm = []
    t = EMAN2_cppwrap.Util.pw_extract(pw, n, iswi, Pixel_size)
    for i in range(len(pw)):
        cur.append(t[i * 2])
        feq.append(t[i * 2 + 1])
    npam = len(t) - 2 * len(pw)
    for i in range(npam):
        k = 2 * len(pw) + i
        parm.append(t[k])
    return [feq, cur, pw, parm]


def adaptive_mask(vol, nsigma=1.0, threshold=-9999.0, ndilation=3, edge_width=5, mode="C"):
    """
		Name
			adaptive_mask - create a mask from a given image.
		Input
			img: input image
			nsigma: value for initial thresholding of the image.
		Output
			mask: The mask will have values one, zero, with cosine smooth transition between two regions.
	"""
    pass  # IMPORTIMPORTIMPORT from utilities  import gauss_edge, model_circle
    pass  # IMPORTIMPORTIMPORT from morphology import binarize, dilation
    nx = vol.get_xsize()
    ny = vol.get_ysize()
    nz = vol.get_zsize()
    mc = sparx_utilities.model_circle( old_div(nx, 2), nx, ny, nz) - sparx_utilities.model_circle( old_div(nx, 3), nx, ny, nz)
    s1 = EMAN2_cppwrap.Util.infomask(vol, mc, True)  # flip true: find statistics under the mask (mask >0.5)
    if threshold <= -9999.0:
        # Use automatic mode
        s1 = [s1[0] + s1[1] * nsigma, s1[0], s1[1], nsigma]
    # new s1[0] is calculated threshold for binarize
    else:
        # use the user-provided threshold
        if s1[1] != 0.0:
            s1 = [threshold, s1[0], s1[1], old_div((threshold - s1[0]), s1[1])]
        else:
            s1 = [threshold, s1[0], s1[1], 0.0]
    # new s1[3] is calculated nsigma corresponding to user-provided threshold
    mask = EMAN2_cppwrap.Util.get_biggest_cluster(binarize(vol, s1[0]))
    for i in range(ndilation):   mask = dilation(mask)
    mask = EMAN2_cppwrap.Util.soft_edge(mask, edge_width, "C")
    return mask


def cosinemask(im, radius=-1, cosine_width=5, bckg=None, s=999999.0):
    """
		Apply mask with a cosine fall-off setting values outside of radius_cosine_width to the average computed outside.
		The fall-off begins from pixel at a distance radius from the center,
		i.e., mask(radius) = 1 and mask(radius+cosine_width)=0.
		if s=999999.0 using average else program takes in user-provided s
	"""
    return EMAN2_cppwrap.Util.cosinemask(im, radius, cosine_width, bckg, s)


"""Multiline Comment2"""


def get_shrink_3dmask(nxinit, mask_file_name):
    pass  # IMPORTIMPORTIMPORT from utilities import get_im
    pass  # IMPORTIMPORTIMPORT from fundamentals import resample
    pass  # IMPORTIMPORTIMPORT from morphology   import binarize
    mask3d = sparx_utilities.get_im(mask_file_name)
    nx2 = nxinit
    nx1 = mask3d.get_xsize()
    if nx1 == nx2:
        return mask3d
    else:
        shrinkage = old_div(float(nx2), nx1)
        mask3d = binarize(sparx_fundamentals.resample(mask3d, shrinkage),
                          0.5)  # added 0.5 here to fix binarization problem
        return mask3d


def get_biggest_cluster(mg):
    """
	  Input: binary image
	  Output: image that contains the largest connected subset in the input image
	  This code was written by Wei in C and put in util_sparx.cpp
	"""
    pass  # IMPORTIMPORTIMPORT from utilities import model_blank
    pass  # IMPORTIMPORTIMPORT from morphology import collapse

    nx = mg.get_xsize()
    ny = mg.get_ysize()
    nz = mg.get_zsize()

    lg = mg.copy()
    s = EMAN2_cppwrap.Util.infomask(lg, None, True)
    nnc = int(s[0] * nx * ny * nz)

    cls = sparx_utilities.model_blank(nx, ny, nz)

    l = []
    grp = 0

    # endless loop
    while (grp < 1):
        grp -= 1

        fif = True
        for k in range(nz):
            if (fif):
                for j in range(ny):
                    if (fif):
                        for i in range(nx):
                            if (fif and lg[i, j, k]):
                                lg[i, j, k] = 0
                                l.append([i, j, k])
                                fif = False

        if (fif):
            #  All points checked, extract the largest group
            grp = -1 - grp
            # check if any groups found
            if (grp == 0):  return cls
            cls *= -1
            # if one group, simply return it
            if (grp == 1):  return cls
            st = 0.0
            for ig in range(1, grp):
                s = EMAN2_cppwrap.Util.infomask(collapse(cls, ig - 0.5, ig + 0.5), None, True)
                if (s[0] > st):
                    st = s[0]
                    tig = ig
            return collapse(cls, tig - 0.5, tig + 0.5)

        while (len(l) > 0):
            cr = l[0]
            del l[0]
            cls[cr[0], cr[1], cr[2]] = grp
            lg[cr[0], cr[1], cr[2]] = 0
            for k in range(-1, 2, 1):
                kq = cr[2] + k
                if (kq > -1 and kq < nz):
                    for j in range(-1, 2, 1):
                        jq = cr[1] + j
                        if (jq > -1 and jq < ny):
                            for i in range(-1, 2, 1):
                                iq = cr[0] + i
                                if (iq > -1 and iq < nx):
                                    if (lg[iq, jq, kq]):
                                        lg[iq, jq, kq] = 0
                                        l.append([iq, jq, kq])


def compute_bfactor(pws, freq_min, freq_max, pixel_size=1.0):
    """
		Estimate B-factor from power spectrum
		pws          : 1D rotational average of power spectrum, length should be half of the image size
		idx_freq_min : the index of the minimum frequency of fitting range
		idx_freq_max : the index of the maximum frequency of fitting range
		pixel_size   :  in A
	"""
    pass  # IMPORTIMPORTIMPORT from math import log, sqrt
    pass  # IMPORTIMPORTIMPORT from statistics import linreg
    nr = len(pws)
    """Multiline Comment4"""

    pws_log = [0.0] * nr
    x = [0.0] * nr
    q = min(pws)
    for i in range(1, nr):
        pws_log[i] = numpy.log(pws[i])  # /q)
        x[i] = (old_div(old_div(float(i), (2 * nr)), pixel_size)) ** 2
    idx_freq_min = 1
    for i in range(1, nr):
        if (x[i] > freq_min ** 2):
            idx_freq_min = i
            break

    idx_freq_max = 1
    for i in range(1, nr):
        idx_freq_max = i
        if (x[i] > freq_max ** 2):
            break

    # Linear regression will crash if min & max frequencies are only one apart
    if idx_freq_max - idx_freq_min <= 1:
        sparx_global_def.ERROR(
            "B_start is too high a resolution! Decrease it (under Advanced) and re-run the program! ",
            "compute_bfactor")

    B, s = sparx_statistics.linreg(x[idx_freq_min:idx_freq_max], pws_log[idx_freq_min:idx_freq_max])
    # print  B,s

    ff = [0.0] * nr
    pass  # IMPORTIMPORTIMPORT from math import exp
    for i in range(nr):  ff[i] = B * x[i] + s

    return -B, [x, ff, pws_log], idx_freq_min, idx_freq_max


################
#
#  CTER code (new version since 2016/03/16)
#
################
#
# NOTE: 2016/03/16 Toshio Moriya
# In this version, the IO-related interface is simplified for sxcter.py and sxgui.py
# Since cter() was used in not only sxcter.py but also e2boxer.py and sxhelixboxer.py,
# This new version is added to avoid breaking e2boxer.py and sxhelixboxer.py
#
# NOTE: 2016/03/16 Toshio Moriya
# To get a single  file name from a GUI application,
# there must be a better way than using guimic...
#
# NOTE: 2016/11/16 Toshio Moriya
# Now, this function assume the MPI setup and clean up is done by caller, such as mpi_init, and mpi_finalize
#
# NOTE: 07/11/2017  PAP
#       This is "exact" copy of mrk version with a switch to amplitudes (square root of PW)
#
def cter_mrk(input_image_path, output_directory, selection_list=None, wn=512, pixel_size=-1.0, \
             Cs=2.0, voltage=300.0, wgh=10.0, f_start=-1.0, f_stop=-1.0, \
             kboot=16, overlap_x=50, overlap_y=50, edge_x=0, edge_y=0, \
             check_consistency=False, stack_mode=False, debug_mode=False, \
             program_name="cter_mrk() in morphology.py", \
             RUNNING_UNDER_MPI=False, main_mpi_proc=0, my_mpi_proc_id=0, n_mpi_procs=1):
    """Multiline Comment5"""
    pass  # IMPORTIMPORTIMPORT from   EMAN2 import periodogram
    pass  # IMPORTIMPORTIMPORT from   EMAN2db import db_check_dict, db_parse_path
    pass  # IMPORTIMPORTIMPORT from   applications import MPI_start_end
    pass  # IMPORTIMPORTIMPORT from   utilities import read_text_file, write_text_file, get_im, model_blank, model_circle, amoeba, generate_ctf
    pass  # IMPORTIMPORTIMPORT from   utilities import if_error_then_all_processes_exit_program
    pass  # IMPORTIMPORTIMPORT from   utilities import wrap_mpi_bcast
    pass  # IMPORTIMPORTIMPORT from   sys import exit
    pass  # IMPORTIMPORTIMPORT import numpy as np
    pass  # IMPORTIMPORTIMPORT import os
    pass  # IMPORTIMPORTIMPORT import glob
    pass  # IMPORTIMPORTIMPORT from   fundamentals import tilemic, rot_avg_table, resample
    pass  # IMPORTIMPORTIMPORT from   morphology   import threshold, bracket_def, bracket, goldsearch_astigmatism
    pass  # IMPORTIMPORTIMPORT from   morphology   import defocus_baseline_fit, simpw1d, movingaverage, localvariance, defocusgett
    pass  # IMPORTIMPORTIMPORT from   morphology   import defocus_guessn, defocusget_from_crf, make_real
    pass  # IMPORTIMPORTIMPORT from   morphology   import fastigmatism, fastigmatism1, fastigmatism2, fastigmatism3, simctf, simctf2, simctf2out, fupw,ctf2_rimg
    pass  # IMPORTIMPORTIMPORT from   alignment    import Numrinit, ringwe
    pass  # IMPORTIMPORTIMPORT from   statistics   import table_stat
    pass  # IMPORTIMPORTIMPORT from   pixel_error  import angle_ave
    pass  # IMPORTIMPORTIMPORT from   inspect      import currentframe, getframeinfo
    pass  # IMPORTIMPORTIMPORT from   global_def   import ERROR
    pass  # IMPORTIMPORTIMPORT import global_def
    pass  # IMPORTIMPORTIMPORT from   time import time
    pass  # IMPORTIMPORTIMPORT from   mpi import MPI_COMM_WORLD, mpi_barrier

    # ====================================================================================
    # Prepare processing
    # ====================================================================================

    # ------------------------------------------------------------------------------------
    # Find the CTER Running Mode before checking error conditions
    # ------------------------------------------------------------------------------------
    i_enum = -1;
    idx_cter_mode_invalid = i_enum;
    i_enum += 1;
    idx_cter_mode_all_mics = i_enum  # All Micrographs Mode - Process all s in a directory
    i_enum += 1;
    idx_cter_mode_selected_mics = i_enum  # Selected Micrographs Mode - Process all s in a selection list file
    i_enum += 1;
    idx_cter_mode_single_mic = i_enum  # Single Micrograph Mode - Process a single
    i_enum += 1;
    idx_cter_mode_stack = i_enum  # Stack Mode - Process a stack (Advanced Option)
    i_enum += 1;
    idx_cter_mode_counts = i_enum

    cter_mode_idx = idx_cter_mode_invalid
    cter_mode_name = None
    if stack_mode == False:
        # One of three Micrograph Modes
        # For any of Micrograph Modes, input image file name must be a file path pattern containing wild card "*"
        if selection_list == None:
            # User did not use selection list option
            # -> All Micrographs Mode
            cter_mode_idx = idx_cter_mode_all_mics
            cter_mode_name = "All Micrographs Mode"
        else:
            if os.path.splitext(selection_list)[1] == ".txt":
                # User specified a selection list text file path containing".txt" extension through selection list option
                # -> Selected Micrographs Mode
                cter_mode_idx = idx_cter_mode_selected_mics
                cter_mode_name = "Selected Micrographs Mode"
            else:
                # User specified an image file path (a non-text file path) through selection list option
                # -> Single Micrograph Mode
                cter_mode_idx = idx_cter_mode_single_mic
                cter_mode_name = "Single Micrograph Mode"
    else:
        # (Particle) Stack Mode
        cter_mode_idx = idx_cter_mode_stack
        cter_mode_name = "Stack Mode"

    if my_mpi_proc_id == main_mpi_proc:
        print(" ")
        print(("----- Running with %s -----" % (cter_mode_name)))

    # ------------------------------------------------------------------------------------
    # Check mode-dependent error conditions of input arguments and options if abort is necessary. All nodes do this checking
    # ------------------------------------------------------------------------------------
    error_message_list = []  # List of error messages. If no error is found, the length should be zero
    if not stack_mode:

        # Check error conditions applicable to any of Micrograph Mode
        if input_image_path.find("*") == -1:
            error_message_list.append(
                "Input image file path (%s) for %s must be a  path pattern containing wild card (*). Please check input_image_path argument." % (
                input_image_path, cter_mode_name))

        if input_image_path[:len("bdb:")].lower() == "bdb:":
            error_message_list.append(
                "BDB file can not be selected as input image file path (%s) for %s. Please check input_image_path argument and convert the image format." % (
                input_image_path, cter_mode_name))

        # Check error conditions applicable to Selected Micrographs Mode
        if cter_mode_idx == idx_cter_mode_selected_mics:
            if not os.path.exists(selection_list):
                error_message_list.append(
                    "Selection list text file specified by selection_list option (%s) for %s does not exists. Please check selection_list option." % (
                    selection_list, cter_mode_name))

        if cter_mode_idx == idx_cter_mode_single_mic:
            if not os.path.exists(os.path.join(os.path.dirname(input_image_path), os.path.basename(selection_list))):
                error_message_list.append(
                    "Micrograph specified by selection_list option (%s) for %s does not exist. Please check selection_list option." % (
                    selection_list, cter_mode_name))
            #
            if RUNNING_UNDER_MPI and n_mpi_procs != 1:
                error_message_list.append(
                    "%s supports only a single processor version. Please change MPI settings." % (cter_mode_name))

    else:
        # Check error conditions
        if input_image_path.find("*") != -1:
            error_message_list.append(
                "Stack file path specified by input_image_path (%s) for %s should not contain wild card (*). Please check input_image_path argument." % (
                input_image_path, cter_mode_name))

        is_not_found_input_image_file = False
        if input_image_path[:len("bdb:")].lower() == "bdb:":
            if not EMAN2db.db_check_dict(input_image_path):
                is_not_found_input_image_file = True
        else:
            if not os.path.exists(input_image_path):
                is_not_found_input_image_file = True
        if is_not_found_input_image_file:
            error_message_list.append(
                "Stack file specified by input_image_path (%s) for %s does not exist. Please check input_image_path argument." % (
                input_image_path, cter_mode_name))

        if RUNNING_UNDER_MPI and n_mpi_procs != 1:
            error_message_list.append(
                "%s supports only a single processor version. Please change MPI settings." % (cter_mode_name))

    # --------------------------------------------------------------------------------
    # check output-related error conditions (mode-independent). All nodes do this checking
    # --------------------------------------------------------------------------------
    if os.path.exists(output_directory):
        error_message_list.append(
            "Output directory (%s) exists already. Please check output_directory argument." % (output_directory))

    # --------------------------------------------------------------------------------
    # Check error conditions of options (mode-independent). All nodes do this checking
    # --------------------------------------------------------------------------------
    if pixel_size <= 0.0:
        error_message_list.append(
            "Pixel size (%f) must not be negative. Please set a pasitive value larger than 0.0 to pixel_size option." % (
                pixel_size))

    if wn <= 0.0:
        error_message_list.append(
            "CTF window size (%d) must not be negative. Please set a valid value larger than 0 to wn option." % (wn))

    # --------------------------------------------------------------------------------
    # Print all error messages and abort the process if necessary.
    # --------------------------------------------------------------------------------
    error_status = None
    if len(error_message_list) > 0:
        # Detected error! Print all error messages
        if my_mpi_proc_id == main_mpi_proc:
            print(" ")
            for error_message in error_message_list:
                print(("ERROR!!! %s" % (error_message)))
        error_status = ("Detected %d error(s) related to arguments and options. Run %s -h for help. Exiting..." % (
        len(error_message_list), program_name), inspect.getframeinfo(inspect.currentframe()))
    sparx_utilities.if_error_then_all_processes_exit_program(error_status)
    if RUNNING_UNDER_MPI:
        # Wait for all mpi processes to check error conditions, especially existence of output directory
        # Without this barrier, main mpi process can create output directory before some child mpi process check this error.
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    del error_message_list  # Don't need this anymore

    # ------------------------------------------------------------------------------------
    # Check warning conditions of options
    # ------------------------------------------------------------------------------------
    if my_mpi_proc_id == main_mpi_proc:
        if stack_mode:
            if selection_list != None:
                print(" ")
                print(("WARNING!!! --selection_list option will be ignored in %s." % (cter_mode_name)))
            if wn != 512:
                print(" ")
                print(("WARNING!!! --wn option will be ignored in %s." % (cter_mode_name)))
            if overlap_x != 50:
                print(" ")
                print(("WARNING!!! --overlap_x option will be ignored in %s." % (cter_mode_name)))
            if overlap_y != 50:
                print(" ")
                print(("WARNING!!! --overlap_y option will be ignored in %s." % (cter_mode_name)))
            if edge_x != 0:
                print(" ")
                print(("WARNING!!! --edge_x option will be ignored in %s." % (cter_mode_name)))
            if edge_y != 0:
                print(" ")
                print(("WARNING!!! --edge_y option will be ignored in %s." % (cter_mode_name)))
            if check_consistency:
                print(" ")
                print(("WARNING!!! --check_consistency option will be ignored in %s." % (cter_mode_name)))

    # ====================================================================================
    # Create the input file path list and also check input-related error conditions if abort is necessary.
    # ====================================================================================
    input_file_path_list = []
    if not stack_mode:
        # --------------------------------------------------------------------------------
        # Prepare the variables for all sections in  mode case
        # --------------------------------------------------------------------------------
        # Micrograph basename pattern (directory path is removed from  path pattern)
        mic_pattern = input_image_path
        mic_basename_pattern = os.path.basename(mic_pattern)

        # Global entry dictionary (all possible entries from all lists) for all mic id substring
        global_entry_dict = {}  # mic id substring is the key
        subkey_input_mic_path = "Input Micrograph Path"
        subkey_selected_mic_basename = "Selected Micrograph Basename"

        # List keeps only id substrings of s whose all necessary information are available
        valid_mic_id_substr_list = []

        # --------------------------------------------------------------------------------
        # Obtain the list of  id sustrings using a single CPU (i.e. main mpi process)
        # --------------------------------------------------------------------------------
        # NOTE: Toshio Moriya 2016/11/15
        # The below is not a real while.
        # It gives if-statements an opportunity to use break when errors need to be reported
        # However, more elegant way is to use 'raise' statement of exception mechanism...
        #
        error_status = None
        while my_mpi_proc_id == main_mpi_proc:

            # --------------------------------------------------------------------------------
            # Prepare variables for this section
            # --------------------------------------------------------------------------------
            # Prefix and suffix of  basename pattern
            # to find the head/tail indices of  id substring
            mic_basename_tokens = mic_basename_pattern.split('*')
            # Find head index of  id substring
            mic_id_substr_head_idx = len(mic_basename_tokens[0])

            # --------------------------------------------------------------------------------
            # Register  id substrings found in the input directory (specified by  path pattern)
            # to the global entry dictionary
            # --------------------------------------------------------------------------------
            # Generate the list of  paths in the input directory
            print(" ")
            print("Checking the input directory...")
            input_mic_path_list = glob.glob(mic_pattern)
            # Check error condition of input  file path list
            print(("Found %d micrographs in %s." % (len(input_mic_path_list), os.path.dirname(mic_pattern))))
            if len(input_mic_path_list) == 0:
                # The result shouldn't be empty if the specified  file name pattern is invalid
                error_status = (
                "There are no micrographs whose paths match with the specified file path pattern (%s) for %s. Please check input_image_path. Run %s -h for help." % (
                mic_pattern, cter_mode_name, program_name), inspect.getframeinfo(inspect.currentframe()))
                break

            # Register  id substrings to the global entry dictionary
            for input_mic_path in input_mic_path_list:
                # Find tail index of  id substring and extract the substring from the  name
                input_mic_basename = os.path.basename(input_mic_path)
                mic_id_substr_tail_idx = input_mic_basename.index(mic_basename_tokens[1])
                mic_id_substr = input_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
                if not mic_id_substr in global_entry_dict:
                    # print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from input_mic_path_list " % (mic_id_substr))
                    global_entry_dict[mic_id_substr] = {}
                global_entry_dict[mic_id_substr][subkey_input_mic_path] = input_mic_path

            # --------------------------------------------------------------------------------
            # Register  id substrings found in the selection list
            # to the global entry dictionary
            # --------------------------------------------------------------------------------
            # Generate the list of selected  paths in the selection file
            selected_mic_path_list = []
            # Generate  lists according to the execution mode
            if cter_mode_idx == idx_cter_mode_all_mics:
                # Treat all s in the input directory as selected ones
                selected_mic_path_list = input_mic_path_list
            else:
                if os.path.splitext(selection_list)[1] == ".txt":
                    print(" ")
                    print("Checking the selection list...")
                    selected_mic_path_list = sparx_utilities.read_text_file(selection_list)

                    # Check error condition of  entry lists
                    print(("Found %d microgarph entries in %s." % (len(selected_mic_path_list), selection_list)))
                    if len(selected_mic_path_list) == 0:
                        error_status = (
                        "The provided  list file (%s) for %s mode contains no entries. Please check selection_list option and make sure the file contains a  list. Run %s -h for help." % (
                        selection_list, cter_mode_name, program_name), inspect.getframeinfo(inspect.currentframe()))
                        break
                else:
                    print(" ")
                    print(("Processing a single micorgprah: %s..." % (selection_list)))
                    selected_mic_path_list = [selection_list]

                selected_mic_directory = os.path.dirname(selected_mic_path_list[0])
                if selected_mic_directory != "":
                    print(("    NOTE: Program disregards the directory paths in the selection list (%s)." % (
                        selected_mic_directory)))

            # Register  id substrings to the global entry dictionary
            for selected_mic_path in selected_mic_path_list:
                # Find tail index of  id substring and extract the substring from the  name
                selected_mic_basename = os.path.basename(selected_mic_path)
                mic_id_substr_tail_idx = selected_mic_basename.index(mic_basename_tokens[1])
                mic_id_substr = selected_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
                if not mic_id_substr in global_entry_dict:
                    # print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from selected_mic_path_list " % (mic_id_substr))
                    global_entry_dict[mic_id_substr] = {}
                global_entry_dict[mic_id_substr][subkey_selected_mic_basename] = selected_mic_basename

            # --------------------------------------------------------------------------------
            # Clean up variables related to registration to the global entry dictionary
            # --------------------------------------------------------------------------------
            del mic_basename_tokens
            del mic_id_substr_head_idx

            # --------------------------------------------------------------------------------
            # Create the list containing only valid  id substrings
            # --------------------------------------------------------------------------------
            # Prepare lists to keep track of invalid (rejected) s
            no_input_mic_id_substr_list = []

            print(" ")
            print("Checking the input datasets consistency...")

            # Loop over substring id list
            for mic_id_substr in global_entry_dict:
                mic_id_entry = global_entry_dict[mic_id_substr]

                warinnig_messages = []
                # selected  basename must have been registed always .
                if subkey_selected_mic_basename in mic_id_entry:
                    # Check if associated input  exists
                    if not subkey_input_mic_path in mic_id_entry:
                        input_mic_path = mic_pattern.replace("*", mic_id_substr)
                        warinnig_messages.append("    associated input  %s." % (input_mic_path))
                        no_input_mic_id_substr_list.append(mic_id_substr)

                    if len(warinnig_messages) > 0:
                        print(("WARNING!!! Micrograph ID %s does not have:" % (mic_id_substr)))
                        for warinnig_message in warinnig_messages:
                            print(warinnig_message)
                        print("    Ignores this as an invalid entry.")
                    else:
                        # print("MRK_DEBUG: adding mic_id_substr := ", mic_id_substr)
                        valid_mic_id_substr_list.append(mic_id_substr)
            # 	# This entry is not in the selection list. Do nothing

            # Check the input dataset consistency and save the result to a text file, if necessary.
            if check_consistency:
                # Create output directory
                os.mkdir(output_directory)

                # Open the consistency check file
                inconsist_mic_list_path = os.path.join(output_directory, "inconsist_mic_id_file.txt")
                print(" ")
                print(("Generating the input datasets consistency report in %s..." % (inconsist_mic_list_path)))
                inconsist_mic_list_file = open(inconsist_mic_list_path, "w")
                inconsist_mic_list_file.write("# The information about inconsistent  IDs\n")
                # Loop over substring id list
                for mic_id_substr in global_entry_dict:
                    mic_id_entry = global_entry_dict[mic_id_substr]

                    consistency_messages = []
                    # Check if associated input  path exists
                    if not subkey_input_mic_path in mic_id_entry:
                        input_mic_path = mic_pattern.replace("*", mic_id_substr)
                        consistency_messages.append("    associated input  %s." % (input_mic_path))

                    # Check if associated selected  basename exists
                    if not subkey_selected_mic_basename in mic_id_entry:
                        input_mic_path = mic_pattern.replace("*", mic_id_substr)
                        consistency_messages.append("    associated selected  %s." % (input_mic_path))

                    if len(consistency_messages) > 0:
                        inconsist_mic_list_file.write("Micrograph ID %s does not have:\n" % (mic_id_substr))
                        for consistency_message in consistency_messages:
                            inconsist_mic_list_file.write(consistency_message)
                            inconsist_mic_list_file.write("\n")

                # Close the consistency check file, if necessary
                inconsist_mic_list_file.flush()
                inconsist_mic_list_file.close()

            # Since mic_id_substr is once stored as the key of global_entry_dict and extracted with the key order
            # we need sort the valid_mic_id_substr_list here
            if debug_mode: print(("BEFORE SORT: valid_mic_id_substr_list := ", valid_mic_id_substr_list))
            valid_mic_id_substr_list.sort(key=str.lower)  # Sort list of  IDs using case insensitive string comparison
            if debug_mode: print(("AFTER SORT: valid_mic_id_substr_list := ", valid_mic_id_substr_list))

            # --------------------------------------------------------------------------------
            # Print out the summary of input consistency
            # --------------------------------------------------------------------------------
            print(" ")
            print("Summary of dataset consistency check...")
            print(("  Detected  IDs               : %6d" % (len(global_entry_dict))))
            print(("  Entries in input directory  : %6d" % (len(input_mic_path_list))))
            print(("  Entries in selection list   : %6d" % (len(selected_mic_path_list))))
            print(("  Rejected by no input        : %6d" % (len(no_input_mic_id_substr_list))))
            print(("  Valid Entries               : %6d" % (len(valid_mic_id_substr_list))))

            # --------------------------------------------------------------------------------
            # Check MPI error condition
            # --------------------------------------------------------------------------------
            if len(valid_mic_id_substr_list) < n_mpi_procs:
                error_status = (
                "Number of MPI processes (%d) supplied by --np in mpirun cannot be greater than %d (number of valid s that satisfy all criteria to be processed). Run %s -h for help." % (
                n_mpi_procs, len(valid_mic_id_substr_list, program_name)), inspect.getframeinfo(inspect.currentframe()))
                break

            # --------------------------------------------------------------------------------
            # Create input file path list
            # --------------------------------------------------------------------------------
            for mic_id_substr in valid_mic_id_substr_list:
                mic_path = global_entry_dict[mic_id_substr][subkey_input_mic_path]
                input_file_path_list.append(mic_path)

            # --------------------------------------------------------------------------------
            # Clean up variables related to tracking of invalid (rejected) s
            # --------------------------------------------------------------------------------
            del input_mic_path_list
            del selected_mic_path_list
            del no_input_mic_id_substr_list

            break

        # --------------------------------------------------------------------------------
        # Clean up the variables for all sections in  mode case
        # --------------------------------------------------------------------------------
        del mic_pattern
        del mic_basename_pattern
        del global_entry_dict
        del subkey_input_mic_path
        del subkey_selected_mic_basename
        del valid_mic_id_substr_list

        # --------------------------------------------------------------------------------
        # Print all error messages and abort the process if necessary.
        # --------------------------------------------------------------------------------
        # NOTE: Toshio Moriya 2016/11/15
        # The following function takes care of the case when an if-statement uses break for occurence of an error.
        # However, more elegant way is to use 'exception' statement of exception mechanism...
        #
        sparx_utilities.if_error_then_all_processes_exit_program(error_status)

    else:
        input_file_path_list.append(input_image_path)

    if RUNNING_UNDER_MPI:
        # Wait for main mpi process to create the input file path list
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

        # All mpi processes should know input file path list
        input_file_path_list = sparx_utilities.wrap_mpi_bcast(input_file_path_list, main_mpi_proc)

    # ====================================================================================
    # Prepare input file path(s)
    # ====================================================================================
    #
    # NOTE: 2016/03/17 Toshio Moriya
    # From here on, stack (and namics) will be used to distinguish stack mode and  mode.
    # However, a single input_file_path_list should be sufficient since we already know the mode.
    # Let's consider this refactoring in the future.
    #
    stack = None  # (particle) stack file name: if it is not None, cter runs with stack mode. Otherwise, runs with  mode
    namics = []  # file name list
    if not stack_mode:
        namics = input_file_path_list
        if debug_mode: print(("BEFORE SORT: namics := ", namics))
        namics.sort(key=str.lower)  # Sort list of s using case insensitive string comparison
        if debug_mode: print(("AFTER SORT: namics := ", namics))
    else:
        stack = input_file_path_list[0]

    del input_file_path_list  # Don't need this anymore

    # Make output directory
    outpwrot = "%s/pwrot" % (output_directory)
    if stack == None:
        outmicthumb = "%s/micthumb" % (output_directory)
    if debug_mode:
        outravg = "%s/ravg" % (output_directory)
    if my_mpi_proc_id == main_mpi_proc:
        # Make output directory
        if not os.path.exists(output_directory):
            os.mkdir(output_directory)
        os.mkdir(outpwrot)
        if stack == None:
            os.mkdir(outmicthumb)
        if debug_mode:
            os.mkdir(outravg)

    if RUNNING_UNDER_MPI:
        # Make all mpi processes wait for main mpi process to create output directory
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    # Set up loop variables depending on the cter mode
    if stack == None:
        if RUNNING_UNDER_MPI:
            set_start, set_end = sparx_applications.MPI_start_end(len(namics), n_mpi_procs, my_mpi_proc_id)
        else:
            set_start = 0
            set_end = len(namics)
    else:
        pw2 = []
        data = EMAN2_cppwrap.EMData.read_images(stack)
        nima = len(data)
        for i in range(nima):
            pw2.append(EMAN2_cppwrap.periodogram(data[i]))
        wn = pw2[0].get_xsize()
        set_start = 0
        set_end = 1

    # Set up progress message
    if my_mpi_proc_id == main_mpi_proc:
        print(" ")
        print("Estimating CTF parameters...")
        if stack == None:
            print("  Micrographs processed by main process (including percent of progress):")
            progress_percent_step =  old_div((
                                                set_end - set_start) , 100.0)  # the number of micrograms for main mpi processer divided by 100

    totresi = []
    missing_img_names = []
    rejected_img_names = []
    for ifi in range(set_start, set_end):
        img_type = ""
        img_name = ""
        img_basename_root = ""

        if stack == None:
            img_type = "Micrograph"
            img_name = namics[ifi]

            if my_mpi_proc_id == main_mpi_proc:
                print(("    Processing %s ---> %6.2f%%" % (img_name, old_div((ifi - set_start), progress_percent_step))))

            if not os.path.exists(img_name):
                missing_img_names.append(img_name)
                print(
                    "    %s %s: Can not find this file. Skipping the estimation and CTF parameters are not stored..." % (
                    img_type, img_name))
                continue

            mic = sparx_utilities.get_im(img_name)
            try:
                pw2 = sparx_fundamentals.tilemic(mic, win_size=wn, overlp_x=overlap_x, overlp_y=overlap_y,
                                                 edge_x=edge_x, edge_y=edge_y)
            except:
                print(
                    "MRK_DEBUG: tilemic() in cter_mrk() raised an exception. The micrographs {} might have a problem. Please check it and remove it if necessary.".format(
                        img_name))
                raise
            del mic
        else:
            img_type = "Stack"
            img_name = stack

            numFM = EMAN2_cppwrap.EMUtil.get_image_count(img_name)
            pw2 = [None] * numFM
            for i in range(numFM):
                pw2.append(EMAN2_cppwrap.periodogram(sparx_utilities.get_im(img_name, i)))

        for i in range(len(pw2)):
            pw2[i] = square_root(pw2[i])
        if debug_mode: print("    %s %s: Process %04d started the processing. Detected %d image(s) in this %s file." % (
        img_type, img_name, ifi, numFM, img_type.lower()))

        if EMAN2db.db_check_dict(img_name) == False:
            img_basename_root = os.path.splitext(os.path.basename(img_name))[0]
        else:
            path, dictname, keys = EMAN2db.db_parse_path(img_name)
            img_basename_root = dictname

        nimi = len(pw2)
        adefocus = [0.0] * kboot
        aamplitu = [0.0] * kboot
        aangle = [0.0] * kboot

        allroo = []
        for imi in range(nimi):
            allroo.append(sparx_fundamentals.rot_avg_table(pw2[imi]))
        lenroo = len(allroo[0])
        # print time(),nimi

        for nboot in range(kboot):
            # at = time()
            if (nboot == 0):
                boot = list(range(nimi))
            else:
                pass  # IMPORTIMPORTIMPORT from random import randint
                for imi in range(nimi): boot[imi] = random.randint(0, nimi - 1)
            qa = sparx_utilities.model_blank(wn, wn)
            roo = np.zeros(lenroo, np.float32)
            sroo = np.zeros(lenroo, np.float32)
            aroo = np.zeros(lenroo, np.float32)

            for imi in range(nimi):
                EMAN2_cppwrap.Util.add_img(qa, pw2[boot[imi]])
                temp1 = np.array(allroo[boot[imi]])
                roo += temp1
                temp2 = movingaverage(temp1, 10)
                aroo += temp2
                sroo += temp2 ** 2
            sroo[0] = sroo[1]
            aroo[0] = aroo[1]
            sroo = old_div((sroo -  old_div(aroo ** 2 , nimi)) , nimi)
            aroo = old_div(aroo, nimi)
            roo =  old_div(roo, nimi)
            qa = old_div(qa,nimi)

            if f_start < 0:
                #  Find a break point
                bp = 1.e23
                for i in range(5, lenroo - 5):
                    # t1 = linreg(sroo[:i])
                    # t2 = linreg(sroo[i:])
                    # tt = t1[1][0] + t2[1][0]
                    xtt = np.array(list(range(i)), np.float32)
                    zet = np.poly1d(np.polyfit(xtt, sroo[:i], 2))
                    t1 = sum((sroo[:i] - zet(xtt)) ** 2)
                    xtt = np.array(list(range(i, lenroo)), np.float32)
                    zet = np.poly1d(np.polyfit(xtt, sroo[i:], 2))
                    tt = t1 + sum((sroo[i:] - zet(xtt)) ** 2)
                    if tt < bp:
                        bp = tt
                        istart = i
                # istart = 25
                # print istart
                f_start = old_div(istart , (pixel_size * wn))
            """Multiline Comment6"""
            # write_text_file([roo.tolist(),aroo.tolist(),sroo.tolist()], "sroo%03d.txt"%ifi)
            rooc = roo.tolist()

            # print namics[ifi],istart,f_start

            defc, subpw, ctf2, baseline, envelope, istart, istop = defocusgett_pap(rooc, wn, voltage=voltage,
                                                                                   Pixel_size=pixel_size, Cs=Cs,
                                                                                   ampcont=wgh, f_start=f_start,
                                                                                   f_stop=f_stop, round_off=1.0, nr1=3,
                                                                                   nr2=6, parent=None, DEBug=debug_mode)
            # defc, subpw, ctf2, baseline, envelope, istart, istop = defocusgett(rooc, wn, voltage = voltage, Pixel_size = pixel_size, Cs = Cs, ampcont = wgh, f_start = f_start, f_stop = f_stop, round_off = 1.0, nr1 = 3, nr2 = 6, parent = None, DEBug = debug_mode)
            if debug_mode:
                print("  RESULT %s" % (img_name), defc, istart, istop)

                freq = list(range(len(subpw)))
                for i in range(len(freq)):  freq[i] =  old_div(old_div(float(i), wn) , pixel_size)
                #				write_text_file([freq, subpw.tolist(), ctf2, envelope.tolist(), baseline.tolist()], "%s/ravg%05d.txt" % (output_directory, ifi))
                fou = os.path.join(outravg, "%s_ravg_%02d.txt" % (img_basename_root, nboot))
                sparx_utilities.write_text_file([freq, subpw.tolist(), ctf2, envelope.tolist(), baseline.tolist()], fou)
            # mpi_barrier(MPI_COMM_WORLD)

            # exit()
            bg = baseline.tolist()
            en = envelope.tolist()

            bckg = sparx_utilities.model_blank(wn, wn, 1, 1)
            envl = sparx_utilities.model_blank(wn, wn, 1, 1)

            pass  # IMPORTIMPORTIMPORT from math import sqrt
            nc = old_div(wn , 2)
            bg.append(bg[-1])
            en.append(en[-1])
            for i in range(wn):
                for j in range(wn):
                    r = numpy.sqrt((i - nc) ** 2 + (j - nc) ** 2)
                    ir = int(r)
                    if (ir < nc):
                        dr = r - ir
                        bckg.set_value_at(i, j, (1. - dr) * bg[ir] + dr * bg[ir + 1])
                        envl.set_value_at(i, j, (1. - dr) * en[ir] + dr * en[ir + 1])

            # qa.write_image("rs1.hdf")

            mask = sparx_utilities.model_circle(istop - 1, wn, wn) * (
                        sparx_utilities.model_blank(wn, wn, 1, 1.0) - sparx_utilities.model_circle(istart, wn, wn))
            qse = threshold((qa - bckg))  # *envl
            # (qse*mask).write_image("rs2.hdf")
            # qse.write_image("rs3.hdf")
            ##  SIMULATION
            # bang = 0.7
            # qse = ctf2_rimg(wn, generate_ctf([defc,Cs,voltage,pixel_size,0.0,wgh, bang, 37.0]) )
            # qse.write_image("rs3.hdf")

            cnx = old_div(wn, 2) + 1
            cny = cnx
            mode = "H"
            istop = min(  old_div(wn, 2) - 2, istop)  # 2-26-2015@ming
            numr = sparx_alignment.Numrinit(istart, istop, 1, mode)
            wr = sparx_alignment.ringwe(numr, mode)

            crefim = EMAN2_cppwrap.Util.Polar2Dm(qse * mask, cnx, cny, numr, mode)
            EMAN2_cppwrap.Util.Frngs(crefim, numr)
            EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)

            # pc = ctf2_rimg(wn,generate_ctf([defc,Cs,voltage,pixel_size,0.0,wgh]))
            # print ccc(pc*envl, subpw, mask)

            bang = 0.0
            bamp = 0.0
            bdef = defc
            bold = 1.e23
            while (True):
                #  in simctf2 data[3] is astigmatism amplitude
                """
				data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
				#astdata = [crefim, numr, wn, bdef, Cs, voltage, pixel_size, wgh, bang]
				for qqq in xrange(200):
					qbdef = 1.0 + qqq*0.001
					print " VALUE AT THE BEGGINING OF while LOOP  ",qbdef,simctf2(qbdef, data)#,fastigmatism3(bamp,astdata)
				"""
                """
				bamp = 0.7
				bang = 37.0

				data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
				astdata = [crefim, numr, wn, bdef, Cs, voltage, Pixel_size, wgh, bang]
				print " VALUE AT THE BEGGINING OF while LOOP  ",bdef,bamp,bang,simctf2(bdef, data),fastigmatism3(bamp,astdata,mask)
				#print  simctf2out(1.568,data)
				#exit()

				for kdef in xrange(14000,17000,10):
					dz = kdef/10000.0
					ard = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
					#print ard
					aqd = [crefim, numr, wn, dz, Cs, voltage, pixel_size, wgh, bang]
					#print aqd
					print  dz,simctf2(dz,ard),fastigmatism3(bamp,aqd,mask)
					#print aqd[-1]
				exit()
				"""
                data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
                h = 0.05 * bdef
                amp1, amp2 = bracket_def(simctf2_pap, data, bdef * 0.9, h)
                # amp1, amp2 = bracket_def(simctf2, data, bdef * 0.9, h)
                # print "bracketing of the defocus  ",amp1, amp2
                # print " ttt ",time()-srtt
                # print "bracketing of the defocus  ",amp1,amp2,simctf2(amp1, data),simctf2(amp2, data),h
                amp1, val2 = goldsearch_astigmatism(simctf2_pap, data, amp1, amp2, tol=1.0e-3)
                # amp1, val2 = goldsearch_astigmatism(simctf2, data, amp1, amp2, tol = 1.0e-3)
                # print "golden defocus ",amp1, val2,simctf2(amp1, data)
                # bdef, bcc = goldsearch_astigmatism(simctf2, data, amp1, amp2, tol=1.0e-3)
                # print "correction of the defocus  ",bdef,bcc
                # print " ttt ",time()-srtt
                """Multiline Comment7"""

                astdata = [crefim, numr, wn, bdef, Cs, voltage, pixel_size, wgh, bang, mask]
                h = 0.01
                amp1, amp2 = bracket(fastigmatism3_pap, astdata, h)
                # amp1, amp2 = bracket(fastigmatism3, astdata, h)
                # print "  astigmatism bracket  ",amp1,amp2,astdata[-1]
                # print " ttt ",time()-srtt
                bamp, bcc = goldsearch_astigmatism(fastigmatism3_pap, astdata, amp1, amp2, tol=1.0e-3)
                junk = fastigmatism3_pap(bamp, astdata)
                # bamp, bcc = goldsearch_astigmatism(fastigmatism3, astdata, amp1, amp2, tol = 1.0e-3)
                # junk = fastigmatism3(bamp,astdata)
                bang = astdata[8]

                # print astdata[8]
                # print  fastigmatism3(0.0,astdata)
                # print astdata[8]
                # temp = 0.0
                # print bdef, Cs, voltage, pixel_size, temp, wgh, bamp, bang, -bcc
                # data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
                # astdata = [crefim, numr, wn, bdef, Cs, voltage, pixel_size, wgh, bang]
                # print " VALUE WITHIN the while LOOP  ",bdef,bamp,bang,simctf2(bdef, data),fastigmatism3(bamp,astdata)
                # print "  golden search ",bamp,data[-1], fastigmatism3(bamp,data), fastigmatism3(0.0,data)
                # print " ttt ",time()-srtt
                # bamp = 0.5
                # bang = 277

                dama = sparx_utilities.amoeba([bdef, bamp], [0.2, 0.2], fupw_pap, 1.e-4, 1.e-4, 500, astdata)
                # dama = amoeba([bdef, bamp], [0.2, 0.2], fupw, 1.e-4, 1.e-4, 500, astdata)
                if debug_mode:  print("AMOEBA    ", dama)
                bdef = dama[0][0]
                bamp = dama[0][1]
                astdata = [crefim, numr, wn, bdef, Cs, voltage, pixel_size, wgh, bang, mask]
                junk = fastigmatism3_pap(bamp, astdata)
                # junk = fastigmatism3(bamp, astdata)
                bang = astdata[8]
                if debug_mode:  print(" after amoeba ", bdef, bamp, bang)
                #  The looping here is blocked as one shot at amoeba is good enough.  To unlock it, remove - from bold.
                if (bcc < -bold):
                    bold = bcc
                else:
                    break

            # data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
            # print " VALUE AFTER the while LOOP  ",bdef,bamp,bang,simctf2(bdef, data),fastigmatism3(bamp,astdata)
            # temp = 0.0
            # print ifi,bdef, Cs, voltage, pixel_size, temp, wgh, bamp, bang, -bcc
            # freq = range(len(subpw))
            # for i in xrange(len(freq)):  freq[i] = float(i)/wn/pixel_size
            # ctf2 = ctf_2(wn, generate_ctf([bdef,Cs,voltage,pixel_size,0.0,wgh]))[:len(freq)]
            # write_text_file([freq, subpw.tolist(), ctf2, envelope.tolist(), baseline.tolist()],"ravg/ravg%05d.txt"%ifi)
            # print " >>>> ",wn, bdef, bamp, Cs, voltage, pixel_size, wgh, bang
            # data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
            # print  simctf2out(bdef, data)
            # exit()
            adefocus[nboot] = bdef
            aamplitu[nboot] = bamp
            aangle[nboot] = bang
        # from sys import exit
        # exit()

        # print " ttt ",time()-srtt
        # from sys import exit
        # exit()
        ad1, ad2, ad3, ad4 = sparx_statistics.table_stat(adefocus)  # return values: average, variance, minimum, maximum
        if ad2 <= 0.0:
            print((
                              "    %s %s: Detected the variance less than zero (defocus statistics: avg = %f, var = %f, min = %f, max = %f)." % (
                      img_type, img_name, ad1, ad2, ad3, ad4)))
            print(("           The program ignores this estimate..."))
            continue

        reject = []
        thr = 3 * numpy.sqrt(ad2)
        for i in range(len(adefocus)):
            if (abs(adefocus[i] - ad1) > thr):
                print((
                                  "    %s %s: Rejected an outlier defocus estimate (defocus = %f, average defocus = %f, threshold = %f)." % (
                          img_type, img_name, adefocus[i], ad1, thr)))
                reject.append(i)

        if (len(reject) > 0):
            print(("    %s %s: Total number of rejects %s" % (img_type, img_name, len(reject))))
            for i in range(len(reject) - 1, -1, -1):
                del adefocus[i]
                del aamplitu[i]
                del aangle[i]

        if (len(adefocus) < 2):
            print((
                              "    %s %s: After rejection of outliers, there is too few estimated defocus values. Skipping the estimation and CTF parameters are not stored..." % (
                      img_type, img_name)))
        else:
            # print "adefocus",adefocus
            # print  "aamplitu",aamplitu
            # print "aangle",aangle
            ad1, ad2, ad3, ad4 = sparx_statistics.table_stat(adefocus)
            bd1, bd2, bd3, bd4 = sparx_statistics.table_stat(aamplitu)
            cd1, cd2 = sparx_pixel_error.angle_ave(
                [2 * q for q in aangle])  # Have to use this trick as the function works for range [0,360]
            cd1 = old_div(cd1,2)
            cd2 = old_div(cd2,2)
            temp = 0.0
            stdavad1 = np.sqrt(kboot * max(0.0, ad2))
            stdavbd1 = np.sqrt(kboot * max(0.0, bd2))
            cd2 *= np.sqrt(kboot)

            # Adjust value ranges of astig. amp. and angle.
            if (bd1 < 0.0):
                bd1 = -bd1
                cd1 += 90.0
            cd1 = cd1 % 180

            #  SANITY CHECK, do not produce anything if defocus abd astigmatism amplitude are out of whack
            reject_img_messages = []
            try:
                pwrot2 = rotavg_ctf(sparx_utilities.model_blank(wn, wn), ad1, Cs, voltage, pixel_size, bd1, cd1)
            except:
                reject_img_messages.append(
                    "    - Astigmatism amplitude (%f) is larger than defocus (%f) or defocus (%f) is negative." % (
                    bd1, ad1, ad1))

            if len(reject_img_messages) > 0:
                rejected_img_names.append(img_name)
                print("    %s %s: Rejected the CTF estimate - " % (img_type, img_name), ad1, Cs, voltage, pixel_size,
                      wgh, bd1, cd1, "(def, Cs, vol, apix, amp_contrast, astig_amp, astig_angle)")
                print("    %s %s: because... " % (img_type, img_name))
                for reject_img_message in reject_img_messages:
                    print(reject_img_message)
                print("    %s %s: Skipping the estimation and CTF parameters are not stored..." % (img_type, img_name))
            else:
                #  Estimate the point at which (sum_errordz ctf_1(dz+errordz))^2 falls to 0.5
                import random as rqt

                supe = sparx_utilities.model_blank(wn, wn)
                niter = 1000
                for it in range(niter):
                    EMAN2_cppwrap.Util.add_img(supe,
                                               EMAN2_cppwrap.Util.ctf_rimg(wn, wn, 1, ad1 + rqt.gauss(0.0, stdavad1),
                                                                           pixel_size, voltage, Cs, 0.0, wgh,
                                                                           bd1 + rqt.gauss(0.0, stdavbd1),
                                                                           cd1 + rqt.gauss(0.0, cd2), 1))
                ni = old_div(wn, 2)
                supe = old_div(supe,niter)
                pwrot2 = rotavg_ctf(supe, ad1, Cs, voltage, pixel_size, bd1, cd1)
                for i in range(ni):  pwrot2[i] = pwrot2[i] ** 2

                ibec = 0
                for it in range(ni - 1, 0, -1):
                    if pwrot2[it] > 0.5:
                        ibec = it
                        break
                pass  # IMPORTIMPORTIMPORT from morphology import ctf_1d
                ct = sparx_utilities.generate_ctf([ad1, Cs, voltage, pixel_size, temp, wgh, 0.0, 0.0])
                cq = ctf_1d(wn, ct)

                supe = [0.0] * ni
                niter = 1000
                for i in range(niter):
                    cq = sparx_utilities.generate_ctf(
                        [ad1 + rqt.gauss(0.0, stdavad1), Cs, voltage, pixel_size, 0.0, wgh, 0.0, 0.0])
                    ci = ctf_1d(wn, cq)[:ni]
                    for l in range(ni):  supe[l] += ci[l]

                for l in range(ni):  supe[l] = (old_div(supe[l] , niter)) ** 2

                ib1 = 0
                for it in range(ni - 1, 0, -1):
                    if supe[it] > 0.5:
                        ib1 = it
                        break
                ibec = old_div(ibec, (pixel_size * wn))  # with astigmatism
                ib1 = old_div(ib1, (pixel_size * wn))  # no astigmatism
                # from utilities import write_text_file
                # write_text_file([range(ni), supe[:ni],pwrot2[:ni]],"fifi.txt")

                # Compute defocus CV and astig. amp. CV (CV: coefficient of variation; ratio of error (SD) relative to average (mean))
                cvavad1 = old_div(stdavad1, ad1) * 100  # use percentage

                cvavbd1 = old_div(stdavbd1, bd1) * 100  # use percentage

                # Compute CTF limit (theoretical resolution limit based on the oscillations of CTF)
                # For output, use ctflim (relative frequency limit [1/A]), not ctflim_abs (absolute frequency limit)
                #
                # NOTE: 2016/03/23 Toshio Moriya
                # xr is limiting frequency [1/A]. Max is Nyquist frequency = 1.0/(2*apix[A/pixel]). <UNIT: [1/(A/pixel)/[pixel])] => [(pixel)/(A*pixel] => [1/A]>
                # 1.0/xr is limiting period (Angstrom resolution) [A]. Min is Nyquist period = (2*apix[A/pixel]). <UNIT: [1/(1/A)] = [A]>
                # fwpix is width of Fourier pixel [pixel/A] := 1.0[pixel]/(2*apix[A/pixel])/box_half[pixel] = 1[pixel]/fullsize[A]). <UNIT: [pixel/(A/pixel)/(pixel)] = [pixel*(pixel/A)*(1/pixel) = [pixel/A]>
                # int(xr/fwpix+0.5) is limiting_absolute_frequency [1/pixel]. <Unit:[(1/A)/(pixel/A)] = [(1/A)*(A/pixel)] = [1/pixel]>
                # return  int(xr/fwpix+0.5),xr, which is limiting_abs_frequency [1/pixel], and Limiting_frequency[1/A]
                #
                ctflim_abs, ctflim = ctflimit(wn, ad1, Cs, voltage, pixel_size)

                """Multiline Comment8"""
                lnsb = len(subpw)
                try:
                    crot1 = rotavg_ctf(
                        ctf_rimg(wn, sparx_utilities.generate_ctf([ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1]),
                                 sign=0), ad1, Cs, voltage, pixel_size, 0.0, 0.0)[:lnsb]
                except:
                    crot1 = [0.0] * lnsb
                try:
                    pwrot1 = rotavg_ctf(threshold(qa - bckg), ad1, Cs, voltage, pixel_size, 0.0, 0.0)[:lnsb]
                except:
                    pwrot1 = [0.0] * lnsb
                try:
                    crot2 = rotavg_ctf(
                        ctf_rimg(wn, sparx_utilities.generate_ctf([ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1]),
                                 sign=0), ad1, Cs, voltage, pixel_size, bd1, cd1)[:lnsb]
                except:
                    crot2 = [0.0] * lnsb
                try:
                    pwrot2 = rotavg_ctf(threshold(qa - bckg), ad1, Cs, voltage, pixel_size, bd1, cd1)[:lnsb]
                except:
                    pwrot2 = [0.0] * lnsb
                #  #1 - rotational averages without astigmatism, #2 - with astigmatism
                lnsb = min(len(crot2), len(pwrot1), len(crot2), len(pwrot2))
                sparx_utilities.write_text_file(
                    [list(range(lnsb)), [old_div(old_div(float(i), wn), pixel_size)    for i in range(lnsb)], pwrot1, crot1, pwrot2,
                     crot2], os.path.join(outpwrot, "%s_rotinf.txt" % (img_basename_root)))

                #
                # NOTE: 2016/03/23 Toshio Moriya
                # Compute mean of extrema differences (differences at peak & trough) between
                # (1) experimental rotational average with astigmatism (pwrot2)
                # (2) experimental rotational average without astigmatism (pwrot1), and
                # as a indication of goodness of astigmatism estimation by cter.
                # The peak & trough detection uses fitted rotational average with astigmatism (crot2)
                # Start from 1st trough while ignoring 1st peak.
                # End at astigmatism frequency limit.
                #
                is_peak_target = True
                pre_crot2_val = crot2[0]
                extremum_counts = 0
                extremum_diff_sum = 0
                for i in range(1, len(crot2)):
                    cur_crot2_val = crot2[i]
                    if is_peak_target == True and pre_crot2_val > cur_crot2_val:
                        # peak search state
                        extremum_i = i - 1
                        extremum_counts += 1
                        extremum_diff_sum += pwrot2[extremum_i] - pwrot1[
                            extremum_i]  # This should be positive if astigmatism estimation is good
                        # print "MRK_DEBUG: Peak Search  : extremum_i = %03d, freq[extremum_i] = %12.5g, extremum_counts = %03d, (pwrot2[extremum_i] - pwrot1[extremum_i]) = %12.5g, extremum_diff_sum = %12.5g " % (extremum_i, freq[extremum_i] , extremum_counts, (pwrot2[extremum_i] - pwrot1[extremum_i]), extremum_diff_sum)
                        is_peak_target = False
                    elif is_peak_target == False and pre_crot2_val < cur_crot2_val:
                        # trough search state
                        extremum_i = i - 1
                        extremum_counts += 1
                        extremum_diff_sum += pwrot1[extremum_i] - pwrot2[
                            extremum_i]  # This should be positive if astigmatism estimation is good
                        # print "MRK_DEBUG: Trough Search: extremum_i = %03d, freq[extremum_i] = %12.5g, extremum_counts = %03d, (pwrot1[extremum_i] - pwrot2[extremum_i]) = %12.5g, extremum_diff_sum = %12.5g " % (extremum_i, freq[extremum_i] , extremum_counts, (pwrot1[extremum_i] - pwrot2[extremum_i]), extremum_diff_sum)
                        is_peak_target = True
                    pre_crot2_val = cur_crot2_val
                #				#if extremum_counts == 0: ERROR("Logical Error: Encountered unexpected zero extremum counts. Consult with the developer." % (bd1), "%s in %s" % (__name__, os.path.basename(__file__))) # MRK_ASSERT
                #				extremum_diff_avg = 1.1#extremum_diff_sum / extremum_counts

                # print "MRK_DEBUG: extremum_avg = %12.5g, extremum_diff_sum = %12.5g, extremum_counts = %03d," % (extremum_avg, extremum_diff_sum, extremum_counts)

                #				if stack == None:     cmd = "echo " + "    " + namics[ifi] + "  >>  " + fou
                #				else:                 cmd = "echo " + "    " + "  >>  " + fou
                #				os.system(cmd)

                ed1 = wgh  # total amplitude contrast. Here, it is also constant amplitude contrast since Volta phase shift is not estimated with this function.
                stdaved1 = 0.0  # dummy value for error of total amplitude contrast estimation
                max_freq = old_div(0.5 , pixel_size)  # dummy value for maximum frequency. set to Nyquist frequency for now. let's add the implementation in near future (Toshio 2017/12/06)
                reserved = 0.0  # dummy value for reserved spot, which might be used for parameter of external programs (e.g. CTFFIND4, GCTF, and etc.)
                # wgh                     # constant amplitude contrast provided by user (default 10%). Here, it is also total amplitude contrast since Volta phase shift is not estimated with this function.
                phase_shift = ampcont2angle(ed1) - ampcont2angle(
                    wgh)  # Volta phase shift [deg] = total amplitude contrast phase shift [deg] (ed1) -  constant amplitude contrast phase shift [deg]; ed1 is boot strap average of total amplitude contrast [%]

                if debug_mode: print((
                                                 "    %s %s: Process %04d finished the processing. Estimated CTF parmaters are stored in %s." % (
                                         img_type, img_name, ifi, os.path.join(output_directory, "partres.txt"))))
                #				if debug_mode: print((ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, cvavad1, cvavbd1, extremum_diff_avg, ib1, ibec, ctflim))
                #				# totresi.append([img_name, ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, cvavad1, cvavbd1, extremum_diff_avg, ib1, ibec, ctflim])
                #				totresi.append([img_name, ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdaved1, stdavbd1, cd2, cvavad1, cvavbd1, extremum_diff_avg, ib1, ibec, ctflim])
                if debug_mode: print((ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdaved1, stdavbd1,
                                      cd2, cvavad1, cvavbd1, ib1, ibec, ctflim, max_freq, reserved, wgh, phase_shift))
                totresi.append(
                    [img_name, ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1, stdavad1, stdaved1, stdavbd1, cd2,
                     cvavad1, cvavbd1, ib1, ibec, ctflim, max_freq, reserved, wgh, phase_shift])

        #				if stack == None:
        #					print  namics[ifi], ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec
        #				else:
        #					print               ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec
        #				if stack == None:
        #					totresi.append( [ namics[ifi], ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec])
        #				else:
        #					totresi.append( [ 0, ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec])
        #				#if ifi == 4 : break

        if stack == None:
            img_mic = sparx_utilities.get_im(namics[ifi])
            # create  thumbnail
            nx = img_mic.get_xsize()
            if nx > 512:
                img_micthumb = sparx_fundamentals.resample(img_mic, old_div(512.0 , nx))
            else:
                img_micthumb = img_mic
            fou = os.path.join(outmicthumb, "%s_thumb.hdf" % (img_basename_root))
            img_micthumb.write_image(fou)

    if RUNNING_UNDER_MPI:
        pass  # IMPORTIMPORTIMPORT from utilities import wrap_mpi_gatherv
        totresi = sparx_utilities.wrap_mpi_gatherv(totresi, 0, mpi.MPI_COMM_WORLD)
        missing_img_names = sparx_utilities.wrap_mpi_gatherv(missing_img_names, 0, mpi.MPI_COMM_WORLD)
        rejected_img_names = sparx_utilities.wrap_mpi_gatherv(rejected_img_names, 0, mpi.MPI_COMM_WORLD)

    if my_mpi_proc_id == main_mpi_proc:
        outf = open(os.path.join(output_directory, "partres.txt"), "w")
        for i in range(len(totresi)):
            for k in range(1, len(totresi[i])):
                outf.write("  %12.5g" % totresi[i][k])
            outf.write("  %s\n" % totresi[i][0])
        outf.close()

        print(" ")
        print(("Summary of %s processing..." % (img_type.lower())))
        missing_counts = len(missing_img_names)
        print(("  Missing  : %d" % (missing_counts)))
        if missing_counts > 0:
            outfile_path = os.path.join(output_directory, "missing_%s_list.txt" % (img_type.lower()))
            print(("    Saving list of missing in %s..." % (outfile_path)))
            outf = open(outfile_path, "w")
            for missing_img_name in missing_img_names:
                outf.write("%s\n" % missing_img_name)
            outf.close()

        rejected_counts = len(rejected_img_names)
        print(("  Rejected : %d" % (rejected_counts)))
        if rejected_counts > 0:
            outfile_path = os.path.join(output_directory, "rejected_%s_list.txt" % (img_type.lower()))
            print(("    Saving list of rejected in %s..." % (outfile_path)))
            outf = open(outfile_path, "w")
            for rejected_img_name in rejected_img_names:
                outf.write("%s\n" % rejected_img_name)
            outf.close()

    if cter_mode_idx == idx_cter_mode_stack:
        return totresi[0][1], totresi[0][7], totresi[0][8], totresi[0][9], totresi[0][10], totresi[0][11]


################
#
#  CTER code (new version since 2016/03/16)
#
################
#
# NOTE: 2016/03/16 Toshio Moriya
# In this version, the IO-related interface is simplified for sxcter.py and sxgui.py
# Since cter() was used in not only sxcter.py but also e2boxer.py and sxhelixboxer.py,
# This new version is added to avoid breaking e2boxer.py and sxhelixboxer.py
#
# NOTE: 2016/03/16 Toshio Moriya
# To get a single  file name from a GUI application,
# there must be a better way than using guimic...
#
# NOTE: 2016/11/16 Toshio Moriya
# Now, this function assume the MPI setup and clean up is done by caller, such as mpi_init, and mpi_finalize
# 07/11/2017  This is power spectrum version of cter_mrk
def cter_pap(input_image_path, output_directory, selection_list=None, wn=512, \
             pixel_size=-1.0, Cs=2.0, voltage=300.0, wgh=10.0, f_start=-1.0, f_stop=-1.0, \
             kboot=16, overlap_x=50, overlap_y=50, edge_x=0, edge_y=0, \
             check_consistency=False, stack_mode=False, debug_mode=False, \
             program_name="cter_pap() in morphology.py", \
             RUNNING_UNDER_MPI=False, main_mpi_proc=0, my_mpi_proc_id=0, n_mpi_procs=1):
    """Multiline Comment9"""
    pass  # IMPORTIMPORTIMPORT from   EMAN2 import periodogram
    pass  # IMPORTIMPORTIMPORT from   EMAN2db import db_check_dict, db_parse_path
    pass  # IMPORTIMPORTIMPORT from   applications import MPI_start_end
    pass  # IMPORTIMPORTIMPORT from   utilities import read_text_file, write_text_file, get_im, model_blank, model_circle, amoeba, generate_ctf
    pass  # IMPORTIMPORTIMPORT from   utilities import if_error_then_all_processes_exit_program
    pass  # IMPORTIMPORTIMPORT from   utilities import wrap_mpi_bcast
    pass  # IMPORTIMPORTIMPORT from   sys import exit
    pass  # IMPORTIMPORTIMPORT import numpy as np
    pass  # IMPORTIMPORTIMPORT import os
    pass  # IMPORTIMPORTIMPORT import glob
    pass  # IMPORTIMPORTIMPORT from   fundamentals import tilemic, rot_avg_table, resample
    pass  # IMPORTIMPORTIMPORT from   morphology   import threshold, bracket_def, bracket, goldsearch_astigmatism
    pass  # IMPORTIMPORTIMPORT from   morphology   import defocus_baseline_fit, simpw1d, movingaverage, localvariance, defocusgett
    pass  # IMPORTIMPORTIMPORT from   morphology   import defocus_guessn, defocusget_from_crf, make_real
    pass  # IMPORTIMPORTIMPORT from   morphology   import fastigmatism, fastigmatism1, fastigmatism2, fastigmatism3, simctf, simctf2, simctf2out, fupw,ctf2_rimg
    pass  # IMPORTIMPORTIMPORT from   alignment    import Numrinit, ringwe
    pass  # IMPORTIMPORTIMPORT from   statistics   import table_stat
    pass  # IMPORTIMPORTIMPORT from   pixel_error  import angle_ave
    pass  # IMPORTIMPORTIMPORT from   inspect      import currentframe, getframeinfo
    pass  # IMPORTIMPORTIMPORT from   global_def   import ERROR
    pass  # IMPORTIMPORTIMPORT import global_def
    pass  # IMPORTIMPORTIMPORT from   time import time
    pass  # IMPORTIMPORTIMPORT from   mpi import MPI_COMM_WORLD, mpi_barrier

    # ====================================================================================
    # Prepare processing
    # ====================================================================================

    # ------------------------------------------------------------------------------------
    # Find the CTER Running Mode before checking error conditions
    # ------------------------------------------------------------------------------------
    i_enum = -1;
    idx_cter_mode_invalid = i_enum;
    i_enum += 1;
    idx_cter_mode_all_mics = i_enum  # All Micrographs Mode - Process all s in a directory
    i_enum += 1;
    idx_cter_mode_selected_mics = i_enum  # Selected Micrographs Mode - Process all s in a selection list file
    i_enum += 1;
    idx_cter_mode_single_mic = i_enum  # Single Micrograph Mode - Process a single
    i_enum += 1;
    idx_cter_mode_stack = i_enum  # Stack Mode - Process a stack (Advanced Option)
    i_enum += 1;
    idx_cter_mode_counts = i_enum

    cter_mode_idx = idx_cter_mode_invalid
    cter_mode_name = None
    if stack_mode == False:
        # One of three Micrograph Modes
        # For any of Micrograph Modes, input image file name must be a file path pattern containing wild card "*"
        if selection_list == None:
            # User did not use selection list option
            # -> All Micrographs Mode
            cter_mode_idx = idx_cter_mode_all_mics
            cter_mode_name = "All Micrographs Mode"
        else:
            if os.path.splitext(selection_list)[1] == ".txt":
                # User specified a selection list text file path containing".txt" extension through selection list option
                # -> Selected Micrographs Mode
                cter_mode_idx = idx_cter_mode_selected_mics
                cter_mode_name = "Selected Micrographs Mode"
            else:
                # User specified an image file path (a non-text file path) through selection list option
                # -> Single Micrograph Mode
                cter_mode_idx = idx_cter_mode_single_mic
                cter_mode_name = "Single Micrograph Mode"
    else:
        # (Particle) Stack Mode
        cter_mode_idx = idx_cter_mode_stack
        cter_mode_name = "Stack Mode"

    if my_mpi_proc_id == main_mpi_proc:
        print(" ")
        print(("----- Running with %s -----" % (cter_mode_name)))

    # ------------------------------------------------------------------------------------
    # Check mode-dependent error conditions of input arguments and options if abort is necessary. All nodes do this checking
    # ------------------------------------------------------------------------------------
    error_message_list = []  # List of error messages. If no error is found, the length should be zero
    if not stack_mode:

        # Check error conditions applicable to any of Micrograph Mode
        if input_image_path.find("*") == -1:
            error_message_list.append(
                "Input image file path (%s) for %s must be a  path pattern containing wild card (*). Please check input_image_path argument." % (
                input_image_path, cter_mode_name))

        if input_image_path[:len("bdb:")].lower() == "bdb:":
            error_message_list.append(
                "BDB file can not be selected as input image file path (%s) for %s. Please check input_image_path argument and convert the image format." % (
                input_image_path, cter_mode_name))

        # Check error conditions applicable to Selected Micrographs Mode
        if cter_mode_idx == idx_cter_mode_selected_mics:
            if not os.path.exists(selection_list):
                error_message_list.append(
                    "Selection list text file specified by selection_list option (%s) for %s does not exists. Please check selection_list option." % (
                    selection_list, cter_mode_name))

        if cter_mode_idx == idx_cter_mode_single_mic:
            if not os.path.exists(os.path.join(os.path.dirname(input_image_path), os.path.basename(selection_list))):
                error_message_list.append(
                    "Micrograph specified by selection_list option (%s) for %s does not exist. Please check selection_list option." % (
                    selection_list, cter_mode_name))
            #
            if RUNNING_UNDER_MPI and n_mpi_procs != 1:
                error_message_list.append(
                    "%s supports only a single processor version. Please change MPI settings." % (cter_mode_name))

    else:
        # Check error conditions
        if input_image_path.find("*") != -1:
            error_message_list.append(
                "Stack file path specified by input_image_path (%s) for %s should not contain wild card (*). Please check input_image_path argument." % (
                input_image_path, cter_mode_name))

        is_not_found_input_image_file = False
        if input_image_path[:len("bdb:")].lower() == "bdb:":
            if not EMAN2db.db_check_dict(input_image_path):
                is_not_found_input_image_file = True
        else:
            if not os.path.exists(input_image_path):
                is_not_found_input_image_file = True
        if is_not_found_input_image_file:
            error_message_list.append(
                "Stack file specified by input_image_path (%s) for %s does not exist. Please check input_image_path argument." % (
                input_image_path, cter_mode_name))

        if RUNNING_UNDER_MPI and n_mpi_procs != 1:
            error_message_list.append(
                "%s supports only a single processor version. Please change MPI settings." % (cter_mode_name))

    # --------------------------------------------------------------------------------
    # check output-related error conditions (mode-independent). All nodes do this checking
    # --------------------------------------------------------------------------------
    if os.path.exists(output_directory):
        error_message_list.append(
            "Output directory (%s) exists already. Please check output_directory argument." % (output_directory))

    # --------------------------------------------------------------------------------
    # Check error conditions of options (mode-independent). All nodes do this checking
    # --------------------------------------------------------------------------------
    if pixel_size <= 0.0:
        error_message_list.append(
            "Pixel size (%f) must not be negative. Please set a pasitive value larger than 0.0 to pixel_size option." % (
                pixel_size))

    if wn <= 0.0:
        error_message_list.append(
            "CTF window size (%d) must not be negative. Please set a valid value larger than 0 to wn option." % (wn))

    # --------------------------------------------------------------------------------
    # Print all error messages and abort the process if necessary.
    # --------------------------------------------------------------------------------
    error_status = None
    if len(error_message_list) > 0:
        # Detected error! Print all error messages
        if my_mpi_proc_id == main_mpi_proc:
            print(" ")
            for error_message in error_message_list:
                print(("ERROR!!! %s" % (error_message)))
        error_status = ("Detected %d error(s) related to arguments and options. Run %s -h for help. Exiting..." % (
        len(error_message_list), program_name), inspect.getframeinfo(inspect.currentframe()))
    sparx_utilities.if_error_then_all_processes_exit_program(error_status)
    if RUNNING_UNDER_MPI:
        # Wait for all mpi processes to check error conditions, especially existence of output directory
        # Without this barrier, main mpi process can create output directory before some child mpi process check this error.
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    del error_message_list  # Don't need this anymore

    # ------------------------------------------------------------------------------------
    # Check warning conditions of options
    # ------------------------------------------------------------------------------------
    if my_mpi_proc_id == main_mpi_proc:
        if stack_mode:
            if selection_list != None:
                print(" ")
                print(("WARNING!!! --selection_list option will be ignored in %s." % (cter_mode_name)))
            if wn != 512:
                print(" ")
                print(("WARNING!!! --wn option will be ignored in %s." % (cter_mode_name)))
            if overlap_x != 50:
                print(" ")
                print(("WARNING!!! --overlap_x option will be ignored in %s." % (cter_mode_name)))
            if overlap_y != 50:
                print(" ")
                print(("WARNING!!! --overlap_y option will be ignored in %s." % (cter_mode_name)))
            if edge_x != 0:
                print(" ")
                print(("WARNING!!! --edge_x option will be ignored in %s." % (cter_mode_name)))
            if edge_y != 0:
                print(" ")
                print(("WARNING!!! --edge_y option will be ignored in %s." % (cter_mode_name)))
            if check_consistency:
                print(" ")
                print(("WARNING!!! --check_consistency option will be ignored in %s." % (cter_mode_name)))

    # ====================================================================================
    # Create the input file path list and also check input-related error conditions if abort is necessary.
    # ====================================================================================
    input_file_path_list = []
    if not stack_mode:
        # --------------------------------------------------------------------------------
        # Prepare the variables for all sections in  mode case
        # --------------------------------------------------------------------------------
        # Micrograph basename pattern (directory path is removed from  path pattern)
        mic_pattern = input_image_path
        mic_basename_pattern = os.path.basename(mic_pattern)

        # Global entry dictionary (all possible entries from all lists) for all mic id substring
        global_entry_dict = {}  # mic id substring is the key
        subkey_input_mic_path = "Input Micrograph Path"
        subkey_selected_mic_basename = "Selected Micrograph Basename"

        # List keeps only id substrings of s whose all necessary information are available
        valid_mic_id_substr_list = []

        # --------------------------------------------------------------------------------
        # Obtain the list of  id sustrings using a single CPU (i.e. main mpi process)
        # --------------------------------------------------------------------------------
        # NOTE: Toshio Moriya 2016/11/15
        # The below is not a real while.
        # It gives if-statements an opportunity to use break when errors need to be reported
        # However, more elegant way is to use 'raise' statement of exception mechanism...
        #
        error_status = None
        while my_mpi_proc_id == main_mpi_proc:

            # --------------------------------------------------------------------------------
            # Prepare variables for this section
            # --------------------------------------------------------------------------------
            # Prefix and suffix of  basename pattern
            # to find the head/tail indices of  id substring
            mic_basename_tokens = mic_basename_pattern.split('*')
            # Find head index of  id substring
            mic_id_substr_head_idx = len(mic_basename_tokens[0])

            # --------------------------------------------------------------------------------
            # Register  id substrings found in the input directory (specified by  path pattern)
            # to the global entry dictionary
            # --------------------------------------------------------------------------------
            # Generate the list of  paths in the input directory
            print(" ")
            print("Checking the input directory...")
            input_mic_path_list = glob.glob(mic_pattern)
            # Check error condition of input  file path list
            print(("Found %d micrographs in %s." % (len(input_mic_path_list), os.path.dirname(mic_pattern))))
            if len(input_mic_path_list) == 0:
                # The result shouldn't be empty if the specified  file name pattern is invalid
                error_status = (
                "There are no micrographs whose paths match with the specified file path pattern (%s) for %s. Please check input_image_path. Run %s -h for help." % (
                mic_pattern, cter_mode_name, program_name), inspect.getframeinfo(inspect.currentframe()))
                break

            # Register  id substrings to the global entry dictionary
            for input_mic_path in input_mic_path_list:
                # Find tail index of  id substring and extract the substring from the  name
                input_mic_basename = os.path.basename(input_mic_path)
                mic_id_substr_tail_idx = input_mic_basename.index(mic_basename_tokens[1])
                mic_id_substr = input_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
                if not mic_id_substr in global_entry_dict:
                    # print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from input_mic_path_list " % (mic_id_substr))
                    global_entry_dict[mic_id_substr] = {}
                global_entry_dict[mic_id_substr][subkey_input_mic_path] = input_mic_path

            # --------------------------------------------------------------------------------
            # Register  id substrings found in the selection list
            # to the global entry dictionary
            # --------------------------------------------------------------------------------
            # Generate the list of selected  paths in the selection file
            selected_mic_path_list = []
            # Generate  lists according to the execution mode
            if cter_mode_idx == idx_cter_mode_all_mics:
                # Treat all s in the input directory as selected ones
                selected_mic_path_list = input_mic_path_list
            else:
                if os.path.splitext(selection_list)[1] == ".txt":
                    print(" ")
                    print("Checking the selection list...")
                    selected_mic_path_list = sparx_utilities.read_text_file(selection_list)

                    # Check error condition of  entry lists
                    print(("Found %d microgarph entries in %s." % (len(selected_mic_path_list), selection_list)))
                    if len(selected_mic_path_list) == 0:
                        error_status = (
                        "The provided  list file (%s) for %s mode contains no entries. Please check selection_list option and make sure the file contains a  list. Run %s -h for help." % (
                        selection_list, cter_mode_name, program_name), inspect.getframeinfo(inspect.currentframe()))
                        break
                else:
                    print(" ")
                    print(("Processing a single micorgprah: %s..." % (selection_list)))
                    selected_mic_path_list = [selection_list]

                selected_mic_directory = os.path.dirname(selected_mic_path_list[0])
                if selected_mic_directory != "":
                    print(("    NOTE: Program disregards the directory paths in the selection list (%s)." % (
                        selected_mic_directory)))

            # Register  id substrings to the global entry dictionary
            for selected_mic_path in selected_mic_path_list:
                # Find tail index of  id substring and extract the substring from the  name
                selected_mic_basename = os.path.basename(selected_mic_path)
                mic_id_substr_tail_idx = selected_mic_basename.index(mic_basename_tokens[1])
                mic_id_substr = selected_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
                if not mic_id_substr in global_entry_dict:
                    # print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from selected_mic_path_list " % (mic_id_substr))
                    global_entry_dict[mic_id_substr] = {}
                global_entry_dict[mic_id_substr][subkey_selected_mic_basename] = selected_mic_basename

            # --------------------------------------------------------------------------------
            # Clean up variables related to registration to the global entry dictionary
            # --------------------------------------------------------------------------------
            del mic_basename_tokens
            del mic_id_substr_head_idx

            # --------------------------------------------------------------------------------
            # Create the list containing only valid  id substrings
            # --------------------------------------------------------------------------------
            # Prepare lists to keep track of invalid (rejected) s
            no_input_mic_id_substr_list = []

            print(" ")
            print("Checking the input datasets consistency...")

            # Loop over substring id list
            for mic_id_substr in global_entry_dict:
                mic_id_entry = global_entry_dict[mic_id_substr]

                warinnig_messages = []
                # selected  basename must have been registed always .
                if subkey_selected_mic_basename in mic_id_entry:
                    # Check if associated input  exists
                    if not subkey_input_mic_path in mic_id_entry:
                        input_mic_path = mic_pattern.replace("*", mic_id_substr)
                        warinnig_messages.append("    associated input  %s." % (input_mic_path))
                        no_input_mic_id_substr_list.append(mic_id_substr)

                    if len(warinnig_messages) > 0:
                        print(("WARNING!!! Micrograph ID %s does not have:" % (mic_id_substr)))
                        for warinnig_message in warinnig_messages:
                            print(warinnig_message)
                        print("    Ignores this as an invalid entry.")
                    else:
                        # print("MRK_DEBUG: adding mic_id_substr := ", mic_id_substr)
                        valid_mic_id_substr_list.append(mic_id_substr)
            # 	# This entry is not in the selection list. Do nothing

            # Check the input dataset consistency and save the result to a text file, if necessary.
            if check_consistency:
                # Create output directory
                os.mkdir(output_directory)

                # Open the consistency check file
                inconsist_mic_list_path = os.path.join(output_directory, "inconsist_mic_id_file.txt")
                print(" ")
                print(("Generating the input datasets consistency report in %s..." % (inconsist_mic_list_path)))
                inconsist_mic_list_file = open(inconsist_mic_list_path, "w")
                inconsist_mic_list_file.write("# The information about inconsistent  IDs\n")
                # Loop over substring id list
                for mic_id_substr in global_entry_dict:
                    mic_id_entry = global_entry_dict[mic_id_substr]

                    consistency_messages = []
                    # Check if associated input  path exists
                    if not subkey_input_mic_path in mic_id_entry:
                        input_mic_path = mic_pattern.replace("*", mic_id_substr)
                        consistency_messages.append("    associated input  %s." % (input_mic_path))

                    # Check if associated selected  basename exists
                    if not subkey_selected_mic_basename in mic_id_entry:
                        input_mic_path = mic_pattern.replace("*", mic_id_substr)
                        consistency_messages.append("    associated selected  %s." % (input_mic_path))

                    if len(consistency_messages) > 0:
                        inconsist_mic_list_file.write("Micrograph ID %s does not have:\n" % (mic_id_substr))
                        for consistency_message in consistency_messages:
                            inconsist_mic_list_file.write(consistency_message)
                            inconsist_mic_list_file.write("\n")

                # Close the consistency check file, if necessary
                inconsist_mic_list_file.flush()
                inconsist_mic_list_file.close()

            # Since mic_id_substr is once stored as the key of global_entry_dict and extracted with the key order
            # we need sort the valid_mic_id_substr_list here
            if debug_mode: print(("BEFORE SORT: valid_mic_id_substr_list := ", valid_mic_id_substr_list))
            valid_mic_id_substr_list.sort(key=str.lower)  # Sort list of  IDs using case insensitive string comparison
            if debug_mode: print(("AFTER SORT: valid_mic_id_substr_list := ", valid_mic_id_substr_list))

            # --------------------------------------------------------------------------------
            # Print out the summary of input consistency
            # --------------------------------------------------------------------------------
            print(" ")
            print("Summary of dataset consistency check...")
            print(("  Detected  IDs               : %6d" % (len(global_entry_dict))))
            print(("  Entries in input directory  : %6d" % (len(input_mic_path_list))))
            print(("  Entries in selection list   : %6d" % (len(selected_mic_path_list))))
            print(("  Rejected by no input        : %6d" % (len(no_input_mic_id_substr_list))))
            print(("  Valid Entries               : %6d" % (len(valid_mic_id_substr_list))))

            # --------------------------------------------------------------------------------
            # Check MPI error condition
            # --------------------------------------------------------------------------------
            if len(valid_mic_id_substr_list) < n_mpi_procs:
                error_status = (
                "Number of MPI processes (%d) supplied by --np in mpirun cannot be greater than %d (number of valid s that satisfy all criteria to be processed). Run %s -h for help." % (
                n_mpi_procs, len(valid_mic_id_substr_list, program_name)), inspect.getframeinfo(inspect.currentframe()))
                break

            # --------------------------------------------------------------------------------
            # Create input file path list
            # --------------------------------------------------------------------------------
            for mic_id_substr in valid_mic_id_substr_list:
                mic_path = global_entry_dict[mic_id_substr][subkey_input_mic_path]
                input_file_path_list.append(mic_path)

            # --------------------------------------------------------------------------------
            # Clean up variables related to tracking of invalid (rejected) s
            # --------------------------------------------------------------------------------
            del input_mic_path_list
            del selected_mic_path_list
            del no_input_mic_id_substr_list

            break

        # --------------------------------------------------------------------------------
        # Clean up the variables for all sections in  mode case
        # --------------------------------------------------------------------------------
        del mic_pattern
        del mic_basename_pattern
        del global_entry_dict
        del subkey_input_mic_path
        del subkey_selected_mic_basename
        del valid_mic_id_substr_list

        # --------------------------------------------------------------------------------
        # Print all error messages and abort the process if necessary.
        # --------------------------------------------------------------------------------
        # NOTE: Toshio Moriya 2016/11/15
        # The following function takes care of the case when an if-statement uses break for occurence of an error.
        # However, more elegant way is to use 'exception' statement of exception mechanism...
        #
        sparx_utilities.if_error_then_all_processes_exit_program(error_status)

    else:
        input_file_path_list.append(input_image_path)

    if RUNNING_UNDER_MPI:
        # Wait for main mpi process to create the input file path list
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

        # All mpi processes should know input file path list
        input_file_path_list = sparx_utilities.wrap_mpi_bcast(input_file_path_list, main_mpi_proc)

    # ====================================================================================
    # Prepare input file path(s)
    # ====================================================================================
    #
    # NOTE: 2016/03/17 Toshio Moriya
    # From here on, stack (and namics) will be used to distinguish stack mode and  mode.
    # However, a single input_file_path_list should be sufficient since we already know the mode.
    # Let's consider this refactoring in the future.
    #
    stack = None  # (particle) stack file name: if it is not None, cter runs with stack mode. Otherwise, runs with  mode
    namics = []  # file name list
    if not stack_mode:
        namics = input_file_path_list
        if debug_mode: print(("BEFORE SORT: namics := ", namics))
        namics.sort(key=str.lower)  # Sort list of s using case insensitive string comparison
        if debug_mode: print(("AFTER SORT: namics := ", namics))
    else:
        stack = input_file_path_list[0]

    del input_file_path_list  # Don't need this anymore

    # Make output directory
    outpwrot = "%s/pwrot" % (output_directory)
    if stack == None:
        outmicthumb = "%s/micthumb" % (output_directory)
    if debug_mode:
        outravg = "%s/ravg" % (output_directory)
    if my_mpi_proc_id == main_mpi_proc:
        # Make output directory
        if not os.path.exists(output_directory):
            os.mkdir(output_directory)
        os.mkdir(outpwrot)
        if stack == None:
            os.mkdir(outmicthumb)
        if debug_mode:
            os.mkdir(outravg)

    if RUNNING_UNDER_MPI:
        # Make all mpi processes wait for main mpi process to create output directory
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    # Set up loop variables depending on the cter mode
    if stack == None:
        if RUNNING_UNDER_MPI:
            set_start, set_end = sparx_applications.MPI_start_end(len(namics), n_mpi_procs, my_mpi_proc_id)
        else:
            set_start = 0
            set_end = len(namics)
    else:
        pw2 = []
        data = EMAN2_cppwrap.EMData.read_images(stack)
        nima = len(data)
        for i in range(nima):
            pw2.append(EMAN2_cppwrap.periodogram(data[i]))
        wn = pw2[0].get_xsize()
        set_start = 0
        set_end = 1

    # Set up progress message
    if my_mpi_proc_id == main_mpi_proc:
        print(" ")
        print("Estimating CTF parameters...")
        if stack == None:
            print("  Micrographs processed by main process (including percent of progress):")
            progress_percent_step = old_div((
                                                set_end - set_start) , 100.0)  # the number of micrograms for main mpi processer divided by 100

    totresi = []
    missing_img_names = []
    rejected_img_names = []
    for ifi in range(set_start, set_end):
        img_type = ""
        img_name = ""
        img_basename_root = ""

        if stack == None:
            img_type = "Micrograph"
            img_name = namics[ifi]

            if my_mpi_proc_id == main_mpi_proc:
                print(("    Processing %s ---> %6.2f%%" % (img_name, old_div((ifi - set_start), progress_percent_step))))

            if not os.path.exists(img_name):
                missing_img_names.append(img_name)
                print(
                    "    %s %s: Can not find this file. Skipping the estimation and CTF parameters are not stored..." % (
                    img_type, img_name))
                continue

            mic = sparx_utilities.get_im(img_name)
            try:
                pw2 = sparx_fundamentals.tilemic(mic, win_size=wn, overlp_x=overlap_x, overlp_y=overlap_y,
                                                 edge_x=edge_x, edge_y=edge_y)
            except:
                print(
                    "MRK_DEBUG: tilemic() in cter_pap() raised an exception. The micrographs {} might have a problem. Please check it and remove it if necessary.".format(
                        img_name))
                raise
            del mic
        else:
            img_type = "Stack"
            img_name = stack

            numFM = EMAN2_cppwrap.EMUtil.get_image_count(img_name)
            pw2 = [None] * numFM
            for i in range(numFM):
                pw2.append(EMAN2_cppwrap.periodogram(sparx_utilities.get_im(img_name, i)))

        if debug_mode: print("    %s %s: Process %04d started the processing. Detected %d image(s) in this %s file." % (
        img_type, img_name, ifi, numFM, img_type.lower()))

        if EMAN2db.db_check_dict(img_name) == False:
            img_basename_root = os.path.splitext(os.path.basename(img_name))[0]
        else:
            path, dictname, keys = EMAN2db.db_parse_path(img_name)
            img_basename_root = dictname

        nimi = len(pw2)
        adefocus = [0.0] * kboot
        aamplitu = [0.0] * kboot
        aangle = [0.0] * kboot

        allroo = []
        for imi in range(nimi):
            allroo.append(sparx_fundamentals.rot_avg_table(pw2[imi]))
        lenroo = len(allroo[0])
        # print time(),nimi

        for nboot in range(kboot):
            # at = time()
            if (nboot == 0):
                boot = list(range(nimi))
            else:
                pass  # IMPORTIMPORTIMPORT from random import randint
                for imi in range(nimi): boot[imi] = random.randint(0, nimi - 1)
            qa = sparx_utilities.model_blank(wn, wn)
            roo = np.zeros(lenroo, np.float32)
            sroo = np.zeros(lenroo, np.float32)
            aroo = np.zeros(lenroo, np.float32)

            for imi in range(nimi):
                EMAN2_cppwrap.Util.add_img(qa, pw2[boot[imi]])
                temp1 = np.array(allroo[boot[imi]])
                roo += temp1
                temp2 = movingaverage(temp1, 10)
                aroo += temp2
                sroo += temp2 ** 2
            sroo[0] = sroo[1]
            aroo[0] = aroo[1]
            sroo = old_div((sroo - old_div(aroo ** 2 , nimi)) , nimi)
            aroo = old_div(aroo,nimi)
            roo = old_div(roo,nimi)
            qa = old_div(qa, nimi)

            if f_start < 0:
                #  Find a break point
                bp = 1.e23
                for i in range(5, lenroo - 5):
                    # t1 = linreg(sroo[:i])
                    # t2 = linreg(sroo[i:])
                    # tt = t1[1][0] + t2[1][0]
                    xtt = np.array(list(range(i)), np.float32)
                    zet = np.poly1d(np.polyfit(xtt, sroo[:i], 2))
                    t1 = sum((sroo[:i] - zet(xtt)) ** 2)
                    xtt = np.array(list(range(i, lenroo)), np.float32)
                    zet = np.poly1d(np.polyfit(xtt, sroo[i:], 2))
                    tt = t1 + sum((sroo[i:] - zet(xtt)) ** 2)
                    if tt < bp:
                        bp = tt
                        istart = i
                # istart = 25
                # print istart
                f_start = old_div(istart , (pixel_size * wn))
            """Multiline Comment10"""
            # write_text_file([roo.tolist(),aroo.tolist(),sroo.tolist()], "sroo%03d.txt"%ifi)
            rooc = roo.tolist()

            # print namics[ifi],istart,f_start

            defc, subpw, ctf2, baseline, envelope, istart, istop = defocusgett(rooc, wn, voltage=voltage,
                                                                               Pixel_size=pixel_size, Cs=Cs,
                                                                               ampcont=wgh, f_start=f_start,
                                                                               f_stop=f_stop, round_off=1.0, nr1=3,
                                                                               nr2=6, parent=None, DEBug=debug_mode)
            if debug_mode:
                print("  RESULT %s" % (img_name), defc, istart, istop)

                freq = list(range(len(subpw)))
                for i in range(len(freq)):  freq[i] = old_div(old_div(float(i) , wn) , pixel_size)
                #				write_text_file([freq, subpw.tolist(), ctf2, envelope.tolist(), baseline.tolist()], "%s/ravg%05d.txt" % (output_directory, ifi))
                fou = os.path.join(outravg, "%s_ravg_%02d.txt" % (img_basename_root, nboot))
                sparx_utilities.write_text_file([freq, subpw.tolist(), ctf2, envelope.tolist(), baseline.tolist()], fou)
            # mpi_barrier(MPI_COMM_WORLD)

            # exit()
            bg = baseline.tolist()
            en = envelope.tolist()

            bckg = sparx_utilities.model_blank(wn, wn, 1, 1)
            envl = sparx_utilities.model_blank(wn, wn, 1, 1)

            pass  # IMPORTIMPORTIMPORT from math import sqrt
            nc = old_div(wn , 2)
            bg.append(bg[-1])
            en.append(en[-1])
            for i in range(wn):
                for j in range(wn):
                    r = numpy.sqrt((i - nc) ** 2 + (j - nc) ** 2)
                    ir = int(r)
                    if (ir < nc):
                        dr = r - ir
                        bckg.set_value_at(i, j, (1. - dr) * bg[ir] + dr * bg[ir + 1])
                        envl.set_value_at(i, j, (1. - dr) * en[ir] + dr * en[ir + 1])

            # qa.write_image("rs1.hdf")

            mask = sparx_utilities.model_circle(istop - 1, wn, wn) * (
                        sparx_utilities.model_blank(wn, wn, 1, 1.0) - sparx_utilities.model_circle(istart, wn, wn))
            qse = threshold((qa - bckg))  # *envl
            # (qse*mask).write_image("rs2.hdf")
            # qse.write_image("rs3.hdf")
            ##  SIMULATION
            # bang = 0.7
            # qse = ctf2_rimg(wn, generate_ctf([defc,Cs,voltage,pixel_size,0.0,wgh, bang, 37.0]) )
            # qse.write_image("rs3.hdf")

            cnx = old_div(wn, 2) + 1
            cny = cnx
            mode = "H"
            istop = min( old_div(wn , 2) - 2, istop)  # 2-26-2015@ming
            numr = sparx_alignment.Numrinit(istart, istop, 1, mode)
            wr = sparx_alignment.ringwe(numr, mode)

            crefim = EMAN2_cppwrap.Util.Polar2Dm(qse * mask, cnx, cny, numr, mode)
            EMAN2_cppwrap.Util.Frngs(crefim, numr)
            EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)

            # pc = ctf2_rimg(wn,generate_ctf([defc,Cs,voltage,pixel_size,0.0,wgh]))
            # print ccc(pc*envl, subpw, mask)

            bang = 0.0
            bamp = 0.0
            bdef = defc
            bold = 1.e23
            while (True):
                #  in simctf2 data[3] is astigmatism amplitude
                """
				data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
				#astdata = [crefim, numr, wn, bdef, Cs, voltage, pixel_size, wgh, bang]
				for qqq in xrange(200):
					qbdef = 1.0 + qqq*0.001
					print " VALUE AT THE BEGGINING OF while LOOP  ",qbdef,simctf2(qbdef, data)#,fastigmatism3(bamp,astdata)
				"""
                """
				bamp = 0.7
				bang = 37.0

				data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
				astdata = [crefim, numr, wn, bdef, Cs, voltage, Pixel_size, wgh, bang]
				print " VALUE AT THE BEGGINING OF while LOOP  ",bdef,bamp,bang,simctf2(bdef, data),fastigmatism3(bamp,astdata,mask)
				#print  simctf2out(1.568,data)
				#exit()

				for kdef in xrange(14000,17000,10):
					dz = kdef/10000.0
					ard = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
					#print ard
					aqd = [crefim, numr, wn, dz, Cs, voltage, pixel_size, wgh, bang]
					#print aqd
					print  dz,simctf2(dz,ard),fastigmatism3(bamp,aqd,mask)
					#print aqd[-1]
				exit()
				"""
                data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
                h = 0.05 * bdef
                amp1, amp2 = bracket_def(simctf2, data, bdef * 0.9, h)
                # print "bracketing of the defocus  ",amp1, amp2
                # print " ttt ",time()-srtt
                # print "bracketing of the defocus  ",amp1,amp2,simctf2(amp1, data),simctf2(amp2, data),h
                amp1, val2 = goldsearch_astigmatism(simctf2, data, amp1, amp2, tol=1.0e-3)
                # print "golden defocus ",amp1, val2,simctf2(amp1, data)
                # bdef, bcc = goldsearch_astigmatism(simctf2, data, amp1, amp2, tol=1.0e-3)
                # print "correction of the defocus  ",bdef,bcc
                # print " ttt ",time()-srtt
                """Multiline Comment11"""

                astdata = [crefim, numr, wn, bdef, Cs, voltage, pixel_size, wgh, bang, mask]
                h = 0.01
                amp1, amp2 = bracket(fastigmatism3, astdata, h)
                # print "  astigmatism bracket  ",amp1,amp2,astdata[-1]
                # print " ttt ",time()-srtt
                bamp, bcc = goldsearch_astigmatism(fastigmatism3, astdata, amp1, amp2, tol=1.0e-3)
                junk = fastigmatism3(bamp, astdata)
                bang = astdata[8]

                # print astdata[8]
                # print  fastigmatism3(0.0,astdata)
                # print astdata[8]
                # temp = 0.0
                # print bdef, Cs, voltage, pixel_size, temp, wgh, bamp, bang, -bcc
                # data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
                # astdata = [crefim, numr, wn, bdef, Cs, voltage, pixel_size, wgh, bang]
                # print " VALUE WITHIN the while LOOP  ",bdef,bamp,bang,simctf2(bdef, data),fastigmatism3(bamp,astdata)
                # print "  golden search ",bamp,data[-1], fastigmatism3(bamp,data), fastigmatism3(0.0,data)
                # print " ttt ",time()-srtt
                # bamp = 0.5
                # bang = 277

                dama = sparx_utilities.amoeba([bdef, bamp], [0.2, 0.2], fupw, 1.e-4, 1.e-4, 500, astdata)
                if debug_mode:  print("AMOEBA    ", dama)
                bdef = dama[0][0]
                bamp = dama[0][1]
                astdata = [crefim, numr, wn, bdef, Cs, voltage, pixel_size, wgh, bang, mask]
                junk = fastigmatism3(bamp, astdata)
                bang = astdata[8]
                if debug_mode:  print(" after amoeba ", bdef, bamp, bang)
                #  The looping here is blocked as one shot at amoeba is good enough.  To unlock it, remove - from bold.
                if (bcc < -bold):
                    bold = bcc
                else:
                    break

            # data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
            # print " VALUE AFTER the while LOOP  ",bdef,bamp,bang,simctf2(bdef, data),fastigmatism3(bamp,astdata)
            # temp = 0.0
            # print ifi,bdef, Cs, voltage, pixel_size, temp, wgh, bamp, bang, -bcc
            # freq = range(len(subpw))
            # for i in xrange(len(freq)):  freq[i] = float(i)/wn/pixel_size
            # ctf2 = ctf_2(wn, generate_ctf([bdef,Cs,voltage,pixel_size,0.0,wgh]))[:len(freq)]
            # write_text_file([freq, subpw.tolist(), ctf2, envelope.tolist(), baseline.tolist()],"ravg/ravg%05d.txt"%ifi)
            # print " >>>> ",wn, bdef, bamp, Cs, voltage, pixel_size, wgh, bang
            # data = [qse, mask, wn, bamp, Cs, voltage, pixel_size, wgh, bang]
            # print  simctf2out(bdef, data)
            # exit()
            adefocus[nboot] = bdef
            aamplitu[nboot] = bamp
            aangle[nboot] = bang
        # from sys import exit
        # exit()

        # print " ttt ",time()-srtt
        # from sys import exit
        # exit()
        ad1, ad2, ad3, ad4 = sparx_statistics.table_stat(adefocus)  # return values: average, variance, minimum, maximum
        if ad2 <= 0.0:
            print((
                              "    %s %s: Detected the variance less than zero (defocus statics: avg = %f, var = %f, min = %f, max = %f)." % (
                      img_type, img_name, ad1, ad2, ad3, ad4)))
            print(("           The program ignores this estimate..."))
            continue

        reject = []
        thr = 3 * numpy.sqrt(ad2)
        for i in range(len(adefocus)):
            if (abs(adefocus[i] - ad1) > thr):
                print((
                                  "    %s %s: Rejected an outlier defocus estimate (defocus = %f, average defocus = %f, threshold = %f)." % (
                          img_type, img_name, adefocus[i], ad1, thr)))
                reject.append(i)

        if (len(reject) > 0):
            print(("    %s %s: Total number of rejects %s" % (img_type, img_name, len(reject))))
            for i in range(len(reject) - 1, -1, -1):
                del adefocus[i]
                del aamplitu[i]
                del aangle[i]

        if (len(adefocus) < 2):
            print((
                              "    %s %s: After rejection of outliers, there is too few estimated defocus values. Skipping the estimation and CTF parameters are not stored..." % (
                      img_type, img_name)))
        else:
            # print "adefocus",adefocus
            # print  "aamplitu",aamplitu
            # print "aangle",aangle
            ad1, ad2, ad3, ad4 = sparx_statistics.table_stat(adefocus)
            bd1, bd2, bd3, bd4 = sparx_statistics.table_stat(aamplitu)
            cd1, cd2 = sparx_pixel_error.angle_ave(
                [2 * q for q in aangle])  # Have to use this trick as the function works for range [0,360]
            cd1 = old_div(cd1,2)
            cd2 = old_div(cd2,2)
            temp = 0.0
            stdavad1 = np.sqrt(kboot * max(0.0, ad2))
            stdavbd1 = np.sqrt(kboot * max(0.0, bd2))
            cd2 *= np.sqrt(kboot)

            # Adjust value ranges of astig. amp. and angle.
            if (bd1 < 0.0):
                bd1 = -bd1
                cd1 = 90.0 + cd1
            cd1 = cd1 % 180

            #  SANITY CHECK, do not produce anything if defocus abd astigmatism amplitude are out of whack
            reject_img_messages = []
            try:
                pwrot2 = rotavg_ctf(sparx_utilities.model_blank(wn, wn), ad1, Cs, voltage, pixel_size, bd1, cd1)
            except:
                reject_img_messages.append(
                    "    - Astigmatism amplitude (%f) is larger than defocus (%f) or defocus (%f) is negative." % (
                    bd1, ad1, ad1))

            if len(reject_img_messages) > 0:
                rejected_img_names.append(img_name)
                print("    %s %s: Rejected the CTF estimate - " % (img_type, img_name), ad1, Cs, voltage, pixel_size,
                      wgh, bd1, cd1, "(def, Cs, vol, apix, amp_contrast, astig_amp, astig_angle)")
                print("    %s %s: because... " % (img_type, img_name))
                for reject_img_message in reject_img_messages:
                    print(reject_img_message)
                print("    %s %s: Skipping the estimation and CTF parameters are not stored..." % (img_type, img_name))
            else:  # assert(len(img_reject_messages) == 0)
                #  Estimate the point at which (sum_errordz ctf_1(dz+errordz))^2 falls to 0.5
                import random as rqt

                supe = sparx_utilities.model_blank(wn, wn)
                niter = 1000
                for it in range(niter):
                    EMAN2_cppwrap.Util.add_img(supe,
                                               EMAN2_cppwrap.Util.ctf_rimg(wn, wn, 1, ad1 + rqt.gauss(0.0, stdavad1),
                                                                           pixel_size, voltage, Cs, 0.0, wgh,
                                                                           bd1 + rqt.gauss(0.0, stdavbd1),
                                                                           cd1 + rqt.gauss(0.0, cd2), 1))
                ni = old_div(wn, 2)
                supe = old_div(supe,niter)
                pwrot2 = rotavg_ctf(supe, ad1, Cs, voltage, pixel_size, bd1, cd1)
                for i in range(ni):  pwrot2[i] = pwrot2[i] ** 2

                ibec = 0
                for it in range(ni - 1, 0, -1):
                    if pwrot2[it] > 0.5:
                        ibec = it
                        break
                pass  # IMPORTIMPORTIMPORT from morphology import ctf_1d
                ct = sparx_utilities.generate_ctf([ad1, Cs, voltage, pixel_size, temp, wgh, 0.0, 0.0])
                cq = ctf_1d(wn, ct)

                supe = [0.0] * ni
                niter = 1000
                for i in range(niter):
                    cq = sparx_utilities.generate_ctf(
                        [ad1 + rqt.gauss(0.0, stdavad1), Cs, voltage, pixel_size, 0.0, wgh, 0.0, 0.0])
                    ci = ctf_1d(wn, cq)[:ni]
                    for l in range(ni):  supe[l] += ci[l]

                for l in range(ni):  supe[l] = (old_div(supe[l] , niter)) ** 2

                ib1 = 0
                for it in range(ni - 1, 0, -1):
                    if supe[it] > 0.5:
                        ib1 = it
                        break
                ibec = old_div(ibec , (pixel_size * wn))  # with astigmatism
                ib1 = old_div(ib1 , (pixel_size * wn))  # no astigmatism
                # from utilities import write_text_file
                # write_text_file([range(ni), supe[:ni],pwrot2[:ni]],"fifi.txt")

                # Compute defocus CV and astig. amp. CV (CV: coefficient of variation; ratio of error (SD) relative to average (mean))
                # if ad1 < max(0.0, valid_min_defocus): ERROR("Logical Error: Encountered unexpected defocus value (%f). Consult with the developer." % (ad1), "%s in %s" % (__name__, os.path.basename(__file__))) # MRK_ASSERT
                #  TOSHIO - no need to ceck, it was computed as sqrt above, so it cannot be <0
                # if stdavad1 < 0.0: ERROR("Logical Error: Encountered unexpected defocus SD value (%f). Consult with the developer." % (stdavad1), "%s in %s" % (__name__, os.path.basename(__file__))) # MRK_ASSERT
                cvavad1 = old_div(stdavad1, ad1) * 100  # use percentage

                if bd1 < 0.0: sparx_global_def.ERROR(
                    "Logical Error: Encountered unexpected astig. amp. value (%f). Consult with the developer." % (bd1),
                    "%s in %s" % (__name__, os.path.basename(__file__)))  # MRK_ASSERT

                bd1 = max(bd1, 1.0e-15)

                cvavbd1 = old_div(stdavbd1, bd1) * 100  # use percentage

                # Compute CTF limit (theoretical resolution limit based on the oscillations of CTF)
                # For output, use ctflim (relative frequency limit [1/A]), not ctflim_abs (absolute frequency limit)
                #
                # NOTE: 2016/03/23 Toshio Moriya
                # xr is limiting frequency [1/A]. Max is Nyquist frequency = 1.0/(2*apix[A/pixel]). <UNIT: [1/(A/pixel)/[pixel])] => [(pixel)/(A*pixel] => [1/A]>
                # 1.0/xr is limiting period (Angstrom resolution) [A]. Min is Nyquist period = (2*apix[A/pixel]). <UNIT: [1/(1/A)] = [A]>
                # fwpix is width of Fourier pixel [pixel/A] := 1.0[pixel]/(2*apix[A/pixel])/box_half[pixel] = 1[pixel]/fullsize[A]). <UNIT: [pixel/(A/pixel)/(pixel)] = [pixel*(pixel/A)*(1/pixel) = [pixel/A]>
                # int(xr/fwpix+0.5) is limiting_absolute_frequency [1/pixel]. <Unit:[(1/A)/(pixel/A)] = [(1/A)*(A/pixel)] = [1/pixel]>
                # return  int(xr/fwpix+0.5),xr, which is limiting_abs_frequency [1/pixel], and Limiting_frequency[1/A]
                #
                ctflim_abs, ctflim = ctflimit(wn, ad1, Cs, voltage, pixel_size)

                """Multiline Comment12"""
                lnsb = len(subpw)
                try:
                    crot1 = rotavg_ctf(
                        ctf_rimg(wn, sparx_utilities.generate_ctf([ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1]),
                                 sign=0), ad1, Cs, voltage, pixel_size, 0.0, 0.0)[:lnsb]
                except:
                    crot1 = [0.0] * lnsb
                try:
                    pwrot1 = rotavg_ctf(threshold(qa - bckg), ad1, Cs, voltage, pixel_size, 0.0, 0.0)[:lnsb]
                except:
                    pwrot1 = [0.0] * lnsb
                try:
                    crot2 = rotavg_ctf(
                        ctf_rimg(wn, sparx_utilities.generate_ctf([ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1]),
                                 sign=0), ad1, Cs, voltage, pixel_size, bd1, cd1)[:lnsb]
                except:
                    crot2 = [0.0] * lnsb
                try:
                    pwrot2 = rotavg_ctf(threshold(qa - bckg), ad1, Cs, voltage, pixel_size, bd1, cd1)[:lnsb]
                except:
                    pwrot2 = [0.0] * lnsb
                #  #1 - rotational averages without astigmatism, #2 - with astigmatism
                lnsb = min(lnsb, len(crot2), len(pwrot1), len(crot2), len(pwrot2))
                sparx_utilities.write_text_file(
                    [list(range(lnsb)), [ old_div(old_div(float(i) , wn) , pixel_size) for i in range(lnsb)], pwrot1[:lnsb], crot1[:lnsb],
                     pwrot2[:lnsb], crot2[:lnsb]], os.path.join(outpwrot, "%s_rotinf.txt" % (img_basename_root)))
                #
                # NOTE: 2016/03/23 Toshio Moriya
                # Compute mean of extrema differences (differences at peak & trough) between
                # (1) experimental rotational average with astigmatism (pwrot2)
                # (2) experimental rotational average without astigmatism (pwrot1), and
                # as a indication of goodness of astigmatism estimation by cter.
                # The peak & trough detection uses fitted rotational average with astigmatism (crot2)
                # Start from 1st trough while ignoring 1st peak.
                # End at astigmatism frequency limit.
                #
                is_peak_target = True
                pre_crot2_val = crot2[0]
                extremum_counts = 0
                extremum_diff_sum = 0
                for i in range(1, len(crot2)):
                    cur_crot2_val = crot2[i]
                    if is_peak_target == True and pre_crot2_val > cur_crot2_val:
                        # peak search state
                        extremum_i = i - 1
                        extremum_counts += 1
                        extremum_diff_sum += pwrot2[extremum_i] - pwrot1[
                            extremum_i]  # This should be positive if astigmatism estimation is good
                        # print "MRK_DEBUG: Peak Search  : extremum_i = %03d, freq[extremum_i] = %12.5g, extremum_counts = %03d, (pwrot2[extremum_i] - pwrot1[extremum_i]) = %12.5g, extremum_diff_sum = %12.5g " % (extremum_i, freq[extremum_i] , extremum_counts, (pwrot2[extremum_i] - pwrot1[extremum_i]), extremum_diff_sum)
                        is_peak_target = False
                    elif is_peak_target == False and pre_crot2_val < cur_crot2_val:
                        # trough search state
                        extremum_i = i - 1
                        extremum_counts += 1
                        extremum_diff_sum += pwrot1[extremum_i] - pwrot2[
                            extremum_i]  # This should be positive if astigmatism estimation is good
                        # print "MRK_DEBUG: Trough Search: extremum_i = %03d, freq[extremum_i] = %12.5g, extremum_counts = %03d, (pwrot1[extremum_i] - pwrot2[extremum_i]) = %12.5g, extremum_diff_sum = %12.5g " % (extremum_i, freq[extremum_i] , extremum_counts, (pwrot1[extremum_i] - pwrot2[extremum_i]), extremum_diff_sum)
                        is_peak_target = True
                    pre_crot2_val = cur_crot2_val
                #				#if extremum_counts == 0: ERROR("Logical Error: Encountered unexpected zero extremum counts. Consult with the developer." % (bd1), "%s in %s" % (__name__, os.path.basename(__file__))) # MRK_ASSERT
                #				extremum_diff_avg = 1.1#extremum_diff_sum / extremum_counts

                # print "MRK_DEBUG: extremum_avg = %12.5g, extremum_diff_sum = %12.5g, extremum_counts = %03d," % (extremum_avg, extremum_diff_sum, extremum_counts)

                #				if stack == None:     cmd = "echo " + "    " + namics[ifi] + "  >>  " + fou
                #				else:                 cmd = "echo " + "    " + "  >>  " + fou
                #				os.system(cmd)

                ed1 = wgh  # total amplitude contrast. Here, it is also constant amplitude contrast since Volta phase shift is not estimated with this function.
                stdaved1 = 0.0  # dummy value for error of total amplitude contrast estimation
                max_freq = old_div(0.5 , pixel_size)  # dummy value for maximum frequency. set to Nyquist frequency for now. let's add the implementation in near future (Toshio 2017/12/06)
                reserved = 0.0  # dummy value for reserved spot, which might be used for parameter of external programs (e.g. CTFFIND4, GCTF, and etc.)
                # wgh                     # constant amplitude contrast provided by user (default 10%). Here, it is also total amplitude contrast since Volta phase shift is not estimated with this function.
                phase_shift = ampcont2angle(ed1) - ampcont2angle(
                    wgh)  # Volta phase shift [deg] = total amplitude contrast phase shift [deg] (ed1) -  constant amplitude contrast phase shift [deg]; ed1 is boot strap average of total amplitude contrast [%]

                if debug_mode: print((
                                                 "    %s %s: Process %04d finished the processing. Estimated CTF parmaters are stored in %s." % (
                                         img_type, img_name, ifi, os.path.join(output_directory, "partres.txt"))))
                #				if debug_mode: print((ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, cvavad1, cvavbd1, extremum_diff_avg, ib1, ibec, ctflim))
                #				# totresi.append([img_name, ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, cvavad1, cvavbd1, extremum_diff_avg, ib1, ibec, ctflim])
                #				totresi.append([img_name, ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdaved1, stdavbd1, cd2, cvavad1, cvavbd1, extremum_diff_avg, ib1, ibec, ctflim])
                if debug_mode: print((ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1, stdavad1, stdavbd1, stdavbd1,
                                      cd2, cvavad1, cvavbd1, ib1, ibec, ctflim, max_freq, reserved, wgh, phase_shift))
                totresi.append(
                    [img_name, ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1, stdavad1, stdaved1, stdavbd1, cd2,
                     cvavad1, cvavbd1, ib1, ibec, ctflim, max_freq, reserved, wgh, phase_shift])

        #				if stack == None:
        #					print  namics[ifi], ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec
        #				else:
        #					print               ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec
        #				if stack == None:
        #					totresi.append( [ namics[ifi], ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec])
        #				else:
        #					totresi.append( [ 0, ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec])
        #				#if ifi == 4 : break

        if stack == None:
            img_mic = sparx_utilities.get_im(namics[ifi])
            # create  thumbnail
            nx = img_mic.get_xsize()
            if nx > 512:
                img_micthumb = sparx_fundamentals.resample(img_mic,old_div(512.0 , nx))
            else:
                img_micthumb = img_mic
            fou = os.path.join(outmicthumb, "%s_thumb.hdf" % (img_basename_root))
            img_micthumb.write_image(fou)

    if RUNNING_UNDER_MPI:
        pass  # IMPORTIMPORTIMPORT from utilities import wrap_mpi_gatherv
        totresi = sparx_utilities.wrap_mpi_gatherv(totresi, 0, mpi.MPI_COMM_WORLD)
        missing_img_names = sparx_utilities.wrap_mpi_gatherv(missing_img_names, 0, mpi.MPI_COMM_WORLD)
        rejected_img_names = sparx_utilities.wrap_mpi_gatherv(rejected_img_names, 0, mpi.MPI_COMM_WORLD)

    if my_mpi_proc_id == main_mpi_proc:
        outf = open(os.path.join(output_directory, "partres.txt"), "w")
        for i in range(len(totresi)):
            for k in range(1, len(totresi[i])):
                outf.write("  %12.5g" % totresi[i][k])
            outf.write("  %s\n" % totresi[i][0])
        outf.close()

        print(" ")
        print(("Summary of %s processing..." % (img_type.lower())))
        missing_counts = len(missing_img_names)
        print(("  Missing  : %d" % (missing_counts)))
        if missing_counts > 0:
            outfile_path = os.path.join(output_directory, "missing_%s_list.txt" % (img_type.lower()))
            print(("    Saving list of missing in %s..." % (outfile_path)))
            outf = open(outfile_path, "w")
            for missing_img_name in missing_img_names:
                outf.write("%s\n" % missing_img_name)
            outf.close()

        rejected_counts = len(rejected_img_names)
        print(("  Rejected : %d" % (rejected_counts)))
        if rejected_counts > 0:
            outfile_path = os.path.join(output_directory, "rejected_%s_list.txt" % (img_type.lower()))
            print(("    Saving list of rejected in %s..." % (outfile_path)))
            outf = open(outfile_path, "w")
            for rejected_img_name in rejected_img_names:
                outf.write("%s\n" % rejected_img_name)
            outf.close()

    if cter_mode_idx == idx_cter_mode_stack:
        return totresi[0][1], totresi[0][7], totresi[0][8], totresi[0][9], totresi[0][10], totresi[0][11]


########################################
# functions used by cter
# Later on make sure these functions don't conflict with those used
# in the cross resolution program getastcrfNOE.py
########################################

def ampcont2angle(A):
    #  convert amplitude contrast to phase shift
    pass  # IMPORTIMPORTIMPORT from math import sqrt, atan, degrees
    if (A == 100.0):
        return 90.0
    elif (A == -100.0):
        return 90.0
    elif (A < 0.0):
        return numpy.degrees(math.atan(old_div(A, numpy.sqrt(1.0e4 - A ** 2)))) + 180.0
    else:
        return numpy.degrees(math.atan(old_div(A, numpy.sqrt(1.0e4 - A ** 2))))


def angle2ampcont(phi):
    #  convert phase shift to amplitude contrast
    pass  # IMPORTIMPORTIMPORT from math import sqrt, tan, radians
    return old_div(numpy.tan(numpy.radians(phi)), numpy.sqrt(1.0 + numpy.tan(numpy.radians(phi)) ** 2)) * 100.0


def bracket_def(f, dat, x1, h):
    c = 1.618033989
    f1 = f(x1, dat)
    x2 = x1 + h
    f2 = f(x2, dat)
    # print x1,f1,x2,f2
    # Determine downhill direction and change sign of h if needed
    if f2 > f1:
        # print  h
        h = -h
        x2 = x1 + h
        f2 = f(x2, dat)
        # Check if minimum between x1 - h and x1 + h
        if f2 > f1: return x2, x1 - h
    # Search loop
    for i in range(100):
        h = c * h
        x3 = x2 + h;
        f3 = f(x3, dat)
        # print i,x1,f1,x2,f2,x3,f3
        if f3 > f2: return x1, x3
        x1 = x2;
        x2 = x3
        f1 = f2;
        f2 = f3
    print("Bracket did not find a mimimum")
    return None, x3


def bracket(f, dat, h):
    c = 1.618033989
    x1 = 0.0
    f1 = f(x1, dat)
    x2 = x1 + h
    f2 = f(x2, dat)
    # Search loop
    for i in range(100):
        h = c * h
        x3 = x2 + h
        f3 = f(x3, dat)
        # print i,x1,f1,x2,f2,x3,f3
        if f3 > f2: return x1, x3
        x1 = x2;
        x2 = x3
        f1 = f2;
        f2 = f3
    print("Bracket did not find a mimimum")


def goldsearch_astigmatism(f, dat, a, b, tol=1.0e-3):
    pass  # IMPORTIMPORTIMPORT from math import log, ceil
    nIter = int(numpy.ceil(-2.078087 * numpy.log(old_div(tol , abs(b - a)))))  # Eq. (10.4)
    R = 0.618033989
    C = 1.0 - R
    # First telescoping
    x1 = R * a + C * b;
    x2 = C * a + R * b
    f1 = f(x1, dat);
    f2 = f(x2, dat)
    # Main loop
    for i in range(nIter):
        if f1 > f2:
            a = x1
            x1 = x2;
            f1 = f2
            x2 = C * a + R * b;
            f2 = f(x2, dat)
        else:
            b = x2
            x2 = x1;
            f2 = f1
            x1 = R * a + C * b;
            f1 = f(x1, dat)
    if f1 < f2:
        return x1, f1
    else:
        return x2, f2


def defocus_baseline_fit(roo, i_start, i_stop, nrank, iswi):
    """
		    iswi = 2 using polynomial n rank to fit envelope function
			iswi = 3 using polynomial n rank to fit background
			The background fit is done between i_start, i_stop, but the entire baseline curve is evaluated and subtracted
	"""
    pass  # IMPORTIMPORTIMPORT import numpy as np
    pass  # IMPORTIMPORTIMPORT from morphology import imf_params_cl1

    TMP = imf_params_cl1(roo[i_start:i_stop], nrank, iswi)
    nroo = len(roo)
    baseline = np.zeros(nroo, np.float32)
    ord = len(TMP[-1])
    freqstep = TMP[0][1]
    if (iswi == 3):
        for i in range(len(roo)):
            freq = freqstep * (i - i_start)
            tmp = TMP[-1][-1]
            for j in range(1, ord): tmp += TMP[-1][j - 1] * freq ** j
            baseline[i] = tmp
    else:
        for i in range(i_start, i_stop):
            freq = freqstep * (i - i_start)
            tmp = TMP[-1][-1]
            for j in range(1, ord): tmp += TMP[-1][j - 1] * freq ** j
            baseline[i] = tmp
    return np.exp(baseline)


def simpw1d(defocus, data):
    pass  # IMPORTIMPORTIMPORT import numpy as np
    pass  # IMPORTIMPORTIMPORT from morphology import ctf_2
    pass  # IMPORTIMPORTIMPORT from utilities import generate_ctf
    pass  # IMPORTIMPORTIMPORT from math import sqrt

    # [defocus, cs, voltage, apix, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle]
    #  data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
    # data[1] - envelope
    # ct = data[1]*np.array( ctf_1d(data[2], generate_ctf([defocus, data[4], data[5], data[6], 0.0, data[7], 0.0, 0.0]), doabs= True)[data[8]:data[9]], np.float32)
    ct = data[1] * np.array(
        ctf_2(data[2], sparx_utilities.generate_ctf([defocus, data[4], data[5], data[6], 0.0, data[7], 0.0, 0.0]))[
        data[8]:data[9]], np.float32)
    # print  " 1d  ",sum(data[0]*ct),np.linalg.norm(ct,2)
    return old_div(-sum(data[0] * ct) , np.linalg.norm(ct, 2))


def simpw1d_pap(defocus, data):
    pass  # IMPORTIMPORTIMPORT import numpy as np
    pass  # IMPORTIMPORTIMPORT from morphology import ctf_1d
    pass  # IMPORTIMPORTIMPORT from utilities import generate_ctf
    pass  # IMPORTIMPORTIMPORT from math import sqrt

    # [defocus, cs, voltage, apix, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle]
    #  data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
    # data[1] - envelope
    # ct = data[1]*np.array( ctf_1d(data[2], generate_ctf([defocus, data[4], data[5], data[6], 0.0, data[7], 0.0, 0.0]), doabs= True)[data[8]:data[9]], np.float32)
    ct = np.array(
        ctf_1d(data[2], sparx_utilities.generate_ctf([defocus, data[4], data[5], data[6], 0.0, data[7], 0.0, 0.0]),
               doabs=True)[data[8]:data[9]], np.float32)
    # print  " 1d  ",sum(data[0]*ct),np.linalg.norm(ct,2)
    return old_div(-sum(old_div(data[0] * ct, data[1])), np.linalg.norm(ct, 2))


def simpw1d_print(defocus, data):
    pass  # IMPORTIMPORTIMPORT import numpy as np
    pass  # IMPORTIMPORTIMPORT from morphology import ctf_1d
    pass  # IMPORTIMPORTIMPORT from utilities import generate_ctf
    pass  # IMPORTIMPORTIMPORT from math import sqrt

    # [defocus, cs, voltage, apix, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle]
    #  data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start, i_stop]
    # data[1] - envelope
    # ct = data[1]*np.array( ctf_1d(data[2], generate_ctf([defocus, data[4], data[5], data[6], 0.0, data[7], 0.0, 0.0]), doabs= True)[data[8]:data[9]], np.float32)
    ct = np.array(
        ctf_1d(data[2], sparx_utilities.generate_ctf([defocus, data[4], data[5], data[6], 0.0, data[7], 0.0, 0.0]),
               doabs=True)[data[8]:data[9]], np.float32)
    # print  " 1d  ",sum(data[0]*ct),np.linalg.norm(ct,2)
    for i in range(len(data[0])):  print(i, i + data[8], data[0][i], ct[i], data[1][i],old_div(data[0][i], data[1][i]))
    return old_div(-sum( old_div(data[0] * ct , data[1])), np.linalg.norm(ct, 2))


def movingaverage(data, window_size, skip=3):
    pass  # IMPORTIMPORTIMPORT import numpy as np
    ld = len(data)
    qt = old_div(sum(data[skip:skip + 4]) , 3.0)
    tt = type(data[0])
    qt = np.concatenate((np.array([qt] * (window_size + skip), tt), data[skip:], np.tile(data[-1], (window_size))))
    out = np.empty(ld, np.float32)
    nc1 = window_size - old_div(window_size , 2)
    nc2 = window_size + old_div(window_size, 2) + 1
    for i in range(ld):   out[i] = sum(qt[i + nc1:i + nc2])
    return out * np.float32(old_div(1.0, window_size))


def defocusgett(roo, nx, voltage=300.0, Pixel_size=1.0, Cs=2.0, ampcont=0.1, f_start=-1.0, f_stop=-1.0, round_off=1.0,
                nr1=3, nr2=6, parent=None, DEBug=False):
    """

		1. Estimate envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get
		   defocus which matches the extracted CTF imprints
	"""
    pass  # IMPORTIMPORTIMPORT from utilities  import generate_ctf
    pass  # IMPORTIMPORTIMPORT import numpy as np
    pass  # IMPORTIMPORTIMPORT from morphology import ctf_2, bracket_def, defocus_baseline_fit, ctflimit, simpw1d, goldsearch_astigmatism

    # print "CTF params:", voltage, Pixel_size, Cs, wgh, f_start, f_stop, round_off, nr1, nr2, parent

    if f_start == 0:
        i_start = 0
    else:
        i_start = int(Pixel_size * nx * f_start + 0.5)
    if f_stop <= f_start:
        i_stop = len(roo)
        adjust_fstop = True
    else:
        i_stop = min(len(roo), int(Pixel_size * nx * f_stop + 0.5))
        adjust_fstop = False

    nroo = len(roo)

    if DEBug:  print("f_start, i_start, f_stop, i_stop:", f_start, i_start, f_stop, i_stop - 1)
    # TE  = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr1), 4)
    # baseline = defocus_baseline_fit(roo, i_start, i_stop, int(nr2), 3)

    baseline = defocus_baseline_fit(roo, i_start, nroo, int(nr2), 3)
    subpw = np.array(roo, np.float32) - baseline
    subpw[0] = subpw[1]
    # write_text_file([roo,baseline,subpw],"dbg.txt")
    # print "IN defocusgett  ",np.min(subpw),np.max(subpw)
    for i in range(len(subpw)):  subpw[i] = max(subpw[i], 0.0)
    # print "IN defocusgett  ",np.min(subpw),np.max(subpw)
    # envelope = movingaverage(  subpw   , nroo//4, 3)
    envelope = np.array([1.0] * len(subpw), np.float32)
    # write_text_file([roo,baseline,subpw],"dbgt.txt")

    # print "IN defocusgett  ",np.min(subpw),np.max(subpw),np.min(envelope)
    # envelope = np.ones(nroo, np.float32)
    defocus = 0.0
    data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start,
            i_stop]
    # for i in xrange(nroo):
    #	print  i,"   ",roo[i],"   ",baseline[i],"   ",subpw[i],"   ",envelope[i]
    h = 0.1
    # def1, def2 = bracket(simpw1d, data, h)
    # if DEBug:  print "first bracket ",def1, def2,simpw1d(def1, data),simpw1d(def2, data)
    # def1=0.1
    ndefs = 18
    defound = []
    for idef in range(ndefs):
        def1 = (idef + 1) * 0.5
        def1, def2 = bracket_def(simpw1d, data, def1, h)
        if DEBug:  print("second bracket ", idef, def1, def2, simpw1d(def1, data), simpw1d(def2, data), h)
        def1, val2 = goldsearch_astigmatism(simpw1d, data, def1, def2, tol=1.0e-3)
        if DEBug:  print("golden ", idef, def1, val2, simpw1d(def1, data))
        if def1 > 0.0:  defound.append([val2, def1])
    defound.sort()
    del defound[3:]
    def1 = defound[0][1]
    if adjust_fstop:
        pass  # IMPORTIMPORTIMPORT from morphology import ctflimit
        newstop, fnewstop = ctflimit(nx, def1, Cs, voltage, Pixel_size)
        if DEBug:  print("newstop  ", int(newstop * 0.7), fnewstop * 0.7, i_stop, newstop)
        if (newstop != i_stop and (newstop - i_start) > min(10, (i_stop - i_start))):
            i_stop = newstop
            data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont,
                    i_start, i_stop]
            """Multiline Comment15"""
            h = 0.05
            for idef in range(3):
                def1, def2 = bracket_def(simpw1d, data, defound[idef][1], h)
                if DEBug:  print(" adjusted def ", def1, def2)
                def1, val2 = goldsearch_astigmatism(simpw1d, data, def1, def2, tol=1.0e-3)
                if DEBug:  print("adjusted golden ", def1, val2, simpw1d(def1, data))
                if def1 > 0.0:  defound[idef] = [val2, def1]
            defound.sort()
            def1 = defound[0][1]
    if DEBug: print(" ultimate defocus", def1, defound)

    # defocus = defocus_guessn(Res_roo, voltage, Cs, Pixel_size, ampcont, i_start, i_stop, 2, round_off)
    # print simpw1d(def1, data),simpw1d(4.372, data)
    """Multiline Comment16"""
    if DEBug:
        qm = 1.e23
        toto = []
        for i in range(1000, 100000, 5):
            dc = old_div(float(i) , 10000.0)
            qt = simpw1d(dc, data)
            toto.append([dc, qt])
            if (qt < qm):
                qm = qt
                defi = dc
        pass  # IMPORTIMPORTIMPORT from utilities import write_text_row
        sparx_utilities.write_text_row(toto, "toto1.txt")
        print(" >>>>>>>>>  ", defi, simpw1d(defi, data))  # ,generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont])
    # def1 = defi
    # exit()
    ctf2 = ctf_2(nx, sparx_utilities.generate_ctf([def1, Cs, voltage, Pixel_size, 0.0, ampcont]))

    return def1, subpw, ctf2, baseline, envelope, i_start, i_stop


def defocusgett_pap(roo, nx, voltage=300.0, Pixel_size=1.0, Cs=2.0, ampcont=0.1, f_start=-1.0, f_stop=-1.0,
                    round_off=1.0, nr1=3, nr2=6, parent=None, DEBug=False):
    """

		1. Estimate envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get
		   defocus which matches the extracted CTF imprints
	"""
    pass  # IMPORTIMPORTIMPORT from utilities  import generate_ctf
    pass  # IMPORTIMPORTIMPORT import numpy as np
    pass  # IMPORTIMPORTIMPORT from morphology import ctf_1d, bracket_def, defocus_baseline_fit, ctflimit, simpw1d_pap, goldsearch_astigmatism

    # print "CTF params:", voltage, Pixel_size, Cs, wgh, f_start, f_stop, round_off, nr1, nr2, parent

    if f_start == 0:
        i_start = 0
    else:
        i_start = int(Pixel_size * nx * f_start + 0.5)
    if f_stop <= f_start:
        i_stop = len(roo)
        adjust_fstop = True
    else:
        i_stop = min(len(roo), int(Pixel_size * nx * f_stop + 0.5))
        adjust_fstop = False

    nroo = len(roo)

    if DEBug:  print("f_start, i_start, f_stop, i_stop:", f_start, i_start, f_stop, i_stop - 1)
    # TE  = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr1), 4)
    # baseline = defocus_baseline_fit(roo, i_start, i_stop, int(nr2), 3)

    baseline = defocus_baseline_fit(roo, i_start, nroo, int(nr2), 3)
    subpw = np.array(roo, np.float32) - baseline
    subpw[0] = subpw[1]
    # write_text_file([roo,baseline,subpw],"dbg.txt")
    # print "IN defocusgett  ",np.min(subpw),np.max(subpw)
    for i in range(len(subpw)):  subpw[i] = max(subpw[i], 0.0)
    # print "IN defocusgett  ",np.min(subpw),np.max(subpw)
    # envelope = movingaverage(  subpw   , nroo//4, 3)
    envelope = np.array([1.0] * len(subpw), np.float32)
    # write_text_file([roo,baseline,subpw],"dbgt.txt")

    # print "IN defocusgett  ",np.min(subpw),np.max(subpw),np.min(envelope)
    # envelope = np.ones(nroo, np.float32)
    defocus = 0.0
    data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start,
            i_stop]
    # for i in xrange(nroo):
    #	print  i,"   ",roo[i],"   ",baseline[i],"   ",subpw[i],"   ",envelope[i]
    h = 0.1
    # def1, def2 = bracket(simpw1d, data, h)
    # if DEBug:  print "first bracket ",def1, def2,simpw1d(def1, data),simpw1d(def2, data)
    # def1=0.1
    ndefs = 18
    defound = []
    for idef in range(ndefs):
        def1 = (idef + 1) * 0.5
        def1, def2 = bracket_def(simpw1d_pap, data, def1, h)
        if DEBug:  print("second bracket ", idef, def1, def2, simpw1d(def1, data), simpw1d_pap(def2, data), h)
        def1, val2 = goldsearch_astigmatism(simpw1d_pap, data, def1, def2, tol=1.0e-3)
        if DEBug:  print("golden ", idef, def1, val2, simpw1d_pap(def1, data))
        if def1 > 0.0:  defound.append([val2, def1])
    defound.sort()
    del defound[3:]
    def1 = defound[0][1]
    if adjust_fstop:
        pass  # IMPORTIMPORTIMPORT from morphology import ctflimit
        newstop, fnewstop = ctflimit(nx, def1, Cs, voltage, Pixel_size)
        if DEBug:  print("newstop  ", int(newstop * 0.7), fnewstop * 0.7, i_stop, newstop)
        if (newstop != i_stop and (newstop - i_start) > min(10, (i_stop - i_start))):
            i_stop = newstop
            data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont,
                    i_start, i_stop]
            """Multiline Comment17"""
            h = 0.05
            for idef in range(3):
                def1, def2 = bracket_def(simpw1d_pap, data, defound[idef][1], h)
                if DEBug:  print(" adjusted def ", def1, def2)
                def1, val2 = goldsearch_astigmatism(simpw1d_pap, data, def1, def2, tol=1.0e-3)
                if DEBug:  print("adjusted golden ", def1, val2, simpw1d_pap(def1, data))
                if (def1 > 0.0):  defound[idef] = [val2, def1]
            defound.sort()
            def1 = defound[0][1]
    if DEBug: print(" ultimate defocus", def1, defound)

    # defocus = defocus_guessn(Res_roo, voltage, Cs, Pixel_size, ampcont, i_start, i_stop, 2, round_off)
    # print simpw1d(def1, data),simpw1d(4.372, data)
    """Multiline Comment18"""
    if DEBug:
        qm = 1.e23
        toto = []
        for i in range(1000, 100000, 5):
            dc = old_div(float(i) , 10000.0)
            qt = simpw1d_pap(dc, data)
            toto.append([dc, qt])
            if (qt < qm):
                qm = qt
                defi = dc
        pass  # IMPORTIMPORTIMPORT from utilities import write_text_row
        sparx_utilities.write_text_row(toto, "toto1.txt")
        print(" >>>>>>>>>  ", defi, simpw1d(defi, data))  # ,generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont])
    # def1 = defi
    # exit()
    ctf2 = ctf_1d(nx, sparx_utilities.generate_ctf([def1, Cs, voltage, Pixel_size, 0.0, ampcont]), doabs=True)

    return def1, subpw, ctf2, baseline, envelope, i_start, i_stop


def fastigmatism3(amp, data):
    pass  # IMPORTIMPORTIMPORT from morphology import ctf2_rimg
    pass  # IMPORTIMPORTIMPORT from utilities  import generate_ctf
    pass  # IMPORTIMPORTIMPORT from alignment  import ornq
    pass  # IMPORTIMPORTIMPORT from math       import sqrt
    #  data[0] - crefim
    #  data[1] - numr
    #  data[2] - nx (image is square)
    #  data[8] - astigmatism amplitude
    #  data[9] - mask defining the region of interest

    cnx = old_div(data[2], 2) + 1
    # qt = 0.5*nx**2
    # B = 0.0
    pc = ctf2_rimg(data[2], sparx_utilities.generate_ctf([data[3], data[4], data[5], data[6], 0.0, data[7], amp, 0.0]))
    # st = Util.infomask(pc, data[9], True)
    # Util.mul_scalar(pc, 1.0/st[0])
    ang, sxs, sys, mirror, peak = sparx_alignment.ornq(pc, data[0], [0.0, 0.0], [0.0, 0.0], 1, "H", data[1], cnx, cnx)
    # print  ang, sxs, sys, mirror, peak
    # exit()
    data[8] = ang
    return -peak


fastigmatism2 = fastigmatism3


def fastigmatism3_pap(amp, data):
    pass  # IMPORTIMPORTIMPORT from morphology import ctf_rimg
    pass  # IMPORTIMPORTIMPORT from utilities  import generate_ctf
    pass  # IMPORTIMPORTIMPORT from math       import sqrt
    #  data[0] - crefim
    #  data[1] - numr
    #  data[2] - nx (image is square)
    #  data[8] - astigmatism amplitude
    #  data[9] - mask defining the region of interest

    cnx = old_div(data[2], 2) + 1
    # qt = 0.5*nx**2
    # B = 0.0
    pc = ctf_rimg(data[2], sparx_utilities.generate_ctf([data[3], data[4], data[5], data[6], 0.0, data[7], amp, 0.0]),
                  sign=0)
    # st = Util.infomask(pc, data[9], True)
    # Util.mul_scalar(pc, 1.0/st[0])
    ang, sxs, sys, mirror, peak = ornq_vpp(pc, data[0], [0.0, 0.0], [0.0, 0.0], 1, "H", data[1], cnx, cnx)
    # print  ang, sxs, sys, mirror, peak
    # exit()
    data[8] = ang
    return -peak


def simctf2(dz, data):
    pass  # IMPORTIMPORTIMPORT from morphology import ctf2_rimg
    pass  # IMPORTIMPORTIMPORT from utilities import generate_ctf

    # nx = data[2]
    # qt = 0.5*nx**2
    # print  data
    # print dz, data[4], data[5], data[6], 0.0, data[7], data[3], data[8]
    pc = ctf_rimg(data[2],
                  sparx_utilities.generate_ctf([dz, data[4], data[5], data[6], 0.0, data[7], data[3], data[8]]))
    bcc = pc.cmp("dot", data[0], {"mask": data[1], "negative": 0, "normalize": 1})
    # print " simctf2   ",amp,-bcc
    return -bcc


def simctf2_pap(dz, data):
    pass  # IMPORTIMPORTIMPORT from morphology import ctf_rimg
    pass  # IMPORTIMPORTIMPORT from utilities import generate_ctf

    # nx = data[2]
    # qt = 0.5*nx**2
    # print  data
    # print dz, data[4], data[5], data[6], 0.0, data[7], data[3], data[8]
    pc = ctf_rimg(data[2],
                  sparx_utilities.generate_ctf([dz, data[4], data[5], data[6], 0.0, data[7], data[3], data[8]]), sign=0)
    bcc = pc.cmp("dot", data[0], {"mask": data[1], "negative": 0, "normalize": 1})
    # print " simctf2   ",amp,-bcc
    return -bcc


def fupw(args, data):
    pass  # IMPORTIMPORTIMPORT from morphology import fastigmatism3
    return -fastigmatism3(args[1],
                          [data[0], data[1], data[2], args[0], data[4], data[5], data[6], data[7], data[8], data[9]])


def fupw_pap(args, data):
    pass  # IMPORTIMPORTIMPORT from morphology import fastigmatism3_pap
    return -fastigmatism3_pap(args[1], [data[0], data[1], data[2], args[0], data[4], data[5], data[6], data[7], data[8],
                                        data[9]])


########################################
# end of functions used by ctfer
########################################

########################################################################
# start of code used for estimation of cross resolution
########################################################################

def cter_vpp(input_image_path, output_directory, selection_list=None, wn=512,
             pixel_size=-1.0, Cs=2.0, voltage=300.0, wgh=10.0, f_start=-1.0, f_stop=-1.0, \
             kboot=16, overlap_x=50, overlap_y=50, edge_x=0, edge_y=0, \
             check_consistency=False, stack_mode=False, debug_mode=False, \
             program_name="cter_vpp() in morphology.py", vpp_options=[], \
             RUNNING_UNDER_MPI=False, main_mpi_proc=0, my_mpi_proc_id=0, n_mpi_procs=1):
    """Multiline Comment28"""
    pass  # IMPORTIMPORTIMPORT from   EMAN2 import periodogram
    pass  # IMPORTIMPORTIMPORT from   EMAN2db import db_check_dict, db_parse_path
    pass  # IMPORTIMPORTIMPORT from   applications import MPI_start_end
    pass  # IMPORTIMPORTIMPORT from   utilities import read_text_file, write_text_file, get_im, model_blank, model_circle, amoeba, generate_ctf
    pass  # IMPORTIMPORTIMPORT from   utilities import if_error_then_all_processes_exit_program
    pass  # IMPORTIMPORTIMPORT from   utilities import wrap_mpi_bcast
    pass  # IMPORTIMPORTIMPORT from   fundamentals import tilemic, rot_avg_table, resample
    pass  # IMPORTIMPORTIMPORT from   morphology   import threshold, bracket_def, bracket, goldsearch_astigmatism
    pass  # IMPORTIMPORTIMPORT from   morphology   import defocus_baseline_fit, simpw1d, movingaverage, localvariance, defocusgett
    pass  # IMPORTIMPORTIMPORT from   morphology   import defocus_guessn, defocusget_from_crf, make_real
    pass  # IMPORTIMPORTIMPORT from   morphology   import fastigmatism, fastigmatism1, fastigmatism2, fastigmatism3, simctf, simctf2, simctf2out, fupw,ctf2_rimg
    pass  # IMPORTIMPORTIMPORT from   alignment    import Numrinit, ringwe
    pass  # IMPORTIMPORTIMPORT from   statistics   import table_stat
    pass  # IMPORTIMPORTIMPORT from   pixel_error  import angle_ave
    pass  # IMPORTIMPORTIMPORT from   global_def   import ERROR
    pass  # IMPORTIMPORTIMPORT import global_def

    pass  # IMPORTIMPORTIMPORT from   sys import exit
    pass  # IMPORTIMPORTIMPORT import numpy as np
    pass  # IMPORTIMPORTIMPORT import os
    pass  # IMPORTIMPORTIMPORT import glob
    pass  # IMPORTIMPORTIMPORT from   time import time
    pass  # IMPORTIMPORTIMPORT from   inspect      import currentframe, getframeinfo
    pass  # IMPORTIMPORTIMPORT from   mpi import MPI_COMM_WORLD, mpi_barrier

    # ====================================================================================
    # Prepare processing
    # ====================================================================================
    #  vpp_options = [defocus_min,  defocus_max,  defocus_step,  phase_min,  phase_max,  phase_step]
    # ------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------
    # Find the CTER Running Mode before checking error conditions
    # ------------------------------------------------------------------------------------
    i_enum = -1;
    idx_cter_mode_invalid = i_enum;
    i_enum += 1;
    idx_cter_mode_all_mics = i_enum  # All Micrographs Mode - Process all s in a directory
    i_enum += 1;
    idx_cter_mode_selected_mics = i_enum  # Selected Micrographs Mode - Process all s in a selection list file
    i_enum += 1;
    idx_cter_mode_single_mic = i_enum  # Single Micrograph Mode - Process a single
    i_enum += 1;
    idx_cter_mode_stack = i_enum  # Stack Mode - Process a stack (Advanced Option)
    i_enum += 1;
    idx_cter_mode_counts = i_enum

    cter_mode_idx = idx_cter_mode_invalid
    cter_mode_name = None
    if stack_mode == False:
        # One of three Micrograph Modes
        # For any of Micrograph Modes, input image file name must be a file path pattern containing wild card "*"
        if selection_list == None:
            # User did not use selection list option
            # -> All Micrographs Mode
            cter_mode_idx = idx_cter_mode_all_mics
            cter_mode_name = "All Micrographs Mode"
        else:
            if os.path.splitext(selection_list)[1] == ".txt":
                # User specified a selection list text file path containing".txt" extension through selection list option
                # -> Selected Micrographs Mode
                cter_mode_idx = idx_cter_mode_selected_mics
                cter_mode_name = "Selected Micrographs Mode"
            else:
                # User specified an image file path (a non-text file path) through selection list option
                # -> Single Micrograph Mode
                cter_mode_idx = idx_cter_mode_single_mic
                cter_mode_name = "Single Micrograph Mode"
    else:
        # (Particle) Stack Mode
        cter_mode_idx = idx_cter_mode_stack
        cter_mode_name = "Stack Mode"

    if my_mpi_proc_id == main_mpi_proc:
        print(" ")
        print(("----- Running with %s -----" % (cter_mode_name)))

    # ------------------------------------------------------------------------------------
    # Check mode-dependent error conditions of input arguments and options if abort is necessary. All nodes do this checking
    # ------------------------------------------------------------------------------------
    error_message_list = []  # List of error messages. If no error is found, the length should be zero
    if not stack_mode:

        # Check error conditions applicable to any of Micrograph Mode
        if input_image_path.find("*") == -1:
            error_message_list.append(
                "Input image file path (%s) for %s must be a  path pattern containing wild card (*). Please check input_image_path argument." % (
                input_image_path, cter_mode_name))

        if input_image_path[:len("bdb:")].lower() == "bdb:":
            error_message_list.append(
                "BDB file can not be selected as input image file path (%s) for %s. Please check input_image_path argument and convert the image format." % (
                input_image_path, cter_mode_name))

        # Check error conditions applicable to Selected Micrographs Mode
        if cter_mode_idx == idx_cter_mode_selected_mics:
            if not os.path.exists(selection_list):
                error_message_list.append(
                    "Selection list text file specified by selection_list option (%s) for %s does not exists. Please check selection_list option." % (
                    selection_list, cter_mode_name))

        if cter_mode_idx == idx_cter_mode_single_mic:
            if not os.path.exists(os.path.join(os.path.dirname(input_image_path), os.path.basename(selection_list))):
                error_message_list.append(
                    "Micrograph specified by selection_list option (%s) for %s does not exist. Please check selection_list option." % (
                    selection_list, cter_mode_name))
            #
            if RUNNING_UNDER_MPI and n_mpi_procs != 1:
                error_message_list.append(
                    "%s supports only a single processor version. Please change MPI settings." % (cter_mode_name))

    else:
        # Check error conditions
        if input_image_path.find("*") != -1:
            error_message_list.append(
                "Stack file path specified by input_image_path (%s) for %s should not contain wild card (*). Please check input_image_path argument." % (
                input_image_path, cter_mode_name))

        is_not_found_input_image_file = False
        if input_image_path[:len("bdb:")].lower() == "bdb:":
            if not EMAN2db.db_check_dict(input_image_path):
                is_not_found_input_image_file = True
        else:
            if not os.path.exists(input_image_path):
                is_not_found_input_image_file = True
        if is_not_found_input_image_file:
            error_message_list.append(
                "Stack file specified by input_image_path (%s) for %s does not exist. Please check input_image_path argument." % (
                input_image_path, cter_mode_name))

        if RUNNING_UNDER_MPI and n_mpi_procs != 1:
            error_message_list.append(
                "%s supports only a single processor version. Please change MPI settings." % (cter_mode_name))

    # --------------------------------------------------------------------------------
    # check output-related error conditions (mode-independent). All nodes do this checking
    # --------------------------------------------------------------------------------
    if os.path.exists(output_directory):
        error_message_list.append(
            "Output directory (%s) exists already. Please check output_directory argument." % (output_directory))

    # --------------------------------------------------------------------------------
    # Check error conditions of options (mode-independent). All nodes do this checking
    # --------------------------------------------------------------------------------
    if pixel_size <= 0.0:
        error_message_list.append(
            "Pixel size (%f) must not be negative. Please set a pasitive value larger than 0.0 to pixel_size option." % (
                pixel_size))

    if wn <= 0.0:
        error_message_list.append(
            "CTF window size (%d) must not be negative. Please set a valid value larger than 0 to wn option." % (wn))

    # --------------------------------------------------------------------------------
    # Print all error messages and abort the process if necessary.
    # --------------------------------------------------------------------------------
    error_status = None
    if len(error_message_list) > 0:
        # Detected error! Print all error messages
        if my_mpi_proc_id == main_mpi_proc:
            print(" ")
            for error_message in error_message_list:
                print(("ERROR!!! %s" % (error_message)))
        error_status = ("Detected %d error(s) related to arguments and options. Run %s -h for help. Exiting..." % (
        len(error_message_list), program_name), inspect.getframeinfo(inspect.currentframe()))
    sparx_utilities.if_error_then_all_processes_exit_program(error_status)
    if RUNNING_UNDER_MPI:
        # Wait for all mpi processes to check error conditions, especially existence of output directory
        # Without this barrier, main mpi process can create output directory before some child mpi process check this error.
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    del error_message_list  # Don't need this anymore

    # ------------------------------------------------------------------------------------
    # Check warning conditions of options
    # ------------------------------------------------------------------------------------
    if my_mpi_proc_id == main_mpi_proc:
        if stack_mode:
            if selection_list != None:
                print(" ")
                print(("WARNING!!! --selection_list option will be ignored in %s." % (cter_mode_name)))
            if wn != 512:
                print(" ")
                print(("WARNING!!! --wn option will be ignored in %s." % (cter_mode_name)))
            if overlap_x != 50:
                print(" ")
                print(("WARNING!!! --overlap_x option will be ignored in %s." % (cter_mode_name)))
            if overlap_y != 50:
                print(" ")
                print(("WARNING!!! --overlap_y option will be ignored in %s." % (cter_mode_name)))
            if edge_x != 0:
                print(" ")
                print(("WARNING!!! --edge_x option will be ignored in %s." % (cter_mode_name)))
            if edge_y != 0:
                print(" ")
                print(("WARNING!!! --edge_y option will be ignored in %s." % (cter_mode_name)))
            if check_consistency:
                print(" ")
                print(("WARNING!!! --check_consistency option will be ignored in %s." % (cter_mode_name)))

    # ====================================================================================
    # Create the input file path list and also check input-related error conditions if abort is necessary.
    # ====================================================================================
    input_file_path_list = []
    if not stack_mode:
        # --------------------------------------------------------------------------------
        # Prepare the variables for all sections in  mode case
        # --------------------------------------------------------------------------------
        # Micrograph basename pattern (directory path is removed from  path pattern)
        mic_pattern = input_image_path
        mic_basename_pattern = os.path.basename(mic_pattern)

        # Global entry dictionary (all possible entries from all lists) for all mic id substring
        global_entry_dict = {}  # mic id substring is the key
        subkey_input_mic_path = "Input Micrograph Path"
        subkey_selected_mic_basename = "Selected Micrograph Basename"

        # List keeps only id substrings of s whose all necessary information are available
        valid_mic_id_substr_list = []

        # --------------------------------------------------------------------------------
        # Obtain the list of  id sustrings using a single CPU (i.e. main mpi process)
        # --------------------------------------------------------------------------------
        # NOTE: Toshio Moriya 2016/11/15
        # The below is not a real while.
        # It gives if-statements an opportunity to use break when errors need to be reported
        # However, more elegant way is to use 'raise' statement of exception mechanism...
        #
        error_status = None
        while my_mpi_proc_id == main_mpi_proc:
            # --------------------------------------------------------------------------------
            # Prepare variables for this section
            # --------------------------------------------------------------------------------
            # Prefix and suffix of  basename pattern
            # to find the head/tail indices of  id substring
            mic_basename_tokens = mic_basename_pattern.split('*')
            # Find head index of  id substring
            mic_id_substr_head_idx = len(mic_basename_tokens[0])

            # --------------------------------------------------------------------------------
            # Register  id substrings found in the input directory (specified by  path pattern)
            # to the global entry dictionary
            # --------------------------------------------------------------------------------
            # Generate the list of  paths in the input directory
            print(" ")
            print("Checking the input directory...")
            input_mic_path_list = glob.glob(mic_pattern)
            # Check error condition of input  file path list
            print(("Found %d micrographs in %s." % (len(input_mic_path_list), os.path.dirname(mic_pattern))))
            if len(input_mic_path_list) == 0:
                # The result shouldn't be empty if the specified  file name pattern is invalid
                error_status = (
                "There are no micrographs whose paths match with the specified file path pattern (%s) for %s. Please check input_image_path. Run %s -h for help." % (
                mic_pattern, cter_mode_name, program_name), inspect.getframeinfo(inspect.currentframe()))
                break

            # Register  id substrings to the global entry dictionary
            for input_mic_path in input_mic_path_list:
                # Find tail index of  id substring and extract the substring from the  name
                input_mic_basename = os.path.basename(input_mic_path)
                mic_id_substr_tail_idx = input_mic_basename.index(mic_basename_tokens[1])
                mic_id_substr = input_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
                if not mic_id_substr in global_entry_dict:
                    # print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from input_mic_path_list " % (mic_id_substr))
                    global_entry_dict[mic_id_substr] = {}
                global_entry_dict[mic_id_substr][subkey_input_mic_path] = input_mic_path

            # --------------------------------------------------------------------------------
            # Register  id substrings found in the selection list
            # to the global entry dictionary
            # --------------------------------------------------------------------------------
            # Generate the list of selected  paths in the selection file
            selected_mic_path_list = []
            # Generate  lists according to the execution mode
            if cter_mode_idx == idx_cter_mode_all_mics:
                # Treat all s in the input directory as selected ones
                selected_mic_path_list = input_mic_path_list
            else:
                if os.path.splitext(selection_list)[1] == ".txt":
                    print(" ")
                    print("Checking the selection list...")
                    selected_mic_path_list = sparx_utilities.read_text_file(selection_list)

                    # Check error condition of  entry lists
                    print(("Found %d microgarph entries in %s." % (len(selected_mic_path_list), selection_list)))
                    if len(selected_mic_path_list) == 0:
                        error_status = (
                        "The provided  list file (%s) for %s mode contains no entries. Please check selection_list option and make sure the file contains a  list. Run %s -h for help." % (
                        selection_list, cter_mode_name, program_name), inspect.getframeinfo(inspect.currentframe()))
                        break
                else:
                    print(" ")
                    print(("Processing a single micorgprah: %s..." % (selection_list)))
                    selected_mic_path_list = [selection_list]

                selected_mic_directory = os.path.dirname(selected_mic_path_list[0])
                if selected_mic_directory != "":
                    print(("    NOTE: Program disregards the directory paths in the selection list (%s)." % (
                        selected_mic_directory)))

            # Register  id substrings to the global entry dictionary
            for selected_mic_path in selected_mic_path_list:
                # Find tail index of  id substring and extract the substring from the  name
                selected_mic_basename = os.path.basename(selected_mic_path)
                mic_id_substr_tail_idx = selected_mic_basename.index(mic_basename_tokens[1])
                mic_id_substr = selected_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
                if not mic_id_substr in global_entry_dict:
                    # print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from selected_mic_path_list " % (mic_id_substr))
                    global_entry_dict[mic_id_substr] = {}
                global_entry_dict[mic_id_substr][subkey_selected_mic_basename] = selected_mic_basename

            # --------------------------------------------------------------------------------
            # Clean up variables related to registration to the global entry dictionary
            # --------------------------------------------------------------------------------
            del mic_basename_tokens
            del mic_id_substr_head_idx

            # --------------------------------------------------------------------------------
            # Create the list containing only valid  id substrings
            # --------------------------------------------------------------------------------
            # Prepare lists to keep track of invalid (rejected) s
            no_input_mic_id_substr_list = []

            print(" ")
            print("Checking the input datasets consistency...")

            # Loop over substring id list
            for mic_id_substr in global_entry_dict:
                mic_id_entry = global_entry_dict[mic_id_substr]

                warinnig_messages = []
                # selected  basename must have been registed always .
                if subkey_selected_mic_basename in mic_id_entry:
                    # Check if associated input  exists
                    if not subkey_input_mic_path in mic_id_entry:
                        input_mic_path = mic_pattern.replace("*", mic_id_substr)
                        warinnig_messages.append("    associated input  %s." % (input_mic_path))
                        no_input_mic_id_substr_list.append(mic_id_substr)

                    if len(warinnig_messages) > 0:
                        print(("WARNING!!! Micrograph ID %s does not have:" % (mic_id_substr)))
                        for warinnig_message in warinnig_messages:
                            print(warinnig_message)
                        print("    Ignores this as an invalid entry.")
                    else:
                        # print("MRK_DEBUG: adding mic_id_substr := ", mic_id_substr)
                        valid_mic_id_substr_list.append(mic_id_substr)
            # 	# This entry is not in the selection list. Do nothing

            # Check the input dataset consistency and save the result to a text file, if necessary.
            if check_consistency:
                # Create output directory
                os.mkdir(output_directory)

                # Open the consistency check file
                inconsist_mic_list_path = os.path.join(output_directory, "inconsist_mic_id_file.txt")
                print(" ")
                print(("Generating the input datasets consistency report in %s..." % (inconsist_mic_list_path)))
                inconsist_mic_list_file = open(inconsist_mic_list_path, "w")
                inconsist_mic_list_file.write("# The information about inconsistent  IDs\n")
                # Loop over substring id list
                for mic_id_substr in global_entry_dict:
                    mic_id_entry = global_entry_dict[mic_id_substr]

                    consistency_messages = []
                    # Check if associated input  path exists
                    if not subkey_input_mic_path in mic_id_entry:
                        input_mic_path = mic_pattern.replace("*", mic_id_substr)
                        consistency_messages.append("    associated input  %s." % (input_mic_path))

                    # Check if associated selected  basename exists
                    if not subkey_selected_mic_basename in mic_id_entry:
                        input_mic_path = mic_pattern.replace("*", mic_id_substr)
                        consistency_messages.append("    associated selected  %s." % (input_mic_path))

                    if len(consistency_messages) > 0:
                        inconsist_mic_list_file.write("Micrograph ID %s does not have:\n" % (mic_id_substr))
                        for consistency_message in consistency_messages:
                            inconsist_mic_list_file.write(consistency_message)
                            inconsist_mic_list_file.write("\n")

                # Close the consistency check file, if necessary
                inconsist_mic_list_file.flush()
                inconsist_mic_list_file.close()

            # Since mic_id_substr is once stored as the key of global_entry_dict and extracted with the key order
            # we need sort the valid_mic_id_substr_list here
            if debug_mode: print(("BEFORE SORT: valid_mic_id_substr_list := ", valid_mic_id_substr_list))
            valid_mic_id_substr_list.sort(key=str.lower)  # Sort list of  IDs using case insensitive string comparison
            if debug_mode: print(("AFTER SORT: valid_mic_id_substr_list := ", valid_mic_id_substr_list))

            # --------------------------------------------------------------------------------
            # Print out the summary of input consistency
            # --------------------------------------------------------------------------------
            print(" ")
            print("Summary of dataset consistency check...")
            print(("  Detected  IDs               : %6d" % (len(global_entry_dict))))
            print(("  Entries in input directory  : %6d" % (len(input_mic_path_list))))
            print(("  Entries in selection list   : %6d" % (len(selected_mic_path_list))))
            print(("  Rejected by no input        : %6d" % (len(no_input_mic_id_substr_list))))
            print(("  Valid Entries               : %6d" % (len(valid_mic_id_substr_list))))

            # --------------------------------------------------------------------------------
            # Check MPI error condition
            # --------------------------------------------------------------------------------
            if len(valid_mic_id_substr_list) < n_mpi_procs:
                error_status = (
                "Number of MPI processes (%d) supplied by --np in mpirun cannot be greater than %d (number of valid s that satisfy all criteria to be processed). Run %s -h for help." % (
                n_mpi_procs, len(valid_mic_id_substr_list, program_name)), inspect.getframeinfo(inspect.currentframe()))
                break

            # --------------------------------------------------------------------------------
            # Create input file path list
            # --------------------------------------------------------------------------------
            for mic_id_substr in valid_mic_id_substr_list:
                mic_path = global_entry_dict[mic_id_substr][subkey_input_mic_path]
                input_file_path_list.append(mic_path)

            # --------------------------------------------------------------------------------
            # Clean up variables related to tracking of invalid (rejected) s
            # --------------------------------------------------------------------------------
            del input_mic_path_list
            del selected_mic_path_list
            del no_input_mic_id_substr_list

            break

        # --------------------------------------------------------------------------------
        # Clean up the variables for all sections in  mode case
        # --------------------------------------------------------------------------------
        del mic_pattern
        del mic_basename_pattern
        del global_entry_dict
        del subkey_input_mic_path
        del subkey_selected_mic_basename
        del valid_mic_id_substr_list

        # --------------------------------------------------------------------------------
        # Print all error messages and abort the process if necessary.
        # --------------------------------------------------------------------------------
        # NOTE: Toshio Moriya 2016/11/15
        # The following function takes care of the case when an if-statement uses break for occurence of an error.
        # However, more elegant way is to use 'exception' statement of exception mechanism...
        #
        sparx_utilities.if_error_then_all_processes_exit_program(error_status)

    else:
        input_file_path_list.append(input_image_path)

    if RUNNING_UNDER_MPI:
        # Wait for main mpi process to create the input file path list
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

        # All mpi processes should know input file path list
        input_file_path_list = sparx_utilities.wrap_mpi_bcast(input_file_path_list, main_mpi_proc)

    # ====================================================================================
    # Prepare input file path(s)
    # ====================================================================================
    #
    # NOTE: 2016/03/17 Toshio Moriya
    # From here on, stack (and namics) will be used to distinguish stack mode and  mode.
    # However, a single input_file_path_list should be sufficient since we already know the mode.
    # Let's consider this refactoring in the future.
    #
    stack = None  # (particle) stack file name: if it is not None, cter runs with stack mode. Otherwise, runs with  mode
    namics = []  # file name list
    if not stack_mode:
        namics = input_file_path_list
        if debug_mode: print(("BEFORE SORT: namics := ", namics))
        namics.sort(key=str.lower)  # Sort list of s using case insensitive string comparison
        if debug_mode: print(("AFTER SORT: namics := ", namics))
    else:
        stack = input_file_path_list[0]

    del input_file_path_list  # Don't need this anymore

    # Make output directory
    outpwrot = "%s/pwrot" % (output_directory)
    if stack == None:
        outmicthumb = "%s/micthumb" % (output_directory)
    if debug_mode:
        outravg = "%s/ravg" % (output_directory)
    if my_mpi_proc_id == main_mpi_proc:
        # Make output directory
        if not os.path.exists(output_directory):
            os.mkdir(output_directory)
        os.mkdir(outpwrot)
        if stack == None:
            os.mkdir(outmicthumb)
        if debug_mode:
            os.mkdir(outravg)

    if RUNNING_UNDER_MPI:
        # Make all mpi processes wait for main mpi process to create output directory
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    # Set up loop variables depending on the cter mode
    if stack == None:
        if RUNNING_UNDER_MPI:
            set_start, set_end = sparx_applications.MPI_start_end(len(namics), n_mpi_procs, my_mpi_proc_id)
        else:
            set_start = 0
            set_end = len(namics)
    else:
        pw2 = []
        data = EMAN2_cppwrap.EMData.read_images(stack)
        nima = len(data)
        for i in range(nima):
            pw2.append(EMAN2_cppwrap.periodogram(data[i]))
        wn = pw2[0].get_xsize()
        set_start = 0
        set_end = 1

    # Set up progress message
    if my_mpi_proc_id == main_mpi_proc:
        print(" ")
        print("Estimating CTF parameters...")
        if stack == None:
            print("  Micrographs processed by main process (including percent of progress):")
            progress_percent_step = ( old_div(
                                                set_end - set_start) , 100.0)  # the number of micrograms for main mpi processer divided by 100

    totresi = []
    missing_img_names = []
    rejected_img_names = []
    for ifi in range(set_start, set_end):
        if stack == None:
            img_type = "Micrograph"
            img_name = namics[ifi]
            img_basename_root = os.path.splitext(os.path.basename(img_name))[0]

            if my_mpi_proc_id == main_mpi_proc:
                print(("    Processing %s ---> %6.2f%%" % (img_name, old_div((ifi - set_start), progress_percent_step))))

            if not os.path.exists(img_name):
                missing_img_names.append(img_name)
                print(
                    "    %s %s: Can not find this file. Skipping the estimation and CTF parameters are not stored..." % (
                    img_type, img_name))
                continue
            mic = sparx_utilities.get_im(img_name)
            try:
                pw2 = sparx_fundamentals.tilemic(mic, win_size=wn, overlp_x=overlap_x, overlp_y=overlap_y,
                                                 edge_x=edge_x, edge_y=edge_y)
            except:
                print(
                    "MRK_DEBUG: tilemic() in cter_vpp() raised an exception. The micrographs {} might have a problem. Please check it and remove it if necessary.".format(
                        img_name))
                raise
            del mic

        else:
            img_type = "Stack"
            img_name = stack

            numFM = EMAN2_cppwrap.EMUtil.get_image_count(img_name)
            pw2 = [None] * numFM
            for i in range(numFM):
                pw2[i] = EMAN2_cppwrap.periodogram(sparx_utilities.get_im(img_name, i))

        for i in range(len(pw2)):
            pw2[i] = square_root(pw2[i])

        if debug_mode: print("    %s %s: Process %04d started the processing. Detected %d image(s) in this %s file." % (
        img_type, img_name, ifi, img_type.lower()))

        if EMAN2db.db_check_dict(img_name) == False:
            img_basename_root = os.path.splitext(os.path.basename(img_name))[0]
        else:
            path, dictname, keys = EMAN2db.db_parse_path(img_name)
            img_basename_root = dictname

        #  VPP code starts here  03/08/2017
        nimi = len(pw2)
        adefocus = [0.0] * kboot
        aampcont = [0.0] * kboot
        aamplitu = [0.0] * kboot
        aangle = [0.0] * kboot

        allroo = []
        for imi in range(nimi):
            allroo.append(sparx_fundamentals.rot_avg_table(pw2[imi]))
        lenroo = len(allroo[0])

        for nboot in range(kboot):
            # at = time()
            if (nboot == 0):
                boot = list(range(nimi))
            else:
                pass  # IMPORTIMPORTIMPORT from random import randint
                for imi in range(nimi): boot[imi] = random.randint(0, nimi - 1)
            qa = sparx_utilities.model_blank(wn, wn)
            roo = np.zeros(lenroo, np.float32)
            sroo = np.zeros(lenroo, np.float32)
            aroo = np.zeros(lenroo, np.float32)

            for imi in range(nimi):
                EMAN2_cppwrap.Util.add_img(qa, pw2[boot[imi]])
                temp1 = np.array(allroo[boot[imi]])
                roo += temp1
                temp2 = movingaverage(temp1, 10)
                aroo += temp2
                sroo += temp2 ** 2
            sroo[0] = sroo[1]
            aroo[0] = aroo[1]
            sroo = old_div((sroo - old_div(aroo ** 2 , nimi)) , nimi)
            aroo = old_div(aroo, nimi)
            roo = old_div(roo, nimi)
            qa = old_div(qa,nimi)

            if f_start < 0:
                #  Find a break point
                bp = 1.e23
                for i in range(5, lenroo - 5):
                    # t1 = linreg(sroo[:i])
                    # t2 = linreg(sroo[i:])
                    # tt = t1[1][0] + t2[1][0]
                    xtt = np.array(list(range(i)), np.float32)
                    zet = np.poly1d(np.polyfit(xtt, sroo[:i], 2))
                    t1 = sum((sroo[:i] - zet(xtt)) ** 2)
                    xtt = np.array(list(range(i, lenroo)), np.float32)
                    zet = np.poly1d(np.polyfit(xtt, sroo[i:], 2))
                    tt = t1 + sum((sroo[i:] - zet(xtt)) ** 2)
                    if tt < bp:
                        bp = tt
                        istart = i
                f_start = old_div(istart , (pixel_size * wn))
            """Multiline Comment29"""
            # write_text_file([roo.tolist(),aroo.tolist(),sroo.tolist()], "sroo%03d.txt"%ifi)
            rooc = roo.tolist()

            defc, ampcont, subpw, baseline, envelope, istart, istop = defocusgett_vpp(rooc, wn, voltage=voltage,
                                                                                      Pixel_size=pixel_size, Cs=Cs, \
                                                                                      f_start=f_start, f_stop=f_stop,
                                                                                      vpp_options=vpp_options, nr1=3,
                                                                                      nr2=6, parent=None,
                                                                                      DEBug=debug_mode)
            if debug_mode:
                print("  RESULT 1 %s" % (img_name), nboot, defc, ampcont, istart, istop)  # , (time()-at)/60.
                """Multiline Comment30"""
            # mpi_barrier(MPI_COMM_WORLD)

            # exit()
            bg = baseline.tolist()
            en = envelope.tolist()

            bckg = sparx_utilities.model_blank(wn, wn, 1, 1)
            envl = sparx_utilities.model_blank(wn, wn, 1, 1)

            pass  # IMPORTIMPORTIMPORT from math import sqrt
            nc = old_div(wn, 2)
            ne = istop
            ns = istart - 1
            bg.append(bg[-1])
            en.append(en[-1])
            for i in range(wn):
                for j in range(wn):
                    r = numpy.sqrt((i - nc) ** 2 + (j - nc) ** 2)
                    ir = int(r + 0.5)
                    if (ir < nc):
                        # This is awkward but it is needed for storing of results
                        dr = r - ir
                        bckg.set_value_at(i, j, (1. - dr) * bg[ir] + dr * bg[ir + 1])
                        if (ir > ns and ir < ne):
                            dr = r - ir
                            envl.set_value_at(i, j, (1. - dr) * en[ir] + dr * en[ir + 1])

            qse = old_div(threshold((qa - bckg)) , envl)
            # print  "  fit1  ", nboot,(time()-at)/60.0
            # at = time()
            # (qse*mask).write_image("rs2.hdf")
            # bckg.write_image("bckg.hdf")
            # envl.write_image("envl.hdf")
            # qse.write_image("qse.hdf")
            del envl, baseline, envelope
            # exit()
            ##  SIMULATION
            # bang = 0.7
            # qse = ctf2_rimg(wn, generate_ctf([defc,Cs,voltage,pixel_size,0.0,wgh, bang, 37.0]) )
            # qse.write_image("rs3.hdf")
            # at = time()
            defc, ampcont, astamp, astang, score = defocusgett_vpp2(qse, wn, defc, ampcont, voltage=voltage,
                                                                    Pixel_size=pixel_size, Cs=Cs, i_start=istart,
                                                                    i_stop=istop, parent=None, DEBug=debug_mode)
            if debug_mode:
                print("  RESULT 2 %s" % (img_name), nboot, defc, ampcont, astamp, astang, score)  # , (time()-at)/60.
                """Multiline Comment31"""
            adefocus[nboot] = defc
            aampcont[nboot] = ampcont
            aamplitu[nboot] = astamp
            aangle[nboot] = astang
        # from sys import exit
        # exit()
        # from morphology import ctf_rimg, ctf_1d
        # cq = ctf_1d(wn, generate_ctf([defc, Cs, voltage, pixel_size, 0.0, ampcont,astamp,astang]), doabs = True)[20:150]
        # write_text_file([subpw[20:150],cq],"pwds%02d.txt"%nboot)
        # print  "  fit2  ", nboot,(time()-at)/60.0

        # print  "  xxx2  ", nboot,(time()-at)/60.0

        ad1, ad2, ad3, ad4 = sparx_statistics.table_stat(adefocus)  # return values: average, variance, minimum, maximum
        if ad2 <= 0.0:
            print((
                              "    %s %s: Detected the variance less than zero (defocus statics: avg = %f, var = %f, min = %f, max = %f)." % (
                      img_type, img_name, ad1, ad2, ad3, ad4)))
            print(("           The program ignores this estimate..."))
            continue

        reject = []
        thr = 3 * numpy.sqrt(ad2)
        for i in range(len(adefocus)):
            if (abs(adefocus[i] - ad1) > thr):
                print((
                                  "    %s %s: Rejected an outlier defocus estimate (defocus = %f, average defocus = %f, threshold = %f)." % (
                          img_type, img_name, adefocus[i], ad1, thr)))
                reject.append(i)

        if (len(reject) > 0):
            print(("    %s %s: Total number of rejects %s" % (img_type, img_name, len(reject))))
            for i in range(len(reject) - 1, -1, -1):
                del adefocus[i]
                del aampcont[i]
                del aamplitu[i]
                del aangle[i]

        # print  "  xxx2  ", nboot,(time()-at)/60.0
        if (len(adefocus) < 2):
            print((
                              "    %s %s: After rejection of outliers, there is too few estimated defocus values. Skipping the estimation and CTF parameters are not stored..." % (
                      img_type, img_name)))
        else:
            # print "adefocus",adefocus
            # print  "aamplitu",aamplitu
            # print "aangle",aangle
            ad1, ad2, ad3, ad4 = sparx_statistics.table_stat(adefocus)
            #  compute statistics of ampcont using phase shifts instead
            ed1, ed2 = sparx_pixel_error.angle_ave([2 * ampcont2angle(q) for q in
                                                    aampcont])  # Have to use this trick as the function works for range [0,360]
            ed1 = old_div(ed1,2)
            ed2 = old_div(ed2,2)
            bd1, bd2, bd3, bd4 = sparx_statistics.table_stat(aamplitu)
            cd1, cd2 = sparx_pixel_error.angle_ave(
                [2 * q for q in aangle])  # Have to use this trick as the function works for range [0,360]
            cd1 = old_div(cd1,2)
            cd2 = old_div(cd2,2)
            temp = 0.0
            stdavad1 = np.sqrt(kboot * max(0.0, ad2))
            ed2 *= np.sqrt(kboot)
            stdavbd1 = np.sqrt(kboot * max(0.0, bd2))
            cd2 *= np.sqrt(kboot)

            # Adjust value ranges of astig. amp. and angle.
            if (bd1 < 0.0):
                bd1 = -bd1
                cd1 += 90.0
            cd1 = cd1 % 180

            #  SANITY CHECK, do not produce anything if defocus abd astigmatism amplitude are out of whack
            reject_img_messages = []
            try:
                pwrot2 = rotavg_ctf(sparx_utilities.model_blank(wn, wn), ad1, Cs, voltage, pixel_size, bd1, cd1)
            except:
                reject_img_messages.append(
                    "    - Astigmatism amplitude (%f) is larger than defocus (%f) or defocus (%f) is negative." % (
                    bd1, ad1, ad1))
            """Multiline Comment32"""

            if len(reject_img_messages) > 0:
                rejected_img_names.append(img_name)
                print("    %s %s: Rejected the CTF estimate - " % (img_type, img_name), ad1, Cs, voltage, pixel_size,
                      angle2ampcont(ed1), bd1, cd1, "(def, Cs, vol, apix, amp_contrast, astig_amp, astig_angle)")
                print("    %s %s: because... " % (img_type, img_name))
                for reject_img_message in reject_img_messages:
                    print(reject_img_message)
                print("    %s %s: Skipping the estimation and CTF parameters are not stored..." % (img_type, img_name))
            else:
                #  Estimate the point at which (sum_errordz ctf_1(dz+errordz))^2 falls to 0.5
                import random as rqt
                #  NOW WE SWITCH PHASE SHIFT TO AMPLITUDE CONTRAST
                ed1 = angle2ampcont(ed1)
                # at = time()
                supe = sparx_utilities.model_blank(wn, wn)
                niter = 1000
                for it in range(niter):
                    EMAN2_cppwrap.Util.add_img(supe,
                                               EMAN2_cppwrap.Util.ctf_rimg(wn, wn, 1, ad1 + rqt.gauss(0.0, stdavad1),
                                                                           pixel_size, voltage, Cs, 0.0,
                                                                           ed1 + angle2ampcont(rqt.gauss(0.0, ed2)),
                                                                           bd1 + rqt.gauss(0.0, stdavbd1),
                                                                           cd1 + rqt.gauss(0.0, cd2), 1))
                ni = old_div(wn , 2)
                supe = old_div(supe, niter)
                pwrot2 = rotavg_ctf(supe, ad1, Cs, voltage, pixel_size, bd1, cd1)
                for i in range(ni):  pwrot2[i] = pwrot2[i] ** 2
                ibec = 0
                for it in range(ni - 1, 0, -1):
                    if pwrot2[it] > 0.5:
                        ibec = it
                        break
                # print  "  uuu1 ",(time()-at)/60.
                #  NOW WE SWITCH PHASE SHIFT TO AMPLITUDE CONTRAST
                ed2 = angle2ampcont(ed2)
                if (ed2 < 0.0): ed2 = 180.0 - ed2

                pass  # IMPORTIMPORTIMPORT from morphology import ctf_1d
                ct = sparx_utilities.generate_ctf([ad1, Cs, voltage, pixel_size, temp, ed1, 0.0, 0.0])
                cq = ctf_1d(wn, ct)
                # at = time()
                supe = [0.0] * ni
                niter = 1000
                for i in range(niter):
                    cq = sparx_utilities.generate_ctf(
                        [ad1 + rqt.gauss(0.0, stdavad1), Cs, voltage, pixel_size, 0.0, ed1, 0.0, 0.0])
                    ci = ctf_1d(wn, cq)[:ni]
                    for l in range(ni):  supe[l] += ci[l]

                for l in range(ni):  supe[l] = (old_div(supe[l] , niter)) ** 2

                ib1 = 0
                for it in range(ni - 1, 0, -1):
                    if supe[it] > 0.5:
                        ib1 = it
                        break
                ibec = old_div(ibec, (pixel_size * wn))  # with astigmatism
                ib1 = old_div(ib1, (pixel_size * wn))  # no astigmatism
                # print  " error est  ",(time()-at)/60.0
                # from utilities import write_text_file
                # write_text_file([range(ni), supe[:ni],pwrot2[:ni]],"fifi.txt")

                # Compute defocus CV and astig. amp. CV (CV: coefficient of variation; ratio of error (SD) relative to average (mean))
                #  I blocked them, make no sense here.
                # if ad1 < max(0.0, valid_min_defocus): ERROR("Logical Error: Encountered unexpected defocus value (%f). Consult with the developer." % (ad1), "%s in %s" % (__name__, os.path.basename(__file__))) # MRK_ASSERT
                # if stdavad1 < 0.0: ERROR("Logical Error: Encountered unexpected defocus SD value (%f). Consult with the developer." % (stdavad1), "%s in %s" % (__name__, os.path.basename(__file__))) # MRK_ASSERT
                cvavad1 = old_div(stdavad1, ad1) * 100  # use percentage

                # if bd1 < 0.0: ERROR("Logical Error: Encountered unexpected astig. amp. value (%f). Consult with the developer." % (bd1), "%s in %s" % (__name__, os.path.basename(__file__))) # MRK_ASSERT
                if stdavbd1 < 0.0: sparx_global_def.ERROR(
                    "Logical Error: Encountered unexpected astig. amp. SD value (%f). Consult with the developer." % (
                        stdavbd1), "%s in %s" % (__name__, os.path.basename(__file__)))  # MRK_ASSERT

                bd1 = max(bd1, 1.0e-15)
                cvavbd1 = old_div(stdavbd1, bd1) * 100  # use percentage

                # Compute CTF limit (theoretical resolution limit based on the oscillations of CTF)
                # For output, use ctflim (relative frequency limit [1/A]), not ctflim_abs (absolute frequency limit)
                #
                # NOTE: 2016/03/23 Toshio Moriya
                # xr is limiting frequency [1/A]. Max is Nyquist frequency = 1.0/(2*apix[A/pixel]). <UNIT: [1/(A/pixel)/[pixel])] => [(pixel)/(A*pixel] => [1/A]>
                # 1.0/xr is limiting period (Angstrom resolution) [A]. Min is Nyquist period = (2*apix[A/pixel]). <UNIT: [1/(1/A)] = [A]>
                # fwpix is width of Fourier pixel [pixel/A] := 1.0[pixel]/(2*apix[A/pixel])/box_half[pixel] = 1[pixel]/fullsize[A]). <UNIT: [pixel/(A/pixel)/(pixel)] = [pixel*(pixel/A)*(1/pixel) = [pixel/A]>
                # int(xr/fwpix+0.5) is limiting_absolute_frequency [1/pixel]. <Unit:[(1/A)/(pixel/A)] = [(1/A)*(A/pixel)] = [1/pixel]>
                # return  int(xr/fwpix+0.5),xr, which is limiting_abs_frequency [1/pixel], and Limiting_frequency[1/A]
                #
                ctflim_abs, ctflim = ctflimit(wn, ad1, Cs, voltage, pixel_size)

                """Multiline Comment33"""
                # print  " error est2  ",(time()-at)/60.0
                # print " ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1 ",ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1
                lnsb = len(subpw)
                try:
                    crot1 = rotavg_ctf(
                        ctf_rimg(wn, sparx_utilities.generate_ctf([ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1]),
                                 sign=0), ad1, Cs, voltage, pixel_size, 0.0, 0.0)[:lnsb]
                except:
                    crot1 = [0.0] * lnsb
                try:
                    pwrot1 = rotavg_ctf(threshold(qa - bckg), ad1, Cs, voltage, pixel_size, 0.0, 0.0)[:lnsb]
                except:
                    pwrot1 = [0.0] * lnsb
                try:
                    crot2 = rotavg_ctf(
                        ctf_rimg(wn, sparx_utilities.generate_ctf([ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1]),
                                 sign=0), ad1, Cs, voltage, pixel_size, bd1, cd1)[:lnsb]
                except:
                    crot2 = [0.0] * lnsb
                try:
                    pwrot2 = rotavg_ctf(threshold(qa - bckg), ad1, Cs, voltage, pixel_size, bd1, cd1)[:lnsb]
                except:
                    pwrot2 = [0.0] * lnsb
                #  #1 - rotational averages without astigmatism, #2 - with astigmatism
                lnsb = min(lnsb, len(crot2), len(pwrot1), len(crot2), len(pwrot2))
                sparx_utilities.write_text_file(
                    [list(range(lnsb)), [old_div(old_div(float(i) , wn), pixel_size) for i in range(lnsb)], pwrot1[:lnsb], crot1[:lnsb],
                     pwrot2[:lnsb], crot2[:lnsb]], os.path.join(outpwrot, "%s_rotinf.txt" % (img_basename_root)))

                #
                # print  " error est3  ",(time()-at)/60.0
                # NOTE: 2016/03/23 Toshio Moriya
                # Compute mean of extrema differences (differences at peak & trough) between
                # (1) experimental rotational average with astigmatism (pwrot2)
                # (2) experimental rotational average without astigmatism (pwrot1), and
                # as a indication of goodness of astigmatism estimation by cter.
                # The peak & trough detection uses fitted rotational average with astigmatism (crot2)
                # Start from 1st trough while ignoring 1st peak.
                # End at astigmatism frequency limit.
                #
                """Multiline Comment34"""
                #				#if extremum_counts == 0: ERROR("Logical Error: Encountered unexpected zero extremum counts. Consult with the developer." % (bd1), "%s in %s" % (__name__, os.path.basename(__file__))) # MRK_ASSERT
                #				extremum_diff_avg = 1.1#extremum_diff_sum / extremum_counts

                #				#print "MRK_DEBUG: extremum_avg = %12.5g, extremum_diff_sum = %12.5g, extremum_counts = %03d," % (extremum_avg, extremum_diff_sum, extremum_counts)
                #				#print "MRK_DEBUG: extremum_diff_avg = %12.5g, extremum_diff_sum = %12.5g, extremum_counts = %03d," % (extremum_diff_avg, extremum_diff_sum, extremum_counts)

                #				if stack == None:     cmd = "echo " + "    " + namics[ifi] + "  >>  " + fou
                #				else:                 cmd = "echo " + "    " + "  >>  " + fou
                #				os.system(cmd)
                cvavbd1 = old_div(stdavbd1 , bd1) * 100  # use percentage

                max_freq = old_div(0.5 , pixel_size)  # dummy value for maximum frequency. set to Nyquist frequency for now. let's add the implementation in near future (Toshio 2017/12/06)
                reserved = 0.0  # dummy value for reserved spot, which might be used for parameter of external programs (e.g. CTFFIND4, GCTF, and etc.)
                # wgh                     # constant amplitude contrast provided by user (default 10%)
                phase_shift = ampcont2angle(ed1) - ampcont2angle(
                    wgh)  # Volta phase shift [deg] = total amplitude contrast phase shift [deg] (ed1) -  constant amplitude contrast phase shift [deg]; ed1 is boot strap average of total amplitude contrast [%]

                if debug_mode: print((
                                                 "    %s %s: Process %04d finished the processing. Estimated CTF parmaters are stored in %s." % (
                                         img_type, img_name, ifi, os.path.join(output_directory, "partres.txt"))))
                #				if debug_mode: print((ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1, stdavad1, ed2, stdavbd1, cd2, cvavad1, cvavbd1, extremum_diff_avg, ib1, ibec, ctflim))
                #				totresi.append([img_name, ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1, stdavad1, ed2, stdavbd1, cd2, cvavad1, cvavbd1, extremum_diff_avg, ib1, ibec, ctflim])
                if debug_mode: print((ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1, stdavad1, ed2, stdavbd1, cd2,
                                      cvavad1, cvavbd1, ib1, ibec, ctflim, max_freq, reserved, wgh, phase_shift))
                totresi.append(
                    [img_name, ad1, Cs, voltage, pixel_size, temp, ed1, bd1, cd1, stdavad1, ed2, stdavbd1, cd2, cvavad1,
                     cvavbd1, ib1, ibec, ctflim, max_freq, reserved, wgh, phase_shift])

        #				if stack == None:
        #					print  namics[ifi], ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec
        #				else:
        #					print               ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec
        #				if stack == None:
        #					totresi.append( [ namics[ifi], ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec])
        #				else:
        #					totresi.append( [ 0, ad1, Cs, voltage, pixel_size, temp, wgh, bd1, cd1, stdavad1, stdavbd1, cd2, ib1, ibec])
        #				#if ifi == 4 : break
        # print  " error est4  ",(time()-at)/60.0

        if stack == None:
            img_mic = sparx_utilities.get_im(namics[ifi])
            # create  thumbnail
            nx = img_mic.get_xsize()
            if nx > 512:
                img_micthumb = sparx_fundamentals.resample(img_mic, old_div(512.0 , nx))
            else:
                img_micthumb = img_mic
            img_basename_root = os.path.splitext(os.path.basename(img_name))[0]
            fou = os.path.join(outmicthumb, "%s_thumb.hdf" % (img_basename_root))
            img_micthumb.write_image(fou)

    if RUNNING_UNDER_MPI:
        pass  # IMPORTIMPORTIMPORT from utilities import wrap_mpi_gatherv
        totresi = sparx_utilities.wrap_mpi_gatherv(totresi, 0, mpi.MPI_COMM_WORLD)
        missing_img_names = sparx_utilities.wrap_mpi_gatherv(missing_img_names, 0, mpi.MPI_COMM_WORLD)
        rejected_img_names = sparx_utilities.wrap_mpi_gatherv(rejected_img_names, 0, mpi.MPI_COMM_WORLD)

    if my_mpi_proc_id == main_mpi_proc:
        outf = open(os.path.join(output_directory, "partres.txt"), "w")
        for i in range(len(totresi)):
            for k in range(1, len(totresi[i])):
                outf.write("  %12.5g" % totresi[i][k])
            outf.write("  %s\n" % totresi[i][0])
        outf.close()

        print(" ")
        print(("Summary of %s processing..." % (img_type.lower())))
        missing_counts = len(missing_img_names)
        print(("  Missing  : %d" % (missing_counts)))
        if missing_counts > 0:
            outfile_path = os.path.join(output_directory, "missing_%s_list.txt" % (img_type.lower()))
            print(("    Saving list of missing in %s..." % (outfile_path)))
            outf = open(outfile_path, "w")
            for missing_img_name in missing_img_names:
                outf.write("%s\n" % missing_img_name)
            outf.close()

        rejected_counts = len(rejected_img_names)
        print(("  Rejected : %d" % (rejected_counts)))
        if rejected_counts > 0:
            outfile_path = os.path.join(output_directory, "rejected_%s_list.txt" % (img_type.lower()))
            print(("    Saving list of rejected in %s..." % (outfile_path)))
            outf = open(outfile_path, "w")
            for rejected_img_name in rejected_img_names:
                outf.write("%s\n" % rejected_img_name)
            outf.close()

    if cter_mode_idx == idx_cter_mode_stack:
        return totresi[0][1], totresi[0][7], totresi[0][8], totresi[0][9], totresi[0][10], totresi[0][11]


########################################
# functions used by cter_vpp
########################################

def defocusgett_vpp(roo, nx, voltage=300.0, Pixel_size=1.0, Cs=2.0, f_start=-1.0, f_stop=-1.0, vpp_options=[], nr1=3,
                    nr2=6, parent=None, DEBug=False):
    """
		1. Estimate envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get
		   defocus which matches the extracted CTF imprints
	"""
    pass  # IMPORTIMPORTIMPORT from utilities  import generate_ctf, write_text_file
    pass  # IMPORTIMPORTIMPORT import numpy as np
    pass  # IMPORTIMPORTIMPORT from morphology import defocus_baseline_fit, simpw1d, angle2ampcont

    # print "CTF params:", voltage, Pixel_size, Cs, wgh, f_start, f_stop, round_off, nr1, nr2, parent

    if f_start == 0:
        i_start = 0
    else:
        i_start = int(Pixel_size * nx * f_start + 0.5)
    if f_stop <= f_start:
        i_stop = len(roo)
        adjust_fstop = True
    else:
        i_stop = min(len(roo), int(Pixel_size * nx * f_stop + 0.5))
        adjust_fstop = False

    i_stop = min(old_div(nx, 2) - 2, i_stop)  # this is for resampling into polar

    nroo = len(roo)

    if DEBug:  print("f_start, i_start, f_stop, i_stop:", f_start, i_start, f_stop, i_stop - 1)
    # TE  = defocus_env_baseline_fit(roo, i_start, i_stop, int(nr1), 4)
    # baseline = defocus_baseline_fit(roo, i_start, i_stop, int(nr2), 3)
    baseline = defocus_baseline_fit(roo, i_start, nroo, int(nr2), 3)
    subpw = np.array(roo, np.float32) - baseline
    subpw[0] = subpw[1]
    # write_text_file([roo,baseline,subpw],"dbg.txt")
    # print "IN defocusgett  ",np.min(subpw),np.max(subpw)
    for i in range(len(subpw)):  subpw[i] = max(subpw[i], 0.0)
    # print "IN defocusgett  ",np.min(subpw),np.max(subpw)
    # envelope = movingaverage(  subpw, nroo//8, 3)

    # envelope = defocus_baseline_fit(roo, i_start, nroo, int(nr2), 2) - baseline
    envelope = defocus_baseline_fit(roo, i_start, min(int(i_stop * 1.45), old_div(nx, 2) - 2), int(nr2), 2) - baseline
    #  Process envelope
    qm = np.max(envelope[5:])
    dc = max(np.min(envelope[5:]), old_div(qm , 1000.))
    for i in range(len(envelope)):
        if (envelope[i] < dc): envelope[i] = qm

    # envelope = np.array([1.0]*len(subpw), np.float32)
    # write_text_file([roo,baseline,subpw,envelope],"dbgt.txt")

    # print "IN defocusgett  ",np.min(subpw),np.max(subpw),np.min(envelope)
    # envelope = np.ones(nroo, np.float32)
    defocus = 0.0
    ampcont = 0.0
    data = [subpw[i_start:i_stop], envelope[i_start:i_stop], nx, defocus, Cs, voltage, Pixel_size, ampcont, i_start,
            i_stop]
    qm = 1.e23
    # toto = []
    #  vpp_options = [defocus_min,  defocus_max,  defocus_step,  phase_min,  phase_max,  phase_step]
    #  This is in degrees
    if (vpp_options[4] < vpp_options[3]): vpp_options[4] += 180.0
    a = vpp_options[3]
    while (a <= vpp_options[4]):
        data[7] = angle2ampcont(a % 180.0)
        dc = vpp_options[0]
        while (dc <= vpp_options[1]):
            qt = simpw1d_pap(dc, data)
            # toto.append([a,data[7],dc,qt])
            if (qt < qm):
                qm = qt
                defi = dc
                ampcont = data[7]
            # print  a,dc,qt
            dc += vpp_options[2]
        a += vpp_options[5]
    # '''
    if DEBug:
        pass  # IMPORTIMPORTIMPORT from utilities import write_text_row
        # write_text_row(toto,"toto1.txt")
        data[7] = ampcont
        print(" >>>>>>>>>  ", defi, data[7], ampcont2angle(data[7]),
              simpw1d_print(defi, data))  # ,generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont])
    # data[7]=10.
    # defi = 4.5
    # print " >>>>>>>>>  ",defi,data[7],simpw1d_print(defi, data)#,generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont])
    # def1 = defi
    # exit()
    # '''
    # ctf2 = ctf_1d(nx, generate_ctf([defi, Cs, voltage, Pixel_size, 0.0, ampcont]), doabs= True)
    """Multiline Comment35"""
    return defi, ampcont, subpw.tolist(), baseline, envelope, i_start, i_stop  # , ctf2


def defocusgett_vpp2(qse, wn, xdefc, xampcont, voltage=300.0, Pixel_size=1.0, Cs=2.0, i_start=0, i_stop=0, parent=None,
                     DEBug=False):
    """
		1. Estimate envelope function and baseline noise using constrained simplex method
		   so as to extract CTF imprints from 1D power spectrum
		2. Based one extracted ctf imprints, perform exhaustive defocus searching to get
		   defocus which matches the extracted CTF imprints
		Switching to phase shift altogether 07/07/2017
	"""
    # from utilities  import generate_ctf
    # import numpy as np
    pass  # IMPORTIMPORTIMPORT from utilities import amoeba
    pass  # IMPORTIMPORTIMPORT from alignment import Numrinit, ringwe

    cnx = old_div(wn , 2) + 1
    cny = cnx
    mode = "H"
    numr = sparx_alignment.Numrinit(i_start, i_stop - 1, 1, mode)
    wr = sparx_alignment.ringwe(numr, mode)

    crefim = EMAN2_cppwrap.Util.Polar2Dm(qse, cnx, cny, numr, mode)
    EMAN2_cppwrap.Util.Frngs(crefim, numr)
    EMAN2_cppwrap.Util.Applyws(crefim, numr, wr)

    bdef = 0.
    bphs = 0.0  # phase_shift
    bamp = 0.0  # initial astigmatism amplitude
    bang = 0.0  # initial astigmatism angle
    astdata = [crefim, numr, wn, bdef, Cs, voltage, Pixel_size, bphs, bamp, bang]
    initial_ast_ang = 0.0
    xphase_shift = ampcont2angle(xampcont)
    dama = sparx_utilities.amoeba([xdefc, xphase_shift, initial_ast_ang], [0.1, 2.0, 0.05], fupw_vpp, 1.e-4, 1.e-4, 500,
                                  astdata)
    qma = -dama[-2]
    if DEBug: print(" amoeba  %7.2f  %7.2f  %12.6g  %12.6g" % (dama[0][0], dama[0][1] % 180.0, dama[0][2], qma))
    dpefi = dama[0][0]
    dphshift = dama[0][1] % 180.0
    dastamp = dama[0][2]
    astdata = [crefim, numr, wn, dpefi, Cs, voltage, Pixel_size, dphshift, dastamp, bang]
    junk = fastigmatism3_vpp(dama[0][2], astdata)
    #  Corrected here PAP 06/16/2018
    dastang = astdata[9]

    """Multiline Comment36"""
    if DEBug:
        # from utilities import write_text_row
        # write_text_row(toto,"toto1.txt")
        print(" repi3  ", dpefi, dphshift, dastamp, dastang, junk)

    return dpefi, angle2ampcont(dphshift), dastamp, dastang, qma  # dp


def fupw_vpp(args, data):
    pass  # IMPORTIMPORTIMPORT from morphology import fastigmatism3_vpp
    #  args = [defocus, phaseshift, astigma-amp]
    #                                   0       1     2   3     4    5         6          7     8     9
    #            (astdata) =          [crefim, numr, wn, bdef, Cs, voltage, Pixel_size, bphs, bamp, bang]
    #
    #   [data[0], data[1], data[2], args[0], data[4], data[5], data[6], args[1], data[8], data[9]]
    #   [crefim,   numr,   wn, (args)defocus, Cs,   voltage, Pixel_size,(a)phshift, (a)astamp, ang, mask]
    #
    # print  " fuw_vpp           ",args[0],args[1],args[2]
    args[0] = max(min(args[0], 6.0), 0.01)
    args[1] = args[1]  # %180.0  #  Phase shift within valid range
    args[2] = max(min(args[2], 3.0), 0.0)
    #                        (a)astamp
    return fastigmatism3_vpp(args[2],
                             [data[0], data[1], data[2], args[0], data[4], data[5], data[6], args[1], data[8], data[9]])


def fastigmatism3_vpp(amp, data):
    pass  # IMPORTIMPORTIMPORT from morphology import ctf2_rimg
    pass  # IMPORTIMPORTIMPORT from utilities  import generate_ctf
    pass  # IMPORTIMPORTIMPORT from alignment  import ornq
    pass  # IMPORTIMPORTIMPORT from math       import sqrt
    #  data[0] - crefim
    #  data[1] - numr
    #  data[2] - nx (image is square)
    #  data[8] - astigmatism amplitude
    #  data[9] - mask defining the region of interest
    #
    #      0        1          2       3        4       5         6         7      8        9
    #   [data[0], data[1], data[2], args[0], data[4], data[5], data[6], args[1], data[8], data[9]]
    #   [crefim,   numr,   wn, (args)defocus, Cs,   voltage, Pixel_size,(a)phshift, (a)astamp, ang]
    #
    #  generate_ctf
    #      0      1    2       3       4        5        6                      7
    #  [defocus, cs, voltage, apix, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle]
    #  [ microns, mm, kV, Angstroms, A^2, microns, radians]

    cnx = old_div(data[2], 2) + 1
    # qt = 0.5*nx**2
    # B = 0.0
    pc = ctf_rimg(data[2], sparx_utilities.generate_ctf(
        [data[3], data[4], data[5], data[6], 0.0, angle2ampcont(data[7] % 180.0), amp, 0.0]), sign=0.)
    # st = pc.cmp("dot", pc, dict(negative = 0, mask = data[10], normalize = 0))
    # Util.mul_scalar(pc, 1.0/st)
    ang, sxs, sys, mirror, peak = ornq_vpp(pc, data[0], [0.0, 0.0], [0.0, 0.0], 1, "H", data[1], cnx, cnx)
    # print  ang, sxs, sys, mirror, peak
    # print  " fastigmatism3_vpp ",round(data[3],3), data[7], amp,round(ang,2),round(peak,3)
    # PAP 06/16/2018 - there seems to be a mistake here as amplitude angle should be #9 on the list
    # So I corrected it
    data[9] = ang
    return peak


def ornq_vpp(image, crefim, xrng, yrng, step, mode, numr, cnx, cny, deltapsi=0.0):
    """Determine shift and rotation between image and reference image (refim)
	   no mirror
		quadratic interpolation
		cnx, cny in FORTRAN convention
	"""
    pass  # IMPORTIMPORTIMPORT from math import pi, cos, sin, radians
    pass  # IMPORTIMPORTIMPORT from alignment import ang_n
    # from utilities import info
    # print "ORNQ"
    peak = -1.0E23

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
            EMAN2_cppwrap.Util.Normalize_ring(cimage, numr, 0)
            retvals = EMAN2_cppwrap.Util.Crosrng_e(crefim, cimage, numr, 0, deltapsi)
            qn = retvals["qn"]
            if qn >= peak:
                sx = -ix
                sy = -iy
                ang = sparx_alignment.ang_n(retvals["tot"], mode, numr[-1])
                peak = qn
    # mirror is returned as zero for consistency
    mirror = 0
    co = numpy.cos(numpy.radians(ang))
    so = -numpy.sin(numpy.radians(ang))
    sxs = sx * co - sy * so
    sys = sx * so + sy * co
    return ang, sxs, sys, mirror, peak

