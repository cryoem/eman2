
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

# Comment by Zhengfan Yang on 06/11/10
# I decided to move everything related to pixel error to this file. Otherwise, they are
# scattered all over the places and there are a lot of duplications and confusions.
#
# This file contains the following functions:
#  1. pixel_error_2D (originally max_pixel_error, also I changed its input format)
#  2. max_3D_pixel_error
#  3. angle_diff
#  4. align_diff_params
#  5. align_diff
#  6. align_diff_textfile (new)
#  7. ave_ali_err_params
#  8. ave_ali_err
#  9. ave_ali_err_textfile (new)
# 10. multi_align_diff_params (new)
# 11. calc_connect_list (new)
# 12. ali_stable_list (new)
#
# Update on 09/01/10, we have decided for multiple alignment function 10 and 11 are not best way
# to determine the stability. We have instead written a new function to to this.
# 13. multi_align_stability (latest)
#
# Here, align_diff_params() and ave_ali_err_params() takes two lists of alignmnet parameters
# in the following format:
# [alpha1, sx1, sy1, mirror1, alpha2, sx2, sy2, mirror2, ...]
# align_diff() and ave_ali_err() takes two lists of images
# align_diff_textfile() and ave_ali_err_textfile() takes two textfiles, which are usually the output
# file of header() or sxheader.py.
# Per previous discussion, I decided not to support two stacks of images.

import EMAN2_cppwrap
import math
import numpy
from . import sp_global_def
from . import sp_statistics
from . import sp_utilities


def pixel_error_2D(ali_params1, ali_params2, r=1.0):
    """
	Compute average squared 2D pixel error
	"""
    return (
        old_div(
            (
                numpy.sin(old_div(numpy.radians(ali_params1[0] - ali_params2[0]), 2))
                * (2 * r + 1)
            )
            ** 2,
            2,
        )
        + (ali_params1[1] - ali_params2[1]) ** 2
        + (ali_params1[2] - ali_params2[2]) ** 2
    )


def max_3D_pixel_error(t1, t2, r=1.0):
    """
	Compute maximum pixel error between two sets of orientation parameters
	assuming object has radius r, t1 is the projection transformation
	of the first projection and t2 of the second one, respectively:
		t = Transform({"type":"spider","phi":phi,"theta":theta,"psi":psi})
		t.set_trans(Vec2f(-tx, -ty))
	Note the function is symmetric in t1, t2.
	"""
    dummy = EMAN2_cppwrap.Transform()
    if type(dummy) != type(t1):
        t = EMAN2_cppwrap.Transform(
            {"type": "spider", "phi": t1[0], "theta": t1[1], "psi": t1[2]}
        )
        t.set_trans(EMAN2_cppwrap.Vec2f(-t1[3], -t1[4]))
    else:
        t = t1
    if type(dummy) != type(t2):
        u = EMAN2_cppwrap.Transform(
            {"type": "spider", "phi": t2[0], "theta": t2[1], "psi": t2[2]}
        )
        u.set_trans(EMAN2_cppwrap.Vec2f(-t2[3], -t2[4]))
    else:
        u = t2

    return EMAN2_cppwrap.EMData().max_3D_pixel_error(t, u, r)


def angle_ave(angle1):
    """
	This function computes average angle of a set of angles.
	It also computes a measure of dispersion (incorrect).
	"""

    nima = len(angle1)

    cosi = 0.0
    sini = 0.0
    for i in range(nima):
        qt = numpy.radians(angle1[i])
        cosi += numpy.cos(qt)
        sini += numpy.sin(qt)
    alphai = numpy.degrees(math.atan2(sini, cosi)) % 360.0
    # what follows is not correct, it is just to give a measure of dispersion
    stdv = 0.0
    for i in range(nima):
        qt = angle1[i] - alphai
        if qt > 180.0:
            qt -= 360.0
        elif qt < -180.0:
            qt += 360.0
        stdv += qt * qt
    stdv = numpy.sqrt(old_div(stdv, nima))

    return alphai, stdv


def angle_diff(angle1, angle2):
    """
	This function determines the relative angle between two sets of angles.
	The resulting angle has to be added (modulo 360) to the first set.
	"""

    nima = len(angle1)
    nima2 = len(angle2)
    if nima2 != nima:
        sp_global_def.ERROR("Error: List lengths do not agree!", "angle_diff", 1)
    else:
        del nima2

    cosi = 0.0
    sini = 0.0
    for i in range(nima):
        qt = numpy.radians(angle2[i] - angle1[i])
        cosi += numpy.cos(qt)
        sini += numpy.sin(qt)
    alphai = numpy.degrees(math.atan2(sini, cosi)) % 360.0

    return alphai


def angle_diff_sym(angle1, angle2, simi=1):
    """
	This function determines the relative angle around Z axis (phi) between two sets of angles
	   taking into account point group symmetry with multiplicity simi.
	The input has to be in the form [[phi0,theta0], [phi1,theta1], ...]
	  Only sets that have theta in the same range (0,90), or (90,180) are included in calculation.
	The resulting angle has to be added (modulo 360/simi) to the first set.
	"""

    nima = len(angle1)
    if len(angle2) != nima:
        sp_global_def.ERROR("List lengths do not agree!", "angle_diff_sym", 1)

    cosi = 0.0
    sini = 0.0
    agree = 0
    for i in range(nima):
        if ((angle2[i][1] < 90.0) and (angle1[i][1] < 90.0)) or (
            (angle2[i][1] > 90.0) and (angle1[i][1] > 90.0)
        ):
            qt = numpy.radians((angle2[i][0] - angle1[i][0]) * simi)
            cosi += numpy.cos(qt)
            sini += numpy.sin(qt)
            agree += 1
    if agree == 0:
        return 0.0
    else:
        return numpy.degrees(old_div(math.atan2(sini, cosi), simi)) % (
            old_div(360.0, simi)
        )


def align_diff_params(ali_params1, ali_params2):
    """
	This function determines the relative angle, shifts and mirrorness between
	two sets of alignment parameters.
	"""

    nima = len(ali_params1)
    nima2 = len(ali_params2)
    if nima2 != nima:
        sp_global_def.sxprint("Error: Number of images do not agree!")
        return 0.0, 0.0, 0.0, 0
    else:
        nima = old_div(nima, 4)
        del nima2

    # Read the alignment parameters and determine the relative mirrorness
    mirror_same = 0
    for i in range(nima):
        if ali_params1[i * 4 + 3] == ali_params2[i * 4 + 3]:
            mirror_same += 1
    if mirror_same > old_div(nima, 2):
        mirror = 0
    else:
        mirror_same = nima - mirror_same
        mirror = 1

    # Determine the relative angle
    cosi = 0.0
    sini = 0.0
    angle1 = []
    angle2 = []
    for i in range(nima):
        mirror1 = ali_params1[i * 4 + 3]
        mirror2 = ali_params2[i * 4 + 3]
        if abs(mirror1 - mirror2) == mirror:
            alpha1 = ali_params1[i * 4]
            alpha2 = ali_params2[i * 4]
            if mirror1 == 1:
                alpha1 = -alpha1
                alpha2 = -alpha2
            angle1.append(alpha1)
            angle2.append(alpha2)
    alphai = angle_diff(angle1, angle2)

    # Determine the relative shift
    sxi = 0.0
    syi = 0.0
    for i in range(nima):
        mirror1 = ali_params1[i * 4 + 3]
        mirror2 = ali_params2[i * 4 + 3]
        if abs(mirror1 - mirror2) == mirror:
            alpha1 = ali_params1[i * 4]
            # alpha2 = ali_params2[i*4]
            sx1 = ali_params1[i * 4 + 1]
            sx2 = ali_params2[i * 4 + 1]
            sy1 = ali_params1[i * 4 + 2]
            sy2 = ali_params2[i * 4 + 2]
            alpha12, sx12, sy12, mirror12 = sp_utilities.combine_params2(
                alpha1, sx1, sy1, int(mirror1), alphai, 0.0, 0.0, 0
            )
            if mirror1 == 0:
                sxi += sx2 - sx12
            else:
                sxi -= sx2 - sx12
            syi += sy2 - sy12

    sxi = old_div(sxi, mirror_same)
    syi = old_div(syi, mirror_same)

    return alphai, sxi, syi, mirror


def multi_align_stability(
    ali_params,
    mir_stab_thld=0.0,
    grp_err_thld=10000.0,
    err_thld=1.732,
    print_individual=False,
    d=64,
):
    def sqerr(a):
        n = len(a)
        avg = sum(a)
        sq = 0.0
        for i in range(n):
            sq += a[i] ** 2
        return old_div((sq - old_div(avg * avg, n)), n)

    # args - G, data - [T, d]
    def func(args, data, return_avg_pixel_error=True):

        ali_params = data[0]
        d = data[1]
        L = len(ali_params)
        N = old_div(len(ali_params[0]), 4)

        args_list = [0.0] * (L * 3)
        for i in range(L * 3 - 3):
            args_list[i] = args[i]

        cosa = [0.0] * L
        sina = [0.0] * L
        for i in range(L):
            cosa[i] = numpy.cos(numpy.radians(args_list[i * 3]))
            sina[i] = numpy.sin(numpy.radians(args_list[i * 3]))
        sqr_pixel_error = [0.0] * N
        ave_params = []
        for i in range(N):
            sum_cosa = 0.0
            sum_sina = 0.0
            sx = [0.0] * L
            sy = [0.0] * L
            for j in range(L):
                if int(ali_params[j][i * 4 + 3]) == 0:
                    sum_cosa += numpy.cos(
                        numpy.radians(args_list[j * 3] + ali_params[j][i * 4])
                    )
                    sum_sina += numpy.sin(
                        numpy.radians(args_list[j * 3] + ali_params[j][i * 4])
                    )
                    sx[j] = (
                        args_list[j * 3 + 1]
                        + ali_params[j][i * 4 + 1] * cosa[j]
                        + ali_params[j][i * 4 + 2] * sina[j]
                    )
                    sy[j] = (
                        args_list[j * 3 + 2]
                        - ali_params[j][i * 4 + 1] * sina[j]
                        + ali_params[j][i * 4 + 2] * cosa[j]
                    )
                else:
                    sum_cosa += numpy.cos(
                        numpy.radians(-args_list[j * 3] + ali_params[j][i * 4])
                    )
                    sum_sina += numpy.sin(
                        numpy.radians(-args_list[j * 3] + ali_params[j][i * 4])
                    )
                    sx[j] = (
                        -args_list[j * 3 + 1]
                        + ali_params[j][i * 4 + 1] * cosa[j]
                        - ali_params[j][i * 4 + 2] * sina[j]
                    )
                    sy[j] = (
                        args_list[j * 3 + 2]
                        + ali_params[j][i * 4 + 1] * sina[j]
                        + ali_params[j][i * 4 + 2] * cosa[j]
                    )
            sqrtP = numpy.sqrt(sum_cosa ** 2 + sum_sina ** 2)
            sqr_pixel_error[i] = max(
                0.0,
                 old_div(d *d , 4.0)  * (1 - old_div(sqrtP, L)) + sqerr(sx) + sqerr(sy),
            )
            # Get ave transform params
            H = EMAN2_cppwrap.Transform({"type": "2D"})
            H.set_matrix(
                [
                    old_div(sum_cosa, sqrtP),
                    old_div(sum_sina, sqrtP),
                    0.0,
                    old_div(sum(sx), L),
                    -old_div(sum_sina, sqrtP),
                    old_div(sum_cosa, sqrtP),
                    0.0,
                    old_div(sum(sy), L),
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                ]
            )
            dd = H.get_params("2D")
            #  We are using here mirror of the LAST SET.
            H = EMAN2_cppwrap.Transform(
                {
                    "type": "2D",
                    "alpha": dd["alpha"],
                    "tx": dd["tx"],
                    "ty": dd["ty"],
                    "mirror": int(ali_params[L - 1][i * 4 + 3]),
                    "scale": 1.0,
                }
            )
            dd = H.get_params("2D")
            ave_params.append([dd["alpha"], dd["tx"], dd["ty"], dd["mirror"]])
        # Warning: Whatever I return here is squared pixel error, this is for the easy expression of derivative
        # Don't forget to square root it after getting the value
        if return_avg_pixel_error:
            return old_div(sum(sqr_pixel_error), N)
        else:
            return sqr_pixel_error, ave_params

    """Multiline Comment0"""

    # I decided not to use scipy in order to reduce the dependency, I wrote the C++ code instead
    # from scipy import array, int32
    # from scipy.optimize.lbfgsb import fmin_l_bfgs_b

    # Find out the subset which is mirror stable over all runs
    all_part = []
    num_ali = len(ali_params)
    nima = old_div(len(ali_params[0]), 4)
    for i in range(num_ali):
        mirror0 = []
        mirror1 = []
        for j in range(nima):
            if ali_params[i][j * 4 + 3] == 0:
                mirror0.append(j)
            else:
                mirror1.append(j)
        mirror0 = numpy.array(mirror0, "int32")
        mirror1 = numpy.array(mirror1, "int32")
        all_part.append([mirror0, mirror1])
    match, stab_part, CT_s, CT_t, ST, st = sp_statistics.k_means_stab_bbenum(
        all_part, T=0, nguesses=1
    )
    mir_stab_part = stab_part[0] + stab_part[1]

    mir_stab_rate = old_div(len(mir_stab_part), float(nima))
    if mir_stab_rate <= mir_stab_thld:
        return [], mir_stab_rate, -1.0
    mir_stab_part.sort()

    del all_part, match, stab_part, CT_s, CT_t, ST, st

    # for i in xrange(len(mir_stab_part)):  print i, mir_stab_part[i]

    nima2 = len(mir_stab_part)

    # print "  mirror stable  ",nima2

    # Keep the alignment parameters of mirror stable particles
    ali_params_mir_stab = [[] for i in range(num_ali)]
    for i in range(num_ali):
        for j in mir_stab_part:
            ali_params_mir_stab[i].extend(ali_params[i][j * 4 : j * 4 + 4])
    # Find out the alignment parameters for each iteration against the last one
    args = []
    for i in range(num_ali - 1):
        alpha, sx, sy, mirror = align_diff_params(
            ali_params_mir_stab[i], ali_params_mir_stab[num_ali - 1]
        )
        args.extend([alpha, sx, sy])

    # Do an initial analysis, purge all outlier particles, whose pixel error are larger than three times the threshold
    data = [ali_params_mir_stab, d]
    pixel_error_before, ave_params = func(
        numpy.array(args), data, return_avg_pixel_error=False
    )
    #  We have to return mir_stab_rate (see above), even if the group does not survive it and the pixel error before cleaning,
    #   see below, and in ISAC print the user the overall statistics (histograms?)
    #   so one can get on overall idea how good/bad data is.  PAP  01/25/2015
    # print  " >>> ",sqrt(sum(pixel_error_before)/nima2)
    ali_params_cleaned = [[] for i in range(num_ali)]
    cleaned_part = []
    for j in range(nima2):
        pixel_error_before[j] = max(0.0, pixel_error_before[j])  # prevent sqrt of 0
        if numpy.sqrt(pixel_error_before[j]) > 3 * err_thld:
            pass  # print "  removed ",3*err_thld,j,sqrt(pixel_error_before[j])
        else:
            cleaned_part.append(mir_stab_part[j])
            for i in range(num_ali):
                ali_params_cleaned[i].extend(ali_params_mir_stab[i][j * 4 : j * 4 + 4])
    nima3 = len(cleaned_part)
    prever = numpy.sqrt(old_div(sum(pixel_error_before), nima2))
    if nima3 <= 1:
        return [], mir_stab_rate, prever

    # print "  cleaned part  ",nima3

    # Use LBFGSB to minimize the sum of pixel errors
    data = [ali_params_cleaned, d]
    # Use Python code
    # ps_lp, val, d = fmin_l_bfgs_b(func, array(args), args=[data], fprime=dfunc, bounds=None, m=10, factr=1e3, pgtol=1e-4, iprint=-1, maxfun=100)
    # Use C++ code
    ali_params_cleaned_list = []
    for params in ali_params_cleaned:
        ali_params_cleaned_list.extend(params)
    results = EMAN2_cppwrap.Util.multi_align_error(args, ali_params_cleaned_list, d)
    ps_lp = results[:-1]

    # Negative val can happen in some rare cases, it should be due to rounding errors,
    # because all results show the val is about 1e-13.
    # print "Strange results"
    # print "args =", args
    # print "ali_params_cleaned_list =", ali_params_cleaned_list
    # print "results = ", results
    val = max(0.0, results[-1])

    del ali_params_cleaned_list

    if numpy.sqrt(val) > grp_err_thld:
        return [], mir_stab_rate, numpy.sqrt(val)

    pixel_error_after, ave_params = func(ps_lp, data, return_avg_pixel_error=False)

    stable_set = []
    val = 0.0
    for i in range(nima):
        if i in cleaned_part:
            j = cleaned_part.index(i)
            err = numpy.sqrt(pixel_error_after[j])
            if err < err_thld:
                stable_set.append([err, i, ave_params[j]])
                val += err
                if print_individual:
                    sp_global_def.sxprint(
                        "Particle %4d :  pixel error = %18.4f" % (i, err)
                    )
            else:
                if print_individual:
                    sp_global_def.sxprint(
                        "Particle %4d :  pixel error = %18.4f  unstable" % (i, err)
                    )
        else:
            if print_individual:
                sp_global_def.sxprint("Particle %4d :  Mirror unstable" % i)
    #  return average pixel error before pruning as it is more informative
    return stable_set, mir_stab_rate, prever  # sqrt(val/len(stable_set))


# args - G, data - [T, d]


"""Multiline Comment1"""


# These are some obsolete codes, we retain them just in case.
"""Multiline Comment2"""


from builtins import range
