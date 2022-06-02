# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import division
from past.utils import old_div

"""
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2019 (toshio.moriya@kek.jp)
Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
Copyright (c) 2000-2006 The University of Texas - Houston Medical School

This software is issued under a joint BSD/GNU license. You may use the source
code in this file under either license. However, note that the complete EMAN2
and SPARX software packages have some GPL dependencies, so you are responsible
for compliance with the licenses of these packages if you opt to use BSD
licensing. The warranty disclaimer below holds in either instance.

This complete copyright notice must be included in any revised version of the
source code. Additional authorship citations may be added, but existing author
citations must be preserved.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 59 Temple
Place, Suite 330, Boston, MA  02111-1307 USA
"""
import EMAN2
import EMAN2_cppwrap
import EMAN2db
import datetime
import errno
import fractions
import inspect
import json
import math
import matplotlib
import matplotlib.pyplot
import mpi
import numpy
import os
import pickle
import re
import scipy.spatial

import six
from . import sp_applications
from . import sp_fundamentals
from . import sp_global_def
from . import sp_morphology
import string
import struct
import sys
import time
import traceback
import zlib


def makerelpath(p1, p2):
    """Takes a pair of paths /a/b/c/d and /a/b/e/f/g and returns a relative path to b from a, ../../e/f/g"""

    p1s = [i for i in os.path.realpath(p1).split("/") if len(i) > 0]
    p2s = [i for i in os.path.realpath(p2).split("/") if len(i) > 0]

    for dv in range(min(len(p1s), len(p2s))):
        if p1s[dv] != p2s[dv]:
            break
    else:
        dv += 1

    p1s = p1s[dv:]
    p2s = p2s[dv:]

    return "../" * len(p1s) + "/".join(p2s)


def make_v_stack_header(path, vstack_path, verbose=False):
    if vstack_path.startswith("bdb:"):
        vspath = EMAN2db.db_parse_path(vstack_path)[0]
        if vspath == "." or vspath == "./":
            vspath = EMAN2db.e2getcwd()
        vspath = os.path.join(vspath, "EMAN2DB/")
    else:
        vspath = vstack_path
        if vspath == "." or vspath == "./":
            vspath = EMAN2db.e2getcwd()

    if path.lower()[:4] == "bdb:" and not "#" in path:
        uu = os.path.split(path)
        if uu[0] == "":
            path = "bdb:.#" + path[4:]
        else:
            path = uu[0] + "#" + uu[1]
    if path.lower()[:4] != "bdb:":
        path = "bdb:" + path
    if "#" in path:
        path, dbs = path.rsplit("#", 1)
        path += "#"
        dbs = [dbs]
    else:
        if not "#" in path and path[-1] != "/":
            path += "#"
        dbs = EMAN2db.db_list_dicts(path)

    dbs.sort()

    step = [0, 1]
    vstack_list = []
    for db in dbs:
        dct, keys = EMAN2db.db_open_dict(path + db, ro=True, with_keys=True)
        if len(step) == 2:
            if keys == None:
                vals = list(range(step[0], len(dct), step[1]))
            else:
                vals = keys[
                    step[0] :: step[1]
                ]  # we apply --step even if we have a list of keys
        else:
            if keys == None:
                vals = list(range(step[0], step[2], step[1]))
            else:
                vals = keys[
                    step[0] : step[2] : step[1]
                ]  # we apply --step even if we have a list of keys

        for n in vals:
            try:
                d = dct.get(n, nodata=1).get_attr_dict()
            except:
                traceback.print_exc()
                print("---\nerror reading ", db, n)
                continue
            # This block converts an absolute path to the actual data to a relative path
            try:
                dpath = os.path.realpath(dct.get_data_path(n))
                if os.name == "nt":
                    vspath = vspath.replace("\\", "/")
                    dpath = dpath.replace("\\", "/")
                rpath = makerelpath(vspath, dpath)
            except Exception as e:
                print(e)
                print("error with data_path ", db, n)
                continue
            d["data_path"] = rpath
            d["data_n"] = n
            d["data_source"] = path + db
            if d["data_path"] == None:
                print("error with data_path ", db, n)
                continue
            vstack_list.append(d)
        dct.close()
    return vstack_list


def amoeba(
    var, scale, func, ftolerance=1.0e-4, xtolerance=1.0e-4, itmax=500, data=None
):
    """Use the simplex method to maximize a function of 1 or more variables.

	   Input:
		  var = the initial guess, a list with one element for each variable
		  scale = the search scale for each variable, a list with one
			  element for each variable.
		  func = the function to maximize.

	   Optional Input:
		  ftolerance = convergence criterion on the function values (default = 1.e-4)
		  xtolerance = convergence criterion on the variable values (default = 1.e-4)
		  itmax = maximum number of iterations allowed (default = 500).
		  data = data to be passed to func (default = None).

	   Output:
		  (varbest,funcvalue,iterations)
		  varbest = a list of the variables at the maximum.
		  funcvalue = the function value at the maximum.
		  iterations = the number of iterations used.

	   - Setting itmax to zero disables the itmax check and the routine will run
		 until convergence, even if it takes forever.
	   - Setting ftolerance or xtolerance to 0.0 turns that convergence criterion
		 off.  But do not set both ftolerance and xtolerance to zero or the routine
		 will exit immediately without finding the maximum.
	   - To check for convergence, check if (iterations < itmax).

	   The function should be defined like func(var,data) where
	   data is optional data to pass to the function.

	   Example:

		   import amoeba
		   def afunc(var,data=None): return 1.0-var[0]*var[0]-var[1]*var[1]
		   print amoeba.amoeba([0.25,0.25],[0.5,0.5],afunc)

	   Version 1.0 2005-March-28 T. Metcalf
		   1.1 2005-March-29 T. Metcalf - Use scale in simsize calculation.
						- Use func convergence *and* x convergence
						  rather than func convergence *or* x
						  convergence.
		   1.2 2005-April-03 T. Metcalf - When contracting, contract the whole
						  simplex.
	   """

    nvar = len(var)  # number of variables in the minimization
    nsimplex = nvar + 1  # number of vertices in the simplex

    # first set up the simplex

    simplex = [0] * (nvar + 1)  # set the initial simplex
    simplex[0] = var[:]
    for i in range(nvar):
        simplex[i + 1] = var[:]
        simplex[i + 1][i] += scale[i]

    fvalue = []
    for i in range(nsimplex):  # set the function values for the simplex
        fvalue.append(func(simplex[i], data=data))

        # Ooze the simplex to the maximum

    iteration = 0

    while 1:
        # find the index of the best and worst vertices in the simplex
        ssworst = 0
        ssbest = 0
        for i in range(nsimplex):
            if fvalue[i] > fvalue[ssbest]:
                ssbest = i
            if fvalue[i] < fvalue[ssworst]:
                ssworst = i

                # get the average of the nsimplex-1 best vertices in the simplex
        pavg = [0.0] * nvar
        for i in range(nsimplex):
            if i != ssworst:
                for j in range(nvar):
                    pavg[j] += simplex[i][j]
        for j in range(nvar):
            pavg[j] = old_div(pavg[j], nvar)  # nvar is nsimplex-1
        simscale = 0.0
        for i in range(nvar):
            simscale += old_div(abs(pavg[i] - simplex[ssworst][i]), scale[i])
        simscale = old_div(simscale, nvar)

        # find the range of the function values
        fscale = old_div((abs(fvalue[ssbest]) + abs(fvalue[ssworst])), 2.0)
        if fscale != 0.0:
            frange = old_div(abs(fvalue[ssbest] - fvalue[ssworst]), fscale)
        else:
            frange = 0.0  # all the fvalues are zero in this case

            # have we converged?
        if (
            (ftolerance <= 0.0 or frange < ftolerance)
            and (xtolerance <= 0.0 or simscale < xtolerance)  # converged to maximum
        ) or (  # simplex contracted enough
            itmax and iteration >= itmax
        ):  # ran out of iterations
            return simplex[ssbest], fvalue[ssbest], iteration

            # reflect the worst vertex
        pnew = [0.0] * nvar
        for i in range(nvar):
            pnew[i] = 2.0 * pavg[i] - simplex[ssworst][i]
        fnew = func(pnew, data=data)
        if fnew <= fvalue[ssworst]:
            # the new vertex is worse than the worst so shrink
            # the simplex.
            for i in range(nsimplex):
                if i != ssbest and i != ssworst:
                    for j in range(nvar):
                        simplex[i][j] = 0.5 * simplex[ssbest][j] + 0.5 * simplex[i][j]
                    fvalue[i] = func(simplex[i], data=data)
            for j in range(nvar):
                pnew[j] = 0.5 * simplex[ssbest][j] + 0.5 * simplex[ssworst][j]
            fnew = func(pnew, data=data)
        elif fnew >= fvalue[ssbest]:
            # the new vertex is better than the best so expand
            # the simplex.
            pnew2 = [0.0] * nvar
            for i in range(nvar):
                pnew2[i] = 3.0 * pavg[i] - 2.0 * simplex[ssworst][i]
            fnew2 = func(pnew2, data=data)
            if fnew2 > fnew:
                # accept the new vertex in the simplexe
                pnew = pnew2
                fnew = fnew2
                # replace the worst vertex with the new vertex
        for i in range(nvar):
            simplex[ssworst][i] = pnew[i]
        fvalue[ssworst] = fnew
        iteration += 1
        # print "Iteration:",iteration,"  ",ssbest,"  ",fvalue[ssbest]


"""Multiline Comment0"""


def compose_transform2(alpha1, sx1, sy1, scale1, alpha2, sx2, sy2, scale2):
    """Print the composition of two transformations  T2*T1
		Here  if v's are vectors:   vnew = T2*T1 vold
			 with T1 described by alpha1, sx1, scale1 etc.

	  Combined parameters correspond to image first transformed by set 1 followed by set 2.

		Usage: compose_transform2(alpha1,sx1,sy1,scale1,alpha2,sx2,sy2,scale2)
		   angles in degrees
	"""

    t1 = EMAN2_cppwrap.Transform(
        {
            "type": "2D",
            "alpha": alpha1,
            "tx": sx1,
            "ty": sy1,
            "mirror": 0,
            "scale": scale1,
        }
    )
    t2 = EMAN2_cppwrap.Transform(
        {
            "type": "2D",
            "alpha": alpha2,
            "tx": sx2,
            "ty": sy2,
            "mirror": 0,
            "scale": scale2,
        }
    )
    tt = t2 * t1
    d = tt.get_params("2D")
    return d["alpha"], d["tx"], d["ty"], d["scale"]


def compose_transform2m(
	alpha1=0.0,
	sx1=0.0,
	sy1=0.0,
	mirror1=0,
	scale1=1.0,
	alpha2=0.0,
	sx2=0.0,
	sy2=0.0,
	mirror2=0,
	scale2=1.0,
):
	"""Print the composition of two transformations  T2*T1
		Here  if v's are vectors:   vnew = T2*T1 vold
			 with T1 described by alpha1, sx1, scale1 etc.

	  Combined parameters correspond to image first transformed by set 1 followed by set 2.

		Usage: compose_transform2(alpha1,sx1,sy1,mirror1,scale1,alpha2,sx2,sy2,mirror2,scale2)
		   angles in degrees
	"""

	t1 = EMAN2_cppwrap.Transform(
		{
			"type": "2D",
			"alpha": alpha1,
			"tx": sx1,
			"ty": sy1,
			"mirror": mirror1,
			"scale": scale1,
		}
	)
	t2 = EMAN2_cppwrap.Transform(
		{
			"type": "2D",
			"alpha": alpha2,
			"tx": sx2,
			"ty": sy2,
			"mirror": mirror2,
			"scale": scale2,
		}
	)
	tt = t2 * t1
	d = tt.get_params("2D")
	return d["alpha"], d["tx"], d["ty"], int(d["mirror"] + 0.1), d["scale"]


def compose_transform3(
    phi1, theta1, psi1, sx1, sy1, sz1, scale1, phi2, theta2, psi2, sx2, sy2, sz2, scale2
):
    """
	  Compute the composition of two transformations  T2*T1
		Here  if v's are vectors:	vnew = T2*T1 vold
		with T1 described by phi1, sx1,  scale1 etc.

		Usage: compose_transform3(phi1,theta1,psi1,sx1,sy1,sz1,scale1,phi2,theta2,psi2,sx2,sy2,sz2,scale2)
		   angles in degrees
	"""

    R1 = EMAN2_cppwrap.Transform(
        {
            "type": "spider",
            "phi": float(phi1),
            "theta": float(theta1),
            "psi": float(psi1),
            "tx": float(sx1),
            "ty": float(sy1),
            "tz": float(sz1),
            "mirror": 0,
            "scale": float(scale1),
        }
    )
    R2 = EMAN2_cppwrap.Transform(
        {
            "type": "spider",
            "phi": float(phi2),
            "theta": float(theta2),
            "psi": float(psi2),
            "tx": float(sx2),
            "ty": float(sy2),
            "tz": float(sz2),
            "mirror": 0,
            "scale": float(scale2),
        }
    )
    Rcomp = R2 * R1
    d = Rcomp.get_params("spider")
    return d["phi"], d["theta"], d["psi"], d["tx"], d["ty"], d["tz"], d["scale"]


def combine_params2(alpha1, sx1, sy1, mirror1, alpha2, sx2, sy2, mirror2):
    """
	  Combine 2D alignment parameters including mirror: Tnew = T2*T1
	  Combined parameters correspond to image first transformed by set 1 followed by set 2.
	"""

    t1 = EMAN2_cppwrap.Transform(
        {
            "type": "2D",
            "alpha": alpha1,
            "tx": sx1,
            "ty": sy1,
            "mirror": mirror1,
            "scale": 1.0,
        }
    )
    t2 = EMAN2_cppwrap.Transform(
        {
            "type": "2D",
            "alpha": alpha2,
            "tx": sx2,
            "ty": sy2,
            "mirror": mirror2,
            "scale": 1.0,
        }
    )
    tt = t2 * t1
    """Multiline Comment3"""
    d = tt.get_params("2D")
    return d["alpha"], d["tx"], d["ty"], int(d["mirror"] + 0.1)


def inverse_transform2(alpha, tx=0.0, ty=0.0, mirror=0):
    """Returns the inverse of the 2d rot and trans matrix

		Usage: nalpha, ntx, nty, mirror = inverse_transform2(alpha,tx,ty,mirror)
	"""

    t = EMAN2_cppwrap.Transform(
        {
            "type": "2D",
            "alpha": alpha,
            "tx": tx,
            "ty": ty,
            "mirror": mirror,
            "scale": 1.0,
        }
    )
    t = t.inverse()
    t = t.get_params("2D")
    return t["alpha"], t["tx"], t["ty"], int(t["mirror"] + 0.1)


def drop_image(imagename, destination, itype="h"):
    """Write an image to the disk.

	Usage:  drop_image(name_of_existing_image, "path/to/image",
			  type = <type>)
	<type> is "h" (hdf) or "s" (spider)
	"""

    if type(destination) == type(""):
        if itype == "h":
            imgtype = EMAN2_cppwrap.EMUtil.ImageType.IMAGE_HDF
        elif itype == "s":
            imgtype = EMAN2_cppwrap.EMUtil.ImageType.IMAGE_SINGLE_SPIDER
        else:
            sp_global_def.ERROR("unknown image type", "drop_image", 1)
        imagename.write_image(destination, 0, imgtype)
    else:
        sp_global_def.ERROR("destination is not a file name", "drop_image", 1)


def drop_spider_doc(filename, data, comment=None):
    """Create a spider-compatible "Doc" file.

	   filename: name of the Doc file
	   data: List of lists, with the inner list being a list of floats
			 and is written as a line into the doc file.
	"""
    outf = open(filename, "w")

    outf.write(
        " ;   %s   %s   %s\n" % (datetime.datetime.now().ctime(), filename, comment)
    )
    count = 1  # start key from 1; otherwise, it is confusing...
    for dat in data:
        try:
            nvals = len(dat)
            if nvals <= 5:
                datstrings = ["%5d %d" % (count, nvals)]
            else:
                datstrings = ["%6d %d" % (count, nvals)]
            for num in dat:
                datstrings.append("%12.5g" % (num))
        except TypeError:
            # dat is a single number
            datstrings = ["%5d 1%12.5g" % (count, dat)]
        datstrings.append("\n")
        outf.write("".join(datstrings))
        count += 1
    outf.close()


def even_angles(
    delta=15.0,
    theta1=0.0,
    theta2=90.0,
    phi1=0.0,
    phi2=359.99,
    method="S",
    phiEqpsi="Minus",
    symmetry="c1",
    ant=0.0,
):
    """Multiline Comment4"""

    angles = []
    symmetryLower = symmetry.lower()
    symmetry_string = symmetry.split()[0]
    if symmetry_string[0] == "c":
        if phi2 == 359.99:
            angles = even_angles_cd(
                delta,
                theta1,
                theta2,
                phi1 - ant,
                old_div(phi2, int(symmetry_string[1:])) + ant,
                method,
                phiEqpsi,
            )
        else:
            angles = even_angles_cd(
                delta, theta1, theta2, phi1 - ant, phi2 + ant, method, phiEqpsi
            )
        if int(symmetry_string[1:]) > 1:
            if int(symmetry_string[1:]) % 2 == 0:
                qt = old_div(360.0, int(symmetry_string[1:]))
            else:
                qt = old_div(180.0, int(symmetry_string[1:]))
            n = len(angles)
            for i in range(n):
                t = n - i - 1
                if angles[t][1] == 90.0:
                    if angles[t][0] >= qt + ant:
                        del angles[t]
    elif symmetry_string[0] == "d":
        if phi2 == 359.99:
            angles = even_angles_cd(
                delta,
                theta1,
                theta2,
                phi1,
                old_div(360.0, int(symmetry_string[1:])),
                method,
                phiEqpsi,
            )
        else:
            angles = even_angles_cd(delta, theta1, theta2, phi1, phi2, method, phiEqpsi)
        n = len(angles)
        badb = old_div(old_div(360.0, int(symmetry_string[1:])), 4) + ant
        bade = 2 * badb - ant
        bbdb = badb + old_div(old_div(360.0, int(symmetry_string[1:])), 2) + ant
        bbde = bbdb + old_div(old_div(360.0, int(symmetry_string[1:])), 4) - ant
        for i in range(n):
            t = n - i - 1
            qt = angles[t][0]
            if (qt >= badb and qt < bade) or (qt >= bbdb and qt < bbde):
                del angles[t]

        if int(symmetry_string[1:]) % 2 == 0:
            qt = old_div(old_div(360.0, 2), int(symmetry_string[1:]))
        else:
            qt = old_div(old_div(180.0, 2), int(symmetry_string[1:]))
        n = len(angles)
        for i in range(n):
            t = n - i - 1
            if angles[t][1] == 90.0:
                if angles[t][0] >= qt + ant:
                    del angles[t]
    elif symmetry_string[0] == "s":

        # if symetry is "s", deltphi=delta, theata intial=theta1, theta end=90, delttheta=theta2
        # for helical, theta1 cannot be 0.0
        if theta1 > 90.0:
            sp_global_def.ERROR(
                "theta1 must be less than 90.0 for helical symmetry", "even_angles", 1
            )
        if theta1 == 0.0:
            theta1 = 90.0
        theta_number = int(old_div((90.0 - theta1), theta2))
        # for helical, symmetry = s or scn
        cn = int(symmetry_string[2:])
        for j in range(theta_number, -1, -1):

            if j == 0:
                if symmetry_string[1] == "c":
                    if cn % 2 == 0:
                        k = int(old_div(old_div(359.99, cn), delta))
                    else:
                        k = int(old_div(old_div(old_div(359.99, 2), cn), delta))
                elif symmetry_string[1] == "d":
                    if cn % 2 == 0:
                        k = int(old_div(old_div(old_div(359.99, 2), cn), delta))
                    else:
                        k = int(old_div(old_div(old_div(359.99, 4), cn), delta))
                else:
                    sp_global_def.ERROR(
                        "For helical strucutre, we only support scn and sdn symmetry",
                        "even_angles",
                        1,
                    )

            else:
                if symmetry_string[1] == "c":
                    k = int(old_div(old_div(359.99, cn), delta))
                elif symmetry_string[1] == "d":
                    k = int(old_div(old_div(old_div(359.99, 2), cn), delta))

            for i in range(k + 1):
                angles.append([i * delta, 90.0 - j * theta2, 90.0])

    else:
        sp_global_def.ERROR(
            "even_angles", "Symmetry not supported: " + symmetry_string, 0
        )
    """Multiline Comment5"""

    return angles


def center_2D(
    image_to_be_centered,
    center_method=1,
    searching_range=-1,
    Gauss_radius_inner=2,
    Gauss_radius_outter=7,
    self_defined_reference=None,
):
    """Multiline Comment1"""
    # MULTILINEMULTILINEMULTILINE 1
    # MULTILINEMULTILINEMULTILINE 1
    # MULTILINEMULTILINEMULTILINE 1
    # MULTILINEMULTILINEMULTILINE 1
    # MULTILINEMULTILINEMULTILINE 1
    # MULTILINEMULTILINEMULTILINE 1
    # MULTILINEMULTILINEMULTILINE 1
    # MULTILINEMULTILINEMULTILINE 1
    # MULTILINEMULTILINEMULTILINE 1
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_utilities import peak_search
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_fundamentals import fshift
    pass  # IMPORTIMPORTIMPORT import types

    if isinstance(image_to_be_centered, (bytes, str)):
        image_to_be_centered = get_im(image_to_be_centered)
    if center_method == 0:
        return image_to_be_centered, 0.0, 0.0
    elif center_method == 1:
        cs = image_to_be_centered.phase_cog()
        if searching_range > 0:
            if abs(cs[0]) > searching_range:
                cs[0] = 0.0
            if abs(cs[1]) > searching_range:
                cs[1] = 0.0
        return (
            sp_fundamentals.fshift(image_to_be_centered, -cs[0], -cs[1]),
            cs[0],
            cs[1],
        )

    elif center_method == 7:
        pass  # IMPORTIMPORTIMPORT from ..libpy.sp_fundamentals import ccf, cyclic_shift
        pass  # IMPORTIMPORTIMPORT from ..libpy.sp_morphology import binarize
        pass  # IMPORTIMPORTIMPORT from ..libpy.sp_utilities import model_blank
        pass  # IMPORTIMPORTIMPORT from EMAN2 import rsconvolution

        p = EMAN2_cppwrap.Util.infomask(image_to_be_centered, None, True)
        cc = sp_morphology.binarize(
            EMAN2_cppwrap.rsconvolution(
                sp_morphology.binarize(image_to_be_centered, p[0] + p[1]),
                model_blank(5, 5, 1, old_div(1.0, (5.0 * 5.0))),
            ),
            0.5,
        )
        c = sp_fundamentals.ccf(cc, self_defined_reference)
        p = EMAN2_cppwrap.Util.infomask(c, None, True)[3]
        nx = c.get_xsize()
        ny = c.get_ysize()
        cx = old_div(nx, 2)
        cy = old_div(ny, 2)
        n = 0
        x = 0
        y = 0
        for i in range(nx):
            for j in range(ny):
                if c.get_value_at(i, j) == p:
                    x += i - cx
                    y += j - cy
                    n += 1
        shiftx = old_div(x, n)
        shifty = old_div(y, n)
        if searching_range > 0:
            if abs(shiftx) > searching_range:
                shiftx = 0
            if abs(shifty) > searching_range:
                shifty = 0
        return (
            sp_fundamentals.cyclic_shift(image_to_be_centered, -shiftx, -shifty),
            shiftx,
            shifty,
        )

    elif center_method == 5:
        pass  # IMPORTIMPORTIMPORT from ..libpy.sp_fundamentals import rot_avg_image, ccf
        pass  # IMPORTIMPORTIMPORT from math import sqrt

        not_centered = True
        tmp_image = image_to_be_centered.copy()
        shiftx = 0
        shifty = 0
        while not_centered:
            reference = sp_fundamentals.rot_avg_image(tmp_image)
            ccmap = sp_fundamentals.ccf(tmp_image, reference)
            if searching_range > 0:
                ccmap = EMAN2_cppwrap.Util.window(
                    ccmap, searching_range, searching_range, 1, 0, 0, 0
                )
            peak = peak_search(ccmap)
            centered_image = sp_fundamentals.fshift(tmp_image, -peak[0][4], -peak[0][5])
            if numpy.sqrt(peak[0][4] ** 2 + peak[0][5] ** 2) < 1.0:
                not_centered = False
            else:
                tmp_image = centered_image.copy()
            shiftx += peak[0][4]
            shifty += peak[0][5]
        return centered_image, shiftx, shifty

    elif center_method == 6:
        nx = image_to_be_centered.get_xsize()
        ny = image_to_be_centered.get_ysize()
        r = old_div(nx, 2) - 2
        mask = model_circle(r, nx, ny)
        [mean, sigma, xmin, xmax] = EMAN2_cppwrap.Util.infomask(
            image_to_be_centered, mask, True
        )
        new_image = sp_morphology.threshold_to_minval(
            image_to_be_centered, mean + sigma
        )
        cs = new_image.phase_cog()
        if searching_range > 0:
            if abs(cs[0]) > searching_range:
                cs[0] = 0.0
            if abs(cs[1]) > searching_range:
                cs[1] = 0.0
        return (
            sp_fundamentals.fshift(image_to_be_centered, -cs[0], -cs[1]),
            cs[0],
            cs[1],
        )

    else:
        nx = image_to_be_centered.get_xsize()
        ny = image_to_be_centered.get_ysize()

        if center_method == 2:
            reference = model_gauss(Gauss_radius_inner, nx, ny)
        if center_method == 3:
            do1 = model_gauss(Gauss_radius_outter, nx, ny)
            do2 = model_gauss(Gauss_radius_inner, nx, ny)
            s = EMAN2_cppwrap.Util.infomask(do1, None, True)
            do1 = old_div(do1, s[3])
            s = EMAN2_cppwrap.Util.infomask(do2, None, True)
            do2 = old_div(do2, s[3])
            reference = do1 - do2
        if center_method == 4:
            reference = self_defined_reference
        ccmap = sp_fundamentals.ccf(image_to_be_centered, reference)
        if searching_range > 1:
            ccmap = EMAN2_cppwrap.Util.window(
                ccmap, searching_range, searching_range, 1, 0, 0, 0
            )
        peak = peak_search(ccmap)
        return (
            sp_fundamentals.fshift(image_to_be_centered, -peak[0][4], -peak[0][5]),
            peak[0][4],
            peak[0][5],
        )


def even_angles_cd(
    delta, theta1=0.0, theta2=90.0, phi1=0.0, phi2=359.99, method="P", phiEQpsi="Minus"
):
    """Create a list of Euler angles suitable for projections.
	   method is either 'S' - for Saff algorithm
					  or   'P' - for Penczek '94 algorithm
			  'S' assumes phi1<phi2 and phi2-phi1>> delta ;
	   phiEQpsi  - set this to 'Minus', if you want psi=-phi;
	"""

    angles = []
    if method == "P":
        temp = EMAN2_cppwrap.Util.even_angles(delta, theta1, theta2, phi1, phi2)
        # 		                                              phi, theta, psi
        for i in range(old_div(len(temp), 3)):
            angles.append([temp[3 * i], temp[3 * i + 1], temp[3 * i + 2]])
    else:  # elif (method == 'S'):
        Deltaz = numpy.cos(old_div(theta2 * numpy.pi, 180.0)) - numpy.cos(
            old_div(theta1 * numpy.pi, 180.0)
        )
        s = old_div(delta * numpy.pi, 180.0)
        NFactor = old_div(3.6, s)
        wedgeFactor = abs(old_div(Deltaz * (phi2 - phi1), 720.0))
        NumPoints = int(NFactor * NFactor * wedgeFactor)
        angles.append([phi1, theta1, 0.0])
        z1 = numpy.cos(old_div(theta1 * numpy.pi, 180.0))
        phi = phi1  # initialize loop
        for k in range(1, (NumPoints - 1)):
            z = z1 + old_div(Deltaz * k, (NumPoints - 1))
            r = numpy.sqrt(1 - z * z)
            phi = phi1 + (phi + old_div(delta, r) - phi1) % (abs(phi2 - phi1))
            # [k, phi,180*acos(z)/pi, 0]
            angles.append([phi, old_div(180 * math.acos(z), numpy.pi), 0.0])
            # angles.append([p2,t2,0])  # This is incorrect, as the last angle is really the border, not the element we need. PAP 01/15/07
    if phiEQpsi == "Minus":
        for k in range(len(angles)):
            angles[k][2] = (720.0 - angles[k][0]) % 360.0
    if theta2 == 180.0 or (theta2 > 180.0 and delta == 180.0):
        angles.append([0.0, 180.0, 0.0])

    return angles


def find(vv, cmp_str, n):
    jFoundVec = []
    for jFound in range(len(vv)):
        if cmp_str == "lt":
            if vv[jFound] < n:
                jFoundVec.append(jFound)
        if cmp_str == "le":
            if vv[jFound] <= n:
                jFoundVec.append(jFound)
        if cmp_str == "eq":
            if vv[jFound] == n:
                jFoundVec.append(jFound)
        if cmp_str == "ge":
            if vv[jFound] >= n:
                jFoundVec.append(jFound)
        if cmp_str == "gt":
            if vv[jFound] > n:
                jFoundVec.append(jFound)
    return jFoundVec


def gauss_edge(sharp_edge_image, kernel_size=7, gauss_standard_dev=3):
    """
		smooth sharp_edge_image with Gaussian function
		1. The sharp-edge image is convoluted with a gassian kernel
		2. The convolution normalized
	"""

    nz = sharp_edge_image.get_ndim()
    if nz == 3:
        kern = model_gauss(gauss_standard_dev, kernel_size, kernel_size, kernel_size)
    elif nz == 2:
        kern = model_gauss(gauss_standard_dev, kernel_size, kernel_size)
    else:
        kern = model_gauss(gauss_standard_dev, kernel_size)
    aves = EMAN2_cppwrap.Util.infomask(kern, None, False)
    nx = kern.get_xsize()
    ny = kern.get_ysize()
    nz = kern.get_zsize()

    kern = old_div(kern, aves[0] * nx * ny * nz)
    return EMAN2_cppwrap.rsconvolution(sharp_edge_image, kern)


def get_image(imagename, nx=0, ny=1, nz=1, im=0):
    """Read an image from the disk or assign existing object to the output.

	Usage: myimage = readImage("path/to/image")
	or     myimage = readImage(name_of_existing_image)
	"""
    if type(imagename) == type(""):
        e = EMAN2_cppwrap.EMData()
        e.read_image(imagename, im)
    elif not imagename:
        e = EMAN2_cppwrap.EMData()
        if nx > 0:
            e.set_size(nx, ny, nz)
    else:
        e = imagename
    return e


def get_im(stackname, im=0):
    """Read an image from the disk stack, or return im's image from the list of images

	Usage: myimage = get_im("path/to/stack", im)
	   or: myimage = get_im( data, im )
	"""
    if type(stackname) == type(""):
        e = EMAN2_cppwrap.EMData()
        e.read_image(stackname, im)
        return e
    else:
        return stackname[im].copy()


def get_image_data(img):
    """
		Return a NumPy array containing the image data.
		Note: The NumPy array and the image data share the same memory,
		so if the NumPy array is altered then the image is altered
		as well (and vice versa).
	"""

    return EMAN2_cppwrap.EMNumPy.em2numpy(img)


def get_symt(symmetry):
    """
	get a list of point-group symmetry transformations, symmetry="c3"
	"""

    scl = sp_fundamentals.symclass(symmetry)
    trans = []
    for q in scl.symangles:
        trans.append(
            EMAN2_cppwrap.Transform(
                {"type": "spider", "phi": q[0], "theta": q[1], "psi": q[2]}
            )
        )
    return trans


def get_input_from_string(str_input):
    """
		Extract input numbers from a given string
	"""

    qq = re.split(" |,", str_input)
    for i in range(len(qq) - 1, -1, -1):
        if qq[i] == "":
            del qq[i]
    o = []
    for i in range(len(qq)):
        if qq[i].find(".") >= 0:
            o.append(float(qq[i]))
        else:
            o.append(int(qq[i]))
    return o


def model_circle(r, nx, ny, nz=1):
    """
	Create a centered circle (or sphere) having radius r.
	"""
    e = EMAN2_cppwrap.EMData()
    e.set_size(nx, ny, nz)
    e.process_inplace("testimage.circlesphere", {"radius": r, "fill": 1})
    return e


def model_gauss(
    xsigma,
    nx,
    ny=1,
    nz=1,
    ysigma=None,
    zsigma=None,
    xcenter=None,
    ycenter=None,
    zcenter=None,
):
    """Multiline Comment6"""
    e = EMAN2_cppwrap.EMData()
    e.set_size(nx, ny, nz)
    if ysigma == None:
        ysigma = xsigma
    if zsigma == None:
        zsigma = xsigma
    if xcenter == None:
        xcenter = old_div(nx, 2)
    if ycenter == None:
        ycenter = old_div(ny, 2)
    if zcenter == None:
        zcenter = old_div(nz, 2)
    e.process_inplace(
        "testimage.puregaussian",
        {
            "x_sigma": xsigma,
            "y_sigma": ysigma,
            "z_sigma": zsigma,
            "x_center": xcenter,
            "y_center": ycenter,
            "z_center": zcenter,
        },
    )
    return e


def model_cylinder(radius, nx, ny, nz):
    """
	 create a cylinder along z axis
	"""
    e = EMAN2_cppwrap.EMData()
    e.set_size(nx, ny, nz)
    e.process_inplace("testimage.cylinder", {"radius": radius})
    return e


def model_rotated_rectangle2D(
    radius_long, radius_short, nx, ny, angle=90, return_numpy=False
):
    """
	Creates a rectangular mask
	:param radius_long: Radius long axis
	:param radius_short: Radius short axis
	:param nx: Mask size x-dim
	:param ny: Mask size y-dim
	:param angle: Rotation angle in degree
	:param return_numpy: Return mask as numpy array instead of as EMObject. [Default: False]
	:return: rotated rectangular mask
	"""

    sizex = nx
    sizey = ny
    if (radius_long * 2) > nx or (radius_long * 2) > ny:
        sizex = radius_long * 2
        sizey = radius_long * 2

    mask = numpy.zeros((sizex, sizey))
    mask[
        (old_div(sizex, 2) - radius_short) : (old_div(sizex, 2) + radius_short),
        (old_div(sizey, 2) - radius_long) : (old_div(sizey, 2) + radius_long),
    ] = 1.0

    mask = scipy.ndimage.rotate(mask, angle, reshape=False)
    offx = 0 if nx % 2 == 0 else 1  # add one column if nx is odd
    offy = 0 if nx % 2 == 0 else 1  # add one row if ny is odd
    mask = mask[  # otherwise this only creates even-sized masks
        (old_div(sizex, 2) - old_div(nx, 2)) : (old_div(sizex, 2) + old_div(nx, 2))
        + offy,
        (old_div(sizey, 2) - old_div(ny, 2)) : (old_div(sizey, 2) + old_div(ny, 2))
        + offx,
    ]

    # eliminate round-off errors introduced by the rotation
    mask[mask <= 0.9] = 0.0
    mask[mask > 0.9] = 1.0

    # deliver the mask
    if return_numpy:
        return mask
    else:
        return_mask = numpy2em_python(mask)
        return return_mask


def model_gauss_noise(sigma, nx, ny=1, nz=1):
    """
	Create an image of noise having standard deviation "sigma",
	and average 0.
	"""
    e = EMAN2_cppwrap.EMData()
    e.set_size(nx, ny, nz)
    if sigma == 0.0:
        e.to_zero()
    else:
        e.process_inplace("testimage.noise.gauss", {"sigma": sigma})
    return e


def model_blank(nx, ny=1, nz=1, bckg=0.0):
    """
	Create a blank image.
	"""
    e = EMAN2_cppwrap.EMData()
    e.set_size(nx, ny, nz)
    e.to_zero()
    if bckg != 0.0:
        e += bckg
    return e


def peak_search(e, npeak=1, invert=1, print_screen=0):
    peaks = e.peak_search(npeak, invert)
    ndim = peaks[0]
    nlist = int(old_div((len(peaks) - 1), ((ndim + 1) * 2)))
    if nlist > 0:
        outpeaks = []
        if print_screen:
            if ndim == 1:
                sp_global_def.sxprint(
                    "%10s%10s%10s%10s%10s"
                    % ("Index  ", " Peak_value", "X   ", "Peak/P_max", "X-NX/2")
                )
                print_list_format(peaks[1:], 4)
            elif ndim == 2:
                sp_global_def.sxprint(
                    "%10s%10s%10s%10s%10s%10s%10s"
                    % (
                        "Index  ",
                        "Peak_value",
                        "X   ",
                        "Y   ",
                        "Peak/P_max",
                        "X-NX/2",
                        "Y-NY/2",
                    )
                )
                print_list_format(peaks[1:], 6)
            elif ndim == 3:
                sp_global_def.sxprint(
                    "%10s%10s%10s%10s%10s%10s%10s%10s%10s"
                    % (
                        "Index  ",
                        "Peak_value",
                        "X   ",
                        "Y   ",
                        "Z   ",
                        "Peak/P_max",
                        "X-NX/2",
                        "Y-NY/2",
                        "Z-NZ/2",
                    )
                )
                print_list_format(peaks[1:], 8)
            else:
                sp_global_def.ERROR(
                    "Image dimension extracted in peak_search is wrong",
                    "Util.peak_search",
                    1,
                )
        for i in range(nlist):
            k = int((ndim + 1) * i * 2)
            if ndim == 1:
                p = [peaks[k + 1], peaks[k + 2], peaks[k + 3], peaks[k + 4]]
            elif ndim == 2:
                p = [
                    peaks[k + 1],
                    peaks[k + 2],
                    peaks[k + 3],
                    peaks[k + 4],
                    peaks[k + 5],
                    peaks[k + 6],
                ]
            elif ndim == 3:
                p = [
                    peaks[k + 1],
                    peaks[k + 2],
                    peaks[k + 3],
                    peaks[k + 4],
                    peaks[k + 5],
                    peaks[k + 6],
                    peaks[k + 7],
                    peaks[k + 8],
                ]
            outpeaks.append(p)
    else:
        ndim = e.get_ndim()
        # ERROR("peak search fails to find any peaks, returns image center as a default peak position","peak_search",0)
        if ndim == 1:
            nx = e.get_xsize()
            outpeaks = [[1.0, float(old_div(nx, 2)), 1.0, 0.0]]
        elif ndim == 2:
            nx = e.get_xsize()
            ny = e.get_ysize()
            outpeaks = [
                [1.0, float(old_div(nx, 2)), float(old_div(ny, 2)), 1.0, 0.0, 0.0]
            ]
        elif ndim == 3:
            nx = e.get_xsize()
            ny = e.get_ysize()
            nz = e.get_ysize()
            outpeaks = [
                [
                    1.0,
                    float(old_div(nx, 2)),
                    float(old_div(ny, 2)),
                    float(old_div(nz, 2)),
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                ]
            ]
    return outpeaks


####--------------------------------------------------------------------------------------------------#########


def print_list_format(m, narray=0):

    """Multiline Comment7"""
    flist = []
    for i in range(len(m)):
        if type(m[i]) is float:
            flist.append("%10.3g" % (m[i]))
        elif type(m[i]) is int:
            flist.append("%10d" % (m[i]))
        else:
            flist.append("%10s" % (m[i]))
    if narray > len(m):
        narray = 0
        sp_global_def.ERROR(
            "improper input narray number, use default value", "print_list_foramt", 0
        )
    if narray == 0:
        num = int(numpy.sqrt(len(m)))
        if len(m) % num != 0:
            lnum = int(old_div(len(m), num)) + 1
        else:
            lnum = int(old_div(len(m), num))
    else:
        num = narray
        if len(m) % num == 0:
            lnum = int(old_div(len(m), num))
        else:
            lnum = int(old_div(len(m), num)) + 1
    ncount = -1
    plist = []
    for i in range(lnum):
        qlist = ""
        for j in range(num):
            ncount += 1
            if ncount <= len(m) - 1:
                qlist = qlist + flist[ncount]
            else:
                break
        plist.append(qlist)
    for i in range(lnum):
        sp_global_def.sxprint("%6d " % (i + 1), plist[i])


def pad(
    image_to_be_padded,
    new_nx,
    new_ny=1,
    new_nz=1,
    background="average",
    off_center_nx=0,
    off_center_ny=0,
    off_center_nz=0,
):

    if not isinstance(background, (bytes, str)):
        background = str(background)
    if background == "average":
        image_padded = EMAN2_cppwrap.Util.pad(
            image_to_be_padded,
            new_nx,
            new_ny,
            new_nz,
            off_center_nx,
            off_center_ny,
            off_center_nz,
            "average",
        )
    elif background == "circumference":
        image_padded = EMAN2_cppwrap.Util.pad(
            image_to_be_padded,
            new_nx,
            new_ny,
            new_nz,
            off_center_nx,
            off_center_ny,
            off_center_nz,
            "circumference",
        )
    else:
        image_padded = EMAN2_cppwrap.Util.pad(
            image_to_be_padded,
            new_nx,
            new_ny,
            new_nz,
            off_center_nx,
            off_center_ny,
            off_center_nz,
            background,
        )
    return image_padded


def chooseformat(t, form_float="  %12.5f"):

    e_form = form_float.replace("f", "e")
    ee = form_float.strip() % t
    if len(ee) > int(form_float.strip().split(".")[0][1:]):
        return e_form
    df1 = float(ee)
    df2 = float(e_form.strip() % t)
    if abs(t - df1) <= abs(t - df2):
        return form_float
    else:
        return e_form


def read_text_row(fnam, format="", skip=";"):
    """
		Read a column-listed txt file.
		INPUT: filename: name of the Doc file
		OUTPUT:
			nc : number of entries in each lines (number of columns)
			len(data)/nc : number of lines (rows)
			data: List of numbers from the doc file
	"""

    inf = open(fnam, "r")
    strg = inf.readline()
    x = []
    data = []
    while len(strg) > 0:
        com_line = False
        for j in range(len(strg)):
            if strg[j] == skip:
                com_line = True
        if com_line == False:
            word = strg.split()
            if format == "s":
                key = int(word[1])
                if key != len(word) - 2:
                    del word
                    word = []
                    word.append(strg[0:5])
                    word.append(strg[6:7])
                    for k in range(key):
                        k_start = 7 + k * 13
                        k_stop = k_start + 13
                        word.append(strg[k_start:k_stop])
            line = []
            for i in range(len(word)):
                try:
                    line.append(int(word[i]))
                except:
                    try:
                        line.append(float(word[i]))
                    except:
                        line.append(word[i])
            data.append(line)
        strg = inf.readline()
    inf.close
    return data


def write_text_row(data, file_name, form_float="  %14.6f", form_int="  %12d"):
    """
	   Write to an ASCII file a list of lists containing floats.

	   filename: name of the text file
	   data: List of lists, with the inner list being a list of floats, i.e., [ [first list], [second list], ...]
			 First list will be written as a first line, second as a second, and so on...
		 If only one list is given, the file will contain one line
	"""
    outf = open(file_name, "w")
    if type(data[0]) == list:
        # It is a list of lists
        for i in range(len(data)):
            for j in range(len(data[i])):
                tpt = data[i][j]
                qtp = type(tpt)
                if qtp == int:
                    outf.write(form_int % tpt)
                elif qtp == float:
                    frmt = chooseformat(tpt, form_float)
                    # if find(frmt, "e") < 0:
                    #     outf.write(frmt % tpt)
                    # else:
                    outf.write(frmt % tpt)
                else:
                    outf.write("  %s" % tpt)
            outf.write("\n")
    else:
        # Single list
        for j in range(len(data)):
            tpt = data[j]
            qtp = type(tpt)
            if qtp == int:
                outf.write(form_int % tpt + "\n")
            elif qtp == float:
                frmt = chooseformat(tpt, form_float)
                # if find(frmt, "e") < 0:
                #     outf.write(frmt % tpt + "\n")
                # else:
                outf.write(frmt % tpt + "\n")
            else:
                outf.write("  %s\n" % tpt)
    outf.flush()
    outf.close()


def read_text_file(file_name, ncol=0):
    """
		Read data from text file, if ncol = -1, read all columns
		if ncol >= 0, just read the (ncol)-th column.
	"""

    inf = open(file_name, "r")
    line = inf.readline()
    data = []
    while len(line) > 0:
        if ncol == -1:
            vdata = line.split()
            if data == []:
                for i in range(len(vdata)):
                    try:
                        data.append([int(vdata[i])])
                    except:
                        try:
                            data.append([float(vdata[i])])
                        except:
                            data.append([vdata[i]])
            else:
                for i in range(len(vdata)):
                    try:
                        data[i].append(int(vdata[i]))
                    except:
                        try:
                            data[i].append(float(vdata[i]))
                        except:
                            data[i].append(vdata[i])
        else:
            vdata = line.split()[ncol]
            try:
                data.append(int(vdata))
            except:
                try:
                    data.append(float(vdata))
                except:
                    data.append(vdata)
        line = inf.readline()
    return data


def write_text_file(data, file_name, form_float="  %14.6f", form_int="  %12d"):
    """
	   Write to an ASCII file a list of lists containing floats.

	   filename: name of the text file
	   data: List of lists, with the inner list being a list of floats, i.e., [ [first list], [second list], ...]
			 First list will be written as a first column, second as a second, and so on...
		 If only one list is given, the file will contain one column
	"""

    if data == []:
        outf = open(file_name, "w")
        outf.close()
        return

    outf = open(file_name, "w")
    if type(data[0]) == list:
        # It is a list of lists
        for i in range(len(data[0])):
            for j in range(len(data)):
                tpt = data[j][i]
                qtp = type(tpt)
                if qtp == int:
                    outf.write(form_int % tpt)
                elif qtp == float:
                    frmt = chooseformat(tpt, form_float)
                    # if find(frmt, "e") < 0:
                    #     outf.write(frmt % tpt)
                    # else:
                    outf.write(frmt % tpt)
                else:
                    outf.write("  %s" % tpt)
            outf.write("\n")
    else:
        # Single list
        for j in range(len(data)):
            tpt = data[j]
            qtp = type(tpt)
            if qtp == int:
                outf.write(form_int % tpt + "\n")
            elif qtp == float:
                frmt = chooseformat(tpt, form_float)
                # if find(frmt, "e") < 0:
                #     outf.write(frmt % tpt + "\n")
                # else:
                outf.write(frmt % tpt + "\n")
            else:
                outf.write("  %s\n" % tpt)
    outf.close()


def rotate_shift_params(paramsin, transf):
    # moved from sxprocess.py
    if len(paramsin[0]) > 3:

        t = EMAN2_cppwrap.Transform(
            {
                "type": "spider",
                "phi": transf[0],
                "theta": transf[1],
                "psi": transf[2],
                "tx": transf[3],
                "ty": transf[4],
                "tz": transf[5],
                "mirror": 0,
                "scale": 1.0,
            }
        )
        t = t.inverse()
        cpar = []
        for params in paramsin:
            d = EMAN2_cppwrap.Transform(
                {
                    "type": "spider",
                    "phi": params[0],
                    "theta": params[1],
                    "psi": params[2],
                }
            )
            d.set_trans(EMAN2_cppwrap.Vec2f(-params[3], -params[4]))
            c = d * t
            u = c.get_params("spider")
            cpar.append([u["phi"], u["theta"], u["psi"], -u["tx"], -u["ty"]])
    else:
        t = EMAN2_cppwrap.Transform(
            {"type": "spider", "phi": transf[0], "theta": transf[1], "psi": transf[2]}
        )
        t = t.inverse()
        cpar = []
        for params in paramsin:
            d = EMAN2_cppwrap.Transform(
                {
                    "type": "spider",
                    "phi": params[0],
                    "theta": params[1],
                    "psi": params[2],
                }
            )
            c = d * t
            u = c.get_params("spider")
            cpar.append([u["phi"], u["theta"], u["psi"]])
            # cpar.append([u["phi"],u["theta"],u["psi"],-u["tx"],-u["ty"]])
    return cpar


def reshape_1d(
    input_object,
    length_current=0,
    length_interpolated=0,
    Pixel_size_current=0.0,
    Pixel_size_interpolated=0.0,
):
    """Multiline Comment9"""

    interpolated = []
    if length_current == 0:
        length_current = len(input_object)
    lt = len(input_object) - 2
    if length_interpolated == 0:
        if Pixel_size_interpolated != Pixel_size_current:
            length_interpolated = int(
                old_div(length_current * Pixel_size_current, Pixel_size_interpolated)
                + 0.5
            )
        else:
            sp_global_def.ERROR("Incorrect input parameters", "reshape_1d", 1)
            return []

    if Pixel_size_current == 0.0:
        Pixel_size_current = 1.0
        Pixel_size_interpolated = old_div(
            Pixel_size_current * float(length_current), float(length_interpolated)
        )
    qt = old_div(Pixel_size_interpolated, Pixel_size_current)

    for i in range(length_interpolated):
        xi = float(i) * qt
        ix = min(int(xi), lt)
        df = xi - ix
        xval = input_object[ix] + df * (input_object[ix + 1] - input_object[ix])
        interpolated.append(xval)
    return interpolated


def estimate_3D_center_MPI(data, nima, myid, number_of_proc, main_node, mpi_comm=None):

    if mpi_comm == None:
        mpi_comm = mpi.MPI_COMM_WORLD

    ali_params_series = []
    for im in data:
        phi, theta, psi, s2x, s2y = get_params_proj(im)
        ali_params_series.append(phi)
        ali_params_series.append(theta)
        ali_params_series.append(psi)
        ali_params_series.append(s2x)
        ali_params_series.append(s2y)

    if myid == main_node:
        for proc in range(number_of_proc):
            if proc != main_node:
                image_start_proc, image_end_proc = sp_applications.MPI_start_end(
                    nima, number_of_proc, proc
                )
                n_params = (image_end_proc - image_start_proc) * 5
                temp = mpi.mpi_recv(n_params, mpi.MPI_FLOAT, proc, proc, mpi_comm)
                for nn in range(n_params):
                    ali_params_series.append(float(temp[nn]))

        ali_params = []
        N = old_div(len(ali_params_series), 5)
        for im in range(N):
            ali_params.append(
                [
                    ali_params_series[im * 5],
                    ali_params_series[im * 5 + 1],
                    ali_params_series[im * 5 + 2],
                    ali_params_series[im * 5 + 3],
                    ali_params_series[im * 5 + 4],
                ]
            )

        A = []
        b = []

        for i in range(N):
            phi_rad = numpy.radians(ali_params[i][0])
            theta_rad = numpy.radians(ali_params[i][1])
            psi_rad = numpy.radians(ali_params[i][2])
            A.append(
                [
                    numpy.cos(psi_rad) * numpy.cos(theta_rad) * numpy.cos(phi_rad)
                    - numpy.sin(psi_rad) * numpy.sin(phi_rad),
                    numpy.cos(psi_rad) * numpy.cos(theta_rad) * numpy.sin(phi_rad)
                    + numpy.sin(psi_rad) * numpy.cos(phi_rad),
                    -numpy.cos(psi_rad) * numpy.sin(theta_rad),
                    1,
                    0,
                ]
            )
            A.append(
                [
                    -numpy.sin(psi_rad) * numpy.cos(theta_rad) * numpy.cos(phi_rad)
                    - numpy.cos(psi_rad) * numpy.sin(phi_rad),
                    -numpy.sin(psi_rad) * numpy.cos(theta_rad) * numpy.sin(phi_rad)
                    + numpy.cos(psi_rad) * numpy.cos(phi_rad),
                    numpy.sin(psi_rad) * numpy.sin(theta_rad),
                    0,
                    1,
                ]
            )
            b.append([ali_params[i][3]])
            b.append([ali_params[i][4]])

        A_matrix = numpy.matrix(A)
        b_matrix = numpy.matrix(b)

        K = numpy.linalg.solve(A_matrix.T * A_matrix, A_matrix.T * b_matrix)
        return (
            float(K[0][0]),
            float(K[1][0]),
            float(K[2][0]),
            float(K[3][0]),
            float(K[4][0]),
        )

    else:
        image_start_proc, image_end_proc = sp_applications.MPI_start_end(
            nima, number_of_proc, myid
        )
        n_params = (image_end_proc - image_start_proc) * 5
        mpi.mpi_send(
            ali_params_series, n_params, mpi.MPI_FLOAT, main_node, myid, mpi_comm
        )

        return 0.0, 0.0, 0.0, 0.0, 0.0


def rotate_3D_shift(data, shift3d):

    t = EMAN2_cppwrap.Transform(
        {
            "type": "spider",
            "phi": 0.0,
            "theta": 0.0,
            "psi": 0.0,
            "tx": -shift3d[0],
            "ty": -shift3d[1],
            "tz": -shift3d[2],
            "mirror": 0,
            "scale": 1.0,
        }
    )

    for i in range(len(data)):
        d = data[i].get_attr("xform.projection")
        c = d * t
        data[i].set_attr("xform.projection", c)


def set_arb_params(img, params, par_str):

    """
		filling arbitary headers
	"""
    for i in range(len(par_str)):
        img.set_attr_dict({par_str[i]: params[i]})


def get_arb_params(img, par_str):

    """
		reading arbitary headers
	"""
    params = []
    for i in range(len(par_str)):
        params.append(img.get_attr(par_str[i]))
    return params


###------------------------------------------------------------------------------------------


"""Multiline Comment10"""


def reduce_EMData_to_root(data, myid, main_node=0, comm=-1):

    if comm == -1 or comm == None:
        comm = mpi.MPI_COMM_WORLD

    array = get_image_data(data)
    n = numpy.shape(array)
    ntot = 1
    for i in n:
        ntot *= i
    count = (75 * 4 + 2) * (75 * 4) ** 2
    array1d = numpy.reshape(array, (ntot,))
    ntime = old_div((ntot - 1), count) + 1
    for i in range(ntime):
        block_begin = i * count
        block_end = min(block_begin + count, ntot)
        block_size = block_end - block_begin
        tmpsum = mpi.mpi_reduce(
            array1d[block_begin : block_begin + block_size],
            block_size,
            mpi.MPI_FLOAT,
            mpi.MPI_SUM,
            main_node,
            comm,
        )
        mpi.mpi_barrier(comm)
        if myid == main_node:
            array1d[block_begin:block_end] = tmpsum[0:block_size]


def bcast_compacted_EMData_all_to_all(list_of_em_objects, myid, comm=-1):

    """
	The assumption in <<bcast_compacted_EMData_all_to_all>> is that each processor
	calculates part of the list of elements and then each processor sends
	its results to the other ones.

	Therefore, each processor has access to the header. If we assume that the
	attributes of interest from the header are the same for all elements then
	we can copy the header and no mpi message is necessary for the
	header.

	"""

    if comm == -1 or comm == None:
        comm = mpi.MPI_COMM_WORLD

    num_ref = len(list_of_em_objects)
    ncpu = mpi.mpi_comm_size(comm)  # Total number of processes, passed by --np option.

    ref_start, ref_end = sp_applications.MPI_start_end(num_ref, ncpu, myid)

    for first_myid_process_that_has_em_elements in range(ncpu):
        sim_start, sim_ref_end = sp_applications.MPI_start_end(
            num_ref, ncpu, first_myid_process_that_has_em_elements
        )
        if sim_start != sim_ref_end:
            break
    else:
        raise ValueError("No processor contains em objects!")

    if myid == first_myid_process_that_has_em_elements:
        # used for copying the header and other info

        reference_em_object = list_of_em_objects[ref_start]
        data = EMAN2_cppwrap.EMNumPy.em2numpy(reference_em_object)
        size_of_one_refring_assumed_common_to_all = data.size

        nx = reference_em_object.get_xsize()
        ny = reference_em_object.get_ysize()
        nz = reference_em_object.get_zsize()

        em_dict = reference_em_object.get_attr_dict()
        dict_to_send = {
            "size_of_one_refring_assumed_common_to_all": size_of_one_refring_assumed_common_to_all,
            "em_dict": em_dict,
            "nx": nx,
            "ny": ny,
            "nz": nz,
        }
    else:
        dict_to_send = None

    dict_received = wrap_mpi_bcast(
        dict_to_send, first_myid_process_that_has_em_elements, comm
    )

    em_dict = dict_received["em_dict"]
    nx = dict_received["nx"]
    ny = dict_received["ny"]
    nz = dict_received["nz"]
    size_of_one_refring_assumed_common_to_all = dict_received[
        "size_of_one_refring_assumed_common_to_all"
    ]

    if size_of_one_refring_assumed_common_to_all * (ref_end - ref_start) > (
        2 ** 31 - 1
    ):
        sp_global_def.sxprint(
            "Sending refrings: size of data to broadcast is greater than 2GB"
        )

    for sender_id in range(ncpu):
        sender_ref_start, sender_ref_end = sp_applications.MPI_start_end(
            num_ref, ncpu, sender_id
        )

        if sender_id == myid:
            if ref_start == ref_end:
                continue
            data = EMAN2_cppwrap.EMNumPy.em2numpy(
                list_of_em_objects[ref_start]
            )  # array([], dtype = 'float32')
            for i in range(ref_start + 1, ref_end):
                data = numpy.concatenate(
                    [data, EMAN2_cppwrap.EMNumPy.em2numpy(list_of_em_objects[i])]
                )
        else:
            if sender_ref_start == sender_ref_end:
                continue
            data = numpy.array([], dtype="float32")

        sender_size_of_refrings = (
            sender_ref_end - sender_ref_start
        ) * size_of_one_refring_assumed_common_to_all

        # size_of_refrings = mpi_bcast(size_of_refrings, 1, MPI_INT, sender_id, comm)
        data = mpi.mpi_bcast(
            data, sender_size_of_refrings, mpi.MPI_FLOAT, sender_id, comm
        )
        # print "Just sent %d float32 elements"%data.size

        if myid != sender_id:
            for i in range(sender_ref_start, sender_ref_end):
                offset_ring = sender_ref_start
                start_p = (i - offset_ring) * size_of_one_refring_assumed_common_to_all
                end_p = (
                    i + 1 - offset_ring
                ) * size_of_one_refring_assumed_common_to_all
                image_data = data[start_p:end_p]

                if int(nz) != 1:
                    image_data = numpy.reshape(image_data, (nz, ny, nx))
                elif ny != 1:
                    image_data = numpy.reshape(image_data, (ny, nx))

                em_object = numpy2em_python(image_data)
                em_object.set_attr_dict(em_dict)
                list_of_em_objects[i] = em_object


def bcast_EMData_to_all(tavg, myid, source_node=0, comm=-1):

    if comm == -1 or comm == None:
        comm = mpi.MPI_COMM_WORLD
    tavg_data = EMAN2_cppwrap.EMNumPy.em2numpy(tavg)
    n = numpy.shape(tavg_data)
    ntot = 1
    for i in n:
        ntot *= i
    tavg_tmp = mpi.mpi_bcast(tavg_data, ntot, mpi.MPI_FLOAT, source_node, comm)
    if myid != source_node:
        tavg_data1d = numpy.reshape(tavg_data, (ntot,))
        tavg_data1d[0:ntot] = tavg_tmp[0:ntot]


"""Multiline Comment12"""


def send_EMData(img, dst, tag, comm=-1):

    if comm == -1:
        comm = mpi.MPI_COMM_WORLD
    img_head = []
    img_head.append(img.get_xsize())
    img_head.append(img.get_ysize())
    img_head.append(img.get_zsize())
    img_head.append(img.is_complex())
    img_head.append(img.is_ri())
    img_head.append(img.get_attr("changecount"))
    img_head.append(img.is_complex_x())
    img_head.append(img.get_attr("is_complex_ri"))
    img_head.append(int(img.get_attr("apix_x") * 10000))
    img_head.append(int(img.get_attr("apix_y") * 10000))
    img_head.append(int(img.get_attr("apix_z") * 10000))

    head_tag = 2 * tag
    mpi.mpi_send(img_head, 11, mpi.MPI_INT, dst, head_tag, comm)

    img_data = get_image_data(img)
    data_tag = 2 * tag + 1
    ntot = img_head[0] * img_head[1] * img_head[2]
    mpi.mpi_send(img_data, ntot, mpi.MPI_FLOAT, dst, data_tag, comm)

    """Multiline Comment13"""


def recv_EMData(src, tag, comm=-1):

    if comm == -1:
        comm = mpi.MPI_COMM_WORLD
    head_tag = 2 * tag
    img_head = mpi.mpi_recv(11, mpi.MPI_INT, src, head_tag, comm)

    nx = int(img_head[0])
    ny = int(img_head[1])
    nz = int(img_head[2])
    is_complex = int(img_head[3])
    is_ri = int(img_head[4])

    data_tag = 2 * tag + 1
    ntot = nx * ny * nz

    img_data = mpi.mpi_recv(ntot, mpi.MPI_FLOAT, src, data_tag, comm)
    if nz != 1:
        img_data = numpy.reshape(img_data, (nz, ny, nx))
    elif ny != 1:
        img_data = numpy.reshape(img_data, (ny, nx))
    else:
        pass

    img = numpy2em_python(img_data)
    img.set_complex(is_complex)
    img.set_ri(is_ri)
    img.set_attr_dict(
        {
            "changecount": int(img_head[5]),
            "is_complex_x": int(img_head[6]),
            "is_complex_ri": int(img_head[7]),
            "apix_x": old_div(int(img_head[8]), 10000.0),
            "apix_y": old_div(int(img_head[9]), 10000.0),
            "apix_z": old_div(int(img_head[10]), 10000.0),
        }
    )
    return img

    """Multiline Comment14"""


def send_string_to_all(str_to_send, source_node=0):
    myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)

    str_to_send_len = len(str_to_send) * int(myid == source_node)
    str_to_send_len = mpi.mpi_bcast(
        str_to_send_len, 1, mpi.MPI_INT, source_node, mpi.MPI_COMM_WORLD
    )[0]
    str_to_send = mpi.mpi_bcast(
        str_to_send, str_to_send_len, mpi.MPI_CHAR, source_node, mpi.MPI_COMM_WORLD
    )

    return b"".join(str_to_send).decode('latin1')


def bcast_number_to_all(number_to_send, source_node=0, mpi_comm=-1):
    """
		number_to_send has to be pre-defined in each node
	"""

    if mpi_comm == -1:
        mpi_comm = mpi.MPI_COMM_WORLD
    if type(number_to_send) is int:
        TMP = mpi.mpi_bcast(number_to_send, 1, mpi.MPI_INT, source_node, mpi_comm)
        return int(TMP[0])
    elif type(number_to_send) is float:
        TMP = mpi.mpi_bcast(number_to_send, 1, mpi.MPI_FLOAT, source_node, mpi_comm)
        return float(TMP[0])
    elif type(number_to_send) is bool:
        if number_to_send:
            number_to_send = 1
        else:
            number_to_send = 0
        TMP = mpi.mpi_bcast(number_to_send, 1, mpi.MPI_INT, source_node, mpi_comm)
        if TMP == 1:
            return True
        else:
            return False
    else:
        sp_global_def.sxprint(" ERROR in bcast_number_to_all")


def bcast_list_to_all(list_to_send, myid, source_node=0, mpi_comm=-1):

    if mpi_comm == -1:
        mpi_comm = mpi.MPI_COMM_WORLD
    if myid == source_node:
        n = len(list_to_send)
        # we will also assume all elements on the list are of the same type
        if type(list_to_send[0]) == int:
            tp = 0
        elif type(list_to_send[0]) == float:
            tp = 1
        else:
            tp = 2
    else:
        n = 0
        tp = 0
    n = bcast_number_to_all(n, source_node=source_node, mpi_comm=mpi_comm)
    tp = bcast_number_to_all(tp, source_node=source_node, mpi_comm=mpi_comm)
    if tp == 2:
        sp_global_def.ERROR(
            "Only list of the same type numbers can be brodcasted",
            "bcast_list_to_all",
            1,
            myid,
        )
    if myid != source_node:
        list_to_send = [0] * n

    if tp == 0:
        list_to_send = mpi.mpi_bcast(
            list_to_send, n, mpi.MPI_INT, source_node, mpi_comm
        )
        return [int(n) for n in list_to_send]
    else:
        list_to_send = mpi.mpi_bcast(
            list_to_send, n, mpi.MPI_FLOAT, source_node, mpi_comm
        )
        return [float(n) for n in list_to_send]


def recv_attr_dict(
    main_node, stack, data, list_params, image_start, image_end, number_of_proc, comm=-1
):

    #   hdf version!
    # This is done on the main node, so for images from the main node, simply write headers

    if comm == -1:
        comm = mpi.MPI_COMM_WORLD

    TransType = type(EMAN2_cppwrap.Transform())
    # prepare keys for float/int
    value = get_arb_params(data[0], list_params)
    ink = []
    len_list = 0
    for il in range(len(list_params)):
        if type(value[il]) is int:
            ink.append(1)
            len_list += 1
        elif type(value[il]) is float:
            ink.append(0)
            len_list += 1
        elif type(value[il]) is TransType:
            ink.append(2)
            len_list += 12
    ldis = []
    headers = []
    for n in range(number_of_proc):
        if n != main_node:
            dis = mpi.mpi_recv(
                2, mpi.MPI_INT, n, sp_global_def.SPARX_MPI_TAG_UNIVERSAL, comm
            )
            value = mpi.mpi_recv(
                len_list * (dis[1] - dis[0]),
                mpi.MPI_FLOAT,
                n,
                sp_global_def.SPARX_MPI_TAG_UNIVERSAL,
                comm,
            )
            ldis.append([dis[0], dis[1]])
            headers.append(value)
            del dis
    del value
    for im in range(image_start, image_end):
        data[im - image_start].write_image(
            stack,
            data[im - image_start].get_attr_default("ID", im),
            EMAN2_cppwrap.EMUtil.ImageType.IMAGE_HDF,
            True,
        )

    for n in range(len(ldis)):
        img_begin = ldis[n][0]
        img_end = ldis[n][1]
        for im in range(img_begin, img_end):
            par_begin = (im - img_begin) * len_list
            nvalue = []
            header = headers[n]
            ilis = 0
            for il in range(len(list_params)):
                if ink[il] == 1:
                    nvalue.append(int(header[par_begin + ilis]))
                    ilis += 1
                elif ink[il] == 0:
                    nvalue.append(float(header[par_begin + ilis]))
                    ilis += 1
                else:
                    assert ink[il] == 2
                    t = EMAN2_cppwrap.Transform()
                    tmp = []
                    for iii in range(par_begin + ilis, par_begin + ilis + 12):
                        tmp.append(float(header[iii]))
                    t.set_matrix(tmp)
                    ilis += 12
                    nvalue.append(t)
            ISID = list_params.count("ID")
            if ISID == 0:
                imm = im
            else:
                imm = nvalue[list_params.index("ID")]
                # read head, set params, and write it
            dummy = EMAN2_cppwrap.EMData()
            dummy.read_image(stack, imm, True)
            set_arb_params(dummy, nvalue, list_params)
            dummy.write_image(
                stack,
                dummy.get_attr_default("ID", im),
                EMAN2_cppwrap.EMUtil.ImageType.IMAGE_HDF,
                True,
            )


def send_attr_dict(main_node, data, list_params, image_start, image_end, comm=-1):

    #  This function is called from a node other than the main node

    if comm == -1:
        comm = mpi.MPI_COMM_WORLD
    TransType = type(EMAN2_cppwrap.Transform())
    mpi.mpi_send(
        [image_start, image_end],
        2,
        mpi.MPI_INT,
        main_node,
        sp_global_def.SPARX_MPI_TAG_UNIVERSAL,
        comm,
    )
    nvalue = []
    for im in range(image_start, image_end):
        value = get_arb_params(data[im - image_start], list_params)
        for il in range(len(value)):
            if type(value[il]) is int:
                nvalue.append(float(value[il]))
            elif type(value[il]) is float:
                nvalue.append(value[il])
            elif type(value[il]) is TransType:
                m = value[il].get_matrix()
                assert len(m) == 12
                for f in m:
                    nvalue.append(f)
    mpi.mpi_send(
        nvalue,
        len(nvalue),
        mpi.MPI_FLOAT,
        main_node,
        sp_global_def.SPARX_MPI_TAG_UNIVERSAL,
        comm,
    )


def recv_attr_dict_bdb(
    main_node, stack, data, list_params, image_start, image_end, number_of_proc, comm=-1
):

    #  bdb version!
    # This is done on the main node, so for images from the main node, simply write headers

    if comm == -1:
        comm = mpi.MPI_COMM_WORLD

    DB = EMAN2db.db_open_dict(stack)
    TransType = type(EMAN2_cppwrap.Transform())
    # prepare keys for float/int
    value = get_arb_params(data[0], list_params)
    ink = []
    len_list = 0
    ISID = -1
    for il in range(len(list_params)):
        if list_params[il] == "ID":
            ISID = il
        if type(value[il]) is int:
            ink.append(1)
            len_list += 1
        elif type(value[il]) is float:
            ink.append(0)
            len_list += 1
        elif type(value[il]) is TransType:
            ink.append(2)
            len_list += 12
    ldis = []
    headers = []
    for n in range(number_of_proc):
        if n != main_node:
            dis = mpi.mpi_recv(
                2, mpi.MPI_INT, n, sp_global_def.SPARX_MPI_TAG_UNIVERSAL, comm
            )
            img_begin = int(dis[0])
            img_end = int(dis[1])
            header = mpi.mpi_recv(
                len_list * (img_end - img_begin),
                mpi.MPI_FLOAT,
                n,
                sp_global_def.SPARX_MPI_TAG_UNIVERSAL,
                comm,
            )
            for im in range(img_begin, img_end):
                par_begin = (im - img_begin) * len_list
                nvalue = []
                ilis = 0
                for il in range(len(list_params)):
                    if ink[il] == 1:
                        nvalue.append(int(header[par_begin + ilis]))
                        ilis += 1
                    elif ink[il] == 0:
                        nvalue.append(float(header[par_begin + ilis]))
                        ilis += 1
                    else:
                        assert ink[il] == 2
                        t = EMAN2_cppwrap.Transform()
                        tmp = []
                        for iii in range(par_begin + ilis, par_begin + ilis + 12):
                            tmp.append(float(header[iii]))
                        t.set_matrix(tmp)
                        ilis += 12
                        nvalue.append(t)
                if ISID == -1:
                    imm = im
                else:
                    imm = nvalue[ISID]
                for i in range(len(list_params)):
                    if list_params[i] != "ID":
                        DB.set_attr(imm, list_params[i], nvalue[i])
        else:
            for n in range(image_start, image_end):
                ID = data[n - image_start].get_attr_default("ID", n)
                for param in list_params:
                    if param != "ID":
                        DB.set_attr(ID, param, data[n - image_start].get_attr(param))
    DB.close()


def print_begin_msg(program_name, onscreen=False):

    t = 100
    stars = "*" * t
    string = (
        "Beginning of the program "
        + program_name
        + ": "
        + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    )
    s = old_div((t - len(string)), 2)
    spacing = " " * s
    if onscreen:
        sp_global_def.sxprint(stars)
        sp_global_def.sxprint(spacing + string)
        sp_global_def.sxprint(stars)
    else:
        print_msg(stars + "\n")
        print_msg(spacing + string + "\n")
        print_msg(stars + "\n")


def print_end_msg(program_name, onscreen=False):

    t = 100
    stars = "*" * t
    string = (
        "End of the program "
        + program_name
        + ": "
        + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    )
    s = old_div((t - len(string)), 2)
    spacing = " " * s
    if onscreen:
        sp_global_def.sxprint(stars)
        sp_global_def.sxprint(spacing + string)
        sp_global_def.sxprint(stars)
    else:
        print_msg(stars + "\n")
        print_msg(spacing + string + "\n")
        print_msg(stars + "\n")


def print_msg(msg):

    if sp_global_def.IS_LOGFILE_OPEN == False:
        sp_global_def.LOGFILE_HANDLE = open(sp_global_def.LOGFILE, "w")
        sp_global_def.IS_LOGFILE_OPEN = True
    if sp_global_def.BATCH:
        sp_global_def.LOGFILE_HANDLE.write(msg)
    else:
        sys.stdout.write(msg)
        sp_global_def.LOGFILE_HANDLE.write(msg)
    sp_global_def.LOGFILE_HANDLE.flush()


def read_fsc(filename):

    f = open(filename, "r")
    fscc = None
    line = f.readline()
    while len(line) > 0:
        items = line.split()
        if fscc is None:
            fscc = [None] * len(items)
            for i in range(len(items)):
                fscc[i] = []

        for i in range(len(items)):
            fscc[i].append(float(items[i]))

        line = f.readline()

    return fscc


"""Multiline Comment15"""


def circumference(img, inner=-1, outer=-1):
    """
	Compute average within a shell between inner and outer radius
	Subtract this shell average from the image
	Multiply image by a spherical mask with inner radius
	"""
    nx = img.get_xsize()
    ny = img.get_ysize()
    nz = img.get_zsize()
    if inner == -1:
        inner = old_div(nx, 2) - 2
        if outer <= inner:
            outer = inner + 1
    else:
        if outer <= inner:
            outer = inner + 1
    inner_sphere = model_circle(inner, nx, ny, nz)

    [mean_a, sigma, imin, imax] = EMAN2_cppwrap.Util.infomask(
        img, model_circle(outer, nx, ny, nz) - inner_sphere, True
    )
    inner_rest = model_blank(nx, ny, nz, 1.0) - inner_sphere
    return EMAN2_cppwrap.Util.muln_img(inner_sphere, img - mean_a)
    # return Util.addn_img(inner_sphere, Util.mult_scalar(inner_rest, mean_a ) )


def write_headers(filename, data, lima):
    """
	  write headers from files in data into a disk file called filename.
	  The filename has to be either hdf or bdb.
	  lima - list with positions in the disk files into which headers will be written,
		i.e., header from data[k] will be written into file number lima[k]
	  WARNING: this function will open and close DB library!
	"""

    ftp = file_type(filename)
    if ftp == "bdb":
        #  For unknown reasons this does not work on Linux, but works on Mac ??? Really?
        DB = EMAN2db.db_open_dict(filename)
        for i in range(len(lima)):
            DB.set_header(lima[i], data[i])
        DB.close()
        # for i in range(len(lima)):
        # 	data[i].write_image(filename, lima[i])
    elif ftp == "hdf":
        for i in range(len(lima)):
            data[i].write_image(
                filename, lima[i], EMAN2_cppwrap.EMUtil.ImageType.IMAGE_HDF, True
            )
    else:
        sp_global_def.ERROR("Unacceptable file format", "write_headers", 1)


def write_header(filename, data, lima):
    """
	  write header from a single file data into a disk file called filename.
	  The filename has to be either hdf or bdb.
	  lima - position in the disk files into which header will be written,
		i.e., header from data will be written into file number lima
	  WARNING: this function assums DB library is opened and will NOT close it!
	"""

    ftp = file_type(filename)
    if ftp == "bdb":
        DB = EMAN2db.db_open_dict(filename)
        DB.set_header(lima, data)
    elif ftp == "hdf":
        data.write_image(filename, lima, EMAN2_cppwrap.EMUtil.ImageType.IMAGE_HDF, True)
    else:
        sp_global_def.ERROR("Unacceptable file format", "write_headers", 1)


def file_type(name):
    if len(name) > 4:
        if name[:4] == "bdb:":
            return "bdb"
        elif name[-4:-3] == ".":
            return name[-3:]
        elif name[-5:] == ".mrcs":
            return "mrcs"
    sp_global_def.ERROR("Unacceptable file format", "file_type", 1)


def get_params2D(ima, xform="xform.align2d"):
    """
	  retrieve 2D alignment parameters from the header
	  alpha, tx, ty, mirror, scale
	"""
    d = EMAN2_cppwrap.Util.get_transform_params(ima, xform, "2D")
    return d["alpha"], d["tx"], d["ty"], d["mirror"], d["scale"]


def set_params2D(ima, p, xform="xform.align2d"):
    """
	  set 2D alignment parameters in the header
	  p = [alpha, tx, ty, mirror, scale]
	"""
    t = EMAN2_cppwrap.Transform(
        {
            "type": "2D",
            "alpha": p[0],
            "tx": p[1],
            "ty": p[2],
            "mirror": p[3],
            "scale": p[4],
        }
    )
    ima.set_attr(xform, t)


def get_params3D(ima, xform="xform.align3d"):
    """
	  retrieve 3D alignment parameters from the header
	  phi,theta, psi, tx, ty, tz, mirror,scale
	"""
    d = EMAN2_cppwrap.Util.get_transform_params(ima, xform, "spider")
    return (
        d["phi"],
        d["theta"],
        d["psi"],
        d["tx"],
        d["ty"],
        d["tz"],
        d["mirror"],
        d["scale"],
    )


def set_params3D(ima, p, xform="xform.align3d"):
    """
	  set 3D alignment parameters in the header
	  p = [phi,theta, psi, tx, ty, tz, mirror,scale]
	"""
    t = EMAN2_cppwrap.Transform(
        {
            "type": "spider",
            "phi": p[0],
            "theta": p[1],
            "psi": p[2],
            "tx": p[3],
            "ty": p[4],
            "tz": p[5],
            "mirror": p[6],
            "scale": p[7],
        }
    )
    ima.set_attr(xform, t)


def get_params_proj(ima, xform="xform.projection"):
    """
	  retrieve projection alignment parameters from the header
	  phi, theta, psi, s2x, s2y
	"""
    d = EMAN2_cppwrap.Util.get_transform_params(ima, xform, "spider")
    return d["phi"], d["theta"], d["psi"], -d["tx"], -d["ty"]


def set_params_proj(ima, p, xform="xform.projection"):
    """
	  set projection alignment parameters in the header
	  p = [phi, theta, psi, s2x, s2y]
	"""

    t = EMAN2_cppwrap.Transform(
        {"type": "spider", "phi": p[0], "theta": p[1], "psi": p[2]}
    )
    t.set_trans(EMAN2_cppwrap.Vec2f(-p[3], -p[4]))
    ima.set_attr(xform, t)


def get_ctf(ima):
    """
	  recover numerical values of CTF parameters from EMAN2 CTF object stored in a header of the input image
	  order of returned parameters:
		defocus, cs, voltage, apix, bfactor, ampcont, astigmatism amplitude, astigmatism angle
	"""

    ctf_params = ima.get_attr("ctf")
    return (
        ctf_params.defocus,
        ctf_params.cs,
        ctf_params.voltage,
        ctf_params.apix,
        ctf_params.bfactor,
        ctf_params.ampcont,
        ctf_params.dfdiff,
        ctf_params.dfang,
    )


def same_ctf(c1, c2):
    """
	  Compare two CTF objects and return True if they are the same
	"""
    return c1.to_string() == c2.to_string()


def generate_ctf(p):
    """
	  generate EMAN2 CTF object using values of CTF parameters given in the list p
	  order of parameters:
		p = [defocus, cs, voltage, apix, bfactor, ampcont, astigmatism_amplitude, astigmatism_angle]
		[ microns, mm, kV, Angstroms, A^2, microns, microns, radians]
	"""

    defocus = p[0]
    cs = p[1]
    voltage = p[2]
    pixel_size = p[3]
    bfactor = p[4]
    amp_contrast = p[5]

    if (
        defocus > 100
    ):  # which means it is very likely in Angstrom, therefore we are using the old convention
        defocus *= 1.0e-4

    ctf = EMAN2_cppwrap.EMAN2Ctf()
    if len(p) == 6:
        ctf.from_dict(
            {
                "defocus": defocus,
                "cs": cs,
                "voltage": voltage,
                "apix": pixel_size,
                "bfactor": bfactor,
                "ampcont": amp_contrast,
            }
        )
    elif len(p) == 8:
        ctf.from_dict(
            {
                "defocus": defocus,
                "cs": cs,
                "voltage": voltage,
                "apix": pixel_size,
                "bfactor": bfactor,
                "ampcont": amp_contrast,
                "dfdiff": p[6],
                "dfang": p[7],
            }
        )
    else:
        sp_global_def.ERROR(
            "Incorrect number of entries on a list, cannot generate CTF",
            "generate_ctf",
            0,
        )
        return None
    return ctf


def delete_bdb(name):
    """
	  Delete bdb stack
	"""

    a = EMAN2db.db_open_dict(name)
    EMAN2db.db_remove_dict(name)


def disable_bdb_cache():

    EMAN2db.BDB_CACHE_DISABLE = True


"""Multiline Comment16"""


def getvec(phi, tht):

    if tht > 180.0:
        tht -= 180.0
        phi += 180.0
    if tht > 90.0:
        tht = 180.0 - tht
        phi += 180.0

    assert tht <= 90.0

    qt = numpy.radians(tht)
    qp = numpy.radians(phi)
    qs = numpy.sin(qt)

    x = qs * numpy.cos(qp)
    y = qs * numpy.sin(qp)
    z = numpy.cos(qt)

    return (x, y, z)


def getfvec(phi, tht):

    qt = numpy.radians(tht)
    qp = numpy.radians(phi)
    qs = numpy.sin(qt)
    x = qs * numpy.cos(qp)
    y = qs * numpy.sin(qp)
    z = numpy.cos(qt)

    return (x, y, z)


def nearest_fang(vecs, phi, tht):
    """
		vecs = [ [x0,y0,z0], [x1,y1,z1], ...]
	"""

    vec = getfvec(phi, tht)
    return EMAN2_cppwrap.Util.nearest_fang(vecs, vec[0], vec[1], vec[2])[0]


"""Multiline Comment17"""


# This is in python, it is very slow, we keep it just for comparison, use Util.assign_projangles instead


def nearest_many_full_k_projangles(
    reference_normals, angles, howmany=1, sym_class=None
):
    #

    # refnormal = normals[:]
    assignments = [-1] * len(angles)
    if sym_class.sym[:2] == "c1":
        for i, q in enumerate(angles):
            ref = getfvec(q[0], q[1])
            assignments[i] = EMAN2_cppwrap.Util.nearest_fang_select(
                reference_normals, ref[0], ref[1], ref[2], howmany
            )
    else:
        for i, q in enumerate(angles):
            ancordir = angles_to_normals(sym_class.symmetry_neighbors([q[:3]]))
            assignments[i] = EMAN2_cppwrap.Util.nearest_fang_sym(
                ancordir, reference_normals, len(ancordir), howmany
            )

    return assignments


"""Multiline Comment18"""


def assign_projdirs_f(projdirs, refdirs, neighbors):
    #  projdirs - data
    #  refdirs  - templates, each template has neighbors related copies
    #  output - list of lists, ofr each of refdirs/neighbors there is a list of projdirs indexes that are closest to it
    """
	qsti = [-1]*len(projdirs)
	for i,q in enumerate(projdirs):
		dn = -2.0
		for l in xrange(len(refdirs)/neighbors):
			sq = -2.0
			for k in xrange(l*neighbors,(l+1)*neighbors):
				sq = max(q[0]*refdirs[k][0] + q[1]*refdirs[k][1] + q[2]*refdirs[k][2], sq)
			if(sq > dn):
				dn = sq
				this = l
		qsti[i] = this
	"""
    #  Create a list that for each projdirs contains an index of the closest refdirs/neighbors
    qsti = EMAN2_cppwrap.Util.assign_projdirs_f(projdirs, refdirs, neighbors)
    assignments = [[] for i in range(old_div(len(refdirs), neighbors))]
    for i in range(len(projdirs)):
        assignments[qsti[i]].append(i)

    return assignments


"""Multiline Comment19"""

"""Multiline Comment20"""


"""Multiline Comment21"""


#  Wrappers for new angular functions


def angles_to_normals(angles):
    temp = EMAN2_cppwrap.Util.angles_to_normals(angles)
    return [[temp[l * 3 + i] for i in range(3)] for l in range(len(angles))]


"""


def symmetry_related(angles, symmetry):  # replace by s.symmetry_related
	if( (symmetry[0] == "c") or (symmetry[0] == "d") ):
		temp = Util.symmetry_related(angles, symmetry)
		nt = len(temp)/3
		return [[temp[l*3+i] for i in xrange(3)] for l in xrange(nt) ]
	else:
		from EMAN2 import Transform
		neighbors = []
		junk = Transform({"type":"spider","phi":angles[0],"theta":angles[1],"psi":angles[2]})
		junk = junk.get_sym_proj(symmetry)
		for p in junk:
			d = p.get_params("spider")
			neighbors.append([d["phi"],d["theta"],d["psi"]])
		return neighbors


def symmetry_related_normals(angles, symmetry):
	from EMAN2 import Transform
	neighbors = []
	junk = Transform({"type":"spider","phi":angles[0],"theta":angles[1],"psi":angles[2]})
	junk = junk.get_sym_proj(symmetry)
	for p in junk:
		neighbors.append(p.get_matrix()[8:11])
	return neighbors
"""


def angular_occupancy(params, angstep=15.0, sym="", method="S", inc_mirror=0):

    # smc = symclass(sym)
    # eah = smc.even_angles(angstep, inc_mirror=inc_mirror, method=method)

    # leah = len(eah)
    # u = []
    # for q in eah:
    # 	# sxprint("q",q)
    # 	m = smc.symmetry_related([(180.0 + q[0]) % 360.0, 180.0 - q[1], 0.0])
    # 	# sxprint("m",m)
    # 	itst = len(u)
    # 	for c in m:
    # 		# sxprint("c",c)
    # 		if smc.is_in_subunit(c[0], c[1], 1):
    # 			# sxprint(" is in 1")
    # 			if not smc.is_in_subunit(c[0], c[1], inc_mirror):
    # 				# sxprint("  outside")
    # 				u.append(c)
    # 				break
    # 	if len(u) != itst + 1:
    # 		u.append(q)  #  This is for exceptions that cannot be easily handled
    # 		"""
    # 		sxprint(q)
    # 		sxprint(m)
    # 		ERROR("balance angles","Fill up upper",1)
    # 		"""
    # seaf = []
    # for q in eah + u:
    # 	seaf += smc.symmetry_related(q)
    #
    # lseaf = len(seaf) / (2 * leah)
    # # sxprint(lseaf)
    # # for i,q in enumerate(seaf):  sxprint(" seaf  ",i,q)
    # # sxprint(seaf)
    # seaf = angles_to_normals(seaf)
    #
    # occupancy = [[] for i in range(leah)]
    #
    # for i, q in enumerate(params):
    # 	l = nearest_fang(seaf, q[0], q[1])
    # 	l = l / lseaf
    # 	if l >= leah:
    # 		l = l - leah
    # 	occupancy[l].append(i)
    # for i,q in enumerate(occupancy):
    # 	if q:
    # 		sxprint("  ",i,q,eah[i])
    # delta = angstep
    # data = params
    if isinstance(sym, sp_fundamentals.symclass):
        sym_class = sym
    else:
        sym_class = sp_fundamentals.symclass(sym)

    eah = sym_class.even_angles(delta=angstep, method=method, inc_mirror=inc_mirror)
    # symclass.even_angles(delta=delta, method=method, inc_mirror=inc_mirror)
    sym_class.build_kdtree()
    occupancy = sym_class.find_k_nearest_neighbors(params, k=1, tolistconv=False)
    # occupancy = indices
    return occupancy, eah


def angular_histogram(params, angstep=15.0, sym="c1", method="S", inc_mirror=0):

    occupancy, eah = angular_occupancy(
        params, angstep, sym, method, inc_mirror=inc_mirror
    )

    # return [len(q) for q in occupancy], eah

    # occupancy = occupancy.astype(int)   # Because of typecasting error

    radius_array = matplotlib.numpy.bincount(occupancy.flatten(),  minlength=len(eah))

    return radius_array, eah


def balance_angular_distribution(params, max_occupy=-1, angstep=15.0, sym="c1"):
    data = sp_fundamentals.symclass(sym).reduce_anglesets(
        matplotlib.numpy.array(params)[:, 0:3], inc_mirror=0, tolistconv=False
    )
    occupancy, eah = angular_occupancy(data, angstep, sym, method="S")

    if max_occupy > 0:
        outo = []
        for i in set(occupancy.flatten()):
            # for l, q in enumerate(occupancy):
            idxs = matplotlib.numpy.where(
                matplotlib.numpy.array(occupancy).flatten() == i
            )[0]
            matplotlib.numpy.random.shuffle(idxs)

            findalidxs = idxs[:max_occupy]
            outo.extend(findalidxs)
            # shuffle(q)
            # q = q[:max_occupy]
            # outo += q
            sp_global_def.sxprint(
                "  %10d   %10d        %6.1f   %6.1f"
                % (i, len(idxs), eah[i][0], eah[i][1])
            )
            # sxprint(l,len(q),q)
        outo.sort()

        # write_text_file(outo,"select.txt")
        return outo, matplotlib.numpy.array(params)[outo]
    else:
        # for l,q in enumerate(occupancy):
        # 	sxprint("  %10d   %10d        %6.1f   %6.1f"%(l,len(q),eah[l][0],eah[l][1]))
        return occupancy


def symmetry_neighbors(angles, symmetry):
    #  input is a list of lists  [[phi0,theta0,psi0],[phi1,theta1,psi1],...]
    #  output is [[phi0,theta0,psi0],[phi0,theta0,psi0]_SYM1,...,[phi1,theta1,psi1],[phi1,theta1,psi1]_SYM1,...]
    temp = EMAN2_cppwrap.Util.symmetry_neighbors(angles, symmetry)
    nt = old_div(len(temp), 3)
    return [[temp[l * 3 + i] for i in range(3)] for l in range(nt)]
    #  We could make it a list of lists
    # mt = len(angles)
    # nt = len(temp)/mt/3
    # return [[ [temp[m*3*nt+l*3+i] for i in xrange(3)] for l in xrange(nt)] for m in xrange(mt) ]


# def nearest_angular_direction(normals, vect, symmetry):


def rotation_between_anglesets(agls1, agls2):
    """
	  Find an overall 3D rotation (phi theta psi) between two sets of Eulerian angles. (psi irrelevant)
	  The two sets have to have the same number of elements and it is assumed that k'th element on the first
	  list corresponds to the k'th element on the second list.
	  Input: two lists [[phi1, theta1, psi1], [phi2, theta2, psi2], ...].  Second list is considered reference.
	  Output: overall rotation phi, theta, psi that has to be applied to the first list (agls1) so resulting
		angles will agree with the second list.
	  Note: all angles have to be in spider convention.
	  For details see: Appendix in Penczek, P., Marko, M., Buttle, K. and Frank, J.:  Double-tilt electron tomography.  Ultramicroscopy 60:393-410, 1995.
	"""

    def ori2xyz(ori):
        if type(ori) == list:
            phi, theta, psi = ori[:3]
        else:
            # it has to be Transformation object
            d = ori.get_params("spider")
            phi = d["phi"]
            theta = d["theta"]
            # psi   = d["psi"]

        phi = numpy.radians(phi)
        theta = numpy.radians(theta)
        sint = numpy.sin(theta)
        x = sint * numpy.sin(phi)
        y = sint * numpy.cos(phi)
        z = numpy.cos(theta)

        return [x, y, z]

    N = len(agls1)
    if N != len(agls2):
        sp_global_def.ERROR(
            "rotation_between_anglesets", "Both lists must have the same length", 1
        )
        return -1
    if N < 2:
        sp_global_def.ERROR(
            "rotation_between_anglesets",
            "At least two orientations are required in each list",
            1,
        )
        return -1

    U1 = [ori2xyz(q) for q in agls1]
    U2 = [ori2xyz(q) for q in agls2]

    # compute all Suv with uv = {xx, xy, xz, yx, ..., zz}
    Suv = [0] * 9

    nbori = len(U1)
    for i in range(3):
        for j in range(3):
            for s in range(nbori):
                Suv[j + 3 * i] += U2[s][i] * U1[s][j]

                # create matrix N
    N = numpy.array(
        [
            [
                Suv[0] + Suv[4] + Suv[8],
                Suv[5] - Suv[7],
                Suv[6] - Suv[2],
                Suv[1] - Suv[3],
            ],
            [
                Suv[5] - Suv[7],
                Suv[0] - Suv[4] - Suv[8],
                Suv[1] + Suv[3],
                Suv[6] + Suv[2],
            ],
            [
                Suv[6] - Suv[2],
                Suv[1] + Suv[3],
                -Suv[0] + Suv[4] - Suv[8],
                Suv[5] + Suv[7],
            ],
            [
                Suv[1] - Suv[3],
                Suv[6] + Suv[2],
                Suv[5] + Suv[7],
                -Suv[0] - Suv[4] + Suv[8],
            ],
        ]
    )

    # eigenvector corresponding to the most positive eigenvalue
    val, vec = numpy.linalg.eig(N)
    q0, qx, qy, qz = vec[:, val.argmax()]
    # create quaternion Rot matrix
    r = [
        [
            q0 * q0 - qx * qx + qy * qy - qz * qz,
            2 * (qy * qx + q0 * qz),
            2 * (qy * qz - q0 * qx),
        ],
        [
            2 * (qx * qy - q0 * qz),
            q0 * q0 + qx * qx - qy * qy - qz * qz,
            2 * (qx * qz + q0 * qy),
        ],
        [
            2 * (qz * qy + q0 * qx),
            2 * (qz * qx - q0 * qy),
            q0 * q0 - qx * qx - qy * qy + qz * qz,
        ],
    ]

    return sp_fundamentals.recmat(r)


"""Multiline Comment22"""


def get_pixel_size(img):
    """
	  Retrieve pixel size from the header.
	  We check attribute Pixel_size and also pixel size from ctf object, if exisits.
	  If the two are different or if the pixel size is not set, return -1.0 and print a warning.
	"""
    p1 = img.get_attr_default("apix_x", -1.0)
    cc = img.get_attr_default("ctf", None)
    if cc == None:
        p2 = -1.0
    else:
        p2 = round(cc.apix, 3)
    if p1 == -1.0 and p2 == -1.0:
        # ERROR("Pixel size not set", "get_pixel_size", 0)
        return -1.0
    elif p1 > -1.0 and p2 > -1.0:
        # if abs(p1-p2) >= 0.001:
        # 	ERROR("Conflict between pixel size in attribute and in ctf object", "get_pixel_size", 0)
        # pixel size is positive, so what follows omits -1 problem
        return max(p1, p2)
    else:
        return max(p1, p2)


def set_pixel_size(img, pixel_size):
    """
	  Set pixel size in the header.
	  Set attribute Pixel_size and also pixel size in ctf object, if exists.
	"""
    nz = img.get_zsize()
    img.set_attr("apix_x", round(pixel_size, 3))
    img.set_attr("apix_y", round(pixel_size, 3))
    img.set_attr("apix_z", round(pixel_size, 3))
    cc = img.get_attr_default("ctf", None)
    if cc:
        cc.apix = pixel_size
        img.set_attr("ctf", cc)


def lacos(x):
    """
		compute acos(x) in degrees after enforcing -1<=x<=1
	"""

    return numpy.degrees(math.acos(max(-1.0, min(1.0, x))))


def nearest_proj(proj_ang, img_per_grp=100, List=[]):
    def ang_diff(v1, v2):
        # The first return value is the angle between two vectors
        # The second return value is whether we need to mirror one of them (0 - no need, 1 - need)

        v = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]
        if v >= 0:
            return lacos(v), 0
        else:
            return lacos(-v), 1

    def get_ref_ang_list(delta, sym):
        ref_ang = even_angles(delta, symmetry=sym)
        ref_ang_list = [0.0] * (len(ref_ang) * 2)
        for i in range(len(ref_ang)):
            ref_ang_list[2 * i] = ref_ang[i][0]
            ref_ang_list[2 * i + 1] = ref_ang[i][1]
        return ref_ang_list, len(ref_ang)

    def binary_search(a, x):
        N = len(a)
        begin = 0
        end = N - 1
        while begin <= end:
            mid = old_div((begin + end), 2)
            if a[mid] == x:
                return mid
            if a[mid] < x:
                begin = mid + 1
            else:
                end = mid - 1
        return -1

    def binary_search_l(a, x):
        # This function returns an index i such that i is the smallest number
        # such that when t >= i, a[t] >= x
        N = len(a)
        t = binary_search(a, x)
        if t != -1:
            while t - 1 >= 0 and a[t - 1] == a[t]:
                t -= 1
            return t
        else:
            if x > a[N - 1]:
                return -1
            if x < a[0]:
                return 0
            begin = 0
            end = N - 2
            while end >= begin:
                mid = old_div((begin + end), 2)
                if x > a[mid] and x < a[mid + 1]:
                    break
                if x < a[mid]:
                    end = mid - 1
                else:
                    begin = mid + 1
            return mid + 1

    def binary_search_r(a, x):
        # This function returns an index i such that i is the largest number
        # such that when t <= i, a[t] <= x
        N = len(a)
        t = binary_search(a, x)
        if t != -1:
            while t + 1 <= N - 1 and a[t + 1] == a[t]:
                t += 1
            return t
        else:
            if x > a[N - 1]:
                return N - 1
            if x < a[0]:
                return -1
            begin = 0
            end = N - 2
            while end >= begin:
                mid = old_div((begin + end), 2)
                if x > a[mid] and x < a[mid + 1]:
                    break
                if x < a[mid]:
                    end = mid - 1
                else:
                    begin = mid + 1
            return mid

    N = len(proj_ang)
    if len(List) == 0:
        List = list(range(N))
    if N < img_per_grp:
        sp_global_def.sxprint(
            "Error: image per group larger than the number of particles!"
        )
        exit()
    phi_list = [[0.0, 0] for i in range(N)]
    theta_list = [[0.0, 0] for i in range(N)]
    vec = [None] * N
    for i in range(N):
        phi = proj_ang[i][0]
        theta = proj_ang[i][1]
        vec[i] = getfvec(phi, theta)
        if theta > 90.0:
            theta = 180.0 - theta
            phi += 180.0
        phi = phi % 360.0
        phi_list[i][0] = phi
        phi_list[i][1] = i
        theta_list[i][0] = theta
        theta_list[i][1] = i
    theta_list.sort()
    phi_list.sort()
    theta_list_l = [0.0] * N
    phi_list_l = [0.0] * N
    for i in range(N):
        theta_list_l[i] = theta_list[i][0]
        phi_list_l[i] = phi_list[i][0]

    g = [[360.0, 0, 0] for i in range(N)]
    proj_list = []
    mirror_list = []
    neighbor = [0] * img_per_grp
    # neighbor2 = [0]*img_per_grp
    dis = [0.0] * img_per_grp
    # dis2      = [0.0]*img_per_grp
    mirror = [0] * img_per_grp
    S = {0}
    T = {0}
    # tt1 = time()
    for i in range(len(List)):
        k = List[i]
        # print "\nCase #%3d: Testing projection %6d"%(i, k)
        # t1 = time()
        phi = proj_ang[k][0]
        theta = proj_ang[k][1]
        if theta > 90.0:
            theta = 180.0 - theta
            phi += 180.0
        phi = phi % 360.0
        delta = 0.01
        while True:
            min_theta = max(0.0, theta - delta)
            max_theta = min(90.0, theta + delta)
            if min_theta == 0.0:
                min_phi = 0.0
                max_phi = 360.0
            else:
                dphi = min(old_div(delta, (2 * min_theta)) * 180.0, 180.0)
                min_phi = phi - dphi
                max_phi = phi + dphi
                if min_phi < 0.0:
                    min_phi += 360.0
                if max_phi > 360.0:
                    max_phi -= 360.0
                if theta + delta > 90.0:
                    phi_mir = (phi + 180.0) % 360.0
                    min_phi_mir = phi_mir - dphi
                    max_phi_mir = phi_mir + dphi
                    if min_phi_mir < 0.0:
                        min_phi_mir += 360.0
                    if max_phi_mir > 360.0:
                        max_phi_mir -= 360.0

            phi_left_bound = binary_search_l(phi_list_l, min_phi)
            phi_right_bound = binary_search_r(phi_list_l, max_phi)
            theta_left_bound = binary_search_l(theta_list_l, min_theta)
            theta_right_bound = binary_search_r(theta_list_l, max_theta)
            if theta + delta > 90.0:
                phi_mir_left_bound = binary_search_l(phi_list_l, min_phi_mir)
                phi_mir_right_bound = binary_search_r(phi_list_l, max_phi_mir)
                # print delta
                # print min_phi, max_phi, min_theta, max_theta
                # print phi_left_bound, phi_right_bound, theta_left_bound, theta_right_bound
            if phi_left_bound < phi_right_bound:
                for j in range(phi_left_bound, phi_right_bound + 1):
                    S.add(phi_list[j][1])
            else:
                for j in range(phi_right_bound + 1):
                    S.add(phi_list[j][1])
                for j in range(phi_left_bound, N):
                    S.add(phi_list[j][1])
            if theta + delta > 90.0:
                if phi_mir_left_bound < phi_mir_right_bound:
                    for j in range(phi_mir_left_bound, phi_mir_right_bound + 1):
                        S.add(phi_list[j][1])
                else:
                    for j in range(phi_mir_right_bound + 1):
                        S.add(phi_list[j][1])
                    for j in range(phi_mir_left_bound, N):
                        S.add(phi_list[j][1])
            for j in range(theta_left_bound, theta_right_bound + 1):
                T.add(theta_list[j][1])
            v = list(T.intersection(S))
            S.clear()
            T.clear()
            if len(v) >= min(1.5 * img_per_grp, N):
                break
            delta *= 2
            del v

        for j in range(len(v)):
            d = ang_diff(vec[v[j]], vec[k])
            g[j][0] = d[0]
            if v[j] == k:
                g[j][
                    0
                ] = -1.0  # To ensure the image itself is always included in the group
            g[j][1] = d[1]
            g[j][2] = v[j]
        g[: len(v)] = sorted(g[: len(v)])
        for j in range(img_per_grp):
            neighbor[j] = g[j][2]
            dis[j] = g[j][0]
            mirror[j] = g[j][1] == 1
        proj_list.append(neighbor[:])
        mirror_list.append(mirror[:])
        # t2 = time()

        """Multiline Comment23"""
        # tt2 = time()
        # print tt2-tt1
    return proj_list, mirror_list


def findall(value, L, start=0):
    """
	 return a list of all indices of a value on the list L beginning from position start
	"""
    positions = []
    lL = len(L) - 1
    i = start - 1
    while i < lL:
        i += 1
        try:
            i = L.index(value, i)
            positions.append(i)
        except:
            break
    return positions


"""Multiline Comment24"""


"""Multiline Comment25"""


# ================ Iterator for list of images


class iterImagesList(object):
    images = []
    imagesIndexes = []
    position = -1

    def __init__(self, list_of_images, list_of_indexes=None):
        if list_of_indexes == None:
            self.images = list_of_images[:]
            self.imagesIndexes = list(range(len(self.images)))
        else:
            for i in list_of_indexes:
                self.images.append(list_of_images[i])
            self.imagesIndexes = list_of_indexes[:]

    def iterNo(self):
        return self.position

    def imageIndex(self):
        return self.imagesIndexes[self.position]

    def image(self):
        return self.images[self.position]

    def goToNext(self):
        if len(self.imagesIndexes) <= self.position:
            return False
        self.position += 1
        return self.position < len(self.imagesIndexes)

    def goToPrev(self):
        if 0 > self.position:
            return False
        self.position -= 1
        return self.position >= 0


# ================ Iterator for stack of images


class iterImagesStack(object):
    stackName = ""
    currentImage = None
    imagesIndexes = []
    position = -1

    def __init__(self, stack_name, list_of_indexes=None):
        if list_of_indexes == None:
            self.imagesIndexes = list(
                range(EMAN2_cppwrap.EMUtil.get_image_count(stack_name))
            )
        else:
            self.imagesIndexes = list_of_indexes[:]
        self.stackName = stack_name

    def iterNo(self):
        return self.position

    def imageIndex(self):
        return self.imagesIndexes[self.position]

    def image(self):
        if self.currentImage == None:
            self.currentImage = EMAN2_cppwrap.EMData()
            self.currentImage.read_image(
                self.stackName, self.imagesIndexes[self.position]
            )
        return self.currentImage

    def goToNext(self):
        self.currentImage = None
        if len(self.imagesIndexes) <= self.position:
            return False
        self.position += 1
        return self.position < len(self.imagesIndexes)

    def goToPrev(self):
        self.currentImage = None
        if 0 > self.position:
            return False
        self.position -= 1
        return self.position >= 0


def pack_message(data):
    """Convert data for transmission efficiently"""
    import sys
    if isinstance(data, str):
        if len(data) > 256:
            return b"C" + zlib.compress(data.encode('utf8'), 1)
        else:
            return b"S" + data.encode('utf8')
    else:
        d2x = pickle.dumps(data, 2)
        if len(d2x) > 256:
            return b"Z" + zlib.compress(d2x, 1)
        else:
            return b"O" + d2x



def unpack_message(msg):
    """Unpack a data payload prepared by pack_message"""

    if msg[0:1] == b"C":
        return zlib.decompress(msg[1:]).decode('latin1')
    elif msg[0:1] == b"S":
        return (msg[1:].decode('latin1'))
    elif msg[0:1] == b"Z":
        return pickle.loads(zlib.decompress((msg[1:])))
    elif msg[0:1] == b"O":
        return pickle.loads((msg[1:]))

    else:
        sp_global_def.ERROR(
            "ERROR: Invalid MPI message. Please contact developers. (%s)"
            % str(msg[:20])  )
        raise Exception("unpack_message")


def update_tag(
    communicator, target_rank
):  # TODO - it doesn't work when communicators are destroyed and recreated
    return 123456
    global statistics_send_recv
    if communicator not in statistics_send_recv:

        statistics_send_recv[communicator] = [0] * mpi.mpi_comm_size(communicator)
    statistics_send_recv[communicator][target_rank] += 1
    return statistics_send_recv[communicator][target_rank]


# ===================================== WRAPPER FOR MPI


def wrap_mpi_send(data, destination, communicator=None):

    if communicator == None:
        communicator = mpi.MPI_COMM_WORLD

    msg = pack_message(data)
    tag = update_tag(communicator, destination)
    # from mpi import mpi_comm_rank
    # print communicator, mpi_comm_rank(communicator), "send to", destination, tag
    mpi.mpi_send(
        msg, len(msg), mpi.MPI_CHAR, destination, tag, communicator
    )  # int MPI_Send( void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm )


def wrap_mpi_recv(source, communicator=None):

    if communicator == None:
        communicator = mpi.MPI_COMM_WORLD

    tag = update_tag(communicator, source)
    # from mpi import mpi_comm_rank
    # print communicator, mpi_comm_rank(communicator), "recv from", source, tag
    mpi.mpi_probe(source, tag, communicator)
    n = mpi.mpi_get_count(mpi.MPI_CHAR)
    msg = mpi.mpi_recv(n, mpi.MPI_CHAR, source, tag, communicator)
    return unpack_message(msg)



'''    print(msg[0:1], msg[1:2], msg[2:3], msg[3:-1])
    print(msg[0:1].decode('latin1') , msg[1:2].decode('latin1')  , msg[2:3].decode('latin1') , msg[3:-1].decode('latin1') )
    print([msg[i:i+1] for i in range(len(msg))])
    print([msg[i:i+1].decode('latin1') for i in range(len(msg))])'''

def wrap_mpi_bcast(data, root, communicator=None):

    if communicator == None:
        communicator = mpi.MPI_COMM_WORLD

    rank = mpi.mpi_comm_rank(communicator)

    if rank == root:
        msg = pack_message(data)
        n = struct.pack("I", len(msg))

    else:
        msg = None
        n = None

    sizeofdata = mpi.mpi_bcast(
        n, 4, mpi.MPI_CHAR, root, communicator
    )  # int MPI_Bcast ( void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm ).

    if mpi.mpi_comm_rank(communicator) == root:
        sizeofdata = n
    else:
        sizeofdata = b''.join([entry if entry != b'' else b'\x00' for entry in sizeofdata])

    n = struct.unpack("I", sizeofdata)[0]
    msgtobcast = mpi.mpi_bcast(
        msg, n, mpi.MPI_CHAR, root, communicator
    )  # int MPI_Bcast ( void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm )

    if mpi.mpi_comm_rank(communicator) == root:
        msgtobcast = msg
    else:
        msgtobcast = b''.join([entry if entry != b'' else b'\x00' for entry in msgtobcast])

    return unpack_message(msgtobcast)


# data must be a python list (numpy array also should be implemented)


def wrap_mpi_gatherv(data, root, communicator=None):

    if communicator == None:
        communicator = mpi.MPI_COMM_WORLD

    rank = mpi.mpi_comm_rank(communicator)
    procs = mpi.mpi_comm_size(communicator)

    out_array = None
    if rank == root:
        if type(data) is list:
            out_array = []
            for p in range(procs):
                if p == rank:
                    out_array.extend(data)
                else:
                    recv_data = wrap_mpi_recv(p, communicator)
                    out_array.extend(recv_data)
        else:
            raise Exception("wrap_mpi_gatherv: type of data not supported")
    else:
        wrap_mpi_send(data, root, communicator)

    return out_array


# def rearrange_ranks_of_processors_to_fit_round_robin_assignment():


def get_colors_and_subsets(
    main_node,
    mpi_comm,
    my_rank,
    shared_comm,
    sh_my_rank,
    masters_from_groups_vs_everything_else_comm,
):
    """Multiline Comment26"""

    sh_my_ranks = [my_rank]
    sh_my_ranks = wrap_mpi_gatherv(sh_my_ranks, main_node, shared_comm)

    group_infos = None
    if sh_my_rank == main_node:
        # executed only by masters from groups
        group_infos = [my_rank, sh_my_ranks]
        group_infos = wrap_mpi_gatherv(
            group_infos, main_node, masters_from_groups_vs_everything_else_comm
        )

    mpi.mpi_barrier(mpi_comm)

    group_infos = wrap_mpi_bcast(group_infos, main_node, mpi_comm)
    number_of_groups = old_div(len(group_infos), 2)

    for i in range(number_of_groups):
        if my_rank in group_infos[2 * i + 1]:
            color = i  # group_infos[2*i]
            break

    number_of_processes_in_each_group = []
    for i in range(number_of_groups):
        number_of_processes_in_each_group.append(len(group_infos[2 * i + 1]))

    balanced_processor_load_on_nodes = len(set(number_of_processes_in_each_group)) == 1

    return color, number_of_groups, balanced_processor_load_on_nodes


def wrap_mpi_split(comm, no_of_groups):
    """

	Takes the processes of a communicator (comm) and splits them in groups (no_of_groups).
	Each subgroup of processes has ids generated from 0 to number of processes per group - 1.
	Consecutive global process ids have consecutive subgroup process ids.

	"""

    nproc = mpi.mpi_comm_size(comm)
    myid = mpi.mpi_comm_rank(comm)

    no_of_proc_per_group = old_div(nproc, no_of_groups)
    color = old_div(myid, no_of_proc_per_group)
    key = myid % no_of_proc_per_group

    return mpi.mpi_comm_split(comm, color, key)


def get_dist(c1, c2):

    d = numpy.sqrt((c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2)
    return d


def eliminate_moons(my_volume, moon_elimination_params):
    """
	moon_elimination_params[0] - mass in KDa
	moon_elimination_params[1] - pixel size in A
	"""

    histogram_threshold = (
        my_volume.find_3d_threshold(
            moon_elimination_params[0], moon_elimination_params[1]
        )
        * 1.1
    )
    # clean11 88x88,  4.84 px/A 750 kDa

    my_volume_binarized = sp_morphology.binarize(my_volume, histogram_threshold)
    # my_volume_binarized.write_image ("my_volume_binarized.hdf")
    my_volume_binarized_with_no_moons = EMAN2_cppwrap.Util.get_biggest_cluster(
        my_volume_binarized
    )
    # my_volume_binarized_with_no_moons.write_image("my_volume_binarized_with_no_moons.hdf")
    volume_difference = my_volume_binarized - my_volume_binarized_with_no_moons
    # volume_difference.write_image("volume_difference.hdf")

    if (
        volume_difference.get_value_at(volume_difference.calc_max_index()) == 0
        and volume_difference.get_value_at(volume_difference.calc_min_index()) == 0
    ):
        return my_volume
    else:

        return gauss_edge(my_volume_binarized_with_no_moons) * my_volume

        # from utilities   import model_blank
        # # mask = model_blank(my_volume_binarized_with_no_moons.get_xsize(), my_volume_binarized_with_no_moons.get_ysize(), my_volume_binarized_with_no_moons.get_zsize())
        # # mask.to_one()
        # this is only in master


def combinations_of_n_taken_by_k(n, k):
    import functools
    return int(
        functools.reduce(
            lambda x, y: x * y, (fractions.Fraction(n - i, i + 1) for i in range(k)), 1
        )
    )


def cmdexecute(cmd, printing_on_success=True):

    # import subprocess  I do not know why this is not used. PAP

    # Remove the variable to avoid problems with subprograms not running with MPI in combination with print_timestamp in global_def
    try:
        mpi_comm_world_rank_value = os.environ["OMPI_COMM_WORLD_RANK"]
        del os.environ["OMPI_COMM_WORLD_RANK"]
    except KeyError:
        mpi_comm_world_rank_value = None

    outcome = os.system(cmd)
    if mpi_comm_world_rank_value is not None:
        os.environ["OMPI_COMM_WORLD_RANK"] = mpi_comm_world_rank_value

    line = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime()) + " =>"
    if outcome != 0:
        sp_global_def.sxprint(
            line,
            "ERROR!!   Command failed:  ",
            cmd,
            " return code of failed command: ",
            outcome,
        )
        return 0
    elif printing_on_success:
        sp_global_def.sxprint(line, "Executed successfully: ", cmd)
        return 1


def string_found_in_file(myregex, filename):

    pattern = re.compile(myregex)
    for line in open(filename):
        if re.findall(pattern, line) != []:
            return True
    return False


def get_latest_directory_increment_value(
    directory_location, directory_name, start_value=1, myformat="%03d"
):

    dir_count = start_value
    while os.path.isdir(directory_location + directory_name + myformat % (dir_count)):
        dir_count += 1
    if dir_count == start_value:
        return start_value
    return dir_count - 1


def if_error_then_all_processes_exit_program(error_status):

    # if "OMPI_COMM_WORLD_SIZE" not in os.environ:
    if (
        len(
            {
                "OMPI_COMM_WORLD_SIZE",
                "PMI_RANK",
                "PMI_ID",
                "LAMRANK",
                "MPI_RANKID",
                "MP_CHILD",
                "MP_CHILD",
                "MP_CHILD",
            }.intersection(set(os.environ))
        )
        == 0
    ):

        def my_mpi_comm_rank(n):
            return 0

        def my_mpi_bcast(*largs):
            return [largs[0]]

        def my_mpi_finalize():
            return None

        MY_MPI_INT, MY_MPI_COMM_WORLD = 0, 0
    else:

        my_mpi_comm_rank = mpi.mpi_comm_rank
        my_mpi_bcast = mpi.mpi_bcast
        my_mpi_finalize = mpi.mpi_finalize
        MY_MPI_INT = mpi.MPI_INT
        MY_MPI_COMM_WORLD = mpi.MPI_COMM_WORLD

    myid = my_mpi_comm_rank(MY_MPI_COMM_WORLD)
    if error_status != None and error_status != 0:
        error_status_info = error_status
        error_status = 1
    else:
        error_status = 0

    error_status = my_mpi_bcast(error_status, 1, MY_MPI_INT, 0, MY_MPI_COMM_WORLD)
    error_status = int(error_status[0])

    if error_status > 0:
        if myid == 0:
            if type(error_status_info) == type((1, 1)):
                if len(error_status_info) == 2:
                    frameinfo = error_status_info[1]
                    print_msg("***********************************\n")
                    print_msg("** Error: %s\n" % error_status_info[0])
                    print_msg("***********************************\n")
                    print_msg(
                        "** Location: %s\n"
                        % (frameinfo.filename + ":" + str(frameinfo.lineno))
                    )
                    print_msg("***********************************\n")
        sys.stdout.flush()
        my_mpi_finalize()
        exit()  # sys.exit(1)


"""Multiline Comment28"""


def getindexdata(stack, partids, partstack, myid, nproc):
    # The function will read from stack a subset of images specified in partids
    #   and assign to them parameters from partstack
    # So, the lengths of partids and partstack are the same.
    #  The read data is properly distributed among MPI threads.

    lpartids = read_text_file(partids)
    ndata = len(lpartids)
    partstack = read_text_row(partstack)

    if ndata < nproc:
        if myid < ndata:
            image_start = myid
            image_end = myid + 1
        else:
            image_start = 0
            image_end = 1
    else:
        image_start, image_end = sp_applications.MPI_start_end(ndata, nproc, myid)
    lpartids = lpartids[image_start:image_end]
    partstack = partstack[image_start:image_end]
    data = EMAN2_cppwrap.EMData.read_images(stack, lpartids)

    for i in range(len(partstack)):
        set_params_proj(data[i], partstack[i])
    return data


def store_value_of_simple_vars_in_json_file(
    filename,
    local_vars,
    exclude_list_of_vars=[],
    write_or_append="w",
    vars_that_will_show_only_size=[],
):
    import collections

    allowed_types = [type(None), bool, int, int, float, complex, str, bytes]

    local_vars_keys = list(local_vars.keys())

    my_vars = dict()
    for key in set(local_vars_keys) - set(exclude_list_of_vars):
        if type(local_vars[key]) in allowed_types:
            my_vars[key] = local_vars[key]
        elif type(local_vars[key]) in [list, tuple, type(set())]:
            if len({type(i) for i in local_vars[key]} - set(allowed_types)) == 0:
                if key in vars_that_will_show_only_size:
                    my_vars[key] = "%s with length: %d" % (
                        str(type(local_vars[key])),
                        len(local_vars[key]),
                    )
                else:
                    if type(local_vars[key]) == type(set()):
                        my_vars[key] = list(local_vars[key])
                    else:
                        my_vars[key] = local_vars[key]
        elif type(local_vars[key]) == dict:
            if (
                len(
                    {type(local_vars[key][i]) for i in local_vars[key]}
                    - set(allowed_types)
                )
                == 0
            ):
                my_vars[key] = local_vars[key]

    ordered_my_vars = collections.OrderedDict(sorted(my_vars.items()))

    with open(filename, write_or_append) as fp:
        json.dump(ordered_my_vars, fp, indent=2)
    fp.close()


def convert_json_fromunicode(data):
    import collections

    if isinstance(data, six.string_types):
        return str(data)
    elif isinstance(data, collections.Mapping):
        return dict(list(map(convert_json_fromunicode, iter(list(data.items())))))
    elif isinstance(data, collections.Iterable):
        return type(data)(list(map(convert_json_fromunicode, data)))
    else:
        return data


"""Multiline Comment30"""


"""Multiline Comment31"""

#### Used in the main programm


"""Multiline Comment32"""


"""Multiline Comment33"""


"""Multiline Comment35"""


def angle_between_projections_directions(proj1, proj2):
    """
	  It returns angle between two projections directions.
	  INPUT: two lists: [phi1, theta1] , [phi2, theta2]
	  OUTPUT: angle (in degrees)
	"""
    pass  # IMPORTIMPORTIMPORT from math import sin, cos, acos, radians, degrees
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_utilities import lacos

    theta1 = numpy.radians(proj1[1])
    theta2 = numpy.radians(proj2[1])
    cp1cp2_sp1sp2 = numpy.cos(numpy.radians(proj1[0]) - numpy.radians(proj2[0]))
    temp = numpy.sin(theta1) * numpy.sin(theta2) * cp1cp2_sp1sp2 + numpy.cos(
        theta1
    ) * numpy.cos(theta2)
    return lacos(temp)


getang3 = angle_between_projections_directions


"""Multiline Comment36"""


def angular_distribution_sabrina(params_file, output_folder, prefix, method, pixel_size, delta, symmetry, box_size, particle_radius, dpi, nth_percentile, do_print=True, exclude=None):
    matplotlib.pyplot.matplotlib.use("Agg")

    matplotlib.pyplot.rc('axes', axisbelow=True)

    def get_color(sorted_array):
        """
        Get the color for the 2D visual representation.

        Arguments:
        sorted_array - Array sorted by size.

        Returns:
        List of colors for each entry in the sorted_array.
        """
        sorted_normalized = numpy.true_divide(sorted_array, numpy.max(sorted_array))
        color_list = []
        for item in sorted_normalized:
            color_list.append((0, item, 1-item))
        return color_list

    COLUMN_Z = 2

    error_template = 'ERROR: {0}'
    error = False
    error_list = []

    if pixel_size <= 0:
        error_list.append('Pixel size cannot be smaller equals 0')
        error = True

    if box_size <= 0:
        error_list.append('Box size cannot be smaller equals 0')
        error = True

    if particle_radius <= 0:
        error_list.append('Particle radius cannot be smaller equals 0')
        error = True

    if delta <= 0:
        error_list.append('delta cannot be smaller equals 0')
        error = True

    if dpi <= 0:
        error_list.append('Dpi cannot be smaller equals 0')
        error = True

    if error:
        for entry in error_list:
            sp_global_def.sxprint(error_template.format(entry))
        return None
    else:
        if do_print:
           sp_global_def.sxprint('Values are valid')

    try:
        os.makedirs(output_folder)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.lexists(output_folder):
            if do_print:
               sp_global_def.sxprint('Output directory already exists: {0}'.format(output_folder))
        else:
            raise
    else:
        if do_print:
           sp_global_def.sxprint('Created output directory: {0}'.format(output_folder))
    sp_global_def.write_command(output_folder)

    COLUMN_X = 0
    COLUMN_Y = 1
    if prefix is None:
        prefix = os.path.basename(os.path.splitext(params_file)[0])

    # Import the parameters, assume the columns 0 (Phi) 1 (Theta) 2 (Psi) are present
    if do_print:
       sp_global_def.sxprint('Import projection parameter')
    data_params = numpy.atleast_2d(numpy.genfromtxt(params_file, usecols=(0, 1, 2)))
    if exclude is not None:
        mask = numpy.array(exclude)
        mask[mask < 0] = 0
        mask = mask.astype(numpy.bool)
        data_params = data_params[~mask]

    # If the symmetry is c0, do not remove mirror projections.
    if do_print:
       sp_global_def.sxprint('Reduce anglesets')
    if symmetry.endswith('_full'):
        symmetry = symmetry.rstrip('_full')
        inc_mirror = 1
    else:
        symmetry = symmetry
        inc_mirror = 0

    # Create 2 symclass objects.
    # One C1 object for the inital reference angles.
    # One related to the actual symmetry, to deal with mirror projections.
    sym_class = sp_fundamentals.symclass(symmetry)


    # data_params = [[0.0, 0.0, 0.0], [22,91,45]]
    # data_params = numpy.atleast_2d(data_params)
    if do_print:
       sp_global_def.sxprint('Reduce data to symmetry - This might take some time for high symmetries')
    # Reduce the parameter data by moving mirror projections into the non-mirror region of the sphere.
    data = sym_class.reduce_anglesets(data_params, inc_mirror=inc_mirror, tolistconv=False)
    # Create cartesian coordinates

    # sym_class.even_angles(delta= delta , method=method, inc_mirror=inc_mirror)
    # sym_class.build_kdtree()
    # indices = sym_class.find_k_nearest_neighbors(data, k=1, tolistconv=False)

    # radius_array = numpy.bincount(indices.flatten(), minlength=sym_class.angles.shape[0])

    radius_array, indices = angular_histogram(data,delta,sym =sym_class,method=method,inc_mirror=inc_mirror)

    angles_no_mirror_cart = sym_class.to_cartesian(sym_class.angles, tolistconv=False)

    nonzero_mask = numpy.nonzero(radius_array)
    radius_array = radius_array[nonzero_mask]


    # Calculate best width and length for the bins in 3D
    width = (pixel_size * particle_radius * numpy.radians(delta) * 2) / float(2 * 3)
    length = particle_radius * 0.2

# Create arrays for 2D plot + normalization
    nonzero_mask = list(nonzero_mask[0])
    sorted_radius = radius_array
    sorted_radius_plot = numpy.sort(radius_array)[::-1]
    array_x = numpy.arange(sorted_radius.shape[0])


    # Normalize the radii so that the 1% value - ordered values - is 1 (robust against outliers)
    quartile_idx = int((len(sorted_radius_plot)-1) * (100-nth_percentile) / 100.0)
    quartile_value = sorted_radius_plot[quartile_idx]
    percentile_idxes = numpy.where(radius_array == quartile_value)[0]
    radius_normalized = numpy.true_divide(radius_array, quartile_value)

    #Normalize for color - max value = yellow
    radius_normalized_color = numpy.true_divide(radius_array, numpy.max(radius_array))

    # Vector pointing to the center of the Box in chimera
    vector_center = 0.5 * numpy.array([
        box_size,
        box_size,
        box_size
        ])

    # Inner and outer vector. The bin starts at inner vector and stops at outer vector
    inner_vector = vector_center + angles_no_mirror_cart[nonzero_mask] * particle_radius

#Length here!!
    outer_vector = vector_center + angles_no_mirror_cart[nonzero_mask] * (particle_radius + radius_normalized.reshape(radius_normalized.shape[0], 1) * length )
    outer_vector_blue = vector_center + angles_no_mirror_cart[nonzero_mask] * (particle_radius + radius_normalized.reshape(radius_normalized.shape[0], 1) * length * 1.01 )

    # Adjust pixel size
    numpy.multiply(inner_vector, pixel_size, out=inner_vector)
    numpy.multiply(outer_vector, pixel_size, out=outer_vector)
    numpy.multiply(outer_vector_blue, pixel_size, out=outer_vector_blue)

    #get cmap for coloring of bins
    a=matplotlib.pyplot.get_cmap('viridis')
    colors = a(radius_normalized_color)

    # Create output bild file
    if do_print:
       sp_global_def.sxprint('Create bild file')
    output_bild_file = os.path.join(output_folder, '{0}.bild'.format(prefix))
    lightblue = [173/255.0, 216/255.0, 230/255.0, 1]
    with open(output_bild_file, 'w') as write:
        for idx, (inner, outer, outer_blue, color) in enumerate(zip(inner_vector, outer_vector, outer_vector_blue, colors)):
            if idx in percentile_idxes:
                write.write('.color {0} {1} {2}\n'.format(*lightblue))
                write.write(
                    '.cylinder {0:.3f} {1:.3f} {2:.3f} {3:.3f} {4:.3f} {5:.3f} {6:.3f}\n'.format(
                        outer[COLUMN_X],
                        outer[COLUMN_Y],
                        outer[COLUMN_Z],
                        outer_blue[COLUMN_X],
                        outer_blue[COLUMN_Y],
                        outer_blue[COLUMN_Z],
                        width
                        )
                    )
            write.write('.color {0} {1} {2}\n'.format(*color))
            write.write(
                '.cylinder {0:.3f} {1:.3f} {2:.3f} {3:.3f} {4:.3f} {5:.3f} {6:.3f}\n'.format(
                    inner[COLUMN_X],
                    inner[COLUMN_Y],
                    inner[COLUMN_Z],
                    outer[COLUMN_X],
                    outer[COLUMN_Y],
                    outer[COLUMN_Z],
                    width
                    )
                )

    """
    lina = numpy.argsort(radius_array)
    sorted_radius = radius_array[lina[::-1]]
    array_x = numpy.arange(sorted_radius.shape[0])
    angles_no_mirror = angles_no_mirror[lina[::-1]]
    nonzero_mask = list(nonzero_mask[0][lina[::-1]])

    """

    #"""
    

    # 2D distribution plot
    if do_print:
       sp_global_def.sxprint('Create 2D legend plot')
    output_bild_legend_png = os.path.join(output_folder, '{0}.png'.format(prefix))
    matplotlib.pyplot.axhline(quartile_value, 0, 1, color=lightblue, label='{0:.2f}th percentile'.format(nth_percentile), zorder=1)
    matplotlib.pyplot.bar(array_x, height=sorted_radius_plot, width=1, color=colors[numpy.argsort(radius_array)[::-1]], zorder=2)
    ax = matplotlib.pyplot.axes()
    ax.yaxis.grid(linestyle='--', linewidth=1, zorder=0)
    matplotlib.pyplot.xlabel('Bin / a.u.')
    matplotlib.pyplot.ylabel('Nr. of Particles')
    matplotlib.pyplot.legend()
    matplotlib.pyplot.savefig(output_bild_legend_png, dpi=dpi)
    matplotlib.pyplot.clf()
    """
    sxprint(array_x)
    sxprint(sorted_radius)
    sxprint(len(angles_no_mirror))
    sxprint(angles_no_mirror)
    """
    # 2D distribution txt file
    if do_print:
       sp_global_def. sxprint('Create 2D legend text file')

    output_bild_legend_txt = os.path.join(output_folder, '{0}.txt'.format(prefix))
    with open(output_bild_legend_txt, 'w') as write:
        for i in range(len(sorted_radius)):
            #	for value_x, value_y in zip(array_x, sorted_radius):
            value_x = '{0:6d}'.format(array_x[i])
            value_y = '{0:6d}'.format(sorted_radius[i])
            phi     = '{0:10f}'.format(sym_class.angles[nonzero_mask[i]][0])
            theta   = '{0:10f}'.format(sym_class.angles[nonzero_mask[i]][1])
            write.write('{0}\n'.format('\t'.join([value_x, value_y, phi, theta])))




def angular_distribution(
    params_file,
    output_folder,
    prefix,
    method,
    pixel_size,
    delta,
    symmetry,
    box_size,
    particle_radius,
    dpi,
    do_print=True,
    exclude=None,
):
    matplotlib.pyplot.matplotlib.use("Agg")

    def get_color(sorted_array):
        """
		Get the color for the 2D visual representation.

		Arguments:
		sorted_array - Array sorted by size.

		Returns:
		List of colors for each entry in the sorted_array.
		"""
        sorted_normalized = matplotlib.numpy.true_divide(
            sorted_array, matplotlib.numpy.max(sorted_array)
        )
        color_list = []
        for item in sorted_normalized:
            color_list.append((0, item, 1 - item))
        return color_list

    COLUMN_Z = 2

    # Use name of the params file as prefix if prefix is None
    # Sanity checks
    # sxprint('Check if values are valid')
    error_template = "ERROR: {0}"
    error = False
    error_list = []

    if pixel_size <= 0:
        error_list.append("Pixel size cannot be smaller equals 0")
        error = True

    if box_size <= 0:
        error_list.append("Box size cannot be smaller equals 0")
        error = True

    if particle_radius <= 0:
        error_list.append("Particle radius cannot be smaller equals 0")
        error = True

    if delta <= 0:
        error_list.append("delta cannot be smaller equals 0")
        error = True

    if dpi <= 0:
        error_list.append("Dpi cannot be smaller equals 0")
        error = True

    if error:
        for entry in error_list:
            sp_global_def.sxprint(error_template.format(entry))
        return None
    else:
        if do_print:
            sp_global_def.sxprint("Values are valid")

    try:
        os.makedirs(output_folder)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.lexists(output_folder):
            if do_print:
                sp_global_def.sxprint(
                    "Output directory already exists: {0}".format(output_folder)
                )
        else:
            raise
    else:
        if do_print:
            sp_global_def.sxprint("Created output directory: {0}".format(output_folder))
    sp_global_def.write_command(output_folder)
    COLUMN_X = 0
    COLUMN_Y = 1
    if prefix is None:
        prefix = os.path.basename(os.path.splitext(params_file)[0])

    # Import the parameters, assume the columns 0 (Phi) 1 (Theta) 2 (Psi) are present
    if do_print:
        sp_global_def.sxprint("Import projection parameter")
    data_params = matplotlib.numpy.atleast_2d(
        matplotlib.numpy.genfromtxt(params_file, usecols=(0, 1, 2))
    )
    if exclude is not None:
        mask = matplotlib.numpy.array(exclude)
        mask[mask < 0] = 0
        mask = mask.astype(matplotlib.numpy.bool)
        data_params = data_params[~mask]

    # If the symmetry is c0, do not remove mirror projections.
    if do_print:
        sp_global_def.sxprint("Reduce anglesets")
    if symmetry.endswith("_full"):
        symmetry = symmetry.rstrip("_full")
        inc_mirror = 1
    else:
        symmetry = symmetry
        inc_mirror = 0

    # Create 2 symclass objects.
    # One C1 object for the inital reference angles.
    # One related to the actual symmetry, to deal with mirror projections.
    sym_class = sp_fundamentals.symclass(symmetry)

    # data_params = [[0.0, 0.0, 0.0], [22,91,45]]
    # data_params = numpy.atleast_2d(data_params)
    if do_print:
        sp_global_def.sxprint(
            "Reduce data to symmetry - This might take some time for high symmetries"
        )
    # Reduce the parameter data by moving mirror projections into the non-mirror region of the sphere.
    data = sym_class.reduce_anglesets(
        data_params, inc_mirror=inc_mirror, tolistconv=False
    )
    # Create cartesian coordinates

    # sym_class.even_angles(delta= delta , method=method, inc_mirror=inc_mirror)
    # sym_class.build_kdtree()
    # indices = sym_class.find_k_nearest_neighbors(data, k=1, tolistconv=False)

    # radius_array = numpy.bincount(indices.flatten(), minlength=sym_class.angles.shape[0])

    radius_array, indices = angular_histogram(
        data, delta, sym=sym_class, method=method, inc_mirror=inc_mirror
    )

    angles_no_mirror_cart = sym_class.to_cartesian(sym_class.angles, tolistconv=False)

    nonzero_mask = matplotlib.numpy.nonzero(radius_array)
    radius_array = radius_array[nonzero_mask]

    # Calculate best width and length for the bins in 3D
    width = old_div(
        (pixel_size * particle_radius * matplotlib.numpy.radians(delta) * 2),
        float(2 * 3),
    )
    length = particle_radius * 0.2

    # Normalize the radii so that the maximum value is 1
    radius_normalized = matplotlib.numpy.true_divide(
        radius_array, matplotlib.numpy.max(radius_array)
    )

    # Vector pointing to the center of the Box in chimera
    vector_center = 0.5 * matplotlib.numpy.array([box_size, box_size, box_size])

    # Inner and outer vector. The bin starts at inner vector and stops at outer vector
    inner_vector = vector_center + angles_no_mirror_cart[nonzero_mask] * particle_radius
    outer_vector = vector_center + angles_no_mirror_cart[nonzero_mask] * (
        particle_radius
        + radius_normalized.reshape(radius_normalized.shape[0], 1) * length
    )

    # Adjust pixel size
    matplotlib.numpy.multiply(inner_vector, pixel_size, out=inner_vector)
    matplotlib.numpy.multiply(outer_vector, pixel_size, out=outer_vector)

    # Create output bild file
    if do_print:
        sp_global_def.sxprint("Create bild file")
    output_bild_file = os.path.join(output_folder, "{0}.bild".format(prefix))
    with open(output_bild_file, "w") as write:
        for inner, outer, radius in zip(inner_vector, outer_vector, radius_normalized):
            write.write(".color 0 {0} {1}\n".format(radius, 1 - radius))
            write.write(
                ".cylinder {0:.3f} {1:.3f} {2:.3f} {3:.3f} {4:.3f} {5:.3f} {6:.3f}\n".format(
                    inner[COLUMN_X],
                    inner[COLUMN_Y],
                    inner[COLUMN_Z],
                    outer[COLUMN_X],
                    outer[COLUMN_Y],
                    outer[COLUMN_Z],
                    width,
                )
            )

    """Multiline Comment37"""

    nonzero_mask = list(nonzero_mask[0])
    sorted_radius = radius_array
    sorted_radius_plot = matplotlib.numpy.sort(radius_array)[::-1]
    array_x = matplotlib.numpy.arange(sorted_radius.shape[0])
    # """

    # 2D distribution plot
    if do_print:
        sp_global_def.sxprint("Create 2D legend plot")
    output_bild_legend_png = os.path.join(output_folder, "{0}.png".format(prefix))
    color = get_color(sorted_radius_plot)
    matplotlib.pyplot.bar(array_x, height=sorted_radius_plot, width=1, color=color)
    matplotlib.pyplot.grid()
    matplotlib.pyplot.xlabel("Bin / a.u.")
    matplotlib.pyplot.ylabel("Nr. of Particles")
    matplotlib.pyplot.savefig(output_bild_legend_png, dpi=dpi)
    matplotlib.pyplot.clf()
    """Multiline Comment38"""
    # 2D distribution txt file
    if do_print:
        sp_global_def.sxprint("Create 2D legend text file")

    output_bild_legend_txt = os.path.join(output_folder, "{0}.txt".format(prefix))
    with open(output_bild_legend_txt, "w") as write:
        for i in range(len(sorted_radius)):
            # 	for value_x, value_y in zip(array_x, sorted_radius):
            value_x = "{0:6d}".format(array_x[i])
            value_y = "{0:6d}".format(sorted_radius[i])
            phi = "{0:10f}".format(sym_class.angles[nonzero_mask[i]][0])
            theta = "{0:10f}".format(sym_class.angles[nonzero_mask[i]][1])
            write.write("{0}\n".format("\t".join([value_x, value_y, phi, theta])))


def numpy2em_python(numpy_array):
    """
	Create an EMData object based on a numpy array by reference.
	The output EMData object will have the reversed order of dimensions.
	x,y,z -> z,y,x

	Arguments:
	numpy_array - Array to convert

	Return:
	EMData object
	"""
    shape = numpy_array.shape[::-1]
    if len(shape) == 1:
        shape = (shape[0], 1)
    return_array = EMAN2.EMData(*shape)
    return_view = EMAN2_cppwrap.EMNumPy.em2numpy(return_array)
    return_view[...] = numpy_array
    return_array.update()
    return return_array


def create_summovie_command(temp_name, micrograph_name, shift_name, frc_name, opt):

    # Handle first and last case events
    if opt["first"] == 0:
        sp_global_def.ERROR(
            "SumMovie indexing starts with 1.\n" + "0 is not a valid entry for --first",
            "sxsummovie.py",
            1,
        )
    elif opt["first"] < 0:
        first = opt["nr_frames"] + opt["first"] + 1
    else:
        first = opt["first"]

    if opt["last"] == 0:
        sp_global_def.ERROR(
            "SumMovie indexing starts with 1.\n" + "0 is not a valid entry for --last",
            "sxsummovie.py",
            1,
        )
    elif opt["last"] < 0:
        last = opt["nr_frames"] + opt["last"] + 1
    else:
        last = opt["last"]

    if first > last:
        sp_global_def.ERROR(
            "First option musst be smaller equals last option!\n"
            + "first: {0}; last: {1}".format(first, last),
            "sxsummovie.py",
            1,
        )

    if opt["nr_frames"] < last or last <= 0:
        sp_global_def.ERROR(
            "--last option {0} is out of range:\n".format(last)
            + "min: 1; max {0}".format(opt["nr_frames"]),
            "sxsummovie.py",
            1,
        )

    if opt["nr_frames"] < first or first <= 0:
        sp_global_def.ERROR(
            "--first option {0} is out of range:\n".format(first)
            + "min: 1; max {0}".format(opt["nr_frames"]),
            "sxsummovie.py",
            1,
        )

    # Command list
    summovie_command = []

    # Input file
    summovie_command.append("{0}".format(temp_name))
    # Number of frames
    summovie_command.append("{0}".format(opt["nr_frames"]))
    # Sum file
    summovie_command.append(micrograph_name)
    # Shift file
    summovie_command.append(shift_name)
    # FRC file
    summovie_command.append(frc_name),
    # First frame
    summovie_command.append("{0}".format(first))
    # Last frame
    summovie_command.append("{0}".format(last))
    # Pixel size
    summovie_command.append("{0}".format(opt["pixel_size"]))
    # Dose correction
    if not opt["apply_dose_filter"]:
        summovie_command.append("NO")
    else:
        summovie_command.append("YES")
        # Exposure per frame
        summovie_command.append("{0}".format(opt["exposure_per_frame"]))
        # Acceleration voltage
        summovie_command.append("{0}".format(opt["voltage"]))
        # Pre exposure
        summovie_command.append("{0}".format(opt["pre_exposure"]))
        # Restore noise power
        if opt["dont_restore_noise"]:
            summovie_command.append("NO")
        else:
            summovie_command.append("YES")

    return summovie_command


def gather_compacted_EMData_to_root(number_of_all_em_objects_distributed_across_processes,
                                    list_of_em_objects_for_myid_process, mpi_rank, comm=-1):
    """
    MPI_Gather that expects all data as a list of EMData objects. Data is transfered by
    collecting all header information, transferring all data as numpy arrays, and re-
    building the EMData objects in the root process.

    Args:
        number_of_all_em_objects_distributed_across_processes (integer): Total number of
            EMData objects to gather across all processes.

        list_of_em_objects_for_myid_process (EMData[]): List of EMData objects that this
            process is tasked with transferring.

        mpi_rank (integer): MPI rank of this process. Process with rank 0 will be used as
            the root process.

        comm: MPI communicator identifier; defaults to mpi.MPI_COMM_WORLD
            [Default: -1]
    """
    # mpi setup
    if comm == -1 or comm == None:
        comm = mpi.MPI_COMM_WORLD

    mpi_size = mpi.mpi_comm_size(comm)  # Total number of processes, passed by --np option. [Unclear what this option is]

    tag_for_send_receive = 123456  # NOTE: we're using the same tag for all messages here, which seems a bit pointless?

    # get data header information
    reference_em_object = list_of_em_objects_for_myid_process[0]

    nx = reference_em_object.get_xsize()
    ny = reference_em_object.get_ysize()
    nz = reference_em_object.get_zsize()
    is_ri = reference_em_object.is_ri()
    changecount = reference_em_object.get_attr("changecount")
    is_complex_x = reference_em_object.is_complex_x()
    is_complex_ri = reference_em_object.get_attr("is_complex_ri")
    apix_x = reference_em_object.get_attr("apix_x")
    apix_y = reference_em_object.get_attr("apix_y")
    apix_z = reference_em_object.get_attr("apix_z")
    is_complex = reference_em_object.get_attr_default("is_complex", 1)
    is_fftpad = reference_em_object.get_attr_default("is_fftpad", 1)
    is_fftodd = reference_em_object.get_attr_default("is_fftodd", nz % 2)

    # process data before transfer
    data = []
    for i in range(len(list_of_em_objects_for_myid_process)):
        data.append(EMAN2_cppwrap.EMNumPy.em2numpy(list_of_em_objects_for_myid_process[i]))
    data = numpy.array(data)

    # transfer: gather all data at the root process (myid==0)
    for sender_id in range(1, mpi_size):

        # transfer/send: sender process w/ rank <sender_id> blocks until data has been sent to root process
        if sender_id == mpi_rank:
            mpi.mpi_send(data,
                     data.size,
                     mpi.MPI_FLOAT,
                     0,
                     tag_for_send_receive,
                     mpi.MPI_COMM_WORLD)

        # transfer/receive: root process blocks until data from process w/ rank <sender_id> has been received
        elif mpi_rank == 0:

            # get transfer parameters of the sender
            ref_start, ref_end = sp_applications.MPI_start_end(number_of_all_em_objects_distributed_across_processes, mpi_size,
                                               sender_id)
            transfer_size = (ref_end - ref_start) * data[0].size  # all data[i] have the same size

            # receive data
            data = mpi.mpi_recv(transfer_size,
                            mpi.MPI_FLOAT,
                            sender_id,
                            tag_for_send_receive,
                            mpi.MPI_COMM_WORLD)

            if int(nz) != 1:
                data = data.reshape([(ref_end - ref_start), nz, ny, nx])
            elif ny != 1:
                data = data.reshape([(ref_end - ref_start), ny, nx])

            # collect received data in EMData objects
            for img_data in data:
                em_object = numpy2em_python(img_data)
                em_object.set_ri(is_ri)
                em_object.set_attr_dict({"changecount": changecount,
                                         "is_complex_x": is_complex_x,
                                         "is_complex_ri": is_complex_ri,
                                         "apix_x": apix_x,
                                         "apix_y": apix_y,
                                         "apix_z": apix_z,
                                         "is_complex": is_complex,
                                         "is_fftodd": is_fftodd,
                                         "is_fftpad": is_fftpad})
                list_of_em_objects_for_myid_process.append(em_object)


def params_2D_3D(alpha, sx, sy, mirror):
	"""
		Convert 2D alignment parameters (alpha, sx, sy, mirror) into
		3D alignment parameters (phi, theta, psi, s2x, s2y)
	"""
	phi = 0
	psi = 0
	theta = 0
	alphan, s2x, s2y, scalen = compose_transform2(0, sx, sy, 1, -alpha, 0, 0, 1)
	if mirror > 0:
		phi = (540.0 + phi) % 360.0
		theta = 180.0 - theta
		psi = (540.0 - psi + alphan) % 360.0
	else:
		psi = (psi + alphan) % 360.0
	return phi, theta, psi, s2x, s2y

def compose_transform2(alpha1, sx1, sy1, scale1, alpha2, sx2, sy2, scale2):
	"""Print the composition of two transformations  T2*T1
		Here  if v's are vectors:   vnew = T2*T1 vold
			 with T1 described by alpha1, sx1, scale1 etc.

	  Combined parameters correspond to image first transformed by set 1 followed by set 2.

		Usage: compose_transform2(alpha1,sx1,sy1,scale1,alpha2,sx2,sy2,scale2)
		   angles in degrees
	"""

	t1 = EMAN2_cppwrap.Transform(
		{
			"type": "2D",
			"alpha": alpha1,
			"tx": sx1,
			"ty": sy1,
			"mirror": 0,
			"scale": scale1,
		}
	)
	t2 = EMAN2_cppwrap.Transform(
		{
			"type": "2D",
			"alpha": alpha2,
			"tx": sx2,
			"ty": sy2,
			"mirror": 0,
			"scale": scale2,
		}
	)
	tt = t2 * t1
	d = tt.get_params("2D")
	return d["alpha"], d["tx"], d["ty"], d["scale"]
# ------------------------------------------------[ import ]

# compatibility
from future import standard_library

standard_library.install_aliases()

# python natives
from builtins import range
from builtins import object

# python commons

# EMAN2 / sparx basics


# EMAN2 / sparx modules

# MPI imports (NOTE: import mpi after EMAN2)


statistics_send_recv = dict()
