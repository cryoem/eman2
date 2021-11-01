
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
from . import sp_filter
from . import sp_fundamentals
from . import sp_global_def
from . import sp_utilities


def insert_slices(reconstructor, proj):
    xforms = [proj.get_attr("xform.projection")]
    weights = [proj.get_attr_default("weight", 1.0)]
    ixform = 0
    while True:
        ixform += 1
        xform_proj = proj.get_attr_default("xform.projection" + str(ixform), None)
        if xform_proj == None:
            break
        # putting params in a list does not seem to be necessary, one could call reconstructor as one goes.
        xforms.append(xform_proj)
        # weights.append(proj.get_attr_default("weight" + str(ixform), 1.0))
        weights.append(1.0)
    for i in range(len(xforms)):
        reconstructor.insert_slice(proj, xforms[i], weights[i])


def recons3d_4nn(
    stack_name,
    list_proj=[],
    symmetry="c1",
    npad=4,
    snr=None,
    weighting=1,
    varsnr=False,
    xysize=-1,
    zsize=-1,
):
    """
	Perform a 3-D reconstruction using Pawel's FFT Back Projection algorithm.

	Input:
		stack_name - name of the file with projection data.

		list_proj -  list of projections to be used in the reconstruction

		symmetry - Point group of the target molecule (defaults to "C1")

		npad -

		Angles and shifts are passed in the file header as set_attr. Keywords are phi, theta, psi, sx, sy

		Return:  3D reconstructed volume image

		Usage:
			vol = recons3d_4nn(filepattern, list_proj, symmetry)
	"""

    if list_proj == []:
        if isinstance(stack_name, (bytes, str)):
            nima = EMAN2_cppwrap.EMUtil.get_image_count(stack_name)
        else:
            nima = len(stack_name)
        list_proj = list(range(nima))
    # read first image to determine the size to use
    if isinstance(stack_name, (bytes, str)):
        proj = EMAN2_cppwrap.EMData()
        proj.read_image(stack_name, list_proj[0])
    else:
        proj = stack_name[list_proj[0]].copy()

    size = proj.get_xsize()
    # sanity check -- image must be square
    if size != proj.get_ysize():
        sp_global_def.ERROR("input data has to be square", "recons3d_4nn", 1)

    # reconstructor
    fftvol = EMAN2_cppwrap.EMData()
    weight = EMAN2_cppwrap.EMData()
    params = {
        "npad": npad,
        "symmetry": symmetry,
        "weighting": weighting,
        "fftvol": fftvol,
        "weight": weight,
    }
    if xysize == -1 and zsize == -1:
        params["size"] = size
        if snr != None:
            params["snr"] = snr
            # params["varsnr"] = int(varsnr)
        r = EMAN2_cppwrap.Reconstructors.get("nn4", params)
    else:
        if xysize != -1 and zsize != -1:
            rx = old_div(float(xysize), size)
            ry = old_div(float(xysize), size)
            rz = old_div(float(zsize), size)
        elif xysize != -1:
            rx = old_div(float(xysize), size)
            ry = old_div(float(xysize), size)
            rz = 1.0
        else:
            rx = 1.0
            ry = 1.0
            rz = old_div(float(zsize), size)

        if snr is None:
            params["sizeprojection"] = size
            params["xratio"] = rx
            params["yratio"] = ry
            params["zratio"] = rz
        else:
            params["sizeprojection"] = size
            params["snr"] = snr
            params["varsnr"] = int(varsnr)
            params["xratio"] = rx
            params["yratio"] = ry
            params["zratio"] = rz
        r = EMAN2_cppwrap.Reconstructors.get("nn4_rect", params)

    r.setup()

    if isinstance(stack_name, (bytes, str)):
        for i in range(len(list_proj)):
            proj.read_image(stack_name, list_proj[i])
            # horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
            # active = proj.get_attr_default('active', 1)
            # if active == 1:
            # 	insert_slices(r, proj)
            insert_slices(r, proj)
    else:
        for i in list_proj:
            # horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
            # active = stack_name[i].get_attr_default('active', 1)
            # if active == 1:
            # 	insert_slices(r, stack_name[i])
            insert_slices(r, stack_name[i])

    dummy = r.finish(True)
    return fftvol


def recons3d_4nn_MPI(
    myid,
    prjlist,
    symmetry="c1",
    finfo=None,
    snr=1.0,
    npad=2,
    xysize=-1,
    zsize=-1,
    mpi_comm=None,
):
    if mpi_comm == None:
        mpi_comm = mpi.MPI_COMM_WORLD

    if type(prjlist) == list:
        prjlist = sp_utilities.iterImagesList(prjlist)

    if not prjlist.goToNext():
        sp_global_def.ERROR("empty input list", "recons3d_4nn_MPI", 1)

    imgsize = prjlist.image().get_xsize()
    if prjlist.image().get_ysize() != imgsize:
        imgsize = max(imgsize, prjlist.image().get_ysize())
        dopad = True
    else:
        dopad = False
    prjlist.goToPrev()

    fftvol = EMAN2_cppwrap.EMData()
    weight = EMAN2_cppwrap.EMData()
    if xysize == -1 and zsize == -1:
        params = {
            "size": imgsize,
            "npad": npad,
            "symmetry": symmetry,
            "fftvol": fftvol,
            "weight": weight,
            "snr": snr,
        }
        r = EMAN2_cppwrap.Reconstructors.get("nn4", params)
    else:
        if xysize != -1 and zsize != -1:
            rx = old_div(float(xysize), imgsize)
            ry = old_div(float(xysize), imgsize)
            rz = old_div(float(zsize), imgsize)
        elif xysize != -1:
            rx = old_div(float(xysize), imgsize)
            ry = old_div(float(xysize), imgsize)
            rz = 1.0
        else:
            rx = 1.0
            ry = 1.0
            rz = old_div(float(zsize), imgsize)
        params = {
            "sizeprojection": imgsize,
            "npad": npad,
            "symmetry": symmetry,
            "fftvol": fftvol,
            "weight": weight,
            "xratio": rx,
            "yratio": ry,
            "zratio": rz,
        }
        r = EMAN2_cppwrap.Reconstructors.get("nn4_rect", params)
    r.setup()

    if not (finfo is None):
        nimg = 0
    while prjlist.goToNext():
        prj = prjlist.image()
        if dopad:
            prj = sp_utilities.pad(prj, imgsize, imgsize, 1, "circumference")
        insert_slices(r, prj)
        if not (finfo is None):
            nimg += 1
            finfo.write("Image %4d inserted.\n" % (nimg))
            finfo.flush()

    if not (finfo is None):
        finfo.write("Begin reducing ...\n")
        finfo.flush()

    sp_utilities.reduce_EMData_to_root(fftvol, myid, comm=mpi_comm)
    sp_utilities.reduce_EMData_to_root(weight, myid, comm=mpi_comm)

    if myid == 0:
        dummy = r.finish(True)
    else:
        if xysize == -1 and zsize == -1:
            fftvol = sp_utilities.model_blank(imgsize, imgsize, imgsize)
        else:
            if zsize == -1:
                fftvol = sp_utilities.model_blank(xysize, xysize, imgsize)
            elif xysize == -1:
                fftvol = sp_utilities.model_blank(imgsize, imgsize, zsize)
            else:
                fftvol = sp_utilities.model_blank(xysize, xysize, zsize)
    return fftvol


"""Multiline Comment0"""

"""Multiline Comment1"""

"""Multiline Comment2"""

"""Multiline Comment3"""


def recons3d_4nnw_MPI(
    myid,
    prjlist,
    bckgdata,
    snr=1.0,
    sign=1,
    symmetry="c1",
    finfo=None,
    npad=2,
    xysize=-1,
    zsize=-1,
    mpi_comm=None,
    smearstep=0.0,
    fsc=None,
):
    """
		recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
		Input
			stack: name of the stack file containing projection data, projections have to be squares
			prjlist: list of projections to be included in the reconstruction or image iterator
			bckgdata = [get_im("tsd.hdf"),read_text_file("data_stamp.txt")]
			snr: Signal-to-Noise Ratio of the data
			sign: sign of the CTF
			symmetry: point-group symmetry to be enforced, each projection will enter the reconstruction in all symmetry-related directions.
	"""
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_utilities  import reduce_EMData_to_root, pad
    pass  # IMPORTIMPORTIMPORT from EMAN2      import Reconstructors
    pass  # IMPORTIMPORTIMPORT from ..libpy.sp_utilities  import iterImagesList, set_params_proj, model_blank
    pass  # IMPORTIMPORTIMPORT from mpi        import MPI_COMM_WORLD
    pass  # IMPORTIMPORTIMPORT import types

    if mpi_comm == None:
        mpi_comm = mpi.MPI_COMM_WORLD

    if type(prjlist) == list:
        prjlist = sp_utilities.iterImagesList(prjlist)
    if not prjlist.goToNext():
        sp_global_def.ERROR("empty input list", "recons3d_4nnw_MPI", 1)
    imgsize = prjlist.image().get_xsize()
    if prjlist.image().get_ysize() != imgsize:
        imgsize = max(imgsize, prjlist.image().get_ysize())
        dopad = True
    else:
        dopad = False
    prjlist.goToPrev()

    #  Do the FSC shtick.
    bnx = old_div(imgsize * npad, 2) + 1
    if fsc:
        pass  # IMPORTIMPORTIMPORT from math import sqrt
        pass  # IMPORTIMPORTIMPORT from ..libpy.sp_utilities import reshape_1d
        t = [0.0] * len(fsc)
        for i in range(len(fsc)):
            t[i] = min(max(fsc[i], 0.0), 0.999)
        t = sp_utilities.reshape_1d(t, len(t), npad * len(t))
        refvol = sp_utilities.model_blank(bnx, 1, 1, 0.0)
        for i in range(len(fsc)):
            refvol.set_value_at(i, t[i])
    else:
        refvol = sp_utilities.model_blank(bnx, 1, 1, 1.0)
    refvol.set_attr("fudge", 1.0)

    fftvol = EMAN2_cppwrap.EMData()
    weight = EMAN2_cppwrap.EMData()

    if smearstep > 0.0:
        # if myid == 0:  print "  Setting smear in prepare_recons_ctf"
        ns = 1
        smear = []
        for j in range(-ns, ns + 1):
            if j != 0:
                for i in range(-ns, ns + 1):
                    for k in range(-ns, ns + 1):
                        smear += [i * smearstep, j * smearstep, k * smearstep, 1.0]
        # Deal with theta = 0.0 cases
        prj = []
        for i in range(-ns, ns + 1):
            for k in range(-ns, ns + 1):
                prj.append(i + k)
        for i in range(-2 * ns, 2 * ns + 1, 1):
            smear += [i * smearstep, 0.0, 0.0, float(prj.count(i))]
        # if myid == 0:  print "  Smear  ",smear
        fftvol.set_attr("smear", smear)

    if xysize == -1 and zsize == -1:
        params = {
            "size": imgsize,
            "npad": npad,
            "snr": snr,
            "sign": sign,
            "symmetry": symmetry,
            "refvol": refvol,
            "fftvol": fftvol,
            "weight": weight,
        }
        r = EMAN2_cppwrap.Reconstructors.get("nn4_ctfw", params)
    else:
        if xysize != -1 and zsize != -1:
            rx = old_div(float(xysize), imgsize)
            ry = old_div(float(xysize), imgsize)
            rz = old_div(float(zsize), imgsize)
        elif xysize != -1:
            rx = old_div(float(xysize), imgsize)
            ry = old_div(float(xysize), imgsize)
            rz = 1.0
        else:
            rx = 1.0
            ry = 1.0
            rz = old_div(float(zsize), imgsize)
        #  There is an error here with sizeprojection  PAP 10/22/2014
        params = {
            "size": sizeprojection,
            "npad": npad,
            "snr": snr,
            "sign": sign,
            "symmetry": symmetry,
            "fftvol": fftvol,
            "weight": weight,
            "xratio": rx,
            "yratio": ry,
            "zratio": rz,
        }
        r = EMAN2_cppwrap.Reconstructors.get("nn4_ctf_rect", params)
    r.setup()

    # from utilities import model_blank, get_im, read_text_file
    # bckgdata = [get_im("tsd.hdf"),read_text_file("data_stamp.txt")]

    nnx = bckgdata[0].get_xsize()
    nny = bckgdata[0].get_ysize()
    bckgnoise = []
    for i in range(nny):
        prj = sp_utilities.model_blank(nnx)
        for k in range(nnx):
            prj[k] = bckgdata[0].get_value_at(k, i)
        bckgnoise.append(prj)

    datastamp = bckgdata[1]
    if not (finfo is None):
        nimg = 0
    while prjlist.goToNext():
        prj = prjlist.image()
        try:
            stmp = old_div(nnx, 0)
            stmp = prj.get_attr("ptcl_source_image")
        except:
            try:
                stmp = prj.get_attr("ctf")
                stmp = round(stmp.defocus, 4)
            except:
                sp_global_def.ERROR(
                    "Either ptcl_source_image or ctf has to be present in the header.",
                    "recons3d_4nnw_MPI",
                    1,
                    myid,
                )
        try:
            indx = datastamp.index(stmp)
        except:
            sp_global_def.ERROR(
                "Problem with indexing ptcl_source_image.", "recons3d_4nnw_MPI", 1, myid
            )

        if dopad:
            prj = sp_utilities.pad(prj, imgsize, imgsize, 1, "circumference")

        prj.set_attr("bckgnoise", bckgnoise[indx])
        insert_slices(r, prj)
        if not (finfo is None):
            nimg += 1
            finfo.write(" %4d inserted\n" % (nimg))
            finfo.flush()
    del sp_utilities.pad
    if not (finfo is None):
        finfo.write("begin reduce\n")
        finfo.flush()

    sp_utilities.reduce_EMData_to_root(fftvol, myid, comm=mpi_comm)
    sp_utilities.reduce_EMData_to_root(weight, myid, comm=mpi_comm)

    if not (finfo is None):
        finfo.write("after reduce\n")
        finfo.flush()

    if myid == 0:
        dummy = r.finish(True)
    else:
        pass  # IMPORTIMPORTIMPORT from ..libpy.sp_utilities import model_blank
        if xysize == -1 and zsize == -1:
            fftvol = sp_utilities.model_blank(imgsize, imgsize, imgsize)
        else:
            if zsize == -1:
                fftvol = sp_utilities.model_blank(xysize, xysize, imgsize)
            elif xysize == -1:
                fftvol = sp_utilities.model_blank(imgsize, imgsize, zsize)
            else:
                fftvol = sp_utilities.model_blank(xysize, xysize, zsize)
    return fftvol


def recons3d_trl_struct_MPI(
    myid,
    main_node,
    prjlist,
    paramstructure,
    refang,
    rshifts_shrank,
    delta,
    upweighted=True,
    mpi_comm=None,
    CTF=True,
    target_size=-1,
    avgnorm=1.0,
    norm_per_particle=None,
):
    """
		recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
		Input
			list_of_prjlist: list of lists of projections to be included in the reconstruction
	"""

    if mpi_comm == None:
        mpi_comm = mpi.MPI_COMM_WORLD

    refvol = sp_utilities.model_blank(target_size)
    refvol.set_attr("fudge", 1.0)

    if CTF:
        do_ctf = 1
    else:
        do_ctf = 0

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

    if prjlist:
        if norm_per_particle == None:
            norm_per_particle = len(prjlist) * [1.0]

        nnx = prjlist[0].get_xsize()
        nny = prjlist[0].get_ysize()
        nshifts = len(rshifts_shrank)
        for im in range(len(prjlist)):
            #  parse projection structure, generate three lists:
            #  [ipsi+iang], [ishift], [probability]
            #  Number of orientations for a given image
            numbor = len(paramstructure[im][2])
            ipsiandiang = [
                old_div(paramstructure[im][2][i][0], 1000) for i in range(numbor)
            ]
            allshifts = [paramstructure[im][2][i][0] % 1000 for i in range(numbor)]
            probs = [paramstructure[im][2][i][1] for i in range(numbor)]
            #  Find unique projection directions
            tdir = list(set(ipsiandiang))
            bckgn = prjlist[im].get_attr("bckgnoise")
            ct = prjlist[im].get_attr("ctf")
            #  For each unique projection direction:
            data = [None] * nshifts
            for ii in range(len(tdir)):
                #  Find the number of times given projection direction appears on the list, it is the number of different shifts associated with it.
                lshifts = sp_utilities.findall(tdir[ii], ipsiandiang)
                toprab = 0.0
                for ki in range(len(lshifts)):
                    toprab += probs[lshifts[ki]]
                recdata = EMAN2_cppwrap.EMData(nny, nny, 1, False)
                recdata.set_attr("is_complex", 0)
                for ki in range(len(lshifts)):
                    lpt = allshifts[lshifts[ki]]
                    if data[lpt] == None:
                        data[lpt] = sp_fundamentals.fshift(
                            prjlist[im], rshifts_shrank[lpt][0], rshifts_shrank[lpt][1]
                        )
                        data[lpt].set_attr("is_complex", 0)
                    EMAN2_cppwrap.Util.add_img(
                        recdata,
                        EMAN2_cppwrap.Util.mult_scalar(
                            data[lpt], old_div(probs[lshifts[ki]], toprab)
                        ),
                    )
                recdata.set_attr_dict(
                    {
                        "padffted": 1,
                        "is_fftpad": 1,
                        "is_fftodd": 0,
                        "is_complex_ri": 1,
                        "is_complex": 1,
                    }
                )
                if not upweighted:
                    recdata = sp_filter.filt_table(recdata, bckgn)
                recdata.set_attr_dict({"bckgnoise": bckgn, "ctf": ct})
                ipsi = tdir[ii] % 100000
                iang = old_div(tdir[ii], 100000)
                r.insert_slice(
                    recdata,
                    EMAN2_cppwrap.Transform(
                        {
                            "type": "spider",
                            "phi": refang[iang][0],
                            "theta": refang[iang][1],
                            "psi": refang[iang][2] + ipsi * delta,
                        }
                    ),
                    old_div(toprab * avgnorm, norm_per_particle[im]),
                )
        #  clean stuff
        del bckgn, recdata, tdir, ipsiandiang, allshifts, probs

    sp_utilities.reduce_EMData_to_root(fftvol, myid, main_node, comm=mpi_comm)
    sp_utilities.reduce_EMData_to_root(weight, myid, main_node, comm=mpi_comm)

    if myid == main_node:
        dummy = r.finish(True)
    mpi.mpi_barrier(mpi_comm)

    if myid == main_node:
        return fftvol, weight, refvol
    else:
        return None, None, None


def recons3d_4nn_ctf(
    stack_name,
    list_proj=[],
    snr=1.0,
    sign=1,
    symmetry="c1",
    verbose=0,
    npad=2,
    xysize=-1,
    zsize=-1,
):
    """Perform a 3-D reconstruction using Pawel's FFT Back Projection algoritm.

	   Input:
	    stack_name - name of the stack file on a disk,
	                 each image has to have the following attributes set:
			 psi, theta, phi, sx, sy, defocus,
	    list_proj - list of images from stack_name to be included in the reconstruction
	    symmetry	 -- Point group of the target molecule (defaults to "C1")

	   Return:  3d reconstructed volume image

	   Usage:

	     anglelist = getAngles("myangles.txt") # not yet written
	     vol = do_reconstruction(filepattern, start, end, anglelist, symmetry)
	"""

    # read first image to determine the size to use
    if list_proj == []:
        if isinstance(stack_name, (bytes,str)) :
            nima = EMAN2_cppwrap.EMUtil.get_image_count(stack_name)
        else:
            nima = len(stack_name)
        list_proj = list(range(nima))
    # read first image to determine the size to use
    if isinstance(stack_name, (bytes,str)):
        proj = EMAN2_cppwrap.EMData()
        proj.read_image(stack_name, list_proj[0])
    else:
        proj = stack_name[list_proj[0]].copy()

    # convert angles to transform (rotation) objects
    # horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
    # active = proj.get_attr_default('active', 1)
    size = proj.get_xsize()
    if proj.get_ysize() != size:
        size = max(size, proj.get_ysize())
        dopad = True
    else:
        dopad = False

    # reconstructor
    fftvol = EMAN2_cppwrap.EMData()
    weight = EMAN2_cppwrap.EMData()
    params = {
        "npad": npad,
        "symmetry": symmetry,
        "snr": snr,
        "sign": sign,
        "fftvol": fftvol,
        "weight": weight,
    }
    if xysize == -1 and zsize == -1:
        params["size"] = size
        r = EMAN2_cppwrap.Reconstructors.get("nn4_ctf", params)
    else:
        if xysize != -1 and zsize != -1:
            rx = old_div(float(xysize), size)
            ry = old_div(float(xysize), size)
            rz = old_div(float(zsize), size)
        elif xysize != -1:
            rx = old_div(float(xysize), size)
            ry = old_div(float(xysize), size)
            rz = 1.0
        else:
            rx = 1.0
            ry = 1.0
            rz = old_div(float(zsize), size)

        params["sizeprojection"] = size
        params["xratio"] = rx
        params["yratio"] = ry
        params["zratio"] = rz
        r = EMAN2_cppwrap.Reconstructors.get("nn4_ctf_rect", params)
    r.setup()

    if isinstance(stack_name, (bytes,str)) :
        for i in range(len(list_proj)):
            proj.read_image(stack_name, list_proj[i])
            if dopad:
                proj = sp_utilities.pad(proj, size, size, 1, "circumference")
            insert_slices(r, proj)
    else:
        for i in range(len(list_proj)):
            insert_slices(r, stack_name[list_proj[i]])
    dummy = r.finish(True)
    return fftvol


def recons3d_4nn_ctf_MPI(
    myid,
    prjlist,
    snr=1.0,
    sign=1,
    symmetry="c1",
    finfo=None,
    npad=2,
    xysize=-1,
    zsize=-1,
    mpi_comm=None,
    smearstep=0.0,
):
    """
		recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
		Input
			stack: name of the stack file containing projection data, projections have to be squares
			list_proj: list of projections to be included in the reconstruction or image iterator
			snr: Signal-to-Noise Ratio of the data
			sign: sign of the CTF
			symmetry: point-group symmetry to be enforced, each projection will enter the reconstruction in all symmetry-related directions.
	"""

    if mpi_comm == None:
        mpi_comm = mpi.MPI_COMM_WORLD

    if type(prjlist) == list:
        prjlist = sp_utilities.iterImagesList(prjlist)
    if not prjlist.goToNext():
        sp_global_def.ERROR("empty input list", "recons3d_4nn_ctf_MPI", 1)
    imgsize = prjlist.image().get_xsize()
    if prjlist.image().get_ysize() != imgsize:
        imgsize = max(imgsize, prjlist.image().get_ysize())
        dopad = True
    else:
        dopad = False
    prjlist.goToPrev()

    fftvol = EMAN2_cppwrap.EMData()

    if smearstep > 0.0:
        # if myid == 0:  print "  Setting smear in prepare_recons_ctf"
        ns = 1
        smear = []
        for j in range(-ns, ns + 1):
            if j != 0:
                for i in range(-ns, ns + 1):
                    for k in range(-ns, ns + 1):
                        smear += [i * smearstep, j * smearstep, k * smearstep, 1.0]
        # Deal with theta = 0.0 cases
        prj = []
        for i in range(-ns, ns + 1):
            for k in range(-ns, ns + 1):
                prj.append(i + k)
        for i in range(-2 * ns, 2 * ns + 1, 1):
            smear += [i * smearstep, 0.0, 0.0, float(prj.count(i))]
        # if myid == 0:  print "  Smear  ",smear
        fftvol.set_attr("smear", smear)

    weight = EMAN2_cppwrap.EMData()
    if xysize == -1 and zsize == -1:
        params = {
            "size": imgsize,
            "npad": npad,
            "snr": snr,
            "sign": sign,
            "symmetry": symmetry,
            "fftvol": fftvol,
            "weight": weight,
        }
        r = EMAN2_cppwrap.Reconstructors.get("nn4_ctf", params)
    else:
        if xysize != -1 and zsize != -1:
            rx = old_div(float(xysize), imgsize)
            ry = old_div(float(xysize), imgsize)
            rz = old_div(float(zsize), imgsize)
        elif xysize != -1:
            rx = old_div(float(xysize), imgsize)
            ry = old_div(float(xysize), imgsize)
            rz = 1.0
        else:
            rx = 1.0
            ry = 1.0
            rz = old_div(float(zsize), imgsize)
        #  There is an error here with sizeprojection  PAP 10/22/2014
        params = {
            "size": sizeprojection,
            "npad": npad,
            "snr": snr,
            "sign": sign,
            "symmetry": symmetry,
            "fftvol": fftvol,
            "weight": weight,
            "xratio": rx,
            "yratio": ry,
            "zratio": rz,
        }
        r = EMAN2_cppwrap.Reconstructors.get("nn4_ctf_rect", params)
    r.setup()

    # if not (finfo is None):
    nimg = 0
    while prjlist.goToNext():
        prj = prjlist.image()
        if dopad:
            prj = sp_utilities.pad(prj, imgsize, imgsize, 1, "circumference")
        # if params:
        insert_slices(r, prj)
        if not (finfo is None):
            nimg += 1
            finfo.write(" %4d inserted\n" % (nimg))
            finfo.flush()
    del sp_utilities.pad
    if not (finfo is None):
        finfo.write("begin reduce\n")
        finfo.flush()

    sp_utilities.reduce_EMData_to_root(fftvol, myid, comm=mpi_comm)
    sp_utilities.reduce_EMData_to_root(weight, myid, comm=mpi_comm)

    if not (finfo is None):
        finfo.write("after reduce\n")
        finfo.flush()

    if myid == 0:
        dummy = r.finish(True)
    else:
        if xysize == -1 and zsize == -1:
            fftvol = sp_utilities.model_blank(imgsize, imgsize, imgsize)
        else:
            if zsize == -1:
                fftvol = sp_utilities.model_blank(xysize, xysize, imgsize)
            elif xysize == -1:
                fftvol = sp_utilities.model_blank(imgsize, imgsize, zsize)
            else:
                fftvol = sp_utilities.model_blank(xysize, xysize, zsize)
    return fftvol




def get_image_size( imgdata, myid ):
	from mpi import mpi_gather, mpi_bcast, MPI_COMM_WORLD, MPI_INT
	nimg = len(imgdata)

	nimgs = mpi_gather( nimg, 1, MPI_INT, 1, MPI_INT, 0, MPI_COMM_WORLD )

	if myid==0:
		src = -1
		for i in range( len(nimgs) ):
			if int(nimgs[i]) > 0 :
				src = i
				break
		if src==-1:
			return 0
	else:
		src = -1

	size_src = mpi_bcast( src, 1, MPI_INT, 0, MPI_COMM_WORLD )

	if myid==int(size_src[0]):
		assert nimg > 0
		size = imgdata[0].get_xsize()
	else:
		size = -1

	nx = mpi_bcast( size, 1, MPI_INT, size_src[0], MPI_COMM_WORLD )
	return int(nx[0])

def prepare_recons_ctf(nx, data, snr, symmetry, myid, main_node_half, half_start, step, finfo=None, npad = 2, mpi_comm=None, smearstep = 0.0):
	from random     import randint
	from sphire.libpy.sp_utilities  import reduce_EMData_to_root
	from mpi        import mpi_barrier, MPI_COMM_WORLD
	from EMAN2 import EMData, Reconstructors

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	fftvol_half = EMData()

	if( smearstep > 0.0 ):
		#if myid == 0:  print "  Setting smear in prepare_recons_ctf"
		ns = 1
		smear = []
		for j in range(-ns,ns+1):
			if( j != 0):
				for i in range(-ns,ns+1):
					for k in range(-ns,ns+1):
						smear += [i*smearstep,j*smearstep,k*smearstep,1.0]
		# Deal with theta = 0.0 cases
		prj = []
		for i in range(-ns,ns+1):
			for k in range(-ns,ns+1):
				prj.append(i+k)
		for i in range(-2*ns,2*ns+1,1):
			smear += [i*smearstep,0.0,0.0,float(prj.count(i))]
		#if myid == 0:  print "  Smear  ",smear
		fftvol_half.set_attr("smear", smear)

	weight_half = EMData()
	half_params = {"size":nx, "npad":npad, "snr":snr, "sign":1, "symmetry":symmetry, "fftvol":fftvol_half, "weight":weight_half}
	half = Reconstructors.get( "nn4_ctf", half_params )
	half.setup()

	for i in range(half_start, len(data), step):
		xform_proj = data[i].get_attr( "xform.projection" )
		half.insert_slice(data[i], xform_proj )

	if not(finfo is None):
		finfo.write( "begin reduce half\n" )
		finfo.flush()

	reduce_EMData_to_root(fftvol_half, myid, main_node_half, mpi_comm)
	reduce_EMData_to_root(weight_half, myid, main_node_half, mpi_comm)

	if not(finfo is None):
		finfo.write( "after reduce half\n" )
		finfo.flush()

	if myid == main_node_half:
		tmpid = randint(0, 1000000)
		fftvol_half_file = ("fftvol_half%d.hdf" % tmpid)
		weight_half_file = ("weight_half%d.hdf" % tmpid)
		fftvol_half.write_image(fftvol_half_file)
		weight_half.write_image(weight_half_file)
	mpi_barrier(mpi_comm)

	fftvol_half = None
	weight_half = None

	if myid == main_node_half:
		return fftvol_half_file, weight_half_file

	return None,None

def recons_ctf_from_fftvol(size, fftvol, weight, snr, symmetry, weighting=1, npad = 2):
	from EMAN2 import Reconstructors

	params = {"size":size, "npad":npad, "snr":snr, "sign":1, "symmetry":symmetry, "fftvol":fftvol, "weight":weight, "weighting":weighting}
	r = Reconstructors.get("nn4_ctf", params)
	r.setup()
	dummy = r.finish(True)
	return fftvol


def rec3D_MPI(data, snr=1.0, symmetry="c1", mask3D=None, fsc_curve=None, \
              myid=0, main_node=0, rstep=1.0, odd_start=0, eve_start=1, finfo=None, \
              index=-1, npad=2, mpi_comm=None, smearstep=0.0):
    '''
	  This function is to be called within an MPI program to do a reconstruction on a dataset kept
	  in the memory, computes reconstruction and through odd-even, in order to get the resolution
	'''
    import os
    from sphire.libpy.sp_statistics import fsc_mask
    from sphire.libpy.sp_utilities import model_blank, model_circle, get_image, send_EMData, recv_EMData
    from mpi import mpi_comm_size, MPI_COMM_WORLD
    from EMAN2 import Util

    if mpi_comm == None:
        mpi_comm = MPI_COMM_WORLD

    nproc = mpi_comm_size(mpi_comm)

    if nproc == 1:
        assert main_node == 0
        main_node_odd = main_node
        main_node_eve = main_node
        main_node_all = main_node
    elif nproc == 2:
        main_node_odd = main_node
        main_node_eve = (main_node + 1) % 2
        main_node_all = main_node

        tag_voleve = 1000
        tag_fftvol_eve = 1001
        tag_weight_eve = 1002
    else:
        # spread CPUs between different nodes to save memory
        main_node_odd = main_node
        main_node_eve = (int(main_node) + nproc - 1) % int(nproc)
        main_node_all = (int(main_node) + nproc // 2) % int(nproc)

        tag_voleve = 1000
        tag_fftvol_eve = 1001
        tag_weight_eve = 1002

        tag_fftvol_odd = 1003
        tag_weight_odd = 1004
        tag_volall = 1005

    if index != -1:
        grpdata = []
        for i in range(len(data)):
            if data[i].get_attr('group') == index:
                grpdata.append(data[i])
        imgdata = grpdata
    else:
        imgdata = data

    nx = get_image_size(imgdata, myid)
    if nx == 0:
        ERROR(
            "Warning: no images were given for reconstruction, this usually means there is an empty group, returning empty volume",
            "rec3D", 0)
        return model_blank(2, 2, 2), None

    fftvol_odd_file, weight_odd_file = prepare_recons_ctf(nx, imgdata, snr, symmetry, myid, main_node_odd, odd_start, 2,
                                                          finfo, npad, mpi_comm=mpi_comm, smearstep=smearstep)
    fftvol_eve_file, weight_eve_file = prepare_recons_ctf(nx, imgdata, snr, symmetry, myid, main_node_eve, eve_start, 2,
                                                          finfo, npad, mpi_comm=mpi_comm, smearstep=smearstep)
    del imgdata

    if nproc == 1:
        fftvol = get_image(fftvol_odd_file)
        weight = get_image(weight_odd_file)
        volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad=npad)

        fftvol = get_image(fftvol_eve_file)
        weight = get_image(weight_eve_file)
        voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad=npad)

        if (not mask3D):
            nx = volodd.get_xsize()
            ny = volodd.get_ysize()
            nz = volodd.get_zsize()
            mask3D = model_circle(min(nx, ny, nz) // 2 - 2, nx, ny, nz)
        fscdat = fsc_mask(volodd, voleve, mask3D, rstep, fsc_curve)
        del volodd, voleve, mask3d

        fftvol = get_image(fftvol_odd_file)
        fftvol_tmp = get_image(fftvol_eve_file)
        fftvol += fftvol_tmp
        fftvol_tmp = None

        weight = get_image(weight_odd_file)
        weight_tmp = get_image(weight_eve_file)
        weight += weight_tmp
        weight_tmp = None

        volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad=npad)
        os.system("rm -f " + fftvol_odd_file + " " + weight_odd_file)
        os.system("rm -f " + fftvol_eve_file + " " + weight_eve_file)

        return volall, fscdat

    if nproc == 2:
        if myid == main_node_odd:
            fftvol = get_image(fftvol_odd_file)
            weight = get_image(weight_odd_file)
            volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad=npad)
            voleve = recv_EMData(main_node_eve, tag_voleve, mpi_comm)

            if (not mask3D):
                nx = volodd.get_xsize()
                ny = volodd.get_ysize()
                nz = volodd.get_zsize()
                mask3D = model_circle(min(nx, ny, nz) // 2 - 2, nx, ny, nz)
            fscdat = fsc_mask(volodd, voleve, mask3D, rstep, fsc_curve)
            del volodd, voleve, mask3D
        else:
            assert myid == main_node_eve
            fftvol = get_image(fftvol_eve_file)
            weight = get_image(weight_eve_file)
            voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad=npad)
            send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)

        if myid == main_node_odd:
            fftvol = get_image(fftvol_odd_file)
            fftvol_tmp = recv_EMData(main_node_eve, tag_fftvol_eve, mpi_comm)
            fftvol += fftvol_tmp
            fftvol_tmp = None

            weight = get_image(weight_odd_file)
            weight_tmp = recv_EMData(main_node_eve, tag_weight_eve, mpi_comm)
            weight += weight_tmp
            weight_tmp = None

            volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad=npad)
            os.system("rm -f " + fftvol_odd_file + " " + weight_odd_file)

            return volall, fscdat
        else:
            assert myid == main_node_eve
            fftvol = get_image(fftvol_eve_file)
            weight = get_image(weight_eve_file)
            send_EMData(fftvol, main_node_odd, tag_fftvol_eve, mpi_comm)
            send_EMData(weight, main_node_odd, tag_weight_eve, mpi_comm)
            os.system("rm -f " + fftvol_eve_file + " " + weight_eve_file)
            return model_blank(nx, nx, nx), None

    # cases from all other number of processors situations
    if myid == main_node_odd:
        fftvol = get_image(fftvol_odd_file)
        send_EMData(fftvol, main_node_eve, tag_fftvol_odd, mpi_comm)

        if not (finfo is None):
            finfo.write("fftvol odd sent\n")
            finfo.flush()

        weight = get_image(weight_odd_file)
        send_EMData(weight, main_node_all, tag_weight_odd, mpi_comm)

        if not (finfo is None):
            finfo.write("weight odd sent\n")
            finfo.flush()

        volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad=npad)
        del fftvol, weight
        voleve = recv_EMData(main_node_eve, tag_voleve, mpi_comm)

        if (not mask3D):
            nx = volodd.get_xsize()
            ny = volodd.get_ysize()
            nz = volodd.get_zsize()
            mask3D = model_circle(min(nx, ny, nz) // 2 - 2, nx, ny, nz)

        fscdat = fsc_mask(volodd, voleve, mask3D, rstep, fsc_curve)
        del volodd, voleve, mask3D
        volall = recv_EMData(main_node_all, tag_volall, mpi_comm)
        os.system("rm -f " + fftvol_odd_file + " " + weight_odd_file)
        return volall, fscdat

    if myid == main_node_eve:
        ftmp = recv_EMData(main_node_odd, tag_fftvol_odd, mpi_comm)
        fftvol = get_image(fftvol_eve_file)
        Util.add_img(ftmp, fftvol)
        send_EMData(ftmp, main_node_all, tag_fftvol_eve, mpi_comm)
        del ftmp

        weight = get_image(weight_eve_file)
        send_EMData(weight, main_node_all, tag_weight_eve, mpi_comm)

        voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad=npad)
        send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)
        os.system("rm -f " + fftvol_eve_file + " " + weight_eve_file);

        return model_blank(nx, nx, nx), None

    if myid == main_node_all:
        fftvol = recv_EMData(main_node_eve, tag_fftvol_eve, mpi_comm)
        if not (finfo is None):
            finfo.write("fftvol odd received\n")
            finfo.flush()

        weight = recv_EMData(main_node_odd, tag_weight_odd, mpi_comm)
        weight_tmp = recv_EMData(main_node_eve, tag_weight_eve, mpi_comm)
        Util.add_img(weight, weight_tmp)
        weight_tmp = None

        volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad=npad)
        send_EMData(volall, main_node_odd, tag_volall, mpi_comm)

        return model_blank(nx, nx, nx), None

    return model_blank(nx, nx, nx), None



def prepare_recons(data, symmetry, myid, main_node_half, half_start, step, index, finfo=None, npad = 2, mpi_comm=None):
	from random     import randint
	from sphire.libpy.sp_utilities  import reduce_EMData_to_root
	from mpi        import mpi_barrier, MPI_COMM_WORLD
	from EMAN2 import EMData, Reconstructors

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	nx = data[0].get_xsize()

	fftvol_half = EMData()
	weight_half = EMData()
	half_params = {"size":nx, "npad":npad, "symmetry":symmetry, "fftvol":fftvol_half, "weight":weight_half}
	half = Reconstructors.get( "nn4", half_params )
	half.setup()

	group = -1
	for i in range(half_start, len(data), step):
		if(index >-1 ):  group = data[i].get_attr('group')
		if(group == index):
			# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
			# if( data[i].get_attr_default('active',1) == 1):
			# 	xform_proj = data[i].get_attr( "xform.projection" )
			# 	half.insert_slice(data[i], xform_proj )
			xform_proj = data[i].get_attr( "xform.projection" )
			half.insert_slice(data[i], xform_proj )

	if not(finfo is None):
		finfo.write( "begin reduce half\n" )
		finfo.flush()

	reduce_EMData_to_root(fftvol_half, myid, main_node_half, mpi_comm)
	reduce_EMData_to_root(weight_half, myid, main_node_half, mpi_comm)

	if not(finfo is None):
		finfo.write( "after reduce half\n" )
		finfo.flush()

	if myid == main_node_half:
		tmpid = randint(0, 1000000)
		fftvol_half_file = ("fftvol_half%d.hdf" % tmpid)
		weight_half_file = ("weight_half%d.hdf" % tmpid)
		fftvol_half.write_image(fftvol_half_file)
		weight_half.write_image(weight_half_file)
	mpi_barrier(mpi_comm)

	fftvol_half = None
	weight_half = None

	if myid == main_node_half:  return fftvol_half_file, weight_half_file

	return None, None


def rec3D_MPI_noCTF(data, symmetry="c1", mask3D=None, fsc_curve=None, myid=2, main_node=0, \
                    rstep=1.0, odd_start=0, eve_start=1, finfo=None, index=-1, npad=2, mpi_comm=None):
    '''
	  This function is to be called within an MPI program to do a reconstruction on a dataset kept in the memory
	  Computes reconstruction and through odd-even, in order to get the resolution
	  if index > -1, projections should have attribute group set and only those whose group matches index will be used in the reconstruction
		this is for multireference alignment
	'''
    import os
    from sphire.libpy.sp_statistics import fsc_mask
    from sphire.libpy.sp_utilities import model_blank, get_image, send_EMData, recv_EMData
    from mpi import mpi_comm_size, MPI_COMM_WORLD
    from EMAN2 import Util

    if mpi_comm == None:
        mpi_comm = MPI_COMM_WORLD

    nproc = mpi_comm_size(mpi_comm)

    if nproc == 1:
        assert main_node == 0
        main_node_odd = main_node
        main_node_eve = main_node
        main_node_all = main_node
    elif nproc == 2:
        main_node_odd = main_node
        main_node_eve = (main_node + 1) % 2
        main_node_all = main_node

        tag_voleve = 1000
        tag_fftvol_eve = 1001
        tag_weight_eve = 1002
    else:
        # spread CPUs between different nodes to save memory
        main_node_odd = main_node
        main_node_eve = (int(main_node) + nproc - 1) % int(nproc)
        main_node_all = (int(main_node) + nproc // 2) % int(nproc)

        tag_voleve = 1000
        tag_fftvol_eve = 1001
        tag_weight_eve = 1002

        tag_fftvol_odd = 1003
        tag_weight_odd = 1004
        tag_volall = 1005

    nx = data[0].get_xsize()

    fftvol_odd_file, weight_odd_file = prepare_recons(data, symmetry, myid, main_node_odd, odd_start, 2, index, finfo,
                                                      npad, mpi_comm=mpi_comm)
    fftvol_eve_file, weight_eve_file = prepare_recons(data, symmetry, myid, main_node_eve, eve_start, 2, index, finfo,
                                                      npad, mpi_comm=mpi_comm)

    if nproc == 1:
        fftvol = get_image(fftvol_odd_file)
        weight = get_image(weight_odd_file)
        volodd = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)

        fftvol = get_image(fftvol_eve_file)
        weight = get_image(weight_eve_file)
        voleve = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)

        fscdat = fsc_mask(volodd, voleve, mask3D, rstep, fsc_curve)
        del volodd, voleve

        fftvol = get_image(fftvol_odd_file)
        Util.add_img(fftvol, get_image(fftvol_eve_file))

        weight = get_image(weight_odd_file)
        Util.add_img(weight, get_image(weight_eve_file))

        volall = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
        os.system("rm -f " + fftvol_odd_file + " " + weight_odd_file);
        os.system("rm -f " + fftvol_eve_file + " " + weight_eve_file);
        return volall, fscdat

    if nproc == 2:
        if myid == main_node_odd:
            fftvol = get_image(fftvol_odd_file)
            weight = get_image(weight_odd_file)
            volodd = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
            voleve = recv_EMData(main_node_eve, tag_voleve, mpi_comm)
            fscdat = fsc_mask(volodd, voleve, mask3D, rstep, fsc_curve)
            del volodd, voleve
        else:
            assert myid == main_node_eve
            fftvol = get_image(fftvol_eve_file)
            weight = get_image(weight_eve_file)
            voleve = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
            send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)

        if myid == main_node_odd:
            fftvol = get_image(fftvol_odd_file)
            fftvol_tmp = recv_EMData(main_node_eve, tag_fftvol_eve, mpi_comm)
            Util.add_img(fftvol, fftvol_tmp)
            fftvol_tmp = None

            weight = get_image(weight_odd_file)
            weight_tmp = recv_EMData(main_node_eve, tag_weight_eve, mpi_comm)
            Util.add_img(weight, weight_tmp)
            weight_tmp = None
            volall = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
            os.system("rm -f " + fftvol_odd_file + " " + weight_odd_file);
            return volall, fscdat
        else:
            assert myid == main_node_eve
            fftvol = get_image(fftvol_eve_file)
            send_EMData(fftvol, main_node_odd, tag_fftvol_eve, mpi_comm)

            weight = get_image(weight_eve_file)
            send_EMData(weight, main_node_odd, tag_weight_eve, mpi_comm)
            os.system("rm -f " + fftvol_eve_file + " " + weight_eve_file);
            return model_blank(nx, nx, nx), None
    # cases from all other number of processors situations
    if myid == main_node_odd:
        fftvol = get_image(fftvol_odd_file)
        send_EMData(fftvol, main_node_eve, tag_fftvol_odd, mpi_comm)

        if not (finfo is None):
            finfo.write("fftvol odd sent\n")
            finfo.flush()

        weight = get_image(weight_odd_file)
        send_EMData(weight, main_node_all, tag_weight_odd, mpi_comm)

        if not (finfo is None):
            finfo.write("weight odd sent\n")
            finfo.flush()

        volodd = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
        del fftvol, weight
        voleve = recv_EMData(main_node_eve, tag_voleve, mpi_comm)
        fscdat = fsc_mask(volodd, voleve, mask3D, rstep, fsc_curve)
        del volodd, voleve
        volall = recv_EMData(main_node_all, tag_volall, mpi_comm)
        os.system("rm -f " + fftvol_odd_file + " " + weight_odd_file);
        return volall, fscdat

    if myid == main_node_eve:
        ftmp = recv_EMData(main_node_odd, tag_fftvol_odd, mpi_comm)
        fftvol = get_image(fftvol_eve_file)
        Util.add_img(ftmp, fftvol)
        send_EMData(ftmp, main_node_all, tag_fftvol_eve, mpi_comm)
        del ftmp

        weight = get_image(weight_eve_file)
        send_EMData(weight, main_node_all, tag_weight_eve, mpi_comm)

        voleve = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
        send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)
        os.system("rm -f " + fftvol_eve_file + " " + weight_eve_file);

        return model_blank(nx, nx, nx), None

    if myid == main_node_all:
        fftvol = recv_EMData(main_node_eve, tag_fftvol_eve, mpi_comm)
        if not (finfo is None):
            finfo.write("fftvol odd received\n")
            finfo.flush()

        weight = recv_EMData(main_node_odd, tag_weight_odd, mpi_comm)
        weight_tmp = recv_EMData(main_node_eve, tag_weight_eve, mpi_comm)
        Util.add_img(weight, weight_tmp)
        weight_tmp = None

        volall = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
        send_EMData(volall, main_node_odd, tag_volall, mpi_comm)

        return model_blank(nx, nx, nx), None

def recons3d_nn_SSNR_MPI(myid, prjlist, mask2D, ring_width=1, npad =1, sign=1, symmetry="c1", CTF = False, random_angles = 0, mpi_comm = None):
	from sphire.libpy.sp_utilities import reduce_EMData_to_root
	from EMAN2 import EMData, Reconstructors, Transform, Util
	from mpi import MPI_COMM_WORLD

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	if( len(prjlist) == 0 ):    ERROR("empty input list","recons3d_nn_SSNR_MPI",1)
	imgsize = prjlist[0].get_xsize()
	if prjlist[0].get_ysize() != imgsize:  ERROR("input data has to be square","recons3d_nn_SSNR_MPI",1)
	fftvol   = EMData()
	weight   = EMData()
	weight2  = EMData()
	SSNR     = EMData()
	vol_ssnr = EMData()
	params = {"size":imgsize, "npad":npad, "symmetry":symmetry, "SSNR":SSNR, "fftvol":fftvol, "weight":weight, "weight2":weight2, "vol_ssnr":vol_ssnr, "w":ring_width }
	if CTF:
		weight3  = EMData()
		params["sign"] = sign
		params["weight3"] = weight3
		r = Reconstructors.get("nnSSNR_ctf", params)
	else:
		r = Reconstructors.get("nnSSNR", params)
	r.setup()

	if prjlist[0].get_xsize() != imgsize or prjlist[0].get_ysize() != imgsize: ERROR("inconsistent image size","recons3d_nn_SSNR_MPI",1)
	for prj in prjlist:
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = prj.get_attr_default('active', 1)
		# if active == 1:
		if random_angles  == 2:
			from  random import  random
			phi	 = 360.0*random()
			theta    = 180.0*random()
			psi	 = 360.0*random()
			xform_proj = Transform( {"type":"spider", "phi":phi, "theta":theta, "psi":psi} )
		elif random_angles  == 3:
			from  random import  random
			phi    = 360.0*random()
			theta  = 180.0*random()
			psi    = 360.0*random()
			tx     = 6.0*(random() - 0.5)
			ty     = 6.0*(random() - 0.5)
			xform_proj = Transform( {"type":"spider", "phi":phi, "theta":theta, "psi":psi, "tx":tx, "ty":ty} )
		elif random_angles  == 1:
			from  random import  random
			old_xform_proj = prj.get_attr( "xform.projection" )
			dict = old_xform_proj.get_rotation( "spider" )
			dict["psi"] = 360.0*random()
			xform_proj = Transform( dict )
		else:
			xform_proj = prj.get_attr( "xform.projection" )
		if mask2D:
			stats = Util.infomask(prj, mask2D, True)
			prj -= stats[0]
			prj *= mask2D
		r.insert_slice(prj, xform_proj )
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1 END

	#from utilities import info
	reduce_EMData_to_root(weight,  myid, 0, comm=mpi_comm)
	reduce_EMData_to_root(fftvol,  myid, 0, comm=mpi_comm)
	reduce_EMData_to_root(weight2, myid, 0, comm=mpi_comm)
	if CTF:
		reduce_EMData_to_root(weight3, myid, 0, comm=mpi_comm)
	if myid == 0 :
		dummy = r.finish(True)
		outlist = [[] for i in range(6)]
		nn = SSNR.get_xsize()
		for i in range(1,nn): outlist[0].append((float(i)-0.5)/(float(nn-1)*2))
		for i in range(1,nn):
			if SSNR(i,1,0) > 0.0:
				outlist[1].append(max(0.0,(SSNR(i,0,0)/SSNR(i,1,0)-1.)))     # SSNR
			else:
				outlist[1].append(0.0)
		for i in range(1,nn):
			if SSNR(i,2,0) > 0.0:
				outlist[2].append(SSNR(i,1,0)/SSNR(i,2,0))	          # variance
			else:
				outlist[2].append(0.0)
		for i in range(1,nn): outlist[3].append(SSNR(i,2,0))				  # number of points in the shell
		for i in range(1,nn): outlist[4].append(SSNR(i,3,0))				  # number of added Fourier points
		for i in range(1,nn): outlist[5].append(SSNR(i,0,0))				  # square of signal
		return [outlist, vol_ssnr]


"""Multiline Comment11"""

from builtins import range
