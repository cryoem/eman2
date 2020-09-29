#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
from past.utils import old_div
"""
Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
Author: Toshio Moriya 2019 (toshio.moriya@kek.jp)
Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)

Copyright (c) 2019 Max Planck Institute of Molecular Physiology
Copyright (c) 2000-2006 The University of Texas - Houston Medical School

This software is issued under a joint BSD/GNU license. You may use the source
code in this file under either license. However, note that the complete EMAN2
and SPARX software packages have some GPL dependencies, so you are responsible
for compliance with the licenses of these packages if you opt to use BSD 
licensing. The warranty disclaimer below holds in either instance.

This complete copyright notice must be included in any revised version of the
source code. Additional authorship citations may be added, but existing author
citations must be preserved.

This program is free software; you can redistribute it and/or modifyit under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later 
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 59 Temple
Place, Suite 330, Boston, MA  02111-1307 USA
"""

import EMAN2
import EMAN2_cppwrap
from EMAN2db import db_open_dict
import ctypes
import mpi
import numpy
import os
import optparse
from ..libpy import sp_statistics
from ..libpy import sp_alignment
from ..libpy import sp_applications
from ..libpy import sp_filter
from ..libpy import sp_fundamentals
from ..libpy import sp_global_def
from ..libpy import sp_logger
from ..libpy import sp_pixel_error
from ..libpy import sp_utilities
from ..libpy import sp_isac
import string
import subprocess
import sys
import time
import random

# ====================================================================[ import ]

# compatibility
from future import standard_library

standard_library.install_aliases()
from builtins import range

# system

# python base packages

# EMAN2 & sparx base


# EMAN2 & sparx general utility


# EMAN2 & sparx specific packages


# =================================================================[ mpi setup ]

mpi.mpi_init(0, [])

# ----------------------------------------------------------[ Blockdata ]

"""
Blockdata holds all the administrative information about running ISAC. This
includes:

	Added in header:
	- nproc: number of processes available to MPI (size of MPI_COMM_WORLD)
	- myid: mpi id of this process
	- main_node: mpi id of the mpi main process (traditionally node_0)
	- shared_comm: communicator to other mpi nodes that share the same memory
	- myid_on_node: mpi id of this process within shared_comms ie, the same node
	- no_of_processes_per_group: number of processes in this process shared_comms group

	Added in main():
	- stack: path to particle images to be run through ISAC
	- masterdir: path to ISAC master directory
	- stack_ali2d: path to .bdb file holding the alignment parameters (this file exists only once)
	- total_nima: total number of images in the stack (see above)
	- 2dalignment: path to the results of the pre-alignment

	Added in create_zero_group():
	- subgroup_comm: mpi communicator for all processes with local/node mpi id zero
	- subgroup_size: size of the zero group
	- subgroup_myid: local process id within the zero group

	NOTE: This list is probably not complete yet !
"""

global Blockdata
Blockdata = {}

Blockdata["nproc"] = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
Blockdata["myid"] = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
Blockdata["main_node"] = 0
Blockdata["shared_comm"] = mpi.mpi_comm_split_type(
    mpi.MPI_COMM_WORLD, mpi.MPI_COMM_TYPE_SHARED, 0, mpi.MPI_INFO_NULL
)
Blockdata["myid_on_node"] = mpi.mpi_comm_rank(Blockdata["shared_comm"])
Blockdata["no_of_processes_per_group"] = mpi.mpi_comm_size(Blockdata["shared_comm"])

masters_from_groups_vs_everything_else_comm = mpi.mpi_comm_split(
    mpi.MPI_COMM_WORLD,
    Blockdata["main_node"] == Blockdata["myid_on_node"],
    Blockdata["myid_on_node"],
)
Blockdata["color"], Blockdata[
    "no_of_groups"
], balanced_processor_load_on_nodes = sp_utilities.get_colors_and_subsets(
    Blockdata["main_node"],
    mpi.MPI_COMM_WORLD,
    Blockdata["myid"],
    Blockdata["shared_comm"],
    Blockdata["myid_on_node"],
    masters_from_groups_vs_everything_else_comm,
)
sp_global_def.BATCH = True
sp_global_def.MPI = True

NAME_OF_JSON_STATE_FILE = "my_state.json"
NAME_OF_ORIGINAL_IMAGE_INDEX = "originalid"
NAME_OF_RUN_DIR = "run"
NAME_OF_MAIN_DIR = "generation_"
DIR_DELIM = os.sep


# create an mpi subgroup containing all nodes whith local/node mpi id zero (?)


def create_zero_group():
    if Blockdata["myid_on_node"] == 0:
        submyids = [Blockdata["myid"]]
    else:
        submyids = []

    submyids = sp_utilities.wrap_mpi_gatherv(
        submyids, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )
    submyids = sp_utilities.wrap_mpi_bcast(
        submyids, Blockdata["main_node"], mpi.MPI_COMM_WORLD
    )

    world_group = mpi.mpi_comm_group(mpi.MPI_COMM_WORLD)
    subgroup = mpi.mpi_group_incl(world_group, len(submyids), submyids)

    Blockdata["subgroup_comm"] = mpi.mpi_comm_create(mpi.MPI_COMM_WORLD, subgroup)
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    Blockdata["subgroup_size"] = -1
    Blockdata["subgroup_myid"] = -1

    if mpi.MPI_COMM_NULL != Blockdata["subgroup_comm"]:
        Blockdata["subgroup_size"] = mpi.mpi_comm_size(Blockdata["subgroup_comm"])
        Blockdata["subgroup_myid"] = mpi.mpi_comm_rank(Blockdata["subgroup_comm"])

    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    return


def checkitem(item, mpi_comm=-1):
    global Blockdata
    if mpi_comm == -1:
        mpi_comm = mpi.MPI_COMM_WORLD
    if Blockdata["myid"] == Blockdata["main_node"]:
        if os.path.exists(item):
            isthere = True
        else:
            isthere = False
    else:
        isthere = False
    isthere = sp_utilities.bcast_number_to_all(
        isthere, source_node=Blockdata["main_node"], mpi_comm=mpi_comm
    )
    mpi.mpi_barrier(mpi_comm)
    return isthere


# ===================================================================[ utility ]


def normalize_particle_images(
    aligned_images,
    shrink_ratio,
    target_radius,
    target_dim,
    align_params,
    filament_width=-1,
    ignore_helical_mask=False,
):
    """
	Function to normalize the images given in <aligned_images>. Note that the
	normalization also includes the shrinking/re-scaling of the particle images.
	The normalization itself is being done by subtracting the mean of the data
	inside a particle mask (signal) and dividing by the variance outsice of the
	mask (noise).

	NOTE: Images will first be shrunken/re-scaled in order to make the particle
	radius match between what ISAC uses for processing and the actual data.
	Afterwards the re-scaled image will be padded/cropped in order to match the
	internal image size used by ISAC.

	Args:
		aligned_images (EMData[]): List of EMData objects holding image data.

		shrink_ratio (float): Ratio by which particles are to be re-scaled.

		target_radius (int): ISAC target radius.

		target_dim (int): ISAC target image size for processing (usually 76).

		align_params (list of lists): Contains (pre-)alignment parameters of
			which we apply the shifts.

		filament_width (int): Filament width when processing helical data. When
			a non-default value is provided this function assumes data to be
			filament images in which case a rectangular mask of the given width
			is applied to all particle images. [Default: -1]

		ignore_helical_mask (bool): Only relevant if filament_width is used. If
			set to False the data will be multiplied with the mask to remove all
			data/noise outside of the mask. If set to True the mask will still
			be used for the normalization but afterwards NOT multiplied with the
			particle images. [Default: False]
	"""

    # particle image dimension after scaling/shrinking
    new_dim = int(aligned_images[0].get_xsize() * shrink_ratio + 0.5)

    # create re-usable mask for non-helical particle images
    if filament_width == -1:
        if new_dim >= target_dim:
            mask = sp_utilities.model_circle(target_radius, target_dim, target_dim)
        else:
            mask = sp_utilities.model_circle(old_div(new_dim, 2) - 2, new_dim, new_dim)

    # pre-process all images
    for im in range(len(aligned_images)):

        # apply any available alignment parameters
        aligned_images[im] = sp_fundamentals.rot_shift2D(
            aligned_images[im], 0, align_params[im][1], align_params[im][2], 0
        )

        # resample if necessary
        if shrink_ratio != 1.0:
            aligned_images[im] = sp_fundamentals.resample(
                aligned_images[im], shrink_ratio
            )

        # crop images if necessary
        if new_dim > target_dim:
            aligned_images[im] = EMAN2.Util.window(
                aligned_images[im], target_dim, target_dim, 1
            )

        current_dim = aligned_images[im].get_xsize()
        # create custom masks for filament particle images
        if filament_width != -1:
            mask = sp_utilities.model_rotated_rectangle2D(
                radius_long=old_div(
                    int(numpy.sqrt(2 * current_dim ** 2)), 2
                ),  # long  edge of the rectangular mask
                radius_short=old_div(
                    int(filament_width * shrink_ratio + 0.5), 2
                ),  # short edge of the rectangular mask
                nx=current_dim,
                ny=current_dim,
                angle=aligned_images[im].get_attr("segment_angle"),
            )

        # normalize using mean of the data and variance of the noise
        p = EMAN2.Util.infomask(aligned_images[im], mask, False)
        aligned_images[im] -= p[0]
        p = EMAN2.Util.infomask(aligned_images[im], mask, True)
        aligned_images[im] = old_div(aligned_images[im], p[1])

        # optional: burn helical mask into particle images
        if filament_width != -1 and not ignore_helical_mask:
            aligned_images[im] *= mask

        # pad images in case they have been shrunken below the target_dim
        if new_dim < target_dim:
            aligned_images[im] = sp_utilities.pad(
                aligned_images[im], target_dim, target_dim, 1, 0.0
            )


# ======================================================================[ ISAC ]


def iter_isac_pap(
    alldata,
    ir,
    ou,
    rs,
    xr,
    yr,
    ts,
    maxit,
    CTF,
    snr,
    dst,
    FL,
    FH,
    FF,
    init_iter,
    main_iter,
    iter_reali,
    match_first,
    max_round,
    match_second,
    stab_ali,
    thld_err,
    indep_run,
    thld_grp,
    img_per_grp,
    generation,
    candidatesexist=False,
    random_seed=None,
    new=False,
):
    """
	Core function to set up the next iteration of ISAC.

	Args:
		alldata (UNKNOWN TYPE): All image data [nima=len(alldata): no. of images]

		ir (int): Inner ring value (in pixels) of the resampling to polar
			coordinates.
			[Default: 1]

		ou (int): Target particle radius used when ISAC processes the data.
			Images will be scaled to conform to this value.
			[Default: 29]

		rs (int): Ring step in pixels used during the resampling of images to
			polar coordinates.
			[Default: 1]

		xr (int): x-range of translational search during alignment.
			[Default: 1]

		yr (int): y-range of translational search during alignment.
			[Default: 1]

		ts (float): Search step size (in pixels) of translational search. (Not
			entirely clear; used in refinement.)
			[Default 1.0]

		maxit (int): Number of iterations for reference-free alignment.
			[Default: 30]

		CTF (bool): If set the data will be phase-flipped using CTF information
			included in image headers.
			[Default: False][UNSUPPORTED]

		snr (float): Signal to noise ratio used if CTF parameter is True.
			[Default: 1.0][UNSUPPORTED]

		dst (float): Discrete angle used during group alignment.
			[Default: 90.0]

		FL (float): Frequency of the lowest stop band used in the tangent filter.
			[Default: 0.2]

		FH (float): Frequency of the highest stop band used in the tangent filter.
			[Default 0.45]

		FF (float): Fall-off value for the tangent filter.
			[Default 0.2]

		init_iter (int): Maximum number of Generation iterations performed for
			a given subset. (Unclear what "subset" refers to.)
			[Default: 7]

		main_iter (int): Number of the current main iteration.

		iter_reali (int): SAC stability check interval. Every iter_reali-th
			iteration of SAC stability checking is performed.
			[Default: 1]

		match_first (int): Number of iterations to run 2-way matching during
			the first phase. (Unclear what "first hpase" refers to.)
			[Default: 1][UNUSED]

		max_round (int): Maximum number of rounds generating candidate class
			averages during the first phase. (Unclear what first phase means.)
			[Default: 20][UNUSED]

		match_second (int): Number of times to run 2-way (or 3-way) matching in
			the second phase. (Unclear what "second phase" refers to.)
			[Default: 5][UNUSED]

		stab_ali (int): Number of alignment iterations when checking stability.
			[Default: 5]

		thld_err (float): Threshold of pixel error when checking stability.
			Equals root mean square of distances between corresponding pixels
			from set of found transformations and theirs average transformation;
			depends linearly on square of radius (parameter target_radius).
			Units are in pixels.
			[Default: 0.7]

		indep_run (int): [Default: 0][UNUSED]

		thld_grp (int): Minimum size of reproducible classes.
			[Default 10/30][UNUSED]

		img_per_grp (int): Number of images per class (maximum group size, also
			defines number of classes: K=(total number of images)/img_per_grp.
			[Default: 200]

		generation (int): Number of iterations in the current generation.

		candidatesexist (bool): Candidate class averages already exist and can
			be used.
			[Default: False][UNUSED]

		random_seed (int): Set random seed manually for testing purposes.
			[Default: None]

		new (bool): Flag to use "new code"; set to False and not used.

	Returns:
		refi (UNKNOWN TYPE): Refinement as returned by isac_MPI_pap()

		all_ali_params (list): List containing 2D alignment parameters for all
			images; entries formatted as [angle, sx, sy, mirror].
	"""

    # ------------------------------------------------------[ mpi ]

    global Blockdata

    number_of_proc = Blockdata["nproc"]
    myid = Blockdata["myid"]
    main_node = Blockdata["main_node"]

    random.seed(myid)
    rand1 = random.randint(1, 1000111222)
    random.seed(random_seed)
    rand2 = random.randint(1, 1000111222)
    random.seed(rand1 + rand2)

    if generation == 0:
        sp_global_def.ERROR(
            "Generation should begin from 1, please reset it and restart the program",
            myid=myid,
        )
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    ali_params_dir = "ali_params_generation_%d" % generation
    if os.path.exists(ali_params_dir):
        sp_global_def.ERROR(
            "Output directory %s for alignment parameters exists, please either change its name or delete it and restart the program"
            % ali_params_dir,
            myid=myid,
        )
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    if new:
        alimethod = "SHC"
    else:
        alimethod = ""

    color = 0  # Blockdata["Ycolor"]
    key = Blockdata["myid"]  # Blockdata["Ykey"]
    group_comm = mpi.MPI_COMM_WORLD  # Blockdata["Ygroup_comm"]
    group_main_node = 0

    nx = alldata[0].get_xsize()
    ndata = len(alldata)
    data = [None] * ndata
    tdummy = EMAN2.Transform({"type": "2D"})

    for im in range(ndata):
        # This is the absolute ID, the only time we use it is
        # when setting the members of 4-way output. All other times, the id in 'members' is
        # the relative ID.
        alldata[im].set_attr_dict({"xform.align2d": tdummy, "ID": im})
        data[im] = alldata[im]

    avg_num = 0
    Iter = 1
    K = old_div(ndata, img_per_grp)

    if myid == main_node:
        sp_global_def.sxprint(
            "     We will process:  %d current images divided equally between %d groups"
            % (ndata, K)
        )
        sp_global_def.sxprint(
            "*************************************************************************************"
        )

    # ------------------------------------------------------[ generate random averages for each group ]

    if key == group_main_node:
        refi = sp_isac.generate_random_averages(data, K, 9023)
    else:
        refi = [sp_utilities.model_blank(nx, nx) for i in range(K)]

    for i in range(K):
        sp_utilities.bcast_EMData_to_all(refi[i], key, group_main_node, group_comm)

    # create d[K*ndata] matrix
    orgsize = K * ndata

    if Blockdata["myid_on_node"] == 0:
        size = orgsize
    else:
        size = 0

    disp_unit = numpy.dtype("f4").itemsize

    win_sm, base_ptr = mpi.mpi_win_allocate_shared(
        size * disp_unit, disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"]
    )
    size = orgsize
    if Blockdata["myid_on_node"] != 0:
        base_ptr, = mpi.mpi_win_shared_query(win_sm, mpi.MPI_PROC_NULL)

    # d = numpy.frombuffer(
    #     numpy.core.multiarray.int_asbuffer(base_ptr, size * disp_unit), dtype="f4"
    # )

    ptr_n = ctypes.cast(base_ptr, ctypes.POINTER(ctypes.c_int * size))
    d = numpy.frombuffer(ptr_n.contents, dtype="f4")

    d = d.reshape(orgsize)
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    # ------------------------------------------------------[ Generate inital averages ]

    refi = isac_MPI_pap(
        data,
        refi,
        d,
        maskfile=None,
        ir=ir,
        ou=ou,
        rs=rs,
        xrng=xr,
        yrng=yr,
        step=ts,
        maxit=maxit,
        isac_iter=init_iter,
        CTF=CTF,
        snr=snr,
        rand_seed=-1,
        color=color,
        comm=group_comm,
        stability=True,
        stab_ali=stab_ali,
        iter_reali=iter_reali,
        thld_err=thld_err,
        FL=FL,
        FH=FH,
        FF=FF,
        dst=dst,
        method=alimethod,
    )

    mpi.mpi_win_free(win_sm)

    del d

    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    if myid == main_node:
        all_ali_params = [None] * len(data)
        for i, im in enumerate(data):
            alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(im)
            all_ali_params[i] = [alpha, sx, sy, mirror]

        sp_global_def.sxprint(
            "****************************************************************************************************"
        )
        sp_global_def.sxprint(
            "*         Generation finished                 "
            + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
            + "                            *"
        )
        sp_global_def.sxprint(
            "****************************************************************************************************"
        )

        return refi, all_ali_params
    else:
        return [], []


def isac_MPI_pap(
    stack,
    refim,
    d,
    maskfile=None,
    ir=1,
    ou=-1,
    rs=1,
    xrng=0,
    yrng=0,
    step=1,
    maxit=30,
    isac_iter=10,
    CTF=False,
    snr=1.0,
    rand_seed=-1,
    color=0,
    comm=-1,
    stability=False,
    stab_ali=5,
    iter_reali=1,
    thld_err=1.732,
    FL=0.1,
    FH=0.3,
    FF=0.2,
    dst=90.0,
    method="",
):
    """
	ISAC core function.

	Args:
		stack (UNKNOWN TYPE): list of images (filename also accepted)

		refim (list OR filename): Class averages. (Providing a filename exits.)

		d (numy.ndarray): Array holding pairwise distances between images.

		maskfile (image OR filename): Image containing mask (filename also
			accepted).

		ir (int): Inner ring value (in pixels) of the resampling to polar
			coordinates.
			[Default: 1]

		ou (int): Target particle radius used when ISAC processes the data.
			Images will be scaled to conform to this value.
			[Default: 29]

		rs (int): Ring step in pixels used during the resampling of images to
			polar coordinates.
			[Default: 1]

		xrng (int): x-range of translational search during alignment.
			[Default: 0]

		yrng (int): y-range of translational search during alignment.
			[Default: 0]

		step (float): Search step size (in pixels) of translational search.
			(Not entirely clear; used in refinement.)
			[Default 1]

		maxit (int): Number of iterations for reference-free alignment.
			[Default: 30]

		isac_iter (UNKNOWN TYPE): Maximum number of Generation iterations
			performed for a given subset. (Unclear what "subset" refers to.)
			[Default: 10]

		CTF (bool): If set the data will be phase-flipped using CTF information
			included in image headers.
			[Default: False][UNSUPPORTED]

		snr (float): Signal to noise ratio used if CTF parameter is True.
			[Default: 1.0][UNSUPPORTED]

		rand_seed (int): Set random seed manually for testing purposes.
			[Default: -1]

		color (mpi color): set to 0; unclear if this is still relevant.
			[Defailt: 0]

		comm (mpi communicator): set to MPI_COMM_WORLD (globally available; redundant parameter)
			[Default: -1]

		stability (bool): If True, ISAC performs stability testing.
			[Default: True]

		stab_ali (bool): Used only when stability testing is performed.
			[Default: 5]

		inter_reali (UNKNOWN TYPE): Used only when stability=True. For each
			iteration i: if (i % iter_reali == 0) stability check is performed.
			[Default: 1]

		thld_err (float): Threshold of pixel error when checking stability.
			Equals root mean square of distances between corresponding pixels
			from set of found transformations and theirs average transformation;
			depends linearly on square of radius (parameter target_radius).
			Units are in pixels.
			[Default: 1.732]

		FL (float): Frequency of the lowest stop band used in the tangent filter.
			[Default: 0.1]

		FH (float): Frequency of the highest stop band used in the tangent filter.
			[Default 0.3]

		FF (float): Fall-off value for the tangent filter.
			[Default 0.2]

		dst (float): Discrete angle used during group alignment.
			[Default: 90.0]

		method (string): Stock method (SHC) for alignment.
			[Default: ""]

	Returns:
		alldata (list): Class averages (format unclear).
	"""

    # ------------------------------------------------------[ initialize mpi ]

    global Blockdata

    if comm == -1:
        comm = mpi.MPI_COMM_WORLD

    number_of_proc = mpi.mpi_comm_size(comm)
    myid = mpi.mpi_comm_rank(comm)
    my_abs_id = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
    main_node = 0

    # ------------------------------------------------------[ initialize ]

    first_ring = int(ir)
    last_ring = int(ou)
    rstep = int(rs)
    max_iter = int(isac_iter)

    if type(stack) == type(""):
        # read all data
        sp_global_def.sxprint("  SHOULD NOT BE HERE")
        sys.exit()
        alldata = EMAN2_cppwrap.EMData.read_images(stack)
    else:
        alldata = stack
    nx = alldata[0].get_xsize()
    ny = alldata[0].get_ysize()

    nima = len(alldata)  # no of images

    # set all parameters to be zero on input
    for im in range(nima):
        sp_utilities.set_params2D(alldata[im], [0.0, 0.0, 0.0, 0, 1.0])

    image_start, image_end = sp_applications.MPI_start_end(nima, number_of_proc, myid)

    if maskfile:
        if isinstance(maskfile, (bytes, str)):
            mask = sp_utilities.get_image(maskfile)
        else:
            mask = maskfile

    else:
        mask = sp_utilities.model_circle(last_ring, nx, nx)

    if type(refim) == type(""):
        refi = EMAN2_cppwrap.EMData.read_images(refim)
    else:
        # It's safer to make a hard copy here. Although I am not sure, I believe a shallow copy
        # has messed up the program.
        #   This is really strange.  It takes much memory without any need.  PAP 01/17/2015
        #      However, later I made changes so refi is deleted from time to time.  All to be checked.
        # refi = refim
        refi = [None for i in range(len(refim))]
        for i in range(len(refim)):
            refi[i] = refim[i].copy()

    numref = len(refi)

    # IMAGES ARE SQUARES! center is in SPIDER convention
    cnx = old_div(nx, 2) + 1
    cny = cnx

    mode = "F"

    # ------------------------------------------------------[ precalculate rings]

    # calculates the necessary information for the 2D polar interpolation (finds the number of rings)
    numr = sp_alignment.Numrinit(first_ring, last_ring, rstep, mode)
    # Calculate ring weights for rotational alignment
    wr = sp_alignment.ringwe(numr, mode)

    if rand_seed > -1:
        random.seed(rand_seed)
    else:
        random.seed(random.randint(1, 2000111222))

    """
    Comment from Adnan
    random.jumpahead function is not available in python 3 . That is why i comment it
    """
    # if myid != main_node:
    #     random.jumpahead(17 * myid + 12345)

    previous_agreement = 0.0
    previous_members = [None] * numref
    for j in range(numref):
        previous_members[j] = set()

    fl = FL
    Iter = -1
    main_iter = 0
    terminate = 0
    while (main_iter < max_iter) and (terminate == 0):

        Iter += 1
        if my_abs_id == main_node:
            sp_global_def.sxprint(
                "Iteration within isac_MPI Iter =>",
                Iter,
                "	main_iter = ",
                main_iter,
                "	len data = ",
                image_end - image_start,
                "   ",
                time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()),
            )
        mashi = cnx - ou - 2

        for j in range(numref):
            refi[j].process_inplace(
                "normalize.mask", {"mask": mask, "no_sigma": 1}
            )  # normalize reference images to N(0,1)
            cimage = EMAN2_cppwrap.Util.Polar2Dm(
                refi[j], cnx, cny, numr, mode
            )  # converting reference images to polar cordinates
            EMAN2_cppwrap.Util.Frngs(
                cimage, numr
            )  # Applying FFT on the reference images
            EMAN2_cppwrap.Util.Applyws(
                cimage, numr, wr
            )  # apply weights to FTs of rings
            refi[j] = cimage.copy()

        # if CTF: ctf2 = [[[0.0]*lctf for k in xrange(2)] for j in xrange(numref)]
        peak_list = [
            numpy.zeros(4 * (image_end - image_start), dtype=numpy.float32)
            for i in range(numref)
        ]
        #  nima is the total number of images, not the one on this node, the latter is (image_end-image_start)
        #    d matrix required by EQ-Kmeans can be huge!!  PAP 01/17/2015
        # d = zeros(numref*nima, dtype=float32)
        if Blockdata["myid_on_node"] == 0:
            d.fill(0.0)

        # --------------------------------------------------[ compute the 2D alignment  ]

        # loop through all images and compute alignment parameters for all references
        for im in range(image_start, image_end):
            alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(
                alldata[im]
            )  # retrieves 2D alignment parameters
            ##  TEST WHETHER PARAMETERS ARE WITHIN RANGE
            alphai, sxi, syi, scalei = sp_utilities.inverse_transform2(
                alpha, sx, sy
            )  # returns the inverse  geometric transformof the 2d rot and trans matrix
            # If shifts are outside of the permissible range, reset them
            if abs(sxi) > mashi or abs(syi) > mashi:
                sxi = 0.0
                syi = 0.0
                sp_utilities.set_params2D(
                    alldata[im], [0.0, 0.0, 0.0, 0, 1.0]
                )  # set 2D alignment parameters (img, alpha, tx, ty, mirror, scale)
            # normalize
            alldata[im].process_inplace(
                "normalize.mask", {"mask": mask, "no_sigma": 0}
            )  # subtract average under the mask
            txrng = sp_alignment.search_range(
                nx, ou, sxi, xrng, "ISAC2"
            )  # Find permissible ranges for translational searches by resampling into polar coordinates
            # nx-> image size; ou -> particle radius, sxi-> particle shift, xrng-> max search rng
            txrng = [txrng[1], txrng[0]]  # for positive and negative shifts
            tyrng = sp_alignment.search_range(ny, ou, syi, yrng, "ISAC2")
            tyrng = [tyrng[1], tyrng[0]]

            # ----------------------------------------------[ align current image to references]

            # compute alignment parameters for single image with index <im> and all references <refi>
            # NOTE: multiref_polar_ali_2d_peaklist includes conversion of image data into polar coordinates AND an FFT
            temp = EMAN2_cppwrap.Util.multiref_polar_ali_2d_peaklist(
                alldata[im], refi, txrng, tyrng, step, mode, numr, cnx + sxi, cny + syi
            )
            for iref in range(numref):
                [alphan, sxn, syn, mn] = sp_utilities.combine_params2(
                    0.0,
                    -sxi,
                    -syi,
                    0,
                    temp[iref * 5 + 1],
                    temp[iref * 5 + 2],
                    temp[iref * 5 + 3],
                    int(temp[iref * 5 + 4]),
                )  # Combine 2D alignment parameters including mirror
                peak_list[iref][(im - image_start) * 4 + 0] = alphan
                peak_list[iref][(im - image_start) * 4 + 1] = sxn
                peak_list[iref][(im - image_start) * 4 + 2] = syn
                peak_list[iref][(im - image_start) * 4 + 3] = mn
                d[iref * nima + im] = temp[
                    iref * 5
                ]  # grab a chunk of the peak list for parallel processing in the next step
            del temp

        del refi

        # --------------------------------------------------[ reduce the results of the correlation computation ]

        if Blockdata["subgroup_myid"] > -1:
            # First check number of nodes, if only one, no reduction necessary.
            if Blockdata["no_of_groups"] > 1:
                # do reduction using numref chunks nima long (it does not matter what is the actual ordering in d)
                at = time.time()
                for j in range(numref):
                    dbuf = numpy.zeros(nima, dtype=numpy.float32)
                    numpy.copyto(dbuf, d[j * nima : (j + 1) * nima])

                    # reduce entries in section <dbuf> of the overall peak_list
                    dbuf = mpi.mpi_reduce(
                        dbuf,
                        nima,
                        mpi.MPI_FLOAT,
                        mpi.MPI_SUM,
                        main_node,
                        Blockdata["subgroup_comm"],
                    )  # RETURNS numpy array
                    if Blockdata["subgroup_myid"] == 0:
                        numpy.copyto(d[j * nima : (j + 1) * nima], dbuf)

                del dbuf
		
        mpi.mpi_barrier(comm)

        # --------------------------------------------------[ find alignment match ]

        if myid == main_node:
            # find maximum in the peak list to determine highest subject/reference match
            id_list_long = EMAN2_cppwrap.Util.assign_groups(
                str(d.__array_interface__["data"][0]), numref, nima
            )  # string with memory address is passed as parameters
            id_list = [[] for i in range(numref)]
            maxasi = old_div(nima, numref)  # max. assignmenx

            for i in range(maxasi * numref):
                id_list[old_div(i, maxasi)].append(id_list_long[i])
            for i in range(nima % maxasi):
                id_list[id_list_long[-1]].append(id_list_long[maxasi * numref + i])
            for iref in range(numref):
                id_list[iref].sort()

            del id_list_long

            belongsto = [0] * nima
            for iref in range(numref):
                for im in id_list[iref]:
                    belongsto[
                        im
                    ] = iref  # image <im> has highest match with reference <iref>

            del id_list

        else:
            belongsto = [0] * nima

        # broadcast assignment result to all mpi nodes
        mpi.mpi_barrier(comm)
        belongsto = mpi.mpi_bcast(belongsto, nima, mpi.MPI_INT, main_node, comm)
        belongsto = list(map(int, belongsto))

        # --------------------------------------------------[ update the averages / create new references ]

        members = [0] * numref  # stores number of images assigned to each reference
        sx_sum = [0.0] * numref
        sy_sum = [0.0] * numref
        refi = [sp_utilities.model_blank(nx, ny) for j in range(numref)]

        for im in range(image_start, image_end):

            # get alignment parameters for the chosen matching reference
            matchref = belongsto[im]
            alphan = float(peak_list[matchref][(im - image_start) * 4 + 0])
            sxn = float(peak_list[matchref][(im - image_start) * 4 + 1])
            syn = float(peak_list[matchref][(im - image_start) * 4 + 2])
            mn = int(peak_list[matchref][(im - image_start) * 4 + 3])

            if mn == 0:
                sx_sum[matchref] += sxn
            else:
                sx_sum[matchref] -= sxn
            sy_sum[matchref] += syn

            # apply current parameters and add to the average
            EMAN2_cppwrap.Util.add_img(
                refi[matchref],
                sp_fundamentals.rot_shift2D(alldata[im], alphan, sxn, syn, mn),
            )
            members[matchref] += 1

        # --------------------------------------------------[ broadcast  ]

        # everyone computes their reduction and then broadcasts it back to the main_node
        sx_sum = mpi.mpi_reduce(
            sx_sum, numref, mpi.MPI_FLOAT, mpi.MPI_SUM, main_node, comm
        )
        sy_sum = mpi.mpi_reduce(
            sy_sum, numref, mpi.MPI_FLOAT, mpi.MPI_SUM, main_node, comm
        )
        members = mpi.mpi_reduce(
            members, numref, mpi.MPI_INT, mpi.MPI_SUM, main_node, comm
        )

        if myid != main_node:
            sx_sum = [0.0] * numref
            sy_sum = [0.0] * numref
            members = [0.0] * numref

        sx_sum = mpi.mpi_bcast(sx_sum, numref, mpi.MPI_FLOAT, main_node, comm)
        sy_sum = mpi.mpi_bcast(sy_sum, numref, mpi.MPI_FLOAT, main_node, comm)
        members = mpi.mpi_bcast(members, numref, mpi.MPI_INT, main_node, comm)
        sx_sum = list(map(float, sx_sum))
        sy_sum = list(map(float, sy_sum))
        members = list(map(int, members))

        for j in range(numref):
            sx_sum[j] = old_div(sx_sum[j], float(members[j]))
            sy_sum[j] = old_div(sy_sum[j], float(members[j]))

        for im in range(image_start, image_end):
            matchref = belongsto[im]
            alphan = float(peak_list[matchref][(im - image_start) * 4 + 0])
            sxn = float(peak_list[matchref][(im - image_start) * 4 + 1])
            syn = float(peak_list[matchref][(im - image_start) * 4 + 2])
            mn = int(peak_list[matchref][(im - image_start) * 4 + 3])

            if mn == 0:
                sp_utilities.set_params2D(
                    alldata[im],
                    [alphan, sxn - sx_sum[matchref], syn - sy_sum[matchref], mn, scale],
                )  # set 2D alignment parameters (img, alpha, tx, ty, mirror, scale)
            else:
                sp_utilities.set_params2D(
                    alldata[im],
                    [alphan, sxn + sx_sum[matchref], syn - sy_sum[matchref], mn, scale],
                )

        del peak_list

        # --------------------------------------------------[ ... we work till here (fabian & adnan)]

        for j in range(numref):
            sp_utilities.reduce_EMData_to_root(refi[j], myid, main_node, comm)

            if myid == main_node:
                # Golden rule when to do within group refinement
                EMAN2_cppwrap.Util.mul_scalar(refi[j], old_div(1.0, float(members[j])))
                refi[j] = sp_filter.filt_tanl(refi[j], fl, FF)
                refi[j] = sp_fundamentals.fshift(refi[j], -sx_sum[j], -sy_sum[j])
                sp_utilities.set_params2D(refi[j], [0.0, 0.0, 0.0, 0, 1.0])

        if myid == main_node:
            # this is most likely meant to center them, if so, it works poorly,
            # it has to be checked and probably a better method used PAP 01/17/2015
            dummy = sp_applications.within_group_refinement(
                refi,
                mask,
                True,
                first_ring,
                last_ring,
                rstep,
                [xrng],
                [yrng],
                [step],
                dst,
                maxit,
                FH,
                FF,
            )
            ref_ali_params = []

            for j in range(numref):
                alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(refi[j])
                refi[j] = sp_fundamentals.rot_shift2D(refi[j], alpha, sx, sy, mirror)
                ref_ali_params.extend([alpha, sx, sy, mirror])
        else:
            ref_ali_params = [0.0] * (numref * 4)

        ref_ali_params = mpi.mpi_bcast(
            ref_ali_params, numref * 4, mpi.MPI_FLOAT, main_node, comm
        )
        ref_ali_params = list(map(float, ref_ali_params))

        for j in range(numref):
            sp_utilities.bcast_EMData_to_all(refi[j], myid, main_node, comm)

        # --------------------------------------------------[ Compensate the centering to averages ]
        for im in range(image_start, image_end):
            matchref = belongsto[im]
            alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(alldata[im])
            alphan, sxn, syn, mirrorn = sp_utilities.combine_params2(
                alpha,
                sx,
                sy,
                mirror,
                ref_ali_params[matchref * 4 + 0],
                ref_ali_params[matchref * 4 + 1],
                ref_ali_params[matchref * 4 + 2],
                int(ref_ali_params[matchref * 4 + 3]),
            )
            sp_utilities.set_params2D(alldata[im], [alphan, sxn, syn, int(mirrorn), 1.0])

        mpi.mpi_barrier( mpi.MPI_COMM_WORLD )
	
        fl += 0.05

        if fl >= FH:
            fl = FL
            do_within_group = 1
        else:
            do_within_group = 0

        # Here stability does not need to be checked for each main iteration, it only needs to
        # be done for every 'iter_reali' iterations. If one really wants it to be checked each time
        # simple set iter_reali to 1, which is the default value right now.
        check_stability = stability and (main_iter % iter_reali == 0)

        # --------------------------------------------------[ within group alignment ]

        if do_within_group == 1:

            # ----------------------------------------------[ Broadcast the alignment parameters to all nodes ]

            for i in range(number_of_proc):
                im_start, im_end = sp_applications.MPI_start_end(
                    nima, number_of_proc, i
                )
                if myid == i:
                    ali_params = []
                    for im in range(image_start, image_end):
                        alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(
                            alldata[im]
                        )
                        ali_params.extend([alpha, sx, sy, mirror])
                else:
                    ali_params = [0.0] * ((im_end - im_start) * 4)
                ali_params = mpi.mpi_bcast(
                    ali_params, len(ali_params), mpi.MPI_FLOAT, i, comm
                )
                ali_params = list(map(float, ali_params))
                for im in range(im_start, im_end):
                    alpha = ali_params[(im - im_start) * 4]
                    sx = ali_params[(im - im_start) * 4 + 1]
                    sy = ali_params[(im - im_start) * 4 + 2]
                    mirror = int(ali_params[(im - im_start) * 4 + 3])
                    sp_utilities.set_params2D(alldata[im], [alpha, sx, sy, mirror, 1.0])

            main_iter += 1
	
            mpi.mpi_barrier( mpi.MPI_COMM_WORLD )

            # There are two approaches to scatter calculations among MPI processes during stability checking.
            # The first one is the original one. I added the second method.
            # Here we try to estimate the calculation time for both approaches.
            stab_calc_time_method_1 = stab_ali * (
                old_div((numref - 1), number_of_proc) + 1
            )
            stab_calc_time_method_2 = (
                old_div((numref * stab_ali - 1), number_of_proc) + 1
            )
            ##if my_abs_id == main_node: print "Times estimation: ", stab_calc_time_method_1, stab_calc_time_method_2

            # When there is no stability checking or estimated calculation time of new method is greater than 80% of estimated calculation time of original method
            # then the original method is used. In other case. the second (new) method is used.
            # if (not check_stability) or (stab_calc_time_method_2 > 0.80 * stab_calc_time_method_1):
            #  For the time being only use this method as the other one is not worked out as far as parameter ranges go.
            # if True :
            ##if my_abs_id == main_node: print "Within group refinement and checking within group stability, original approach .......", check_stability, "  ",time.localtime()[0:5]
            # ====================================== standard approach is used, calculations are parallelized by scatter groups (averages) among MPI processes
            gpixer = []
            mpi.mpi_barrier(comm)
            for j in range(myid, numref, number_of_proc):
                assign = []
                for im in range(nima):
                    if j == belongsto[im]:
                        assign.append(im)

                randomize = True  # I think there is no reason not to be True
                class_data = [alldata[im] for im in assign]
                refi[j] = sp_applications.within_group_refinement(
                    class_data,
                    mask,
                    randomize,
                    first_ring,
                    last_ring,
                    rstep,
                    [xrng],
                    [yrng],
                    [step],
                    dst,
                    maxit,
                    FH,
                    FF,
                    method=method,
                )

                if check_stability:
                    ###if my_abs_id == main_node: print "Checking within group stability, original approach .......", check_stability, "  ",time.localtime()[0:5]
                    ali_params = [[] for qq in range(stab_ali)]
                    for ii in range(stab_ali):
                        if ii > 0:  # The first one does not have to be repeated
                            dummy = sp_applications.within_group_refinement(
                                class_data,
                                mask,
                                randomize,
                                first_ring,
                                last_ring,
                                rstep,
                                [xrng],
                                [yrng],
                                [step],
                                dst,
                                maxit,
                                FH,
                                FF,
                                method=method,
                            )
                        for im in range(len(class_data)):
                            alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(
                                class_data[im]
                            )
                            ali_params[ii].extend([alpha, sx, sy, mirror])

                    stable_set, mirror_consistent_rate, err = sp_pixel_error.multi_align_stability(
                        ali_params, 0.0, 10000.0, thld_err, False, last_ring * 2
                    )
                    gpixer.append(err)

                    # print  "Class %4d ...... Size of the group = %4d and of the stable subset = %4d  Mirror consistent rate = %5.3f  Average pixel error prior to class pruning = %10.2f"\
                    # 				%(j, len(class_data), len(stable_set), mirror_consistent_rate, err)

                    # If the size of stable subset is too small (say 1, 2), it will cause many problems, so we manually increase it to 5
                    while len(stable_set) < 5:
                        duplicate = True
                        while duplicate:
                            duplicate = False
                            p = random.randint(0, len(class_data) - 1)
                            for ss in stable_set:
                                if p == ss[1]:
                                    duplicate = True
                        stable_set.append([100.0, p, [0.0, 0.0, 0.0, 0]])
                    stable_data = []
                    stable_members = []
                    for err in stable_set:
                        im = err[1]
                        stable_members.append(assign[im])
                        stable_data.append(class_data[im])
                        sp_utilities.set_params2D(
                            class_data[im],
                            [err[2][0], err[2][1], err[2][2], int(err[2][3]), 1.0],
                        )
                    stable_members.sort()

                    refi[j] = sp_filter.filt_tanl(
                        sp_statistics.ave_series(stable_data), FH, FF
                    )
                    refi[j].set_attr("members", stable_members)
                    refi[j].set_attr("n_objects", len(stable_members))
                    # print  "Class %4d ...... Size of the stable subset = %4d  "%(j, len(stable_members))
                    del stable_members
                # end of stability
                del assign
            mpi.mpi_barrier(comm)

            for im in range(nima):
                done_on_node = belongsto[im] % number_of_proc
                if myid == done_on_node:
                    alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(alldata[im])
                    ali_params = [alpha, sx, sy, mirror]
                else:
                    ali_params = [0.0] * 4
                ali_params = mpi.mpi_bcast(
                    ali_params, 4, mpi.MPI_FLOAT, done_on_node, comm
                )
                ali_params = list(map(float, ali_params))
                sp_utilities.set_params2D(
                    alldata[im],
                    [
                        ali_params[0],
                        ali_params[1],
                        ali_params[2],
                        int(ali_params[3]),
                        1.0,
                    ],
                )

            for j in range(numref):
                sp_utilities.bcast_EMData_to_all(refi[j], myid, j % number_of_proc, comm)

            terminate = 0
            if check_stability:
                # In this case, we need to set the 'members' attr using stable members from the stability test
                for j in range(numref):
                    # print " SSS ",Blockdata["myid"],j
                    done_on_node = j % number_of_proc
                    if done_on_node != main_node:
                        if myid == main_node:
                            mem_len = mpi.mpi_recv(
                                1,
                                mpi.MPI_INT,
                                done_on_node,
                                sp_global_def.SPARX_MPI_TAG_UNIVERSAL,
                                comm,
                            )
                            mem_len = int(mem_len[0])
                            members = mpi.mpi_recv(
                                mem_len,
                                mpi.MPI_INT,
                                done_on_node,
                                sp_global_def.SPARX_MPI_TAG_UNIVERSAL,
                                comm,
                            )
                            members = list(map(int, members))
                            refi[j].set_attr_dict(
                                {"members": members, "n_objects": mem_len}
                            )
                        elif myid == done_on_node:
                            members = refi[j].get_attr("members")
                            mpi.mpi_send(
                                len(members),
                                1,
                                mpi.MPI_INT,
                                main_node,
                                sp_global_def.SPARX_MPI_TAG_UNIVERSAL,
                                comm,
                            )
                            mpi.mpi_send(
                                members,
                                len(members),
                                mpi.MPI_INT,
                                main_node,
                                sp_global_def.SPARX_MPI_TAG_UNIVERSAL,
                                comm,
                            )

                if myid == main_node:
                    #  compare with previous
                    totprevious = 0.0
                    totcurrent = 0.0
                    common = 0.0
                    for j in range(numref):
                        totprevious += len(previous_members[j])
                        members = set(refi[j].get_attr("members"))
                        totcurrent += len(members)
                        common += len(previous_members[j].intersection(members))
                        previous_members[j] = members
                        # print "  memebers  ",j,len(members)
                    agreement = old_div(
                        common, float(totprevious + totcurrent - common)
                    )
                    j = agreement - previous_agreement
                    if (agreement > 0.5) and (j > 0.0) and (j < 0.05):
                        terminate = 1
                    previous_agreement = agreement
                    sp_global_def.sxprint(
                        ">>>  Assignment agreement with previous iteration  %5.1f"
                        % (agreement * 100),
                        "   ",
                        time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()),
                    )
                terminate = sp_utilities.bcast_number_to_all(
                    terminate, source_node=main_node
                )

            if check_stability and ((main_iter == max_iter) or (terminate == 1)):
                #  gather all pixers and print a histogram
                gpixer = sp_utilities.wrap_mpi_gatherv(gpixer, main_node, comm)
                if my_abs_id == main_node and color == 0:
                    lhist = 12
                    region, histo = sp_statistics.hist_list(gpixer, lhist)
                    sp_global_def.sxprint(
                        "\n=== Histogram of average within-class pixel errors prior to class pruning ==="
                    )
                    for lhx in range(lhist):
                        sp_global_def.sxprint(
                            "     %10.3f     %7d" % (region[lhx], histo[lhx])
                        )
                    sp_global_def.sxprint(
                        "=============================================================================\n"
                    )
                del gpixer

        # end of within group alignment
        mpi.mpi_barrier(comm)

    if myid == main_node:
        i = [refi[j].get_attr("n_objects") for j in range(numref)]
        lhist = max(12, old_div(numref, 2))
        region, histo = sp_statistics.hist_list(i, lhist)
        sp_global_def.sxprint(
            "\n=== Histogram of group sizes ================================================"
        )
        for lhx in range(lhist):
            sp_global_def.sxprint("     %10.1f     %7d" % (region[lhx], histo[lhx]))
        sp_global_def.sxprint(
            "=============================================================================\n"
        )

    mpi.mpi_barrier(comm)

    return refi


def do_generation(
    main_iter, generation_iter, target_nx, target_xr, target_yr, target_radius, options
):
    """
	Perform one iteration of ISAC processing within current generation.

	Args:
		main_iter (int): Number of SAC main iterations, i.e., the number of
			runs of cluster alignment for stability evaluation in SAC.
			[Default: 3]

		generation_iter (int): Number of iterations in the current generation.

		target_nx (int): Target particle image size. This is the actual image
			size on which ISAC will process data. Images will be scaled
			according to target particle radius and pruned/padded to achieve
			target_nx size. If xr > 0 (see below), the final image size for
			ISAC processing equals target_nx + xr -1.
			[Default: 76]

		target_xr (int): x-range of translational search during alignment.
			[Default: 1]

		target_yr (int): y-range of translational search during alignment.
			[Default: 1]

		target_radius (int): Target particle radius used when ISAC processes
			the data. Images will be scaled to conform to this value.
			[Default: 29]

		options (options object): Provided by the Python OptionParser. This
			structure contains all command line options (option "--value" is
			accessed by options.value).

	Returns:
		keepdoing_main (bool): Indicates the main ISAC iteration should stop.

		keepdoing_generation (bool): Indicates the iterations within this
			generation should stop.
	"""

    global Blockdata

    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    if Blockdata["myid"] == Blockdata["main_node"]:
        plist = sp_utilities.read_text_file(
            os.path.join(
                Blockdata["masterdir"],
                "main%03d" % main_iter,
                "generation%03d" % (generation_iter - 1),
                "to_process_next_%03d_%03d.txt" % (main_iter, generation_iter - 1),
            )
        )
        nimastack = len(plist)

    else:
        plist = 0
        nimastack = 0

    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    nimastack = sp_utilities.bcast_number_to_all(
        nimastack, source_node=Blockdata["main_node"], mpi_comm=mpi.MPI_COMM_WORLD
    )

    # Bcast plist to all zero CPUs
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    if Blockdata["subgroup_myid"] > -1:
        # First check number of nodes, if only one, no reduction necessary.
        # print  "  subgroup_myid   ",Blockdata["subgroup_myid"],Blockdata["no_of_groups"],nimastack
        if Blockdata["no_of_groups"] > 1:
            plist = sp_utilities.bcast_list_to_all(
                plist,
                Blockdata["subgroup_myid"],
                source_node=0,
                mpi_comm=Blockdata["subgroup_comm"],
            )
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    # reserve buffers
    disp_unit = numpy.dtype("f4").itemsize
    size_of_one_image = target_nx * target_nx
    orgsize = (
        nimastack * size_of_one_image
    )  # This is number of projections to be computed simultaneously times their size

    if Blockdata["myid_on_node"] == 0:
        size = orgsize
    else:
        size = 0

    win_sm, base_ptr = mpi.mpi_win_allocate_shared(
        size * disp_unit, disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"]
    )
    size = orgsize

    if Blockdata["myid_on_node"] != 0:
        base_ptr, = mpi.mpi_win_shared_query(win_sm, mpi.MPI_PROC_NULL)


    # buffer = numpy.frombuffer(
    #     numpy.core.multiarray.int_asbuffer(base_ptr, size * disp_unit), dtype="f4"
    # )

    """
    Comments from Adnan:
    numpy.core.mutliarray.int_asbuffer function is not available anymore . That is why a different alternative is tried here.
    """
    ptr = ctypes.cast(base_ptr, ctypes.POINTER(ctypes.c_int * size))
    buffer = numpy.frombuffer(ptr.contents, dtype="f4")

    buffer = buffer.reshape(nimastack, target_nx, target_nx)

    emnumpy2 = EMAN2_cppwrap.EMNumPy()
    bigbuffer = emnumpy2.register_numpy_to_emdata(buffer)

    if Blockdata["myid_on_node"] == 0:
        #  read data on process 0 of each node
        # print "  READING DATA FIRST :",Blockdata["myid"],Blockdata["stack_ali2d"],len(plist)
        for i in range(nimastack):
            bigbuffer.insert_clip(
                sp_utilities.get_im(Blockdata["stack_ali2d"], plist[i]), (0, 0, i)
            )
        del plist

    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    alldata = [None] * nimastack
    emnumpy3 = [None] * nimastack

    msk = sp_utilities.model_blank(target_nx, target_nx, 1, 1)
    for i in range(nimastack):
        newpoint = base_ptr + i * size_of_one_image * disp_unit
        pointer_location = ctypes.cast(newpoint, ctypes.POINTER(ctypes.c_int * size_of_one_image))
        img_buffer = numpy.frombuffer(pointer_location.contents, dtype="f4")
        img_buffer = img_buffer.reshape(target_nx, target_nx)
        emnumpy3[i] = EMAN2_cppwrap.EMNumPy()
        alldata[i] = emnumpy3[i].register_numpy_to_emdata(img_buffer)
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    # former options
    indep_run = 0
    match_first = 0
    thld_grp = 0
    match_second = 0
    max_round = 0
    dummy_main_iter = 0

    if Blockdata["myid"] == 0:
        sp_global_def.sxprint(
            "*************************************************************************************"
        )
        sp_global_def.sxprint(
            "     Main iteration: %3d,  Generation: %3d. "
            % (main_iter, generation_iter)
            + "   "
            + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        )

    ave, all_params = iter_isac_pap(
        alldata,
        options.ir,
        target_radius,
        options.rs,
        target_xr,
        target_yr,
        options.ts,
        options.maxit,
        False,
        1.0,
        options.dst,
        options.FL,
        options.FH,
        options.FF,
        options.init_iter,
        dummy_main_iter,
        options.iter_reali,
        match_first,
        max_round,
        match_second,
        options.stab_ali,
        options.thld_err,
        indep_run,
        thld_grp,
        options.img_per_grp,
        generation_iter,
        False,
        random_seed=options.rand_seed,
        new=False,
    )  # options.new)

    #  Clean the stack
    # delay(Blockdata["myid"],"  PROCESSING DONE")
    mpi.mpi_win_free(win_sm)
    emnumpy2.unregister_numpy_from_emdata()
    del emnumpy2
    for i in range(nimastack):
        emnumpy3[i].unregister_numpy_from_emdata()
    del alldata
    # delay(Blockdata["myid"],"  CLEANED")

    if Blockdata["myid"] == Blockdata["main_node"]:
        #  How many averages alreay exist
        if os.path.exists(os.path.join(Blockdata["masterdir"], "class_averages.hdf")):
            nave_exist = EMAN2_cppwrap.EMUtil.get_image_count(
                os.path.join(Blockdata["masterdir"], "class_averages.hdf")
            )
        else:
            nave_exist = 0
        #  Read all parameters table from masterdir
        all_parameters = sp_utilities.read_text_row(
            os.path.join(Blockdata["masterdir"], "all_parameters.txt")
        )
        plist = sp_utilities.read_text_file(
            os.path.join(
                Blockdata["masterdir"],
                "main%03d" % main_iter,
                "generation%03d" % (generation_iter - 1),
                "to_process_next_%03d_%03d.txt" % (main_iter, generation_iter - 1),
            )
        )
        # print "****************************************************************************************************",os.path.join(Blockdata["masterdir"], "main%03d"%main_iter, \
        # 	"generation%03d"%(generation_iter-1), "to_process_next_%03d_%03d.txt"%(main_iter,generation_iter-1))
        j = 0
        good = []
        bad = []
        new_good_classes = 0
        for i, q in enumerate(ave):
            #  Convert local numbering to absolute numbering of images
            local_members = q.get_attr("members")
            members = [plist[l] for l in local_members]
            q.write_image(
                os.path.join(
                    Blockdata["masterdir"],
                    "main%03d" % main_iter,
                    "generation%03d" % generation_iter,
                    "original_class_averages_%03d_%03d.hdf"
                    % (main_iter, generation_iter),
                ),
                i,
            )
            q.set_attr("members", members)
            q.write_image(
                os.path.join(
                    Blockdata["masterdir"],
                    "main%03d" % main_iter,
                    "generation%03d" % generation_iter,
                    "class_averages_%03d_%03d.hdf" % (main_iter, generation_iter),
                ),
                i,
            )
            if len(members) > options.minimum_grp_size:
                new_good_classes += 1
                good += members
                q.write_image(
                    os.path.join(
                        Blockdata["masterdir"],
                        "main%03d" % main_iter,
                        "generation%03d" % generation_iter,
                        "good_class_averages_%03d_%03d.hdf"
                        % (main_iter, generation_iter),
                    ),
                    j,
                )
                q.write_image(
                    os.path.join(Blockdata["masterdir"], "class_averages.hdf"),
                    j + nave_exist,
                )
                j += 1

                # We have to update all parameters table
                for l, m in enumerate(members):
                    #  I had to remove it as in case of restart there will be conflicts
                    # if( all_parameters[m][-1] > -1):
                    # 	print "  CONFLICT !!!"
                    # 	exit()
                    all_parameters[m] = all_params[local_members[l]]
            else:
                bad += members

        if len(good) > 0:
            sp_utilities.write_text_row(
                all_parameters,
                os.path.join(Blockdata["masterdir"], "all_parameters.txt"),
            )
            good.sort()
            #  Add currently assigned images to the overall list
            if os.path.exists(
                os.path.join(Blockdata["masterdir"], "processed_images.txt")
            ):
                lprocessed = good + sp_utilities.read_text_file(
                    os.path.join(Blockdata["masterdir"], "processed_images.txt")
                )
                lprocessed.sort()
                sp_utilities.write_text_file(
                    lprocessed,
                    os.path.join(Blockdata["masterdir"], "processed_images.txt"),
                )
            else:
                sp_utilities.write_text_file(
                    good, os.path.join(Blockdata["masterdir"], "processed_images.txt")
                )

        if len(bad) > 0:
            bad.sort()
            sp_utilities.write_text_file(
                bad,
                os.path.join(
                    Blockdata["masterdir"],
                    "main%03d" % main_iter,
                    "generation%03d" % (generation_iter),
                    "to_process_next_%03d_%03d.txt" % (main_iter, generation_iter),
                ),
            )

            if (int(len(bad) * 1.2) < 2 * options.img_per_grp) or (
                (generation_iter == 1) and (new_good_classes <= options.delta_good)
            ):
                #  Insufficient number of images to keep processing bad set
                #    or
                #  Program cannot produce any good averages from what is left at the beginning of new main
                try:
                    lprocessed = sp_utilities.read_text_file(
                        os.path.join(Blockdata["masterdir"], "processed_images.txt")
                    )
                except:
                    lprocessed = []
                nprocessed = len(lprocessed)
                leftout = sorted(
                    list(set(range(Blockdata["total_nima"])) - set(lprocessed))
                )
                sp_utilities.write_text_file(
                    leftout,
                    os.path.join(Blockdata["masterdir"], "not_processed_images.txt"),
                )
                # Check whether what remains can be still processed in a new main interation
                if (len(leftout) < 2 * options.img_per_grp) or (
                    (generation_iter == 1) and (new_good_classes <= options.delta_good)
                ):
                    #    if the the number of remaining all bad too low full stop
                    keepdoing_main = False
                    keepdoing_generation = False
                    cmd = "{} {}".format(
                        "touch",
                        os.path.join(
                            Blockdata["masterdir"],
                            "main%03d" % main_iter,
                            "generation%03d" % generation_iter,
                            "finished",
                        ),
                    )
                    junk = sp_utilities.cmdexecute(cmd)
                    cmd = "{} {}".format(
                        "touch",
                        os.path.join(
                            Blockdata["masterdir"], "main%03d" % main_iter, "finished"
                        ),
                    )
                    junk = sp_utilities.cmdexecute(cmd)
                    cmd = "{} {}".format(
                        "touch", os.path.join(Blockdata["masterdir"], "finished")
                    )
                    junk = sp_utilities.cmdexecute(cmd)
                    sp_global_def.sxprint(
                        "*         There are no more images to form averages, program finishes     "
                        + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
                        + "     *"
                    )
                else:
                    #  Will have to increase main, which means putting all bad left as new good,
                    keepdoing_main = True
                    keepdoing_generation = False
                    #  Will have to increase main, which means putting all bad left as new good,
                    cmd = "{} {}".format(
                        "touch",
                        os.path.join(
                            Blockdata["masterdir"],
                            "main%03d" % main_iter,
                            "generation%03d" % generation_iter,
                            "finished",
                        ),
                    )
                    junk = sp_utilities.cmdexecute(cmd)
                    cmd = "{} {}".format(
                        "touch",
                        os.path.join(
                            Blockdata["masterdir"], "main%03d" % main_iter, "finished"
                        ),
                    )
                    junk = sp_utilities.cmdexecute(cmd)
            else:
                keepdoing_main = True
                keepdoing_generation = True
        else:
            keepdoing_main = False
            keepdoing_generation = False
        # print "****************************************************************************************************",keepdoing_main,keepdoing_generation

    else:
        keepdoing_main = False
        keepdoing_generation = False
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    keepdoing_main = sp_utilities.bcast_number_to_all(
        keepdoing_main, source_node=Blockdata["main_node"], mpi_comm=mpi.MPI_COMM_WORLD
    )
    keepdoing_generation = sp_utilities.bcast_number_to_all(
        keepdoing_generation,
        source_node=Blockdata["main_node"],
        mpi_comm=mpi.MPI_COMM_WORLD,
    )

    return keepdoing_main, keepdoing_generation


# ======================================================================[ main ]


def parse_parameters(prog_name, usage, args):
    prog_name = os.path.basename(prog_name)
    parser = optparse.OptionParser(usage, version=sp_global_def.SPARXVERSION)

    # ISAC command line parameters (public)
    parser.add_option(
        "--radius",
        type="int",
        help="particle radius: there is no default, a sensible number has to be provided, units - pixels (default required int)",
    )
    parser.add_option(
        "--target_radius",
        type="int",
        default=29,
        help="target particle radius: actual particle radius on which isac will process data. Images will be shrinked/enlarged to achieve this radius (default 29)",
    )
    parser.add_option(
        "--target_nx",
        type="int",
        default=76,
        help="target particle image size: actual image size on which isac will process data. Images will be shrinked/enlarged according to target particle radius and then cut/padded to achieve target_nx size. When xr > 0, the final image size for isac processing is 'target_nx + xr - 1'  (default 76)",
    )
    parser.add_option(
        "--img_per_grp",
        type="int",
        default=200,
        help="number of images per class (maximum group size, also defines number of classes K=(total number of images)/img_per_grp (default 200)",
    )
    parser.add_option(
        "--minimum_grp_size",
        type="int",
        default=60,
        help="minimum size of class (default 60)",
    )
    parser.add_option(
        "--CTF",
        action="store_true",
        default=False,
        help="apply phase-flip for CTF correction: if set the data will be phase-flipped using CTF information included in image headers (default False)",
    )
    parser.add_option(
        "--VPP",
        action="store_true",
        default=False,
        help="Phase Plate data (default False)",
    )
    parser.add_option(
        "--ir",
        type="int",
        default=1,
        help="inner ring: of the resampling to polar coordinates. units - pixels (default 1)",
    )
    parser.add_option(
        "--rs",
        type="int",
        default=1,
        help="ring step: of the resampling to polar coordinates. units - pixels (default 1)",
    )
    parser.add_option(
        "--xr",
        type="int",
        default=1,
        help="x range: of translational search. By default, set by the program. (default 1)",
    )
    parser.add_option(
        "--yr",
        type="int",
        default=-1,
        help="y range: of translational search. By default, same as xr. (default -1)",
    )
    parser.add_option(
        "--ts",
        type="float",
        default=1.0,
        help="search step: of translational search: units - pixels (default 1.0)",
    )
    parser.add_option(
        "--maxit",
        type="int",
        default=30,
        help="number of iterations for reference-free alignment (default 30)",
    )
    parser.add_option(
        "--center_method",
        type="int",
        default=-1,
        help="method for centering: of global 2D average during initial prealignment of data (0 : no centering; -1 : average shift method; please see center_2D in utilities.py for methods 1-7) (default -1)",
    )
    parser.add_option(
        "--dst",
        type="float",
        default=90.0,
        help="discrete angle used in within group alignment (default 90.0)",
    )
    parser.add_option(
        "--FL",
        type="float",
        default=0.2,
        help="lowest stopband: frequency used in the tangent filter (default 0.2)",
    )
    parser.add_option(
        "--FH",
        type="float",
        default=0.45,
        help="highest stopband: frequency used in the tangent filter (default 0.45)",
    )
    parser.add_option(
        "--FF",
        type="float",
        default=0.2,
        help="fall-off of the tangent filter (default 0.2)",
    )
    parser.add_option(
        "--init_iter",
        type="int",
        default=7,
        help="Maximum number of Generation iterations performed for a given subset (default 7)",
    )
    parser.add_option(
        "--iter_reali",
        type="int",
        default=1,
        help="SAC stability check interval: every iter_reali iterations of SAC stability checking is performed (default 1)",
    )
    parser.add_option(
        "--stab_ali",
        type="int",
        default=5,
        help="number of alignments when checking stability (default 5)",
    )
    parser.add_option(
        "--thld_err",
        type="float",
        default=0.7,
        help="threshold of pixel error when checking stability: equals root mean square of distances between corresponding pixels from set of found transformations and theirs average transformation, depends linearly on square of radius (parameter target_radius). units - pixels. (default 0.7)",
    )
    parser.add_option(
        "--restart",
        type="int",
        default="-1",
        help="0: restart ISAC2 after last completed main iteration (meaning there is file >finished< in it.  k: restart ISAC2 after k'th main iteration (It has to be completed, meaning there is file >finished< in it. Higer iterations will be removed.)  Default: no restart",
    )
    parser.add_option(
        "--rand_seed",
        type="int",
        help="random seed set before calculations: useful for testing purposes. By default, total randomness (type int)",
    )

    parser.add_option(
        "--skip_prealignment",
        action="store_true",
        default=False,
        help="skip pre-alignment step: to be used if images are already centered. 2dalignment directory will still be generated but the parameters will be zero. (default False)",
    )

    parser.add_option(
        "--filament_width",
        type="int",
        default=-1,
        help="When this is set to a non-default value helical data is assumed in which case particle images will be masked with a rectangular mask. (Default: -1)",
    )
    parser.add_option(
        "--filament_mask_ignore",
        action="store_true",
        default=False,
        help="Only relevant when parameter '--filament_width' is set. When set to False a rectangular mask is used to (a) normalize and (b) to mask the particle images. The latter can be disabled by setting this flag to True. (Default: False)",
    )
    parser.add_option(
        "--skip_ordering",
        action="store_true",
        default=False,
        help="Skip the ordered class_averages creation. (Default: False)",
    )
    parser.add_option(
        "--delta_good",
        type="int",
        default=0,
        help="Convergence criteria for ISAC. As soon as the number of new good classes is lower that this value ISAC stops.")

    return parser.parse_args(args)


def run(args):
    # ------------------------------------------------------[ command line parameters ]

    usage = (
        sys.argv[0]
        + " stack_file  [output_directory] --radius=particle_radius"
        + " --img_per_grp=img_per_grp --CTF <The remaining parameters are"
        + " optional --ir=ir --rs=rs --xr=xr --yr=yr --ts=ts --maxit=maxit"
        + " --dst=dst --FL=FL --FH=FH --FF=FF --init_iter=init_iter"
        + " --iter_reali=iter_reali --stab_ali=stab_ali --thld_err=thld_err"
        + " --rand_seed=rand_seed>"
    )

    options, args = parse_parameters(
        sys.argv[0], usage, args
    )  # NOTE: output <args> != input <args>

    # after parsing, the only remaining args should be path to input & output folders
    if len(args) > 2:
        print("usage: " + usage)
        print("Please run '" + sys.argv[0] + " -h' for detailed options")
        sys.exit()
    elif len(args) == 2:
        Blockdata["stack"] = args[0]
        Blockdata["masterdir"] = args[1]
    elif len(args) == 1:
        Blockdata["stack"] = args[0]
        Blockdata["masterdir"] = ""

    options.new = False

    # check required options
    required_option_list = ["radius"]

    for required_option in required_option_list:
        if not options.__dict__[required_option]:
            print("\n ==%s== mandatory option is missing.\n" % required_option)
            print("Please run '" + sys.argv[0] + " -h' for detailed options")
            return 1

    # sanity check: make sure the minimum group size is smaller than the actual group size
    if options.minimum_grp_size > options.img_per_grp:
        if Blockdata["myid"] == Blockdata["main_node"]:
            print(
                "\nERROR! Minimum group size ("
                + str(options.minimum_grp_size)
                + ") is larger than the actual group size ("
                + str(options.img_per_grp)
                + "). Oh dear :(\n"
            )
        return 1

    # TODO: what does this do?
    if sp_global_def.CACHE_DISABLE:
        sp_utilities.disable_bdb_cache()
    sp_global_def.BATCH = True

    # ------------------------------------------------------[ master directory setup ]

    # get mpi id values (NOTE: these are not used consistently throughout the code)
    main_node = Blockdata["main_node"]
    myid = Blockdata["myid"]
    nproc = Blockdata["nproc"]

    # main process creates the master directory
    if Blockdata["myid"] == Blockdata["main_node"]:
        if Blockdata["masterdir"] == "":
            timestring = time.strftime("_%d_%b_%Y_%H_%M_%S", time.localtime())
            Blockdata["masterdir"] = "isac_directory" + timestring
            li = len(Blockdata["masterdir"])
            cmd = "{} {}".format("mkdir -p", Blockdata["masterdir"])
            junk = sp_utilities.cmdexecute(cmd)
        else:
            if not os.path.exists(Blockdata["masterdir"]):
                cmd = "{} {}".format("mkdir -p", Blockdata["masterdir"])
                junk = sp_utilities.cmdexecute(cmd)
            li = 0
        sp_global_def.write_command(Blockdata["masterdir"])
    else:
        li = 0

    # main process broadcasts length of master directory string
    li = mpi.mpi_bcast(li, 1, mpi.MPI_INT, Blockdata["main_node"], mpi.MPI_COMM_WORLD)[
        0
    ]

    # main process broadcasts path to ISAC master directory
    if li > 0:
        Blockdata["masterdir"] = mpi.mpi_bcast(
            Blockdata["masterdir"],
            li,
            mpi.MPI_CHAR,
            Blockdata["main_node"],
            mpi.MPI_COMM_WORLD,
        )
        Blockdata["masterdir"] = string.join(Blockdata["masterdir"], "")

    # add stack_ali2d path to blockdata
    Blockdata["stack_ali2d"] = "bdb:" + os.path.join(
        Blockdata["masterdir"], "stack_ali2d"
    )

    if myid == main_node:
        sp_global_def.sxprint(
            "****************************************************************************************************"
        )
        sp_global_def.sxprint(
            "*                                                                                                  *"
        )
        sp_global_def.sxprint(
            "* ISAC (Iterative Stable Alignment and Clustering)   "
            + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
            + "                     *"
        )
        sp_global_def.sxprint(
            "* By Zhengfan Yang, Jia Fang, Francisco Asturias and Pawel A. Penczek                              *"
        )
        sp_global_def.sxprint(
            "*                                                                                                  *"
        )
        sp_global_def.sxprint(
            '* REFERENCE: Z. Yang, J. Fang, J. Chittuluru, F. J. Asturias and P. A. Penczek, "Iterative Stable  *'
        )
        sp_global_def.sxprint(
            '*            Alignment and Clustering of 2D Transmission Electron Microscope Images",              *'
        )
        sp_global_def.sxprint(
            "*            Structure 20, 237-247, February 8, 2012.                                              *"
        )
        sp_global_def.sxprint(
            "*                                                                                                  *"
        )
        sp_global_def.sxprint(
            "* Last updated: 05/30/2017 PAP                                                                     *"
        )
        sp_global_def.sxprint(
            "****************************************************************************************************"
        )
        EMAN2.Util.version()
        sp_global_def.sxprint(
            "****************************************************************************************************"
        )
        sys.stdout.flush()

        i = "  "
        for a in sys.argv:
            i += a + "  "
        sp_global_def.sxprint("* shell line command: ")
        sp_global_def.sxprint(i)
        sp_global_def.sxprint(
            "*                                                                                                  *"
        )
        sp_global_def.sxprint("* Master directory: %s" % Blockdata["masterdir"])
        sp_global_def.sxprint(
            "****************************************************************************************************"
        )
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    # ------------------------------------------------------[ zero group ]

    create_zero_group()

    # ------------------------------------------------------[ gather image parameters ]

    if options.CTF and options.VPP:
        sp_global_def.ERROR("Options CTF and VPP cannot be used together", myid=myid)

    # former options
    indep_run = 0
    match_first = 0
    match_second = 0
    max_round = 0
    main_iter = 0

    radi = options.radius
    target_radius = options.target_radius
    target_nx = options.target_nx  # dimension of particle images after downscaling
    center_method = options.center_method

    if radi < 1:
        sp_global_def.ERROR("Particle radius has to be provided!", myid=myid)

    target_xr = options.xr
    target_nx += target_xr - 1  # subtract one, which is default

    if options.yr == -1:
        target_yr = options.xr
    else:
        target_yr = options.yr

    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    # get total number of images (nima) and broadcast
    if myid == main_node:
        Blockdata["total_nima"] = EMAN2_cppwrap.EMUtil.get_image_count(
            Blockdata["stack"]
        )
    else:
        Blockdata["total_nima"] = 0

    Blockdata["total_nima"] = sp_utilities.bcast_number_to_all(
        Blockdata["total_nima"], source_node=main_node
    )

    nxrsteps = 4

    # ------------------------------------------------------[ prepare ISAC loop to run from scratch or continue ]

    error = 0
    if Blockdata["myid"] == Blockdata["main_node"]:

        # fresh start
        if not os.path.exists(
            os.path.join(Blockdata["masterdir"], "main001", "generation000")
        ):
            #  NOTE: we do not create processed_images.txt selection file as it has to be initially empty
            #  we do, however, initialize all parameters with empty values
            sp_utilities.write_text_row(
                [[0.0, 0.0, 0.0, -1] for i in range(Blockdata["total_nima"])],
                os.path.join(Blockdata["masterdir"], "all_parameters.txt"),
            )
            if options.restart > -1:
                error = 1

        # continue ISAC from a previous run
        else:

            # restart ISAC from the last finished iteration
            if options.restart == 0:
                keepdoing_main = True
                main_iter = 0
                # go through the main iteration folders and find the first not to contain the "finished" file, tehn start from there
                while keepdoing_main:
                    main_iter += 1
                    if os.path.exists(
                        os.path.join(Blockdata["masterdir"], "main%03d" % main_iter)
                    ):
                        if not os.path.exists(
                            os.path.join(
                                Blockdata["masterdir"],
                                "main%03d" % main_iter,
                                "finished",
                            )
                        ):
                            sp_utilities.cmdexecute(
                                "{} {}".format(
                                    "rm -rf",
                                    os.path.join(
                                        Blockdata["masterdir"], "main%03d" % main_iter
                                    ),
                                )
                            )
                            keepdoing_main = False
                    else:
                        keepdoing_main = False

            # restart ISAC from a specfied iteration number
            else:
                main_iter = options.restart
                if not os.path.exists(
                    os.path.join(
                        Blockdata["masterdir"], "main%03d" % main_iter, "finished"
                    )
                ):
                    error = 2
                else:
                    keepdoing_main = True
                    main_iter += 1
                    # when continuing from iteration <n> remove all existing directories for iterations <n+i>
                    while keepdoing_main:
                        if os.path.exists(
                            os.path.join(Blockdata["masterdir"], "main%03d" % main_iter)
                        ):
                            sp_utilities.cmdexecute(
                                "{} {}".format(
                                    "rm -rf",
                                    os.path.join(
                                        Blockdata["masterdir"], "main%03d" % main_iter
                                    ),
                                )
                            )
                            main_iter += 1
                        else:
                            keepdoing_main = False

            # if we're re-starting from an iteration folder containing the "finished" file, remove that file; we're just getting started!
            if os.path.exists(os.path.join(Blockdata["masterdir"], "finished")):
                sp_utilities.cmdexecute(
                    "{} {}".format(
                        "rm -rf", os.path.join(Blockdata["masterdir"], "finished")
                    )
                )

    error = sp_utilities.bcast_number_to_all(error, source_node=Blockdata["main_node"])
    if error == 1:
        sp_global_def.ERROR(
            "Good Sir or Madam, for thine restart value %d there exists no ISAC directory for thine continuation. Please doth choose a restart value of -1 to start a fresh run or 0 to continue the iteration last completed."
            % options.restart,
            myid=Blockdata["myid"],
        )
    elif error == 2:
        sp_global_def.ERROR(
            "Cannot restart from unfinished main iteration number %d. To automatically detect the latest ISAC directory for continuation please choose a restart value of 0."
            % main_iter,
            myid=Blockdata["myid"],
        )

    # everyone waits for the root process to run the above checks
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    # ------------------------------------------------------[ initial 2D alignment (centering) ]

    init2dir = os.path.join(Blockdata["masterdir"], "2dalignment")

    if not checkitem(os.path.join(init2dir, "Finished_initial_2d_alignment.txt")):

        if myid == 0:

            #  Create output directory
            log2d = sp_logger.Logger(sp_logger.BaseLogger_Files())
            log2d.prefix = os.path.join(init2dir)
            cmd = "mkdir -p " + log2d.prefix
            outcome = subprocess.call(cmd, shell=True)
            log2d.prefix += "/"
        else:
            outcome = 0
            log2d = None

        if myid == main_node:
            a = sp_utilities.get_im(Blockdata["stack"])
            nnxo = a.get_xsize()
        else:
            nnxo = 0
        nnxo = sp_utilities.bcast_number_to_all(nnxo, source_node=main_node)

        image_start, image_end = sp_applications.MPI_start_end(
            Blockdata["total_nima"], nproc, myid
        )

        original_images = EMAN2_cppwrap.EMData.read_images(
            Blockdata["stack"], list(range(image_start, image_end))
        )

        if options.VPP:
            ntp = len(sp_fundamentals.rops_table(original_images[0]))
            rpw = [0.0] * ntp
            for q in original_images:
                tpw = sp_fundamentals.rops_table(q)
                for i in range(ntp):
                    rpw[i] += numpy.sqrt(tpw[i])
            del tpw
            rpw = mpi.mpi_reduce(
                rpw, ntp, mpi.MPI_FLOAT, mpi.MPI_SUM, main_node, mpi.MPI_COMM_WORLD
            )
            if myid == 0:
                rpw = [float(old_div(Blockdata["total_nima"], q)) for q in rpw]
                rpw[0] = 1.0
                sp_utilities.write_text_file(
                    rpw, os.path.join(Blockdata["masterdir"], "rpw.txt")
                )
            else:
                rpw = []
            rpw = sp_utilities.bcast_list_to_all(
                rpw, myid, source_node=main_node, mpi_comm=mpi.MPI_COMM_WORLD
            )
            for i in range(len(original_images)):
                original_images[i] = sp_filter.filt_table(original_images[i], rpw)
            del rpw

        else:
            if myid == 0:
                ntp = len(sp_fundamentals.rops_table(original_images[0]))
                sp_utilities.write_text_file(
                    [0.0] * ntp, os.path.join(Blockdata["masterdir"], "rpw.txt")
                )

        if options.skip_prealignment:
            params2d = [[0.0, 0.0, 0.0, 0] for i in range(image_start, image_end)]
        else:
            #  We assume the target radius will be 29, and xr = 1.
            shrink_ratio = old_div(float(target_radius), float(radi))

            for im in range(len(original_images)):
                if shrink_ratio != 1.0:
                    original_images[im] = sp_fundamentals.resample(
                        original_images[im], shrink_ratio
                    )

            nx = original_images[0].get_xsize()
            # nx = int(nx*shrink_ratio + 0.5)

            txrm = old_div((nx - 2 * (target_radius + 1)), 2)
            if txrm < 0:
                sp_global_def.ERROR(
                    "Radius of the structure larger than the window data size permits %d"
                    % (radi),
                    myid=myid,
                )
            if old_div(txrm, nxrsteps) > 0:
                tss = ""
                txr = ""
                while old_div(txrm, nxrsteps) > 0:
                    tts = old_div(txrm, nxrsteps)
                    tss += "  %d" % tts
                    txr += "  %d" % (tts * nxrsteps)
                    txrm = old_div(txrm, 2)
            else:
                tss = "1"
                txr = "%d" % txrm

            # perform the 2D alignment
            if Blockdata["myid"] == 0:
                sp_global_def.sxprint(
                    "* 2D alignment   "
                    + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
                )

            params2d = sp_applications.ali2d_base(
                original_images,
                init2dir,
                None,
                1,
                target_radius,
                1,
                txr,
                txr,
                tss,
                False,
                90.0,
                center_method,
                14,
                options.CTF,
                1.0,
                False,
                "ref_ali2d",
                "",
                log2d,
                nproc,
                myid,
                main_node,
                mpi.MPI_COMM_WORLD,
                write_headers=False,
            )

            del original_images

            for i in range(len(params2d)):
                alpha, sx, sy, mirror = sp_utilities.combine_params2(
                    0, params2d[i][1], params2d[i][2], 0, -params2d[i][0], 0, 0, 0
                )
                sx = old_div(sx, shrink_ratio)
                sy = old_div(sy, shrink_ratio)
                params2d[i][0] = 0.0
                params2d[i][1] = sx
                params2d[i][2] = sy
                params2d[i][3] = 0
                # util.set_params2D(aligned_images[i],[0.0, sx,sy,0.,1.0])

        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        tmp = params2d[:]
        tmp = sp_utilities.wrap_mpi_gatherv(tmp, main_node, mpi.MPI_COMM_WORLD)
        if myid == main_node:
            if options.skip_prealignment:
                sp_global_def.sxprint("=========================================")
                sp_global_def.sxprint(
                    "There is no alignment step, '%s' params are set to zero for later use."
                    % os.path.join(init2dir, "initial2Dparams.txt")
                )
                sp_global_def.sxprint("=========================================")
            sp_utilities.write_text_row(
                tmp, os.path.join(init2dir, "initial2Dparams.txt")
            )
        del tmp
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

        #  We assume the target image size will be target_nx, radius will be 29, and xr = 1.
        #  Note images can be also padded, in which case shrink_ratio > 1.
        shrink_ratio = old_div(float(target_radius), float(radi))

        aligned_images = EMAN2_cppwrap.EMData.read_images(
            Blockdata["stack"], list(range(image_start, image_end))
        )
        nx = aligned_images[0].get_xsize()
        nima = len(aligned_images)
        newx = int(nx * shrink_ratio + 0.5)

        while not os.path.exists(os.path.join(init2dir, "initial2Dparams.txt")):
            time.sleep(1)
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

        params = sp_utilities.read_text_row(
            os.path.join(init2dir, "initial2Dparams.txt")
        )
        params = params[image_start:image_end]

        # --------------------------------------------------[ defocus value correction ]

        # for phase plate data: read a given rotational power spectrum (1D amplitude profile)
        if options.VPP:
            if myid == 0:
                rpw = sp_utilities.read_text_file(
                    os.path.join(Blockdata["masterdir"], "rpw.txt")
                )
            else:
                rpw = [0.0]
            rpw = sp_utilities.bcast_list_to_all(
                rpw, myid, source_node=main_node, mpi_comm=mpi.MPI_COMM_WORLD
            )

        # if we're not looking at filament images we create a circular mask that we can use for all particles
        if options.filament_width == -1:
            mask = sp_utilities.model_circle(radi, nx, nx)

        # defocus value correction for all images
        if Blockdata["myid"] == main_node:
            sp_global_def.sxprint("Apply CTF")
        for im in range(nima):
            # create custom mask per particle in case we're processing filament images
            if options.filament_width != -1:
                mask = sp_utilities.model_rotated_rectangle2D(
                    radius_long=old_div(
                        int(numpy.sqrt(2 * nx ** 2)), 2
                    ),  # use a length that will be guaranteed to cross the whole image
                    radius_short=int(options.filament_width * 0.5),
                    # use a conservative, slightly larger mask since filaments might not be aligned as well as they could be
                    nx=nx,
                    ny=nx,
                    angle=aligned_images[im].get_attr("segment_angle"),
                )
            # subtract mean of the within-mask area
            st = EMAN2_cppwrap.Util.infomask(aligned_images[im], mask, False)
            aligned_images[im] -= st[0]
            # for normal data: CTF correction
            if options.CTF:
                aligned_images[im] = sp_filter.filt_ctf(
                    aligned_images[im], aligned_images[im].get_attr("ctf"), binary=True
                )
            # for phase plate data we force images to match the rotational power spectrum (provided to us from the outside)
            elif options.VPP:
                aligned_images[im] = sp_fundamentals.fft(
                    sp_filter.filt_table(
                        sp_filter.filt_ctf(
                            sp_fundamentals.fft(aligned_images[im]),
                            aligned_images[im].get_attr("ctf"),
                            binary=True,
                        ),
                        rpw,
                    )
                )

        if options.VPP:
            del rpw

        # NOTE: the defocus value correction above is done with full-sized images, while
        # the normalization below is done on the shrunken / re-scaled images. Because of
        # this we cannot integrate the above into the normalization function below.

        # --------------------------------------------------[ shrinking / re-scaling ]

        # normalize all particle images after applying ctf correction (includes shrinking/re-scaling)
        if Blockdata["myid"] == main_node:
            sp_global_def.sxprint("Create normalized particles")
        normalize_particle_images(
            aligned_images,
            shrink_ratio,
            target_radius,
            target_nx,
            params,
            filament_width=options.filament_width,
            ignore_helical_mask=options.filament_mask_ignore,
        )

        if Blockdata["myid"] == main_node:
            sp_global_def.sxprint("Gather EMData")

        # gather normalized particles at the root node
        sp_utilities.gather_compacted_EMData_to_root(
            Blockdata["total_nima"], aligned_images, myid
        )

        if Blockdata["myid"] == main_node:
            sp_global_def.sxprint("Write aligned stack")
            for i in range(Blockdata["total_nima"]):
                aligned_images[i].write_image(Blockdata["stack_ali2d"], i)
            del aligned_images
            #  It has to be explicitly closed
            DB = db_open_dict(Blockdata["stack_ali2d"])
            DB.close()

            fp = open(
                os.path.join(Blockdata["masterdir"], "README_shrink_ratio.txt"), "w"
            )
            output_text = """
			Since for processing purposes isac changes the image dimensions,
			adjustment of pixel size needs to be made in subsequent steps, (e.g.
			running sxviper.py). The shrink ratio and radius used for this particular isac run is
			--------
			%.5f
			%.5f
			--------
			To get the pixel size for the isac output the user needs to divide
			the original pixel size by the above value. This info is saved in
			the following file: README_shrink_ratio.txt
			""" % (
                shrink_ratio,
                radi,
            )
            fp.write(output_text)
            fp.flush()
            fp.close()
            sp_global_def.sxprint(output_text)
            junk = sp_utilities.cmdexecute(
                "sp_header.py  --consecutive  --params=originalid   %s"
                % Blockdata["stack_ali2d"]
            )

            fp = open(os.path.join(init2dir, "Finished_initial_2d_alignment.txt"), "w")
            fp.flush()
            fp.close()

    else:
        if Blockdata["myid"] == Blockdata["main_node"]:
            sp_global_def.sxprint("Skipping 2D alignment since it was already done!")

    # ------------------------------------------------------[ ISAC main loop ]

    keepdoing_main = True
    main_iter = 0

    while keepdoing_main:

        main_iter += 1
        if checkitem(os.path.join(Blockdata["masterdir"], "finished")):
            keepdoing_main = False

        else:
            if not checkitem(
                os.path.join(Blockdata["masterdir"], "main%03d" % main_iter)
            ):
                #  CREATE masterdir
                #  Create generation000 and put files in it
                generation_iter = 0
                if Blockdata["myid"] == 0:
                    cmd = "{} {}".format(
                        "mkdir",
                        os.path.join(Blockdata["masterdir"], "main%03d" % main_iter),
                    )
                    junk = sp_utilities.cmdexecute(cmd)
                    cmd = "{} {}".format(
                        "mkdir",
                        os.path.join(
                            Blockdata["masterdir"],
                            "main%03d" % main_iter,
                            "generation%03d" % generation_iter,
                        ),
                    )
                    junk = sp_utilities.cmdexecute(cmd)
                    if main_iter > 1:
                        #  It may be restart from unfinished main, so replace files in master
                        cmd = "{} {} {} {}".format(
                            "cp -Rp",
                            os.path.join(
                                Blockdata["masterdir"],
                                "main%03d" % (main_iter - 1),
                                "processed_images.txt",
                            ),
                            os.path.join(
                                Blockdata["masterdir"],
                                "main%03d" % (main_iter - 1),
                                "class_averages.hdf",
                            ),
                            os.path.join(Blockdata["masterdir"]),
                        )
                        junk = sp_utilities.cmdexecute(cmd)
                        junk = os.path.join(
                            Blockdata["masterdir"],
                            "main%03d" % (main_iter - 1),
                            "not_processed_images.txt",
                        )
                        if os.path.exists(junk):
                            cmd = "{} {} {}".format(
                                "cp -Rp", junk, os.path.join(Blockdata["masterdir"])
                            )
                            junk = sp_utilities.cmdexecute(cmd)

                    if os.path.exists(
                        os.path.join(Blockdata["masterdir"], "not_processed_images.txt")
                    ):
                        cmd = "{} {} {}".format(
                            "cp -Rp",
                            os.path.join(
                                Blockdata["masterdir"], "not_processed_images.txt"
                            ),
                            os.path.join(
                                Blockdata["masterdir"],
                                "main%03d" % main_iter,
                                "generation%03d" % generation_iter,
                                "to_process_next_%03d_%03d.txt"
                                % (main_iter, generation_iter),
                            ),
                        )
                        junk = sp_utilities.cmdexecute(cmd)
                    else:
                        sp_utilities.write_text_file(
                            list(range(Blockdata["total_nima"])),
                            os.path.join(
                                Blockdata["masterdir"],
                                "main%03d" % main_iter,
                                "generation%03d" % generation_iter,
                                "to_process_next_%03d_%03d.txt"
                                % (main_iter, generation_iter),
                            ),
                        )
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            if not checkitem(
                os.path.join(Blockdata["masterdir"], "main%03d" % main_iter, "finished")
            ):
                keepdoing_generation = True
                generation_iter = 0

                while keepdoing_generation:
                    generation_iter += 1
                    if checkitem(
                        os.path.join(
                            Blockdata["masterdir"],
                            "main%03d" % main_iter,
                            "generation%03d" % generation_iter,
                        )
                    ):
                        if checkitem(
                            os.path.join(
                                Blockdata["masterdir"],
                                "main%03d" % main_iter,
                                "generation%03d" % generation_iter,
                                "finished",
                            )
                        ):
                            okdo = False
                        else:
                            #  rm -f THIS GENERATION
                            if Blockdata["myid"] == 0:
                                cmd = "{} {}".format(
                                    "rm -rf",
                                    os.path.join(
                                        Blockdata["masterdir"],
                                        "main%03d" % main_iter,
                                        "generation%03d" % generation_iter,
                                    ),
                                )
                                junk = sp_utilities.cmdexecute(cmd)
                            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
                            okdo = True
                    else:
                        okdo = True

                    if okdo:
                        if Blockdata["myid"] == 0:
                            cmd = "{} {}".format(
                                "mkdir",
                                os.path.join(
                                    Blockdata["masterdir"],
                                    "main%03d" % main_iter,
                                    "generation%03d" % generation_iter,
                                ),
                            )
                            junk = sp_utilities.cmdexecute(cmd)
                        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

                        # DO THIS GENERATION
                        keepdoing_main, keepdoing_generation = do_generation(
                            main_iter,
                            generation_iter,
                            target_nx,
                            target_xr,
                            target_yr,
                            target_radius,
                            options,
                        )
                        # Preserve results obtained so far
                        if not keepdoing_generation:
                            if Blockdata["myid"] == 0:
                                cmd = "{} {} {} {}".format(
                                    "cp -Rp",
                                    os.path.join(
                                        Blockdata["masterdir"], "processed_images.txt"
                                    ),
                                    os.path.join(
                                        Blockdata["masterdir"], "class_averages.hdf"
                                    ),
                                    os.path.join(
                                        Blockdata["masterdir"], "main%03d" % main_iter
                                    ),
                                )
                                junk = sp_utilities.cmdexecute(cmd)
                                junk = os.path.join(
                                    Blockdata["masterdir"], "not_processed_images.txt"
                                )
                                if os.path.exists(junk):
                                    cmd = "{} {} {}".format(
                                        "cp -Rp",
                                        junk,
                                        os.path.join(
                                            Blockdata["masterdir"],
                                            "main%03d" % main_iter,
                                        ),
                                    )
                                    junk = sp_utilities.cmdexecute(cmd)

    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    if Blockdata["myid"] == 0:
        if options.skip_ordering:
            pass
        elif os.path.exists(os.path.join(Blockdata["masterdir"], "class_averages.hdf")):
            cmd = "{} {} {} {} {} {} {} {} {} {}".format(
                "sp_chains.py",
                os.path.join(Blockdata["masterdir"], "class_averages.hdf"),
                os.path.join(Blockdata["masterdir"], "junk.hdf"),
                os.path.join(Blockdata["masterdir"], "ordered_class_averages.hdf"),
                "--circular",
                "--radius=%d" % target_radius,
                "--xr=%d" % (target_xr + 1),
                "--yr=%d" % (target_yr + 1),
                "--align",
                ">/dev/null",
            )
            junk = sp_utilities.cmdexecute(cmd)
            cmd = "{} {}".format(
                "rm -rf", os.path.join(Blockdata["masterdir"], "junk.hdf")
            )
            junk = sp_utilities.cmdexecute(cmd)
        else:
            print("ISAC could not find any stable class averaging, terminating...")
    return

def main():
    sp_global_def.print_timestamp("Start")
    run(sys.argv[1:])
    sp_global_def.print_timestamp("Finish")
    mpi.mpi_finalize()

if __name__ == "__main__":
    main()

