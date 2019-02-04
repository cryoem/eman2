#!/home/schoenf/applications/sphire/v1.1/bin/python

#Command to run the code
# mpirun sxisac2.py 'bdb:/home/adnan/applications/isac2data/isac_dummy_data_64#faces' '/home/adnan/applications/isac2results/out/' --radius=32 --img_per_grp=8 --minimum_grp_size=4 --skip_prealignment
"""
Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
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

#====================================================================[ import ]

# compatibility
from __future__ import print_function

# system
import os
import sys
import time

import ctypes
import numpy as np

import mpi

from EMAN2 import EMNumPy

import sparx as spx

from future import standard_library
standard_library.install_aliases()
from builtins import range

from logger import Logger, BaseLogger_Files

import global_def
from global_def   import EMData, Util, ERROR, Transform

from applications import  ali2d_base, MPI_start_end, within_group_refinement

from isac import *

from   optparse import OptionParser, SUPPRESS_HELP
import configparser
from inspect import currentframe, getframeinfo


from random       import randint, seed, jumpahead

from alignment    import Numrinit, ringwe, search_range
from filter       import filt_tanl
from fundamentals import rot_shift2D, fshift, fft, resample
from pixel_error  import multi_align_stability
from statistics   import ave_series

import subprocess
from logger import Logger, BaseLogger_Files

from statistics   import hist_list

from EMAN2db import db_open_dict

#=================================================================[ mpi setup ]

mpi.mpi_init(0, [])

global Blockdata
Blockdata = {}
Blockdata["nproc"]        = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
Blockdata["myid"]         = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
Blockdata["main_node"]    = 0
Blockdata["shared_comm"]  = mpi.mpi_comm_split_type(mpi.MPI_COMM_WORLD, mpi.MPI_COMM_TYPE_SHARED,  0, mpi.MPI_INFO_NULL)
Blockdata["myid_on_node"] = mpi.mpi_comm_rank(Blockdata["shared_comm"])
Blockdata["no_of_processes_per_group"] = mpi.mpi_comm_size(Blockdata["shared_comm"])

masters_from_groups_vs_everything_else_comm = ( mpi.mpi_comm_split(mpi.MPI_COMM_WORLD, Blockdata["main_node"] == Blockdata["myid_on_node"], Blockdata["myid_on_node"]) )
Blockdata["color"], Blockdata["no_of_groups"], balanced_processor_load_on_nodes = spx.get_colors_and_subsets( Blockdata["main_node"], 
                                                                                                              mpi.MPI_COMM_WORLD, Blockdata["myid"],
                                                                                                              Blockdata["shared_comm"], 
                                                                                                              Blockdata["myid_on_node"], 
                                                                                                              masters_from_groups_vs_everything_else_comm )
global_def.BATCH = True
global_def.MPI   = True

NAME_OF_JSON_STATE_FILE = "my_state.json"
NAME_OF_ORIGINAL_IMAGE_INDEX = "originalid"
NAME_OF_RUN_DIR  = "run"
NAME_OF_MAIN_DIR = "generation_"
DIR_DELIM = os.sep

def create_zero_group():

    # select a subset of mpi ids to be in subdivision
    if( Blockdata["myid_on_node"] == 0 ): 
        submyids = [ Blockdata["myid"] ]
    else:  
        submyids = []

    submyids = spx.wrap_mpi_gatherv( submyids, Blockdata["main_node"], mpi.MPI_COMM_WORLD )
    submyids = spx.wrap_mpi_bcast  ( submyids, Blockdata["main_node"], mpi.MPI_COMM_WORLD )
    
    world_group = mpi.mpi_comm_group( mpi.MPI_COMM_WORLD )
    subgroup    = mpi.mpi_group_incl( world_group, len(submyids), submyids )
    
    Blockdata["subgroup_comm"] = mpi.mpi_comm_create( mpi.MPI_COMM_WORLD, subgroup )
    mpi.mpi_barrier( mpi.MPI_COMM_WORLD )

    Blockdata["subgroup_size"] = -1
    Blockdata["subgroup_myid"] = -1
    
    if( mpi.MPI_COMM_NULL != Blockdata["subgroup_comm"] ):
        Blockdata["subgroup_size"] = mpi.mpi_comm_size( Blockdata["subgroup_comm"] )
        Blockdata["subgroup_myid"] = mpi.mpi_comm_rank( Blockdata["subgroup_comm"] )
    
    mpi.mpi_barrier( mpi.MPI_COMM_WORLD )
    return

def checkitem( item, mpi_comm = -1 ):
    global Blockdata
    if mpi_comm == -1:  
        mpi_comm = mpi.MPI_COMM_WORLD
    if( Blockdata["myid"] == Blockdata["main_node"] ):
        if( os.path.exists(item) ): 
            isthere = True
        else: 
            isthere = False
    else: 
        isthere = False
    isthere = spx.bcast_number_to_all( isthere, source_node=Blockdata["main_node"], mpi_comm=mpi_comm )
    mpi.mpi_barrier( mpi_comm )
    return isthere

#======================================================================[ ISAC ]

def iter_isac_pap(alldata, ir, ou, rs, xr, yr, ts, maxit, CTF, snr, dst, FL, FH, FF, init_iter, main_iter, iter_reali, match_first, max_round, match_second, stab_ali, thld_err, indep_run, thld_grp, img_per_grp, generation, candidatesexist = False, random_seed=None, new = False):
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

    #------------------------------------------------------[ mpi related part of the code]

    number_of_proc = Blockdata["nproc"]
    myid = Blockdata["myid"]
    main_node = Blockdata["main_node"]

    seed(myid)
    rand1 = randint(1,1000111222)
    seed(random_seed)
    rand2 = randint(1,1000111222)
    seed(rand1 + rand2)

    if generation == 0:
        ERROR("Generation should begin from 1, please reset it and restart the program", "iter_isac", 1, myid)
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    ali_params_dir = "ali_params_generation_%d"%generation
    if os.path.exists(ali_params_dir):  
        ERROR('Output directory %s for alignment parameters exists, please either change its name or delete it and restart the program'%ali_params_dir, "iter_isac", 1, myid)
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    if new: 
        alimethod = "SHC"
    else:  
        alimethod = ""

    color = 0                        # Blockdata["Ycolor"]
    key = Blockdata["myid"]          # Blockdata["Ykey"]
    group_comm = mpi.MPI_COMM_WORLD  # Blockdata["Ygroup_comm"]
    group_main_node = 0

    nx     = alldata[0].get_xsize()
    ndata  = len(alldata)
    data   = [None]*ndata
    tdummy = Transform({"type":"2D"})

    for im in range(ndata):
        # This is the absolute ID, the only time we use it is
        # when setting the members of 4-way output. All other times, the id in 'members' is 
        # the relative ID.
        alldata[im].set_attr_dict({"xform.align2d": tdummy, "ID": im})
        data[im] = alldata[im]

    avg_num = 0
    Iter = 1
    K = ndata/img_per_grp

    if myid == main_node:
        print("     We will process:  %d current images divided equally between %d groups"%(ndata, K))
        print("*************************************************************************************")

    #------------------------------------------------------[ generate random averages for each group ]

    if key == group_main_node:
        refi = generate_random_averages(data, K, 9023)
    else:
        refi = [spx.model_blank(nx, nx) for i in range(K)]

    for i in range(K):
        spx.bcast_EMData_to_all(refi[i], key, group_main_node, group_comm)

    # create d[K*ndata] matrix 
    orgsize = K*ndata

    if( Blockdata["myid_on_node"] == 0 ): size = orgsize
    else:  size = 0

    disp_unit = np.dtype("f4").itemsize

    win_sm, base_ptr  = mpi.mpi_win_allocate_shared( size*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
    size = orgsize
    if( Blockdata["myid_on_node"] != 0 ):
        base_ptr, = mpi.mpi_win_shared_query(win_sm, mpi.MPI_PROC_NULL)

    d = np.frombuffer(np.core.multiarray.int_asbuffer(base_ptr, size*disp_unit), dtype = 'f4')
    d = d.reshape(orgsize)
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    #------------------------------------------------------[ Generate inital averages ] 

    refi = isac_MPI_pap(data, refi, d, maskfile=None, ir=ir, ou=ou, rs=rs, xrng=xr, yrng=yr, step=ts, 
            maxit=maxit, isac_iter=init_iter, CTF=CTF, snr=snr, rand_seed=-1, color=color, comm=group_comm, 
            stability=True, stab_ali=stab_ali, iter_reali=iter_reali, thld_err=thld_err,
            FL=FL, FH=FH, FF=FF, dst=dst, method = alimethod)

    mpi.mpi_win_free(win_sm)
    del d
    # print "  AFTER FIRST isac_MPI_pap"

    ## broadcast current_refim to all nodes
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    if myid == main_node:
        all_ali_params = [None]*len(data)
        for i,im in enumerate(data):
            alpha, sx, sy, mirror, scale = spx.get_params2D(im)
            all_ali_params[i] = [alpha, sx, sy, mirror]

        print("****************************************************************************************************")
        print("*         Generation finished                 "+time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())+"                            *")
        print("****************************************************************************************************")
        return refi, all_ali_params

    else:  
        return [], []

def isac_MPI_pap(stack, refim, d, maskfile = None, ir=1, ou=-1, rs=1, xrng=0, yrng=0, step=1, maxit=30, isac_iter=10, CTF=False, snr=1.0, rand_seed=-1, color=0, comm=-1, stability=False, stab_ali=5, iter_reali=1, thld_err=1.732, FL=0.1, FH=0.3, FF=0.2, dst=90.0, method = ""):
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

    #------------------------------------------------------[ initialize mpi ]

    if comm == -1: comm = mpi.MPI_COMM_WORLD    

    number_of_proc = mpi.mpi_comm_size(comm)
    myid = mpi.mpi_comm_rank(comm)
    my_abs_id = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
    main_node = 0

    #------------------------------------------------------[ initialize ]

    first_ring = int(ir)
    last_ring  = int(ou)
    rstep      = int(rs)
    max_iter   = int(isac_iter)

    if type(stack) == type(""):
        print("  SHOULD NOT BE HERE")
        sys.exit()
        alldata = EMData.read_images(stack)
    else:
        alldata = stack

    nx = alldata[0].get_xsize()
    ny = alldata[0].get_ysize()

    nima = len(alldata)  # no of images

    #  setall parameters to be zero on input
    for im in range(nima):
        spx.set_params2D(alldata[im], [0.,0.,0.,0, 1.0])

    image_start, image_end = MPI_start_end(nima, number_of_proc, myid)

    if maskfile:
        if type(maskfile) is bytes:
            mask = get_image(maskfile)
        else: 
            mask = maskfile

    else : 
        mask = spx.model_circle(last_ring, nx, nx)

    if type(refim) == type(""):
        refi = EMData.read_images(refim)
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
    cnx = nx/2+1
    cny = cnx

    mode = "F"

    #------------------------------------------------------[ precalculate rings]

    numr = Numrinit(first_ring, last_ring, rstep, mode)    # calculates the necessary information for the 2D polar interpolation (finds the number of rings)
    wr = ringwe(numr, mode)   # Calculate ring weights for rotational alignment

    if rand_seed > -1:
        seed(rand_seed)
    else:
        seed(randint(1,2000111222))
    if myid != main_node:   
        jumpahead(17*myid + 12345)

    previous_agreement = 0.0
    previous_members = [None]*numref
    for j in range(numref):
        previous_members[j] = set()

    fl = FL
    Iter = -1
    main_iter = 0
    terminate = 0
    while( (main_iter < max_iter) and (terminate == 0) ):
        
        Iter += 1
        if my_abs_id == main_node: 
            print("Iteration within isac_MPI Iter =>", Iter, "   main_iter = ", main_iter, " len data = ", image_end-image_start,"   ",time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
        mashi = cnx-ou-2
        
        for j in range(numref):
            refi[j].process_inplace("normalize.mask", {"mask":mask, "no_sigma":1}) # normalize reference images to N(0,1)
            cimage = Util.Polar2Dm(refi[j] , cnx, cny, numr, mode)    # converting reference images to polar cordinates 
            Util.Frngs(cimage, numr)    # Applying FFT on the reference images 
            Util.Applyws(cimage, numr, wr)     # apply weights to FTs of rings
            refi[j] = cimage.copy()

        # inializing the peak list (stores the results of the cross correlations done in 2D alignment)
        peak_list = [np.zeros(4*(image_end-image_start), dtype=np.float32) for i in range(numref)]
        #  nima is the total number of images, not the one on this node, the latter is (image_end-image_start)
        #    d matrix required by EQ-Kmeans can be huge!!  PAP 01/17/2015
        #d = zeros(numref*nima, dtype=float32)
        if( Blockdata["myid_on_node"] == 0 ): 
            d.fill(0.0)

        #--------------------------------------------------[ compute the 2D alignment  ]

        # loop through all images and compute alignment parameters for all references
        for im in range(image_start, image_end):
            alpha, sx, sy, mirror, scale = spx.get_params2D(alldata[im])  # retrieves 2D alignment parameters
            ##  TEST WHETHER PARAMETERS ARE WITHIN RANGE
            alphai, sxi, syi, scalei = spx.inverse_transform2(alpha, sx, sy)  # returns the inverse  geometric transformof the 2d rot and trans matrix
            # If shifts are outside of the permissible range, reset them
            if(abs(sxi)>mashi or abs(syi)>mashi):
                sxi = 0.0
                syi = 0.0
                spx.set_params2D(alldata[im],[0.0,0.0,0.0,0,1.0])  # set 2D alignment parameters (img, alpha, tx, ty, mirror, scale)
            # normalize
            alldata[im].process_inplace("normalize.mask", {"mask":mask, "no_sigma":0}) # subtract average under the mask
            txrng = search_range(nx, ou, sxi, xrng, "ISAC2")   # Find permissible ranges for translational searches by resampling into polar coordinates
                                                               # nx-> image size; ou -> particle radius, sxi-> particle shift, xrng-> max search rng 
            txrng = [txrng[1],txrng[0]]   # for positive and negative shifts
            tyrng = search_range(ny, ou, syi, yrng, "ISAC2")
            tyrng = [tyrng[1],tyrng[0]]

            #----------------------------------------------[ align current image to references]

            # compute alignment parameters for single image with index <im> and all references <refi>
            # NOTE: multiref_polar_ali_2d_peaklist includes conversion of image data into polar coordinates AND an FFT
            temp = Util.multiref_polar_ali_2d_peaklist(alldata[im], refi, txrng, tyrng, step, mode, numr, cnx+sxi, cny+syi)
            for iref in range(numref):
                [alphan, sxn, syn, mn] = \
                   spx.combine_params2(0.0, -sxi, -syi, 0, temp[iref*5+1], temp[iref*5+2], temp[iref*5+3], int(temp[iref*5+4]))  #  Combine 2D alignment parameters including mirror
                peak_list[iref][(im-image_start)*4+0] = alphan
                peak_list[iref][(im-image_start)*4+1] = sxn
                peak_list[iref][(im-image_start)*4+2] = syn
                peak_list[iref][(im-image_start)*4+3] = mn
                d[iref*nima+im] = temp[iref*5] # grab a chunk of the peak list for parallel processing in the next step
            del temp

        del refi

        #--------------------------------------------------[ reduce the results of the correlation computation ]

        if(Blockdata["subgroup_myid"] > -1 ):
            # First check number of nodes, if only one, no reduction necessary.
            if(Blockdata["no_of_groups"] > 1):
                # do reduction using numref chunks nima long (it does not matter what is the actual ordering in d)
                at = time()
                for j in range(numref):
                    dbuf = np.zeros(nima, dtype=np.float32)
                    np.copyto(dbuf,d[j*nima:(j+1)*nima])

                    # reduce entries in section <dbuf> of the overall peak_list 
                    dbuf = mpi.mpi_reduce(dbuf, nima, MPI_FLOAT, MPI_SUM, main_node, Blockdata["subgroup_comm"])  #  RETURNS numpy array
                    if( Blockdata["subgroup_myid"] == 0 ):  
                        np.copyto(d[j*nima:(j+1)*nima],dbuf)

                del dbuf

        #--------------------------------------------------[ find alignment match ]

        if myid == main_node:
            # find maximum in the peak list to determine highest subject/reference match
            id_list_long = Util.assign_groups(str(d.__array_interface__['data'][0]), numref, nima) # string with memory address is passed as parameters
            id_list = [[] for i in range(numref)]
            maxasi = nima/numref  # max. assignmenx

            for i in range(maxasi*numref): id_list[i/maxasi].append(id_list_long[i])
            for i in range(nima%maxasi):   id_list[id_list_long[-1]].append(id_list_long[maxasi*numref+i])
            for iref in range(numref):     id_list[iref].sort()

            del id_list_long

            belongsto = [0]*nima
            for iref in range(numref):
                for im in id_list[iref]:
                    belongsto[im] = iref # image <im> has highest match with reference <iref>
            
            del id_list

        else:
            belongsto = [0]*nima

        # broadcast assignment result to all mpi nodes
        mpi.mpi_barrier(comm)
        belongsto = mpi.mpi_bcast(belongsto, nima, mpi.MPI_INT, main_node, comm)
        belongsto = list(map(int, belongsto))

        #--------------------------------------------------[ update the averages / create new references ]

        members = [0]*numref   # stores number of images assigned to each reference
        sx_sum = [0.0]*numref
        sy_sum = [0.0]*numref
        refi = [ spx.model_blank(nx,ny) for j in range(numref) ]

        for im in range(image_start, image_end):
            
            # get alignment parameters for the chosen matching reference        
            matchref = belongsto[im]
            alphan = float(peak_list[matchref][(im-image_start)*4+0])
            sxn    = float(peak_list[matchref][(im-image_start)*4+1])
            syn    = float(peak_list[matchref][(im-image_start)*4+2])
            mn     =   int(peak_list[matchref][(im-image_start)*4+3])

            if mn == 0: sx_sum[matchref] += sxn
            else:       sx_sum[matchref] -= sxn
            sy_sum[matchref] += syn

            # apply current parameters and add to the average
            Util.add_img( refi[matchref], rot_shift2D(alldata[im], alphan, sxn, syn, mn) )
            members[matchref] += 1

        #--------------------------------------------------[ broadcast  ]

        # everyone computes their reduction and then broadcasts it back to the main_node
        sx_sum  = mpi.mpi_reduce( sx_sum,  numref, mpi.MPI_FLOAT, mpi.MPI_SUM, main_node, comm )
        sy_sum  = mpi.mpi_reduce( sy_sum,  numref, mpi.MPI_FLOAT, mpi.MPI_SUM, main_node, comm )
        members = mpi.mpi_reduce( members, numref, mpi.MPI_INT,   mpi.MPI_SUM, main_node, comm )

        if myid != main_node:
            sx_sum  = [0.0]*numref
            sy_sum  = [0.0]*numref
            members = [0.0]*numref

        sx_sum  = mpi.mpi_bcast( sx_sum,  numref, mpi.MPI_FLOAT, main_node, comm )
        sy_sum  = mpi.mpi_bcast( sy_sum,  numref, mpi.MPI_FLOAT, main_node, comm )
        members = mpi.mpi_bcast( members, numref, mpi.MPI_INT,   main_node, comm )
        sx_sum  = list( map(float, sx_sum) )
        sy_sum  = list( map(float, sy_sum) )
        members = list( map(int, members) )

        for j in range(numref):
            sx_sum[j] /= float(members[j])
            sy_sum[j] /= float(members[j])

        for im in range(image_start, image_end):
            matchref = belongsto[im]
            alphan = float(peak_list[matchref][(im-image_start)*4+0])
            sxn    = float(peak_list[matchref][(im-image_start)*4+1])
            syn    = float(peak_list[matchref][(im-image_start)*4+2])
            mn     =   int(peak_list[matchref][(im-image_start)*4+3])

            if mn == 0:
                spx.set_params2D(alldata[im], [alphan, sxn-sx_sum[matchref], syn-sy_sum[matchref], mn, scale])  # set 2D alignment parameters (img, alpha, tx, ty, mirror, scale)
            else:
                spx.set_params2D(alldata[im], [alphan, sxn+sx_sum[matchref], syn-sy_sum[matchref], mn, scale])

        del peak_list

        #--------------------------------------------------[ ... we work till here (fabian & adnan)]

        for j in range(numref):
            spx.reduce_EMData_to_root(refi[j], myid, main_node, comm)

            if myid == main_node:
                # Golden rule when to do within group refinement
                Util.mul_scalar(refi[j], 1.0/float(members[j]))
                refi[j] = filt_tanl(refi[j], fl, FF)
                refi[j] = fshift(refi[j], -sx_sum[j], -sy_sum[j])
                spx.set_params2D(refi[j], [0.0, 0.0, 0.0, 0, 1.0])

        if myid == main_node:
            # this is most likely meant to center them, if so, it works poorly,
            # it has to be checked and probably a better method used PAP 01/17/2015
            dummy = within_group_refinement(refi, mask, True, first_ring, last_ring, rstep, [xrng], [yrng], [step], dst, maxit, FH, FF)
            ref_ali_params = []

            for j in range(numref):
                alpha, sx, sy, mirror, scale = spx.get_params2D(refi[j])
                refi[j] = rot_shift2D(refi[j], alpha, sx, sy, mirror)
                ref_ali_params.extend([alpha, sx, sy, mirror])
        else:
            ref_ali_params = [0.0]*(numref*4)

        ref_ali_params = mpi.mpi_bcast( ref_ali_params, numref*4, mpi.MPI_FLOAT, main_node, comm )
        ref_ali_params = list( map(float, ref_ali_params) )

        for j in range(numref):
            spx.bcast_EMData_to_all( refi[j], myid, main_node, comm )

        #--------------------------------------------------------------------[ Compensate the centering to averages ]
        for im in range(image_start, image_end):
            matchref = belongsto[im]
            alpha, sx, sy, mirror, scale = spx.get_params2D(alldata[im])
            alphan, sxn, syn, mirrorn = spx.combine_params2(alpha, sx, sy, mirror, ref_ali_params[matchref*4+0], 
                                                                                   ref_ali_params[matchref*4+1],
                                                                                   ref_ali_params[matchref*4+2], 
                                                                                   int(ref_ali_params[matchref*4+3]))
            spx.set_params2D( alldata[im], [alphan, sxn, syn, int(mirrorn), 1.0] )

        fl += 0.05
        ##if my_abs_id == main_node: print "Increased fl .......", fl,localtime()[0:5]
        if fl >= FH:
            fl = FL
            do_within_group = 1
        else:  
            do_within_group = 0

        # Here stability does not need to be checked for each main iteration, it only needs to
        # be done for every 'iter_reali' iterations. If one really wants it to be checked each time
        # simple set iter_reali to 1, which is the default value right now.
        check_stability = (stability and (main_iter%iter_reali==0))

        if do_within_group == 1:

            #--------------------------------------------------------["Doing within group alignment ......."]

            #--------------------------------------------------------[ Broadcast the alignment parameters to all nodes ]
            
            for i in range(number_of_proc):
                im_start, im_end = MPI_start_end(nima, number_of_proc, i)
                if myid == i:
                    ali_params = []
                    for im in range(image_start, image_end):
                        alpha, sx, sy, mirror, scale = spx.get_params2D(alldata[im])
                        ali_params.extend([alpha, sx, sy, mirror])
                else:
                    ali_params = [0.0]*((im_end-im_start)*4)
                ali_params = mpi.mpi_bcast(ali_params, len(ali_params), mpi.MPI_FLOAT, i, comm)
                ali_params = list(map(float, ali_params))
                for im in range(im_start, im_end):
                    alpha = ali_params[(im-im_start)*4]
                    sx = ali_params[(im-im_start)*4+1]
                    sy = ali_params[(im-im_start)*4+2]
                    mirror = int(ali_params[(im-im_start)*4+3])
                    spx.set_params2D(alldata[im], [alpha, sx, sy, mirror, 1.0])

            main_iter += 1

            # There are two approaches to scatter calculations among MPI processes during stability checking.
            # The first one is the original one. I added the second method.
            # Here we try to estimate the calculation time for both approaches.
            stab_calc_time_method_1 = stab_ali * ((numref-1) // number_of_proc + 1)
            stab_calc_time_method_2 = (numref * stab_ali - 1) // number_of_proc + 1
            ##if my_abs_id == main_node: print "Times estimation: ", stab_calc_time_method_1, stab_calc_time_method_2

            # When there is no stability checking or estimated calculation time of new method is greater than 80% of estimated calculation time of original method 
            # then the original method is used. In other case. the second (new) method is used.
            #if (not check_stability) or (stab_calc_time_method_2 > 0.80 * stab_calc_time_method_1):
            #  For the time being only use this method as the other one is not worked out as far as parameter ranges go.
            #if True :
            ##if my_abs_id == main_node: print "Within group refinement and checking within group stability, original approach .......", check_stability, "  ",localtime()[0:5]
            # ====================================== standard approach is used, calculations are parallelized by scatter groups (averages) among MPI processes
            gpixer = []
            for j in range(myid, numref, number_of_proc):
                assign = []
                for im in range(nima):
                    if j == belongsto[im]:
                        assign.append(im)

                randomize = True  # I think there is no reason not to be True
                class_data = [alldata[im] for im in assign]
                refi[j] = within_group_refinement(class_data, mask, randomize, first_ring, last_ring, rstep,
                                                  [xrng], [yrng], [step], dst, maxit, FH, FF, method=method)

                if check_stability:
                    ###if my_abs_id == main_node: print "Checking within group stability, original approach .......", check_stability, "  ",localtime()[0:5]
                    ali_params = [[] for qq in range(stab_ali)]
                    for ii in range(stab_ali):
                        if ii > 0:  # The first one does not have to be repeated
                            dummy = within_group_refinement(class_data, mask, randomize, first_ring, last_ring, rstep,
                                                            [xrng], [yrng], [step], dst, maxit, FH, FF, method=method)
                        for im in range(len(class_data)):
                            alpha, sx, sy, mirror, scale = spx.get_params2D(class_data[im])
                            ali_params[ii].extend([alpha, sx, sy, mirror])

                    stable_set, mirror_consistent_rate, err = multi_align_stability(ali_params, 0.0, 10000.0, thld_err, False, last_ring*2)
                    gpixer.append(err)

                    #print  "Class %4d ...... Size of the group = %4d and of the stable subset = %4d  Mirror consistent rate = %5.3f  Average pixel error prior to class pruning = %10.2f"\
                    #               %(j, len(class_data), len(stable_set), mirror_consistent_rate, err)

                    # If the size of stable subset is too small (say 1, 2), it will cause many problems, so we manually increase it to 5
                    while len(stable_set) < 5:
                        duplicate = True
                        while duplicate:
                            duplicate = False
                            p = randint(0, len(class_data)-1)
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
                        spx.set_params2D( class_data[im], [err[2][0], err[2][1], err[2][2], int(err[2][3]), 1.0] )
                    
                    stable_members.sort()

                    refi[j] = filt_tanl(ave_series(stable_data), FH, FF)
                    refi[j].set_attr('members', stable_members)
                    refi[j].set_attr('n_objects', len(stable_members))
                    #print  "Class %4d ...... Size of the stable subset = %4d  "%(j, len(stable_members))
                    del stable_members
                # end of stability
                del assign

            mpi.mpi_barrier(comm)

            for im in range(nima):
                done_on_node = belongsto[im]%number_of_proc
                if myid == done_on_node:
                    alpha, sx, sy, mirror, scale = spx.get_params2D(alldata[im])
                    ali_params = [alpha, sx, sy, mirror]
                else:
                    ali_params = [0.0]*4

                ali_params = mpi.mpi_bcast(ali_params, 4, mpi.MPI_FLOAT, done_on_node, comm)
                ali_params = list(map(float, ali_params))
                spx.set_params2D(alldata[im], [ali_params[0], ali_params[1], ali_params[2], int(ali_params[3]), 1.0])

            for j in range(numref):
                spx.bcast_EMData_to_all(refi[j], myid, j%number_of_proc, comm)

            terminate = 0
            if check_stability:
                # In this case, we need to set the 'members' attr using stable members from the stability test
                for j in range(numref):
                    #print " SSS ",Blockdata["myid"],j
                    done_on_node = j%number_of_proc
                    if done_on_node != main_node:
                        if myid == main_node:
                            mem_len = mpi.mpi_recv(1, mpi.MPI_INT, done_on_node, SPARX_MPI_TAG_UNIVERSAL, comm)
                            mem_len = int(mem_len[0])
                            members = mpi.mpi_recv(mem_len, mpi.MPI_INT, done_on_node, SPARX_MPI_TAG_UNIVERSAL, comm)
                            members = list(map(int, members))
                            refi[j].set_attr_dict({'members': members,'n_objects': mem_len})

                        elif myid == done_on_node:
                            members = refi[j].get_attr('members')
                            mpi.mpi_send(len(members), 1, mpi.MPI_INT, main_node, SPARX_MPI_TAG_UNIVERSAL, comm)
                            mpi.mpi_send(members, len(members), mpi.MPI_INT, main_node, SPARX_MPI_TAG_UNIVERSAL, comm)

                if myid == main_node:
                    #  compare with previous
                    totprevious = 0.0
                    totcurrent = 0.0
                    common = 0.0
                    for j in range(numref):
                        totprevious += len(previous_members[j])
                        members = set(refi[j].get_attr('members'))
                        totcurrent += len(members)
                        common += len(previous_members[j].intersection(members))
                        previous_members[j] = members
                        #print "  memebers  ",j,len(members)

                    agreement = common/float(totprevious + totcurrent - common)
                    j = agreement - previous_agreement
                    
                    if( (agreement>0.5) and (j > 0.0) and (j < 0.05) ): 
                        terminate = 1
                    
                    previous_agreement = agreement
                    print(">>>  Assignment agreement with previous iteration  %5.1f"%(agreement*100),"   ",time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
                
                terminate = spx.bcast_number_to_all(terminate, source_node = main_node)


            if( check_stability and ( (main_iter == max_iter) or (terminate == 1) ) ):
                #  gather all pixers and print a histogram
                gpixer = spx.wrap_mpi_gatherv(gpixer, main_node, comm)
                if my_abs_id == main_node and color == 0:
                    
                    lhist = 12
                    region, histo = hist_list(gpixer, lhist)
                    print("\n=== Histogram of average within-class pixel errors prior to class pruning ===")
                    for lhx in range(lhist): print("     %10.3f     %7d"%(region[lhx], histo[lhx]))
                    print("=============================================================================\n")

                del gpixer


            ###if myid == main_node:
            ### print  "  WRITING refrealigned  for color:",color
            ### for j in xrange(numref):
            ###     refi[j].write_image("refrealigned%02d_round%02d.hdf"%(color, Iter), j)
            #if myid == main_node:  print "within group alignment done. ", localtime()[0:5]

        # end of do_within_group
        mpi.mpi_barrier(comm)

    if myid == main_node:
        #i = [len(q) for q in id_list]
        #for j in xrange(numref):
        #   refi[j].set_attr_dict({'members': id_list[j], 'n_objects': i[j]})
        #del id_list
        i = [refi[j].get_attr("n_objects") for j in range(numref)]
        lhist = max(12, numref/2)
        region, histo = hist_list(i, lhist)
        print("\n=== Histogram of group sizes ================================================")
        for lhx in range(lhist):  print("     %10.1f     %7d"%(region[lhx], histo[lhx]))
        print("=============================================================================\n")

    mpi.mpi_barrier(comm)

    return refi

def do_generation(main_iter, generation_iter, target_nx, target_xr, target_yr, target_radius, options):
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

    if( Blockdata["myid"] == Blockdata["main_node"] ):
        plist = spx.read_text_file( os.path.join(Blockdata["masterdir"], 
                                                 "main%03d"%main_iter,
                                                 "generation%03d"%(generation_iter-1), 
                                                 "to_process_next_%03d_%03d.txt"%(main_iter,generation_iter-1)) )
        nimastack = len(plist)

    else:
        plist = 0
        nimastack = 0

    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    #print  " willbcast  ",Blockdata["myid"],nimastack,Blockdata["main_node"]
    nimastack = spx.bcast_number_to_all(nimastack, source_node = Blockdata["main_node"], mpi_comm=mpi.MPI_COMM_WORLD)

    #if( Blockdata["myid"] == Blockdata["main_node"] ):  print  "  returned  ",nimastack
    #print  " nimastack2  ",Blockdata["myid"],nimastack

    # Bcast plist to all zero CPUs
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    
    #mpi_finalize()
    #exit()
    
    if(Blockdata["subgroup_myid"] > -1 ):
        # First check number of nodes, if only one, no reduction necessary.
        #print  "  subgroup_myid   ",Blockdata["subgroup_myid"],Blockdata["no_of_groups"],nimastack
        if(Blockdata["no_of_groups"] > 1):          
            plist = bcast_list_to_all(plist, Blockdata["subgroup_myid"], source_node = 0, mpi_comm = Blockdata["subgroup_comm"])
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    #mpi_finalize()
    #exit()

    #delay( Blockdata["myid"],"  AFTER  ")
    
    #mpi_finalize()
    #exit()
    
    # reserve buffers
    disp_unit = np.dtype("f4").itemsize
    size_of_one_image = target_nx*target_nx
    orgsize = nimastack*size_of_one_image #  This is number of projections to be computed simultaneously times their size

    if( Blockdata["myid_on_node"] == 0 ): 
        size = orgsize
    else:  
        size = 0

    win_sm, base_ptr = mpi.mpi_win_allocate_shared( size*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])
    size = orgsize

    #if( Blockdata["myid_on_node"] == 0 ): print " LOHI ",Blockdata["myid"],nimastack,size_of_one_image,size,disp_unit,win_sm,base_ptr

    if( Blockdata["myid_on_node"] != 0 ):
        base_ptr, = mpi.mpi_win_shared_query(win_sm, mpi.MPI_PROC_NULL)

    buffer = np.frombuffer(np.core.multiarray.int_asbuffer(base_ptr, size*disp_unit), dtype = 'f4')
    buffer = buffer.reshape(nimastack, target_nx, target_nx)

    emnumpy2 = spx.EMNumPy()
    bigbuffer = emnumpy2.register_numpy_to_emdata(buffer)

    if( Blockdata["myid_on_node"] == 0 ):
        #  read data on process 0 of each node
        #print "  READING DATA FIRST :",Blockdata["myid"],Blockdata["stack_ali2d"],len(plist)
        for i in range(nimastack):
            bigbuffer.insert_clip(spx.get_im(Blockdata["stack_ali2d"],plist[i]),(0,0,i))
        del plist
        #print "  ALL READ  ",Blockdata["myid"],Blockdata["myid_on_node"],target_nx

    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    #delay(Blockdata["myid"],"  ALL READ DONE")
    alldata = [None]*nimastack
    #img_buffer = [None]*nimastack
    emnumpy3 = [None]*nimastack

    msk = spx.model_blank(target_nx, target_nx,1,1)
    for i in range(nimastack):
        pointer_location = base_ptr + i*size_of_one_image*disp_unit
        img_buffer  = np.frombuffer(np.core.multiarray.int_asbuffer(pointer_location, size_of_one_image*disp_unit), dtype = 'f4')
        img_buffer  = img_buffer.reshape(target_nx, target_nx)
        emnumpy3[i] = spx.EMNumPy()
        alldata[i]  = emnumpy3[i].register_numpy_to_emdata(img_buffer)
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    #  former options
    indep_run    = 0
    match_first  = 0
    thld_grp     = 0
    match_second = 0
    max_round    = 0
    dummy_main_iter = 0

    if( Blockdata["myid"] == 0 ):
        print("*************************************************************************************")
        print("     Main iteration: %3d,  Generation: %3d. "%(main_iter,generation_iter)+"   "+time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))

    ave, all_params = iter_isac_pap(alldata, options.ir, target_radius, options.rs, target_xr, target_yr, options.ts, options.maxit, False, 1.0,\
        options.dst, options.FL, options.FH, options.FF, options.init_iter, dummy_main_iter, options.iter_reali, match_first, \
        max_round, match_second, options.stab_ali, options.thld_err, indep_run, thld_grp, \
        options.img_per_grp, generation_iter, False, random_seed=options.rand_seed, new=False)#options.new)

    #  Clean the stack
    #delay(Blockdata["myid"],"  PROCESSING DONE")
    mpi.mpi_win_free(win_sm)
    emnumpy2.unregister_numpy_from_emdata()
    del emnumpy2
    for i in range(nimastack):  
        emnumpy3[i].unregister_numpy_from_emdata()
    del alldata
    #delay(Blockdata["myid"],"  CLEANED")

    if( Blockdata["myid"] == Blockdata["main_node"] ):
        #  How many averages alreay exist
        if( os.path.exists(os.path.join(Blockdata["masterdir"],"class_averages.hdf")) ):
            nave_exist = global_def.EMUtil.get_image_count(os.path.join(Blockdata["masterdir"],"class_averages.hdf"))
        else: nave_exist = 0
        #  Read all parameters table from masterdir
        all_parameters = spx.read_text_row( os.path.join(Blockdata["masterdir"],"all_parameters.txt"))
        plist = spx.read_text_file(os.path.join(Blockdata["masterdir"], "main%03d"%main_iter, \
            "generation%03d"%(generation_iter-1), "to_process_next_%03d_%03d.txt"%(main_iter,generation_iter-1)) )
        #print "****************************************************************************************************",os.path.join(Blockdata["masterdir"], "main%03d"%main_iter, \
        #   "generation%03d"%(generation_iter-1), "to_process_next_%03d_%03d.txt"%(main_iter,generation_iter-1))
        j = 0
        good = []
        bad = []
        for i,q in enumerate(ave):
            #  Convert local numbering to absolute numbering of images
            local_members = q.get_attr("members")
            members = [plist[l] for l in local_members]
            q.write_image(os.path.join(Blockdata["masterdir"], "main%03d"%main_iter, "generation%03d"%generation_iter,"original_class_averages_%03d_%03d.hdf"%(main_iter,generation_iter)),i)
            q.set_attr("members",members)
            q.write_image(os.path.join(Blockdata["masterdir"], "main%03d"%main_iter, "generation%03d"%generation_iter,"class_averages_%03d_%03d.hdf"%(main_iter,generation_iter)),i)
            if(len(members)> options.minimum_grp_size):
                good += members
                q.write_image(os.path.join(Blockdata["masterdir"], "main%03d"%main_iter, "generation%03d"%generation_iter,"good_class_averages_%03d_%03d.hdf"%(main_iter,generation_iter)),j)
                q.write_image(os.path.join(Blockdata["masterdir"],"class_averages.hdf"),j+nave_exist)
                j += 1

                # We have to update all parameters table
                for l,m in enumerate(members):
                    #  I had to remove it as in case of restart there will be conflicts
                    #if( all_parameters[m][-1] > -1):
                    #   print "  CONFLICT !!!"
                    #   exit()
                    all_parameters[m] = all_params[local_members[l]]
            else:
                bad += members

        if(len(good)> 0):
            spx.write_text_row( all_parameters, os.path.join(Blockdata["masterdir"],"all_parameters.txt"))
            good.sort()
            #  Add currently assigned images to the overall list
            if( os.path.exists( os.path.join(Blockdata["masterdir"], "processed_images.txt") ) ):
                lprocessed = good + spx.read_text_file(os.path.join(Blockdata["masterdir"], "processed_images.txt" ))
                lprocessed.sort()
                spx.write_text_file(lprocessed, os.path.join(Blockdata["masterdir"], "processed_images.txt" ))
            else:
                spx.write_text_file(good, os.path.join(Blockdata["masterdir"], "processed_images.txt" ))
        
        if(len(bad)> 0):
            bad.sort()
            spx.write_text_file(bad, os.path.join(Blockdata["masterdir"], "main%03d"%main_iter, \
                "generation%03d"%(generation_iter), "to_process_next_%03d_%03d.txt"%(main_iter,generation_iter)) )
            
            if( (int(len(bad)*1.2) < 2*options.img_per_grp) or ( (len(good) == 0) and (generation_iter == 1) ) ):
                #  Insufficient number of images to keep processing bad set
                #    or 
                #  Program cannot produce any good averages from what is left at the beginning of new main
                try:  lprocessed = spx.read_text_file(os.path.join(Blockdata["masterdir"], "processed_images.txt" ))
                except:  lprocessed = []
                nprocessed = len(lprocessed)
                leftout = sorted(list(set(range(Blockdata["total_nima"])) - set(lprocessed)))
                spx.write_text_file(leftout, os.path.join(Blockdata["masterdir"], "not_processed_images.txt" ))
                # Check whether what remains can be still processed in a new main interation
                if( ( len(leftout) < 2*options.img_per_grp) or ( (len(good) == 0) and (generation_iter == 1) ) ):
                    #    if the the number of remaining all bad too low full stop
                    keepdoing_main = False
                    keepdoing_generation = False
                    cmd = "{} {}".format("touch", os.path.join(Blockdata["masterdir"], "main%03d"%main_iter, "generation%03d"%generation_iter, "finished"))
                    junk = spx.cmdexecute(cmd)
                    cmd = "{} {}".format("touch", os.path.join(Blockdata["masterdir"], "main%03d"%main_iter, "finished"))
                    junk = spx.cmdexecute(cmd)
                    cmd = "{} {}".format("touch", os.path.join(Blockdata["masterdir"], "finished"))
                    junk = spx.cmdexecute(cmd)
                    print("*         There are no more images to form averages, program finishes     "+time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())+"     *")
                else:
                    #  Will have to increase main, which means putting all bad left as new good, 
                    keepdoing_main = True
                    keepdoing_generation = False
                    #  Will have to increase main, which means putting all bad left as new good, 
                    cmd = "{} {}".format("touch", os.path.join(Blockdata["masterdir"], "main%03d"%main_iter, "generation%03d"%generation_iter, "finished"))
                    junk = spx.cmdexecute(cmd)
                    cmd = "{} {}".format("touch", os.path.join(Blockdata["masterdir"], "main%03d"%main_iter, "finished"))
                    junk = spx.cmdexecute(cmd)
            else:
                keepdoing_main = True
                keepdoing_generation = True
        else:
            keepdoing_main = False
            keepdoing_generation = False
        #print "****************************************************************************************************",keepdoing_main,keepdoing_generation

    else:
        keepdoing_main = False
        keepdoing_generation = False
        
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    keepdoing_main = spx.bcast_number_to_all(keepdoing_main, source_node = Blockdata["main_node"], mpi_comm = mpi.MPI_COMM_WORLD)
    keepdoing_generation = spx.bcast_number_to_all(keepdoing_generation, source_node = Blockdata["main_node"], mpi_comm = mpi.MPI_COMM_WORLD)

    return keepdoing_main, keepdoing_generation

#================================================================[ parameters ]

def parse_parameters( prog_name, usage, args ):

    prog_name = os.path.basename( prog_name )
    parser = OptionParser( usage, version=spx.SPARXVERSION )

    # ISAC command line parameters (public)
    parser.add_option( "--radius",           type="int",                         help="particle radius: there is no default, a sensible number has to be provided, units - pixels (default required int)" )
    parser.add_option( "--target_radius",    type="int",          default=29,    help="target particle radius: actual particle radius on which isac will process data. Images will be shrinked/enlarged to achieve this radius (default 29)" )
    parser.add_option( "--target_nx",        type="int",          default=76,    help="target particle image size: actual image size on which isac will process data. Images will be shrinked/enlarged according to target particle radius and then cut/padded to achieve target_nx size. When xr > 0, the final image size for isac processing is 'target_nx + xr - 1'  (default 76)" )
    parser.add_option( "--img_per_grp",      type="int",          default=200,   help="number of images per class (maximum group size, also defines number of classes K=(total number of images)/img_per_grp (default 200)" )
    parser.add_option( "--minimum_grp_size", type="int",          default=60,    help="minimum size of class (default 60)" )
    parser.add_option( "--CTF",              action="store_true", default=False, help="apply phase-flip for CTF correction: if set the data will be phase-flipped using CTF information included in image headers (default False)" )
    parser.add_option( "--VPP",              action="store_true", default=False, help="Phase Plate data (default False)" )
    parser.add_option( "--ir",               type="int",          default=1,     help="inner ring: of the resampling to polar coordinates. units - pixels (default 1)" )
    parser.add_option( "--rs",               type="int",          default=1,     help="ring step: of the resampling to polar coordinates. units - pixels (default 1)" )
    parser.add_option( "--xr",               type="int",          default=1,     help="x range: of translational search. By default, set by the program. (default 1)" )
    parser.add_option( "--yr",               type="int",          default=-1,    help="y range: of translational search. By default, same as xr. (default -1)" )
    parser.add_option( "--ts",               type="float",        default=1.0,   help="search step: of translational search: units - pixels (default 1.0)" )
    parser.add_option( "--maxit",            type="int",          default=30,    help="number of iterations for reference-free alignment (default 30)" )
    parser.add_option( "--center_method",    type="int",          default=-1,    help="method for centering: of global 2D average during initial prealignment of data (0 : no centering; -1 : average shift method; please see center_2D in utilities.py for methods 1-7) (default -1)" )
    parser.add_option( "--dst",              type="float",        default=90.0,  help="discrete angle used in within group alignment (default 90.0)" )
    parser.add_option( "--FL",               type="float",        default=0.2,   help="lowest stopband: frequency used in the tangent filter (default 0.2)" )
    parser.add_option( "--FH",               type="float",        default=0.45,  help="highest stopband: frequency used in the tangent filter (default 0.45)" )
    parser.add_option( "--FF",               type="float",        default=0.2,   help="fall-off of the tangent filter (default 0.2)" )
    parser.add_option( "--init_iter",        type="int",          default=7,     help="Maximum number of Generation iterations performed for a given subset (default 7)" )
    parser.add_option( "--iter_reali",       type="int",          default=1,     help="SAC stability check interval: every iter_reali iterations of SAC stability checking is performed (default 1)" )
    parser.add_option( "--stab_ali",         type="int",          default=5,     help="number of alignments when checking stability (default 5)" )
    parser.add_option( "--thld_err",         type="float",        default=0.7,   help="threshold of pixel error when checking stability: equals root mean square of distances between corresponding pixels from set of found transformations and theirs average transformation, depends linearly on square of radius (parameter target_radius). units - pixels. (default 0.7)" )
    parser.add_option( "--restart",          type="int",          default='-1',  help="0: restart ISAC2 after last completed main iteration (meaning there is file >finished< in it.  k: restart ISAC2 after k'th main iteration (It has to be completed, meaning there is file >finished< in it. Higer iterations will be removed.)  Default: no restart" )
    parser.add_option( "--rand_seed",        type="int",                         help="random seed set before calculations: useful for testing purposes. By default, total randomness (type int)" )

    # developer parameters (not listed in docs)
    parser.add_option("--skip_prealignment", action="store_true", default=False, help="skip pre-alignment step: to be used if images are already centered. 2dalignment directory will still be generated but the parameters will be zero. (default False)")

    return parser.parse_args(args)

#======================================================================[ main ]

def compile_gpu_code( cuda_path ):
    print( "Compiling GPU code.." )
    subprocess.Popen( ("nvcc "+cuda_path+"multiref_ccf.cu -o "+cuda_path+"multiref.so -shared -Xcompiler -fPIC -lcufft -lcublas -lopencv_core -lopencv_imgcodecs").split() ).wait()  

def get_ptr_array( emdata_list ):
    float_ptr = ctypes.POINTER(ctypes.c_float)

    ptr_list = []
    for img in emdata_list:
        img_np = EMNumPy.em2numpy( img )
        assert img_np.flags['C_CONTIGUOUS'] == True
        assert img_np.dtype == np.float32
        img_ptr = img_np.ctypes.data_as(float_ptr)
        ptr_list.append(img_ptr)
    
    return (float_ptr*len(emdata_list))(*ptr_list)


def run_gpu_testing():

    #---------------------------------------------[ param ]

    DEFAULT_IMG_DIM       = 76
    DEFAULT_RING_NUM      = 32
    DEFAULT_RING_LEN      = 256

    DEFAULT_SHIFT_STEP    = 1
    DEFAULT_SHIFT_RNG_X   = 1
    DEFAULT_SHIFT_RNG_Y   = 1

    CUDA_DEVICE_ID        = 0
    CUDA_STREAM_POOL_SIZE = 8

    CUDA_PATH = "/home/schoenf/work/code/cuISAC/cuda/"

    class CcfBatchConfiguration( ctypes.Structure ):
        _fields_ = [ # data param
                     ("sbj_num", ctypes.c_uint),
                     ("ref_num", ctypes.c_uint),
                     ("img_dim", ctypes.c_uint),
                     # polar sampling parameters
                     ("ring_num", ctypes.c_uint),
                     ("ring_len", ctypes.c_uint),
                     # shift parameters
                     ("shift_step", ctypes.c_uint),
                     ("shift_rng_x", ctypes.c_uint),
                     ("shift_rng_y", ctypes.c_uint),
                     # cuda parameters
                     ("stream_pool_size", ctypes.c_uint) ]

    sbj_num = 4096
    ref_num =    8

    cfg = CcfBatchConfiguration( # data param
                                 sbj_num, 
                                 ref_num,
                                 DEFAULT_IMG_DIM,
                                 # polar sampling parameters
                                 DEFAULT_RING_NUM,
                                 DEFAULT_RING_LEN,
                                 # shift parameters
                                 DEFAULT_SHIFT_STEP,
                                 DEFAULT_SHIFT_RNG_X,
                                 DEFAULT_SHIFT_RNG_Y,
                                 # cuda parameters
                                 CUDA_STREAM_POOL_SIZE )

    #----------------------------------------------[ data ]

    random_data = False

    if random_data:
        ref_stack = np.random.random( [DEFAULT_IMG_DIM * ref_num, DEFAULT_IMG_DIM] ).astype( np.float32 )
        sbj_stack = np.random.random( [DEFAULT_IMG_DIM * ref_num, DEFAULT_IMG_DIM] ).astype( np.float32 )

    else:
        ref_stack = EMData.read_images( "/home/schoenf/work/code/cuISAC/cuda/aln_calibration_test/ref_stack_spx_8.hdf" )
        sbj_stack = EMData.read_images( "/home/schoenf/work/code/cuISAC/cuda/aln_calibration_test/sbj_stack_spx_4096.hdf" )
            
    #---------------------------------------[ c interface ]

    # compile gpu code (optional)
    do_compile = False
    if do_compile:
        compile_gpu_code( CUDA_PATH )

    # c interface setup
    cu_module = ctypes.CDLL( CUDA_PATH + "multiref.so" )
    
    # call c code like any other module function
    cu_module.multiref_alignment_2D( get_ptr_array(sbj_stack), get_ptr_array(ref_stack), ctypes.byref(cfg) )

def main(args):

    #------------------------------------------------------[ gpu code integration testing area ]

    if "--gpu_test" in sys.argv:

        if Blockdata["myid"] == Blockdata["main_node"]:
            print( "Running GPU testing.." )
            run_gpu_testing()
            print( "All done." )

        mpi.mpi_barrier( mpi.MPI_COMM_WORLD )
        mpi.mpi_finalize()
        sys.exit()

    #------------------------------------------------------[ command line parameters ]

    usage = ( sys.argv[0] + " stack_file  [output_directory] --radius=particle_radius" 
              + " --img_per_grp=img_per_grp --CTF <The remaining parameters are" 
              + " optional --ir=ir --rs=rs --xr=xr --yr=yr --ts=ts --maxit=maxit"
              + " --dst=dst --FL=FL --FH=FH --FF=FF --init_iter=init_iter" 
              + " --iter_reali=iter_reali --stab_ali=stab_ali --thld_err=thld_err"
              + " --rand_seed=rand_seed>" )

    options, args = parse_parameters( sys.argv[0], usage, args )

    # only remaining args should be path to input & output folders
    if len(args) > 2:
        print("usage: " + usage)
        print("Please run '" +  sys.argv[0] + " -h' for detailed options")
        sys.exit()
    elif( len(args) == 2):
        Blockdata["stack"]  = args[0]
        Blockdata["masterdir"] = args[1]
    elif( len(args) == 1):
        Blockdata["stack"]  = args[0]
        Blockdata["masterdir"] = ""

    options.new = False
    
    # check required options
    required_option_list = ['radius']

    for required_option in required_option_list:
        if not options.__dict__[required_option]:
            print( "\n ==%s== mandatory option is missing.\n" % required_option )
            print( "Please run '" +  sys.argv[0] + " -h' for detailed options" )
            return 1
    
    # other parameter checks
    if options.minimum_grp_size > options.img_per_grp:
        if Blockdata["myid"] == Blockdata["main_node"]:
            print( "\nERROR! Minimum group size (" + str(options.minimum_grp_size) + ") is larger than the actual group size (" + str(options.img_per_grp) + "). Oh dear :(\n" )
        return 1

    #------------------------------------------------------[ :: ]

    if global_def.CACHE_DISABLE:
        spx.disable_bdb_cache()
    global_def.BATCH = True
    
    main_node = Blockdata["main_node"]
    myid  = Blockdata["myid"]
    nproc = Blockdata["nproc"]

    #  MASTER DIRECTORY
    if(Blockdata["myid"] == Blockdata["main_node"]):
        if( Blockdata["masterdir"] == ""):
            timestring = time.strftime("_%d_%b_%Y_%H_%M_%S", time.localtime())
            Blockdata["masterdir"] = "isac_directory"+timestring
            li = len(Blockdata["masterdir"])
            cmd = "{} {}".format("mkdir", Blockdata["masterdir"])
            junk = spx.cmdexecute(cmd)
        else:
            if not os.path.exists(Blockdata["masterdir"]):
                cmd = "{} {}".format("mkdir", Blockdata["masterdir"])
                junk = spx.cmdexecute(cmd)
            li = 0
    else:
        li = 0
    li = mpi.mpi_bcast(li,1,mpi.MPI_INT,Blockdata["main_node"],mpi.MPI_COMM_WORLD)[0]

    if( li > 0 ):
        Blockdata["masterdir"] = mpi.mpi_bcast(Blockdata["masterdir"],li,MPI_CHAR,Blockdata["main_node"],mpi.MPI_COMM_WORLD)
        Blockdata["masterdir"] = string.join(Blockdata["masterdir"],"")

    Blockdata["stack_ali2d"] = "bdb:" + os.path.join(Blockdata["masterdir"], "stack_ali2d" )
    
    if(myid == main_node):
        print( "****************************************************************************************************" )
        print( "*                                                                                                  *" )
        print( "* ISAC (Iterative Stable Alignment and Clustering)   "+time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())+"                     *" )
        print( "* By Zhengfan Yang, Jia Fang, Francisco Asturias and Pawel A. Penczek                              *" )
        print( "*                                                                                                  *" )
        print( "* REFERENCE: Z. Yang, J. Fang, J. Chittuluru, F. J. Asturias and P. A. Penczek, \"Iterative Stable  *" )
        print( "*            Alignment and Clustering of 2D Transmission Electron Microscope Images\",              *" ) 
        print( "*            Structure 20, 237-247, February 8, 2012.                                              *" )
        print( "*                                                                                                  *" )
        print( "* Last updated: 05/30/2017 PAP                                                                     *" )
        print( "****************************************************************************************************" )
        Util.version()
        print( "****************************************************************************************************" )
        sys.stdout.flush()

        i = "  "
        for a in sys.argv:
            i +=a+"  "
        print(" * shell line command: " )
        print( i )
        print( "*                                                                                                  *" )
        print( "* Master directory: %s"%Blockdata["masterdir"])
        print( "****************************************************************************************************" )
    mpi.mpi_barrier( mpi.MPI_COMM_WORLD )

    create_zero_group()

    if options.CTF and options.VPP: ERROR("Options CTF and VPP cannot be used together", "isac2", 1, myid)

    #  former options
    indep_run    = 0
    match_first  = 0
    match_second = 0
    max_round    = 0
    main_iter    = 0

    radi          = options.radius
    target_radius = options.target_radius
    target_nx     = options.target_nx
    center_method = options.center_method
    if( radi < 1 ):
        ERROR( "Particle radius has to be provided!", "sxisac", 1, myid )

    target_xr  = options.xr
    target_nx += target_xr - 1 # subtract one, which is default
    
    if (options.yr == -1): 
        target_yr = options.xr
    else: 
        target_yr = options.yr

    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    #--------------------------------------[ Initialization of stacks ]
    if(myid == main_node): Blockdata["total_nima"] = global_def.EMUtil.get_image_count(Blockdata["stack"])
    else: Blockdata["total_nima"] = 0

    Blockdata["total_nima"] = spx.bcast_number_to_all(Blockdata["total_nima"], source_node = main_node)

    nxrsteps = 4

    #------------------------------------------------------[ initial 2D alignment (centering) ]

    init2dir = os.path.join(Blockdata["masterdir"],"2dalignment")

    if not checkitem(os.path.join(init2dir, "Finished_initial_2d_alignment.txt")):

        if(myid == 0):

            #  Create output directory
            log2d = Logger(BaseLogger_Files())
            log2d.prefix = os.path.join(init2dir)
            cmd = "mkdir -p "+log2d.prefix
            outcome = subprocess.call(cmd, shell=True)
            log2d.prefix += "/"
        else:
            outcome = 0
            log2d = None

        if(myid == main_node):
            a = spx.get_im(Blockdata["stack"])
            nnxo = a.get_xsize()
        else:
            nnxo = 0
        nnxo = spx.bcast_number_to_all(nnxo, source_node = main_node)

        image_start, image_end = spx.MPI_start_end(Blockdata["total_nima"], nproc, myid)

        original_images = EMData.read_images(Blockdata["stack"], list(range(image_start,image_end)))

        if options.VPP:
            ntp = len(rops_table(original_images[0]))
            rpw = [0.0]*ntp
            for q in original_images:
                tpw = rops_table(q)
                for i in range(ntp):  rpw[i] += sqrt(tpw[i])
            del tpw
            rpw = mpi.mpi_reduce(rpw, ntp, MPI_FLOAT, MPI_SUM, main_node, mpi.MPI_COMM_WORLD)
            if(myid == 0):
                rpw = [float(Blockdata["total_nima"]/q) for q in rpw]
                rpw[0] = 1.0
                spx.write_text_file(rpw,os.path.join(Blockdata["masterdir"], "rpw.txt"))
            else:  rpw = []
            rpw = bcast_list_to_all(rpw, myid, source_node = main_node, mpi_comm = mpi.MPI_COMM_WORLD)
            for i in range(len(original_images)):  original_images[i] = filt_table(original_images[i], rpw)
            del rpw
        else:
            if(myid == 0):
                ntp = len( spx.rops_table(original_images[0]) )
                spx.write_text_file([0.0]*ntp,os.path.join(Blockdata["masterdir"], "rpw.txt"))

        if options.skip_prealignment:
            params2d = [[0.0,0.0,0.0,0] for i in range(image_start, image_end)]
        else:
            #  We assume the target radius will be 29, and xr = 1.
            shrink_ratio = float(target_radius)/float(radi)

            for im in range(len(original_images)):
                if(shrink_ratio != 1.0):
                    original_images[im] = resample(original_images[im], shrink_ratio)

            nx = original_images[0].get_xsize()
            # nx = int(nx*shrink_ratio + 0.5)

            txrm = (nx - 2*(target_radius+1))//2
            if(txrm < 0):           ERROR( "ERROR!!   Radius of the structure larger than the window data size permits   %d"%(radi), "sxisac",1, myid)
            if(txrm/nxrsteps>0):
                tss = ""
                txr = ""
                while(txrm/nxrsteps>0):
                    tts=txrm/nxrsteps
                    tss += "  %d"%tts
                    txr += "  %d"%(tts*nxrsteps)
                    txrm =txrm//2
            else:
                tss = "1"
                txr = "%d"%txrm
            
            # 2D alignment
            if(Blockdata["myid"] == 0): print("* 2D alignment   "+time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
            params2d = ali2d_base(original_images, init2dir, None, 1, target_radius, 1, txr, txr, tss, \
                False, 90.0, center_method, 14, options.CTF, 1.0, False, \
                "ref_ali2d", "", log2d, nproc, myid, main_node, mpi.MPI_COMM_WORLD, write_headers = False)

            del original_images

            for i in range(len(params2d)):
                alpha, sx, sy, mirror = spx.combine_params2(0, params2d[i][1],params2d[i][2], 0, -params2d[i][0], 0, 0, 0)
                sx /= shrink_ratio
                sy /= shrink_ratio
                params2d[i][0] = 0.0
                params2d[i][1] = sx
                params2d[i][2] = sy
                params2d[i][3] = 0

        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        tmp = params2d[:]
        tmp = spx.wrap_mpi_gatherv(tmp, main_node, mpi.MPI_COMM_WORLD)
        if( myid == main_node ):        
            if options.skip_prealignment:
                print("=========================================")
                print("There is no alignment step, '%s' params are set to zero for later use."%os.path.join(init2dir, "initial2Dparams.txt"))
                print("=========================================")
            spx.write_text_row(tmp,os.path.join(init2dir, "initial2Dparams.txt"))
        del tmp
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    
        #  We assume the target image size will be target_nx, radius will be 29, and xr = 1.  
        #  Note images can be also padded, in which case shrink_ratio > 1.
        shrink_ratio = float(target_radius)/float(radi)
        
        aligned_images = EMData.read_images(Blockdata["stack"], list(range(image_start,image_end)))
        nx = aligned_images[0].get_xsize()
        nima = len(aligned_images)
        newx = int(nx*shrink_ratio + 0.5)
       
        while not os.path.exists(os.path.join(init2dir, "initial2Dparams.txt")):
            time.sleep(1)
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        
        params = spx.read_text_row(os.path.join(init2dir, "initial2Dparams.txt"))
        params = params[image_start:image_end]

        msk = spx.model_circle(radi, nx, nx)
        if options.VPP:
            if myid == 0:  rpw = spx.read_text_file(os.path.join(Blockdata["masterdir"], "rpw.txt"))
            else:  rpw = [0.0]
            rpw = bcast_list_to_all(rpw, myid, source_node = main_node, mpi_comm = mpi.MPI_COMM_WORLD)
        for im in range(nima):
            st = Util.infomask(aligned_images[im], msk, False)
            aligned_images[im] -= st[0]
            if options.CTF:
                aligned_images[im] = filt_ctf(aligned_images[im], aligned_images[im].get_attr("ctf"), binary = True)
            elif options.VPP:
                aligned_images[im] = fft(filt_table(filt_ctf(fft(aligned_images[im]), aligned_images[im].get_attr("ctf"), binary = True), rpw))
        if options.VPP: del rpw
        if(shrink_ratio < 1.0):
            if    newx > target_nx  :
                msk = spx.model_circle(target_radius, target_nx, target_nx)
                for im in range(nima):
                    #  Here we should use only shifts
                    #alpha, sx, sy, mirror, scale = get_params2D(aligned_images[im])
                    #alpha, sx, sy, mirror = combine_params2(0, sx,sy, 0, -alpha, 0, 0, 0)
                    #aligned_images[im] = rot_shift2D(aligned_images[im], 0, sx, sy, 0)
                    aligned_images[im] = rot_shift2D(aligned_images[im], 0, params[im][1], params[im][2], 0)
                    aligned_images[im] = resample(aligned_images[im], shrink_ratio)
                    aligned_images[im] = Util.window(aligned_images[im], target_nx, target_nx, 1)
                    p = Util.infomask(aligned_images[im], msk, False)
                    aligned_images[im] -= p[0]
                    p = Util.infomask(aligned_images[im], msk, True)
                    aligned_images[im] /= p[1]
            elif  newx == target_nx :
                msk = spx.model_circle(target_radius, target_nx, target_nx)
                for im in range(nima):
                    #  Here we should use only shifts
                    #alpha, sx, sy, mirror, scale = get_params2D(aligned_images[im])
                    #alpha, sx, sy, mirror = combine_params2(0, sx,sy, 0, -alpha, 0, 0, 0)
                    aligned_images[im] = rot_shift2D(aligned_images[im], 0, params[im][1], params[im][2], 0)
                    aligned_images[im] = resample(aligned_images[im], shrink_ratio)
                    p = Util.infomask(aligned_images[im], msk, False)
                    aligned_images[im] -= p[0]
                    p = Util.infomask(aligned_images[im], msk, True)
                    aligned_images[im] /= p[1]
            elif  newx < target_nx  :   
                msk = spx.model_circle(newx//2-2, newx,  newx)
                for im in range(nima):
                    #  Here we should use only shifts
                    #alpha, sx, sy, mirror, scale = get_params2D(aligned_images[im])
                    #alpha, sx, sy, mirror = combine_params2(0, sx,sy, 0, -alpha, 0, 0, 0)
                    aligned_images[im] = rot_shift2D(aligned_images[im], 0, params[im][1], params[im][2], 0)
                    aligned_images[im] = resample(aligned_images[im], shrink_ratio)
                    p = Util.infomask(aligned_images[im], msk, False)
                    aligned_images[im] -= p[0]
                    p = Util.infomask(aligned_images[im], msk, True)
                    aligned_images[im] /= p[1]
                    aligned_images[im] = spx.pad(aligned_images[im], target_nx, target_nx, 1, 0.0)
        elif(shrink_ratio == 1.0):
            if    newx > target_nx  :
                msk = spx.model_circle(target_radius, target_nx, target_nx)
                for im in range(nima):
                    #  Here we should use only shifts
                    #alpha, sx, sy, mirror, scale = get_params2D(aligned_images[im])
                    #alpha, sx, sy, mirror = combine_params2(0, sx,sy, 0, -alpha, 0, 0, 0)
                    aligned_images[im] = rot_shift2D(aligned_images[im], 0, params[im][1], params[im][2], 0)
                    aligned_images[im] = Util.window(aligned_images[im], target_nx, target_nx, 1)
                    p = Util.infomask(aligned_images[im], msk, False)
                    aligned_images[im] -= p[0]
                    p = Util.infomask(aligned_images[im], msk, True)
                    aligned_images[im] /= p[1]
            elif  newx == target_nx :
                msk = spx.model_circle(target_radius, target_nx, target_nx)
                for im in range(nima):
                    #  Here we should use only shifts
                    #alpha, sx, sy, mirror, scale = get_params2D(aligned_images[im])
                    #alpha, sx, sy, mirror = combine_params2(0, sx,sy, 0, -alpha, 0, 0, 0)
                    aligned_images[im] = rot_shift2D(aligned_images[im], 0, params[im][1], params[im][2], 0)
                    p = Util.infomask(aligned_images[im], msk, False)
                    aligned_images[im] -= p[0]
                    p = Util.infomask(aligned_images[im], msk, True)
                    aligned_images[im] /= p[1]
            elif  newx < target_nx  :           
                msk = spx.model_circle(newx//2-2, newx,  newx)
                for im in range(nima):
                    #  Here we should use only shifts
                    #alpha, sx, sy, mirror, scale = get_params2D(aligned_images[im])
                    #alpha, sx, sy, mirror = combine_params2(0, sx,sy, 0, -alpha, 0, 0, 0)
                    aligned_images[im] = rot_shift2D(aligned_images[im], 0, params[im][1], params[im][2], 0)
                    #aligned_images[im]  = resample(aligned_images[im], shrink_ratio)
                    p = Util.infomask(aligned_images[im], msk, False)
                    aligned_images[im] -= p[0]
                    p = Util.infomask(aligned_images[im], msk, True)
                    aligned_images[im] /= p[1]
                    aligned_images[im] = spx.pad(aligned_images[im], target_nx, target_nx, 1, 0.0)
        elif(shrink_ratio > 1.0):
            if    newx > target_nx  :
                msk = spx.model_circle(target_radius, target_nx, target_nx)
                for im in range(nima):
                    #  Here we should use only shifts
                    #alpha, sx, sy, mirror, scale = get_params2D(aligned_images[im])
                    #alpha, sx, sy, mirror = combine_params2(0, sx,sy, 0, -alpha, 0, 0, 0)
                    aligned_images[im] = rot_shift2D(aligned_images[im], 0, params[im][1], params[im][2], 0)
                    aligned_images[im] = resample(aligned_images[im], shrink_ratio)
                    aligned_images[im] = Util.window(aligned_images[im], target_nx, target_nx, 1)
                    p = Util.infomask(aligned_images[im], msk, False)
                    aligned_images[im] -= p[0]
                    p = Util.infomask(aligned_images[im], msk, True)
                    aligned_images[im] /= p[1]
            elif  newx == target_nx :
                msk = spx.model_circle(target_radius, target_nx, target_nx)
                for im in range(nima):
                    #  Here we should use only shifts
                    #alpha, sx, sy, mirror, scale = get_params2D(aligned_images[im])
                    #alpha, sx, sy, mirror = combine_params2(0, sx,sy, 0, -alpha, 0, 0, 0)
                    aligned_images[im] = rot_shift2D(aligned_images[im], 0, params[im][1], params[im][2], 0)
                    aligned_images[im] = resample(aligned_images[im], shrink_ratio)
                    p = Util.infomask(aligned_images[im], msk, False)
                    aligned_images[im] -= p[0]
                    p = Util.infomask(aligned_images[im], msk, True)
                    aligned_images[im] /= p[1]
            elif  newx < target_nx  :
                msk = spx.model_circle(newx//2-2, newx,  newx)
                for im in range(nima):
                    #  Here we should use only shifts
                    #alpha, sx, sy, mirror, scale = get_params2D(aligned_images[im])
                    #alpha, sx, sy, mirror = combine_params2(0, sx,sy, 0, -alpha, 0, 0, 0)
                    aligned_images[im] = rot_shift2D(aligned_images[im], 0, params[im][1], params[im][2], 0)
                    aligned_images[im] = resample(aligned_images[im], shrink_ratio)
                    p = Util.infomask(aligned_images[im], msk, False)
                    aligned_images[im] -= p[0]
                    p = Util.infomask(aligned_images[im], msk, True)
                    aligned_images[im] /= p[1]
                    aligned_images[im] = spx.pad(aligned_images[im], target_nx, target_nx, 1, 0.0)
        del msk
    
        spx.gather_compacted_EMData_to_root(Blockdata["total_nima"], aligned_images, myid)
    
        if( Blockdata["myid"] == main_node ):
            for i in range( Blockdata["total_nima"] ):  
                aligned_images[i].write_image( Blockdata["stack_ali2d"],i )
            del aligned_images
            
            DB = db_open_dict(Blockdata["stack_ali2d"])
            DB.close() # has to be explicitly closed

            fp = open( os.path.join(Blockdata["masterdir"], "README_shrink_ratio.txt"), "w" )
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
            """%(shrink_ratio, radi)
            fp.write(output_text); fp.flush() ;fp.close()
            print(output_text)
            junk = spx.cmdexecute("sxheader.py  --consecutive  --params=originalid   %s"%Blockdata["stack_ali2d"])

            fp = open(os.path.join(init2dir, "Finished_initial_2d_alignment.txt"), "w"); fp.flush() ;fp.close()
    else:
        if( Blockdata["myid"] == Blockdata["main_node"] ):
            print("Skipping 2D alignment since it was already done!")

    #------------------------------------------------------[ pre-alignment finished ]

    error = 0
    if( Blockdata["myid"] == Blockdata["main_node"] ):
        if( not os.path.exists( os.path.join(Blockdata["masterdir"], "main001", "generation000") ) ):
            #  We do not create processed_images.txt selection file as it has to be initially empty
            #  We do create all params filled with zeroes.
            spx.write_text_row( [[0.0,0.0,0.0,-1] for i in range(Blockdata["total_nima"])], os.path.join(Blockdata["masterdir"], "all_parameters.txt") )
            if(options.restart > -1):
                error = 1
        else:
            if(options.restart == 0):
                keepdoing_main = True
                main_iter = 0
                while(keepdoing_main):
                    main_iter += 1
                    if( os.path.exists( os.path.join(Blockdata["masterdir"], "main%03d"%main_iter ) ) ):
                        if( not  os.path.exists( os.path.join(Blockdata["masterdir"], "main%03d"%main_iter, "finished" ) ) ):
                            cmd = "{} {}".format( "rm -rf", 
                                                  os.path.join(Blockdata["masterdir"], "main%03d"%main_iter) )
                            junk = spx.cmdexecute(cmd)
                            keepdoing_main = False
                    else: 
                        keepdoing_main = False

            else:
                main_iter = options.restart
                if( not os.path.exists( os.path.join(Blockdata["masterdir"], "main%03d"%main_iter, "finished")) ):
                    error = 1
                else:
                    keepdoing_main = True
                    main_iter += 1
                    while( keepdoing_main ):
                        if( os.path.exists(os.path.join(Blockdata["masterdir"], "main%03d"%main_iter)) ):
                            cmd = "{} {}".format( "rm -rf", 
                                                  os.path.join(Blockdata["masterdir"], "main%03d"%main_iter) )
                            junk = spx.cmdexecute( cmd )
                            main_iter += 1
                        else: 
                            keepdoing_main = False

            if( os.path.exists(os.path.join(Blockdata["masterdir"], "finished")) ):
                cmd = "{} {}".format( "rm -rf", 
                                      os.path.join(Blockdata["masterdir"], "finished") )
                junk = spx.cmdexecute( cmd )

    error = spx.bcast_number_to_all( error, source_node = Blockdata["main_node"] )
    if( error == 1 ):
        ERROR( "isac2","cannot restart from unfinished main iteration %d" % main_iter )

    mpi.mpi_barrier( mpi.MPI_COMM_WORLD )

    #------------------------------------------------------[ ISAC main loop ]

    keepdoing_main = True
    main_iter = 0

    while( keepdoing_main ):

        main_iter += 1
        if( checkitem(os.path.join(Blockdata["masterdir"], "finished")) ):
            keepdoing_main = False

        else:
            if( not checkitem(os.path.join(Blockdata["masterdir"], "main%03d"%main_iter )) ):
                #  CREATE masterdir
                #  Create generation000 and put files in it
                generation_iter = 0
                if( Blockdata["myid"] == 0 ):
                    cmd = "{} {}".format( "mkdir", 
                                          os.path.join(Blockdata["masterdir"], "main%03d"%main_iter) )
                    junk = spx.cmdexecute(cmd)
                    cmd = "{} {}".format( "mkdir", 
                                          os.path.join(Blockdata["masterdir"], "main%03d"%main_iter, "generation%03d"%generation_iter) )
                    junk = spx.cmdexecute(cmd)
                    if( main_iter > 1 ):
                        #  It may be restart from unfinished main, so replace files in master
                        cmd = "{} {} {} {}".format( "cp -Rp",
                                                    os.path.join(Blockdata["masterdir"], "main%03d"%(main_iter-1), "processed_images.txt"),
                                                    os.path.join(Blockdata["masterdir"], "main%03d"%(main_iter-1), "class_averages.hdf"),
                                                    os.path.join(Blockdata["masterdir"]) )
                        junk = spx.cmdexecute( cmd )
                        junk = os.path.join( Blockdata["masterdir"], "main%03d"%(main_iter-1), "not_processed_images.txt" )
                        if( os.path.exists(junk) ):
                            cmd = "{} {} {}".format( "cp -Rp", 
                                                     junk, 
                                                     os.path.join(Blockdata["masterdir"]) )
                            junk = spx.cmdexecute( cmd )

                    if( os.path.exists( os.path.join(Blockdata["masterdir"], "not_processed_images.txt")) ):
                        cmd = "{} {} {}".format( "cp -Rp", 
                                                 os.path.join(Blockdata["masterdir"], "not_processed_images.txt"),
                                                 os.path.join(Blockdata["masterdir"], 
                                                              "main%03d"%main_iter, 
                                                              "generation%03d"%generation_iter, 
                                                              "to_process_next_%03d_%03d.txt"%(main_iter,generation_iter)) )
                        junk = spx.cmdexecute( cmd )
                    else:
                        spx.write_text_file( list(range(Blockdata["total_nima"])),
                                             os.path.join(Blockdata["masterdir"], 
                                                          "main%03d"%main_iter, 
                                                          "generation%03d"%generation_iter, 
                                                          "to_process_next_%03d_%03d.txt"%(main_iter,generation_iter)) )
                mpi.mpi_barrier( mpi.MPI_COMM_WORLD )

            if( not checkitem(os.path.join(Blockdata["masterdir"], "main%03d"%main_iter, "finished")) ):
                keepdoing_generation = True
                generation_iter = 0
    
                while( keepdoing_generation ):
                    generation_iter += 1
                    if checkitem( os.path.join(Blockdata["masterdir"], "main%03d"%main_iter, "generation%03d"%generation_iter) ):
                        if checkitem( os.path.join(Blockdata["masterdir"], "main%03d"%main_iter, "generation%03d"%generation_iter, "finished") ):
                            okdo = False
                        else:
                            #  rm -f THIS GENERATION
                            if( Blockdata["myid"] == 0 ):
                                cmd = "{} {}".format( "rm -rf", 
                                                      os.path.join(Blockdata["masterdir"], "main%03d"%main_iter, "generation%03d"%generation_iter) )
                                junk = spx.cmdexecute( cmd )
                            mpi.mpi_barrier( mpi.MPI_COMM_WORLD )
                            okdo = True
                    else:
                        okdo = True

                    if okdo:
                        if( Blockdata["myid"] == 0 ):
                            cmd = "{} {}".format( "mkdir", 
                                                  os.path.join(Blockdata["masterdir"], "main%03d"%main_iter, "generation%03d"%generation_iter) )
                            junk = spx.cmdexecute( cmd )
                        mpi.mpi_barrier( mpi.MPI_COMM_WORLD )

                        # DO THIS GENERATION
                        keepdoing_main, keepdoing_generation = do_generation(main_iter, generation_iter, target_nx, target_xr, target_yr, target_radius, options)
                        # Preserve results obtained so far
                        if( not keepdoing_generation ):
                            if( Blockdata["myid"] == 0 ):
                                cmd = "{} {} {} {}".format( "cp -Rp",
                                                            os.path.join(Blockdata["masterdir"], "processed_images.txt"),
                                                            os.path.join(Blockdata["masterdir"], "class_averages.hdf"),
                                                            os.path.join(Blockdata["masterdir"], "main%03d"%main_iter) )
                                junk = spx.cmdexecute( cmd )
                                junk = os.path.join( Blockdata["masterdir"], "not_processed_images.txt" )
                                if( os.path.exists(junk) ):
                                    cmd = "{} {} {}".format( "cp -Rp",
                                                             junk,
                                                             os.path.join(Blockdata["masterdir"], "main%03d"%main_iter) )
                                    junk = spx.cmdexecute( cmd )

    mpi.mpi_barrier( mpi.MPI_COMM_WORLD )
    if( Blockdata["myid"] == 0 ):
        if( os.path.exists(os.path.join(Blockdata["masterdir"],"class_averages.hdf")) ):
            cmd = "{} {} {} {} {} {} {} {} {} {}".format( "sxchains.py", 
                                                          os.path.join(Blockdata["masterdir"],"class_averages.hdf"),
                                                          os.path.join(Blockdata["masterdir"],"junk.hdf"),
                                                          os.path.join(Blockdata["masterdir"],"ordered_class_averages.hdf"),
                                                          "--circular",
                                                          "--radius=%d"%target_radius ,
                                                          "--xr=%d"%(target_xr+1),
                                                          "--yr=%d"%(target_yr+1),
                                                          "--align", 
                                                          ">/dev/null" )
            junk = spx.cmdexecute( cmd )
            cmd = "{} {}".format( "rm -rf", 
                                  os.path.join(Blockdata["masterdir"], "junk.hdf") )
            junk = spx.cmdexecute( cmd )
        else:
            print( "ISAC could not find any stable class averaging, terminating..." )

    mpi.mpi_finalize()
    sys.exit()

if __name__=="__main__":
    main(sys.argv[1:])
