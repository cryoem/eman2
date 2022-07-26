#! /usr/bin/env python
from __future__ import print_function

#
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
#

from past.utils import old_div
import os
import sys
import os
import types
import mpi
import EMAN2
import optparse
import random
import numpy
import EMAN2_cppwrap
from sphire.libpy import sp_global_def
from sphire.libpy import sp_user_functions
from sphire.libpy import sp_utilities
from sphire.libpy import sp_statistics
from sphire.libpy import sp_filter
from sphire.libpy import sp_alignment
from sphire.libpy import sp_fundamentals
from sphire.libpy import sp_morphology
from sphire.libpy import sp_applications

mpi.mpi_init( 0, [] )

def mref_ali2d_MPI(stack, refim, outdir, maskfile = None, ir=1, ou=-1, rs=1, 
                   xrng=0, yrng=0, step=1, center=1, maxit=10, CTF=False, snr=1.0, 
                   user_func_name="ref_ali2d", rand_seed=1000, upscale_isac_dir=None):
    """
    2D multi-reference alignment using rotational ccf in polar coordinates and quadratic interpolation.
    
    Adapted from sp_applications.py
    
    Changelog:
        2021-02-18 (trs) -- Updated for SPHIRE 1.4
        2020-10-23 (trs) -- Added option to upscale reference images based on ISAC SHRINK_RATIO
        
    Parameters:
        stack
        refim
        outdir
        maskfile
        ir
        ou
        rs
        xrng
        yrng
        step
        center
        maxit
        CTF
        snr
        user_func_name
        rand_seed
        upscale_isac_dir
    """

    number_of_proc = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
    myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
    main_node = 0
    
    # create the output directory, if it does not exist
    if os.path.exists(outdir):  
        sp_global_def.ERROR(
            'Output directory exists, please change the name and restart the program', 
            "mref_ali2d_MPI ", 
            1, 
            myid
            )
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    if myid == main_node:
        os.makedirs(outdir)
        sp_global_def.LOGFILE =  os.path.join(outdir, sp_global_def.LOGFILE)
        sp_global_def.sxprint("mref_ali2d_MPI")

    nima = EMAN2.EMUtil.get_image_count(stack)
    
    image_start, image_end = sp_applications.MPI_start_end(nima, number_of_proc, myid)

    nima = EMAN2.EMUtil.get_image_count(stack)
    ima = EMAN2.EMData()
    ima.read_image(stack, image_start)

    first_ring=int(ir); last_ring=int(ou); rstep=int(rs); max_iter=int(maxit)

    if max_iter == 0:
        max_iter  = 10
        auto_stop = True
    else:
        auto_stop = False

    if myid == main_node:
        sp_global_def.sxprint("Input stack                 : %s"%(stack))
        sp_global_def.sxprint("Reference stack             : %s"%(refim))	
        sp_global_def.sxprint("Output directory            : %s"%(outdir))
        sp_global_def.sxprint("Maskfile                    : %s"%(maskfile))
        sp_global_def.sxprint("Inner radius                : %i"%(first_ring))

    nx = ima.get_xsize()
    # default value for the last ring
    if last_ring == -1: last_ring= old_div(nx, 2) - 2
    
    if myid == main_node:
        sp_global_def.sxprint("Outer radius                : %i"%(last_ring))
        sp_global_def.sxprint("Ring step                   : %i"%(rstep))
        sp_global_def.sxprint("X search range              : %f"%(xrng))
        sp_global_def.sxprint("Y search range              : %f"%(yrng))
        sp_global_def.sxprint("Translational step          : %f"%(step))
        sp_global_def.sxprint("Center type                 : %i"%(center))
        sp_global_def.sxprint("Maximum iteration           : %i"%(max_iter))
        sp_global_def.sxprint("CTF correction              : %s"%(CTF))
        sp_global_def.sxprint("Signal-to-Noise Ratio       : %f"%(snr))
        sp_global_def.sxprint("Random seed                 : %i"%(rand_seed))	
        sp_global_def.sxprint("User function               : %s\n"%(user_func_name))
    user_func = sp_user_functions.factory[user_func_name]

    if maskfile:
        if type(maskfile) is str:  
            mask = sp_utilities.get_image(maskfile)
        else: 
            mask = maskfile
    else: 
        mask = sp_utilities.model_circle(last_ring, nx, nx)
    
    #  references, do them on all processors...
    refi = []
    numref = EMAN2.EMUtil.get_image_count(refim)
    
    #  CTF stuff
    if CTF:
        ctf_params = ima.get_attr("ctf")
        data_had_ctf = ima.get_attr("ctf_applied")
        ctm = sp_morphology.ctf_2(nx, ctf_params)
        lctf = len(ctm)

    # IMAGES ARE SQUARES! center is in SPIDER convention
    cnx = old_div(nx, 2) + 1
    cny = cnx

    mode = "F"
    #precalculate rings
    numr = sp_alignment.Numrinit(first_ring, last_ring, rstep, mode)
    wr = sp_alignment.ringwe(numr, mode)
    # reference images
    again = True
    params = []
    
    # prepare reference images on all nodes
    ima.to_zero()
    
    # Upscale references if necessary
    refdim = sp_utilities.get_im(refim, 0).get_xsize()
    if upscale_isac_dir!=None:
        isac_shrink_path= os.path.join(upscale_isac_dir, "README_shrink_ratio.txt")
        isac_shrink_file= open(isac_shrink_path, "r")
        isac_shrink_lines= isac_shrink_file.readlines()
        isac_shrink_ratio= float(isac_shrink_lines[5])
        
        if myid == main_node : 
            sp_global_def.sxprint("Upscale ISAC directory      : %s" % upscale_isac_dir)
        
    for j in range(numref):
        #  even, odd, number of even, number of images.  After frc, totav
        if upscale_isac_dir==None:
            refi.append([sp_utilities.get_im(refim,j), ima.copy(), 0])
        
        # Resize if necessary
        else:
            rescaled = sp_fundamentals.resample(
                sp_utilities.get_im(refim,j), old_div(1., isac_shrink_ratio) 
                )
            
            if rescaled['nx'] > nx:
                rewindowed = EMAN2.Util.window(rescaled, nx, nx, 1)
            elif rescaled['nx'] < nx:
                sp_utilities.pad(rescaled, nx, nx, background='circumference')
            else:
                rewindowed = rescaled
                
            refi.append([rewindowed, ima.copy(), 0])
        # End resize if-then

    #  for each node read its share of data
    data = EMAN2.EMData.read_images(stack, list(range(image_start, image_end)))
    for im in range(image_start, image_end):
        data[im-image_start].set_attr('ID', im)
        if CTF:
            ctf_params = data[im-image_start].get_attr( "ctf" )
            if data[im-image_start].get_attr("ctf_applied") == 0:
                st = EMAN2_cppwrap.Util.infomask(data[im-image_start], mask, False)
                data[im-image_start] -= st[0]
                data[im-image_start] = sp_filter.filt_ctf(data[im-image_start], ctf_params)
                data[im-image_start].set_attr('ctf_applied', 1)
    if myid == main_node:  random.seed(rand_seed)

    a0 = -1.0
    again = True
    Iter = 0

    ref_data = [mask, center, None, None]

    while Iter < max_iter and again:
        ringref = []
        mashi = cnx-last_ring-2
        for j in range(numref):
            assert refi[j][0]['nx'] == mask['nx'], "ERROR!! Reference dimensions %s unequal to subject dimensions %s" % (
                refi[j][0]['nx'], mask['nx']
                )
            
            # normalize reference images to N(0,1)
            refi[j][0].process_inplace("normalize.mask", {"mask":mask, "no_sigma":1}) 
            
            cimage = EMAN2_cppwrap.Util.Polar2Dm(refi[j][0] , cnx, cny, numr, mode)
            EMAN2_cppwrap.Util.Frngs(cimage, numr)
            EMAN2_cppwrap.Util.Applyws(cimage, numr, wr)
            ringref.append(cimage)
            
            # zero refi
            refi[j][0].to_zero()
            refi[j][1].to_zero()
            refi[j][2] = 0
        
        if CTF: ctf2 = [[[0.0]*lctf for k in range(2)] for j in range(numref)]
        assign = [[] for i in range(numref)]
        
        # begin MPI section
        for im in range(image_start, image_end):
            alpha, sx, sy, mirror, scale = sp_utilities.get_params2D(data[im-image_start])
            #  Why inverse?  07/11/2015 PAP
            alphai, sxi, syi, scalei = sp_utilities.inverse_transform2(alpha, sx, sy)
            
            # normalize
            # subtract average under the mask
            data[im-image_start].process_inplace("normalize.mask", {"mask":mask, "no_sigma":0}) 
            # If shifts are outside of the permissible range, reset them
            if(abs(sxi)>mashi or abs(syi)>mashi):
                sxi = 0.0
                syi = 0.0
                sp_utilities.set_params2D(data[im-image_start],[0.0,0.0,0.0,0,1.0])
            
            ny = nx
            txrng = sp_alignment.search_range(nx, last_ring, sxi, xrng, "mref_ali2d_MPI")
            txrng = [txrng[1],txrng[0]]
            tyrng = sp_alignment.search_range(ny, last_ring, syi, yrng, "mref_ali2d_MPI")
            tyrng = [tyrng[1],tyrng[0]]
            
            # align current image to the reference
            [angt, sxst, syst, mirrort, xiref, peakt] = EMAN2_cppwrap.Util.multiref_polar_ali_2d(
                data[im-image_start], 
                ringref, 
                txrng, 
                tyrng, 
                step, 
                mode, 
                numr, 
                cnx+sxi, 
                cny+syi
                )
            
            iref = int(xiref)
            
            # combine parameters and set them to the header, ignore previous angle and mirror
            [alphan, sxn, syn, mn] = sp_utilities.combine_params2(
                0.0, -sxi, -syi, 0, 
                angt, sxst, syst, int(mirrort)
                )
            sp_utilities.set_params2D(data[im-image_start], [alphan, sxn, syn, int(mn), scale])
            data[im-image_start].set_attr('assign',iref)
            
            # apply current parameters and add to the average
            temp = sp_fundamentals.rot_shift2D(data[im-image_start], alphan, sxn, syn, mn)
            it = im%2
            EMAN2_cppwrap.Util.add_img( refi[iref][it], temp)
            assign[iref].append(im)
            if CTF:
                #  I wonder whether params are still there....
                ctf_params = data[im-image_start].get_attr("ctf")
                ctm = sp_morphology.ctf_2(nx, ctf_params)
                for i in range(lctf):  ctf2[iref][it][i] += ctm[i]
            
            refi[iref][2] += 1.0
        del ringref
        
        # end MPI section, bring partial things together, calculate new reference images, broadcast them back
        
        if CTF:
            # bring ctf2 together on main node
            s = numpy.shape(ctf2)
            ctf2  = mpi.mpi_reduce(
                ctf2, 2*lctf*numref, 
                mpi.MPI_FLOAT, 
                mpi.MPI_SUM, 
                main_node, 
                mpi.MPI_COMM_WORLD
                )
            if myid == main_node: ctf2 = numpy.reshape(ctf2, s)
        for j in range(numref):
            sp_utilities.reduce_EMData_to_root(refi[j][0], myid, main_node)
            sp_utilities.reduce_EMData_to_root(refi[j][1], myid, main_node)
            refi[j][2] = mpi.mpi_reduce(
                refi[j][2], 1, 
                mpi.MPI_FLOAT, 
                mpi.MPI_SUM, 
                main_node, 
                mpi.MPI_COMM_WORLD
                )
            if(myid == main_node): refi[j][2] = int(refi[j][2][0])
        
        # gather assignements
        for j in range(numref):
            if myid == main_node:
                for n in range(number_of_proc):
                    if n != main_node:
                        ln =  mpi.mpi_recv(
                            1,
                            mpi.MPI_INT, 
                            n, 
                            sp_global_def.SPARX_MPI_TAG_UNIVERSAL, 
                            mpi.MPI_COMM_WORLD
                            )
                        lis = mpi.mpi_recv(
                            ln[0], 
                            mpi.MPI_INT, 
                            n, 
                            sp_global_def.SPARX_MPI_TAG_UNIVERSAL, 
                            mpi.MPI_COMM_WORLD
                            )
                        for l in range(ln[0]): assign[j].append(int(lis[l]))
            else:
                mpi.mpi_send(
                    len(assign[j]), 
                    1, 
                    mpi.MPI_INT, 
                    main_node, 
                    sp_global_def.SPARX_MPI_TAG_UNIVERSAL, 
                    mpi.MPI_COMM_WORLD
                    )
                mpi.mpi_send(
                    assign[j], 
                    len(assign[j]), 
                    mpi.MPI_INT, 
                    main_node, 
                    sp_global_def.SPARX_MPI_TAG_UNIVERSAL, 
                    mpi.MPI_COMM_WORLD
                    )

        if myid == main_node:
            # replace the name of the stack with reference with the current one
            refim = os.path.join(outdir,"aqm%03d.hdf"%Iter)
            a1 = 0.0
            ave_fsc = []
            for j in range(numref):
                if refi[j][2] < 4:
                    #  if vanished, put a random image (only from main node!) there
                    assign[j] = []
                    assign[j].append( random.randint(image_start, image_end-1) - image_start )
                    refi[j][0] = data[assign[j][0]].copy()
                else:
                    if CTF:
                        for i in range(lctf):  
                            ctm[i] = old_div( 1.0, (ctf2[j][0][i] + old_div(1.0, snr) ) )
                        av1 = sp_filter.filt_table( refi[j][0], ctm)
                        for i in range(lctf):  
                            ctm[i] = old_div( 1.0, (ctf2[j][1][i] + old_div(1.0, snr) ) )
                        av2 = sp_filter.filt_table( refi[j][1], ctm)
                        frsc = sp_statistics.fsc(
                            av1, 
                            av2, 
                            1.0, 
                            os.path.join(outdir,"drm%03d_class%04d.txt"%(Iter, j) )
                            )
                        
                        # Now the total average
                        for i in range(lctf):  
                            ctm[i] = old_div( 1.0, (ctf2[j][0][i] + ctf2[j][1][i] + old_div( 1.0, snr ) ) )
                        refi[j][0] = sp_filter.filt_table( 
                            EMAN2_cppwrap.Util.addn_img( refi[j][0], refi[j][1] ), ctm
                            )
                    else:
                        frsc = sp_statistics.fsc(
                            refi[j][0], 
                            refi[j][1], 
                            1.0, 
                            os.path.join(outdir,"drm%03d_class%04d.txt"%(Iter,j) )
                            )
                        EMAN2_cppwrap.Util.add_img( refi[j][0], refi[j][1] )
                        EMAN2_cppwrap.Util.mul_scalar( refi[j][0], old_div(1.0, float(refi[j][2]) ) )
                            
                    if ave_fsc == []:
                        for i in range(len(frsc[1])): ave_fsc.append(frsc[1][i])
                        c_fsc = 1
                    else:
                        for i in range(len(frsc[1])): ave_fsc[i] += frsc[1][i]
                        c_fsc += 1

            if sum(ave_fsc) != 0:		
                for i in range(len(ave_fsc)):
                    ave_fsc[i]= old_div( ave_fsc[i], float(c_fsc) )
                    frsc[1][i]  = ave_fsc[i]
            
            for j in range(numref):
                ref_data[2]    = refi[j][0]
                ref_data[3]    = frsc
                refi[j][0], cs = user_func(ref_data)	

                # write the current average
                member_list = []
                for i_tmp in range(len(assign[j])): member_list.append(assign[j][i_tmp])
                member_list.sort()
                refi[j][0].set_attr_dict(
                    {'ave_n': refi[j][2],  'members': member_list, 'n_objects': len(member_list)}
                    )
                del member_list
                refi[j][0].process_inplace("normalize.mask", {"mask":mask, "no_sigma":1})
                refi[j][0].write_image(refim, j)
                a1 += refi[j][0].cmp("dot", refi[j][0], {"negative":0, "mask":mask})
                # send refi[j][0]  back!

            Iter += 1
            msg = "ITERATION #%3d        criterion = %20.7e\n"%(Iter,a1)
            sp_global_def.sxprint(msg)
            for j in range(numref):
                msg = "   group #%3d   number of particles = %7d"%(j, refi[j][2])
                sp_global_def.sxprint(msg)
            
            if a1 < a0:
                if (auto_stop == True):	again = False
            else:
                a0 = a1
        Iter  = sp_utilities.bcast_number_to_all(Iter, main_node)
        if CTF:  del  ctf2
        if again:
            for j in range(numref):
                sp_utilities.bcast_EMData_to_all(refi[j][0], myid, main_node)

    #  clean up
    del assign
    
    # write out headers  and STOP, under MPI writing has to be done sequentially
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    if CTF and data_had_ctf == 0:
        for im in range(len(data)): data[im].set_attr('ctf_applied', 0)
    par_str = ['xform.align2d', 'assign', 'ID']
    if myid == main_node:
        if(sp_utilities.file_type(stack) == "bdb"):
            sp_utilities.recv_attr_dict_bdb(
                main_node, stack, data, par_str, image_start, image_end, number_of_proc
                )
        else:
            sp_utilities.recv_attr_dict(
                main_node, stack, data, par_str, image_start, image_end, number_of_proc
                )
    else:           
        sp_utilities.send_attr_dict(main_node, data, par_str, image_start, image_end)

    if myid == main_node:
        newrefim = os.path.join(outdir, "multi_ref.hdf")
        for j in range(numref):
            refi[j][0].write_image(newrefim, j)
        sp_applications.header( stack, 'xform.align2d', fexport=os.path.join(outdir, 'docmrefparams.txt') )
        sp_global_def.sxprint("mref_ali2d_MPI")

def mref_ali2d():
    arglist = []
    for arg in sys.argv:
        arglist.append( arg )
    progname = os.path.basename(sys.argv[0])
    usage = progname + " data_stack reference_stack outdir <maskfile> --ir=inner_radius --ou=outer_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translation_step --center=center_type --maxit=max_iteration --CTF --snr=SNR --function=user_function_name --rand_seed=random_seed --MPI"
    parser = optparse.OptionParser(usage,version=sp_global_def.SPARXVERSION)
    parser.add_option(
        "--ir", 
        type="float", 
        default=1, 
        help="  inner radius for rotational correlation > 0 (set to 1)"
        )
    parser.add_option(
        "--ou", 
        type="float", 
        default=-1, 
        help="  outer radius for rotational correlation < nx/2-1 (set to the radius of the particle)"
        )
    parser.add_option(
        "--rs", 
        type="float", 
        default=1, 
        help="  step between rings in rotational correlation > 0 (set to 1)" )
    parser.add_option(
        "--xr", 
        type="float", 
        default=0, 
        help="  range for translation search in x direction, search is +/-xr "
        )
    parser.add_option(
        "--yr", 
        type="float",
        default=0, 
        help="  range for translation search in y direction, search is +/-yr "
        )
    parser.add_option(
        "--ts", 
        type="float", 
        default=1, 
        help="  step of translation search in both directions"
        )
    parser.add_option(
        "--center", 
        type="float", 
        default=1, 
        help="  0 - if you do not want the average to be centered, 1 - center the average (default=1)"
        )
    parser.add_option(
        "--maxit", 
        type="float", 
        default=10, 
        help="  maximum number of iterations (set to 10) "
        )
    parser.add_option(
        "--CTF", 
        action="store_true", 
        default=False, 
        help=" Consider CTF correction during multiple reference alignment"
        )
    parser.add_option(
        "--upscale_isac_dir", 
        type="str", 
        default=None, 
        help=" Upscale reference images from ISAC class-average directory"
        )
    parser.add_option(
        "--snr", 
        type="float",  
        default= 1.0, 
        help="  signal-to-noise ratio of the data (set to 1.0)"
        )
    parser.add_option(
        "--function", 
        type="string", 
        default="ref_ali2d", 
        help="  name of the reference preparation function"
        )
    parser.add_option(
        "--rand_seed", 
        type="int", 
        default=1000, 
        help=" random seed of initial (set to 1000)" 
        )
    parser.add_option(
        "--MPI", 
        action="store_true", 
        default=False,     
        help="  whether to use MPI version "
        )
    
    (options, args) = parser.parse_args(arglist[1:])
    
    if len(args) < 3 or len(args) > 4:
            sp_global_def.sxprint( "Usage: " + usage )
            sp_global_def.sxprint( "Please run '" + progname + " -h' for detailed options" )
            sp_global_def.ERROR( "Invalid number of parameters used. Please see usage information above." )
            return
    else:
    
        if len(args) == 3:
            mask = None
        else:
            mask = args[3]

        if sp_global_def.CACHE_DISABLE:
            sp_utilities.disable_bdb_cache()
        
        sp_global_def.BATCH = True
        
        mref_ali2d_MPI(
            args[0], 
            args[1], 
            args[2], 
            mask, 
            options.ir, 
            options.ou, 
            options.rs, 
            options.xr, 
            options.yr, 
            options.ts, 
            options.center, 
            options.maxit, 
            options.CTF, 
            options.snr, 
            options.function, 
            options.rand_seed, 
            options.upscale_isac_dir
            )
        sp_global_def.BATCH = False

def main():
    sp_global_def.print_timestamp( "Start" )
    sp_global_def.write_command()
    mref_ali2d()
    sp_global_def.print_timestamp( "Finish" )
    mpi.mpi_finalize()

if __name__ == "__main__":
    main()
