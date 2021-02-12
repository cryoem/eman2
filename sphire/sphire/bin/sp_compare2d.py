#!/usr/bin/env python
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
# Copyright (c) 2012 The University of Texas - Houston Medical School
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
from __future__ import division
from past.utils import old_div
import os
import EMAN2
import EMAN2_cppwrap
import datetime
import argparse
import math
try:
    from sphire.libpy import sp_global_def
    from sphire.libpy import sp_utilities 
    from sphire.libpy import sp_statistics 
    from sphire.libpy import sp_logger 
    from sphire.libpy import sp_alignment
    from sphire.libpy import sp_fundamentals
    from sphire.libpy import sp_applications
    from sphire.libpy import sp_filter
except ImportError:
    import sp_global_def
    import sp_utilities 
    import sp_statistics 
    import sp_logger 
    import sp_alignment
    import sp_fundamentals
    import sp_applications
    import sp_filter

# Set default values (global variables written in ALL CAPS)
TIMESTAMP_LENGTH= 23  # chars

REFANGLES= 'refangles.txt'
ALIGNDOC=  'docalign2d.txt'
MASKFILE=  'maskccc.hdf'
STACKFILE= 'stkcompare2d.hdf'

USAGE= """
PURPOSE:
Compare two image stacks and display the best match side-by-side.

SYNTAX:
%s <experimental_images> <reference_images> {output_directory>
General, required command-line parameters:
1. Input experimental image stack
2. Input reference image stack, same xy-dimensions as first stack
3. Output directory
Outputs:
docalign2d.txt : Alignment parameters
comparison2d.hdf : Stack of experimental images next to aligned reference images

More options:
%s <experimental_images> <reference_images> {output_directory> --maxshift=1 --outerradius=29 --ringstep=1
Parameters:
--maxshift : Maximum shift to allow during translation alignment (default 2)
--outerradius : Outer alignment radius (default Automatic)
--ringstep : Alignment radius step size (default 1)
--verbosity : verbosity level (0..2)

""" % ((__file__,)*2)
    
MODIFIED= "Modified 2021-01-19"

"""
CHANGELOG:
    2021-01-19 (trs) -- added normalization option
    2021-01-19 (trs) -- added CCC (correctly) and best reference to comparison stack
    2021-01-19 (trs) -- updated for Python 3
    2020-06-18 (trs) -- fixed alignment normalization bug
    2020-04-01 (trs & dr) -- updated usage
    2020-04-01 (trs) -- applies (inverse of) alignment parameters to reference images
"""

def check(file, log=None, verbose=False):
    """
    Checks whether file exists.
    
    Arguments:
        file : File to look for
        log : Instance of Logger object
        verbose : (boolean) Whether to write to screen
    """
    
    if os.path.exists(file):
        if verbose: print_log_msg("Found '%s'" % file, log)
    else:
        print_log_msg("ERROR!! '%s' doesn't exist!\n" % log)
        exit()

def prepare_outdir_log(outdir='.', verbose=False, is_main=True):
    """
    Prepares output directory and sets up log file.
    
    Arguments:
        outdir : output directory
        verbose : (boolean) whether to write to screen
        is_main : (boolean) if using multiple cores, some tasks only need to be performed once, not by all cores
    Returns:
        log : instance of Logger class
        verbose : (boolean) new version of Logger can write to screen simultaneously
    """
    
    # Create directory if it doesn't exist
    if is_main:
        if os.path.isdir(outdir):
            print("Writing to output directory: %s" % outdir)
        else:
            print("Created output directory: %s" % outdir)
            os.makedirs(outdir)  # os.mkdir() can only operate one directory deep
            
    logname= "logfile.txt"
    logname= os.path.join(outdir, logname)
    
    # May be using old version of logger.py
    try:
        if verbose:
            log= sp_logger.Logger(
                base_logger=sp_logger.BaseLogger_Files(), 
                base_logger2=sp_logger.BaseLogger_Print(), 
                file_name=logname
                )
            verbose= False  # logger output will be echoed to screen
        else:
            log= sp_logger.Logger(base_logger=sp_logger.BaseLogger_Files(), file_name=logname)
    except TypeError:
        if is_main: print("WARNING: Using old logger.py library")
        log= sp_logger.Logger(base_logger=sp_logger.BaseLogger_Files())#, file_name=logname)
        logname= 'log.txt'
        
    if is_main:
        print("Writing log file to %s" % logname)
        progbase= os.path.basename(__file__).split('.')[0].upper()
        length= len(progbase) + 4
        
        log.add("\n" +
                " "*TIMESTAMP_LENGTH + "*"*length + "\n" +
                " "*TIMESTAMP_LENGTH + "* " + progbase + " *\n" +
                " "*TIMESTAMP_LENGTH + "*"*length)
    
    return log, verbose

def print_log_msg(mesg, log=None, verbose=False, is_main=True):
    """
    Prints messages to log file and, optionally, to the screen.
    
    Arguments:
        mesg : Message to write
        log : Instance of Logger class
        verbose : (boolean) Whether to write to screen
        is_main : (boolean) If using MPI, some tasks only need to be performed once, not by all cores
    """
    
    if is_main:
        if verbose: sp_global_def.sxprint(mesg)
        if log: log.add(mesg)

def compare2d(expstack, refstack, outdir, options):
    """
    Main function overseeing various projection-comparison modes.
    
    Arguments:
        expstack : Input experimental image stack
        refstack : Input reference image stack
        outdir : Output directory
        options : (list) Command-line options, run 'sxproj_compare.py -h' for an exhaustive list
    """
    
    # Expand path for outputs
    refanglesdoc= os.path.join(outdir, REFANGLES)
    outaligndoc= os.path.join(outdir, ALIGNDOC)

    # Check that dimensions of images and volume agree (maybe rescale volume)
    refdim= EMAN2.EMData(refstack, 0).get_xsize()
    imgdim= EMAN2.EMData(expstack,0).get_xsize()
    if refdim != imgdim:
        sp_global_def.ERROR("\nERROR!! Dimensions of input stacks don't match: %s vs. %s" % 
                (refdim, imgdim), __file__, 1)
        exit()
    
    # Set output directory and log file name
    log, _= prepare_outdir_log(outdir)
    sp_global_def.write_command(outdir)

    if options.verbosity>=1:
        print_log_msg("Finding best match between '%s' and '%s'" % (expstack, refstack), log)
        
    # Check if inputs exist
    check(expstack, log=log, verbose=options.verbosity>=1)
    check(refstack, log=log, verbose=options.verbosity>=1)
    
    if options.normalize==None:
        mesg= "Not normalizing images"
    elif options.normalize=='minmax' or options.normalize=='rops' or options.normalize=='sigmaone':
        mesg= "Normalization mode: '%s'" % options.normalize
    else:
        sp_global_def.ERROR("\nNormalization mode '%s' unknown!" % options.normalize, __file__, 1)
        exit()
    if options.verbosity>=1 : print_log_msg(mesg, log)
    
    # Perform multi-reference alignment
    align_params= apsh(
        refstack, 
        expstack, 
        outangles=None, 
        refanglesdoc=None, 
        outaligndoc=outaligndoc, 
        outerradius=options.outerradius, 
        maxshift=options.maxshift, 
        ringstep=options.ringstep, 
        log=log, 
        verbose=options.verbosity>=1
        )
    
    # Set filenames
    maskfile= os.path.join(outdir, MASKFILE)
    comparison_stack= os.path.join(outdir, STACKFILE)
    
    # Make comparison stack between two sets of images, e.g., experimental images and class averages
    compare_imgs(
        refstack, 
        expstack, 
        align_params, 
        comparison_stack, 
        maskfile= maskfile,
        norm_mode=options.normalize, 
        log=log, 
        verbosity=options.verbosity
        )

def apsh(refstack, imgs2align, outangles=None, refanglesdoc=None, outaligndoc=None, 
        outerradius=-1, maxshift=0, ringstep=1, mode="F", log=None, verbose=False):
    """
    Generates polar representations of a series of images to be used as alignment references.
    
    Arguments:
        refstack : Input reference image stack (filename or EMData object)
        imgs2align : Image stack to be aligned (filename or EMData object)
        outangles : Output Euler angles doc file
        refanglesdoc : Input Euler angles for reference projections
        outaligndoc : Output 2D alignment doc file
        outerradius : Outer alignment radius
        maxshift : Maximum shift allowed
        ringstep : Alignment radius step size
        mode : Mode, full circle ("F") vs. half circle ("H")
        log : Logger object
        verbose : (boolean) whether to write to screen
    """
    
    # Generate polar representation(s) of reference(s)
    alignrings, polarrefs= mref2polar(refstack, outerradius=outerradius, ringstep=ringstep, 
            log=log, verbose=verbose)
            
    # Read image stack (as a filename or already an EMDataobject)
    if isinstance(imgs2align, str): 
        imagestacklist= EMAN2.EMData.read_images(imgs2align)
    else:
        imagestacklist= [imgs2align]
    
    # Get number of images
    numimg= len(imagestacklist)
    
    # Get image dimensions (assuming square, and that images and references have the same dimension)
    idim= imagestacklist[0]['nx']

    # Calculate image center
    halfdim= idim//2 + 1
    
    # Set constants
    currshift= 0
    scale= 1
    
    # Initialize output angles
    outangleslist= []
    outalignlist= []
    
    if outerradius <= 0:
        outerradius= halfdim - 3
        
    if refanglesdoc:
        refangleslist= sp_utilities.read_text_row(refanglesdoc)
        
    # Set search range
    txrng= tyrng= sp_alignment.search_range(idim, outerradius, currshift, maxshift)
    
    if verbose : 
        print_log_msg('Running multireference alignment allowing a maximum shift of %s' % maxshift, log, False)
    
    # Loop through images
    for imgindex in range(numimg):
        currimg= imagestacklist[imgindex]
        
        # Perform multi-reference alignment (adapted from alignment.mref_ali2d)
        [angt, sxst, syst, mirrorfloat, bestreffloat, peakt]= EMAN2_cppwrap.Util.multiref_polar_ali_2d(
            currimg, polarrefs, txrng, tyrng, ringstep, mode, alignrings, halfdim, halfdim)
        bestref= int(bestreffloat)
        mirrorflag= int(mirrorfloat)
        
        # Store parameters
        params2dlist= [angt, sxst, syst, mirrorflag, scale, bestref]
        outalignlist.append(params2dlist)
    
        if refanglesdoc:
            besteulers= refangleslist[bestref]
        else:
            besteulers= [0]*5
        
        # Check for mirroring
        if mirrorflag== 1:
            tempeulers= list(sp_utilities.compose_transform3(
                besteulers[0],besteulers[1],besteulers[2], besteulers[3],besteulers[4],0, 1,
                0,180,0, 0,0,0, 1)
            )
            combinedparams= list(sp_utilities.compose_transform3(
                tempeulers[0],tempeulers[1],tempeulers[2], tempeulers[3],tempeulers[4],0, 1, 
                0,0,-angt, 0,0,0, 1)
            )
        else:
            combinedparams= list(sp_utilities.compose_transform3(
                besteulers[0],besteulers[1],besteulers[2], besteulers[3],besteulers[4],0, 1, 
                0,0,-angt, 0,0,0, 1)
            )
        # (compose_transform3: returns phi,theta,psi, tx,ty,tz, scale)
        
        outangleslist.append(combinedparams)
        
    if outangles or outaligndoc:
        mesg= ''
        if outangles : 
            sp_utilities.write_text_row(outangleslist, outangles)
            if verbose:
                mesg += "Wrote alignment angles to '%s'\n" % outangles
                print_log_msg(mesg, log, False)
        if outaligndoc : 
            sp_utilities.write_text_row(outalignlist, outaligndoc)
            if verbose:
                mesg += "Wrote 2D alignment parameters to '%s'\n" % outaligndoc
                print_log_msg(mesg, log, False)
    else:
        if verbose : print_log_msg('Finished running multireference alignment\n', log, False)
    
    return outalignlist
    
def mref2polar(refstack, firstring=1, outerradius=-1, ringstep=1, mode="F", normbysquare=0, 
               log=None, verbose=False):
    """
    Generates polar representations of a series of images to be used as alignment references.
    
    Arguments:
        refstack : Input reference image stack (filename or EMData object)
        firstring : Inner alignment radius
        outerradius : Outer alignment radius
        ringstep : Alignment radius step size
        mode : Mode, full circle ("F") vs. half circle ("H)
        normbysquare : If other than 0, normalization by setting the norm to 1
        log : Logger object
        verbose : (boolean) whether to write to screen
    Returns:
        alignringlist : List of alignment-ring data
        polarreflist : List of polar representation of refernences
    """
    
    # Read reference stack
    if isinstance(refstack, str):
        referencelist= EMAN2.EMData.read_images(refstack)
    else:
        referencelist= [refstack]  # For single image
    
    numrefs= len(referencelist)
    polarreflist= []
    
    # Get image dimensions (assuming square, and that images and references have the same dimension)
    idim= referencelist[0]['nx']

    # Calculate image center
    halfdim= idim//2
    
    if outerradius <= 0:
        outerradius= halfdim - 2
        
    # Prepare alignment rings
    alignringlist= sp_alignment.Numrinit(firstring, outerradius, ringstep, mode)
    
    # Calculate ring weights
    ringweightlist= sp_alignment.ringwe(alignringlist, mode)
    
    if verbose : 
        mesg= 'Converting %s references to polar coordinates from radius %s to %s with step %s and mode "%s"' % \
            (numrefs, firstring, outerradius, ringstep, mode)
        print_log_msg(mesg, log, verbose=False)#, log)
    
    # Loop through reference images (adapted from sxisac2)
    for refindex in range(numrefs):
        # Convert to polar
        cimage= EMAN2_cppwrap.Util.Polar2Dm(
            referencelist[refindex], 
            halfdim, 
            halfdim, 
            alignringlist, 
            mode
            )
        
        # Normalize
        EMAN2_cppwrap.Util.Normalize_ring(cimage, alignringlist, normbysquare)  
        #normbysquare: if other than 0, normalizes by setting the norm to 1
        
        # Fourier transform of rings
        EMAN2_cppwrap.Util.Frngs(cimage, alignringlist)

        # Apply weights to rings
        EMAN2_cppwrap.Util.Applyws(cimage, alignringlist, ringweightlist)

        # Copy to reference stack
        polarreflist.append(cimage.copy())
        
    return alignringlist, polarreflist

def compare_imgs(refstack, expstack, align_params, compstack, maskfile=None, norm_mode=None, 
                 log=None, verbosity=0):
    """
    Make comparison stack between class averages and re-projections.
    
    Arguments:
        refstack : Input reference image stack
        expstack ; Input experimental image stack
        align_params : Input alignment parameter list
        compstack : Output comparison stack
        maskfile : Mask file for computing correlation
        norm_mode : Normalization mode for comparison images
        log : Logger object
        verbosity : Verbosity level (0..2)
    """
    
    nimg1  = EMAN2.EMUtil.get_image_count(expstack)
    
    # Get image dimension
    recondata= EMAN2.EMData(refstack, 0)
    nx= recondata.get_xsize()
    
    # Initialize list of cross-correlation coefficients
    ccclist= []
    
    #  Here you need actual radius to compute proper ccc's, but if you do, you have to deal with translations, PAP
    mask_obj= sp_utilities.model_circle(nx//2-2,nx,nx)
    if maskfile : mask_obj.write_image(maskfile)
    
    if verbosity>=1 : print_log_msg("Creating montage of images", log)
    
    for imgnum in range(nimg1):
        # Get experimental image
        expimg= sp_utilities.get_im(expstack, imgnum)
        
        # Get reference image number (last column in alignment parameters)
        [angt, sxst, syst, mirrorflag, scale, bestref]= align_params[imgnum]
        
        # Get reference image
        refimg= sp_utilities.get_im(refstack, bestref)
        
        # apply 2D alignment
        if mirrorflag== 1 : refimg= sp_fundamentals.mirror(refimg, axis='x')
        angt, sxst, syst, _= sp_utilities.combine_params2(0,-sxst,-syst,0, -angt,0,0,0)
        refimg= sp_fundamentals.rot_shift2D(
            refimg, 
            angle=angt, 
            sx=sxst, 
            sy=syst, 
            interpolation_method='linear'
            )
        
        cccoeff= sp_statistics.ccc(refimg, expimg, mask_obj)
        if verbosity>= 2 : 
            mesg= "  Image %s, bestref %s mirror %s, CCC %.4f" % (imgnum, bestref, mirrorflag, cccoeff)
            print_log_msg(mesg, log)
            
        montagestack= [expimg, refimg]
        comparison_pair= montage_norm(montagestack, ncol=2, marginwidth=1, norm_mode=norm_mode)
        comparison_pair.set_attr_dict({'compare2d.ccc': cccoeff})
        comparison_pair.set_attr_dict({'compare2d.best': bestref})
        comparison_pair.write_image(compstack, imgnum)
        
        ccclist.append(cccoeff)
        
    meanccc= sum(ccclist)/nimg1
    if verbosity>=1 : print_log_msg("Average CCC is %.5f" % meanccc, log)
    
    nimg2= EMAN2.EMUtil.get_image_count(compstack)
    
    for imgnum in range(nimg2):  # xrange will be deprecated in Python3
        refimg= sp_utilities.get_im(compstack,imgnum)
        meanccc1= refimg.get_attr_default('mean-cross-corr', -1.0)
        refimg.set_attr_dict({'mean-cross-corr':meanccc})
        sp_utilities.write_header(compstack,refimg,imgnum)
    
    if verbosity>=1 : print_log_msg("Wrote %s images to '%s'\n" % (nimg2, compstack), log, False)
    
def montage_norm(input_stack_or_list, ncol=None, marginwidth=0, bkgd=0, outfile=None, norm_mode=None):
    """
    Generates montage of images into one image.
    Adapted from sxmontage.py and SPIDER's MN S
    Modified from in sp_utilities in SPHIRE 1.3
    2021-01-19 -- Added normalization option: avg->0, sigma->1
    2021-01-19 -- Added rotational power-spectrum normalization (adapted from sp_proj_compare)
    
    Arguments:
        input_stack_or_list : Stack of input images to merge into montage
        ncol : Number of images per row (default: all on one row)
        marginwidth : Margin width, pixels
        bkgd : Background value of montage
        outfile : Optional output file with montage output
        norm_mode: Normalization mode: None, minmax, rops, sigmaone
    Returns:
        montage : EMData object of image montage
    """
    
    if isinstance(input_stack_or_list, str): 
        image_list= EMAN2.EMData.read_images(input_stack_or_list)
    else:
        assert isinstance(input_stack_or_list, list), "Uh oh! Don't know type %s" % type(input_stack_or_list)
        image_list= input_stack_or_list
    
    # Get single-image dimensions
    nx= image_list[0].get_xsize()
    ny= image_list[0].get_ysize()
    
    # Get number of images and calculate montage dimensions
    numimgs= len(image_list)
    if ncol != None:
        numrows= old_div( (numimgs-1), ncol) + 1
    else:
        ncol= numimgs
        numrows= 1
    
    # Create blank image
    montage_xdim= (nx + marginwidth)*ncol
    montage_ydim= (ny + marginwidth)*numrows
    montage= sp_utilities.model_blank(montage_xdim, montage_ydim, 1, bkgd)
    mask= sp_utilities.model_blank(nx, ny, bckg=1.0)
    
    # Loop through images
    for img_num in range(numimgs):
        # Horizontal grid position is image# modulo NCOL
        col_num= img_num % ncol
        
        # Montage is numbered from the top down
        row_num= numrows - 1 - old_div(img_num, ncol)
        
        xoffset= col_num*(nx+marginwidth)
        yoffset= row_num*(ny+marginwidth)
        
        if norm_mode== None:
            img_norm= image_list[img_num]
        elif norm_mode== 'minmax':
            currentx= image_list[img_num].get_xsize()
            assert currentx== nx, "ERROR!! Mismatched dimensions %s != %s" % (currentx, nx)
            [avg,var,fmin,fmax]= EMAN2.Util.infomask(image_list[img_num], mask, True)
            
            img_norm= (image_list[img_num] - fmin) * old_div(2., fmax - fmin)
        elif norm_mode== 'rops':
            if img_num==0:
                rops_src= sp_fundamentals.rops_table(image_list[img_num]*mask)
                img_norm= image_list[img_num]
            else:
                rops_dst= sp_fundamentals.rops_table(image_list[img_num]*mask)  
                
                #  Set power spectrum of reprojection to the data.
                #  Since data has an envelope, it would make more sense to set data to reconstruction,
                #  but to do it one would have to know the actual resolution of the data. 
                #  you can check sxprocess.py --adjpw to see how this is done properly  PAP
                roo_list= [0.0]*len(rops_dst)  # initialize table
                for j in range( len(rops_dst) ):
                    roo_list[j]= math.sqrt( rops_dst[j]/rops_src[j] )
                
                # match FFT amplitudes
                img_norm= sp_filter.filt_table(image_list[img_num], roo_list)
        elif norm_mode== 'sigmaone':
            img_stat= EMAN2_cppwrap.Util.infomask(image_list[img_num]*mask, mask, False)
            img_mask= image_list[img_num]*mask
            img_zero= img_mask - img_stat[0]
            img_norm= old_div(img_zero, img_stat[1])
        else:
            sp_global_def.ERROR("\nNormalization mode '%s' unknown!" % norm_mode, __file__, 1)
            exit()
        insert_image(img_norm, montage, xoffset, yoffset)
    # End image-loop
    
    if outfile: montage.write_image(outfile)
    
    return montage
    
def insert_image(smallimg, largeimage, xoffset, yoffset):
    """
    Inserts small image into large image.
    Adapted from sxmontage.py
    
    Arguments:
        smallimg : Small image to insert into large image
        largeimage : Large image (OVERWRITTEN!) into which small image will be inserted
        xoffset : Top-left x-coordinate of large image where small image will be inserted
        yoffset : Top-left y-coordinate of large image where small image will be inserted
    """
    
    # Get small-image dimensions
    nx= smallimg.get_xsize()
    ny= smallimg.get_ysize()
    
    for xcoord in range(nx):
        for ycoord in range(ny):
            getpixel= smallimg.get_value_at(xcoord, ycoord)
            largeimage.set_value_at(xoffset+xcoord, yoffset+ycoord, getpixel)

def parse_command_line():
    """
    Parse the command line.  Adapted from sxmask.py

    Arguments:
    None:

    Returns:
    Parsed arguments object
    """

    parser= argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
        usage=USAGE, 
        epilog=MODIFIED
        )

    parser.add_argument(
        'expimgs', 
        type=str, 
        help='Input experimental image stack')
    
    parser.add_argument(
        'refimgs', 
        type=str, 
        help='Input reference image stack')
    
    parser.add_argument(
        'outdir', 
        type=str, 
        default='.', 
        help='Output directory (default .)')
    
    parser.add_argument(
        '--maxshift', 
        type=int, 
        default=2, 
        help='Maximum translational shift during multireference alignment (pixels, default 2)')
    
    parser.add_argument(
        '--outerradius', 
        type=int, 
        default=-1, 
        help='Outer radius for alignment (pixels, default automatic)')
    
    parser.add_argument(
        '--ringstep', 
        type=int, 
        default=1, 
        help='Radius increment during alignment (pixels, default 1)')

    parser.add_argument(
        '--normalize', '--norm',
        type=str, 
        default=None, 
        help='Normalization mode: None, minmax, rops, sigmaone')

    parser.add_argument(
        '--verbosity', '-v', 
        type=int, 
        default=1, 
        help='Extent to which to write information to the screen (0..2)')

    return parser.parse_args()

def main():
    sp_global_def.print_timestamp( "Start" )
    options= parse_command_line()
    
    ##print args, options  # (Everything is in options.)
    #print('options', options)
    #print('LOGFILE', sp_global_def.LOGFILE)
    #exit()
    
    main(options.expimgs, options.refimgs, options.outdir, options)
    sp_global_def.print_timestamp( "Finish" )

if __name__== "__main__":
    main()
