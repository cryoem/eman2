from past.utils import old_div  #### old_div(
import numpy as np
import sys
import argparse
import os

import EMAN2
import datetime
import textwrap

from sphire.libpy import sp_global_def
from sphire.libpy import sp_utilities
from sphire.libpy import sp_fundamentals
from sphire.libpy import sp_morphology
from sphire.libpy import sp_projection
from sphire.libpy import sp_logger
from sphire.libpy import sp_applications

# Set default values (global variables written in ALL CAPS)
TIMESTAMP_LENGTH= 23  # chars
HEADER_LABEL= 'xform.center2d'

# Output filenames
SHIFTONLY_DOC= 'docshifts.txt'
MASK3D_CENT= 'mask_centeredbox.hdf'
VOL_CENT_BY_BOX= 'vol_box_centered.hdf'
VOL_CENT_BY_BLOB= 'vol_blob_centered.hdf'
MASK_STACK= 'stkmask.hdf'
MASK_APPLIED= 'stkappliedmask.hdf'
MASK_MONTAGE= 'stkcompare.hdf'
MASK2D_CENT= 'stkmask_centered.hdf'
AVG_SHIFT_BY_ACF= 'stkavgs_centered.hdf'
AVG_SHIFT_BY_RT180= 'stkcentered_rt180.hdf'

# Additional outputs in DEBUG mode
MASK3D_PADDED= 'tst3d_1pad.hdf'
MASK3D_PATTERSON= 'tst3d_2patterson.hdf'
MASK3D_ROTATED= 'tst3d_3rot.hdf'
MASK3D_ZSHIFT= 'tst3d_99shiftz.hdf'

PROJ2D_XY= 'tst2d_0proj_xy.hdf'
PROJ2D_PAD= 'tst2d_1proj_pad.hdf'
PROJ2D_ACF= 'tst2d_2acf.hdf'
PROJ2D_ACF_THRESH= 'tst2d_3acf_thresh.hdf'
PROJ2D_ROT= 'tst2d_4proj_rot.hdf'
PROJ2D_ROT_THRESH= 'tst2d_5proj_rot_thresh.hdf'
PROJ2D_ROTSHIFT= 'tst2d_6proj_rotshift.hdf'
PROJ2D_SHIFTONLY= 'tst2d_99proj_shiftxy.hdf'

BLOB_CENT= 'tstvol1_blob.hdf'

MASK2D_PAD= 'tstmask1_pad.hdf'
MASK2D_ACF= 'tstmask2_acf.hdf'
MASK2D_ACF_THRESH= 'tstmask3_acf_thresh.hdf'
MASK2D_ROT= 'tstmask4_rot.hdf'
MASK2D_ROT_THRESH= 'tstmask5_rot_thresh.hdf'
MASK2D_ROTSHIFT= 'tstmask6_rotshift.hdf'
MASK2D_SHIFTONLY= 'tstmsk7_shiftxy.hdf'

RT180_ROT= 'tst1_rot.hdf'
RT180_CCF= 'tst2_ccf.hdf'

USAGE= """
PURPOSE:
Various centering tools for 2D and 3D.

General options:
%s <input> <output_directory> --mode=<mode>
Required command-line parameters:
1. Input volume/stack (sometimes binarized, see detailed usage below)
2. Output directory
Parameters:
--mode :  mode, allowed options are: acf (default, 2D or 3D), rt180 (2D), boundingbox (3D), blob (3D)
Advanced parameters:
--dimensions : in case it is ambiguous whether the input is 2D or 3D (e.g., MRC format), enter 2 or 3, respectively
--verbosity : verbosity level
--debug : writes various intermediate diagostic files

2D options:

To center a series of binarized images using the ACF:
%s <input_image_stack> <output_directory> --mode=acf --threshsigma=0.7 --ndilation=3 --nerosion=0
Parameters:
--threshsigma : binarization threshold, in units of sigma above the mean (default 0)
Advanced parameters:
--ndilation : number of dilations of binary mask (default 3)
--nerosion : number of erosions of binary mask (default 0)

To center a series of images using 180-degree rotation:
%s <input_image_stack> <output_directory> --mode=rt180

3D options:

To center a binarized volume using an autocorrelation function:
%s <input_binary_mask> <output_directory> --mode=acf --applyshift <input_map>
Parameters:
--applyshift : volume to which to apply shifts 
Advanced parameters:
--threshint : intermediate threshold for binarization (0..1, default 0.1)
--pad : padding factor, to prevent wraparound artifacts in ACF (default 2)
--stem : stem to embed in output volume names (default acf)

To center a binarized volume using a bounding box:
%s <input_binary_mask> <output_directory> --mode=boundingbox
Parameters:
--applyshift : volume to which to apply shifts 

To center a volume by alignment to a Gaussian blob:
%s <input_volume> <output_directory> --mode=blob --volradius=<radius>
--volradius : Radius of blob to use for alignment

""" % ( (__file__,)*6) 

MODIFIED= "Modified 2021-02-18"

"""
Modifications log:
    2021-02-18 (trs) -- Updated for Python 3
    2020-05-22 (trs) -- Renamed the header tag xform.center2d, for use with other programs
    2020-04-21 (trs & dq) -- In 2D ACF mode, uses scipy masking to disallow disconnected mask
    2020-04-16 (trs) -- In masking mode, checks for empty images
    2020-04-10 (trs) -- Writes shifts to header, in 2D modes
    2020-03-22 (trs) -- Adapted blob-alignment from signal-subtraction
    2020-03-22 (trs) -- Input 3D mask is binarized (at 0.9) in case it was an adaptive mask
    TODO: Check whether inputs are binary when they're supposed to be
    TODO: In 3D mode, give dimensions in Angstroms, which may be informative
"""

def center_2d3d(inputfile, mode, outdir, options=None):
    """
    Main function overseeing various centering options.
    
    Arguments:
        inputfile : Stack of input images to merge into montage
        mode : Mode: allowed options are acf (2D or 3D), boundingbox (3D), acf (2D), or rt180 (2D)
        outdir : Output directory
        options : (list) Command-line options, run 'sp_center_2d3d.py -h' for an exhaustive list
    """
    
    # Set output directory and log file name
    log, _= prepare_outdir_log(outdir)

    # Read input
    input_obj= EMAN2.EMData(inputfile)
    
    # Check if volume or image(s), choices are cubic volume, non-MRC stack, and MRC stack
    xdim= input_obj['nx']
    ydim= input_obj['ny']
    zdim= input_obj['nz']
    
    if xdim != ydim:
        sp_global_def.ERROR("Program expects square/cubic inputs", __file__, 1)
        sp_global_def.print_timestamp( "Finish" )
        exit()
    
    # If a volume or MRC stack
    if zdim > 1:
        if zdim == xdim:
            if options.dimensions == 3:
                if options.verbosity >= 1 : print_log_msg("Reading '%s' as a volume\n" % inputfile, log)
            else:
                options.dimensions= 2
        # If not a cube, then assuming a stack
        else:
            options.dimensions= 2
    elif zdim == 1:
        options.dimensions= 2
    else:
        sp_global_def.ERROR("Reading unrecognizable z-dimension of %s" % zdim, __file__, 1)
        sp_global_def.print_timestamp( "Finish" )
        exit()
        
    if options.dimensions == 3:
        if options.volradius:
            if mode != 'blob':
                mode= 'blob'
                print_log_msg('WARNING: volradius specified, so assuming --mode=blob', log)
                    
        if mode == 'acf':
            if options.verbosity >=1 :
                print_log_msg('Centering using autocorrelation method', log)
                if options.applyshift: 
                    assert os.path.exists(options.applyshift), "ERROR!! %s doesn't exist!" % options.applyshift
                    print_log_msg("Shifts will be applied to '%s'" % options.applyshift, log)
                print_log_msg('Intermediate binary files will be thresholded at %s' % options.threshint, log)
                print_log_msg('Will pad by a factor of %s before computing autocorrelation' % options.pad, log)
                print_log_msg("Output volume names will be named '%s'" % ('vol_' + options.stem), log)
                print_log_msg('Verbosity level to set to: %s' % options.verbosity, log)
                print_log_msg('Debugging diagnostics set to: %s\n' % options.debug, log)
            
            center_patterson_3d(
                input_obj, 
                outdir=outdir, 
                vol2shift=options.applyshift, 
                threshold=options.threshint, 
                pad=options.pad, 
                stem=options.stem, 
                verbosity=options.verbosity, 
                debug=options.debug, 
                log=log
                )
        
        elif mode == 'boundingbox':
            if options.verbosity >=1 :
                print_log_msg('Centering using bounding-box method', log)
                if options.applyshift: 
                    assert os.path.exists(options.applyshift), "ERROR!! %s doesn't exist!" % options.applyshift
                    print_log_msg("Shifts will be applied to '%s'" % options.applyshift, log)
                print_log_msg('Intermediate binary files will be thresholded at %s' % options.threshint, log)
                print_log_msg('Verbosity level to set to: %s' % options.verbosity, log)
                print_log_msg('Debugging diagnostics set to: %s\n' % options.debug, log)
            
            center_boundingbox(
                input_obj, 
                outdir=outdir, 
                vol2shift=options.applyshift, 
                threshold=options.threshint, 
                verbosity=options.verbosity, 
                debug=options.debug, 
                log=log
                )
        
        elif mode == 'blob':
            if not options.volradius:
                sp_global_def.ERROR("Volume radius '--volradius' expected for blob mode", __file__, 1)
            else:
                if options.verbosity >=1 :
                    print_log_msg('Input volume will be aligned to Gaussian blob of radius %s' % \
                        options.volradius, log)
                    print_log_msg('Verbosity level to set to: %s' % options.verbosity, log)
                    print_log_msg('Debugging diagnostics set to: %s\n' % options.debug, log)
                
                center_blob(
                    input_obj, 
                    volradius=options.volradius, 
                    outdir=outdir, 
                    verbosity=options.verbosity, 
                    debug=options.debug, 
                    log=log
                    )
        
        else:
            msg= "Mode %s not supported. Allowed 3D modes are: acf, boundingbox" % mode
            sp_global_def.ERROR(msg, __file__, 1)
            sp_global_def.print_timestamp( "Finish" )
            exit()
    elif options.dimensions == 2:
        if options.verbosity >= 1 : print_log_msg("Reading '%s' as a stack\n" % inputfile, log)
    
        if mode == 'acf':
            if options.verbosity >=1 :
                print_log_msg("Binarizing images in stack '%s'" % inputfile, log)
                print_log_msg('Will binarize at %s sigma above mean' % options.threshsigma, log)
                print_log_msg('Number of dilations set to: %s' % options.ndilation, log)
                print_log_msg('Number of erosions set to: %s' % options.nerosion, log)
                print_log_msg('Verbosity level to set to: %s' % options.verbosity, log)
                print_log_msg('Debugging diagnostics set to: %s\n' % options.debug, log)
            
            mask_stack(
                inputfile, 
                outdir=outdir, 
                threshold_sd=options.threshsigma, 
                ndilation=options.ndilation, 
                nerosion=options.nerosion, 
                verbosity=options.verbosity, 
                log=log
                )
        
            if options.verbosity >=1 :
                print_log_msg('Centering using autocorrelation method', log)
                print_log_msg("Shifts will be applied to '%s'" % inputfile, log)
                print_log_msg('Will pad by a factor of %s before computing autocorrelation' % options.pad, log)
                print_log_msg('Intermediate binary files will be thresholded at %s\n' % options.threshint, log)
            
            center_stack_acf(
                os.path.join(outdir, MASK_STACK), 
                outdir=outdir, 
                imgs2shift=inputfile, 
                pad=options.pad, 
                threshold=options.threshint, 
                verbosity=options.verbosity, 
                debug=options.debug, 
                log=log
                )
            
            # Message to the user on what to look for
            if options.verbosity >=1 :
                msg=  "DONE!\n\n"
                mask_montage= os.path.join(outdir, MASK_MONTAGE)
                msg+= "HINT:\nCheck mask thresholding in output '%s'\n\n" % (mask_montage)
                shifted_avgs= os.path.join(outdir, AVG_SHIFT_BY_ACF)
                centered_masks= os.path.join(outdir, MASK2D_CENT)
                msg+= "*If* the masking is reasonable, then the subsequent centering of \n'%s' and '%s' \nshould also be reasonable.\n" % (shifted_avgs, centered_masks)
                print_log_msg(msg, log)
        
        elif mode == 'rt180':
            if options.verbosity >=1 :
                print_log_msg('Centering by alignment to 180-degree-rotated version', log)
                print_log_msg('Verbosity level to set to: %s' % options.verbosity, log)
                print_log_msg('Debugging diagnostics set to: %s\n' % options.debug, log)
            
            center_rt180(
                inputfile, 
                outdir=outdir, 
                verbosity=options.verbosity, 
                debug=options.debug, 
                log=log
                )
            
        else:
            msg= "Mode %s not supported. Allowed 2D modes are: acf, rt180" % mode
            sp_global_def.ERROR(msg, __file__, 1)
            sp_global_def.print_timestamp( "Finish" )
            exit()
    else:
        msg= "Dimension %s not supported. Allowed modes are: 3 (volume), 2 (image)" % options.dimensions
        sp_global_def.ERROR(msg, __file__, 1)
        sp_global_def.print_timestamp( "Finish" )
        exit()
        
def prepare_outdir_log(outdir='.', verbose=False, is_main=True):
    """
    Prepares output directory and sets up log file.
    
    Arguments:
        outdir : Output directory
        verbose : (boolean) Whether to write to screen
        is_main : (boolean) If using multiple cores, some tasks only need to be performed once, not by all cores
    Returns:
        log : instance of Logger class
        verbose : (boolean) New version of Logger can write to screen simultaneously
    """
    
    # Create directory if it doesn't exist
    if is_main:
        if os.path.isdir(outdir):
            sp_global_def.sxprint("Writing to output directory: %s" % outdir)
        else:
            sp_global_def.sxprint("Created output directory: %s" % outdir)
            os.makedirs(outdir)  # os.mkdir() can only operate one directory deep
        sp_global_def.write_command(outdir)

    logname= "log_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S") +  ".txt"
    logname= os.path.join(outdir, logname)
    
    # May be using old version of logger.py
    try:
        if verbose:
            log= sp_logger.Logger(base_logger=sp_logger.BaseLogger_Files(), base_logger2=sp_logger.BaseLogger_Print(), file_name=logname)
            verbose= False  # logger output will be echoed to screen
        else:
            log= sp_logger.Logger(base_logger=sp_logger.BaseLogger_Files(), file_name=logname)
    except TypeError:
        if is_main: sp_global_def.sxprint("WARNING: Using old sp_logger.py library")
        log= sp_logger.Logger(base_logger=sp_logger.BaseLogger_Files())#, file_name=logname)
        logname= 'log.txt'
        
    if is_main: sp_global_def.sxprint("Writing log file to %s\n" % logname)
    
    if is_main:
        progbase= os.path.basename(__file__).split('.')[0].upper()
        length= len(progbase) + 4
        
        log.add("\n" +
                " "*TIMESTAMP_LENGTH + "*"*length + "\n" +
                " "*TIMESTAMP_LENGTH + "* " + progbase + " *\n" +
                " "*TIMESTAMP_LENGTH + "*"*length)
    
    return log, verbose

def print_log_msg(msg, log=None, verbose=False, is_main=True):
    """
    Prints messages to log file and, optionally, to the screen.
    
    Arguments:
        msg : Message to write
        log : Instance of Logger class
        verbose : (boolean) Whether to write to screen
        is_main : (boolean) If using MPI, some tasks only need to be performed once, not by all cores
    """
    
    if is_main:
        if verbose: sp_global_def.sxprint(msg)
        if log: log.add(msg)

def center_patterson_3d(mask_obj, outdir='.', vol2shift=None, threshold=0.1, pad=2, stem='cent', verbosity=1, debug=False, log=None):
    """
    Centers a volume by finding the longest dimension using a Patterson map. 
    The long dimension will then be rotated to the z-axis. 
    The rotated volume is projected onto xy, and the longest dimension is then determined in xy.
    The projection is rotated, and the boundary in the last dimension will be determined.
    I've tested it only on binary volumes thus far.
    
    Arguments:
        mask_obj : input EMData object
        outdir : output directory
        vol2shift : volume to which to apply shift
        threshold : binarization threshold (0..1, default 0.1)
        pad : padding factor, to prevent wraparound artifacts in ACF
        stem : stem for output volumes
        verbosity : verbosity level (0..3)
        debug : writes intermediate diagnostic files
        log : instance of Logger class
    Returns:
        transform_instance : transform object of rotation + shift parameters
        mask_rotshift : input volume after rotationg and translation
        mask_shiftonly : input volume after x, y, and z shifts only (no rotation)
        shift_list : 3D shifts (no rotation) applied to mask_shiftonly
    """
    
    # Get volume dimension (assuming cube)
    idim= mask_obj['nx']

    # Binarize object (in case input was adaptive)
    temp_vol= sp_morphology.binarize(mask_obj, 0.9)
    
    # Pad volume (or else voxels near box edge may find voxels in next unit cell over)
    temp_vol= padvolby2(temp_vol, factor=pad)
    if debug : 
        filename= os.path.join(outdir, MASK3D_PADDED)
        temp_vol.write_image(filename)

    # Compute autocorrelation
    if verbosity >= 1 : print_log_msg('Computing autocorrelation map', log)
    patterson_3d= sp_fundamentals.scf(temp_vol)
    
    if verbosity >= 1 : 
        print_log_msg('Finished computing autocorrelation map', log)
        print_log_msg('Binarizing and finding longest dimension', log)

    # Binarize
    patterson_3d= sp_morphology.binarize(patterson_3d, minval=threshold)
    if debug : 
        filename= os.path.join(outdir, MASK3D_PATTERSON)
        patterson_3d.write_image(filename)

    # Find maximum dimension
    longest3d_dim, longest3d_vector, _= findMaxRadius3dFast(patterson_3d)
    if verbosity >= 2 : 
        msg= '  Longest dimension %.1f, longest vector %s' % (longest3d_dim, longest3d_vector)
        print_log_msg(msg, log)
    del patterson_3d

    # Rotate longest dimension to z-axis
    angle1= float( np.degrees( np.arctan2(longest3d_vector[1], longest3d_vector[0]) ) )
    angle2= float( np.degrees( np.arctan2(longest3d_vector[2], longest3d_vector[0]) ) )
    temp_vol= sp_fundamentals.rot_shift3D(mask_obj, phi=angle1, theta=angle2)
    if verbosity >= 1 : 
        print_log_msg( '  Rotated by phi=%.1f, theta=%.1f' % (angle1, angle2), log)
    if debug : 
        filename= os.path.join(outdir, MASK3D_ROTATED)
        temp_vol.write_image(filename)

    # Binarize
    temp_vol= sp_morphology.binarize(temp_vol, minval=threshold)

    # Find first & last nonzero voxels
    firstz, lastz= findextrema3d(temp_vol)
    centerz= old_div( (firstz+lastz), 2)
    shiftz= float( old_div(-(firstz+lastz-idim), 2) + 1 )

    # Shift midpoint along z to center
    temp_vol= sp_fundamentals.rot_shift3D(temp_vol, sz=shiftz)
    if verbosity >= 1 : 
        print_log_msg('  Center at z=%s, shifted by (0,0,%s)' % (centerz, shiftz), log)
    if debug : 
        filename= os.path.join(outdir, MASK3D_ZSHIFT)
        temp_vol.write_image(filename)

    # Project along 0,0,0
    proj_xy= sp_projection.project(temp_vol, [0,0,0, 0,0])

    # Binarize
    proj_xy= sp_morphology.binarize(proj_xy, minval=threshold)
    if debug : 
        filename= os.path.join(outdir, PROJ2D_XY)
        proj_xy.write_image(filename)
    
    if debug:
        anglexy, shiftx, shifty, _= center_patterson_2d(
            proj_xy,
            pad_factor=pad,
            bin_thresh=threshold, 
            verbosity=verbosity, 
            filenumber=None, 
            log=log, 
            outdir=outdir,
            padded_fn=PROJ2D_PAD, 
            acf_fn=PROJ2D_ACF, 
            acf_thresh_fn=PROJ2D_ACF_THRESH, 
            rot_fn=PROJ2D_ROT, 
            rot_thresh_fn=PROJ2D_ROT_THRESH, 
            shift_rot_fn=PROJ2D_ROTSHIFT, 
            shift_fn=PROJ2D_SHIFTONLY
            )
    else:
        anglexy, shiftx, shifty, _= center_patterson_2d(
            proj_xy,
            pad_factor=pad,
            bin_thresh=threshold, 
            verbosity=verbosity, 
            filenumber=None, 
            log=log, 
            )
    
    # Combined rotation + shift
    transform_instance= EMAN2.Transform({'type': 'SPIDER', 'phi': angle1, 'theta': angle2, 'psi': anglexy, 
                    'tx': shiftx, 'ty': shifty, 'tz': shiftz})
    mask_rotshift= mask_obj.rot_scale_trans_background(transform_instance, None)
    filename= os.path.join(outdir, 'mask_' + stem + '_rotshift.hdf')
    mask_rotshift.write_image(filename)
    if verbosity >= 1 : 
        msg= "  Wrote rotated + shifted volume to '%s'" % (filename)
        print_log_msg(msg, log)

    # "Undo" rotation, and only consider shifts
    phi, theta, psi, onlyx, onlyy, onlyz, scale= sp_utilities.compose_transform3(
        angle1, angle2, anglexy, shiftx, shifty, shiftz, 1, -anglexy, -angle2, -angle1, 0,0,0,1)
    shift_list= [onlyx, onlyy, onlyz]
    mask_shiftonly= sp_fundamentals.rot_shift3D(mask_obj, sx=shift_list[0], sy=shift_list[1], sz=shift_list[2])
    filename= os.path.join(outdir, 'mask_' + stem + '_shiftonly.hdf')
    mask_shiftonly.write_image(filename)
    if verbosity >= 1 : 
        msg= "  Shifted by (%.1f,%.1f,%.1f) and wrote to '%s'" % \
            (shift_list[0], shift_list[1], shift_list[2], filename)
        print_log_msg(msg, log)

    # Optionally repeat shift on a different volume
    if vol2shift:
        filename= os.path.join(outdir, 'vol_' + stem + '_shifted.hdf')
        if verbosity >= 1 : 
            msg= "  Also shifting %s by (%.1f,%.1f,%.1f) and writing to '%s'" % \
                (vol2shift, shift_list[0], shift_list[1], shift_list[2], filename)
            print_log_msg(msg, log)
            msg= "  e2proc3d.py %s %s --trans='%.1f,%.1f,%.1f'" % \
                (vol2shift, filename, shift_list[0], shift_list[1], shift_list[2])
            print_log_msg(msg, log)
            print()

        vol_shifted= sp_fundamentals.rot_shift3D(EMAN2.EMData(vol2shift), sx=shift_list[0], sy=shift_list[1], sz=shift_list[2])
        vol_shifted.write_image(filename)
    
    return transform_instance, mask_shiftonly, mask_rotshift, shift_list

def center_patterson_2d(input_img_obj, pad_factor=2, bin_thresh=0.1, 
                        verbosity=1, filenumber=None, log=None, outdir='.',
                        padded_fn=None, acf_fn=None, acf_thresh_fn=None, 
                        rot_fn=None, rot_thresh_fn=None, shift_rot_fn=None, 
                        shift_fn=None):

    """
    Centers a binary image by finding zero columns & rows around it, and centering the remaining box.
    The longest dimension is determined by using a Pattern map (ACF), and rotates the long dimension to the x-axis.
    
    Arguments:
        input_img_obj : input EMData object
        pad_factor : padding factor, to prevent wraparound artifacts in ACF
        bin_thresh : binarization threshold, for when interpolation gives non-binary values 
        verbosity : verbosity level (0..3)
        filenumber : prepends counter in prints statements, perhaps useful within a loop
        log : instance of Logger class
        outdir : output directory
        padded_fn : filename for option padded image
        acf_fn : filename for optional ACF
        acf_thresh_fn : filename for optional thresholded ACF
        rot_fn : filename for optional rotated image
        rot_thresh_fn : filename for optional rotated, thresholded image
        shift_rot_fn : filename for optional shifted + rotated image
        shift_fn : filename for optional shifted-only image
    Returns:
        transform_instance : transform object of rotation + shift parameters
        vol_shiftonly : input volume after x, y, and z shifts only (no rotation)
        vol_rotshift : input volume after rotationg and translation
        anglexy : rotation angle
        shiftx : x-shift
        shifty : y-shift
        shifted_img : shifted image (no rotation)
    """
    
    idim= input_img_obj['nx']
    
    # Pad
    intermed_img= padimgby2(input_img_obj, factor=pad_factor, verbosity=verbosity, log=log)
    if padded_fn : 
        filename= os.path.join(outdir, padded_fn)
        intermed_img.write_image(filename)

    # Compute autocorrelation function
    patterson_2d= sp_fundamentals.scf(intermed_img)
    if acf_fn :
        filename= os.path.join(outdir, acf_fn)
        patterson_2d.write_image(filename)

    # Binarize
    patterson_2d= sp_morphology.adaptive_mask_scipy(patterson_2d, threshold=bin_thresh, allow_disconnected=False, edge_width=0)
    if acf_thresh_fn :
        filename= os.path.join(outdir, acf_thresh_fn)
        patterson_2d.write_image(filename)

    # Find maximum dimension
    longest2d_dim, longest2d_vector, _= findMaxRadius2d(patterson_2d)
    del patterson_2d
    if verbosity >= 3 : 
        msg= 'Longest dimension %.1f, longest vector %s' % (longest2d_dim, longest2d_vector)
        if filenumber != None: msg= 'Image %s: ' % filenumber + msg
        print_log_msg(msg, log)

    # Rotate longest dimension to x-axis
    anglexy= float( np.degrees( np.arctan2(longest2d_vector[1], longest2d_vector[0]) ) )
    intermed_img= sp_fundamentals.rot_shift2D(input_img_obj, angle=anglexy, 
                                mode='background', interpolation_method='linear')
    if verbosity >= 2 : 
        msg= 'Rotated by psi=%.1f' % anglexy
        if filenumber != None : msg= 'Image %s: ' % filenumber + msg
        print_log_msg(msg, log)
    if rot_fn : 
        filename= os.path.join(outdir, rot_fn)
        intermed_img.write_image(filename)

    # Binarize (TODO: there were wraparound artifacts using rot_shift2d, maybe pad it?)
    intermed_img= sp_morphology.adaptive_mask_scipy(intermed_img, threshold=bin_thresh, allow_disconnected=False, edge_width=0)
    if rot_thresh_fn : 
        filename= os.path.join(outdir, rot_thresh_fn)
        intermed_img.write_image(filename)

    # Find first & last nonzero pixels along x
    firstx, lastx= findextrema2d(intermed_img, coord='x')
    if verbosity >= 3 : 
        msg= 'Nonzero from x {%s..%s}' % (firstx, lastx)
        if filenumber != None : msg= 'Image %s: ' % filenumber + msg
        print_log_msg(msg, log)
    centerx= old_div( (firstx+lastx), 2 )
    shiftx= float( old_div(-(firstx+lastx-idim), 2) )

    # Find first & last nonzero pixels along y
    firsty, lasty= findextrema2d(intermed_img, coord='y')
    if verbosity >= 3 : 
        msg= 'Nonzero from y {%s..%s}' % (firsty, lasty)
        if filenumber != None : msg= 'Image %s: ' % filenumber + msg
        print_log_msg(msg, log)
    centery= old_div( (firsty+lasty), 2 )
    shifty= float( old_div(-(firsty+lasty-idim), 2) )

    # Apply rotation + shifts to original image
    intermed_img= sp_fundamentals.rot_shift2D(input_img_obj, angle=anglexy, sx=shiftx, sy=shifty, 
                                mode='cyclic', interpolation_method='linear')
    if shift_rot_fn : 
        filename= os.path.join(outdir, shift_rot_fn)
        intermed_img.write_image(filename)
    del intermed_img
    
    # "Undo" rotation, and only consider shifts
    alpha, onlyx, onlyy, scale= sp_utilities.compose_transform2(anglexy, shiftx, shifty, 1, -anglexy, 0, 0, 1)
    shifted_img= sp_fundamentals.rot_shift2D(input_img_obj, sx=onlyx, sy=onlyy, 
                                mode='cyclic', interpolation_method='linear')
    if shift_fn : 
        filename= os.path.join(outdir, shift_fn)
        shifted_img.write_image(filename)
    if verbosity >= 2 :
        msg= 'Center at (%s,%s), shifted by (%s,%s)\n' % (centerx,centery, shiftx,shifty)
        if filenumber != None : msg= 'Image #%s: ' % filenumber + msg
        print_log_msg(msg, log)
        
    return anglexy, shiftx, shifty, shifted_img

def center_stack_acf(mask2center, outdir='.', imgs2shift=None, pad=2, 
                    threshold=0.1, verbosity=1, debug=False, log=None):

    """
    Centers an image by finding the longest dimension using a Patterson map. 
    The long dimension will then be rotated to the x-axis. 
    The boundary in the last dimension will be determined.
    
    Arguments:
        mask2center : stack of images (e.g., binary masks) to center
        outdir : output directory
        imgs2shift : stack of images (e.g., class averages) to which to apply shifts
        pad : padding factor, to prevent wraparound artifacts in ACF
        threshold : binarization threshold (0..1, default 0.1)
        verbosity : verbosity level (0..3)
        debug : writes intermediate diagnostic files
        log : instance of Logger class
    """
    # Read stack as list
    mask_list= read_stack(mask2center)
    
    if imgs2shift:
        avg_list= read_stack(imgs2shift)
        
        if len(mask_list) != len(avg_list):
            msg= "Inputs %s and %s have different numbers of images, %s and %s" % ( mask2center, imgs2shift, len(mask_list), len(avg_list) )
            sp_global_def.ERROR(msg, __file__, 1)
            sp_global_def.print_timestamp( "Finish" )
            exit()
        else: print_log_msg("Will apply shifts to second stack file '%s'" % imgs2shift, log)
            
    if verbosity >= 1 : print_log_msg('Centering images', log)
    
    centered_mask_file= os.path.join(outdir, MASK2D_CENT)
    centered_avg_file= os.path.join(outdir, AVG_SHIFT_BY_ACF)
    shift_params_file= os.path.join(outdir, SHIFTONLY_DOC)
    
    # Loop through images
    for idx, mask_orig in enumerate(mask_list):
        if debug and idx == 0:
            anglexy, shiftx, shifty, mask_aligned= center_patterson_2d(
                mask_orig,
                pad_factor=pad,
                bin_thresh=threshold, 
                verbosity=verbosity, 
                filenumber=idx, 
                log=log, 
                outdir=outdir,
                padded_fn=MASK2D_PAD, 
                acf_fn=MASK2D_ACF, 
                acf_thresh_fn=MASK2D_ACF_THRESH, 
                rot_fn=MASK2D_ROT, 
                rot_thresh_fn=MASK2D_ROT_THRESH, 
                shift_rot_fn=MASK2D_ROTSHIFT, 
                shift_fn=MASK2D_SHIFTONLY
                )
        else:
            anglexy, shiftx, shifty, mask_aligned= center_patterson_2d(
                mask_orig,
                pad_factor=pad,
                bin_thresh=threshold, 
                verbosity=verbosity, 
                filenumber=idx, 
                log=log, 
                outdir=outdir
                )
        
        # Write alignment parameters to header
        sp_utilities.set_params2D(mask_aligned, [anglexy, shiftx, shifty, 0, 1.0], xform='xform.align2d')
        mask_aligned.write_image(centered_mask_file, idx)
        
        if imgs2shift:
            # Reported shifts presume rotation first, so first we need to "undo" the rotation
            alpha, onlyx, onlyy, scale= sp_utilities.compose_transform2(anglexy, shiftx, shifty, 1, -anglexy, 0, 0, 1)
            
            # Apply shift
            avg_shifted= sp_fundamentals.rot_shift2D(avg_list[idx], sx=onlyx, sy=onlyy, 
                                mode='cyclic', interpolation_method='linear')
            
            # Write alignment parameters to header to be used by other programs
            sp_utilities.set_params2D(avg_shifted, [alpha, onlyx, onlyy, 0, scale], xform=HEADER_LABEL)
            
            # Also write to xform.align2d so that sp_header writes in a recognizable format
            sp_utilities.set_params2D(avg_shifted, [alpha, onlyx, onlyy, 0, scale], xform='xform.align2d')
            
            # Write to disk
            avg_shifted.write_image(centered_avg_file, idx)
    # End image-loop

    if verbosity >= 1 : 
        msg= "Wrote %s centered-mask images to '%s'" % (len(mask_list), centered_mask_file)
        print_log_msg(msg, log)
        
    if imgs2shift:
        # Write alignment parameters to text file
        sp_applications.header(centered_avg_file, 'xform.align2d', fexport=shift_params_file)
        
        if verbosity >= 1 : 
            msg= "Wrote %s shifted images to '%s'" % (len(avg_list), centered_avg_file)
            print_log_msg(msg, log)
            
            print_log_msg("Exported shift parameters from '%s' into '%s'" % (centered_avg_file, shift_params_file), log)
            cmd= "  sp_header.py %s --params='%s' --export=%s\n" % (centered_avg_file, HEADER_LABEL, shift_params_file)
            print_log_msg(cmd, log)
            
def center_boundingbox(mask_obj, outdir='.', vol2shift=None, threshold=0.1, verbosity=1, debug=False, log=None):
    """
    Centers a binary volume by finding zero slices around it, and centering the remaining box.
    In comparison to the Patterson-map method above, this method doesn't look for the longest dimension. 
    For irregularly-shaped objects, the result will be somewhat different.
    
    Arguments:
        mask_obj : input EMData object
        outdir : output directory
        vol2shift : volume to which to apply shift
        threshold : binarization threshold (0..1, default 0.1)
        verbosity : verbosity level (0..3)
        debug : writes intermediate diagnostic files
        log : instance of Logger class
    Returns:
        transform_instance : transform object of rotation + shift parameters
        mask_shifted : shifted volume
        shift_list : list of applied 3D shifts
    """

    # Get volume dimension (assuming cube)
    idim= mask_obj['nx']

    # Binarize object (in case input was adaptive)
    temp_vol= sp_morphology.binarize(mask_obj, 0.9)
    
    # Find first & last nonzero voxels
    firstz, lastz= findextrema3d(temp_vol)
    centerz= old_div( (firstz+lastz), 2 )
    shiftz= float( old_div(-(firstz+lastz-idim), 2) + 1 )

    # Shift midpoint along z to center
    temp_vol= sp_fundamentals.rot_shift3D(temp_vol, sz=shiftz)
    if verbosity >= 1 :
        print_log_msg('Center at z=%s, shifted by (0,0,%s)' % (centerz, shiftz), log)
    if debug : 
        filename= os.path.join(outdir, MASK3D_ZSHIFT)
        temp_vol.write_image(filename)

    # Project along 0,0,0
    proj_xy= sp_projection.project(temp_vol, [0,0,0, 0,0])

    # Binarize
    proj_xy= sp_morphology.binarize(proj_xy, minval=threshold)
    if debug : 
        filename= os.path.join(outdir, PROJ2D_XY)
        proj_xy.write_image(filename)

    # Find first & last nonzero voxels along x
    firstx, lastx= findextrema2d(proj_xy, coord='x')
    if verbosity >= 2 : 
        print_log_msg('Nonzero from x {%s..%s}' % (firstx, lastx), log)
    centerx= old_div( (firstx+lastx), 2)
    shiftx= float( old_div(-(firstx+lastx-idim), 2) + 1 )

    # Find first & last nonzero voxels along y
    firsty, lasty= findextrema2d(proj_xy, coord='y')
    if verbosity >= 2 : 
        print_log_msg('Nonzero from y {%s..%s}' % (firsty, lasty), log)
    centery= old_div( (firsty+lasty), 2 )
    shifty= float( old_div(-(firsty+lasty-idim), 2) + 1 )

    proj_xy= sp_fundamentals.rot_shift2D(proj_xy, sx=shiftx, sy=shifty, 
                                mode='cyclic', interpolation_method='linear')
    if debug : 
        filename= os.path.join(outdir, 'PROJ2D_SHIFTONLY')
        proj_xy.write_image(filename)
    if verbosity >= 1 : 
        msg= 'Center at (%s,%s), shifted by (%s,%s,0)' % (centerx,centery, shiftx,shifty)
        print_log_msg(msg, log)

    # Combined shift
    transform_instance= EMAN2.Transform({'type': 'SPIDER', 'tx': shiftx, 'ty': shifty, 'tz': shiftz})
    shift_list= [shiftx, shifty, shiftz]
    mask_shifted= temp_vol.rot_scale_trans_background(transform_instance, None)
    filename= os.path.join(outdir, MASK3D_CENT)
    mask_shifted.write_image(filename)
    if verbosity >= 1 : 
        msg= 'Wrote rotated + shifted volume to %s' % (filename)
        print_log_msg(msg, log)

    # Optionally repeat shift on a different volume
    if vol2shift:
        filename= os.path.join(outdir, VOL_CENT_BY_BOX)
        if verbosity >= 1 : 
            msg= 'Also shifting %s by (%.1f,%.1f,%.1f) and writing to %s' % \
                (vol2shift, shift_list[0], shift_list[1], shift_list[2], filename)
            print_log_msg(msg, log)
            msg= "  e2proc3d.py %s %s --trans='%.1f,%.1f,%.1f'" % \
                (vol2shift, filename, shift_list[0], shift_list[1], shift_list[2])
            print_log_msg(msg, log)

        vol_shifted= sp_fundamentals.rot_shift3D(EMAN2.EMData(vol2shift), sx=shift_list[0], sy=shift_list[1], sz=shift_list[2])
        vol_shifted.write_image(filename)
    
    return transform_instance, mask_shifted, shift_list

def center_blob(input_vol_obj, volradius, outdir='.', verbosity=1, debug=False, log=None):
    """
    Centers a volume by aligning it to a Gaussian blob. 
    
    Arguments:
        input_vol_obj : input stack, EMData object or filename
        volradius : (float) radius of Gaussian blob to align reconstruction, pixels
        outdir : output directory
        verbosity : verbosity level (0..3)
        debug : writes intermediate diagnostic files
        log : instance of Logger class
    """
    
    # Make blob
    if debug:
        blobfile= os.path.join(outdir, BLOB_CENT)
        if verbosity >=1 : print_log_msg('Writing blob file %s' % blobfile, log)
    else:
        blobfile= None
    blob_obj= make_blob(input_vol_obj, volradius, outfile=blobfile)  # gets dimensions from input volume
    
    # Align experimental volume to blob
    bestshift= align_vols(input_vol_obj, blob_obj, numpeaks=2, log=log, verbose=verbosity>=2)
    
    # Apply shifts
    input_vol_obj.transform(bestshift)
    filename= os.path.join(outdir, VOL_CENT_BY_BLOB)
    input_vol_obj.write_image(filename)
    xshift= bestshift.get_params('spider').get('tx')
    yshift= bestshift.get_params('spider').get('ty')
    zshift= bestshift.get_params('spider').get('tz')
    if verbosity >=1 : 
        msg= "Shifted by (%s, %s, %s) to %s\n" % (xshift, yshift, zshift, filename)
        print_log_msg(msg, log)
    
    return

def make_blob(dimensions, radiuspx, outfile=None):
    """
    Generates a Gaussian blob.
    Adapted from sp_signalsubtract.py, 2020-03-22
    Input parameter 'dimensions' can be an integer, 2020-03-22
    
    Arguments:
        dimensions : list of 1 or 3 elements of volume dimension, or EMData object from which dimensions will be copied
        radiuspx : (float) radius of Gaussian blob to align reconstruction, pixels
        outfile : output volume of Gaussian blob (optional)
    Returns:
        blob : EMData object of Gaussian blob
    """
    
    # Dimensions can be either from a volume or a list of one element [i] or three elements [x,y,z]
    if isinstance(dimensions, list):
        if len(dimensions) == 3:
            xdim= dimensions[0]
            ydim= dimensions[1]
            zdim= dimensions[2]
        elif len(dimensions) == 1:
            xdim= dimensions[0]
            ydim= dimensions[0]
            zdim= dimensions[0]
        else:
            sp_global_def.ERROR("One or three dimensions expected", __file__, 1)
            return
    elif isinstance(dimensions, EMAN2.EMData):
        xdim= dimensions['nx']
        ydim= dimensions['ny']
        zdim= dimensions['nz']
    elif isinstance(dimensions, int):
        xdim= dimensions
        ydim= dimensions
        zdim= dimensions
    else:
        sp_global_def.ERROR("Expecting either a volume or a list of length 3 or 1", __file__, 1)
        return
        
    blob= EMAN2.EMData()
    blob.set_size(xdim, ydim, zdim)
    
    xcenter= old_div(xdim, 2)
    ycenter= old_div(ydim, 2)
    zcenter= old_div(zdim, 2)
    blob.process_inplace('testimage.puregaussian', {'x_center':xcenter, 'x_sigma':radiuspx, 'y_center':ycenter,
                                                'y_sigma':radiuspx, 'z_center':zcenter, 'z_sigma':radiuspx })
    
    if outfile: blob.write_image(outfile)
    
    return blob

def align_vols(expvol, refvol, numpeaks=1, log=None, verbose=False, is_main=True):
    """
    Translationally aligns two volumes.
    Adapted from sp_signalsubtract.py, 2020-03-22
    
    Arguments:
        expvol : input experimental volume to be aligned
        refvol : input reference volume to be aligned to
        numpeaks : (int) number of cross-correlation peaks desired
        log : instance of Logger class
        verbose : (boolean) whether to write to screen
        is_main : (boolean) if using multiple cores, some tasks only need to be performed once, not by all cores
    Returns:
        shift : instance of Transform class containing shifts necessary for alignment
    """
    
    ccmap= sp_fundamentals.ccf(expvol, refvol)
    peaks= sp_utilities.peak_search(ccmap, numpeaks)  # e.g., [109.94314575195312, 187.0, 185.0, 192.0, 1.0, -3.0, -5.0, 2.0]
    
    if verbose:
        msg= "   PEAK  XPEAK  YPEAK  ZPEAK  RPEAK  XORIG  YORIG  ZORIG"
        print_log_msg(msg, log)
        
        for line in peaks: print_log_msg("".join(format(item, "7.1f") for item in line), log)

    xshift= -peaks[0][5]
    yshift= -peaks[0][6]
    zshift= -peaks[0][7]
    shift= EMAN2.Transform({"type":"spider","phi":0,"theta":0,"psi":0,"tx":xshift,"ty":yshift,"tz":zshift})
    
    return shift

def center_rt180(input, outdir='.', verbosity=1, debug=False, log=None):
    """
    Centers an image by rotating an image by 180 degrees and aligning with the original. 
    
    Arguments:
        input : input stack, EMData object or filename
        outdir : output directory
        verbosity : verbosity level (0..3)
        debug : writes intermediate diagnostic files
        log : instance of Logger class
    """
    
    # Read stack as list
    image_list= read_stack(input)
    
    if verbosity >= 1 : print_log_msg('Centering images', log)
    
    centered_file= os.path.join(outdir, AVG_SHIFT_BY_RT180)
    shift_params_file= os.path.join(outdir, SHIFTONLY_DOC)
    
    num_peaks= 2
    
    # Loop through images
    for img_idx, img_orig in enumerate(image_list):
        # Rotate image by 180 degrees
        rot_img= sp_fundamentals.rot_shift2D(img_orig, angle=180, 
                                mode='background', interpolation_method='linear')
        if debug and img_idx == 0:
            filename= os.path.join(outdir, RT180_ROT)
            rot_img.write_image(filename)
            
        # Compute cross-correlation
        ccmap= sp_fundamentals.ccf(rot_img, img_orig)
        if debug and img_idx == 0:
            filename= os.path.join(outdir, RT180_CCF)
            ccmap.write_image(filename)

        # Find peak
        peak_list= sp_utilities.peak_search(ccmap, npeak=num_peaks)
        
        # Print peak list
        if verbosity >= 3 : 
            headers= " Img Peak#  Value    X     Y   Ratio  Xctr  Yctr"
            print_log_msg(headers, log)
            
            for peak_idx, peak in enumerate(peak_list):
                peakdata= "%3i    %1i   %7.3f %5.1f %5.1f %5.3f %5.1f %5.1f" % \
                    (img_idx, peak_idx, peak[0], peak[1], peak[2], peak[3], peak[4], peak[5])
                print_log_msg(peakdata, log)
        
        # Shift by half the distance
        xshift= old_div(peak_list[0][4], 2)
        yshift= old_div(peak_list[0][5], 2)
        centered_img= sp_fundamentals.rot_shift2D(img_orig, sx=xshift, sy=yshift, 
                                mode='background', interpolation_method='linear')

        if verbosity >= 2 : 
            msg= "Shifting image #%s by (%s,%s)" % (img_idx, xshift, yshift)
            print_log_msg(msg, log)
        
        # Write alignment parameters to header
        sp_utilities.set_params2D(centered_img, [0.0, xshift, yshift, 0, 1.0], xform=HEADER_LABEL)
        
        # Write centered file
        centered_img.write_image(centered_file, img_idx)
    # End image-loop

    # Write alignment parameters to text file
    sp_applications.header(centered_file, HEADER_LABEL, fexport=shift_params_file)
    
    if verbosity >= 1 : 
        msg= "Wrote %s centered images to '%s'" % (len(image_list), centered_file)
        print_log_msg(msg, log)
        
        print_log_msg("Exported shift parameters from '%s' into '%s'" % (centered_file, shift_params_file), log)
        cmd= "sp_header.py %s --params='%s' --export=%s\n" % (centered_file, HEADER_LABEL, shift_params_file)
        print_log_msg(cmd, log)
        
def padvolby2(input_vol_obj, factor=2):
    """
    Pads volume. 
    Original volume will be centered in the padded one.
    
    Arguments:
        input_vol_obj : input EMData object
        factor : (int) multiplication factor
    Returns:
        padvol : padded volume
    """
    
    # Assume cube
    idim= input_vol_obj['nx']

    # Pad by specified padding factor
    paddim= idim*factor

    # Pad (will be centered in new box)
    padvol= sp_utilities.pad(input_vol_obj, paddim, paddim, paddim)

    return padvol

def padimgby2(input_img_obj, factor=2, verbosity=1, log=None):
    """
    Pads image. 
    Original image will be centered in the padded one.
    
    Arguments:
        input_img_obj : input EMData object
        factor : (int) multiplication factor
    Returns:
        padimg : padded image
    """
    
    # Assume square
    idim= input_img_obj['nx']
    if verbosity >= 3 : print_log_msg('Starting dimension: %s' % idim, log)

    # Pad by specified padding factor
    paddim= idim*factor
    if verbosity >= 3 : print_log_msg('Computed padded dimension: %s' % paddim, log)

    # Pad (will be centered in new box)
    padimg= sp_utilities.pad(input_img_obj, paddim, paddim, background='circumference')
    if verbosity >= 3 : print_log_msg('Confirmed dimension: %s' % padimg['nx'], log)

    return padimg

def findMaxRadius3dFast(input_vol_obj, threshold=0):
    """
    Finds voxel above a given threshold with the highest distance from the volume center. 
    
    Arguments:
        input_vol_obj : input EMData object
        threshold : only voxels above this threshold will be considered
    Returns:
        sqrt_maximum_value : longest distance vector, voxels
        vector_centered : position of most distant voxel from the center
        num_possible_indices : longest vector should come in pairs (positive and negative), but there could be more
    """
    
    # Convert the volume to a NumPy array
    data_numpy= input_vol_obj.get_3dview()

    # Get volume dimension (assuming square)
    nx= input_vol_obj['nx']

    # Make an array of the Pythagorean radii for each position in the matrix
    indices= ((np.mgrid[0:nx]-nx//2)**2).astype(np.dtype('int32'))
    squared_distance= np.add.outer(np.add.outer(indices, indices), indices)

    # Make a boolean mask such that values above the threshold are true
    mask= data_numpy > threshold

    # Find the maximum of the masked distance matrix
    maximum_value= np.max(squared_distance[mask])
    sqrt_maximum_value= np.sqrt(maximum_value)

    possible_indices= list( zip( *np.where( np.logical_and((squared_distance == maximum_value), mask) ) ) )
    num_possible_indices= len(possible_indices)
    # There should be two solutions, the second being the negative of the first
    # (I don't understand the notation. Markus wrote it.)

    # Distance from center
    vector_centered= [x - nx//2 for x in possible_indices[0]][::-1]
    # (The notation [::-1] reorders the axes.)

    return sqrt_maximum_value, vector_centered, num_possible_indices

def findMaxRadius2d(inputimgobj, threshold=0):
    """
    Finds pixel above a given threshold with the highest distance from the volume center. 
    
    Arguments:
        inputimgobj : input EMData object
        threshold : only pixels above this threshold will be considered
    Returns:
        sqrt_maximum_value : longest distance vector, pixels
        vector_centered : position of most distant pixel from the center
        num_possible_indices : longest vector should come in pairs (positive and negative), but there could be more
    """
    
    # Convert the volume to a NumPy array
    data_numpy= inputimgobj.get_2dview()

    # Get volume dimension (assuming square)
    nx= inputimgobj['nx']

    # Make an array of the Pythagorean radii for each position in the matrix
    indices= ((np.mgrid[0:nx]-nx//2)**2).astype(np.dtype('int32'))
    squared_distance= np.add.outer(indices, indices)

    # Make a boolean mask such that values above the threshold are true
    mask= data_numpy > threshold

    # Find the maximum of the masked distance matrix
    maximum_value= np.max(squared_distance[mask])
    sqrt_maximum_value= np.sqrt(maximum_value)

    possible_indices_3d= zip(*np.where(np.logical_and((squared_distance == maximum_value), mask)))
    # There should be two solutions, the second being the negative of the first
    # (I don't understand the notation. Markus wrote the 3d version.)

    # Convert the 3D coordinates back to 2D (first coordinate will be 0)
    possible_indices= [t for t in possible_indices_3d]
    num_possible_indices= len(possible_indices)

    # Distance from center
    vector_centered= [x - nx//2 for x in possible_indices[0]][::-1]
    # The notation [::-1] reorders the axes

    return sqrt_maximum_value, vector_centered, num_possible_indices

def findextrema3d(input_vol_obj):
    """
    Finds first and last nonzero z-slices in a volume. 
    Reshapes 3D array into a 2D array, and then looks for nonzero rows.
    
    Arguments:
        input_vol_obj : input EMData object
    Returns:
        first_index : top nonzero z-slice
        last_index : bottom nonzero z-slice
    """
    
    # Convert the volume to a NumPy array
    data_numpy= input_vol_obj.get_3dview()

    # Flatten 3D array to 2D
    threed_to_twod= data_numpy.reshape(data_numpy.shape[0], -1)

    # Find nonzero vectors (reshaped 2D arrays)
    nonzero_boolean= np.any(threed_to_twod, axis=1)

    # Find indices of nonzero slices
    nonzero_indices= np.where(nonzero_boolean)[0]
    first_index= nonzero_indices[0]  + 1
    last_index = nonzero_indices[-1] + 1

    # Return first & last nonzero slice numbers
    return first_index, last_index

def findextrema2d(inputimgobj, coord='x'):
    """
    Finds first and last nonzero columns in an image. 
    Reshapes 2D array into a 1D array, and then looks for nonzero elements.
    
    Arguments:
        inputimgobj : input EMData object
    Returns:
        first_index : top nonzero column
        last_index : bottom nonzero column
    """
    
    # Convert the image to a NumPy array
    data_numpy= inputimgobj.get_2dview()

    if coord == 'x':
        axis= 0
    elif coord == 'y':
        axis= 1
    else:
        sp_global_def.ERROR("Coordinate axis %s not supported" % coord, __file__, 1)
        sp_global_def.print_timestamp( "Finish" )
        exit()

    # Find nonzero columns
    nonzero_boolean= np.any(data_numpy, axis=axis)

    # Find indices of nonzero slices
    nonzero_indices= np.where(nonzero_boolean)[0]
    first_index= nonzero_indices[0]
    last_index = nonzero_indices[-1]

    # Return first & last nonzero slice numbers
    return first_index, last_index

def mask_stack(input, outdir='.', threshold_sd=0, ndilation=3, nerosion=0, verbosity=1, log=None):
    """
    Masks a stack of images.
    
        input : input stack, EMData object or filename
        outdir : output directory
        threshold_sd : binarization threshold, in units of sigma above the mean
        ndilation : number of dilations
        nerosion : number of erosions
        verbosity : verbosity level (0..3)
        log : instance of Logger class
    """

    # Read stack as list
    image_list= read_stack(input)
    
    if verbosity >= 1 : 
        msg= 'Using a threshold of %s sigma above the mean, with %s dilations and %s erosions' % (
            threshold_sd, ndilation, nerosion
            )
        print_log_msg(msg, log)
    
    mask_hdf_stack= os.path.join(outdir, MASK_STACK)
    masked_stack= os.path.join(outdir, MASK_APPLIED)
    comparison_stack= os.path.join(outdir, MASK_MONTAGE)
    
    # Loop through images
    for idx, img in enumerate(image_list):
        # Compute threshold as a function of the avg and sigma
        threshold= img['mean'] + threshold_sd*img['sigma']
        
        if verbosity >= 2 : 
            print_log_msg("Image #%s: mean %.4f, sigma %.4f, threshold %.4f" % \
                (idx, img['mean'], img['sigma'], threshold), log=log )
        
        # Make binary mask
        mask_img= sp_morphology.adaptive_mask_scipy(
            img, 
            threshold=threshold, 
            allow_disconnected=True, 
            do_fill=True, 
            edge_width=0, 
            ndilation=ndilation,
            nerosion=nerosion
            )
        
        # The option allow_disconnected seems to be applied before dilation, so I'm going to do it last, separately.
        mask_img= sp_morphology.adaptive_mask_scipy(
            mask_img, 
            threshold=0.5, 
            allow_disconnected=False, 
            edge_width=0, 
            ndilation=0,
            nerosion=0
            )
        
        # Write mask and masked images
        mask_img.write_image(mask_hdf_stack, idx)
        masked_img= mask_img*img
        masked_img.write_image(masked_stack, idx)
        
        # Make side-by-side montage
        montagestack= []
        montagestack.append(img)
        montagestack.append(mask_img)
        montagestack.append(masked_img)
        comparison_pair= mn_s(
            montagestack, 
            ncol=3,
            marginwidth=2, 
            bkgd=0.25, 
            scale=True, 
            debuginfo=idx, 
            log=log, 
            verbose=verbosity>=1
            )
        comparison_pair.write_image(comparison_stack, idx)
    # End image-loop
    
    if verbosity >= 1 : 
        print_log_msg("Wrote mask stack to '%s'" % mask_hdf_stack, log)
        print_log_msg("Wrote stack of original images multipled by mask to '%s'" % masked_stack, log)
        print_log_msg("Wrote stack of montages of original, mask, and masked images to '%s'\n" % comparison_stack, log)
    
def read_stack(input):
    """
    Read stack. 
    Can be EMData object or filename.
    
    Arguments:
        input : input EMData object or filename
    Returns:
        image_list : list of EMData objects
    """
    
    if isinstance(input, list):
        # Do nothing, is already a list
        image_list= input
    # The following should work for a filename and for MRC stacks
    elif isinstance(input, str):
        # Read from file
        image_list= EMAN2.EMData.read_images(input)
    else:
        msg= "Input %s is not a recognized type %s" % ( input, type(input) )
        sp_global_def.ERROR(msg, __file__, 1)
        sp_global_def.print_timestamp( "Finish" )
        exit()
    
    return image_list

def mn_s(inputstack, ncol, marginwidth=0, bkgd=0, outfile=None, scale=False, debuginfo=None, log=None, verbose=False):
    """
    Generates montage of images into one image.
    Adapted from sp_projcompare.montage2, 2020-03-20
    Scale flag scales each image from 0..1, similar to SPIDER's MN S operation
    
    Arguments:
        inputstack : Stack of input images to merge into montage
        ncol : Number of images per row
        marginwidth : Margin width, pixels
        bkgd : Background value of montage
        outfile : Optional output file with montage output
        scale : Flag to scale each image's intensity from 0..1 (cf. MN S from SPIDER)
        debuginfo : Information to print in case of error
        log : instance of Logger class
        verbose : Flag to print messages to screen
    Returns:
        montage : EMData object of image montage
    """
    
    if isinstance(inputstack, str): inputstack= EMAN2.EMData.read_images(inputstack)
    
    # Get single-image dimensions
    nx= inputstack[0].get_xsize()
    ny= inputstack[0].get_ysize()
    
    # Get number of images and calculate montage dimensions
    numimgs= len(inputstack)
    numrows= old_div( (numimgs-1), ncol) + 1
    
    # Create blank image
    montage_xdim= (nx + marginwidth)*ncol
    montage_ydim= (ny + marginwidth)*numrows
    montage= sp_utilities.model_blank(montage_xdim, montage_ydim, 1, bkgd)
    
    # Loop through images
    for imgnum in range(numimgs):
        # Horizontal grid position is image# modulo NCOL
        colnum= imgnum % ncol
        
        # Montage is numbered from the top down
        rownum= numrows - 1 - old_div(imgnum, ncol)
        
        xoffset= colnum*(nx+marginwidth)
        yoffset= rownum*(ny+marginwidth)
        img= inputstack[imgnum]
        
        # Scale, if requested
        if scale:
            # Trap for blank images
            if img['maximum'] != img['minimum']: 
                img= old_div( (img-img['minimum']), (img['maximum']-img['minimum']) )
            else:
                if verbose and debuginfo: print_log_msg('Empty image found at %s' % debuginfo, log)
                img= sp_utilities.model_blank(nx, ny, 1, bkgd)
        
        insert_image(img, montage, xoffset, yoffset)
    
    if outfile: montage.write_image(outfile)
    
    return montage
    
def insert_image(smallimg, largeimage, xoffset, yoffset):
    """
    Inserts small image into large image.
    Copied from sp_projcompare.py, 2020-03-20
    
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
        'inputfile', 
        type=str, 
        help='Input file to center (3D volume or 2D stack)')
    
    parser.add_argument(
        'outdir', 
        type=str, 
        help='Output directory')
    
    parser.add_argument(
        '--mode', 
        type=str, 
        default='acf', 
        help='Centering options: acf, boundingbox, blob, rt180')
    
    parser.add_argument(
        '--applyshift', 
        type=str, 
        help='File to which to apply shifts, if mask used as input')
    
    parser.add_argument(
        '--dimensions', 
        type=int, 
        default=3, 
        help='Number of dimensions for input(s)')
    
    parser.add_argument(
        '--threshint', 
        type=float, 
        default=0.1, 
        help='Binarization threshold for intermediate masks')
    
    parser.add_argument(
        '--volradius', 
        type=int, 
        help='Volume radius, pixels, used only in 3D blob mode')
    
    parser.add_argument(
        '--stem', 
        type=str, 
        default='acf', 
        help='File stem for output files, used only in 3D ACF mode')
    
    parser.add_argument(
        '--verbosity', '-v', 
        type=int,
        default=1, 
        help='Verbosity level (0..3)')

    parser.add_argument(
        '--debug', 
        action="store_true", 
        help='Write additional files')

    group_acf= parser.add_argument_group(
        title='ACF options',
        description='Options using autocorrelation function (Patterson map), 2D and 3D.')
    
    group_acf.add_argument(
        '--threshsigma', 
        type=float, 
        default=1., 
        help='Binarization threshold for initial masks, units of sigma above mean')
    
    group_acf.add_argument(
        '--ndilation', '-nd', 
        type=int, 
        default=3, 
        help='Number of iterations to dilate the binarized volume. Each iteration adds about two pixels to the mask.')
    
    group_acf.add_argument(
        '--nerosion', '-ne', 
        type=int, 
        default=0, 
        help='Number of iterations to erode the binarized volume. Each iteration removes about two pixels from the mask.')
    
    group_acf.add_argument(
        '--pad', 
        type=int, 
        default=2, 
        help='Padding factor')
    
    return parser.parse_args()
    
def main():
    sp_global_def.print_timestamp( "Start" )
    options= parse_command_line()

    ##print args, options  # (Everything is in options.)
    #print('options', options)
    #print('LOGFILE', sp_global_def.LOGFILE)
    #exit()

    center_2d3d(
        options.inputfile, 
        options.mode, 
        options.outdir, 
        options=options
        )

    sp_global_def.print_timestamp( "Finish" )

if __name__ == "__main__":
    main()
