#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
from past.utils import old_div

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
import os
import math
import datetime
import shutil
import glob
import EMAN2
from sphire.libpy import sp_applications
from sphire.libpy import sp_utilities
from sphire.libpy import sp_statistics
from sphire.libpy import sp_fundamentals
from sphire.libpy import sp_projection
from sphire.libpy import sp_filter
from sphire.libpy import sp_logger
from sphire.libpy import sp_alignment
from sphire.libpy import sp_global_def
from sphire.libpy import sp_pixel_error
import numpy as np
import argparse
import tqdm
import inspect
import six
import EMAN2_cppwrap

# Set default values and filenames (global variables written in ALL CAPS)
TIMESTAMP_LENGTH= 23  # chars
REFPROJS= 'refproj.hdf'
REFANGLES= 'refangles.txt'
OUTALIGN= 'docalign2d.txt'
PROJPARAMS= 'docprojparams.txt'
CLASSDIR= 'Byclass'
CLASSMAP= 'classmap.txt'
CLASSDOCSTEM= 'docclass'
MASIMG= 'maskalign.hdf'
COMPARISON_STACK= 'stkcompare.hdf'
GOODCLASSPARTS= 'goodpartsclass{{0:0{0}d}}.txt'  # Meridien mode w/outliers

USAGE= """
PURPOSE:
Compare projections of a 3D model to a stack of 2D projections (averages or raw images).

To use angles stored in the image header:
%s <input_imgs> <input_volume>
General, required command-line parameters:
1. Input image stack
2. Input volume to be projected, same dimension as images
Outputs:
docprojection.txt : Projection angles applied to input model
stkcompare.hdf : Stack of reprojections and averages 

General options:
%s <input_imgs> <input_volume> <output_directory> --projparams <angles_file> --select <img_selection_file> --projmethod <interpolation_method> --display
Parameters:
--projmethod : Interpolation method : trilinear (default), gridding, nn (nearest neighbor)
--display : Automatically open montage in e2display

To use angles from a VIPER or RVIPER run:
%s <input_imgs> <input_volume> <output_directory> --mode viper --projparams <angles_file> --projselect <img_selection_file>
Parameters:
--projparams : Projection angles (For RVIPER, this file will have a name like main003/run002/rotated_reduced_params.txt)
--projselect : Image selection file (For RVIPER, not all images will necessarily be assigned alignment angles, in which case the number of angles in the doc file above will differ from the number of images in the input stack. This file will have a name like main003/index_keep_image.txt)

To apply VIPER angles to the member particles' headers:
%s <input_imgs> <input_volume> <output_directory> --mode viper --projparams <angles_file> --projselect <img_selection_file> --fullstack <bdb_stack> 
Parameters:
--fullstack : Input particle stack to which VIPER angles will be added

To perform an internal round of projection matching against the input model:
%s <input_imgs> <input_volume> <output_directory> --mode projmatch --delta <angular_increment> --matchshift <shift_range> --matchrad <outer_radius> --matchstep <ring_step> --symmetry <optional_symmetry>
Parameters:
--delta : Angular-sampling increment
--matchshift : Maximum shift to allow during translation alignment (default 0)
--matchrad : Outer alignment radius (default Automatic)
--matchstep : Alignment radius step size (default 1)
--symmetry : To limit angular projections (default c1)

To use the average projection angles from MERIDIEN refinement:
%s <input_imgs> <input_volume> <output_directory> --mode meridien --refineparams <refinement_params> --partselect <substack_select> --refineshift <shift_range> --outliers <max_angle_diff> --refinerad <outer_radius> --refinestep <ring_step> --alignmethod <2d_alignment_method>
Parameters:
--refineparams : Input refinement parameters
--partselect : Input substack selection file if particles removed before refinement (e.g., Substack/isac_substack_particle_id_list.txt)
--outliers : Particles differing from average Euler angle by more than this threshold (in degrees) will be excluded from average calculation (default keep all)
--refinerad : Outer alignment radius (default Automatic)
--refineshift : Maximum shift to allow during translation alignment (default 2)
--refinestep : Alignment radius step size (default 1)
--alignmethod : Alignment method: apsh (default) or scf

""" % ((__file__,)*6)
    
MODIFIED="Modified 2021-02-19"

"""
Modifications log:
    2021-02-19 (trs) -- updated for SPHIRE 1.4
    2020-06-22 (trs) -- writes command.txt
    2020-06-18 (trs) -- fixed normalization
    2020-04-22 (trs & sp) -- can be used on single images
    2020-04-09 (trs & ab) -- added class membership information to comparison stack
"""


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
        'classavgs', 
        type=str, 
        help='Input class averages')
    
    parser.add_argument(
        'vol3d', 
        type=str, 
        help='Input 3D reconstruction')
    
    parser.add_argument(
        'outdir', 
        type=str, 
        default='.', 
        help='Output directory (default .)')
    
    parser.add_argument(
        '--mode', 
        type=str, 
        default='viper', 
        help='Options: viper (default), projmatch, meridien')
    
    parser.add_argument(
        '--projmethod', 
        type=str, 
        default='trilinear', 
        help='Projection method: trilinear (default), gridding, nn (nearest neighbor)')
    
    parser.add_argument(
        '--display', 
        action="store_true", 
        help='Automatically open montage in e2display')
    
    parser.add_argument(
        '--verbosity', "-v", 
        type=int, 
        default=2,
        help='Increase verbosity (0..2)')
    
    group_viper= parser.add_argument_group(
        title='VIPER options (--mode=viper)',
        description='Projection angles for input images will be extracted from initial-reconstruction results.')
    
    group_viper.add_argument(
        '--projparams', 
        type=str, 
        help='Mode viper, angles file, which will be imported into the header of the class-average stack')
    
    group_viper.add_argument(
        '--projselect', 
        type=str, 
        help='Selection file for included classes.')
    
    group_viper.add_argument(
        '--fullstack', 
        type=str, 
        help='Input particle stack to which VIPER angles will be added.')
    
    group_projmatch= parser.add_argument_group(
        title='Projection-matching options (--mode=projmatch)',
        description='Input image will be aligned to re-projections of the input volume.')
    
    group_projmatch.add_argument(
        '--delta', 
        type=float, 
        default=7.5, 
        help='Angular search increment for reference projections (degrees, default 7.5)')
    
    group_projmatch.add_argument(
        '--symmetry', 
        type=str, 
        default='c1', 
        help='Symmetry (default c1)')
    
    group_projmatch.add_argument(
        '--matchshift', 
        type=int, 
        default=2, 
        help='Maximum translational shift during multireference alignment (pixels, default 2)')
    
    group_projmatch.add_argument(
        '--matchrad', 
        type=int, 
        default=-1, 
        help='Outer radius for alignment (pixels, default automatic)')
    
    group_projmatch.add_argument(
        '--matchstep', 
        type=int, 
        default=1, 
        help='Radius increment during alignment (pixels, default 1)')
    
    group_meridien= parser.add_argument_group(
        title='MERIDIEN options (--mode=meridien)',
        description='Projection angles from existing refinement will be averages. Outliers can optionally be omitted.')
    
    group_meridien.add_argument(
        '--refineparams', 
        type=str, 
        help='Mode meridien, alignment parameters for experimental images')
    
    group_meridien.add_argument(
        '--partselect', 
        type=str, 
        help='Selection file for included particles.')
    
    group_meridien.add_argument(
        '--classdocs', 
        type=str, 
        help='Mode meridien, file pattern for class-to-particle lists')
    
    group_meridien.add_argument(
        '--outliers', 
        type=float, 
        help='Mode meridien, maximum angular difference from initial estimate')
    
    group_meridien.add_argument(
        '--refinerad', 
        type=int, 
        default=-1, 
        help='Outer radius for alignment (pixels, default automatic)')
    
    group_meridien.add_argument(
        '--refineshift', 
        type=int, 
        default=2, 
        help='Maximum translational shift during multireference alignment (pixels, default 2)')
    
    group_meridien.add_argument(
        '--refinestep', 
        type=int, 
        default=1, 
        help='Radius increment during alignment (pixels, default 1)')
    
    group_meridien.add_argument(
        '--alignmethod', 
        type=str, 
        help='Mode meridien, 2D alignment method: apsh (default) or scf')
    
    return parser.parse_args()

def compare_projections(classavgstack, reconfile, outdir, 
            mode='viper', projmethod='trilinear', 
            projparamfile=None, refineparams=None, selectdoc=None, 
            displayYN=False, verbosity=0, options=None):
    """
    Main function overseeing various projection-comparison modes.
    
    Arguments:
        classavgstack : Input image stack
        reconfile : Map of which to generate projections (an optionally perform alignment)
        outdir : Output directory
        mode : Mode, viper (pre-existing angles for each input image), projmatch (angles from internal projection-matching)
        ####verbose : (boolean) Whether to write additional information to screen
        projparamfile : Angles and shifts for each input class average
        refineparams : Angles and shifts for each particle (mode meridien)
        selectdoc : Selection file for included images
        projmethod : Interpolation method to use
        displayYN : (boolean) Whether to automatically open montage
        verbosity : Verbosity level (0..2)
        options : (list) Command-line options
    """
    
    # Check if inputs exist
    check(classavgstack, verbose=verbosity>=1)
    check(reconfile, verbose=verbosity>=1)
    
    # Set output directory and log file name
    prepare_outdir(outdir, verbose=verbosity>=1)
    sp_global_def.write_command(outdir)
    log, _= prepare_log(outdir)

    # Expand path for outputs
    refprojstack= os.path.join(outdir, REFPROJS)
    refanglesdoc= os.path.join(outdir, REFANGLES)
    out2daligndoc= os.path.join(outdir, OUTALIGN)
    goodclassparttemplate= None  # will be overridden in Meridien/outliers mode

    # If not an input, will create an output, in modes projmatch
    exportAnglesTF= False
    if projparamfile == None:
        exportAnglesTF= True
        projparamfile= os.path.join(outdir, PROJPARAMS)
        
        # You need either input angles (mode viper) or to calculate them on the fly (mode projmatch)
        if mode=='viper':
            assert sp_utilities.get_im(classavgstack,0).has_attr('xform.projection'), "ERROR! Projection angles provided neither in header nor file"
            if verbosity>=1: 
                print_log_msg('Parameter document not provided, will extract from header...', log)
    
    # Check that dimensions of images and volume agree (maybe rescale volume)
    voldim= EMAN2.EMData(reconfile).get_xsize()
    imgdim= EMAN2.EMData(classavgstack,0).get_xsize()
    if voldim != imgdim:
            sp_global_def.ERROR("\nERROR!! Dimension of input volume doesn't match that of image stack: %s vs. %s" % 
                    (voldim, imgdim), __file__, 1)
            
            scale= old_div(float(imgdim), voldim)  # only approximate, since full-sized particle radius is arbitrary
            msg = 'The command to resize the volume will be of the form:\n'
            msg += 'e2proc3d.py %s resized_vol.hdf --scale=%1.5f --clip=%s,%s,%s\n' % (
                reconfile, scale, imgdim, imgdim, imgdim
                )
            msg += 'Check the file in the ISAC directory named "README_shrink_ratio.txt" for confirmation.\n'
            print_log_msg(msg, log)
            exit()
    
    # Check whether dimension is odd (prep_vol and/or prgl does something weird if so. --Tapu)
    if voldim % 2 != 0:
        voldim= voldim + 1
        if verbosity>=1: 
            print_log_msg("WARNING! Inputs have odd dimension, will pad to %s" % voldim, log)
        
    #  Here if you want to be fancy, there should be an option to chose the projection method,
    #  the mechanism can be copied from sxproject3d.py  PAP
    if projmethod=='trilinear':
        method_num= 1
    elif projmethod=='gridding':
        method_num= -1
    elif projmethod=='nn':
        method_num= 0
    else:
        sp_global_def.ERROR("\nERROR!! Valid projection methods are: trilinear (default), gridding, and nn (nearest neighbor).", __file__, 1)
        exit()
    
    if verbosity >=1 : print_log_msg('Using mode %s\n' % mode, log)
    
    # In case class averages include discarded images, apply selection file
    if mode == 'viper':
        if selectdoc:
            goodavgs, extension= os.path.splitext(os.path.basename(classavgstack))
            newclasses= os.path.join(outdir, goodavgs + "_kept" + extension)
            
            # e2proc2d appends to existing files, so rename existing output
            if os.path.exists(newclasses):
                renamefile= newclasses + '.bak'
                if verbosity >=1 : 
                    mesg= "Selected-classes stack %s exists, renaming to %s" % (newclasses, renamefile)
                    print_log_msg(mesg, log)
                    print_log_msg("mv %s %s\n" % (newclasses, renamefile), log)
                os.rename(newclasses, renamefile)
            
            cmd= "e2proc2d.py %s %s --list=%s" % (classavgstack, newclasses, selectdoc)
            if verbosity >=1 : 
                print_log_msg('Creating subset of %s to %s based on selection list %s' % (
                    classavgstack, newclasses, selectdoc
                    ), log)
                print_log_msg(cmd + '\n', log)
            os.system(cmd)
            
            # Update class-averages
            classavgstack= newclasses
        
        # Create parameter doc if needed
        if exportAnglesTF:
            projparamfile= os.path.join(outdir, PROJPARAMS)
            if verbosity >=1 : 
                mesg= "Exporting parameter information from %s into %s" % (classavgstack, projparamfile)
                print_log_msg(mesg, log)
                cmd= "sp_header.py %s --params=xform.projection --export=%s\n" % (
                    classavgstack, projparamfile
                    )
                print_log_msg(cmd, log)
            sp_applications.header(classavgstack, 'xform.projection', fexport=projparamfile)
            
        
        # (TODO: Move to function)
        if options.fullstack:
            if verbosity >=1 : print_log_msg('Will copy angles to %s' % options.fullstack, log)
    
            # Make sure enough digits
            num_imgs  = EMAN2.EMUtil.get_image_count(classavgstack)
            num_digits= math.ceil(math.log10(num_imgs) )
            
            # Filenames
            classdir= os.path.join(outdir, CLASSDIR)
            if not os.path.isdir(classdir): 
                os.makedirs(classdir)
            classmap= os.path.join(classdir, CLASSMAP)
            classdoc= os.path.join(classdir, CLASSDOCSTEM + '{{0:0{0}d}}.txt'.format(num_digits) )
            classdocs= os.path.join(classdir, CLASSDOCSTEM + '*.txt')
            
            # Separate particles by class
            vomq(classavgstack, classmap, classdoc, log=log, verbose=verbosity>=1)
            
            # Read class lists
            classdoclist= glob.glob(classdocs)
            
            # Read Euler angles
            angleslist= sp_utilities.read_text_row(projparamfile)
            
            # Read stack
            if verbosity >=1 : print_log_msg("Reading %s" % options.fullstack, log)
            stack_obj= EMAN2.EMData.read_images(options.fullstack)
            num_tot_imgs= len(stack_obj)
            if verbosity>=1: 
                mesg= "Finished reading %s images from %s\n" % (num_tot_imgs, options.fullstack)
                print_log_msg(mesg, log)
            
            # Loop through class lists
            for classidx, classdoc in enumerate(classdoclist):
                # Strip out filenumber
                classexample= os.path.splitext(classdoc)
                classnum= int(classexample[0][-num_digits:])  # number of digits in class number
                if verbosity>=1: 
                    print_log_msg('Index %s, class number %s' % (classidx, classnum), log)
                
                # Read angles
                curr_angles= angleslist[classidx]
                
                # Read class list
                partlist= sp_utilities.read_text_row(classdoc)
                
                # Loop through particles
                for totpartnum in partlist:
                    img_obj= sp_utilities.get_im(stack_obj, totpartnum[0])  
                    # totpartnum is a list of length 1
                    
                    sp_utilities.set_params_proj(img_obj, curr_angles, xform="xform.projection")
                    
    # align class averages de novo to projections of a reference map
    if mode=='projmatch':
        projparamfile= os.path.join(outdir, PROJPARAMS)
        
        mode_projmatch(
            reconfile, 
            refprojstack, 
            classavgstack, 
            outdir, 
            refanglesdoc, 
            projparamfile, 
            out2daligndoc, 
            delta=options.delta, 
            symmetry=options.symmetry, 
            projmethod=options.projmethod, 
            align_method=options.alignmethod, 
            align_radius=options.matchrad, 
            align_shift=options.matchshift, 
            align_step=options.matchstep, 
            log=log, 
            verbosity=verbosity
            )
        
    # Get alignment parameters from MERIDIEN 
    if mode=='meridien':
        projparamfile= os.path.join(outdir, PROJPARAMS)
        
        continueTF= True  # Will proceed unless some information is missing
        
        if not refineparams:
            sp_global_def.ERROR("\nERROR!! Input alignment parameters not provided.", __file__, 1)
            continueTF= False
    
        if not continueTF:
            print_log_msg('Type %s --help to see available options\n' % os.path.basename(__file__), log)
            exit()
        
        if not options.classdocs or options.outliers:
            classdir= os.path.join(outdir, CLASSDIR)
            if not os.path.isdir(classdir): os.makedirs(classdir)
            
            # Make sure enough digits
            num_imgs  = EMAN2.EMUtil.get_image_count(classavgstack)
            num_digits= math.ceil(math.log10(num_imgs) )
            
            if options.outliers: 
                goodclassparttemplate= os.path.join(classdir, GOODCLASSPARTS.format(num_digits) )
        
            if not options.classdocs:
                classmap= os.path.join(classdir, CLASSMAP)
                classdoc= os.path.join(classdir, CLASSDOCSTEM + '{{0:0{0}d}}.txt'.format(num_digits) )
                options.classdocs= os.path.join(classdir, CLASSDOCSTEM + '*.txt')
                
                # Separate particles by class
                vomq(classavgstack, classmap, classdoc, log=log, verbose=verbosity>=1)
            
        mode_meridien(
            reconfile, 
            classavgstack, 
            options.classdocs, 
            refineparams, 
            selectdoc, 
            options.refineshift, 
            options.refinerad, 
            projparamfile, 
            out2daligndoc, 
            interpolation_method=method_num, 
            outliers=options.outliers,
            goodclassparttemplate=goodclassparttemplate, 
            alignopt=options.alignmethod, 
            ringstep=options.refinestep, 
            log=log, 
            verbosity=verbosity
            )
        
    # Import Euler angles (Don't need to import them if they were exported from the same place)
    if not exportAnglesTF:
        # Sanity check
        params_length= len( open(projparamfile).readlines() )
        num_imgs= EMAN2.EMUtil.get_image_count(classavgstack)
        if num_imgs != params_length:
            msg= "ERROR!! Number of images in '%s' (%s) does not equal number of lines in the parameter file '%s' (%s)" % (
                classavgstack, num_imgs, projparamfile, params_length
                )
            msg+= "\nYou probably need to supply a class-selection file (flag '--projselect') from the RVIPER directory, with a name like 'index_keep_images.txt'\n"
            print(msg)
            exit()
        
        if verbosity >=1 : 
            msg= "Importing parameter information into %s from %s" % (classavgstack, projparamfile)
            print_log_msg(msg, log)
            cmd= "sp_header.py %s --params=xform.projection --import=%s\n" % (
                classavgstack, projparamfile
                )
            print_log_msg(cmd, log)
        sp_applications.header(classavgstack, 'xform.projection', fimport=projparamfile)
    
    # Make comparison stack between class averages (images 0,2,4,...) and re-projections (images 1,3,5,...)
    comparison_stack= compare_projs(
        reconfile, 
        classavgstack, 
        outdir, 
        inputprojparams=projparamfile, 
        interpolation_method=method_num, 
        classdoctemplate=goodclassparttemplate, 
        log=log, 
        verbosity=verbosity
        )

    # Optionally pop up e2display
    if displayYN:
        cmd= "e2display.py %s" % comparison_stack
        if verbosity >=1 : 
            print_log_msg('Opening montage', log)
            print_log_msg(cmd + '\n', log)
        os.system(cmd)
    
def check(file, verbose=True):
    """
    Checks whether file exists.
    
    Arguments:
        file : File to look for
        verbose : (boolean) Whether to write to screen
    """
    
    if os.path.exists(file):
        if verbose: sp_global_def.sxprint("Found %s" % file)
    else:
        sp_global_def.sxprint("ERROR!! %s doesn't exist!\n" % file)
        exit()

def prepare_outdir(outdir='.', verbose=False, main=True):
    """
    Create directory if it doesn't exist.
    
    Arguments:
        outdir : Output directory
        verbose : (boolean) Whether to write additional information to screen
        main : (boolean) Whether main MPI process
    """
    
    # If using MPI, only need to check once
    if main:
        if os.path.isdir(outdir):
            if verbose: print("Writing to output directory: %s" % outdir)
        else:
            if verbose: print("Created output directory: %s" % outdir)
            os.makedirs(outdir)  # os.mkdir() can only operate one directory deep
            
def prepare_log(outdir='.', verbose=False, main=True):
    """
    Prepare log file.
    
    Arguments:
        outdir : Output directory
        verbose : (boolean) Whether to write additional information to screen
        main : (boolean) Whether main MPI process
    Returns:
        log : Instance of Logger class
        verbose : (boolean) Updates flag to False if Logger class can mirror output to screen
    """
    
    logname= "log_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S") +  ".txt"
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
        print("WARNING: Using old sp_logger.py library")
        log= sp_logger.Logger(base_logger=sp_logger.BaseLogger_Files())#, file_name=logname)
        logname= 'log.txt'
        
    print("Writing log file to %s" % logname)
    
    if main:
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
    
    ####printvars([
        ####'is_main',
        ####'verbose',
        ####'log',
        ####'msg',
        ####], quitTF=True, typeTF=True)
    
    if is_main:
        if verbose: print(msg)
        if log: log.add(msg)

def printvars(variables, quitTF=False, typeTF=False):
    """Print the local variables in the caller's frame.
    
    Adapted from https://stackoverflow.com/questions/6618795/get-locals-from-calling-namespace-in-python
    """
    
    if type(variables) is list:
        # Weird things happen if 
        assert isinstance(variables[0], six.string_types), "UH OH!! Passed non-string %s instead of variable name" % variables[0]
        
        variable_list= variables
    elif isinstance(variables, six.string_types):  # six works with both Pythons 2 & 3
        variable_list= [variables]
    else:
        print("ERROR!! Don't know how to deal with type %s" % type(variables) )
        exit()
    
    frame= inspect.currentframe()
    dictionary= frame.f_back.f_locals
    
    print("")
    for variable in variable_list : 
        try: 
            msg= "%s: %s" % (variable, dictionary[variable])
            if typeTF: msg+= " %s" % type(dictionary[variable])
            print(msg)
        except KeyError: 
            print('ERROR!! printvars')
            print(dictionary)
            print('EXITING!!')
            exit()
        
    del frame
    
    if quitTF: 
        print('\nExiting printvars...')  # reminder in case I forget to take out the quit flag
        exit()

def mode_projmatch(reconfile, refprojstack, classavgstack, outdir, refanglesdoc, projparamfile, 
                   out2daligndoc, delta=7.5, symmetry='c1', projmethod='trilinear', 
                   align_method='apsh', align_radius=-1, align_shift=2, align_step=1, 
                   log=None, verbosity=0):
    """
    Makes projections of a map and runs multireference alignment between projections and images.
    
    Arguments:
        reconfile : Input volume from which to generate re-projections
        refprojstack : Output stack of reference projections
        classavgstack ; Input image stack
        outdir ; Output directory
        refanglesdoc : Output reference-projection angles
        projparamfile : Output experimental-image projection parameters 
        out2daligndoc : Output 2D alignment parameters
        delta : anglar increment for reference projections
        symmetry : Symmetry of input reconstruction
        projmethod : Projection method (default trilinear)
        align_method : Alignment method (defaul apsh)
        align_radius : Alignment radius, pixels
        align_shift : Maximum translation tested during alignment
        align_step : Translation step size
        log : Logger object
        ####verbose : (boolean) Whether to write additional information to screen
        verbosity : Verbosity level (0..2)
    """
    
    # Check if volume needs to be padded
    vol_obj= EMAN2.EMData(reconfile)
    nx= vol_obj.get_xsize()
    if nx % 2 != 0:
        nx= nx +1
        
        # Volume
        if verbosity >=1 : print_log_msg("Padding volume and images to %s" % nx, log)
        vol_obj= sp_utilities.pad(vol_obj, nx, nx, nx, background='circumference')
        reconfile= os.path.join(outdir, 'padvol.hdf')
        vol_obj.write_image(reconfile)
        # (TODO: Padding is performed multiple times.)
        
        # Images
        padstack= os.path.join(outdir, 'padimgs.hdf')
        cmd= 'e2proc2d.py %s %s --clip=%s,%s' % (classavgstack, padstack, nx, nx)
        os.system(cmd)
        if verbosity >=1 : print_log_msg(cmd + '\n', log)
        
        # Update class-average stack name
        classavgstack= padstack
        
        # Sanity check
        img_obj= sp_utilities.get_im(classavgstack, 0)
        imgx= img_obj.get_xsize()
        assert imgx == nx, "MRK_DEBUG"
    
    # Generate reference projections
    if verbosity >=1 : 
        msg= 'Projecting %s to output %s using an increment of %s degrees using %s symmetry' % (
            reconfile, refprojstack, delta, symmetry
            )
        print_log_msg(msg, log)
        cmd= 'sp_project3d.py %s %s --delta=%s --method=S --phiEqpsi=Minus --symmetry=%s' % (
            reconfile, refprojstack, delta, symmetry
            )
        if projmethod == 'trilinear':
            cmd += ' --trilinear'
            #trilinear= True
        #else: 
            #trilinear= False
        cmd += '\n'
        print_log_msg(cmd, log)
    sp_applications.project3d(
        reconfile, 
        refprojstack, 
        delta=delta, 
        symmetry=symmetry, 
        trillinear=projmethod=='trilinear'
        )  # 'trillinear' is spelled wrong
    
    # Export projection angles
    if verbosity >=1 : 
        print_log_msg("Exporting projection angles from %s to %s" % (refprojstack, refanglesdoc), log)
        cmd= "sp_header.py %s --params=xform.projection --export=%s\n" % (refprojstack, refanglesdoc)
        print_log_msg(cmd, log)
    sp_applications.header(refprojstack, 'xform.projection', fexport=refanglesdoc)
    
    # Perform multi-reference alignment
    if align_method=='ali2d':
        projdir= os.path.join(outdir, 'Projdir')  # used if input angles no provided
        if os.path.isdir(projdir):
                if verbosity >=1 : 
                    print_log_msg('Removing pre-existing directory %s' % projdir, log)
                    print_log_msg('rm -r %s\n' % projdir, log)
                shutil.rmtree(projdir)  # os.rmdir only removes empty directories
        
        # Zero out alignment parameters in header
        if verbosity >=1 : 
            print_log_msg('Zeroing out alignment parameters in header of %s' % classavgstack, log)
            cmd= 'sp_header.py %s --params xform.align2d --zero\n' % classavgstack
            print_log_msg(cmd, log)
        sp_applications.header(classavgstack, 'xform.align2d', zero=True)
        
        # Perform multi-reference alignment
        if verbosity >=1 : 
            msg= 'Aligning images in %s to projections %s with a radius of %s and a maximum allowed shift of %s' % (classavgstack, refprojstack, align_radius, align_shift)
            print_log_msg(msg, log)
            cmd= 'sp_mref_ali2d.py %s %s %s --ou=%s --xr=%s --yr=%s\n' % (
                classavgstack, refprojstack, projdir, align_radius, align_shift, align_shift
                )
            print_log_msg(cmd, log)
        sp_applications.mref_ali2d(
            classavgstack, 
            refprojstack, 
            projdir, 
            ou=align_radius, 
            xrng=align_shift, 
            yrng=align_shift
            )
        
        # Export alignment parameters
        if verbosity >=1 : 
            print_log_msg('Exporting angles from %s into %s' % (classavgstack, projparamfile), log)
            cmd= "sp_header.py %s --params=xform.align2d --export=%s\n" % (classavgstack, projparamfile)
            print_log_msg(cmd, log)
        sp_applications.header(classavgstack, 'xform.align2d', fexport=projparamfile)	
    
    # By default, use AP SH
    else:
        apsh(
            refprojstack, 
            classavgstack, 
            outprojparamfile=projparamfile, 
            refanglesdoc=refanglesdoc, 
            out2daligndoc=out2daligndoc, 
            outerradius=align_radius, 
            maxshift=align_shift, 
            ringstep=align_step, 
            log=log, 
            verbosity=verbosity
            )
    
    # Diagnostic
    paramlist= sp_utilities.read_text_row(projparamfile)  # contain 2D alignment parameters
    nimg1  = EMAN2.EMUtil.get_image_count(classavgstack)
    assert len(paramlist) == nimg1, "UH OH!! Number of images in %s and %s are different!" % (
        projparamfile, classavgstack
        )
    
    
def apsh(refimgs, imgs2align, outprojparamfile=None, refanglesdoc=None, out2daligndoc=None, 
         outerradius=-1, maxshift=0, ringstep=1, mode="F", log=None, verbosity=0):
    """
    Generates polar representations of a series of images to be used as alignment references.
    
    Arguments:
        refimgs : Input reference image stack (filename or EMData object)
        imgs2align : Image stack to be aligned (filename or EMData object)
        outprojparamfile : Output Euler angles doc file
        refanglesdoc : Input Euler angles for reference projections
        out2daligndoc : Output 2D alignment doc file
        outerradius : Outer alignment radius
        maxshift : Maximum shift allowed
        ringstep : Alignment radius step size
        mode : Mode, full circle ("F") vs. half circle ("H")
        log : Logger object
        ####verbose : (boolean) Whether to write additional information to screen
        verbosity : Verbosity level (0..2)
    """
    
    # Generate polar representation(s) of reference(s)
    alignrings, polarrefs= mref2polar(refimgs, outerradius=outerradius, ringstep=ringstep, 
            log=log, verbose=verbosity>=1)
            
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
    halfdim= old_div(idim, 2) + 1
    
    # Set constants
    currshift= 0
    scale= 1
    
    # Initialize output angles
    outprojparamslist= []
    out2dalignlist= []
    
    if outerradius <= 0:
        outerradius= halfdim - 3
        
    if refanglesdoc:
        refangleslist= sp_utilities.read_text_row(refanglesdoc)
        
    # Set search range
    txrng= tyrng= sp_alignment.search_range(idim, outerradius, currshift, maxshift)
    
    if verbosity>=1: 
        print_log_msg('Running multireference alignment allowing a maximum shift of %s' % maxshift, log)
    
    disableTF= verbosity<1
    
    # Loop through images
    for imgindex in tqdm.tqdm(range(numimg), unit='img', disable=disableTF):
        currimg= imagestacklist[imgindex]
        
        # Perform multi-reference alignment (adapted from alignment.mref_ali2d)
        [angt, sxst, syst, mirrorfloat, bestreffloat, peakt]= EMAN2_cppwrap.Util.multiref_polar_ali_2d(
            currimg, polarrefs, txrng, tyrng, ringstep, mode, alignrings, halfdim, halfdim)
        bestref= int(bestreffloat)
        mirrorflag= int(mirrorfloat)
        
        # Store parameters
        params2dlist= [angt, sxst, syst, mirrorflag, scale]
        out2dalignlist.append(params2dlist)
    
        if refanglesdoc:
            besteulers= refangleslist[bestref]
        else:
            besteulers= [0]*5
        
        # Check for mirroring
        if mirrorflag == 1:
            tempeulers= list(
                sp_utilities.compose_transform3(
                    besteulers[0],besteulers[1],besteulers[2], besteulers[3],besteulers[4],0, 1,
                    0,            180,          0,             0,            0,            0, 1
                    )
                )
            combinedparams= list(
                sp_utilities.compose_transform3(
                    tempeulers[0],tempeulers[1],tempeulers[2], tempeulers[3],tempeulers[4],0, 1, 
                    0,            0,            -angt,         0,            0,            0, 1
                    )
                )
        else:
            combinedparams= list(
                sp_utilities.compose_transform3(
                    besteulers[0],besteulers[1],besteulers[2], besteulers[3],besteulers[4],0, 1, 
                    0,            0,            -angt,         0,            0,            0, 1
                    )
                )
        # compose_transform3: returns phi,theta,psi, tx,ty,tz, scale
        
        outprojparamslist.append(combinedparams)
        
        # Set transformations as image attribute
        sp_utilities.set_params2D(currimg, params2dlist, xform="xform.align2d")  # sometimes I get a vector error with sxheader
        sp_utilities.set_params_proj(currimg, besteulers, xform="xform.projection")  # use shifts
    # End image loop
        
    if verbosity >=1 : print_log_msg('Finished running multireference alignment', log)

    if outprojparamfile or out2daligndoc:
        if outprojparamfile : 
            sp_utilities.write_text_row(outprojparamslist, outprojparamfile)
            if verbosity >=1 : 
                msg= 'Wrote alignment angles to %s' % outprojparamfile
                print_log_msg(msg, log)
        if out2daligndoc : 
            sp_utilities.write_text_row(out2dalignlist, out2daligndoc)
            if verbosity >=1 : 
                msg= 'Wrote 2D alignment parameters to %s\n' % out2daligndoc
                print_log_msg(msg, log)
        
    return out2dalignlist
    
def mref2polar(refimgs, firstring=1, outerradius=-1, ringstep=1, mode="F", normbysquare=0, 
               log=None, verbose=False):
    """
    Generates polar representations of a series of images to be used as alignment references.
    
    Arguments:
        refimgs : Input reference image stack (filename or EMData object)
        firstring : Inner alignment radius
        outerradius : Outer alignment radius
        ringstep : Alignment radius step size
        mode : Mode, full circle ("F") vs. half circle ("H)
        normbysquare : If other than 0, normalization by setting the norm to 1
        log : Logger object
        verbose : (boolean) Whether to write additional information to screen
    Returns:
        alignringlist : List of alignment-ring data
        polarreflist : List of polar representation of refernences
    """
    
    # Read reference stack
    if isinstance(refimgs, str):
        referencelist= EMAN2.EMData.read_images(refimgs)
    else:
        referencelist= [refimgs]  # For single image
    
    numrefs= len(referencelist)
    polarreflist= []
    
    # Get image dimensions (assuming square, and that images and references have the same dimension)
    idim= referencelist[0]['nx']

    # Calculate image center
    halfdim= old_div(idim, 2) + 1
    
    if outerradius <= 0: outerradius= halfdim - 3
        
    # Prepare alignment rings
    alignringlist= sp_alignment.Numrinit(firstring, outerradius, ringstep, mode)
    
    # Calculate ring weights
    ringweightlist= sp_alignment.ringwe(alignringlist, mode)
    
    if verbose : 
        msg= 'Converting %s references to polar coordinates from radius %s to %s with step %s and mode "%s"' % (
            numrefs, firstring, outerradius, ringstep, mode)
        print_log_msg(msg, log)
    
    # Loop through reference images (adapted from sxisac2)
    for refindex in range(numrefs):
        # Convert to polar
        cimage= EMAN2_cppwrap.Util.Polar2Dm(referencelist[refindex] , halfdim, halfdim, alignringlist, mode)
        
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

def vomq(classavgstack, classmap, classdoc, log=None, verbose=False):
    """
    Separate particles according to class assignment.
    
    Arguments:
        classavgstack : Input image stack
        classmap : Output class-to-particle lookup table. Each (long) line contains particles assigned to a class, one file for all classes
        classdoc : Output lists of particles assigned to a class, one file per class
        mode : Mode, viper (pre-existing angles for each input image), projmatch (angles from internal projection-matching)
        log : instance of Logger class
        verbose : (boolean) Whether to write additional information to screen
    """
    
    # Generate class-to-particle lookup table
    if verbose : 
        mesg= "Exporting members of stack %s to class map %s" % (classavgstack, classmap)
        print_log_msg(mesg, log, verbose)
        cmd= "sp_header.py %s --params=members --export=%s" % (classavgstack, classmap) 
        print_log_msg(cmd, log)
    sp_applications.header(classavgstack, 'members', fexport=classmap)
    
    counter= 0
    
    # Loop through classes
    with open(classmap) as r:
        for idx, line in enumerate(r.readlines()):
            with open(classdoc.format(idx), 'w') as w:
                w.write('\n'.join(line[1:-3].split(', ')))
            counter += 1

    if verbose : print_log_msg("Wrote %s class selection files\n" % counter, log)

def mode_meridien(reconfile, classavgstack, classdocs, refineparams, selectdoc, maxshift, outerrad, 
                    outprojparamsfile, out2daligndoc, interpolation_method=1, outliers=None, 
                    goodclassparttemplate=None, alignopt='apsh', ringstep=1, 
                    log=None, verbosity=0):
    """
    Averages projection angles within a class, computes projection, and aligns to input images. 
    Optionally, an angular threshold can be applied, and outliers will be excluded.
    
    Arguments:
        reconfile : Input volume from which to generate re-projections
        classavgstack ; Input image stack
        classdocs : Input class-membership files
        refineparams : Projection parameters from refinement
        selectdoc : Selection file for included images
        maxshift : Maximum translation to be tested, pixels
        outerrad : Alignment radius, pixels
        outprojparamsfile : Output projection parameter file
        out2daligndoc : Output 2D alignment file
        interpolation_method : Interpolation method
        outliers : Threshold for allowed angular deviation from average, degrees
        goodclassparttemplate : Particles within outlier threshold
        alignopt : Alignment method
        ringstep : Alignment radius increment, pixels
        log : Logger object
        ####verbose : (boolean) Whether to write additional information to screen
        verbosity : Verbosity level (0..2)
        ####debug : Prints additional information to screen
    """
    
    # Resample reference
    recondata= EMAN2.EMData(reconfile)
    idim= recondata['nx']
    reconprep= sp_projection.prep_vol(recondata, npad=2, interpolation_method=interpolation_method)
    
    # Initialize output angles
    outprojparamslist= []
    out2dalignlist= []

    # Read class lists
    classdoclist= glob.glob(classdocs)
    refineparamslist= sp_utilities.read_text_row(refineparams)
    
    # Allocate number of digits for class number
    num_imgs  = EMAN2.EMUtil.get_image_count(classavgstack)
    num_digits= math.ceil(math.log10(num_imgs) )
    
    if verbosity >=1 : print_log_msg('Averaging Euler angles in each class', log)
    
    disableTF= verbosity<1 or verbosity>1

    # Loop through class lists
    for classidx in tqdm.tqdm(range( len(classdoclist) ), unit='class', disable=disableTF):
        classdoc= classdoclist[classidx]
        
        # Strip out filenumber
        classexample= os.path.splitext(classdoc)
        classnum= int(classexample[0][-num_digits:])  # number of digits in class number
        if verbosity >=2 : print_log_msg('Index %s, class number %s' % (classidx, classnum), log)
        
        # Initial average
        [avg_phi_init, avg_theta_init]= average_angles(
            refineparamslist, 
            classdoc, 
            selectdoc=selectdoc, 
            log=log, 
            verbose=verbosity>=2
            )
        
        # Look for outliers
        if outliers:
            goodpartdoc= goodclassparttemplate.format(classnum)
            [avg_phi_final, avg_theta_final]= average_angles(
                refineparamslist, 
                classdoc, 
                selectdoc=selectdoc,
                init_angles=[avg_phi_init, avg_theta_init], 
                threshold=outliers, 
                goodpartdoc=goodpartdoc, 
                log=log, 
                verbose=verbosity>=2
                )
        else:
            [avg_phi_final, avg_theta_final]= [avg_phi_init, avg_theta_init]
        
        # Compute re-projection
        refprjreal= sp_projection.prgl(
            reconprep, [avg_phi_final,avg_theta_final,0,0,0], interpolation_method=1, return_real=True
            )
        
        # Align to class average
        classavg= sp_utilities.get_im(classavgstack, classnum)
        
        # Alignment using self-correlation function
        if alignopt == 'scf':
            ang_align2d, sxs, sys, mirrorflag, peak= sp_alignment.align2d_scf(
                classavg, refprjreal, maxshift, maxshift, ou=outerrad
                )
        
        # Weird results
        elif alignopt == 'align2d':
            # Calculate outer radius (default of -1 not OK)
            if outerrad <= 0 : outerrad= old_div(idim, 2) - 2
            
            # Set search range
            currshift= 0
            txrng= tyrng= sp_alignment.search_range(idim, outerrad, currshift, maxshift)
            
            # Perform alignment
            ang_align2d, sxs, sys, mirrorflag, peak= sp_alignment.align2d(
                classavg, refprjreal, txrng, tyrng, last_ring=outerrad
                )  
        
        # Direct3 (angles seemed to be quantized)
        elif alignopt == 'direct3':
            [[ang_align2d, sxs, sys, mirrorflag, peak]]= sp_alignment.align2d_direct3(
                [classavg], refprjreal, maxshift, maxshift, ou=outerrad
                )
        
        # APSH-like alignment (default)
        else:
            [[ang_align2d, sxs, sys, mirrorflag, scale]]= apsh(refprjreal, classavg, 
                outerradius=outerrad, maxshift=maxshift, ringstep=ringstep)
            
        out2dalignlist.append([ang_align2d, sxs, sys, mirrorflag, 1])
        msg= "ang_align2d=%s sx=%s sy=%s mirror=%s" % (ang_align2d, sxs, sys, mirrorflag)
        if outliers != None and verbosity < 2: 
            msg += "\n"
        if verbosity >=2 : print_log_msg(msg, log)
        
        # Check for mirroring
        if mirrorflag == 1:
            tempeulers= list(
                sp_utilities.compose_transform3(
                    avg_phi_final,avg_theta_final,0, 0,0,0, 1, 
                    0,            180,            0, 0,0,0, 1
                    )
                )
            combinedparams= list(
                sp_utilities.compose_transform3(
                    tempeulers[0],tempeulers[1],tempeulers[2], tempeulers[3],tempeulers[4],0, 1, 
                    0,            0,            -ang_align2d,  0,            0,            0, 1
                    )
                )
        else:
            combinedparams= list(
                sp_utilities.compose_transform3(
                    avg_phi_final,avg_theta_final,0,            0,0,0, 1,
                    0,            0,              -ang_align2d, 0,0, 0,1
                    )
                )
        # compose_transform3: returns phi,theta,psi, tx,ty,tz, scale
        
        if verbosity >=2 : 
            print_log_msg("combinedparams %s\n" % (combinedparams), log)
        outprojparamslist.append(combinedparams)
    # End class-loop
    
    sp_utilities.write_text_row(outprojparamslist, outprojparamsfile)
    sp_utilities.write_text_row(out2dalignlist, out2daligndoc)
    if verbosity>=1:
        mesg= 'Wrote alignment parameters to %s and %s\n' % (outprojparamsfile, out2daligndoc)
        print_log_msg(mesg, log)
    
    del recondata  # Clean up
        
    
def average_angles(alignlist, partdoc, 
                   selectdoc=None, init_angles=None, threshold=None, 
                   goodpartdoc=None, log=None, verbose=False):
    """
    Computes a vector average of a set of particles' Euler angles phi and theta.
    
    Arguments:
        alignlist : Alignment parameter doc file, i.e., from MERIDIEN refinement
        partdoc : List of particle indices whose angles should be averaged
        selectdoc : Input substack selection file if particles removed before refinement (e.g., Substack/isac_substack_particle_id_list.txt)
        init_angles : List (2 elements) with initial phi and theta angles, for excluding outliers
        threshold : Angular threshold (degrees) beyond which particles exceeding this angular difference from init_angles will be excluded
        goodpartdoc : Output list of retained particles if a threshold was specified
        log : Logger object
        verbose : (boolean) Whether to write additional information to screen
    Returns:
        list of 2 elements:
            avg_phi
            avg_theta
    """
    
    # Read alignment parameters
    if isinstance(alignlist, str): alignlist= sp_utilities.read_text_row(out2daligndoc)
    # (If loading the same parameters repeatedly, better to read the file once externally and pass only the list.)
    
    # Read class list
    partlist= sp_utilities.read_text_row(partdoc)
            
    if selectdoc: 
        selectlist= sp_utilities.read_text_row(selectdoc)
    else:
        selectlist= None
        
    sum_phi= np.array([0.0,0.0])
    sum_theta= np.array([0.0,0.0])
    totparts= 0
    num_outliers= 0
    goodpartlist= []
    goodpartcounter= 0
    
    # Loop through particles
    for totpartnum in partlist:
        if selectlist:
            goodpartnum= selectlist.index(totpartnum)
        else:
            goodpartnum= totpartnum[0]

        try:
            phi_deg= alignlist[goodpartnum][0]
            theta_deg= alignlist[goodpartnum][1]
            phi_rad= np.deg2rad(phi_deg)
            theta_rad=  np.deg2rad(theta_deg)
        except IndexError:
            msg= "\nERROR!! %s tries to access particle #%s" % (partdoc, goodpartnum)
            numalignparts= len(alignlist)
            msg += "\nAlignment doc file has only %s entries" % (numalignparts)
            msg += "\nMaybe try substack selection file with flag '--select <substack_select>'?"
            sp_global_def.ERROR(msg, __file__, 1)
            exit()
        
        if init_angles:
            angdif= sp_pixel_error.angle_diff(init_angles, [phi_deg, theta_deg])
            if angdif > 180: angdif= 360.0 - angdif
        
        totparts += 1
        
        # Exclude particles exceeding optional threshold
        if threshold==None or angdif < threshold:
            sum_phi += (np.cos(phi_rad), np.sin(phi_rad))
            sum_theta += (np.cos(theta_rad), np.sin(theta_rad))
            goodpartlist.append(totpartnum)
            goodpartcounter += 1
        else:
            num_outliers += 1
        
    # Compute final average
    avg_phi= math.degrees(math.atan2(sum_phi[1],sum_phi[0]))
    avg_theta= math.degrees(math.atan2(sum_theta[1],sum_theta[0]))
    
    # Clean up, might reuse
    del alignlist
    del partlist
    del selectlist
    
    if verbose : 
        ####msg= "Particle list %s: average angles (%s, %s)" % (partdoc, avg_phi, avg_theta)
        msg= "Average angles (%s, %s)" % (avg_phi, avg_theta)
        print_log_msg(msg, log)

    if threshold:
        if verbose : 
            msg= "Found %s out of %s outliers exceeding an angle difference of %s degrees from initial estimate" % (num_outliers, totparts, threshold)
            print_log_msg(msg, log)
        
        if goodpartdoc:
            if goodpartcounter > 0:
                sp_utilities.write_text_row(goodpartlist, goodpartdoc)
                if verbose : 
                    msg= "Wrote %s particles to %s" % (goodpartcounter, goodpartdoc)
                    print_log_msg(msg, log)
            else:
                if verbose : 
                    msg= "WARNING!! Kept 0 particles from class %s" % partdoc
                    print_log_msg(msg, log)
                [avg_phi, avg_theta]= init_angles
                
    return [avg_phi, avg_theta]
    
def compare_projs(reconfile, classavgstack, outdir, inputprojparams=None, interpolation_method=1, 
                classdoctemplate=None, log=None, verbosity=0):
    """
    Make comparison stack between class averages (even-numbered (starts from 0)) and re-projections (odd-numbered).
    
    Arguments:
        reconfile : Input volume from which to generate re-projections
        classavgstack ; Input image stack
        outdir ; Output directory
        inputprojparams : Input Euler angles doc
        interpolation_method : Interpolation method: nearest neighbor (nn, 0), trilinear (1, default), gridding (-1)
        classdoctemplate : Template for class-particle template (if Meridien/outlier mode)
        log : Logger object
        verbosity : Verbosity level (0..2)
    Returns:
        compstack : Stack of comparisons between input image stack 
    """
    
    nimg1  = EMAN2.EMUtil.get_image_count(classavgstack)
    if inputprojparams : angleslist= sp_utilities.read_text_row(inputprojparams)  # if not provided will try to read from header
    
    recondata= EMAN2.EMData(reconfile)
    nx= recondata.get_xsize()
    
    # Check whether dimension is odd (prgl does something weird if so. --Tapu)
    padYN= False
    if nx % 2 != 0:
        padYN= True
        nx= nx + 1
        if verbosity>=1 : print_log_msg("Padding volume and images to %s" % nx, log)
        recondata= sp_utilities.pad(recondata, nx, nx, nx, background='circumference')
    
    # Check whether images have class-membership information
    classimg= sp_utilities.get_im(classavgstack, 0)
    if not classimg.has_attr('members'):
        msg= 'No class membership information found'
        print_log_msg(msg, log)
    
    # Resample reference
    reconprep= sp_projection.prep_vol(recondata, npad=2, interpolation_method=interpolation_method)

    ccclist= []
    
    #  Here you need actual radius to compute proper ccc's, but if you do, you have to deal with translations, PAP
    mask= sp_utilities.model_circle(nx//2-2,nx,nx)
    mask.write_image(os.path.join(outdir, MASIMG) )
    compstack= os.path.join(outdir, COMPARISON_STACK)
    
    if verbosity >=1 : print_log_msg('Comparing images', log)
    
    disableTF= verbosity<1
    
    # Loop through images
    for classnum in tqdm.tqdm(range(nimg1), unit=' class', disable=disableTF):
        # Get class average
        classimg= sp_utilities.get_im(classavgstack, classnum)
        
        # Pad if necessary
        if padYN:
            classimg= pad(classimg, nx, nx, background='circumference')
        
        # Get angles either from file or header
        if inputprojparams: 
            params3d= angleslist[classnum]
        else:
            assert classimg.has_attr('xform.projection'), "ERROR! Projection angles provided neither in header nor file"
            params3d= sp_utilities.get_params_proj(classimg)
            
        # Compute re-projection
        prjimg= sp_projection.prgl(reconprep, params3d, interpolation_method=1, return_real=False)
        
        # Calculate 1D power spectra
        rops_dst= sp_fundamentals.rops_table(classimg*mask)  
        rops_src= sp_fundamentals.rops_table(prjimg)
        
        #  Set power spectrum of reprojection to the data.
        #  Since data has an envelope, it would make more sense to set data to reconstruction,
        #  but to do it one would have to know the actual resolution of the data. 
        #  you can check sxprocess.py --adjpw to see how this is done properly  PAP
        table= [0.0]*len(rops_dst)  # initialize table
        for j in range( len(rops_dst) ):
            table[j]= math.sqrt( old_div( rops_dst[j], rops_src[j]) )
        
        # match FFT amplitudes of re-projection and class average
        prjimg= sp_fundamentals.fft(sp_filter.filt_table(prjimg, table))  

        cccoeff= sp_statistics.ccc(prjimg, classimg, mask)
        #print classnum, cccoeff
        classimg.set_attr_dict({'cross-corr':cccoeff})
        prjimg.set_attr_dict({'cross-corr':cccoeff})
        
        montagestack= []
        montagestack.append(prjimg)
        montagestack.append(classimg)
        
        # Copy membership information
        if classdoctemplate:
            goodpartdoc= classdoctemplate.format(classnum)
            
            # After outlier-removal, may be empty
            if os.path.exists(goodpartdoc): 
                members= [item for sublist in sp_utilities.read_text_row(goodpartdoc) for item in sublist]
                # (read_text_row returns a list of N lists of length 1.)
            else:
                members= None
        else:
            if classimg.has_attr('members'):
                members= classimg.get_attr('members')
            else:
                members= None
        
        comparison_pair= montage2(montagestack, ncol=2, marginwidth=1)
        if members : comparison_pair.set_attr('members', members)
        comparison_pair.write_image(compstack,classnum)
        
        ccclist.append(cccoeff)
    del angleslist
    meanccc= old_div(sum(ccclist), nimg1)
    if verbosity>=1 : print_log_msg("Average CCC is %s\n" % meanccc, log)
    
    nimg2= EMAN2.EMUtil.get_image_count(compstack)
    
    for classnum in range(nimg2):  # xrange will be deprecated in Python3
        prjimg= sp_utilities.get_im(compstack,classnum)
        meanccc1= prjimg.get_attr_default('mean-cross-corr', -1.0)
        prjimg.set_attr_dict({'mean-cross-corr':meanccc})
        sp_utilities.write_header(compstack,prjimg,classnum)
    
    return compstack
    
def montage2(inputstack, ncol, marginwidth=0, bkgd=0, outfile=None):
    """
    Generates montage of images into one image.
    Adapted from sxmontage.py
    
    Arguments:
        inputstack : Stack of input images to merge into montage
        ncol : Number of images per row
        marginwidth : Margin width, pixels
        bkgd : Background value of montage
        outfile : Optional output file with montage output
    Returns:
        montage : EMData object of image montage
    """
    
    if isinstance(inputstack, str): inputstack= EMAN2.EMData.read_images(inputstack)
    
    # Get single-image dimensions
    nx= inputstack[0].get_xsize()
    ny= inputstack[0].get_ysize()
    
    # Get number of images and calculate montage dimensions
    numimgs= len(inputstack)
    numrows= old_div( (numimgs-1), ncol ) + 1
    
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
        insert_image(inputstack[imgnum], montage, xoffset, yoffset)
    
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

def main():
    options= parse_command_line()
    
    ##print args, options  # (Everything is in options.)
    #print('options', options)
    #print('LOGFILE', sp_global_def.LOGFILE)
    #exit()
    
    # If output directory not specified, write to same directory as class averages
    if not options.outdir:
        outdir= os.path.dirname(os.path.realpath(options.classavgs))
    else:
        outdir= options.outdir

    if options.mode == 'viper':
        selectdoc= options.projselect
    elif options.mode == 'projmatch':
        selectdoc= None
    elif options.mode == 'meridien':
        selectdoc= options.partselect
    else:
        sp_global_def.ERROR("\nERROR!! Valid mode not specified. Valid modes are: viper, projmatch, and meridien.", __file__, 1)
        sp_global_def.sxprint('Type %s --help to see available options\n' % os.path.basename(__file__))
        exit()

    compare_projections(
        options.classavgs, 
        options.vol3d, 
        outdir, 
        mode=options.mode, 
        projmethod=options.projmethod, 
        projparamfile=options.projparams, 
        refineparams=options.refineparams, 
        selectdoc=selectdoc, 
        displayYN=options.display, 
        verbosity=options.verbosity, 
        options=options 
        )
    sp_global_def.print_timestamp( "Finish" )

if __name__ == "__main__":
    main()
