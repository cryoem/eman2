#!/usr/bin/env python
from __future__ import print_function
#
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2019 (toshio.moriya@kek.jp)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
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
#

from __future__ import division
from past.utils import old_div
import os
import EMAN2
import datetime
import operator
import EMAN2db
import argparse
import math
import tqdm
import sys
import subprocess
import numpy as np
from matplotlib import pyplot as plt
import mpi
import random
import time

try:
    from sphire.libpy import sp_global_def
    from sphire.libpy import sp_utilities
    from sphire.libpy import sp_logger
    from sphire.libpy import sp_applications
    from sphire.libpy import sp_fundamentals
    from sphire.libpy import sp_filter
    from sphire.libpy import sp_statistics
    from sphire.libpy import sp_morphology
    ##from sp_utilities import montage_scale  #### Moved to sp_utilities in SPHIRE 1.3
except ImportError:
    import sp_global_def
    import sp_utilities
    import sp_logger
    import sp_applications
    import sp_fundamentals
    import sp_filter
    import sp_statistics
    import sp_morphology
    ##from sp_utilities import montage_scale  #### Moved to sp_utilities in SPHIRE 1.3
    
# Set default values (global variables written in ALL CAPS)
TIMESTAMP_LENGTH=  23               # Padding length for timestamps
MAX_VERBOSITY=     6                # Maximum verbosity
GLOBAL_NUM_TAG=   'global_num'      # Header tag for fixed, consecutive particle number
ROTATION_TAG=     'rot.align2d'     # Header tag for rotation angle
CENTERING_TAG=    'xform.center2d'  # Header tag for centering (from sp_center_2d3d.py)

# Inputs
ISAC_ALIGN_DIR=      "2dalignment"                 # ISAC 2d alignment directory
ISAC_INIT_PARAMS=    "initial2Dparams.txt"         # ISAC initial alignment parameters
ISAC_FINAL_PARAMS=   "all_parameters.txt"          # ISAC final alignment parameters
ISAC_SHRINK_RATIO=   "README_shrink_ratio.txt"     # ISAC shrink ratio file
ISAC_PROCESSED_IMGS= 'processed_images.txt'        # ISAC processed images
ORDERED_CLASS_AVGS=  'ordered_class_averages.hdf'  # ISAC ordered class averages
BEAUT_FINAL_PARAMS=  'ali2d_local_params.txt'      # Beautifier final parameters
BEAUT_INIT_PARAMS=   'init_isac_params.txt'        # Beautifier initial parameters (from ISAC)
ORDERED_PARAMS=      'chains_params.txt'           # Ordered alignment parameters (requires new sp_chains.py)

# Outputs
IMGFORMAT=            ".mrcs"                 # Stack-file format
DOCFILEDIR=           "ClassDocs"             # Class doc file directory
CLASSMAPFILE=         "classmap.txt"          # Class-to-particle lookup table
BIGSTACKCOPY=         'all_copy'              # Copy of input stack
CLASSDOCPREFIX=       'sel_class_'            # Class-membership files
OLD2NEW_PARTICLES=    'sel_particles.txt'     # List of includes particles, using original numbering
STACKFILEDIR=         "Stacks"                # Stack file directory
LOCALMASKS=           "stkmask.hdf"           # Mask stack
CLASSORIGBDBPREFIX=   'original_'             # Original particles, separated into classes
CLASS_STACK_PREFIX=   'stkclass_'             # Filtered stack prefix
CLASS_BDB_STEM=       'class_'                # Stem of output BDB stack
FILTCLASSBDBPREFIX=   'filtered_'             # Filtered BDB prefix
ALIGNBDBPREFIX=       'aligned_'              # Aligned BDB prefix
PARAMS_PREFIX=        'params_'               # Parameter file prefix
MPI_STACK_PREFIX=     'ave_mpi_'              # Stacks with average, variance, and eigenimages (one per MPI process)
AVGVARSTACK=          'stkavgvar.hdf'         # Average + variance stack
EIGENIMGSTACK=        'stkeigen.hdf'          # Eigenimage stack
AVGVARBDB=            'avgvariance'           # Average + variance BDB
EIGENIMGBDB=          'eigenimages'           # Eigenimage BDB
CONSECPARTICLES=      'listparticles.txt'     # List of all particles, consecutively numbered
BIG_ALIGN_BDB_PREFIX= 'allclasses_'           # Combined, aligned BDB stack
CLASS_CENTER_DOC=     'center_params.txt'     # Class-centering parameters
PWS_1D_PREFIX=        'docpws'                # (CTF phase-flip mode) Prefix for 1D power spectra
BANDPASS_AVGS=        'stkbandpass.hdf'       # Band-pass filtered averages
BANDPASS_FILTER_DOC=  'docbandpass.txt'       # Band-pass filter profile text file
BANDPASS_FILTER_PLOT= 'plotbandpass.png'      # Band-pass filter profile text file

USAGE= """ 
PURPOSE:
Separate particles according to class assignment.

General usage:
%s <input_avgs> <input_images> <output_directory> 
Required command-line parameters:
1. Input class averages
2. Input image stack
3. Output directory
Outputs:
%s/%s : Class-to-particle lookup table, one file for all classes 
%s/%s???.txt : List of particles for each class, one file per class
EMAN2DB/%s???.bdb : Virtual stacks of particles for each class
%s/%s???.mrcs : (Optional) stacks of aligned particles for each class
%s/%s???.mrcs : (Optional) stacks of filtered/shrunken particles for each class

Advanced usage:
%s <input_avgs> <input_images> <output_directory> --verbosity=<verbosity>
Advanced parameters:
--verbosity : Adjust verbosity (0..4)

To apply a low-pass (tangent) filter:
%s <input_avgs> <input_images> <output_directory> --filtrad=<filter_radius> --apix=<pixel_size>
Parameters:
--filtrad : Low-pass filter radius, Angstroms (or if, apix not provided, pixels^-1)
--apix : Pixel size of input images (NOT class averages), Angstroms (optional, if not provided, filter radius assumed to be pixels^-1)

To downsample the output images:
%s <input_avgs> <input_images> <output_directory> --filtrad=<filter_radius> --shrink=<shrink_factor> --apix=<pixel_size>
Parameter:
--shrink : Downsampling factor (e.g., 4 -> 1/4 original size)

To apply alignments from ISAC to output image stacks, also creates 2D averages + variances:
%s <input_avgs> <input_images> <output_directory> --align_isac_dir <ISAC_directory> --ctf=wiener --chains_radius=<alignment_radius>
Parameter:
--align_isac_dir : If applying alignments, directory for ISAC output
Advanced Parameters:
--ctf : CTF-correction options: flip, wiener
--chains_radius : Radius (pixels) for internal alignment of averages, overrides existing ordered class-average alignments

To compute eigenimages:
%s <input_avgs> <input_images> <output_directory> --align_isac_dir <ISAC_directory> --nvec=<number_of_factors> --pca_radius=<radius> --mask_drawn <drawn_mask_file> --mask_binary <binary_mask>
Parameter:
--nvec : Number of (factors) eigenimages to compute
Advanced Parameters:
--pca_radius : Radius for masking for PCA, for full-size images, pixels (default: half image dimension)
--mask_drawn : Mask file for PCA (draw using e2display.py in single-image mode, with pen intensity greater than image maximum)
--mask_binary : Binary mask file for PCA

To apply class-centering parameters (i.e., from sp_center_2d3d.py):
%s <input_avgs> <input_images> <output_directory> --align_isac_dir <ISAC_directory>  --write_centered --applyparams=<align_option> --debug
Parameter:
--write_centered : Write centered images (assuming images centered using sp_center_2d3d.py)
Advanced Parameters:
--applyparams : Option to apply alignments to outputs: None (default), combined (rotation + shift), intshifts (integer shifts)
--debug : Applies alignment option and computes average, for testing

""" % (__file__, DOCFILEDIR, CLASSMAPFILE, DOCFILEDIR, CLASSDOCPREFIX, CLASSORIGBDBPREFIX, 
    STACKFILEDIR, CLASS_STACK_PREFIX, STACKFILEDIR, CLASS_STACK_PREFIX, 
    __file__,__file__,__file__,__file__,__file__,__file__,)

MODIFIED="Modified 2021-01-01"

"""
Modifications log:
    TODO: bcast average
    2020-11-06 (trs) -- can use binary mask for PCA, e.g., from sp_center_2d3d
    2020-11-06 (trs & bv) -- without '--particle' flag, simply writes class doc files
    2020-11-02 (trs & dq) -- added custom mask option, drawn over class average (using e2display)
    2020-10-20 (trs) -- uses MPI
    2020-10-06 (trs) -- if ISAC/beautifier files are missing, will list them
    2020-07-23 (trs) -- input particle stack now is preceded by '--particles' flag
    2020-07-23 (trs) -- added band-pass filter option
    2020-07-02 (trs) -- added class size (n_objects) to class averages
    2020-05-23 (trs & ac) -- bug fix when running PCA with no shrink factor
    2020-05-20 (trs) -- can use subsets of class averages
    2020-05-16 (trs) -- reads ordered class-average parameters from header from new sp_chains.py
    2020-05-15 (trs) -- can read alignment parameters from beautifier directory
    2020-05-04 (trs & ac) -- optional CTF flag applied to averages
    2020-04-25 (trs & dq) -- converts 2D parameters to 3D (to use for restacking)
    2020-04-23 (trs) --  sorts BDB stack with original numbering
    2020-01-21 (trs & dq) -- added CTF-correction to particles
    2020-01-17 (trs & dq) -- added radius option for PCA
"""

class BdbNames:
    """
    Adapted from sp_signalsubtract.py
    
    Class containing the following attributes:
        bdb_stem : Stem, no directory, no leading "bdb:", no trailing ".bdb"
        bdb_dir : Directory, not counting "EMAN2DB"
        eman2db_dir : Directory, including "EMAN2DB"
        bdb_name : Of the form "bdb:DirectoryName#stem"
        bdb_path : Of the form "DirectoryName/EMAN2DB/stem.bdb"
    """

    def __init__(self, bdb_stem, dir_above_eman2db=None, quiet=False):
        """
        Parameters:
            bdb_stem : Name of BDB not including preceding 'bdb:' or trailing '.bdb'
            dir_above_eman2db : Directory below which EMAN2DB/<bdb_stem>.bdb will be written
            quiet : (boolean) Flag to suppress warnings
        """
        
        # TODO: sanity check for input stem, make sure it doesn't start with 'bdb:', etc.
        self.bdb_stem= bdb_stem
        if dir_above_eman2db != None:
            self.bdb_dir= dir_above_eman2db
        else:
            self.bdb_dir= '.'
        self.eman2db_dir= os.path.join(self.bdb_dir, 'EMAN2DB')
        self.bdb_name= self.stem2name()
        self.bdb_path= self.stem2path()

    def stem2name(self):
        if self.bdb_dir != None:
            name= 'bdb:' + self.bdb_dir + "#" + self.bdb_stem
        else:
            name= 'bdb:' + self.bdb_stem
        
        return name

    def stem2path(self):
        if self.bdb_dir != None:
            path= os.path.join(self.bdb_dir, 'EMAN2DB', self.bdb_stem + '.bdb')
        else:
            path= os.path.join('EMAN2DB', self.bdb_stem + '.bdb')
            
        return path
    
    def list_info(self):
        info= []
        info.append(self.bdb_stem)
        info.append(self.bdb_dir)
        info.append(self.eman2db_dir)
        info.append(self.bdb_name)
        info.append(self.bdb_path)
        
        return info
        
def setup_mpi(usempi, verbose=False):
    """
    Set up MPI.
    
    Argument:
        usempi : (boolean) Whether to use Message Passing Interface
        
    Returns:
        mpidict : dictionary with keys:
            'use' : (boolean) whether to use MPI
            'myid' : MPI process number
            'main_process' : main MPI process number
            'is_main' : (boolean) whether current process is main
            'size' : number of processes
    """
    
    main_process= 0  # default

    ####usempi= "OMPI_COMM_WORLD_SIZE" in os.environ

    if usempi:
        number_of_proc= mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
        myid= mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
        
        #  choose a random node as main
        if myid == 0: main_process= random.randint(0, number_of_proc-1)
        main_process= mpi.mpi_bcast(main_process, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD)
    else:
        myid= 0
        number_of_proc= 1
        # (IF-THEN statements will be simpler if I pretend non-MPI runs are the main node.)
    
    is_main= myid==int(main_process)  # boolean: whether current process is the main one
    mpidict= {'use':usempi, 'myid':myid, 'main_process':main_process, 'is_main':is_main, 'size': number_of_proc}
    
    if verbose: 
        for key in mpidict.keys() : print("%s %s" % (key, mpidict[key]) )
        print()

    return mpidict
    
def quick_barrier(mesg=None):
    """
    Synchronizes parallel processes before continuing
    Safe for non-MPI calls also
    """
    
    if "OMPI_COMM_WORLD_SIZE" in os.environ:
        if mesg!=None : print("Synchronizing: %s" % mesg)
    
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        
        if mesg!=None : print("Synchronized: %s" % mesg)
    
def safe_exit(verbosity=0):
        """
        Properly closes MPI before exiting, preventing ugly warnings
        Safe for non-MPI calls also

        Argument:
                verbosity : controls how much information to write to screen (0..3)
        """
        
        if "OMPI_COMM_WORLD_SIZE" in os.environ: 
                my_rank= mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD) 
                mpi.mpi_finalize()
        else:
                my_rank= 0
        
        sp_global_def.print_timestamp( "Finish" )
        
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
            log= sp_logger.Logger(base_logger=sp_logger.BaseLogger_Files(), base_logger2=sp_logger.BaseLogger_Print(), file_name=logname)
            verbose= False  # logger output will be echoed to screen
        else:
            log= sp_logger.Logger(base_logger=sp_logger.BaseLogger_Files(), file_name=logname)
    except TypeError:
        print("WARNING: Using old sp_logger.py library")
        log= sp_logger.Logger(base_logger=sp_logger.BaseLogger_Files())#, file_name=logname)
        logname= 'log.txt'
        
    if main:
        print("Writing log file to %s" % logname)
        
        progbase= os.path.basename(__file__).split('.')[0].upper()
        length= len(progbase) + 4
        
        log.add("\n" +
                " "*TIMESTAMP_LENGTH + "*"*length + "\n" +
                " "*TIMESTAMP_LENGTH + "* " + progbase + " *\n" +
                " "*TIMESTAMP_LENGTH + "*"*length)
    
    return log, verbose

def print_log_msg(mesg, log=None, is_main=True):
    """
    Prints messages to log file and, optionally, to the screen.
    
    Arguments:
        mesg : message to write
        log : instance of Logger class
        is_main : (boolean) if using multiple cores, some tasks only need to be performed once, not by all cores
    """
    
    if is_main:
        if log: 
            log.add(mesg)
        else:
            print(mesg)

def printvars(variables, quitTF=False, typeTF=False):
    """Print the local variables in the caller's frame.
    
    Adapted from https://stackoverflow.com/questions/6618795/get-locals-from-calling-namespace-in-python
    """
    
    import inspect
    import six

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
        mesg= "%s: %s" % (variable, dictionary[variable])
        if typeTF: mesg+= " %s" % type(dictionary[variable])
        print(mesg)
        
    del frame
    
    if quitTF: 
        print('\nExiting printvars...')  # reminder in case I forget to take out the quit flag
        safe_exit()

def setup_subprocess():
    """
    Returns
        Python version number
        allowed instance of DEVNULL (specific to Python version)
        
    Adapted from https://stackoverflow.com/questions/11269575/how-to-hide-output-of-subprocess-in-python-2-7
    """
    
    python_version= sys.version.split()[0]
    
    # Keep up to one decimal
    version_float= float('.'.join(python_version.split('.', 2)[:2]))
    
    try:
        from subprocess import DEVNULL
    except ImportError:
        DEVNULL = open(os.devnull, 'wb')
    
    return version_float, DEVNULL

def system_call_23(cmd, args, lenient=False, stdout=None, stderr=None, usempi=False, log=None, verbose=False):
    """
    Runs subprocess safely.
    
    Arguments:
        cmd : Executable
        args : Command-line arguments
        lenient : (boolean) Will simply print warning in case of error
        stdout : Where to direct standard out
        stderr : Where to direct standard error
        usempi : (boolean) Whether using MPI, in order to remove "OMPI_COMM_WORLD_RANK" from os.environ if necessary
        log : Logger instance
        verbose : (boolean) Whether to write to screen
    """
    
    # Check Python version
    python_version= sys.version.split()[0]
    
    # Keep up to one decimal
    version_float= float('.'.join(python_version.split('.', 2)[:2]))
    
    # Build command line
    cmdline= "%s %s" % (cmd, args)
    
    if verbose : print_log_msg(cmdline, log)
    
    # Some SPHIRE programs call sp_global_def which checks for pre-existing "OMPI_COMM_WORLD_RANK" values
    if usempi:
        mpi_rank= os.environ["OMPI_COMM_WORLD_RANK"]
        del os.environ["OMPI_COMM_WORLD_RANK"]
    
    try:
        if version_float < 3.5:
            #subprocess.check_call([cmd, args], stdout=stdout, stderr=stderr, shell=True)
            subprocess.check_call(cmdline, stdout=stdout, stderr=stderr, shell=True)
            # (shell=True is risky, but I couldn't get it to work without)
        else:
            subprocess.run([cmd] + args.split(), stdout=stdout, stderr=stderr)
    except subprocess.CalledProcessError:
        if not lenient:
            print("\nERROR!! Cannot execute '%s'." % cmdline)
            if "OMPI_COMM_WORLD_RANK" in os.environ:
                print("Maybe try to remove 'OMPI_COMM_WORLD_RANK' from os.environ by using 'usempi=True' in 'system_call_23'.")
            print()
            exit()
        else:
            print("\nWARNING! Cannot execute '%s'\n" % cmdline)

    if usempi:
        os.environ["OMPI_COMM_WORLD_RANK"]= mpi_rank
    
def eval_isac(classavgstack, options, outdir='.', verbosity=0, usempi=False):
    """
    Main function overseeing various ISAC-evaluation modes.
    
    Arguments:
        classavgstack : input class-average stack
        options : (namespace) command-line options, run 'sp_eval_isac.py -h' for an exhaustive list
        outdir : output directory
        verbosity : verbosity level
        usempi : (boolean) whether to use parallelization
    """
    
    #Set up MPI
    mpioptions= setup_mpi(usempi)
    is_main= mpioptions['is_main']
    
    # Set output directory and log file name
    if is_main:
        prepare_outdir(outdir, verbose=verbosity>=1)
        prepare_outdir(os.path.join(outdir, DOCFILEDIR), verbose=verbosity>=1)
        if options.align_isac_dir or options.filtrad or options.shrink!=1:
            prepare_outdir(os.path.join(outdir, STACKFILEDIR) )
    
    quick_barrier()
    log, _= prepare_log(outdir, main=is_main)
    sp_global_def.write_command(outdir)  # expects mpi_comm_rank to be 0, which isn't necessarily main process, which is randomly assigned
    
    if verbosity>=1 : 
        print_log_msg('Using verbosity level %s (maximum %s)' % (verbosity, MAX_VERBOSITY), log, is_main)
        
    if options.bandpass_radius:
        if is_main:
            bandpass_filter(
                classavgstack, 
                options, 
                outdir, 
                log=log, 
                verbosity=verbosity
                )
    else:
        separate_classes(
            classavgstack, 
            options, 
            outdir, 
            verbosity=verbosity,
            log=log, 
            mpioptions=mpioptions
            )
    
def separate_classes(classavgstack, options, outdir='.', verbosity=0, log=None, mpioptions=None):
    """
    Oversees various class-sepataion modes.
    
    Arguments:
        classavgstack : Input class-average stack
        options : (namespace) Command-line options, run 'sp_eval_isac.py -h' for an exhaustive list
        outdir : Output directory
        verbosity : Verbosity level
        log : Instance of Logger class
        mpioptions : (dict) MPI options
    """
    
    # Parallelization stuff
    is_main= mpioptions['is_main']
    main_mpi_proc= mpioptions['main_process']
    myid= mpioptions['myid']
    usempi= mpioptions['use']
    
    # Make sure enough digits
    num_classes= EMAN2.EMUtil.get_image_count(classavgstack)
    num_digits= math.ceil( math.log10(num_classes) )
    digit_pattern= '{{0:0{0}d}}'.format(num_digits)

    if options.particles==None:
        if verbosity>=1:
            mesg= "Simply splitting '%s' into text files with lists of members" % classavgstack
            print_log_msg(mesg, log, is_main)
        
        # Print warning if extraneous flags were provided
        if options.align_isac_dir or options.applyparams or options.chains_radius or options.filtrad or options.nvec!=0:
            print_log_msg("WARNING! Further options require the particle stack, i.e., '--particles PARTICLE_STACK'", log, is_main)
    else:
        partstack= options.particles
        
        # Sanity check
        try:
            num_parts= EMAN2.EMUtil.get_image_count(partstack)
        except:
            print_log_msg("ERROR!! Couldn't open '%s'!" % partstack, log, is_main)
            safe_exit()
        
        if verbosity>=1 : print_log_msg("Using particle stack '%s'" % partstack, log, is_main)
    
    if options.filtrad: 
        if options.apix: 
            filter_radius= old_div(options.apix, options.filtrad)
            mesg= f"Will low-pass filter to {filter_radius: .1f} px^-1 (={options.apix}/{options.filtrad})"
        else:
            filter_radius= options.filtrad
            mesg= f"Will low-pass filter to {filter_radius: .1f} px^-1"
        if verbosity>=1 : print_log_msg(mesg, log, is_main)
    else:
        filter_radius= None
    
    if options.shrink!=1 and verbosity>=1: 
        print_log_msg("Will downsample stacks by a factor of %s" % options.shrink, log, is_main)
    
    if options.align_isac_dir!=None and verbosity>=1: 
        print_log_msg("Will apply alignments from '%s'" % options.align_isac_dir, log, is_main)
    
    if options.nvec>0 and options.align_isac_dir == None:
        print()
        print_log_msg("ERROR!! To compute eigenimages, need to specify --align_isac_dir", log, is_main)
        safe_exit()
    
    print()
    if options.ctf != None:
        if options.align_isac_dir:
            if options.ctf == 'flip':
                if verbosity>=1:
                    print_log_msg('Phase-flipping images to compute averages', log, is_main)
                    mesg= 'Program will also compute 1D power spectra profiles with and w/o CTF-correction\n'
                    print_log_msg(mesg, log, is_main)
            elif options.ctf == 'wiener':
                if verbosity>=1 : print_log_msg('Wiener-filtering images to compute averages', log, is_main)
            else:
                mesg= "ERROR!! Option --ctf='%s' not recognized. Allowed options: flip, wiener" % options.ctf
                print_log_msg(mesg, log, is_main)
                safe_exit()
        else:
            mesg=  "CTF flag specified, but program doesn't write images without alignment option."
            mesg+= "\nRe-run with the '--align_isac_dir=<isac_directory>' option."
            print_log_msg(mesg, log, is_main)
            safe_exit()
    else:
        if options.align_isac_dir and verbosity>=1:
            print_log_msg("Not applying CTF correction", log, is_main)
        
    centeredTF= False
    if options.write_centered:
        img= sp_utilities.get_im(classavgstack, 0)
        
        # Check that centering parameter (global CENTERING_TAG) exists 
        if not img.has_attr(CENTERING_TAG):
            mesg= "ERROR!! '%s' is missing header tag %s from sp_center_2d3d.py" % (classavgstack, CENTERING_TAG)
            print_log_msg(mesg, log, is_main)
            print_log_msg("Run sp_center_2d3d.py or remove '--write_centered' flag.", log, is_main)
            safe_exit()
        
        # Average has centering tag
        else:
            if options.align_isac_dir:
                centeredTF= True
                if verbosity>=1:
                    mesg= "Applying centering from '%s' to parameters" % classavgstack
                    print_log_msg(mesg, log, is_main)
                
                if options.applyparams:
                    if options.applyparams == 'combined' and options.debug:
                        mesg= "ERROR!! Option --debug implemented only for option --applyparams=intshifts."
                        print_log_msg(mesg, log, is_main)
                        safe_exit()
                        
                    if options.shrink != 1: 
                        mesg= "ERROR!! Option --shrink=%s not implemented with --applyparams"
                        print_log_msg(mesg, log, is_main)
                        safe_exit()
                        
                    if options.applyparams == 'combined' or options.applyparams == 'intshifts':
                        selected_params_basename= PARAMS_PREFIX + options.applyparams + '.txt'
                        if verbosity>=1:
                            mesg= "Will center using parameter file '%s'" % selected_params_basename
                            print_log_msg(mesg, log, is_main)
                    
                        if options.applyparams == 'intshifts' and verbosity>=1:
                            if options.debug:
                                mesg= "DEBUG mode: Rotation applied to integer-shifted images AFTER writing, for computing averages"
                                print_log_msg(mesg, log, is_main)
                            else:
                                mesg= "NOTE: This option does not compute averages, since rotation is not applied. (See '--debug' flag.)"
                                print_log_msg(mesg, log, is_main)
                    else:
                        mesg= "ERROR!! Option --applyparams='%s' unknown." % options.applyparams
                        mesg+= "Available options are: combined and intshifts"
                        print_log_msg(mesg, log, is_main)
                        safe_exit()
                    
                else:
                    mesg= "WARNING! Will not apply alignment parameters to stack, will simply write parameter file (see parameter '--applyparams')"
                    print_log_msg(mesg, log, is_main)
                    selected_params_basename= PARAMS_PREFIX + 'combined.txt'
                
            else:
                mesg=  "Centering specified, but program doesn't write images without alignment option."
                mesg+= "\nRe-run with the '--align_isac_dir=<isac_directory>' option."
                print_log_msg(mesg, log, is_main)
                safe_exit()
            # End alignment IF-THEN
        # End centering IF-THEN
    
    chains_dir= None
    if options.align_isac_dir:
        # Initialize mask info
        mask2count= None  
        mask_type= None
        
        # Check whether to use ordered class averages
        order_method= check_ordered_avgs(options.align_isac_dir, classavgstack)
        
        if verbosity>=1:
            if order_method == 'header':
                mesg= "Combining using ordered class-average alignment parameters stored in header"
            elif order_method == 'textfile':
                chains_dir= options.align_isac_dir
                chains_params_file= os.path.join(chains_dir, ORDERED_PARAMS)
                mesg= "Combining ordered class-average alignment parameters from '%s'" % (chains_params_file)
            else:
                mesg= "Not using existing ordered class-average alignments"
                if not options.chains_radius : mesg+= ", use '--chains_radius=<alignment_radius>' flag to force internal alignment"
        if verbosity>=1 : print_log_msg(mesg, log, is_main)
        
        if options.chains_radius:
            if usempi: 
                print_log_msg("ERROR!! Option '--chains_radius' doesn't work under MPI! Sorry!\n", log, is_main)
                safe_exit()
            
            if order_method != None:
                if order_method == 'header':
                    if verbosity>=1:
                        mesg= "WARNING! Alignment information for input averages already present. Overriding..."
                        print_log_msg(mesg, log, is_main)
            
            if options.chains_exe:
                if not os.path.exists(options.chains_exe):
                    print_log_msg("ERROR!! Path '%s' does not exist!" % options.chains_exe, log, is_main)
                    safe_exit()
                else:
                    if verbosity>=1:
                        mesg= "Using non-default executable at '%s'" % options.chains_exe
                        print_log_msg(mesg, log, is_main)

                chains_exe= 'python ' + options.chains_exe
            else:
                chains_exe= "sp_chains.py"
            
            if verbosity>=1:
                mesg= "Running internal alignment of class averages using a radius of %s pixels" % options.chains_radius
                print_log_msg(mesg, log, is_main)
                print_log_msg("WARNING: Might be slow for large stacks or image sizes\n", log, is_main)
                
            new_avgs= os.path.join(outdir, ORDERED_CLASS_AVGS)
            
            # Run sp_chains
            jnk= os.path.join(outdir, 'junk.hdf')
            cmd= '%s' % chains_exe
            args= '%s %s %s %s --circular --radius=%s --align' % (
                chains_exe, classavgstack, jnk, new_avgs, options.chains_radius
                )
            system_call_23(cmd, args, log=log, verbose=verbosity>=1)
            
            # Sanity checks
            if not os.path.exists(new_avgs):
                print_log_msg("ERROR!! sp_chains.py did not produce '%s'" % new_avgs, log, is_main)
                safe_exit()
            
            prms= os.path.join(outdir, ORDERED_PARAMS)
            num_new_classes= EMAN2.EMUtil.get_image_count(new_avgs)
            os.remove(jnk)
            
            if not num_new_classes == num_classes:
                mesg= "ERROR!! sp_chains.py produced a different number of averages: %s vs %s" % (num_new_classes, num_classes)
                print_log_msg(mesg, log, is_main)
                safe_exit()
            
            # Need for later
            chains_dir= outdir  
            classavgstack= new_avgs
            order_method= check_ordered_avgs(outdir, classavgstack)
        # End chains_radius IF-THEN
        
    if options.nvec:
        if options.write_centered:
            print_log_msg("ERROR!! Option --nvec not available with option --write_centered", log, is_main)
            print_log_msg("(sp_applications.prepare_2d_forPCA() doesn't work consistently with all options.)", log, is_main)
            safe_exit()
                        
        mesg= 'Writing %s eigenimages per class' % options.nvec
        if options.pca_radius >=1 and (options.mask_drawn or options.mask_binary):
            print()
            print_log_msg("ERROR!! Option '--pca_radius' not available with option '--mask_drawn' or '--mask_binary'", log, is_main)
            print_log_msg("Pick one and restart\n", log, is_main)
            safe_exit()
        
        elif options.mask_drawn and options.mask_binary:
            print()
            print_log_msg("ERROR!! Options '--mask_drawn' and '--mask_binary' cannot be used simulatenously", log, is_main)
            print_log_msg("Pick one and restart\n", log, is_main)
            safe_exit()
            
        elif options.pca_radius >= 1:
            mesg+= ', using a mask radius of %s pixels (based on the full-sized images)' % options.pca_radius
        
        elif options.mask_drawn:
            mask2count= options.mask_drawn
            if not os.path.exists(mask2count):
                print_log_msg("ERROR!! Mask file '%s' doesn't exist!" % mask2count, log, is_main)
                safe_exit()
            mesg+= ", using manually-drawn mask '%s'" % mask2count
            mask_type= 'drawn'
            
        elif options.mask_binary:
            mask2count= options.mask_binary
            if not os.path.exists(mask2count):
                print_log_msg("ERROR!! Mask file '%s' doesn't exist!" % mask2count, log, is_main)
                safe_exit()
            mesg+= ", using binary mask '%s'" % mask2count
            mask_type= 'binary'
            
        else:
            mesg+= ", using no mask"
        
        if verbosity>=1 : print_log_msg(mesg, log, is_main)
    # End PCA if-then
        
    # Set filenames
    classmap= os.path.join(outdir, DOCFILEDIR, CLASSMAPFILE)
    classdoc_template= os.path.join(outdir, DOCFILEDIR, CLASSDOCPREFIX + digit_pattern + '.txt')
    processed_vomq_imgs_file= os.path.join(outdir, OLD2NEW_PARTICLES)
    
    if is_main:
        if verbosity>=1 : print()
        num_vomq_classes, num_vomq_parts= vomq(
            classavgstack, 
            classmap, 
            classdoc=classdoc_template, 
            combined_parts_file=processed_vomq_imgs_file, 
            log=log, 
            verbose=verbosity>=1
            )
    else:
        num_vomq_classes= 0
        num_vomq_parts= 0
    
    # If SPIDER-like 'VO MQ' mode, we're already done
    if options.particles==None : return
    
    # Broadcast to all processes
    if usempi:
        num_vomq_classes= sp_utilities.bcast_number_to_all(num_vomq_classes, source_node= main_mpi_proc, mpi_comm= mpi.MPI_COMM_WORLD)
        num_vomq_parts=   sp_utilities.bcast_number_to_all(num_vomq_parts,   source_node= main_mpi_proc, mpi_comm= mpi.MPI_COMM_WORLD)
    
    # If not aligning images
    if not options.align_isac_dir:
        stack2split= partstack
    
    # If aligning images
    else:  # if options.align_isac_dir
        # Check whether ISAC or beautified
        if is_main:
            params_type, isac_shrink_ratio, combined_params_file, idim_parts= check_isac_or_beautify(
                processed_vomq_imgs_file, 
                options.align_isac_dir, 
                partstack, 
                classavgstack, 
                outdir, 
                processed_beautified_imgs_file= os.path.join(outdir, ISAC_PROCESSED_IMGS), 
                subsetTF=(order_method != None or centeredTF), 
                log=log, 
                verbose=verbosity>=1
                )
            
            # Resize mask, or simply print maximum if PCA won't be performed
            pca_mask= prepare_mask(
                classavgstack, 
                mask2count, 
                mask_type, 
                outstack=os.path.join(outdir, STACKFILEDIR, LOCALMASKS), 
                idim_parts=idim_parts, 
                params_type=params_type, 
                isac_shrink_ratio=isac_shrink_ratio, 
                shrink_factor=options.shrink, 
                is_main=is_main, 
                log=log, 
                verbose=verbosity>=1
                )
        else:
            params_type= ''
            isac_shrink_ratio= 0.0
            combined_params_file= ''
            pca_mask=''
            idim_parts= 0
        
        if usempi:
            params_type= sp_utilities.send_string_to_all(params_type, source_node= main_mpi_proc)
            isac_shrink_ratio= sp_utilities.bcast_number_to_all(isac_shrink_ratio, source_node= main_mpi_proc, mpi_comm= mpi.MPI_COMM_WORLD)
            combined_params_file= sp_utilities.send_string_to_all(combined_params_file, source_node= main_mpi_proc)
            pca_mask= sp_utilities.send_string_to_all(pca_mask, source_node= main_mpi_proc)
            idim_parts= sp_utilities.bcast_number_to_all(idim_parts, source_node= main_mpi_proc, mpi_comm= mpi.MPI_COMM_WORLD)
            quick_barrier()
            
            # Sanity check
            assert idim_parts!=0, "ERROR!! idim_parts %s wasn't broadcasted" % idim_parts
        
        combined_bdb_obj= BdbNames(BIGSTACKCOPY, outdir)
        if os.path.exists(combined_bdb_obj.bdb_path) : os.remove(combined_bdb_obj.bdb_path)  # will otherwise merge with pre-existing file
        
        # Need to make a subset with processed images
        if order_method != None or centeredTF:
            if is_main:
                # Make substack with processed images
                if verbosity>=1: 
                    mesg= "Making substack of original stack '%s' using subset in '%s'" % (partstack, processed_vomq_imgs_file)
                    print_log_msg(mesg, log, is_main)
                
                sane_bdb_subset(
                    partstack, 
                    combined_bdb_obj.bdb_name, 
                    processed_vomq_imgs_file, 
                    log=log, 
                    verbose=verbosity>=1
                    )
            
                if centeredTF:
                    # Write a simple list of particles.
                    # After splitting into stacks, original numbering will otherwise be lost.
                    particle_listdoc= os.path.join(outdir, CONSECPARTICLES)
                    if verbosity>=1: 
                        print_log_msg("Writing consecutive particle number to '%s'" % particle_listdoc, log, is_main)
                    sp_utilities.write_text_row(range(EMAN2.EMUtil.get_image_count(combined_bdb_obj.bdb_name) ), particle_listdoc) 
                
                    # Import original number into a new label that will be preserved (e2bdb + sanity check)
                    sane_update_header(
                        combined_bdb_obj.bdb_name, 
                        particle_listdoc, 
                        header_param=GLOBAL_NUM_TAG, 
                        log=log,
                        verbose=verbosity>=1
                        )
            # End is_main if-then
        
        # Simply copy image stack
        else:
            if verbosity>=1 : print_log_msg("Copying '%s' to virtual stack '%s'" % (partstack, combined_bdb_obj.bdb_name), log, is_main)
            cmd= "e2bdb.py %s --makevstack %s" % (partstack, combined_bdb_obj.bdb_name)
            if verbosity>=1 : print_log_msg("  " + cmd, log, is_main)
            if is_main: os.system(cmd)
            quick_barrier()
        # End subset if-then
        
        # Don't need to combine parameters if beautified and not further centered
        if params_type=='isac' or order_method != None or options.write_centered:
            if is_main:
                combine_isac_params(
                    classavgstack, 
                    options.align_isac_dir, 
                    outdir, 
                    params_type=params_type, 
                    isac_shrink_ratio=isac_shrink_ratio, 
                    beautify_params_file=combined_params_file,
                    old_combined_parts=processed_vomq_imgs_file, 
                    part_list_template=classdoc_template, 
                    chains_dir=chains_dir, 
                    centerTF=options.write_centered, 
                    log=log, 
                    verbosity=verbosity
                    )
            quick_barrier()
            
            # Beautified parameters have their own combined_params_file, which has been superceded
            combined_params_file= os.path.join(outdir, PARAMS_PREFIX + 'combined.txt')
        # End combine-parameters if-then
        
        stack2split= combined_bdb_obj.bdb_name
    # End alignment if-then
    
    # BDB where 2D parameters will be converted to 3D, might be updated depending on options
    final_bdb_name= stack2split
    
    tot_parts= 0  # initialize particle-counter
    class_init_bdb_template= os.path.join(outdir + "#" + CLASSORIGBDBPREFIX + digit_pattern)
    
    # Not centered averages is the default
    if not centeredTF:
        # If no ISAC alignments
        if not options.align_isac_dir:
            if verbosity>=1: 
                print_log_msg("Writing %s class BDB stacks" % num_vomq_classes, log, is_main)
                if options.filtrad != None or options.shrink != 1:
                    fn= os.path.join(outdir, STACKFILEDIR, CLASS_STACK_PREFIX)
                    print_log_msg("Writing stacks to '%s*%s'" % (fn, IMGFORMAT), log, is_main)
                
            class_stack_basename_template= CLASS_STACK_PREFIX + digit_pattern + IMGFORMAT
            
            separate_by_class(
                num_vomq_classes,
                stack2split, 
                classdoc_template, 
                class_init_bdb_template,
                outdir, 
                os.path.join(outdir, STACKFILEDIR), 
                class_stack_basename_template, 
                filter_radius=options.filtrad,
                shrink_factor=options.shrink,
                log=log, 
                verbosity=verbosity, 
                mpioptions=mpioptions
                )
                
            if is_main: 
                tot_parts= count_total_particles(
                    num_vomq_classes, 
                    class_init_bdb_template
                    )
            else:
                tot_parts= 0
            
            if usempi: 
                quick_barrier()
                tot_parts= sp_utilities.bcast_number_to_all(tot_parts, source_node= main_mpi_proc, mpi_comm= mpi.MPI_COMM_WORLD)
                assert tot_parts!=0, "ERROR!! tot_parts %s wasn't broadcasted" % tot_parts
                
        # Applying ISAC alignments
        else:
            # Import combined alignment parameters
            if is_main:
                sane_update_header(
                    stack2split, 
                    combined_params_file, 
                    header_param='xform.align2d', 
                    log=log,
                    verbose=verbosity>=1
                    )
            quick_barrier()
            
            # Build BDB name
            class_bdb_template= CLASS_BDB_STEM
            if options.shrink != 1 : class_bdb_template= ALIGNBDBPREFIX + class_bdb_template
            if options.filtrad : class_bdb_template= FILTCLASSBDBPREFIX + class_bdb_template
            class_bdb_template= outdir + "#" + class_bdb_template + digit_pattern
                
            if options.nvec>0:
                montage_file= EIGENIMGBDB
            else:
                montage_file= AVGVARBDB
            
            class_stack_basename_template= CLASS_STACK_PREFIX + digit_pattern + IMGFORMAT
            if verbosity>=1: 
                print_log_msg("Writing %s class BDB stacks" % num_vomq_classes, log, is_main)
                if options.align_isac_dir : print_log_msg("Writing stacks for each class", log, is_main)
                
            # Only used if phase-flipping
            pws_doc_template= os.path.join(outdir, DOCFILEDIR, PWS_1D_PREFIX + digit_pattern + '.txt')
            
            separate_by_class(
                num_vomq_classes,
                stack2split, 
                classdoc_template, 
                class_init_bdb_template,
                outdir, 
                os.path.join(outdir, STACKFILEDIR), 
                class_stack_basename_template, 
                do_align=True, 
                output_class_bdb_template=class_bdb_template,
                filter_radius=options.filtrad,
                shrink_factor=options.shrink,
                num_factors=options.nvec, 
                pca_radius=options.pca_radius, 
                pca_mask=pca_mask, 
                ctf_method=options.ctf, 
                montage_file=montage_file, 
                pws_doc_template=pws_doc_template, 
                log=log, 
                verbosity=verbosity, 
                mpioptions=mpioptions
                )
            
            if is_main: 
                tot_parts= count_total_particles(
                    num_vomq_classes, 
                    class_init_bdb_template
                    )
            else:
                tot_parts= 0
            
            if usempi: 
                quick_barrier()
                tot_parts= sp_utilities.bcast_number_to_all(tot_parts, source_node= main_mpi_proc, mpi_comm= mpi.MPI_COMM_WORLD)
                assert tot_parts!=0, "ERROR!! tot_parts %s wasn't broadcasted" % tot_parts
        # End alignment if-then
    
    # If averages were centered
    else:  # if centeredTF:
        # Read xform.align2d
        if is_main : sane_update_header(
            stack2split, 
            os.path.join(outdir, selected_params_basename), 
            header_param='xform.align2d', 
            log=log,
            verbose=verbosity>=1
            )
        quick_barrier()
        
        # If debugging, will write averages + variance, which require rotation angle, stored in a special ROTATION_TAG
        if options.debug and options.applyparams == 'intshifts':
            params_file_rotonly= os.path.join(outdir, PARAMS_PREFIX + 'rotonly.txt')
            selected_params_file= params_file_rotonly
            mesg= "DEBUG MODE: Importing rotation from '%s'" % selected_params_file
            print_log_msg(mesg, log, is_main)
            
            if is_main : sane_update_header(
                stack2split, 
                selected_params_file, 
                header_param=ROTATION_TAG, 
                log=log,
                verbose=verbosity>=1
                )
            quick_barrier()
            
        if options.applyparams:  # a subset of use of centering parameters
            # Filename
            class_stack_basename_template= CLASS_STACK_PREFIX + digit_pattern + IMGFORMAT
            montage_file= AVGVARBDB
            pws_doc_template= os.path.join(outdir, DOCFILEDIR, PWS_1D_PREFIX + digit_pattern + '.txt')
            aligned_class_bdb_template= os.path.join(outdir, ALIGNBDBPREFIX + digit_pattern)
            
            separate_by_class(
                num_vomq_classes,
                stack2split, 
                classdoc_template, 
                class_init_bdb_template,
                outdir, 
                os.path.join(outdir, STACKFILEDIR), 
                class_stack_basename_template, 
                do_align=True, 
                output_class_bdb_template=aligned_class_bdb_template,
                filter_radius=options.filtrad,
                shrink_factor=options.shrink,
                num_factors=options.nvec, 
                pca_radius=options.pca_radius, 
                pca_mask=pca_mask, 
                ctf_method=options.ctf, 
                montage_file=montage_file, 
                pws_doc_template=pws_doc_template, 
                center_params_type=options.applyparams, 
                center_debugTF=options.debug, 
                log=log, 
                verbosity=verbosity, 
                mpioptions=mpioptions
                )
            
            if is_main: 
                tot_parts= count_total_particles(
                    num_vomq_classes, 
                    class_init_bdb_template
                    )
            else:
                tot_parts= 0
            
            if usempi: 
                quick_barrier()
                tot_parts= sp_utilities.bcast_number_to_all(tot_parts, source_node= main_mpi_proc, mpi_comm= mpi.MPI_COMM_WORLD)
                assert tot_parts!=0, "ERROR!! tot_parts %s wasn't broadcasted" % tot_parts
                
            # Merge class stacks
            combined_aligned_bdb_name= "bdb:" + os.path.join(outdir, BIG_ALIGN_BDB_PREFIX + options.applyparams)
            if verbosity>=1 : print_log_msg("Merging class BDBs into '%s'" % (combined_aligned_bdb_name), log, is_main)
            
            if is_main:
                sort_merge_bdbs(
                    aligned_class_bdb_template, 
                    num_vomq_classes, 
                    outdir, 
                    BIG_ALIGN_BDB_PREFIX + options.applyparams, 
                    classdoc_template, 
                    log=log, 
                    verbosity=verbosity
                    )
            
                # Sanity check
                combined_length= EMAN2.EMUtil.get_image_count(combined_aligned_bdb_name)
                if not combined_length == tot_parts :
                    print_log_msg("WARNING! %s != %s" % (combined_length, tot_parts), log, is_main)
                
            quick_barrier()
            final_bdb_name= combined_aligned_bdb_name
        
        # If not applying shifts
        else:
            quick_barrier()
        
        # End apply-centering IF-THEN
    # End centering if-then
        
    # Convert 2D parameters to 3D in case of future restacking
    if options.align_isac_dir:
        # BDB seems to have temporary gaps when using MPI
        if mpioptions['size'] == 1:
            if verbosity>=1: 
                mesg= 'Converting 2D parameters to 3D (for possible restacking, will work better if averages are centered)'
                print_log_msg(mesg, log, is_main)
                print_log_msg('  sp_params_2D_to_3D.py %s' % final_bdb_name, log, is_main)
            sp_applications.wrapper_params_2D_to_3D(final_bdb_name)
            
            # Write "3D" parameters
            projection_params_file= os.path.join(outdir, PARAMS_PREFIX + 'projection3d.txt')
            if verbosity>=1: 
                cmd= "sp_header.py %s --params=xform.projection --export=%s" % (final_bdb_name, projection_params_file) 
                print_log_msg("  "+cmd+"\n", log, is_main)
            sp_applications.header(final_bdb_name, 'xform.projection', fexport=projection_params_file)
        
    if verbosity>=1 and is_main: 
        mesg= "Done! Separated %s particles from %s classes, outputs written to '%s'" % (num_vomq_parts, num_vomq_classes, outdir)
        print_log_msg(mesg, log, is_main)
        
        if options.align_isac_dir:
            mesg= 'To convert 2D parameters to 3D (for possible restacking, will work better if averages are centered), type:'
            print_log_msg(mesg, log, is_main)
            print_log_msg('  sp_params_2D_to_3D.py %s' % final_bdb_name, log, is_main)
    
        num_final_parts= EMAN2.EMUtil.get_image_count(final_bdb_name)
        if num_vomq_parts < num_final_parts :
            print_log_msg("To generate a substack containing the selected %s particles (out of %s) in '%s', type:" % (
                num_vomq_parts, num_final_parts, os.path.basename(classavgstack)
                ), log, is_main)
            print_log_msg('  e2bdb.py %s --makevstack <output_bdb> --list %s' % (final_bdb_name, processed_vomq_imgs_file), log, is_main)
        
        print()
    
    quick_barrier()
    
def vomq(classavgstack, classmap, classdoc=None, combined_parts_file=None, log=None, verbose=True):
    """
    Separate particles according to class assignment.
    
    Arguments:
        classavgstack : Input image stack
        classmap : Output class-to-particle lookup table. Each (long) line contains particles assigned to a class, one file for all classes
        classdoc : Filename template for output lists of particles assigned to a class, of the form 'classdoc{0:03d}.txt'
        combined_parts_file : Output combined particle list
        log : instance of Logger class
        verbose : (boolean) Whether to write additional information to screen
        
    Returns:
        class_counter : Number of classes
        part_counter : Number of particles
    """
    
    if verbose: 
        print_log_msg("Exporting class members of stack '%s'" % classavgstack, log)
        print_log_msg("  sp_header.py %s --params=members --export=%s" % (classavgstack, classmap), log)
    
    # Generate class-to-particle lookup table
    sp_applications.header(classavgstack, 'members', fexport=classmap)
    
    # Generate class-to-particle lookup table and class-selection lists
    if verbose: 
        print_log_msg("Finished exporting class members to class map '%s'\n" % classmap, log)
        print_log_msg("Splitting '%s' into class selection files" % classavgstack, log)
    
    # Initialize
    combined_parts_list= []  # unsorted
    class_counter= 0
    part_counter= 0
    
    # Loop through classes
    with open(classmap) as r:
        for class_num, line in enumerate(r.readlines() ):
            class_list= line[1:-3].split(', ')
            class_counter += 1
            part_counter += len(class_list)
            combined_parts_list+= class_list
            if classdoc : sp_utilities.write_text_row(class_list, classdoc.format(class_num) )

    # Sort particle list numerically
    combined_parts_list.sort(key=int)
    if combined_parts_file : sp_utilities.write_text_row(combined_parts_list, combined_parts_file)
    
    if verbose: 
        if classdoc:
            doc_dir= os.path.dirname(classdoc)
            print_log_msg("Wrote %s class selection files to '%s/'" % (class_counter, doc_dir), log)
        else:
            print_log_msg("Found assignments for %s classes" % class_counter, log)
        
        if combined_parts_file:
            print_log_msg("Wrote list of %s particles to '%s'" % (part_counter, combined_parts_file), log)
        else:
            print_log_msg("Found %s total particle assignments" % part_counter, log)
        
        print()
    # End verbose IF-THEN
    
    return class_counter, part_counter

def check_ordered_avgs(isac_dir, classavgstack):
    """
    Checked whether there are ordered alignment parameters for class averages.
    The two options are:
        1) Whether sp_chains wrote to xform.chains2d, 
        2) File is called ordered_class_averages.hdf, and chains.params.txt is present.
    Either option requires a 2020 version of sp_chains.py
    
    Arguments:
        isac_dir : ISAC directory
        classavgstack : class average stack file
    Returns:
        order_method : where the alignment parameters are stored: 'header', 'textfile', or None
    """
    
    # Look for params file in ISAC directory
    chains_params_file= os.path.join(isac_dir, ORDERED_PARAMS)
    
    # Will check image header
    img= sp_utilities.get_im(classavgstack, 0)
    
    if img.has_attr('xform.chains2d'):
        order_method= 'header'
    elif os.path.basename(classavgstack) == ORDERED_CLASS_AVGS and os.path.exists(chains_params_file): 
        order_method= 'textfile'
    else:
        order_method= None
        
    return order_method
    
def check_isac_or_beautify(processed_imgs_file, isac_dir, partstack, classavgstack, outdir, \
        processed_beautified_imgs_file=None, subsetTF=False, log=None, verbose=False):
    """
    Determines whether a directory has ISAC parameters or beautifier parameters.
    
    Arguments:
        processed_imgs_file : list of classified images
        isac_dir : directory to be tested
        partstack : particle stack, to confirm that image dimensions are valid
        classavgstack ; class-average stack
        outdir ; directory where outputs will be written
        processed_beautified_imgs_file : list of beautified images
        subsetTF (boolean) : whether all particles are accounted for or only the processed ones
        log : instance of Logger class
        verbose : (boolean) whether to write additional information to screen
        
    Returns:
        params_type : either 'isac' or 'beautified'
        isac_shrink_ratio : downsampling factor for class averages
        combined_params_file : name of parameter file
        idim_parts : dimension of particle images
    """

    # Initialize
    isacTF= False
    beautifiedTF= False
    params_type= None
    isac_shrink_ratio= 1.
    
    init_isac_params_file= os.path.join(isac_dir, ISAC_ALIGN_DIR, ISAC_INIT_PARAMS)
    all_isac_params_file= os.path.join(isac_dir, ISAC_FINAL_PARAMS)
    isac_shrink_path= os.path.join(isac_dir, ISAC_SHRINK_RATIO)
    if processed_beautified_imgs_file == None:
        processed_beautified_imgs_file= os.path.join(outdir, ISAC_PROCESSED_IMGS)
    
    # Check whether necessary files exist
    classed_imglist_exists= os.path.exists(processed_imgs_file)
    init_params_exists= os.path.exists(init_isac_params_file)
    all_params_exists= os.path.exists(all_isac_params_file)
    isac_shrink_exists= os.path.exists(isac_shrink_path)
    
    # If all three files exist, then proceed
    if classed_imglist_exists and init_params_exists and all_params_exists and isac_shrink_exists:
        isacTF= True
        params_type='isac'
        combined_params_file= os.path.join(outdir, PARAMS_PREFIX + 'combined.txt')
        
        isac_shrink_file= open(isac_shrink_path, "r")
        isac_shrink_lines= isac_shrink_file.readlines()
        isac_shrink_ratio= float(isac_shrink_lines[5])
        
    # Check for beautifier files
    beautify_final_file= os.path.join(isac_dir, BEAUT_FINAL_PARAMS)
    beautify_init_file= os.path.join(isac_dir, BEAUT_INIT_PARAMS)
    
    beautify_final_exists=  os.path.exists(beautify_final_file)
    beautify_init_exists= os.path.exists(beautify_init_file)
    
    # Check image dimensions
    idim_parts= sp_utilities.get_im(partstack, 0)['nx']
    idim_classavg= sp_utilities.get_im(classavgstack, 0)['nx']
        
    if beautify_final_exists or beautify_init_exists:
        # Make sure that the dimensions are the same
        if idim_parts == idim_classavg:
            beautifiedTF= True 
            params_type='beautified'
            
            # Final parameters would be better, but can use initial parameters
            if beautify_final_exists:
                params2import_file= beautify_final_file
            else:
                params2import_file= beautify_init_file
            params2import_list= sp_utilities.read_text_row(params2import_file)
            num_import_parts= len(params2import_list)
            
            # Need to remove first column of parameter file to make selection file
            processed_img_list= [row[0] for row in params2import_list]
            sp_utilities.write_text_row(processed_img_list, processed_beautified_imgs_file)
            
            # If using ordered class averages, only processed images should be represented in the parameters doc
            if subsetTF:
                combined_params_file= os.path.join(outdir, PARAMS_PREFIX + "selected.txt")
                
                # Read 2nd column onward from beautifier parameter file
                beautified_params_list= [row[1:] for row in params2import_list]
            
                # Write parameters file of selected particles
                sel_part_list= sp_utilities.read_text_row(processed_imgs_file)
                sel_params_list= []
                
                # Loop through particles (sel_part_list is a list of lists of length 1)
                for newidx, partnum in enumerate(sel_part_list):
                    try:
                        intermed_idx= processed_img_list.index(partnum[0])
                        part_params= beautified_params_list[intermed_idx]
                    except IndexError:
                        print('last processed_img_list', processed_img_list[-1], type(processed_img_list[-1]) )
                        print('newidx %s, partnum %s, beautified_params_list %s' % (newidx, partnum[0], len(beautified_params_list) ) )
                        exit()
                    sel_params_list.append(part_params)
                
            # Beautifier parameter doc includes particle number & skips not-processed images, contrary to ISAC's
            else:
                combined_params_file= os.path.join(outdir, PARAMS_PREFIX + "beautified.txt")
                
                good_part_counter= 0
                num_tot_parts= EMAN2.EMUtil.get_image_count(partstack)
                sel_params_list= []
                
                for global_idx in range(num_tot_parts):
                    # Don't let program proceed part the end of the file.
                    if good_part_counter < num_import_parts : 
                        try:
                            row= params2import_list[good_part_counter]
                        except IndexError:
                            print("%s: good_part_counter %s, num_tot_parts %s, params2import_file %s" % (
                                global_idx, good_part_counter, num_tot_parts, params2import_file
                                ) )
                            exit()
                        
                        # Look for skipped particles
                        if global_idx == processed_img_list[good_part_counter]: 
                            sel_params_list.append(row[1:])
                            good_part_counter+= 1
                        else:
                            sel_params_list.append([0.,0.,0.,-1])
                    else:
                        sel_params_list.append([0.,0.,0.,-1])
                        
            sp_utilities.write_text_row(sel_params_list, combined_params_file)
        
        else:
            err= "ERROR!! Input image size %s is different from class average size %s" % (idim_parts, idim_classavg)
            print_log_msg(err, log)
            print_log_msg('Fix and re-start', log)
            exit()
        # End image-size IF-THEN
    # End beautify IF-THEN
    
    # Neither ISAC nor beautified
    if not isacTF and not beautifiedTF:
        mesg=  "\nERROR!! Directory '%s' is missing files that ISAC or beautified directories should have!" % options.align_isac_dir
        
        mesg+= "\n\nISAC directory should have:"
        mesg+= "\n  %s %s" % (processed_imgs_file, classed_imglist_exists)
        mesg+= "\n  %s %s" % (init_isac_params_file, init_params_exists)
        mesg+= "\n  %s %s" % (all_isac_params_file, all_params_exists)
        mesg+= "\n  %s %s" % (isac_shrink_path, isac_shrink_exists)
        
        mesg+= "\n\nBeautifier directory should have at least one of:"
        mesg+= "\n  %s %s" % (beautify_init_file, beautify_init_exists)
        mesg+= "\n  %s %s" % (beautify_final_file, beautify_final_exists)
        
        mesg+= "\n\nExiting..."

        print(mesg)
        exit()
    
    # Both types of files are present
    elif isacTF and beautifiedTF:
        mesg= "WARNING!! Directory '%s' has files of both ISAC or beautified directories" % options.align_isac_dir
        print_log_msg(mesg, log)
        print_log_msg('Will proceed without applying alignments...', log)
        options.align_isac_dir= None
        
    # One or the other
    elif isacTF and not beautifiedTF:
        combined_params_file= os.path.join(outdir, PARAMS_PREFIX + 'combined.txt')
        
        if verbose:
            print_log_msg("Found necessary files in ISAC directory '%s/'" % isac_dir, log)
            
            mesg= "Combining alignment parameters from '%s' and '%s', dividing by %s" % \
                (os.path.join(ISAC_ALIGN_DIR, ISAC_INIT_PARAMS), ISAC_FINAL_PARAMS, isac_shrink_ratio)
            print_log_msg(mesg, log)
            
            print_log_msg("Will write combined alignment parameters to '%s'" % combined_params_file, log)
    
    elif beautifiedTF and not isacTF:
        if verbose:
            mesg= "Found necessary files in beautifier directory '%s/'\n" % isac_dir
            print_log_msg(mesg, log)
            mesg= "Wrote list of %s particles to '%s'" % (len(processed_img_list), processed_beautified_imgs_file)
            print_log_msg(mesg, log)
            mesg= "Wrote %s sets of beautified parameters to '%s'" % (len(sel_params_list), combined_params_file)
            print_log_msg(mesg, log)
        
    else:
        print("ERROR!! Unknown state %s %s" % (isacTF, beautifiedTF) )
        exit()
    
    print()
    
    return params_type, isac_shrink_ratio, combined_params_file, idim_parts
    
def prepare_mask(classavgstack, mask_file, mask_type, outstack=LOCALMASKS, idim_parts=0, params_type=None, 
                 isac_shrink_ratio=1.0, shrink_factor=1, is_main=True, log=None, verbose=False):
    """
    Prepare mask for principal component analysis.
    
    Arguments:
        log : instance of Logger class
        verbosity : How much information to write to screen
    """
    
    num_classes= EMAN2.EMUtil.get_image_count(classavgstack)
    
    if mask_file:
        # Sanity check
        assert idim_parts>0, "\nERROR!! Need dimensions for particle images"
        assert params_type!=None, "\nERROR!! Need ISAC vs. beautifier directory structure"
        num_mask_imgs= EMAN2.EMUtil.get_image_count(mask_file)
        
        # If particle images and mask images have the same dimension
        idim_mask= sp_utilities.get_im(mask_file)['nx']
        idim_classavg= sp_utilities.get_im(classavgstack)['nx']        
        
        if verbose:
            if num_mask_imgs==1:
                print_log_msg("Mask file '%s' has single image, will use same mask for all classes" % mask_file, log, is_main)
            else:
                if num_classes!=num_mask_imgs:
                    print_log_msg("ERROR!! Mask file '%s' has %s images" % (mask_file, num_mask_imgs), log, is_main)
                    mesg= "Mask (n=%s) must match number of class averages (%s) or must use the same mask for all classes" % (num_mask_imgs, num_classes)
                    print_log_msg(mesg, log, is_main)
                    safe_exit()
            
            if idim_mask==idim_parts:
                print_log_msg("Mask dimension is same as particle dimension: %s" % idim_mask, log, is_main)
                
                # Print warning if class averages have a different dimension
                if idim_mask!=idim_classavg:
                    print()
                    print_log_msg("WARNING! Mask dimension %s is not same as class-average dimension %s" % (idim_mask, idim_classavg), log, is_main)
                    print_log_msg("Binarization threshold for masked images may not be appropriate\n", log, is_main)
                
                # Print warning if not beautified parameters
                if idim_parts==idim_classavg and params_type!='beautified':
                    print()
                    print_log_msg("WARNING! Averages appear to be full-sized, but parameters do not correspond to beautifier\n", log, is_main)
            
            # If particle images and mask images have different dimensions, will rescale mask
            else:
                assert params_type=='isac', "\nERROR!! Was expecting parameters from ISAC-compatible directory"
                mesg= "Mask dimension %s is different from particle-image dimension %s, will divide by %s" % (idim_mask, idim_parts, isac_shrink_ratio)
                print_log_msg(mesg, log, is_main)
                
            if shrink_factor>1:
                print_log_msg("Will downsample mask by a factor of %s" % shrink_factor, log, is_main)
            
            print_log_msg("Preparing mask for each class....", log, is_main)
        # End verbose if-then
            
        upscaled_dim, upscale_rate= get_rescaled_dimension(mask_file, isac_shrink_ratio)
        
        # Loop through classes
        for class_num in range(num_classes):
            if num_mask_imgs==1:
                mask_img= sp_utilities.get_im(mask_file)
            elif num_mask_imgs==num_classes:
                mask_img= sp_utilities.get_im(mask_file, class_num)
            else:
                print_log_msg("ERROR!! Mask file '%s' has %s images" % (mask_file, num_mask_imgs), log_is_main)
                mesg= "Mask must match number of class averages (%s) or must use the same mask for all classes" % (num_classes, num_mask_imgs)
                print_log_msg(mesg, log_is_main)
                safe_exit()
                
            # If mask is small, enlarge it
            mask_upscaled= sp_fundamentals.resample(mask_img, upscale_rate)
            
            # Pad or crop as necessary
            if upscaled_dim<idim_parts:
                mask_part_sized= sp_utilities.pad(mask_upscaled, idim_parts, idim_parts, 1, background="circumference")
            elif upscaled_dim>idim_parts:
                mask_part_sized= EMAN2.Util.window(mask_upscaled, idim_parts, idim_parts, 1)
            else:
                mask_part_sized= mask_upscaled
            
            # Apply optional downsampling factor
            if shrink_factor!=1:
                shrunk_dim, shrink_rate= get_rescaled_dimension(mask_part_sized, shrink_factor)
                mask_shrunken= sp_fundamentals.resample(mask_part_sized, shrink_rate)
            else:
                mask_shrunken= mask_part_sized
            
            # Binarize mask
            if mask_type=='drawn':
                classavg_max= sp_utilities.get_im(classavgstack, class_num)['maximum']
                mask_thresh= sp_morphology.binarize(mask_shrunken, minval=classavg_max)
                # (Threshold is assumed to be greater than the original maximum of the class average.)
                
                # Sanity check
                if mask_thresh['maximum'] == 0.0:
                    print()
                    print_log_msg("ERROR!! Mask for class #%s is empty" % class_num, log, is_main)
                    print_log_msg("Maximum in '%s' is %s (after resampling)," % (
                        mask_file, mask_shrunken['maximum']
                        ), log, is_main)
                    print_log_msg("  whereas maximum in '%s' used for threshold was %s" % (
                        classavgstack, classavg_max
                        ), log, is_main)
                    print_log_msg("Fix mask and try again.\n", log, is_main)
                    safe_exit()
                
            elif mask_type=='binary':
                # Maybe the once-binary mask is no longer binary if it has been interpolated
                mask_thresh= sp_morphology.binarize(mask_shrunken, minval=0.5)
            else:
                print_log_msg("ERROR!! Mask type '%s' unknown!" % mask_type, log_is_main)
                safe_exit()
            
            # Write image
            mask_thresh.write_image(outstack, class_num)
        # End class loop
        
        if verbose : print_log_msg("Wrote %s mask images to '%s'\n" % (num_classes, outstack), log, is_main)
    
    # If not using custom mask, simply write out statistics
    else:
        outstack= ''
        if verbose: 
            for class_num in range(num_classes):
                classavg_max= sp_utilities.get_im(classavgstack, class_num)['maximum']
                
                if class_num==0:
                    img_max= classavg_max
                else:
                    if classavg_max>img_max : img_max= classavg_max
            # End class loop
            
            print_log_msg("Maximum pixel value in class averages: %.4f" % img_max, log, is_main)
            print_log_msg("If you draw a custom mask using e2display.py, use a pen value higher than this number.\n", log, is_main)
    # End drawn-mask if-then
    
    return outstack
        
def combine_isac_params(classavgstack, isac_dir, outdir, params_type='isac', isac_shrink_ratio=1.0, 
                beautify_params_file=None, old_combined_parts=None, part_list_template=None, 
                chains_dir=None, centerTF=False, log=None, verbosity=0):
    """
    Combines initial and final params from ISAC.
    
    Arguments:
        classavgstack : Input image stack
        isac_dir : Input ISAC directory
        outdir : Output directory
        params_type : ISAC or Beautified
        isac_shrink_ratio : Ratio of downsampled ISAC size to original size (<=1.0)
        chains_dir : Directory to look for chains_params.txt
        old_combined_parts : List of images (e.g., processed_images.txt)
        part_list_template : Particle-membership doc file
        centerTF : (boolean) Apply centering from sp_center_2d3d
        log : instance of Logger class
        verbosity : How much information to write to screen
    """
    
    # Check whether directory is ISAC or Beautified
    
    # File-handling
    outparams_prefix= os.path.join(outdir, PARAMS_PREFIX)
    combined_params_file= outparams_prefix + 'combined.txt'
    params_file_shiftonly= outparams_prefix + 'intshifts.txt'
    params_file_rotonly= outparams_prefix + 'rotonly.txt'
    center_docfile= os.path.join(outdir, CLASS_CENTER_DOC)
    
    if verbosity>=1 : print_log_msg("Collecting alignment parameters", log)
        
    if params_type == 'isac':
        init_params_file= os.path.join(isac_dir, ISAC_ALIGN_DIR, ISAC_INIT_PARAMS)
        init_params_list= sp_utilities.read_text_row(init_params_file)
        
        all_params_file= os.path.join(isac_dir, ISAC_FINAL_PARAMS)
        all_params_list= sp_utilities.read_text_row(all_params_file)
    elif params_type == 'beautified':
        beautify_params_list= sp_utilities.read_text_row(beautify_params_file)
    else:
        print("ERROR!! Unknown params_type %s" % params_type)
        exit()
    
    # Check if using ordered class-average alignment parameters
    if chains_dir: 
        order_method= check_ordered_avgs(chains_dir, classavgstack)
    else:
        order_method= check_ordered_avgs(isac_dir, classavgstack)
    
    if order_method == 'textfile': 
        if chains_dir == None : chains_dir= isac_dir
        chains_params_file= os.path.join(chains_dir, ORDERED_PARAMS)
        chains_params_list= sp_utilities.read_text_row(chains_params_file)
    
    # Check if class-centering file will be used (TODO: Move to function)
    if centerTF: 
        # Check that centering parameter (global CENTERING_TAG) exists 
        img= sp_utilities.get_im(classavgstack, 0)
        if not img.has_attr(CENTERING_TAG):
            mesg= "ERROR!! '%s' is missing header tag %s from sp_center_2d3d.py" % (classavgstack, CENTERING_TAG)
            print_log_msg(mesg, log)
            print_log_msg("Run sp_center_2d3d.py or remove '--write_centered' flag.", log)
            exit()
        
        centerTF= True
        center_param_list= []  # will write centering parameters to a file
        
        if verbosity>=1: 
            mesg= "Exporting centering parameters from stack '%s' to '%s'" % (classavgstack, center_docfile)
            print_log_msg(mesg, log)
    
    # TODO: Using single source for parameters (i.e., not combining), e.g., from sp_mref_ali2d.py
    
    # SCENARIO #1: Simply applying ISAC alignments and not centering.
    # We can simply loop through the particles all at once and not class-wise.
    if not centerTF and not order_method:
        combined_params= []
        disableTF= verbosity>2 or verbosity<2
        
        # Loop through images
        for img in tqdm.tqdm(range(len(all_params_list) ), unit='part', disable=disableTF, file=sys.stdout): 
            # Get final ISAC parameters
            P= sp_utilities.combine_params2(
                init_params_list[img][0], 
                init_params_list[img][1], 
                init_params_list[img][2], 
                init_params_list[img][3],
                all_params_list[img][0], 
                old_div(all_params_list[img][1], isac_shrink_ratio), 
                old_div(all_params_list[img][2], isac_shrink_ratio), 
                all_params_list[img][3] 
                )
            
            combined_params.append([P[0], P[1], P[2], P[3], 1.0])
    # End simple if-then
    
    # SCENARIO #2: Using ordered_class_averages (w/chain_params) and/or class-centering
    else:
        old_combined_list= sp_utilities.read_text_row(old_combined_parts)
        
        num_classes= EMAN2.EMUtil.get_image_count(classavgstack)
        unsorted_combined= []
        
        if centerTF:
            unsort_shiftonly= []
            unsort_rotonly= []
        
        # Debug tools
        debug_class_num= None  # 13  # 4  # 
        if debug_class_num:
            first= debug_class_num
            last= first + 1
            disableTF= True
        else:
            first= 0
            last= num_classes
            
            # Turn off progress bar unless verbosity 2
            disableTF= verbosity>2 or verbosity<2
        # End debug IF-THEN
        
        # Loop through classes
        if verbosity>=1 : print_log_msg('Combining parameters for each class', log)
        for class_num in tqdm.tqdm(range(first, last), unit='class', disable=disableTF, file=sys.stdout):
            # Extract members
            classavg= sp_utilities.get_im(classavgstack, class_num)
            members= sorted(classavg.get_attr("members") )
            old_class_list= sp_utilities.read_text_row(part_list_template.format(class_num) )
            new_class_list= []
            
            # Get optional centering parameters
            if centerTF: 
                centered_avg= sp_utilities.get_im(classavgstack, class_num)
                d= centered_avg[CENTERING_TAG].get_params("2D")
                
                if verbosity>=3:
                    mesg= "  Class #%s: %16.6f %16.6f %16.6f %10d %10.3f "%(class_num, d["alpha"],d["tx"],d["ty"],d["mirror"],d["scale"])
                    print_log_msg(mesg, log)
                
                classctr_params= [ d["alpha"],d["tx"],d["ty"],d["mirror"],d["scale"] ]
                center_param_list.append(classctr_params)
                
            if order_method != None:
                if order_method == 'textfile': 
                    class_chains_params= chains_params_list[class_num][2:]  # first two columns are indices
                if order_method == 'header':
                    img= sp_utilities.get_im(classavgstack, class_num)
                    class_chains_params= list( sp_utilities.get_params2D(img, 'xform.chains2d') )
            
            # Loop through particles in class
            for idx, partnum in enumerate(members):
                if params_type == 'isac':
                    cmbprm= sp_utilities.combine_params2(
                        init_params_list[partnum][0], 
                        init_params_list[partnum][1], 
                        init_params_list[partnum][2], 
                        init_params_list[partnum][3],
                        all_params_list[partnum][0], 
                        old_div(all_params_list[partnum][1], isac_shrink_ratio), 
                        old_div(all_params_list[partnum][2], isac_shrink_ratio), 
                        all_params_list[partnum][3] 
                        )
                else:
                    try:
                        new_part_num= old_combined_list.index([partnum])
                        cmbprm= beautify_params_list[new_part_num]
                    except IndexError:
                        print('\npartnum %s, len(beautify_params_list) %s, old_combined_parts %s' % (
                            partnum, len(beautify_params_list), old_combined_parts
                            ) )
                        new_part_num= old_combined_list.index([partnum])
                        print('\npartnum %s, new_part_num %s, old_combined_list[new_part_num] %s' % (
                            partnum, new_part_num, old_combined_list[new_part_num]
                            ) )
                        exit()
                
                if debug_class_num : isacprm= cmbprm
                
                # If ordered_class_averages w/chains_params
                if order_method != None:
                    cmbprm= sp_utilities.combine_params2(
                            cmbprm[0], 
                            cmbprm[1], 
                            cmbprm[2], 
                            cmbprm[3], 
                            class_chains_params[0], 
                            class_chains_params[1], 
                            class_chains_params[2], 
                            class_chains_params[3]
                            )
                
                # Combine with class-centering parameters
                if centerTF:
                    cmbprm= sp_utilities.combine_params2(
                        cmbprm[0], 
                        cmbprm[1], 
                        cmbprm[2], 
                        cmbprm[3], 
                        classctr_params[0], 
                        classctr_params[1], 
                        classctr_params[2], 
                        classctr_params[3]
                        )
                    
                    # Now, "undo" the rotation
                    if cmbprm[3] == 0:
                        unrot_params= sp_utilities.combine_params2(cmbprm[0],cmbprm[1],cmbprm[2],cmbprm[3], -cmbprm[0],0,0,0)
                        unsort_rotonly.append([partnum, cmbprm[0], 0.,0.,0])
                    else:
                        unrot_params= sp_utilities.combine_params2(cmbprm[0],cmbprm[1],cmbprm[2],cmbprm[3], cmbprm[0],0,0,0)
                        unsort_rotonly.append([partnum, -cmbprm[0], 0.,0.,0])
                    
                    # Sanity check for angle, and round off shifts
                    assert abs( math.sin( math.radians(unrot_params[0]) ) ) < 0.01, "MRK_DEBUG: Non-zero rotation angle"
                    
                    shiftx= round(unrot_params[1])
                    shifty= round(unrot_params[2])
                    
                    unsort_shiftonly.append([ partnum, 0.0, shiftx, shifty, unrot_params[3] ])
                
                # Pre-pend the particle number so that we can sort later
                unsorted_combined.append([ partnum, cmbprm[0], cmbprm[1], cmbprm[2], cmbprm[3] ])
                
                # Need to update class number in class docs
                try:
                    old_part_num= old_class_list[idx]
                except IndexError:
                    print('\npart_list_template.format(class_num)', part_list_template.format(class_num))
                    printvars(['class_num', 'idx'], quitTF=True)
                
                try:
                    new_part_num= old_combined_list.index(old_part_num)
                except ValueError:
                    print("Couldn't find particle: class_num %s, old_part_num %s, new_part_num %s" % (class_num, old_part_num[0], new_part_num) )
                
                new_class_list.append(new_part_num)
                
                if debug_class_num:
                    print(class_num, idx, old_part_num, new_part_num, 
                            isacprm[0], isacprm[1], isacprm[2], isacprm[3], 
                            cmbprm[0], cmbprm[1], cmbprm[2], cmbprm[3])
            # End particle-loop
        
            # Overwrite pre-existing class doc
            sp_utilities.write_text_row(new_class_list, part_list_template.format(class_num) )  
        # End class-loop
        
        if debug_class_num : exit()
        
        # Sort by particle number
        combined_params= sorted( unsorted_combined, key=operator.itemgetter(0) )
        if centerTF:
            shiftonly_params= sorted( unsort_shiftonly, key=operator.itemgetter(0) )
            rotonly_params= sorted( unsort_rotonly, key=operator.itemgetter(0) )
            
        # Remove first column
        for row in combined_params: del row[0]
        if centerTF:
            for row in shiftonly_params: del row[0]
            for row in rotonly_params: del row[0]
    # End simple/combined if-then
    
    sp_utilities.write_text_row(combined_params, combined_params_file)
    if verbosity>=1 : print_log_msg("Wrote %s sets of alignment parameters to '%s'" % (len(combined_params), combined_params_file), log)

    if centerTF:
        sp_utilities.write_text_row(center_param_list, center_docfile)
        
        if params_file_shiftonly:
            sp_utilities.write_text_row(shiftonly_params, params_file_shiftonly)
            if verbosity>=1 : print_log_msg("Wrote %s integer shifts to '%s'" % (len(shiftonly_params), params_file_shiftonly), log)
        
        if params_file_rotonly:
            sp_utilities.write_text_row(rotonly_params, params_file_rotonly)
            if verbosity>=1 : print_log_msg("Wrote %s rotation angles to '%s'" % (len(rotonly_params), params_file_rotonly), log)
            
    if verbosity>=1 : print()
        
def sane_update_header(stack_file, params_file, header_param='xform.align2d', log=None, verbose=False):
    """
    Imports header information after performing a sanity check.
    
    Arguments:
        stack_file stack file to update
        params_file : parameter file to import
        header_param='xform.align2d'
        log : instance of Logger class
        verbose : (boolean) Whether to write additional information to screen
    """

    num_stack_parts= EMAN2.EMUtil.get_image_count(stack_file)
    num_params_parts= len( sp_utilities.read_text_row(params_file) )
    assert num_params_parts == num_stack_parts, "ERROR!! Mismatch: len(%s) %s != len(%s) %s" % (
            params_file, num_params_parts, stack_file, num_stack_parts
            )
    
    # Import combined alignment parameters
    if verbose: 
        mesg= "Importing combined parameters '%s'" % params_file
        print_log_msg(mesg, log)
        cmd= "sp_header.py %s --params=%s --import=%s" % (stack_file, header_param, params_file) 
        print_log_msg("  " + cmd, log)
    sp_applications.header(stack_file, header_param, fimport=params_file)
    
    if verbose : print_log_msg("Finished importing parameters into '%s'\n" % stack_file, log)
        
def sane_bdb_subset(input_bdb_stack, output_bdb_stack, particle_list, prepend="", log=None, verbose=False):
    """
    Makes a subset of a BDB stack after performing a sanity check.
    
    Arguments:
        input_bdb_stack : should include preceding 'bdb:'
        output_bdb_stack : should include preceding 'bdb:'
        particle_list : selection file
        prepend : string to prepend to printed command
        log : instance of Logger class
        verbose : (boolean) Whether to write additional information to screen
    
    Returns:
        number of images in output stack
    """
    
    # Read subset list
    last_particle= sp_utilities.read_text_row(particle_list)[-1][0]  # read_text_row row generates a list of lists
    
    # Check whether list extend past length of stack
    num_bdb_parts= EMAN2.EMUtil.get_image_count(input_bdb_stack)
    assert last_particle <= num_bdb_parts, "ERROR!! %s specifies particle number %s beyond length of %s (%s)" % (
            particle_list, last_particle, input_bdb_stack, num_bdb_parts
            )
    
    cmd= "e2bdb.py %s --makevstack %s --list %s" % (input_bdb_stack, output_bdb_stack, particle_list)
    if verbose: 
        print_log_msg(prepend + cmd, log)
    else:
        cmd+= ' >> ' + os.devnull
    os.system(cmd)
    num_class_imgs= EMAN2.EMUtil.get_image_count(output_bdb_stack)
    
    return num_class_imgs
    
def separate_by_class(num_classes, stack2split, classdoc_template, class_init_bdb_template, 
                outdir, stack_dir, class_stack_basename_template, do_align=False, 
                output_class_bdb_template=None, filter_radius=None, shrink_factor=1, num_factors=0, 
                pca_radius=-1, pca_mask="", ctf_method=None, montage_file=None, pws_doc_template=None, 
                center_params_type=None, center_debugTF=False, log=None, verbosity=0, mpioptions=None):
    """
    Separates particles by class.
    
    Arguments:
        num_classes : number of classes
        stack2split : input BDB stack which will be separated
        classdoc_template : template for class selection filenames
        class_init_bdb_template : template for BDB stacks for each class, otherwise unchanged from original stack
        outdir : output directory
        stack_dir : directory where MRC stacks will be written
        class_stack_basename_template : template for MRC stack filenames
        do_align : (boolean) whether to apply aligments
        output_class_bdb_template : template for BDBs after alignment/filtration/downsampling/etc.
        filter_radius : low-pass filter radius
        shrink_factor : downsampling factor
        num_factors : number of eigenimages to calculate
        pca_radius : radius for principal component analysis
        pca_mask : mask file for principal component analysis
        ctf_method : method for CTF correctioni: None, 'flip', 'wiener'
        montage_file : output stack of average + variance + eigenimages
        pws_doc_template : if phase-flipping, 1D profiles of amplitudes before & after CTF-correction
        center_params_type : method for centering: None, 'intshifts', 'combined'
        center_debugTF : debug flag when using 'intshifts' above, rotates before shifting by integer values
        log : instance of Logger class
        verbose : (boolean) Whether to write additional information to screen
        mpioptions : (dict) MPI options
    """
    
    # For subprocesses
    python_version, devnull= setup_subprocess()
    
    tot_parts= 0  # initialize particle-counter
    
    if pca_radius > 0:
        pca_radius= pca_radius//shrink_factor
    else:
        pca_radius= -1
    
    # Parallelization stuff
    is_main= mpioptions['is_main']
    myid= mpioptions['myid']
    usempi= mpioptions['use']
    main_mpi_proc= mpioptions['main_process']
    num_proc= mpioptions['size']
    mpi_start, mpi_end= sp_applications.MPI_start_end(num_classes, mpioptions['size'], mpioptions['myid'] )
    
    if verbosity>=2:
        if usempi:
            print_log_msg("Splitting %s classes into %s chunks of about %s classes each" % (
                num_classes, mpioptions['size'], mpi_end-mpi_start
                ), log, is_main)
        else:
            print_log_msg("Separating '%s' into %s classes" % (stack2split, num_classes), log, is_main)
    
    # Turn off progress bar unless verbosity 2
    disable_class= verbosity>2 or verbosity<2 or not is_main
    
    # Get image dimension (TODO: Bcast using MPI)
    shrunk_size, _= get_rescaled_dimension(stack2split, shrink_factor)
    montage_xdim= shrunk_size*(2 + num_factors)
    
    if do_align:
        mpi_digits= len( str(num_proc) )
        digit_pattern= '{{0:0{0}d}}'.format(mpi_digits)
        mpi_ave_stack_path= os.path.join(stack_dir, 'stk' + MPI_STACK_PREFIX + digit_pattern.format(myid) + '.mrcs')
        mpi_ave_stack_obj= EMAN2.EMData(montage_xdim, shrunk_size, (mpi_end-mpi_start) )
        mpi_ave_bdb_template= outdir + "/" + MPI_STACK_PREFIX + digit_pattern  # The '/' will make makerelpath easier to use
        mpi_ave_bdb_name= "bdb:" + mpi_ave_bdb_template.format(myid)
        mpi_ave_bdb_dict= EMAN2db.db_open_dict(mpi_ave_bdb_name)
    
    for class_num in tqdm.tqdm(range(mpi_start, mpi_end), unit='class', disable=disable_class, file=sys.stdout):
        class_mpi_idx= class_num - mpi_start
        
        class_init_bdb_name= 'bdb:' + class_init_bdb_template.format(class_num)
        classdoc_current= classdoc_template.format(class_num)
        
        # Write class stack (e2bdb + sanity check)
        tot_parts+= sane_bdb_subset(
            stack2split, 
            class_init_bdb_name, 
            classdoc_current, 
            prepend="  ", 
            log=log, 
            verbose=verbosity>=4
            )
        
        # If not applying ISAC alignments
        if not do_align:
            # Optional filtered stack
            if filter_radius != None or shrink_factor != 1:
                output_class_stack= os.path.join(stack_dir, class_stack_basename_template.format(class_num))
                
                # So that EMAN output doesn't break up progress bar
                if verbosity<=2:
                    stdout= devnull
                else:
                    stdout=None
                
                cmd= "e2proc2d.py"
                args= "%s %s --inplace" % (class_init_bdb_name, output_class_stack)
                # --inplace overwrites existing images
                
                if filter_radius : args= args + " --process=filter.lowpass.gauss:cutoff_freq=%s" % filter_radius
                if shrink_factor != 1 : args= args + " --meanshrink=%s" % shrink_factor
                
                system_call_23(cmd, args, stdout=stdout, log=log, verbose=verbosity>=4)
        
        # If applying ISAC alignments:
        # If we're not computing eigenimages and not centering, we: 1) compute averages, and 2) apply alignments.
        # If we are computing eigenimages, we: 1) filter the stack, 2) compute averages, and 3) apply alignments.
        # If we are centering, we 1) apply alignment, 2) compute averages
        # Computing eigenimages and centering is not implemented yet
        else: # if do_align:
            # Set filenames
            output_class_stack= os.path.join(stack_dir, class_stack_basename_template.format(class_num))
            aligned_bdb_name= "bdb:" + output_class_bdb_template.format(class_num)
            
            if verbosity>=3:
                if is_main : print_log_msg("  Working on class '%s'" % class_init_bdb_name, log, is_main)
            
            # If not centering
            if not center_params_type:
                avg_var_eig_list= []
                
                # If computing eigenimages (and not centering)
                if num_factors != 0:
                    assert center_params_type == None, "ERROR!! Option '--nvec' not compatible with option '--applyparams'"
                    
                    # Low-pass filter and shrink images
                    num_class_imgs= align_filter_shrink(
                        class_init_bdb_name, 
                        output_class_stack, 
                        aligned_bdb_name, 
                        alignYN=False, 
                        filter_radius=filter_radius, 
                        shrink_factor=shrink_factor, 
                        usectf=ctf_method, 
                        method='quadratic', 
                        log=log, 
                        verbosity=0, 
                        is_main=is_main
                        )
                    # (Will re-run align_filter_shrink below, so don't repeat screen output.)

                    tmp_classavg, class_stack_list= prepare_2d_forPCA(aligned_bdb_name, mode='a', CTF=False)
                    # prepare_2d_forPCA w/CTF uses a Wiener filter, and 
                    # prepare_2d_forPCA w/o CTF always applies alignment regardless of 'mode' parameter.
                    # So we need to filter+shrink the images twice, first without alignment & later with alignment.
                    
                    if pca_mask!="":
                        assert pca_radius<=0, "ERROR!! Mask radius (%s) and mask file (%s) are mutually exclusive for PCA!" % (pca_radius, pca_mask)
                        
                        curr_mask= sp_utilities.get_im(pca_mask, class_num)
                        avg_var_eig_list= pca(class_stack_list, nvec=num_factors, maskfile=curr_mask)
                    else:
                        avg_var_eig_list= pca(class_stack_list, nvec=num_factors, mask_radius=pca_radius)
                    
                    ####avg_img, var_img= avg_optional_ctf(
                        ####aligned_bdb_name, 
                        ####ctf_method=ctf_method, 
                        ####pws_docfile=pws_doc_template.format(class_num), 
                        ####do_align=True
                        ####)
                # End PCA if-then
                    
                # Apply final alignments
                num_class_imgs= align_filter_shrink(
                    class_init_bdb_name, 
                    output_class_stack, 
                    aligned_bdb_name, 
                    alignYN=True, 
                    filter_radius=filter_radius, 
                    shrink_factor=shrink_factor, 
                    usectf=ctf_method, 
                    method='quadratic', 
                    log=log, 
                    verbosity=verbosity, 
                    is_main=is_main
                    )
                
                avg_img, var_img= avg_optional_ctf(
                    aligned_bdb_name, 
                    ctf_method=ctf_method, 
                    pws_docfile=pws_doc_template.format(class_num), 
                    do_align=False
                    )
                num_class_imgs= EMAN2.EMUtil.get_image_count(class_init_bdb_name)
                
                image_list= [avg_img] + [var_img] + avg_var_eig_list
                montage_img= montage_scale(image_list, scale=True)  # moved to sp_utilities in v1.3
                ####montage_img.set_attr('n_objects', num_class_imgs)
                ####montage_list[class_num]= montage_img
                ####if usempi : sp_utilities.bcast_EMData_to_all(montage_list[class_num], myid, main_mpi_proc, mpi.MPI_COMM_WORLD)
                mpi_ave_stack_obj.insert_clip(montage_img, (0, 0, class_mpi_idx) )
                avg_dict= avg_img.get_attr_dict()
                avg_dict['data_path']= sp_utilities.makerelpath(mpi_ave_bdb_name, mpi_ave_stack_path)
                avg_dict['n_objects']= num_class_imgs
                avg_dict["ptcl_source_coord_id"]= class_mpi_idx
                
                # Add memership list to header
                class_members= sp_utilities.read_text_row(classdoc_current)
                avg_dict['members']= [member[0] for member in class_members]
                # (Do we want to use the rest of the class-average header also?)
                
                # New dictionary entry
                mpi_ave_bdb_dict[class_mpi_idx]= avg_dict
                
            # If applying centering
            else:  # if center_params_type:
                if center_params_type=='combined':
                    interpolation_method= 'quadratic'
                else:
                    interpolation_method= 'linear'
                
                # Sanity checks
                assert shrink_factor == 1, "ERROR!! Shrink factor %s assumed to be 1 with centering option" % shrink_factor
                assert filter_radius == None, "ERROR!! Filter radius %s assumed to be unused with centering option" % filter_radius
                
                num_class_imgs= align_filter_shrink(
                        class_init_bdb_name, 
                        output_class_stack, 
                        aligned_bdb_name, 
                        alignYN=True, 
                        filter_radius=None, 
                        shrink_factor=1, 
                        usectf=ctf_method, 
                        method=interpolation_method, 
                        log=log, 
                        verbosity=verbosity, 
                        is_main=is_main
                        )
                
                if center_params_type=='combined' or center_debugTF:
                    # For more intshifts, update xform.align2d from my arbitrary header
                    if center_params_type=='intshifts':
                        # Copy rotation angle to header position used by averaging functions ('xform.align2d')
                        copy_header_attr(aligned_bdb_name, ROTATION_TAG, 'xform.align2d')
                        
                        align_mode=True
                    elif center_params_type=='combined':
                        assert not center_debugTF, "ERROR!! Flag '--debug' only used with --applyparams='intshifts'"
                        
                        align_mode=False
                    else:
                        print_log_msg('ERROR!! Disallowed --applyparams=%s' % center_params_type, log, is_main)
                        exit()
                    
                    avg_img, var_img= avg_optional_ctf(
                        aligned_bdb_name, 
                        ctf_method=ctf_method, 
                        pws_docfile=pws_doc_template.format(class_num), 
                        do_align=align_mode
                        )
                    
                    image_list= [avg_img] + [var_img]  # + avg_var_eig_list
                    montage_img= montage_scale(image_list, scale=True)  # moved to sp_utilities in v1.3
                    ####montage_img.set_attr('n_objects', num_class_imgs)
                    ####montage_list[class_num]= montage_img
                    ####if usempi : sp_utilities.bcast_EMData_to_all(montage_list[class_num], myid, main_mpi_proc, mpi.MPI_COMM_WORLD)
                    mpi_ave_stack_obj.insert_clip(montage_img, (0, 0, class_mpi_idx) )
                    avg_dict= avg_img.get_attr_dict()
                    avg_dict['data_path']= sp_utilities.makerelpath(mpi_ave_bdb_name, mpi_ave_stack_path)
                    avg_dict['n_objects']= num_class_imgs
                    avg_dict["ptcl_source_coord_id"]= class_mpi_idx
                    mpi_ave_bdb_dict[class_mpi_idx]= avg_dict
                # End combined if-then
            # End centering if-then
        # End ISAC if-then
    # End class loop
    
    if do_align:
        mpi_ave_stack_obj.write_image(mpi_ave_stack_path)
        EMAN2db.db_close_dict(mpi_ave_bdb_name)
    
    ####quick_barrier(myid)
    ####if usempi:#### and not is_main: 
        ####for class_num in range(mpi_start, mpi_end):
            ####printvars(['myid','class_num'])
            ####sp_utilities.bcast_EMData_to_all(montage_list[class_num], myid, main_mpi_proc, mpi.MPI_COMM_WORLD)

    ##### Write average/variance/eigenimages
    ####if is_main:
        ########if usempi: 
            ########for class_num in range(mpi_start, mpi_end):
                ########printvars(['myid','class_num'])
                ########sp_utilities.bcast_EMData_to_all(montage_list[class_num], myid, main_mpi_proc, mpi.MPI_COMM_WORLD)
        
        ####for class_num in range(num_classes):
            ####montage_list[class_num].write_image(montage_file, class_num)
    
    quick_barrier()
    
    if do_align and is_main:
        # Merge BDBs
        sort_merge_bdbs(
            mpi_ave_bdb_template, 
            num_proc, 
            outdir, 
            montage_file, 
            log=log, 
            verbosity=verbosity
            )
    
    return tot_parts
    
def prepare_2d_forPCA(data, mode = "a", output_stack = None, CTF = False):
    """
        Prepare 2D images for PCA
        
        Adapted from SPHIRE 1.2 sp_applications
        
        Average of all images is calculated using header alignment information, 
        subtracted from each image and the difference is written to the output_stack
        If CTF, the average is calculated as
        Av = sum(CTF_k*Im_k)/sum(CTF_k^2)
        and the difference as
        CTF_k(Im_k - CTF_k*Av)/sum(CTF_k^2)
        average outside of a circle r = nx//2-1 is subtracted from each image
    """
    dopa = True
    if type(data) == type(""):
        inmem = False
    else:
        inmem = True

    if inmem:
        n = len(data)
    else:
        n = EMAN2.EMUtil.get_image_count(data)

    if inmem:
        img = data[0]
    else:
        img = sp_utilities.get_im(data,0)

    nx = img.get_xsize()
    ny = img.get_ysize()
    nz = img.get_zsize()
    
    if( output_stack == None):  outstack = [None]*n

    mask = sp_utilities.model_circle( nx//2-2, nx, ny)
    if  CTF:
        if(img.get_attr_default('ctf_applied', 0) > 0):
            ERROR("data cannot be ctf-applied","prepare_2d_forPCA",1)

        nx2 = 2*nx
        ny2 = 2*ny
        ave       = EMData(nx2, ny2, 1, False)
        ctf_2_sum = EMData(nx2, ny2, 1, False)

        for i in range(n):
            if inmem:
                img = data[i].copy()
            else:
                img = sp_utilities.get_im(data, i)
            ctf_params = img.get_attr("ctf")
            if (mode == 'a'):
                angle, sx, sy, mirror, scale = sp_utilities.get_params2D(img)
                img = sp_fundamentals.rot_shift2D(img, angle, sx, sy, mirror, scale)
            st = EMAN2.Util.infomask(img, mask, False)
            img -= st[0]
            img = sp_utilities.pad(img, nx2,ny2, 1, background = "circumference")
            sp_fundamentals.fftip(img)
            EMAN2.Util.add_img(ave, sp_filter.filt_ctf(img, ctf_params))
            EMAN2.Util.add_img2(ctf_2_sum, sp_morphology.ctf_img(nx2, ctf_params))
        EMAN2.Util.div_filter(ave, ctf_2_sum)
        for i in range(n):
            if inmem:
                img = data[i].copy()
            else:
                img = sp_utilities.get_im(data, i)
            ctf_params = img.get_attr("ctf")
            if (mode == 'a'):
                angle, sx, sy, mirror, scale = sp_utilities.get_params2D(img)
                img = sp_fundamentals.rot_shift2D(img, angle, sx, sy, mirror, scale)
            st = EMAN2.Util.infomask(img, mask, False)
            img -= st[0]
            img = sp_utilities.pad(img, nx2,ny2, 1, background = "circumference")
            sp_fundamentals.fftip(img)
            img = sp_filter.filt_ctf(img - sp_filter.filt_ctf(ave, ctf_params, dopa), ctf_params, dopa)
            EMAN2.Util.div_filter(img, ctf_2_sum)
            img = sp_fundamentals.window2d(sp_fundamentals.fft(img),nx,ny)
            sp_utilities.set_params2D(img, [0.0,0.0,0.0,0,1.0])
            if( output_stack == None):  outstack[i] = img
            else:                       img.write_image(output_stack, i)
    else:
        ave  = sp_utilities.model_blank( nx, ny)
        for i in range(n):
            if inmem:
                img = data[i].copy()
            else:
                img = sp_utilities.get_im(data, i)
            angle, sx, sy, mirror, scale = sp_utilities.get_params2D(img)
            img = sp_fundamentals.rot_shift2D(img, angle, sx, sy, mirror, scale)
            st = EMAN2.Util.infomask(img, mask, False)
            img -= st[0]
            EMAN2.Util.add_img(ave, img)
        ####ave /= float(n)    #### (TODO: check whether float() makes a difference)
        ave = old_div(ave, n)
        for i in range(n):
            if inmem:
                img = data[i].copy()
            else:
                img = sp_utilities.get_im(data, i)
            angle, sx, sy, mirror, scale = sp_utilities.get_params2D(img)
            img = sp_fundamentals.rot_shift2D(img, angle, sx, sy, mirror, scale)
            st = EMAN2.Util.infomask(img, mask, False)
            img -= st[0]
            EMAN2.Util.sub_img(img, ave)
            sp_utilities.set_params2D(img, [0.0,0.0,0.0,0,1.0])
            if( output_stack == None):  outstack[i] = img
            else:                       img.write_image(output_stack, i)
    if( output_stack == None):  return ave, outstack
    else:                       return None

def pca(input_stacks, subavg="", mask_radius=-1, nvec=3, incore=False, shuffle=False, genbuf=True, maskfile="", MPI=False, verbose=False):
    """
        PCA of a set of images (can be 1-2-3-D).
        Adapted from SPHIRE 1.2 sp_applications
        
        input_stacks - 
        subavg       - file name containing the average of the input stack.  If None, average will not be subtracted
        mask_radius  - radius of a spherical mask, cannot be specified if maskfile provided
        nvec         - number of egeinimages to be computed
        incore       - do in-core calculations, preferable for small datasets (default False)
        shuffle      - Shuffle test (default False)
        genbuf       - generate disk buffer (default True), to use the disk buffer with data set to False
        maskfile     - name of the mask file 
    """

    if type(input_stacks[0]) is str: data_on_disk = True	 # input_stacks is a file name
    else:
        data_on_disk = False # input_stacks is a list of images not a file name
        if MPI:
            ERROR('MPI version for data in memory version is not implemented', "pca", 1)

    if mask_radius > 0 and maskfile !="":
        ERROR('Error: mask radius and mask file cannot be used at the same time', "pca", 1)

    if mask_radius >0:

        if(verbose): sxprint("Using spherical mask, rad=", mask_radius)

        if maskfile!="":   ERROR('mask radius and mask file cannot be used at the same time', "pca", 1)
        if data_on_disk:
            data = sp_utilities.get_im( input_stacks[0] )
        else:
            data = input_stacks[0]
        mask = sp_utilities.model_circle(mask_radius, data.get_xsize(), data.get_ysize(), data.get_zsize())

    elif(maskfile!="") :
        if(verbose): sxprint("Using mask: ", maskfile)
        mask = sp_utilities.get_image( maskfile )
    else:
        data = EMAN2.EMData()
        if data_on_disk:
            data.read_image( input_stacks[0], 0, True)
        else:
            data = input_stacks[0]
        mask = sp_utilities.model_blank(data.get_xsize(), data.get_ysize(), data.get_zsize(), bckg=1.0)

    pca = pcanalyzer(mask, nvec, incore, MPI)

    if subavg != "":
        if(verbose): sxprint("Subtracting ", subavg, " from each image")
        avg = sp_utilities.get_image( subavg )
        pca.setavg( avg )

    if data_on_disk:
        files = file_set( input_stacks )
    if MPI:
        myid = mpi.mpi_comm_rank( mpi.MPI_COMM_WORLD )
        ncpu = mpi.mpi_comm_size( mpi.MPI_COMM_WORLD )
    else:
        myid = 0
        ncpu = 1

    if genbuf:
        if shuffle: ERROR('Shuffle works only with usebuf', "pca", 1)

        if data_on_disk:
            bgn,end = sp_applications.MPI_start_end( files.nimg(), ncpu, myid )
        else:
            bgn,end = sp_applications.MPI_start_end( len(input_stacks), ncpu, myid )
        for i in range(bgn,end):
            if data_on_disk:
                fname, imgid = files.get( i )
                data = sp_utilities.get_im( fname, imgid)
                if(verbose):  sxprint("Inserting image %s, %4d" % (fname, imgid))
            else:
                data = input_stacks[i]
            pca.insert( data )

    else:
        pca.usebuf( )
        if(verbose):  sxprint(myid, "using existing buff, nimg: ", pca.nimg)
        if shuffle:
            pca.shuffle()

    vecs = pca.analyze()
    return vecs
    
class pcanalyzer(object):
    # Adapted from sp_statistics
    
    def __init__(self, mask, nvec=3, incore=False, MPI=False, scratch=None):
        self.mask = mask.copy()
       
        if MPI:
            self.myid = mpi.mpi_comm_rank( mpi.MPI_COMM_WORLD )
            if( scratch == None):
                self.file = os.path.join("." , "maskedimg%04d.bin" % self.myid )
            else:
                self.file = os.path.join(scratch , "maskedimg%04d.bin" % self.myid )
            self.MPI  = True
        else:
            if( scratch == None):
                self.file = os.path.join("." , "maskedimg.bin" )
            else:
                self.file = os.path.join(scratch , "maskedimg.bin" )
            self.MPI  = False
            self.myid = 0
        
        self.nimg   = 0
        self.nvec   = nvec
        self.fw     = None
        self.fr     = None
        self.avgdat = None
        self.myBuff = []
        self.myBuffPos = 0
        self.incore = incore

    def writedat( self, data ):
        if not self.incore:
            if self.fw is None:
                self.fw = open( self.file, "wb" )
            data.tofile( self.fw )
        else:
            if len(self.myBuff) <= self.myBuffPos:
                self.myBuff.append(data.copy())
                self.myBuffPos = len(self.myBuff)
            else:
                self.myBuff[self.myBuffPos] = data.copy()
                self.myBuffPos += 1

    def read_dat( self, data ):
        if not self.incore:
            if not(self.fw is None) and not( self.fw.closed ):
                self.fw.close()
            if self.fr is None:
                self.fr = open( self.file, "rb" )
            assert not(self.fr is None) and not self.fr.closed
            EMAN2.Util.readarray( self.fr, data, self.ncov )
        else:
            data[:] = self.myBuff[self.myBuffPos]
            self.myBuffPos += 1
        if not(self.avgdat is None):
            data -= self.avgdat

    def close_dat( self ):
        if not self.incore:
            if not(self.fw is None) and not( self.fw.closed ):
                self.fw.close()
            self.fw = None
            if not(self.fr is None) and not( self.fr.closed ):
                self.fr.close()
            self.fr = None	
        else:
            self.myBuffPos = 0

    def usebuf( self ):
        nx = self.mask.get_xsize()
        ny = self.mask.get_ysize()
        nz = self.mask.get_zsize()

        self.ncov = 0
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    if( self.mask.get_value_at(ix,iy,iz) >= 0.5 ):
                        self.ncov += 1
        size = os.stat( self.file )[6]
        self.nimg = old_div( size, (self.ncov*4) )
        assert self.nimg * self.ncov*4 == size
        self.bufused = True

    def shuffle( self ):
        assert self.bufused
        random.seed( 10000 + 10*self.myid )

        shfflfile = string.replace( self.file, "masked", "shuffled" )

        sumdata = np.zeros( (self.ncov), np.float32 )
        imgdata = np.zeros( (self.ncov), np.float32 )
        if not self.incore: 
            self.fr = open( self.file, "rb" )
        self.avgdata = None

        if not self.incore: 
            fw = open( shfflfile, "wb" )
        for i in range(self.nimg):
            self.read_dat( imgdata )
            random.shuffle( imgdata )
            sumdata += imgdata
            if not self.incore:
                imgdata.tofile( fw )
            else:
                self.myBuff[self.myBuffPos-1] = imgdata.copy()

        if not self.incore: 
            self.fr.close()
            fw.close()
        else:
            self.close_dat()
        
        if self.MPI:
            sumdata = mpi.mpi_reduce( sumdata, self.ncov, mpi.MPI_FLOAT, mpi.MPI_SUM, 0, mpi.MPI_COMM_WORLD )
            sumdata = mpi.mpi_bcast(  sumdata, self.ncov, mpi.MPI_FLOAT, 0, mpi.MPI_COMM_WORLD )
            sumdata = np.array(sumdata, np.float32)

            sumnimg = mpi.mpi_reduce( self.nimg, 1, mpi.MPI_INT, mpi.MPI_SUM, 0, mpi.MPI_COMM_WORLD )
            sumnimg = mpi.mpi_bcast(  sumnimg,   1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD )
        else:
            sumnimg = self.nimg

        self.file = shfflfile
        self.avgdat = old_div( sumdata[:], float(sumnimg) )

    def setavg( self, avg ):
        tmpimg = EMAN2.Util.compress_image_mask( avg, self.mask )
        avgdat = sp_utilities.get_image_data(tmpimg)
        self.avgdat = np.zeros( (len(avgdat)), np.float32 )
        self.avgdat[:] = avgdat[:]

    def insert( self, img ):
        assert self.mask.get_xsize()==img.get_xsize()
        assert self.mask.get_ysize()==img.get_ysize()
        assert self.mask.get_zsize()==img.get_zsize()

        tmpimg = EMAN2.Util.compress_image_mask( img, self.mask )
        tmpdat = sp_utilities.get_image_data(tmpimg)
        if self.incore:
            self.myBuffPos = len(self.myBuff)
        self.writedat( tmpdat )
        if self.incore:
            self.close_dat()                                   #   WRITEDAT
        self.nimg +=1
        self.ncov = tmpimg.get_xsize()

    def analyze( self ):
        ncov = self.ncov
        kstep = self.nvec + 20 # the choice of kstep is purely heuristic

        diag    = np.zeros( (kstep), np.float32 )
        subdiag = np.zeros( (kstep), np.float32 )
        vmat    = np.zeros( (kstep, ncov), np.float32 )

        lanczos_start = time.time()
        kstep = self.lanczos( kstep, diag, subdiag, vmat )
        if not self.MPI or self.myid==0:
            qmat = np.zeros( (kstep,kstep), np.float32 )
            lfwrk = 100 + 4*kstep + kstep*kstep
            liwrk =   3 + 5*kstep

            fwork = np.zeros( (lfwrk), np.float32 )
            iwork = np.zeros( (liwrk), np.int32 )
            info = EMAN2.Util.sstevd( "V", kstep, diag, subdiag, qmat, kstep, fwork, lfwrk, iwork, liwrk)

            eigval = np.zeros( (self.nvec), np.float32 )
            for j in range(self.nvec):
                eigval[j] = diag[kstep-j-1]

            eigimgs = []
            for j in range(self.nvec):
                tmpimg = sp_utilities.model_blank(ncov, 1, 1)
                eigvec = sp_utilities.get_image_data( tmpimg )
                trans = 'N'
                EMAN2.Util.sgemv( trans, ncov, kstep, 1.0, vmat, ncov, qmat[kstep-j-1], 1, 0.0, eigvec, 1 );

                eigimg = EMAN2.Util.reconstitute_image_mask(tmpimg, self.mask)
                eigimg.set_attr( "eigval", old_div( float(eigval[j]), (self.nimg - 1) ) )
                eigimgs.append( eigimg )

            return eigimgs


    def lanczos( self, kstep, diag, subdiag, V ):
        all_start = time.time()

        ncov = self.ncov
        v0 = np.zeros( (ncov), np.float32)
        Av = np.zeros( (ncov), np.float32)

        hvec = np.zeros( (kstep), np.float32 )
        htmp = np.zeros( (kstep), np.float32 )
        imgdata = np.zeros( (ncov), np.float32 )

        for i in range(ncov):
            v0[i] = 1.0

        beta = EMAN2.Util.snrm2(ncov, v0, 1)
        for i in range(ncov):
            V[0][i] = old_div(v0[i], beta)

        for i in range(self.nimg):
            self.read_dat(imgdata)                                     #  READ_DAT			
            alpha = EMAN2.Util.sdot( ncov, imgdata, 1, V[0], 1 )
            EMAN2.Util.saxpy( ncov, alpha, imgdata, 1, Av, 1 )
        self.close_dat()

        if self.MPI:
            Av = mpi.mpi_reduce( Av, ncov, mpi.MPI_FLOAT, mpi.MPI_SUM, 0, mpi.MPI_COMM_WORLD )
            Av = mpi.mpi_bcast(  Av, ncov, mpi.MPI_FLOAT, 0, mpi.MPI_COMM_WORLD )
            Av = np.array(Av, np.float32)

        diag[0] = EMAN2.Util.sdot( ncov, V[0], 1, Av, 1 )
        alpha = -diag[0]
        EMAN2.Util.saxpy( ncov, float(alpha), V[0], 1, Av, 1 )

        TOL = 1.0e-7
        for iter in range(1, kstep):
            iter_start = time.time()
            beta = EMAN2.Util.snrm2(ncov, Av, 1)
            if( beta < TOL ):
                kstep = iter+1
                break

            subdiag[iter-1] = beta
            for i in range(ncov):
                V[iter][i] = old_div(Av[i], beta)

            Av[:] = 0.0

            for i in range(self.nimg):
                self.read_dat( imgdata )                                #READ_DAT
                alpha = EMAN2.Util.sdot( ncov, imgdata, 1, V[iter], 1 )
                EMAN2.Util.saxpy( ncov, float(alpha), imgdata, 1, Av, 1 )
            self.close_dat()


            if self.MPI:
                Av = mpi.mpi_reduce( Av, ncov, mpi.MPI_FLOAT, mpi.MPI_SUM, 0, mpi.MPI_COMM_WORLD )
                Av = mpi.mpi_bcast(  Av, ncov, mpi.MPI_FLOAT, 0, mpi.MPI_COMM_WORLD )
                Av = np.array(Av, np.float32)

            trans = 'T'
            EMAN2.Util.sgemv( trans, ncov, iter+1,  1.0, V, ncov, Av, 1,
                        0.0, hvec, 1 )

            trans = 'N'
            EMAN2.Util.sgemv( trans, ncov, iter+1, -1.0, V, ncov, hvec, 1,
                        1.0,     Av, 1 )

            trans = 'T'
            EMAN2.Util.sgemv( trans, ncov, iter+1,  1.0, V, ncov, Av, 1,
                        0.0,   htmp, 1 )

            EMAN2.Util.saxpy(iter+1, 1.0, htmp, 1, hvec, 1)

            trans = 'N'
            EMAN2.Util.sgemv( trans, ncov, iter+1, -1.0, V, ncov, htmp, 1,
                        1.0,     Av, 1 )

            diag[iter] = hvec[iter]
        
        return kstep

def count_total_particles(num_classes, bdb_template):
    """
    Gets total number of particles from a series of BDB stacks.
    
    Arguments:
        num_classes : number of classes
        bdb_template : template for BDB stacks, without leading 'bdb:' or trailing '.bdb'
    
    Returns:
        number of total particles
    """
    
    tot_parts= 0
    
    for class_num in range(num_classes):
        bdb_name= 'bdb:' + bdb_template.format(class_num)
        class_parts= EMAN2.EMUtil.get_image_count(bdb_name)
        tot_parts+= class_parts
        
    return tot_parts
        
def copy_header_attr(bdb_name, source_tag, destination_tag):
    """
    Copies header values from one location to another.
    Hardwired for rotation -> xform.align2d, but can be generalized.
    
    Arguments:
        bdb_name : BDB to be modified
        source_tag : source header tag
        destination_tag : destination header tag
    """
    
    rotation_list= EMAN2.EMUtil.get_all_attributes(bdb_name, source_tag)
    num_imgs= EMAN2.EMUtil.get_image_count(bdb_name)
    local_bdb_stack= EMAN2db.db_open_dict(bdb_name)
    
    # Loop through particles
    for idx in range(num_imgs):
        img= sp_utilities.get_im(bdb_name, idx)
        part_header= img.get_attr_dict()
        part_header[destination_tag]= EMAN2.Transform({"type": "2D", "alpha": rotation_list[idx]})
        local_bdb_stack[idx]= part_header
    # End particle loop
    
    # Close BDB
    EMAN2db.db_close_dict(local_bdb_stack)

def avg_optional_ctf(bdb_or_list, ctf_method=None, pws_docfile=None, do_align=True):
    """
    
    Arguments:
        bdb_or_list : input image BDB stack or list of EMData objects
        ctf_method : method for CTF correctioni: None, 'flip', 'wiener'
        pws_docfile : if phase-flipping, 1D profile of amplitudes before & after CTF-correction
        do_align : (boolean) whether to apply alignment to images in input stack
    
    Returns:
        average image
        variance image
    """
    
    # Alignment is applied by default
    if do_align: 
        align_mode= 'a'
    else:
        align_mode= 'n'
    
    if ctf_method == 'wiener':
        avg_img, var_img= sp_statistics.avgvar_ctf(bdb_or_list, dopa=True, mode=align_mode)
    elif ctf_method == 'flip':
        avg_img, var_img, pxsz= avgvar_flipctf(bdb_or_list, dopa=False, mode=align_mode)
        avg_noctf, var_noctf=   avgvar(bdb_or_list, mode=align_mode)  # move to sp_statistics
        
        # Compute 1D power spectra
        rops_ctf= sp_fundamentals.rops_table(avg_img)
        rops_noctf= sp_fundamentals.rops_table(avg_noctf)
        spfreq= [ old_div( x, 2.0 * pxsz * (len(rops_ctf)-1) ) for x in range(len(rops_ctf))]
        pws_tuple_list= list( zip(spfreq, rops_ctf, rops_noctf) )  # list() should work in both python 2 & 3
        ####pws_list= map(list, pws_tuple_list)  # wrote_text_row won't like list of tuples
        pws_list = [list(x) for x in pws_tuple_list]  # wrote_text_row won't like list of tuples
        
        # Write 1D power spectra
        if pws_docfile : sp_utilities.write_text_row(pws_list, pws_docfile)
    elif ctf_method == None:
        avg_img, var_img= avgvar(bdb_or_list, mode=align_mode)  # move to sp_statistics
    else:
        print("ERROR!! ctf_method %s not known" % ctf_method)
        exit()
    
    return avg_img, var_img

def avgvar_flipctf(data, mode='a', interp='quadratic', i1=0, i2=0, use_odd=True, use_even=True, snr=1.0, dopa = True):
    '''
    Averages data while applying phase-flip.
    Adapted from sp_statistics.avgvar_ctf
    
    Parameters
        data: image stack, must be 2D, must be in real space
        mode: whether to apply alignment parameters. Default mode='a' means apply parameters
        rot_method: specifies the function by which images are rotated/shifted if alignment parameters are to be applied. 
            This is only relevant for the case where images are 2D, 
            in which case rot_method can be either rot_shift2D or combined2dg, with the default being rot_shift2D. 
            If images are 3D, rot_shift3D will be used to rotate/shift the images.
        interp: interpolation method to use for rot_method when applying alignment parameters.
        i1: index of first image to be used.
        i2: index of last image to be used. If i2 = 0, then i2 defaults to one less than number of images in the data
        use_odd: images with indices between i1 and i2 which are odd are used if and only if use_odd is set to True. Default is True.
        use_even: images with indices between i1 and i2 which are even are used if and only if use_even is set to True. Default is True.
        snr: signal to noise ratio, default 1.0
    
    Returns:
        average
        variance
        pixel size
    '''

    inmem = True
    if type(data) == type(""):
        inmem = False

    if inmem:
        img = data[0]
    else:
        img = sp_utilities.get_im(data,0)
    nx = img.get_xsize()
    ny = img.get_ysize()
    nz = img.get_zsize()
    if nz > 1:
        sp_global_def.ERROR("images must be 2D for CTF correction.....exiting","avgvar_ctf",1)

    if img.get_attr_default('ctf_applied', 0) == 1:
        sp_global_def.ERROR("data cannot be ctf-applied....exiting","avgvar_ctf",1)

    if inmem:
        data_nima = len(data)
    else:
        data_nima = EMAN2.EMUtil.get_image_count(data)

    if i2 == 0: i2 = data_nima-1
    if dopa:
        nx2 = nx*2
        ny2 = ny*2
    else:
        nx2 = nx
        ny2 = ny
    ave = EMAN2.EMData(nx2, ny2, 1, False)
    var = EMAN2.EMData(nx2, ny2, 1, True)

    nima = 0
    for i in range(i1, i2+1):
        if not(use_odd) and i%2 == 1: continue
        if not(use_even) and i%2 == 0: continue
        nima += 1
        if inmem: img = data[i].copy()
        else: img = sp_utilities.get_im(data, i)

        ctf_params = img.get_attr("ctf")

        if(mode == 'a'):
            angle, sx, sy, mirror, scale = sp_utilities.get_params2D(img)
            img = sp_fundamentals.rot_shift2D(img, angle, sx, sy, mirror, scale, interp)
            ctf_params.dfang += angle
            if mirror == 1:  ctf_params.dfang = 270.0-ctf_params.dfang

        img = sp_utilities.pad(img, nx2, ny2, 1, background="circumference")
        sp_fundamentals.fftip(img)
        img= sp_filter.filt_ctf(img, ctf_params, binary=1)
        
        # Update average and variance
        EMAN2.Util.add_img(ave, img)
        EMAN2.Util.add_img2( var, sp_fundamentals.fft(img) )

    EMAN2.Util.mul_scalar(ave, old_div(1.0, float(nima) ) )  # The "1.0" and "float()" are probably unnecessary
    
    # Inverse FT
    avg_padded= sp_fundamentals.fft(ave)
    var_padded= (var - avg_padded*avg_padded*nima)/(nima-1)
    
    # Back to original size
    avg_unpad= sp_fundamentals.window2d(avg_padded, nx, ny)
    var_unpad= sp_fundamentals.window2d(var_padded, nx, ny)
    
    return avg_unpad, var_unpad, ctf_params.apix
    
def avgvar(data, mode='a', interp='quadratic', i1=0, i2=0, use_odd=True, use_even=True):
    '''
    
    INPUT
    
    data: image stack, can be 2D or 3D, must be in real space
    mode: whether to apply alignment parameters. Default mode='a' means apply parameters
    rot_method: specifies the function by which images are rotated/shifted if alignment parameters are to be applied. This is only relevant for the case where images are 2D, in which case rot_method can be either rot_shift2D or rotshift2dg, with the default being rot_shift2D. If images are 3D, rot_shift3D will be used to rotate/shift the images.
    interp: interpolation method to use for rot_method when applying alignment parameters.
    i1: index of first image to be used.
    i2: index of last image to be used. If i2 = 0, then i2 defaults to one less than number of images in the data
    use_odd: images with indices between i1 and i2 which are odd are used if and only if use_odd is set to True. Default is True.
    use_even: images with indices between i1 and i2 which are even are used if and only if use_even is set to True. Default is True.
    
    OUTPUT
        
    ave: the average of the image series in real space
    var: the variance of the image series in real space

    '''
    
    #### TODO: modernize and return to sp_statistics
    
    inmem = True
    if type(data) == type(""):
        inmem = False
        from sp_utilities    import get_im	

    img2D = True
    if inmem:
        img = data[0]
    else:
        img = get_im(data,0)
    nx = img.get_xsize()
    ny = img.get_ysize()
    nz = img.get_zsize()
    if nz > 1:
        img2D = False

    if mode == 'a':
        if img2D:
            from sp_utilities import get_params2D
            from sp_fundamentals import rot_shift2D
        else:
            from sp_utilities import get_params3D
            from sp_fundamentals import rot_shift3D

    if inmem:
        data_nima = len(data)
    else:
        data_nima = EMAN2.EMUtil.get_image_count(data)
    if i2 == 0: i2 = data_nima-1

    ave = sp_utilities.model_blank(nx,ny,nz)
    var = sp_utilities.model_blank(nx,ny,nz)
    nima = 0
    for i in range(i1, i2+1):
        if not(use_odd) and i%2 == 1:
            continue
        if not(use_even) and i%2 == 0:
            continue
        nima += 1
        if inmem:
            img = data[i]
        else:
            img = get_im(data, i)
        if (mode == 'a'):
            if img2D:
                angle, sx, sy, mirror, scale = get_params2D(img)
                img = rot_shift2D(img, angle, sx, sy, mirror, scale, interp)
            else:
                phi, theta, psi, s3x, s3y, s3z, mirror, scale = get_params3D(img)
                img = rot_shift3D(img, phi, theta, psi, s3x, s3y, s3z, scale)
        EMAN2.Util.add_img(ave, img)
        EMAN2.Util.add_img2(var, img)

    EMAN2.Util.mul_scalar(ave, 1.0 /float(nima) )
    return ave, (var - ave*ave*nima)/(nima-1)

def montage_scale(inputstack, ncol=None, marginwidth=0, bkgd=0, outfile=None, scale=False):
    """
    Generates montage of images into one image.
    Adapted from sxmontage.py and SPIDER's MN S
    Modified from in sp_utilities in SPHIRE 1.3
    2021-01-04 -- Updated for Python 3
    
    Arguments:
        inputstack : Stack of input images to merge into montage
        ncol : Number of images per row (default: all on one row)
        marginwidth : Margin width, pixels
        bkgd : Background value of montage
        outfile : Optional output file with montage output
        scale : Normalize images from 0..2
    Returns:
        montage : EMData object of image montage
    """
    
    if isinstance(inputstack, str): inputstack= EMAN2.EMData.read_images(inputstack)
    
    # Get single-image dimensions
    nx= inputstack[0].get_xsize()
    ny= inputstack[0].get_ysize()
    
    # Get number of images and calculate montage dimensions
    numimgs= len(inputstack)
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
    for imgnum in range(numimgs):
        # Horizontal grid position is image# modulo NCOL
        colnum= imgnum % ncol
        
        # Montage is numbered from the top down
        rownum= numrows - 1 - old_div(imgnum, ncol)
        
        xoffset= colnum*(nx+marginwidth)
        yoffset= rownum*(ny+marginwidth)
        
        if scale == True:
            currentx= inputstack[imgnum].get_xsize()
            assert currentx == nx, "ERROR!! Mismatched dimensions %s != %s" % (currentx, nx)
            [avg,var,fmin,fmax]= EMAN2.Util.infomask(inputstack[imgnum], mask, True)
            
            try:
                img_norm= (inputstack[imgnum] - fmin) * old_div(2., fmax - fmin)
            except ZeroDivisionError:
                print(f"WARNING! Blank image in #{imgnum} out of {numimgs}")
                img_norm= inputstack[imgnum]
        else:
            img_norm= inputstack[imgnum]
        insert_image(img_norm, montage, xoffset, yoffset)
    
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
            try:
                largeimage.set_value_at(xoffset+xcoord, int(yoffset+ycoord), getpixel)
            except:
                print("ERROR!! xoffset+xcoord %s, yoffset+ycoord %s" % (
                    xoffset+xcoord, yoffset+ycoord
                    ))
                printvars(['xcoord','xoffset','ycoord','yoffset','largeimage','getpixel'], True, True)

def align_filter_shrink(input_bdb_name, output_stack_path, output_bdb_name, alignYN=False, \
                filter_radius=None, shrink_factor=1, usectf=None, method='quadratic', log=None, verbosity=0, is_main=True):
    """
    Filters and shrinks image stack.
    
    Arguments:
        input_bdb_name : input BDB stack (in the form bdb:DIRECTORY#STACK)
        output_stack_path : name for output image stack (MRCS/HDF/etc.)
        output_bdb_name : name for output BDB stack (in the form bdb:DIRECTORY#STACK)
        alignYN : (boolean) whether to apply alignments
        filter_radius : filter radius, reciprocal pixels
        shrink_factor : downsampling factor
        usectf: CTF options
        method : interpolation method
        log : instance of Logger class
        verbosity : how much information to write to screen (0..6)
        is_main : (boolean) if main MPI process
    """
    
    input_stack= EMAN2.EMData.read_images(input_bdb_name)
    num_imgs= EMAN2.EMUtil.get_image_count(input_bdb_name)
    
    box_size, sub_rate= get_rescaled_dimension(input_stack, shrink_factor, log=log, verbose=verbosity>=6)
    
    # Initialize stack & BDB
    aligned_stack_obj= EMAN2.EMData(box_size, box_size, num_imgs) 
    assert output_bdb_name[:4] == 'bdb:', "align_filter_shrink, %s does not start with 'bdb:'"
    new_bdb_dict= EMAN2db.db_open_dict(output_bdb_name)
    
    # If we pretend the BDB file is in a subdirectory by substituting a '/' for the "#", then makerelpath is easy
    output_bdb_stem= output_bdb_name[4:].replace("#",'/')
    
    # Turn off progress bar unless verbosity 5
    disable_part= verbosity<5 or not is_main
    
    # Loop through images
    for img_num in tqdm.tqdm(range( len(input_stack) ), unit='part', disable=disable_part, file=sys.stdout):
        img_orig= input_stack[img_num]
        
        alpha, sx, sy, mirror, scale= sp_utilities.get_params2D(img_orig)
        
        try:
            alpha, sx, sy, mirror, scale= sp_utilities.get_params2D(img_orig)
        
        # Probably xform.align2d is missing
        except RuntimeError:  
            print('\nERROR! Exiting with RuntimeError')
            img_prev= input_stack[img_num-1]
            print('\nPrevious particle: %s %s' % (img_num-1, img_prev.get_attr_dict() ) )
            print('\nCurrent particle: %s %s' % (img_num, img_orig.get_attr_dict() ) )
            exit()
            
        # Optionally apply alignment parameters
        if alignYN:
            img_ali= sp_fundamentals.rot_shift2D(img_orig, alpha, sx, sy, mirror, scale, method)
        else:
            img_ali= img_orig
        
        if usectf:
            ctf_params= img_orig.get_attr("ctf")
            img_ali= sp_filter.filt_ctf(img_ali, ctf_params, dopad=False)
        
        if filter_radius:
            img_ali= sp_filter.filt_gaussl(img_ali, filter_radius)
        
        if shrink_factor != 1:
            img_ali= sp_fundamentals.resample(img_ali, sub_rate)
            # (Maybe update resample_ratio)
        
        img_ali_dict= img_ali.get_attr_dict()
        img_ali_dict["data_path"]= sp_utilities.makerelpath(output_bdb_stem, output_stack_path)
        img_ali_dict["ptcl_source_coord_id"]= img_num
        new_bdb_dict[img_num]= img_ali_dict
        
        aligned_stack_obj.insert_clip(img_ali, (0, 0, img_num) )
    # End image-loop
        
    aligned_stack_obj.write_image(output_stack_path)
    EMAN2db.db_close_dict(output_bdb_name)
    if verbosity>=3 : 
        mesg= "  Wrote %s images to '%s' from '%s'" % (num_imgs, output_bdb_name, input_bdb_name)
        print_log_msg(mesg, log, is_main)
    
    return num_imgs

def get_rescaled_dimension(image_or_stack, shrink_factor=1, log=None, verbose=False):
    """
    Gets image dimension and subsampling rate.
    
    Arguments:
        image_or_stack : input stack or image
        shrink_factor : downsampling factor, assumed to be integer
        log : instance of Logger class
        verbose : (boolean) whether to write to screen
    """
    
    if type(image_or_stack) is str:  # used to be 'bytes' in python2
        img_obj= sp_utilities.get_im(image_or_stack)
    elif isinstance(image_or_stack, list): 
        img_obj= image_or_stack[0]
    else:
        img_obj= image_or_stack
    
    try:
        box_size= img_obj.get_attr('nx')
    except AttributeError:
        print("\nAttributeError: img_obj", img_obj, type(img_obj) )
        print("AttributeError: image_or_stack", image_or_stack, type(image_or_stack) )
        exit()
    
    if shrink_factor != 1:
        sub_rate= old_div(1.0, shrink_factor)
        box_size= int( old_div(float(box_size), shrink_factor) + 0.5 )
    else:
        sub_rate= 1.0
        
    if verbose: 
        mesg= '  sub_rate: %s %s, box_size: %s, shrink_factor: %s' % (sub_rate, type(sub_rate), box_size, shrink_factor)
        print_log_msg(mesg, log)
    
    return box_size, sub_rate

def sort_merge_bdbs(input_bdb_stem_template, num_bdbs, output_bdb_dir, output_bdb_stem, 
        doc_template=None, log=None, verbosity=0):
    """
    Merges a series of BDBs. 
    Sorts according to GLOBAL_NUM_TAG, if it is present.
    Doesn't look for gaps in the final BDB.
    
    Modified 2020-12-29
    
    Arguments:
        input_bdb_stem_template : Template for BDB name, w/o leading "bdb:" or trailing ".bdb"
        num_bdbs  : number of BDBs to loop through
        output_bdb_dir : target directory for merged BDB (not including EMAN2DB subdirectory)
        output_bdb_stem : stem of output BDB (that is, without preceding "bdb:" or trailing ".bdb")
        doc_template : optional template of particle lists for each BDB, to make sure of correct number
        verbosity : verbosity level
    """

    if verbosity>=2 : print_log_msg("Merging %s BDBs" % num_bdbs, log)
    
    # Output BDB
    output_bdb_obj= BdbNames(output_bdb_stem, output_bdb_dir)
    if os.path.exists(output_bdb_obj.bdb_path) : os.remove(output_bdb_obj.bdb_path)  # will otherwise merge with pre-existing file

    # We'll need this variable when updating the data_path in the header
    out_eman2db_dir= os.path.abspath(output_bdb_obj.eman2db_dir)
    
    # Open new database
    new_bdb_dict= EMAN2db.db_open_dict(output_bdb_obj.bdb_name)

    img_counter= 0
    disableTF= verbosity>=3 or verbosity<2 or num_bdbs==1
    
    # Loop through BDBs
    for bdb_num in tqdm.tqdm(range(num_bdbs), unit='bdb', disable=disableTF, file=sys.stdout):
        input_bdb_name= "bdb:" + input_bdb_stem_template.format(bdb_num)
        input_bdb_dict= EMAN2db.db_open_dict(input_bdb_name, ro=True)  
        
        # Get number of particles
        bdb_numparts= len(input_bdb_dict)
        if verbosity>=3 : print_log_msg("  Opened %s images from '%s'" % (bdb_numparts, input_bdb_name), log)
        
        # Sanity check
        if doc_template:
            docfile= doc_template.format(bdb_num)
            num_doc_parts= len(sp_utilities.read_text_row(docfile) )
            err= "%s %s particles, %s %s particles" % (docfile, num_doc_parts, input_bdb_name, bdb_numparts)
            assert num_doc_parts == bdb_numparts, err
            
        # Loop through particles
        for partnum in range(bdb_numparts):
            # Read image header
            img_header= input_bdb_dict.get_header(partnum)
            
            # Update relative path between new BDB and image stack
            img_header['data_path']= update_data_path(
                    img_header['data_path'], 
                    input_bdb_stem_template, 
                    out_eman2db_dir)
            
            # (When are the following used?)
            img_header['data_source']= input_bdb_name
            img_header['data_n']= partnum
            
            if GLOBAL_NUM_TAG in img_header:
                # Extract global particle number (what we will sort)
                glonum= img_header[GLOBAL_NUM_TAG]
                
                # Write in new database
                new_bdb_dict[glonum]= img_header
            else:
                new_bdb_dict[img_counter]= img_header
            # End global if-then
            
            img_counter+= 1
        # End particle loop
        
        # Close database
        EMAN2db.db_close_dict(input_bdb_name)
    # End BDB loop
    
    # Close new database
    EMAN2db.db_close_dict(output_bdb_obj.bdb_name)
    numoutimgs= EMAN2.EMUtil.get_image_count(output_bdb_obj.bdb_name)
    assert numoutimgs == img_counter, "Uh oh!! %s: %s != %s" % (output_bdb_obj.bdb_name, numoutimgs, img_counter)
    if verbosity>=2 : print_log_msg("Wrote %s images to '%s'\n" % (numoutimgs, output_bdb_obj.bdb_name), log)
    
def update_data_path(curr_path, in_eman2db_dir, out_eman2db_dir):
    """
    Update relative data path between BDB and particle stack. 
    
    Arguments:
        curr_path : current data_path, from image header
        in_eman2db_dir : EMAN2DB directory (or something at the same depth in the directory tree)
        out_eman2db_dir : output EMAN2DB directory
    Returns:
        rel_path : relative path between source and target
    """
    
    # Get full directory
    bdb2curr_path= os.path.join(in_eman2db_dir, curr_path)
    
    # Get absolute directory
    abs_img_path= os.path.abspath(bdb2curr_path)
    # (abspath removes '..'s from the path.)
    
    # Sanity check
    if not os.path.exists(abs_img_path):
        print("ERROR!! update_data_path: '%s' does not exist!" % abs_img_path)
        printvars(['in_eman2db_dir', 'curr_path', 'out_eman2db_dir', 'abs_img_path'], quitTF=True)
    
    # Get relative path between source (BDB) and target (image stack)
    rel_path= sp_utilities.makerelpath(out_eman2db_dir, abs_img_path)
    
    return rel_path

def bandpass_filter(input_stack, options, outdir, log=None, verbosity=0):
    """
    Applies band-pass filter.
    
    Arguments:
        input_stack : input class-average stack
        options : (namespace) command-line options, run 'sp_eval_isac.py -h' for an exhaustive list
        outdir : output directory
        log : instance of Logger class
    """
    
    # Check whether ISAC or beautifier directory (to confirm pixel size)
    
    # Check pixel size
    if options.apix: 
        pxsz= options.apix
        filter_center= old_div(pxsz, options.bandpass_radius)
        mesg= "Band-pass filter centered at %.3f Angstroms (=%s/%.3f px^-1) and width %.3f px^-1" % (
            options.bandpass_radius, pxsz, filter_center, options.bandpass_width
            )
        if verbosity>=1 : print_log_msg(mesg, log)
    else:
        mesg= "ERROR!! Provide pixel size ('--apix=pixel_size') with this option"
        print_log_msg(mesg, log)
        exit()
    
    output_stack= os.path.join(outdir, BANDPASS_AVGS)
    if os.path.exists(output_stack) : os.remove(output_stack)
    
    # Get low-pass filter radius
    lowpass_radius_absfreq= filter_center + options.bandpass_width
    lowpass_radius_angstroms= old_div(pxsz, lowpass_radius_absfreq)
    if verbosity>=1 : 
        mesg= "Gaussian low-pass filter at %.3f Angstroms (=%s/%.3f px^-1)" % (
            lowpass_radius_angstroms, pxsz, lowpass_radius_absfreq
            )
        print_log_msg(mesg, log)
        cmd= "Command: \ne2proc2d.py %s %s --process=filter.lowpass.gauss:cutoff_abs=%.3f" % (
            input_stack, output_stack, lowpass_radius_absfreq
            )
    
    # Get high-pass filter radius
    highpass_radius_absfreq= filter_center - options.bandpass_width
    if highpass_radius_absfreq <= 0:
        print_log_msg("WARNING! High-pass filter beyond origin. Low-pass filtering only...\n", log)
    else:
        highpass_radius_angstroms= old_div(pxsz, highpass_radius_absfreq)
        if verbosity>=1 : 
            mesg= "Gaussian high-pass filter at %.3f Angstroms (=%s/%.3f px^-1)" % (
                highpass_radius_angstroms, pxsz, highpass_radius_absfreq
                )
            print_log_msg(mesg, log)
            cmd+= " --process=filter.highpass.gauss:cutoff_abs=%.3f" % (highpass_radius_absfreq)
    
    # Set up filter profile
    filter_profile_doc= os.path.join(outdir, BANDPASS_FILTER_DOC)
    filter_profile_plot= os.path.join(outdir, BANDPASS_FILTER_PLOT)
    if verbosity>=1 : 
        mesg= "Will write filter profile to '%s' and '%s'" % (
            filter_profile_doc, filter_profile_plot
            )
        print_log_msg(mesg, log)
    
    if verbosity>=2 : cmd += "\n"
    if verbosity>=1 : print_log_msg(cmd, log)
    
    num_classes= EMAN2.EMUtil.get_image_count(input_stack)
    
    disableTF= verbosity>2 or verbosity<2
    
    # Loop through images
    for class_num in tqdm.tqdm(range( num_classes ), unit='class', disable=disableTF, file=sys.stdout):
        im0= sp_utilities.get_im(input_stack, class_num)
        
        if highpass_radius_absfreq > 0:
            # Apply low-pass and high-pass filters in series
            flt= sp_filter.filt_gaussh(sp_filter.filt_gaussl(im0, lowpass_radius_absfreq), highpass_radius_absfreq)
        else:
            flt= sp_filter.filt_gaussl(im0, lowpass_radius_absfreq)
            flt-= flt['mean']
            
        # Write example !D filter profile + plot to disk
        if class_num == 0:
            rooim0= np.array(sp_fundamentals.rops_table(im0), dtype=np.float)
            rooflt= np.array(sp_fundamentals.rops_table(flt), dtype=np.float)
            filter_profile_list= old_div(rooflt, rooim0)
            sp_utilities.write_text_row(filter_profile_list, filter_profile_doc)
            
            try:
                idim= im0['ny']
                spfreq_list= old_div( np.array(range( len(filter_profile_list) ), dtype=np.float), idim*pxsz )  # does this work with arrays
                plt.plot(spfreq_list, filter_profile_list)
                plt.xlabel(r'Spatial frequency, $\mathrm{\AA}$$^{-1}$')
                plt.title('Filter profile')
                plt.savefig(filter_profile_plot)
                
                plotFailed= False
            except RuntimeError:
                plotFailed= True
        
        if verbosity>=3: 
            print_log_msg("%s: unfiltered['mean'] %s, filtered['mean'] %s" % (
                class_num, im0['mean'], flt['mean']
                ), log)
            
        flt.write_image(output_stack, class_num)
    # End image-loop
    
    num_filt= EMAN2.EMUtil.get_image_count(output_stack)
    assert num_classes == num_filt, "ERROR!!, Number of starting %s and final images %s don't agree!" % (
        num_classes, num_filt
        )
            
    if verbosity>=1 :
        print()
        
        assert os.path.exists(filter_profile_doc), "ERROR!! '%s' not written!" % filter_profile_doc
        print_log_msg("Wrote profile as text to '%s'" % filter_profile_doc, log)
    
        if plotFailed:
            print_log_msg("WARNING! Can't create profile plot (a machine-specific problem). Skipping...", log)
        else:
            assert os.path.exists(filter_profile_plot), "ERROR!! '%s' not written!" % filter_profile_plot
            print_log_msg("Wrote profile plot to '%s'" % filter_profile_plot, log)
        
        assert os.path.exists(output_stack), "ERROR!! '%s' not written!" % output_stack
        print_log_msg("Wrote %s images to '%s'\n" % (num_filt, output_stack), log)
    
    print()
    
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
        help='Input class averages')
    
    parser.add_argument(
        'outdir', type=str, 
        help='Output directory')
    
    parser.add_argument(
        '--particles', 
        help='Input image stack')
    
    parser.add_argument(
        '--align_isac_dir', 
        type=str, 
        default=None, 
        help='ISAC directory, for aligning images')
    
    parser.add_argument(
        '--filtrad', 
        type=float, 
        default=None, 
        help='For optional filtered images, low-pass filter radius (1/px or, if pixel size specified, Angstroms)')
    
    parser.add_argument(
        '--apix', 
        type=float, 
        default=None, 
        help='Pixel size, Angstroms (might be downsampled by ISAC)')
    
    parser.add_argument(
        '--shrink', 
        type=int, 
        default=1, 
        help='Optional downsampling factor')
    
    parser.add_argument(
        '--ctf', 
        type=str, 
        default=None, 
        help='CTF-correction options for averages: flip, wiener')
    
    parser.add_argument(
        '--chains_radius', 
        type=int, 
        help='Alignment radius for internal ordering of averages (Hint: Use 29 for ISAC averages)')
    
    parser.add_argument(
        '--chains_exe', 
        type=str, 
        default=None, 
        help='Path for sp_chains.py executable if not default')
    
    parser.add_argument(
        '--verbosity', "-v", 
        type=int, 
        default=2,
        help='Increase verbosity (0..%s)' % MAX_VERBOSITY)
    """
    Verbosity levels:
        0: Only warnings
        1: Basic
        2: Progress bar over classes
        3: Merged classes
        3: Summary in align_filter_shrink()
        3: Stats in bandpass_filter()
        4: e2bdb command line for each class in sane_bdb_subset()
        5: Progress bar over particles
        6: Shrink factor, etc., in align_filter_shrink()
    """
    
    group_pca= parser.add_argument_group(
        title='PCA options',
        description='Options when running principal component analysis.')
    
    group_pca.add_argument(
        '--nvec', 
        type=int, 
        default=0, 
        help='Number of eigenimages to compute')
    
    group_pca.add_argument(
        '--pca_radius', 
        type=int, 
        default=-1, 
        help='Radius for PCA (pixels, for full-sized images)')
    
    group_pca.add_argument(
        '--mask_binary', 
        type=str, 
        default=None, 
        help='Binary mask to use for PCA')
    
    group_pca.add_argument(
        '--mask_drawn', 
        type=str, 
        default=None, 
        help='Manually-drawm mask to use for PCA, drawn on top of image')
    
    group_center= parser.add_argument_group(
        title='Class-centering options',
        description='Options if class averages were centered separately.')
    
    group_center.add_argument(
        '--applyparams', 
        type=str, 
        default=None, 
        help='Parameters file to apply to images: none (default), combined, intshifts')
    
    group_center.add_argument(
        '--write_centered', 
        action="store_true", 
        help='Center using parameters from sp_center_2d3d.py')
    
    group_center.add_argument(
        '--debug', 
        action="store_true", 
        help='Debug mode')
    
    group_bandpass= parser.add_argument_group(
        title='Band-pass filtering options',
        description='Options if applying a band-pass filter to existing averages.')
    
    group_bandpass.add_argument(
        '--bandpass_radius', 
        type=float, 
        default=None, 
        help='Band-pass filter radius, Angstroms')
    
    group_bandpass.add_argument(
        '--bandpass_width', 
        type=float, 
        default=0.03, 
        help='Band-pass filter width (1/px, i.e., absolute frequency)')
    
    return parser.parse_args()

def main():
    options= parse_command_line()

    #print('options', options)
    #exit()
    
    RUNNING_UNDER_MPI= "OMPI_COMM_WORLD_SIZE" in os.environ
    if RUNNING_UNDER_MPI:
            mpi.mpi_init( 0, [] )  # On OS X, there is an error if MPI is initialized and not finalized, hence the conditional
    sp_global_def.print_timestamp( "Start" )
    
    eval_isac(options.classavgs, options, outdir=options.outdir, verbosity=options.verbosity, usempi=RUNNING_UNDER_MPI)
    
    sp_global_def.print_timestamp( "Finish" )
    if RUNNING_UNDER_MPI:
            mpi.mpi_finalize()

if __name__ == "__main__":
    main()
